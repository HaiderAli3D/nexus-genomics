"""Audit the email's strict table requirements without overstating scientific evidence."""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from typing import Any

from nexus_genomics.adapters.base import Adapter
from nexus_genomics.common import canonical_json
from nexus_genomics.nexus_csv import EXTENSION, read_sidecar_and_shape
from nexus_genomics.validation import ValidationResult, validate_file

__all__ = [
    "ComplianceCheck",
    "ComplianceReport",
    "ComplianceStatus",
    "audit_email_compliance",
    "render_markdown",
]


class ComplianceStatus(StrEnum):
    PASS = "PASS"
    BLOCKED_EXTERNAL = "BLOCKED_EXTERNAL"
    UNSUPPORTED_BY_SOURCE = "UNSUPPORTED_BY_SOURCE"
    FAIL = "FAIL"


@dataclass(frozen=True, slots=True)
class ComplianceCheck:
    requirement: str
    status: ComplianceStatus
    detail: str


@dataclass(slots=True)
class ComplianceReport:
    checks: Sequence[ComplianceCheck]

    @property
    def ok(self) -> bool:
        return all(check.status is not ComplianceStatus.FAIL for check in self.checks)

    @property
    def status_counts(self) -> dict[str, int]:
        return {
            status.value: sum(1 for check in self.checks if check.status is status)
            for status in ComplianceStatus
        }

    def to_json(self) -> str:
        """The report's canonical machine-readable representation."""
        return canonical_json(
            {
                "checks": [
                    {
                        "detail": check.detail,
                        "requirement": check.requirement,
                        "status": check.status.value,
                    }
                    for check in self.checks
                ],
                "ok": self.ok,
                "status_counts": self.status_counts,
            }
        )


@dataclass(frozen=True, slots=True)
class _Finding:
    passed: bool
    detail: str


@dataclass(frozen=True, slots=True)
class _FileEvidence:
    path: Path
    validates: _Finding
    primary_key: _Finding
    single_target: _Finding
    features: _Finding
    encoding: _Finding


def _validation_finding(
    result: ValidationResult | None,
    required_names: tuple[str, ...],
    error: str,
) -> _Finding:
    if result is None:
        return _Finding(False, error)
    by_name = {check.name: check for check in result.checks}
    missing = [name for name in required_names if name not in by_name]
    failed = [name for name in required_names if name in by_name and not by_name[name].passed]
    if missing or failed:
        detail = []
        if missing:
            detail.append(f"checks did not run: {', '.join(missing)}")
        if failed:
            detail.append(f"failed: {', '.join(failed)}")
        return _Finding(False, "; ".join(detail))
    return _Finding(True, "required validator checks passed")


def _inspect_file(path: Path) -> _FileEvidence:
    validation: ValidationResult | None = None
    validation_error = "validator did not return a result"
    try:
        validation = validate_file(path)
    except Exception as exc:
        validation_error = f"validator raised {type(exc).__name__}: {exc}"

    primary_key = _validation_finding(
        validation,
        ("primary_key_column_exists", "primary_key_non_null", "primary_key_unique"),
        validation_error,
    )
    integer_features = _validation_finding(
        validation,
        ("all_position_columns_integer",),
        validation_error,
    )

    sidecar_error = "manifest could not be read"
    target_finding = _Finding(False, sidecar_error)
    feature_finding = _Finding(False, sidecar_error)
    encoding_finding = _Finding(False, sidecar_error)
    sidecar_ok = False
    try:
        payload, header, _n_rows, _widths = read_sidecar_and_shape(path)
        columns: dict[str, Any] = payload["columns"]
        manifest: dict[str, Any] = payload["manifest"]
        index_columns = list(columns["index_columns"])
        target_columns = list(columns["target_columns"])
        feature_columns = list(columns["feature_columns"])
        sidecar_ok = True

        target_finding = _Finding(
            len(target_columns) == 1,
            f"manifest declares {len(target_columns)} target column(s)",
        )

        feature_start = len(index_columns) + len(target_columns)
        names_ok = all(
            isinstance(name, str) and name.startswith("position_") for name in feature_columns
        )
        order_ok = header[feature_start:] == feature_columns
        feature_finding = _Finding(
            names_ok and order_ok and integer_features.passed,
            f"{len(feature_columns)} declared feature column(s); "
            f"position names={names_ok}, manifest order={order_ok}, "
            f"integer={integer_features.passed}",
        )

        encoding = manifest.get("encoding")
        encoding_map = encoding if isinstance(encoding, dict) else {}
        padding = encoding_map.get("padding_token")
        encoding_finding = _Finding(
            manifest.get("encoding_mode") == "letter"
            and encoding_map.get("cipher") == "A=1 .. Z=26"
            and type(padding) is int
            and padding == 0,
            f"mode={manifest.get('encoding_mode')!r}, cipher="
            f"{encoding_map.get('cipher')!r}, padding={padding!r}",
        )
    except Exception as exc:
        sidecar_error = f"manifest read raised {type(exc).__name__}: {exc}"
        target_finding = _Finding(False, sidecar_error)
        feature_finding = _Finding(False, sidecar_error)
        encoding_finding = _Finding(False, sidecar_error)

    if validation is None:
        validates = _Finding(False, validation_error)
    elif not validation.ok:
        failed = [check.name for check in validation.checks if not check.passed]
        validates = _Finding(False, f"validator failed: {', '.join(failed)}")
    elif not sidecar_ok:
        validates = _Finding(False, sidecar_error)
    else:
        validates = _Finding(True, f"all {len(validation.checks)} validator checks passed")

    return _FileEvidence(
        path,
        validates,
        primary_key,
        target_finding,
        feature_finding,
        encoding_finding,
    )


def _discover(out_dir: Path, adapter_name: str, profile: str) -> list[Path]:
    if not out_dir.is_dir():
        return []
    native = f"{adapter_name}_{profile}{EXTENSION}"
    projected_prefix = f"{adapter_name}_{profile}__"
    return sorted(
        path
        for path in out_dir.iterdir()
        if path.is_file()
        and (
            path.name == native
            or (path.name.startswith(projected_prefix) and path.name.endswith(EXTENSION))
        )
    )


def _aggregate(
    requirement: str,
    evidence: Sequence[_FileEvidence],
    finding_of: Callable[[_FileEvidence], _Finding],
    empty_status: ComplianceStatus,
) -> ComplianceCheck:
    if not evidence:
        return ComplianceCheck(
            requirement,
            empty_status,
            "0 file(s) checked; no exact native or projected output exists for this profile",
        )
    failed = [
        (item.path.name, finding_of(item).detail)
        for item in evidence
        if not finding_of(item).passed
    ]
    if failed:
        examples = "; ".join(f"{name} ({detail})" for name, detail in failed[:3])
        return ComplianceCheck(
            requirement,
            ComplianceStatus.FAIL,
            f"{len(evidence)} file(s) checked; {len(failed)} failed: {examples}",
        )
    return ComplianceCheck(
        requirement,
        ComplianceStatus.PASS,
        f"{len(evidence)} file(s) checked; requirement established for every file",
    )


def audit_email_compliance(
    adapters: Mapping[str, Adapter],
    raw_dir: Path,
    strict_output_dir: Path,
    profile: str,
) -> ComplianceReport:
    """Audit exact strict outputs using injected adapters so the audit remains offline."""
    files_by_source = {
        registry_key: _discover(strict_output_dir, adapter.name, profile)
        for registry_key, adapter in adapters.items()
    }
    evidence = [
        _inspect_file(path) for registry_key in adapters for path in files_by_source[registry_key]
    ]

    coverage: list[ComplianceCheck] = []
    missing_ready_source = False
    for registry_key, adapter in adapters.items():
        paths = files_by_source[registry_key]
        if paths:
            coverage.append(
                ComplianceCheck(
                    f"source coverage: {registry_key}",
                    ComplianceStatus.PASS,
                    f"{len(paths)} strict output file(s) present for {adapter.name}_{profile}",
                )
            )
            continue
        try:
            availability = adapter.availability(raw_dir / registry_key)
        except Exception as exc:
            missing_ready_source = True
            coverage.append(
                ComplianceCheck(
                    f"source coverage: {registry_key}",
                    ComplianceStatus.FAIL,
                    f"output absent and availability raised {type(exc).__name__}: {exc}",
                )
            )
            continue
        if availability.ready:
            missing_ready_source = True
            coverage.append(
                ComplianceCheck(
                    f"source coverage: {registry_key}",
                    ComplianceStatus.FAIL,
                    f"output absent although the source is ready: {availability.detail}",
                )
            )
            continue
        context = availability.detail
        if availability.manual_steps:
            context += f" Manual steps: {availability.manual_steps}"
        coverage.append(
            ComplianceCheck(
                f"source coverage: {registry_key}",
                ComplianceStatus.BLOCKED_EXTERNAL,
                f"output absent because the source is not ready: {context}",
            )
        )

    empty_status = (
        ComplianceStatus.FAIL if missing_ready_source else ComplianceStatus.BLOCKED_EXTERNAL
    )
    checks = [
        _aggregate(
            "CSV loads and internally validates",
            evidence,
            lambda item: item.validates,
            empty_status,
        ),
        _aggregate(
            "sample_id is present, non-null, and unique",
            evidence,
            lambda item: item.primary_key,
            empty_status,
        ),
        _aggregate(
            "exactly one target column exists",
            evidence,
            lambda item: item.single_target,
            empty_status,
        ),
        _aggregate(
            "all declared features are position_*, use the manifest's literal order, "
            "and are integer",
            evidence,
            lambda item: item.features,
            empty_status,
        ),
        _aggregate(
            "encoding mode is letter, cipher is exactly A=1 .. Z=26, and padding token is 0",
            evidence,
            lambda item: item.encoding,
            empty_status,
        ),
        *coverage,
        ComplianceCheck(
            "natural-versus-engineered target",
            ComplianceStatus.UNSUPPORTED_BY_SOURCE,
            "None of the four sources is a general natural-versus-engineered benchmark. "
            "CodonTransformer compares natural NCBI coding sequences only with output generated "
            "by one model, not with engineered sequences in general.",
        ),
        ComplianceCheck(
            "verified Nexus compatibility",
            ComplianceStatus.BLOCKED_EXTERNAL,
            "No authoritative Nexus fixture or schema has been supplied. Matching the requested "
            "strict CSV shape is not evidence that Nexus accepts the files.",
        ),
    ]
    return ComplianceReport(tuple(checks))


def render_markdown(report: ComplianceReport) -> str:
    """Render the same report object used by the canonical JSON output."""
    lines = [
        "# Email compliance report",
        "",
        f"**Result**: {'PASS' if report.ok else 'FAIL'}",
        "",
        "## Status counts",
        "",
        "| Status | Count |",
        "|---|---:|",
    ]
    lines.extend(f"| {status} | {count} |" for status, count in report.status_counts.items())
    lines += [
        "",
        "## Checks",
        "",
        "| Requirement | Status | Detail |",
        "|---|---|---|",
    ]
    for check in report.checks:
        requirement = check.requirement.replace("|", "\\|").replace("\n", " ")
        detail = check.detail.replace("|", "\\|").replace("\n", " ")
        lines.append(f"| {requirement} | {check.status.value} | {detail} |")
    return "\n".join(lines) + "\n"
