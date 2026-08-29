"""Audit strict exports, their run freshness, and their scientific provenance as data."""

from __future__ import annotations

import json
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path, PurePosixPath, PureWindowsPath
from typing import Any

import pandas as pd

from nexus_genomics.adapters.base import Adapter, Availability
from nexus_genomics.common import (
    CONTENT_HASH_ALGORITHM,
    PACKAGE_VERSION,
    canonical_json,
    hash_file,
)
from nexus_genomics.nexus_csv import EXTENSION, INDEX_COLUMN, read_sidecar_and_shape
from nexus_genomics.projection import ordered_target_hash
from nexus_genomics.run_receipt import RECEIPT_FORMAT, RECEIPT_VERSION, receipt_path
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
    def has_failures(self) -> bool:
        return any(check.status is ComplianceStatus.FAIL for check in self.checks)

    @property
    def exit_ok(self) -> bool:
        return not self.has_failures

    @property
    def all_requirements_satisfied(self) -> bool:
        return all(check.status is ComplianceStatus.PASS for check in self.checks)

    @property
    def overall_status(self) -> str:
        if self.has_failures:
            return "FAIL"
        if self.all_requirements_satisfied:
            return "PASS"
        return "INCOMPLETE"

    @property
    def status_counts(self) -> dict[str, int]:
        return {
            status.value: sum(1 for check in self.checks if check.status is status)
            for status in ComplianceStatus
        }

    def to_json(self) -> str:
        """Publish distinct execution and requirement meanings; a bare ``ok`` is ambiguous."""
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
                "has_failures": self.has_failures,
                "exit_ok": self.exit_ok,
                "all_requirements_satisfied": self.all_requirements_satisfied,
                "overall_status": self.overall_status,
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
    manifest: Mapping[str, Any] | None
    validates: _Finding
    primary_key: _Finding
    single_target: _Finding
    features: _Finding
    encoding: _Finding
    provenance: _Finding
    projection_error: str | None


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
        detail: list[str] = []
        if missing:
            detail.append(f"checks did not run: {', '.join(missing)}")
        if failed:
            detail.append(f"failed: {', '.join(failed)}")
        return _Finding(False, "; ".join(detail))
    return _Finding(True, "required validator checks passed")


def _trusted_provenance(
    manifest: Mapping[str, Any],
    adapter: Adapter,
    validation: ValidationResult | None,
) -> _Finding:
    failures: list[str] = []
    expected = {
        "source": adapter.name,
        "source_url": adapter.source_url,
        "licence": adapter.licence,
    }
    for field, value in expected.items():
        if manifest.get(field) != value:
            failures.append(f"{field} differs from the adapter declaration")
    description = manifest.get("target_description")
    if not isinstance(description, str) or not description.startswith(adapter.target_description):
        failures.append("target_description does not begin with the adapter declaration")
    semantics = manifest.get("label_semantics")
    if not isinstance(semantics, Mapping) or not semantics:
        failures.append("label_semantics is empty or not an object")
    documented = _validation_finding(
        validation,
        ("labels_only_documented_values",),
        "validator did not return a result",
    )
    if not documented.passed:
        failures.append(
            f"label_semantics validator evidence is absent or failed ({documented.detail})"
        )
    return _Finding(
        not failures,
        "manifest and documented target values match the trusted adapter declaration"
        if not failures
        else "; ".join(failures),
    )


def _projection_shape_error(projection: object) -> str | None:
    if not isinstance(projection, Mapping):
        return (
            "target_projection metadata is missing"
            if projection is None
            else "target_projection metadata must be an object"
        )
    count = projection.get("original_target_count")
    index = projection.get("original_target_index")
    if type(count) is not int or count <= 0:
        return "target_projection must contain a positive integer original_target_count"
    if type(index) is not int:
        return "target_projection must contain an integer original_target_index"
    if not 0 <= index < count:
        return f"target_projection index {index} is outside declared range 0..{count - 1}"
    for field in ("source_target_names", "source_target_slugs", "source_target_hashes"):
        value = projection.get(field)
        if (
            not isinstance(value, list)
            or len(value) != count
            or not all(isinstance(item, str) and item for item in value)
        ):
            return f"target_projection {field} must contain {count} nonempty strings"
    if projection.get("source_target_hash_algorithm") != CONTENT_HASH_ALGORITHM:
        return f"target_projection source_target_hash_algorithm must be {CONTENT_HASH_ALGORITHM}"
    row_count = projection.get("source_row_count")
    if type(row_count) is not int or row_count <= 0:
        return "target_projection source_row_count must be a positive integer"
    name = projection.get("original_target_name")
    slug = projection.get("original_target_slug")
    if not isinstance(name, str) or not isinstance(slug, str):
        return "target_projection selected name and slug must be strings"
    return None


def _inspect_file(path: Path, adapter: Adapter) -> _FileEvidence:
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
    manifest: Mapping[str, Any] | None = None
    target_finding = _Finding(False, "manifest could not be read")
    feature_finding = _Finding(False, "manifest could not be read")
    encoding_finding = _Finding(False, "manifest could not be read")
    provenance_finding = _Finding(False, "manifest could not be read")
    projection_error: str | None = None
    sidecar_error = "manifest could not be read"
    sidecar_ok = False
    try:
        payload, header, _n_rows, _widths = read_sidecar_and_shape(path)
        columns: Mapping[str, Any] = payload["columns"]
        manifest = payload["manifest"]
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
            f"{len(feature_columns)} declared feature column(s); position names={names_ok}, "
            f"manifest order={order_ok}, integer={integer_features.passed}",
        )
        encoding = manifest.get("encoding")
        encoding_map = encoding if isinstance(encoding, Mapping) else {}
        padding = encoding_map.get("padding_token")
        encoding_finding = _Finding(
            manifest.get("encoding_mode") == "letter"
            and encoding_map.get("cipher") == "A=1 .. Z=26"
            and type(padding) is int
            and padding == 0,
            f"mode={manifest.get('encoding_mode')!r}, cipher={encoding_map.get('cipher')!r}, "
            f"padding={padding!r}",
        )
        provenance_finding = _trusted_provenance(manifest, adapter, validation)
        if "__" in path.name:
            projection_error = _projection_shape_error(manifest.get("target_projection"))
    except Exception as exc:
        sidecar_error = f"manifest read raised {type(exc).__name__}: {exc}"
        target_finding = _Finding(False, sidecar_error)
        feature_finding = _Finding(False, sidecar_error)
        encoding_finding = _Finding(False, sidecar_error)
        provenance_finding = _Finding(False, sidecar_error)
        projection_error = sidecar_error

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
        manifest,
        validates,
        primary_key,
        target_finding,
        feature_finding,
        encoding_finding,
        provenance_finding,
        projection_error,
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


def _receipt_freshness(
    out_dir: Path,
    profile: str,
    config_digest: str,
    enabled_source_keys: Sequence[str],
) -> ComplianceCheck:
    """Bind the receipt to disk: same config, same sources, same tables, same bytes.

    This re-reads and re-hashes every table independently of the per-file evidence gathered
    elsewhere in one ``compliance`` invocation, catching any output left over from an older or
    interrupted run. It does not hold a lock across the two reads, so a concurrent writer
    landing strictly between them could in principle let the two views disagree -- a narrow
    window (it needs a second process writing during this one invocation, not merely a stale
    file), accepted rather than solved with file locking for a single-user local tool.
    """
    path = receipt_path(out_dir)
    try:
        raw_receipt: object = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        return ComplianceCheck(
            "conversion receipt and output freshness",
            ComplianceStatus.FAIL,
            f"complete matching {path.name} is required: {type(exc).__name__}: {exc}",
        )
    if not isinstance(raw_receipt, Mapping):
        return ComplianceCheck(
            "conversion receipt and output freshness",
            ComplianceStatus.FAIL,
            f"complete matching {path.name} must contain a JSON object",
        )
    receipt: Mapping[str, Any] = raw_receipt

    errors: list[str] = []
    if receipt.get("format") != RECEIPT_FORMAT or receipt.get("version") != RECEIPT_VERSION:
        errors.append("format/version does not match the conversion-run contract")
    if receipt.get("generator") != f"nexus-genomics {PACKAGE_VERSION}":
        errors.append("generator version does not match this executable")
    if receipt.get("complete") is not True:
        errors.append("receipt is incomplete")
    if receipt.get("profile") != profile:
        errors.append(f"profile differs ({receipt.get('profile')!r} != {profile!r})")
    if receipt.get("target_layout") != "single":
        errors.append("target_layout is not the required single-target layout")
    if receipt.get("config_hash") != config_digest:
        errors.append("config_hash differs from the current raw configuration")
    enabled = receipt.get("enabled_source_keys")
    if enabled != list(enabled_source_keys):
        errors.append("enabled_source_keys differs from the current enabled registry set")

    converted = receipt.get("successfully_converted_source_keys")
    gated = receipt.get("externally_gated_source_keys")
    if not isinstance(converted, list) or not all(isinstance(key, str) for key in converted):
        errors.append("successfully_converted_source_keys is malformed")
        converted = []
    if not isinstance(gated, list) or not all(isinstance(key, str) for key in gated):
        errors.append("externally_gated_source_keys is malformed")
        gated = []
    if set(converted) & set(gated):
        errors.append("a source is both converted and externally gated")
    if set(converted) | set(gated) != set(enabled_source_keys):
        errors.append("converted and gated source keys do not exactly partition enabled sources")

    raw_entries = receipt.get("generated_tables")
    entries: list[tuple[str, str]] = []
    if not isinstance(raw_entries, list):
        errors.append("generated_tables is not a list")
    else:
        for entry in raw_entries:
            if not isinstance(entry, Mapping):
                errors.append("generated_tables contains a non-object entry")
                continue
            relative = entry.get("path")
            digest = entry.get("content_hash")
            if not isinstance(relative, str) or not isinstance(digest, str):
                errors.append("generated table entry lacks string path/content_hash")
                continue
            posix = PurePosixPath(relative)
            windows = PureWindowsPath(relative)
            # Checked under both path flavours, not just the host platform's own Path: a
            # Windows drive-letter path (C:/...) or UNC path is not "absolute" by POSIX rules,
            # so a POSIX-only check would wave it through wherever this audit happens to run.
            if (
                posix.is_absolute()
                or windows.is_absolute()
                or ".." in posix.parts
                or ".." in windows.parts
                or not relative.endswith(EXTENSION)
            ):
                errors.append(
                    f"generated table path is not a safe relative .ml.csv path: {relative}"
                )
                continue
            entries.append((relative, digest))

    receipt_paths = [relative for relative, _digest in entries]
    if receipt_paths != sorted(receipt_paths) or len(receipt_paths) != len(set(receipt_paths)):
        errors.append("generated table paths are not unique and sorted")
    discovered = sorted(
        path.relative_to(out_dir).as_posix() for path in out_dir.rglob(f"*{EXTENSION}")
    )
    if receipt_paths != discovered:
        missing = sorted(set(receipt_paths) - set(discovered))
        unexpected = sorted(set(discovered) - set(receipt_paths))
        errors.append(
            f"receipt inventory differs from disk; missing={missing[:3]}, "
            f"unexpected={unexpected[:3]}"
        )
    for relative, expected_hash in entries:
        table = out_dir / Path(relative)
        if table.is_file() and hash_file(table) != expected_hash:
            errors.append(f"content hash differs for {relative}")

    return ComplianceCheck(
        "conversion receipt and output freshness",
        ComplianceStatus.FAIL if errors else ComplianceStatus.PASS,
        "; ".join(errors)
        if errors
        else (
            "complete receipt matches config, enabled sources, and "
            f"{len(entries)} exact table hash(es)"
        ),
    )


def _projection_coverage(
    registry_key: str,
    adapter: Adapter,
    profile: str,
    projected: Sequence[_FileEvidence],
    config_digest: str,
) -> ComplianceCheck:
    requirement = f"source coverage: {registry_key}"
    malformed = [item for item in projected if item.projection_error is not None]
    if malformed:
        details = "; ".join(f"{item.path.name} ({item.projection_error})" for item in malformed[:3])
        return ComplianceCheck(
            requirement,
            ComplianceStatus.FAIL,
            f"{len(projected)} projected file(s); invalid projection metadata: {details}",
        )

    projections = [item.manifest["target_projection"] for item in projected if item.manifest]
    counts = sorted({projection["original_target_count"] for projection in projections})
    if len(counts) != 1:
        return ComplianceCheck(
            requirement,
            ComplianceStatus.FAIL,
            f"{len(projected)} projected file(s) have inconsistent "
            f"original_target_count values {counts}",
        )
    count = counts[0]
    indexes = [projection["original_target_index"] for projection in projections]
    duplicate_indexes = sorted({index for index in indexes if indexes.count(index) > 1})
    if duplicate_indexes:
        return ComplianceCheck(
            requirement,
            ComplianceStatus.FAIL,
            f"{len(projected)} projected file(s) have duplicate indexes {duplicate_indexes}",
        )
    expected_indexes = set(range(count))
    observed_indexes = set(indexes)
    if observed_indexes != expected_indexes:
        return ComplianceCheck(
            requirement,
            ComplianceStatus.FAIL,
            f"{len(projected)} projected file(s) do not cover 0..{count - 1}; missing indexes "
            f"{sorted(expected_indexes - observed_indexes)}, unexpected indexes "
            f"{sorted(observed_indexes - expected_indexes)}",
        )

    binding_contracts: set[str] = set()
    mapping_errors: list[str] = []
    hash_errors: list[str] = []
    for item, projection in zip(projected, projections, strict=True):
        manifest = item.manifest or {}
        index = projection["original_target_index"]
        names = projection["source_target_names"]
        slugs = projection["source_target_slugs"]
        hashes = projection["source_target_hashes"]
        if projection["original_target_name"] != names[index]:
            mapping_errors.append(f"{item.path.name} name/index mapping differs")
        if projection["original_target_slug"] != slugs[index]:
            mapping_errors.append(f"{item.path.name} slug/index mapping differs")
        expected_name = f"{adapter.name}_{profile}__{slugs[index]}{EXTENSION}"
        if item.path.name != expected_name:
            mapping_errors.append(f"{item.path.name} does not match selected slug {slugs[index]!r}")
        binding_contracts.add(
            json.dumps(
                {
                    "names": names,
                    "slugs": slugs,
                    "hashes": hashes,
                    "hash_algorithm": projection["source_target_hash_algorithm"],
                    "source_row_count": projection["source_row_count"],
                    "config_hash": manifest.get("config_hash"),
                    "raw_input_hashes": manifest.get("raw_input_hashes"),
                    "samples_content_hash": manifest.get("samples_content_hash"),
                    "samples_content_hash_algorithm": manifest.get(
                        "samples_content_hash_algorithm"
                    ),
                    "samples_n_rows": manifest.get("samples_n_rows"),
                },
                sort_keys=True,
                separators=(",", ":"),
            )
        )
        if manifest.get("config_hash") != config_digest:
            mapping_errors.append(f"{item.path.name} config_hash differs from current config")
        if (
            manifest.get("samples_n_rows") != projection["source_row_count"]
            or manifest.get("n_samples") != projection["source_row_count"]
        ):
            mapping_errors.append(f"{item.path.name} source/sample row counts differ")
        try:
            frame = pd.read_csv(
                item.path,
                usecols=[INDEX_COLUMN, "target"],
                dtype={INDEX_COLUMN: "string"},
            )
            if not pd.api.types.is_integer_dtype(frame["target"]):
                raise TypeError("target dtype is not integer")
            current_target_hash = ordered_target_hash(
                frame[INDEX_COLUMN].astype(str).tolist(), frame["target"].tolist()
            )
            if current_target_hash != hashes[index]:
                hash_errors.append(f"{item.path.name} selected target hash differs")
        except Exception as exc:
            hash_errors.append(f"{item.path.name} selected target hash failed: {exc}")

    if len(binding_contracts) != 1:
        mapping_errors.append("projected files do not share one source-level contract")
    if mapping_errors or hash_errors:
        return ComplianceCheck(
            requirement,
            ComplianceStatus.FAIL,
            "; ".join([*mapping_errors, *hash_errors][:6]),
        )
    return ComplianceCheck(
        requirement,
        ComplianceStatus.PASS,
        f"{len(projected)} projected file(s) provide the bound indexes, names, slugs, rows, "
        f"and target hashes 0..{count - 1}",
    )


def _present_source_coverage(
    registry_key: str,
    adapter: Adapter,
    profile: str,
    evidence: Sequence[_FileEvidence],
    config_digest: str,
) -> ComplianceCheck:
    requirement = f"source coverage: {registry_key}"
    native_name = f"{adapter.name}_{profile}{EXTENSION}"
    native = [item for item in evidence if item.path.name == native_name]
    projected = [item for item in evidence if item.path.name != native_name]
    if native and projected:
        return ComplianceCheck(
            requirement,
            ComplianceStatus.FAIL,
            f"source mixes native and projected outputs across {len(evidence)} file(s)",
        )
    if native:
        return ComplianceCheck(
            requirement,
            ComplianceStatus.PASS,
            f"1 native strict output file present for {adapter.name}_{profile}",
        )
    return _projection_coverage(registry_key, adapter, profile, projected, config_digest)


def audit_email_compliance(
    adapters: Mapping[str, Adapter],
    raw_dir: Path,
    strict_output_dir: Path,
    profile: str,
    *,
    config_digest: str,
    enabled_source_keys: Sequence[str],
) -> ComplianceReport:
    """Audit only a complete current strict run, while keeping external blockers explicit."""
    enabled_set = set(enabled_source_keys)
    files_by_source = {
        registry_key: _discover(strict_output_dir, adapter.name, profile)
        for registry_key, adapter in adapters.items()
        if registry_key in enabled_set
    }
    evidence_by_source = {
        registry_key: [_inspect_file(path, adapters[registry_key]) for path in paths]
        for registry_key, paths in files_by_source.items()
    }
    evidence = [
        item
        for registry_key in enabled_source_keys
        for item in evidence_by_source.get(registry_key, [])
    ]

    availability_by_source: dict[str, Availability | Exception] = {}
    for registry_key, adapter in adapters.items():
        try:
            availability_by_source[registry_key] = adapter.availability(raw_dir / registry_key)
        except Exception as exc:
            availability_by_source[registry_key] = exc

    coverage: list[ComplianceCheck] = []
    missing_ready_source = False
    for registry_key in enabled_source_keys:
        adapter = adapters[registry_key]
        paths = files_by_source.get(registry_key, [])
        availability = availability_by_source[registry_key]
        if isinstance(availability, Exception):
            missing_ready_source = True
            coverage.append(
                ComplianceCheck(
                    f"source coverage: {registry_key}",
                    ComplianceStatus.FAIL,
                    f"availability raised {type(availability).__name__}: {availability}",
                )
            )
            continue
        if not availability.ready:
            context = availability.detail
            if availability.manual_steps:
                context += f" Manual steps: {availability.manual_steps}"
            stale = (
                " Existing output is stale and is not counted as current coverage." if paths else ""
            )
            coverage.append(
                ComplianceCheck(
                    f"source coverage: {registry_key}",
                    ComplianceStatus.BLOCKED_EXTERNAL,
                    f"source is not currently ready: {context}.{stale}",
                )
            )
            continue
        if not paths:
            missing_ready_source = True
            coverage.append(
                ComplianceCheck(
                    f"source coverage: {registry_key}",
                    ComplianceStatus.FAIL,
                    f"output absent although the source is ready: {availability.detail}",
                )
            )
            continue
        coverage.append(
            _present_source_coverage(
                registry_key,
                adapter,
                profile,
                evidence_by_source[registry_key],
                config_digest,
            )
        )

    empty_status = (
        ComplianceStatus.FAIL if missing_ready_source else ComplianceStatus.BLOCKED_EXTERNAL
    )
    checks = [
        _receipt_freshness(
            strict_output_dir,
            profile,
            config_digest,
            enabled_source_keys,
        ),
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
        _aggregate(
            "label provenance matches trusted adapter declarations",
            evidence,
            lambda item: item.provenance,
            empty_status,
        ),
        *coverage,
        ComplianceCheck(
            "natural-versus-engineered target",
            ComplianceStatus.UNSUPPORTED_BY_SOURCE,
            "None of the four original sources supplies a general natural-versus-engineered "
            "target. The fifth felix_signatures source supplies only a narrower recorded-event "
            "contrast (host context not introduced by that event versus sequence introduced by "
            "that event); it does not establish a general benchmark.",
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
        f"**Result**: {report.overall_status}",
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
