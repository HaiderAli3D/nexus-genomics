"""Executable evidence for every requirement in the email compliance report."""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

import pytest
from typer.testing import CliRunner

from nexus_genomics.adapters.base import Availability, LoadedSource
from nexus_genomics.cleaning import Alphabet
from nexus_genomics.cli import app
from nexus_genomics.common import SampleRecord, canonical_json, hash_file, hash_text
from nexus_genomics.pipeline import ConvertOptions, convert

if TYPE_CHECKING:
    from nexus_genomics.compliance import ComplianceReport

runner = CliRunner()


class _Adapter:
    source_url = "https://example.invalid/compliance"
    licence = "Synthetic compliance fixture; no redistribution restrictions."
    target_description = "A documented synthetic target used only for offline compliance tests."

    def __init__(
        self,
        name: str,
        *,
        ready: bool,
        detail: str = "synthetic source is ready",
        manual_steps: str = "",
    ) -> None:
        self.name = name
        self._availability = Availability(ready, detail, manual_steps)
        self.availability_calls = 0

    def availability(self, raw_dir: Path) -> Availability:
        _ = raw_dir
        self.availability_calls += 1
        return self._availability

    def load(self, raw_dir: Path, options: Mapping[str, Any]) -> LoadedSource:
        _ = raw_dir, options
        raise AssertionError("the compliance audit must not load source data")


def _write_table(
    path: Path,
    *,
    n_targets: int = 1,
    labels: list[tuple[int, ...]] | None = None,
    manifest_extra: Mapping[str, Any] | None = None,
    sample_ids: tuple[str, str] = ("s1", "s2"),
) -> None:
    resolved_labels = labels or ([(0,), (1,)] if n_targets == 1 else [(1, 0), (0, 1)])
    convert(
        LoadedSource(
            records=[
                SampleRecord(sample_ids[0], "ACGT", resolved_labels[0], "fixture", {}),
                SampleRecord(sample_ids[1], "TGCA", resolved_labels[1], "fixture", {}),
            ],
            n_targets=n_targets,
            label_semantics={0: "absent", 1: "present"},
            alphabet=Alphabet.DNA_NUCLEOTIDE,
        ),
        path,
        ConvertOptions(width=4),
        source_name="fixture_output",
        source_url=_Adapter.source_url,
        licence=_Adapter.licence,
        target_description=_Adapter.target_description,
        notes="Synthetic compliance fixture",
        config_digest="fixture",
        extra_manifest=manifest_extra,
    )


def _write_projection(
    path: Path,
    *,
    index: int,
    count: int,
    sample_ids: tuple[str, str] = ("s1", "s2"),
) -> None:
    names = ["Target Zero", "Target One", "Target Two"][:count]
    slugs = ["target_zero", "target_one", "target_two"][:count]
    source_labels = [
        tuple(1 if target_index == row_index else 0 for target_index in range(count))
        for row_index in range(2)
    ]

    def target_hash(target_index: int) -> str:
        pairs = [
            [sample_ids[row_index], row[target_index]]
            for row_index, row in enumerate(source_labels)
        ]
        encoded = json.dumps(pairs, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
        return hashlib.blake2b(encoded, digest_size=32).hexdigest()

    _write_table(
        path,
        labels=[(source_labels[0][index],), (source_labels[1][index],)],
        manifest_extra={
            "target_layout": "single",
            "target_projection": {
                "original_target_count": count,
                "original_target_index": index,
                "original_target_name": names[index],
                "original_target_slug": slugs[index],
                "source_target_names": names,
                "source_target_slugs": slugs,
                "source_target_hashes": [target_hash(i) for i in range(count)],
                "source_target_hash_algorithm": "blake2b-256",
                "source_row_count": 2,
            },
        },
        sample_ids=sample_ids,
    )


def _write_receipt(
    out_dir: Path,
    adapters: Mapping[str, _Adapter],
    *,
    config_hash: str = "fixture",
    complete: bool = True,
    generated_paths: list[Path] | None = None,
    gated_source_keys: list[str] | None = None,
) -> None:
    paths = generated_paths
    if paths is None:
        paths = sorted(out_dir.glob("*.ml.csv"))
    generated = []
    for path in paths:
        manifest_path = path.with_suffix(".manifest.json")
        if not manifest_path.exists():
            continue
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
        generated.append(
            {
                "path": path.relative_to(out_dir).as_posix(),
                "content_hash": payload["manifest"]["content_hash"],
            }
        )
    receipt = {
        "format": "nexus-genomics-conversion-run",
        "version": 1,
        "generator": "nexus-genomics 0.1.0",
        "complete": complete,
        "profile": "demo",
        "target_layout": "single",
        "config_hash": config_hash,
        "enabled_source_keys": list(adapters),
        "successfully_converted_source_keys": (
            [
                key
                for key, adapter in adapters.items()
                if any(path.name.startswith(adapter.name) for path in paths)
            ]
        ),
        "externally_gated_source_keys": gated_source_keys or [],
        "generated_tables": generated,
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / ".conversion-run.json").write_text(canonical_json(receipt), encoding="utf-8")


def _audit(
    adapters: Mapping[str, _Adapter],
    raw_dir: Path,
    out_dir: Path,
    profile: str = "demo",
    config_digest: str = "fixture",
) -> ComplianceReport:
    from nexus_genomics.compliance import audit_email_compliance

    if not (out_dir / ".conversion-run.json").exists():
        _write_receipt(out_dir, adapters)
    return audit_email_compliance(
        adapters,
        raw_dir,
        out_dir,
        profile,
        config_digest=config_digest,
        enabled_source_keys=tuple(adapters),
    )


def _write_config(path: Path, adapters: Mapping[str, _Adapter]) -> str:
    payload = {
        "encoding": {"mode": "letter"},
        "output": {"target_layout": "native"},
        "window": {
            "length": 4,
            "policy": "pad_or_truncate",
            "stride": 4,
            "max_windows_per_sample": 1,
        },
        "validation": {"round_trip_sample_size": 2, "round_trip_seed": 0},
        "sources": {key: {"enabled": True, "demo": {}, "full": {}} for key in adapters},
    }
    text = json.dumps(payload)
    path.write_text(text, encoding="utf-8")
    return hash_text(text)


def _checks(report: ComplianceReport) -> dict[str, Any]:
    return {check.requirement: check for check in report.checks}


def test_compliance_models_separate_failures_exit_policy_and_requirement_completion() -> None:
    """Collapsing blockers into failures would make the audit misstate external constraints."""
    from nexus_genomics.compliance import ComplianceCheck, ComplianceReport, ComplianceStatus

    assert [status.value for status in ComplianceStatus] == [
        "PASS",
        "BLOCKED_EXTERNAL",
        "UNSUPPORTED_BY_SOURCE",
        "FAIL",
    ]
    checks = tuple(
        ComplianceCheck(status.value, status, f"{status.value} detail")
        for status in ComplianceStatus
    )
    report = ComplianceReport(checks)
    assert report.has_failures is True
    assert report.exit_ok is False
    assert report.all_requirements_satisfied is False
    assert report.overall_status == "FAIL"
    assert report.status_counts == {
        "PASS": 1,
        "BLOCKED_EXTERNAL": 1,
        "UNSUPPORTED_BY_SOURCE": 1,
        "FAIL": 1,
    }
    blockers_only = ComplianceReport(checks[1:3])
    assert blockers_only.has_failures is False
    assert blockers_only.exit_ok is True
    assert blockers_only.all_requirements_satisfied is False
    assert blockers_only.overall_status == "INCOMPLETE"
    passes_only = ComplianceReport((checks[0],))
    assert passes_only.all_requirements_satisfied is True
    assert passes_only.overall_status == "PASS"


def test_audit_discovers_only_exact_native_and_projected_profile_names(tmp_path: Path) -> None:
    """A loose glob could grade stale profiles or unrelated tables as current strict output."""
    out_dir = tmp_path / "strict"
    _write_table(out_dir / "fixture_output_demo.ml.csv")
    out_dir.joinpath("fixture_output_full.ml.csv").write_text("broken", encoding="utf-8")
    out_dir.joinpath("fixture_output_full__selected_target.ml.csv").write_text(
        "broken", encoding="utf-8"
    )
    out_dir.joinpath("fixture_output_demo_extra.ml.csv").write_text("broken", encoding="utf-8")
    out_dir.joinpath("unrelated_demo.ml.csv").write_text("broken", encoding="utf-8")

    report = _audit({"fixture": _Adapter("fixture_output", ready=True)}, tmp_path / "raw", out_dir)

    assert report.has_failures is True
    assert _checks(report)["conversion receipt and output freshness"].status.value == "FAIL"
    structural = report.checks[1:6]
    assert all(check.status.value == "PASS" for check in structural)
    assert all("1 file(s)" in check.detail for check in structural)


def test_complete_two_target_projection_set_passes_source_coverage(tmp_path: Path) -> None:
    """A complete projection fan-out must remain accepted after completeness enforcement."""
    out_dir = tmp_path / "strict"
    _write_projection(out_dir / "fixture_output_demo__target_zero.ml.csv", index=0, count=2)
    _write_projection(out_dir / "fixture_output_demo__target_one.ml.csv", index=1, count=2)

    report = _audit({"fixture": _Adapter("fixture_output", ready=True)}, tmp_path / "raw", out_dir)
    coverage = _checks(report)["source coverage: fixture"]

    assert coverage.status.value == "PASS"
    assert "target hashes 0..1" in coverage.detail
    assert report.has_failures is False


def test_projection_hashing_preserves_numeric_looking_sample_id_text(tmp_path: Path) -> None:
    """Pandas inference must not turn source key ``001`` into ``1`` before target hashing."""
    out_dir = tmp_path / "strict"
    sample_ids = ("001", "002")
    _write_projection(
        out_dir / "fixture_output_demo__target_zero.ml.csv",
        index=0,
        count=2,
        sample_ids=sample_ids,
    )
    _write_projection(
        out_dir / "fixture_output_demo__target_one.ml.csv",
        index=1,
        count=2,
        sample_ids=sample_ids,
    )

    report = _audit({"fixture": _Adapter("fixture_output", ready=True)}, tmp_path / "raw", out_dir)

    assert _checks(report)["source coverage: fixture"].status.value == "PASS"


def test_deleted_projection_makes_source_coverage_fail(tmp_path: Path) -> None:
    """One surviving FunSoC projection must not represent a complete 32-target strict export."""
    out_dir = tmp_path / "strict"
    first = out_dir / "fixture_output_demo__target_zero.ml.csv"
    second = out_dir / "fixture_output_demo__target_one.ml.csv"
    _write_projection(first, index=0, count=2)
    _write_projection(second, index=1, count=2)
    second.unlink()

    report = _audit({"fixture": _Adapter("fixture_output", ready=True)}, tmp_path / "raw", out_dir)
    coverage = _checks(report)["source coverage: fixture"]

    assert coverage.status.value == "FAIL"
    assert "missing indexes [1]" in coverage.detail
    assert report.has_failures is True


def test_projected_output_without_projection_metadata_fails_coverage(tmp_path: Path) -> None:
    """A projected-looking filename alone cannot establish which target it preserves."""
    out_dir = tmp_path / "strict"
    _write_table(out_dir / "fixture_output_demo__unknown.ml.csv")

    report = _audit({"fixture": _Adapter("fixture_output", ready=True)}, tmp_path / "raw", out_dir)
    coverage = _checks(report)["source coverage: fixture"]

    assert coverage.status.value == "FAIL"
    assert "target_projection metadata is missing" in coverage.detail
    assert report.has_failures is True


def test_duplicate_projection_indexes_fail_coverage(tmp_path: Path) -> None:
    """Two filenames for one target must not conceal that another target was omitted."""
    out_dir = tmp_path / "strict"
    _write_projection(out_dir / "fixture_output_demo__first.ml.csv", index=0, count=2)
    _write_projection(out_dir / "fixture_output_demo__duplicate.ml.csv", index=0, count=2)

    report = _audit({"fixture": _Adapter("fixture_output", ready=True)}, tmp_path / "raw", out_dir)
    coverage = _checks(report)["source coverage: fixture"]

    assert coverage.status.value == "FAIL"
    assert "duplicate indexes [0]" in coverage.detail
    assert report.has_failures is True


def test_inconsistent_projection_counts_fail_coverage(tmp_path: Path) -> None:
    """Files that disagree on source target width cannot define a complete projection set."""
    out_dir = tmp_path / "strict"
    _write_projection(out_dir / "fixture_output_demo__first.ml.csv", index=0, count=2)
    _write_projection(out_dir / "fixture_output_demo__second.ml.csv", index=1, count=3)

    report = _audit({"fixture": _Adapter("fixture_output", ready=True)}, tmp_path / "raw", out_dir)
    coverage = _checks(report)["source coverage: fixture"]

    assert coverage.status.value == "FAIL"
    assert "inconsistent original_target_count values [2, 3]" in coverage.detail
    assert report.has_failures is True


def test_malformed_projection_metadata_fails_without_crashing(tmp_path: Path) -> None:
    """Boolean and string metadata must not be coerced into apparently valid indexes."""
    out_dir = tmp_path / "strict"
    _write_table(
        out_dir / "fixture_output_demo__malformed.ml.csv",
        manifest_extra={
            "target_projection": {
                "original_target_count": True,
                "original_target_index": "0",
                "original_target_name": "Malformed",
            }
        },
    )

    report = _audit({"fixture": _Adapter("fixture_output", ready=True)}, tmp_path / "raw", out_dir)
    coverage = _checks(report)["source coverage: fixture"]

    assert coverage.status.value == "FAIL"
    assert "positive integer original_target_count" in coverage.detail
    assert report.has_failures is True


def test_mixing_native_and_projected_outputs_fails_source_coverage(tmp_path: Path) -> None:
    """A stale native file beside strict fan-out must not count as one coherent export mode."""
    out_dir = tmp_path / "strict"
    _write_table(out_dir / "fixture_output_demo.ml.csv")
    _write_projection(out_dir / "fixture_output_demo__target_zero.ml.csv", index=0, count=1)

    report = _audit({"fixture": _Adapter("fixture_output", ready=True)}, tmp_path / "raw", out_dir)
    coverage = _checks(report)["source coverage: fixture"]

    assert coverage.status.value == "FAIL"
    assert "mixes native and projected outputs" in coverage.detail
    assert report.has_failures is True


def test_audit_fails_a_native_multitarget_table_even_when_internal_validation_passes(
    tmp_path: Path,
) -> None:
    """Native FunSoC shape is valid internally but violates the email's one-target contract."""
    out_dir = tmp_path / "strict"
    _write_table(out_dir / "fixture_output_demo.ml.csv", n_targets=2)

    report = _audit({"fixture": _Adapter("fixture_output", ready=True)}, tmp_path / "raw", out_dir)
    checks = _checks(report)

    assert checks["CSV loads and internally validates"].status.value == "PASS"
    assert checks["exactly one target column exists"].status.value == "FAIL"
    assert report.has_failures is True


def test_audit_turns_an_invalid_discovered_table_into_a_failure_without_crashing(
    tmp_path: Path,
) -> None:
    """Corruption is the audit's input, so it must become evidence rather than an exception."""
    out_dir = tmp_path / "strict"
    table = out_dir / "fixture_output_demo.ml.csv"
    _write_table(table)
    lines = table.read_text(encoding="utf-8").splitlines()
    table.write_text("\n".join([*lines, lines[1]]) + "\n", encoding="utf-8")

    report = _audit({"fixture": _Adapter("fixture_output", ready=True)}, tmp_path / "raw", out_dir)
    checks = _checks(report)

    assert checks["CSV loads and internally validates"].status.value == "FAIL"
    assert checks["sample_id is present, non-null, and unique"].status.value == "FAIL"
    assert "fixture_output_demo.ml.csv" in checks["CSV loads and internally validates"].detail


def test_audit_requires_declared_features_to_use_position_names(tmp_path: Path) -> None:
    """A manifest and table can agree on a wrong feature name while the validator stays green."""
    out_dir = tmp_path / "strict"
    table = out_dir / "fixture_output_demo.ml.csv"
    _write_table(table)
    lines = table.read_text(encoding="utf-8").splitlines()
    lines[0] = lines[0].replace("position_1", "variable_1")
    table.write_text("\n".join(lines) + "\n", encoding="utf-8")
    manifest_path = table.with_suffix(".manifest.json")
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    payload["columns"]["feature_columns"][0] = "variable_1"
    payload["manifest"]["content_hash"] = hash_file(table)
    manifest_path.write_text(canonical_json(payload), encoding="utf-8")

    report = _audit({"fixture": _Adapter("fixture_output", ready=True)}, tmp_path / "raw", out_dir)
    checks = _checks(report)

    requirement = (
        "all declared features are position_*, use the manifest's literal order, and are integer"
    )
    assert checks["CSV loads and internally validates"].status.value == "PASS"
    assert checks[requirement].status.value == "FAIL"


def test_audit_requires_the_exact_letter_cipher_metadata(tmp_path: Path) -> None:
    """Strict table bytes are insufficient if the manifest describes a different encoding."""
    out_dir = tmp_path / "strict"
    table = out_dir / "fixture_output_demo.ml.csv"
    _write_table(table)
    manifest_path = table.with_suffix(".manifest.json")
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    payload["manifest"]["encoding"]["cipher"] = "A=0 .. Z=25"
    manifest_path.write_text(canonical_json(payload), encoding="utf-8")

    report = _audit({"fixture": _Adapter("fixture_output", ready=True)}, tmp_path / "raw", out_dir)
    checks = _checks(report)

    requirement = "encoding mode is letter, cipher is exactly A=1 .. Z=26, and padding token is 0"
    assert checks["CSV loads and internally validates"].status.value == "PASS"
    assert checks[requirement].status.value == "FAIL"


def test_missing_output_distinguishes_ready_sources_from_gated_sources(tmp_path: Path) -> None:
    """Missing authorized input is a blocker; missing output from ready input is a defect."""
    report = _audit(
        {
            "ready_key": _Adapter("ready_output", ready=True, detail="raw fixture exists"),
            "gated_key": _Adapter(
                "gated_output",
                ready=False,
                detail="license acceptance required",
                manual_steps="Sign in and download gated.csv",
            ),
        },
        tmp_path / "raw",
        tmp_path / "strict",
    )
    checks = _checks(report)

    assert checks["source coverage: ready_key"].status.value == "FAIL"
    blocked = checks["source coverage: gated_key"]
    assert blocked.status.value == "BLOCKED_EXTERNAL"
    assert "license acceptance required" in blocked.detail
    assert "Sign in and download gated.csv" in blocked.detail
    assert report.has_failures is True


def test_a_stale_config_hash_fails_the_conversion_receipt_check(tmp_path: Path) -> None:
    """A coherent table from yesterday is not evidence for today's byte-affecting config."""
    out_dir = tmp_path / "strict"
    adapters = {"fixture": _Adapter("fixture_output", ready=True)}
    _write_table(out_dir / "fixture_output_demo.ml.csv")
    _write_receipt(out_dir, adapters, config_hash="old-config")

    report = _audit(adapters, tmp_path / "raw", out_dir, config_digest="current-config")

    freshness = _checks(report)["conversion receipt and output freshness"]
    assert freshness.status.value == "FAIL"
    assert "config_hash" in freshness.detail


def test_an_incomplete_conversion_receipt_fails_even_when_old_tables_validate(
    tmp_path: Path,
) -> None:
    """A failed rebuild must not certify the valid triplets left by its predecessor."""
    out_dir = tmp_path / "strict"
    adapters = {"fixture": _Adapter("fixture_output", ready=True)}
    _write_table(out_dir / "fixture_output_demo.ml.csv")
    _write_receipt(out_dir, adapters, complete=False)

    report = _audit(adapters, tmp_path / "raw", out_dir)

    freshness = _checks(report)["conversion receipt and output freshness"]
    assert freshness.status.value == "FAIL"
    assert "incomplete" in freshness.detail


def test_a_windows_drive_letter_path_in_the_receipt_is_rejected_as_unsafe(tmp_path: Path) -> None:
    """A POSIX-only absolute-path check would wave a drive-letter path through as "relative".

    ``PurePosixPath("C:/Windows/...").is_absolute()`` is ``False``, so a receipt entry naming
    an absolute Windows path used to slip past the safety gate and reach the join/hash step,
    even though it plainly does not live inside ``out_dir``. Checked under ``PureWindowsPath``
    too, since this audit is not guaranteed to run on the same platform the path string names.
    """
    out_dir = tmp_path / "strict"
    adapters = {"fixture": _Adapter("fixture_output", ready=True)}
    _write_table(out_dir / "fixture_output_demo.ml.csv")
    _write_receipt(out_dir, adapters)
    receipt_path = out_dir / ".conversion-run.json"
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    receipt["generated_tables"] = [
        {"path": "C:/Windows/System32/drivers/etc/hosts.ml.csv", "content_hash": "0" * 64}
    ]
    receipt_path.write_text(canonical_json(receipt), encoding="utf-8")

    report = _audit(adapters, tmp_path / "raw", out_dir)

    freshness = _checks(report)["conversion receipt and output freshness"]
    assert freshness.status.value == "FAIL"
    assert "not a safe relative" in freshness.detail


def test_a_non_object_conversion_receipt_is_reported_without_crashing(tmp_path: Path) -> None:
    """JSON syntax alone does not make a receipt an object with enforceable fields."""
    out_dir = tmp_path / "strict"
    adapters = {"fixture": _Adapter("fixture_output", ready=True)}
    _write_table(out_dir / "fixture_output_demo.ml.csv")
    out_dir.joinpath(".conversion-run.json").write_text("[]\n", encoding="utf-8")

    report = _audit(adapters, tmp_path / "raw", out_dir)

    freshness = _checks(report)["conversion receipt and output freshness"]
    assert freshness.status.value == "FAIL"
    assert "object" in freshness.detail


def test_a_receipt_table_hash_mismatch_fails_freshness_without_relying_on_validation(
    tmp_path: Path,
) -> None:
    """The receipt must bind the exact table hashes produced by that bulk invocation."""
    out_dir = tmp_path / "strict"
    adapters = {"fixture": _Adapter("fixture_output", ready=True)}
    _write_table(out_dir / "fixture_output_demo.ml.csv")
    _write_receipt(out_dir, adapters)
    receipt_path = out_dir / ".conversion-run.json"
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    receipt["generated_tables"][0]["content_hash"] = "0" * 64
    receipt_path.write_text(canonical_json(receipt), encoding="utf-8")

    report = _audit(adapters, tmp_path / "raw", out_dir)

    freshness = _checks(report)["conversion receipt and output freshness"]
    assert freshness.status.value == "FAIL"
    assert "content hash" in freshness.detail


def test_a_stale_gated_output_is_blocked_external_and_never_counts_as_current_coverage(
    tmp_path: Path,
) -> None:
    """A prior Addgene table cannot bypass today's missing authorized source files."""
    out_dir = tmp_path / "strict"
    adapter = _Adapter(
        "fixture_output",
        ready=False,
        detail="authorized input is absent",
        manual_steps="Download after accepting the rules",
    )
    adapters = {"gated": adapter}
    _write_table(out_dir / "fixture_output_demo.ml.csv")
    _write_receipt(out_dir, adapters, generated_paths=[], gated_source_keys=["gated"])

    report = _audit(adapters, tmp_path / "raw", out_dir)
    checks = _checks(report)

    assert checks["conversion receipt and output freshness"].status.value == "FAIL"
    coverage = checks["source coverage: gated"]
    assert coverage.status.value == "BLOCKED_EXTERNAL"
    assert "stale" in coverage.detail.lower()
    assert adapter.availability_calls == 1


def test_label_provenance_must_match_every_trusted_adapter_declaration(tmp_path: Path) -> None:
    """Self-consistent manifests cannot replace the adapter declarations trusted by the audit."""
    cases = {
        "source": "wrong_source",
        "source_url": "https://wrong.invalid",
        "licence": "wrong licence",
        "target_description": "wrong target",
        "label_semantics": {},
    }
    for field, replacement in cases.items():
        out_dir = tmp_path / field
        adapters = {"fixture": _Adapter("fixture_output", ready=True)}
        table = out_dir / "fixture_output_demo.ml.csv"
        _write_table(table)
        sidecar = table.with_suffix(".manifest.json")
        payload = json.loads(sidecar.read_text(encoding="utf-8"))
        payload["manifest"][field] = replacement
        sidecar.write_text(canonical_json(payload), encoding="utf-8")

        report = _audit(adapters, tmp_path / "raw", out_dir)
        provenance = _checks(report)["label provenance matches trusted adapter declarations"]

        assert provenance.status.value == "FAIL", field
        assert field in provenance.detail, field


def test_projection_name_index_and_slug_must_match_the_ordered_source_contract(
    tmp_path: Path,
) -> None:
    """Swapping two human names while retaining indexes silently relabels model outputs."""
    out_dir = tmp_path / "strict"
    adapters = {"fixture": _Adapter("fixture_output", ready=True)}
    first = out_dir / "fixture_output_demo__target_zero.ml.csv"
    second = out_dir / "fixture_output_demo__target_one.ml.csv"
    _write_projection(first, index=0, count=2)
    _write_projection(second, index=1, count=2)
    sidecar = first.with_suffix(".manifest.json")
    payload = json.loads(sidecar.read_text(encoding="utf-8"))
    payload["manifest"]["target_projection"]["original_target_name"] = "Target One"
    sidecar.write_text(canonical_json(payload), encoding="utf-8")

    report = _audit(adapters, tmp_path / "raw", out_dir)
    coverage = _checks(report)["source coverage: fixture"]

    assert coverage.status.value == "FAIL"
    assert "name/index" in coverage.detail


def test_projection_contracts_must_match_config_raw_inputs_samples_and_row_count(
    tmp_path: Path,
) -> None:
    """Mixing independently generated projections must not look like one coherent fan-out."""
    out_dir = tmp_path / "strict"
    adapters = {"fixture": _Adapter("fixture_output", ready=True)}
    first = out_dir / "fixture_output_demo__target_zero.ml.csv"
    second = out_dir / "fixture_output_demo__target_one.ml.csv"
    _write_projection(first, index=0, count=2)
    _write_projection(second, index=1, count=2)
    sidecar = second.with_suffix(".manifest.json")
    payload = json.loads(sidecar.read_text(encoding="utf-8"))
    payload["manifest"]["config_hash"] = "another-generation"
    payload["manifest"]["raw_input_hashes"] = {"raw.fa": "changed"}
    payload["manifest"]["samples_content_hash"] = "f" * 64
    payload["manifest"]["samples_n_rows"] = 3
    sidecar.write_text(canonical_json(payload), encoding="utf-8")

    report = _audit(adapters, tmp_path / "raw", out_dir)
    coverage = _checks(report)["source coverage: fixture"]

    assert coverage.status.value == "FAIL"
    assert "source-level contract" in coverage.detail


def test_a_projected_target_value_swap_fails_the_selected_target_hash(tmp_path: Path) -> None:
    """Valid binary values can still belong to the wrong target; column hashes bind meaning."""
    out_dir = tmp_path / "strict"
    adapters = {"fixture": _Adapter("fixture_output", ready=True)}
    first = out_dir / "fixture_output_demo__target_zero.ml.csv"
    second = out_dir / "fixture_output_demo__target_one.ml.csv"
    _write_projection(first, index=0, count=2)
    _write_projection(second, index=1, count=2)
    lines = first.read_text(encoding="utf-8").splitlines()
    first_row = lines[1].split(",")
    second_row = lines[2].split(",")
    first_row[1], second_row[1] = second_row[1], first_row[1]
    lines[1], lines[2] = ",".join(first_row), ",".join(second_row)
    first.write_text("\n".join(lines) + "\n", encoding="utf-8")
    sidecar = first.with_suffix(".manifest.json")
    payload = json.loads(sidecar.read_text(encoding="utf-8"))
    payload["manifest"]["content_hash"] = hash_file(first)
    sidecar.write_text(canonical_json(payload), encoding="utf-8")

    report = _audit(adapters, tmp_path / "raw", out_dir)
    coverage = _checks(report)["source coverage: fixture"]

    assert coverage.status.value == "FAIL"
    assert "selected target hash" in coverage.detail


def test_scientific_limitation_names_the_fifth_source_as_narrower_not_general(
    tmp_path: Path,
) -> None:
    """The fifth source corrects a narrower event question, not the original benchmark request."""
    adapters = {"fixture": _Adapter("fixture_output", ready=False)}
    report = _audit(adapters, tmp_path / "raw", tmp_path / "strict")
    check = _checks(report)["natural-versus-engineered target"]

    assert "four original sources" in check.detail
    assert "fifth" in check.detail
    assert "recorded-event" in check.detail
    assert "general benchmark" in check.detail


def test_json_and_markdown_render_the_same_status_counts() -> None:
    """Two separately computed reports could disagree while both looked plausible."""
    from nexus_genomics.compliance import (
        ComplianceCheck,
        ComplianceReport,
        ComplianceStatus,
        render_markdown,
    )

    report = ComplianceReport(
        (
            ComplianceCheck("pass one", ComplianceStatus.PASS, "established"),
            ComplianceCheck("pass two", ComplianceStatus.PASS, "established"),
            ComplianceCheck("blocked", ComplianceStatus.BLOCKED_EXTERNAL, "awaiting input"),
            ComplianceCheck("unsupported", ComplianceStatus.UNSUPPORTED_BY_SOURCE, "label absent"),
        )
    )
    payload = json.loads(report.to_json())
    markdown = render_markdown(report)

    assert payload["status_counts"] == {
        "PASS": 2,
        "BLOCKED_EXTERNAL": 1,
        "UNSUPPORTED_BY_SOURCE": 1,
        "FAIL": 0,
    }
    assert "ok" not in payload
    assert payload["has_failures"] is False
    assert payload["exit_ok"] is True
    assert payload["all_requirements_satisfied"] is False
    assert payload["overall_status"] == "INCOMPLETE"
    assert "**Result**: INCOMPLETE" in markdown
    for status, count in payload["status_counts"].items():
        assert f"| {status} | {count} |" in markdown


def test_compliance_cli_writes_reports_and_exits_zero_with_only_blockers(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """External and scientific blockers must stay visible without failing automation."""
    import nexus_genomics.cli as cli_module

    adapters = {
        "gated": _Adapter(
            "gated_output",
            ready=False,
            detail="account approval required",
            manual_steps="Download after approval",
        )
    }
    monkeypatch.setattr(cli_module, "ADAPTERS", adapters)
    reports_dir = tmp_path / "reports"
    strict_dir = tmp_path / "strict"
    config_path = tmp_path / "config.json"
    config_digest = _write_config(config_path, adapters)
    _write_receipt(
        strict_dir,
        adapters,
        config_hash=config_digest,
        gated_source_keys=["gated"],
    )

    result = runner.invoke(
        app,
        [
            "compliance",
            "--config",
            str(config_path),
            "--raw-dir",
            str(tmp_path / "raw"),
            "--out-dir",
            str(strict_dir),
            "--reports-dir",
            str(reports_dir),
        ],
    )

    assert result.exit_code == 0, result.output
    payload = json.loads((reports_dir / "email_compliance_report.json").read_text("utf-8"))
    markdown = (reports_dir / "email_compliance_report.md").read_text("utf-8")
    assert "ok" not in payload
    assert payload["exit_ok"] is True
    assert payload["all_requirements_satisfied"] is False
    assert payload["overall_status"] == "INCOMPLETE"
    assert payload["status_counts"]["FAIL"] == 0
    assert "BLOCKED_EXTERNAL" in result.output
    assert "UNSUPPORTED_BY_SOURCE" in markdown


def test_compliance_cli_exits_one_when_a_ready_source_has_no_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A real compliance defect must make the command usable as a CI gate."""
    import nexus_genomics.cli as cli_module

    adapters = {"ready": _Adapter("ready_output", ready=True)}
    monkeypatch.setattr(cli_module, "ADAPTERS", adapters)
    reports_dir = tmp_path / "reports"
    strict_dir = tmp_path / "strict"
    config_path = tmp_path / "config.json"
    config_digest = _write_config(config_path, adapters)
    _write_receipt(strict_dir, adapters, config_hash=config_digest)

    result = runner.invoke(
        app,
        [
            "compliance",
            "--config",
            str(config_path),
            "--raw-dir",
            str(tmp_path / "raw"),
            "--out-dir",
            str(strict_dir),
            "--reports-dir",
            str(reports_dir),
        ],
    )

    assert result.exit_code == 1
    payload = cast(
        dict[str, Any],
        json.loads((reports_dir / "email_compliance_report.json").read_text("utf-8")),
    )
    assert payload["has_failures"] is True
    assert payload["exit_ok"] is False
    assert payload["overall_status"] == "FAIL"
    assert payload["status_counts"]["FAIL"] > 0


def test_compliance_rejects_a_misspelled_profile_before_availability_checks(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Auditing `_ful` defaults would certify artifacts outside the configured profiles."""
    import nexus_genomics.cli as cli_module

    adapter = _Adapter("fixture_output", ready=True)
    monkeypatch.setattr(cli_module, "ADAPTERS", {"fixture": adapter})

    result = runner.invoke(app, ["compliance", "--profile", "ful"])

    assert result.exit_code == 2
    assert adapter.availability_calls == 0


def test_compliance_requires_the_selected_profile_block_before_availability_checks(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A valid profile name still needs explicit source settings in the audited config."""
    import nexus_genomics.cli as cli_module

    adapter = _Adapter("fixture_output", ready=True)
    adapters = {"fixture": adapter}
    monkeypatch.setattr(cli_module, "ADAPTERS", adapters)
    config_path = tmp_path / "config.json"
    digest = _write_config(config_path, adapters)
    parsed = json.loads(config_path.read_text(encoding="utf-8"))
    del parsed["sources"]["fixture"]["full"]
    config_path.write_text(json.dumps(parsed), encoding="utf-8")
    _write_receipt(tmp_path / "strict", adapters, config_hash=digest)

    result = runner.invoke(
        app,
        [
            "compliance",
            "--profile",
            "full",
            "--config",
            str(config_path),
            "--out-dir",
            str(tmp_path / "strict"),
        ],
    )

    assert result.exit_code == 1
    assert adapter.availability_calls == 0
