"""Executable evidence for every requirement in the email compliance report."""

from __future__ import annotations

import json
from collections.abc import Mapping
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

import pytest
from typer.testing import CliRunner

from nexus_genomics.adapters.base import Availability, LoadedSource
from nexus_genomics.cleaning import Alphabet
from nexus_genomics.cli import app
from nexus_genomics.common import SampleRecord, canonical_json, hash_file
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

    def availability(self, raw_dir: Path) -> Availability:
        _ = raw_dir
        return self._availability

    def load(self, raw_dir: Path, options: Mapping[str, Any]) -> LoadedSource:
        _ = raw_dir, options
        raise AssertionError("the compliance audit must not load source data")


def _write_table(path: Path, *, n_targets: int = 1) -> None:
    labels = [(0,), (1,)] if n_targets == 1 else [(1, 0), (0, 1)]
    convert(
        LoadedSource(
            records=[
                SampleRecord("s1", "ACGT", labels[0], "fixture", {}),
                SampleRecord("s2", "TGCA", labels[1], "fixture", {}),
            ],
            n_targets=n_targets,
            label_semantics={0: "absent", 1: "present"},
            alphabet=Alphabet.DNA_NUCLEOTIDE,
        ),
        path,
        ConvertOptions(width=4),
        source_name="fixture",
        source_url="https://example.invalid/compliance",
        licence="Synthetic fixture",
        target_description="Synthetic target",
        notes="Synthetic compliance fixture",
        config_digest="fixture",
    )


def _audit(
    adapters: Mapping[str, _Adapter], raw_dir: Path, out_dir: Path, profile: str = "demo"
) -> ComplianceReport:
    from nexus_genomics.compliance import audit_email_compliance

    return audit_email_compliance(adapters, raw_dir, out_dir, profile)


def _checks(report: ComplianceReport) -> dict[str, Any]:
    return {check.requirement: check for check in report.checks}


def test_compliance_models_expose_all_four_statuses_and_fail_only_for_fail() -> None:
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
    assert report.ok is False
    assert report.status_counts == {
        "PASS": 1,
        "BLOCKED_EXTERNAL": 1,
        "UNSUPPORTED_BY_SOURCE": 1,
        "FAIL": 1,
    }
    blockers_only = ComplianceReport(checks[1:3])
    assert blockers_only.ok is True


def test_audit_discovers_only_exact_native_and_projected_profile_names(tmp_path: Path) -> None:
    """A loose glob could grade stale profiles or unrelated tables as current strict output."""
    out_dir = tmp_path / "strict"
    _write_table(out_dir / "fixture_output_demo.ml.csv")
    _write_table(out_dir / "fixture_output_demo__selected_target.ml.csv")
    out_dir.joinpath("fixture_output_full.ml.csv").write_text("broken", encoding="utf-8")
    out_dir.joinpath("fixture_output_demo_extra.ml.csv").write_text("broken", encoding="utf-8")
    out_dir.joinpath("unrelated_demo.ml.csv").write_text("broken", encoding="utf-8")

    report = _audit({"fixture": _Adapter("fixture_output", ready=True)}, tmp_path / "raw", out_dir)

    assert report.ok is True
    structural = report.checks[:5]
    assert all(check.status.value == "PASS" for check in structural)
    assert all("2 file(s)" in check.detail for check in structural)


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
    assert report.ok is False


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
    assert report.ok is False


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
    for status, count in payload["status_counts"].items():
        assert f"| {status} | {count} |" in markdown


def test_compliance_cli_writes_reports_and_exits_zero_with_only_blockers(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """External and scientific blockers must stay visible without failing automation."""
    import nexus_genomics.cli as cli_module

    monkeypatch.setattr(
        cli_module,
        "ADAPTERS",
        {
            "gated": _Adapter(
                "gated_output",
                ready=False,
                detail="account approval required",
                manual_steps="Download after approval",
            )
        },
    )
    reports_dir = tmp_path / "reports"

    result = runner.invoke(
        app,
        [
            "compliance",
            "--raw-dir",
            str(tmp_path / "raw"),
            "--out-dir",
            str(tmp_path / "strict"),
            "--reports-dir",
            str(reports_dir),
        ],
    )

    assert result.exit_code == 0, result.output
    payload = json.loads((reports_dir / "email_compliance_report.json").read_text("utf-8"))
    markdown = (reports_dir / "email_compliance_report.md").read_text("utf-8")
    assert payload["ok"] is True
    assert payload["status_counts"]["FAIL"] == 0
    assert "BLOCKED_EXTERNAL" in result.output
    assert "UNSUPPORTED_BY_SOURCE" in markdown


def test_compliance_cli_exits_one_when_a_ready_source_has_no_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A real compliance defect must make the command usable as a CI gate."""
    import nexus_genomics.cli as cli_module

    monkeypatch.setattr(cli_module, "ADAPTERS", {"ready": _Adapter("ready_output", ready=True)})
    reports_dir = tmp_path / "reports"

    result = runner.invoke(
        app,
        [
            "compliance",
            "--raw-dir",
            str(tmp_path / "raw"),
            "--out-dir",
            str(tmp_path / "strict"),
            "--reports-dir",
            str(reports_dir),
        ],
    )

    assert result.exit_code == 1
    payload = cast(
        dict[str, Any],
        json.loads((reports_dir / "email_compliance_report.json").read_text("utf-8")),
    )
    assert payload["ok"] is False
    assert payload["status_counts"]["FAIL"] > 0
