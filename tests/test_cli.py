"""Focused CLI contracts for command discovery, orchestration, reports, and exit codes.

The tests use offline synthetic adapters and temporary configuration. They establish the
documented command names and selected option behavior; they do not claim to execute every
literal README or runbook command against the live public sources.
"""

from __future__ import annotations

import csv
import json
import subprocess
from collections.abc import Mapping
from pathlib import Path
from typing import Any, cast

import pytest
import yaml
from typer.testing import CliRunner

from nexus_genomics.adapters.base import LoadedSource
from nexus_genomics.cleaning import Alphabet
from nexus_genomics.cli import app
from nexus_genomics.common import SampleRecord
from nexus_genomics.pipeline import ConvertOptions, convert

runner = CliRunner()


class _TwoTargetAdapter:
    name = "fake_two_target"
    source_url = "https://example.invalid/two-target"
    licence = "Synthetic test fixture; not a distributable dataset."
    target_description = "Two documented binary scientific targets from a synthetic fixture."

    def __init__(self) -> None:
        self.load_calls = 0
        self.loaded: LoadedSource | None = None

    def load(self, raw_dir: Path, options: Mapping[str, Any]) -> LoadedSource:
        _ = raw_dir, options
        self.load_calls += 1
        self.loaded = LoadedSource(
            records=[
                SampleRecord("s1", "ACGT", (1, 0), "fixture", {}),
                SampleRecord("s2", "TGCA", (0, 1), "fixture", {}),
            ],
            n_targets=2,
            label_semantics={0: "absent", 1: "present"},
            alphabet=Alphabet.DNA_NUCLEOTIDE,
            target_names=("First mechanism", "Second mechanism"),
            manifest_extra={"adapter_metadata": "preserved"},
        )
        return self.loaded


def _configure_fake_source(config: Path, *, target_layout: str) -> None:
    parsed = yaml.safe_load(config.read_text(encoding="utf-8"))
    parsed["output"] = {"target_layout": target_layout}
    parsed["sources"]["fake"] = {"enabled": True, "demo": {}}
    config.write_text(yaml.safe_dump(parsed), encoding="utf-8")


def _configure_only_fake_source(config: Path, *, target_layout: str) -> None:
    parsed = yaml.safe_load(config.read_text(encoding="utf-8"))
    parsed["output"] = {"target_layout": target_layout}
    for source in parsed["sources"].values():
        source["enabled"] = False
    parsed["sources"]["fake"] = {"enabled": True, "demo": {}, "full": {}}
    config.write_text(yaml.safe_dump(parsed), encoding="utf-8")


def _read_table(path: Path) -> list[list[str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.reader(stream))


def _read_manifest(path: Path) -> dict[str, Any]:
    payload = json.loads(path.with_suffix(".manifest.json").read_text(encoding="utf-8"))
    return cast(dict[str, Any], payload["manifest"])


@pytest.fixture
def config(tmp_path: Path) -> Path:
    path = tmp_path / "c.yaml"
    path.write_text(
        yaml.safe_dump(
            {
                "encoding": {"mode": "letter"},
                "output": {"target_layout": "native"},
                "window": {
                    "length": 12,
                    "policy": "pad_or_truncate",
                    "stride": 12,
                    "max_windows_per_sample": 4,
                },
                "validation": {"round_trip_sample_size": 4, "round_trip_seed": 0},
                "sources": {
                    "codontransformer": {"enabled": False},
                    "seqscreen": {"enabled": False},
                    "felix": {"enabled": False},
                    "addgene": {"enabled": True, "demo": {"max_records": 5}},
                },
            }
        ),
        encoding="utf-8",
    )
    return path


@pytest.fixture
def table(tmp_path: Path) -> Path:
    out = tmp_path / "t.ml.csv"
    convert(
        LoadedSource(
            records=[
                SampleRecord("s1", "ACGT", (0,), "t", {}),
                SampleRecord("s2", "TTTTAA", (1,), "t", {}),
            ],
            n_targets=1,
            label_semantics={0: "natural", 1: "engineered"},
            alphabet=Alphabet.DNA_NUCLEOTIDE,
        ),
        out,
        ConvertOptions(width=12),
        source_name="t",
        source_url="https://example.invalid",
        licence="t",
        target_description="t",
        notes="t",
        config_digest="d",
    )
    return out


def test_every_documented_command_exists() -> None:
    """A renamed command silently breaks the README and the demo script."""
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
    for command in (
        "audit",
        "compliance",
        "convert",
        "convert-all",
        "validate",
        "reproduce",
        "formats",
    ):
        assert command in result.output, command


def test_audit_reports_every_source_and_writes_a_file(tmp_path: Path) -> None:
    """The first thing the demo runs. It must work with an empty data directory."""
    out = tmp_path / "audit.md"
    result = runner.invoke(app, ["audit", "--raw-dir", str(tmp_path / "raw"), "--out", str(out)])
    assert result.exit_code == 0, result.output
    text = out.read_text(encoding="utf-8")
    for source in ("codontransformer", "seqscreen", "felix", "addgene"):
        assert f"`{source}`" in text


def test_audit_prints_the_gated_download_steps(tmp_path: Path) -> None:
    """If this text goes missing, a blocked source looks like a broken one."""
    out = tmp_path / "audit.md"
    runner.invoke(app, ["audit", "--raw-dir", str(tmp_path / "raw"), "--out", str(out)])
    text = out.read_text(encoding="utf-8")
    assert "drivendata.org/accounts/signup" in text
    assert "pubs.acs.org" in text


def test_validate_passes_a_good_file_and_writes_both_reports(
    table: Path, config: Path, tmp_path: Path
) -> None:
    reports = tmp_path / "reports"
    result = runner.invoke(
        app, ["validate", str(table), "--config", str(config), "--reports-dir", str(reports)]
    )
    assert result.exit_code == 0, result.output
    payload = json.loads((reports / "validation_report.json").read_text(encoding="utf-8"))
    assert payload["all_pass"] is True
    assert payload["files"][0]["n_failed"] == 0
    markdown = (reports / "validation_report.md").read_text(encoding="utf-8")
    assert "1 of 1 file(s) pass every check." in markdown


def test_validate_exits_non_zero_on_a_broken_file(
    table: Path, config: Path, tmp_path: Path
) -> None:
    """A validator that always exits 0 cannot be used in CI or trusted in a demo."""
    lines = table.read_text(encoding="utf-8").splitlines()
    lines.append(lines[1])  # duplicate primary key
    table.write_text("\n".join(lines) + "\n", encoding="utf-8")
    result = runner.invoke(
        app,
        [
            "validate",
            str(table),
            "--config",
            str(config),
            "--reports-dir",
            str(tmp_path / "r"),
        ],
    )
    assert result.exit_code == 1


def test_reproduce_compares_two_runs(table: Path, tmp_path: Path) -> None:
    twin = tmp_path / "twin.ml.csv"
    twin.write_bytes(table.read_bytes())
    assert runner.invoke(app, ["reproduce", str(table), str(twin)]).exit_code == 0

    twin.write_bytes(table.read_bytes() + b"extra\n")
    assert runner.invoke(app, ["reproduce", str(table), str(twin)]).exit_code == 1


def test_convert_all_reports_a_gated_source_instead_of_aborting(
    config: Path, tmp_path: Path
) -> None:
    """One blocked source must not stop the other three from being built.

    This is the behaviour the Monday demo depends on: `convert-all` has to finish.
    """
    result = runner.invoke(
        app,
        [
            "convert-all",
            "--config",
            str(config),
            "--raw-dir",
            str(tmp_path / "raw"),
            "--out-dir",
            str(tmp_path / "out"),
        ],
    )
    assert result.exit_code == 0, result.output
    assert "addgene" in result.output
    assert "skipping codontransformer" in result.output
    receipt = json.loads((tmp_path / "out" / ".conversion-run.json").read_text("utf-8"))
    assert receipt["complete"] is True
    assert receipt["successfully_converted_source_keys"] == []
    assert receipt["externally_gated_source_keys"] == ["addgene"]
    assert receipt["generated_tables"] == []


def test_convert_all_exits_non_zero_on_a_failure_that_is_not_a_gated_download(
    config: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A crash must not look like a source waiting on a manual download.

    Both exiting 0 means a scheduled run that ran out of memory goes green, and the next
    step consumes the previous run's stale outputs believing they are fresh.
    """
    import nexus_genomics.cli as cli_module

    def boom(*_args: object, **_kwargs: object) -> None:
        raise MemoryError("out of memory")

    monkeypatch.setattr(cli_module, "convert_source", boom)
    result = runner.invoke(
        app,
        [
            "convert-all",
            "--config",
            str(config),
            "--raw-dir",
            str(tmp_path / "raw"),
            "--out-dir",
            str(tmp_path / "out"),
        ],
    )
    assert result.exit_code == 1
    assert "MemoryError" in result.output


def test_convert_all_treats_an_unexpected_file_not_found_as_a_real_failure(
    config: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Missing open-source or output files must not masquerade as an authorized-data gate."""
    import nexus_genomics.cli as cli_module

    def missing_open_file(*_args: object, **_kwargs: object) -> list[Path]:
        raise FileNotFoundError("open-source extraction vanished")

    monkeypatch.setattr(cli_module, "convert_source", missing_open_file)
    result = runner.invoke(
        app,
        [
            "convert-all",
            "--config",
            str(config),
            "--raw-dir",
            str(tmp_path / "raw"),
            "--out-dir",
            str(tmp_path / "out"),
        ],
    )

    assert result.exit_code == 1
    assert "FileNotFoundError" in result.output
    assert "gated download" not in result.output
    receipt = json.loads((tmp_path / "out" / ".conversion-run.json").read_text("utf-8"))
    assert receipt["complete"] is False


def test_convert_all_writes_a_complete_receipt_bound_to_exact_tables_and_hashes(
    config: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A terminal log cannot prove which files survived; the atomic receipt must do so."""
    import nexus_genomics.cli as cli_module

    adapter = _TwoTargetAdapter()
    _configure_only_fake_source(config, target_layout="native")
    monkeypatch.setattr(cli_module, "ADAPTERS", {"fake": adapter})
    monkeypatch.setattr(cli_module, "get_adapter", lambda _name: adapter)
    out_dir = tmp_path / "strict"

    result = runner.invoke(
        app,
        [
            "convert-all",
            "--config",
            str(config),
            "--out-dir",
            str(out_dir),
            "--strict-single-target",
        ],
    )

    assert result.exit_code == 0, result.output
    receipt = json.loads((out_dir / ".conversion-run.json").read_text("utf-8"))
    assert receipt["format"] == "nexus-genomics-conversion-run"
    assert receipt["version"] == 1
    assert receipt["generator"] == "nexus-genomics 0.1.0"
    assert receipt["complete"] is True
    assert receipt["profile"] == "demo"
    assert receipt["target_layout"] == "single"
    assert receipt["enabled_source_keys"] == ["fake"]
    assert receipt["successfully_converted_source_keys"] == ["fake"]
    assert receipt["externally_gated_source_keys"] == []
    assert [entry["path"] for entry in receipt["generated_tables"]] == [
        "fake_two_target_demo__first_mechanism.ml.csv",
        "fake_two_target_demo__second_mechanism.ml.csv",
    ]
    for entry in receipt["generated_tables"]:
        table = out_dir / entry["path"]
        manifest = _read_manifest(table)
        assert entry["content_hash"] == manifest["content_hash"]


def test_a_failed_rebuild_atomically_invalidates_a_prior_complete_receipt(
    config: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Coherent old triplets may remain, but no failed run may leave them certified current."""
    import nexus_genomics.cli as cli_module

    adapter = _TwoTargetAdapter()
    _configure_only_fake_source(config, target_layout="native")
    monkeypatch.setattr(cli_module, "ADAPTERS", {"fake": adapter})
    monkeypatch.setattr(cli_module, "get_adapter", lambda _name: adapter)
    out_dir = tmp_path / "strict"
    command = [
        "convert-all",
        "--config",
        str(config),
        "--out-dir",
        str(out_dir),
        "--strict-single-target",
    ]
    first = runner.invoke(app, command)
    assert first.exit_code == 0, first.output
    assert json.loads((out_dir / ".conversion-run.json").read_text("utf-8"))["complete"]

    def fail_rebuild(*_args: object, **_kwargs: object) -> list[Path]:
        raise RuntimeError("synthetic rebuild failed")

    monkeypatch.setattr(cli_module, "convert_source", fail_rebuild)
    second = runner.invoke(app, command)

    assert second.exit_code == 1
    receipt = json.loads((out_dir / ".conversion-run.json").read_text("utf-8"))
    assert receipt["complete"] is False
    assert (out_dir / "fake_two_target_demo__first_mechanism.ml.csv").exists()


def test_convert_native_mode_writes_one_table_with_every_target(
    config: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Fan-out in native mode would change established files and multilabel table behavior."""
    import nexus_genomics.cli as cli_module

    adapter = _TwoTargetAdapter()
    _configure_fake_source(config, target_layout="native")
    monkeypatch.setattr(cli_module, "get_adapter", lambda _name: adapter)
    out_dir = tmp_path / "out"

    result = runner.invoke(
        app,
        ["convert", "fake", "--config", str(config), "--out-dir", str(out_dir)],
    )

    assert result.exit_code == 0, result.output
    assert "target_layout" in result.output
    assert "native" in result.output
    assert "2 rows x 15 columns" in result.output
    path = out_dir / "fake_two_target_demo.ml.csv"
    rows = _read_table(path)
    assert rows[0][:3] == ["sample_id", "target_0", "target_1"]
    assert [row[1:3] for row in rows[1:]] == [["1", "0"], ["0", "1"]]
    assert _read_manifest(path)["target_layout"] == "native"
    assert adapter.load_calls == 1


def test_convert_strict_mode_writes_one_named_table_per_target(
    config: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Wrong projection indexes or reused paths silently swap or overwrite target labels."""
    import nexus_genomics.cli as cli_module

    adapter = _TwoTargetAdapter()
    _configure_fake_source(config, target_layout="native")
    monkeypatch.setattr(cli_module, "get_adapter", lambda _name: adapter)
    out_dir = tmp_path / "out"

    result = runner.invoke(
        app,
        [
            "convert",
            "fake",
            "--config",
            str(config),
            "--out-dir",
            str(out_dir),
            "--strict-single-target",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "single" in result.output
    assert result.output.count("2 rows x 14 columns") == 2
    expected = {
        "first_mechanism": (["1", "0"], 0, "First mechanism"),
        "second_mechanism": (["0", "1"], 1, "Second mechanism"),
    }
    for slug, (values, index, name) in expected.items():
        path = out_dir / f"fake_two_target_demo__{slug}.ml.csv"
        rows = _read_table(path)
        assert rows[0][:2] == ["sample_id", "target"]
        assert [row[1] for row in rows[1:]] == values
        manifest = _read_manifest(path)
        assert manifest["target_layout"] == "single"
        assert manifest["adapter_metadata"] == "preserved"
        projection = manifest["target_projection"]
        assert projection["original_target_count"] == 2
        assert projection["original_target_index"] == index
        assert projection["original_target_name"] == name
        assert projection["original_target_slug"] == slug
        assert projection["source_target_names"] == ["First mechanism", "Second mechanism"]
        assert projection["source_target_slugs"] == ["first_mechanism", "second_mechanism"]
        assert len(projection["source_target_hashes"]) == 2
        assert projection["source_target_hash_algorithm"] == "blake2b-256"
        assert projection["source_row_count"] == 2
        assert adapter.target_description in manifest["target_description"]
        assert name in manifest["target_description"]
    assert adapter.load_calls == 1
    assert adapter.loaded is not None
    assert adapter.loaded.manifest_extra == {"adapter_metadata": "preserved"}


def test_configured_single_layout_is_used_when_the_convert_flag_is_omitted(
    config: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Collapsing an omitted flag to false would ignore the byte-affecting config default."""
    import nexus_genomics.cli as cli_module

    adapter = _TwoTargetAdapter()
    _configure_fake_source(config, target_layout="single")
    monkeypatch.setattr(cli_module, "get_adapter", lambda _name: adapter)
    out_dir = tmp_path / "out"

    result = runner.invoke(
        app,
        ["convert", "fake", "--config", str(config), "--out-dir", str(out_dir)],
    )

    assert result.exit_code == 0, result.output
    assert (out_dir / "fake_two_target_demo__first_mechanism.ml.csv").is_file()
    assert (out_dir / "fake_two_target_demo__second_mechanism.ml.csv").is_file()


def test_native_targets_flag_overrides_a_configured_single_layout(
    config: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """An explicit native request must win over the configured strict export default."""
    import nexus_genomics.cli as cli_module

    adapter = _TwoTargetAdapter()
    _configure_fake_source(config, target_layout="single")
    monkeypatch.setattr(cli_module, "get_adapter", lambda _name: adapter)
    out_dir = tmp_path / "out"

    result = runner.invoke(
        app,
        [
            "convert",
            "fake",
            "--config",
            str(config),
            "--out-dir",
            str(out_dir),
            "--native-targets",
        ],
    )

    assert result.exit_code == 0, result.output
    assert (out_dir / "fake_two_target_demo.ml.csv").is_file()
    assert not (out_dir / "fake_two_target_demo__first_mechanism.ml.csv").exists()


def test_convert_refuses_an_unknown_configured_target_layout(
    config: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Guessing an unknown layout could emit an unexpected number of training tables."""
    import nexus_genomics.cli as cli_module

    adapter = _TwoTargetAdapter()
    _configure_fake_source(config, target_layout="wide")
    monkeypatch.setattr(cli_module, "get_adapter", lambda _name: adapter)

    result = runner.invoke(
        app,
        ["convert", "fake", "--config", str(config), "--out-dir", str(tmp_path / "out")],
    )

    assert result.exit_code == 1
    assert isinstance(result.exception, ValueError)
    assert "native" in str(result.exception)
    assert "single" in str(result.exception)
    assert adapter.load_calls == 0


def test_convert_refuses_a_boolean_window_cap_from_yaml_before_source_io(
    config: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """YAML parses ``true`` as bool, which must not silently become a one-window cap."""
    import nexus_genomics.cli as cli_module

    adapter = _TwoTargetAdapter()
    _configure_fake_source(config, target_layout="native")
    parsed = yaml.safe_load(config.read_text(encoding="utf-8"))
    parsed["window"]["max_windows_per_sample"] = True
    config.write_text(yaml.safe_dump(parsed), encoding="utf-8")
    monkeypatch.setattr(cli_module, "get_adapter", lambda _name: adapter)

    result = runner.invoke(
        app,
        ["convert", "fake", "--config", str(config), "--out-dir", str(tmp_path / "out")],
    )

    assert result.exit_code == 1
    assert isinstance(result.exception, TypeError)
    assert "non-boolean integer" in str(result.exception)
    assert adapter.load_calls == 0


def test_convert_rejects_a_misspelled_profile_before_source_io(
    config: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """``ful`` must not run with empty adapter defaults and create plausible `_ful` files."""
    import nexus_genomics.cli as cli_module

    adapter = _TwoTargetAdapter()
    _configure_fake_source(config, target_layout="native")
    monkeypatch.setattr(cli_module, "get_adapter", lambda _name: adapter)

    result = runner.invoke(
        app,
        [
            "convert",
            "fake",
            "--profile",
            "ful",
            "--config",
            str(config),
            "--out-dir",
            str(tmp_path / "out"),
        ],
    )

    assert result.exit_code == 2
    assert adapter.load_calls == 0


def test_convert_all_rejects_a_misspelled_profile_before_any_source_io(
    config: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Bulk mode must share the same closed profile vocabulary as single-source conversion."""
    import nexus_genomics.cli as cli_module

    calls = 0

    def record_call(*_args: object, **_kwargs: object) -> list[Path]:
        nonlocal calls
        calls += 1
        return []

    monkeypatch.setattr(cli_module, "convert_source", record_call)
    result = runner.invoke(
        app,
        [
            "convert-all",
            "--profile",
            "ful",
            "--config",
            str(config),
            "--out-dir",
            str(tmp_path / "out"),
        ],
    )

    assert result.exit_code == 2
    assert calls == 0


def test_convert_requires_the_selected_profile_block_before_source_io(
    config: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A valid profile name without source settings must not fall back to adapter defaults."""
    import nexus_genomics.cli as cli_module

    adapter = _TwoTargetAdapter()
    _configure_fake_source(config, target_layout="native")
    monkeypatch.setattr(cli_module, "get_adapter", lambda _name: adapter)

    result = runner.invoke(
        app,
        [
            "convert",
            "fake",
            "--profile",
            "full",
            "--config",
            str(config),
            "--out-dir",
            str(tmp_path / "out"),
        ],
    )

    assert result.exit_code == 1
    assert isinstance(result.exception, ValueError)
    assert "has no 'full' profile block" in str(result.exception)
    assert adapter.load_calls == 0


def test_convert_all_forwards_the_resolved_strict_layout(
    config: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Dropping the mode in convert-all would make bulk output differ from convert output."""
    import nexus_genomics.cli as cli_module

    calls: list[bool | None] = []

    def record_call(
        source: str,
        config_path: Path,
        profile: str,
        raw_dir: Path,
        out_dir: Path,
        strict_single_target: bool | None,
    ) -> list[Path]:
        _ = source, config_path, profile, raw_dir, out_dir
        calls.append(strict_single_target)
        return []

    monkeypatch.setattr(cli_module, "convert_source", record_call)
    result = runner.invoke(
        app,
        [
            "convert-all",
            "--config",
            str(config),
            "--raw-dir",
            str(tmp_path / "raw"),
            "--out-dir",
            str(tmp_path / "out"),
            "--strict-single-target",
        ],
    )

    assert result.exit_code == 0, result.output
    assert calls == [True]


def test_formats_states_that_nexus_compatibility_is_not_claimed() -> None:
    """The standing rule inherited from qecgen. It must be visible from the CLI."""
    result = runner.invoke(app, ["formats"])
    assert result.exit_code == 0
    assert "not claimed" in result.output
    assert "CATEGORICAL TOKENS" in result.output


def test_the_shipped_config_parses_and_names_every_registered_source() -> None:
    """config/default.yaml drifting from ADAPTERS would make convert-all silently skip one."""
    from nexus_genomics.adapters import ADAPTERS

    parsed = yaml.safe_load(Path("config/default.yaml").read_text(encoding="utf-8"))
    assert set(parsed["sources"]) == set(ADAPTERS)
    for name, entry in parsed["sources"].items():
        assert "demo" in entry, name
        assert entry["enabled"] is True, name
    assert parsed["window"]["length"] == 1000
    assert parsed["encoding"]["mode"] == "letter"
    assert parsed["output"]["target_layout"] == "native"


def test_every_default_or_documented_generated_report_destination_is_git_ignored() -> None:
    """Source-derived class counts and gated-input hashes must not become commit candidates."""
    generated = {
        "reports/validation_report.json",
        "reports/validation_report.md",
        "reports/strict/validation_report.json",
        "reports/strict/validation_report.md",
        "reports/email_compliance_report.json",
        "reports/email_compliance_report.md",
        "reports/source_audit_status.md",
    }
    result = subprocess.run(
        [
            "git",
            "-c",
            "core.excludesFile=NUL",
            "check-ignore",
            "--no-index",
            *sorted(generated),
        ],
        cwd=Path.cwd(),
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert set(result.stdout.splitlines()) == generated
