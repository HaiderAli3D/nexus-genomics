"""The command line, exercised the way the README tells someone to use it.

Every command documented in README.md and reports/monday_demo_checklist.md is invoked here,
so a renamed command or a changed flag breaks the suite rather than the demo.
"""

from __future__ import annotations

import csv
import json
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
        assert manifest["target_projection"] == {
            "original_target_count": 2,
            "original_target_index": index,
            "original_target_name": name,
        }
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
    ) -> None:
        _ = source, config_path, profile, raw_dir, out_dir
        calls.append(strict_single_target)

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
