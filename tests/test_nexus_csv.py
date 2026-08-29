"""The writer's conventions, and the traps they exist to close."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from nexus_genomics.nexus_csv import (
    COLUMN_ORDER,
    FORMAT,
    NotANexusGenomicsDatasetError,
    Row,
    manifest_path,
    position_columns,
    read_nexus_csv,
    samples_path,
    target_columns,
    write_nexus_csv,
)


def write(tmp_path: Path, rows: list[Row], n_features: int, n_targets: int = 1) -> Path:
    path = tmp_path / "d.ml.csv"
    write_nexus_csv(
        path,
        rows,
        n_features=n_features,
        n_targets=n_targets,
        manifest={"source": "test", "label_semantics": {"0": "a", "1": "b"}},
        sample_metadata=[{"sample_id": r.sample_id, "row_sequence_blake2b": "x"} for r in rows],
    )
    return path


def test_column_names_are_zero_padded_to_this_files_own_width() -> None:
    """`sorted(df.columns)` must not permute the feature matrix.

    Unpadded, a column sort -- inside a join, a concat, or a feature-store schema -- puts
    position_10 before position_2 and silently reorders every sequence. The model still
    trains and converges; it just means nothing.
    """
    assert position_columns(9) == [f"position_{i}" for i in range(1, 10)]
    assert position_columns(10)[:2] == ["position_01", "position_02"]
    assert position_columns(1000)[0] == "position_0001"
    assert position_columns(1000)[-1] == "position_1000"


def test_sorting_the_columns_cannot_permute_the_feature_matrix() -> None:
    """The property the padding exists for, asserted directly."""
    names = position_columns(1000)
    assert sorted(names) == names


def test_position_columns_are_one_based() -> None:
    """Sequence coordinates are one-based in biology; a position_0 invites an off-by-one."""
    assert position_columns(5)[0] == "position_1"
    assert "position_0" not in position_columns(5)


def test_multilabel_targets_are_padded_too() -> None:
    """32 FunSoC targets unpadded would sort target_10 before target_2.

    qecgen leaves its multi-target spelling unpadded, but it only ever writes one
    observable, so that path was never exercised past nine. Here it is exercised at 32.
    """
    assert target_columns(1) == ["target"]
    assert target_columns(32)[:2] == ["target_00", "target_01"]
    assert sorted(target_columns(32)) == target_columns(32)


def test_the_layout_is_key_then_target_then_variables(tmp_path: Path) -> None:
    """The consumer asked for `Primary Key, Target, Variable_1, ...`."""
    assert COLUMN_ORDER == ("index", "target", "feature")
    path = write(tmp_path, [Row("s1", (1,), (1, 3, 7))], 3)
    header = path.read_text(encoding="utf-8").splitlines()[0].split(",")
    assert header == ["sample_id", "target", "position_1", "position_2", "position_3"]


def test_the_table_has_one_header_row_and_no_comment_lines(tmp_path: Path) -> None:
    """The whole point of the format: `pandas.read_csv(path)` with no arguments works.

    A `#` line here would make that call return a single garbage column.
    """
    path = write(tmp_path, [Row("s1", (0,), (1, 3, 7))], 3)
    text = path.read_text(encoding="utf-8")
    assert not any(line.startswith("#") for line in text.splitlines())
    assert text.endswith("\n")
    assert "\r" not in text


def test_the_sidecar_names_every_column_so_no_consumer_builds_one(tmp_path: Path) -> None:
    """The pad width is a property of the file, so the names must be published, not described."""
    path = write(tmp_path, [Row("s1", (0,), (1, 3, 7))], 3)
    payload = json.loads(manifest_path(path).read_text(encoding="utf-8"))
    assert payload["format"] == FORMAT
    assert payload["columns"]["feature_columns"] == ["position_1", "position_2", "position_3"]
    assert payload["columns"]["index_columns"] == ["sample_id"]
    assert payload["columns"]["target_columns"] == ["target"]


def test_a_missing_sidecar_is_not_a_nexus_genomics_dataset(tmp_path: Path) -> None:
    """A bare .ml.csv is shape-identical to any other one-header-row CSV.

    The named sidecar is what makes a file provably ours rather than merely plausible, and
    the group is committed atomically so a table can never legitimately appear without it.
    """
    path = write(tmp_path, [Row("s1", (0,), (1, 3, 7))], 3)
    manifest_path(path).unlink()
    with pytest.raises(NotANexusGenomicsDatasetError, match="plain CSV"):
        read_nexus_csv(path)


def test_a_future_sidecar_version_is_refused_rather_than_half_understood(tmp_path: Path) -> None:
    """A later revision may add a required key; reading it as v1 would drop it silently."""
    path = write(tmp_path, [Row("s1", (0,), (1, 3, 7))], 3)
    side = manifest_path(path)
    payload = json.loads(side.read_text(encoding="utf-8"))
    payload["version"] = 2
    side.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="understands only 1"):
        read_nexus_csv(path)


def test_duplicate_primary_keys_are_refused(tmp_path: Path) -> None:
    """A duplicate key fans out on every join and silently reweights the dataset."""
    with pytest.raises(ValueError, match="duplicated id"):
        write(tmp_path, [Row("s1", (0,), (1,)), Row("s1", (1,), (3,))], 1)


def test_a_row_of_the_wrong_width_is_refused_rather_than_padded(tmp_path: Path) -> None:
    """Padding a short row here would fabricate residues that were never in the source."""
    with pytest.raises(ValueError, match="fabricate residues"):
        write(tmp_path, [Row("s1", (0,), (1, 3))], 3)


def test_an_empty_table_is_refused(tmp_path: Path) -> None:
    """An empty output is almost always a silently failed adapter, and it passes most checks."""
    with pytest.raises(ValueError, match="no rows"):
        write(tmp_path, [], 3)


def test_the_sidecars_share_the_stem_of_the_table(tmp_path: Path) -> None:
    """`with_suffix` replaces only the final suffix, so `d.ml.csv` keeps the `d.ml` stem."""
    path = tmp_path / "d.ml.csv"
    assert manifest_path(path).name == "d.ml.manifest.json"
    assert samples_path(path).name == "d.ml.samples.csv"


def test_nothing_is_written_when_the_group_cannot_be_completed(tmp_path: Path) -> None:
    """The atomic-group promise: a table must never appear without its sidecar.

    Every guard runs before the first byte is staged, so a refusal leaves no partial file.
    """
    path = tmp_path / "d.ml.csv"
    with pytest.raises(ValueError):
        write_nexus_csv(
            path,
            [Row("s1", (0,), (1, 3))],
            n_features=3,
            n_targets=1,
            manifest={},
        )
    assert not path.exists()
    assert not manifest_path(path).exists()
    assert list(tmp_path.iterdir()) == []


def test_a_free_text_field_cannot_break_the_samples_sidecar(tmp_path: Path) -> None:
    """Descriptions contain commas and quotes; csv must quote them, not corrupt the file."""
    path = tmp_path / "d.ml.csv"
    write_nexus_csv(
        path,
        [Row("s1", (0,), (1, 3, 7))],
        n_features=3,
        n_targets=1,
        manifest={},
        sample_metadata=[
            {
                "sample_id": "s1",
                "row_sequence_blake2b": "x",
                "description": 'a, b "quoted", c\nnewline',
            }
        ],
    )
    with samples_path(path).open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["description"] == 'a, b "quoted", c\nnewline'
