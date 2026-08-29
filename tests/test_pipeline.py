"""End to end: records in, a validated Nexus CSV out."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from nexus_genomics.adapters.base import LoadedSource
from nexus_genomics.cleaning import Alphabet
from nexus_genomics.common import SampleRecord
from nexus_genomics.encoding import EncodingMode, LengthPolicy
from nexus_genomics.nexus_csv import manifest_path, samples_path
from nexus_genomics.pipeline import ConvertOptions, convert
from nexus_genomics.validation import validate_file


def source(records: list[SampleRecord], **kwargs: object) -> LoadedSource:
    defaults: dict[str, object] = {
        "n_targets": 1,
        "label_semantics": {0: "natural", 1: "engineered"},
        "alphabet": Alphabet.DNA_NUCLEOTIDE,
    }
    defaults.update(kwargs)
    return LoadedSource(records=records, **defaults)  # type: ignore[arg-type]


def run(tmp_path: Path, src: LoadedSource, options: ConvertOptions) -> Path:
    out = tmp_path / "t.ml.csv"
    convert(
        src,
        out,
        options,
        source_name="test",
        source_url="https://example.invalid",
        licence="test",
        target_description="test",
        notes="test",
        config_digest="deadbeef",
    )
    return out


def test_a_converted_file_passes_every_validation_check(tmp_path: Path) -> None:
    """The suite's spine: if this breaks, none of the other guarantees hold."""
    records = [SampleRecord(f"s{i}", "ACGT" * (i + 1), (i % 2,), "test", {}) for i in range(1, 9)]
    out = run(tmp_path, source(records), ConvertOptions(width=40))
    result = validate_file(out, sample_size=8, seed=1)
    failed = [c.name for c in result.checks if not c.passed]
    assert failed == [], failed
    assert result.stats["row_count"] == 8


def test_a_random_row_decodes_back_to_its_source_sequence(tmp_path: Path) -> None:
    """The end-to-end check the brief asks for: source -> clean -> encode -> pad -> CSV -> decode.

    Comparing digests rather than re-reading the source means this works without the raw
    download present, which matters for the two gated sources.
    """
    records = [SampleRecord("s1", "acgt nnrt\n", (0,), "test", {})]
    out = run(tmp_path, source(records), ConvertOptions(width=12))
    result = validate_file(out, sample_size=1, seed=0)
    check = next(c for c in result.checks if c.name.startswith("sampled_rows_decode"))
    assert check.passed, check.detail

    row = out.read_text(encoding="utf-8").splitlines()[1].split(",")
    assert row[0] == "s1"
    # "ACGTNNRT" -> A=1 C=3 G=7 T=20 N=14 N=14 R=18 T=20, then padding.
    assert row[2:] == ["1", "3", "7", "20", "14", "14", "18", "20", "0", "0", "0", "0"]


def test_padding_is_the_only_source_of_zero(tmp_path: Path) -> None:
    """No letter maps to 0, so a zero always means 'no residue here'."""
    records = [SampleRecord("s1", "ACGT", (0,), "test", {})]
    out = run(tmp_path, source(records), ConvertOptions(width=10))
    cells = out.read_text(encoding="utf-8").splitlines()[1].split(",")[2:]
    assert cells.count("0") == 6
    assert cells[:4] == ["1", "3", "7", "20"]


def test_truncation_is_recorded_per_sample_not_silently_dropped(tmp_path: Path) -> None:
    """Truncating without recording the loss is called out by name in the brief."""
    records = [SampleRecord("s1", "ACGT" * 10, (0,), "test", {})]
    out = run(tmp_path, source(records), ConvertOptions(width=5))
    meta = samples_path(out).read_text(encoding="utf-8").splitlines()
    assert "residues_removed_by_truncation" in meta[0]
    row = dict(zip(meta[0].split(","), meta[1].split(","), strict=True))
    assert row["original_length"] == "40"
    assert row["residues_removed_by_truncation"] == "35"

    manifest = json.loads(manifest_path(out).read_text(encoding="utf-8"))
    assert manifest["manifest"]["window"]["residues_removed_by_truncation"] == 35


def test_window_ids_from_one_sequence_never_collide(tmp_path: Path) -> None:
    """Colliding window ids would be refused by the writer, but only after the damage.

    Each window carries its parent id and its offsets so windows can be grouped back
    together -- which is what a leakage-free split needs.
    """
    records = [SampleRecord("s1", "ACGT" * 10, (0,), "test", {})]
    out = run(
        tmp_path,
        source(records),
        ConvertOptions(width=8, policy=LengthPolicy.WINDOW, stride=8, max_windows_per_sample=None),
    )
    lines = out.read_text(encoding="utf-8").splitlines()[1:]
    ids = [line.split(",")[0] for line in lines]
    assert len(ids) == len(set(ids)) == 5
    assert ids[0] == "s1__w00"

    meta = samples_path(out).read_text(encoding="utf-8").splitlines()
    header = meta[0].split(",")
    rows = [dict(zip(header, line.split(","), strict=True)) for line in meta[1:]]
    assert {r["parent_sample_id"] for r in rows} == {"s1"}
    assert [int(r["window_start"]) for r in rows] == [1, 9, 17, 25, 33]
    assert [int(r["window_end"]) for r in rows] == [8, 16, 24, 32, 40]


def test_a_split_grouped_by_parent_keeps_every_window_on_one_side(tmp_path: Path) -> None:
    """Windows of one sequence on both sides of a split is textbook leakage.

    The sidecar has to make the grouping possible without re-deriving it from the id string.
    """
    records = [SampleRecord(f"s{i}", "ACGT" * 6, (0,), "test", {}) for i in range(3)]
    out = run(
        tmp_path,
        source(records),
        ConvertOptions(width=8, policy=LengthPolicy.WINDOW, stride=8, max_windows_per_sample=None),
    )
    meta = samples_path(out).read_text(encoding="utf-8").splitlines()
    header = meta[0].split(",")
    rows = [dict(zip(header, line.split(","), strict=True)) for line in meta[1:]]
    groups: dict[str, set[str]] = {}
    for row in rows:
        groups.setdefault(row["parent_sample_id"], set()).add(row["sample_id"])
    assert len(groups) == 3
    assert all(len(v) == 3 for v in groups.values())


def test_the_window_cap_is_counted_rather_than_silently_applied(tmp_path: Path) -> None:
    """A cap that is not reported reads as 'this sequence was only that long'."""
    records = [SampleRecord("s1", "ACGT" * 20, (0,), "test", {})]
    out = run(
        tmp_path,
        source(records),
        ConvertOptions(width=8, policy=LengthPolicy.WINDOW, stride=8, max_windows_per_sample=3),
    )
    manifest = json.loads(manifest_path(out).read_text(encoding="utf-8"))
    assert manifest["manifest"]["window"]["windows_dropped_by_cap"] == 7
    assert len(out.read_text(encoding="utf-8").splitlines()) == 4

    # The RESIDUES the cap threw away must be counted too. `place` reports removed=0 for a
    # window because windowing discards nothing; the cap is what discards, and a sequence
    # that lost 70% of its length must not be recorded as having lost none.
    assert manifest["manifest"]["window"]["residues_dropped_by_cap"] == 56
    meta = samples_path(out).read_text(encoding="utf-8").splitlines()
    header = meta[0].split(",")
    rows = [dict(zip(header, line.split(","), strict=True)) for line in meta[1:]]
    assert sum(int(r["residues_dropped_by_window_cap"]) for r in rows) == 56


def test_a_record_with_no_label_is_refused(tmp_path: Path) -> None:
    """A missing label must never default to the negative class.

    That single silent default is how a dataset acquires fabricated 'natural' rows.
    """
    records = [SampleRecord("s1", "ACGT", None, "test", {})]
    with pytest.raises(ValueError, match="never default to the negative class"):
        run(tmp_path, source(records), ConvertOptions(width=4))


def test_an_unimplemented_encoding_mode_is_refused_with_an_explanation(tmp_path: Path) -> None:
    """Silently approximating codons as nucleotides is the confusion the brief warns about."""
    records = [SampleRecord("s1", "ACGT", (0,), "test", {})]
    with pytest.raises(NotImplementedError, match="reading frame"):
        run(tmp_path, source(records), ConvertOptions(width=4, mode=EncodingMode.CODON))
    with pytest.raises(NotImplementedError, match="protein that does not exist"):
        run(
            tmp_path,
            source(records),
            ConvertOptions(width=4, mode=EncodingMode.AMINO_ACID_TRANSLATED),
        )


def test_quarantined_records_are_excluded_from_the_table_and_kept_in_the_manifest(
    tmp_path: Path,
) -> None:
    """Dropping a bad record without recording it makes the row count unexplainable."""
    records = [
        SampleRecord("good", "ACGT", (0,), "test", {}),
        SampleRecord("bad", "ACGTXX", (1,), "test", {}),
    ]
    out = run(tmp_path, source(records), ConvertOptions(width=8))
    assert len(out.read_text(encoding="utf-8").splitlines()) == 2
    manifest = json.loads(manifest_path(out).read_text(encoding="utf-8"))
    quarantined = manifest["manifest"]["quarantined_samples"]
    assert [q["sample_id"] for q in quarantined] == ["bad"]
    assert manifest["manifest"]["cleaning"]["quarantined"] == 1


def test_rerunning_with_the_same_config_reproduces_the_content_hash(tmp_path: Path) -> None:
    """Reproducibility is a promise the sidecar makes; the timestamp must not break it."""
    records = [SampleRecord(f"s{i}", "ACGT" * (i + 1), (i % 2,), "test", {}) for i in range(4)]
    first = run(tmp_path / "a", source(records), ConvertOptions(width=20))
    second = run(tmp_path / "b", source(records), ConvertOptions(width=20))
    assert first.read_bytes() == second.read_bytes()

    hashes = {
        json.loads(manifest_path(p).read_text(encoding="utf-8"))["manifest"]["content_hash"]
        for p in (first, second)
    }
    assert len(hashes) == 1

    stamps = {
        json.loads(manifest_path(p).read_text(encoding="utf-8"))["manifest"]["generated_at"]
        for p in (first, second)
    }
    assert len(stamps) >= 1  # the timestamp may differ; the table bytes may not
