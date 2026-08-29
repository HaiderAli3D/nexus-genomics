"""Generic single-target projection behavior and its SeqScreen integration."""

from __future__ import annotations

from pathlib import Path

import pytest

from nexus_genomics.adapters.base import LoadedSource
from nexus_genomics.adapters.seqscreen import SeqScreenFunSoCAdapter
from nexus_genomics.cleaning import Alphabet
from nexus_genomics.common import SampleRecord


def _source(
    *,
    labels: list[tuple[int, ...] | None],
    n_targets: int,
    target_names: tuple[str, ...] = (),
    manifest_extra: dict[str, object] | None = None,
) -> LoadedSource:
    records = [
        SampleRecord(f"sample-{index}", "ACGT", label, "fixture", {"row": index})
        for index, label in enumerate(labels)
    ]
    return LoadedSource(
        records=records,
        n_targets=n_targets,
        label_semantics={0: "absent", 1: "present"},
        alphabet=Alphabet.DNA_NUCLEOTIDE,
        target_names=target_names,
        manifest_extra=manifest_extra or {},
        raw_input_hashes={"fixture.txt": "abc123"},
        sample_metadata_fields=("row",),
    )


def test_a_single_target_source_uses_the_original_source_and_an_empty_slug() -> None:
    """Copying a one-target source could needlessly change its established output bytes."""
    from nexus_genomics.projection import single_target_projections

    source = _source(labels=[(1,)], n_targets=1)

    projections = single_target_projections(source)

    assert len(projections) == 1
    assert projections[0].index == 0
    assert projections[0].name == ""
    assert projections[0].slug == ""
    assert projections[0].source is source
    assert source.target_names == ()


def test_multilabel_values_are_split_into_named_single_target_sources() -> None:
    """Selecting the wrong vector index silently assigns one mechanism's bits to another."""
    from nexus_genomics.projection import single_target_projections

    source = _source(
        labels=[(1, 0), (0, 1)],
        n_targets=2,
        target_names=("Beta / Toxin++", "Tox\u00edn"),
    )

    projections = single_target_projections(source)

    assert [(item.index, item.name, item.slug) for item in projections] == [
        (0, "Beta / Toxin++", "beta_toxin"),
        (1, "Tox\u00edn", "tox_n"),
    ]
    assert [[record.label for record in item.source.records] for item in projections] == [
        [(1,), (0,)],
        [(0,), (1,)],
    ]
    assert all(item.source.n_targets == 1 for item in projections)
    assert all(item.source.target_names == () for item in projections)


def test_projection_preserves_provenance_without_mutating_source_metadata() -> None:
    """Projection metadata must enrich each output without contaminating later projections."""
    from nexus_genomics.projection import single_target_projections

    source = _source(
        labels=[(1, 0)],
        n_targets=2,
        target_names=("First", "Second"),
        manifest_extra={"existing": {"kept": True}},
    )

    first = single_target_projections(source)[0].source

    assert first.records[0].sample_id == "sample-0"
    assert first.records[0].sequence == "ACGT"
    assert first.records[0].source == "fixture"
    assert first.records[0].metadata == {"row": 0}
    assert first.label_semantics is source.label_semantics
    assert first.alphabet is source.alphabet
    assert first.raw_input_hashes is source.raw_input_hashes
    assert first.sample_metadata_fields == source.sample_metadata_fields
    assert first.manifest_extra == {
        "existing": {"kept": True},
        "target_projection": {
            "original_target_index": 0,
            "original_target_name": "First",
            "original_target_count": 2,
        },
    }
    assert first.manifest_extra is not source.manifest_extra
    assert source.manifest_extra == {"existing": {"kept": True}}


def test_multilabel_projection_requires_one_name_for_every_target() -> None:
    """An unnamed index cannot produce a traceable or deterministic output filename."""
    from nexus_genomics.projection import single_target_projections

    source = _source(labels=[(1, 0)], n_targets=2, target_names=("Only one",))

    with pytest.raises(ValueError, match="exactly one target name"):
        single_target_projections(source)


def test_a_target_name_that_has_no_ascii_alphanumeric_slug_is_refused() -> None:
    """An empty suffix would collide with the unsuffixed source output."""
    from nexus_genomics.projection import single_target_projections

    source = _source(labels=[(1, 0)], n_targets=2, target_names=("!!!", "Valid"))

    with pytest.raises(ValueError, match="empty slug"):
        single_target_projections(source)


def test_target_names_that_normalize_to_the_same_slug_are_refused() -> None:
    """Duplicate projected filenames would overwrite one target with another target's data."""
    from nexus_genomics.projection import single_target_projections

    source = _source(
        labels=[(1, 0)],
        n_targets=2,
        target_names=("Toxin Activity", "toxin-activity"),
    )

    with pytest.raises(ValueError, match="duplicate slug"):
        single_target_projections(source)


def test_a_missing_record_label_is_refused_before_any_projection_is_returned() -> None:
    """Defaulting a missing bit to zero would fabricate a negative class."""
    from nexus_genomics.projection import single_target_projections

    source = _source(labels=[(1, 0), None], n_targets=2, target_names=("First", "Second"))

    with pytest.raises(ValueError, match="missing label"):
        single_target_projections(source)


def test_a_record_label_with_the_wrong_width_is_refused() -> None:
    """A short vector shifts or drops target meaning while still containing valid integers."""
    from nexus_genomics.projection import single_target_projections

    source = _source(labels=[(1,)], n_targets=2, target_names=("First", "Second"))

    with pytest.raises(ValueError, match="exactly 2 values"):
        single_target_projections(source)


def test_an_empty_source_is_refused_instead_of_producing_empty_tables() -> None:
    """An empty projection usually means loading failed and must not look successful."""
    from nexus_genomics.projection import single_target_projections

    source = _source(labels=[], n_targets=2, target_names=("First", "Second"))

    with pytest.raises(ValueError, match="no records"):
        single_target_projections(source)


def test_an_empty_slug_leaves_the_base_path_unchanged() -> None:
    """Single-target sources retain their established path, even before suffix validation."""
    from nexus_genomics.projection import projected_path

    base_path = Path("output.csv")

    assert projected_path(base_path, "") == base_path


def test_a_projection_slug_is_inserted_before_the_exact_ml_csv_extension() -> None:
    """Appending after .csv would break sidecar discovery from the shared filename stem."""
    from nexus_genomics.projection import projected_path

    assert projected_path(Path("out/funsoc_demo.ml.csv"), "toxin") == Path(
        "out/funsoc_demo__toxin.ml.csv"
    )


def test_a_suffixed_projection_refuses_a_path_without_the_exact_ml_csv_extension() -> None:
    """Guessing at another suffix could write a file the existing reader cannot recognize."""
    from nexus_genomics.projection import projected_path

    with pytest.raises(ValueError, match=r"\.ml\.csv"):
        projected_path(Path("funsoc_demo.csv"), "toxin")


def test_seqscreen_publishes_target_names_in_the_same_order_as_its_vectors(
    tmp_path: Path,
) -> None:
    """Projection names must use the exact index order used to set each FunSoC bit."""
    raw = tmp_path / "seqscreen"
    funsoc_dir = raw / "funsoc"
    funsoc_dir.mkdir(parents=True)
    expected = tuple(f"FunSoC{index:02d}" for index in range(32))
    for index, name in enumerate(expected):
        (funsoc_dir / f"{name}_v1.0_2026.txt").write_text(f"P{index:05d}\n", encoding="utf-8")
    (raw / "uniprot_sequences.tsv").write_text(
        "".join(f"P{index:05d}\tMKT\n" for index in range(32)), encoding="utf-8"
    )

    loaded = SeqScreenFunSoCAdapter().load(raw, {})

    assert loaded.target_names == expected
    assert loaded.manifest_extra["target_index_to_funsoc"] == dict(enumerate(expected))
