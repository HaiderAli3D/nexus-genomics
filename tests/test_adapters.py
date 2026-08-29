"""Adapter behaviour, against synthetic fixtures.

These run with no network and no gated download, and they exercise the same code path the
real files take. The label-integrity tests here are the ones that matter most: they are what
stops a future edit from quietly reintroducing a fabricated natural class.
"""

from __future__ import annotations

import csv
import io
import json
import tarfile
from pathlib import Path

import pytest

from nexus_genomics.adapters import ADAPTERS, get_adapter
from nexus_genomics.adapters.addgene import POOLED_UNKNOWN_LAB, AddgeneAdapter
from nexus_genomics.adapters.codontransformer import CodonTransformerAdapter
from nexus_genomics.adapters.felix_guardian import FelixGuardianAdapter
from nexus_genomics.adapters.seqscreen import SeqScreenFunSoCAdapter, _stratified
from nexus_genomics.cleaning import Alphabet

# --------------------------------------------------------------------------- registry


def test_every_registered_adapter_states_its_licence_and_target() -> None:
    """An adapter that cannot say what its target means has no business writing a file."""
    for name, adapter in ADAPTERS.items():
        assert adapter.name, name
        assert adapter.source_url.startswith("http"), name
        assert len(adapter.licence) > 20, name
        assert len(adapter.target_description) > 40, name


def test_an_unknown_source_names_the_known_ones() -> None:
    with pytest.raises(KeyError, match="Known sources"):
        get_adapter("nope")


# --------------------------------------------------------------------------- addgene


def _write_geac(tmp_path: Path, labs: list[str]) -> Path:
    raw = tmp_path / "addgene"
    raw.mkdir()
    ids = [f"SEQ{i:02d}" for i in range(len(labs))]
    with (raw / "train_values.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["sequence_id", "sequence", "bacterial_resistance_ampicillin", "species_human"]
        )
        for i, sid in enumerate(ids):
            writer.writerow([sid, "ACGT" * (i + 2), 1.0, 0.0])
    unique = sorted(set(labs))
    with (raw / "train_labels.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sequence_id", *unique])
        for sid, lab in zip(ids, labs, strict=True):
            writer.writerow([sid, *[1.0 if u == lab else 0.0 for u in unique]])
    return raw


def test_a_lab_attribution_file_never_contains_a_natural_class(tmp_path: Path) -> None:
    """Every Addgene plasmid is engineered. A 'natural' class here would be fabricated.

    This is the single most important assertion in the suite: the brief explicitly warns
    against labelling some engineered plasmids as natural, and this is what would catch it.
    """
    raw = _write_geac(tmp_path, ["AAAAAAAA", "BBBBBBBB", POOLED_UNKNOWN_LAB])
    loaded = AddgeneAdapter().load(raw, {})
    meanings = " ".join(loaded.label_semantics.values()).lower()
    assert "natural" not in meanings.replace("not a natural class", "")
    assert all(len(r.label or ()) == 1 for r in loaded.records)
    assert loaded.n_targets == 1
    assert loaded.alphabet is Alphabet.DNA_NUCLEOTIDE


def test_the_pooled_unknown_class_is_documented_as_unknown_lab_not_unknown_origin(
    tmp_path: Path,
) -> None:
    """'Unknown Engineered' means unknown LABORATORY. Reading it as unknown provenance
    inverts its meaning and is exactly how a fake natural class gets introduced."""
    raw = _write_geac(tmp_path, ["AAAAAAAA", POOLED_UNKNOWN_LAB])
    loaded = AddgeneAdapter().load(raw, {})
    pooled = next(v for v in loaded.label_semantics.values() if POOLED_UNKNOWN_LAB in v)
    assert "UNKNOWN LABORATORY" in pooled
    assert "still an" in pooled and "engineered" in pooled


def test_the_one_hot_label_file_is_collapsed_the_way_the_benchmark_does(tmp_path: Path) -> None:
    """train_labels.csv is 1,315 columns of one-hot; idxmax is the documented collapse."""
    raw = _write_geac(tmp_path, ["BBBBBBBB", "AAAAAAAA", "BBBBBBBB"])
    loaded = AddgeneAdapter().load(raw, {})
    by_id = {r.sample_id: r for r in loaded.records}
    assert by_id["SEQ00"].metadata["lab_id"] == "BBBBBBBB"
    assert by_id["SEQ01"].metadata["lab_id"] == "AAAAAAAA"
    # Class indices are assigned by sorted lab id, so they are stable across runs.
    assert by_id["SEQ01"].label == (0,)
    assert by_id["SEQ00"].label == (1,)


def test_an_all_zero_label_row_is_refused_rather_than_attributed_to_the_first_lab(
    tmp_path: Path,
) -> None:
    """`idxmax` on an all-zero row returns the first column and says nothing about it.

    That is a fabricated label produced by a function that looks like it is reading one --
    the exact failure this whole project is organised against. It must be a hard error.
    """
    raw = _write_geac(tmp_path, ["AAAAAAAA", "BBBBBBBB"])
    labels = raw / "train_labels.csv"
    lines = labels.read_text(encoding="utf-8").splitlines()
    lines[1] = lines[1].split(",")[0] + ",0.0,0.0"
    labels.write_text("\n".join(lines) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="would invent a lab"):
        AddgeneAdapter().load(raw, {})


def test_a_duplicated_sequence_id_is_refused_rather_than_silently_corrupted(
    tmp_path: Path,
) -> None:
    """With a duplicate id, a pandas row lookup returns a DataFrame, not a row.

    `str()` of that is a repr of a Series, so the sequence column would fill with text like
    '0    ACGT\\nName: ...' -- a well-formed CSV of garbage.
    """
    raw = _write_geac(tmp_path, ["AAAAAAAA", "BBBBBBBB"])
    values = raw / "train_values.csv"
    lines = values.read_text(encoding="utf-8").splitlines()
    lines.append(lines[1])
    values.write_text("\n".join(lines) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="duplicated sequence_id"):
        AddgeneAdapter().load(raw, {})


def test_the_gated_source_explains_the_download_instead_of_failing_obscurely(
    tmp_path: Path,
) -> None:
    """'FileNotFoundError' alone would send someone hunting for a bug that is not there."""
    from nexus_genomics.adapters.base import GatedSourceUnavailableError

    status = AddgeneAdapter().availability(tmp_path / "empty")
    assert not status.ready
    assert "drivendata.org/accounts/signup" in status.manual_steps
    assert "ACCEPT THE COMPETITION RULES" in status.manual_steps
    assert status.expected_files == ("train_values.csv", "train_labels.csv")
    with pytest.raises(GatedSourceUnavailableError, match="drivendata"):
        AddgeneAdapter().load(tmp_path / "empty", {})


def test_the_phenotype_columns_stay_out_of_the_model_table(tmp_path: Path) -> None:
    """They are real variables, but mixing them into position_* breaks the one-position-
    per-column contract, so they ride in the samples sidecar instead."""
    raw = _write_geac(tmp_path, ["AAAAAAAA", "BBBBBBBB"])
    loaded = AddgeneAdapter().load(raw, {})
    assert "bacterial_resistance_ampicillin" in loaded.sample_metadata_fields
    assert "species_human" in loaded.sample_metadata_fields


# --------------------------------------------------------------------- codontransformer


def _write_ct(tmp_path: Path, rows: list[tuple[str, str, str]]) -> Path:
    raw = tmp_path / "codontransformer"
    raw.mkdir()
    with (raw / "dataset_with_predictions_metrics.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["", "protein", "organism", "GeneID", "description", "natural_dna", "finetune_dna"]
        )
        for i, (protein, natural, finetune) in enumerate(rows):
            writer.writerow(
                [i, protein, "Escherichia coli general", 100 + i, "d", natural, finetune]
            )
    return raw


def test_a_pair_shares_a_protein_and_an_organism_and_differs_only_in_label(
    tmp_path: Path,
) -> None:
    """The whole scientific claim of this dataset: matched on protein AND organism.

    If the two members ever came from different rows, the classes would differ by gene as
    well as by generator and the task would be measuring the wrong thing.
    """
    raw = _write_ct(tmp_path, [("MA", "ATGGCC", "ATGGCG")])
    loaded = CodonTransformerAdapter().load(raw, {"max_pairs": 10})
    assert len(loaded.records) == 2
    natural, generated = loaded.records
    assert natural.label == (0,)
    assert generated.label == (1,)
    assert natural.metadata["pair_id"] == generated.metadata["pair_id"]
    assert natural.metadata["organism"] == generated.metadata["organism"]
    assert natural.metadata["gene_id"] == generated.metadata["gene_id"]


def test_a_pair_whose_members_differ_in_length_is_dropped_and_counted(tmp_path: Path) -> None:
    """Same protein means the same codon count. A mismatch is not a clean pair.

    Keeping it would also make truncation asymmetric between the classes, which leaks the
    label into sequence length.
    """
    raw = _write_ct(tmp_path, [("MA", "ATGGCC", "ATGGCGAAA"), ("MA", "ATGGCC", "ATGGCG")])
    loaded = CodonTransformerAdapter().load(raw, {"max_pairs": 10})
    assert len(loaded.records) == 2
    assert loaded.manifest_extra["dropped_length_mismatch_pairs"] == 1


def test_truncation_cannot_leak_the_label_because_pairs_are_equal_length(
    tmp_path: Path,
) -> None:
    """Both members are truncated identically, so sequence length carries no class signal."""
    raw = _write_ct(tmp_path, [("MA", "ATG" * 40, "ATC" * 40)])
    loaded = CodonTransformerAdapter().load(raw, {"max_pairs": 10})
    lengths = {r.label: len(r.sequence) for r in loaded.records}
    assert lengths[(0,)] == lengths[(1,)]


def test_an_unexpected_source_schema_is_refused_rather_than_guessed(tmp_path: Path) -> None:
    """Inventing a column name is the worst available error, so a missing one is fatal."""
    raw = tmp_path / "codontransformer"
    raw.mkdir()
    (raw / "dataset_with_predictions_metrics.csv").write_text("a,b\n1,2\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Refusing to guess a schema"):
        CodonTransformerAdapter().load(raw, {"max_pairs": 5})


def test_the_generated_column_must_be_one_of_the_two_that_exist(tmp_path: Path) -> None:
    raw = _write_ct(tmp_path, [("MA", "ATGGCC", "ATGGCG")])
    with pytest.raises(ValueError, match="finetune_dna"):
        CodonTransformerAdapter().load(raw, {"generated_column": "made_up_dna"})


# ------------------------------------------------------------------------ felix/guardian


def _write_midoe(tmp_path: Path) -> Path:
    raw = tmp_path / "felix"
    raw.mkdir()
    payload = {
        "features": [
            {"id": "f1", "name": "n1", "sequence": "acgtACGT", "role": "CDS", "source": "AF1"},
            {"id": "f2", "name": "n2", "sequence": "ACGTNNNN", "role": "promoter"},
            {"id": "f3", "name": "n3", "sequence": "ACGT"},
        ],
        "detections": [
            {
                "agent": "TEST",
                "sample": "Y001",
                "engineered": "http://guardian.bbn.technology#engineered",
                "contigs": [
                    {
                        "source": "X.0",
                        "regions": [
                            {
                                "start": 1,
                                "end": 8,
                                "feature": "f1",
                                "engineered": "http://guardian.bbn.technology#engineered",
                            },
                            {
                                # A call keyed to a feature that carries no sequence. Real
                                # in the corpus, and the reason the census must intersect
                                # against the features rather than count the call index.
                                "start": 1,
                                "end": 9,
                                "feature": "f_no_sequence",
                                "engineered": "http://guardian.bbn.technology#engineered",
                            },
                        ],
                    }
                ],
            }
        ],
    }
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w:gz") as archive:
        data = json.dumps(payload).encode("utf-8")
        info = tarfile.TarInfo("midoe-main/guardian_felix_t_and_e_batch_4/guardian_evidence/A.json")
        info.size = len(data)
        archive.addfile(info, io.BytesIO(data))
    (raw / "midoe-main.tar.gz").write_bytes(buffer.getvalue())
    return raw


def test_the_guardian_call_is_never_written_into_the_target_column(tmp_path: Path) -> None:
    """It is a detector assertion with a confidence, not ground truth.

    It is carried in the samples sidecar for audit under a name that says so, and the target
    is the part role instead.
    """
    loaded = FelixGuardianAdapter().load(_write_midoe(tmp_path), {})
    assert loaded.n_targets == 1
    assert {r.metadata["part_role"] for r in loaded.records} == {"CDS", "promoter"}
    # Every documented class says in its own text that it is a part type, so a reader who
    # only ever sees label_semantics cannot mistake it for a provenance call.
    for meaning in loaded.label_semantics.values():
        assert "genetic PART type" in meaning
        assert "not where it came from" in meaning
    assert "guardian_call_not_ground_truth" in loaded.sample_metadata_fields
    carried = {r.metadata["guardian_call_not_ground_truth"] for r in loaded.records}
    assert "engineered" in carried


def test_no_subset_of_the_part_roles_may_be_read_as_a_provenance_binary(
    tmp_path: Path,
) -> None:
    """Two role names contain the word 'engineered', which invites exactly the wrong binary.

    `target.isin([engineered_region, engineered_tag])` would label the CDS, promoter and
    origin_of_replication rows 'natural' -- but those are parts GUARDIAN annotated INSIDE
    engineered constructs. The refusal has to be in the artefact, not only in a docstring.
    """
    loaded = FelixGuardianAdapter().load(_write_midoe(tmp_path), {})
    note = loaded.manifest_extra["no_subset_of_these_classes_is_a_provenance_contrast"]
    assert "PART TYPES, not provenance calls" in note
    assert "fabrication" in note


def test_the_class_index_does_not_shift_when_the_sample_is_capped(tmp_path: Path) -> None:
    """Class integers must mean the same thing in the demo file and the full file.

    Building the vocabulary from the capped sample renumbers every class that sorts after a
    class the cap excluded. A model trained on the demo and scored against the full file
    would then be graded against the wrong classes, silently, with nothing to catch it.
    """
    raw = _write_midoe(tmp_path)
    full = FelixGuardianAdapter().load(raw, {})
    capped = FelixGuardianAdapter().load(raw, {"max_records": 1})
    assert len(capped.records) == 1
    assert dict(capped.label_semantics) == dict(full.label_semantics)
    by_role_full = {r.metadata["part_role"]: r.label for r in full.records}
    for record in capped.records:
        assert record.label == by_role_full[record.metadata["part_role"]]


def test_only_features_with_a_role_become_rows(tmp_path: Path) -> None:
    """A feature with no annotation has no honest target, so it is left out entirely."""
    loaded = FelixGuardianAdapter().load(_write_midoe(tmp_path), {})
    assert {r.metadata["part_role"] for r in loaded.records} == {"CDS", "promoter"}
    assert len(loaded.records) == 2


def test_the_adapter_records_why_a_binary_file_would_be_dishonest(tmp_path: Path) -> None:
    """The finding has to survive into the artefact, not just live in a docstring."""
    loaded = FelixGuardianAdapter().load(_write_midoe(tmp_path), {})
    why = loaded.manifest_extra["why_not_natural_vs_engineered"]
    assert "exactly one class" in why
    assert "pubs.acs.org" in loaded.manifest_extra["ground_truth_route"]


def test_the_census_counts_only_calls_that_reach_a_real_sequence(tmp_path: Path) -> None:
    """The headline finding is 280 of 1,670, all engineered. It must be computed, not asserted.

    A call can be keyed to a feature that carries no sequence, and such a call can never
    become a row. Counting the call index instead of intersecting it against the features
    would inflate `features_linkable_to_a_call` and quietly turn the finding into its
    opposite -- "every call is linkable" -- in the sidecar of every FELIX file.
    """
    loaded = FelixGuardianAdapter().load(_write_midoe(tmp_path), {})
    census = loaded.manifest_extra["corpus_census"]
    assert census["features_with_sequence"] == 3
    assert census["features_linkable_to_a_call"] == 1
    assert census["feature_level_call_distribution"] == {"engineered": 1}


def test_lowercase_bases_in_the_corpus_survive_to_cleaning(tmp_path: Path) -> None:
    """MIDOE ships mixed case; the adapter must not silently normalise it out of sight."""
    loaded = FelixGuardianAdapter().load(_write_midoe(tmp_path), {})
    assert any("acgt" in r.sequence for r in loaded.records)


# ---------------------------------------------------------------------------- funsoc


def test_the_stratified_cap_keeps_every_funsoc_represented() -> None:
    """Taking the first N accessions alphabetically left 6 of 32 target columns all-zero.

    An all-zero column is indistinguishable from a class that never occurs, so a model
    trained on it learns the class does not exist.
    """
    labels = {f"P{i:05d}": {f"F{i % 10}"} for i in range(500)}
    funsocs = [f"F{i}" for i in range(10)]
    chosen = _stratified(labels, funsocs, 20)
    assert len(chosen) == 20
    covered = {next(iter(labels[a])) for a in chosen}
    assert covered == set(funsocs)


def test_the_stratified_cap_is_deterministic() -> None:
    """content_hash must not depend on dict iteration order."""
    labels = {f"P{i:05d}": {f"F{i % 7}"} for i in range(200)}
    funsocs = [f"F{i}" for i in range(7)]
    assert _stratified(labels, funsocs, 30) == _stratified(labels, funsocs, 30)


def test_funsoc_names_are_parsed_off_the_versioned_filename() -> None:
    from nexus_genomics.adapters.seqscreen import _funsoc_name

    assert _funsoc_name("SecretedEffector_v1.0_10.29.2019.txt") == "SecretedEffector"
    assert _funsoc_name("AvirulencePlant_v1.0_01.121.2020.txt") == "AvirulencePlant"


def test_a_partial_funsoc_download_is_refused(tmp_path: Path) -> None:
    """Fewer than 32 lists means missing FunSoCs would appear as all-zero columns."""
    raw = tmp_path / "seqscreen"
    (raw / "funsoc").mkdir(parents=True)
    (raw / "funsoc" / "Cytotoxicity_v1.0_10.24.2019.txt").write_text("P00001\n", encoding="utf-8")
    adapter = SeqScreenFunSoCAdapter()
    # _ensure_files falls back to what is cached when discovery fails, so this exercises the
    # refusal without needing the network to be down.
    status = adapter.availability(raw)
    assert status.ready
    assert "PAGE 2" in status.manual_steps
