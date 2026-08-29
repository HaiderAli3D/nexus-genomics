"""Proof that the validator fails when it should.

A validator that only ever passes is worse than none: it converts "we did not look" into
"we checked". Every check with real teeth gets a test that deliberately breaks the file and
asserts the check catches it.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from nexus_genomics.adapters.base import LoadedSource
from nexus_genomics.cleaning import Alphabet
from nexus_genomics.common import SampleRecord, canonical_json, hash_file
from nexus_genomics.encoding import LengthPolicy
from nexus_genomics.nexus_csv import manifest_path, samples_path
from nexus_genomics.pipeline import ConvertOptions, convert
from nexus_genomics.validation import compare_runs, render_markdown, validate_file


def build(tmp_path: Path, records: list[SampleRecord], width: int = 12, n_targets: int = 1) -> Path:
    out = tmp_path / "t.ml.csv"
    convert(
        LoadedSource(
            records=records,
            n_targets=n_targets,
            label_semantics={0: "natural", 1: "engineered"},
            alphabet=Alphabet.DNA_NUCLEOTIDE,
        ),
        out,
        ConvertOptions(width=width, policy=LengthPolicy.PAD_OR_TRUNCATE),
        source_name="t",
        source_url="https://example.invalid",
        licence="t",
        target_description="t",
        notes="t",
        config_digest="d",
    )
    return out


def default_records() -> list[SampleRecord]:
    return [
        SampleRecord("s1", "ACGT", (0,), "t", {}),
        SampleRecord("s2", "ACGTAC", (1,), "t", {}),
        SampleRecord("s3", "TTTT", (0,), "t", {}),
        SampleRecord("s4", "GGGGCC", (1,), "t", {}),
    ]


def failed(path: Path) -> list[str]:
    return [c.name for c in validate_file(path, sample_size=4).checks if not c.passed]


def rewrite_samples(path: Path, rows: list[str], *, update_hash: bool) -> None:
    """Rewrite a synthetic samples sidecar and optionally make its byte hash self-consistent."""
    side = samples_path(path)
    side.write_text("\n".join(rows) + "\n", encoding="utf-8")
    if update_hash:
        companion = manifest_path(path)
        payload = json.loads(companion.read_text(encoding="utf-8"))
        payload["manifest"]["samples_content_hash"] = hash_file(side)
        companion.write_text(canonical_json(payload), encoding="utf-8")


def test_a_clean_file_passes_everything(tmp_path: Path) -> None:
    assert failed(build(tmp_path, default_records())) == []


def test_padding_in_the_middle_of_a_row_is_caught(tmp_path: Path) -> None:
    """The invariant is residues-then-padding. A 0 before a residue means a spliced gap."""
    path = build(tmp_path, default_records())
    lines = path.read_text(encoding="utf-8").splitlines()
    cells = lines[1].split(",")
    cells[2] = "0"  # blank out position_1, leaving non-zero values after it
    lines[1] = ",".join(cells)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    assert "padding_only_ever_trails_the_sequence" in failed(path)


def test_an_out_of_range_token_is_caught(tmp_path: Path) -> None:
    """27 is not a letter and 0 is padding; anything else did not come from this cipher."""
    path = build(tmp_path, default_records())
    lines = path.read_text(encoding="utf-8").splitlines()
    cells = lines[1].split(",")
    cells[2] = "99"
    lines[1] = ",".join(cells)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    assert "tokens_within_configured_range" in failed(path)


def test_a_duplicated_primary_key_is_caught(tmp_path: Path) -> None:
    path = build(tmp_path, default_records())
    lines = path.read_text(encoding="utf-8").splitlines()
    lines.append(lines[1])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    names = failed(path)
    assert "primary_key_unique" in names


def test_an_edited_table_no_longer_matches_its_content_hash(tmp_path: Path) -> None:
    """The digest is what makes 'this file is the one the manifest describes' checkable."""
    path = build(tmp_path, default_records())
    text = path.read_text(encoding="utf-8")
    path.write_text(text.replace("s1", "sX", 1), encoding="utf-8")
    assert "content_hash_matches_the_table_on_disk" in failed(path)


def test_an_all_padding_row_is_caught(tmp_path: Path) -> None:
    """A row with no residues is well-formed and meaningless."""
    path = build(tmp_path, default_records())
    lines = path.read_text(encoding="utf-8").splitlines()
    width = len(lines[0].split(",")) - 2
    lines.append(",".join(["s9", "0", *(["0"] * width)]))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    assert "no_row_is_silently_empty" in failed(path)


def test_an_undocumented_label_value_is_caught(tmp_path: Path) -> None:
    """label_semantics is the file's own statement of what may appear in the target."""
    path = build(tmp_path, default_records())
    lines = path.read_text(encoding="utf-8").splitlines()
    cells = lines[1].split(",")
    cells[1] = "7"
    lines[1] = ",".join(cells)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    names = failed(path)
    assert "labels_only_documented_values" in names


def test_an_empty_present_label_map_fails_documented_value_validation(tmp_path: Path) -> None:
    """An empty declared vocabulary cannot document any observed categorical target."""
    path = build(tmp_path, default_records())
    companion = manifest_path(path)
    payload = json.loads(companion.read_text(encoding="utf-8"))
    payload["manifest"]["label_semantics"] = {}
    companion.write_text(canonical_json(payload), encoding="utf-8")

    names = failed(path)
    assert "labels_only_documented_values" in names
    assert "label_dependent_checks_skipped" in names


def test_an_all_zero_target_column_is_caught(tmp_path: Path) -> None:
    """A class with no positives cannot be learned and reads as 'never occurs'."""
    records = [SampleRecord(f"s{i}", "ACGT", (0,), "t", {}) for i in range(4)]
    path = build(tmp_path, records)
    assert "every_target_column_carries_at_least_one_positive" in failed(path)


def test_a_missing_sidecar_stops_validation_immediately(tmp_path: Path) -> None:
    """Without the sidecar there is nothing to validate the table against."""
    path = build(tmp_path, default_records())
    manifest_path(path).unlink()
    assert "sidecar_present_and_matches_header" in failed(path)


def test_a_tampered_row_fails_the_round_trip(tmp_path: Path) -> None:
    """Changing a residue must be caught by the per-sample digest, not go unnoticed."""
    path = build(tmp_path, default_records())
    lines = path.read_text(encoding="utf-8").splitlines()
    cells = lines[1].split(",")
    cells[2] = "3"  # A -> C
    lines[1] = ",".join(cells)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    names = failed(path)
    assert "sampled_rows_decode_back_to_their_source_sequence" in names


def test_the_round_trip_catches_an_off_by_one_in_the_encoder(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The check must compare against the SOURCE, not against the encoder's own output.

    An earlier version hashed `decode(placed.tokens)` on both sides, which made the check a
    tautology: this exact mutation -- dropping the first residue of every row -- passed on
    every row of every file. Injecting the bug and asserting the validator now fails is the
    only way to prove the check has teeth.
    """
    import nexus_genomics.pipeline as pipeline_module
    from nexus_genomics.encoding import PAD_TOKEN, Placed

    def off_by_one(tokens: list[int], width: int, _policy: object, _stride: int) -> list[Placed]:
        kept = tokens[1 : width + 1]  # the bug: skips the first residue
        pad = width - len(kept)
        return [
            Placed(
                tokens=tuple(kept) + (PAD_TOKEN,) * pad,
                start=1,
                end=len(kept),
                pad_count=pad,
                removed=max(len(tokens) - len(kept), 0),
            )
        ]

    monkeypatch.setattr(pipeline_module, "place", off_by_one)
    path = build(tmp_path, default_records())
    assert "sampled_rows_decode_back_to_their_source_sequence" in failed(path)


def test_a_non_integer_feature_cell_is_reported_not_crashed_on(tmp_path: Path) -> None:
    """A blank cell is what a spreadsheet round-trip leaves behind.

    It makes pandas read the column as float64 with NaN, and the token-range check does
    int(values.min()). Raising there would turn the validator into a crash on exactly the
    corruption it exists to detect.
    """
    path = build(tmp_path, default_records())
    lines = path.read_text(encoding="utf-8").splitlines()
    cells = lines[1].split(",")
    cells[3] = ""
    lines[1] = ",".join(cells)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    names = failed(path)
    assert "all_position_columns_integer" in names
    assert "remaining_checks_skipped" in names


def test_a_non_integer_target_cell_is_reported_without_coercion_or_a_crash(
    tmp_path: Path,
) -> None:
    """``int(1.5)`` either crashes or fabricates class 1; neither is validation."""
    path = build(tmp_path, default_records())
    lines = path.read_text(encoding="utf-8").splitlines()
    cells = lines[1].split(",")
    cells[1] = "1.5"
    lines[1] = ",".join(cells)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    result = validate_file(path, sample_size=4)
    names = [check.name for check in result.checks]
    failed_names = [check.name for check in result.checks if not check.passed]

    assert "all_target_columns_integer" in failed_names
    assert "label_dependent_checks_skipped" in failed_names
    assert "content_hash_matches_the_table_on_disk" in names
    assert "sampled_rows_decode_back_to_their_source_sequence" in names


def test_a_padding_failure_stops_the_decode_dependent_checks(tmp_path: Path) -> None:
    """decode() skips a 0 wherever it sits, so interleaved padding yields a plausible lie.

    Reporting 'round trip: pass' underneath 'padding invariant: FAIL' would tell a reader
    the file is mostly fine when in fact the sequence is wrong.
    """
    path = build(tmp_path, default_records())
    lines = path.read_text(encoding="utf-8").splitlines()
    cells = lines[1].split(",")
    cells[2] = "0"
    lines[1] = ",".join(cells)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    names = [c.name for c in validate_file(path, sample_size=4).checks]
    assert "padding_only_ever_trails_the_sequence" in names
    assert "remaining_checks_skipped" in names
    assert "sampled_rows_decode_back_to_their_source_sequence" not in names


def test_a_sidecar_that_contradicts_its_own_manifest_is_caught(tmp_path: Path) -> None:
    path = build(tmp_path, default_records())
    side = manifest_path(path)
    payload = json.loads(side.read_text(encoding="utf-8"))
    payload["manifest"]["n_samples"] = 999
    side.write_text(json.dumps(payload), encoding="utf-8")
    assert "sidecar_counts_match_manifest" in failed(path)


def test_duplicate_sequences_and_label_conflicts_are_reported_not_failed(
    tmp_path: Path,
) -> None:
    """Both are facts about the data a reader must know, not defects in the conversion.

    They are reported with detail so nobody has to grep the CSV to discover them.
    """
    records = [
        SampleRecord("a", "ACGT", (0,), "t", {}),
        SampleRecord("b", "ACGT", (1,), "t", {}),
        SampleRecord("c", "TTTT", (1,), "t", {}),
    ]
    result = validate_file(build(tmp_path, records), sample_size=3)
    assert result.ok
    dupes = next(c for c in result.checks if c.name == "duplicate_sequences_reported")
    conflict = next(
        c for c in result.checks if c.name == "conflicting_labels_on_identical_sequences_reported"
    )
    assert "1 distinct sequence" in dupes.detail
    assert "1 sequence" in conflict.detail


def test_the_report_records_class_counts_and_balance(tmp_path: Path) -> None:
    result = validate_file(build(tmp_path, default_records()), sample_size=4)
    assert result.stats["class_counts"]["target"] == {"0": 2, "1": 2}
    assert result.stats["class_balance_min_over_max"] == 1.0
    assert result.stats["row_count"] == 4
    assert result.stats["column_count"] == 14


def test_the_markdown_and_json_reports_come_from_the_same_object(tmp_path: Path) -> None:
    """Two renderers over one result cannot disagree; two separate walks of the data can."""
    result = validate_file(build(tmp_path, default_records()), sample_size=4)
    markdown = render_markdown([result], {"config_hash": "abc"})
    payload = result.to_json_dict()
    assert payload["n_checks"] == len(result.checks)
    assert "1 of 1 file(s) pass every check." in markdown
    for check in result.checks:
        assert check.name in markdown


def test_two_identical_runs_compare_equal_and_a_changed_one_does_not(tmp_path: Path) -> None:
    first = build(tmp_path / "a", default_records())
    second = build(tmp_path / "b", default_records())
    assert compare_runs(first, second).passed

    third = build(tmp_path / "c", default_records()[:3])
    assert not compare_runs(first, third).passed


def test_the_samples_sidecar_covers_every_row(tmp_path: Path) -> None:
    """A row with no provenance entry cannot be round-tripped or audited."""
    path = build(tmp_path, default_records())
    table = path.read_text(encoding="utf-8").splitlines()[1:]
    side = samples_path(path).read_text(encoding="utf-8").splitlines()[1:]
    assert len(table) == len(side) == 4


def test_numeric_looking_sample_ids_keep_their_exact_text_during_validation(
    tmp_path: Path,
) -> None:
    """Default pandas inference must not collapse ``001`` and break an exact sidecar join."""
    path = build(
        tmp_path,
        [
            SampleRecord("001", "ACGT", (0,), "t", {}),
            SampleRecord("002", "TGCA", (1,), "t", {}),
        ],
    )

    checks = {check.name: check for check in validate_file(path).checks}

    assert checks["samples_sidecar_sample_ids_match_table_order"].passed
    assert checks["primary_key_unique"].passed


def test_numeric_looking_sample_ids_remain_distinct_in_the_round_trip(
    tmp_path: Path,
) -> None:
    """The semantic spot check must not compare ``001`` against metadata for ``1``.

    The complete sidecar hash protects the bytes, but pandas inference in the separate
    sampled round trip used to collapse both keys to integer 1 and select the wrong source
    digest. That creates a false corruption report for a valid pair of exact textual IDs.
    """
    path = build(
        tmp_path,
        [
            SampleRecord("001", "ACGT", (0,), "t", {}),
            SampleRecord("1", "TGCA", (1,), "t", {}),
        ],
    )

    checks = {check.name: check for check in validate_file(path, sample_size=2).checks}

    assert checks["sampled_rows_decode_back_to_their_source_sequence"].passed


def test_a_duplicate_samples_sidecar_id_is_rejected_even_with_a_matching_hash(
    tmp_path: Path,
) -> None:
    """A digest proves bytes, not that the sidecar key remains a one-to-one join key."""
    path = build(tmp_path, default_records())
    rows = samples_path(path).read_text(encoding="utf-8").splitlines()
    second = rows[2].split(",")
    second[0] = "s1"
    rows[2] = ",".join(second)
    rewrite_samples(path, rows, update_hash=True)

    names = failed(path)
    assert "samples_sidecar_sample_ids_unique" in names
    assert "samples_sidecar_sample_ids_match_table_order" in names


def test_a_missing_samples_sidecar_row_is_rejected_even_with_a_matching_hash(
    tmp_path: Path,
) -> None:
    """A table row without provenance cannot participate in a deterministic audit."""
    path = build(tmp_path, default_records())
    rows = samples_path(path).read_text(encoding="utf-8").splitlines()
    rewrite_samples(path, rows[:-1], update_hash=True)

    names = failed(path)
    assert "samples_sidecar_row_count_matches_manifest" in names
    assert "samples_sidecar_sample_ids_match_table_order" in names


def test_an_extra_samples_sidecar_row_is_rejected_even_with_a_matching_hash(
    tmp_path: Path,
) -> None:
    """Extra provenance silently fans out a downstream join just like a duplicate table key."""
    path = build(tmp_path, default_records())
    rows = samples_path(path).read_text(encoding="utf-8").splitlines()
    extra = rows[-1].split(",")
    extra[0] = "s-extra"
    rows.append(",".join(extra))
    rewrite_samples(path, rows, update_hash=True)

    names = failed(path)
    assert "samples_sidecar_row_count_matches_manifest" in names
    assert "samples_sidecar_sample_ids_match_table_order" in names


def test_reordered_samples_sidecar_ids_are_rejected_even_with_a_matching_hash(
    tmp_path: Path,
) -> None:
    """The contract is ordered equality, so row-position metadata cannot drift between IDs."""
    path = build(tmp_path, default_records())
    rows = samples_path(path).read_text(encoding="utf-8").splitlines()
    rows[1], rows[2] = rows[2], rows[1]
    rewrite_samples(path, rows, update_hash=True)

    assert "samples_sidecar_sample_ids_match_table_order" in failed(path)


def test_a_stale_samples_sidecar_from_another_table_is_rejected(tmp_path: Path) -> None:
    """A coherent old sidecar beside a new table is the partial-replacement failure mode."""
    current = build(tmp_path / "current", default_records())
    stale = build(
        tmp_path / "stale",
        [
            SampleRecord("x1", "ACGT", (0,), "t", {}),
            SampleRecord("x2", "ACGTAC", (1,), "t", {}),
            SampleRecord("x3", "TTTT", (0,), "t", {}),
            SampleRecord("x4", "GGGGCC", (1,), "t", {}),
        ],
    )
    samples_path(current).write_bytes(samples_path(stale).read_bytes())

    names = failed(current)
    assert "samples_sidecar_content_hash_matches" in names
    assert "samples_sidecar_sample_ids_match_table_order" in names


def test_tampering_an_unsampled_metadata_row_is_always_caught_by_the_full_hash(
    tmp_path: Path,
) -> None:
    """Random semantic spot checks cannot probabilistically protect the complete sidecar."""
    path = build(tmp_path, default_records())
    rows = samples_path(path).read_text(encoding="utf-8").splitlines()
    header = rows[0].split(",")
    original_length = header.index("original_length")
    first = rows[1].split(",")
    first[original_length] = "999"
    rows[1] = ",".join(first)
    rewrite_samples(path, rows, update_hash=False)

    result = validate_file(path, sample_size=1, seed=0)
    checks = {check.name: check for check in result.checks}
    assert not checks["samples_sidecar_content_hash_matches"].passed
    assert checks["sampled_rows_decode_back_to_their_source_sequence"].passed
