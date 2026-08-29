"""The cipher and the length policies.

Each test names the trap it exists to stop, following the house style in
``C:\\Projects\\qecgen\\tests``: a test whose docstring only restates its assertion tells a
future reader nothing about why deleting it would be dangerous.
"""

from __future__ import annotations

import pytest

from nexus_genomics.encoding import (
    LETTER_TO_TOKEN,
    MAX_TOKEN,
    PAD_TOKEN,
    EncodingMode,
    LengthPolicy,
    decode,
    encode,
    place,
)


def test_the_cipher_matches_the_brief_worked_example() -> None:
    """`1,3,7,20` must decode to ACGT.

    The brief gives this row as its only concrete specification of the encoding, and it is
    also the evidence that the intended representation is nucleotides rather than amino
    acids -- `A,C,G,T` is valid DNA and no plausible peptide. If this breaks, the whole
    representation decision behind the project is no longer supported by its own example.
    """
    assert encode("ACGT") == [1, 3, 7, 20]
    assert decode([1, 3, 7, 20]) == "ACGT"
    assert decode([1, 20, 7, 20]) == "ATGT"


def test_zero_is_never_emitted_for_a_biological_character() -> None:
    """0 is reserved for padding, so no letter may map to it.

    If any letter encoded to 0, a padding cell and a residue would be indistinguishable and
    every downstream check that counts padding would be silently wrong.
    """
    assert PAD_TOKEN == 0
    assert 0 not in LETTER_TO_TOKEN.values()
    assert min(LETTER_TO_TOKEN.values()) == 1
    assert max(LETTER_TO_TOKEN.values()) == MAX_TOKEN == 26
    assert len(LETTER_TO_TOKEN) == 26


def test_an_unmapped_character_is_refused_not_substituted() -> None:
    """Folding an unexpected character to N (or to 0) would corrupt the sequence quietly.

    The encode step is the last place that can still tell the difference, so it raises
    rather than choosing a value.
    """
    with pytest.raises(ValueError, match="no token"):
        encode("AC-GT")
    with pytest.raises(ValueError, match="no token"):
        encode("acgt")


def test_decoding_a_token_outside_the_cipher_is_refused() -> None:
    """A 27 in a feature column means the file was not written by this encoder."""
    with pytest.raises(ValueError, match="outside 0-26"):
        decode([1, 27])


def test_padding_is_appended_and_counted_never_interleaved() -> None:
    """Padding in the middle of a row would splice a gap into the sequence.

    The validation suite relies on "residues then padding" being a hard invariant to detect
    corruption, so the placer must never produce anything else.
    """
    (placed,) = place(encode("ACGT"), 8, LengthPolicy.PAD_OR_TRUNCATE, 8)
    assert placed.tokens == (1, 3, 7, 20, 0, 0, 0, 0)
    assert placed.pad_count == 4
    assert placed.removed == 0
    nonzero = [i for i, t in enumerate(placed.tokens) if t != 0]
    assert nonzero == list(range(len(nonzero)))


def test_truncation_records_how_much_was_removed() -> None:
    """Truncating without recording the loss is the failure the brief calls out by name."""
    (placed,) = place(encode("ACGTACGT"), 3, LengthPolicy.PAD_OR_TRUNCATE, 3)
    assert placed.tokens == (1, 3, 7)
    assert placed.removed == 5
    assert placed.pad_count == 0


def test_windowing_keeps_the_tail_instead_of_dropping_it() -> None:
    """A final short window must be padded, not discarded.

    Dropping it would silently delete the end of every sequence whose length is not a
    multiple of the stride -- a systematic, invisible truncation across the whole dataset.
    """
    placements = place(encode("ACGTACGTA"), 4, LengthPolicy.WINDOW, 4)
    assert [p.tokens for p in placements] == [
        (1, 3, 7, 20),
        (1, 3, 7, 20),
        (1, 0, 0, 0),
    ]
    assert [(p.start, p.end) for p in placements] == [(1, 4), (5, 8), (9, 9)]
    assert sum(p.removed for p in placements) == 0


def test_window_offsets_are_one_based_and_inclusive() -> None:
    """Sequence coordinates are one-based everywhere in biology.

    A zero-based offset here would be off by one against every GenBank coordinate anyone
    compares this table to, and the error would look like a real biological discrepancy.
    """
    placements = place(encode("A" * 10), 4, LengthPolicy.WINDOW, 4)
    assert placements[0].start == 1
    assert placements[0].end == 4
    assert placements[1].start == 5


def test_an_empty_sequence_still_produces_one_row_under_pad_policy() -> None:
    """Guards the loop condition in `place`: `while offset < length or not out`."""
    placements = place([], 4, LengthPolicy.WINDOW, 4)
    assert len(placements) == 1
    assert placements[0].tokens == (0, 0, 0, 0)


def test_unimplemented_modes_are_named_rather_than_missing() -> None:
    """Codon and amino-acid-translated must be distinguishable from a typo.

    A `KeyError` would tell a reader the mode does not exist; the enum plus the refusal in
    the pipeline tells them it exists, is different from nucleotides, and is not built.
    """
    assert EncodingMode.CODON.value == "codon"
    assert EncodingMode.AMINO_ACID_TRANSLATED.value == "amino_acid_translated"
    assert EncodingMode.LETTER.value == "letter"


def test_a_stride_wider_than_the_window_is_refused() -> None:
    """It would drop the residues between consecutive windows and report removed=0.

    With width=1000 and stride=1500, a 3,000-residue sequence emits two rows covering
    1-1000 and 1501-2500: a third of the sequence appears in no row, and every Placed still
    says `removed=0`, so nothing in the manifest or the sidecar records the loss.
    """
    with pytest.raises(ValueError, match="greater than the window width"):
        place(encode("A" * 3000), 1000, LengthPolicy.WINDOW, 1500)


def test_overlapping_windows_stop_once_the_tail_is_covered() -> None:
    """Every window past the end of the sequence is a strict suffix of an earlier one.

    Emitting them adds no residue, weights the tail of every sequence far above its head,
    and produces final rows that are almost entirely padding.
    """
    placements = place(encode("A" * 1050), 1000, LengthPolicy.WINDOW, 100)
    assert len(placements) == 2
    assert [(p.start, p.end) for p in placements] == [(1, 1000), (101, 1050)]
    assert placements[-1].pad_count == 50


def test_windows_still_cover_the_whole_sequence_when_they_do_not_overlap() -> None:
    """The early stop must not lose the tail it was added to protect."""
    placements = place(encode("A" * 25), 10, LengthPolicy.WINDOW, 10)
    assert [(p.start, p.end) for p in placements] == [(1, 10), (11, 20), (21, 25)]
    covered = set()
    for p in placements:
        covered |= set(range(p.start, p.end + 1))
    assert covered == set(range(1, 26))
