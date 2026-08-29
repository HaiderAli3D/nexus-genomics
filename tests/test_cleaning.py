"""Cleaning rules, and the refusals that keep them honest."""

from __future__ import annotations

import pytest

from nexus_genomics.cleaning import (
    Alphabet,
    CleaningReport,
    CleanResult,
    Quarantined,
    clean_sequence,
)


def clean(
    raw: str, alphabet: Alphabet = Alphabet.DNA_NUCLEOTIDE, *, rna_to_dna: bool = False
) -> tuple[CleanResult | Quarantined, CleaningReport]:
    report = CleaningReport()
    return clean_sequence(raw, alphabet, report=report, rna_to_dna=rna_to_dna), report


def test_whitespace_and_case_are_normalised_and_counted() -> None:
    """FASTA arrives wrapped at 60 or 70 columns and MIDOE ships lowercase bases.

    Both are cosmetic, but a cleaning step that is not counted cannot be audited, so the
    report must show that something happened.
    """
    result, report = clean("ac gt\nAC\tGT\r\n")
    assert isinstance(result, CleanResult)
    assert result.sequence == "ACGTACGT"
    assert result.original_length == 8
    assert report.whitespace_stripped == 1
    assert report.uppercased == 1


def test_an_invalid_character_is_quarantined_not_coerced() -> None:
    """Rewriting an unexpected base to N would look like real ambiguity in the output."""
    result, report = clean("ACGTXZ")
    assert isinstance(result, Quarantined)
    assert result.reason == "invalid_characters"
    assert "'X'" in result.detail or "X" in result.detail
    assert report.quarantined == 1


def test_a_u_in_a_dna_source_is_refused_rather_than_guessed() -> None:
    """A U means the source is RNA or the alphabet is wrong. Converting hides which.

    The error text has to say so, because the fix differs: one is a config change, the other
    is a bug in the adapter.
    """
    result, _ = clean("ACGU")
    assert isinstance(result, Quarantined)
    assert "rna_to_dna" in result.detail


def test_u_to_t_is_applied_only_when_the_source_is_declared_rna() -> None:
    """Opt-in per source, never inferred from the data."""
    result, report = clean("ACGU", rna_to_dna=True)
    assert isinstance(result, CleanResult)
    assert result.sequence == "ACGT"
    assert report.u_to_t_converted == 1


def test_u_to_t_is_refused_outright_for_a_protein_source() -> None:
    """U is selenocysteine in a protein. Converting it would rewrite a real residue.

    This is the sharpest case of the nucleotide/amino-acid confusion the brief warns about:
    the conversion is correct for one alphabet and silently destructive for the other.
    """
    with pytest.raises(ValueError, match="selenocysteine"):
        clean("MACU", Alphabet.AMINO_ACID, rna_to_dna=True)


def test_selenocysteine_and_pyrrolysine_are_valid_residues() -> None:
    """U and O are real amino acids, so a protein containing them must not be quarantined."""
    result, _ = clean("MACUO", Alphabet.AMINO_ACID)
    assert isinstance(result, CleanResult)
    assert result.sequence == "MACUO"


def test_iupac_ambiguity_codes_are_kept_and_counted_not_flattened_to_n() -> None:
    """R means purine and N means any base. Collapsing R to N discards real information.

    Each keeps its own cipher value, which is also what makes the decode round-trip exact.
    """
    result, report = clean("ACGTNRY")
    assert isinstance(result, CleanResult)
    assert result.sequence == "ACGTNRY"
    assert report.ambiguity_counts == {"N": 1, "R": 1, "Y": 1}


def test_a_protein_terminator_is_stripped_only_from_the_end() -> None:
    """CodonTransformer writes `_` and UniProt writes `*` as a stop marker.

    Stripping one from the middle would splice two open reading frames into one sequence.
    """
    result, report = clean("MACD*", Alphabet.AMINO_ACID)
    assert isinstance(result, CleanResult)
    assert result.sequence == "MACD"
    assert report.terminator_stripped == 1

    interior, _ = clean("MA*CD", Alphabet.AMINO_ACID)
    assert isinstance(interior, Quarantined)


def test_an_empty_sequence_is_quarantined_rather_than_padded() -> None:
    """An all-padding row is a well-formed row that means nothing.

    Left in, it trains the model on a sample with no content but a real label.
    """
    result, report = clean("   \n\t ")
    assert isinstance(result, Quarantined)
    assert result.reason == "empty_after_cleaning"
    assert report.quarantine_reasons["empty_after_cleaning"] == 1


def test_a_protein_declared_as_dna_is_caught_by_the_alphabet() -> None:
    """E, Q and I are not bases, so this direction of the mix-up is detectable."""
    result, _ = clean("MEEIQ", Alphabet.DNA_NUCLEOTIDE)
    assert isinstance(result, Quarantined)
    assert result.reason == "invalid_characters"


def test_dna_declared_as_protein_is_undetectable_by_alphabet_so_it_is_counted() -> None:
    """The asymmetry that matters, asserted so nobody assumes the check is symmetric.

    Every nucleotide letter is also a valid amino-acid letter -- A alanine, C cysteine,
    G glycine, T threonine, N asparagine -- so a DNA source mislabelled `amino_acid` passes
    cleaning untouched and no letter-based check can catch it. Refusing on a heuristic could
    reject a genuine peptide, so the counter carries the suspicion into the sidecar instead.
    """
    result, report = clean("ACGTACGTACGT", Alphabet.AMINO_ACID)
    assert isinstance(result, CleanResult)
    assert report.protein_sequences_using_only_nucleotide_letters == 1
    assert report.to_json_dict()["protein_sequences_using_only_nucleotide_letters"] == 1

    real_protein, report2 = clean("MEEIQKACGT", Alphabet.AMINO_ACID)
    assert isinstance(real_protein, CleanResult)
    assert report2.protein_sequences_using_only_nucleotide_letters == 0
