"""Deterministic sequence cleaning, with every operation counted.

The rule this module exists to enforce: **a sequence is either cleaned by a documented rule
or quarantined**. Nothing is silently coerced. Folding an unexpected character to ``N``, or
an empty sequence to a row of padding, produces a file that parses perfectly and means
something else -- and the damage only surfaces much later as a model that trains, converges
and is wrong.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from enum import StrEnum
from typing import Final

__all__ = [
    "AMINO_ACID_AMBIGUITY",
    "AMINO_ACID_CANONICAL",
    "DNA_AMBIGUITY",
    "DNA_CANONICAL",
    "Alphabet",
    "CleanResult",
    "CleaningReport",
    "Quarantined",
    "clean_sequence",
]


class Alphabet(StrEnum):
    """What the letters in a sequence mean.

    This is *not* an encoding mode. The cipher is always A=1..Z=26 regardless; the alphabet
    records which letters are legitimate and how to read them, so that a protein sequence is
    never validated against a nucleotide letter set or vice versa. Mixing the two up is the
    single most damaging error available here, because both are strings of capital letters
    and neither looks wrong.
    """

    DNA_NUCLEOTIDE = "dna_nucleotide"
    AMINO_ACID = "amino_acid"


DNA_CANONICAL: Final = frozenset("ACGT")
DNA_AMBIGUITY: Final = frozenset("RYSWKMBDHVN")
"""IUPAC nucleotide ambiguity codes. Kept, not dropped, and never rewritten to N.

Each takes its own cipher value (N=14, R=18, ...), so the encoding stays reversible and the
distinction between "any base" and "purine" survives into the table. Rewriting them all to
N would discard real information and make the decode check in `validation` pass on a
sequence that is no longer the source sequence.
"""

AMINO_ACID_CANONICAL: Final = frozenset("ACDEFGHIKLMNPQRSTVWY")
AMINO_ACID_AMBIGUITY: Final = frozenset("BZXJUO")
"""B=Asx, Z=Glx, X=any, J=Leu/Ile, U=selenocysteine, O=pyrrolysine.

U and O are real residues, not ambiguity -- which is exactly why `U -> T` conversion must
never be applied to a protein sequence. It would turn every selenocysteine into a threonine
and leave a sequence that still looks like a protein.
"""

_TERMINATORS: Final = frozenset("*_")
"""Stop-codon markers seen in the wild. CodonTransformer writes `_`, UniProt/EMBOSS write `*`.

Stripped only from the END of a protein sequence, and counted. Stripped from the middle
would silently splice two open reading frames together.
"""

_NUCLEOTIDE_LOOKALIKE: Final = frozenset("ACGTN")
"""Letters that are simultaneously valid bases and valid residues.

A=alanine, C=cysteine, G=glycine, T=threonine, N=asparagine. See
:class:`CleaningReport.protein_sequences_using_only_nucleotide_letters`.
"""

_WHITESPACE: Final = frozenset(" \t\r\n\v\f")


@dataclass(slots=True)
class CleaningReport:
    """Counts of every transformation applied across a whole source.

    Serialised into the sidecar. A cleaning step that is not counted here is a cleaning step
    that cannot be audited, so every branch in :func:`clean_sequence` increments something.
    """

    seen: int = 0
    uppercased: int = 0
    whitespace_stripped: int = 0
    terminator_stripped: int = 0
    u_to_t_converted: int = 0
    ambiguity_counts: Counter[str] = field(default_factory=Counter)
    quarantined: int = 0
    quarantine_reasons: Counter[str] = field(default_factory=Counter)
    protein_sequences_using_only_nucleotide_letters: int = 0
    """Sequences declared amino-acid that contain nothing but A, C, G, T and N.

    The alphabet check catches a protein mistakenly declared DNA -- E, Q, I and friends are
    not bases. It CANNOT catch the reverse, because every nucleotide letter is also a valid
    amino-acid letter (A=alanine, C=cysteine, G=glycine, T=threonine, N=asparagine). So a
    DNA source mislabelled `amino_acid` would sail through cleaning untouched.

    Since a heuristic here could be wrong -- a genuine peptide of only those five residues
    is possible, if vanishingly unlikely at length -- this counts rather than refuses. The
    count reaches the sidecar and the validation report, which turns an otherwise invisible
    class of error into one a reader can see.
    """

    def to_json_dict(self) -> dict[str, object]:
        return {
            "seen": self.seen,
            "uppercased": self.uppercased,
            "whitespace_stripped": self.whitespace_stripped,
            "terminator_stripped": self.terminator_stripped,
            "u_to_t_converted": self.u_to_t_converted,
            "ambiguity_counts": dict(sorted(self.ambiguity_counts.items())),
            "quarantined": self.quarantined,
            "quarantine_reasons": dict(sorted(self.quarantine_reasons.items())),
            "protein_sequences_using_only_nucleotide_letters": (
                self.protein_sequences_using_only_nucleotide_letters
            ),
        }


@dataclass(frozen=True, slots=True)
class CleanResult:
    """A sequence that passed cleaning, with what happened to it."""

    sequence: str
    original_length: int
    """Length after whitespace removal, before any truncation. Recorded per sample."""


@dataclass(frozen=True, slots=True)
class Quarantined:
    """A sequence that did not pass, and precisely why."""

    reason: str
    detail: str


def clean_sequence(
    raw: str,
    alphabet: Alphabet,
    *,
    report: CleaningReport,
    rna_to_dna: bool = False,
) -> CleanResult | Quarantined:
    """Apply the documented cleaning rules, or refuse.

    ``rna_to_dna`` is opt-in per source and only legal for a nucleotide alphabet. It is not
    inferred from the presence of a U: a U in a sequence the source documents as DNA is a
    defect in the input or a mislabelled alphabet, and converting it would paper over
    whichever it is.
    """
    report.seen += 1

    if rna_to_dna and alphabet is not Alphabet.DNA_NUCLEOTIDE:
        raise ValueError(
            "rna_to_dna was requested for an amino-acid source. U is selenocysteine in a "
            "protein, not uracil; converting it would rewrite a real residue as threonine."
        )

    stripped = "".join(c for c in raw if c not in _WHITESPACE)
    if len(stripped) != len(raw):
        report.whitespace_stripped += 1

    upper = stripped.upper()
    if upper != stripped:
        report.uppercased += 1

    if alphabet is Alphabet.AMINO_ACID:
        trimmed = upper.rstrip("".join(_TERMINATORS))
        if trimmed != upper:
            report.terminator_stripped += 1
        upper = trimmed

    if not upper:
        report.quarantined += 1
        report.quarantine_reasons["empty_after_cleaning"] += 1
        return Quarantined("empty_after_cleaning", "sequence is empty once whitespace is removed")

    if rna_to_dna and "U" in upper:
        upper = upper.replace("U", "T")
        report.u_to_t_converted += 1

    canonical, ambiguity = (
        (DNA_CANONICAL, DNA_AMBIGUITY)
        if alphabet is Alphabet.DNA_NUCLEOTIDE
        else (AMINO_ACID_CANONICAL, AMINO_ACID_AMBIGUITY)
    )
    allowed = canonical | ambiguity

    invalid = sorted({c for c in upper if c not in allowed})
    if invalid:
        report.quarantined += 1
        reason = "invalid_characters"
        report.quarantine_reasons[reason] += 1
        hint = ""
        if alphabet is Alphabet.DNA_NUCLEOTIDE and "U" in invalid:
            hint = (
                " A U in a source documented as DNA means either the source is RNA (set "
                "rna_to_dna for it) or the alphabet is wrong. It is not converted by "
                "guesswork."
            )
        return Quarantined(
            reason,
            f"characters {invalid!r} are not valid for alphabet {alphabet.value}.{hint}",
        )

    for char in upper:
        if char in ambiguity:
            report.ambiguity_counts[char] += 1

    if alphabet is Alphabet.AMINO_ACID and set(upper) <= _NUCLEOTIDE_LOOKALIKE:
        report.protein_sequences_using_only_nucleotide_letters += 1

    return CleanResult(sequence=upper, original_length=len(upper))
