"""The alphabet cipher, and the length policies built on it.

    A -> 1,  B -> 2,  ...,  Z -> 26,  padding -> 0

**These integers are categorical tokens.** Their numerical distance carries no biological
meaning whatsoever: T=20 is not "ten times" C=3, and G=7 does not sit "between" C and T in
any chemical sense. A model must treat these columns as categories -- an embedding layer, or
one-hot -- and a model that treats them as magnitudes is learning from an artefact of the
Latin alphabet. The warning is repeated in the README and in every ``.ml.manifest.json``,
because this is the single easiest way to misuse these files. It is deliberately noted that
it is NOT in the ``.ml.samples.csv``: handing someone the two CSVs alone does not carry it.

``0`` is reserved for padding and nothing else. No letter maps to it, so a zero in a feature
column always means "no residue here", and :mod:`nexus_genomics.validation` checks that
every zero sits where the sidecar says padding begins.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from typing import Final

__all__ = [
    "LETTER_TO_TOKEN",
    "MAX_TOKEN",
    "PAD_TOKEN",
    "TOKEN_TO_LETTER",
    "EncodingMode",
    "LengthPolicy",
    "Placed",
    "decode",
    "encode",
    "place",
]

PAD_TOKEN: Final = 0
MAX_TOKEN: Final = 26

LETTER_TO_TOKEN: Final[dict[str, int]] = {chr(ord("A") + i): i + 1 for i in range(26)}
TOKEN_TO_LETTER: Final[dict[int, str]] = {v: k for k, v in LETTER_TO_TOKEN.items()}


class EncodingMode(StrEnum):
    """How a sequence becomes columns.

    Only ``LETTER`` is implemented. The other two are named rather than omitted so that a
    configuration asking for them fails with an explanation instead of an unhelpful
    ``KeyError`` -- and so that nobody reads the absence of a codon mode as "codons are the
    same thing as nucleotides".
    """

    LETTER = "letter"
    """One column per letter of the sequence, as written. The default and the only one built.

    For a nucleotide source this is one column per base, which is what the brief's own worked
    example encodes: `1,3,7,20` decodes to A,C,G,T. It preserves synonymous codon choice,
    which a translated representation would destroy, and it needs no reading frame.
    """

    CODON = "codon"
    """Three nucleotides per column, 1..64. Defined, not implemented."""

    AMINO_ACID_TRANSLATED = "amino_acid_translated"
    """Translate nucleotides to residues, then encode. Defined, not implemented, and refused.

    Translation requires knowing where the coding regions are and in which frame. A plasmid
    or a genome is not one open reading frame, so translating it end to end would emit a
    protein that does not exist. A source whose sequences are ALREADY amino acids does not
    need this mode -- it sets `alphabet: amino_acid` and uses LETTER.
    """


class LengthPolicy(StrEnum):
    """What to do when a sequence is not exactly the configured width."""

    PAD_OR_TRUNCATE = "pad_or_truncate"
    WINDOW = "window"


def encode(sequence: str) -> list[int]:
    """Letters to tokens. Raises on anything outside A-Z, never substitutes."""
    try:
        return [LETTER_TO_TOKEN[c] for c in sequence]
    except KeyError as exc:
        raise ValueError(
            f"character {exc.args[0]!r} has no token. The cipher covers A-Z only, and 0 is "
            f"reserved for padding, so there is no value this could be mapped to without "
            f"inventing one. Clean or quarantine the sequence first."
        ) from None


def decode(tokens: list[int] | tuple[int, ...]) -> str:
    """Tokens back to letters, dropping padding. The inverse used by the round-trip check."""
    out: list[str] = []
    for token in tokens:
        if token == PAD_TOKEN:
            continue
        try:
            out.append(TOKEN_TO_LETTER[token])
        except KeyError:
            raise ValueError(
                f"token {token!r} is outside 0-26, so this row was not written by this encoder."
            ) from None
    return "".join(out)


@dataclass(frozen=True, slots=True)
class Placed:
    """One row's tokens, plus what had to happen to make them fit."""

    tokens: tuple[int, ...]
    start: int
    """1-based inclusive offset of the first residue, in the original sequence."""
    end: int
    """1-based inclusive offset of the last real residue. Equals ``start - 1`` if empty."""
    pad_count: int
    removed: int
    """Residues discarded by truncation. Always recorded, never silently zero."""


def place(tokens: list[int], width: int, policy: LengthPolicy, stride: int) -> list[Placed]:
    """Fit a token list into rows of exactly ``width`` columns.

    ``PAD_OR_TRUNCATE`` yields exactly one row and records how much was thrown away.
    ``WINDOW`` yields one row per window and throws nothing away -- the final window is
    padded rather than dropped, because dropping it would silently delete the tail of every
    sequence whose length is not a multiple of the stride.
    """
    if width <= 0:
        raise ValueError(f"width must be positive, got {width}")
    length = len(tokens)

    if policy is LengthPolicy.PAD_OR_TRUNCATE:
        kept = tokens[:width]
        pad = width - len(kept)
        return [
            Placed(
                tokens=tuple(kept) + (PAD_TOKEN,) * pad,
                start=1,
                end=len(kept),
                pad_count=pad,
                removed=length - len(kept),
            )
        ]

    if stride <= 0:
        raise ValueError(f"stride must be positive, got {stride}")
    if stride > width:
        # A stride wider than the window is not a windowing scheme, it is undocumented
        # subsampling: the `stride - width` residues between consecutive windows appear in
        # no row at all, and every Placed would still report `removed=0`, so the loss would
        # be invisible in both the manifest and the per-sample sidecar.
        raise ValueError(
            f"stride ({stride}) is greater than the window width ({width}), which would "
            f"drop the {stride - width} residues between consecutive windows without "
            f"recording them anywhere. Set stride <= width, or use pad_or_truncate if you "
            f"genuinely want only the first {width} residues."
        )
    out: list[Placed] = []
    offset = 0
    while offset < length or not out:
        chunk = tokens[offset : offset + width]
        pad = width - len(chunk)
        out.append(
            Placed(
                tokens=tuple(chunk) + (PAD_TOKEN,) * pad,
                start=offset + 1,
                end=offset + len(chunk),
                pad_count=pad,
                removed=0,
            )
        )
        # Stop as soon as a window has reached the end. With stride < width every further
        # window is a strict suffix of this one -- it contains no residue this one does not
        # -- so continuing would emit progressively emptier duplicate rows and weight the
        # tail of every sequence far above its head.
        if offset + width >= length:
            break
        offset += stride
    return out
