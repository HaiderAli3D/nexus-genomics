"""FELIX engineering signatures: event-relative host context versus introduced sequence.

This narrower recorded-event contrast comes from the ACS Supporting Information of the
GUARDIAN paper, which had to be fetched by hand through a Cloudflare challenge;
:mod:`nexus_genomics.adapters.felix_guardian` explains why nothing in the public MIDOE
repository can supply a general sequence-provenance target.

The ZIP holds three files, whose real schemas are recorded here because nobody had seen them
before this adapter was written:

``guardian_felix_t_and_e_batch_4_signature_sequences.fa``
    335 FASTA records, 562,114 bp. IDs look like
    ``LBNL_IF00015_AmpR_861bp_plasmid.element_Staphylococcus.aureus``: a lab prefix, an
    ``IF#####`` signature number, the element name, its length, **its type**, and the
    organism the element came from. Records suffixed ``_flank1`` / ``_flank2`` are the
    sequence on either side of the change, which the paper describes explicitly.

``guardian_felix_t_and_e_batch_4_sample_signatures.csv``
    ``Sample_ID,Master_Sample_ID,IF_Element`` -- 1,004 rows tying 75 samples to 340 elements.

``guardian_felix_t_and_e_batch_4_sample_metadata.xlsx``
    100 samples with ``Engineering_Key`` (yes 75 / no 25), ``Insert_Over_10_bp``
    (yes 60 / no 40 -- exactly the paper's Table 1), ``Organism``, ``BioSample``, ``SRA``,
    ``Sample_Type`` and ``Fraction_Engineered``. Read for context and recorded in the
    manifest; it is *sample*-level, so it cannot label a sequence.

**The label is event-relative.** A signature's operation type supports only what this
recorded event introduced, removed, or bordered:

* For an **insertion**, ``plasmid.element``, or ``plasmid``, the sequence is cargo the
  recorded event introduced. Class 1.
* For a **deletion**, the sequence was present before this recorded event and removed by it.
  That does not establish whether an earlier event introduced the sequence. Class 0.
* A **flank** is host-locus context bordering the recorded change. Its earlier provenance is
  likewise not established. Class 0.
* ``insertion.site`` names a locus rather than establishing sequence provenance. The source does
  not show that the supplied sequence was introduced by this event or was pre-existing host
  context, so these records are excluded.

So: ``1`` = sequence introduced by this recorded event, ``0`` = host-context sequence not
introduced by this recorded event (a removed region or flank). The two classes come from the
same source collection, but they do not establish a global natural-versus-engineered benchmark.

**Rearrangements are excluded, not guessed at.** For an inversion, a reassortment, a
translocation or a point mutation, the operation name alone does not place the sequence in
either event-relative class. Those records are dropped and counted.

**Length leaks the label.** Introduced elements are generally much longer than host-context
sequence. A single length threshold is therefore a strong baseline, and because the encoder pads
to a fixed width, length is plainly visible as the position where padding starts. This is a
property of the source, not of the conversion; the current baseline is measured from every loaded
profile and recorded in its manifest rather than copied from a stale documentation constant.
"""

from __future__ import annotations

import csv
import io
import re
import zipfile
from collections.abc import Mapping
from pathlib import Path
from typing import Any, Final

from nexus_genomics.adapters.base import Availability, GatedSourceUnavailableError, LoadedSource
from nexus_genomics.cleaning import Alphabet
from nexus_genomics.common import SampleRecord, hash_file

__all__ = ["SI_ZIP_NAME", "FelixSignaturesAdapter"]

SI_ZIP_NAME: Final = "sb3c00398_si_001.zip"
_FASTA: Final = "guardian_felix_t_and_e_batch_4_signature_sequences.fa"
_SIGNATURES: Final = "guardian_felix_t_and_e_batch_4_sample_signatures.csv"
_METADATA: Final = "guardian_felix_t_and_e_batch_4_sample_metadata.xlsx"

INTRODUCED_TYPES: Final = frozenset(
    {"insertion", "compound.insertion", "plasmid.element", "plasmid"}
)
"""Element types whose sequence this recorded event introduced. Class 1."""

REMOVED_HOST_CONTEXT_TYPES: Final = frozenset({"deletion", "compound.deletion"})
"""Element types whose sequence this recorded event removed. Class 0."""

_UNSUPPORTED_TYPE_REASONS: Final = {
    "insertion.site": (
        "insertion.site names a host locus rather than establishing sequence provenance; the "
        "source does not show whether the supplied sequence was introduced by this recorded "
        "event or was pre-existing host context"
    )
}

_KNOWN_TYPE_TOKENS: Final = (
    INTRODUCED_TYPES | REMOVED_HOST_CONTEXT_TYPES | frozenset(_UNSUPPORTED_TYPE_REASONS)
)
"""Every operation-type token the label logic can recognise, for the multi-token scan below."""

_NUCLEOTIDE_LINE: Final = re.compile(r"^[ACGTNRYSWKMBDHV]+$", re.IGNORECASE)
_FLANK: Final = re.compile(r"_flank\d$", re.IGNORECASE)
_LENGTH_TOKEN: Final = re.compile(r"^\d+bp$", re.IGNORECASE)

_STEPS: Final = f"""\
This narrower recorded-event sequence contrast is FREE, but sits behind a Cloudflare bot challenge
(HTTP 403 with 'Cf-Mitigated: challenge' to any script), so it needs a browser:

  1. Open https://pubs.acs.org/doi/10.1021/acssynbio.3c00398?goto=supporting-info
     and let the Cloudflare interstitial clear.
  2. Under 'Supporting Information', click the link rendered as 'sifile1 - zip file'.
  3. Save it as data/raw/felix_signatures/{SI_ZIP_NAME}   (about 175 KB)\
"""


class FelixSignaturesAdapter:
    """Recorded-event host context versus introduced sequence from FELIX signatures."""

    name = "felix_natural_vs_engineered"
    source_url = "https://doi.org/10.1021/acssynbio.3c00398"
    licence = (
        "ACS Supporting Information, free to access. (c) 2024 American Chemical Society; the "
        "underlying files were produced by the IARPA FELIX test-and-evaluation team. The "
        "Acknowledgments grant reproduction rights to the U.S. Government only, so "
        "REDISTRIBUTION BY THIRD PARTIES IS NOT ESTABLISHED -- keep derived outputs local "
        "and cite Adler et al., ACS Synth. Biol. 2024, 13, 1105-1115."
    )
    target_description = (
        "0 = host-context sequence not introduced by this recorded event (a flank or a "
        "region this event removed); 1 = sequence introduced by this recorded event. This "
        "event-relative distinction does not establish global sequence origin."
    )

    def availability(self, raw_dir: Path) -> Availability:
        zip_path = raw_dir / SI_ZIP_NAME
        if zip_path.exists():
            return Availability(True, f"{SI_ZIP_NAME} present ({zip_path.stat().st_size:,} bytes)")
        return Availability(
            False,
            f"missing {SI_ZIP_NAME} in {raw_dir}",
            manual_steps=_STEPS,
            expected_files=(SI_ZIP_NAME,),
        )

    def load(self, raw_dir: Path, options: Mapping[str, Any]) -> LoadedSource:
        status = self.availability(raw_dir)
        if not status.ready:
            raise GatedSourceUnavailableError(f"{status.detail}\n\n{status.manual_steps}")
        zip_path = raw_dir / SI_ZIP_NAME

        with zipfile.ZipFile(zip_path) as archive:
            present = set(archive.namelist())
            missing = [n for n in (_FASTA, _SIGNATURES) if n not in present]
            if missing:
                raise ValueError(
                    f"{SI_ZIP_NAME} does not contain {missing}. Observed members: "
                    f"{sorted(present)}. Refusing to guess which file is which."
                )
            sequences, malformed = _parse_fasta(archive.read(_FASTA).decode("utf-8", "replace"))
            element_samples = _read_signature_map(archive.read(_SIGNATURES).decode("utf-8"))
            sample_context = (
                _read_sample_metadata(archive.read(_METADATA)) if _METADATA in present else {}
            )

        records: list[SampleRecord] = []
        excluded: dict[str, int] = {}
        for element, sequence in sorted(sequences.items()):
            kind = _element_type(element)
            is_flank = bool(_FLANK.search(element))
            if kind in _UNSUPPORTED_TYPE_REASONS:
                excluded[kind] = excluded.get(kind, 0) + 1
                continue
            if is_flank:
                label = 0
                basis = "flank -- host context not introduced by this recorded event"
            elif kind in INTRODUCED_TYPES:
                label, basis = 1, f"{kind} -- sequence introduced by this recorded event"
            elif kind in REMOVED_HOST_CONTEXT_TYPES:
                label = 0
                basis = f"{kind} -- host context removed, not introduced, by this recorded event"
            else:
                # A rearrangement or point mutation does not identify a sequence as cargo or
                # host context under this recorded event, so neither class is supported.
                excluded[kind] = excluded.get(kind, 0) + 1
                continue

            samples = element_samples.get(element, set())
            records.append(
                SampleRecord(
                    sample_id=element,
                    sequence=sequence,
                    label=(label,),
                    source=self.name,
                    metadata={
                        "if_element": element,
                        "element_type": kind,
                        "is_flank": int(is_flank),
                        "label_basis": basis,
                        "element_origin_organism": _origin_organism(element),
                        "n_samples_containing": len(samples),
                        "samples_containing": "|".join(sorted(samples)[:8]),
                    },
                )
            )

        limit = options.get("max_records")
        if limit:
            # Sorted by element id above, so a cap is deterministic -- but it would also skew
            # the class balance, because ids cluster by type. Take a balanced interleave.
            records = _balanced_head(records, int(limit))

        counts = {
            0: sum(1 for r in records if r.label == (0,)),
            1: sum(1 for r in records if r.label == (1,)),
        }
        return LoadedSource(
            records=records,
            n_targets=1,
            label_semantics={
                0: (
                    "host-context sequence not introduced by this recorded event -- either "
                    "a flank beside the change or a region this event removed; earlier "
                    "provenance is not established"
                ),
                1: (
                    "sequence introduced by this recorded event (insertion or plasmid "
                    "element cargo)"
                ),
            },
            alphabet=Alphabet.DNA_NUCLEOTIDE,
            manifest_extra={
                "label_provenance_scope": (
                    "The operation type supports an event-relative distinction: sequence "
                    "introduced by this recorded event versus host context not introduced by "
                    "this recorded event. It does not establish global sequence origin."
                ),
                "why_class_zero_is_event_relative": (
                    "A deletion signature records sequence this event removed, and a flank "
                    "records host-locus context beside this event. Both are labelled 0 because "
                    "this recorded event did not introduce them; earlier provenance remains "
                    "unknown."
                ),
                "class_counts": {str(k): v for k, v in counts.items()},
                "excluded_element_types": excluded,
                "excluded_element_type_reasons": {
                    kind: _UNSUPPORTED_TYPE_REASONS.get(
                        kind,
                        "the operation type does not identify the sequence as introduced cargo "
                        "or host context under this recorded event",
                    )
                    for kind in excluded
                },
                "excluded_reason": (
                    "Unsupported provenance is excluded rather than guessed. Rearrangements "
                    "and point mutations do not identify a sequence-level event-relative "
                    "class; insertion.site names a locus without establishing whether the "
                    "supplied sequence is introduced cargo or pre-existing host context."
                ),
                "length_only_baseline_accuracy": _length_baseline(records),
                "length_confounder_warning": (
                    "Introduced elements are much longer than host-context flanks, so a single "
                    "length threshold already classifies a large share of this dataset. The "
                    "encoder pads to a fixed width, which makes length plainly visible as the "
                    "position where padding starts. A model must beat the baseline above to "
                    "be learning sequence content rather than length."
                ),
                "signatures_with_no_sample_linkage": [
                    r.sample_id for r in records if not r.metadata["n_samples_containing"]
                ],
                "sample_link_note": (
                    "One signature carries no sample linkage, and it is not a join bug. "
                    "Its FASTA header is upper-cased and drops the LBNL_ prefix the other "
                    "334 carry, which looks like a naming mismatch -- but the "
                    "sample-signature CSV lists only that element's _flank1 and _flank2, "
                    "never the element itself. So n_samples_containing = 0 is the truth "
                    "about the source, not an artefact of matching. A case-folded join was "
                    "tried, rescued nothing, and was removed rather than left in as "
                    "machinery that could conflate two genuinely different elements."
                ),
                "length_token_disagreements": _length_token_check(sequences),
                "length_token_note": (
                    "Every non-flank signature name states its own length (the '861bp' in "
                    "LBNL_IF00015_AmpR_861bp_plasmid.element_...). Comparing that against the "
                    "sequence actually parsed is a free integrity check on both the file and "
                    "this parser -- it is what would have caught the malformed records below "
                    "even if nobody had looked. 211 of 222 agree exactly; the disagreements "
                    "listed here are defects in the published file, not in the conversion."
                ),
                "malformed_fasta_records_quarantined": malformed,
                "malformed_reason": (
                    "Two records in the published FASTA lost the '>' from their header line, "
                    "so a naive parser concatenates the header text and the following "
                    "sequence onto the PREVIOUS record -- which is how that record's length "
                    "reaches 202,036 bp. Their names do not appear in the sample-signature "
                    "CSV, so the repair could not be verified against anything and the two "
                    "records were dropped rather than guessed at."
                ),
                "elements_in_csv_without_a_sequence": sorted(set(element_samples) - set(sequences))[
                    :10
                ],
                "sample_level_context": sample_context,
                "sample_level_is_not_the_target": (
                    "Engineering_Key (75 yes / 25 no) and Insert_Over_10_bp (60 yes / 40 no) "
                    "label whole SAMPLES, and the sample sequences are not in this ZIP -- they "
                    "are in SRA and GenBank. They are recorded as context; the target here is "
                    "per-signature."
                ),
                "grouping_key_for_splits": "if_element_number",
                "leakage_warning": (
                    "An element and its _flank1/_flank2 records share an IF number and sit "
                    "beside each other in the genome, and they are in OPPOSITE classes. Split "
                    "grouped by the IF number (the IF##### in sample_id) or a model will see "
                    "one side of a locus in training and the other in test."
                ),
            },
            raw_input_hashes={SI_ZIP_NAME: hash_file(zip_path)},
            sample_metadata_fields=(
                "if_element",
                "element_type",
                "is_flank",
                "label_basis",
                "element_origin_organism",
                "n_samples_containing",
                "samples_containing",
            ),
        )


def _parse_fasta(text: str) -> tuple[dict[str, str], list[str]]:
    """Parse FASTA, refusing to append anything that is not a sequence line.

    The published file contains two header lines that lost their ``>``. Appending them would
    splice a header and its sequence onto the previous record; the record's length is what
    gives the defect away (202,036 bp for a signature whose name says 16 bp). Anything that is
    not a run of nucleotide characters ends the current record and is quarantined by name.
    """
    sequences: dict[str, str] = {}
    malformed: list[str] = []
    current: str | None = None
    buffer: list[str] = []

    def flush() -> None:
        if current is not None:
            sequences[current] = "".join(buffer)

    for raw_line in text.split("\n"):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            flush()
            current, buffer[:] = line[1:].strip(), []
        elif _NUCLEOTIDE_LINE.match(line):
            buffer.append(line)
        else:
            # Not a header and not a sequence. Close the current record so the junk cannot
            # reach it, and record what was seen.
            flush()
            malformed.append(line[:120])
            current, buffer[:] = None, []
    flush()
    return sequences, malformed


def _length_token_check(sequences: Mapping[str, str]) -> dict[str, str]:
    """Compare each signature's self-declared length against the sequence actually parsed.

    A name like ``..._861bp_plasmid.element_...`` states its own length, which makes the file
    self-checking: a parser that concatenated two records would show a large disagreement here.
    Flanks are skipped, because the length token belongs to the element, not to its flank.
    """
    out: dict[str, str] = {}
    for name, sequence in sequences.items():
        if _FLANK.search(name):
            continue
        match = re.search(r"_(\d+)bp_", name, re.IGNORECASE) or re.search(
            r"_(\d+)bp$", name, re.IGNORECASE
        )
        if match and int(match.group(1)) != len(sequence):
            out[name] = f"name says {match.group(1)} bp, sequence is {len(sequence)} bp"
    return out


def _read_signature_map(text: str) -> dict[str, set[str]]:
    """``IF_Element -> {Sample_ID}`` from the sample-signature CSV.

    Matched on the literal name. A case-folded, prefix-stripped fallback was tried, because
    one FASTA header is upper-cased and drops the ``LBNL_`` prefix the other 334 carry -- but
    it rescued nothing: the CSV genuinely does not list that element, only its two flanks.
    Normalising away real distinctions is how two different elements silently become one, so
    the fallback was removed rather than kept for a case that does not arise.
    """
    reader = csv.DictReader(io.StringIO(text))
    expected = {"Sample_ID", "IF_Element"}
    if not expected <= set(reader.fieldnames or []):
        raise ValueError(
            f"{_SIGNATURES} does not have the columns {sorted(expected)}; it has "
            f"{reader.fieldnames}. Refusing to guess a schema."
        )
    out: dict[str, set[str]] = {}
    for row in reader:
        out.setdefault(row["IF_Element"], set()).add(row["Sample_ID"])
    return out


def _read_sample_metadata(blob: bytes) -> dict[str, Any]:
    """Sample-level context. Never a target -- the sample sequences are not in this ZIP."""
    try:
        import openpyxl
    except ImportError:  # pragma: no cover - openpyxl is a declared dependency
        return {"unavailable": "openpyxl is not installed, so the XLSX was not read"}

    sheet = openpyxl.load_workbook(io.BytesIO(blob), read_only=True).worksheets[0]
    rows = list(sheet.values)
    header = [str(h) for h in rows[0]]
    body = [r for r in rows[1:] if any(v is not None for v in r)]

    def tally(column: str) -> dict[str, int]:
        if column not in header:
            return {}
        index = header.index(column)
        counts: dict[str, int] = {}
        for row in body:
            value = row[index]
            if value is not None:
                counts[str(value)] = counts.get(str(value), 0) + 1
        return dict(sorted(counts.items()))

    return {
        "n_samples": len(body),
        "columns": [h for h in header if h != "None"],
        "engineering_key": tally("Engineering_Key"),
        "insert_over_10_bp": tally("Insert_Over_10_bp"),
        "sample_type": tally("Sample_Type"),
        "n_distinct_organisms": len({str(r[header.index("Organism")]) for r in body}),
    }


def _element_type(element: str) -> str:
    """The recognised operation-type token following a signature's ``<n>bp`` length token.

    ``LBNL_IF00015_AmpR_861bp_plasmid.element_Staphylococcus.aureus`` -> ``plasmid.element``.
    That is the common case, but some names carry an extra descriptor between the length and
    the type -- ``lacZalpha_168bp_CDSpartial_insertion_...``, or an inline mutation/indel
    description such as ``voz1_9bp_607-615delACGACCTCC_DDLQE29-33EE_deletion_...``. Returning
    the token immediately after the length unconditionally would misparse those as an
    unrecognised type and drop a genuinely classifiable record into "excluded" -- confirmed
    against the published file for 6 of 335 elements. Every token after the length marker is
    checked in order against the known type vocabulary, and the first exact match wins; only
    when nothing after the length token is recognised does this fall back to the
    immediately-following token, so a true rearrangement or point mutation (e.g. ``inversion``,
    ``reassortment``) is still reported under its own name for the excluded-type tally.
    """
    core = _FLANK.sub("", element)
    parts = core.split("_")
    for index, part in enumerate(parts):
        if _LENGTH_TOKEN.match(part):
            tail = parts[index + 1 :]
            for token in tail:
                if token.lower() in _KNOWN_TYPE_TOKENS:
                    return token.lower()
            return tail[0].lower() if tail else "unparsed"
    return "unparsed"


def _origin_organism(element: str) -> str:
    """The trailing organism token, where the name carries one."""
    parts = _FLANK.sub("", element).split("_")
    tail = parts[-1] if parts else ""
    return tail.replace(".", " ") if tail and tail[0].isupper() else ""


def _balanced_head(records: list[SampleRecord], limit: int) -> list[SampleRecord]:
    """Take ``limit`` records while keeping both classes represented.

    A plain head would skew the balance badly, because element ids cluster by type: the first
    N sorted ids are dominated by whichever type happens to sort first. Interleaving the two
    classes keeps the demo file's balance close to the full file's, deterministically.
    """
    host_context = [r for r in records if r.label == (0,)]
    introduced = [r for r in records if r.label == (1,)]
    out: list[SampleRecord] = []
    for index in range(max(len(host_context), len(introduced))):
        for bucket in (host_context, introduced):
            if index < len(bucket) and len(out) < limit:
                out.append(bucket[index])
        if len(out) >= limit:
            break
    return out


def _length_baseline(records: list[SampleRecord]) -> float:
    """Best accuracy obtainable from sequence length alone, rounded.

    Published so nobody mistakes a length detector for an engineering detector.
    """
    host_context = sorted(len(r.sequence) for r in records if r.label == (0,))
    introduced = sorted(len(r.sequence) for r in records if r.label == (1,))
    if not host_context or not introduced:
        return 0.0
    total = len(host_context) + len(introduced)
    thresholds = sorted({*host_context, *introduced})
    best = max(
        (
            sum(1 for length in introduced if length >= threshold)
            + sum(1 for length in host_context if length < threshold)
        )
        / total
        for threshold in thresholds
    )
    return round(best, 4)
