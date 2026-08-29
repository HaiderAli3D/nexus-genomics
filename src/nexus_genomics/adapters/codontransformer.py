"""CodonTransformer: natural coding sequence versus the model's codon-optimised version.

**What this dataset is.** Each pair is one protein, in one organism, encoded two ways: the
DNA nature actually uses, and the DNA CodonTransformer generates for it. Both members
translate to the *same protein*, so the only thing distinguishing them is which synonymous
codon was chosen at each position. That makes this a **model-artefact detection** task.

**What it is not.** It is not a natural-versus-engineered benchmark, and nothing here should
be presented as one. The positive class is the output of exactly one generator, restricted
to the 17 organisms that model was fine-tuned on. It will not transfer to GeneArt, IDT,
JCat, DNAWorks or Twist output, and it has no bearing on biosecurity screening.

**No model is run and no weights are downloaded.** Zenodo record 13262517 already publishes
``natural_dna``, ``base_dna`` and ``finetune_dna`` side by side in the same row for 622,915
proteins, so the matched pairs are precomputed. This also sidesteps a real hazard in the
package: ``predict_dna_sequence`` calls ``tokenize`` without passing ``max_len``, so the
2048-token ceiling is unreachable through the public API and long proteins would be
silently truncated.

**Do not join the two Zenodo records.** 13262517 has 622,915 rows; the training corpus
(12509224) has 1,001,197 rows in total but only 470,375 for these same 17 organisms. They are
not row-aligned. Pairs are taken from *within* one row of 13262517, which is safe by
construction.

**Sampling is honest about its bias.** The file is 6.08 GB and Zenodo honours HTTP Range, so
this reads a handful of megabyte windows at fixed, evenly spaced offsets rather than
downloading the whole thing. Rows are grouped by organism, so a byte prefix would be one
organism only; strided offsets spread the sample across the file. It is still not a uniform
random sample, and the exact offsets are recorded in the sidecar so the bias is inspectable
rather than hidden.
"""

from __future__ import annotations

import csv
import urllib.request
from collections.abc import Mapping
from pathlib import Path
from typing import Any, Final

from nexus_genomics.adapters.base import Availability, LoadedSource
from nexus_genomics.cleaning import Alphabet
from nexus_genomics.common import SampleRecord, blake2b_256, hash_file

__all__ = ["CodonTransformerAdapter"]

RECORD_URL: Final = "https://zenodo.org/records/13262517"
FILE_URL: Final = (
    "https://zenodo.org/api/records/13262517/files/dataset_with_predictions_metrics.csv/content"
)
LOCAL_NAME: Final = "dataset_with_predictions_metrics.csv"
PUBLISHED_MD5: Final = "9d6571275fa7818f0a2eca04a02fed19"
PUBLISHED_SIZE: Final = 6_083_009_710
TOTAL_ROWS: Final = 622_915

FINE_TUNE_NOTE: Final = (
    "Class 1 is restricted to the organisms CodonTransformer was fine-tuned on, so this is "
    "not a 164-organism sample."
)


class CodonTransformerAdapter:
    """Matched natural / model-generated coding sequences, streamed by HTTP Range."""

    name = "codontransformer_natural_vs_optimized"
    source_url = RECORD_URL
    licence = "CC-BY-4.0 (Zenodo 13262517). Redistribution permitted with attribution."
    target_description = (
        "0 = the natural coding sequence from NCBI; 1 = the CodonTransformer-generated "
        "codon-optimised sequence for the SAME protein in the SAME organism. Model-artefact "
        "detection, not a general natural-versus-engineered benchmark."
    )

    def availability(self, raw_dir: Path) -> Availability:
        local = raw_dir / LOCAL_NAME
        if local.exists():
            return Availability(True, f"local copy present ({local.stat().st_size:,} bytes)")
        return Availability(
            True,
            "no local copy; the adapter reads bounded HTTP Range windows from Zenodo "
            "(a few MB, open access, no credentials)",
            manual_steps=(
                "Optional. For the full 622,915-pair dataset rather than a sampled demo, "
                f"download {FILE_URL} (6.08 GB, md5 {PUBLISHED_MD5}) to "
                f"data/raw/codontransformer/{LOCAL_NAME} and re-run with "
                "--profile full. The demo needs no download."
            ),
            expected_files=(LOCAL_NAME,),
        )

    def load(self, raw_dir: Path, options: Mapping[str, Any]) -> LoadedSource:
        max_pairs = int(options.get("max_pairs", 1000))
        n_windows = int(options.get("n_windows", 8))
        window_bytes = int(options.get("window_bytes", 1_000_000))
        generated = str(options.get("generated_column", "finetune_dna"))
        if generated not in {"finetune_dna", "base_dna"}:
            raise ValueError(
                f"generated_column must be 'finetune_dna' or 'base_dna', got {generated!r}"
            )

        local = raw_dir / LOCAL_NAME
        if local.exists():
            header, rows, provenance = _read_local(local, max_pairs)
            raw_hashes = {LOCAL_NAME: hash_file(local)}
        else:
            header, rows, provenance = _read_ranges(n_windows, window_bytes, max_pairs)
            raw_hashes = {
                f"{FILE_URL} (published md5)": PUBLISHED_MD5,
                "bytes actually read (blake2b-256)": provenance["read_digest"],
            }

        for required in ("protein", "organism", "GeneID", "natural_dna", generated):
            if required not in header:
                raise ValueError(
                    f"the source header does not contain {required!r}. Observed columns: "
                    f"{header}. Refusing to guess a schema."
                )

        records: list[SampleRecord] = []
        length_mismatch = 0
        skipped_missing = 0
        seen_ids: dict[str, int] = {}
        for row in rows:
            natural = row["natural_dna"].strip()
            optimised = row[generated].strip()
            if not natural or not optimised:
                skipped_missing += 1
                continue
            if len(natural) != len(optimised):
                # Both encode the same protein, so the codon counts must match. A mismatch
                # means the row is not a clean pair; counted and dropped rather than used.
                length_mismatch += 1
                continue
            # Derived from the source row's own content, never from its position in the
            # sample. The demo profile reads 8 Range windows and the full profile reads 40,
            # so an enumeration index names a different gene in each -- and concatenating
            # the two files would produce hundreds of duplicate primary keys pointing at
            # different sequences, with nothing to catch it (`require_unique` only checks
            # within one file).
            stem = blake2b_256(
                "\x1f".join(
                    (row["organism"], str(row["GeneID"]), row.get("description", ""), natural)
                ).encode("utf-8")
            )[:12]
            repeat = seen_ids.get(stem, 0)
            seen_ids[stem] = repeat + 1
            pair_id = f"ct_{stem}" if repeat == 0 else f"ct_{stem}_{repeat}"
            shared = {
                "pair_id": pair_id,
                "organism": row["organism"],
                "gene_id": row["GeneID"],
                "protein_length": len(row["protein"]),
                "description": row.get("description", "")[:120],
            }
            records.append(
                SampleRecord(
                    sample_id=f"{pair_id}__natural",
                    sequence=natural,
                    label=(0,),
                    source=self.name,
                    metadata={**shared, "member": "natural"},
                )
            )
            records.append(
                SampleRecord(
                    sample_id=f"{pair_id}__{generated.removesuffix('_dna')}",
                    sequence=optimised,
                    label=(1,),
                    source=self.name,
                    metadata={**shared, "member": generated.removesuffix("_dna")},
                )
            )

        return LoadedSource(
            records=records,
            n_targets=1,
            label_semantics={
                0: "natural coding sequence (NCBI RefSeq CDS)",
                1: f"CodonTransformer-generated codon-optimised sequence ({generated})",
            },
            alphabet=Alphabet.DNA_NUCLEOTIDE,
            manifest_extra={
                "pairing": (
                    "Matched on protein AND organism by construction: both members come "
                    "from the same row of the source file, so they encode the same protein "
                    "for the same organism and differ only in synonymous codon choice."
                ),
                "generated_column": generated,
                "n_pairs": len(records) // 2,
                "source_total_rows": TOTAL_ROWS,
                "sampling": provenance,
                "dropped_length_mismatch_pairs": length_mismatch,
                "dropped_missing_sequence_rows": skipped_missing,
                "grouping_key_for_splits": "gene_id",
                "row_order_predicts_the_class": (
                    "Rows strictly alternate: every even row is natural and every odd row is "
                    "the generated partner of the row above it. So row POSITION predicts the "
                    "target perfectly. Retaining the row index as a feature (a bare "
                    "reset_index() before dropping columns will do it), or using "
                    "KFold(shuffle=False), or taking df.head(n) with odd n, all leak or skew. "
                    "Shuffle with a fixed seed and split grouped by pair_id and gene_id."
                ),
                "sample_id_encodes_the_class": (
                    "sample_id ends in '__natural' or '__finetune', and the .samples.csv "
                    "sidecar carries a 'member' column saying the same thing. Both are "
                    "deliberate -- provenance in the primary key is worth having -- but "
                    "sample_id is a PRIMARY KEY and the sidecar is not the model table. "
                    "Feeding either to a model as a feature would leak the label perfectly. "
                    "The feature block is the position_* columns and nothing else."
                ),
                "leakage_warning": (
                    "Rows in the source are per-transcript, not per-gene: one GeneID recurs "
                    "with near-identical isoform sequences. Any train/test split must group "
                    "by gene_id (and by pair_id) or the two classes of one pair, and the "
                    "isoforms of one gene, will straddle the split."
                ),
                "scope_warning": FINE_TUNE_NOTE,
                "not_a_general_benchmark": (
                    "The positive class is one generator's output. This measures whether a "
                    "model can recognise CodonTransformer, not whether a sequence was "
                    "engineered."
                ),
            },
            raw_input_hashes=raw_hashes,
            sample_metadata_fields=(
                "pair_id",
                "member",
                "organism",
                "gene_id",
                "protein_length",
                "description",
            ),
        )


def _parse_chunk(text: str, header: list[str], *, drop_first: bool) -> list[dict[str, str]]:
    """Parse whole CSV rows out of a byte window, discarding the partial lines at each end.

    The first line of a mid-file window starts mid-record and the last line ends mid-record;
    both are dropped rather than padded, because a half-row silently becomes a wrong row.
    """
    lines = text.split("\n")
    if drop_first:
        lines = lines[1:]
    lines = lines[:-1]
    out: list[dict[str, str]] = []
    malformed = 0
    for row in csv.reader(lines):
        if len(row) == len(header):
            out.append(dict(zip(header, row, strict=True)))
        elif row:
            malformed += 1
    _LAST_CHUNK_STATS["malformed_lines"] = malformed
    return out


_LAST_CHUNK_STATS: dict[str, int] = {"malformed_lines": 0}
"""Field-count mismatches from the most recent chunk, surfaced in the sidecar.

A line whose field count disagrees with the header is a record split across the window
boundary or one containing an embedded newline. Dropping it is correct -- half a row would
become a wrong row -- but dropping it silently is not, because a badly chosen
``window_bytes`` would then look like a sparse region of the file rather than a
misconfiguration.
"""


def _read_ranges(
    n_windows: int, window_bytes: int, max_pairs: int
) -> tuple[list[str], list[dict[str, str]], dict[str, Any]]:
    """Read the header plus evenly spaced windows, without downloading 6.08 GB."""
    head = _range_bytes(0, 65_536)
    header_line = head.decode("utf-8", "replace").split("\n", 1)[0]
    header = next(csv.reader([header_line]))

    offsets = [round(i * PUBLISHED_SIZE / (n_windows + 1)) for i in range(1, n_windows + 1)]
    rows: list[dict[str, str]] = []
    digest_material = bytearray()
    used: list[dict[str, int]] = []
    empty_windows = 0
    malformed = 0
    for offset in offsets:
        if len(rows) >= max_pairs:
            break
        chunk = _range_bytes(offset, window_bytes)
        digest_material.extend(chunk)
        recovered = _parse_chunk(chunk.decode("utf-8", "replace"), header, drop_first=True)
        malformed += _LAST_CHUNK_STATS["malformed_lines"]
        used.append({"start": offset, "length": len(chunk), "rows": len(recovered)})
        if not recovered:
            # A window narrower than a single record yields no complete row. The source has
            # rows up to roughly 324 kB (three sequences of ~108 kB), so this is reachable
            # with a small window_bytes -- and it would otherwise look like a sparse file.
            empty_windows += 1
        rows.extend(recovered)

    if empty_windows == len(used) and used:
        raise RuntimeError(
            f"every one of the {len(used)} Range windows of {window_bytes:,} bytes yielded "
            f"no complete row. Rows in this file reach roughly 324 kB, so window_bytes is "
            f"too small; raise it rather than shipping an empty sample."
        )

    provenance = {
        "method": "http_range_strided",
        "explanation": (
            "Rows are grouped by organism in the source, so a byte prefix would sample one "
            "organism. These offsets are evenly spaced and fixed, so the sample is "
            "reproducible -- but it is NOT a uniform random sample of the 622,915 rows."
        ),
        "windows": used,
        "windows_yielding_no_complete_row": empty_windows,
        "lines_dropped_for_field_count_mismatch": malformed,
        "read_digest": blake2b_256(bytes(digest_material)),
        "rows_recovered": len(rows),
    }
    return header, rows[:max_pairs], provenance


def _range_bytes(start: int, length: int) -> bytes:
    request = urllib.request.Request(
        FILE_URL,
        headers={
            "Range": f"bytes={start}-{start + length - 1}",
            "User-Agent": "nexus-genomics",
        },
    )
    with urllib.request.urlopen(request, timeout=240) as response:
        if response.status != 206:
            raise RuntimeError(
                f"Zenodo returned HTTP {response.status} for a Range request; this adapter "
                f"depends on Range support to avoid a 6.08 GB download."
            )
        return bytes(response.read())


def _read_local(
    path: Path, max_pairs: int
) -> tuple[list[str], list[dict[str, str]], dict[str, Any]]:
    """Stream a local full copy, stopping once enough pairs are collected."""
    rows: list[dict[str, str]] = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        header = list(reader.fieldnames or [])
        for row in reader:
            rows.append({k: (v or "") for k, v in row.items()})
            if len(rows) >= max_pairs:
                break
    return (
        header,
        rows,
        {
            "method": "local_file_prefix",
            "explanation": (
                "Read from the head of a local copy. Rows are grouped by organism, so a "
                "prefix is organism-biased; use the streaming sampler or raise max_pairs."
            ),
            "rows_recovered": len(rows),
        },
    )
