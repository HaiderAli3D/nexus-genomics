"""The writer: one header row, one row per sample, metadata in a JSON sidecar.

The layout and every convention here is copied from ``qecgen/qecgen/exporters/ml_csv.py``,
which was built for a consumer that asked for ``Primary Key, Target, Variable_1, ...`` --
the same phrasing the Nexus brief uses. Three files share a stem::

    d.ml.csv              one header row, then one row per sample. No '#' line anywhere.
    d.ml.manifest.json    the manifest, plus this file's literal column names
    d.ml.samples.csv      per-sample provenance that must NOT enter the model table

**Nexus compatibility is not claimed.** ``qecgen`` records the Nexus input format as
unknown and refuses to claim compatibility until an exporter passes a fixture supplied by
the Nexus team. Nothing here changes that; this format matches the shape the brief
describes, and that is a different statement.

**This is not a ``qecgen-ml-csv`` file and must not pretend to be one.** That format's
reader compares every feature cell literally against ``"0"`` and ``"1"``, and its sidecar
declares ``bit_values: ["0", "1"]``. A table of A=1..Z=26 tokens would be refused by it, so
this format carries its own marker and its own ``token_values``.
"""

from __future__ import annotations

import csv
import io
import json
import warnings
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import Any, Final

from nexus_genomics.common import (
    CONTENT_HASH_ALGORITHM,
    StagedWrite,
    blake2b_256,
    canonical_json,
    require_unique,
)
from nexus_genomics.encoding import MAX_TOKEN, PAD_TOKEN

__all__ = [
    "COLUMN_ORDER",
    "EXTENSION",
    "FORMAT",
    "INDEX_COLUMN",
    "TARGET_COLUMN",
    "VERSION",
    "NotANexusGenomicsDatasetError",
    "Row",
    "manifest_path",
    "pad_width",
    "position_columns",
    "read_nexus_csv",
    "samples_path",
    "target_columns",
    "write_nexus_csv",
]

FORMAT: Final = "nexus-genomics-csv"
VERSION: Final = 1
EXTENSION: Final = ".ml.csv"

INDEX_COLUMN: Final = "sample_id"
TARGET_COLUMN: Final = "target"
POSITION_PREFIX: Final = "position_"

COLUMN_ORDER: Final = ("index", "target", "feature")
"""The order the header builder, the row writer and the read offsets all derive from.

One tuple, three consumers -- the same discipline as ``qecgen.exporters.bit_columns``. When
these were three independent hardcodings in that project, changing two of the three stored
feature values under the target's column name, and the round trip stayed green because both
sides shared the header builder. The damage showed up only as a model that trained on
nothing.
"""

_MANIFEST_SUFFIX: Final = ".manifest.json"
_SAMPLES_SUFFIX: Final = ".samples.csv"
_LINE_TERMINATOR: Final = "\n"


SIZE_WARNING_CELLS: Final = 50_000_000
"""Cell count above which a one-position-per-column text table is an actively bad idea.

Borrowed from ``qecgen``'s SIZE_WARNING_THRESHOLD, expressed in cells rather than rows
because width matters as much as height here: 1,000 position columns turn a modest row count
into a very large file. The writer buffers the whole table to compute its content hash, so
peak memory runs several times the output size.
"""


def warn_if_large(n_rows: int, n_columns: int, name: str) -> None:
    """Warn before writing a text table that will be unreasonably large."""
    cells = n_rows * n_columns
    if cells <= SIZE_WARNING_CELLS:
        return
    warnings.warn(
        f"{name}: {n_rows:,} rows x {n_columns:,} columns is {cells:,} cells. This format "
        f"writes one column per position and buffers the table to hash it, so expect peak "
        f"memory of several times the output size. Reduce window.length or the record cap.",
        UserWarning,
        stacklevel=2,
    )


class NotANexusGenomicsDatasetError(Exception):
    """Raised when a file is provably not one of ours, as opposed to corrupt."""


def pad_width(count: int, *, one_based: bool) -> int:
    """Digits needed for the largest index in a block of ``count`` columns."""
    largest = count if one_based else max(count - 1, 0)
    return len(str(max(largest, 0)))


def position_columns(count: int) -> list[str]:
    """``position_0001 ... position_1000``. One-based, zero-padded to this file's own width.

    **One-based** because sequence positions are one-based everywhere in biology; a
    ``position_0`` would be an off-by-one waiting to happen when anyone compares this table
    against a GenBank coordinate. This deliberately differs from ``qecgen``'s zero-based
    ``detector_00``, where the index is an array offset rather than a biological coordinate.

    **Zero-padded** for the reason ``qecgen`` documents: ``sorted(df.columns)`` -- and every
    column sort hiding inside a join, a concat or a feature-store schema -- puts
    ``position_10`` before ``position_2``, which permutes the feature matrix silently and
    leaves a model that trains, converges and means nothing. The width comes from this
    file's own count, so it is a property of the file rather than a constant to hardcode,
    which is exactly why the sidecar publishes the literal names.
    """
    digits = pad_width(count, one_based=True)
    return [f"{POSITION_PREFIX}{i:0{digits}d}" for i in range(1, count + 1)]


def target_columns(count: int) -> list[str]:
    """``target`` for a single target, ``target_00 ... target_31`` for a multilabel one.

    The bare name for the single-target case follows ``qecgen.dataset.target_columns``: the
    column holding the answer says so in its own name, so a consumer never has to read the
    manifest to find it.

    The multi-target spelling is zero-padded, which ``qecgen`` does not do -- but that
    project only ever writes one observable, so its unpadded ``target_{i}`` was never
    exercised past nine. FunSoC has 32, where the unpadded spelling would sort
    ``target_10`` before ``target_2`` and permute the label matrix: the same defect the
    feature padding exists to prevent, and no less damaging for being on the target side.
    """
    if count < 1:
        raise ValueError("a Nexus table must have at least one target column")
    if count == 1:
        return [TARGET_COLUMN]
    digits = pad_width(count, one_based=False)
    return [f"{TARGET_COLUMN}_{i:0{digits}d}" for i in range(count)]


def manifest_path(path: Path) -> Path:
    """``d.ml.csv`` -> ``d.ml.manifest.json``. ``with_suffix`` replaces only the last suffix."""
    return path.with_suffix(_MANIFEST_SUFFIX)


def samples_path(path: Path) -> Path:
    """``d.ml.csv`` -> ``d.ml.samples.csv``."""
    return path.with_suffix(_SAMPLES_SUFFIX)


class Row:
    """One output row: primary key, targets, features."""

    __slots__ = ("features", "sample_id", "targets")

    def __init__(self, sample_id: str, targets: Sequence[int], features: Sequence[int]) -> None:
        self.sample_id = sample_id
        self.targets = tuple(targets)
        self.features = tuple(features)


def _columns_block(n_features: int, n_targets: int) -> dict[str, Any]:
    """The literal, ordered column names this file carries, split by role.

    Split rather than given as one flat list because the split *is* the contract: a consumer
    reads ``feature_columns`` and ``target_columns`` and never reconstructs either. That is
    what lets the pad width be a property of the file.
    """
    return {
        "index_columns": [INDEX_COLUMN],
        "target_columns": target_columns(n_targets),
        "feature_columns": position_columns(n_features),
        "token_values": {
            "padding": PAD_TOKEN,
            "letters": f"A=1 .. Z={MAX_TOKEN}",
            "note": (
                "Categorical tokens. Numerical distance between them has no biological "
                "meaning; treat these columns as categories, not magnitudes."
            ),
        },
        "row_order": (
            f"{INDEX_COLUMN} is unique; rows are in emission order. Sorting is safe because "
            f"every row is self-describing, unlike qecgen's positional 'shot' column."
        ),
    }


def _header(columns: Mapping[str, Any]) -> list[str]:
    blocks = {
        "index": list(columns["index_columns"]),
        "target": list(columns["target_columns"]),
        "feature": list(columns["feature_columns"]),
    }
    return [name for block in COLUMN_ORDER for name in blocks[block]]


def write_nexus_csv(
    path: Path,
    rows: Sequence[Row],
    *,
    n_features: int,
    n_targets: int,
    manifest: Mapping[str, Any],
    sample_metadata: Sequence[Mapping[str, Any]] | None = None,
) -> list[Path]:
    """Write the table and its sidecars as one atomic group."""
    if not rows:
        raise ValueError(
            f"{path.name}: refusing to write a table with no rows. An empty dataset is "
            f"almost always a silently failed adapter, and it would pass most validation."
        )
    require_unique([r.sample_id for r in rows], path.name)

    columns = _columns_block(n_features, n_targets)
    header = _header(columns)
    warn_if_large(len(rows), len(header), path.name)

    for row in rows:
        if len(row.features) != n_features:
            raise ValueError(
                f"{path.name}: row {row.sample_id!r} has {len(row.features)} features but "
                f"the header declares {n_features}; padding it here would fabricate residues"
            )
        if len(row.targets) != n_targets:
            raise ValueError(
                f"{path.name}: row {row.sample_id!r} has {len(row.targets)} targets but the "
                f"header declares {n_targets}"
            )

    buffer = io.StringIO(newline="")
    writer = csv.writer(buffer, lineterminator=_LINE_TERMINATOR)
    writer.writerow(header)
    for row in rows:
        cells = {
            "index": [row.sample_id],
            "target": [str(t) for t in row.targets],
            "feature": [str(t) for t in row.features],
        }
        writer.writerow([cell for block in COLUMN_ORDER for cell in cells[block]])
    table = buffer.getvalue()

    # The digest covers the table only. `generated_at` lives in the manifest and changes on
    # every run, so folding the manifest in would make the reproducibility check vacuous.
    content_hash = blake2b_256(table.encode("utf-8"))

    sidecar = {
        "format": FORMAT,
        "version": VERSION,
        "columns": columns,
        "manifest": {
            **dict(manifest),
            "n_samples": len(rows),
            "n_features": n_features,
            "n_targets": n_targets,
            "content_hash": content_hash,
            "content_hash_algorithm": CONTENT_HASH_ALGORITHM,
        },
    }

    staged = StagedWrite()
    staged.add_bytes(path, table.encode("utf-8"))
    staged.add_text(manifest_path(path), canonical_json(sidecar))
    if sample_metadata:
        staged.add_text(samples_path(path), _samples_csv(sample_metadata))
    return staged.commit()


def _samples_csv(records: Sequence[Mapping[str, Any]]) -> str:
    """Per-sample provenance, kept out of the model table.

    Richer metadata in the Nexus table itself risks an ingestion failure over a column the
    consumer did not expect, and a free-text field is a CSV-quoting hazard. It lives here
    instead, keyed by the same ``sample_id``.
    """
    fields: list[str] = []
    for record in records:
        for key in record:
            if key not in fields:
                fields.append(key)
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields, lineterminator=_LINE_TERMINATOR)
    writer.writeheader()
    for record in records:
        writer.writerow({k: record.get(k, "") for k in fields})
    return buffer.getvalue()


def _load_sidecar(path: Path) -> dict[str, Any]:
    """The manifest sidecar, refusing a foreign or future one. The format's magic line."""
    companion = manifest_path(path)
    if not companion.exists():
        raise NotANexusGenomicsDatasetError(
            f"{path.name} has no {companion.name} beside it, so it is a plain CSV rather "
            f"than a nexus-genomics dataset."
        )
    payload: dict[str, Any] = json.loads(companion.read_text(encoding="utf-8"))
    if payload.get("format") != FORMAT:
        raise NotANexusGenomicsDatasetError(
            f"{companion.name} is not a {FORMAT} sidecar (format={payload.get('format')!r})."
        )
    if payload.get("version") != VERSION:
        raise ValueError(
            f"{companion.name} declares version {payload.get('version')!r}, but this reader "
            f"understands only {VERSION}. Refused rather than read as this version, because "
            f"a later revision may add a required key."
        )
    return payload


def read_nexus_csv(path: Path) -> tuple[dict[str, Any], list[str], list[list[str]]]:
    """Read a table back, refusing anything that is provably not ours.

    Returns the sidecar payload, the header and the raw rows as strings. Deliberately does
    not coerce to int here: :mod:`nexus_genomics.validation` needs to see the literal cells
    to tell a genuine ``0`` apart from an empty field a spreadsheet left behind.
    """
    companion = manifest_path(path)
    if not companion.exists():
        raise NotANexusGenomicsDatasetError(
            f"{path.name} has no {companion.name} beside it, so it is a plain CSV rather "
            f"than a nexus-genomics dataset. This format keeps its manifest in that "
            f"sidecar and commits the whole set in one move, so a table standing alone was "
            f"not written by this tool."
        )
    payload = _load_sidecar(path)

    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration:
            raise NotANexusGenomicsDatasetError(f"{path.name} is empty.") from None
        rows = [row for row in reader if row]

    expected = _header(payload["columns"])
    if header != expected:
        raise ValueError(
            f"{path}: the header row disagrees with its sidecar. The file has {len(header)} "
            f"columns and the sidecar names {len(expected)}. Resolving that in favour of "
            f"either one would be a guess about which half is corrupt."
        )
    return payload, header, rows


def read_sidecar_and_shape(path: Path) -> tuple[dict[str, Any], list[str], int, set[int]]:
    """The sidecar, the header, the row count and the set of row widths -- without the rows.

    :func:`read_nexus_csv` materialises every row as a list of strings, which for the 70 MB
    CodonTransformer table is a second whole copy on top of the pandas frame the validator
    already holds. The validator only ever needs the count and the widths, so it uses this
    streaming form instead and the memory stays bounded by one row.
    """
    payload = _load_sidecar(path)
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration:
            raise NotANexusGenomicsDatasetError(f"{path.name} is empty.") from None
        widths: set[int] = set()
        count = 0
        for row in reader:
            if not row:
                continue
            count += 1
            widths.add(len(row))

    expected = _header(payload["columns"])
    if header != expected:
        raise ValueError(
            f"{path}: the header row disagrees with its sidecar. The file has {len(header)} "
            f"columns and the sidecar names {len(expected)}."
        )
    return payload, header, count, widths


def iter_rows_as_ints(header: Sequence[str], rows: Iterable[Sequence[str]]) -> Iterable[list[int]]:
    """Parse feature cells strictly, refusing anything a spreadsheet may have rewritten."""
    del header
    for row in rows:
        yield [int(cell) for cell in row]
