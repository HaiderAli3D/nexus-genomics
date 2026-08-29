"""Addgene / Genetic Engineering Attribution Challenge.

**Every sequence in this corpus is an engineered Addgene plasmid.** There is no natural
class, there never was one, and this adapter will not manufacture one. The real ground truth
is **lab of origin** -- a 1,314-way multiclass target -- and that is what it emits.

One trap deserves naming. The label set contains a pooled class, ``I7FXTVDP``, holding the
8,286 training sequences (13.15%) from the 2,438 labs that deposited fewer than ten
plasmids. The competition calls it **"Unknown Engineered"**. That means *unknown laboratory*,
not *unknown provenance*: everything in it is still an engineered plasmid. Reading it as an
"unknown or natural" class would invert its meaning, so
:func:`~nexus_genomics.adapters.addgene.AddgeneAdapter.load` records it explicitly in
``label_semantics`` and a test asserts that no class in the emitted file is called natural.

``train_labels.csv`` is **one-hot**: 63,017 rows by 1,315 columns, one column per lab ID,
exactly one ``1.0`` per row. It is collapsed with ``idxmax(axis=1)``, which is what
DrivenData's own benchmark does. Read in chunks, because materialising 63,017 x 1,314
float64 cells at once is roughly 660 MB for a column index.

**Access and redistribution.** The files sit behind a DrivenData account *and* acceptance of
the competition rules, which state that participants "agree not to transmit, duplicate,
publish, redistribute or otherwise provide or make available the Data to any party not
participating in the Competition". So the raw files cannot be committed here and neither can
a derived CSV containing the sequences: outputs from this adapter stay local. An unofficial
Kaggle mirror exists; it confers no rights and this adapter deliberately does not use it.
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Any, Final

import pandas as pd

from nexus_genomics.adapters.base import Availability, LoadedSource, stable_class_index
from nexus_genomics.cleaning import Alphabet
from nexus_genomics.common import SampleRecord, hash_file

__all__ = ["POOLED_UNKNOWN_LAB", "AddgeneAdapter"]

VALUES_NAME: Final = "train_values.csv"
LABELS_NAME: Final = "train_labels.csv"
COMPETITION_URL: Final = (
    "https://www.drivendata.org/competitions/63/genetic-engineering-attribution/"
)

_PHENOTYPE_PREFIXES: Final = (
    "bacterial_resistance_",
    "copy_number_",
    "growth_",
    "selectable_markers_",
    "species_",
)
"""The 39 binary metadata columns, by prefix. Kept out of the model table; see load()."""

POOLED_UNKNOWN_LAB: Final = "I7FXTVDP"
"""The pooled 'Unknown Engineered' class. Unknown LAB, not unknown provenance."""

_STEPS: Final = f"""\
The competition data is gated. Exact steps:

  1. Create a free account at https://www.drivendata.org/accounts/signup/ and verify it.
  2. Log in, open {COMPETITION_URL} and ACCEPT THE COMPETITION RULES -- accepting is the
     act that grants data access.
  3. Open the 'Data' tab and download {VALUES_NAME} and {LABELS_NAME}.
     (The links are per-user pre-signed S3 URLs that expire after about 24 hours, so if you
     script the download, quote the URL and do it promptly.)
  4. Place both files in data/raw/addgene/

The rules forbid redistributing the data or anything derived from it, so the CSV this
adapter produces must stay on your machine. Do not commit it.\
"""


class AddgeneAdapter:
    """Lab-of-origin attribution. Multiclass, and emphatically not a binary."""

    name = "addgene_lab_attribution"
    source_url = COMPETITION_URL
    licence = (
        "DrivenData competition rules; redistribution of the Data or derived data is "
        "PROHIBITED. Outputs must stay local."
    )
    target_description = (
        "lab_id: which laboratory deposited this plasmid (multiclass, up to 1,314 classes). "
        "There is NO natural class -- 100% of this corpus is engineered plasmids."
    )

    def availability(self, raw_dir: Path) -> Availability:
        values, labels = raw_dir / VALUES_NAME, raw_dir / LABELS_NAME
        missing = [p.name for p in (values, labels) if not p.exists()]
        if not missing:
            return Availability(
                True,
                f"{VALUES_NAME} ({values.stat().st_size:,} B) and "
                f"{LABELS_NAME} ({labels.stat().st_size:,} B) present",
            )
        return Availability(
            False,
            f"missing {', '.join(missing)} in {raw_dir}",
            manual_steps=_STEPS,
            expected_files=(VALUES_NAME, LABELS_NAME),
        )

    def load(self, raw_dir: Path, options: Mapping[str, Any]) -> LoadedSource:
        availability = self.availability(raw_dir)
        if not availability.ready:
            raise FileNotFoundError(f"{availability.detail}\n\n{availability.manual_steps}")

        values_path, labels_path = raw_dir / VALUES_NAME, raw_dir / LABELS_NAME
        limit = options.get("max_records")
        n = int(limit) if limit else None

        lab_of = _collapse_one_hot(labels_path, n)
        frame = pd.read_csv(
            values_path,
            usecols=lambda c: c in {"sequence_id", "sequence"} or c.startswith(_PHENOTYPE_PREFIXES),
            nrows=n,
        )

        ids = [str(s) for s in frame["sequence_id"]]
        duplicated = len(ids) - len(set(ids))
        if duplicated:
            raise ValueError(
                f"{VALUES_NAME} has {duplicated} duplicated sequence_id value(s). Every "
                f"lookup below assumes the id is unique; with a duplicate, pandas returns a "
                f"DataFrame rather than a row and the sequence would be silently replaced by "
                f"a repr of a Series. Refused rather than corrupted."
            )
        aligned = [sid for sid in ids if sid in lab_of]

        # The class vocabulary is the HEADER of train_labels.csv -- every lab column the
        # competition defines -- not the labs that happen to survive `max_records`. Deriving
        # it from the sample would renumber classes whenever the cap changes, so a model
        # trained on the demo file and scored against the full one would be graded against
        # the wrong labs. Reading one header line is cheap and is the authoritative list.
        classes = stable_class_index(_label_vocabulary(labels_path))
        unknown = sorted({lab_of[sid] for sid in aligned} - set(classes))
        if unknown:
            raise ValueError(
                f"{LABELS_NAME}: {len(unknown)} lab id(s) appear in the collapsed labels but "
                f"not in the header (e.g. {unknown[:3]}). The header is the class vocabulary; "
                f"a mismatch means the two files disagree."
            )

        metadata_columns = [c for c in frame.columns if c not in {"sequence_id", "sequence"}]
        # Columns pulled out once and indexed by row position, rather than `.loc` per row.
        # `.loc` inside a 63,017-iteration loop constructs a Series each time, which turns a
        # few seconds into minutes on the full file -- the archetypal defect that never shows
        # up on a 500-row demo.
        column_values = {c: frame[c].tolist() for c in frame.columns}
        position = {sid: i for i, sid in enumerate(ids)}

        records: list[SampleRecord] = []
        for sid in aligned:
            at = position[sid]
            lab = lab_of[sid]
            records.append(
                SampleRecord(
                    sample_id=sid,
                    sequence=str(column_values["sequence"][at]),
                    label=(classes[lab],),
                    source=self.name,
                    metadata={
                        "lab_id": lab,
                        "is_pooled_unknown_lab": int(lab == POOLED_UNKNOWN_LAB),
                        **{c: column_values[c][at] for c in metadata_columns},
                    },
                )
            )

        semantics: dict[int, str] = {}
        for lab, index in classes.items():
            if lab == POOLED_UNKNOWN_LAB:
                semantics[index] = (
                    f"lab_id {lab} -- the pooled 'Unknown Engineered' class: labs that "
                    f"deposited fewer than 10 plasmids. UNKNOWN LABORATORY, still an "
                    f"engineered plasmid. This is NOT a natural class."
                )
            else:
                semantics[index] = f"lab_id {lab} (obfuscated 8-character identifier)"

        return LoadedSource(
            records=records,
            n_targets=1,
            label_semantics=semantics,
            alphabet=Alphabet.DNA_NUCLEOTIDE,
            manifest_extra={
                "no_natural_class": (
                    "Every sequence here is an engineered plasmid deposited at Addgene. A "
                    "natural-versus-engineered target does not exist for this source and "
                    "was not fabricated."
                ),
                "pooled_unknown_lab_id": POOLED_UNKNOWN_LAB,
                "pooled_unknown_meaning": (
                    "Unknown LABORATORY, not unknown provenance. 8,286 of 63,017 training "
                    "rows (13.15%), pooled from 2,438 labs with fewer than 10 deposits."
                ),
                "class_map": {str(v): k for k, v in classes.items()},
                "n_classes": len(classes),
                "phenotype_columns_location": (
                    "The 39 binary phenotype/metadata columns are real numeric variables but "
                    "are kept in the .samples.csv sidecar rather than mixed into the "
                    "position_* block, so the feature block stays one-position-per-column. "
                    "Join them on sample_id if you want them."
                ),
                "redistribution": "PROHIBITED by the competition rules. Keep outputs local.",
                "known_duplicates_upstream": (
                    "The official train_values.csv contains 374 duplicated sequence strings "
                    "over 63,017 rows; the validation report will surface them."
                ),
            },
            raw_input_hashes={
                VALUES_NAME: hash_file(values_path),
                LABELS_NAME: hash_file(labels_path),
            },
            sample_metadata_fields=("lab_id", "is_pooled_unknown_lab", *metadata_columns),
        )


def _label_vocabulary(path: Path) -> list[str]:
    """Every lab id the competition defines, read from the header of ``train_labels.csv``.

    One line of I/O against a 331 MB file, and it is the only list that does not shift when
    the sampling cap changes.
    """
    header = pd.read_csv(path, nrows=0).columns.tolist()
    return [str(c) for c in header[1:]]


def _collapse_one_hot(path: Path, limit: int | None) -> dict[str, str]:
    """One-hot 63,017 x 1,315 -> {sequence_id: lab_id}, in chunks.

    Chunked because the whole frame as float64 is roughly 660 MB, and the demo path should
    not need a machine that can hold it. ``idxmax`` per chunk gives the same answer as
    ``idxmax`` over the whole frame, since exactly one column is set per row.
    """
    out: dict[str, str] = {}
    for chunk in pd.read_csv(path, index_col=0, chunksize=2000):
        # `idxmax` on an all-zero row returns the FIRST column and says nothing about it.
        # That would silently attribute a plasmid to whichever lab happens to sort first --
        # a fabricated label, produced by a function that looks like it is reading one.
        # Every row must have exactly one positive, and a row that does not is refused.
        row_max = chunk.max(axis=1)
        row_sum = chunk.sum(axis=1)
        bad = chunk.index[(row_max != 1) | (row_sum != 1)]
        if len(bad):
            raise ValueError(
                f"{path.name}: {len(bad)} row(s) are not one-hot (e.g. {list(bad[:3])}). "
                f"train_labels.csv must have exactly one 1.0 per row; idxmax on an all-zero "
                f"row would invent a lab, so this is refused rather than collapsed."
            )
        for sequence_id, lab in chunk.idxmax(axis=1).items():
            out[str(sequence_id)] = str(lab)
        if limit and len(out) >= limit:
            break
    return out
