"""Records in, Nexus CSV out. The one path every source travels.

Cleaning, encoding, padding, windowing and writing happen here and only here, so that all
five sources are subject to identical rules and a fix lands in one place.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from nexus_genomics.adapters.base import LoadedSource
from nexus_genomics.cleaning import CleaningReport, Quarantined, clean_sequence
from nexus_genomics.common import (
    PACKAGE_VERSION,
    SampleRecord,
    blake2b_256,
    utc_now_iso,
    versions_block,
)
from nexus_genomics.encoding import EncodingMode, LengthPolicy, encode, place
from nexus_genomics.nexus_csv import Row, write_nexus_csv

__all__ = ["ConvertOptions", "ConvertOutcome", "convert"]


@dataclass(frozen=True, slots=True)
class ConvertOptions:
    width: int = 1000
    policy: LengthPolicy = LengthPolicy.PAD_OR_TRUNCATE
    stride: int = 1000
    mode: EncodingMode = EncodingMode.LETTER
    max_windows_per_sample: int | None = 8
    """Cap on windows from one sequence, so one 52 kb feature cannot dominate a demo file.

    Applied deterministically (the first N windows) and COUNTED, never silently. A capped
    sample records ``windows_dropped`` so the tail it lost is visible in the sidecar.
    """

    def __post_init__(self) -> None:
        cap = self.max_windows_per_sample
        if cap is None:
            return
        if type(cap) is not int:
            raise TypeError(
                "max_windows_per_sample must be None or a non-boolean integer greater "
                f"than zero, got {cap!r} ({type(cap).__name__})"
            )
        if cap <= 0:
            raise ValueError(
                f"max_windows_per_sample must be greater than zero when configured, got {cap}"
            )


@dataclass(slots=True)
class ConvertOutcome:
    path: Path
    written: list[Path]
    n_rows: int
    n_features: int
    cleaning: CleaningReport
    quarantine: list[tuple[str, str, str]]
    windows_dropped: int
    residues_dropped_by_cap: int
    truncated_residues: int


def _require_supported(mode: EncodingMode) -> None:
    if mode is EncodingMode.LETTER:
        return
    raise NotImplementedError(
        f"encoding mode {mode.value!r} is defined but not implemented. "
        + (
            "A codon representation needs a reading frame, which a plasmid or a genome does "
            "not have a single one of."
            if mode is EncodingMode.CODON
            else "Translating nucleotides to residues requires knowing where the coding "
            "regions are and in which frame; translating a whole plasmid end to end would "
            "emit a protein that does not exist. A source whose sequences are already amino "
            "acids should set alphabet=amino_acid and use the 'letter' mode."
        )
    )


def convert(
    source: LoadedSource,
    out_path: Path,
    options: ConvertOptions,
    *,
    source_name: str,
    source_url: str,
    licence: str,
    target_description: str,
    notes: str,
    config_digest: str,
    extra_manifest: Mapping[str, Any] | None = None,
) -> ConvertOutcome:
    """Clean, encode, place and write. Every discarded residue is counted."""
    _require_supported(options.mode)

    report = CleaningReport()
    quarantine: list[tuple[str, str, str]] = []
    rows: list[Row] = []
    sample_meta: list[dict[str, Any]] = []
    windows_dropped = 0
    residues_dropped_by_cap = 0
    truncated = 0

    for record in source.records:
        cleaned = clean_sequence(
            record.sequence,
            source.alphabet,
            report=report,
            rna_to_dna=source.rna_to_dna,
        )
        if isinstance(cleaned, Quarantined):
            quarantine.append((record.sample_id, cleaned.reason, cleaned.detail))
            continue
        if record.label is None:
            raise ValueError(
                f"{record.sample_id}: adapter produced a record with no label. A missing "
                f"label must never default to the negative class; fix the adapter or drop "
                f"the record explicitly."
            )

        tokens = encode(cleaned.sequence)
        placements = place(tokens, options.width, options.policy, options.stride)

        cap = options.max_windows_per_sample
        dropped_by_cap = 0
        if options.policy is LengthPolicy.WINDOW and cap is not None and len(placements) > cap:
            windows_dropped += len(placements) - cap
            # The residues in the dropped windows have to be counted here. `place` reports
            # removed=0 for every window because windowing itself discards nothing; the cap
            # is what discards, and if it does not say so the sidecar reports a sequence
            # that lost 70% of its length as having lost none.
            kept_to = placements[cap - 1].end
            dropped_by_cap = max(cleaned.original_length - kept_to, 0)
            placements = placements[:cap]
        residues_dropped_by_cap += dropped_by_cap

        multi = len(placements) > 1 or options.policy is LengthPolicy.WINDOW
        digits = len(str(max(len(placements) - 1, 0)))
        for i, placed in enumerate(placements):
            sample_id = (
                f"{record.sample_id}__w{i:0{max(digits, 2)}d}" if multi else record.sample_id
            )
            rows.append(Row(sample_id, record.label, placed.tokens))
            truncated += placed.removed
            # Hashed from the CLEANED SOURCE, never from `decode(placed.tokens)`.
            #
            # This is the difference between a round-trip check and a tautology. Hashing the
            # decoded tokens would compare the encoder's output against itself: introduce an
            # off-by-one in `place` -- `tokens[1 : width + 1]` instead of `tokens[:width]` --
            # and every row would silently lose its first residue while the digest matched
            # perfectly and the validator reported a pass on every file.
            source_slice = cleaned.sequence[placed.start - 1 : placed.end]
            sample_meta.append(
                {
                    "sample_id": sample_id,
                    "parent_sample_id": record.sample_id,
                    "source": record.source,
                    "original_length": cleaned.original_length,
                    "window_start": placed.start,
                    "window_end": placed.end,
                    "residues_in_row": len(source_slice),
                    "padding_cells": placed.pad_count,
                    "residues_removed_by_truncation": placed.removed,
                    "residues_dropped_by_window_cap": (
                        dropped_by_cap if i == len(placements) - 1 else 0
                    ),
                    "source_slice_blake2b": blake2b_256(source_slice.encode("utf-8")),
                    **{k: record.metadata.get(k, "") for k in source.sample_metadata_fields},
                }
            )

    manifest: dict[str, Any] = {
        "source": source_name,
        "source_url": source_url,
        "licence": licence,
        "target_description": target_description,
        "label_semantics": {str(k): v for k, v in source.label_semantics.items()},
        "alphabet": source.alphabet.value,
        "encoding_mode": options.mode.value,
        "encoding": {
            "cipher": "A=1 .. Z=26",
            "padding_token": 0,
            "reversible": True,
            "warning": (
                "These integers are categorical tokens. Their numerical distance has no "
                "biological meaning; a model must treat them as categories, not magnitudes."
            ),
        },
        "window": {
            "policy": options.policy.value,
            "length": options.width,
            "stride": options.stride if options.policy is LengthPolicy.WINDOW else None,
            "max_windows_per_sample": options.max_windows_per_sample,
            "windows_dropped_by_cap": windows_dropped,
            "residues_dropped_by_cap": residues_dropped_by_cap,
            "residues_removed_by_truncation": truncated,
        },
        "cleaning": report.to_json_dict(),
        "quarantined_samples": [
            {"sample_id": s, "reason": r, "detail": d} for s, r, d in quarantine
        ],
        "raw_input_hashes": source.raw_input_hashes,
        "config_hash": config_digest,
        "generated_at": utc_now_iso(),
        "generator": f"nexus-genomics {PACKAGE_VERSION}",
        "versions": versions_block(),
        "nexus_compatibility": (
            "NOT VERIFIED. This file matches the CSV shape described in the brief and "
            "mirrors the conventions of qecgen's ml_csv exporter. No fixture supplied by "
            "the Nexus team has been run against it."
        ),
        "notes": notes,
        **(dict(extra_manifest) if extra_manifest else {}),
    }

    written = write_nexus_csv(
        out_path,
        rows,
        n_features=options.width,
        n_targets=source.n_targets,
        manifest=manifest,
        sample_metadata=sample_meta,
    )
    return ConvertOutcome(
        path=out_path,
        written=written,
        n_rows=len(rows),
        n_features=options.width,
        cleaning=report,
        quarantine=quarantine,
        windows_dropped=windows_dropped,
        residues_dropped_by_cap=residues_dropped_by_cap,
        truncated_residues=truncated,
    )


def records_summary(records: list[SampleRecord]) -> dict[str, Any]:
    """Cheap shape report used by the CLI before any conversion happens."""
    lengths = [len(r.sequence) for r in records]
    return {
        "n_records": len(records),
        "min_length": min(lengths) if lengths else 0,
        "max_length": max(lengths) if lengths else 0,
        "median_length": sorted(lengths)[len(lengths) // 2] if lengths else 0,
    }
