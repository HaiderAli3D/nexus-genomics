"""Project a labelled source into deterministic single-target sources."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from nexus_genomics.adapters.base import LoadedSource
from nexus_genomics.common import SampleRecord

__all__ = ["TargetProjection", "projected_path", "single_target_projections"]

_OUTPUT_SUFFIX = ".ml.csv"


@dataclass(frozen=True, slots=True)
class TargetProjection:
    """One target's index, provenance name, filename slug, and projected source."""

    index: int
    name: str
    slug: str
    source: LoadedSource


def single_target_projections(source: LoadedSource) -> list[TargetProjection]:
    """Return one validated source per target without changing the loaded source."""
    if not source.records:
        raise ValueError("target projection source has no records")

    validated: list[tuple[SampleRecord, tuple[int, ...]]] = []
    for record in source.records:
        if record.label is None:
            raise ValueError(f"record {record.sample_id!r} has a missing label")
        if len(record.label) != source.n_targets:
            raise ValueError(
                f"record {record.sample_id!r} label must have exactly "
                f"{source.n_targets} values, got {len(record.label)}"
            )
        validated.append((record, record.label))

    if source.n_targets == 1:
        return [TargetProjection(index=0, name="", slug="", source=source)]

    if len(source.target_names) != source.n_targets:
        raise ValueError(
            "multilabel projection requires exactly one target name per target: "
            f"expected {source.n_targets}, got {len(source.target_names)}"
        )

    projections: list[TargetProjection] = []
    slugs: set[str] = set()
    for index, name in enumerate(source.target_names):
        slug = re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")
        if not slug:
            raise ValueError(f"target name {name!r} produces an empty slug")
        if slug in slugs:
            raise ValueError(f"target name {name!r} produces duplicate slug {slug!r}")
        slugs.add(slug)

        manifest_extra = dict(source.manifest_extra)
        manifest_extra["target_projection"] = {
            "original_target_index": index,
            "original_target_name": name,
            "original_target_count": source.n_targets,
        }
        records = [
            SampleRecord(
                sample_id=record.sample_id,
                sequence=record.sequence,
                label=(label[index],),
                source=record.source,
                metadata=record.metadata,
            )
            for record, label in validated
        ]
        projected_source = LoadedSource(
            records=records,
            n_targets=1,
            label_semantics=source.label_semantics,
            alphabet=source.alphabet,
            target_names=(),
            rna_to_dna=source.rna_to_dna,
            manifest_extra=manifest_extra,
            raw_input_hashes=source.raw_input_hashes,
            sample_metadata_fields=source.sample_metadata_fields,
        )
        projections.append(TargetProjection(index, name, slug, projected_source))
    return projections


def projected_path(base_path: Path, slug: str) -> Path:
    """Insert a non-empty target slug before the recognized table extension."""
    if not slug:
        return base_path
    if not base_path.name.endswith(_OUTPUT_SUFFIX):
        raise ValueError(f"projected output path must end with {_OUTPUT_SUFFIX!r}: {base_path}")
    stem = base_path.name[: -len(_OUTPUT_SUFFIX)]
    return base_path.with_name(f"{stem}__{slug}{_OUTPUT_SUFFIX}")
