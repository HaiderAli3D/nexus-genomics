"""Project a labelled source into deterministic single-target sources."""

from __future__ import annotations

import json
import re
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from nexus_genomics.adapters.base import LoadedSource
from nexus_genomics.common import CONTENT_HASH_ALGORITHM, SampleRecord, blake2b_256

__all__ = [
    "TargetProjection",
    "ordered_target_hash",
    "projected_path",
    "single_target_projections",
]

_OUTPUT_SUFFIX = ".ml.csv"


@dataclass(frozen=True, slots=True)
class TargetProjection:
    """One target's index, provenance name, filename slug, and projected source."""

    index: int
    name: str
    slug: str
    source: LoadedSource


def ordered_target_hash(sample_ids: Sequence[str], target_values: Sequence[int]) -> str:
    """Bind one target to its ordered primary keys without accepting typed-label drift."""
    if len(sample_ids) != len(target_values):
        raise ValueError(
            f"target hash needs one value per sample id, got {len(sample_ids)} and "
            f"{len(target_values)}"
        )
    pairs: list[list[str | int]] = []
    for sample_id, value in zip(sample_ids, target_values, strict=True):
        if type(value) is not int:
            raise TypeError(
                f"target hash values must be exact integers, got {value!r} "
                f"({type(value).__name__}) for {sample_id!r}"
            )
        pairs.append([sample_id, value])
    encoded = json.dumps(
        pairs,
        ensure_ascii=False,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return blake2b_256(encoded)


def _slug(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")


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
        invalid = [value for value in record.label if type(value) is not int]
        if invalid:
            value = invalid[0]
            raise TypeError(
                f"record {record.sample_id!r} target values must be exact integers; got "
                f"{value!r} ({type(value).__name__})"
            )
        validated.append((record, record.label))

    if source.n_targets == 1:
        return [TargetProjection(index=0, name="", slug="", source=source)]

    if len(source.target_names) != source.n_targets:
        raise ValueError(
            "multilabel projection requires exactly one target name per target: "
            f"expected {source.n_targets}, got {len(source.target_names)}"
        )

    target_slugs: list[str] = []
    seen_slugs: set[str] = set()
    for name in source.target_names:
        slug = _slug(name)
        if not slug:
            raise ValueError(f"target name {name!r} produces an empty slug")
        if slug in seen_slugs:
            raise ValueError(f"target name {name!r} produces duplicate slug {slug!r}")
        seen_slugs.add(slug)
        target_slugs.append(slug)

    sample_ids = [record.sample_id for record, _label in validated]
    target_hashes = [
        ordered_target_hash(sample_ids, [label[index] for _record, label in validated])
        for index in range(source.n_targets)
    ]
    source_contract = {
        "source_target_names": list(source.target_names),
        "source_target_slugs": target_slugs,
        "source_target_hashes": target_hashes,
        "source_target_hash_algorithm": CONTENT_HASH_ALGORITHM,
        "source_row_count": len(validated),
    }

    projections: list[TargetProjection] = []
    for index, (name, slug) in enumerate(zip(source.target_names, target_slugs, strict=True)):
        manifest_extra = dict(source.manifest_extra)
        manifest_extra["target_projection"] = {
            "original_target_index": index,
            "original_target_name": name,
            "original_target_slug": slug,
            "original_target_count": source.n_targets,
            **source_contract,
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
