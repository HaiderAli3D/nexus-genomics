"""What every adapter must provide, and nothing more.

An adapter's whole job is to turn one source's native format into
:class:`~nexus_genomics.common.SampleRecord` objects and to state, honestly, what its labels
mean. It does no cleaning, no encoding and no padding -- those are shared and must behave
identically for all four sources, or a defect fixed in one place survives in three others.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Protocol, runtime_checkable

from nexus_genomics.cleaning import Alphabet
from nexus_genomics.common import SampleRecord

__all__ = ["Adapter", "Availability", "LoadedSource"]


@dataclass(frozen=True, slots=True)
class Availability:
    """Whether this adapter can run right now, and if not, exactly what a human must do.

    ``manual_steps`` is not optional politeness. Two of the four sources are gated, and the
    difference between "this pipeline is broken" and "this pipeline needs one file you have
    to fetch in a browser" is the difference between a wasted afternoon and a five-minute
    task.
    """

    ready: bool
    detail: str
    manual_steps: str = ""
    expected_files: tuple[str, ...] = ()


@dataclass(slots=True)
class LoadedSource:
    """Everything the shared pipeline needs, once an adapter has done its part."""

    records: list[SampleRecord]
    n_targets: int
    label_semantics: Mapping[int, str]
    """Every value that may legally appear in a target column, and what it means.

    Written into the sidecar and checked by validation. A source whose labels are not
    ground truth says so here, in the text a reader will actually see.
    """
    alphabet: Alphabet
    target_names: tuple[str, ...] = ()
    rna_to_dna: bool = False
    manifest_extra: dict[str, Any] = field(default_factory=dict)
    raw_input_hashes: dict[str, str] = field(default_factory=dict)
    sample_metadata_fields: tuple[str, ...] = ()
    """Keys of ``SampleRecord.metadata`` to copy into the ``.samples.csv`` sidecar."""


@runtime_checkable
class Adapter(Protocol):
    """The one interface. Adding a source is one module plus one line in ``ADAPTERS``."""

    @property
    def name(self) -> str: ...

    @property
    def source_url(self) -> str: ...

    @property
    def licence(self) -> str: ...

    @property
    def target_description(self) -> str:
        """What the target column actually is. Never a claim the source does not support."""
        ...

    def availability(self, raw_dir: Path) -> Availability: ...

    def load(self, raw_dir: Path, options: Mapping[str, Any]) -> LoadedSource: ...


def one_hot(index: int, width: int) -> tuple[int, ...]:
    """A multilabel target vector with a single positive, used by the multiclass adapters."""
    if not 0 <= index < width:
        raise ValueError(f"class index {index} outside 0..{width - 1}")
    return tuple(1 if i == index else 0 for i in range(width))


def stable_class_index(values: Sequence[str]) -> dict[str, int]:
    """Sorted, deterministic string -> integer class map.

    Sorted rather than first-seen so that two runs over the same data -- in any row order --
    assign the same integer to the same class. A first-seen map would make ``content_hash``
    depend on file ordering, and the reproducibility check would fail for a reason that has
    nothing to do with the data.
    """
    return {value: i for i, value in enumerate(sorted(set(values)))}
