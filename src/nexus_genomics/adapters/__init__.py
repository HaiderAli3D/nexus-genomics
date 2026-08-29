"""Source-specific readers. Adding a source is one module plus one line in ``ADAPTERS``."""

from __future__ import annotations

from typing import Final

from nexus_genomics.adapters.addgene import AddgeneAdapter
from nexus_genomics.adapters.base import Adapter, Availability, LoadedSource
from nexus_genomics.adapters.codontransformer import CodonTransformerAdapter
from nexus_genomics.adapters.felix_guardian import FelixGuardianAdapter
from nexus_genomics.adapters.seqscreen import SeqScreenFunSoCAdapter

__all__ = [
    "ADAPTERS",
    "Adapter",
    "Availability",
    "LoadedSource",
    "get_adapter",
]

ADAPTERS: Final[dict[str, Adapter]] = {
    "codontransformer": CodonTransformerAdapter(),
    "seqscreen": SeqScreenFunSoCAdapter(),
    "felix": FelixGuardianAdapter(),
    "addgene": AddgeneAdapter(),
}
"""Keyed by the short name the CLI takes, ordered by how ready each one is to run.

``codontransformer`` and ``seqscreen`` need no credentials. ``felix`` downloads an
Apache-2.0 repository but cannot produce a natural-versus-engineered file from it (see that
module). ``addgene`` needs a gated manual download.
"""


def get_adapter(name: str) -> Adapter:
    try:
        return ADAPTERS[name]
    except KeyError:
        raise KeyError(
            f"unknown source {name!r}. Known sources: {', '.join(sorted(ADAPTERS))}"
        ) from None
