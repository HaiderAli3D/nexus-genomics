"""The internal representation every adapter produces, and the primitives around it.

One dataclass and a handful of functions, so that the four adapters cannot each invent
their own idea of what a sample is. Everything downstream of :class:`SampleRecord` is
source-agnostic; everything upstream of it is source-specific and lives in ``adapters/``.
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Final

__all__ = [
    "CONTENT_HASH_ALGORITHM",
    "PACKAGE_VERSION",
    "SampleRecord",
    "StagedWrite",
    "blake2b_256",
    "canonical_json",
    "hash_file",
    "hash_text",
    "utc_now_iso",
    "versions_block",
]

PACKAGE_VERSION: Final = "0.1.0"

CONTENT_HASH_ALGORITHM: Final = "blake2b-256"
"""BLAKE2b with a 32-byte digest, named exactly as ``qecgen``'s DATA_CONTRACT.md names it.

Recorded alongside every digest rather than left implicit. The sibling project learned this
the hard way: the field was once called ``content_sha256`` while holding a BLAKE2b digest,
so anyone who verified it with SHA-256 concluded the file was corrupt.
"""


@dataclass(frozen=True, slots=True)
class SampleRecord:
    """One row's worth of source data, before any cleaning or encoding.

    ``label`` is a tuple rather than a bare int so that a multilabel source (FunSoC, 32
    binary targets) and a single-target source use one type. ``None`` means this source has
    no ground truth at all, which is a state the writer refuses rather than defaults to 0 --
    a missing label silently becoming the negative class is the exact failure this whole
    project is organised against.
    """

    sample_id: str
    sequence: str
    label: tuple[int, ...] | None
    source: str
    metadata: Mapping[str, Any] = field(default_factory=dict)


def blake2b_256(data: bytes) -> str:
    """The project's one digest. 32-byte BLAKE2b, hex."""
    return hashlib.blake2b(data, digest_size=32).hexdigest()


def hash_text(text: str) -> str:
    """Digest of text as UTF-8, so a digest never depends on the platform encoding."""
    return blake2b_256(text.encode("utf-8"))


def hash_file(path: Path) -> str:
    """Digest a file in chunks, so a multi-gigabyte raw input never enters memory."""
    digest = hashlib.blake2b(digest_size=32)
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_json(payload: Mapping[str, Any]) -> str:
    """JSON with sorted keys and no NaN, the one spelling used for every sidecar.

    ``allow_nan=False`` so a non-finite value can never be written as the bare token ``NaN``,
    which is not JSON and which ``json.load`` would nonetheless read back. ``sort_keys`` so
    that two runs of the same configuration produce byte-identical sidecars, which is what
    makes the reproducibility check in ``validation`` meaningful.
    """
    return json.dumps(payload, sort_keys=True, allow_nan=False, indent=2) + "\n"


def utc_now_iso() -> str:
    """Timestamp for the manifest's ``generated_at`` and nothing else.

    This is the only wall-clock value that reaches a written file. It is deliberately kept
    out of ``content_hash``, so that re-running the pipeline reproduces the digest even
    though the timestamp differs.
    """
    return datetime.now(UTC).isoformat()


def versions_block() -> dict[str, str]:
    """Versions of everything that can change an output byte."""
    import numpy
    import pandas

    return {
        "nexus_genomics": PACKAGE_VERSION,
        "python": sys.version.split()[0],
        "numpy": numpy.__version__,
        "pandas": pandas.__version__,
    }


class StagedWrite:
    """Write a group of files, or none of them.

    Mirrors ``qecgen.run.staged``. The table and its sidecar are one artefact: a ``.ml.csv``
    with no ``.ml.manifest.json`` beside it is reported as *not one of ours* rather than as
    corrupt, and that claim is only sound if a half-finished run can never produce one. Each
    file is written to a ``.tmp`` sibling and then moved into place with ``os.replace``,
    which is atomic per file; the group is committed in one pass at the end, so a crash before
    the commit leaves nothing. The commit loop itself is not a single transaction -- a kill
    between two replaces can leave one new file beside one old one -- but that state is
    detectable, because the sidecar's ``content_hash`` will not match the table.
    """

    def __init__(self) -> None:
        self._pending: list[tuple[Path, bytes]] = []

    def add_bytes(self, path: Path, data: bytes) -> None:
        self._pending.append((path, data))

    def add_text(self, path: Path, text: str) -> None:
        # newline handling is the caller's business: the CSV writer has already chosen its
        # terminator, and re-encoding here would turn every "\n" into "\r\n" on Windows.
        self._pending.append((path, text.encode("utf-8")))

    def commit(self) -> list[Path]:
        """Stage every file, then move them all into place."""
        staged: list[tuple[Path, Path]] = []
        try:
            for path, data in self._pending:
                path.parent.mkdir(parents=True, exist_ok=True)
                tmp = path.with_name(path.name + ".tmp")
                tmp.write_bytes(data)
                staged.append((tmp, path))
            for tmp, path in staged:
                os.replace(tmp, path)
        except BaseException:
            for tmp, _ in staged:
                tmp.unlink(missing_ok=True)
            raise
        return [path for _, path in staged]


def require_unique(ids: Sequence[str], what: str) -> None:
    """Refuse duplicate primary keys rather than let a join silently fan out.

    A duplicate ``sample_id`` is not a cosmetic problem: downstream it turns one row into
    several on every merge, and the resulting model trains on a reweighted dataset that
    nothing in the file records.
    """
    seen: set[str] = set()
    duplicates: list[str] = []
    for value in ids:
        if value in seen:
            duplicates.append(value)
        seen.add(value)
    if duplicates:
        sample = ", ".join(sorted(set(duplicates))[:5])
        raise ValueError(
            f"{what} contains {len(set(duplicates))} duplicated id(s) (e.g. {sample}). "
            f"Primary keys must be unique: a duplicate fans out on every join and silently "
            f"reweights the dataset."
        )
