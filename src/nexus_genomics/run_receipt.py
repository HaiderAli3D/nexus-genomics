"""Atomic evidence that one bulk conversion invocation completed and produced exact tables."""

from __future__ import annotations

import json
from collections.abc import Sequence
from pathlib import Path
from typing import Any, Final

from nexus_genomics.common import PACKAGE_VERSION, StagedWrite, canonical_json
from nexus_genomics.nexus_csv import manifest_path

__all__ = [
    "RECEIPT_FILENAME",
    "RECEIPT_FORMAT",
    "RECEIPT_VERSION",
    "build_conversion_run_receipt",
    "receipt_path",
    "write_conversion_run_receipt",
]

RECEIPT_FILENAME: Final = ".conversion-run.json"
RECEIPT_FORMAT: Final = "nexus-genomics-conversion-run"
RECEIPT_VERSION: Final = 1


def receipt_path(out_dir: Path) -> Path:
    return out_dir / RECEIPT_FILENAME


def _generated_tables(out_dir: Path, table_paths: Sequence[Path]) -> list[dict[str, str]]:
    root = out_dir.resolve()
    entries: list[dict[str, str]] = []
    for table_path in table_paths:
        resolved = table_path.resolve()
        try:
            relative = resolved.relative_to(root)
        except ValueError as exc:
            raise ValueError(
                f"generated table {table_path} is outside output directory {out_dir}"
            ) from exc
        payload: dict[str, Any] = json.loads(manifest_path(table_path).read_text(encoding="utf-8"))
        content_hash = payload.get("manifest", {}).get("content_hash")
        if not isinstance(content_hash, str):
            raise ValueError(f"{manifest_path(table_path)} has no string content_hash")
        entries.append({"path": relative.as_posix(), "content_hash": content_hash})
    return sorted(entries, key=lambda entry: entry["path"])


def build_conversion_run_receipt(
    out_dir: Path,
    *,
    complete: bool,
    profile: str,
    target_layout: str,
    config_hash: str,
    enabled_source_keys: Sequence[str],
    successfully_converted_source_keys: Sequence[str],
    externally_gated_source_keys: Sequence[str],
    table_paths: Sequence[Path],
) -> dict[str, Any]:
    """Describe one invocation; incomplete receipts deliberately certify no prior run."""
    return {
        "format": RECEIPT_FORMAT,
        "version": RECEIPT_VERSION,
        "generator": f"nexus-genomics {PACKAGE_VERSION}",
        "complete": complete,
        "profile": profile,
        "target_layout": target_layout,
        "config_hash": config_hash,
        "enabled_source_keys": list(enabled_source_keys),
        "successfully_converted_source_keys": list(successfully_converted_source_keys),
        "externally_gated_source_keys": list(externally_gated_source_keys),
        "generated_tables": _generated_tables(out_dir, table_paths),
    }


def write_conversion_run_receipt(out_dir: Path, payload: dict[str, Any]) -> Path:
    """Replace the receipt atomically so a reader sees either the old or complete new JSON."""
    destination = receipt_path(out_dir)
    staged = StagedWrite()
    staged.add_text(destination, canonical_json(payload))
    staged.commit()
    return destination
