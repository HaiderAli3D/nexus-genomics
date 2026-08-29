"""Every check the brief asks for, expressed as data rather than as prose.

One :class:`ValidationResult` is produced per file. ``reports/validation_report.json`` is
that object serialised and ``reports/validation_report.md`` is rendered from the same
object, so the machine-readable and human-readable reports cannot drift apart.

A check that merely *reports* something (class balance, duplicate sequences) is recorded
with ``passed=True`` and its findings in ``detail``; only a check that establishes a defect
fails. The distinction matters: a class imbalance is a fact the reader must know, not a bug.
"""

from __future__ import annotations

import random
from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from nexus_genomics.common import CONTENT_HASH_ALGORITHM, blake2b_256, hash_file
from nexus_genomics.encoding import MAX_TOKEN, PAD_TOKEN, TOKEN_TO_LETTER
from nexus_genomics.nexus_csv import (
    INDEX_COLUMN,
    manifest_path,
    read_sidecar_and_shape,
    samples_path,
)

__all__ = ["Check", "ValidationResult", "validate_file"]


@dataclass(slots=True)
class Check:
    name: str
    passed: bool
    detail: str


@dataclass(slots=True)
class ValidationResult:
    path: str
    checks: list[Check] = field(default_factory=list)
    stats: dict[str, Any] = field(default_factory=dict)

    @property
    def ok(self) -> bool:
        return all(c.passed for c in self.checks)

    def add(self, name: str, passed: bool, detail: str) -> None:
        self.checks.append(Check(name, passed, detail))

    def to_json_dict(self) -> dict[str, Any]:
        return {
            "path": self.path,
            "ok": self.ok,
            "n_checks": len(self.checks),
            "n_failed": sum(1 for c in self.checks if not c.passed),
            "checks": [asdict(c) for c in self.checks],
            "stats": self.stats,
        }


_LETTER_LUT = np.array([""] + [TOKEN_TO_LETTER[i] for i in range(1, MAX_TOKEN + 1)], dtype="<U1")
"""Token -> letter, with index 0 mapping to the empty string so padding vanishes on join.

Only used after `tokens_within_configured_range` has passed, so every index is in bounds.
"""


def _decode_matrix(values: np.ndarray[Any, Any]) -> list[str]:
    """Decode a whole table at once.

    The per-row Python loop this replaces did 31.7 million dict lookups on the full
    CodonTransformer table and dominated the runtime of the entire validator.
    """
    return ["".join(row) for row in _LETTER_LUT[values]]


def validate_file(path: Path, *, sample_size: int = 25, seed: int = 0) -> ValidationResult:
    """Run every check against one written table."""
    result = ValidationResult(path=str(path))

    # --- loads at all -----------------------------------------------------------------
    try:
        frame = pd.read_csv(path)
        result.add(
            "loads_with_pandas_no_arguments",
            True,
            f"pandas.read_csv returned {frame.shape[0]} rows x {frame.shape[1]} columns "
            f"with no extra arguments",
        )
    except Exception as exc:
        result.add("loads_with_pandas_no_arguments", False, f"{type(exc).__name__}: {exc}")
        return result

    # --- sidecar agreement ------------------------------------------------------------
    try:
        payload, header, n_rows, widths = read_sidecar_and_shape(path)
    except Exception as exc:
        result.add("sidecar_present_and_matches_header", False, f"{type(exc).__name__}: {exc}")
        return result
    columns = payload["columns"]
    manifest = payload["manifest"]
    feature_cols: list[str] = list(columns["feature_columns"])
    target_cols: list[str] = list(columns["target_columns"])
    result.add(
        "sidecar_present_and_matches_header",
        True,
        f"{manifest_path(path).name} names {len(header)} columns and the header has {len(header)}",
    )

    result.add(
        "sidecar_counts_match_manifest",
        len(feature_cols) == manifest["n_features"]
        and len(target_cols) == manifest["n_targets"]
        and n_rows == manifest["n_samples"],
        f"sidecar lists {len(feature_cols)} feature / {len(target_cols)} target columns and "
        f"{n_rows} rows; manifest declares {manifest['n_features']} / "
        f"{manifest['n_targets']} / {manifest['n_samples']}",
    )

    # --- primary key ------------------------------------------------------------------
    has_pk = INDEX_COLUMN in frame.columns
    result.add("primary_key_column_exists", has_pk, f"looked for {INDEX_COLUMN!r}")
    if not has_pk:
        return result
    # The initial no-arguments load proves ordinary pandas compatibility, but its inference turns
    # an exact textual key such as ``001`` into integer ``1``. Re-read just the key as text before
    # any identity or sidecar comparison so validation does not manufacture a mismatch.
    pk = pd.read_csv(path, usecols=[INDEX_COLUMN], dtype={INDEX_COLUMN: "string"})[INDEX_COLUMN]
    pk_values = pk.astype(str).tolist()
    result.add("primary_key_non_null", bool(pk.notna().all()), f"{int(pk.isna().sum())} null")
    n_dupe = int(len(pk) - pk.nunique())
    result.add(
        "primary_key_unique",
        n_dupe == 0,
        "all unique" if n_dupe == 0 else f"{n_dupe} duplicated primary key(s)",
    )

    # --- complete samples-sidecar integrity -------------------------------------------
    sample_sidecar = samples_path(path)
    sample_frame: pd.DataFrame | None = None
    sample_error = ""
    try:
        sample_frame = pd.read_csv(sample_sidecar, dtype={INDEX_COLUMN: str})
    except Exception as exc:
        sample_error = f"{type(exc).__name__}: {exc}"

    declared_samples_hash = manifest.get("samples_content_hash")
    declared_samples_algorithm = manifest.get("samples_content_hash_algorithm")
    current_samples_hash = hash_file(sample_sidecar) if sample_frame is not None else None
    samples_hash_ok = (
        isinstance(declared_samples_hash, str)
        and declared_samples_algorithm == CONTENT_HASH_ALGORITHM
        and current_samples_hash == declared_samples_hash
    )
    result.add(
        "samples_sidecar_content_hash_matches",
        samples_hash_ok,
        (
            f"recomputed {current_samples_hash[:16]}... vs manifest "
            f"{declared_samples_hash[:16]}... using {declared_samples_algorithm}"
            if current_samples_hash is not None and isinstance(declared_samples_hash, str)
            else (
                "required samples sidecar hash is unavailable: "
                f"{sample_error or 'manifest fields missing'}"
            )
        ),
    )

    sample_ids: list[str] = []
    if sample_frame is not None and INDEX_COLUMN in sample_frame.columns:
        sample_ids = sample_frame[INDEX_COLUMN].astype(str).tolist()
    declared_sample_rows = manifest.get("samples_n_rows")
    samples_count_ok = (
        sample_frame is not None
        and len(sample_frame) == declared_sample_rows
        and declared_sample_rows == manifest.get("n_samples")
    )
    result.add(
        "samples_sidecar_row_count_matches_manifest",
        samples_count_ok,
        f"sidecar has {len(sample_frame) if sample_frame is not None else 'unreadable'} rows; "
        f"samples_n_rows={declared_sample_rows!r}, n_samples={manifest.get('n_samples')!r}",
    )
    samples_unique = bool(sample_ids) and len(sample_ids) == len(set(sample_ids))
    result.add(
        "samples_sidecar_sample_ids_unique",
        samples_unique,
        "all unique"
        if samples_unique
        else f"{len(sample_ids) - len(set(sample_ids))} duplicate(s), or sample_id missing",
    )
    sample_ids_match = sample_ids == pk_values
    result.add(
        "samples_sidecar_sample_ids_match_table_order",
        sample_ids_match,
        "exact ordered equality"
        if sample_ids_match
        else (
            f"table has {len(pk_values)} IDs and sidecar has {len(sample_ids)}; "
            "order/content differs"
        ),
    )

    # --- targets ----------------------------------------------------------------------
    missing_targets = [c for c in target_cols if c not in frame.columns]
    result.add(
        "target_columns_exist",
        not missing_targets,
        "present" if not missing_targets else f"missing {missing_targets}",
    )
    non_integer_targets = [
        column
        for column in target_cols
        if column in frame.columns and not pd.api.types.is_integer_dtype(frame[column])
    ]
    targets_are_integer = not missing_targets and not non_integer_targets
    result.add(
        "all_target_columns_integer",
        targets_are_integer,
        "all integer dtype"
        if targets_are_integer
        else f"non-integer or missing target columns: {non_integer_targets or missing_targets}",
    )
    documented = manifest.get("label_semantics")
    allowed: set[int] | None = None
    semantics_valid = True
    if documented is not None and (not isinstance(documented, Mapping) or not documented):
        semantics_valid = False
        result.add(
            "labels_only_documented_values",
            False,
            "label_semantics is present but does not contain a nonempty class mapping",
        )
    elif isinstance(documented, Mapping):
        try:
            allowed = {int(key) for key in documented}
        except (TypeError, ValueError):
            semantics_valid = False
            result.add(
                "labels_only_documented_values",
                False,
                "label_semantics contains a key that is not an integer class",
            )
    labels_usable = targets_are_integer and semantics_valid
    if allowed is not None and labels_usable:
        seen_values: set[int] = set()
        for col in target_cols:
            seen_values.update(frame[col].unique().tolist())
        undocumented = sorted(seen_values - allowed)
        result.add(
            "labels_only_documented_values",
            not undocumented,
            f"documented {sorted(allowed)}; found {sorted(seen_values)}"
            + (f"; UNDOCUMENTED {undocumented}" if undocumented else ""),
        )
    if not labels_usable:
        result.add(
            "label_dependent_checks_skipped",
            False,
            "target columns are not exact integer dtype or label semantics are malformed, "
            "so label comparisons, conflicts, class counts, and balance were not coerced",
        )

    # --- feature columns --------------------------------------------------------------
    missing_features = [c for c in feature_cols if c not in frame.columns]
    result.add(
        "all_position_columns_present",
        not missing_features,
        f"{len(feature_cols)} declared"
        + ("" if not missing_features else f"; missing {len(missing_features)}"),
    )
    if missing_features:
        return result

    features = frame[feature_cols]
    non_numeric = [c for c in feature_cols if not pd.api.types.is_integer_dtype(frame[c])]
    result.add(
        "all_position_columns_integer",
        not non_numeric,
        "all integer dtype"
        if not non_numeric
        else f"{len(non_numeric)} non-integer, e.g. {non_numeric[:5]}",
    )

    if non_numeric:
        # An empty cell makes pandas read the column as float64 with NaN, and the token-range
        # check below does int(values.min()) -- which raises rather than reporting. The
        # validator must never crash on exactly the corruption it exists to detect.
        result.add(
            "remaining_checks_skipped",
            False,
            f"{len(non_numeric)} feature column(s) are not integer (an empty or non-numeric "
            f"cell), so the token-range and decode-dependent checks could not run",
        )
        return result

    result.add(
        "every_row_has_the_same_width",
        widths == {len(header)},
        f"row widths observed: {sorted(widths)}; header has {len(header)}",
    )

    values = features.to_numpy()
    lo, hi = int(values.min()), int(values.max())
    in_range = lo >= PAD_TOKEN and hi <= MAX_TOKEN
    result.add(
        "tokens_within_configured_range",
        in_range,
        f"observed {lo}..{hi}, allowed {PAD_TOKEN}..{MAX_TOKEN}",
    )
    if not in_range:
        # Everything below decodes tokens back to letters, which is undefined for a value
        # outside the cipher. Returning here reports the defect; letting `decode` raise
        # would turn the validator into a crash on exactly the file it exists to catch.
        result.add(
            "remaining_checks_skipped",
            False,
            "token range check failed, so the decode-dependent checks (padding invariant, "
            "duplicates, label conflicts, round trip) could not run on this file",
        )
        return result

    # --- the padding invariant --------------------------------------------------------
    # Padding is appended, never interleaved. So every row must be a run of 1..26 followed
    # by a run of 0. A zero with a non-zero after it means either the encoder mapped a
    # biological character to 0, or padding landed in the middle -- both silently corrupt
    # the sequence while leaving a well-formed table.
    bad_rows: list[str] = []
    empty_rows: list[str] = []
    pad_total = 0
    for idx in range(values.shape[0]):
        row = values[idx]
        nonzero = row != PAD_TOKEN
        n_real = int(nonzero.sum())
        pad_total += int(row.shape[0] - n_real)
        if n_real == 0:
            empty_rows.append(pk_values[idx])
        if not bool(nonzero[:n_real].all()):
            bad_rows.append(pk_values[idx])
    result.add(
        "padding_only_ever_trails_the_sequence",
        not bad_rows,
        "every row is residues then padding"
        if not bad_rows
        else f"{len(bad_rows)} row(s) have a 0 before a non-zero, e.g. {bad_rows[:5]}",
    )
    if bad_rows:
        # decode() skips a PAD_TOKEN wherever it sits, so on a row with interleaved padding
        # it returns a wrong-but-plausible sequence rather than raising. Continuing would
        # print confirmatory passes for the duplicate, conflict and round-trip checks
        # underneath a failure that invalidates all three.
        result.add(
            "remaining_checks_skipped",
            False,
            "the padding invariant failed, so the decode-dependent checks (duplicates, "
            "label conflicts, round trip) could not run on this file",
        )
        return result
    result.add(
        "no_row_is_silently_empty",
        not empty_rows,
        "every row carries at least one residue"
        if not empty_rows
        else f"{len(empty_rows)} all-padding row(s), e.g. {empty_rows[:5]}",
    )
    result.add(
        "zero_is_reserved_for_padding",
        True,
        f"{pad_total} padding cells across {values.shape[0]} rows; no letter maps to 0 by "
        f"construction (the cipher starts at A=1)",
    )

    # --- reported, not judged ---------------------------------------------------------
    decoded = _decode_matrix(values)
    seq_counts = Counter(decoded)
    dupes = {s: n for s, n in seq_counts.items() if n > 1}
    result.add(
        "duplicate_sequences_reported",
        True,
        f"{len(dupes)} distinct sequence(s) appear more than once "
        f"({sum(dupes.values())} rows of {len(decoded)})",
    )

    # Extracted once, not per row per column. `frame.iloc[i]` inside the loop rebuilt a
    # 1,033-element object Series for every target column of every row: 87 seconds on the
    # 13,164-row FunSoC file against 0.02 seconds for the vectorised form.
    if labels_usable:
        label_matrix = frame[target_cols].to_numpy()
        label_of: dict[str, set[tuple[int, ...]]] = {}
        for i, seq in enumerate(decoded):
            key = tuple(label_matrix[i].tolist())
            label_of.setdefault(seq, set()).add(key)
        conflicting = {s: v for s, v in label_of.items() if len(v) > 1}
        result.add(
            "conflicting_labels_on_identical_sequences_reported",
            True,
            f"{len(conflicting)} sequence(s) carry more than one distinct label vector"
            + (
                ""
                if not conflicting
                else f"; first example has labels {sorted(next(iter(conflicting.values())))}"
            ),
        )

        # An all-zero target column cannot be learned and is indistinguishable, to a reader,
        # from a class that genuinely never occurs. It is almost always a capped sample that
        # dropped a class, so it fails rather than being reported.
        all_zero = [c for c in target_cols if int(frame[c].abs().sum()) == 0]
        result.add(
            "every_target_column_carries_at_least_one_positive",
            not all_zero,
            "every target column has a non-zero value"
            if not all_zero
            else f"{len(all_zero)} target column(s) are entirely zero: {all_zero[:8]}",
        )

        counts: dict[str, Any] = {}
        for col in target_cols:
            counts[col] = {str(k): int(v) for k, v in frame[col].value_counts().items()}
        result.stats["class_counts"] = counts
        if len(target_cols) == 1:
            vc = frame[target_cols[0]].value_counts()
            result.stats["n_distinct_classes"] = len(vc)
            # A constant target is a degenerate file. The all-zero check above only catches it
            # when the constant happens to be 0, and reporting balance = 1.0 for a single class
            # would describe the most useless possible dataset as perfectly balanced.
            result.add(
                "a_single_target_column_has_at_least_two_classes",
                len(vc) >= 2,
                f"{len(vc)} distinct value(s) in {target_cols[0]!r}: {sorted(vc.index.tolist())}"
                + (
                    ""
                    if len(vc) >= 2
                    else " -- a constant target cannot be learned and is almost always a "
                    "sampling cap that dropped every other class"
                ),
            )
            result.stats["class_balance_min_over_max"] = (
                round(float(vc.min() / vc.max()), 4) if len(vc) > 1 else None
            )

    result.stats["row_count"] = int(frame.shape[0])
    result.stats["column_count"] = int(frame.shape[1])
    result.stats["padding_cells"] = pad_total
    result.stats["content_hash"] = manifest["content_hash"]
    result.stats["content_hash_algorithm"] = manifest["content_hash_algorithm"]
    result.stats["file_digest"] = hash_file(path)
    for stat_key in ("cleaning", "window", "raw_input_hashes", "licence", "source_url", "alphabet"):
        if stat_key in manifest:
            result.stats[stat_key] = manifest[stat_key]

    # --- content hash actually matches the bytes on disk ------------------------------
    # Reuses the digest streamed in 1 MB chunks above. `path.read_bytes()` here would
    # allocate a second whole-file bytes object -- 4.38 GB at the largest size this
    # convention has shipped -- inside the checker whose job is to report problems rather
    # than die on them.
    recomputed = result.stats["file_digest"]
    result.add(
        "content_hash_matches_the_table_on_disk",
        recomputed == manifest["content_hash"],
        f"recomputed {recomputed[:16]}... vs manifest {manifest['content_hash'][:16]}...",
    )

    # --- round trip against the per-sample digests ------------------------------------
    result.add(
        *_round_trip(path, pk_values, decoded, sample_size=sample_size, seed=seed),
    )
    return result


def _round_trip(
    path: Path,
    sample_ids: Sequence[str],
    decoded: Sequence[str],
    *,
    sample_size: int,
    seed: int,
) -> tuple[str, bool, str]:
    """Decode randomly sampled rows and compare against the digest the adapter recorded.

    The samples sidecar stores ``source_slice_blake2b`` -- the digest of the CLEANED SOURCE
    slice that row was supposed to carry, taken before encoding. That is what makes this a
    real check rather than a tautology: hashing ``decode(tokens)`` on both sides would
    compare the encoder against itself, and an off-by-one in ``place`` would pass on every
    row of every file. Comparing digests rather than re-reading the source also means the
    check works without the raw download present, which matters for the gated sources.
    """
    name = "sampled_rows_decode_back_to_their_source_sequence"
    side = samples_path(path)
    if not side.exists():
        return (name, False, f"{side.name} is missing, so no round trip is possible")
    meta = pd.read_csv(side, dtype={INDEX_COLUMN: "string"})
    if "source_slice_blake2b" not in meta.columns:
        return (name, False, f"{side.name} has no source_slice_blake2b column")
    expected = dict(zip(meta[INDEX_COLUMN].astype(str), meta["source_slice_blake2b"], strict=True))

    rng = random.Random(seed)
    n = min(sample_size, len(decoded))
    picks = rng.sample(range(len(decoded)), n)
    bad: list[str] = []
    for i in picks:
        sid = sample_ids[i]
        want = expected.get(sid)
        got = blake2b_256(decoded[i].encode("utf-8"))
        if want != got:
            bad.append(sid)
    return (
        name,
        not bad,
        f"{n} row(s) sampled with seed {seed}; all decode to their recorded digest"
        if not bad
        else f"{len(bad)} of {n} sampled rows do not match, e.g. {bad[:5]}",
    )


def compare_runs(first: Path, second: Path) -> Check:
    """Two runs of the same configuration must produce the same table bytes."""
    a, b = hash_file(first), hash_file(second)
    return Check(
        "rerunning_with_the_same_config_is_byte_identical",
        a == b,
        f"{first.name}={a[:16]}...  {second.name}={b[:16]}...",
    )


def render_markdown(results: Sequence[ValidationResult], extra: Mapping[str, Any]) -> str:
    """The human-readable report, rendered from the same objects as the JSON one."""
    lines: list[str] = ["# Validation report", ""]
    total = len(results)
    passing = sum(1 for r in results if r.ok)
    lines += [f"{passing} of {total} file(s) pass every check.", ""]
    for extra_key, extra_value in extra.items():
        lines.append(f"- **{extra_key}**: {extra_value}")
    if extra:
        lines.append("")
    for r in results:
        lines += [f"## `{Path(r.path).name}`", ""]
        lines.append(
            f"**{'PASS' if r.ok else 'FAIL'}** — {len(r.checks)} checks, "
            f"{sum(1 for c in r.checks if not c.passed)} failed."
        )
        lines += ["", "| Check | Result | Detail |", "|---|---|---|"]
        for c in r.checks:
            detail = c.detail.replace("|", "\\|")
            lines.append(f"| `{c.name}` | {'pass' if c.passed else '**FAIL**'} | {detail} |")
        lines += ["", "<details><summary>Statistics</summary>", "", "```json"]
        import json

        lines.append(json.dumps(r.stats, indent=2, sort_keys=True))
        lines += ["```", "", "</details>", ""]
    return "\n".join(lines) + "\n"
