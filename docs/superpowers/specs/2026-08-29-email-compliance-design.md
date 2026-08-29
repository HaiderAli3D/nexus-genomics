# Email Compliance Remediation Design

## Goal

Make every code-addressable deliverable conform to the email's tabular contract without
inventing natural-versus-engineered labels or claiming Nexus acceptance without evidence.
The result must distinguish implemented compliance from external blockers.

## Non-Negotiable Boundaries

- Preserve each source's actual label semantics.
- Do not classify Addgene, GUARDIAN, or FunSoC records as natural or engineered when the
  source does not provide that ground truth.
- Keep `nexus_compatibility` as `NOT VERIFIED` until an authoritative Nexus fixture passes.
- Do not redistribute the gated Addgene data or derivatives.

## Export Contract

Existing validated outputs remain unchanged. A generic projection path will additionally
emit one single-target table for every target in a multilabel source. For FunSoC this means
32 binary files, each with `sample_id,target,position_0001,...,position_1000`. Projection
preserves every target bit and records the original target index, name, and source manifest
in each projected manifest. Single-target and multiclass sources continue to emit one table.

Projection is generic pipeline behavior rather than SeqScreen-specific adapter behavior.
Adapters remain responsible only for native-format loading and honest label semantics.

## Compliance Audit

Add a CLI command that evaluates the email requirements as data and writes Markdown and JSON.
Each requirement receives one of:

- `PASS`: established by executable evidence.
- `BLOCKED_EXTERNAL`: implementation exists but required gated data or Nexus fixture is absent.
- `UNSUPPORTED_BY_SOURCE`: the requested scientific label does not exist in the named source.
- `FAIL`: an actionable implementation requirement is violated.

The command checks output shape, primary-key integrity, a single target column in strict
exports, A=1 through Z=26 encoding metadata, one position per feature column, source coverage,
label provenance, and Nexus verification status. It exits nonzero only for `FAIL`; explicit
external or scientific blockers remain visible without masquerading as code failures.

## CLI And Configuration

Add an opt-in strict export flag to `convert` and `convert-all`. The default preserves current
filenames and output behavior. Under strict export, multilabel datasets produce one file per
target with deterministic, filesystem-safe names; other datasets use the current filename.
Configuration remains the source of byte-affecting defaults, and manifests record the mode.

## Error Handling

Refuse malformed target vectors, undocumented target indexes, duplicate projected filenames,
and projections with no rows. Missing Addgene files remain a gated availability result.
Absence of a Nexus fixture remains `BLOCKED_EXTERNAL`, never an inferred pass.

## Testing And Verification

Use test-driven development. Add focused tests for projection values, names, manifests,
single-target strict shape, deterministic output, CLI exit behavior, and the four compliance
statuses. Re-run the complete pytest, Ruff, format, strict mypy, validation, and reproducibility
checks. Generate strict demo outputs for all available sources and both compliance reports.

## Independent Review

After implementation, run parallel reviews for:

1. Email requirement traceability and Nexus-shape conservatism.
2. Scientific validity and label provenance.
3. Data integrity, validation, and reproducibility.
4. Operational readiness, CLI behavior, and missing-test analysis.

Resolve every verified code finding. The final report must state separately what passes, what
is unsupported by the source data, and what remains blocked on inputs unavailable in this pass.

## Success Criteria

- Every available strict export has one primary key, one target, and one token per position.
- No source receives a fabricated natural-versus-engineered label.
- FunSoC's 32 labels are preserved losslessly across 32 strict binary tables.
- Addgene's adapter remains ready for the authorized files and is reported as externally blocked.
- Nexus acceptance is reported as externally blocked until a fixture exists.
- All automated quality gates pass and independent reviews find no unresolved actionable issue.
