# Email Compliance Remediation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Produce conservative single-target exports and an executable compliance report for every requirement that can be established without unavailable private data or a Nexus fixture.

**Architecture:** A source-agnostic projection module converts an `n_targets > 1` `LoadedSource` into deterministic one-target views before the existing shared conversion path. A separate compliance module inspects manifests, validation results, source availability, and explicit scientific limitations; the CLI only orchestrates both units.

**Tech Stack:** Python 3.13, Typer, dataclasses, pytest, pandas-backed existing validator, Ruff, mypy strict.

**Spec:** `docs/superpowers/specs/2026-08-29-email-compliance-design.md`

**Status (2026-08-29): IMPLEMENTATION COMPLETE; ARTIFACT REGENERATION BLOCKED.** All code,
tests, final-review remediation, and locally executable quality gates are complete. The required
CodonTransformer/FunSoC regeneration could not finish because this environment prohibits outbound
Zenodo/UniProt sockets, so receipts remain incomplete and artifact validation/compliance correctly
fail closed. Separate standing external conditions remain: Addgene requires authorized gated
files, Nexus requires an authoritative fixture, and the requested general target is unsupported
by the four original sources. The fifth FELIX-signatures adapter supplies only a narrower
recorded-event contrast.

## Global Constraints

- Never fabricate natural-versus-engineered labels.
- Never claim Nexus compatibility without an authoritative passing fixture.
- Preserve existing native outputs and their filenames.
- Addgene data and derivatives remain local and are never redistributed.
- New byte-affecting defaults belong in `config/default.yaml` and manifests record resolved layout.
- Commit and push steps are outside this implementation plan.

---

### Task 1: Generic Single-Target Projection

**Files:**
- Create: `src/nexus_genomics/projection.py`
- Modify: `src/nexus_genomics/adapters/base.py`
- Modify: `src/nexus_genomics/adapters/seqscreen.py`
- Test: `tests/test_projection.py`

**Interfaces:**
- Produces: `TargetProjection(index: int, name: str, slug: str, source: LoadedSource)`.
- Produces: `single_target_projections(source: LoadedSource) -> list[TargetProjection]`.
- Produces: `projected_path(base_path: Path, slug: str) -> Path`.
- Extends: `LoadedSource.target_names: tuple[str, ...] = ()`.

- [x] **Step 1: Write failing projection tests**

Test that a two-target source yields two `n_targets == 1` sources, preserves each bit exactly,
records `target_projection` metadata, creates deterministic safe slugs, and refuses missing or
duplicate target names. Test that a one-target source is returned as one unsuffixed projection.

- [x] **Step 2: Verify the tests fail**

Run: `.venv\Scripts\python -m pytest tests\test_projection.py -q`
Expected: import failure because `nexus_genomics.projection` does not exist.

- [x] **Step 3: Implement the minimal projection module**

Slug names with lowercase ASCII alphanumerics joined by `_`; reject empty and duplicate slugs.
For each projected record replace `(v0, ..., vn)` with `(vi,)`, preserving all other fields.
Require every label vector to have exactly `source.n_targets` values. Copy manifest metadata and
add original index, name, and target count without mutating the original source.

- [x] **Step 4: Publish target names from SeqScreen**

Set `target_names=tuple(funsocs)` alongside `target_index_to_funsoc`; leave the latter in the
manifest for compatibility. Other adapters use the dataclass default because they have one target.

- [x] **Step 5: Run focused and regression tests**

Run: `.venv\Scripts\python -m pytest tests\test_projection.py tests\test_adapters.py -q`
Expected: all pass.

### Task 2: Strict Export CLI

**Files:**
- Modify: `config/default.yaml`
- Modify: `src/nexus_genomics/cli.py`
- Test: `tests/test_cli.py`

**Interfaces:**
- Adds config: `output.target_layout: native`.
- Adds options: `--strict-single-target/--native-targets` to `convert` and `convert-all`.
- Consumes: `single_target_projections` and `projected_path` from Task 1.

- [x] **Step 1: Write failing CLI tests**

Use a monkeypatched two-target adapter and assert strict conversion writes two tables whose
headers each contain only `target`. Assert native conversion remains one two-target table. Assert
`convert-all --strict-single-target` forwards the mode, and config includes the native default.

- [x] **Step 2: Verify the tests fail**

Run: `.venv\Scripts\python -m pytest tests\test_cli.py -q`
Expected: unrecognized strict option or missing output configuration.

- [x] **Step 3: Implement strict conversion orchestration**

Load each adapter once. In native mode call `convert` exactly as before. In strict mode project
multilabel sources and call the same `convert` once per projection, using a deterministic suffix.
Merge `target_layout: single` into extra manifest data and report each written path. Do not alter
existing native filenames or default behavior.

- [x] **Step 4: Wire `convert-all` and resolved output logging**

Forward the strict flag to `convert_source`; print the resolved target layout in `_echo_config`.
Preserve exit 0 for gated downloads and exit 1 for all other failures.

- [x] **Step 5: Run focused tests**

Run: `.venv\Scripts\python -m pytest tests\test_cli.py tests\test_pipeline.py -q`
Expected: all pass.

### Task 3: Executable Email Compliance Audit

**Files:**
- Create: `src/nexus_genomics/compliance.py`
- Modify: `src/nexus_genomics/cli.py`
- Create: `tests/test_compliance.py`
- Modify: `tests/test_cli.py`

**Interfaces:**
- Produces: `ComplianceStatus` string enum with `PASS`, `BLOCKED_EXTERNAL`, `UNSUPPORTED_BY_SOURCE`, `FAIL`.
- Produces: `ComplianceCheck`, `ComplianceReport`, `audit_email_compliance(...)` and Markdown rendering.
- Adds CLI: `compliance --profile demo --raw-dir data/raw --out-dir outputs/strict --reports-dir reports`.

- [x] **Step 1: Write failing report tests**

Build valid temporary strict tables through `convert`. Assert structural requirements pass, a
ready source with no output fails, unavailable Addgene is externally blocked, the unsupported
general binary target is explicit, and Nexus acceptance is externally blocked. Assert JSON and
Markdown carry identical status counts.

- [x] **Step 2: Verify the tests fail**

Run: `.venv\Scripts\python -m pytest tests\test_compliance.py -q`
Expected: import failure because `nexus_genomics.compliance` does not exist.

- [x] **Step 3: Implement data-first compliance checks**

Inspect strict manifests and call the existing validator. Check one unique `sample_id`, exactly one
target, all declared `position_*` columns, letter encoding, `A=1 .. Z=26`, and internal validation.
Match outputs by each adapter's long output stem. Missing output is `BLOCKED_EXTERNAL` only when
adapter availability says gated; it is `FAIL` when a ready source was not exported.

- [x] **Step 4: Encode non-computable requirement truthfully**

Report the general natural-versus-engineered request as `UNSUPPORTED_BY_SOURCE`: none of the four
original sources supplies it, CodonTransformer is a one-generator model-artifact contrast, and the
fifth FELIX-signatures source is only a narrower recorded-event contrast. Report Nexus acceptance
as `BLOCKED_EXTERNAL` until an authoritative fixture exists. These statuses produce an
`INCOMPLETE` overall status, never a false green `PASS`.

- [x] **Step 5: Add CLI rendering and exit behavior**

Write `reports/email_compliance_report.json` and `.md` from one `ComplianceReport`. Print a compact
table. Exit 1 only when any status is `FAIL`; external and source limitations remain visible but do
not turn an otherwise correct implementation into a software failure.

- [x] **Step 6: Run focused tests**

Run: `.venv\Scripts\python -m pytest tests\test_compliance.py tests\test_cli.py -q`
Expected: all pass.

### Task 4: Artifacts, Documentation, And Adversarial Review

**Files:**
- Modify: `README.md`
- Modify: `reports/monday_demo_checklist.md`
- Generate: `outputs/strict/*.ml.csv` and sidecars
- Generate: `reports/email_compliance_report.json`
- Generate: `reports/email_compliance_report.md`

**Interfaces:**
- Consumes the strict CLI and compliance command from Tasks 2 and 3.
- Produces reviewed demo artifacts and final requirement evidence.

- [x] **Step 1: Update operational documentation**

Document native versus strict layouts, the 32 FunSoC projections, compliance status meanings,
commands, and the unchanged external blockers. Do not call either format Nexus-compatible.

- [ ] **Step 2: Generate strict demo artifacts — BLOCKED by prohibited outbound sockets**

Run: `.venv\Scripts\python -m nexus_genomics.cli convert-all --profile demo --out-dir outputs\strict --strict-single-target`
Expected: CodonTransformer, FELIX part roles, and FELIX signatures each write one table, FunSoC
writes 32, and Addgene is reported gated.

- [ ] **Step 3: Validate all strict artifacts and write compliance reports — BLOCKED on Step 2**

Run: `.venv\Scripts\python -m nexus_genomics.cli validate outputs\strict\*.ml.csv`
Run: `.venv\Scripts\python -m nexus_genomics.cli compliance --profile demo --out-dir outputs\strict`
Expected: every generated file passes validation; compliance is `INCOMPLETE` with no `FAIL` and
retains explicit blockers.

- [x] **Step 4: Run complete deterministic quality gates**

Run: `.venv\Scripts\python -m pytest -q`
Run: `.venv\Scripts\python -m ruff check .`
Run: `.venv\Scripts\python -m ruff format --check .`
Run: `.venv\Scripts\python -m mypy src tests`
Expected: all commands exit 0.

- [x] **Step 5: Run four independent reviews in parallel**

Review requirement traceability, scientific label validity, data integrity/reproducibility, and
operational/test readiness. Verify each reported finding against code and artifacts; fix every
confirmed actionable issue with a focused regression test.

- [ ] **Step 6: Re-run all affected checks after review fixes — artifact gates BLOCKED on Step 2**

Repeat the complete quality gates, strict validation, and compliance command. Report external
blockers separately from implementation completion.
