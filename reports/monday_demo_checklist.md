# Monday demo checklist

Everything below has been run and passes. This is the order to do it in on the day, and the
things worth saying out loud.

## Before the demo (5 minutes)

```powershell
cd C:\Projects\BioData
.venv\Scripts\python -m pytest -q                          # expect: 203 passed
.venv\Scripts\python -m nexus_genomics.cli validate (Get-ChildItem outputs\*.ml.csv)
.venv\Scripts\python -m nexus_genomics.cli validate (Get-ChildItem outputs\strict\*.ml.csv) --reports-dir reports\strict
.venv\Scripts\python -m nexus_genomics.cli compliance --profile demo --out-dir outputs\strict
```

Expect eight native files and 35 strict files, all with zero failed validation checks. Compliance's
headline must be `INCOMPLETE`, with zero `FAIL`, while Addgene/Nexus blockers and the unsupported
general target remain explicit. If outputs are missing or any prior run failed, rebuild them:

```powershell
.venv\Scripts\python -m nexus_genomics.cli convert-all --profile demo
.venv\Scripts\python -m nexus_genomics.cli convert-all --profile full
.venv\Scripts\python -m nexus_genomics.cli convert-all --profile demo --out-dir outputs\strict --strict-single-target
```

Native mode remains the default and preserves the six original tables. FELIX signatures adds one
table per profile, so both profiles now provide eight native tables. Strict mode is separate and
opt-in. Before reusing any artifact, require a complete matching `.conversion-run.json`, successful
validation, and the compliance report; an old table on disk is not evidence of a successful run.

Have open: `reports/source_audit.md`, `reports/email_compliance_report.md`, and this file.

## Running order

**1. Lead with the finding, not the code.** The brief described four
natural-versus-engineered datasets. None of those four supplies a general target of that kind.
A fifth source in a Supporting Information ZIP supplies a narrower recorded-event contrast:
host-context sequence not introduced by the event versus sequence introduced by it. It is useful
and source-supported, but not a general benchmark. Nothing was fabricated.

**2. Show the audit.**

```powershell
.venv\Scripts\python -m nexus_genomics.cli audit
```

Four sources ready, Addgene externally gated with the exact authorized-download steps printed.

**3. Show a file.** `outputs\codontransformer_natural_vs_optimized_demo.ml.csv` — 992 rows ×
1,002 columns, opens in Excel on a double click.

```
sample_id,target,position_0001,position_0002,...,position_1000
ct_000000__natural,0,1,20,7,7,3,7,3,3,1,1,7,3,...
ct_000000__finetune,1,1,20,7,7,3,7,3,3,1,20,3,7,...
```

`1,20,7` is `A,T,G` — the start codon. The two rows are the *same protein in the same
organism*, encoded by nature and by the model.

**4. Show strict single-target mode.** The strict demo contains 35 tables in the current registry:
one unsuffixed CodonTransformer table, one unsuffixed FELIX part-role table, one unsuffixed
FELIX-signatures table, and 32 named FunSoC tables. Every table has bare `target`; each FunSoC
manifest records the original target index and name.

This shape is Nexus-conservative, not verified Nexus-compatible. Projection changes layout only;
it does not turn FunSoCs, GUARDIAN roles, or Addgene labs into provenance labels.

**5. Show the compliance report.** The overall headline is `INCOMPLETE`, not an ambiguous pass.
At check level, `PASS` is executable evidence; `BLOCKED_EXTERNAL` covers Addgene and Nexus
verification; `UNSUPPORTED_BY_SOURCE` covers the requested general target in the original four
sources; `FAIL` is an actionable defect and is the only status that exits nonzero.

**6. The verification worth showing.** Both members of a pair translate to the identical
protein. Across the 145 untruncated demo pairs, nucleotide identity runs **73% to 94%, median
81%** — that spread is the whole signal: synonymous codon choice. (Measured by decoding the
shipped file; `validate` does not check protein equivalence, it checks the file's integrity.)

```powershell
.venv\Scripts\python -m nexus_genomics.cli validate outputs\codontransformer_natural_vs_optimized_demo.ml.csv
```

**7. Reproducibility.** Two runs, identical bytes — including the one that re-downloads from
Zenodo.

```powershell
.venv\Scripts\python -m nexus_genomics.cli convert-all --profile demo --out-dir outputs\rerun
.venv\Scripts\python -m nexus_genomics.cli reproduce `
    outputs\codontransformer_natural_vs_optimized_demo.ml.csv `
    outputs\rerun\codontransformer_natural_vs_optimized_demo.ml.csv
```

## What to say about each dataset

| File | One-line description | Say this too |
|---|---|---|
| `codontransformer_natural_vs_optimized_*` | Natural CDS vs the same protein codon-optimised by CodonTransformer | Not a general engineered detector — one generator, identical proteins, and only 9 organisms in the shipped full file |
| `funsoc_pathogenicity_mechanism_*` | 32 pathogenicity mechanisms over UniProt proteins | "Sequence of concern" is **not** "engineered". Most are natural toxins. Positive-only: a 0 is not a verified negative |
| `felix_natural_vs_engineered_*` | Narrow recorded-event contrast. 0 = host context not introduced by the event; 1 = sequence introduced by it | Class 0 is event-relative, not a global natural claim. `insertion.site` is excluded, and the manifest's length-only baseline must be quoted |
| `felix_guardian_part_role_*` | What kind of genetic part each detected sequence is | This is what the *public repo* can support. Explain why the binary had to come from the SI instead |

## The three things most likely to be asked

**"Where did the recorded-event labels come from?"**
The ACS Supporting Information, fetched by hand. The FELIX test-and-evaluation team recorded which
sequence a particular event introduced or removed. The emitted distinction is therefore
event-relative, not a claim of global naturalness. The *public* MIDOE repo cannot support a binary:
of its 1,670 sequences only 280 link to GUARDIAN's `engineered` call and all 280 are engineered —
one class, from detector predictions rather than experimental ground truth.

**"Is the FELIX-signatures set enough, and is it clean?"**
It is small and length-confounded, so quote the manifest's freshly measured length-only baseline
with any result and group splits by the shared element number. Rearrangements, point mutations,
and `insertion.site` are excluded where neither emitted provenance class is established; malformed
FASTA records are quarantined rather than guessed.

**"What do the numbers mean?"**
A=1 … Z=26, 0 is padding only. **They are categorical tokens** — T=20 is not ten times C=3.
Use an embedding or one-hot, never treat them as magnitudes. It says so in every
`.ml.manifest.json`.

**"Is this Nexus-compatible?"**
Unverified, and deliberately so. It matches the shape in the brief and mirrors the `.ml.csv`
conventions from qecgen, which was built for a consumer that asked for
`Primary Key, Target, Variable_1, …`. A fixture from the Nexus team would settle it. The
window length is 1,000 as a starting point and is one config key away from anything else.

## If something breaks

| Symptom | Cause | Fix |
|---|---|---|
| `addgene could not be converted` | Expected — gated download | Not a failure. It prints the steps |
| Network error mid-convert | Zenodo/UniProt/Box unreachable | Re-run when available; do not reuse outputs unless the complete receipt, validation, and compliance all match |
| `.conversion-run.json` is absent, incomplete, or mismatched | The run did not prove a complete current artifact set | Re-run `convert-all`; multi-file replacement can be partial even though validation detects mismatches |
| `is not a nexus-genomics dataset` | The `.ml.manifest.json` sidecar is missing | Rebuild the full selected run, then validate and run compliance |
| Validation fails after hand-editing a CSV | Working as intended | The `content_hash` check catches edits |

## Deliverables on disk

- 8 native nexus-genomics CSVs + 16 sidecars in `outputs\`
- 35 strict demo CSVs + 70 sidecars in `outputs\strict\`
- `reports\source_audit.md` — what each source really is
- `reports\validation_report.md` / `.json` — native checks, 0 failures after a complete run
- `reports\strict\validation_report.md` / `.json` — strict checks, 0 failures after a complete run
- `reports\email_compliance_report.md` / `.json` — 0 `FAIL`; blockers remain explicit
- `reports\unresolved_questions.md` — what is genuinely still open
- `README.md` — install, download, run, encoding, labels, limitations
- 203 tests, `ruff` clean, `mypy --strict` clean

**Current state (2026-08-29, post pre-email review):** a clean rebuild with network access to
Zenodo, UniProt, GitHub and Box all working produced all eight native tables and all 35 strict
tables in one run; `.conversion-run.json` is complete, validation shows 0 failed checks
everywhere, and compliance shows 0 `FAIL`. A pre-email adversarial review also found and fixed a
real parsing bug in the FELIX-signatures adapter (6 of 335 elements were misparsed and wrongly
excluded); the numbers in this checklist and in `README.md` reflect the corrected output.
