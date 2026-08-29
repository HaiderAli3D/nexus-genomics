# Monday demo checklist

Everything below has been run and passes. This is the order to do it in on the day, and the
things worth saying out loud.

## Before the demo (5 minutes)

```powershell
cd C:\Projects\BioData
.venv\Scripts\python -m pytest -q                          # expect: 104 passed
.venv\Scripts\python -m nexus_genomics.cli validate outputs\*.ml.csv
```

Expect **6 files, 20–21 checks each, 0 failed** (124 checks in total). If the outputs are missing, rebuild them —
three of the four sources need no credentials:

```powershell
.venv\Scripts\python -m nexus_genomics.cli convert-all --profile demo
.venv\Scripts\python -m nexus_genomics.cli convert-all --profile full
```

Each profile writes three files; both are needed for the six the validate step expects.

Have open: `reports/source_audit.md` and this file.

## Running order

**1. Lead with the finding, not the code.** The brief described four
natural-versus-engineered datasets. None of them is one. Three still yield genuinely labelled
CSVs — just not with that label. Nothing was fabricated to fill the gap.

**2. Show the audit.**

```powershell
.venv\Scripts\python -m nexus_genomics.cli audit
```

Three sources ready, one blocked with the exact steps printed.

**3. Show a file.** `outputs\codontransformer_natural_vs_optimized_demo.ml.csv` — 992 rows ×
1,002 columns, opens in Excel on a double click.

```
sample_id,target,position_0001,position_0002,...,position_1000
ct_000000__natural,0,1,20,7,7,3,7,3,3,1,1,7,3,...
ct_000000__finetune,1,1,20,7,7,3,7,3,3,1,20,3,7,...
```

`1,20,7` is `A,T,G` — the start codon. The two rows are the *same protein in the same
organism*, encoded by nature and by the model.

**4. The verification worth showing.** Both members of a pair translate to the identical
protein. Across the 145 untruncated demo pairs, nucleotide identity runs **73% to 94%, median
81%** — that spread is the whole signal: synonymous codon choice. (Measured by decoding the
shipped file; `validate` does not check protein equivalence, it checks the file's integrity.)

```powershell
.venv\Scripts\python -m nexus_genomics.cli validate outputs\codontransformer_natural_vs_optimized_demo.ml.csv
```

**5. Reproducibility.** Two runs, identical bytes — including the one that re-downloads from
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
| `felix_guardian_part_role_*` | What kind of genetic part each detected sequence is | This is where the natural-vs-engineered dataset was *supposed* to come from. Explain why it isn't |

## The three things most likely to be asked

**"Why isn't there a natural-versus-engineered file?"**
Because at sequence level the FELIX corpus has one class. 1,670 sequences ship in the public
repo; only 280 can be linked to GUARDIAN's `engineered` call, and all 280 are *engineered* —
a detector that finds nothing emits no feature, so there is nothing for a negative to attach
to. The balanced 257/250 figure is per whole sample, and it is a detector prediction rather
than ground truth. The real ground truth is one browser download away; see
`reports/unresolved_questions.md` §4.

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
| Network error mid-convert | Zenodo/UniProt/Box unreachable | Cached raw files under `data\raw\` mean re-running is cheap; the outputs already on disk are valid |
| `is not a nexus-genomics dataset` | The `.ml.manifest.json` sidecar is missing | Rebuild that source; the group is staged, so a failed run leaves nothing half-written |
| Validation fails after hand-editing a CSV | Working as intended | The `content_hash` check catches edits |

## Deliverables on disk

- 6 Nexus CSVs + 12 sidecars in `outputs\`
- `reports\source_audit.md` — what each source really is
- `reports\validation_report.md` / `.json` — 124 checks, 0 failures
- `reports\unresolved_questions.md` — what is genuinely still open
- `README.md` — install, download, run, encoding, labels, limitations
- 104 tests, `ruff` clean, `mypy --strict` clean
