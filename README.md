# nexus-genomics

Convert public genomic datasets into a Nexus-shaped CSV (compatibility unverified), **without inventing labels**.

> **Nexus compatibility is not claimed.** The Nexus input format is not published. These
> files match the CSV shape described in the brief and mirror the conventions of the
> `ml_csv` exporter in the sibling `qecgen` project, which was built for a consumer that
> asked for `Primary Key, Target, Variable_1, …`. No fixture supplied by the Nexus team has
> been run against them. That claim will not be made until one passes.

## The task, and what it turned into

Four sources were handed over described as "natural versus engineered" genomic datasets.
**None of them is one.** Details and evidence are in [`reports/source_audit.md`](reports/source_audit.md);
the short version:

| Source | Reality | Honest target |
|---|---|---|
| Addgene / GEAC | 100% engineered plasmids. No natural class exists. | lab of origin, 1,314 classes |
| FELIX / GUARDIAN | Detector *predictions*, single-class at sequence level | genetic part role |
| SeqScreen / FunSoC | An analysis pipeline; "sequence of concern" ≠ "engineered" | 32 pathogenicity mechanisms |
| CodonTransformer | A codon-optimisation model over natural genes | natural vs model-generated |

Three of the four produce a genuinely labelled CSV. The fourth needs one gated download.
Nothing was fabricated to fill the gap.

## Install

Requires Python 3.13+. A project-local virtualenv is **required**, not merely tidy: the
optional `generate` extra pins `numpy<2.0.0`, and installing it globally would downgrade
numpy underneath anything else on the machine.

```powershell
python -m venv .venv
.venv\Scripts\python -m pip install -e ".[dev]"
```

## Quick start

```powershell
.venv\Scripts\python -m nexus_genomics.cli audit                    # what is available now
.venv\Scripts\python -m nexus_genomics.cli convert-all --profile demo
.venv\Scripts\python -m nexus_genomics.cli convert-all --profile full   # the three _full files
.venv\Scripts\python -m nexus_genomics.cli validate outputs\*.ml.csv
```

Each profile writes its own three files (`*_demo.ml.csv` and `*_full.ml.csv`), so run both to
reproduce all six.

Three of the four sources download what they need automatically — about 11 MB for the demo
profile, about 166 MB for `--profile full` — and need no credentials. The fourth reports
exactly what you must supply.

## Where each dataset comes from

| Source | How to get it | Size | Credentials |
|---|---|---|---|
| `codontransformer` | Streamed by HTTP Range from Zenodo 13262517 (6.08 GB file, never downloaded whole) | demo 5 MB · **full 160 MB** | none |
| `seqscreen` | Rice Box accession lists + UniProt REST, cached under `data\raw\` | 5 MB | none |
| `felix` | Repository tarball from GitHub | 1.3 MB | none |
| `addgene` | **You must download it** — see below | 647 MB | free DrivenData account |

Measured from the shipped sidecars, not estimated. The `full` CodonTransformer profile reads
40 windows of 4 MB; on a metered link that is the number to watch.

### You must supply: the GEAC files

1. Create a free account at <https://www.drivendata.org/accounts/signup/> and verify it.
2. Log in, open the [competition page](https://www.drivendata.org/competitions/63/genetic-engineering-attribution/)
   and **accept the competition rules** — accepting is the act that grants data access.
3. From the *Data* tab, download `train_values.csv` and `train_labels.csv`.
   The links are per-user pre-signed S3 URLs that expire after roughly 24 hours.
4. Put both in `data\raw\addgene\`.
5. `.venv\Scripts\python -m nexus_genomics.cli convert addgene --profile demo`

**The competition rules forbid redistributing the data or anything derived from it.** The
output CSV must stay on your machine; `.gitignore` already excludes `outputs/`.

### Optional: the full CodonTransformer file

The demo and full profiles both stream what they need. If you want the complete 622,915-pair
set, download
<https://zenodo.org/api/records/13262517/files/dataset_with_predictions_metrics.csv/content>
(6.08 GB, md5 `9d6571275fa7818f0a2eca04a02fed19`) to
`data\raw\codontransformer\dataset_with_predictions_metrics.csv` and raise `max_pairs`.

## Pipeline stages

```
adapter.load()  →  clean  →  encode  →  place  →  write  →  validate  →  report
```

| Stage | Module | What it does |
|---|---|---|
| Read the native format | `adapters/*.py` | One module per source; emits `SampleRecord`s |
| Clean | `cleaning.py` | Uppercase, strip whitespace, U→T, IUPAC policy, quarantine |
| Encode | `encoding.py` | A=1…Z=26, 0 for padding |
| Place | `encoding.place` | Pad, truncate or window to the configured length |
| Write | `nexus_csv.py` | The table plus two sidecars, committed atomically |
| Validate | `validation.py` | 20–21 checks per file |
| Report | `validation.render_markdown` / `cli.py` | Markdown and JSON, rendered from one object |

Everything after `adapter.load()` is source-agnostic, so all four sources obey identical
rules and a fix lands in one place.

## Encoding rules

```
A → 1    B → 2    C → 3   …   Z → 26
padding → 0
```

**These integers are categorical tokens. Their numerical distance has no biological meaning.**
`T=20` is not "ten times" `C=3`, and `G=7` does not sit "between" C and T in any chemical
sense. A model must treat these columns as categories — an embedding layer, or one-hot. A
model that treats them as magnitudes is learning from an artefact of the Latin alphabet. This
warning is repeated in every `.ml.manifest.json` because it is the easiest way to misuse
these files. Note it is *not* in the `.ml.samples.csv`, so if you hand someone only the two
CSVs, hand them this paragraph too.

`0` is reserved for padding and nothing else, so a zero in a feature column always means "no
residue here".

### Nucleotide, codon or amino acid?

The brief's own worked example settles the default: `1,3,7,20` decodes to `A,C,G,T` and
`1,20,7,20` to `A,T,G,T`. Both are valid DNA and neither is a plausible peptide. So
**nucleotide-position encoding is the default**.

But FunSoC sequences are proteins, so the model is *not* a translation mode:

- The **cipher is always** `A=1…Z=26`. Letters are letters.
- An `alphabet` field per dataset — `dna_nucleotide` or `amino_acid` — records what the
  letters *mean*, and is validated against the expected letter set.
- `codon` (three nucleotides per column) is **defined and refused** with an explanation: it
  needs a reading frame, which a plasmid does not have a single one of.
- `amino_acid_translated` is **defined and refused**: translating a whole plasmid end to end
  would emit a protein that does not exist.

**A known limitation.** The alphabet check catches a protein mislabelled as DNA — E, Q and I
are not bases. It *cannot* catch the reverse, because every nucleotide letter is also a valid
amino-acid letter (A alanine, C cysteine, G glycine, T threonine, N asparagine). So the
cleaning report counts `protein_sequences_using_only_nucleotide_letters` and surfaces it in
the sidecar, turning an otherwise invisible error into a visible one.

### Ambiguity codes

IUPAC codes take their own cipher value (`N`→14, `R`→18, …) rather than being flattened to N.
That keeps the encoding reversible and preserves the difference between "any base" and
"purine". `U→T` is opt-in per source and only legal for a nucleotide alphabet — U is
selenocysteine in a protein, and converting it would rewrite a real residue as threonine.
Anything outside A–Z is quarantined, never coerced.

## Label mappings

| Dataset | Target | Values |
|---|---|---|
| `codontransformer_natural_vs_optimized` | `target` | `0` natural NCBI CDS · `1` CodonTransformer-generated |
| `funsoc_pathogenicity_mechanism` | `target_00 … target_31` | `1` curated positive · `0` **not curated positive** (not a verified negative) |
| `felix_guardian_part_role` | `target` | 8 classes: CDS, engineered_region, engineered_tag, five_prime_UTR, origin_of_replication, promoter, replicon, terminator |
| `addgene_lab_attribution` | `target` | one class per lab ID; `I7FXTVDP` is the pooled **"Unknown Engineered"** class — unknown *laboratory*, still engineered |

The full class map and the meaning of every value are in each file's
`.ml.manifest.json` under `label_semantics`.

## Padding, truncation and windowing

Configured in `config/default.yaml`:

- `window.length` — 1,000 positions by default, giving 1,002 columns. Raise it when Nexus's
  real column ceiling is known.
- `window.policy` — `pad_or_truncate` gives exactly one row per sequence and **records how
  much was discarded, per sample**; `window` gives one row per window and discards nothing,
  padding the final window rather than dropping it.
- `window.stride`, `window.max_windows_per_sample` — the cap is applied deterministically and
  counted as `windows_dropped_by_cap`.

Windows get `sample_id` of the form `parent__w00`, and every one carries `parent_sample_id`,
`window_start` and `window_end` (one-based, inclusive) in the samples sidecar — so windows of
one sequence can be kept on one side of a train/test split.

## Output files

Each dataset is three files sharing a stem:

```
<name>.ml.csv             one header row, no '#' lines; sample_id, target(s), position_*
<name>.ml.manifest.json   format/version/columns/manifest — the file's magic line
<name>.ml.samples.csv     per-sample provenance, deliberately kept out of the model table
```

**The manifest sidecar is what makes a file provably ours.** A bare `.ml.csv` is
shape-identical to any other one-header-row CSV. The group is staged and then committed with
one `os.replace` per file, so a failure *before* the commit leaves nothing behind. The
replaces themselves are not one transaction: a kill between them can leave a new table beside
the previous run's sidecar, which the `content_hash` check then catches.

Produced and validated (2026-08-29):

| File | Rows | Cols | Size | Classes |
|---|---|---|---|---|
| `codontransformer_natural_vs_optimized_demo.ml.csv` | 992 | 1,002 | 2.2 MB | 496 / 496 |
| `codontransformer_natural_vs_optimized_full.ml.csv` | 31,730 | 1,002 | 70.4 MB | 15,865 / 15,865 |
| `funsoc_pathogenicity_mechanism_demo.ml.csv` | 583 | 1,033 | 1.4 MB | 32 labels, 12–70 positives each |
| `funsoc_pathogenicity_mechanism_full.ml.csv` | 13,164 | 1,033 | 29.8 MB | 32 labels, 23–4,549 positives each |
| `felix_guardian_part_role_demo.ml.csv` | 400 | 1,002 | 0.9 MB | 8 classes |
| `felix_guardian_part_role_full.ml.csv` | 481 | 1,002 | 1.0 MB | 8 classes |

`addgene_lab_attribution_*` is not listed because it needs the gated download.

Column naming: `position_0001 … position_1000`, **one-based** (sequence coordinates are
one-based in biology) and **zero-padded** (unpadded, `sorted(df.columns)` puts `position_10`
before `position_2` and silently permutes the feature matrix). The width is a property of
each file, which is why the sidecar publishes the literal names rather than describing them:

```python
import json, pandas as pd

side = json.load(open("outputs/codontransformer_natural_vs_optimized_demo.ml.manifest.json"))
df = pd.read_csv("outputs/codontransformer_natural_vs_optimized_demo.ml.csv")
X = df[side["columns"]["feature_columns"]]
y = df[side["columns"]["target_columns"][0]]
```

## How to verify the outputs

```powershell
.venv\Scripts\python -m pytest -q                        # 104 tests
.venv\Scripts\python -m ruff check . ; .venv\Scripts\python -m ruff format --check .
.venv\Scripts\python -m mypy src tests                   # strict
.venv\Scripts\python -m nexus_genomics.cli validate outputs\*.ml.csv
```

Reproducibility — two runs of the same configuration must produce identical bytes:

```powershell
.venv\Scripts\python -m nexus_genomics.cli convert-all --profile demo --out-dir outputs\rerun
.venv\Scripts\python -m nexus_genomics.cli reproduce `
    outputs\codontransformer_natural_vs_optimized_demo.ml.csv `
    outputs\rerun\codontransformer_natural_vs_optimized_demo.ml.csv
```

Each file gets 20–21 checks (21 when it has a single target column), including several that are easy to skip: the padding invariant
(every row must be residues then padding, never a zero before a non-zero), a decode round
trip of randomly sampled rows against a per-sample digest, sidecar-versus-table agreement,
recomputation of the `content_hash`, and a refusal to ship an all-zero target column.
Duplicate sequences, conflicting labels on identical sequences and class balance are
*reported* rather than failed — they are facts a reader must know, not defects.

Reports land in `reports/validation_report.md` and `reports/validation_report.json`, both
rendered from the same object so they cannot disagree.

## Self-review

After everything was working, the project was reviewed as if it were someone else's work:
four independent reviewers by failure mode, then an adversarial verifier per finding whose
job was to refute it. 43 findings, 28 confirmed, 10 partial, 5 refuted; all confirmed ones
fixed. The most valuable one was that the round-trip check had been comparing the encoder
against itself and so could not have caught an off-by-one. See
[`reports/self_review.md`](reports/self_review.md).

## Known limitations

1. **No general natural-versus-engineered benchmark exists in these four sources.** The
   closest thing, the FELIX engineering-signature FASTA, needs a manual browser download.
2. **CodonTransformer's positive class is one generator.** Both members of a pair encode the
   same protein, so the only signal is synonymous codon choice. The *source* is restricted to
   the 17 fine-tuned organisms; the shipped files contain fewer still — 3 in the demo and 9 in
   the full file, with Homo sapiens alone at 36% of rows. It measures whether a model can
   recognise CodonTransformer, on a narrow slice of life.
3. **FunSoC is positive-unlabelled.** A `0` means "not curated positive". Treating the zeros
   as true negatives will overstate precision. Class sizes range from 23 to 4,549.
4. **FELIX part roles are pipeline annotations**, not experimentally verified ground truth,
   and two classes have one or two members.
5. **Truncation at 1,000 positions is substantial** — 23.6 M residues in the full
   CodonTransformer file. It is recorded per sample and is symmetric between the two classes
   there (paired sequences are equal length), so it leaks no label; but raise `window.length`
   or switch to `window` policy if you need whole sequences.
6. **Sampling from the big Zenodo file is strided, not uniform random.** Rows are grouped by
   organism, so a byte prefix would sample one organism; the offsets used are fixed and
   recorded in the sidecar, but the sample is not statistically uniform.
7. **Nexus compatibility is unverified**, as stated at the top.

## Repository layout

```
config/default.yaml              every setting; its digest is recorded in each output
src/nexus_genomics/
  cli.py          Typer; every command prints its resolved config first
  common.py       SampleRecord, blake2b-256, atomic staged write
  cleaning.py     the cleaning rules and the refusals
  encoding.py     the cipher and the length policies
  nexus_csv.py    the writer, the sidecar, and the reader that refuses foreign files
  validation.py   20-21 checks, as data rather than prose
  adapters/       one module per source, registered in ADAPTERS
tests/                           104 tests, no network required
data/raw/                        downloads; never committed
outputs/                         the CSVs and sidecars
reports/                         source audit, validation, self-review, open questions, demo checklist
```

Conventions follow the sibling `qecgen` project: exact-pinned dependencies, `ruff` at line
length 100, `mypy --strict`, flat `tests/`, and docstrings that explain **why a trap exists**
rather than what the code does.
