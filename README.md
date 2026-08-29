# nexus-genomics

Convert public genomic datasets into a Nexus-shaped CSV (compatibility unverified), **without inventing labels**.

New here? [`GUIDE.md`](GUIDE.md) is a five-minute, plain-language walkthrough of what this project
found and how to run it. This README is the complete technical reference underneath that.

> **Nexus compatibility is not claimed.** The Nexus input format is not published. These
> files match the CSV shape described in the brief and mirror the conventions of the
> `ml_csv` exporter in the sibling `qecgen` project, which was built for a consumer that
> asked for `Primary Key, Target, Variable_1, …`. No fixture supplied by the Nexus team has
> been run against them. That claim will not be made until one passes.

## The task, and what it turned into

Four sources were handed over described as "natural versus engineered" genomic datasets.
**None of those four supplies a general target of that kind.** A fifth adapter reads a narrower,
recorded-event contrast from a Supporting Information ZIP that has to be fetched by hand. Details
and evidence are in
[`reports/source_audit.md`](reports/source_audit.md); the short version:

| Source | Reality | Honest target |
|---|---|---|
| **FELIX SI signatures** | Event-relative provenance recorded by the test-and-evaluation team; not a general benchmark | **0 host-context sequence not introduced by the event / 1 introduced by the event** |
| FELIX / GUARDIAN repo | Detector *predictions*, single-class at sequence level | genetic part role |
| Addgene / GEAC | 100% engineered plasmids. No natural class exists. | lab of origin, 1,314 classes |
| SeqScreen / FunSoC | An analysis pipeline; "sequence of concern" ≠ "engineered" | 32 pathogenicity mechanisms |
| CodonTransformer | A codon-optimisation model over natural genes | natural vs model-generated |

Four of the five produce CSVs with labels supported by their sources. Only Addgene is still
blocked, on a gated download. Nothing was fabricated to fill the gap.

**The event-relative label that had to be got right.** A deletion's "engineering signature" is
DNA the recorded event removed, so it belongs in class 0: host-context sequence not introduced by
that event. Class 1 is reserved for sequence introduced by the recorded event. `insertion.site`
records are excluded because the source does not establish either sequence-level provenance class.
Treating every record in a file called "engineering signatures" as introduced would fabricate the
meaning of the target. See
[`adapters/felix_signatures.py`](src/nexus_genomics/adapters/felix_signatures.py).

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
.venv\Scripts\python -m nexus_genomics.cli convert-all --profile full
.venv\Scripts\python -m nexus_genomics.cli validate (Get-ChildItem outputs\*.ml.csv)
.venv\Scripts\python -m nexus_genomics.cli convert-all --profile demo --out-dir outputs\strict --strict-single-target
.venv\Scripts\python -m nexus_genomics.cli validate (Get-ChildItem outputs\strict\*.ml.csv) --reports-dir reports\strict
.venv\Scripts\python -m nexus_genomics.cli compliance --profile demo --out-dir outputs\strict
```

Native target layout remains the default. The six original demo/full tables retain their names and
target blocks; FELIX signatures adds one demo and one full table, for eight current native tables.
Strict layout is opt-in and writes to the separate directory shown above.

Three sources download what they need automatically — about 11 MB for the demo profile and about
166 MB for `--profile full` — and need no credentials. FELIX signatures needs the free browser
download described below. Addgene reports exactly what authorized files you must supply.

## Where each dataset comes from

| Source | How to get it | Size | Credentials |
|---|---|---|---|
| `codontransformer` | Streamed by HTTP Range from Zenodo 13262517 (6.08 GB file, never downloaded whole) | demo 5 MB · **full 160 MB** | none |
| `seqscreen` | Rice Box accession lists + UniProt REST, cached under `data\raw\` | 5 MB | none |
| `felix` | Repository tarball from GitHub | 1.3 MB | none |
| `felix_signatures` | **You must download it** — see below | 175 KB | none (browser only) |
| `addgene` | **You must download it** — see below | 647 MB | free DrivenData account |

Measured from the shipped sidecars, not estimated. The `full` CodonTransformer profile reads
40 windows of 4 MB; on a metered link that is the number to watch.

### You must supply: the FELIX Supporting Information

This is the fifth source's narrower recorded-event contrast, not general
natural-versus-engineered ground truth. It is **free** — the HTTP 403 a script gets is a
Cloudflare bot challenge, not a paywall — but it needs a browser:

1. Open <https://pubs.acs.org/doi/10.1021/acssynbio.3c00398?goto=supporting-info> and let the
   Cloudflare interstitial clear.
2. Under *Supporting Information*, click the link rendered as **sifile1 - zip file**.
3. Save it to `data\raw\felix_signatures\sb3c00398_si_001.zip` (about 175 KB).
4. `.venv\Scripts\python -m nexus_genomics.cli convert felix_signatures --profile full`

Redistribution by third parties is not established for these files (the paper grants
reproduction rights to the U.S. Government only), so keep the outputs local and cite
Adler et al., *ACS Synth. Biol.* 2024, **13**, 1105–1115.

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
adapter.load()  →  optional project  →  clean  →  encode  →  place  →  write  →  validate
```

| Stage | Module | What it does |
|---|---|---|
| Read the native format | `adapters/*.py` | One module per source; emits `SampleRecord`s |
| Project targets | `projection.py` | Losslessly selects one named target per strict table |
| Clean | `cleaning.py` | Uppercase, strip whitespace, U→T, IUPAC policy, quarantine |
| Encode | `encoding.py` | A=1…Z=26, 0 for padding |
| Place | `encoding.place` | Pad, truncate or window to the configured length |
| Write | `nexus_csv.py` | The table plus two sidecars, with hashes binding both CSVs |
| Record run | `run_receipt.py` | Atomic `convert-all` receipt binding config, layout, sources, paths, and hashes |
| Validate | `validation.py` | Integrity and semantic checks expressed as data |
| Report | `validation.render_markdown` / `cli.py` | Markdown and JSON, rendered from one object |
| Audit compliance | `compliance.py` | Email-contract evidence with four explicit statuses |

Everything after `adapter.load()` is source-agnostic, so all five registered sources obey identical
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
| `felix_natural_vs_engineered` | `target` | `0` host-context sequence not introduced by the recorded event (including flanks or sequence it removed) · `1` sequence introduced by that event |
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

The manifest records a BLAKE2b-256 hash for the complete table and the complete samples sidecar,
and validation requires unique sidecar IDs in exact table order. `convert-all` also writes
`.conversion-run.json` atomically in the selected output directory. That receipt binds a complete
run to its profile, raw config hash, target layout, enabled/converted/gated source sets, exact table
inventory, and table hashes.

### Native and strict target layouts

`output.target_layout: native` is the default. It preserves the six established tables from the
original source set: CodonTransformer and FELIX remain unsuffixed single-target tables, while
FunSoC remains one lossless 32-target table. FELIX signatures is also single-target and
unsuffixed.

`--strict-single-target` is opt-in. CodonTransformer, FELIX, and FELIX signatures each remain one
unsuffixed table. FunSoC is split losslessly into 32 deterministically named binary tables; every
table has bare `target`, and every projected manifest records the full ordered source target names,
deterministic slugs, row population, samples hash, and an ordered hash of every target column.
These outputs are conservative about the email's single-target shape, but are still **not verified
Nexus-compatible**. Projection changes layout, never scientific label semantics.

The current five-adapter registry produces 35 available strict demo tables: the original strict
set's 1 CodonTransformer, 1 FELIX, and 32 FunSoC tables, plus 1 FELIX-signatures table. Addgene
remains absent until its authorized DrivenData files are supplied.

**The manifest sidecar is what makes a file provably ours.** A bare `.ml.csv` is
shape-identical to any other one-header-row CSV. The group is staged and then committed with
one `os.replace` per file, so a failure *before* the commit leaves nothing behind. The
replaces themselves are not one transaction: a kill between them can leave a new table beside
the previous run's sidecar, which the `content_hash` check then catches.

Freshly regenerated from locally available raw inputs (2026-08-29):

| File | Rows | Cols | Size | Classes |
|---|---|---|---|---|
| `codontransformer_natural_vs_optimized_demo.ml.csv` | 992 | 1,002 | 2.2 MB | 496 / 496 |
| `codontransformer_natural_vs_optimized_full.ml.csv` | 31,730 | 1,002 | 70.6 MB | 15,865 / 15,865 |
| `felix_natural_vs_engineered_demo.ml.csv` | 200 | 1,002 | 0.4 MB | 100 / 100 |
| `felix_natural_vs_engineered_full.ml.csv` | 315 | 1,002 | 0.7 MB | 164 class 0 / 151 class 1 |
| `felix_guardian_part_role_demo.ml.csv` | 400 | 1,002 | 0.9 MB | 8 classes |
| `felix_guardian_part_role_full.ml.csv` | 481 | 1,002 | 1.0 MB | 8 classes |
| `funsoc_pathogenicity_mechanism_demo.ml.csv` | 583 | 1,033 | 1.4 MB | 32 labels, 12–70 positives each |
| `funsoc_pathogenicity_mechanism_full.ml.csv` | 13,164 | 1,033 | 29.8 MB | 32 labels, 23–4,549 positives each |

`addgene_lab_attribution_*` is not listed because it needs the gated download. All eight native
tables above come from one clean rebuild with network access to Zenodo, UniProt, GitHub and Box
all working; the run receipt is complete and every table's hash matches it.

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
.venv\Scripts\python -m pytest -q                        # 203 tests
.venv\Scripts\python -m ruff check . ; .venv\Scripts\python -m ruff format --check .
.venv\Scripts\python -m mypy src tests                   # strict
.venv\Scripts\python -m nexus_genomics.cli validate (Get-ChildItem outputs\*.ml.csv)
.venv\Scripts\python -m nexus_genomics.cli validate (Get-ChildItem outputs\strict\*.ml.csv) --reports-dir reports\strict
.venv\Scripts\python -m nexus_genomics.cli compliance --profile demo --out-dir outputs\strict
```

A successful run leaves eight native files and 35 strict demo files with zero failed validator
checks. Compliance then reports overall `INCOMPLETE` (correctly, since two requirements remain
externally blocked), zero `FAIL`, explicit Addgene and Nexus blockers, and `UNSUPPORTED_BY_SOURCE`
for the requested general target. That is the current state: `.conversion-run.json` is complete,
`compliance` exits 0, and reproducibility has been re-confirmed after the FELIX-signatures label
fix described in Known limitations.

Reproducibility — two runs of the same configuration must produce identical bytes:

```powershell
.venv\Scripts\python -m nexus_genomics.cli convert-all --profile demo --out-dir outputs\rerun
.venv\Scripts\python -m nexus_genomics.cli reproduce `
    outputs\codontransformer_natural_vs_optimized_demo.ml.csv `
    outputs\rerun\codontransformer_natural_vs_optimized_demo.ml.csv
```

Each file gets the applicable integrity and semantic checks, including several that are easy to
skip: the padding invariant (every row must be residues then padding, never a zero before a
non-zero), a decode round trip of randomly sampled rows against a per-sample digest, the full
samples-sidecar hash and exact ordered IDs, recomputation of the table `content_hash`, strict
integer targets, documented label values, and a refusal to ship an all-zero target column.
Duplicate sequences, conflicting labels on identical sequences and class balance are
*reported* rather than failed — they are facts a reader must know, not defects.

Native reports land in `reports/validation_report.md` and `.json`; strict reports land under
`reports/strict/`. Compliance writes `reports/email_compliance_report.md` and `.json`. Each pair
is rendered from the same report object so its formats cannot disagree. Compliance distinguishes
`has_failures`, `exit_ok`, and `all_requirements_satisfied`; its headline is `INCOMPLETE` whenever
an external blocker or unsupported requirement remains, even though only `FAIL` makes the command
exit nonzero.

After an interrupted or failed conversion, do not infer that older tables are current. Reuse is
allowed only when `.conversion-run.json` is complete and matches the current profile, config,
single-target layout, registry source set, exact table inventory, and hashes, and when validation
and compliance also pass without `FAIL`. Multi-file table/sidecar replacement is detectable but is
not a single filesystem transaction. The receipt cannot prove that a raw file stayed unchanged
after a completed run; rerun conversion if that boundary matters.

## Self-review

After everything was working, the project was reviewed as if it were someone else's work:
four independent reviewers by failure mode, then an adversarial verifier per finding whose
job was to refute it. 43 findings, 28 confirmed, 10 partial, 5 refuted; all confirmed ones
fixed. The most valuable one was that the round-trip check had been comparing the encoder
against itself and so could not have caught an off-by-one. See
[`reports/self_review.md`](reports/self_review.md).

## Known limitations

The requested general natural-versus-engineered target remains unsupported by the four original
sources. Strict projection cannot change that scientific fact. FELIX signatures provides the
narrower provenance contrast documented below, not a general benchmark over arbitrary sequences.

1. **The FELIX-signatures contrast is small, event-relative, and length-confounded.** Class 0 is
   host-context sequence not introduced by the recorded event; it is not a global claim that the
   sequence is natural. Class 1 is sequence introduced by that event. Padding makes the strong
   length signal plainly visible, so `length_only_baseline_accuracy` is recorded in the manifest.
   Rearrangements, point mutations, and `insertion.site` records are excluded because neither
   emitted sequence-level provenance class is established for them; malformed FASTA records are
   quarantined rather than repaired by guessing. A pre-release adversarial review found that the
   element-type parser read only the token immediately after a signature's length marker, which
   misread 6 of 335 elements whose name carries an extra descriptor there (`..._168bp_CDSpartial_
   insertion_...`, or an inline indel description before the true type) and wrongly excluded them
   as unsupported. Fixed by scanning every token after the length marker for a recognised type
   keyword; the full file grew from 309 to 315 rows (163/146 to 164/151) once those 6 were
   correctly classified. See `adapters/felix_signatures.py`'s `_element_type` and
   `tests/test_felix_signatures.py`'s two regression tests for the exact case.
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
  projection.py   generic lossless single-target projection for strict outputs
  run_receipt.py  atomic complete-run evidence for convert-all
  validation.py   integrity and semantic checks, as data rather than prose
  compliance.py   executable email-contract audit and Markdown/JSON renderers
  adapters/       one module per source, registered in ADAPTERS
tests/                           203 tests, no network required
data/raw/                        downloads; never committed
outputs/                         the CSVs and sidecars
reports/                         source audit, validation, self-review, open questions, demo checklist
```

Conventions follow the sibling `qecgen` project: exact-pinned dependencies, `ruff` at line
length 100, `mypy --strict`, flat `tests/`, and docstrings that explain **why a trap exists**
rather than what the code does.
