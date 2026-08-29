# AGENTS.md

This file provides guidance to Codex when working with code in this repository.

## What this is

`nexus-genomics` converts public genomic datasets into a supervised-learning CSV: one column
per sequence position, integer tokens, a manifest sidecar carrying the provenance. Four
adapters read four native formats; everything after that is one shared path.

**The load-bearing finding, which the original brief got wrong: none of the four sources is a
natural-versus-engineered dataset.** The evidence is in `reports/source_audit.md`. Each adapter
emits the label its source actually supports:

| Registry key | `adapter.name` (the output filename stem) | Honest target |
|---|---|---|
| `codontransformer` | `codontransformer_natural_vs_optimized` | natural NCBI CDS vs CodonTransformer-generated |
| `seqscreen` | `funsoc_pathogenicity_mechanism` | 32 FunSoC pathogenicity mechanisms, multilabel |
| `felix` | `felix_guardian_part_role` | 8 genetic part roles |
| `addgene` | `addgene_lab_attribution` | lab of origin, 1,314 classes (gated download) |

Two standing refusals govern everything below:

- **Never add a natural/engineered column to a source that lacks one.** Such a set *could* be
  constructed by pairing Addgene against RefSeq; `reports/unresolved_questions.md:112-120`
  records why it deliberately was not — the classes would differ by database, assembly pipeline
  and GC content as much as by engineering, so a model would score well by learning the source.
- **Never claim Nexus compatibility** in code, docs or commit messages until an exporter passes
  a fixture supplied by the Nexus team. The input format has never been published. This is
  enforced in the data, not just in prose: every manifest carries
  `"nexus_compatibility": "NOT VERIFIED"` (`pipeline.py:200-204`), and `cli formats` prints the
  disclaimer, asserted by `tests/test_cli.py:203-208`.

## Commands

Run from the repo root — `cli.py:37-40` uses relative defaults and `tests/test_cli.py:215`
reads `Path("config/default.yaml")` relative to the working directory.

```powershell
.venv\Scripts\python -m pip install -e ".[dev]"

.venv\Scripts\python -m pytest -q                                   # 104 tests, no network
.venv\Scripts\python -m pytest tests\test_adapters.py               # one file
.venv\Scripts\python -m pytest tests\test_adapters.py::test_the_stratified_cap_is_deterministic
.venv\Scripts\python -m pytest -k stratified_cap -v                 # by keyword

.venv\Scripts\python -m ruff check . ; .venv\Scripts\python -m ruff format --check .
.venv\Scripts\python -m mypy src tests                              # strict

.venv\Scripts\python -m nexus_genomics.cli audit
.venv\Scripts\python -m nexus_genomics.cli convert <key> --profile demo
.venv\Scripts\python -m nexus_genomics.cli convert-all --profile demo   # or full
.venv\Scripts\python -m nexus_genomics.cli validate outputs\*.ml.csv
.venv\Scripts\python -m nexus_genomics.cli reproduce <a.ml.csv> <b.ml.csv>
.venv\Scripts\python -m nexus_genomics.cli formats
```

- **The virtualenv is required, not merely tidy.** The optional `generate` extra pins
  `numpy<2.0.0` against this project's `numpy==2.3.3` (`pyproject.toml:36-46`), so installing it
  alongside the core dependencies downgrades numpy underneath everything else on the machine.
  It exists only for the optional `verify-generation` cross-check and must never become a
  silent fallback — separate environment, always.
- pip only. No `uv`, no poetry, and **no lockfiles**: reproducibility comes from exact `==` pins
  (`pyproject.toml:13-16`).
- `filterwarnings = ["error::DeprecationWarning"]` (`pyproject.toml:80`) — a deprecation warning
  fails the suite.
- A `network` marker is declared (`pyproject.toml:77-79`) but **used by zero tests**. The whole
  suite is offline and needs no gated download; keep it that way.
- This is PowerShell: chain with `;`, not `&&`. The repo is not a git repository.

## Architecture

Dependencies flow one way. Everything downstream of `SampleRecord` is source-agnostic;
everything upstream of it is source-specific and lives in `adapters/`.

```
adapters/base.py   Adapter Protocol, Availability, LoadedSource, stable_class_index, one_hot
adapters/*.py      native format -> SampleRecord. The ONLY source-specific code. It does no
                   cleaning, no encoding and no padding (base.py:5-7) -- those are shared, or
                   a defect fixed in one place survives in three others
common.py          SampleRecord, blake2b_256, canonical_json, StagedWrite, require_unique
cleaning.py        Alphabet; clean_sequence -> CleanResult | Quarantined; every op counted
encoding.py        the A=1..Z=26 cipher, EncodingMode, LengthPolicy, place()
pipeline.py        convert(): clean -> encode -> place -> manifest -> write. The one path.
nexus_csv.py       the writer, the two sidecars, and a reader that refuses foreign files
validation.py      20-21 checks expressed as data; render_markdown from the same object
cli.py             Typer; every command prints its fully resolved config before doing work
```

`cli.convert_source` (`cli.py:110-165`) loads `config/default.yaml`, hashes it as **raw text**
into `config_hash`, resolves the adapter from `ADAPTERS`, calls `adapter.load(raw_dir / key,
options)`, then `pipeline.convert`, then `nexus_csv.write_nexus_csv`.

Two naming facts, which have already caused one bug (`reports/self_review.md:59`): raw inputs
live under the **registry key** (`data/raw/felix/`), outputs under the adapter's **long name**
(`felix_guardian_part_role_demo.ml.csv`, `cli.py:136`).

Each dataset is three files sharing a stem — `.ml.csv`, `.ml.manifest.json`, `.ml.samples.csv`.
The table is `sample_id, target[…], position_0001 … position_1000`.

## Invariants that produce silently wrong data if broken

Breaking one of these yields a well-formed file containing wrong data, which passes casual
inspection and often passes the tests too.

**Labels**

- **A missing label is refused, never defaulted to 0** (`pipeline.py:107-112`,
  `common.py:48-52`). A missing label quietly becoming the negative class is the exact failure
  this project is organised against.
- **Class maps come from the whole corpus, never the capped sample** (`felix_guardian.py:126-131`,
  `addgene.py:133-138`). Otherwise `target == 2` means `origin_of_replication` in the demo file
  and `engineered_tag` in the full one, and a model trained on one is silently graded against
  the wrong classes on the other. Addgene reads the vocabulary from the *header* of
  `train_labels.csv` with `nrows=0` for this reason.
- **`stable_class_index` sorts rather than using first-seen** (`base.py:88-96`), or
  `content_hash` depends on input row order and the reproducibility check fails for a reason
  that has nothing to do with the data.
- **An all-zero one-hot row is refused** (`addgene.py:241-255`): `idxmax` returns the first
  column and says nothing about it — "a fabricated label, produced by a function that looks like
  it is reading one".

**Encoding**

- **The round-trip digest is taken from the cleaned source slice, never from
  `decode(placed.tokens)`** (`pipeline.py:138-145`, restated at `validation.py:377-383`, guarded
  by `tests/test_validation.py:149`). This is the difference between a round-trip check and a
  tautology: hashing the decoded tokens compares the encoder against itself, so an off-by-one in
  `place` would pass on every row of every file. The project's own review found this, and it is
  the single most load-bearing comment in the tree.
- **Tokens are categorical, not ordinal** (`encoding.py:5-11`). `T=20` is not ten times `C=3`.
  The warning is repeated in every manifest because it is the easiest way to misuse these files.
- **`0` is reserved for padding, and padding only ever trails.** No letter maps to 0; the
  validator refuses a zero before a non-zero (`validation.py:221-244`).
- **`stride > width` is refused** (`encoding.py:149-159`) — it would drop the residues between
  consecutive windows while every `Placed` still reported `removed=0`, making the loss invisible
  in both the manifest and the sidecar. The window loop also breaks as soon as one window
  reaches the end (`encoding.py:174-179`), or the tail of every sequence is weighted far above
  its head.
- **`codon` and `amino_acid_translated` are defined and refused, not omitted**
  (`encoding.py:61-71`, `pipeline.py:56-70`), so nobody reads the absence of a codon mode as
  "codons are the same thing as nucleotides". A protein source sets `alphabet: amino_acid` and
  still uses `letter`.

**Cleaning**

- **Refuse, never coerce** (`cleaning.py:3-8`). A sequence is either cleaned by a documented rule
  or quarantined. IUPAC codes keep their own cipher value rather than collapsing to N;
  terminators are stripped from the *end* only, because mid-sequence would splice two reading
  frames together; `U→T` is opt-in per source and raises for a protein alphabet, since U is
  selenocysteine (`cleaning.py:56-61`, `151-163`).
- **The alphabet check is one-directional and the code knows it.** A protein mislabelled DNA is
  caught — E, Q and I are not bases — but DNA mislabelled `amino_acid` cannot be, because every
  nucleotide letter is a valid residue. So
  `protein_sequences_using_only_nucleotide_letters` **counts** rather than guessing, and the
  count reaches the sidecar (`cleaning.py:96-108`).

**The writer**

- **`COLUMN_ORDER` has one definition and three consumers** — header builder, row writer, read
  offsets (`nexus_csv.py:67-75`). When these were independent hardcodings in the sibling project,
  feature values were stored under the target's column name and the round trip stayed green.
- **`position_*` is one-based and zero-padded to the file's own width** (`nexus_csv.py:110-132`).
  One-based because biology is, and a `position_0` is an off-by-one waiting to happen against a
  GenBank coordinate. Zero-padded because `sorted(df.columns)` — and every column sort hiding
  inside a join or a feature-store schema — puts `position_10` before `position_2` and silently
  permutes the feature matrix. **`target_00 … target_31` is zero-padded too**
  (`nexus_csv.py:135-153`), which qecgen does not do because it never wrote more than one
  observable; FunSoC has 32, where the unpadded spelling permutes the *label* matrix. Because the
  width is a property of each file, the sidecar publishes the literal names — read
  `columns.feature_columns`, never reconstruct it.
- **The manifest sidecar is what makes a file provably ours** (`nexus_csv.py:303-322`): a bare
  `.ml.csv` is reported *not one of ours* rather than corrupt. That claim is sound only because
  `StagedWrite` commits the group in one pass (`common.py:115-154`). The commit loop is not a
  single transaction — a kill between two replaces can leave a new table beside an old sidecar —
  but that state is detectable, because `content_hash` will not match.
- **`content_hash` covers the table only** (`nexus_csv.py:257-259`). `generated_at` changes every
  run, so folding the manifest in would make the reproducibility check vacuous. Digests are
  BLAKE2b-256 and say so in `content_hash_algorithm` (`common.py:35-41`) — the sibling project
  once shipped a field named `content_sha256` holding a BLAKE2b digest.
- **An empty table is refused** (`nexus_csv.py:222-226`) — it is almost always a silently failed
  adapter and would pass most validation. A short row is refused too, never padded
  (`nexus_csv.py:233-243`), because padding it there would fabricate residues.

**The validator**

- **It must never crash on the corruption it exists to detect.** A non-integer feature column, an
  out-of-range token and a failed padding invariant each bail out with an explicit *failing*
  `remaining_checks_skipped` (`validation.py:183-193`, `209-219`, `245-256`). The last matters
  most: `decode()` skips a pad wherever it sits, so on an interleaved row it returns a
  wrong-but-plausible sequence, and continuing would print confirmatory passes underneath a
  failure that invalidates all three.
- **Reported and failed are different things** (`validation.py:6-9`). Duplicate sequences,
  conflicting labels and class balance are recorded with `passed=True` — they are facts a reader
  must know, not defects. Only an established defect fails.
- **The vectorised paths are load-bearing, not cosmetic.** `_decode_matrix` replaced 31.7 M dict
  lookups (`validation.py:74-80`); hoisting `label_matrix` out of the row loop took the FunSoC
  file from 87 s to 0.02 s (`validation.py:282-285`).

**Exit codes**

- **`convert-all` exits 0 for a missing gated file and 1 for anything else** (`cli.py:185-207`).
  Otherwise a scheduled run that ran out of memory looks identical to one waiting on a download,
  and the next step consumes stale outputs believing they are fresh.

## Per-adapter traps

Each module's docstring is the full statement; these are the ones that bite from outside.

- **codontransformer** — model-artefact detection, not a benchmark. Rows strictly alternate
  natural/generated, so row *position* predicts the target perfectly: a bare `reset_index()`,
  `KFold(shuffle=False)` or `df.head(odd_n)` all leak or skew. `sample_id` ends `__natural` or
  `__finetune`, so **the feature block is the `position_*` columns and nothing else**
  (`codontransformer.py:195-211`). Primary keys derive from row content, not position, because
  the demo reads 8 Range windows and the full reads 40 — an enumeration index would name a
  different gene in each, and `require_unique` only checks within one file
  (`codontransformer.py:134-147`). Do not join Zenodo 13262517 with 12509224; they are not
  row-aligned.
- **seqscreen / FunSoC** — a FunSoC is a pathogenicity mechanism, usually in a *natural* toxin;
  `vbeo == 's'` is a decoy that flags cloning vectors. The Box listing **paginates at 20**, so
  fetching page 1 silently loses 12 FunSoCs, and Box answers a throttled download with HTTP 200
  and an HTML challenge page that would become a FunSoC's accession list
  (`seqscreen.py:306-320`). Batch-splitting is limited to `{400, 414}` (`seqscreen.py:64-71`):
  halving on a 429 turns one rate limit into a ban. The cap walks the FunSoCs round-robin
  because alphabetical-first-N left 6 target columns all-zero (`seqscreen.py:217-228`); a capped
  file is therefore **not** distribution-representative — weight loss with
  `positives_per_funsoc_in_this_file`, never the source-list counts.
- **felix / GUARDIAN** — the adapter measured the corpus and refuses the binary: only 280 of
  1,670 sequences link to a detector call, and all 280 are `engineered`, so a binary file would
  contain exactly one class. Two class names contain the word "engineered"
  (`engineered_region`, `engineered_tag`) but they are **part types**, and
  `target.isin([...])` on them is exactly the fabrication this adapter exists to avoid
  (`felix_guardian.py:186-195`). Two coordinate conventions coexist in the sidecar:
  `source_start`/`source_end` are 0-based half-open in GenBank accession coordinates, while
  `window_start`/`window_end` are 1-based inclusive in the feature's own sequence
  (`felix_guardian.py:203-210`). Conflicting detector calls are recorded as `"CONFLICT"`, not
  resolved — two detectors disagreeing is a fact about the corpus, not a tie to break.
- **addgene** — 100% engineered plasmids. `I7FXTVDP` is the pooled "Unknown Engineered" class,
  meaning unknown *laboratory*, still engineered; reading it as "unknown or natural" inverts it.
  **The DrivenData competition rules prohibit redistributing the data or anything derived from
  it** — the raw files and the output CSV stay on this machine, `.gitignore` excludes
  `outputs/`, and the unofficial Kaggle mirror confers no rights and is deliberately not used.

## Extension points

**A new source is four things, not one.** The module and the registry line are the advertised
pair; the test suite enforces two more.

1. One module in `adapters/` with class attributes `name`, `source_url`, `licence`,
   `target_description`, plus `availability()` and `load()`. It satisfies the `Adapter` Protocol
   structurally — nothing is subclassed. The registry contract test
   (`tests/test_adapters.py:28-34`) asserts `source_url` starts with `http`,
   `len(licence) > 20` and `len(target_description) > 40`.
2. One line in `ADAPTERS` (`adapters/__init__.py:21-26`).
3. A `sources.<key>` block in `config/default.yaml` with `enabled: true` and a `demo:` block.
   `tests/test_cli.py:211-221` asserts `set(parsed["sources"]) == set(ADAPTERS)`, so an adapter
   without a config entry **breaks the suite** — which is the point, since `convert-all` would
   otherwise skip it in silence.
4. A banner-comment section in `tests/test_adapters.py` with a `_write_<source>(tmp_path)`
   synthetic-fixture writer and trap-named tests. If the source is gated, add an
   `availability()` test asserting the `manual_steps` text: `FileNotFoundError` alone sends
   someone hunting for a bug that is not there.

`load()` returns a `LoadedSource` whose `label_semantics` names every value that may legally
appear in a target column and what it means — validation checks the table against it, so a
source whose labels are not ground truth says so there, in text a reader will actually see. Set
the right `Alphabet`, record `raw_input_hashes`, and use `stable_class_index` / `one_hot` from
`base.py` rather than re-deriving a class map. `Availability.manual_steps` is not optional
politeness: it is the difference between "this pipeline is broken" and "this pipeline needs one
file you fetch in a browser".

**A new validation check** is one `result.add(...)` in `validate_file`. Decide deliberately
whether it *fails* or merely *reports*, and place it after the bail-outs it depends on.

**A new config key** belongs in `config/default.yaml`, not a CLI default: that file's text digest
is `config_hash` in every output, so a setting that changes bytes must be inside it.

## Conventions

The fuller statement is `C:\Projects\qecgen\AGENTS.md`, which this project follows.

- Dependencies exact-pinned with `==`. Python ≥3.13. Typer CLI. Flat `tests/` with **no
  `conftest.py`, no `parametrize`, and no shared fixture module** — synthetic fixtures are built
  by module-level `_write_*` helpers into `tmp_path`.
- ruff `E,F,I,N,UP,B,A,C4,SIM,RUF` at line length 100, with `E501` ignored under `tests/`.
  `mypy --strict` clean over `src` and `tests`.
- **Docstrings explain why a trap exists**, not what the code does. That is the house style
  throughout and the reason these invariants survive refactoring — preserve the explanations when
  editing, and if you correct a claim, correct it in the docstring, `README.md` and this file
  together.
- Test names are full descriptive sentences naming the trap, and each docstring says why deleting
  the test would be dangerous.
- **Refuse, never coerce.** Error text states what the coercion would have fabricated.
- `reports/source_audit.md`, `reports/unresolved_questions.md` and
  `reports/monday_demo_checklist.md` are hand-written and kept. `reports/validation_report.md`,
  `reports/validation_report.json` and `reports/source_audit_status.md` are generated by the CLI,
  gitignored, and must not be hand-edited.
- **`AGENTS.md` is a deliberate copy of this file** — Codex reads that name, Claude Code reads
  this one, and both are kept self-contained rather than one pointing at the other. Every edit
  here must be made there too. The only lines that may differ are the title and the one sentence
  naming the tool; `diff CLAUDE.md AGENTS.md` should report exactly those two hunks.

## Known gaps

- **Blocked on the Nexus team:** the real column and row ceiling (1,000 positions was chosen as a
  safe starting point, not established), and whether a 32-column multilabel target block is
  acceptable at all — the brief says "column", singular. If it is not, FunSoC must become 32
  binary files or one collapsed column, and both lose information.
- **Blocked on a download only the user can perform:** the GEAC competition files (free
  DrivenData account plus rules acceptance), and `sb3c00398_si_001.zip`, which is the **only**
  genuine sequence-level natural-versus-engineered ground truth found across all four sources.
  Its HTTP 403 is a Cloudflare challenge, not a paywall, so a browser gets it and a script does
  not. Even once it lands, what the honest *negative* class should be is an open question that
  `reports/unresolved_questions.md:69` marks as one not to decide unilaterally.
- **Recorded, non-blocking:** `fine_tuned` in the CodonTransformer corpus is undocumented and
  must never be used as a target; the strided sample from the 6.08 GB Zenodo file is reproducible
  but not uniform random, because rows are grouped by organism.
- **Pending in code:** `read_nexus_csv` (`nexus_csv.py:325`) still materialises the whole table as
  strings, which is why the streaming `read_sidecar_and_shape` exists beside it for the validator.
