# Source audit

What each of the four original sources **actually is**, plus the narrower fifth source added from
FELIX Supporting Information.

Every claim below was produced by one researcher agent working from primary sources and then
re-checked by an independent adversarial verifier whose job was to refute it. Where the two
disagreed, the verifier's finding is recorded — several of the researcher's claims were wrong
and are corrected here.

## Headline

**None of the four original sources supplies a general natural-versus-engineered benchmark.**
The fifth source, ACS Supporting Information for the GUARDIAN paper, records the narrower
event-relative distinction between sequence introduced by a recorded engineering event and
host-context sequence not introduced by that event. That supports a useful contrast, but not a
general benchmark over arbitrary sequences. The original sources retain only the labels they
actually support. No labels have been fabricated to fill any gap.

| Source | What it actually is | Honest target | Status |
|---|---|---|---|
| **FELIX Supporting Information** | Per-signature, event-relative provenance recorded by the test-and-evaluation team | **`0` host context not introduced by the event / `1` introduced by the event** | **Converted** |
| Addgene / GEAC | A real labelled corpus — of **lab of origin**, 1,314 classes. 100% engineered; no natural class exists. | `lab_id` | **Blocked**: gated download |
| IARPA FELIX / GUARDIAN | Detector *output* plus 1.27 Mbp of real DNA. Its `engineered` field is a prediction, and at sequence level it is single-class. | `part_role` (annotation) | Converted, narrowly |
| SeqScreen / FunSoC | An analysis **pipeline**. No provenance labels at all. But a real biocurated pathogenicity-mechanism label set exists. | 32 FunSoC binaries | Converted |
| CodonTransformer | A codon-optimisation **model** plus an open corpus of natural CDS, and a second record with matched model output. | natural vs generated | Converted |

---

## 1. Addgene / Genetic Engineering Attribution Challenge

| Field | Value |
|---|---|
| Official URL | <https://www.drivendata.org/competitions/63/genetic-engineering-attribution/> |
| Papers | Nielsen & Voigt 2018 ([10.1038/s41467-018-05378-z](https://doi.org/10.1038/s41467-018-05378-z)); Alley et al. 2020 ([10.1038/s41467-020-19612-0](https://doi.org/10.1038/s41467-020-19612-0)) |
| Raw file location | `data/raw/addgene/train_values.csv`, `train_labels.csv` — **you must supply these** |
| File format | CSV. `train_values.csv` 63,017 × 41; `train_labels.csv` 63,017 × 1,315 |
| Unit of one sample | One engineered plasmid |
| Sequence field | `sequence` — DNA nucleotide string |
| Ground-truth field | `train_labels.csv`, **one-hot** across 1,314 lab-ID columns, exactly one `1.0` per row |
| Label semantics | 8-character obfuscated lab identifiers. `I7FXTVDP` is the pooled **"Unknown Engineered"** class: 8,286 training rows (13.15%) from 2,438 labs with fewer than ten deposits |
| Licence | DrivenData competition rules. **Redistribution of the data or anything derived from it is prohibited.** |
| Scientifically valid to convert? | Yes — **as multiclass lab attribution only** |
| Blockers | Free DrivenData account **and** acceptance of the competition rules |

**The trap.** Every sequence here is an engineered Addgene plasmid. There is no natural class
and there never was one. "Unknown Engineered" means *unknown laboratory*, not *unknown
provenance* — reading it the other way is how a fabricated natural class would enter the
dataset. `tests/test_adapters.py::test_a_lab_attribution_file_never_contains_a_natural_class`
exists to catch exactly that.

Neither paper released its corpus (2018: "available from the corresponding author upon
request"; 2020: "Then request the data"). An unofficial Kaggle mirror is live; it confers no
rights, so this project does not use it.

One more property, **taken from the competition's own data description rather than measured
here** (the gated files have never been present on this machine): the official training file
is reported to contain 374 duplicated sequence strings over 63,017 rows. The validator
reports duplicate sequences on every file it sees, so this will be confirmed or corrected the
moment the download arrives.

## 1b. FELIX Supporting Information — a recorded-event sequence contrast

| Field | Value |
|---|---|
| Official URL | [10.1021/acssynbio.3c00398](https://doi.org/10.1021/acssynbio.3c00398), *Supporting Information* |
| Raw file location | `data/raw/felix_signatures/sb3c00398_si_001.zip` — **you must fetch this in a browser** (175 KB, free) |
| File format | ZIP of one FASTA, one CSV and one XLSX |
| Unit of one sample | One engineering signature (or one of its flanks) |
| Sequence field | The FASTA — DNA nucleotide, 335 records, 562,114 bp |
| Ground-truth field | The element **type**, encoded in each signature's name |
| Label semantics | `0` host-context sequence not introduced by the recorded event · `1` sequence introduced by that event |
| Licence | ACS SI, free to access. Reproduction rights granted to the U.S. Government only, so third-party redistribution is **not** established. Keep outputs local. |
| Scientifically valid to convert? | **Yes, as this narrower event-relative contrast; not as a general benchmark** |
| Blockers | None, once the ZIP is in place |

**The three files, whose schemas nobody had seen before this project read them:**

- `..._signature_sequences.fa` — 335 records. Names carry a lab prefix, an `IF#####` number,
  the element, **its length**, **its type** and the origin organism, e.g.
  `LBNL_IF00015_AmpR_861bp_plasmid.element_Staphylococcus.aureus`. `_flank1` / `_flank2`
  records are the sequence either side of the change.
- `..._sample_signatures.csv` — `Sample_ID,Master_Sample_ID,IF_Element`, 1,004 rows over 75
  samples and 340 elements.
- `..._sample_metadata.xlsx` — 100 samples with `Engineering_Key` (**yes 75 / no 25**),
  `Insert_Over_10_bp` (**yes 60 / no 40** — exactly the paper's Table 1), `Organism`,
  `BioSample`, `SRA`, `Sample_Type`, `Fraction_Engineered`. Sample-level, so it cannot label a
  sequence; recorded in the manifest as context.

**The label distinction that decides everything.** A signature file name does not establish
that every listed sequence was introduced by the recorded event:

- **insertion / plasmid.element / plasmid** → sequence introduced by the recorded event.
  **Class 1.**
- **deletion / compound.deletion** → host-context sequence the recorded event removed, so it
  was not introduced by that event. **Class 0.** This is event-relative provenance, not a
  global assertion that the sequence is natural.
- **flank** → host-context sequence bordering the recorded change and not introduced by it.
  **Class 0.**
- **insertion.site** → excluded: the source does not establish that the supplied sequence is
  the inserted material rather than host context at the site.
- **inversion, reassortment, translocation, point mutation** → excluded: position or base-level
  change does not establish either emitted sequence-level provenance class.

The emitted records come from the same samples, laboratories, and sequencing run. That reduces
the obvious source-batch confounder of pairing them against unrelated reference genomes, but it
does not turn this recorded-event contrast into a general natural-versus-engineered benchmark.
Fresh full output contains **315 rows: 164 class 0 and 151 class 1**. The adapter excludes 20
unsupported records, including all 8 `insertion.site` records, and quarantines 2 malformed FASTA
records.

**A parsing defect found by a pre-email adversarial review, not a defect in the source.** The
element-type parser originally read only the token immediately after a signature's length marker.
Six elements carry an extra descriptor in that slot before the real type keyword — a partial-CDS
marker or an inline indel description, e.g. `..._9bp_607-615delACGACCTCC_..._deletion_...` — so
the parser returned the descriptor, found no recognised type, and wrongly excluded a genuinely
classifiable record as unsupported (309 rows, 26 excluded, instead of the correct 315 rows, 20
excluded). It never assigned the wrong class to a record; the failure was exclusion, not
mislabelling. Fixed by scanning every token after the length marker for a recognised keyword.

**Two defects in the published file, both surfaced rather than absorbed:**

1. **Two FASTA records are missing the `>` on their header line.** A naive parser appends the
   header text and the following sequence to the *previous* record — which is how one 16 bp
   signature reads as 202,036 bp. Their names do not appear in the sample-signature CSV, so
   the repair could not be verified against anything and both were **quarantined, not
   guessed at**.
2. **11 of 222 non-flank records disagree with the length their own name states** (e.g.
   `adh1_24bp_deletion` is 1,023 bp). Every signature name declares its own length, which
   makes the file self-checking; the disagreements are recorded in the manifest as
   `length_token_disagreements`.

**The confounder you must quote alongside any result.** Introduced elements are much longer
than much of the host-context class. The encoder pads to a fixed width, so length is visible as
the point where padding starts. The freshly measured full-output
`length_only_baseline_accuracy` is **0.7587**, recorded in the manifest for exactly this reason.

## 2. IARPA FELIX / GUARDIAN / MIDOE (the public repository)

| Field | Value |
|---|---|
| Official URL | <https://github.com/raytheonbbn/midoe>; paper [10.1021/acssynbio.3c00398](https://doi.org/10.1021/acssynbio.3c00398) |
| Raw file location | `data/raw/felix/midoe-main.tar.gz` — downloaded automatically (1.2 MB) |
| File format | 340 JSON evidence files against a published JSON Schema |
| Unit of one sample | One detected sequence **feature** |
| Sequence field | `features[].sequence` — DNA nucleotide, alphabet `{A,C,G,N,T}` plus lowercase |
| Ground-truth field | **None that is ground truth.** `detections[].engineered` is a *detector assertion* |
| Licence | Apache-2.0. Redistribution permitted with attribution |
| Scientifically valid to convert? | **Not as natural-vs-engineered.** Yes as part-role annotation |
| Blockers | The repository itself has no general sequence-provenance ground truth; the separate SI source supports only the narrower recorded-event contrast above |

**Measured, by parsing all 340 files rather than reading the docs:** 1,672 features, of which
**1,670 carry real DNA** — 1,266,189 bp, 11 to 52,182 bp, median 288. 627 detections across
exactly 100 blinded `Y###` samples, matching the paper's 100 test-and-evaluation samples.

**Why no binary file was produced.** The MIDOE JSON Schema does define precisely the enum the
brief predicted:

```json
"enum": ["http://guardian.bbn.technology#natural",
         "http://guardian.bbn.technology#engineered",
         "http://guardian.bbn.technology#indeterminate_origin"]
```

described as *"whether this **evidence** is natural, engineered, or of indeterminate origin"*.
Three things disqualify it as a target:

1. It is **detector output**, not ground truth. Each detection names an `agent` — there are
   seven, and they disagree — and most carry a `confidence`. Training on it would distil
   GUARDIAN's 2022 calls including its documented false negatives.
2. At *sample* level it looks balanced (natural 257 / engineered 250). At *sequence* level it
   collapses: only **280 of the 1,670 sequences** can be linked to a call at all, and **all
   280 are `engineered`**. A binary file built here would contain one class. A detector that
   finds nothing engineered emits no feature, so there is nothing for a negative to attach to.
3. `indeterminate_origin` **never appears** in the committed data, so the
   exclude-and-audit configuration the brief asked for has nothing to act on. The
   configuration hook exists anyway; it is simply never exercised by this corpus.

**Two alternatives considered and rejected**, both documented in the output sidecar:

- *NCBI BioProject PRJNA607328* carries a per-BioSample `biol_stat` (179 engineered / 52
  natural / **59 unresolved**, and the source vocabulary contains the misspelling
  `enginered mutant`). Rejected because the unit is a multi-megabase genome that is >99.9%
  natural sequence: windowing it would stamp "engineered" onto overwhelmingly natural windows.
- *Pairing the 280 engineered features against sampled natural reference genomes.* Rejected
  because the classes would then differ by origin, assembly pipeline and length distribution
  as much as by engineering, so a model would learn the batch effect.

**What was produced instead:** `felix_guardian_part_role`, an 8-class annotation of what kind
of genetic part each sequence is. Those annotations are pipeline output too, but *"this is a
promoter"* is a claim the annotation actually makes, whereas *"this organism was engineered"*
is not one any per-sequence field here supports.

## 3. SeqScreen / FunSoC

| Field | Value |
|---|---|
| Official URL | <https://gitlab.com/treangenlab/seqscreen>; labels at <https://rice.box.com/s/60380nrh1mqfib5yic8x4damycbje838> |
| Raw file location | `data/raw/seqscreen/funsoc/*.txt` + `uniprot_sequences.tsv` — downloaded automatically |
| File format | 32 newline-delimited UniProt accession lists; sequences from the UniProt REST API |
| Unit of one sample | One UniProt protein |
| Sequence field | Fetched from UniProt — **amino acid**, not DNA |
| Ground-truth field | Membership of each of the 32 curated FunSoC lists |
| Label semantics | `1` = curated positive for that FunSoC; `0` = **not curated positive**, which is *not* a verified negative |
| Licence | SeqScreen code GPLv3; UniProt CC-BY-4.0; no explicit licence on the accession lists |
| Scientifically valid to convert? | Yes — as a **pathogenicity-mechanism** multilabel set |
| Blockers | None |

**SeqScreen is a pipeline, not a dataset.** It needs a ~169 GiB database, 32–256 GB RAM and
Linux, and it *produces* annotations rather than shipping labelled data.

**It has no natural-versus-engineered ground truth of any kind.** "Synthetic construct",
"wild-type" and "provenance" appear **zero times** in the Genome Biology paper.

**"Sequence of concern" does not mean "engineered".** A FunSoC is a pathogenicity mechanism;
the overwhelming majority of sequences carrying one are entirely natural toxins, adhesins and
virulence factors. The one field that looks like an engineering signal — output column
`vbeo == 's'` — is a decoy: it fires when a read's taxid falls under NCBI "synthetic
construct" (taxid 32630), which flags cloning vectors and artefact GenBank records. An
engineered pathogen gene hits its natural homolog and is never flagged.

**Two operational traps, both handled in code.** The Box listing **paginates at 20**; fetching
only page 1 silently loses 12 of the 32 FunSoCs, so the adapter reads both pages and refuses
to build a matrix from fewer than 32 lists. And UniProt returns a bare `HTTP 400` once the
accession query grows too long — 100 per request worked, 200 failed on the full run — so the
adapter halves the batch and retries.

**Could SeqScreen outputs become additional Nexus variables?** Yes in principle: its
per-query outputs (taxid, UniProt hit, GO terms, `bsat_hit`, `vfdb_hit`, the FunSoC columns)
are numeric or categorical and could join to another dataset on a query ID. But producing
them means running the pipeline, which needs the 169 GiB database and Linux. Not attempted.

## 4. CodonTransformer

| Field | Value |
|---|---|
| Official URL | <https://zenodo.org/records/13262517> (pairs); <https://zenodo.org/records/12509224> (training corpus) |
| Raw file location | Streamed by HTTP Range; no local file required |
| File format | One CSV, 6,083,009,710 bytes, 622,915 rows |
| Unit of one sample | One coding sequence (each source row yields a matched pair) |
| Sequence field | `natural_dna`, `base_dna`, `finetune_dna` — DNA nucleotide |
| Ground-truth field | Which column the sequence came from — provenance is exact by construction |
| Label semantics | `0` = natural NCBI RefSeq CDS; `1` = CodonTransformer-generated sequence for the **same protein in the same organism** |
| Licence | CC-BY-4.0. Redistribution permitted with attribution |
| Scientifically valid to convert? | Yes, **as model-artefact detection**, clearly named as such |
| Blockers | None |

**The training corpus is natural, as the brief expected.** Zenodo 12509224 is 1,001,197 NCBI
CDS rows; its header, read directly, is
`,dna,protein,organism,GeneID,description,tokenized,fine_tuned,len`, and sampling the
`description` field at four widely separated byte offsets returned only RefSeq accessions.

**The matched pairs already exist and did not need generating.** Zenodo 13262517 publishes
`natural_dna`, `base_dna` and `finetune_dna` in the *same row*. So no model is run, no 358 MB
of weights are downloaded, and a real hazard is sidestepped: `predict_dna_sequence` calls
`tokenize` **without** passing `max_len`, making the 2048-token ceiling unreachable through
the public API, so bulk generation would silently truncate long proteins.

**What this dataset is not.** Both members of a pair encode the identical protein, so the only
available signal is synonymous codon choice. The positive class is one generator's output,
restricted to the 17 organisms the model was fine-tuned on. It will not transfer to GeneArt,
IDT, JCat, DNAWorks or Twist output and has no bearing on biosecurity screening.

**Hazards respected by the adapter:**

- The two Zenodo records are **not row-aligned**: 13262517 has 622,915 rows, while 12509224
  has 1,001,197 rows in total of which only 470,375 belong to the 17 fine-tuned organisms.
  So pairs are taken from within a single row of 13262517 and never joined across records.
- Rows are **per-transcript, not per-gene** — one `GeneID` recurs with near-identical isoforms.
  `gene_id` is carried into the sidecar and named as the grouping key for any split.
- `fine_tuned` in the training corpus is an **undocumented** boolean and is never used as a
  target.
- The two members of a pair are asserted to be equal length, so truncation is symmetric and
  sequence length carries no class signal.

---

## Sources deliberately not used

- **The Kaggle mirror of the GEAC data.** Live, but almost certainly an unauthorised
  re-upload; it confers no redistribution rights.
- **Zenodo 13130104** (CodonTransformer ncRNA, ~45,000 human sequences). `access_right` is
  `restricted`; the files API returns `HTTP 403 Permission denied`. Access must be requested
  from the depositor.
- **The Addgene API.** Requires a token and a per-scope data licence, and even then
  "Principal Investigator/Head of Laboratory names are not included", so lab-of-origin labels
  cannot be rebuilt from it.
- **IARPA's FELIX programme page** advertises "Extensive training data sets for machine
  learning algorithms" under Accomplishments with no link, no repository and no request
  procedure. There is no public FELIX benchmark.
