# Data dictionary

What every column in every output file means, and what value each encoded number represents.
This is a human-readable index over information that already exists, per file, in that file's
own `.ml.manifest.json` — if this document and a manifest ever disagree, **the manifest is
correct**; it is generated fresh by the code every run, this file is hand-maintained.

## Every file has the same three-column shape

```
sample_id, target(s), position_0001 ... position_NNNN
```

| Column | Type | Meaning |
|---|---|---|
| `sample_id` | string | Primary key. Unique within the file. For windowed sequences, a suffix like `__w00` marks which window of a longer parent sequence this row is. |
| `position_0001` … `position_NNNN` | integer, 1–26 or 0 | One DNA/protein letter per column, ciphered `A=1, B=2, C=3 … Z=26`. `0` means "no residue here" (padding) — it is never a valid letter. Column count (`NNNN`) is fixed per file and recorded in the manifest; currently 1,000 everywhere. |
| `target` (or `target_00` … `target_31` for FunSoC) | integer | What's actually being predicted. **Meaning is different for every dataset** — see below. Never treat these as the same variable across files. |

**The encoded numbers are categories, not quantities.** `T=20` is not "ten times" `C=3`, and `G=7`
does not sit "between" C and T. A model must treat `position_*` columns as categorical (one-hot
or an embedding layer), never as a magnitude. This warning is repeated in every manifest for
exactly this reason — it's the easiest way to misuse these files.

**Whether the letters are DNA or protein depends on the file** (recorded as `alphabet` in each
manifest):

| Dataset | Alphabet |
|---|---|
| `codontransformer_natural_vs_optimized` | `dna_nucleotide` |
| `felix_natural_vs_engineered` | `dna_nucleotide` |
| `felix_guardian_part_role` | `dna_nucleotide` |
| `funsoc_pathogenicity_mechanism` | `amino_acid` — this one is protein, not DNA |
| `addgene_lab_attribution` | `dna_nucleotide` |

## What `target` actually means, per dataset

### `codontransformer_natural_vs_optimized` — `target`

| Value | Meaning |
|---|---|
| `0` | Natural coding sequence (NCBI RefSeq CDS) |
| `1` | The same protein, re-coded by the CodonTransformer model (synonymous codon substitution only — identical protein, different DNA spelling) |

### `felix_natural_vs_engineered` — `target` *(the FELIX Supporting Information dataset — the one with real natural-vs-engineered ground truth)*

| Value | Meaning |
|---|---|
| `0` | Host-context sequence **not introduced** by the recorded engineering test event — either a flank beside the change, or a region the event *removed* (a deletion's signature is the DNA that was taken out, so it counts as pre-existing/natural, not engineered) |
| `1` | Sequence **introduced** by the recorded event (an insertion or a plasmid-element cargo sequence) |

This is event-relative, not a universal natural/engineered claim — see `README.md` and
`GUIDE.md` for the full caveat (small sample, length-confounded, and a specific labelling trap
this project caught). Sourced from the Supporting Information of
[Adler et al., *ACS Synthetic Biology* 2024, 13, 1105–1115](https://doi.org/10.1021/acssynbio.3c00398).

### `felix_guardian_part_role` — `target` *(the public FELIX/GUARDIAN/MIDOE repository — a different dataset from the one above)*

| Value | Part role |
|---|---|
| `0` | CDS |
| `1` | engineered_region |
| `2` | engineered_tag |
| `3` | five_prime_UTR |
| `4` | origin_of_replication |
| `5` | promoter |
| `6` | replicon |
| `7` | terminator |

This is **what kind of genetic part the sequence is** (a promoter, a coding region, etc.), not
where it came from. `engineered_region`/`engineered_tag` are part *type* names, not an
engineered-vs-natural signal — see the caveat in `README.md`. These are GUARDIAN pipeline
annotations, not experimentally verified ground truth.

### `funsoc_pathogenicity_mechanism` — `target_00` … `target_31` (32 columns, multilabel)

Each column is an independent 0/1 flag for one specific pathogenicity mechanism (FunSoC = Function
of Sequence of Concern); a protein sequence can carry several. `0` means "not curated as a positive
for this FunSoC" — **not** a verified negative. This is a disease/pathogenicity-mechanism signal,
not an engineered-vs-natural label; most flagged sequences are natural toxins, adhesins, or
virulence factors.

| Column | FunSoC |
|---|---|
| `target_00` | AntibioticResistance |
| `target_01` | AvirulencePlant |
| `target_02` | BacterialCounterSignaling |
| `target_03` | CounterImmunoglobin |
| `target_04` | Cytotoxicity |
| `target_05` | DegradeECM |
| `target_06` | DevelopmentInHost |
| `target_07` | DisableOrgan |
| `target_08` | HostCellCycle |
| `target_09` | HostCellDeath |
| `target_10` | HostCytoskeleton |
| `target_11` | HostGTPase |
| `target_12` | HostTranscription |
| `target_13` | HostTranslation |
| `target_14` | HostUbiquitin |
| `target_15` | HostXenophagy |
| `target_16` | InduceInflammation |
| `target_17` | InvasionHostCell_Viral |
| `target_18` | NonviralAdhesion |
| `target_19` | NonviralInvasion |
| `target_20` | PlantRNASilencing_Viral |
| `target_21` | ResistHostComplement |
| `target_22` | ResistOxidative |
| `target_23` | SecretedEffector |
| `target_24` | Secretion |
| `target_25` | SuppressDetection |
| `target_26` | ToxinSynthase |
| `target_27` | ViralAdhesion |
| `target_28` | ViralCounterSignaling |
| `target_29` | ViralMovement |
| `target_30` | VirulenceActivity |
| `target_31` | VirulenceRegulator |

In `--strict-single-target` mode, each of these becomes its own file with a bare `target` column;
the manifest's `target_projection` block records which of the 32 it holds.

### `addgene_lab_attribution` — `target` *(not yet built — needs the gated GEAC files)*

| Value | Meaning |
|---|---|
| one integer per lab | Which laboratory deposited this plasmid (up to 1,314 classes) |
| `I7FXTVDP` | The pooled **"Unknown Engineered"** class — unknown *laboratory*, still engineered. Every sequence in this dataset is engineered; there is no natural class. |

The full 1,314-entry integer-to-lab-ID map is read from the (currently missing) `train_labels.csv`
header at conversion time, so it can't be published here in advance — it will appear in
`addgene_lab_attribution_*.ml.manifest.json` the moment the gated files are supplied and the
dataset is built.

## The `.ml.samples.csv` sidecar — per-row provenance, kept out of the model table

Every dataset's samples file starts with these columns, generated by the shared pipeline
(so their meaning never varies by source):

| Column | Meaning |
|---|---|
| `sample_id` | Matches the `.ml.csv` row |
| `parent_sample_id` | The un-windowed sequence this row came from (same as `sample_id` unless the sequence was split into windows) |
| `source` | Which dataset this row belongs to |
| `original_length` | Length of the cleaned sequence before padding/truncation |
| `window_start`, `window_end` | 1-based, inclusive — where in the original sequence this row's residues sit |
| `residues_in_row` | How many real (non-padding) residues this row holds |
| `padding_cells` | How many `position_*` columns are padding (`0`) in this row |
| `residues_removed_by_truncation` | How much of the sequence was cut off to fit the fixed width |
| `residues_dropped_by_window_cap` | Residues lost because a sequence produced more windows than the configured cap |
| `source_slice_blake2b` | A hash of the exact residues encoded in this row — what the validator's decode round-trip checks against |

Then each source appends its own extra columns:

| Dataset | Extra sample columns |
|---|---|
| CodonTransformer | `pair_id`, `member` (natural/finetune), `organism`, `gene_id`, `protein_length`, `description` |
| FELIX natural-vs-engineered | `if_element`, `element_type`, `is_flank`, `label_basis`, `element_origin_organism`, `n_samples_containing`, `samples_containing` |
| FELIX part-role | `feature_name`, `part_role`, `genbank_source`, `source_start`, `source_end`, `source_orientation`, `guardian_call_not_ground_truth` |
| FunSoC | `uniprot_accession`, `funsocs`, `n_funsocs` |
| Addgene | `lab_id`, `is_pooled_unknown_lab`, plus ~39 phenotype columns (`bacterial_resistance_*`, `copy_number_*`, `growth_*`, …) |

## The authoritative, machine-readable version

Every `.ml.manifest.json` carries the same information this document summarises, always in sync
with the file it ships beside, under:

```json
{
  "columns": {
    "index_columns": ["sample_id"],
    "target_columns": ["target"],
    "feature_columns": ["position_0001", "..."],
    "token_values": {"letters": "A=1 .. Z=26", "padding": 0, "note": "..."}
  },
  "manifest": {
    "alphabet": "dna_nucleotide",
    "label_semantics": {"0": "...", "1": "..."}
  }
}
```

If you're building anything automated against these files, read `columns.feature_columns` and
`manifest.label_semantics` from the sidecar rather than hard-coding column names or class counts —
the manifest is generated fresh from the real data every run and is the one place these can never
drift out of date.
