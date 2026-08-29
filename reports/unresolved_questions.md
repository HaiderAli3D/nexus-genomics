# Unresolved questions

Things that are genuinely unknown, as distinct from things that were decided. Each says who
can answer it and what changes when they do.

## Answered during planning

These were asked and settled before implementation, and are recorded here so the reasoning
survives.

| Question | Answer | Consequence |
|---|---|---|
| How many position columns can Nexus ingest? | **1,000** to start, "scale up in the future if needed" | `window.length` is a config key, not a constant |
| What should the target column be called? | **`target`** | Matches `qecgen`, whose own comment says the consumer asked for `Primary Key, Target, Variable_1, …`. The brief's illustrative `label` was not used |
| Zero-pad the position column numbers? | **Yes** | `position_0001 … position_1000`, so `sorted(df.columns)` cannot permute the feature matrix |

## Blocking — needs the Nexus team

### 1. The real column and row ceiling

1,000 positions was chosen as a safe starting point, not because anything establishes it.
Nothing on this machine records a Nexus limit; the sibling `qecgen` project states outright
that *"the Nexus input format is unknown"* and refuses to claim compatibility. For scale,
that project has shipped a 21,866-column × 100,000-row CSV (4.38 GB), so very wide tables are
not unprecedented — but nothing says Nexus will take one.

*Changes if answered:* `window.length` in `config/default.yaml`, and whether long sequences
need windowing at all. Everything else is unaffected.

### 2. Does Nexus accept a multilabel target block?

`funsoc_pathogenicity_mechanism` emits 32 target columns, `target_00 … target_31`, because a
protein genuinely can carry several FunSoCs and collapsing them to one would discard real
labels. The brief describes "an annotated target/category column", singular.

*Changes if answered:* if Nexus needs exactly one target, the FunSoC set has to be emitted
either as 32 separate binary files or as a single most-specific-label column, both of which
lose information. Worth confirming before Monday.

### 3. Will a fixture be supplied?

`qecgen`'s standing rule is that Nexus compatibility is not claimed until an exporter passes a
fixture from the Nexus team. That rule is inherited here and no fixture exists. Every sidecar
carries `"nexus_compatibility": "NOT VERIFIED"`.

## Blocking — needs a file only you can fetch

### 4. The FELIX engineering-signature FASTA

This is the **only genuine sequence-level natural-versus-engineered ground truth** found
across all four sources. It is free; the `HTTP 403` is a Cloudflare bot challenge
(`Cf-Mitigated: challenge`), not a paywall, so a browser gets it and a script does not.

1. Open <https://pubs.acs.org/doi/10.1021/acssynbio.3c00398?goto=supporting-info>
2. Let the Cloudflare interstitial clear; click the link rendered as *sifile1 - zip file*
3. Save to `data\raw\felix\sb3c00398_si_001.zip`

It contains an XLSX of sample metadata, a CSV mapping sample IDs to engineering-signature IDs,
and a FASTA of the 1,004 engineering-signature sequences (234 unique elements).

**The XLSX and CSV column headers have not been observed by anyone.** The adapter will read
the FASTA — whose shape is fixed by the format — and refuse the other two until their real
headers are seen. Guessing a schema is the worst available error here.

Open question even once the ZIP arrives: what the honest **negative** class is. The signatures
are all engineered elements. A defensible negative would be sequence sampled from the matched
natural host reference genomes, but that needs a documented sampling rule (matched length
distribution, matched GC content, matched assembly source) or the model learns a batch effect
rather than engineering. **This decision should not be made unilaterally.**

### 5. The GEAC competition files

Gated behind a free DrivenData account and acceptance of the competition rules. Steps are in
the README and printed by `nexus_genomics.cli audit`. The adapter is written and tested
against a synthetic fixture; it will run unchanged the moment the files appear.

## Non-blocking — recorded so they are not forgotten

### 6. `indeterminate_origin` has no data to exercise it

The brief asked for a documented configuration handling indeterminate samples: exclude them
from the binary output by default, preserve them in an audit file, never silently treat them
as natural. The GUARDIAN schema does define the value — but it **never appears** in the
committed corpus (natural 257, engineered 250, indeterminate 0). So the requirement is
satisfied vacuously. If the SI ZIP turns out to contain indeterminate rows, this needs
implementing for real rather than assumed to work.

### 7. `fine_tuned` in the CodonTransformer training corpus

An undocumented boolean (955,322 False / 45,875 True). The token appears nowhere in the
package source, model card or Zenodo description. Empirically consistent with a top-10%-CSI
flag over the 17 fine-tune organisms (10% of 470,375 ≈ 47,038 vs 45,875 observed), but that is
inference. **It is never used as a target** and should not be, until the authors confirm it.

### 8. Whether the strided sample is representative

Sampling from the 6.08 GB Zenodo file uses fixed, evenly spaced byte offsets. It is
reproducible and its offsets are recorded, but it is not a uniform random sample: rows are
grouped by organism, so the organism mix in the demo reflects where the offsets land. A
uniform sample would require either the full download or a row index that Zenodo does not
provide.

### 9. Should the 39 GEAC phenotype columns enter the model table?

They are genuine numeric variables (`bacterial_resistance_*`, `copy_number_*`, `species_*`, …)
and would likely help lab attribution. They currently sit in the `.samples.csv` sidecar so the
feature block stays strictly one-position-per-column, and can be joined on `sample_id`. Moving
them into the Nexus table is a one-line change if that is preferred — but it makes the feature
block heterogeneous, which the brief's "one sequence position per feature column" preference
argues against.

### 10. A natural-vs-engineered set could be *constructed*, but should it be?

Addgene plasmids are engineered by construction; RefSeq bacterial genomes and PLSDB natural
plasmids are natural by construction. Pairing them would produce a large binary dataset
quickly. It would also be **methodologically weak**: the classes would differ by database,
assembly pipeline, length distribution and GC content as much as by engineering, so a model
would score well by learning the source rather than the biology. This was not built. If it is
ever wanted, it needs an explicit name that says it is a constructed contrast, and a
documented matching procedure.
