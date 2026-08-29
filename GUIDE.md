# Guide: what this is, and how to try it

A five-minute walkthrough. For every technical detail, [`README.md`](README.md) is the full
reference — this page is just the friendly on-ramp.

## The short version

We were handed four public genomic datasets, each described as "natural versus engineered" —
the kind of label you'd want to train a biosecurity detector on. We dug in, and **none of the
four actually supports that label.** Rather than force one onto the data anyway, we figured out
what each dataset *genuinely* measures, converted all of them into a clean, consistent CSV
format, and then went and found a fifth source that actually does contain real natural-vs-
engineered sequence data.

| Source | What it's really good for |
|---|---|
| Addgene / GEAC | Which laboratory built a given plasmid (1,314 labs) |
| FELIX / GUARDIAN (public repo) | What role a detected DNA part plays — promoter, gene, etc. |
| SeqScreen / FunSoC | Which of 32 disease-causing mechanisms a protein carries |
| CodonTransformer | A natural gene vs. that same gene rewritten by an AI model |
| **FELIX Supporting Information** (new) | **Real natural-vs-engineered labels**, from one specific paper's test data |

That last one is the good part: buried in the supplementary files of a 2024 paper
([Adler et al., *ACS Synthetic Biology*](https://doi.org/10.1021/acssynbio.3c00398)), there's a
small set of DNA sequences with genuine, source-verified labels for whether a specific
engineering test introduced them or not. It's small (315 sequences) and it answers a narrower
question than "is this natural in general," but it's the real thing — everything else here is an
honest re-framing, not a fabrication.

## Why this took care

The tempting shortcut was to just slap a 0/1 "natural/engineered" column on all five datasets and
call it done. We didn't, for two reasons:

1. **Most of the sources can't honestly support that label.** Addgene is 100% engineered
   plasmids — there's no natural class to contrast against. FunSoC flags *disease risk*, which
   isn't the same thing as *engineered*. Inventing a label the data doesn't back up would produce
   a dataset that looks fine and trains a model on nonsense.
2. **Even the one source that DOES have real labels has a trap in it.** In that fifth dataset, a
   "deletion" record describes DNA that was *removed* by an experiment — meaning it was part of
   the original organism, so it's natural, not engineered. Read the filename wrong and you'd
   flip that label backwards. We built tests specifically to catch this.

## Try it yourself

```powershell
python -m venv .venv
.venv\Scripts\python -m pip install -e ".[dev]"
.venv\Scripts\python -m nexus_genomics.cli audit                    # what's ready to convert right now
.venv\Scripts\python -m nexus_genomics.cli convert-all --profile demo
.venv\Scripts\python -m nexus_genomics.cli validate (Get-ChildItem outputs\*.ml.csv)
```

That downloads a small sample from each available source (a few MB total, no account needed for
four of the five), converts each into a CSV, and checks the result against ~26 automated
integrity rules per file. Open any `outputs\*.ml.csv` in Excel — one row per DNA/protein
sequence, one column per position, a label column, done.

## What you'll actually see

Every output file comes as a set of three:

```
some_dataset.ml.csv             the data itself — spreadsheet-ready
some_dataset.ml.manifest.json   where it came from, what was excluded, and why
some_dataset.ml.samples.csv     extra per-row detail, kept separate from the model table
```

DNA/protein letters are turned into numbers the simple way — A=1, B=2, C=3, and so on — with 0
meaning "no data here" (padding). That's a deliberate, boring choice: it's exactly what the
project asking for this data requested, and boring is good when the point is that a stranger can
read the manifest and know precisely what they're looking at.

## Confidence, not just claims

Before calling any of this finished, we had it reviewed adversarially — separate passes
specifically trying to find mistakes, then separate passes trying to disprove each thing the
first pass found. That process has caught real bugs twice now, including once catching that a
verification check was quietly checking nothing at all. See
[`reports/self_review.md`](reports/self_review.md) for the full history, and
[`reports/source_audit.md`](reports/source_audit.md) for the detailed evidence behind the "none
of the four" finding above.

## What's still open

- Two of the five sources need something only a human can fetch: Addgene needs a free account
  and to accept a competition's terms; the FELIX Supporting Information needs a browser to get
  past a bot-check (both are free, neither needs special access).
- Nobody outside this project has confirmed the exact output shape is accepted by the downstream
  system it's built for — every file says so honestly in its own manifest.

Full detail on both, and everything else, is in [`README.md`](README.md).
