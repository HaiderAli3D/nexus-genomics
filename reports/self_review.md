# Self-review

After the pipeline was working and all outputs validated, the project was reviewed as if it
were another engineer's work: four independent reviewers, one per failure-mode dimension
(label integrity; encoding and off-by-one; robustness and scale; documentation accuracy),
then an independent adversarial verifier per finding whose job was to *refute* it.

**43 findings, 28 confirmed, 10 partially confirmed, 5 refuted.** 47 agents, 979 tool calls.
Everything confirmed or partially confirmed has been fixed. The refutations are kept below,
because a finding that turned out to be wrong is also a result.

The single most valuable finding is the first one: the round-trip check was a tautology.

## Fixed — label integrity

| Finding | Why it mattered | Fix |
|---|---|---|
| **The round-trip check compared the encoder against itself.** `row_sequence_blake2b` hashed `decode(placed.tokens)`, and the validator recomputed exactly that. | An off-by-one in `place` — `tokens[1:width+1]` — would drop the first residue of every row in every file while all 25 sampled rows matched and the check reported *pass*. | Hash the **cleaned source slice** (`cleaned.sequence[start-1:end]`) instead, as `source_slice_blake2b`. `test_the_round_trip_catches_an_off_by_one_in_the_encoder` injects that exact mutation and asserts the validator now fails. |
| Class indices were assigned **after** the `max_records` cap. | `target == 2` meant `origin_of_replication` in the demo file and `engineered_tag` in the full one. A model trained on one and scored against the other is graded against the wrong classes, silently. | Build the vocabulary from the whole corpus before capping (FELIX), and from the `train_labels.csv` header (Addgene). Regression test asserts demo and full agree. |
| Demo and full CodonTransformer files reused `ct_000000` for different genes. | `pd.concat([demo, full])` produced hundreds of duplicate primary keys pointing at different sequences. `require_unique` only checks within one file. | Derive the pair id from the source row's own content (organism + GeneID + description + sequence digest). Verified: 0 colliding ids across the two shipped files. |
| An all-zero label row would have been attributed to the first lab by `idxmax`. | A fabricated label, produced by a function that looks like it is reading one. | Refuse any row that is not exactly one-hot. |
| A duplicated `sequence_id` would make a pandas row lookup return a DataFrame. | `str()` of that is a repr of a Series, so the sequence column fills with `'0    ACGT\nName: ...'` — a well-formed CSV of garbage. | Refuse duplicates explicitly. |
| A single-class target passed validation and was reported as `class_balance = 1.0`. | Byte-identical to a perfectly balanced binary in the report. | Fail when a single-target column has fewer than two distinct values; emit `null` balance and `n_distinct_classes`. |
| FunSoC's `positives_per_funsoc` counted the source lists, not the emitted table. | A reader weighting a loss from it would be out by 7× on 32 columns at once (4,558 vs the file's 38). | Renamed to `..._in_source_lists`, added `..._in_this_file`, and recorded that the capped profile is round-robin stratified and therefore *not* distribution-representative. |
| Two FELIX class names contain the word "engineered". | `target.isin([engineered_region, engineered_tag])` looks like it recovers the promised binary. It would label CDS, promoters and origins inside engineered constructs as "natural" — the exact fabrication the project exists to avoid. | Every class meaning now says it is a genetic part type; a manifest field states that no subset of the eight is a provenance contrast. Test asserts both. |
| Row order in the CodonTransformer files predicts the class perfectly (rows alternate). | A bare `reset_index()` before dropping columns leaks the label completely; `KFold(shuffle=False)` skews. | Documented in the manifest as `row_order_predicts_the_class`, alongside the note that `sample_id` and the sidecar's `member` column also name the target and are not features. |

## Fixed — encoding and accounting

| Finding | Fix |
|---|---|
| `stride > width` silently dropped the residues between windows and reported `removed=0`. | Refused, with an error saying how much would have been lost. |
| The `max_windows_per_sample` cap discarded residues but every row recorded `residues_removed_by_truncation = 0`. | New `residues_dropped_by_cap` in the manifest and `residues_dropped_by_window_cap` per sample. A 80-nt sequence capped to 3 windows of 8 now reports the 56 residues it lost. |
| With `stride < width`, every window past the end was a strict suffix of an earlier one — pure duplication, weighting sequence tails far above their heads. | Stop once a window reaches the end. |
| `validate_file` **crashed** with `ValueError: cannot convert float NaN to integer` on an empty feature cell — the exact corruption a spreadsheet round-trip produces. | Report and return, like the token-range branch. |
| The validator ran decode-dependent checks after the padding invariant failed, printing confirmatory passes beneath a failure that invalidates them. | Return early with `remaining_checks_skipped`. |
| FELIX `source_start`/`source_end` (0-based half-open, GenBank coordinates) sat beside `window_start`/`window_end` (1-based inclusive, feature coordinates) with nothing distinguishing them, and `sourceOrientation` was discarded. | Carry `source_orientation`; document both conventions in the manifest. |

## Fixed — robustness

| Finding | Fix |
|---|---|
| The UniProt retry split the batch on **any** non-2xx. On a 429 a 200-accession batch becomes 399 failing requests; the full run would issue tens of thousands back to back. | Split only on 400/414 (query too long). Anything else raises with advice. |
| Raw caches were written non-atomically. A truncated `uniprot_sequences.tsv` line still contains a tab, so the next run reads a protein cut off mid-sequence and encodes it. | All cache writes go through a temp sibling + `os.replace`. |
| Box answers a throttled download with HTTP 200 and an HTML challenge page, which would be cached as a FunSoC's accession list — the file count still reaches 32, so that FunSoC silently loses every positive. | Validate the body is an accession list before caching. |
| `except Exception` around the Box listing collapsed a 403, a changed page and a regex bug into the same "download them yourself" message. | Catch only network errors, and carry the cause into the error the user sees. |
| Validation took **87 s** and **1.33 GB** on the full files: a 1,033-element Series rebuilt per target column per row, the table materialised four times, and the whole file slurped twice for a digest it had already streamed. | Vectorised the label and decode paths, replaced the row materialisation with a streaming shape read, reused the streamed digest. Now **13.6 s / 275 MB** (FunSoC full) and **32.6 s / 669 MB** (CodonTransformer full). |
| `convert-all` exited 0 on any failure, so an out-of-memory crash looked identical to a source awaiting a manual download. | A missing gated file still exits 0; anything else exits 1. Tested both ways. |
| No size guard on the writer, which buffers the whole table to hash it. | Warn above 50 M cells. |

## Fixed — documentation

Every number in the docs was checked against the artefacts. Corrections:

- Nucleotide identity of CodonTransformer pairs is **73–94%, median 81%** (was written as 77–87%), measured over the 145 untruncated demo pairs.
- The shipped files cover **3 organisms (demo) and 9 (full)**, not the source's 17. Homo sapiens alone is 36% of the full file.
- Download volume is **5 MB demo / 160 MB full** for CodonTransformer and 5 MB for FunSoC, not "~8 MB" and "~20 MB".
- `--full` is not a CLI option; the adapter was telling users to run `--profile full`.
- The FELIX manual-download instruction pointed at `data/raw/felix_guardian/`, which the pipeline never reads. It is `data/raw/felix/`.
- Nothing documented how to build the three `_full` files the README advertises.
- `reports/validation_report.*` and `reports/source_audit_status.md` are now gitignored: once the GEAC files are present they carry per-class counts derived from `train_labels.csv` and digests of both gated files, which the competition rules forbid redistributing.
- The "374 duplicated sequences" figure is now labelled as coming from the competition's documentation, not from this project's validator — the gated files have never been on this machine.
- The categorical-token warning is in every `.ml.manifest.json`, **not** in the `.ml.samples.csv`; the docs said "every sidecar".
- The atomicity claim was overstated: the group is *staged* atomically, but the commit is one `os.replace` per file, so a kill between two replaces can leave a new table beside an old sidecar. That state is detectable via `content_hash`, and the docs now say so.
- Zenodo 12509224 has 1,001,197 rows in total, of which 470,375 belong to the 17 fine-tuned organisms. Two places quoted these as if they were the same number.
- The one-line description in the README, `pyproject.toml` and `--help` asserted "Nexus-ingestable", which is the one claim the project says it will not make.
- `nexus_csv.py` said "Four files share a stem" and listed three; the README named two modules that do not exist (`windows`, `report`).

## Refuted

Kept because they were checked and found not to be defects:

- *"CodonTransformer silently delivers 35% fewer pairs than `max_pairs`"* — `max_pairs` is a ceiling, and the shortfall is already visible in `sampling.rows_recovered`.
- *"FELIX class names invite an inverted binary with nothing marking them as part roles"* — the target description and adapter docstring already said so; the fix strengthened it anyway.
- *"`_parse_chunk` uses a module-level mutable global"* — real, but single-threaded and read on the next statement; no reachable wrong result.
- *"Addgene `.tolist()` duplication blows memory"* — the profile that would trigger it needs a gated file, and the frame dominates either way.
- *"`ruff format --check .` fails on the README's own code block"* — it does not; both pass.

One finding was downgraded rather than refuted: the single-class validation gap is real, but no shipped adapter or config setting can reach it, because a surviving lone class always gets index 0 and trips the all-zero guard. It is fixed regardless, since `validate` is also usable on files this pipeline did not write.

## After the fixes

```
104 tests passed
ruff check   All checks passed!
ruff format  29 files left unchanged
mypy --strict  Success: no issues found in 21 source files
validate     6 files, 124 checks, 0 failed
reproduce    3/3 demo files byte-identical across independent runs
```

Manual spot-checks, re-run against the sources rather than against the pipeline's own output:
6 of 6 CodonTransformer pairs translate to the identical protein; 6 of 6 FELIX rows decode
exactly to their source feature in the MIDOE tarball; 4 of 4 FunSoC rows decode exactly to
the sequence UniProt returns for that accession today.
