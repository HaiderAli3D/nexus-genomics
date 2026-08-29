"""SeqScreen / FunSoC.

**SeqScreen is an analysis pipeline, not a labelled dataset, and it has no
natural-versus-engineered ground truth of any kind.** The words "synthetic construct",
"wild-type" and "provenance" appear zero times in the Genome Biology paper. The one field
that looks like an engineering signal -- the output column ``vbeo == 's'`` -- is a decoy: it
fires when a read's assigned taxid falls under NCBI "synthetic construct" (taxid 32630),
which flags cloning vectors and artefact GenBank records. An engineered pathogen gene hits
its natural homolog and is never flagged.

**"Sequence of concern" does not mean "engineered".** A FunSoC -- Function of Sequence of
Concern -- is a *pathogenicity mechanism*. The overwhelming majority of sequences carrying
one are entirely natural: toxins, adhesins and virulence factors that evolved. Reading the
FunSoC label as an engineering label would be a category error, and this adapter therefore
names its output ``funsoc_pathogenicity_mechanism`` and never emits a binary
natural/engineered column.

What *does* exist, openly and without credentials, is a genuinely human-biocurated label
set: the Rice Box ``FunSoCTrainingData`` folder, 32 files, one per FunSoC, each a bare list
of UniProt accessions. This adapter reads all 32 (the listing is paginated -- 20 files on
page 1 and 12 on page 2, and fetching only page 1 silently loses 12 FunSoCs), fetches the
**amino-acid** sequences from the UniProt REST API, and emits a 32-column multilabel matrix.

Three caveats travel with every file it writes, and they are in the sidecar as well as here:

* **Positive-only.** A ``0`` means "not curated as a positive for this FunSoC", not
  "verified negative". This is positive-unlabelled data and evaluating it as though the
  zeros were true negatives will overstate precision.
* **Severely imbalanced**, from 23 positives (CounterImmunoglobin) to 4,558 (DisableOrgan).
* **Annotation-biased.** The curation was seeded from UniProt keyword and GO queries over
  reviewed Swiss-Prot entries, and the published model was trained on UniProt *annotation*
  fields rather than on residues -- so its reported accuracies are not a baseline for any
  sequence-based model trained on this table.
"""

from __future__ import annotations

import os
import re
import urllib.error
import urllib.parse
import urllib.request
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any, Final

from nexus_genomics.adapters.base import Availability, LoadedSource
from nexus_genomics.cleaning import Alphabet
from nexus_genomics.common import SampleRecord, hash_text

__all__ = ["SeqScreenFunSoCAdapter"]

SHARED_NAME: Final = "60380nrh1mqfib5yic8x4damycbje838"
FOLDER_URL: Final = f"https://rice.box.com/s/{SHARED_NAME}"
DOWNLOAD_TEMPLATE: Final = (
    "https://rice.app.box.com/index.php?rm=box_download_shared_file"
    "&shared_name={shared}&file_id={file_id}"
)
UNIPROT_STREAM: Final = "https://rest.uniprot.org/uniprotkb/stream"
_BROWSER_UA: Final = "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"
_ACCESSION: Final = re.compile(r"[A-Z0-9][A-Z0-9-]{3,}")
_EXPECTED_FILES: Final = 32

_SPLITTABLE_STATUSES: Final = frozenset({400, 414})
"""Statuses that mean "your query is too long", and only those.

UniProt answers an over-long ``accession:A OR accession:B OR ...`` query with a bare 400 and
publishes no limit, so the batch size has to adapt. Every other non-2xx status -- 429, 500,
502, 503, 504 -- is transient or a server fault, and halving on those turns one failure into
an exponential burst of them.
"""


class SeqScreenFunSoCAdapter:
    """FunSoC multilabel over UniProt proteins. Pathogenicity mechanism, not provenance."""

    name = "funsoc_pathogenicity_mechanism"
    source_url = "https://gitlab.com/treangenlab/seqscreen"
    licence = (
        "SeqScreen code is GPLv3. The FunSoC curation is distributed openly by the Treangen "
        "Lab via Rice Box; no explicit licence is attached to the accession lists. UniProt "
        "sequence data is CC-BY-4.0. Cite Balaji et al., Genome Biology 2022."
    )
    target_description = (
        "32 binary targets, one per FunSoC (Function of Sequence of Concern). Multilabel, "
        "positive-only. This is a pathogenicity-mechanism annotation and is NOT an "
        "engineered-versus-natural label."
    )

    def availability(self, raw_dir: Path) -> Availability:
        cached = raw_dir / "funsoc"
        n = len(list(cached.glob("*.txt"))) if cached.exists() else 0
        if n >= _EXPECTED_FILES:
            return Availability(True, f"{n} cached FunSoC accession lists in {cached}")
        return Availability(
            True,
            f"{n}/{_EXPECTED_FILES} cached; the adapter will fetch the rest from Rice Box "
            f"and UniProt (open, no credentials)",
            manual_steps=(
                f"Only needed if the network is unavailable. Open {FOLDER_URL} in a browser, "
                f"download all {_EXPECTED_FILES} .txt files (REMEMBER PAGE 2 -- the listing "
                f"paginates at 20 and the remaining 12 FunSoCs are easy to miss), and place "
                f"them in data/raw/seqscreen/funsoc/."
            ),
        )

    def load(self, raw_dir: Path, options: Mapping[str, Any]) -> LoadedSource:
        cache = raw_dir / "funsoc"
        cache.mkdir(parents=True, exist_ok=True)
        files = _ensure_files(cache)
        if len(files) < _EXPECTED_FILES:
            raise RuntimeError(
                f"only {len(files)} of {_EXPECTED_FILES} FunSoC lists are available. "
                f"Refusing to build a partial label matrix: the missing FunSoCs would "
                f"appear as all-zero columns, which is indistinguishable from a FunSoC that "
                f"genuinely has no positives."
                + (f" Listing failed: {_DISCOVERY_ERROR[-1]}." if _DISCOVERY_ERROR else "")
                + f" {self.availability(raw_dir).manual_steps}"
            )

        funsocs = sorted(files)
        index = {name: i for i, name in enumerate(funsocs)}
        labels: dict[str, set[str]] = {}
        per_class: dict[str, int] = {}
        for funsoc in funsocs:
            accessions = _read_accessions(files[funsoc])
            per_class[funsoc] = len(accessions)
            for accession in accessions:
                labels.setdefault(accession, set()).add(funsoc)

        limit = options.get("max_records")
        wanted = _stratified(labels, funsocs, int(limit)) if limit else sorted(labels)

        sequences, unresolved = _fetch_sequences(wanted, raw_dir, int(options.get("batch", 100)))

        records: list[SampleRecord] = []
        for accession in wanted:
            sequence = sequences.get(accession)
            if not sequence:
                continue
            vector = [0] * len(funsocs)
            for funsoc in labels[accession]:
                vector[index[funsoc]] = 1
            records.append(
                SampleRecord(
                    sample_id=accession,
                    sequence=sequence,
                    label=tuple(vector),
                    source=self.name,
                    metadata={
                        "uniprot_accession": accession,
                        "funsocs": "|".join(sorted(labels[accession])),
                        "n_funsocs": len(labels[accession]),
                    },
                )
            )

        return LoadedSource(
            records=records,
            n_targets=len(funsocs),
            label_semantics={
                0: "not curated as a positive for this FunSoC (NOT a verified negative)",
                1: "curated positive for this FunSoC",
            },
            alphabet=Alphabet.AMINO_ACID,
            target_names=tuple(funsocs),
            manifest_extra={
                "target_index_to_funsoc": dict(enumerate(funsocs)),
                "positives_per_funsoc_in_source_lists": per_class,
                "positives_per_funsoc_in_this_file": {
                    name: sum(1 for r in records if r.label and r.label[index])
                    for index, name in enumerate(funsocs)
                },
                "sampling": (
                    "The capped profiles select accessions ROUND-ROBIN across the 32 FunSoCs "
                    "so that no target column ends up all-zero. That deliberately flattens "
                    "the class distribution, so a capped file is NOT distribution-"
                    "representative of the source: use positives_per_funsoc_in_this_file for "
                    "loss weighting, never the source-list counts."
                    if options.get("max_records")
                    else "No cap applied; every curated accession that resolved is included."
                ),
                "total_curated_accessions": len(labels),
                "unresolved_accessions": len(unresolved),
                "unresolved_examples": sorted(unresolved)[:20],
                "sequence_source": "UniProt REST (rest.uniprot.org), fields=accession,sequence",
                "positive_unlabelled_warning": (
                    "A 0 means 'not curated positive', NOT 'verified negative'. This is "
                    "positive-unlabelled data; treating the zeros as true negatives will "
                    "overstate precision."
                ),
                "not_an_engineering_label": (
                    "A Function of Sequence of Concern is a pathogenicity mechanism. Most "
                    "sequences carrying one are entirely natural toxins or virulence "
                    "factors. This is not, and must not be presented as, a "
                    "natural-versus-engineered label."
                ),
                "annotation_bias": (
                    "Curation was seeded from UniProt keyword/GO queries over reviewed "
                    "Swiss-Prot entries; the published model used annotation fields rather "
                    "than residues as features, so its accuracies are not a baseline here."
                ),
                "seqscreen_outputs_as_extra_variables": (
                    "SeqScreen's own per-query outputs (taxid, UniProt hit, GO terms, "
                    "bsat_hit, vfdb_hit, funsoc columns) could become additional numeric "
                    "variables for another dataset, but only by running the pipeline: it "
                    "needs a ~169 GiB database, 32-256 GB RAM and Linux. Not attempted here."
                ),
            },
            raw_input_hashes={
                name: hash_text(path.read_text(encoding="utf-8", errors="replace"))
                for name, path in sorted(files.items())
            },
            sample_metadata_fields=("uniprot_accession", "funsocs", "n_funsocs"),
        )


def _stratified(labels: Mapping[str, set[str]], funsocs: Sequence[str], limit: int) -> list[str]:
    """Cap the sample without letting any FunSoC drop out of the file entirely.

    Taking the first N accessions alphabetically left 6 of the 32 target columns all-zero,
    which is the worst kind of quiet defect: an all-zero column is indistinguishable from a
    FunSoC that genuinely has no positives, and a model trained on it learns that the class
    never occurs. This walks the FunSoCs round-robin instead, so every column carries at
    least one positive as long as the cap is at least the number of FunSoCs.

    Deterministic: both the FunSoC order and the accessions within each are sorted, so the
    same cap always yields the same sample and ``content_hash`` stays reproducible.
    """
    by_funsoc: dict[str, list[str]] = {f: [] for f in funsocs}
    for accession in sorted(labels):
        for funsoc in sorted(labels[accession]):
            by_funsoc[funsoc].append(accession)

    chosen: list[str] = []
    seen: set[str] = set()
    depth = 0
    while len(chosen) < limit:
        added = False
        for funsoc in funsocs:
            bucket = by_funsoc[funsoc]
            if depth >= len(bucket):
                continue
            added = True
            accession = bucket[depth]
            if accession not in seen:
                seen.add(accession)
                chosen.append(accession)
                if len(chosen) >= limit:
                    break
        if not added:
            break
        depth += 1
    return chosen


def _funsoc_name(filename: str) -> str:
    """``SecretedEffector_v1.0_10.29.2019.txt`` -> ``SecretedEffector``."""
    return filename.split("_v")[0]


def _ensure_files(cache: Path) -> dict[str, Path]:
    """Return {funsoc: path}, downloading anything not already cached."""
    found = {_funsoc_name(p.name): p for p in cache.glob("*.txt")}
    if len(found) >= _EXPECTED_FILES:
        return found
    try:
        listing = _discover()
    except (urllib.error.URLError, OSError) as exc:
        # Only the network is tolerated here, and the cause is carried forward. A bare
        # `except Exception` collapsed a 403 from Box, a changed listing page and a genuine
        # regex bug into the same silent "download them yourself" message, which sends
        # someone to a browser to fix a problem that is not on their end.
        _DISCOVERY_ERROR.append(repr(exc))
        return found
    for file_id, filename in listing.items():
        target = cache / filename
        if target.exists():
            found[_funsoc_name(filename)] = target
            continue
        url = DOWNLOAD_TEMPLATE.format(shared=SHARED_NAME, file_id=file_id)
        request = urllib.request.Request(url, headers={"User-Agent": _BROWSER_UA})
        with urllib.request.urlopen(request, timeout=120) as response:
            body = response.read()
        _require_accession_list(body, filename)
        _atomic_write(target, body)
        found[_funsoc_name(filename)] = target
    return found


_DISCOVERY_ERROR: list[str] = []
"""Why listing the Box folder failed, if it did. Surfaced in the RuntimeError from load()."""


def _atomic_write(path: Path, data: bytes) -> None:
    """Write via a temp sibling and os.replace, so a killed run leaves no half file.

    A truncated cache is worse than a missing one: the next run reads it, gets a protein cut
    off mid-sequence, and encodes it into the table where nothing downstream can tell it from
    a genuinely short protein.
    """
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_bytes(data)
    os.replace(tmp, path)


def _require_accession_list(body: bytes, filename: str) -> None:
    """Refuse anything that is not a list of UniProt accessions.

    Box answers a throttled shared-link download with HTTP 200 and an HTML challenge page.
    Cached blindly, that page becomes a FunSoC's accession list: the file count still reaches
    32 so the guard in `load` does not fire, and that FunSoC silently loses every positive.
    """
    lines = [ln.strip() for ln in body.decode("utf-8", "replace").splitlines() if ln.strip()]
    offender = next((ln for ln in lines if not _ACCESSION.fullmatch(ln)), None)
    if not lines or offender is not None:
        raise RuntimeError(
            f"the download for {filename} is not an accession list -- it contains "
            f"{(offender or '<empty file>')[:60]!r}. Box serves an HTML challenge page with "
            f"HTTP 200 when it throttles; caching it would silently empty that FunSoC."
        )


def _discover() -> dict[str, str]:
    """Both pages of the Box listing. Page 2 exists and is easy to miss."""
    out: dict[str, str] = {}
    for page in (1, 2):
        url = FOLDER_URL if page == 1 else f"{FOLDER_URL}?page={page}"
        request = urllib.request.Request(url, headers={"User-Agent": _BROWSER_UA})
        with urllib.request.urlopen(request, timeout=120) as response:
            html = response.read().decode("utf-8", "replace")
        ids = re.findall(r"f_\d+", html)
        names = re.findall(r'"([A-Za-z0-9_.]+_v[\d.]+_[\d.]+\.txt)"', html)
        # The listing renders ids and names in the same order; both lists must be the same
        # length or the pairing is a guess, so a mismatch is skipped rather than zipped.
        seen_ids = list(dict.fromkeys(ids))
        seen_names = list(dict.fromkeys(names))
        if len(seen_ids) == len(seen_names):
            out.update(zip(seen_ids, seen_names, strict=True))
    return out


def _read_accessions(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8", errors="replace")
    return [line.strip() for line in text.splitlines() if line.strip()]


def _fetch_batch(chunk: Sequence[str], cached: dict[str, str]) -> None:
    """One UniProt request, halving the batch and retrying if the server rejects it.

    UniProt answers a long ``accession:A OR accession:B OR ...`` query with a bare
    ``HTTP 400`` once the query grows past a limit it does not publish. A fixed batch size
    is therefore not safe: 100 worked for the demo and 200 failed on the full 13,203-accession
    run, which is precisely the shape of bug that only appears at full scale. Halving on
    failure adapts to whatever the real limit is, and recursing to a single accession means
    one genuinely bad accession costs one failed request rather than the whole batch.
    """
    if not chunk:
        return
    query = " OR ".join(f"accession:{a}" for a in chunk)
    url = (
        f"{UNIPROT_STREAM}?query={urllib.parse.quote(f'({query})')}"
        f"&format=tsv&fields=accession,sequence"
    )
    request = urllib.request.Request(url, headers={"User-Agent": "nexus-genomics"})
    try:
        with urllib.request.urlopen(request, timeout=240) as response:
            body = response.read().decode("utf-8", "replace")
    except urllib.error.HTTPError as exc:
        if exc.code not in _SPLITTABLE_STATUSES:
            # Halving only helps when the query is too long. On a 429 or a 5xx it doubles the
            # request count at every level: a 200-accession batch becomes 399 failing
            # requests, and the full 13,203-accession run would issue tens of thousands of
            # them back to back -- which is how a rate limit becomes a ban.
            raise RuntimeError(
                f"UniProt returned HTTP {exc.code} for a batch of {len(chunk)} accessions. "
                f"That is not a query-length problem, so the batch is not split. Retry "
                f"later, reduce the `batch` option, or supply "
                f"data/raw/seqscreen/uniprot_sequences.tsv yourself."
            ) from exc
        if len(chunk) == 1:
            # A single accession the server will not answer for. Left unresolved and
            # counted, never substituted with an empty sequence.
            return
        middle = len(chunk) // 2
        _fetch_batch(chunk[:middle], cached)
        _fetch_batch(chunk[middle:], cached)
        return
    for line in body.splitlines()[1:]:
        if "\t" in line:
            accession, sequence = line.split("\t", 1)
            cached[accession.strip()] = sequence.strip()


def _fetch_sequences(
    accessions: Sequence[str], raw_dir: Path, batch: int
) -> tuple[dict[str, str], set[str]]:
    """Amino-acid sequences from UniProt, cached to disk so a re-run costs nothing."""
    cache_path = raw_dir / "uniprot_sequences.tsv"
    cached: dict[str, str] = {}
    if cache_path.exists():
        for line in cache_path.read_text(encoding="utf-8").splitlines():
            if "\t" in line:
                acc, seq = line.split("\t", 1)
                cached[acc] = seq

    missing = [a for a in accessions if a not in cached]
    for start in range(0, len(missing), batch):
        _fetch_batch(missing[start : start + batch], cached)

    if missing:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        _atomic_write(
            cache_path,
            "".join(f"{a}\t{s}\n" for a, s in sorted(cached.items())).encode("utf-8"),
        )

    resolved = {a: cached[a] for a in accessions if a in cached}
    unresolved = {a for a in accessions if a not in cached}
    return resolved, unresolved
