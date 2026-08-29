"""IARPA FELIX / GUARDIAN / MIDOE.

**This source does not yield a natural-versus-engineered sequence dataset, and this adapter
refuses to pretend otherwise.** The brief called it the strongest candidate. It is not.

What the public Apache-2.0 repository actually contains, measured by parsing all 340
evidence JSONs rather than read off the documentation:

* 1,672 ``features``, of which **1,670 carry a real DNA ``sequence``** -- 1,266,189 bp,
  11 to 52,182 bp, alphabet ``{A,C,G,N,T}`` plus lowercase. So the repo does ship sequence.
* 627 ``detections`` over exactly 100 blinded ``Y###`` samples, 507 of which carry the
  ``engineered`` enum. Sample-level, that enum looks nicely balanced: natural 257,
  engineered 250. ``indeterminate_origin`` never appears in the committed data at all.
* But that enum is **detector output, not ground truth**. The schema describes it as
  "whether this *evidence* is natural, engineered, or of indeterminate origin"; every
  detection names an ``agent`` (seven of them, and they disagree) and most carry a
  ``confidence``. Training on it distils GUARDIAN's 2022 calls, false negatives included.
* And the balance evaporates at sequence level. Only **280 of the 1,670 sequences** can be
  linked to a call at all (via ``region.feature`` or ``alteration.derivedFeature``), and
  **all 280 are ``engineered``**. There is no negative class, because a detector that finds
  nothing engineered emits no feature to attach a sequence to. A binary file built here
  would have one class in it.

Two constructions were considered and rejected:

* **NCBI BioProject PRJNA607328** gives a per-BioSample ``biol_stat`` (179 engineered /
  52 natural / 59 unresolved, and the source vocabulary contains the misspelling
  ``enginered mutant``). Rejected because the unit is a multi-megabase genome that is
  >99.9% natural sequence: windowing it would stamp "engineered" onto overwhelmingly
  natural windows.
* **Pairing the 280 engineered features against sampled natural reference genomes.**
  Rejected because the two classes would then differ by origin, assembly pipeline and
  length distribution as much as by engineering, so a model would learn the batch effect.

What *is* emitted is a deliberately narrow, honestly named multiclass set: the annotated
**genetic part role** of each feature (CDS, promoter, terminator, origin_of_replication,
...). Those annotations come from GUARDIAN's targeted search, so they are pipeline output
too -- but "what kind of genetic part is this sequence" is a claim the annotation actually
makes, whereas "this organism was engineered" is not one any per-sequence field here
supports.

The genuine sequence-level ground truth exists, but not here: it is the ACS Supporting
Information ZIP, which needs one manual browser download. See :func:`availability`.
"""

from __future__ import annotations

import json
import tarfile
import urllib.request
from collections.abc import Mapping
from pathlib import Path
from typing import Any, Final

from nexus_genomics.adapters.base import Availability, LoadedSource, stable_class_index
from nexus_genomics.cleaning import Alphabet
from nexus_genomics.common import SampleRecord, hash_file

__all__ = ["FelixGuardianAdapter"]

TARBALL_URL: Final = "https://codeload.github.com/raytheonbbn/midoe/tar.gz/refs/heads/main"
TARBALL_NAME: Final = "midoe-main.tar.gz"
SI_ZIP_NAME: Final = "sb3c00398_si_001.zip"

_SI_STEPS: Final = f"""\
The genuine sequence-level natural-versus-engineered ground truth is NOT in this repository.
It lives in the ACS Supporting Information and is handled by a separate adapter,
`felix_signatures` (nexus_genomics/adapters/felix_signatures.py), which reads all three of
its files and explains why a deletion's sequence is NATURAL. The SI is FREE but sits behind
a Cloudflare bot challenge (HTTP 403 with 'Cf-Mitigated: challenge' to any script), so it
needs a browser:

  1. Open https://pubs.acs.org/doi/10.1021/acssynbio.3c00398?goto=supporting-info
     in a normal desktop browser and let the Cloudflare interstitial clear.
  2. Under 'Supporting Information', click the link rendered as 'sifile1 - zip file'.
  3. Save it as data/raw/felix_signatures/{SI_ZIP_NAME}

Then run:  python -m nexus_genomics.cli convert felix_signatures --profile full\
"""


class FelixGuardianAdapter:
    """Genetic part roles from the MIDOE evidence corpus. Not an engineering label."""

    name = "felix_guardian_part_role"
    source_url = "https://github.com/raytheonbbn/midoe"
    licence = "Apache-2.0 (code and committed evidence). Redistribution permitted with attribution."
    target_description = (
        "Annotated genetic part role of a detected sequence feature (multiclass). "
        "GUARDIAN pipeline annotation, NOT experimentally verified ground truth, and NOT "
        "an engineered-versus-natural label."
    )

    def availability(self, raw_dir: Path) -> Availability:
        tar = raw_dir / TARBALL_NAME
        if tar.exists():
            return Availability(True, f"{tar.name} present ({tar.stat().st_size:,} bytes)")
        return Availability(
            True,
            "not downloaded yet; the adapter will fetch it (1.2 MB, Apache-2.0, no credentials)",
            manual_steps=_SI_STEPS,
            expected_files=(TARBALL_NAME,),
        )

    def fetch(self, raw_dir: Path) -> Path:
        """Download the repository tarball once, rather than 340 separate API calls."""
        raw_dir.mkdir(parents=True, exist_ok=True)
        tar = raw_dir / TARBALL_NAME
        if not tar.exists():
            request = urllib.request.Request(TARBALL_URL, headers={"User-Agent": "nexus-genomics"})
            with urllib.request.urlopen(request, timeout=180) as response:
                tar.write_bytes(response.read())
        return tar

    def load(self, raw_dir: Path, options: Mapping[str, Any]) -> LoadedSource:
        tar_path = self.fetch(raw_dir)
        features, call_index, stats = _parse(tar_path)

        roled = [f for f in features if f.get("role")]
        # Sorted by id so a cap takes a deterministic subset rather than whatever order the
        # tar happened to yield.
        roled.sort(key=lambda f: (f["_file"], f["id"]))

        # The class vocabulary comes from the WHOLE corpus, before any cap. Building it from
        # the capped sample would renumber every class that sorts after a class the cap
        # happened to exclude, so `target == 2` would mean origin_of_replication in the demo
        # file and engineered_tag in the full one -- and a model trained on one and scored
        # against the other would be silently graded against the wrong classes.
        classes = stable_class_index([str(f["role"]) for f in roled])

        limit = options.get("max_records")
        if limit:
            roled = roled[: int(limit)]
        records: list[SampleRecord] = []
        for feature in roled:
            fid = f"{feature['_file'].removesuffix('.json')}::{feature['id']}"
            call = call_index.get((feature["_file"], feature["id"]))
            records.append(
                SampleRecord(
                    sample_id=fid,
                    sequence=str(feature["sequence"]),
                    label=(classes[str(feature["role"])],),
                    source=self.name,
                    metadata={
                        "feature_name": str(feature.get("name", ""))[:120],
                        "part_role": str(feature["role"]),
                        "genbank_source": str(feature.get("source", "")),
                        "source_start": feature.get("sourceStart", ""),
                        "source_end": feature.get("sourceEnd", ""),
                        "source_orientation": str(feature.get("sourceOrientation", "")).split("#")[
                            -1
                        ],
                        # Carried for audit only. It is a detector assertion and is single
                        # class across every feature it covers, which is precisely why it
                        # is not the target.
                        "guardian_call_not_ground_truth": call or "",
                    },
                )
            )

        return LoadedSource(
            records=records,
            n_targets=1,
            label_semantics={
                v: (
                    f"{k} -- a genetic PART type annotated by GUARDIAN's targeted search. "
                    f"This is what KIND of element the sequence is, not where it came from."
                )
                for k, v in classes.items()
            },
            alphabet=Alphabet.DNA_NUCLEOTIDE,
            manifest_extra={
                "corpus_census": stats,
                "why_not_natural_vs_engineered": (
                    "Only 280 of the 1,670 sequences in this repository can be linked to "
                    "the GUARDIAN 'engineered' enum at all, and all 280 of those are the "
                    "engineered value. A binary file built from this source would contain "
                    "exactly one class. The 257-natural / 250-engineered balance visible in "
                    "the corpus is per-detection (per whole sample), not per-sequence, and "
                    "is in any case detector output rather than ground truth."
                ),
                "ground_truth_route": _SI_STEPS,
                "class_map": {str(v): k for k, v in classes.items()},
                "no_subset_of_these_classes_is_a_provenance_contrast": (
                    "Two class names contain the word 'engineered' (engineered_region, "
                    "engineered_tag). They are PART TYPES, not provenance calls. Deriving a "
                    "binary such as target.isin([engineered_region, engineered_tag]) would "
                    "be wrong: the CDS, promoter, origin_of_replication and terminator rows "
                    "are parts that GUARDIAN annotated INSIDE engineered constructs, so "
                    "labelling them 'natural' would be exactly the fabrication this adapter "
                    "exists to avoid."
                ),
                "rows_with_a_guardian_call": sum(
                    1 for f in roled if call_index.get((f["_file"], f["id"]))
                ),
                "guardian_call_coverage_note": (
                    "The call-linked features and the role-annotated features are disjoint "
                    "sets in this corpus, so guardian_call_not_ground_truth is empty on "
                    "every emitted row. That is the corpus topology, not a broken join."
                ),
                "source_coordinate_convention": (
                    "source_start/source_end are copied verbatim from MIDOE sourceStart/"
                    "sourceEnd and are 0-based half-open in the coordinates of the GenBank "
                    "accession in genbank_source; source_orientation says whether the "
                    "sequence is the reverse complement of that span. window_start/"
                    "window_end are a DIFFERENT convention: 1-based inclusive, in the "
                    "coordinates of this feature's own sequence."
                ),
            },
            raw_input_hashes={TARBALL_NAME: hash_file(tar_path)},
            sample_metadata_fields=(
                "feature_name",
                "part_role",
                "genbank_source",
                "source_start",
                "source_end",
                "source_orientation",
                "guardian_call_not_ground_truth",
            ),
        )


_Feature = dict[str, Any]
_CallIndex = dict[tuple[str, str], str]


def _parse(tar_path: Path) -> tuple[list[_Feature], _CallIndex, dict[str, Any]]:
    """Read every evidence JSON out of the tarball in one pass."""
    features: list[dict[str, Any]] = []
    calls: dict[tuple[str, str], set[str]] = {}
    n_detections = 0
    detection_calls: dict[str, int] = {}
    samples: set[str] = set()

    with tarfile.open(tar_path) as archive:
        members = sorted(
            (
                m
                for m in archive.getmembers()
                if "guardian_evidence/" in m.name and m.name.endswith(".json")
            ),
            key=lambda m: m.name,
        )
        for member in members:
            handle = archive.extractfile(member)
            if handle is None:
                continue
            payload = json.loads(handle.read().decode("utf-8"))
            filename = member.name.split("/")[-1]

            for feature in payload.get("features", []):
                if "sequence" in feature:
                    features.append({**feature, "_file": filename})

            for detection in payload.get("detections", []):
                n_detections += 1
                if detection.get("sample"):
                    samples.add(str(detection["sample"]))
                if "engineered" in detection:
                    key = str(detection["engineered"]).split("#")[-1]
                    detection_calls[key] = detection_calls.get(key, 0) + 1
                for contig in detection.get("contigs", []) or []:
                    for region in contig.get("regions", []) or []:
                        if "feature" in region and "engineered" in region:
                            calls.setdefault((filename, str(region["feature"])), set()).add(
                                str(region["engineered"]).split("#")[-1]
                            )
                for alteration in detection.get("alterations", []) or []:
                    if "derivedFeature" in alteration and "engineered" in alteration:
                        calls.setdefault((filename, str(alteration["derivedFeature"])), set()).add(
                            str(alteration["engineered"]).split("#")[-1]
                        )

    # A feature with two conflicting calls is recorded as CONFLICT rather than resolved by
    # picking one: two detectors disagreeing is a fact about the corpus, not a tie to break.
    flat = {k: (next(iter(v)) if len(v) == 1 else "CONFLICT") for k, v in calls.items()}
    stats = {
        "evidence_files": len(members),
        "features_with_sequence": len(features),
        "total_bp": sum(len(f["sequence"]) for f in features),
        "detections": n_detections,
        "distinct_samples": len(samples),
        "detection_level_calls": detection_calls,
        "features_linkable_to_a_call": _linkable(flat, features),
        "feature_level_call_distribution": _distribution(flat, features),
    }
    return features, flat, stats


def _sequence_keys(features: list[_Feature]) -> set[tuple[str, str]]:
    """Keys of the features that actually carry a sequence.

    A call keyed to a feature that has no sequence cannot become a row, so the census must
    intersect against this set rather than counting the call index. Getting that wrong is
    how the headline number -- 280 of 1,670, all of them `engineered` -- would silently
    become "every call is linkable", which is the opposite of the finding this adapter
    exists to record.
    """
    return {(f["_file"], f["id"]) for f in features}


def _linkable(flat: _CallIndex, features: list[_Feature]) -> int:
    have = _sequence_keys(features)
    return sum(1 for key in flat if key in have)


def _distribution(flat: _CallIndex, features: list[_Feature]) -> dict[str, int]:
    have = _sequence_keys(features)
    out: dict[str, int] = {}
    for key, value in flat.items():
        if key in have:
            out[value] = out.get(value, 0) + 1
    return out
