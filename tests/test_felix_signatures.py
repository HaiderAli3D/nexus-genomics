"""The FELIX recorded-event provenance contrast, and the traps in its labelling."""

from __future__ import annotations

import zipfile
from pathlib import Path

import pytest

from nexus_genomics.adapters.felix_signatures import (
    SI_ZIP_NAME,
    FelixSignaturesAdapter,
    _element_type,
    _length_baseline,
    _parse_fasta,
)
from nexus_genomics.cleaning import Alphabet

# One record of each shape that matters, mirroring the real file's naming exactly.
FASTA = """\
>LBNL_IF00015_AmpR_861bp_plasmid.element_Staphylococcus.aureus
ATGAGTATTCAACATTTCCGTGTCGCC
>LBNL_IF00016_AmpR_861bp_insertion_Staphylococcus.aureus
ATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATT
>LBNL_IF00275_PureD_promoter_16bp_deletion_Pseudomonas.putida
TTGGCCCGCAACCTGC
>LBNL_IF00275_PureD_promoter_16bp_deletion_Pseudomonas.putida_flank2
CGATGCCAGAGCGGGCACCGCACAGGTT
>LBNL_IF00500_gene_900bp_inversion_Escherichia.coli
ACGTACGTACGTACGTACGT
"""

SIGNATURES = """\
Sample_ID,Master_Sample_ID,IF_Element
Y302,MS818,LBNL_IF00015_AmpR_861bp_plasmid.element_Staphylococcus.aureus
Y304,MS819,LBNL_IF00016_AmpR_861bp_insertion_Staphylococcus.aureus
Y302,MS818,LBNL_IF00275_PureD_promoter_16bp_deletion_Pseudomonas.putida
"""


def make_zip(tmp_path: Path, fasta: str = FASTA) -> Path:
    raw = tmp_path / "felix_signatures"
    raw.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(raw / SI_ZIP_NAME, "w") as archive:
        archive.writestr("guardian_felix_t_and_e_batch_4_signature_sequences.fa", fasta)
        archive.writestr("guardian_felix_t_and_e_batch_4_sample_signatures.csv", SIGNATURES)
    return raw


def labels(tmp_path: Path, fasta: str = FASTA) -> dict[str, int]:
    loaded = FelixSignaturesAdapter().load(make_zip(tmp_path, fasta), {})
    return {r.sample_id: (r.label or (-1,))[0] for r in loaded.records}


def test_a_deleted_sequence_is_host_context_not_introduced_by_this_event(tmp_path: Path) -> None:
    """A deletion says what this event removed, not whether an earlier event introduced it."""
    got = labels(tmp_path)
    assert got["LBNL_IF00275_PureD_promoter_16bp_deletion_Pseudomonas.putida"] == 0


def test_an_insertion_or_plasmid_element_is_introduced_by_the_recorded_event(
    tmp_path: Path,
) -> None:
    """These element types identify the event's introduced cargo without a global claim."""
    got = labels(tmp_path)
    assert got["LBNL_IF00015_AmpR_861bp_plasmid.element_Staphylococcus.aureus"] == 1
    assert got["LBNL_IF00016_AmpR_861bp_insertion_Staphylococcus.aureus"] == 1


def test_a_flank_is_host_context_not_introduced_by_the_recorded_event(tmp_path: Path) -> None:
    """A flank is locus-side context; its earlier provenance is not established here."""
    got = labels(tmp_path)
    assert got["LBNL_IF00275_PureD_promoter_16bp_deletion_Pseudomonas.putida_flank2"] == 0


def test_a_rearrangement_is_excluded_rather_than_assigned_a_class(tmp_path: Path) -> None:
    """For an inversion only the orientation is described by the operation type.

    Neither event-relative binary class identifies the rearranged sequence itself, so it is
    dropped and counted rather than guessed.
    """
    loaded = FelixSignaturesAdapter().load(make_zip(tmp_path), {})
    assert "LBNL_IF00500_gene_900bp_inversion_Escherichia.coli" not in {
        r.sample_id for r in loaded.records
    }
    assert loaded.manifest_extra["excluded_element_types"] == {"inversion": 1}


def test_an_insertion_site_is_excluded_when_sequence_provenance_is_unsupported(
    tmp_path: Path,
) -> None:
    """A locus name alone proves neither introduced cargo nor pre-existing host context."""
    name = "LBNL_IF00501_amyE_CDSpartial-5prime_30bp_insertion.site_Bacillus.subtilis"
    fasta = f">{name}\nACGTACGTACGT\n"

    loaded = FelixSignaturesAdapter().load(make_zip(tmp_path, fasta), {})

    assert loaded.records == []
    assert loaded.manifest_extra["excluded_element_types"] == {"insertion.site": 1}
    reason = loaded.manifest_extra["excluded_element_type_reasons"]["insertion.site"]
    assert "recorded event" in reason
    assert "does not show" in reason


def test_the_element_type_is_read_from_the_slot_after_the_length_token() -> None:
    """`AmpR_861bp_insertion` and `AmpR_861bp_plasmid.element` differ only in that slot.

    Searching the whole name for a keyword would also match the gene name or the organism --
    and the type is what decides the label.
    """
    assert _element_type("LBNL_IF00016_AmpR_861bp_insertion_Staphylococcus.aureus") == "insertion"
    assert (
        _element_type("LBNL_IF00015_AmpR_861bp_plasmid.element_Staphylococcus.aureus")
        == "plasmid.element"
    )
    # The flank suffix must not become the type.
    assert (
        _element_type("LBNL_IF00275_PureD_promoter_16bp_deletion_Pseudomonas.putida_flank1")
        == "deletion"
    )


def test_the_element_type_is_found_past_an_intervening_descriptor_token() -> None:
    """A confirmed real-file bug: some names carry a descriptor between length and type.

    ``lacZalpha_168bp_CDSpartial_insertion_...`` and ``voz1_9bp_607-615delACGACCTCC_..._deletion_...``
    put a partial-CDS or inline indel descriptor in the slot right after the length token, with
    the real type further along. Reading only the immediately-following slot returns the
    descriptor and silently excludes a record whose real type is recognised -- confirmed against
    the published SI for 6 of 335 elements, including this exact deletion. Getting this specific
    one wrong is worse than the others: a deletion signature is host context the recorded event
    removed (class 0, natural), so misreading its type as unrecognised would drop a genuinely
    natural record from the dataset entirely rather than merely mislabel it -- but the failure
    mode to guard hardest against is the reverse, a natural sequence ending up labelled 1.
    """
    assert (
        _element_type("LBNL_IF00390_lacZalpha_168bp_CDSpartial_insertion_Escherichia.coli")
        == "insertion"
    )
    assert (
        _element_type(
            "LBNL_IF00265_voz1_9bp_607-615delACGACCTCC_DDLQE29-33EE_deletion_Oryza.sativa"
        )
        == "deletion"
    )


def test_an_element_with_an_intervening_descriptor_still_gets_its_real_label(
    tmp_path: Path,
) -> None:
    """The parsing fix must reach the emitted label, not just the helper function."""
    fasta = (
        ">LBNL_IF00390_lacZalpha_168bp_CDSpartial_insertion_Escherichia.coli\n"
        "ATGACCATGATTACGCCAAGCTT\n"
        ">LBNL_IF00265_voz1_9bp_607-615delACGACCTCC_DDLQE29-33EE_deletion_Oryza.sativa\n"
        "TTGGCCCGC\n"
    )
    got = labels(tmp_path, fasta)
    assert got["LBNL_IF00390_lacZalpha_168bp_CDSpartial_insertion_Escherichia.coli"] == 1
    assert got["LBNL_IF00265_voz1_9bp_607-615delACGACCTCC_DDLQE29-33EE_deletion_Oryza.sativa"] == 0
    loaded = FelixSignaturesAdapter().load(make_zip(tmp_path, fasta), {})
    assert loaded.manifest_extra["excluded_element_types"] == {}


def test_a_header_that_lost_its_marker_cannot_leak_into_the_previous_sequence() -> None:
    """Two records in the published SI are missing the '>' on their header line.

    A naive parser appends the header text and the following sequence to the PREVIOUS record.
    In the real file that inflates one 16 bp signature to 202,036 bp -- a well-formed FASTA
    record that is not the sequence it claims to be.
    """
    broken = (
        ">LBNL_IF00275_PureD_promoter_16bp_deletion_Pseudomonas.putida\n"
        "TTGGCCCGCAACCTGC\n"
        "IF00275_PURED_PROMOTER_16BP_DELETION_PSEUDOMONAS.PUTIDA_FLANK1\n"
        "CGATGCCAGAGCGGGCACCGCACAGGTT\n"
    )
    sequences, malformed = _parse_fasta(broken)
    assert sequences["LBNL_IF00275_PureD_promoter_16bp_deletion_Pseudomonas.putida"] == (
        "TTGGCCCGCAACCTGC"
    )
    assert len(malformed) == 1
    assert "FLANK1" in malformed[0]
    # The orphaned sequence is dropped, not attached to something.
    assert len(sequences) == 1


def test_the_quarantined_records_are_named_in_the_manifest(tmp_path: Path) -> None:
    """Dropping data silently is the failure; dropping it and saying so is the fix."""
    broken = FASTA + "IF99999_ORPHAN_100BP_INSERTION_SYNTHETIC\nACGTACGT\n"
    loaded = FelixSignaturesAdapter().load(make_zip(tmp_path, broken), {})
    quarantined = loaded.manifest_extra["malformed_fasta_records_quarantined"]
    assert any("IF99999" in q for q in quarantined)


def test_the_length_only_baseline_is_published(tmp_path: Path) -> None:
    """Introduced elements are far longer than flanks, so length alone scores highly.

    Publishing the number is what stops someone reporting a length detector as an
    engineering detector.
    """
    loaded = FelixSignaturesAdapter().load(make_zip(tmp_path), {})
    baseline = loaded.manifest_extra["length_only_baseline_accuracy"]
    assert 0.0 < baseline <= 1.0
    assert "beat the baseline" in loaded.manifest_extra["length_confounder_warning"]


def test_the_baseline_is_1_0_when_length_separates_the_classes_perfectly() -> None:
    """A sanity check on the statistic itself, so the reassuring number cannot be wrong."""
    from nexus_genomics.common import SampleRecord

    records = [
        SampleRecord("a", "A" * 10, (0,), "t", {}),
        SampleRecord("b", "A" * 12, (0,), "t", {}),
        SampleRecord("c", "A" * 100, (1,), "t", {}),
        SampleRecord("d", "A" * 120, (1,), "t", {}),
    ]
    assert _length_baseline(records) == 1.0


def test_both_classes_survive_a_record_cap(tmp_path: Path) -> None:
    """Element ids cluster by type, so a plain head would skew the demo file's balance."""
    loaded = FelixSignaturesAdapter().load(make_zip(tmp_path), {"max_records": 2})
    assert len(loaded.records) == 2
    assert {(r.label or ())[0] for r in loaded.records} == {0, 1}


def test_every_emitted_label_claim_is_explicitly_event_relative(tmp_path: Path) -> None:
    """Global natural/wild-type wording would exceed what operation names establish."""
    loaded = FelixSignaturesAdapter().load(make_zip(tmp_path), {})
    assert loaded.alphabet is Alphabet.DNA_NUCLEOTIDE
    assert loaded.n_targets == 1
    claims = [
        FelixSignaturesAdapter.target_description,
        *loaded.label_semantics.values(),
        *(str(record.metadata["label_basis"]) for record in loaded.records),
        str(loaded.manifest_extra["label_provenance_scope"]),
    ]
    assert all("recorded event" in claim.lower() for claim in claims)
    combined = " ".join(claims).lower()
    assert "wild type" not in combined
    assert "natural" not in combined
    assert "untouched" not in combined


def test_a_missing_zip_explains_the_browser_download(tmp_path: Path) -> None:
    from nexus_genomics.adapters.base import GatedSourceUnavailableError

    status = FelixSignaturesAdapter().availability(tmp_path / "empty")
    assert not status.ready
    assert "pubs.acs.org" in status.manual_steps
    assert status.expected_files == (SI_ZIP_NAME,)
    with pytest.raises(GatedSourceUnavailableError, match="Cloudflare"):
        FelixSignaturesAdapter().load(tmp_path / "empty", {})


def test_an_unexpected_zip_layout_is_refused_rather_than_guessed(tmp_path: Path) -> None:
    raw = tmp_path / "felix_signatures"
    raw.mkdir(parents=True)
    with zipfile.ZipFile(raw / SI_ZIP_NAME, "w") as archive:
        archive.writestr("something_else.txt", "hello")
    with pytest.raises(ValueError, match="Refusing to guess"):
        FelixSignaturesAdapter().load(raw, {})


def test_a_signature_csv_with_the_wrong_columns_is_refused(tmp_path: Path) -> None:
    raw = tmp_path / "felix_signatures"
    raw.mkdir(parents=True)
    with zipfile.ZipFile(raw / SI_ZIP_NAME, "w") as archive:
        archive.writestr("guardian_felix_t_and_e_batch_4_signature_sequences.fa", FASTA)
        archive.writestr("guardian_felix_t_and_e_batch_4_sample_signatures.csv", "a,b\n1,2\n")
    with pytest.raises(ValueError, match="Refusing to guess a schema"):
        FelixSignaturesAdapter().load(raw, {})


def test_an_element_and_its_flanks_are_in_opposite_classes_and_share_an_if_number(
    tmp_path: Path,
) -> None:
    """The leakage the manifest warns about, asserted so the warning stays true.

    A deletion and its flanks sit beside each other in the genome and here they are class 0
    and class 0 -- but an INSERTION and its flanks straddle the boundary, which is why a
    split has to group by IF number rather than by row.
    """
    fasta = (
        ">LBNL_IF00016_AmpR_861bp_insertion_Staphylococcus.aureus\nACGTACGTAC\n"
        ">LBNL_IF00016_AmpR_861bp_insertion_Staphylococcus.aureus_flank1\nTTTTTTTT\n"
    )
    got = labels(tmp_path, fasta)
    assert got["LBNL_IF00016_AmpR_861bp_insertion_Staphylococcus.aureus"] == 1
    assert got["LBNL_IF00016_AmpR_861bp_insertion_Staphylococcus.aureus_flank1"] == 0
    loaded = FelixSignaturesAdapter().load(make_zip(tmp_path, fasta), {})
    assert "group" in loaded.manifest_extra["leakage_warning"].lower()
