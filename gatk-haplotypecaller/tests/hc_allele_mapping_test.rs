//! CIGAR-aware allele support at padded-ref loci (N2 / cluster genotyping).

use gatk_haplotypecaller::bio_ids::HaplotypeIndex;
use gatk_haplotypecaller::cigar::{Cigar, CigarOperator};
use gatk_haplotypecaller::event_map::VariationEvent;
use gatk_haplotypecaller::genome_loc::{GenomeLoc, GenomePosition};
use gatk_haplotypecaller::haplotype::Haplotype;
use gatk_haplotypecaller::hc_allele_mapping::{
    create_allele_mapper, hap_base_at_ref_locus, haplotype_supports_allele_at,
};

#[test]
fn hap_base_at_ref_locus_skips_deletion() {
    let pad = 100u64;
    let mut h = Haplotype::new(b"ACGT", false);
    let mut cigar = Cigar::new();
    cigar.push(2, CigarOperator::Match);
    cigar.push(2, CigarOperator::Deletion);
    cigar.push(2, CigarOperator::Match);
    h.cigar = Some(cigar);
    h.alignment_start_hap_wrt_ref = 0;
    assert_eq!(hap_base_at_ref_locus(&h, pad, 100), Some(b'A'));
    assert_eq!(hap_base_at_ref_locus(&h, pad, 102), None);
    assert_eq!(hap_base_at_ref_locus(&h, pad, 104), Some(b'G'));
}

#[test]
fn snp_support_uses_cigar_not_linear_offset() {
    let pad = 0u64;
    let loc = 37u64;
    let mut ref_hap = Haplotype::new(vec![b'T'; 50], true);
    let mut ref_cigar = Cigar::new();
    ref_cigar.push(50, CigarOperator::Match);
    ref_hap.cigar = Some(ref_cigar);

    let mut alt = Haplotype::new(vec![b'T'; 48], false);
    let mut alt_cigar = Cigar::new();
    alt_cigar.push(35, CigarOperator::Match);
    alt_cigar.push(2, CigarOperator::Deletion);
    alt_cigar.push(13, CigarOperator::Match);
    alt.bases[35] = b'C';
    alt.cigar = Some(alt_cigar);

    assert!(haplotype_supports_allele_at(
        &alt, &ref_hap, loc, pad, "T", "C"
    ));
    assert!(!haplotype_supports_allele_at(
        &ref_hap, &ref_hap, loc, pad, "T", "C"
    ));
}

#[test]
fn indel_support_uses_event_map_not_linear_offset() {
    let pad = 92307249u64;
    let loc = 92307324u64;
    let mut ref_bytes = vec![b'A'; 120];
    ref_bytes[75] = b'T';
    ref_bytes[76] = b'T';
    ref_bytes[77] = b'C';
    let mut ref_hap = Haplotype::new(ref_bytes.clone(), true);
    let mut ref_cigar = Cigar::new();
    ref_cigar.push(ref_bytes.len(), CigarOperator::Match);
    ref_hap.cigar = Some(ref_cigar);

    let mut alt_bases = ref_bytes[0..76].to_vec();
    alt_bases.extend_from_slice(&ref_bytes[79..]);
    let mut alt = Haplotype::new(alt_bases, false);
    let mut alt_cigar = Cigar::new();
    alt_cigar.push(76, CigarOperator::Match);
    alt_cigar.push(2, CigarOperator::Deletion);
    alt_cigar.push(alt.bases.len().saturating_sub(76), CigarOperator::Match);
    alt.cigar = Some(alt_cigar);

    assert!(haplotype_supports_allele_at(
        &alt, &ref_hap, loc, pad, "TTC", "T"
    ));
    assert!(!haplotype_supports_allele_at(
        &ref_hap, &ref_hap, loc, pad, "TTC", "T"
    ));
}

#[test]
fn coupled_atg_support_survives_preceding_deletion_offset() {
    let pad = 92307319u64;
    let ttc_loc = 92307324u64;
    let atg_loc = 92307327u64;
    let mut ref_bytes = vec![b'A'; 24];
    ref_bytes[(ttc_loc - pad) as usize..=(atg_loc - pad) as usize].copy_from_slice(b"TTCA");
    let mut ref_hap = Haplotype::new(ref_bytes.clone(), true);
    let mut ref_cigar = Cigar::new();
    ref_cigar.push(ref_bytes.len(), CigarOperator::Match);
    ref_hap.cigar = Some(ref_cigar);

    let mut alt_bases = ref_bytes[..(ttc_loc - pad + 1) as usize].to_vec();
    alt_bases.push(b'A');
    alt_bases.extend_from_slice(b"TG");
    alt_bases.extend_from_slice(&ref_bytes[(atg_loc - pad + 1) as usize..]);
    let mut alt = Haplotype::new(alt_bases, false);
    let mut alt_cigar = Cigar::new();
    alt_cigar.push((ttc_loc - pad + 1) as usize, CigarOperator::Match);
    alt_cigar.push(2, CigarOperator::Deletion);
    alt_cigar.push(1, CigarOperator::Match);
    alt_cigar.push(2, CigarOperator::Insertion);
    alt_cigar.push(
        ref_bytes.len().saturating_sub((atg_loc - pad + 1) as usize),
        CigarOperator::Match,
    );
    alt.cigar = Some(alt_cigar);

    assert!(haplotype_supports_allele_at(
        &alt, &ref_hap, atg_loc, pad, "A", "ATG",
    ));
}

#[test]
fn create_allele_mapper_snp_pools_by_base_without_eventmap() {
    let pad = 92305524u64;
    let loc = 92305634u64;
    let ref_bytes = vec![b'G'; 200];
    let mut ref_hap = Haplotype::new(ref_bytes.clone(), true);
    let mut ref_cigar = gatk_haplotypecaller::cigar::Cigar::new();
    ref_cigar.push(
        ref_bytes.len(),
        gatk_haplotypecaller::cigar::CigarOperator::Match,
    );
    ref_hap.cigar = Some(ref_cigar);

    let mut alt_t = Haplotype::new(ref_bytes.clone(), false);
    alt_t.bases[(loc - pad) as usize] = b'T';
    let mut alt_cigar = gatk_haplotypecaller::cigar::Cigar::new();
    alt_cigar.push(
        ref_bytes.len(),
        gatk_haplotypecaller::cigar::CigarOperator::Match,
    );
    alt_t.cigar = Some(alt_cigar);

    let haps = vec![ref_hap, alt_t];
    let merged = VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(loc),
        end_1based: GenomePosition::new_1based(loc),
        ref_allele: "G".into(),
        alt_allele: "T".into(),
    };
    let m = create_allele_mapper(&merged, loc, &haps, pad, &ref_bytes, 1, true);
    assert_eq!(m.ref_haplotype_indices, vec![HaplotypeIndex::new(0)]);
    assert_eq!(m.alt_haplotype_indices, vec![HaplotypeIndex::new(1)]);
}

#[test]
fn allele_mapper_ref_and_alt_pools_are_disjoint() {
    let pad = 100u64;
    let loc = 110u64;
    let mut ref_bytes = vec![b'A'; 40];
    ref_bytes[(loc - pad) as usize] = b'G';
    let mut ref_hap = Haplotype::new(ref_bytes.clone(), true);
    let mut ref_cigar = Cigar::new();
    ref_cigar.push(ref_bytes.len(), CigarOperator::Match);
    ref_hap.cigar = Some(ref_cigar);

    let mut alt = Haplotype::new(ref_bytes.clone(), false);
    alt.bases[(loc - pad) as usize] = b'T';
    let mut alt_cigar = Cigar::new();
    alt_cigar.push(ref_bytes.len(), CigarOperator::Match);
    alt.cigar = Some(alt_cigar);

    // Ambiguous third hap matches both REF and ALT paths in linear offset; mapper must not
    // keep it in both pools (Java createAlleleMapper disjoint membership).
    let mut dual = Haplotype::new(ref_bytes.clone(), false);
    dual.bases[(loc - pad) as usize] = b'T';
    let mut dual_cigar = Cigar::new();
    dual_cigar.push(ref_bytes.len(), CigarOperator::Match);
    dual.cigar = Some(dual_cigar);

    let haps = vec![ref_hap, alt, dual];
    let merged = VariationEvent {
        contig: "20".into(),
        start_1based: GenomePosition::new_1based(loc),
        end_1based: GenomePosition::new_1based(loc),
        ref_allele: "G".into(),
        alt_allele: "T".into(),
    };
    let m = create_allele_mapper(&merged, loc, &haps, pad, &ref_bytes, 1, true);
    let ref_set: std::collections::HashSet<_> = m.ref_haplotype_indices.iter().copied().collect();
    let alt_set: std::collections::HashSet<_> = m.alt_haplotype_indices.iter().copied().collect();
    assert!(
        ref_set.is_disjoint(&alt_set),
        "ref={ref_set:?} alt={alt_set:?}"
    );
    assert!(ref_set.contains(&HaplotypeIndex::new(0)));
    assert!(alt_set.contains(&HaplotypeIndex::new(1)));
}

/// Trim-window genotyping pad vs full-pad `alignment_start_hap_wrt_ref` (Class-A2 root cause).
/// Hap carries REF at the SNP but ALT at a nearby site — must not be read as ALT at the SNP.
#[test]
fn hap_base_reconciles_trim_pad_with_full_pad_alignment_start() {
    let full_pad = 10008590u64;
    let trim_pad = 10009207u64;
    let snp_loc = 10009227u64;
    let align0 = (trim_pad - full_pad) as usize; // 617
    let off = (snp_loc - trim_pad) as usize; // 20

    let mut bases = vec![b'A'; 80];
    bases[off] = b'A'; // REF at SNP
    bases[off + 19] = b'G'; // nearby alt (10009246)

    let mut hap = Haplotype::new(bases, false);
    let mut cigar = Cigar::new();
    cigar.push(80, CigarOperator::Match);
    hap.cigar = Some(cigar);
    hap.alignment_start_hap_wrt_ref = align0;
    hap.genome_loc = Some(GenomeLoc::new(trim_pad, trim_pad + 79));

    // Caller passes trim pad (production genotyping) — must still see REF A.
    assert_eq!(hap_base_at_ref_locus(&hap, trim_pad, snp_loc), Some(b'A'));
    // Caller passes full pad — same answer.
    assert_eq!(hap_base_at_ref_locus(&hap, full_pad, snp_loc), Some(b'A'));

    let ref_bases = vec![b'A'; 80];
    let mut ref_hap = Haplotype::new(ref_bases.clone(), true);
    let mut ref_cigar = Cigar::new();
    ref_cigar.push(80, CigarOperator::Match);
    ref_hap.cigar = Some(ref_cigar);
    ref_hap.alignment_start_hap_wrt_ref = align0;
    ref_hap.genome_loc = Some(GenomeLoc::new(trim_pad, trim_pad + 79));

    let mut alt_bases = vec![b'A'; 80];
    alt_bases[off] = b'G';
    let mut alt = Haplotype::new(alt_bases, false);
    let mut alt_cigar = Cigar::new();
    alt_cigar.push(80, CigarOperator::Match);
    alt.cigar = Some(alt_cigar);
    alt.alignment_start_hap_wrt_ref = align0;
    alt.genome_loc = Some(GenomeLoc::new(trim_pad, trim_pad + 79));

    let merged = VariationEvent {
        contig: "20".into(),
        start_1based: GenomePosition::new_1based(snp_loc),
        end_1based: GenomePosition::new_1based(snp_loc),
        ref_allele: "A".into(),
        alt_allele: "G".into(),
    };
    let m = create_allele_mapper(
        &merged,
        snp_loc,
        &[ref_hap, hap, alt],
        trim_pad,
        &ref_bases,
        1,
        false,
    );
    assert!(
        m.ref_haplotype_indices.contains(&HaplotypeIndex::new(1)),
        "nearby-only alt hap must stay in REF pool at SNP; ref={:?} alt={:?}",
        m.ref_haplotype_indices,
        m.alt_haplotype_indices
    );
    assert!(
        m.alt_haplotype_indices.contains(&HaplotypeIndex::new(2)),
        "true SNP alt hap must be in ALT pool; alt={:?}",
        m.alt_haplotype_indices
    );
    assert!(
        !m.alt_haplotype_indices.contains(&HaplotypeIndex::new(1)),
        "nearby-only alt hap must not enter ALT pool; alt={:?}",
        m.alt_haplotype_indices
    );
}

/// 6R.66: Java `createAlleleMapper` with no overlapping EventMap events at `loc`
/// assigns the hap to REF. Rust must not pad-slice an apparent indel (`TG/T`).
fn six_r65_empty_overlap_indel_fixture() -> (VariationEvent, Vec<Haplotype>, Vec<u8>, u64, u64) {
    let pad = 100u64;
    let loc = 110u64;
    let mut ref_bytes = vec![b'A'; 40];
    ref_bytes[(loc - pad) as usize] = b'T';
    ref_bytes[(loc - pad + 1) as usize] = b'G';
    let mut ref_hap = Haplotype::new(ref_bytes.clone(), true);
    let mut ref_cigar = Cigar::new();
    ref_cigar.push(ref_bytes.len(), CigarOperator::Match);
    ref_hap.cigar = Some(ref_cigar);

    // All-M hap: match at `loc` (T), SNP G>A at loc+1. EventMap overlaps loc+1, not loc.
    // Pad slice at loc is `TA…`, which matches merged ALT `T` vs REF `TG`.
    let mut alt = Haplotype::new(ref_bytes.clone(), false);
    alt.bases[(loc - pad + 1) as usize] = b'A';
    let mut alt_cigar = Cigar::new();
    alt_cigar.push(ref_bytes.len(), CigarOperator::Match);
    alt.cigar = Some(alt_cigar);

    let merged = VariationEvent {
        contig: "20".into(),
        start_1based: GenomePosition::new_1based(loc),
        end_1based: GenomePosition::new_1based(loc + 1),
        ref_allele: "TG".into(),
        alt_allele: "T".into(),
    };
    (merged, vec![ref_hap, alt], ref_bytes, pad, loc)
}

#[test]
fn empty_overlapping_eventmap_does_not_assign_indel_via_pad_slice() {
    let (merged, haps, ref_bytes, pad, loc) = six_r65_empty_overlap_indel_fixture();
    let m = create_allele_mapper(&merged, loc, &haps, pad, &ref_bytes, 0, true);
    assert!(
        m.ref_haplotype_indices.contains(&HaplotypeIndex::new(1)),
        "Java createAlleleMapper: empty overlapping EventMap → REF; ref={:?} alt={:?}",
        m.ref_haplotype_indices,
        m.alt_haplotype_indices
    );
    assert!(
        !m.alt_haplotype_indices.contains(&HaplotypeIndex::new(1)),
        "6R.66: no pad-slice indel fallback for TG/T; alt={:?}",
        m.alt_haplotype_indices
    );
}
