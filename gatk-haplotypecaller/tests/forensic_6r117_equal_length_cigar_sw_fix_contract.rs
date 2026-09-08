//! 6R.117 coordinate-free: equal-length haplotype CIGAR assignment matches
//! GATK 4.4 `CigarUtils.calculateCigar`.
//!
//! GATK 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`) short-circuits SW
//! **only** when `Arrays.equals(refSeq, altSeq)`. Non-identical equal-length
//! sequences take padded SW. A coupled `TTCATGA` vs `TATGTGA` window therefore
//! gets the EventMap-compatible `2D1M2I` CIGAR, not leftover-Match `{len}M`.
//!
//! Does not pin genomic `92307327`. Sequence relationship only.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r117_equal_length_cigar_sw_fix_contract
//! ```

use gatk_haplotypecaller::cigar::{Cigar, CigarOperator};
use gatk_haplotypecaller::event_map::{variation_events_for_haplotype, VariationEvent};
use gatk_haplotypecaller::genome_loc::GenomeLoc;
use gatk_haplotypecaller::haplotype::Haplotype;
use gatk_haplotypecaller::haplotype_cigar::{
    calculate_haplotype_cigar_for_assembly_with_offset, calculate_haplotype_cigar_java_padded_sw,
};
use gatk_haplotypecaller::smith_waterman::{SwOverhangStrategy, SwParameters};

const REF_MOTIF: &[u8] = b"TTCATGA";
const HAP_MOTIF: &[u8] = b"TATGTGA";
const PREFIX: usize = 79;
const REF_LEN: usize = 191;

fn fill_window(prefix: usize, motif: &[u8], len: usize) -> Vec<u8> {
    let mut v = vec![b'C'; len];
    assert!(prefix + motif.len() <= len);
    v[prefix..prefix + motif.len()].copy_from_slice(motif);
    v
}

fn hap_with_cigar(bases: Vec<u8>, cigar: Cigar, is_ref: bool, pad: u64) -> Haplotype {
    let len = bases.len() as u64;
    let mut h = Haplotype::new(bases, is_ref);
    h.cigar = Some(cigar);
    h.alignment_start_hap_wrt_ref = 0;
    h.genome_loc = Some(GenomeLoc::new(
        pad,
        pad.saturating_add(len).saturating_sub(1),
    ));
    h
}

fn has_allele(events: &[VariationEvent], start: u64, r: &str, a: &str) -> bool {
    events
        .iter()
        .any(|e| e.start_1based.get() == start && e.ref_allele == r && e.alt_allele == a)
}

fn cigar_is_all_m(c: &Cigar) -> bool {
    c.elements.len() == 1 && matches!(c.elements[0].operator, CigarOperator::Match)
}

fn cigar_has_2d1m2i(c: &Cigar) -> bool {
    c.to_gatk_string().contains("2D1M2I")
}

/// A. Exact equality remains exact-match / all-`M`.
#[test]
fn forensic_6r117_identical_equal_length_remains_all_match() {
    let bases = fill_window(PREFIX, REF_MOTIF, REF_LEN);
    let sw = SwParameters::gatk_haplotype_to_reference();
    let prod = calculate_haplotype_cigar_for_assembly_with_offset(&bases, &bases, REF_LEN, &sw)
        .expect("identical assembly CIGAR");
    assert!(
        cigar_is_all_m(&prod.cigar),
        "identical sequences must remain all-M, got {}",
        prod.cigar.to_gatk_string()
    );
    assert_eq!(prod.cigar.to_gatk_string(), format!("{REF_LEN}M"));
}

/// B+C. Non-identical equal-length haplotypes take padded SW, not `{len}M`.
#[test]
fn forensic_6r117_non_identical_equal_length_uses_padded_sw() {
    let ref_bases = fill_window(PREFIX, REF_MOTIF, REF_LEN);
    let hap_bases = fill_window(PREFIX, HAP_MOTIF, REF_LEN);
    assert_eq!(ref_bases.len(), hap_bases.len());
    assert_ne!(ref_bases, hap_bases);
    let sw = SwParameters::gatk_haplotype_to_reference();
    let prod =
        calculate_haplotype_cigar_for_assembly_with_offset(&ref_bases, &hap_bases, REF_LEN, &sw)
            .expect("production assembly CIGAR");
    let java = calculate_haplotype_cigar_java_padded_sw(
        &ref_bases,
        &hap_bases,
        &sw,
        SwOverhangStrategy::SoftClip,
    )
    .expect("Java padded SW");
    assert!(
        !cigar_is_all_m(&prod.cigar),
        "non-identical equal-length must not take {{len}}M shortcut, got {}",
        prod.cigar.to_gatk_string()
    );
    assert!(
        cigar_has_2d1m2i(&prod.cigar),
        "Java-compatible indel CIGAR must contain 2D1M2I, got {}",
        prod.cigar.to_gatk_string()
    );
    assert_eq!(
        prod.cigar.to_gatk_string(),
        java.cigar.to_gatk_string(),
        "production must match Java calculateCigar"
    );
    assert_eq!(prod.cigar.reference_length(), REF_LEN);
}

/// D. Production CIGAR EventMap is the indel, not leftover-Match SNP.
#[test]
fn forensic_6r117_production_cigar_eventmap_is_indel_not_leftover_snp() {
    let pad: u64 = 1;
    let ref_bases = fill_window(PREFIX, REF_MOTIF, REF_LEN);
    let hap_bases = fill_window(PREFIX, HAP_MOTIF, REF_LEN);
    let sw = SwParameters::gatk_haplotype_to_reference();
    let mut ref_m = Cigar::new();
    ref_m.push(REF_LEN, CigarOperator::Match);
    let ref_hap = hap_with_cigar(ref_bases.clone(), ref_m, true, pad);
    let atg_start = pad + PREFIX as u64 + 3;

    let prod =
        calculate_haplotype_cigar_for_assembly_with_offset(&ref_bases, &hap_bases, REF_LEN, &sw)
            .expect("production assembly CIGAR");
    let sw_hap = hap_with_cigar(hap_bases, prod.cigar, false, pad);
    let sw_ev = variation_events_for_haplotype(&sw_hap, &ref_hap, &ref_bases, pad, 0, "synth");
    assert!(
        has_allele(&sw_ev, atg_start, "A", "ATG"),
        "production SW CIGAR must replay A/ATG, got {:?}",
        sw_ev
            .iter()
            .filter(|e| e.start_1based.get() == atg_start)
            .map(|e| format!("{}/{}", e.ref_allele, e.alt_allele))
            .collect::<Vec<_>>()
    );
    assert!(
        !has_allele(&sw_ev, atg_start, "A", "G"),
        "production SW CIGAR must not leftover-Match A/G"
    );
}
