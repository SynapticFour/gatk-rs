//! 6R.116 coordinate-free: equal-length haplotype-to-reference CIGAR must run
//! Java padded SW, not an all-`M` shortcut.
//!
//! GATK 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//! `CigarUtils.calculateCigar` short-circuits **only** when `Arrays.equals(refSeq, altSeq)`.
//! A 2-mismatch shortcut was reverted because SW parameters can prefer indel over
//! substitutions. Equal-length `TTCATGA` vs `TATGTGA` therefore gets padded SW and
//! the EventMap-compatible `2D1M2I` CIGAR, not leftover-Match SNPs.
//!
//! Production `calculate_haplotype_cigar_sw` skips SW only on exact sequence
//! equality (6R.117). 6R.116 proved the Java `Arrays.equals` contract and that
//! the former equal-length `{len}M` shortcut leftover-Matched `A/G`.
//!
//! Does not pin genomic `92307327`. Sequence relationship only.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r116_equal_length_cigar_sw_contract
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
    let s = c.to_gatk_string();
    s.contains("2D1M2I")
}

/// Java `CigarUtils.calculateCigar` on equal-length coupled motif is indel, not all-`M`.
#[test]
fn forensic_6r116_java_padded_sw_encodes_coupled_indel_not_all_match() {
    let ref_bases = fill_window(PREFIX, REF_MOTIF, REF_LEN);
    let hap_bases = fill_window(PREFIX, HAP_MOTIF, REF_LEN);
    let sw = SwParameters::gatk_haplotype_to_reference();
    let java = calculate_haplotype_cigar_java_padded_sw(
        &ref_bases,
        &hap_bases,
        &sw,
        SwOverhangStrategy::SoftClip,
    )
    .expect("Java padded SW");
    assert!(
        cigar_has_2d1m2i(&java.cigar),
        "Java calculateCigar must encode 2D1M2I, got {}",
        java.cigar.to_gatk_string()
    );
    assert!(
        !cigar_is_all_m(&java.cigar),
        "Java must not all-M equal-length TATGTGA vs TTCATGA"
    );
    assert_eq!(java.cigar.reference_length(), REF_LEN);
}

/// Production assembly CIGAR matches Java padded SW (`2D1M2I`), not `{len}M`.
/// 6R.117 removed the equal-length shortcut 6R.116 identified.
#[test]
fn forensic_6r116_production_assembly_cigar_is_all_match_for_equal_length() {
    let ref_bases = fill_window(PREFIX, REF_MOTIF, REF_LEN);
    let hap_bases = fill_window(PREFIX, HAP_MOTIF, REF_LEN);
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
        cigar_has_2d1m2i(&prod.cigar),
        "production must encode 2D1M2I after 6R.117, got {}",
        prod.cigar.to_gatk_string()
    );
    assert!(
        !cigar_is_all_m(&prod.cigar),
        "equal-length shortcut must not write all-M, got {}",
        prod.cigar.to_gatk_string()
    );
    assert_eq!(
        prod.cigar.to_gatk_string(),
        java.cigar.to_gatk_string(),
        "production extract CIGAR must match Java calculateCigar"
    );
}

/// All-`M` leftover-Match EventMap is `A/G`; Java SW CIGAR EventMap is `A/ATG`.
#[test]
fn forensic_6r116_all_match_vs_java_sw_eventmap_at_insertion() {
    let pad: u64 = 1;
    let ref_bases = fill_window(PREFIX, REF_MOTIF, REF_LEN);
    let hap_bases = fill_window(PREFIX, HAP_MOTIF, REF_LEN);
    let sw = SwParameters::gatk_haplotype_to_reference();
    let mut ref_m = Cigar::new();
    ref_m.push(REF_LEN, CigarOperator::Match);
    let ref_hap = hap_with_cigar(ref_bases.clone(), ref_m, true, pad);
    let atg_start = pad + PREFIX as u64 + 3;

    let mut all_m = Cigar::new();
    all_m.push(REF_LEN, CigarOperator::Match);
    let leftover = hap_with_cigar(hap_bases.clone(), all_m, false, pad);
    let leftover_ev =
        variation_events_for_haplotype(&leftover, &ref_hap, &ref_bases, pad, 0, "synth");
    assert!(
        has_allele(&leftover_ev, atg_start, "A", "G"),
        "all-M leftover-Matches A/G"
    );
    assert!(!has_allele(&leftover_ev, atg_start, "A", "ATG"));

    let java = calculate_haplotype_cigar_java_padded_sw(
        &ref_bases,
        &hap_bases,
        &sw,
        SwOverhangStrategy::SoftClip,
    )
    .expect("Java padded SW");
    let sw_hap = hap_with_cigar(hap_bases, java.cigar, false, pad);
    let sw_ev = variation_events_for_haplotype(&sw_hap, &ref_hap, &ref_bases, pad, 0, "synth");
    assert!(
        has_allele(&sw_ev, atg_start, "A", "ATG"),
        "Java SW CIGAR must replay A/ATG, got {:?}",
        sw_ev
            .iter()
            .filter(|e| e.start_1based.get() == atg_start)
            .map(|e| format!("{}/{}", e.ref_allele, e.alt_allele))
            .collect::<Vec<_>>()
    );
    assert!(!has_allele(&sw_ev, atg_start, "A", "G"));
}
