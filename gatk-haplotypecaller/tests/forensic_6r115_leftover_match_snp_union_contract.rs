//! 6R.115 coordinate-free: leftover-Match SNP and coupled insertion are
//! independent EventMap alleles across haplotypes.
//!
//! GATK 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//! `EventMap.processCigarForInitialEvents` is per haplotype. An all-`M` CIGAR on
//! a coupled `TTCATGA`→`TATGTGA` sequence leftover-Matches `A/G` at the insertion
//! locus. The EventMap-compatible `2D1M2I` CIGAR on the same bytes replays `A/ATG`.
//! `prefer_indel_over_colocated_snps` runs per haplotype, so it cannot drop the
//! SNP on a haplotype that has no indel at that start.
//!
//! `getVariantContextsFromActiveHaplotypes` unions unique `(start, alleles)`.
//! When both encodings are present, the union is `{A/ATG, A/G}` and
//! `makeMergedVariantContext` / [`merged_alleles_for_genotyping`] is colocated
//! SNP+indel `[A, ATG, G]` (6R.61 class), not same-REF multi-indel (6R.110).
//! A population with only REF + the compatible-CIGAR hap has `{A/ATG}` only.
//!
//! Does not pin genomic `92307327` as the contract. Sequence relationship only.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r115_leftover_match_snp_union_contract
//! ```

use gatk_haplotypecaller::cigar::{Cigar, CigarOperator};
use gatk_haplotypecaller::event_map::{
    collect_variation_events, is_colocated_snp_indel_merged_site,
    is_same_ref_multi_indel_merged_site, merged_alleles_for_genotyping,
    variation_events_for_haplotype, VariationEvent,
};
use gatk_haplotypecaller::genome_loc::GenomeLoc;
use gatk_haplotypecaller::haplotype::Haplotype;

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

fn cigar_all_match(len: usize) -> Cigar {
    let mut c = Cigar::new();
    c.push(len, CigarOperator::Match);
    c
}

fn cigar_coupled(lead: usize, ref_len: usize) -> Cigar {
    let tail = ref_len.saturating_sub(lead.saturating_add(3));
    let mut c = Cigar::new();
    if lead > 0 {
        c.push(lead, CigarOperator::Match);
    }
    c.push(2, CigarOperator::Deletion);
    c.push(1, CigarOperator::Match);
    c.push(2, CigarOperator::Insertion);
    if tail > 0 {
        c.push(tail, CigarOperator::Match);
    }
    c
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

fn alleles_at(events: &[VariationEvent], start: u64) -> Vec<String> {
    events
        .iter()
        .filter(|e| e.start_1based.get() == start)
        .map(|e| format!("{}/{}", e.ref_allele, e.alt_allele))
        .collect()
}

/// Compatible CIGAR + REF only: EventMap at the insertion loc is `A/ATG`, not leftover `A/G`.
#[test]
fn forensic_6r115_compatible_population_eventmap_is_insertion_only() {
    let pad: u64 = 1;
    let ref_bases = fill_window(PREFIX, REF_MOTIF, REF_LEN);
    let hap_bases = fill_window(PREFIX, HAP_MOTIF, REF_LEN);
    let ref_hap = hap_with_cigar(ref_bases.clone(), cigar_all_match(REF_LEN), true, pad);
    let alt = hap_with_cigar(hap_bases, cigar_coupled(PREFIX + 1, REF_LEN), false, pad);
    let haps = vec![ref_hap, alt];
    let union = collect_variation_events(&haps, &ref_bases, pad, "synth", 0);
    let atg_start = pad + PREFIX as u64 + 3;
    assert_eq!(
        alleles_at(&union, atg_start),
        vec!["A/ATG".to_string()],
        "REF + EventMap-compatible coupled hap must not inject leftover A/G"
    );
    assert!(merged_alleles_for_genotyping(&union, atg_start).is_none());
}

/// All-M encoding of the same alt bytes leftover-Matches `A/G` at the insertion loc.
#[test]
fn forensic_6r115_all_match_encoding_replays_leftover_snp() {
    let pad: u64 = 1;
    let ref_bases = fill_window(PREFIX, REF_MOTIF, REF_LEN);
    let hap_bases = fill_window(PREFIX, HAP_MOTIF, REF_LEN);
    let ref_hap = hap_with_cigar(ref_bases.clone(), cigar_all_match(REF_LEN), true, pad);
    let alt = hap_with_cigar(hap_bases, cigar_all_match(REF_LEN), false, pad);
    let ev = variation_events_for_haplotype(&alt, &ref_hap, &ref_bases, pad, 0, "synth");
    let atg_start = pad + PREFIX as u64 + 3;
    assert!(
        has_allele(&ev, atg_start, "A", "G"),
        "all-M TATGTGA vs TTCATGA leftover-Matches A/G"
    );
    assert!(
        !has_allele(&ev, atg_start, "A", "ATG"),
        "all-M encoding must not reconstruct A/ATG"
    );
}

/// Union of both encodings is `{A/ATG, A/G}` and Java `simpleMerge` of SNP+indel.
#[test]
fn forensic_6r115_union_of_both_encodings_is_colocated_snp_indel_merge() {
    let pad: u64 = 1;
    let ref_bases = fill_window(PREFIX, REF_MOTIF, REF_LEN);
    let hap_bases = fill_window(PREFIX, HAP_MOTIF, REF_LEN);
    let ref_hap = hap_with_cigar(ref_bases.clone(), cigar_all_match(REF_LEN), true, pad);
    let coupled = hap_with_cigar(
        hap_bases.clone(),
        cigar_coupled(PREFIX + 1, REF_LEN),
        false,
        pad,
    );
    let leftover = hap_with_cigar(hap_bases, cigar_all_match(REF_LEN), false, pad);
    let haps = vec![ref_hap, coupled, leftover];
    let union = collect_variation_events(&haps, &ref_bases, pad, "synth", 0);
    let atg_start = pad + PREFIX as u64 + 3;
    let keys = alleles_at(&union, atg_start);
    assert!(
        keys.contains(&"A/ATG".to_string()) && keys.contains(&"A/G".to_string()),
        "per-hap EventMap union retains both encodings, got {keys:?}"
    );
    let (long_ref, alts) =
        merged_alleles_for_genotyping(&union, atg_start).expect("SNP+indel simpleMerge");
    assert_eq!(long_ref, "A");
    assert_eq!(alts, vec!["ATG".to_string(), "G".to_string()]);
    assert!(is_colocated_snp_indel_merged_site(&long_ref, &alts));
    assert!(!is_same_ref_multi_indel_merged_site(
        &union, atg_start, &long_ref, &alts
    ));
}
