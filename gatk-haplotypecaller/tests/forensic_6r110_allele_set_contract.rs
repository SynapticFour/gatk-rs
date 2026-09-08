//! 6R.110 coordinate-free: Java `simpleMerge` of same-REF multi-indel alts.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//! `AssemblyBasedCallerUtils.makeMergedVariantContext` / `simpleMerge` builds one
//! genotyping VC from every EventMap event at the loc. Same-REF insertions
//! `G/GTT` + `G/GTTT` become `[G, GTT, GTTT]`. `calculateOutputAlleleSubset`
//! then drops an ALT that fails the AFC emit threshold; unused-ALT subset after
//! GT `0/2` is the same remaining-allele class.
//!
//! Rust previously returned `None` from [`merged_alleles_for_genotyping`] for
//! that geometry (biallelic walk + emit-time coalesce), which is how the extra
//! insertion reached VCF. 6R.104 SNP+`*` (one non-star alt) stays off this path.
//! Nested-STR dels with unequal REF lengths remain the remapped biallelic walk.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r110_allele_set_contract
//! HOLDOUT_6R110=1 cargo test -p gatk-haplotypecaller --test holdout_6r110_allele_set -- --nocapture
//! ```

use gatk_haplotypecaller::event_map::{
    is_colocated_snp_indel_merged_site, is_same_ref_multi_indel_merged_site,
    merged_alleles_for_genotyping, merged_site_uses_joint_gls, VariationEvent,
};
use gatk_haplotypecaller::hc_allele_mapping::SPAN_DEL_ALLELE;
use gatk_haplotypecaller::subset_unused_alts_after_merged_genotyping;

fn ev(start: u64, r: &str, a: &str) -> VariationEvent {
    VariationEvent::from_alleles("chr", start, r, a)
}

#[test]
fn forensic_6r110_same_ref_insertions_are_java_simple_merge() {
    let loc = 50u64;
    let a = ev(loc, "G", "GTT");
    let b = ev(loc, "G", "GTTT");
    let events = [a, b];
    let (long_ref, alts) = merged_alleles_for_genotyping(&events, loc).expect("simpleMerge");
    assert_eq!(long_ref, "G");
    assert_eq!(alts, vec!["GTT".to_string(), "GTTT".to_string()]);
    assert!(!is_colocated_snp_indel_merged_site(&long_ref, &alts));
    assert!(is_same_ref_multi_indel_merged_site(
        &events, loc, &long_ref, &alts
    ));
    assert!(merged_site_uses_joint_gls(&events, loc, &long_ref, &alts));
}

#[test]
fn forensic_6r110_snp_star_stays_biallelic_walk() {
    let loc = 200u64;
    let snp = ev(loc, "C", "T");
    let star = ev(loc, "C", SPAN_DEL_ALLELE);
    assert_eq!(
        merged_alleles_for_genotyping(&[snp, star], loc),
        None,
        "6R.104: same-REF SNP+* is one non-star alt, not this merge"
    );
}

#[test]
fn forensic_6r110_nested_str_unequal_ref_is_not_this_class() {
    let loc = 100u64;
    let short = ev(loc, "ATGTGTGTG", "A");
    let long = ev(loc, "ATGTGTGTGTGTGTGTGTG", "A");
    let events = [short, long];
    let (long_ref, alts) = merged_alleles_for_genotyping(&events, loc).expect("remap");
    assert!(!is_colocated_snp_indel_merged_site(&long_ref, &alts));
    assert!(
        !is_same_ref_multi_indel_merged_site(&events, loc, &long_ref, &alts),
        "unequal REF lengths are not the 6R.110 same-REF insertion class"
    );
}

#[test]
fn forensic_6r110_gt_0_2_drops_unused_insertion() {
    let alts = vec!["GTT".to_string(), "GTTT".to_string()];
    let gls = vec![0.0, -8.0, -40.0, -1.0, -20.0, -3.0];
    let ad = vec![10, 4, 8];
    let subset =
        subset_unused_alts_after_merged_genotyping(&alts, &[0, 2], &gls, &ad).expect("subset");
    assert_eq!(
        subset.alt_alleles,
        vec!["GTTT".to_string()],
        "GT 0/2 keeps allele 2; emit-time coalesce of two hets is not this contract"
    );
}
