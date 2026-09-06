//! 6R.104 coordinate-free: overlapping spanning deletion + SNP after `replaceSpanDels`.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//! `HaplotypeCallerGenotypingEngine.assignGenotypeLikelihoods`:
//!
//! ```text
//! getVariantContextsFromActiveHaplotypes(loc, haps, includeSpanning=true)
//!   → events whose EventMap interval overlaps loc
//! replaceSpanDels
//!   → start == loc kept; start < loc becomes (single-base REF, *)
//! makeMergedVariantContext / simpleMerge
//!   → one genotyping VC; * is a real allele when a spanning event existed
//! ```
//!
//! Live HOLDOUT_6R53 (`20:29455388 C/T`) does **not** exercise that merge: Java
//! EventMap at the loc is solely `C/T` (no overlapping indel). Rust's union
//! also carries a long deletion that starts earlier; after `replaceSpanDels`
//! that becomes `C/*` beside `C/T`. `merged_alleles_for_genotyping` returns
//! `None` for same-REF SNP+* (it only merges when a shorter REF remaps), so
//! production then genotypes `C/T` as an independent biallelic.
//!
//! That same-REF SNP+* geometry is **not** the cause of the holdout VCF extra:
//! Java already has `C/T` in EventMap and still does not emit it. Production
//! change: NONE. This file pins the merge-input contract so the next arrow does
//! not confuse EventMap presence with VCF emission.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r104_holdout_allele_set_contract
//! HOLDOUT_6R104=1 cargo test -p gatk-haplotypecaller --test holdout_6r104_holdout_allele_set -- --nocapture
//! ```

use gatk_haplotypecaller::event_map::{merged_alleles_for_genotyping, VariationEvent};
use gatk_haplotypecaller::hc_allele_mapping::{replace_span_del_events, SPAN_DEL_ALLELE};
use std::collections::HashSet;

fn ev(start: u64, r: &str, a: &str) -> VariationEvent {
    VariationEvent::from_alleles("chr", start, r, a)
}

/// Java `simpleMerge` allele list after `replaceSpanDels` (same-REF SNP + *).
fn java_simple_merge_same_ref(events: &[VariationEvent]) -> Vec<String> {
    let long_ref = events
        .iter()
        .map(|e| e.ref_allele.as_str())
        .max_by_key(|r| r.len())
        .unwrap_or("")
        .to_string();
    let mut alleles = Vec::new();
    let mut seen = HashSet::new();
    seen.insert(long_ref.clone());
    alleles.push(long_ref.clone());
    for e in events {
        if e.ref_allele == long_ref && e.alt_allele != long_ref && seen.insert(e.alt_allele.clone())
        {
            alleles.push(e.alt_allele.clone());
        }
    }
    alleles
}

/// Overlapping spanning del (earlier start) + SNP: `replaceSpanDels` yields same-REF SNP and *.
#[test]
fn forensic_6r104_overlapping_span_del_becomes_star_beside_snp() {
    let loc = 200u64;
    let spanning = ev(180, "AAAAAAAAAAAAAAAAAAAAA", "A");
    let snp = ev(200, "C", "T");
    assert!(spanning.start_1based.get() < loc);
    assert!(
        spanning.end_1based.get() >= loc,
        "fixture spanning del must overlap loc (end={})",
        spanning.end_1based.get()
    );
    let replaced = replace_span_del_events(&[spanning, snp], loc, 100, &vec![b'N'; 120]);
    assert_eq!(replaced.len(), 2);
    assert!(replaced.iter().any(|e| e.alt_allele == SPAN_DEL_ALLELE));
    assert!(replaced
        .iter()
        .any(|e| e.ref_allele == "C" && e.alt_allele == "T"));
    assert!(replaced.iter().all(|e| e.start_1based.get() == loc));
    assert!(replaced.iter().all(|e| e.ref_allele.len() == 1));
}

/// Java simpleMerge of that pair is `[C, *, T]` (order: REF, then encounter).
/// Rust `merged_alleles_for_genotyping` requires a shorter REF to remap, so it is
/// `None` — the production walk then genotypes the SNP independently.
#[test]
fn forensic_6r104_same_ref_snp_and_star_is_not_colocated_indel_merge() {
    let loc = 200u64;
    let star = ev(200, "C", SPAN_DEL_ALLELE);
    let snp = ev(200, "C", "T");
    let java = java_simple_merge_same_ref(&[star.clone(), snp.clone()]);
    assert_eq!(java[0], "C");
    assert!(
        java.iter().any(|a| a == SPAN_DEL_ALLELE),
        "Java simpleMerge keeps * when a spanning event was replaced: {java:?}"
    );
    assert!(
        java.iter().any(|a| a == "T"),
        "Java simpleMerge keeps SNP T"
    );
    assert_eq!(
        merged_alleles_for_genotyping(&[star, snp], loc),
        None,
        "same-REF SNP+* is not the 6R.61 shorter-REF colocated merge"
    );
}

/// When EventMap at loc is only the SNP (Java live at HOLDOUT_6R53), replaceSpanDels
/// is a no-op and simpleMerge is biallelic `[REF, SNP]`.
#[test]
fn forensic_6r104_snp_only_eventmap_stays_biallelic() {
    let loc = 200u64;
    let snp = ev(200, "C", "T");
    let replaced = replace_span_del_events(&[snp.clone()], loc, 100, &vec![b'N'; 120]);
    assert_eq!(replaced.len(), 1);
    assert_eq!(replaced[0].alt_allele, "T");
    assert!(!replaced.iter().any(|e| e.alt_allele == SPAN_DEL_ALLELE));
    assert_eq!(merged_alleles_for_genotyping(&replaced, loc), None);
    let java = java_simple_merge_same_ref(&replaced);
    assert_eq!(java, vec!["C".to_string(), "T".to_string()]);
}

/// EventMap SNP presence is not VCF emission: unused-ALT / emit still apply.
/// Coordinate-free: unused-ALT of a 0/0-like GT drops the unused SNP alt.
#[test]
fn forensic_6r104_eventmap_snp_is_not_automatic_vcf_emit() {
    let alts = vec!["T".to_string()];
    let gls = vec![0.0, -8.1, -3.6];
    let ad = vec![44, 4];
    let subset =
        gatk_haplotypecaller::subset_unused_alts_after_merged_genotyping(&alts, &[0, 0], &gls, &ad)
            .expect("subset");
    assert!(
        subset.alt_alleles.is_empty(),
        "hom-ref GT drops the EventMap SNP before VCF: {:?}",
        subset.alt_alleles
    );
}
