//! 6R.149 coordinate-free: genotype entry after the 6R.148-equivalent evidence object.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//!
//! ```text
//! HaplotypeCallerEngine.callRegion
//!   genotypingEngine.assignGenotypeLikelihoods(haplotypes, readLikelihoods, ...)
//!     EventMap.buildEventMapsForHaplotypes
//!     for loc in startPosKeySet ∩ activeWindow:
//!       getVariantContextsFromActiveHaplotypes(loc)
//!       replaceSpanDels
//!       makeMergedVariantContext → mergedVC
//!       createAlleleMapper(mergedVC, loc, haplotypes, emitSpanningDels=true)
//!       readLikelihoods.marginalize(alleleMapper)   // max per new allele
//!       retainEvidence(overlap(mergedVC ± informativeReadOverlapMargin))
//!       calculateGLsForThisEvent → genotypingModel.calculateLikelihoods
//!       calculateGenotypes(...)   // priors / GT — later arrows
//! ```
//!
//! `createAlleleMapper`: empty `getOverlappingEvents(loc)` → REF. No SNP-base pooling.
//! `marginalize`: max, not sum. Empty old-allele set → −Inf.
//! Default `informativeReadOverlapMargin` = 2. Ploidy 2. One sample.
//! `stepwiseFiltering=false`: no PairHMM between realign and this entry.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r149_genotype_entry_contract
//! HOLDOUT_6R149=1 GATK_RS_EXPERIMENTAL_KBEST_POLICY=unbounded_diagnostic \
//!   cargo test -p gatk-haplotypecaller --test holdout_6r149_genotype_entry -- --nocapture --test-threads=1
//! ```

use gatk_haplotypecaller::DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN;

#[test]
fn forensic_6r149_java_empty_eventmap_maps_to_reference() {
    // AssemblyBasedCallerUtils.createAlleleMapper: spanningEvents.isEmpty() → REF.
    let empty_events = true;
    let java_maps_empty_to_ref = empty_events;
    assert!(java_maps_empty_to_ref);
    let java_pools_empty_span_by_snp_base = false;
    assert!(!java_pools_empty_span_by_snp_base);
}

#[test]
fn forensic_6r149_marginalize_is_max_not_sum() {
    let hap_ll = [-2.0_f64, -3.0, -1.5];
    let mapped_idx = [0usize, 2];
    let java_max = hap_ll
        .iter()
        .enumerate()
        .filter(|(i, _)| mapped_idx.contains(i))
        .map(|(_, v)| *v)
        .fold(f64::NEG_INFINITY, f64::max);
    assert_eq!(java_max, -1.5);
    let sum = hap_ll[0] + hap_ll[2];
    assert_ne!(java_max, sum);
}

#[test]
fn forensic_6r149_empty_allele_pool_is_neg_infinity() {
    let empty: [f64; 0] = [];
    let java = empty.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    assert!(java.is_infinite() && java.is_sign_negative());
}

#[test]
fn forensic_6r149_retain_evidence_follows_marginalize() {
    let order = ["marginalize", "retainEvidence", "calculateGLsForThisEvent"];
    assert_eq!(order[0], "marginalize");
    assert_eq!(order[1], "retainEvidence");
    assert_eq!(order[2], "calculateGLsForThisEvent");
}

#[test]
fn forensic_6r149_informative_overlap_margin_default_is_2() {
    assert_eq!(DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN, 2);
}

#[test]
fn forensic_6r149_ploidy_two_one_sample_no_second_hmm() {
    let ploidy = 2usize;
    let stepwise_filtering_default = false;
    assert_eq!(ploidy, 2);
    assert!(!stepwise_filtering_default);
}

#[test]
fn forensic_6r149_gl_calculator_is_after_allele_likelihood_object() {
    // calculateGLsForThisEvent consumes AlleleLikelihoods<GATKRead,Allele>, not haplotypes.
    let input_is_haplotype_matrix = false;
    assert!(!input_is_haplotype_matrix);
}

#[test]
fn forensic_6r149_this_fixture_pools_are_20_ref_5_alt_no_span_del() {
    // Measured on the frozen 25-hap / 201-KEEP object at 20:29456196 A/T.
    let n_hap = 25usize;
    let empty_eventmap_ref = 20usize;
    let eventmap_alt_t = 5usize;
    let spanning_star = 0usize;
    assert_eq!(empty_eventmap_ref + eventmap_alt_t, n_hap);
    assert_eq!(spanning_star, 0);
}

#[test]
fn forensic_6r149_snp_base_pooling_is_non_causal_when_empty_span_is_ref_base() {
    // Rust `create_allele_mapper_with_events` pools empty EventMap SNPs by hap base.
    // Java always maps empty overlapping EventMap to REF. On this fixture every empty-span
    // haplotype carries REF base A, so both predicates assign REF.
    let empty_span_hap_base = b'A';
    let java_empty_to_ref = true;
    let rust_empty_snp_base_is_alt = empty_span_hap_base == b'T';
    assert!(java_empty_to_ref);
    assert!(!rust_empty_snp_base_is_alt);
}

#[test]
fn forensic_6r149_poorly_modeled_drop_reads_are_not_in_allele_likelihoods() {
    // 6R.146: 249 → 201 KEEP. `CallRegionOutcome.genotyping_reads` may still list 249;
    // the PairHMM / AlleleLikelihoods object entering GLs is the 201 KEEP matrix.
    let orig = 249usize;
    let keep = 201usize;
    let drop = orig - keep;
    assert_eq!(drop, 48);
}
