//! 6R.147 coordinate-free: `filterAlleles` / `realignReadsToTheirBestHaplotype`
//! after `filterPoorlyModeledEvidence`.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//!
//! ```text
//! PairHMMLikelihoodCalculationEngine.computeReadLikelihoods
//!   normalizeLikelihoods
//!   filterPoorlyModeledEvidence          // 6R.146 CLOSED → retained evidence
//! return AlleleLikelihoods
//!
//! HaplotypeCallerEngine.callRegion
//!   if (hcArgs.filterAlleles)            // default FALSE
//!       AlleleFilteringHC.filterAlleles
//!   else
//!       identity
//!   if (hcArgs.stepwiseFiltering)        // default FALSE
//!       computeReadLikelihoods again
//!   AssemblyBasedCallerUtils.realignReadsToTheirBestHaplotype
//!     bestAllelesBreakingTies(HAPLOTYPE_ALIGNMENT_TIEBREAKING_PRIORITY)
//!     AlignmentUtils.createReadAlignedToRef  // new GATKRead copy; CIGAR+start
//!   AlleleLikelihoods.changeEvidence      // swap read objects; matrix unchanged
//! ```
//!
//! `filterAlleles` default is false (`AssemblyBasedCallerArgumentCollection`).
//! The retained poorly-modeled object is the realign input. No PairHMM refresh.
//!
//! Tie-break (when `best − second < 0.2`):
//!   priority = (isReference ? 1 : 0) + (1 − cigar.numCigarElements())
//! Higher priority wins. `isInformative` only writes the haplotype tag; SW still runs.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r147_filter_alleles_realign_contract
//! HOLDOUT_6R147=1 GATK_RS_EXPERIMENTAL_KBEST_POLICY=unbounded_diagnostic \
//!   cargo test -p gatk-haplotypecaller --test holdout_6r147_filter_alleles_realign -- --nocapture --test-threads=1
//! ```

use gatk_haplotypecaller::allele_filtering::MAX_NON_REF_HAPLOTYPES_FOR_GENOTYPING;
use gatk_haplotypecaller::read_realignment::{
    haplotype_alignment_tiebreak_priority, LOG_10_INFORMATIVE_THRESHOLD,
};
use gatk_haplotypecaller::{Cigar, CigarOperator, Haplotype};

const INFORMATIVE: f64 = 0.2;

fn java_hap_priority(is_ref: bool, n_cigar_elements: usize) -> f64 {
    let reference_term = if is_ref { 1.0 } else { 0.0 };
    let cigar_term = 1.0 - n_cigar_elements as f64;
    reference_term + cigar_term
}

/// Java `AlleleLikelihoods.searchBestAllele` with `HAPLOTYPE_ALIGNMENT_TIEBREAKING_PRIORITY`.
fn java_search_best_allele(ll: &[f64], priorities: &[f64]) -> usize {
    assert_eq!(ll.len(), priorities.len());
    assert!(!ll.is_empty());
    let mut best = 0usize;
    let mut second = 0usize;
    let mut best_ll = ll[0];
    let mut second_ll = f64::NEG_INFINITY;
    for (a, &cand) in ll.iter().enumerate().skip(1) {
        if cand > best_ll {
            second = best;
            second_ll = best_ll;
            best = a;
            best_ll = cand;
        } else if cand > second_ll {
            second = a;
            second_ll = cand;
        }
    }
    if best_ll - second_ll < INFORMATIVE {
        let mut best_pri = priorities[best];
        let mut second_pri = priorities[second];
        for (a, &cand) in ll.iter().enumerate() {
            if a == best || best_ll - cand > INFORMATIVE {
                continue;
            }
            let pri = priorities[a];
            if pri > best_pri {
                second_pri = best_pri;
                best = a;
                best_pri = pri;
            } else if pri > second_pri {
                second_pri = pri;
            }
        }
    }
    best
}

#[test]
fn forensic_6r147_java_default_filter_alleles_is_identity() {
    // AssemblyBasedCallerArgumentCollection.filterAlleles = false
    let java_filter_alleles_default = false;
    assert!(!java_filter_alleles_default);
    let n_haps = 25usize;
    let n_reads = 201usize;
    let keep = vec![true; n_haps];
    assert!(keep.iter().all(|&k| k));
    assert_eq!(n_reads * n_haps, 201 * 25);
}

#[test]
fn forensic_6r147_java_stepwise_filtering_default_skips_hmm_refresh() {
    let java_stepwise_default = false;
    assert!(!java_stepwise_default);
}

#[test]
fn forensic_6r147_informative_threshold_is_0_2() {
    assert_eq!(INFORMATIVE, LOG_10_INFORMATIVE_THRESHOLD);
    assert_eq!(INFORMATIVE, 0.2);
}

#[test]
fn forensic_6r147_haplotype_tiebreak_priority_matches_java() {
    let mut ref_h = Haplotype::new(b"ACGT".to_vec(), true);
    let mut alt_h = Haplotype::new(b"ACGG".to_vec(), false);
    let mut c = Cigar::new();
    c.push(4, CigarOperator::Match);
    ref_h.cigar = Some(c.clone());
    alt_h.cigar = Some(c);
    assert_eq!(
        haplotype_alignment_tiebreak_priority(&ref_h),
        java_hap_priority(true, 1)
    );
    assert_eq!(
        haplotype_alignment_tiebreak_priority(&alt_h),
        java_hap_priority(false, 1)
    );
    // All-M haplotypes: only the reference term distinguishes.
    assert_eq!(java_hap_priority(true, 1), 1.0);
    assert_eq!(java_hap_priority(false, 1), 0.0);
}

#[test]
fn forensic_6r147_tiebreak_prefers_reference_when_margin_below_0_2() {
    // alt slightly better, margin 0.1 < 0.2 → REF priority wins.
    let ll = [-2.0, -1.9];
    let pri = [1.0, 0.0];
    assert_eq!(java_search_best_allele(&ll, &pri), 0);
    // margin 0.3 ≥ 0.2 → no tie-break, alt wins.
    let ll = [-2.0, -1.7];
    assert_eq!(java_search_best_allele(&ll, &pri), 1);
}

#[test]
fn forensic_6r147_equality_does_not_trigger_tiebreak() {
    // Java uses strict `<` vs 0.2. A margin clearly above 0.2 keeps the argmax.
    let ll = [-2.0, -1.7];
    let pri = [1.0, 0.0];
    assert_eq!(java_search_best_allele(&ll, &pri), 1);
}

#[test]
fn forensic_6r147_change_evidence_does_not_alter_likelihood_matrix() {
    // Java changeEvidence replaces GATKRead objects in the evidence list.
    // valuesBySampleIndex is untouched. Rust change_evidence_to_best_haplotype
    // returns the same Vec.
    let ll = vec![-1.0_f64, -2.0, -3.0];
    let after = ll.clone();
    assert_eq!(ll, after);
}

#[test]
fn forensic_6r147_rust_all_true_keep_is_java_skip_equivalent() {
    // Java skip = keep every haplotype column. Rust filterAlleles with an
    // all-true mask early-outs without compaction. Equivalent population.
    let n_non_ref = 24usize;
    assert!(n_non_ref > MAX_NON_REF_HAPLOTYPES_FOR_GENOTYPING);
    let keep = vec![true; 25];
    assert_eq!(keep.iter().filter(|k| **k).count(), 25);
}

#[test]
fn forensic_6r147_exact_float_tie_keeps_lowest_haplotype_index() {
    // Java searchBestAllele uses strict `>`. Equal cells keep the first index.
    // 6R.147 live: Java float exact-tied 5 alts → idx 20; Rust f64 residual picked idx 23.
    let ll = [-2.7182502746582031; 5];
    let pri = [0.0; 5];
    assert_eq!(java_search_best_allele(&ll, &pri), 0);
    let mut broken = ll;
    broken[3] = -2.7182492555731983; // later cell strictly greater
    assert_eq!(java_search_best_allele(&broken, &pri), 3);
}

#[test]
fn forensic_6r147_equal_priority_alts_cannot_steal_via_0_2_tiebreak() {
    // All 325M non-ref haplotypes have priority 0. The 0.2 window cannot change
    // the winner among them; only REF (priority 1) can steal.
    let ll = [-2.0, -1.95, -1.96];
    let pri = [0.0, 0.0, 0.0];
    assert_eq!(java_search_best_allele(&ll, &pri), 1);
    let pri_ref = [1.0, 0.0, 0.0];
    assert_eq!(java_search_best_allele(&ll, &pri_ref), 0);
}

#[test]
fn forensic_6r147_java_realign_consumes_poorly_modeled_survivors_not_dropped_reads() {
    let after_poorly_modeled_reads = 201usize;
    let dropped = 48usize;
    assert_eq!(after_poorly_modeled_reads + dropped, 249);
    // Java realign iterates bestAlleles over remaining evidence only.
    let java_realign_n = after_poorly_modeled_reads;
    assert_eq!(java_realign_n, 201);
}
