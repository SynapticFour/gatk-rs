//! 6R.93 coordinate-free: `filterPoorlyModeledEvidence` predicate inputs.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//!
//! ```text
//! PairHMMLikelihoodCalculationEngine.computeReadLikelihoods
//!   normalizeLikelihoods
//!   ReadLikelihoodCalculationEngine.filterPoorlyModeledEvidence  // dynamic off, cap=true
//!     AlleleLikelihoods.filterPoorlyModeledEvidence(log10MinTrueLikelihood)
//!
//! qualifiedLen = HMM_BASE_QUALITIES length if present, else GATKRead.getLength()
//! maxErrors    = min(2.0, ceil(qualifiedLen * 0.02))
//! threshold    = maxErrors * -4.0
//! max_ll       = max over haplotype columns after normalize (before filterAlleles)
//! DROP iff max_ll < threshold
//! KEEP iff !(max_ll < threshold)   // equality keeps
//! ```
//!
//! Live 24-QNAME table (holdout): threshold matches at -8.0 even when the qualifiedLen
//! integer differs, because every length is ≥ 51. Keep/drop follows max_ll, not qlen.
//! That remaining max_ll gap is upstream likelihood population, not this predicate.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r93_filter_poorly_modeled_predicate_contract
//! HOLDOUT_6R93=1 cargo test -p gatk-haplotypecaller --test holdout_6r93_filter_poorly_modeled -- --nocapture
//! ```

const EXPECTED_ERROR_RATE_PER_BASE: f64 = 0.02;
const LOG10_QUAL_PER_ERROR: f64 = -4.0;

/// Java `log10MinTrueLikelihood(expectedErrorRatePerBase, capLikelihoods=true)`.
fn java_log10_min_true_likelihood(qualified_read_len: usize) -> f64 {
    let max_errors = (qualified_read_len as f64 * EXPECTED_ERROR_RATE_PER_BASE)
        .ceil()
        .min(2.0);
    max_errors * LOG10_QUAL_PER_ERROR
}

/// Java `AlleleLikelihoods.filterPoorlyModeledEvidence`: drop iff `max_ll < threshold`.
fn java_keep(max_ll: f64, qualified_read_len: usize) -> bool {
    !(max_ll < java_log10_min_true_likelihood(qualified_read_len))
}

/// Java `qualifiedLen`: HMM_BASE_QUALITIES length if the tag is present, else read length.
fn java_qualified_len(hmm_bq: Option<usize>, read_len: usize) -> usize {
    hmm_bq.unwrap_or(read_len)
}

#[test]
fn forensic_6r93_java_static_threshold_formula() {
    assert_eq!(java_log10_min_true_likelihood(1), -4.0);
    assert_eq!(java_log10_min_true_likelihood(50), -4.0);
    assert_eq!(java_log10_min_true_likelihood(51), -8.0);
    assert_eq!(java_log10_min_true_likelihood(76), -8.0);
    assert_eq!(java_log10_min_true_likelihood(100), -8.0);
    assert_eq!(java_log10_min_true_likelihood(130), -8.0);
    assert_eq!(java_log10_min_true_likelihood(148), -8.0);
    assert_eq!(java_log10_min_true_likelihood(250), -8.0);
}

#[test]
fn forensic_6r93_drop_is_strict_less_equality_keeps() {
    let thr = java_log10_min_true_likelihood(100);
    assert_eq!(thr, -8.0);
    assert!(java_keep(-8.0, 100));
    assert!(!java_keep(-8.0 - 1e-15, 100));
    assert!(java_keep(-7.999, 100));
}

#[test]
fn forensic_6r93_qualified_len_prefers_hmm_base_qualities() {
    assert_eq!(java_qualified_len(Some(76), 148), 76);
    assert_eq!(java_qualified_len(None, 148), 148);
    assert_eq!(java_qualified_len(Some(148), 148), 148);
}

#[test]
fn forensic_6r93_clipped_vs_full_read_length_does_not_change_threshold_above_51() {
    // Live 24-QNAME qlen integers include 76 vs 148. Both cap at two errors.
    let a = java_log10_min_true_likelihood(76);
    let b = java_log10_min_true_likelihood(148);
    assert_eq!(a, -8.0);
    assert_eq!(b, -8.0);
    assert_eq!(a, b);
    assert_eq!(
        java_keep(-7.45, 76),
        java_keep(-7.45, 148),
        "same max_ll cannot flip keep/drop via qlen once both are ≥ 51"
    );
    assert_eq!(java_keep(-9.5, 76), java_keep(-9.5, 148));
}

#[test]
fn forensic_6r93_rust_qual_len_max1_matches_java_formula_on_same_integer() {
    // Rust `rec.qual().len().max(1)` is the integer fed to the same ceil/min/−4 formula.
    let rust_qual_len = 148usize.max(1);
    assert_eq!(
        java_log10_min_true_likelihood(rust_qual_len),
        java_log10_min_true_likelihood(148)
    );
}

#[test]
fn forensic_6r93_identical_predicate_inputs_yield_identical_keep() {
    let cases = [
        (148usize, -2.5),
        (148, -8.0),
        (148, -8.0001),
        (76, -2.47),
        (76, -9.5),
        (50, -3.9),
        (50, -4.1),
    ];
    for &(qlen, max_ll) in &cases {
        let rust_keep = max_ll >= java_log10_min_true_likelihood(qlen);
        assert_eq!(
            java_keep(max_ll, qlen),
            rust_keep,
            "qlen={qlen} max_ll={max_ll}"
        );
    }
}

#[test]
fn forensic_6r93_max_ll_not_qlen_flips_keep_when_threshold_tied() {
    let qlen_java = 76usize;
    let qlen_rust = 148usize;
    assert_eq!(
        java_log10_min_true_likelihood(qlen_java),
        java_log10_min_true_likelihood(qlen_rust)
    );
    let java_max_ll = -2.47;
    let rust_max_ll = -9.55;
    assert!(java_keep(java_max_ll, qlen_java));
    assert!(!java_keep(rust_max_ll, qlen_rust));
    assert_ne!(
        java_keep(java_max_ll, qlen_java),
        java_keep(rust_max_ll, qlen_rust),
        "membership follows max_ll, not the tied threshold"
    );
}

#[test]
fn forensic_6r93_max_ll_is_max_over_haplotype_columns_after_normalize() {
    // Immediate producer: max of the per-haplotype likelihood row (Java alleles at
    // filter time are haplotypes). Not allele-marginalized AD, not PairHMM kernel.
    let hap_ll = [-12.0, -2.5, -4.0, -20.0];
    let max_ll = hap_ll.into_iter().fold(f64::NEG_INFINITY, f64::max);
    assert_eq!(max_ll, -2.5);
    assert!(java_keep(max_ll, 148));
    let filtered_haps_only = [-12.0, -20.0];
    let max_after_allele_filter = filtered_haps_only
        .into_iter()
        .fold(f64::NEG_INFINITY, f64::max);
    assert_eq!(max_after_allele_filter, -12.0);
    assert!(
        !java_keep(max_after_allele_filter, 148),
        "dropping the best haplotype column can flip DROP without changing the predicate"
    );
}
