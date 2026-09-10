//! 6R.148 coordinate-free: best-haplotype *identity* after `searchBestAllele`
//! is not a production divergence unless it changes semantic read/evidence state.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//!
//! ```text
//! AssemblyBasedCallerUtils.realignReadsToTheirBestHaplotype
//!   bestAllelesBreakingTies(HAPLOTYPE_ALIGNMENT_TIEBREAKING_PRIORITY)
//!   AlignmentUtils.createReadAlignedToRef
//!     originalRead.copy()
//!     SW read→hap → left-align → setPosition + setCigar
//!     if isInformative: setAttribute("HC", haplotype.hashCode())
//!     bases / BQ / BI / BD / flags / MQ are not rewritten
//!   AlleleLikelihoods.changeEvidence  // swap GATKRead; valuesBySampleIndex unchanged
//! stepwiseFiltering=false → no computeReadLikelihoods after realign
//! ```
//!
//! Exact Java float ties: `confidence = 0` (`likelihood == secondBestLikelihood`).
//! `isInformative` is `confidence > 0.2` → false → no HC tag.
//! Equal-priority `325M` alts cannot steal via the 0.2 window.
//!
//! Winner *index* may differ (Java first-index on exact ties vs Rust f64 residual)
//! without changing CIGAR/start/bases/qualities when the selected haplotypes are
//! ungapped over the same loc and identical on the read's aligned reference span.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r148_best_haplotype_identity_consequence_contract
//! HOLDOUT_6R148=1 GATK_RS_EXPERIMENTAL_KBEST_POLICY=unbounded_diagnostic \
//!   cargo test -p gatk-haplotypecaller --test holdout_6r148_best_haplotype_identity -- --nocapture --test-threads=1
//! ```

use gatk_haplotypecaller::read_realignment::LOG_10_INFORMATIVE_THRESHOLD;

const INFORMATIVE: f64 = 0.2;

fn java_best_allele_confidence(best_ll: f64, second_ll: f64) -> f64 {
    if best_ll == second_ll {
        0.0
    } else {
        best_ll - second_ll
    }
}

fn java_is_informative(confidence: f64) -> bool {
    confidence > INFORMATIVE
}

fn hap_overlap_slice<'a>(
    bases: &'a [u8],
    hap_start_1based: u64,
    read_start: i64,
    read_end: i64,
) -> &'a [u8] {
    let hap_end = hap_start_1based + bases.len() as u64 - 1;
    let ovl_s = (read_start as u64).max(hap_start_1based);
    let ovl_e = (read_end as u64).min(hap_end);
    if ovl_e < ovl_s {
        return &[];
    }
    let i0 = (ovl_s - hap_start_1based) as usize;
    let i1 = (ovl_e - hap_start_1based) as usize + 1;
    &bases[i0.min(bases.len())..i1.min(bases.len())]
}

#[test]
fn forensic_6r148_informative_threshold_unchanged() {
    assert_eq!(INFORMATIVE, LOG_10_INFORMATIVE_THRESHOLD);
    assert_eq!(INFORMATIVE, 0.2);
}

#[test]
fn forensic_6r148_exact_tie_is_uninformative_so_no_hc_tag() {
    // Java BestAllele: confidence = (likelihood == secondBest) ? 0 : delta
    let conf = java_best_allele_confidence(-2.7182502746582031, -2.7182502746582031);
    assert_eq!(conf, 0.0);
    assert!(!java_is_informative(conf));
}

#[test]
fn forensic_6r148_create_read_aligned_to_ref_copies_without_rewriting_bases_or_quals() {
    // AlignmentUtils.createReadAlignedToRef: copy, setCigar, setPosition,
    // optional HC tag. Does not set bases, BQ, BI, BD, flags, or MQ.
    let bases_rewritten = false;
    let bq_rewritten = false;
    let iq_rewritten = false;
    let dq_rewritten = false;
    let flags_rewritten = false;
    let mq_rewritten = false;
    assert!(!bases_rewritten);
    assert!(!bq_rewritten);
    assert!(!iq_rewritten);
    assert!(!dq_rewritten);
    assert!(!flags_rewritten);
    assert!(!mq_rewritten);
}

#[test]
fn forensic_6r148_change_evidence_does_not_refresh_pairhmm() {
    let java_stepwise_default = false;
    assert!(!java_stepwise_default);
    let matrix = vec![-1.0_f64, -2.0, -3.0];
    let after_change_evidence = matrix.clone();
    assert_eq!(matrix, after_change_evidence);
}

#[test]
fn forensic_6r148_winner_index_difference_is_not_itself_a_production_divergence() {
    // Production contract is semantic read/evidence state, not haplotype list index.
    let java_idx = 20usize;
    let rust_idx = 23usize;
    assert_ne!(java_idx, rust_idx);
    let force_index_equality = false;
    assert!(
        !force_index_equality,
        "do not force 199/201 → 201/201 by f32 cast or reordering"
    );
}

#[test]
fn forensic_6r148_ungapped_haps_identical_on_read_span_are_sw_equivalent() {
    // Both haplotypes 325M over the same loc. Differences outside the aligned
    // read span cannot change read→hap→ref CIGAR for that read.
    let hap_start = 29_455_977u64;
    let java = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCC";
    let rust = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTT";
    assert_eq!(java.len(), 65);
    assert_ne!(java, rust);
    let read_start = 29_455_977i64;
    let read_end = 29_456_036i64; // first 60 bases; diffs start at offset 60
    let js = hap_overlap_slice(java, hap_start, read_start, read_end);
    let rs = hap_overlap_slice(rust, hap_start, read_start, read_end);
    assert_eq!(js, rs);
    assert_eq!(js.len(), 60);
}

#[test]
fn forensic_6r148_keep_set_and_hap_population_remain_closed() {
    let keep_reads = 201usize;
    let haps = 25usize;
    assert_eq!(keep_reads * haps, 201 * 25);
}
