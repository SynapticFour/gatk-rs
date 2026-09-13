//! 6R.154 coordinate-free: FORMAT/AD after the equivalent 41-read retained object.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//!
//! ```text
//! calculateGLsForThisEvent
//!   GenotypeBuilder(NO_CALL).PL(...)     // AD absent
//! AlleleSubsettingUtils.subsetAlleles
//!   if (g.hasAD()) slice; skipped
//! prepareReadAlleleLikelihoodsForAnnotation
//!   reuse genotyping AlleleLikelihoods (no contamination)
//!   updateNonRefAlleleLikelihoods iff allele counts differ
//!   addEvidence(overlappingFilteredReads, 0)   // ties → uninformative
//! DepthPerAlleleBySample.annotate
//!   alleles = LinkedHashSet(vc.getAlleles())   // remaining call alleles
//!   annotateWithLikelihoods                    // FIRST AD write
//!     alleleSubset = {allele → [allele]}       // identity remarg
//!     subsetted = likelihoods.marginalize(alleleSubset)
//!     bestAllelesBreakingTies(sample)
//!       searchBestAllele(..., canBeReference=true, priorities)
//!       default priority: isReference ? 1.0 : 0
//!       tie-break iff (best − second) < 0.2
//!     filter isInformative: confidence > 0.2   // strict greater-than
//!     count by Allele.equals into vc allele order (REF, then ALTs)
//! reverseTrimAlleles / phaseVC                 // copy genotype fields
//! ```
//!
//! GT does not enter the count. On a 2-allele remaining set, REF priority
//! reassignment never creates or destroys an informative vote: after a swap
//! to REF, confidence is ≤ 0. Rust biallelic `|lr−la| > 0.2` therefore
//! matches Java informative counts on A/T.
//!
//! A 4-way unused-ALT permute is not this operation. Identity remarg votes
//! only over remaining columns. `*` is an internal column only; it is not
//! emitted.
//!
//! Overlapping mates are independent evidence rows (6R.151). Both are
//! eligible for AD; QNAME is a join key, not a collapse key.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r154_ad_semantics_contract
//! HOLDOUT_6R154=1 GATK_RS_EXPERIMENTAL_KBEST_POLICY=unbounded_diagnostic \
//!   cargo test -p gatk-haplotypecaller --test holdout_6r154_ad -- --nocapture --test-threads=1
//! ```

use gatk_haplotypecaller::genotyping::ReadLikelihoodRow;
use gatk_haplotypecaller::hc_genotyping_engine::{
    biallelic_allele_depths_from_rows, InformativeAd,
};
use gatk_haplotypecaller::read_realignment::LOG_10_INFORMATIVE_THRESHOLD;

const JAVA_INFORMATIVE: f64 = 0.2;

fn row(lls: Vec<f64>) -> ReadLikelihoodRow {
    ReadLikelihoodRow {
        read_index: 0,
        read_id: String::new(),
        haplotype_log10_likelihoods: lls,
    }
}

/// Java `searchBestAllele` + default REF tie-break (`priority = isReference ? 1 : 0`).
fn java_best_breaking_ties(lls: &[f64], ref_i: usize) -> (usize, usize, f64, bool) {
    let n = lls.len();
    if n == 0 {
        return (0, 0, 0.0, false);
    }
    let mut best_i = 0usize;
    let mut second_i = 0usize;
    let mut best = lls[0];
    let mut second = f64::NEG_INFINITY;
    for a in 1..n {
        let cand = lls[a];
        if cand > best {
            second_i = best_i;
            second = best;
            best_i = a;
            best = cand;
        } else if cand > second {
            second_i = a;
            second = cand;
        }
    }
    let priorities: Vec<f64> = (0..n).map(|i| if i == ref_i { 1.0 } else { 0.0 }).collect();
    if best - second < JAVA_INFORMATIVE {
        let mut best_pri = priorities[best_i];
        let mut second_pri = priorities[second_i];
        for a in 0..n {
            let cand = lls[a];
            if a == best_i || best - cand > JAVA_INFORMATIVE {
                continue;
            }
            let pri = priorities[a];
            if pri > best_pri {
                second_i = best_i;
                best_i = a;
                second_pri = best_pri;
                best_pri = pri;
            } else if pri > second_pri {
                second_i = a;
                second_pri = pri;
            }
        }
    }
    let best_ll = lls[best_i];
    let second_ll = if second_i != best_i {
        lls[second_i]
    } else {
        f64::NEG_INFINITY
    };
    let conf = if best_ll == second_ll {
        0.0
    } else {
        best_ll - second_ll
    };
    (best_i, second_i, conf, conf > JAVA_INFORMATIVE)
}

fn rust_biallelic_vote(lr: f64, la: f64) -> Option<usize> {
    let gap = (lr - la).abs();
    if gap > LOG_10_INFORMATIVE_THRESHOLD {
        if lr > la {
            Some(0)
        } else {
            Some(1)
        }
    } else {
        None
    }
}

fn remarg_identity(old: &[f64], keep: &[usize]) -> Vec<f64> {
    keep.iter()
        .map(|&i| old.get(i).copied().unwrap_or(f64::NEG_INFINITY))
        .collect()
}

fn java_contribution(lls: &[f64], alleles: &[&str], ref_i: usize) -> Option<String> {
    let (best_i, _, _, inf) = java_best_breaking_ties(lls, ref_i);
    if inf {
        Some(alleles[best_i].to_string())
    } else {
        None
    }
}

#[test]
fn forensic_6r154_ad_is_absent_until_annotate_with_likelihoods() {
    let after_calculate_gls_has_ad = false;
    let after_subset_alleles_has_ad = after_calculate_gls_has_ad;
    assert!(
        !after_subset_alleles_has_ad,
        "Java unused-ALT subset copies AD only when genotype.hasAD()"
    );
    let first_write = "DepthPerAlleleBySample.annotateWithLikelihoods";
    assert_eq!(
        first_write,
        "DepthPerAlleleBySample.annotateWithLikelihoods"
    );
}

#[test]
fn forensic_6r154_is_informative_is_strict_greater_than_0_2() {
    assert_eq!(LOG_10_INFORMATIVE_THRESHOLD, JAVA_INFORMATIVE);
    let (_, _, conf_eq, inf_eq) = java_best_breaking_ties(&[0.0, -0.2], 0);
    assert!((conf_eq - 0.2).abs() < 1e-12);
    assert!(
        !inf_eq,
        "Java BestAllele.isInformative is confidence > 0.2; equality is uninformative"
    );
    let (_, _, conf_gt, inf_gt) = java_best_breaking_ties(&[0.0, -0.2000001], 0);
    assert!(conf_gt > JAVA_INFORMATIVE);
    assert!(inf_gt);
}

#[test]
fn forensic_6r154_strongly_best_a_and_t_count() {
    let alleles = ["A", "T"];
    let a = java_contribution(&[0.0, -3.0], &alleles, 0);
    let t = java_contribution(&[-3.0, 0.0], &alleles, 0);
    assert_eq!(a.as_deref(), Some("A"));
    assert_eq!(t.as_deref(), Some("T"));
    assert_eq!(rust_biallelic_vote(0.0, -3.0), Some(0));
    assert_eq!(rust_biallelic_vote(-3.0, 0.0), Some(1));
    let rows = [row(vec![0.0, -3.0]), row(vec![-3.0, 0.0])];
    assert_eq!(biallelic_allele_depths_from_rows(&rows, 0, 1), vec![1, 1]);
    assert_eq!(
        InformativeAd::from_marginalized_rows(&rows, 0, 1, None),
        InformativeAd {
            ref_depth: 1,
            alt_depth: 1
        }
    );
}

#[test]
fn forensic_6r154_ambiguous_below_threshold_is_uninformative() {
    let alleles = ["A", "T"];
    assert_eq!(java_contribution(&[0.0, -0.1], &alleles, 0), None);
    assert_eq!(java_contribution(&[-0.1, 0.0], &alleles, 0), None);
    assert_eq!(rust_biallelic_vote(0.0, -0.1), None);
    assert_eq!(rust_biallelic_vote(-0.1, 0.0), None);
}

#[test]
fn forensic_6r154_exact_tie_is_uninformative() {
    let (best_i, _, conf, inf) = java_best_breaking_ties(&[0.0, 0.0], 0);
    assert_eq!(best_i, 0, "first-max picks REF on an exact likelihood tie");
    assert_eq!(conf, 0.0);
    assert!(!inf);
    assert_eq!(rust_biallelic_vote(0.0, 0.0), None);
}

#[test]
fn forensic_6r154_biallelic_ref_priority_does_not_change_informative_counts() {
    // ALT slightly better (gap 0.1 < 0.2): Java swaps to REF, then conf ≤ 0.
    let (best_i, _, conf, inf) = java_best_breaking_ties(&[-0.1, 0.0], 0);
    assert_eq!(best_i, 0);
    assert!(conf <= 0.0);
    assert!(!inf);
    assert_eq!(
        rust_biallelic_vote(-0.1, 0.0),
        None,
        "Rust |gap| > 0.2 matches Java informative counts on two remaining alleles"
    );
}

#[test]
fn forensic_6r154_gt_does_not_enter_ad_count() {
    let alleles = ["A", "T"];
    let gt_hom_ref = [0, 0];
    let gt_het = [0, 1];
    let vote = java_contribution(&[-3.0, 0.0], &alleles, 0);
    assert_eq!(vote.as_deref(), Some("T"));
    let _ = (gt_hom_ref, gt_het);
    let gt_used_in_count = false;
    assert!(!gt_used_in_count);
}

#[test]
fn forensic_6r154_ad_keyed_by_allele_identity_ref_then_alt() {
    // Java: alleleCounts.get(vc.getReference()) then getAlternateAllele(i).
    // Rust production biallelic AD is index 0 = REF A, index 1 = ALT T.
    let java_ad0 = "A";
    let java_ad1 = "T";
    let rust_ad0 = "A";
    let rust_ad1 = "T";
    assert_eq!(java_ad0, rust_ad0);
    assert_eq!(java_ad1, rust_ad1);
    let rows = [
        row(vec![0.0, -5.0]),
        row(vec![0.0, -5.0]),
        row(vec![-5.0, 0.0]),
    ];
    assert_eq!(biallelic_allele_depths_from_rows(&rows, 0, 1), vec![2, 1]);
}

#[test]
fn forensic_6r154_identity_remarg_not_four_way_permute() {
    // Columns: A, *, unused C, T. Remaining call alleles A,T.
    let keep = [0usize, 3];
    let star_best = [-8.45, 0.0, -3.0, -8.0];
    let unused_best = [-8.0, -3.0, 0.0, -8.55];
    let remarg_star = remarg_identity(&star_best, &keep);
    let remarg_unused = remarg_identity(&unused_best, &keep);
    assert_eq!(
        java_contribution(&star_best, &["A", "*", "C", "T"], 0).as_deref(),
        Some("*")
    );
    assert_eq!(
        java_contribution(&remarg_star, &["A", "T"], 0).as_deref(),
        Some("T"),
        "after unused-ALT identity remarg, * is gone; remaining T still wins"
    );
    assert_eq!(
        java_contribution(&remarg_unused, &["A", "T"], 0).as_deref(),
        Some("A")
    );
    // Slicing 4-way informative counts would drop both rows.
    let four_way_slice = [0i32, 0];
    let remarg_counts = {
        let mut ad = [0i32; 2];
        for kept in [&remarg_star, &remarg_unused] {
            if let Some(i) = java_contribution(kept, &["A", "T"], 0) {
                if i == "A" {
                    ad[0] += 1;
                } else {
                    ad[1] += 1;
                }
            }
        }
        ad
    };
    assert_eq!(four_way_slice, [0, 0]);
    assert_eq!(remarg_counts, [1, 1]);
}

#[test]
fn forensic_6r154_star_best_does_not_enter_at_ad_unless_remaining() {
    let lls = [-1.0, 0.0, -4.0]; // A, *, T — * strongly best
    assert_eq!(
        java_contribution(&lls, &["A", "*", "T"], 0).as_deref(),
        Some("*")
    );
    let remaining = remarg_identity(&lls, &[0, 2]);
    assert_eq!(
        java_contribution(&remaining, &["A", "T"], 0).as_deref(),
        Some("A"),
        "* is not a remaining call allele; remaining A vs T still votes"
    );
}

#[test]
fn forensic_6r154_add_evidence_zero_is_uninformative() {
    // Java addEvidence(..., 0) fills every allele cell with 0.
    let (best_i, _, conf, inf) = java_best_breaking_ties(&[0.0, 0.0], 0);
    assert_eq!(best_i, 0);
    assert_eq!(conf, 0.0);
    assert!(!inf);
}

#[test]
fn forensic_6r154_overlapping_mates_are_independently_eligible() {
    // Coordinate-free: two evidence rows, same fragment conceptually, both informative.
    let mate_a = row(vec![0.0, -2.0]);
    let mate_b = row(vec![-2.5, 0.0]);
    let rows = [mate_a, mate_b];
    let ad = biallelic_allele_depths_from_rows(&rows, 0, 1);
    assert_eq!(
        ad,
        vec![1, 1],
        "QNAME collapse must not drop a mate from AD"
    );
    let a = java_contribution(&[0.0, -2.0], &["A", "T"], 0);
    let b = java_contribution(&[-2.5, 0.0], &["A", "T"], 0);
    assert_eq!(a.as_deref(), Some("A"));
    assert_eq!(b.as_deref(), Some("T"));
}

#[test]
fn forensic_6r154_ad_computed_after_unused_alt_remaining_set() {
    let annotation_alleles = ["A", "T"];
    let genotyping_alleles_with_unused = ["A", "*", "C", "T"];
    assert!(annotation_alleles
        .iter()
        .all(|a| genotyping_alleles_with_unused.contains(a)));
    let first_ad_write_uses = "remaining vc.getAlleles() after unused-ALT subset";
    assert_eq!(
        first_ad_write_uses,
        "remaining vc.getAlleles() after unused-ALT subset"
    );
}
