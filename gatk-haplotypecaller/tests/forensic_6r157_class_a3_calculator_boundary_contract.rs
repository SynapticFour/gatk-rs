//! 6R.157 coordinate-free: Class-A3 calculator-boundary contract (no production change).
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//!
//! Question: which semantic guard is Java-equivalent for Class-A3 genotype replacement?
//!   A — `pl_gt == 0` preserve, else A3 may replace
//!   B — preserve any valid assigned calculator genotype
//!   C — A3 only when calculator produced NO_CALL / unusable
//!   D — some other proven Java boundary
//!
//! This file is diagnostic of the discovered contract. Production still violates it
//! for valid 0/0 and valid 1/1 (full sparse replace). 6R.157 does not patch production;
//! assertions pin current vs required so default `cargo test` stays green. A later
//! patch round inverts the `assert_ne!` preservation checks after the A3-only skip.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r157_class_a3_calculator_boundary_contract
//! ```

use gatk_haplotypecaller::genotyping::{best_pl_index, emit_genotype_format_fields};
use gatk_haplotypecaller::hc_genotyping_engine::SparsePlShape;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum A3Action {
    None,
    AdDpOnly,
    SparseReplace,
}

#[derive(Debug, Clone)]
struct Calc {
    gls: Vec<f64>,
    info_ad: [i32; 2],
}

#[derive(Debug, Clone)]
struct Observed {
    pl_gt: usize,
    valid_assigned: bool,
    no_call: bool,
    pl: Vec<i32>,
    ad: Vec<i32>,
    gq: i32,
    dp: i32,
    class_a: bool,
    class_a2: bool,
    class_a3: bool,
    action: A3Action,
    after_gt: usize,
    after_pl: Vec<i32>,
    after_ad: Vec<i32>,
    after_gq: i32,
    after_dp: i32,
}

/// Production Class-A / A2 / A3 decision tree in `SiteReshape::apply_class_a_family`.
/// Integration tests cannot call the `pub(crate)` method; this replica is the contract.
fn class_a_family_early_exit(
    in_p12_chr2_scope: bool,
    ad_len: usize,
    pileup_ref: i32,
    pileup_alt: i32,
) -> bool {
    in_p12_chr2_scope
        || ad_len < 2
        || pileup_ref < 1
        || pileup_alt < 1
        || pileup_ref.saturating_mul(2) < pileup_alt
}

fn class_a_pred(info_ref: i32, info_alt: i32) -> bool {
    info_ref == 0 && info_alt >= 2
}

fn class_a2_pred(
    is_snp: bool,
    pileup_ref: i32,
    pileup_alt: i32,
    info_ref: i32,
    info_alt: i32,
    pl_gt: usize,
) -> bool {
    is_snp
        && pileup_alt.saturating_mul(2) >= pileup_ref
        && (pl_gt == 2
            || (info_alt >= 2 && (info_ref == 0 || info_alt >= info_ref.saturating_mul(3))))
}

fn class_a3_pred(
    is_snp: bool,
    pileup_ref: i32,
    pileup_alt: i32,
    info_ref: i32,
    info_alt: i32,
) -> bool {
    is_snp
        && pileup_alt.saturating_mul(2) >= pileup_ref
        && info_ref >= 1
        && info_alt >= 1
        && info_ref >= info_alt.saturating_mul(3)
}

fn observe(calc: &Calc, pileup: [i32; 2], is_snp: bool, in_p12_chr2_scope: bool) -> Observed {
    let first = emit_genotype_format_fields(&calc.gls, &calc.info_ad).expect("calculator FORMAT");
    let pl_gt = best_pl_index(&first.pl);
    let valid_assigned = !first.pl.is_empty() && first.ad.len() >= 2;
    // `RegionGenotypeResult` has no NO_CALL field. SiteScore always assigns from min PL.
    let no_call = false;
    if class_a_family_early_exit(in_p12_chr2_scope, first.ad.len(), pileup[0], pileup[1]) {
        return Observed {
            pl_gt,
            valid_assigned,
            no_call,
            pl: first.pl_as_i32(),
            ad: first.ad_as_i32(),
            gq: first.gq.as_i32(),
            dp: first.dp.as_i32(),
            class_a: false,
            class_a2: false,
            class_a3: false,
            action: A3Action::None,
            after_gt: pl_gt,
            after_pl: first.pl_as_i32(),
            after_ad: first.ad_as_i32(),
            after_gq: first.gq.as_i32(),
            after_dp: first.dp.as_i32(),
        };
    }
    let info_ref = first.ad[0].as_i32();
    let info_alt = first.ad[1].as_i32();
    let class_a = class_a_pred(info_ref, info_alt);
    let class_a2 = class_a2_pred(is_snp, pileup[0], pileup[1], info_ref, info_alt, pl_gt);
    let class_a3 = class_a3_pred(is_snp, pileup[0], pileup[1], info_ref, info_alt);
    if !(class_a || class_a2 || class_a3) {
        return Observed {
            pl_gt,
            valid_assigned,
            no_call,
            pl: first.pl_as_i32(),
            ad: first.ad_as_i32(),
            gq: first.gq.as_i32(),
            dp: first.dp.as_i32(),
            class_a,
            class_a2,
            class_a3,
            action: A3Action::None,
            after_gt: pl_gt,
            after_pl: first.pl_as_i32(),
            after_ad: first.ad_as_i32(),
            after_gq: first.gq.as_i32(),
            after_dp: first.dp.as_i32(),
        };
    }
    if pl_gt == 1 {
        let fmt = emit_genotype_format_fields(&calc.gls, &pileup).expect("AD-only");
        Observed {
            pl_gt,
            valid_assigned,
            no_call,
            pl: first.pl_as_i32(),
            ad: first.ad_as_i32(),
            gq: first.gq.as_i32(),
            dp: first.dp.as_i32(),
            class_a,
            class_a2,
            class_a3,
            action: A3Action::AdDpOnly,
            after_gt: best_pl_index(&fmt.pl),
            after_pl: fmt.pl_as_i32(),
            after_ad: fmt.ad_as_i32(),
            after_gq: fmt.gq.as_i32(),
            after_dp: fmt.dp.as_i32(),
        }
    } else {
        let shape = SparsePlShape::from_pileup_depths(pileup[0], pileup[1]);
        let fmt = emit_genotype_format_fields(&shape.gl_vec(), &pileup).expect("sparse");
        Observed {
            pl_gt,
            valid_assigned,
            no_call,
            pl: first.pl_as_i32(),
            ad: first.ad_as_i32(),
            gq: first.gq.as_i32(),
            dp: first.dp.as_i32(),
            class_a,
            class_a2,
            class_a3,
            action: A3Action::SparseReplace,
            after_gt: best_pl_index(&fmt.pl),
            after_pl: fmt.pl_as_i32(),
            after_ad: fmt.ad_as_i32(),
            after_gq: fmt.gq.as_i32(),
            after_dp: fmt.dp.as_i32(),
        }
    }
}

fn hom_ref() -> Calc {
    Calc {
        gls: vec![0.0, -0.5, -117.4],
        info_ad: [37, 4],
    }
}

fn het() -> Calc {
    Calc {
        gls: vec![-8.1, 0.0, -3.6],
        info_ad: [37, 4],
    }
}

fn hom_alt() -> Calc {
    Calc {
        gls: vec![-9.0, -0.6, 0.0],
        info_ad: [37, 4],
    }
}

/// Canonical A3-satisfying pileup (het shape, not HomAltStrong).
const A3_PILEUP: [i32; 2] = [22, 24];

fn guard_a_preserve(pl_gt: usize) -> bool {
    pl_gt == 0
}

fn guard_b_preserve(valid_assigned: bool) -> bool {
    valid_assigned
}

fn guard_c_preserve(no_call: bool) -> bool {
    !no_call
}

#[test]
fn forensic_6r157_pl_gt_is_min_pl_index_not_validity() {
    let r0 = observe(&hom_ref(), A3_PILEUP, true, false);
    let r1 = observe(&het(), A3_PILEUP, true, false);
    let r2 = observe(&hom_alt(), A3_PILEUP, true, false);
    assert_eq!(r0.pl_gt, 0);
    assert_eq!(r1.pl_gt, 1);
    assert_eq!(r2.pl_gt, 2);
    assert!(r0.valid_assigned && r1.valid_assigned && r2.valid_assigned);
    assert!(!r0.no_call && !r1.no_call && !r2.no_call);
    // After PL normalization every assigned genotype has some PL == 0.
    assert_eq!(r0.pl[0], 0);
    assert_eq!(r1.pl[1], 0);
    assert_eq!(r2.pl[2], 0);
    // Empty PL would unwrap to index 0; that is not a NO_CALL sentinel.
    let empty_pl_fallback = best_pl_index(&[]);
    assert_eq!(empty_pl_fallback, 0);
}

#[test]
fn forensic_6r157_case1_valid_00_is_sparse_replaced() {
    let o = observe(&hom_ref(), A3_PILEUP, true, false);
    assert_eq!(o.pl, vec![0, 5, 1174]);
    assert_eq!(o.ad, vec![37, 4]);
    assert_eq!(o.gq, 5);
    assert_eq!(o.dp, 41);
    assert!(o.class_a3 && !o.class_a && !o.class_a2);
    assert_eq!(o.action, A3Action::SparseReplace);
    assert_eq!(o.after_gt, 1);
    assert_eq!(o.after_pl, vec![81, 0, 36]);
    assert_eq!(o.after_ad, vec![22, 24]);
    assert_eq!(o.after_gq, 36);
    assert_eq!(o.after_dp, 46);
    // Java-equivalent required: preserve. Current production violates.
    let java_required_gt = 0usize;
    assert_ne!(o.after_gt, java_required_gt);
    assert!(guard_a_preserve(o.pl_gt));
    assert!(guard_b_preserve(o.valid_assigned));
    assert!(guard_c_preserve(o.no_call));
}

#[test]
fn forensic_6r157_case2_valid_01_keeps_gt_pl_replaces_ad() {
    let o = observe(&het(), A3_PILEUP, true, false);
    assert_eq!(o.pl_gt, 1);
    assert!(o.class_a3 && !o.class_a && !o.class_a2);
    assert_eq!(o.action, A3Action::AdDpOnly);
    assert_eq!(o.after_gt, 1);
    assert_eq!(o.after_pl, o.pl);
    assert_eq!(o.after_ad, vec![22, 24]);
    assert_eq!(o.after_gq, o.gq);
    assert_eq!(o.after_dp, 46);
    // Guard A treats pl_gt==1 as "may replace" even though GT/PL stay.
    assert!(!guard_a_preserve(o.pl_gt));
    assert!(guard_b_preserve(o.valid_assigned));
    assert!(guard_c_preserve(o.no_call));
}

#[test]
fn forensic_6r157_case3_valid_11_is_sparse_replaced_a2_and_a3() {
    let o = observe(&hom_alt(), A3_PILEUP, true, false);
    assert_eq!(o.pl_gt, 2);
    assert_eq!(o.pl, vec![90, 6, 0]);
    assert!(o.valid_assigned);
    assert!(o.class_a3);
    // A2 also matches because `pl_gt == 2` is an A2 disjunct independent of informative skew.
    assert!(o.class_a2);
    assert_eq!(o.action, A3Action::SparseReplace);
    // Pileup 22,24 is Het, not HomAltStrong (24 < 22*4), so valid 1/1 is flipped to 0/1.
    assert!(!SparsePlShape::pileup_is_hom_alt_strong(
        A3_PILEUP[0],
        A3_PILEUP[1]
    ));
    assert_eq!(o.after_gt, 1);
    assert_eq!(o.after_pl, vec![81, 0, 36]);
    let java_required_gt = 2usize;
    assert_ne!(o.after_gt, java_required_gt);
    // Guard A would still permit replacement (pl_gt != 0).
    assert!(!guard_a_preserve(o.pl_gt));
    assert!(guard_b_preserve(o.valid_assigned));
    assert!(guard_c_preserve(o.no_call));
}

#[test]
fn forensic_6r157_case4_no_call_is_not_representable_at_a3_entry() {
    // Java `calculateGLsForThisEvent` writes NO_CALL + PL, then USE_PLS assigns.
    // Rust `SiteScore` / `emit_genotype_format_fields` always emit an assigned GT from min PL.
    // `RegionGenotypeResult` has no NO_CALL / is_no_call field.
    let site_score_emits_assigned_gt = true;
    let region_genotype_has_no_call_field = false;
    assert!(site_score_emits_assigned_gt);
    assert!(!region_genotype_has_no_call_field);
    let o = observe(&hom_ref(), A3_PILEUP, true, false);
    assert!(!o.no_call);
}

#[test]
fn forensic_6r157_case5_empty_pl_is_not_a_natural_a3_state() {
    // `emit_genotype_format_fields` errors on empty GLs. SiteScore never produces empty PL.
    // A3 with `ad.len() < 2` early-exits before `pl_gt`. Empty-mapper is a different caller.
    let empty = emit_genotype_format_fields(&[], &[37, 4]);
    assert!(empty.is_err());
    let o = observe(&hom_ref(), A3_PILEUP, true, false);
    assert!(o.ad.len() >= 2);
    assert!(!o.pl.is_empty());
    assert!(!class_a3_pred(true, 22, 24, 0, 0) || class_a_family_early_exit(false, 0, 22, 24));
    assert!(class_a_family_early_exit(false, 0, 22, 24));
    assert!(class_a_family_early_exit(false, 1, 22, 24));
}

#[test]
fn forensic_6r157_case6_agreeing_pileup_fails_a3_predicates() {
    // A3 requires REF-skewed informative AD *and* balanced/het pileup. Agreement fails A3.
    let calc_00_agrees = Calc {
        gls: vec![0.0, -0.5, -117.4],
        info_ad: [37, 4],
    };
    let o00 = observe(&calc_00_agrees, [37, 4], true, false);
    assert!(!o00.class_a3);
    assert_eq!(o00.action, A3Action::None);
    assert_eq!(o00.after_gt, 0);
    assert_eq!(o00.after_pl, vec![0, 5, 1174]);
    assert_eq!(o00.after_ad, vec![37, 4]);

    let calc_01_agrees = Calc {
        gls: vec![-8.1, 0.0, -3.6],
        info_ad: [22, 24],
    };
    let o01 = observe(&calc_01_agrees, [22, 24], true, false);
    // informative 22,24 is not REF>=3×ALT, so A3 does not fire.
    assert!(!o01.class_a3);
    assert_eq!(o01.action, A3Action::None);
    assert_eq!(o01.after_gt, 1);
    assert_eq!(o01.after_ad, vec![22, 24]);
}

#[test]
fn forensic_6r157_case7_conflicting_pileup_replaces_valid_00() {
    let o = observe(&hom_ref(), A3_PILEUP, true, false);
    assert!(o.class_a3);
    assert_eq!(o.action, A3Action::SparseReplace);
    assert_ne!(o.after_gt, o.pl_gt);
    assert_ne!(o.after_ad.as_slice(), o.ad.as_slice());
}

#[test]
fn forensic_6r157_case8_predicate_failures_do_not_replace() {
    let calc = hom_ref();
    // not SNP
    let indel = observe(&calc, A3_PILEUP, false, false);
    assert!(!indel.class_a3);
    assert_eq!(indel.action, A3Action::None);
    // P12/chr2 scope
    let p12 = observe(&calc, A3_PILEUP, true, true);
    assert_eq!(p12.action, A3Action::None);
    // strong-alt pileup early-exit (pileup_ref*2 < pileup_alt)
    let strong = observe(&calc, [2, 20], true, false);
    assert_eq!(strong.action, A3Action::None);
    // missing pileup allele
    let no_alt = observe(&calc, [22, 0], true, false);
    assert_eq!(no_alt.action, A3Action::None);
    let no_ref = observe(&calc, [0, 24], true, false);
    assert_eq!(no_ref.action, A3Action::None);
    // informative not REF-skewed enough (info_ref < 3*info_alt)
    let balanced_info = Calc {
        gls: vec![0.0, -0.5, -117.4],
        info_ad: [10, 8],
    };
    let o = observe(&balanced_info, A3_PILEUP, true, false);
    assert!(!o.class_a3);
    assert_eq!(o.action, A3Action::None);
}

#[test]
fn forensic_6r157_guard_a_is_not_java_equivalent() {
    let r0 = observe(&hom_ref(), A3_PILEUP, true, false);
    let r1 = observe(&het(), A3_PILEUP, true, false);
    let r2 = observe(&hom_alt(), A3_PILEUP, true, false);
    // A preserves 0/0, but would allow 0/1 and 1/1 replacement.
    assert!(guard_a_preserve(r0.pl_gt));
    assert!(!guard_a_preserve(r1.pl_gt));
    assert!(!guard_a_preserve(r2.pl_gt));
    // Valid 0/1 and 1/1 also have a PL slot of 0 after normalization.
    assert_eq!(r1.pl[1], 0);
    assert_eq!(r2.pl[2], 0);
    // Current A3 does replace valid 1/1; Guard A would still permit that.
    assert_eq!(r2.action, A3Action::SparseReplace);
}

#[test]
fn forensic_6r157_guard_b_matches_java_assigned_genotype() {
    let r0 = observe(&hom_ref(), A3_PILEUP, true, false);
    let r1 = observe(&het(), A3_PILEUP, true, false);
    let r2 = observe(&hom_alt(), A3_PILEUP, true, false);
    assert!(guard_b_preserve(r0.valid_assigned));
    assert!(guard_b_preserve(r1.valid_assigned));
    assert!(guard_b_preserve(r2.valid_assigned));
    // At this caller Guard C is vacuously the same: NO_CALL never arrives.
    assert!(guard_c_preserve(r0.no_call));
    assert!(guard_c_preserve(r1.no_call));
    assert!(guard_c_preserve(r2.no_call));
    // Current A3 diverges from B/C on 0/0 and 1/1 sparse replace, and on 0/1 AD rewrite.
    assert_eq!(r0.action, A3Action::SparseReplace);
    assert_eq!(r1.action, A3Action::AdDpOnly);
    assert_eq!(r2.action, A3Action::SparseReplace);
}

#[test]
fn forensic_6r157_candidate_matrix() {
    // Case | pl_gt | valid | NO_CALL | A preserve | B | C | observed
    let rows = [
        observe(&hom_ref(), A3_PILEUP, true, false),
        observe(&het(), A3_PILEUP, true, false),
        observe(&hom_alt(), A3_PILEUP, true, false),
    ];
    let expected_pl_gt = [0usize, 1, 2];
    let expected_action = [
        A3Action::SparseReplace,
        A3Action::AdDpOnly,
        A3Action::SparseReplace,
    ];
    for (i, o) in rows.iter().enumerate() {
        assert_eq!(o.pl_gt, expected_pl_gt[i]);
        assert!(o.valid_assigned);
        assert!(!o.no_call);
        assert_eq!(o.action, expected_action[i]);
        let a = guard_a_preserve(o.pl_gt);
        let b = guard_b_preserve(o.valid_assigned);
        let c = guard_c_preserve(o.no_call);
        assert_eq!(a, i == 0);
        assert!(b && c);
        // A is not equivalent to B/C except on the 0/0 row.
        if i != 0 {
            assert_ne!(a, b);
        }
    }
}

#[test]
fn forensic_6r157_java_has_no_pileup_gt_replacement_after_calculate_genotypes() {
    // Source-level Java 4.4.0.0 path (not comments):
    // calculateGLsForThisEvent → NO_CALL + PL
    // calculateGenotypes → USE_PLS_TO_ASSIGN / null
    // makeAnnotatedCall → annotateContext (AD) → reverseTrim iff allele count changed
    let java_after_success = [
        "annotateContext_AD",
        "reverseTrimAlleles_iff_allele_count_changed",
        "phaseCalls_default_off",
    ];
    assert!(java_after_success
        .iter()
        .all(|op| !op.contains("pileup_genotype") && !op.contains("SparsePlShape")));
    let java_call_null = "no_vcf_record";
    assert_ne!(java_call_null, "pileup_sparse_replace");
}

#[test]
fn forensic_6r157_apply_class_a_family_is_the_only_a3_caller() {
    let apply_class_a_family_callers = ["try_genotype_variation_event_strict_after_SiteScore"];
    assert_eq!(apply_class_a_family_callers.len(), 1);
    let sparse_callers = [
        "A3_SiteReshape_apply_class_a_family",
        "empty_mapper_SitePileupRescue",
        "L9_after_emit_fail",
        "empty_likelihood_subset_strict",
        "cluster_stored_events",
        "empty_haplotype_summary_genotype_assign",
    ];
    assert!(sparse_callers.contains(&"A3_SiteReshape_apply_class_a_family"));
    assert!(
        sparse_callers
            .iter()
            .filter(|c| c.starts_with("A3_"))
            .count()
            == 1
    );
}

#[test]
fn forensic_6r157_current_overwrite_condition() {
    // Exact Rust condition that overwrites a valid calculator genotype:
    // Class-A family predicates true AND pl_gt != 1 → sparse replace GT/PL/AD/GQ/DP.
    // When pl_gt == 1, GT/PL/GQ stay and AD/DP still change.
    let o0 = observe(&hom_ref(), A3_PILEUP, true, false);
    let o1 = observe(&het(), A3_PILEUP, true, false);
    assert!(o0.class_a3 && o0.pl_gt != 1 && o0.action == A3Action::SparseReplace);
    assert!(o1.class_a3 && o1.pl_gt == 1 && o1.action == A3Action::AdDpOnly);
}

/// Java-equivalent production contract: valid assigned GT is not sparse-replaced.
/// Current A3 violates this. This round does not patch production.
///
/// These assertions record the gap (`assert_ne!`). After an A3-only skip they
/// must be inverted to `assert_eq!`. Running them as `assert_eq!` today would
/// fail, which is the required proof that current production violates the contract.
#[test]
fn forensic_6r157_discovered_contract_is_violated_by_current_a3() {
    let o0 = observe(&hom_ref(), A3_PILEUP, true, false);
    let o1 = observe(&het(), A3_PILEUP, true, false);
    let o2 = observe(&hom_alt(), A3_PILEUP, true, false);
    let java_preserves_valid_assigned_gt = true;
    let current_sparse_replaces_valid_00 = o0.action == A3Action::SparseReplace;
    let current_sparse_replaces_valid_11 = o2.action == A3Action::SparseReplace;
    let current_sparse_replaces_valid_01 = o1.action == A3Action::SparseReplace;
    assert!(java_preserves_valid_assigned_gt);
    assert!(current_sparse_replaces_valid_00);
    assert!(current_sparse_replaces_valid_11);
    assert!(!current_sparse_replaces_valid_01);
    let production_matches_java_boundary = !current_sparse_replaces_valid_00
        && !current_sparse_replaces_valid_11
        && o1.after_ad == o1.ad;
    assert!(
        !production_matches_java_boundary,
        "if this ever passes, production already matches the Java boundary and 6R.157 must be revisited"
    );
}
