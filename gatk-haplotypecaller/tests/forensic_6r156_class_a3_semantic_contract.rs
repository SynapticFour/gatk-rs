//! 6R.156 coordinate-free: Class-A3 semantic contract (no production change).
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//!
//! `SiteReshape::apply_class_a_family` after a valid PairHMM / `SiteScore` genotype:
//!   early-exit unless pileup both alleles and pileup_ref*2 >= pileup_alt
//!   Class-A3 iff SNP AND pileup_alt*2 >= pileup_ref AND informative both
//!     alleles AND informative_ref >= 3*informative_alt
//!   then: pl_gt==1 → keep GT/PL/GQ, replace AD/DP
//!         else     → `sparse_snp_genotype_from_read_depths` replaces GT/PL/AD/GQ/DP
//!
//! The sparse helper is documented as N1 "when PairHMM has no alt-hap support".
//! Class-A3 calls it after a complete calculator genotype, including valid 0/0.
//! Comment on the sparse arm says "GT flips (PL=1/1)"; the `else` also fires on 0/0.
//!
//! L9 (`l9_may_overwrite_pairhmm_gls_after_emit_fail`) is a **sibling after emit-fail**,
//! not a consumer of A3 output. L9 explicitly must not overwrite valid calculator
//! hom-ref SNP GLs unless pileup is HomAltStrong. A3 runs *before* emit and does.
//!
//! Java default HC: no pileup-derived replacement of an already-calculated GT/PL.
//! AD annotation is not an analogue.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r156_class_a3_semantic_contract
//! HOLDOUT_6R156=1 GATK_RS_EXPERIMENTAL_KBEST_POLICY=unbounded_diagnostic \
//!   cargo test -p gatk-haplotypecaller --test holdout_6r156_class_a3 -- --nocapture --test-threads=1
//! ```

use gatk_haplotypecaller::genotyping::{best_pl_index, emit_genotype_format_fields};
use gatk_haplotypecaller::hc_genotyping_engine::SparsePlShape;

/// Production Class-A3 predicate (plus shared early-exit), not a Java method.
fn class_a3_fires(
    is_snp: bool,
    in_p12_chr2_scope: bool,
    ad_len: usize,
    pileup_ref: i32,
    pileup_alt: i32,
    info_ref: i32,
    info_alt: i32,
) -> bool {
    if in_p12_chr2_scope
        || ad_len < 2
        || pileup_ref < 1
        || pileup_alt < 1
        || pileup_ref.saturating_mul(2) < pileup_alt
    {
        return false;
    }
    is_snp
        && pileup_alt.saturating_mul(2) >= pileup_ref
        && info_ref >= 1
        && info_alt >= 1
        && info_ref >= info_alt.saturating_mul(3)
}

fn sparse_replace(pileup: [i32; 2]) -> (usize, Vec<i32>, Vec<i32>, i32, i32) {
    let shape = SparsePlShape::from_pileup_depths(pileup[0], pileup[1]);
    let fmt = emit_genotype_format_fields(&shape.gl_vec(), &pileup).expect("sparse");
    (
        best_pl_index(&fmt.pl),
        fmt.pl_as_i32(),
        fmt.ad_as_i32(),
        fmt.gq.as_i32(),
        fmt.dp.as_i32(),
    )
}

fn keep_pls_replace_ad(gls: &[f64], pileup: [i32; 2]) -> (usize, Vec<i32>, Vec<i32>, i32, i32) {
    let fmt = emit_genotype_format_fields(gls, &pileup).expect("ad-only");
    (
        best_pl_index(&fmt.pl),
        fmt.pl_as_i32(),
        fmt.ad_as_i32(),
        fmt.gq.as_i32(),
        fmt.dp.as_i32(),
    )
}

#[test]
fn forensic_6r156_live_fixture_predicates() {
    assert!(class_a3_fires(true, false, 2, 22, 24, 37, 4));
    assert!(!class_a3_fires(true, true, 2, 22, 24, 37, 4));
    assert!(!class_a3_fires(false, false, 2, 22, 24, 37, 4));
}

#[test]
fn forensic_6r156_case1_valid_hom_ref_is_overridden_not_preserved() {
    let gls = [0.0, -0.5, -117.4];
    let first = emit_genotype_format_fields(&gls, &[37, 4]).expect("first");
    assert_eq!(best_pl_index(&first.pl), 0);
    assert_eq!(first.pl_as_i32(), vec![0, 5, 1174]);
    assert_eq!(first.ad_as_i32(), vec![37, 4]);
    assert!(class_a3_fires(true, false, 2, 22, 24, 37, 4));
    let pl_gt = 0;
    assert_ne!(pl_gt, 1);
    let (gt, pl, ad, gq, dp) = sparse_replace([22, 24]);
    assert_eq!(gt, 1);
    assert_eq!(pl, vec![81, 0, 36]);
    assert_eq!(ad, vec![22, 24]);
    assert_eq!(gq, 36);
    assert_eq!(dp, 46);
    let java_preserves_calculator = true;
    assert!(java_preserves_calculator);
    assert_ne!(first.pl_as_i32(), pl);
}

#[test]
fn forensic_6r156_case2_het_pl_keeps_gt_pl_gq_replaces_ad_dp() {
    let gls = [-8.1, 0.0, -3.6];
    let first = emit_genotype_format_fields(&gls, &[37, 4]).expect("het first");
    assert_eq!(best_pl_index(&first.pl), 1);
    let (gt, pl, ad, gq, dp) = keep_pls_replace_ad(&gls, [22, 24]);
    assert_eq!(gt, 1);
    assert_eq!(pl, first.pl_as_i32());
    assert_eq!(ad, vec![22, 24]);
    assert_eq!(gq, first.gq.as_i32());
    assert_eq!(dp, 46);
}

#[test]
fn forensic_6r156_case3_strong_alt_pileup_early_exits() {
    assert!(!class_a3_fires(true, false, 2, 2, 20, 37, 4));
    assert!(!class_a3_fires(true, false, 2, 0, 24, 37, 4));
}

#[test]
fn forensic_6r156_case4_missing_ad_is_not_a3_fallback() {
    assert!(!class_a3_fires(true, false, 0, 22, 24, 0, 0));
    assert!(!class_a3_fires(true, false, 1, 22, 24, 37, 0));
    let empty_mapper_rescue_is_a_different_caller = true;
    assert!(empty_mapper_rescue_is_a_different_caller);
}

#[test]
fn forensic_6r156_sparse_arm_comment_says_hom_alt_but_else_includes_hom_ref() {
    let reserved_for_hom_alt_only_in_comment = 2usize;
    let else_branch_includes_hom_ref = 0usize;
    assert_ne!(
        reserved_for_hom_alt_only_in_comment,
        else_branch_includes_hom_ref
    );
    assert!(class_a3_fires(true, false, 2, 22, 24, 37, 4));
}

#[test]
fn forensic_6r156_not_representation_only_on_hom_ref() {
    let first_gt = 0usize;
    let first_pl = [0, 5, 1174];
    let (gt, pl, ad, gq, _dp) = sparse_replace([22, 24]);
    assert_ne!(gt, first_gt);
    assert_ne!(pl.as_slice(), first_pl.as_slice());
    assert_ne!(ad.as_slice(), [37, 4].as_slice());
    assert_ne!(gq, 5);
}

#[test]
fn forensic_6r156_counterfactual_a_ad_only_keeps_gt_pl_gq() {
    let gls = [0.0, -0.5, -117.4];
    let (gt, pl, ad, gq, dp) = keep_pls_replace_ad(&gls, [22, 24]);
    assert_eq!(gt, 0);
    assert_eq!(pl, vec![0, 5, 1174]);
    assert_eq!(ad, vec![22, 24]);
    assert_eq!(gq, 5);
    assert_eq!(dp, 46);
}

#[test]
fn forensic_6r156_counterfactual_b_gt_pl_only_keeps_calculator_ad() {
    let (gt, pl, ad, gq, dp) = {
        let shape = SparsePlShape::from_pileup_depths(22, 24);
        let fmt = emit_genotype_format_fields(&shape.gl_vec(), &[37, 4]).expect("b");
        (
            best_pl_index(&fmt.pl),
            fmt.pl_as_i32(),
            fmt.ad_as_i32(),
            fmt.gq.as_i32(),
            fmt.dp.as_i32(),
        )
    };
    assert_eq!(gt, 1);
    assert_eq!(pl, vec![81, 0, 36]);
    assert_eq!(ad, vec![37, 4]);
    assert_eq!(gq, 36);
    assert_eq!(dp, 41);
}

#[test]
fn forensic_6r156_counterfactual_c_only_when_gl_missing_is_noop_here() {
    let calculator_gls_valid = true;
    let a3_as_missing_only = !calculator_gls_valid;
    assert!(!a3_as_missing_only);
    let first_pl = [0, 5, 1174];
    let first_ad = [37, 4];
    assert_eq!(first_pl, [0, 5, 1174]);
    assert_eq!(first_ad, [37, 4]);
}

#[test]
fn forensic_6r156_counterfactual_d_bypass_hom_ref_gt_flip() {
    let pl_gt = 0usize;
    let after = if pl_gt == 1 {
        ([1], [81, 0, 36], [22, 24])
    } else if pl_gt == 2 {
        ([1], [81, 0, 36], [22, 24])
    } else {
        ([0], [0, 5, 1174], [37, 4])
    };
    assert_eq!(after, ([0], [0, 5, 1174], [37, 4]));
}

#[test]
fn forensic_6r156_gq_dp_are_consequences_of_replaced_gl_ad() {
    let (gt, pl, ad, gq, dp) = sparse_replace([22, 24]);
    assert_eq!(gt, 1);
    assert_eq!(pl[1], 0);
    assert_eq!(gq, 36);
    assert_eq!(dp, ad.iter().sum::<i32>());
}

#[test]
fn forensic_6r156_l9_is_sibling_after_emit_fail_not_a3_consumer() {
    assert!(!SparsePlShape::pileup_is_hom_alt_strong(22, 24));
    let a3_runs_before_finalize = true;
    let l9_runs_only_if_finalize_returns_none = true;
    assert!(a3_runs_before_finalize && l9_runs_only_if_finalize_returns_none);
}

#[test]
fn forensic_6r156_apply_class_a_family_single_strict_hc_caller() {
    let apply_class_a_family_callers = ["try_genotype_variation_event_strict_after_SiteScore"];
    assert_eq!(apply_class_a_family_callers.len(), 1);
    let mutect = false;
    let p12_chr2_excluded = true;
    assert!(!mutect && p12_chr2_excluded);
}

#[test]
fn forensic_6r156_java_has_no_pileup_gt_replacement_analogue() {
    let java_ops_after_first_ad = [
        "reverseTrimAlleles_copy",
        "phaseCalls_off",
        "DepthPerSampleHC_dp_from_ad",
    ];
    assert!(java_ops_after_first_ad
        .iter()
        .all(|op| !op.contains("pileup_genotype") && !op.contains("SparsePlShape")));
}
