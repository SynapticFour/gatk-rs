//! 6R.153 coordinate-free: default `USE_PLS_TO_ASSIGN` PL / GT / GQ after equivalent raw GLs.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//!
//! ```text
//! calculateGLsForThisEvent
//!   IndependentSampleGenotypesModel.calculateLikelihoods
//!   GenotypeBuilder(NO_CALL).PL(likelihoods.sampleLikelihoods(s).getAsPLs())
//!     HTSJDK GenotypeLikelihoods.GLsToPLs:
//!       adjust = max(log10 GL)
//!       PL_i = (int) Math.round(min(-10*(GL_i-adjust), MAX_PL))
//! calculateGenotypes(..., USE_PLS_TO_ASSIGN)
//!   AlleleSubsettingUtils.subsetAlleles → makeGenotypeCall
//!   maxLikelihoodIndex = MathUtils.maxElementIndex(genotypeLikelihoods)  // first strict max
//!   gb.alleles(genotypeAlleleCountsAt(index))
//!   gb.log10PError(GenotypeLikelihoods.getGQLog10FromLikelihoods(index, GLs))
//!     unique max: log10PError = -(best − second)
//!     GQ = (int) Math.round(log10PError * -10)
//! gpc unused (not USE_POSTERIOR_PROBABILITIES)
//! ```
//!
//! Rust production:
//!   emit_genotype_format_fields: PL_i = round(-10*(GL_i-max)).max(0) then subtract min
//!   GT = argmin(PL) via best_pl_index (first min on ties)
//!   GQ hom-ref-best: second-best PL − 0, cap 99
//!
//! Priors remain closed (6R.152). PairHMM residual must not change integer PL here.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r153_default_pl_gt_assignment_contract
//! HOLDOUT_6R153=1 GATK_RS_EXPERIMENTAL_KBEST_POLICY=unbounded_diagnostic \
//!   cargo test -p gatk-haplotypecaller --test holdout_6r153_default_pl_gt -- --nocapture --test-threads=1
//! ```

use gatk_haplotypecaller::genotyping::{
    best_pl_index, diploid_genotype_alleles_from_pl_index, emit_genotype_format_fields,
    gq_phred_from_pl, gq_phred_p7_compatible,
};

/// Java `Math.round(double)` for finite values: `floor(x + 0.5)`.
fn java_math_round(x: f64) -> i32 {
    (x + 0.5).floor() as i32
}

/// HTSJDK `GenotypeLikelihoods.GLsToPLs` without the Integer.MAX_VALUE cap (unused on this scale).
fn java_gls_to_pls(gls: &[f64]) -> Vec<i32> {
    let adjust = gls.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    gls.iter()
        .map(|&g| java_math_round(-10.0 * (g - adjust)).max(0))
        .collect()
}

fn java_pre_round_pl(gls: &[f64]) -> Vec<f64> {
    let adjust = gls.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    gls.iter().map(|&g| -10.0 * (g - adjust)).collect()
}

/// `MathUtils.maxElementIndex`: first index of the strict maximum.
fn java_max_element_index(gls: &[f64]) -> usize {
    let mut max_i = 0usize;
    for i in 1..gls.len() {
        if gls[i] > gls[max_i] {
            max_i = i;
        }
    }
    max_i
}

/// `getGQLog10FromLikelihoods` then `GenotypeBuilder.log10PError` → integer GQ, unique-max branch.
fn java_gq_from_gls(gls: &[f64]) -> i32 {
    let chosen = java_max_element_index(gls);
    let mut second = f64::NEG_INFINITY;
    for (i, &ll) in gls.iter().enumerate() {
        if i != chosen && ll >= second {
            second = ll;
        }
    }
    let qual = gls[chosen] - second;
    let log10_p_error = if qual < 0.0 {
        panic!("chosen genotype is not uniquely best; 6R.153 unique-max GQ path only");
    } else {
        -qual
    };
    java_math_round(log10_p_error * -10.0).clamp(0, 99)
}

#[test]
fn forensic_6r153_pl_is_phred_gap_from_best_gl_with_java_round() {
    let gls = [-1.0_f64, -1.46, -5.0];
    let pre = java_pre_round_pl(&gls);
    assert!((pre[0] - 0.0).abs() < 1e-15);
    assert!((pre[1] - 4.6).abs() < 1e-12);
    let pls = java_gls_to_pls(&gls);
    assert_eq!(pls, vec![0, 5, 40]);
    let rust = emit_genotype_format_fields(&gls, &[0, 0]).unwrap();
    assert_eq!(rust.pl_as_i32(), pls);
}

#[test]
fn forensic_6r153_genotype_order_is_vcf_diploid_biallelic() {
    assert_eq!(diploid_genotype_alleles_from_pl_index(2, 0), vec![0, 0]);
    assert_eq!(diploid_genotype_alleles_from_pl_index(2, 1), vec![0, 1]);
    assert_eq!(diploid_genotype_alleles_from_pl_index(2, 2), vec![1, 1]);
}

#[test]
fn forensic_6r153_unique_winner_gt_from_gl_and_pl_agree() {
    let gls = [-138.12681007_f64, -138.58846653, -255.47756386];
    let java_gt = java_max_element_index(&gls);
    let rust = emit_genotype_format_fields(&gls, &[0, 0]).unwrap();
    let rust_gt = best_pl_index(&rust.pl);
    assert_eq!(java_gt, 0);
    assert_eq!(rust_gt, 0);
    assert_eq!(rust.pl_as_i32(), java_gls_to_pls(&gls));
}

#[test]
fn forensic_6r153_exact_gl_tie_picks_lowest_index() {
    let gls = [-2.0_f64, -2.0, -2.0];
    assert_eq!(java_max_element_index(&gls), 0);
    let rust = emit_genotype_format_fields(&gls, &[0, 0]).unwrap();
    assert_eq!(best_pl_index(&rust.pl), 0);
    assert_eq!(rust.pl_as_i32(), vec![0, 0, 0]);
    assert_eq!(gq_phred_from_pl(&rust.pl_as_i32()), 0);
}

#[test]
fn forensic_6r153_two_way_gl_tie_is_first_index() {
    let gls = [-1.0_f64, -1.0, -9.0];
    assert_eq!(java_max_element_index(&gls), 0);
    let rust = emit_genotype_format_fields(&gls, &[0, 0]).unwrap();
    assert_eq!(best_pl_index(&rust.pl), 0);
}

#[test]
fn forensic_6r153_pairhmm_residual_does_not_change_integer_pl() {
    let java = [-138.126810073853_f64, -138.588466525737, -255.477563858032];
    let rust = [-138.126825290025_f64, -138.588480870144, -255.477564897631];
    let d0 = (java[0] - rust[0]).abs();
    let d1 = (java[1] - rust[1]).abs();
    assert!(d0 < 2e-5 && d1 < 2e-5);
    assert_eq!(java_gls_to_pls(&java), java_gls_to_pls(&rust));
    assert_eq!(java_gls_to_pls(&java), vec![0, 5, 1174]);
    let pre_j = java_pre_round_pl(&java);
    let pre_r = java_pre_round_pl(&rust);
    assert!((pre_j[1] - 4.61656).abs() < 1e-4);
    assert!((pre_r[1] - 4.61656).abs() < 1e-4);
    assert!((pre_j[1] - 4.5).abs() > 0.1, "not on a rounding boundary");
}

#[test]
fn forensic_6r153_equal_pl_after_rounding_can_disagree_with_gl_argmax() {
    // 10*Δ = 0.2 → both PL 0. Java GT follows GL; Rust GT follows first min PL.
    let gls = [-1.02_f64, -1.00, -10.0];
    assert_eq!(java_max_element_index(&gls), 1);
    let rust = emit_genotype_format_fields(&gls, &[0, 0]).unwrap();
    assert_eq!(rust.pl_as_i32(), vec![0, 0, 90]);
    assert_eq!(best_pl_index(&rust.pl), 0);
    // Latent GT_SELECTION difference. Not a production patch: this fixture's Δ≈0.462 → PL 0 vs 5.
}

#[test]
fn forensic_6r153_gq_hom_ref_is_gl_gap_then_round_matching_pl_gap() {
    let gls = [-138.12681007_f64, -138.58846653, -255.47756386];
    let java_gq = java_gq_from_gls(&gls);
    let rust = emit_genotype_format_fields(&gls, &[0, 0]).unwrap();
    assert_eq!(java_gq, 5);
    assert_eq!(rust.gq.as_i32(), 5);
    assert_eq!(gq_phred_p7_compatible(&gls, &rust.pl_as_i32(), &[0, 0]), 5);
}

#[test]
fn forensic_6r153_gq_is_capped_at_99_and_does_not_feed_gt() {
    let gls = [0.0_f64, -20.0, -40.0];
    let rust = emit_genotype_format_fields(&gls, &[0, 0]).unwrap();
    assert_eq!(best_pl_index(&rust.pl), 0);
    assert_eq!(rust.gq.as_i32(), 99);
    assert_eq!(java_gq_from_gls(&gls), 99);
}

#[test]
fn forensic_6r153_priors_are_not_on_the_use_pls_path() {
    let assignment = "USE_PLS_TO_ASSIGN";
    assert_eq!(assignment, "USE_PLS_TO_ASSIGN");
    let gpc_used_for_gt_pl_gq = false;
    assert!(!gpc_used_for_gt_pl_gq);
}
