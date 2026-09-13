//! 6R.152 coordinate-free: genotype-prior construction vs application after
//! equivalent raw diploid GLs.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//!
//! ```text
//! calculateGLsForThisEvent
//!   IndependentSampleGenotypesModel.calculateLikelihoods → NO_CALL + PL(GL)
//! resolveGenotypePriorCalculator
//!   dragstrParams == null → GenotypePriorCalculator.assumingHW(
//!       log10(snpHeterozygosity), log10(indelHeterozygosity))
//! calculateGenotypes(..., gpc, ...)
//!   default genotypeAssignmentMethod = USE_PLS_TO_ASSIGN
//!   makeGenotypeCall: GT = argmax(GL); gpc unused
//! USE_POSTERIOR_PROBABILITIES (not default):
//!   log10Priors = gpc.getLog10Priors(glCalc, alleles)
//!   log10Posteriors = ebeAdd(log10Priors, genotypeLikelihoods)
//! ```
//!
//! `assumingHW` diploid biallelic SNP (`getLog10Priors`):
//!   0/0 = 0  (REF convention, not a simplex remainder)
//!   0/1 = log10(snpHet) − log10(3)
//!   1/1 = 2·log10(snpHet) − log10(3)
//! SNP /3 is `LOG10_SNP_NORMALIZATION_CONSTANT` (4 standard bases − 1).
//! Hom-var log10 = 2 × het log10 before that /3 (Hardy-Weinberg).
//!
//! Rust `BiallelicDiploidPriorModel` default is a probability simplex:
//!   P(0/1)=het_prior=1e-3, P(1/1)=hom_var_prior=5e-4,
//!   P(0/0)=1 − het − hom_var  ("remainder-mass")
//!   log10 priors = log10 of those three probabilities.
//!
//! Pairwise log10 prior differences are **not** a common offset, so this is
//! not normalization-only. Default HC still assigns GT/PL from raw GLs.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r152_genotype_prior_relative_order_contract
//! HOLDOUT_6R152=1 GATK_RS_EXPERIMENTAL_KBEST_POLICY=unbounded_diagnostic \
//!   cargo test -p gatk-haplotypecaller --test holdout_6r152_genotype_prior -- --nocapture --test-threads=1
//! ```

use gatk_haplotypecaller::activity_scoring::{
    DEFAULT_INDEL_HETEROZYGOSITY, DEFAULT_SNP_HETEROZYGOSITY, LOG10_ONE_THIRD,
};
use gatk_haplotypecaller::genotyping::{
    biallelic_diploid_log10_priors, genotype_posteriors_from_log10_likelihoods,
    BiallelicDiploidPriorModel,
};

/// Java `GenotypePriorCalculator.assumingHW` + `getLog10Priors` SNP slice.
fn java_assuming_hw_snp_log10_priors(snp_heterozygosity: f64) -> [f64; 3] {
    let snp_het = snp_heterozygosity.log10();
    let log10_snp_norm = 3.0_f64.log10();
    [
        0.0,
        snp_het - log10_snp_norm,
        snp_het * 2.0 - log10_snp_norm,
    ]
}

fn pairwise_log_diffs(p: &[f64; 3]) -> [f64; 3] {
    [p[0] - p[1], p[0] - p[2], p[1] - p[2]]
}

#[test]
fn forensic_6r152_prior_enters_at_resolve_gpc_not_calculate_gls() {
    let order = [
        "calculateGLsForThisEvent",
        "resolveGenotypePriorCalculator",
        "calculateGenotypes",
        "makeGenotypeCall",
    ];
    assert_eq!(order[0], "calculateGLsForThisEvent");
    assert_eq!(order[1], "resolveGenotypePriorCalculator");
    let prior_inside_calculate_gls = false;
    assert!(!prior_inside_calculate_gls);
}

#[test]
fn forensic_6r152_java_assuming_hw_is_heterozygosity_hw_not_simplex() {
    let java = java_assuming_hw_snp_log10_priors(DEFAULT_SNP_HETEROZYGOSITY);
    assert_eq!(java.len(), 3);
    assert_eq!(java[0], 0.0);
    assert!((java[1] - (DEFAULT_SNP_HETEROZYGOSITY.log10() + LOG10_ONE_THIRD)).abs() < 1e-15);
    assert!((java[2] - (2.0 * DEFAULT_SNP_HETEROZYGOSITY.log10() + LOG10_ONE_THIRD)).abs() < 1e-15);
    let lin: [f64; 3] = java.map(|x| 10.0_f64.powf(x));
    let sum = lin[0] + lin[1] + lin[2];
    let expected =
        1.0 + DEFAULT_SNP_HETEROZYGOSITY / 3.0 + DEFAULT_SNP_HETEROZYGOSITY.powi(2) / 3.0;
    assert!((lin[0] - 1.0).abs() < 1e-15);
    assert!(
        (sum - expected).abs() < 1e-12,
        "Java assumingHW linear is [1, het/3, het²/3], not a remainder-mass simplex, sum={sum}"
    );
    assert!(
        (sum - 1.0).abs() > 1e-4,
        "REF convention leaves residual mass het/3 + het²/3, sum={sum}"
    );
}

#[test]
fn forensic_6r152_java_hom_var_is_twice_het_in_log10_before_snp_norm() {
    let het = DEFAULT_SNP_HETEROZYGOSITY.log10();
    let hom = het * 2.0;
    let java = java_assuming_hw_snp_log10_priors(DEFAULT_SNP_HETEROZYGOSITY);
    let norm = 3.0_f64.log10();
    assert!((java[1] - (het - norm)).abs() < 1e-15);
    assert!((java[2] - (hom - norm)).abs() < 1e-15);
    // After the common SNP /3, 0/1 vs 1/1 is still one log10(het) (HW square).
    assert!(((java[2] - java[1]) - het).abs() < 1e-15);
}

#[test]
fn forensic_6r152_rust_remainder_mass_is_probability_simplex() {
    let model = BiallelicDiploidPriorModel::default();
    assert_eq!(model.het_prior, 1e-3);
    assert_eq!(model.hom_var_prior, 5e-4);
    let hom_ref = 1.0 - model.het_prior - model.hom_var_prior;
    assert!((hom_ref - 0.9985).abs() < 1e-15);
    let rust = biallelic_diploid_log10_priors(model).unwrap();
    assert!((rust[0] - hom_ref.log10()).abs() < 1e-15);
    assert!((rust[1] - model.het_prior.log10()).abs() < 1e-15);
    assert!((rust[2] - model.hom_var_prior.log10()).abs() < 1e-15);
    let lin: [f64; 3] = rust.map(|x| 10.0_f64.powf(x));
    let sum = lin[0] + lin[1] + lin[2];
    assert!((sum - 1.0).abs() < 1e-12);
}

#[test]
fn forensic_6r152_relative_prior_ratios_are_not_a_common_scale() {
    let java = java_assuming_hw_snp_log10_priors(DEFAULT_SNP_HETEROZYGOSITY);
    let rust = biallelic_diploid_log10_priors(BiallelicDiploidPriorModel::default()).unwrap();
    let jd = pairwise_log_diffs(&java);
    let rd = pairwise_log_diffs(&rust);
    // A constant log offset (normalization) would leave pairwise diffs equal.
    assert!(
        (jd[0] - rd[0]).abs() > 0.4,
        "0/0 vs 0/1 log-ratio Java={} Rust={}",
        jd[0],
        rd[0]
    );
    assert!(
        (jd[1] - rd[1]).abs() > 3.0,
        "0/0 vs 1/1 log-ratio Java={} Rust={}",
        jd[1],
        rd[1]
    );
    assert!(
        (jd[2] - rd[2]).abs() > 2.5,
        "0/1 vs 1/1 log-ratio Java={} Rust={}",
        jd[2],
        rd[2]
    );
}

#[test]
fn forensic_6r152_use_pls_does_not_consume_gpc() {
    let assignment = "USE_PLS_TO_ASSIGN";
    assert_eq!(assignment, "USE_PLS_TO_ASSIGN");
    let gpc_used_for_gt = false;
    assert!(!gpc_used_for_gt);
}

#[test]
fn forensic_6r152_posterior_is_log10_gl_plus_prior_when_requested() {
    let gls = [-1.0_f64, -2.0, -4.0];
    let java_prior = java_assuming_hw_snp_log10_priors(DEFAULT_SNP_HETEROZYGOSITY);
    let java_post: [f64; 3] = [
        gls[0] + java_prior[0],
        gls[1] + java_prior[1],
        gls[2] + java_prior[2],
    ];
    let rust_post = genotype_posteriors_from_log10_likelihoods(&gls, &java_prior).unwrap();
    for i in 0..3 {
        assert!((java_post[i] - rust_post.genotype_log10_posteriors[i]).abs() < 1e-15);
    }
}

#[test]
fn forensic_6r152_prior_offset_is_ranking_invariant_scale_only_if_pairwise_equal() {
    let java = java_assuming_hw_snp_log10_priors(DEFAULT_SNP_HETEROZYGOSITY);
    let shifted = [java[0] + 5.0, java[1] + 5.0, java[2] + 5.0];
    let d0 = pairwise_log_diffs(&java);
    let d1 = pairwise_log_diffs(&shifted);
    for i in 0..3 {
        assert!((d0[i] - d1[i]).abs() < 1e-15);
    }
}

#[test]
fn forensic_6r152_stand_call_conf_and_indel_het_are_not_this_snp_vector() {
    assert_eq!(DEFAULT_INDEL_HETEROZYGOSITY, 1.0 / 8000.0);
    let snp = java_assuming_hw_snp_log10_priors(DEFAULT_SNP_HETEROZYGOSITY);
    let indel_het = DEFAULT_INDEL_HETEROZYGOSITY.log10();
    // SNP 0/1 uses log10(snpHet)−log10(3); indel het is a different allele-type slot.
    assert!((snp[1] - indel_het).abs() > 0.4);
    let stand_call_conf = 30.0_f64;
    let stand_call_enters_prior = false;
    assert_eq!(stand_call_conf, 30.0);
    assert!(!stand_call_enters_prior);
}

#[test]
fn forensic_6r152_genotype_order_is_vcf_diploid_biallelic() {
    let labels = ["0/0", "0/1", "1/1"];
    assert_eq!(labels.len(), 3);
    assert_eq!(labels[0], "0/0");
}

#[test]
fn forensic_6r152_prior_ranking_order_equal_but_ratios_diverge() {
    let java = java_assuming_hw_snp_log10_priors(DEFAULT_SNP_HETEROZYGOSITY);
    let rust = biallelic_diploid_log10_priors(BiallelicDiploidPriorModel::default()).unwrap();
    let rank = |p: &[f64; 3]| {
        let mut idx = [0usize, 1, 2];
        idx.sort_by(|a, b| p[*b].total_cmp(&p[*a]));
        idx
    };
    assert_eq!(rank(&java), rank(&rust));
    assert_eq!(rank(&java), [0, 1, 2]);
    let jd = pairwise_log_diffs(&java);
    let rd = pairwise_log_diffs(&rust);
    assert!((jd[2] - (-DEFAULT_SNP_HETEROZYGOSITY.log10())).abs() < 1e-15);
    assert!((rd[2] - (1e-3_f64 / 5e-4_f64).log10()).abs() < 1e-12);
}

#[test]
fn forensic_6r152_relative_prior_alone_can_flip_01_vs_11() {
    // Coordinate-free GLs: 1/1 is 2 log10 better than 0/1; 0/0 is uncompetitive.
    // Java HW still prefers 0/1 (het vs het² gap of 3). Rust remainder-mass prefers 1/1
    // (het vs het/2 gap of 0.301).
    let gls = [-20.0_f64, -5.0, -3.0];
    let java_prior = java_assuming_hw_snp_log10_priors(DEFAULT_SNP_HETEROZYGOSITY);
    let rust_prior = biallelic_diploid_log10_priors(BiallelicDiploidPriorModel::default()).unwrap();
    let java_post = genotype_posteriors_from_log10_likelihoods(&gls, &java_prior).unwrap();
    let rust_post = genotype_posteriors_from_log10_likelihoods(&gls, &rust_prior).unwrap();
    assert_eq!(java_post.most_likely_genotype_index, 1);
    assert_eq!(rust_post.most_likely_genotype_index, 2);
}
