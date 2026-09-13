//! 6R.159 coordinate-free: hom-ref calculator result → emission eligibility.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//! `HaplotypeCallerGenotypingEngine.assignGenotypeLikelihoods` (default HC VCF):
//!
//! ```text
//! calculateGLsForThisEvent          → NO_CALL + PL
//! USE_PLS_TO_ASSIGN                  → GT = argmin(PL)
//! AlleleFrequencyCalculator.calculate
//! calculateOutputAlleleSubset:
//!   siteIsMonomorphic iff no ALT passesThreshold(stand-call-conf)
//! passesEmitThreshold(QUAL, siteIsMonomorphic):
//!   (EMIT_ALL_CONFIDENT_SITES || !siteIsMonomorphic)
//!   && QUAL >= standardConfidenceForCalling
//! if !passesEmitThreshold && !emitAllActiveSites() → calculateGenotypes = null
//! if (call != null) returnCalls.add(makeAnnotatedCall(...))
//! ```
//!
//! Default HC: `OutputMode.EMIT_VARIANTS_ONLY`, `emitReferenceConfidence=NONE`,
//! `standardConfidenceForCalling=30.0`. GATK 4.4 has **no** separate
//! `standardConfidenceForEmitting` on this path. Sample GQ is not an emit
//! predicate. `bestGuessIsRef` is AF `siteIsMonomorphic`, not sample GT==0/0.
//!
//! Rust production (strict Java, biallelic SNP, not coupled/CTC):
//!
//! ```text
//! SiteScore → SiteReshape (A3 identity) → calculator_is_hom_ref
//!   → GenotypeFinalize::finalize_strict_java_variation_genotype_java
//!     → java_emit_would_pass
//!          1. coupled/CTC special-case (SNP: skip)
//!          2. passes_java_emit_not_hom_ref  // sample GT after PL round-trip
//!          3. passes_hc_variant_emit_biallelic (AF + QUAL vs stand_emit)
//!     → if None: L9 overwrite iff HomAltStrong; else Ok(None)
//! ```
//!
//! Frozen calculator (6R.158, after A3 skip):
//! GT=0/0 PL=0,5,1174 AD=37,4 GQ=5 DP=41 GLs=[0.0,-0.5,-117.4]
//!
//! Production change: NONE.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r159_homref_emission_boundary_contract -- --nocapture --test-threads=1
//! HOLDOUT_6R159=1 GATK_RS_EXPERIMENTAL_KBEST_POLICY=unbounded_diagnostic \
//!   cargo test -p gatk-haplotypecaller --test holdout_6r159_homref_emission -- --nocapture --test-threads=1
//! ```

use gatk_haplotypecaller::emit_gates::{
    java_emit_would_pass, passes_hc_variant_emit_biallelic, passes_java_emit_not_hom_ref,
};
use gatk_haplotypecaller::event_map::VariationEvent;
use gatk_haplotypecaller::genotyping::{
    best_biallelic_diploid_genotype_index, best_pl_index, diploid_genotype_alleles_from_pl_index,
    emit_genotype_format_fields, GenotypeFormatFields, HaplotypeLikelihoodAggregation,
};
use gatk_haplotypecaller::hc_genotyping_engine::{
    java_emit_af_decision, l9_may_overwrite_pairhmm_gls_after_emit_fail, HcGenotypingConfig,
    JavaEmitAfDecision, RegionGenotypeResult, SiteReshape, SparsePlShape,
    DEFAULT_STAND_EMIT_CONFIDENCE,
};
use gatk_haplotypecaller::GenomePosition;

/// Java 4.4 `GenotypeCalculationArgumentCollection.DEFAULT_STANDARD_CONFIDENCE_FOR_CALLING`.
const JAVA_STAND_CALL_CONF: f64 = 30.0;

/// Canonical calculator GLs after Java PL round-trip (`PL=0,5,1174`).
const CANONICAL_GL: [f64; 3] = [0.0, -0.5, -117.4];
const CANONICAL_AD: [i32; 2] = [37, 4];
const A3_PILEUP: [i32; 2] = [22, 24];

/// Test-only trace of the production `java_emit_would_pass` evaluation order.
/// Not a production type.
#[derive(Debug, Clone, PartialEq)]
struct HomRefEmitTrace {
    gt_alleles: Vec<i32>,
    pl: Vec<i32>,
    ad: Vec<i32>,
    gq: i32,
    dp: i32,
    best_genotype_index: usize,
    is_coupled_or_ctc_snp_path: bool,
    passes_java_emit_not_hom_ref: bool,
    gq_below_stand_emit_10: bool,
    gq_gate_evaluated_in_java_emit_would_pass: bool,
    af_at_rust_10: JavaEmitAfDecision,
    af_at_java_30: JavaEmitAfDecision,
    passes_hc_variant_emit_biallelic_at_10: bool,
    passes_hc_variant_emit_biallelic_at_30: bool,
    java_emit_would_pass_at_10: bool,
    java_emit_would_pass_at_30: bool,
    first_failing_at_rust_10: Option<&'static str>,
    first_failing_at_java_30: Option<&'static str>,
    l9_overwrite: bool,
    a3_mutated: bool,
}

fn snp_event() -> VariationEvent {
    VariationEvent {
        contig: "synth".into(),
        start_1based: GenomePosition::new_1based(1),
        end_1based: GenomePosition::new_1based(1),
        ref_allele: "A".into(),
        alt_allele: "T".into(),
    }
}

fn gl_from_pl(pl: &[i32]) -> Vec<f64> {
    pl.iter().map(|&p| (p as f64) / -10.0).collect()
}

fn assigned(gls: &[f64], ad: [i32; 2]) -> RegionGenotypeResult {
    let format = emit_genotype_format_fields(gls, &ad).expect("FORMAT");
    RegionGenotypeResult {
        aggregation: HaplotypeLikelihoodAggregation {
            haplotype_log10_sums: vec![0.0, -1.0],
            read_count: ad.iter().sum::<i32>().max(0) as usize,
        },
        best_haplotype_index: 0,
        ref_haplotype_index: 0,
        alt_haplotype_index: 1,
        genotype_log10_likelihoods: gls.to_vec(),
        format,
    }
}

fn first_failing_java_emit_would_pass(
    event: &VariationEvent,
    gls: &[f64],
    format: &GenotypeFormatFields,
    stand: f64,
) -> (bool, Option<&'static str>) {
    let would = java_emit_would_pass(event, gls, format, stand, &[]).expect("emit");
    if would {
        return (true, None);
    }
    // Production order inside java_emit_would_pass for a biallelic SNP with empty
    // region_events: coupled/CTC skip → passes_java_emit_not_hom_ref → AF.
    if !passes_java_emit_not_hom_ref(gls, format) {
        return (false, Some("passes_java_emit_not_hom_ref"));
    }
    if !passes_hc_variant_emit_biallelic(gls, stand).expect("af") {
        return (false, Some("passes_hc_variant_emit_biallelic"));
    }
    (false, Some("java_emit_would_pass_other"))
}

fn trace_case(name: &str, gls: &[f64], ad: [i32; 2]) -> HomRefEmitTrace {
    trace_case_with_a3(name, gls, ad, false)
}

fn trace_canonical_through_a3(name: &str, gls: &[f64], ad: [i32; 2]) -> HomRefEmitTrace {
    trace_case_with_a3(name, gls, ad, true)
}

fn trace_case_with_a3(name: &str, gls: &[f64], ad: [i32; 2], apply_a3: bool) -> HomRefEmitTrace {
    let event = snp_event();
    let incoming = assigned(gls, ad);
    let gt = if apply_a3 {
        SiteReshape::apply_class_a_family(
            incoming.clone(),
            &event,
            A3_PILEUP[0],
            A3_PILEUP[1],
            &HcGenotypingConfig::strict_java(),
        )
        .expect("A3")
    } else {
        incoming.clone()
    };
    let a3_mutated = apply_a3
        && (gt.format != incoming.format
            || gt.genotype_log10_likelihoods != incoming.genotype_log10_likelihoods);
    let pl = gt.format.pl_as_i32();
    let ad_v = gt.format.ad_as_i32();
    let gq = gt.format.gq.as_i32();
    let dp = gt.format.dp.as_i32();
    let best = best_biallelic_diploid_genotype_index(&gt.genotype_log10_likelihoods, &ad_v);
    let gt_alleles = diploid_genotype_alleles_from_pl_index(2, best_pl_index(&gt.format.pl));
    let not_hom_ref = passes_java_emit_not_hom_ref(&gt.genotype_log10_likelihoods, &gt.format);
    let af10 = java_emit_af_decision(
        &gt.genotype_log10_likelihoods,
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("af10");
    let af30 =
        java_emit_af_decision(&gt.genotype_log10_likelihoods, JAVA_STAND_CALL_CONF).expect("af30");
    let af_emit_10 = passes_hc_variant_emit_biallelic(
        &gt.genotype_log10_likelihoods,
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("af emit 10");
    let af_emit_30 =
        passes_hc_variant_emit_biallelic(&gt.genotype_log10_likelihoods, JAVA_STAND_CALL_CONF)
            .expect("af emit 30");
    let (pass10, fail10) = first_failing_java_emit_would_pass(
        &event,
        &gt.genotype_log10_likelihoods,
        &gt.format,
        DEFAULT_STAND_EMIT_CONFIDENCE,
    );
    let (pass30, fail30) = first_failing_java_emit_would_pass(
        &event,
        &gt.genotype_log10_likelihoods,
        &gt.format,
        JAVA_STAND_CALL_CONF,
    );
    let l9 = l9_may_overwrite_pairhmm_gls_after_emit_fail(
        &event,
        A3_PILEUP[0],
        A3_PILEUP[1],
        best_pl_index(&gt.format.pl) == 0,
    );
    eprintln!(
        "6R159\tcase={name}\tGT={}/{} PL={:?} AD={:?} GQ={gq} DP={dp} best={best} \
         not_hom_ref={not_hom_ref} af10_mono={} af10_plausible={} af10_phred={:.4} af10_emit={} \
         af30_mono={} af30_plausible={} af30_phred={:.4} af30_emit={} \
         java_emit_10={pass10} fail10={fail10:?} java_emit_30={pass30} fail30={fail30:?} \
         gq_lt_10={} gq_in_java_emit_would_pass=false a3_mutated={a3_mutated} l9={l9}",
        gt_alleles[0],
        gt_alleles[1],
        pl,
        ad_v,
        af10.site_is_monomorphic,
        af10.alt_plausible,
        af10.phred_scaled,
        af10.passes_emit,
        af30.site_is_monomorphic,
        af30.alt_plausible,
        af30.phred_scaled,
        af30.passes_emit,
        (gq as f64) < DEFAULT_STAND_EMIT_CONFIDENCE,
    );
    HomRefEmitTrace {
        gt_alleles,
        pl,
        ad: ad_v,
        gq,
        dp,
        best_genotype_index: best,
        is_coupled_or_ctc_snp_path: false,
        passes_java_emit_not_hom_ref: not_hom_ref,
        gq_below_stand_emit_10: (gq as f64) < DEFAULT_STAND_EMIT_CONFIDENCE,
        gq_gate_evaluated_in_java_emit_would_pass: false,
        af_at_rust_10: af10,
        af_at_java_30: af30,
        passes_hc_variant_emit_biallelic_at_10: af_emit_10,
        passes_hc_variant_emit_biallelic_at_30: af_emit_30,
        java_emit_would_pass_at_10: pass10,
        java_emit_would_pass_at_30: pass30,
        first_failing_at_rust_10: fail10,
        first_failing_at_java_30: fail30,
        l9_overwrite: l9,
        a3_mutated,
    }
}

#[test]
fn forensic_6r159_canonical_calculator_reaches_emit_decision_unchanged() {
    let t = trace_canonical_through_a3("A_canonical_00_gq5", &CANONICAL_GL, CANONICAL_AD);
    assert_eq!(t.gt_alleles, vec![0, 0]);
    assert_eq!(t.pl, vec![0, 5, 1174]);
    assert_eq!(t.ad, vec![37, 4]);
    assert_eq!(t.gq, 5);
    assert_eq!(t.dp, 41);
    assert_eq!(t.best_genotype_index, 0);
    assert!(!t.a3_mutated, "6R.158 A3 skip must still hold");
    assert!(!t.l9_overwrite, "pileup 22,24 is Het, not HomAltStrong");
    assert!(!t.is_coupled_or_ctc_snp_path);
    assert!(!t.gq_gate_evaluated_in_java_emit_would_pass);
}

#[test]
fn forensic_6r159_first_failing_predicate_is_passes_java_emit_not_hom_ref() {
    let t = trace_canonical_through_a3("A_canonical_00_gq5", &CANONICAL_GL, CANONICAL_AD);
    assert!(
        !t.passes_java_emit_not_hom_ref,
        "canonical best genotype index is 0"
    );
    assert_eq!(
        t.first_failing_at_rust_10,
        Some("passes_java_emit_not_hom_ref")
    );
    assert_eq!(
        t.first_failing_at_java_30,
        Some("passes_java_emit_not_hom_ref")
    );
    assert!(!t.java_emit_would_pass_at_10);
    assert!(!t.java_emit_would_pass_at_30);
    // AF is not consulted after the GT veto. Record it; do not treat it as the first fail.
    assert_eq!(
        t.passes_hc_variant_emit_biallelic_at_10,
        t.af_at_rust_10.passes_emit
    );
    assert_eq!(
        t.passes_hc_variant_emit_biallelic_at_30,
        t.af_at_java_30.passes_emit
    );
}

#[test]
fn forensic_6r159_gq5_is_not_the_first_emit_predicate() {
    // Case E: same 0/0 genotype with GQ well above Rust stand_emit=10.
    let high_gq_00 = gl_from_pl(&[0, 40, 1000]);
    let e = trace_case("E_00_gq40", &high_gq_00, CANONICAL_AD);
    assert_eq!(e.gt_alleles, vec![0, 0]);
    assert_eq!(e.gq, 40);
    assert!(!e.gq_below_stand_emit_10);
    assert!(!e.passes_java_emit_not_hom_ref);
    assert_eq!(
        e.first_failing_at_rust_10,
        Some("passes_java_emit_not_hom_ref")
    );
    assert!(
        !e.java_emit_would_pass_at_10,
        "raising GQ on 0/0 does not change the first Rust emit predicate"
    );

    // Case A still has GQ=5, but that gate is not inside java_emit_would_pass.
    let a = trace_canonical_through_a3("A_canonical_00_gq5", &CANONICAL_GL, CANONICAL_AD);
    assert!(a.gq_below_stand_emit_10);
    assert!(!a.gq_gate_evaluated_in_java_emit_would_pass);
}

#[test]
fn forensic_6r159_non_homref_low_gq_is_not_dropped_by_the_gt_veto() {
    // Case C: valid 0/1 with the same GQ=5 as the canonical site.
    let het_low_gq = gl_from_pl(&[5, 0, 1174]);
    let c = trace_case("C_01_gq5", &het_low_gq, CANONICAL_AD);
    assert_eq!(c.gt_alleles, vec![0, 1]);
    assert_eq!(c.gq, 5);
    assert!(c.passes_java_emit_not_hom_ref);
    assert_ne!(
        c.first_failing_at_rust_10,
        Some("passes_java_emit_not_hom_ref")
    );

    // Case D: valid 1/1 with GQ=5. Do not run A3 — that family is not this round's
    // emit predicate and would reshape AD [0,4] against the canonical pileup.
    let hom_alt_low_gq = gl_from_pl(&[90, 5, 0]);
    let d = trace_case("D_11_gq5", &hom_alt_low_gq, [0, 4]);
    assert_eq!(d.gt_alleles, vec![1, 1]);
    assert_eq!(d.gq, 5);
    assert!(!d.a3_mutated);
    assert!(d.passes_java_emit_not_hom_ref);
    assert_ne!(
        d.first_failing_at_rust_10,
        Some("passes_java_emit_not_hom_ref")
    );
}

#[test]
fn forensic_6r159_java_emit_is_af_site_monomorphic_not_sample_gt() {
    let a = trace_canonical_through_a3("A_canonical_00_gq5", &CANONICAL_GL, CANONICAL_AD);
    // Java 4.4 default EMIT_VARIANTS_ONLY: hom-ref sample GT can still be emitted
    // iff the AF calculator keeps an ALT (siteIsMonomorphic=false) and QUAL>=30.
    // For this calculator vector, record the AF truth and the Java null-call rule.
    let java_would_return_non_null_at_30 = a.af_at_java_30.passes_emit;
    let java_would_return_non_null_at_10 = a.af_at_rust_10.passes_emit;
    eprintln!(
        "6R159\tjava_calculateGenotypes_non_null_at_30={java_would_return_non_null_at_30} \
         java_calculateGenotypes_non_null_at_10={java_would_return_non_null_at_10} \
         java_siteIsMonomorphic_30={} java_siteIsMonomorphic_10={}",
        a.af_at_java_30.site_is_monomorphic, a.af_at_rust_10.site_is_monomorphic
    );
    assert!(
        a.af_at_java_30.site_is_monomorphic,
        "canonical PL 0,5,1174 is AF-monomorphic at Java stand-call-conf=30"
    );
    assert!(
        a.af_at_java_30.phred_scaled >= JAVA_STAND_CALL_CONF,
        "monomorphic-site QUAL can still be >=30; Java still drops because bestGuessIsRef"
    );
    assert!(
        !java_would_return_non_null_at_30,
        "Java calculateGenotypes is null at default stand-call-conf=30"
    );
    assert!(
        a.af_at_rust_10.site_is_monomorphic,
        "canonical PL 0,5,1174 is also AF-monomorphic at Rust stand_emit=10"
    );
    assert!(!java_would_return_non_null_at_10);
}

#[test]
fn forensic_6r159_config_calling_vs_emitting() {
    assert_eq!(DEFAULT_STAND_EMIT_CONFIDENCE, 10.0);
    assert_eq!(JAVA_STAND_CALL_CONF, 30.0);
    let cfg = HcGenotypingConfig::strict_java();
    assert_eq!(cfg.stand_emit_confidence, 10.0);
    // GQ=5 vs 10 is a later/unreached Rust gate, not the first drop.
    let a = trace_canonical_through_a3("A_canonical_00_gq5", &CANONICAL_GL, CANONICAL_AD);
    assert!(a.gq_below_stand_emit_10);
    assert_eq!(
        a.first_failing_at_rust_10,
        Some("passes_java_emit_not_hom_ref")
    );
}

#[test]
fn forensic_6r159_no_call_not_representable_at_finalize() {
    // Case F: SiteScore always assigns via USE_PLS_TO_ASSIGN. RegionGenotypeResult
    // has no NO_CALL allele vector at GenotypeFinalize.
    let t = trace_canonical_through_a3("A_canonical_00_gq5", &CANONICAL_GL, CANONICAL_AD);
    assert_ne!(t.gt_alleles, vec![-1, -1]);
    assert!(!t.pl.is_empty());
}

#[test]
fn forensic_6r159_case_b_same_gt_higher_gq_still_gt_veto() {
    // Case B: same 0/0 PL shape with a larger PL gap (GQ=40 from PL 0,40,*).
    let b = trace_case(
        "B_00_gq40_same_shape",
        &gl_from_pl(&[0, 40, 1174]),
        CANONICAL_AD,
    );
    assert_eq!(b.gt_alleles, vec![0, 0]);
    assert_eq!(b.gq, 40);
    assert_eq!(
        b.first_failing_at_rust_10,
        Some("passes_java_emit_not_hom_ref")
    );
    assert!(!b.java_emit_would_pass_at_10);
}

#[test]
fn forensic_6r159_het_emits_when_af_passes() {
    // Existing 6R.105 het vector: PL 81,0,36 emits at 10 and 30. Proves the
    // machinery is not a universal GQ rule.
    let het = trace_case("het_81_0_36", &gl_from_pl(&[81, 0, 36]), [44, 4]);
    assert_eq!(het.gt_alleles, vec![0, 1]);
    assert!(het.passes_java_emit_not_hom_ref);
    assert!(het.java_emit_would_pass_at_10);
    assert!(het.java_emit_would_pass_at_30);
    assert!(het.first_failing_at_rust_10.is_none());
}

#[test]
fn forensic_6r159_sparse_shape_het_is_not_l9() {
    assert_eq!(
        SparsePlShape::from_pileup_depths(A3_PILEUP[0], A3_PILEUP[1]),
        SparsePlShape::Het
    );
    assert!(!SparsePlShape::pileup_is_hom_alt_strong(
        A3_PILEUP[0],
        A3_PILEUP[1]
    ));
}
