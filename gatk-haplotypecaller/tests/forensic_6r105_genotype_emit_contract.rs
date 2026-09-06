//! 6R.105 coordinate-free: genotype assignment / emit / output-allele subset.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//! `GenotypingEngine.calculateGenotypes` (default HC VCF, `EMIT_VARIANTS_ONLY`):
//!
//! ```text
//! calculateGLsForThisEvent → NO_CALL genotype with PLs
//! USE_PLS_TO_ASSIGN         → GT = argmin(PL)
//! AlleleFrequencyCalculator.calculate
//! calculateOutputAlleleSubset:
//!   keep ALT a iff AFresult.passesThreshold(a, stand-call-conf=30)
//!   siteIsMonomorphic iff no ALT is plausible
//! passesEmitThreshold(QUAL, siteIsMonomorphic):
//!   (EMIT_ALL_CONFIDENT_SITES || !siteIsMonomorphic) && QUAL >= 30
//! if !passesEmitThreshold → calculateGenotypes returns null (no VCF record)
//! ```
//!
//! Live HOLDOUT_6R53 (`20:29455388 C/T`):
//! Java GLs/PLs are hom-ref `PL=0,6,1780` → GT 0/0, T fails `passesThreshold(30)`,
//! `calculateGenotypes` is null.
//! Rust GLs/PLs are het `PL=81,0,36` → GT 0/1, QUAL 73.64, emit true at conf 10 and 30.
//!
//! The emit / subset **predicates agree** given the same GLs. The first unequal
//! object is the GL/PL vector itself. Production change: NONE.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r105_genotype_emit_contract
//! HOLDOUT_6R105=1 cargo test -p gatk-haplotypecaller --test holdout_6r105_genotype_emit -- --nocapture
//! ```

use gatk_haplotypecaller::emit_gates::{
    java_emit_would_pass, passes_hc_variant_emit_biallelic, passes_java_emit_not_hom_ref,
};
use gatk_haplotypecaller::event_map::VariationEvent;
use gatk_haplotypecaller::genotyping::{
    best_pl_index, diploid_genotype_alleles_from_pl_index, emit_genotype_format_fields,
};
use gatk_haplotypecaller::hc_genotyping_engine::{
    java_emit_af_decision, DEFAULT_STAND_EMIT_CONFIDENCE,
};

/// Java 4.4 `GenotypeCalculationArgumentCollection.DEFAULT_STANDARD_CONFIDENCE_FOR_CALLING`.
const JAVA_STAND_CALL_CONF: f64 = 30.0;

fn snp(r: &str, a: &str) -> VariationEvent {
    VariationEvent::from_alleles("chr", 100, r, a)
}

fn gl_from_pl(pl: &[i32]) -> Vec<f64> {
    pl.iter().map(|&p| (p as f64) / -10.0).collect()
}

fn keep_from_gt(n_alleles: usize, gt: &[i32]) -> Vec<usize> {
    let mut keep = vec![0];
    for i in 1..n_alleles {
        if gt.iter().any(|&g| g >= 0 && (g as usize) == i) {
            keep.push(i);
        }
    }
    keep
}

/// Live Java `calculateGLsForThisEvent` PLs at the holdout (hom-ref). Coordinate-free vector.
fn java_holdout_homref_pl() -> [i32; 3] {
    [0, 6, 1780]
}

/// Live Rust genotyped PLs at the holdout (het). Coordinate-free vector.
fn rust_holdout_het_pl() -> [i32; 3] {
    [81, 0, 36]
}

#[test]
fn forensic_6r105_use_pls_to_assign_is_argmin_pl() {
    let event = snp("C", "T");
    let _ = event;
    let java_pl = java_holdout_homref_pl();
    let rust_pl = rust_holdout_het_pl();
    let java_gl = gl_from_pl(&java_pl);
    let rust_gl = gl_from_pl(&rust_pl);
    let java_fmt = emit_genotype_format_fields(&java_gl, &[40, 0]).expect("fmt java");
    let rust_fmt = emit_genotype_format_fields(&rust_gl, &[44, 4]).expect("fmt rust");
    assert_eq!(
        diploid_genotype_alleles_from_pl_index(2, best_pl_index(&java_fmt.pl)),
        [0, 0]
    );
    assert_eq!(
        diploid_genotype_alleles_from_pl_index(2, best_pl_index(&rust_fmt.pl)),
        [0, 1]
    );
}

#[test]
fn forensic_6r105_homref_gls_fail_java_emit_at_10_and_30() {
    let event = snp("C", "T");
    let gl = gl_from_pl(&java_holdout_homref_pl());
    let fmt = emit_genotype_format_fields(&gl, &[40, 0]).expect("fmt");
    assert!(
        !passes_java_emit_not_hom_ref(&gl, &fmt),
        "PL 0,6,1780 is hom-ref after USE_PLS_TO_ASSIGN"
    );
    let af30 = java_emit_af_decision(&gl, JAVA_STAND_CALL_CONF).expect("af30");
    assert!(af30.site_is_monomorphic);
    assert!(!af30.alt_plausible);
    assert!(!af30.passes_emit);
    assert!(!java_emit_would_pass(&event, &gl, &fmt, JAVA_STAND_CALL_CONF, &[]).unwrap());
    assert!(!java_emit_would_pass(&event, &gl, &fmt, DEFAULT_STAND_EMIT_CONFIDENCE, &[]).unwrap());
}

#[test]
fn forensic_6r105_het_gls_pass_java_emit_at_30() {
    let event = snp("C", "T");
    let gl = gl_from_pl(&rust_holdout_het_pl());
    let fmt = emit_genotype_format_fields(&gl, &[44, 4]).expect("fmt");
    assert!(passes_java_emit_not_hom_ref(&gl, &fmt));
    assert!(passes_hc_variant_emit_biallelic(&gl, JAVA_STAND_CALL_CONF).unwrap());
    let af30 = java_emit_af_decision(&gl, JAVA_STAND_CALL_CONF).expect("af30");
    assert!(!af30.site_is_monomorphic);
    assert!(af30.alt_plausible);
    assert!(af30.phred_scaled >= JAVA_STAND_CALL_CONF);
    assert!(java_emit_would_pass(&event, &gl, &fmt, JAVA_STAND_CALL_CONF, &[]).unwrap());
    assert!(java_emit_would_pass(&event, &gl, &fmt, DEFAULT_STAND_EMIT_CONFIDENCE, &[]).unwrap());
}

#[test]
fn forensic_6r105_output_allele_subset_follows_assigned_gt() {
    assert_eq!(
        keep_from_gt(2, &[0, 0]),
        vec![0],
        "hom-ref drops the SNP ALT"
    );
    assert_eq!(
        keep_from_gt(2, &[0, 1]),
        vec![0, 1],
        "het keeps the SNP ALT"
    );
}

#[test]
fn forensic_6r105_emit_threshold_not_the_first_divergence() {
    let event = snp("C", "T");
    let java_gl = gl_from_pl(&java_holdout_homref_pl());
    let rust_gl = gl_from_pl(&rust_holdout_het_pl());
    let java_fmt = emit_genotype_format_fields(&java_gl, &[40, 0]).expect("j");
    let rust_fmt = emit_genotype_format_fields(&rust_gl, &[44, 4]).expect("r");
    let rust_on_java =
        java_emit_would_pass(&event, &java_gl, &java_fmt, JAVA_STAND_CALL_CONF, &[]).unwrap();
    let java_on_rust =
        java_emit_would_pass(&event, &rust_gl, &rust_fmt, JAVA_STAND_CALL_CONF, &[]).unwrap();
    assert!(
        !rust_on_java,
        "Rust emit predicate on Java GLs must not emit (same as calculateGenotypes=null)"
    );
    assert!(
        java_on_rust,
        "Java emit predicate on Rust GLs would emit (QUAL>30, not hom-ref)"
    );
}
