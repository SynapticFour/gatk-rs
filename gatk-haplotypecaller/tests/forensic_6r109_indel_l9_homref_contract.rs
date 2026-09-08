//! 6R.109 coordinate-free: valid hom-ref calculator GLs are not replaced by indel L9.
//!
//! Java 4.4 preserves a valid `calculateGLsForThisEvent` result after emit failure.
//! Post-`finalize_site` indel L9 must not replace that result with [`SparsePlShape::Het`].
//! Empty / missing calculator GLs keep the existing genome-wide indel fallback.
//! SNP HomAltStrong (6R.107) is unchanged.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r109_indel_l9_homref_contract
//! HOLDOUT_6R109=1 cargo test -p gatk-haplotypecaller --test holdout_6r109_indel_l9 -- --nocapture
//! ```

use gatk_haplotypecaller::event_map::VariationEvent;
use gatk_haplotypecaller::genotyping::{
    best_pl_index, diploid_genotype_alleles_from_pl_index, emit_genotype_format_fields,
};
use gatk_haplotypecaller::hc_genotyping_engine::{
    java_emit_would_pass, l9_may_overwrite_pairhmm_gls_after_emit_fail, SparsePlShape,
    DEFAULT_STAND_EMIT_CONFIDENCE,
};

const JAVA_STAND_CALL_CONF: f64 = 30.0;

fn snp() -> VariationEvent {
    VariationEvent::from_alleles("20", 100, "C", "T")
}

fn indel() -> VariationEvent {
    VariationEvent::from_alleles("20", 100, "A", "AT")
}

fn gl_from_pl(pl: &[i32]) -> Vec<f64> {
    pl.iter().map(|&p| (p as f64) / -10.0).collect()
}

fn calculator_is_hom_ref(pl: &[gatk_haplotypecaller::bio_ids::PhredLikelihood]) -> bool {
    !pl.is_empty() && best_pl_index(pl) == 0
}

#[test]
fn forensic_6r109_valid_homref_calculator_is_not_overwritten_for_indel() {
    let event = indel();
    let calc = gl_from_pl(&[0, 43, 4393]);
    let fmt = emit_genotype_format_fields(&calc, &[8, 4]).expect("fmt");
    assert!(calculator_is_hom_ref(&fmt.pl));
    assert_eq!(
        diploid_genotype_alleles_from_pl_index(2, best_pl_index(&fmt.pl)).as_slice(),
        [0, 0]
    );
    assert!(!java_emit_would_pass(&event, &calc, &fmt, JAVA_STAND_CALL_CONF, &[]).unwrap());
    assert!(
        !java_emit_would_pass(&event, &calc, &fmt, DEFAULT_STAND_EMIT_CONFIDENCE, &[]).unwrap()
    );
    assert_eq!(SparsePlShape::from_pileup_depths(8, 4), SparsePlShape::Het);
    assert!(
        !l9_may_overwrite_pairhmm_gls_after_emit_fail(&event, 8, 4, true),
        "valid hom-ref calculator GLs must not be replaced by indel L9"
    );
}

#[test]
fn forensic_6r109_missing_calculator_gls_keep_indel_fallback() {
    let event = indel();
    assert!(!calculator_is_hom_ref(&[]));
    assert!(
        l9_may_overwrite_pairhmm_gls_after_emit_fail(&event, 0, 4, false),
        "empty-mapper / missing calculator GLs keep the existing indel fallback"
    );
    assert!(!l9_may_overwrite_pairhmm_gls_after_emit_fail(
        &event, 43, 2, false
    ));
}

#[test]
fn forensic_6r109_snp_hom_alt_strong_contract_unchanged() {
    let event = snp();
    assert!(
        !l9_may_overwrite_pairhmm_gls_after_emit_fail(&event, 44, 4, true),
        "SNP non-HomAltStrong still does not overwrite"
    );
    assert!(
        l9_may_overwrite_pairhmm_gls_after_emit_fail(&event, 0, 4, true),
        "SNP HomAltStrong overwrite remains permitted"
    );
    assert!(!l9_may_overwrite_pairhmm_gls_after_emit_fail(
        &VariationEvent::from_alleles("2", 100, "C", "T"),
        0,
        4,
        true
    ));
}
