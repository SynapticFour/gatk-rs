//! 6R.107 coordinate-free: L9 post-emit-fail SparsePlShape overwrite gate.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`) does not replace a
//! valid `calculateGLsForThisEvent` result with a pileup shape because `alt_ad >= 1`.
//! Hom-ref calculator GLs fail `passesEmitThreshold` → `calculateGenotypes` is null.
//!
//! Rust L9 comment (post-`finalize_site` None): overwrite when PairHMM looks
//! non-variant **and** BAM pileup is hom-alt / strong-alt. That class is
//! [`SparsePlShape::from_pileup_depths`] → [`SparsePlShape::HomAltStrong`], already
//! defined as `alt >= 2 && (ref == 0 || alt >= 4*ref)`.
//!
//! `genome_wide_genotype_read_support` for SNPs remains `alt_ad >= 1` (empty-mapper /
//! empty-subset fallback). Post-emit-fail overwrite of existing PairHMM GLs uses
//! [`l9_may_overwrite_pairhmm_gls_after_emit_fail`].
//!
//! Holdout phenotype `20:29455388 C/T` is REF 44 / ALT 4 → `Het`, not HomAltStrong.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r107_l9_overwrite_contract
//! HOLDOUT_6R107=1 cargo test -p gatk-haplotypecaller --test holdout_6r107_l9_overwrite -- --nocapture
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

fn snp_on(contig: &str, r: &str, a: &str) -> VariationEvent {
    VariationEvent::from_alleles(contig, 100, r, a)
}

fn indel() -> VariationEvent {
    VariationEvent::from_alleles("20", 100, "AT", "A")
}

fn gl_from_pl(pl: &[i32]) -> Vec<f64> {
    pl.iter().map(|&p| (p as f64) / -10.0).collect()
}

#[test]
fn forensic_6r107_calculator_homref_fails_java_emit() {
    let event = snp_on("20", "C", "T");
    let java_gl = gl_from_pl(&[0, 6, 1780]);
    let java_fmt = emit_genotype_format_fields(&java_gl, &[43, 4]).expect("fmt");
    assert_eq!(
        diploid_genotype_alleles_from_pl_index(2, best_pl_index(&java_fmt.pl)).as_slice(),
        [0, 0]
    );
    assert!(!java_emit_would_pass(&event, &java_gl, &java_fmt, JAVA_STAND_CALL_CONF, &[]).unwrap());
    assert!(!java_emit_would_pass(
        &event,
        &java_gl,
        &java_fmt,
        DEFAULT_STAND_EMIT_CONFIDENCE,
        &[]
    )
    .unwrap());
}

#[test]
fn forensic_6r107_ref44_alt4_is_het_not_hom_alt_strong() {
    assert_eq!(SparsePlShape::from_pileup_depths(44, 4), SparsePlShape::Het);
    assert_eq!(SparsePlShape::Het.pl(), [81, 0, 36]);
    assert!(!SparsePlShape::pileup_is_hom_alt_strong(44, 4));
}

#[test]
fn forensic_6r107_genome_wide_alt_ge_1_is_not_post_emit_overwrite() {
    let event = snp_on("20", "C", "T");
    assert!(
        !l9_may_overwrite_pairhmm_gls_after_emit_fail(&event, 44, 4),
        "valid PairHMM GLs + REF-majority pileup must not be replaced"
    );
}

#[test]
fn forensic_6r107_hom_alt_strong_pileup_may_overwrite() {
    let event = snp_on("20", "C", "T");
    assert!(SparsePlShape::pileup_is_hom_alt_strong(0, 4));
    assert!(l9_may_overwrite_pairhmm_gls_after_emit_fail(&event, 0, 4));
    assert!(l9_may_overwrite_pairhmm_gls_after_emit_fail(&event, 1, 8));
}

#[test]
fn forensic_6r107_p12_scope_never_takes_this_l9() {
    let event = snp_on("2", "C", "T");
    assert!(!l9_may_overwrite_pairhmm_gls_after_emit_fail(&event, 0, 4));
    assert!(!l9_may_overwrite_pairhmm_gls_after_emit_fail(&event, 44, 4));
}

#[test]
fn forensic_6r107_no_alt_does_not_overwrite() {
    let event = snp_on("20", "C", "T");
    assert!(!l9_may_overwrite_pairhmm_gls_after_emit_fail(&event, 44, 0));
}

#[test]
fn forensic_6r107_indel_keeps_genome_wide_gate() {
    let event = indel();
    assert!(l9_may_overwrite_pairhmm_gls_after_emit_fail(&event, 0, 4));
    assert!(!l9_may_overwrite_pairhmm_gls_after_emit_fail(&event, 43, 2));
}

#[test]
fn forensic_6r107_preserving_calculator_does_not_emit_het_shape() {
    let event = snp_on("20", "C", "T");
    let calc = gl_from_pl(&[0, 6, 1780]);
    let calc_fmt = emit_genotype_format_fields(&calc, &[44, 4]).expect("c");
    let het = SparsePlShape::Het.gl_vec();
    let het_fmt = emit_genotype_format_fields(&het, &[44, 4]).expect("h");
    assert!(!l9_may_overwrite_pairhmm_gls_after_emit_fail(&event, 44, 4));
    assert!(!java_emit_would_pass(&event, &calc, &calc_fmt, JAVA_STAND_CALL_CONF, &[]).unwrap());
    assert!(java_emit_would_pass(&event, &het, &het_fmt, JAVA_STAND_CALL_CONF, &[]).unwrap());
}
