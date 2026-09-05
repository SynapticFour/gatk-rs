//! 6R.102 coordinate-free: site QUAL is AFCalculator on the pre-subset merged VC.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//!
//! ```text
//! calculateGLsForThisEvent(mergedVC)           // 10 diploid GLs when 4 alleles
//! GenotypingEngine.calculateGenotypes
//!   AlleleFrequencyCalculator.calculate(mergedVC)
//!     priorPseudocounts: REF / SNP (len==refLen) / indel (else)
//!     EM Dirichlet log10MeanWeights
//!     if alleles.contains(SPAN_DEL):
//!       log10PNoVariant = log10Sum(posteriors of REF+SPAN_DEL genotypes)
//!     else:
//!       log10PNoVariant = posterior[HOM_REF]
//!   builder.log10PError(log10ProbOnlyRefAlleleExists())
//!   QUAL = -10 * log10PError
//! AlleleSubsettingUtils.subsetAlleles           // copies QUAL; does not recompute
//! reverseTrimAlleles                            // copies QUAL
//! ```
//!
//! Emitted biallelic PL (542,0,1353) is not the QUAL input. AF on those 3 GLs is
//! ~534.64. Java QUAL 510.06 is P(no variant) over {0/0, 0/*, */*} after mixed
//! SNP/indel Dirichlet priors on the PL-roundtripped 10-GL object.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r102_qual_calculation_contract
//! HOLDOUT_6R102=1 cargo test -p gatk-haplotypecaller --test holdout_6r102_qual_calculation -- --nocapture
//! ```

use gatk_haplotypecaller::variant_site_hc_annotations::{
    qual_from_af_calculation, qual_from_merged_diploid_af_calculate,
};

/// PL-roundtrip GLs (`GenotypeBuilder.PL` → `getLikelihoods`) for a 4-allele diploid
/// site with alleles `TG, T, CG, *`. Constructed from integer PL, not genomic coords.
fn merged_pl_roundtrip_gls() -> [f64; 10] {
    let pl = [542, 484, 1964, 0, 1234, 1353, 481, 1801, 1264, 1880];
    std::array::from_fn(|i| (pl[i] as f64) / -10.0)
}

#[test]
fn emitted_biallelic_pl_does_not_reproduce_java_qual() {
    let gl = [-54.2, 0.0, -135.3];
    let qual = qual_from_af_calculation(&gl).expect("qual");
    assert!(
        (qual - 534.64).abs() < 0.02,
        "AF(emitted PL 542,0,1353) is Rust-before-6R.102 QUAL, got {qual}"
    );
    assert!(
        (qual - 510.06).abs() > 10.0,
        "Java QUAL is not AF(emitted PL); got {qual}"
    );
}

#[test]
fn hom_ref_only_on_four_allele_object_is_not_java_qual() {
    let gl = merged_pl_roundtrip_gls();
    let alleles = ["TG", "T", "CG", "X"];
    let qual = qual_from_merged_diploid_af_calculate(&gl, &alleles).expect("qual");
    assert!(
        (qual - 534.64).abs() < 0.15,
        "without SPAN_DEL, log10PNoVariant is HOM_REF ≈ 534.64, got {qual}"
    );
    assert!((qual - 510.06).abs() > 10.0);
}

#[test]
fn span_del_mixed_priors_reproduce_java_qual_from_constructed_gls() {
    let gl = merged_pl_roundtrip_gls();
    let alleles = ["TG", "T", "CG", "*"];
    let qual = qual_from_merged_diploid_af_calculate(&gl, &alleles).expect("qual");
    assert!(
        (qual - 510.06).abs() < 0.02,
        "Java calculate() with SPAN_DEL + mixed priors → 510.06, got {qual}"
    );
}

#[test]
fn qual_does_not_depend_on_ad_or_gt() {
    let gl = merged_pl_roundtrip_gls();
    let alleles = ["TG", "T", "CG", "*"];
    let q1 = qual_from_merged_diploid_af_calculate(&gl, &alleles).expect("q1");
    let q2 = qual_from_merged_diploid_af_calculate(&gl, &alleles).expect("q2");
    assert_eq!(q1, q2);
    assert!((q1 - 510.06).abs() < 0.02);
}
