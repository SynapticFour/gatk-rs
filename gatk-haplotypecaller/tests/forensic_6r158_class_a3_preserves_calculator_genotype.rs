//! 6R.158 coordinate-free: Class-A3 preserves an assigned calculator genotype.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//!
//! Production change (this round): `SiteReshape::apply_class_a_family` returns the
//! incoming calculator genotype when Class-A3 predicates hold. The generic helper
//! `sparse_snp_genotype_from_read_depths` is unchanged.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r158_class_a3_preserves_calculator_genotype
//! cargo test -p gatk-haplotypecaller --lib forensic_6r158 -- --test-threads=1
//! HOLDOUT_6R158=1 GATK_RS_EXPERIMENTAL_KBEST_POLICY=unbounded_diagnostic \
//!   cargo test -p gatk-haplotypecaller --test holdout_6r158_class_a3_preserve -- --nocapture --test-threads=1
//! ```

use gatk_haplotypecaller::event_map::VariationEvent;
use gatk_haplotypecaller::genotyping::{
    best_pl_index, emit_genotype_format_fields, HaplotypeLikelihoodAggregation,
};
use gatk_haplotypecaller::hc_genotyping_engine::{
    HcGenotypingConfig, RegionGenotypeResult, SiteReshape, SparsePlShape,
};
use gatk_haplotypecaller::GenomePosition;

fn snp_event() -> VariationEvent {
    VariationEvent {
        contig: "synth".into(),
        start_1based: GenomePosition::new_1based(1),
        end_1based: GenomePosition::new_1based(1),
        ref_allele: "A".into(),
        alt_allele: "T".into(),
    }
}

fn assigned(gls: &[f64], info_ad: [i32; 2]) -> RegionGenotypeResult {
    let format = emit_genotype_format_fields(gls, &info_ad).expect("assigned FORMAT");
    RegionGenotypeResult {
        aggregation: HaplotypeLikelihoodAggregation {
            haplotype_log10_sums: vec![0.0, -1.0],
            read_count: 2,
        },
        best_haplotype_index: 0,
        ref_haplotype_index: 0,
        alt_haplotype_index: 1,
        genotype_log10_likelihoods: gls.to_vec(),
        format,
    }
}

fn fmt_tuple(gt: &RegionGenotypeResult) -> (usize, Vec<i32>, Vec<i32>, i32, i32) {
    (
        best_pl_index(&gt.format.pl),
        gt.format.pl_as_i32(),
        gt.format.ad_as_i32(),
        gt.format.gq.as_i32(),
        gt.format.dp.as_i32(),
    )
}

const A3_PILEUP: [i32; 2] = [22, 24];

#[test]
fn forensic_6r158_case_a_valid_00_preserved() {
    let incoming = assigned(&[0.0, -0.5, -117.4], [37, 4]);
    let before = fmt_tuple(&incoming);
    assert_eq!(before, (0, vec![0, 5, 1174], vec![37, 4], 5, 41));
    let out = SiteReshape::apply_class_a_family(
        incoming.clone(),
        &snp_event(),
        A3_PILEUP[0],
        A3_PILEUP[1],
        &HcGenotypingConfig::strict_java(),
    )
    .expect("A3 skip");
    assert_eq!(fmt_tuple(&out), before);
    assert_eq!(out.format, incoming.format);
    assert_eq!(
        out.genotype_log10_likelihoods,
        incoming.genotype_log10_likelihoods
    );
    assert_eq!(out.best_haplotype_index, incoming.best_haplotype_index);
}

#[test]
fn forensic_6r158_case_b_valid_01_preserved() {
    let incoming = assigned(&[-8.1, 0.0, -3.6], [37, 4]);
    let before = fmt_tuple(&incoming);
    assert_eq!(before.0, 1);
    let out = SiteReshape::apply_class_a_family(
        incoming.clone(),
        &snp_event(),
        A3_PILEUP[0],
        A3_PILEUP[1],
        &HcGenotypingConfig::strict_java(),
    )
    .expect("A3 skip");
    assert_eq!(fmt_tuple(&out), before);
    assert_eq!(out.format, incoming.format);
}

#[test]
fn forensic_6r158_case_c_valid_11_preserved() {
    let incoming = assigned(&[-9.0, -0.6, 0.0], [37, 4]);
    let before = fmt_tuple(&incoming);
    assert_eq!(before.0, 2);
    assert_eq!(before.1, vec![90, 6, 0]);
    let out = SiteReshape::apply_class_a_family(
        incoming.clone(),
        &snp_event(),
        A3_PILEUP[0],
        A3_PILEUP[1],
        &HcGenotypingConfig::strict_java(),
    )
    .expect("A3 skip");
    assert_eq!(fmt_tuple(&out), before);
    assert_eq!(out.format, incoming.format);
}

#[test]
fn forensic_6r158_case_d_a3_false_non_a3_unchanged() {
    let no_class = assigned(&[0.0, -0.5, -117.4], [10, 8]);
    let before = fmt_tuple(&no_class);
    let out = SiteReshape::apply_class_a_family(
        no_class.clone(),
        &snp_event(),
        A3_PILEUP[0],
        A3_PILEUP[1],
        &HcGenotypingConfig::strict_java(),
    )
    .expect("no-op");
    assert_eq!(fmt_tuple(&out), before);

    let a2_only = assigned(&[-9.0, -0.6, 0.0], [4, 12]);
    assert_eq!(fmt_tuple(&a2_only).0, 2);
    let a2_out = SiteReshape::apply_class_a_family(
        a2_only,
        &snp_event(),
        A3_PILEUP[0],
        A3_PILEUP[1],
        &HcGenotypingConfig::strict_java(),
    )
    .expect("A2 still runs");
    assert_eq!(
        fmt_tuple(&a2_out),
        (1, vec![81, 0, 36], vec![22, 24], 36, 46)
    );
}

#[test]
fn forensic_6r158_case_e_sparse_helper_shape_unchanged() {
    let shape = SparsePlShape::from_pileup_depths(22, 24);
    assert_eq!(shape, SparsePlShape::Het);
    let fmt = emit_genotype_format_fields(&shape.gl_vec(), &[22, 24]).expect("helper shape");
    assert_eq!(best_pl_index(&fmt.pl), 1);
    assert_eq!(fmt.pl_as_i32(), vec![81, 0, 36]);
    assert_eq!(fmt.ad_as_i32(), vec![22, 24]);
    assert_eq!(fmt.gq.as_i32(), 36);
    assert_eq!(fmt.dp.as_i32(), 46);
}
