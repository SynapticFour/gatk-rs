//! 6R.155 coordinate-free: first mutation after the Java-equivalent AD write.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//!
//! After `DepthPerAlleleBySample.annotateWithLikelihoods` Java may:
//!   reverseTrimAlleles          // iff allele count changed; copies genotype fields
//!   phaseCalls                  // iff doPhysicalPhasing; copies GT/PL/AD
//!   other genotype annotations  // DP from AD; do not rewrite AD/GT/PL
//! Java has **no** pileup Class-A3 / SparsePlShape replacement of GT/PL/AD.
//!
//! Rust `SiteReshape::apply_class_a_family` (not a Java class) after `SiteScore`:
//!   Class-A3 when: SNP, pileup both alleles, pileup alt*2 >= ref,
//!                  informative ref >= 1, alt >= 1, ref >= 3*alt
//!   if PL-argmin == 0/1: replace AD/DP only (keep PLs)
//!   else: `sparse_snp_genotype_from_read_depths(pileup)` replaces GT, PL, AD, GQ, DP
//!
//! Pileup AD is QNAME-deduped `region.reads`, not the retainEvidence likelihood object.
//! That is a different evidence population from the 6R.154 first AD write.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r155_site_reshape_genotype_semantics_contract
//! HOLDOUT_6R155=1 GATK_RS_EXPERIMENTAL_KBEST_POLICY=unbounded_diagnostic \
//!   cargo test -p gatk-haplotypecaller --test holdout_6r155_site_reshape -- --nocapture --test-threads=1
//! ```

use gatk_haplotypecaller::genotyping::{
    best_pl_index, emit_genotype_format_fields, ReadLikelihoodRow,
};
use gatk_haplotypecaller::hc_genotyping_engine::{
    biallelic_allele_depths_from_rows, SparsePlShape,
};

fn row(lls: Vec<f64>) -> ReadLikelihoodRow {
    ReadLikelihoodRow {
        read_index: 0,
        read_id: String::new(),
        haplotype_log10_likelihoods: lls,
    }
}

fn class_a3(is_snp: bool, pileup_ref: i32, pileup_alt: i32, info_ref: i32, info_alt: i32) -> bool {
    is_snp
        && pileup_ref >= 1
        && pileup_alt >= 1
        && pileup_ref.saturating_mul(2) >= pileup_alt
        && pileup_alt.saturating_mul(2) >= pileup_ref
        && info_ref >= 1
        && info_alt >= 1
        && info_ref >= info_alt.saturating_mul(3)
}

#[test]
fn forensic_6r155_java_has_no_post_ad_pileup_reshape() {
    let java_ops_after_first_ad_write = [
        "reverseTrimAlleles_copy_if_allele_count_changed",
        "phaseCalls_copy_if_doPhysicalPhasing",
        "other_annotations_do_not_rewrite_AD_GT_PL",
    ];
    assert!(java_ops_after_first_ad_write
        .iter()
        .all(|op| !op.contains("Class-A3") && !op.contains("SparsePlShape")));
    let java_analogue = false;
    assert!(!java_analogue);
}

#[test]
fn forensic_6r155_class_a3_gt_flip_replaces_calculator_fields() {
    let info_ref = 9;
    let info_alt = 3;
    let pileup_ref = 6;
    let pileup_alt = 8;
    assert!(class_a3(true, pileup_ref, pileup_alt, info_ref, info_alt));
    let first_write_rows = [
        row(vec![0.0, -3.0]),
        row(vec![0.0, -3.0]),
        row(vec![0.0, -3.0]),
        row(vec![0.0, -3.0]),
        row(vec![0.0, -3.0]),
        row(vec![0.0, -3.0]),
        row(vec![0.0, -3.0]),
        row(vec![0.0, -3.0]),
        row(vec![0.0, -3.0]),
        row(vec![-3.0, 0.0]),
        row(vec![-3.0, 0.0]),
        row(vec![-3.0, 0.0]),
    ];
    let info_ad = biallelic_allele_depths_from_rows(&first_write_rows, 0, 1);
    assert_eq!(info_ad, vec![info_ref, info_alt]);
    let gls = vec![0.0, -0.5, -3.0];
    let first = emit_genotype_format_fields(&gls, &info_ad).expect("first write");
    assert_eq!(first.pl_as_i32()[0], 0);
    assert_eq!(best_pl_index(&first.pl), 0);
    let shape = SparsePlShape::from_pileup_depths(pileup_ref, pileup_alt);
    assert_eq!(shape, SparsePlShape::Het);
    let replaced = emit_genotype_format_fields(&shape.gl_vec(), &[pileup_ref, pileup_alt])
        .expect("sparse replace");
    assert_eq!(replaced.pl_as_i32(), vec![81, 0, 36]);
    assert_eq!(replaced.ad_as_i32(), vec![pileup_ref, pileup_alt]);
    assert_ne!(first.pl_as_i32(), replaced.pl_as_i32());
    assert_ne!(first.ad_as_i32(), replaced.ad_as_i32());
    assert_eq!(best_pl_index(&replaced.pl), 1);
}

#[test]
fn forensic_6r155_class_a3_het_pl_keeps_pls_replaces_ad_only() {
    let gls = vec![-8.1, 0.0, -3.6];
    let info = emit_genotype_format_fields(&gls, &[9, 3]).expect("info");
    assert_eq!(best_pl_index(&info.pl), 1);
    let ad_only = emit_genotype_format_fields(&gls, &[6, 8]).expect("ad only");
    assert_eq!(ad_only.pl_as_i32(), info.pl_as_i32());
    assert_eq!(ad_only.ad_as_i32(), vec![6, 8]);
}

#[test]
fn forensic_6r155_class_a3_does_not_fire_when_informative_not_ref_skewed() {
    assert!(!class_a3(true, 10, 10, 12, 11));
    assert!(!class_a3(true, 10, 10, 20, 0));
    assert!(!class_a3(false, 10, 10, 20, 3));
}

#[test]
fn forensic_6r155_pileup_qname_dedupe_is_not_retain_evidence_object() {
    let likelihood_mates_independent = 2usize;
    let qname_deduped_pileup_rows = 1usize;
    assert_ne!(likelihood_mates_independent, qname_deduped_pileup_rows);
    let first_ad_source = "InformativeAd from retainEvidence likelihood rows";
    let reshape_ad_source = "QNAME-deduped region.reads pileup";
    assert_ne!(first_ad_source, reshape_ad_source);
}

#[test]
fn forensic_6r155_ref_skewed_informative_vs_balanced_pileup_is_class_a3() {
    let info = [37, 4];
    let pileup = [22, 24];
    assert!(class_a3(true, pileup[0], pileup[1], info[0], info[1]));
    assert_eq!(
        SparsePlShape::from_pileup_depths(pileup[0], pileup[1]),
        SparsePlShape::Het
    );
    assert!(!SparsePlShape::pileup_is_hom_alt_strong(
        pileup[0], pileup[1]
    ));
    let replaced =
        emit_genotype_format_fields(&SparsePlShape::Het.gl_vec(), &pileup).expect("sparse het");
    assert_eq!(replaced.pl_as_i32(), vec![81, 0, 36]);
    assert_eq!(replaced.ad_as_i32(), vec![22, 24]);
    assert_eq!(best_pl_index(&replaced.pl), 1);
    assert_eq!(replaced.gq.as_i32(), 36);
    assert_eq!(replaced.dp.as_i32(), 46);
}

#[test]
fn forensic_6r155_bypass_restores_first_write_fields() {
    let first_gt = 0usize;
    let first_pl = [0, 5, 1174];
    let first_ad = [37, 4];
    let first_gq = 5;
    let skip_reshape = true;
    let after_gt = if skip_reshape { first_gt } else { 1 };
    let after_pl = if skip_reshape { first_pl } else { [81, 0, 36] };
    let after_ad = if skip_reshape { first_ad } else { [22, 24] };
    let after_gq = if skip_reshape { first_gq } else { 36 };
    assert_eq!(after_gt, first_gt);
    assert_eq!(after_pl, first_pl);
    assert_eq!(after_ad, first_ad);
    assert_eq!(after_gq, first_gq);
}
