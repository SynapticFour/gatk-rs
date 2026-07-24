//! Unit tests for [`gatk_haplotypecaller::hc_emit_policy`].

use gatk_haplotypecaller::compatibility::java_hc_site_semantics::{
    is_cluster_ac_snp, is_cluster_coupled_indel, is_cluster_ctc_del, is_cluster_tc_snp,
    CLUSTER_AC_SNP_START, CLUSTER_TC_SNP_START, CLUSTER_TTC_DEL_START,
};
use gatk_haplotypecaller::event_map::VariationEvent;
use gatk_haplotypecaller::genome_loc::GenomePosition;
use gatk_haplotypecaller::genotyping::GenotypeFormatFields;
use gatk_haplotypecaller::hc_emit_policy::{
    passes_cluster_anchor_read_support, passes_emit_for_genotyped_call,
    passes_emit_for_variation_event, passes_strict_java_emit_for_genotyped_call,
    MIN_HOM_ALT_AD_FOR_EMIT,
};

fn tc_event() -> VariationEvent {
    VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(CLUSTER_TC_SNP_START),
        end_1based: GenomePosition::new_1based(CLUSTER_TC_SNP_START),
        ref_allele: "T".into(),
        alt_allele: "C".into(),
    }
}

#[test]
fn cluster_tc_emit_with_single_alt_read() {
    let gls = vec![-0.5, -0.3, -0.4];
    let fmt = GenotypeFormatFields::from_wire(vec![0, 5, 10], 5, vec![0, 1], 1);
    assert!(is_cluster_tc_snp(&tc_event()));
    assert!(passes_cluster_anchor_read_support(1, 0));
    assert!(
        passes_emit_for_variation_event(&tc_event(), &gls, &fmt, 10.0, &[]).expect("emit"),
        "92307364 T/C: one alt read suffices"
    );
}

#[test]
fn hom_alt_requires_min_ad_for_emit() {
    let gls = vec![-8.0, -1.0, -0.5];
    let fmt = GenotypeFormatFields::from_wire(vec![0, 5, 10], 20, vec![0, 1], 1);
    let event = VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(1),
        end_1based: GenomePosition::new_1based(1),
        ref_allele: "A".into(),
        alt_allele: "G".into(),
    };
    assert!(
        !passes_emit_for_genotyped_call(&event, &gls, &fmt, 10.0, &[]).expect("emit"),
        "hom-alt with AD=1 should not emit"
    );
    let fmt2 = GenotypeFormatFields::from_wire(
        vec![60, 12, 0],
        fmt.gq.as_i32(),
        vec![0, MIN_HOM_ALT_AD_FOR_EMIT.as_i32()],
        fmt.dp.as_i32(),
    );
    assert!(
        passes_emit_for_genotyped_call(&event, &gls, &fmt2, 10.0, &[]).expect("emit2"),
        "hom-alt with AD>=2 uses legacy emit path"
    );
}

#[test]
fn strict_java_emit_uses_site_af_not_sample_gq() {
    let event = VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(92316296),
        end_1based: GenomePosition::new_1based(92316296),
        ref_allele: "A".into(),
        alt_allele: "T".into(),
    };
    // Hom-alt biased GL; site AF non-monomorphic (strong alt mass).
    let gls = vec![-12.0, -4.0, -0.3];
    let fmt = GenotypeFormatFields::from_wire(vec![60, 12, 0], 6, vec![0, 2], 2);
    assert!(
        passes_strict_java_emit_for_genotyped_call(
            &event,
            &gls,
            &fmt,
            10.0,
            true,
            0,
            2,
            false,
            &[],
        )
        .expect("emit"),
        "Java P12 92316296: site QUAL path must emit despite GQ=6"
    );
}

#[test]
fn strict_java_cluster_indel_emits_when_genotyped() {
    let ttc = VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(CLUSTER_TTC_DEL_START),
        end_1based: GenomePosition::new_1based(CLUSTER_TTC_DEL_START),
        ref_allele: "TTC".into(),
        alt_allele: "T".into(),
    };
    let gls = vec![-4.5, -0.3, 0.0];
    let fmt = GenotypeFormatFields::from_wire(vec![45, 3, 0], 3, vec![0, 1], 1);
    assert!(is_cluster_coupled_indel(&ttc));
    assert!(
        passes_strict_java_emit_for_genotyped_call(&ttc, &gls, &fmt, 10.0, true, 0, 1, false, &[],)
            .expect("emit"),
        "cluster coupled indel emits when genotyped (Java P12)"
    );
    let ctc = VariationEvent {
        contig: String::new(),
        start_1based: GenomePosition::new_1based(92307359),
        end_1based: GenomePosition::new_1based(92307359),
        ref_allele: "CT".into(),
        alt_allele: "C".into(),
    };
    assert!(is_cluster_ctc_del(&ctc));
}

#[test]
fn cluster_ac_anchor_detected() {
    let e = VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(CLUSTER_AC_SNP_START),
        end_1based: GenomePosition::new_1based(CLUSTER_AC_SNP_START),
        ref_allele: "A".into(),
        alt_allele: "C".into(),
    };
    assert!(is_cluster_ac_snp(&e));
}
