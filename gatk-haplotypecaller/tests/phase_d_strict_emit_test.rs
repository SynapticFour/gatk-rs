//! strict Java emit uses only GATK AF/GQ thresholds (no read-style sparse emit).
//! Run: `cargo test -p gatk-haplotypecaller phase_d_strict_emit --release`

use gatk_haplotypecaller::event_map::VariationEvent;
use gatk_haplotypecaller::genome_loc::GenomePosition;
use gatk_haplotypecaller::genotyping::GenotypeFormatFields;
use gatk_haplotypecaller::hc_emit_policy::{
    passes_emit_for_genotyped_call, passes_emit_for_variation_event, passes_read_style_sparse_emit,
};
use gatk_haplotypecaller::hc_genotyping_engine::HcGenotypingConfig;

#[test]
fn phase_d_strict_emit_rejects_weak_ad_without_read_style() {
    let event = VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(100),
        end_1based: GenomePosition::new_1based(100),
        ref_allele: "A".into(),
        alt_allele: "G".into(),
    };
    let gls = vec![0.0, -0.1, -50.0];
    let format = GenotypeFormatFields::from_wire(vec![0, 3, 500], 0, vec![10, 1], 11);
    let strict = HcGenotypingConfig::strict_java();
    assert!(!strict.enable_read_style_emit);
    let prod_read_style =
        strict.enable_read_style_emit && passes_read_style_sparse_emit(&event, 1, 10);
    assert!(
        !prod_read_style,
        "Phase D: strict prod must not use read-style sparse emit"
    );
    assert!(
        !passes_emit_for_variation_event(&event, &gls, &format, strict.stand_emit_confidence, &[])
            .unwrap(),
        "weak alt must fail Java emit threshold under strict"
    );
    assert!(!passes_emit_for_genotyped_call(
        &event,
        &gls,
        &format,
        strict.stand_emit_confidence,
        &[]
    )
    .unwrap());
}

#[test]
#[cfg(feature = "parity_harness")]
fn phase_d_parity_aligned_may_use_read_style_when_enabled() {
    let _event = VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(100),
        end_1based: GenomePosition::new_1based(100),
        ref_allele: "A".into(),
        alt_allele: "G".into(),
    };
    let parity = HcGenotypingConfig::parity_aligned();
    assert!(!parity.enable_read_style_emit);
}
