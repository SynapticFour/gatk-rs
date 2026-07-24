//! `call_region` uses Java order (filter → realign, CIGAR-only events in parity mode).

#![allow(deprecated)]
use gatk_haplotypecaller::event_map::VariationEvent;
use gatk_haplotypecaller::genome_loc::GenomePosition;
use gatk_haplotypecaller::genotyping::GenotypeFormatFields;
use gatk_haplotypecaller::hc_emit_policy::passes_emit_for_variation_event;
use gatk_haplotypecaller::hc_genotyping_engine::HcGenotypingConfig;
use gatk_haplotypecaller::CallRegionArgs;

#[test]
fn parity_call_region_args_disable_supplements_and_sparse_genotype() {
    let args = CallRegionArgs::strict_java();
    assert!(!args.enable_read_event_supplement);
    assert!(!args.genotyping.enable_sparse_read_genotype);
    assert!(!args.genotyping.enable_read_style_emit);
    assert!(!args.genotyping.genotype_stored_events_only);
    assert!(!args.genotyping.enable_l4_emit_gl_rescue);
    assert!(args.genotyping.enable_java_strict());
    assert_eq!(args.mode, gatk_haplotypecaller::CallRegionMode::StrictJava);
    assert!(args.is_strict_java());
    assert!(args.enable_allele_filtering);
    assert!(!args.enable_assembly_cluster_indel_inject);
    assert!(args.assemble.strict_java_assembly);
    assert!(args.assemble.assembler.dangling_java_exact);
}

#[test]
fn strict_java_genotyping_config_has_no_bridges() {
    let g = HcGenotypingConfig::strict_java();
    assert!(g.enable_java_strict());
    assert!(!g.enable_sparse_read_genotype);
    assert!(!g.enable_read_style_emit);
    assert!(!g.genotype_stored_events_only);
    assert!(!g.enable_l4_emit_gl_rescue);
}

#[test]
fn strict_java_emit_never_uses_read_style_sparse() {
    use gatk_haplotypecaller::hc_emit_policy::passes_read_style_sparse_emit;
    let event = VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(1),
        end_1based: GenomePosition::new_1based(1),
        ref_allele: "A".into(),
        alt_allele: "G".into(),
    };
    let g = HcGenotypingConfig::strict_java();
    assert!(!g.enable_read_style_emit);
    let prod_read_style = g.enable_read_style_emit && passes_read_style_sparse_emit(&event, 5, 0);
    assert!(
        !prod_read_style,
        "Phase D: strict prod must not use read-style sparse emit"
    );
}

#[cfg(feature = "parity_harness")]
#[test]
fn parity_aligned_keeps_transitional_assembly_hooks() {
    let args = CallRegionArgs::parity_aligned();
    assert!(!args.is_strict_java());
    assert!(!args.assemble.strict_java_assembly);
}

#[cfg(feature = "parity_harness")]
#[test]
fn phase_c_parity_emit_rejects_cluster_ad_bypass() {
    let event = VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(1),
        end_1based: GenomePosition::new_1based(1),
        ref_allele: "A".into(),
        alt_allele: "G".into(),
    };
    let gls = vec![0.0, -1.0, -50.0];
    let format = GenotypeFormatFields::from_wire(vec![], 0, vec![10, 1], 11);
    let config = HcGenotypingConfig::parity_aligned();
    let passes =
        passes_emit_for_variation_event(&event, &gls, &format, config.stand_emit_confidence, &[])
            .unwrap();
    assert!(!passes, "weak alt with high GQ bypass removed (Phase C)");
}
