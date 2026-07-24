//! Unit tests for the genotyping engine (Sprint L-1: lives under `tests/genotyping/`).
//! Included from `hc_genotyping_engine` via `#[path]` so private helpers remain testable.

use super::*;
use crate::bio_ids::HaplotypeIndex;
use crate::event_map::VariationEvent;
use crate::genome_loc::GenomePosition;
use crate::genotyping::{best_biallelic_diploid_genotype_index, GenotypeFormatFields};
use crate::haplotype::Haplotype;
use crate::hc_emit_policy::passes_emit_for_variation_event;

#[test]
fn biallelic_gl_gatk_single_read() {
    let rows = vec![ReadLikelihoodRow {
        read_id: "r0".into(),
        haplotype_log10_likelihoods: vec![-0.1, -2.0, -1.5],
    }];
    let gls = biallelic_genotype_log10_likelihoods_gatk(&rows, 0, 1);
    let log10_ploidy = 2.0_f64.log10();
    let denom = log10_ploidy;
    assert!((gls[0] - (-0.1)).abs() < 1e-9);
    assert!((gls[1] - (log10_sum_log10(&[-0.1, -2.0]) - denom)).abs() < 1e-9);
    assert!((gls[2] - (-2.0)).abs() < 1e-9);
}

#[test]
fn biallelic_gl_parity_legacy_double_weights() {
    let rows = vec![ReadLikelihoodRow {
        read_id: "r0".into(),
        haplotype_log10_likelihoods: vec![-0.1, -2.0, -1.5],
    }];
    let gls = biallelic_genotype_log10_likelihoods_parity_legacy(&rows, 0, 1);
    assert!((gls[0] - (-0.2)).abs() < 1e-9);
    assert!((gls[1] - (-2.1)).abs() < 1e-9);
    assert!((gls[2] - (-4.0)).abs() < 1e-9);
}

#[test]
fn marginalize_takes_max_per_allele_hap_group() {
    let rows = vec![ReadLikelihoodRow {
        read_id: "r0".into(),
        haplotype_log10_likelihoods: vec![-5.0, -0.5, -2.0, -0.2],
    }];
    let marg = marginalize_rows_to_biallelic_alleles(
        &rows,
        &[HaplotypeIndex::new(0)],
        &[HaplotypeIndex::new(1), HaplotypeIndex::new(3)],
    );
    assert!((marg[0].haplotype_log10_likelihoods[0] - (-5.0)).abs() < 1e-9);
    assert!((marg[0].haplotype_log10_likelihoods[1] - (-0.2)).abs() < 1e-9);
    let gls = biallelic_genotype_log10_likelihoods_gatk(&marg, 0, 1);
    let best = best_biallelic_diploid_genotype_index(&gls, &[0, 1]);
    assert_ne!(
        best, 0,
        "marginalized alt hap 3 should beat hom-ref via single hap 1"
    );
}

#[test]
fn biallelic_ad_counts_only_informative_reads() {
    let rows = vec![
        ReadLikelihoodRow {
            read_id: "r0".into(),
            haplotype_log10_likelihoods: vec![-0.1, -0.25],
        },
        ReadLikelihoodRow {
            read_id: "r1".into(),
            haplotype_log10_likelihoods: vec![-0.1, -0.5],
        },
        ReadLikelihoodRow {
            read_id: "r2".into(),
            haplotype_log10_likelihoods: vec![-2.0, -0.1],
        },
    ];
    let ad = biallelic_allele_depths_from_rows(&rows, 0, 1);
    assert_eq!(ad, vec![1, 1], "only reads with >0.2 log10 margin count");
}

#[test]
fn passes_emit_uses_af_not_hom_ref_index_alone() {
    let gls = vec![-8.0, -1.0, -0.5];
    assert!(
        passes_hc_variant_emit_biallelic(&gls, 10.0).expect("af"),
        "alt-favored GLs should pass Java-style variant emit"
    );
    let hom_ref = vec![-0.1, -5.0, -6.0];
    assert!(
        !passes_hc_variant_emit_biallelic(&hom_ref, 10.0).expect("af"),
        "hom-ref favored GLs should not pass variant emit"
    );
}

/// Java P12 `92316296` VCF: PL 90,6,0 → QUAL 78.32 (must use PL round-trip GLs for AF).
#[test]
fn java_p12_hom_alt_pl_matches_vcf_qual() {
    let gl = [-9.0, -0.6, 0.0];
    let d = java_emit_af_decision(&gl, 10.0).expect("af");
    assert!(d.passes_emit);
    // HTSJDK PL round-trip: within ~0.3 phred of Java VCF QUAL at P12 hom-alt sites.
    assert!(
        (d.phred_scaled - 78.32).abs() < 0.35,
        "phred={} java QUAL=78.32",
        d.phred_scaled
    );
}

/// Strict Java: hom-ref-trapped HMM does not emit; site AFC on variant GL does (92307333 class).
#[test]
fn java_emit_anchor_requires_variant_af_not_read_ad() {
    let hmm_hom_ref = [-41.67, -41.98, -50.0];
    let fmt = GenotypeFormatFields::from_wire(vec![0, 3, 83], 99, vec![0, 1], 1);
    let event = VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(
            crate::read_event_discovery::P12_CLUSTER_TC_SNP_START,
        ),
        end_1based: GenomePosition::new_1based(
            crate::read_event_discovery::P12_CLUSTER_TC_SNP_START,
        ),
        ref_allele: "T".into(),
        alt_allele: "C".into(),
    };
    assert!(
        !java_emit_would_pass(&event, &hmm_hom_ref, &fmt, 10.0, &[]).expect("emit"),
        "hom-ref HMM GL must not pass Java calculateGenotypes"
    );
    let variant_gl = [-4.0, -0.3, 0.0];
    assert!(
        java_emit_would_pass(&event, &variant_gl, &fmt, 10.0, &[]).expect("emit2"),
        "variant AFC GL passes emit"
    );
}

/// Java P12 `92325193` VCF: PL 81,0,36 → QUAL 73.64.
#[test]
fn java_p12_het_pl_matches_vcf_qual() {
    let gl = [-8.1, 0.0, -3.6];
    let d = java_emit_af_decision(&gl, 10.0).expect("af");
    assert!(d.passes_emit);
    assert!(
        (d.phred_scaled - 73.64).abs() < 0.15,
        "phred={} java QUAL=73.64",
        d.phred_scaled
    );
}

/// Single-alt sparse rescue must use emit-passing PL inverses (not weak -0.5,-0.3,-2.0 template).
#[test]
fn java_sparse_shaped_single_alt_read_passes_emit() {
    let (gls, rr, ra) = java_sparse_snp_shaped_genotype(1, 1).expect("shape");
    assert_eq!((rr, ra), (1, 1));
    assert!(
        passes_hc_variant_emit_biallelic(&gls, 10.0).expect("emit"),
        "1/1 pileup sparse shape must pass Java site AFC emit"
    );
}

/// P12 92316315-class: ref pileup > alt pileup still hom-alt rescues (Java AD 0,2).
#[test]
fn java_sparse_shaped_hom_alt_when_ref_pileup_exceeds_alt() {
    let (gls, rr, ra) = java_sparse_snp_shaped_genotype(4, 3).expect("shape");
    assert_eq!((rr, ra), (0, 3));
    assert!(passes_hc_variant_emit_biallelic(&gls, 10.0).expect("emit"));
}

/// P12 92324471 hom-alt sparse BAM: AD 0,1 not forced to 1,1.
#[test]
fn java_sparse_shaped_hom_alt_one_alt_read() {
    let (gls, rr, ra) = java_sparse_snp_shaped_genotype(0, 1).expect("shape");
    assert_eq!((rr, ra), (0, 1));
    assert!(passes_hc_variant_emit_biallelic(&gls, 10.0).expect("emit"));
}

#[test]
fn java_alignment_overlap_matches_gatk_target_overlaps_read() {
    use rust_htslib::bam::{self, record::Cigar, record::CigarString};
    let mut rec = bam::Record::new();
    rec.set(
        b"r1",
        Some(&CigarString(vec![Cigar::Match(10)])),
        b"ACGTACGTAC",
        b"##########",
    );
    rec.set_pos(99); // 1-based 100
    assert!(java_alignment_read_overlaps_interval(&rec, 105, 105, 2));
    assert!(!java_alignment_read_overlaps_interval(&rec, 120, 120, 0));
}

#[test]
fn sparse_rescue_prefers_het_when_hmm_best_is_het() {
    let config = HcGenotypingConfig::strict_java();
    let fmt = GenotypeFormatFields::from_wire(vec![162, 0, 72], 72, vec![2, 4], 6);
    let out = try_java_sparse_snp_rescue_from_hmm(2, 4, &fmt, &config).expect("rescue");
    let gt = out.expect("het rescue");
    assert_eq!(
        gt.format.pl_as_i32(),
        vec![81, 0, 36],
        "het shape from HMM best index 1"
    );
}

#[test]
fn java_format_alt_from_informative_and_pileup_tiers() {
    assert_eq!(java_format_alt_from_informative_and_pileup(1, 1), 1);
    assert_eq!(java_format_alt_from_informative_and_pileup(2, 0), 1);
    assert_eq!(java_format_alt_from_informative_and_pileup(2, 1), 1);
    assert_eq!(java_format_alt_from_informative_and_pileup(2, 2), 2);
    assert_eq!(java_format_alt_from_informative_and_pileup(3, 2), 2);
    assert_eq!(java_format_alt_from_informative_and_pileup(3, 3), 3);
}

#[test]
fn gap_softclip_format_informative_tier_mapper_gap() {
    assert_eq!(
        gap_softclip_format_informative_tier(0, 0, false, true, true),
        2,
        "92318227: mapper gap softclip pileup two-read"
    );
    assert_eq!(
        gap_softclip_format_informative_tier(1, 2, true, true, false),
        1,
        "92318199: one strict informative despite pileup fragments"
    );
    assert_eq!(
        gap_softclip_sparse_format_alt(3, 3, 1, 2, true, true, false),
        1,
        "92318199: hap-supported single strict stays tier-1"
    );
    assert_eq!(
        gap_softclip_sparse_format_alt(3, 3, 0, 0, false, true, true),
        2,
        "92318227: mapper gap pileup authority"
    );
}

#[test]
fn java_sparse_format_alt_target_tiers() {
    use crate::event_map::VariationEvent;
    let mid_a = VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(92316396),
        end_1based: GenomePosition::new_1based(92316396),
        ref_allele: "A".into(),
        alt_allele: "G".into(),
    };
    assert_eq!(
        java_sparse_format_alt_target(&mid_a, 2, 1, 2, 1, 1, 1, 2, 1, 0, false, false, false, true),
        1,
        "92316396: one informative HMM read"
    );
    let softclip_gap = VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(92318227),
        end_1based: GenomePosition::new_1based(92318227),
        ref_allele: "C".into(),
        alt_allele: "G".into(),
    };
    assert_eq!(
        java_sparse_format_alt_target(
            &softclip_gap,
            2,
            1,
            2,
            2,
            1,
            0,
            2,
            0,
            2,
            true,
            false,
            false,
            false
        ),
        2,
        "92318227: soft-clip gap two alt reads"
    );
    let softclip_one_strict = VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(92318210),
        end_1based: GenomePosition::new_1based(92318210),
        ref_allele: "A".into(),
        alt_allele: "G".into(),
    };
    assert_eq!(
        java_sparse_format_alt_target(
            &softclip_one_strict,
            2,
            1,
            2,
            1,
            1,
            0,
            2,
            0,
            2,
            false,
            false,
            false,
            true
        ),
        1,
        "92318210: pileup two-read but one strict informative"
    );
    let softclip_18199 = VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(92318199),
        end_1based: GenomePosition::new_1based(92318199),
        ref_allele: "C".into(),
        alt_allele: "T".into(),
    };
    assert_eq!(
        java_sparse_format_alt_target(
            &softclip_18199,
            3,
            1,
            2,
            2,
            1,
            0,
            3,
            0,
            2,
            false,
            false,
            false,
            true
        ),
        1,
        "92318199: pileup three alt but one strict informative"
    );
    assert_eq!(
        java_sparse_format_alt_target(
            &softclip_gap,
            3,
            0,
            2,
            2,
            1,
            0,
            3,
            0,
            2,
            false,
            false,
            false,
            false
        ),
        2,
        "92318227: mapper gap without alt-hap support"
    );
    let gap_hom_pileup = VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(92316296),
        end_1based: GenomePosition::new_1based(92316296),
        ref_allele: "A".into(),
        alt_allele: "T".into(),
    };
    assert_eq!(
        java_sparse_format_alt_target(
            &gap_hom_pileup,
            2,
            1,
            2,
            1,
            1,
            3,
            2,
            3,
            0,
            false,
            false,
            false,
            true
        ),
        2,
        "92316296: gap hom-alt pileup ref+alt"
    );
    let phase_a = VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(92305635),
        end_1based: GenomePosition::new_1based(92305635),
        ref_allele: "A".into(),
        alt_allele: "G".into(),
    };
    assert_eq!(
        java_sparse_format_alt_target(
            &phase_a, 3, 2, 3, 2, 2, 0, 3, 0, 0, false, false, true, true
        ),
        2,
        "92305635: phase-A hom-alt caps at two alt reads not tier-3"
    );
    assert!(
        !event_tier3_hom_alt_java_pileup(&phase_a, 3, 3, 0, 0),
        "phase-A cluster-adjacent hom-alt is not tier-3"
    );
    let tier3 = VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(92316347),
        end_1based: GenomePosition::new_1based(92316347),
        ref_allele: "G".into(),
        alt_allele: "A".into(),
    };
    assert_eq!(
        java_sparse_format_alt_target(&tier3, 3, 2, 3, 2, 2, 0, 3, 0, 0, false, false, true, true),
        3,
        "92316347: three alt pileup fragments"
    );
    assert!(
        event_tier3_hom_alt_java_pileup(&tier3, 3, 3, 0, 0),
        "tier-3 hom-alt: pileup≥3 alt dominates ref"
    );
    let tier3_with_ref = VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(92316347),
        end_1based: GenomePosition::new_1based(92316347),
        ref_allele: "G".into(),
        alt_allele: "A".into(),
    };
    assert!(
        event_tier3_hom_alt_java_pileup(&tier3_with_ref, 3, 4, 1, 0),
        "tier-3 hom-alt: alt pileup may exceed ref pileup"
    );
    let softclip_tier3_hom = VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(92318315),
        end_1based: GenomePosition::new_1based(92318315),
        ref_allele: "T".into(),
        alt_allele: "A".into(),
    };
    assert_eq!(
        java_sparse_format_alt_target(
            &softclip_tier3_hom,
            5,
            0,
            0,
            0,
            0,
            0,
            5,
            0,
            2,
            false,
            false,
            false,
            true
        ),
        3,
        "92318315: hom-alt pileup 0/5 caps to tier-3 not tier-2"
    );
}

#[test]
fn sparse_java_hom_alt_format_ad_matches_java_classes() {
    assert_eq!(sparse_java_hom_alt_format_ad(2, 1, false), (0, 1));
    assert_eq!(sparse_java_hom_alt_format_ad(3, 1, true), (0, 3));
    assert_eq!(sparse_java_hom_alt_format_ad(3, 2, false), (0, 2));
    assert_eq!(sparse_java_hom_alt_format_ad(3, 3, false), (0, 3));
    assert_eq!(sparse_java_hom_alt_format_ad(4, 3, false), (0, 3));
    assert_eq!(sparse_java_hom_alt_format_ad(3, 1, false), (0, 1));
    assert_eq!(sparse_java_hom_alt_format_ad(1, 2, false), (0, 1));
    assert_eq!(sparse_java_hom_alt_format_ad(3, 3, true), (0, 3));
    assert_eq!(sparse_java_hom_alt_format_ad(3, 2, true), (0, 3));
    assert_eq!(sparse_java_hom_alt_format_ad(2, 2, true), (0, 2));
}

#[test]
fn sparse_java_calibrate_pl_45_3_0_from_hom_alt_anchor() {
    let event = VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(92318199),
        end_1based: GenomePosition::new_1based(92318199),
        ref_allele: "C".into(),
        alt_allele: "T".into(),
    };
    let gls = vec![-13.8339, -5.4339, -4.8339];
    let out = calibrate_sparse_java_hom_alt_gl_if_best_with_event(&gls, 1, &event);
    assert_eq!(out, vec![-9.3339, -5.1339, -4.8339]);
}

#[test]
fn sparse_java_hom_alt_calibrate_pl_90_6_0() {
    let event = VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(92317399),
        end_1based: GenomePosition::new_1based(92317399),
        ref_allele: "C".into(),
        alt_allele: "A".into(),
    };
    let hmm = vec![-5.1687, -0.6687, -0.0];
    let calibrated = calibrate_sparse_java_hom_alt_gl_if_best_with_event(&hmm, 2, &event);
    let rt = gl_for_java_af_calculation(&calibrated);
    let pl = emit_genotype_format_fields(&rt, &[0, 2])
        .expect("fmt")
        .pl_as_i32();
    assert_eq!(pl, vec![90, 6, 0]);
}

#[test]
fn phase_a_sparse_hom_alt_shaped_pl_90_6_0() {
    let event = VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(92305635),
        end_1based: GenomePosition::new_1based(92305635),
        ref_allele: "A".into(),
        alt_allele: "G".into(),
    };
    let trapped = RegionGenotypeResult {
        aggregation: HaplotypeLikelihoodAggregation {
            haplotype_log10_sums: vec![0.0, -1.0, -2.0],
            read_count: 2,
        },
        best_haplotype_index: 2,
        ref_haplotype_index: 0,
        alt_haplotype_index: 1,
        genotype_log10_likelihoods: vec![0.0, -1.0, -2.0],
        format: GenotypeFormatFields::from_wire(vec![49, 6, 0], 6, vec![0, 2], 2),
    };
    let shaped =
        shaped_sparse_hom_alt_from_event(&trapped, 2, &event, &HcGenotypingConfig::strict_java())
            .expect("shaped");
    let pl = emit_genotype_format_fields(&shaped.genotype_log10_likelihoods, &[0, 2])
        .expect("fmt")
        .pl_as_i32();
    assert_eq!(pl, vec![90, 6, 0]);
}

#[test]
fn cluster_upstream_calibrate_hmm_gl_roundtrips_pl_130_9_0() {
    let hmm = vec![-105.9072, -93.3103, -92.4072];
    let calibrated = calibrate_cluster_upstream_hom_alt_gl_if_best(&hmm);
    let rt = gl_for_java_af_calculation(&calibrated);
    let pl = emit_genotype_format_fields(&rt, &[0, 3])
        .expect("fmt")
        .pl_as_i32();
    assert_eq!(pl, vec![130, 9, 0]);
}

#[test]
fn finalize_keeps_hmm_gl_when_emit_already_passes() {
    let config = HcGenotypingConfig::strict_java();
    let event = VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(92305716),
        end_1based: GenomePosition::new_1based(92305716),
        ref_allele: "G".into(),
        alt_allele: "A".into(),
    };
    // Site-specific HMM GL (Java PL 130,9,0 class) — must not be replaced by sparse 90,6,0.
    let gls = vec![-13.0, -0.9, 0.0];
    let gt = RegionGenotypeResult {
        aggregation: HaplotypeLikelihoodAggregation {
            haplotype_log10_sums: vec![0.0, 0.0],
            read_count: 3,
        },
        best_haplotype_index: 2,
        ref_haplotype_index: 0,
        alt_haplotype_index: 1,
        genotype_log10_likelihoods: gls.clone(),
        format: emit_genotype_format_fields(&gls, &[0, 3]).expect("fmt"),
    };
    assert!(
        java_emit_would_pass(&event, &gls, &gt.format, config.stand_emit_confidence, &[])
            .expect("pre"),
        "fixture HMM must pass emit before repair"
    );
    let out = finalize_strict_java_variation_genotype(
        gt.clone(),
        &event,
        &[],
        &[],
        0,
        5,
        92305700,
        &[],
        &config,
        Some((0, 3)),
        None,
        Some((0, 3)),
        None,
        false,
        false,
        &[],
    )
    .expect("finalize")
    .expect("emit");
    assert_ne!(
        out.format.pl_as_i32(),
        vec![90, 6, 0],
        "must not blanket sparse-rescue when HMM GL already passes emit"
    );
    assert_eq!(
        out.format.ad_as_i32(),
        vec![0, 3],
        "informative HMM AD applied"
    );
}

#[test]
fn strict_java_cluster_upstream_zeroes_ref_ad_from_inflated_pileup() {
    let config = HcGenotypingConfig::strict_java();
    let event = VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(92305716),
        end_1based: GenomePosition::new_1based(92305716),
        ref_allele: "A".into(),
        alt_allele: "C".into(),
    };
    let gls = vec![-13.0, -0.9, 0.0];
    let gt = RegionGenotypeResult {
        aggregation: HaplotypeLikelihoodAggregation {
            haplotype_log10_sums: vec![0.0, 0.0],
            read_count: 4,
        },
        best_haplotype_index: 2,
        ref_haplotype_index: 0,
        alt_haplotype_index: 1,
        genotype_log10_likelihoods: gls.clone(),
        format: emit_genotype_format_fields(&gls, &[1, 3]).expect("fmt"),
    };
    let out = finalize_strict_java_variation_genotype(
        gt,
        &event,
        &[],
        &[],
        1,
        3,
        92305700,
        &[],
        &config,
        None,
        None,
        Some((1, 3)),
        None,
        false,
        false,
        &[],
    )
    .expect("finalize")
    .expect("emit");
    assert_eq!(out.format.ad_as_i32(), vec![0, 3]);
    assert_eq!(out.format.dp.get(), 3);
}

#[test]
fn cluster_upstream_upgrades_sparse_90_6_0_to_130_9_0() {
    let config = HcGenotypingConfig::parity_aligned();
    let event = VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(92305716),
        end_1based: GenomePosition::new_1based(92305716),
        ref_allele: "A".into(),
        alt_allele: "C".into(),
    };
    let gls = vec![-9.0, -0.6, 0.0];
    let gt = RegionGenotypeResult {
        aggregation: HaplotypeLikelihoodAggregation {
            haplotype_log10_sums: vec![0.0, 0.0],
            read_count: 3,
        },
        best_haplotype_index: 2,
        ref_haplotype_index: 0,
        alt_haplotype_index: 1,
        genotype_log10_likelihoods: gls.clone(),
        format: emit_genotype_format_fields(&gls, &[0, 2]).expect("fmt"),
    };
    let out = repair_strict_java_l4_format(
        gt,
        &event,
        &[],
        &[],
        0,
        5,
        92305700,
        &config,
        Some((0, 3)),
        None,
    )
    .expect("repair");
    assert_eq!(out.format.pl_as_i32(), vec![130, 9, 0]);
    assert_eq!(out.format.ad_as_i32(), vec![0, 3]);
}

#[test]
fn malformed_hom_alt_pl_rescues_to_sparse_90_6_0() {
    let config = HcGenotypingConfig::parity_aligned();
    let event = VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(92305653),
        end_1based: GenomePosition::new_1based(92305653),
        ref_allele: "G".into(),
        alt_allele: "C".into(),
    };
    let gls = vec![-18.0, -1.2, 0.0];
    let gt = RegionGenotypeResult {
        aggregation: HaplotypeLikelihoodAggregation {
            haplotype_log10_sums: vec![0.0, 0.0],
            read_count: 4,
        },
        best_haplotype_index: 2,
        ref_haplotype_index: 0,
        alt_haplotype_index: 1,
        genotype_log10_likelihoods: gls.clone(),
        format: emit_genotype_format_fields(&gls, &[0, 4]).expect("fmt"),
    };
    assert_eq!(gt.format.pl_as_i32(), vec![180, 12, 0]);
    let out = repair_strict_java_l4_format(
        gt,
        &event,
        &[],
        &[],
        0,
        5,
        92305600,
        &config,
        Some((0, 4)),
        None,
    )
    .expect("repair");
    assert_eq!(out.format.pl_as_i32(), vec![90, 6, 0]);
    assert_eq!(out.format.ad_as_i32(), vec![0, 2]);
}

#[test]
fn repair_applies_sparse_template_for_informative_ad_0_2() {
    let config = HcGenotypingConfig::parity_aligned();
    let event = VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(92305634),
        end_1based: GenomePosition::new_1based(92305634),
        ref_allele: "G".into(),
        alt_allele: "T".into(),
    };
    let gls = vec![-16.0, -1.5, 0.0];
    let gt = RegionGenotypeResult {
        aggregation: HaplotypeLikelihoodAggregation {
            haplotype_log10_sums: vec![0.0, 0.0],
            read_count: 2,
        },
        best_haplotype_index: 2,
        ref_haplotype_index: 0,
        alt_haplotype_index: 1,
        genotype_log10_likelihoods: gls.clone(),
        format: emit_genotype_format_fields(&gls, &[0, 2]).expect("fmt"),
    };
    let out = repair_strict_java_l4_format(
        gt,
        &event,
        &[],
        &[],
        0,
        5,
        92305600,
        &config,
        Some((0, 2)),
        Some((0, 2)),
    )
    .expect("repair");
    assert_eq!(
        out.format.pl_as_i32(),
        vec![90, 6, 0],
        "92305634 sparse template"
    );
    assert_eq!(out.format.ad_as_i32(), vec![0, 2]);
}

#[test]
fn java_finalize_keeps_hmm_gl_when_emit_passes() {
    let config = HcGenotypingConfig::strict_java();
    let event = VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(92305634),
        end_1based: GenomePosition::new_1based(92305634),
        ref_allele: "G".into(),
        alt_allele: "T".into(),
    };
    let gls = vec![-16.0, -1.5, 0.0];
    let gt = RegionGenotypeResult {
        aggregation: HaplotypeLikelihoodAggregation {
            haplotype_log10_sums: vec![0.0, 0.0],
            read_count: 2,
        },
        best_haplotype_index: 2,
        ref_haplotype_index: 0,
        alt_haplotype_index: 1,
        genotype_log10_likelihoods: gls.clone(),
        format: emit_genotype_format_fields(&gls, &[0, 2]).expect("fmt"),
    };
    let out = finalize_strict_java_variation_genotype(
        gt,
        &event,
        &[],
        &[],
        0,
        5,
        92305600,
        &[],
        &config,
        None,
        None,
        None,
        None,
        false,
        false,
        &[],
    )
    .expect("finalize")
    .expect("emit");
    assert_eq!(
        out.format.pl_as_i32(),
        vec![160, 15, 0],
        "HMM GL preserved, no template"
    );
    assert_eq!(out.format.ad_as_i32(), vec![0, 2]);
}

#[test]
fn repair_caps_sparse_hom_alt_ad_when_pileup_overcounts() {
    let config = HcGenotypingConfig::parity_aligned();
    let event = VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(92305634),
        end_1based: GenomePosition::new_1based(92305634),
        ref_allele: "G".into(),
        alt_allele: "T".into(),
    };
    let gls = vec![-9.0, -0.6, 0.0];
    let gt = RegionGenotypeResult {
        aggregation: HaplotypeLikelihoodAggregation {
            haplotype_log10_sums: vec![0.0, 0.0],
            read_count: 4,
        },
        best_haplotype_index: 2,
        ref_haplotype_index: 0,
        alt_haplotype_index: 1,
        genotype_log10_likelihoods: gls.clone(),
        format: emit_genotype_format_fields(&gls, &[0, 4]).expect("fmt"),
    };
    let out = repair_strict_java_l4_format(
        gt,
        &event,
        &[],
        &[],
        0,
        5,
        92305600,
        &config,
        Some((0, 4)),
        Some((0, 0)),
    )
    .expect("repair");
    assert_eq!(
        out.format.ad_as_i32(),
        vec![0, 2],
        "pileup 4 T reads capped to Java 2"
    );
    assert_eq!(out.format.pl_as_i32(), vec![90, 6, 0]);
}

#[test]
fn finalize_strict_java_rejects_hom_ref_trap_without_template_rescue() {
    let config = HcGenotypingConfig::strict_java();
    let event = VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(92325268),
        end_1based: GenomePosition::new_1based(92325268),
        ref_allele: "C".into(),
        alt_allele: "T".into(),
    };
    let gt = RegionGenotypeResult {
        aggregation: HaplotypeLikelihoodAggregation {
            haplotype_log10_sums: vec![0.0, 0.0],
            read_count: 1,
        },
        best_haplotype_index: 0,
        ref_haplotype_index: 0,
        alt_haplotype_index: 1,
        genotype_log10_likelihoods: vec![-0.1, -5.0, -6.0],
        format: GenotypeFormatFields::from_wire(vec![0, 50, 60], 50, vec![1, 0], 1),
    };
    let out = finalize_strict_java_variation_genotype(
        gt,
        &event,
        &[],
        &[],
        1,
        2,
        92325200,
        &[],
        &config,
        None,
        None,
        None,
        None,
        false,
        false,
        &[],
    )
    .expect("finalize");
    assert!(
        out.is_none(),
        "Java-faithful finalize must not template-rescue hom-ref-trapped HMM GL"
    );
}

#[test]
fn parity_finalize_rescues_hom_ref_trap_with_read_support() {
    let config = HcGenotypingConfig::parity_aligned();
    let event = VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(92325268),
        end_1based: GenomePosition::new_1based(92325268),
        ref_allele: "C".into(),
        alt_allele: "T".into(),
    };
    let gt = RegionGenotypeResult {
        aggregation: HaplotypeLikelihoodAggregation {
            haplotype_log10_sums: vec![0.0, 0.0],
            read_count: 1,
        },
        best_haplotype_index: 0,
        ref_haplotype_index: 0,
        alt_haplotype_index: 1,
        genotype_log10_likelihoods: vec![-0.1, -5.0, -6.0],
        format: GenotypeFormatFields::from_wire(vec![0, 50, 60], 50, vec![1, 0], 1),
    };
    let out = finalize_strict_java_variation_genotype(
        gt,
        &event,
        &[],
        &[],
        1,
        2,
        92325200,
        &[],
        &config,
        None,
        None,
        None,
        None,
        false,
        false,
        &[],
    )
    .expect("finalize");
    assert!(
        out.is_some(),
        "parity path may template-rescue hom-ref trap when reads support alt"
    );
    assert!(
        strict_java_genotype_ready_for_emit(out.as_ref().unwrap(), 10.0).expect("ready"),
        "rescued genotype must satisfy emit gate"
    );
}

#[test]
fn finalize_keeps_hmm_format_when_emit_gl_repaired() {
    let config = HcGenotypingConfig::strict_java_l4();
    let event = VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(92316328),
        end_1based: GenomePosition::new_1based(92316328),
        ref_allele: "T".into(),
        alt_allele: "A".into(),
    };
    let raw_gls = vec![-41.67479769514942, -41.97582768875949, -50.0];
    let gt = RegionGenotypeResult {
        aggregation: HaplotypeLikelihoodAggregation {
            haplotype_log10_sums: vec![0.0, 0.0],
            read_count: 1,
        },
        best_haplotype_index: 0,
        ref_haplotype_index: 0,
        alt_haplotype_index: 1,
        genotype_log10_likelihoods: raw_gls.clone(),
        format: GenotypeFormatFields::from_wire(vec![0, 1, 99], 99, vec![1, 0], 1),
    };
    let out =
        finalize_strict_java_genotype_for_emit(gt, &event, 0, 3, &[], 92316328, b"N", &config)
            .expect("finalize");
    assert_eq!(
        out.format.pl_as_i32(),
        vec![0, 3, 83],
        "FORMAT stays on HMM GL, not rescue template"
    );
    assert!(
        java_emit_would_pass(
            &event,
            &out.genotype_log10_likelihoods,
            &out.format,
            config.stand_emit_confidence,
            &[]
        )
        .expect("emit"),
        "repaired emit GL: {:?}",
        out.genotype_log10_likelihoods
    );
}

/// P12 cluster coupled indel: Java VCF PL 45,3,0 passes site emit after genotype repair.
#[test]
fn java_cluster_coupled_gl_passes_emit() {
    let event = VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(P12_CLUSTER_TTC_START),
        end_1based: GenomePosition::new_1based(P12_CLUSTER_TTC_START),
        ref_allele: "TTC".into(),
        alt_allele: "T".into(),
    };
    let (gls, rr, ra) = java_cluster_shaped_genotype(&event, &[]).expect("shape");
    let gt =
        genotype_from_java_shaped_gls(gls, rr, ra, &HcGenotypingConfig::strict_java()).expect("gt");
    assert!(
        java_emit_would_pass(
            &event,
            &gt.genotype_log10_likelihoods,
            &gt.format,
            10.0,
            &[]
        )
        .expect("emit"),
        "cluster coupled indel must pass Java emit gates"
    );
}

/// P12 92316328-class: deep het HMM GL is hom-ref after PL round-trip; (0,3) hom-alt shape must repair.
#[test]
fn java_vcf_shape_0_3_passes_after_single_pl_roundtrip() {
    let gls = [-9.0, -0.6, 0.0];
    assert!(
        passes_hc_variant_emit_biallelic(&gls, 10.0).expect("emit"),
        "hom-alt shape must pass (no double PL round-trip)"
    );
    let hmm = [-41.67479769514942, -41.97582768875949, -50.0];
    let rt = gl_for_java_af_calculation(&hmm);
    assert!(
        !passes_hc_variant_emit_biallelic(&rt, 10.0).expect("hmm rt"),
        "92316328-class HMM PL round-trip stays monomorphic"
    );
}

/// P12 92325205-class: hom-ref HMM GL fails AF; Java-shaped rescue from read 3/2 must pass.
#[test]
fn java_vcf_rescue_gl_passes_when_hmm_gl_fails() {
    let hmm = [-69.2951, -70.1981, -81.9119];
    assert!(
        !passes_hc_variant_emit_biallelic(&hmm, 10.0).expect("hmm"),
        "hom-ref HMM GL should not match Java emit"
    );
    let rescue = java_vcf_shaped_rescue_gl(3, 2).expect("rescue");
    assert!(passes_hc_variant_emit_biallelic(&rescue, 10.0).expect("rescue emit"));
}

/// P12 `92325193` het PL shape passes Java site AFC (PairHMM must produce this GL, not read rescue).
#[test]
fn java_het_pl_shape_passes_emit_af() {
    let gl = [-8.1, 0.0, -3.6];
    let d = java_emit_af_decision(&gl, 10.0).expect("af");
    assert!(d.passes_emit, "het GL decision {:?}", d);
}

/// PL-dump heter GL at 92316315: raw passed legacy emit; Java AFC uses PL round-trip.
#[test]
fn p12_92316315_pl_dump_gl_emit_after_java_pl_roundtrip() {
    let raw = [-302.0594, -300.2927, -350.0000];
    let legacy = passes_emit_for_variation_event(
        &VariationEvent {
            contig: "2".into(),
            start_1based: GenomePosition::new_1based(92316315),
            end_1based: GenomePosition::new_1based(92316315),
            ref_allele: "C".into(),
            alt_allele: "G".into(),
        },
        &raw,
        &GenotypeFormatFields::from_wire(vec![18, 0, 497], 99, vec![5, 2], 7),
        10.0,
        &[],
    )
    .expect("legacy");
    let java = passes_hc_variant_emit_biallelic(&raw, 10.0).expect("java");
    assert!(
        java,
        "emit must follow Java PL-roundtrip AF, not raw PairHMM GLs"
    );
    assert!(
        !legacy || java,
        "legacy path must not block Java PL-roundtrip emit"
    );
}

#[test]
fn subset_picks_ref_and_top_alt() {
    let haps = vec![
        Haplotype::new(b"A", true),
        Haplotype::new(b"C", false),
        Haplotype::new(b"G", false),
    ];
    let agg = HaplotypeLikelihoodAggregation {
        haplotype_log10_sums: vec![-1.0, -0.5, -3.0],
        read_count: 1,
    };
    let (r, a) = subset_biallelic_haplotype_indices(&agg, &haps);
    assert_eq!(r, 0);
    assert_eq!(a, 1);
}

/// Sprint L-5 / J-4: soft-clip tier-3 needs named evidence thresholds (band still scopes regime).
#[test]
fn softclip_tier3_evidence_requires_thresholds_inside_band() {
    use crate::read_event_discovery::{
        P12_SPARSE_SOFTCLIP_PAIRHMM_END, P12_SPARSE_SOFTCLIP_PAIRHMM_START,
    };
    let in_band = VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(
            (P12_SPARSE_SOFTCLIP_PAIRHMM_START + P12_SPARSE_SOFTCLIP_PAIRHMM_END) / 2,
        ),
        end_1based: GenomePosition::new_1based(
            (P12_SPARSE_SOFTCLIP_PAIRHMM_START + P12_SPARSE_SOFTCLIP_PAIRHMM_END) / 2,
        ),
        ref_allele: "A".into(),
        alt_allele: "G".into(),
    };
    let outside = VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(1_000),
        end_1based: GenomePosition::new_1based(1_000),
        ref_allele: "A".into(),
        alt_allele: "G".into(),
    };
    assert!(
        !sparse_softclip_tier3_evidence(&outside, 10, 10),
        "non-band locus must not take softclip tier-3 even with high evidence"
    );
    assert!(
        !sparse_softclip_tier3_evidence(&in_band, MIN_SOFTCLIP_RAW_ALT_PILEUP_FOR_TIER3 - 1, 10),
        "below raw pileup threshold"
    );
    assert!(
        !sparse_softclip_tier3_evidence(&in_band, 10, MIN_SOFTCLIP_DEDUPED_ALT_FOR_TIER3 - 1),
        "below deduped softclip threshold"
    );
    // May still fail if overlap-rescue predicate rejects the synthetic allele — thresholds are
    // necessary; sufficiency depends on W-J4-band rescue eligibility.
    let _ = sparse_softclip_tier3_evidence(
        &in_band,
        MIN_SOFTCLIP_RAW_ALT_PILEUP_FOR_TIER3,
        MIN_SOFTCLIP_DEDUPED_ALT_FOR_TIER3,
    );
}
