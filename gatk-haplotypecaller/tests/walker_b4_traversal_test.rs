//! Walker traversal API smoke (multi-span shard + apply stats).

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::{
    traverse_assembly_region_walker, ReadFilterParams, WalkerApplyStats, WalkerTraversalConfig,
    GATK_DEFAULT_ASSEMBLY_REGION_PADDING, GATK_DEFAULT_MAX_READS_PER_ALIGNMENT_START,
};

#[test]
fn walker_traversal_disjoint_spans_matches_b3_apply_contract() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../parity/fixtures");
    let ref_fa = root.join("reference.fa");
    let bam = root.join("sample.bam");
    let dict = SequenceDictionary::from_fasta_path(&ref_fa).unwrap();
    let specs = parse_intervals_cli_string(&dict, "chr1:1-5;chr1:20-25").unwrap();
    let filters = ReadFilterParams::default();
    let walk = traverse_assembly_region_walker(
        &dict,
        &specs,
        &ref_fa,
        &bam,
        &filters,
        &WalkerTraversalConfig::gatk_haplotype_caller_production(5),
    )
    .unwrap();
    assert_eq!(walk.shards.len(), 1);
    assert_eq!(walk.shards[0].regions.len(), 2);
    assert_eq!(
        walk.apply_stats,
        WalkerApplyStats {
            total_apply: 2,
            inactive_fast_path: 2,
            active_full: 0,
        }
    );
}

#[test]
fn production_traversal_config_enables_positional_ds_cap_50() {
    let cfg = WalkerTraversalConfig::gatk_haplotype_caller_production(
        GATK_DEFAULT_ASSEMBLY_REGION_PADDING,
    );
    assert_eq!(
        cfg.shard_pipeline.downsample.max_reads_per_alignment_start,
        GATK_DEFAULT_MAX_READS_PER_ALIGNMENT_START
    );
    assert_eq!(
        cfg.downsample.max_reads_per_alignment_start,
        GATK_DEFAULT_MAX_READS_PER_ALIGNMENT_START
    );
    assert!(cfg.shard_pipeline.apply_iupac_pre_transform);
}

#[test]
fn walker_traversal_default_padding_chr1_5_15() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../parity/fixtures");
    let ref_fa = root.join("reference.fa");
    let bam = root.join("sample.bam");
    let dict = SequenceDictionary::from_fasta_path(&ref_fa).unwrap();
    let specs = parse_intervals_cli_string(&dict, "chr1:5-15").unwrap();
    let walk = traverse_assembly_region_walker(
        &dict,
        &specs,
        &ref_fa,
        &bam,
        &ReadFilterParams::default(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(
            GATK_DEFAULT_ASSEMBLY_REGION_PADDING,
        ),
    )
    .unwrap();
    assert_eq!(walk.apply_stats.total_apply, 1);
    assert_eq!(walk.apply_stats.inactive_fast_path, 1);
}
