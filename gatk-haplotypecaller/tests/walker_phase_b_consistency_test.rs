//! B.3 apply-summary == B.4 traversal; B.2 region list via traversal.

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::{
    collect_assembly_regions, traverse_assembly_region_walker, ReadFilterParams,
    WalkerTraversalConfig, GATK_DEFAULT_ASSEMBLY_REGION_PADDING,
};
use std::path::PathBuf;

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../parity/fixtures")
        .join(name)
}

#[test]
fn apply_summary_matches_traversal_on_b3_disjoint_fixture() {
    let ref_fa = fixture("reference.fa");
    let bam = fixture("sample.bam");
    let dict = SequenceDictionary::from_fasta_path(&ref_fa).unwrap();
    let specs = parse_intervals_cli_string(&dict, "chr1:1-5;chr1:20-25").unwrap();
    let filters = ReadFilterParams::default();
    let cfg = WalkerTraversalConfig::gatk_haplotype_caller_defaults(5);
    let walk =
        traverse_assembly_region_walker(&dict, &specs, &ref_fa, &bam, &filters, &cfg).unwrap();
    let flat = collect_assembly_regions(&dict, &specs, &ref_fa, &bam, &filters, &cfg).unwrap();
    assert_eq!(flat.len(), walk.apply_stats.total_apply);
    assert_eq!(walk.apply_stats.total_apply, 2);
    assert_eq!(walk.apply_stats.inactive_fast_path, 2);
}

#[test]
fn b2_region_count_matches_golden_chr1_5_15() {
    let ref_fa = fixture("reference.fa");
    let bam = fixture("sample.bam");
    let dict = SequenceDictionary::from_fasta_path(&ref_fa).unwrap();
    let specs = parse_intervals_cli_string(&dict, "chr1:5-15").unwrap();
    let regions = collect_assembly_regions(
        &dict,
        &specs,
        &ref_fa,
        &bam,
        &ReadFilterParams::default(),
        &WalkerTraversalConfig::gatk_haplotype_caller_defaults(
            GATK_DEFAULT_ASSEMBLY_REGION_PADDING,
        ),
    )
    .unwrap();
    assert_eq!(regions.len(), 1);
    assert!(!regions[0].is_active);
}
