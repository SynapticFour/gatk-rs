//! j2-vcf L2 fixtures: p5 must not emit; p11 must emit (Java golden).

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::read_model::ReadFilterParams;
use gatk_haplotypecaller::walker_traversal::{
    flatten_assembly_regions, traverse_assembly_region_walker,
};
use gatk_haplotypecaller::{
    call_disposition, region_vcf_emit::DEFAULT_STAND_EMIT_CONFIDENCE, try_emit_call_region_variant,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, WalkerTraversalConfig,
};
use std::path::Path;

fn first_active_outcome(
    _repo: &Path,
    ref_fa: &Path,
    bam: &Path,
    interval: &str,
) -> (
    gatk_haplotypecaller::assembly_region_iterator::AssemblyRegion,
    gatk_haplotypecaller::engine::CallRegionOutcome,
) {
    let dict = SequenceDictionary::from_fasta_path(ref_fa).expect("dict");
    let specs = parse_intervals_cli_string(&dict, interval).expect("interval");
    let walk = traverse_assembly_region_walker(
        &dict,
        &specs,
        ref_fa,
        bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(0),
    )
    .expect("walk");
    let regions = flatten_assembly_regions(&walk);
    let region = regions
        .into_iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            )
        })
        .expect("active region");
    let outcome =
        HaplotypeCallerEngine::call_region(&region, &dict, ref_fa, &CallRegionArgs::strict_java())
            .expect("call_region")
            .expect("outcome");
    (region, outcome)
}

#[test]
fn j2_vcf_p5_sparse_snp_does_not_emit() {
    let repo = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
    let ref_fa = repo.join("parity/fixtures/p5_live_reference.fa");
    let bam = repo.join("parity/build/sam-indexed-bam/p5_live_case_snp.bam");
    assert!(bam.is_file());
    let (region, outcome) = first_active_outcome(&repo, &ref_fa, &bam, "chrLive:1-24");
    assert!(
        try_emit_call_region_variants(&region, &outcome, "S", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("multi")
            .is_empty()
    );
    assert!(
        try_emit_call_region_variant(&region, &outcome, "S", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("single")
            .is_none()
    );
}

#[test]
fn j2_vcf_p11_java_positive_emits() {
    let repo = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
    let ref_fa = repo.join("parity/fixtures/p5_live_reference.fa");
    let bam = repo.join("parity/build/sam-indexed-bam/p11_java_positive.bam");
    assert!(bam.is_file());
    let (region, outcome) = first_active_outcome(&repo, &ref_fa, &bam, "chrLive:1-63");
    let multi =
        try_emit_call_region_variants(&region, &outcome, "S", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("multi");
    assert_eq!(multi.len(), 1);
    assert_eq!(multi[0].position, 15);
}
