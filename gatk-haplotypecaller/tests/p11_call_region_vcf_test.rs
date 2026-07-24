//! Live p11 fixture: strict `call_region` must emit chrLive:15 (Java j2-vcf row).

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::read_model::ReadFilterParams;
use gatk_haplotypecaller::walker_traversal::{
    flatten_assembly_regions, traverse_assembly_region_walker,
};
use gatk_haplotypecaller::{
    call_disposition, region_vcf_emit::DEFAULT_STAND_EMIT_CONFIDENCE,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, WalkerTraversalConfig,
};
use std::path::Path;

#[test]
fn p11_java_positive_strict_call_region_emits_snp() {
    let repo = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
    let ref_fa = repo.join("parity/fixtures/p5_live_reference.fa");
    let bam = repo.join("parity/build/sam-indexed-bam/p11_java_positive.bam");
    assert!(
        bam.is_file(),
        "run hc_full_parity_gate_dump once to build indexed BAM"
    );
    let dict = SequenceDictionary::from_fasta_path(&ref_fa).expect("dict");
    let specs = parse_intervals_cli_string(&dict, "chrLive:1-63").expect("interval");
    let filters = ReadFilterParams::gatk_standard_hc();
    let cfg = WalkerTraversalConfig::gatk_haplotype_caller_production(0);
    let walk = traverse_assembly_region_walker(&dict, &specs, &ref_fa, &bam, &filters, &cfg)
        .expect("walk");
    let regions = flatten_assembly_regions(&walk);
    let region = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            )
        })
        .expect("active region");
    let args = CallRegionArgs::strict_java();
    let outcome = HaplotypeCallerEngine::call_region(region, &dict, &ref_fa, &args)
        .expect("call_region")
        .expect("outcome");
    let rec =
        try_emit_call_region_variants(region, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("emit")
            .into_iter()
            .next()
            .expect("variant at chrLive:15");
    assert_eq!(rec.chromosome, "chrLive");
    assert_eq!(rec.position, 15);
    assert_eq!(rec.reference, "T");
    assert_eq!(rec.alternate[0], "A");
    assert!(rec.quality.unwrap() > 1000.0);
}
