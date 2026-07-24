//! P12 L5.2: active-region RCM must use post-realign genotyping reads (Java getPileupsOverReference).
//! Run: `P12_REFERENCE=parity/realworld/assets/hs37d5.simple.fa cargo test -p gatk-haplotypecaller p12_rcm_fringe --release -- --ignored --nocapture`

use gatk_core::reference::ReferenceWindowCache;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::ref_confidence::{
    reference_confidence_loci_for_active_region, reference_confidence_loci_for_region,
    ClusterRcmEvidenceMode, ReferenceConfidenceConfig,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions,
    reference_vcf_emit::{emitted_variant_starts_in_region, first_emitted_variant_start_in_region},
    traverse_assembly_region_walker, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
};
use rust_htslib::bam::Read;
use std::path::Path;

const FRINGE_POS: u64 = 92305554;

fn locus_gq(loci: &[gatk_haplotypecaller::genotyping::ReferenceConfidenceLocus], pos: u64) -> i32 {
    loci.iter()
        .find(|l| l.position_1based as u64 == pos)
        .expect("locus")
        .gq
}

fn p12_paths() -> Option<(std::path::PathBuf, std::path::PathBuf)> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
    let ref_path = std::env::var("P12_REFERENCE")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|_| root.join("parity/realworld/assets/hs37d5.simple.fa"));
    let ref_path = if ref_path.is_absolute() {
        ref_path
    } else {
        root.join(ref_path)
    };
    let bam = root.join("parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam");
    if ref_path.is_file() && bam.is_file() {
        Some((ref_path, bam))
    } else {
        None
    }
}

#[test]
#[ignore = "P12 BAM: fringe RCM pileup source (~2 min)"]
fn p12_rcm_fringe_uses_genotyping_reads_for_zero_gq() {
    let Some((ref_fasta, bam)) = p12_paths() else {
        eprintln!("skip: P12_REFERENCE / BAM");
        return;
    };

    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, "2:92305500-92305640").expect("interval");
    let walk = traverse_assembly_region_walker(
        &dict,
        &specs,
        &ref_fasta,
        &bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
    )
    .expect("walk");
    let regions = flatten_assembly_regions(&walk);
    let region = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= FRINGE_POS
                && r.end.get() >= FRINGE_POS
        })
        .expect("active region covering fringe");

    let header = rust_htslib::bam::Reader::from_path(&bam)
        .expect("bam")
        .header()
        .clone();
    let outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("outcome");

    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
    let config = ReferenceConfidenceConfig::default();
    let filters = ReadFilterParams::gatk_standard_hc();

    let raw_loci = reference_confidence_loci_for_region(
        region,
        &header,
        &config,
        &filters,
        &mut ref_cache,
        &dict,
    )
    .expect("raw loci");
    let first_variant = first_emitted_variant_start_in_region(region, &outcome, 10.0)
        .expect("first emitted variant");
    let emitted_variants =
        emitted_variant_starts_in_region(region, &outcome, 10.0).expect("emitted variants");
    let gt_loci = reference_confidence_loci_for_active_region(
        region,
        &outcome.genotyping_reads,
        first_variant,
        &emitted_variants,
        &header,
        &config,
        &filters,
        &mut ref_cache,
        &dict,
        ClusterRcmEvidenceMode::Production,
    )
    .expect("hybrid loci");
    let raw_gq = locus_gq(&raw_loci, FRINGE_POS);
    let gt_gq = locus_gq(&gt_loci, FRINGE_POS);

    eprintln!("fringe {FRINGE_POS}: raw gq={raw_gq}, hybrid active gq={gt_gq}");
    assert!(
        raw_gq > 0,
        "fixture: raw pileup non-zero at fringe (got {raw_gq})"
    );
    assert_eq!(
        gt_gq, 0,
        "hybrid active pileup GQ=0 at fringe (got {gt_gq})"
    );
}
