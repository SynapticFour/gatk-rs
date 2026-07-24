//! mid-region 923162 strict genotyping — EventMap events must yield genotyped calls.
//! Run: `P12_REFERENCE=… cargo test -p gatk-haplotypecaller p12_region_923162_genotyping --release -- --nocapture`

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    AssemblyRegionCallDisposition, CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams,
    WalkerTraversalConfig,
};
use std::path::Path;

const REGION_LO: u64 = 92316227;
const REGION_HI: u64 = 92316475;

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
    if !ref_path.is_file() || !bam.is_file() {
        return None;
    }
    Some((ref_path, bam))
}

#[test]
fn p12_region_923162_genotyping_probe() {
    let Some((ref_fasta, bam)) = p12_paths() else {
        eprintln!("skip: set P12_REFERENCE");
        return;
    };
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, "2:92300000-92350000").expect("interval");
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
            ) && r.start.get() <= REGION_LO
                && r.end.get() >= REGION_HI
        })
        .expect("active region covering 92316227-92316475");
    let args = CallRegionArgs::strict_java();
    let outcome = HaplotypeCallerEngine::call_region(region, &dict, &ref_fasta, &args)
        .expect("call")
        .expect("strict_java call_region");
    eprintln!(
        "events={} genotyped_calls={} read_ll={}",
        outcome.assembly.variation_events().len(),
        outcome.genotyped_calls.len(),
        outcome.read_likelihoods.len()
    );
    for call in &outcome.genotyped_calls {
        eprintln!(
            "call\t{} {}/{} GQ={}",
            call.event.start_1based.get(),
            call.event.ref_allele,
            call.event.alt_allele,
            call.genotype.format.gq
        );
    }
    assert!(
        outcome.assembly.variation_events().len() >= 2,
        "Phase C: need EventMap events from ASM"
    );
    assert!(
        !outcome.read_likelihoods.is_empty(),
        "Phase C: PairHMM likelihoods required for genotyping"
    );
    assert!(
        !outcome.genotyped_calls.is_empty(),
        "Phase C: expected genotyped_calls from stored events + GLs"
    );
    let span_calls = outcome
        .genotyped_calls
        .iter()
        .filter(|c| {
            c.event.start_1based.get() >= REGION_LO && c.event.start_1based.get() <= REGION_HI
        })
        .count();
    assert!(
        span_calls >= 1,
        "Phase C: expected ≥1 genotyped call in mid-A span, got {span_calls}"
    );
}
