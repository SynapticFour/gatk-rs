//! Graph-only ASM-8: five historical `no_event` sites must appear on the EventMap.
//! Run: `P12_REFERENCE=… cargo test -p gatk-haplotypecaller p12_graph_only_gap_events --release -- --nocapture`

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::genome_loc::GenomePosition;
use gatk_haplotypecaller::read_event_discovery::P12_PHASE_E_GAP_SNPS;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    AssemblyRegionCallDisposition, CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams,
    WalkerTraversalConfig,
};
use std::path::Path;

const MID_A_LO: u64 = 92316227;
const MID_A_HI: u64 = 92316475;
const TAIL_LO: u64 = 92325071;
const TAIL_HI: u64 = 92325332;

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

fn has_event(
    events: &[gatk_haplotypecaller::event_map::VariationEvent],
    pos: u64,
    ref_a: &str,
    alt_a: &str,
) -> bool {
    events.iter().any(|e| {
        e.start_1based == GenomePosition::new_1based(pos)
            && e.ref_allele == ref_a
            && e.alt_allele == alt_a
    })
}

#[test]
fn p12_graph_only_gap_events() {
    std::env::remove_var("GATK_RS_P12_EVENT_REGISTRY");
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
    let args = CallRegionArgs::strict_java();

    let mid_a = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= MID_A_LO
                && r.end.get() >= MID_A_HI
        })
        .expect("mid-A active region");
    let tail = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= TAIL_LO
                && r.end.get() >= TAIL_HI
        })
        .expect("tail active region");

    for (label, region) in [("mid_a", mid_a), ("tail", tail)] {
        let outcome = HaplotypeCallerEngine::call_region(region, &dict, &ref_fasta, &args)
            .expect("call")
            .expect("call_region outcome");
        let events = outcome.assembly.variation_events();
        for &(pos, ref_a, alt_a) in P12_PHASE_E_GAP_SNPS {
            if pos < region.start.get() || pos > region.end.get() {
                continue;
            }
            let ok = has_event(events, pos, ref_a, alt_a);
            eprintln!("graph_only_gap\t{label}\t{pos}\t{ref_a}/{alt_a}\tevent={ok}");
            assert!(
                ok,
                "graph-only gap site missing from EventMap: {pos} {ref_a}/{alt_a} ({label})"
            );
        }
    }
}
