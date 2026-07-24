//! mid-A no_event sites 92316416–92316458 must appear in EventMap + genotyped_calls.
//! Run: `P12_REFERENCE=… cargo test -p gatk-haplotypecaller p12_region_923162_mid_a --release -- --nocapture`

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::genome_loc::GenomePosition;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
};
use std::path::Path;

const REGION_LO: u64 = 92316227;
const REGION_HI: u64 = 92316475;

use gatk_haplotypecaller::read_event_discovery::P12_MID_A_JAVA_SNPS;

const MID_A_SITES: &[(u64, &str, &str)] = P12_MID_A_JAVA_SNPS;
/// Gap-backfill anchor in same active region (hom-alt shaped GL regression).
const MID_A_GAP_ANCHORS: &[(u64, &str, &str)] = &[(92316315, "C", "G")];

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
fn p12_region_923162_mid_a_genotyping() {
    // ASM-8 graph-only: no list inject.
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
    let events = outcome.assembly.variation_events();
    for &(pos, ref_a, alt_a) in MID_A_SITES.iter().chain(MID_A_GAP_ANCHORS.iter()) {
        let in_events = has_event(events, pos, ref_a, alt_a);
        let in_calls = outcome.genotyped_calls.iter().any(|c| {
            c.event.start_1based == GenomePosition::new_1based(pos)
                && c.event.ref_allele == ref_a
                && c.event.alt_allele == alt_a
        });
        eprintln!("mid_a\t{pos}\t{ref_a}/{alt_a}\tevent={in_events}\tcall={in_calls}");
        assert!(
            in_events,
            "Phase C mid-A: missing EventMap event {pos} {ref_a}/{alt_a}"
        );
        assert!(
            in_calls,
            "Phase C mid-A: missing genotyped call {pos} {ref_a}/{alt_a}"
        );
        let emitted = try_emit_call_region_variants(region, &outcome, "SAMPLE", 10.0)
            .expect("emit")
            .iter()
            .any(|r| {
                r.position == pos
                    && r.reference == ref_a
                    && r.alternate.first().map(String::as_str) == Some(alt_a)
            });
        assert!(
            emitted,
            "Phase C mid-A: missing VCF emit {pos} {ref_a}/{alt_a}"
        );
    }
}
