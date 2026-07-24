//! Probe P12 region 92316227–92316475 (dense java-only `no_event` bucket).
//! Under `strict_java`, `call_region` may return `None` until ASM/EventMap matches Java (CR-2).
//! Run: `P12_REFERENCE=…/hs37d5.simple.fa cargo test -p gatk-haplotypecaller p12_region_923162_probe --release -- --nocapture`

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
fn p12_region_923162_probe() {
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
    eprintln!("region\t{}-{}", region.start.get(), region.end.get());
    let args = CallRegionArgs::strict_java();
    let outcome =
        HaplotypeCallerEngine::call_region(region, &dict, &ref_fasta, &args).expect("call");
    let outcome = outcome.expect("strict_java: call_region must run when reads support variation");
    eprintln!(
        "haps={} variation_present={} events={} genotyped_calls={} read_ll={}",
        outcome.assembly.haplotypes.len(),
        outcome.assembly.is_variation_present(),
        outcome.assembly.variation_events().len(),
        outcome.genotyped_calls.len(),
        outcome.read_likelihoods.len()
    );
    for (i, h) in outcome.assembly.haplotypes.iter().enumerate() {
        eprintln!(
            "hap{i}\tref={}\tlen={}\talign={}\tcigar={}",
            h.is_reference,
            h.bases.len(),
            h.alignment_start_hap_wrt_ref,
            h.cigar
                .as_ref()
                .map(|c| c.to_gatk_string())
                .unwrap_or_default()
        );
    }
    let mut in_span = 0usize;
    for e in outcome.assembly.variation_events() {
        if e.start_1based.get() >= REGION_LO && e.start_1based.get() <= REGION_HI {
            in_span += 1;
            eprintln!(
                "event\t{} {}/{}",
                e.start_1based.get(),
                e.ref_allele,
                e.alt_allele
            );
        }
    }
    eprintln!("events_in_span\t{in_span}");
    assert!(
        in_span >= 2,
        "strict_java: expected EventMap events in span when call_region runs; got {in_span}"
    );
    assert!(
        outcome.assembly.haplotypes.iter().any(|h| {
            !h.is_reference
                && h.cigar
                    .as_ref()
                    .is_some_and(|c| c.elements.iter().any(|e| e.operator.is_indel()))
        }),
        "ASM-8: expect ≥1 alt hap with I/D CIGAR in mid-region strict assembly"
    );
}
