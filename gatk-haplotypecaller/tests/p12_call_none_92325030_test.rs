//! N1: `92325030` must not return `call_region` = None when cluster/read variation exists.
//! Run: `P12_REFERENCE=…/hs37d5.simple.fa cargo test -p gatk-haplotypecaller p12_call_none_92325030 --release`

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::genome_loc::GenomePosition;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
};
use std::path::Path;

const TARGET: u64 = 92325030;

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
fn p12_call_none_92325030() {
    let Some((ref_fasta, bam)) = p12_paths() else {
        eprintln!("skip: set P12_REFERENCE");
        return;
    };
    let interval = format!("2:{TARGET}-{TARGET}");
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, &interval).expect("interval");
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
            ) && r.start.get() <= TARGET
                && r.end.get() >= TARGET
        })
        .expect("active region covering 92325030");
    let outcome =
        HaplotypeCallerEngine::call_region(region, &dict, &ref_fasta, &CallRegionArgs::default())
            .expect("call");
    let outcome = outcome.expect("call_region outcome");
    let has_at = outcome.assembly.variation_events().iter().any(|e| {
        e.start_1based == GenomePosition::new_1based(TARGET)
            && e.ref_allele == "A"
            && e.alt_allele == "T"
    });
    assert!(has_at, "92325030 A/T must appear in variation_events");
    let genotyped = outcome.genotyped_calls.iter().any(|c| {
        c.event.start_1based == GenomePosition::new_1based(TARGET)
            && c.event.ref_allele == "A"
            && c.event.alt_allele == "T"
    });
    assert!(
        genotyped,
        "92325030 A/T must be genotyped; calls={:?}",
        outcome
            .genotyped_calls
            .iter()
            .map(|c| (
                c.event.start_1based,
                &c.event.ref_allele,
                &c.event.alt_allele
            ))
            .collect::<Vec<_>>()
    );
    let emitted = try_emit_call_region_variants(region, &outcome, "SAMPLE", 10.0).expect("emit");
    assert!(
        emitted.iter().any(|r| {
            r.position == TARGET
                && r.reference == "A"
                && r.alternate.first().map(String::as_str) == Some("T")
        }),
        "92325030 A/T must emit VCF row; rows={}",
        emitted.len()
    );
}
