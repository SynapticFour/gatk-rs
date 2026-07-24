//! Mid-B `call_none` sites (92317399–92317412, 92318593): `call_region` must return `Some`.
//! Run: `P12_REFERENCE=…/hs37d5.simple.fa cargo test -p gatk-haplotypecaller p12_call_none_mid_b --release`

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::genome_loc::GenomePosition;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    AssemblyRegionCallDisposition, CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams,
    WalkerTraversalConfig,
};
use std::path::Path;

const SITES: &[(u64, &str, &str)] = &[
    (92317399, "C", "A"),
    (92317407, "T", "C"),
    (92317412, "G", "C"),
    (92318593, "T", "A"),
];

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
fn p12_call_none_mid_b() {
    let Some((ref_fasta, bam)) = p12_paths() else {
        eprintln!("skip: set P12_REFERENCE");
        return;
    };
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, "2:92317000-92319000").expect("interval");
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

    for &(pos, ref_a, alt_a) in SITES {
        let region = regions
            .iter()
            .find(|r| {
                matches!(
                    call_disposition(r),
                    AssemblyRegionCallDisposition::ActiveFull
                ) && r.start.get() <= pos
                    && r.end.get() >= pos
            })
            .unwrap_or_else(|| panic!("no ActiveFull region for {pos}"));
        let outcome =
            HaplotypeCallerEngine::call_region(region, &dict, &ref_fasta, &args).expect("call");
        let Some(outcome) = outcome else {
            panic!(
                "call_region None for {pos} in region {}-{}",
                region.start.get(),
                region.end.get()
            );
        };
        let has_event = outcome.assembly.variation_events().iter().any(|e| {
            e.start_1based == GenomePosition::new_1based(pos)
                && e.ref_allele == ref_a
                && e.alt_allele == alt_a
        });
        assert!(
            has_event,
            "{pos} {ref_a}/{alt_a} missing from variation_events"
        );
    }
}
