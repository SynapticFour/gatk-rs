//! R4-2 / L7: dense GIAB insertion `20:10001436 A→AAGGCT` must genotype.

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::genome_loc::GenomePosition;
use gatk_haplotypecaller::read_model::ReadFilterParams;
use gatk_haplotypecaller::walker_traversal::{
    flatten_assembly_regions, traverse_assembly_region_walker,
};
use gatk_haplotypecaller::{
    call_disposition, AssemblyRegionCallDisposition, CallRegionArgs, HaplotypeCallerEngine,
    WalkerTraversalConfig,
};
use std::path::Path;

#[test]
fn dense_giab_insertion_reaches_event_map_and_call() {
    let repo = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
    let ref_fa = repo.join("parity/realworld/assets/hs37d5.simple.fa");
    let bam = repo.join("parity/realworld/na12878_giab_window_b37/NA12878_giab_window.b37.bam");
    // Realworld assets are gitignored / fetched locally — skip on bare CI checkouts.
    if !ref_fa.is_file() || !bam.is_file() {
        return;
    }
    let dict = SequenceDictionary::from_fasta_path(&ref_fa).expect("dict");
    let specs = parse_intervals_cli_string(&dict, "20:10001400-10001500").expect("interval");
    let filters = ReadFilterParams::gatk_standard_hc();
    let cfg = WalkerTraversalConfig::gatk_haplotype_caller_production(100);
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
    let events = outcome.assembly.variation_events();
    let event = events
        .iter()
        .find(|e| {
            e.start_1based == GenomePosition::new_1based(10001436)
                && e.ref_allele == "A"
                && e.alt_allele == "AAGGCT"
        })
        .expect("EventMap must contain A>AAGGCT");
    let (rr, ra) = gatk_haplotypecaller::read_event_discovery::read_allele_depths_at_locus(
        &region.reads,
        event,
        region.start.get(),
    );
    assert!(ra >= 2, "pileup alt AD expected >=2, got {rr},{ra}");
    assert!(
        outcome.genotyped_calls.iter().any(|c| {
            c.event.start_1based == GenomePosition::new_1based(10001436)
                && c.event.ref_allele == "A"
                && c.event.alt_allele == "AAGGCT"
        }),
        "expected genotyped call A>AAGGCT at 10001436"
    );
}
