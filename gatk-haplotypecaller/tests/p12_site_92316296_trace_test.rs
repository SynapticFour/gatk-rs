//! Trace FORMAT tier path for P12 gap hom-alt pileup site 92316296.

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::VariationEvent;
use gatk_haplotypecaller::genome_loc::GenomePosition;
use gatk_haplotypecaller::hc_genotyping_engine::{
    diagnose_genotype_variation_event, read_allele_depths_for_strict_emit,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    AssemblyRegionCallDisposition, CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams,
    WalkerTraversalConfig,
};
use std::path::Path;

const POS: u64 = 92316296;

#[test]
#[ignore = "P12 BAM: 92316296 trace"]
fn p12_site_92316296_trace() {
    if std::env::var("P12_PHASE_E").is_err() {
        return;
    }
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
    let ref_path = root.join("parity/realworld/assets/hs37d5.simple.fa");
    let bam = root.join("parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam");
    let dict = SequenceDictionary::from_fasta_path(&ref_path).expect("dict");
    let specs = parse_intervals_cli_string(&dict, "2:92300000-92350000").expect("interval");
    let walk = traverse_assembly_region_walker(
        &dict,
        &specs,
        &ref_path,
        &bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
    )
    .expect("walk");
    let args = CallRegionArgs::strict_java();
    let event = VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(POS),
        end_1based: GenomePosition::new_1based(POS),
        ref_allele: "A".into(),
        alt_allele: "T".into(),
    };
    for region in flatten_assembly_regions(&walk) {
        if !matches!(
            call_disposition(&region),
            AssemblyRegionCallDisposition::ActiveFull
        ) {
            continue;
        }
        if region.start.get() > POS || region.end.get() < POS {
            continue;
        }
        let Some(outcome) =
            HaplotypeCallerEngine::call_region(&region, &dict, &ref_path, &args).expect("call")
        else {
            continue;
        };
        let ref_hap = outcome
            .assembly
            .haplotypes
            .iter()
            .find(|h| h.is_reference)
            .expect("ref");
        let pad = ref_hap
            .genome_loc
            .map(|g| g.start_1based())
            .unwrap_or_else(|| outcome.assembly.padded_reference_start_1based());
        let full_ref = outcome.assembly.reference_bases();
        let full_pad = outcome.assembly.padded_reference_start_1based();
        let (rr, ra) = read_allele_depths_for_strict_emit(
            &outcome.genotyping_reads,
            Some(region.reads.as_slice()),
            &event,
            pad,
            &args.genotyping,
            &ref_hap.bases,
            full_ref,
            full_pad,
        );
        eprintln!("emit_AD={rr}/{ra}");
        let diag = diagnose_genotype_variation_event(
            &event,
            &outcome.read_likelihoods,
            &outcome.genotyping_reads,
            &outcome.genotyping_reads,
            Some(region.reads.as_slice()),
            &outcome.assembly.haplotypes,
            &ref_hap.bases,
            pad,
            full_ref,
            full_pad,
            region.start.get(),
            region.end.get(),
            outcome.assembly.max_mnp_distance(),
            &args.genotyping,
        )
        .expect("diag");
        if let Ok(call) = diag {
            eprintln!(
                "PL={:?} AD={:?} GL={:.4?}",
                call.genotype.format.pl,
                call.genotype.format.ad,
                call.genotype.genotype_log10_likelihoods
            );
        } else {
            eprintln!("reject={diag:?}");
        }
        if let Some(c) = outcome
            .genotyped_calls
            .iter()
            .find(|c| c.event.start_1based == GenomePosition::new_1based(POS))
        {
            eprintln!(
                "region_call PL={:?} AD={:?}",
                c.genotype.format.pl, c.genotype.format.ad
            );
        }
    }
}
