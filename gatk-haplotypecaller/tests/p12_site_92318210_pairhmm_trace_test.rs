//! L4: PairHMM trace for P12 site 92318210 (A/G het-trap → hom-alt, Java PL 45,3,0 AD 0,1).

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::VariationEvent;
use gatk_haplotypecaller::genome_loc::GenomePosition;
use gatk_haplotypecaller::hc_genotyping_engine::read_allele_depths_for_strict_emit;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
};
use std::path::Path;

const POS: u64 = 92318210;

#[test]
#[ignore = "P12 BAM: 92318210 trace"]
fn p12_site_92318210_pairhmm_trace() {
    if std::env::var("P12_PHASE_E").is_err() {
        return;
    }
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
    let region = flatten_assembly_regions(&walk)
        .into_iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= 92318129
                && r.end.get() >= 92318289
        })
        .expect("region");
    let args = CallRegionArgs::strict_java();
    let outcome = HaplotypeCallerEngine::call_region(&region, &dict, &ref_path, &args)
        .expect("call")
        .expect("outcome");
    let ref_hap = outcome
        .assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .unwrap();
    let pad = ref_hap.genome_loc.map(|g| g.start_1based()).unwrap_or(1);
    let event = VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(POS),
        end_1based: GenomePosition::new_1based(POS),
        ref_allele: "A".into(),
        alt_allele: "G".into(),
    };
    let (rr, ra) = read_allele_depths_for_strict_emit(
        &outcome.genotyping_reads,
        Some(region.reads.as_slice()),
        &event,
        pad,
        &args.genotyping,
        &ref_hap.bases,
        outcome.assembly.reference_bases(),
        outcome.assembly.padded_reference_start_1based(),
    );
    eprintln!("pileup_AD={rr}/{ra}");
    for c in &outcome.genotyped_calls {
        if c.event.start_1based == GenomePosition::new_1based(POS) {
            eprintln!(
                "genotyped PL={:?} AD={:?} GL={:.4?}",
                c.genotype.format.pl, c.genotype.format.ad, c.genotype.genotype_log10_likelihoods
            );
        }
    }
    let emitted = try_emit_call_region_variants(&region, &outcome, "SAMPLE", 10.0).expect("emit");
    if let Some(rec) = emitted.iter().find(|r| r.position == POS) {
        let s = rec.samples.first().unwrap();
        eprintln!("VCF PL={:?} AD={:?} DP={:?}", s.pl, s.ad, s.dp);
    }
}
