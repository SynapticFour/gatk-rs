//! L4: PairHMM trace for P12 site 92318244 (T/C sparse hom-alt, Java PL 90,6,0 AD 0,2).
//! Run: `P12_PHASE_E=1 P12_REFERENCE=parity/realworld/assets/hs37d5.simple.fa cargo test -p gatk-haplotypecaller p12_site_92318244 --release -- --ignored --nocapture`

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::VariationEvent;
use gatk_haplotypecaller::genome_loc::GenomePosition;
use gatk_haplotypecaller::hc_genotyping_engine::read_allele_depths_for_strict_emit;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, pairhmm_locus_trace_dump,
    traverse_assembly_region_walker, try_emit_call_region_variants, AssemblyRegionCallDisposition,
    CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
};
use std::path::Path;

const POS: u64 = 92318244;
const REGION_LO: u64 = 92318129;
const REGION_HI: u64 = 92318289;

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
#[ignore = "P12 BAM: 92318244 PairHMM trace"]
fn p12_site_92318244_pairhmm_trace() {
    if std::env::var("P12_PHASE_E").is_err() {
        eprintln!("skip: set P12_PHASE_E=1");
        return;
    }
    let Some((ref_fasta, bam)) = p12_paths() else {
        eprintln!("skip: P12_REFERENCE / BAM");
        return;
    };

    eprintln!("=== Java target ===\nPOS={POS} PL=90,6,0 GQ=6 AD=0,2 DP=2");

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
        .expect("active region");
    let args = CallRegionArgs::strict_java();
    let gt_cfg = args.genotyping.clone();
    let outcome = HaplotypeCallerEngine::call_region(region, &dict, &ref_fasta, &args)
        .expect("call")
        .expect("outcome");
    let ref_hap = outcome
        .assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .expect("ref hap");
    let pad = ref_hap
        .genome_loc
        .map(|g| g.start_1based())
        .unwrap_or_else(|| outcome.assembly.padded_reference_start_1based());
    let full_ref = outcome.assembly.reference_bases();
    let full_pad = outcome.assembly.padded_reference_start_1based();
    let event = VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(POS),
        end_1based: GenomePosition::new_1based(POS),
        ref_allele: "T".into(),
        alt_allele: "C".into(),
    };
    let (geno_rr, geno_ra) = read_allele_depths_for_strict_emit(
        &outcome.genotyping_reads,
        Some(region.reads.as_slice()),
        &event,
        pad,
        &gt_cfg,
        &ref_hap.bases,
        full_ref,
        full_pad,
    );
    eprintln!("pileup_AD genotyping={geno_rr}/{geno_ra}");
    let dump = pairhmm_locus_trace_dump(
        &event,
        &outcome.read_likelihoods,
        &outcome.genotyping_reads,
        &outcome.assembly.haplotypes,
        &ref_hap.bases,
        pad,
        region.start.get(),
        region.end.get(),
        outcome.assembly.max_mnp_distance(),
        &gt_cfg,
    )
    .expect("trace");
    eprintln!("\n=== PairHMM trace ===\n{dump}");
    for c in &outcome.genotyped_calls {
        if c.event.start_1based == GenomePosition::new_1based(POS) {
            eprintln!(
                "genotyped GL={:.4?} PL={:?} AD={:?} DP={}",
                c.genotype.genotype_log10_likelihoods,
                c.genotype.format.pl,
                c.genotype.format.ad,
                c.genotype.format.dp
            );
        }
    }
    let emitted = try_emit_call_region_variants(region, &outcome, "SAMPLE", 10.0).expect("emit");
    if let Some(rec) = emitted.iter().find(|r| r.position == POS) {
        let s = rec.samples.first().expect("sample");
        eprintln!(
            "VCF PL={:?} AD={:?} DP={:?} GQ={:?}",
            s.pl, s.ad, s.dp, s.gq
        );
    }
}
