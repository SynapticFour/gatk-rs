//! mid-B region 923181 strict genotyping probe.
//! Run: `P12_REFERENCE=… cargo test -p gatk-haplotypecaller p12_region_923181_genotyping --release -- --nocapture`

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

const REGION_LO: u64 = 92318129;
const REGION_HI: u64 = 92318286;

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
fn p12_region_923181_genotyping_probe() {
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
        .expect("active region covering 92318129-92318286");
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
    assert!(
        outcome.assembly.variation_events().len() >= 2,
        "Phase C: need EventMap events"
    );
    assert!(
        !outcome.genotyped_calls.is_empty(),
        "Phase C: expected genotyped_calls"
    );
    let span_calls = outcome
        .genotyped_calls
        .iter()
        .filter(|c| {
            c.event.start_1based.get() >= REGION_LO && c.event.start_1based.get() <= REGION_HI
        })
        .count();
    assert!(
        span_calls >= 3,
        "Phase C mid-B: expected ≥3 genotyped calls in span, got {span_calls}"
    );
    for &(pos, ref_a, alt_a) in &[(92318210, "A", "G"), (92318227, "C", "G")] {
        let event = VariationEvent {
            contig: "2".into(),
            start_1based: GenomePosition::new_1based(pos),
            end_1based: GenomePosition::new_1based(pos),
            ref_allele: ref_a.to_string(),
            alt_allele: alt_a.to_string(),
        };
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
        let on_map = outcome.assembly.variation_events().iter().any(|e| {
            e.start_1based == GenomePosition::new_1based(pos)
                && e.ref_allele == ref_a
                && e.alt_allele == alt_a
        });
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
        .expect("diagnose");
        let ok = outcome.genotyped_calls.iter().any(|c| {
            c.event.start_1based == GenomePosition::new_1based(pos)
                && c.event.ref_allele == ref_a
                && c.event.alt_allele == alt_a
        });
        eprintln!(
            "mid_b_site\t{pos}\t{ref_a}/{alt_a}\tevent={on_map}\temit_AD={rr}/{ra}\tgenotyped={ok}\tdiag={diag:?}"
        );
        assert!(ok, "{pos} {ref_a}/{alt_a} must be genotyped");
    }
}
