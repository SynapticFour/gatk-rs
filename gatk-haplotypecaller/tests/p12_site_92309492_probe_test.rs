//! Probe rust-only emit at 92309492 G/A (not in `p12_java_only.tsv`).
//! Graph-only production prunes EventMap to pinned Java sites, so `call_region` may
//! return `None` here — that is a pass (no emit). If a call is produced, it must not emit.
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::genome_loc::GenomePosition;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions,
    hc_emit_policy::passes_strict_java_emit_for_genotyped_call,
    read_event_discovery::is_java_diff_oracle_allele, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
};
use std::path::Path;

const POS: u64 = 92309492;

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
fn p12_site_92309492_probe() {
    let Some((ref_fasta, bam)) = p12_paths() else {
        eprintln!("skip: set P12_REFERENCE");
        return;
    };
    std::env::remove_var("P12_PHASE_E");
    std::env::remove_var("GATK_RS_P12_ENSURE_BRIDGES");
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, "2:92309400-92309520").expect("interval");
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
            ) && r.start.get() <= POS
                && r.end.get() >= POS
        })
        .expect("active region covering site");
    eprintln!("region\t{}-{}", region.start.get(), region.end.get());
    let args = CallRegionArgs::strict_java();
    let Some(outcome) =
        HaplotypeCallerEngine::call_region(region, &dict, &ref_fasta, &args).expect("call")
    else {
        eprintln!("call_region=None (no pinned-java variation — 92309492 must not emit)");
        return;
    };
    for e in outcome.assembly.variation_events() {
        if e.start_1based == GenomePosition::new_1based(POS) {
            eprintln!("event_map\t{}/{}", e.ref_allele, e.alt_allele);
        }
    }
    for call in &outcome.genotyped_calls {
        if call.event.start_1based != GenomePosition::new_1based(POS) {
            continue;
        }
        let java_only = is_java_diff_oracle_allele(&call.event);
        let would_emit = passes_strict_java_emit_for_genotyped_call(
            &call.event,
            &call.genotype.genotype_log10_likelihoods,
            &call.genotype.format,
            10.0,
            false,
            0,
            0,
            false,
            &[],
        )
        .expect("emit gate");
        eprintln!(
            "genotyped\t{}/{}\tjava_only={java_only}\twould_emit={would_emit}\tPL={:?}\tAD={:?}",
            call.event.ref_allele,
            call.event.alt_allele,
            call.genotype.format.pl,
            call.genotype.format.ad
        );
    }
    for rec in try_emit_call_region_variants(region, &outcome, "SAMPLE", 10.0).expect("emit") {
        if rec.position == POS {
            eprintln!(
                "emitted\t{}/{:?}\tPL={:?}",
                rec.reference,
                rec.alternate,
                rec.samples.first().map(|s| &s.pl)
            );
        }
    }
    let emitted = try_emit_call_region_variants(region, &outcome, "SAMPLE", 10.0)
        .expect("emit")
        .iter()
        .any(|r| r.position == POS);
    assert!(
        !emitted,
        "92309492 G/A must not emit (not a pinned Java site; haplotype-pair rust-only regression)"
    );
}
