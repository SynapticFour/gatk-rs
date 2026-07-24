//! Dump emit gates for the six P12 graph-only sites that were `genotyped_not_emitted`.
//! Run: `P12_PHASE_E=1 P12_REFERENCE=… cargo test -p gatk-haplotypecaller --features parity_harness p12_emit_six_java_sites --release -- --ignored --nocapture`
#![cfg(feature = "parity_harness")]

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::genome_loc::GenomePosition;
use gatk_haplotypecaller::hc_emit_policy::{
    explain_strict_java_emit_gates, passes_strict_java_emit_for_genotyped_call,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, read_event_discovery::read_allele_depths_at_locus,
    traverse_assembly_region_walker, try_emit_call_region_variants, AssemblyRegionCallDisposition,
    CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
};
use std::path::Path;

const SIX: &[(u64, &str, &str)] = &[
    (92316296, "A", "T"),
    (92316315, "C", "G"),
    (92316328, "T", "A"),
    (92324471, "C", "T"),
    (92325193, "C", "T"),
    (92325205, "G", "A"),
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
    if ref_path.is_file() && bam.is_file() {
        Some((ref_path, bam))
    } else {
        None
    }
}

#[test]
#[ignore = "Phase E: six-site emit gate dump (~10+ min)"]
fn p12_emit_six_java_sites() {
    if std::env::var("P12_PHASE_E").is_err() {
        eprintln!("skip: set P12_PHASE_E=1");
        return;
    }
    let Some((ref_fasta, bam)) = p12_paths() else {
        eprintln!("skip: P12_REFERENCE / BAM");
        return;
    };
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, "2:92300000-92350000").expect("interval");
    let filters = ReadFilterParams::gatk_standard_hc();
    let cfg = WalkerTraversalConfig::gatk_haplotype_caller_production(100);
    let walk = traverse_assembly_region_walker(&dict, &specs, &ref_fasta, &bam, &filters, &cfg)
        .expect("walk");
    let regions = flatten_assembly_regions(&walk);
    let args = CallRegionArgs::strict_java();
    let stand = args.genotyping.stand_emit_confidence;

    let mut found = 0usize;
    for region in &regions {
        if !matches!(
            call_disposition(region),
            AssemblyRegionCallDisposition::ActiveFull
        ) {
            continue;
        }
        let Some(outcome) =
            HaplotypeCallerEngine::call_region(region, &dict, &ref_fasta, &args).expect("call")
        else {
            continue;
        };
        let pad = outcome
            .assembly
            .haplotypes
            .iter()
            .find(|h| h.is_reference)
            .and_then(|h| h.genome_loc.map(|g| g.start_1based()))
            .unwrap_or_else(|| outcome.assembly.padded_reference_start_1based());
        for (pos, ref_a, alt_a) in SIX {
            let Some(call) = outcome.genotyped_calls.iter().find(|c| {
                c.event.start_1based == GenomePosition::new_1based(*pos)
                    && c.event.ref_allele == *ref_a
                    && c.event.alt_allele == *alt_a
            }) else {
                continue;
            };
            found += 1;
            let (read_ref, read_alt) = read_allele_depths_at_locus(&region.reads, &call.event, pad);
            let gates = explain_strict_java_emit_gates(
                &call.event,
                &call.genotype.genotype_log10_likelihoods,
                &call.genotype.format,
                stand,
                args.genotyping.genotype_stored_events_only,
                read_ref,
                read_alt,
                &[],
            )
            .expect("gates");
            let strict = passes_strict_java_emit_for_genotyped_call(
                &call.event,
                &call.genotype.genotype_log10_likelihoods,
                &call.genotype.format,
                stand,
                args.genotyping.genotype_stored_events_only,
                read_ref,
                read_alt,
                false,
                &[],
            )
            .expect("strict");
            let recs =
                try_emit_call_region_variants(region, &outcome, "SAMPLE", stand).expect("emit");
            let emitted = recs.iter().any(|r| r.position == *pos);
            eprintln!("SITE\t{pos}\t{ref_a}/{alt_a}\tstrict={strict}\temitted={emitted}\t{gates}");
        }
    }
    eprintln!("found_genotyped\t{found}/{}", SIX.len());
    assert_eq!(found, SIX.len(), "expected all six genotyped on P12");
}
