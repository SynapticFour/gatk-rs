//! Regression gate for two formerly divergent P12 sites (92316328, 92325205).
//! Run: `P12_PHASE_E=1 P12_REFERENCE=… cargo test -p gatk-haplotypecaller p12_two_remaining_sites --release -- --ignored --nocapture`

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::genome_loc::GenomePosition;
use gatk_haplotypecaller::hc_emit_policy::explain_strict_java_emit_gates;
use gatk_haplotypecaller::hc_genotyping_engine::java_emit_af_decision;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, read_event_discovery::read_allele_depths_at_locus,
    traverse_assembly_region_walker, try_emit_call_region_variants, AssemblyRegionCallDisposition,
    CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
};
use std::path::Path;

const SITES: &[(u64, &str, &str, &str)] = &[
    (92316328, "T", "A", "2:92316227-92316475"),
    (92325205, "G", "A", "2:92325071-92325332"),
];

#[test]
#[ignore = "Phase E: two-site emit gate (~5 min each)"]
fn p12_two_remaining_sites() {
    if std::env::var("P12_PHASE_E").is_err() {
        eprintln!("skip: set P12_PHASE_E=1");
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
    if !ref_path.is_file() || !bam.is_file() {
        eprintln!("skip: missing P12 assets");
        return;
    }
    let dict = SequenceDictionary::from_fasta_path(&ref_path).expect("dict");
    let args = CallRegionArgs::strict_java();

    for &(pos, ref_a, alt_a, interval) in SITES {
        let specs = parse_intervals_cli_string(&dict, interval).expect("interval");
        let walk = traverse_assembly_region_walker(
            &dict,
            &specs,
            &ref_path,
            &bam,
            &ReadFilterParams::gatk_standard_hc(),
            &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
        )
        .expect("walk");
        let regions = flatten_assembly_regions(&walk);
        let mut emitted = false;
        for region in &regions {
            if !matches!(
                call_disposition(region),
                AssemblyRegionCallDisposition::ActiveFull
            ) {
                continue;
            }
            let Some(outcome) =
                HaplotypeCallerEngine::call_region(region, &dict, &ref_path, &args).expect("call")
            else {
                continue;
            };
            let Some(call) = outcome.genotyped_calls.iter().find(|c| {
                c.event.start_1based == GenomePosition::new_1based(pos)
                    && c.event.ref_allele == ref_a
                    && c.event.alt_allele == alt_a
            }) else {
                continue;
            };
            let pad = outcome
                .assembly
                .haplotypes
                .iter()
                .find(|h| h.is_reference)
                .and_then(|h| h.genome_loc.map(|g| g.start_1based()))
                .unwrap_or_else(|| outcome.assembly.padded_reference_start_1based());
            let (read_ref, read_alt) = read_allele_depths_at_locus(&region.reads, &call.event, pad);
            let gl = &call.genotype.genotype_log10_likelihoods;
            let af = java_emit_af_decision(gl, 10.0).expect("af");
            let gates = explain_strict_java_emit_gates(
                &call.event,
                gl,
                &call.genotype.format,
                10.0,
                true,
                read_ref,
                read_alt,
                &[],
            )
            .expect("gates");
            if std::env::var("GATK_RS_P12_CLUSTER_DEBUG")
                .ok()
                .is_some_and(|v| matches!(v.as_str(), "1" | "true" | "TRUE"))
            {
                eprintln!(
                    "SITE {pos}\t{ref_a}/{alt_a}\tread_AD={read_ref}/{read_alt}\t{af:?}\t{gates}"
                );
            }
            for rec in
                try_emit_call_region_variants(region, &outcome, "SAMPLE", 10.0).expect("emit")
            {
                if rec.position == pos
                    && rec.reference == ref_a
                    && rec.alternate.first().map(|s| s.as_str()) == Some(alt_a)
                {
                    emitted = true;
                }
            }
        }
        assert!(emitted, "{pos} {ref_a}/{alt_a} must emit on {interval}");
    }
}
