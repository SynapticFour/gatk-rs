//! P12 NA12878: why full-interval HC emits 0 VCF rows while Java emits ~66.
//! Run: `P12_REFERENCE=…/hs37d5.simple.fa cargo test -p gatk-haplotypecaller p12_emit_audit -- --nocapture`

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
};

const STAND_EMIT: f64 = 10.0;
use std::path::Path;

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
fn p12_emit_audit() {
    let Some((ref_fasta, bam)) = p12_paths() else {
        eprintln!("skip p12_emit_audit: set P12_REFERENCE to hs37d5.simple.fa");
        return;
    };
    let interval = "2:92300000-92350000";
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, interval).expect("interval");
    let filters = ReadFilterParams::gatk_standard_hc();
    let cfg = WalkerTraversalConfig::gatk_haplotype_caller_production(100);
    let walk = traverse_assembly_region_walker(&dict, &specs, &ref_fasta, &bam, &filters, &cfg)
        .expect("walk");
    let regions = flatten_assembly_regions(&walk);
    let args = CallRegionArgs::default();

    let mut active = 0usize;
    let mut call_none = 0usize;
    let mut no_variation = 0usize;
    let mut no_genotype = 0usize;
    let mut gq_filtered = 0usize;
    let mut ref_only_haps = 0usize;
    let mut emitted = 0usize;
    let mut genotyped_calls_total = 0usize;

    for region in &regions {
        if !matches!(
            call_disposition(region),
            AssemblyRegionCallDisposition::ActiveFull
        ) {
            continue;
        }
        active += 1;
        let Some(outcome) =
            HaplotypeCallerEngine::call_region(region, &dict, &ref_fasta, &args).expect("call")
        else {
            call_none += 1;
            continue;
        };
        if !outcome.assembly.is_variation_present() {
            no_variation += 1;
        }
        let non_ref = outcome
            .assembly
            .haplotypes
            .iter()
            .filter(|h| !h.is_reference)
            .count();
        if non_ref == 0 {
            ref_only_haps += 1;
        }
        genotyped_calls_total += outcome.genotyped_calls.len();
        if outcome.genotype.is_none() && outcome.genotyped_calls.is_empty() {
            no_genotype += 1;
        }
        let recs =
            try_emit_call_region_variants(region, &outcome, "SAMPLE", STAND_EMIT).expect("emit");
        if recs.is_empty() {
            if let Some(gt) = &outcome.genotype {
                if (gt.format.gq.as_i32() as f64) < STAND_EMIT {
                    gq_filtered += 1;
                }
            }
        } else {
            emitted += recs.len();
            eprintln!(
                "EMIT\t{}:{}-{}\t{}",
                region.contig,
                region.start.get(),
                region.end.get(),
                recs.len()
            );
        }
    }

    eprintln!("=== P12 emit audit ({interval}) ===");
    eprintln!("active_regions\t{active}");
    eprintln!("call_region_none\t{call_none}");
    eprintln!("call_some_no_variation\t{no_variation}");
    eprintln!("ref_only_haplotypes\t{ref_only_haps}");
    eprintln!("no_genotype\t{no_genotype}");
    eprintln!("genotyped_calls_total\t{genotyped_calls_total}");
    eprintln!("gq_filtered_empty_emit\t{gq_filtered}");
    eprintln!("vcf_rows_emitted\t{emitted}");
    assert!(active > 0, "expected active regions on P12 interval");
}
