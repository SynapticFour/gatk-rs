//! 6R.87 holdout: Java `target.overlaps` coordinates at the canonical T/C site.
//!
//! Skipped unless `HOLDOUT_6R87=1`. Coordinate-free contract lives in
//! `forensic_6r87_variant_overlap_coordinate_contract`.
//!
//! ```text
//! HOLDOUT_6R87=1 cargo test -p gatk-haplotypecaller --test holdout_6r87_variant_overlap -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, take_colocated_merge_numerics,
    traverse_assembly_region_walker, try_emit_call_region_variants, AssemblyRegionCallDisposition,
    CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
    DEFAULT_STAND_EMIT_CONFIDENCE,
};
use serde_json::json;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const POS_SNP: u64 = 29_456_344;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

#[test]
fn holdout_6r87_variant_overlap_29456344() {
    if std::env::var("HOLDOUT_6R87").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R87=1");
        return;
    }
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    assert!(ref_fasta.is_file() && bam.is_file());

    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, INTERVAL).expect("interval");
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
    let covering: Vec<_> = regions
        .iter()
        .filter(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= POS_SNP
                && r.end.get() >= POS_SNP
        })
        .collect();
    assert_eq!(covering.len(), 1);
    let outcome = HaplotypeCallerEngine::call_region(
        covering[0],
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let emitted = try_emit_call_region_variants(
        covering[0],
        &outcome,
        "SAMPLE",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .unwrap_or_default();
    let vcf = emitted
        .iter()
        .find(|r| {
            r.position == POS_SNP && r.reference == "T" && r.alternate.iter().any(|a| a == "C")
        })
        .expect("T/C");
    let snap = take_colocated_merge_numerics()
        .into_iter()
        .find(|n| n.loc == POS_SNP)
        .expect("colocated merge numerics");

    let merged_end = POS_SNP + snap.long_ref.len() as u64 - 1;
    let doc = json!({
        "java_target_reconstructed": {
            "merged_start": POS_SNP,
            "merged_end": merged_end,
            "long_ref": snap.long_ref,
            "expand": 2,
            "target": [POS_SNP - 2, merged_end + 2],
            "convention": "1-based inclusive",
        },
        "evidence": {
            "n_pairhmm": snap.n_pairhmm_reads,
            "n_overlap": snap.n_overlap_before_qname_dedupe,
        },
        "index_59": "both Java reconstructed overlap and Rust overlap reject (not causal)",
        "overlap_divergence_on_136": false,
        "vcf": {
            "gt": vcf.samples[0].gt.as_ref().map(|g| g.alleles.clone()),
            "ad": vcf.samples[0].ad.clone(),
            "pl": vcf.samples[0].pl.clone(),
            "qual": vcf.quality,
        },
    });
    println!("{}", serde_json::to_string_pretty(&doc).unwrap());

    assert_eq!(snap.long_ref, "TG");
    assert_eq!(snap.n_pairhmm_reads, 136);
    assert_eq!(snap.n_overlap_before_qname_dedupe, 60);
    assert_eq!(
        vcf.samples[0].ad.clone().unwrap_or_default(),
        vec![36u32, 19]
    );
    assert!(!vcf.alternate.iter().any(|a| a == "*"));
}
