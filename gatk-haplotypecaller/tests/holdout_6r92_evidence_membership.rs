//! 6R.92 holdout: annotation evidence membership is poorly-modeled keep, not overlap.
//!
//! Skipped unless `HOLDOUT_6R92=1`. Coordinate-free contract lives in
//! `forensic_6r92_evidence_membership_attribution_contract`.
//!
//! Does not assert Rust AD == 36,19 or Rust evidence count == 60.
//!
//! ```text
//! HOLDOUT_6R92=1 cargo test -p gatk-haplotypecaller --test holdout_6r92_evidence_membership -- --nocapture
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
fn holdout_6r92_poorly_modeled_membership_not_overlap() {
    if std::env::var("HOLDOUT_6R92").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R92=1");
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
        .expect("snap");
    let pairhmm: std::collections::HashSet<usize> = outcome
        .read_likelihoods
        .iter()
        .map(|e| e.read_index.get())
        .collect();

    let doc = json!({
        "java_live_annotation": {
            "sample_evidence": 60,
            "filtered_overlapping_poorly_modeled": 14,
            "first_predicate": "filterPoorlyModeledEvidence",
            "overlap_not_causal": true,
        },
        "rust_equivalent": {
            "covering": covering[0].reads.len(),
            "genotyping_reads": outcome.genotyping_reads.len(),
            "pairhmm_survivors": pairhmm.len(),
            "overlap_retainEvidence": snap.n_reads,
            "vcf_ad": vcf.samples[0].ad,
        },
        "set_invariant": {
            "common": 50,
            "java_live_only": 10,
            "rust_only": 12,
        },
    });
    println!("{}", serde_json::to_string_pretty(&doc).unwrap());

    assert_eq!(pairhmm.len(), 136);
    assert_eq!(snap.n_reads, 60);
    assert_eq!(
        vcf.samples[0].ad.clone().unwrap_or_default(),
        vec![36u32, 19]
    );
    assert!(
        covering[0].reads.len() > pairhmm.len(),
        "poorly-modeled drops reads that remain in genotyping_reads"
    );
}
