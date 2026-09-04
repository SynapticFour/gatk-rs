//! 6R.86 holdout: AD evidence population at the canonical T/C site.
//!
//! Skipped unless `HOLDOUT_6R86=1`. Coordinate-free proof lives in
//! `forensic_6r86_ad_evidence_contract`.
//!
//! ```text
//! HOLDOUT_6R86=1 cargo test -p gatk-haplotypecaller --test holdout_6r86_ad_evidence -- --nocapture
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
fn holdout_6r86_ad_evidence_29456344() {
    if std::env::var("HOLDOUT_6R86").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R86=1");
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

    let doc = json!({
        "vcf": {
            "gt": vcf.samples[0].gt.as_ref().map(|g| g.alleles.clone()),
            "ad": vcf.samples[0].ad.clone(),
            "pl": vcf.samples[0].pl.clone(),
            "qual": vcf.quality,
        },
        "java_oracle": {"gt": [0, 1], "ad": [36, 19], "pl": [542, 0, 1353], "qual": 510.06},
        "evidence": {
            "n_pairhmm": snap.n_pairhmm_reads,
            "n_overlap": snap.n_overlap_before_qname_dedupe,
            "n_after_qname": snap.n_reads,
            "n_multi_qname": snap.n_qnames_with_multiple_overlapping_reads,
            "merged_ad_4way": snap.merged_ad,
            "remarg": snap.subset_ad_remarginalized,
            "permute": snap.subset_ad_permuted,
        },
        "first_divergence": "read eligibility: overlap retainEvidence 136→62; QNAME collapse not causal",
    });
    println!("{}", serde_json::to_string_pretty(&doc).unwrap());

    assert_eq!(
        vcf.samples[0].ad.clone().unwrap_or_default(),
        vec![26u32, 9]
    );
    assert_eq!(snap.n_qnames_with_multiple_overlapping_reads, 0);
    assert_eq!(snap.n_overlap_before_qname_dedupe, snap.n_reads);
    assert_eq!(snap.subset_ad_remarginalized, vec![27, 9]);
    assert_eq!(snap.subset_ad_permuted, vec![26, 9]);
    assert!(!vcf.alternate.iter().any(|a| a == "*"));
}
