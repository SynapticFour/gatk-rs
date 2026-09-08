//! 6R.88 holdout: AlleleLikelihoods object at DepthPerAlleleBySample for the canonical T/C site.
//!
//! Skipped unless `HOLDOUT_6R88=1`. Coordinate-free contract lives in
//! `forensic_6r88_ad_likelihoods_evidence_contract`.
//!
//! ```text
//! HOLDOUT_6R88=1 cargo test -p gatk-haplotypecaller --test holdout_6r88_ad_likelihoods_evidence -- --nocapture
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
fn holdout_6r88_ad_likelihoods_evidence_29456344() {
    if std::env::var("HOLDOUT_6R88").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R88=1");
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
        "java_ad_input_object": {
            "construction": "marginalize(alleleMapper) → retainEvidence → reuse for annotation + addEvidence(filtered,0)",
            "real_likelihood_evidence": snap.n_overlap_before_qname_dedupe,
            "alleles_at_annotateWithLikelihoods": ["TG", "*", "T", "CG"],
            "updateNonRef": "no-op without <NON_REF>",
            "inside_annotateWithLikelihoods": "marginalize remaining call alleles {TG, CG}; same evidence",
        },
        "rust_ad_input_object": {
            "construction": "overlap retainEvidence subset → 4-way pool_max → informative_ad_n_alleles → unused-ALT permute",
            "evidence": snap.n_reads,
            "alleles_4way": snap.alts,
            "long_ref": snap.long_ref,
        },
        "cardinalities": {
            "A_pairhmm": snap.n_pairhmm_reads,
            "B_retainEvidence": snap.n_overlap_before_qname_dedupe,
            "C_rust_ad_object": snap.n_reads,
            "E_rust_vcf_ad": vcf.samples[0].ad.clone(),
            "E_java_vcf_ad": [36, 19],
        },
        "b_equals_c_real_likelihood": snap.n_overlap_before_qname_dedupe == snap.n_reads,
        "vcf": {
            "gt": vcf.samples[0].gt.as_ref().map(|g| g.alleles.clone()),
            "ad": vcf.samples[0].ad.clone(),
            "pl": vcf.samples[0].pl.clone(),
            "qual": vcf.quality,
        },
        "diagnostic_remarg": snap.subset_ad_remarginalized,
        "diagnostic_permute": snap.subset_ad_permuted,
    });
    println!("{}", serde_json::to_string_pretty(&doc).unwrap());

    assert_eq!(snap.long_ref, "TG");
    assert_eq!(snap.n_pairhmm_reads, 136);
    assert_eq!(snap.n_overlap_before_qname_dedupe, 60);
    assert_eq!(
        snap.n_reads, snap.n_overlap_before_qname_dedupe,
        "Rust AD input C equals shared retainEvidence B"
    );
    assert_eq!(
        snap.n_qnames_with_multiple_overlapping_reads, 0,
        "QNAME collapse is not a C reconstruction"
    );
    assert_eq!(
        vcf.samples[0].ad.clone().unwrap_or_default(),
        vec![36u32, 19]
    );
    assert_eq!(
        snap.subset_ad_remarginalized,
        vec![36, 19],
        "6R.100: remarg of C matches Java live annotation AD"
    );
    assert!(!vcf.alternate.iter().any(|a| a == "*"));
}
