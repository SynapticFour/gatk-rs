//! 6R.91 holdout: live Java annotation-call object at the canonical T/C site.
//!
//! Skipped unless `HOLDOUT_6R91=1`. Coordinate-free contract lives in
//! `forensic_6r91_live_ad_annotation_contract`.
//!
//! Live Java (GATK 4.4.0.0 `2dbc0258`) `DepthPerAlleleBySample.annotate`:
//! 60×4 likelihoods `TG,*,T,CG`, remaining call alleles `TG,CG`, identity remarg AD `[36,19]`.
//! This is not the 6R.88 reconstructed 62×4 remarg `[27,9]`.
//!
//! ```text
//! HOLDOUT_6R91=1 cargo test -p gatk-haplotypecaller --test holdout_6r91_live_ad_annotation -- --nocapture
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
fn holdout_6r91_live_java_annotation_state() {
    if std::env::var("HOLDOUT_6R91").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R91=1");
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
    let vcf_ad = vcf.samples[0].ad.clone().unwrap_or_default();

    let doc = json!({
        "java_live_annotation": {
            "call_alleles": ["TG", "CG"],
            "likelihoods_alleles": ["TG", "*", "T", "CG"],
            "evidence_count": 60,
            "sample": "NA12878",
            "variant_calling_subset": "20:29456342-29456347",
            "remaining_map": {"TG": ["TG"], "CG": ["CG"]},
            "independent_remarg_ad": [36, 19],
            "annotateWithLikelihoods_ad": [36, 19],
            "first_ad_write": true,
            "vcf_after_reverseTrim": {"alleles": ["T", "C"], "ad": [36, 19]},
        },
        "rust_equivalent": {
            "n_c": snap.n_reads,
            "diagnostic_remarg": snap.subset_ad_remarginalized,
            "permute": snap.subset_ad_permuted,
            "vcf_ad": vcf_ad,
        },
        "reconstructed_62x4_remarg": [27, 9],
        "live_object_explains_java_ad": true,
        "reconstructed_c_is_not_live_annotation_input": snap.n_reads != 60
            || snap.subset_ad_remarginalized != vec![36, 19],
    });
    println!("{}", serde_json::to_string_pretty(&doc).unwrap());

    assert_eq!(snap.n_reads, 60);
    assert_eq!(snap.subset_ad_remarginalized, vec![36, 19]);
    assert_eq!(vcf_ad, vec![36u32, 19]);
}
