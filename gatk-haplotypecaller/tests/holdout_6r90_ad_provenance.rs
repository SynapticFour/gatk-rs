//! 6R.90 holdout: AD provenance after unused-ALT subset / reverseTrim.
//!
//! Skipped unless `HOLDOUT_6R90=1`. Coordinate-free contract lives in
//! `forensic_6r90_ad_provenance_contract`.
//!
//! ```text
//! HOLDOUT_6R90=1 cargo test -p gatk-haplotypecaller --test holdout_6r90_ad_provenance -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, reverse_trim_alleles,
    take_colocated_merge_numerics, traverse_assembly_region_walker, try_emit_call_region_variants,
    AssemblyRegionCallDisposition, CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams,
    WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
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
fn holdout_6r90_ad_provenance_after_annotation() {
    if std::env::var("HOLDOUT_6R90").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R90=1");
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

    let (trim_ref, trim_alts) = reverse_trim_alleles(&snap.long_ref, &["CG".to_string()]);
    let vcf_ad = vcf.samples[0].ad.clone().unwrap_or_default();

    let doc = json!({
        "java_lifecycle": {
            "calculateGLs": "PL only; AD absent",
            "subsetAlleles": "slice AD iff hasAD; therefore AD still absent",
            "DepthPerAlleleBySample": "FIRST AD write on remaining call alleles",
            "reverseTrimAlleles": "copy AD; remap TG/CG → T/C",
            "phaseVC": "copy AD",
        },
        "rust_lifecycle": {
            "four_way_informative": snap.merged_ad,
            "unused_alt_slice": snap.subset_ad_permuted,
            "diagnostic_remarg": snap.subset_ad_remarginalized,
            "reverse_trim": format!("{trim_ref}/{}", trim_alts.join(",")),
            "vcf_ad": vcf_ad,
        },
        "java_vcf_ad": [36, 19],
        "reverse_trim_preserves_remarg": snap.subset_ad_remarginalized == vec![36, 19],
        "vcf_equals_remarg_not_permute": vcf_ad.iter().map(|&x| x as i32).collect::<Vec<_>>()
            == snap.subset_ad_remarginalized,
        "n_c": snap.n_reads,
    });
    println!("{}", serde_json::to_string_pretty(&doc).unwrap());

    assert_eq!(snap.n_reads, 60);
    assert_eq!(trim_ref, "T");
    assert_eq!(trim_alts, vec!["C".to_string()]);
    assert_eq!(snap.subset_ad_remarginalized, vec![36, 19]);
    assert_eq!(vcf_ad, vec![36u32, 19]);
    assert_eq!(snap.subset_ad_permuted, vec![34, 17]);
    assert_eq!(
        vcf_ad.iter().map(|&x| x as i32).collect::<Vec<_>>(),
        snap.subset_ad_remarginalized
    );
}
