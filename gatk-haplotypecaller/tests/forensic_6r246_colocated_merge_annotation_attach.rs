//! 6R.246: production. Colocated merge at `20:29455649 T/TGTTTG` propagates
//! the existing retainEvidence `subset` onto `Call.annotation_likelihoods`
//! (`subset.into_owned()` after `hap_rows`). INFO DP 230 → 123.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! Do not retune PL. Do not change Coverage or fallback semantics.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r246_colocated_merge_annotation_attach -- --nocapture --test-threads=1
//! HOLDOUT_6R246=1 cargo test -p gatk-haplotypecaller --test holdout_6r246_colocated_merge_annotation_attach -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::genotyping::best_pl_index;
use gatk_haplotypecaller::hc_genotyping_engine::take_colocated_merge_numerics;
use gatk_haplotypecaller::read_threading_assembler::DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH;
use gatk_haplotypecaller::variant_site_hc_annotations::coverage_evidence_count;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, ReadFilterParams, RegionReadLikelihood, WalkerTraversalConfig,
    DEFAULT_STAND_EMIT_CONFIDENCE,
};
use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_VCF_REL: &str = "parity/reports/6r43/chr20_tiny/java.vcf";
const TARGET: u64 = 29_455_649;
const TARGET_REF: &str = "T";
const TARGET_ALT: &str = "TGTTTG";
const JAVA_INFO_DP: i32 = 123;
const HAP_N: usize = 24;
const CELL_N: usize = 2952;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R246\t{key}\t{}", value.as_ref());
}

fn info_i32(info: &[InfoValue], key: &str) -> Option<i32> {
    for v in info {
        if let InfoValue::Integer(k, xs) = v {
            if k == key {
                return xs.first().copied();
            }
        }
    }
    None
}

fn unique_indices(likelihoods: &[RegionReadLikelihood]) -> BTreeSet<usize> {
    likelihoods.iter().map(|c| c.read_index.get()).collect()
}

fn fn_body<'a>(src: &'a str, sig: &str) -> &'a str {
    let start = src.find(sig).unwrap_or(0);
    let rest = &src[start..];
    let rel = rest[sig.len()..]
        .find("\nfn ")
        .map(|i| sig.len() + i)
        .unwrap_or(rest.len());
    &rest[..rel]
}

#[test]
fn forensic_6r246_source_is_subset_into_owned_after_hap_rows() {
    kv("java_pin", JAVA_PIN);
    kv(
        "production_change",
        "try_genotype_colocated_snp_indel_merge Call annotation_likelihoods: subset.into_owned()",
    );
    assert_eq!(DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, 128);
    assert_eq!(JAVA_INFO_DP as usize * HAP_N, CELL_N);

    let assign = fs::read_to_string(
        repo_root().join("gatk-haplotypecaller/src/hc_genotyping_engine/genotype_assign.rs"),
    )
    .unwrap();
    let merge = fn_body(&assign, "fn try_genotype_colocated_snp_indel_merge");
    assert!(merge.contains("let hap_rows = region_likelihoods_to_rows(subset.as_ref()"));
    assert!(merge.contains("annotation_likelihoods: subset.into_owned()"));
    assert!(
        !merge.contains("annotation_likelihoods: Vec::new()"),
        "empty sentinel must not remain on the Call constructor"
    );
    let hap_pos = merge
        .find("let hap_rows = region_likelihoods_to_rows(subset.as_ref()")
        .expect("hap_rows");
    let owned_pos = merge
        .find("annotation_likelihoods: subset.into_owned()")
        .expect("into_owned");
    assert!(hap_pos < owned_pos, "into_owned is after hap_rows");
    assert_eq!(
        merge
            .matches("ColocatedMergeGenotype::Call(GenotypedSiteCall")
            .count(),
        1
    );

    let emit = fs::read_to_string(repo_root().join("gatk-haplotypecaller/src/region_vcf_emit.rs"))
        .unwrap();
    let for_call = fn_body(&emit, "fn annotation_likelihoods_for_call");
    assert!(for_call.contains("if call.annotation_likelihoods.is_empty()"));
    assert!(for_call.contains("likelihoods: &call.annotation_likelihoods"));
    kv(
        "fallback_semantics",
        "unchanged; non-empty selects Call-local object",
    );
}

#[test]
fn forensic_6r246_live_info_dp_is_123_not_region_wide_230() {
    kv("java_pin", JAVA_PIN);
    let root = repo_root();
    let java_line = fs::read_to_string(root.join(JAVA_VCF_REL))
        .unwrap()
        .lines()
        .find(|l| {
            let f: Vec<_> = l.split('\t').collect();
            f.len() >= 10
                && f[0] == "20"
                && f[1] == "29455649"
                && f[3] == TARGET_REF
                && f[4] == TARGET_ALT
        })
        .expect("Java covering VCF has T/TGTTTG")
        .to_string();
    kv("java_vcf_line", &java_line);
    assert!(java_line.contains("DP=123"));
    assert!(java_line.contains("570,0,3517"));
    let java_fields: Vec<_> = java_line.split('\t').collect();
    assert_eq!(java_fields[5], "562.60");
    assert!(java_fields[9].starts_with("0/1:88,22:110:99:570,0,3517"));

    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
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
    let region = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= TARGET
                && r.end.get() >= TARGET
        })
        .expect("ActiveFull");
    let outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    assert_eq!(outcome.assembly.haplotypes.len(), HAP_N);

    let call = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based.get() == TARGET
                && c.event.ref_allele == TARGET_REF
                && c.event.alt_allele == TARGET_ALT
        })
        .expect("T/TGTTTG");
    assert!(call.post_merge_unused_alt_subset);
    let unique = unique_indices(&call.annotation_likelihoods);
    kv(
        "call_annotation",
        format!(
            "cells={} unique={} hap_n={}",
            call.annotation_likelihoods.len(),
            unique.len(),
            HAP_N
        ),
    );
    assert_eq!(
        unique.len(),
        JAVA_INFO_DP as usize,
        "unique read_index is Java annotation evidence 123, not Vec len"
    );
    assert_eq!(
        call.annotation_likelihoods.len(),
        CELL_N,
        "propagated object is 123×24 haplotype cells, not region-wide 230"
    );
    assert_eq!(JAVA_INFO_DP as usize * HAP_N, CELL_N);
    assert_ne!(
        unique.len(),
        unique_indices(&outcome.read_likelihoods).len(),
        "must not be the region-wide 230-read matrix"
    );

    let attached_dp = coverage_evidence_count(
        &outcome.genotyping_reads,
        &call.annotation_likelihoods,
        TARGET,
        TARGET,
        2,
    );
    kv("attached_coverage_evidence_count", attached_dp.to_string());
    assert_eq!(attached_dp, JAVA_INFO_DP);

    let region_wide = coverage_evidence_count(
        &outcome.genotyping_reads,
        &outcome.read_likelihoods,
        TARGET,
        TARGET,
        2,
    );
    kv("region_wide_still_230", region_wide.to_string());
    assert_eq!(region_wide, 230);

    let snap = take_colocated_merge_numerics()
        .into_iter()
        .find(|s| s.loc == TARGET)
        .expect("merge snapshot");
    assert_eq!(snap.n_overlap_before_qname_dedupe, JAVA_INFO_DP as usize);
    assert_eq!(snap.n_pairhmm_reads, 230);

    let pl = call.genotype.format.pl_as_i32();
    assert_eq!(best_pl_index(&call.genotype.format.pl), 1);
    assert_eq!(call.genotype.format.ad_as_i32(), vec![88, 22]);
    assert_eq!(&pl[..2], &[570, 0]);
    assert!((pl[2] - 3517).abs() <= 1, "PL ±1 retained; got {}", pl[2]);
    assert_eq!(call.genotype.format.gq.as_i32(), 99);
    assert_eq!(call.genotype.format.dp.as_i32(), 110);

    let recs =
        try_emit_call_region_variants(region, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    let emitted = recs
        .iter()
        .find(|r| {
            r.position == TARGET
                && r.reference == TARGET_REF
                && r.alternate.iter().any(|a| a == TARGET_ALT)
        })
        .expect("emit");
    let live_dp = info_i32(&emitted.info, "DP").expect("INFO DP");
    kv("live_emit_info_dp", live_dp.to_string());
    kv("live_emit_qual", format!("{:?}", emitted.quality));
    assert_eq!(live_dp, JAVA_INFO_DP);
    let qual = emitted.quality.expect("QUAL");
    assert!((qual - 562.60).abs() < 0.01, "QUAL unchanged; got {qual}");
    let sample = emitted.samples.first().expect("sample");
    assert_eq!(
        sample.gt.as_ref().map(|g| g.alleles.as_slice()),
        Some(&[0, 1][..])
    );
    assert_eq!(sample.ad.as_deref(), Some(&[88u32, 22][..]));
    assert_eq!(sample.dp, Some(110));
    assert_eq!(sample.gq, Some(99.0));
    let emit_pl = sample.pl.as_deref().expect("PL");
    assert_eq!(&emit_pl[..2], &[570, 0]);
    assert!((emit_pl[2] as i32 - 3517).abs() <= 1);

    kv(
        "classification",
        "ANNOTATION_LIKELIHOOD_LIFECYCLE_DIVERGENCE",
    );
    kv(
        "first_divergent_operation",
        "ColocatedMergeGenotype::Call construction (now into_owned)",
    );
}
