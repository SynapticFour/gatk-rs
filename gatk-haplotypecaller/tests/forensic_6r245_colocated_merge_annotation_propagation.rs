//! 6R.245 specified the propagation contract at `20:29455649 T/TGTTTG`.
//! 6R.246 implements it: `annotation_likelihoods: subset.into_owned()` after
//! `hap_rows`. This file locks that closed contract. Do not retune PL.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r245_colocated_merge_annotation_propagation -- --nocapture --test-threads=1
//! HOLDOUT_6R245=1 cargo test -p gatk-haplotypecaller --test holdout_6r245_colocated_merge_annotation_propagation -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::hc_genotyping_engine::{
    java_alignment_read_overlaps_interval, take_colocated_merge_numerics,
    DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN,
};
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
const RUST_FALLBACK_DP: i32 = 230;
const MARGIN: i32 = DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R245\t{key}\t{}", value.as_ref());
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
fn forensic_6r245_source_propagation_contract() {
    kv("java_pin", JAVA_PIN);
    kv(
        "production_change",
        "6R.246 annotation_likelihoods: subset.into_owned()",
    );
    assert_eq!(DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, 128);

    let assign = fs::read_to_string(
        repo_root().join("gatk-haplotypecaller/src/hc_genotyping_engine/genotype_assign.rs"),
    )
    .unwrap();
    let merge = fn_body(&assign, "fn try_genotype_colocated_snp_indel_merge");
    assert_eq!(
        merge
            .matches("ColocatedMergeGenotype::Call(GenotypedSiteCall")
            .count(),
        1,
        "exactly one successful Call constructor"
    );
    assert_eq!(
        merge
            .matches("annotation_likelihoods: subset.into_owned()")
            .count(),
        1
    );
    assert!(
        !merge.contains("annotation_likelihoods: Vec::new()"),
        "empty sentinel must not remain on the Call constructor"
    );
    assert!(merge.contains("let subset = likelihood_subset_for_event"));
    assert!(merge.contains("let hap_rows = region_likelihoods_to_rows(subset.as_ref()"));
    let hap_rows_pos = merge
        .find("let hap_rows = region_likelihoods_to_rows(subset.as_ref()")
        .expect("hap_rows");
    let owned_pos = merge
        .find("annotation_likelihoods: subset.into_owned()")
        .expect("into_owned");
    assert!(
        hap_rows_pos < owned_pos,
        "hap_rows borrows subset.as_ref() before into_owned moves it"
    );
    let hap_line_end = hap_rows_pos
        + merge[hap_rows_pos..]
            .find('\n')
            .expect("hap_rows statement");
    let after_hap = &merge[hap_line_end..owned_pos];
    assert!(
        !after_hap.contains("subset.as_ref()") && !after_hap.contains("let subset ="),
        "GL/AD/QUAL consume hap_rows; into_owned is the next use of subset"
    );
    assert!(merge.contains("return Ok(ColocatedMergeGenotype::NotApplicable)"));
    assert!(merge.contains("return Ok(ColocatedMergeGenotype::MergedNoEmit)"));
    assert!(
        !merge.contains("with_annotation_likelihoods"),
        "Call uses the struct literal, not the SiteScore helper"
    );
    kv(
        "call_constructors",
        "one Call(GenotypedSiteCall { … subset.into_owned() }); NotApplicable and MergedNoEmit do not produce a Call",
    );

    let subset_fn = fs::read_to_string(
        repo_root().join("gatk-haplotypecaller/src/hc_genotyping_engine/mod.rs"),
    )
    .unwrap();
    assert!(subset_fn.contains("fn likelihood_subset_for_event<'a>("));
    assert!(subset_fn.contains("-> Cow<'a, [RegionReadLikelihood]>"));
    assert!(
        subset_fn.contains("return Cow::Owned(subset);")
            && subset_fn.contains("if config.enable_java_strict()")
    );
    assert!(subset_fn.contains("pub annotation_likelihoods: Vec<RegionReadLikelihood>"));
    assert!(subset_fn.contains(
        "pub fn with_annotation_likelihoods(mut self, likelihoods: Vec<RegionReadLikelihood>)"
    ));
    assert!(
        subset_fn.contains("Empty → annotations keep the region-wide PairHMM matrix."),
        "empty Vec is the emit-fallback sentinel"
    );
    kv(
        "annotation_row_type",
        "Vec<RegionReadLikelihood> (owned cells: read_index, haplotype_index, log10_likelihood)",
    );
    kv(
        "subset_type",
        "Cow<[RegionReadLikelihood]>; strict_java returns Cow::Owned from filter_likelihoods_for_variant (.cloned() already paid for GLs)",
    );
    kv(
        "ownership_after_hap_rows",
        "hap_rows is Vec<ReadLikelihoodRow> copied from subset.as_ref(); subset.into_owned() is a Vec move on Cow::Owned, not a second matrix clone",
    );

    let rrl =
        fs::read_to_string(repo_root().join("gatk-haplotypecaller/src/region_read_likelihood.rs"))
            .unwrap();
    assert!(rrl.contains("pub struct RegionReadLikelihood"));
    assert!(rrl.contains("pub read_index: ReadIndex"));
    assert!(rrl.contains("pub haplotype_index: HaplotypeIndex"));
    assert!(rrl.contains("pub log10_likelihood: f64"));

    let emit = fs::read_to_string(repo_root().join("gatk-haplotypecaller/src/region_vcf_emit.rs"))
        .unwrap();
    let for_call = fn_body(&emit, "fn annotation_likelihoods_for_call");
    assert!(for_call.contains("if call.annotation_likelihoods.is_empty()"));
    assert!(for_call.contains("likelihoods: region_wide.likelihoods"));
    assert!(for_call.contains("likelihoods: &call.annotation_likelihoods"));
    kv(
        "sentinel",
        "empty → region-wide; non-empty → &call.annotation_likelihoods. Replacing Vec::new() with subset.into_owned() is sufficient to suppress fallback",
    );

    let coverage = fs::read_to_string(
        repo_root().join("gatk-haplotypecaller/src/variant_site_hc_annotations.rs"),
    )
    .unwrap();
    assert!(coverage.contains("pub fn coverage_evidence_count"));
    assert!(coverage.contains("seen.insert(cell.read_index.get())"));
    assert!(coverage.contains("Some(ev) => coverage_evidence_count("));
    kv(
        "coverage",
        "INFO DP = unique read_index of the annotation_likelihoods_for_call matrix; no second overlap filter",
    );
    kv(
        "proposed_6r246",
        "try_genotype_colocated_snp_indel_merge Call: annotation_likelihoods: subset.into_owned() after hap_rows; do not change merged_handled_locs, Coverage, fallback, PairHMM, PL",
    );
    kv(
        "contamination",
        "off; merge has no contaminationDownsampling",
    );
}

#[test]
fn forensic_6r245_live_subset_is_sufficient_for_info_dp_123() {
    if std::env::var("FORENSIC_6R245_LIVE").ok().as_deref() != Some("1") {
        eprintln!("skip: pre-6R.246 live empty-attach counterfactual");
        return;
    }
    kv("java_pin", JAVA_PIN);
    kv(
        "production_change",
        "6R.246 annotation_likelihoods: subset.into_owned()",
    );
    let root = repo_root();
    if !root.join(JAVA_VCF_REL).is_file()
        || !root.join(REF_REL).is_file()
        || !root.join(BAM_REL).is_file()
    {
        eprintln!("skip: missing chr20_tiny BAM/REF");
        return;
    }
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
    assert!(java_line.contains("DP=123"));

    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing chr20_tiny BAM/REF");
        return;
    }
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
    assert!(
        call.annotation_likelihoods.is_empty(),
        "production still empty; this arrow does not attach"
    );
    kv(
        "current_call_ann_n",
        call.annotation_likelihoods.len().to_string(),
    );

    let snap = take_colocated_merge_numerics()
        .into_iter()
        .find(|s| s.loc == TARGET)
        .expect("merge snapshot is the live subset identity");
    kv(
        "live_merge_subset",
        format!(
            "n_reads={} n_overlap={} n_overlap_qnames={} n_pairhmm={} subset_ad={:?}",
            snap.n_reads,
            snap.n_overlap_before_qname_dedupe,
            snap.n_overlap_unique_qnames,
            snap.n_pairhmm_reads,
            snap.subset_ad_remarginalized
        ),
    );
    assert_eq!(snap.n_reads, JAVA_INFO_DP as usize);
    assert_eq!(snap.n_overlap_before_qname_dedupe, JAVA_INFO_DP as usize);
    assert_eq!(snap.n_overlap_unique_qnames, JAVA_INFO_DP as usize);
    assert_eq!(snap.n_pairhmm_reads, RUST_FALLBACK_DP as usize);
    assert_eq!(
        snap.subset_ad_remarginalized,
        vec![88, 22],
        "GL/AD identities are the same 123-read subset"
    );
    assert_eq!(call.genotype.format.ad_as_i32(), vec![88, 22]);

    let keep: BTreeSet<usize> = unique_indices(&outcome.read_likelihoods)
        .into_iter()
        .filter(|&i| {
            outcome
                .genotyping_reads
                .get(i)
                .is_some_and(|r| java_alignment_read_overlaps_interval(r, TARGET, TARGET, MARGIN))
        })
        .collect();
    assert_eq!(keep.len(), JAVA_INFO_DP as usize);
    let subset_identity: Vec<RegionReadLikelihood> = outcome
        .read_likelihoods
        .iter()
        .filter(|c| keep.contains(&c.read_index.get()))
        .cloned()
        .collect();
    kv(
        "subset_identity_cells_unique",
        format!(
            "cells={} unique={}",
            subset_identity.len(),
            unique_indices(&subset_identity).len()
        ),
    );
    assert_eq!(
        unique_indices(&subset_identity).len(),
        JAVA_INFO_DP as usize
    );
    let desired_dp = coverage_evidence_count(
        &outcome.genotyping_reads,
        &subset_identity,
        TARGET,
        TARGET,
        MARGIN,
    );
    kv("desired_coverage_evidence_count", desired_dp.to_string());
    assert_eq!(
        desired_dp, JAVA_INFO_DP,
        "Coverage.evidenceCount of the subset identity is Java INFO DP 123"
    );
    let current_attached = coverage_evidence_count(
        &outcome.genotyping_reads,
        &call.annotation_likelihoods,
        TARGET,
        TARGET,
        MARGIN,
    );
    assert_eq!(current_attached, 0);
    let region_wide = coverage_evidence_count(
        &outcome.genotyping_reads,
        &outcome.read_likelihoods,
        TARGET,
        TARGET,
        MARGIN,
    );
    assert_eq!(region_wide, RUST_FALLBACK_DP);

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
    assert_eq!(live_dp, RUST_FALLBACK_DP);

    kv(
        "propagation_sufficient",
        "yes: non-empty annotation_likelihoods selects &call.annotation_likelihoods; Coverage unique read_index of subset identity is 123; fallback is not required",
    );
    kv("expected_info_dp_after_6r246", "230 → 123");
    kv(
        "no_second_pairhmm",
        "subset is already the post-retainEvidence haplotype-cell matrix; into_owned moves it",
    );
    kv(
        "classification",
        "ANNOTATION_LIKELIHOOD_LIFECYCLE_DIVERGENCE",
    );
    kv("next_arrow", "6R.246");
}
