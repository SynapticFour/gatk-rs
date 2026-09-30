//! 6R.244 showed the 123-read subset was in scope and not attached at
//! `20:29455649 T/TGTTTG`. 6R.246 attaches it with
//! `annotation_likelihoods: subset.into_owned()`. This file locks that
//! closed contract. Do not retune PL.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r244_colocated_merge_annotation_lifecycle -- --nocapture --test-threads=1
//! HOLDOUT_6R244=1 cargo test -p gatk-haplotypecaller --test holdout_6r244_colocated_merge_annotation -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{merged_alleles_for_genotyping, merged_site_uses_joint_gls};
use gatk_haplotypecaller::genotyping::best_pl_index;
use gatk_haplotypecaller::hc_genotyping_engine::{
    java_emit_af_decision, take_colocated_merge_numerics, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use gatk_haplotypecaller::read_threading_assembler::DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH;
use gatk_haplotypecaller::variant_site_hc_annotations::coverage_evidence_count;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig, DEFAULT_STAND_CALL_CONF,
};
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

/// Java 4.4 `HaplotypeCallerGenotypingEngine.prepareReadAlleleLikelihoodsForAnnotation`
/// (pinned SHA). Contamination off / `useFilteredReadMapForAnnotations` reuses
/// the genotyping object; it does not construct a new AlleleLikelihoods.
const JAVA_PREPARE_ANN_REUSE: &str = "\
if (hcArgs.useFilteredReadMapForAnnotations || !configuration.isSampleContaminationPresent()) {\n\
readAlleleLikelihoodsForAnnotations = readAlleleLikelihoodsForGenotyping;\n\
}";

/// Java loc loop after `calculateGenotypes` non-null: prepare then annotate.
const JAVA_CALL_THEN_ANNOTATE: &str = "\
if( call != null ) {\n\
readAlleleLikelihoods = prepareReadAlleleLikelihoodsForAnnotation(";

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R244\t{key}\t{}", value.as_ref());
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
fn forensic_6r244_source_merge_drops_in_scope_subset() {
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
    kv(
        "merge_fn_has_subset",
        format!(
            "{}",
            merge.contains("let subset = likelihood_subset_for_event")
        ),
    );
    kv(
        "merge_fn_annotation",
        format!(
            "{}",
            merge.contains("annotation_likelihoods: subset.into_owned()")
        ),
    );
    kv(
        "merge_fn_with_annotation",
        format!("{}", merge.contains("with_annotation_likelihoods")),
    );

    assert!(
        merge.contains("let subset = likelihood_subset_for_event"),
        "merge receives/builds the genotyping AlleleLikelihoods as local `subset`"
    );
    assert!(
        merge.contains("if subset.is_empty()"),
        "empty subset is NotApplicable, not an empty-annotation Call"
    );
    let subset_pos = merge
        .find("let subset = likelihood_subset_for_event")
        .expect("subset");
    let owned_pos = merge
        .find("annotation_likelihoods: subset.into_owned()")
        .expect("into_owned");
    assert!(
        subset_pos < owned_pos,
        "123-read subset is in scope when the Call attaches it"
    );
    assert!(
        !merge.contains("with_annotation_likelihoods"),
        "merge never calls GenotypedSiteCall::with_annotation_likelihoods"
    );
    assert_eq!(
        merge
            .matches("annotation_likelihoods: subset.into_owned()")
            .count(),
        1,
        "every successful colocated merge Call attaches the same subset"
    );
    assert!(
        !merge.contains("annotation_likelihoods: Vec::new()"),
        "empty sentinel must not remain on the Call constructor"
    );
    assert!(
        merge.contains("ColocatedMergeGenotype::Call(GenotypedSiteCall"),
        "result type is GenotypedSiteCall, which can carry annotation_likelihoods"
    );
    assert!(
        !merge.contains("subset.clone()"),
        "into_owned moves the Cow; it does not clone the matrix"
    );

    let walker = fn_body(&assign, "pub fn assign_genotype_likelihoods_for_region");
    assert!(
        walker.contains("merged_handled_locs.insert(loc)"),
        "Call and MergedNoEmit both record the loc"
    );
    assert!(
        walker.contains("if !merged_site_handled"),
        "biallelic SiteScore walk runs only when merge is NotApplicable"
    );
    let handled_if = walker
        .find("if !merged_site_handled")
        .expect("merged_site_handled guard");
    let first_variation = walker[handled_if..]
        .find("try_genotype_variation_event")
        .expect("SiteScore after guard");
    kv(
        "site_score_is_behind_merged_handled_guard",
        format!("{}", first_variation > 0),
    );
    assert!(
        walker.contains("merged_handled_locs.contains(&loc)")
            || walker.contains("merged_handled_locs.contains(&event.start_1based.get())"),
        "later supplement walks also skip handled locs"
    );

    let site_call = fs::read_to_string(
        repo_root().join("gatk-haplotypecaller/src/hc_genotyping_engine/mod.rs"),
    )
    .unwrap();
    assert!(
        site_call.contains("pub annotation_likelihoods: Vec<RegionReadLikelihood>"),
        "GenotypedSiteCall can carry annotation evidence"
    );
    assert!(
        site_call.contains("Empty → annotations keep the region-wide PairHMM matrix."),
        "empty Vec is the emit-fallback sentinel, not a merge-specific clear"
    );
    assert!(
        site_call.contains("pub fn with_annotation_likelihoods"),
        "attach helper exists; merge never uses it"
    );

    let pipeline = fs::read_to_string(
        repo_root().join("gatk-haplotypecaller/src/hc_genotyping_engine/genotype_site_pipeline.rs"),
    )
    .unwrap();
    assert!(
        pipeline
            .contains("GenotypedSiteCall::new(event, gt).with_annotation_likelihoods(annotation)"),
        "SiteScore attach path exists and is the Java prepare-ann analogue"
    );

    let emit = fs::read_to_string(repo_root().join("gatk-haplotypecaller/src/region_vcf_emit.rs"))
        .unwrap();
    assert!(
        emit.contains("if call.annotation_likelihoods.is_empty()")
            && emit.contains("fn annotation_likelihoods_for_call"),
        "empty annotation_likelihoods selects region-wide fallback"
    );

    kv(
        "empty_is_sentinel",
        "empty Vec means emit uses region-wide PairHMM; it is not computed from retainEvidence==0",
    );
    kv(
        "intentional_scope",
        "all successful ColocatedMergeGenotype::Call constructors; not a T/TGTTTG-only branch",
    );
}

#[test]
fn forensic_6r244_java_lifecycle_reuses_genotyping_object() {
    kv("java_pin", JAVA_PIN);
    kv(
        "java_source_sha",
        "2dbc025821bc5f686c423ff332a41e6cef892a77",
    );
    kv(
        "contamination",
        "off (default HC; isSampleContaminationPresent=false)",
    );
    kv(
        "java_unified_loc_loop",
        "makeMergedVariantContext → marginalize → retainEvidence(±2) → calculateGLsForThisEvent → calculateGenotypes → prepareReadAlleleLikelihoodsForAnnotation → makeAnnotatedCall → annotateContext → Coverage.annotate",
    );
    kv(
        "java_new_vs_reuse",
        "default HC reuses readAlleleLikelihoodsForGenotyping; does not construct a new AlleleLikelihoods",
    );
    kv(
        "java_merged_vc_carries_evidence",
        "no; VariantContext holds alleles/genotypes; evidence stays on AlleleLikelihoods passed to annotateContext",
    );
    kv(
        "java_colocated_vs_biallelic",
        "Java has one loc loop. makeMergedVariantContext always runs. There is no separate merge-vs-SiteScore split.",
    );
    assert!(JAVA_PREPARE_ANN_REUSE
        .contains("readAlleleLikelihoodsForAnnotations = readAlleleLikelihoodsForGenotyping"));
    assert!(JAVA_CALL_THEN_ANNOTATE.contains("prepareReadAlleleLikelihoodsForAnnotation"));
    kv(
        "java_prepare_ann_reuse_quote",
        JAVA_PREPARE_ANN_REUSE.replace('\n', " "),
    );
}

#[test]
fn forensic_6r244_live_merge_result_loses_123_then_handled_locs_skip_sitescore() {
    if std::env::var("FORENSIC_6R244_LIVE").ok().as_deref() != Some("1") {
        eprintln!("skip: pre-6R.246 live empty-attach counterfactual");
        return;
    }
    kv("java_pin", JAVA_PIN);
    kv(
        "production_change",
        "6R.246 annotation_likelihoods: subset.into_owned()",
    );
    let root = repo_root();
    let java_vcf = root.join(JAVA_VCF_REL);
    if !java_vcf.is_file() || !root.join(REF_REL).is_file() || !root.join(BAM_REL).is_file() {
        eprintln!("skip: missing chr20_tiny BAM/REF");
        return;
    }
    let java_line = fs::read_to_string(&java_vcf)
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
    assert!(java_line.contains("DP=123"), "Java INFO DP=123");
    assert!(java_line.contains("570,0,3517"), "Java PL 570,0,3517");
    let java_fields: Vec<_> = java_line.split('\t').collect();
    kv("java_qual", java_fields[5]);
    kv("java_gt_ad_pl", java_fields[9]);

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
        .expect("ActiveFull containing target");

    let outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");

    let events_at: Vec<_> = outcome
        .assembly
        .variation_events()
        .iter()
        .filter(|e| e.start_1based.get() == TARGET)
        .cloned()
        .collect();
    assert_eq!(events_at.len(), 2);
    assert!(events_at.iter().any(|e| e.alt_allele == TARGET_ALT));
    assert!(events_at.iter().any(|e| e.alt_allele == "TTTG"));
    let merged = merged_alleles_for_genotyping(&events_at, TARGET).expect("merged alleles");
    kv(
        "merged_objects",
        format!(
            "EventMap T/TGTTTG + T/TTTG → merged REF={} ALTS={}",
            merged.0,
            merged.1.join(",")
        ),
    );
    assert!(merged_site_uses_joint_gls(
        &events_at, TARGET, &merged.0, &merged.1
    ));

    let calls_at_loc: Vec<_> = outcome
        .genotyped_calls
        .iter()
        .filter(|c| c.event.start_1based.get() == TARGET)
        .collect();
    kv("genotyped_calls_at_loc", calls_at_loc.len().to_string());
    assert_eq!(
        calls_at_loc.len(),
        1,
        "merged_handled_locs skipped SiteScore; a second biallelic call was not produced"
    );

    let call = calls_at_loc[0];
    assert_eq!(call.event.ref_allele, TARGET_REF);
    assert_eq!(call.event.alt_allele, TARGET_ALT);
    let pl = call.genotype.format.pl_as_i32();
    let ad = call.genotype.format.ad_as_i32();
    kv(
        "merge_result_object",
        format!(
            "VC={}/{} GT_idx={} AD={:?} PL={:?} GQ={} QUAL_log10_p_error={:?} extra_alt={:?} post_merge={} ann_n={} ann_is_empty={}",
            call.event.ref_allele,
            call.event.alt_allele,
            best_pl_index(&call.genotype.format.pl),
            ad,
            pl,
            call.genotype.format.gq.as_i32(),
            call.qual_log10_p_error,
            call.extra_alt_alleles,
            call.post_merge_unused_alt_subset,
            call.annotation_likelihoods.len(),
            call.annotation_likelihoods.is_empty()
        ),
    );
    assert_eq!(best_pl_index(&call.genotype.format.pl), 1);
    assert_eq!(ad, vec![88, 22]);
    assert_eq!(&pl[..2], &[570, 0]);
    assert!(
        (pl[2] - 3517).abs() <= 1,
        "PL ±1 is representation; got {}",
        pl[2]
    );
    assert!(call.post_merge_unused_alt_subset);
    assert!(
        call.extra_alt_alleles.is_empty(),
        "no other field holds leftover alleles or evidence"
    );
    assert!(
        call.annotation_likelihoods.is_empty(),
        "merge result annotation object is empty; do not inject a fake one"
    );
    assert!(
        call.qual_log10_p_error.is_none(),
        "no SPAN_DEL; QUAL comes from AF on unused-ALT GLs like Java"
    );

    let snaps = take_colocated_merge_numerics();
    let snap = snaps
        .iter()
        .find(|s| s.loc == TARGET)
        .expect("colocated merge snapshot at target");
    kv(
        "pre_merge_annotation_evidence",
        format!(
            "n_reads={} n_overlap={} n_overlap_qnames={} n_pairhmm={}",
            snap.n_reads,
            snap.n_overlap_before_qname_dedupe,
            snap.n_overlap_unique_qnames,
            snap.n_pairhmm_reads
        ),
    );
    assert_eq!(
        snap.n_overlap_before_qname_dedupe, JAVA_INFO_DP as usize,
        "pre-merge retainEvidence unique n=123 exists inside the merge function"
    );
    assert_eq!(snap.n_reads, JAVA_INFO_DP as usize);
    assert_eq!(snap.n_overlap_unique_qnames, JAVA_INFO_DP as usize);
    assert_eq!(snap.n_pairhmm_reads, RUST_FALLBACK_DP as usize);
    kv(
        "subset_available_but_not_on_call",
        "local subset / n_overlap=123 is used for GLs/AD then dropped; snapshot is forensic-only TLS, not annotation_likelihoods",
    );

    let attached_cov = coverage_evidence_count(
        &outcome.genotyping_reads,
        &call.annotation_likelihoods,
        TARGET,
        TARGET,
        2,
    );
    kv("merge_result_annotation_evidence", attached_cov.to_string());
    assert_eq!(attached_cov, 0);

    let region_wide = coverage_evidence_count(
        &outcome.genotyping_reads,
        &outcome.read_likelihoods,
        TARGET,
        TARGET,
        2,
    );
    kv("region_wide_fallback_count", region_wide.to_string());
    assert_eq!(region_wide, RUST_FALLBACK_DP);

    let gls = &call.genotype.genotype_log10_likelihoods;
    let af = java_emit_af_decision(gls, DEFAULT_STAND_CALL_CONF).expect("af");
    kv(
        "emit_pred",
        format!(
            "stand30_pass={} phred={:.2} mono={}",
            af.passes_emit, af.phred_scaled, af.site_is_monomorphic
        ),
    );
    assert!(af.passes_emit);

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
        .expect("both engines emit");
    let live_dp = info_i32(&emitted.info, "DP").expect("INFO DP");
    kv("live_emit_info_dp", live_dp.to_string());
    kv("live_emit_qual", format!("{:?}", emitted.quality));
    assert_eq!(live_dp, RUST_FALLBACK_DP);
    if let Some(sample) = emitted.samples.first() {
        kv(
            "live_emit_gt_ad_pl_gq_dp",
            format!(
                "GT={:?} AD={:?} PL={:?} GQ={:?} FORMAT_DP={:?}",
                sample.gt.as_ref().map(|g| &g.alleles),
                sample.ad,
                sample.pl,
                sample.gq,
                sample.dp
            ),
        );
        assert_eq!(sample.ad.as_deref(), Some(&[88u32, 22][..]));
    }

    kv(
        "merged_handled_locs_behavior",
        "Call inserts loc 29455649; if !merged_site_handled is false; SiteScore with_annotation_likelihoods never runs; genotyped_calls_at_loc=1",
    );
    kv(
        "ordering",
        "1) merge constructs empty annotation while subset n=123 is in scope; 2) merged_handled_locs.insert(loc); 3) SiteScore skipped; 4) emit empty→region-wide 230",
    );
    kv(
        "candidate_a",
        "TRUE: merge fails to propagate in-scope subset onto GenotypedSiteCall.annotation_likelihoods",
    );
    kv(
        "candidate_b",
        "FALSE: empty is the generic emit-fallback default; no later merge-specific annotation path exists",
    );
    kv(
        "candidate_c",
        "TRUE: merge path has no prepareReadAlleleLikelihoodsForAnnotation analogue",
    );
    kv(
        "candidate_d",
        "PARTIAL: merged_handled_locs skips SiteScore attach (correct Java one-loc merged VC), but does not wipe a populated field",
    );
    kv(
        "candidate_e",
        "TRUE as interaction: missing propagate/prepare AND handled-locs skip SiteScore recovery; first divergent op is Vec::new()",
    );
    kv(
        "smallest_semantic_gap",
        "Java unified loc loop reuses post-retainEvidence AlleleLikelihoods for annotation; Rust merge never copies local subset onto the Call, then skips the only attach helper",
    );
    kv(
        "pl_plus_minus_1",
        "Java 570,0,3517 vs Rust 570,0,3518; independent of empty annotation (PL is written before Vec::new())",
    );
    kv(
        "first_divergent_operation",
        "try_genotype_colocated_snp_indel_merge: annotation_likelihoods = Vec::new() while subset (n=123) is in local scope",
    );
    kv(
        "classification",
        "ANNOTATION_LIKELIHOOD_LIFECYCLE_DIVERGENCE",
    );
    kv("next_arrow", "6R.245");
}
