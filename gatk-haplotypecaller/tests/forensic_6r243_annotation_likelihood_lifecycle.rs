//! 6R.243 located the empty merge constructor at `20:29455649 T/TGTTTG`.
//! 6R.246 replaced it with `annotation_likelihoods: subset.into_owned()`.
//! This file locks that closed contract. Do not retune PL.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r243_annotation_likelihood_lifecycle -- --nocapture --test-threads=1
//! HOLDOUT_6R243=1 cargo test -p gatk-haplotypecaller --test holdout_6r243_annotation_likelihood_lifecycle -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{
    build_per_haplotype_variation_events, merged_alleles_for_genotyping,
    merged_site_uses_joint_gls, variation_events_at_position_from_cache,
};
use gatk_haplotypecaller::genotyping::best_pl_index;
use gatk_haplotypecaller::hc_genotyping_engine::{
    java_alignment_read_overlaps_interval, java_emit_af_decision, take_colocated_merge_numerics,
    DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use gatk_haplotypecaller::read_threading_assembler::DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH;
use gatk_haplotypecaller::variant_site_hc_annotations::coverage_evidence_count;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, region_likelihoods_to_rows,
    traverse_assembly_region_walker, try_emit_call_region_variants, AssemblyRegionCallDisposition,
    CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams, RegionReadLikelihood,
    WalkerTraversalConfig, DEFAULT_STAND_CALL_CONF,
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
const MARGIN: i32 = DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN;
const JAVA_INFO_DP: i32 = 123;
const RUST_FALLBACK_DP: i32 = 230;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R243\t{key}\t{}", value.as_ref());
}

fn unique_indices(likelihoods: &[RegionReadLikelihood]) -> BTreeSet<usize> {
    likelihoods.iter().map(|c| c.read_index.get()).collect()
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

fn src_contains(rel: &str, needle: &str) -> bool {
    let text = fs::read_to_string(repo_root().join(rel)).unwrap_or_default();
    text.contains(needle)
}

#[test]
fn forensic_6r243_source_lifecycle_is_merge_empty_attach() {
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
    let merge_at = assign
        .find("fn try_genotype_colocated_snp_indel_merge")
        .expect("merge fn");
    let merge_rest = &assign[merge_at..];
    let merge_end = merge_rest[1..]
        .find("\nfn ")
        .map(|i| i + 1)
        .unwrap_or(merge_rest.len());
    let merge = &merge_rest[..merge_end];
    kv(
        "merge_annotation_field",
        "annotation_likelihoods: subset.into_owned()",
    );
    assert!(
        merge.contains("annotation_likelihoods: subset.into_owned()"),
        "6R.246 attaches retainEvidence; colocated merge does not construct an empty annotation object"
    );
    assert!(
        !merge.contains("annotation_likelihoods: Vec::new()"),
        "empty annotation_likelihoods is not the live merge constructor"
    );
    assert!(
        !merge.contains("with_annotation_likelihoods"),
        "merge path writes the struct field; it does not call with_annotation_likelihoods"
    );
    assert!(
        assign.contains("merged_handled_locs.insert(loc)"),
        "successful merge skips the per-event SiteScore loc-loop"
    );

    let pipeline = fs::read_to_string(
        repo_root().join("gatk-haplotypecaller/src/hc_genotyping_engine/genotype_site_pipeline.rs"),
    )
    .unwrap();
    assert!(
        pipeline.contains("let annotation = annotation_likelihoods_from_stored_haplotypes"),
        "SiteScore loc-loop still attaches stored-hap retainEvidence"
    );
    assert!(
        src_contains(
            "gatk-haplotypecaller/src/hc_genotyping_engine/genotype_site_early_template.rs",
            "fn annotation_likelihoods_from_stored_haplotypes",
        ),
        "loc-loop helper exists but is not on the merge constructor"
    );

    let emit = fs::read_to_string(repo_root().join("gatk-haplotypecaller/src/region_vcf_emit.rs"))
        .unwrap();
    assert!(
        emit.contains("if call.annotation_likelihoods.is_empty()")
            && emit.contains("fn annotation_likelihoods_for_call"),
        "empty annotation_likelihoods selects region-wide fallback"
    );

    let coverage = fs::read_to_string(
        repo_root().join("gatk-haplotypecaller/src/variant_site_hc_annotations.rs"),
    )
    .unwrap();
    assert!(
        coverage.contains("6R.174: INFO DP is Java `Coverage.evidenceCount`"),
        "Coverage formula is not this arrow"
    );
    kv(
        "java_prepare_ann",
        "contamination off: readAlleleLikelihoodsForAnnotations = readAlleleLikelihoodsForGenotyping; addEvidence(overlappingFilteredReads,0)",
    );
    kv(
        "java_coverage",
        "Coverage.annotate = likelihoods.evidenceCount() of that annotation object",
    );
}

#[test]
fn forensic_6r243_live_empty_is_merge_skip_not_retain_zero() {
    if std::env::var("FORENSIC_6R243_LIVE").ok().as_deref() != Some("1") {
        eprintln!("skip: pre-6R.246 live empty-attach counterfactual");
        return;
    }
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
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
    kv(
        "live_active_region",
        format!(
            "{}:{}-{} extended={}-{} reads={}",
            region.contig,
            region.start.get(),
            region.end.get(),
            region.extended_start.get(),
            region.extended_end.get(),
            region.reads.len()
        ),
    );

    let outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let hap_n = outcome.assembly.haplotypes.len();
    let alt_hap_n = outcome
        .assembly
        .haplotypes
        .iter()
        .filter(|h| !h.is_reference)
        .count();
    kv("hap_n", hap_n.to_string());
    kv("alt_hap_n", alt_hap_n.to_string());

    let events_at: Vec<_> = outcome
        .assembly
        .variation_events()
        .iter()
        .filter(|e| e.start_1based.get() == TARGET)
        .cloned()
        .collect();
    kv(
        "eventmap_at_loc",
        format!(
            "n={}\t{}",
            events_at.len(),
            events_at
                .iter()
                .map(|e| format!("{}/{}", e.ref_allele, e.alt_allele))
                .collect::<Vec<_>>()
                .join(",")
        ),
    );
    assert_eq!(events_at.len(), 2);
    assert!(events_at.iter().any(|e| e.alt_allele == TARGET_ALT));
    assert!(events_at.iter().any(|e| e.alt_allele == "TTTG"));

    let merged = merged_alleles_for_genotyping(&events_at, TARGET).expect("merged alleles");
    kv(
        "merged_alleles",
        format!("REF={} ALTS={}", merged.0, merged.1.join(",")),
    );
    assert!(merged_site_uses_joint_gls(
        &events_at, TARGET, &merged.0, &merged.1
    ));
    kv("merged_site_uses_joint_gls", "true");

    let call = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based.get() == TARGET
                && c.event.ref_allele == TARGET_REF
                && c.event.alt_allele == TARGET_ALT
        })
        .expect("genotyped T/TGTTTG");
    let pl = call.genotype.format.pl_as_i32();
    kv(
        "live_genotype",
        format!(
            "GT_idx={} AD={:?} PL={:?} GQ={} extra_alt={} post_merge={} ann_n={}",
            best_pl_index(&call.genotype.format.pl),
            call.genotype.format.ad_as_i32(),
            pl,
            call.genotype.format.gq.as_i32(),
            call.extra_alt_alleles.len(),
            call.post_merge_unused_alt_subset,
            call.annotation_likelihoods.len()
        ),
    );
    assert_eq!(best_pl_index(&call.genotype.format.pl), 1);
    assert_eq!(call.genotype.format.ad_as_i32(), vec![88, 22]);
    assert_eq!(&pl[..2], &[570, 0]);
    assert!(
        (pl[2] - 3517).abs() <= 1,
        "PL ±1 is representation; got {}",
        pl[2]
    );
    assert!(
        call.post_merge_unused_alt_subset,
        "site is the colocated unused-ALT subset product"
    );
    assert!(
        call.extra_alt_alleles.is_empty(),
        "emitted VC is biallelic after unused-ALT subset"
    );
    assert!(
        call.annotation_likelihoods.is_empty(),
        "production annotation_likelihoods is empty on the merge path"
    );

    kv(
        "early_template_chr20_indel",
        "P12 cluster/gap/sparse early-templates are contig-2 coordinate classes; this is chr20 T/TGTTTG",
    );
    let early = fs::read_to_string(
        root.join("gatk-haplotypecaller/src/hc_genotyping_engine/genotype_site_early_template.rs"),
    )
    .unwrap();
    assert!(
        early.contains("if is_cluster_downstream_snp(&event)")
            && early.contains("if is_p12_phase_e_gap_het_event(&event)")
            && early.contains("if is_cluster_tg_snp(&event)"),
        "early-template arms exist but are skipped because merge handles loc first"
    );
    assert_eq!(call.event.contig, "20");
    assert_eq!(call.event.ref_allele, TARGET_REF);
    assert_eq!(call.event.alt_allele, TARGET_ALT);

    let snaps = take_colocated_merge_numerics();
    let snap = snaps
        .iter()
        .find(|s| s.loc == TARGET)
        .expect("colocated merge snapshot at target");
    kv(
        "merge_snapshot",
        format!(
            "long_ref={}\talts={:?}\tn_reads={}\tn_pairhmm={}\tn_overlap={}\tn_overlap_qnames={}\tpool_sizes={:?}\tmerged_pl={:?}\tsubset_pl={:?}\tsubset_ad_remarg={:?}\tn_haps={}\thap_events={:?}",
            snap.long_ref,
            snap.alts,
            snap.n_reads,
            snap.n_pairhmm_reads,
            snap.n_overlap_before_qname_dedupe,
            snap.n_overlap_unique_qnames,
            snap.pool_sizes,
            snap.merged_pl,
            snap.subset_pl,
            snap.subset_ad_remarginalized,
            snap.n_haps,
            snap.hap_event_signatures_at_loc
        ),
    );
    assert_eq!(
        snap.n_overlap_before_qname_dedupe, 123,
        "merge retainEvidence unique evidence is Java INFO DP 123"
    );
    assert_eq!(snap.n_reads, 123);
    assert_eq!(
        snap.n_pairhmm_reads, 230,
        "stored hap / PairHMM unique is the region-wide fallback 230"
    );
    assert_eq!(snap.n_overlap_unique_qnames, 123);
    kv(
        "merge_overlap_vs_java_info_dp",
        format!(
            "overlap={} java_INFO_DP={} pairhmm={}",
            snap.n_overlap_before_qname_dedupe, JAVA_INFO_DP, snap.n_pairhmm_reads
        ),
    );

    let reads = &outcome.genotyping_reads;
    let stored_unique = unique_indices(&outcome.read_likelihoods);
    kv("stored_hap_unique_n", stored_unique.len().to_string());
    let retain: BTreeSet<usize> = stored_unique
        .iter()
        .copied()
        .filter(|&i| {
            reads
                .get(i)
                .is_some_and(|r| java_alignment_read_overlaps_interval(r, TARGET, TARGET, MARGIN))
        })
        .collect();
    kv(
        "diagnostic_loc_loop_retainEvidence_n",
        retain.len().to_string(),
    );
    assert_eq!(
        retain.len(),
        123,
        "diagnostic loc-loop retainEvidence ±2 is the same 123-read object merge already built"
    );
    let mut retain_ids: Vec<(String, u16)> = retain
        .iter()
        .filter_map(|&i| {
            reads
                .get(i)
                .map(|r| (String::from_utf8_lossy(r.qname()).into_owned(), r.flags()))
        })
        .collect();
    retain_ids.sort();
    kv(
        "retain_head",
        format!("{:?}", &retain_ids[..5.min(retain_ids.len())]),
    );
    kv(
        "retain_tail",
        format!("{:?}", &retain_ids[retain_ids.len().saturating_sub(3)..]),
    );

    let region_wide =
        coverage_evidence_count(reads, &outcome.read_likelihoods, TARGET, TARGET, MARGIN);
    kv(
        "region_wide_coverage_evidence_count",
        region_wide.to_string(),
    );
    assert_eq!(
        region_wide, RUST_FALLBACK_DP,
        "empty attach → region-wide Coverage 230"
    );

    let attached_cov =
        coverage_evidence_count(reads, &call.annotation_likelihoods, TARGET, TARGET, MARGIN);
    kv("attached_coverage_evidence_count", attached_cov.to_string());
    assert_eq!(attached_cov, 0);

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
    assert!(!af.site_is_monomorphic);

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
    assert_eq!(live_dp, RUST_FALLBACK_DP);

    let (full_ref, full_pad) = outcome.assembly.event_map_reference();
    let hap_events = build_per_haplotype_variation_events(
        &outcome.assembly.haplotypes,
        full_ref,
        full_pad,
        outcome.assembly.max_mnp_distance(),
        &region.contig,
    );
    let cache_at = variation_events_at_position_from_cache(&hap_events, TARGET, true);
    kv("hap_eventmap_at_loc_n", cache_at.len().to_string());
    kv(
        "rows_n",
        region_likelihoods_to_rows(&outcome.read_likelihoods, hap_n)
            .len()
            .to_string(),
    );

    kv(
        "pl_plus_minus_1",
        "Java 570,0,3517 vs Rust 570,0,3518; both GT 0/1; does not empty annotation_likelihoods",
    );
    kv(
        "first_divergent_operation",
        "try_genotype_colocated_snp_indel_merge constructs GenotypedSiteCall with annotation_likelihoods=Vec::new(); merged_handled_locs skips SiteScore annotation_likelihoods_from_stored_haplotypes; annotation_likelihoods_for_call then selects region-wide fallback",
    );
    kv(
        "classification",
        "ANNOTATION_LIKELIHOOD_LIFECYCLE_DIVERGENCE",
    );
    assert_ne!(
        snap.n_overlap_before_qname_dedupe as i32, RUST_FALLBACK_DP,
        "merge retainEvidence is not the region-wide fallback object"
    );
}
