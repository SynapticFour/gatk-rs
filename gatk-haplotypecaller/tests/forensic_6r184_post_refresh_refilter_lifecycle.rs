//! 6R.184: `read_likelihoods` lifecycle at `2:92307333 T/G`.
//! After 6R.185 the last P12 refresh is followed by the existing
//! Java-order normalize+filter (stored n=2).
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! PRODUCTION CHANGE: NONE.
//!
//! Does not investigate annotation formulas, INFO, or FORMAT construction.
//! Does not attach any candidate object.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r184_post_refresh_refilter_lifecycle -- --nocapture --test-threads=1
//! HOLDOUT_6R184=1 cargo test -p gatk-haplotypecaller --test holdout_6r184_post_refresh_refilter -- --nocapture --test-threads=1
//! ```

use gatk_haplotypecaller::event_map::{build_per_haplotype_variation_events, VariationEvent};
use gatk_haplotypecaller::hc_allele_mapping::create_allele_mapper_with_events;
use gatk_haplotypecaller::hc_genotyping_engine::java_alignment_read_overlaps_interval;
use gatk_haplotypecaller::{
    begin_likelihood_pipeline_observe, begin_poorly_modeled_observe, call_disposition,
    flatten_assembly_regions, marginalize_rows_to_biallelic_alleles, region_likelihoods_to_rows,
    take_likelihood_pipeline_cells, take_likelihood_pipeline_snaps, take_poorly_modeled_observe,
    traverse_assembly_region_walker, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, LikelihoodPipelineCell, LikelihoodPipelineSnap, ReadFilterParams,
    RegionReadLikelihood, WalkerTraversalConfig,
};
use rust_htslib::bam::Record;
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "2:92307200-92307550";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_307_333;
const MERGED_REF: &str = "T";
const MERGED_ALT: &str = "G";
const CLOSED_INDEL: u64 = 92_307_324;
const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const MARGIN: i32 = 2;
const P12_CLUSTER_UPSTREAM_START: u64 = 92_305_716;
const P12_CLUSTER_UPSTREAM_END: u64 = 92_305_728;
const P12_CLUSTER_TTC_START: u64 = 92_307_324;
const P12_CLUSTER_AC_SNP_START: u64 = 92_307_383;
const P12_CLUSTER_UPSTREAM_MARGINAL_READ_KEEP_LOG10: f64 = -13.5;
const LOG10_GLOBAL_READ_MISMATCHING_RATE: f64 = -4.5;
const JAVA_MQ44_QNAME: &str = "H06JUADXX130110:1:1101:10052:88682";
const JAVA_SURVIVOR2_QNAME: &str = "H06HDADXX130110:1:1101:10061:17286";

#[derive(Clone, Copy)]
#[allow(dead_code)]
struct JavaOrigRead {
    qname: &'static str,
    flag: u16,
    mapq: u8, // identity dump / Java pin
    start: u64,
    orig_len: usize, // Java threshold uses qualified length; pin ≥51 → −8
    norm_max_ll: f64,
    java_keep: bool,
}

/// Live Java `genotype-emit-at-loc` 92307333 (`/tmp/gatk-rs-6r182/genotype-emit-92307333.tsv`).
const JAVA_ORIG: &[JavaOrigRead] = &[
    JavaOrigRead {
        qname: "H06HDADXX130110:1:1101:10061:17286",
        flag: 163,
        mapq: 40,
        start: 92_307_248,
        orig_len: 120,
        norm_max_ll: -13.056375503540,
        java_keep: false,
    },
    JavaOrigRead {
        qname: "H06JUADXX130110:1:1101:10011:51168",
        flag: 81,
        mapq: 22,
        start: 92_307_248,
        orig_len: 190,
        norm_max_ll: -16.238531112671,
        java_keep: false,
    },
    JavaOrigRead {
        qname: JAVA_MQ44_QNAME,
        flag: 163,
        mapq: 50,
        start: 92_307_260,
        orig_len: 179,
        norm_max_ll: -19.377605438232,
        java_keep: false,
    },
    JavaOrigRead {
        qname: "H06HDADXX130110:1:1101:10073:74407",
        flag: 163,
        mapq: 41,
        start: 92_307_272,
        orig_len: 167,
        norm_max_ll: -10.475610733032,
        java_keep: false,
    },
    JavaOrigRead {
        qname: "H06HDADXX130110:1:1101:10073:74407",
        flag: 83,
        mapq: 41,
        start: 92_307_272,
        orig_len: 167,
        norm_max_ll: -10.475610733032,
        java_keep: false,
    },
    JavaOrigRead {
        qname: JAVA_MQ44_QNAME,
        flag: 83,
        mapq: 44,
        start: 92_307_292,
        orig_len: 147,
        norm_max_ll: -5.817768096924,
        java_keep: true,
    },
    JavaOrigRead {
        qname: JAVA_SURVIVOR2_QNAME,
        flag: 83,
        mapq: 40,
        start: 92_307_338,
        orig_len: 100,
        norm_max_ll: -2.494571685791,
        java_keep: true,
    },
    JavaOrigRead {
        qname: "H06HDADXX130110:2:1101:10097:72839",
        flag: 147,
        mapq: 23,
        start: 92_307_367,
        orig_len: 72,
        norm_max_ll: -18.344545364380,
        java_keep: false,
    },
];

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R184\t{key}\t{}", value.as_ref());
}

fn qname(rec: &Record) -> String {
    String::from_utf8_lossy(rec.qname()).into_owned()
}

fn java_log10_min_true_likelihood(qualified_read_len: usize) -> f64 {
    let max_errors = (qualified_read_len as f64 * 0.02).ceil().min(2.0);
    max_errors * -4.0
}

fn java_keep(max_ll: f64, qualified_read_len: usize) -> bool {
    max_ll.is_finite() && !(max_ll < java_log10_min_true_likelihood(qualified_read_len))
}

fn unique_likelihood_indices(likelihoods: &[RegionReadLikelihood]) -> BTreeSet<usize> {
    likelihoods.iter().map(|c| c.read_index.get()).collect()
}

fn best_ll_by_index(likelihoods: &[RegionReadLikelihood]) -> BTreeMap<usize, f64> {
    let mut best = BTreeMap::new();
    for cell in likelihoods {
        let idx = cell.read_index.get();
        let ll = cell.log10_likelihood;
        best.entry(idx)
            .and_modify(|m| {
                if ll > *m {
                    *m = ll;
                }
            })
            .or_insert(ll);
    }
    best
}

fn match_java_orig(rec: &Record) -> Option<&'static JavaOrigRead> {
    let qn = qname(rec);
    let flag = rec.flags();
    let start = (rec.pos() + 1) as u64;
    JAVA_ORIG
        .iter()
        .find(|j| j.qname == qn && j.flag == flag && j.start == start)
        .or_else(|| JAVA_ORIG.iter().find(|j| j.qname == qn && j.flag == flag))
}

fn p12_cluster_span(active_start: u64, active_end: u64) -> bool {
    active_end >= P12_CLUSTER_TTC_START.saturating_sub(50)
        && active_start <= P12_CLUSTER_AC_SNP_START.saturating_add(50)
}

fn unique_from_cells(cells: &[LikelihoodPipelineCell], seq: u32) -> BTreeSet<usize> {
    cells
        .iter()
        .filter(|c| c.seq == seq)
        .map(|c| c.read_index)
        .collect()
}

/// Diagnostic replica of `normalize_region_read_likelihoods` (engine.rs).
/// Integration tests cannot call the private production fn without exporting it.
fn diagnostic_normalize_all_columns(ll: &mut [RegionReadLikelihood], n_haps: usize) {
    if ll.is_empty() {
        return;
    }
    let eligible: Vec<usize> = (0..n_haps).collect();
    let max_read = ll.iter().map(|e| e.read_index.get()).max().unwrap_or(0);
    let max_hap = ll
        .iter()
        .map(|e| e.haplotype_index.get())
        .max()
        .unwrap_or(0);
    let mut eligible_mask = vec![false; max_hap + 1];
    for &i in &eligible {
        if i < eligible_mask.len() {
            eligible_mask[i] = true;
        }
    }
    let mut best = vec![f64::NEG_INFINITY; max_read + 1];
    for entry in ll.iter() {
        if eligible_mask
            .get(entry.haplotype_index.get())
            .copied()
            .unwrap_or(false)
            && entry.log10_likelihood > best[entry.read_index.get()]
        {
            best[entry.read_index.get()] = entry.log10_likelihood;
        }
    }
    for entry in ll.iter() {
        let i = entry.read_index.get();
        if !best[i].is_finite() && entry.log10_likelihood > best[i] {
            best[i] = entry.log10_likelihood;
        }
    }
    for entry in ll.iter_mut() {
        let b = best[entry.read_index.get()];
        if !b.is_finite() {
            continue;
        }
        let floor = b + LOG10_GLOBAL_READ_MISMATCHING_RATE;
        if entry.log10_likelihood < floor {
            entry.log10_likelihood = floor;
        }
    }
}

fn extras_keep(best_ll: f64, rec: &Record, active: (u64, u64)) -> bool {
    if !best_ll.is_finite() {
        return false;
    }
    let qual_len = rec.qual().len().max(1);
    let p12_up = qual_len >= 100
        && best_ll >= P12_CLUSTER_UPSTREAM_MARGINAL_READ_KEEP_LOG10
        && java_alignment_read_overlaps_interval(
            rec,
            P12_CLUSTER_UPSTREAM_START,
            P12_CLUSTER_UPSTREAM_END,
            0,
        );
    if p12_up {
        return true;
    }
    use gatk_haplotypecaller::hc_genotyping_engine::soft_unclipped_read_overlaps_interval;
    let (active_start, active_end) = active;
    let soft_overlaps = soft_unclipped_read_overlaps_interval(rec, active_start, active_end, 2);
    let align_overlaps = java_alignment_read_overlaps_interval(rec, active_start, active_end, 2);
    if !soft_overlaps {
        return false;
    }
    if !align_overlaps && qual_len >= 20 && best_ll >= P12_CLUSTER_UPSTREAM_MARGINAL_READ_KEEP_LOG10
    {
        return true;
    }
    qual_len >= 50 && best_ll >= P12_CLUSTER_UPSTREAM_MARGINAL_READ_KEEP_LOG10 && !align_overlaps
}

/// Diagnostic replica of `filter_normalized_region_read_likelihoods` →
/// `filter_poorly_modeled_region_read_likelihoods` (engine.rs).
fn diagnostic_filter_existing(
    ll: &[RegionReadLikelihood],
    reads: &[impl std::borrow::Borrow<Record>],
    active: (u64, u64),
) -> (
    Vec<RegionReadLikelihood>,
    BTreeMap<usize, (f64, f64, bool, bool)>,
) {
    let best = best_ll_by_index(ll);
    let mut keep = BTreeSet::new();
    let mut dump = BTreeMap::new();
    for (&idx, &max_ll) in &best {
        let Some(slot) = reads.get(idx) else {
            continue;
        };
        let rec = slot.borrow();
        let qlen = rec.qual().len().max(1);
        let thresh = java_log10_min_true_likelihood(qlen);
        let java_equiv = java_keep(max_ll, qlen);
        let extra = extras_keep(max_ll, rec, active);
        let retain = java_equiv || extra;
        dump.insert(idx, (max_ll, thresh, java_equiv, extra));
        if retain {
            keep.insert(idx);
        }
    }
    let filtered: Vec<RegionReadLikelihood> = ll
        .iter()
        .filter(|c| keep.contains(&c.read_index.get()))
        .cloned()
        .collect();
    (filtered, dump)
}

#[test]
fn forensic_6r184_source_contracts_no_production_change() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    let engine = include_str!("../src/engine.rs");
    let early = include_str!("../src/hc_genotyping_engine/genotype_site_early_template.rs");
    let fin = include_str!("../src/hc_genotyping_engine/genotype_finalize.rs");
    let pipe = include_str!("../src/hc_genotyping_engine/genotype_site_pipeline.rs");
    let emit = include_str!("../src/region_vcf_emit.rs");
    let ann = include_str!("../src/variant_site_hc_annotations.rs");

    assert!(
        engine.contains("fn normalize_region_read_likelihoods")
            && engine.contains("fn filter_normalized_region_read_likelihoods")
            && engine.contains("fn filter_poorly_modeled_region_read_likelihoods"),
        "existing Java-order normalize/filter still present"
    );
    assert!(
        engine.contains("let floor = b + LOG10_GLOBAL_READ_MISMATCHING_RATE;")
            && engine.contains("if entry.log10_likelihood < floor"),
        "normalize still floors losing cells; max_ll unchanged"
    );
    assert!(
        engine.contains("let ll_normalize = !args.is_strict_java();"),
        "first PairHMM still skips normalize+filter under strict_java"
    );
    assert!(
        engine.contains("fn apply_java_order_normalize_and_filter")
            && engine.contains("normalize_region_read_likelihoods(read_likelihoods, &norm_haps)")
            && engine.matches("apply_java_order_normalize_and_filter(").count() >= 2,
        "Java-order normalize+filter helper still used after allele filter and after last P12 refresh"
    );

    let outcome_move = engine
        .find("Ok(Some(CallRegionOutcome {")
        .expect("outcome move");
    let last_refresh = engine[..outcome_move]
        .rfind("read_likelihoods = ll;")
        .expect("at least one refresh assignment before outcome");
    let after_last = &engine[last_refresh..outcome_move];
    assert!(
        after_last.contains("apply_java_order_normalize_and_filter"),
        "6R.185 restores Java-order normalize+filter after the last P12 refresh"
    );
    assert!(
        after_last.contains("likelihoods: &read_likelihoods"),
        "after last assignment, likelihoods are borrowed for genotyping"
    );

    assert!(
        early.contains("if is_cluster_tg_snp(&event)")
            && early.contains("finish_strict_java_shaped_site_call"),
        "cluster-TG still selected in SiteEarlyTemplate::try_shaped"
    );
    let fin_fn = fin
        .split("fn finish_strict_java_shaped_site_call")
        .nth(1)
        .expect("finish fn");
    let fin_body = fin_fn
        .split("fn finalize_strict_java_variation_genotype")
        .next()
        .unwrap_or(fin_fn);
    assert!(
        !fin_body.contains("annotation_likelihoods") && !fin_body.contains("read_likelihoods ="),
        "finish_strict_java_shaped_site_call still does not write likelihoods"
    );
    assert!(
        pipe.contains("SiteEarlyTemplate::try_shaped") && pipe.contains("return Ok(Some(call));"),
        "try_shaped still returns without reconstructing likelihoods"
    );
    assert!(
        engine.contains("change_evidence_to_best_haplotype")
            && include_str!("../src/read_realignment.rs")
                .contains("pub fn change_evidence_to_best_haplotype")
            && include_str!("../src/read_realignment.rs").contains("likelihoods"),
        "change_evidence remains identity on the likelihood matrix"
    );
    assert!(
        !ann.contains("92307333") && !emit.contains("92307333") && !engine.contains("92307333"),
        "no locus-specific 6R.184 production patch"
    );
    assert!(
        !engine.contains(JAVA_MQ44_QNAME) && !ann.contains(JAVA_MQ44_QNAME),
        "no QNAME exception"
    );
}

#[test]
fn forensic_6r184_lifecycle_refresh_is_last_membership_write() {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing P12 ref/BAM");
        return;
    }

    kv("java_pin", JAVA_PIN);
    kv("target", "2:92307333 T/G");
    kv("production_change", "NONE");
    kv(
        "scope",
        "read_likelihoods lifecycle + diagnostic post-refresh re-filter; not annotation attach",
    );

    let dict = gatk_core::reference::SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs =
        gatk_core::reference::parse_intervals_cli_string(&dict, INTERVAL).expect("interval");
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
    let covering = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= TARGET
                && r.end.get() >= TARGET
        })
        .expect("covering");
    kv(
        "active_region",
        format!(
            "{}:{}-{}",
            covering.contig,
            covering.start.get(),
            covering.end.get()
        ),
    );
    let in_cluster = p12_cluster_span(covering.start.get(), covering.end.get());
    kv("strict_java_p12_cluster_span", in_cluster.to_string());
    assert!(in_cluster);

    begin_poorly_modeled_observe();
    begin_likelihood_pipeline_observe();
    let outcome = HaplotypeCallerEngine::call_region(
        covering,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("outcome");
    let observed = take_poorly_modeled_observe();
    let snaps: Vec<LikelihoodPipelineSnap> = take_likelihood_pipeline_snaps();
    let cells = take_likelihood_pipeline_cells();

    assert!(
        !snaps.is_empty(),
        "pipeline snaps must record every kernel/stage"
    );
    for s in &snaps {
        let uniq = unique_from_cells(&cells, s.seq);
        kv(
            "LIFE_SNAP",
            format!(
                "seq={} stage={} n_reads={} n_haps={} n_ll={} unique_reads={}",
                s.seq,
                s.stage,
                s.n_reads,
                s.n_haps,
                s.n_ll_entries,
                uniq.len()
            ),
        );
    }

    let last_pass = observed.iter().map(|r| r.pass).max().unwrap_or(0);
    let last_rows: Vec<_> = observed.iter().filter(|r| r.pass == last_pass).collect();
    let java_equiv_keep_n = last_rows.iter().filter(|r| r.java_equiv_keep).count();
    let rust_keep_n = last_rows.iter().filter(|r| r.rust_keep).count();
    let extra_n = last_rows.iter().filter(|r| r.extra_retain).count();
    kv("pre_refresh_filter_pass", last_pass.to_string());
    kv(
        "pre_refresh_filter",
        format!("java_equiv_keep={java_equiv_keep_n} rust_keep={rust_keep_n} extra={extra_n}"),
    );
    assert_eq!(java_equiv_keep_n, 2);
    assert_eq!(rust_keep_n, 2);
    assert_eq!(extra_n, 0);

    let normalize_snaps: Vec<_> = snaps.iter().filter(|s| s.stage == "normalize").collect();
    assert!(
        normalize_snaps.len() >= 2,
        "Java-order normalize runs after allele filter and again after the last P12 refresh"
    );
    let last = snaps.last().expect("last snap");
    kv("last_snap_stage", last.stage);
    assert_eq!(
        last.stage, "filter",
        "last captured likelihood write is the restored poorly-modeled filter"
    );
    let last_unique = unique_from_cells(&cells, last.seq);
    kv("last_filter_unique_n", last_unique.len().to_string());
    kv("last_filter_n_haps", last.n_haps.to_string());

    let stored_unique = unique_likelihood_indices(&outcome.read_likelihoods);
    kv("stored_unique_n", stored_unique.len().to_string());
    kv(
        "stored_hap_n",
        outcome.assembly.haplotypes.len().to_string(),
    );
    assert_eq!(stored_unique.len(), 2);
    assert_eq!(
        last_unique.len(),
        2,
        "last filter unique evidence is the stored n=2"
    );
    assert_eq!(
        last.n_haps,
        outcome.assembly.haplotypes.len(),
        "last filter haplotype population is the stored haplotype set"
    );

    let refresh_after_normalize = snaps
        .iter()
        .filter(|s| s.seq > normalize_snaps[0].seq && s.stage == "refresh")
        .count();
    kv(
        "refresh_snaps_after_normalize",
        refresh_after_normalize.to_string(),
    );
    assert!(
        refresh_after_normalize >= 1,
        "at least one unfiltered refresh after the Java-order filter"
    );

    let call = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(TARGET)
                && c.event.ref_allele == MERGED_REF
                && c.event.alt_allele == MERGED_ALT
        })
        .expect("T/G");
    assert!(
        call.annotation_likelihoods.is_empty(),
        "cluster-TG still does not construct/subset a likelihood object"
    );
    let indel = outcome.genotyped_calls.iter().find(|c| {
        c.event.start_1based == GenomePosition::new_1based(CLOSED_INDEL)
            && c.event.ref_allele == "TTC"
    });
    assert!(indel.is_some());

    let active = (covering.start.get(), covering.end.get());
    let n_haps = outcome.assembly.haplotypes.len();
    let stored_best = best_ll_by_index(&outcome.read_likelihoods);

    let (raw_filtered, raw_dump) =
        diagnostic_filter_existing(&outcome.read_likelihoods, &outcome.genotyping_reads, active);
    let raw_unique = unique_likelihood_indices(&raw_filtered);

    let mut normalized = outcome.read_likelihoods.clone();
    diagnostic_normalize_all_columns(&mut normalized, n_haps);
    let norm_best = best_ll_by_index(&normalized);
    let (norm_filtered, norm_dump) =
        diagnostic_filter_existing(&normalized, &outcome.genotyping_reads, active);
    let norm_unique = unique_likelihood_indices(&norm_filtered);

    kv("pre_refresh_filtered_n", "2");
    kv("post_refresh_raw_n", stored_unique.len().to_string());
    kv("post_refresh_refilter_raw_n", raw_unique.len().to_string());
    kv(
        "post_refresh_refilter_after_normalize_n",
        norm_unique.len().to_string(),
    );

    kv(
        "TABLE",
        "READ\tJAVA_BEST_LL\tJAVA_FILTER\tRUST_STORED_BEST\tTHRESH\tRAW_KEEP\tNORM_MAX\tNORM_KEEP\tEXTRA\tMAX_LL_CHANGED",
    );
    for &idx in &stored_unique {
        let rec = outcome.genotyping_reads.get(idx).expect("read");
        let orig = match_java_orig(rec);
        let stored = stored_best.get(&idx).copied().unwrap_or(f64::NEG_INFINITY);
        let after_norm = norm_best.get(&idx).copied().unwrap_or(f64::NEG_INFINITY);
        let (raw_max, thresh, java_equiv, extra) = raw_dump.get(&idx).copied().unwrap_or((
            stored,
            java_log10_min_true_likelihood(rec.qual().len().max(1)),
            false,
            false,
        ));
        let norm_keep = norm_unique.contains(&idx);
        let raw_keep = raw_unique.contains(&idx);
        let max_changed = (stored - after_norm).abs() > 1e-15;
        kv(
            "ROW",
            format!(
                "{} FLAG={} MAPQ={} start={}\t{:.6}\t{}\t{:.6}\t{:.1}\t{}\t{:.6}\t{}\t{}\t{}",
                qname(rec),
                rec.flags(),
                rec.mapq(),
                rec.pos() + 1,
                orig.map(|j| j.norm_max_ll).unwrap_or(f64::NAN),
                orig.map(|j| if j.java_keep { "KEEP" } else { "DROP" })
                    .unwrap_or("?"),
                stored,
                thresh,
                if raw_keep { "KEEP" } else { "DROP" },
                after_norm,
                if norm_keep { "KEEP" } else { "DROP" },
                extra,
                max_changed,
            ),
        );
        assert_eq!(raw_keep, java_equiv || extra);
        assert!(
            !max_changed,
            "normalize floors losers; max_ll must be unchanged for idx={idx}"
        );
        let _ = (raw_max, java_equiv);
    }

    assert_eq!(raw_unique, norm_unique, "normalize does not flip KEEP/DROP");
    assert_eq!(
        raw_unique.len(),
        2,
        "STOP: post-refresh re-filter is not Java n=2"
    );
    assert_eq!(norm_unique.len(), 2);

    let mut keep_ids = Vec::new();
    for &idx in &norm_unique {
        let rec = &outcome.genotyping_reads[idx];
        keep_ids.push((qname(rec), rec.flags(), rec.mapq(), rec.pos() + 1));
        kv(
            "REFILTER_KEEP",
            format!(
                "{} FLAG={} MAPQ={} start={} best_ll={:.6} thresh={:.1}",
                qname(rec),
                rec.flags(),
                rec.mapq(),
                rec.pos() + 1,
                stored_best[&idx],
                java_log10_min_true_likelihood(rec.qual().len().max(1)),
            ),
        );
    }
    keep_ids.sort();
    let has_mq44 = keep_ids
        .iter()
        .any(|(q, f, m, s)| q == JAVA_MQ44_QNAME && *f == 83 && *m == 44 && *s == 92_307_292);
    let has_survivor2 = keep_ids
        .iter()
        .any(|(q, f, m, s)| q == JAVA_SURVIVOR2_QNAME && *f == 83 && *m == 40 && *s == 92_307_338);
    assert!(
        has_mq44,
        "re-filter must keep Java survivor FLAG=83 MAPQ=44"
    );
    assert!(
        has_survivor2,
        "re-filter must keep Java survivor FLAG=83 MAPQ=40 start=92307338"
    );
    assert_eq!(
        extra_n, 0,
        "extras KEEP still does not fire on the stored matrix"
    );
    let extra_on_refilter = norm_dump.values().filter(|(_, _, _, extra)| *extra).count();
    assert_eq!(extra_on_refilter, 0);

    let (full_ref, full_pad) = outcome.assembly.event_map_reference();
    let apply_pad = outcome
        .assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .and_then(|h| h.genome_loc.as_ref().map(|g| g.start.get()))
        .unwrap_or(full_pad);
    let hap_cache = build_per_haplotype_variation_events(
        &outcome.assembly.haplotypes,
        full_ref,
        full_pad,
        outcome.assembly.max_mnp_distance(),
        &covering.contig,
    );
    let mapping = create_allele_mapper_with_events(
        &VariationEvent::from_alleles(&covering.contig, TARGET, MERGED_REF, MERGED_ALT),
        TARGET,
        &outcome.assembly.haplotypes,
        apply_pad,
        outcome.assembly.reference_bases(),
        outcome.assembly.max_mnp_distance(),
        true,
        Some(&hap_cache),
    );
    let rows = region_likelihoods_to_rows(&norm_filtered, n_haps);
    let marg = marginalize_rows_to_biallelic_alleles(
        &rows,
        &mapping.ref_haplotype_indices,
        &mapping.alt_haplotype_indices,
    );
    kv("diagnostic_marginalize_n", marg.len().to_string());
    assert_eq!(marg.len(), 2, "marginalize of re-filtered n=2 stays n=2");

    let retained: Vec<_> = marg
        .iter()
        .filter(|row| {
            outcome
                .genotyping_reads
                .get(row.read_index)
                .is_some_and(|r| java_alignment_read_overlaps_interval(r, TARGET, TARGET, MARGIN))
        })
        .collect();
    kv("diagnostic_retainEvidence_n", retained.len().to_string());
    assert_eq!(retained.len(), 1);
    let rec = &outcome.genotyping_reads[retained[0].read_index];
    let lls = &retained[0].haplotype_log10_likelihoods;
    let ll_ref = lls.first().copied().unwrap_or(f64::NEG_INFINITY);
    let ll_alt = lls.get(1).copied().unwrap_or(f64::NEG_INFINITY);
    let best_allele = if ll_alt > ll_ref { "G" } else { "T" };
    kv(
        "diagnostic_n1",
        format!(
            "{} FLAG={} MAPQ={} start={} best_allele={} ll_ref={:.6} ll_alt={:.6}",
            qname(rec),
            rec.flags(),
            rec.mapq(),
            rec.pos() + 1,
            best_allele,
            ll_ref,
            ll_alt,
        ),
    );
    assert_eq!(qname(rec), JAVA_MQ44_QNAME);
    assert_eq!(rec.flags(), 83);
    assert_eq!(rec.mapq(), 44);
    assert_eq!(rec.pos() + 1, 92_307_292);
    assert_eq!(best_allele, "G");

    kv(
        "first_causal_arrow",
        "filtered n=2 → unfiltered refresh → restored Java-order filter → stored n=2",
    );
    kv(
        "candidate",
        "preserve refresh; restore existing normalize+filter on the last refreshed matrix before it becomes stored evidence",
    );
    kv("classification", "D — WRONG UPSTREAM GENOTYPING SUBSET");
}
