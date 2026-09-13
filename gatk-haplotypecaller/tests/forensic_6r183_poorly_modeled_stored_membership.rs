//! 6R.183: why Rust stored hap likelihoods are n=8 where Java
//! `filterPoorlyModeledEvidence` leaves n=2 at `2:92307333 T/G` (proof-only).
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! PRODUCTION CHANGE: NONE.
//!
//! Does not investigate annotation, INFO, or FORMAT.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r183_poorly_modeled_stored_membership -- --nocapture --test-threads=1
//! HOLDOUT_6R183=1 cargo test -p gatk-haplotypecaller --test holdout_6r183_poorly_modeled -- --nocapture --test-threads=1
//! ```

use gatk_haplotypecaller::hc_genotyping_engine::java_alignment_read_overlaps_interval;
use gatk_haplotypecaller::{
    begin_likelihood_pipeline_observe, begin_poorly_modeled_observe, call_disposition,
    flatten_assembly_regions, take_likelihood_pipeline_snaps, take_poorly_modeled_observe,
    traverse_assembly_region_walker, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
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
const P12_CLUSTER_UPSTREAM_START: u64 = 92_305_716;
const P12_CLUSTER_UPSTREAM_END: u64 = 92_305_728;
const P12_CLUSTER_TTC_START: u64 = 92_307_324;
const P12_CLUSTER_AC_SNP_START: u64 = 92_307_383;
const JAVA_MQ44_QNAME: &str = "H06JUADXX130110:1:1101:10052:88682";
const JAVA_SURVIVOR2_QNAME: &str = "H06HDADXX130110:1:1101:10061:17286";

#[derive(Clone, Copy)]
struct JavaOrigRead {
    qname: &'static str,
    flag: u16,
    mapq: u8,
    start: u64,
    end: u64,
    cigar: &'static str,
    orig_len: usize,
    best_hap: &'static str,
    prim_max_ll: f64,
    norm_max_ll: f64,
    java_keep: bool,
}

/// Live Java `genotype-emit-at-loc` 92307333 (`/tmp/gatk-rs-6r182/genotype-emit-92307333.tsv`).
/// `stored_n_ev=2` is the post-filter haplotype AlleleLikelihoods.
const JAVA_ORIG: &[JavaOrigRead] = &[
    JavaOrigRead {
        qname: "H06HDADXX130110:1:1101:10061:17286",
        flag: 163,
        mapq: 40,
        start: 92_307_248,
        end: 92_307_368,
        cigar: "56H112M1D8M74H",
        orig_len: 120,
        best_hap: "cc52c0e6e338ead1",
        prim_max_ll: -13.056375503540,
        norm_max_ll: -13.056375503540,
        java_keep: false,
    },
    JavaOrigRead {
        qname: "H06JUADXX130110:1:1101:10011:51168",
        flag: 81,
        mapq: 22,
        start: 92_307_248,
        end: 92_307_438,
        cigar: "4H112M1D78M56H",
        orig_len: 190,
        best_hap: "cc52c0e6e338ead1",
        prim_max_ll: -16.238531112671,
        norm_max_ll: -16.238531112671,
        java_keep: false,
    },
    JavaOrigRead {
        qname: JAVA_MQ44_QNAME,
        flag: 163,
        mapq: 50,
        start: 92_307_260,
        end: 92_307_438,
        cigar: "179M71H",
        orig_len: 179,
        best_hap: "d926aa50c2c117aa",
        prim_max_ll: -19.377605438232,
        norm_max_ll: -19.377605438232,
        java_keep: false,
    },
    JavaOrigRead {
        qname: "H06HDADXX130110:1:1101:10073:74407",
        flag: 163,
        mapq: 41,
        start: 92_307_272,
        end: 92_307_438,
        cigar: "167M83H",
        orig_len: 167,
        best_hap: "d926aa50c2c117aa",
        prim_max_ll: -10.475610733032,
        norm_max_ll: -10.475610733032,
        java_keep: false,
    },
    JavaOrigRead {
        qname: "H06HDADXX130110:1:1101:10073:74407",
        flag: 83,
        mapq: 41,
        start: 92_307_272,
        end: 92_307_438,
        cigar: "31H167M52H",
        orig_len: 167,
        best_hap: "d926aa50c2c117aa",
        prim_max_ll: -10.475610733032,
        norm_max_ll: -10.475610733032,
        java_keep: false,
    },
    JavaOrigRead {
        qname: JAVA_MQ44_QNAME,
        flag: 83,
        mapq: 44,
        start: 92_307_292,
        end: 92_307_438,
        cigar: "14H147M89H",
        orig_len: 147,
        best_hap: "d926aa50c2c117aa",
        prim_max_ll: -5.817768096924,
        norm_max_ll: -5.817768096924,
        java_keep: true,
    },
    JavaOrigRead {
        qname: JAVA_SURVIVOR2_QNAME,
        flag: 83,
        mapq: 40,
        start: 92_307_338,
        end: 92_307_438,
        cigar: "1H22M1D78M149H",
        orig_len: 100,
        best_hap: "9188a56b92ef6932/cc52c0e6e338ead1",
        prim_max_ll: -2.494571685791,
        norm_max_ll: -2.494571685791,
        java_keep: true,
    },
    JavaOrigRead {
        qname: "H06HDADXX130110:2:1101:10097:72839",
        flag: 147,
        mapq: 23,
        start: 92_307_367,
        end: 92_307_438,
        cigar: "83H72M95H",
        orig_len: 72,
        best_hap: "9188a56b92ef6932/cc52c0e6e338ead1",
        prim_max_ll: -18.344545364380,
        norm_max_ll: -18.344545364380,
        java_keep: false,
    },
];

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R183\t{key}\t{}", value.as_ref());
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

fn unique_likelihood_indices(
    likelihoods: &[gatk_haplotypecaller::RegionReadLikelihood],
) -> BTreeSet<usize> {
    likelihoods.iter().map(|c| c.read_index.get()).collect()
}

fn best_ll_by_index(
    likelihoods: &[gatk_haplotypecaller::RegionReadLikelihood],
) -> BTreeMap<usize, f64> {
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

fn overlaps_p12_cluster_upstream(rec: &Record) -> bool {
    java_alignment_read_overlaps_interval(
        rec,
        P12_CLUSTER_UPSTREAM_START,
        P12_CLUSTER_UPSTREAM_END,
        0,
    )
}

fn p12_cluster_span(active_start: u64, active_end: u64) -> bool {
    active_end >= P12_CLUSTER_TTC_START.saturating_sub(50)
        && active_start <= P12_CLUSTER_AC_SNP_START.saturating_add(50)
}

#[test]
fn forensic_6r183_java_filter_contract_from_pin_and_dump() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    kv(
        "java_call_chain",
        "PairHMMLikelihoodCalculationEngine.computeReadLikelihoods → AlleleLikelihoods.normalizeLikelihoods → ReadLikelihoodCalculationEngine.filterPoorlyModeledEvidence(dynamic=false) → AlleleLikelihoods.filterPoorlyModeledEvidence(log10MinTrueLikelihood(cap=true))",
    );
    kv(
        "java_predicate",
        "DROP iff maximumLikelihoodOverAllAlleles(sample,i) < log10MinTrueLikelihood(evidence[i]); KEEP iff !(max_ll < thresh). Mutates via removeEvidenceByIndex. Not relative-to-best, not overlap.",
    );
    kv(
        "java_threshold",
        "qualifiedLen=HMM_BASE_QUALITIES.len else GATKRead.getLength(); maxErrors=min(2, ceil(qlen*0.02)); thresh=maxErrors*-4.0. All 8 reads qlen≥51 → thresh=-8.0",
    );
    kv("java_hap_n", "6");
    kv("java_orig_n", "8");
    kv("java_stored_n", "2");
    let keep: Vec<_> = JAVA_ORIG.iter().filter(|r| r.java_keep).collect();
    let drop: Vec<_> = JAVA_ORIG.iter().filter(|r| !r.java_keep).collect();
    assert_eq!(keep.len(), 2);
    assert_eq!(drop.len(), 6);
    assert_eq!(keep[0].qname, JAVA_MQ44_QNAME);
    assert_eq!(keep[0].flag, 83);
    assert_eq!(keep[0].mapq, 44);
    assert_eq!(keep[1].qname, JAVA_SURVIVOR2_QNAME);
    assert_eq!(keep[1].flag, 83);
    assert_eq!(keep[1].mapq, 40);
    assert_eq!(keep[1].start, 92_307_338);
    kv(
        "JAVA_FILTER_SURVIVORS",
        format!(
            "{} FLAG={} MAPQ={} start={} max_ll={:.6}; {} FLAG={} MAPQ={} start={} max_ll={:.6}",
            keep[0].qname,
            keep[0].flag,
            keep[0].mapq,
            keep[0].start,
            keep[0].norm_max_ll,
            keep[1].qname,
            keep[1].flag,
            keep[1].mapq,
            keep[1].start,
            keep[1].norm_max_ll,
        ),
    );
    for r in drop {
        kv(
            "JAVA_FILTER_DROP",
            format!(
                "{} FLAG={} MAPQ={} start={} CIGAR={} orig_len={} best_hap={} prim_max={:.6} norm_max={:.6} thresh=-8 DROP because {:.6} < -8",
                r.qname,
                r.flag,
                r.mapq,
                r.start,
                r.cigar,
                r.orig_len,
                r.best_hap,
                r.prim_max_ll,
                r.norm_max_ll,
                r.norm_max_ll,
            ),
        );
        assert!(r.norm_max_ll < -8.0);
        assert!(!java_keep(r.norm_max_ll, r.orig_len));
    }
    for r in JAVA_ORIG.iter().filter(|r| r.java_keep) {
        assert!(java_keep(r.norm_max_ll, r.orig_len));
        kv(
            "JAVA_FILTER_KEEP",
            format!(
                "{} FLAG={} MAPQ={} start={} CIGAR={} orig_len={} best_hap={} prim_max={:.6} norm_max={:.6} thresh=-8 KEEP because {:.6} >= -8",
                r.qname,
                r.flag,
                r.mapq,
                r.start,
                r.cigar,
                r.orig_len,
                r.best_hap,
                r.prim_max_ll,
                r.norm_max_ll,
                r.norm_max_ll,
            ),
        );
    }
}

#[test]
fn forensic_6r183_source_contracts_no_production_change() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    let engine = include_str!("../src/engine.rs");
    let early = include_str!("../src/hc_genotyping_engine/genotype_site_early_template.rs");
    let emit = include_str!("../src/region_vcf_emit.rs");
    let ann = include_str!("../src/variant_site_hc_annotations.rs");
    assert!(
        engine.contains("fn filter_poorly_modeled_region_read_likelihoods")
            && engine.contains("retain_marginal_p12_cluster_upstream_read")
            && engine.contains("retain_marginal_sparse_softclip_read"),
        "Rust extras KEEP still inside poorly-modeled filter"
    );
    assert!(
        engine.contains("let ll_normalize = !args.is_strict_java();"),
        "strict_java first PairHMM still skips normalize+filter"
    );
    assert!(
        engine.contains("filter_normalized_region_read_likelihoods")
            && engine.contains("refresh_region_read_likelihoods"),
        "strict_java later filter + refresh still present"
    );
    let refresh_false = engine.matches("refresh_region_read_likelihoods").count();
    assert!(
        refresh_false >= 1,
        "refresh helper must exist for post-filter rescore"
    );
    assert!(
        engine.contains("apply_normalize: bool")
            && engine.contains("if !apply_normalize")
            && engine.contains("return ll;"),
        "apply_normalize=false still skips filterPoorlyModeledEvidence"
    );
    assert!(
        early.contains("if is_cluster_tg_snp(&event)"),
        "cluster-TG path unchanged"
    );
    assert!(
        !ann.contains("92307333") && !emit.contains("92307333") && !engine.contains("92307333"),
        "no locus-specific 6R.183 production patch"
    );
}

#[test]
fn forensic_6r183_rust_filter_vs_java_stored_membership() {
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
        "stored haplotype AlleleLikelihoods membership only; not annotation",
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
    let in_cluster_span = p12_cluster_span(covering.start.get(), covering.end.get());
    kv("strict_java_p12_cluster_span", in_cluster_span.to_string());
    assert!(
        in_cluster_span,
        "this ActiveFull is inside P12 cluster span; post-filter refresh applies"
    );

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
    let snaps = take_likelihood_pipeline_snaps();

    kv("observe_filter_rows", observed.len().to_string());
    let passes: BTreeSet<u32> = observed.iter().map(|r| r.pass).collect();
    kv("observe_filter_passes", format!("{passes:?}"));
    for s in &snaps {
        kv(
            "PIPE_SNAP",
            format!(
                "seq={} stage={} n_reads={} n_haps={} n_ll={}",
                s.seq, s.stage, s.n_reads, s.n_haps, s.n_ll_entries
            ),
        );
    }

    let last_pass = observed.iter().map(|r| r.pass).max().unwrap_or(0);
    let last_rows: Vec<_> = observed.iter().filter(|r| r.pass == last_pass).collect();
    kv("last_filter_pass", last_pass.to_string());
    kv("last_filter_n_rows", last_rows.len().to_string());
    let last_java_keep_n = last_rows.iter().filter(|r| r.java_equiv_keep).count();
    let last_rust_keep_n = last_rows.iter().filter(|r| r.rust_keep).count();
    let last_extra_n = last_rows.iter().filter(|r| r.extra_retain).count();
    kv("last_pass_java_equiv_keep_n", last_java_keep_n.to_string());
    kv("last_pass_rust_keep_n", last_rust_keep_n.to_string());
    kv("last_pass_extra_retain_n", last_extra_n.to_string());

    let stored_unique = unique_likelihood_indices(&outcome.read_likelihoods);
    let stored_best = best_ll_by_index(&outcome.read_likelihoods);
    kv("stored_unique_n", stored_unique.len().to_string());
    kv(
        "stored_hap_n",
        outcome.assembly.haplotypes.len().to_string(),
    );
    kv(
        "java_vs_rust_hap_n",
        format!("java=6 rust={}", outcome.assembly.haplotypes.len()),
    );
    assert_eq!(
        stored_unique.len(),
        2,
        "6R.185 stored hap matrix is Java n=2"
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
        "annotation object is downstream; must stay empty this round"
    );
    let indel = outcome.genotyped_calls.iter().find(|c| {
        c.event.start_1based == GenomePosition::new_1based(CLOSED_INDEL)
            && c.event.ref_allele == "TTC"
    });
    assert!(indel.is_some(), "TTC/T still genotyped");

    kv(
        "TABLE",
        "READ\tJAVA_BEST_LL\tJAVA_FILTER\tRUST_BEST_LL\tRUST_FILTER_LAST_PASS\tRUST_STORED\tTHRESHOLD\tFIRST_DIFFERING_DECISION",
    );

    let extras_keep_caused_filter_n8 = last_extra_n > 0 && last_rust_keep_n == 8;
    let last_pass_would_be_java_n2 = last_java_keep_n == 2;
    let mut stored_includes_java_drops = 0usize;

    for &idx in &stored_unique {
        let rec = outcome.genotyping_reads.get(idx).expect("read");
        let orig = match_java_orig(rec);
        let rust_best = stored_best.get(&idx).copied().unwrap_or(f64::NEG_INFINITY);
        let qlen = rec.qual().len().max(1);
        let thresh = java_log10_min_true_likelihood(qlen);
        let rust_java_pred = java_keep(rust_best, qlen);
        let obs = last_rows
            .iter()
            .find(|r| r.qname == qname(rec) && r.flags == rec.flags());
        let java_best = orig.map(|j| j.norm_max_ll);
        let java_filter = orig.map(|j| if j.java_keep { "KEEP" } else { "DROP" });
        let rust_filter = obs.map(|r| {
            if r.rust_keep {
                if r.extra_retain {
                    "KEEP_EXTRA"
                } else {
                    "KEEP"
                }
            } else {
                "DROP"
            }
        });
        let ll_gap = java_best.map(|j| (rust_best - j).abs()).unwrap_or(f64::NAN);
        let residual_causal = java_best
            .map(|j| java_keep(j, qlen) != rust_java_pred)
            .unwrap_or(false);
        let p12_up = overlaps_p12_cluster_upstream(rec);
        let align_active =
            java_alignment_read_overlaps_interval(rec, covering.start.get(), covering.end.get(), 2);
        if orig.is_some_and(|j| !j.java_keep) {
            stored_includes_java_drops += 1;
        }
        let reason = if orig.is_some_and(|j| j.java_keep) && rust_java_pred {
            "both KEEP on Java static threshold"
        } else if residual_causal {
            "numerical LL residual flips Java predicate"
        } else if obs.is_some_and(|r| r.extra_retain) {
            "Rust extras KEEP (java_equiv_keep=false, rust_keep=true)"
        } else if orig.is_some_and(|j| !j.java_keep) && rust_java_pred {
            "Rust max_ll passes Java threshold (hap/LL population)"
        } else if orig.is_some_and(|j| !j.java_keep)
            && !rust_java_pred
            && obs.is_some_and(|r| !r.rust_keep)
        {
            "last filter pass DROPS like Java; stored KEEP is post-filter overwrite"
        } else {
            "see dumps"
        };
        kv(
            "ROW",
            format!(
                "{} FLAG={} MAPQ={} start={} CIGAR={}\tjava_norm_max={}\t{}\trust_best={:.6}\t{}\tSTORED=yes\tthresh={thresh}\tp12_up={p12_up} align_active={align_active} ll_gap={} extra={} java_equiv={} {}",
                qname(rec),
                rec.flags(),
                rec.mapq(),
                rec.pos() + 1,
                rec.cigar(),
                java_best.map(|v| format!("{v:.6}")).unwrap_or_else(|| "?".into()),
                java_filter.unwrap_or("?"),
                rust_best,
                rust_filter.unwrap_or("?"),
                if ll_gap.is_finite() { format!("{ll_gap:.6}") } else { "?".into() },
                obs.map(|r| r.extra_retain).unwrap_or(false),
                obs.map(|r| r.java_equiv_keep).unwrap_or(false),
                reason,
            ),
        );
        if let Some(o) = orig {
            kv(
                "JAVA_DETAIL",
                format!(
                    "{} FLAG={} end={} orig_len={} best_hap={} prim_max={:.6} java_keep={}",
                    o.qname, o.flag, o.end, o.orig_len, o.best_hap, o.prim_max_ll, o.java_keep
                ),
            );
        }
        if let Some(r) = obs {
            kv(
                "RUST_FILTER_PASS",
                format!(
                    "pass={} {} FLAG={} start={} end={} qlen={} thresh={:.1} max_ll={:.6} n_hap_cells={} n_columns={} java_equiv={} rust_keep={} extra={}",
                    r.pass,
                    r.qname,
                    r.flags,
                    r.start_1based,
                    r.end_1based,
                    r.qual_len,
                    r.threshold,
                    r.max_ll,
                    r.n_hap_cells,
                    r.n_columns,
                    r.java_equiv_keep,
                    r.rust_keep,
                    r.extra_retain,
                ),
            );
        }
        assert!(
            !residual_causal,
            "LL residual must not flip KEEP/DROP vs Java at this site"
        );
        assert!(
            !p12_up,
            "P12 cluster-upstream interval is 92305716-728; these reads must not overlap it"
        );
    }

    kv(
        "stored_includes_java_drop_n",
        stored_includes_java_drops.to_string(),
    );
    assert_eq!(stored_includes_java_drops, 0);

    kv(
        "filter_pass_java_equiv_matches_java_n2",
        last_pass_would_be_java_n2.to_string(),
    );
    kv(
        "filter_pass_extras_keep_n8",
        extras_keep_caused_filter_n8.to_string(),
    );

    // First causal arrow: either extras KEEP at the filter, or a later apply_normalize=false
    // refresh that restores Java-dropped reads after a Java-equivalent drop.
    let first_arrow = if extras_keep_caused_filter_n8 {
        "filter_poorly_modeled_region_read_likelihoods extras KEEP (predicate), not LL residual"
    } else if last_pass_would_be_java_n2 && stored_unique.len() == 2 {
        "6R.185 restored Java-order normalize+filter after last P12 refresh; stored n=2"
    } else if last_pass_would_be_java_n2 && stored_unique.len() == 8 {
        "post-filter refresh_region_read_likelihoods(apply_normalize=false) overwrites Java-equivalent n=2 with raw n=8"
    } else {
        "see observe dumps"
    };
    kv("first_arrow", first_arrow);
    kv(
        "classification",
        if extras_keep_caused_filter_n8 {
            "D — WRONG UPSTREAM GENOTYPING SUBSET (Rust poorly-modeled extras KEEP predicate)"
        } else if last_pass_would_be_java_n2 {
            "D — WRONG UPSTREAM GENOTYPING SUBSET (filter result overwritten by unfiltered refresh)"
        } else {
            "D — WRONG UPSTREAM GENOTYPING SUBSET"
        },
    );

    assert!(
        last_pass_would_be_java_n2 || extras_keep_caused_filter_n8,
        "must pin either extras KEEP or post-filter overwrite"
    );
}
