//! 6R.182: construct the Java-equivalent per-variant annotation `AlleleLikelihoods`
//! at `2:92307333 T/G` (proof-only).
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! PRODUCTION CHANGE: NONE.
//!
//! Live Java object dump (`genotype-emit-at-loc` 92307333) replaces 6R.181 INFO
//! inference. Sequence:
//!   stored hap `AlleleLikelihoods<GATKRead,Haplotype>` n=2 after
//!     `filterPoorlyModeledEvidence` (8→2)
//!   NEW `AlleleLikelihoods<GATKRead,Allele>` = `marginalize(alleleMapper)` n=2
//!   `retainEvidence(target.overlaps)` 2→1
//!   `prepareReadAlleleLikelihoodsForAnnotation` reuses that object
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r182_per_variant_annotation_object -- --nocapture --test-threads=1
//! HOLDOUT_6R182=1 cargo test -p gatk-haplotypecaller --test holdout_6r182_per_variant_object -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{build_per_haplotype_variation_events, VariationEvent};
use gatk_haplotypecaller::hc_allele_mapping::create_allele_mapper_with_events;
use gatk_haplotypecaller::hc_genotyping_engine::java_alignment_read_overlaps_interval;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, marginalize_rows_to_biallelic_alleles,
    region_likelihoods_to_rows, traverse_assembly_region_walker, try_emit_call_region_variants,
    AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition, HaplotypeCallerEngine,
    ReadFilterParams, RegionReadLikelihood, SharedBamRecord, WalkerTraversalConfig,
    DEFAULT_STAND_EMIT_CONFIDENCE,
};
use rust_htslib::bam::Record;
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "2:92307200-92307550";
const CLOSED_SNP: u64 = 92_305_634;
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_VCF_REL: &str = "parity/reports/6r43/p12_indel_mix/java.vcf";
const TARGET: u64 = 92_307_333;
const MERGED_REF: &str = "T";
const MERGED_ALT: &str = "G";
const CLOSED_INDEL: u64 = 92_307_324;
const CLOSED_INDEL_REF: &str = "TTC";
const CLOSED_INDEL_ALT: &str = "T";
const FLAG_REVERSE: u16 = 0x10;
const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const MARGIN: i32 = 2;
const JAVA_MQ44_QNAME: &str = "H06JUADXX130110:1:1101:10052:88682";
const EXPECTED_ERROR_RATE_PER_BASE: f64 = 0.02;
const LOG10_QUAL_PER_ERROR: f64 = -4.0;

/// Live Java `genotype-emit-at-loc` 92307333 (`/tmp/gatk-rs-6r182/genotype-emit-92307333.tsv`).
/// Identity is (QNAME, FLAG, start), not QNAME collapse.
#[derive(Clone, Copy)]
struct JavaOrigRead {
    qname: &'static str,
    flag: u16,
    mapq: u8,
    start: u64,
    cigar: &'static str,
    orig_len: usize,
    norm_max_ll: f64,
    java_poorly_modeled_keep: bool,
    java_hap_ll_at_loc: bool,
    java_pre_retain: bool,
    java_post_retain: bool,
}

const JAVA_ORIG: &[JavaOrigRead] = &[
    JavaOrigRead {
        qname: "H06HDADXX130110:1:1101:10061:17286",
        flag: 163,
        mapq: 40,
        start: 92_307_248,
        cigar: "56H112M1D8M74H",
        orig_len: 120,
        norm_max_ll: -13.056375503540,
        java_poorly_modeled_keep: false,
        java_hap_ll_at_loc: false,
        java_pre_retain: false,
        java_post_retain: false,
    },
    JavaOrigRead {
        qname: "H06JUADXX130110:1:1101:10011:51168",
        flag: 81,
        mapq: 22,
        start: 92_307_248,
        cigar: "4H112M1D78M56H",
        orig_len: 190,
        norm_max_ll: -16.238531112671,
        java_poorly_modeled_keep: false,
        java_hap_ll_at_loc: false,
        java_pre_retain: false,
        java_post_retain: false,
    },
    JavaOrigRead {
        qname: "H06JUADXX130110:1:1101:10052:88682",
        flag: 163,
        mapq: 50,
        start: 92_307_260,
        cigar: "179M71H",
        orig_len: 179,
        norm_max_ll: -19.377605438232,
        java_poorly_modeled_keep: false,
        java_hap_ll_at_loc: false,
        java_pre_retain: false,
        java_post_retain: false,
    },
    JavaOrigRead {
        qname: "H06HDADXX130110:1:1101:10073:74407",
        flag: 163,
        mapq: 41,
        start: 92_307_272,
        cigar: "167M83H",
        orig_len: 167,
        norm_max_ll: -10.475610733032,
        java_poorly_modeled_keep: false,
        java_hap_ll_at_loc: false,
        java_pre_retain: false,
        java_post_retain: false,
    },
    JavaOrigRead {
        qname: "H06HDADXX130110:1:1101:10073:74407",
        flag: 83,
        mapq: 41,
        start: 92_307_272,
        cigar: "31H167M52H",
        orig_len: 167,
        norm_max_ll: -10.475610733032,
        java_poorly_modeled_keep: false,
        java_hap_ll_at_loc: false,
        java_pre_retain: false,
        java_post_retain: false,
    },
    JavaOrigRead {
        qname: JAVA_MQ44_QNAME,
        flag: 83,
        mapq: 44,
        start: 92_307_292,
        cigar: "14H147M89H",
        orig_len: 147,
        norm_max_ll: -5.817768096924,
        java_poorly_modeled_keep: true,
        java_hap_ll_at_loc: true,
        java_pre_retain: true,
        java_post_retain: true,
    },
    JavaOrigRead {
        qname: "H06HDADXX130110:1:1101:10061:17286",
        flag: 83,
        mapq: 40,
        start: 92_307_338,
        cigar: "1H22M1D78M149H",
        orig_len: 100,
        norm_max_ll: -2.494571685791,
        java_poorly_modeled_keep: true,
        java_hap_ll_at_loc: true,
        java_pre_retain: true,
        java_post_retain: false,
    },
    JavaOrigRead {
        qname: "H06HDADXX130110:2:1101:10097:72839",
        flag: 147,
        mapq: 23,
        start: 92_307_367,
        cigar: "83H72M95H",
        orig_len: 72,
        norm_max_ll: -18.344545364380,
        java_poorly_modeled_keep: false,
        java_hap_ll_at_loc: false,
        java_pre_retain: false,
        java_post_retain: false,
    },
];

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R182\t{key}\t{}", value.as_ref());
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

fn info_f64(info: &[InfoValue], key: &str) -> Option<f64> {
    for v in info {
        if let InfoValue::Float(k, xs) = v {
            if k == key {
                return xs.first().copied();
            }
        }
    }
    None
}

fn info_has(info: &[InfoValue], key: &str) -> bool {
    info.iter().any(|v| match v {
        InfoValue::Flag(k)
        | InfoValue::Integer(k, _)
        | InfoValue::Float(k, _)
        | InfoValue::String(k, _)
        | InfoValue::Character(k, _) => k == key,
    })
}

fn qname(rec: &Record) -> String {
    String::from_utf8_lossy(rec.qname()).into_owned()
}

fn strand_label(rec: &Record) -> &'static str {
    if rec.flags() & FLAG_REVERSE != 0 {
        "rev"
    } else {
        "fwd"
    }
}

fn read_id(rec: &Record) -> String {
    format!(
        "{} FLAG={} start={}",
        qname(rec),
        rec.flags(),
        rec.pos() + 1
    )
}

fn java_log10_min_true_likelihood(qualified_read_len: usize) -> f64 {
    let max_errors = (qualified_read_len as f64 * EXPECTED_ERROR_RATE_PER_BASE)
        .ceil()
        .min(2.0);
    max_errors * LOG10_QUAL_PER_ERROR
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

fn parse_info_map(info: &str) -> BTreeMap<String, String> {
    let mut out = BTreeMap::new();
    for piece in info.split(';') {
        if piece.is_empty() {
            continue;
        }
        match piece.split_once('=') {
            Some((k, v)) => {
                out.insert(k.to_string(), v.to_string());
            }
            None => {
                out.insert(piece.to_string(), String::new());
            }
        }
    }
    out
}

fn java_target_record(path: &Path) -> (BTreeMap<String, String>, String, String, String) {
    let text = std::fs::read_to_string(path).expect("java.vcf");
    for line in text.lines() {
        if line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 10 {
            continue;
        }
        if f[0] == "2" && f[1] == "92307333" && f[3] == MERGED_REF && f[4] == MERGED_ALT {
            return (
                parse_info_map(f[7]),
                f[5].to_string(),
                f[8].to_string(),
                f[9].to_string(),
            );
        }
    }
    panic!("java.vcf missing 2:92307333 T/G");
}

fn match_java_orig(rec: &Record) -> Option<&'static JavaOrigRead> {
    let qn = qname(rec);
    let flag = rec.flags();
    let start = (rec.pos() + 1) as u64;
    if let Some(exact) = JAVA_ORIG
        .iter()
        .find(|j| j.qname == qn && j.flag == flag && j.start == start)
    {
        return Some(exact);
    }
    // Soft-clip/CIGAR can shift start by 1 (Java 92307248 vs Rust 92307249 on FLAG=163).
    JAVA_ORIG.iter().find(|j| j.qname == qn && j.flag == flag)
}

/// Copy of production `per_variant_annotation_likelihoods` (test-only reconstruction).
fn hypothetical_annotation_from_subset(
    subset: &[RegionReadLikelihood],
    reads: &[SharedBamRecord],
) -> Vec<RegionReadLikelihood> {
    if subset.is_empty() {
        return Vec::new();
    }
    let mut best_ll: BTreeMap<usize, f64> = BTreeMap::new();
    for cell in subset {
        let idx = cell.read_index.get();
        let ll = cell.log10_likelihood;
        best_ll
            .entry(idx)
            .and_modify(|m| {
                if ll > *m {
                    *m = ll;
                }
            })
            .or_insert(ll);
    }
    let mut qnames: BTreeSet<Vec<u8>> = BTreeSet::new();
    for &idx in best_ll.keys() {
        if let Some(rec) = reads.get(idx) {
            qnames.insert(rec.qname().to_owned());
        }
    }
    if qnames.len() == 1 && best_ll.len() > 1 {
        let Some((&keep_idx, _)) = best_ll.iter().max_by(|a, b| a.1.total_cmp(b.1)) else {
            return subset.to_vec();
        };
        return subset
            .iter()
            .filter(|cell| cell.read_index.get() == keep_idx)
            .cloned()
            .collect();
    }
    subset.to_vec()
}

#[test]
fn forensic_6r182_java_object_from_live_dump() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    kv(
        "java_dump",
        "genotype-emit-at-loc 92307333 /tmp/gatk-rs-6r182/genotype-emit-92307333.tsv",
    );
    kv("java_hap_count", "6");
    kv("java_orig_n", "8");
    kv("java_stored_hap_ll_n", "2");
    kv("java_pre_retain_n", "2");
    kv("java_post_retain_n", "1");
    kv(
        "java_annotation_reuse",
        "prepareReadAlleleLikelihoodsForAnnotation contamination-off reuses genotyping object",
    );
    kv("java_merged_alleles", "[T*, G]");
    kv(
        "java_mapper",
        "T* ← 3 haps (cc52c0e6e338ead1, af6e0f986fc324c8 ref, 77c687de18319e97); G ← 3 haps (9188a56b92ef6932, 9415a7310df771a1, d926aa50c2c117aa)",
    );
    kv(
        "java_retained_read",
        "H06JUADXX130110:1:1101:10052:88682 FLAG=83 MAPQ=44 start=92307292 CIGAR=14H147M89H best=G=-5.817768 (T*=-10.317768)",
    );
    kv(
        "java_dropped_by_retain",
        "H06HDADXX130110:1:1101:10061:17286 FLAG=83 MAPQ=40 start=92307338 (hap_ll n=2; starts after 92307333±2)",
    );

    let java_keep_n = JAVA_ORIG
        .iter()
        .filter(|r| r.java_poorly_modeled_keep)
        .count();
    let java_post = JAVA_ORIG.iter().filter(|r| r.java_post_retain).count();
    assert_eq!(java_keep_n, 2);
    assert_eq!(java_post, 1);
    let retained = JAVA_ORIG
        .iter()
        .find(|r| r.java_post_retain)
        .expect("retained");
    assert_eq!(retained.qname, JAVA_MQ44_QNAME);
    assert_eq!(retained.flag, 83);
    assert_eq!(retained.mapq, 44);
    assert!(retained.java_poorly_modeled_keep);
    assert!(java_keep(-5.817768096924, 147));
    assert!(java_keep(-2.494571685791, 100));
    assert!(!java_keep(-10.475610733032, 167));
    assert!(!java_keep(-13.056375503540, 120));
    // Direct object, not INFO inference: Coverage.evidenceCount of post-retain n=1.
    kv(
        "java_n1_source",
        "DIRECT AlleleLikelihoods post-retainEvidence read_ll row, not 6R.181 INFO reconstruction",
    );
}

#[test]
fn forensic_6r182_source_contracts_no_production_change() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    let early = include_str!("../src/hc_genotyping_engine/genotype_site_early_template.rs");
    let fin = include_str!("../src/hc_genotyping_engine/genotype_finalize.rs");
    let pipe = include_str!("../src/hc_genotyping_engine/genotype_site_pipeline.rs");
    let ge = include_str!("../src/hc_genotyping_engine/mod.rs");
    let emit = include_str!("../src/region_vcf_emit.rs");
    let ann = include_str!("../src/variant_site_hc_annotations.rs");
    let ll = include_str!("../src/engine_likelihoods.rs");
    let engine = include_str!("../src/engine.rs");
    assert!(
        early.contains("if is_cluster_tg_snp(&event)")
            && early.contains("finish_strict_java_shaped_site_call"),
        "cluster-TG still selected in SiteEarlyTemplate::try_shaped"
    );
    assert!(
        early.contains("likelihoods: &[RegionReadLikelihood]"),
        "try_shaped already receives the region hap matrix"
    );
    let tg_arm = early
        .split("if is_cluster_tg_snp(&event)")
        .nth(1)
        .expect("tg arm");
    let tg_body = tg_arm
        .split("if is_cluster_tc_snp(&event)")
        .next()
        .expect("tg body");
    assert!(
        !tg_body.contains("likelihood_subset_for_event")
            && !tg_body.contains("with_annotation_likelihoods")
            && !tg_body.contains("per_variant_annotation_likelihoods"),
        "cluster-TG arm still does not construct a per-variant likelihood object"
    );
    assert!(
        fin.contains("fn finish_strict_java_shaped_site_call")
            && fin.contains("GenotypedSiteCall::new(event, genotype)"),
        "shaped path still constructs GenotypedSiteCall::new"
    );
    let fin_fn = fin
        .split("fn finish_strict_java_shaped_site_call")
        .nth(1)
        .expect("finish fn");
    let fin_body = fin_fn
        .split("fn finalize_strict_java_variation_genotype")
        .next()
        .expect("body");
    assert!(
        !fin_body.contains("with_annotation_likelihoods"),
        "early-template finish must not attach annotation_likelihoods"
    );
    assert!(
        pipe.contains("SiteEarlyTemplate::try_shaped") && pipe.contains("return Ok(Some(call));"),
        "early template still returns before subset construction"
    );
    assert!(
        emit.contains("if call.annotation_likelihoods.is_empty()")
            && emit.contains("likelihoods: region_wide.likelihoods"),
        "empty annotation_likelihoods still falls back to region-wide"
    );
    assert!(
        ge.contains("fn per_variant_annotation_likelihoods")
            && pipe.contains("with_annotation_likelihoods"),
        "6R.180 helper still exists on the later SiteScore path only"
    );
    assert!(
        engine.contains("retain_marginal_p12_cluster_upstream_read")
            && engine.contains("retain_marginal_sparse_softclip_read"),
        "poorly-modeled extras KEEP helpers still present (not this round's production target)"
    );
    assert!(
        ann.contains("6R.174: INFO DP is Java `Coverage.evidenceCount`")
            && ann.contains("mq_rms_of_sample_evidence")
            && ann.contains("strand_bias_sample_counts"),
        "DP/MQ/SOR formulas must stay closed"
    );
    assert!(
        ll.contains("score_pairhmm_from_records_java_mate_contig"),
        "6R.176 mate-contig gate must stay closed"
    );
    assert!(
        emit.contains("INBREEDING_COEFF_MIN_SAMPLES"),
        "6R.178 InbreedingCoeff gate must stay closed"
    );
    assert!(
        !ann.contains("92307333") && !emit.contains("92307333") && !pipe.contains("92307333"),
        "no locus-specific 6R.182 production patch"
    );
}

#[test]
fn forensic_6r182_construct_candidate_on_cluster_tg_path() {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    let java_vcf = root.join(JAVA_VCF_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing P12 ref/BAM");
        return;
    }

    kv("java_pin", JAVA_PIN);
    kv("target", "2:92307333 T/G");
    kv("production_change", "NONE");
    kv(
        "classification_target",
        "ANNOTATION_INPUT_MEMBERSHIP_DIVERGENCE",
    );

    if java_vcf.is_file() {
        let (jinfo, jqual, jfmt, jsample) = java_target_record(&java_vcf);
        assert_eq!(jinfo.get("DP").map(String::as_str), Some("1"));
        assert_eq!(jinfo.get("MQ").map(String::as_str), Some("44.00"));
        assert_eq!(jinfo.get("SOR").map(String::as_str), Some("1.609"));
        let fields: Vec<&str> = jfmt.split(':').collect();
        let vals: Vec<&str> = jsample.split(':').collect();
        let get = |name: &str| {
            fields
                .iter()
                .zip(vals.iter())
                .find(|(k, _)| **k == name)
                .map(|(_, v)| *v)
                .unwrap_or("")
        };
        assert_eq!(get("GT"), "1/1");
        assert_eq!(get("AD"), "0,1");
        assert_eq!(get("DP"), "1");
        assert_eq!(get("GQ"), "3");
        assert_eq!(get("PL"), "45,3,0");
        assert!((jqual.parse::<f64>().unwrap_or(0.0) - 35.48).abs() < 0.02);
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
    let outcome = HaplotypeCallerEngine::call_region(
        covering,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("outcome");

    let call = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(TARGET)
                && c.event.ref_allele == MERGED_REF
                && c.event.alt_allele == MERGED_ALT
        })
        .expect("T/G call");
    assert_eq!(call.genotype.format.ad_as_i32(), vec![0, 1]);
    assert_eq!(call.genotype.format.dp.as_i32(), 1);
    assert_eq!(call.genotype.format.gq.as_i32(), 3);
    assert_eq!(call.genotype.format.pl_as_i32(), vec![45, 3, 0]);
    assert!(
        call.annotation_likelihoods.is_empty(),
        "cluster-TG path must still leave annotation_likelihoods empty"
    );

    let indel = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_INDEL)
                && c.event.ref_allele == CLOSED_INDEL_REF
                && c.event.alt_allele == CLOSED_INDEL_ALT
        })
        .expect("TTC/T call");
    let indel_unique = unique_likelihood_indices(&indel.annotation_likelihoods);
    assert_eq!(indel_unique.len(), 1, "6R.180 TTC/T attachment stays");

    kv(
        "rust_hap_count",
        outcome.assembly.haplotypes.len().to_string(),
    );
    kv(
        "java_vs_rust_hap_count",
        format!("java=6 rust={}", outcome.assembly.haplotypes.len()),
    );

    let region_unique = unique_likelihood_indices(&outcome.read_likelihoods);
    let best_ll = best_ll_by_index(&outcome.read_likelihoods);
    kv("region_wide_unique_n", region_unique.len().to_string());
    kv(
        "genotyping_reads_n",
        outcome.genotyping_reads.len().to_string(),
    );
    assert_eq!(region_unique.len(), 2);

    let overlap: BTreeSet<usize> = outcome
        .genotyping_reads
        .iter()
        .enumerate()
        .filter(|(_, r)| java_alignment_read_overlaps_interval(r, TARGET, TARGET, MARGIN))
        .map(|(i, _)| i)
        .filter(|i| region_unique.contains(i))
        .collect();
    kv("naive_overlap_on_stored_n", overlap.len().to_string());
    assert_eq!(
        overlap.len(),
        1,
        "naive retainEvidence on stored n=2 is n=1; still not attached"
    );

    let overlap_cells: Vec<RegionReadLikelihood> = outcome
        .read_likelihoods
        .iter()
        .filter(|c| overlap.contains(&c.read_index.get()))
        .cloned()
        .collect();
    let hypo = hypothetical_annotation_from_subset(&overlap_cells, &outcome.genotyping_reads);
    let hypo_unique = unique_likelihood_indices(&hypo);
    kv(
        "naive_overlap_then_6r180_helper_n",
        hypo_unique.len().to_string(),
    );
    assert_eq!(hypo_unique.len(), 1);

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
    kv(
        "rust_mapper",
        format!(
            "ref_pool={:?} alt_pool={:?}",
            mapping
                .ref_haplotype_indices
                .iter()
                .map(|i| i.get())
                .collect::<Vec<_>>(),
            mapping
                .alt_haplotype_indices
                .iter()
                .map(|i| i.get())
                .collect::<Vec<_>>()
        ),
    );
    let rows =
        region_likelihoods_to_rows(&outcome.read_likelihoods, outcome.assembly.haplotypes.len());
    let marg = marginalize_rows_to_biallelic_alleles(
        &rows,
        &mapping.ref_haplotype_indices,
        &mapping.alt_haplotype_indices,
    );
    let allele_by_idx: BTreeMap<usize, (String, f64, f64)> = marg
        .iter()
        .map(|row| {
            let lls = &row.haplotype_log10_likelihoods;
            let ll_ref = lls.first().copied().unwrap_or(f64::NEG_INFINITY);
            let ll_alt = lls.get(1).copied().unwrap_or(f64::NEG_INFINITY);
            let best = if ll_alt > ll_ref { "G" } else { "T" };
            (row.read_index, (best.to_string(), ll_ref, ll_alt))
        })
        .collect();

    let mut java_threshold_keep: BTreeSet<usize> = BTreeSet::new();
    let mut java_recipe: BTreeSet<usize> = BTreeSet::new();
    let mut extras_fail_java_threshold = true;

    kv(
        "TABLE",
        "READ\tJAVA_PER_VARIANT\tRUST_CANDIDATE\tDIFFERENCE_REASON",
    );
    for &idx in &region_unique {
        let rec = outcome.genotyping_reads.get(idx).expect("read");
        let orig = match_java_orig(rec);
        let ov = overlap.contains(&idx);
        let max_ll = best_ll.get(&idx).copied().unwrap_or(f64::NEG_INFINITY);
        let qlen = rec.qual().len().max(1);
        let thresh = java_log10_min_true_likelihood(qlen);
        let rust_java_keep = java_keep(max_ll, qlen);
        if rust_java_keep {
            java_threshold_keep.insert(idx);
            if ov {
                java_recipe.insert(idx);
            }
        }
        let java_post = orig.map(|j| j.java_post_retain).unwrap_or(false);
        let java_hap = orig.map(|j| j.java_hap_ll_at_loc).unwrap_or(false);
        let java_pm = orig.map(|j| j.java_poorly_modeled_keep).unwrap_or(false);
        if ov && !java_post && rust_java_keep {
            extras_fail_java_threshold = false;
        }
        let (best_allele, ll_t, ll_g) = allele_by_idx
            .get(&idx)
            .cloned()
            .unwrap_or_else(|| ("?".into(), f64::NEG_INFINITY, f64::NEG_INFINITY));
        let reason = if java_post && ov && rust_java_keep {
            "both: Java post-retain and Rust java-threshold+overlap"
        } else if java_post {
            "Java post-retainEvidence n=1"
        } else if java_hap && !ov {
            "Java hap_ll n=2; dropped by retainEvidence (no overlap 92307333±2)"
        } else if !java_pm && ov {
            "Java filterPoorlyModeledEvidence DROP (norm max_ll < -8); Rust extras KEEP still in n=8; naive overlap keeps"
        } else if !java_pm {
            "Java filterPoorlyModeledEvidence DROP; not in hap_ll at loc; no overlap"
        } else {
            "other"
        };
        kv(
            "ROW",
            format!(
                "{}\tJAVA_POST={}\tNAIVE_OV={}\tJAVA_THRESH_KEEP={}\tJAVA_RECIPE={}\t{}",
                read_id(rec),
                java_post,
                ov,
                rust_java_keep,
                rust_java_keep && ov,
                reason
            ),
        );
        kv(
            "LL_READ",
            format!(
                "idx={idx} {} MAPQ={} strand={} CIGAR={} java_cigar={} java_orig_len={} rust_best_ll={max_ll:.6} thresh={thresh} rust_java_keep={rust_java_keep} java_norm_max={} java_pm_keep={} java_hap_ll={} java_pre_retain={} java_post_retain={} overlap={} best_allele={} T={ll_t:.6} G={ll_g:.6}",
                read_id(rec),
                rec.mapq(),
                strand_label(rec),
                rec.cigar(),
                orig.map(|j| j.cigar).unwrap_or("?"),
                orig.map(|j| j.orig_len).unwrap_or(0),
                orig.map(|j| format!("{:.6}", j.norm_max_ll)).unwrap_or_else(|| "?".into()),
                java_pm,
                java_hap,
                orig.map(|j| j.java_pre_retain).unwrap_or(false),
                java_post,
                ov,
                best_allele,
            ),
        );
    }

    kv(
        "java_threshold_keep_on_rust_n",
        java_threshold_keep.len().to_string(),
    );
    kv(
        "java_threshold_then_overlap_n",
        java_recipe.len().to_string(),
    );
    kv(
        "extras_fail_java_threshold_on_rust_lls",
        extras_fail_java_threshold.to_string(),
    );

    for &idx in &java_recipe {
        let rec = &outcome.genotyping_reads[idx];
        kv("JAVA_RECIPE_READ", read_id(rec));
        assert_eq!(qname(rec), JAVA_MQ44_QNAME);
        assert_eq!(rec.flags(), 83);
        assert_eq!(rec.mapq(), 44);
    }

    // The Java-equivalent candidate on Rust likelihood *values* (static poorly-modeled
    // then overlap) is n=1 FLAG=83 MAPQ=44. Stored n=2 overlap is also n=1; still not attached.
    assert_eq!(
        java_recipe.len(),
        1,
        "Java-recipe (threshold then overlap) on Rust LLs must be the MAPQ=44 FLAG=83 read"
    );
    assert!(
        extras_fail_java_threshold,
        "no extra overlap read may pass Java poorly-modeled on Rust best_ll"
    );

    let mq44 = region_unique.iter().filter(|&&idx| {
        outcome
            .genotyping_reads
            .get(idx)
            .is_some_and(|r| qname(r) == JAVA_MQ44_QNAME && r.flags() == 83 && r.mapq() == 44)
    });
    assert_eq!(mq44.count(), 1);

    kv(
        "earliest_construction_point",
        "SiteEarlyTemplate::try_shaped already has likelihoods + mapping + event + reads; finish_strict_java_shaped_site_call does not use them for annotation",
    );
    kv(
        "java_semantic_recipe",
        "stored_hap_ll.marginalize(alleleMapper) [NEW] -> retainEvidence(overlap) [in-place] -> reuse for annotation",
    );
    kv(
        "naive_candidate_is_not_java",
        "overlap on stored n=2 = 1; Java object is still a NEW per-variant AlleleLikelihoods (6R.181 B); do not attach here",
    );
    kv(
        "first_arrow",
        "AlleleLikelihoods.filterPoorlyModeledEvidence membership is now stored n=2 (6R.185). Cluster-TG missing construction is later.",
    );
    kv(
        "classification",
        "D — WRONG UPSTREAM GENOTYPING SUBSET (stored hap AlleleLikelihoods). 6R.181 B remains true as a later missing construction.",
    );

    let emitted =
        try_emit_call_region_variants(covering, &outcome, "NA12878", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("emit");
    let rec = emitted
        .iter()
        .find(|r| {
            r.position == TARGET
                && r.reference == MERGED_REF
                && r.alternate.first().map(String::as_str) == Some(MERGED_ALT)
        })
        .expect("T/G record");
    let sample = rec.samples.first().expect("sample");
    assert_eq!(
        sample.gt.as_ref().map(|g| g.to_string()).as_deref(),
        Some("1/1")
    );
    assert_eq!(sample.ad.as_deref(), Some(&[0u32, 1][..]));
    assert_eq!(sample.dp.map(|d| d as i32), Some(1));
    assert_eq!(sample.gq.map(|g| g as i32), Some(3));
    assert_eq!(sample.pl.as_deref(), Some(&[45u32, 3, 0][..]));
    assert!((rec.quality.unwrap_or(0.0) - 35.48).abs() < 0.05);
    assert_eq!(
        info_i32(&rec.info, "DP"),
        Some(1),
        "region-wide fallback on stored n=2"
    );
    let mq = info_f64(&rec.info, "MQ").unwrap_or(-1.0);
    let sor = info_f64(&rec.info, "SOR").unwrap_or(-1.0);
    assert!((mq - 42.05).abs() < 0.005, "MQ={mq}");
    assert!((sor - 0.6931471805599453).abs() < 1e-9, "SOR={sor}");

    let closed_indel = emitted
        .iter()
        .find(|r| r.position == CLOSED_INDEL && r.reference == CLOSED_INDEL_REF)
        .expect("TTC/T");
    assert_eq!(info_i32(&closed_indel.info, "DP"), Some(1));
    let closed_mq = info_f64(&closed_indel.info, "MQ").unwrap_or(-1.0);
    let closed_sor = info_f64(&closed_indel.info, "SOR").unwrap_or(-1.0);
    assert!((closed_mq - 44.0).abs() < 0.005);
    assert!((closed_sor - 1.6094379124341003).abs() < 1e-3 || (closed_sor - 1.609).abs() < 0.002);

    let closed_specs = parse_intervals_cli_string(&dict, "2:92305500-92305850").expect("closed");
    let closed_walk = traverse_assembly_region_walker(
        &dict,
        &closed_specs,
        &ref_fasta,
        &bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
    )
    .expect("closed walk");
    let closed_regions = flatten_assembly_regions(&closed_walk);
    let closed_covering = closed_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_SNP
                && r.end.get() >= CLOSED_SNP
        })
        .expect("closed covering");
    let closed_outcome = HaplotypeCallerEngine::call_region(
        closed_covering,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("closed call")
    .expect("closed outcome");
    let closed_emitted = try_emit_call_region_variants(
        closed_covering,
        &closed_outcome,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("closed emit");
    let closed_rec = closed_emitted
        .iter()
        .find(|r| r.position == CLOSED_SNP)
        .expect("closed G/T");
    assert_eq!(info_i32(&closed_rec.info, "DP"), Some(3), "6R.174 stays");
    assert!(
        !info_has(&closed_rec.info, "InbreedingCoeff"),
        "6R.178 stays"
    );
    let closed_snp_sor = info_f64(&closed_rec.info, "SOR").unwrap_or(-1.0);
    assert!(
        (closed_snp_sor - 0.6931471805599453).abs() < 1e-6
            || (closed_snp_sor - 0.693).abs() < 0.002
    );
}
