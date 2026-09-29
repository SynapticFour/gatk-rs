//! 6R.308: proof that `java_read_coordinate_compare` matches Java
//! `ReadCoordinateComparator`. Production uses that comparator by default.
//! This test still sets it explicitly and replays the pre-6R.312 name sort
//! in memory. Likelihood values are not injected.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r308_java_read_coordinate_comparator -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::assembly_region_finalize::{
    java_read_coordinate_compare, set_forensic_6r308_java_read_coordinate_order,
};
use gatk_haplotypecaller::hc_genotyping_engine::diploid_genotype_log10_likelihoods_from_allele_rows;
use gatk_haplotypecaller::read_assembly_filter::{passes_assembly_read, AssemblyReadFilterConfig};
use gatk_haplotypecaller::{
    begin_forensic_6r307_order_observe, begin_hap_list_observe, begin_trim_calc_observe,
    call_disposition, flatten_assembly_regions, forensic_6r306_excluded_count,
    set_forensic_6r294_padded_span, set_forensic_6r306_exclude_failed_mate,
    take_colocated_merge_numerics, take_forensic_6r307_order, take_hap_list_snaps,
    take_hap_list_trim_span, traverse_assembly_region_walker, try_emit_call_region_variants,
    AssemblyRegionCallDisposition, CallRegionArgs, Forensic6r307Read, HapListSnap,
    HaplotypeCallerEngine, ReadFilterParams, ReadLikelihoodRow, WalkerTraversalConfig,
    DEFAULT_STAND_EMIT_CONFIDENCE,
};
use rust_htslib::bam::record::{Cigar, CigarString};
use rust_htslib::bam::{HeaderView, Record};
use std::collections::{BTreeSet, HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::sync::Arc;

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_TSV: &str = "gatk-haplotypecaller/tests/6r297_java_hap_ll.tsv";
const LOCUS: u64 = 29_455_649;
const ALT: &str = "TGTTTG";
const JAVA_SPAN: (u64, u64) = (29_455_560, 29_455_728);
const DIAG_QNAME: &str = "HISEQ1:11:H8GV6ADXX:1:2116:18670:99941";
const DIAG_FLAGS: u16 = 99;
const MATE_QNAME: &str = "HISEQ1:11:H8GV6ADXX:2:1103:14252:55237";
const MATE_FLAGS: u16 = 97;
const FORWARD_QNAME: &str = "HISEQ1:11:H8GV6ADXX:1:1116:4033:65919";
const FORWARD_FLAGS: u16 = 99;
const REVERSE_QNAME: &str = "HISEQ1:11:H8GV6ADXX:1:1109:9994:7054";
const REVERSE_FLAGS: u16 = 83;

type Key = (String, u16);

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R308\t{key}\t{}", value.as_ref());
}

fn fmt_f(x: f64) -> String {
    format!("{x:.17}")
}

fn cell_bits(cell: &str) -> f64 {
    let hex = cell.split("bits=").nth(1).unwrap().trim();
    f64::from_bits(u64::from_str_radix(hex.trim_start_matches("0x"), 16).unwrap())
}

fn field_cells(text: &str, key: &str) -> Vec<f64> {
    for line in text.lines() {
        let mut p = line.split('\t');
        if p.next() == Some("6R297") && p.next() == Some(key) {
            return p.map(cell_bits).collect();
        }
    }
    panic!("missing {key}");
}

fn header() -> Arc<HeaderView> {
    Arc::new(HeaderView::from_bytes(
        b"@HD\tVN:1.0\n@SQ\tSN:20\tLN:100000000\n@SQ\tSN:7\tLN:100000000\n",
    ))
}

fn rec(name: &str, flags: u16, pos_1based: i64) -> Record {
    let mut r = Record::new();
    r.set_header(header());
    r.set(
        name.as_bytes(),
        Some(&CigarString::from(vec![Cigar::Match(10)])),
        b"AAAAAAAAAA",
        &vec![30u8; 10],
    );
    r.set_tid(0);
    r.set_pos(pos_1based - 1);
    r.set_flags(flags);
    r.set_mapq(60);
    r.set_insert_size(0);
    r
}

fn antisym(a: &Record, b: &Record) {
    assert_eq!(
        java_read_coordinate_compare(a, b),
        -java_read_coordinate_compare(b, a)
    );
}

#[test]
fn comparator_forward_before_reverse_at_same_start() {
    let forward = rec("fwd", 99, 29_455_560);
    let reverse = rec("rev", 83, 29_455_560);
    assert!(!forward.is_reverse());
    assert!(reverse.is_reverse());
    assert!(java_read_coordinate_compare(&forward, &reverse) < 0);
    antisym(&forward, &reverse);
}

#[test]
fn comparator_name_after_start_and_orientation() {
    let a = rec("AAA", 99, 29_455_560);
    let b = rec("BBB", 99, 29_455_560);
    assert!(java_read_coordinate_compare(&a, &b) < 0);
    antisym(&a, &b);
}

#[test]
fn comparator_later_java_tie_breakers() {
    let mut low_flags = rec("SAME", 99, 100);
    let mut high_flags = rec("SAME", 163, 100);
    assert_eq!(java_read_coordinate_compare(&low_flags, &high_flags), -1);
    antisym(&low_flags, &high_flags);

    low_flags.set_flags(99);
    high_flags.set_flags(99);
    low_flags.set_mapq(20);
    high_flags.set_mapq(60);
    assert_eq!(java_read_coordinate_compare(&low_flags, &high_flags), -1);

    low_flags.set_mapq(60);
    high_flags.set_mapq(60);
    low_flags.set_mtid(0);
    high_flags.set_mtid(1);
    low_flags.set_mpos(10);
    high_flags.set_mpos(10);
    assert_eq!(java_read_coordinate_compare(&low_flags, &high_flags), -1);

    high_flags.set_mtid(0);
    low_flags.set_mpos(10);
    high_flags.set_mpos(40);
    assert_eq!(java_read_coordinate_compare(&low_flags, &high_flags), -1);

    high_flags.set_mpos(10);
    low_flags.set_insert_size(10);
    high_flags.set_insert_size(50);
    assert_eq!(java_read_coordinate_compare(&low_flags, &high_flags), -1);
    antisym(&low_flags, &high_flags);
    assert_eq!(java_read_coordinate_compare(&low_flags, &low_flags), 0);
}

fn load_java(text: &str) -> (Vec<Key>, HashMap<Key, [f64; 3]>) {
    let mut order = Vec::new();
    let mut alleles = HashMap::new();
    for line in text.lines() {
        let p: Vec<_> = line.split('\t').collect();
        if p.len() >= 7 && p[0] == "6R297" && p[1] == "allele_ll" {
            let key = (p[2].to_string(), p[3].parse().unwrap());
            let vals: Vec<f64> = p[4..].iter().copied().map(cell_bits).collect();
            order.push(key.clone());
            alleles.insert(key, [vals[0], vals[1], vals[2]]);
        }
    }
    (order, alleles)
}

fn arm(on: bool) {
    gatk_haplotypecaller::hc_genotyping_engine::set_forensic_6r287_global_hap_floor(None, None);
    gatk_haplotypecaller::hc_genotyping_engine::set_forensic_6r288_haplotype_substitute(
        None, None, 0, None,
    );
    gatk_haplotypecaller::likelihood_engine::set_forensic_6r289_target(None, 0);
    gatk_haplotypecaller::likelihood_engine::set_forensic_6r289_substitute(None, None);
    if on {
        set_forensic_6r294_padded_span(Some(JAVA_SPAN.0), Some(JAVA_SPAN.1));
        set_forensic_6r306_exclude_failed_mate(true);
        set_forensic_6r308_java_read_coordinate_order(true);
        gatk_haplotypecaller::likelihood_engine::set_forensic_6r289_target(
            Some(DIAG_QNAME),
            DIAG_FLAGS,
        );
    } else {
        set_forensic_6r294_padded_span(None, None);
        set_forensic_6r306_exclude_failed_mate(false);
        set_forensic_6r308_java_read_coordinate_order(false);
    }
}

fn ulp_distance(a: f64, b: f64) -> u64 {
    if a.to_bits() == b.to_bits() {
        return 0;
    }
    fn ordered(x: f64) -> u64 {
        let bits = x.to_bits();
        if bits & (1 << 63) != 0 {
            !bits
        } else {
            bits | (1 << 63)
        }
    }
    ordered(a).abs_diff(ordered(b))
}

fn rows_gl(order: &[Key], alleles: &HashMap<Key, [f64; 3]>) -> Vec<f64> {
    let rows: Vec<_> = order
        .iter()
        .enumerate()
        .map(|(i, k)| ReadLikelihoodRow {
            read_index: i,
            read_id: String::new(),
            haplotype_log10_likelihoods: alleles[k].to_vec(),
        })
        .collect();
    diploid_genotype_log10_likelihoods_from_allele_rows(&rows, 3)
}

fn continuous_pl(gls: &[f64]) -> Vec<f64> {
    let best = gls.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    gls.iter().map(|g| -10.0 * (g - best)).collect()
}

struct Pass {
    geno: Vec<Key>,
    alleles: HashMap<Key, [f64; 3]>,
    gl: Vec<f64>,
    integer: String,
    subset: String,
    pairhmm_122: Vec<Key>,
    by_key: HashMap<Key, Forensic6r307Read>,
    excluded: usize,
    active: usize,
    passing: usize,
    geometry_ok: bool,
}

fn one_pass(root: &Path, with_geometry: bool) -> Pass {
    arm(true);
    begin_forensic_6r307_order_observe();
    if with_geometry {
        begin_trim_calc_observe();
        begin_hap_list_observe();
    }
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
            ) && r.start.get() <= LOCUS
                && r.end.get() >= LOCUS
        })
        .expect("ActiveFull");
    let filter_cfg = AssemblyReadFilterConfig::gatk_defaults();
    let active = region.reads.len();
    let passing = region
        .reads
        .iter()
        .filter(|r| passes_assembly_read(r, &filter_cfg))
        .count();
    let outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let excluded = forensic_6r306_excluded_count();
    let stages = take_forensic_6r307_order();
    let geometry_ok = if with_geometry {
        let applied = take_hap_list_trim_span().expect("span");
        let snaps = take_hap_list_snaps();
        let read = gatk_haplotypecaller::likelihood_engine::take_forensic_6r289_snaps();
        let hap = snaps
            .iter()
            .find(|s: &&HapListSnap| s.stage == "pairhmm_input")
            .expect("haps");
        let read = read.first().expect("diag");
        let ref_h = hap.columns.iter().find(|c| c.is_reference).expect("ref");
        let lens: BTreeSet<_> = hap.columns.iter().map(|c| c.len).collect();
        let tg: Vec<_> = (19..24)
            .filter_map(|i| hap.columns.get(i).map(|c| c.len))
            .collect();
        applied.trim_start == JAVA_SPAN.0
            && applied.trim_end == JAVA_SPAN.1
            && format!(
                "{} @ {} len={}",
                read.cigar,
                read.pos_1based,
                read.read_bases.len()
            ) == "51H97M @ 29455560 len=97"
            && ref_h.cigar == "169M"
            && ref_h.len == 169
            && lens == BTreeSet::from([169, 170, 172, 174])
            && tg == vec![174, 174, 174, 174, 174]
    } else {
        true
    };
    arm(false);
    let numerics = take_colocated_merge_numerics();
    let site = numerics
        .iter()
        .find(|n| n.loc == LOCUS && n.alts.iter().any(|a| a == ALT))
        .expect("site");
    let col = |name: &str| -> usize {
        if site.long_ref == name {
            0
        } else {
            site.alts.iter().position(|a| a == name).unwrap() + 1
        }
    };
    let cols = [col("T"), col("TTTG"), col(ALT)];
    let mut geno = Vec::new();
    let mut alleles = HashMap::new();
    for i in 0..site.ad_row_qname.len() {
        let key = (site.ad_row_qname[i].clone(), site.ad_row_flags[i]);
        let row = &site.ad_row_lls[i];
        alleles.insert(key.clone(), [row[cols[0]], row[cols[1]], row[cols[2]]]);
        geno.push(key);
    }
    let pairhmm = stages
        .iter()
        .find(|(s, _)| s == "pairhmm_input")
        .expect("pairhmm")
        .1
        .clone();
    let keep: HashSet<Key> = geno.iter().cloned().collect();
    let mut by_key = HashMap::new();
    let mut pairhmm_122 = Vec::new();
    for row in &pairhmm {
        let key = (row.qname.clone(), row.flags);
        if keep.contains(&key) {
            pairhmm_122.push(key.clone());
        }
        by_key.insert(key, row.clone());
    }
    let integer = site
        .merged_pl
        .iter()
        .map(|p| p.to_string())
        .collect::<Vec<_>>()
        .join(",");
    let recs =
        try_emit_call_region_variants(region, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    let emitted = recs
        .iter()
        .find(|r| r.position == LOCUS && r.alternate.iter().any(|a| a == ALT))
        .expect("emit");
    let subset = emitted.samples[0]
        .pl
        .clone()
        .unwrap_or_default()
        .iter()
        .map(|p| p.to_string())
        .collect::<Vec<_>>()
        .join(",");
    Pass {
        geno,
        alleles,
        gl: site.merged_gls.clone(),
        integer,
        subset,
        pairhmm_122,
        by_key,
        excluded,
        active,
        passing,
        geometry_ok,
    }
}

#[test]
fn forensic_6r308_java_read_coordinate_comparator() {
    let root = repo_root();
    let text = std::fs::read_to_string(root.join(JAVA_TSV)).unwrap();
    let (java, java_alleles) = load_java(&text);
    let java_gl = field_cells(&text, "java_gl_vector");
    let java_pl = field_cells(&text, "java_continuous_pl_vector");
    assert_eq!(java.len(), 122);

    let a = one_pass(&root, true);
    let b = one_pass(&root, false);
    assert!(a.geometry_ok);
    assert_eq!(a.geno, b.geno);
    assert_eq!(a.active, 450);
    assert_eq!(a.passing, 449);
    assert_eq!(a.excluded, 1);
    assert!(!a
        .geno
        .iter()
        .any(|k| k.0 == MATE_QNAME && k.1 == MATE_FLAGS));
    assert_eq!(a.geno.len(), 122);
    assert_eq!(a.geno, java);
    assert_eq!(a.pairhmm_122, java);

    let fwd = &a.by_key[&(FORWARD_QNAME.to_string(), FORWARD_FLAGS)];
    let rev = &a.by_key[&(REVERSE_QNAME.to_string(), REVERSE_FLAGS)];
    assert_eq!(fwd.start_1based, 29_455_560);
    assert_eq!(rev.start_1based, 29_455_560);
    let pair = java_read_coordinate_compare_reads(fwd, rev);
    kv("compare_forward_vs_reverse", pair.to_string());
    assert!(pair < 0);
    assert_eq!(a.geno[0].0, FORWARD_QNAME);
    assert_eq!(a.geno[0].1, FORWARD_FLAGS);
    assert_eq!(a.geno[1].0, "HISEQ1:11:H8GV6ADXX:1:2106:20088:100260");
    assert_eq!(a.geno[1].1, 163);
    assert_eq!(a.geno[2].0, DIAG_QNAME);

    let mut compared = 0usize;
    let mut differing = 0usize;
    let mut max_abs = 0.0_f64;
    let mut max_ulp = 0u64;
    for key in &java {
        let rust = a.alleles[key];
        let java_row = java_alleles[key];
        for i in 0..3 {
            compared += 1;
            let abs = (rust[i] - java_row[i]).abs();
            let ulp = ulp_distance(rust[i], java_row[i]);
            max_abs = max_abs.max(abs);
            if ulp != u64::MAX {
                max_ulp = max_ulp.max(ulp);
            }
            if rust[i].to_bits() != java_row[i].to_bits() {
                differing += 1;
            }
        }
    }
    assert_eq!(compared, 366);
    assert_eq!(differing, 0);
    assert_eq!(max_ulp, 0);

    let gl_delta: Vec<f64> =
        a.gl.iter()
            .zip(java_gl.iter())
            .map(|(r, j)| r - j)
            .collect();
    let pl = continuous_pl(&a.gl);
    let pl_delta: Vec<f64> = pl.iter().zip(java_pl.iter()).map(|(r, j)| r - j).collect();
    assert!(a
        .gl
        .iter()
        .zip(java_gl.iter())
        .all(|(r, j)| r.to_bits() == j.to_bits()));
    assert!(pl
        .iter()
        .zip(java_pl.iter())
        .all(|(r, j)| r.to_bits() == j.to_bits()));
    assert_eq!(a.integer, "570,148,3485,0,2762,3517");
    assert_eq!(a.subset, "570,0,3517");

    let mut prefix: Vec<Key> = a.geno.clone();
    prefix.sort_by(|x, y| {
        let p = &a.by_key[x];
        let q = &a.by_key[y];
        p.start_1based
            .cmp(&q.start_1based)
            .then_with(|| p.qname.cmp(&q.qname))
    });
    let prefix_gl = rows_gl(&prefix, &a.alleles);
    let java_replay = rows_gl(&java, &a.alleles);
    let live_replay = rows_gl(&a.geno, &a.alleles);
    assert!(java_replay
        .iter()
        .zip(java_gl.iter())
        .all(|(r, j)| r.to_bits() == j.to_bits()));
    assert!(live_replay
        .iter()
        .zip(a.gl.iter())
        .all(|(r, j)| r.to_bits() == j.to_bits()));
    assert_ne!(prefix, java);
    assert_eq!(ulp_distance(prefix_gl[3], java_gl[3]), 2);
    assert_eq!(ulp_distance(prefix_gl[5], java_gl[5]), 1);

    for i in 0..20 {
        kv(
            "pos",
            format!(
                "{i}\t{}\t{}\t{}\t{}",
                java[i].0, java[i].1, a.geno[i].0, a.geno[i].1
            ),
        );
    }
    kv("positional_differences", "0");
    kv(
        "gl_rust_minus_java",
        gl_delta
            .iter()
            .map(|d| fmt_f(*d))
            .collect::<Vec<_>>()
            .join(","),
    );
    kv(
        "pl_rust_minus_java",
        pl_delta
            .iter()
            .map(|d| fmt_f(*d))
            .collect::<Vec<_>>()
            .join(","),
    );
    kv("integer_pl", &a.integer);
    kv("emitted_subset_pl", &a.subset);
    kv("deterministic", "YES");
    kv(
        "classification",
        "JAVA_READ_COORDINATE_COMPARATOR_REPRODUCES_JAVA",
    );
}

fn java_read_coordinate_compare_reads(a: &Forensic6r307Read, b: &Forensic6r307Read) -> i32 {
    let mut left = rec(&a.qname, a.flags, a.start_1based);
    let mut right = rec(&b.qname, b.flags, b.start_1based);
    left.set_mapq(a.mapq);
    right.set_mapq(b.mapq);
    left.set_mtid(a.mate_tid);
    right.set_mtid(b.mate_tid);
    if a.mate_start_1based > 0 {
        left.set_mpos(a.mate_start_1based - 1);
    }
    if b.mate_start_1based > 0 {
        right.set_mpos(b.mate_start_1based - 1);
    }
    left.set_insert_size(a.tlen);
    right.set_insert_size(b.tlen);
    java_read_coordinate_compare(&left, &right)
}
