//! 6R.98 holdout: Java default float vs Java double vs Rust NeonF64 primitive buffers.
//!
//! Skipped unless `HOLDOUT_6R98=1`. Coordinate-free contract lives in
//! `forensic_6r98_java_double_precision_result_contract`.
//!
//! ```text
//! HOLDOUT_6R98=1 cargo test -p gatk-haplotypecaller --test holdout_6r98_java_double_precision -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::{
    begin_likelihood_pipeline_observe, call_disposition, flatten_assembly_regions,
    take_likelihood_pipeline_cells, take_likelihood_pipeline_snaps,
    traverse_assembly_region_walker, try_emit_call_region_variants, AssemblyRegionCallDisposition,
    CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
    DEFAULT_STAND_EMIT_CONFIDENCE,
};
use serde_json::json;
use std::collections::{BTreeSet, HashMap, HashSet};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const POS_SNP: u64 = 29_456_344;
const JAVA_FLOAT: &str = include_str!("6r96_java_seq6_post_kernel.tsv");
const JAVA_DOUBLE: &str = include_str!("6r98_java_seq6_double_post_kernel.tsv");
const JAVA_ONLY: [u64; 2] = [0xfa2d2442dde7f8ff, 0x48eb4b18de00d4fd];

const READS: &[(&str, u16)] = &[
    ("HISEQ1:11:H8GV6ADXX:2:2216:2203:76921", 147),
    ("HISEQ1:13:H8G92ADXX:1:1111:12251:89078", 83),
    ("HISEQ1:9:H8962ADXX:1:1112:19265:60083", 99),
    ("HWI-D00360:5:H814YADXX:1:1202:11051:34179", 147),
    ("HWI-D00360:5:H814YADXX:1:2207:10890:76583", 147),
    ("HWI-D00360:5:H814YADXX:2:1102:2154:52493", 163),
    ("HWI-D00360:6:H81VLADXX:2:1104:15554:2818", 83),
    ("HWI-D00360:6:H81VLADXX:2:1202:18367:85709", 163),
    ("HWI-D00360:7:H88WKADXX:1:1116:9273:30844", 83),
    ("HWI-D00360:8:H88U0ADXX:1:2108:16806:75328", 163),
    ("HISEQ1:13:H8G92ADXX:1:1205:16330:83279", 163),
    ("HISEQ1:9:H8962ADXX:2:1212:17767:73796", 83),
    ("HWI-D00360:5:H814YADXX:2:2103:4936:45407", 83),
    ("HWI-D00360:6:H81VLADXX:1:1103:1948:22968", 147),
    ("HWI-D00360:6:H81VLADXX:1:1210:4156:72506", 83),
    ("HWI-D00360:7:H88WKADXX:1:2111:4466:65743", 147),
    ("HWI-D00360:7:H88WKADXX:1:2203:20480:101193", 163),
    ("HWI-D00360:7:H88WKADXX:2:1214:6938:52704", 83),
    ("HWI-D00360:8:H88U0ADXX:1:1205:11075:4786", 147),
    ("HWI-D00360:8:H88U0ADXX:1:1213:18559:65935", 163),
    ("HWI-D00360:8:H88U0ADXX:1:2213:15618:11579", 163),
    ("HWI-D00360:8:H88U0ADXX:2:1213:15376:17578", 163),
];

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn parse_kv(s: &str) -> HashMap<&str, &str> {
    let mut d = HashMap::new();
    for part in s.split('\t') {
        if let Some((k, v)) = part.split_once('=') {
            d.insert(k, v);
        }
    }
    d
}

fn is_exact_f32_widened(x: f64) -> bool {
    (f64::from(x as f32)).to_bits() == x.to_bits()
}

fn ordered_f64_bits(x: f64) -> i64 {
    let bits = x.to_bits() as i64;
    if bits < 0 {
        i64::MIN - bits
    } else {
        bits
    }
}

fn ulp_distance(a: f64, b: f64) -> u64 {
    ordered_f64_bits(a).abs_diff(ordered_f64_bits(b))
}

struct JavaDump {
    hap_fnv: Vec<u64>,
    bits: HashMap<(String, u16), Vec<u64>>,
}

fn parse_java_dump(text: &str, prefix: &str) -> JavaDump {
    let mut hap_fnv = Vec::new();
    let mut bits = HashMap::new();
    for line in text.lines() {
        let p: Vec<&str> = line.split('\t').collect();
        if p.len() < 4 || p[0] != prefix {
            continue;
        }
        let key = p[2];
        let val = p[3..].join("\t");
        let kv = parse_kv(&val);
        if key.starts_with("hap_") {
            let fnv = u64::from_str_radix(kv.get("fnv").copied().unwrap_or("0"), 16).unwrap_or(0);
            hap_fnv.push(fnv);
        } else if key.starts_with("rowbits_") {
            let q = kv.get("qname").copied().unwrap_or("").to_string();
            let flags: u16 = kv.get("flags").copied().unwrap_or("0").parse().unwrap_or(0);
            let b = kv
                .get("bits")
                .copied()
                .unwrap_or("")
                .split(',')
                .filter(|s| !s.is_empty())
                .map(|s| u64::from_str_radix(s, 16).unwrap_or(0))
                .collect();
            bits.insert((q, flags), b);
        }
    }
    JavaDump { hap_fnv, bits }
}

fn cells_from_dump(dump: &JavaDump, reads: &[(&str, u16)]) -> HashMap<(String, u16, u64), f64> {
    let mut map = HashMap::new();
    for &(q, flags) in reads.iter().take(22) {
        let bits = dump
            .bits
            .get(&(q.to_string(), flags))
            .unwrap_or_else(|| panic!("missing {q} flags={flags}"));
        assert_eq!(bits.len(), dump.hap_fnv.len());
        for (col, &fnv) in dump.hap_fnv.iter().enumerate() {
            map.insert((q.to_string(), flags, fnv), f64::from_bits(bits[col]));
        }
    }
    map
}

struct PairStats {
    n: usize,
    eq: usize,
    diff: usize,
    max_abs: f64,
    mean_abs: f64,
    median_abs: f64,
    max_ulp: u64,
    mean_ulp: f64,
    f32_wide_a: usize,
    f32_wide_b: usize,
}

fn pair_stats(
    a: &HashMap<(String, u16, u64), f64>,
    b: &HashMap<(String, u16, u64), f64>,
    reads: &[(&str, u16)],
    common: &BTreeSet<u64>,
) -> PairStats {
    let mut n = 0usize;
    let mut eq = 0usize;
    let mut diff = 0usize;
    let mut absds = Vec::new();
    let mut ulps = Vec::new();
    let mut f32_a = 0usize;
    let mut f32_b = 0usize;
    for &(q, flags) in reads.iter().take(22) {
        for &h in common {
            let k = (q.to_string(), flags, h);
            let Some(&av) = a.get(&k) else {
                continue;
            };
            let Some(&bv) = b.get(&k) else {
                continue;
            };
            n += 1;
            if is_exact_f32_widened(av) {
                f32_a += 1;
            }
            if is_exact_f32_widened(bv) {
                f32_b += 1;
            }
            if av.to_bits() == bv.to_bits() {
                eq += 1;
            } else {
                diff += 1;
            }
            absds.push((av - bv).abs());
            ulps.push(ulp_distance(av, bv));
        }
    }
    absds.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
    let median = if absds.is_empty() {
        0.0
    } else {
        absds[absds.len() / 2]
    };
    PairStats {
        n,
        eq,
        diff,
        max_abs: absds.iter().copied().fold(0.0, f64::max),
        mean_abs: if n == 0 {
            0.0
        } else {
            absds.iter().sum::<f64>() / n as f64
        },
        median_abs: median,
        max_ulp: ulps.iter().copied().max().unwrap_or(0),
        mean_ulp: if n == 0 {
            0.0
        } else {
            ulps.iter().sum::<u64>() as f64 / n as f64
        },
        f32_wide_a: f32_a,
        f32_wide_b: f32_b,
    }
}

fn rust_post(
    cells: &[gatk_haplotypecaller::LikelihoodPipelineCell],
) -> (HashMap<(String, u16, u64), f64>, usize, usize, HashSet<u64>) {
    let seq = cells
        .iter()
        .filter(|c| c.stage == "post_kernel")
        .map(|c| c.seq)
        .min()
        .unwrap_or(0);
    let slice: Vec<_> = cells
        .iter()
        .filter(|c| c.stage == "post_kernel" && c.seq == seq)
        .collect();
    let mut map = HashMap::new();
    let mut haps = HashSet::new();
    let mut n_reads = 0usize;
    let mut n_haps = 0usize;
    for c in &slice {
        map.insert((c.qname.clone(), c.flags, c.hap_fnv), c.log10_likelihood);
        haps.insert(c.hap_fnv);
        n_reads = c.n_reads;
        n_haps = c.n_haps;
    }
    (map, n_reads, n_haps, haps)
}

fn classify(fd: &PairStats, dr: &PairStats) -> &'static str {
    if dr.eq == dr.n && dr.n > 0 {
        "JAVA_DOUBLE_MATCHES_RUST"
    } else if fd.eq == fd.n && fd.n > 0 && fd.f32_wide_b == fd.n {
        "PRECISION_SWITCH_NOT_ACTIVE"
    } else if fd.diff == fd.n && dr.diff == dr.n && fd.f32_wide_b == 0 {
        "JAVA_DOUBLE_STILL_DIVERGES"
    } else if dr.eq == dr.n {
        "JAVA_FLOAT_RESULT_PATH"
    } else {
        "NO_PROVEN_RELATIONSHIP"
    }
}

#[test]
fn holdout_6r98_java_double_precision_three_way() {
    if std::env::var("HOLDOUT_6R98").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R98=1");
        return;
    }
    let java_f = parse_java_dump(JAVA_FLOAT, "6R96");
    let java_d = parse_java_dump(JAVA_DOUBLE, "6R98");
    assert_eq!(java_f.hap_fnv, java_d.hap_fnv);
    assert_eq!(java_f.hap_fnv.len(), 70);

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

    begin_likelihood_pipeline_observe();
    let outcome = HaplotypeCallerEngine::call_region(
        covering[0],
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let cells = take_likelihood_pipeline_cells();
    let _snaps = take_likelihood_pipeline_snaps();
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
        .expect("canonical T/C");

    let (rust, rust_r, rust_c, rust_haps) = rust_post(&cells);
    let java_only: HashSet<u64> = JAVA_ONLY.into_iter().collect();
    let java_set: HashSet<u64> = java_f.hap_fnv.iter().copied().collect();
    let common: BTreeSet<u64> = java_set
        .intersection(&rust_haps)
        .copied()
        .filter(|h| !java_only.contains(h))
        .collect();

    let jf = cells_from_dump(&java_f, READS);
    let jd = cells_from_dump(&java_d, READS);
    let fd = pair_stats(&jf, &jd, READS, &common);
    let dr = pair_stats(&jd, &rust, READS, &common);
    let fr = pair_stats(&jf, &rust, READS, &common);
    let classification = classify(&fd, &dr);

    let first_hap = *common.iter().min().unwrap_or(&0);
    let live = READS[0];
    let k = (live.0.to_string(), live.1, first_hap);
    eprintln!(
        "{}",
        json!({
            "classification": classification,
            "common_haps": common.len(),
            "rust_rows": rust_r,
            "rust_cols": rust_c,
            "j_float_vs_j_double": {
                "n": fd.n, "eq": fd.eq, "diff": fd.diff,
                "max_abs": fd.max_abs, "mean_abs": fd.mean_abs, "median_abs": fd.median_abs,
                "max_ulp": fd.max_ulp, "mean_ulp": fd.mean_ulp,
                "f32_wide_float": fd.f32_wide_a, "f32_wide_double": fd.f32_wide_b,
            },
            "j_double_vs_r_f64": {
                "n": dr.n, "eq": dr.eq, "diff": dr.diff,
                "max_abs": dr.max_abs, "mean_abs": dr.mean_abs, "median_abs": dr.median_abs,
                "max_ulp": dr.max_ulp, "mean_ulp": dr.mean_ulp,
            },
            "j_float_vs_r_f64": {
                "n": fr.n, "eq": fr.eq, "diff": fr.diff,
                "max_abs": fr.max_abs, "mean_abs": fr.mean_abs, "median_abs": fr.median_abs,
                "max_ulp": fr.max_ulp, "mean_ulp": fr.mean_ulp,
            },
            "first_hap": format!("{:x}", first_hap),
            "first_j_float": jf.get(&k),
            "first_j_double": jd.get(&k),
            "first_r_f64": rust.get(&k),
            "first_j_float_bits": jf.get(&k).map(|v| format!("{:016x}", v.to_bits())),
            "first_j_double_bits": jd.get(&k).map(|v| format!("{:016x}", v.to_bits())),
            "first_r_bits": rust.get(&k).map(|v| format!("{:016x}", v.to_bits())),
            "vcf_ad": vcf.samples.first().map(|s| s.ad.clone()),
        })
    );

    assert_eq!(fd.n, 1496);
    assert_eq!(dr.n, 1496);
    assert_eq!(fr.n, 1496);
    assert_eq!(fd.eq, 0);
    assert_eq!(fd.diff, 1496);
    assert_eq!(fr.eq, 0);
    assert_eq!(fr.diff, 1496);
    assert_eq!(fd.f32_wide_a, 1496);
    assert_eq!(fd.f32_wide_b, 0);
    assert_eq!(rust_r, 153);
    assert_eq!(rust_c, 70);
    assert_eq!(classification, "JAVA_DOUBLE_STILL_DIVERGES");
    assert_eq!(dr.eq, 0);
    assert_eq!(dr.diff, 1496);
    let _ = outcome;
}
