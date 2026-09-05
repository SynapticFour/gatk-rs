//! 6R.95 holdout: common-68 likelihood cells at poorly-modeled (live vs Java dump).
//!
//! Skipped unless `HOLDOUT_6R95=1`. Coordinate-free contract lives in
//! `forensic_6r95_common_haplotype_likelihood_contract`.
//!
//! Java cells: `6r95_java_seq6_filter_object.tsv` from GATK 4.4.0.0
//! `filterPoorlyModeledEvidence` seq=6 (`20:29456294-29456500`).
//!
//! ```text
//! HOLDOUT_6R95=1 cargo test -p gatk-haplotypecaller --test holdout_6r95_common_haplotype_likelihood -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::{
    begin_poorly_modeled_observe, call_disposition, flatten_assembly_regions,
    take_poorly_modeled_cells, take_poorly_modeled_haplotypes, take_poorly_modeled_observe,
    traverse_assembly_region_walker, try_emit_call_region_variants, AssemblyRegionCallDisposition,
    CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
    DEFAULT_STAND_EMIT_CONFIDENCE,
};
use serde_json::json;
use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const POS_SNP: u64 = 29_456_344;
const JAVA_DUMP: &str = include_str!("6r95_java_seq6_filter_object.tsv");

const JAVA_ONLY: [u64; 2] = [0xfa2d2442dde7f8ff, 0x48eb4b18de00d4fd];
const JAVA_LIVE_ONLY_N: usize = 10;
const RUST_ONLY_N: usize = 12;

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
    ("HISEQ1:11:H8GV6ADXX:2:2105:12137:22761", 163),
    ("HISEQ1:13:H8G92ADXX:1:1201:11859:45984", 147),
];

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn tag(i: usize) -> &'static str {
    if i < JAVA_LIVE_ONLY_N {
        "JAVA_LIVE_ONLY"
    } else if i < JAVA_LIVE_ONLY_N + RUST_ONLY_N {
        "RUST_ONLY"
    } else {
        "BOTH_DROP"
    }
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

struct JavaDump {
    hap_fnv: Vec<u64>,
    row_index: HashMap<(String, u16), usize>,
    bits: HashMap<(String, u16), Vec<u64>>,
}

fn parse_java_dump(text: &str) -> JavaDump {
    let mut hap_fnv = Vec::new();
    let mut row_index = HashMap::new();
    let mut bits = HashMap::new();
    for line in text.lines() {
        let p: Vec<&str> = line.split('\t').collect();
        if p.len() < 4 {
            continue;
        }
        let key = p[2];
        let val = p[3..].join("\t");
        let kv = parse_kv(&val);
        if key.starts_with("hap_") {
            let fnv = u64::from_str_radix(kv.get("fnv").copied().unwrap_or("0"), 16).unwrap_or(0);
            hap_fnv.push(fnv);
        } else if key.starts_with("row_") && !key.starts_with("rowbits_") {
            let q = kv.get("qname").copied().unwrap_or("").to_string();
            let flags: u16 = kv.get("flags").copied().unwrap_or("0").parse().unwrap_or(0);
            let row: usize = kv.get("row").copied().unwrap_or("0").parse().unwrap_or(0);
            row_index.insert((q, flags), row);
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
    JavaDump {
        hap_fnv,
        row_index,
        bits,
    }
}

fn max_and_argmax(vals: &[(u64, f64)]) -> (f64, u64) {
    let mut best = f64::NEG_INFINITY;
    let mut arg = 0u64;
    for &(h, v) in vals {
        if v > best {
            best = v;
            arg = h;
        }
    }
    (best, arg)
}

#[test]
fn holdout_6r95_common_68_likelihood_cells() {
    if std::env::var("HOLDOUT_6R95").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R95=1");
        return;
    }
    let java = parse_java_dump(JAVA_DUMP);
    assert_eq!(java.hap_fnv.len(), 70);
    assert_eq!(java.bits.len(), 153);
    let java_only: HashSet<u64> = JAVA_ONLY.into_iter().collect();
    assert_eq!(
        java.hap_fnv
            .iter()
            .copied()
            .filter(|h| java_only.contains(h))
            .collect::<HashSet<_>>(),
        java_only
    );

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

    begin_poorly_modeled_observe();
    let outcome = HaplotypeCallerEngine::call_region(
        covering[0],
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let observed = take_poorly_modeled_observe();
    let hap_obs = take_poorly_modeled_haplotypes();
    let cells = take_poorly_modeled_cells();
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

    let last_pass = observed.iter().map(|r| r.pass).max().unwrap_or(0);
    let last: Vec<_> = observed.iter().filter(|r| r.pass == last_pass).collect();
    let rust_rows: HashMap<(String, u16), usize> = last
        .iter()
        .map(|r| ((r.qname.clone(), r.flags), r.row_index))
        .collect();
    let rust_haps: Vec<u64> = hap_obs
        .iter()
        .filter(|h| h.pass == last_pass)
        .map(|h| h.fnv1a)
        .collect();
    let java_set: HashSet<u64> = java.hap_fnv.iter().copied().collect();
    let rust_set: HashSet<u64> = rust_haps.iter().copied().collect();
    let common: HashSet<u64> = java_set.intersection(&rust_set).copied().collect();
    let rust_only: HashSet<u64> = rust_set.difference(&java_set).copied().collect();
    assert_eq!(common.len(), 68);
    assert_eq!(rust_haps.len(), 68);
    assert!(rust_only.is_empty());
    assert_eq!(last.len(), 153);

    let rust_cell: HashMap<(String, u16, u64), f64> = cells
        .iter()
        .filter(|c| c.pass == last_pass)
        .map(|c| ((c.qname.clone(), c.flags, c.hap_fnv), c.log10_likelihood))
        .collect();

    eprintln!(
        "6R.95 java={}x{} rust={}x{} common={} JAVA_ONLY=2 RUST_ONLY={}",
        java.bits.len(),
        java.hap_fnv.len(),
        last.len(),
        rust_haps.len(),
        common.len(),
        rust_only.len()
    );
    eprintln!(
        "TAG\tQNAME\tJrow\tRrow\tJmax70\tJmax68\tRmax68\tsame68max\tJ70arg\tJ68arg\tR68arg\tncmp\tndiff\tmaxabs\tmeanabs"
    );

    let mut n_cells = 0usize;
    let mut n_diff = 0usize;
    let mut n_eq = 0usize;
    let mut sum_abs = 0.0f64;
    let mut max_abs = 0.0f64;
    let mut n_same68 = 0usize;
    let mut n_java70_gt_68 = 0usize;
    let mut n_java70_eq_68 = 0usize;
    let mut n_flips_explained_by_java_only = 0usize;

    for (i, &(q, flags)) in READS.iter().enumerate() {
        if i >= 22 {
            continue;
        }
        let key = (q.to_string(), flags);
        let jrow = *java.row_index.get(&key).expect("java row");
        let rrow = *rust_rows.get(&key).expect("rust row");
        let jbits = java.bits.get(&key).expect("java bits");
        assert_eq!(jbits.len(), 70);
        let mut j70 = Vec::with_capacity(70);
        let mut j68 = Vec::with_capacity(68);
        let mut r68 = Vec::with_capacity(68);
        let mut nd = 0usize;
        let mut ncmp = 0usize;
        let mut sum = 0.0;
        let mut mxd = 0.0;
        for (col, &fnv) in java.hap_fnv.iter().enumerate() {
            let jv = f64::from_bits(jbits[col]);
            j70.push((fnv, jv));
            if !common.contains(&fnv) {
                continue;
            }
            ncmp += 1;
            n_cells += 1;
            let rv = rust_cell.get(&(q.to_string(), flags, fnv)).copied();
            match rv {
                Some(rv) => {
                    j68.push((fnv, jv));
                    r68.push((fnv, rv));
                    if jv.to_bits() == rv.to_bits() {
                        n_eq += 1;
                    } else {
                        nd += 1;
                        n_diff += 1;
                    }
                    let d = (jv - rv).abs();
                    sum += d;
                    sum_abs += d;
                    if d > mxd {
                        mxd = d;
                    }
                    if d > max_abs {
                        max_abs = d;
                    }
                }
                None => {
                    nd += 1;
                    n_diff += 1;
                    j68.push((fnv, jv));
                }
            }
        }
        let (jmax70, jarg70) = max_and_argmax(&j70);
        let (jmax68, jarg68) = max_and_argmax(&j68);
        let (rmax68, rarg68) = max_and_argmax(&r68);
        let same68 = jmax68.to_bits() == rmax68.to_bits();
        if same68 {
            n_same68 += 1;
        }
        if jmax70.to_bits() == jmax68.to_bits() {
            n_java70_eq_68 += 1;
        } else if jmax70 > jmax68 {
            n_java70_gt_68 += 1;
        }
        let jkeep70 = jmax70 >= -8.0;
        let jkeep68 = jmax68 >= -8.0;
        let rkeep68 = rmax68 >= -8.0;
        if jkeep70 != rkeep68 && jkeep68 == rkeep68 {
            n_flips_explained_by_java_only += 1;
        }
        let mean = if ncmp == 0 { 0.0 } else { sum / ncmp as f64 };
        eprintln!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:x}\t{:x}\t{:x}\t{}\t{}\t{}\t{}",
            tag(i),
            q,
            jrow,
            rrow,
            jmax70,
            jmax68,
            rmax68,
            same68,
            jarg70,
            jarg68,
            rarg68,
            ncmp,
            nd,
            mxd,
            mean
        );
        assert_eq!(ncmp, 68);
    }

    let mean_abs = if n_cells == 0 {
        0.0
    } else {
        sum_abs / n_cells as f64
    };
    let classification = if n_diff > 0 {
        "PRE_FILTER_LIKELIHOOD_VALUE"
    } else if n_same68 < 22 {
        "MAX_LL_REDUCTION"
    } else {
        "LIKELIHOOD_COLUMN_POPULATION_ONLY"
    };
    let java_only_explain = if n_flips_explained_by_java_only == 22 {
        "YES"
    } else {
        "NO"
    };

    eprintln!(
        "{}",
        json!({
            "classification": classification,
            "common_cells": n_cells,
            "differing_cells": n_diff,
            "exactly_equal": n_eq,
            "max_abs_delta": max_abs,
            "mean_abs_delta": mean_abs,
            "n_same68_max": n_same68,
            "n_java70_eq_68": n_java70_eq_68,
            "n_java70_gt_68": n_java70_gt_68,
            "n_flips_explained_by_java_only": n_flips_explained_by_java_only,
            "java_only_explain_flips": java_only_explain,
            "vcf_ad": vcf.samples[0].ad,
        })
    );
    eprintln!("6R.95 classification={classification} java_only_explain={java_only_explain}");

    assert_eq!(n_cells, 22 * 68);
    assert_eq!(n_diff, 1496);
    assert_eq!(n_eq, 0);
    assert_eq!(n_java70_eq_68, 22);
    assert_eq!(n_flips_explained_by_java_only, 0);
    assert_eq!(classification, "PRE_FILTER_LIKELIHOOD_VALUE");
}
