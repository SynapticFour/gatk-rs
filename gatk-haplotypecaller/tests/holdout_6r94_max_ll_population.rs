//! 6R.94 holdout: live max_ll row/column/argmax attribution at poorly-modeled.
//!
//! Skipped unless `HOLDOUT_6R94=1`. Coordinate-free contract lives in
//! `forensic_6r94_max_ll_population_contract`.
//!
//! Does not assert Rust AD/PL/QUAL. Does not inspect PairHMM kernel arithmetic.
//!
//! ```text
//! HOLDOUT_6R94=1 cargo test -p gatk-haplotypecaller --test holdout_6r94_max_ll_population -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::{
    begin_poorly_modeled_observe, call_disposition, flatten_assembly_regions,
    take_colocated_merge_numerics, take_poorly_modeled_haplotypes, take_poorly_modeled_observe,
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

const JAVA_LIVE_ONLY_N: usize = 10;
const RUST_ONLY_N: usize = 12;

/// Java GATK 4.4.0.0 `filterPoorlyModeledEvidence` seq=6 covering region
/// `20:29456294-29456500`: haplotype-column FNV-1a 64 of bases (70 unique).
const JAVA_HAP_FNV: &[u64] = &[
    0xef45e6aca9014ffd,
    0xccc49a253783bf31,
    0xea70e110402e3bf1,
    0x501e24eda83c4dbd,
    0x85cf71706504b7db,
    0x2e8fd35c0fdb0acf,
    0xac368db7301b1796,
    0x857a39e95c5a92d2,
    0x1e02cf3b25f742ae,
    0xd9dab6ee47f0cd8a,
    0x9f575d82d4ab089c,
    0x51c5482bfeee2058,
    0xf099bba97cf7f4ff,
    0x9f1f32d1edf15f69,
    0xd199f5b1361bdab5,
    0xeabec86dbbe561b9,
    0xd139639be3a02545,
    0x5004750d727fef75,
    0x3e7496f001438729,
    0xfa2d2442dde7f8ff,
    0xc1c6bd92fffa3a4b,
    0xd82092d53640bbcd,
    0x7b572bed34509241,
    0xa6462ff765962201,
    0xe47d5607d425922b,
    0xb1b46e27187945df,
    0xc00bcb2e11106653,
    0x5c748fb6602799e7,
    0x091bc3562a42f13b,
    0xa96dac94841a9aaf,
    0xfdb01091639a3b09,
    0x9b1ebb3aa3acd2d5,
    0x754395d544564246,
    0x2484a2e756a6d042,
    0x3942a66dc14c62e0,
    0xd1a76b206bf16e04,
    0xde70499491a88f58,
    0x2906a7b3bfd0acab,
    0xf63dbfd30424605f,
    0xa3f513a51e07057e,
    0x1c9b32f38c58b60b,
    0xcd49423e5a6000da,
    0xf67585e483f200c2,
    0x473478d271a172c6,
    0x6daa42dccac5e33c,
    0xa10be7920919fbfa,
    0x1a8f682df141de2f,
    0x16baa899e50d6b9e,
    0x1f87fc79947b6c0c,
    0x464e661c6d90f0ca,
    0x1539e5b00300965f,
    0xf250501a49babdee,
    0xd377923e0a5a286e,
    0x7a003b854891edb8,
    0x48eb4b18de00d4fd,
    0x273373fc53c8617c,
    0xa940f3f22a5fb678,
    0x6e24cdcf3fac039e,
    0xf8760cc763b893fa,
    0x1eb2931d452797b0,
    0x7c6040e7a654e1b6,
    0x6c7d2390bb99f7f2,
    0xde4cc7844e4759a2,
    0xfde0a00ce139e404,
    0x0225f0666632f2e1,
    0x657bdb5a3694d8e0,
    0xe648106add71ca97,
    0x80366a53c7c41f83,
    0xa721bf5f47a71821,
    0xd97975790509bf2d,
];

/// qname, flags, row, start, end, cigar, max_ll, argmax_col, argmax_fnv, keep.
const JAVA: &[(&str, u16, usize, i64, i64, &str, f64, usize, u64, bool)] = &[
    (
        "HISEQ1:11:H8GV6ADXX:2:2216:2203:76921",
        147,
        32,
        29456268,
        29456398,
        "18H76M1D54M",
        -7.450010299682617,
        0,
        0xef45e6aca9014ffd,
        true,
    ),
    (
        "HISEQ1:13:H8G92ADXX:1:1111:12251:89078",
        83,
        98,
        29456301,
        29456448,
        "148M",
        -2.5224151611328125,
        12,
        0xf099bba97cf7f4ff,
        true,
    ),
    (
        "HISEQ1:9:H8962ADXX:1:1112:19265:60083",
        99,
        8,
        29456268,
        29456343,
        "72H76M",
        -2.4746932983398438,
        36,
        0xde70499491a88f58,
        true,
    ),
    (
        "HWI-D00360:5:H814YADXX:1:1202:11051:34179",
        147,
        52,
        29456268,
        29456363,
        "52H96M",
        -2.5021400451660156,
        32,
        0x754395d544564246,
        true,
    ),
    (
        "HWI-D00360:5:H814YADXX:1:2207:10890:76583",
        147,
        104,
        29456319,
        29456466,
        "148M",
        -2.6972694396972656,
        13,
        0x9f1f32d1edf15f69,
        true,
    ),
    (
        "HWI-D00360:5:H814YADXX:2:1102:2154:52493",
        163,
        109,
        29456335,
        29456482,
        "148M",
        -2.8041610717773438,
        19,
        0xfa2d2442dde7f8ff,
        true,
    ),
    (
        "HWI-D00360:6:H81VLADXX:2:1104:15554:2818",
        83,
        65,
        29456268,
        29456410,
        "5H143M",
        -2.522064208984375,
        52,
        0xd377923e0a5a286e,
        true,
    ),
    (
        "HWI-D00360:6:H81VLADXX:2:1202:18367:85709",
        163,
        16,
        29456268,
        29456401,
        "14H134M",
        -2.6402854919433594,
        36,
        0xde70499491a88f58,
        true,
    ),
    (
        "HWI-D00360:7:H88WKADXX:1:1116:9273:30844",
        83,
        66,
        29456268,
        29456399,
        "16H132M",
        -2.5336990356445312,
        56,
        0xa940f3f22a5fb678,
        true,
    ),
    (
        "HWI-D00360:8:H88U0ADXX:1:2108:16806:75328",
        163,
        21,
        29456268,
        29456379,
        "36H112M",
        -2.5473976135253906,
        56,
        0xa940f3f22a5fb678,
        true,
    ),
    (
        "HISEQ1:13:H8G92ADXX:1:1205:16330:83279",
        163,
        87,
        29456284,
        29456431,
        "148M",
        -9.549158096313477,
        19,
        0xfa2d2442dde7f8ff,
        false,
    ),
    (
        "HISEQ1:9:H8962ADXX:2:1212:17767:73796",
        83,
        96,
        29456296,
        29456443,
        "148M",
        -9.387346267700195,
        26,
        0xc00bcb2e11106653,
        false,
    ),
    (
        "HWI-D00360:5:H814YADXX:2:2103:4936:45407",
        83,
        57,
        29456268,
        29456351,
        "64H84M",
        -9.047658920288086,
        12,
        0xf099bba97cf7f4ff,
        false,
    ),
    (
        "HWI-D00360:6:H81VLADXX:1:1103:1948:22968",
        147,
        86,
        29456283,
        29456430,
        "148M",
        -12.48884391784668,
        34,
        0x3942a66dc14c62e0,
        false,
    ),
    (
        "HWI-D00360:6:H81VLADXX:1:1210:4156:72506",
        83,
        62,
        29456268,
        29456360,
        "56H76M1D16M",
        -9.251577377319336,
        0,
        0xef45e6aca9014ffd,
        false,
    ),
    (
        "HWI-D00360:7:H88WKADXX:1:2111:4466:65743",
        147,
        68,
        29456268,
        29456411,
        "4H144M",
        -8.373727798461914,
        34,
        0x3942a66dc14c62e0,
        false,
    ),
    (
        "HWI-D00360:7:H88WKADXX:1:2203:20480:101193",
        163,
        92,
        29456290,
        29456438,
        "54M1D94M",
        -13.432897567749023,
        2,
        0xea70e110402e3bf1,
        false,
    ),
    (
        "HWI-D00360:7:H88WKADXX:2:1214:6938:52704",
        83,
        70,
        29456268,
        29456392,
        "23H125M",
        -8.664169311523438,
        34,
        0x3942a66dc14c62e0,
        false,
    ),
    (
        "HWI-D00360:8:H88U0ADXX:1:1205:11075:4786",
        147,
        97,
        29456296,
        29456443,
        "148M",
        -11.753900527954102,
        34,
        0x3942a66dc14c62e0,
        false,
    ),
    (
        "HWI-D00360:8:H88U0ADXX:1:1213:18559:65935",
        163,
        100,
        29456312,
        29456460,
        "33M1D115M",
        -11.865598678588867,
        6,
        0xac368db7301b1796,
        false,
    ),
    (
        "HWI-D00360:8:H88U0ADXX:1:2213:15618:11579",
        163,
        22,
        29456268,
        29456372,
        "43H105M",
        -9.34518814086914,
        26,
        0xc00bcb2e11106653,
        false,
    ),
    (
        "HWI-D00360:8:H88U0ADXX:2:1213:15376:17578",
        163,
        24,
        29456268,
        29456371,
        "44H104M",
        -18.327007293701172,
        6,
        0xac368db7301b1796,
        false,
    ),
    (
        "HISEQ1:11:H8GV6ADXX:2:2105:12137:22761",
        163,
        1,
        29456268,
        29456377,
        "38H110M",
        -9.151548385620117,
        26,
        0xc00bcb2e11106653,
        false,
    ),
    (
        "HISEQ1:13:H8G92ADXX:1:1201:11859:45984",
        147,
        39,
        29456268,
        29456364,
        "51H97M",
        -9.154520034790039,
        12,
        0xf099bba97cf7f4ff,
        false,
    ),
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

#[test]
fn holdout_6r94_max_ll_row_and_column_population() {
    if std::env::var("HOLDOUT_6R94").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R94=1");
        return;
    }
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
    let _snap = take_colocated_merge_numerics()
        .into_iter()
        .find(|n| n.loc == POS_SNP);

    let last_pass = observed.iter().map(|r| r.pass).max().unwrap_or(0);
    let last: Vec<_> = observed.iter().filter(|r| r.pass == last_pass).collect();
    let last_n = last.len();
    let rust_by: HashMap<(String, u16), _> = last
        .into_iter()
        .map(|r| ((r.qname.clone(), r.flags), r))
        .collect();
    let rust_haps: Vec<_> = hap_obs.iter().filter(|h| h.pass == last_pass).collect();
    let java_cols: HashSet<u64> = JAVA_HAP_FNV.iter().copied().collect();
    let rust_cols: HashSet<u64> = rust_haps.iter().map(|h| h.fnv1a).collect();
    let common: HashSet<u64> = java_cols.intersection(&rust_cols).copied().collect();
    let java_only: HashSet<u64> = java_cols.difference(&rust_cols).copied().collect();
    let rust_only: HashSet<u64> = rust_cols.difference(&java_cols).copied().collect();
    eprintln!("JAVA_ONLY hashes: {}", {
        let mut v: Vec<String> = java_only.iter().map(|h| format!("{h:x}")).collect();
        v.sort();
        v.join(",")
    });

    let java_row_keys: HashSet<(String, u16)> =
        JAVA.iter().map(|&(q, f, ..)| (q.to_string(), f)).collect();
    let rust_row_keys: HashSet<(String, u16)> = JAVA
        .iter()
        .filter_map(|&(q, f, ..)| rust_by.get(&(q.to_string(), f)).map(|_| (q.to_string(), f)))
        .collect();

    eprintln!(
        "6R.94 last_pass={} rust_rows={} rust_cols={} java_cols={} covering_reads={}",
        last_pass,
        last_n,
        rust_cols.len(),
        java_cols.len(),
        covering[0].reads.len()
    );
    eprintln!(
        "column set: common={} JAVA_ONLY={} RUST_ONLY={} same_set={} same_order={}",
        common.len(),
        java_only.len(),
        rust_only.len(),
        java_cols == rust_cols,
        JAVA_HAP_FNV.iter().eq(rust_haps.iter().map(|h| &h.fnv1a))
    );

    eprintln!(
        "TAG\tQNAME\tJ_row\tR_row\tJ_cols\tR_cols\tJ_max\tR_max\tJ_argfnv\tR_argfnv\tsame_argmax\targ_in"
    );

    let mut n_row_absent = 0usize;
    let mut n_same_argmax = 0usize;
    let mut n_argmax_java_only = 0usize;
    let mut n_argmax_rust_only = 0usize;
    let mut n_argmax_common_value_diff = 0usize;
    let mut n_cigar_eq = 0usize;
    let mut n_start_eq = 0usize;

    for (i, &(q, flags, jrow, jstart, _je, jcig, jll, _jcol, jfnv, _jkeep)) in
        JAVA.iter().enumerate()
    {
        let rr = rust_by.get(&(q.to_string(), flags));
        let (rrow, rcols, rll, rfnv) = match rr {
            Some(r) => (
                format!("{}", r.row_index),
                format!("{}", r.n_columns),
                format!("{}", r.max_ll),
                format!("{:x}", r.argmax_fnv),
            ),
            None => (
                "ABSENT".into(),
                "ABSENT".into(),
                "ABSENT".into(),
                "ABSENT".into(),
            ),
        };
        let same_argmax = rr.is_some_and(|r| r.argmax_fnv == jfnv);
        let arg_in = if same_argmax {
            "COMMON"
        } else if rust_only.contains(&jfnv) {
            "J_ARG_RUST_ONLY?"
        } else if java_only.contains(&jfnv) && rr.is_some_and(|r| rust_only.contains(&r.argmax_fnv))
        {
            "BOTH_SIDE_ONLY"
        } else if java_only.contains(&jfnv) {
            "JAVA_ONLY"
        } else if rr.is_some_and(|r| rust_only.contains(&r.argmax_fnv)) {
            "RUST_ONLY"
        } else if rr.is_some_and(|r| common.contains(&r.argmax_fnv) && common.contains(&jfnv)) {
            "COMMON_VALUE"
        } else {
            "OTHER"
        };
        eprintln!(
            "{}\t{}\t{}\t{}\t70\t{}\t{}\t{}\t{:x}\t{}\t{}\t{}",
            tag(i),
            q,
            jrow,
            rrow,
            rcols,
            jll,
            rll,
            jfnv,
            rfnv,
            same_argmax,
            arg_in
        );
        match rr {
            None => n_row_absent += 1,
            Some(r) => {
                if r.cigar == jcig {
                    n_cigar_eq += 1;
                }
                if r.start_1based == jstart {
                    n_start_eq += 1;
                }
                if r.argmax_fnv == jfnv {
                    n_same_argmax += 1;
                    if (r.max_ll - jll).abs() > 1e-9 {
                        n_argmax_common_value_diff += 1;
                    }
                } else if java_only.contains(&jfnv) {
                    n_argmax_java_only += 1;
                } else if rust_only.contains(&r.argmax_fnv) {
                    n_argmax_rust_only += 1;
                } else if common.contains(&jfnv) && common.contains(&r.argmax_fnv) {
                    n_argmax_common_value_diff += 1;
                }
            }
        }
    }

    let row_pop_same = n_row_absent == 0 && java_row_keys.len() == rust_row_keys.len();
    let classification = if !row_pop_same {
        "LIKELIHOOD_ROW_POPULATION"
    } else if java_cols != rust_cols {
        "LIKELIHOOD_COLUMN_POPULATION"
    } else if JAVA.iter().any(|&(q, f, _, _, _, _, jll, _, _, _)| {
        rust_by
            .get(&(q.to_string(), f))
            .is_some_and(|r| (r.max_ll - jll).abs() > 1e-12)
    }) {
        "PRE_FILTER_LIKELIHOOD_VALUE"
    } else {
        "AFTER_MAX_LL"
    };

    eprintln!(
        "{}",
        json!({
            "classification": classification,
            "java_columns": java_cols.len(),
            "rust_columns": rust_cols.len(),
            "common_columns": common.len(),
            "java_only_columns": java_only.len(),
            "rust_only_columns": rust_only.len(),
            "java_filter_rows": 153,
            "rust_filter_rows": last_n,
            "n_cigar_eq": n_cigar_eq,
            "n_start_eq": n_start_eq,
            "n_same_argmax": n_same_argmax,
            "n_argmax_java_only": n_argmax_java_only,
            "n_argmax_rust_only": n_argmax_rust_only,
            "n_argmax_common_value_diff": n_argmax_common_value_diff,
            "vcf_ad": vcf.samples[0].ad,
        })
    );

    assert_eq!(JAVA_HAP_FNV.len(), 70);
    assert_eq!(JAVA.len(), 24);
    assert_eq!(n_row_absent, 0, "all 24 QNAME+flags exist as filter rows");
    assert_eq!(last_n, 153);
    assert_eq!(rust_cols.len(), 68);
    assert_eq!(java_only.len(), 2);
    assert_eq!(rust_only.len(), 0);
    assert_eq!(classification, "LIKELIHOOD_COLUMN_POPULATION");
    eprintln!("6R.94 classification={classification}");
}
