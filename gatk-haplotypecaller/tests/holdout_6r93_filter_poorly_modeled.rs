//! 6R.93 holdout: live poorly-modeled predicate inputs for 24 named QNAMEs.
//!
//! Skipped unless `HOLDOUT_6R93=1`. Coordinate-free contract lives in
//! `forensic_6r93_filter_poorly_modeled_predicate_contract`.
//!
//! Does not assert Rust AD == 36,19 or equal evidence membership.
//!
//! ```text
//! HOLDOUT_6R93=1 cargo test -p gatk-haplotypecaller --test holdout_6r93_filter_poorly_modeled -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::{
    begin_poorly_modeled_observe, call_disposition, flatten_assembly_regions,
    take_colocated_merge_numerics, take_poorly_modeled_observe, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use serde_json::json;
use std::collections::HashMap;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const POS_SNP: u64 = 29_456_344;

/// Live Java GATK 4.4.0.0 `filterPoorlyModeledEvidence` at the covering region
/// `20:29456294-29456500` (dump seq=6): qname, flags, qualifiedLen, threshold, max_ll, keep.
const JAVA: &[(&str, u16, usize, f64, f64, bool)] = &[
    (
        "HISEQ1:11:H8GV6ADXX:2:2216:2203:76921",
        147,
        130,
        -8.0,
        -7.450010299682617,
        true,
    ),
    (
        "HISEQ1:13:H8G92ADXX:1:1111:12251:89078",
        83,
        148,
        -8.0,
        -2.5224151611328125,
        true,
    ),
    (
        "HISEQ1:9:H8962ADXX:1:1112:19265:60083",
        99,
        76,
        -8.0,
        -2.4746932983398438,
        true,
    ),
    (
        "HWI-D00360:5:H814YADXX:1:1202:11051:34179",
        147,
        96,
        -8.0,
        -2.5021400451660156,
        true,
    ),
    (
        "HWI-D00360:5:H814YADXX:1:2207:10890:76583",
        147,
        148,
        -8.0,
        -2.6972694396972656,
        true,
    ),
    (
        "HWI-D00360:5:H814YADXX:2:1102:2154:52493",
        163,
        148,
        -8.0,
        -2.8041610717773438,
        true,
    ),
    (
        "HWI-D00360:6:H81VLADXX:2:1104:15554:2818",
        83,
        143,
        -8.0,
        -2.522064208984375,
        true,
    ),
    (
        "HWI-D00360:6:H81VLADXX:2:1202:18367:85709",
        163,
        134,
        -8.0,
        -2.6402854919433594,
        true,
    ),
    (
        "HWI-D00360:7:H88WKADXX:1:1116:9273:30844",
        83,
        132,
        -8.0,
        -2.5336990356445312,
        true,
    ),
    (
        "HWI-D00360:8:H88U0ADXX:1:2108:16806:75328",
        163,
        112,
        -8.0,
        -2.5473976135253906,
        true,
    ),
    (
        "HISEQ1:13:H8G92ADXX:1:1205:16330:83279",
        163,
        148,
        -8.0,
        -9.549158096313477,
        false,
    ),
    (
        "HISEQ1:9:H8962ADXX:2:1212:17767:73796",
        83,
        148,
        -8.0,
        -9.387346267700195,
        false,
    ),
    (
        "HWI-D00360:5:H814YADXX:2:2103:4936:45407",
        83,
        84,
        -8.0,
        -9.047658920288086,
        false,
    ),
    (
        "HWI-D00360:6:H81VLADXX:1:1103:1948:22968",
        147,
        148,
        -8.0,
        -12.48884391784668,
        false,
    ),
    (
        "HWI-D00360:6:H81VLADXX:1:1210:4156:72506",
        83,
        92,
        -8.0,
        -9.251577377319336,
        false,
    ),
    (
        "HWI-D00360:7:H88WKADXX:1:2111:4466:65743",
        147,
        144,
        -8.0,
        -8.373727798461914,
        false,
    ),
    (
        "HWI-D00360:7:H88WKADXX:1:2203:20480:101193",
        163,
        148,
        -8.0,
        -13.432897567749023,
        false,
    ),
    (
        "HWI-D00360:7:H88WKADXX:2:1214:6938:52704",
        83,
        125,
        -8.0,
        -8.664169311523438,
        false,
    ),
    (
        "HWI-D00360:8:H88U0ADXX:1:1205:11075:4786",
        147,
        148,
        -8.0,
        -11.753900527954102,
        false,
    ),
    (
        "HWI-D00360:8:H88U0ADXX:1:1213:18559:65935",
        163,
        148,
        -8.0,
        -11.865598678588867,
        false,
    ),
    (
        "HWI-D00360:8:H88U0ADXX:1:2213:15618:11579",
        163,
        105,
        -8.0,
        -9.34518814086914,
        false,
    ),
    (
        "HWI-D00360:8:H88U0ADXX:2:1213:15376:17578",
        163,
        104,
        -8.0,
        -18.327007293701172,
        false,
    ),
    (
        "HISEQ1:11:H8GV6ADXX:2:2105:12137:22761",
        163,
        110,
        -8.0,
        -9.151548385620117,
        false,
    ),
    (
        "HISEQ1:13:H8G92ADXX:1:1201:11859:45984",
        147,
        97,
        -8.0,
        -9.154520034790039,
        false,
    ),
];

const JAVA_LIVE_ONLY_N: usize = 10;
const RUST_ONLY_N: usize = 12;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

#[test]
fn holdout_6r93_filter_poorly_modeled_24_qnames() {
    if std::env::var("HOLDOUT_6R93").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R93=1");
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
    let snap = take_colocated_merge_numerics()
        .into_iter()
        .find(|n| n.loc == POS_SNP);

    let last_pass = observed.iter().map(|r| r.pass).max().unwrap_or(0);
    let last: Vec<_> = observed.iter().filter(|r| r.pass == last_pass).collect();
    let rust_by: HashMap<(String, u16), _> = last
        .iter()
        .map(|r| ((r.qname.clone(), r.flags), *r))
        .collect();

    eprintln!(
        "6R.93 observe_passes={} last_pass={} last_n={} covering_reads={}",
        observed.iter().map(|r| r.pass).max().unwrap_or(0),
        last_pass,
        last.len(),
        covering[0].reads.len()
    );
    eprintln!(
        "QNAME\tJava_qlen\tRust_qlen\tJava_thr\tRust_thr\tJava_max_ll\tRust_max_ll\tJava_dec\tRust_dec\textra_retain\tn_hap_cells"
    );

    let mut n_qlen_diff = 0usize;
    let mut n_thr_diff = 0usize;
    let mut n_ll_diff = 0usize;
    let mut n_ll_absent = 0usize;
    let mut n_extra = 0usize;
    let mut n_dec_diff = 0usize;

    for (i, &(q, flags, jqlen, jthr, jll, jkeep)) in JAVA.iter().enumerate() {
        let tag = if i < JAVA_LIVE_ONLY_N {
            "JAVA_LIVE_ONLY"
        } else if i < JAVA_LIVE_ONLY_N + RUST_ONLY_N {
            "RUST_ONLY"
        } else {
            "BOTH_DROP"
        };
        let rr = rust_by.get(&(q.to_string(), flags));
        let (rqlen, rthr, rll, rkeep, extra, cells) = match rr {
            Some(r) => (
                format!("{}", r.qual_len),
                format!("{}", r.threshold),
                format!("{}", r.max_ll),
                r.rust_keep,
                r.extra_retain,
                format!("{}", r.n_hap_cells),
            ),
            None => (
                "ABSENT".to_string(),
                "ABSENT".to_string(),
                "ABSENT".to_string(),
                false,
                false,
                "0".to_string(),
            ),
        };
        eprintln!(
            "{tag}\t{q}\t{jqlen}\t{rqlen}\t{jthr}\t{rthr}\t{jll}\t{rll}\t{jkeep}\t{rkeep}\t{extra}\t{cells}"
        );
        match rr {
            None => {
                n_ll_absent += 1;
                n_ll_diff += 1;
                if jkeep != false {
                    n_dec_diff += 1;
                }
            }
            Some(r) => {
                if r.qual_len != jqlen {
                    n_qlen_diff += 1;
                }
                if r.threshold != jthr {
                    n_thr_diff += 1;
                }
                if r.max_ll != jll {
                    n_ll_diff += 1;
                }
                if r.extra_retain {
                    n_extra += 1;
                }
                if r.rust_keep != jkeep {
                    n_dec_diff += 1;
                }
            }
        }
    }

    // Thresholds match for all 24 (−8.0): qlen integers may differ but cannot
    // flip keep/drop once every length is ≥ 51. Causal remaining input is max_ll.
    let classification = if n_thr_diff > 0 {
        "QUALIFIED_LENGTH"
    } else if n_ll_diff > 0 {
        "LIKELIHOOD_POPULATION"
    } else if n_dec_diff > 0 {
        "NOT_IN_FILTER_PREDICATE"
    } else {
        "NOT_IN_FILTER_PREDICATE"
    };

    eprintln!(
        "{}",
        json!({
            "classification": classification,
            "n_qlen_diff": n_qlen_diff,
            "n_thr_diff": n_thr_diff,
            "n_ll_diff": n_ll_diff,
            "n_ll_absent": n_ll_absent,
            "n_extra_retain": n_extra,
            "n_dec_diff": n_dec_diff,
            "last_pass": last_pass,
            "snap_n_reads": snap.as_ref().map(|s| s.n_reads),
            "vcf_ad": vcf.samples[0].ad,
        })
    );

    assert_eq!(JAVA.len(), 24);
    for &(q, flags, _, _, jll, jkeep) in JAVA.iter().take(JAVA_LIVE_ONLY_N) {
        assert!(jkeep, "{q}:{flags} Java KEEP");
        assert!(jll >= -8.0, "{q} Java max_ll above threshold");
    }
    for &(q, flags, _, _, jll, jkeep) in JAVA.iter().skip(JAVA_LIVE_ONLY_N).take(RUST_ONLY_N + 2) {
        assert!(!jkeep, "{q}:{flags} Java DROP");
        assert!(jll < -8.0, "{q} Java max_ll below threshold");
    }
    assert_eq!(n_thr_diff, 0, "all 24 thresholds must match at -8.0");
    assert_eq!(
        n_extra, 0,
        "P12 extra-retain bands must not fire on these 24"
    );
    assert_eq!(n_dec_diff, 0);
    assert_eq!(n_ll_diff, 24);
    assert_eq!(classification, "LIKELIHOOD_POPULATION");
    eprintln!("6R.93 classification={classification}");
}
