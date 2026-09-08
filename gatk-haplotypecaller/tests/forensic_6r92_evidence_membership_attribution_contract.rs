//! 6R.92 coordinate-free: attribute Java vs Rust annotation evidence membership.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//! `assignGenotypeLikelihoods`:
//!
//! ```text
//! mergedVC interval.expandWithinContig(margin=2)
//! retainEvidence(target.overlaps(read))   // GATKRead getStart/getEnd = post-realign alignment
//! prepareReadAlleleLikelihoodsForAnnotation
//!   default HC: reuse genotyping likelihoods
//!   addEvidence(overlapping filterNonPassingReads, 0)
//! DepthPerAlleleBySample.annotateWithLikelihoods
//! ```
//!
//! 6R.91: live Java annotation evidence is 60 unique QNAMEs; Rust overlap is 62.
//! Symmetric difference is 10 JAVA_LIVE_ONLY + 12 RUST_ONLY.
//!
//! First membership predicate: `filterPoorlyModeledEvidence`, not `retainEvidence` overlap.
//! Java `retainEvidence` also applies the overlap predicate to the poorly-modeled
//! filtered list, so overlapping disqualified reads remain in `filteredSampleEvidence`
//! (14 at this site) and are **not** annotation evidence. Likelihood cells are not compared.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r92_evidence_membership_attribution_contract
//! HOLDOUT_6R92=1 cargo test -p gatk-haplotypecaller --test forensic_6r92_evidence_membership_attribution_contract live_ -- --nocapture
//! ```

use gatk_haplotypecaller::hc_genotyping_engine::java_alignment_read_overlaps_interval;
use gatk_haplotypecaller::read_unclip::alignment_end_1based;
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::Record;
use std::collections::{HashMap, HashSet};

/// Live Java seq=52 `AlleleLikelihoods.sampleEvidence`: qname, flags, start, end, cigar.
const JAVA_LIVE: &[(&str, u16, i64, i64, &str)] = &[
    (
        "HISEQ1:11:H8GV6ADXX:1:1213:10338:76652",
        163,
        29456268,
        29456411,
        "4H144M",
    ),
    (
        "HISEQ1:12:H8GVUADXX:2:1216:17904:51162",
        163,
        29456268,
        29456350,
        "65H83M",
    ),
    (
        "HISEQ1:13:H8G92ADXX:2:1205:1862:61782",
        99,
        29456268,
        29456361,
        "54H94M",
    ),
    (
        "HISEQ1:9:H8962ADXX:1:1112:19265:60083",
        99,
        29456268,
        29456343,
        "72H76M",
    ),
    (
        "HISEQ1:9:H8962ADXX:1:2109:11968:22409",
        99,
        29456268,
        29456397,
        "18H130M",
    ),
    (
        "HWI-D00360:5:H814YADXX:2:2211:1401:25869",
        99,
        29456268,
        29456351,
        "64H84M",
    ),
    (
        "HWI-D00360:6:H81VLADXX:1:1102:1483:21701",
        163,
        29456268,
        29456365,
        "50H98M",
    ),
    (
        "HWI-D00360:6:H81VLADXX:2:1116:2361:8950",
        163,
        29456268,
        29456355,
        "60H88M",
    ),
    (
        "HWI-D00360:6:H81VLADXX:2:1202:18367:85709",
        163,
        29456268,
        29456401,
        "14H134M",
    ),
    (
        "HWI-D00360:6:H81VLADXX:2:2103:8562:24900",
        99,
        29456268,
        29456393,
        "23H76M1D49M",
    ),
    (
        "HWI-D00360:7:H88WKADXX:2:2106:13368:29882",
        163,
        29456268,
        29456363,
        "52H96M",
    ),
    (
        "HWI-D00360:8:H88U0ADXX:1:1114:5563:73517",
        99,
        29456268,
        29456358,
        "57H91M",
    ),
    (
        "HWI-D00360:8:H88U0ADXX:1:2108:16806:75328",
        163,
        29456268,
        29456379,
        "36H112M",
    ),
    (
        "HWI-D00360:8:H88U0ADXX:2:1206:5788:95711",
        163,
        29456268,
        29456342,
        "73H75M",
    ),
    (
        "HISEQ1:11:H8GV6ADXX:1:1108:6005:20931",
        83,
        29456268,
        29456388,
        "27H121M",
    ),
    (
        "HISEQ1:11:H8GV6ADXX:2:2216:2203:76921",
        147,
        29456268,
        29456398,
        "18H76M1D54M",
    ),
    (
        "HISEQ1:12:H8GVUADXX:1:1106:3080:10032",
        147,
        29456268,
        29456349,
        "66H82M",
    ),
    (
        "HISEQ1:13:H8G92ADXX:1:1103:8983:35134",
        83,
        29456268,
        29456348,
        "67H81M",
    ),
    (
        "HISEQ1:13:H8G92ADXX:1:1214:8757:76066",
        83,
        29456268,
        29456369,
        "47H77M1D24M",
    ),
    (
        "HISEQ1:9:H8962ADXX:2:2111:4053:81460",
        83,
        29456268,
        29456376,
        "39H109M",
    ),
    (
        "HWI-D00360:5:H814YADXX:1:1202:11051:34179",
        147,
        29456268,
        29456363,
        "52H96M",
    ),
    (
        "HWI-D00360:5:H814YADXX:1:1208:17759:33171",
        83,
        29456268,
        29456375,
        "40H108M",
    ),
    (
        "HWI-D00360:5:H814YADXX:2:1109:4475:85583",
        147,
        29456268,
        29456398,
        "18H77M1D53M",
    ),
    (
        "HWI-D00360:5:H814YADXX:2:2107:1177:54891",
        83,
        29456268,
        29456360,
        "55H93M",
    ),
    (
        "HWI-D00360:5:H814YADXX:2:2209:16484:91770",
        147,
        29456268,
        29456397,
        "18H130M",
    ),
    (
        "HWI-D00360:6:H81VLADXX:2:1104:15554:2818",
        83,
        29456268,
        29456410,
        "5H143M",
    ),
    (
        "HWI-D00360:7:H88WKADXX:1:1116:9273:30844",
        83,
        29456268,
        29456399,
        "16H132M",
    ),
    (
        "HWI-D00360:7:H88WKADXX:2:1210:8920:81214",
        147,
        29456268,
        29456364,
        "51H97M",
    ),
    (
        "HWI-D00360:7:H88WKADXX:2:2105:16900:85825",
        83,
        29456268,
        29456345,
        "71H76M1D1M",
    ),
    (
        "HWI-D00360:8:H88U0ADXX:1:1101:14447:3965",
        147,
        29456268,
        29456405,
        "11H77M1D60M",
    ),
    (
        "HWI-D00360:8:H88U0ADXX:1:1106:15654:72455",
        83,
        29456268,
        29456389,
        "26H122M",
    ),
    (
        "HWI-D00360:8:H88U0ADXX:1:1202:9544:72702",
        147,
        29456268,
        29456402,
        "13H135M",
    ),
    (
        "HWI-D00360:8:H88U0ADXX:1:2201:20293:70383",
        83,
        29456268,
        29456362,
        "54H77M1D17M",
    ),
    (
        "HISEQ1:9:H8962ADXX:2:2101:14564:85219",
        163,
        29456272,
        29456419,
        "148M",
    ),
    (
        "HISEQ1:9:H8962ADXX:2:1206:11075:50547",
        147,
        29456273,
        29456420,
        "148M",
    ),
    (
        "HWI-D00360:6:H81VLADXX:1:2212:1317:49311",
        163,
        29456275,
        29456422,
        "148M",
    ),
    (
        "HISEQ1:9:H8962ADXX:1:1213:5437:90850",
        147,
        29456276,
        29456423,
        "148M",
    ),
    (
        "HWI-D00360:6:H81VLADXX:2:2205:5601:38112",
        147,
        29456279,
        29456426,
        "148M",
    ),
    (
        "HISEQ1:11:H8GV6ADXX:1:1104:4186:41128",
        83,
        29456283,
        29456430,
        "148M",
    ),
    (
        "HISEQ1:9:H8962ADXX:2:1104:19714:91101",
        147,
        29456286,
        29456433,
        "148M",
    ),
    (
        "HISEQ1:9:H8962ADXX:2:1105:19831:29366",
        83,
        29456286,
        29456433,
        "148M",
    ),
    (
        "HWI-D00360:7:H88WKADXX:1:1112:12463:50733",
        147,
        29456286,
        29456433,
        "148M",
    ),
    (
        "HISEQ1:13:H8G92ADXX:1:2103:3021:81022",
        163,
        29456289,
        29456436,
        "148M",
    ),
    (
        "HWI-D00360:8:H88U0ADXX:1:1106:10867:71386",
        83,
        29456290,
        29456437,
        "148M",
    ),
    (
        "HWI-D00360:5:H814YADXX:1:1212:12402:73354",
        83,
        29456294,
        29456441,
        "148M",
    ),
    (
        "HISEQ1:13:H8G92ADXX:2:1209:13449:57024",
        163,
        29456296,
        29456443,
        "148M",
    ),
    (
        "HISEQ1:13:H8G92ADXX:1:1111:12251:89078",
        83,
        29456301,
        29456448,
        "148M",
    ),
    (
        "HWI-D00360:8:H88U0ADXX:2:1204:18013:40432",
        83,
        29456306,
        29456453,
        "148M",
    ),
    (
        "HISEQ1:13:H8G92ADXX:2:1101:15512:24344",
        99,
        29456316,
        29456463,
        "148M",
    ),
    (
        "HWI-D00360:7:H88WKADXX:1:2101:19221:59478",
        147,
        29456316,
        29456463,
        "148M",
    ),
    (
        "HWI-D00360:6:H81VLADXX:2:2204:18162:74438",
        99,
        29456319,
        29456466,
        "148M",
    ),
    (
        "HWI-D00360:5:H814YADXX:1:2207:10890:76583",
        147,
        29456319,
        29456466,
        "148M",
    ),
    (
        "HWI-D00360:8:H88U0ADXX:2:1116:13828:55178",
        147,
        29456319,
        29456466,
        "148M",
    ),
    (
        "HWI-D00360:6:H81VLADXX:1:2101:16091:54483",
        83,
        29456323,
        29456470,
        "148M",
    ),
    (
        "HISEQ1:11:H8GV6ADXX:2:2205:10648:31409",
        147,
        29456324,
        29456471,
        "148M",
    ),
    (
        "HISEQ1:9:H8962ADXX:2:2115:11042:98655",
        99,
        29456333,
        29456480,
        "148M",
    ),
    (
        "HWI-D00360:5:H814YADXX:2:1102:2154:52493",
        163,
        29456335,
        29456482,
        "148M",
    ),
    (
        "HWI-D00360:7:H88WKADXX:2:1106:19349:74372",
        163,
        29456345,
        29456492,
        "148M",
    ),
    (
        "HISEQ1:13:H8G92ADXX:1:2209:4644:71639",
        147,
        29456345,
        29456492,
        "148M",
    ),
    (
        "HISEQ1:9:H8962ADXX:2:1209:5235:78510",
        99,
        29456347,
        29456494,
        "148M",
    ),
];

/// Live Java seq=52 `filteredSampleEvidence` after `retainEvidence` (overlapping poorly-modeled).
const JAVA_FILTERED_OVERLAP: &[&str] = &[
    "HISEQ1:11:H8GV6ADXX:2:2105:12137:22761",
    "HISEQ1:13:H8G92ADXX:1:1201:11859:45984",
    "HISEQ1:13:H8G92ADXX:1:1205:16330:83279",
    "HISEQ1:9:H8962ADXX:2:1212:17767:73796",
    "HWI-D00360:5:H814YADXX:2:2103:4936:45407",
    "HWI-D00360:6:H81VLADXX:1:1103:1948:22968",
    "HWI-D00360:6:H81VLADXX:1:1210:4156:72506",
    "HWI-D00360:7:H88WKADXX:1:2111:4466:65743",
    "HWI-D00360:7:H88WKADXX:1:2203:20480:101193",
    "HWI-D00360:7:H88WKADXX:2:1214:6938:52704",
    "HWI-D00360:8:H88U0ADXX:1:1205:11075:4786",
    "HWI-D00360:8:H88U0ADXX:1:1213:18559:65935",
    "HWI-D00360:8:H88U0ADXX:1:2213:15618:11579",
    "HWI-D00360:8:H88U0ADXX:2:1213:15376:17578",
];

fn java_simple_interval_overlaps(tstart: i64, tend: i64, rstart: i64, rend: i64) -> bool {
    tstart <= rend && rstart <= tend
}

fn format_cigar(rec: &Record) -> String {
    rec.cigar()
        .iter()
        .map(|c| {
            let (n, op) = match c {
                Cigar::Match(n) => (*n, 'M'),
                Cigar::Ins(n) => (*n, 'I'),
                Cigar::Del(n) => (*n, 'D'),
                Cigar::SoftClip(n) => (*n, 'S'),
                Cigar::HardClip(n) => (*n, 'H'),
                Cigar::Equal(n) => (*n, '='),
                Cigar::Diff(n) => (*n, 'X'),
                Cigar::RefSkip(n) => (*n, 'N'),
                Cigar::Pad(n) => (*n, 'P'),
            };
            format!("{n}{op}")
        })
        .collect()
}

fn rec_align(rec: &Record) -> (i64, i64) {
    (
        (rec.pos() + 1).max(1),
        i64::from(alignment_end_1based(rec).max(1)),
    )
}

fn flag_bits(flags: u16) -> String {
    let mut bits = Vec::new();
    if flags & 0x1 != 0 {
        bits.push("PAIRED");
    }
    if flags & 0x40 != 0 {
        bits.push("R1");
    }
    if flags & 0x80 != 0 {
        bits.push("R2");
    }
    if flags & 0x10 != 0 {
        bits.push("REV");
    }
    if flags & 0x100 != 0 {
        bits.push("SECONDARY");
    }
    if flags & 0x200 != 0 {
        bits.push("QCFAIL");
    }
    if flags & 0x400 != 0 {
        bits.push("DUP");
    }
    if flags & 0x800 != 0 {
        bits.push("SUPP");
    }
    bits.join(",")
}

/// Java `SimpleInterval(mergedVC).expandWithinContig(2)` for untrimmed TG at 29456344-29456345.
#[test]
fn forensic_6r92_merged_tg_interval_is_two_bp_plus_margin() {
    let merged_start = 29_456_344i64;
    let merged_end = 29_456_345i64;
    let (tstart, tend) = (merged_start - 2, merged_end + 2);
    assert_eq!((tstart, tend), (29_456_342, 29_456_347));
    let snp_end = merged_start;
    let snp_tend = snp_end + 2;
    assert_ne!(
        tend, snp_tend,
        "1bp SNP ±2 is not the live Java variantCallingSubset"
    );
}

/// A read whose alignment starts on the extra right-edge base is kept only by the 2bp interval.
#[test]
fn forensic_6r92_right_edge_start_depends_on_merged_end() {
    let read_start = 29_456_347i64;
    let read_end = 29_456_494i64;
    assert!(java_simple_interval_overlaps(
        29_456_342, 29_456_347, read_start, read_end
    ));
    assert!(
        !java_simple_interval_overlaps(29_456_342, 29_456_346, read_start, read_end),
        "1bp-ref ±2 drops a read that starts at 29456347"
    );
}

/// A left-hard-clipped realigned read ending on the expanded start is kept; one base shorter is not.
#[test]
fn forensic_6r92_left_edge_uses_post_realign_alignment_end() {
    assert!(java_simple_interval_overlaps(
        29_456_342, 29_456_347, 29_456_268, 29_456_342
    ));
    assert!(!java_simple_interval_overlaps(
        29_456_342, 29_456_347, 29_456_268, 29_456_341
    ));
}

/// Same QNAME with different FLAGS is a different SAM record (mates).
#[test]
fn forensic_6r92_qname_alone_does_not_identify_a_record() {
    let a = ("q", 99u16, 100i64);
    let b = ("q", 147u16, 200i64);
    assert_eq!(a.0, b.0);
    assert_ne!((a.1, a.2), (b.1, b.2));
}

/// `retainEvidence` keeps overlapping poorly-modeled reads in `filteredSampleEvidence` only.
#[test]
fn forensic_6r92_overlapping_poorly_modeled_is_filtered_not_annotation() {
    let annotation = 60usize;
    let overlapping_poorly_modeled = 14usize;
    let annotation_plus_filtered = annotation + overlapping_poorly_modeled;
    assert_ne!(
        annotation_plus_filtered, annotation,
        "overlapping poorly-modeled reads are not sampleEvidence"
    );
    assert_eq!(JAVA_FILTERED_OVERLAP.len(), overlapping_poorly_modeled);
    let annot: HashSet<&str> = JAVA_LIVE.iter().map(|r| r.0).collect();
    for q in JAVA_FILTERED_OVERLAP {
        assert!(
            !annot.contains(q),
            "filtered overlapping poorly-modeled must not appear in annotation evidence"
        );
    }
}

/// Overlap on genotyping_reads is not annotation membership; poorly-modeled must keep the read.
#[test]
fn forensic_6r92_annotation_requires_poorly_modeled_keep_not_only_overlap() {
    let rust_overlap = 62usize;
    let java_annotation = 60usize;
    assert_ne!(
        rust_overlap, java_annotation,
        "Rust overlap ∩ PairHMM survivors is not Java annotation evidence"
    );
    assert_eq!(50 + 10, java_annotation);
    assert_eq!(50 + 12, rust_overlap);
}

/// Live Java evidence is 60 unique QNAMEs; ordinary primary paired flags only.
#[test]
fn forensic_6r92_java_live_is_sixty_primary_records() {
    assert_eq!(JAVA_LIVE.len(), 60);
    let q: HashSet<&str> = JAVA_LIVE.iter().map(|r| r.0).collect();
    assert_eq!(q.len(), 60);
    for &(_, flags, _, _, _) in JAVA_LIVE {
        assert_eq!(flags & 0x100, 0, "no secondary");
        assert_eq!(flags & 0x800, 0, "no supplementary");
        assert_eq!(flags & 0x400, 0, "no duplicate");
        assert_eq!(flags & 0x200, 0, "no QC-fail");
        assert_ne!(flags & 0x1, 0, "paired");
    }
}

#[test]
fn live_java_vs_rust_membership_attribution() {
    if std::env::var("HOLDOUT_6R92").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R92=1");
        return;
    }
    use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
    use gatk_haplotypecaller::read_assembly_filter::{
        passes_assembly_read, AssemblyReadFilterConfig,
    };
    use gatk_haplotypecaller::{
        call_disposition, flatten_assembly_regions, take_colocated_merge_numerics,
        traverse_assembly_region_walker, AssemblyRegionCallDisposition, CallRegionArgs,
        HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
        DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN,
    };
    use std::path::Path;

    const INTERVAL: &str = "20:29455000-29456500";
    const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
    const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
    const POS_SNP: u64 = 29_456_344;
    const JAVA_TSTART: i64 = 29_456_342;
    const JAVA_TEND: i64 = 29_456_347;

    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: live BAM/ref missing");
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
    let orig_by_qname: HashMap<String, Vec<(u16, i64, i64, String, u8, bool)>> = {
        let mut m: HashMap<String, Vec<_>> = HashMap::new();
        for rec in &covering[0].reads {
            let q = String::from_utf8_lossy(rec.qname()).into_owned();
            let (s, e) = rec_align(rec);
            m.entry(q).or_default().push((
                rec.flags(),
                s,
                e,
                format_cigar(rec),
                rec.mapq(),
                rec.flags() & 0xF00 != 0,
            ));
        }
        m
    };
    let outcome = HaplotypeCallerEngine::call_region(
        covering[0],
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let snap = take_colocated_merge_numerics()
        .into_iter()
        .find(|n| n.loc == POS_SNP)
        .expect("snap");

    let loc = POS_SNP;
    let end_long = loc.saturating_add(snap.long_ref.len().saturating_sub(1) as u64);
    let end_snp = loc;
    let margin = DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN;
    let pairhmm: HashSet<usize> = outcome
        .read_likelihoods
        .iter()
        .map(|e| e.read_index.get())
        .collect();

    #[derive(Clone)]
    struct RustRec {
        idx: usize,
        qname: String,
        flags: u16,
        start: i64,
        end: i64,
        cigar: String,
        mapq: u8,
        in_pairhmm: bool,
        overlap_long: bool,
        overlap_snp: bool,
        overlap_java_interval: bool,
        passing: bool,
    }
    let filter_cfg = AssemblyReadFilterConfig::gatk_defaults();
    let rust_all: Vec<RustRec> = outcome
        .genotyping_reads
        .iter()
        .enumerate()
        .map(|(idx, rec)| {
            let (start, end) = rec_align(rec);
            let in_pairhmm = pairhmm.contains(&idx);
            RustRec {
                idx,
                qname: String::from_utf8_lossy(rec.qname()).into_owned(),
                flags: rec.flags(),
                start,
                end,
                cigar: format_cigar(rec),
                mapq: rec.mapq(),
                in_pairhmm,
                overlap_long: java_alignment_read_overlaps_interval(rec, loc, end_long, margin),
                overlap_snp: java_alignment_read_overlaps_interval(rec, loc, end_snp, margin),
                overlap_java_interval: java_simple_interval_overlaps(
                    JAVA_TSTART,
                    JAVA_TEND,
                    start,
                    end,
                ),
                passing: passes_assembly_read(rec, &filter_cfg),
            }
        })
        .collect();
    let rust_overlap: Vec<&RustRec> = rust_all
        .iter()
        .filter(|r| r.in_pairhmm && r.overlap_long)
        .collect();
    let rust_q: HashSet<&str> = rust_overlap.iter().map(|r| r.qname.as_str()).collect();
    let java_q: HashSet<&str> = JAVA_LIVE.iter().map(|r| r.0).collect();
    let mut common: Vec<&str> = java_q.intersection(&rust_q).copied().collect();
    let mut java_only: Vec<&str> = java_q.difference(&rust_q).copied().collect();
    let mut rust_only: Vec<&str> = rust_q.difference(&java_q).copied().collect();
    common.sort_unstable();
    java_only.sort_unstable();
    rust_only.sort_unstable();

    eprintln!(
        "6R.92 source covering={} genotyping={} pairhmm={} rust_overlap_long={} rust_overlap_snp={} java_live={} long_ref={} loc={} end_long={} snap_n={}",
        covering[0].reads.len(),
        outcome.genotyping_reads.len(),
        pairhmm.len(),
        rust_overlap.len(),
        rust_all.iter().filter(|r| r.in_pairhmm && r.overlap_snp).count(),
        JAVA_LIVE.len(),
        snap.long_ref,
        loc,
        end_long,
        snap.n_reads
    );
    eprintln!(
        "6R.92 sets common={} JAVA_LIVE_ONLY={} RUST_ONLY={} invariant_java={} invariant_rust={}",
        common.len(),
        java_only.len(),
        rust_only.len(),
        common.len() + java_only.len(),
        common.len() + rust_only.len()
    );

    let rust_by_q: HashMap<&str, Vec<&RustRec>> = {
        let mut m: HashMap<&str, Vec<&RustRec>> = HashMap::new();
        for r in &rust_all {
            m.entry(r.qname.as_str()).or_default().push(r);
        }
        m
    };
    let rust_ov_by_q: HashMap<&str, &RustRec> = rust_overlap
        .iter()
        .map(|r| (r.qname.as_str(), *r))
        .collect();
    let java_by_q: HashMap<&str, &(&str, u16, i64, i64, &str)> =
        JAVA_LIVE.iter().map(|r| (r.0, r)).collect();

    let dump_one = |tag: &str, q: &str| {
        let j = java_by_q.get(q);
        let ro = rust_ov_by_q.get(q);
        let r_all = rust_by_q.get(q).cloned().unwrap_or_default();
        let orig = orig_by_qname.get(q);
        let (jf, js, je, jc) = j.map(|t| (t.1, t.2, t.3, t.4)).unwrap_or((0, 0, 0, "."));
        let j_ov = j
            .map(|t| java_simple_interval_overlaps(JAVA_TSTART, JAVA_TEND, t.2, t.3))
            .unwrap_or(false);
        eprintln!(
            "6R.92 {tag}\tq={q}\tmate={}\tjava_flags={jf}\tjava={js}-{je}\tjava_cigar={jc}\tjava_interval_ov={j_ov}",
            flag_bits(jf)
        );
        if let Some(rs) = orig {
            for (i, (f, s, e, c, mq, special)) in rs.iter().enumerate() {
                let ov = java_simple_interval_overlaps(JAVA_TSTART, JAVA_TEND, *s, *e);
                eprintln!(
                    "6R.92 {tag}\t  covering[{i}] flags={f} {} {s}-{e} cigar={c} mq={mq} special={special} java_interval_ov={ov}",
                    flag_bits(*f)
                );
            }
        } else {
            eprintln!("6R.92 {tag}\t  covering=MISSING");
        }
        if r_all.is_empty() {
            eprintln!("6R.92 {tag}\t  genotyping=MISSING");
        }
        for r in &r_all {
            eprintln!(
                "6R.92 {tag}\t  gt idx={} flags={} {} {}-{} cigar={} mq={} pairhmm={} ov_long={} ov_snp={} ov_java_iv={} passing={} in_rust_overlap={}",
                r.idx,
                r.flags,
                flag_bits(r.flags),
                r.start,
                r.end,
                r.cigar,
                r.mapq,
                r.in_pairhmm,
                r.overlap_long,
                r.overlap_snp,
                r.overlap_java_interval,
                r.passing,
                ro.is_some()
            );
        }
        if let (Some(j), Some(r)) = (j, ro) {
            let coord_eq = j.2 == r.start && j.3 == r.end && j.4 == r.cigar && j.1 == r.flags;
            eprintln!("6R.92 {tag}\t  live_vs_rust_coord_eq={coord_eq}");
        }
    };

    eprintln!("6R.92 --- JAVA_LIVE_ONLY ---");
    for q in &java_only {
        dump_one("JAVA_LIVE_ONLY", q);
    }
    eprintln!("6R.92 --- RUST_ONLY ---");
    for q in &rust_only {
        dump_one("RUST_ONLY", q);
    }

    let mut n_coord_mismatch_common = 0usize;
    for q in &common {
        if let (Some(j), Some(r)) = (java_by_q.get(q), rust_ov_by_q.get(q)) {
            if j.2 != r.start || j.3 != r.end || j.4 != r.cigar {
                n_coord_mismatch_common += 1;
            }
        }
    }
    eprintln!(
        "6R.92 common_coord_mismatch={n_coord_mismatch_common}/{}",
        common.len()
    );

    let java_filtered: HashSet<&str> = JAVA_FILTERED_OVERLAP.iter().copied().collect();
    for q in &java_only {
        let recs = rust_by_q.get(q).cloned().unwrap_or_default();
        let r = recs
            .iter()
            .find(|x| x.flags == java_by_q.get(q).map(|j| j.1).unwrap_or(0))
            .or(recs.first())
            .expect("JAVA_LIVE_ONLY in genotyping_reads");
        assert!(
            r.overlap_long && r.passing,
            "{q}: Rust overlap would keep JAVA_LIVE_ONLY"
        );
        assert!(
            !r.in_pairhmm,
            "{q}: Rust filterPoorlyModeledEvidence dropped JAVA_LIVE_ONLY"
        );
    }
    for q in &rust_only {
        let r = rust_ov_by_q.get(q).expect("RUST_ONLY in rust overlap");
        assert!(r.in_pairhmm && r.overlap_long && r.passing);
        assert!(
            java_filtered.contains(q),
            "{q}: RUST_ONLY must be Java overlapping poorly-modeled filtered"
        );
    }
    let rust_only_set: HashSet<&str> = rust_only.iter().copied().collect();
    let extra_java_filtered: Vec<&str> = JAVA_FILTERED_OVERLAP
        .iter()
        .copied()
        .filter(|q| !rust_only_set.contains(q))
        .collect();
    eprintln!(
        "6R.92 cause=filterPoorlyModeledEvidence JAVA_LIVE_ONLY_rust_drop={} RUST_ONLY_java_poorly_modeled={} extra_java_filtered_both_drop={}",
        java_only.len(),
        rust_only.len(),
        extra_java_filtered.len()
    );
    for q in &extra_java_filtered {
        eprintln!("6R.92 extra_java_filtered {q}");
    }
    assert_eq!(extra_java_filtered.len(), 2);

    assert_eq!(JAVA_LIVE.len(), 60);
    assert_eq!(common.len() + java_only.len(), 60);
    assert_eq!(common.len() + rust_only.len(), rust_overlap.len());
    assert_eq!(snap.n_reads, rust_overlap.len());
    assert_eq!(common.len(), 50);
    assert_eq!(java_only.len(), 10);
    assert_eq!(rust_only.len(), 12);
}
