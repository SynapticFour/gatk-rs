use super::*;

#[test]
fn p11_pl_yields_java_af_qual() {
    let pl = [2238_i32, 151, 0];
    let gl: Vec<f64> = pl.iter().map(|p| -((p - pl[2]) as f64) / 10.0).collect();
    let qual = qual_from_af_calculation(&gl).expect("qual");
    assert!((qual - 2224.06).abs() < 0.15, "qual={qual}");
}

#[test]
fn six_r64_java_emitted_pl_does_not_reproduce_java_qual() {
    // Java VCF PL=542,0,1353 QUAL=510.06. Rust AFCalculator on those PLs is ~534.64,
    // so Java QUAL is not AF(emitted biallelic PL). GenotypingEngine computes QUAL
    // from AF on the pre-subset merged VC (6 GLs), then copies it through subset+trim.
    let gl = [-54.2, 0.0, -135.3];
    let qual = qual_from_af_calculation(&gl).expect("qual");
    assert!(
        (qual - 534.64).abs() < 0.15,
        "Rust AF on Java emitted PL: got {qual}"
    );
    assert!(
        (qual - 510.06).abs() > 10.0,
        "must not confuse AF(emitted PL) with Java site QUAL 510.06, got {qual}"
    );
}

#[test]
fn six_r64_rust_pl_298_0_1103_qual_matches_emitted() {
    // Rust 6R.63 VCF: PL=298,0,1103 QUAL=290.64.
    let gl = [-29.8, 0.0, -110.3];
    let qual = qual_from_af_calculation(&gl).expect("qual");
    assert!(
        (qual - 290.64).abs() < 0.15,
        "Rust PL→QUAL: got {qual}, expected ~290.64"
    );
}

#[test]
fn six_r64_rust_merged_6gl_qual_is_not_java_510() {
    // AF on merged 6-state GLs (Java QUAL source stage). If this were ~510, QUAL
    // would be a post-subset AF-order issue. It is not — 6-state inputs already differ.
    let gl6 = [-29.8, -33.7, -162.0, 0.0, -105.8, -110.3];
    let af =
        crate::af_calc::calculate_multiallelic_af_em(&[&gl6], 3, &AfCalculatorConfig::default())
            .expect("af6");
    let qual = (-10.0 * af.log10_posterior_no_variant) + 0.0;
    // Measured: same as AF on the subsetted 3-GLs (290.64). Java lifecycle F
    // (AF before unused-ALT subset) would not produce 510.06 from these 6-GLs.
    assert!(
        (qual - 290.64).abs() < 0.15,
        "Rust 6-GL AF QUAL must match subsetted QUAL ~290.64, not Java 510.06; got {qual}"
    );
}

#[test]
fn six_r40_pl_90_6_0_qual_mleac_match_java_af() {
    let gl = [-9.0, -0.6, 0.0];
    let af = calculate_biallelic_af_em(&[&gl], &AfCalculatorConfig::default()).expect("af");
    let qual = (-10.0 * af.log10_posterior_no_variant) + 0.0;
    assert_eq!(af.alt_allele_count, 1);
    assert!((qual - 78.32).abs() < 0.02, "qual={qual}");
    let mleaf = (af.alt_allele_count as f64 / 2.0).min(1.0);
    assert!((mleaf - 0.5).abs() < 1e-9);
    let _qd = crate::annotator::plugins::qual_by_depth::hold_process_qd_rng_for_test();
    crate::annotator::plugins::qual_by_depth::reset_gatk_qual_by_depth_rng();
    let qd = crate::annotator::plugins::qual_by_depth::qual_by_depth(qual, 2);
    assert!(
        (qd - 25.36).abs() < 0.005,
        "QD first GATK-seed gaussian, qd={qd}"
    );
}

fn bam_rec(
    qname: &[u8],
    pos0: i64,
    mapq: u8,
    paired: bool,
    unmapped: bool,
    mate_unmapped: bool,
    mtid: i32,
) -> SharedBamRecord {
    use rust_htslib::bam::record::{Cigar, CigarString};
    use rust_htslib::bam::HeaderView;
    use std::sync::Arc;
    let mut r = Record::new();
    r.set_header(Arc::new(HeaderView::from_bytes(
        b"@HD\tVN:1.0\n@SQ\tSN:chr1\tLN:10000\n@SQ\tSN:chr2\tLN:10000\n",
    )));
    r.set(
        qname,
        Some(&CigarString::from(vec![Cigar::Match(10)])),
        b"AAAAAAAAAA",
        &vec![30u8; 10],
    );
    r.set_tid(0);
    r.set_pos(pos0);
    r.set_mapq(mapq);
    r.set_mtid(mtid);
    if paired {
        r.set_paired();
    }
    if unmapped {
        r.set_unmapped();
    } else {
        r.unset_unmapped();
    }
    if mate_unmapped {
        r.set_mate_unmapped();
    } else {
        r.unset_mate_unmapped();
    }
    crate::shared_bam::share_record(r)
}

fn ll_cell(read_index: usize) -> RegionReadLikelihood {
    RegionReadLikelihood {
        read_index: crate::bio_ids::ReadIndex::new(read_index),
        haplotype_index: crate::bio_ids::HaplotypeIndex::new(0),
        log10_likelihood: -1.0,
    }
}

#[test]
fn mq_adds_overlapping_coverage_only_failed_mate_without_a_likelihood_row() {
    use crate::engine::engine_observe::{
        clear_coverage_only_failed_mates, forensic_6r306_note_pre_filter_reads,
        production_failed_mate_key, replace_coverage_only_failed_mates,
        with_coverage_only_failed_mates,
    };
    clear_coverage_only_failed_mates();
    let same = bam_rec(b"same-contig", 99, 53, true, false, false, 0);
    let ordinary = bam_rec(b"ordinary", 99, 43, false, false, false, 0);
    let mate_unmapped = bam_rec(b"mate-unmapped", 99, 30, true, false, true, 1);
    let other = bam_rec(b"other-contig", 99, 25, true, false, false, 1);
    let other_far = bam_rec(b"other-far", 5000, 10, true, false, false, 1);
    let unavailable = bam_rec(b"other-unavail", 99, 255, true, false, false, 1);
    let originals = [
        same.clone(),
        ordinary.clone(),
        mate_unmapped.clone(),
        other.clone(),
        other_far.clone(),
        unavailable.clone(),
    ];
    forensic_6r306_note_pre_filter_reads(&originals);
    let mut scored = Vec::new();
    let mut coverage_only = Vec::new();
    for rec in &originals {
        if production_failed_mate_key(rec.qname(), rec.flags()) {
            coverage_only.push(rec.clone());
        } else {
            scored.push(rec.clone());
        }
    }
    assert!(scored.iter().any(|r| r.qname() == b"same-contig"));
    assert!(scored.iter().any(|r| r.qname() == b"ordinary"));
    assert!(scored.iter().any(|r| r.qname() == b"mate-unmapped"));
    assert!(coverage_only.iter().all(|r| {
        r.qname() == b"other-contig" || r.qname() == b"other-far" || r.qname() == b"other-unavail"
    }));
    assert!(!coverage_only
        .iter()
        .any(|r| scored.iter().any(|s| s.qname() == r.qname())));
    replace_coverage_only_failed_mates(coverage_only);
    let cells = vec![ll_cell(0), ll_cell(1), ll_cell(2)];
    let with_coverage = mq_rms_of_sample_evidence(&scored, &cells, 100, 100, 0);
    let likelihood_only = {
        let mqs = rms_mapping_quality_sample_mapqs(&scored, &cells);
        java_finalize_mq(rms_mapping_quality_raw(&mqs).expect("raw").2)
    };
    let expected = java_finalize_mq(rms_mapping_quality_raw(&[53, 43, 30, 25]).expect("raw").2);
    assert!(
        (with_coverage - expected).abs() < 1e-9,
        "overlapping failed mate enters MQ once, got {with_coverage} expected {expected}"
    );
    assert!(
        (with_coverage - likelihood_only).abs() > 1.0,
        "coverage-only MAPQ must change the RMS"
    );
    assert!(
        (likelihood_only - java_finalize_mq(rms_mapping_quality_raw(&[53, 43, 30]).unwrap().2))
            .abs()
            < 1e-9
    );
    assert_eq!(cells.len(), 3, "failed mate has no likelihood row");
    assert!(with_coverage_only_failed_mates(|mates| {
        mates.iter().any(|r| r.qname() == b"other-contig")
    }));
    clear_coverage_only_failed_mates();
    let cleared = mq_rms_of_sample_evidence(&scored, &cells, 100, 100, 0);
    assert!(
        (cleared - likelihood_only).abs() < 1e-9,
        "without the coverage-only list, MQ stays on genotyping MAPQs"
    );
}

fn allele_confidence(best_l: f64, second_l: f64) -> f64 {
    if best_l == second_l {
        0.0
    } else {
        best_l - second_l
    }
}

#[test]
fn search_best_allele_reference_priority_breaks_near_tie() {
    // The second allele is 0.1 above REF. That gap is below 0.2, so REF priority 1.0 wins.
    let (best, second, best_l, second_l) = search_best_allele(&[1.0, 1.1, 0.0], &[1.0, 0.0, 0.0]);
    assert_eq!((best, second), (0, 1));
    let confidence = allele_confidence(best_l, second_l);
    assert!(confidence < 0.0);
    assert!(!(confidence > LOG_10_INFORMATIVE_THRESHOLD));
}

#[test]
fn search_best_allele_requires_confidence_strictly_above_threshold() {
    let (best, _, best_l, second_l) = search_best_allele(&[0.0, 0.2], &[1.0, 0.0]);
    assert_eq!(best, 1, "a gap of exactly 0.2 is not a near-tie");
    assert_eq!(allele_confidence(best_l, second_l), 0.2);
    assert!(!(allele_confidence(best_l, second_l) > LOG_10_INFORMATIVE_THRESHOLD));

    let (best, _, best_l, second_l) = search_best_allele(&[0.0, 0.2000001], &[1.0, 0.0]);
    assert_eq!(best, 1);
    assert!(allele_confidence(best_l, second_l) > LOG_10_INFORMATIVE_THRESHOLD);
}

#[test]
fn search_best_allele_far_reference_does_not_take_a_distant_tie() {
    // Two non-reference alleles tie at the top. REF is far below that maximum,
    // so its priority does not replace the earlier tied allele.
    let (best, second, best_l, second_l) = search_best_allele(&[0.0, 5.0, 5.0], &[1.0, 0.0, 0.0]);
    assert_eq!((best, second), (1, 2));
    assert_eq!(allele_confidence(best_l, second_l), 0.0);
}

#[test]
fn search_best_allele_equal_likelihood_keeps_earlier_allele() {
    let (best, second, best_l, second_l) = search_best_allele(&[1.0, 1.0, 0.5], &[1.0, 0.0, 0.0]);
    assert_eq!((best, second), (0, 1));
    assert_eq!(allele_confidence(best_l, second_l), 0.0);
}

#[test]
fn non_record_best_allele_is_not_reassigned_to_emitted_alt() {
    assert_eq!(
        annotation_allele_class("TTTG", "T", "TGTTTG"),
        AnnotationAlleleClass::Other
    );
    assert_eq!(
        annotation_allele_class("T", "T", "TGTTTG"),
        AnnotationAlleleClass::EmittedRef
    );
    assert_eq!(
        annotation_allele_class("TGTTTG", "T", "TGTTTG"),
        AnnotationAlleleClass::EmittedAlt
    );
    let (best, _, best_l, second_l) = search_best_allele(&[0.0, 5.0, 0.1], &[1.0, 0.0, 0.0]);
    assert_eq!(best, 1);
    assert!(allele_confidence(best_l, second_l) > LOG_10_INFORMATIVE_THRESHOLD);
}

fn three_decimal(value: f64) -> String {
    format!("{value:.3}")
}

fn optional_three_decimal(value: Option<f64>) -> String {
    value.map(three_decimal).unwrap_or_else(|| ".".to_string())
}

/// Frozen GATK 4.4.0.0 read classifications at 20:29455942 (7R.72).
/// Multiset of (qname, best allele, class, informative, counted, reverse).
const JAVA_29455942_READS: &[(&str, &str, &str, bool, bool, bool)] = &[
    (
        "HWI-D00360:8:H88U0ADXX:1:2104:6643:85211",
        "AAT",
        "ALT",
        true,
        true,
        true,
    ),
    (
        "HWI-D00360:6:H81VLADXX:2:2205:5601:38112",
        "AAT",
        "ALT",
        true,
        true,
        false,
    ),
    (
        "HISEQ1:9:H8962ADXX:2:1212:4076:61181",
        "AAT",
        "ALT",
        true,
        true,
        true,
    ),
    (
        "HISEQ1:9:H8962ADXX:2:1104:8730:7170",
        "ATT",
        "IGNORED",
        true,
        false,
        true,
    ),
    (
        "HWI-D00360:8:H88U0ADXX:2:1206:3145:74907",
        "AAT",
        "ALT",
        true,
        true,
        true,
    ),
    (
        "HWI-D00360:8:H88U0ADXX:2:2116:18111:55930",
        "T",
        "IGNORED",
        true,
        false,
        true,
    ),
    (
        "HWI-D00360:6:H81VLADXX:2:1210:13801:21690",
        "A",
        "REF",
        true,
        true,
        false,
    ),
    (
        "HWI-D00360:6:H81VLADXX:1:1216:15328:67949",
        "AAT",
        "ALT",
        true,
        true,
        false,
    ),
    (
        "HISEQ1:9:H8962ADXX:1:1107:11064:51936",
        "AAT",
        "ALT",
        true,
        true,
        true,
    ),
    (
        "HWI-D00360:8:H88U0ADXX:2:1109:7902:95502",
        "T",
        "IGNORED",
        true,
        false,
        true,
    ),
    (
        "HISEQ1:11:H8GV6ADXX:1:1204:16745:38126",
        "A",
        "REF",
        true,
        true,
        false,
    ),
    (
        "HWI-D00360:5:H814YADXX:1:1212:17430:32056",
        "ATT",
        "IGNORED",
        true,
        false,
        true,
    ),
    (
        "HWI-D00360:5:H814YADXX:2:2107:1177:54891",
        "AAT",
        "ALT",
        true,
        true,
        false,
    ),
    (
        "HISEQ1:12:H8GVUADXX:2:1208:7611:58296",
        "A",
        "REF",
        true,
        true,
        true,
    ),
    (
        "HISEQ1:12:H8GVUADXX:2:2215:6274:41210",
        "AAT",
        "ALT",
        true,
        true,
        true,
    ),
    (
        "HWI-D00360:5:H814YADXX:1:2211:10796:37074",
        "A",
        "REF",
        true,
        true,
        true,
    ),
    (
        "HISEQ1:13:H8G92ADXX:2:1207:9142:14976",
        "ATT",
        "IGNORED",
        true,
        false,
        true,
    ),
    (
        "HWI-D00360:7:H88WKADXX:1:1105:11855:14490",
        "AAT",
        "ALT",
        true,
        true,
        true,
    ),
    (
        "HISEQ1:9:H8962ADXX:1:1213:5437:90850",
        "A",
        "REF",
        true,
        true,
        false,
    ),
    (
        "HISEQ1:13:H8G92ADXX:1:2104:10433:94988",
        "A",
        "REF",
        true,
        true,
        true,
    ),
    (
        "HISEQ1:9:H8962ADXX:1:2214:3946:13940",
        "A",
        "REF",
        true,
        true,
        false,
    ),
    (
        "HWI-D00360:6:H81VLADXX:1:1102:19603:47839",
        "AAT",
        "ALT",
        true,
        true,
        true,
    ),
    (
        "HWI-D00360:8:H88U0ADXX:1:2215:12793:48329",
        "AAT",
        "ALT",
        true,
        true,
        false,
    ),
    (
        "HISEQ1:12:H8GVUADXX:2:2104:18947:36190",
        "AAT",
        "ALT",
        true,
        true,
        true,
    ),
    (
        "HWI-D00360:7:H88WKADXX:2:1105:6275:37858",
        "T",
        "IGNORED",
        true,
        false,
        true,
    ),
    (
        "HISEQ1:9:H8962ADXX:2:1104:19714:91101",
        "A",
        "REF",
        true,
        true,
        false,
    ),
    (
        "HISEQ1:13:H8G92ADXX:2:1204:11285:24052",
        "A",
        "REF",
        true,
        true,
        true,
    ),
    (
        "HISEQ1:9:H8962ADXX:2:1206:9048:57469",
        "AAT",
        "ALT",
        true,
        true,
        true,
    ),
    (
        "HWI-D00360:8:H88U0ADXX:1:1216:6344:27329",
        "ATT",
        "IGNORED",
        true,
        false,
        true,
    ),
    (
        "HWI-D00360:5:H814YADXX:2:2215:4392:84484",
        "A",
        "REF",
        true,
        true,
        false,
    ),
    (
        "HWI-D00360:5:H814YADXX:2:1207:4181:62077",
        "A",
        "REF",
        true,
        true,
        true,
    ),
    (
        "HWI-D00360:6:H81VLADXX:2:1103:4397:32757",
        "ATT",
        "IGNORED",
        true,
        false,
        true,
    ),
    (
        "HISEQ1:9:H8962ADXX:2:1106:13373:82431",
        "A",
        "REF",
        true,
        true,
        true,
    ),
    (
        "HISEQ1:9:H8962ADXX:1:2207:8490:5517",
        "AAT",
        "ALT",
        true,
        true,
        true,
    ),
    (
        "HWI-D00360:5:H814YADXX:1:2207:10890:76583",
        "A",
        "REF",
        true,
        true,
        false,
    ),
    (
        "HISEQ1:12:H8GVUADXX:2:2208:10531:71035",
        "T",
        "IGNORED",
        true,
        false,
        true,
    ),
    (
        "HWI-D00360:5:H814YADXX:1:2115:10943:3938",
        "AAT",
        "ALT",
        true,
        true,
        true,
    ),
    (
        "HWI-D00360:5:H814YADXX:2:1102:14713:45312",
        "A",
        "REF",
        true,
        true,
        true,
    ),
    (
        "HWI-D00360:8:H88U0ADXX:2:1106:13399:2229",
        "ATT",
        "IGNORED",
        true,
        false,
        true,
    ),
    (
        "HWI-D00360:7:H88WKADXX:2:1206:6658:94818",
        "ATT",
        "IGNORED",
        true,
        false,
        true,
    ),
    (
        "HISEQ1:12:H8GVUADXX:1:1110:16456:8464",
        "A",
        "REF",
        true,
        true,
        true,
    ),
    (
        "HWI-D00360:6:H81VLADXX:2:2115:16939:20047",
        "AAT",
        "ALT",
        true,
        true,
        true,
    ),
    (
        "HWI-D00360:5:H814YADXX:2:1205:3134:10685",
        "A",
        "REF",
        true,
        true,
        true,
    ),
    (
        "HISEQ1:13:H8G92ADXX:2:2102:1356:32269",
        "A",
        "REF",
        true,
        true,
        false,
    ),
    (
        "HISEQ1:11:H8GV6ADXX:1:2205:11286:28668",
        "ATT",
        "IGNORED",
        true,
        false,
        true,
    ),
    (
        "HWI-D00360:8:H88U0ADXX:2:2105:9169:69257",
        "T",
        "IGNORED",
        true,
        false,
        true,
    ),
    (
        "HISEQ1:12:H8GVUADXX:1:2114:19696:40833",
        "ATT",
        "IGNORED",
        true,
        false,
        true,
    ),
    (
        "HWI-D00360:8:H88U0ADXX:1:1207:4297:19722",
        "ATT",
        "IGNORED",
        true,
        false,
        true,
    ),
    (
        "HISEQ1:11:H8GV6ADXX:1:2208:19582:60012",
        "AAT",
        "ALT",
        true,
        true,
        true,
    ),
    (
        "HISEQ1:9:H8962ADXX:1:1116:1737:10125",
        "AAT",
        "ALT",
        true,
        true,
        true,
    ),
    (
        "HWI-D00360:6:H81VLADXX:2:2111:17451:73433",
        "ATT",
        "IGNORED",
        true,
        false,
        true,
    ),
    (
        "HISEQ1:9:H8962ADXX:1:2202:18336:27346",
        "A",
        "REF",
        true,
        true,
        false,
    ),
    (
        "HWI-D00360:5:H814YADXX:1:2106:9068:58347",
        "AAT",
        "ALT",
        true,
        true,
        false,
    ),
    (
        "HWI-D00360:6:H81VLADXX:2:1103:9648:12748",
        "AAT",
        "ALT",
        true,
        true,
        false,
    ),
    (
        "HWI-D00360:6:H81VLADXX:1:2203:12713:95012",
        "AAT",
        "ALT",
        true,
        true,
        true,
    ),
    (
        "HISEQ1:9:H8962ADXX:2:1206:11075:50547",
        "A",
        "REF",
        true,
        true,
        false,
    ),
    (
        "HISEQ1:12:H8GVUADXX:1:2113:19174:50109",
        "A",
        "REF",
        true,
        true,
        true,
    ),
    (
        "HWI-D00360:5:H814YADXX:2:1113:5150:90271",
        "AAT",
        "ALT",
        true,
        true,
        true,
    ),
    (
        "HISEQ1:11:H8GV6ADXX:2:2209:12379:22070",
        "A",
        "REF",
        true,
        true,
        false,
    ),
    (
        "HWI-D00360:5:H814YADXX:2:1102:8133:87063",
        "A",
        "REF",
        true,
        true,
        false,
    ),
    (
        "HISEQ1:13:H8G92ADXX:2:1207:5295:91465",
        "ATT",
        "IGNORED",
        true,
        false,
        true,
    ),
    (
        "HWI-D00360:5:H814YADXX:1:1107:10132:8763",
        "AAT",
        "ALT",
        true,
        true,
        true,
    ),
    (
        "HWI-D00360:5:H814YADXX:1:2216:10487:39079",
        "ATT",
        "IGNORED",
        true,
        false,
        true,
    ),
    (
        "HISEQ1:11:H8GV6ADXX:1:1104:4186:41128",
        "A",
        "REF",
        true,
        true,
        false,
    ),
    (
        "HISEQ1:9:H8962ADXX:2:2201:14338:36616",
        "AT",
        "IGNORED",
        true,
        false,
        true,
    ),
    (
        "HISEQ1:11:H8GV6ADXX:2:2209:12379:22070",
        "A",
        "REF",
        true,
        true,
        true,
    ),
    (
        "HISEQ1:11:H8GV6ADXX:1:2116:18670:99941",
        "ATT",
        "IGNORED",
        true,
        false,
        true,
    ),
    (
        "HWI-D00360:6:H81VLADXX:1:1216:15328:67949",
        "AAT",
        "ALT",
        true,
        true,
        true,
    ),
    (
        "HISEQ1:11:H8GV6ADXX:1:2104:1513:33781",
        "A",
        "REF",
        true,
        true,
        true,
    ),
    (
        "HWI-D00360:5:H814YADXX:1:1107:5408:90614",
        "A",
        "REF",
        true,
        true,
        true,
    ),
    (
        "HWI-D00360:7:H88WKADXX:2:1206:11064:69408",
        "A",
        "REF",
        true,
        true,
        true,
    ),
    (
        "HWI-D00360:5:H814YADXX:1:2106:9068:58347",
        "AAT",
        "ALT",
        true,
        true,
        true,
    ),
    (
        "HWI-D00360:6:H81VLADXX:1:2101:16091:54483",
        "A",
        "REF",
        true,
        true,
        false,
    ),
    (
        "HISEQ1:11:H8GV6ADXX:1:1203:17733:55670",
        "A",
        "REF",
        true,
        true,
        true,
    ),
    (
        "HWI-D00360:6:H81VLADXX:1:1102:16926:25984",
        "A",
        "REF",
        true,
        true,
        true,
    ),
    (
        "HWI-D00360:8:H88U0ADXX:2:2101:7566:46476",
        "AAT",
        "ALT",
        true,
        true,
        true,
    ),
    (
        "HISEQ1:9:H8962ADXX:1:1106:10350:18258",
        "AAT",
        "ALT",
        true,
        true,
        true,
    ),
    (
        "HWI-D00360:5:H814YADXX:2:2107:14652:60528",
        "A",
        "REF",
        true,
        true,
        true,
    ),
    (
        "HWI-D00360:7:H88WKADXX:1:2103:7800:34534",
        "AAT",
        "ALT",
        true,
        true,
        true,
    ),
    (
        "HISEQ1:9:H8962ADXX:1:1110:12402:16069",
        "AAT",
        "ALT",
        true,
        true,
        true,
    ),
    (
        "HWI-D00360:8:H88U0ADXX:2:1105:14270:89425",
        "AAT",
        "ALT",
        true,
        true,
        true,
    ),
    (
        "HISEQ1:13:H8G92ADXX:1:2211:12679:32928",
        "AT",
        "IGNORED",
        true,
        false,
        true,
    ),
    (
        "HWI-D00360:8:H88U0ADXX:1:1202:9544:72702",
        "A",
        "REF",
        true,
        true,
        false,
    ),
    (
        "HWI-D00360:6:H81VLADXX:1:1104:11023:44992",
        "A",
        "REF",
        true,
        true,
        true,
    ),
    (
        "HISEQ1:12:H8GVUADXX:1:2115:4673:74518",
        "A",
        "REF",
        true,
        true,
        false,
    ),
    (
        "HWI-D00360:6:H81VLADXX:1:2105:13991:32173",
        "AAT",
        "ALT",
        true,
        true,
        true,
    ),
    (
        "HWI-D00360:5:H814YADXX:2:1211:2735:67553",
        "A",
        "REF",
        true,
        true,
        false,
    ),
    (
        "HWI-D00360:6:H81VLADXX:2:1215:11178:56038",
        "ATT",
        "IGNORED",
        true,
        false,
        true,
    ),
    (
        "HISEQ1:13:H8G92ADXX:1:2215:3048:31329",
        "A",
        "REF",
        true,
        true,
        false,
    ),
    (
        "HWI-D00360:5:H814YADXX:2:1202:19283:18749",
        "A",
        "REF",
        true,
        true,
        false,
    ),
    (
        "HISEQ1:13:H8G92ADXX:2:2106:15182:80281",
        "AAT",
        "ALT",
        true,
        true,
        true,
    ),
    (
        "HISEQ1:12:H8GVUADXX:1:2102:3262:64140",
        "A",
        "REF",
        true,
        true,
        false,
    ),
    (
        "HISEQ1:11:H8GV6ADXX:1:1116:4033:65919",
        "AAT",
        "ALT",
        true,
        true,
        true,
    ),
    (
        "HWI-D00360:7:H88WKADXX:1:2111:7437:40887",
        "ATT",
        "IGNORED",
        true,
        false,
        true,
    ),
    (
        "HISEQ1:12:H8GVUADXX:1:1201:10428:98142",
        "AAT",
        "ALT",
        true,
        true,
        true,
    ),
    (
        "HWI-D00360:8:H88U0ADXX:2:2214:8541:27804",
        "ATT",
        "IGNORED",
        true,
        false,
        true,
    ),
    (
        "HISEQ1:9:H8962ADXX:2:2212:1195:72188",
        "A",
        "REF",
        true,
        true,
        true,
    ),
    (
        "HWI-D00360:5:H814YADXX:1:2215:1427:64469",
        "ATT",
        "IGNORED",
        true,
        false,
        true,
    ),
    (
        "HISEQ1:11:H8GV6ADXX:1:1210:3128:10108",
        "A",
        "REF",
        true,
        true,
        true,
    ),
    (
        "HWI-D00360:6:H81VLADXX:2:1104:15554:2818",
        "A",
        "REF",
        true,
        true,
        false,
    ),
    (
        "HWI-D00360:7:H88WKADXX:2:1105:8086:83865",
        "AAT",
        "ALT",
        true,
        true,
        true,
    ),
    (
        "HWI-D00360:7:H88WKADXX:2:2209:12675:63734",
        "AAT",
        "ALT",
        true,
        true,
        false,
    ),
    (
        "HWI-D00360:7:H88WKADXX:1:1209:8005:92936",
        "T",
        "UNINFO",
        false,
        false,
        true,
    ),
    (
        "HWI-D00360:7:H88WKADXX:1:2211:9082:22140",
        "A",
        "REF",
        true,
        true,
        true,
    ),
    (
        "HISEQ1:9:H8962ADXX:1:1105:6235:77999",
        "ATT",
        "IGNORED",
        true,
        false,
        true,
    ),
    (
        "HISEQ1:13:H8G92ADXX:2:1209:2537:66911",
        "A",
        "REF",
        true,
        true,
        false,
    ),
    (
        "HWI-D00360:7:H88WKADXX:2:2102:14743:6741",
        "A",
        "REF",
        true,
        true,
        false,
    ),
    (
        "HWI-D00360:7:H88WKADXX:1:1202:17988:48058",
        "AAT",
        "ALT",
        true,
        true,
        true,
    ),
    (
        "HISEQ1:11:H8GV6ADXX:1:1201:12241:62342",
        "A",
        "REF",
        true,
        true,
        false,
    ),
    (
        "HISEQ1:13:H8G92ADXX:2:1210:1858:43133",
        "A",
        "REF",
        true,
        true,
        false,
    ),
    (
        "HWI-D00360:5:H814YADXX:2:1108:10112:66562",
        "AAT",
        "ALT",
        true,
        true,
        true,
    ),
    (
        "HISEQ1:12:H8GVUADXX:1:2215:18324:5571",
        "T",
        "UNINFO",
        false,
        false,
        true,
    ),
    (
        "HISEQ1:13:H8G92ADXX:2:1112:5418:11026",
        "A",
        "REF",
        true,
        true,
        false,
    ),
];

#[test]
fn shared_annotation_selection_matches_java_at_two_loci() {
    let root = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("..");
    let ref_fasta = root.join("parity/realworld/assets/hs37d5.simple.fa");
    let bam = root.join("parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam");
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing chr20 tiny reference or BAM");
        return;
    }
    let dict =
        gatk_core::reference::SequenceDictionary::from_fasta_path(&ref_fasta).expect("dictionary");
    let specs = gatk_core::reference::parse_intervals_cli_string(&dict, "20:29455000-29456500")
        .expect("interval");
    let walk = crate::traverse_assembly_region_walker(
        &dict,
        &specs,
        &ref_fasta,
        &bam,
        &crate::ReadFilterParams::gatk_standard_hc(),
        &crate::WalkerTraversalConfig::gatk_haplotype_caller_production(100),
    )
    .expect("walk");
    let regions = crate::flatten_assembly_regions(&walk);
    check_locus_29455649(&dict, &ref_fasta, &regions);
    check_locus_29455942(&dict, &ref_fasta, &regions);
    check_locus_29456344(&dict, &ref_fasta, &regions);
}

fn active_region<'a>(regions: &'a [crate::AssemblyRegion], loc: u64) -> &'a crate::AssemblyRegion {
    regions
        .iter()
        .find(|region| {
            matches!(
                crate::call_disposition(region),
                crate::AssemblyRegionCallDisposition::ActiveFull
            ) && region.start.get() <= loc
                && region.end.get() >= loc
        })
        .unwrap_or_else(|| panic!("no ActiveFull region covers {loc}"))
}

fn evidence_for<'a>(
    region: &'a crate::AssemblyRegion,
    outcome: &'a crate::CallRegionOutcome,
    call: &'a crate::hc_genotyping_engine::GenotypedSiteCall,
) -> HcStrandBiasLikelihoods<'a> {
    let (full_ref, full_pad) = outcome.assembly.event_map_reference();
    let apply_pad = outcome
        .assembly
        .haplotypes
        .iter()
        .find(|haplotype| haplotype.is_reference)
        .and_then(|haplotype| haplotype.genome_loc.as_ref().map(|loc| loc.start.get()))
        .unwrap_or(full_pad);
    let likelihoods = if call.annotation_likelihoods.is_empty() {
        outcome.read_likelihoods.as_slice()
    } else {
        call.annotation_likelihoods.as_slice()
    };
    HcStrandBiasLikelihoods {
        reads: &outcome.genotyping_reads,
        likelihoods,
        haplotypes: &outcome.assembly.haplotypes,
        contig: &region.contig,
        ref_bytes: outcome.assembly.reference_bases(),
        pad_start_1based: apply_pad,
        full_ref_bytes: full_ref,
        full_pad_1based: full_pad,
        max_mnp_distance: outcome.assembly.max_mnp_distance(),
        emit_spanning_dels: !HcGenotypingConfig::strict_java().disable_spanning_event_genotyping,
    }
}

fn call_emitted<'a>(
    outcome: &'a crate::CallRegionOutcome,
    loc: u64,
    ref_allele: &str,
    alt_allele: &str,
) -> &'a crate::hc_genotyping_engine::GenotypedSiteCall {
    outcome
        .genotyped_calls
        .iter()
        .find(|call| {
            call.event.start_1based.get() == loc
                && call.event.ref_allele == ref_allele
                && call.event.alt_allele == alt_allele
        })
        .unwrap_or_else(|| {
            let nearby = outcome
                .genotyped_calls
                .iter()
                .filter(|call| call.event.start_1based.get().abs_diff(loc) < 30)
                .map(|call| {
                    format!(
                        "{} {}/{}",
                        call.event.start_1based.get(),
                        call.event.ref_allele,
                        call.event.alt_allele
                    )
                })
                .collect::<Vec<_>>();
            panic!("missing {loc} {ref_allele}/{alt_allele}; nearby {nearby:?}")
        })
}

fn assert_candidate_pools(
    evidence: &HcStrandBiasLikelihoods<'_>,
    loc: u64,
    expected_ref: &str,
    expected_alts: &[&str],
    expected_pool_sizes: &[usize],
) {
    let (long_ref, alts) = annotation_candidate_alleles(evidence, loc)
        .unwrap_or_else(|| panic!("{loc} did not recover a multi-allele candidate set"));
    assert_eq!(long_ref, expected_ref, "{loc} reference allele");
    assert_eq!(alts, expected_alts, "{loc} candidate alt order");
    let cache = build_per_haplotype_variation_events(
        evidence.haplotypes,
        evidence.full_ref_bytes,
        evidence.full_pad_1based,
        evidence.max_mnp_distance,
        evidence.contig,
    );
    let (pools, _, _) = java_create_allele_mapper_pools(
        evidence.haplotypes.len(),
        &cache,
        loc,
        &long_ref,
        &alts,
        evidence.emit_spanning_dels,
    );
    let sizes = pools.iter().map(Vec::len).collect::<Vec<_>>();
    assert_eq!(sizes, expected_pool_sizes, "{loc} haplotype pool sizes");
}

fn assert_genotype_contract(
    call: &crate::hc_genotyping_engine::GenotypedSiteCall,
    pl: &[i32],
    ad: &[i32],
    dp: i32,
    qual: &str,
    ann_qual: f64,
) {
    assert_eq!(call.genotype.format.pl_as_i32(), pl);
    assert_eq!(call.genotype.format.ad_as_i32(), ad);
    assert_eq!(call.genotype.format.dp.as_i32(), dp);
    let best = call
        .genotype
        .format
        .pl_as_i32()
        .iter()
        .enumerate()
        .min_by_key(|(_, value)| *value)
        .map(|(index, _)| index);
    assert_eq!(best, Some(1), "GT 0/1 is the minimum PL");
    assert_eq!(format!("{ann_qual:.2}"), qual);
}

fn attached_baseq_lists(
    evidence: &HcStrandBiasLikelihoods<'_>,
    loc: u64,
    ref_allele: &str,
    alt_allele: &str,
) -> (Vec<f64>, Vec<f64>) {
    let mut ref_quals = Vec::new();
    let mut alt_quals = Vec::new();
    for decision in annotation_allele_decisions(evidence, loc, ref_allele, alt_allele, None) {
        if !decision.informative {
            continue;
        }
        let Some(rec) = evidence.reads.get(decision.read_index) else {
            continue;
        };
        if !java_ranksum_mq_usable(rec) {
            continue;
        }
        let Some(quality) = read_base_quality_at_ref_coord_1based(rec, loc as i32) else {
            continue;
        };
        match decision.class {
            AnnotationAlleleClass::EmittedRef => ref_quals.push(f64::from(quality)),
            AnnotationAlleleClass::EmittedAlt => alt_quals.push(f64::from(quality)),
            AnnotationAlleleClass::Other => {}
        }
    }
    (ref_quals, alt_quals)
}

fn check_locus_29455649(
    dict: &gatk_core::reference::SequenceDictionary,
    ref_fasta: &std::path::Path,
    regions: &[crate::AssemblyRegion],
) {
    const LOC: u64 = 29_455_649;
    let region = active_region(regions, LOC);
    let outcome = crate::HaplotypeCallerEngine::call_region(
        region,
        dict,
        ref_fasta,
        &crate::CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("outcome");
    let call = call_emitted(&outcome, LOC, "T", "TGTTTG");
    let evidence = evidence_for(region, &outcome, call);
    assert_candidate_pools(&evidence, LOC, "T", &["TTTG", "TGTTTG"], &[14, 5, 5]);
    let table = strand_bias_sample_counts(&evidence, LOC, "T", "TGTTTG", None);
    assert_eq!(table, (60, 23, 10, 9));
    let fs = fisher_strand::fisher_strand_statistic(table.0, table.1, table.2, table.3);
    assert_eq!(fs.to_bits(), 9.677494113791163_f64.to_bits());
    let ann = annotate_hc_variant_site(
        Some(region),
        LOC,
        "T",
        "TGTTTG",
        &call.genotype,
        &HcGenotypingConfig::strict_java(),
        call.qual_log10_p_error,
        Some(&evidence),
        None,
    )
    .expect("annotate");
    assert_eq!(ann.fs.to_bits(), fs.to_bits());
    assert_eq!(three_decimal(ann.fs), "9.677");
    let (diag_ref, diag_alt) = attached_baseq_lists(&evidence, LOC, "T", "TGTTTG");
    assert_eq!((diag_ref.len(), diag_alt.len()), (83, 18));
    let diag_z = rank_sum_baseq::base_quality_rank_sum(&diag_ref, &diag_alt).expect("diag z");
    assert_eq!(diag_z.to_bits(), (-0.1729717140331544_f64).to_bits());
    let (prod_ref, prod_alt) = baseq_rank_sum_quals_from_likelihoods(&evidence, LOC, "T", "TGTTTG");
    assert_eq!(
        (prod_ref.len(), prod_alt.len()),
        (83, 18),
        "production BaseQ is the attached alignment lookup"
    );
    let prod_z = rank_sum_baseq::base_quality_rank_sum(&prod_ref, &prod_alt).expect("prod z");
    assert_eq!(prod_z.to_bits(), (-0.1729717140331544_f64).to_bits());
    assert_eq!(three_decimal(prod_z), "-0.173");
    assert_eq!(optional_three_decimal(ann.base_q_rank_sum), "-0.173");
    assert_genotype_contract(call, &[570, 0, 3517], &[88, 22], 110, "562.60", ann.qual);
    assert_eq!(
        ann.read_pos_rank_sum.map(f64::to_bits),
        Some(0.8124396993400359_f64.to_bits())
    );
    assert_eq!(optional_three_decimal(ann.read_pos_rank_sum), "0.812");
    assert_eq!(optional_three_decimal(ann.mq_rank_sum), "-3.705");
    assert_eq!(three_decimal(ann.sor), "0.172");
}

fn check_locus_29455942(
    dict: &gatk_core::reference::SequenceDictionary,
    ref_fasta: &std::path::Path,
    regions: &[crate::AssemblyRegion],
) {
    const LOC: u64 = 29_455_942;
    let region = active_region(regions, LOC);
    let outcome = crate::HaplotypeCallerEngine::call_region(
        region,
        dict,
        ref_fasta,
        &crate::CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("outcome");
    let call = call_emitted(&outcome, LOC, "A", "AAT");
    let evidence = evidence_for(region, &outcome, call);
    assert_candidate_pools(
        &evidence,
        LOC,
        "A",
        &["ATT", "T", "AT", "AAT"],
        &[32, 10, 5, 7, 22],
    );
    let table = strand_bias_sample_counts(&evidence, LOC, "A", "AAT", None);
    assert_eq!(table, (26, 21, 7, 31));
    let fs = fisher_strand::fisher_strand_statistic(table.0, table.1, table.2, table.3);
    assert_eq!(fs.to_bits(), 31.505493249597485_f64.to_bits());
    let (long_ref, alts) = annotation_candidate_alleles(&evidence, LOC).expect("alleles");
    let mut allele_names = vec![long_ref];
    allele_names.extend(alts);
    let mut actual = Vec::new();
    for decision in annotation_allele_decisions(&evidence, LOC, "A", "AAT", None) {
        let rec = evidence
            .reads
            .get(decision.read_index)
            .unwrap_or_else(|| panic!("missing read {}", decision.read_index));
        let qname = String::from_utf8_lossy(rec.qname()).into_owned();
        let best = allele_names
            .get(decision.best_index)
            .map(String::as_str)
            .unwrap_or("?")
            .to_string();
        let class = if !decision.informative {
            "UNINFO"
        } else {
            match decision.class {
                AnnotationAlleleClass::EmittedRef => "REF",
                AnnotationAlleleClass::EmittedAlt => "ALT",
                AnnotationAlleleClass::Other => "IGNORED",
            }
        };
        let counted = decision.informative
            && matches!(
                decision.class,
                AnnotationAlleleClass::EmittedRef | AnnotationAlleleClass::EmittedAlt
            );
        let reverse = rec.flags() & FLAG_REVERSE != 0;
        actual.push((
            qname,
            best,
            class.to_string(),
            decision.informative,
            counted,
            reverse,
        ));
    }
    let mut expected = JAVA_29455942_READS
        .iter()
        .map(|(qname, best, class, info, counted, reverse)| {
            (
                (*qname).to_string(),
                (*best).to_string(),
                (*class).to_string(),
                *info,
                *counted,
                *reverse,
            )
        })
        .collect::<Vec<_>>();
    expected.sort();
    actual.sort();
    assert_eq!(actual.len(), 113);
    assert_eq!(actual, expected);
    let ann = annotate_hc_variant_site(
        Some(region),
        LOC,
        "A",
        "AAT",
        &call.genotype,
        &HcGenotypingConfig::strict_java(),
        call.qual_log10_p_error,
        Some(&evidence),
        None,
    )
    .expect("annotate");
    assert_eq!(ann.fs.to_bits(), fs.to_bits());
    assert_eq!(three_decimal(ann.fs), "31.505");
    assert_genotype_contract(call, &[1558, 0, 1842], &[48, 44], 92, "1550.60", ann.qual);
    let (diag_ref, diag_alt) = attached_baseq_lists(&evidence, LOC, "A", "AAT");
    assert_eq!((diag_ref.len(), diag_alt.len()), (46, 36));
    let diag_z = rank_sum_baseq::base_quality_rank_sum(&diag_ref, &diag_alt).expect("diag z");
    assert_eq!(diag_z.to_bits(), (-0.22673512728192455_f64).to_bits());
    assert_eq!(three_decimal(diag_z), "-0.227");
    let (prod_ref, prod_alt) = baseq_rank_sum_quals_from_likelihoods(&evidence, LOC, "A", "AAT");
    let prod_z = rank_sum_baseq::base_quality_rank_sum(&prod_ref, &prod_alt).expect("prod z");
    assert_eq!(
        prod_z.to_bits(),
        ann.base_q_rank_sum.expect("prod baseq").to_bits()
    );
    assert_eq!((prod_ref.len(), prod_alt.len()), (46, 36));
    assert_eq!(prod_z.to_bits(), (-0.22673512728192455_f64).to_bits());
    assert_eq!(three_decimal(prod_z), "-0.227");
    assert_eq!(
        ann.read_pos_rank_sum.map(f64::to_bits),
        Some((-0.3803560137447075_f64).to_bits())
    );
    assert_eq!(optional_three_decimal(ann.read_pos_rank_sum), "-0.380");
    assert_eq!(optional_three_decimal(ann.mq_rank_sum), "-6.667");
    assert_eq!(three_decimal(ann.sor), "2.813");
}

fn check_locus_29456344(
    dict: &gatk_core::reference::SequenceDictionary,
    ref_fasta: &std::path::Path,
    regions: &[crate::AssemblyRegion],
) {
    const LOC: u64 = 29_456_344;
    let region = active_region(regions, LOC);
    let outcome = crate::HaplotypeCallerEngine::call_region(
        region,
        dict,
        ref_fasta,
        &crate::CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("outcome");
    let call = call_emitted(&outcome, LOC, "T", "C");
    assert!(call.extra_alt_alleles.is_empty(), "VCF emits T/C");
    let (pre_ref, kept_alts) = call
        .annotation_output_alleles
        .as_ref()
        .expect("pre-trim TG/CG retained");
    assert_eq!(pre_ref, "TG");
    assert_eq!(kept_alts.as_slice(), ["CG".to_string()]);
    let evidence = evidence_for(region, &outcome, call);
    let output = Some((pre_ref.as_str(), kept_alts.as_slice()));
    let decisions = annotation_allele_decisions(&evidence, LOC, "T", "C", output);
    let mut ref_n = 0usize;
    let mut alt_n = 0usize;
    let mut other_n = 0usize;
    let mut uninformative = 0usize;
    for decision in &decisions {
        if !decision.informative {
            uninformative += 1;
            continue;
        }
        match decision.class {
            AnnotationAlleleClass::EmittedRef => ref_n += 1,
            AnnotationAlleleClass::EmittedAlt => alt_n += 1,
            AnnotationAlleleClass::Other => other_n += 1,
        }
    }
    assert_eq!(decisions.len(), 60);
    assert_eq!((ref_n, alt_n, other_n, uninformative), (34, 17, 7, 2));
    let ann = annotate_hc_variant_site(
        Some(region),
        LOC,
        "T",
        "C",
        &call.genotype,
        &HcGenotypingConfig::strict_java(),
        call.qual_log10_p_error,
        Some(&evidence),
        output,
    )
    .expect("annotate");
    assert_eq!(
        ann.base_q_rank_sum.map(f64::to_bits),
        Some((-3.998014443084962_f64).to_bits())
    );
    assert_eq!(optional_three_decimal(ann.base_q_rank_sum), "-3.998");
    assert_eq!(
        ann.mq_rank_sum.map(f64::to_bits),
        Some((-4.8940920751716872_f64).to_bits())
    );
    assert_eq!(optional_three_decimal(ann.mq_rank_sum), "-4.894");
    assert_eq!(
        ann.read_pos_rank_sum.map(f64::to_bits),
        Some(0.67837048764109575_f64.to_bits())
    );
    assert_eq!(optional_three_decimal(ann.read_pos_rank_sum), "0.678");
    assert_eq!(ann.sor.to_bits(), 0.7929925155296615_f64.to_bits());
    assert_eq!(three_decimal(ann.sor), "0.793");
    assert_genotype_contract(call, &[542, 0, 1353], &[36, 19], 55, "510.06", ann.qual);
}

fn record_with_cigar(pos0: i64, ops: &[Cigar], qual: &[u8]) -> Record {
    let qlen: usize = ops
        .iter()
        .map(|op| match op {
            Cigar::Match(n)
            | Cigar::Equal(n)
            | Cigar::Diff(n)
            | Cigar::Ins(n)
            | Cigar::SoftClip(n) => *n as usize,
            _ => 0,
        })
        .sum();
    let seq = vec![b'A'; qlen];
    let qualities = if qual.len() == qlen {
        qual.to_vec()
    } else {
        vec![30u8; qlen]
    };
    let mut rec = Record::new();
    rec.set(
        b"readpos-edge",
        Some(&CigarString::from(ops.to_vec())),
        &seq,
        &qualities,
    );
    rec.set_pos(pos0);
    rec.set_mapq(60);
    rec
}

#[test]
fn read_pos_leading_insertion_returns_zero_only_on_the_variant_end_edge() {
    let adjacent_start0 = 99i64;
    let cases: &[(&str, &[Cigar], i64, i32, i32, Option<f64>)] = &[
        (
            "5I79M64H immediately after the variant",
            &[Cigar::Ins(5), Cigar::Match(79), Cigar::HardClip(64)],
            adjacent_start0,
            99,
            99,
            Some(0.0),
        ),
        (
            "2I79M67H immediately after the variant",
            &[Cigar::Ins(2), Cigar::Match(79), Cigar::HardClip(67)],
            adjacent_start0,
            99,
            99,
            Some(0.0),
        ),
        (
            "leading insertion whose start is the variant end, not the next base",
            &[Cigar::Ins(5), Cigar::Match(79), Cigar::HardClip(64)],
            adjacent_start0,
            100,
            100,
            Some(5.0),
        ),
        (
            "leading insertion that does not meet the variant",
            &[Cigar::Ins(5), Cigar::Match(79), Cigar::HardClip(64)],
            adjacent_start0,
            50,
            50,
            None,
        ),
        (
            "hard clip then insertion on the next base",
            &[Cigar::HardClip(3), Cigar::Ins(5), Cigar::Match(10)],
            adjacent_start0,
            99,
            99,
            Some(0.0),
        ),
        (
            "soft clip then insertion on the next base",
            &[Cigar::SoftClip(3), Cigar::Ins(5), Cigar::Match(10)],
            adjacent_start0,
            99,
            99,
            Some(0.0),
        ),
        (
            "soft clip then match on the next base",
            &[Cigar::SoftClip(4), Cigar::Match(10)],
            adjacent_start0,
            100,
            99,
            Some(4.0),
        ),
        (
            "match is the first non-clipping operator",
            &[Cigar::Match(10)],
            adjacent_start0,
            105,
            99,
            Some(4.0),
        ),
        (
            "insertion after a match does not take the edge",
            &[Cigar::Match(10), Cigar::Ins(5)],
            adjacent_start0,
            105,
            99,
            Some(5.0),
        ),
    ];
    for (label, ops, pos0, vc_start, vc_end, expected) in cases {
        let rec = record_with_cigar(*pos0, ops, &[]);
        assert_eq!(
            java_read_pos_rank_element(&rec, *vc_start, *vc_end),
            *expected,
            "{label}"
        );
    }

    let long_ref_end = java_variant_context_end_1based(100, "AAAAAA");
    assert_eq!(long_ref_end, 105);
    let after_end = record_with_cigar(105, &[Cigar::Ins(2), Cigar::Match(10)], &[]);
    assert_eq!(
        java_read_pos_rank_element(&after_end, 100, long_ref_end),
        Some(0.0)
    );
    let after_start_only = record_with_cigar(100, &[Cigar::Ins(2), Cigar::Match(10)], &[]);
    assert_eq!(
        java_read_pos_rank_element(&after_start_only, 100, long_ref_end),
        None,
        "start + 1 is not the variant end when the reference allele is longer than one base"
    );
}

#[test]
fn read_pos_edge_does_not_change_baseq_pre_realign_fallback() {
    let ops = [Cigar::Ins(5), Cigar::Match(10)];
    let mut qual = vec![30u8; 15];
    qual[9] = 42;
    let mut rec = record_with_cigar(99, &ops, &qual);
    let original_cigar = "15M".to_string();
    rec.push_aux(
        crate::fragment_overlap::RANK_SUM_ORIGINAL_POS_TAG,
        rust_htslib::bam::record::Aux::I32(90),
    )
    .expect("OP");
    rec.push_aux(
        crate::fragment_overlap::RANK_SUM_ORIGINAL_CIGAR_TAG,
        rust_htslib::bam::record::Aux::String(&original_cigar),
    )
    .expect("OC");

    assert_eq!(
        java_read_pos_rank_element(&rec, 99, 99),
        Some(0.0),
        "adjacent leading insertion returns 0 before the OP/OC position"
    );
    assert_eq!(
        java_read_pos_rank_element(&rec, 99, 100),
        Some(5.0),
        "without the edge, ReadPos still uses the pre-realign CIGAR"
    );
    assert_eq!(
        java_baseq_rank_element(&rec, 99),
        None,
        "BaseQ uses the attached alignment and does not recover OP/OC when that alignment misses the coordinate"
    );
    assert_eq!(
        java_baseq_rank_element(&rec, 104),
        Some(42.0),
        "BaseQ returns the attached quality when the attached alignment covers the coordinate"
    );
}
