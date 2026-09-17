//! 6R.200: proof-only (frozen). At `2:92307403 C/A` the attached object is n=6.
//! Java `Coverage.annotate` is `likelihoods.evidenceCount()` with no overlap.
//! The second `java_alignment_read_overlaps_interval(±2)` filter that used to
//! drop two stored rows (6 → 4) was proven here and removed by 6R.201.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! PRODUCTION CHANGE: NONE (6R.200). Coverage cardinality is 6R.201.
//!
//! Live pins keep the overlap *diagnostic* (those two rows still fail ±2) and
//! expect Coverage / INFO DP = 6 after 6R.201.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r200_coverage_no_secondary_overlap_filter -- --nocapture --test-threads=1
//! HOLDOUT_6R200=1 cargo test -p gatk-haplotypecaller --test holdout_6r200_coverage_secondary_overlap -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::hc_genotyping_engine::java_alignment_read_overlaps_interval;
use gatk_haplotypecaller::variant_site_hc_annotations::{
    coverage_evidence_count, rms_mapping_quality_raw, rms_mapping_quality_sample_mapqs,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use rust_htslib::bam::Record;
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
/// Pinned GATK 4.4.0.0 `Coverage.annotate` (SHA `2dbc0258`). No overlap predicate.
const JAVA_COVERAGE_ANNOTATE: &str = r#"
if (likelihoods == null || likelihoods.evidenceCount() == 0) {
    return Collections.emptyMap();
}
final int depth = likelihoods.evidenceCount();
return Collections.singletonMap(getKeyNames().get(0), String.format("%d", depth));
"#;
/// Pinned `AlleleLikelihoods.evidenceCount()` (SHA `2dbc0258`).
const JAVA_EVIDENCE_COUNT: &str =
    "return evidenceBySampleIndex.stream().mapToInt(List::size).sum();";
const GAP_INTERVAL: &str = "2:92305500-92305850";
const TG_INTERVAL: &str = "2:92307200-92307550";
const MID_INTERVAL: &str = "2:92316200-92316580";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_307_403;
const TARGET_REF: &str = "C";
const TARGET_ALT: &str = "A";
const CLOSED_GT: u64 = 92_305_634;
const CLOSED_AG: u64 = 92_305_635;
const CLOSED_AC: u64 = 92_305_716;
const CLOSED_INDEL: u64 = 92_307_324;
const CLOSED_TG: u64 = 92_307_333;
const CLOSED_QUAL: u64 = 92_307_359;
const CLOSED_MID: u64 = 92_316_347;
const MARGIN: i32 = 2;
const EXTRA_QNAME_A: &str = "H06HDADXX130110:1:1101:10061:17286";
const EXTRA_FLAG_A: u16 = 83;
const EXTRA_QNAME_B: &str = "H06JUADXX130110:1:1101:10011:51168";
const EXTRA_FLAG_B: u16 = 81;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R200\t{key}\t{}", value.as_ref());
}

fn qname(rec: &Record) -> String {
    String::from_utf8_lossy(rec.qname()).into_owned()
}

fn unique_indices(likelihoods: &[gatk_haplotypecaller::RegionReadLikelihood]) -> BTreeSet<usize> {
    likelihoods.iter().map(|c| c.read_index.get()).collect()
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

fn has_member(
    reads: &[gatk_haplotypecaller::SharedBamRecord],
    idxs: &BTreeSet<usize>,
    name: &str,
    flag: u16,
) -> bool {
    idxs.iter().any(|idx| {
        let rec: &Record = reads[*idx].as_ref();
        rec.flags() == flag && qname(rec) == name
    })
}

fn rec_at<'a>(
    emitted: &'a [gatk_core::io::vcf::VcfRecord],
    pos: u64,
    r: &str,
) -> &'a gatk_core::io::vcf::VcfRecord {
    emitted
        .iter()
        .find(|rec| rec.position == pos && rec.reference == r)
        .unwrap_or_else(|| panic!("emit {pos} {r}"))
}

fn dump_reads(
    label: &str,
    reads: &[gatk_haplotypecaller::SharedBamRecord],
    idxs: &BTreeSet<usize>,
) {
    for idx in idxs {
        let rec: &Record = reads[*idx].as_ref();
        kv(
            label,
            format!(
                "idx={idx} qname={} flag={} mapq={} start={} end={} cigar={}",
                qname(rec),
                rec.flags(),
                rec.mapq(),
                rec.pos() + 1,
                gatk_haplotypecaller::read_unclip::alignment_end_1based(rec),
                rec.cigar()
            ),
        );
    }
}

#[test]
fn forensic_6r200_java_vs_rust_coverage_contracts() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    kv(
        "classification",
        "ANNOTATION_INPUT_MEMBERSHIP_DIVERGENCE / WRONG SECONDARY FILTER / DOUBLE OVERLAP FILTER",
    );
    assert!(
        JAVA_COVERAGE_ANNOTATE.contains("likelihoods.evidenceCount()")
            && !JAVA_COVERAGE_ANNOTATE.contains("overlaps")
            && !JAVA_COVERAGE_ANNOTATE.contains("retainEvidence"),
        "Java Coverage.annotate is evidenceCount only"
    );
    assert!(
        JAVA_EVIDENCE_COUNT.contains("List::size") && !JAVA_EVIDENCE_COUNT.contains("overlaps"),
        "Java AlleleLikelihoods.evidenceCount is remaining-list cardinality"
    );

    let ann = include_str!("../src/variant_site_hc_annotations.rs");
    let emit = include_str!("../src/region_vcf_emit.rs");
    let early = include_str!("../src/hc_genotyping_engine/genotype_site_early_template.rs");
    let cov = ann
        .split("pub fn coverage_evidence_count")
        .nth(1)
        .expect("coverage fn")
        .split("fn mq_rms_of_sample_evidence")
        .next()
        .expect("coverage body");
    assert!(
        !cov.contains("java_alignment_read_overlaps_interval")
            && !cov.contains("seen.remove(&idx)"),
        "6R.200 proved the second overlap filter; 6R.201 removed it from Coverage"
    );
    let annotate = ann
        .split("pub fn annotate_hc_variant_site")
        .nth(1)
        .expect("annotate")
        .split("let qd_depth")
        .next()
        .expect("dp assign");
    assert!(
        annotate.contains("coverage_evidence_count(")
            && annotate.contains("position_1based")
            && annotate.contains("config.informative_read_overlap_margin"),
        "INFO DP is assigned from coverage_evidence_count with locus+margin"
    );
    assert!(
        early.contains("fn annotation_likelihoods_from_stored_unique_evidence")
            && early.contains("with_annotation_likelihoods"),
        "6R.199 stored-unique attach stays; this round does not reopen it"
    );
    assert!(
        emit.contains("if call.annotation_likelihoods.is_empty()")
            && emit.contains("likelihoods: &call.annotation_likelihoods"),
        "emit still consumes the attached object when present"
    );
    assert!(
        !ann.contains("92307403") && !early.contains("92307403"),
        "no locus-specific Coverage patch"
    );
}

#[test]
fn forensic_6r200_coverage_no_secondary_overlap_filter() {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing P12 ref/BAM");
        return;
    }
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    kv("target", "2:92307403 C/A");
    kv(
        "classification",
        "ANNOTATION_INPUT_MEMBERSHIP_DIVERGENCE / WRONG SECONDARY FILTER / DOUBLE OVERLAP FILTER",
    );

    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let walk = |interval: &str| {
        let specs = parse_intervals_cli_string(&dict, interval).expect("interval");
        traverse_assembly_region_walker(
            &dict,
            &specs,
            &ref_fasta,
            &bam,
            &ReadFilterParams::gatk_standard_hc(),
            &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
        )
        .expect("walk")
    };
    let gap_regions = flatten_assembly_regions(&walk(GAP_INTERVAL));
    let tg_regions = flatten_assembly_regions(&walk(TG_INTERVAL));
    let mid_regions = flatten_assembly_regions(&walk(MID_INTERVAL));
    let covering_gap = |pos: u64| {
        gap_regions
            .iter()
            .find(|r| {
                matches!(
                    call_disposition(r),
                    AssemblyRegionCallDisposition::ActiveFull
                ) && r.start.get() <= pos
                    && r.end.get() >= pos
            })
            .unwrap_or_else(|| panic!("gap covering {pos}"))
    };
    let covering_tg = |pos: u64| {
        tg_regions
            .iter()
            .find(|r| {
                matches!(
                    call_disposition(r),
                    AssemblyRegionCallDisposition::ActiveFull
                ) && r.start.get() <= pos
                    && r.end.get() >= pos
            })
            .unwrap_or_else(|| panic!("tg covering {pos}"))
    };
    let covering_mid = mid_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_MID
                && r.end.get() >= CLOSED_MID
        })
        .expect("mid covering");
    let covering_gt = covering_gap(CLOSED_GT);
    let covering_ag = covering_gap(CLOSED_AG);
    let covering_ac = covering_gap(CLOSED_AC);
    let covering_indel = covering_tg(CLOSED_INDEL);
    let covering_tg_snp = covering_tg(CLOSED_TG);
    let covering_qual = covering_tg(CLOSED_QUAL);
    let covering_target = covering_tg(TARGET);

    let call_one = |region: &gatk_haplotypecaller::AssemblyRegion| {
        HaplotypeCallerEngine::call_region(
            region,
            &dict,
            &ref_fasta,
            &CallRegionArgs::strict_java(),
        )
        .expect("call")
        .expect("outcome")
    };
    let reuse = |a: &gatk_haplotypecaller::AssemblyRegion,
                 b: &gatk_haplotypecaller::AssemblyRegion,
                 oa: &gatk_haplotypecaller::CallRegionOutcome| {
        if std::ptr::eq(a, b) {
            oa.clone()
        } else {
            call_one(b)
        }
    };
    let gap_outcome = call_one(covering_gt);
    let ag_outcome = reuse(covering_gt, covering_ag, &gap_outcome);
    let ac_outcome = reuse(covering_gt, covering_ac, &gap_outcome);
    let indel_outcome = call_one(covering_indel);
    let tg_snp_outcome = reuse(covering_indel, covering_tg_snp, &indel_outcome);
    let qual_outcome = reuse(covering_indel, covering_qual, &indel_outcome);
    let outcome = reuse(covering_indel, covering_target, &indel_outcome);
    let mid_outcome = call_one(covering_mid);

    let find = |o: &gatk_haplotypecaller::CallRegionOutcome, pos: u64, r: &str, a: &str| {
        o.genotyped_calls
            .iter()
            .find(|c| {
                c.event.start_1based == GenomePosition::new_1based(pos)
                    && c.event.ref_allele == r
                    && c.event.alt_allele == a
            })
            .unwrap_or_else(|| panic!("{pos} {r}/{a}"))
            .clone()
    };
    let gt_snp = find(&gap_outcome, CLOSED_GT, "G", "T");
    let ag = find(&ag_outcome, CLOSED_AG, "A", "G");
    let ac = find(&ac_outcome, CLOSED_AC, "A", "C");
    let indel = find(&indel_outcome, CLOSED_INDEL, "TTC", "T");
    let tg = find(&tg_snp_outcome, CLOSED_TG, "T", "G");
    let qual = find(&qual_outcome, CLOSED_QUAL, "CT", "C");
    let mid = find(&mid_outcome, CLOSED_MID, "G", "A");
    let target = find(&outcome, TARGET, TARGET_REF, TARGET_ALT);

    let attached = unique_indices(&target.annotation_likelihoods);
    kv("attached_n", attached.len().to_string());
    assert_eq!(
        attached.len(),
        6,
        "6R.199 annotation object remains stored unique n=6"
    );
    assert_eq!(
        attached,
        unique_indices(&outcome.read_likelihoods),
        "Coverage receives the attached stored-unique object, not a later subset"
    );
    dump_reads("attached_all", &outcome.genotyping_reads, &attached);

    let mut kept = BTreeSet::new();
    let mut removed = BTreeSet::new();
    for idx in &attached {
        let rec = &outcome.genotyping_reads[*idx];
        if java_alignment_read_overlaps_interval(rec, TARGET, TARGET, MARGIN) {
            kept.insert(*idx);
        } else {
            removed.insert(*idx);
        }
    }
    kv("overlap_kept_n", kept.len().to_string());
    kv("overlap_removed_n", removed.len().to_string());
    dump_reads("coverage_kept", &outcome.genotyping_reads, &kept);
    dump_reads("coverage_removed", &outcome.genotyping_reads, &removed);
    assert_eq!(kept.len(), 4, "overlap predicate retains four rows");
    assert_eq!(removed.len(), 2, "overlap predicate removes two rows");
    assert!(
        has_member(
            &outcome.genotyping_reads,
            &removed,
            EXTRA_QNAME_A,
            EXTRA_FLAG_A
        ) && has_member(
            &outcome.genotyping_reads,
            &removed,
            EXTRA_QNAME_B,
            EXTRA_FLAG_B
        ),
        "removed rows are the two extra stored unique reads"
    );
    assert!(
        has_member(
            &outcome.genotyping_reads,
            &attached,
            EXTRA_QNAME_A,
            EXTRA_FLAG_A
        ) && has_member(
            &outcome.genotyping_reads,
            &attached,
            EXTRA_QNAME_B,
            EXTRA_FLAG_B
        ),
        "those two rows are still on the annotation object Coverage receives"
    );

    let java_equiv_n = attached.len() as i32;
    let rust_cov = coverage_evidence_count(
        &outcome.genotyping_reads,
        &target.annotation_likelihoods,
        TARGET,
        TARGET,
        MARGIN,
    );
    kv("java_equivalent_evidenceCount", java_equiv_n.to_string());
    kv("rust_coverage_evidence_count", rust_cov.to_string());
    assert_eq!(
        java_equiv_n, 6,
        "Java-equivalent evidenceCount is unique attached cardinality"
    );
    assert_eq!(
        rust_cov, 6,
        "6R.201: Coverage consumes attached unique cardinality (6R.200 proved the old 6→4 filter)"
    );
    assert_eq!(
        java_equiv_n, rust_cov,
        "Coverage evidenceCount matches unique attached n=6"
    );

    assert_eq!(unique_indices(&gt_snp.annotation_likelihoods).len(), 3);
    assert_eq!(unique_indices(&ag.annotation_likelihoods).len(), 3);
    assert_eq!(unique_indices(&ac.annotation_likelihoods).len(), 4);
    assert_eq!(unique_indices(&indel.annotation_likelihoods).len(), 1);
    assert_eq!(unique_indices(&tg.annotation_likelihoods).len(), 1);
    assert_eq!(unique_indices(&qual.annotation_likelihoods).len(), 0);
    assert_eq!(unique_indices(&mid.annotation_likelihoods).len(), 3);

    assert_eq!(target.genotype.format.ad_as_i32(), vec![2, 4]);
    assert_eq!(target.genotype.format.dp.as_i32(), 6);
    assert_eq!(target.genotype.format.gq.as_i32(), 72);
    assert_eq!(target.genotype.format.pl_as_i32(), vec![162, 0, 72]);

    let mqs =
        rms_mapping_quality_sample_mapqs(&outcome.genotyping_reads, &target.annotation_likelihoods);
    let rms = rms_mapping_quality_raw(&mqs).map(|t| t.2).expect("rms");
    assert_eq!(mqs.len(), 6);
    assert!((rms - 40.58).abs() < 0.01);

    let emit = |region: &gatk_haplotypecaller::AssemblyRegion,
                o: &gatk_haplotypecaller::CallRegionOutcome| {
        try_emit_call_region_variants(region, o, "NA12878", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("emit")
    };
    let gap_emitted = emit(covering_ag, &ag_outcome);
    let ac_emitted = emit(covering_ac, &ac_outcome);
    let indel_emitted = emit(covering_indel, &indel_outcome);
    let tg_snp_emitted = emit(covering_tg_snp, &tg_snp_outcome);
    let qual_emitted = emit(covering_qual, &qual_outcome);
    let target_emitted = emit(covering_target, &outcome);
    let mid_emitted = emit(covering_mid, &mid_outcome);

    assert_eq!(
        info_i32(&rec_at(&gap_emitted, CLOSED_GT, "G").info, "DP"),
        Some(3)
    );
    assert_eq!(
        info_i32(&rec_at(&gap_emitted, CLOSED_AG, "A").info, "DP"),
        Some(3)
    );
    assert_eq!(
        info_i32(&rec_at(&ac_emitted, CLOSED_AC, "A").info, "DP"),
        Some(4)
    );
    assert_eq!(
        info_i32(&rec_at(&indel_emitted, CLOSED_INDEL, "TTC").info, "DP"),
        Some(1)
    );
    assert_eq!(
        info_i32(&rec_at(&tg_snp_emitted, CLOSED_TG, "T").info, "DP"),
        Some(1)
    );
    let pin = rec_at(&qual_emitted, CLOSED_QUAL, "CT");
    assert_eq!(info_i32(&pin.info, "DP"), Some(2));
    assert!((pin.quality.expect("qual") - 31.60).abs() < 0.005);
    assert!((info_f64(&pin.info, "QD").expect("QD") - 15.80).abs() < 0.005);
    assert_eq!(
        (info_f64(&pin.info, "BaseQRankSum").expect("BaseQ") * 1000.0).round(),
        0.0
    );
    assert_eq!(
        (info_f64(&pin.info, "MQRankSum").expect("MQRS") * 1000.0).round(),
        -674.0
    );
    assert!(info_has(&pin.info, "ReadPosRankSum"));
    assert_eq!(
        info_i32(&rec_at(&mid_emitted, CLOSED_MID, "G").info, "DP"),
        Some(3)
    );

    let rec = rec_at(&target_emitted, TARGET, TARGET_REF);
    assert!((rec.quality.expect("QUAL") - 154.64).abs() < 0.005);
    assert_eq!(
        info_i32(&rec.info, "DP"),
        Some(6),
        "INFO DP is attached unique evidenceCount after 6R.201"
    );
    assert!((info_f64(&rec.info, "MQ").expect("MQ") - 40.58).abs() < 0.005);
    let sample = rec.samples.first().expect("sample");
    assert_eq!(
        sample.gt.as_ref().map(|g| g.to_string()).as_deref(),
        Some("0/1")
    );
    assert_eq!(sample.ad.as_deref(), Some(&[2u32, 4][..]));
    assert_eq!(sample.dp, Some(6));
    assert_eq!(sample.pl.as_deref(), Some(&[162u32, 0, 72][..]));
    kv(
        "first_arrow",
        "6R.200 proved Java Coverage.annotate = evidenceCount(); the second overlap filter 6→4 was removed by 6R.201",
    );
    kv(
        "recommendation_6r201",
        "closed by 6R.201: Coverage consumes attached unique cardinality; 6R.199 object construction unchanged",
    );
}
