//! 6R.202: proof-only. First remaining INFO at `2:92307403 C/A` is Java-only
//! BaseQRankSum / ReadPosRankSum. PRODUCTION CHANGE: NONE.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! 6R.199 object (n=6) and 6R.201 Coverage cardinality stay closed.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r202_ranksum_annotation_boundary -- --nocapture --test-threads=1
//! HOLDOUT_6R202=1 cargo test -p gatk-haplotypecaller --test holdout_6r202_ranksum_boundary -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::fragment_overlap::read_base_quality_at_ref_coord_1based;
use gatk_haplotypecaller::hc_genotyping_engine::java_alignment_read_overlaps_interval;
use gatk_haplotypecaller::variant_site_hc_annotations::{
    baseq_rank_sum_quals_from_likelihoods, coverage_evidence_count,
    mq_rank_sum_quals_from_likelihoods, read_pos_rank_sum_quals_from_likelihoods,
    rms_mapping_quality_raw, rms_mapping_quality_sample_mapqs, HcStrandBiasLikelihoods,
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
const JAVA_RANKSUM_ANNOTATE: &str = r#"
if( likelihoods != null) {
    fillQualsFromLikelihood(vc, likelihoods, refQuals, altQuals);
}
if ( refQuals.isEmpty() && altQuals.isEmpty() ) {
    return Collections.emptyMap();
}
final double zScore = mannWhitneyU.test(...).getZ();
if (Double.isNaN(zScore)) return Collections.emptyMap();
return singletonMap(key, String.format("%.3f", zScore));
"#;
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
    eprintln!("6R202\t{key}\t{}", value.as_ref());
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

#[test]
fn forensic_6r202_ranksum_boundary_contract() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    kv(
        "classification",
        "ANNOTATION_INPUT_MEMBERSHIP_DIVERGENCE / fillQuals getElementForRead empty REF for BaseQ and ReadPos",
    );
    assert!(
        JAVA_RANKSUM_ANNOTATE.contains("fillQualsFromLikelihood")
            && JAVA_RANKSUM_ANNOTATE.contains("isNaN"),
        "Java RankSumTest.annotate is fillQuals then Mann-Whitney; NaN/empty → no key"
    );

    let ann = include_str!("../src/variant_site_hc_annotations.rs");
    let emit = include_str!("../src/region_vcf_emit.rs");
    let early = include_str!("../src/hc_genotyping_engine/genotype_site_early_template.rs");
    let annotate = ann
        .split("pub fn annotate_hc_variant_site")
        .nth(1)
        .expect("annotate");
    assert!(
        annotate.contains("baseq_rank_sum_quals_from_likelihoods")
            && annotate.contains("read_pos_rank_sum_quals_from_likelihoods")
            && annotate.contains("mq_rank_sum_quals_from_likelihoods")
            && annotate.contains("rank_sum_baseq::base_quality_rank_sum")
            && annotate.contains("read_pos_rank_sum::read_pos_rank_sum"),
        "RankSums are invoked inside annotate_hc_variant_site, not skipped on this branch"
    );
    assert!(
        early.contains("fn annotation_likelihoods_from_stored_unique_evidence"),
        "6R.199 stored-unique attach stays"
    );
    let cov = ann
        .split("pub fn coverage_evidence_count")
        .nth(1)
        .expect("coverage fn")
        .split("fn mq_rms_of_sample_evidence")
        .next()
        .expect("coverage body");
    assert!(
        !cov.contains("java_alignment_read_overlaps_interval"),
        "6R.201 Coverage cardinality stays"
    );
    assert!(
        emit.contains("if let Some(z) = ann.base_q_rank_sum")
            && emit.contains("if let Some(z) = ann.read_pos_rank_sum")
            && emit.contains("if let Some(z) = ann.mq_rank_sum"),
        "hc_info_values still inserts finite RankSum z including 0.0"
    );
    assert!(
        !ann.contains("92307403") && !early.contains("92307403") && !emit.contains("92307403"),
        "no locus-specific RankSum patch"
    );
}

#[test]
fn forensic_6r202_ranksum_annotation_boundary() {
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

    let java_vcf = root.join("parity/reports/6r43/p12_indel_mix/java.vcf");
    if java_vcf.is_file() {
        let text = std::fs::read_to_string(&java_vcf).expect("java.vcf");
        let line = text
            .lines()
            .find(|l| l.contains("\t92307403\t") && l.contains("\tC\tA\t"))
            .expect("frozen Java 4.4.0.0 p12_indel_mix 2:92307403 C/A");
        let info = line.split('\t').nth(7).expect("INFO");
        kv("java_vcf_info", info);
        assert!(
            info.contains("BaseQRankSum=-1.834"),
            "pinned Java execution emits BaseQRankSum=-1.834, got {info}"
        );
        assert!(
            info.contains("ReadPosRankSum=1.282"),
            "pinned Java execution emits ReadPosRankSum=1.282, got {info}"
        );
        assert!(
            info.contains("MQRankSum=1.834") && info.contains("DP=6") && info.contains("MQ=40.58"),
            "pinned Java MQ/DP stay, got {info}"
        );
    }

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
    kv("annotation_evidence_n", attached.len().to_string());
    assert_eq!(attached.len(), 6, "6R.199 object remains n=6");
    assert_eq!(
        attached,
        unique_indices(&outcome.read_likelihoods),
        "attached object is stored unique haplotype evidence"
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
        )
    );
    assert_eq!(
        coverage_evidence_count(
            &outcome.genotyping_reads,
            &target.annotation_likelihoods,
            TARGET,
            TARGET,
            MARGIN,
        ),
        6,
        "6R.201 Coverage evidenceCount stays 6"
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

    let (full_ref, full_pad) = outcome.assembly.event_map_reference();
    let apply_pad = outcome
        .assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .and_then(|h| h.genome_loc.as_ref().map(|g| g.start.get()))
        .unwrap_or(full_pad);
    let evidence = HcStrandBiasLikelihoods {
        reads: &outcome.genotyping_reads,
        likelihoods: &target.annotation_likelihoods,
        haplotypes: &outcome.assembly.haplotypes,
        contig: covering_target.contig.as_str(),
        ref_bytes: outcome.assembly.reference_bases(),
        pad_start_1based: apply_pad,
        full_ref_bytes: full_ref,
        full_pad_1based: full_pad,
        max_mnp_distance: outcome.assembly.max_mnp_distance(),
        emit_spanning_dels: true,
    };

    for idx in &attached {
        let rec: &Record = outcome.genotyping_reads[*idx].as_ref();
        let overlap = java_alignment_read_overlaps_interval(rec, TARGET, TARGET, MARGIN);
        let bq = read_base_quality_at_ref_coord_1based(rec, TARGET as i32);
        kv(
            "attached_row",
            format!(
                "idx={idx} qname={} flag={} mapq={} start={} end={} cigar={} overlap_pm2={} baseq_at_locus={:?}",
                qname(rec),
                rec.flags(),
                rec.mapq(),
                rec.pos() + 1,
                gatk_haplotypecaller::read_unclip::alignment_end_1based(rec),
                rec.cigar(),
                overlap,
                bq
            ),
        );
        let extra = (qname(rec) == EXTRA_QNAME_A && rec.flags() == EXTRA_FLAG_A)
            || (qname(rec) == EXTRA_QNAME_B && rec.flags() == EXTRA_FLAG_B);
        if extra {
            assert!(
                !overlap,
                "extra stored unique REF rows miss ±2 overlap at 92307403 on the attached CIGAR"
            );
            assert_eq!(
                bq, None,
                "getReadBaseQualityAtReferenceCoordinate is empty on extra REF rows"
            );
        }
    }

    let (ref_bq, alt_bq) =
        baseq_rank_sum_quals_from_likelihoods(&evidence, TARGET, TARGET_REF, TARGET_ALT);
    let (ref_rp, alt_rp) =
        read_pos_rank_sum_quals_from_likelihoods(&evidence, TARGET, TARGET_REF, TARGET_ALT);
    let (ref_mq, alt_mq) =
        mq_rank_sum_quals_from_likelihoods(&evidence, TARGET, TARGET_REF, TARGET_ALT);
    kv("baseq_lists", format!("REF={:?} ALT={:?}", ref_bq, alt_bq));
    kv(
        "readpos_lists",
        format!("REF={:?} ALT={:?}", ref_rp, alt_rp),
    );
    kv("mqrs_lists", format!("REF={:?} ALT={:?}", ref_mq, alt_mq));
    kv(
        "baseq_n",
        format!("REF={} ALT={}", ref_bq.len(), alt_bq.len()),
    );
    kv(
        "readpos_n",
        format!("REF={} ALT={}", ref_rp.len(), alt_rp.len()),
    );
    kv(
        "mqrs_n",
        format!("REF={} ALT={}", ref_mq.len(), alt_mq.len()),
    );
    assert_eq!(
        ref_bq,
        vec![30.0, 30.0],
        "6R.203: extra REF rows yield original-alignment BaseQ=30"
    );
    assert_eq!(alt_bq, vec![20.0, 20.0, 20.0, 20.0]);
    assert_eq!(
        ref_rp,
        vec![65.0, 91.0],
        "6R.203: extra REF rows yield original-alignment ReadPos"
    );
    assert_eq!(alt_rp, vec![118.0, 87.0, 106.0, 124.0]);
    let mut mq_ref = ref_mq.clone();
    mq_ref.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut mq_alt = alt_mq.clone();
    mq_alt.sort_by(|a, b| a.partial_cmp(b).unwrap());
    assert_eq!(
        mq_ref,
        vec![22.0, 40.0],
        "MQ REF is MAPQ of the two extra stored rows"
    );
    assert_eq!(mq_alt, vec![41.0, 41.0, 44.0, 50.0]);
    kv(
        "shared_empty_ref",
        "attached CIGAR still misses the locus; 6R.203 retries pre-realign CIGAR so REF lists are non-empty",
    );

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
    assert_eq!(info_i32(&rec.info, "DP"), Some(6));
    let mqs =
        rms_mapping_quality_sample_mapqs(&outcome.genotyping_reads, &target.annotation_likelihoods);
    let rms = rms_mapping_quality_raw(&mqs).map(|t| t.2).expect("rms");
    assert!((rms - 40.58).abs() < 0.01);
    assert!((info_f64(&rec.info, "MQ").expect("MQ") - 40.58).abs() < 0.005);
    let mqrs = info_f64(&rec.info, "MQRankSum").expect("MQRankSum still emits");
    assert_eq!((mqrs * 1000.0).round(), 1834.0);
    let bq = info_f64(&rec.info, "BaseQRankSum").expect("6R.203 emits BaseQRankSum");
    assert_eq!((bq * 1000.0).round(), -1834.0);
    let rprs = info_f64(&rec.info, "ReadPosRankSum").expect("6R.203 emits ReadPosRankSum");
    assert_eq!((rprs * 1000.0).round(), 1282.0);
    let sample = rec.samples.first().expect("sample");
    assert_eq!(
        sample.gt.as_ref().map(|g| g.to_string()).as_deref(),
        Some("0/1")
    );
    assert_eq!(sample.ad.as_deref(), Some(&[2u32, 4][..]));
    assert_eq!(sample.dp, Some(6));
    assert_eq!(sample.gq, Some(72.0));
    assert_eq!(sample.pl.as_deref(), Some(&[162u32, 0, 72][..]));
    kv(
        "first_arrow",
        "invocation+n=6 object reach RankSum; attached CIGAR still misses the locus; 6R.203 pre-realign getElementForRead fills REF and emits Java z",
    );
    kv("java_baseq", "-1.834");
    kv("java_readpos", "1.282");
    kv("java_mqrs", "1.834");
}
