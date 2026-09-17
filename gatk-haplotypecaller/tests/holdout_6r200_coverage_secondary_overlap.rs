//! 6R.200 live: frozen proof that ±2 overlap would drop two of the n=6
//! attached rows at `2:92307403 C/A`. After 6R.201 Coverage / INFO DP = 6.
//! PRODUCTION CHANGE: NONE (6R.200).
//! Skipped unless `HOLDOUT_6R200=1`.
//!
//! ```text
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

const GAP_INTERVAL: &str = "2:92305500-92305850";
const TG_INTERVAL: &str = "2:92307200-92307550";
const MID_INTERVAL: &str = "2:92316200-92316580";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_307_403;
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
    println!("6R200\t{key}\t{}", value.as_ref());
}

fn unique_indices(likelihoods: &[gatk_haplotypecaller::RegionReadLikelihood]) -> BTreeSet<usize> {
    likelihoods.iter().map(|c| c.read_index.get()).collect()
}

fn qname(rec: &Record) -> String {
    String::from_utf8_lossy(rec.qname()).into_owned()
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

fn covering_full(
    regions: &[gatk_haplotypecaller::AssemblyRegion],
    pos: u64,
) -> &gatk_haplotypecaller::AssemblyRegion {
    regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= pos
                && r.end.get() >= pos
        })
        .unwrap_or_else(|| panic!("covering {pos}"))
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

fn find_call(
    outcome: &gatk_haplotypecaller::CallRegionOutcome,
    pos: u64,
    r: &str,
    a: &str,
) -> gatk_haplotypecaller::hc_genotyping_engine::GenotypedSiteCall {
    outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(pos)
                && c.event.ref_allele == r
                && c.event.alt_allele == a
        })
        .unwrap_or_else(|| panic!("{pos} {r}/{a}"))
        .clone()
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

#[test]
fn holdout_6r200_coverage_secondary_overlap() {
    if std::env::var("HOLDOUT_6R200").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R200=1");
        return;
    }
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing P12 ref/BAM");
        return;
    }
    kv(
        "java_pin",
        "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77",
    );
    kv("production_change", "NONE");
    kv("target", "2:92307403 C/A");

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
    let covering_gt = covering_full(&gap_regions, CLOSED_GT);
    let covering_ag = covering_full(&gap_regions, CLOSED_AG);
    let covering_ac = covering_full(&gap_regions, CLOSED_AC);
    let covering_indel = covering_full(&tg_regions, CLOSED_INDEL);
    let covering_tg_snp = covering_full(&tg_regions, CLOSED_TG);
    let covering_qual = covering_full(&tg_regions, CLOSED_QUAL);
    let covering_tg = covering_full(&tg_regions, TARGET);
    let covering_mid = covering_full(&mid_regions, CLOSED_MID);

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
    let gap_outcome = call_one(covering_ag);
    let ac_outcome = reuse(covering_ag, covering_ac, &gap_outcome);
    let indel_outcome = call_one(covering_indel);
    let tg_snp_outcome = reuse(covering_indel, covering_tg_snp, &indel_outcome);
    let qual_outcome = reuse(covering_indel, covering_qual, &indel_outcome);
    let tg_outcome = reuse(covering_indel, covering_tg, &indel_outcome);
    let mid_outcome = call_one(covering_mid);

    let closed_gt = find_call(&gap_outcome, CLOSED_GT, "G", "T");
    let closed_ag = find_call(&gap_outcome, CLOSED_AG, "A", "G");
    let closed_ac = find_call(&ac_outcome, CLOSED_AC, "A", "C");
    let closed_indel = find_call(&indel_outcome, CLOSED_INDEL, "TTC", "T");
    let closed_tg = find_call(&tg_snp_outcome, CLOSED_TG, "T", "G");
    let closed_qual = find_call(&qual_outcome, CLOSED_QUAL, "CT", "C");
    let closed_mid = find_call(&mid_outcome, CLOSED_MID, "G", "A");
    let target = find_call(&tg_outcome, TARGET, "C", "A");

    assert_eq!(unique_indices(&closed_gt.annotation_likelihoods).len(), 3);
    assert_eq!(unique_indices(&closed_ag.annotation_likelihoods).len(), 3);
    assert_eq!(unique_indices(&closed_ac.annotation_likelihoods).len(), 4);
    assert_eq!(
        unique_indices(&closed_indel.annotation_likelihoods).len(),
        1
    );
    assert_eq!(unique_indices(&closed_tg.annotation_likelihoods).len(), 1);
    assert_eq!(unique_indices(&closed_qual.annotation_likelihoods).len(), 0);
    assert_eq!(unique_indices(&closed_mid.annotation_likelihoods).len(), 3);

    let attached = unique_indices(&target.annotation_likelihoods);
    assert_eq!(attached.len(), 6);
    let removed: BTreeSet<usize> = attached
        .iter()
        .copied()
        .filter(|idx| {
            !java_alignment_read_overlaps_interval(
                &tg_outcome.genotyping_reads[*idx],
                TARGET,
                TARGET,
                MARGIN,
            )
        })
        .collect();
    assert_eq!(removed.len(), 2);
    assert!(
        has_member(
            &tg_outcome.genotyping_reads,
            &removed,
            EXTRA_QNAME_A,
            EXTRA_FLAG_A
        ) && has_member(
            &tg_outcome.genotyping_reads,
            &removed,
            EXTRA_QNAME_B,
            EXTRA_FLAG_B
        )
    );
    assert_eq!(
        coverage_evidence_count(
            &tg_outcome.genotyping_reads,
            &target.annotation_likelihoods,
            TARGET,
            TARGET,
            MARGIN,
        ),
        6,
        "6R.201: Coverage consumes the six attached rows directly"
    );
    let mqs = rms_mapping_quality_sample_mapqs(
        &tg_outcome.genotyping_reads,
        &target.annotation_likelihoods,
    );
    let rms = rms_mapping_quality_raw(&mqs).map(|t| t.2).expect("rms");
    assert_eq!(mqs.len(), 6);
    assert!((rms - 40.58).abs() < 0.01);

    assert_eq!(target.genotype.format.ad_as_i32(), vec![2, 4]);
    assert_eq!(target.genotype.format.dp.as_i32(), 6);
    assert_eq!(target.genotype.format.gq.as_i32(), 72);
    assert_eq!(target.genotype.format.pl_as_i32(), vec![162, 0, 72]);

    let gap_emitted = try_emit_call_region_variants(
        covering_ag,
        &gap_outcome,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("gap emit");
    let ac_emitted = try_emit_call_region_variants(
        covering_ac,
        &ac_outcome,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("ac emit");
    let indel_emitted = try_emit_call_region_variants(
        covering_indel,
        &indel_outcome,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("indel emit");
    let tg_snp_emitted = try_emit_call_region_variants(
        covering_tg_snp,
        &tg_snp_outcome,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("tg snp emit");
    let qual_emitted = try_emit_call_region_variants(
        covering_qual,
        &qual_outcome,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("qual emit");
    let tg_emitted = try_emit_call_region_variants(
        covering_tg,
        &tg_outcome,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("tg emit");
    let mid_emitted = try_emit_call_region_variants(
        covering_mid,
        &mid_outcome,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("mid emit");
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
    let target_rec = rec_at(&tg_emitted, TARGET, "C");
    assert!((target_rec.quality.expect("QUAL") - 154.64).abs() < 0.005);
    assert_eq!(info_i32(&target_rec.info, "DP"), Some(6));
    assert!((info_f64(&target_rec.info, "MQ").expect("MQ") - 40.58).abs() < 0.005);
    let sample = target_rec.samples.first().expect("sample");
    assert_eq!(
        sample.gt.as_ref().map(|g| g.to_string()).as_deref(),
        Some("0/1")
    );
    assert_eq!(sample.ad.as_deref(), Some(&[2u32, 4][..]));
    assert_eq!(sample.dp, Some(6));
    assert_eq!(sample.pl.as_deref(), Some(&[162u32, 0, 72][..]));
    kv("attached_n", "6");
    kv("coverage_n", "6");
    kv("info_dp", "6");
}
