//! 6R.210 live: `2:92318199 C/T` gap-sparse shaped-early attaches the
//! stored-haplotype loc-loop annotation object (n=1). Coverage consumes
//! that object (INFO DP=1, MQ=24.00, SOR=1.609). FORMAT/QUAL stay
//! Java-equivalent. MQRankSum is omitted (empty REF fillQuals).
//! Skipped unless `HOLDOUT_6R210=1`.
//!
//! ```text
//! HOLDOUT_6R210=1 cargo test -p gatk-haplotypecaller --test holdout_6r210_annotation_source_boundary -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::hc_genotyping_engine::java_alignment_read_overlaps_interval;
use gatk_haplotypecaller::variant_site_hc_annotations::coverage_evidence_count;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use rust_htslib::bam::Record;
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

const TG_INTERVAL: &str = "2:92307200-92307550";
const MID_INTERVAL: &str = "2:92316200-92316580";
const MID_B_INTERVAL: &str = "2:92317000-92319000";
const P12_POST_INTERVAL: &str = "2:92318150-92319220";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_318_199;
const CLOSED_TG: u64 = 92_307_333;
const CLOSED_CA: u64 = 92_307_403;
const CLOSED_HOM_ALT: u64 = 92_316_296;
const CLOSED_ONE_READ: u64 = 92_316_416;
const CLOSED_MID_B: u64 = 92_317_399;
const MARGIN: i32 = 2;
const JAVA_QNAME: &str = "H06JUADXX130110:1:1101:10002:20630";
const JAVA_FLAG: u16 = 99;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R210\t{key}\t{}", value.as_ref());
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

#[test]
fn holdout_6r210_annotation_source_boundary() {
    if std::env::var("HOLDOUT_6R210").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R210=1");
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
    kv(
        "production_change",
        "one annotation-object construction arrow on gap_sparse_shaped_early",
    );
    kv("target", "2:92318199 C/T");
    kv(
        "classification",
        "ANNOTATION_SOURCE_OBJECT_DIVERGENCE / loc-loop retainEvidence n=1 attached",
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
    let tg_regions = flatten_assembly_regions(&walk(TG_INTERVAL));
    let mid_regions = flatten_assembly_regions(&walk(MID_INTERVAL));
    let mid_b_regions = flatten_assembly_regions(&walk(MID_B_INTERVAL));
    let post_regions = flatten_assembly_regions(&walk(P12_POST_INTERVAL));
    let covering_tg = tg_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_TG
                && r.end.get() >= CLOSED_TG
        })
        .expect("tg covering T/G");
    let covering_ca = tg_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_CA
                && r.end.get() >= CLOSED_CA
        })
        .expect("tg covering C/A");
    let covering_hom = mid_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_HOM_ALT
                && r.end.get() >= CLOSED_HOM_ALT
        })
        .expect("mid covering A/T");
    let covering_one = mid_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_ONE_READ
                && r.end.get() >= CLOSED_ONE_READ
        })
        .expect("mid covering C/A");
    let covering_mid_b = mid_b_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_MID_B
                && r.end.get() >= CLOSED_MID_B
        })
        .expect("mid-B covering C/A");
    let covering_target = post_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= TARGET
                && r.end.get() >= TARGET
        })
        .expect("p12_post covering 92318199");
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
    let tg_outcome = call_one(covering_tg);
    let ca_outcome = if std::ptr::eq(covering_tg, covering_ca) {
        tg_outcome.clone()
    } else {
        call_one(covering_ca)
    };
    let hom_outcome = call_one(covering_hom);
    let one_outcome = if std::ptr::eq(covering_hom, covering_one) {
        hom_outcome.clone()
    } else {
        call_one(covering_one)
    };
    let mid_b_outcome = call_one(covering_mid_b);
    let outcome = call_one(covering_target);
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
    let tg = find(&tg_outcome, CLOSED_TG, "T", "G");
    let ca = find(&ca_outcome, CLOSED_CA, "C", "A");
    let hom = find(&hom_outcome, CLOSED_HOM_ALT, "A", "T");
    let one = find(&one_outcome, CLOSED_ONE_READ, "C", "A");
    let mid_b = find(&mid_b_outcome, CLOSED_MID_B, "C", "A");
    let target = find(&outcome, TARGET, "C", "T");
    assert_eq!(unique_indices(&tg.annotation_likelihoods).len(), 1);
    assert_eq!(unique_indices(&ca.annotation_likelihoods).len(), 6);
    assert_eq!(unique_indices(&hom.annotation_likelihoods).len(), 2);
    assert_eq!(unique_indices(&one.annotation_likelihoods).len(), 1);
    assert_eq!(unique_indices(&mid_b.annotation_likelihoods).len(), 2);
    let attached = unique_indices(&target.annotation_likelihoods);
    let stored = unique_indices(&outcome.read_likelihoods);
    assert_eq!(stored.len(), 4);
    assert_eq!(attached.len(), 1);
    assert!(!target.annotation_likelihoods.is_empty());
    let java_idx = stored
        .iter()
        .copied()
        .find(|&idx| {
            let rec: &Record = outcome.genotyping_reads[idx].as_ref();
            rec.flags() == JAVA_FLAG && qname(rec) == JAVA_QNAME
        })
        .expect("Java-counted MAPQ24 row");
    assert!(attached.contains(&java_idx));
    assert!(java_alignment_read_overlaps_interval(
        &outcome.genotyping_reads[java_idx],
        TARGET,
        TARGET,
        MARGIN,
    ));
    assert_eq!(
        coverage_evidence_count(
            &outcome.genotyping_reads,
            &target.annotation_likelihoods,
            TARGET,
            TARGET,
            MARGIN,
        ),
        1
    );
    assert_eq!(
        coverage_evidence_count(
            &outcome.genotyping_reads,
            &outcome.read_likelihoods,
            TARGET,
            TARGET,
            MARGIN,
        ),
        4
    );
    assert_eq!(target.genotype.format.ad_as_i32(), vec![0, 1]);
    assert_eq!(target.genotype.format.dp.as_i32(), 1);
    assert_eq!(target.genotype.format.pl_as_i32(), vec![45, 3, 0]);

    let emit = |region: &gatk_haplotypecaller::AssemblyRegion,
                o: &gatk_haplotypecaller::CallRegionOutcome| {
        try_emit_call_region_variants(region, o, "NA12878", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("emit")
    };
    let tg_rec = emit(covering_tg, &tg_outcome)
        .into_iter()
        .find(|r| r.position == CLOSED_TG && r.reference == "T")
        .expect("T/G");
    let ca_rec = emit(covering_ca, &ca_outcome)
        .into_iter()
        .find(|r| r.position == CLOSED_CA && r.reference == "C")
        .expect("C/A");
    let hom_rec = emit(covering_hom, &hom_outcome)
        .into_iter()
        .find(|r| r.position == CLOSED_HOM_ALT && r.reference == "A")
        .expect("A/T");
    let one_rec = emit(covering_one, &one_outcome)
        .into_iter()
        .find(|r| r.position == CLOSED_ONE_READ && r.reference == "C")
        .expect("C/A");
    let mid_b_rec = emit(covering_mid_b, &mid_b_outcome)
        .into_iter()
        .find(|r| r.position == CLOSED_MID_B && r.reference == "C")
        .expect("C/A");
    let target_rec = emit(covering_target, &outcome)
        .into_iter()
        .find(|r| r.position == TARGET && r.reference == "C")
        .expect("C/T");
    assert_eq!(info_i32(&tg_rec.info, "DP"), Some(1));
    assert!((info_f64(&tg_rec.info, "MQ").expect("mq") - 44.00).abs() < 0.005);
    assert!((info_f64(&tg_rec.info, "SOR").expect("sor") - 1.609).abs() < 0.002);
    assert_eq!(info_i32(&ca_rec.info, "DP"), Some(6));
    assert!((info_f64(&ca_rec.info, "MQ").expect("mq") - 40.58).abs() < 0.005);
    assert!((info_f64(&ca_rec.info, "SOR").expect("sor") - 0.105).abs() < 0.002);
    assert_eq!(info_i32(&hom_rec.info, "DP"), Some(2));
    assert!((info_f64(&hom_rec.info, "MQ").expect("mq") - 47.00).abs() < 0.005);
    assert!((info_f64(&hom_rec.info, "SOR").expect("sor") - 2.303).abs() < 0.002);
    assert_eq!(info_i32(&one_rec.info, "DP"), Some(1));
    assert!((info_f64(&one_rec.info, "MQ").expect("mq") - 21.00).abs() < 0.005);
    assert!((info_f64(&one_rec.info, "SOR").expect("sor") - 1.609).abs() < 0.002);
    assert_eq!(info_i32(&mid_b_rec.info, "DP"), Some(2));
    assert!((info_f64(&mid_b_rec.info, "MQ").expect("mq") - 27.00).abs() < 0.005);
    assert!((info_f64(&mid_b_rec.info, "SOR").expect("sor") - 0.693).abs() < 0.002);
    assert_eq!(info_i32(&target_rec.info, "DP"), Some(1));
    assert!((target_rec.quality.expect("QUAL") - 35.48).abs() < 0.005);
    let mq = info_f64(&target_rec.info, "MQ").expect("MQ");
    let sor = info_f64(&target_rec.info, "SOR").expect("SOR");
    assert!((mq - 24.00).abs() < 0.005);
    assert!((sor - 1.609).abs() < 0.002);
    assert!(!info_has(&target_rec.info, "MQRankSum"));
    kv("java_info_dp", "1");
    kv("rust_info_dp", "1");
    kv("mq_sor", "Java MQ=24.00 SOR=1.609; loc-loop n=1");
}
