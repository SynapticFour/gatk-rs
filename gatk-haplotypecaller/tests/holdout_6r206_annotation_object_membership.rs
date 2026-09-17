//! 6R.206 live: `2:92316416 C/A` one-read hom-alt now attaches retainEvidence n=1
//! (closed by 6R.207). Coverage consumes that object (INFO DP=1).
//! Skipped unless `HOLDOUT_6R206=1`.
//!
//! ```text
//! HOLDOUT_6R206=1 cargo test -p gatk-haplotypecaller --test holdout_6r206_annotation_object_membership -- --nocapture --test-threads=1
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

const GAP_INTERVAL: &str = "2:92305500-92305850";
const TG_INTERVAL: &str = "2:92307200-92307550";
const MID_INTERVAL: &str = "2:92316200-92316580";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_316_416;
const CLOSED_AG: u64 = 92_305_635;
const CLOSED_QUAL: u64 = 92_307_359;
const CLOSED_CA: u64 = 92_307_403;
const CLOSED_HOM_ALT: u64 = 92_316_296;
const CLOSED_MID: u64 = 92_316_347;
const MARGIN: i32 = 2;
const JAVA_QNAME: &str = "H06HDADXX130110:2:1101:10046:78083";
const JAVA_FLAG: u16 = 147;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R206\t{key}\t{}", value.as_ref());
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
fn holdout_6r206_annotation_object_membership() {
    if std::env::var("HOLDOUT_6R206").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R206=1");
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
    kv("production_change", "CLOSED by 6R.207 object attach");
    kv("target", "2:92316416 C/A");
    kv(
        "classification",
        "ANNOTATION_SOURCE_OBJECT_DIVERGENCE / one-read hom-alt now attaches loc-loop retainEvidence n=1",
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
    let covering_ag = gap_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_AG
                && r.end.get() >= CLOSED_AG
        })
        .expect("covering A/G");
    let covering_qual = tg_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_QUAL
                && r.end.get() >= CLOSED_QUAL
        })
        .expect("covering CT/C");
    let covering_ca = tg_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_CA
                && r.end.get() >= CLOSED_CA
        })
        .expect("covering C/A");
    let covering_target = mid_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= TARGET
                && r.end.get() >= TARGET
        })
        .expect("covering C/A target");
    let covering_hom = mid_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_HOM_ALT
                && r.end.get() >= CLOSED_HOM_ALT
        })
        .expect("covering A/T");
    let covering_mid = mid_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_MID
                && r.end.get() >= CLOSED_MID
        })
        .expect("covering G/A");
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
    let ag_outcome = call_one(covering_ag);
    let qual_outcome = call_one(covering_qual);
    let ca_outcome = if std::ptr::eq(covering_qual, covering_ca) {
        qual_outcome.clone()
    } else {
        call_one(covering_ca)
    };
    let outcome = call_one(covering_target);
    let hom_outcome = if std::ptr::eq(covering_target, covering_hom) {
        outcome.clone()
    } else {
        call_one(covering_hom)
    };
    let mid_outcome = if std::ptr::eq(covering_target, covering_mid) {
        outcome.clone()
    } else {
        call_one(covering_mid)
    };
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
    let ag = find(&ag_outcome, CLOSED_AG, "A", "G");
    let qual = find(&qual_outcome, CLOSED_QUAL, "CT", "C");
    let ca = find(&ca_outcome, CLOSED_CA, "C", "A");
    let hom = find(&hom_outcome, CLOSED_HOM_ALT, "A", "T");
    let mid = find(&mid_outcome, CLOSED_MID, "G", "A");
    let target = find(&outcome, TARGET, "C", "A");
    assert_eq!(unique_indices(&ag.annotation_likelihoods).len(), 3);
    assert_eq!(unique_indices(&qual.annotation_likelihoods).len(), 0);
    assert_eq!(unique_indices(&ca.annotation_likelihoods).len(), 6);
    assert_eq!(unique_indices(&hom.annotation_likelihoods).len(), 2);
    assert_eq!(unique_indices(&mid.annotation_likelihoods).len(), 3);
    let attached = unique_indices(&target.annotation_likelihoods);
    assert_eq!(attached.len(), 1);
    assert!(!target.annotation_likelihoods.is_empty());
    let stored = unique_indices(&outcome.read_likelihoods);
    assert_eq!(stored.len(), 3);
    let java_idx = stored.iter().copied().find(|&idx| {
        let rec: &Record = outcome.genotyping_reads[idx].as_ref();
        rec.flags() == JAVA_FLAG && qname(rec) == JAVA_QNAME
    });
    let java_idx = java_idx.expect("Java-counted MAPQ21 row");
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
        3
    );
    assert_eq!(target.genotype.format.ad_as_i32(), vec![0, 1]);
    assert_eq!(target.genotype.format.dp.as_i32(), 1);
    assert_eq!(target.genotype.format.pl_as_i32(), vec![45, 3, 0]);

    let emit = |region: &gatk_haplotypecaller::AssemblyRegion,
                o: &gatk_haplotypecaller::CallRegionOutcome| {
        try_emit_call_region_variants(region, o, "NA12878", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("emit")
    };
    let ag_rec = emit(covering_ag, &ag_outcome)
        .into_iter()
        .find(|r| r.position == CLOSED_AG && r.reference == "A")
        .expect("A/G");
    let pin = emit(covering_qual, &qual_outcome)
        .into_iter()
        .find(|r| r.position == CLOSED_QUAL && r.reference == "CT")
        .expect("CT/C");
    let ca_rec = emit(covering_ca, &ca_outcome)
        .into_iter()
        .find(|r| r.position == CLOSED_CA && r.reference == "C")
        .expect("C/A");
    let hom_rec = emit(covering_hom, &hom_outcome)
        .into_iter()
        .find(|r| r.position == CLOSED_HOM_ALT && r.reference == "A")
        .expect("A/T");
    let target_rec = emit(covering_target, &outcome)
        .into_iter()
        .find(|r| r.position == TARGET && r.reference == "C")
        .expect("C/A");
    assert_eq!(info_i32(&ag_rec.info, "DP"), Some(3));
    assert!((pin.quality.expect("qual") - 31.60).abs() < 0.005);
    assert!(info_has(&pin.info, "BaseQRankSum"));
    assert_eq!(info_i32(&ca_rec.info, "DP"), Some(6));
    assert!((info_f64(&ca_rec.info, "MQ").expect("mq") - 40.58).abs() < 0.005);
    assert!((info_f64(&ca_rec.info, "SOR").expect("sor") - 0.105).abs() < 0.002);
    assert_eq!(info_i32(&hom_rec.info, "DP"), Some(2));
    assert!((info_f64(&hom_rec.info, "MQ").expect("mq") - 47.00).abs() < 0.005);
    assert!((info_f64(&hom_rec.info, "SOR").expect("sor") - 2.303).abs() < 0.002);
    assert_eq!(info_i32(&target_rec.info, "DP"), Some(1));
    assert!((target_rec.quality.expect("QUAL") - 35.48).abs() < 0.005);
    kv("java_info_dp", "1");
    kv("rust_info_dp", "1");
    kv(
        "mq_sor",
        "Java MQ=21.00 SOR=1.609; Rust MQ/SOR follow attached n=1 object",
    );
}
