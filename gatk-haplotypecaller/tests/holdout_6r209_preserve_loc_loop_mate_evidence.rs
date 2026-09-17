//! 6R.209 live: `2:92317399 C/A` SiteScore loc-loop retainEvidence n=2
//! (same-QNAME mates FLAG 99/147) is preserved. Coverage consumes n=2
//! (INFO DP=2, SOR=0.693). FORMAT/QUAL unchanged. 6R.180 collapse stays
//! off this object. Skipped unless `HOLDOUT_6R209=1`.
//!
//! ```text
//! HOLDOUT_6R209=1 cargo test -p gatk-haplotypecaller --test holdout_6r209_preserve_loc_loop_mate_evidence -- --nocapture --test-threads=1
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
const MID_B_INTERVAL: &str = "2:92317000-92319000";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_317_399;
const CLOSED_AG: u64 = 92_305_635;
const CLOSED_QUAL: u64 = 92_307_359;
const CLOSED_CA: u64 = 92_307_403;
const CLOSED_HOM_ALT: u64 = 92_316_296;
const CLOSED_ONE_READ: u64 = 92_316_416;
const CLOSED_MID: u64 = 92_316_347;
const MARGIN: i32 = 2;
const JAVA_QNAME: &str = "H06HDADXX130110:1:1101:10040:45938";
const JAVA_FWD_FLAG: u16 = 99;
const JAVA_REV_FLAG: u16 = 147;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R209\t{key}\t{}", value.as_ref());
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
fn holdout_6r209_preserve_loc_loop_mate_evidence() {
    if std::env::var("HOLDOUT_6R209").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R209=1");
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
        "one lifecycle-scoped removal of same-QNAME collapse from the Java loc-loop annotation object",
    );
    kv("target", "2:92317399 C/A");
    kv(
        "classification",
        "ANNOTATION_SOURCE_OBJECT_DIVERGENCE / loc-loop retainEvidence n=2 preserved",
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
    let mid_b_regions = flatten_assembly_regions(&walk(MID_B_INTERVAL));
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
    let covering_one = mid_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_ONE_READ
                && r.end.get() >= CLOSED_ONE_READ
        })
        .expect("covering one-read C/A");
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
    let covering_target = mid_b_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= TARGET
                && r.end.get() >= TARGET
        })
        .expect("covering mid-B C/A");

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
    let hom_outcome = call_one(covering_hom);
    let one_outcome = if std::ptr::eq(covering_hom, covering_one) {
        hom_outcome.clone()
    } else {
        call_one(covering_one)
    };
    let mid_outcome = if std::ptr::eq(covering_hom, covering_mid) {
        hom_outcome.clone()
    } else {
        call_one(covering_mid)
    };
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
    let ag = find(&ag_outcome, CLOSED_AG, "A", "G");
    let qual = find(&qual_outcome, CLOSED_QUAL, "CT", "C");
    let ca = find(&ca_outcome, CLOSED_CA, "C", "A");
    let hom = find(&hom_outcome, CLOSED_HOM_ALT, "A", "T");
    let one = find(&one_outcome, CLOSED_ONE_READ, "C", "A");
    let mid = find(&mid_outcome, CLOSED_MID, "G", "A");
    let target = find(&outcome, TARGET, "C", "A");
    assert_eq!(unique_indices(&ag.annotation_likelihoods).len(), 3);
    assert_eq!(unique_indices(&qual.annotation_likelihoods).len(), 0);
    assert_eq!(unique_indices(&ca.annotation_likelihoods).len(), 6);
    assert_eq!(unique_indices(&hom.annotation_likelihoods).len(), 2);
    assert_eq!(unique_indices(&one.annotation_likelihoods).len(), 1);
    assert_eq!(unique_indices(&mid.annotation_likelihoods).len(), 3);
    let attached = unique_indices(&target.annotation_likelihoods);
    let stored = unique_indices(&outcome.read_likelihoods);
    assert_eq!(stored.len(), 2);
    assert_eq!(attached.len(), 2);
    assert!(!target.annotation_likelihoods.is_empty());
    let mut retain = BTreeSet::new();
    for idx in &stored {
        if java_alignment_read_overlaps_interval(
            &outcome.genotyping_reads[*idx],
            TARGET,
            TARGET,
            MARGIN,
        ) {
            retain.insert(*idx);
        }
    }
    assert_eq!(retain.len(), 2);
    assert_eq!(attached, retain);
    let rec_at = |idx: usize| {
        let rec: &Record = outcome.genotyping_reads[idx].as_ref();
        (qname(rec), rec.flags())
    };
    assert!(retain
        .iter()
        .any(|&i| rec_at(i) == (JAVA_QNAME.to_string(), JAVA_FWD_FLAG)));
    assert!(retain
        .iter()
        .any(|&i| rec_at(i) == (JAVA_QNAME.to_string(), JAVA_REV_FLAG)));
    assert!(attached
        .iter()
        .any(|&i| rec_at(i) == (JAVA_QNAME.to_string(), JAVA_FWD_FLAG)));
    assert!(attached
        .iter()
        .any(|&i| rec_at(i) == (JAVA_QNAME.to_string(), JAVA_REV_FLAG)));
    assert_eq!(
        coverage_evidence_count(
            &outcome.genotyping_reads,
            &target.annotation_likelihoods,
            TARGET,
            TARGET,
            MARGIN,
        ),
        2
    );
    assert_eq!(
        coverage_evidence_count(
            &outcome.genotyping_reads,
            &outcome.read_likelihoods,
            TARGET,
            TARGET,
            MARGIN,
        ),
        2
    );
    assert_eq!(target.genotype.format.ad_as_i32(), vec![0, 2]);
    assert_eq!(target.genotype.format.dp.as_i32(), 2);
    assert_eq!(target.genotype.format.pl_as_i32(), vec![90, 6, 0]);

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
    let one_rec = emit(covering_one, &one_outcome)
        .into_iter()
        .find(|r| r.position == CLOSED_ONE_READ && r.reference == "C")
        .expect("one-read C/A");
    let target_rec = emit(covering_target, &outcome)
        .into_iter()
        .find(|r| r.position == TARGET && r.reference == "C")
        .expect("mid-B C/A");
    assert_eq!(info_i32(&ag_rec.info, "DP"), Some(3));
    assert!((pin.quality.expect("qual") - 31.60).abs() < 0.005);
    assert!(info_has(&pin.info, "BaseQRankSum"));
    assert_eq!(info_i32(&ca_rec.info, "DP"), Some(6));
    assert!((info_f64(&ca_rec.info, "MQ").expect("mq") - 40.58).abs() < 0.005);
    assert!((info_f64(&ca_rec.info, "SOR").expect("sor") - 0.105).abs() < 0.002);
    assert_eq!(info_i32(&hom_rec.info, "DP"), Some(2));
    assert!((info_f64(&hom_rec.info, "MQ").expect("mq") - 47.00).abs() < 0.005);
    assert!((info_f64(&hom_rec.info, "SOR").expect("sor") - 2.303).abs() < 0.002);
    assert_eq!(info_i32(&one_rec.info, "DP"), Some(1));
    assert!((info_f64(&one_rec.info, "MQ").expect("mq") - 21.00).abs() < 0.005);
    assert!((info_f64(&one_rec.info, "SOR").expect("sor") - 1.609).abs() < 0.002);
    assert_eq!(info_i32(&target_rec.info, "DP"), Some(2));
    assert!((info_f64(&target_rec.info, "MQ").expect("mq") - 27.00).abs() < 0.005);
    assert!((info_f64(&target_rec.info, "SOR").expect("sor") - 0.693).abs() < 0.002);
    assert!((target_rec.quality.expect("QUAL") - 78.32).abs() < 0.005);
    kv("java_info_dp", "2");
    kv("rust_info_dp", "2");
    kv("java_sor", "0.693");
    kv("rust_sor", "0.693");
    kv(
        "first_arrow",
        "loc-loop retainEvidence n=2 same-QNAME mates preserved",
    );
}
