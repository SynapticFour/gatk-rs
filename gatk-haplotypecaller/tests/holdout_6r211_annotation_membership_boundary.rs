//! 6R.211 live: `2:92325193 C/T` gap-tail het early-template attaches the
//! stored-haplotype loc-loop annotation object (n=3). Coverage consumes
//! that object (INFO DP=3, MQ=28.03, SOR=0.223). FORMAT/QUAL stay
//! Java-equivalent. The non-overlapping same-QNAME mate is not retained.
//! Skipped unless `HOLDOUT_6R211=1`.
//!
//! ```text
//! HOLDOUT_6R211=1 cargo test -p gatk-haplotypecaller --test holdout_6r211_annotation_membership_boundary -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::variant_site_hc_annotations::coverage_evidence_count;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

const TG_INTERVAL: &str = "2:92307200-92307550";
const MID_INTERVAL: &str = "2:92316200-92316580";
const MID_B_INTERVAL: &str = "2:92317000-92319000";
const P12_POST_INTERVAL: &str = "2:92318150-92319220";
const HET_TAIL_INTERVAL: &str = "2:92324900-92325400";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_325_193;
const CLOSED_TG: u64 = 92_307_333;
const CLOSED_CA: u64 = 92_307_403;
const CLOSED_HOM_ALT: u64 = 92_316_296;
const CLOSED_ONE_READ: u64 = 92_316_416;
const CLOSED_MID_B: u64 = 92_317_399;
const CLOSED_POST: u64 = 92_318_199;
const MARGIN: i32 = 2;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R211\t{key}\t{}", value.as_ref());
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

#[test]
fn holdout_6r211_annotation_membership_boundary() {
    if std::env::var("HOLDOUT_6R211").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R211=1");
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
        "one annotation-object construction arrow on is_p12_phase_e_gap_het_event",
    );
    kv("target", "2:92325193 C/T");

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
    let het_regions = flatten_assembly_regions(&walk(HET_TAIL_INTERVAL));
    let covering_tg = tg_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_TG
                && r.end.get() >= CLOSED_TG
        })
        .expect("tg");
    let covering_ca = tg_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_CA
                && r.end.get() >= CLOSED_CA
        })
        .expect("ca");
    let covering_hom = mid_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_HOM_ALT
                && r.end.get() >= CLOSED_HOM_ALT
        })
        .expect("hom");
    let covering_one = mid_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_ONE_READ
                && r.end.get() >= CLOSED_ONE_READ
        })
        .expect("one");
    let covering_mid_b = mid_b_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_MID_B
                && r.end.get() >= CLOSED_MID_B
        })
        .expect("mid-B");
    let covering_post = post_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_POST
                && r.end.get() >= CLOSED_POST
        })
        .expect("post");
    let covering_target = het_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= TARGET
                && r.end.get() >= TARGET
        })
        .expect("het-tail");
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
    let post_outcome = call_one(covering_post);
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
    let post = find(&post_outcome, CLOSED_POST, "C", "T");
    let target = find(&outcome, TARGET, "C", "T");
    assert_eq!(unique_indices(&tg.annotation_likelihoods).len(), 1);
    assert_eq!(unique_indices(&ca.annotation_likelihoods).len(), 6);
    assert_eq!(unique_indices(&hom.annotation_likelihoods).len(), 2);
    assert_eq!(unique_indices(&one.annotation_likelihoods).len(), 1);
    assert_eq!(unique_indices(&mid_b.annotation_likelihoods).len(), 2);
    assert_eq!(unique_indices(&post.annotation_likelihoods).len(), 1);
    let attached = unique_indices(&target.annotation_likelihoods);
    let stored = unique_indices(&outcome.read_likelihoods);
    assert_eq!(stored.len(), 4);
    assert_eq!(attached.len(), 3);
    assert_eq!(
        coverage_evidence_count(
            &outcome.genotyping_reads,
            &target.annotation_likelihoods,
            TARGET,
            TARGET,
            MARGIN,
        ),
        3
    );
    assert_eq!(target.genotype.format.ad_as_i32(), vec![1, 2]);
    assert_eq!(target.genotype.format.pl_as_i32(), vec![81, 0, 36]);

    let emit = |region: &gatk_haplotypecaller::AssemblyRegion,
                o: &gatk_haplotypecaller::CallRegionOutcome| {
        try_emit_call_region_variants(region, o, "NA12878", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("emit")
    };
    let pin = |emitted: Vec<gatk_core::io::vcf::VcfRecord>, pos: u64, r: &str| {
        emitted
            .into_iter()
            .find(|rec| rec.position == pos && rec.reference == r)
            .unwrap_or_else(|| panic!("emit {pos}"))
    };
    let tg_rec = pin(emit(covering_tg, &tg_outcome), CLOSED_TG, "T");
    let ca_rec = pin(emit(covering_ca, &ca_outcome), CLOSED_CA, "C");
    let hom_rec = pin(emit(covering_hom, &hom_outcome), CLOSED_HOM_ALT, "A");
    let one_rec = pin(emit(covering_one, &one_outcome), CLOSED_ONE_READ, "C");
    let mid_b_rec = pin(emit(covering_mid_b, &mid_b_outcome), CLOSED_MID_B, "C");
    let post_rec = pin(emit(covering_post, &post_outcome), CLOSED_POST, "C");
    let target_rec = pin(emit(covering_target, &outcome), TARGET, "C");
    assert_eq!(info_i32(&tg_rec.info, "DP"), Some(1));
    assert!((info_f64(&tg_rec.info, "MQ").expect("mq") - 44.00).abs() < 0.005);
    assert_eq!(info_i32(&ca_rec.info, "DP"), Some(6));
    assert!((info_f64(&ca_rec.info, "MQ").expect("mq") - 40.58).abs() < 0.005);
    assert_eq!(info_i32(&hom_rec.info, "DP"), Some(2));
    assert!((info_f64(&hom_rec.info, "MQ").expect("mq") - 47.00).abs() < 0.005);
    assert_eq!(info_i32(&one_rec.info, "DP"), Some(1));
    assert!((info_f64(&one_rec.info, "MQ").expect("mq") - 21.00).abs() < 0.005);
    assert_eq!(info_i32(&mid_b_rec.info, "DP"), Some(2));
    assert!((info_f64(&mid_b_rec.info, "MQ").expect("mq") - 27.00).abs() < 0.005);
    assert_eq!(info_i32(&post_rec.info, "DP"), Some(1));
    assert!((info_f64(&post_rec.info, "MQ").expect("mq") - 24.00).abs() < 0.005);
    assert_eq!(info_i32(&target_rec.info, "DP"), Some(3));
    assert!((target_rec.quality.expect("QUAL") - 73.64).abs() < 0.005);
    let mq = info_f64(&target_rec.info, "MQ").expect("MQ");
    let sor = info_f64(&target_rec.info, "SOR").expect("SOR");
    assert!((mq - 28.03).abs() < 0.005);
    assert!((sor - 0.223).abs() < 0.002);
    kv("java_info_dp", "3");
    kv("rust_info_dp", "3");
    kv("mq_sor", "Java MQ=28.03 SOR=0.223; loc-loop n=3");
}
