//! 6R.69 forensic: walk `modifyReadQualities` on the canonical first genotyping read.
//!
//! Skipped unless `HOLDOUT_6R69=1`. Coordinate-free proof is
//! `forensic_6r69_quality_contract`. This harness records integer Phred states.
//!
//! ```text
//! HOLDOUT_6R69=1 cargo test -p gatk-haplotypecaller --test holdout_6r69_modify_read_qualities -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::pairhmm_qual::{set_to_fixed_value_if_too_low, MIN_USABLE_Q_SCORE};
use gatk_haplotypecaller::pcr_error_model::{
    apply_pcr_error_model, tandem_repeat_units, PcrErrorModel,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, indel_gop_from_optional_tag,
    prepare_read_quals_for_pairhmm_inplace, take_colocated_merge_numerics,
    traverse_assembly_region_walker, try_emit_call_region_variants, AssemblyRegionCallDisposition,
    CallRegionArgs, HaplotypeCallerEngine, HcLikelihoodEngineConfig, ReadFilterParams,
    WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE, GATK_PARITY_DEFAULT_GCP,
    GATK_PARITY_DEFAULT_INS_QUAL,
};
use rust_htslib::bam::record::Aux;
use serde_json::json;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const POS_SNP: u64 = 29_456_344;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn aux_phred_z(rec: &rust_htslib::bam::Record, tag: &[u8]) -> Option<Vec<u8>> {
    match rec.aux(tag) {
        Ok(Aux::String(s)) => Some(s.bytes().map(|b| b.saturating_sub(33)).collect()),
        _ => None,
    }
}

fn java_fast_round(d: f64) -> i32 {
    if d > 0.0 {
        (d + 0.5) as i32
    } else {
        (d - 0.5) as i32
    }
}

fn java_44_adjusted_qual(repeat_length: usize, rate_factor: f64) -> u8 {
    let q = 40.0 - (repeat_length as f64 / (rate_factor * std::f64::consts::PI)).exp() + 1.0;
    java_fast_round(q).max(10) as u8
}

fn apply_java_44_conservative_pcr(read_bases: &[u8], ins: &mut [u8], del: &mut [u8]) {
    const MAX_REPEAT: usize = 20;
    let mut cache = [0u8; MAX_REPEAT + 1];
    for (i, slot) in cache.iter_mut().enumerate() {
        *slot = java_44_adjusted_qual(i, 3.0);
    }
    for i in 1..read_bases.len() {
        let repeat = tandem_repeat_units(read_bases, i - 1).min(MAX_REPEAT);
        let cap = cache[repeat];
        let idx = i - 1;
        ins[idx] = ins[idx].min(cap);
        del[idx] = del[idx].min(cap);
    }
}

#[test]
fn holdout_6r69_first_genotyping_read_quality_pipeline() {
    if std::env::var("HOLDOUT_6R69").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R69=1");
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
    let region = covering[0];
    let outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");

    let emitted =
        try_emit_call_region_variants(region, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    let vcf = emitted
        .iter()
        .find(|r| {
            r.position == POS_SNP && r.reference == "T" && r.alternate.iter().any(|a| a == "C")
        })
        .expect("lifecycle T/C");
    let pl = vcf.samples[0].pl.clone().expect("PL");
    assert_eq!(pl, vec![542, 0, 1353]);
    let live = take_colocated_merge_numerics();
    let numerics = live.iter().find(|n| n.loc == POS_SNP).expect("merge");
    assert_eq!(numerics.pool_sizes, vec![35, 6, 21, 6]);

    let rec = outcome
        .genotyping_reads
        .iter()
        .find(|r| aux_phred_z(r, b"BI").is_some() && aux_phred_z(r, b"BD").is_some())
        .expect("BI/BD genotyping read");
    let bases = rec.seq().as_bytes();
    let raw_bq = rec.qual().to_vec();
    let mapq = rec.mapq();
    let bi = aux_phred_z(rec, b"BI").expect("BI");
    let bd = aux_phred_z(rec, b"BD").expect("BD");
    assert_eq!(bi.len(), bases.len());
    assert_eq!(bd.len(), bases.len());

    let mut rust_bq = raw_bq.clone();
    prepare_read_quals_for_pairhmm_inplace(
        &mut rust_bq,
        mapq,
        &HcLikelihoodEngineConfig::default(),
    );

    let mut rust_ins = indel_gop_from_optional_tag(Some(&bi), bases.len()).unwrap();
    let mut rust_del = indel_gop_from_optional_tag(Some(&bd), bases.len()).unwrap();
    apply_pcr_error_model(
        &bases,
        &mut rust_ins,
        &mut rust_del,
        PcrErrorModel::Conservative,
    );

    let mut rust_ins_from_bi = bi.clone();
    let mut rust_del_from_bd = bd.clone();
    apply_pcr_error_model(
        &bases,
        &mut rust_ins_from_bi,
        &mut rust_del_from_bd,
        PcrErrorModel::Conservative,
    );

    let mut java_ins = bi.clone();
    let mut java_del = bd.clone();
    apply_java_44_conservative_pcr(&bases, &mut java_ins, &mut java_del);
    let java_ins_pre_floor = java_ins[0];
    for q in java_ins.iter_mut().chain(java_del.iter_mut()) {
        *q = set_to_fixed_value_if_too_low(*q, MIN_USABLE_Q_SCORE, MIN_USABLE_Q_SCORE);
    }

    let repeat0 = tandem_repeat_units(&bases, 0);
    let java_cap = java_44_adjusted_qual(repeat0.min(20), 3.0);
    let rust_cap = gatk_haplotypecaller::pcr_error_model::error_model_adjusted_qual(
        repeat0.min(20),
        PcrErrorModel::Conservative.rate_factor().unwrap(),
    );

    let doc = json!({
        "read_index": 0,
        "read_len": bases.len(),
        "base_0": (bases[0] as char).to_string(),
        "raw": {
            "bq0": raw_bq[0],
            "mapq": mapq,
            "bi0": bi[0],
            "bd0": bd[0],
        },
        "after_bq_cap": { "bq0": rust_bq[0] },
        "indel_source": {
            "q45_fallback": GATK_PARITY_DEFAULT_INS_QUAL,
            "java_bi0": bi[0],
            "java_bd0": bd[0],
        },
        "pcr": {
            "repeat_len_at_0": repeat0,
            "java_cache_slot": java_cap,
            "rust_cache_slot": rust_cap,
            "rust_production_ins0": rust_ins[0],
            "rust_after_bi_sub_ins0": rust_ins_from_bi[0],
            "java_after_pcr_ins0": java_ins_pre_floor,
            "java_after_iq_floor_ins0": java_ins[0],
        },
        "gcp0": GATK_PARITY_DEFAULT_GCP,
        "first_divergence_after_bi_bd": "6R.73: production GOP source is BI/BD; canonical GOP 40 is still PCR-capped",
    });
    println!("{}", serde_json::to_string_pretty(&doc).unwrap());

    assert_eq!(
        GATK_PARITY_DEFAULT_INS_QUAL, 45,
        "Q45 remains the Java/Rust fallback when the tag is absent"
    );
    assert_eq!(
        rust_ins_from_bi[0],
        bi[0].min(rust_cap),
        "after BI sub, Rust PCR min(BI, rust cache)"
    );
    assert_eq!(
        java_ins_pre_floor,
        bi[0].min(java_cap),
        "after BI, Java PCR min(BI, java cache)"
    );
    assert_eq!(java_ins[0], java_ins_pre_floor, "IQ floor 6 does not fire");
    assert_eq!(
        rust_ins_from_bi[0], java_ins[0],
        "PCR cache matches after BI substitution"
    );
    assert_eq!(
        rust_ins[bases.len() - 1],
        bi[bases.len() - 1],
        "last base keeps BI (never PCR-written)"
    );
}
