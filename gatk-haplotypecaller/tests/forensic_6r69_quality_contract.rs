//! 6R.69 coordinate-free: complete `modifyReadQualities` contract vs Rust PairHMM prep.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//! `PairHMMLikelihoodCalculationEngine.modifyReadQualities` order:
//! keep/hard-clip softclips → clone BQ + `BI`/`BD` (else Q45) → PCR min-cap →
//! cap BQ by MAPQ/threshold and floor IQ/DQ at `MIN_USABLE_Q_SCORE` (6).
//! GCP is assigned later by the imputator (`gcpHMM=10`), not inside
//! `modifyReadQualities`.
//!
//! 6R.68 first divergence (BI/BD vs Q45) is closed in production by 6R.73.
//! After BI/BD substitution the PCR min-cap matches Java (6R.72).
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r69_quality_contract
//! ```

use gatk_haplotypecaller::pairhmm_log10::{
    GATK_PARITY_DEFAULT_DEL_QUAL, GATK_PARITY_DEFAULT_GCP, GATK_PARITY_DEFAULT_INS_QUAL,
};
use gatk_haplotypecaller::pairhmm_qual::{
    cap_read_base_qualities, set_to_fixed_value_if_too_low, MIN_USABLE_Q_SCORE,
};
use gatk_haplotypecaller::pcr_error_model::{
    apply_pcr_error_model, error_model_adjusted_qual, tandem_repeat_units, PcrErrorModel,
};
use gatk_haplotypecaller::{
    indel_gop_from_optional_tag, prepare_read_quals_for_pairhmm_inplace, HcLikelihoodEngineConfig,
};

fn java_fast_round(d: f64) -> i32 {
    if d > 0.0 {
        (d + 0.5) as i32
    } else {
        (d - 0.5) as i32
    }
}

/// GATK 4.4 `getErrorModelAdjustedQual` (`INITIAL_QSCORE=40`, min 10).
fn java_44_adjusted_qual(repeat_length: usize, rate_factor: f64) -> u8 {
    let q = 40.0 - (repeat_length as f64 / (rate_factor * std::f64::consts::PI)).exp() + 1.0;
    java_fast_round(q).max(10) as u8
}

/// Java `applyPCRErrorModel` with CONSERVATIVE rate 3.0, using the production
/// finder so the *only* isolated difference is the cache.
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

fn java_floor_indel_quals(ins: &mut [u8], del: &mut [u8]) {
    for q in ins.iter_mut().chain(del.iter_mut()) {
        *q = set_to_fixed_value_if_too_low(*q, MIN_USABLE_Q_SCORE, MIN_USABLE_Q_SCORE);
    }
}

/// Java `ReadUtils.getBaseInsertionQualities`: BI FastQ → Phred, else fill Q45.
fn java_indel_gop_from_bi_or_q45(bi: Option<&[u8]>, read_len: usize) -> Vec<u8> {
    match bi {
        Some(tag) if tag.len() == read_len => tag.to_vec(),
        Some(_) => {
            // Length mismatch: Java `createQualityModifiedRead` throws. No silent fill.
            panic!("Java contract: BI/BD length must equal read length")
        }
        None => vec![GATK_PARITY_DEFAULT_INS_QUAL; read_len],
    }
}

#[test]
fn missing_bi_bd_falls_back_to_q45() {
    let ins = java_indel_gop_from_bi_or_q45(None, 8);
    assert_eq!(ins, vec![45u8; 8]);
    assert_eq!(GATK_PARITY_DEFAULT_INS_QUAL, 45);
    assert_eq!(GATK_PARITY_DEFAULT_DEL_QUAL, 45);
}

#[test]
fn bi_present_is_preferred_over_q45() {
    let bi = [44u8, 42, 43, 41, 45, 40, 46, 38];
    let java = java_indel_gop_from_bi_or_q45(Some(&bi), bi.len());
    let rust = indel_gop_from_optional_tag(Some(&bi), bi.len()).unwrap();
    let q45 = vec![GATK_PARITY_DEFAULT_INS_QUAL; bi.len()];
    assert_eq!(java[0], 44);
    assert_eq!(rust, java);
    assert_ne!(java.as_slice(), q45.as_slice());
}

#[test]
#[should_panic(expected = "BI/BD length must equal read length")]
fn java_rejects_length_mismatched_bi() {
    let _ = java_indel_gop_from_bi_or_q45(Some(&[44u8, 42]), 8);
}

/// After substituting BI/BD, PCR is the next Java operation. Same STR length
/// (6R.71), same starting GOP (44), matching cache → matching min-cap (6R.72).
#[test]
fn after_bi_substitution_pcr_cache_matches_java() {
    let read = b"ACGTACGTACGTACGT";
    let bi = vec![44u8; read.len()];
    assert_eq!(tandem_repeat_units(read, 0), 3);
    assert_eq!(java_44_adjusted_qual(3, 3.0), 40);
    assert_eq!(
        error_model_adjusted_qual(3, PcrErrorModel::Conservative.rate_factor().unwrap()),
        40
    );

    let mut rust_after_bi = bi.clone();
    let mut rust_del = bi.clone();
    apply_pcr_error_model(
        read,
        &mut rust_after_bi,
        &mut rust_del,
        PcrErrorModel::Conservative,
    );

    let mut java_after_bi = bi.clone();
    let mut java_del = bi.clone();
    apply_java_44_conservative_pcr(read, &mut java_after_bi, &mut java_del);

    assert_eq!(
        rust_after_bi[0], 40,
        "Rust cache[3]=40 caps BI 44 to 40: {}",
        rust_after_bi[0]
    );
    assert_eq!(
        java_after_bi[0], 40,
        "Java cache[3]=40 caps BI 44 to 40: {}",
        java_after_bi[0]
    );
    java_floor_indel_quals(&mut java_after_bi, &mut java_del);
    assert_eq!(
        java_after_bi[0], 40,
        "Java IQ floor at 6 does not fire on Q40"
    );
}

/// PCR is a per-position MIN against a repeat-length cache, not a transform of BI
/// through an error probability.
#[test]
fn java_pcr_is_a_min_cap_from_repeat_cache() {
    let bi = 44u8;
    let java_cap = java_44_adjusted_qual(1, 3.0);
    let rust_cap = error_model_adjusted_qual(1, PcrErrorModel::Conservative.rate_factor().unwrap());
    assert_eq!(java_cap, 40);
    assert_eq!(rust_cap, 40);
    assert_eq!(bi.min(java_cap), 40);
    assert_eq!(bi.min(rust_cap), 40);
}

#[test]
fn last_base_is_not_pcr_capped() {
    let read = b"ACGT";
    let mut ins = vec![44u8; 4];
    let mut del = vec![44u8; 4];
    apply_java_44_conservative_pcr(read, &mut ins, &mut del);
    assert_eq!(ins[3], 44, "Java PCR loop never writes the last index");
    apply_pcr_error_model(read, &mut ins, &mut del, PcrErrorModel::Conservative);
    assert_eq!(ins[3], 44);
}

#[test]
fn java_indel_floor_fires_only_below_min_usable() {
    let mut ins = vec![5u8, 40];
    let mut del = vec![5u8, 40];
    java_floor_indel_quals(&mut ins, &mut del);
    assert_eq!(ins, vec![MIN_USABLE_Q_SCORE, 40]);
    assert_eq!(del, vec![MIN_USABLE_Q_SCORE, 40]);
}

#[test]
fn bq_cap_is_independent_of_indel_gop() {
    let mut quals = vec![5u8, 18, 30, 40];
    let cfg = HcLikelihoodEngineConfig::default();
    assert_eq!(cfg.base_quality_score_threshold, 18);
    prepare_read_quals_for_pairhmm_inplace(&mut quals, 25, &cfg);
    // MAPQ 25 caps 30 and 40; 5 < threshold 18 → 6; 18 stays.
    assert_eq!(quals, vec![6, 18, 25, 25]);
    let mut q2 = vec![5u8, 18, 30, 40];
    cap_read_base_qualities(&mut q2, 25, 18, false);
    assert_eq!(q2, quals);
}

#[test]
fn gcp_is_constant_ten_not_from_bi() {
    assert_eq!(GATK_PARITY_DEFAULT_GCP, 10);
    let gcp = vec![GATK_PARITY_DEFAULT_GCP; 8];
    assert!(gcp.iter().all(|&q| q == 10));
}

#[test]
fn rust_production_q45_fill_is_masked_by_q40_cache() {
    let read = b"ACGTACGT";
    let bi = vec![44u8; read.len()];
    assert_eq!(GATK_PARITY_DEFAULT_INS_QUAL, 45);
    assert_ne!(bi[0], GATK_PARITY_DEFAULT_INS_QUAL);
    let mut rust = vec![GATK_PARITY_DEFAULT_INS_QUAL; read.len()];
    let mut rust_del = rust.clone();
    apply_pcr_error_model(read, &mut rust, &mut rust_del, PcrErrorModel::Conservative);
    let mut java = bi.clone();
    let mut java_del = bi.clone();
    apply_java_44_conservative_pcr(read, &mut java, &mut java_del);
    assert_eq!(rust[0], 40, "Q45 PCR-capped to cache");
    assert_eq!(java[0], 40, "BI 44 PCR-capped to cache");
    assert_eq!(
        rust[read.len() - 1],
        45,
        "last base keeps production Q45 fill"
    );
    assert_eq!(java[read.len() - 1], 44, "last base keeps BI");
}
