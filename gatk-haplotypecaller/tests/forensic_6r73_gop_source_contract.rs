//! 6R.73 coordinate-free: production indel GOP source is BAM BI/BD, else Q45.
//!
//! SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`
//! `ReadUtils.getBaseInsertionQualities` / `getBaseDeletionQualities`
//! then `applyPCRErrorModel`.
//!
//! Canonical `TAAGAAAA[0]` BI=44 / cache[2]=40 is **masked** (min=40 either source).
//! The primary proof is an anti-masking cell where BI/BD < cache[r].
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r73_gop_source_contract
//! ```

use gatk_haplotypecaller::indel_gop_from_optional_tag;
use gatk_haplotypecaller::pairhmm_log10::{
    GATK_PARITY_DEFAULT_DEL_QUAL, GATK_PARITY_DEFAULT_GCP, GATK_PARITY_DEFAULT_INS_QUAL,
};
use gatk_haplotypecaller::pairhmm_qual::{
    cap_read_base_qualities, set_to_fixed_value_if_too_low, MIN_USABLE_Q_SCORE,
};
use gatk_haplotypecaller::pcr_error_model::{
    apply_pcr_error_model, error_model_adjusted_qual, find_tandem_repeat_units,
    tandem_repeat_units, PcrErrorModel,
};

const LEGACY_Q45: u8 = 45;

fn cache_at(repeat: usize) -> u8 {
    error_model_adjusted_qual(repeat, PcrErrorModel::Conservative.rate_factor().unwrap())
}

fn pcr(seq: &[u8], ins: &mut [u8], del: &mut [u8]) {
    apply_pcr_error_model(seq, ins, del, PcrErrorModel::Conservative);
}

fn floor_iq_dq(ins: &mut [u8], del: &mut [u8]) {
    for q in ins.iter_mut().chain(del.iter_mut()) {
        *q = set_to_fixed_value_if_too_low(*q, MIN_USABLE_Q_SCORE, MIN_USABLE_Q_SCORE);
    }
}

/// Unique tetranucleotide: finder length 1, cache[1]=40. PCR writes all but last.
const ANTI_SEQ: &[u8] = b"ACGT";
const ANTI_BI: u8 = 30;
const ANTI_BD: u8 = 25;

#[test]
fn anti_masking_sequence_has_cache_above_bi_and_bd() {
    assert_eq!(ANTI_SEQ.len(), 4);
    let r0 = tandem_repeat_units(ANTI_SEQ, 0);
    assert!(r0 >= 1);
    assert_eq!(cache_at(r0), 40);
    assert!(ANTI_BI < 40 && ANTI_BD < 40);
    assert!(ANTI_BI >= MIN_USABLE_Q_SCORE && ANTI_BD >= MIN_USABLE_Q_SCORE);
    assert_ne!(ANTI_BI.min(cache_at(r0)), LEGACY_Q45.min(cache_at(r0)));
    assert_ne!(ANTI_BD.min(cache_at(r0)), LEGACY_Q45.min(cache_at(r0)));
}

#[test]
fn legacy_q45_fill_is_the_pre_fix_divergence() {
    let n = ANTI_SEQ.len();
    let bi = vec![ANTI_BI; n];
    let bd = vec![ANTI_BD; n];
    let mut rust_legacy_ins = vec![LEGACY_Q45; n];
    let mut rust_legacy_del = vec![LEGACY_Q45; n];
    pcr(ANTI_SEQ, &mut rust_legacy_ins, &mut rust_legacy_del);
    let mut java_ins = bi.clone();
    let mut java_del = bd.clone();
    pcr(ANTI_SEQ, &mut java_ins, &mut java_del);
    assert_eq!(java_ins[0], ANTI_BI, "Java min(30, cache[1]=40)=30");
    assert_eq!(java_del[0], ANTI_BD, "Java min(25, cache[1]=40)=25");
    assert_eq!(rust_legacy_ins[0], 40, "legacy min(45, 40)=40");
    assert_eq!(rust_legacy_del[0], 40);
    assert_ne!(rust_legacy_ins[0], java_ins[0]);
    assert_ne!(rust_legacy_del[0], java_del[0]);
}

#[test]
fn production_helper_matches_java_when_both_tags_present() {
    let n = ANTI_SEQ.len();
    let bi = vec![ANTI_BI; n];
    let bd = vec![ANTI_BD; n];
    let mut ins = indel_gop_from_optional_tag(Some(&bi), n).unwrap();
    let mut del = indel_gop_from_optional_tag(Some(&bd), n).unwrap();
    assert_eq!(ins[0], ANTI_BI, "pre-PCR insertion GOP is BI");
    assert_eq!(del[0], ANTI_BD, "pre-PCR deletion GOP is BD");
    assert_ne!(ins[0], del[0], "BI and BD are independent");
    pcr(ANTI_SEQ, &mut ins, &mut del);
    assert_eq!(ins[0], ANTI_BI);
    assert_eq!(del[0], ANTI_BD);
    assert_eq!(ins[n - 1], ANTI_BI, "last base is never PCR-written");
    assert_eq!(del[n - 1], ANTI_BD);
}

#[test]
fn case_b_bi_present_bd_absent() {
    let n = ANTI_SEQ.len();
    let bi = vec![ANTI_BI; n];
    let mut ins = indel_gop_from_optional_tag(Some(&bi), n).unwrap();
    let mut del = indel_gop_from_optional_tag(None, n).unwrap();
    assert_eq!(ins, vec![ANTI_BI; n]);
    assert_eq!(del, vec![GATK_PARITY_DEFAULT_DEL_QUAL; n]);
    pcr(ANTI_SEQ, &mut ins, &mut del);
    assert_eq!(ins[0], ANTI_BI);
    assert_eq!(del[0], 40, "Q45 deletion is PCR-capped to cache[1]");
}

#[test]
fn case_c_bi_absent_bd_present() {
    let n = ANTI_SEQ.len();
    let bd = vec![ANTI_BD; n];
    let mut ins = indel_gop_from_optional_tag(None, n).unwrap();
    let mut del = indel_gop_from_optional_tag(Some(&bd), n).unwrap();
    assert_eq!(ins, vec![GATK_PARITY_DEFAULT_INS_QUAL; n]);
    assert_eq!(del, vec![ANTI_BD; n]);
    pcr(ANTI_SEQ, &mut ins, &mut del);
    assert_eq!(ins[0], 40);
    assert_eq!(del[0], ANTI_BD);
}

#[test]
fn case_d_both_absent_falls_back_to_q45() {
    let n = ANTI_SEQ.len();
    let ins = indel_gop_from_optional_tag(None, n).unwrap();
    let del = indel_gop_from_optional_tag(None, n).unwrap();
    assert_eq!(ins, vec![45u8; n]);
    assert_eq!(del, vec![45u8; n]);
}

#[test]
fn length_mismatch_is_not_silent_q45() {
    let err = indel_gop_from_optional_tag(Some(&[30u8, 31]), 4).unwrap_err();
    let msg = format!("{err}");
    assert!(
        msg.contains("length"),
        "Java createQualityModifiedRead rejects length mismatch: {msg}"
    );
}

#[test]
fn first_differing_cell_anti_masking_pipeline() {
    let seq = ANTI_SEQ;
    let n = seq.len();
    let bi = vec![ANTI_BI; n];
    let bd = vec![ANTI_BD; n];
    let bq = vec![31u8; n];
    let mapq = 60u8;
    let (unit, repeat) = find_tandem_repeat_units(seq, 0);
    let cap = cache_at(repeat);

    let rust_pre_ins = indel_gop_from_optional_tag(Some(&bi), n).unwrap();
    let rust_pre_del = indel_gop_from_optional_tag(Some(&bd), n).unwrap();
    let mut rust_ins = rust_pre_ins.clone();
    let mut rust_del = rust_pre_del.clone();
    pcr(seq, &mut rust_ins, &mut rust_del);

    let mut java_ins = bi.clone();
    let mut java_del = bd.clone();
    pcr(seq, &mut java_ins, &mut java_del);
    let java_ins_pre_floor = java_ins.clone();
    let java_del_pre_floor = java_del.clone();
    floor_iq_dq(&mut java_ins, &mut java_del);

    let mut final_bq = bq.clone();
    cap_read_base_qualities(&mut final_bq, mapq, 18, false);

    eprintln!(
        "cell0 unit={} repeat={repeat} cache={cap} BI={ANTI_BI} BD={ANTI_BD} \
         rust_pre_ins={} rust_pre_del={} rust_ins={} rust_del={} \
         java_ins={} java_del={} bq={} gcp={}",
        std::str::from_utf8(&unit).unwrap_or("?"),
        rust_pre_ins[0],
        rust_pre_del[0],
        rust_ins[0],
        rust_del[0],
        java_ins[0],
        java_del[0],
        final_bq[0],
        GATK_PARITY_DEFAULT_GCP,
    );

    assert_eq!(rust_pre_ins[0], ANTI_BI);
    assert_eq!(rust_pre_del[0], ANTI_BD);
    assert_eq!(rust_ins[0], java_ins_pre_floor[0]);
    assert_eq!(rust_del[0], java_del_pre_floor[0]);
    assert_eq!(
        java_ins[0], java_ins_pre_floor[0],
        "IQ floor 6 does not fire"
    );
    assert_eq!(java_del[0], java_del_pre_floor[0]);
    assert_eq!(final_bq[0], 31);
    assert_eq!(GATK_PARITY_DEFAULT_GCP, 10);
    // First differing value vs legacy Q45 fill is insertion GOP before PCR: 30 vs 45.
    assert_ne!(rust_pre_ins[0], LEGACY_Q45);
}

#[test]
fn canonical_taagaaaa_is_masked_but_sources_differed() {
    let seq = b"TAAGAAAA";
    let (unit, nrep) = find_tandem_repeat_units(seq, 0);
    assert_eq!(unit, b"A");
    assert_eq!(nrep, 2);
    let cap = cache_at(2);
    assert_eq!(cap, 40);
    let bi = 44u8;
    assert_eq!(bi.min(cap), 40);
    assert_eq!(LEGACY_Q45.min(cap), 40);
    assert_eq!(
        bi.min(cap),
        LEGACY_Q45.min(cap),
        "canonical final GOP is masked"
    );
    assert_ne!(bi, LEGACY_Q45, "pre-PCR sources are not equivalent");
    let mut from_bi = vec![bi; seq.len()];
    let mut from_q45 = vec![LEGACY_Q45; seq.len()];
    let mut d1 = from_bi.clone();
    let mut d2 = from_q45.clone();
    pcr(seq, &mut from_bi, &mut d1);
    pcr(seq, &mut from_q45, &mut d2);
    assert_eq!(from_bi[0], 40);
    assert_eq!(from_q45[0], 40);
    let prod = indel_gop_from_optional_tag(Some(&[bi; 8]), 8).unwrap();
    assert_eq!(prod[0], 44, "production source after 6R.73 is BI");
}
