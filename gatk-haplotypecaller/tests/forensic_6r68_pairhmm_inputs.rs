//! 6R.68 coordinate-free: indel GOP source is BAM BI/BD, else Q45.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//! `ReadUtils.getBaseInsertionQualities` / `getBaseDeletionQualities`.
//! 6R.73: production `score_read_against_haplotypes` consumes those tags.
//! Layer B (PairHMM kernel) is not entered.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r68_pairhmm_inputs
//! ```

use gatk_haplotypecaller::bio_ids::{HaplotypeIndex, ReadIndex};
use gatk_haplotypecaller::pairhmm_log10::{
    GATK_PARITY_DEFAULT_DEL_QUAL, GATK_PARITY_DEFAULT_GCP, GATK_PARITY_DEFAULT_INS_QUAL,
};
use gatk_haplotypecaller::pairhmm_qual::MIN_USABLE_Q_SCORE;
use gatk_haplotypecaller::pcr_error_model::{
    apply_pcr_error_model, error_model_adjusted_qual, tandem_repeat_units, PcrErrorModel,
};
use gatk_haplotypecaller::region_read_likelihood::RegionReadLikelihood;
use gatk_haplotypecaller::{
    indel_gop_from_optional_tag, logless_pairhmm_likelihood,
    prepare_read_quals_for_pairhmm_inplace, HcLikelihoodEngineConfig,
};
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

/// GATK 4.4 `MathUtils.fastRound`.
fn java_fast_round(d: f64) -> i32 {
    if d > 0.0 {
        (d + 0.5) as i32
    } else {
        (d - 0.5) as i32
    }
}

/// GATK 4.4.0.0 `getErrorModelAdjustedQual`.
fn java_44_adjusted_qual(repeat_length: usize, rate_factor: f64) -> u8 {
    let q = 40.0 - (repeat_length as f64 / (rate_factor * std::f64::consts::PI)).exp() + 1.0;
    java_fast_round(q).max(10) as u8
}

/// Same loop as Java `applyPCRErrorModel` / Rust `apply_pcr_error_model`, with
/// the pinned Java CONSERVATIVE cache. Repeat lengths use the production
/// finder so the isolated difference is the cache constants.
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

fn rust_conservative_gop(read_bases: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let mut ins = vec![GATK_PARITY_DEFAULT_INS_QUAL; read_bases.len()];
    let mut del = vec![GATK_PARITY_DEFAULT_DEL_QUAL; read_bases.len()];
    apply_pcr_error_model(read_bases, &mut ins, &mut del, PcrErrorModel::Conservative);
    (ins, del)
}

fn java_conservative_gop(read_bases: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let mut ins = vec![GATK_PARITY_DEFAULT_INS_QUAL; read_bases.len()];
    let mut del = vec![GATK_PARITY_DEFAULT_DEL_QUAL; read_bases.len()];
    apply_java_44_conservative_pcr(read_bases, &mut ins, &mut del);
    (ins, del)
}

fn first_gop_diff(a: &[u8], b: &[u8]) -> Option<usize> {
    a.iter().zip(b.iter()).position(|(x, y)| x != y)
}

/// Stable diagnostic row: order, identity hashes, lengths, likelihood.
fn matrix_row_id(
    read_index: usize,
    hap_index: usize,
    read_bases: &[u8],
    hap_bases: &[u8],
    log10_likelihood: f64,
) -> String {
    let mut rh = DefaultHasher::new();
    read_bases.hash(&mut rh);
    let mut hh = DefaultHasher::new();
    hap_bases.hash(&mut hh);
    format!(
        "r{read_index}\th{hap_index}\tread_len={}\thap_len={}\tread_hash={:016x}\thap_hash={:016x}\tll={:.17}",
        read_bases.len(),
        hap_bases.len(),
        rh.finish(),
        hh.finish(),
        log10_likelihood
    )
}

#[test]
fn java_and_production_use_bi_when_present() {
    let bi = vec![42u8, 43, 41, 45, 40, 46, 38, 44];
    let rust = indel_gop_from_optional_tag(Some(&bi), bi.len()).unwrap();
    let legacy = vec![GATK_PARITY_DEFAULT_INS_QUAL; bi.len()];
    assert_eq!(GATK_PARITY_DEFAULT_INS_QUAL, 45);
    assert_eq!(rust, bi);
    assert_ne!(bi, legacy, "legacy constant Q45 remains a distinct source");
    assert_eq!(first_gop_diff(&bi, &legacy), Some(0));
    assert_eq!(bi[0], 42);
    assert_eq!(legacy[0], 45);
}

#[test]
fn same_read_hap_bq_gcp_yields_different_ll_when_only_indel_gop_source_changes() {
    let read = b"ACGTACGTACGTACGT";
    let hap = b"ACGTACGTACGTACGTACGT";
    let mut quals = vec![30u8; read.len()];
    let cfg = HcLikelihoodEngineConfig::default();
    prepare_read_quals_for_pairhmm_inplace(&mut quals, 60, &cfg);
    let gcp = vec![GATK_PARITY_DEFAULT_GCP; read.len()];
    // Pre-PCR GOP: BI-like 42 vs constant Q45 still moves the kernel (sensitivity, not production).
    let rust_gop = vec![GATK_PARITY_DEFAULT_INS_QUAL; read.len()];
    let java_bi_gop = vec![42u8; read.len()];
    let rust_ll =
        logless_pairhmm_likelihood(read, &quals, hap, &rust_gop, &rust_gop, &gcp).expect("q45");
    let java_ll = logless_pairhmm_likelihood(read, &quals, hap, &java_bi_gop, &java_bi_gop, &gcp)
        .expect("bi");
    assert!(
        (rust_ll - java_ll).abs() > 1e-9,
        "BI GOP vs Q45 must move PairHMM output: rust={rust_ll:.17} java_bi={java_ll:.17}"
    );
}

/// 6R.72: CONSERVATIVE cache matches Java. BI/BD vs Q45 remains first.
#[test]
fn java_4_4_conservative_cache_matches_rust_at_every_short_repeat() {
    let rust_rate = PcrErrorModel::Conservative.rate_factor().unwrap();
    assert_eq!(rust_rate, 3.0);
    for repeat in 0..=10 {
        let rust = error_model_adjusted_qual(repeat, rust_rate);
        let java = java_44_adjusted_qual(repeat, 3.0);
        assert_eq!(
            rust, java,
            "repeat {repeat}: Rust CONSERVATIVE {rust} vs Java 4.4 {java}"
        );
    }
    assert_eq!(java_44_adjusted_qual(0, 3.0), 40);
    assert_eq!(error_model_adjusted_qual(0, rust_rate), 40);
}

#[test]
fn pcr_gop_matches_after_shared_q45_fill() {
    // Shared Q45 fill + matching cache → matching post-PCR GOP. That masks BI vs Q45.
    let read = b"ACGTACGTACGTACGT";
    let (rust_ins, rust_del) = rust_conservative_gop(read);
    let (java_ins, java_del) = java_conservative_gop(read);
    assert_eq!(first_gop_diff(&rust_ins, &java_ins), None);
    assert_eq!(first_gop_diff(&rust_del, &java_del), None);
    assert_eq!(rust_ins[0], 40);
    assert_eq!(java_ins[0], 40);
    assert_eq!(rust_del[0], 40);
    assert_eq!(java_del[0], 40);
    assert_eq!(
        rust_ins[read.len() - 1],
        GATK_PARITY_DEFAULT_INS_QUAL,
        "last base is never PCR-written; production fill stays Q45"
    );
}

#[test]
fn bq_cap_matches_java_min_usable_before_gop() {
    let mut quals = vec![5u8, 18, 30];
    let cfg = HcLikelihoodEngineConfig::default();
    prepare_read_quals_for_pairhmm_inplace(&mut quals, 60, &cfg);
    assert_eq!(quals, vec![MIN_USABLE_Q_SCORE, 18, 30]);
}

#[test]
fn likelihood_matrix_dump_is_ordered_by_read_then_haplotype_index() {
    let reads: [&[u8]; 2] = [b"ACGT", b"TGCA"];
    let haps: [&[u8]; 2] = [b"ACGTAC", b"TGCATG"];
    let mut cells = vec![
        RegionReadLikelihood {
            read_index: ReadIndex::new(1),
            haplotype_index: HaplotypeIndex::new(0),
            log10_likelihood: -2.0,
        },
        RegionReadLikelihood {
            read_index: ReadIndex::new(0),
            haplotype_index: HaplotypeIndex::new(1),
            log10_likelihood: -3.0,
        },
        RegionReadLikelihood {
            read_index: ReadIndex::new(0),
            haplotype_index: HaplotypeIndex::new(0),
            log10_likelihood: -1.0,
        },
        RegionReadLikelihood {
            read_index: ReadIndex::new(1),
            haplotype_index: HaplotypeIndex::new(1),
            log10_likelihood: -4.0,
        },
    ];
    cells.sort_by_key(|c| (c.read_index.get(), c.haplotype_index.get()));
    let dump: Vec<String> = cells
        .iter()
        .map(|c| {
            matrix_row_id(
                c.read_index.get(),
                c.haplotype_index.get(),
                reads[c.read_index.get()],
                haps[c.haplotype_index.get()],
                c.log10_likelihood,
            )
        })
        .collect();
    for window in dump.windows(2) {
        let left = window[0].as_str();
        let right = window[1].as_str();
        assert!(
            left < right,
            "matrix dump must be lexicographic in r then h:\n{left}\n{right}"
        );
    }
    assert!(dump[0].starts_with("r0\th0\t"));
    assert!(dump[1].starts_with("r0\th1\t"));
    assert!(dump[2].starts_with("r1\th0\t"));
    assert!(dump[3].starts_with("r1\th1\t"));
}
