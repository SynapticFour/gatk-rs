//! 6R.278: measurement-only. GKL AVX last-stripe `sumX = VEC_ADD(sumX, X_t.d)`
//! vs the Rust f64 fold of final-row X at `20:29455649 T/TGTTTG`.
//! PRODUCTION CHANGE: NONE.
//!
//! `X_t` is the completed 6R.267 f32 add. The M/X/Y recurrence and the
//! 6R.277 sumM fold stay Rust. This round injects only the first both-nonzero
//! sumX add on the extracted last-row lane. Native AVX GKL is not executed
//! (oracle x86_64; this host aarch64).
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r278_gkl_avx_last_stripe_sumx_add -- --nocapture --test-threads=1
//! HOLDOUT_6R278=1 cargo test -p gatk-haplotypecaller --test holdout_6r278_gkl_avx_last_stripe_sumx_add -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::assembly_region_finalize::{
    clip_finalized_reads_to_region, finalize_region_reads_for_assembly,
    gatk_min_tail_quality_for_assembly,
};
use gatk_haplotypecaller::hc_genotyping_engine::take_colocated_merge_numerics;
use gatk_haplotypecaller::pairhmm_log10::GATK_PARITY_DEFAULT_GCP;
use gatk_haplotypecaller::pcr_error_model::apply_pcr_error_model;
use gatk_haplotypecaller::read_threading_assembler::DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH;
use gatk_haplotypecaller::{
    begin_hap_list_observe, biallelic_genotype_log10_likelihoods_gatk, call_disposition,
    flatten_assembly_regions, indel_gop_from_optional_tag, logless_pairhmm_likelihood,
    prepare_read_quals_for_pairhmm_inplace, region_likelihoods_to_rows,
    score_read_against_haplotypes, take_hap_list_trim_span, traverse_assembly_region_walker,
    AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition, HaplotypeCallerEngine,
    HcLikelihoodEngineConfig, PairHmmBackend, ReadFilterParams, ReadLikelihoodRow,
    WalkerTraversalConfig, INITIAL_CONDITION, INITIAL_CONDITION_LOG10,
};
use rust_htslib::bam::record::Aux;
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const GKL_PIN: &str = "0.8.8 IntelLabs/GKL avx-pairhmm-template.h:357 sumX = VEC_ADD(sumX, X_t.d)";
const GKL_INITIAL_F32: f32 = f32::from_bits((127 + 120) << 23);
const Q20_GKL_PH2PR_F32: u32 = 0x3c23d70a;
const Q20_GKL_MATCH_F32: u32 = 0x3f7d70a4;
const Q20_GKL_MISMATCH_F32: u32 = 0x3b5a740d;
const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 29_455_649;
const TARGET_REF: &str = "T";
const TARGET_ALT: &str = "TGTTTG";
const FROZEN_IDX: [usize; 5] = [19, 20, 21, 22, 23];
const RUST_FNV: [&str; 5] = [
    "55012fcf3b430591",
    "341b2e070ccb5846",
    "6a65e4c02733c2ed",
    "fa07750bf228b3c2",
    "79451c576721a729",
];
const JAVA_FNV: [&str; 5] = [
    "fcd72c6ce600dd16",
    "c1a4e8204522f645",
    "eb03271fa7548f26",
    "7dc2a8ae5da116e1",
    "343e7c443c1d9732",
];
const JAVA_PREFIX: &[u8] = b"CAAAGAGTA";
const JAVA_SUFFIX: &[u8] = b"TAAA";
const FLOOR: f64 = -4.5;
const JOINT_GL: f64 = -351.91571812977724676;
const JAVA_CLIP: (u64, u64) = (29_455_560, 29_455_728);
const RUST_CLIP: (u64, u64) = (29_455_569, 29_455_724);
const JAVA_TSV: &str = include_str!("forensic_6r252_java_hap_trim.tsv");
const BOUNDARY: f64 = 3517.5;
const REQUIRED_MOVE: f64 = 1.65718129777269496;
const GKL_MAX_QUAL: usize = 254;
const JACOBIAN_STEP: f64 = 0.0001;
const JACOBIAN_TOL: f32 = 8.0;
const JACOBIAN_INV_STEP: f32 = 10_000.0;
const JACOBIAN_SIZE: usize = 80_001;
const GKL_INV_LN10: f64 = 0.434294;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R278\t{key}\t{}", value.as_ref());
}

fn fmt_f64(x: f64) -> String {
    format!("{x:.17} bits=0x{:016x}", x.to_bits())
}

fn fmt_f32(x: f32) -> String {
    format!(
        "{x:.9} bits=0x{:08x} widened={}",
        x.to_bits(),
        fmt_f64(f64::from(x))
    )
}

fn fnv1a64_hex(bases: &[u8]) -> String {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in bases {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    format!("{h:016x}")
}

fn java_seq(hash: &str) -> Vec<u8> {
    for line in JAVA_TSV.lines() {
        if line.contains("stage=trimmed") && line.contains(&format!("hash={hash}")) {
            if let Some(seq) = line.split("\tseq=").nth(1) {
                return seq.as_bytes().to_vec();
            }
        }
    }
    panic!("missing Java trimmed seq for {hash}");
}

fn bam_indel_phred(rec: &rust_htslib::bam::Record, tag: &[u8]) -> Option<Vec<u8>> {
    match rec.aux(tag) {
        Ok(Aux::String(s)) => Some(s.bytes().map(|b| b.saturating_sub(33)).collect()),
        _ => None,
    }
}

fn floor_pair(l0: f64, l2: f64) -> (f64, f64) {
    let best = l0.max(l2);
    if !best.is_finite() {
        return (l0, l2);
    }
    let floor = best + FLOOR;
    let lift = |v: f64| {
        if v.is_finite() && v < floor {
            floor
        } else {
            v
        }
    };
    (lift(l0), lift(l2))
}

fn biallelic_gls(l0: &[f64], l2: &[f64]) -> Vec<f64> {
    let rows: Vec<ReadLikelihoodRow> = l0
        .iter()
        .zip(l2.iter())
        .enumerate()
        .map(|(i, (a, b))| {
            let (fa, fb) = floor_pair(*a, *b);
            ReadLikelihoodRow {
                read_index: i,
                read_id: String::new(),
                haplotype_log10_likelihoods: vec![fa, fb],
            }
        })
        .collect();
    biallelic_genotype_log10_likelihoods_gatk(&rows, 0, 1)
}

fn emitted_hom_alt(gls: &[f64]) -> (f64, f64, i32) {
    let best = gls.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let rel = gls[2] - best;
    let cont = -10.0 * rel;
    let pl = (cont + 0.5).floor() as i32;
    (rel, cont, pl)
}

fn floor5(raw: [f64; 5], other_best: f64) -> [f64; 5] {
    let raw_max = raw.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let best = other_best.max(raw_max);
    let floor = best + FLOOR;
    std::array::from_fn(|k| {
        if raw[k].is_finite() && raw[k] < floor {
            floor
        } else {
            raw[k]
        }
    })
}

fn clip_map(
    finalized: &[rust_htslib::bam::Record],
    region: &gatk_haplotypecaller::assembly_region_iterator::AssemblyRegion,
    start: u64,
    end: u64,
) -> std::collections::HashMap<(Vec<u8>, u16), rust_htslib::bam::Record> {
    let mut clip_region = region.clone();
    clip_region.extended_start = GenomePosition::new_1based(start);
    clip_region.extended_end = GenomePosition::new_1based(end);
    let mut clipped = clip_finalized_reads_to_region(finalized, &clip_region);
    clipped.retain(|r| r.seq().len() >= 10);
    let mut by_id = std::collections::HashMap::new();
    for rec in clipped {
        by_id.insert((rec.qname().to_vec(), rec.flags()), rec);
    }
    by_id
}

struct Planes {
    bases: Vec<u8>,
    bq: Vec<u8>,
    iq: Vec<u8>,
    dq: Vec<u8>,
    gcp: Vec<u8>,
}

fn planes(rec: &rust_htslib::bam::Record, cfg: &HcLikelihoodEngineConfig) -> Planes {
    let bases = rec.seq().as_bytes();
    let n = bases.len();
    let mut bq = rec.qual().to_vec();
    prepare_read_quals_for_pairhmm_inplace(&mut bq, rec.mapq(), cfg);
    let mut iq = indel_gop_from_optional_tag(bam_indel_phred(rec, b"BI").as_deref(), n).unwrap();
    let mut dq = indel_gop_from_optional_tag(bam_indel_phred(rec, b"BD").as_deref(), n).unwrap();
    apply_pcr_error_model(&bases, &mut iq, &mut dq, cfg.pcr_error_model);
    Planes {
        bases,
        bq,
        iq,
        dq,
        gcp: vec![GATK_PARITY_DEFAULT_GCP; n],
    }
}

fn rust_err(q: u8) -> f64 {
    10f64.powf(-(q as f64) / 10.0)
}

fn rust_match(q: u8) -> f64 {
    1.0 - rust_err(q)
}

fn rust_approx_log10_sum(a: f64, b: f64) -> f64 {
    let (x, y) = if a > b { (b, a) } else { (a, b) };
    if x.is_infinite() && x.is_sign_negative() {
        return y;
    }
    y + (1.0 + 10f64.powf(x - y)).log10()
}

fn rust_m2m(ins: u8, del: u8) -> f64 {
    let (min_q, max_q) = if ins <= del {
        (ins as usize, del as usize)
    } else {
        (del as usize, ins as usize)
    };
    let log10_sum = rust_approx_log10_sum(-0.1 * min_q as f64, -0.1 * max_q as f64);
    let log10_m2m = (1.0 - 10f64.powf(log10_sum).min(1.0)).log10();
    10f64.powf(log10_m2m)
}

fn rust_trans(ins: u8, del: u8, gcp: u8) -> [f64; 6] {
    let gcp_err = rust_err(gcp);
    [
        rust_m2m(ins, del),
        1.0 - rust_err(gcp),
        rust_err(ins),
        gcp_err,
        rust_err(del),
        gcp_err,
    ]
}

fn libc_powf(base: f32, exp: f32) -> f32 {
    extern "C" {
        fn powf(x: f32, y: f32) -> f32;
    }
    unsafe { powf(base, exp) }
}

fn gkl_exponent(q: u8) -> f32 {
    -((q as usize & 127) as f32) / 10.0f32
}

fn gkl_ph2pr(q: u8) -> f32 {
    libc_powf(10.0, gkl_exponent(q))
}

/// GKL AVX `_1_distm = VEC_SUB(1.0, distm)` with `NUMBER=float`.
fn gkl_match(q: u8) -> f32 {
    1.0f32 - gkl_ph2pr(q)
}

/// GKL AVX `distm = VEC_DIV(distm, 3.0)` with `NUMBER=float`.
fn gkl_mismatch(q: u8) -> f32 {
    gkl_ph2pr(q) / 3.0f32
}

fn rust_mismatch(q: u8) -> f64 {
    rust_err(q) / 3.0
}

fn bases_equal(x: u8, y: u8) -> bool {
    x == y || x == b'N' || y == b'N'
}

/// GKL `computeDistVec`: `VEC_BLENDV(distmChosen, distm, _1_distm, mask)`.
/// AVX: `_mm256_blendv_ps(distm, _1_distm, mask)` — MSB of mask selects
/// `_1_distm` (match) else `distm` (mismatch). Bit-preserving; no arithmetic.
fn gkl_distm_f32(q: u8, equal: bool) -> f32 {
    if equal {
        gkl_match(q)
    } else {
        gkl_mismatch(q)
    }
}

fn rust_distm(q: u8, equal: bool) -> f64 {
    if equal {
        rust_match(q)
    } else {
        rust_mismatch(q)
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum XCf {
    Baseline,
    FirstBothNonzeroSumxAdd,
}

/// GKL `pGAPM = ctx._(1.0) - ctx.ph2pr[_c]` already f32 (`_c = tc->c[r-1] & 127`).
fn gkl_pgapm(gcp: u8) -> f32 {
    1.0f32 - gkl_ph2pr(gcp)
}

fn rust_pgapm(gcp: u8) -> f64 {
    1.0 - rust_err(gcp)
}

/// GKL `Context<float>::set_mm_prob` via `matchToMatchProb` table (NUMBER=float).
fn gkl_pmm_table() -> &'static [f32] {
    static T: OnceLock<Vec<f32>> = OnceLock::new();
    T.get_or_init(|| {
        let mut jacobian = vec![0.0f32; JACOBIAN_SIZE];
        for k in 0..JACOBIAN_SIZE {
            jacobian[k] = (1.0 + 10f64.powf(-(k as f64) * JACOBIAN_STEP)).log10() as f32;
        }
        let approx = |small: f32, big: f32| -> f32 {
            let (mut s, mut b) = (small, big);
            if s > b {
                std::mem::swap(&mut s, &mut b);
            }
            if s.is_infinite() || b.is_infinite() {
                return b;
            }
            let diff = b - s;
            if diff >= JACOBIAN_TOL {
                return b;
            }
            let d = diff * JACOBIAN_INV_STEP;
            let ind = if d > 0.0 {
                (d + 0.5) as i32
            } else {
                (d - 0.5) as i32
            } as usize;
            b + jacobian[ind.min(JACOBIAN_SIZE - 1)]
        };
        let n = ((GKL_MAX_QUAL + 1) * (GKL_MAX_QUAL + 2)) >> 1;
        let mut m2m = vec![0.0f32; n];
        let mut offset = 0usize;
        for i in 0..=GKL_MAX_QUAL {
            for j in 0..=i {
                let log10_sum = approx(-0.1f32 * i as f32, -0.1f32 * j as f32) as f64;
                let match_log10 =
                    ((-1.0f64 * (10f64.powf(log10_sum)).min(1.0)).ln_1p()) * GKL_INV_LN10;
                m2m[offset + j] = 10f64.powf(match_log10) as f32;
            }
            offset += i + 1;
        }
        m2m
    })
}

fn gkl_pmm(ins_q: u8, del_q: u8) -> f32 {
    let ins = (ins_q as usize) & 127;
    let del = (del_q as usize) & 127;
    let (min_q, max_q) = if ins <= del { (ins, del) } else { (del, ins) };
    gkl_pmm_table()[((max_q * (max_q + 1)) >> 1) + min_q]
}

/// GKL f32 `M_t_2 * pMM` (6R.271) in INITIAL-normalized space.
fn gkl_mm_f32(m_diag: f64, ins_q: u8, del_q: u8) -> f32 {
    (m_diag / INITIAL_CONDITION) as f32 * gkl_pmm(ins_q, del_q)
}

/// GKL f32 `X_t_2 * pGAPM` (6R.272) in INITIAL-normalized space.
fn gkl_xgapm_f32(x_diag: f64, gcp: u8) -> f32 {
    (x_diag / INITIAL_CONDITION) as f32 * gkl_pgapm(gcp)
}

/// GKL f32 `Y_t_2 * pGAPM` (6R.273) in INITIAL-normalized space.
fn gkl_ygapm_f32(y_diag: f64, gcp: u8) -> f32 {
    (y_diag / INITIAL_CONDITION) as f32 * gkl_pgapm(gcp)
}

/// GKL f32 6R.274 partial: `VEC_ADD(M*pMM, X*pGAPM)`.
fn gkl_partial_m_f32(m_diag: f64, x_diag: f64, ins_q: u8, del_q: u8, gcp: u8) -> f32 {
    gkl_mm_f32(m_diag, ins_q, del_q) + gkl_xgapm_f32(x_diag, gcp)
}

/// Exact 6R.275 f32 inner sum `(M*pMM + X*pGAPM) + Y*pGAPM` in ratio space.
fn gkl_inner_f32(m_diag: f64, x_diag: f64, y_diag: f64, ins_q: u8, del_q: u8, gcp: u8) -> f32 {
    gkl_partial_m_f32(m_diag, x_diag, ins_q, del_q, gcp) + gkl_ygapm_f32(y_diag, gcp)
}

/// GKL `IntelPairHmm.cc` sets `_MM_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON)` and does not set DAZ.
/// Subnormal inputs participate; a subnormal add result is stored as +0.
fn gkl_add_ftz(a: f32, b: f32) -> f32 {
    let sum = a + b;
    if sum.is_subnormal() {
        0.0
    } else {
        sum
    }
}

/// GKL `pMX = ctx.ph2pr[ins]` (initializeVectors).
fn gkl_pmx(ins_q: u8) -> f32 {
    gkl_ph2pr(ins_q)
}

/// GKL `pXX = ctx.ph2pr[gcp]` (initializeVectors).
fn gkl_pxx(gcp: u8) -> f32 {
    gkl_ph2pr(gcp)
}

/// IEEE-754 f32 `M_t_1 * pMX` in INITIAL-normalized space (6R.265).
fn gkl_mx_f32(m_up: f64, ins_q: u8) -> f32 {
    (m_up / INITIAL_CONDITION) as f32 * gkl_pmx(ins_q)
}

/// IEEE-754 f32 `X_t_1 * pXX` in INITIAL-normalized space (6R.266).
fn gkl_xx_f32(x_up: f64, gcp: u8) -> f32 {
    (x_up / INITIAL_CONDITION) as f32 * gkl_pxx(gcp)
}

/// Completed `X_t` after 6R.267: `VEC_ADD(M_t_1*pMX, X_t_1*pXX)`.
fn gkl_closed_x_f32(m_up: f64, x_up: f64, ins_q: u8, gcp: u8) -> f32 {
    gkl_mx_f32(m_up, ins_q) + gkl_xx_f32(x_up, gcp)
}

/// Completed `M_t` after 6R.275 and 6R.276, in INITIAL-normalized f32.
fn gkl_closed_m_f32(
    m_diag: f64,
    x_diag: f64,
    y_diag: f64,
    ins_q: u8,
    del_q: u8,
    gcp: u8,
    bq: u8,
    equal: bool,
) -> f32 {
    gkl_inner_f32(m_diag, x_diag, y_diag, ins_q, del_q, gcp) * gkl_distm_f32(bq, equal)
}

fn ulp_f64(a: f64, b: f64) -> u64 {
    a.to_bits().abs_diff(b.to_bits())
}

#[derive(Clone, Copy)]
struct AccCell {
    retain_row: usize,
    hap_index: usize,
    read_len: usize,
    hap_len: usize,
    hap_j: usize,
    prior_nonzero_j: usize,
    prior_gkl_x: f32,
    ins_q: u8,
    gcp: u8,
    m_up: f64,
    x_up: f64,
    gkl_sum_before: f32,
    gkl_x: f32,
    rust_sum_before: f64,
    rust_x: f64,
}

fn rust_diag_ll(
    read: &[u8],
    quals: &[u8],
    hap: &[u8],
    iq: &[u8],
    dq: &[u8],
    gcp: &[u8],
    mode: XCf,
) -> f64 {
    let rn = read.len();
    let hn = hap.len();
    let cols = hn + 1;
    let init_del = INITIAL_CONDITION / hn as f64;
    let mut m = vec![0.0f64; (rn + 1) * cols];
    let mut ins = vec![0.0f64; (rn + 1) * cols];
    let mut del = vec![0.0f64; (rn + 1) * cols];
    for j in 0..=hn {
        del[j] = init_del;
    }
    let mut gkl_last_x = vec![0.0f32; hn + 1];
    for i in 1..=rn {
        let t = rust_trans(iq[i - 1], dq[i - 1], gcp[i - 1]);
        let x = read[i - 1];
        let q = quals[i - 1];
        let row = i * cols;
        let prev = (i - 1) * cols;
        m[row] = 0.0;
        ins[row] = 0.0;
        del[row] = 0.0;
        for j in 1..=hn {
            let y = hap[j - 1];
            let equal = bases_equal(x, y);
            let p = rust_distm(q, equal);
            let m_diag = m[prev + j - 1];
            let x_diag = ins[prev + j - 1];
            let y_diag = del[prev + j - 1];
            m[row + j] = p * (m_diag * t[0] + x_diag * t[1] + y_diag * t[1]);
            if i == rn {
                gkl_last_x[j] = gkl_closed_x_f32(m[prev + j], ins[prev + j], iq[i - 1], gcp[i - 1]);
            }
            ins[row + j] = m[prev + j] * t[2] + ins[prev + j] * t[3];
            del[row + j] = m[row + j - 1] * t[4] + del[row + j - 1] * t[5];
        }
    }
    let end = rn * cols;
    let mut sum = 0.0;
    let mut sum_x = 0.0;
    let mut gkl_acc = 0.0f32;
    let mut injected = false;
    for j in 1..=hn {
        let m_j = m[end + j];
        let x_j = ins[end + j];
        let gkl_x = gkl_last_x[j];
        let gkl_before = gkl_acc;
        let rust_before = sum_x;
        let mut term = m_j + x_j;
        if mode == XCf::FirstBothNonzeroSumxAdd && !injected && gkl_before != 0.0 && gkl_x != 0.0 {
            injected = true;
            let gkl_add = gkl_before + gkl_x;
            let rescaled = INITIAL_CONDITION * f64::from(gkl_add);
            term += rescaled - (rust_before + x_j);
            sum_x = rescaled;
        } else {
            sum_x = rust_before + x_j;
            gkl_acc = gkl_before + gkl_x;
        }
        sum += term;
    }
    if sum <= 0.0 || !sum.is_finite() {
        return f64::NEG_INFINITY;
    }
    sum.log10() - INITIAL_CONDITION_LOG10
}

fn collect_acc(
    read: &[u8],
    quals: &[u8],
    hap: &[u8],
    iq: &[u8],
    dq: &[u8],
    gcp: &[u8],
    retain_row: usize,
    hap_index: usize,
    first_both: &mut Option<AccCell>,
) {
    if first_both.is_some() {
        return;
    }
    let rn = read.len();
    let hn = hap.len();
    let cols = hn + 1;
    let init_del = INITIAL_CONDITION / hn as f64;
    let mut m = vec![0.0f64; (rn + 1) * cols];
    let mut ins = vec![0.0f64; (rn + 1) * cols];
    let mut del = vec![0.0f64; (rn + 1) * cols];
    for j in 0..=hn {
        del[j] = init_del;
    }
    let mut gkl_last_x = vec![0.0f32; hn + 1];
    let mut m_up_last = vec![0.0f64; hn + 1];
    let mut x_up_last = vec![0.0f64; hn + 1];
    for i in 1..=rn {
        let t = rust_trans(iq[i - 1], dq[i - 1], gcp[i - 1]);
        let x = read[i - 1];
        let q = quals[i - 1];
        let row = i * cols;
        let prev = (i - 1) * cols;
        m[row] = 0.0;
        ins[row] = 0.0;
        del[row] = 0.0;
        for j in 1..=hn {
            let y = hap[j - 1];
            let equal = bases_equal(x, y);
            let p = rust_distm(q, equal);
            let m_diag = m[prev + j - 1];
            let x_diag = ins[prev + j - 1];
            let y_diag = del[prev + j - 1];
            m[row + j] = p * (m_diag * t[0] + x_diag * t[1] + y_diag * t[1]);
            if i == rn {
                m_up_last[j] = m[prev + j];
                x_up_last[j] = ins[prev + j];
                gkl_last_x[j] = gkl_closed_x_f32(m_up_last[j], x_up_last[j], iq[i - 1], gcp[i - 1]);
            }
            ins[row + j] = m[prev + j] * t[2] + ins[prev + j] * t[3];
            del[row + j] = m[row + j - 1] * t[4] + del[row + j - 1] * t[5];
        }
    }
    assert_eq!(quals.len(), rn);
    let end = rn * cols;
    let mut gkl_acc = 0.0f32;
    let mut rust_acc = 0.0f64;
    let mut prior_nonzero_j = 0usize;
    let mut prior_gkl_x = 0.0f32;
    for j in 1..=hn {
        let gkl_x = gkl_last_x[j];
        let rust_x = ins[end + j];
        if gkl_acc != 0.0 && gkl_x != 0.0 {
            *first_both = Some(AccCell {
                retain_row,
                hap_index,
                read_len: rn,
                hap_len: hn,
                hap_j: j,
                prior_nonzero_j,
                prior_gkl_x,
                ins_q: iq[rn - 1],
                gcp: gcp[rn - 1],
                m_up: m_up_last[j],
                x_up: x_up_last[j],
                gkl_sum_before: gkl_acc,
                gkl_x,
                rust_sum_before: rust_acc,
                rust_x,
            });
            return;
        }
        if gkl_x != 0.0 && prior_nonzero_j == 0 {
            prior_nonzero_j = j;
            prior_gkl_x = gkl_x;
        }
        gkl_acc += gkl_x;
        rust_acc += rust_x;
    }
}

fn dump_acc(cell: &AccCell) {
    let gkl_mx = gkl_mx_f32(cell.m_up, cell.ins_q);
    let gkl_xx = gkl_xx_f32(cell.x_up, cell.gcp);
    let gkl_x = gkl_mx + gkl_xx;
    let ieee_add = cell.gkl_sum_before + cell.gkl_x;
    let ftz_add = gkl_add_ftz(cell.gkl_sum_before, cell.gkl_x);
    let gkl_add_w = INITIAL_CONDITION * f64::from(ieee_add);
    let gkl_ops_f64 = INITIAL_CONDITION * (f64::from(cell.gkl_sum_before) + f64::from(cell.gkl_x));
    let rust_add = cell.rust_sum_before + cell.rust_x;
    let new_f32_add = gkl_add_w.to_bits() != gkl_ops_f64.to_bits();
    let sum_before_is_prior_x = cell.gkl_sum_before.to_bits() == cell.prior_gkl_x.to_bits();
    let ftz_changes_store = ieee_add.to_bits() != ftz_add.to_bits();
    kv(
        "acc_cell",
        format!(
            "retain_row={}\thap={}\tread_len={}\thap_len={}\thap_j={}\tprior_nonzero_j={}\tins_q={}\tgcp={}\tremaining_rows={}\tresult_lane={}\tfirst_add_one_zero_operand=true\tgkl_sum_before={}\tgkl_x={}\tgkl_mx_bits=0x{:08x}\tgkl_xx_bits=0x{:08x}\tieee_add={}\tgkl_ftz_stored={}\tieee_add_rescaled={}\trust_sum_before={}\trust_x={}\trust_add={}\tgkl_ops_f64_add={}\tulp_vs_rust={}\tulp_ieee_vs_f64_add_of_gkl_ops={}\tnew_f32_add_rounding={}\tx_is_6r267_add={}\tsum_before_equals_prior_x_bits={}\tsum_before_subnormal={}\tx_subnormal={}\tieee_add_subnormal={}\tftz_changes_stored_value={}\tshift_changes_sumx=false\tftz_mode_set=true\tdaz_enabled=false\tcf1_value=ieee_add",
            cell.retain_row,
            cell.hap_index,
            cell.read_len,
            cell.hap_len,
            cell.hap_j,
            cell.prior_nonzero_j,
            cell.ins_q,
            cell.gcp,
            if cell.read_len % 8 == 0 { 8 } else { cell.read_len % 8 },
            (cell.read_len - 1) % 8,
            fmt_f32(cell.gkl_sum_before),
            fmt_f32(cell.gkl_x),
            gkl_mx.to_bits(),
            gkl_xx.to_bits(),
            fmt_f32(ieee_add),
            fmt_f32(ftz_add),
            fmt_f64(gkl_add_w),
            fmt_f64(cell.rust_sum_before),
            fmt_f64(cell.rust_x),
            fmt_f64(rust_add),
            fmt_f64(gkl_ops_f64),
            ulp_f64(gkl_add_w, rust_add),
            ulp_f64(gkl_add_w, gkl_ops_f64),
            new_f32_add,
            gkl_x.to_bits() == cell.gkl_x.to_bits(),
            sum_before_is_prior_x,
            cell.gkl_sum_before.is_subnormal(),
            cell.gkl_x.is_subnormal(),
            ieee_add.is_subnormal(),
            ftz_changes_store,
        ),
    );
}

#[test]
fn forensic_6r278_gkl_avx_last_stripe_sumx_add() {
    kv("java_pin", JAVA_PIN);
    kv("gkl_pin", GKL_PIN);
    kv(
        "gkl_initial_f32_bits",
        format!("0x{:08x}", GKL_INITIAL_F32.to_bits()),
    );
    assert_eq!(gkl_ph2pr(20).to_bits(), Q20_GKL_PH2PR_F32);
    assert_eq!(gkl_match(20).to_bits(), Q20_GKL_MATCH_F32);
    assert_eq!(gkl_mismatch(20).to_bits(), Q20_GKL_MISMATCH_F32);
    kv("production_change", "NONE");
    kv("target", "20:29455649 T/TGTTTG");
    kv(
        "frozen_input",
        "6R.257: Java174-mers + Java clip 20:29455560-29455728; 123 retainEvidence",
    );
    kv(
        "gkl_x_op",
        "compute_full_prob last stripe: sumX = VEC_ADD(sumX, X_t.d) at avx-pairhmm-template.h:357; f32 add; not FMA; NUMBER=float; X_t is the completed post-computeMXY value from line 219",
    );
    kv(
        "gkl_source",
        "IntelLabs/GKL 0.8.8 src/main/native/pairhmm/avx-pairhmm-template.h:357 sumX = VEC_ADD(sumX, X_t.d); avx-functions-float.h:96-97 VEC_ADD=_mm256_add_ps; avx512-functions-float.h:120-121 VEC_ADD=_mm512_add_ps; sumX is zero-initialized by VEC_SET1_VAL(zero) at line 336; result lane is sumMX.f[remainingRows-1] after line 369",
    );
    kv(
        "oracle_vs_host",
        "Java oracle is x86_64 (AVX or AVX-512 if supported, not Apple); this host is aarch64 Darwin; native GKL AVX is not executed",
    );
    kv(
        "first_cell_init",
        "sumX starts at +0. The first add on the extracted last-row lane is 0+X and has one zero operand. The first meaningful add is the first last-row column where both the f32 sumX and the completed f32 X_t are nonzero. Only lane (read_len-1)%8 is returned",
    );
    kv(
        "vector_shift",
        "Line 355 _vector_shift_last(M_t) runs before this add and shifts M_t, not X_t and not sumX. X_t was already produced by computeMXY line 219. Line 358 shifts X_t after the add. avx-vector-shift.h:56-78 is a byte-lane move. Neither shift changes the arithmetic value entering sumX = VEC_ADD(sumX, X_t.d)",
    );
    kv(
        "rust_formula",
        "final_sum += m[rn,j] + ins[rn,j] in f64; this arrow adjusts only the X running sum at the first both-nonzero last-row add. The M fold and the M/X/Y recurrence stay Rust",
    );
    kv(
        "rust_source",
        "gatk-haplotypecaller/src/pairhmm_logless.rs 358-361 f64 fold of (m+ins); pairhmm recurrence itself is unchanged",
    );
    kv(
        "gkl_primitive_kind",
        "f32 add: VEC_ADD=_mm256_add_ps (avx-functions-float.h:96-97) and _mm512_add_ps (avx512-functions-float.h:120-121); not FMA",
    );
    kv(
        "counterfactual_rule",
        "CF1 first-both-nonzero-sumX-add only: on the last read row, the first column with GKL sumX!=0 and GKL X_t!=0 replaces only that add with INITIAL*f32(gkl_sumX_before + gkl_closed_x_6r267). sumM, recurrence, probabilities, geometry, normalization, and later adds stay Rust. CF1 stores the IEEE f32 sum; FTZ is recorded separately",
    );
    assert_eq!(DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, 128);

    let java_haps: Vec<Vec<u8>> = JAVA_FNV.iter().map(|h| java_seq(h)).collect();
    for (i, seq) in java_haps.iter().enumerate() {
        assert_eq!(fnv1a64_hex(seq), JAVA_FNV[i]);
        assert_eq!(seq.len(), 174);
        assert_eq!(&seq[..9], JAVA_PREFIX);
        assert_eq!(&seq[seq.len() - 4..], JAVA_SUFFIX);
    }

    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing chr20_tiny BAM/REF");
        return;
    }
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, INTERVAL).expect("interval");
    begin_hap_list_observe();
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
    let region = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= TARGET
                && r.end.get() >= TARGET
        })
        .expect("ActiveFull");
    let outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let trim = take_hap_list_trim_span().expect("trim");
    assert_eq!((trim.trim_start, trim.trim_end), RUST_CLIP);
    let snap = take_colocated_merge_numerics()
        .into_iter()
        .find(|s| s.loc == TARGET)
        .expect("snap");
    assert_eq!(snap.n_reads, 123);
    assert_eq!(TARGET_REF, "T");
    assert_eq!(TARGET_ALT, "TGTTTG");

    let rust_haps: Vec<Vec<u8>> = FROZEN_IDX
        .iter()
        .map(|&idx| outcome.assembly.haplotypes[idx].bases.clone())
        .collect();
    for (i, (java, rust)) in java_haps.iter().zip(rust_haps.iter()).enumerate() {
        assert_eq!(fnv1a64_hex(rust), RUST_FNV[i]);
        assert_eq!(&java[9..java.len() - 4], rust.as_slice());
    }

    let cfg = HcLikelihoodEngineConfig::gatk_haplotype_caller_production();
    assert_eq!(cfg.resolved_pair_hmm_backend(), PairHmmBackend::NeonF64);
    kv(
        "pairhmm_backend",
        format!(
            "configured={} resolved={}",
            cfg.primary_engine_label(),
            cfg.resolved_pair_hmm_backend().label()
        ),
    );

    let retain: std::collections::HashSet<usize> = snap.ad_row_read_index.iter().copied().collect();
    let subset: Vec<_> = outcome
        .read_likelihoods
        .iter()
        .filter(|c| retain.contains(&c.read_index.get()))
        .cloned()
        .collect();
    let hap_rows = region_likelihoods_to_rows(&subset, outcome.assembly.haplotypes.len());
    let mut by_read = std::collections::HashMap::new();
    for row in &hap_rows {
        by_read.entry(row.read_index).or_insert(row);
    }
    let l0: Vec<f64> = snap.ad_row_lls.iter().map(|ll| ll[0]).collect();

    let finalized = finalize_region_reads_for_assembly(
        &region.reads,
        region,
        true,
        gatk_min_tail_quality_for_assembly(10),
        false,
    );
    let java_map = clip_map(&finalized, region, JAVA_CLIP.0, JAVA_CLIP.1);
    let java_refs: Vec<&[u8]> = java_haps.iter().map(|s| s.as_slice()).collect();

    let mut unique_q = BTreeSet::new();
    let mut n_bq_total = 0usize;
    let mut rust_pooled = vec![0.0f64; 123];
    let mut cf1_pooled = vec![0.0f64; 123];
    let mut first_both: Option<AccCell> = None;
    let mut n_cf1_hap_diff = 0u64;
    let mut first_score_dumped = false;

    for (ri, &read_idx) in snap.ad_row_read_index.iter().enumerate() {
        let key = (
            snap.ad_row_qname[ri].as_bytes().to_vec(),
            snap.ad_row_flags[ri],
        );
        let rec = java_map
            .get(&key)
            .unwrap_or_else(|| panic!("missing Java-clip read row={ri}"));
        let p = planes(rec, &cfg);
        n_bq_total += p.bq.len();
        unique_q.extend(p.bq.iter().copied());

        for k in 0..5 {
            collect_acc(
                &p.bases,
                &p.bq,
                &java_haps[k],
                &p.iq,
                &p.dq,
                &p.gcp,
                ri,
                k,
                &mut first_both,
            );
        }

        let prod_row = by_read.get(&read_idx).expect("prod");
        let p5: [f64; 5] =
            std::array::from_fn(|k| prod_row.haplotype_log10_likelihoods[FROZEN_IDX[k]]);
        if p5.iter().all(|&v| v == 0.0) {
            continue;
        }
        let other_best = prod_row
            .haplotype_log10_likelihoods
            .iter()
            .enumerate()
            .filter(|(i, _)| !FROZEN_IDX.contains(i))
            .map(|(_, v)| *v)
            .fold(f64::NEG_INFINITY, f64::max);

        let scores = score_read_against_haplotypes(
            &cfg,
            &p.bases,
            rec.qual(),
            rec.mapq(),
            &java_refs,
            bam_indel_phred(rec, b"BI").as_deref(),
            bam_indel_phred(rec, b"BD").as_deref(),
        )
        .expect("neon");
        let raw: [f64; 5] = std::array::from_fn(|k| scores[k]);
        rust_pooled[ri] = floor5(raw, other_best)
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);

        let mut cf1_raw = [0.0f64; 5];
        for k in 0..5 {
            cf1_raw[k] = rust_diag_ll(
                &p.bases,
                &p.bq,
                &java_haps[k],
                &p.iq,
                &p.dq,
                &p.gcp,
                XCf::FirstBothNonzeroSumxAdd,
            );
            let scalar =
                logless_pairhmm_likelihood(&p.bases, &p.bq, &java_haps[k], &p.iq, &p.dq, &p.gcp)
                    .expect("scalar");
            assert_eq!(
                scalar.to_bits(),
                raw[k].to_bits(),
                "NEON must remain bit-identical to scalar (6R.258)"
            );
            let baseline = rust_diag_ll(
                &p.bases,
                &p.bq,
                &java_haps[k],
                &p.iq,
                &p.dq,
                &p.gcp,
                XCf::Baseline,
            );
            assert_eq!(
                baseline.to_bits(),
                raw[k].to_bits(),
                "diagnostic baseline must match NEON/scalar"
            );
            if cf1_raw[k].to_bits() != raw[k].to_bits() {
                n_cf1_hap_diff += 1;
            }
            if !first_score_dumped && k == 0 && ri == 0 {
                first_score_dumped = true;
                kv(
                    "first_matrix_ll",
                    format!(
                        "baseline={} cf1={} bits_equal={}",
                        fmt_f64(raw[k]),
                        fmt_f64(cf1_raw[k]),
                        raw[k].to_bits() == cf1_raw[k].to_bits(),
                    ),
                );
            }
        }
        cf1_pooled[ri] = floor5(cf1_raw, other_best)
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
    }

    kv(
        "unique_bq_values",
        unique_q
            .iter()
            .map(|q| q.to_string())
            .collect::<Vec<_>>()
            .join(","),
    );
    kv("n_unique_bq", unique_q.len().to_string());
    kv("n_bq_total", n_bq_total.to_string());
    assert!(unique_q.contains(&20));

    let cell = first_both.expect("first both-nonzero last-stripe sumX add");
    assert_ne!(
        cell.gkl_sum_before, 0.0,
        "sumX before the add must be nonzero"
    );
    assert_ne!(cell.gkl_x, 0.0, "X_t must be nonzero");
    assert_eq!(
        cell.gkl_sum_before.to_bits(),
        cell.prior_gkl_x.to_bits(),
        "prior adds have a zero operand, so sumX equals the single prior nonzero X_t"
    );
    assert!(cell.hap_j > cell.prior_nonzero_j);
    assert!(
        cell.prior_nonzero_j > 0,
        "a trivial 0+X add precedes this one"
    );
    let gkl_mx = gkl_mx_f32(cell.m_up, cell.ins_q);
    let gkl_xx = gkl_xx_f32(cell.x_up, cell.gcp);
    assert_eq!(
        (gkl_mx + gkl_xx).to_bits(),
        cell.gkl_x.to_bits(),
        "X_t must be the 6R.267 f32 add"
    );
    kv(
        "first_x_primitive",
        "sumX = VEC_ADD(sumX, X_t.d) = _mm256_add_ps / _mm512_add_ps",
    );
    kv("gkl_op_kind", "f32 add; not FMA");
    kv(
        "first_both_accumulation",
        format!(
            "retain_row={} hap={} read_len={} hap_len={} hap_j={} prior_nonzero_j={} ins_q={} gcp={} remaining_rows={} result_lane={} first_add_has_one_zero_operand=true",
            cell.retain_row,
            cell.hap_index,
            cell.read_len,
            cell.hap_len,
            cell.hap_j,
            cell.prior_nonzero_j,
            cell.ins_q,
            cell.gcp,
            if cell.read_len % 8 == 0 { 8 } else { cell.read_len % 8 },
            (cell.read_len - 1) % 8,
        ),
    );
    dump_acc(&cell);

    kv(
        "counterfactual_note",
        "CF1 injects only the first both-nonzero last-row VEC_ADD(sumX, X_t) as INITIAL*f32(sumX_before + X_t). sumM and the M/X/Y recurrence stay Rust f64. CF1 uses the IEEE add; the FTZ store is reported separately",
    );

    let rust_gls = biallelic_gls(&l0, &rust_pooled);
    let (rr_rel, rr_cont, rr_pl) = emitted_hom_alt(&rust_gls);
    kv("baseline_emitted_GL", fmt_f64(rr_rel));
    kv("baseline_continuous_PL", fmt_f64(rr_cont));
    kv("baseline_integer_PL", rr_pl.to_string());
    assert!((rr_rel - JOINT_GL).abs() < 1e-9);
    assert_eq!(rr_pl, 3519);

    let cf1_gls = biallelic_gls(&l0, &cf1_pooled);
    let (cf1_rel, cf1_cont, cf1_pl) = emitted_hom_alt(&cf1_gls);
    kv("first_divergent_operation", "GKL_AVX_LAST_STRIPE_SUMX_ADD");
    kv("cf1_label", "first both-nonzero last-stripe sumX f32 add");
    kv("cf1_emitted_GL", fmt_f64(cf1_rel));
    kv("cf1_continuous_PL", fmt_f64(cf1_cont));
    kv("cf1_integer_PL", cf1_pl.to_string());
    kv("cf1_delta_GL", fmt_f64(cf1_rel - rr_rel));
    kv("cf1_delta_continuous_PL", fmt_f64(cf1_cont - rr_cont));
    kv("cf1_delta_integer_PL", (cf1_pl - rr_pl).to_string());
    kv("n_cf1_hap_score_diff", n_cf1_hap_diff.to_string());
    kv("required_pl_movement", format!("{REQUIRED_MOVE:.17}"));
    kv("cf1_crosses_3517_5", (cf1_cont < BOUNDARY).to_string());
    let toward_java = rr_cont - cf1_cont;
    kv("cf1_movement_toward_java", fmt_f64(toward_java));

    let ieee_add = cell.gkl_sum_before + cell.gkl_x;
    let ftz_add = gkl_add_ftz(cell.gkl_sum_before, cell.gkl_x);
    let gkl_add = ieee_add;
    let isolated_add = ieee_add;
    assert_eq!(isolated_add.to_bits(), gkl_add.to_bits());
    let gkl_add_w = INITIAL_CONDITION * f64::from(gkl_add);
    let gkl_ops_f64 = INITIAL_CONDITION * (f64::from(cell.gkl_sum_before) + f64::from(cell.gkl_x));
    let rust_add = cell.rust_sum_before + cell.rust_x;
    let gkl_sum_w = INITIAL_CONDITION * f64::from(cell.gkl_sum_before);
    let gkl_x_w = INITIAL_CONDITION * f64::from(cell.gkl_x);
    kv(
        "selected_add_bits",
        format!(
            "gkl_sum_before_f32=0x{:08x} gkl_x_f32=0x{:08x} ieee_add_f32=0x{:08x} ftz_stored_f32=0x{:08x} cf1_f32=0x{:08x} gkl_sum_before_rescaled={} gkl_x_rescaled={} gkl_add_rescaled={} rust_sum_before={} rust_x={} rust_add={} gkl_ops_f64_add={} isolated_matches={} new_f32_add_rounding={} gkl_is_f32_add_not_fma=true ftz_mode_set=true daz_enabled=false cf1_applies_ftz=false ftz_changes_stored_value={} sum_before_subnormal={} x_subnormal={} ieee_add_subnormal={} shift_changes_arithmetic=false",
            cell.gkl_sum_before.to_bits(),
            cell.gkl_x.to_bits(),
            ieee_add.to_bits(),
            ftz_add.to_bits(),
            isolated_add.to_bits(),
            fmt_f64(gkl_sum_w),
            fmt_f64(gkl_x_w),
            fmt_f64(gkl_add_w),
            fmt_f64(cell.rust_sum_before),
            fmt_f64(cell.rust_x),
            fmt_f64(rust_add),
            fmt_f64(gkl_ops_f64),
            isolated_add.to_bits() == gkl_add.to_bits(),
            gkl_add_w.to_bits() != gkl_ops_f64.to_bits(),
            ieee_add.to_bits() != ftz_add.to_bits(),
            cell.gkl_sum_before.is_subnormal(),
            cell.gkl_x.is_subnormal(),
            ieee_add.is_subnormal(),
        ),
    );
    kv(
        "inherited_divergence",
        format!(
            "sumX_ulp_vs_rust={} x_t_ulp_vs_rust={} (X_t is the 6R.267 closed f32 add of Rust M_t_1 and X_t_1; not new 6R.278)",
            ulp_f64(gkl_sum_w, cell.rust_sum_before),
            ulp_f64(gkl_x_w, cell.rust_x),
        ),
    );
    kv(
        "new_6r278_rounding",
        format!(
            "ieee_f32_add_vs_f64_add_of_same_gkl_ops ulp={} bits_differ={}",
            ulp_f64(gkl_add_w, gkl_ops_f64),
            gkl_add_w.to_bits() != gkl_ops_f64.to_bits(),
        ),
    );

    let classification = if toward_java > REQUIRED_MOVE {
        "PAIRHMM_LAST_STRIPE_SUMX_ADD_DIVERGENCE_CAUSAL"
    } else {
        "PAIRHMM_LAST_STRIPE_SUMX_ADD_DIVERGENCE_NOT_CAUSAL"
    };
    kv("classification", classification);
    kv("production_change", "NONE");
    kv(
        "next_arrow",
        if classification == "PAIRHMM_LAST_STRIPE_SUMX_ADD_DIVERGENCE_CAUSAL" {
            "STOP — causal arrow identified; no production patch applied yet"
        } else {
            "GKL_AVX_LAST_STRIPE_SUMMX_ADD: the next arithmetic line is avx-pairhmm-template.h:368 sumMX.d = VEC_ADD(sumM, sumX) (_mm256_add_ps). Lines 358 and 360 are _vector_shift_last of X_t and Y_t_1 and are byte-lane moves; there is no sumY accumulation; diagnostic-only"
        },
    );
    assert!(
        classification == "PAIRHMM_LAST_STRIPE_SUMX_ADD_DIVERGENCE_CAUSAL"
            || classification == "PAIRHMM_LAST_STRIPE_SUMX_ADD_DIVERGENCE_NOT_CAUSAL"
    );
}
