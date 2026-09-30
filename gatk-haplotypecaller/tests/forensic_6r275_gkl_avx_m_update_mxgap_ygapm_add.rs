//! 6R.275: measurement-only. GKL AVX `computeMXY` inner second M add
//! `VEC_ADD(VEC_ADD(M*pMM, X*pGAPM), VEC_MUL(Y_t_2, pGAPM))` vs Rust f64
//! `(m_diag*t_m2m + ins_diag*t_i2m) + del_diag*t_i2m` at
//! `20:29455649 T/TGTTTG`. PRODUCTION CHANGE: NONE.
//!
//! 6R.274 closed `(M*pMM)+(X*pGAPM)`; 6R.273 closed `Y*pGAPM`. This round
//! injects only that final inner `VEC_ADD`. Native AVX GKL is not executed
//! (oracle x86_64; this host aarch64).
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r275_gkl_avx_m_update_mxgap_ygapm_add -- --nocapture --test-threads=1
//! HOLDOUT_6R275=1 cargo test -p gatk-haplotypecaller --test holdout_6r275_gkl_avx_m_update_mxgap_ygapm_add -- --nocapture --test-threads=1
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
const GKL_PIN: &str =
    "0.8.8 IntelLabs/GKL avx-pairhmm-template.h:213 VEC_ADD(VEC_ADD(M*pMM,X*pGAPM), VEC_MUL(Y_t_2,pGAPM))";
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
    eprintln!("6R275\t{key}\t{}", value.as_ref());
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

fn classify_mismatch(q: u8) -> &'static str {
    let gkl_w = f64::from(gkl_mismatch(q));
    let rust = rust_mismatch(q);
    let f64_div = f64::from(gkl_ph2pr(q)) / 3.0;
    if gkl_w.to_bits() == rust.to_bits() {
        "exact-match"
    } else if gkl_w.to_bits() != f64_div.to_bits() {
        "new-f32-division-rounding"
    } else {
        "inherited-only"
    }
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
enum MCf {
    Baseline,
    FirstBothMxgapYgapmAdd,
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

/// IEEE-754 f32 `VEC_ADD(partial_M, Y*pGAPM)`, rescaled to Rust INITIAL.
fn gkl_mxgap_ygapm_add_rescaled(
    m_diag: f64,
    x_diag: f64,
    y_diag: f64,
    ins_q: u8,
    del_q: u8,
    gcp: u8,
) -> f64 {
    INITIAL_CONDITION
        * f64::from(
            gkl_partial_m_f32(m_diag, x_diag, ins_q, del_q, gcp) + gkl_ygapm_f32(y_diag, gcp),
        )
}

fn ulp_f64(a: f64, b: f64) -> u64 {
    a.to_bits().abs_diff(b.to_bits())
}

#[derive(Clone, Copy)]
struct MCell {
    category: &'static str,
    retain_row: usize,
    hap_index: usize,
    read_i: usize,
    hap_j: usize,
    bq: u8,
    ins_q: u8,
    del_q: u8,
    gcp: u8,
    m_diag: f64,
    x_diag: f64,
    y_diag: f64,
}

fn rust_diag_ll(
    read: &[u8],
    quals: &[u8],
    hap: &[u8],
    iq: &[u8],
    dq: &[u8],
    gcp: &[u8],
    mode: MCf,
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
    let mut injected = false;
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
            let p = rust_distm(q, bases_equal(x, y));
            let m_diag = m[prev + j - 1];
            let x_diag = ins[prev + j - 1];
            let y_diag = del[prev + j - 1];
            let rust_partial = m_diag * t[0] + x_diag * t[1];
            let rust_yg = y_diag * t[1];
            let mx = if mode == MCf::FirstBothMxgapYgapmAdd
                && !injected
                && rust_partial != 0.0
                && rust_yg != 0.0
            {
                injected = true;
                gkl_mxgap_ygapm_add_rescaled(
                    m_diag,
                    x_diag,
                    y_diag,
                    iq[i - 1],
                    dq[i - 1],
                    gcp[i - 1],
                )
            } else {
                rust_partial + rust_yg
            };
            m[row + j] = p * mx;
            ins[row + j] = m[prev + j] * t[2] + ins[prev + j] * t[3];
            del[row + j] = m[row + j - 1] * t[4] + del[row + j - 1] * t[5];
        }
    }
    let end = rn * cols;
    let mut sum = 0.0;
    for j in 1..=hn {
        sum += m[end + j] + ins[end + j];
    }
    if sum <= 0.0 || !sum.is_finite() {
        return f64::NEG_INFINITY;
    }
    sum.log10() - INITIAL_CONDITION_LOG10
}

fn collect_m_cells(
    read: &[u8],
    quals: &[u8],
    hap: &[u8],
    iq: &[u8],
    dq: &[u8],
    gcp: &[u8],
    retain_row: usize,
    hap_index: usize,
    first_i1j1: &mut Option<MCell>,
    first_i2j2: &mut Option<MCell>,
    first_i3j3: &mut Option<MCell>,
    first_both: &mut Option<MCell>,
) {
    if first_i1j1.is_some() && first_i2j2.is_some() && first_i3j3.is_some() && first_both.is_some()
    {
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
            let p = rust_distm(q, bases_equal(x, y));
            let m_diag = m[prev + j - 1];
            let x_diag = ins[prev + j - 1];
            let y_diag = del[prev + j - 1];
            let gkl_partial = gkl_partial_m_f32(m_diag, x_diag, iq[i - 1], dq[i - 1], gcp[i - 1]);
            let gkl_yg = gkl_ygapm_f32(y_diag, gcp[i - 1]);
            let cell = MCell {
                category: "",
                retain_row,
                hap_index,
                read_i: i,
                hap_j: j,
                bq: q,
                ins_q: iq[i - 1],
                del_q: dq[i - 1],
                gcp: gcp[i - 1],
                m_diag,
                x_diag,
                y_diag,
            };
            if i == 1 && j == 1 && first_i1j1.is_none() {
                *first_i1j1 = Some(MCell {
                    category: "mxgap_ygap_i1j1",
                    ..cell
                });
            }
            if i == 2 && j == 2 && first_i2j2.is_none() {
                *first_i2j2 = Some(MCell {
                    category: "sixr271_i2j2",
                    ..cell
                });
            }
            if i == 3 && j == 3 && first_i3j3.is_none() {
                *first_i3j3 = Some(MCell {
                    category: "sixr274_i3j3",
                    ..cell
                });
            }
            if gkl_partial != 0.0 && gkl_yg != 0.0 && first_both.is_none() {
                *first_both = Some(MCell {
                    category: "first_both_nonzero_mxgap_ygapm_add",
                    ..cell
                });
            }
            m[row + j] = p * (m_diag * t[0] + x_diag * t[1] + y_diag * t[1]);
            ins[row + j] = m[prev + j] * t[2] + ins[prev + j] * t[3];
            del[row + j] = m[row + j - 1] * t[4] + del[row + j - 1] * t[5];
            if first_i1j1.is_some()
                && first_i2j2.is_some()
                && first_i3j3.is_some()
                && first_both.is_some()
            {
                return;
            }
        }
    }
}

fn dump_m_cell(cell: &MCell) {
    let rust_mm = cell.m_diag * rust_m2m(cell.ins_q, cell.del_q);
    let rust_xg = cell.x_diag * rust_pgapm(cell.gcp);
    let rust_partial = rust_mm + rust_xg;
    let rust_yg = cell.y_diag * rust_pgapm(cell.gcp);
    let rust_add = rust_partial + rust_yg;
    let gkl_mm = gkl_mm_f32(cell.m_diag, cell.ins_q, cell.del_q);
    let gkl_xg = gkl_xgapm_f32(cell.x_diag, cell.gcp);
    let gkl_partial = gkl_mm + gkl_xg;
    let gkl_yg = gkl_ygapm_f32(cell.y_diag, cell.gcp);
    let gkl_add = gkl_partial + gkl_yg;
    let isolated_add = gkl_partial + gkl_yg;
    let gkl_add_w = INITIAL_CONDITION * f64::from(gkl_add);
    let gkl_ops_f64_add = INITIAL_CONDITION * (f64::from(gkl_partial) + f64::from(gkl_yg));
    let gkl_partial_w = INITIAL_CONDITION * f64::from(gkl_partial);
    let gkl_yg_w = INITIAL_CONDITION * f64::from(gkl_yg);
    let abs = (gkl_add_w - rust_add).abs();
    let rel = if rust_add != 0.0 {
        abs / rust_add.abs()
    } else {
        0.0
    };
    let new_f32_add_rounding = gkl_add_w.to_bits() != gkl_ops_f64_add.to_bits();
    let y_below_one_ulp = gkl_add.to_bits() == gkl_partial.to_bits();
    kv(
        "m_cell",
        format!(
            "cat={}\tretain_row={}\thap={}\tread_i={}\thap_j={}\tbq={}\tins_q={}\tdel_q={}\tgcp={}\tm_diag_ratio={}\tx_diag_ratio={}\ty_diag_ratio={}\tgkl_mm_f32_bits=0x{:08x}\tgkl_xg_f32_bits=0x{:08x}\tgkl_partial_m_f32_bits=0x{:08x}\tgkl_yg_f32_bits=0x{:08x}\tgkl_add_f32_bits=0x{:08x}\tgkl_partial_rescaled={}\tgkl_yg_rescaled={}\tgkl_add_rescaled={}\trust_partial={}\trust_yg={}\trust_add={}\tgkl_ops_f64_add={}\tabs_delta={:.17e}\trel_delta={:.17e}\tulp_vs_rust={}\tulp_f32_vs_f64_gkl_ops={}\tnew_f32_add_rounding={}\tisolated_matches_gkl={}\ty_below_one_ulp_of_partial={}\tadd_trivial={}",
            cell.category,
            cell.retain_row,
            cell.hap_index,
            cell.read_i,
            cell.hap_j,
            cell.bq,
            cell.ins_q,
            cell.del_q,
            cell.gcp,
            fmt_f64(cell.m_diag / INITIAL_CONDITION),
            fmt_f64(cell.x_diag / INITIAL_CONDITION),
            fmt_f64(cell.y_diag / INITIAL_CONDITION),
            gkl_mm.to_bits(),
            gkl_xg.to_bits(),
            gkl_partial.to_bits(),
            gkl_yg.to_bits(),
            gkl_add.to_bits(),
            fmt_f64(gkl_partial_w),
            fmt_f64(gkl_yg_w),
            fmt_f64(gkl_add_w),
            fmt_f64(rust_partial),
            fmt_f64(rust_yg),
            fmt_f64(rust_add),
            fmt_f64(gkl_ops_f64_add),
            abs,
            rel,
            ulp_f64(gkl_add_w, rust_add),
            ulp_f64(gkl_add_w, gkl_ops_f64_add),
            new_f32_add_rounding,
            isolated_add.to_bits() == gkl_add.to_bits(),
            y_below_one_ulp,
            rust_partial == 0.0 || rust_yg == 0.0,
        ),
    );
}

#[test]
fn forensic_6r275_gkl_avx_m_update_mxgap_ygapm_add() {
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
        "gkl_m_op",
        "computeMXY: M_t.d=VEC_MUL(VEC_ADD(VEC_ADD(VEC_MUL(M_t_2.d, pMM), VEC_MUL(X_t_2.d, pGAPM)), VEC_MUL(Y_t_2.d, pGAPM)), distmSel); this arrow is inner VEC_ADD(partial_M, Y*pGAPM); not FMA; NUMBER=float",
    );
    kv(
        "gkl_source",
        "IntelLabs/GKL 0.8.8 src/main/native/pairhmm/avx-pairhmm-template.h computeMXY line 213 VEC_ADD(VEC_ADD(M*pMM, X*pGAPM), VEC_MUL(Y_t_2.d, pGAPM)); avx-functions-float.h:96-97 VEC_ADD=_mm256_add_ps; avx512-functions-float.h:120-121 VEC_ADD=_mm512_add_ps; left is 6R.274 f32 partial_M; right is 6R.273 f32 Y*pGAPM; evaluation order is (MM+XGAPM)+YGAPM",
    );
    kv(
        "oracle_vs_host",
        "Java oracle is x86_64 (AVX or AVX-512 if supported, not Apple); this host is aarch64 Darwin; native GKL AVX is not executed",
    );
    kv(
        "first_cell_init",
        "stripe i==0: M_t_2=X_t_2=0, Y_t_2=init_Y so (1,1) is 0+YGAPM (partial_M=0); 6R.271 (2,2) has MM nonzero and Y_t_2=del[1,1]=0 (trivial Y add); 6R.274 (3,3) has MM+X both nonzero and Y_t_2=del[2,2]=0. First meaningful is first cell with GKL partial_M!=0 and GKL Y*pGAPM!=0",
    );
    kv(
        "rust_formula",
        "m[i,j]=distm*(m[i-1,j-1]*t_m2m + ins[i-1,j-1]*t_i2m + del[i-1,j-1]*t_i2m); this arrow is f64 add of the first two-term partial plus the Y term; production NEON vaddq_f64(vaddq_f64(vmul m, vmul ins), vmul del)",
    );
    kv(
        "rust_source",
        "gatk-haplotypecaller/src/pairhmm_logless.rs 347-350 (m*t_m2m + ins*t_i2m + del*t_i2m); pairhmm_simd/neon.rs 304-307 vaddq_f64 of the MM+XGAPM vaddq_f64 plus vmulq_f64(d_diag, ti2m); same left-assoc order as GKL",
    );
    kv(
        "counterfactual_rule",
        "CF1 first-both-nonzero-only: first cell with rust (MM+XGAPM)!=0 and rust Y*pGAPM!=0 replaces only that final add with INITIAL*f32(gkl_partial_m_f32 + gkl_ygapm_f32); 6R.271-274 products, outer distm mul, and later cells stay Rust",
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
    let mut first_i1j1: Option<MCell> = None;
    let mut first_i2j2: Option<MCell> = None;
    let mut first_i3j3: Option<MCell> = None;
    let mut first_both: Option<MCell> = None;
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
            collect_m_cells(
                &p.bases,
                &p.bq,
                &java_haps[k],
                &p.iq,
                &p.dq,
                &p.gcp,
                ri,
                k,
                &mut first_i1j1,
                &mut first_i2j2,
                &mut first_i3j3,
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
                MCf::FirstBothMxgapYgapmAdd,
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
                MCf::Baseline,
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

    let i1 = first_i1j1.expect("mxgap+ygap i=1 j=1");
    let i2j2 = first_i2j2.expect("6R.271 i=2 j=2");
    let i3j3 = first_i3j3.expect("6R.274 i=3 j=3");
    let both = first_both.expect("first both-nonzero MXGAP+YGAPM add");
    assert_eq!(i1.read_i, 1);
    assert_eq!(i1.hap_j, 1);
    assert_eq!(i1.m_diag, 0.0, "(1,1) M_t_2 must be 0");
    assert_eq!(i1.x_diag, 0.0, "(1,1) X_t_2 must be 0");
    assert_ne!(i1.y_diag, 0.0, "(1,1) Y_t_2 is init_Y");
    assert_eq!((i2j2.read_i, i2j2.hap_j), (2, 2));
    assert_ne!(i2j2.m_diag, 0.0, "6R.271 cell M_t_2 is nonzero");
    assert_eq!(
        i2j2.y_diag, 0.0,
        "6R.271 cell Y_t_2=del[1,1] is still 0; Y add trivial"
    );
    assert_eq!((i3j3.read_i, i3j3.hap_j), (3, 3));
    assert_ne!(i3j3.m_diag, 0.0, "6R.274 cell M_t_2 is nonzero");
    assert_ne!(i3j3.x_diag, 0.0, "6R.274 cell X_t_2 is nonzero");
    assert_eq!(
        i3j3.y_diag, 0.0,
        "6R.274 first both-MM+X cell still has Y_t_2=0"
    );
    let gkl_partial_both =
        gkl_partial_m_f32(both.m_diag, both.x_diag, both.ins_q, both.del_q, both.gcp);
    let gkl_yg_both = gkl_ygapm_f32(both.y_diag, both.gcp);
    assert_ne!(
        gkl_partial_both, 0.0,
        "first-both partial_M must be nonzero"
    );
    assert_ne!(gkl_yg_both, 0.0, "first-both Y*pGAPM must be nonzero");
    kv(
        "first_m_primitive",
        "VEC_ADD(VEC_ADD(M*pMM, X*pGAPM), VEC_MUL(Y_t_2, pGAPM)) = _mm256_add_ps / _mm512_add_ps",
    );
    kv(
        "evaluation_order",
        "(M*pMM + X*pGAPM) + Y*pGAPM; left operand is the 6R.274 f32 result; right operand is the 6R.273 f32 product; not (M*pMM + Y*pGAPM) + X*pGAPM",
    );
    kv(
        "i1j1_cell",
        format!(
            "retain_row={} hap={} i={} j={} m_zero={} x_zero={} y_zero={} add_trivial={}",
            i1.retain_row,
            i1.hap_index,
            i1.read_i,
            i1.hap_j,
            i1.m_diag == 0.0,
            i1.x_diag == 0.0,
            i1.y_diag == 0.0,
            true,
        ),
    );
    kv(
        "sixr271_cell",
        format!(
            "retain_row={} hap={} i={} j={} m_zero={} x_zero={} y_zero={} add_trivial={}",
            i2j2.retain_row,
            i2j2.hap_index,
            i2j2.read_i,
            i2j2.hap_j,
            i2j2.m_diag == 0.0,
            i2j2.x_diag == 0.0,
            i2j2.y_diag == 0.0,
            i2j2.y_diag == 0.0,
        ),
    );
    kv(
        "sixr274_cell",
        format!(
            "retain_row={} hap={} i={} j={} m_zero={} x_zero={} y_zero={} add_trivial={}",
            i3j3.retain_row,
            i3j3.hap_index,
            i3j3.read_i,
            i3j3.hap_j,
            i3j3.m_diag == 0.0,
            i3j3.x_diag == 0.0,
            i3j3.y_diag == 0.0,
            i3j3.y_diag == 0.0,
        ),
    );
    kv(
        "first_both_cell",
        format!(
            "retain_row={} hap={} i={} j={} bq={} ins_q={} del_q={} gcp={}",
            both.retain_row,
            both.hap_index,
            both.read_i,
            both.hap_j,
            both.bq,
            both.ins_q,
            both.del_q,
            both.gcp
        ),
    );
    dump_m_cell(&i1);
    dump_m_cell(&i2j2);
    dump_m_cell(&i3j3);
    dump_m_cell(&both);

    kv(
        "counterfactual_note",
        "CF1 injects only the first both-nonzero VEC_ADD of GKL f32 6R.274 partial_M and 6R.273 Y*pGAPM; 6R.271-274 products themselves, outer distm mul, and later cells stay Rust",
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
    kv(
        "first_divergent_operation",
        "GKL_AVX_M_UPDATE_MXGAP_YGAPM_ADD",
    );
    kv("cf1_label", "first both-nonzero MXGAP+YGAPM add only");
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

    let gkl_mm = gkl_mm_f32(both.m_diag, both.ins_q, both.del_q);
    let gkl_xg = gkl_xgapm_f32(both.x_diag, both.gcp);
    let gkl_partial = gkl_mm + gkl_xg;
    let gkl_yg = gkl_ygapm_f32(both.y_diag, both.gcp);
    let gkl_add = gkl_partial + gkl_yg;
    let isolated_add = gkl_partial + gkl_yg;
    assert_eq!(
        isolated_add.to_bits(),
        gkl_add.to_bits(),
        "isolated f32 add must match GKL bit-for-bit"
    );
    let rust_mm = both.m_diag * rust_m2m(both.ins_q, both.del_q);
    let rust_xg = both.x_diag * rust_pgapm(both.gcp);
    let rust_partial = rust_mm + rust_xg;
    let rust_yg = both.y_diag * rust_pgapm(both.gcp);
    let rust_add = rust_partial + rust_yg;
    let gkl_partial_w = INITIAL_CONDITION * f64::from(gkl_partial);
    let gkl_yg_w = INITIAL_CONDITION * f64::from(gkl_yg);
    let gkl_add_w = INITIAL_CONDITION * f64::from(gkl_add);
    let gkl_ops_f64_add = INITIAL_CONDITION * (f64::from(gkl_partial) + f64::from(gkl_yg));
    kv(
        "selected_add_bits",
        format!(
            "gkl_partial_m_f32=0x{:08x} gkl_yg_f32=0x{:08x} gkl_add_f32=0x{:08x} isolated_f32=0x{:08x} gkl_partial_rescaled={} gkl_yg_rescaled={} gkl_add_rescaled={} rust_partial={} rust_yg={} rust_add={} gkl_ops_f64_add={} isolated_matches={} y_below_one_ulp_of_partial={} new_f32_add_rounding={} equals_6r274_partial={} equals_6r273_yg={}",
            gkl_partial.to_bits(),
            gkl_yg.to_bits(),
            gkl_add.to_bits(),
            isolated_add.to_bits(),
            fmt_f64(gkl_partial_w),
            fmt_f64(gkl_yg_w),
            fmt_f64(gkl_add_w),
            fmt_f64(rust_partial),
            fmt_f64(rust_yg),
            fmt_f64(rust_add),
            fmt_f64(gkl_ops_f64_add),
            isolated_add.to_bits() == gkl_add.to_bits(),
            gkl_add.to_bits() == gkl_partial.to_bits(),
            gkl_add_w.to_bits() != gkl_ops_f64_add.to_bits(),
            gkl_partial.to_bits() == (gkl_mm + gkl_xg).to_bits(),
            gkl_yg.to_bits() == gkl_ygapm_f32(both.y_diag, both.gcp).to_bits(),
        ),
    );
    kv(
        "inherited_divergence",
        format!(
            "partial_M_ulp_vs_rust={} yg_ulp_vs_rust={} (6R.271+272+274 in left; 6R.273 in right; not new 6R.275)",
            ulp_f64(gkl_partial_w, rust_partial),
            ulp_f64(gkl_yg_w, rust_yg),
        ),
    );
    kv(
        "new_6r275_rounding",
        format!(
            "f32_add_vs_f64_add_of_same_gkl_ops ulp={} bits_differ={}",
            ulp_f64(gkl_add_w, gkl_ops_f64_add),
            gkl_add_w.to_bits() != gkl_ops_f64_add.to_bits(),
        ),
    );

    let classification = if toward_java > REQUIRED_MOVE {
        "PAIRHMM_M_UPDATE_MXGAP_YGAPM_ADD_CAUSAL"
    } else {
        "PAIRHMM_M_UPDATE_MXGAP_YGAPM_ADD_DIVERGENCE_NOT_CAUSAL"
    };
    kv("classification", classification);
    kv("production_change", "NONE");
    kv(
        "next_arrow",
        if classification == "PAIRHMM_M_UPDATE_MXGAP_YGAPM_ADD_CAUSAL" {
            "STOP — causal arrow identified; no production patch applied yet"
        } else {
            "GKL_AVX_M_UPDATE_SUM_DISTM_MUL: outer M multiply VEC_MUL(closed_6R275_f32_inner_sum, distmSel) in computeMXY line 213 using the exact GKL f32 inner sum as left operand; diagnostic-only"
        },
    );
    assert!(
        classification == "PAIRHMM_M_UPDATE_MXGAP_YGAPM_ADD_CAUSAL"
            || classification == "PAIRHMM_M_UPDATE_MXGAP_YGAPM_ADD_DIVERGENCE_NOT_CAUSAL"
    );
}
