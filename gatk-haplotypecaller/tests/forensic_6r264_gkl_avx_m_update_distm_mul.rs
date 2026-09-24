//! 6R.264: measurement-only. GKL AVX `computeMXY` outer `VEC_MUL` of
//! `distmSel` vs Rust f64 `p * sum` at `20:29455649 T/TGTTTG` under the
//! frozen 6R.257 joint input plane. PRODUCTION CHANGE: NONE.
//!
//! 6R.259–263 closed ph2pr / match / mismatch / distm blend (not causal).
//! This round injects only the GKL f32 result of `float(M_operand)*float(distm)`.
//! Prefer first-cell-only; also record labeled all-M-multiply CF.
//! Native AVX GKL is not executed (oracle x86_64; this host aarch64).
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r264_gkl_avx_m_update_distm_mul -- --nocapture --test-threads=1
//! HOLDOUT_6R264=1 cargo test -p gatk-haplotypecaller --test holdout_6r264_gkl_avx_m_update_distm_mul -- --nocapture --test-threads=1
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

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const GKL_PIN: &str =
    "0.8.8 Intel-HLS/GKL avx-pairhmm-template.h computeMXY VEC_MUL(..., distmSel)";
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

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R264\t{key}\t{}", value.as_ref());
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
enum MulCf {
    Baseline,
    FirstCellGklF32,
    AllMGklF32,
}

/// IEEE-754 f32 multiply of the INITIAL-normalized M-operand and GKL distm,
/// rescaled back to Rust INITIAL. Direct `sum as f32` overflows (2^1020).
fn gkl_m_distm_mul_rescaled(sum_f64: f64, q: u8, equal: bool) -> f64 {
    let ratio = sum_f64 / INITIAL_CONDITION;
    INITIAL_CONDITION * f64::from((ratio as f32) * gkl_distm_f32(q, equal))
}

#[derive(Clone, Copy, PartialEq)]
struct CellDump {
    category: &'static str,
    retain_row: usize,
    hap_index: usize,
    read_i: usize,
    hap_j: usize,
    q: u8,
    gcp: u8,
    read_base: u8,
    hap_base: u8,
    equal: bool,
    rust_m_operand: f64,
    rust_distm: f64,
}

fn take_first(slot: &mut Option<CellDump>, cell: CellDump) {
    if slot.is_none() {
        *slot = Some(cell);
    }
}

fn rust_diag_ll(
    read: &[u8],
    quals: &[u8],
    hap: &[u8],
    iq: &[u8],
    dq: &[u8],
    gcp: &[u8],
    mode: MulCf,
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
            let sum = m[prev + j - 1] * t[0] + ins[prev + j - 1] * t[1] + del[prev + j - 1] * t[1];
            let prod = match mode {
                MulCf::Baseline => p * sum,
                MulCf::FirstCellGklF32 if i == 1 && j == 1 => {
                    gkl_m_distm_mul_rescaled(sum, q, equal)
                }
                MulCf::FirstCellGklF32 => p * sum,
                MulCf::AllMGklF32 => gkl_m_distm_mul_rescaled(sum, q, equal),
            };
            m[row + j] = prod;
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

fn collect_m_mul_cells(
    read: &[u8],
    quals: &[u8],
    hap: &[u8],
    iq: &[u8],
    dq: &[u8],
    gcp: &[u8],
    retain_row: usize,
    hap_index: usize,
    first_op: &mut Option<CellDump>,
    first_q20_match: &mut Option<CellDump>,
    first_q20_mismatch: &mut Option<CellDump>,
    first_q6: &mut Option<CellDump>,
    first_inherited: &mut Option<CellDump>,
    first_new_f32: &mut Option<CellDump>,
) {
    if first_op.is_some()
        && first_q20_match.is_some()
        && first_q20_mismatch.is_some()
        && first_q6.is_some()
        && first_inherited.is_some()
        && first_new_f32.is_some()
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
        let gcp_i = gcp[i - 1];
        let row = i * cols;
        let prev = (i - 1) * cols;
        m[row] = 0.0;
        ins[row] = 0.0;
        del[row] = 0.0;
        for j in 1..=hn {
            let y = hap[j - 1];
            let equal = bases_equal(x, y);
            let p = rust_distm(q, equal);
            let sum = m[prev + j - 1] * t[0] + ins[prev + j - 1] * t[1] + del[prev + j - 1] * t[1];
            let cell = || CellDump {
                category: "",
                retain_row,
                hap_index,
                read_i: i,
                hap_j: j,
                q,
                gcp: gcp_i,
                read_base: x,
                hap_base: y,
                equal,
                rust_m_operand: sum,
                rust_distm: p,
            };
            if i == 1 && j == 1 {
                take_first(
                    first_op,
                    CellDump {
                        category: "first_op_i1j1",
                        ..cell()
                    },
                );
            }
            if q == 20 && equal && sum != 0.0 {
                take_first(
                    first_q20_match,
                    CellDump {
                        category: "q20_match",
                        ..cell()
                    },
                );
            }
            if q == 20 && !equal && sum != 0.0 {
                take_first(
                    first_q20_mismatch,
                    CellDump {
                        category: "q20_mismatch",
                        ..cell()
                    },
                );
            }
            if q == 6 && sum != 0.0 {
                take_first(
                    first_q6,
                    CellDump {
                        category: "q6",
                        ..cell()
                    },
                );
            }
            if sum != 0.0 {
                match classify_mismatch(q) {
                    "inherited-only" => take_first(
                        first_inherited,
                        CellDump {
                            category: "inherited_only",
                            ..cell()
                        },
                    ),
                    "new-f32-division-rounding" => take_first(
                        first_new_f32,
                        CellDump {
                            category: "new_f32_rounding",
                            ..cell()
                        },
                    ),
                    _ => {}
                }
            }
            m[row + j] = p * sum;
            ins[row + j] = m[prev + j] * t[2] + ins[prev + j] * t[3];
            del[row + j] = m[row + j - 1] * t[4] + del[row + j - 1] * t[5];
            if first_op.is_some()
                && first_q20_match.is_some()
                && first_q20_mismatch.is_some()
                && first_q6.is_some()
                && first_inherited.is_some()
                && first_new_f32.is_some()
            {
                return;
            }
        }
    }
}

fn dump_cell(cell: &CellDump) {
    let q = cell.q;
    let gkl_distm = gkl_distm_f32(q, cell.equal);
    let rust_sum = cell.rust_m_operand;
    let rust_d = cell.rust_distm;
    let rust_prod = rust_sum * rust_d;
    let rust_ratio = rust_sum / INITIAL_CONDITION;
    let gkl_m_f32 = rust_ratio as f32;
    let gkl_prod_f32 = gkl_m_f32 * gkl_distm;
    let gkl_prod_w = INITIAL_CONDITION * f64::from(gkl_prod_f32);
    let isolated = (rust_ratio as f32) * (rust_d as f32);
    let isolated_w = INITIAL_CONDITION * f64::from(isolated);
    let abs = (gkl_prod_w - rust_prod).abs();
    let rel = if rust_prod != 0.0 {
        abs / rust_prod.abs()
    } else {
        0.0
    };
    let isolated_abs = (isolated_w - rust_prod).abs();
    kv(
        "m_mul_cell",
        format!(
            "cat={}\tretain_row={}\thap={}\tread_i={}\thap_j={}\tq={}\tgcp={}\tread_base={}\thap_base={}\tequal={}\trust_m_ratio={}\tgkl_m_operand_f32_bits=0x{:08x}\tgkl_distm_f32_bits=0x{:08x}\trust_distm={}\tgkl_f32_mul_bits=0x{:08x}\tgkl_f32_mul_rescaled={}\trust_f64_mul={}\tisolated_f32_same_operands_bits=0x{:08x}\tisolated_rescaled={}\tabs_delta={:.17e}\trel_delta={:.17e}\tisolated_abs_delta={:.17e}\tnew_f32_mul_rounding={}",
            cell.category,
            cell.retain_row,
            cell.hap_index,
            cell.read_i,
            cell.hap_j,
            q,
            cell.gcp,
            cell.read_base as char,
            cell.hap_base as char,
            cell.equal,
            fmt_f64(rust_ratio),
            gkl_m_f32.to_bits(),
            gkl_distm.to_bits(),
            fmt_f64(rust_d),
            gkl_prod_f32.to_bits(),
            fmt_f64(gkl_prod_w),
            fmt_f64(rust_prod),
            isolated.to_bits(),
            fmt_f64(isolated_w),
            abs,
            rel,
            isolated_abs,
            isolated_w.to_bits() != rust_prod.to_bits(),
        ),
    );
    let inherit_key = format!("m_mul_{}_inherited_distm_vs_new_mul", cell.category);
    kv(
        &inherit_key,
        format!(
            "gkl_prod_vs_isolated_differs={} new_f32_mul_vs_rust={}",
            gkl_prod_f32.to_bits() != isolated.to_bits(),
            isolated_w.to_bits() != rust_prod.to_bits()
        ),
    );
}

#[test]
fn forensic_6r264_gkl_avx_m_update_distm_mul() {
    kv("java_pin", JAVA_PIN);
    kv("gkl_pin", GKL_PIN);
    kv("production_change", "NONE");
    kv("target", "20:29455649 T/TGTTTG");
    kv(
        "frozen_input",
        "6R.257: Java174-mers + Java clip 20:29455560-29455728; 123 retainEvidence",
    );
    kv(
        "gkl_m_mul_op",
        "computeMXY: M_t.d=VEC_MUL(VEC_ADD(VEC_ADD(VEC_MUL(M_t_2,pMM), VEC_MUL(X_t_2,pGAPM)), VEC_MUL(Y_t_2,pGAPM)), distmSel); outer VEC_MUL is the first distm multiply; NUMBER=float",
    );
    kv(
        "gkl_source",
        "Intel-HLS/GKL 0.8.8 src/main/native/pairhmm/avx-pairhmm-template.h computeMXY; avx-functions-float.h VEC_MUL=_mm256_mul_ps; avx512-functions-float.h VEC_MUL=_mm512_mul_ps; IntelPairHmm.cc g_compute_full_prob_float",
    );
    kv(
        "oracle_vs_host",
        "Java oracle is x86_64 (AVX or AVX-512 if supported, not Apple); this host is aarch64 Darwin; native GKL AVX is not executed",
    );
    kv(
        "first_cell_init",
        "stripe i==0: M_t_2=0 X_t_2=0 Y_t_2=VEC_SET_LSE(init_Y); init_Y=INITIAL_CONSTANT/(NUMBER)haplen; first M operand of outer mul is (0*pMM+0*pGAPM+init_Y*pGAPM) in float before distmSel mul",
    );
    kv(
        "rust_formula",
        "m[i,j]=p*(m[i-1,j-1]*t_mm + ins[i-1,j-1]*t_gapm + del[i-1,j-1]*t_gapm); first cell i=1,j=1: m=p*(init_del*t_gapm); f64 mul",
    );
    kv(
        "rust_source",
        "gatk-haplotypecaller/src/pairhmm_logless.rs ~347-350; pairhmm_simd/neon.rs ~304-308 vmulq_f64(p, sum)",
    );
    kv(
        "counterfactual_rule",
        "CF1 first-cell-only: only (i=1,j=1) M product replaced with INITIAL*(f32(sum/INITIAL)*gkl_distm); CF2 labeled all-M-multiply GKL-f32 of the same ratio-space mul; transitions/init/I/D/log10 stay Rust; never narrow 2^1020 to f32",
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
    let mut cf2_pooled = vec![0.0f64; 123];
    let mut first_op: Option<CellDump> = None;
    let mut first_q20_match: Option<CellDump> = None;
    let mut first_q20_mismatch: Option<CellDump> = None;
    let mut first_q6: Option<CellDump> = None;
    let mut first_inherited: Option<CellDump> = None;
    let mut first_new_f32: Option<CellDump> = None;

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
            collect_m_mul_cells(
                &p.bases,
                &p.bq,
                &java_haps[k],
                &p.iq,
                &p.dq,
                &p.gcp,
                ri,
                k,
                &mut first_op,
                &mut first_q20_match,
                &mut first_q20_mismatch,
                &mut first_q6,
                &mut first_inherited,
                &mut first_new_f32,
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
        let mut cf2_raw = [0.0f64; 5];
        for k in 0..5 {
            cf1_raw[k] = rust_diag_ll(
                &p.bases,
                &p.bq,
                &java_haps[k],
                &p.iq,
                &p.dq,
                &p.gcp,
                MulCf::FirstCellGklF32,
            );
            cf2_raw[k] = rust_diag_ll(
                &p.bases,
                &p.bq,
                &java_haps[k],
                &p.iq,
                &p.dq,
                &p.gcp,
                MulCf::AllMGklF32,
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
                MulCf::Baseline,
            );
            assert_eq!(
                baseline.to_bits(),
                raw[k].to_bits(),
                "diagnostic baseline must match NEON/scalar"
            );
        }
        cf1_pooled[ri] = floor5(cf1_raw, other_best)
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        cf2_pooled[ri] = floor5(cf2_raw, other_best)
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

    assert!(unique_q.contains(&20), "frozen object must include Q=20");
    assert!(unique_q.contains(&6), "frozen object must include Q=6");
    assert_eq!(gkl_ph2pr(20).to_bits(), Q20_GKL_PH2PR_F32);
    assert_eq!(gkl_match(20).to_bits(), Q20_GKL_MATCH_F32);
    assert_eq!(gkl_mismatch(20).to_bits(), Q20_GKL_MISMATCH_F32);

    let first_op = first_op.expect("first M*distm cell");
    assert_eq!(first_op.read_i, 1);
    assert_eq!(first_op.hap_j, 1);
    let hn = 174usize;
    let gkl_init = GKL_INITIAL_F32 / hn as f32;
    let gkl_pgapm = 1.0f32 - gkl_ph2pr(first_op.gcp);
    let gkl_native_m_op = gkl_init * gkl_pgapm;
    kv("gkl_initial_f32", fmt_f32(GKL_INITIAL_F32));
    kv("gkl_init_del_f32", fmt_f32(gkl_init));
    kv("gkl_native_first_m_operand", fmt_f32(gkl_native_m_op));
    let gkl_native_prod = gkl_native_m_op * gkl_distm_f32(first_op.q, first_op.equal);
    kv("gkl_native_first_m_distm_mul", fmt_f32(gkl_native_prod));
    kv(
        "gkl_native_first_mul_over_initial",
        fmt_f64(f64::from(gkl_native_prod) / f64::from(GKL_INITIAL_F32)),
    );
    kv(
        "rust_first_mul_over_initial",
        fmt_f64((first_op.rust_m_operand * first_op.rust_distm) / INITIAL_CONDITION),
    );
    kv(
        "gkl_native_scale_note",
        "GKL INITIAL is 2^120 f32; Rust INITIAL is 2^1020 f64. Direct sum-as-f32 overflows. CF multiplies INITIAL-normalized M-operand in f32, then rescales.",
    );

    for cell in [
        &Some(first_op.clone()),
        &first_q20_match,
        &first_q20_mismatch,
        &first_q6,
        &first_inherited,
        &first_new_f32,
    ] {
        dump_cell(cell.as_ref().expect("missing representative cell"));
    }

    kv(
        "counterfactual_note",
        "CF1=first-cell-only (i=1,j=1) INITIAL*(f32(sum/INITIAL)*gkl_distm); CF2=all-M-multiply GKL-f32 of the same ratio-space mul; I/D/transitions/init/log10 stay Rust NEON-equivalent scalar",
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
    kv("first_divergent_operation", "GKL_AVX_M_UPDATE_DISTM_MUL");
    kv(
        "first_divergent_detail",
        format!(
            "outer VEC_MUL=_mm256_mul_ps(sum, distmSel) float; first_op q={} equal={} rust_sum_bits=0x{:016x}",
            first_op.q,
            first_op.equal,
            first_op.rust_m_operand.to_bits()
        ),
    );
    kv("cf1_label", "first-cell-only");
    kv("cf1_emitted_GL", fmt_f64(cf1_rel));
    kv("cf1_continuous_PL", fmt_f64(cf1_cont));
    kv("cf1_integer_PL", cf1_pl.to_string());
    kv("cf1_delta_GL", fmt_f64(cf1_rel - rr_rel));
    kv("cf1_delta_continuous_PL", fmt_f64(cf1_cont - rr_cont));
    kv("cf1_delta_integer_PL", (cf1_pl - rr_pl).to_string());
    kv("required_pl_movement", format!("{REQUIRED_MOVE:.17}"));
    kv("cf1_crosses_3517_5", (cf1_cont < BOUNDARY).to_string());

    let cf2_gls = biallelic_gls(&l0, &cf2_pooled);
    let (cf2_rel, cf2_cont, cf2_pl) = emitted_hom_alt(&cf2_gls);
    kv("cf2_label", "all-M-multiply GKL-f32 counterfactual");
    kv("cf2_emitted_GL", fmt_f64(cf2_rel));
    kv("cf2_continuous_PL", fmt_f64(cf2_cont));
    kv("cf2_integer_PL", cf2_pl.to_string());
    kv("cf2_delta_GL", fmt_f64(cf2_rel - rr_rel));
    kv("cf2_delta_continuous_PL", fmt_f64(cf2_cont - rr_cont));
    kv("cf2_delta_integer_PL", (cf2_pl - rr_pl).to_string());
    kv("cf2_crosses_3517_5", (cf2_cont < BOUNDARY).to_string());

    let classification = if cf1_cont.to_bits() == rr_cont.to_bits()
        && cf1_pl == rr_pl
        && first_op.rust_distm.to_bits()
            == f64::from(gkl_distm_f32(first_op.q, first_op.equal)).to_bits()
    {
        "PAIRHMM_M_MULTIPLY_RESULT_MATCH"
    } else if cf1_cont < BOUNDARY && cf1_pl == 3517 {
        "PAIRHMM_M_MULTIPLY_DIVERGENCE_CAUSAL"
    } else {
        "PAIRHMM_M_MULTIPLY_DIVERGENCE_NOT_CAUSAL"
    };
    kv("classification", classification);
    kv("classification_basis", "CF1 first-cell-only");
    kv("production_change", "NONE");
    kv(
        "next_arrow",
        if classification == "PAIRHMM_M_MULTIPLY_DIVERGENCE_CAUSAL" {
            "M*distm f32 mul is PL-causal; production still unchanged this round"
        } else {
            "GKL_AVX_X_UPDATE: computeMXY X_t=VEC_ADD(VEC_MUL(M_t_1,pMX), VEC_MUL(X_t_1,pXX)); first arithmetic after distm mul; diagnostic-only"
        },
    );
    assert!(
        classification == "PAIRHMM_M_MULTIPLY_DIVERGENCE_CAUSAL"
            || classification == "PAIRHMM_M_MULTIPLY_DIVERGENCE_NOT_CAUSAL"
            || classification == "PAIRHMM_M_MULTIPLY_RESULT_MATCH"
    );
}
