//! GATK `Log10PairHMM` (exact log10) — parity with `HcParityNativePairHmm` / **GAP-F-02**.
//! R4-1: flat DP matrices + thread-local scratch (same numerics, far fewer allocations).

use gatk_common::{GatkError, GatkResult};
use std::cell::RefCell;
use std::sync::OnceLock;

const LOG10_3: f64 = 0.47712125471966244; // log10(3)
/// GATK `PairHMMModel.matchToMatchProbLog10`: log10(1 - min(1, 10^log10Sum)).
fn log10_match_to_match_from_error_log10_sum(log10_sum: f64) -> f64 {
    (1.0 - 10f64.powf(log10_sum).min(1.0)).log10()
}
const MAX_QUAL: usize = 127;

/// Default insertion/deletion/GCP qualities used by `HcParityNativePairHmm` parity dumps.
pub const GATK_PARITY_DEFAULT_INS_QUAL: u8 = 45;
pub const GATK_PARITY_DEFAULT_DEL_QUAL: u8 = 45;
pub const GATK_PARITY_DEFAULT_GCP: u8 = 10;

const MATCH_TO_MATCH: usize = 0;
const INDEL_TO_MATCH: usize = 1;
const MATCH_TO_INSERTION: usize = 2;
const INSERTION_TO_INSERTION: usize = 3;
const MATCH_TO_DELETION: usize = 4;
const DELETION_TO_DELETION: usize = 5;

fn qual_to_prob_log10(qual: u8) -> f64 {
    let q = qual as usize;
    if q > MAX_QUAL {
        let err = 10f64.powf(-(qual as f64) / 10.0);
        (1.0 - err).log10()
    } else {
        qual_to_prob_log10_table()[q]
    }
}

#[inline]
fn qual_to_error_prob_log10(qual: u8) -> f64 {
    -(qual as f64) * 0.1
}

#[inline]
fn log10_sum2(a: f64, b: f64) -> f64 {
    if a > b {
        return log10_sum2(b, a);
    }
    if a.is_infinite() && a.is_sign_negative() {
        return b;
    }
    if b.is_infinite() && b.is_sign_negative() {
        return a;
    }
    b + (1.0 + 10f64.powf(a - b)).log10()
}

#[cfg(test)]
#[inline]
fn exact_log10_sum_log10(values: &[f64]) -> f64 {
    let mut max_v = f64::NEG_INFINITY;
    let mut max_idx = 0usize;
    for (idx, &v) in values.iter().enumerate() {
        if v > max_v {
            max_v = v;
            max_idx = idx;
        }
    }
    if max_v.is_infinite() && max_v.is_sign_negative() {
        return max_v;
    }
    let mut sum = 1.0;
    for (idx, &v) in values.iter().enumerate() {
        if idx == max_idx || (v.is_infinite() && v.is_sign_negative()) {
            continue;
        }
        sum += 10f64.powf(v - max_v);
    }
    max_v + if sum != 1.0 { sum.log10() } else { 0.0 }
}

/// Same numerics as [`exact_log10_sum_log10`] on three values (no slice walk).
#[inline]
fn exact_log10_sum3(a: f64, b: f64, c: f64) -> f64 {
    let (max_v, x, y) = if a >= b && a >= c {
        (a, b, c)
    } else if b >= a && b >= c {
        (b, a, c)
    } else {
        (c, a, b)
    };
    if max_v.is_infinite() && max_v.is_sign_negative() {
        return max_v;
    }
    let mut sum = 1.0;
    if !(x.is_infinite() && x.is_sign_negative()) {
        sum += 10f64.powf(x - max_v);
    }
    if !(y.is_infinite() && y.is_sign_negative()) {
        sum += 10f64.powf(y - max_v);
    }
    max_v + if sum != 1.0 { sum.log10() } else { 0.0 }
}

fn approximate_log10_sum_log10(a: f64, b: f64) -> f64 {
    log10_sum2(a, b)
}

fn match_to_match_prob_log10(ins_qual: u8, del_qual: u8) -> f64 {
    let (min_q, max_q) = if ins_qual <= del_qual {
        (ins_qual as usize, del_qual as usize)
    } else {
        (del_qual as usize, ins_qual as usize)
    };
    if max_q > MAX_QUAL {
        let log10_sum = approximate_log10_sum_log10(-0.1 * min_q as f64, -0.1 * max_q as f64);
        log10_match_to_match_from_error_log10_sum(log10_sum)
    } else {
        match_to_match_log10_table()[((max_q * (max_q + 1)) >> 1) + min_q]
    }
}

fn qual_to_trans_probs_log10(ins_qual: u8, del_qual: u8, gcp: u8) -> [f64; 6] {
    let gcp_err = qual_to_error_prob_log10(gcp);
    [
        match_to_match_prob_log10(ins_qual, del_qual),
        qual_to_prob_log10(gcp),
        qual_to_error_prob_log10(ins_qual),
        gcp_err,
        qual_to_error_prob_log10(del_qual),
        gcp_err,
    ]
}

fn build_match_to_match_log10_table() -> [f64; ((MAX_QUAL + 1) * (MAX_QUAL + 2)) >> 1] {
    let mut table = [0.0_f64; ((MAX_QUAL + 1) * (MAX_QUAL + 2)) >> 1];
    let mut offset = 0usize;
    for i in 0..=MAX_QUAL {
        for j in 0..=i {
            let log10_sum = approximate_log10_sum_log10(-0.1 * j as f64, -0.1 * i as f64);
            table[offset + j] = log10_match_to_match_from_error_log10_sum(log10_sum);
        }
        offset += i + 1;
    }
    table
}

fn match_to_match_log10_table() -> &'static [f64; ((MAX_QUAL + 1) * (MAX_QUAL + 2)) >> 1] {
    static TABLE: OnceLock<[f64; ((MAX_QUAL + 1) * (MAX_QUAL + 2)) >> 1]> = OnceLock::new();
    TABLE.get_or_init(build_match_to_match_log10_table)
}

fn qual_to_prob_log10_table() -> &'static [f64; MAX_QUAL + 1] {
    static TABLE: OnceLock<[f64; MAX_QUAL + 1]> = OnceLock::new();
    TABLE.get_or_init(|| {
        let mut t = [0.0_f64; MAX_QUAL + 1];
        for (i, slot) in t.iter_mut().enumerate() {
            let err = 10f64.powf(-(i as f64) / 10.0);
            *slot = (1.0 - err).log10();
        }
        t
    })
}

struct PairHmmScratch {
    transition: Vec<[f64; 6]>,
    prior: Vec<f64>,
    m: Vec<f64>,
    ins: Vec<f64>,
    del: Vec<f64>,
}

impl PairHmmScratch {
    fn new() -> Self {
        Self {
            transition: Vec::new(),
            prior: Vec::new(),
            m: Vec::new(),
            ins: Vec::new(),
            del: Vec::new(),
        }
    }

    fn ensure(&mut self, rn: usize, hn: usize) {
        let rows = rn + 1;
        let cols = hn + 1;
        let cells = rows * cols;
        if self.transition.len() < rows {
            self.transition.resize(rows, [0.0; 6]);
        }
        if self.prior.len() < cells {
            self.prior.resize(cells, 0.0);
            self.m.resize(cells, f64::NEG_INFINITY);
            self.ins.resize(cells, f64::NEG_INFINITY);
            self.del.resize(cells, f64::NEG_INFINITY);
        }
    }

    /// Drop oversized TLS capacity so Peak-RSS can fall after a deep region.
    /// Keeps a modest high-water mark to avoid thrashing on typical region sizes.
    fn shrink_to_budget(&mut self, max_keep_cells: usize) {
        if self.prior.capacity() <= max_keep_cells.saturating_mul(2) {
            return;
        }
        *self = Self::new();
    }
}

thread_local! {
    static PAIRHMM_SCRATCH: RefCell<PairHmmScratch> = RefCell::new(PairHmmScratch::new());
}

/// Soft ceiling retained after [`release_pairhmm_tls_scratch`] (cells across DP planes).
const PAIRHMM_TLS_KEEP_CELLS: usize = 256 * 1024;

/// Release PairHMM Log10 TLS scratch when it grew past a modest region-scale budget.
/// Call after finishing an assembly region (or when Peak-RSS pressure is high).
pub fn release_pairhmm_tls_scratch() {
    PAIRHMM_SCRATCH.with(|cell| {
        cell.borrow_mut().shrink_to_budget(PAIRHMM_TLS_KEEP_CELLS);
    });
}

/// Log10 P(read | haplotype) using GATK `Log10PairHMM` semantics (exact log10 sums).
pub fn log10_pairhmm_likelihood(
    read_bases: &[u8],
    read_quals: &[u8],
    haplotype_bases: &[u8],
    insertion_gop: &[u8],
    deletion_gop: &[u8],
    overall_gcp: &[u8],
) -> GatkResult<f64> {
    if read_bases.len() != read_quals.len()
        || read_bases.len() != insertion_gop.len()
        || read_bases.len() != deletion_gop.len()
        || read_bases.len() != overall_gcp.len()
    {
        return Err(GatkError::argument(
            "PairHMM read arrays must have equal length",
        ));
    }
    if haplotype_bases.is_empty() {
        return Err(GatkError::argument("haplotype must be non-empty"));
    }

    let rn = read_bases.len();
    let hn = haplotype_bases.len();
    if rn == 0 {
        return Ok(0.0);
    }
    // Guard: contig-scale inputs (e.g. full chr20) make DP matrices tens of GiB and OOM the process.
    // GATK assembly regions are hundreds of bp; keep Peak-RSS fail-closed well below 16 GiB hosts.
    // 8e6 cells × 4 f64 planes ≈ 256 MiB TLS high-water (was 50e6 ≈ 1.6 GiB before refuse).
    const MAX_PAIRHMM_DIM: usize = 100_000;
    const MAX_PAIRHMM_CELLS: usize = 8_000_000;
    let cells = (rn + 1).saturating_mul(hn + 1);
    if rn > MAX_PAIRHMM_DIM || hn > MAX_PAIRHMM_DIM || cells > MAX_PAIRHMM_CELLS {
        return Err(GatkError::algorithm(format!(
            "PairHMM Log10 refused oversized DP (read_len={rn}, hap_len={hn}, cells={cells}); \
             inputs must be assembly-region scale, not contig scale"
        )));
    }

    PAIRHMM_SCRATCH.with(|cell| {
        let mut scratch = cell.borrow_mut();
        scratch.ensure(rn, hn);
        Ok(log10_pairhmm_likelihood_into(
            read_bases,
            read_quals,
            haplotype_bases,
            insertion_gop,
            deletion_gop,
            overall_gcp,
            &mut scratch,
        ))
    })
}

fn log10_pairhmm_likelihood_into(
    read_bases: &[u8],
    read_quals: &[u8],
    haplotype_bases: &[u8],
    insertion_gop: &[u8],
    deletion_gop: &[u8],
    overall_gcp: &[u8],
    scratch: &mut PairHmmScratch,
) -> f64 {
    let rn = read_bases.len();
    let hn = haplotype_bases.len();
    let cols = hn + 1;

    for i in 0..rn {
        scratch.transition[i + 1] =
            qual_to_trans_probs_log10(insertion_gop[i], deletion_gop[i], overall_gcp[i]);
    }

    let cells = (rn + 1) * cols;
    for slot in &mut scratch.m[..cells] {
        *slot = f64::NEG_INFINITY;
    }
    for slot in &mut scratch.ins[..cells] {
        *slot = f64::NEG_INFINITY;
    }
    for slot in &mut scratch.del[..cells] {
        *slot = f64::NEG_INFINITY;
    }

    for i in 0..rn {
        let x = read_bases[i];
        let qual = read_quals[i];
        let match_p = qual_to_prob_log10(qual);
        let mismatch_p = qual_to_error_prob_log10(qual) - LOG10_3;
        let row = (i + 1) * cols;
        for j in 0..hn {
            let y = haplotype_bases[j];
            scratch.prior[row + j + 1] = if x == y || x == b'N' || y == b'N' {
                match_p
            } else {
                mismatch_p
            };
        }
    }

    let init_del = (1.0 / hn as f64).log10();
    for j in 0..=hn {
        scratch.del[j] = init_del;
    }

    for i in 1..=rn {
        let t = scratch.transition[i];
        let row = i * cols;
        let prev = (i - 1) * cols;
        for j in 1..=hn {
            let p = scratch.prior[row + j];
            scratch.m[row + j] = p + exact_log10_sum3(
                scratch.m[prev + j - 1] + t[MATCH_TO_MATCH],
                scratch.ins[prev + j - 1] + t[INDEL_TO_MATCH],
                scratch.del[prev + j - 1] + t[INDEL_TO_MATCH],
            );
            scratch.ins[row + j] = log10_sum2(
                scratch.m[prev + j] + t[MATCH_TO_INSERTION],
                scratch.ins[prev + j] + t[INSERTION_TO_INSERTION],
            );
            // Two-operand sum (same as historical `exact_log10_sum_log10(&[a, b])`).
            scratch.del[row + j] = log10_sum2(
                scratch.m[row + j - 1] + t[MATCH_TO_DELETION],
                scratch.del[row + j - 1] + t[DELETION_TO_DELETION],
            );
        }
    }

    let end_row = rn * cols;
    let mut final_sum = log10_sum2(scratch.m[end_row + 1], scratch.ins[end_row + 1]);
    for j in 2..=hn {
        final_sum = exact_log10_sum3(final_sum, scratch.m[end_row + j], scratch.ins[end_row + j]);
    }
    final_sum
}

/// Parity-path likelihood with GATK default indel/GCP qualities (45/45/10).
pub fn log10_pairhmm_likelihood_parity_defaults(
    read_bases: &[u8],
    read_quals: &[u8],
    haplotype_bases: &[u8],
) -> GatkResult<f64> {
    let n = read_bases.len();
    let ins = vec![GATK_PARITY_DEFAULT_INS_QUAL; n];
    let del = vec![GATK_PARITY_DEFAULT_DEL_QUAL; n];
    let gcp = vec![GATK_PARITY_DEFAULT_GCP; n];
    log10_pairhmm_likelihood(read_bases, read_quals, haplotype_bases, &ins, &del, &gcp)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn log10_sum_single_finite_matches_java() {
        let s = exact_log10_sum_log10(&[f64::NEG_INFINITY, f64::NEG_INFINITY, -1.045757490560675]);
        assert!((s - (-1.045757490560675)).abs() < 1e-12, "s={s}");
    }

    #[test]
    fn exact_log10_sum3_matches_slice() {
        let a = -2.5;
        let b = -1.0;
        let c = -3.25;
        let via3 = exact_log10_sum3(a, b, c);
        let via_slice = exact_log10_sum_log10(&[a, b, c]);
        assert!((via3 - via_slice).abs() < 1e-15);
    }

    #[test]
    fn first_cell_m11_matches_java() {
        let read = b"ACGTACGTAC";
        let quals = vec![32u8; 10];
        let hap = b"ACGTACGTAC";
        let n = read.len();
        let ins = vec![GATK_PARITY_DEFAULT_INS_QUAL; n];
        let del = vec![GATK_PARITY_DEFAULT_DEL_QUAL; n];
        let gcp = vec![GATK_PARITY_DEFAULT_GCP; n];
        let rn = read.len();
        let hn = hap.len();
        let mut scratch = PairHmmScratch::new();
        scratch.ensure(rn, hn);
        for i in 0..rn {
            scratch.transition[i + 1] = qual_to_trans_probs_log10(ins[i], del[i], gcp[i]);
        }
        let cols = hn + 1;
        for i in 0..rn {
            for j in 0..hn {
                let x = read[i];
                let y = hap[j];
                scratch.prior[(i + 1) * cols + j + 1] = if x == y || x == b'N' || y == b'N' {
                    qual_to_prob_log10(quals[i])
                } else {
                    qual_to_error_prob_log10(quals[i]) - LOG10_3
                };
            }
        }
        for slot in &mut scratch.m {
            *slot = f64::NEG_INFINITY;
        }
        for slot in &mut scratch.ins {
            *slot = f64::NEG_INFINITY;
        }
        for slot in &mut scratch.del {
            *slot = f64::NEG_INFINITY;
        }
        let init_del = (1.0 / hn as f64).log10();
        for j in 0..=hn {
            scratch.del[j] = init_del;
        }
        let i = 1;
        let j = 1;
        let t = scratch.transition[i];
        let p = scratch.prior[i * cols + j];
        let m11 = p + exact_log10_sum3(
            scratch.m[(i - 1) * cols + j - 1] + t[MATCH_TO_MATCH],
            scratch.ins[(i - 1) * cols + j - 1] + t[INDEL_TO_MATCH],
            scratch.del[(i - 1) * cols + j - 1] + t[INDEL_TO_MATCH],
        );
        assert!((m11 - (-1.0460315983379533)).abs() < 1e-9, "m[1][1]={m11}");
    }

    #[test]
    fn parity_fixture_case1_near_java_native() {
        let read = b"ACGTACGTAC";
        let quals = vec![32u8; 10];
        let hap = b"ACGTACGTAC";
        let ll = log10_pairhmm_likelihood_parity_defaults(read, &quals, hap).unwrap();
        assert!((ll - (-1.0487303647971473)).abs() < 1e-6, "ll={ll}");
    }

    #[test]
    fn parity_fixture_case2_near_java_native() {
        let read = b"ACGTACGTAC";
        let quals = vec![32u8; 10];
        let hap = b"ACGTTCGTAC";
        let ll = log10_pairhmm_likelihood_parity_defaults(read, &quals, hap).unwrap();
        assert!((ll - (-4.7255750676296255)).abs() < 1e-6, "ll={ll}");
    }

    #[test]
    fn parity_fixture_case3_near_java_native() {
        let read = b"ACGTACGTAC";
        let quals = vec![32u8; 10];
        let hap = b"ACGTACGTTC";
        let ll = log10_pairhmm_likelihood_parity_defaults(read, &quals, hap).unwrap();
        assert!((ll - (-4.7180722504710590)).abs() < 1e-6, "ll={ll}");
    }
}
