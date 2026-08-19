//! Architecture C phase-1: row-wavefront Logless PairHMM.
//!
//! Rolling 2-row DP with optional AVX2/NEON column striping for Match/Insertion.
//! Deletion stays serial (left-cell dependence). f32 primary + f64 retry.
//! Scalar [`crate::pairhmm_logless`] remains the correctness oracle.
//!
//! Opt-in via [`super::dispatch::PairHmmImpl::Wavefront`] — not wired into
//! `FASTEST_AVAILABLE` until gates pass.

mod prep;
mod rolling_f32;
mod rolling_f64;

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
mod avx2_f32;
#[cfg(target_arch = "aarch64")]
mod neon_f32;

pub use prep::ReadPrep;
pub use rolling_f64::score_one_rolling_f64;

use gatk_common::{GatkError, GatkResult};
use prep::ReadPrep as Prep;
use rolling_f32::score_one_rolling_f32_with_retry;
use std::cell::RefCell;

const MAX_PAIRHMM_DIM: usize = 100_000;
const MAX_PAIRHMM_CELLS: usize = 8_000_000;

/// Contiguous rolling DP scratch (2 rows × cols for m/ins/del).
pub struct WavefrontScratch {
    pub m: Vec<f64>,
    pub ins: Vec<f64>,
    pub del: Vec<f64>,
    pub m32: Vec<f32>,
    pub ins32: Vec<f32>,
    pub del32: Vec<f32>,
}

impl WavefrontScratch {
    fn empty() -> Self {
        Self {
            m: Vec::new(),
            ins: Vec::new(),
            del: Vec::new(),
            m32: Vec::new(),
            ins32: Vec::new(),
            del32: Vec::new(),
        }
    }

    pub fn ensure_rolling_f64(&mut self, cols: usize) {
        let need = cols.saturating_mul(2);
        if self.m.len() < need {
            self.m.resize(need, 0.0);
            self.ins.resize(need, 0.0);
            self.del.resize(need, 0.0);
        }
    }

    pub fn ensure_rolling_f32(&mut self, cols: usize) {
        let need = cols.saturating_mul(2);
        if self.m32.len() < need {
            self.m32.resize(need, 0.0);
            self.ins32.resize(need, 0.0);
            self.del32.resize(need, 0.0);
        }
    }
}

thread_local! {
    static WF_SCRATCH: RefCell<WavefrontScratch> = RefCell::new(WavefrontScratch::empty());
}

/// Keep wavefront TLS high-water (Peak path).
pub fn release_wavefront_tls_scratch() {
    // High-water retention — no munmap (matches other PairHMM TLS).
}

/// Score one read against many haplotypes with shared [`ReadPrep`] (HC shape).
pub fn score_haps_wavefront_f32(
    read_bases: &[u8],
    read_quals: &[u8],
    haplotypes: &[&[u8]],
    insertion_gop: &[u8],
    deletion_gop: &[u8],
    overall_gcp: &[u8],
) -> GatkResult<Vec<f64>> {
    if haplotypes.is_empty() {
        return Ok(Vec::new());
    }
    let rn = read_bases.len();
    if rn == 0 {
        return Ok(vec![0.0; haplotypes.len()]);
    }
    validate_read_arrays(
        read_bases,
        read_quals,
        insertion_gop,
        deletion_gop,
        overall_gcp,
    )?;
    let max_hn = haplotypes.iter().map(|h| h.len()).max().unwrap_or(0);
    let max_cells = (rn + 1).saturating_mul(max_hn + 1);
    if rn > MAX_PAIRHMM_DIM || max_hn > MAX_PAIRHMM_DIM || max_cells > MAX_PAIRHMM_CELLS {
        return Err(GatkError::algorithm(format!(
            "PairHMM wavefront refused oversized DP (read_len={rn}, max_hap_len={max_hn}, cells={max_cells})"
        )));
    }
    for hap in haplotypes {
        if hap.is_empty() {
            return Err(GatkError::argument("haplotype must be non-empty"));
        }
    }

    let prep = Prep::build(
        read_bases,
        read_quals,
        insertion_gop,
        deletion_gop,
        overall_gcp,
    );
    let mut out = Vec::with_capacity(haplotypes.len());
    WF_SCRATCH.with(|cell| {
        let mut scratch = cell.borrow_mut();
        for &hap in haplotypes {
            out.push(score_one_with_prep(&prep, read_bases, hap, &mut scratch));
        }
    });
    Ok(out)
}

/// Score with a pre-built [`ReadPrep`] (amortize prep across tiles/benches).
pub fn score_haps_wavefront_f32_with_prep(
    prep: &Prep,
    read_bases: &[u8],
    haplotypes: &[&[u8]],
) -> GatkResult<Vec<f64>> {
    if haplotypes.is_empty() {
        return Ok(Vec::new());
    }
    if prep.read_len() != read_bases.len() {
        return Err(GatkError::argument(
            "ReadPrep length must match read_bases length",
        ));
    }
    let mut out = Vec::with_capacity(haplotypes.len());
    WF_SCRATCH.with(|cell| {
        let mut scratch = cell.borrow_mut();
        for &hap in haplotypes {
            if hap.is_empty() {
                // Collect error outside — use sentinel then rewrite.
                out.clear();
                return;
            }
            out.push(score_one_with_prep(prep, read_bases, hap, &mut scratch));
        }
    });
    if out.len() != haplotypes.len() {
        return Err(GatkError::argument("haplotype must be non-empty"));
    }
    Ok(out)
}

/// Rolling f64 only (oracle layout for tests / retry).
pub fn score_haps_wavefront_rolling_f64(
    read_bases: &[u8],
    read_quals: &[u8],
    haplotypes: &[&[u8]],
    insertion_gop: &[u8],
    deletion_gop: &[u8],
    overall_gcp: &[u8],
) -> GatkResult<Vec<f64>> {
    if haplotypes.is_empty() {
        return Ok(Vec::new());
    }
    validate_read_arrays(
        read_bases,
        read_quals,
        insertion_gop,
        deletion_gop,
        overall_gcp,
    )?;
    let prep = Prep::build(
        read_bases,
        read_quals,
        insertion_gop,
        deletion_gop,
        overall_gcp,
    );
    let mut out = Vec::with_capacity(haplotypes.len());
    WF_SCRATCH.with(|cell| {
        let mut scratch = cell.borrow_mut();
        for &hap in haplotypes {
            if hap.is_empty() {
                out.clear();
                return;
            }
            out.push(score_one_rolling_f64(&prep, read_bases, hap, &mut scratch));
        }
    });
    if out.len() != haplotypes.len() {
        return Err(GatkError::argument("haplotype must be non-empty"));
    }
    Ok(out)
}

fn score_one_with_prep(
    prep: &Prep,
    read_bases: &[u8],
    hap: &[u8],
    scratch: &mut WavefrontScratch,
) -> f64 {
    score_one_rolling_f32_with_retry(prep, read_bases, hap, scratch)
}

fn validate_read_arrays(
    read_bases: &[u8],
    read_quals: &[u8],
    insertion_gop: &[u8],
    deletion_gop: &[u8],
    overall_gcp: &[u8],
) -> GatkResult<()> {
    if read_bases.len() != read_quals.len()
        || read_bases.len() != insertion_gop.len()
        || read_bases.len() != deletion_gop.len()
        || read_bases.len() != overall_gcp.len()
    {
        return Err(GatkError::argument(
            "PairHMM read arrays must have equal length",
        ));
    }
    Ok(())
}

/// Score forcing the portable f32 kernel (no AVX2/NEON) — for A/B vs host SIMD.
pub fn score_haps_wavefront_portable_f32(
    read_bases: &[u8],
    read_quals: &[u8],
    haplotypes: &[&[u8]],
    insertion_gop: &[u8],
    deletion_gop: &[u8],
    overall_gcp: &[u8],
) -> GatkResult<Vec<f64>> {
    if haplotypes.is_empty() {
        return Ok(Vec::new());
    }
    validate_read_arrays(
        read_bases,
        read_quals,
        insertion_gop,
        deletion_gop,
        overall_gcp,
    )?;
    let prep = Prep::build(
        read_bases,
        read_quals,
        insertion_gop,
        deletion_gop,
        overall_gcp,
    );
    let mut out = Vec::with_capacity(haplotypes.len());
    WF_SCRATCH.with(|cell| {
        let mut scratch = cell.borrow_mut();
        for &hap in haplotypes {
            if hap.is_empty() {
                out.clear();
                return;
            }
            let sum = rolling_f32::score_one_portable_f32(&prep, read_bases, hap, &mut scratch);
            if !sum.is_finite() || sum < 1e-20 {
                out.push(score_one_rolling_f64(&prep, read_bases, hap, &mut scratch));
            } else {
                let s = sum as f64;
                out.push(s.log10() - rolling_f32::INITIAL_CONDITION_F32_LOG10);
            }
        }
    });
    if out.len() != haplotypes.len() {
        return Err(GatkError::argument("haplotype must be non-empty"));
    }
    Ok(out)
}

/// Which wavefront kernel this host will use for M/I striping.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WavefrontKernel {
    PortableF32,
    Avx2F32,
    NeonF32,
}

pub fn select_wavefront_kernel() -> WavefrontKernel {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if is_x86_feature_detected!("avx2") {
            return WavefrontKernel::Avx2F32;
        }
    }
    #[cfg(target_arch = "aarch64")]
    {
        if std::arch::is_aarch64_feature_detected!("neon") {
            return WavefrontKernel::NeonF32;
        }
    }
    WavefrontKernel::PortableF32
}
