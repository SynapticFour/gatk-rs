//! Runtime PairHMM implementation selection.

use super::pack::{score_haps_logless_packed_f32, score_haps_logless_packed_f64};
use crate::pairhmm_log10::log10_pairhmm_likelihood;
use crate::pairhmm_logless::logless_pairhmm_likelihood;
use gatk_common::{GatkError, GatkResult};

/// User / CLI-selected PairHMM implementation family.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairHmmImpl {
    /// Exact scalar `Log10PairHMM` (F.1/F.2 dumps; current production default).
    Log10PairHmm,
    /// Scalar linear-space `LoglessPairHMM`.
    LoglessPairHmm,
    /// Best host SIMD (f64 lanes; optional f32+retry via `SimdF32`).
    Simd,
    /// SIMD f32 with f64 retry on underflow.
    SimdF32,
    /// Resolve at score time: SIMD if available, else Logless, else Log10.
    /// Until validation gates close, production config still defaults to [`Self::Log10PairHmm`].
    FastestAvailable,
}

/// Concrete backend chosen after runtime feature detection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairHmmBackend {
    Log10Scalar,
    LoglessScalar,
    PackedF64,
    Avx2F64,
    NeonF64,
    PackedF32Retry,
}

impl PairHmmImpl {
    pub fn label(self) -> &'static str {
        match self {
            Self::Log10PairHmm => "LOG10_PAIRHMM",
            Self::LoglessPairHmm => "LOGLESS_HMM",
            Self::Simd => "SIMD",
            Self::SimdF32 => "SIMD_F32",
            Self::FastestAvailable => "FASTEST_AVAILABLE",
        }
    }
}

/// Parse GATK-shaped `--pair-hmm` / `--pairhmm-impl` values (case-insensitive).
pub fn parse_pair_hmm_impl(s: &str) -> GatkResult<PairHmmImpl> {
    let u = s.trim().to_ascii_uppercase();
    match u.as_str() {
        "LOG10_PAIRHMM" | "LOG10" | "SCALAR" => Ok(PairHmmImpl::Log10PairHmm),
        "LOGLESS_HMM" | "LOGLESS" => Ok(PairHmmImpl::LoglessPairHmm),
        "AVX" | "SIMD" => Ok(PairHmmImpl::Simd),
        "SIMD_F32" | "AVX_F32" => Ok(PairHmmImpl::SimdF32),
        "FASTEST_AVAILABLE" | "FASTEST" => Ok(PairHmmImpl::FastestAvailable),
        _ => Err(GatkError::argument(format!(
            "unknown --pair-hmm value '{s}' (expected LOG10_PAIRHMM|LOGLESS_HMM|AVX|SIMD|SIMD_F32|FASTEST_AVAILABLE)"
        ))),
    }
}

/// Resolve [`PairHmmImpl`] to a concrete backend for this host.
pub fn resolve_pair_hmm_impl(imp: PairHmmImpl) -> PairHmmBackend {
    match imp {
        PairHmmImpl::Log10PairHmm => PairHmmBackend::Log10Scalar,
        PairHmmImpl::LoglessPairHmm => PairHmmBackend::LoglessScalar,
        PairHmmImpl::SimdF32 => PairHmmBackend::PackedF32Retry,
        PairHmmImpl::Simd => best_simd_backend(),
        PairHmmImpl::FastestAvailable => {
            // Until GIAB sign-off, callers should still pass Log10 as production default.
            // When FastestAvailable is explicitly requested, prefer SIMD → Logless.
            let simd = best_simd_backend();
            if matches!(
                simd,
                PairHmmBackend::Avx2F64 | PairHmmBackend::NeonF64 | PairHmmBackend::PackedF64
            ) {
                simd
            } else {
                PairHmmBackend::LoglessScalar
            }
        }
    }
}

pub fn best_simd_backend() -> PairHmmBackend {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if is_x86_feature_detected!("avx2") {
            return PairHmmBackend::Avx2F64;
        }
    }
    #[cfg(target_arch = "aarch64")]
    {
        if std::arch::is_aarch64_feature_detected!("neon") {
            return PairHmmBackend::NeonF64;
        }
    }
    PairHmmBackend::PackedF64
}

/// Score one read against many haplotypes with the selected implementation.
pub fn score_read_haps_logless(
    backend: PairHmmBackend,
    read_bases: &[u8],
    read_quals: &[u8],
    haplotypes: &[&[u8]],
    insertion_gop: &[u8],
    deletion_gop: &[u8],
    overall_gcp: &[u8],
) -> GatkResult<Vec<f64>> {
    match backend {
        PairHmmBackend::Log10Scalar => haplotypes
            .iter()
            .map(|hap| {
                log10_pairhmm_likelihood(
                    read_bases,
                    read_quals,
                    hap,
                    insertion_gop,
                    deletion_gop,
                    overall_gcp,
                )
            })
            .collect(),
        PairHmmBackend::LoglessScalar => haplotypes
            .iter()
            .map(|hap| {
                logless_pairhmm_likelihood(
                    read_bases,
                    read_quals,
                    hap,
                    insertion_gop,
                    deletion_gop,
                    overall_gcp,
                )
            })
            .collect(),
        PairHmmBackend::PackedF64 => score_haps_logless_packed_f64(
            read_bases,
            read_quals,
            haplotypes,
            insertion_gop,
            deletion_gop,
            overall_gcp,
        ),
        PairHmmBackend::PackedF32Retry => score_haps_logless_packed_f32(
            read_bases,
            read_quals,
            haplotypes,
            insertion_gop,
            deletion_gop,
            overall_gcp,
        ),
        PairHmmBackend::Avx2F64 => {
            #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
            {
                return super::avx2::score_haps_avx2_f64(
                    read_bases,
                    read_quals,
                    haplotypes,
                    insertion_gop,
                    deletion_gop,
                    overall_gcp,
                );
            }
            #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
            {
                score_haps_logless_packed_f64(
                    read_bases,
                    read_quals,
                    haplotypes,
                    insertion_gop,
                    deletion_gop,
                    overall_gcp,
                )
            }
        }
        PairHmmBackend::NeonF64 => {
            #[cfg(target_arch = "aarch64")]
            {
                super::neon::score_haps_neon_f64(
                    read_bases,
                    read_quals,
                    haplotypes,
                    insertion_gop,
                    deletion_gop,
                    overall_gcp,
                )
            }
            #[cfg(not(target_arch = "aarch64"))]
            {
                score_haps_logless_packed_f64(
                    read_bases,
                    read_quals,
                    haplotypes,
                    insertion_gop,
                    deletion_gop,
                    overall_gcp,
                )
            }
        }
    }
}
