//! SIMD / pairs-in-lanes Logless PairHMM (runtime feature dispatch).
//! Parallelization axis: **independent haplotypes in SIMD lanes** for one shared
//! read (Endeavor-style), not GKL anti-diagonal striping. Scalar
//! [`crate::pairhmm_logless`] remains the correctness reference.
//!
//! Architecture C (opt-in): [`wavefront`] — rolling-row f32 + column striping.

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
mod avx2;
mod dispatch;
#[cfg(target_arch = "aarch64")]
mod neon;
mod pack;
pub mod wavefront;

pub use dispatch::{
    best_simd_backend, parse_pair_hmm_impl, resolve_pair_hmm_impl, score_read_haps_logless,
    PairHmmBackend, PairHmmImpl,
};
pub use pack::{
    score_haps_logless_packed_f32, score_haps_logless_packed_f64, take_pack_dp_cell_stats,
};
pub use wavefront::{
    score_haps_wavefront_f32, score_haps_wavefront_portable_f32, score_haps_wavefront_rolling_f64,
    select_wavefront_kernel, ReadPrep, WavefrontKernel,
};

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
pub use avx2::{release_pairhmm_avx2_tls_scratch, take_avx2_pack_stats};
#[cfg(target_arch = "aarch64")]
pub use neon::{release_pairhmm_neon_tls_scratch, take_neon_pack_stats};

/// Release SIMD PairHMM TLS scratch when present.
pub fn release_pairhmm_simd_tls_scratch() {
    #[cfg(target_arch = "aarch64")]
    release_pairhmm_neon_tls_scratch();
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    release_pairhmm_avx2_tls_scratch();
    wavefront::release_wavefront_tls_scratch();
}
