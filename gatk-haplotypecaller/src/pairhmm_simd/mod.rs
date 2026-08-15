//! SIMD / pairs-in-lanes Logless PairHMM (runtime feature dispatch).
//! Parallelization axis: **independent haplotypes in SIMD lanes** for one shared
//! read (Endeavor-style), not GKL anti-diagonal striping. Scalar
//! [`crate::pairhmm_logless`] remains the correctness reference.

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
mod avx2;
mod dispatch;
#[cfg(target_arch = "aarch64")]
mod neon;
mod pack;

pub use dispatch::{
    best_simd_backend, parse_pair_hmm_impl, resolve_pair_hmm_impl, score_read_haps_logless,
    PairHmmBackend, PairHmmImpl,
};
pub use pack::{score_haps_logless_packed_f32, score_haps_logless_packed_f64};

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
pub use avx2::release_pairhmm_avx2_tls_scratch;
#[cfg(target_arch = "aarch64")]
pub use neon::{release_pairhmm_neon_tls_scratch, take_neon_pack_stats};

/// Release SIMD PairHMM TLS scratch when present.
pub fn release_pairhmm_simd_tls_scratch() {
    #[cfg(target_arch = "aarch64")]
    release_pairhmm_neon_tls_scratch();
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    release_pairhmm_avx2_tls_scratch();
}
