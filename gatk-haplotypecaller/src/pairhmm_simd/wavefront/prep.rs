//! Amortized read-side PairHMM preparation (transitions + emission priors).

use crate::pairhmm_logless::{logless_match_mismatch_prior, logless_qual_to_trans_probs};

/// Precomputed planes for one read, shared across all haplotypes.
#[derive(Debug, Clone)]
pub struct ReadPrep {
    pub transitions_f64: Vec<[f64; 6]>,
    pub transitions_f32: Vec<[f32; 6]>,
    /// Per read base: (match_p, mismatch_p) as f32.
    pub match_mm_f32: Vec<(f32, f32)>,
    /// Per read base: (match_p, mismatch_p) as f64 (retry / rolling f64).
    pub match_mm_f64: Vec<(f64, f64)>,
}

impl ReadPrep {
    pub fn build(
        read_bases: &[u8],
        read_quals: &[u8],
        insertion_gop: &[u8],
        deletion_gop: &[u8],
        overall_gcp: &[u8],
    ) -> Self {
        let rn = read_bases.len();
        let mut transitions_f64 = vec![[0.0f64; 6]; rn + 1];
        let mut transitions_f32 = vec![[0.0f32; 6]; rn + 1];
        let mut match_mm_f64 = Vec::with_capacity(rn);
        let mut match_mm_f32 = Vec::with_capacity(rn);
        for i in 0..rn {
            let t = logless_qual_to_trans_probs(insertion_gop[i], deletion_gop[i], overall_gcp[i]);
            transitions_f64[i + 1] = t;
            transitions_f32[i + 1] = [
                t[0] as f32,
                t[1] as f32,
                t[2] as f32,
                t[3] as f32,
                t[4] as f32,
                t[5] as f32,
            ];
            let (m, mm) = logless_match_mismatch_prior(read_quals[i]);
            match_mm_f64.push((m, mm));
            match_mm_f32.push((m as f32, mm as f32));
        }
        Self {
            transitions_f64,
            transitions_f32,
            match_mm_f32,
            match_mm_f64,
        }
    }

    #[inline]
    pub fn read_len(&self) -> usize {
        self.match_mm_f64.len()
    }
}
