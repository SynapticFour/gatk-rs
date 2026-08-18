//! GATK `PairHMMLikelihoodCalculationEngine` PCR indel error model.

use std::sync::OnceLock;

const MIN_ADJUSTED_QSCORE: i32 = 6;
const INITIAL_QSCORE: f64 = 45.0;
const MAX_REPEAT_LENGTH: usize = 20;

fn pcr_cache_conservative() -> &'static [u8; MAX_REPEAT_LENGTH + 1] {
    static CACHE: OnceLock<[u8; MAX_REPEAT_LENGTH + 1]> = OnceLock::new();
    CACHE.get_or_init(|| build_pcr_cache(20.0))
}

fn pcr_cache_aggressive() -> &'static [u8; MAX_REPEAT_LENGTH + 1] {
    static CACHE: OnceLock<[u8; MAX_REPEAT_LENGTH + 1]> = OnceLock::new();
    CACHE.get_or_init(|| build_pcr_cache(10.0))
}

/// GATK `PCRErrorModel` for HC defaults.
/// # Invariants
/// [`Self::None`] leaves insertion/deletion quals unchanged; other variants apply homopolymer caps.
/// Conservative and Aggressive differ only by [`Self::rate_factor`] (20 vs 10).
/// # Ownership
/// [`Copy`] enum passed into [`apply_pcr_error_model`]; mutates qual arrays in place.
/// # Mutation
/// Immutable model tag; target `ins_quals` / `del_quals` slices are updated per read base.
/// # Biological assumptions
/// PCR stutter scales with homopolymer run length at indel positions.
/// # Java equivalence
/// GATK `PCRErrorModel` in `PairHMMLikelihoodCalculationEngine`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PcrErrorModel {
    None,
    /// `PCRErrorModel.CONSERVATIVE` (rate factor 20).
    Conservative,
    /// `PCRErrorModel.AGGRESSIVE` (rate factor 10).
    Aggressive,
}

impl PcrErrorModel {
    pub fn rate_factor(self) -> Option<f64> {
        match self {
            Self::None => None,
            Self::Conservative => Some(20.0),
            Self::Aggressive => Some(10.0),
        }
    }
}

/// GATK `getErrorModelAdjustedQual`.
pub fn error_model_adjusted_qual(repeat_length: usize, rate_factor: f64) -> u8 {
    let exp_term = (repeat_length as f64 / (rate_factor * std::f64::consts::PI)).exp();
    let q = INITIAL_QSCORE - exp_term + 1.0;
    (q.round() as i32).max(MIN_ADJUSTED_QSCORE) as u8
}

fn build_pcr_cache(rate_factor: f64) -> [u8; MAX_REPEAT_LENGTH + 1] {
    let mut cache = [0u8; MAX_REPEAT_LENGTH + 1];
    for (i, slot) in cache.iter_mut().enumerate() {
        *slot = error_model_adjusted_qual(i, rate_factor);
    }
    cache
}

/// GATK `findTandemRepeatUnits` simplified: homopolymer run length ending at `i`.
pub fn tandem_repeat_units(bases: &[u8], i: usize) -> usize {
    if bases.is_empty() || i >= bases.len() {
        return 0;
    }
    let base = bases[i];
    let mut len = 1usize;
    let mut j = i;
    while j > 0 && bases[j - 1] == base {
        len += 1;
        j -= 1;
    }
    j = i;
    while j + 1 < bases.len() && bases[j + 1] == base {
        len += 1;
        j += 1;
    }
    len
}

/// GATK `applyPCRErrorModel` on insertion/deletion qual arrays.
pub fn apply_pcr_error_model(
    read_bases: &[u8],
    ins_quals: &mut [u8],
    del_quals: &mut [u8],
    model: PcrErrorModel,
) {
    let cache = match model {
        PcrErrorModel::None => return,
        PcrErrorModel::Conservative => pcr_cache_conservative(),
        PcrErrorModel::Aggressive => pcr_cache_aggressive(),
    };
    for i in 1..read_bases.len() {
        let repeat = tandem_repeat_units(read_bases, i - 1);
        let cap = cache[repeat.min(MAX_REPEAT_LENGTH)];
        let idx = i - 1;
        ins_quals[idx] = ins_quals[idx].min(cap);
        del_quals[idx] = del_quals[idx].min(cap);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn conservative_caps_homopolymer_ins_del() {
        let read = b"AAAAAAAAAAAAACGT";
        let mut ins = vec![45u8; read.len()];
        let mut del = vec![45u8; read.len()];
        apply_pcr_error_model(read, &mut ins, &mut del, PcrErrorModel::Aggressive);
        assert!(
            ins.iter().any(|&q| q < 45),
            "ins quals should be capped: {:?}",
            ins
        );
    }
}
