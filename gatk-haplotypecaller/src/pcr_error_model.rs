//! GATK `PairHMMLikelihoodCalculationEngine` PCR indel error model.

use std::sync::OnceLock;

/// GATK 4.4 `PairHMMLikelihoodCalculationEngine.MIN_ADJUSTED_QSCORE`.
const MIN_ADJUSTED_QSCORE: i32 = 10;
/// GATK 4.4 `PairHMMLikelihoodCalculationEngine.INITIAL_QSCORE`.
const INITIAL_QSCORE: f64 = 40.0;
/// GATK `ReadLikelihoodCalculationEngine.MAX_REPEAT_LENGTH`. Clip on **unit count**.
pub const MAX_REPEAT_LENGTH: usize = 20;
/// GATK `ReadLikelihoodCalculationEngine.MAX_STR_UNIT_LENGTH`.
pub const MAX_STR_UNIT_LENGTH: usize = 8;

fn pcr_cache_conservative() -> &'static [u8; MAX_REPEAT_LENGTH + 1] {
    static CACHE: OnceLock<[u8; MAX_REPEAT_LENGTH + 1]> = OnceLock::new();
    CACHE.get_or_init(|| build_pcr_cache(3.0))
}

fn pcr_cache_aggressive() -> &'static [u8; MAX_REPEAT_LENGTH + 1] {
    static CACHE: OnceLock<[u8; MAX_REPEAT_LENGTH + 1]> = OnceLock::new();
    CACHE.get_or_init(|| build_pcr_cache(2.0))
}

/// GATK `PCRErrorModel` for HC defaults.
/// # Invariants
/// [`Self::None`] leaves insertion/deletion quals unchanged; other variants apply
/// tandem-repeat caps (unit width 1..=8, including homopolymers).
/// Conservative and Aggressive differ only by [`Self::rate_factor`] (3 vs 2).
/// # Ownership
/// [`Copy`] enum passed into [`apply_pcr_error_model`]; mutates qual arrays in place.
/// # Mutation
/// Immutable model tag; target `ins_quals` / `del_quals` slices are updated per read base.
/// # Biological assumptions
/// PCR stutter scales with tandem-repeat unit count at indel positions.
/// # Java equivalence
/// GATK `PCRErrorModel` in `PairHMMLikelihoodCalculationEngine`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PcrErrorModel {
    None,
    /// `PCRErrorModel.CONSERVATIVE` (rate factor 3).
    Conservative,
    /// `PCRErrorModel.AGGRESSIVE` (rate factor 2).
    Aggressive,
}

impl PcrErrorModel {
    pub fn rate_factor(self) -> Option<f64> {
        match self {
            Self::None => None,
            Self::Conservative => Some(3.0),
            Self::Aggressive => Some(2.0),
        }
    }
}

/// GATK 4.4 `MathUtils.fastRound`.
fn fast_round(d: f64) -> i32 {
    if d > 0.0 {
        (d + 0.5) as i32
    } else {
        (d - 0.5) as i32
    }
}

/// GATK 4.4 `getErrorModelAdjustedQual`.
pub fn error_model_adjusted_qual(repeat_length: usize, rate_factor: f64) -> u8 {
    let exp_term = (repeat_length as f64 / (rate_factor * std::f64::consts::PI)).exp();
    let q = INITIAL_QSCORE - exp_term + 1.0;
    fast_round(q).max(MIN_ADJUSTED_QSCORE) as u8
}

fn build_pcr_cache(rate_factor: f64) -> [u8; MAX_REPEAT_LENGTH + 1] {
    let mut cache = [0u8; MAX_REPEAT_LENGTH + 1];
    for (i, slot) in cache.iter_mut().enumerate() {
        *slot = error_model_adjusted_qual(i, rate_factor);
    }
    cache
}

fn equal_range(a: &[u8], a_off: usize, b: &[u8], b_off: usize, len: usize) -> bool {
    a.get(a_off..a_off + len) == b.get(b_off..b_off + len)
}

/// GATK 4.4 `GATKVariantContextUtils.findNumberOfRepetitions` (subarray form).
/// Counts **whole units** only; a leftover partial unit is not counted.
fn find_number_of_repetitions(
    unit: &[u8],
    unit_off: usize,
    unit_len: usize,
    test: &[u8],
    test_off: usize,
    test_len: usize,
    leading: bool,
) -> usize {
    if unit_len == 0 || test_len == 0 {
        return 0;
    }
    let length_difference = test_len as isize - unit_len as isize;
    if leading {
        let mut n = 0usize;
        let mut start = 0isize;
        while start <= length_difference {
            if equal_range(test, start as usize + test_off, unit, unit_off, unit_len) {
                n += 1;
                start += unit_len as isize;
            } else {
                return n;
            }
        }
        n
    } else {
        let mut n = 0usize;
        let mut start = length_difference;
        while start >= 0 {
            if equal_range(test, start as usize + test_off, unit, unit_off, unit_len) {
                n += 1;
                start -= unit_len as isize;
            } else {
                return n;
            }
        }
        n
    }
}

/// GATK 4.4 `ReadLikelihoodCalculationEngine.findTandemRepeatUnits`.
///
/// Returns `(repeat unit, unit count)`. The count is the PCR cache index.
/// Widths are tried **1..=8**; the first width with count `> 1` wins (not longest).
/// Backward is trailing copies in `bases[0..=offset]`; forward is leading copies
/// in the suffix. Mismatching FW/BW units keep the **forward** unit and re-count
/// trailing FW copies in the prefix. Final count is clipped to [`MAX_REPEAT_LENGTH`].
///
/// The last-base PCR skip lives in [`apply_pcr_error_model`], not here.
pub fn find_tandem_repeat_units(bases: &[u8], offset: usize) -> (Vec<u8>, usize) {
    if bases.is_empty() || offset >= bases.len() {
        return (Vec::new(), 0);
    }
    let mut max_bw = 0usize;
    let mut best_bw: Vec<u8> = vec![bases[offset]];
    for str_len in 1..=MAX_STR_UNIT_LENGTH {
        if offset + 1 < str_len {
            break;
        }
        max_bw = find_number_of_repetitions(
            bases,
            offset + 1 - str_len,
            str_len,
            bases,
            0,
            offset + 1,
            false,
        );
        if max_bw > 1 {
            best_bw = bases[offset + 1 - str_len..=offset].to_vec();
            break;
        }
    }
    let mut best_unit = best_bw.clone();
    let mut max_rl = max_bw;

    if offset < bases.len() - 1 {
        let mut best_fw: Vec<u8> = vec![bases[offset + 1]];
        let mut max_fw = 0usize;
        for str_len in 1..=MAX_STR_UNIT_LENGTH {
            if offset + str_len + 1 > bases.len() {
                break;
            }
            max_fw = find_number_of_repetitions(
                bases,
                offset + 1,
                str_len,
                bases,
                offset + 1,
                bases.len() - offset - 1,
                true,
            );
            if max_fw > 1 {
                best_fw = bases[offset + 1..offset + 1 + str_len].to_vec();
                break;
            }
        }
        if best_fw == best_bw {
            max_rl = max_bw + max_fw;
            best_unit = best_fw;
        } else {
            let test_len = offset + 1;
            max_bw =
                find_number_of_repetitions(&best_fw, 0, best_fw.len(), bases, 0, test_len, false);
            max_rl = max_fw + max_bw;
            best_unit = best_fw;
        }
    }
    if max_rl > MAX_REPEAT_LENGTH {
        max_rl = MAX_REPEAT_LENGTH;
    }
    (best_unit, max_rl)
}

/// PCR cache index: unit count from [`find_tandem_repeat_units`].
pub fn tandem_repeat_units(bases: &[u8], i: usize) -> usize {
    find_tandem_repeat_units(bases, i).1
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

    /// 6R.72: GATK 4.4.0.0 `getErrorModelAdjustedQual` CONSERVATIVE cache[0] is Q40.
    #[test]
    fn conservative_repeat_zero_is_java_4_4_q40() {
        let rust = error_model_adjusted_qual(0, PcrErrorModel::Conservative.rate_factor().unwrap());
        let java = java_44_error_model_adjusted_qual(0, 3.0);
        assert_eq!(rust, 40, "Rust CONSERVATIVE cache[0]");
        assert_eq!(java, 40, "Java 4.4 CONSERVATIVE cache[0]");
        assert_eq!(rust, java);
    }

    #[test]
    fn find_tandem_repeat_units_canonical_prefix() {
        let (u, n) = find_tandem_repeat_units(b"TAAGAAAA", 0);
        assert_eq!(u, b"A");
        assert_eq!(n, 2);
        assert_eq!(tandem_repeat_units(b"TAAGAAAA", 0), 2);
    }

    /// GATK 4.4.0.0 `PairHMMLikelihoodCalculationEngine.getErrorModelAdjustedQual`.
    fn java_44_error_model_adjusted_qual(repeat_length: usize, rate_factor: f64) -> u8 {
        let q = 40.0 - (repeat_length as f64 / (rate_factor * std::f64::consts::PI)).exp() + 1.0;
        let rounded = if q > 0.0 {
            (q + 0.5) as i32
        } else {
            (q - 0.5) as i32
        };
        rounded.max(10) as u8
    }
}
