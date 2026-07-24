//! Apache Commons Math `Well19937c` + `MathArrays.shuffle` used by GATK
//! `MathUtils.sampleIndicesWithoutReplacement` / `AlleleBiasedDownsamplingUtils`.

/// GATK `Utils.GATK_RANDOM_SEED` — also seeds `RandomDataGenerator(Well19937c)`.
pub const GATK_WELL19937C_SEED: i64 = 47_382_911;

/// `org.apache.commons.math3.random.Well19937c` (WELL19937c, K=19937).
/// # Invariants
/// Internal state matches Apache Commons Math WELL19937c with GATK seed [`GATK_WELL19937C_SEED`] when reset.
/// Used for deterministic without-replacement sampling / shuffles.
/// # Ownership
/// Owns large internal state vectors; cloneable for independent streams.
/// # Mutation
/// Advances internal state on each random draw / shuffle.
/// # Biological assumptions
/// None — RNG for allele-biased downsampling parity.
/// # Java equivalence
/// Apache Commons Math `Well19937c` via GATK `MathUtils` / `AlleleBiasedDownsamplingUtils`.
#[derive(Debug, Clone)]
pub struct Well19937c {
    v: Vec<i32>,
    index: usize,
    i_rm1: Vec<usize>,
    i_rm2: Vec<usize>,
    i1: Vec<usize>,
    i2: Vec<usize>,
    i3: Vec<usize>,
}

impl Well19937c {
    pub fn new(seed: i64) -> Self {
        Self::from_seed_array(&[(seed >> 32) as i32, seed as i32])
    }

    pub fn reset_gatk_default() -> Self {
        Self::new(GATK_WELL19937C_SEED)
    }

    fn from_seed_array(seed: &[i32]) -> Self {
        const K: usize = 19937;
        const M1: usize = 70;
        const M2: usize = 179;
        const M3: usize = 449;
        const W: usize = 32;
        let r = K.div_ceil(W);
        let mut v = vec![0i32; r];
        let n = seed.len().min(v.len());
        v[..n].copy_from_slice(&seed[..n]);
        if seed.len() < v.len() {
            for i in n..v.len() {
                let l = v[i - seed.len()] as i64;
                v[i] = ((1812433253i64 * (l ^ (l >> 30)) + i as i64) & 0xffff_ffff) as i32;
            }
        }
        let mut i_rm1 = vec![0usize; r];
        let mut i_rm2 = vec![0usize; r];
        let mut i1 = vec![0usize; r];
        let mut i2 = vec![0usize; r];
        let mut i3 = vec![0usize; r];
        for j in 0..r {
            i_rm1[j] = (j + r - 1) % r;
            i_rm2[j] = (j + r - 2) % r;
            i1[j] = (j + M1) % r;
            i2[j] = (j + M2) % r;
            i3[j] = (j + M3) % r;
        }
        Self {
            v,
            index: 0,
            i_rm1,
            i_rm2,
            i1,
            i2,
            i3,
        }
    }

    fn next_bits(&mut self, bits: u32) -> u32 {
        let index_rm1 = self.i_rm1[self.index];
        let index_rm2 = self.i_rm2[self.index];

        let v0 = self.v[self.index];
        let v_m1 = self.v[self.i1[self.index]];
        let v_m2 = self.v[self.i2[self.index]];
        let v_m3 = self.v[self.i3[self.index]];

        let z0 = (self.v[index_rm1] & i32::MIN) ^ (self.v[index_rm2] & i32::MAX);
        let z1 = (v0 ^ v0.wrapping_shl(25)) ^ (v_m1 ^ (v_m1 as u32 >> 27) as i32);
        let z2 = (v_m2 as u32 >> 9) as i32 ^ (v_m3 ^ (v_m3 as u32 >> 1) as i32);
        let z3 = z1 ^ z2;
        let mut z4 = z0
            ^ (z1 ^ z1.wrapping_shl(9))
            ^ (z2 ^ z2.wrapping_shl(21))
            ^ (z3 ^ (z3 as u32 >> 21) as i32);

        self.v[self.index] = z3;
        self.v[index_rm1] = z4;
        self.v[index_rm2] &= i32::MIN;
        self.index = index_rm1;

        z4 ^= z4.wrapping_shl(7) & 0xe46e1700u32 as i32;
        z4 ^= z4.wrapping_shl(15) & 0x9b868000u32 as i32;

        (z4 as u32) >> (32 - bits)
    }

    /// `RandomGenerator.nextInt(bound)` — Commons Math `BitsStreamGenerator`.
    pub fn next_int(&mut self, bound: u32) -> u32 {
        assert!(bound > 0, "bound must be positive");
        loop {
            let bits = self.next_bits(31);
            let val = bits % bound;
            if (bits as i64) - (val as i64) + (bound as i64) > 0 {
                return val;
            }
        }
    }

    /// `MathArrays.shuffle` (TAIL): at index `j`, swap with `k` uniform in `[0, j]` inclusive.
    pub fn shuffle_usize(&mut self, list: &mut [usize]) {
        for j in (1..list.len()).rev() {
            let k = self.next_int((j + 1) as u32) as usize;
            list.swap(j, k);
        }
    }

    /// `MathUtils.sampleIndicesWithoutReplacement` → `RandomDataGenerator.nextPermutation`.
    pub fn sample_indices_without_replacement(&mut self, n: usize, k: usize) -> Vec<usize> {
        assert!(k > 0 && k <= n, "invalid sample size");
        let mut index: Vec<usize> = (0..n).collect();
        self.shuffle_usize(&mut index);
        index.truncate(k);
        index
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn permutation_deterministic_at_gatk_seed() {
        let mut rng = Well19937c::reset_gatk_default();
        let p = rng.sample_indices_without_replacement(10, 3);
        assert_eq!(p.len(), 3);
        let mut rng2 = Well19937c::reset_gatk_default();
        assert_eq!(p, rng2.sample_indices_without_replacement(10, 3));
    }
}
