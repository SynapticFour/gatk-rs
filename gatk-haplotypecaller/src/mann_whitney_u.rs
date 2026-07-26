//! GATK `org.broadinstitute.hellbender.utils.MannWhitneyU` (rank-sum / BaseQRankSum).

use statrs::distribution::{ContinuousCDF, Normal};
use std::cmp::Ordering;
use std::collections::HashMap;

const NORMAL_MEAN: f64 = 0.0;
const NORMAL_SD: f64 = 1.0;

/// Which tail of the Mann–Whitney U distribution to use (one- or two-sided).
/// # Invariants
/// Drives U statistic selection and p-value halving in [`MannWhitneyU::test`].
/// # Ownership
/// [`Copy`] enum.
/// # Mutation
/// Immutable test configuration per call.
/// # Biological assumptions
/// Used to compare two numeric samples (e.g., ref vs alt base qualities).
/// # Java equivalence
/// GATK `MannWhitneyU.TestType`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TestType {
    FirstDominates,
    #[allow(dead_code)] // GATK API surface; callers may flip series order
    SecondDominates,
    TwoSided,
}

/// Mann–Whitney U test statistics (U, z, p, median shift).
/// # Invariants
/// When either input series is empty, `u`, `z`, `p`, and `median_shift` are NaN.
/// # Ownership
/// [`Copy`] result bundle returned from [`MannWhitneyU::test`].
/// # Mutation
/// Immutable statistics snapshot.
/// # Biological assumptions
/// Compares two independent numeric samples (rank-sum style QC annotations).
/// # Java equivalence
/// GATK `MannWhitneyU` result object (`org.broadinstitute.hellbender.utils.MannWhitneyU`).
#[derive(Debug, Clone, Copy)]
pub struct MannWhitneyResult {
    #[allow(dead_code)] // public stats bundle; BaseQ path currently consumes `z`
    pub u: f64,
    pub z: f64,
    #[allow(dead_code)]
    pub p: f64,
    #[allow(dead_code)]
    pub median_shift: f64,
}

struct Rank {
    value: f64,
    rank: f32,
    series: u8,
}

struct RankedData {
    ranks: Vec<Rank>,
    num_of_ties: Vec<i32>,
}

struct TestStatistic {
    u1: f64,
    u2: f64,
    true_u: f64,
    num_of_ties: f64,
}

/// Mann–Whitney U rank-sum test engine (BaseQ rank sum and related annotations).
/// # Invariants
/// `minimum_normal_n` controls normal approximation vs exact permutation path.
/// # Ownership
/// Stateless calculator except for `minimum_normal_n`; input series are copied/sorted internally.
/// # Mutation
/// [`Self::test`] may sort input copies; the engine itself is not mutated.
/// # Biological assumptions
/// Two-sample comparison of per-read or per-base numeric scores at a variant site.
/// # Java equivalence
/// GATK `org.broadinstitute.hellbender.utils.MannWhitneyU`.
pub struct MannWhitneyU {
    minimum_normal_n: usize,
}

impl Default for MannWhitneyU {
    fn default() -> Self {
        Self {
            minimum_normal_n: 10,
        }
    }
}

impl MannWhitneyU {
    pub fn test(
        &self,
        series1: &[f64],
        series2: &[f64],
        which_side: TestType,
    ) -> MannWhitneyResult {
        let mut s1 = series1.to_vec();
        let mut s2 = series2.to_vec();
        let n1 = s1.len();
        let n2 = s2.len();
        if n1 == 0 || n2 == 0 {
            return MannWhitneyResult {
                u: f64::NAN,
                z: f64::NAN,
                p: f64::NAN,
                median_shift: f64::NAN,
            };
        }
        let (u, nties) = match which_side {
            TestType::TwoSided => {
                let stat = self.calculate_u1_and_u2(&mut s1, &mut s2);
                (stat.u1.min(stat.u2), stat.num_of_ties)
            }
            TestType::FirstDominates => {
                let stat = self.calculate_one_sided_u(&mut s1, &mut s2, TestType::FirstDominates);
                (stat.true_u, stat.num_of_ties)
            }
            TestType::SecondDominates => {
                let stat = self.calculate_one_sided_u(&mut s1, &mut s2, TestType::SecondDominates);
                (stat.true_u, stat.num_of_ties)
            }
        };
        let (z, p) = if n1 >= self.minimum_normal_n || n2 >= self.minimum_normal_n {
            let z = self.calculate_z(u, n1, n2, nties, which_side);
            let normal = Normal::new(NORMAL_MEAN, NORMAL_SD).expect("standard normal");
            let mut p = 2.0 * normal.cdf(NORMAL_MEAN + z * NORMAL_SD);
            if which_side != TestType::TwoSided {
                p /= 2.0;
            }
            (z, p)
        } else {
            let p = self.permutation_test(&mut s1, &mut s2, u);
            let normal = Normal::new(NORMAL_MEAN, NORMAL_SD).expect("standard normal");
            let z = normal.inverse_cdf(p);
            (z, p)
        };
        MannWhitneyResult {
            u,
            z,
            p,
            median_shift: (median(&s1) - median(&s2)).abs(),
        }
    }

    fn calculate_rank(&self, series1: &mut [f64], series2: &mut [f64]) -> RankedData {
        series1.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
        series2.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
        let mut ranks = Vec::with_capacity(series1.len() + series2.len());
        let mut i = 0usize;
        let mut j = 0usize;
        while ranks.len() < series1.len() + series2.len() {
            let r = ranks.len() + 1;
            if i >= series1.len() {
                ranks.push(Rank {
                    value: series2[j],
                    rank: r as f32,
                    series: 2,
                });
                j += 1;
            } else if j >= series2.len() {
                ranks.push(Rank {
                    value: series1[i],
                    rank: r as f32,
                    series: 1,
                });
                i += 1;
            } else if series1[i] <= series2[j] {
                ranks.push(Rank {
                    value: series1[i],
                    rank: r as f32,
                    series: 1,
                });
                i += 1;
            } else {
                ranks.push(Rank {
                    value: series2[j],
                    rank: r as f32,
                    series: 2,
                });
                j += 1;
            }
        }
        let mut num_of_ties = Vec::new();
        let mut idx = 0usize;
        while idx < ranks.len() {
            let mut rank_sum = ranks[idx].rank;
            let mut count = 1usize;
            let value = ranks[idx].value;
            let mut j = idx + 1;
            while j < ranks.len() && ranks[j].value == value {
                rank_sum += ranks[j].rank;
                count += 1;
                j += 1;
            }
            if count > 1 {
                let avg = rank_sum / count as f32;
                for r in &mut ranks[idx..idx + count] {
                    r.rank = avg;
                }
                num_of_ties.push(count as i32);
            }
            idx += count;
        }
        RankedData { ranks, num_of_ties }
    }

    fn calculate_u1_and_u2(&self, series1: &mut [f64], series2: &mut [f64]) -> TestStatistic {
        let ranked = self.calculate_rank(series1, series2);
        let length = ranked.ranks.len();
        let nties = transform_ties(length, &ranked.num_of_ties);
        let mut r1 = 0.0f32;
        let mut r2 = 0.0f32;
        for rank in &ranked.ranks {
            if rank.series == 1 {
                r1 += rank.rank;
            } else {
                r2 += rank.rank;
            }
        }
        let n1 = series1.len() as f64;
        let n2 = series2.len() as f64;
        let u1 = f64::from(r1) - (n1 * (n1 + 1.0)) / 2.0;
        let u2 = f64::from(r2) - (n2 * (n2 + 1.0)) / 2.0;
        TestStatistic {
            u1,
            u2,
            true_u: f64::NAN,
            num_of_ties: nties,
        }
    }

    fn calculate_one_sided_u(
        &self,
        series1: &mut [f64],
        series2: &mut [f64],
        which: TestType,
    ) -> TestStatistic {
        let stat = self.calculate_u1_and_u2(series1, series2);
        let true_u = if which == TestType::FirstDominates {
            stat.u1
        } else {
            stat.u2
        };
        TestStatistic {
            u1: stat.u1,
            u2: stat.u2,
            true_u,
            num_of_ties: stat.num_of_ties,
        }
    }

    fn calculate_z(&self, u: f64, n1: usize, n2: usize, nties: f64, which_side: TestType) -> f64 {
        let n1f = n1 as f64;
        let n2f = n2 as f64;
        let m = (n1f * n2f) / 2.0;
        let correction = match which_side {
            TestType::TwoSided => {
                if (u - m) >= 0.0 {
                    0.5
                } else {
                    -0.5
                }
            }
            TestType::FirstDominates => -0.5,
            TestType::SecondDominates => 0.5,
        };
        let correction = if nties == 0.0 { 0.0 } else { correction };
        let sigma = ((n1f * n2f / 12.0)
            * ((n1f + n2f + 1.0) - nties / ((n1f + n2f) * (n1f + n2f - 1.0))))
            .sqrt();
        (u - m - correction) / sigma
    }

    fn permutation_test(&self, series1: &mut [f64], series2: &mut [f64], test_stat_u: f64) -> f64 {
        let n1 = series1.len();
        let n2 = series2.len();
        let ranked = self.calculate_rank(series1, series2);
        let ranks = ranked.ranks;
        let mut first_perm: Vec<u8> = vec![0; n1 + n2];
        for i in n1..n1 + n2 {
            first_perm[i] = 1;
        }
        let perms = get_permutations(&first_perm);
        let mut histo: HashMap<i64, f64> = HashMap::new();
        let mut new_series1 = vec![0.0; n1];
        let mut new_series2 = vec![0.0; n2];
        for perm in perms {
            let mut s1_end = 0usize;
            let mut s2_end = 0usize;
            for (i, &grouping) in perm.iter().enumerate() {
                if grouping == 0 {
                    new_series1[s1_end] = f64::from(ranks[i].rank);
                    s1_end += 1;
                } else {
                    new_series2[s2_end] = f64::from(ranks[i].rank);
                    s2_end += 1;
                }
            }
            let new_u = new_series1.iter().sum::<f64>() - (n1 as f64 * (n1 as f64 + 1.0)) / 2.0;
            let key = (new_u * 2.0).round() as i64;
            *histo.entry(key).or_insert(0.0) += 1.0;
        }
        let test_key = (test_stat_u * 2.0).round() as i64;
        let mut sum_smaller = histo.get(&test_key).copied().unwrap_or(0.0) / 2.0;
        for (key, val) in &histo {
            if *key < test_key {
                sum_smaller += val;
            }
        }
        let total: f64 = histo.values().sum();
        sum_smaller / total
    }
}

fn transform_ties(num_of_ranks: usize, num_of_ties: &[i32]) -> f64 {
    let mut total = 0.0;
    for &count in num_of_ties {
        if count as usize != num_of_ranks {
            total += f64::from(count).powi(3) - f64::from(count);
        }
    }
    total
}

fn median(data: &[f64]) -> f64 {
    let len = data.len();
    let mid = len / 2;
    if len % 2 == 0 {
        (data[mid] + data[mid - 1]) / 2.0
    } else {
        data[mid]
    }
}

fn get_permutations(first_perm: &[u8]) -> Vec<Vec<u8>> {
    let mut temp: Vec<u8> = first_perm.to_vec();
    let mut all = Vec::new();
    // CLONE: needed because owned element into collection.
    all.push(temp.clone());
    loop {
        let mut k: Option<usize> = None;
        for i in (0..temp.len().saturating_sub(1)).rev() {
            if temp[i] < temp[i + 1] {
                k = Some(i);
                break;
            }
        }
        let k = match k {
            Some(v) => v,
            None => break,
        };
        let mut l: Option<usize> = None;
        for i in (k + 1..temp.len()).rev() {
            if temp[k] < temp[i] {
                l = Some(i);
                break;
            }
        }
        let l = l.expect("l exists");
        temp.swap(k, l);
        let mut end = temp.len() - 1;
        let mut begin = k + 1;
        while begin < end {
            temp.swap(begin, end);
            begin += 1;
            end -= 1;
        }
        // CLONE: needed because owned element into collection.
        all.push(temp.clone());
    }
    all
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gatk_baseq_unit_test_vector() {
        let u = MannWhitneyU::default();
        let r = u.test(&[10.0, 20.0], &[50.0, 60.0], TestType::FirstDominates);
        assert!(!r.z.is_nan(), "z={}", r.z);
    }

    #[test]
    fn parity_strand_bias_bqs() {
        let u = MannWhitneyU::default();
        let r = u.test(&[28.0, 29.0], &[30.0, 31.0, 32.0], TestType::FirstDominates);
        assert!(!r.z.is_nan(), "z={}", r.z);
    }
}
