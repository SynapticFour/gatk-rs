//! PairHMM workload counters for production profiling.

use std::collections::BTreeMap;

#[derive(Debug, Clone, Default)]
pub struct PairHmmCallSample {
    pub reads: u64,
    pub haplotypes: u64,
    pub read_len_sum: u64,
    pub hap_len_sum: u64,
    pub read_lens: Vec<u32>,
    pub hap_lens: Vec<u32>,
    /// SIMD pack hits (NEON pack2 / AVX2 pack4 units).
    pub simd_packs: u64,
    /// Haplotypes scored via hapStartIndex prefix reuse path.
    pub prefix_reuse_haps: u64,
    /// Scalar leftover singles.
    pub leftover_haps: u64,
    /// DP cells evaluated (approximate).
    pub dp_cells_evaluated: u64,
    /// DP cells skipped by hapStartIndex (approximate).
    pub dp_cells_avoided_prefix: u64,
    pub wall_ns: u64,
}

#[derive(Debug, Default)]
pub struct PairHmmAgg {
    pub calls: u64,
    pub reads_scored: u64,
    pub haplotypes_scored: u64,
    pub read_hap_pairs: u64,
    pub read_len_sum: u64,
    pub hap_len_sum: u64,
    pub simd_packs: u64,
    pub prefix_reuse_haps: u64,
    pub leftover_haps: u64,
    pub dp_cells_evaluated: u64,
    pub dp_cells_avoided_prefix: u64,
    pub wall_ns: u64,
    pub read_len_hist: BTreeMap<u32, u64>,
    pub hap_len_hist: BTreeMap<u32, u64>,
}

impl PairHmmAgg {
    pub fn add(&mut self, s: PairHmmCallSample) {
        self.calls += 1;
        self.reads_scored += s.reads;
        self.haplotypes_scored += s.haplotypes;
        self.read_hap_pairs += s.reads.saturating_mul(s.haplotypes);
        self.read_len_sum += s.read_len_sum;
        self.hap_len_sum += s.hap_len_sum;
        self.simd_packs += s.simd_packs;
        self.prefix_reuse_haps += s.prefix_reuse_haps;
        self.leftover_haps += s.leftover_haps;
        self.dp_cells_evaluated += s.dp_cells_evaluated;
        self.dp_cells_avoided_prefix += s.dp_cells_avoided_prefix;
        self.wall_ns += s.wall_ns;
        for len in s.read_lens {
            let bucket = bucket_len(len);
            *self.read_len_hist.entry(bucket).or_default() += 1;
        }
        for len in s.hap_lens {
            let bucket = bucket_len(len);
            *self.hap_len_hist.entry(bucket).or_default() += 1;
        }
    }

    pub fn mean_haps_per_read(&self) -> f64 {
        if self.reads_scored == 0 {
            0.0
        } else {
            self.haplotypes_scored as f64 / self.reads_scored as f64
        }
    }

    pub fn simd_occupancy(&self) -> PairHmmOccupancy {
        let pack_lanes = simd_pack_lanes();
        let pack_hap_est = self.simd_packs.saturating_mul(pack_lanes);
        let denom = self
            .prefix_reuse_haps
            .saturating_add(self.leftover_haps)
            .saturating_add(pack_hap_est)
            .max(1);
        PairHmmOccupancy {
            pack_units: self.simd_packs,
            pack_lanes,
            pack_hap_est,
            pack_occupancy_pct: 100.0 * pack_hap_est as f64 / denom as f64,
            prefix_reuse_pct: 100.0 * self.prefix_reuse_haps as f64 / denom as f64,
            leftover_pct: 100.0 * self.leftover_haps as f64 / denom as f64,
            hap_oriented_denom: denom,
        }
    }
}

fn simd_pack_lanes() -> u64 {
    #[cfg(target_arch = "aarch64")]
    {
        2
    }
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        4
    }
    #[cfg(not(any(target_arch = "aarch64", target_arch = "x86", target_arch = "x86_64")))]
    {
        2
    }
}

#[derive(Debug, Clone)]
pub struct PairHmmOccupancy {
    pub pack_units: u64,
    pub pack_lanes: u64,
    pub pack_hap_est: u64,
    pub pack_occupancy_pct: f64,
    pub prefix_reuse_pct: f64,
    pub leftover_pct: f64,
    pub hap_oriented_denom: u64,
}

fn bucket_len(len: u32) -> u32 {
    // 25 bp buckets for readable histograms.
    (len / 25) * 25
}
