//! Deterministic scenario decoding from arbitrary bytes (libFuzzer / campaign seeds).

use serde::{Deserialize, Serialize};

/// Compact HC differential scenario. Kept small for M4 / CI survivability.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct Scenario {
    pub seed: u64,
    /// Reference length (bp), inclusive of padding.
    pub ref_len: u32,
    pub n_reads: u32,
    pub read_len: u32,
    /// Mean coverage target (approximate).
    pub coverage: u32,
    /// Probability×255 of inserting a 1–indel_max indel event per read.
    pub indel_p: u8,
    pub indel_max: u8,
    /// Probability×255 of soft-clipping 1–softclip_max bases.
    pub softclip_p: u8,
    pub softclip_max: u8,
    /// Probability×255 of emitting an overlapping mate pair.
    pub mate_overlap_p: u8,
    /// Floor MAPQ (0–60).
    pub mapq_min: u8,
    /// MAPQ span above floor.
    pub mapq_span: u8,
    /// Mean Phred base quality (ASCII offset applied at emit).
    pub bq_mean: u8,
    /// BQ jitter (0–20).
    pub bq_jitter: u8,
    /// Plant a SNP near mid-reference.
    pub plant_snp: bool,
    /// Plant a small indel near mid-reference.
    pub plant_indel: bool,
}

impl Default for Scenario {
    fn default() -> Self {
        Self {
            seed: 1,
            ref_len: 200,
            n_reads: 12,
            read_len: 75,
            coverage: 10,
            indel_p: 20,
            indel_max: 3,
            softclip_p: 30,
            softclip_max: 8,
            mate_overlap_p: 40,
            mapq_min: 10,
            mapq_span: 40,
            bq_mean: 30,
            bq_jitter: 8,
            plant_snp: true,
            plant_indel: false,
        }
    }
}

fn take_u8(data: &[u8], i: &mut usize) -> u8 {
    if *i < data.len() {
        let v = data[*i];
        *i += 1;
        v
    } else {
        // Stable pad when fuzz input is short.
        let v = (*i as u8).wrapping_mul(31).wrapping_add(17);
        *i += 1;
        v
    }
}

fn take_u64(data: &[u8], i: &mut usize) -> u64 {
    let mut b = [0u8; 8];
    for slot in &mut b {
        *slot = take_u8(data, i);
    }
    u64::from_le_bytes(b)
}

fn clamp_u32(v: u8, lo: u32, hi: u32) -> u32 {
    lo + (u32::from(v) % (hi - lo + 1))
}

/// Decode a scenario from arbitrary bytes (libFuzzer-compatible).
pub fn scenario_from_bytes(data: &[u8]) -> Scenario {
    let mut i = 0usize;
    let seed = if data.is_empty() {
        0xC0FFEE
    } else {
        take_u64(data, &mut i)
    };
    let mut s = Scenario {
        seed,
        ref_len: clamp_u32(take_u8(data, &mut i), 120, 360),
        n_reads: clamp_u32(take_u8(data, &mut i), 4, 28),
        read_len: clamp_u32(take_u8(data, &mut i), 40, 120),
        coverage: clamp_u32(take_u8(data, &mut i), 4, 20),
        indel_p: take_u8(data, &mut i),
        indel_max: 1 + (take_u8(data, &mut i) % 6),
        softclip_p: take_u8(data, &mut i),
        softclip_max: 1 + (take_u8(data, &mut i) % 12),
        mate_overlap_p: take_u8(data, &mut i),
        mapq_min: take_u8(data, &mut i) % 41,
        mapq_span: 1 + (take_u8(data, &mut i) % 40),
        bq_mean: 15 + (take_u8(data, &mut i) % 26),
        bq_jitter: take_u8(data, &mut i) % 16,
        plant_snp: take_u8(data, &mut i) % 2 == 0,
        plant_indel: take_u8(data, &mut i) % 3 == 0,
    };
    // Keep reads inside reference.
    if s.read_len >= s.ref_len {
        s.read_len = (s.ref_len / 2).max(30);
    }
    s
}

/// Encode a scenario back to a compact seed (for fixture metadata / re-runs).
pub fn scenario_to_seed_bytes(s: &Scenario) -> Vec<u8> {
    let mut out = Vec::with_capacity(24);
    out.extend_from_slice(&s.seed.to_le_bytes());
    out.push((s.ref_len.saturating_sub(120) % 256) as u8);
    out.push((s.n_reads.saturating_sub(4) % 256) as u8);
    out.push((s.read_len.saturating_sub(40) % 256) as u8);
    out.push((s.coverage.saturating_sub(4) % 256) as u8);
    out.push(s.indel_p);
    out.push(s.indel_max.saturating_sub(1));
    out.push(s.softclip_p);
    out.push(s.softclip_max.saturating_sub(1));
    out.push(s.mate_overlap_p);
    out.push(s.mapq_min);
    out.push(s.mapq_span.saturating_sub(1));
    out.push(s.bq_mean.saturating_sub(15));
    out.push(s.bq_jitter);
    out.push(u8::from(s.plant_snp));
    out.push(u8::from(s.plant_indel));
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_bytes_yield_valid_scenario() {
        let s = scenario_from_bytes(&[]);
        assert!(s.ref_len >= 120);
        assert!(s.n_reads >= 4);
        assert!(s.read_len < s.ref_len);
    }

    #[test]
    fn deterministic() {
        let a = scenario_from_bytes(&[1, 2, 3, 4, 5, 6, 7, 8, 9]);
        let b = scenario_from_bytes(&[1, 2, 3, 4, 5, 6, 7, 8, 9]);
        assert_eq!(a, b);
    }
}
