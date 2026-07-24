//! `DepthPerSampleHC` vs FORMAT DP reconciliation scaffold (I-D04).

/// HC depth: sum of AD entries (parity scaffold; Java uses read-depth model).
pub fn depth_per_sample_hc(ad: &[i32]) -> i32 {
    ad.iter().sum()
}

/// FORMAT DP from the same AD vector (parity v1 uses sum for both).
pub fn format_dp_from_ad(ad: &[i32]) -> i32 {
    depth_per_sample_hc(ad)
}
