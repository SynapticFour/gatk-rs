//! Property tests for `normalize_from_log10_to_linear_space` (Rust-native R2).
use gatk_haplotypecaller::normalize_from_log10_to_linear_space;
use proptest::prelude::*;

fn argmax(xs: &[f64]) -> usize {
    xs.iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.total_cmp(b))
        .map(|(i, _)| i)
        .unwrap_or(0)
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(128))]

    /// Linear probabilities sum to 1, stay non-negative, and preserve argmax.
    #[test]
    fn normalize_log10_to_linear_preserves_argmax_and_sums(
        vals in prop::collection::vec(-80.0f64..0.0, 2..8)
    ) {
        prop_assume!(vals.iter().all(|x| x.is_finite()));
        let mut sorted = vals.clone();
        sorted.sort_by(|a, b| b.total_cmp(a));
        prop_assume!(sorted[0] - sorted[1] >= 1e-6);
        let before = argmax(&vals);
        let linear = normalize_from_log10_to_linear_space(&vals);
        prop_assert_eq!(linear.len(), vals.len());
        prop_assert!(linear.iter().all(|x| x.is_finite() && *x >= 0.0));
        let sum: f64 = linear.iter().sum();
        prop_assert!(
            (sum - 1.0).abs() < 1e-9,
            "expected sum≈1, got {sum}"
        );
        prop_assert_eq!(argmax(&linear), before);
    }
}
