//! Property tests for Java PL→GL round-trip used in AF / emit (Rust-native R1).
use gatk_haplotypecaller::hc_genotyping_engine::gl_for_java_af_calculation;
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

    /// After Java PL round-trip scaling, the best genotype index is unchanged and
    /// the max log10-GL is near zero (stability scale).
    #[test]
    fn java_af_gl_roundtrip_preserves_argmax_and_scales(
        g0 in -80.0f64..0.0,
        g1 in -80.0f64..0.0,
        g2 in -80.0f64..0.0,
    ) {
        let raw = [g0, g1, g2];
        // PL round-trip is integer-phred; near-ties can flip argmax — require a clear winner.
        let mut sorted = raw;
        sorted.sort_by(|a, b| b.total_cmp(a));
        prop_assume!(sorted[0] - sorted[1] >= 0.05);
        let before = argmax(&raw);
        let rt = gl_for_java_af_calculation(&raw);
        prop_assume!(rt.len() == 3);
        prop_assume!(rt.iter().all(|x| x.is_finite()));
        prop_assert_eq!(argmax(&rt), before);
        let max_v = rt.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        prop_assert!(
            (max_v - 0.0).abs() < 1e-6,
            "expected max GL ~ 0 after scale, got {max_v}"
        );
        let rt2 = gl_for_java_af_calculation(&rt);
        prop_assert_eq!(argmax(&rt2), before);
        for (a, b) in rt.iter().zip(rt2.iter()) {
            prop_assert!((a - b).abs() < 1e-9, "second pass not idempotent: {a} vs {b}");
        }
    }
}
