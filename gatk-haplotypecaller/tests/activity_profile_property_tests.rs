//! Property-style tests for activity profile boundaries (Phase 4, step 60).
use gatk_haplotypecaller::{
    ActivityProfileState, BandPassActivityProfile, BandPassActivityProfileParams, PositiveSigma,
};
use proptest::prelude::*;

fn no_smoothing_profile(
    threshold: f64,
    max_prob_propagation_distance: u32,
    contig_len: u64,
) -> BandPassActivityProfile {
    BandPassActivityProfile::new(
        "chr1",
        contig_len,
        BandPassActivityProfileParams {
            max_prob_propagation_distance,
            active_prob_threshold: threshold,
            max_filter_size: 0,
            sigma: PositiveSigma::try_new(1.0).unwrap(),
            adaptive_filter_size: false,
        },
    )
}

#[test]
fn empty_profile_pop_ready_regions_is_empty() {
    let mut p = no_smoothing_profile(0.2, 0, 1000);
    assert!(p.pop_ready_regions(0, 1, 50, true).expect("pop").is_empty());
}

#[test]
fn pop_ready_rejects_zero_min_or_max_region_size() {
    let mut p = no_smoothing_profile(0.2, 0, 1000);
    p.add(ActivityProfileState::new("chr1", 1, 0.5))
        .expect("add");
    assert!(p.pop_ready_regions(0, 0, 10, true).is_err());
    assert!(p.pop_ready_regions(0, 1, 0, true).is_err());
}

#[test]
fn non_contiguous_add_returns_error() {
    let mut p = no_smoothing_profile(0.2, 0, 1000);
    p.add(ActivityProfileState::new("chr1", 5, 0.1)).unwrap();
    let err = p
        .add(ActivityProfileState::new("chr1", 7, 0.1))
        .expect_err("gap");
    assert!(err.to_string().contains("not immediately after"), "{err}");
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    #[test]
    fn contiguous_walk_pop_ready_never_panics_and_regions_ordered(
        start in 1u64..500u64,
        probs in prop::collection::vec(0.0..=1.0f64, 5..60)
    ) {
        let contig_len = start + probs.len() as u64 + 200;
        let mut p = no_smoothing_profile(0.15, 2, contig_len);
        for (i, pr) in probs.iter().enumerate() {
            p.add(ActivityProfileState::new("chr1", start + i as u64, *pr)).unwrap();
        }
        let regions = p.pop_ready_regions(0, 1, 100, true).expect("pop");
        let mut last_end: Option<u64> = None;
        for r in regions {
            prop_assert!(r.start <= r.end);
            if let Some(le) = last_end {
                prop_assert!(r.start > le, "regions should advance after drain");
            }
            last_end = Some(r.end);
        }
    }

    #[test]
    fn hq_soft_clip_expansion_respects_contig_bounds(
        pos in 1u64..20u64,
        clip in 0u32..10u32,
        contig_len in 20u64..80u64
    ) {
        let mut p = no_smoothing_profile(0.05, 5, contig_len);
        p.add(ActivityProfileState::high_quality_soft_clips("chr1", pos, 0.9, clip))
            .expect("hq add");
        let _ = p.pop_ready_regions(0, 1, 500, true).expect("pop");
    }
}
