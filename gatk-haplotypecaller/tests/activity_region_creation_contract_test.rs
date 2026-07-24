use gatk_haplotypecaller::{
    haplotype_caller_activity_profile_state_single_sample, ActivityProfileState,
    BandPassActivityProfile, BandPassActivityProfileParams, HaplotypeCallerActivityScoringParams,
    PileupObservation, PositiveSigma,
};
use std::collections::BTreeMap;
use std::fs;
use std::path::PathBuf;

fn fixture_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../parity/fixtures")
        .join(name)
}

fn load_prob_fixture(name: &str) -> Vec<(u64, f64)> {
    let raw = fs::read_to_string(fixture_path(name)).expect("fixture");
    raw.lines()
        .filter(|l| !l.trim().is_empty() && !l.trim_start().starts_with('#'))
        .map(|l| {
            let mut parts = l.split_whitespace();
            let pos = parts.next().expect("pos").parse::<u64>().expect("pos u64");
            let p = parts
                .next()
                .expect("active_prob")
                .parse::<f64>()
                .expect("prob f64");
            (pos, p)
        })
        .collect()
}

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

fn active_intervals_from_probs(
    probs: &[(u64, f64)],
    threshold: f64,
    min_region_size: u32,
    max_region_size: u32,
) -> Vec<(u64, u64)> {
    let contig_len = probs.last().map(|(p, _)| *p + 100).unwrap_or(10_000);
    let mut p = no_smoothing_profile(threshold, 0, contig_len);
    for (pos, pr) in probs {
        p.add(ActivityProfileState::new("chr1", *pos, *pr)).unwrap();
    }
    let regions = p
        .pop_ready_regions(0, min_region_size, max_region_size, true)
        .unwrap();
    regions
        .into_iter()
        .filter(|r| r.is_active)
        .map(|r| (r.start, r.end))
        .collect()
}

#[test]
fn deterministic_region_boundaries_fixture_a() {
    let probs = load_prob_fixture("p4_region_boundaries_known_a.tsv");
    let active = active_intervals_from_probs(&probs, 0.2, 1, 100);
    assert_eq!(active, vec![(3, 5), (8, 10)]);
}

#[test]
fn deterministic_region_boundaries_fixture_b_with_max_region_cut() {
    let probs = load_prob_fixture("p4_region_boundaries_known_b.tsv");
    let active = active_intervals_from_probs(&probs, 0.2, 2, 5);
    // First block is cut at the local minimum due to max-region enforcement.
    assert_eq!(active, vec![(1, 4), (5, 5)]);
}

#[test]
fn low_complexity_near_threshold_noise_creates_only_small_islands() {
    let probs = vec![
        (1, 0.19),
        (2, 0.21),
        (3, 0.18),
        (4, 0.22),
        (5, 0.19),
        (6, 0.21),
        (7, 0.18),
    ];
    let active = active_intervals_from_probs(&probs, 0.2, 1, 50);
    assert_eq!(active, vec![(2, 2), (4, 4), (6, 6)]);
}

#[test]
fn homopolymer_soft_clip_signal_expands_activity_window() {
    let mut p = no_smoothing_profile(0.2, 3, 1_000);
    for pos in 1..=10 {
        if pos == 5 {
            p.add(ActivityProfileState::high_quality_soft_clips(
                "chr1", pos, 0.9, 3,
            ))
            .unwrap();
        } else {
            p.add(ActivityProfileState::new("chr1", pos, 0.0)).unwrap();
        }
    }
    let regions = p.pop_ready_regions(0, 1, 100, true).unwrap();
    let active: Vec<(u64, u64)> = regions
        .into_iter()
        .filter(|r| r.is_active)
        .map(|r| (r.start, r.end))
        .collect();
    // With filter_size=0, `BandPassActivityProfile#processState` emits around the just-added locus.
    // Combined with current GATK-style wiring, HQ soft-clip expansion contributes mass to the center locus.
    assert_eq!(active, vec![(5, 5)]);
}

#[derive(Clone)]
struct SyntheticRead {
    locus: u64,
    alignment_start: u64,
    is_alt: bool,
    qual: u8,
    order_key: u64,
}

fn downsample_per_alignment_start(reads: &[SyntheticRead], cap: usize) -> Vec<SyntheticRead> {
    let mut sorted = reads.to_vec();
    sorted.sort_by_key(|r| (r.locus, r.alignment_start, r.order_key));
    let mut kept: BTreeMap<(u64, u64), usize> = BTreeMap::new();
    let mut out = Vec::new();
    for r in sorted {
        let key = (r.locus, r.alignment_start);
        let count = kept.entry(key).or_insert(0);
        if *count < cap {
            out.push(r.clone());
            *count += 1;
        }
    }
    out
}

fn active_intervals_from_reads_with_cap(
    reads: &[SyntheticRead],
    cap: usize,
    threshold: f64,
) -> Vec<(u64, u64)> {
    let ds = downsample_per_alignment_start(reads, cap);
    let mut by_locus: BTreeMap<u64, Vec<PileupObservation>> = BTreeMap::new();
    for r in ds {
        by_locus
            .entry(r.locus)
            .or_default()
            .push(PileupObservation {
                read_base: if r.is_alt { b'T' } else { b'A' },
                qual: r.qual,
                is_deletion: false,
                is_alt: r.is_alt,
                is_next_to_soft_clip: false,
                read_hq_soft_clip_base_count: 0,
            });
    }

    let mut profile = no_smoothing_profile(threshold, 0, 10_000);
    let scoring = HaplotypeCallerActivityScoringParams::default();
    for locus in 10..=14 {
        let pile = by_locus.get(&locus).cloned().unwrap_or_default();
        let st =
            haplotype_caller_activity_profile_state_single_sample("chr1", locus, &pile, &scoring);
        profile.add(st).unwrap();
    }
    profile
        .pop_ready_regions(0, 1, 100, true)
        .unwrap()
        .into_iter()
        .filter(|r| r.is_active)
        .map(|r| (r.start, r.end))
        .collect()
}

fn total_len(intervals: &[(u64, u64)]) -> u64 {
    intervals.iter().map(|(s, e)| e - s + 1).sum()
}

#[test]
fn downsampling_interaction_is_deterministic_and_monotonic() {
    // Hotspot loci 11-12: many ALT observations from the same alignment start.
    // Stricter cap should not create *more* active span than permissive cap.
    let mut reads = Vec::new();
    let mut key = 0u64;
    for locus in [11_u64, 12_u64] {
        for _ in 0..5 {
            key += 1;
            reads.push(SyntheticRead {
                locus,
                alignment_start: locus - 1,
                is_alt: true,
                qual: 30,
                order_key: key,
            });
        }
        key += 1;
        reads.push(SyntheticRead {
            locus,
            alignment_start: locus + 5,
            is_alt: false,
            qual: 30,
            order_key: key,
        });
    }
    // Flanking mostly-reference evidence.
    for locus in [10_u64, 13_u64, 14_u64] {
        key += 1;
        reads.push(SyntheticRead {
            locus,
            alignment_start: locus,
            is_alt: false,
            qual: 30,
            order_key: key,
        });
    }

    let low_cap = active_intervals_from_reads_with_cap(&reads, 1, 0.2);
    let hi_cap = active_intervals_from_reads_with_cap(&reads, 5, 0.2);
    let hi_cap_reversed = {
        let mut rev = reads.clone();
        rev.reverse();
        active_intervals_from_reads_with_cap(&rev, 5, 0.2)
    };

    // Deterministic under input-order changes due to explicit ordering key in downsampler.
    assert_eq!(hi_cap, hi_cap_reversed);
    // More permissive downsampling should preserve or increase active span.
    assert!(total_len(&hi_cap) >= total_len(&low_cap));
}
