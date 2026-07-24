//! Concurrent reproducibility for pure activity-profile math (Phase 4, step 61).
use gatk_haplotypecaller::{
    haplotype_caller_activity_profile_state_single_sample, BandPassActivityProfile,
    BandPassActivityProfileParams, HaplotypeCallerActivityScoringParams, PileupObservation,
    PositiveSigma,
};
use std::collections::BTreeMap;
use std::sync::mpsc;
use std::thread;

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

fn hotspot_reads() -> Vec<SyntheticRead> {
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
    reads
}

#[test]
fn concurrent_identical_runs_match_bit_for_bit() {
    let reads = hotspot_reads();
    let baseline = active_intervals_from_reads_with_cap(&reads, 3, 0.2);

    let (tx, rx) = mpsc::channel();
    let mut handles = Vec::new();
    for _ in 0..16 {
        let tx = tx.clone();
        let r = reads.clone();
        handles.push(thread::spawn(move || {
            let got = active_intervals_from_reads_with_cap(&r, 3, 0.2);
            tx.send(got).expect("send");
        }));
    }
    drop(tx);
    for h in handles {
        h.join().expect("thread");
    }
    for got in rx {
        assert_eq!(got, baseline);
    }
}
