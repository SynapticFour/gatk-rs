//! Allocation probe for the activity-scoring hot path.
//! Compares the pre-audit pattern (clone nonempty strata every locus) against the
//! post-audit borrow path on the same synthetic multi-sample pileups.
//! ```text
//! cargo run -p gatk-haplotypecaller --release --example allocation_rss_probe
//! /usr/bin/time -l cargo run -p gatk-haplotypecaller --release --example allocation_rss_probe
//! ```

use gatk_haplotypecaller::{
    haplotype_caller_activity_profile_state_multi_sample, HaplotypeCallerActivityScoringParams,
    PileupObservation,
};
use std::time::Instant;

fn obs(alt: bool) -> PileupObservation {
    PileupObservation {
        read_base: if alt { b'T' } else { b'A' },
        qual: 30,
        is_deletion: false,
        is_alt: alt,
        is_next_to_soft_clip: false,
        read_hq_soft_clip_base_count: 0,
    }
}

fn main() {
    let mode = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "both".to_string());
    let params = HaplotypeCallerActivityScoringParams::default();
    // ~3 samples × 80 obs — sized to make per-locus clones visible in time / RSS.
    let sample_pileups: Vec<Vec<PileupObservation>> = (0..3)
        .map(|s| (0..80).map(|i| obs((i + s) % 3 == 0)).collect())
        .collect();
    const LOCI: usize = 50_000;

    let bytes_per_locus_clone = sample_pileups.iter().map(|p| p.len()).sum::<usize>()
        * std::mem::size_of::<PileupObservation>();
    println!("allocation_rss_probe mode={mode}");
    println!("loci={LOCI}");
    println!("approx_clone_bytes_per_locus={bytes_per_locus_clone}");
    println!(
        "approx_clone_bytes_total={}",
        bytes_per_locus_clone.saturating_mul(LOCI)
    );

    let mut sink = 0.0_f64;
    if mode == "clone" || mode == "both" {
        let t0 = Instant::now();
        for locus in 0..LOCI {
            // Pre-audit pattern: clone nonempty strata every locus.
            let nonempty: Vec<Vec<PileupObservation>> = sample_pileups
                .iter()
                .filter(|p| !p.is_empty())
                .cloned()
                .collect();
            let st = haplotype_caller_activity_profile_state_multi_sample(
                "chr21",
                (locus as u64) + 1,
                &nonempty,
                &params,
            );
            sink += st.active_prob;
        }
        println!("clone_path_ms={:.3}", t0.elapsed().as_secs_f64() * 1000.0);
    }
    if mode == "borrow" || mode == "both" {
        let t1 = Instant::now();
        for locus in 0..LOCI {
            let st = haplotype_caller_activity_profile_state_multi_sample(
                "chr21",
                (locus as u64) + 1,
                &sample_pileups,
                &params,
            );
            sink += st.active_prob;
        }
        println!("borrow_path_ms={:.3}", t1.elapsed().as_secs_f64() * 1000.0);
    }
    println!("sink={sink}");
}
