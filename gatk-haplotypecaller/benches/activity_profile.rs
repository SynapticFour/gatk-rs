//! Microbenchmarks for Phase-4 activity profile hot paths (step 59).
#![allow(clippy::result_large_err)]

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use gatk_haplotypecaller::{
    ActivityProfileState, BandPassActivityProfile, BandPassActivityProfileParams, PositiveSigma,
};

fn bench_gaussian_kernel_gatk_defaults(c: &mut Criterion) {
    let params = BandPassActivityProfileParams::gatk_haplotype_caller_defaults();
    c.bench_function("band_pass_normalized_kernel_gatk_defaults", |b| {
        b.iter(|| black_box(params.normalized_kernel()))
    });
}

fn bench_add_and_pop_ready_regions(c: &mut Criterion) {
    c.bench_function("band_pass_add_then_pop_ready_256_loci", |b| {
        b.iter(|| {
            let contig: std::sync::Arc<str> = std::sync::Arc::from("chr1");
            let mut p = BandPassActivityProfile::new(
                contig.clone(),
                10_000,
                BandPassActivityProfileParams {
                    max_prob_propagation_distance: 12,
                    active_prob_threshold: 0.002,
                    max_filter_size: 8,
                    sigma: PositiveSigma::try_new(2.0).unwrap(),
                    adaptive_filter_size: false,
                },
            );
            for pos in 1u64..=256 {
                let pr = ((pos.wrapping_mul(17)) % 100) as f64 * 0.001;
                p.add(ActivityProfileState::new(contig.clone(), pos, pr))
                    .unwrap();
            }
            black_box(p.pop_ready_regions(50, 50, 300, true).unwrap())
        });
    });
}

criterion_group!(
    benches,
    bench_gaussian_kernel_gatk_defaults,
    bench_add_and_pop_ready_regions
);
criterion_main!(benches);
