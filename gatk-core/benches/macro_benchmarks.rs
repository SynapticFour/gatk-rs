use criterion::{criterion_group, criterion_main, Bencher};
use gatk_core::bench_test;

bench_test!(macro_sum_bench, |b: &mut Bencher| {
    b.iter(|| {
        let values: Vec<u64> = (0..10_000).collect();
        let _sum: u64 = values.iter().sum();
    });
});

criterion_group!(macro_benches, macro_sum_bench);
criterion_main!(macro_benches);
