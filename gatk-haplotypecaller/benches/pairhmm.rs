use criterion::{criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};
use gatk_haplotypecaller::{
    pairhmm_log10_likelihoods_vectorized, resolve_pair_hmm_impl, score_read_haps_logless,
    PairHmmImpl, PairHmmParams, GATK_PARITY_DEFAULT_DEL_QUAL, GATK_PARITY_DEFAULT_GCP,
    GATK_PARITY_DEFAULT_INS_QUAL,
};

fn make_batch(
    hap_count: usize,
    read_len: usize,
) -> (Vec<u8>, Vec<u8>, Vec<Vec<u8>>, Vec<u8>, Vec<u8>, Vec<u8>) {
    let read = b"ACGT".repeat(read_len / 4);
    let quals = vec![30u8; read.len()];
    let ins = vec![GATK_PARITY_DEFAULT_INS_QUAL; read.len()];
    let del = vec![GATK_PARITY_DEFAULT_DEL_QUAL; read.len()];
    let gcp = vec![GATK_PARITY_DEFAULT_GCP; read.len()];
    let haps = (0..hap_count)
        .map(|k| {
            let mut s = read.clone();
            if k % 3 != 0 {
                let pos = k % s.len();
                s[pos] = b'T';
            }
            s
        })
        .collect::<Vec<_>>();
    (read, quals, haps, ins, del, gcp)
}

fn bench_pairhmm_scaffold_vectorized(c: &mut Criterion) {
    let mut group = c.benchmark_group("pairhmm_scaffold_vectorized");
    let params = PairHmmParams::default();

    for hap_count in [4usize, 16usize, 64usize] {
        group.bench_with_input(BenchmarkId::new("batch", hap_count), &hap_count, |b, &n| {
            b.iter_batched(
                || {
                    let read = "ACGT".repeat(50);
                    let quals = vec![30u8; read.len()];
                    let haps = (0..n)
                        .map(|k| {
                            if k % 3 == 0 {
                                read.clone()
                            } else {
                                let mut s = read.clone().into_bytes();
                                let pos = k % s.len();
                                s[pos] = b'T';
                                String::from_utf8(s).expect("utf8")
                            }
                        })
                        .collect::<Vec<_>>();
                    (read, quals, haps)
                },
                |(read, quals, haps)| {
                    let out =
                        pairhmm_log10_likelihoods_vectorized(&read, &quals, 60, &haps, &params)
                            .unwrap();
                    assert_eq!(out.len(), haps.len());
                },
                BatchSize::SmallInput,
            )
        });
    }
    group.finish();
}

fn bench_pairhmm_logless_simd(c: &mut Criterion) {
    let mut group = c.benchmark_group("pairhmm_logless_simd");
    let backends = [
        (
            "logless_scalar",
            resolve_pair_hmm_impl(PairHmmImpl::LoglessPairHmm),
        ),
        ("simd", resolve_pair_hmm_impl(PairHmmImpl::Simd)),
        ("simd_f32", resolve_pair_hmm_impl(PairHmmImpl::SimdF32)),
        (
            "log10_scalar",
            resolve_pair_hmm_impl(PairHmmImpl::Log10PairHmm),
        ),
    ];

    for hap_count in [8usize, 32usize, 64usize] {
        for (name, backend) in backends {
            group.bench_with_input(
                BenchmarkId::new(format!("{name}_haps"), hap_count),
                &hap_count,
                |b, &n| {
                    b.iter_batched(
                        || make_batch(n, 200),
                        |(read, quals, haps, ins, del, gcp)| {
                            let refs: Vec<&[u8]> = haps.iter().map(|h| h.as_slice()).collect();
                            let out = score_read_haps_logless(
                                backend, &read, &quals, &refs, &ins, &del, &gcp,
                            )
                            .unwrap();
                            assert_eq!(out.len(), n);
                        },
                        BatchSize::SmallInput,
                    )
                },
            );
        }
    }
    group.finish();
}

criterion_group!(
    benches,
    bench_pairhmm_scaffold_vectorized,
    bench_pairhmm_logless_simd
);
criterion_main!(benches);
