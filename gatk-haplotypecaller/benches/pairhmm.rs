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

/// Equal-length haps with long shared prefixes (GIAB assembly-like).
fn make_realistic_prefix_batch(
    hap_count: usize,
    read_len: usize,
) -> (Vec<u8>, Vec<u8>, Vec<Vec<u8>>, Vec<u8>, Vec<u8>, Vec<u8>) {
    let read = b"ACGT".repeat(read_len / 4);
    let quals = vec![30u8; read.len()];
    let ins = vec![GATK_PARITY_DEFAULT_INS_QUAL; read.len()];
    let del = vec![GATK_PARITY_DEFAULT_DEL_QUAL; read.len()];
    let gcp = vec![GATK_PARITY_DEFAULT_GCP; read.len()];
    let hap_len = read_len + 20;
    let mut base = b"ACGT".repeat(hap_len / 4);
    base.resize(hap_len, b'N');
    let haps = (0..hap_count)
        .map(|k| {
            let mut h = base.clone();
            if k > 0 {
                let pos = hap_len.saturating_sub(1 + (k % 8));
                h[pos] = b'T';
            }
            h
        })
        .collect::<Vec<_>>();
    (read, quals, haps, ins, del, gcp)
}

fn make_read_tile(
    read_count: usize,
    hap_count: usize,
    read_len: usize,
) -> (
    Vec<(Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>)>,
    Vec<Vec<u8>>,
) {
    let haps = (0..hap_count)
        .map(|k| {
            let mut s = b"ACGT".repeat(read_len / 4);
            if k % 3 != 0 {
                let pos = k % s.len();
                s[pos] = b'T';
            }
            s
        })
        .collect::<Vec<_>>();
    let reads = (0..read_count)
        .map(|r| {
            let mut read = b"ACGT".repeat(read_len / 4);
            if r % 5 != 0 {
                let pos = r % read.len();
                read[pos] = b'G';
            }
            let quals = vec![30u8; read.len()];
            let ins = vec![GATK_PARITY_DEFAULT_INS_QUAL; read.len()];
            let del = vec![GATK_PARITY_DEFAULT_DEL_QUAL; read.len()];
            let gcp = vec![GATK_PARITY_DEFAULT_GCP; read.len()];
            (read, quals, ins, del, gcp)
        })
        .collect();
    (reads, haps)
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

    // Phenotype classes (read length × hap count), not locus pins.
    for read_len in [100usize, 200, 300] {
        for hap_count in [8usize, 32, 64, 128] {
            for (name, backend) in backends {
                group.bench_with_input(
                    BenchmarkId::new(format!("{name}_r{read_len}_h"), hap_count),
                    &hap_count,
                    |b, &n| {
                        b.iter_batched(
                            || make_batch(n, read_len),
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
    }
    group.finish();
}

/// Architecture C wavefront vs scalar / hap-axis SIMD.
fn bench_pairhmm_wavefront(c: &mut Criterion) {
    let mut group = c.benchmark_group("pairhmm_wavefront");
    let backends = [
        (
            "logless_scalar",
            resolve_pair_hmm_impl(PairHmmImpl::LoglessPairHmm),
        ),
        ("simd", resolve_pair_hmm_impl(PairHmmImpl::Simd)),
        ("wavefront", resolve_pair_hmm_impl(PairHmmImpl::Wavefront)),
    ];

    // 1×N matrix (read 200).
    let read_len = 200usize;
    for hap_count in [1usize, 8, 16, 32, 64] {
        for (name, backend) in backends {
            group.bench_with_input(
                BenchmarkId::new(format!("{name}_1x{hap_count}"), hap_count),
                &hap_count,
                |b, &n| {
                    b.iter_batched(
                        || make_batch(n, read_len),
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

    // Multi-read tiles (prep amortized per read; kernels sequential in tile).
    for (reads_n, haps_n) in [(8usize, 32usize), (16usize, 32usize)] {
        for (name, backend) in backends {
            group.bench_with_input(
                BenchmarkId::new(format!("{name}_tile_{reads_n}x{haps_n}"), reads_n),
                &reads_n,
                |b, &_rn| {
                    b.iter_batched(
                        || make_read_tile(reads_n, haps_n, read_len),
                        |(reads, haps)| {
                            let hrefs: Vec<&[u8]> = haps.iter().map(|h| h.as_slice()).collect();
                            for (read, quals, ins, del, gcp) in &reads {
                                let out = score_read_haps_logless(
                                    backend, read, quals, &hrefs, ins, del, gcp,
                                )
                                .unwrap();
                                assert_eq!(out.len(), haps_n);
                            }
                        },
                        BatchSize::SmallInput,
                    )
                },
            );
        }
    }

    // Realistic shared-prefix phenotype.
    for hap_count in [8usize, 32, 64] {
        for (name, backend) in backends {
            group.bench_with_input(
                BenchmarkId::new(format!("{name}_realistic_h"), hap_count),
                &hap_count,
                |b, &n| {
                    b.iter_batched(
                        || make_realistic_prefix_batch(n, read_len),
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
    bench_pairhmm_logless_simd,
    bench_pairhmm_wavefront
);
criterion_main!(benches);
