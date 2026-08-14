//! Benchmark tests for GATK-RS core components
#![allow(clippy::result_large_err)]

use criterion::{criterion_group, criterion_main, Criterion};
use gatk_core::io::FastaSequence;
use gatk_core::memory::MemoryMappedFile;
use gatk_core::tests::TestData;
use gatk_core::types::{Allele, Base, BaseQuality, GenomicPosition, VariantContext};
use gatk_core::utils;
use std::hint::black_box;

fn benchmark_fasta_sequence_ops(c: &mut Criterion) {
    let sequence = FastaSequence::new("seq".to_string(), vec![b'A'; 1_000_000]);
    c.bench_function("fasta_gc_content", |b| {
        b.iter(|| black_box(sequence.gc_content()))
    });
    c.bench_function("fasta_revcomp", |b| {
        b.iter(|| black_box(sequence.reverse_complement()))
    });
    c.bench_function("fasta_subsequence", |b| {
        b.iter(|| black_box(sequence.subsequence(100_000, 10_000)))
    });
}

fn benchmark_variant_ops(c: &mut Criterion) {
    let variants: Vec<VariantContext> = (0..1000)
        .map(|i| {
            VariantContext::new(
                GenomicPosition {
                    contig: 1,
                    position: (i * 100) as u64,
                },
                Allele::new(vec![Base::A]),
                vec![Allele::new(vec![if i % 2 == 0 {
                    Base::T
                } else {
                    Base::G
                }])],
            )
        })
        .collect();
    c.bench_function("variant_type", |b| {
        b.iter(|| {
            for v in &variants {
                black_box(v.variant_type());
            }
        })
    });

    let qualities: Vec<BaseQuality> = (0..=93).map(BaseQuality::new).collect();
    c.bench_function("quality_error_prob", |b| {
        b.iter(|| {
            for q in &qualities {
                black_box(q.error_probability());
                black_box(q.phred_score());
            }
        })
    });
}

fn benchmark_mmap(c: &mut Criterion) {
    let data = TestData::new();
    let file_path = data.create_file("benchmark.txt", &"A".repeat(1_000_000));
    c.bench_function("memory_mapped_access", |b| {
        b.iter(|| {
            let mmap_file = MemoryMappedFile::open(file_path.to_str().unwrap()).unwrap();
            black_box(mmap_file.as_bytes());
            black_box(mmap_file.slice(100_000, 1000).unwrap());
        })
    });
}

fn benchmark_hamming_distance(c: &mut Criterion) {
    let seq1 = vec![Base::A, Base::T, Base::C, Base::G];
    let seq2 = vec![Base::A, Base::G, Base::C, Base::T];
    c.bench_function("hamming_distance", |b| {
        b.iter(|| black_box(utils::hamming_distance(&seq1, &seq2)))
    });
}

criterion_group!(
    benches,
    benchmark_fasta_sequence_ops,
    benchmark_variant_ops,
    benchmark_mmap,
    benchmark_hamming_distance
);

criterion_main!(benches);
