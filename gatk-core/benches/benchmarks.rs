//! Benchmark tests for GATK-RS core components
#![allow(clippy::result_large_err)]

use criterion::{criterion_group, criterion_main, Criterion};
use gatk_core::io::FastaSequence;
use gatk_core::memory::{
    GenomicCache, GenomicInterval, IntervalTree, MemoryMappedFile, MemoryPool, StreamProcessor,
};
use gatk_core::tests::{mocks, TestData};
use gatk_core::types::{Allele, Base, BaseQuality, GenomicPosition, VariantContext};
use gatk_core::{math, utils};
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

fn benchmark_memory_ops(c: &mut Criterion) {
    let pool = MemoryPool::new(100);
    c.bench_function("memory_pool_allocation", |b| {
        b.iter(|| {
            let buffer = pool.get_buffer(1024);
            black_box(buffer.len());
            pool.return_buffer(buffer);
        })
    });

    let cache: GenomicCache<String, String> = GenomicCache::new(1000);
    c.bench_function("cache_put_get", |b| {
        b.iter(|| {
            for i in 0..1000 {
                cache.put(format!("key_{i}"), format!("value_{i}"));
            }
            black_box(cache.get(&"key_500".to_string()));
        })
    });
}

fn benchmark_interval_tree(c: &mut Criterion) {
    let mut tree: IntervalTree<i32> = IntervalTree::new();
    for i in 0..10_000 {
        tree.insert(GenomicInterval {
            chromosome: format!("chr{}", i % 10),
            start: i * 10,
            end: i * 10 + 5,
            data: i as i32,
        });
    }
    tree.sort();
    c.bench_function("interval_query", |b| {
        b.iter(|| {
            for i in 0..1000 {
                black_box(tree.find_overlapping("chr1", i * 10 + 2));
            }
        })
    });
}

fn benchmark_variant_and_math(c: &mut Criterion) {
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

    c.bench_function("log_addition", |b| {
        b.iter(|| black_box(math::likelihood::log_add(-10.0, -12.0)))
    });
}

fn benchmark_stream_and_mmap(c: &mut Criterion) {
    let pool = std::sync::Arc::new(MemoryPool::new(100));
    let processor = StreamProcessor::new(1024, pool);
    c.bench_function("stream_processing", |b| {
        b.iter(|| {
            let mock_reader = mocks::MockFileReader::new(&"A".repeat(10_000));
            let _ = processor.process_chunks(mock_reader, |chunk: &[u8]| {
                black_box(chunk.len());
                Ok(())
            });
        })
    });

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
    benchmark_memory_ops,
    benchmark_interval_tree,
    benchmark_variant_and_math,
    benchmark_stream_and_mmap,
    benchmark_hamming_distance
);

criterion_main!(benches);
