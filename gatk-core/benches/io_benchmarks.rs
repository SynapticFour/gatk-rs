//! I/O benchmarks for GATK-RS FASTA/FASTQ parsers
//! Performance benchmarks to ensure GATK-RS parsers meet or exceed
//! the performance of the original GATK implementation.

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use gatk_core::io::*;
use gatk_core::tests::*;
use std::time::Duration;

fn benchmark_fasta_parsing(c: &mut Criterion) {
    let test_data = TestData::new();

    // Create test files of different sizes
    let sizes = vec![100, 1000, 10000];

    let mut group = c.benchmark_group("fasta_parsing");
    group.measurement_time(Duration::from_secs(10));

    for size in sizes {
        // Create FASTA file
        let mut sequences = Vec::new();
        for i in 0..size {
            sequences.push((
                format!("seq{}", i),
                "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
            ));
        }

        let content: String = sequences
            .iter()
            .map(|(name, seq)| format!(">{}\n{}\n", name, seq))
            .collect();

        let fasta_path = test_data.create_file(format!("test_{}.fasta", size), &content);

        // Benchmark buffered reading
        group.bench_with_input(
            BenchmarkId::new("buffered", size),
            &fasta_path,
            |b, path| {
                b.iter(|| {
                    let mut reader = FastaReader::from_file_buffered(black_box(path)).unwrap();
                    let sequences = reader.read_all_sequences().unwrap();
                    black_box(sequences)
                })
            },
        );

        // Benchmark memory-mapped reading
        group.bench_with_input(
            BenchmarkId::new("memory_mapped", size),
            &fasta_path,
            |b, path| {
                b.iter(|| {
                    let mut reader = FastaReader::from_file_memory_mapped(black_box(path)).unwrap();
                    let sequences = reader.read_all_sequences().unwrap();
                    black_box(sequences)
                })
            },
        );

        // Benchmark iterator
        group.bench_with_input(
            BenchmarkId::new("iterator", size),
            &fasta_path,
            |b, path| {
                b.iter(|| {
                    let mut reader = FastaReader::from_file(black_box(path)).unwrap();
                    let count: usize = reader.iter().count();
                    black_box(count)
                })
            },
        );
    }

    group.finish();
}

fn benchmark_fastq_parsing(c: &mut Criterion) {
    let test_data = TestData::new();

    // Create test files of different sizes
    let sizes = vec![100, 1000, 10000];

    let mut group = c.benchmark_group("fastq_parsing");
    group.measurement_time(Duration::from_secs(10));

    for size in sizes {
        // Create FASTQ file
        let mut content = String::new();
        for i in 0..size {
            content.push_str(&format!("@read{} Test read {}\n", i, i));
            // 68bp read; quality line must match sequence length.
            content
                .push_str("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n");
            content.push_str("+\n");
            content
                .push_str("IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n");
        }

        let fastq_path = test_data.create_file(format!("test_{}.fastq", size), &content);

        // Benchmark buffered reading
        group.bench_with_input(
            BenchmarkId::new("buffered", size),
            &fastq_path,
            |b, path| {
                b.iter(|| {
                    let mut reader = FastqReader::from_file_buffered(black_box(path)).unwrap();
                    let reads = reader.read_all_reads().unwrap();
                    black_box(reads)
                })
            },
        );

        // Benchmark memory-mapped reading
        group.bench_with_input(
            BenchmarkId::new("memory_mapped", size),
            &fastq_path,
            |b, path| {
                b.iter(|| {
                    let mut reader = FastqReader::from_file_memory_mapped(black_box(path)).unwrap();
                    let reads = reader.read_all_reads().unwrap();
                    black_box(reads)
                })
            },
        );

        // Benchmark iterator
        group.bench_with_input(
            BenchmarkId::new("iterator", size),
            &fastq_path,
            |b, path| {
                b.iter(|| {
                    let mut reader = FastqReader::from_file(black_box(path)).unwrap();
                    let count: usize = reader.iter().count();
                    black_box(count)
                })
            },
        );
    }

    group.finish();
}

fn benchmark_sequence_operations(c: &mut Criterion) {
    // Create test sequence
    let sequence = FastaSequence::new(
        "test".to_string(),
        b"ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG".to_vec(),
    );

    let mut group = c.benchmark_group("sequence_operations");

    // Benchmark GC content calculation
    group.bench_function("gc_content", |b| {
        b.iter(|| {
            let gc_content = black_box(&sequence).gc_content();
            black_box(gc_content)
        })
    });

    // Benchmark reverse complement
    group.bench_function("reverse_complement", |b| {
        b.iter(|| {
            let rev_comp = black_box(&sequence).reverse_complement();
            black_box(rev_comp)
        })
    });

    // Benchmark subsequence extraction
    group.bench_function("subsequence", |b| {
        b.iter(|| {
            let subseq = black_box(&sequence).subsequence(10, 20);
            black_box(subseq)
        })
    });

    // Benchmark pattern counting
    group.bench_function("count_pattern", |b| {
        b.iter(|| {
            let count = black_box(&sequence).count_pattern(b"ATCG");
            black_box(count)
        })
    });

    group.finish();
}

fn benchmark_fastq_operations(c: &mut Criterion) {
    // Create test read
    let read = FastqRead::new(
        "test_read".to_string(),
        b"ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG".to_vec(),
        b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII".to_vec(),
    );

    let mut group = c.benchmark_group("fastq_operations");

    // Benchmark quality calculations
    group.bench_function("average_quality", |b| {
        b.iter(|| {
            let avg_quality = black_box(&read).average_quality();
            black_box(avg_quality)
        })
    });

    group.bench_function("min_quality", |b| {
        b.iter(|| {
            let min_quality = black_box(&read).min_quality();
            black_box(min_quality)
        })
    });

    group.bench_function("max_quality", |b| {
        b.iter(|| {
            let max_quality = black_box(&read).max_quality();
            black_box(max_quality)
        })
    });

    group.bench_function("median_quality", |b| {
        b.iter(|| {
            let median_quality = black_box(&read).median_quality();
            black_box(median_quality)
        })
    });

    // Benchmark quality trimming
    group.bench_function("trim_quality", |b| {
        b.iter(|| {
            let trimmed = black_box(&read).trim_quality(35);
            black_box(trimmed)
        })
    });

    // Benchmark reverse complement
    group.bench_function("reverse_complement", |b| {
        b.iter(|| {
            let rev_comp = black_box(&read).reverse_complement();
            black_box(rev_comp)
        })
    });

    // Benchmark validation
    group.bench_function("is_valid", |b| {
        b.iter(|| {
            let is_valid = black_box(&read).is_valid();
            black_box(is_valid)
        })
    });

    group.finish();
}

fn benchmark_fasta_indexing(c: &mut Criterion) {
    let test_data = TestData::new();

    // Create large FASTA file for indexing
    let mut sequences = Vec::new();
    for i in 0..1000 {
        sequences.push((
            format!("chr{}", i),
            "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
        ));
    }

    let content: String = sequences
        .iter()
        .map(|(name, seq)| format!(">{}\n{}\n", name, seq))
        .collect();

    let fasta_path = test_data.create_file("large_index.fasta", &content);

    let mut group = c.benchmark_group("fasta_indexing");
    group.measurement_time(Duration::from_secs(10));

    // Benchmark index building
    group.bench_function("build_index", |b| {
        b.iter(|| {
            let index = FastaIndex::build_from_file(black_box(&fasta_path)).unwrap();
            black_box(index)
        })
    });

    // Build index once for random access tests
    let index = FastaIndex::build_from_file(&fasta_path).unwrap();

    // Benchmark random access
    group.bench_function("random_access", |b| {
        b.iter(|| {
            let seq = black_box(&index).get_sequence("chr500", 10, 20).unwrap();
            black_box(seq)
        })
    });

    // Benchmark save/load index
    let index_path = test_data.path().join("test_index.fai");

    group.bench_function("save_index", |b| {
        b.iter(|| {
            black_box(&index)
                .save_to_file(black_box(&index_path))
                .unwrap();
        })
    });

    group.bench_function("load_index", |b| {
        b.iter(|| {
            let loaded_index = FastaIndex::load_from_file(black_box(&index_path)).unwrap();
            black_box(loaded_index)
        })
    });

    group.finish();
}

fn benchmark_fastq_filtering(c: &mut Criterion) {
    let test_data = TestData::new();

    // Create mixed quality FASTQ file
    let mut content = String::new();
    for i in 0..1000 {
        let quality = if i % 3 == 0 {
            "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
        // High
        } else if i % 3 == 1 {
            "HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH"
        // Medium
        } else {
            "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            // Low
        };

        content.push_str(&format!("@read{}\n", i));
        content.push_str("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n");
        content.push_str("+\n");
        content.push_str(quality);
        content.push('\n');
    }

    let fastq_path = test_data.create_file("mixed_quality.fastq", &content);

    let mut group = c.benchmark_group("fastq_filtering");
    group.measurement_time(Duration::from_secs(10));

    // Benchmark quality filtering
    group.bench_function("filter_by_quality", |b| {
        b.iter(|| {
            let mut reader = FastqReader::from_file(black_box(&fastq_path)).unwrap();
            let filtered = reader.filter_by_quality(35).unwrap();
            black_box(filtered)
        })
    });

    // Benchmark length filtering
    group.bench_function("filter_by_length", |b| {
        b.iter(|| {
            let mut reader = FastqReader::from_file(black_box(&fastq_path)).unwrap();
            let filtered = reader.filter_by_length(60, 60).unwrap();
            black_box(filtered)
        })
    });

    // Benchmark sampling
    group.bench_function("sample_reads", |b| {
        b.iter(|| {
            let mut reader = FastqReader::from_file(black_box(&fastq_path)).unwrap();
            let sampled = reader.sample_reads(100).unwrap();
            black_box(sampled)
        })
    });

    group.finish();
}

fn benchmark_writing(c: &mut Criterion) {
    let test_data = TestData::new();

    let mut group = c.benchmark_group("io_writing");
    group.measurement_time(Duration::from_secs(10));

    // Prepare test data
    let fasta_sequences: Vec<FastaSequence> = (0..1000)
        .map(|i| {
            FastaSequence::new(
                format!("seq{}", i),
                b"ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG".to_vec(),
            )
        })
        .collect();

    let fastq_reads: Vec<FastqRead> = (0..1000)
        .map(|i| {
            FastqRead::new(
                format!("read{}", i),
                b"ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG".to_vec(),
                b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII".to_vec(),
            )
        })
        .collect();

    // Benchmark FASTA writing
    group.bench_function("write_fasta", |b| {
        let fasta_path = test_data.path().join("bench_write.fasta");
        b.iter(|| {
            let mut writer = FastaWriter::new(black_box(&fasta_path)).unwrap();
            writer.write_sequences(black_box(&fasta_sequences)).unwrap();
            writer.finish().unwrap();
        })
    });

    // Benchmark FASTQ writing
    group.bench_function("write_fastq", |b| {
        let fastq_path = test_data.path().join("bench_write.fastq");
        b.iter(|| {
            let mut writer = FastqWriter::new(black_box(&fastq_path)).unwrap();
            writer.write_reads(black_box(&fastq_reads)).unwrap();
            writer.finish().unwrap();
        })
    });

    group.finish();
}

fn benchmark_memory_usage(c: &mut Criterion) {
    let test_data = TestData::new();

    // Create large FASTA file
    let mut sequences = Vec::new();
    for i in 0..10000 {
        sequences.push((
            format!("seq{}", i),
            "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
        ));
    }

    let content: String = sequences
        .iter()
        .map(|(name, seq)| format!(">{}\n{}\n", name, seq))
        .collect();

    let fasta_path = test_data.create_file("memory_test.fasta", &content);

    let mut group = c.benchmark_group("memory_usage");

    // Benchmark memory usage during parsing
    group.bench_function("parse_memory_usage", |b| {
        b.iter(|| {
            let mut reader = FastaReader::from_file(black_box(&fasta_path)).unwrap();
            let sequences = reader.read_all_sequences().unwrap();

            // Calculate total memory used
            let total_memory: usize = sequences.iter().map(|s| s.sequence.len()).sum();
            black_box(total_memory)
        })
    });

    group.finish();
}

criterion_group!(
    benches,
    benchmark_fasta_parsing,
    benchmark_fastq_parsing,
    benchmark_sequence_operations,
    benchmark_fastq_operations,
    benchmark_fasta_indexing,
    benchmark_fastq_filtering,
    benchmark_writing,
    benchmark_memory_usage
);

criterion_main!(benches);
