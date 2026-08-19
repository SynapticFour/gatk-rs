//! Microbenchmarks for assembly-graph / read-threading hot paths.
#![allow(clippy::result_large_err)]

use std::collections::HashMap;
use std::hint::black_box;
use std::sync::Arc;
use std::time::Duration;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use gatk_haplotypecaller::kmer_key::{KmerKey, RollingKmer};
use gatk_haplotypecaller::{
    assembly_graph_from_reads_threading, AssemblyGraph, AssemblyGraphParams,
    AssemblyGraphPruningParams, AssemblyRead, KmerSize,
};

fn smoke_criterion() -> Criterion {
    Criterion::default()
        .warm_up_time(Duration::from_millis(200))
        .measurement_time(Duration::from_secs(2))
        .sample_size(20)
}

fn mk_read(seq: &str, q: u8) -> AssemblyRead {
    AssemblyRead {
        bases: seq.as_bytes().to_vec(),
        base_quals: vec![q; seq.len()],
    }
}

/// ~120 bp ACGT backbone so both k=10 and k=25 form many overlapping windows.
fn backbone_120() -> String {
    "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        .to_string()
}

/// Shallow coverage (sparse region).
fn low_depth_reads() -> Vec<AssemblyRead> {
    let mut reads = Vec::with_capacity(16);
    let base = backbone_120();
    for i in 0..16 {
        let mut seq = base.clone();
        let cut = 40 + (i % 8);
        if i % 4 == 0 {
            seq.insert_str(cut, "TTAA");
        } else {
            seq.insert_str(cut, "TGAA");
        }
        reads.push(mk_read(&seq, 30));
    }
    reads
}

fn medium_depth_reads() -> Vec<AssemblyRead> {
    let mut reads = Vec::with_capacity(64);
    let base = backbone_120();
    for i in 0..64 {
        let mut seq = base.clone();
        let cut = 40 + (i % 16);
        match i % 3 {
            0 => seq.insert_str(cut, "TTAA"),
            1 => seq.insert_str(cut, "TCAA"),
            _ => seq.insert_str(cut, "TGAA"),
        }
        reads.push(mk_read(&seq, 30));
    }
    reads
}

fn high_depth_reads() -> Vec<AssemblyRead> {
    let mut reads = Vec::with_capacity(512);
    let base = backbone_120();
    for i in 0..512 {
        let mut seq = base.clone();
        let cut = 36 + (i % 24);
        match i % 7 {
            0 => seq.insert_str(cut, "TTAA"),
            1 => seq.insert_str(cut, "TCAA"),
            _ => seq.insert_str(cut, "TGAA"),
        }
        // Local homopolymer to stress non-unique / pruning edges without emptying k=25.
        seq.push_str("AAAAAAAA");
        reads.push(mk_read(&seq, 30));
    }
    reads
}

/// A/B: sliding-window map fill — Arc<[u8]> keys vs packed [`KmerKey`].
fn fill_arc_map(seq: &[u8], k: usize) -> usize {
    let mut map: HashMap<Arc<[u8]>, u32> = HashMap::new();
    if seq.len() < k {
        return 0;
    }
    for start in 0..=seq.len() - k {
        let key: Arc<[u8]> = Arc::from(&seq[start..start + k]);
        *map.entry(key).or_insert(0) += 1;
    }
    map.len()
}

fn fill_packed_map(seq: &[u8], k: usize) -> usize {
    let mut map: HashMap<KmerKey, u32> = HashMap::new();
    if seq.len() < k {
        return 0;
    }
    let mut roll = RollingKmer::new(k);
    for start in 0..=seq.len() - k {
        let key = roll.key_at(seq, start);
        *map.entry(key).or_insert(0) += 1;
    }
    map.len()
}

fn threading_params(k: u16) -> AssemblyGraphParams {
    AssemblyGraphParams {
        kmer_size: KmerSize::try_new(k).expect("k"),
        min_base_quality: 10,
        min_edge_weight: 2,
        dangling_path_max_nodes: 4,
        max_haplotypes: 32,
        max_haplotype_bases: 128,
        ..Default::default()
    }
}

fn assemble_pipeline(reads: &[AssemblyRead], k: u16) {
    let params = threading_params(k);
    let mut graph = AssemblyGraph::from_reads(reads, &params).expect("graph");
    let mut pruning = AssemblyGraphPruningParams::gatk_haplotype_caller_defaults();
    pruning.min_prune_factor = params.min_edge_weight;
    graph.apply_pruning(&pruning);
    graph.remove_dangling_paths(params.dangling_path_max_nodes);
    graph.cleanup_isolated_nodes();
    black_box(
        graph.extract_candidate_haplotypes(params.max_haplotypes, params.max_haplotype_bases),
    );
}

/// Read-threading build only (packed k-mer path).
fn threading_build_only(reads: &[AssemblyRead], k: u16) {
    let params = threading_params(k);
    let g = assembly_graph_from_reads_threading(reads, &params).expect("rt");
    black_box(g.edges_sorted().len());
}

fn bench_depth_matrix(c: &mut Criterion) {
    let mut group = c.benchmark_group("assembly_graph_depth");
    let suites = [
        ("low", low_depth_reads()),
        ("medium", medium_depth_reads()),
        ("high", high_depth_reads()),
    ];
    for (label, reads) in &suites {
        for k in [10u16, 25] {
            group.bench_with_input(
                BenchmarkId::new(format!("full_pipeline_{label}_k{k}"), reads.len()),
                reads,
                |b, reads| b.iter(|| assemble_pipeline(black_box(reads), k)),
            );
            group.bench_with_input(
                BenchmarkId::new(format!("threading_build_{label}_k{k}"), reads.len()),
                reads,
                |b, reads| b.iter(|| threading_build_only(black_box(reads), k)),
            );
        }
    }
    group.finish();
}

fn bench_medium_depth_k10(c: &mut Criterion) {
    let reads = medium_depth_reads();
    c.bench_function("assembly_graph_medium_depth_k10", |b| {
        b.iter(|| assemble_pipeline(black_box(&reads), 10))
    });
}

fn bench_high_depth_k10(c: &mut Criterion) {
    let reads = high_depth_reads();
    c.bench_function("assembly_graph_high_depth_k10", |b| {
        b.iter(|| assemble_pipeline(black_box(&reads), 10))
    });
}

fn bench_kmer_key_ab(c: &mut Criterion) {
    let mut group = c.benchmark_group("kmer_key_representation");
    // Concatenate medium-depth reads into one long stream (many overlapping windows).
    let reads = medium_depth_reads();
    let mut seq = Vec::with_capacity(reads.len() * 130);
    for r in &reads {
        seq.extend_from_slice(&r.bases);
    }
    for k in [10usize, 25] {
        group.bench_with_input(BenchmarkId::new("A_arc_bytes", k), &seq, |b, seq| {
            b.iter(|| fill_arc_map(black_box(seq), k))
        });
        group.bench_with_input(BenchmarkId::new("B_packed_key", k), &seq, |b, seq| {
            b.iter(|| fill_packed_map(black_box(seq), k))
        });
    }
    group.finish();
}

criterion_group!(
    name = benches;
    config = smoke_criterion();
    targets = bench_medium_depth_k10, bench_high_depth_k10, bench_depth_matrix, bench_kmer_key_ab
);
criterion_main!(benches);
