//! Microbenchmarks for Phase-5 assembly-graph hot paths (step 73).
#![allow(clippy::result_large_err)]

use std::hint::black_box;
use std::time::Duration;

use criterion::{criterion_group, criterion_main, Criterion};
use gatk_haplotypecaller::{
    AssemblyGraph, AssemblyGraphParams, AssemblyGraphPruningParams, AssemblyRead, KmerSize,
};

fn smoke_criterion() -> Criterion {
    Criterion::default()
        .warm_up_time(Duration::from_millis(200))
        .measurement_time(Duration::from_secs(1))
        .sample_size(10)
}

fn mk_read(seq: &str, q: u8) -> AssemblyRead {
    AssemblyRead {
        bases: seq.as_bytes().to_vec(),
        base_quals: vec![q; seq.len()],
    }
}

fn medium_depth_reads() -> Vec<AssemblyRead> {
    let mut reads = Vec::with_capacity(64);
    for i in 0..64 {
        let mut seq = String::from("ACGTACGTAC");
        if i % 3 == 0 {
            seq.push_str("TTAA");
        } else if i % 3 == 1 {
            seq.push_str("TCAA");
        } else {
            seq.push_str("TGAA");
        }
        seq.push_str("GCTAGCTAGC");
        reads.push(mk_read(&seq, 30));
    }
    reads
}

fn high_depth_reads() -> Vec<AssemblyRead> {
    let mut reads = Vec::with_capacity(512);
    for i in 0..512 {
        let mut seq = String::from("ACGTACGT");
        if i % 7 == 0 {
            seq.push_str("TTAA");
        } else if i % 7 == 1 {
            seq.push_str("TCAA");
        } else {
            seq.push_str("TGAA");
        }
        seq.push_str("AAAAAAAAAAAAAAAA");
        reads.push(mk_read(&seq, 30));
    }
    reads
}

fn assemble_pipeline(reads: &[AssemblyRead], k: u16) {
    let params = AssemblyGraphParams {
        kmer_size: KmerSize::try_new(k).expect("k"),
        min_base_quality: 10,
        min_edge_weight: 2,
        dangling_path_max_nodes: 4,
        max_haplotypes: 32,
        max_haplotype_bases: 128,
        ..Default::default()
    };
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

criterion_group!(
    name = benches;
    config = smoke_criterion();
    targets = bench_medium_depth_k10, bench_high_depth_k10
);
criterion_main!(benches);
