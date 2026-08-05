//! Phase-5 Step-70: frozen Java-export candidate sets vs Rust local-assembly output.
use gatk_haplotypecaller::{
    AssemblyGraph, AssemblyGraphParams, AssemblyGraphPruningParams, AssemblyRead, KmerSize,
};
use std::collections::BTreeSet;
use std::fs;
use std::path::PathBuf;

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("workspace root")
        .to_path_buf()
}

fn load_reads_tsv(path: &PathBuf) -> Vec<AssemblyRead> {
    fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("read {}: {e}", path.display()))
        .lines()
        .filter(|l| {
            let t = l.trim();
            !t.is_empty() && !t.starts_with('#')
        })
        .map(|l| {
            let mut parts = l.split_whitespace();
            let bases = parts
                .next()
                .unwrap_or_else(|| panic!("missing sequence in {}", path.display()))
                .as_bytes()
                .to_vec();
            let q = parts
                .next()
                .unwrap_or_else(|| panic!("missing qual in {}", path.display()))
                .parse::<u8>()
                .unwrap_or_else(|_| panic!("invalid qual in {}", path.display()));
            let n = bases.len();
            AssemblyRead {
                bases,
                base_quals: vec![q; n],
            }
        })
        .collect()
}

fn load_candidates(path: &PathBuf) -> BTreeSet<(String, u32)> {
    fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("read {}: {e}", path.display()))
        .lines()
        .filter(|l| {
            let t = l.trim();
            !t.is_empty() && !t.starts_with('#')
        })
        .map(|l| {
            let mut parts = l.split_whitespace();
            let seq = parts
                .next()
                .unwrap_or_else(|| panic!("missing sequence in {}", path.display()))
                .to_string();
            let support = parts
                .next()
                .unwrap_or_else(|| panic!("missing support in {}", path.display()))
                .parse::<u32>()
                .unwrap_or_else(|_| panic!("invalid support in {}", path.display()));
            (seq, support)
        })
        .collect()
}

#[test]
fn p5_assembly_case1_candidates_match_frozen_java_export() {
    let root = repo_root();
    let reads_path = root.join("parity/fixtures/p5_assembly_case1_reads.tsv");
    let expected_path = root.join("parity/expected/p5_assembly_case1.java_candidates.tsv");
    let reads = load_reads_tsv(&reads_path);
    assert!(
        !reads.is_empty(),
        "p5_assembly_case1_reads.tsv must be non-empty"
    );

    // Same k=3 / min-qual=10 contract as the frozen Java-export fixture (support ACGTT=8, ACGTA=7).
    let params = AssemblyGraphParams {
        kmer_size: KmerSize::try_new(3).expect("k=3"),
        min_base_quality: 10,
        min_edge_weight: 1,
        max_haplotypes: 16,
        max_haplotype_bases: 64,
        dangling_path_max_nodes: 0,
        ..Default::default()
    };
    let mut graph = AssemblyGraph::from_reads(&reads, &params).expect("build graph");
    let mut pruning = AssemblyGraphPruningParams::gatk_haplotype_caller_defaults();
    pruning.min_prune_factor = params.min_edge_weight;
    graph.apply_pruning(&pruning);
    graph.remove_dangling_paths(params.dangling_path_max_nodes);
    graph.cleanup_isolated_nodes();
    let got = graph.extract_candidate_haplotypes(params.max_haplotypes, params.max_haplotype_bases);

    let actual: BTreeSet<(String, u32)> = got
        .into_iter()
        .map(|h| (String::from_utf8_lossy(&h.sequence).into_owned(), h.support))
        .collect();
    let expected = load_candidates(&expected_path);

    assert_eq!(
        actual, expected,
        "Rust local-assembly candidates must match frozen Java-export TSV"
    );
}
