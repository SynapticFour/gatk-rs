//! Phase 5 assembly regression + determinism contracts (Steps 74–75 / P0 matrix).
//!
//! Invoked by:
//! - `cargo test -p gatk-haplotypecaller --test p5_assembly_regression_test`
//! - `run_p5_determinism_matrix.sh` / `run_p5_assembly_stability_contract.sh`
//!   filter `outputs_are_stable_across_repeated_runs_and_input_order`

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

fn case1_params() -> AssemblyGraphParams {
    // Same contract as `p5_haplotype_candidate_diff_test` / frozen Java export.
    AssemblyGraphParams {
        kmer_size: KmerSize::try_new(3).expect("k=3"),
        min_base_quality: 10,
        min_edge_weight: 1,
        max_haplotypes: 16,
        max_haplotype_bases: 64,
        dangling_path_max_nodes: 0,
        ..Default::default()
    }
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
                .to_string();
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

fn mk_read(seq: &str, q: u8) -> AssemblyRead {
    AssemblyRead {
        bases: seq.to_string(),
        base_quals: vec![q; seq.len()],
    }
}

fn assemble_candidates(reads: &[AssemblyRead], params: &AssemblyGraphParams) -> Vec<(String, u32)> {
    let mut graph = AssemblyGraph::from_reads(reads, params).expect("build graph");
    let mut pruning = AssemblyGraphPruningParams::gatk_haplotype_caller_defaults();
    pruning.min_prune_factor = params.min_edge_weight;
    graph.apply_pruning(&pruning);
    graph.remove_dangling_paths(params.dangling_path_max_nodes);
    graph.cleanup_isolated_nodes();
    graph
        .extract_candidate_haplotypes(params.max_haplotypes, params.max_haplotype_bases)
        .into_iter()
        .map(|h| (h.sequence, h.support))
        .collect()
}

fn case1_reads() -> Vec<AssemblyRead> {
    let path = repo_root().join("parity/fixtures/p5_assembly_case1_reads.tsv");
    let reads = load_reads_tsv(&path);
    assert!(!reads.is_empty(), "case1 reads fixture must be non-empty");
    reads
}

/// P0 / Step-75: identical candidates across repeats, input order, and Rayon thread env.
#[test]
fn outputs_are_stable_across_repeated_runs_and_input_order() {
    let params = case1_params();
    let reads = case1_reads();
    let baseline = assemble_candidates(&reads, &params);
    assert!(!baseline.is_empty(), "case1 must produce candidates");

    for _ in 0..3 {
        assert_eq!(assemble_candidates(&reads, &params), baseline);
    }

    let mut reversed = reads.clone();
    reversed.reverse();
    assert_eq!(
        assemble_candidates(&reversed, &params),
        baseline,
        "candidate set must be stable under input reversal"
    );

    let mut rotated = reads.clone();
    rotated.rotate_left(reads.len() / 2);
    assert_eq!(
        assemble_candidates(&rotated, &params),
        baseline,
        "candidate set must be stable under input rotation"
    );

    // Harness varies RAYON_NUM_THREADS; assembly output must not depend on it.
    let _ = std::env::var("RAYON_NUM_THREADS");
    assert_eq!(assemble_candidates(&reads, &params), baseline);
}

#[test]
fn case1_matches_frozen_java_candidate_set() {
    let params = case1_params();
    let reads = case1_reads();
    let got: BTreeSet<(String, u32)> = assemble_candidates(&reads, &params).into_iter().collect();
    let expected_path = repo_root().join("parity/expected/p5_assembly_case1.java_candidates.tsv");
    let expected: BTreeSet<(String, u32)> = fs::read_to_string(&expected_path)
        .unwrap_or_else(|e| panic!("read {}: {e}", expected_path.display()))
        .lines()
        .filter(|l| {
            let t = l.trim();
            !t.is_empty() && !t.starts_with('#')
        })
        .map(|l| {
            let mut parts = l.split_whitespace();
            let seq = parts.next().unwrap().to_string();
            let support = parts.next().unwrap().parse::<u32>().unwrap();
            (seq, support)
        })
        .collect();
    assert_eq!(got, expected);
}

#[test]
fn alt_path_and_homopolymer_motifs_remain_recoverable() {
    let snp_params = AssemblyGraphParams {
        kmer_size: KmerSize::try_new(3).unwrap(),
        min_base_quality: 10,
        min_edge_weight: 2,
        dangling_path_max_nodes: 0,
        max_haplotypes: 16,
        ..Default::default()
    };
    let snp_reads = vec![
        mk_read("ACGTT", 30),
        mk_read("ACGTT", 30),
        mk_read("ACGTT", 30),
        mk_read("ACGTA", 30),
        mk_read("ACGTA", 30),
    ];
    let snp = assemble_candidates(&snp_reads, &snp_params);
    assert!(snp.iter().any(|(s, _)| s.ends_with('T')));
    assert!(snp.iter().any(|(s, _)| s.ends_with('A')));

    let hp_params = AssemblyGraphParams {
        kmer_size: KmerSize::try_new(4).unwrap(),
        min_base_quality: 10,
        min_edge_weight: 1,
        dangling_path_max_nodes: 0,
        max_haplotypes: 8,
        ..Default::default()
    };
    let hp_reads = vec![
        mk_read("AAAACAAA", 30),
        mk_read("AAAACAAA", 30),
        mk_read("AAAAGAAA", 30),
        mk_read("AAAAGAAA", 30),
    ];
    let hp = assemble_candidates(&hp_reads, &hp_params);
    assert!(hp.iter().any(|(s, _)| s.contains("AAAAC")));
    assert!(hp.iter().any(|(s, _)| s.contains("AAAAG")));
}

#[test]
fn equal_support_tie_break_is_lexicographic() {
    let params = AssemblyGraphParams {
        kmer_size: KmerSize::try_new(3).unwrap(),
        min_base_quality: 10,
        min_edge_weight: 1,
        dangling_path_max_nodes: 0,
        max_haplotypes: 16,
        ..Default::default()
    };
    let reads = vec![
        mk_read("ACGTA", 30),
        mk_read("ACGTA", 30),
        mk_read("ACGTC", 30),
        mk_read("ACGTC", 30),
    ];
    let hs = assemble_candidates(&reads, &params);
    let eq = hs
        .windows(2)
        .find(|w| w[0].1 == w[1].1)
        .expect("equal-support pair");
    assert!(eq[0].0 <= eq[1].0);
}
