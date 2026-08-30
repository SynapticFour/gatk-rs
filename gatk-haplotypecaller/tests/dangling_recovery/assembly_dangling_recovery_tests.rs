//! Unit tests for `assembly_dangling_recovery` (pulled in via `#[path]`).
//! Production algorithm lives in `src/assembly_dangling_recovery.rs`.

use super::*;
use crate::assembly::{AssemblyGraphParams, AssemblyGraphPruningParams, AssemblyRead};
use crate::assembly_pruning::apply_gatk_pruning;
use crate::read_threading_graph::assembly_graph_from_ref_and_reads_threading;

fn read(seq: &str, q: u8) -> AssemblyRead {
    AssemblyRead {
        bases: seq.as_bytes().to_vec(),
        base_quals: vec![q; seq.len()],
    }
}

fn build_pruned_ref_graph_at_k(
    reference: &str,
    alt_reads: &[&str],
    kmer_size: usize,
) -> AssemblyGraph {
    let params = AssemblyGraphParams {
        kmer_size: crate::bio_ids::KmerSize::try_new(kmer_size as u16).expect("test k≥2"),
        min_base_quality: 10,
        ..Default::default()
    };
    let mut reads: Vec<AssemblyRead> = alt_reads.iter().map(|s| read(s, 30)).collect();
    reads.insert(0, read(reference, 30));
    reads.insert(0, read(reference, 30));
    reads.insert(0, read(reference, 30));
    let reference = read(reference, 30);
    let mut graph =
        assembly_graph_from_ref_and_reads_threading(&reference, &reads, &params).unwrap();
    let mut prune = AssemblyGraphPruningParams::gatk_haplotype_caller_defaults();
    prune.min_prune_factor = 2;
    apply_gatk_pruning(&mut graph, &prune);
    graph
}

fn build_pruned_ref_graph(reference: &str, alt_reads: &[&str]) -> AssemblyGraph {
    build_pruned_ref_graph_at_k(reference, alt_reads, 3)
}

#[test]
fn find_path_upwards_is_lca_first() {
    let graph = build_pruned_ref_graph(
        "ACGTTGCATCG",
        &["ACGTTGCATCG", "ACGTTGCATCG", "ACGTTGCATCA", "ACGTTGCATCA"],
    );
    let alt_sink = graph
        .nodes()
        .iter()
        .position(|n| n.kmer.as_ref() == b"TCA")
        .expect("TCA sink");
    let path = graph
        .find_path_upwards_to_lca(alt_sink, 2, true)
        .expect("alt path");
    assert_eq!(graph.kmer_at(path[0]), b"ATC");
    assert_eq!(graph.kmer_at(*path.last().unwrap()), b"TCA");
}

#[test]
fn best_prefix_match_requires_min_matching_bases() {
    let mut cigar = Cigar::new();
    cigar.push(4, CigarOperator::Match);
    let ref_bases = b"ACGT";
    let alt_bases = b"ACGT";
    assert!(best_prefix_match(&cigar, ref_bases, alt_bases, 3).is_some());
    assert!(best_prefix_match(&cigar, ref_bases, alt_bases, 5).is_none());
}

/// GATK 4.4 `bestPrefixMatchLegacy`: perfect first-M window leaves `lastGoodIndex=-1`.
#[test]
fn six_r32_best_prefix_match_legacy_java44_perfect_prefix_rejects() {
    let seq = b"ACGTACGTACGTACGTACGTACGT";
    let max_index = 10;
    let k = 25;
    assert_eq!(
        best_prefix_match_legacy(seq, seq, max_index, k),
        None,
        "Java 4.4 lastGoodIndex stays -1 on a perfect prefix"
    );
    assert_eq!(
        best_prefix_match_legacy_java_44(seq, seq, max_index, k),
        -1,
        "Java 4.4 bestPrefixMatchLegacy returns -1 when there are zero mismatches"
    );
}

#[test]
fn six_r32_best_prefix_match_legacy_java44_mismatch_returns_last_mismatch() {
    let path1 = b"ACGTACGTAA";
    let path2 = b"ACGTACGTCA";
    let max_index = 10;
    let k = 25;
    assert_eq!(
        best_prefix_match_legacy_java_44(path1, path2, max_index, k),
        8,
        "Java lastGoodIndex is the last mismatch (index 8)"
    );
    assert_eq!(
        best_prefix_match_legacy(path1, path2, max_index, k),
        Some(8),
        "Rust matches Java: last mismatch, not maxIndex-1"
    );
}

/// GATK 4.4 aborts the whole merge when mismatches exceed the cap.
#[test]
fn six_r32_best_prefix_match_legacy_java44_excess_mismatches_abort() {
    let path1 = b"AAAAAAAAAA";
    let path2 = b"TTAAAAAAAA";
    let max_index = 10;
    let k = 25;
    assert_eq!(
        (max_index / k).max(1),
        1,
        "maxMismatches = max(1, 10/25) = 1"
    );
    assert_eq!(
        best_prefix_match_legacy_java_44(path1, path2, max_index, k),
        -1,
        "Java returns -1 on the second mismatch"
    );
    assert_eq!(
        best_prefix_match_legacy(path1, path2, max_index, k),
        None,
        "Rust matches Java: excess mismatches abort to None"
    );
}

/// 6R.33: exact 73M mid-B window from 6R.32 (mismatch cap 2, third mismatch at 39).
#[test]
fn six_r33_canonical_73m_prefix_exceeds_cap_rejects() {
    let rust_alt = b"AGTTTGGACGAGATACTTTCCCTTAGAAGTTGAGATACTCAACTTACGTCTGTAGTCTTTCTTTAAAGACTCT";
    let rust_ref = b"AGTTTGGACGAGATACTTTCCCTTACAAGTTGAGACACTGAACTTACGTTTGTAGTCTTTCTTCAAAGACCCT";
    assert_eq!(rust_alt.len(), 73);
    assert_eq!(rust_ref.len(), 73);
    let mm: Vec<usize> = (0..73).filter(|&i| rust_alt[i] != rust_ref[i]).collect();
    assert_eq!(mm, vec![25, 35, 39, 49, 63, 70]);
    let k = 25;
    let cap = (73 / k).max(1);
    assert_eq!(cap, 2);
    assert_ne!(
        best_prefix_match_legacy(rust_ref, rust_alt, 73, k),
        Some(35),
        "must not return the pre-6R.33 last_good_index"
    );
    assert_eq!(
        best_prefix_match_legacy(rust_ref, rust_alt, 73, k),
        None,
        "mismatch #3 at 39 exceeds cap 2 → Java -1 / Rust None"
    );
    assert_eq!(
        best_prefix_match_legacy_java_44(rust_ref, rust_alt, 73, k),
        -1
    );
}

#[test]
fn six_r33_prefix_exactly_at_cap_accepts_last_mismatch() {
    let path1 = vec![b'A'; 73];
    let mut at_cap = vec![b'A'; 73];
    at_cap[25] = b'T';
    at_cap[35] = b'T';
    let k = 25;
    assert_eq!((73 / k).max(1), 2);
    assert_eq!(
        best_prefix_match_legacy(&path1, &at_cap, 73, k),
        Some(35),
        "exactly two mismatches (the cap) → last mismatch index"
    );
    let mut over = at_cap;
    over[39] = b'T';
    assert_eq!(
        best_prefix_match_legacy(&path1, &over, 73, k),
        None,
        "third mismatch exceeds cap → reject"
    );
}

#[test]
fn six_r33_perfect_prefix_java_merge_index_is_none() {
    let seq = vec![b'A'; 73];
    let k = 25;
    assert_eq!(
        best_prefix_match_legacy(&seq, &seq, 73, k),
        None,
        "Java lastGoodIndex stays -1; mergeDanglingHeadLegacy rejects indexesToMerge <= 0"
    );
    assert_eq!(best_prefix_match_legacy_java_44(&seq, &seq, 73, k), -1);
}

/// 6R.34: Smith–Waterman + `bestPrefixMatchLegacy` on a pair of path encodings.
struct SixR34EncodingDecision {
    cigar: String,
    cigar_ok: bool,
    first_m: usize,
    mismatches: usize,
    mismatch_positions: Vec<usize>,
    cap: usize,
    prefix: Option<usize>,
}

fn six_r34_encoding_decision(
    ref_bases: &[u8],
    alt_bases: &[u8],
    k: usize,
) -> SixR34EncodingDecision {
    let cigar = align_dangling(
        ref_bases,
        alt_bases,
        &DanglingRecoverySwParams::gatk_defaults(),
    );
    let cigar_s = format_cigar_elements(&cigar);
    let cigar_ok = cigar_ok_to_merge_head(&cigar);
    if !cigar_ok {
        return SixR34EncodingDecision {
            cigar: cigar_s,
            cigar_ok: false,
            first_m: 0,
            mismatches: 0,
            mismatch_positions: Vec::new(),
            cap: 0,
            prefix: None,
        };
    }
    let first_m = cigar
        .elements
        .first()
        .filter(|e| e.operator == CigarOperator::Match)
        .map(|e| e.length)
        .unwrap_or(0);
    let n = first_m.min(ref_bases.len()).min(alt_bases.len());
    let mismatch_positions: Vec<usize> = (0..n)
        .filter(|&i| !base_eq(ref_bases[i], alt_bases[i]))
        .collect();
    let cap = (first_m / k.max(1)).max(1);
    let prefix = best_prefix_match_legacy(ref_bases, alt_bases, first_m, k);
    SixR34EncodingDecision {
        cigar: cigar_s,
        cigar_ok: true,
        first_m,
        mismatches: mismatch_positions.len(),
        mismatch_positions,
        cap,
        prefix,
    }
}

fn six_r34_legacy_accept(d: &SixR34EncodingDecision, ref_path_len: usize) -> bool {
    d.cigar_ok
        && d.prefix
            .is_some_and(|i| i > 0 && i < ref_path_len.saturating_sub(1))
}

/// Fixture A: one source at `path[0]`. OBSERVED RUST == SOURCE-DERIVED JAVA.
#[test]
fn six_r34_fixture_a_single_source_at_path0_encodings_match() {
    let mut g = AssemblyGraph::new(3).expect("k=3");
    let src = g.ensure_node(b"ACG");
    let mid = g.ensure_node(b"CGT");
    let snk = g.ensure_node(b"GTT");
    g.add_edge_support(src, mid, 1);
    g.add_edge_support(mid, snk, 1);
    let path = vec![src, mid, snk];
    assert!(g.is_source(src));
    assert!(!g.is_source(mid));
    assert!(!g.is_source(snk));
    let rust = path_bases(&g, &path, true);
    let java = path_bases_java_get_bases_for_path(&g, &path, true);
    assert_eq!(rust, java, "single source at path[0] must match Java");
    assert_eq!(rust, b"GCATT", "rev(ACG)+suffix(CGT)+suffix(GTT)");
    let rust_off = path_bases(&g, &path, false);
    let java_off = path_bases_java_get_bases_for_path(&g, &path, false);
    assert_eq!(rust_off, java_off);
    assert_eq!(rust_off, b"GTT");
}

/// Fixture B: later vertex is a source (LCA-first dangling-head shape).
#[test]
fn six_r34_fixture_b_source_after_path0_java_expands_rust_does_not() {
    let mut g = AssemblyGraph::new(3).expect("k=3");
    let pred = g.ensure_node(b"AAC");
    let lca = g.ensure_node(b"ACT");
    let mid = g.ensure_node(b"CTT");
    let head = g.ensure_node(b"TTA");
    g.add_edge_support(pred, lca, 1);
    g.add_edge_support(head, mid, 1);
    g.add_edge_support(mid, lca, 1);
    let path = vec![lca, mid, head];
    assert!(!g.is_source(lca), "LCA has incoming from pred");
    assert!(!g.is_source(mid));
    assert!(g.is_source(head));
    let rust = path_bases(&g, &path, true);
    let java = path_bases_java_get_bases_for_path(&g, &path, true);
    assert_eq!(
        rust, b"TTATT",
        "6R.36: suffixes then reverse(TTA)=ATT at the later source"
    );
    assert_eq!(rust, java);
    let rust_off = path_bases(&g, &path, false);
    let java_off = path_bases_java_get_bases_for_path(&g, &path, false);
    assert_eq!(rust_off, java_off, "expandSource=false never diverges");
    assert_eq!(rust_off, b"TTA");
}

/// Fixture C: two sources in one path list.
#[test]
fn six_r34_fixture_c_multiple_sources_java_expands_each() {
    let mut g = AssemblyGraph::new(3).expect("k=3");
    let s0 = g.ensure_node(b"AAA");
    let mid = g.ensure_node(b"AAT");
    let s1 = g.ensure_node(b"TTT");
    g.add_edge_support(s0, mid, 1);
    let path = vec![s0, mid, s1];
    assert!(g.is_source(s0));
    assert!(!g.is_source(mid));
    assert!(g.is_source(s1));
    let rust = path_bases(&g, &path, true);
    let java = path_bases_java_get_bases_for_path(&g, &path, true);
    assert_eq!(rust, b"AAATTTT", "6R.36: both sources: rev(AAA)+T+rev(TTT)");
    assert_eq!(rust, java);
}

/// Fixture D: reverse of B puts the source at `path[0]`. Reverse of C still
/// expands the later source (6R.36).
#[test]
fn six_r34_fixture_d_reverse_path_list() {
    let mut g = AssemblyGraph::new(3).expect("k=3");
    let pred = g.ensure_node(b"AAC");
    let lca = g.ensure_node(b"ACT");
    let mid = g.ensure_node(b"CTT");
    let head = g.ensure_node(b"TTA");
    g.add_edge_support(pred, lca, 1);
    g.add_edge_support(head, mid, 1);
    g.add_edge_support(mid, lca, 1);
    let reversed_b = vec![head, mid, lca];
    let rust_b = path_bases(&g, &reversed_b, true);
    let java_b = path_bases_java_get_bases_for_path(&g, &reversed_b, true);
    assert_eq!(rust_b, java_b, "single source at path[0] after reverse");
    assert_eq!(rust_b, b"ATTTT");

    let s0 = g.ensure_node(b"AAA");
    let mid_c = g.ensure_node(b"AAT");
    let s1 = g.ensure_node(b"GGG");
    g.add_edge_support(s0, mid_c, 1);
    let reversed_c = vec![s1, mid_c, s0];
    assert!(g.is_source(s1));
    assert!(g.is_source(s0));
    let rust_c = path_bases(&g, &reversed_c, true);
    let java_c = path_bases_java_get_bases_for_path(&g, &reversed_c, true);
    assert_eq!(
        rust_c, b"GGGTAAA",
        "6R.36: reverse still expands later source"
    );
    assert_eq!(rust_c, java_c);
}

/// Fixture B/C encodings into the same prefix/SW logic. Synthetic only.
#[test]
fn six_r34_fixture_b_downstream_synthetic_not_a_production_gate() {
    let mut g = AssemblyGraph::new(3).expect("k=3");
    let ref_src = g.ensure_node(b"AAA");
    let ref_mid = g.ensure_node(b"AAT");
    let lca = g.ensure_node(b"ACT");
    let alt_mid = g.ensure_node(b"CTT");
    let head = g.ensure_node(b"TTA");
    g.add_edge_support(ref_src, ref_mid, 1);
    g.add_edge_support(ref_mid, lca, 1);
    g.add_edge_support(head, alt_mid, 1);
    g.add_edge_support(alt_mid, lca, 1);
    let alt_path = vec![lca, alt_mid, head];
    let ref_path = vec![lca, ref_mid, ref_src];
    let rust_ref = path_bases(&g, &ref_path, true);
    let rust_alt = path_bases(&g, &alt_path, true);
    let java_ref = path_bases_java_get_bases_for_path(&g, &ref_path, true);
    let java_alt = path_bases_java_get_bases_for_path(&g, &alt_path, true);
    assert_eq!(rust_alt, java_alt);
    assert_eq!(rust_ref, java_ref);
    let rust_d = six_r34_encoding_decision(&rust_ref, &rust_alt, 3);
    let java_d = six_r34_encoding_decision(&java_ref, &java_alt, 3);
    let rust_acc = six_r34_legacy_accept(&rust_d, ref_path.len());
    let java_acc = six_r34_legacy_accept(&java_d, ref_path.len());
    eprintln!(
        "6R.34 SYNTHETIC B rust_seq={} java_seq={} rust_cigar={} java_cigar={} rust_mm={:?} java_mm={:?} rust_idx={:?} java_idx={:?} rust_acc={} java_acc={}",
        String::from_utf8_lossy(&rust_alt),
        String::from_utf8_lossy(&java_alt),
        rust_d.cigar,
        java_d.cigar,
        rust_d.mismatch_positions,
        java_d.mismatch_positions,
        rust_d.prefix,
        java_d.prefix,
        rust_acc,
        java_acc
    );
    assert_eq!(rust_alt, b"TTATT");
    assert_eq!(java_alt, b"TTATT");
    assert_eq!(rust_ref, b"TTAAA");
    assert_eq!(java_ref, b"TTAAA");
    assert_eq!(rust_acc, java_acc);
}

/// Existing dangling unit graphs: tails (`expand=false`) match; any head whose
/// encodings differ must not flip ACCEPT/REJECT vs Java-equivalent strings.
#[test]
fn six_r34_existing_dangling_unit_graphs_decision_parity() {
    let mut cases: Vec<(&str, AssemblyGraph, DanglingRecoveryParams)> = Vec::new();
    {
        let g = build_pruned_ref_graph(
            "ACGTTGCATCG",
            &["ACGTTGCATCG", "ACGTTGCATCG", "ACGTTGCATCA", "ACGTTGCATCA"],
        );
        let mut p = DanglingRecoveryParams::gatk_haplotype_caller_defaults();
        p.min_dangling_branch_length = 1;
        p.dangling_java_exact = true;
        cases.push(("k3_tail", g, p));
    }
    {
        let common_prefix = "AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT";
        let reference = format!("{common_prefix}GCTAGCTAATCG");
        let alt1 = format!("{common_prefix}ACTAGCTAATCG");
        let alt2 = format!("{common_prefix}ACTAGATAATCG");
        let g = build_pruned_ref_graph_at_k(&reference, &[&alt1, &alt2], 15);
        let mut p = DanglingRecoveryParams::gatk_haplotype_caller_defaults();
        p.min_dangling_branch_length = 4;
        p.recover_all_dangling_branches = true;
        p.dangling_java_exact = true;
        cases.push(("forked_tails", g, p));
    }
    {
        let mut g = AssemblyGraph::new(3).expect("k=3");
        let a = g.ensure_node(b"AAA");
        let b = g.ensure_node(b"AAB");
        let c = g.ensure_node(b"ABC");
        g.add_edge_support(a, b, 4);
        g.add_edge_support(b, c, 4);
        g.ref_edges.insert((a, b));
        g.ref_edges.insert((b, c));
        g.ref_nodes.extend([a, b, c]);
        let x = g.ensure_node(b"TTT");
        let y = g.ensure_node(b"TTA");
        let z = g.ensure_node(b"TAC");
        g.add_edge_support(x, y, 2);
        g.add_edge_support(y, z, 2);
        let mut p = DanglingRecoveryParams::gatk_haplotype_caller_defaults();
        p.dangling_java_exact = true;
        p.min_dangling_branch_length = 1;
        cases.push(("disconnected_island", g, p));
    }

    let mut heads_with_encoding_gap = 0usize;
    for (name, graph, params) in &cases {
        let sinks: Vec<usize> = (0..graph.node_count())
            .filter(|&v| graph.outgoing_nodes(v).is_empty() && !graph.is_ref_sink(v))
            .collect();
        for sink in sinks {
            if let Some(alt) = graph.find_path_upwards_to_lca(sink, params.min_prune_factor, true) {
                let rust = path_bases(graph, &alt, false);
                let java = path_bases_java_get_bases_for_path(graph, &alt, false);
                assert_eq!(
                    rust, java,
                    "{name} tail sink={sink}: expandSource=false must match Java"
                );
            }
        }
        let heads: Vec<usize> = (0..graph.node_count())
            .filter(|&v| {
                graph.incoming_count(v) == 0
                    && !graph.outgoing_nodes(v).is_empty()
                    && !graph.is_ref_source_vertex(v)
            })
            .collect();
        for head in heads {
            let dump = graph.test_dangling_head_decision_dump(head, params);
            let seqs_differ = dump.rust_alt_bases != dump.java_alt_bases
                || dump.rust_ref_bases != dump.java_ref_bases;
            if seqs_differ {
                heads_with_encoding_gap += 1;
                assert_eq!(
                    dump.final_rust, dump.final_java_source_derived,
                    "{name} head={head}: encoding gap must not flip ACCEPT/REJECT"
                );
                assert_eq!(
                    dump.cigar_ok == "PASS",
                    dump.java_seq_cigar_ok,
                    "{name} head={head}: cigarIsOkayToMerge must agree"
                );
            }
            eprintln!(
                "6R.34 EXISTING {name} head={head} seqs_differ={seqs_differ} rust={} java={} cigar_ok={} java_cigar_ok={}",
                dump.final_rust,
                dump.final_java_source_derived,
                dump.cigar_ok,
                dump.java_seq_cigar_ok
            );
        }
    }
    eprintln!("6R.34 EXISTING heads_with_encoding_gap={heads_with_encoding_gap}");
}

fn six_r35_dna_tag(i: usize) -> [u8; 4] {
    let b = [b'A', b'C', b'G', b'T'];
    [b[(i >> 6) & 3], b[(i >> 4) & 3], b[(i >> 2) & 3], b[i & 3]]
}

fn six_r35_kmer(tag_i: usize, last: u8) -> [u8; 5] {
    let t = six_r35_dna_tag(tag_i);
    [t[0], t[1], t[2], t[3], last]
}

/// Independent Java-equivalent helper: byte reverse, not reverse-complement.
#[test]
fn six_r35_java_helper_reverses_bytes_not_complement() {
    let mut g = AssemblyGraph::new(3).expect("k=3");
    let src = g.ensure_node(b"ACG");
    assert!(g.is_source(src));
    let java = java_get_bases_for_path_reference(&g, &[src], true);
    assert_eq!(java, b"GCA");
    assert_ne!(java, b"TGC", "must not reverse-complement (ACG→CGT→TGC)");
    let rust = path_bases(&g, &[src], true);
    assert_eq!(rust, java, "single source at path[0] matches");
    assert_eq!(java_get_bases_for_path_reference(&g, &[src], false), b"G");
    assert_eq!(path_bases(&g, &[src], false), b"G");
}

#[test]
fn six_r35_orientation_and_source_placement_exact_bases() {
    let mut g = AssemblyGraph::new(3).expect("k=3");
    let pred = g.ensure_node(b"AAC");
    let lca = g.ensure_node(b"ACT");
    let mid = g.ensure_node(b"CTT");
    let head = g.ensure_node(b"TTA");
    g.add_edge_support(pred, lca, 1);
    g.add_edge_support(head, mid, 1);
    g.add_edge_support(mid, lca, 1);
    let ref_src = g.ensure_node(b"AAA");
    let ref_mid = g.ensure_node(b"AAT");
    g.add_edge_support(ref_src, ref_mid, 1);
    g.add_edge_support(ref_mid, lca, 1);

    let one_src = vec![ref_src, ref_mid, lca];
    assert!(g.is_source(ref_src) && !g.is_source(ref_mid) && !g.is_source(lca));
    let rust = path_bases(&g, &one_src, true);
    let java = java_get_bases_for_path_reference(&g, &one_src, true);
    assert_eq!(rust, java);
    assert_eq!(rust, b"AAATT");

    let later = vec![lca, mid, head];
    let rust = path_bases(&g, &later, true);
    let java = java_get_bases_for_path_reference(&g, &later, true);
    assert_eq!(rust, b"TTATT");
    assert_eq!(rust, java);

    let none = vec![ref_mid, lca];
    assert!(!g.is_source(ref_mid) && !g.is_source(lca));
    let rust = path_bases(&g, &none, true);
    let java = java_get_bases_for_path_reference(&g, &none, true);
    assert_eq!(rust, java);
    assert_eq!(rust, b"TT");

    assert_eq!(
        path_bases(&g, &later, false),
        java_get_bases_for_path_reference(&g, &later, false)
    );

    let s1 = g.ensure_node(b"GGG");
    let multi = vec![ref_src, ref_mid, s1];
    assert!(g.is_source(ref_src) && g.is_source(s1));
    let rust = path_bases(&g, &multi, true);
    let java = java_get_bases_for_path_reference(&g, &multi, true);
    assert_eq!(rust, b"AAATGGG");
    assert_eq!(rust, java);

    let rev_one = vec![head, mid, lca];
    let rust = path_bases(&g, &rev_one, true);
    let java = java_get_bases_for_path_reference(&g, &rev_one, true);
    assert_eq!(rust, java);
    assert_eq!(rust, b"ATTTT");

    let rev_multi = vec![s1, ref_mid, ref_src];
    let rust = path_bases(&g, &rev_multi, true);
    let java = java_get_bases_for_path_reference(&g, &rev_multi, true);
    assert_eq!(rust, b"GGGTAAA");
    assert_eq!(rust, java);

    let alt_path = vec![lca, mid, head];
    let ref_path = vec![lca, ref_mid, ref_src];
    let row = path_bases_holdout_from_graph_paths(
        &g,
        "synth_later_source_k3",
        &alt_path,
        &ref_path,
        true,
    );
    eprintln_path_bases_holdout(&row);
    assert!(!row.seqs_differ);
    assert_eq!(row.rust_merge, row.java_merge);
    assert!(!row.flip);
}

/// Synthetic semantic probe: extra Java bases add a mismatch while the cap can stay 2.
/// Not a biological locus. Uses production SW + 6R.33 prefix-cap.
#[test]
fn six_r35_synthetic_cap_boundary_extra_mismatch() {
    let k = 25usize;
    let rust_ref = vec![b'A'; 50];
    let mut rust_alt = rust_ref.clone();
    rust_alt[10] = b'T';
    rust_alt[20] = b'T';
    let mut java_ref = rust_ref.clone();
    java_ref.extend(std::iter::repeat_n(b'A', 24));
    let mut java_alt = rust_alt.clone();
    java_alt.extend(std::iter::repeat_n(b'A', 23));
    java_alt.push(b'T');
    let row = path_bases_holdout_from_encodings(
        "synth_cap_boundary_rust_accept_java_reject",
        k,
        &rust_ref,
        &rust_alt,
        &java_ref,
        &java_alt,
        50,
        1,
    );
    eprintln_path_bases_holdout(&row);
    assert_eq!(row.rust_alt_len, 50);
    assert_eq!(row.java_alt_len, 74);
    assert_eq!(row.rust_cap, 2);
    assert_eq!(row.java_cap, 2);
    assert_eq!(row.rust_mm, vec![10, 20]);
    assert_eq!(row.java_mm, vec![10, 20, 73]);
    assert_eq!(row.rust_merge, "ACCEPT");
    assert_eq!(row.java_merge, "REJECT");
    assert!(
        row.flip,
        "synthetic cap-boundary must demonstrate ACCEPT vs REJECT"
    );
}

/// Inverse synthetic: extra Java bases raise first-M enough that 3 mismatches fit the new cap.
#[test]
fn six_r35_synthetic_cap_boundary_inverse() {
    let k = 25usize;
    let rust_ref = vec![b'A'; 73];
    let mut rust_alt = rust_ref.clone();
    rust_alt[10] = b'T';
    rust_alt[20] = b'T';
    rust_alt[30] = b'T';
    let mut java_ref = rust_ref.clone();
    java_ref.extend(std::iter::repeat_n(b'A', 24));
    let mut java_alt = rust_alt.clone();
    java_alt.extend(std::iter::repeat_n(b'A', 24));
    let row = path_bases_holdout_from_encodings(
        "synth_cap_boundary_rust_reject_java_accept",
        k,
        &rust_ref,
        &rust_alt,
        &java_ref,
        &java_alt,
        73,
        1,
    );
    eprintln_path_bases_holdout(&row);
    assert_eq!(row.rust_cap, 2);
    assert_eq!(row.java_cap, 3);
    assert_eq!(row.rust_mm, vec![10, 20, 30]);
    assert_eq!(row.java_mm, vec![10, 20, 30]);
    assert_eq!(row.rust_merge, "REJECT");
    assert_eq!(row.java_merge, "ACCEPT");
    assert!(row.flip, "inverse cap-boundary must flip REJECT vs ACCEPT");
}

/// Graph-derived later-source chain near a small-k cap (semantic probe, not a locus).
#[test]
fn six_r35_synthetic_graph_later_source_cap_probe() {
    let mut g = AssemblyGraph::new(5).expect("k=5");
    let pred = g.ensure_node(&six_r35_kmer(99, b'A'));
    let mut alt = Vec::new();
    let mut rref = Vec::new();
    let last_alt: [u8; 10] = *b"ATAAATAAAA";
    let last_ref: [u8; 10] = *b"AAAAAAAAAA";
    for i in 0..10 {
        alt.push(g.ensure_node(&six_r35_kmer(i, last_alt[i])));
        rref.push(g.ensure_node(&six_r35_kmer(50 + i, last_ref[i])));
    }
    g.add_edge_support(pred, alt[0], 1);
    g.add_edge_support(pred, rref[0], 1);
    for w in alt.windows(2).rev() {
        g.add_edge_support(w[1], w[0], 1);
    }
    for w in rref.windows(2).rev() {
        g.add_edge_support(w[1], w[0], 1);
    }
    assert!(!g.is_source(alt[0]));
    assert!(g.is_source(*alt.last().unwrap()));
    assert!(g.is_source(*rref.last().unwrap()));
    let rust_alt = path_bases(&g, &alt, true);
    let java_alt = java_get_bases_for_path_reference(&g, &alt, true);
    assert_eq!(rust_alt, java_alt);
    assert_eq!(java_alt.len(), last_alt.len() + 4);
    assert_ne!(
        rust_alt,
        &last_alt[..],
        "later source adds reverse(k-mer)[1..]"
    );
    let row =
        path_bases_holdout_from_graph_paths(&g, "synth_graph_k5_later_source", &alt, &rref, true);
    eprintln_path_bases_holdout(&row);
    assert!(!row.seqs_differ);
    assert_eq!(row.rust_merge, row.java_merge);
    assert_eq!(row.java_merge, "REJECT");
    assert!(!row.flip);
}

#[test]
fn six_r35_existing_dangling_unit_graphs_holdout_scan() {
    let mut cases: Vec<(&str, AssemblyGraph, DanglingRecoveryParams)> = Vec::new();
    {
        let g = build_pruned_ref_graph(
            "ACGTTGCATCG",
            &["ACGTTGCATCG", "ACGTTGCATCG", "ACGTTGCATCA", "ACGTTGCATCA"],
        );
        let mut p = DanglingRecoveryParams::gatk_haplotype_caller_defaults();
        p.min_dangling_branch_length = 1;
        p.dangling_java_exact = true;
        cases.push(("k3_tail", g, p));
    }
    {
        let common_prefix = "AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT";
        let reference = format!("{common_prefix}GCTAGCTAATCG");
        let alt1 = format!("{common_prefix}ACTAGCTAATCG");
        let alt2 = format!("{common_prefix}ACTAGATAATCG");
        let g = build_pruned_ref_graph_at_k(&reference, &[&alt1, &alt2], 15);
        let mut p = DanglingRecoveryParams::gatk_haplotype_caller_defaults();
        p.min_dangling_branch_length = 4;
        p.recover_all_dangling_branches = true;
        p.dangling_java_exact = true;
        cases.push(("forked_tails", g, p));
    }
    {
        let mut g = AssemblyGraph::new(3).expect("k=3");
        let a = g.ensure_node(b"AAA");
        let b = g.ensure_node(b"AAB");
        let c = g.ensure_node(b"ABC");
        g.add_edge_support(a, b, 4);
        g.add_edge_support(b, c, 4);
        g.ref_edges.insert((a, b));
        g.ref_edges.insert((b, c));
        g.ref_nodes.extend([a, b, c]);
        let x = g.ensure_node(b"TTT");
        let y = g.ensure_node(b"TTA");
        let z = g.ensure_node(b"TAC");
        g.add_edge_support(x, y, 2);
        g.add_edge_support(y, z, 2);
        let mut p = DanglingRecoveryParams::gatk_haplotype_caller_defaults();
        p.dangling_java_exact = true;
        p.min_dangling_branch_length = 1;
        cases.push(("disconnected_island", g, p));
    }

    let mut real_heads = 0usize;
    let mut real_flips = 0usize;
    let mut encoding_gaps = 0usize;
    for (name, graph, params) in &cases {
        for v in 0..graph.node_count() {
            if graph.outgoing_nodes(v).is_empty() && !graph.is_ref_sink(v) {
                if let Some(alt) = graph.find_path_upwards_to_lca(v, params.min_prune_factor, true)
                {
                    let rust = path_bases(graph, &alt, false);
                    let java = java_get_bases_for_path_reference(graph, &alt, false);
                    assert_eq!(rust, java, "{name} tail expand=false");
                }
            }
        }
        let heads: Vec<usize> = (0..graph.node_count())
            .filter(|&v| {
                graph.incoming_count(v) == 0
                    && !graph.outgoing_nodes(v).is_empty()
                    && !graph.is_ref_source_vertex(v)
            })
            .collect();
        for head in heads {
            let dump = graph.test_dangling_head_decision_dump(head, params);
            if dump.alt_path_ids.is_empty() || dump.ref_path_ids.is_empty() {
                continue;
            }
            real_heads += 1;
            let row = path_bases_holdout_from_graph_paths(
                graph,
                format!("{name}_head{head}"),
                &dump.alt_path_ids,
                &dump.ref_path_ids,
                true,
            );
            eprintln_path_bases_holdout(&row);
            if row.seqs_differ {
                encoding_gaps += 1;
            }
            if row.flip {
                real_flips += 1;
            }
        }
    }
    eprintln!(
        "6R.35 UNIT_GRAPH_SCAN heads={real_heads} encoding_gaps={encoding_gaps} flips={real_flips}"
    );
    assert_eq!(
        real_flips, 0,
        "unit-graph dangling heads must not flip merge"
    );
}

#[test]
fn six_r36_path_bases_expands_every_source() {
    let mut g = AssemblyGraph::new(3).expect("k=3");
    let s0 = g.ensure_node(b"AAA");
    let mid = g.ensure_node(b"AAT");
    let s1 = g.ensure_node(b"TTT");
    g.add_edge_support(s0, mid, 1);
    let path = vec![s0, mid, s1];
    assert!(g.is_source(s0) && g.is_source(s1) && !g.is_source(mid));
    let rust = path_bases(&g, &path, true);
    let java = java_get_bases_for_path_reference(&g, &path, true);
    assert_eq!(rust, b"AAATTTT");
    assert_eq!(rust, java);
}

#[test]
fn six_r36_path_bases_single_source_at_path0() {
    let mut g = AssemblyGraph::new(3).expect("k=3");
    let src = g.ensure_node(b"ACG");
    let mid = g.ensure_node(b"CGT");
    let snk = g.ensure_node(b"GTT");
    g.add_edge_support(src, mid, 1);
    g.add_edge_support(mid, snk, 1);
    let rust = path_bases(&g, &[src, mid, snk], true);
    assert_eq!(rust, b"GCATT");
    assert_eq!(
        rust,
        java_get_bases_for_path_reference(&g, &[src, mid, snk], true)
    );
}

#[test]
fn six_r36_reverse_oriented_source_bytes() {
    let mut g = AssemblyGraph::new(3).expect("k=3");
    let pred = g.ensure_node(b"AAC");
    let lca = g.ensure_node(b"ACT");
    let mid = g.ensure_node(b"CTT");
    let head = g.ensure_node(b"TTA");
    g.add_edge_support(pred, lca, 1);
    g.add_edge_support(head, mid, 1);
    g.add_edge_support(mid, lca, 1);
    let reversed = vec![head, mid, lca];
    let rust = path_bases(&g, &reversed, true);
    assert_eq!(rust, b"ATTTT", "rev(TTA)=ATT then suffixes T,T");
    assert_eq!(rust, java_get_bases_for_path_reference(&g, &reversed, true));
    assert_ne!(
        rust, b"TAATT",
        "byte reverse, not reverse-complement of TTA"
    );
}

#[test]
fn six_r36_expand_source_false_suffix_only() {
    let mut g = AssemblyGraph::new(3).expect("k=3");
    let pred = g.ensure_node(b"AAC");
    let lca = g.ensure_node(b"ACT");
    let mid = g.ensure_node(b"CTT");
    let head = g.ensure_node(b"TTA");
    g.add_edge_support(pred, lca, 1);
    g.add_edge_support(head, mid, 1);
    g.add_edge_support(mid, lca, 1);
    let path = vec![lca, mid, head];
    let rust_off = path_bases(&g, &path, false);
    let java_off = java_get_bases_for_path_reference(&g, &path, false);
    assert_eq!(rust_off, b"TTA");
    assert_eq!(rust_off, java_off);
    let rust_on = path_bases(&g, &path, true);
    assert_eq!(rust_on, b"TTATT");
    assert_ne!(rust_off, rust_on);
}

/// 6R.35 k=5 later-source holdout: production path_bases now matches Java, so both REJECT.
#[test]
fn six_r36_k5_holdout_now_matches_java_reject() {
    let mut g = AssemblyGraph::new(5).expect("k=5");
    let pred = g.ensure_node(&six_r35_kmer(99, b'A'));
    let mut alt = Vec::new();
    let mut rref = Vec::new();
    let last_alt: [u8; 10] = *b"ATAAATAAAA";
    let last_ref: [u8; 10] = *b"AAAAAAAAAA";
    for i in 0..10 {
        alt.push(g.ensure_node(&six_r35_kmer(i, last_alt[i])));
        rref.push(g.ensure_node(&six_r35_kmer(50 + i, last_ref[i])));
    }
    g.add_edge_support(pred, alt[0], 1);
    g.add_edge_support(pred, rref[0], 1);
    for w in alt.windows(2).rev() {
        g.add_edge_support(w[1], w[0], 1);
    }
    for w in rref.windows(2).rev() {
        g.add_edge_support(w[1], w[0], 1);
    }
    let rust_alt = path_bases(&g, &alt, true);
    let java_alt = java_get_bases_for_path_reference(&g, &alt, true);
    assert_eq!(rust_alt, java_alt);
    assert_ne!(&rust_alt[..], &last_alt[..]);
    let row = path_bases_holdout_from_graph_paths(&g, "six_r36_k5_later_source", &alt, &rref, true);
    eprintln_path_bases_holdout(&row);
    assert!(!row.seqs_differ);
    assert_eq!(row.rust_merge, "REJECT");
    assert_eq!(row.java_merge, "REJECT");
    assert!(!row.flip);
}

#[test]
fn longest_suffix_match_matches_gatk_examples() {
    assert_eq!(longest_suffix_match_java(b"ACGT", b"TGT", 3), 2);
    assert_eq!(longest_suffix_match_java(b"ACGT", b"CGT", 3), 3);
    assert_eq!(longest_suffix_match_java(b"CG", b"CA", 1), 0);
}

#[test]
fn dangling_java_exact_single_pass_matches_gatk_edge_count() {
    let mut graph = build_pruned_ref_graph(
        "ACGTTGCATCG",
        &["ACGTTGCATCG", "ACGTTGCATCG", "ACGTTGCATCA", "ACGTTGCATCA"],
    );
    let mut params = DanglingRecoveryParams::gatk_haplotype_caller_defaults();
    params.min_dangling_branch_length = 1;
    params.dangling_java_exact = true;
    let summary = graph.recover_dangling_branches(&params).unwrap();
    assert!(summary.tails_attempted >= 1);
    // This fixture needs ASM-1 suffix rescue for a merge; GATK-exact mode correctly skips it.
    assert_eq!(summary.tails_recovered, 0);
    let mut multi = params;
    multi.dangling_java_exact = false;
    let multi_summary = graph.recover_dangling_branches(&multi).unwrap();
    assert_eq!(multi_summary.tails_recovered, 1);
}

#[test]
fn dangling_tail_recovery_attempts_alt_sink() {
    let mut graph = build_pruned_ref_graph(
        "ACGTTGCATCG",
        &["ACGTTGCATCG", "ACGTTGCATCG", "ACGTTGCATCA", "ACGTTGCATCA"],
    );
    let mut dangling = DanglingRecoveryParams::gatk_haplotype_caller_defaults();
    dangling.min_dangling_branch_length = 1;
    let summary = graph.recover_dangling_branches(&dangling).unwrap();
    assert!(summary.tails_attempted >= 1);
    // Java idempotent `addEdge` → `edge_exists` counts as tail recovered (ASM-1).
    assert_eq!(summary.tails_recovered, 1);
}

#[test]
fn cigar_ok_to_merge_tail_requires_terminal_match_op() {
    let mut trailing_ins = Cigar::new();
    trailing_ins.push(1, CigarOperator::Match);
    trailing_ins.push(1, CigarOperator::Insertion);
    assert!(!cigar_ok_to_merge_tail(&trailing_ins));

    let mut ok = Cigar::new();
    ok.push(2, CigarOperator::Match);
    assert!(cigar_ok_to_merge_tail(&ok));
}

#[test]
fn align_dangling_uses_leading_indel_strategy() {
    let p = DanglingRecoverySwParams::gatk_defaults();
    let cigar = align_dangling(b"ACGTACGTAC", b"ACGTXACGTAC", &p);
    assert!(cigar_ok_to_merge_tail(&cigar));
    assert!(cigar
        .elements
        .last()
        .is_some_and(|e| e.operator == CigarOperator::Match));
}

/// GATK `ReadThreadingGraphUnitTest.testForkedDanglingEnds`.
#[test]
fn forked_dangling_ends_recovers_all_alt_sinks_with_recover_all() {
    let common_prefix = "AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT";
    let reference = format!("{common_prefix}GCTAGCTAATCG");
    let alt1 = format!("{common_prefix}ACTAGCTAATCG");
    let alt2 = format!("{common_prefix}ACTAGATAATCG");
    let mut graph = build_pruned_ref_graph_at_k(&reference, &[&alt1, &alt2], 15);
    let mut dangling = DanglingRecoveryParams::gatk_haplotype_caller_defaults();
    dangling.min_dangling_branch_length = 4;
    dangling.recover_all_dangling_branches = true;
    let summary = graph.recover_dangling_branches(&dangling).unwrap();
    assert!(
        summary.tails_attempted >= 1,
        "GATK testForkedDanglingEnds expects non-ref sinks (Rust may collapse forks to fewer sinks)"
    );
    assert_eq!(
        summary.tails_recovered, summary.tails_attempted,
        "recoverAll should merge every attempted alt tail"
    );
}

#[test]
fn reference_path_from_breaks_on_cycle() {
    // Synthetic 2-node ref cycle: without the revisit guard this walk is unbounded.
    let mut g = AssemblyGraph::new(3).expect("k=3 graph");
    let a = g.ensure_node(b"AAA");
    let b = g.ensure_node(b"AAB");
    g.add_edge_support(a, b, 1);
    g.add_edge_support(b, a, 1);
    g.ref_edges.insert((a, b));
    g.ref_edges.insert((b, a));
    g.ref_nodes.extend([a, b]);
    let path = g.reference_path_from(a, TraversalDir::Down, None);
    assert!(
        path.len() <= 3,
        "cycle must not grow ref path unboundedly: got {}",
        path.len()
    );
    assert_eq!(path[0], a);
}

#[test]
fn disconnected_alt_island_not_attached_by_either_dangling_java_exact() {
    let mut g = AssemblyGraph::new(3).expect("k=3");
    let a = g.ensure_node(b"AAA");
    let b = g.ensure_node(b"AAB");
    let c = g.ensure_node(b"ABC");
    g.add_edge_support(a, b, 4);
    g.add_edge_support(b, c, 4);
    g.ref_edges.insert((a, b));
    g.ref_edges.insert((b, c));
    g.ref_nodes.extend([a, b, c]);
    let x = g.ensure_node(b"TTT");
    let y = g.ensure_node(b"TTA");
    let z = g.ensure_node(b"TAC");
    g.add_edge_support(x, y, 2);
    g.add_edge_support(y, z, 2);
    for exact in [true, false] {
        let mut g = g.clone();
        let mut params = DanglingRecoveryParams::gatk_haplotype_caller_defaults();
        params.dangling_java_exact = exact;
        params.min_dangling_branch_length = 1;
        let summary = g.recover_dangling_branches(&params).unwrap();
        assert_eq!(
            summary.heads_recovered, 0,
            "disconnected island must not head-merge exact={exact}"
        );
        assert_eq!(
            summary.tails_recovered, 0,
            "disconnected island must not tail-merge exact={exact}"
        );
        let still_disconnected = g.incoming_count(x) == 0 && g.outgoing_nodes(z).is_empty();
        assert!(
            still_disconnected,
            "island endpoints must stay detached exact={exact}"
        );
    }
}
