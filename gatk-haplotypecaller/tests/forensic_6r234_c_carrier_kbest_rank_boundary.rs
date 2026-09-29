//! 6R.234: why Java ranks carrier `c7acc50dfb9f9ecc` at 116 while Rust
//! ranks the same FNV at 129. Target site remains `20:29455379 G/A`.
//!
//! Proof-only. Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! PRODUCTION CHANGE: NONE. Do not raise K.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r234_c_carrier_kbest_rank_boundary -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::assembly_based_caller::assemble_reads_with_finalized;
use gatk_haplotypecaller::assembly_region_finalize::{
    assembly_reference_read, create_graph_reference_read, records_to_assembly_reads,
};
use gatk_haplotypecaller::read_threading_assembler::{
    build_threading_graph_for_seq_assembly, DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH,
};
use gatk_haplotypecaller::seq_graph::SeqGraph;
use gatk_haplotypecaller::seq_kbest_haplotype::{
    find_best_haplotypes_seq_graph, seq_kbest_path_score_terms, SeqKbestEdgeTerm,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    AssemblyRegionCallDisposition, CallRegionArgs, ReadFilterParams, WalkerTraversalConfig,
};
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 29_455_379;
const EVENT_POS: u64 = 29_455_314;
const JAVA_CARRIER: &str = "c7acc50dfb9f9ecc";
const JAVA_CARRIER_RANK: usize = 116;
const JAVA_CARRIER_SCORE: f64 = -2.647_867_02;
const RUST_DIAG_RANK: usize = 129;
const RUST_DIAG_SCORE: f64 = -2.740_621_07;
const JAVA_KBEST_K: usize = 128;
const JAVA_TSV: &str = include_str!("forensic_6r234_java_kbest.tsv");

#[derive(Clone)]
struct JavaEdge {
    i: usize,
    from: usize,
    to: usize,
    mult: u32,
    out: u32,
    is_ref: bool,
    penalty: f64,
    to_seq: String,
}

#[derive(Clone)]
struct RankRow {
    rank: usize,
    hash: String,
    score: f64,
    is_ref: bool,
    len: usize,
}

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R234\t{key}\t{}", value.as_ref());
}

fn fnv1a64_hex(data: &[u8]) -> String {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in data {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    format!("{h:016x}")
}

fn parse_kv_fields(line: &str) -> BTreeMap<String, String> {
    let mut out = BTreeMap::new();
    for p in line.split('\t').skip(2) {
        if let Some((k, v)) = p.split_once('=') {
            out.insert(k.to_string(), v.to_string());
        }
    }
    out
}

fn parse_rank_row(line: &str) -> RankRow {
    let f = parse_kv_fields(line);
    RankRow {
        rank: f["rank"].parse().unwrap(),
        hash: f["hash"].clone(),
        score: f["score"].parse().unwrap(),
        is_ref: f["isRef"] == "true",
        len: f["len"].parse().unwrap(),
    }
}

fn log_penalty(mult: u32, out: u32) -> f64 {
    if out == 0 {
        return 0.0;
    }
    (mult.max(1) as f64).log10() - (out.max(1) as f64).log10()
}

fn term_to_seq(graph: &SeqGraph, t: &SeqKbestEdgeTerm) -> String {
    String::from_utf8_lossy(&graph.vertices()[t.to].sequence).into_owned()
}

#[test]
fn forensic_6r234_c_carrier_kbest_rank_boundary() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    kv(
        "kbest_cap",
        DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH.to_string(),
    );
    assert_eq!(DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, JAVA_KBEST_K);

    let mut java128 = Vec::new();
    let mut java256 = Vec::new();
    let mut java_edges = Vec::new();
    let mut java_score_sum = None;
    for line in JAVA_TSV.lines() {
        if line.contains("\tkbest128\t") {
            java128.push(parse_rank_row(line));
        } else if line.contains("\tkbest256\t") {
            java256.push(parse_rank_row(line));
        } else if line.contains("\tcarrier_edge\t") {
            let f = parse_kv_fields(line);
            java_edges.push(JavaEdge {
                i: f["i"].parse().unwrap(),
                from: f["from"].parse().unwrap(),
                to: f["to"].parse().unwrap(),
                mult: f["mult"].parse().unwrap(),
                out: f["out"].parse().unwrap(),
                is_ref: f["is_ref"] == "true",
                penalty: f["penalty"].parse().unwrap(),
                to_seq: f["to_seq"].clone(),
            });
        } else if line.contains("\tcarrier_score_sum\t") {
            java_score_sum = line.split('\t').last().and_then(|s| s.parse().ok());
        }
    }
    assert_eq!(java128.len(), 128);
    assert_eq!(java_edges.len(), 18);
    let java_carrier = java128
        .iter()
        .find(|r| r.hash == JAVA_CARRIER)
        .expect("java K=128 carrier");
    assert_eq!(java_carrier.rank, JAVA_CARRIER_RANK);
    assert!((java_carrier.score - JAVA_CARRIER_SCORE).abs() < 1e-8);
    let java_sum: f64 = java_edges.iter().map(|e| e.penalty).sum();
    let pinned_sum: f64 = java_score_sum.expect("sum");
    assert!((java_sum - pinned_sum).abs() < 1e-9);
    assert!((java_sum - JAVA_CARRIER_SCORE).abs() < 1e-8);
    for e in &java_edges {
        let expect = log_penalty(e.mult, e.out);
        assert!(
            (e.penalty - expect).abs() < 1e-9,
            "Java carrier edge {} formula is log10(mult/out)",
            e.i
        );
    }
    kv(
        "java_carrier",
        format!(
            "K=128 rank={} score={:.8} n_edges={} formula=sum log10(mult/out) using SeqGraph getMultiplicity",
            java_carrier.rank, java_carrier.score, java_edges.len()
        ),
    );
    kv(
        "java_first_c_edge",
        format!(
            "0->12 mult={} out={} penalty={:.10} to_seq_len={}",
            java_edges[0].mult,
            java_edges[0].out,
            java_edges[0].penalty,
            java_edges[0].to_seq.len()
        ),
    );

    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing chr20_tiny BAM/REF");
        return;
    }
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, INTERVAL).expect("interval");
    let walk = traverse_assembly_region_walker(
        &dict,
        &specs,
        &ref_fasta,
        &bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
    )
    .expect("walk");
    let regions = flatten_assembly_regions(&walk);
    let region = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= TARGET
                && r.end.get() >= TARGET
        })
        .expect("ActiveFull");
    let args = CallRegionArgs::strict_java();
    let mut region_for_assemble = region.clone();
    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
    let assembled = assemble_reads_with_finalized(
        &mut region_for_assemble,
        &dict,
        &mut ref_cache,
        &args.assemble,
    )
    .expect("assemble");
    let padded_ref = assembly_reference_read(&dict, &mut ref_cache, region).expect("pad ref");
    let graph_ref = create_graph_reference_read(&padded_ref, region, &dict);
    let graph_reads = records_to_assembly_reads(&assembled.finalized_reads);
    let assembler = args.assemble.assembler.clone();
    let g = build_threading_graph_for_seq_assembly(
        &graph_ref,
        &graph_reads,
        25,
        &assembler,
        false,
        false,
    )
    .expect("seq asm")
    .expect("k=25 graph");
    let mut seq = SeqGraph::from_assembly_graph(&g);
    seq.clean_non_ref_paths();
    let _ = seq.cleanup_seq_graph();
    kv(
        "seq_graph",
        format!(
            "java_nodes=37 java_edges=52 rust_nodes={} rust_edges={}",
            seq.node_count(),
            seq.edge_count()
        ),
    );
    let off = (EVENT_POS - region.extended_start.get()) as usize;

    let paths128 = find_best_haplotypes_seq_graph(&seq, JAVA_KBEST_K).expect("k128");
    let rust128: Vec<RankRow> = paths128
        .iter()
        .enumerate()
        .map(|(rank, p)| {
            let bases = seq.path_bases_bytes(p.start, &p.edges);
            RankRow {
                rank,
                hash: fnv1a64_hex(&bases),
                score: p.score,
                is_ref: p.is_reference,
                len: bases.len(),
            }
        })
        .collect();
    assert_eq!(rust128.len(), 128);
    let k128_carrier = rust128.iter().find(|r| r.hash == JAVA_CARRIER);
    assert!(
        k128_carrier.is_some(),
        "6R.239: production K=128 includes carrier"
    );
    kv("k128_carrier_selected", "true");
    kv(
        "k128_carrier_rank",
        k128_carrier.map(|r| r.rank.to_string()).unwrap(),
    );

    let paths256 = find_best_haplotypes_seq_graph(&seq, 256).expect("k256");
    let rust256: Vec<RankRow> = paths256
        .iter()
        .enumerate()
        .map(|(rank, p)| {
            let bases = seq.path_bases_bytes(p.start, &p.edges);
            RankRow {
                rank,
                hash: fnv1a64_hex(&bases),
                score: p.score,
                is_ref: p.is_reference,
                len: bases.len(),
            }
        })
        .collect();
    let rust_carrier_row = rust256
        .iter()
        .find(|r| r.hash == JAVA_CARRIER)
        .expect("diagnostic K=256 carrier");
    assert_eq!(rust_carrier_row.rank, JAVA_CARRIER_RANK);
    assert!((rust_carrier_row.score - JAVA_CARRIER_SCORE).abs() < 1e-7);
    assert_eq!(
        paths256[rust_carrier_row.rank].score,
        rust_carrier_row.score
    );
    let cpath = &paths256[rust_carrier_row.rank];
    let cbases = seq.path_bases_bytes(cpath.start, &cpath.edges);
    assert_eq!(cbases.get(off).copied(), Some(b'C'));
    kv(
        "k256_carrier",
        format!(
            "rank={} score={:.8} len={} isRef={} n_edges={}",
            rust_carrier_row.rank,
            rust_carrier_row.score,
            rust_carrier_row.len,
            rust_carrier_row.is_ref,
            cpath.edges.len()
        ),
    );
    kv(
        "score_vs_k",
        "K=128 omitted; K=256 score computed from the same path terms (no K-dependent mutation)",
    );

    let rust_terms = seq_kbest_path_score_terms(&seq, cpath);
    let rust_sum: f64 = rust_terms.iter().map(|t| t.penalty).sum();
    assert!(
        (rust_sum - cpath.score).abs() < 1e-12,
        "Rust k-best score must equal sum of log10(mult/out) terms"
    );
    kv(
        "rust_score_sum",
        format!("{rust_sum:.10} n_edges={}", rust_terms.len()),
    );
    kv(
        "formula",
        "both sides: path.score = sum_i (log10(edgeMult_i) - log10(totalOutgoing_i)); multiplicity = SeqGraph support / Java BaseEdge.getMultiplicity(); not pruning multiplicity",
    );

    kv(
        "path_edge_n",
        format!(
            "java={} rust={} same_fnv=true same_bases_len={}",
            java_edges.len(),
            rust_terms.len(),
            rust_carrier_row.len == 460
        ),
    );
    for (i, rt) in rust_terms.iter().enumerate() {
        let to_seq = term_to_seq(&seq, rt);
        kv(
            &format!("rust_carrier_edge_{i}"),
            format!(
                "{}->{} mult={} out={} ref={} pen={:.10} to_len={} to_seq={to_seq}",
                rt.from,
                rt.to,
                rt.edge_support,
                rt.total_outgoing,
                rt.is_ref,
                rt.penalty,
                to_seq.len()
            ),
        );
    }
    for je in &java_edges {
        kv(
            &format!("java_carrier_edge_{}", je.i),
            format!(
                "{}->{} mult={} out={} ref={} pen={:.10} to_len={} to_seq={}",
                je.from,
                je.to,
                je.mult,
                je.out,
                je.is_ref,
                je.penalty,
                je.to_seq.len(),
                je.to_seq
            ),
        );
    }
    let rust_to: Vec<String> = rust_terms.iter().map(|t| term_to_seq(&seq, t)).collect();
    let java_to: Vec<String> = java_edges.iter().map(|e| e.to_seq.clone()).collect();
    let mut first_path_diff = None;
    for i in 0..java_to.len().max(rust_to.len()) {
        if java_to.get(i) != rust_to.get(i) {
            first_path_diff = Some(i);
            kv(
                "first_path_edge_mismatch",
                format!(
                    "i={i} java_to_len={:?} rust_to_len={:?} java_to={:?} rust_to={:?}",
                    java_to.get(i).map(|s| s.len()),
                    rust_to.get(i).map(|s| s.len()),
                    java_to.get(i).cloned().unwrap_or_default(),
                    rust_to.get(i).cloned().unwrap_or_default()
                ),
            );
            break;
        }
    }
    let first_i = first_path_diff;
    kv("first_path_edge_mismatch", format!("{first_i:?}"));
    assert!(
        first_i.is_none(),
        "6R.239: Java/Rust carrier k-best edge sequences match"
    );
    assert_eq!(rust_terms.len(), java_edges.len());
    for i in 0..java_edges.len() {
        assert_eq!(
            rust_to[i], java_to[i],
            "carrier edge {i} additional sequence"
        );
        assert_eq!(rust_terms[i].edge_support, java_edges[i].mult);
        assert_eq!(rust_terms[i].total_outgoing, java_edges[i].out);
        assert!(
            (rust_terms[i].penalty - java_edges[i].penalty).abs() < 1e-9,
            "penalty {i}"
        );
    }
    kv(
        "score_delta",
        format!(
            "path_d_score={:.10} (Java {java_sum:.8} vs Rust {rust_sum:.8}) after 6R.239 sentinel fix",
            rust_sum - java_sum
        ),
    );

    kv(
        "java_ranks_112_128",
        "see fixture; K=128 max rank 127; rank128 is first excluded",
    );
    for r in rust128.iter().filter(|r| (112..=127).contains(&r.rank)) {
        kv(
            &format!("rust128_rank_{}", r.rank),
            format!(
                "hash={} score={:.8} isRef={} len={}",
                r.hash, r.score, r.is_ref, r.len
            ),
        );
    }
    for r in rust256.iter().filter(|r| (112..=131).contains(&r.rank)) {
        kv(
            &format!("rust256_rank_{}", r.rank),
            format!(
                "hash={} score={:.8} isRef={} len={} java128={}",
                r.hash,
                r.score,
                r.is_ref,
                r.len,
                java128
                    .iter()
                    .find(|j| j.hash == r.hash)
                    .map(|j| j.rank.to_string())
                    .unwrap_or_else(|| "absent".into())
            ),
        );
    }

    let java128_by_hash: BTreeMap<&str, &RankRow> =
        java128.iter().map(|r| (r.hash.as_str(), r)).collect();
    let rust_top128: BTreeSet<&str> = rust128.iter().map(|r| r.hash.as_str()).collect();
    let java_top128: BTreeSet<&str> = java128.iter().map(|r| r.hash.as_str()).collect();
    let rust_only_top: Vec<&RankRow> = rust128
        .iter()
        .filter(|r| !java_top128.contains(r.hash.as_str()))
        .collect();
    let java_only_top: Vec<&RankRow> = java128
        .iter()
        .filter(|r| !rust_top128.contains(r.hash.as_str()))
        .collect();
    kv(
        "top128_symmetric_diff",
        format!(
            "rust_only={} java_only={}",
            rust_only_top.len(),
            java_only_top.len()
        ),
    );
    for r in &rust_only_top {
        kv(
            "rust_only_top128",
            format!("rank={} hash={} score={:.8}", r.rank, r.hash, r.score),
        );
    }
    for r in &java_only_top {
        kv(
            "java_only_top128",
            format!("rank={} hash={} score={:.8}", r.rank, r.hash, r.score),
        );
    }

    kv(
        "carrier_cross",
        format!(
            "hash={JAVA_CARRIER} java_rank=116 java_score={:.8} rust_rank={} rust_score={:.8}",
            JAVA_CARRIER_SCORE, rust_carrier_row.rank, rust_carrier_row.score
        ),
    );
    let outrank: Vec<&RankRow> = rust256
        .iter()
        .filter(|r| r.rank < JAVA_CARRIER_RANK)
        .filter(|r| {
            java128_by_hash
                .get(r.hash.as_str())
                .map(|j| j.rank > JAVA_CARRIER_RANK)
                .unwrap_or(true)
        })
        .collect();
    kv(
        "paths_outranking_carrier_in_rust_not_java",
        outrank.len().to_string(),
    );
    for r in &outrank {
        let jr = java128_by_hash
            .get(r.hash.as_str())
            .map(|j| format!("{}:{:.8}", j.rank, j.score))
            .unwrap_or_else(|| "absent_from_java_K128".into());
        kv(
            "rank_cross",
            format!(
                "hash={} rust_rank={} rust_score={:.8} java={}",
                r.hash, r.rank, r.score, jr
            ),
        );
    }

    let classification = "INDEX_SENTINEL_SEMANTICS_DIVERGENCE";
    kv("classification", classification);
    kv(
        "first_causal_operation",
        format!(
            "6R.239 closed the 6R.234 path-state split: Java/Rust carrier now share 18 edges, score {JAVA_CARRIER_SCORE}, rank 116. first_path_diff={first_i:?}"
        ),
    );
    kv(
        "why_k128_is_downstream",
        "K=128 was never the cause; after the sentinel fix the carrier is inside K=128. Production K is unchanged.",
    );
    kv(
        "not_investigated",
        "EventMap/trim/clip/PairHMM consequences of the now-present carrier — next arrow if FORMAT still diverges",
    );
    assert_eq!(classification, "INDEX_SENTINEL_SEMANTICS_DIVERGENCE");
    assert!((rust_carrier_row.score - java_carrier.score).abs() < 1e-7);
    assert!(
        outrank.is_empty(),
        "no Rust-only outrankers once the carrier is Java-equivalent"
    );
}
