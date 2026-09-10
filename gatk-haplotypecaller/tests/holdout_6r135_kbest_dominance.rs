//! 6R.135 forensic: pop-1494 expand-refuse dominance / heap-state analysis.
//! Skipped unless `HOLDOUT_6R135=1`. Production change: NONE.
//!
//! ```text
//! HOLDOUT_6R135=1 cargo test -p gatk-haplotypecaller --test holdout_6r135_kbest_dominance -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::assembly_region_finalize::{
    create_graph_reference_read, records_to_assembly_reads,
};
use gatk_haplotypecaller::kbest_haplotype::KBestPath;
use gatk_haplotypecaller::read_threading_assembler::build_threading_graph_for_seq_assembly;
use gatk_haplotypecaller::seq_graph::{SeqGraph, SeqGraphCleanupStatus};
use gatk_haplotypecaller::seq_kbest_haplotype::{
    find_best_haplotypes_seq_graph_forensic, SeqKbestCapPolicy, SeqKbestExpandRefuse,
    SeqKbestFirstRefuseSnapshot, SEQ_KBEST_PRODUCTION_MAX_EXPANSIONS,
    SEQ_KBEST_PRODUCTION_MAX_HEAP, SEQ_KBEST_PRODUCTION_MAX_PATH_EDGES,
};
use gatk_haplotypecaller::{
    assemble_reads_with_finalized, call_disposition, flatten_assembly_regions,
    traverse_assembly_region_walker, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    ReadFilterParams, WalkerTraversalConfig,
};
use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_KBEST_DUMP_REL: &str = "parity/giab/runs/local-pairhmm-diff/6r133_java.txt";
const TARGET: u64 = 29_456_196;
const JAVA_ACTIVE_START: u64 = 29_455_995;
const JAVA_ACTIVE_END: u64 = 29_456_293;
const JAVA_PAD_START: u64 = 29_455_895;
const JAVA_PAD_END: u64 = 29_456_393;
const EXPECTED_REF_HASH: &str = "659741f99b7f78a5";
const EXPECTED_FINALIZED: usize = 379;
const JAVA_RANK56: &str = "28b8f48b3a0c4954";

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R135\t{key}\t{}", value.as_ref());
}

fn fnv1a64(bases: &[u8]) -> u64 {
    let mut h = 0xcbf29ce484222325u64;
    for &b in bases {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

fn hex(h: u64) -> String {
    format!("{h:016x}")
}

fn path_hash(graph: &SeqGraph, start: usize, edges: &[(usize, usize)]) -> String {
    hex(fnv1a64(&graph.path_bases_bytes(start, edges)))
}

fn unique_path_hashes(graph: &SeqGraph, paths: &[KBestPath]) -> BTreeSet<String> {
    paths
        .iter()
        .map(|p| path_hash(graph, p.start, &p.edges))
        .collect()
}

fn edges_have_prefix(
    full: &[(usize, usize)],
    prefix: &[(usize, usize)],
    extra: (usize, usize),
) -> bool {
    full.len() >= prefix.len() + 1
        && full[..prefix.len()] == prefix[..]
        && full[prefix.len()] == extra
}

fn load_java_k128(path: &Path) -> (BTreeSet<String>, BTreeMap<String, (usize, f64)>) {
    let mut set = BTreeSet::new();
    let mut ranks = BTreeMap::new();
    let Ok(text) = std::fs::read_to_string(path) else {
        return (set, ranks);
    };
    for line in text.lines() {
        let mut parts = line.split('\t');
        if parts.next() != Some("6R133") {
            continue;
        }
        if parts.next() != Some("kbest") {
            continue;
        }
        let mut kmer = None;
        let mut kcap = None;
        let mut rank = None;
        let mut hash = None;
        let mut score = None;
        for tok in parts {
            if let Some(v) = tok.strip_prefix("kmer=") {
                kmer = v.parse().ok();
            }
            if let Some(v) = tok.strip_prefix("K=") {
                kcap = v.parse().ok();
            }
            if let Some(v) = tok.strip_prefix("rank=") {
                rank = v.parse().ok();
            }
            if let Some(v) = tok.strip_prefix("hash=") {
                hash = Some(v.to_string());
            }
            if let Some(v) = tok.strip_prefix("score=") {
                score = v.parse().ok();
            }
        }
        if kmer == Some(25) && kcap == Some(128) {
            if let (Some(h), Some(r), Some(s)) = (hash, rank, score) {
                set.insert(h.clone());
                ranks.insert(h, (r, s));
            }
        }
    }
    (set, ranks)
}

fn outs_of(graph: &SeqGraph, from: usize) -> Vec<(usize, u32)> {
    graph
        .edges()
        .iter()
        .filter(|e| e.from == from)
        .map(|e| (e.to, e.support))
        .collect()
}

fn can_reach_sink(graph: &SeqGraph, start: usize, sink: usize) -> bool {
    if start == sink {
        return true;
    }
    let mut seen = vec![false; graph.node_count()];
    let mut q = VecDeque::new();
    seen[start] = true;
    q.push_back(start);
    while let Some(v) = q.pop_front() {
        for (to, _) in outs_of(graph, v) {
            if !seen[to] {
                if to == sink {
                    return true;
                }
                seen[to] = true;
                q.push_back(to);
            }
        }
    }
    false
}

fn per_dest_policy(per_dest: usize) -> SeqKbestCapPolicy {
    SeqKbestCapPolicy {
        max_heap_paths: Some(SEQ_KBEST_PRODUCTION_MAX_HEAP),
        max_expansions: Some(SEQ_KBEST_PRODUCTION_MAX_EXPANSIONS),
        max_path_edges: Some(SEQ_KBEST_PRODUCTION_MAX_PATH_EDGES),
        max_heap_paths_per_dest: Some(per_dest),
    }
}

fn build_seqgraph() -> SeqGraph {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
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
    let covering = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= TARGET
                && r.end.get() >= TARGET
        })
        .expect("ActiveFull")
        .clone();
    let mut java_bounds = covering;
    java_bounds.start = GenomePosition::new_1based(JAVA_ACTIVE_START);
    java_bounds.end = GenomePosition::new_1based(JAVA_ACTIVE_END);
    java_bounds.extended_start = GenomePosition::new_1based(JAVA_PAD_START);
    java_bounds.extended_end = GenomePosition::new_1based(JAVA_PAD_END);
    let args = CallRegionArgs::strict_java();
    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
    let mut owned = java_bounds.clone();
    let assembled =
        assemble_reads_with_finalized(&mut owned, &dict, &mut ref_cache, &args.assemble)
            .expect("assemble");
    let graph_ref_hash = hex(fnv1a64(assembled.assembly.reference_bases()));
    kv(
        "matched_input",
        format!(
            "graph_ref={graph_ref_hash}\tlen={}\tfinalized={}",
            assembled.assembly.reference_bases().len(),
            assembled.finalized_reads.len()
        ),
    );
    assert_eq!(graph_ref_hash, EXPECTED_REF_HASH);
    assert_eq!(assembled.finalized_reads.len(), EXPECTED_FINALIZED);
    let padded_ref = gatk_haplotypecaller::assembly_region_finalize::assembly_reference_read(
        &dict,
        &mut ref_cache,
        &java_bounds,
    )
    .expect("pad ref");
    let graph_ref = create_graph_reference_read(&padded_ref, &java_bounds, &dict);
    let graph_reads = records_to_assembly_reads(&assembled.finalized_reads);
    let mut assembler = args.assemble.assembler.clone();
    assembler.dangling_java_exact = true;
    let rt = build_threading_graph_for_seq_assembly(
        &graph_ref,
        &graph_reads,
        25,
        &assembler,
        false,
        false,
    )
    .expect("rt")
    .expect("k=25 graph");
    let mut seq = SeqGraph::from_assembly_graph(&rt);
    seq.clean_non_ref_paths();
    assert_eq!(
        seq.cleanup_seq_graph(),
        SeqGraphCleanupStatus::AssembledSomeVariation
    );
    kv(
        "seqgraph",
        format!("nodes={}\tedges={}", seq.node_count(), seq.edge_count()),
    );
    seq
}

fn successor_java_hits(
    graph: &SeqGraph,
    rec: &SeqKbestExpandRefuse,
    to: usize,
    unbounded: &[KBestPath],
    java: &BTreeSet<String>,
    ranks: &BTreeMap<String, (usize, f64)>,
) -> Vec<(String, usize, f64)> {
    let extra = (rec.parent_last, to);
    let mut hits = Vec::new();
    for p in unbounded {
        if !edges_have_prefix(&p.edges, &rec.parent_edges, extra) {
            continue;
        }
        let h = path_hash(graph, p.start, &p.edges);
        if java.contains(&h) {
            let (r, s) = ranks.get(&h).copied().unwrap_or((usize::MAX, p.score));
            hits.push((h, r, s));
        }
    }
    hits
}

fn dump_heap_by_dest(snap: &SeqKbestFirstRefuseSnapshot, k: usize) {
    let mut by: BTreeMap<usize, Vec<(f64, usize, bool, u64, u64)>> = BTreeMap::new();
    for e in &snap.heap {
        by.entry(e.last).or_default().push((
            e.score,
            e.n_edges,
            e.is_sink,
            e.edges_fnv,
            e.bases_fnv,
        ));
    }
    let mut n_gt1 = 0usize;
    let mut n_gt8 = 0usize;
    let mut n_gt32 = 0usize;
    let mut n_gt_k = 0usize;
    let mut max_n = 0usize;
    let mut sink_n = 0usize;
    for (dest, rows) in &by {
        let n = rows.len();
        max_n = max_n.max(n);
        if n > 1 {
            n_gt1 += 1;
        }
        if n > 8 {
            n_gt8 += 1;
        }
        if n > 32 {
            n_gt32 += 1;
        }
        if n > k {
            n_gt_k += 1;
        }
        let n_sink = rows.iter().filter(|r| r.2).count();
        sink_n += n_sink;
        let scores: Vec<f64> = rows.iter().map(|r| r.0).collect();
        let min_s = scores.iter().copied().fold(f64::INFINITY, f64::min);
        let max_s = scores.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let uniq_edges: BTreeSet<u64> = rows.iter().map(|r| r.3).collect();
        let uniq_bases: BTreeSet<u64> = rows.iter().map(|r| r.4).collect();
        kv(
            "heap_dest",
            format!(
                "last={dest}\tn={n}\tn_sink={n_sink}\tmin_score={min_s:.9}\tmax_score={max_s:.9}\tuniq_edges={}\tuniq_bases={}\tn_edges_min={}\tn_edges_max={}",
                uniq_edges.len(),
                uniq_bases.len(),
                rows.iter().map(|r| r.1).min().unwrap_or(0),
                rows.iter().map(|r| r.1).max().unwrap_or(0)
            ),
        );
    }
    let uniq_edges: BTreeSet<u64> = snap.heap.iter().map(|e| e.edges_fnv).collect();
    let uniq_bases: BTreeSet<u64> = snap.heap.iter().map(|e| e.bases_fnv).collect();
    let uniq_p8: BTreeSet<u64> = snap.heap.iter().map(|e| e.prefix8_fnv).collect();
    let uniq_p16: BTreeSet<u64> = snap.heap.iter().map(|e| e.prefix16_fnv).collect();
    let total_edges: usize = snap.heap.iter().map(|e| e.n_edges).sum();
    kv(
        "heap_summary",
        format!(
            "n={}\tdestinct_dest={}\tmax_per_dest={max_n}\tdest_gt1={n_gt1}\tdest_gt8={n_gt8}\tdest_gt32={n_gt32}\tdest_gt_k={n_gt_k}\tsink_entries={sink_n}\tuniq_full_edges={}\tuniq_bases={}\tuniq_prefix8={}\tuniq_prefix16={}\ttotal_owned_edge_pairs={}\tedge_bytes={}",
            snap.heap.len(),
            by.len(),
            uniq_edges.len(),
            uniq_bases.len(),
            uniq_p8.len(),
            uniq_p16.len(),
            total_edges,
            total_edges * std::mem::size_of::<(usize, usize)>()
        ),
    );
}

#[test]
fn holdout_6r135_kbest_dominance() {
    if std::env::var("HOLDOUT_6R135").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R135=1");
        return;
    }
    kv(
        "sizes",
        format!(
            "KBestPath={}\tpair_usize={}\theap_cap={SEQ_KBEST_PRODUCTION_MAX_HEAP}",
            std::mem::size_of::<KBestPath>(),
            std::mem::size_of::<(usize, usize)>()
        ),
    );
    let root = repo_root();
    let (java, ranks) = load_java_k128(&root.join(JAVA_KBEST_DUMP_REL));
    kv("java_k128_loaded", java.len().to_string());
    let seq = build_seqgraph();
    let sink = seq.reference_sink_vertex().expect("sink");

    let unbounded = find_best_haplotypes_seq_graph_forensic(
        &seq,
        128,
        128,
        SeqKbestCapPolicy::unbounded(),
        &[],
    )
    .expect("unbounded");
    let unb_set = unique_path_hashes(&seq, &unbounded.paths);
    kv(
        "unbounded",
        format!(
            "n={}\tmax_heap={}\tCOMMON={}",
            unbounded.paths.len(),
            unbounded.max_heap,
            unb_set.intersection(&java).count()
        ),
    );

    let prod = find_best_haplotypes_seq_graph_forensic(
        &seq,
        128,
        128,
        SeqKbestCapPolicy::production(),
        &[],
    )
    .expect("prod");
    let prod_set = unique_path_hashes(&seq, &prod.paths);
    let java_only: BTreeSet<_> = java.difference(&prod_set).cloned().collect();
    kv(
        "prod",
        format!(
            "n={}\tmax_heap={}\tskip_exp={}\tCOMMON={}\tJAVA_ONLY={}\tRUST_ONLY={}",
            prod.paths.len(),
            prod.max_heap,
            prod.skip_heap_full_at_expand,
            prod_set.intersection(&java).count(),
            java_only.len(),
            prod_set.difference(&java).count()
        ),
    );

    let rec = prod.first_expand_refuse.as_ref().expect("first refuse");
    let snap = prod.first_refuse_snapshot.as_ref().expect("snapshot");
    let parent_h = path_hash(&seq, rec.parent_start, &rec.parent_edges);
    let parent_bases = seq.path_bases_bytes(rec.parent_start, &rec.parent_edges);
    kv(
        "pop1494_parent",
        format!(
            "pop={}\tresult_len={}\theap_len={}\tscore={:.9}\tlast={}\tstart={}\tn_edges={}\thash={parent_h}\tbases_len={}\tvisit_parent={}\tn_outs={}\tinserted={}\trefused={}",
            rec.pop_count,
            rec.result_len,
            rec.heap_len_at_refuse,
            rec.parent_score,
            rec.parent_last,
            rec.parent_start,
            rec.parent_n_edges,
            parent_bases.len(),
            snap.vertex_counts.get(rec.parent_last).copied().unwrap_or(0),
            rec.n_outs,
            rec.inserted_tos.len(),
            rec.refused_tos.len()
        ),
    );

    for (kind, tos) in [
        ("accepted", rec.inserted_tos.as_slice()),
        ("refused", rec.refused_tos.as_slice()),
    ] {
        for &(to, score, sup) in tos {
            let mut edges = rec.parent_edges.clone();
            edges.push((rec.parent_last, to));
            let sh = path_hash(&seq, rec.parent_start, &edges);
            let blen = seq.path_bases_bytes(rec.parent_start, &edges).len();
            let n_heap_at = snap.heap.iter().filter(|e| e.last == to).count();
            let visit = snap.vertex_counts.get(to).copied().unwrap_or(0);
            let reach = can_reach_sink(&seq, to, sink);
            let out_n = outs_of(&seq, to).len();
            let hits = successor_java_hits(&seq, rec, to, &unbounded.paths, &java, &ranks);
            let in_prod = successor_java_hits(&seq, rec, to, &prod.paths, &java, &ranks);
            kv(
                "pop1494_succ",
                format!(
                    "kind={kind}\tto={to}\tsupport={sup}\tscore={score:.9}\thash={sh}\tbases_len={blen}\treach_sink={reach}\tout_deg={out_n}\theap_at_dest={n_heap_at}\tvisit_dest={visit}\tjava_top128={}\tprod_java_hits={}\tjava_ranks={}",
                    hits.len(),
                    in_prod.len(),
                    hits.iter()
                        .map(|(h, r, s)| format!("{h}:rank{r}:{s:.8}"))
                        .collect::<Vec<_>>()
                        .join(",")
                ),
            );
        }
    }

    kv(
        "siblings_same_dest",
        format!(
            "{}",
            rec.inserted_tos.first().map(|t| t.0) == rec.refused_tos.first().map(|t| t.0)
        ),
    );

    dump_heap_by_dest(snap, 128);

    let refused_to = rec.refused_tos[0].0;
    let refused_score = rec.refused_tos[0].1;
    let at_refused: Vec<f64> = snap
        .heap
        .iter()
        .filter(|e| e.last == refused_to)
        .map(|e| e.score)
        .collect();
    let n_better_or_eq = at_refused.iter().filter(|&&s| s >= refused_score).count();
    kv(
        "refused_dest_dominance",
        format!(
            "to={refused_to}\theap_already_at_dest={}\tn_score_ge_refused={n_better_or_eq}\trefused_score={refused_score:.9}\tvisit={}",
            at_refused.len(),
            snap.vertex_counts.get(refused_to).copied().unwrap_or(0)
        ),
    );

    let accepted_to = rec.inserted_tos[0].0;
    let accepted_outs: BTreeSet<usize> = outs_of(&seq, accepted_to)
        .into_iter()
        .map(|(t, _)| t)
        .collect();
    let refused_outs: BTreeSet<usize> = outs_of(&seq, refused_to)
        .into_iter()
        .map(|(t, _)| t)
        .collect();
    kv(
        "future_outs",
        format!(
            "accepted_to={accepted_to}\trefused_to={refused_to}\tacc_outs={:?}\tref_outs={:?}\tsame_out_set={}",
            accepted_outs, refused_outs, accepted_outs == refused_outs
        ),
    );

    let mut n_java_only_refuses = 0usize;
    let mut emitted = 0usize;
    for (ri, r) in prod.expand_refuses.iter().enumerate() {
        for &(to, score, _) in &r.refused_tos {
            let hits = successor_java_hits(&seq, r, to, &unbounded.paths, &java, &ranks);
            let n_only = hits
                .iter()
                .filter(|(h, _, _)| java_only.contains(h))
                .count();
            if n_only == 0 {
                continue;
            }
            n_java_only_refuses += 1;
            if emitted < 8 {
                kv(
                    "later_java_only_refuse",
                    format!(
                        "idx={ri}\tpop={}\tresult_len={}\tparent_score={:.9}\tto={to}\tsucc_score={score:.9}\tparent_last={}\theap={}\tjava_only_hits={n_only}\tranks={}",
                        r.pop_count,
                        r.result_len,
                        r.parent_score,
                        r.parent_last,
                        r.heap_len_at_refuse,
                        hits.iter()
                            .filter(|(h, _, _)| java_only.contains(h))
                            .map(|(_, rank, _)| rank.to_string())
                            .collect::<Vec<_>>()
                            .join(",")
                    ),
                );
                emitted += 1;
            }
        }
    }
    kv("later_java_only_refuse_n", n_java_only_refuses.to_string());

    for per in [128usize, 1] {
        let r = find_best_haplotypes_seq_graph_forensic(&seq, 128, 128, per_dest_policy(per), &[])
            .expect("per-dest");
        let set = unique_path_hashes(&seq, &r.paths);
        kv(
            "counterfactual_per_dest",
            format!(
                "per_dest={per}\tmax_heap={}\trefused_exp={}\tn={}\tCOMMON={}\tJAVA_ONLY={}\tRUST_ONLY={}\thas_rank56={}",
                r.max_heap,
                r.skip_heap_full_at_expand,
                r.paths.len(),
                set.intersection(&java).count(),
                java.difference(&set).count(),
                set.difference(&java).count(),
                set.contains(JAVA_RANK56)
            ),
        );
    }
}
