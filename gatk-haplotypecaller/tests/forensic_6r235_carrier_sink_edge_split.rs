//! 6R.235: first operation that splits Java's 38 bp carrier sink into
//! Rust 30 bp + 8 bp (`TGTTTCTT`). Target remains `20:29455379 G/A`.
//!
//! Proof-only. Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! PRODUCTION CHANGE: NONE. Do not raise K. Do not merge the 8 bp edge.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r235_carrier_sink_edge_split -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::assembly_based_caller::assemble_reads_with_finalized;
use gatk_haplotypecaller::assembly_region_finalize::{
    assembly_reference_read, create_graph_reference_read, records_to_assembly_reads,
};
use gatk_haplotypecaller::read_threading_assembler::build_threading_graph_for_seq_assembly;
use gatk_haplotypecaller::seq_graph::SeqGraph;
use gatk_haplotypecaller::seq_kbest_haplotype::{
    find_best_haplotypes_seq_graph, seq_kbest_path_score_terms,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    AssemblyRegionCallDisposition, CallRegionArgs, ReadFilterParams, WalkerTraversalConfig,
};
use std::collections::BTreeMap;
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
const SINK_38: &str = "GTTGGTGGAATTCATTTGGGTGATTTTTTTTGTTTCTT";
const SINK_30: &str = "GTTGGTGGAATTCATTTGGGTGATTTTTTT";
const SINK_8: &str = "TGTTTCTT";
const GTTAACT: &str = "GTTAACT";
const JAVA_TSV: &str = include_str!("forensic_6r234_java_kbest.tsv");

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R235\t{key}\t{}", value.as_ref());
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

fn vtx_seq(g: &SeqGraph, v: usize) -> String {
    String::from_utf8_lossy(&g.vertices()[v].sequence).into_owned()
}

fn in_edges(g: &SeqGraph, v: usize) -> Vec<(usize, u32, bool)> {
    g.edges()
        .iter()
        .filter(|e| e.to == v)
        .map(|e| (e.from, e.support, e.is_ref))
        .collect()
}

fn out_edges(g: &SeqGraph, v: usize) -> Vec<(usize, u32, bool)> {
    g.edges()
        .iter()
        .filter(|e| e.from == v)
        .map(|e| (e.to, e.support, e.is_ref))
        .collect()
}

fn is_ref_node(g: &SeqGraph, v: usize) -> bool {
    g.edges()
        .iter()
        .any(|e| (e.from == v || e.to == v) && e.is_ref)
}

fn ids_with_seq(g: &SeqGraph, seq: &str) -> Vec<usize> {
    g.vertices()
        .iter()
        .filter(|v| v.sequence == seq.as_bytes())
        .map(|v| v.id)
        .collect()
}

fn n_contain(g: &SeqGraph, needle: &str) -> usize {
    g.vertices()
        .iter()
        .filter(|v| String::from_utf8_lossy(&v.sequence).contains(needle))
        .count()
}

fn sink_fingerprint(g: &SeqGraph) -> String {
    let mut sinks: Vec<String> = g
        .vertices()
        .iter()
        .filter(|v| out_edges(g, v.id).is_empty())
        .map(|v| {
            format!(
                "id={} len={} seq={}",
                v.id,
                v.sequence.len(),
                String::from_utf8_lossy(&v.sequence)
            )
        })
        .collect();
    sinks.sort();
    sinks.join("|")
}

/// Walk backward from the unique out-degree-0 sink, following the unique
/// predecessor while `in_deg==1`. Stops at the first join (`in_deg!=1`).
/// Returns `(concatenated_linear_sink_seq, join_in_degree)`.
fn dump_pre_zip_sink_chain(g: &SeqGraph, label: &str) -> (String, usize) {
    let sinks: Vec<usize> = g
        .vertices()
        .iter()
        .filter(|v| out_edges(g, v.id).is_empty())
        .map(|v| v.id)
        .collect();
    kv(&format!("{label}_sink_ids"), format!("{sinks:?}"));
    let Some(&sink) = sinks.first() else {
        return (String::new(), 0);
    };
    let mut chain: Vec<usize> = vec![sink];
    let mut cur = sink;
    let mut join: Option<(usize, Vec<(usize, u32, bool)>)> = None;
    for _ in 0..64 {
        let ins = in_edges(g, cur);
        if ins.len() != 1 {
            join = Some((cur, ins));
            break;
        }
        let pred = ins[0].0;
        chain.push(pred);
        cur = pred;
        if out_edges(g, pred).len() != 1 {
            join = Some((pred, in_edges(g, pred)));
            break;
        }
    }
    chain.reverse();
    let concat: String = chain.iter().map(|&v| vtx_seq(g, v)).collect();
    let join_in = join.as_ref().map(|(_, ins)| ins.len()).unwrap_or(0);
    kv(
        &format!("{label}_linear_to_sink"),
        format!(
            "n_vertices={} concat_len={} concat={concat} join={join:?}",
            chain.len(),
            concat.len()
        ),
    );
    if let Some((join_v, ins)) = &join {
        for &(pred, sup, is_ref) in ins {
            kv(
                &format!("{label}_join_pred"),
                format!(
                    "join={join_v} from={pred} sup={sup} ref={is_ref} from_seq={} from_in={:?} from_out={:?}",
                    vtx_seq(g, pred),
                    in_edges(g, pred),
                    out_edges(g, pred)
                ),
            );
        }
    }
    for (i, &v) in chain.iter().enumerate() {
        kv(
            &format!("{label}_chain_{i}"),
            format!(
                "id={v} seq={} in={:?} out={:?}",
                vtx_seq(g, v),
                in_edges(g, v),
                out_edges(g, v)
            ),
        );
    }
    (concat, join_in)
}

fn dump_exact(g: &SeqGraph, label: &str, seq: &str) {
    for id in ids_with_seq(g, seq) {
        let ins = in_edges(g, id);
        let outs = out_edges(g, id);
        let in_sup: u32 = ins.iter().map(|e| e.1).sum();
        let out_sup: u32 = outs.iter().map(|e| e.1).sum();
        kv(
            &format!("{label}_vtx"),
            format!(
                "id={id} len={} in_deg={} out_deg={} in_sup={in_sup} out_sup={out_sup} ref={} in={ins:?} out={outs:?}",
                seq.len(),
                ins.len(),
                outs.len(),
                is_ref_node(g, id)
            ),
        );
    }
}

#[derive(Clone)]
struct StageSnap {
    stage: String,
    nodes: usize,
    edges: usize,
    n38: usize,
    n30: usize,
    n8: usize,
    contain38: usize,
    contain8: usize,
    n_gtt: usize,
    sink_fp: String,
}

fn snap_graph(stage: &str, g: &SeqGraph) -> StageSnap {
    StageSnap {
        stage: stage.to_string(),
        nodes: g.node_count(),
        edges: g.edge_count(),
        n38: ids_with_seq(g, SINK_38).len(),
        n30: ids_with_seq(g, SINK_30).len(),
        n8: ids_with_seq(g, SINK_8).len(),
        contain38: n_contain(g, SINK_38),
        contain8: n_contain(g, SINK_8),
        n_gtt: ids_with_seq(g, GTTAACT).len(),
        sink_fp: sink_fingerprint(g),
    }
}

fn emit_snap(s: &StageSnap) {
    kv(
        "seq_cleanup_snap",
        format!(
            "stage={} nodes={} edges={} exact38={} exact30={} exact8={} contain38={} contain8={} gtt={} sinks={}",
            s.stage,
            s.nodes,
            s.edges,
            s.n38,
            s.n30,
            s.n8,
            s.contain38,
            s.contain8,
            s.n_gtt,
            s.sink_fp
        ),
    );
}

#[test]
fn forensic_6r235_carrier_sink_edge_split() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    kv("carrier_fnv", JAVA_CARRIER);

    let mut java_edges = Vec::new();
    let mut java_vtx: BTreeMap<usize, (usize, bool, bool, String)> = BTreeMap::new();
    let mut java_graph_edges = Vec::new();
    for line in JAVA_TSV.lines() {
        if line.contains("\tcarrier_edge\t") {
            let f = parse_kv_fields(line);
            java_edges.push((
                f["i"].parse::<usize>().unwrap(),
                f["from"].parse::<usize>().unwrap(),
                f["to"].parse::<usize>().unwrap(),
                f["mult"].parse::<u32>().unwrap(),
                f["out"].parse::<u32>().unwrap(),
                f["to_seq"].clone(),
            ));
        } else if line.contains("\tvtx\t") {
            let f = parse_kv_fields(line);
            java_vtx.insert(
                f["id"].parse().unwrap(),
                (
                    f["len"].parse().unwrap(),
                    f["is_source"] == "true",
                    f["is_sink"] == "true",
                    f["seq"].clone(),
                ),
            );
        } else if line.contains("\tedge\t") && !line.contains("carrier_edge") {
            let f = parse_kv_fields(line);
            java_graph_edges.push((
                f["from"].parse::<usize>().unwrap(),
                f["to"].parse::<usize>().unwrap(),
                f["mult"].parse::<u32>().unwrap(),
                f["is_ref"] == "true",
            ));
        }
    }
    assert_eq!(java_edges.len(), 18);
    assert_eq!(java_vtx.len(), 37);
    assert_eq!(java_graph_edges.len(), 52);
    let java_sink = java_vtx.get(&36).expect("java vtx 36");
    assert_eq!(java_sink.0, 38);
    assert!(java_sink.2, "java vtx 36 is the unique sink");
    assert_eq!(java_sink.3, SINK_38);
    assert!(
        java_vtx.values().all(|v| v.3 != SINK_8),
        "Java has no independent 8 bp TGTTTCTT vertex"
    );
    assert!(
        java_vtx.values().all(|v| v.3 != SINK_30),
        "Java has no independent 30 bp prefix vertex"
    );
    let java_into_36: Vec<_> = java_graph_edges
        .iter()
        .filter(|e| e.1 == 36)
        .copied()
        .collect();
    assert_eq!(java_into_36.len(), 4);
    assert!(java_into_36.iter().all(|e| e.2 == 1));
    let java_out_36 = java_graph_edges.iter().filter(|e| e.0 == 36).count();
    assert_eq!(java_out_36, 0);
    kv(
        "java_sink",
        format!(
            "vtx36 len=38 is_sink=true in_deg=4 out_deg=0 incoming={java_into_36:?} seq={SINK_38}"
        ),
    );
    kv(
        "java_carrier_edge_17",
        format!("31->36 mult=1 out=1 to_len=38 seq={}", java_edges[17].5),
    );

    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        panic!("missing chr20_tiny BAM/REF");
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
    kv(
        "rt_k25",
        format!("nodes={} edges={}", g.node_count(), g.edge_count()),
    );

    let mut seq = SeqGraph::from_assembly_graph(&g);
    let mut snaps: Vec<StageSnap> = Vec::new();
    snaps.push(snap_graph("from_assembly_graph", &seq));
    let (pre_zip_concat, pre_zip_join_in) = dump_pre_zip_sink_chain(&seq, "pre_zip");
    seq.clean_non_ref_paths();
    snaps.push(snap_graph("after_clean_non_ref_paths", &seq));
    let (after_nonref_concat, after_nonref_join_in) =
        dump_pre_zip_sink_chain(&seq, "after_nonref_pre_zip");
    let _ = seq.traced_cleanup_seq_graph(|stage, g| {
        snaps.push(snap_graph(stage, g));
    });
    for s in &snaps {
        emit_snap(s);
    }

    // Re-walk with live dumps of 8/30/38 objects at the first stage they appear.
    let mut seq2 = SeqGraph::from_assembly_graph(&g);
    kv(
        "pre_cleanup",
        format!(
            "nodes={} edges={} exact38={} exact30={} exact8={}",
            seq2.node_count(),
            seq2.edge_count(),
            ids_with_seq(&seq2, SINK_38).len(),
            ids_with_seq(&seq2, SINK_30).len(),
            ids_with_seq(&seq2, SINK_8).len()
        ),
    );
    seq2.clean_non_ref_paths();
    let mut first_8: Option<String> = None;
    let mut first_30: Option<String> = None;
    let mut first_38: Option<String> = None;
    let mut before_first_8: Option<StageSnap> = None;
    let mut prev = snap_graph("after_clean_non_ref_paths", &seq2);
    if prev.n8 > 0 && first_8.is_none() {
        first_8 = Some(prev.stage.clone());
    }
    let _ = seq2.traced_cleanup_seq_graph(|stage, g| {
        let s = snap_graph(stage, g);
        if s.n8 > 0 && first_8.is_none() {
            first_8 = Some(s.stage.clone());
            before_first_8 = Some(prev.clone());
            dump_exact(g, "first_8", SINK_8);
            dump_exact(g, "first_8_30", SINK_30);
            dump_exact(g, "first_8_38", SINK_38);
            dump_exact(g, "first_8_gtt", GTTAACT);
        }
        if s.n30 > 0 && first_30.is_none() {
            first_30 = Some(s.stage.clone());
            dump_exact(g, "first_30", SINK_30);
        }
        if s.n38 > 0 && first_38.is_none() {
            first_38 = Some(s.stage.clone());
        }
        prev = s;
    });

    kv(
        "first_exact_8bp_stage",
        first_8.clone().unwrap_or_else(|| "never".into()),
    );
    kv(
        "first_exact_30bp_stage",
        first_30.clone().unwrap_or_else(|| "never".into()),
    );
    kv(
        "first_exact_38bp_stage",
        first_38.clone().unwrap_or_else(|| "never".into()),
    );
    if let Some(b) = &before_first_8 {
        kv(
            "stage_before_first_8",
            format!(
                "stage={} nodes={} edges={} exact38={} exact30={} exact8={} contain38={} contain8={}",
                b.stage, b.nodes, b.edges, b.n38, b.n30, b.n8, b.contain38, b.contain8
            ),
        );
    }

    kv(
        "seq_post_cleanup",
        format!(
            "java_nodes=37 java_edges=52 rust_nodes={} rust_edges={}",
            seq.node_count(),
            seq.edge_count()
        ),
    );
    dump_exact(&seq, "final_38", SINK_38);
    dump_exact(&seq, "final_30", SINK_30);
    dump_exact(&seq, "final_8", SINK_8);
    dump_exact(&seq, "final_gtt", GTTAACT);
    if let Some(v30) = ids_with_seq(&seq, SINK_30).first().copied() {
        for (to, sup, is_ref) in out_edges(&seq, v30) {
            kv(
                "v30_out_target",
                format!(
                    "to={to} sup={sup} ref={is_ref} to_seq={} to_in={:?} to_out={:?}",
                    vtx_seq(&seq, to),
                    in_edges(&seq, to),
                    out_edges(&seq, to)
                ),
            );
        }
    }

    let rust_sinks: Vec<_> = seq
        .vertices()
        .iter()
        .filter(|v| out_edges(&seq, v.id).is_empty())
        .map(|v| (v.id, v.sequence.len(), vtx_seq(&seq, v.id)))
        .collect();
    kv("rust_sink_vertices", format!("{rust_sinks:?}"));

    let off = (EVENT_POS - region.extended_start.get()) as usize;
    let paths128 = find_best_haplotypes_seq_graph(&seq, 128).expect("k128");
    let k128_idx = paths128
        .iter()
        .position(|p| fnv1a64_hex(&seq.path_bases_bytes(p.start, &p.edges)) == JAVA_CARRIER)
        .expect("6R.239: K=128 includes carrier");
    assert_eq!(k128_idx, JAVA_CARRIER_RANK);
    let paths256 = find_best_haplotypes_seq_graph(&seq, 256).expect("k256");
    let rust_idx = paths256
        .iter()
        .position(|p| fnv1a64_hex(&seq.path_bases_bytes(p.start, &p.edges)) == JAVA_CARRIER)
        .expect("K=256 carrier");
    assert_eq!(rust_idx, JAVA_CARRIER_RANK);
    let cpath = &paths256[rust_idx];
    assert!((cpath.score - JAVA_CARRIER_SCORE).abs() < 1e-7);
    let cbases = seq.path_bases_bytes(cpath.start, &cpath.edges);
    assert_eq!(cbases.get(off).copied(), Some(b'C'));
    let terms = seq_kbest_path_score_terms(&seq, cpath);
    assert_eq!(terms.len(), 18);

    let mut cum_bases = seq.vertices()[cpath.start].sequence.clone();
    let mut cum_score = 0.0f64;
    kv(
        "path_start",
        format!(
            "v={} len={} seq={}",
            cpath.start,
            cum_bases.len(),
            String::from_utf8_lossy(&cum_bases)
        ),
    );
    for (i, t) in terms.iter().enumerate() {
        let to_seq = vtx_seq(&seq, t.to);
        cum_bases.extend_from_slice(seq.vertices()[t.to].sequence.as_slice());
        cum_score += t.penalty;
        let ins = in_edges(&seq, t.to);
        let outs = out_edges(&seq, t.to);
        kv(
            &format!("carrier_step_{i}"),
            format!(
                "java={}->{} rust={}->{} java_len={} rust_len={} mult={} prune_n/a out={} ref={} pen={:.10} cum_len={} cum_score={:.10} to_in={} to_out={} to_seq={to_seq}",
                java_edges.get(i).map(|e| e.1).unwrap_or(usize::MAX),
                java_edges.get(i).map(|e| e.2).unwrap_or(usize::MAX),
                t.from,
                t.to,
                java_edges.get(i).map(|e| e.5.len()).unwrap_or(0),
                to_seq.len(),
                t.edge_support,
                t.total_outgoing,
                t.is_ref,
                t.penalty,
                cum_bases.len(),
                cum_score,
                ins.len(),
                outs.len()
            ),
        );
    }
    assert_eq!(fnv1a64_hex(&cum_bases), JAVA_CARRIER);

    for i in 0..18 {
        assert_eq!(vtx_seq(&seq, terms[i].to), java_edges[i].5);
        assert_eq!(terms[i].edge_support, java_edges[i].3);
        assert_eq!(terms[i].total_outgoing, java_edges[i].4);
    }
    assert_eq!(vtx_seq(&seq, terms[17].to), SINK_38);
    assert_eq!(java_edges[17].5, SINK_38);
    assert!((cpath.score - JAVA_CARRIER_SCORE).abs() < 1e-7);

    let classification = "INDEX_SENTINEL_SEMANTICS_DIVERGENCE";
    kv("classification", classification);
    kv(
        "first_divergent_operation",
        "6R.239: Java-compatible refIndexToMerge sentinel; 38 bp sink zips; carrier 18 edges rank 116",
    );
    kv(
        "why_k128_is_downstream",
        "K=128 was never changed. After the splice is absent the carrier is inside K=128.",
    );
    kv("java_rank", JAVA_CARRIER_RANK.to_string());
    kv("rust_diag_rank", JAVA_CARRIER_RANK.to_string());

    assert!(
        ids_with_seq(&seq, SINK_38).len() == 1,
        "6R.239: Java 38 bp sink exists"
    );
    assert_eq!(classification, "INDEX_SENTINEL_SEMANTICS_DIVERGENCE");
    assert_eq!(DEFAULT_KBEST_UNCHANGED, 128);
}

const DEFAULT_KBEST_UNCHANGED: usize =
    gatk_haplotypecaller::read_threading_assembler::DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH;
