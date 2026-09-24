//! 6R.236: first SeqGraph/RT operation that creates the Rust-only
//! 1 bp sink join (`427 → 428` ref + `971 → 428` non-ref support 1).
//! Target remains `20:29455379 G/A`. Carrier FNV `c7acc50dfb9f9ecc`
//! is identity only — this gate does not use k-best rank or VCF.
//!
//! Proof-only. Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! PRODUCTION CHANGE: NONE. Do not raise K. Do not merge the 8 bp edge.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r236_seqgraph_sink_topology -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::assembly::{
    AssemblyGraph, AssemblyGraphParams, AssemblyGraphPruningParams,
};
use gatk_haplotypecaller::assembly_based_caller::assemble_reads_with_finalized;
use gatk_haplotypecaller::assembly_dangling_recovery::DanglingRecoveryParams;
use gatk_haplotypecaller::assembly_region_finalize::{
    assembly_reference_read, create_graph_reference_read, records_to_assembly_reads,
};
use gatk_haplotypecaller::read_threading_assembler::build_threading_graph_for_seq_assembly;
use gatk_haplotypecaller::seq_graph::SeqGraph;
use gatk_haplotypecaller::{
    assembly_graph_from_ref_and_reads_threading, call_disposition, flatten_assembly_regions,
    traverse_assembly_region_walker, AssemblyRegionCallDisposition, CallRegionArgs, KmerSize,
    ReadFilterParams, WalkerTraversalConfig,
};
use rust_htslib::bam::record::Cigar as BamCigar;
use std::path::{Path, PathBuf};

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 29_455_379;
const SINK_8: &str = "TGTTTCTT";
const JAVA_TSV: &str = include_str!("forensic_6r234_java_kbest.tsv");

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R236\t{key}\t{}", value.as_ref());
}

fn vtx_seq(g: &SeqGraph, v: usize) -> String {
    String::from_utf8_lossy(&g.vertices()[v].sequence).into_owned()
}

fn in_edges_seq(g: &SeqGraph, v: usize) -> Vec<(usize, u32, bool)> {
    g.edges()
        .iter()
        .filter(|e| e.to == v)
        .map(|e| (e.from, e.support, e.is_ref))
        .collect()
}

fn out_edges_seq(g: &SeqGraph, v: usize) -> Vec<(usize, u32, bool)> {
    g.edges()
        .iter()
        .filter(|e| e.from == v)
        .map(|e| (e.to, e.support, e.is_ref))
        .collect()
}

fn rt_find_kmer(g: &AssemblyGraph, kmer: &[u8]) -> Option<usize> {
    g.nodes().iter().position(|n| n.kmer.as_ref() == kmer)
}

fn rt_edge(g: &AssemblyGraph, from_k: &[u8], to_k: &[u8]) -> Option<(usize, usize, u32, bool)> {
    let from = rt_find_kmer(g, from_k)?;
    let to = rt_find_kmer(g, to_k)?;
    let e = g
        .edges_sorted()
        .into_iter()
        .find(|e| e.from == from && e.to == to)?;
    Some((from, to, e.support, g.edge_is_ref(from, to)))
}

fn rt_outs(g: &AssemblyGraph, id: usize) -> Vec<(usize, u32, bool)> {
    g.edges_sorted()
        .into_iter()
        .filter(|e| e.from == id)
        .map(|e| (e.to, e.support, g.edge_is_ref(e.from, e.to)))
        .collect()
}

fn rt_ins(g: &AssemblyGraph, id: usize) -> Vec<(usize, u32, bool)> {
    g.edges_sorted()
        .into_iter()
        .filter(|e| e.to == id)
        .map(|e| (e.from, e.support, g.edge_is_ref(e.from, e.to)))
        .collect()
}

fn dump_rt_edge(stage: &str, g: &AssemblyGraph, from_k: &[u8], to_k: &[u8], label: &str) {
    match rt_edge(g, from_k, to_k) {
        Some((from, to, sup, is_ref)) => kv(
            &format!("{stage}_{label}"),
            format!(
                "present=true from={from} to={to} sup={sup} ref={is_ref} from_out={} to_in={} nodes={} edges={}",
                rt_outs(g, from).len(),
                rt_ins(g, to).len(),
                g.node_count(),
                g.edge_count()
            ),
        ),
        None => kv(
            &format!("{stage}_{label}"),
            format!(
                "present=false from_kmer={} to_kmer={} nodes={} edges={}",
                rt_find_kmer(g, from_k).is_some(),
                rt_find_kmer(g, to_k).is_some(),
                g.node_count(),
                g.edge_count()
            ),
        ),
    }
}

fn kmer_hex(k: &[u8]) -> String {
    String::from_utf8_lossy(k).into_owned()
}

fn threading_params(
    k: usize,
    assembler: &gatk_haplotypecaller::ReadThreadingAssemblerArgs,
) -> AssemblyGraphParams {
    AssemblyGraphParams {
        kmer_size: KmerSize::try_from_usize(k).expect("k"),
        min_base_quality: assembler.min_base_quality,
        min_edge_weight: 1,
        dangling_path_max_nodes: 0,
        max_haplotypes: assembler.num_best_haplotypes_per_graph,
        max_haplotype_bases: 4096,
        start_threading_only_at_existing_vertex: !assembler.recover_dangling_branches,
    }
}

fn cigar_str(rec: &rust_htslib::bam::Record) -> String {
    rec.cigar().to_string()
}

fn alignment_end_1based(rec: &rust_htslib::bam::Record) -> u64 {
    let mut pos = rec.pos() + 1;
    for c in rec.cigar().iter() {
        match c {
            BamCigar::Match(n)
            | BamCigar::Equal(n)
            | BamCigar::Diff(n)
            | BamCigar::Del(n)
            | BamCigar::RefSkip(n) => pos += *n as i64,
            _ => {}
        }
    }
    pos as u64 - 1
}

/// Unique-predecessor walk from the SeqGraph sink; first `in_deg!=1` is the join.
fn pre_zip_join(seq: &SeqGraph) -> (usize, Vec<(usize, u32, bool)>, String) {
    let sinks: Vec<usize> = seq
        .vertices()
        .iter()
        .filter(|v| out_edges_seq(seq, v.id).is_empty())
        .map(|v| v.id)
        .collect();
    assert_eq!(sinks.len(), 1, "unique SeqGraph sink before zip");
    let mut chain = vec![sinks[0]];
    let mut cur = sinks[0];
    let mut join_ins = Vec::new();
    let mut join_v = cur;
    for _ in 0..64 {
        let ins = in_edges_seq(seq, cur);
        if ins.len() != 1 {
            join_v = cur;
            join_ins = ins;
            break;
        }
        chain.push(ins[0].0);
        cur = ins[0].0;
        if out_edges_seq(seq, cur).len() != 1 {
            join_v = cur;
            join_ins = in_edges_seq(seq, cur);
            break;
        }
    }
    chain.reverse();
    let concat: String = chain.iter().map(|&v| vtx_seq(seq, v)).collect();
    (join_v, join_ins, concat)
}

#[test]
fn forensic_6r236_seqgraph_sink_topology() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");

    let mut java_has_gt2 = false;
    let mut java_has_8 = false;
    for line in JAVA_TSV.lines() {
        if !line.contains("\tvtx\t") {
            continue;
        }
        let mut seq = None;
        for p in line.split('\t').skip(2) {
            if let Some(("seq", v)) = p.split_once('=') {
                seq = Some(v);
            }
        }
        if let Some(s) = seq {
            if s == "GT" {
                java_has_gt2 = true;
            }
            if s == SINK_8 {
                java_has_8 = true;
            }
        }
    }
    assert!(!java_has_gt2, "Java post-cleanup SeqGraph has no GT vertex");
    assert!(!java_has_8, "Java post-cleanup SeqGraph has no 8 bp sink");
    kv(
        "java_post_cleanup_sink",
        "38bp vertex 36; no TGTTTCTT vertex; no GT vertex",
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

    let seq_asm = build_threading_graph_for_seq_assembly(
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
        "rt_seq_asm",
        format!(
            "nodes={} edges={}",
            seq_asm.node_count(),
            seq_asm.edge_count()
        ),
    );

    let seq = SeqGraph::from_assembly_graph(&seq_asm);
    let (join_v, join_ins, concat) = pre_zip_join(&seq);
    kv(
        "pre_zip_join",
        format!(
            "v={join_v} seq={} in_deg={} concat={concat} ins={join_ins:?}",
            vtx_seq(&seq, join_v),
            join_ins.len()
        ),
    );
    kv(
        "sixr239_zip_note",
        "6R.239: 971→428 splice is absent so zip is not blocked at TGTTTCTT",
    );
    let k427 = b"TGGAATTCATTTGGGTGATTTTTTT".to_vec();
    let k428 = b"GGAATTCATTTGGGTGATTTTTTTT".to_vec();
    let k970 = b"GGAATTCATTTGGGTGATTTTTTTG".to_vec();
    let k971 = b"GAATTCATTTGGGTGATTTTTTTGT".to_vec();
    kv("kmer_427", kmer_hex(&k427));
    kv("kmer_428", kmer_hex(&k428));
    kv("kmer_970", kmer_hex(&k970));
    kv("kmer_971", kmer_hex(&k971));
    assert!(
        graph_ref.bases.windows(25).any(|w| w == k428.as_slice()),
        "join vertex 428 is a reference k-mer"
    );
    assert!(
        !graph_ref.bases.windows(25).any(|w| w == k971.as_slice()),
        "971 is not a reference k-mer"
    );
    assert!(
        !graph_ref.bases.windows(25).any(|w| w == k970.as_slice()),
        "970 GT-start is not a reference k-mer"
    );

    // Vertex identity: k-mer intern. 428 exists because the reference path
    // created that k-mer; 971→428 reuses it.
    kv(
        "vertex_identity_rule",
        "AssemblyGraph unique_kmers / kmer_to_id keyed by k-mer bytes; SeqGraph ids are 1:1 with RT node index",
    );

    let params = threading_params(25, &assembler);
    let raw = assembly_graph_from_ref_and_reads_threading(&graph_ref, &graph_reads, &params)
        .expect("raw thread");
    dump_rt_edge("raw_thread", &raw, &k427, &k428, "ref_427_428");
    dump_rt_edge("raw_thread", &raw, &k971, &k428, "alt_971_428");
    dump_rt_edge("raw_thread", &raw, &k427, &k970, "gt_427_970");
    dump_rt_edge("raw_thread", &raw, &k970, &k971, "gt_970_971");
    if let Some(id971) = rt_find_kmer(&raw, &k971) {
        kv(
            "raw_971_is_sink",
            format!(
                "id={id971} out_deg={} in={:?}",
                rt_outs(&raw, id971).len(),
                rt_ins(&raw, id971)
            ),
        );
    } else {
        kv("raw_971_is_sink", "absent");
    }

    let mut pruned = raw.clone();
    let mut pruning = AssemblyGraphPruningParams::gatk_haplotype_caller_defaults();
    pruning.min_prune_factor = assembler.min_prune_factor;
    pruning.use_adaptive_pruning = assembler.use_adaptive_pruning;
    let _ = pruned.apply_pruning(&pruning);
    dump_rt_edge("after_prune", &pruned, &k427, &k428, "ref_427_428");
    dump_rt_edge("after_prune", &pruned, &k971, &k428, "alt_971_428");
    dump_rt_edge("after_prune", &pruned, &k427, &k970, "gt_427_970");
    dump_rt_edge("after_prune", &pruned, &k970, &k971, "gt_970_971");
    if let Some(id971) = rt_find_kmer(&pruned, &k971) {
        kv(
            "prune_971_is_sink",
            format!(
                "id={id971} out_deg={} in={:?}",
                rt_outs(&pruned, id971).len(),
                rt_ins(&pruned, id971)
            ),
        );
    }

    let mut dangled = pruned.clone();
    let dangling = DanglingRecoveryParams::from_assembler_args(&assembler);
    let summary = dangled
        .recover_dangling_branches(&dangling)
        .expect("dangling");
    kv(
        "dangling_summary",
        format!(
            "{summary:?} dangling_java_exact={}",
            dangling.dangling_java_exact
        ),
    );
    dump_rt_edge("after_dangling", &dangled, &k427, &k428, "ref_427_428");
    dump_rt_edge("after_dangling", &dangled, &k971, &k428, "alt_971_428");
    dump_rt_edge("after_dangling", &dangled, &k427, &k970, "gt_427_970");
    dump_rt_edge("after_dangling", &dangled, &k970, &k971, "gt_970_971");
    if let Some(id971) = rt_find_kmer(&dangled, &k971) {
        kv(
            "dangle_971_degree",
            format!(
                "id={id971} out={:?} in={:?}",
                rt_outs(&dangled, id971),
                rt_ins(&dangled, id971)
            ),
        );
    }

    let mut connected = dangled.clone();
    connected
        .remove_paths_not_connected_to_ref()
        .expect("removePaths");
    dump_rt_edge(
        "after_remove_paths",
        &connected,
        &k427,
        &k428,
        "ref_427_428",
    );
    dump_rt_edge(
        "after_remove_paths",
        &connected,
        &k971,
        &k428,
        "alt_971_428",
    );

    let seq2 = SeqGraph::from_assembly_graph(&connected);
    dump_rt_edge(
        "seqgraph_conversion",
        &connected,
        &k971,
        &k428,
        "alt_971_428",
    );
    let mut seq_clean = seq2.clone();
    seq_clean.clean_non_ref_paths();
    let (join2, ins2, concat2) = pre_zip_join(&seq_clean);
    kv(
        "after_clean_non_ref",
        format!("join={join2} in_deg={} concat={concat2}", ins2.len()),
    );
    assert_ne!(concat2, SINK_8, "6R.239: zip is not blocked at TGTTTCTT");

    let raw_has = rt_edge(&raw, &k971, &k428).is_some();
    let prune_has = rt_edge(&pruned, &k971, &k428).is_some();
    let dangle_has = rt_edge(&dangled, &k971, &k428).is_some();
    let connected_has = rt_edge(&connected, &k971, &k428).is_some();
    kv(
        "alt_edge_lifecycle",
        format!(
            "raw={raw_has} prune={prune_has} dangling={dangle_has} remove_paths={connected_has}"
        ),
    );

    let stage_class = if raw_has {
        "RAW_GRAPH_CONSTRUCTION"
    } else if prune_has {
        "EDGE_PRUNING"
    } else if dangle_has {
        "GRAPH_CLEANUP"
    } else if connected_has {
        "NON_REF_PATH_CLEANUP"
    } else {
        "LATER_GRAPH_TRANSFORMATION"
    };
    kv("earliest_stage_class", stage_class);
    assert!(!raw_has, "971→428 is absent at raw threading");
    assert!(!prune_has, "971→428 is absent after prune");
    assert!(!dangle_has, "6R.239: Java sentinel, no addEdge(971,428)");
    assert!(
        dangling.dangling_java_exact,
        "production dangling recovery is Java-exact addEdge weight=1"
    );
    assert_eq!(stage_class, "LATER_GRAPH_TRANSFORMATION");

    let first_op = if !raw_has && dangle_has {
        "recover_dangling_branches / recover_dangling_tail addEdge weight=1"
    } else if raw_has {
        "thread_sequence / extend_chain_by_one"
    } else {
        "unknown"
    };
    kv("first_creating_function", first_op);

    // Read provenance for kmer 971.
    let mut supporting = Vec::new();
    for rec in &assembled.finalized_reads {
        let seq = rec.seq().as_bytes();
        if seq.windows(25).any(|w| w == k971.as_slice()) {
            let q = String::from_utf8_lossy(rec.qname()).into_owned();
            supporting.push((
                q,
                rec.flags(),
                rec.pos() + 1,
                alignment_end_1based(rec),
                cigar_str(rec),
            ));
        }
    }
    kv("k971_supporting_finalized_n", supporting.len().to_string());
    for (q, fl, st, en, cg) in &supporting {
        kv(
            "k971_read",
            format!("qname={q} flag={fl} start={st} end={en} cigar={cg}"),
        );
    }
    let mut orig_hits = Vec::new();
    for rec in &region.reads {
        let seq = rec.seq().as_bytes();
        if seq.windows(25).any(|w| w == k971.as_slice()) {
            orig_hits.push((
                String::from_utf8_lossy(rec.qname()).into_owned(),
                rec.flags(),
                cigar_str(rec),
            ));
        }
    }
    kv("k971_supporting_original_n", orig_hits.len().to_string());
    for (q, fl, cg) in &orig_hits {
        kv("k971_orig_read", format!("qname={q} flag={fl} cigar={cg}"));
    }

    let dangle_edge = rt_edge(&dangled, &k971, &k428);
    let classification = "INDEX_SENTINEL_SEMANTICS_DIVERGENCE";
    kv("classification", classification);
    kv(
        "vertex_reuse",
        "6R.239: no support-1 edge onto 428; in_degree stays 1; zip proceeds.",
    );
    kv(
        "multiplicity_causality",
        "6R.239: splice absent; 63/78 zip boundary does not form",
    );
    kv(
        "first_divergent_operation",
        format!(
            "6R.239 closed EDGE_INSERTION: recover_dangling_tail no longer addEdge(971,428). stage={stage_class} edge={dangle_edge:?}"
        ),
    );
    kv(
        "why_k128_is_downstream",
        "K=128 unchanged. After the splice is absent the carrier is inside K=128.",
    );

    assert!(!java_has_8);
    assert!(!connected_has, "6R.239: 971→428 never enters SeqGraph");
    assert_eq!(classification, "INDEX_SENTINEL_SEMANTICS_DIVERGENCE");
    assert_eq!(stage_class, "LATER_GRAPH_TRANSFORMATION");
    assert!(!dangle_has);
}
