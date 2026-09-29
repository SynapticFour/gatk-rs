//! 6R.237: first dangling-tail decision that inserts Rust-only
//! `971 → 428` (weight 1). Target `20:29455379 G/A`. Carrier FNV
//! `c7acc50dfb9f9ecc` is identity only — this gate does not use k-best
//! rank or VCF.
//!
//! Proof-only. Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! PRODUCTION CHANGE: NONE. Do not raise K. Do not delete the splice.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r237_dangling_tail_splice_decision -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::assembly::{AssemblyGraph, AssemblyGraphPruningParams};
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
use std::path::{Path, PathBuf};

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 29_455_379;
const K971: &[u8] = b"GAATTCATTTGGGTGATTTTTTTGT";
const K428: &[u8] = b"GGAATTCATTTGGGTGATTTTTTTT";
const JAVA_TSV: &str = include_str!("forensic_6r234_java_kbest.tsv");

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R237\t{key}\t{}", value.as_ref());
}

fn kmer_str(k: &[u8]) -> String {
    String::from_utf8_lossy(k).into_owned()
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

fn threading_params(
    k: usize,
    assembler: &gatk_haplotypecaller::ReadThreadingAssemblerArgs,
) -> gatk_haplotypecaller::assembly::AssemblyGraphParams {
    gatk_haplotypecaller::assembly::AssemblyGraphParams {
        kmer_size: KmerSize::try_from_usize(k).expect("k"),
        min_base_quality: assembler.min_base_quality,
        min_edge_weight: 1,
        dangling_path_max_nodes: 0,
        max_haplotypes: assembler.num_best_haplotypes_per_graph,
        max_haplotype_bases: 4096,
        start_threading_only_at_existing_vertex: !assembler.recover_dangling_branches,
    }
}

fn dump_one(
    tag: &str,
    d: &gatk_haplotypecaller::assembly_dangling_recovery::DanglingTailDecisionDump,
) {
    kv(
        tag,
        format!(
            "sink={} kmer={} classified={} reject={:?} alt_len={} ref_len={} cigar={} cigar_ok={} match_suf={} alt_idx={:?} ref_idx={:?} from_k={} to_k={} mm={} indel={} offset={} java_pred={}",
            d.sink,
            kmer_str(&d.sink_kmer),
            d.classified_dangling_tail,
            d.reject_reason,
            d.alt_path_ids.len(),
            d.ref_path_ids.len(),
            d.cigar,
            d.cigar_ok,
            d.matching_suffix,
            d.alt_index_to_merge,
            d.ref_index_to_merge,
            kmer_str(&d.from_kmer),
            kmer_str(&d.to_kmer),
            d.mismatch_count,
            d.indel_count,
            d.alignment_offset,
            d.java_merge_predicate
        ),
    );
    if !d.alt_path_ids.is_empty() {
        kv(&format!("{tag}_alt_path"), format!("{:?}", d.alt_path_ids));
        kv(&format!("{tag}_alt_bases"), kmer_str(&d.alt_bases));
    }
    if !d.ref_path_ids.is_empty() {
        kv(
            &format!("{tag}_ref_path_head"),
            format!(
                "ids_head={:?} n={} bases_head={}",
                &d.ref_path_ids[..d.ref_path_ids.len().min(8)],
                d.ref_path_ids.len(),
                kmer_str(&d.ref_bases[..d.ref_bases.len().min(40)])
            ),
        );
    }
}

#[test]
fn forensic_6r237_dangling_tail_splice_decision() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");

    let mut java_has_8 = false;
    for line in JAVA_TSV.lines() {
        if line.contains("seq=") && line.contains("TGTTTCTT") && line.contains("len=8") {
            java_has_8 = true;
        }
    }
    kv(
        "java_post_cleanup_sink",
        "38bp vertex 36; no TGTTTCTT vertex (6R.235/6R.236)",
    );
    kv(
        "java_rt_after_prune",
        "k=25 nodes=1064 edges=1079 (forensic_6r233_java_seqgraph.tsv; same as Rust prune)",
    );
    assert!(!java_has_8);

    let root = repo_root();
    let bam = root.join(BAM_REL);
    let ref_fasta = root.join(REF_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing chr20_tiny BAM/REF");
        return;
    }
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, INTERVAL).expect("interval");
    let walked = traverse_assembly_region_walker(
        &dict,
        &specs,
        &ref_fasta,
        &bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
    )
    .expect("walk");
    let regions = flatten_assembly_regions(&walked);
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
    kv(
        "assembler_knobs",
        format!(
            "min_dangling={} recover_all={} java_exact={} min_prune={} min_match={}",
            assembler.min_dangling_branch_length,
            assembler.recover_all_dangling_branches,
            assembler.dangling_java_exact,
            assembler.min_prune_factor,
            assembler.min_matching_bases_to_dangling_end_recovery
        ),
    );
    assert!(assembler.dangling_java_exact);
    assert_eq!(assembler.min_dangling_branch_length, 4);
    assert!(!assembler.recover_all_dangling_branches);
    assert_eq!(assembler.min_matching_bases_to_dangling_end_recovery, -1);

    let params = threading_params(25, &assembler);
    let raw = assembly_graph_from_ref_and_reads_threading(&graph_ref, &graph_reads, &params)
        .expect("raw thread");
    let mut pruned = raw.clone();
    let mut pruning = AssemblyGraphPruningParams::gatk_haplotype_caller_defaults();
    pruning.min_prune_factor = assembler.min_prune_factor;
    pruning.use_adaptive_pruning = assembler.use_adaptive_pruning;
    let _ = pruned.apply_pruning(&pruning);
    kv(
        "after_prune",
        format!(
            "nodes={} edges={}",
            pruned.node_count(),
            pruned.edge_count()
        ),
    );
    assert_eq!(pruned.node_count(), 1064);
    assert_eq!(pruned.edge_count(), 1079);
    assert!(rt_edge(&pruned, K971, K428).is_none());

    let dangling = DanglingRecoveryParams::from_assembler_args(&assembler);
    assert!(dangling.dangling_java_exact);
    let probes = pruned.probe_dangling_tail_failures(&dangling);
    kv("dangling_sinks_n", probes.len().to_string());
    let mut splice_dump = None;
    for (i, (sink, kmer, reason)) in probes.iter().enumerate() {
        kv(
            &format!("sink_{i}"),
            format!("id={sink} kmer={kmer} probe={reason}"),
        );
        let d = pruned.dangling_tail_decision_dump(*sink, &dangling);
        dump_one(&format!("decision_{i}"), &d);
        if d.cigar == "3I9M" && d.matching_suffix == 9 && d.last_ref_idx == 8 {
            splice_dump = Some(d);
        }
    }
    let splice = splice_dump.expect("canonical 3I9M lastRef=8 matchSuf=9 must exist");
    dump_one("splice_candidate", &splice);
    kv("splice_alt_bases", kmer_str(&splice.alt_bases));
    kv("splice_ref_bases_len", splice.ref_bases.len().to_string());
    kv(
        "splice_alt_kmers",
        splice
            .alt_path_ids
            .iter()
            .map(|&id| kmer_str(pruned.nodes()[id].kmer.as_ref()))
            .collect::<Vec<_>>()
            .join(","),
    );
    kv(
        "support_weight_is_standard",
        "Java mergeDanglingTail always addEdge(..., createEdge(false, 1)); support=1 is not a reject predicate",
    );

    // Java `mergeDanglingTail` uses signed int:
    //   refIndexToMerge = lastRefIndex - matchingSuffix + 1
    //   if (refIndexToMerge == 0) return 0;  // would cycle to LCA
    // Rust uses saturating_sub, so last_ref_idx < matching_suffix becomes 1, not 0.
    let java_ref_index = splice.last_ref_idx as i32 - splice.matching_suffix as i32 + 1;
    kv("java_considers_candidate", "true");
    kv("java_path_ok", "true");
    kv("java_cigar", splice.cigar.clone());
    kv("java_matching_suffix", splice.matching_suffix.to_string());
    kv("last_ref_idx", splice.last_ref_idx.to_string());
    kv("java_ref_index_to_merge", java_ref_index.to_string());
    kv(
        "rust_ref_index_to_merge",
        format!("{:?}", splice.ref_index_to_merge),
    );
    kv(
        "java_merge_predicate",
        "mergeDanglingTail: if (refIndexToMerge == 0) return 0",
    );
    assert_eq!(splice.cigar, "3I9M");
    assert_eq!(splice.matching_suffix, 9);
    assert_eq!(splice.last_ref_idx, 8);
    assert_eq!(splice.alt_index_to_merge, Some(2));
    assert_eq!(java_ref_index, 0, "Java signed arithmetic rejects at LCA");
    assert_eq!(
        splice.ref_index_to_merge,
        Some(0),
        "6R.239: checked arithmetic preserves Java path-index 0"
    );

    let mut dangled = pruned.clone();
    let summary = dangled
        .recover_dangling_branches(&dangling)
        .expect("dangling");
    kv("dangling_summary", format!("{summary:?}"));
    let after = rt_edge(&dangled, K971, K428);
    kv("after_dangling_971_428", format!("{after:?}"));
    assert!(
        after.is_none(),
        "6R.239: 971→428 must not appear after recover_dangling_tail"
    );

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
    let seq = SeqGraph::from_assembly_graph(&seq_asm);
    let join = seq
        .edges()
        .iter()
        .any(|e| e.support == 1 && !e.is_ref && seq.vertices()[e.to].sequence == b"T");
    kv("seqgraph_keeps_support1_join", join.to_string());
    assert!(!join, "6R.239: no support-1 splice onto T");

    let classification = "INDEX_SENTINEL_SEMANTICS_DIVERGENCE";
    kv("classification", classification);
    kv(
        "first_divergent_operation",
        "6R.239: Java mergeDanglingTail 8-9+1=0 now preserved; no addEdge(971,428).",
    );
    kv(
        "why_k128_is_downstream",
        "K=128 unchanged. After the splice is absent the carrier is inside K=128.",
    );
    kv(
        "support1_causal",
        "false — Java/Rust both use dangling addEdge weight 1; support is not a reject predicate",
    );

    assert_eq!(splice.reject_reason, Some("ref_index_zero_cycle"));
    assert!(splice.from_kmer.is_empty());
    assert!(splice.to_kmer.is_empty());
    assert_eq!(classification, "INDEX_SENTINEL_SEMANTICS_DIVERGENCE");
}
