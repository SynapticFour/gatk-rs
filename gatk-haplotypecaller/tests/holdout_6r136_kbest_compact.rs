//! 6R.136 forensic: compact ancestry vs lazy heap occupancy vs Java top-128.
//! Skipped unless `HOLDOUT_6R136=1`. Production change: NONE.
//!
//! ```text
//! HOLDOUT_6R136=1 cargo test -p gatk-haplotypecaller --test holdout_6r136_kbest_compact -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::assembly_region_finalize::{
    create_graph_reference_read, records_to_assembly_reads,
};
use gatk_haplotypecaller::kbest_haplotype::KBestPath;
use gatk_haplotypecaller::read_threading_assembler::build_threading_graph_for_seq_assembly;
use gatk_haplotypecaller::seq_graph::{SeqGraph, SeqGraphCleanupStatus};
use gatk_haplotypecaller::seq_kbest_compact_forensic::{
    compact_heap_item_size, compact_node_size, copied_edge_payload_bytes,
    find_best_haplotypes_seq_graph_compact,
};
use gatk_haplotypecaller::seq_kbest_haplotype::{
    find_best_haplotypes_seq_graph_forensic, seq_kbest_production_layout, SeqKbestCapPolicy,
    SEQ_KBEST_PRODUCTION_MAX_HEAP,
};
use gatk_haplotypecaller::{
    assemble_reads_with_finalized, call_disposition, flatten_assembly_regions,
    traverse_assembly_region_walker, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    ReadFilterParams, WalkerTraversalConfig,
};
use std::collections::BTreeSet;
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
    println!("6R136\t{key}\t{}", value.as_ref());
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

fn load_java_k128(path: &Path) -> BTreeSet<String> {
    let mut set = BTreeSet::new();
    let Ok(text) = std::fs::read_to_string(path) else {
        return set;
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
        let mut hash = None;
        for tok in parts {
            if let Some(v) = tok.strip_prefix("kmer=") {
                kmer = v.parse().ok();
            }
            if let Some(v) = tok.strip_prefix("K=") {
                kcap = v.parse().ok();
            }
            if let Some(v) = tok.strip_prefix("hash=") {
                hash = Some(v.to_string());
            }
        }
        if kmer == Some(25) && kcap == Some(128) {
            if let Some(h) = hash {
                set.insert(h);
            }
        }
    }
    set
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

fn dump_compact(
    label: &str,
    graph: &SeqGraph,
    java: &BTreeSet<String>,
    report: &gatk_haplotypecaller::seq_kbest_compact_forensic::SeqKbestCompactReport,
) {
    let set = unique_path_hashes(graph, &report.paths);
    kv(
        "run",
        format!(
            "label={label}\tmode={}\tn={}\tpeak_heap={}\tpeak_live={}\tpeak_deferred={}\tpeak_logical={}\tarena={}\trecon={}\tpops={}\texpansions={}\tskip_heap={}\tCOMMON={}\tJAVA_ONLY={}\tRUST_ONLY={}\thas_rank56={}",
            report.mode,
            report.paths.len(),
            report.peak_heap,
            report.peak_live_unpopped,
            report.peak_deferred_sibs,
            report.peak_logical_unpopped,
            report.arena_nodes,
            report.reconstructions,
            report.pops,
            report.expansions,
            report.skip_heap_full,
            set.intersection(java).count(),
            java.difference(&set).count(),
            set.difference(java).count(),
            set.contains(JAVA_RANK56)
        ),
    );
}

#[test]
fn holdout_6r136_kbest_compact() {
    if std::env::var("HOLDOUT_6R136").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R136=1");
        return;
    }
    let lay = seq_kbest_production_layout();
    kv(
        "layout",
        format!(
            "PathState={}\tHeapItem={}\tedge_pair={}\tCompactNode={}\tCompactHeapItem={}\tVec_hdr=24\tcopied_29_edges={}\tcopied_4096_edges={}\tcompact_1812={}",
            lay.path_state_bytes,
            lay.heap_item_bytes,
            lay.edge_pair_bytes,
            compact_node_size(),
            compact_heap_item_size(),
            copied_edge_payload_bytes(29),
            copied_edge_payload_bytes(4096),
            compact_node_size() * 1812
        ),
    );
    kv(
        "rss_proxy",
        format!(
            "prod_1024x29_payload={}\tprod_1024x4096_payload={}\tcompact_1812_nodes={}\theap_items_1812={}",
            1024 * copied_edge_payload_bytes(29) + 1024 * lay.heap_item_bytes,
            1024 * copied_edge_payload_bytes(4096) + 1024 * lay.heap_item_bytes,
            compact_node_size() * 1812,
            compact_heap_item_size() * 1812
        ),
    );

    let java = load_java_k128(&repo_root().join(JAVA_KBEST_DUMP_REL));
    kv("java_k128_loaded", java.len().to_string());
    let seq = build_seqgraph();

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
        "copied_unbounded",
        format!(
            "n={}\tmax_heap={}\tCOMMON={}\thas_rank56={}",
            unbounded.paths.len(),
            unbounded.max_heap,
            unb_set.intersection(&java).count(),
            unb_set.contains(JAVA_RANK56)
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
    kv(
        "copied_production",
        format!(
            "n={}\tmax_heap={}\tskip_exp={}\tCOMMON={}\thas_rank56={}",
            prod.paths.len(),
            prod.max_heap,
            prod.skip_heap_full_at_expand,
            prod_set.intersection(&java).count(),
            prod_set.contains(JAVA_RANK56)
        ),
    );

    let eager = find_best_haplotypes_seq_graph_compact(&seq, 128, 128, false, None).expect("eager");
    dump_compact("compact_eager_unbounded", &seq, &java, &eager);

    let lazy = find_best_haplotypes_seq_graph_compact(&seq, 128, 128, true, None).expect("lazy");
    dump_compact("compact_lazy_unbounded", &seq, &java, &lazy);

    let eager_cap = find_best_haplotypes_seq_graph_compact(
        &seq,
        128,
        128,
        false,
        Some(SEQ_KBEST_PRODUCTION_MAX_HEAP),
    )
    .expect("eager cap");
    dump_compact("compact_eager_cap1024", &seq, &java, &eager_cap);

    let lazy_cap = find_best_haplotypes_seq_graph_compact(
        &seq,
        128,
        128,
        true,
        Some(SEQ_KBEST_PRODUCTION_MAX_HEAP),
    )
    .expect("lazy cap");
    dump_compact("compact_lazy_cap1024", &seq, &java, &lazy_cap);

    kv(
        "eager_vs_unbounded_set",
        format!(
            "equal={}",
            unique_path_hashes(&seq, &eager.paths) == unb_set
        ),
    );
    kv(
        "lazy_vs_unbounded_set",
        format!("equal={}", unique_path_hashes(&seq, &lazy.paths) == unb_set),
    );
}
