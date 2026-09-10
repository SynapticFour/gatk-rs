//! 6R.140: experimental SeqGraph k-best resource policies on the canonical graph.
//! Skipped unless `HOLDOUT_6R140=1`. Production default remains `legacy_1024`.
//!
//! ```text
//! HOLDOUT_6R140=1 cargo test -p gatk-haplotypecaller --test holdout_6r140_kbest_resource_policy -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::assembly_region_finalize::{
    create_graph_reference_read, records_to_assembly_reads,
};
use gatk_haplotypecaller::kbest_haplotype::KBestPath;
use gatk_haplotypecaller::read_threading_assembler::build_threading_graph_for_seq_assembly;
use gatk_haplotypecaller::seq_graph::{SeqGraph, SeqGraphCleanupStatus};
use gatk_haplotypecaller::seq_kbest_haplotype::find_best_haplotypes_seq_graph_with_policy;
use gatk_haplotypecaller::seq_kbest_resource_policy::{
    SeqKbestResourcePolicy, DIAGNOSTIC_GENEROUS_BYTE_BUDGET,
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

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R140\t{key}\t{}", value.as_ref());
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

fn dump_run(
    label: &str,
    graph: &SeqGraph,
    java: &BTreeSet<String>,
    policy: SeqKbestResourcePolicy,
) -> (usize, bool, usize, usize, usize) {
    let report = find_best_haplotypes_seq_graph_with_policy(graph, 128, policy).expect(label);
    let set = unique_path_hashes(graph, &report.paths);
    let common = set.intersection(java).count();
    let java_only = java.difference(&set).count();
    let rust_only = set.difference(java).count();
    let d = &report.diagnostics;
    kv(
        "run",
        format!(
            "label={label}\tpolicy={}\trequested_k={}\tcompleted_paths={}\tresource_limit={}\tresource_limit_hit={}\tpeak_logical={}\tpeak_heap={}\tpeak_frontier_B={}\tpeak_pathstate_B={}\tpeak_payload_cap_B={}\tpaths_refused={}\texpansions={}\tjava_top_k_match={}\tjava_common={}\tjava_only={}\trust_only={}\trss_before={:?}\trss_peak={:?}\trss_after={:?}\taccounting={}",
            d.policy,
            d.requested_k,
            d.completed_paths,
            d.resource_limit,
            d.resource_limit_hit,
            d.peak_logical_frontier,
            d.peak_heap_entries,
            d.peak_frontier_bytes,
            d.peak_pathstate_bytes,
            d.peak_edge_payload_bytes,
            d.paths_refused,
            d.expansions,
            common == java.len() && rust_only == 0,
            common,
            java_only,
            rust_only,
            d.rss_before_mib,
            d.rss_peak_mib,
            d.rss_after_mib,
            d.accounting
        ),
    );
    (
        common,
        d.resource_limit_hit,
        d.peak_logical_frontier,
        d.peak_heap_entries,
        d.paths_refused,
    )
}

#[test]
fn holdout_6r140_kbest_resource_policy() {
    if std::env::var("HOLDOUT_6R140").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R140=1");
        return;
    }
    let java = load_java_k128(&repo_root().join(JAVA_KBEST_DUMP_REL));
    kv("java_k128_loaded", java.len().to_string());
    assert_eq!(java.len(), 128);
    let seq = build_seqgraph();

    let (legacy_common, legacy_hit, legacy_logical, _, _) = dump_run(
        "legacy_1024",
        &seq,
        &java,
        SeqKbestResourcePolicy::Legacy1024,
    );
    let (unb_common, unb_hit, unb_logical, _, _) = dump_run(
        "unbounded_diagnostic",
        &seq,
        &java,
        SeqKbestResourcePolicy::UnboundedDiagnostic,
    );
    let (byte_common, byte_hit, byte_logical, _, _) = dump_run(
        "generous_byte_budget",
        &seq,
        &java,
        SeqKbestResourcePolicy::ByteBudget {
            max_frontier_bytes: DIAGNOSTIC_GENEROUS_BYTE_BUDGET,
        },
    );

    kv(
        "summary",
        format!(
            "legacy_common={legacy_common}/128 hit={legacy_hit} logical={legacy_logical}\tunbounded_common={unb_common}/128 hit={unb_hit} logical={unb_logical}\tgenerous_byte_common={byte_common}/128 hit={byte_hit} logical={byte_logical}\tgenerous_budget_B={DIAGNOSTIC_GENEROUS_BYTE_BUDGET}\tnote=generous_budget_is_diagnostic_not_a_product_contract"
        ),
    );

    assert_eq!(legacy_common, 22, "legacy_1024 must reproduce COMMON=22");
    assert!(legacy_hit, "legacy_1024 must hit the resource limit");
    assert_eq!(unb_common, 128, "unbounded must match Java top-128");
    assert!(!unb_hit);
    assert_eq!(
        unb_logical, 1812,
        "unbounded logical peak is the 6R.136 Java-equivalent frontier"
    );
    assert_eq!(byte_common, 128, "generous byte budget must not constrain");
    assert!(!byte_hit);
    assert_eq!(byte_logical, 1812);
}
