//! 6R.238: Java `mergeDanglingTail` `refIndexToMerge == 0` is a path-index
//! sentinel (no merge / LCA cycle), not a graph vertex id. Rust
//! `saturating_sub` maps the canonical underflow onto path index 1.
//! Proof-only. Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! PRODUCTION CHANGE: NONE. Do not patch saturating_sub.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r238_dangling_index_sentinel_contract -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::assembly::{AssemblyGraph, AssemblyGraphPruningParams};
use gatk_haplotypecaller::assembly_based_caller::assemble_reads_with_finalized;
use gatk_haplotypecaller::assembly_dangling_recovery::DanglingRecoveryParams;
use gatk_haplotypecaller::assembly_region_finalize::{
    assembly_reference_read, create_graph_reference_read, records_to_assembly_reads,
};
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
const P11_REF_REL: &str = "parity/fixtures/p5_live_reference.fa";
const P11_BAM_REL: &str = "parity/build/sam-indexed-bam/p11_java_positive.bam";
const INDEL4_REF_REL: &str = "parity/fixtures/g2_subset_live_indel4.fa";
const INDEL4_BAM_REL: &str = "parity/build/sam-indexed-bam/g2_subset_live_indel4.bam";
const P12_BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const P12_INTERVAL: &str = "2:92307200-92307550";
const P12_TTC: u64 = 92_307_324;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R238\t{key}\t{}", value.as_ref());
}

/// Pinned Java `mergeDanglingTail` (SHA `2dbc0258`):
/// `refIndexToMerge = lastRefIndex - matchingSuffix + 1 + (leadingDel ? 1 : 0)`.
fn java_ref_index(last_ref_idx: i32, matching_suffix: i32, leading_del: bool) -> i32 {
    last_ref_idx - matching_suffix + 1 + i32::from(leading_del)
}

fn java_merges(last_ref_idx: i32, matching_suffix: i32, leading_del: bool) -> bool {
    java_ref_index(last_ref_idx, matching_suffix, leading_del) != 0
}

/// Current Rust `plan_dangling_tail_merge` arithmetic (not patched).
fn rust_ref_index(last_ref_idx: usize, matching_suffix: usize, leading_del: bool) -> usize {
    last_ref_idx.saturating_sub(matching_suffix) + 1 + usize::from(leading_del)
}

/// Java `longestSuffixMatch`: max return is `seqStart + 1` when the whole prefix of
/// `seq[0..=seqStart]` matches the k-mer suffix (or k-mer length if shorter).
fn java_longest_suffix_match_max(seq_start: i32) -> i32 {
    seq_start + 1
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

struct Scan {
    sinks: usize,
    underflow: usize,
    rust_path_idx_zero: usize,
}

fn prune_k25(
    region: &gatk_haplotypecaller::AssemblyRegion,
    root: &Path,
) -> Option<(AssemblyGraph, DanglingRecoveryParams)> {
    let ref_fasta = root.join(REF_REL);
    if !ref_fasta.is_file() {
        return None;
    }
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).ok()?;
    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
    let args = CallRegionArgs::strict_java();
    let mut owned = region.clone();
    let assembled =
        assemble_reads_with_finalized(&mut owned, &dict, &mut ref_cache, &args.assemble).ok()?;
    let padded = assembly_reference_read(&dict, &mut ref_cache, region).ok()?;
    let graph_ref = create_graph_reference_read(&padded, region, &dict);
    let graph_reads = records_to_assembly_reads(&assembled.finalized_reads);
    let assembler = args.assemble.assembler.clone();
    let params = threading_params(25, &assembler);
    let raw =
        assembly_graph_from_ref_and_reads_threading(&graph_ref, &graph_reads, &params).ok()?;
    let mut pruned = raw;
    let mut pruning = AssemblyGraphPruningParams::gatk_haplotype_caller_defaults();
    pruning.min_prune_factor = assembler.min_prune_factor;
    pruning.use_adaptive_pruning = assembler.use_adaptive_pruning;
    let _ = pruned.apply_pruning(&pruning);
    Some((
        pruned,
        DanglingRecoveryParams::from_assembler_args(&assembler),
    ))
}

fn scan_underflow(g: &AssemblyGraph, dangling: &DanglingRecoveryParams) -> Scan {
    let mut sinks = 0usize;
    let mut underflow = 0usize;
    let mut rust_path_idx_zero = 0usize;
    for (sink, _, _) in g.probe_dangling_tail_failures(dangling) {
        sinks += 1;
        let d = g.dangling_tail_decision_dump(sink, dangling);
        if d.matching_suffix > d.last_ref_idx {
            underflow += 1;
        }
        if d.ref_index_to_merge == Some(0) {
            rust_path_idx_zero += 1;
        }
    }
    Scan {
        sinks,
        underflow,
        rust_path_idx_zero,
    }
}

fn walk_active(
    root: &Path,
    ref_rel: &str,
    bam_rel: &str,
    interval: &str,
    target: Option<u64>,
) -> Option<gatk_haplotypecaller::AssemblyRegion> {
    let ref_fasta = root.join(ref_rel);
    let bam = root.join(bam_rel);
    if !ref_fasta.is_file() || !bam.is_file() {
        return None;
    }
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).ok()?;
    let specs = parse_intervals_cli_string(&dict, interval).ok()?;
    let walked = traverse_assembly_region_walker(
        &dict,
        &specs,
        &ref_fasta,
        &bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
    )
    .ok()?;
    flatten_assembly_regions(&walked).into_iter().find(|r| {
        matches!(
            call_disposition(r),
            AssemblyRegionCallDisposition::ActiveFull
        ) && target.map_or(true, |t| r.start.get() <= t && r.end.get() >= t)
    })
}

#[test]
fn forensic_6r238_java_refindex_arithmetic_contract() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    kv(
        "java_contract",
        "refIndexToMerge = lastRefIndex - matchingSuffix + 1 + (mustHandleLeadingDeletionCase ? 1 : 0); if (refIndexToMerge == 0) return 0; else addEdge(..., weight=1); return 1",
    );

    let rows: &[(i32, i32, i32, bool, usize)] = &[
        // lastRef, matchSuf, java_raw, java_merge, rust_idx
        (8, 9, 0, false, 1),
        (9, 9, 1, true, 1),
        (10, 9, 2, true, 2),
        (8, 8, 1, true, 1),
        (0, 1, 0, false, 1),
    ];
    kv(
        "boundary_header",
        "lastRefIndex\tmatchingSuffix\tJava raw\tJava merge?\tRust saturating",
    );
    for &(last, suf, java_raw, java_ok, rust) in rows {
        let got = java_ref_index(last, suf, false);
        assert_eq!(got, java_raw, "java {last},{suf}");
        assert_eq!(java_merges(last, suf, false), java_ok, "merge {last},{suf}");
        assert_eq!(
            rust_ref_index(last as usize, suf as usize, false),
            rust,
            "rust {last},{suf}"
        );
        kv(
            "boundary_row",
            format!("{last}\t{suf}\t{java_raw}\t{java_ok}\t{rust}"),
        );
    }
    let neg = java_ref_index(0, 2, false);
    assert_eq!(neg, -1);
    assert_ne!(neg, 0, "negative is not the ==0 sentinel");
    assert_eq!(rust_ref_index(0, 2, false), 1);
    kv(
        "boundary_row_unreachable_negative",
        "0\t2\t-1\tList.get would throw if reachable\t1",
    );

    // Domain: matchingSuffix <= lastRefIndex+1 ⇒ Java raw >= 0 (no leading del).
    for last in 0..16 {
        let max_suf = java_longest_suffix_match_max(last);
        assert_eq!(max_suf, last + 1);
        let min_raw = java_ref_index(last, max_suf, false);
        assert_eq!(min_raw, 0);
        assert!(!java_merges(last, max_suf, false));
        assert_eq!(
            rust_ref_index(last as usize, max_suf as usize, false),
            1,
            "Rust underflow at last={last} max_suf={max_suf} becomes path index 1"
        );
    }
    kv(
        "negative_raw",
        "UNREACHABLE from longestSuffixMatch: max matchingSuffix = lastRefIndex+1 so lastRefIndex-matchingSuffix+1 >= 0. Java if (==0) does not treat negative as sentinel; List.get(negative) would throw. lastRef=0 matchSuf=2 is outside the CIGAR/suffix domain.",
    );
    kv(
        "zero_meaning",
        "PATH INDEX 0 of referencePath = LCA = do not addEdge (cycle). Not a graph vertex id. Java recoverDanglingTail returns 0 (no merge) vs 1 (merged).",
    );
    kv(
        "established_sentinel_convention",
        "best_prefix_match_legacy maps Java -1 → Option::None. refIndexToMerge==0 is a second sentinel: no-merge, not Option graph-id.",
    );
}

#[test]
fn forensic_6r238_dangling_index_sentinel_contract() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    kv("classification", "INDEX_SENTINEL_SEMANTICS_DIVERGENCE");
    kv("closed_by", "6R.239 production checked arithmetic");

    let root = repo_root();
    let region =
        walk_active(&root, REF_REL, BAM_REL, INTERVAL, Some(TARGET)).expect("carrier ActiveFull");
    let (pruned, dangling) = prune_k25(&region, &root).expect("k=25 prune");
    kv(
        "carrier_prune",
        format!(
            "nodes={} edges={}",
            pruned.node_count(),
            pruned.edge_count()
        ),
    );
    let src = pruned.reference_source_vertex();
    kv(
        "carrier_graph_vertex_0_exists",
        (pruned.node_count() > 0).to_string(),
    );
    kv("carrier_ref_source_graph_id", format!("{src:?}"));
    assert!(src.is_some());
    // Graph id 0 can exist; the sentinel is PATH index 0, not graph id 0.
    let scan = scan_underflow(&pruned, &dangling);
    kv(
        "carrier_scan",
        format!(
            "sinks={} matchingSuffix_gt_lastRef={} rust_path_idx_0={}",
            scan.sinks, scan.underflow, scan.rust_path_idx_zero
        ),
    );
    assert!(
        scan.underflow >= 1,
        "canonical 3I9M must exercise matchingSuffix > lastRefIndex"
    );
    assert!(
        scan.rust_path_idx_zero >= 1,
        "6R.239: checked arithmetic yields path index 0 (Java sentinel)"
    );

    let mut saw_canonical = false;
    for (sink, _, _) in pruned.probe_dangling_tail_failures(&dangling) {
        let d = pruned.dangling_tail_decision_dump(sink, &dangling);
        if d.cigar == "3I9M" && d.matching_suffix == 9 && d.last_ref_idx == 8 {
            saw_canonical = true;
            let java = java_ref_index(8, 9, false);
            let rust = rust_ref_index(8, 9, false);
            kv(
                "canonical",
                format!(
                    "sink={} java_path_idx={java} rust_path_idx={rust} rust_opt={:?} lca={:?} dest={:?}",
                    d.sink, d.ref_index_to_merge, d.ref_path_ids.first(), d.to
                ),
            );
            assert_eq!(java, 0);
            assert!(!java_merges(8, 9, false));
            assert_eq!(rust, 1, "historical saturating helper still maps 8,9 → 1");
            assert_eq!(d.ref_index_to_merge, Some(0));
            assert_eq!(d.reject_reason, Some("ref_index_zero_cycle"));
            assert!(d.to.is_none(), "6R.239: sentinel does not pick ref_path[1]");
        }
    }
    assert!(
        saw_canonical,
        "3I9M lastRef=8 matchSuf=9 must still be present"
    );

    // Negative controls: do they rely on matchingSuffix > lastRefIndex?
    for (label, ref_rel, bam_rel, interval, tgt) in [
        ("l2_p11", P11_REF_REL, P11_BAM_REL, "chrLive:1-63", None),
        (
            "l2_indel4",
            INDEL4_REF_REL,
            INDEL4_BAM_REL,
            "chrIndel4:1-68",
            None,
        ),
        ("p12_ttc", REF_REL, P12_BAM_REL, P12_INTERVAL, Some(P12_TTC)),
    ] {
        let Some(reg) = walk_active(&root, ref_rel, bam_rel, interval, tgt) else {
            kv(label, "skip missing fixture");
            continue;
        };
        // p11/indel4 use their own reference fasta; reuse prune helper only for hs37d5.
        if ref_rel != REF_REL {
            let ref_fasta = root.join(ref_rel);
            let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
            let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
            let args = CallRegionArgs::strict_java();
            let mut owned = reg.clone();
            let assembled =
                assemble_reads_with_finalized(&mut owned, &dict, &mut ref_cache, &args.assemble)
                    .expect("assemble");
            let padded = assembly_reference_read(&dict, &mut ref_cache, &reg).expect("pad");
            let graph_ref = create_graph_reference_read(&padded, &reg, &dict);
            let graph_reads = records_to_assembly_reads(&assembled.finalized_reads);
            let assembler = args.assemble.assembler.clone();
            let k = if graph_ref.bases.len() < 50 { 10 } else { 25 };
            let params = threading_params(k, &assembler);
            let Ok(raw) =
                assembly_graph_from_ref_and_reads_threading(&graph_ref, &graph_reads, &params)
            else {
                kv(label, "skip graph build");
                continue;
            };
            let mut pruned = raw;
            let mut pruning = AssemblyGraphPruningParams::gatk_haplotype_caller_defaults();
            pruning.min_prune_factor = assembler.min_prune_factor;
            let _ = pruned.apply_pruning(&pruning);
            let dangling = DanglingRecoveryParams::from_assembler_args(&assembler);
            let s = scan_underflow(&pruned, &dangling);
            kv(
                label,
                format!(
                    "sinks={} underflow={} rust_path_idx_0={}",
                    s.sinks, s.underflow, s.rust_path_idx_zero
                ),
            );
        } else {
            let Some((pruned, dangling)) = prune_k25(&reg, &root) else {
                kv(label, "skip prune");
                continue;
            };
            let s = scan_underflow(&pruned, &dangling);
            kv(
                label,
                format!(
                    "sinks={} underflow={} rust_path_idx_0={}",
                    s.sinks, s.underflow, s.rust_path_idx_zero
                ),
            );
        }
    }

    kv(
        "callers",
        "plan_dangling_tail_merge ← recover_dangling_tail ← recover_dangling_branches ← build_threading_graph_core / audit_threading_dangling_recovery / dumps",
    );
    kv(
        "safe_candidate_unimplemented",
        "6R.239 production: checked (last+1).checked_sub(matching) so Ok(0) feeds ref_index_zero_cycle. Historical saturating helper remains in this file as the pre-fix contract.",
    );
    kv(
        "why_k128_is_downstream",
        "Index sentinel only; K=128 unchanged",
    );
}
