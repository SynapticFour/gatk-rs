//! 6R.138 forensic: k-best frontier bytes vs process RSS (canonical + bomb window).
//! Skipped unless `HOLDOUT_6R138=1`. Production change: NONE.
//!
//! ```text
//! HOLDOUT_6R138=1 cargo test -p gatk-haplotypecaller --test holdout_6r138_kbest_rss -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::assembly_region_finalize::{
    create_graph_reference_read, records_to_assembly_reads,
};
use gatk_haplotypecaller::kbest_haplotype::KBestPath;
use gatk_haplotypecaller::read_threading_assembler::build_threading_graph_for_seq_assembly;
use gatk_haplotypecaller::runtime_config::current_rss_mib;
use gatk_haplotypecaller::seq_graph::{SeqGraph, SeqGraphCleanupStatus};
use gatk_haplotypecaller::seq_kbest_compact_forensic::{
    compact_heap_item_size, compact_node_size, find_best_haplotypes_seq_graph_compact,
};
use gatk_haplotypecaller::seq_kbest_haplotype::{
    find_best_haplotypes_seq_graph_forensic, seq_kbest_production_layout, SeqKbestCapPolicy,
    SeqKbestFrontierMem,
};
use gatk_haplotypecaller::{
    assemble_reads_with_finalized, call_disposition, flatten_assembly_regions,
    traverse_assembly_region_walker, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    ReadFilterParams, WalkerTraversalConfig,
};
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

const CANON_INTERVAL: &str = "20:29455000-29456500";
const CANON_BAM: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_KBEST_DUMP_REL: &str = "parity/giab/runs/local-pairhmm-diff/6r133_java.txt";
const TARGET: u64 = 29_456_196;
const JAVA_ACTIVE_START: u64 = 29_455_995;
const JAVA_ACTIVE_END: u64 = 29_456_293;
const JAVA_PAD_START: u64 = 29_455_895;
const JAVA_PAD_END: u64 = 29_456_393;
const EXPECTED_REF_HASH: &str = "659741f99b7f78a5";
const EXPECTED_FINALIZED: usize = 379;

const BOMB_INTERVAL: &str = "20:10098000-10099600";
const BOMB_BAM: &str =
    "parity/realworld/na12878_giab_window_mem_500kb_b37/NA12878_giab_window.b37.bam";

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R138\t{key}\t{}", value.as_ref());
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

fn fmt_rss(v: Option<f64>) -> String {
    v.map(|x| format!("{x:.3}")).unwrap_or_else(|| "NA".into())
}

fn rss_delta(before: Option<f64>, peak: Option<f64>) -> String {
    match (before, peak) {
        (Some(b), Some(p)) => format!("{:.3}", p - b),
        _ => "NA".into(),
    }
}

fn dump_mem(label: &str, mem: &SeqKbestFrontierMem) {
    kv(
        "mem",
        format!(
            "label={label}\tpeak_logical={}\tpeak_heap_len={}\tpeak_heap_cap={}\tpeak_edge_len={}\tpeak_edge_cap={}\tpeak_sum_edge_len={}\tpeak_sum_edge_cap={}\tpeak_inline_B={}\tpeak_payload_len_B={}\tpeak_payload_cap_B={}\tpeak_frontier_B={}\tpeak_result_n={}\tpeak_result_payload_cap_B={}\tpeak_with_results_B={}\trss_before_MiB={}\trss_peak_MiB={}\trss_after_MiB={}\trss_delta_peak_minus_before_MiB={}\trss_samples={}\trss_aborted={}",
            mem.peak_logical,
            mem.peak_heap_len,
            mem.peak_heap_capacity,
            mem.peak_edge_len,
            mem.peak_edge_cap,
            mem.peak_sum_edge_len,
            mem.peak_sum_edge_cap,
            mem.peak_inline_bytes,
            mem.peak_payload_len_bytes,
            mem.peak_payload_cap_bytes,
            mem.peak_frontier_bytes,
            mem.peak_result_n,
            mem.peak_result_payload_cap_bytes,
            mem.peak_frontier_plus_results_bytes,
            fmt_rss(mem.rss_before_mib),
            fmt_rss(mem.rss_peak_mib),
            fmt_rss(mem.rss_after_mib),
            rss_delta(mem.rss_before_mib, mem.rss_peak_mib),
            mem.rss_samples,
            mem.rss_aborted
        ),
    );
}

fn build_canon_seqgraph() -> SeqGraph {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(CANON_BAM);
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, CANON_INTERVAL).expect("interval");
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
        "canon_matched_input",
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
        "canon_seqgraph",
        format!("nodes={}\tedges={}", seq.node_count(), seq.edge_count()),
    );
    seq
}

fn rss_calibration() {
    let before = current_rss_mib();
    let blob = vec![0xABu8; 16 * 1024 * 1024];
    std::hint::black_box(&blob);
    let after = current_rss_mib();
    kv(
        "rss_calibration_16MiB",
        format!(
            "before_MiB={}\tafter_MiB={}\tdelta_MiB={}\tnote=process_RSS_not_allocator",
            fmt_rss(before),
            fmt_rss(after),
            rss_delta(before, after)
        ),
    );
    drop(blob);
}

fn run_copied(
    label: &str,
    graph: &SeqGraph,
    java: Option<&BTreeSet<String>>,
    caps: SeqKbestCapPolicy,
) -> SeqKbestFrontierMem {
    let report = find_best_haplotypes_seq_graph_forensic(graph, 128, 128, caps, &[]).expect(label);
    let set = unique_path_hashes(graph, &report.paths);
    let common = java.map(|j| set.intersection(j).count());
    kv(
        "copied_run",
        format!(
            "label={label}\tn={}\tmax_heap={}\texpansions={}\tCOMMON={}\trss_aborted={}",
            report.paths.len(),
            report.max_heap,
            report.expansions,
            common
                .map(|c| c.to_string())
                .unwrap_or_else(|| "n/a".into()),
            report.mem.rss_aborted
        ),
    );
    dump_mem(label, &report.mem);
    report.mem
}

fn run_compact(
    label: &str,
    graph: &SeqGraph,
    java: Option<&BTreeSet<String>>,
    lazy: bool,
    cap: Option<usize>,
) -> SeqKbestFrontierMem {
    let report = find_best_haplotypes_seq_graph_compact(graph, 128, 128, lazy, cap).expect(label);
    let set = unique_path_hashes(graph, &report.paths);
    let common = java.map(|j| set.intersection(j).count());
    kv(
        "compact_run",
        format!(
            "label={label}\tmode={}\tn={}\tpeak_heap={}\tpeak_logical={}\tarena={}\tCOMMON={}\trss_aborted={}",
            report.mode,
            report.paths.len(),
            report.peak_heap,
            report.peak_logical_unpopped,
            report.arena_nodes,
            common
                .map(|c| c.to_string())
                .unwrap_or_else(|| "n/a".into()),
            report.mem.rss_aborted
        ),
    );
    dump_mem(label, &report.mem);
    report.mem
}

fn bomb_graphs() -> Vec<(String, SeqGraph)> {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BOMB_BAM);
    if !bam.exists() || !ref_fasta.exists() {
        kv(
            "bomb_skip",
            format!("missing bam={} ref={}", bam.display(), ref_fasta.display()),
        );
        return Vec::new();
    }
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, BOMB_INTERVAL).expect("interval");
    let walk = traverse_assembly_region_walker(
        &dict,
        &specs,
        &ref_fasta,
        &bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
    )
    .expect("bomb walk");
    let regions = flatten_assembly_regions(&walk);
    kv("bomb_regions", format!("n={}", regions.len()));
    let args = CallRegionArgs::strict_java();
    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
    let mut out = Vec::new();
    for (ri, region) in regions.iter().enumerate() {
        if !matches!(
            call_disposition(region),
            AssemblyRegionCallDisposition::ActiveFull
        ) {
            continue;
        }
        let mut owned = region.clone();
        let assembled = match assemble_reads_with_finalized(
            &mut owned,
            &dict,
            &mut ref_cache,
            &args.assemble,
        ) {
            Ok(a) => a,
            Err(e) => {
                kv("bomb_assemble_err", format!("ri={ri}\terr={e}"));
                continue;
            }
        };
        let padded_ref =
            match gatk_haplotypecaller::assembly_region_finalize::assembly_reference_read(
                &dict,
                &mut ref_cache,
                &owned,
            ) {
                Ok(r) => r,
                Err(e) => {
                    kv("bomb_pad_err", format!("ri={ri}\terr={e}"));
                    continue;
                }
            };
        let graph_ref = create_graph_reference_read(&padded_ref, &owned, &dict);
        let graph_reads = records_to_assembly_reads(&assembled.finalized_reads);
        let mut assembler = args.assemble.assembler.clone();
        assembler.dangling_java_exact = true;
        for kmer in [10usize, 25] {
            let rt = match build_threading_graph_for_seq_assembly(
                &graph_ref,
                &graph_reads,
                kmer,
                &assembler,
                false,
                false,
            ) {
                Ok(Some(g)) => g,
                Ok(None) => continue,
                Err(e) => {
                    kv("bomb_rt_err", format!("ri={ri}\tk={kmer}\terr={e}"));
                    continue;
                }
            };
            let mut seq = SeqGraph::from_assembly_graph(&rt);
            seq.clean_non_ref_paths();
            let status = seq.cleanup_seq_graph();
            kv(
                "bomb_seqgraph",
                format!(
                    "ri={ri}\tspan={}:{}-{}\tk={kmer}\tnodes={}\tedges={}\tstatus={status:?}\treads={}",
                    owned.contig,
                    owned.start.get(),
                    owned.end.get(),
                    seq.node_count(),
                    seq.edge_count(),
                    assembled.finalized_reads.len()
                ),
            );
            if seq.node_count() == 0 {
                continue;
            }
            out.push((
                format!(
                    "ri{ri}_k{kmer}_{}_{}-{}",
                    owned.contig,
                    owned.start.get(),
                    owned.end.get()
                ),
                seq,
            ));
        }
    }
    kv("bomb_graphs_built", format!("n={}", out.len()));
    out
}

#[test]
fn holdout_6r138_kbest_rss() {
    if std::env::var("HOLDOUT_6R138").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R138=1");
        return;
    }

    let lay = seq_kbest_production_layout();
    kv(
        "layout",
        format!(
            "PathState={}\tHeapItem={}\tedge_pair={}\tCompactNode={}\tCompactHeapItem={}\tVec_hdr={}",
            lay.path_state_bytes,
            lay.heap_item_bytes,
            lay.edge_pair_bytes,
            compact_node_size(),
            compact_heap_item_size(),
            std::mem::size_of::<Vec<(usize, usize)>>()
        ),
    );
    kv(
        "tightness_proxy",
        format!(
            "1024x4096x16B={}\tthreshold=128MiB\ttest=kbest_bounds_are_finite_and_tight",
            1024usize.saturating_mul(4096).saturating_mul(16)
        ),
    );

    rss_calibration();

    let java = load_java_k128(&repo_root().join(JAVA_KBEST_DUMP_REL));
    kv("java_k128_loaded", java.len().to_string());
    let seq = build_canon_seqgraph();

    let mut copied_deltas = Vec::new();
    for trial in 0..3 {
        let mem = run_copied(
            &format!("canon_copied_unbounded_t{trial}"),
            &seq,
            Some(&java),
            SeqKbestCapPolicy::unbounded(),
        );
        copied_deltas.push(rss_delta(mem.rss_before_mib, mem.rss_peak_mib));
    }
    kv("canon_copied_unbounded_rss_deltas", copied_deltas.join(","));

    run_copied(
        "canon_copied_production",
        &seq,
        Some(&java),
        SeqKbestCapPolicy::production(),
    );
    run_compact(
        "canon_compact_eager_unbounded",
        &seq,
        Some(&java),
        false,
        None,
    );
    run_compact(
        "canon_compact_lazy_unbounded",
        &seq,
        Some(&java),
        true,
        None,
    );

    let graphs = bomb_graphs();
    for (name, g) in &graphs {
        let prod = run_copied(
            &format!("bomb_copied_production_{name}"),
            g,
            None,
            SeqKbestCapPolicy::production(),
        );
        run_compact(
            &format!("bomb_compact_eager_productioncap_{name}"),
            g,
            None,
            false,
            Some(1024),
        );
        // Unbounded copied on a graph that already saturates 1024 can be the
        // historical multi-GiB failure. Only run it when production stayed small.
        if prod.peak_heap_len < 1024 && g.node_count() <= 500 {
            run_copied(
                &format!("bomb_copied_unbounded_{name}"),
                g,
                None,
                SeqKbestCapPolicy::unbounded(),
            );
            run_compact(
                &format!("bomb_compact_eager_unbounded_{name}"),
                g,
                None,
                false,
                None,
            );
        } else {
            kv(
                "bomb_unbounded_skipped",
                format!(
                    "name={name}\treason=production_heap_or_nodes_too_large\tpeak_heap={}\tnodes={}",
                    prod.peak_heap_len,
                    g.node_count()
                ),
            );
        }
    }
}
