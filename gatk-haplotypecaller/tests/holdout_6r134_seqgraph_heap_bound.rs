//! 6R.134 forensic: first expansion refused by `MAX_KBEST_HEAP_PATHS` on the matched k=25 graph.
//! Skipped unless `HOLDOUT_6R134=1`. Production change: NONE.
//!
//! ```text
//! HOLDOUT_6R134=1 cargo test -p gatk-haplotypecaller --test holdout_6r134_seqgraph_heap_bound -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::assembly_region_finalize::{
    create_graph_reference_read, records_to_assembly_reads,
};
use gatk_haplotypecaller::kbest_haplotype::KBestPath;
use gatk_haplotypecaller::read_threading_assembler::build_threading_graph_for_seq_assembly;
use gatk_haplotypecaller::seq_graph::{SeqGraph, SeqGraphCleanupStatus};
use gatk_haplotypecaller::seq_kbest_haplotype::{
    find_best_haplotypes_seq_graph, find_best_haplotypes_seq_graph_forensic, SeqKbestCapPolicy,
    SeqKbestExpandRefuse, SEQ_KBEST_PRODUCTION_MAX_EXPANSIONS, SEQ_KBEST_PRODUCTION_MAX_HEAP,
    SEQ_KBEST_PRODUCTION_MAX_PATH_EDGES,
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
const HEAP_BOUNDS: &[usize] = &[1024, 1280, 1536, 2048, 4096];

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R134\t{key}\t{}", value.as_ref());
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

fn compare(label: &str, rust: &BTreeSet<String>, java: &BTreeSet<String>) {
    kv(
        "compare",
        format!(
            "label={label}\trust={}\tjava={}\tCOMMON={}\tJAVA_ONLY={}\tRUST_ONLY={}",
            rust.len(),
            java.len(),
            rust.intersection(java).count(),
            java.difference(rust).count(),
            rust.difference(java).count()
        ),
    );
}

fn load_java_k128(path: &Path) -> BTreeSet<String> {
    let mut kbest = BTreeSet::new();
    let mut assemble = BTreeSet::new();
    let Ok(text) = std::fs::read_to_string(path) else {
        return kbest;
    };
    for line in text.lines() {
        let mut parts = line.split('\t');
        let Some(tag) = parts.next() else { continue };
        if tag != "6R133" {
            continue;
        }
        let Some(kind) = parts.next() else { continue };
        match kind {
            "kbest" => {
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
                        kbest.insert(h);
                    }
                }
            }
            "assemble" => {
                for tok in parts {
                    if let Some(v) = tok.strip_prefix("hash=") {
                        assemble.insert(v.to_string());
                    }
                }
            }
            _ => {}
        }
    }
    if !kbest.is_empty() {
        kbest
    } else {
        assemble
    }
}

fn heap_policy(max_heap: usize) -> SeqKbestCapPolicy {
    SeqKbestCapPolicy {
        max_heap_paths: Some(max_heap),
        max_expansions: Some(SEQ_KBEST_PRODUCTION_MAX_EXPANSIONS),
        max_path_edges: Some(SEQ_KBEST_PRODUCTION_MAX_PATH_EDGES),
        max_heap_paths_per_dest: None,
    }
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

fn refuse_hits<'a>(rec: &SeqKbestExpandRefuse, paths: &'a [KBestPath]) -> Vec<(usize, f64, usize)> {
    rec.refused_tos
        .iter()
        .map(|&(to, score, _)| {
            let extra = (rec.parent_last, to);
            let n = paths
                .iter()
                .filter(|p| edges_have_prefix(&p.edges, &rec.parent_edges, extra))
                .count();
            (to, score, n)
        })
        .collect()
}

fn build_seqgraph() -> (SeqGraph, BTreeSet<String>) {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    let java = load_java_k128(&root.join(JAVA_KBEST_DUMP_REL));
    kv("java_k128_loaded", java.len().to_string());

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
    (seq, java)
}

#[test]
fn holdout_6r134_seqgraph_heap_bound() {
    if std::env::var("HOLDOUT_6R134").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R134=1");
        return;
    }
    kv(
        "prod_caps",
        format!(
            "heap={SEQ_KBEST_PRODUCTION_MAX_HEAP}\texpansions={SEQ_KBEST_PRODUCTION_MAX_EXPANSIONS}\tpath_edges={SEQ_KBEST_PRODUCTION_MAX_PATH_EDGES}"
        ),
    );
    let (seq, java) = build_seqgraph();

    let unbounded = find_best_haplotypes_seq_graph_forensic(
        &seq,
        128,
        128,
        SeqKbestCapPolicy::unbounded(),
        &[],
    )
    .expect("unbounded");
    let unb_set = unique_path_hashes(&seq, &unbounded.paths);
    compare("unbounded_k128_vs_java", &unb_set, &java);
    kv(
        "unbounded",
        format!(
            "n={}\tmax_heap={}\tpops={}\texpansions={}\tskip_heap_pop={}\tskip_heap_exp={}",
            unbounded.paths.len(),
            unbounded.max_heap,
            unbounded.pop_count,
            unbounded.expansions,
            unbounded.skip_heap_full_at_pop,
            unbounded.skip_heap_full_at_expand
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
    let prod_fn = find_best_haplotypes_seq_graph(&seq, 128).expect("prod fn");
    let prod_set = unique_path_hashes(&seq, &prod.paths);
    let prod_fn_set = unique_path_hashes(&seq, &prod_fn);
    compare("prod_forensic_vs_prod_fn", &prod_set, &prod_fn_set);
    compare("prod_k128_vs_java", &prod_set, &java);
    kv(
        "prod",
        format!(
            "n={}\tmax_heap={}\tpops={}\texpansions={}\tskip_heap_pop={}\tskip_heap_exp={}\theap_first_hit_cap_pop={:?}\trefuses={}",
            prod.paths.len(),
            prod.max_heap,
            prod.pop_count,
            prod.expansions,
            prod.skip_heap_full_at_pop,
            prod.skip_heap_full_at_expand,
            prod.heap_first_hit_cap_pop,
            prod.expand_refuses.len()
        ),
    );

    let java_only: BTreeSet<_> = java.difference(&prod_set).cloned().collect();
    let rust_only: BTreeSet<_> = prod_set.difference(&java).cloned().collect();
    kv(
        "private",
        format!(
            "JAVA_ONLY={}\tRUST_ONLY={}",
            java_only.len(),
            rust_only.len()
        ),
    );

    if let Some(rec) = prod.first_expand_refuse.as_ref() {
        let parent_h = path_hash(&seq, rec.parent_start, &rec.parent_edges);
        kv(
            "first_refuse",
            format!(
                "pop={}\texpansions={}\tresult_len={}\theap_len={}\tparent_score={:.9}\tparent_edges={}\tparent_last={}\tparent_hash={parent_h}\tn_outs={}\tinserted_this={}\trefused={}",
                rec.pop_count,
                rec.expansions,
                rec.result_len,
                rec.heap_len_at_refuse,
                rec.parent_score,
                rec.parent_n_edges,
                rec.parent_last,
                rec.n_outs,
                rec.n_inserted_this_expand,
                rec.refused_tos.len()
            ),
        );
        let hits_unb = refuse_hits(rec, &unbounded.paths);
        let hits_prod = refuse_hits(rec, &prod.paths);
        for (i, ((to, score, n_unb), (_, _, n_prod))) in
            hits_unb.iter().zip(hits_prod.iter()).enumerate()
        {
            let mut edges = rec.parent_edges.clone();
            edges.push((rec.parent_last, *to));
            let sh = path_hash(&seq, rec.parent_start, &edges);
            let in_java = unbounded
                .paths
                .iter()
                .filter(|p| edges_have_prefix(&p.edges, &rec.parent_edges, (rec.parent_last, *to)))
                .filter_map(|p| {
                    let h = path_hash(&seq, p.start, &p.edges);
                    java.contains(&h).then_some(h)
                })
                .collect::<BTreeSet<_>>();
            let in_java_only = in_java.iter().filter(|h| java_only.contains(*h)).count();
            kv(
                "first_refuse_succ",
                format!(
                    "idx={i}\tto={to}\tscore={score:.9}\thash={sh}\tunbounded_completions={n_unb}\tprod_completions={n_prod}\tjava_top128={}\tjava_only={in_java_only}",
                    in_java.len()
                ),
            );
        }
    } else {
        kv("first_refuse", "none");
    }

    let mut causal: Option<(usize, String)> = None;
    for (ri, rec) in prod.expand_refuses.iter().enumerate() {
        for &(to, score, _) in &rec.refused_tos {
            extra_loop(
                &seq,
                rec,
                to,
                score,
                &unbounded.paths,
                &java,
                &java_only,
                ri,
                &mut causal,
            );
        }
        if causal.is_some() {
            break;
        }
    }
    match &causal {
        Some((ri, detail)) => kv("causal_refuse", format!("refuse_index={ri}\t{detail}")),
        None => kv("causal_refuse", "none_in_expand_refuses"),
    }

    for &bound in HEAP_BOUNDS {
        let r = find_best_haplotypes_seq_graph_forensic(&seq, 128, 128, heap_policy(bound), &[])
            .expect("bound");
        let set = unique_path_hashes(&seq, &r.paths);
        kv(
            "counterfactual",
            format!(
                "bound={bound}\tmax_heap={}\trefused_exp={}\tskip_pop={}\tn={}\tCOMMON={}\tJAVA_ONLY={}\tRUST_ONLY={}",
                r.max_heap,
                r.skip_heap_full_at_expand,
                r.skip_heap_full_at_pop,
                r.paths.len(),
                set.intersection(&java).count(),
                java.difference(&set).count(),
                set.difference(&java).count()
            ),
        );
    }
}

fn extra_loop(
    seq: &SeqGraph,
    rec: &SeqKbestExpandRefuse,
    to: usize,
    score: f64,
    unbounded: &[KBestPath],
    java: &BTreeSet<String>,
    java_only: &BTreeSet<String>,
    ri: usize,
    causal: &mut Option<(usize, String)>,
) {
    if causal.is_some() {
        return;
    }
    let extra = (rec.parent_last, to);
    let mut java_hits = Vec::new();
    for p in unbounded {
        if !edges_have_prefix(&p.edges, &rec.parent_edges, extra) {
            continue;
        }
        let h = path_hash(seq, p.start, &p.edges);
        if java.contains(&h) {
            java_hits.push(h);
        }
    }
    let n_only = java_hits.iter().filter(|h| java_only.contains(*h)).count();
    if n_only > 0 {
        let parent_h = path_hash(seq, rec.parent_start, &rec.parent_edges);
        *causal = Some((
            ri,
            format!(
                "pop={}\texpansions={}\theap={}\tparent_hash={parent_h}\tparent_score={:.9}\tto={to}\tsucc_score={score:.9}\tjava_top128={}\tjava_only={n_only}\thit0={}",
                rec.pop_count,
                rec.expansions,
                rec.heap_len_at_refuse,
                rec.parent_score,
                java_hits.len(),
                java_hits.first().map(|s| s.as_str()).unwrap_or(".")
            ),
        ));
    }
}
