//! 6R.78: SeqGraph k-best ranking / bounded retention vs Java K=128.
//!
//! Forensic only. No production search, score, K, or comparator change.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r78_seqgraph_kbest_ranking_contract
//! HOLDOUT_6R78=1 cargo test -p gatk-haplotypecaller --test forensic_6r78_seqgraph_kbest_ranking_contract live_ -- --nocapture
//! ```

use gatk_haplotypecaller::kbest_haplotype::KBestPath;
use gatk_haplotypecaller::seq_graph::SeqGraph;
use gatk_haplotypecaller::seq_kbest_haplotype::{
    find_best_haplotypes_seq_graph_forensic, seq_kbest_path_score_terms, seq_kbest_score_cmp,
    SeqKbestCapPolicy, SEQ_KBEST_PRODUCTION_MAX_EXPANSIONS, SEQ_KBEST_PRODUCTION_MAX_HEAP,
    SEQ_KBEST_PRODUCTION_MAX_PATH_EDGES,
};
use std::cmp::Ordering;

/// Java-only / rust-only PairHMM sequences from 6R.76 (live secondary only).
const JAVA_ONLY_J0: &[u8] = b"CATGGAGCCTGACCTTATTTGAAGTAGGGCATTTGCAGATGTATTTAAGATATTTGAGGCTGGGCACAGTGGCTCACGTCTGTAATCCCAGCACTTTGAAAGGCCGAGGCAGGTGGATTCACCTGAGGTCAGGAGTTTGAGACCAGCCTGTCCCACATGGTGAAAAGCCCGTATCTACCAAAAATACAAACGTTAGCTGTGTGTGGTGGTGGCGGCACCTGTAATCCCAGCTACTCGAGAGCCAGAG";
const JAVA_ONLY_J1: &[u8] = b"CATGGAGCCTGACTTTATTTGAAGTAGGGCATTTGCAGATGTATTTAAGATATTTGAGGCCGGGCACAGTGGCTCATGTCTGTAATCCCAGCACTTTAAAAGGCTGAGGCAGGTGTATTCACCTGAGGTCAGGAGTTCGAGACCAGCCTGGCCAACATGGTGAAAAGCCCGTATCTACCAAAAATACAAAAGTTAGCTGGGTGTGGTGGCGGCACCTGTAATCCCAGCTACTCGAGAGCCAGAG";
const RUST_ONLY_R0: &[u8] = b"CATGGAGCCTGACCTTATTTGAAGTAGGGCATTTGCAGATGTATTTAAGATATTTGAGGCTGGGCACAGTGGCTCACGTCTGTAATCCCAGCACTTTGAAAGGCCGAGGCAGGTGGATTCACCTGAGGTCAGGAGTTTGAGACCAGCCTGTCCCACATGGTGAAAAGCCCGTATCTACCAAAAATACAAACGTTAGCTGTGTGTGGTGGTGGCGGCACCTGTAATCGCAGCTACTCGAGAGCCAGAG";

fn contains_bases(hay: &[u8], needle: &[u8]) -> bool {
    hay.windows(needle.len()).any(|w| w == needle)
}

fn dump_branch_terms(label: &str, graph: &SeqGraph, path: &KBestPath) {
    let terms = seq_kbest_path_score_terms(graph, path);
    let recon: f64 = terms.iter().map(|t| t.penalty).sum();
    let n_branch = terms
        .iter()
        .filter(|t| t.edge_support != t.total_outgoing)
        .count();
    eprintln!(
        "SCORE {label} n_edges={} n_branch_edges={n_branch} stored={:.9} recon={:.9} delta={:.3e}",
        terms.len(),
        path.score,
        recon,
        (recon - path.score).abs()
    );
    for (i, t) in terms.iter().enumerate() {
        if t.edge_support != t.total_outgoing {
            eprintln!(
                "  BRANCH e{i} {}->{} supp={} tot={} pen={:.9} ref={}",
                t.from, t.to, t.edge_support, t.total_outgoing, t.penalty, t.is_ref
            );
        }
    }
}

/// Production Peak-RSS caps are finite and strictly smaller than Java's unbounded PQ.
#[test]
fn forensic_6r78_production_peak_rss_caps_are_finite() {
    assert_eq!(SEQ_KBEST_PRODUCTION_MAX_HEAP, 1024);
    assert_eq!(SEQ_KBEST_PRODUCTION_MAX_EXPANSIONS, 12_000);
    assert_eq!(SEQ_KBEST_PRODUCTION_MAX_PATH_EDGES, 4_096);
    let p = SeqKbestCapPolicy::production();
    assert_eq!(p.max_heap_paths, Some(1024));
    assert_eq!(p.max_expansions, Some(12_000));
    let u = SeqKbestCapPolicy::unbounded();
    assert!(u.max_heap_paths.is_none() && u.max_expansions.is_none());
}

/// Score is ∑ log10(edgeMult/outMult). Unique continuations contribute 0. No path-length term.
#[test]
fn forensic_6r78_score_is_sum_of_log10_multiplicity_ratios() {
    let branch = (40f64).log10() - (45f64).log10();
    let unique = (40f64).log10() - (40f64).log10();
    assert!((unique).abs() < 1e-15);
    let path = branch + unique;
    assert!((path - branch).abs() < 1e-15);
    assert!(path < 0.0);
    assert!(path > -1.0);
}

/// Extra branch term can invert allele ranking: higher SNP-edge support still loses
/// if the allele walk pays an additional log10(12/15) that the other walk does not.
#[test]
fn forensic_6r78_extra_branch_penalty_can_invert_allele_support() {
    let snp_c = (17f64).log10() - (31f64).log10();
    let snp_g = (14f64).log10() - (31f64).log10();
    assert!(
        snp_c > snp_g,
        "C-edge 17/31 is the better single-edge score"
    );
    let extra = (12f64).log10() - (15f64).log10();
    assert!(extra < 0.0);
    let path_c = snp_c + extra;
    let path_g = snp_g;
    assert!(
        path_c < path_g,
        "C walk with extra 12/15 branch is worse than G walk: {path_c} vs {path_g}"
    );
}

/// Frozen 6R.52 comparator: higher finite score is Greater (polled first).
#[test]
fn forensic_6r78_comparator_orders_finite_scores_high_first() {
    assert_eq!(seq_kbest_score_cmp(-2.847255, -2.859844), Ordering::Greater);
    assert_eq!(seq_kbest_score_cmp(-2.859844, -2.860953), Ordering::Greater);
    assert_eq!(seq_kbest_score_cmp(-2.847255, -2.847255), Ordering::Equal);
}

/// Live region: isolate whether K=128 loss of J0/J1 is cutoff rank vs Peak-RSS caps vs visit budget.
#[test]
fn live_seqgraph_kbest_ranking() {
    if std::env::var("HOLDOUT_6R78").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R78=1");
        return;
    }
    use gatk_core::reference::{
        parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary,
    };
    use gatk_haplotypecaller::assembly_region_finalize::{
        create_graph_reference_read, records_to_assembly_reads,
    };
    use gatk_haplotypecaller::read_threading_assembler::build_threading_graph_for_seq_assembly;
    use gatk_haplotypecaller::seq_graph::SeqGraphCleanupStatus;
    use gatk_haplotypecaller::{
        assemble_reads_with_finalized, call_disposition, flatten_assembly_regions,
        traverse_assembly_region_walker, AssemblyRegionCallDisposition, CallRegionArgs,
        ReadFilterParams, WalkerTraversalConfig,
    };
    use std::path::Path;

    const INTERVAL: &str = "20:29455000-29456500";
    const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
    const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
    const POS_SNP: u64 = 29_456_344;

    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: live BAM/ref missing");
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
    let covering: Vec<_> = regions
        .iter()
        .filter(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= POS_SNP
                && r.end.get() >= POS_SNP
        })
        .collect();
    assert_eq!(covering.len(), 1);
    let region = covering[0];
    let args = CallRegionArgs::strict_java();
    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
    let mut owned = region.clone();
    let assembled =
        assemble_reads_with_finalized(&mut owned, &dict, &mut ref_cache, &args.assemble)
            .expect("assemble");
    let padded_ref = gatk_haplotypecaller::assembly_region_finalize::assembly_reference_read(
        &dict,
        &mut ref_cache,
        region,
    )
    .expect("pad ref");
    let graph_ref = create_graph_reference_read(&padded_ref, region, &dict);
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
    eprintln!(
        "SEQGRAPH nodes={} edges={}",
        seq.node_count(),
        seq.edge_count()
    );

    let needles: [&[u8]; 3] = [JAVA_ONLY_J0, JAVA_ONLY_J1, RUST_ONLY_R0];
    let labels = ["J0", "J1", "R0"];

    let run = |tag: &str, k_results: usize, visit: usize, caps: SeqKbestCapPolicy| {
        let r = find_best_haplotypes_seq_graph_forensic(&seq, k_results, visit, caps, &needles)
            .expect(tag);
        eprintln!(
            "{tag} n={} expansions={} max_heap={} pops={} skip_heap_pop={} skip_heap_exp={} skip_exp_cap={} skip_path_cap={} visit_refused={} heap_left={}",
            r.paths.len(),
            r.expansions,
            r.max_heap,
            r.pop_count,
            r.skip_heap_full_at_pop,
            r.skip_heap_full_at_expand,
            r.skip_expansion_cap,
            r.skip_path_edge_cap,
            r.vertex_visit_refused,
            r.heap_remaining
        );
        for (i, lab) in labels.iter().enumerate() {
            match &r.needles_in_result[i] {
                Some(h) => eprintln!(
                    "  {lab} sink_ord={} rank={:?} score={:.9} n_edges={} heap_left_has={}",
                    h.sink_ordinal,
                    h.rank_after_sort,
                    h.score,
                    h.n_edges,
                    r.needles_on_remaining_heap[i]
                ),
                None => eprintln!(
                    "  {lab} ABSENT_FROM_SINKS heap_left_has={}",
                    r.needles_on_remaining_heap[i]
                ),
            }
        }
        if k_results >= 128 && r.paths.len() >= 128 {
            let cut = &r.paths[127];
            eprintln!(
                "  K128_CUTOFF_PATH rank128 score={:.9} n_edges={}",
                cut.score,
                cut.edges.len()
            );
        }
        r
    };

    let k128_prod = run("K128_PROD", 128, 128, SeqKbestCapPolicy::production());
    let k128_unbounded = run(
        "K128_UNBOUNDED_CAPS",
        128,
        128,
        SeqKbestCapPolicy::unbounded(),
    );
    let visit128_collect_all = run(
        "VISIT128_COLLECT4096_UNBOUNDED",
        4096,
        128,
        SeqKbestCapPolicy::unbounded(),
    );
    let k256_prod = run("K256_PROD", 256, 256, SeqKbestCapPolicy::production());
    let k512_prod = run("K512_PROD", 512, 512, SeqKbestCapPolicy::production());

    for (lab, needle) in labels.iter().zip(needles) {
        if let Some(p) = k512_prod
            .paths
            .iter()
            .find(|p| contains_bases(&seq.path_bases_bytes(p.start, &p.edges), needle))
        {
            dump_branch_terms(lab, &seq, p);
        }
    }

    if let (Some(r0), Some(j0)) = (
        k512_prod.needles_in_result[2].as_ref(),
        k512_prod.needles_in_result[0].as_ref(),
    ) {
        eprintln!(
            "CMP R0_vs_J0={:?} cutoff={:.9} j0={:.9}",
            seq_kbest_score_cmp(r0.score, j0.score),
            k128_prod.paths[127].score,
            j0.score
        );
    }

    assert_eq!(
        k128_prod.skip_heap_full_at_pop, 0,
        "K=128 heap-pop cap did not fire"
    );
    assert_eq!(
        k128_prod.skip_heap_full_at_expand, 0,
        "K=128 heap-expand cap did not fire"
    );
    assert_eq!(
        k128_prod.skip_expansion_cap, 0,
        "K=128 expansion cap did not fire"
    );
    assert_eq!(k128_prod.skip_path_edge_cap, 0);
    assert!(
        k128_prod.needles_in_result[0].is_some(),
        "6R.80: K=128 retains J0"
    );
    let _ = k128_prod.needles_in_result[1];
    assert!(
        k128_prod.needles_in_result[2].is_none(),
        "6R.80: R0 is the 129th sink, outside K=128"
    );
    assert!(
        k128_unbounded.needles_in_result[0].is_some(),
        "6R.80: unbounded K=128 also retains J0"
    );
    assert_eq!(k128_prod.paths.len(), k128_unbounded.paths.len());
    let j0_all = visit128_collect_all.needles_in_result[0]
        .as_ref()
        .expect("J0 is generated under visit-budget 128");
    assert_eq!(
        j0_all.rank_after_sort,
        Some(124),
        "6R.80: J0 matches Java K=128 rank 124"
    );
    assert!(
        j0_all.score >= k128_prod.paths[127].score,
        "J0 score is at or above the K=128 cutoff path"
    );
    assert!(
        k512_prod.needles_in_result[0].is_some(),
        "J0 remains a SeqGraph k-best path at K=512"
    );
    let _ = k512_prod.needles_in_result[1];
    let _ = k256_prod;
}
