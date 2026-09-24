//! 6R.239: Java-compatible dangling-tail `refIndexToMerge` sentinel.
//! Production change: replace saturating unsigned subtraction with checked
//! arithmetic so path-index 0 remains the LCA / no-splice sentinel.
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! Do not raise K.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r239_dangling_index_sentinel_java_semantics -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::assembly::{AssemblyGraph, AssemblyGraphPruningParams};
use gatk_haplotypecaller::assembly_based_caller::assemble_reads_with_finalized;
use gatk_haplotypecaller::assembly_dangling_recovery::{
    dangling_tail_ref_index_to_merge, DanglingRecoveryParams,
};
use gatk_haplotypecaller::assembly_region_finalize::{
    assembly_reference_read, create_graph_reference_read, records_to_assembly_reads,
};
use gatk_haplotypecaller::hc_genotyping_engine::java_alignment_read_overlaps_interval;
use gatk_haplotypecaller::read_threading_assembler::{
    build_threading_graph_for_seq_assembly, DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH,
};
use gatk_haplotypecaller::seq_graph::SeqGraph;
use gatk_haplotypecaller::seq_kbest_haplotype::{
    find_best_haplotypes_seq_graph, seq_kbest_path_score_terms,
};
use gatk_haplotypecaller::{
    assembly_graph_from_ref_and_reads_threading, call_disposition, flatten_assembly_regions,
    traverse_assembly_region_walker, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, KmerSize, ReadFilterParams, WalkerTraversalConfig,
};
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
const JAVA_CARRIER_EDGES: usize = 18;
const K971: &[u8] = b"GAATTCATTTGGGTGATTTTTTTGT";
const K428: &[u8] = b"GGAATTCATTTGGGTGATTTTTTTT";
const SINK_38: &str = "GTTGGTGGAATTCATTTGGGTGATTTTTTTTGTTTCTT";
const SINK_8: &str = "TGTTTCTT";

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R239\t{key}\t{}", value.as_ref());
}

fn fnv1a64_hex(data: &[u8]) -> String {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in data {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    format!("{h:016x}")
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

fn rt_in_degree(g: &AssemblyGraph, node: usize) -> usize {
    g.edges_sorted()
        .into_iter()
        .filter(|e| e.to == node)
        .count()
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

fn merge_decision(last: usize, suf: usize, lead: bool) -> Result<usize, &'static str> {
    let idx = dangling_tail_ref_index_to_merge(last, suf, lead)?;
    if idx == 0 {
        Err("ref_index_zero_cycle")
    } else {
        Ok(idx)
    }
}

#[test]
fn forensic_6r239_dangling_index_sentinel_java_semantics() {
    kv("java_pin", JAVA_PIN);
    kv(
        "production_change",
        "plan_dangling_tail_merge: saturating_sub → checked (last+1).checked_sub(matching_suffix)",
    );
    kv("classification", "INDEX_SENTINEL_SEMANTICS_DIVERGENCE");
    assert_eq!(DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, 128);

    let rows: &[(usize, usize, bool, Result<usize, &str>)] = &[
        (8, 9, false, Err("ref_index_zero_cycle")),
        (9, 9, false, Ok(1)),
        (10, 9, false, Ok(2)),
        (8, 8, false, Ok(1)),
        (0, 1, false, Err("ref_index_zero_cycle")),
        (8, 9, true, Ok(1)),
        (8, 8, true, Ok(2)),
        (0, 1, true, Ok(1)),
    ];
    kv(
        "boundary_header",
        "lastRef\tmatchingSuffix\tleadingDel\texpected",
    );
    for &(last, suf, lead, expect) in rows {
        let got = merge_decision(last, suf, lead);
        assert_eq!(got, expect, "last={last} suf={suf} lead={lead}");
        kv("boundary_row", format!("{last}\t{suf}\t{lead}\t{expect:?}"));
    }
    assert_eq!(
        dangling_tail_ref_index_to_merge(0, 2, false),
        Err("ref_index_to_merge_underflow")
    );
    kv(
        "underflow_0_2",
        "unreachable from longestSuffixMatch; must not become path index 0 or 1",
    );

    let root = repo_root();
    let bam = root.join(BAM_REL);
    let ref_fasta = root.join(REF_REL);
    assert!(bam.is_file() && ref_fasta.is_file());
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
    let region = flatten_assembly_regions(&walked)
        .into_iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= TARGET
                && r.end.get() >= TARGET
        })
        .expect("ActiveFull");
    let args = CallRegionArgs::strict_java();
    let mut owned = region.clone();
    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
    let assembled =
        assemble_reads_with_finalized(&mut owned, &dict, &mut ref_cache, &args.assemble)
            .expect("assemble");
    let padded = assembly_reference_read(&dict, &mut ref_cache, &region).expect("pad");
    let graph_ref = create_graph_reference_read(&padded, &region, &dict);
    let graph_reads = records_to_assembly_reads(&assembled.finalized_reads);
    let assembler = args.assemble.assembler.clone();
    let params = threading_params(25, &assembler);
    let raw =
        assembly_graph_from_ref_and_reads_threading(&graph_ref, &graph_reads, &params).expect("rt");
    let mut pruned = raw;
    let mut pruning = AssemblyGraphPruningParams::gatk_haplotype_caller_defaults();
    pruning.min_prune_factor = assembler.min_prune_factor;
    pruning.use_adaptive_pruning = assembler.use_adaptive_pruning;
    let _ = pruned.apply_pruning(&pruning);
    let dangling = DanglingRecoveryParams::from_assembler_args(&assembler);
    assert!(dangling.dangling_java_exact);

    assert!(rt_edge(&pruned, K971, K428).is_none());
    let id428 = rt_find_kmer(&pruned, K428).expect("k428");
    kv(
        "prune_428_in_degree",
        rt_in_degree(&pruned, id428).to_string(),
    );

    let mut saw_canonical = false;
    for (sink, _, _) in pruned.probe_dangling_tail_failures(&dangling) {
        let d = pruned.dangling_tail_decision_dump(sink, &dangling);
        if d.cigar == "3I9M" && d.matching_suffix == 9 && d.last_ref_idx == 8 {
            saw_canonical = true;
            kv(
                "canonical_dump",
                format!(
                    "sink={} ref_idx={:?} reject={:?} from={:?} to={:?}",
                    d.sink, d.ref_index_to_merge, d.reject_reason, d.from, d.to
                ),
            );
            assert_eq!(d.ref_index_to_merge, Some(0));
            assert_eq!(d.reject_reason, Some("ref_index_zero_cycle"));
            assert!(d.from.is_none());
            assert!(d.to.is_none());
            assert_eq!(
                merge_decision(d.last_ref_idx, d.matching_suffix, false),
                Err("ref_index_zero_cycle")
            );
        }
    }
    assert!(saw_canonical, "3I9M lastRef=8 matchSuf=9 must still exist");

    let mut dangled = pruned.clone();
    let _ = dangled
        .recover_dangling_branches(&dangling)
        .expect("dangling");
    let splice = rt_edge(&dangled, K971, K428);
    kv("after_dangling_971_428", format!("{splice:?}"));
    assert!(splice.is_none(), "Java sentinel: no addEdge(971,428)");
    let id428 = rt_find_kmer(&dangled, K428).expect("k428 after dangling");
    let in_deg = rt_in_degree(&dangled, id428);
    kv("after_dangling_428_in_degree", in_deg.to_string());
    assert_eq!(in_deg, 1, "vertex 428 must not gain the Rust-only incoming");

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
    let mut seq = SeqGraph::from_assembly_graph(&seq_asm);
    seq.clean_non_ref_paths();
    let _ = seq.cleanup_seq_graph();
    let sink_seqs: Vec<String> = seq
        .vertices()
        .iter()
        .filter(|v| !seq.edges().iter().any(|e| e.from == v.id))
        .map(|v| String::from_utf8_lossy(&v.sequence).into_owned())
        .collect();
    kv("seq_sinks", sink_seqs.join("|"));
    assert!(
        sink_seqs.iter().any(|s| s == SINK_38),
        "Java-equivalent 38 bp sink must appear after zip"
    );
    assert!(
        !sink_seqs.iter().any(|s| s == SINK_8),
        "8 bp TGTTTCTT sink is the zip-blocked artifact of the splice"
    );

    let off = (EVENT_POS - region.extended_start.get()) as usize;
    let paths128 = find_best_haplotypes_seq_graph(&seq, 128).expect("k128");
    let k128_idx = paths128
        .iter()
        .position(|p| fnv1a64_hex(&seq.path_bases_bytes(p.start, &p.edges)) == JAVA_CARRIER);
    kv("k128_carrier_rank", format!("{k128_idx:?}"));
    let k128_idx = k128_idx.expect("K=128 must include carrier after sentinel fix");
    assert_eq!(k128_idx, JAVA_CARRIER_RANK);
    assert_eq!(DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, 128);

    let paths256 = find_best_haplotypes_seq_graph(&seq, 256).expect("k256");
    let k256_idx = paths256
        .iter()
        .position(|p| fnv1a64_hex(&seq.path_bases_bytes(p.start, &p.edges)) == JAVA_CARRIER)
        .expect("K=256 carrier");
    let cpath = &paths256[k256_idx];
    let terms = seq_kbest_path_score_terms(&seq, cpath);
    kv(
        "k256_carrier",
        format!(
            "rank={} score={:.8} n_edges={} k128_rank={}",
            k256_idx,
            cpath.score,
            terms.len(),
            k128_idx
        ),
    );
    assert_eq!(k256_idx, JAVA_CARRIER_RANK);
    assert!((cpath.score - JAVA_CARRIER_SCORE).abs() < 1e-7);
    assert_eq!(terms.len(), JAVA_CARRIER_EDGES);
    let cbases = seq.path_bases_bytes(cpath.start, &cpath.edges);
    assert_eq!(cbases.get(off).copied(), Some(b'C'));
    kv("k_unchanged", "K=128 production cap not modified");

    let outcome = HaplotypeCallerEngine::call_region(
        &region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let ga = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(TARGET)
                && c.event.ref_allele == "G"
                && c.event.alt_allele == "A"
        })
        .expect("G/A");
    kv("format_ad", format!("{:?}", ga.genotype.format.ad_as_i32()));
    kv("format_pl", format!("{:?}", ga.genotype.format.pl_as_i32()));
    kv("format_dp", format!("{:?}", ga.genotype.format.dp.as_i32()));
    let pairhmm_idx: std::collections::BTreeSet<usize> = outcome
        .read_likelihoods
        .iter()
        .map(|c| c.read_index.get())
        .collect();
    let rust_n = outcome
        .genotyping_reads
        .iter()
        .enumerate()
        .filter(|(i, r)| {
            pairhmm_idx.contains(i) && java_alignment_read_overlaps_interval(r, TARGET, TARGET, 2)
        })
        .count();
    kv("n_likelihoods", outcome.read_likelihoods.len().to_string());
    kv("rust_overlap_membership", rust_n.to_string());
    assert_eq!(ga.genotype.format.ad_as_i32(), vec![42, 5]);
    assert_eq!(ga.genotype.format.pl_as_i32(), vec![84, 0, 1738]);
    assert_eq!(ga.genotype.format.dp.as_i32(), 47);
    assert_eq!(rust_n, 47);
}
