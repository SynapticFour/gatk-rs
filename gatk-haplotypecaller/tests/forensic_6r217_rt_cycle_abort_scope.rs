//! 6R.217: diagnostic Java `createGraph` cycle-abort scope for post-SeqGraph
//! RT extract at `20:29455015 G/T`.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! Proof-only. PRODUCTION CHANGE: NONE.
//!
//! Java `createGraph`: prune, then `if (generateSeqGraph && hasCycles()) return
//! null` **before** dangling recovery. Rust `extract_rt_haplotypes_before_remove_paths`
//! calls `build_threading_graph_core(..., abort_cyclic_before_dangling=false)`,
//! so a cyclic graph still reaches dangling recovery and k-best.
//!
//! Diagnostic (test-only): skip RT extract when the SeqGraph `createGraph`
//! builder returns `None`. Do not change production merge_rt / extract.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r217_rt_cycle_abort_scope -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::assembly::{AssemblyGraphParams, AssemblyGraphPruningParams};
use gatk_haplotypecaller::assembly_region_finalize::{
    assembly_reference_read, create_graph_reference_read, records_to_assembly_reads,
};
use gatk_haplotypecaller::bio_ids::KmerSize;
use gatk_haplotypecaller::read_threading_assembler::{
    assemble_from_ref_and_reads, build_threading_graph_for_haplotype_dump,
    build_threading_graph_for_seq_assembly, extract_haplotypes_from_seq_kbest_paths,
    extract_rt_haplotypes_before_remove_paths, merge_rt_kbest_pre_remove_paths,
};
use gatk_haplotypecaller::read_threading_graph::assembly_graph_from_ref_and_reads_threading_with_summary;
use gatk_haplotypecaller::seq_graph::SeqGraph;
use gatk_haplotypecaller::seq_kbest_haplotype::find_best_haplotypes_seq_graph;
use gatk_haplotypecaller::{
    assemble_reads_with_finalized, call_disposition, flatten_assembly_regions,
    traverse_assembly_region_walker, AssemblyRead, AssemblyRegionCallDisposition, CallRegionArgs,
    Cigar, Haplotype, ReadFilterParams, ReadThreadingAssemblerArgs, WalkerTraversalConfig,
};
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const P11_REF_REL: &str = "parity/fixtures/p5_live_reference.fa";
const P11_BAM_REL: &str = "parity/build/sam-indexed-bam/p11_java_positive.bam";
const INDEL4_REF_REL: &str = "parity/fixtures/g2_subset_live_indel4.fa";
const INDEL4_BAM_REL: &str = "parity/build/sam-indexed-bam/g2_subset_live_indel4.bam";
const P12_BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const TARGET: u64 = 29_455_015;
const JAVA_ACTIVE: &str = "20:29455000-29455145";
const JAVA_PADDED: &str = "20:29454900-29455245";
const JAVA_UNTRIMMED_N: usize = 78;
const P12_TTC: u64 = 92_307_324;
const P12_TG: u64 = 92_307_333;
const TG_INTERVAL: &str = "2:92307200-92307550";
const MID_B_INTERVAL: &str = "2:92317000-92319000";
const CLOSED_MID_B: u64 = 92_317_399;

const JAVA_UNTRIMMED_HASHES: &[&str] = &[
    "72369f70c97564fb",
    "9b4130f578a5c3b1",
    "0760164d9ef6dd05",
    "a411db2e5639776a",
    "c3ca8558d183ddd0",
    "ad1727667b42fe9e",
    "c599fb57c27f5a17",
    "31f12d6e2e7365a1",
    "b88de6b98d11d0d7",
    "ed01a60a7f2d5970",
    "b55bcced2fefd423",
    "bfe3c2fb2c2a240d",
    "3990d39c635fde1a",
    "9177506879ad0a9c",
    "810d9fc18079ce30",
    "44d25f40711fe1fa",
    "54afee1dd89e373b",
    "649c345e483e5443",
    "53cb3d3d6311786a",
    "a2f24e4d4a6e8d9e",
    "b9c46fd1377b9041",
    "cdef5e7202584ed0",
    "d1be531684070ff1",
    "eeaca65161357903",
    "5c88c58d266f75d4",
    "cc3f9cadeaeeb4c7",
    "833d4f233780f9a8",
    "9825a1930286fe3e",
    "4afa9304857dc4b1",
    "4685cd1bf0436e72",
    "2cf89addc6b88d5c",
    "8bea11740b6e4204",
    "3bfc7896af17a76f",
    "a14f1b575e96b1fe",
    "e114fb78453ccc5b",
    "3bfbd484c5d0ed45",
    "8aca083d6ff2a5cd",
    "40e4247f4b1a91ec",
    "71e588c6222f0401",
    "8e95b473d340af19",
    "cceb12a3a0c7b29e",
    "20bc439d334fc9ff",
    "de388cc436a18b40",
    "1c39b8c66923c1d2",
    "3fc7de0954270452",
    "c48d5100bcad6748",
    "ae7714c59dafc75e",
    "188da0dd30b021c1",
    "793815145ccc35bb",
    "391503e0ef3650da",
    "dfb7ae3851aed80d",
    "15f1bf3597e17af4",
    "0994a4aaad72cff3",
    "b10c816794fadac2",
    "97cf358727e398f9",
    "284791f3257fae56",
    "93b1d9eda04c9065",
    "985a04a85c946f9a",
    "6ed6c60d101efa23",
    "a7c57619bbaef828",
    "7f563f731e72eaf7",
    "23afacd2bee260a7",
    "edc349af27dd02d0",
    "ac3b46130edb998e",
    "b7ac973e561ef465",
    "8016ead01172b854",
    "9a86c939b608f767",
    "8dea0ce0dc47c9fc",
    "6083fb7b49feebc5",
    "fe10cb506100d9cb",
    "add30f91e861475a",
    "3a3a91f45eb372c8",
    "fdfe44dbde6f93d1",
    "1dc3e207f92421d8",
    "e6a394771b59fa81",
    "2a2b6e4d256f5695",
    "550fe41efa3bd69c",
    "7a4a7948ede22768",
];

const REPRESENTATIVE_EXTRAS: &[&str] = &[
    "70825a603033cd15",
    "20734c56283366c5",
    "e0c1092ca29b8eb6",
    "00bda5d67ed6ce26",
];

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R217\t{key}\t{}", value.as_ref());
}

fn fnv1a64_hex(bases: &[u8]) -> String {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in bases {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    format!("{h:016x}")
}

fn hap_hashes(haps: &[Haplotype]) -> BTreeSet<String> {
    haps.iter().map(|h| fnv1a64_hex(&h.bases)).collect()
}

fn unique_insert(haps: &mut Vec<Haplotype>, seen: &mut BTreeSet<String>, h: Haplotype) {
    if seen.insert(fnv1a64_hex(&h.bases)) {
        haps.push(h);
    }
}

/// Java `createGraph` cycle test: prune, then `hasCycles()`, before dangling.
fn cyclic_before_dangling(
    reference: &AssemblyRead,
    reads: &[AssemblyRead],
    assembler: &ReadThreadingAssemblerArgs,
    kmer: usize,
) -> Option<bool> {
    let kmer_size = KmerSize::try_from_usize(kmer).ok()?;
    let params = AssemblyGraphParams {
        kmer_size,
        min_base_quality: assembler.min_base_quality,
        min_edge_weight: 1,
        dangling_path_max_nodes: 0,
        max_haplotypes: assembler.num_best_haplotypes_per_graph,
        max_haplotype_bases: 4096,
        start_threading_only_at_existing_vertex: !assembler.recover_dangling_branches,
    };
    let (mut graph, summary) =
        assembly_graph_from_ref_and_reads_threading_with_summary(reference, reads, &params).ok()?;
    if !assembler.allow_low_complexity_graphs && summary.is_low_complexity {
        return None;
    }
    let mut pruning = AssemblyGraphPruningParams::gatk_haplotype_caller_defaults();
    pruning.min_prune_factor = assembler.min_prune_factor;
    pruning.use_adaptive_pruning = assembler.use_adaptive_pruning;
    if assembler.prune_before_cycle_counting {
        graph.apply_pruning(&pruning);
    }
    Some(graph.has_cycle())
}

fn seq_kbest_haps(
    graph_ref: &AssemblyRead,
    graph_reads: &[AssemblyRead],
    assembler: &ReadThreadingAssemblerArgs,
    kmer: usize,
) -> Vec<Haplotype> {
    let Some(graph) = build_threading_graph_for_seq_assembly(
        graph_ref,
        graph_reads,
        kmer,
        assembler,
        false,
        false,
    )
    .ok()
    .flatten() else {
        return Vec::new();
    };
    let mut seq = SeqGraph::from_assembly_graph(&graph);
    seq.clean_non_ref_paths();
    let _ = seq.cleanup_seq_graph();
    let Ok(paths) = find_best_haplotypes_seq_graph(&seq, assembler.num_best_haplotypes_per_graph)
    else {
        return Vec::new();
    };
    let mut ref_hap = Haplotype::new(graph_ref.bases.as_slice(), true);
    let mut ref_cigar = Cigar::new();
    ref_cigar.push(
        ref_hap.bases.len(),
        gatk_haplotypecaller::CigarOperator::Match,
    );
    ref_hap.cigar = Some(ref_cigar);
    let ref_cigar_len = ref_hap.cigar.as_ref().unwrap().reference_length();
    extract_haplotypes_from_seq_kbest_paths(
        &paths,
        &seq,
        kmer,
        &ref_hap,
        ref_cigar_len,
        &assembler.haplotype_to_reference_sw,
    )
    .unwrap_or_default()
}

/// Test-only: Java `createGraph(generateSeqGraph=true)` → None means no RT extract.
fn diagnostic_extract_honoring_create_graph(
    graph_ref: &AssemblyRead,
    graph_reads: &[AssemblyRead],
    assembler: &ReadThreadingAssemblerArgs,
    kmer: usize,
) -> Vec<Haplotype> {
    let seq_g = build_threading_graph_for_seq_assembly(
        graph_ref,
        graph_reads,
        kmer,
        assembler,
        false,
        false,
    )
    .ok()
    .flatten();
    if seq_g.is_none() {
        return Vec::new();
    }
    extract_rt_haplotypes_before_remove_paths(graph_ref, graph_reads, assembler, kmer, false, false)
        .unwrap_or_default()
}

fn diagnostic_merge(
    seq_haps: &[Haplotype],
    graph_ref: &AssemblyRead,
    graph_reads: &[AssemblyRead],
    assembler: &ReadThreadingAssemblerArgs,
    kmers: &[usize],
) -> Vec<Haplotype> {
    let mut haps = seq_haps.to_vec();
    let mut seen = hap_hashes(&haps);
    for &k in kmers {
        let batch = diagnostic_extract_honoring_create_graph(graph_ref, graph_reads, assembler, k);
        for h in batch {
            unique_insert(&mut haps, &mut seen, h);
        }
        if haps.iter().any(|h| !h.is_reference) && haps.len() > 1 {
            break;
        }
    }
    haps
}

#[test]
fn forensic_6r217_caller_matrix_pin() {
    kv("java_pin", JAVA_PIN);
    kv(
        "extract_abort_flag",
        "extract_rt_haplotypes_from_built_graph always passes abort_cyclic_before_dangling=false",
    );
    kv(
        "seq_assembly_abort_flag",
        "build_threading_graph_for_seq_assembly passes true (Java createGraph)",
    );
    kv(
        "dump_abort_flag",
        "build_threading_graph_for_haplotype_dump passes false (parity dumps / RT fallback)",
    );
    kv(
        "caller_merge_rt",
        "merge_rt_kbest_pre_remove_paths → extract_rt_before_remove; NO Java counterpart; cycle abort OFF; candidate YES",
    );
    kv(
        "caller_supplement_p12",
        "supplement_p12_cluster_coupled_haplotypes → extract_rt_before_remove; P12/strict scoring; cycle abort OFF; candidate YES (same extract)",
    );
    kv(
        "caller_seq_createGraph",
        "assemble_from_ref_and_reads_seq_graph → build_threading_graph_for_seq_assembly; Java createGraph; cycle abort ON; candidate N/A (already correct)",
    );
    kv(
        "caller_rt_fallback",
        "assemble_from_ref_and_reads_rt_graph / use_seq_graph=false zip-collapse; Java generateSeqGraph=false does not abort cycles; candidate NO",
    );
    kv(
        "flag_reason",
        "abort=false is the haplotype-dump builder so cyclic graphs still k-best; L2 comments about before_remove vs after_remove are path-connectivity, not cycle abort",
    );
    assert!(ReadThreadingAssemblerArgs::default().use_seq_graph);
    assert!(ReadThreadingAssemblerArgs::default().abort_seq_graph_on_cycles);
    assert_eq!(JAVA_UNTRIMMED_HASHES.len(), JAVA_UNTRIMMED_N);
}

#[test]
fn forensic_6r217_live_cycle_abort_diagnostic() {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing chr20_tiny BAM/REF");
        return;
    }
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE — proof-only");

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
        .expect("covering");
    assert_eq!(
        format!(
            "{}:{}-{}",
            region.contig,
            region.start.get(),
            region.end.get()
        ),
        JAVA_ACTIVE
    );
    assert_eq!(
        format!(
            "{}:{}-{}",
            region.contig,
            region.extended_start.get(),
            region.extended_end.get()
        ),
        JAVA_PADDED
    );

    let mut owned = region.clone();
    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
    let args = CallRegionArgs::strict_java();
    let assembled =
        assemble_reads_with_finalized(&mut owned, &dict, &mut ref_cache, &args.assemble)
            .expect("assemble");
    let padded_ref = assembly_reference_read(&dict, &mut ref_cache, region).expect("pad ref");
    let graph_ref = create_graph_reference_read(&padded_ref, region, &dict);
    let graph_reads = records_to_assembly_reads(&assembled.finalized_reads);
    let assembler = args.assemble.assembler.clone();
    assert!(assembler.use_seq_graph);
    assert!(assembler.abort_seq_graph_on_cycles);
    assert!(assembler.dangling_java_exact);

    let cyclic10 = cyclic_before_dangling(&graph_ref, &graph_reads, &assembler, 10);
    let cyclic25 = cyclic_before_dangling(&graph_ref, &graph_reads, &assembler, 25);
    let seq_g10 = build_threading_graph_for_seq_assembly(
        &graph_ref,
        &graph_reads,
        10,
        &assembler,
        false,
        false,
    )
    .expect("seq k10");
    let seq_g25 = build_threading_graph_for_seq_assembly(
        &graph_ref,
        &graph_reads,
        25,
        &assembler,
        false,
        false,
    )
    .expect("seq k25");
    let dump_g10 = build_threading_graph_for_haplotype_dump(
        &graph_ref,
        &graph_reads,
        10,
        &assembler,
        false,
        false,
    )
    .expect("dump k10");
    kv("k10_cyclic_before_dangling", format!("{cyclic10:?}"));
    kv("k25_cyclic_before_dangling", format!("{cyclic25:?}"));
    kv(
        "java_createGraph_k10",
        "null (generateSeqGraph && hasCycles before dangling)",
    );
    kv(
        "rust_seq_createGraph_k10",
        if seq_g10.is_none() { "None" } else { "Some" },
    );
    kv(
        "rust_dump_graph_k10_after_dangling",
        if dump_g10.is_some() { "Some" } else { "None" },
    );
    kv(
        "first_divergent_op",
        "cycle detection result ignored in extract_rt (abort_cyclic_before_dangling=false); dangling recovery is downstream, not the first split",
    );
    assert_eq!(cyclic10, Some(true));
    assert_eq!(cyclic25, Some(false));
    assert!(
        seq_g10.is_none(),
        "Java createGraph k=10 is None, not empty k-best"
    );
    assert!(seq_g25.is_some());
    assert!(
        dump_g10.is_some(),
        "dump/extract builder still materializes the cyclic graph then dangling/k-best"
    );

    let java_untrimmed: BTreeSet<&str> = JAVA_UNTRIMMED_HASHES.iter().copied().collect();
    let seq25 = seq_kbest_haps(&graph_ref, &graph_reads, &assembler, 25);
    let seq10 = seq_kbest_haps(&graph_ref, &graph_reads, &assembler, 10);
    let seq25_h = hap_hashes(&seq25);
    let seq25_refs: BTreeSet<&str> = seq25_h.iter().map(String::as_str).collect();
    kv("seq_kbest_k10_n", seq10.len().to_string());
    kv("seq_kbest_k25_n", seq25.len().to_string());
    assert_eq!(seq10.len(), 0);
    assert_eq!(seq25.len(), JAVA_UNTRIMMED_N);
    assert_eq!(java_untrimmed.difference(&seq25_refs).count(), 0);
    assert_eq!(seq25_refs.difference(&java_untrimmed).count(), 0);

    let rt10_seqgraph = extract_rt_haplotypes_before_remove_paths(
        &graph_ref,
        &graph_reads,
        &assembler,
        10,
        false,
        false,
    )
    .unwrap_or_default();
    let mut dump_args = assembler.clone();
    dump_args.use_seq_graph = false;
    let rt10 = extract_rt_haplotypes_before_remove_paths(
        &graph_ref,
        &graph_reads,
        &dump_args,
        10,
        false,
        false,
    )
    .unwrap_or_default();
    let rt25 = extract_rt_haplotypes_before_remove_paths(
        &graph_ref,
        &graph_reads,
        &assembler,
        25,
        false,
        false,
    )
    .unwrap_or_default();
    let rt10_h = hap_hashes(&rt10);
    let diag10 = diagnostic_extract_honoring_create_graph(&graph_ref, &graph_reads, &assembler, 10);
    let diag25 = diagnostic_extract_honoring_create_graph(&graph_ref, &graph_reads, &assembler, 25);
    kv("rt_k10_seqgraph_extract_n", rt10_seqgraph.len().to_string());
    kv("rt_k10_fallback_extract_n", rt10.len().to_string());
    kv("rt_k10_after_cycle_abort_n", diag10.len().to_string());
    kv("rt_k25_before_n", rt25.len().to_string());
    kv("rt_k25_after_cycle_abort_n", diag25.len().to_string());
    assert!(
        rt10_seqgraph.is_empty(),
        "6R.218 SeqGraph-path extract aborts cyclic k=10"
    );
    assert_eq!(
        rt10.len(),
        54,
        "use_seq_graph=false extract still recovers the cyclic k=10 paths Java never materializes"
    );
    assert_eq!(diag10.len(), 0, "createGraph None, not empty findBestPaths");
    assert_eq!(diag25.len(), rt25.len());

    let assemble_h = hap_hashes(&assembled.assembly.haplotypes);
    let assemble_refs: BTreeSet<&str> = assemble_h.iter().map(String::as_str).collect();
    let rust_only: Vec<&str> = assemble_refs.difference(&java_untrimmed).copied().collect();
    kv("extras_before_n", rust_only.len().to_string());
    assert_eq!(
        rust_only.len(),
        0,
        "6R.218 production assemble no longer admits the 44 extras"
    );
    let extras_in_rt10 = rust_only.iter().filter(|h| rt10_h.contains(**h)).count();
    assert_eq!(extras_in_rt10, 0);
    for h in rust_only.iter().take(8) {
        let rec = rt10.iter().find(|x| fnv1a64_hex(&x.bases) == *h);
        if let Some(hap) = rec {
            kv(
                "extra_provenance",
                format!(
                    "hash={}\tlen={}\tscore={}\tkmer={}\tsource=cyclic_RT_k10_before_remove\tcyclic=true\tstage=extract_rt_abort_false\tjava_createGraph=null",
                    h,
                    hap.bases.len(),
                    hap.score,
                    hap.kmer_size
                ),
            );
        }
    }
    for pin in REPRESENTATIVE_EXTRAS {
        assert!(
            rt10_h.contains(*pin),
            "representative extra {pin} must still exist in dump extract of the cyclic graph"
        );
        assert!(
            !assemble_refs.contains(pin),
            "representative extra {pin} must be absent from production assemble after cycle abort, not hash-deleted from the dump graph"
        );
    }

    let mut merge_current = seq25.clone();
    merge_rt_kbest_pre_remove_paths(
        &graph_ref,
        &graph_reads,
        &assembler,
        &[10, 25],
        &mut merge_current,
    )
    .expect("merge current");
    let merge_current_h = hap_hashes(&merge_current);
    let merge_diag = diagnostic_merge(&seq25, &graph_ref, &graph_reads, &assembler, &[10, 25]);
    let merge_diag_h = hap_hashes(&merge_diag);
    let extras_after: Vec<&str> = rust_only
        .iter()
        .copied()
        .filter(|h| merge_diag_h.contains(*h))
        .collect();
    kv("merged_before_unique_n", merge_current_h.len().to_string());
    kv("merged_after_unique_n", merge_diag_h.len().to_string());
    kv("extras_after_n", extras_after.len().to_string());
    kv(
        "prod_assemble_n",
        assembled.assembly.haplotypes.len().to_string(),
    );
    let prod = assemble_from_ref_and_reads(&graph_ref, &graph_reads, &assembler).expect("prod");
    kv(
        "prod_assemble_from_ref_n",
        prod.haplotypes.len().to_string(),
    );

    let merge_then_supp_diag = {
        let mut haps = merge_diag.clone();
        let mut seen = hap_hashes(&haps);
        for &k in &[10usize, 25] {
            let batch =
                diagnostic_extract_honoring_create_graph(&graph_ref, &graph_reads, &assembler, k);
            for h in batch {
                unique_insert(&mut haps, &mut seen, h);
            }
            if haps.iter().any(|h| !h.is_reference) && haps.len() > 1 {
                break;
            }
        }
        haps
    };
    kv(
        "merge_plus_supplement_diag_unique_n",
        hap_hashes(&merge_then_supp_diag).len().to_string(),
    );

    let merge_diag_refs: BTreeSet<&str> = merge_diag_h.iter().map(String::as_str).collect();
    assert_eq!(extras_after.len(), 0);
    assert_eq!(merge_diag_h.len(), JAVA_UNTRIMMED_N);
    assert_eq!(java_untrimmed.difference(&merge_diag_refs).count(), 0);
    assert_eq!(hap_hashes(&merge_then_supp_diag).len(), JAVA_UNTRIMMED_N);
    assert_eq!(
        prod.haplotypes.len(),
        JAVA_UNTRIMMED_N,
        "6R.218 production extract now aborts cyclic k=10"
    );
    kv("classification", "GRAPH_STATE_DIVERGENCE");
    kv(
        "first_causal_arrow",
        "6R.218: SeqGraph-path extract_rt abort_cyclic_before_dangling=true → k=10 contributes 0 → extras gone → unique merge = 78",
    );
    kv(
        "why_not_lifecycle_alone",
        "merge_rt still runs; k=10 empty + SeqGraph already has alts → early-stop skips k=25 RT; extras vanish from cycle abort, not from deleting merge_rt",
    );
    kv("kbest_policy", "legacy_1024 unchanged; K=128 unchanged");
}

#[test]
fn forensic_6r217_l2_p11_acyclic_alt_survives() {
    let root = repo_root();
    let ref_fasta = root.join(P11_REF_REL);
    let bam = root.join(P11_BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing p11 fixture");
        return;
    }
    kv("neg_control", "L2 p11_java_positive");
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, "chrLive:1-63").expect("interval");
    let walk = traverse_assembly_region_walker(
        &dict,
        &specs,
        &ref_fasta,
        &bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(0),
    )
    .expect("walk");
    let region = flatten_assembly_regions(&walk)
        .into_iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            )
        })
        .expect("p11 active");
    let mut owned = region.clone();
    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
    let args = CallRegionArgs::strict_java();
    let assembled =
        assemble_reads_with_finalized(&mut owned, &dict, &mut ref_cache, &args.assemble)
            .expect("assemble");
    let padded_ref = assembly_reference_read(&dict, &mut ref_cache, &region).expect("pad");
    let graph_ref = create_graph_reference_read(&padded_ref, &region, &dict);
    let graph_reads = records_to_assembly_reads(&assembled.finalized_reads);
    let assembler = args.assemble.assembler.clone();
    let cyclic10 = cyclic_before_dangling(&graph_ref, &graph_reads, &assembler, 10);
    let seq10 = seq_kbest_haps(&graph_ref, &graph_reads, &assembler, 10);
    let seq25 = seq_kbest_haps(&graph_ref, &graph_reads, &assembler, 25);
    let diag10 = diagnostic_extract_honoring_create_graph(&graph_ref, &graph_reads, &assembler, 10);
    let prod = assemble_from_ref_and_reads(&graph_ref, &graph_reads, &assembler).expect("prod");
    kv("p11_k10_cyclic", format!("{cyclic10:?}"));
    kv("p11_seq10_n", seq10.len().to_string());
    kv("p11_seq25_n", seq25.len().to_string());
    kv("p11_diag_rt10_n", diag10.len().to_string());
    kv("p11_assemble_n", prod.haplotypes.len().to_string());
    kv(
        "p11_alt_n",
        prod.haplotypes
            .iter()
            .filter(|h| !h.is_reference)
            .count()
            .to_string(),
    );
    assert_eq!(cyclic10, Some(false));
    assert_eq!(seq10.len(), 2);
    assert!(seq10.iter().any(|h| !h.is_reference));
    assert!(
        !diag10.is_empty(),
        "cycle-only abort must still extract acyclic k=10"
    );
    assert!(prod.haplotypes.iter().any(|h| !h.is_reference));
}

#[test]
fn forensic_6r217_p12_ttc_and_acyclic_rt_control() {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let p12_bam = root.join(P12_BAM_REL);
    if !ref_fasta.is_file() || !p12_bam.is_file() {
        eprintln!("skip: missing P12 BAM/REF");
        return;
    }
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let args = CallRegionArgs::strict_java();

    let probe = |interval: &str, pos: u64, label: &str| {
        let specs = parse_intervals_cli_string(&dict, interval).expect("interval");
        let walk = traverse_assembly_region_walker(
            &dict,
            &specs,
            &ref_fasta,
            &p12_bam,
            &ReadFilterParams::gatk_standard_hc(),
            &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
        )
        .expect("walk");
        let region = flatten_assembly_regions(&walk)
            .into_iter()
            .find(|r| {
                matches!(
                    call_disposition(r),
                    AssemblyRegionCallDisposition::ActiveFull
                ) && r.start.get() <= pos
                    && r.end.get() >= pos
            })
            .unwrap_or_else(|| panic!("{label} region"));
        let mut owned = region.clone();
        let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
        let assembled =
            assemble_reads_with_finalized(&mut owned, &dict, &mut ref_cache, &args.assemble)
                .expect("assemble");
        let padded_ref = assembly_reference_read(&dict, &mut ref_cache, &region).expect("pad");
        let graph_ref = create_graph_reference_read(&padded_ref, &region, &dict);
        let graph_reads = records_to_assembly_reads(&assembled.finalized_reads);
        let assembler = args.assemble.assembler.clone();
        let cyclic10 = cyclic_before_dangling(&graph_ref, &graph_reads, &assembler, 10);
        let seq10 = seq_kbest_haps(&graph_ref, &graph_reads, &assembler, 10);
        let seq25 = seq_kbest_haps(&graph_ref, &graph_reads, &assembler, 25);
        let seq_union = {
            let mut s = hap_hashes(&seq10);
            s.extend(hap_hashes(&seq25));
            s
        };
        let rt10 = extract_rt_haplotypes_before_remove_paths(
            &graph_ref,
            &graph_reads,
            &assembler,
            10,
            false,
            false,
        )
        .unwrap_or_default();
        let diag10 =
            diagnostic_extract_honoring_create_graph(&graph_ref, &graph_reads, &assembler, 10);
        let prod = assemble_from_ref_and_reads(&graph_ref, &graph_reads, &assembler).expect("prod");
        let prod_h = hap_hashes(&prod.haplotypes);
        let rt10_only: BTreeSet<_> = hap_hashes(&rt10).difference(&seq_union).cloned().collect();
        let acyclic_rt_needed = cyclic10 == Some(false)
            && seq10.iter().any(|h| !h.is_reference)
            && !rt10_only.is_empty();
        kv(
            label,
            format!(
                "cyclic10={cyclic10:?}\tseq10={}\tseq25={}\trt10={}\tdiag10={}\tassemble={}\trt10_not_in_seq={}\tacyclic_rt_needed={acyclic_rt_needed}",
                seq10.len(),
                seq25.len(),
                rt10.len(),
                diag10.len(),
                prod.haplotypes.len(),
                rt10_only.len()
            ),
        );
        (
            cyclic10,
            seq25.len(),
            prod_h.len(),
            acyclic_rt_needed,
            seq10.len(),
        )
    };

    let (ttc_cyc, ttc_seq25, ttc_asm, _, _) = probe(TG_INTERVAL, P12_TTC, "p12_ttc");
    let _ = probe(TG_INTERVAL, P12_TG, "p12_tg");
    let (mid_cyc, _, _, mid_needed, mid_seq10) = probe(MID_B_INTERVAL, CLOSED_MID_B, "p12_mid_b");
    assert_eq!(ttc_cyc, Some(true));
    assert_eq!(ttc_seq25, 6);
    assert_eq!(
        ttc_asm, 6,
        "cycle abort must not drop Java-materialized P12 TTC haplotypes"
    );
    kv("p12_ttc_java_n", "6");

    let indel4_ref = root.join(INDEL4_REF_REL);
    let indel4_bam = root.join(INDEL4_BAM_REL);
    let mut found_acyclic_required = mid_needed;
    if indel4_ref.is_file() && indel4_bam.is_file() {
        let dict4 = SequenceDictionary::from_fasta_path(&indel4_ref).expect("indel4 dict");
        let specs = parse_intervals_cli_string(&dict4, "chrIndel4:1-68").expect("indel4");
        let walk = traverse_assembly_region_walker(
            &dict4,
            &specs,
            &indel4_ref,
            &indel4_bam,
            &ReadFilterParams::gatk_standard_hc(),
            &WalkerTraversalConfig::gatk_haplotype_caller_production(0),
        )
        .expect("indel4 walk");
        let region = flatten_assembly_regions(&walk)
            .into_iter()
            .find(|r| {
                matches!(
                    call_disposition(r),
                    AssemblyRegionCallDisposition::ActiveFull
                )
            })
            .expect("indel4 active");
        let mut owned = region.clone();
        let mut ref_cache = ReferenceWindowCache::new(indel4_ref.clone(), 4);
        let assembled =
            assemble_reads_with_finalized(&mut owned, &dict4, &mut ref_cache, &args.assemble)
                .expect("indel4 assemble");
        let padded_ref = assembly_reference_read(&dict4, &mut ref_cache, &region).expect("pad");
        let graph_ref = create_graph_reference_read(&padded_ref, &region, &dict4);
        let graph_reads = records_to_assembly_reads(&assembled.finalized_reads);
        let assembler = args.assemble.assembler.clone();
        let cyclic10 = cyclic_before_dangling(&graph_ref, &graph_reads, &assembler, 10);
        let seq10 = seq_kbest_haps(&graph_ref, &graph_reads, &assembler, 10);
        let seq25 = seq_kbest_haps(&graph_ref, &graph_reads, &assembler, 25);
        let mut seq_union = hap_hashes(&seq10);
        seq_union.extend(hap_hashes(&seq25));
        let rt10 = extract_rt_haplotypes_before_remove_paths(
            &graph_ref,
            &graph_reads,
            &assembler,
            10,
            false,
            false,
        )
        .unwrap_or_default();
        let prod = assemble_from_ref_and_reads(&graph_ref, &graph_reads, &assembler).expect("prod");
        let rt10_only = hap_hashes(&rt10).difference(&seq_union).count();
        let needed = cyclic10 == Some(false) && rt10_only > 0;
        kv(
            "l2_indel4",
            format!(
                "cyclic10={cyclic10:?}\tseq10={}\tseq25={}\trt10={}\tassemble={}\trt10_not_in_seq={rt10_only}\tacyclic_rt_needed={needed}",
                seq10.len(),
                seq25.len(),
                rt10.len(),
                prod.haplotypes.len()
            ),
        );
        found_acyclic_required = true;
        assert_eq!(cyclic10, Some(false), "indel4 k=10 must stay acyclic");
        assert!(
            prod.haplotypes.iter().any(|h| !h.is_reference),
            "indel4 alt must survive cycle-only abort"
        );
    } else {
        kv("l2_indel4", "skip missing BAM");
    }

    kv(
        "acyclic_rt_positive",
        format!(
            "found={found_acyclic_required}\tmid_b_k10_cyclic={mid_cyc:?}\tmid_b_seq10={mid_seq10}"
        ),
    );
    kv(
        "acyclic_rt_rule",
        "candidate is abort-cyclic, not disable-RT: p11/indel4 acyclic k=10 remains extractable; P12 zip-collapse fallback uses use_seq_graph=false and is out of createGraph scope",
    );
    assert!(
        found_acyclic_required,
        "need a live acyclic lower-k control (indel4 k=10)"
    );
}
