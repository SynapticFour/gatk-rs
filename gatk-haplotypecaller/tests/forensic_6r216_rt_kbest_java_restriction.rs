//! 6R.216: Java-equivalent restriction on post-SeqGraph RT k=10 merge
//! at `20:29455015 G/T`.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! Default HC: `generateSeqGraph = !useLinkedDeBruijnGraph` (false → true).
//! Proof-only. PRODUCTION CHANGE: NONE.
//!
//! Java `ReadThreadingAssembler` (pinned SHA):
//! - `runLocalAssembly` → `assembleKmerGraphsAndHaplotypeCall` when
//!   `generateSeqGraph`
//! - `assemble()` tries **every** configured kmer via `createGraph` (does not
//!   stop after a successful higher k)
//! - `createGraph`: `if (generateSeqGraph && rtgraph.hasCycles()) return null`
//!   **before** `recoverDanglingTails`
//! - successful SeqGraphs with variation → `findBestPaths` only
//! - no post-SeqGraph ReadThreadingGraph k-best merge
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r216_rt_kbest_java_restriction -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::assembly_region_finalize::{
    assembly_reference_read, create_graph_reference_read, records_to_assembly_reads,
};
use gatk_haplotypecaller::read_threading_assembler::{
    assemble_from_ref_and_reads, build_threading_graph_for_seq_assembly,
    extract_haplotypes_from_seq_kbest_paths, extract_rt_haplotypes_before_remove_paths,
    merge_rt_kbest_pre_remove_paths, DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH,
};
use gatk_haplotypecaller::seq_graph::SeqGraph;
use gatk_haplotypecaller::seq_kbest_haplotype::find_best_haplotypes_seq_graph;
use gatk_haplotypecaller::{
    assemble_reads_with_finalized, call_disposition, flatten_assembly_regions,
    probe_seq_graph_kmer_attempts, traverse_assembly_region_walker, AssemblyRegionCallDisposition,
    CallRegionArgs, Cigar, Haplotype, ReadFilterParams, WalkerTraversalConfig,
};
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const P11_REF_REL: &str = "parity/fixtures/p5_live_reference.fa";
const P11_BAM_REL: &str = "parity/build/sam-indexed-bam/p11_java_positive.bam";
const P12_BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const TARGET: u64 = 29_455_015;
const JAVA_ACTIVE: &str = "20:29455000-29455145";
const JAVA_PADDED: &str = "20:29454900-29455245";
const JAVA_UNTRIMMED_N: usize = 78;
const JAVA_KBEST_K: usize = 128;
const P12_TTC: u64 = 92_307_324;
const TG_INTERVAL: &str = "2:92307200-92307550";

/// Java `assembleReads` / untrimmed FNV hashes (`hap-trim-at-loc` stage=untrimmed).
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

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R216\t{key}\t{}", value.as_ref());
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

fn seq_kbest_haps(
    graph_ref: &gatk_haplotypecaller::AssemblyRead,
    graph_reads: &[gatk_haplotypecaller::AssemblyRead],
    assembler: &gatk_haplotypecaller::ReadThreadingAssemblerArgs,
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

/// Java 4.4 `ReadThreadingAssembler` call-graph pins (source SHA in JAVA_PIN).
#[test]
fn forensic_6r216_java_call_graph_pin() {
    kv("java_pin", JAVA_PIN);
    kv(
        "java_runLocalAssembly",
        "if (generateSeqGraph) assembleKmerGraphsAndHaplotypeCall else assembleGraphsAndExpandKmersGivenHaplotypes",
    );
    kv(
        "java_generateSeqGraph",
        "!useLinkedDeBruijnGraph; default HC useLinkedDeBruijnGraph=false so generateSeqGraph=true",
    );
    kv(
        "java_assembleKmerGraphsAndHaplotypeCall",
        "assemble() then findBestPaths on SeqGraphs with ASSEMBLED_SOME_VARIATION only",
    );
    kv(
        "java_assemble_loop",
        "for (kmerSize : kmerSizes) createGraph(...); expand only if results.isEmpty() && !dontIncreaseKmerSizesForCycles",
    );
    kv(
        "java_createGraph_cycles",
        "if (generateSeqGraph && rtgraph.hasCycles()) return null; BEFORE recoverDanglingTails",
    );
    kv(
        "java_findBestPaths",
        "GraphBasedKBestHaplotypeFinder on SeqGraph only; no ReadThreadingGraph k-best after toSequenceGraph",
    );
    kv(
        "java_no_merge_rt",
        "no counterpart of merge_rt_kbest_pre_remove_paths on the default SeqGraph path",
    );
    kv(
        "too_broad_rule",
        "if higher_k_assembly_succeeded: do not merge lower_k_rt_haplotypes — REJECTED: Java still calls createGraph for every configured kmer",
    );
    assert_eq!(JAVA_UNTRIMMED_HASHES.len(), JAVA_UNTRIMMED_N);
    assert_eq!(JAVA_KBEST_K, DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH);
    assert!(
        gatk_haplotypecaller::ReadThreadingAssemblerArgs::default().use_seq_graph,
        "Rust default must follow Java generateSeqGraph=true"
    );
    assert!(
        gatk_haplotypecaller::ReadThreadingAssemblerArgs::default().abort_seq_graph_on_cycles,
        "Rust SeqGraph createGraph aborts cycles"
    );
}

#[test]
fn forensic_6r216_live_rt_kbest_java_restriction() {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing chr20_tiny BAM/REF");
        return;
    }
    kv("java_pin", JAVA_PIN);
    kv("kbest_requested", JAVA_KBEST_K.to_string());
    kv("kbest_policy", "legacy_1024 (frontier, not K)");
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
    let covering: Vec<_> = regions
        .iter()
        .filter(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= TARGET
                && r.end.get() >= TARGET
        })
        .collect();
    assert_eq!(covering.len(), 1);
    let region = covering[0];
    let active = format!(
        "{}:{}-{}",
        region.contig,
        region.start.get(),
        region.end.get()
    );
    let padded = format!(
        "{}:{}-{}",
        region.contig,
        region.extended_start.get(),
        region.extended_end.get()
    );
    kv("active", &active);
    kv("padded", &padded);
    assert_eq!(active, JAVA_ACTIVE);
    assert_eq!(padded, JAVA_PADDED);

    let mut owned = region.clone();
    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
    let args = CallRegionArgs::strict_java();
    let assembled =
        assemble_reads_with_finalized(&mut owned, &dict, &mut ref_cache, &args.assemble)
            .expect("assemble");
    kv(
        "assemble_reads_n",
        assembled.assembly.haplotypes.len().to_string(),
    );

    let padded_ref = assembly_reference_read(&dict, &mut ref_cache, region).expect("pad ref");
    let graph_ref = create_graph_reference_read(&padded_ref, region, &dict);
    let graph_reads = records_to_assembly_reads(&assembled.finalized_reads);
    let assembler = args.assemble.assembler.clone();
    kv(
        "assembler",
        format!(
            "kmers={:?}\tuse_seq_graph={}\tabort_seq_graph_on_cycles={}\tnum_best={}\tdangling_java_exact={}",
            assembler.kmer_sizes,
            assembler.use_seq_graph,
            assembler.abort_seq_graph_on_cycles,
            assembler.num_best_haplotypes_per_graph,
            assembler.dangling_java_exact
        ),
    );
    assert!(assembler.use_seq_graph);
    assert!(assembler.abort_seq_graph_on_cycles);
    assert!(assembler.dangling_java_exact);
    assert_eq!(assembler.kmer_sizes, vec![10, 25]);
    assert_eq!(
        assembler.num_best_haplotypes_per_graph,
        DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH
    );

    let probes =
        probe_seq_graph_kmer_attempts(&graph_ref, &graph_reads, &assembler).expect("probe");
    for p in &probes {
        if p.kmer_size > 35 && p.phase == "expanded" {
            continue;
        }
        kv(
            "seq_probe",
            format!(
                "phase={}\tkmer={}\toutcome={}\trt_nodes={}\trt_edges={}\tcleanup={}\tkbest_n={}\textracted={}\tnon_ref={}",
                p.phase,
                p.kmer_size,
                p.outcome,
                p.thread_nodes,
                p.thread_edges,
                p.cleanup_status,
                p.kbest_paths,
                p.extracted_haps,
                p.non_ref_haps
            ),
        );
    }
    kv(
        "java_k10_createGraph",
        "return null (cycles_before_dangling)",
    );
    kv(
        "probe_note",
        "probe_seq_graph_kmer_attempts uses haplotype-dump builder (abort_cyclic=false); createGraph abort is seq_assembly graph",
    );

    let seq_graph_k10 = build_threading_graph_for_seq_assembly(
        &graph_ref,
        &graph_reads,
        10,
        &assembler,
        false,
        false,
    )
    .expect("seq k10 build");
    kv(
        "seq_assembly_graph_k10",
        if seq_graph_k10.is_some() {
            "Some (would violate Java createGraph)"
        } else {
            "None (Java createGraph null)"
        },
    );
    assert!(
        seq_graph_k10.is_none(),
        "build_threading_graph_for_seq_assembly(k=10) must match Java createGraph null"
    );

    let java_untrimmed: BTreeSet<&str> = JAVA_UNTRIMMED_HASHES.iter().copied().collect();
    let seq25 = seq_kbest_haps(&graph_ref, &graph_reads, &assembler, 25);
    let seq10 = seq_kbest_haps(&graph_ref, &graph_reads, &assembler, 10);
    let seq25_h = hap_hashes(&seq25);
    let seq25_refs: BTreeSet<&str> = seq25_h.iter().map(String::as_str).collect();
    kv("seq_kbest_k10_n", seq10.len().to_string());
    kv("seq_kbest_k25_n", seq25.len().to_string());
    kv(
        "seq25_vs_java",
        format!(
            "common={}\tjava_only={}\trust_only={}",
            java_untrimmed.intersection(&seq25_refs).count(),
            java_untrimmed.difference(&seq25_refs).count(),
            seq25_refs.difference(&java_untrimmed).count()
        ),
    );
    assert_eq!(seq10.len(), 0, "Java never SeqGraph-k-bests cyclic k=10");
    assert_eq!(seq25.len(), JAVA_UNTRIMMED_N);
    assert_eq!(java_untrimmed.difference(&seq25_refs).count(), 0);
    assert_eq!(seq25_refs.difference(&java_untrimmed).count(), 0);

    let rt10_before = extract_rt_haplotypes_before_remove_paths(
        &graph_ref,
        &graph_reads,
        &assembler,
        10,
        false,
        false,
    )
    .unwrap_or_default();
    let rt25_before = extract_rt_haplotypes_before_remove_paths(
        &graph_ref,
        &graph_reads,
        &assembler,
        25,
        false,
        false,
    )
    .unwrap_or_default();
    let rt10_h = hap_hashes(&rt10_before);
    let rt25_h = hap_hashes(&rt25_before);
    let rt10_refs: BTreeSet<&str> = rt10_h.iter().map(String::as_str).collect();
    let rt25_refs: BTreeSet<&str> = rt25_h.iter().map(String::as_str).collect();
    kv("rt_before_k10_n", rt10_before.len().to_string());
    kv("rt_before_k25_n", rt25_before.len().to_string());
    kv(
        "rt_k10_reason",
        "extract_rt_haplotypes_before_remove_paths now honors use_seq_graph && abort_seq_graph_on_cycles (6R.218)",
    );

    let rust_only: Vec<String> = seq25_h
        .symmetric_difference(&hap_hashes(&assembled.assembly.haplotypes))
        .cloned()
        .collect();
    let assemble_h = hap_hashes(&assembled.assembly.haplotypes);
    let assemble_refs: BTreeSet<&str> = assemble_h.iter().map(String::as_str).collect();
    let rust_only_untrimmed: Vec<&str> =
        assemble_refs.difference(&java_untrimmed).copied().collect();
    let java_only: Vec<&str> = java_untrimmed.difference(&assemble_refs).copied().collect();
    kv(
        "untrimmed_set",
        format!(
            "common={}\tjava_only={}\trust_only={}",
            java_untrimmed.intersection(&assemble_refs).count(),
            java_only.len(),
            rust_only_untrimmed.len()
        ),
    );
    assert!(java_only.is_empty());
    assert_eq!(
        rust_only_untrimmed.len(),
        0,
        "6R.218 SeqGraph-path cycle abort removes the 44 extras upstream"
    );
    for h in &rust_only_untrimmed {
        kv("rust_only_untrimmed", *h);
    }
    let _ = rust_only;

    let extras_in_rt10: Vec<_> = rust_only_untrimmed
        .iter()
        .copied()
        .filter(|h| rt10_refs.contains(h))
        .collect();
    let extras_in_rt25: Vec<_> = rust_only_untrimmed
        .iter()
        .copied()
        .filter(|h| rt25_refs.contains(h))
        .collect();
    let extras_in_seq25: Vec<_> = rust_only_untrimmed
        .iter()
        .copied()
        .filter(|h| seq25_refs.contains(h))
        .collect();
    kv(
        "extras_in_rt_before_k10",
        format!("{}/{}", extras_in_rt10.len(), rust_only_untrimmed.len()),
    );
    kv(
        "extras_in_rt_before_k25",
        format!("{}/{}", extras_in_rt25.len(), rust_only_untrimmed.len()),
    );
    kv(
        "extras_in_seq_kbest_k25",
        format!("{}/{}", extras_in_seq25.len(), rust_only_untrimmed.len()),
    );
    assert_eq!(extras_in_rt10.len(), 0);
    assert_eq!(extras_in_rt25.len(), 0);
    assert_eq!(extras_in_seq25.len(), 0);

    for (i, h) in rt10_before.iter().enumerate() {
        let hash = fnv1a64_hex(&h.bases);
        if !java_untrimmed.contains(hash.as_str()) {
            kv(
                "k10_candidate",
                format!(
                    "idx={i}\thash={hash}\tisRef={}\tlen={}\tscore={}\tkmer={}\tcigar={}\tsource=RT_before_remove_k10\taccepted=merge_rt unique insert",
                    h.is_reference,
                    h.bases.len(),
                    h.score,
                    h.kmer_size,
                    h.cigar
                        .as_ref()
                        .map(|c| c.to_gatk_string())
                        .unwrap_or_else(|| ".".into())
                ),
            );
        }
    }

    let mut merge_in = seq25.clone();
    let merge_in_n = merge_in.len();
    merge_rt_kbest_pre_remove_paths(
        &graph_ref,
        &graph_reads,
        &assembler,
        &[10, 25],
        &mut merge_in,
    )
    .expect("merge_rt");
    let merge_out_h = hap_hashes(&merge_in);
    let merge_out_refs: BTreeSet<&str> = merge_out_h.iter().map(String::as_str).collect();
    let merge_added: Vec<&str> = merge_out_refs.difference(&seq25_refs).copied().collect();
    kv("merge_input_n", merge_in_n.to_string());
    kv("merge_output_n", merge_in.len().to_string());
    kv("merge_added_n", merge_added.len().to_string());
    kv(
        "merge_source_k25",
        "SeqGraph findBestPaths k=25 = 78 Java hashes",
    );
    kv(
        "merge_source_k10",
        "RT before_remove k=10 (cyclic graph Java createGraph rejected)",
    );
    assert_eq!(merge_in_n, JAVA_UNTRIMMED_N);
    assert_eq!(
        merge_added.len(),
        0,
        "6R.218: merge_rt on cyclic k=10 now extracts nothing"
    );

    let prod = assemble_from_ref_and_reads(&graph_ref, &graph_reads, &assembler).expect("prod");
    kv(
        "assemble_from_ref_and_reads_n",
        prod.haplotypes.len().to_string(),
    );

    kv("q1_java_materializes_k10_rt", "NO");
    kv("q2_java_runs_lower_k_rt_after_k25", "createGraph(k=10) is still invoked independently; it returns null on cycles. Java does not k-best that RT graph.");
    kv(
        "q3_java_merges_rt_paths",
        "NO — findBestPaths is SeqGraph-only",
    );
    kv(
        "q4_java_condition",
        "assembleKmerGraphsAndHaplotypeCall materializes haplotypes only from SeqGraph findBestPaths of createGraph survivors; createGraph aborts generateSeqGraph&&hasCycles before dangling",
    );
    kv(
        "q5_condition_class",
        "lifecycle (no post-SeqGraph RT k-best) AND graph-state (k=10 cyclic abort). Not kmer-fallback: assemble() still tries every configured k.",
    );
    kv("classification", "GRAPH_STATE_DIVERGENCE");
    kv(
        "first_divergent_arrow",
        "closed by 6R.218: SeqGraph-path extract now aborts cyclic k=10 before dangling/k-best; extras 0; assemble 78",
    );
    kv(
        "proposed_production",
        "APPLIED in 6R.218: abort_cyclic_before_dangling = use_seq_graph && abort_seq_graph_on_cycles at extract_rt",
    );
    kv("global_rt_disable", "INCORRECT — L2 g2-subset-live / P12 still need merge_rt on acyclic graphs where SeqGraph zip drops alts");
}

#[test]
fn forensic_6r216_l2_p11_merge_rt_still_required() {
    let root = repo_root();
    let ref_fasta = root.join(P11_REF_REL);
    let bam = root.join(P11_BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing p11 L2 fixture");
        return;
    }
    kv(
        "neg_control",
        "L2 g2-subset-live p11_java_positive chrLive:1-63",
    );
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
    let regions = flatten_assembly_regions(&walk);
    let region = regions
        .iter()
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
    let padded_ref = assembly_reference_read(&dict, &mut ref_cache, region).expect("pad ref");
    let graph_ref = create_graph_reference_read(&padded_ref, region, &dict);
    let graph_reads = records_to_assembly_reads(&assembled.finalized_reads);
    let assembler = args.assemble.assembler.clone();

    let seq10 = seq_kbest_haps(&graph_ref, &graph_reads, &assembler, 10);
    let seq25 = seq_kbest_haps(&graph_ref, &graph_reads, &assembler, 25);
    let seq_graph_k10 = build_threading_graph_for_seq_assembly(
        &graph_ref,
        &graph_reads,
        10,
        &assembler,
        false,
        false,
    )
    .expect("p11 seq k10");
    let seq_graph_k25 = build_threading_graph_for_seq_assembly(
        &graph_ref,
        &graph_reads,
        25,
        &assembler,
        false,
        false,
    )
    .expect("p11 seq k25");
    let prod = assemble_from_ref_and_reads(&graph_ref, &graph_reads, &assembler).expect("prod");
    kv("p11_seq_k10_n", seq10.len().to_string());
    kv("p11_seq_k25_n", seq25.len().to_string());
    kv(
        "p11_seq_graph_k10",
        if seq_graph_k10.is_some() {
            "Some (acyclic — createGraph would proceed)"
        } else {
            "None (cyclic abort)"
        },
    );
    kv(
        "p11_seq_graph_k25",
        if seq_graph_k25.is_some() {
            "Some"
        } else {
            "None"
        },
    );
    kv("p11_assemble_n", prod.haplotypes.len().to_string());
    kv(
        "p11_alt_n",
        prod.haplotypes
            .iter()
            .filter(|h| !h.is_reference)
            .count()
            .to_string(),
    );
    kv(
        "p11_java_haplotype_count",
        "2 (g2-subset-live p11_java_positive_chrlive.tsv)",
    );
    assert!(
        prod.haplotypes.iter().any(|h| !h.is_reference),
        "L2 p11 must keep an alt; global RT-merge disable is not Java-equivalent"
    );
    kv(
        "p11_rule",
        "k=10 is acyclic so Java createGraph proceeds; SeqGraph k=10 already n=2 vs k=25 n=1. Skipping lower-k after higher-k success is not Java and would drop the alt. Cyclic-k=10 abort in merge_rt would not hit this fixture.",
    );
}

#[test]
fn forensic_6r216_p12_ttc_merge_rt_still_required() {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(P12_BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing P12 BAM/REF");
        return;
    }
    kv("neg_control", "P12 TTC cluster 2:92307324");
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, TG_INTERVAL).expect("interval");
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
            ) && r.start.get() <= P12_TTC
                && r.end.get() >= P12_TTC
        })
        .expect("p12 ttc");
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
    let seq_graph_k10 = build_threading_graph_for_seq_assembly(
        &graph_ref,
        &graph_reads,
        10,
        &assembler,
        false,
        false,
    )
    .expect("p12 seq k10");
    let seq25 = seq_kbest_haps(&graph_ref, &graph_reads, &assembler, 25);
    let prod = assemble_from_ref_and_reads(&graph_ref, &graph_reads, &assembler).expect("prod");
    kv(
        "p12_seq_graph_k10",
        if seq_graph_k10.is_some() {
            "Some (acyclic)"
        } else {
            "None (cyclic abort)"
        },
    );
    kv("p12_seq_k25_n", seq25.len().to_string());
    kv("p12_assemble_n", prod.haplotypes.len().to_string());
    kv(
        "p12_rule",
        "this TTC snapshot: k=10 cyclic (createGraph null), SeqGraph k=25 n=6 = assemble n=6. Coupled-event recovery on other P12 windows still uses merge_rt/supplement on acyclic kmers; do not disable RT after any SeqGraph success.",
    );
    assert!(
        !prod.haplotypes.is_empty(),
        "P12 TTC assembly must remain populated"
    );
}
