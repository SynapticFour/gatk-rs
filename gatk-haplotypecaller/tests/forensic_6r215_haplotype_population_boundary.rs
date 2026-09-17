//! 6R.215: first causal haplotype-population arrow at `20:29455015 G/T`.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! 6R.214 closed proof-only: PairHMM hap_n Java 30 vs Rust 84/89. This round
//! locates the assembly/trim boundary that first admits a Rust-only haplotype.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r215_haplotype_population_boundary -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::assembly_region_finalize::{
    assembly_reference_read, create_graph_reference_read, records_to_assembly_reads,
};
use gatk_haplotypecaller::read_threading_assembler::{
    assemble_from_ref_and_reads, build_threading_graph_for_seq_assembly,
    extract_haplotypes_from_seq_kbest_paths, extract_rt_haplotypes_after_remove_paths,
    extract_rt_haplotypes_before_remove_paths, DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH,
};
use gatk_haplotypecaller::seq_graph::SeqGraph;
use gatk_haplotypecaller::seq_kbest_haplotype::find_best_haplotypes_seq_graph;
use gatk_haplotypecaller::{
    assemble_reads_with_finalized, begin_hap_list_observe, begin_trim_io_observe, call_disposition,
    flatten_assembly_regions, probe_seq_graph_kmer_attempts, take_hap_list_snaps,
    take_hap_list_trim_span, take_trim_io, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs, Cigar,
    HapListColumn, Haplotype, HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
    DEFAULT_STAND_EMIT_CONFIDENCE,
};
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 29_455_015;
const SIB: u64 = 29_455_019;
const JAVA_TRIM_START: u64 = 29_454_995;
const JAVA_TRIM_END: u64 = 29_455_124;
const JAVA_ACTIVE: &str = "20:29455000-29455145";
const JAVA_PADDED: &str = "20:29454900-29455245";
const JAVA_UNTRIMMED_N: usize = 78;
const JAVA_POST_TRIM_N: usize = 30;
const JAVA_KBEST_K: usize = 128;
const JAVA_KBEST_RETURNED: usize = 78;

/// Java PairHMM / post-trim FNV hashes (`hap-trim-at-loc` stage=trimmed).
const JAVA_POST_TRIM_HASHES: &[&str] = &[
    "a9a6a6a3c7e334fe",
    "17b253b3caefbc24",
    "8cdf90fd00055da9",
    "94152d17f55d2c0f",
    "504af452f2b4a369",
    "421f004a656474e4",
    "f61271afcbc2fb1d",
    "dd22d427216a7c13",
    "a075ca46f85f60ba",
    "66751242bb88dd8c",
    "f5d1e8d60cfb9666",
    "1ec62235b5a6f477",
    "ba566a88f4a37b49",
    "f9fb0ec02bc3a8ef",
    "e24a3d9b7270ab1e",
    "a288c9e831abc8c3",
    "65773746fc8fc8cd",
    "80b1e7bc60742470",
    "8a328239b600c7bd",
    "cf7d9939fecc0fb3",
    "fad3d8d2e57bba5a",
    "a1c7814baa41b868",
    "04cd4d695bcfcef2",
    "a911d3c74c6ac1ab",
    "7722ec426dd6775e",
    "8183d80853176284",
    "405bde68757e1909",
    "7c0149c997f1c2d1",
    "85fb4ea1d9d0a2b7",
    "fdfee8cd394ccda6",
];

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

const RUST_ONLY_WINDOW_HASHES: &[&str] = &[
    "0103d498d05ea41d",
    "0ff22ea407fcd612",
    "34844af71ff38140",
    "352750bcca900ff9",
    "377d32f330bb29c6",
    "4242c8e361353c25",
    "79c9e31c50bf2ad4",
    "8f89eb8e7e09f008",
    "92d5943204d3f9f6",
    "a2d31d7854bf7ae5",
    "c7f943e17b4f3367",
    "fca50f1c73d7f7be",
    "fe7b7c9bc838c98f",
];

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R215\t{key}\t{}", value.as_ref());
}

fn hex(h: u64) -> String {
    format!("{h:016x}")
}

fn fnv1a64_hex(bases: &[u8]) -> String {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in bases {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    format!("{h:016x}")
}

fn unique_hashes(cols: &[HapListColumn]) -> BTreeSet<String> {
    cols.iter().map(|c| hex(c.fnv1a)).collect()
}

fn lens(cols: &[HapListColumn]) -> String {
    let mut m: BTreeMap<usize, usize> = BTreeMap::new();
    for c in cols {
        *m.entry(c.len).or_insert(0) += 1;
    }
    m.iter()
        .map(|(k, v)| format!("{k}={v}"))
        .collect::<Vec<_>>()
        .join(",")
}

fn hap_hashes(haps: &[Haplotype]) -> BTreeSet<String> {
    haps.iter().map(|h| fnv1a64_hex(&h.bases)).collect()
}

fn window_bases(h: &Haplotype, win_start: u64, win_end: u64) -> Option<Vec<u8>> {
    let gl = h.genome_loc?;
    let gs = gl.start_1based();
    if win_start < gs || win_end < win_start {
        return None;
    }
    let off = (win_start - gs) as usize;
    let len = (win_end - win_start + 1) as usize;
    h.bases.get(off..off.checked_add(len)?).map(|s| s.to_vec())
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

#[test]
fn forensic_6r215_java_population_pin() {
    kv("java_pin", JAVA_PIN);
    kv("java_active", JAVA_ACTIVE);
    kv("java_padded", JAVA_PADDED);
    kv("java_untrimmed_n", JAVA_UNTRIMMED_N.to_string());
    kv("java_kbest_requested", JAVA_KBEST_K.to_string());
    kv("java_kbest_returned", JAVA_KBEST_RETURNED.to_string());
    kv("java_k10", "skip cycles_before_dangling");
    kv("java_k25_kbest_n", "78");
    kv("java_post_trim_n", JAVA_POST_TRIM_N.to_string());
    kv("java_trim_span", "20:29454995-29455124");
    kv("java_pairhmm_hap_n", "30");
    kv(
        "rust_only_window_n",
        RUST_ONLY_WINDOW_HASHES.len().to_string(),
    );
    assert_eq!(JAVA_UNTRIMMED_HASHES.len(), 78);
    assert_eq!(JAVA_POST_TRIM_HASHES.len(), 30);
    assert_eq!(RUST_ONLY_WINDOW_HASHES.len(), 13);
    assert_eq!(JAVA_KBEST_K, DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH);
    assert_eq!(JAVA_KBEST_RETURNED, JAVA_UNTRIMMED_N);
}

#[test]
fn forensic_6r215_live_haplotype_population_boundary() {
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

    begin_hap_list_observe();
    begin_trim_io_observe();
    let outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    if let Some(span) = take_hap_list_trim_span() {
        kv(
            "trim_span",
            format!(
                "active={}-{} trim={}-{}",
                span.active_start, span.active_end, span.trim_start, span.trim_end
            ),
        );
        kv(
            "java_trim_span",
            format!("{JAVA_TRIM_START}-{JAVA_TRIM_END}"),
        );
    }
    let snaps = take_hap_list_snaps();
    let io = take_trim_io();

    let mut by_stage: BTreeMap<&str, BTreeSet<String>> = BTreeMap::new();
    let mut n_by_stage: BTreeMap<&str, usize> = BTreeMap::new();
    for s in &snaps {
        let uniq = unique_hashes(&s.columns);
        let n_ref = s.columns.iter().filter(|c| c.is_reference).count();
        let mut counts: BTreeMap<String, usize> = BTreeMap::new();
        for c in &s.columns {
            *counts.entry(hex(c.fnv1a)).or_insert(0) += 1;
        }
        let dup = counts.values().filter(|n| **n > 1).count();
        kv(
            "snap",
            format!(
                "stage={}\tn_cols={}\tunique={}\tn_ref={}\tdup_seqs={}\tlens={}",
                s.stage,
                s.n,
                uniq.len(),
                n_ref,
                dup,
                lens(&s.columns)
            ),
        );
        n_by_stage.insert(s.stage, s.n);
        by_stage.insert(s.stage, uniq);
    }

    let mut n_dropped = 0usize;
    let mut n_kept = 0usize;
    let mut n_collapsed = 0usize;
    for row in &io {
        match row.outcome {
            "dropped" => n_dropped += 1,
            "kept" => n_kept += 1,
            "collapsed" => n_collapsed += 1,
            _ => {}
        }
    }
    kv(
        "trim_io_summary",
        format!("kept={n_kept} dropped={n_dropped} collapsed={n_collapsed}"),
    );

    let java_untrimmed: BTreeSet<&str> = JAVA_UNTRIMMED_HASHES.iter().copied().collect();
    let assemble_uniq = by_stage.get("after_assemble").cloned().unwrap_or_default();
    let assemble_refs: BTreeSet<&str> = assemble_uniq.iter().map(String::as_str).collect();
    let rust_only_untrimmed: Vec<_> = assemble_refs.difference(&java_untrimmed).copied().collect();
    let java_only_untrimmed: Vec<_> = java_untrimmed.difference(&assemble_refs).copied().collect();
    let common_untrimmed = java_untrimmed.intersection(&assemble_refs).count();
    kv(
        "untrimmed_set",
        format!(
            "common={common_untrimmed}\tjava_only={}\trust_only={}",
            java_only_untrimmed.len(),
            rust_only_untrimmed.len()
        ),
    );
    kv("java_only_untrimmed", java_only_untrimmed.join(","));
    kv(
        "rust_only_untrimmed_n",
        rust_only_untrimmed.len().to_string(),
    );
    for h in &rust_only_untrimmed {
        kv("rust_only_untrimmed", *h);
    }

    let assemble_n = *n_by_stage.get("after_assemble").unwrap_or(&0);
    let before_trim_n = *n_by_stage.get("before_trim").unwrap_or(&0);
    let after_trim_n = *n_by_stage.get("after_trim_to").unwrap_or(&0);
    let after_preserve_n = *n_by_stage.get("after_preserve_untrimmed").unwrap_or(&0);
    let after_prune_n = *n_by_stage.get("after_prune_spillover").unwrap_or(&0);
    let pairhmm_n = *n_by_stage.get("pairhmm_input").unwrap_or(&0);
    kv("assemble_n", assemble_n.to_string());
    kv("before_trim_n", before_trim_n.to_string());
    kv("after_trim_to_n", after_trim_n.to_string());
    kv("after_preserve_n", after_preserve_n.to_string());
    kv("after_prune_n", after_prune_n.to_string());
    kv("pairhmm_n", pairhmm_n.to_string());
    kv(
        "call_region_hap_n",
        outcome.assembly.haplotypes.len().to_string(),
    );
    kv("java_untrimmed_n", JAVA_UNTRIMMED_N.to_string());
    kv("java_pairhmm_n", JAVA_POST_TRIM_N.to_string());

    let extras: BTreeSet<&str> = RUST_ONLY_WINDOW_HASHES.iter().copied().collect();
    let mut window_present = 0usize;
    let mut extras_at_pairhmm = Vec::new();
    for h in &outcome.assembly.haplotypes {
        if let Some(b) = window_bases(h, JAVA_TRIM_START, JAVA_TRIM_END) {
            let wh = fnv1a64_hex(&b);
            if extras.contains(wh.as_str()) {
                extras_at_pairhmm.push(format!(
                    "{} full={} len={} kmer={} score={} cigar={}",
                    wh,
                    fnv1a64_hex(&h.bases),
                    h.bases.len(),
                    h.kmer_size,
                    h.score,
                    h.cigar
                        .as_ref()
                        .map(|c| c.to_gatk_string())
                        .unwrap_or_else(|| ".".into())
                ));
            }
            window_present += 1;
        }
    }
    kv(
        "java_window_coverable",
        format!("{window_present}/{}", outcome.assembly.haplotypes.len()),
    );
    let extra_window_unique: BTreeSet<String> = extras_at_pairhmm
        .iter()
        .filter_map(|row| row.split_whitespace().next().map(str::to_string))
        .collect();
    kv(
        "extras_at_pairhmm_window",
        format!(
            "rows={}/unique={}/pinned=13",
            extras_at_pairhmm.len(),
            extra_window_unique.len()
        ),
    );
    for row in &extras_at_pairhmm {
        kv("rust_only_window_provenance", row);
    }
    assert_eq!(extra_window_unique.len(), 0);

    assert_eq!(common_untrimmed, JAVA_UNTRIMMED_N);
    assert!(java_only_untrimmed.is_empty());
    assert_eq!(assemble_n, before_trim_n);
    assert_eq!(assemble_n, JAVA_UNTRIMMED_N);
    assert_eq!(after_trim_n, after_preserve_n);
    assert_eq!(after_preserve_n, pairhmm_n);
    assert_eq!(pairhmm_n, JAVA_POST_TRIM_N);

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
    for (i, h) in assembled.assembly.haplotypes.iter().enumerate() {
        let hash = fnv1a64_hex(&h.bases);
        if java_untrimmed.contains(hash.as_str()) {
            continue;
        }
        kv(
            "assemble_extra",
            format!(
                "idx={i}\thash={hash}\tisRef={}\tlen={}\tscore={}\tkmer={}\tcigar={}",
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

    let padded_ref = assembly_reference_read(&dict, &mut ref_cache, region).expect("pad ref");
    let graph_ref = create_graph_reference_read(&padded_ref, region, &dict);
    let graph_reads = records_to_assembly_reads(&assembled.finalized_reads);
    kv(
        "graph_ref",
        format!(
            "len={}\thash={}",
            graph_ref.bases.len(),
            fnv1a64_hex(&graph_ref.bases)
        ),
    );
    kv("assembly_reads", graph_reads.len().to_string());
    let assembler = args.assemble.assembler.clone();
    kv(
        "assembler",
        format!(
            "kmers={:?}\tuse_seq_graph={}\tnum_best={}\tdangling_java_exact={}",
            assembler.kmer_sizes,
            assembler.use_seq_graph,
            assembler.num_best_haplotypes_per_graph,
            assembler.dangling_java_exact
        ),
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

    let seq25 = seq_kbest_haps(&graph_ref, &graph_reads, &assembler, 25);
    let seq10 = seq_kbest_haps(&graph_ref, &graph_reads, &assembler, 10);
    let seq25_h = hap_hashes(&seq25);
    let seq10_h = hap_hashes(&seq10);
    let seq25_refs: BTreeSet<&str> = seq25_h.iter().map(String::as_str).collect();
    let seq10_refs: BTreeSet<&str> = seq10_h.iter().map(String::as_str).collect();
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
    kv(
        "seq10_vs_java",
        format!(
            "common={}\tjava_only={}\trust_only={}",
            java_untrimmed.intersection(&seq10_refs).count(),
            java_untrimmed.difference(&seq10_refs).count(),
            seq10_refs.difference(&java_untrimmed).count()
        ),
    );

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
    let rt10_after = extract_rt_haplotypes_after_remove_paths(
        &graph_ref,
        &graph_reads,
        &assembler,
        10,
        false,
        false,
    )
    .unwrap_or_default();
    let rt25_after = extract_rt_haplotypes_after_remove_paths(
        &graph_ref,
        &graph_reads,
        &assembler,
        25,
        false,
        false,
    )
    .unwrap_or_default();
    let rt10b_h = hap_hashes(&rt10_before);
    let rt25b_h = hap_hashes(&rt25_before);
    let rt10a_h = hap_hashes(&rt10_after);
    let rt25a_h = hap_hashes(&rt25_after);
    let rt10b_refs: BTreeSet<&str> = rt10b_h.iter().map(String::as_str).collect();
    let rt25b_refs: BTreeSet<&str> = rt25b_h.iter().map(String::as_str).collect();
    kv("rt_before_k10_n", rt10_before.len().to_string());
    kv("rt_before_k25_n", rt25_before.len().to_string());
    kv("rt_after_k10_n", rt10_after.len().to_string());
    kv("rt_after_k25_n", rt25_after.len().to_string());

    let extras_in_rt10: Vec<_> = rust_only_untrimmed
        .iter()
        .copied()
        .filter(|h| rt10b_refs.contains(h))
        .collect();
    let extras_in_rt25: Vec<_> = rust_only_untrimmed
        .iter()
        .copied()
        .filter(|h| rt25b_refs.contains(h))
        .collect();
    let extras_in_seq25: Vec<_> = rust_only_untrimmed
        .iter()
        .copied()
        .filter(|h| seq25_refs.contains(h))
        .collect();
    let extras_in_seq10: Vec<_> = rust_only_untrimmed
        .iter()
        .copied()
        .filter(|h| seq10_refs.contains(h))
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
        "extras_in_seq_kbest_k10",
        format!("{}/{}", extras_in_seq10.len(), rust_only_untrimmed.len()),
    );
    kv(
        "extras_in_seq_kbest_k25",
        format!("{}/{}", extras_in_seq25.len(), rust_only_untrimmed.len()),
    );

    let prod = assemble_from_ref_and_reads(&graph_ref, &graph_reads, &assembler).expect("prod");
    kv(
        "assemble_from_ref_and_reads_n",
        prod.haplotypes.len().to_string(),
    );

    let seq25_matches_java = seq25.len() == JAVA_UNTRIMMED_N
        && java_untrimmed.difference(&seq25_refs).count() == 0
        && seq25_refs.difference(&java_untrimmed).count() == 0;
    kv(
        "seq25_matches_java_untrimmed",
        seq25_matches_java.to_string(),
    );

    let first_boundary;
    let classification;
    let arrow;
    if seq25_matches_java
        && rust_only_untrimmed.is_empty()
        && assemble_n == JAVA_UNTRIMMED_N
        && after_trim_n == JAVA_POST_TRIM_N
    {
        first_boundary = "none at 20:29455015";
        classification = "CLOSED_AFTER_6R218";
        arrow = "6R.218 SeqGraph-path cycle abort: SeqGraph k=25 = 78, cyclic k=10 contributes 0, trim = Java 30";
    } else if !seq25_matches_java
        && extras_in_seq25.len() == rust_only_untrimmed.len()
        && !rust_only_untrimmed.is_empty()
    {
        first_boundary = "SeqGraph findBestPaths (k=25)";
        classification = "KBEST_HAPLOTYPE_POPULATION_DIVERGENCE";
        arrow = "Rust SeqGraph k=25 k-best already contains the Rust-only untrimmed haplotypes; Java k=25 returns 78";
    } else if seq25_matches_java
        && extras_in_rt10.len() == rust_only_untrimmed.len()
        && !rust_only_untrimmed.is_empty()
    {
        first_boundary = "merge_rt_kbest_pre_remove_paths after SeqGraph k=25 findBestPaths";
        classification = "HAPLOTYPE_MATERIALIZATION_DIVERGENCE";
        arrow = "Java SeqGraph k=25 findBestPaths (K=128) returns 78 haplotypes and skips k=10 (cycles_before_dangling). Rust SeqGraph k=25 matches those 78 exactly. merge_rt_kbest_pre_remove_paths then extracts RT-before-remove k=10 and inserts 44 unique kmer=10 haplotypes Java never materializes";
    } else if seq25_matches_java && extras_in_rt25.len() == rust_only_untrimmed.len() {
        first_boundary = "RT k-best before remove_paths (k=25) merged after SeqGraph";
        classification = "HAPLOTYPE_MATERIALIZATION_DIVERGENCE";
        arrow = "SeqGraph k=25 matches Java 78; RT-before-remove k=25 then admits the Rust-only haplotypes";
    } else if assemble_n != JAVA_UNTRIMMED_N {
        first_boundary = "after_assemble (before trim)";
        classification = "HAPLOTYPE_MATERIALIZATION_DIVERGENCE";
        arrow = "Rust after_assemble already differs from Java untrimmed 78; trim is not the first split";
    } else if after_trim_n != JAVA_POST_TRIM_N {
        first_boundary = "trim_to";
        classification = "HAPLOTYPE_TRIM_DIVERGENCE";
        arrow = "untrimmed populations match; trim changes unique sequences";
    } else {
        first_boundary = "post-trim";
        classification = "HAPLOTYPE_MATERIALIZATION_DIVERGENCE";
        arrow = "trim matched Java 30; a later supplement added haplotypes";
    }
    kv("first_boundary", first_boundary);
    kv("classification", classification);
    kv("first_divergent_arrow", arrow);
    kv(
        "production_change",
        "NONE in 6R.215; closed upstream by 6R.218 GRAPH_STATE_DIVERGENCE",
    );
    kv("kbest_resource_policy", "unchanged legacy_1024");
    let _ = (rt10a_h, rt25a_h);

    assert_eq!(assemble_n, JAVA_UNTRIMMED_N);
    assert_eq!(n_dropped, 0, "trim drops none; it collapses");
    assert_eq!(classification, "CLOSED_AFTER_6R218");

    let emitted =
        try_emit_call_region_variants(region, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    let sib = emitted
        .iter()
        .find(|r| r.position == SIB && r.reference == "G" && r.alternate.iter().any(|a| a == "A"))
        .expect("sib");
    assert_eq!(
        sib.samples[0].pl.as_deref(),
        Some([645, 0, 1147].as_slice())
    );
    let tgt = emitted
        .iter()
        .find(|r| r.position == TARGET && r.reference == "G")
        .expect("target");
    assert_eq!(tgt.samples[0].pl.as_deref(), Some([69, 0, 2140].as_slice()));
}
