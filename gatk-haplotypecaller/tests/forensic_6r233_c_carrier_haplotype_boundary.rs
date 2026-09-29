//! 6R.233: first assembly operation that materializes Java C-bearing
//! haplotype `c7acc50dfb9f9ecc` while Rust does not. Target site remains
//! `20:29455379 G/A`. Missing event: `20:29455314 G>C`.
//!
//! Proof-only. Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! PRODUCTION CHANGE: NONE.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r233_c_carrier_haplotype_boundary -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::assembly::{
    AssemblyGraph, AssemblyGraphParams, AssemblyGraphPruningParams,
};
use gatk_haplotypecaller::assembly_based_caller::assemble_reads_with_finalized;
use gatk_haplotypecaller::assembly_region_finalize::{
    assembly_reference_read, create_graph_reference_read, records_to_assembly_reads,
};
use gatk_haplotypecaller::cigar::CigarOperator;
use gatk_haplotypecaller::event_map::build_per_haplotype_variation_events;
use gatk_haplotypecaller::haplotype::Haplotype;
use gatk_haplotypecaller::read_threading_assembler::{
    build_threading_graph_for_seq_assembly, extract_haplotypes_from_seq_kbest_paths,
    probe_seq_graph_kmer_attempts, DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH,
};
use gatk_haplotypecaller::read_unclip::alignment_end_1based;
use gatk_haplotypecaller::seq_graph::{SeqGraph, SeqGraphCleanupStatus};
use gatk_haplotypecaller::seq_kbest_haplotype::find_best_haplotypes_seq_graph;
use gatk_haplotypecaller::{
    assembly_graph_from_ref_and_reads_threading, call_disposition, flatten_assembly_regions,
    traverse_assembly_region_walker, AssemblyRegionCallDisposition, CallRegionArgs, Cigar,
    KmerSize, ReadFilterParams, WalkerTraversalConfig,
};
use rust_htslib::bam::record::Cigar as BamCigar;
use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::path::{Path, PathBuf};

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 29_455_379;
const EVENT_POS: u64 = 29_455_314;
const JAVA_CARRIER: &str = "c7acc50dfb9f9ecc";
const JAVA_CARRIER_SCORE: &str = "-2.64786702";
const JAVA_CARRIER_RANK: usize = 116;
const JAVA_KBEST_K: usize = 128;
const JAVA_C_VERTEX: &str = "ACCTCTGTTTTGTCATCTCTAATAAATACAGACACAATATTC";
const JAVA_G_VERTEX: &str = "TCCTCTGTTTTGTCATCTGTAATAAATACAGATACAATATTT";
const JAVA_C_EDGE_MULT: u32 = 3;
const JAVA_TSV: &str = include_str!("forensic_6r233_java_seqgraph.tsv");
const C_READS: &[(&str, u16)] = &[
    ("HWI-D00360:8:H88U0ADXX:2:1203:13289:20746", 83),
    ("HISEQ1:9:H8962ADXX:2:2114:7354:86795", 147),
    ("HWI-D00360:6:H81VLADXX:2:1106:13656:97971", 147),
    ("HWI-D00360:6:H81VLADXX:2:1211:6676:53019", 83),
];

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R233\t{key}\t{}", value.as_ref());
}

fn fnv1a64_hex(data: &[u8]) -> String {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in data {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    format!("{h:016x}")
}

fn parse_kv_fields(line: &str) -> BTreeMap<String, String> {
    let mut out = BTreeMap::new();
    for p in line.split('\t').skip(2) {
        if let Some((k, v)) = p.split_once('=') {
            out.insert(k.to_string(), v.to_string());
        }
    }
    out
}

fn hap_base_at(h: &Haplotype, pos: u64) -> Option<u8> {
    let gl = h.genome_loc?;
    if pos < gl.start_1based() || pos > gl.end_1based() {
        return None;
    }
    let cigar = h.cigar.as_ref()?;
    let mut ref_pos = gl.start_1based();
    let mut hap_pos = 0usize;
    for el in &cigar.elements {
        match el.operator {
            CigarOperator::Match => {
                for _ in 0..el.length {
                    if ref_pos == pos {
                        return h.bases.get(hap_pos).copied();
                    }
                    ref_pos += 1;
                    hap_pos += 1;
                }
            }
            CigarOperator::Insertion | CigarOperator::SoftClip => {
                hap_pos += el.length;
            }
            CigarOperator::Deletion => {
                for _ in 0..el.length {
                    if ref_pos == pos {
                        return None;
                    }
                    ref_pos += 1;
                }
            }
            CigarOperator::HardClip => {}
        }
    }
    None
}

fn snp_kmers(ref_bases: &[u8], k: usize, off: usize, alt: u8) -> HashSet<Vec<u8>> {
    let mut out = HashSet::new();
    if off >= ref_bases.len() || k == 0 {
        return out;
    }
    let start0 = off.saturating_sub(k - 1);
    for s in start0..=off {
        if s + k > ref_bases.len() {
            break;
        }
        let mut km = ref_bases[s..s + k].to_vec();
        km[off - s] = alt;
        out.insert(km);
    }
    out
}

fn graph_has_kmers(graph: &AssemblyGraph, kmers: &HashSet<Vec<u8>>) -> (usize, u32) {
    let mut n = 0usize;
    let mut support = 0u32;
    for node in graph.nodes() {
        if kmers.contains(node.kmer.as_ref()) {
            n += 1;
            support = support.saturating_add(node.support);
        }
    }
    (n, support)
}

fn seq_has_c_vertex(seq: &SeqGraph) -> bool {
    seq.vertices().iter().any(|v| {
        v.sequence
            .windows(JAVA_C_VERTEX.len())
            .any(|w| w == JAVA_C_VERTEX.as_bytes())
    })
}

fn seq_has_g_vertex(seq: &SeqGraph) -> bool {
    seq.vertices().iter().any(|v| {
        v.sequence
            .windows(JAVA_G_VERTEX.len())
            .any(|w| w == JAVA_G_VERTEX.as_bytes())
    })
}

fn seq_c_edge_support(seq: &SeqGraph) -> Option<(usize, usize, u32, bool)> {
    let c_ids: Vec<usize> = seq
        .vertices()
        .iter()
        .filter(|v| {
            v.sequence
                .windows(JAVA_C_VERTEX.len())
                .any(|w| w == JAVA_C_VERTEX.as_bytes())
        })
        .map(|v| v.id)
        .collect();
    if c_ids.is_empty() {
        return None;
    }
    seq.edges()
        .iter()
        .find(|e| c_ids.contains(&e.to) || c_ids.contains(&e.from))
        .map(|e| (e.from, e.to, e.support, e.is_ref))
}

fn read_base_at(rec: &rust_htslib::bam::Record, pos: u64) -> Option<(u8, usize)> {
    if rec.tid() < 0 || rec.is_unmapped() {
        return None;
    }
    let mut ref_pos = (rec.pos() + 1) as u64;
    let seq = rec.seq();
    let mut read_pos = 0usize;
    for c in rec.cigar().iter() {
        match c {
            BamCigar::Match(n) | BamCigar::Equal(n) | BamCigar::Diff(n) => {
                for _ in 0..*n {
                    if ref_pos == pos {
                        return Some((seq.as_bytes()[read_pos], read_pos));
                    }
                    ref_pos += 1;
                    read_pos += 1;
                }
            }
            BamCigar::Ins(n) | BamCigar::SoftClip(n) => read_pos += *n as usize,
            BamCigar::Del(n) | BamCigar::RefSkip(n) => ref_pos += u64::from(*n),
            BamCigar::HardClip(_) | BamCigar::Pad(_) => {}
        }
    }
    None
}

fn read_c_kmers(rec: &rust_htslib::bam::Record, k: usize, pos: u64) -> HashSet<Vec<u8>> {
    let Some((b, qidx)) = read_base_at(rec, pos) else {
        return HashSet::new();
    };
    if b != b'C' {
        return HashSet::new();
    }
    let bases = rec.seq().as_bytes();
    let mut out = HashSet::new();
    if qidx >= bases.len() || k == 0 {
        return out;
    }
    let start0 = qidx.saturating_sub(k - 1);
    for s in start0..=qidx {
        if s + k > bases.len() {
            break;
        }
        out.insert(bases[s..s + k].to_vec());
    }
    out
}

fn threading_params(
    k: usize,
    assembler: &gatk_haplotypecaller::ReadThreadingAssemblerArgs,
) -> AssemblyGraphParams {
    AssemblyGraphParams {
        kmer_size: KmerSize::try_from_usize(k).expect("k"),
        min_base_quality: assembler.min_base_quality,
        min_edge_weight: 1,
        dangling_path_max_nodes: 0,
        max_haplotypes: assembler.num_best_haplotypes_per_graph,
        max_haplotype_bases: 4096,
        start_threading_only_at_existing_vertex: !assembler.recover_dangling_branches,
    }
}

#[test]
fn forensic_6r233_c_carrier_haplotype_boundary() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    kv(
        "kbest_cap",
        DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH.to_string(),
    );
    assert_eq!(DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, JAVA_KBEST_K);

    let mut saw_carrier = false;
    let mut saw_c_vertex = false;
    let mut saw_c_edge = false;
    let mut java_k10_skip = false;
    let mut java_k25_rt_nodes = 0usize;
    for line in JAVA_TSV.lines() {
        if line.contains("\tassemble\t") && line.contains(JAVA_CARRIER) {
            let f = parse_kv_fields(line);
            assert_eq!(f.get("idx").map(String::as_str), Some("115"));
            assert_eq!(f.get("len").map(String::as_str), Some("460"));
            assert_eq!(f.get("isRef").map(String::as_str), Some("false"));
            assert_eq!(f.get("score").map(String::as_str), Some(JAVA_CARRIER_SCORE));
            saw_carrier = true;
            kv(
                "java_carrier",
                format!(
                    "hash={JAVA_CARRIER} idx=115 len=460 score={JAVA_CARRIER_SCORE} isRef=false cigar=460M"
                ),
            );
        }
        if line.contains("\tvtx\tkmer=25\tid=12\t") {
            let f = parse_kv_fields(line);
            assert_eq!(f.get("seq").map(String::as_str), Some(JAVA_C_VERTEX));
            saw_c_vertex = true;
            kv("java_c_vertex", format!("id=12 seq={JAVA_C_VERTEX}"));
        }
        if line.contains("\tedge\tkmer=25\tfrom=0\tto=12\t") {
            let f = parse_kv_fields(line);
            assert_eq!(f.get("mult").map(String::as_str), Some("3"));
            assert_eq!(f.get("is_ref").map(String::as_str), Some("false"));
            saw_c_edge = true;
            kv(
                "java_c_edge",
                format!("0->12 mult={JAVA_C_EDGE_MULT} is_ref=false"),
            );
        }
        if line.contains("\tkmer_skip\tkmer=10\treason=cycles_before_dangling") {
            java_k10_skip = true;
        }
        if line.contains("\trt\tkmer=25\t") {
            let f = parse_kv_fields(line);
            java_k25_rt_nodes = f.get("nodes").and_then(|s| s.parse().ok()).unwrap_or(0);
            kv(
                "java_rt_k25",
                line.split('\t').skip(2).collect::<Vec<_>>().join("\t"),
            );
        }
        if line.contains("\tkbest\tkmer=25\tK=128\t") && line.contains(JAVA_CARRIER) {
            let f = parse_kv_fields(line);
            assert_eq!(
                f.get("rank").and_then(|s| s.parse::<usize>().ok()),
                Some(JAVA_CARRIER_RANK)
            );
            kv(
                "java_kbest_carrier",
                format!("K=128 rank={JAVA_CARRIER_RANK} score={JAVA_CARRIER_SCORE}"),
            );
        }
        if line.contains("\tseq\tkmer=25\t") {
            kv(
                "java_seq_k25",
                line.split('\t').skip(2).collect::<Vec<_>>().join("\t"),
            );
        }
    }
    assert!(saw_carrier, "Java fixture must pin carrier haplotype");
    assert!(
        saw_c_vertex,
        "Java fixture must pin C-bearing SeqGraph vertex 12"
    );
    assert!(
        saw_c_edge,
        "Java fixture must pin C-bearing SeqGraph edge 0->12"
    );
    assert!(
        java_k10_skip,
        "Java k=10 must abort on cycles before dangling"
    );
    assert_eq!(java_k25_rt_nodes, 1064);
    kv(
        "java_carrier_events",
        "29455296 T>A, 29455314 G>C, 29455328 T>C, 29455337 T>C on the same 42bp SeqGraph vertex 12",
    );
    kv(
        "java_graph_before_kbest",
        "C-bearing SeqGraph vertex+edge exist after cleanup, before k-best; k-best rank 116/128 selects the path",
    );

    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing chr20_tiny BAM/REF");
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
    let region = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= TARGET
                && r.end.get() >= TARGET
        })
        .expect("ActiveFull");
    kv(
        "active_region",
        format!("20:{}-{}", region.start.get(), region.end.get()),
    );
    assert_eq!(region.start.get(), 29_455_300);
    assert_eq!(region.end.get(), 29_455_559);

    let args = CallRegionArgs::strict_java();
    let mut region_for_assemble = region.clone();
    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
    let assembled = assemble_reads_with_finalized(
        &mut region_for_assemble,
        &dict,
        &mut ref_cache,
        &args.assemble,
    )
    .expect("assemble");
    let haps = &assembled.assembly.haplotypes;
    kv("rust_untrimmed_hap_n", haps.len().to_string());
    let rust_hashes: BTreeSet<String> = haps.iter().map(|h| fnv1a64_hex(&h.bases)).collect();
    kv(
        "rust_has_java_carrier_hash",
        rust_hashes.contains(JAVA_CARRIER).to_string(),
    );
    assert!(
        rust_hashes.contains(JAVA_CARRIER),
        "6R.239: Rust K=128 includes Java carrier {JAVA_CARRIER}"
    );

    let mut hap_base_counts: BTreeMap<char, usize> = BTreeMap::new();
    let mut n_cover = 0usize;
    for h in haps {
        if let Some(b) = hap_base_at(h, EVENT_POS) {
            n_cover += 1;
            *hap_base_counts.entry(b as char).or_insert(0) += 1;
        }
    }
    kv("rust_haps_covering_29455314", n_cover.to_string());
    kv(
        "rust_hap_base_29455314",
        hap_base_counts
            .iter()
            .map(|(b, n)| format!("{b}:{n}"))
            .collect::<Vec<_>>()
            .join(","),
    );
    assert!(hap_base_counts.get(&'C').copied().unwrap_or(0) >= 1);
    assert_eq!(n_cover, 128);

    let (full_ref, full_pad) = assembled.assembly.event_map_reference();
    let per_hap = build_per_haplotype_variation_events(
        haps,
        full_ref,
        full_pad,
        assembled.assembly.max_mnp_distance(),
        region.contig.as_str(),
    );
    let mut rust_eventmap_c = false;
    for i in 0..per_hap.hap_count() {
        for e in per_hap.events_for(i) {
            if e.start_1based.get() == EVENT_POS && e.ref_allele == "G" && e.alt_allele == "C" {
                rust_eventmap_c = true;
            }
        }
    }
    kv("rust_eventmap_has_G_to_C", rust_eventmap_c.to_string());
    assert!(rust_eventmap_c, "6R.239: EventMap G>C from K=128 carrier");

    let padded_ref = assembly_reference_read(&dict, &mut ref_cache, region).expect("pad ref");
    let graph_ref = create_graph_reference_read(&padded_ref, region, &dict);
    let graph_reads = records_to_assembly_reads(&assembled.finalized_reads);
    let assembler = args.assemble.assembler.clone();
    let pad = region.extended_start.get();
    let off = (EVENT_POS - pad) as usize;
    kv("graph_ref_len", graph_ref.bases.len().to_string());
    kv("graph_ref_hash", fnv1a64_hex(&graph_ref.bases));
    kv("snp_offset", off.to_string());
    assert_eq!(graph_ref.bases.len(), 460);
    assert_eq!(graph_ref.bases[off], b'G');
    assert_eq!(fnv1a64_hex(&graph_ref.bases), "442d1a5b4b11b07e");

    kv(
        "finalized_reads_n",
        assembled.finalized_reads.len().to_string(),
    );

    let mut orig_ok = 0usize;
    let mut fin_ok = 0usize;
    for (qname, flag) in C_READS {
        let orig = region
            .reads
            .iter()
            .find(|r| {
                let rec = r.as_ref();
                rec.qname() == qname.as_bytes() && rec.flags() == *flag
            })
            .map(|r| r.as_ref());
        let Some(orig) = orig else {
            kv("c_read_orig_missing", format!("{qname} FLAG{flag}"));
            continue;
        };
        orig_ok += 1;
        let orig_base = read_base_at(orig, EVENT_POS);
        kv(
            &format!("c_read_orig_{}", qname.rsplit(':').nth(1).unwrap_or(qname)),
            format!(
                "qname={qname} FLAG={flag} start={} end={} cigar={} base={:?}",
                orig.pos() + 1,
                alignment_end_1based(orig),
                orig.cigar().to_string(),
                orig_base.map(|(b, _)| b as char)
            ),
        );
        let fin = assembled
            .finalized_reads
            .iter()
            .find(|r| r.qname() == qname.as_bytes() && r.flags() == *flag);
        match fin {
            Some(fin) => {
                fin_ok += 1;
                let fin_base = read_base_at(fin, EVENT_POS);
                let ck = read_c_kmers(fin, 25, EVENT_POS);
                kv(
                    &format!(
                        "c_read_finalized_{}",
                        qname.rsplit(':').nth(1).unwrap_or(qname)
                    ),
                    format!(
                        "start={} end={} cigar={} base={:?} c_kmers_25={}",
                        fin.pos() + 1,
                        alignment_end_1based(fin),
                        fin.cigar().to_string(),
                        fin_base.map(|(b, _)| b as char),
                        ck.len()
                    ),
                );
            }
            None => kv(
                &format!(
                    "c_read_finalized_missing_{}",
                    qname.rsplit(':').nth(1).unwrap_or(qname)
                ),
                format!("{qname} FLAG{flag}"),
            ),
        }
    }
    kv("c_reads_orig", format!("{orig_ok}/4"));
    kv("c_reads_finalized", format!("{fin_ok}/4"));
    assert_eq!(orig_ok, 4);
    assert_eq!(fin_ok, 4);

    let probes =
        probe_seq_graph_kmer_attempts(&graph_ref, &graph_reads, &assembler).expect("probe");
    for row in &probes {
        kv(
            "kmer_probe",
            format!(
                "phase={} k={} outcome={} thread_n/e={}/{} src={} sink={} kbest={} haps={} nonref={}",
                row.phase,
                row.kmer_size,
                row.outcome,
                row.thread_nodes,
                row.thread_edges,
                row.has_ref_source,
                row.has_ref_sink,
                row.kbest_paths,
                row.extracted_haps,
                row.non_ref_haps
            ),
        );
    }

    let c_kmers_25 = snp_kmers(&graph_ref.bases, 25, off, b'C');
    let g_kmers_25 = snp_kmers(&graph_ref.bases, 25, off, b'G');
    kv("c_kmers_25_n", c_kmers_25.len().to_string());
    kv("g_kmers_25_n", g_kmers_25.len().to_string());

    let mut read_c_kmer_in_graph = BTreeMap::new();
    let params25 = threading_params(25, &assembler);
    let thread_graph =
        assembly_graph_from_ref_and_reads_threading(&graph_ref, &graph_reads, &params25)
            .expect("thread k=25");
    kv(
        "rt_thread_k25",
        format!(
            "nodes={} edges={} cycles={}",
            thread_graph.node_count(),
            thread_graph.edge_count(),
            thread_graph.has_cycle()
        ),
    );
    let (c_thread_n, c_thread_sup) = graph_has_kmers(&thread_graph, &c_kmers_25);
    let (g_thread_n, g_thread_sup) = graph_has_kmers(&thread_graph, &g_kmers_25);
    kv(
        "rt_thread_c_kmers",
        format!("nodes={c_thread_n} support_sum={c_thread_sup}"),
    );
    kv(
        "rt_thread_g_kmers",
        format!("nodes={g_thread_n} support_sum={g_thread_sup}"),
    );

    let mut pruned = thread_graph.clone();
    let mut pruning = AssemblyGraphPruningParams::gatk_haplotype_caller_defaults();
    pruning.min_prune_factor = assembler.min_prune_factor;
    pruning.use_adaptive_pruning = assembler.use_adaptive_pruning;
    let pruned_n = pruned.apply_pruning(&pruning);
    kv(
        "rt_pruned_k25",
        format!(
            "removed={pruned_n} nodes={} edges={} cycles={}",
            pruned.node_count(),
            pruned.edge_count(),
            pruned.has_cycle()
        ),
    );
    let (c_prune_n, c_prune_sup) = graph_has_kmers(&pruned, &c_kmers_25);
    let (g_prune_n, g_prune_sup) = graph_has_kmers(&pruned, &g_kmers_25);
    kv(
        "rt_pruned_c_kmers",
        format!("nodes={c_prune_n} support_sum={c_prune_sup}"),
    );
    kv(
        "rt_pruned_g_kmers",
        format!("nodes={g_prune_n} support_sum={g_prune_sup}"),
    );

    let seq_asm = build_threading_graph_for_seq_assembly(
        &graph_ref,
        &graph_reads,
        25,
        &assembler,
        false,
        false,
    )
    .expect("seq asm");
    kv("rt_seq_asm_present", seq_asm.is_some().to_string());
    let mut c_seq_asm_n = 0usize;
    let mut g_seq_asm_n = 0usize;
    let mut seq_pre_c = false;
    let mut seq_post_c = false;
    let mut seq_post_g = false;
    let mut seq_c_edge: Option<(usize, usize, u32, bool)> = None;
    let mut kbest_c_n = 0usize;
    let mut kbest_g_n = 0usize;
    let mut kbest_c_rank: Option<usize> = None;
    let mut kbest_c_score: Option<f64> = None;
    let mut extract_c_n = 0usize;
    let mut seq_nodes = 0usize;
    let mut seq_edges = 0usize;
    let mut cleanup_snaps: Vec<(String, bool, bool, usize, usize)> = Vec::new();
    if let Some(g) = seq_asm.as_ref() {
        kv(
            "rt_seq_asm_k25",
            format!(
                "nodes={} edges={} cycles={}",
                g.node_count(),
                g.edge_count(),
                g.has_cycle()
            ),
        );
        (c_seq_asm_n, _) = graph_has_kmers(g, &c_kmers_25);
        (g_seq_asm_n, _) = graph_has_kmers(g, &g_kmers_25);
        kv("rt_seq_asm_c_kmers", format!("nodes={c_seq_asm_n}"));
        kv("rt_seq_asm_g_kmers", format!("nodes={g_seq_asm_n}"));

        let mut seq = SeqGraph::from_assembly_graph(g);
        seq_pre_c = seq_has_c_vertex(&seq);
        kv(
            "seq_pre_cleanup",
            format!(
                "nodes={} edges={} c_vertex={seq_pre_c} g_vertex={}",
                seq.node_count(),
                seq.edge_count(),
                seq_has_g_vertex(&seq)
            ),
        );
        seq.clean_non_ref_paths();
        let status = seq.traced_cleanup_seq_graph(|stage, g| {
            cleanup_snaps.push((
                stage.to_string(),
                seq_has_c_vertex(g),
                seq_has_g_vertex(g),
                g.node_count(),
                g.edge_count(),
            ));
        });
        kv("seq_cleanup_status", format!("{status:?}"));
        for (stage, c, g_v, n, e) in &cleanup_snaps {
            kv(
                "seq_cleanup_snap",
                format!("stage={stage} c_vertex={c} g_vertex={g_v} nodes={n} edges={e}"),
            );
        }
        seq_post_c = seq_has_c_vertex(&seq);
        seq_post_g = seq_has_g_vertex(&seq);
        seq_nodes = seq.node_count();
        seq_edges = seq.edge_count();
        seq_c_edge = seq_c_edge_support(&seq);
        kv(
            "seq_post_cleanup",
            format!(
                "nodes={seq_nodes} edges={seq_edges} c_vertex={seq_post_c} g_vertex={seq_post_g}"
            ),
        );
        if let Some((from, to, sup, is_ref)) = seq_c_edge {
            kv(
                "seq_c_edge",
                format!("from={from} to={to} support={sup} is_ref={is_ref}"),
            );
        } else {
            kv("seq_c_edge", "absent");
        }
        assert_eq!(status, SeqGraphCleanupStatus::AssembledSomeVariation);

        let src_out: Vec<_> = seq
            .edges()
            .iter()
            .filter(|e| e.from == 0)
            .map(|e| format!("{}:sup={}:ref={}", e.to, e.support, e.is_ref))
            .collect();
        kv("seq_source_out_edges", src_out.join(","));
        kv(
            "seq_vs_java",
            format!("java_nodes=37 java_edges=52 rust_nodes={seq_nodes} rust_edges={seq_edges}"),
        );

        let paths = find_best_haplotypes_seq_graph(&seq, JAVA_KBEST_K).expect("kbest");
        kv("kbest_n", paths.len().to_string());
        let mut n_contain_c_vertex = 0usize;
        for (rank, p) in paths.iter().enumerate() {
            let bases = seq.path_bases_bytes(p.start, &p.edges);
            if bases
                .windows(JAVA_C_VERTEX.len())
                .any(|w| w == JAVA_C_VERTEX.as_bytes())
            {
                n_contain_c_vertex += 1;
                kv(
                    "kbest_contains_c_vertex",
                    format!(
                        "rank={rank} score={:.8} hash={} len={}",
                        p.score,
                        fnv1a64_hex(&bases),
                        bases.len()
                    ),
                );
            }
            match bases.get(off).copied() {
                Some(b'C') => {
                    kbest_c_n += 1;
                    if kbest_c_rank.is_none() {
                        kbest_c_rank = Some(rank);
                        kbest_c_score = Some(p.score);
                    }
                    kv(
                        "kbest_c_path",
                        format!(
                            "rank={rank} score={:.8} isRef={} hash={}",
                            p.score,
                            p.is_reference,
                            fnv1a64_hex(&bases)
                        ),
                    );
                }
                Some(b'G') => kbest_g_n += 1,
                Some(other) => kv("kbest_other_base", format!("{} rank={rank}", other as char)),
                None => {}
            }
            if rank == 116 || rank == 127 {
                kv(
                    &format!("kbest_rank_{rank}"),
                    format!(
                        "hash={} score={:.8} len={} base114={}",
                        fnv1a64_hex(&bases),
                        p.score,
                        bases.len(),
                        bases.get(off).copied().map(|b| b as char).unwrap_or('?')
                    ),
                );
            }
        }
        kv(
            "kbest_base_29455314",
            format!("C={kbest_c_n} G={kbest_g_n} n={}", paths.len()),
        );
        kv(
            "kbest_paths_containing_java_c_vertex",
            n_contain_c_vertex.to_string(),
        );
        kv(
            "kbest_c_rank",
            kbest_c_rank
                .map(|r| r.to_string())
                .unwrap_or_else(|| "absent".into()),
        );
        assert!(
            n_contain_c_vertex >= 1,
            "6R.239: production K=128 takes the C-bearing SeqGraph vertex"
        );
        let mut diag256_rank = None;
        let mut diag256_hash = String::new();
        let mut diag256_score = 0.0f64;
        for &kdiag in &[256usize, 512, 1024] {
            let extra = find_best_haplotypes_seq_graph(&seq, kdiag).expect("kbest extra");
            let mut first_c: Option<(usize, f64, String)> = None;
            let mut c_n = 0usize;
            let mut contain = 0usize;
            for (rank, p) in extra.iter().enumerate() {
                let bases = seq.path_bases_bytes(p.start, &p.edges);
                if bases
                    .windows(JAVA_C_VERTEX.len())
                    .any(|w| w == JAVA_C_VERTEX.as_bytes())
                {
                    contain += 1;
                }
                if bases.get(off) == Some(&b'C') {
                    c_n += 1;
                    if first_c.is_none() {
                        first_c = Some((rank, p.score, fnv1a64_hex(&bases)));
                    }
                }
            }
            kv(
                &format!("kbest_diagnostic_K{kdiag}"),
                format!(
                    "n={} C_at_114={c_n} contain_c_vertex={contain} first={:?}",
                    extra.len(),
                    first_c
                ),
            );
            if kdiag == 256 {
                if let Some((rank, score, hash)) = first_c {
                    diag256_rank = Some(rank);
                    diag256_score = score;
                    diag256_hash = hash;
                }
            }
        }
        kv(
            "diagnostic_k256_carrier",
            format!(
                "rank={} score={:.8} hash={} (production K=128 cutoff is rank<128)",
                diag256_rank.unwrap_or(usize::MAX),
                diag256_score,
                diag256_hash
            ),
        );
        assert_eq!(diag256_rank, Some(116));
        assert_eq!(diag256_hash, JAVA_CARRIER);

        let mut ref_hap = Haplotype::new(graph_ref.bases.as_slice(), true);
        let mut ref_cigar = Cigar::new();
        ref_cigar.push(ref_hap.bases.len(), CigarOperator::Match);
        ref_hap.cigar = Some(ref_cigar);
        let ref_cigar_len = ref_hap.cigar.as_ref().unwrap().reference_length();
        let extracted = extract_haplotypes_from_seq_kbest_paths(
            &paths,
            &seq,
            25,
            &ref_hap,
            ref_cigar_len,
            &assembler.haplotype_to_reference_sw,
        )
        .unwrap_or_default();
        extract_c_n = extracted
            .iter()
            .filter(|h| h.bases.get(off) == Some(&b'C'))
            .count();
        kv(
            "extract_haps",
            format!("n={} C_at_offset={extract_c_n}", extracted.len()),
        );
        for h in &extracted {
            if h.bases.get(off) == Some(&b'C') {
                kv(
                    "extract_c_hap",
                    format!(
                        "hash={} len={} score={} isRef={}",
                        fnv1a64_hex(&h.bases),
                        h.bases.len(),
                        h.score,
                        h.is_reference
                    ),
                );
            }
        }
    }

    for (qname, flag) in C_READS {
        let fin = assembled
            .finalized_reads
            .iter()
            .find(|r| r.qname() == qname.as_bytes() && r.flags() == *flag)
            .expect("finalized C read");
        let ck = read_c_kmers(fin, 25, EVENT_POS);
        let in_thread = ck
            .iter()
            .filter(|k| {
                thread_graph
                    .nodes()
                    .iter()
                    .any(|n| n.kmer.as_ref() == k.as_slice())
            })
            .count();
        let in_pruned = ck
            .iter()
            .filter(|k| {
                pruned
                    .nodes()
                    .iter()
                    .any(|n| n.kmer.as_ref() == k.as_slice())
            })
            .count();
        let in_seq_asm = seq_asm.as_ref().map(|g| {
            ck.iter()
                .filter(|k| g.nodes().iter().any(|n| n.kmer.as_ref() == k.as_slice()))
                .count()
        });
        let short = qname.rsplit(':').nth(1).unwrap_or(qname);
        kv(
            &format!("c_read_graph_{short}"),
            format!(
                "c_kmers={} thread={} pruned={} seq_asm={:?}",
                ck.len(),
                in_thread,
                in_pruned,
                in_seq_asm
            ),
        );
        read_c_kmer_in_graph.insert(*qname, (ck.len(), in_thread, in_pruned));
    }

    let first = if c_thread_n == 0 {
        "ASSEMBLY_GRAPH_CONSTRUCTION_DIVERGENCE"
    } else if c_prune_n == 0 {
        "ASSEMBLY_GRAPH_PRUNING_DIVERGENCE"
    } else if c_seq_asm_n == 0 {
        "ASSEMBLY_DANGLING_RECOVERY_DIVERGENCE"
    } else if !seq_post_c {
        "ASSEMBLY_GRAPH_CONSTRUCTION_DIVERGENCE"
    } else if kbest_c_n == 0 {
        "ASSEMBLY_KBEST_DIVERGENCE"
    } else if extract_c_n == 0 {
        "HAPLOTYPE_MATERIALIZATION_DIVERGENCE"
    } else if hap_base_counts.get(&'C').copied().unwrap_or(0) == 0 {
        "HAPLOTYPE_DEDUPLICATION_DIVERGENCE"
    } else {
        "NO_DIVERGENCE_AT_THIS_BOUNDARY"
    };
    kv("classification", first);
    kv(
        "first_divergent_operation",
        format!(
            "thread_c={c_thread_n} prune_c={c_prune_n} seq_asm_c={c_seq_asm_n} seq_pre_c={seq_pre_c} seq_post_c={seq_post_c} kbest_c={kbest_c_n} extract_c={extract_c_n} assemble_c={}",
            hap_base_counts.get(&'C').copied().unwrap_or(0)
        ),
    );
    kv(
        "why_eventmap_downstream",
        "EventMap G>C is CIGAR M-mismatch on the Java carrier; Rust never materializes a C-bearing haplotype, so EventMap cannot emit the SNP",
    );
    kv(
        "why_clipping_downstream",
        "trim pad 29455355 vs Java 29455294 depends on leftmost EventMap SNP; that event is missing because the carrier haplotype is missing",
    );
    kv(
        "why_pairhmm_downstream",
        "PairHMM sequences are hard-clipped to the trim interval derived from EventMap; assembly of the carrier is earlier",
    );
    kv("java_kbest_rank", JAVA_CARRIER_RANK.to_string());
    kv(
        "rust_kbest_rank",
        kbest_c_rank
            .map(|r| r.to_string())
            .unwrap_or_else(|| "absent".into()),
    );
    kv(
        "rust_kbest_score",
        kbest_c_score
            .map(|s| format!("{s:.8}"))
            .unwrap_or_else(|| "absent".into()),
    );
    kv("java_kbest_score", JAVA_CARRIER_SCORE);

    assert_eq!(first, "NO_DIVERGENCE_AT_THIS_BOUNDARY");
    assert!(c_thread_n > 0, "C-bearing k-mers exist before prune");
    assert!(c_prune_n > 0, "C-bearing k-mers survive prune");
    assert!(
        c_seq_asm_n > 0,
        "C-bearing k-mers survive dangling/remove_paths"
    );
    assert!(seq_post_c, "C-bearing SeqGraph vertex exists after cleanup");
    assert!(kbest_c_n >= 1, "6R.239: K=128 selects the C-bearing path");
    assert_eq!(
        seq_c_edge.map(|e| e.2),
        Some(JAVA_C_EDGE_MULT),
        "C-bearing SeqGraph edge multiplicity must match Java 3"
    );
    let _ = (
        seq_c_edge,
        kbest_g_n,
        seq_post_g,
        read_c_kmer_in_graph,
        seq_nodes,
        seq_edges,
        g_thread_n,
        g_prune_n,
        g_seq_asm_n,
        c_thread_sup,
        c_prune_sup,
        g_thread_sup,
        g_prune_sup,
    );
}
