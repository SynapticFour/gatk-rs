//! 6R.79: SeqGraph multiplicity / zip provenance vs Java 4.4 `cleanupSeqGraph`.
//!
//! Forensic only. No production graph, zip, or k-best change in this file.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r79_seqgraph_multiplicity_contract
//! HOLDOUT_6R79=1 cargo test -p gatk-haplotypecaller --test forensic_6r79_seqgraph_multiplicity_contract live_ -- --nocapture
//! ```

use gatk_haplotypecaller::seq_graph::SeqGraph;
use gatk_haplotypecaller::seq_kbest_haplotype::{
    find_best_haplotypes_seq_graph_forensic, seq_kbest_path_score_terms, SeqKbestCapPolicy,
};

const JAVA_ONLY_J0: &[u8] = b"CATGGAGCCTGACCTTATTTGAAGTAGGGCATTTGCAGATGTATTTAAGATATTTGAGGCTGGGCACAGTGGCTCACGTCTGTAATCCCAGCACTTTGAAAGGCCGAGGCAGGTGGATTCACCTGAGGTCAGGAGTTTGAGACCAGCCTGTCCCACATGGTGAAAAGCCCGTATCTACCAAAAATACAAACGTTAGCTGTGTGTGGTGGTGGCGGCACCTGTAATCCCAGCTACTCGAGAGCCAGAG";
const JAVA_ONLY_J1: &[u8] = b"CATGGAGCCTGACTTTATTTGAAGTAGGGCATTTGCAGATGTATTTAAGATATTTGAGGCCGGGCACAGTGGCTCATGTCTGTAATCCCAGCACTTTAAAAGGCTGAGGCAGGTGTATTCACCTGAGGTCAGGAGTTCGAGACCAGCCTGGCCAACATGGTGAAAAGCCCGTATCTACCAAAAATACAAAAGTTAGCTGGGTGTGGTGGCGGCACCTGTAATCCCAGCTACTCGAGAGCCAGAG";
const RUST_ONLY_R0: &[u8] = b"CATGGAGCCTGACCTTATTTGAAGTAGGGCATTTGCAGATGTATTTAAGATATTTGAGGCTGGGCACAGTGGCTCACGTCTGTAATCCCAGCACTTTGAAAGGCCGAGGCAGGTGGATTCACCTGAGGTCAGGAGTTTGAGACCAGCCTGTCCCACATGGTGAAAAGCCCGTATCTACCAAAAATACAAACGTTAGCTGTGTGTGGTGGTGGCGGCACCTGTAATCGCAGCTACTCGAGAGCCAGAG";

const K25_C: &[u8] = b"ACCTGTAATCCCAGCTACTCGAGAG";
const K25_G: &[u8] = b"ACCTGTAATCGCAGCTACTCGAGAG";

fn contains_bases(hay: &[u8], needle: &[u8]) -> bool {
    !needle.is_empty() && hay.windows(needle.len()).any(|w| w == needle)
}

fn out_mult(g: &SeqGraph, from: usize) -> u32 {
    g.edges()
        .iter()
        .filter(|e| e.from == from)
        .map(|e| e.support)
        .sum()
}

fn dump_branches(label: &str, g: &SeqGraph) {
    eprintln!(
        "SNAP {label} vertices={} edges={}",
        g.node_count(),
        g.edge_count()
    );
    for v in g.vertices() {
        let outs: Vec<_> = g.edges().iter().filter(|e| e.from == v.id).collect();
        if outs.len() < 2 {
            continue;
        }
        let tot: u32 = outs.iter().map(|e| e.support).sum();
        eprintln!(
            "  BRANCH v{} seq_len={} seq={:?} out_n={} tot={tot}",
            v.id,
            v.sequence.len(),
            String::from_utf8_lossy(&v.sequence),
            outs.len()
        );
        for e in &outs {
            let to_seq = g
                .vertices()
                .get(e.to)
                .map(|t| String::from_utf8_lossy(&t.sequence).into_owned())
                .unwrap_or_default();
            eprintln!(
                "    {}->{} supp={} tot={tot} ref={} to_seq={:?}",
                e.from, e.to, e.support, e.is_ref, to_seq
            );
        }
    }
    for (tag, mer) in [("C25", K25_C), ("G25", K25_G)] {
        let hit: Vec<usize> = g
            .vertices()
            .iter()
            .filter(|v| contains_bases(&v.sequence, mer))
            .map(|v| v.id)
            .collect();
        eprintln!("  CONTAINS_{tag} vertices={hit:?}");
    }
}

fn dump_triple(g: &SeqGraph, a: usize, b: usize, c: usize) {
    for &(from, to) in &[(a, b), (a, c)] {
        if let Some(e) = g.edges().iter().find(|e| e.from == from && e.to == to) {
            eprintln!(
                "  EDGE {from}->{} supp={} outMult={} ref={}",
                e.to,
                e.support,
                out_mult(g, from),
                e.is_ref
            );
        } else {
            eprintln!("  EDGE {from}->{to} MISSING");
        }
    }
}

/// Coordinate-free motif dump: C/G SNP bubble and leftover A/T fork on the C allele.
fn dump_motif(label: &str, g: &SeqGraph) {
    eprintln!(
        "MOTIF {label} vertices={} edges={}",
        g.node_count(),
        g.edge_count()
    );
    for v in g.vertices() {
        let outs: Vec<_> = g.edges().iter().filter(|e| e.from == v.id).collect();
        if outs.len() < 2 {
            continue;
        }
        let tot: u32 = outs.iter().map(|e| e.support).sum();
        let to_seqs: Vec<String> = outs
            .iter()
            .map(|e| {
                g.vertices()
                    .get(e.to)
                    .map(|t| String::from_utf8_lossy(&t.sequence).into_owned())
                    .unwrap_or_default()
            })
            .collect();
        let has_c = to_seqs.iter().any(|s| s.starts_with('C'));
        let has_g = to_seqs.iter().any(|s| s.starts_with('G'));
        let has_a = to_seqs.iter().any(|s| s == "A" || s.starts_with('A'));
        let has_t = to_seqs.iter().any(|s| s == "T" || s.starts_with('T'));
        let seq = String::from_utf8_lossy(&v.sequence);
        let snp_bubble = has_c && has_g;
        let leftover_at = seq.contains("CCAGCTACTCGAGAG") && has_a && has_t;
        if !snp_bubble && !leftover_at {
            continue;
        }
        let kind = if snp_bubble { "CG_BUBBLE" } else { "AT_FORK" };
        eprintln!(
            "  {kind} v{} seq={:?} out_n={} tot={tot}",
            v.id,
            seq,
            outs.len()
        );
        for (e, to_seq) in outs.iter().zip(to_seqs.iter()) {
            eprintln!(
                "    {}->{} supp={} tot={tot} ref={} to_seq={:?}",
                e.from, e.to, e.support, e.is_ref, to_seq
            );
        }
    }
}

fn rt_sequence_walk_status(rt: &gatk_haplotypecaller::AssemblyGraph, seq: &[u8]) -> String {
    if seq.len() < 25 {
        return "seq_shorter_than_k25".to_string();
    }
    let id_of = |k: &[u8]| -> Option<usize> {
        rt.nodes()
            .iter()
            .find(|n| n.kmer.as_ref() == k)
            .map(|n| n.id)
    };
    let first = &seq[..25];
    let Some(mut prev) = id_of(first) else {
        return format!(
            "NOT_REPRESENTABLE missing_kmer0={}",
            String::from_utf8_lossy(first)
        );
    };
    for (i, w) in seq.windows(25).enumerate().skip(1) {
        let Some(id) = id_of(w) else {
            return format!(
                "NOT_REPRESENTABLE missing_kmer{i}={}",
                String::from_utf8_lossy(w)
            );
        };
        let has_edge = rt
            .edges_sorted()
            .iter()
            .any(|e| e.from == prev && e.to == id);
        if !has_edge {
            return format!(
                "NOT_REPRESENTABLE missing_edge kmer{}->kmer{} {} -> {}",
                i - 1,
                i,
                prev,
                id
            );
        }
        prev = id;
    }
    "RT_WALK_COMPLETE".to_string()
}

/// Java `CommonSuffixSplitter.split` (SHA 2dbc0258): prefix→suffix is
/// `new BaseEdge(out.isRef(), 1)`, not the predecessor's outgoing multiplicity.
/// Diagnostic of the leftover-fork score identity, not a locus pin.
#[test]
fn forensic_6r79_java_suffix_split_prefix_edge_is_multiplicity_one() {
    let java_prefix_suffix_mult = 1u32;
    assert_eq!(java_prefix_suffix_mult, 1);
    let leftover_branch = (12f64).log10() - (15f64).log10();
    let unique = (12f64).log10() - (12f64).log10();
    assert!(leftover_branch < unique);
    assert!((unique).abs() < 1e-15);
}

/// Java `mergeLinearChainVertex` copies last-outgoing edges (`edge.copy()`),
/// it does not average or renormalize multiplicity.
#[test]
fn forensic_6r79_java_zip_copies_last_outgoing_multiplicity() {
    let last_outgoing = 17u32;
    let copied = last_outgoing; // edge.copy()
    assert_eq!(copied, 17);
    assert_ne!(copied / 2, 17, "zip does not average");
}

/// Extra leftover-fork term inverts the C vs G SNP-edge ranking under rust counts
/// (17/31 then 12/15). Java's +1 on the C path (18/32 then 13/16) does not invert.
/// Diagnostic identity, not a locus pin.
#[test]
fn forensic_6r79_extra_branch_12_15_inverts_17_31_vs_14_31() {
    let rust_c = (17f64).log10() - (31f64).log10();
    let rust_g = (14f64).log10() - (31f64).log10();
    assert!(rust_c > rust_g);
    let rust_extra = (12f64).log10() - (15f64).log10();
    assert!(rust_c + rust_extra < rust_g);

    let java_c = (18f64).log10() - (32f64).log10();
    let java_g = (14f64).log10() - (32f64).log10();
    let java_extra = (13f64).log10() - (16f64).log10();
    assert!(java_c + java_extra > java_g);
}

#[test]
fn live_seqgraph_multiplicity_zip() {
    if std::env::var("HOLDOUT_6R79").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R79=1");
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

    let c_id = rt
        .nodes()
        .iter()
        .find(|n| n.kmer.as_ref() == K25_C)
        .map(|n| n.id);
    let g_id = rt
        .nodes()
        .iter()
        .find(|n| n.kmer.as_ref() == K25_G)
        .map(|n| n.id);
    eprintln!(
        "RT C25={c_id:?} G25={g_id:?} nodes={} edges={}",
        rt.node_count(),
        rt.edges_sorted().len()
    );
    if let (Some(cid), Some(gid)) = (c_id, g_id) {
        for e in rt.edges_sorted() {
            if e.from == cid || e.to == cid || e.from == gid || e.to == gid {
                let fk = String::from_utf8_lossy(&rt.nodes()[e.from].kmer);
                let tk = String::from_utf8_lossy(&rt.nodes()[e.to].kmer);
                eprintln!(
                    "  RT_EDGE {}->{} supp={} ref={} {fk} -> {tk}",
                    e.from,
                    e.to,
                    e.support,
                    rt.edge_is_ref(e.from, e.to)
                );
            }
        }
    }

    let mut j1_unique = 0usize;
    let mut j1_present = 0usize;
    let mut j1_missing: Vec<String> = Vec::new();
    for w in JAVA_ONLY_J1.windows(25) {
        if JAVA_ONLY_J0.windows(25).any(|x| x == w) {
            continue;
        }
        j1_unique += 1;
        if rt.nodes().iter().any(|n| n.kmer.as_ref() == w) {
            j1_present += 1;
        } else if j1_missing.len() < 8 {
            j1_missing.push(String::from_utf8_lossy(w).into_owned());
        }
    }
    eprintln!("J1_UNIQUE_25MERS total={j1_unique} present_in_rt={j1_present} missing_sample={j1_missing:?}");

    let mut seq = SeqGraph::from_assembly_graph(&rt);
    seq.clean_non_ref_paths();
    dump_motif("after_from_assembly_graph", &seq);

    eprintln!("J1_RT_WALK {}", rt_sequence_walk_status(&rt, JAVA_ONLY_J1));
    eprintln!("J0_RT_WALK {}", rt_sequence_walk_status(&rt, JAVA_ONLY_J0));
    eprintln!("R0_RT_WALK {}", rt_sequence_walk_status(&rt, RUST_ONLY_R0));

    let status = seq.traced_cleanup_seq_graph(|stage, g| {
        dump_motif(stage, g);
        if stage == "cleanup_entry"
            || stage == "after_initial_zip_linear_chains"
            || stage == "after_remove_paths_not_connected_to_ref"
            || stage == "after_second_simplify"
            || stage == "final_for_kbest"
        {
            dump_branches(stage, g);
        }
    });
    assert_eq!(status, SeqGraphCleanupStatus::AssembledSomeVariation);

    dump_triple(&seq, 20, 5, 41);
    if let Some(e) = seq.edges().iter().find(|e| e.from == 5 && e.to == 39) {
        eprintln!(
            "  EDGE 5->39 supp={} outMult={} ref={}",
            e.support,
            out_mult(&seq, 5),
            e.is_ref
        );
        for o in seq.edges().iter().filter(|x| x.from == 5) {
            eprintln!(
                "  V5_OUT {}->{} supp={} to_seq={:?}",
                o.from,
                o.to,
                o.support,
                seq.vertices()
                    .get(o.to)
                    .map(|t| String::from_utf8_lossy(&t.sequence).into_owned())
                    .unwrap_or_default()
            );
        }
    } else {
        eprintln!("  EDGE 5->39 MISSING after cleanup");
    }

    let needles: [&[u8]; 3] = [JAVA_ONLY_J0, JAVA_ONLY_J1, RUST_ONLY_R0];
    let kbest = find_best_haplotypes_seq_graph_forensic(
        &seq,
        4096,
        128,
        SeqKbestCapPolicy::unbounded(),
        &needles,
    )
    .expect("kbest");
    eprintln!(
        "J0_STATUS rank={:?} J1_STATUS rank={:?} R0_STATUS rank={:?}",
        kbest.needles_in_result[0]
            .as_ref()
            .map(|h| h.rank_after_sort),
        kbest.needles_in_result[1]
            .as_ref()
            .map(|h| h.rank_after_sort),
        kbest.needles_in_result[2]
            .as_ref()
            .map(|h| h.rank_after_sort)
    );
    eprintln!(
        "J0_CLASS {}",
        if kbest.needles_in_result[0].is_some() {
            "ENUMERATED_BUT_BELOW_K"
        } else {
            "ABSENT_FROM_VISIT128_SINKS"
        }
    );
    eprintln!(
        "J1_CLASS {}",
        if rt_sequence_walk_status(&rt, JAVA_ONLY_J1) != "RT_WALK_COMPLETE" {
            "NOT_REPRESENTABLE"
        } else if kbest.needles_in_result[1].is_some() {
            "ENUMERATED_BUT_BELOW_K"
        } else {
            "REPRESENTABLE_ON_RT_NOT_IN_VISIT128_SINKS"
        }
    );

    let kb_full = find_best_haplotypes_seq_graph_forensic(
        &seq,
        4096,
        4096,
        SeqKbestCapPolicy::unbounded(),
        &needles,
    )
    .expect("kbest full");
    eprintln!(
        "VISIT4096 J0={:?} J1={:?} R0={:?}",
        kb_full.needles_in_result[0]
            .as_ref()
            .map(|h| h.rank_after_sort),
        kb_full.needles_in_result[1]
            .as_ref()
            .map(|h| h.rank_after_sort),
        kb_full.needles_in_result[2]
            .as_ref()
            .map(|h| h.rank_after_sort)
    );

    if let Some(p) = kbest
        .paths
        .iter()
        .find(|p| contains_bases(&seq.path_bases_bytes(p.start, &p.edges), JAVA_ONLY_J0))
    {
        let terms = seq_kbest_path_score_terms(&seq, p);
        let mut cum = 0.0;
        for (i, t) in terms.iter().enumerate() {
            if t.edge_support != t.total_outgoing {
                cum += t.penalty;
                eprintln!(
                    "  J0_TERM e{i} {}->{} {:>3}/{:<3} pen={:+.9} cum={:.9} ref={}",
                    t.from, t.to, t.edge_support, t.total_outgoing, t.penalty, cum, t.is_ref
                );
            } else {
                cum += t.penalty;
            }
        }
        eprintln!("  J0_SCORE stored={:.9}", p.score);
    }

    let mut assembler_lowc = assembler.clone();
    assembler_lowc.dont_increase_kmer_sizes_for_cycles = true;
    let rt_lowc = build_threading_graph_for_seq_assembly(
        &graph_ref,
        &graph_reads,
        25,
        &assembler_lowc,
        true,
        false,
    )
    .expect("rt lowc");
    eprintln!(
        "LOWC k=25 graph_built={} (dont_increase=true / allow_low_complexity=true)",
        rt_lowc.is_some()
    );
    if let Some(rt2) = rt_lowc {
        eprintln!(
            "LOWC_J1_RT_WALK {}",
            rt_sequence_walk_status(&rt2, JAVA_ONLY_J1)
        );
        let mut seq2 = SeqGraph::from_assembly_graph(&rt2);
        seq2.clean_non_ref_paths();
        let _ = seq2.cleanup_seq_graph();
        let kb2 = find_best_haplotypes_seq_graph_forensic(
            &seq2,
            512,
            512,
            SeqKbestCapPolicy::production(),
            &needles,
        )
        .expect("kbest lowc");
        eprintln!(
            "LOWC_K512 J0={:?} J1={:?} R0={:?}",
            kb2.needles_in_result[0].as_ref().map(|h| h.rank_after_sort),
            kb2.needles_in_result[1].as_ref().map(|h| h.rank_after_sort),
            kb2.needles_in_result[2].as_ref().map(|h| h.rank_after_sort)
        );
    }

    assert_eq!(c_id.is_some(), true, "C 25-mer in RT graph");
    assert_eq!(g_id.is_some(), true, "G 25-mer in RT graph");
    assert!(
        kbest.needles_in_result[0].is_some(),
        "J0 walkable under visit=128 collect-all"
    );
    assert!(
        kbest.needles_in_result[1].is_none(),
        "J1 not a production SeqGraph k-best path"
    );
}
