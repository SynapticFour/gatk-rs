//! 6R.52 re-measure of the 6R.51 ActiveFull covering 20:29455379 G/A.
//!
//! After the SeqGraph k-best heap-order fix, production K=128 must include the
//! A-bearing path. Downstream stages are recorded, not patched.
//!
//! Skipped unless `HOLDOUT_6R51=1`.
//!
//! ```text
//! HOLDOUT_6R51=1 cargo test -p gatk-haplotypecaller --test holdout_6r51_chr20_tiny_test -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::assembly_region_finalize::{
    assembly_reference_read, create_graph_reference_read, records_to_assembly_reads,
};
use gatk_haplotypecaller::event_map::variation_events_for_haplotype;
use gatk_haplotypecaller::hc_allele_mapping::hap_base_at_ref_locus;
use gatk_haplotypecaller::hc_genotyping_engine::DEFAULT_STAND_EMIT_CONFIDENCE;
use gatk_haplotypecaller::read_threading_assembler::{
    build_threading_graph_for_seq_assembly, extract_rt_haplotypes_after_remove_paths,
    extract_rt_haplotypes_before_remove_paths, DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH,
};
use gatk_haplotypecaller::seq_graph::SeqGraph;
use gatk_haplotypecaller::seq_kbest_haplotype::find_best_haplotypes_seq_graph;
use gatk_haplotypecaller::{
    assemble_from_ref_and_reads, assemble_reads_with_finalized,
    assembly_graph_from_ref_and_reads_threading, calculate_haplotype_cigar_for_assembly,
    call_disposition, flatten_assembly_regions, probe_seq_graph_kmer_attempts,
    query_index_at_reference_position, reference_has_non_unique_kmers,
    traverse_assembly_region_walker, try_emit_call_region_variants, AssemblyGraph,
    AssemblyGraphParams, AssemblyRegionCallDisposition, CallRegionArgs, Haplotype,
    HaplotypeCallerEngine, KmerSize, ReadFilterParams, SwParameters, WalkerTraversalConfig,
};
use rust_htslib::bam::record::CigarString;
use serde_json::{json, Value};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET_POS: u64 = 29_455_379;
const TARGET_REF: u8 = b'G';
const TARGET_ALT: u8 = b'A';
const NEARBY_POS: u64 = 29_455_389;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn cigar_str(h: &Haplotype) -> String {
    h.cigar
        .as_ref()
        .map(|c| c.to_gatk_string())
        .unwrap_or_else(|| "-".to_string())
}

fn hap_row(h: &Haplotype, pad: u64, loc: u64) -> Value {
    let base = hap_base_at_ref_locus(h, pad, loc)
        .map(|b| (b as char).to_string())
        .unwrap_or_else(|| ".".to_string());
    json!({
        "is_ref": h.is_reference,
        "len": h.bases.len(),
        "k": h.kmer_size,
        "score": h.score,
        "align0": h.alignment_start_hap_wrt_ref,
        "cigar": cigar_str(h),
        "loc": h.genome_loc.map(|g| [g.start.get(), g.end.get()]),
        "base_at_target": base,
        "eq_target_alt": hap_base_at_ref_locus(h, pad, loc) == Some(TARGET_ALT),
        "eq_target_ref": hap_base_at_ref_locus(h, pad, loc) == Some(TARGET_REF),
    })
}

fn hap_has_alt(haps: &[Haplotype], pad: u64, loc: u64) -> bool {
    haps.iter()
        .any(|h| hap_base_at_ref_locus(h, pad, loc) == Some(TARGET_ALT))
}

fn event_hits(events: &[gatk_haplotypecaller::event_map::VariationEvent], loc: u64) -> Vec<Value> {
    events
        .iter()
        .filter(|e| e.start_1based.get() == loc)
        .map(|e| {
            json!({
                "pos": e.start_1based.get(),
                "ref": e.ref_allele,
                "alt": e.alt_allele,
            })
        })
        .collect()
}

fn pileup_at(reads: &[rust_htslib::bam::Record], loc_1based: u64) -> BTreeMap<char, usize> {
    let pos0 = (loc_1based - 1) as i64;
    let mut counts: BTreeMap<char, usize> = BTreeMap::new();
    for rec in reads {
        if rec.is_unmapped() {
            continue;
        }
        let cigar = CigarString(rec.cigar().iter().copied().collect());
        let Some(qi) = query_index_at_reference_position(rec.pos(), &cigar, pos0) else {
            continue;
        };
        let seq = rec.seq();
        if qi >= seq.len() {
            continue;
        }
        let b = (seq.as_bytes()[qi] as char).to_ascii_uppercase();
        *counts.entry(b).or_insert(0) += 1;
    }
    counts
}

fn kmer_windows_spanning(seq: &[u8], offset: usize, k: usize) -> Vec<Vec<u8>> {
    if seq.len() < k || offset >= seq.len() {
        return Vec::new();
    }
    let start_min = offset.saturating_sub(k - 1);
    let start_max = offset.min(seq.len() - k);
    (start_min..=start_max)
        .map(|i| seq[i..i + k].to_vec())
        .collect()
}

fn count_reads_containing_kmer(reads: &[gatk_haplotypecaller::AssemblyRead], kmer: &[u8]) -> usize {
    reads
        .iter()
        .filter(|r| r.bases.windows(kmer.len()).any(|w| w == kmer))
        .count()
}

fn alt_only_windows(ref_seq: &[u8], alt_seq: &[u8], offset: usize, k: usize) -> Vec<Vec<u8>> {
    let ref_w = kmer_windows_spanning(ref_seq, offset, k);
    kmer_windows_spanning(alt_seq, offset, k)
        .into_iter()
        .filter(|km| !ref_w.iter().any(|r| r == km))
        .collect()
}

fn seq_contains_any_kmer(seq: &[u8], kmers: &[Vec<u8>]) -> bool {
    kmers
        .iter()
        .any(|km| seq.windows(km.len()).any(|w| w == km.as_slice()))
}

fn graph_alt_kmer_stats(graph: &AssemblyGraph, alt_windows: &[Vec<u8>]) -> Value {
    let mut nodes_hit = 0usize;
    let mut max_node_support = 0u32;
    for n in graph.nodes() {
        if alt_windows.iter().any(|km| &*n.kmer == km.as_slice()) {
            nodes_hit += 1;
            max_node_support = max_node_support.max(n.support);
        }
    }
    let alt_ids: std::collections::HashSet<usize> = graph
        .nodes()
        .iter()
        .enumerate()
        .filter(|(_, n)| alt_windows.iter().any(|km| &*n.kmer == km.as_slice()))
        .map(|(i, _)| i)
        .collect();
    let mut max_edge_support = 0u32;
    let mut edges_touching = 0usize;
    for e in graph.edges_sorted() {
        if alt_ids.contains(&e.from) || alt_ids.contains(&e.to) {
            edges_touching += 1;
            max_edge_support = max_edge_support.max(e.support);
        }
    }
    json!({
        "nodes": graph.node_count(),
        "edges": graph.edge_count(),
        "alt_nodes": nodes_hit,
        "max_alt_node_support": max_node_support,
        "alt_touching_edges": edges_touching,
        "max_alt_edge_support": max_edge_support,
    })
}

fn hap_compact(
    haps: &[Haplotype],
    pad: u64,
    loc: u64,
    alt_windows: &[Vec<u8>],
    ref_len: usize,
) -> Value {
    let mut cigars: BTreeMap<String, usize> = BTreeMap::new();
    let mut mapped_alt = 0usize;
    let mut substr_alt = 0usize;
    let mut eq_len_idx_alt = 0usize;
    for h in haps {
        *cigars.entry(cigar_str(h)).or_insert(0) += 1;
        if hap_base_at_ref_locus(h, pad, loc) == Some(TARGET_ALT) {
            mapped_alt += 1;
        }
        if seq_contains_any_kmer(&h.bases, alt_windows) {
            substr_alt += 1;
        }
        if h.bases.len() == ref_len {
            let off = loc.saturating_sub(pad) as usize;
            if h.bases.get(off) == Some(&TARGET_ALT) {
                eq_len_idx_alt += 1;
            }
        }
    }
    json!({
        "n": haps.len(),
        "mapped_alt": mapped_alt,
        "substr_alt_kmer": substr_alt,
        "eq460_index_alt": eq_len_idx_alt,
        "cigars": cigars,
    })
}

fn hap_base_from_seq(seq: &[u8], ref_seq: &[u8], pad: u64, loc: u64) -> Option<u8> {
    if seq.len() == ref_seq.len() {
        let off = loc.saturating_sub(pad) as usize;
        return seq.get(off).copied();
    }
    let sw = SwParameters::gatk_haplotype_to_reference();
    let cigar = calculate_haplotype_cigar_for_assembly(ref_seq, seq, ref_seq.len(), &sw)?;
    let mut h = Haplotype::new(seq, false);
    h.cigar = Some(cigar);
    hap_base_at_ref_locus(&h, pad, loc)
}

#[test]
fn holdout_6r51_chr20_tiny_covering_activefull() {
    if std::env::var("HOLDOUT_6R51").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R51=1");
        return;
    }
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    assert!(ref_fasta.is_file(), "missing {}", ref_fasta.display());
    assert!(bam.is_file(), "missing {}", bam.display());

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
            ) && r.start.get() <= TARGET_POS
                && r.end.get() >= TARGET_POS
        })
        .collect();
    assert!(!covering.is_empty(), "no ActiveFull covers {TARGET_POS}");
    let region = covering[0];
    let args = CallRegionArgs::strict_java();
    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);

    let raw_reads: Vec<_> = region.reads.iter().map(|r| r.as_ref().clone()).collect();
    let pileup_raw = pileup_at(&raw_reads, TARGET_POS);
    let pileup_raw_near = pileup_at(&raw_reads, NEARBY_POS);

    let mut owned = region.clone();
    let assembled =
        assemble_reads_with_finalized(&mut owned, &dict, &mut ref_cache, &args.assemble)
            .expect("assemble");
    let pad = assembled.assembly.padded_reference_start_1based();
    let k_used = assembled.assembly.kmer_size_for_dump();
    let untrimmed = &assembled.assembly;
    let ref_hap = untrimmed
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .expect("ref hap");
    let ref_bytes = untrimmed.apply_bases_shared();

    let mut per_hap_events = Vec::new();
    for (i, h) in untrimmed.haplotypes.iter().enumerate() {
        let ev =
            variation_events_for_haplotype(h, ref_hap, ref_bytes.as_ref(), pad, 0, &region.contig);
        per_hap_events.push(json!({
            "i": i,
            "target_events": event_hits(&ev, TARGET_POS),
            "nearby_events": event_hits(&ev, NEARBY_POS),
            "n_events": ev.len(),
        }));
    }

    let padded_ref = assembly_reference_read(&dict, &mut ref_cache, region).expect("pad ref");
    let graph_ref = create_graph_reference_read(&padded_ref, region, &dict);
    let graph_reads = records_to_assembly_reads(&assembled.finalized_reads);
    let pileup_final = pileup_at(&assembled.finalized_reads, TARGET_POS);
    let assembler = args.assemble.assembler.clone();
    let off = TARGET_POS.saturating_sub(pad) as usize;
    assert!(off < graph_ref.bases.len(), "target outside graph ref");
    let mut alt_ref = graph_ref.bases.clone();
    alt_ref[off] = TARGET_ALT;
    let alt10 = alt_only_windows(&graph_ref.bases, &alt_ref, off, 10);
    let alt25 = alt_only_windows(&graph_ref.bases, &alt_ref, off, 25);

    let probe = probe_seq_graph_kmer_attempts(&graph_ref, &graph_reads, &assembler).expect("probe");
    let probe_json: Vec<Value> = probe
        .iter()
        .map(|r| {
            json!({
                "phase": r.phase,
                "k": r.kmer_size,
                "allow_lc": r.allow_low_complexity,
                "allow_nu": r.allow_non_unique_ref,
                "outcome": r.outcome,
                "thread_nodes": r.thread_nodes,
                "thread_edges": r.thread_edges,
                "cleanup": r.cleanup_status,
                "kbest_paths": r.kbest_paths,
                "extracted_haps": r.extracted_haps,
                "non_ref_haps": r.non_ref_haps,
                "path_bases_len": r.path_bases_len,
                "path_eq_ref": r.path_eq_ref_bases,
            })
        })
        .collect();

    let mut rt_rows = Vec::new();
    for k in [10usize, 25, 35] {
        let unique = !reference_has_non_unique_kmers(&graph_ref, k);
        let before = extract_rt_haplotypes_before_remove_paths(
            &graph_ref,
            &graph_reads,
            &assembler,
            k,
            false,
            false,
        )
        .unwrap_or_default();
        let after = extract_rt_haplotypes_after_remove_paths(
            &graph_ref,
            &graph_reads,
            &assembler,
            k,
            false,
            false,
        )
        .unwrap_or_default();
        rt_rows.push(json!({
            "k": k,
            "ref_unique": unique,
            "before": hap_compact(&before, pad, TARGET_POS, if k == 10 { &alt10 } else { &alt25 }, graph_ref.bases.len()),
            "after": hap_compact(&after, pad, TARGET_POS, if k == 10 { &alt10 } else { &alt25 }, graph_ref.bases.len()),
        }));
    }

    let seq_raw = assemble_from_ref_and_reads(&graph_ref, &graph_reads, &assembler).expect("seq");
    let seq_has_alt = hap_has_alt(&seq_raw.haplotypes, pad, TARGET_POS);

    let mut params10 = AssemblyGraphParams::default();
    params10.kmer_size = KmerSize::try_new(10).unwrap();
    params10.min_base_quality = assembler.min_base_quality;
    let mut params25 = AssemblyGraphParams::default();
    params25.kmer_size = KmerSize::try_new(25).unwrap();
    params25.min_base_quality = assembler.min_base_quality;
    let raw10 = assembly_graph_from_ref_and_reads_threading(&graph_ref, &graph_reads, &params10)
        .expect("raw10");
    let raw25 = assembly_graph_from_ref_and_reads_threading(&graph_ref, &graph_reads, &params25)
        .expect("raw25");
    let pruned10 = build_threading_graph_for_seq_assembly(
        &graph_ref,
        &graph_reads,
        10,
        &assembler,
        false,
        false,
    )
    .expect("pruned10");
    let pruned25 = build_threading_graph_for_seq_assembly(
        &graph_ref,
        &graph_reads,
        25,
        &assembler,
        false,
        false,
    )
    .expect("pruned25");

    fn seq_kbest_alt_trace(
        graph: &AssemblyGraph,
        ref_bases: &[u8],
        alt_windows: &[Vec<u8>],
        off: usize,
        ks: &[usize],
    ) -> Value {
        let mut seq = SeqGraph::from_assembly_graph(graph);
        seq.clean_non_ref_paths();
        let status = seq.cleanup_seq_graph();
        let ref_path = seq.reference_path_bytes();
        let mut by_k = Vec::new();
        for &k in ks {
            let paths = find_best_haplotypes_seq_graph(&seq, k).unwrap_or_default();
            let mut substr = 0usize;
            let mut idx_alt = 0usize;
            let mut eq_len = 0usize;
            let mut first_index_eq_alt: Option<usize> = None;
            let mut first_index_substr: Option<usize> = None;
            for (i, p) in paths.iter().enumerate() {
                let b = seq.path_bases_bytes(p.start, &p.edges);
                if b.len() == ref_bases.len() {
                    eq_len += 1;
                    if b.get(off) == Some(&TARGET_ALT) {
                        idx_alt += 1;
                        if first_index_eq_alt.is_none() {
                            first_index_eq_alt = Some(i);
                        }
                    }
                }
                if seq_contains_any_kmer(&b, alt_windows) {
                    substr += 1;
                    if first_index_substr.is_none() {
                        first_index_substr = Some(i);
                    }
                }
            }
            by_k.push(json!({
                "k_best": k,
                "n_paths": paths.len(),
                "eq_ref_len": eq_len,
                "index_alt": idx_alt,
                "first_index_eq_alt": first_index_eq_alt,
                "first_index_substr": first_index_substr,
                "substr_alt_kmer": substr,
            }));
        }
        json!({
            "seq_nodes": seq.node_count(),
            "seq_edges": seq.edge_count(),
            "cleanup": format!("{status:?}"),
            "ref_path_len": ref_path.as_ref().map(|b| b.len()),
            "ref_path_eq": ref_path.as_ref().map(|b| b.as_slice() == ref_bases),
            "kbest": by_k,
        })
    }

    let stage_e = json!({
        "offset_in_graph_ref": off,
        "graph_ref_len": graph_ref.bases.len(),
        "alt10_windows": alt10.len(),
        "alt25_windows": alt25.len(),
        "reads_with_alt10": alt10.iter().map(|km| count_reads_containing_kmer(&graph_reads, km)).max().unwrap_or(0),
        "reads_with_alt25": alt25.iter().map(|km| count_reads_containing_kmer(&graph_reads, km)).max().unwrap_or(0),
        "raw_k10": graph_alt_kmer_stats(&raw10, &alt10),
        "raw_k25": graph_alt_kmer_stats(&raw25, &alt25),
        "pruned_dangling_k10": pruned10.as_ref().map(|g| graph_alt_kmer_stats(g, &alt10)),
        "pruned_dangling_k25": pruned25.as_ref().map(|g| graph_alt_kmer_stats(g, &alt25)),
        "seq_kbest_k10": pruned10.as_ref().map(|g| seq_kbest_alt_trace(g, graph_ref.bases.as_slice(), &alt10, off, &[128, 256, 512])),
        "seq_kbest_k25": pruned25.as_ref().map(|g| seq_kbest_alt_trace(g, graph_ref.bases.as_slice(), &alt25, off, &[DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, 256, 512])),
        "rt_before_k10": hap_compact(&extract_rt_haplotypes_before_remove_paths(&graph_ref, &graph_reads, &assembler, 10, false, false).unwrap_or_default(), pad, TARGET_POS, &alt10, graph_ref.bases.len()),
        "rt_after_k10": hap_compact(&extract_rt_haplotypes_after_remove_paths(&graph_ref, &graph_reads, &assembler, 10, false, false).unwrap_or_default(), pad, TARGET_POS, &alt10, graph_ref.bases.len()),
        "rt_before_k25": hap_compact(&extract_rt_haplotypes_before_remove_paths(&graph_ref, &graph_reads, &assembler, 25, false, false).unwrap_or_default(), pad, TARGET_POS, &alt25, graph_ref.bases.len()),
        "rt_after_k25": hap_compact(&extract_rt_haplotypes_after_remove_paths(&graph_ref, &graph_reads, &assembler, 25, false, false).unwrap_or_default(), pad, TARGET_POS, &alt25, graph_ref.bases.len()),
        "assemble_from_ref_and_reads": hap_compact(&seq_raw.haplotypes, pad, TARGET_POS, &alt25, graph_ref.bases.len()),
        "assemble_untrimmed": hap_compact(&untrimmed.haplotypes, pad, TARGET_POS, &alt25, graph_ref.bases.len()),
    });

    let call =
        HaplotypeCallerEngine::call_region(region, &dict, &ref_fasta, &args).expect("call_region");
    let (trim_haps, trim_events, trim_has_alt) = match &call {
        Some(outcome) => {
            let tpad = outcome.assembly.padded_reference_start_1based();
            let th = outcome
                .assembly
                .haplotypes
                .iter()
                .map(|h| hap_row(h, tpad, TARGET_POS))
                .collect::<Vec<_>>();
            let ev: Vec<Value> = outcome
                .assembly
                .variation_events()
                .iter()
                .map(|e| {
                    json!({
                        "pos": e.start_1based.get(),
                        "ref": e.ref_allele,
                        "alt": e.alt_allele,
                    })
                })
                .collect();
            let has = hap_has_alt(&outcome.assembly.haplotypes, tpad, TARGET_POS);
            (th, ev, has)
        }
        None => (Vec::new(), Vec::new(), false),
    };

    let (gt_at_target, gt_has_ga, vcf_rows, vcf_has_ga) = match &call {
        Some(outcome) => {
            let gt: Vec<Value> = outcome
                .genotyped_calls
                .iter()
                .filter(|c| c.event.start_1based.get() == TARGET_POS)
                .map(|c| {
                    json!({
                        "pos": c.event.start_1based.get(),
                        "ref": c.event.ref_allele,
                        "alt": c.event.alt_allele,
                    })
                })
                .collect();
            let gt_has = outcome.genotyped_calls.iter().any(|c| {
                c.event.start_1based.get() == TARGET_POS
                    && c.event.ref_allele == "G"
                    && c.event.alt_allele == "A"
            });
            let emitted = try_emit_call_region_variants(
                region,
                outcome,
                "SAMPLE",
                DEFAULT_STAND_EMIT_CONFIDENCE,
            )
            .unwrap_or_default();
            let vcf: Vec<Value> = emitted
                .iter()
                .map(|r| {
                    json!({
                        "pos": r.position,
                        "ref": r.reference,
                        "alt": r.alternate,
                        "qual": r.quality,
                    })
                })
                .collect();
            let vcf_has = emitted.iter().any(|r| {
                r.position == TARGET_POS
                    && r.reference == "G"
                    && r.alternate.iter().any(|a| a == "A")
            });
            (gt, gt_has, vcf, vcf_has)
        }
        None => (Vec::new(), false, Vec::new(), false),
    };

    let seq_path_alt_count = seq_raw
        .haplotypes
        .iter()
        .filter(|h| {
            hap_base_from_seq(&h.bases, graph_ref.bases.as_slice(), pad, TARGET_POS)
                == Some(TARGET_ALT)
        })
        .count();

    let k128_row = stage_e["seq_kbest_k25"]["kbest"].as_array().and_then(|a| {
        a.iter()
            .find(|r| r["k_best"].as_u64() == Some(DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH as u64))
            .cloned()
    });
    let k128_index_alt = k128_row
        .as_ref()
        .and_then(|r| r["index_alt"].as_u64())
        .unwrap_or(0);
    let k128_substr = k128_row
        .as_ref()
        .and_then(|r| r["substr_alt_kmer"].as_u64())
        .unwrap_or(0);
    let kbest_k128_has_alt = k128_index_alt > 0 || k128_substr > 0;
    let eventmap_has = !event_hits(untrimmed.variation_events(), TARGET_POS).is_empty()
        || !event_hits_from_values(&trim_events, TARGET_POS).is_empty();
    let assemble_has = hap_has_alt(&untrimmed.haplotypes, pad, TARGET_POS);

    let first_loss = first_loss_label(
        pileup_final.get(&'A').copied().unwrap_or(0) > 0,
        kbest_k128_has_alt,
        rt_rows
            .iter()
            .any(|r| r["before"]["mapped_alt"].as_u64().unwrap_or(0) > 0),
        rt_rows
            .iter()
            .any(|r| r["after"]["mapped_alt"].as_u64().unwrap_or(0) > 0),
        seq_has_alt,
        assemble_has,
        trim_has_alt,
        eventmap_has,
        gt_has_ga,
        vcf_has_ga,
    );

    let doc = json!({
        "id": "chr20_tiny",
        "target": {"contig": "20", "pos": TARGET_POS, "ref": "G", "alt": "A"},
        "covering_activefull": {
            "active": [region.start.get(), region.end.get()],
            "extended": [region.extended_start.get(), region.extended_end.get()],
            "n_walker_reads": region.reads.len(),
            "n_finalized": assembled.finalized_reads.len(),
            "graph_ref_len": graph_ref.bases.len(),
            "graph_ref_start": pad,
            "kmer_used": k_used,
            "num_best_haplotypes_per_graph": assembler.num_best_haplotypes_per_graph,
            "kmer_sizes": assembler.kmer_sizes,
            "use_seq_graph": assembler.use_seq_graph,
            "min_prune_factor": assembler.min_prune_factor,
        },
        "pileup_raw_target": pileup_raw,
        "pileup_finalized_target": pileup_final,
        "pileup_raw_nearby": pileup_raw_near,
        "assemble_untrimmed": hap_compact(
            &untrimmed.haplotypes,
            pad,
            TARGET_POS,
            &alt25,
            graph_ref.bases.len(),
        ),
        "assembly_events_at_target": event_hits(untrimmed.variation_events(), TARGET_POS),
        "n_haps_with_target_eventmap": per_hap_events
            .iter()
            .filter(|e| e["target_events"].as_array().is_some_and(|a| !a.is_empty()))
            .count(),
        "seqgraph_assemble_from_ref_and_reads": {
            "status": format!("{:?}", seq_raw.status),
            "kmer": seq_raw.kmer_size,
            "has_alt_at_target": seq_has_alt,
            "seq_path_alt_count": seq_path_alt_count,
            "compact": hap_compact(
                &seq_raw.haplotypes,
                pad,
                TARGET_POS,
                &alt25,
                graph_ref.bases.len(),
            ),
        },
        "seqgraph_kmer_probe": probe_json,
        "rt_extract": rt_rows,
        "stage_e_chain": stage_e,
        "call_region_trimmed": {
            "n_haps": trim_haps.len(),
            "has_alt_at_target": trim_has_alt,
            "haps": trim_haps,
            "events_at_target": event_hits_from_values(&trim_events, TARGET_POS),
            "n_events": trim_events.len(),
        },
        "genotyping": {
            "n_calls": call.as_ref().map(|o| o.genotyped_calls.len()).unwrap_or(0),
            "at_target": gt_at_target,
            "has_ga": gt_has_ga,
        },
        "vcf_emit": {
            "n_records": vcf_rows.len(),
            "records": vcf_rows,
            "has_ga_at_target": vcf_has_ga,
        },
        "pipeline_survival": {
            "kbest_k128": k128_row,
            "assemble_from_ref_and_reads": seq_has_alt,
            "untrimmed_hap": assemble_has,
            "trim": trim_has_alt,
            "eventmap": eventmap_has,
            "genotyping": gt_has_ga,
            "vcf": vcf_has_ga,
        },
        "first_loss": first_loss,
    });
    println!("{}", serde_json::to_string_pretty(&doc).unwrap());
    assert!(
        covering.len() == 1,
        "expected a single covering ActiveFull, got {}",
        covering.len()
    );
    assert!(
        kbest_k128_has_alt,
        "6R.52: production K={} must contain the A-bearing path; k128={:?}",
        DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, k128_row
    );
}

fn event_hits_from_values(events: &[Value], loc: u64) -> Vec<Value> {
    events
        .iter()
        .filter(|e| e["pos"].as_u64() == Some(loc))
        .cloned()
        .collect()
}

fn first_loss_label(
    reads_have_alt: bool,
    kbest_k128: bool,
    rt_before: bool,
    rt_after: bool,
    seq_extract: bool,
    assemble: bool,
    trimmed: bool,
    eventmap: bool,
    genotyped: bool,
    vcf: bool,
) -> &'static str {
    if !reads_have_alt {
        "absent_from_finalized_reads"
    } else if rt_before && !rt_after {
        "lost_at_remove_paths_not_connected_to_ref"
    } else if !kbest_k128 {
        "lost_at_seqgraph_kbest_k128"
    } else if !seq_extract {
        "lost_at_haplotype_construction_or_sw_extract"
    } else if seq_extract && !assemble {
        "lost_at_findBestPaths_retention_or_normalize"
    } else if assemble && !trimmed {
        "lost_at_trim"
    } else if trimmed && !eventmap {
        "lost_at_eventmap"
    } else if eventmap && !genotyped {
        "lost_at_genotyping"
    } else if genotyped && !vcf {
        "lost_at_vcf_emit"
    } else if vcf {
        "present_in_vcf"
    } else {
        "present_through_eventmap_downstream_unknown"
    }
}
