//! 6R.57: earliest Java/Rust loss of `20:29455019 G/A`.
//!
//! After 6R.56 SeqGraph control flow, the remaining Java-only SNP is present on
//! retained haplotypes and per-haplotype EventMaps, then dropped by the Rust-only
//! union filter `prefer_indel_over_colocated_snps` when another haplotype has an
//! insertion at the same start. Production change: do not apply that filter on
//! EventMap *union* (Java `getAllVariantContexts`).
//!
//! Skipped unless `HOLDOUT_6R57=1`.
//!
//! ```text
//! HOLDOUT_6R57=1 cargo test -p gatk-haplotypecaller --test holdout_6r57_chr20_tiny_test -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::assembly_region_finalize::{
    assembly_reference_read, create_graph_reference_read, records_to_assembly_reads,
};
use gatk_haplotypecaller::event_map::variation_events_for_haplotype;
use gatk_haplotypecaller::haplotype::prune_fragment_non_reference_haplotypes;
use gatk_haplotypecaller::haplotype_cigar::trace_find_best_paths_gates;
use gatk_haplotypecaller::hc_allele_mapping::hap_base_at_ref_locus;
use gatk_haplotypecaller::hc_genotyping_engine::DEFAULT_STAND_EMIT_CONFIDENCE;
use gatk_haplotypecaller::read_threading_assembler::{
    assemble_from_ref_and_reads, build_threading_graph_for_seq_assembly,
    extract_haplotypes_from_seq_kbest_paths, extract_rt_haplotypes_after_remove_paths,
    extract_rt_haplotypes_before_remove_paths, AssemblyScoringContext,
    DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, MIN_HAPLOTYPE_REFERENCE_LENGTH,
};
use gatk_haplotypecaller::seq_graph::SeqGraph;
use gatk_haplotypecaller::seq_kbest_haplotype::find_best_haplotypes_seq_graph;
use gatk_haplotypecaller::{
    assemble_reads_with_finalized, assembly_graph_from_ref_and_reads_threading,
    calculate_haplotype_cigar_for_assembly, call_disposition, flatten_assembly_regions,
    probe_seq_graph_kmer_attempts, query_index_at_reference_position,
    reference_has_non_unique_kmers, traverse_assembly_region_walker, try_emit_call_region_variants,
    AssemblyGraph, AssemblyGraphParams, AssemblyRegionCallDisposition, CallRegionArgs, Cigar,
    CigarOperator, Haplotype, HaplotypeCallerEngine, KmerSize, ReadFilterParams, SwParameters,
    WalkerTraversalConfig, GATK_HC_DEFAULT_MIN_MAPPING_QUALITY,
};
use rust_htslib::bam::record::CigarString;
use rust_htslib::bam::Read;
use serde_json::{json, Value};
use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_VCF_REL: &str = "parity/reports/6r43/chr20_tiny/java.vcf";
const TARGET: u64 = 29_455_019;
const TARGET_REF: u8 = b'G';
const TARGET_ALT: u8 = b'A';
const K: usize = DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn cigar_str(h: &Haplotype) -> String {
    h.cigar
        .as_ref()
        .map(|c| c.to_gatk_string())
        .unwrap_or_else(|| "-".to_string())
}

fn hap_row(h: &Haplotype, pad: u64) -> Value {
    json!({
        "is_ref": h.is_reference,
        "len": h.bases.len(),
        "k": h.kmer_size,
        "score": h.score,
        "cigar": cigar_str(h),
        "cigar_ref_len": h.cigar.as_ref().map(|c| c.reference_length()),
        "align0": h.alignment_start_hap_wrt_ref,
        "base_at_target": hap_base_at_ref_locus(h, pad, TARGET)
            .map(|b| (b as char).to_string())
            .unwrap_or_else(|| ".".to_string()),
        "eq_alt_a": hap_base_at_ref_locus(h, pad, TARGET) == Some(TARGET_ALT),
        "eq_ref_g": hap_base_at_ref_locus(h, pad, TARGET) == Some(TARGET_REF),
    })
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

fn mq_hist(reads: &[rust_htslib::bam::Record]) -> BTreeMap<String, usize> {
    let mut h: BTreeMap<String, usize> = BTreeMap::new();
    for rec in reads {
        let mq = rec.mapq();
        let bucket = if mq == 255 {
            "255_unavailable".to_string()
        } else if mq < GATK_HC_DEFAULT_MIN_MAPPING_QUALITY {
            format!("lt_{GATK_HC_DEFAULT_MIN_MAPPING_QUALITY}")
        } else {
            "pass".to_string()
        };
        *h.entry(bucket).or_insert(0) += 1;
    }
    h
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

fn hap_compact(haps: &[Haplotype], pad: u64, alt_windows: &[Vec<u8>], ref_len: usize) -> Value {
    let mut cigars: BTreeMap<String, usize> = BTreeMap::new();
    let mut mapped_alt = 0usize;
    let mut substr_alt = 0usize;
    let mut eq_len_idx_alt = 0usize;
    let mut a_rows = Vec::new();
    for h in haps {
        *cigars.entry(cigar_str(h)).or_insert(0) += 1;
        let mapped = hap_base_at_ref_locus(h, pad, TARGET) == Some(TARGET_ALT);
        if mapped {
            mapped_alt += 1;
            if a_rows.len() < 12 {
                a_rows.push(hap_row(h, pad));
            }
        }
        if seq_contains_any_kmer(&h.bases, alt_windows) {
            substr_alt += 1;
        }
        if h.bases.len() == ref_len {
            let off = TARGET.saturating_sub(pad) as usize;
            if h.bases.get(off) == Some(&TARGET_ALT) {
                eq_len_idx_alt += 1;
            }
        }
    }
    json!({
        "n": haps.len(),
        "mapped_alt": mapped_alt,
        "substr_alt_kmer": substr_alt,
        "eq_len_index_alt": eq_len_idx_alt,
        "cigars": cigars,
        "a_bearing_sample": a_rows,
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

fn java_vcf_has_ga(path: &Path) -> Option<String> {
    if !path.is_file() {
        return None;
    }
    fs::read_to_string(path)
        .ok()?
        .lines()
        .find(|l| {
            let f: Vec<_> = l.split('\t').collect();
            f.len() >= 5 && f[0] == "20" && f[1] == "29455019" && f[3] == "G" && f[4] == "A"
        })
        .map(str::to_string)
}

fn java_bamout_observe(root: &Path) -> Value {
    let candidates = [
        "parity/reports/6r51/chr20_tiny/java_full.bamout.bam",
        "parity/reports/6r51/chr20_tiny/java_covering.bamout.bam",
    ];
    for rel in candidates {
        let p = root.join(rel);
        if !p.is_file() {
            continue;
        }
        let Ok(mut bam) = rust_htslib::bam::Reader::from_path(&p) else {
            continue;
        };
        let mut hc: BTreeMap<String, (String, usize, char, i64)> = BTreeMap::new();
        for rec in bam.records().flatten() {
            if rec.is_unmapped() {
                continue;
            }
            let q = String::from_utf8_lossy(rec.qname()).into_owned();
            if !q.starts_with("HC_") {
                continue;
            }
            let cigar = CigarString(rec.cigar().iter().copied().collect());
            let pos0 = (TARGET - 1) as i64;
            let base = query_index_at_reference_position(rec.pos(), &cigar, pos0)
                .map(|qi| rec.seq()[qi])
                .map(|b| (b as char).to_ascii_uppercase())
                .unwrap_or('.');
            let cstr = cigar.to_string();
            hc.entry(q)
                .or_insert((cstr, rec.seq_len(), base, rec.pos()));
        }
        let n_a = hc.values().filter(|v| v.2 == 'A').count();
        let n_g = hc.values().filter(|v| v.2 == 'G').count();
        let cigars: BTreeMap<String, usize> = {
            let mut m = BTreeMap::new();
            for v in hc.values() {
                *m.entry(v.0.clone()).or_insert(0) += 1;
            }
            m
        };
        return json!({
            "source": rel,
            "n_hc_qnames": hc.len(),
            "n_a": n_a,
            "n_g": n_g,
            "cigars": cigars,
            "align_starts_0based": hc.values().map(|v| v.3).collect::<std::collections::BTreeSet<_>>(),
        });
    }
    json!({"status": "UNKNOWN", "reason": "no java bamout"})
}

fn seq_kbest_alt_trace(
    graph: &AssemblyGraph,
    ref_bases: &[u8],
    alt_windows: &[Vec<u8>],
    pad: u64,
    off: usize,
    ks: &[usize],
    sw: &SwParameters,
) -> Value {
    let mut seq_before = SeqGraph::from_assembly_graph(graph);
    let nodes_before = seq_before.node_count();
    let edges_before = seq_before.edge_count();
    seq_before.clean_non_ref_paths();
    let status = seq_before.cleanup_seq_graph();
    let seq = seq_before;
    let mut ref_hap = Haplotype::new(ref_bases, true);
    let mut rc = Cigar::new();
    rc.push(ref_hap.bases.len(), CigarOperator::Match);
    ref_hap.cigar = Some(rc);
    let ref_cigar_len = ref_bases.len();

    let mut by_k = Vec::new();
    for &k in ks {
        let paths = find_best_haplotypes_seq_graph(&seq, k).unwrap_or_default();
        let mut substr = 0usize;
        let mut idx_alt = 0usize;
        let mut eq_len = 0usize;
        let mut mapped_alt = 0usize;
        let mut first_index_eq_alt: Option<usize> = None;
        let mut first_index_substr: Option<usize> = None;
        let mut first_index_mapped: Option<usize> = None;
        let mut a_bearing = Vec::new();
        let mut seen: Vec<(Vec<u8>, bool)> = Vec::new();
        for (i, p) in paths.iter().enumerate() {
            let b = seq.path_bases_bytes(p.start, &p.edges);
            let eq = b.len() == ref_bases.len();
            if eq {
                eq_len += 1;
            }
            let index_hit = eq && b.get(off) == Some(&TARGET_ALT);
            if index_hit {
                idx_alt += 1;
                if first_index_eq_alt.is_none() {
                    first_index_eq_alt = Some(i);
                }
            }
            let substr_hit = seq_contains_any_kmer(&b, alt_windows);
            if substr_hit {
                substr += 1;
                if first_index_substr.is_none() {
                    first_index_substr = Some(i);
                }
            }
            let mapped = hap_base_from_seq(&b, ref_bases, pad, TARGET) == Some(TARGET_ALT);
            if mapped {
                mapped_alt += 1;
                if first_index_mapped.is_none() {
                    first_index_mapped = Some(i);
                }
            }
            if index_hit || mapped || (substr_hit && a_bearing.len() < 8) {
                let t = trace_find_best_paths_gates(
                    ref_bases,
                    &b,
                    p.is_reference,
                    ref_cigar_len,
                    sw,
                    &seen,
                );
                if t.rust_extract_keep {
                    seen.push((b.clone(), p.is_reference));
                }
                if index_hit || mapped {
                    a_bearing.push(json!({
                        "kbest_ordinal": i,
                        "score": p.score,
                        "finite": p.score.is_finite(),
                        "is_reference_flag": p.is_reference,
                        "seq_len": t.seq_len,
                        "eq_ref_len": eq,
                        "base_at_offset": b.get(off).copied().map(|c| (c as char).to_string()),
                        "rust_prod_cigar": t.rust_prod_cigar,
                        "rust_prod_ref_len": t.rust_prod_ref_len,
                        "java_would_retain": t.java_would_retain,
                        "rust_extract_keep": t.rust_extract_keep,
                        "first_rust_reject": t.first_rust_reject,
                        "first_java_reject": t.first_java_reject,
                        "duplicate": t.duplicate,
                    }));
                }
            }
        }
        let extracted = extract_haplotypes_from_seq_kbest_paths(
            &paths,
            &seq,
            graph.kmer_size,
            &ref_hap,
            ref_cigar_len,
            sw,
        )
        .unwrap_or_default();
        let extract_has = extracted
            .iter()
            .any(|h| hap_base_at_ref_locus(h, pad, TARGET) == Some(TARGET_ALT));
        let mut pruned = extracted.clone();
        prune_fragment_non_reference_haplotypes(
            &mut pruned,
            &ref_hap,
            MIN_HAPLOTYPE_REFERENCE_LENGTH,
        );
        let prune_has = pruned
            .iter()
            .any(|h| hap_base_at_ref_locus(h, pad, TARGET) == Some(TARGET_ALT));
        by_k.push(json!({
            "k_best": k,
            "n_paths": paths.len(),
            "eq_ref_len": eq_len,
            "index_alt": idx_alt,
            "mapped_alt": mapped_alt,
            "first_index_eq_alt": first_index_eq_alt,
            "first_index_mapped": first_index_mapped,
            "first_index_substr": first_index_substr,
            "substr_alt_kmer": substr,
            "a_bearing": a_bearing,
            "extract_n": extracted.len(),
            "extract_has_a": extract_has,
            "after_prune_fragment_has_a": prune_has,
        }));
    }
    json!({
        "seq_nodes_before_cleanup": nodes_before,
        "seq_edges_before_cleanup": edges_before,
        "seq_nodes": seq.node_count(),
        "seq_edges": seq.edge_count(),
        "cleanup": format!("{status:?}"),
        "kbest": by_k,
    })
}

#[test]
fn holdout_6r57_chr20_tiny_29455019_ga() {
    if std::env::var("HOLDOUT_6R57").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R57=1");
        return;
    }
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    let java_vcf = root.join(JAVA_VCF_REL);
    assert!(ref_fasta.is_file(), "missing {}", ref_fasta.display());
    assert!(bam.is_file(), "missing {}", bam.display());

    let java_line = java_vcf_has_ga(&java_vcf);
    let java_bam = java_bamout_observe(&root);

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
    assert!(!covering.is_empty(), "no ActiveFull covers {TARGET}");
    let region = covering[0];
    let args = CallRegionArgs::strict_java();
    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);

    let raw_reads: Vec<_> = region.reads.iter().map(|r| r.as_ref().clone()).collect();
    let pileup_raw = pileup_at(&raw_reads, TARGET);
    let mq_raw = mq_hist(&raw_reads);

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
        .cloned();
    let ref_bytes = untrimmed.apply_bases_shared();

    let mut per_hap_events = Vec::new();
    if let Some(ref_h) = ref_hap.as_ref() {
        for (i, h) in untrimmed.haplotypes.iter().enumerate() {
            let ev = variation_events_for_haplotype(
                h,
                ref_h,
                ref_bytes.as_ref(),
                pad,
                0,
                &region.contig,
            );
            let hits = event_hits(&ev, TARGET);
            if !hits.is_empty() {
                per_hap_events.push(json!({
                    "i": i,
                    "is_ref": h.is_reference,
                    "cigar": cigar_str(h),
                    "target_events": hits,
                }));
            }
        }
    }

    let padded_ref = assembly_reference_read(&dict, &mut ref_cache, region).expect("pad ref");
    let graph_ref = create_graph_reference_read(&padded_ref, region, &dict);
    let graph_reads = records_to_assembly_reads(&assembled.finalized_reads);
    let pileup_final = pileup_at(&assembled.finalized_reads, TARGET);
    let mq_final = mq_hist(&assembled.finalized_reads);

    let mut assembler = args.assemble.assembler.clone();
    assembler.dangling_java_exact = true;
    assembler.scoring = Some(AssemblyScoringContext {
        padded_reference_start_1based: region.extended_start.get(),
        active_start_1based: region.start.get(),
        active_end_1based: region.end.get(),
        contig: region.contig.clone(),
    });
    let off = TARGET.saturating_sub(pad) as usize;
    assert!(off < graph_ref.bases.len(), "target outside graph ref");
    let mut alt_ref = graph_ref.bases.clone();
    alt_ref[off] = TARGET_ALT;
    let alt10 = alt_only_windows(&graph_ref.bases, &alt_ref, off, 10);
    let alt25 = alt_only_windows(&graph_ref.bases, &alt_ref, off, 25);
    let sw = assembler.haplotype_to_reference_sw;

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
            })
        })
        .collect();

    let seq_raw = assemble_from_ref_and_reads(&graph_ref, &graph_reads, &assembler).expect("seq");

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

    let rt_before_25 = extract_rt_haplotypes_before_remove_paths(
        &graph_ref,
        &graph_reads,
        &assembler,
        25,
        false,
        false,
    )
    .unwrap_or_default();
    let rt_after_25 = extract_rt_haplotypes_after_remove_paths(
        &graph_ref,
        &graph_reads,
        &assembler,
        25,
        false,
        false,
    )
    .unwrap_or_default();
    let rt_before_10 = extract_rt_haplotypes_before_remove_paths(
        &graph_ref,
        &graph_reads,
        &assembler,
        10,
        false,
        false,
    )
    .unwrap_or_default();
    let rt_after_10 = extract_rt_haplotypes_after_remove_paths(
        &graph_ref,
        &graph_reads,
        &assembler,
        10,
        false,
        false,
    )
    .unwrap_or_default();

    let seq_kbest_k25 = pruned25.as_ref().map(|g| {
        seq_kbest_alt_trace(
            g,
            graph_ref.bases.as_slice(),
            &alt25,
            pad,
            off,
            &[K, 256, 512],
            &sw,
        )
    });
    let seq_kbest_k10 = pruned10.as_ref().map(|g| {
        seq_kbest_alt_trace(
            g,
            graph_ref.bases.as_slice(),
            &alt10,
            pad,
            off,
            &[K, 256],
            &sw,
        )
    });

    let call =
        HaplotypeCallerEngine::call_region(region, &dict, &ref_fasta, &args).expect("call_region");
    let (trim_has, trim_events, trim_compact, gt_has, vcf_has, vcf_rows) = match &call {
        Some(outcome) => {
            let tpad = outcome.assembly.padded_reference_start_1based();
            let has = outcome
                .assembly
                .haplotypes
                .iter()
                .any(|h| hap_base_at_ref_locus(h, tpad, TARGET) == Some(TARGET_ALT));
            let tev = event_hits(outcome.assembly.variation_events(), TARGET);
            let compact = hap_compact(
                &outcome.assembly.haplotypes,
                tpad,
                &alt25,
                graph_ref.bases.len(),
            );
            let gt_has = outcome.genotyped_calls.iter().any(|c| {
                c.event.start_1based.get() == TARGET
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
                .filter(|r| r.position == TARGET)
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
                r.position == TARGET && r.reference == "G" && r.alternate.iter().any(|a| a == "A")
            });
            (has, tev, compact, gt_has, vcf_has, vcf)
        }
        None => (false, Vec::new(), json!(null), false, false, Vec::new()),
    };

    let k128_row = seq_kbest_k25.as_ref().and_then(|v| {
        v["kbest"].as_array().and_then(|a| {
            a.iter()
                .find(|r| r["k_best"].as_u64() == Some(K as u64))
                .cloned()
        })
    });
    let k256_row = seq_kbest_k25.as_ref().and_then(|v| {
        v["kbest"].as_array().and_then(|a| {
            a.iter()
                .find(|r| r["k_best"].as_u64() == Some(256))
                .cloned()
        })
    });
    let k128_has = k128_row
        .as_ref()
        .map(|r| {
            r["index_alt"].as_u64().unwrap_or(0) > 0
                || r["mapped_alt"].as_u64().unwrap_or(0) > 0
                || r["substr_alt_kmer"].as_u64().unwrap_or(0) > 0
        })
        .unwrap_or(false);
    let k256_has = k256_row
        .as_ref()
        .map(|r| {
            r["index_alt"].as_u64().unwrap_or(0) > 0
                || r["mapped_alt"].as_u64().unwrap_or(0) > 0
                || r["substr_alt_kmer"].as_u64().unwrap_or(0) > 0
        })
        .unwrap_or(false);
    let extract_has = k128_row
        .as_ref()
        .and_then(|r| r["extract_has_a"].as_bool())
        .unwrap_or(false);
    let assemble_has = hap_base_at_ref_locus_any(&untrimmed.haplotypes, pad);
    let seq_extract_has = hap_base_at_ref_locus_any(&seq_raw.haplotypes, pad);
    let eventmap_union = event_hits(untrimmed.variation_events(), TARGET);
    let eventmap_has_ga = eventmap_union
        .iter()
        .any(|e| e["ref"] == "G" && e["alt"] == "A");
    let per_hap_ga = per_hap_events.iter().any(|row| {
        row["target_events"]
            .as_array()
            .is_some_and(|a| a.iter().any(|e| e["ref"] == "G" && e["alt"] == "A"))
    });
    let pileup_has_a = pileup_final.get(&'A').copied().unwrap_or(0) > 0;
    let rt_before_has = hap_base_at_ref_locus_any(&rt_before_25, pad);
    let rt_after_has = hap_base_at_ref_locus_any(&rt_after_25, pad);
    let graph_has = pruned25.as_ref().is_some_and(|g| {
        graph_alt_kmer_stats(g, &alt25)["alt_nodes"]
            .as_u64()
            .unwrap_or(0)
            > 0
    });
    let raw_has = graph_alt_kmer_stats(&raw25, &alt25)["alt_nodes"]
        .as_u64()
        .unwrap_or(0)
        > 0;

    let first_loss = if !pileup_has_a {
        "A_finalized_pileup"
    } else if !raw_has && !rt_before_has {
        "B_rt_graph"
    } else if raw_has && !graph_has {
        "B_rt_prune"
    } else if graph_has && !k256_has && !k128_has {
        "C_or_D_seqgraph_or_kbest_absent_even_at_k256"
    } else if !k128_has && k256_has {
        "D_kbest_k128_absent_present_in_larger_search"
    } else if k128_has && !extract_has {
        "E_findBestPaths_retention"
    } else if extract_has && !seq_extract_has && !assemble_has {
        "E_assemble_post_extract"
    } else if assemble_has && !eventmap_has_ga && !per_hap_ga {
        "F_A_event_never_constructed"
    } else if per_hap_ga && !eventmap_has_ga {
        "F_B_lost_at_union"
    } else if eventmap_has_ga && !trim_has {
        "G_trim"
    } else if eventmap_has_ga && trim_has && !gt_has {
        "H_genotyping"
    } else if gt_has && !vcf_has {
        "I_vcf_emit"
    } else if vcf_has {
        "present_in_vcf"
    } else {
        "unclassified"
    };

    let doc = json!({
        "id": "chr20_tiny",
        "target": {"contig": "20", "pos": TARGET, "ref": "G", "alt": "A"},
        "k_production": K,
        "min_mapping_quality": GATK_HC_DEFAULT_MIN_MAPPING_QUALITY,
        "java_oracle": {
            "vcf_line_present": java_line.is_some(),
            "vcf_line": java_line,
            "bamout": java_bam,
            "graph_dot": "UNKNOWN",
        },
        "covering_activefull": {
            "n_covering": covering.len(),
            "active": [region.start.get(), region.end.get()],
            "extended": [region.extended_start.get(), region.extended_end.get()],
            "n_walker_reads": region.reads.len(),
            "n_finalized": assembled.finalized_reads.len(),
            "graph_ref_len": graph_ref.bases.len(),
            "graph_ref_start": pad,
            "target_offset": off,
            "kmer_used": k_used,
            "num_best_haplotypes_per_graph": assembler.num_best_haplotypes_per_graph,
            "kmer_sizes": assembler.kmer_sizes,
            "use_seq_graph": assembler.use_seq_graph,
            "min_prune_factor": assembler.min_prune_factor,
            "ref_unique_k10": !reference_has_non_unique_kmers(&graph_ref, 10),
            "ref_unique_k25": !reference_has_non_unique_kmers(&graph_ref, 25),
        },
        "stage_a_pileup": {
            "raw": pileup_raw,
            "finalized": pileup_final,
            "mq_raw": mq_raw,
            "mq_finalized": mq_final,
        },
        "stage_b_rt": {
            "reads_with_alt10": alt10.iter().map(|km| count_reads_containing_kmer(&graph_reads, km)).max().unwrap_or(0),
            "reads_with_alt25": alt25.iter().map(|km| count_reads_containing_kmer(&graph_reads, km)).max().unwrap_or(0),
            "alt10_windows": alt10.len(),
            "alt25_windows": alt25.len(),
            "raw_k10": graph_alt_kmer_stats(&raw10, &alt10),
            "raw_k25": graph_alt_kmer_stats(&raw25, &alt25),
            "pruned_dangling_k10": pruned10.as_ref().map(|g| graph_alt_kmer_stats(g, &alt10)),
            "pruned_dangling_k25": pruned25.as_ref().map(|g| graph_alt_kmer_stats(g, &alt25)),
            "rt_before_k10": hap_compact(&rt_before_10, pad, &alt10, graph_ref.bases.len()),
            "rt_after_k10": hap_compact(&rt_after_10, pad, &alt10, graph_ref.bases.len()),
            "rt_before_k25": hap_compact(&rt_before_25, pad, &alt25, graph_ref.bases.len()),
            "rt_after_k25": hap_compact(&rt_after_25, pad, &alt25, graph_ref.bases.len()),
        },
        "stage_c_d_seqgraph": {
            "probe": probe_json,
            "k10": seq_kbest_k10,
            "k25": seq_kbest_k25,
        },
        "stage_e_assemble": {
            "assemble_from_ref_and_reads": hap_compact(&seq_raw.haplotypes, pad, &alt25, graph_ref.bases.len()),
            "assemble_untrimmed": hap_compact(&untrimmed.haplotypes, pad, &alt25, graph_ref.bases.len()),
            "status": format!("{:?}", seq_raw.status),
            "kmer": seq_raw.kmer_size,
        },
        "stage_f_eventmap": {
            "union_at_target": eventmap_union,
            "per_hap_with_target_event": per_hap_events,
        },
        "stage_g_trim": {
            "has_alt": trim_has,
            "events_at_target": trim_events,
            "compact": trim_compact,
        },
        "stage_h_i": {
            "genotyped_ga": gt_has,
            "vcf_ga": vcf_has,
            "vcf_at_target": vcf_rows,
        },
        "first_loss": first_loss,
        "survival": {
            "pileup_a": pileup_has_a,
            "raw_graph_alt_nodes": raw_has,
            "pruned_graph_alt_nodes": graph_has,
            "rt_before_k25": rt_before_has,
            "rt_after_k25": rt_after_has,
            "kbest_k128": k128_has,
            "kbest_k256": k256_has,
            "extract_k128": extract_has,
            "assemble_from_ref": seq_extract_has,
            "untrimmed": assemble_has,
            "per_hap_eventmap_ga": per_hap_ga,
            "eventmap_union_ga": eventmap_has_ga,
            "trim": trim_has,
            "genotyping": gt_has,
            "vcf": vcf_has,
        },
    });
    println!("{}", serde_json::to_string_pretty(&doc).unwrap());
    assert_eq!(K, 128, "do not raise K");
    assert_eq!(covering.len(), 1, "expected one covering ActiveFull");
    assert!(
        java_line.is_some(),
        "Java 4.4 oracle VCF must still contain 20:29455019 G/A"
    );
    assert!(
        per_hap_ga,
        "G/A must exist on a per-haplotype EventMap (6R.57 F-A would be CIGAR construction)"
    );
    assert!(
        eventmap_has_ga,
        "6R.57: EventMap union must keep G/A when another haplotype has an insertion at the same start"
    );
    assert!(
        vcf_has,
        "6R.57: G/A must reach VCF after EventMap union keeps the SNP (allele returned; QUAL not Java-equivalent)"
    );
}

fn hap_base_at_ref_locus_any(haps: &[Haplotype], pad: u64) -> bool {
    haps.iter()
        .any(|h| hap_base_at_ref_locus(h, pad, TARGET) == Some(TARGET_ALT))
}
