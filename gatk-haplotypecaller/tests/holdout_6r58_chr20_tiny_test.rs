//! 6R.58 forensic: `20:29456344` Java `T/C` vs Rust `TG/T`.
//!
//! Representation-first, then stage A–I. No production algorithm change in this test.
//!
//! Skipped unless `HOLDOUT_6R58=1`.
//!
//! ```text
//! HOLDOUT_6R58=1 cargo test -p gatk-haplotypecaller --test holdout_6r58_chr20_tiny_test -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::assembly_region_finalize::{
    assembly_reference_read, create_graph_reference_read, records_to_assembly_reads,
};
use gatk_haplotypecaller::event_map::{
    merged_biallelic_sites_at_position, prefer_indel_over_colocated_snps,
    variation_events_for_haplotype, VariationEvent,
};
use gatk_haplotypecaller::genotyping::biallelic_genotype_index_from_pl;
use gatk_haplotypecaller::haplotype_cigar::trace_find_best_paths_gates;
use gatk_haplotypecaller::hc_allele_mapping::hap_base_at_ref_locus;
use gatk_haplotypecaller::hc_emit_policy::{
    explain_strict_java_emit_gates, passes_strict_java_emit_for_genotyped_call,
};
use gatk_haplotypecaller::hc_genotyping_engine::DEFAULT_STAND_EMIT_CONFIDENCE;
use gatk_haplotypecaller::read_threading_assembler::{
    assemble_from_ref_and_reads, build_threading_graph_for_seq_assembly,
    extract_haplotypes_from_seq_kbest_paths, extract_rt_haplotypes_after_remove_paths,
    extract_rt_haplotypes_before_remove_paths, AssemblyScoringContext,
    DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH,
};
use gatk_haplotypecaller::seq_graph::SeqGraph;
use gatk_haplotypecaller::seq_kbest_haplotype::find_best_haplotypes_seq_graph;
use gatk_haplotypecaller::{
    assemble_reads_with_finalized, assembly_graph_from_ref_and_reads_threading, call_disposition,
    flatten_assembly_regions, query_index_at_reference_position, reference_has_non_unique_kmers,
    traverse_assembly_region_walker, try_emit_call_region_variants, AssemblyGraph,
    AssemblyGraphParams, AssemblyRegionCallDisposition, CallRegionArgs, Cigar, CigarOperator,
    Haplotype, HaplotypeCallerEngine, KmerSize, ReadFilterParams, WalkerTraversalConfig,
    GATK_HC_DEFAULT_MIN_MAPPING_QUALITY,
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
const RUST_VCF_REL: &str = "parity/reports/6r43/chr20_tiny/rust.vcf";
/// Java SNP POS.
const POS_SNP: u64 = 29_456_344;
/// Next ref base (G); deleted in Rust `TG/T`.
const POS_G: u64 = 29_456_345;
const POS_AT_DEL: u64 = 29_456_343;
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
    let b344 = hap_base_at_ref_locus(h, pad, POS_SNP);
    let b345 = hap_base_at_ref_locus(h, pad, POS_G);
    json!({
        "is_ref": h.is_reference,
        "len": h.bases.len(),
        "k": h.kmer_size,
        "score": h.score,
        "cigar": cigar_str(h),
        "cigar_ref_len": h.cigar.as_ref().map(|c| c.reference_length()),
        "align0": h.alignment_start_hap_wrt_ref,
        "base_29456344": b344.map(|b| (b as char).to_string()).unwrap_or_else(|| "-".into()),
        "base_29456345": b345.map(|b| (b as char).to_string()).unwrap_or_else(|| "-".into()),
        "is_java_tc_snp": b344 == Some(b'C') && b345 == Some(b'G'),
        "is_ref_tg": b344 == Some(b'T') && b345 == Some(b'G'),
        "deletes_t_344": b344.is_none(),
        "deletes_g_345": b345.is_none(),
        "snippet": snippet_at(h, pad, POS_SNP),
    })
}

fn snippet_at(h: &Haplotype, pad: u64, loc: u64) -> String {
    if h.bases.is_empty() {
        return String::new();
    }
    // Equal-length vs padded ref: offset in bases.
    let off = loc.saturating_sub(pad) as usize;
    if off < h.bases.len() {
        let lo = off.saturating_sub(8);
        let hi = (off + 10).min(h.bases.len());
        return String::from_utf8_lossy(&h.bases[lo..hi]).into_owned();
    }
    String::from("uneq_len")
}

fn pileup_at(reads: &[rust_htslib::bam::Record], loc_1based: u64) -> BTreeMap<String, usize> {
    let pos0 = (loc_1based - 1) as i64;
    let mut counts: BTreeMap<String, usize> = BTreeMap::new();
    for rec in reads {
        if rec.is_unmapped() {
            continue;
        }
        let cigar = CigarString(rec.cigar().iter().copied().collect());
        match query_index_at_reference_position(rec.pos(), &cigar, pos0) {
            Some(qi) => {
                let seq = rec.seq();
                if qi < seq.len() {
                    let b = (seq[qi] as char).to_ascii_uppercase();
                    *counts.entry(b.to_string()).or_insert(0) += 1;
                }
            }
            None => {
                *counts.entry("gap_or_del".to_string()).or_insert(0) += 1;
            }
        }
    }
    counts
}

/// Strip a shared suffix (VCF left-trim of TG/CG → T/C). Not a production normalizer.
fn strip_common_suffix(ref_a: &str, alt_a: &str) -> (String, String) {
    let mut r: Vec<char> = ref_a.chars().collect();
    let mut a: Vec<char> = alt_a.chars().collect();
    while r.len() > 1 && a.len() > 1 && r.last() == a.last() {
        r.pop();
        a.pop();
    }
    (r.into_iter().collect(), a.into_iter().collect())
}

fn event_hits(events: &[VariationEvent], lo: u64, hi: u64) -> Vec<Value> {
    events
        .iter()
        .filter(|e| e.start_1based.get() >= lo && e.start_1based.get() <= hi)
        .map(|e| {
            json!({
                "pos": e.start_1based.get(),
                "ref": e.ref_allele,
                "alt": e.alt_allele,
            })
        })
        .collect()
}

fn hap_compact(haps: &[Haplotype], pad: u64) -> Value {
    let mut cigars: BTreeMap<String, usize> = BTreeMap::new();
    let mut n_tc = 0usize;
    let mut n_ref_tg = 0usize;
    let mut n_del_t = 0usize;
    let mut n_del_g = 0usize;
    let mut tc_sample = Vec::new();
    let mut del_sample = Vec::new();
    for h in haps {
        *cigars.entry(cigar_str(h)).or_insert(0) += 1;
        let b344 = hap_base_at_ref_locus(h, pad, POS_SNP);
        let b345 = hap_base_at_ref_locus(h, pad, POS_G);
        if b344 == Some(b'C') && b345 == Some(b'G') {
            n_tc += 1;
            if tc_sample.len() < 6 {
                tc_sample.push(hap_row(h, pad));
            }
        }
        if b344 == Some(b'T') && b345 == Some(b'G') {
            n_ref_tg += 1;
        }
        if b344.is_none() {
            n_del_t += 1;
            if del_sample.len() < 4 {
                del_sample.push(hap_row(h, pad));
            }
        }
        if b345.is_none() {
            n_del_g += 1;
            if del_sample.len() < 8 {
                del_sample.push(hap_row(h, pad));
            }
        }
    }
    json!({
        "n": haps.len(),
        "n_java_T_to_C_with_G": n_tc,
        "n_ref_T_G": n_ref_tg,
        "n_deletes_T_29456344": n_del_t,
        "n_deletes_G_29456345": n_del_g,
        "cigars": cigars,
        "tc_snp_sample": tc_sample,
        "deletion_sample": del_sample,
    })
}

fn vcf_keys_at(path: &Path, lo: u64, hi: u64) -> Vec<String> {
    if !path.is_file() {
        return Vec::new();
    }
    fs::read_to_string(path)
        .unwrap_or_default()
        .lines()
        .filter_map(|l| {
            let f: Vec<_> = l.split('\t').collect();
            if f.len() < 5 || f[0] != "20" {
                return None;
            }
            let p: u64 = f[1].parse().ok()?;
            if p < lo || p > hi {
                return None;
            }
            Some(format!("{} {}/{}", p, f[3], f[4]))
        })
        .collect()
}

fn java_bamout_observe(root: &Path) -> Value {
    let p = root.join("parity/reports/6r51/chr20_tiny/java_full.bamout.bam");
    if !p.is_file() {
        return json!({"status": "UNKNOWN"});
    }
    let Ok(mut bam) = rust_htslib::bam::Reader::from_path(&p) else {
        return json!({"status": "UNKNOWN", "reason": "open failed"});
    };
    let mut n_tc = 0usize;
    let mut n_tg = 0usize;
    let mut n_other = 0usize;
    let mut cigars: BTreeMap<String, usize> = BTreeMap::new();
    let mut starts: BTreeMap<i64, usize> = BTreeMap::new();
    for rec in bam.records().flatten() {
        if rec.is_unmapped() {
            continue;
        }
        let q = String::from_utf8_lossy(rec.qname());
        if !q.starts_with("HC_") {
            continue;
        }
        let cigar = CigarString(rec.cigar().iter().copied().collect());
        let cstr = cigar.to_string();
        *cigars.entry(cstr.clone()).or_insert(0) += 1;
        *starts.entry(rec.pos() + 1).or_insert(0) += 1;
        let pos0 = (POS_SNP - 1) as i64;
        let pos0g = (POS_G - 1) as i64;
        let b344 = query_index_at_reference_position(rec.pos(), &cigar, pos0)
            .map(|qi| rec.seq()[qi] as char)
            .map(|c| c.to_ascii_uppercase());
        let b345 = query_index_at_reference_position(rec.pos(), &cigar, pos0g)
            .map(|qi| rec.seq()[qi] as char)
            .map(|c| c.to_ascii_uppercase());
        match (b344, b345) {
            (Some('C'), Some('G')) => n_tc += 1,
            (Some('T'), Some('G')) => n_tg += 1,
            _ => n_other += 1,
        }
    }
    json!({
        "source": "parity/reports/6r51/chr20_tiny/java_full.bamout.bam",
        "n_tc_snp_CG": n_tc,
        "n_ref_TG": n_tg,
        "n_other_or_indel_or_noncovering": n_other,
        "cigars_covering_window": cigars,
        "align_starts_1based": starts,
    })
}

fn graph_alt_kmer_stats(graph: &AssemblyGraph, alt_windows: &[Vec<u8>]) -> Value {
    let mut nodes_hit = 0usize;
    let mut max_edge = 0u32;
    let alt_ids: std::collections::HashSet<usize> = graph
        .nodes()
        .iter()
        .enumerate()
        .filter(|(_, n)| alt_windows.iter().any(|km| &*n.kmer == km.as_slice()))
        .map(|(i, _)| {
            nodes_hit += 1;
            i
        })
        .collect();
    let mut edges_touching = 0usize;
    for e in graph.edges_sorted() {
        if alt_ids.contains(&e.from) || alt_ids.contains(&e.to) {
            edges_touching += 1;
            max_edge = max_edge.max(e.support);
        }
    }
    json!({
        "nodes": graph.node_count(),
        "edges": graph.edge_count(),
        "alt_nodes": nodes_hit,
        "alt_touching_edges": edges_touching,
        "max_alt_edge_support": max_edge,
    })
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

fn seq_kbest_trace(
    graph: &AssemblyGraph,
    ref_bases: &[u8],
    pad: u64,
    off: usize,
    sw: &gatk_haplotypecaller::SwParameters,
) -> Value {
    let mut seq = SeqGraph::from_assembly_graph(graph);
    let before = (seq.node_count(), seq.edge_count());
    seq.clean_non_ref_paths();
    let status = seq.cleanup_seq_graph();
    let mut ref_hap = Haplotype::new(ref_bases, true);
    let mut rc = Cigar::new();
    rc.push(ref_hap.bases.len(), CigarOperator::Match);
    ref_hap.cigar = Some(rc);
    let paths = find_best_haplotypes_seq_graph(&seq, K).unwrap_or_default();
    let mut n_tc = 0usize;
    let mut n_eq = 0usize;
    let mut first_tc = None;
    let mut seen: Vec<(Vec<u8>, bool)> = Vec::new();
    let mut tc_rows = Vec::new();
    for (i, p) in paths.iter().enumerate() {
        let b = seq.path_bases_bytes(p.start, &p.edges);
        let eq = b.len() == ref_bases.len();
        if eq {
            n_eq += 1;
        }
        let is_tc = eq && b.get(off) == Some(&b'C') && b.get(off + 1) == Some(&b'G');
        if is_tc {
            n_tc += 1;
            if first_tc.is_none() {
                first_tc = Some(i);
            }
            let t = trace_find_best_paths_gates(
                ref_bases,
                &b,
                p.is_reference,
                ref_bases.len(),
                sw,
                &seen,
            );
            if t.rust_extract_keep {
                seen.push((b.clone(), p.is_reference));
            }
            if tc_rows.len() < 8 {
                tc_rows.push(json!({
                    "kbest_ordinal": i,
                    "score": p.score,
                    "finite": p.score.is_finite(),
                    "java_would_retain": t.java_would_retain,
                    "rust_extract_keep": t.rust_extract_keep,
                    "first_rust_reject": t.first_rust_reject,
                    "rust_prod_cigar": t.rust_prod_cigar,
                    "snippet": String::from_utf8_lossy(&b[off.saturating_sub(8)..(off+10).min(b.len())]).into_owned(),
                }));
            }
        }
    }
    let extracted = extract_haplotypes_from_seq_kbest_paths(
        &paths,
        &seq,
        graph.kmer_size,
        &ref_hap,
        ref_bases.len(),
        sw,
    )
    .unwrap_or_default();
    json!({
        "seq_nodes_before": before.0,
        "seq_edges_before": before.1,
        "seq_nodes": seq.node_count(),
        "seq_edges": seq.edge_count(),
        "cleanup": format!("{status:?}"),
        "n_paths": paths.len(),
        "eq_ref_len": n_eq,
        "n_tc_snp_paths": n_tc,
        "first_tc_ordinal": first_tc,
        "tc_sample": tc_rows,
        "extract": hap_compact(&extracted, pad),
    })
}

#[test]
fn holdout_6r58_chr20_tiny_29456344_tc_vs_tgt() {
    if std::env::var("HOLDOUT_6R58").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R58=1");
        return;
    }
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    assert!(ref_fasta.is_file() && bam.is_file());

    let java_near = vcf_keys_at(&root.join(JAVA_VCF_REL), 29_456_330, 29_456_360);
    let rust_near = vcf_keys_at(&root.join(RUST_VCF_REL), 29_456_330, 29_456_360);
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
            ) && r.start.get() <= POS_SNP
                && r.end.get() >= POS_SNP
        })
        .collect();
    assert!(!covering.is_empty(), "no ActiveFull covers {POS_SNP}");
    let region = covering[0];
    let args = CallRegionArgs::strict_java();
    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);

    let raw_reads: Vec<_> = region.reads.iter().map(|r| r.as_ref().clone()).collect();
    let mut owned = region.clone();
    let assembled =
        assemble_reads_with_finalized(&mut owned, &dict, &mut ref_cache, &args.assemble)
            .expect("assemble");
    let pad = assembled.assembly.padded_reference_start_1based();
    let untrimmed = &assembled.assembly;
    let ref_hap = untrimmed
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .cloned();
    let ref_bytes = untrimmed.apply_bases_shared();
    let off = POS_SNP.saturating_sub(pad) as usize;

    let padded_ref = assembly_reference_read(&dict, &mut ref_cache, region).expect("pad ref");
    let graph_ref = create_graph_reference_read(&padded_ref, region, &dict);
    let graph_reads = records_to_assembly_reads(&assembled.finalized_reads);
    let mut assembler = args.assemble.assembler.clone();
    assembler.dangling_java_exact = true;
    assembler.scoring = Some(AssemblyScoringContext {
        padded_reference_start_1based: region.extended_start.get(),
        active_start_1based: region.start.get(),
        active_end_1based: region.end.get(),
        contig: region.contig.clone(),
    });

    let ref_win = if off + 10 < graph_ref.bases.len() {
        String::from_utf8_lossy(&graph_ref.bases[off.saturating_sub(8)..off + 10]).into_owned()
    } else {
        String::new()
    };
    let mut alt_snp = graph_ref.bases.clone();
    if off < alt_snp.len() {
        alt_snp[off] = b'C';
    }
    let alt25: Vec<Vec<u8>> = {
        let ref_w = kmer_windows_spanning(&graph_ref.bases, off, 25);
        kmer_windows_spanning(&alt_snp, off, 25)
            .into_iter()
            .filter(|km| !ref_w.iter().any(|r| r == km))
            .collect()
    };

    let mut params25 = AssemblyGraphParams::default();
    params25.kmer_size = KmerSize::try_new(25).unwrap();
    params25.min_base_quality = assembler.min_base_quality;
    let raw25 = assembly_graph_from_ref_and_reads_threading(&graph_ref, &graph_reads, &params25)
        .expect("raw25");
    let pruned25 = build_threading_graph_for_seq_assembly(
        &graph_ref,
        &graph_reads,
        25,
        &assembler,
        false,
        false,
    )
    .expect("pruned25");
    let sw = assembler.haplotype_to_reference_sw;
    let seq_kbest = pruned25
        .as_ref()
        .map(|g| seq_kbest_trace(g, graph_ref.bases.as_slice(), pad, off, &sw));

    let rt_before = extract_rt_haplotypes_before_remove_paths(
        &graph_ref,
        &graph_reads,
        &assembler,
        25,
        false,
        false,
    )
    .unwrap_or_default();
    let rt_after = extract_rt_haplotypes_after_remove_paths(
        &graph_ref,
        &graph_reads,
        &assembler,
        25,
        false,
        false,
    )
    .unwrap_or_default();
    let seq_raw = assemble_from_ref_and_reads(&graph_ref, &graph_reads, &assembler).expect("seq");

    let mut per_hap = Vec::new();
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
            let hits = event_hits(&ev, POS_AT_DEL, POS_G);
            if !hits.is_empty() {
                per_hap.push(json!({
                    "i": i,
                    "is_ref": h.is_reference,
                    "cigar": cigar_str(h),
                    "base_344": hap_base_at_ref_locus(h, pad, POS_SNP).map(|b| (b as char).to_string()),
                    "base_345": hap_base_at_ref_locus(h, pad, POS_G).map(|b| (b as char).to_string()),
                    "events": hits,
                }));
            }
        }
    }

    let call = HaplotypeCallerEngine::call_region(region, &dict, &ref_fasta, &args).expect("call");
    let (trim_compact, trim_events, gt_rows, vcf_rows, leftover_events, emit_audit, merged_at_snp) =
        match &call {
            Some(outcome) => {
                let tpad = outcome.assembly.padded_reference_start_1based();
                let compact = hap_compact(&outcome.assembly.haplotypes, tpad);
                let leftover = outcome.assembly.variation_events();
                let tev = event_hits(leftover, POS_AT_DEL, POS_G);
                let pileup_snp = pileup_at(&assembled.finalized_reads, POS_SNP);
                let read_ref_ad = *pileup_snp.get("T").unwrap_or(&0) as i32;
                let read_alt_ad = *pileup_snp.get("C").unwrap_or(&0) as i32;
                let gt: Vec<Value> = outcome
                    .genotyped_calls
                    .iter()
                    .filter(|c| {
                        c.event.start_1based.get() >= POS_AT_DEL
                            && c.event.start_1based.get() <= POS_G
                    })
                    .map(|c| {
                        json!({
                            "pos": c.event.start_1based.get(),
                            "ref": c.event.ref_allele,
                            "alt": c.event.alt_allele,
                            "extra_alts": c.extra_alt_alleles,
                            "n_gl": c.genotype.genotype_log10_likelihoods.len(),
                            "gq": c.genotype.format.gq.as_i32(),
                            "dp": c.genotype.format.dp.as_i32(),
                            "ad": c.genotype.format.ad_as_i32(),
                            "pl": c.genotype.format.pl_as_i32(),
                            "gt_idx": biallelic_genotype_index_from_pl(&c.genotype.format.pl).get(),
                        })
                    })
                    .collect();
                let emit_rows: Vec<Value> = outcome
                    .genotyped_calls
                    .iter()
                    .filter(|c| {
                        c.event.start_1based.get() >= POS_AT_DEL
                            && c.event.start_1based.get() <= POS_G
                    })
                    .map(|c| {
                        let explain = explain_strict_java_emit_gates(
                            &c.event,
                            &c.genotype.genotype_log10_likelihoods,
                            &c.genotype.format,
                            DEFAULT_STAND_EMIT_CONFIDENCE,
                            false,
                            read_ref_ad,
                            read_alt_ad,
                            leftover,
                        )
                        .unwrap_or_else(|_| "err".into());
                        let strict_pass = passes_strict_java_emit_for_genotyped_call(
                            &c.event,
                            &c.genotype.genotype_log10_likelihoods,
                            &c.genotype.format,
                            DEFAULT_STAND_EMIT_CONFIDENCE,
                            false,
                            read_ref_ad,
                            read_alt_ad,
                            false,
                            leftover,
                        )
                        .unwrap_or(false);
                        json!({
                            "pos": c.event.start_1based.get(),
                            "ref": c.event.ref_allele,
                            "alt": c.event.alt_allele,
                            "strict_java_emit_pass": strict_pass,
                            "explain": explain,
                        })
                    })
                    .collect();
                let emitted = try_emit_call_region_variants(
                    region,
                    outcome,
                    "SAMPLE",
                    DEFAULT_STAND_EMIT_CONFIDENCE,
                )
                .unwrap_or_default();
                let vcf: Vec<Value> = emitted
                    .iter()
                    .filter(|r| r.position >= POS_AT_DEL && r.position <= POS_G)
                    .map(|r| {
                        json!({
                            "pos": r.position,
                            "ref": r.reference,
                            "alt": r.alternate,
                            "qual": r.quality,
                            "gt": r.samples.first().and_then(|s| s.gt.as_ref().map(|g| g.alleles.clone())),
                            "ad": r.samples.first().and_then(|s| s.ad.clone()),
                            "pl": r.samples.first().and_then(|s| s.pl.clone()),
                        })
                    })
                    .collect();
                let mut at_snp: Vec<VariationEvent> = leftover
                    .iter()
                    .filter(|e| {
                        e.start_1based.get() == POS_SNP
                            || (e.start_1based.get() <= POS_SNP && e.end_1based.get() >= POS_SNP)
                    })
                    .cloned()
                    .collect();
                let merged_before = merged_biallelic_sites_at_position(&at_snp, POS_SNP);
                prefer_indel_over_colocated_snps(&mut at_snp);
                let merged_after_prefer = merged_biallelic_sites_at_position(&at_snp, POS_SNP);
                let remap =
                    gatk_haplotypecaller::event_map::remap_alt_onto_longer_ref("T", "C", "TG");
                let stripped = remap.as_ref().map(|alt| strip_common_suffix("TG", alt));
                (
                    compact,
                    tev,
                    gt,
                    vcf,
                    event_hits(leftover, POS_AT_DEL, POS_G),
                    emit_rows,
                    json!({
                        "merged_biallelic_from_leftover": merged_before.iter().map(|e| json!({
                            "pos": e.start_1based.get(),
                            "ref": e.ref_allele,
                            "alt": e.alt_allele,
                        })).collect::<Vec<_>>(),
                        "merged_after_prefer_indel": merged_after_prefer.iter().map(|e| json!({
                            "pos": e.start_1based.get(),
                            "ref": e.ref_allele,
                            "alt": e.alt_allele,
                        })).collect::<Vec<_>>(),
                        "remap_T_C_onto_TG": remap,
                        "left_trim_TG_remap": stripped.map(|(r, a)| format!("{r}/{a}")),
                        "java_equivalent_of_TG_CG": "T/C after stripping shared suffix G",
                    }),
                )
            }
            None => (
                json!(null),
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
                json!(null),
            ),
        };

    let union = event_hits(untrimmed.variation_events(), POS_AT_DEL, POS_G);
    let union_has_tc = union
        .iter()
        .any(|e| e["pos"] == POS_SNP && e["ref"] == "T" && e["alt"] == "C");
    let union_has_tgt = union
        .iter()
        .any(|e| e["pos"] == POS_SNP && e["ref"] == "TG" && e["alt"] == "T");
    let mut union_for_merge: Vec<VariationEvent> = untrimmed
        .variation_events()
        .iter()
        .filter(|e| e.start_1based.get() == POS_SNP)
        .cloned()
        .collect();
    let untrimmed_merged = merged_biallelic_sites_at_position(&union_for_merge, POS_SNP);
    prefer_indel_over_colocated_snps(&mut union_for_merge);
    let untrimmed_merged_after_prefer =
        merged_biallelic_sites_at_position(&union_for_merge, POS_SNP);
    let gt_has_tc = gt_rows
        .iter()
        .any(|e| e["pos"] == POS_SNP && e["ref"] == "T" && e["alt"] == "C");
    let gt_has_merged_tg_cg = gt_rows.iter().any(|e| {
        if e["pos"] != POS_SNP || e["ref"] != "TG" {
            return false;
        }
        let extra = e["extra_alts"]
            .as_array()
            .map(|a| {
                a.iter()
                    .filter_map(|x| x.as_str().map(String::from))
                    .collect::<Vec<_>>()
            })
            .unwrap_or_default();
        extra.iter().any(|a| a == "CG") || e["alt"] == "CG"
    });
    let gt_merged_pl: Option<Vec<i64>> = gt_rows.iter().find_map(|e| {
        if e["pos"] != POS_SNP {
            return None;
        }
        // After 6R.63 the genotyped event is reverse-trimmed T/C; PL still from merged subset.
        let is_trim_tc = e["ref"] == "T" && e["alt"] == "C";
        let is_pre_trim_tg = e["ref"] == "TG"
            && (e["alt"] == "CG"
                || e["extra_alts"]
                    .as_array()
                    .is_some_and(|a| a.iter().any(|x| x.as_str() == Some("CG"))));
        if !is_trim_tc && !is_pre_trim_tg {
            return None;
        }
        e["pl"]
            .as_array()
            .map(|a| a.iter().filter_map(|v| v.as_i64()).collect::<Vec<_>>())
    });
    let gt_merged_ad: Option<Vec<i64>> = gt_rows.iter().find_map(|e| {
        if e["pos"] != POS_SNP || e["ref"] != "T" || e["alt"] != "C" {
            return None;
        }
        e["ad"]
            .as_array()
            .map(|a| a.iter().filter_map(|v| v.as_i64()).collect::<Vec<_>>())
    });
    let vcf_tc = vcf_rows.iter().find(|e| {
        e["pos"] == POS_SNP && e["ref"] == "T" && vcf_alt_strings(e).iter().any(|a| a == "C")
    });
    let vcf_has_tc = vcf_rows.iter().any(|e| {
        e["pos"] == POS_SNP && e["ref"] == "T" && vcf_alt_strings(e).iter().any(|a| a == "C")
    });
    let vcf_has_tgt = vcf_rows.iter().any(|e| {
        e["pos"] == POS_SNP && e["ref"] == "TG" && vcf_alt_strings(e).iter().any(|a| a == "T")
    });
    let vcf_has_tg_cg = vcf_rows.iter().any(|e| {
        e["pos"] == POS_SNP && e["ref"] == "TG" && vcf_alt_strings(e).iter().any(|a| a == "CG")
    });
    let per_hap_tc = per_hap.iter().any(|r| {
        r["events"].as_array().is_some_and(|a| {
            a.iter()
                .any(|e| e["pos"] == POS_SNP && e["ref"] == "T" && e["alt"] == "C")
        })
    });
    let assemble_tc = hap_compact(&untrimmed.haplotypes, pad)["n_java_T_to_C_with_G"]
        .as_u64()
        .unwrap_or(0)
        > 0;
    let k128_tc = seq_kbest
        .as_ref()
        .and_then(|v| v["n_tc_snp_paths"].as_u64())
        .unwrap_or(0)
        > 0;
    let pileup_c = pileup_at(&assembled.finalized_reads, POS_SNP)
        .get("C")
        .copied()
        .unwrap_or(0);

    let first = if pileup_c == 0 {
        "A_no_C_in_pileup"
    } else if !k128_tc && !assemble_tc {
        "D_or_E_no_T_to_C_haplotype_in_kbest_or_assemble"
    } else if assemble_tc && !per_hap_tc {
        "F_A_snp_never_constructed"
    } else if per_hap_tc && !union_has_tc {
        "F_B_snp_lost_at_union"
    } else if union_has_tc && !gt_has_tc {
        "H_snp_in_eventmap_not_genotyped"
    } else if gt_has_tc && !vcf_has_tc {
        "I_snp_genotyped_not_emitted"
    } else if vcf_has_tc {
        "emitted_T_C"
    } else {
        "unclassified"
    };

    let doc = json!({
        "k_production": K,
        "representation": {
            "ref_window_around_29456344_8_plus_10": ref_win,
            "ref_29456343_to_46": graph_ref.bases.get(off.saturating_sub(1)..off.saturating_add(3))
                .map(|s| String::from_utf8_lossy(s).into_owned()),
            "java_T_C_means": "substitution T->C at 29456344 keeping G at 29456345 (haplotype ...CACGTCT... vs ref ...CATGTCT...)",
            "rust_TG_T_means": "1bp deletion of G at 29456345, VCF-anchored as TG/T at 29456344",
            "rust_AT_A_nearby": "1bp deletion of T at 29456344, VCF-anchored as AT/A at 29456343",
            "equivalent_representation": false,
            "tg_cg_left_trim_is_T_C": true,
            "tg_t_is_not_T_C": true,
        },
        "java_oracle": {
            "vcf_330_360": java_near,
            "bamout": java_bam,
            "graph_dot": "UNKNOWN",
        },
        "rust_vcf_330_360": rust_near,
        "covering_activefull": {
            "n_covering": covering.len(),
            "active": [region.start.get(), region.end.get()],
            "extended": [region.extended_start.get(), region.extended_end.get()],
            "n_walker_reads": region.reads.len(),
            "n_finalized": assembled.finalized_reads.len(),
            "graph_ref_start": pad,
            "graph_ref_len": graph_ref.bases.len(),
            "target_offset": off,
            "kmer_used": assembled.assembly.kmer_size_for_dump(),
            "kmer_sizes": assembler.kmer_sizes,
            "use_seq_graph": assembler.use_seq_graph,
            "min_prune_factor": assembler.min_prune_factor,
            "num_best_haplotypes": assembler.num_best_haplotypes_per_graph,
            "ref_unique_k25": !reference_has_non_unique_kmers(&graph_ref, 25),
            "min_mq": GATK_HC_DEFAULT_MIN_MAPPING_QUALITY,
        },
        "stage_a_pileup": {
            "raw_29456344": pileup_at(&raw_reads, POS_SNP),
            "finalized_29456344": pileup_at(&assembled.finalized_reads, POS_SNP),
            "finalized_29456345": pileup_at(&assembled.finalized_reads, POS_G),
        },
        "stage_b_rt": {
            "alt25_windows_for_T_to_C": alt25.len(),
            "raw_k25": graph_alt_kmer_stats(&raw25, &alt25),
            "pruned_k25": pruned25.as_ref().map(|g| graph_alt_kmer_stats(g, &alt25)),
            "rt_before_k25": hap_compact(&rt_before, pad),
            "rt_after_k25": hap_compact(&rt_after, pad),
        },
        "stage_c_d_seqgraph": seq_kbest,
        "stage_e_assemble": {
            "assemble_from_ref_and_reads": hap_compact(&seq_raw.haplotypes, pad),
            "untrimmed": hap_compact(&untrimmed.haplotypes, pad),
        },
        "stage_f_eventmap": {
            "union_29456343_to_345": union,
            "per_hap_events_in_window": per_hap,
            "union_has_T_C": union_has_tc,
            "union_has_TG_T": union_has_tgt,
            "per_hap_has_T_C": per_hap_tc,
        },
        "stage_g_trim": {
            "compact": trim_compact,
            "events": trim_events,
        },
        "stage_h_i": {
            "genotyped": gt_rows,
            "vcf": vcf_rows,
            "leftover_variation_events": leftover_events,
            "emit_audit": emit_audit,
            "merged_from_leftover": merged_at_snp,
            "untrimmed_merged_biallelic": untrimmed_merged.iter().map(|e| json!({
                "pos": e.start_1based.get(),
                "ref": e.ref_allele,
                "alt": e.alt_allele,
            })).collect::<Vec<_>>(),
            "untrimmed_merged_after_prefer_indel": untrimmed_merged_after_prefer.iter().map(|e| json!({
                "pos": e.start_1based.get(),
                "ref": e.ref_allele,
                "alt": e.alt_allele,
            })).collect::<Vec<_>>(),
            "gt_has_T_C": gt_has_tc,
            "gt_has_merged_TG_CG": gt_has_merged_tg_cg,
            "vcf_has_T_C": vcf_has_tc,
            "vcf_has_TG_T": vcf_has_tgt,
            "vcf_has_TG_CG": vcf_has_tg_cg,
        },
        "first_loss_of_java_T_C": first,
    });
    println!("{}", serde_json::to_string_pretty(&doc).unwrap());
    assert_eq!(K, 128, "do not raise K");
    assert_eq!(covering.len(), 1);
    assert!(
        java_near.iter().any(|s| s.contains("29456344 T/C")),
        "Java oracle must still contain T/C"
    );
    assert!(
        vcf_has_tc,
        "6R.63: reverseTrimAlleles of subsetted TG/CG must emit T/C"
    );
    assert!(
        !vcf_has_tg_cg,
        "6R.63: untrimmed TG/CG must not remain after reverse trim"
    );
    assert!(
        !vcf_has_tgt,
        "6R.62/6R.63: unused deletion T must be removed; VCF must not keep TG/T"
    );
    assert!(
        gt_has_tc,
        "6R.63: genotyped event is reverse-trimmed T/C (not an independent SNP genotype)"
    );
    let pl = gt_merged_pl.expect("merged-then-trimmed site PL");
    assert_eq!(
        pl.len(),
        3,
        "after unused-ALT subset diploid PL has 3 states, got {pl:?}"
    );
    assert_ne!(
        pl,
        vec![90, 30, 60, 30, 0, 60],
        "must not be the 6R.60 fabricated emit-merge PL: {pl:?}"
    );
    assert_eq!(
        pl,
        vec![266, 0, 1018],
        "6R.84 SPAN_DEL haplotypes no longer dumped into REF; unused-ALT subset + reverse-trim still emit T/C, got {pl:?}"
    );
    let ad = gt_merged_ad.expect("merged-then-trimmed site AD");
    assert_eq!(
        ad,
        vec![26, 9],
        "reverse-trim must not recalculate AD, got {ad:?}"
    );
    let vcf = vcf_tc.expect("VCF T/C row");
    assert_eq!(
        vcf["gt"],
        json!([0, 1]),
        "GT 0/1 must remain associated with the SNP allele after reverse-trim"
    );
    assert_eq!(
        vcf["pl"],
        json!([266, 0, 1018]),
        "VCF PL must still be the unused-ALT subset vector"
    );
    assert_eq!(
        vcf["ad"],
        json!([26, 9]),
        "VCF AD must still be the unused-ALT subset vector unless Java reverse-trim itself changes it"
    );
}

fn vcf_alt_strings(row: &Value) -> Vec<String> {
    match &row["alt"] {
        Value::Array(a) => a
            .iter()
            .filter_map(|x| x.as_str().map(String::from))
            .collect(),
        Value::String(s) => vec![s.clone()],
        _ => Vec::new(),
    }
}
