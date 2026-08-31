//! 6R.54 forensic: follow A-bearing SeqGraph k-best candidates through every
//! GATK 4.4.0.0 `findBestPaths` gate (no production algorithm change).
//!
//! Phenotype: `20:29455902 G/A` on ActiveFull `20:29455745–29455993`.
//! 6R.53 proved the A path is inside production K=128 and absent from untrimmed
//! haplotypes. This dump classifies each candidate against Java retention.
//!
//! Skipped unless `HOLDOUT_6R54=1`.
//!
//! ```text
//! HOLDOUT_6R54=1 cargo test -p gatk-haplotypecaller --test holdout_6r54_chr20_tiny_test -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::assembly_region_finalize::{
    assembly_reference_read, create_graph_reference_read, records_to_assembly_reads,
};
use gatk_haplotypecaller::haplotype::prune_fragment_non_reference_haplotypes;
use gatk_haplotypecaller::haplotype_cigar::trace_find_best_paths_gates;
use gatk_haplotypecaller::hc_allele_mapping::hap_base_at_ref_locus;
use gatk_haplotypecaller::read_threading_assembler::{
    assemble_from_ref_and_reads, build_threading_graph_for_seq_assembly,
    extract_haplotypes_from_seq_kbest_paths, extract_rt_haplotypes_after_remove_paths,
    extract_rt_haplotypes_before_remove_paths, AssemblyScoringContext,
    DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, MIN_HAPLOTYPE_REFERENCE_LENGTH,
};
use gatk_haplotypecaller::seq_graph::SeqGraph;
use gatk_haplotypecaller::seq_kbest_haplotype::find_best_haplotypes_seq_graph;
use gatk_haplotypecaller::{
    assemble_reads_with_finalized, call_disposition, flatten_assembly_regions,
    query_index_at_reference_position, traverse_assembly_region_walker, AssemblyRegion,
    AssemblyRegionCallDisposition, CallRegionArgs, Cigar, CigarOperator, Haplotype,
    ReadFilterParams, WalkerTraversalConfig,
};
use rust_htslib::bam::record::CigarString;
use rust_htslib::bam::Read;
use serde_json::{json, Value};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const ACTIVE: (u64, u64) = (29_455_745, 29_455_993);
const TARGET: u64 = 29_455_902;
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
        let b = rec.seq()[qi] as char;
        *counts.entry(b.to_ascii_uppercase()).or_default() += 1;
    }
    counts
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
        "sha8": format!("{:08x}", {
            use std::hash::{Hash, Hasher};
            let mut hasher = std::collections::hash_map::DefaultHasher::new();
            h.bases.hash(&mut hasher);
            hasher.finish() as u32
        }),
    })
}

fn seq_kbest_trace(
    seq: &SeqGraph,
    ref_hap: &Haplotype,
    pad: u64,
    off: usize,
    kmer: usize,
    sw: &gatk_haplotypecaller::SwParameters,
) -> Value {
    let ref_bytes = ref_hap.bases.as_slice();
    let ref_cigar_len = ref_hap
        .cigar
        .as_ref()
        .map(|c| c.reference_length())
        .unwrap_or(ref_bytes.len());
    let paths = find_best_haplotypes_seq_graph(seq, K).unwrap_or_default();
    let mut seen: Vec<(Vec<u8>, bool)> = Vec::new();
    let mut a_bearing = Vec::new();
    let mut retained_non_a = Vec::new();
    for (i, p) in paths.iter().enumerate() {
        let bases = seq.path_bases_bytes(p.start, &p.edges);
        let eq_len = bases.len() == ref_bytes.len();
        let has_a = eq_len && bases.get(off) == Some(&TARGET_ALT);
        let t = trace_find_best_paths_gates(
            ref_bytes,
            &bases,
            p.is_reference,
            ref_cigar_len,
            sw,
            &seen,
        );
        if t.rust_extract_keep {
            seen.push((bases.clone(), p.is_reference));
        }
        let row = json!({
            "kbest_ordinal": i,
            "score": p.score,
            "is_reference_flag": p.is_reference,
            "seq_len": t.seq_len,
            "ref_hap_len": t.ref_hap_len,
            "eq_ref_len": eq_len,
            "base_at_offset": bases.get(off).copied().map(|b| (b as char).to_string()),
            "rust_prod_cigar": t.rust_prod_cigar,
            "rust_prod_ref_len": t.rust_prod_ref_len,
            "java_softclip_cigar": t.java_softclip_cigar,
            "java_softclip_ref_len": t.java_softclip_ref_len,
            "java_indel_cigar": t.java_indel_cigar,
            "java_indel_ref_len": t.java_indel_ref_len,
            "cigar_contains_n": t.cigar_contains_n,
            "rust_prod_spans_required": t.rust_prod_spans_required,
            "java_softclip_spans_required": t.java_softclip_spans_required,
            "duplicate": t.duplicate,
            "rust_extract_keep": t.rust_extract_keep,
            "java_would_retain": t.java_would_retain,
            "first_rust_reject": t.first_rust_reject,
            "first_java_reject": t.first_java_reject,
        });
        if has_a {
            a_bearing.push(row);
        } else if t.rust_extract_keep && !p.is_reference {
            if retained_non_a.len() < 8 {
                retained_non_a.push(row);
            }
        }
    }

    let extracted =
        extract_haplotypes_from_seq_kbest_paths(&paths, seq, kmer, ref_hap, ref_cigar_len, sw)
            .unwrap_or_default();
    let extract_has_a = extracted
        .iter()
        .any(|h| h.bases.len() == ref_bytes.len() && h.bases.get(off) == Some(&TARGET_ALT));
    let mut pruned = extracted.clone();
    prune_fragment_non_reference_haplotypes(&mut pruned, ref_hap, MIN_HAPLOTYPE_REFERENCE_LENGTH);
    let prune_has_a = pruned
        .iter()
        .any(|h| h.bases.len() == ref_bytes.len() && h.bases.get(off) == Some(&TARGET_ALT));

    json!({
        "kmer": kmer,
        "n_paths": paths.len(),
        "n_a_bearing_eq_len": a_bearing.len(),
        "first_a_ordinal": a_bearing.first().and_then(|v| v.get("kbest_ordinal").cloned()),
        "a_bearing": a_bearing,
        "retained_non_a_sample": retained_non_a,
        "extract": {
            "n": extracted.len(),
            "has_a": extract_has_a,
            "haps": extracted.iter().map(|h| hap_row(h, pad)).collect::<Vec<_>>(),
        },
        "after_prune_fragment": {
            "n": pruned.len(),
            "has_a": prune_has_a,
        },
    })
}

fn rt_trace(
    graph_ref: &gatk_haplotypecaller::assembly::AssemblyRead,
    graph_reads: &[gatk_haplotypecaller::assembly::AssemblyRead],
    assembler: &gatk_haplotypecaller::ReadThreadingAssemblerArgs,
    pad: u64,
    kmer: usize,
) -> Value {
    let before = extract_rt_haplotypes_before_remove_paths(
        graph_ref,
        graph_reads,
        assembler,
        kmer,
        false,
        false,
    )
    .unwrap_or_default();
    let after = extract_rt_haplotypes_after_remove_paths(
        graph_ref,
        graph_reads,
        assembler,
        kmer,
        false,
        false,
    )
    .unwrap_or_default();
    json!({
        "kmer": kmer,
        "before_remove": {
            "n": before.len(),
            "has_a": before.iter().any(|h| hap_base_at_ref_locus(h, pad, TARGET) == Some(TARGET_ALT)),
            "haps": before.iter().map(|h| hap_row(h, pad)).collect::<Vec<_>>(),
        },
        "after_remove": {
            "n": after.len(),
            "has_a": after.iter().any(|h| hap_base_at_ref_locus(h, pad, TARGET) == Some(TARGET_ALT)),
            "haps": after.iter().map(|h| hap_row(h, pad)).collect::<Vec<_>>(),
        },
    })
}

fn assemble_summary(label: &str, haps: &[Haplotype], pad: u64) -> Value {
    json!({
        "label": label,
        "n": haps.len(),
        "has_a": haps.iter().any(|h| hap_base_at_ref_locus(h, pad, TARGET) == Some(TARGET_ALT)),
        "haps": haps.iter().map(|h| hap_row(h, pad)).collect::<Vec<_>>(),
    })
}

fn java_bamout_observe(root: &Path) -> Value {
    // Prefer a region-specific bamout if a prior round left one; otherwise UNKNOWN.
    let candidates = [
        "parity/reports/6r54/java_active.bamout.bam",
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
        let mut hc: BTreeMap<String, (String, usize, bool)> = BTreeMap::new();
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
            let has_a = query_index_at_reference_position(rec.pos(), &cigar, pos0)
                .map(|qi| rec.seq()[qi])
                .is_some_and(|b| b.eq_ignore_ascii_case(&TARGET_ALT));
            let cstr = cigar.to_string();
            hc.entry(q)
                .and_modify(|e| {
                    e.2 |= has_a;
                })
                .or_insert((cstr, rec.seq_len(), has_a));
        }
        let a_haps: Vec<_> = hc
            .iter()
            .filter(|(_, v)| v.2)
            .map(|(k, v)| json!({"qname": k, "cigar": v.0, "len": v.1}))
            .collect();
        return json!({
            "source": rel,
            "n_hc_qnames": hc.len(),
            "n_a_bearing_hc": a_haps.len(),
            "a_bearing": a_haps,
            "all_hc": hc.iter().map(|(k, v)| json!({
                "qname": k, "cigar": v.0, "len": v.1, "has_a": v.2
            })).collect::<Vec<_>>(),
        });
    }
    json!({"status": "UNKNOWN", "reason": "no java bamout at expected paths"})
}

fn trace_region(region: &AssemblyRegion, dict: &SequenceDictionary, ref_fasta: &Path) -> Value {
    let args = CallRegionArgs::strict_java();
    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.to_path_buf(), 4);
    let mut owned = region.clone();
    let assembled = assemble_reads_with_finalized(&mut owned, dict, &mut ref_cache, &args.assemble)
        .expect("assemble");
    let pad = assembled.assembly.padded_reference_start_1based();
    let padded_ref = assembly_reference_read(dict, &mut ref_cache, region).expect("pad ref");
    let graph_ref = create_graph_reference_read(&padded_ref, region, dict);
    let graph_reads = records_to_assembly_reads(&assembled.finalized_reads);

    let mut assembler = args.assemble.assembler.clone();
    assembler.dangling_java_exact = true;
    assembler.scoring = Some(AssemblyScoringContext {
        padded_reference_start_1based: region.extended_start.get(),
        active_start_1based: region.start.get(),
        active_end_1based: region.end.get(),
        contig: region.contig.clone(),
    });

    let prod_like =
        assemble_from_ref_and_reads(&graph_ref, &graph_reads, &assembler).expect("prod");
    let mut seq_only_args = assembler.clone();
    seq_only_args.scoring = None;
    let seq_only =
        assemble_from_ref_and_reads(&graph_ref, &graph_reads, &seq_only_args).expect("seq-only");

    let mut ref_hap = Haplotype::new(graph_ref.bases.as_slice(), true);
    let mut rc = Cigar::new();
    rc.push(ref_hap.bases.len(), CigarOperator::Match);
    ref_hap.cigar = Some(rc);
    let off = TARGET.saturating_sub(pad) as usize;
    let sw = assembler.haplotype_to_reference_sw;

    let mut seq_by_k = Vec::new();
    for kmer in [10usize, 25] {
        if let Ok(Some(g)) = build_threading_graph_for_seq_assembly(
            &graph_ref,
            &graph_reads,
            kmer,
            &assembler,
            false,
            false,
        ) {
            let mut seq = SeqGraph::from_assembly_graph(&g);
            let _ = seq.cleanup_seq_graph();
            seq_by_k.push(json!({
                "kmer": kmer,
                "nodes": seq.node_count(),
                "edges": seq.edge_count(),
                "kbest_gates": seq_kbest_trace(&seq, &ref_hap, pad, off, kmer, &sw),
            }));
        } else {
            seq_by_k.push(json!({"kmer": kmer, "built": false}));
        }
    }

    let rt10 = rt_trace(&graph_ref, &graph_reads, &assembler, pad, 10);
    let rt25 = rt_trace(&graph_ref, &graph_reads, &assembler, pad, 25);

    json!({
        "active": [region.start.get(), region.end.get()],
        "extended": [region.extended_start.get(), region.extended_end.get()],
        "graph_ref_start": pad,
        "graph_ref_len": graph_ref.bases.len(),
        "target_offset": off,
        "pileup_finalized": pileup_at(&assembled.finalized_reads, TARGET),
        "production_assemble_reads_with_finalized": assemble_summary(
            "assemble_reads_with_finalized",
            &assembled.assembly.haplotypes,
            pad,
        ),
        "assemble_from_ref_scoring_on": assemble_summary(
            "scoring_on_rt_first_eligible",
            &prod_like.haplotypes,
            pad,
        ),
        "assemble_from_ref_scoring_off": assemble_summary(
            "scoring_off_seqgraph_path",
            &seq_only.haplotypes,
            pad,
        ),
        "seqgraph_kbest": seq_by_k,
        "rt_kbest": [rt10, rt25],
    })
}

#[test]
fn holdout_6r54_find_best_paths_a_bearing_gates() {
    if std::env::var("HOLDOUT_6R54").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R54=1");
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
    let region = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() == ACTIVE.0
                && r.end.get() == ACTIVE.1
        })
        .expect("ActiveFull 29455745-29455993");

    let rust = trace_region(region, &dict, &ref_fasta);
    let java_obs = java_bamout_observe(&root);
    let doc = json!({
        "k_production": K,
        "java_observability": java_obs,
        "rust": rust,
    });
    println!("{}", serde_json::to_string_pretty(&doc).unwrap());
    assert_eq!(K, 128, "do not raise K");
}
