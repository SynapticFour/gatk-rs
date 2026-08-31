//! 6R.53 forensic: remaining chr20_tiny allele-set divergences after 6R.52.
//!
//! Covering ActiveFull `20:29455300–29455559` VCF allele set already matches Java.
//! This dump traces leftover EventMap extras there, then the first Java-only
//! ActiveFull (`20:29455745–29455993`).
//!
//! Skipped unless `HOLDOUT_6R53=1`. No production algorithm change.
//!
//! ```text
//! HOLDOUT_6R53=1 cargo test -p gatk-haplotypecaller --test holdout_6r53_chr20_tiny_test -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::assembly_region_finalize::{
    assembly_reference_read, create_graph_reference_read, records_to_assembly_reads,
};
use gatk_haplotypecaller::hc_allele_mapping::hap_base_at_ref_locus;
use gatk_haplotypecaller::hc_genotyping_engine::DEFAULT_STAND_EMIT_CONFIDENCE;
use gatk_haplotypecaller::read_threading_assembler::{
    build_threading_graph_for_seq_assembly, DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH,
};
use gatk_haplotypecaller::seq_graph::SeqGraph;
use gatk_haplotypecaller::seq_kbest_haplotype::find_best_haplotypes_seq_graph;
use gatk_haplotypecaller::{
    assemble_reads_with_finalized, call_disposition, flatten_assembly_regions,
    query_index_at_reference_position, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegion, AssemblyRegionCallDisposition, CallRegionArgs,
    Haplotype, HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
};
use rust_htslib::bam::record::CigarString;
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_VCF_REL: &str = "parity/reports/6r43/chr20_tiny/java.vcf";
const RUST_VCF_REL: &str = "parity/reports/6r43/chr20_tiny/rust.vcf";
const COVERING: (u64, u64) = (29_455_300, 29_455_559);
const JAVA_ONLY_REGION: (u64, u64) = (29_455_745, 29_455_993);
const K: usize = DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn parse_vcf_keys(path: &Path) -> BTreeSet<(u64, String, String)> {
    let mut out = BTreeSet::new();
    if !path.is_file() {
        return out;
    }
    for line in fs::read_to_string(path).unwrap_or_default().lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<_> = line.split('\t').collect();
        if f.len() < 5 || f[0] != "20" {
            continue;
        }
        out.insert((f[1].parse().unwrap(), f[3].to_string(), f[4].to_string()));
    }
    out
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

fn hap_base_counts(haps: &[Haplotype], pad: u64, loc: u64) -> BTreeMap<String, usize> {
    let mut c: BTreeMap<String, usize> = BTreeMap::new();
    for h in haps {
        let key = hap_base_at_ref_locus(h, pad, loc)
            .map(|b| (b as char).to_string())
            .unwrap_or_else(|| ".".to_string());
        *c.entry(key).or_default() += 1;
    }
    c
}

fn seq_kbest_alt_at(seq: &SeqGraph, ref_bases: &[u8], off: usize, k: usize) -> Value {
    let paths = find_best_haplotypes_seq_graph(seq, k).unwrap_or_default();
    let mut eq_alt = 0usize;
    let mut first = None;
    let mut eq_len = 0usize;
    for (i, p) in paths.iter().enumerate() {
        let b = seq.path_bases_bytes(p.start, &p.edges);
        if b.len() != ref_bases.len() {
            continue;
        }
        eq_len += 1;
        if b.get(off) == Some(&b'A') {
            eq_alt += 1;
            if first.is_none() {
                first = Some(i);
            }
        }
    }
    json!({
        "k_best": k,
        "n_paths": paths.len(),
        "eq_ref_len": eq_len,
        "index_alt_a": eq_alt,
        "first_index_alt_a": first,
    })
}

fn trace_region(
    region: &AssemblyRegion,
    dict: &SequenceDictionary,
    ref_fasta: &Path,
    args: &CallRegionArgs,
    target: u64,
    target_ref: u8,
    target_alt: u8,
    label: &str,
) -> Value {
    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.to_path_buf(), 4);
    let mut owned = region.clone();
    let assembled = assemble_reads_with_finalized(&mut owned, dict, &mut ref_cache, &args.assemble)
        .expect("assemble");
    let pad = assembled.assembly.padded_reference_start_1based();
    let untrimmed = &assembled.assembly;
    let raw_reads: Vec<_> = region.reads.iter().map(|r| r.as_ref().clone()).collect();
    let padded_ref = assembly_reference_read(dict, &mut ref_cache, region).expect("pad ref");
    let graph_ref = create_graph_reference_read(&padded_ref, region, dict);
    let graph_reads = records_to_assembly_reads(&assembled.finalized_reads);
    let assembler = args.assemble.assembler.clone();
    let off = target.saturating_sub(pad) as usize;

    let mut seq_kbest = Value::Null;
    let mut seq_nodes = 0usize;
    let mut seq_edges = 0usize;
    if off < graph_ref.bases.len() {
        if let Ok(Some(g)) = build_threading_graph_for_seq_assembly(
            &graph_ref,
            &graph_reads,
            25,
            &assembler,
            false,
            false,
        ) {
            let mut seq = SeqGraph::from_assembly_graph(&g);
            let _ = seq.cleanup_seq_graph();
            seq_nodes = seq.node_count();
            seq_edges = seq.edge_count();
            seq_kbest = json!({
                "k128": seq_kbest_alt_at(&seq, graph_ref.bases.as_slice(), off, K),
                "k256": seq_kbest_alt_at(&seq, graph_ref.bases.as_slice(), off, 256),
            });
        }
    }

    let events: Vec<Value> = untrimmed
        .variation_events()
        .iter()
        .filter(|e| e.start_1based.get() == target)
        .map(|e| {
            json!({
                "pos": e.start_1based.get(),
                "ref": e.ref_allele,
                "alt": e.alt_allele,
            })
        })
        .collect();

    let call = HaplotypeCallerEngine::call_region(region, dict, ref_fasta, args).expect("call");
    let (trim_n, trim_has, trim_events, vcf_has, vcf_n) = match &call {
        Some(outcome) => {
            let tpad = outcome.assembly.padded_reference_start_1based();
            let has = outcome
                .assembly
                .haplotypes
                .iter()
                .any(|h| hap_base_at_ref_locus(h, tpad, target) == Some(target_alt));
            let tev: Vec<Value> = outcome
                .assembly
                .variation_events()
                .iter()
                .filter(|e| e.start_1based.get() == target)
                .map(|e| {
                    json!({
                        "pos": e.start_1based.get(),
                        "ref": e.ref_allele,
                        "alt": e.alt_allele,
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
            let vcf_hit = emitted.iter().any(|r| {
                r.position == target
                    && r.reference.as_bytes() == [target_ref]
                    && r.alternate.iter().any(|a| a.as_bytes() == [target_alt])
            });
            (
                outcome.assembly.haplotypes.len(),
                has,
                tev,
                vcf_hit,
                emitted.len(),
            )
        }
        None => (0, false, Vec::new(), false, 0),
    };

    json!({
        "label": label,
        "active": [region.start.get(), region.end.get()],
        "extended": [region.extended_start.get(), region.extended_end.get()],
        "n_reads": region.reads.len(),
        "n_finalized": assembled.finalized_reads.len(),
        "kmer": assembled.assembly.kmer_size_for_dump(),
        "graph_ref_start": pad,
        "graph_ref_len": graph_ref.bases.len(),
        "target": {"pos": target, "ref": (target_ref as char).to_string(), "alt": (target_alt as char).to_string(), "offset": off},
        "pileup_raw": pileup_at(&raw_reads, target),
        "pileup_finalized": pileup_at(&assembled.finalized_reads, target),
        "seqgraph": {"nodes": seq_nodes, "edges": seq_edges},
        "kbest": seq_kbest,
        "untrimmed": {
            "n_haps": untrimmed.haplotypes.len(),
            "bases_at_target": hap_base_counts(&untrimmed.haplotypes, pad, target),
            "eventmap_at_target": events,
        },
        "trimmed": {
            "n_haps": trim_n,
            "has_alt": trim_has,
            "eventmap_at_target": trim_events,
        },
        "vcf_emit_has_snp": vcf_has,
        "vcf_n_records": vcf_n,
    })
}

#[test]
fn holdout_6r53_chr20_tiny_remaining_allele_sets() {
    if std::env::var("HOLDOUT_6R53").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R53=1");
        return;
    }
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    let java_vcf = root.join(JAVA_VCF_REL);
    let rust_vcf = root.join(RUST_VCF_REL);
    assert!(ref_fasta.is_file(), "missing {}", ref_fasta.display());
    assert!(bam.is_file(), "missing {}", bam.display());

    let jk = parse_vcf_keys(&java_vcf);
    let rk = parse_vcf_keys(&rust_vcf);
    let in_span = |lo: u64, hi: u64, k: &(u64, String, String)| k.0 >= lo && k.0 <= hi;
    let cov_j: BTreeSet<_> = jk
        .iter()
        .filter(|k| in_span(COVERING.0, COVERING.1, k))
        .cloned()
        .collect();
    let cov_r: BTreeSet<_> = rk
        .iter()
        .filter(|k| in_span(COVERING.0, COVERING.1, k))
        .cloned()
        .collect();
    let jo_j: BTreeSet<_> = jk
        .iter()
        .filter(|k| in_span(JAVA_ONLY_REGION.0, JAVA_ONLY_REGION.1, k))
        .cloned()
        .collect();
    let jo_r: BTreeSet<_> = rk
        .iter()
        .filter(|k| in_span(JAVA_ONLY_REGION.0, JAVA_ONLY_REGION.1, k))
        .cloned()
        .collect();

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
    let args = CallRegionArgs::strict_java();

    let covering = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() == COVERING.0
                && r.end.get() == COVERING.1
        })
        .expect("covering ActiveFull");
    let java_only_reg = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() == JAVA_ONLY_REGION.0
                && r.end.get() == JAVA_ONLY_REGION.1
        })
        .expect("java-only ActiveFull");

    let covering_json = trace_region(
        covering,
        &dict,
        &ref_fasta,
        &args,
        29_455_379,
        b'G',
        b'A',
        "covering_29455379",
    );
    let extra_375 = trace_region(
        covering,
        &dict,
        &ref_fasta,
        &args,
        29_455_375,
        b'T',
        b'A',
        "covering_eventmap_extra_29455375",
    );
    let java_only_902 = trace_region(
        java_only_reg,
        &dict,
        &ref_fasta,
        &args,
        29_455_902,
        b'G',
        b'A',
        "java_only_29455902",
    );

    let doc = json!({
        "k_production": K,
        "vcf_full_tiny": {
            "java_n": jk.len(),
            "rust_n": rk.len(),
            "shared": jk.intersection(&rk).count(),
            "java_only": jk.difference(&rk).cloned().collect::<Vec<_>>(),
            "rust_only": rk.difference(&jk).cloned().collect::<Vec<_>>(),
        },
        "covering_vcf": {
            "java_n": cov_j.len(),
            "rust_n": cov_r.len(),
            "java_only": cov_j.difference(&cov_r).cloned().collect::<Vec<_>>(),
            "rust_only": cov_r.difference(&cov_j).cloned().collect::<Vec<_>>(),
        },
        "first_java_only_span_vcf": {
            "active": JAVA_ONLY_REGION,
            "java": jo_j.iter().cloned().collect::<Vec<_>>(),
            "rust": jo_r.iter().cloned().collect::<Vec<_>>(),
            "java_only": jo_j.difference(&jo_r).cloned().collect::<Vec<_>>(),
            "rust_only": jo_r.difference(&jo_j).cloned().collect::<Vec<_>>(),
        },
        "covering_29455379": covering_json,
        "covering_29455375": extra_375,
        "java_only_29455902": java_only_902,
    });
    println!("{}", serde_json::to_string_pretty(&doc).unwrap());

    assert!(
        cov_j == cov_r,
        "6R.52: covering VCF allele set must still match Java; java_only={:?} rust_only={:?}",
        cov_j.difference(&cov_r).collect::<Vec<_>>(),
        cov_r.difference(&cov_j).collect::<Vec<_>>()
    );
    assert_eq!(K, 128, "do not raise K");
}
