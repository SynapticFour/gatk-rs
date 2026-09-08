//! 6R.123 forensic dump: haplotype likelihoods → `marginalize` at `20:29456196`.
//! Does not recompute PairHMM as a separate investigation; consumes `call_region`
//! likelihoods already produced for genotyping. Skipped unless `HOLDOUT_6R123=1`.
//!
//! ```text
//! HOLDOUT_6R123=1 cargo test -p gatk-haplotypecaller --test holdout_6r123_marginalize -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{build_per_haplotype_variation_events, VariationEvent};
use gatk_haplotypecaller::hc_allele_mapping::create_allele_mapper_with_events;
use gatk_haplotypecaller::hc_genotyping_engine::{
    marginalize_rows_to_biallelic_alleles, region_likelihoods_to_rows,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    AssemblyRegionCallDisposition, CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams,
    WalkerTraversalConfig,
};
use serde_json::{json, Value};
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 29_456_196;
const EMIT_SPANNING_DELS: bool = true;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn fnv1a64_hex(bases: &[u8]) -> String {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in bases {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    format!("{h:016x}")
}

fn read_key(qname: &str, flags: u16, start_1based: i64) -> String {
    format!("{qname}\tflags={flags}\tstart={start_1based}")
}

#[test]
fn holdout_6r123_marginalize() {
    if std::env::var("HOLDOUT_6R123").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R123=1");
        return;
    }
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
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
    let covering = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= TARGET
                && r.end.get() >= TARGET
        })
        .expect("ActiveFull");
    let args = CallRegionArgs::strict_java();
    let outcome = HaplotypeCallerEngine::call_region(covering, &dict, &ref_fasta, &args)
        .expect("call")
        .expect("outcome");
    let pad = outcome.assembly.padded_reference_start_1based();
    let ref_bytes = outcome.assembly.reference_bases();
    let haps = &outcome.assembly.haplotypes;
    let hap_cache = build_per_haplotype_variation_events(
        haps,
        ref_bytes,
        pad,
        outcome.assembly.max_mnp_distance(),
        "20",
    );
    let merged = VariationEvent::from_alleles("20", TARGET, "A", "T");
    let mapping = create_allele_mapper_with_events(
        &merged,
        TARGET,
        haps,
        pad,
        ref_bytes,
        outcome.assembly.max_mnp_distance(),
        EMIT_SPANNING_DELS,
        Some(&hap_cache),
    );
    let ref_set: BTreeSet<usize> = mapping
        .ref_haplotype_indices
        .iter()
        .map(|i| i.get())
        .collect();
    let alt_set: BTreeSet<usize> = mapping
        .alt_haplotype_indices
        .iter()
        .map(|i| i.get())
        .collect();
    let hap_hashes: Vec<String> = haps.iter().map(|h| fnv1a64_hex(&h.bases)).collect();
    let mut rows = region_likelihoods_to_rows(&outcome.read_likelihoods, haps.len());
    for row in &mut rows {
        if let Some(rec) = outcome.genotyping_reads.get(row.read_index) {
            let qname = String::from_utf8_lossy(rec.qname()).into_owned();
            row.read_id = read_key(&qname, rec.flags(), rec.pos() + 1);
        }
    }
    let marg = marginalize_rows_to_biallelic_alleles(
        &rows,
        &mapping.ref_haplotype_indices,
        &mapping.alt_haplotype_indices,
    );
    let hap_cols: Vec<Value> = haps
        .iter()
        .enumerate()
        .map(|(i, h)| {
            json!({
                "idx": i,
                "hash": hap_hashes[i],
                "len": h.bases.len(),
                "is_reference": h.is_reference,
                "mapped": if alt_set.contains(&i) { "ALT" } else if ref_set.contains(&i) { "REF" } else { "unmapped" },
            })
        })
        .collect();
    let pre: Vec<Value> = rows
        .iter()
        .map(|r| {
            json!({
                "read": r.read_id,
                "ll": r.haplotype_log10_likelihoods,
            })
        })
        .collect();
    let post: Vec<Value> = marg
        .iter()
        .map(|r| {
            json!({
                "read": r.read_id,
                "A": r.haplotype_log10_likelihoods.first().copied(),
                "T": r.haplotype_log10_likelihoods.get(1).copied(),
            })
        })
        .collect();
    let doc = json!({
        "target": TARGET,
        "active": [covering.start.get(), covering.end.get()],
        "input": {
            "hap_count": haps.len(),
            "evidence_count": rows.len(),
            "likelihood_cells": outcome.read_likelihoods.len(),
            "mapper_ref": mapping.ref_haplotype_indices.len(),
            "mapper_alt": mapping.alt_haplotype_indices.len(),
            "reduction": "max_log10_per_allele_pool",
            "empty_pool": -50.0,
        },
        "hap_cols": hap_cols,
        "pre_marginalize": pre,
        "post_marginalize": post,
    });
    println!("{}", serde_json::to_string_pretty(&doc).unwrap());
}
