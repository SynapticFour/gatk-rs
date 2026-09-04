//! 6R.68 forensic: live read–allele likelihood matrix at the canonical T/C site.
//!
//! Skipped unless `HOLDOUT_6R68=1`. The coordinate-free proof lives in
//! `forensic_6r68_pairhmm_inputs.rs`; this harness only records region facts.
//!
//! ```text
//! HOLDOUT_6R68=1 cargo test -p gatk-haplotypecaller --test holdout_6r68_read_allele_likelihood -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, take_colocated_merge_numerics,
    traverse_assembly_region_walker, try_emit_call_region_variants, AssemblyRegionCallDisposition,
    CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
    DEFAULT_STAND_EMIT_CONFIDENCE, GATK_PARITY_DEFAULT_INS_QUAL,
};
use rust_htslib::bam::record::Aux;
use serde_json::json;
use std::collections::HashSet;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const POS_SNP: u64 = 29_456_344;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn aux_phred_z(rec: &rust_htslib::bam::Record, tag: &[u8]) -> Option<Vec<u8>> {
    match rec.aux(tag) {
        Ok(Aux::String(s)) => Some(s.bytes().map(|b| b.saturating_sub(33)).collect()),
        _ => None,
    }
}

#[test]
fn holdout_6r68_read_allele_likelihood_29456344() {
    if std::env::var("HOLDOUT_6R68").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R68=1");
        return;
    }
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    assert!(ref_fasta.is_file() && bam.is_file());

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
    let outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");

    let emitted =
        try_emit_call_region_variants(region, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    let vcf = emitted
        .iter()
        .find(|r| {
            r.position == POS_SNP && r.reference == "T" && r.alternate.iter().any(|a| a == "C")
        })
        .expect("lifecycle T/C");
    let pl = vcf.samples[0].pl.clone().expect("PL");
    assert_eq!(pl, vec![266, 0, 1018], "6R.84 SPAN_DEL not dumped into REF");

    let live = take_colocated_merge_numerics();
    let numerics = live
        .iter()
        .find(|n| n.loc == POS_SNP)
        .cloned()
        .expect("merge numerics");
    assert_eq!(
        numerics.alts,
        vec!["T".to_string(), "CG".to_string(), "*".to_string()]
    );
    assert_eq!(numerics.pool_sizes.len(), 4);
    assert_eq!(numerics.pool_sizes[1], 6);
    assert_eq!(numerics.pool_sizes[2], 21);
    assert_eq!(numerics.pool_sizes[3], 6);

    let n_pairhmm_reads: HashSet<usize> = outcome
        .read_likelihoods
        .iter()
        .map(|c| c.read_index.get())
        .collect();
    let n_haps: HashSet<usize> = outcome
        .read_likelihoods
        .iter()
        .map(|c| c.haplotype_index.get())
        .collect();

    let mut n_bi = 0usize;
    let mut n_bd = 0usize;
    let mut n_bi_ne_45_reads = 0usize;
    let mut n_bd_ne_45_reads = 0usize;
    let mut n_pre_pcr_ins_mismatch_reads = 0usize;
    let mut n_pre_pcr_ins_mismatch_bases = 0usize;
    let mut n_gop_bases = 0usize;
    let mut first_mismatch_read: Option<usize> = None;
    let mut first_mismatch_base: Option<(usize, u8, u8)> = None;
    for (ri, rec) in outcome.genotyping_reads.iter().enumerate() {
        let bi = aux_phred_z(rec, b"BI");
        let bd = aux_phred_z(rec, b"BD");
        if bi.is_some() {
            n_bi += 1;
        }
        if bd.is_some() {
            n_bd += 1;
        }
        let bases = rec.seq().as_bytes();
        if bases.is_empty() {
            continue;
        }
        n_gop_bases += bases.len();
        if let Some(ref bi) = bi {
            if bi.iter().any(|&q| q != 45) {
                n_bi_ne_45_reads += 1;
            }
            let rust_ins = vec![GATK_PARITY_DEFAULT_INS_QUAL; bases.len()];
            let n = rust_ins.len().min(bi.len());
            let mismatches = rust_ins[..n]
                .iter()
                .zip(bi[..n].iter())
                .filter(|(a, b)| a != b)
                .count();
            if mismatches > 0 {
                n_pre_pcr_ins_mismatch_reads += 1;
                n_pre_pcr_ins_mismatch_bases += mismatches;
                if first_mismatch_read.is_none() {
                    first_mismatch_read = Some(ri);
                    if let Some(i) = rust_ins[..n]
                        .iter()
                        .zip(bi[..n].iter())
                        .position(|(a, b)| a != b)
                    {
                        first_mismatch_base = Some((i, rust_ins[i], bi[i]));
                    }
                }
            }
        }
        if let Some(ref bd) = bd {
            if bd.iter().any(|&q| q != 45) {
                n_bd_ne_45_reads += 1;
            }
        }
    }

    let doc = json!({
        "locus": "20:29456344 T/C",
        "vcf": {
            "ref": vcf.reference,
            "alt": vcf.alternate,
            "gt": vcf.samples[0].gt.as_ref().map(|g| g.alleles.clone()),
            "pl": pl,
        },
        "java_oracle_pl": [542, 0, 1353],
        "matrix": {
            "n_genotyping_reads": outcome.genotyping_reads.len(),
            "n_pairhmm_read_indices": n_pairhmm_reads.len(),
            "n_haplotype_indices": n_haps.len(),
            "n_assembly_haps": outcome.assembly.haplotypes.len(),
            "n_likelihood_cells": outcome.read_likelihoods.len(),
            "n_overlap_reads": numerics.n_reads,
            "n_pairhmm_reads_snap": numerics.n_pairhmm_reads,
            "n_overlap_before_qname_dedupe": numerics.n_overlap_before_qname_dedupe,
            "pool_sizes": numerics.pool_sizes,
            "merged_pl": numerics.merged_pl,
        },
        "layer_a_indel_gop_source": {
            "n_reads_with_bi_tag": n_bi,
            "n_reads_with_bd_tag": n_bd,
            "n_bi_not_all_45": n_bi_ne_45_reads,
            "n_bd_not_all_45": n_bd_ne_45_reads,
            "n_pre_pcr_ins_mismatch_reads": n_pre_pcr_ins_mismatch_reads,
            "n_pre_pcr_ins_mismatch_bases": n_pre_pcr_ins_mismatch_bases,
            "n_gop_bases": n_gop_bases,
            "first_mismatch_read_index": first_mismatch_read,
            "first_mismatch_base_rust_vs_bi": first_mismatch_base,
        },
        "first_divergence": "PairHMM input construction (indel GOP source: BAM BI/BD vs constant Q45)",
    });
    println!("{}", serde_json::to_string_pretty(&doc).unwrap());

    assert!(n_bi > 0, "canonical BAM carries BI tags");
    assert!(n_bd > 0, "canonical BAM carries BD tags");
    assert!(
        n_pre_pcr_ins_mismatch_reads > 0,
        "live genotyping reads must show BI vs Q45 GOP mismatch before PCR"
    );
}
