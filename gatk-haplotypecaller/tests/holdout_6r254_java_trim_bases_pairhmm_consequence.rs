//! 6R.254 live: Java 174-mers on frozen Rust 123-read evidence at
//! `20:29455649 T/TGTTTG` move 2/2 GL *away* from 3517.5 (diagnostic PL 3519).
//! Proof-only. Skipped unless `HOLDOUT_6R254=1`.
//!
//! ```text
//! HOLDOUT_6R254=1 cargo test -p gatk-haplotypecaller --test holdout_6r254_java_trim_bases_pairhmm_consequence -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::assembly_region_finalize::{
    clip_finalized_reads_to_region, finalize_region_reads_for_assembly,
    gatk_min_tail_quality_for_assembly,
};
use gatk_haplotypecaller::hc_genotyping_engine::{
    take_colocated_merge_numerics, HcGenotypingConfig,
};
use gatk_haplotypecaller::read_threading_assembler::DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH;
use gatk_haplotypecaller::{
    begin_hap_list_observe, biallelic_genotype_log10_likelihoods_gatk, call_disposition,
    flatten_assembly_regions, region_likelihoods_to_rows, score_read_against_haplotypes,
    take_hap_list_trim_span, traverse_assembly_region_walker, try_emit_call_region_variants,
    AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition, HaplotypeCallerEngine,
    HcLikelihoodEngineConfig, PairHmmBackend, ReadFilterParams, ReadLikelihoodRow,
    WalkerTraversalConfig, DEFAULT_STAND_CALL_CONF, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use rust_htslib::bam::record::Aux;
use std::collections::HashMap;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_VCF_REL: &str = "parity/reports/6r43/chr20_tiny/java.vcf";
const TARGET: u64 = 29_455_649;
const TARGET_REF: &str = "T";
const TARGET_ALT: &str = "TGTTTG";
const CLOSED_GC: u64 = 29_455_314;
const FROZEN_IDX: [usize; 5] = [19, 20, 21, 22, 23];
const RUST_FNV: [&str; 5] = [
    "55012fcf3b430591",
    "341b2e070ccb5846",
    "6a65e4c02733c2ed",
    "fa07750bf228b3c2",
    "79451c576721a729",
];
const JAVA_FNV: [&str; 5] = [
    "fcd72c6ce600dd16",
    "c1a4e8204522f645",
    "eb03271fa7548f26",
    "7dc2a8ae5da116e1",
    "343e7c443c1d9732",
];
const FLOOR: f64 = -4.5;
const JAVA_TSV: &str = include_str!("forensic_6r252_java_hap_trim.tsv");
const TRIMMER_SRC: &str = include_str!("../src/assembly_region_trimmer.rs");

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R254\t{key}\t{}", value.as_ref());
}

fn fnv1a64_hex(bases: &[u8]) -> String {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in bases {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    format!("{h:016x}")
}

fn java_seq(hash: &str) -> Vec<u8> {
    for line in JAVA_TSV.lines() {
        if line.contains("stage=trimmed") && line.contains(&format!("hash={hash}")) {
            if let Some(seq) = line.split("\tseq=").nth(1) {
                return seq.as_bytes().to_vec();
            }
        }
    }
    panic!("missing Java trimmed seq for {hash}");
}

fn bam_indel_phred(rec: &rust_htslib::bam::Record, tag: &[u8]) -> Option<Vec<u8>> {
    match rec.aux(tag) {
        Ok(Aux::String(s)) => Some(s.bytes().map(|b| b.saturating_sub(33)).collect()),
        _ => None,
    }
}

fn floor_pair(l0: f64, l2: f64) -> (f64, f64) {
    let best = l0.max(l2);
    if !best.is_finite() {
        return (l0, l2);
    }
    let floor = best + FLOOR;
    let lift = |v: f64| {
        if v.is_finite() && v < floor {
            floor
        } else {
            v
        }
    };
    (lift(l0), lift(l2))
}

fn biallelic_gls(l0: &[f64], l2: &[f64]) -> Vec<f64> {
    let rows: Vec<ReadLikelihoodRow> = l0
        .iter()
        .zip(l2.iter())
        .enumerate()
        .map(|(i, (a, b))| {
            let (fa, fb) = floor_pair(*a, *b);
            ReadLikelihoodRow {
                read_index: i,
                read_id: String::new(),
                haplotype_log10_likelihoods: vec![fa, fb],
            }
        })
        .collect();
    biallelic_genotype_log10_likelihoods_gatk(&rows, 0, 1)
}

fn emitted_hom_alt(gls: &[f64]) -> (f64, f64, i32) {
    let best = gls.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let rel = gls[2] - best;
    let cont = -10.0 * rel;
    let pl = (cont + 0.5).floor() as i32;
    (rel, cont, pl)
}

fn vcf_has(path: &Path, pos: u64, r: &str, a: &str) -> bool {
    for line in std::fs::read_to_string(path).unwrap_or_default().lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<_> = line.split('\t').collect();
        if f.len() < 5 || f[0] != "20" {
            continue;
        }
        if f[1].parse::<u64>().ok() == Some(pos) && f[3] == r && f[4] == a {
            return true;
        }
    }
    false
}

#[test]
fn holdout_6r254_java_trim_bases_pairhmm_consequence() {
    if std::env::var("HOLDOUT_6R254").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R254=1");
        return;
    }
    assert_eq!(DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, 128);
    assert_eq!(
        HcGenotypingConfig::strict_java().stand_emit_confidence,
        DEFAULT_STAND_CALL_CONF
    );
    let cfg = HcLikelihoodEngineConfig::gatk_haplotype_caller_production();
    assert_eq!(cfg.resolved_pair_hmm_backend(), PairHmmBackend::NeonF64);
    assert!(TRIMMER_SRC.contains("let slice = &bases[rel_start..=rel_end];"));
    assert!(TRIMMER_SRC.contains("if slice.len() < 2"));

    let root = repo_root();
    assert!(vcf_has(
        &root.join(JAVA_VCF_REL),
        TARGET,
        TARGET_REF,
        TARGET_ALT
    ));
    assert!(!vcf_has(&root.join(JAVA_VCF_REL), CLOSED_GC, "G", "C"));

    let java_haps: Vec<Vec<u8>> = JAVA_FNV.iter().map(|h| java_seq(h)).collect();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, INTERVAL).expect("interval");
    begin_hap_list_observe();
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
    let outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let trim = take_hap_list_trim_span().expect("trim");
    assert_eq!(trim.trim_start, 29_455_569);
    assert_eq!(trim.trim_end, 29_455_724);
    let snap = take_colocated_merge_numerics()
        .into_iter()
        .find(|s| s.loc == TARGET)
        .expect("snap");
    assert_eq!(snap.n_reads, 123);
    assert_eq!(snap.long_ref, TARGET_REF);
    assert_eq!(snap.assigned_gt, vec![0, 2]);
    assert_eq!(snap.alts[1], TARGET_ALT);

    let rust_haps: Vec<Vec<u8>> = FROZEN_IDX
        .iter()
        .map(|&idx| outcome.assembly.haplotypes[idx].bases.clone())
        .collect();
    for (i, (java, rust)) in java_haps.iter().zip(rust_haps.iter()).enumerate() {
        assert_eq!(fnv1a64_hex(java), JAVA_FNV[i]);
        assert_eq!(fnv1a64_hex(rust), RUST_FNV[i]);
        assert_eq!(java.len(), 174);
        assert_eq!(rust.len(), 161);
        assert_eq!(&java[9..java.len() - 4], rust.as_slice());
        assert_eq!(&java[..9], b"CAAAGAGTA");
        assert_eq!(&java[java.len() - 4..], b"TAAA");
    }

    let recs =
        try_emit_call_region_variants(region, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    let emitted = recs
        .iter()
        .find(|r| {
            r.position == TARGET
                && r.reference == TARGET_REF
                && r.alternate.iter().any(|a| a == TARGET_ALT)
        })
        .expect("emit");
    assert_eq!(
        emitted.samples[0].pl.as_deref(),
        Some(&[570, 0, 3518][..]),
        "production PL must remain 3518"
    );

    let retain: std::collections::HashSet<usize> = snap.ad_row_read_index.iter().copied().collect();
    let subset: Vec<_> = outcome
        .read_likelihoods
        .iter()
        .filter(|c| retain.contains(&c.read_index.get()))
        .cloned()
        .collect();
    let hap_rows = region_likelihoods_to_rows(&subset, outcome.assembly.haplotypes.len());
    let mut by_read = HashMap::new();
    for row in &hap_rows {
        by_read.entry(row.read_index).or_insert(row);
    }
    let l0: Vec<f64> = snap.ad_row_lls.iter().map(|ll| ll[0]).collect();
    let l2_prod: Vec<f64> = snap.ad_row_lls.iter().map(|ll| ll[2]).collect();
    let baseline = biallelic_gls(&l0, &l2_prod);
    let (base_rel, base_cont, base_pl) = emitted_hom_alt(&baseline);
    assert_eq!(base_pl, 3518);
    assert!(base_cont > 3517.5);

    let finalized = finalize_region_reads_for_assembly(
        &region.reads,
        region,
        true,
        gatk_min_tail_quality_for_assembly(10),
        false,
    );
    let mut clip_region = region.clone();
    clip_region.extended_start = GenomePosition::new_1based(trim.trim_start);
    clip_region.extended_end = GenomePosition::new_1based(trim.trim_end);
    let mut pairhmm_reads = clip_finalized_reads_to_region(&finalized, &clip_region);
    pairhmm_reads.retain(|r| r.seq().len() >= 10);
    let mut by_id: HashMap<(Vec<u8>, u16), &rust_htslib::bam::Record> = HashMap::new();
    for rec in &pairhmm_reads {
        by_id.insert((rec.qname().to_vec(), rec.flags()), rec);
    }
    let rust_refs: Vec<&[u8]> = rust_haps.iter().map(|s| s.as_slice()).collect();
    let java_refs: Vec<&[u8]> = java_haps.iter().map(|s| s.as_slice()).collect();
    let mut all_haps = rust_refs.clone();
    all_haps.extend_from_slice(&java_refs);
    let mut java_diag = vec![0.0f64; 123];
    let mut rust_diag = vec![0.0f64; 123];
    for (ri, &read_idx) in snap.ad_row_read_index.iter().enumerate() {
        let rec = *by_id
            .get(&(
                snap.ad_row_qname[ri].as_bytes().to_vec(),
                snap.ad_row_flags[ri],
            ))
            .expect("PairHMM-time read");
        let prod_row = by_read.get(&read_idx).expect("prod row");
        let p5: [f64; 5] =
            std::array::from_fn(|k| prod_row.haplotype_log10_likelihoods[FROZEN_IDX[k]]);
        if p5.iter().all(|&v| v == 0.0) {
            continue;
        }
        let scores = score_read_against_haplotypes(
            &cfg,
            &rec.seq().as_bytes(),
            rec.qual(),
            rec.mapq(),
            &all_haps,
            bam_indel_phred(rec, b"BI").as_deref(),
            bam_indel_phred(rec, b"BD").as_deref(),
        )
        .expect("score");
        let r5: [f64; 5] = std::array::from_fn(|k| scores[k]);
        let j5: [f64; 5] = std::array::from_fn(|k| scores[5 + k]);
        let other_best = prod_row
            .haplotype_log10_likelihoods
            .iter()
            .enumerate()
            .filter(|(i, _)| !FROZEN_IDX.contains(i))
            .map(|(_, v)| *v)
            .fold(f64::NEG_INFINITY, f64::max);
        let floor5 = |raw: [f64; 5]| -> [f64; 5] {
            let raw_max = raw.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            let best = other_best.max(raw_max);
            let floor = best + FLOOR;
            std::array::from_fn(|k| {
                if raw[k].is_finite() && raw[k] < floor {
                    floor
                } else {
                    raw[k]
                }
            })
        };
        rust_diag[ri] = floor5(r5).iter().copied().fold(f64::NEG_INFINITY, f64::max);
        java_diag[ri] = floor5(j5).iter().copied().fold(f64::NEG_INFINITY, f64::max);
    }
    let rust_gls = biallelic_gls(&l0, &rust_diag);
    let (_, _, r_pl) = emitted_hom_alt(&rust_gls);
    assert_eq!(r_pl, 3518);
    let java_gls = biallelic_gls(&l0, &java_diag);
    let (j_rel, j_cont, j_pl) = emitted_hom_alt(&java_gls);
    assert_eq!(j_pl, 3519);
    assert!(j_cont > 3517.5);
    assert!(j_rel < base_rel);
    kv("production_integer_PL", "3518");
    kv("java174_integer_PL", "3519");
    kv("crosses_3517_5", "false");
    kv("classification", "TRIM_SPAN_MOVES_GL_BUT_NOT_PL_BOUNDARY");
    kv("production_change", "NONE");
}
