//! 6R.254: measurement-only. Do Java 174-mer TGTTTG haplotypes move Rust
//! 2/2 GL across the 3517.5 PL boundary at `20:29455649 T/TGTTTG`?
//!
//! Frozen 123-read retainEvidence. Diagnostic PairHMM only — production
//! haplotypes/trim/STR are unchanged.
//! PRODUCTION CHANGE: NONE.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r254_java_trim_bases_pairhmm_consequence -- --nocapture --test-threads=1
//! HOLDOUT_6R254=1 cargo test -p gatk-haplotypecaller --test holdout_6r254_java_trim_bases_pairhmm_consequence -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::assembly_region_finalize::{
    clip_finalized_reads_to_region, finalize_region_reads_for_assembly,
    gatk_min_tail_quality_for_assembly,
};
use gatk_haplotypecaller::hc_genotyping_engine::take_colocated_merge_numerics;
use gatk_haplotypecaller::read_threading_assembler::DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH;
use gatk_haplotypecaller::{
    begin_hap_list_observe, biallelic_genotype_log10_likelihoods_gatk, call_disposition,
    flatten_assembly_regions, region_likelihoods_to_rows, score_read_against_haplotypes,
    take_hap_list_trim_span, traverse_assembly_region_walker, AssemblyRegionCallDisposition,
    CallRegionArgs, GenomePosition, HaplotypeCallerEngine, HcLikelihoodEngineConfig,
    PairHmmBackend, ReadFilterParams, ReadLikelihoodRow, WalkerTraversalConfig,
};
use rust_htslib::bam::record::Aux;
use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 29_455_649;
const TARGET_REF: &str = "T";
const TARGET_ALT: &str = "TGTTTG";
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
const JAVA_PREFIX: &[u8] = b"CAAAGAGTA";
const JAVA_SUFFIX: &[u8] = b"TAAA";
const FLOOR: f64 = -4.5;
const RUST_EMITTED_HOM_ALT_GL: f64 = -351.75162674083281900;
const JAVA_TSV: &str = include_str!("forensic_6r252_java_hap_trim.tsv");

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R254\t{key}\t{}", value.as_ref());
}

fn fmt_f64(x: f64) -> String {
    format!("{x:.17} bits=0x{:016x}", x.to_bits())
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

#[test]
fn forensic_6r254_java_trim_bases_pairhmm_consequence() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    kv("target", "20:29455649 T/TGTTTG");
    kv(
        "experiment",
        "diagnostic PairHMM: frozen 123 reads × Java 174-mers",
    );
    assert_eq!(DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, 128);
    assert_eq!(JAVA_PREFIX.len(), 9);
    assert_eq!(JAVA_SUFFIX.len(), 4);

    let java_haps: Vec<Vec<u8>> = JAVA_FNV.iter().map(|h| java_seq(h)).collect();
    for (i, seq) in java_haps.iter().enumerate() {
        assert_eq!(fnv1a64_hex(seq), JAVA_FNV[i]);
        assert_eq!(seq.len(), 174);
        assert_eq!(&seq[..9], JAVA_PREFIX);
        assert_eq!(&seq[seq.len() - 4..], JAVA_SUFFIX);
        kv(
            &format!("java_H{i}"),
            format!(
                "fnv={} len=174 cigar=90M5I79M prefix={} suffix={}",
                JAVA_FNV[i],
                std::str::from_utf8(JAVA_PREFIX).unwrap(),
                std::str::from_utf8(JAVA_SUFFIX).unwrap()
            ),
        );
        kv(
            &format!("java_H{i}_bases"),
            String::from_utf8_lossy(seq).into_owned(),
        );
    }

    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing chr20_tiny BAM/REF");
        return;
    }
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
    let trim = take_hap_list_trim_span().expect("trim span");
    assert_eq!(trim.trim_start, 29_455_569);
    assert_eq!(trim.trim_end, 29_455_724);
    kv(
        "pairhmm_clip_span",
        format!("20:{}-{}", trim.trim_start, trim.trim_end),
    );
    let snap = take_colocated_merge_numerics()
        .into_iter()
        .find(|s| s.loc == TARGET)
        .expect("snap");
    assert_eq!(snap.n_reads, 123);
    assert_eq!(snap.long_ref, TARGET_REF);
    assert_eq!(snap.assigned_gt, vec![0, 2]);
    assert_eq!(snap.alts[1], TARGET_ALT);
    kv("n_reads_retainEvidence", "123");

    let rust_haps: Vec<Vec<u8>> = FROZEN_IDX
        .iter()
        .map(|&idx| outcome.assembly.haplotypes[idx].bases.clone())
        .collect();
    let mut stripped: Vec<Vec<u8>> = Vec::new();
    for (i, (java, rust)) in java_haps.iter().zip(rust_haps.iter()).enumerate() {
        assert_eq!(fnv1a64_hex(rust), RUST_FNV[i]);
        assert_eq!(rust.len(), 161);
        kv(
            &format!("rust_H{i}_bases"),
            String::from_utf8_lossy(rust).into_owned(),
        );
        let interior = &java[9..java.len() - 4];
        assert_eq!(
            interior,
            rust.as_slice(),
            "Java 174-mer interior must be the Rust 161-mer for H{i}"
        );
        let left_aligned_mismatch = java
            .iter()
            .zip(rust.iter())
            .position(|(a, b)| a != b)
            .unwrap_or(rust.len());
        kv(
            &format!("pair_H{i}"),
            format!(
                "rust_fnv={} java_fnv={} rust_len=161 java_len=174 prefix={} suffix={} left_aligned_first_mismatch={} structural=Java[9:-4]==Rust interior_match=YES",
                RUST_FNV[i],
                JAVA_FNV[i],
                std::str::from_utf8(JAVA_PREFIX).unwrap(),
                std::str::from_utf8(JAVA_SUFFIX).unwrap(),
                left_aligned_mismatch
            ),
        );
        stripped.push(interior.to_vec());
        assert_eq!(stripped[i], *rust);
    }

    let cfg = HcLikelihoodEngineConfig::gatk_haplotype_caller_production();
    assert_eq!(cfg.resolved_pair_hmm_backend(), PairHmmBackend::NeonF64);
    kv(
        "pairhmm_backend",
        format!(
            "configured={} resolved={}",
            cfg.primary_engine_label(),
            cfg.resolved_pair_hmm_backend().label()
        ),
    );

    let retain: HashSet<usize> = snap.ad_row_read_index.iter().copied().collect();
    assert_eq!(retain.len(), 123);
    let subset: Vec<_> = outcome
        .read_likelihoods
        .iter()
        .filter(|c| retain.contains(&c.read_index.get()))
        .cloned()
        .collect();
    let hap_rows = region_likelihoods_to_rows(&subset, outcome.assembly.haplotypes.len());
    let mut by_read = std::collections::HashMap::new();
    for row in &hap_rows {
        by_read.entry(row.read_index).or_insert(row);
    }

    let mut l0 = Vec::with_capacity(123);
    let mut l2_prod = Vec::with_capacity(123);
    for ll in &snap.ad_row_lls {
        l0.push(ll[0]);
        l2_prod.push(ll[2]);
    }
    assert_eq!(l0.len(), 123);

    let baseline = biallelic_gls(&l0, &l2_prod);
    let (base_rel, base_cont, base_pl) = emitted_hom_alt(&baseline);
    kv("baseline_GL_00", fmt_f64(baseline[0]));
    kv("baseline_GL_02", fmt_f64(baseline[1]));
    kv("baseline_GL_22", fmt_f64(baseline[2]));
    kv("baseline_emitted_hom_alt_GL", fmt_f64(base_rel));
    kv("baseline_continuous_PL", fmt_f64(base_cont));
    kv("baseline_integer_PL", base_pl.to_string());
    assert!((base_rel - RUST_EMITTED_HOM_ALT_GL).abs() < 1e-9);
    assert!(base_cont > 3517.5);
    assert_eq!(base_pl, 3518);

    let rust_refs: Vec<&[u8]> = rust_haps.iter().map(|s| s.as_slice()).collect();
    let java_refs: Vec<&[u8]> = java_haps.iter().map(|s| s.as_slice()).collect();
    let strip_refs: Vec<&[u8]> = stripped.iter().map(|s| s.as_slice()).collect();
    let mut all_haps: Vec<&[u8]> = Vec::with_capacity(15);
    all_haps.extend_from_slice(&rust_refs);
    all_haps.extend_from_slice(&java_refs);
    all_haps.extend_from_slice(&strip_refs);

    let mut rust_diag = vec![0.0f64; 123];
    let mut java_diag = vec![0.0f64; 123];
    let mut strip_diag = vec![0.0f64; 123];
    let mut n_winner_change = 0usize;
    let mut n_h3_still = 0usize;
    let mut n_unique_prod = 0usize;
    let mut n_unique_java = 0usize;
    let mut n_unique_prod_h3 = 0usize;
    let mut n_unique_java_h3 = 0usize;
    let mut unique_winner_change = 0usize;
    let mut col_d_min = [f64::INFINITY; 5];
    let mut col_d_max = [f64::NEG_INFINITY; 5];
    let mut col_d_sum = [0.0f64; 5];
    let mut pooled_d_min = f64::INFINITY;
    let mut pooled_d_max = f64::NEG_INFINITY;
    let mut pooled_d_sum = 0.0f64;
    let mut n_control_exact = 0usize;

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
    kv(
        "pairhmm_time_reads",
        format!(
            "finalized={} clipped={} mapped_ids={}",
            finalized.len(),
            pairhmm_reads.len(),
            by_id.len()
        ),
    );

    for (ri, &read_idx) in snap.ad_row_read_index.iter().enumerate() {
        let qname = snap.ad_row_qname[ri].as_bytes();
        let flags = snap.ad_row_flags[ri];
        let rec = *by_id
            .get(&(qname.to_vec(), flags))
            .unwrap_or_else(|| panic!("missing PairHMM-time read row={ri} qname flags"));
        let bi = bam_indel_phred(rec, b"BI");
        let bd = bam_indel_phred(rec, b"BD");
        let scores = score_read_against_haplotypes(
            &cfg,
            &rec.seq().as_bytes(),
            rec.qual(),
            rec.mapq(),
            &all_haps,
            bi.as_deref(),
            bd.as_deref(),
        )
        .expect("score");
        assert_eq!(scores.len(), 15);
        let r5: [f64; 5] = std::array::from_fn(|k| scores[k]);
        let j5: [f64; 5] = std::array::from_fn(|k| scores[5 + k]);
        let s5: [f64; 5] = std::array::from_fn(|k| scores[10 + k]);
        for k in 0..5 {
            assert_eq!(
                s5[k].to_bits(),
                r5[k].to_bits(),
                "control strip H{k} row {ri} must match Rust 161-mer LL"
            );
            let d = j5[k] - r5[k];
            col_d_min[k] = col_d_min[k].min(d);
            col_d_max[k] = col_d_max[k].max(d);
            col_d_sum[k] += d;
        }
        n_control_exact += 1;
        let prod_row = by_read.get(&read_idx).expect("prod row");
        let p5: [f64; 5] =
            std::array::from_fn(|k| prod_row.haplotype_log10_likelihoods[FROZEN_IDX[k]]);
        let pmax = p5.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        // Java `addEvidence(..., 0)` mate-contig / skipped PairHMM cells. Do not
        // invent a kernel score for those rows; haplotype bases never entered PairHMM.
        if p5.iter().all(|&v| v == 0.0) {
            rust_diag[ri] = 0.0;
            java_diag[ri] = 0.0;
            strip_diag[ri] = 0.0;
            continue;
        }
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
        let r5f = floor5(r5);
        let j5f = floor5(j5);
        let s5f = floor5(s5);
        let rmax = r5f.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let jmax = j5f.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        rust_diag[ri] = rmax;
        java_diag[ri] = jmax;
        strip_diag[ri] = s5f.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let pd = jmax - rmax;
        pooled_d_min = pooled_d_min.min(pd);
        pooled_d_max = pooled_d_max.max(pd);
        pooled_d_sum += pd;
        assert!(
            (rmax - pmax).abs() < 1e-9,
            "diagnostic hap-floored Rust 161-mer max must match production five-column max row={ri} {rmax} vs {pmax}"
        );

        let pw: Vec<usize> = p5
            .iter()
            .enumerate()
            .filter(|(_, v)| **v == pmax)
            .map(|(k, _)| k)
            .collect();
        let jw: Vec<usize> = j5f
            .iter()
            .enumerate()
            .filter(|(_, v)| **v == jmax)
            .map(|(k, _)| k)
            .collect();
        if pw.len() == 1 {
            n_unique_prod += 1;
            if pw[0] == 3 {
                n_unique_prod_h3 += 1;
            }
            kv(
                &format!("unique_max_row_{ri}"),
                format!(
                    "qname={} flags={} rust_winner=H{} java_winner=H{} floored_pooled_delta={} rust_floored_max={} java_floored_max={} rust_raw_max={} java_raw_max={}",
                    snap.ad_row_qname[ri],
                    flags,
                    pw[0],
                    jw.first().copied().unwrap_or(99),
                    fmt_f64(pd),
                    fmt_f64(rmax),
                    fmt_f64(jmax),
                    fmt_f64(r5.iter().copied().fold(f64::NEG_INFINITY, f64::max)),
                    fmt_f64(j5.iter().copied().fold(f64::NEG_INFINITY, f64::max))
                ),
            );
        }
        if jw.len() == 1 {
            n_unique_java += 1;
            if jw[0] == 3 {
                n_unique_java_h3 += 1;
            }
        }
        let p_set: HashSet<usize> = pw.iter().copied().collect();
        let j_set: HashSet<usize> = jw.iter().copied().collect();
        if p_set != j_set {
            n_winner_change += 1;
            if pw.len() == 1 && jw.len() == 1 {
                unique_winner_change += 1;
            }
        }
        if j_set.contains(&3) {
            n_h3_still += 1;
        }
    }
    assert_eq!(n_control_exact, 123);
    kv(
        "control_strip_reproduces_rust_161",
        "YES 123/123 exact bits",
    );
    kv("n_unique_max_production", n_unique_prod.to_string());
    kv("n_unique_max_java174", n_unique_java.to_string());
    kv("n_unique_max_prod_H3", n_unique_prod_h3.to_string());
    kv("n_unique_max_java_H3", n_unique_java_h3.to_string());
    kv("n_winner_set_changed", n_winner_change.to_string());
    kv("n_unique_winner_changed", unique_winner_change.to_string());
    kv("n_rows_H3_in_java_winner_set", n_h3_still.to_string());
    kv(
        "pooled_delta_java_minus_rust",
        format!(
            "min={} max={} mean={}",
            fmt_f64(pooled_d_min),
            fmt_f64(pooled_d_max),
            fmt_f64(pooled_d_sum / 123.0)
        ),
    );
    for k in 0..5 {
        kv(
            &format!("col_H{k}_delta"),
            format!(
                "min={} max={} mean={}",
                fmt_f64(col_d_min[k]),
                fmt_f64(col_d_max[k]),
                fmt_f64(col_d_sum[k] / 123.0)
            ),
        );
    }

    let strip_gls = biallelic_gls(&l0, &strip_diag);
    let (s_rel, s_cont, s_pl) = emitted_hom_alt(&strip_gls);
    kv("control_emitted_hom_alt_GL", fmt_f64(s_rel));
    kv("control_continuous_PL", fmt_f64(s_cont));
    kv("control_integer_PL", s_pl.to_string());
    assert!(
        (s_rel - base_rel).abs() < 1e-9,
        "stripped Java 174-mers must reproduce baseline GL"
    );
    assert_eq!(s_pl, 3518);

    let rust_rescore_gls = biallelic_gls(&l0, &rust_diag);
    let (_, r_cont, r_pl) = emitted_hom_alt(&rust_rescore_gls);
    assert_eq!(r_pl, 3518);
    assert!((r_cont - base_cont).abs() < 1e-6);

    let java_gls = biallelic_gls(&l0, &java_diag);
    let (j_rel, j_cont, j_pl) = emitted_hom_alt(&java_gls);
    kv("java174_GL_00", fmt_f64(java_gls[0]));
    kv("java174_GL_02", fmt_f64(java_gls[1]));
    kv("java174_GL_22", fmt_f64(java_gls[2]));
    kv("java174_emitted_hom_alt_GL", fmt_f64(j_rel));
    kv("java174_continuous_PL", fmt_f64(j_cont));
    kv("java174_integer_PL", j_pl.to_string());
    kv("displacement_emitted_GL", fmt_f64(j_rel - base_rel));
    kv("displacement_continuous_PL", fmt_f64(j_cont - base_cont));
    kv("crosses_3517_5", (j_cont <= 3517.5).to_string());
    kv(
        "java_pairhmm_oracle",
        "NOT collected: GKL native is x86_64; dump uses assembly-only on this host",
    );

    // A: Java 174-mers cross 3517.5 toward Java PL 3517.
    // B: they move GL but remain above 3517.5 (this round: they move *away*, to 3519).
    // C: they leave the 3518 boundary essentially unchanged.
    let classification = if j_cont <= 3517.5 {
        "TRIM_SPAN_CAUSALLY_EXPLAINS_PL_BOUNDARY"
    } else if (j_rel - base_rel).abs() > 1e-6 {
        "TRIM_SPAN_MOVES_GL_BUT_NOT_PL_BOUNDARY"
    } else {
        "TRIM_SPAN_INSUFFICIENT_FOR_PL"
    };
    assert_eq!(j_pl, 3519);
    assert!(j_cont > 3517.5);
    assert!(
        j_rel < base_rel,
        "Java 174-mers make 2/2 more negative, not less"
    );
    kv("classification", classification);
    kv(
        "next_arrow",
        if classification == "TRIM_SPAN_CAUSALLY_EXPLAINS_PL_BOUNDARY" {
            "prove whether Java STR-aware padding is the correct production change (do not patch this round)"
        } else if classification == "TRIM_SPAN_MOVES_GL_BUT_NOT_PL_BOUNDARY" {
            "measure remaining PL displacement; do not assume GKL residual. Next isolated input: Java-span hardClipToRegion on the same 123 reads (not STR padding)"
        } else {
            "trim difference remains real but another downstream difference is required"
        },
    );
    assert!(
        classification == "TRIM_SPAN_CAUSALLY_EXPLAINS_PL_BOUNDARY"
            || classification == "TRIM_SPAN_MOVES_GL_BUT_NOT_PL_BOUNDARY"
            || classification == "TRIM_SPAN_INSUFFICIENT_FOR_PL"
            || classification == "TRIM_SPAN_MEASURED_OTHER"
    );
}
