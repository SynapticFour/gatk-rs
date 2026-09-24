//! 6R.255: measurement-only. Does the Java hard-clip interval
//! `20:29455560-29455728` vs Rust `20:29455569-29455724` move 2/2 GL
//! across 3517.5 at `20:29455649 T/TGTTTG`?
//!
//! Frozen 123-read retainEvidence and Rust 161-mer TGTTTG haplotypes.
//! Only the read clip span changes. PRODUCTION CHANGE: NONE.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r255_read_clip_interval_pairhmm_consequence -- --nocapture --test-threads=1
//! HOLDOUT_6R255=1 cargo test -p gatk-haplotypecaller --test holdout_6r255_read_clip_interval_pairhmm_consequence -- --nocapture --test-threads=1
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
use std::collections::HashMap;
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
const FLOOR: f64 = -4.5;
const RUST_EMITTED_HOM_ALT_GL: f64 = -351.75162674083281900;
const RUST_CLIP: (u64, u64) = (29_455_569, 29_455_724);
const JAVA_CLIP: (u64, u64) = (29_455_560, 29_455_728);
const C1_LEFT: (u64, u64) = (29_455_560, 29_455_724);
const C2_RIGHT: (u64, u64) = (29_455_569, 29_455_728);

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R255\t{key}\t{}", value.as_ref());
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

fn floor5(raw: [f64; 5], other_best: f64) -> [f64; 5] {
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
}

fn winner_mask(five: &[f64; 5], mx: f64) -> u8 {
    let mut m = 0u8;
    for (k, &v) in five.iter().enumerate() {
        if v == mx {
            m |= 1 << k;
        }
    }
    m
}

fn clip_map(
    finalized: &[rust_htslib::bam::Record],
    region: &gatk_haplotypecaller::assembly_region_iterator::AssemblyRegion,
    start: u64,
    end: u64,
) -> HashMap<(Vec<u8>, u16), rust_htslib::bam::Record> {
    let mut clip_region = region.clone();
    clip_region.extended_start = GenomePosition::new_1based(start);
    clip_region.extended_end = GenomePosition::new_1based(end);
    let mut clipped = clip_finalized_reads_to_region(finalized, &clip_region);
    clipped.retain(|r| r.seq().len() >= 10);
    let mut by_id = HashMap::new();
    for rec in clipped {
        by_id.insert((rec.qname().to_vec(), rec.flags()), rec);
    }
    by_id
}

#[allow(dead_code)]
struct ExpOut {
    pooled: Vec<f64>,
    raw_max: Vec<f64>,
    winners: Vec<u8>,
    n_bases: usize,
    n_len_changed: usize,
    n_missing: usize,
    n_above_floor: usize,
}

#[test]
fn forensic_6r255_read_clip_interval_pairhmm_consequence() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    kv("target", "20:29455649 T/TGTTTG");
    kv(
        "experiment",
        "diagnostic PairHMM: frozen 123 reads × Rust 161-mers; only clip interval varies",
    );
    kv(
        "java_clip_site",
        "GATK 4.4 AssemblyBasedCallerUtils.finalizeRegion → ReadClipper.hardClipToRegion(padded variant span); PairHMM uses those clipped GATKRead bases unchanged",
    );
    kv("rust_clip_interval", "20:29455569-29455724");
    kv("java_clip_interval", "20:29455560-29455728");
    assert_eq!(DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, 128);

    let root = repo_root();
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
    let trim = take_hap_list_trim_span().expect("trim span");
    assert_eq!((trim.trim_start, trim.trim_end), RUST_CLIP);
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
    for (i, seq) in rust_haps.iter().enumerate() {
        assert_eq!(fnv1a64_hex(seq), RUST_FNV[i]);
        assert_eq!(seq.len(), 161);
    }
    kv(
        "haplotypes",
        "Rust 161-mers idx=19..23 frozen; Java 174-mers NOT used",
    );

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

    let retain: std::collections::HashSet<usize> = snap.ad_row_read_index.iter().copied().collect();
    assert_eq!(retain.len(), 123);
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
    kv("snap_GL_00", fmt_f64(baseline[0]));
    kv("snap_GL_02", fmt_f64(baseline[1]));
    kv("snap_GL_22", fmt_f64(baseline[2]));
    kv("snap_emitted_hom_alt_GL", fmt_f64(base_rel));
    kv("snap_continuous_PL", fmt_f64(base_cont));
    kv("snap_integer_PL", base_pl.to_string());
    assert!((base_rel - RUST_EMITTED_HOM_ALT_GL).abs() < 1e-9);
    assert!(base_cont > 3517.5);
    assert_eq!(base_pl, 3518);

    let finalized = finalize_region_reads_for_assembly(
        &region.reads,
        region,
        true,
        gatk_min_tail_quality_for_assembly(10),
        false,
    );
    kv("n_finalized_preclip", finalized.len().to_string());
    kv(
        "java_hardclip_contract",
        "ReadClipper.hardClipToRegion on regionForGenotyping padded span; Rust clip_finalized_reads_in_place/hard_clip_to_region is the same op",
    );

    let hap_refs: Vec<&[u8]> = rust_haps.iter().map(|s| s.as_slice()).collect();
    let experiments = [
        ("A_rust", RUST_CLIP),
        ("B_java", JAVA_CLIP),
        ("C1_left", C1_LEFT),
        ("C2_right", C2_RIGHT),
    ];
    let mut maps = Vec::new();
    for (name, (start, end)) in experiments {
        maps.push((name, start, end, clip_map(&finalized, region, start, end)));
    }

    let mut outs: Vec<ExpOut> = Vec::new();
    for (name, start, end, map) in &maps {
        let mut pooled = vec![0.0f64; 123];
        let mut raw_max = vec![0.0f64; 123];
        let mut winners = vec![0u8; 123];
        let mut n_bases = 0usize;
        let mut n_len_changed = 0usize;
        let mut n_missing = 0usize;
        let mut n_above_floor = 0usize;
        let rust_map = &maps[0].3;
        for (ri, &read_idx) in snap.ad_row_read_index.iter().enumerate() {
            let key = (
                snap.ad_row_qname[ri].as_bytes().to_vec(),
                snap.ad_row_flags[ri],
            );
            let Some(rec) = map.get(&key) else {
                n_missing += 1;
                kv(
                    &format!("{name}_missing_row_{ri}"),
                    format!(
                        "qname={} flags={}",
                        snap.ad_row_qname[ri], snap.ad_row_flags[ri]
                    ),
                );
                continue;
            };
            n_bases += rec.seq().len();
            if let Some(r0) = rust_map.get(&key) {
                if rec.seq().len() != r0.seq().len() {
                    n_len_changed += 1;
                }
            }
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
                &hap_refs,
                bam_indel_phred(rec, b"BI").as_deref(),
                bam_indel_phred(rec, b"BD").as_deref(),
            )
            .expect("score");
            let r5: [f64; 5] = std::array::from_fn(|k| scores[k]);
            let other_best = prod_row
                .haplotype_log10_likelihoods
                .iter()
                .enumerate()
                .filter(|(i, _)| !FROZEN_IDX.contains(i))
                .map(|(_, v)| *v)
                .fold(f64::NEG_INFINITY, f64::max);
            let r5f = floor5(r5, other_best);
            let rmax = r5f.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            let raw = r5.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            pooled[ri] = rmax;
            raw_max[ri] = raw;
            winners[ri] = winner_mask(&r5f, rmax);
            if (raw - rmax).abs() > 1e-15 {
                // floored
            } else {
                n_above_floor += 1;
            }
            if *name == "A_rust" {
                let pmax = p5.iter().copied().fold(f64::NEG_INFINITY, f64::max);
                assert!(
                    (rmax - pmax).abs() < 1e-9,
                    "Experiment A must match production five-column max row={ri} {rmax} vs {pmax}"
                );
            }
        }
        kv(
            &format!("{name}_clip"),
            format!(
                "span=20:{start}-{end} n_mapped={} n_missing={n_missing} n_bases={n_bases} n_len_changed_vs_A={n_len_changed} n_above_hap_floor={n_above_floor}",
                map.len()
            ),
        );
        assert_eq!(
            n_missing, 0,
            "{name}: frozen 123 QNAME+FLAG must survive clip 20:{start}-{end}"
        );
        outs.push(ExpOut {
            pooled,
            raw_max,
            winners,
            n_bases,
            n_len_changed,
            n_missing,
            n_above_floor,
        });
    }

    let mut gls_pl = Vec::new();
    for (i, (name, _, _, _)) in maps.iter().enumerate() {
        let gls = biallelic_gls(&l0, &outs[i].pooled);
        let (rel, cont, pl) = emitted_hom_alt(&gls);
        kv(&format!("{name}_GL_00"), fmt_f64(gls[0]));
        kv(&format!("{name}_GL_02"), fmt_f64(gls[1]));
        kv(&format!("{name}_GL_22"), fmt_f64(gls[2]));
        kv(&format!("{name}_emitted_hom_alt_GL"), fmt_f64(rel));
        kv(&format!("{name}_continuous_PL"), fmt_f64(cont));
        kv(&format!("{name}_integer_PL"), pl.to_string());
        kv(
            &format!("{name}_total_read_bases"),
            outs[i].n_bases.to_string(),
        );
        gls_pl.push((rel, cont, pl, gls));
    }

    let (a_rel, a_cont, a_pl, _) = &gls_pl[0];
    assert!(
        (a_rel - base_rel).abs() < 1e-9,
        "Experiment A must reproduce 6R.254 baseline GL"
    );
    assert_eq!(*a_pl, 3518);
    assert!((*a_cont - base_cont).abs() < 1e-6);

    let a = &outs[0];
    let b = &outs[1];
    let mut deltas: Vec<(usize, f64)> = Vec::new();
    let mut n_ident = 0usize;
    let mut n_winner_change = 0usize;
    let mut n_floor_cross = 0usize;
    for ri in 0..123 {
        let d = b.pooled[ri] - a.pooled[ri];
        if d == 0.0 {
            n_ident += 1;
        } else {
            deltas.push((ri, d));
        }
        if a.winners[ri] != b.winners[ri] {
            n_winner_change += 1;
        }
        let a_floored = (a.raw_max[ri] - a.pooled[ri]).abs() > 1e-15 && a.pooled[ri] != 0.0;
        let b_floored = (b.raw_max[ri] - b.pooled[ri]).abs() > 1e-15 && b.pooled[ri] != 0.0;
        if a_floored != b_floored {
            n_floor_cross += 1;
        }
    }
    deltas.sort_by(|x, y| x.1.abs().partial_cmp(&y.1.abs()).unwrap());
    let n_changed = deltas.len();
    let dvals: Vec<f64> = (0..123).map(|i| b.pooled[i] - a.pooled[i]).collect();
    let mut sorted = dvals.clone();
    sorted.sort_by(|x, y| x.partial_cmp(y).unwrap());
    let dmin = sorted[0];
    let dmax = sorted[122];
    let dmean = dvals.iter().sum::<f64>() / 123.0;
    let dmed = sorted[61];
    kv("A_vs_B_n_identical_pooled", n_ident.to_string());
    kv("A_vs_B_n_changed_pooled", n_changed.to_string());
    kv("A_vs_B_delta_min", fmt_f64(dmin));
    kv("A_vs_B_delta_max", fmt_f64(dmax));
    kv("A_vs_B_delta_mean", fmt_f64(dmean));
    kv("A_vs_B_delta_median", fmt_f64(dmed));
    kv("A_vs_B_n_winner_mask_changed", n_winner_change.to_string());
    kv("A_vs_B_n_hap_floor_cross", n_floor_cross.to_string());
    kv("A_vs_B_n_len_changed", outs[1].n_len_changed.to_string());
    for &(ri, d) in deltas.iter().rev().take(16) {
        kv(
            &format!("A_vs_B_shifted_row_{ri}"),
            format!(
                "qname={} flags={} delta={} A_pooled={} B_pooled={} A_raw={} B_raw={} A_win={:05b} B_win={:05b} A_len={} B_len={}",
                snap.ad_row_qname[ri],
                snap.ad_row_flags[ri],
                fmt_f64(d),
                fmt_f64(a.pooled[ri]),
                fmt_f64(b.pooled[ri]),
                fmt_f64(a.raw_max[ri]),
                fmt_f64(b.raw_max[ri]),
                a.winners[ri],
                b.winners[ri],
                maps[0].3.get(&(snap.ad_row_qname[ri].as_bytes().to_vec(), snap.ad_row_flags[ri])).map(|r| r.seq().len()).unwrap_or(0),
                maps[1].3.get(&(snap.ad_row_qname[ri].as_bytes().to_vec(), snap.ad_row_flags[ri])).map(|r| r.seq().len()).unwrap_or(0),
            ),
        );
    }

    let (b_rel, b_cont, b_pl, _) = &gls_pl[1];
    assert_eq!(*b_pl, 4280);
    assert!(*b_cont > 3517.5);
    assert!(*b_rel < *a_rel);
    kv("displacement_emitted_GL", fmt_f64(b_rel - a_rel));
    kv("displacement_continuous_PL", fmt_f64(b_cont - a_cont));
    kv("crosses_3517_5", (*b_cont < 3517.5).to_string());
    kv(
        "B_total_read_bases_minus_A",
        (outs[1].n_bases as i64 - outs[0].n_bases as i64).to_string(),
    );

    for (label, idx) in [("C1_minus_A", 2usize), ("C2_minus_A", 3)] {
        let d = gls_pl[idx].1 - a_cont;
        kv(
            label,
            format!(
                "emitted_GL={} continuous_PL={} integer_PL={} dPL={}",
                fmt_f64(gls_pl[idx].0),
                fmt_f64(gls_pl[idx].1),
                gls_pl[idx].2,
                fmt_f64(d)
            ),
        );
    }

    let classification = if *b_cont < 3517.5 {
        "READ_CLIP_INTERVAL_CAUSALLY_MOVES_PL_BOUNDARY"
    } else if (*b_rel - a_rel).abs() > 1e-6 {
        "READ_CLIP_INTERVAL_MOVES_GL_BUT_NOT_PL_BOUNDARY"
    } else {
        "READ_CLIP_INTERVAL_NOT_CAUSAL"
    };
    kv("classification", classification);
    kv(
        "next_arrow",
        if classification == "READ_CLIP_INTERVAL_CAUSALLY_MOVES_PL_BOUNDARY" {
            "do not patch this round; next would prove whether Java STR-aware padding (and its clip span) is the correct production change"
        } else if classification == "READ_CLIP_INTERVAL_MOVES_GL_BUT_NOT_PL_BOUNDARY" {
            "clip interval is consequential but not the PL ±1 cause; remaining arrow is not STR padding by itself"
        } else {
            "clip interval is not causal; remaining first input divergence is elsewhere"
        },
    );
    assert!(
        classification == "READ_CLIP_INTERVAL_CAUSALLY_MOVES_PL_BOUNDARY"
            || classification == "READ_CLIP_INTERVAL_MOVES_GL_BUT_NOT_PL_BOUNDARY"
            || classification == "READ_CLIP_INTERVAL_NOT_CAUSAL"
    );
}
