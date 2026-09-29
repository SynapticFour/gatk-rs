//! 6R.256: measurement-only. Joint Java 174-mer haplotypes + Java clip
//! interval vs frozen 123-read retainEvidence at `20:29455649 T/TGTTTG`.
//! PRODUCTION CHANGE: NONE.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r256_joint_java_input_plane -- --nocapture --test-threads=1
//! HOLDOUT_6R256=1 cargo test -p gatk-haplotypecaller --test holdout_6r256_joint_java_input_plane -- --nocapture --test-threads=1
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
const CTRL_RR_GL: f64 = -351.75162674083281900;
const CTRL_JR_GL: f64 = -351.89646448035699700;
const CTRL_RJ_GL: f64 = -428.04642955956700234;
const RUST_CLIP: (u64, u64) = (29_455_569, 29_455_724);
const JAVA_CLIP: (u64, u64) = (29_455_560, 29_455_728);
const UNIQUE_MAX_254: [usize; 12] = [10, 34, 35, 65, 67, 69, 73, 83, 88, 113, 120, 121];
const JAVA_TSV: &str = include_str!("forensic_6r252_java_hap_trim.tsv");

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R256\t{key}\t{}", value.as_ref());
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

struct RowScore {
    pooled: f64,
    raw_max: f64,
    winner: u8,
    floored: bool,
}

fn score_row(
    rec: &rust_htslib::bam::Record,
    haps: &[&[u8]],
    cfg: &HcLikelihoodEngineConfig,
    other_best: f64,
    skip: bool,
) -> RowScore {
    if skip {
        return RowScore {
            pooled: 0.0,
            raw_max: 0.0,
            winner: 0,
            floored: false,
        };
    }
    let scores = score_read_against_haplotypes(
        cfg,
        &rec.seq().as_bytes(),
        rec.qual(),
        rec.mapq(),
        haps,
        bam_indel_phred(rec, b"BI").as_deref(),
        bam_indel_phred(rec, b"BD").as_deref(),
    )
    .expect("score");
    let raw: [f64; 5] = std::array::from_fn(|k| scores[k]);
    let raw_max = raw.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let floored = floor5(raw, other_best);
    let pooled = floored.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    RowScore {
        pooled,
        raw_max,
        winner: winner_mask(&floored, pooled),
        floored: (raw_max - pooled).abs() > 1e-15,
    }
}

#[test]
fn forensic_6r256_joint_java_input_plane() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    kv("target", "20:29455649 T/TGTTTG");
    kv(
        "experiment",
        "2x2: {Rust161, Java174} × {RustClip, JavaClip}; primary = Java174+JavaClip",
    );
    kv(
        "java_gkl",
        "NOT collected: GKL native is x86_64; optional oracle skipped",
    );
    assert_eq!(DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, 128);

    let java_haps: Vec<Vec<u8>> = JAVA_FNV.iter().map(|h| java_seq(h)).collect();
    for (i, seq) in java_haps.iter().enumerate() {
        assert_eq!(fnv1a64_hex(seq), JAVA_FNV[i]);
        assert_eq!(seq.len(), 174);
        assert_eq!(&seq[..9], JAVA_PREFIX);
        assert_eq!(&seq[seq.len() - 4..], JAVA_SUFFIX);
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
    let trim = take_hap_list_trim_span().expect("trim");
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
    for (i, (java, rust)) in java_haps.iter().zip(rust_haps.iter()).enumerate() {
        assert_eq!(fnv1a64_hex(rust), RUST_FNV[i]);
        assert_eq!(rust.len(), 161);
        assert_eq!(&java[9..java.len() - 4], rust.as_slice());
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

    let finalized = finalize_region_reads_for_assembly(
        &region.reads,
        region,
        true,
        gatk_min_tail_quality_for_assembly(10),
        false,
    );
    let rust_map = clip_map(&finalized, region, RUST_CLIP.0, RUST_CLIP.1);
    let java_map = clip_map(&finalized, region, JAVA_CLIP.0, JAVA_CLIP.1);
    kv("n_mapped_rust_clip", rust_map.len().to_string());
    kv("n_mapped_java_clip", java_map.len().to_string());

    let rust_refs: Vec<&[u8]> = rust_haps.iter().map(|s| s.as_slice()).collect();
    let java_refs: Vec<&[u8]> = java_haps.iter().map(|s| s.as_slice()).collect();

    let mut rr = Vec::with_capacity(123);
    let mut jr = Vec::with_capacity(123);
    let mut rj = Vec::with_capacity(123);
    let mut jj = Vec::with_capacity(123);
    let mut n_bases_rust = 0usize;
    let mut n_bases_java = 0usize;

    for (ri, &read_idx) in snap.ad_row_read_index.iter().enumerate() {
        let key = (
            snap.ad_row_qname[ri].as_bytes().to_vec(),
            snap.ad_row_flags[ri],
        );
        let rrec = rust_map
            .get(&key)
            .unwrap_or_else(|| panic!("missing Rust-clip read row={ri}"));
        let jrec = java_map
            .get(&key)
            .unwrap_or_else(|| panic!("missing Java-clip read row={ri}"));
        n_bases_rust += rrec.seq().len();
        n_bases_java += jrec.seq().len();
        let prod_row = by_read.get(&read_idx).expect("prod");
        let p5: [f64; 5] =
            std::array::from_fn(|k| prod_row.haplotype_log10_likelihoods[FROZEN_IDX[k]]);
        let skip = p5.iter().all(|&v| v == 0.0);
        let other_best = prod_row
            .haplotype_log10_likelihoods
            .iter()
            .enumerate()
            .filter(|(i, _)| !FROZEN_IDX.contains(i))
            .map(|(_, v)| *v)
            .fold(f64::NEG_INFINITY, f64::max);
        rr.push(score_row(rrec, &rust_refs, &cfg, other_best, skip));
        jr.push(score_row(rrec, &java_refs, &cfg, other_best, skip));
        rj.push(score_row(jrec, &rust_refs, &cfg, other_best, skip));
        jj.push(score_row(jrec, &java_refs, &cfg, other_best, skip));
        if !skip {
            let pmax = p5.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            assert!(
                (rr[ri].pooled - pmax).abs() < 1e-9,
                "RR must match production five-col max row={ri}"
            );
        }
    }
    kv("total_read_bases_rust_clip", n_bases_rust.to_string());
    kv("total_read_bases_java_clip", n_bases_java.to_string());
    kv("pairhmm_dims_per_config", "123 reads × 5 TGTTTG haps");

    let pack = |rows: &[RowScore]| -> (f64, f64, i32, Vec<f64>) {
        let pooled: Vec<f64> = rows.iter().map(|r| r.pooled).collect();
        let gls = biallelic_gls(&l0, &pooled);
        let (rel, cont, pl) = emitted_hom_alt(&gls);
        (rel, cont, pl, gls)
    };
    let (rr_rel, rr_cont, rr_pl, rr_gls) = pack(&rr);
    let (jr_rel, jr_cont, jr_pl, jr_gls) = pack(&jr);
    let (rj_rel, rj_cont, rj_pl, rj_gls) = pack(&rj);
    let (jj_rel, jj_cont, jj_pl, jj_gls) = pack(&jj);

    for (name, rel, cont, pl, gls) in [
        ("RR", rr_rel, rr_cont, rr_pl, &rr_gls),
        ("JR", jr_rel, jr_cont, jr_pl, &jr_gls),
        ("RJ", rj_rel, rj_cont, rj_pl, &rj_gls),
        ("JJ", jj_rel, jj_cont, jj_pl, &jj_gls),
    ] {
        kv(&format!("{name}_GL_00"), fmt_f64(gls[0]));
        kv(&format!("{name}_GL_02"), fmt_f64(gls[1]));
        kv(&format!("{name}_GL_22"), fmt_f64(gls[2]));
        kv(&format!("{name}_emitted_hom_alt_GL"), fmt_f64(rel));
        kv(&format!("{name}_continuous_PL"), fmt_f64(cont));
        kv(&format!("{name}_integer_PL"), pl.to_string());
    }

    assert!((rr_rel - CTRL_RR_GL).abs() < 1e-9, "Control 1 RR");
    assert_eq!(rr_pl, 3518);
    assert!((jr_rel - CTRL_JR_GL).abs() < 1e-9, "Control 2 JR");
    assert_eq!(jr_pl, 3519);
    assert!((rj_rel - CTRL_RJ_GL).abs() < 1e-9, "Control 3 RJ");
    assert_eq!(rj_pl, 4280);
    kv("controls", "RR=3518 JR=3519 RJ=4280 reproduced");

    let d_hap = jr_rel - rr_rel;
    let d_clip = rj_rel - rr_rel;
    let d_joint = jj_rel - rr_rel;
    let d_int = d_joint - d_hap - d_clip;
    kv("delta_hap_GL", fmt_f64(d_hap));
    kv("delta_clip_GL", fmt_f64(d_clip));
    kv("delta_joint_GL", fmt_f64(d_joint));
    kv("delta_interaction_GL", fmt_f64(d_int));
    kv("delta_hap_PL", fmt_f64(jr_cont - rr_cont));
    kv("delta_clip_PL", fmt_f64(rj_cont - rr_cont));
    kv("delta_joint_PL", fmt_f64(jj_cont - rr_cont));
    kv(
        "delta_interaction_PL",
        fmt_f64((jj_cont - rr_cont) - (jr_cont - rr_cont) - (rj_cont - rr_cont)),
    );
    kv("crosses_3517_5", (jj_cont < 3517.5).to_string());

    let pair_stats = |tag: &str, a: &[RowScore], b: &[RowScore]| {
        let mut n_changed = 0usize;
        let mut n_win = 0usize;
        let mut n_floor = 0usize;
        let mut ds: Vec<f64> = Vec::with_capacity(123);
        for ri in 0..123 {
            let d = b[ri].pooled - a[ri].pooled;
            ds.push(d);
            if d != 0.0 {
                n_changed += 1;
            }
            if b[ri].winner != a[ri].winner {
                n_win += 1;
            }
            if b[ri].floored != a[ri].floored {
                n_floor += 1;
            }
        }
        let mut sorted = ds.clone();
        sorted.sort_by(|x, y| x.partial_cmp(y).unwrap());
        kv(&format!("{tag}_n_changed"), n_changed.to_string());
        kv(&format!("{tag}_n_winner_changed"), n_win.to_string());
        kv(&format!("{tag}_n_floor_transition"), n_floor.to_string());
        kv(&format!("{tag}_delta_min"), fmt_f64(sorted[0]));
        kv(&format!("{tag}_delta_max"), fmt_f64(sorted[122]));
        kv(
            &format!("{tag}_delta_mean"),
            fmt_f64(sorted.iter().sum::<f64>() / 123.0),
        );
        kv(&format!("{tag}_delta_median"), fmt_f64(sorted[61]));
    };
    pair_stats("JR_vs_RR", &rr, &jr);
    pair_stats("RJ_vs_RR", &rr, &rj);
    pair_stats("JJ_vs_JR", &jr, &jj);
    pair_stats("JJ_vs_RJ", &rj, &jj);

    let mut n_changed_jj = 0usize;
    let mut n_win_jj = 0usize;
    let mut n_floor_jj = 0usize;
    let mut n_interact = 0usize;
    let mut n_unique254_changed = 0usize;
    let mut n_unique254_jj_ne_jr = 0usize;
    let mut deltas_jj = Vec::new();
    for ri in 0..123 {
        let d = jj[ri].pooled - rr[ri].pooled;
        if d != 0.0 {
            n_changed_jj += 1;
            deltas_jj.push((ri, d));
        }
        if jj[ri].winner != rr[ri].winner {
            n_win_jj += 1;
        }
        if jj[ri].floored != rr[ri].floored {
            n_floor_jj += 1;
        }
        let pred = (jr[ri].pooled - rr[ri].pooled) + (rj[ri].pooled - rr[ri].pooled);
        if (d - pred).abs() > 1e-12 {
            n_interact += 1;
        }
    }
    for &ri in &UNIQUE_MAX_254 {
        if jj[ri].pooled != rr[ri].pooled {
            n_unique254_changed += 1;
        }
        if jj[ri].pooled != jr[ri].pooled {
            n_unique254_jj_ne_jr += 1;
        }
    }
    let mut all_d: Vec<f64> = (0..123).map(|i| jj[i].pooled - rr[i].pooled).collect();
    all_d.sort_by(|a, b| a.partial_cmp(b).unwrap());
    kv("JJ_vs_RR_n_changed", n_changed_jj.to_string());
    kv("JJ_vs_RR_n_winner_changed", n_win_jj.to_string());
    kv("JJ_vs_RR_n_floor_transition", n_floor_jj.to_string());
    kv("JJ_vs_RR_n_nonadditive_rows", n_interact.to_string());
    kv(
        "unique254_n_changed_JJ_vs_RR",
        n_unique254_changed.to_string(),
    );
    kv("unique254_n_JJ_ne_JR", n_unique254_jj_ne_jr.to_string());
    kv("JJ_vs_RR_delta_min", fmt_f64(all_d[0]));
    kv("JJ_vs_RR_delta_max", fmt_f64(all_d[122]));
    kv(
        "JJ_vs_RR_delta_mean",
        fmt_f64(all_d.iter().sum::<f64>() / 123.0),
    );
    kv("JJ_vs_RR_delta_median", fmt_f64(all_d[61]));
    deltas_jj.sort_by(|a, b| a.1.abs().partial_cmp(&b.1.abs()).unwrap());
    for &(ri, d) in deltas_jj.iter().rev().take(22) {
        kv(
            &format!("JJ_shifted_row_{ri}"),
            format!(
                "qname={} flags={} dJJ={} RR={} JR={} RJ={} JJ={} winRR={:05b} winJJ={:05b} flRR={} flJJ={} rawRR={} rawJJ={}",
                snap.ad_row_qname[ri],
                snap.ad_row_flags[ri],
                fmt_f64(d),
                fmt_f64(rr[ri].pooled),
                fmt_f64(jr[ri].pooled),
                fmt_f64(rj[ri].pooled),
                fmt_f64(jj[ri].pooled),
                rr[ri].winner,
                jj[ri].winner,
                rr[ri].floored,
                jj[ri].floored,
                fmt_f64(rr[ri].raw_max),
                fmt_f64(jj[ri].raw_max),
            ),
        );
    }
    for &ri in &UNIQUE_MAX_254 {
        kv(
            &format!("unique254_row_{ri}"),
            format!(
                "RR={} JR={} RJ={} JJ={} winRR={:05b} winJJ={:05b}",
                fmt_f64(rr[ri].pooled),
                fmt_f64(jr[ri].pooled),
                fmt_f64(rj[ri].pooled),
                fmt_f64(jj[ri].pooled),
                rr[ri].winner,
                jj[ri].winner,
            ),
        );
    }

    let classification = if jj_cont < 3517.5 && jj_pl == 3517 {
        "JOINT_JAVA_INPUT_PL_BOUNDARY_REPRODUCED"
    } else if (jj_rel - rr_rel).abs() > 1e-6 && jj_cont >= 3517.5 && (jj_cont - 3517.5).abs() < 5.0
    {
        "JOINT_JAVA_INPUTS_MOVE_GL_BUT_NOT_PL_BOUNDARY"
    } else {
        "JOINT_JAVA_INPUTS_NOT_SUFFICIENT"
    };
    kv("classification", classification);
    kv(
        "next_arrow",
        if classification == "JOINT_JAVA_INPUT_PL_BOUNDARY_REPRODUCED" {
            "joint sufficiency only; next proof: which Java transform (STR pad vs clip) is required before any production change"
        } else if classification == "JOINT_JAVA_INPUTS_MOVE_GL_BUT_NOT_PL_BOUNDARY" {
            "joint plane moves GL but not across 3517.5; remaining PairHMM input difference"
        } else {
            "joint Java hap+clip plane is not sufficient; first remaining PairHMM input difference"
        },
    );
    assert!(
        classification == "JOINT_JAVA_INPUT_PL_BOUNDARY_REPRODUCED"
            || classification == "JOINT_JAVA_INPUTS_MOVE_GL_BUT_NOT_PL_BOUNDARY"
            || classification == "JOINT_JAVA_INPUTS_NOT_SUFFICIENT"
    );
}
