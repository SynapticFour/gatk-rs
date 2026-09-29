//! 6R.258: measurement-only. First PairHMM backend initialization /
//! probability / recurrence divergence at `20:29455649 T/TGTTTG` under the
//! frozen 6R.257 joint input plane. PRODUCTION CHANGE: NONE.
//!
//! Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77` + GKL 0.8.8
//! `VectorLoglessPairHMM` / `IntelPairHmm` (`useDoublePrecision=false`).
//! Native GKL is not executed (x86_64). Diagnostic scalar reconstructs
//! `Context<float>` + `LoglessPairHMM` recurrence from source.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r258_pairhmm_recurrence -- --nocapture --test-threads=1
//! HOLDOUT_6R258=1 cargo test -p gatk-haplotypecaller --test holdout_6r258_pairhmm_recurrence -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::assembly_region_finalize::{
    clip_finalized_reads_to_region, finalize_region_reads_for_assembly,
    gatk_min_tail_quality_for_assembly,
};
use gatk_haplotypecaller::hc_genotyping_engine::take_colocated_merge_numerics;
use gatk_haplotypecaller::pairhmm_log10::GATK_PARITY_DEFAULT_GCP;
use gatk_haplotypecaller::pcr_error_model::apply_pcr_error_model;
use gatk_haplotypecaller::read_threading_assembler::DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH;
use gatk_haplotypecaller::{
    begin_hap_list_observe, biallelic_genotype_log10_likelihoods_gatk, call_disposition,
    flatten_assembly_regions, indel_gop_from_optional_tag, logless_pairhmm_likelihood,
    prepare_read_quals_for_pairhmm_inplace, region_likelihoods_to_rows,
    score_read_against_haplotypes, take_hap_list_trim_span, traverse_assembly_region_walker,
    AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition, HaplotypeCallerEngine,
    HcLikelihoodEngineConfig, PairHmmBackend, ReadFilterParams, ReadLikelihoodRow,
    WalkerTraversalConfig, INITIAL_CONDITION, INITIAL_CONDITION_LOG10,
};
use rust_htslib::bam::record::Aux;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const GKL_PIN: &str = "0.8.8 IntelLabs/GKL Context.h + avx-pairhmm-template.h";
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
const JOINT_GL: f64 = -351.91571812977724676;
const JAVA_CLIP: (u64, u64) = (29_455_560, 29_455_728);
const RUST_CLIP: (u64, u64) = (29_455_569, 29_455_724);
const JAVA_TSV: &str = include_str!("forensic_6r252_java_hap_trim.tsv");
const BOUNDARY: f64 = 3517.5;
const REQUIRED_MOVE: f64 = 1.65718129777269496;

/// GKL `Context<float>::INITIAL_CONSTANT = ldexpf(1.f, 120.f)`.
const GKL_INITIAL_F32: f32 = f32::from_bits((127 + 120) << 23);
const GKL_MAX_QUAL: usize = 254;
const JACOBIAN_STEP: f64 = 0.0001;
const JACOBIAN_TOL: f32 = 8.0;
const JACOBIAN_INV_STEP: f32 = 10_000.0;
const JACOBIAN_SIZE: usize = 80_001;
const GKL_INV_LN10: f64 = 0.434294;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R258\t{key}\t{}", value.as_ref());
}

fn fmt_f64(x: f64) -> String {
    format!("{x:.17} bits=0x{:016x}", x.to_bits())
}

fn fmt_f32(x: f32) -> String {
    format!(
        "{x:.9} bits=0x{:08x} widened={}",
        x.to_bits(),
        fmt_f64(f64::from(x))
    )
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

fn clip_map(
    finalized: &[rust_htslib::bam::Record],
    region: &gatk_haplotypecaller::assembly_region_iterator::AssemblyRegion,
    start: u64,
    end: u64,
) -> std::collections::HashMap<(Vec<u8>, u16), rust_htslib::bam::Record> {
    let mut clip_region = region.clone();
    clip_region.extended_start = GenomePosition::new_1based(start);
    clip_region.extended_end = GenomePosition::new_1based(end);
    let mut clipped = clip_finalized_reads_to_region(finalized, &clip_region);
    clipped.retain(|r| r.seq().len() >= 10);
    let mut by_id = std::collections::HashMap::new();
    for rec in clipped {
        by_id.insert((rec.qname().to_vec(), rec.flags()), rec);
    }
    by_id
}

struct Planes {
    bases: Vec<u8>,
    bq: Vec<u8>,
    iq: Vec<u8>,
    dq: Vec<u8>,
    gcp: Vec<u8>,
}

fn planes(rec: &rust_htslib::bam::Record, cfg: &HcLikelihoodEngineConfig) -> Planes {
    let bases = rec.seq().as_bytes();
    let n = bases.len();
    let mut bq = rec.qual().to_vec();
    prepare_read_quals_for_pairhmm_inplace(&mut bq, rec.mapq(), cfg);
    let mut iq = indel_gop_from_optional_tag(bam_indel_phred(rec, b"BI").as_deref(), n).unwrap();
    let mut dq = indel_gop_from_optional_tag(bam_indel_phred(rec, b"BD").as_deref(), n).unwrap();
    apply_pcr_error_model(&bases, &mut iq, &mut dq, cfg.pcr_error_model);
    Planes {
        bases,
        bq,
        iq,
        dq,
        gcp: vec![GATK_PARITY_DEFAULT_GCP; n],
    }
}

fn rust_err(q: u8) -> f64 {
    10f64.powf(-(q as f64) / 10.0)
}

fn rust_match_mm(q: u8) -> (f64, f64) {
    let err = rust_err(q);
    (1.0 - err, err / 3.0)
}

fn rust_approx_log10_sum(a: f64, b: f64) -> f64 {
    let (x, y) = if a > b { (b, a) } else { (a, b) };
    if x.is_infinite() && x.is_sign_negative() {
        return y;
    }
    y + (1.0 + 10f64.powf(x - y)).log10()
}

fn rust_m2m(ins: u8, del: u8) -> f64 {
    let (min_q, max_q) = if ins <= del {
        (ins as usize, del as usize)
    } else {
        (del as usize, ins as usize)
    };
    let log10_sum = rust_approx_log10_sum(-0.1 * min_q as f64, -0.1 * max_q as f64);
    let log10_m2m = (1.0 - 10f64.powf(log10_sum).min(1.0)).log10();
    10f64.powf(log10_m2m)
}

fn rust_trans(ins: u8, del: u8, gcp: u8) -> [f64; 6] {
    let gcp_err = rust_err(gcp);
    [
        rust_m2m(ins, del),
        1.0 - rust_err(gcp),
        rust_err(ins),
        gcp_err,
        rust_err(del),
        gcp_err,
    ]
}

struct GklFloatTables {
    ph2pr: [f32; 128],
    m2m: Vec<f32>,
    log10_initial: f32,
}

fn gkl_tables() -> &'static GklFloatTables {
    static T: OnceLock<GklFloatTables> = OnceLock::new();
    T.get_or_init(|| {
        let mut jacobian = vec![0.0f32; JACOBIAN_SIZE];
        for k in 0..JACOBIAN_SIZE {
            jacobian[k] = (1.0 + 10f64.powf(-(k as f64) * JACOBIAN_STEP)).log10() as f32;
        }
        let approx = |small: f32, big: f32| -> f32 {
            let (mut s, mut b) = (small, big);
            if s > b {
                std::mem::swap(&mut s, &mut b);
            }
            if s.is_infinite() || b.is_infinite() {
                return b;
            }
            let diff = b - s;
            if diff >= JACOBIAN_TOL {
                return b;
            }
            let d = diff * JACOBIAN_INV_STEP;
            let ind = if d > 0.0 {
                (d + 0.5) as i32
            } else {
                (d - 0.5) as i32
            } as usize;
            b + jacobian[ind.min(JACOBIAN_SIZE - 1)]
        };
        let mut ph2pr = [0.0f32; 128];
        for x in 0..128 {
            ph2pr[x] = 10f32.powf(-(x as f32) / 10.0);
        }
        let n_m2m = ((GKL_MAX_QUAL + 1) * (GKL_MAX_QUAL + 2)) >> 1;
        let mut m2m = vec![0.0f32; n_m2m];
        let mut offset = 0usize;
        for i in 0..=GKL_MAX_QUAL {
            for j in 0..=i {
                let log10_sum = approx(-0.1f32 * i as f32, -0.1f32 * j as f32) as f64;
                let match_log10 =
                    ((-1.0f64 * (10f64.powf(log10_sum)).min(1.0)).ln_1p()) * GKL_INV_LN10;
                m2m[offset + j] = 10f64.powf(match_log10) as f32;
            }
            offset += i + 1;
        }
        GklFloatTables {
            ph2pr,
            m2m,
            log10_initial: GKL_INITIAL_F32.log10(),
        }
    })
}

fn gkl_m2m(ins: u8, del: u8) -> f32 {
    let t = gkl_tables();
    let (min_q, max_q) = if ins <= del {
        (ins as usize, del as usize)
    } else {
        (del as usize, ins as usize)
    };
    t.m2m[((max_q * (max_q + 1)) >> 1) + min_q]
}

fn gkl_trans(ins: u8, del: u8, gcp: u8) -> [f32; 6] {
    let t = gkl_tables();
    let i = (ins as usize) & 127;
    let d = (del as usize) & 127;
    let c = (gcp as usize) & 127;
    [
        gkl_m2m(ins, del),
        1.0 - t.ph2pr[c],
        t.ph2pr[i],
        t.ph2pr[c],
        t.ph2pr[d],
        t.ph2pr[c],
    ]
}

fn gkl_match_mm(q: u8) -> (f32, f32) {
    let err = gkl_tables().ph2pr[(q as usize) & 127];
    (1.0 - err, err / 3.0)
}

#[derive(Clone, Copy)]
struct Cell {
    m: f64,
    ins: f64,
    del: f64,
}

fn rust_first_cell(read: &[u8], hap: &[u8], bq: u8, iq: u8, dq: u8, gcp: u8, hn: usize) -> Cell {
    let init = INITIAL_CONDITION / hn as f64;
    let t = rust_trans(iq, dq, gcp);
    let (match_p, mismatch_p) = rust_match_mm(bq);
    let p = if read[0] == hap[0] || read[0] == b'N' || hap[0] == b'N' {
        match_p
    } else {
        mismatch_p
    };
    let m = p * (0.0 * t[0] + 0.0 * t[1] + init * t[1]);
    let ins = 0.0 * t[2] + 0.0 * t[3];
    let del = m * t[4] + init * t[5];
    Cell { m, ins, del }
}

fn gkl_first_cell(read: &[u8], hap: &[u8], bq: u8, iq: u8, dq: u8, gcp: u8, hn: usize) -> Cell {
    let init = GKL_INITIAL_F32 / hn as f32;
    let t = gkl_trans(iq, dq, gcp);
    let (match_p, mismatch_p) = gkl_match_mm(bq);
    let p = if read[0] == hap[0] || read[0] == b'N' || hap[0] == b'N' {
        match_p
    } else {
        mismatch_p
    };
    let m = p * (0.0 * t[0] + 0.0 * t[1] + init * t[1]);
    let ins = 0.0 * t[2] + 0.0 * t[3];
    let del = m * t[4] + init * t[5];
    Cell {
        m: f64::from(m),
        ins: f64::from(ins),
        del: f64::from(del),
    }
}

fn rust_diag_ll(
    read: &[u8],
    quals: &[u8],
    hap: &[u8],
    iq: &[u8],
    dq: &[u8],
    gcp: &[u8],
    init_del: f64,
    match_mm: impl Fn(u8) -> (f64, f64),
    trans: impl Fn(u8, u8, u8) -> [f64; 6],
) -> f64 {
    let rn = read.len();
    let hn = hap.len();
    let cols = hn + 1;
    let mut m = vec![0.0f64; (rn + 1) * cols];
    let mut ins = vec![0.0f64; (rn + 1) * cols];
    let mut del = vec![0.0f64; (rn + 1) * cols];
    for j in 0..=hn {
        del[j] = init_del;
    }
    for i in 1..=rn {
        let t = trans(iq[i - 1], dq[i - 1], gcp[i - 1]);
        let x = read[i - 1];
        let (match_p, mismatch_p) = match_mm(quals[i - 1]);
        let row = i * cols;
        let prev = (i - 1) * cols;
        m[row] = 0.0;
        ins[row] = 0.0;
        del[row] = 0.0;
        for j in 1..=hn {
            let y = hap[j - 1];
            let p = if x == y || x == b'N' || y == b'N' {
                match_p
            } else {
                mismatch_p
            };
            m[row + j] =
                p * (m[prev + j - 1] * t[0] + ins[prev + j - 1] * t[1] + del[prev + j - 1] * t[1]);
            ins[row + j] = m[prev + j] * t[2] + ins[prev + j] * t[3];
            del[row + j] = m[row + j - 1] * t[4] + del[row + j - 1] * t[5];
        }
    }
    let end = rn * cols;
    let mut sum = 0.0;
    for j in 1..=hn {
        sum += m[end + j] + ins[end + j];
    }
    if sum <= 0.0 || !sum.is_finite() {
        return f64::NEG_INFINITY;
    }
    sum.log10() - INITIAL_CONDITION_LOG10
}

fn gkl_float_ll(read: &[u8], quals: &[u8], hap: &[u8], iq: &[u8], dq: &[u8], gcp: &[u8]) -> f64 {
    let rn = read.len();
    let hn = hap.len();
    let cols = hn + 1;
    let init = GKL_INITIAL_F32 / hn as f32;
    let mut m = vec![0.0f32; (rn + 1) * cols];
    let mut ins = vec![0.0f32; (rn + 1) * cols];
    let mut del = vec![0.0f32; (rn + 1) * cols];
    for j in 0..=hn {
        del[j] = init;
    }
    for i in 1..=rn {
        let t = gkl_trans(iq[i - 1], dq[i - 1], gcp[i - 1]);
        let x = read[i - 1];
        let (match_p, mismatch_p) = gkl_match_mm(quals[i - 1]);
        let row = i * cols;
        let prev = (i - 1) * cols;
        m[row] = 0.0;
        ins[row] = 0.0;
        del[row] = 0.0;
        for j in 1..=hn {
            let y = hap[j - 1];
            let p = if x == y || x == b'N' || y == b'N' {
                match_p
            } else {
                mismatch_p
            };
            m[row + j] =
                p * (m[prev + j - 1] * t[0] + ins[prev + j - 1] * t[1] + del[prev + j - 1] * t[1]);
            ins[row + j] = m[prev + j] * t[2] + ins[prev + j] * t[3];
            del[row + j] = m[row + j - 1] * t[4] + del[row + j - 1] * t[5];
        }
    }
    let end = rn * cols;
    let mut sum = 0.0f32;
    for j in 1..=hn {
        sum += m[end + j] + ins[end + j];
    }
    if sum <= 0.0 || !sum.is_finite() {
        return f64::NEG_INFINITY;
    }
    f64::from(sum.log10() - gkl_tables().log10_initial)
}

#[test]
fn forensic_6r258_pairhmm_recurrence() {
    kv("java_pin", JAVA_PIN);
    kv("gkl_pin", GKL_PIN);
    kv("production_change", "NONE");
    kv("target", "20:29455649 T/TGTTTG");
    kv(
        "frozen_input",
        "6R.257: Java174-mers + Java clip 20:29455560-29455728; 123 retainEvidence",
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

    let rust_haps: Vec<Vec<u8>> = FROZEN_IDX
        .iter()
        .map(|&idx| outcome.assembly.haplotypes[idx].bases.clone())
        .collect();
    for (i, (java, rust)) in java_haps.iter().zip(rust_haps.iter()).enumerate() {
        assert_eq!(fnv1a64_hex(rust), RUST_FNV[i]);
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
    kv(
        "java_gkl_path",
        "PairHMM.FASTEST_AVAILABLE → VectorLoglessPairHMM(AVX) → IntelPairHmm.initialize(useDoublePrecision=false) → g_compute_full_prob_float",
    );
    kv(
        "rust_path",
        "FastestAvailable → NeonF64 → score_haps_neon_f64 / logless_fill_transitions / score_pack2",
    );

    let retain: std::collections::HashSet<usize> = snap.ad_row_read_index.iter().copied().collect();
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
    let l0: Vec<f64> = snap.ad_row_lls.iter().map(|ll| ll[0]).collect();

    let finalized = finalize_region_reads_for_assembly(
        &region.reads,
        region,
        true,
        gatk_min_tail_quality_for_assembly(10),
        false,
    );
    let java_map = clip_map(&finalized, region, JAVA_CLIP.0, JAVA_CLIP.1);
    let java_refs: Vec<&[u8]> = java_haps.iter().map(|s| s.as_slice()).collect();
    let hn = java_haps[0].len();
    assert_eq!(hn, 174);

    kv("gkl_initial_f32", fmt_f32(GKL_INITIAL_F32));
    kv("rust_initial_f64", fmt_f64(INITIAL_CONDITION));
    kv("gkl_log10_initial", fmt_f32(gkl_tables().log10_initial));
    kv("rust_log10_initial", fmt_f64(INITIAL_CONDITION_LOG10));
    let gkl_init = GKL_INITIAL_F32 / hn as f32;
    let rust_init = INITIAL_CONDITION / hn as f64;
    let gkl_quot = f64::from(1.0f32 / hn as f32);
    let rust_quot = 1.0 / hn as f64;
    kv("gkl_init_del_raw", fmt_f32(gkl_init));
    kv("rust_init_del_raw", fmt_f64(rust_init));
    kv("gkl_init_del_over_initial", fmt_f64(gkl_quot));
    kv("rust_init_del_over_initial", fmt_f64(rust_quot));
    kv(
        "init_policy",
        "MATCH: M[0]=0 I[0]=0 D[0][j]=INITIAL/haplen col0=0 free leading deletions",
    );

    let mut first_op = "none";
    let mut first_detail = String::new();
    if GKL_INITIAL_F32.to_bits() != 0
        && INITIAL_CONDITION.to_bits() != GKL_INITIAL_F32 as f64 as u64
    {
        // raw INITIAL always differs by design; comparable layer is haplen quotient
    }
    if gkl_quot.to_bits() != rust_quot.to_bits() {
        first_op = "GKL_FLOAT_HAPLEN_QUOTIENT";
        first_detail = format!(
            "init_del/INITIAL gkl_f32_widened={} rust_f64={} abs_delta={:.17} haplen=174",
            fmt_f64(gkl_quot),
            fmt_f64(rust_quot),
            (gkl_quot - rust_quot).abs()
        );
    }

    let mut rust_pooled = vec![0.0f64; 123];
    let mut cf_pooled = vec![0.0f64; 123];
    let mut gkl_pooled = vec![0.0f64; 123];
    let mut traced = false;
    let mut n_ph2pr_diff = 0usize;
    let mut n_m2m_diff = 0usize;
    let mut n_trans_diff = 0usize;
    let mut n_neon_vs_scalar = 0usize;
    let mut max_neon_vs_scalar = 0.0f64;
    let mut max_gkl_vs_neon = 0.0f64;
    let mut first_material_row = None;

    for (ri, &read_idx) in snap.ad_row_read_index.iter().enumerate() {
        let key = (
            snap.ad_row_qname[ri].as_bytes().to_vec(),
            snap.ad_row_flags[ri],
        );
        let rec = java_map
            .get(&key)
            .unwrap_or_else(|| panic!("missing Java-clip read row={ri}"));
        let p = planes(rec, &cfg);
        let prod_row = by_read.get(&read_idx).expect("prod");
        let p5: [f64; 5] =
            std::array::from_fn(|k| prod_row.haplotype_log10_likelihoods[FROZEN_IDX[k]]);
        if p5.iter().all(|&v| v == 0.0) {
            continue;
        }
        let other_best = prod_row
            .haplotype_log10_likelihoods
            .iter()
            .enumerate()
            .filter(|(i, _)| !FROZEN_IDX.contains(i))
            .map(|(_, v)| *v)
            .fold(f64::NEG_INFINITY, f64::max);

        let scores = score_read_against_haplotypes(
            &cfg,
            &p.bases,
            rec.qual(),
            rec.mapq(),
            &java_refs,
            bam_indel_phred(rec, b"BI").as_deref(),
            bam_indel_phred(rec, b"BD").as_deref(),
        )
        .expect("neon");
        let raw: [f64; 5] = std::array::from_fn(|k| scores[k]);
        rust_pooled[ri] = floor5(raw, other_best)
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);

        let mut gkl_raw = [0.0f64; 5];
        let mut cf_raw = [0.0f64; 5];
        let cf_init = INITIAL_CONDITION * gkl_quot;
        for k in 0..5 {
            gkl_raw[k] = gkl_float_ll(&p.bases, &p.bq, &java_haps[k], &p.iq, &p.dq, &p.gcp);
            cf_raw[k] = rust_diag_ll(
                &p.bases,
                &p.bq,
                &java_haps[k],
                &p.iq,
                &p.dq,
                &p.gcp,
                cf_init,
                rust_match_mm,
                rust_trans,
            );
            let scalar =
                logless_pairhmm_likelihood(&p.bases, &p.bq, &java_haps[k], &p.iq, &p.dq, &p.gcp)
                    .expect("scalar");
            let d_ns = (scalar - raw[k]).abs();
            if d_ns > 0.0 {
                n_neon_vs_scalar += 1;
            }
            max_neon_vs_scalar = max_neon_vs_scalar.max(d_ns);
            max_gkl_vs_neon = max_gkl_vs_neon.max((gkl_raw[k] - raw[k]).abs());
        }
        gkl_pooled[ri] = floor5(gkl_raw, other_best)
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        cf_pooled[ri] = floor5(cf_raw, other_best)
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);

        for i in 0..p.bq.len() {
            let (rm, rmm) = rust_match_mm(p.bq[i]);
            let (gm, gmm) = gkl_match_mm(p.bq[i]);
            if rm.to_bits() != f64::from(gm).to_bits() || rmm.to_bits() != f64::from(gmm).to_bits()
            {
                n_ph2pr_diff += 1;
            }
            let rt = rust_trans(p.iq[i], p.dq[i], p.gcp[i]);
            let gt = gkl_trans(p.iq[i], p.dq[i], p.gcp[i]);
            if rust_m2m(p.iq[i], p.dq[i]).to_bits()
                != f64::from(gkl_m2m(p.iq[i], p.dq[i])).to_bits()
            {
                n_m2m_diff += 1;
            }
            for c in 0..6 {
                if rt[c].to_bits() != f64::from(gt[c]).to_bits() {
                    n_trans_diff += 1;
                    break;
                }
            }
        }

        if !traced {
            traced = true;
            let hap = &java_haps[0];
            kv("trace_read_index", ri.to_string());
            kv("trace_qname", &snap.ad_row_qname[ri]);
            kv("trace_flags", snap.ad_row_flags[ri].to_string());
            kv("trace_haplotype_index", "19");
            kv("trace_haplotype_fnv", JAVA_FNV[0]);
            kv("trace_read_length", p.bases.len().to_string());
            kv("trace_haplotype_length", hap.len().to_string());
            kv("trace_read_base0", format!("{}", p.bases[0] as char));
            kv("trace_hap_base0", format!("{}", hap[0] as char));
            kv("trace_bq0", p.bq[0].to_string());
            kv("trace_iq0", p.iq[0].to_string());
            kv("trace_dq0", p.dq[0].to_string());
            kv("trace_gcp0", p.gcp[0].to_string());

            let (rm, rmm) = rust_match_mm(p.bq[0]);
            let (gm, gmm) = gkl_match_mm(p.bq[0]);
            kv("rust_match_p0", fmt_f64(rm));
            kv("gkl_match_p0", fmt_f32(gm));
            kv("rust_mismatch_p0", fmt_f64(rmm));
            kv("gkl_mismatch_p0", fmt_f32(gmm));
            kv("rust_ph2pr_q0", fmt_f64(rust_err(p.bq[0])));
            kv(
                "gkl_ph2pr_q0",
                fmt_f32(gkl_tables().ph2pr[(p.bq[0] as usize) & 127]),
            );
            let rt = rust_trans(p.iq[0], p.dq[0], p.gcp[0]);
            let gt = gkl_trans(p.iq[0], p.dq[0], p.gcp[0]);
            for (name, i) in [
                ("matchToMatch", 0),
                ("indelToMatch", 1),
                ("matchToInsertion", 2),
                ("insertionToInsertion", 3),
                ("matchToDeletion", 4),
                ("deletionToDeletion", 5),
            ] {
                kv(&format!("rust_trans0_{name}"), fmt_f64(rt[i]));
                kv(&format!("gkl_trans0_{name}"), fmt_f32(gt[i]));
            }

            let rc = rust_first_cell(&p.bases, hap, p.bq[0], p.iq[0], p.dq[0], p.gcp[0], hn);
            let gc = gkl_first_cell(&p.bases, hap, p.bq[0], p.iq[0], p.dq[0], p.gcp[0], hn);
            kv("rust_first_cell_M", fmt_f64(rc.m));
            kv("gkl_first_cell_M", fmt_f64(gc.m));
            kv("rust_first_cell_I", fmt_f64(rc.ins));
            kv("gkl_first_cell_I", fmt_f64(gc.ins));
            kv("rust_first_cell_D", fmt_f64(rc.del));
            kv("gkl_first_cell_D", fmt_f64(gc.del));
            kv(
                "first_cell_M_scale_normalized_rust",
                fmt_f64(rc.m / INITIAL_CONDITION),
            );
            kv(
                "first_cell_M_scale_normalized_gkl",
                fmt_f64(gc.m / f64::from(GKL_INITIAL_F32)),
            );
            kv("neon_ll_hap0", fmt_f64(raw[0]));
            kv(
                "scalar_logless_ll_hap0",
                fmt_f64(
                    logless_pairhmm_likelihood(&p.bases, &p.bq, hap, &p.iq, &p.dq, &p.gcp).unwrap(),
                ),
            );
            kv("gkl_float_diag_ll_hap0", fmt_f64(gkl_raw[0]));
            kv("cf_haplen_quot_ll_hap0", fmt_f64(cf_raw[0]));

            if first_op == "none" {
                if f64::from(gm).to_bits() != rm.to_bits() {
                    first_op = "GKL_FLOAT_PH2PR_MATCH_PRIOR";
                    first_detail = format!(
                        "read_pos=0 bq={} gkl_match={} rust_match={}",
                        p.bq[0],
                        fmt_f32(gm),
                        fmt_f64(rm)
                    );
                } else if f64::from(gt[0]).to_bits() != rt[0].to_bits() {
                    first_op = "GKL_FLOAT_MATCH_TO_MATCH_TABLE";
                    first_detail = format!(
                        "read_pos=0 iq={} dq={} gkl={} rust={}",
                        p.iq[0],
                        p.dq[0],
                        fmt_f32(gt[0]),
                        fmt_f64(rt[0])
                    );
                } else if (rc.m / INITIAL_CONDITION).to_bits()
                    != (gc.m / f64::from(GKL_INITIAL_F32)).to_bits()
                {
                    first_op = "FIRST_DP_CELL_MATCH_STATE";
                    first_detail = format!(
                        "i=1 j=1 scale-normalized M gkl={} rust={}",
                        fmt_f64(gc.m / f64::from(GKL_INITIAL_F32)),
                        fmt_f64(rc.m / INITIAL_CONDITION)
                    );
                }
            }
        }

        if first_material_row.is_none() {
            let d = (gkl_raw[0] - raw[0]).abs();
            if d > 1e-8 {
                first_material_row = Some((ri, d));
            }
        }
    }

    kv("n_ph2pr_positions_bit_unequal", n_ph2pr_diff.to_string());
    kv("n_m2m_positions_bit_unequal", n_m2m_diff.to_string());
    kv(
        "n_trans_rows_any_component_unequal",
        n_trans_diff.to_string(),
    );
    kv(
        "n_neon_vs_scalar_cells_unequal",
        n_neon_vs_scalar.to_string(),
    );
    kv("max_abs_neon_vs_scalar", fmt_f64(max_neon_vs_scalar));
    kv("max_abs_gkl_float_vs_neon", fmt_f64(max_gkl_vs_neon));
    if let Some((ri, d)) = first_material_row {
        kv(
            "first_material_gkl_vs_neon_row",
            format!("ri={ri} abs_delta_ll={:.17}", d),
        );
    }

    let rust_gls = biallelic_gls(&l0, &rust_pooled);
    let (rr_rel, rr_cont, rr_pl) = emitted_hom_alt(&rust_gls);
    kv("joint_control_emitted_GL", fmt_f64(rr_rel));
    kv("joint_control_continuous_PL", fmt_f64(rr_cont));
    kv("joint_control_integer_PL", rr_pl.to_string());
    assert!((rr_rel - JOINT_GL).abs() < 1e-9);
    assert_eq!(rr_pl, 3519);

    let cf_gls = biallelic_gls(&l0, &cf_pooled);
    let (cf_rel, cf_cont, cf_pl) = emitted_hom_alt(&cf_gls);
    kv("counterfactual_component", first_op);
    kv("counterfactual_emitted_GL", fmt_f64(cf_rel));
    kv("counterfactual_continuous_PL", fmt_f64(cf_cont));
    kv("counterfactual_integer_PL", cf_pl.to_string());
    kv(
        "counterfactual_delta_continuous_PL",
        fmt_f64(cf_cont - rr_cont),
    );
    kv(
        "counterfactual_crosses_3517_5",
        (cf_cont < BOUNDARY).to_string(),
    );

    let gkl_gls = biallelic_gls(&l0, &gkl_pooled);
    let (gkl_rel, gkl_cont, gkl_pl) = emitted_hom_alt(&gkl_gls);
    kv("gkl_float_diag_full_kernel_emitted_GL", fmt_f64(gkl_rel));
    kv(
        "gkl_float_diag_full_kernel_continuous_PL",
        fmt_f64(gkl_cont),
    );
    kv("gkl_float_diag_full_kernel_integer_PL", gkl_pl.to_string());
    kv(
        "note_full_gkl_kernel",
        "informational only; not the single-component counterfactual",
    );

    kv("first_divergent_operation", first_op);
    kv("first_divergent_detail", first_detail);
    kv("required_pl_movement", format!("{REQUIRED_MOVE:.17}"));

    let classification = if first_op == "none" {
        "PAIRHMM_ALGORITHMIC_SEMANTICS_MATCH"
    } else if cf_cont < BOUNDARY && cf_pl == 3517 {
        "PAIRHMM_RECURRENCE_DIVERGENCE_CAUSALLY_MOVES_PL"
    } else {
        "PAIRHMM_RECURRENCE_DIVERGENCE_NOT_SUFFICIENT"
    };
    kv("classification", classification);
    kv("production_change", "NONE");
    kv(
        "next_arrow",
        if classification == "PAIRHMM_ALGORITHMIC_SEMANTICS_MATCH" {
            "accumulated numerical error of identical recurrence without changing production precision"
        } else if classification == "PAIRHMM_RECURRENCE_DIVERGENCE_CAUSALLY_MOVES_PL" {
            "first init/prob operation is PL-causal; production still unchanged this round"
        } else {
            "first init/prob operation is present but does not cross 3517.5; do not patch production"
        },
    );
    assert!(
        classification == "PAIRHMM_RECURRENCE_DIVERGENCE_CAUSALLY_MOVES_PL"
            || classification == "PAIRHMM_RECURRENCE_DIVERGENCE_NOT_SUFFICIENT"
            || classification == "PAIRHMM_ALGORITHMIC_SEMANTICS_MATCH"
    );
    assert_ne!(first_op, "none");
}
