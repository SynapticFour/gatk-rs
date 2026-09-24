//! 6R.262: measurement-only. GKL float `ph2pr / 3.f` mismatch prior vs
//! Rust `err / 3.0` at `20:29455649 T/TGTTTG` under the frozen 6R.257
//! joint input plane. PRODUCTION CHANGE: NONE.
//!
//! 6R.259–261 closed exponent / powf / 1-ph2pr (not causal). This round
//! injects only the GKL f32 mismatch `ph2pr/3.f` widened to f64. Match
//! stays Rust `1.0-err`. Java ph2pr is not substituted as match.
//!
//! Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77` + GKL 0.8.8
//! `avx-pairhmm-template.h` `distm = VEC_DIV(distm, 3.0)` after loading
//! `distm = ph2pr[q]` (`Context<float>`).
//! Native AVX GKL is not executed (x86_64).
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r262_gkl_float_ph2pr_div3 -- --nocapture --test-threads=1
//! HOLDOUT_6R262=1 cargo test -p gatk-haplotypecaller --test holdout_6r262_gkl_float_ph2pr_div3 -- --nocapture --test-threads=1
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
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const GKL_PIN: &str = "0.8.8 Intel-HLS/GKL avx-pairhmm-template.h distm=VEC_DIV(distm, 3.0)";
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

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R262\t{key}\t{}", value.as_ref());
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

fn rust_match(q: u8) -> f64 {
    1.0 - rust_err(q)
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

fn libc_powf(base: f32, exp: f32) -> f32 {
    extern "C" {
        fn powf(x: f32, y: f32) -> f32;
    }
    unsafe { powf(base, exp) }
}

fn gkl_exponent(q: u8) -> f32 {
    -((q as usize & 127) as f32) / 10.0f32
}

fn gkl_ph2pr(q: u8) -> f32 {
    libc_powf(10.0, gkl_exponent(q))
}

/// GKL AVX `distm = VEC_DIV(distm, 3.0)` with `NUMBER=float`.
fn gkl_mismatch(q: u8) -> f32 {
    gkl_ph2pr(q) / 3.0f32
}

fn rust_mismatch(q: u8) -> f64 {
    rust_err(q) / 3.0
}

fn classify_mismatch(q: u8) -> &'static str {
    let gkl_w = f64::from(gkl_mismatch(q));
    let rust = rust_mismatch(q);
    let f64_div = f64::from(gkl_ph2pr(q)) / 3.0;
    if gkl_w.to_bits() == rust.to_bits() {
        "exact-match"
    } else if gkl_w.to_bits() != f64_div.to_bits() {
        "new-f32-division-rounding"
    } else {
        "inherited-only"
    }
}

/// CF1: Rust match, GKL f32 mismatch widened. Java ph2pr not used as match.
fn cf1_gkl_mismatch_only(q: u8) -> (f64, f64) {
    (rust_match(q), f64::from(gkl_mismatch(q)))
}

/// CF2: Rust match, f64(GKL ph2pr)/3.0. Isolates inherited ph2pr from /3.f.
fn cf2_gkl_ph2pr_f64_div3(q: u8) -> (f64, f64) {
    (rust_match(q), f64::from(gkl_ph2pr(q)) / 3.0)
}

fn rust_diag_ll(
    read: &[u8],
    quals: &[u8],
    hap: &[u8],
    iq: &[u8],
    dq: &[u8],
    gcp: &[u8],
    match_mm: impl Fn(u8) -> (f64, f64),
) -> f64 {
    let rn = read.len();
    let hn = hap.len();
    let cols = hn + 1;
    let init_del = INITIAL_CONDITION / hn as f64;
    let mut m = vec![0.0f64; (rn + 1) * cols];
    let mut ins = vec![0.0f64; (rn + 1) * cols];
    let mut del = vec![0.0f64; (rn + 1) * cols];
    for j in 0..=hn {
        del[j] = init_del;
    }
    for i in 1..=rn {
        let t = rust_trans(iq[i - 1], dq[i - 1], gcp[i - 1]);
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

fn dump_mismatch_row(q: u8) {
    let gkl_p = gkl_ph2pr(q);
    let gkl_mm = gkl_mismatch(q);
    let gkl_mm_w = f64::from(gkl_mm);
    let rust_mm = rust_mismatch(q);
    let f64_div = f64::from(gkl_p) / 3.0;
    let abs = (gkl_mm_w - rust_mm).abs();
    let rel = if rust_mm != 0.0 {
        abs / rust_mm.abs()
    } else {
        0.0
    };
    kv(
        "mismatch_row",
        format!(
            "q={q}\tclass={}\tgkl_ph2pr_f32_bits=0x{:08x}\tgkl_ph2pr_widened={}\tgkl_mismatch_f32_bits=0x{:08x}\tgkl_mismatch_widened={}\trust_mismatch_f64_bits=0x{:016x}\trust_mismatch={}\tabs_delta={:.17e}\trel_delta={:.17e}\tf64_div_same_gkl_ph2pr={}",
            classify_mismatch(q),
            gkl_p.to_bits(),
            fmt_f64(f64::from(gkl_p)),
            gkl_mm.to_bits(),
            fmt_f64(gkl_mm_w),
            rust_mm.to_bits(),
            fmt_f64(rust_mm),
            abs,
            rel,
            fmt_f64(f64_div),
        ),
    );
}

#[test]
fn forensic_6r262_gkl_float_ph2pr_div3() {
    kv("java_pin", JAVA_PIN);
    kv("gkl_pin", GKL_PIN);
    kv("production_change", "NONE");
    kv("target", "20:29455649 T/TGTTTG");
    kv(
        "frozen_input",
        "6R.257: Java174-mers + Java clip 20:29455560-29455728; 123 retainEvidence",
    );
    kv(
        "gkl_mismatch_op",
        "Context<float> stores only ph2pr; stripeINITIALIZATION: distm=ph2pr[q]; distm=VEC_DIV(distm, 3.0); NUMBER=float; packed3=VEC_SET1_VAL(3.0); no widen before div; no wrapper; no later recompute",
    );
    kv(
        "gkl_source",
        "Intel-HLS/GKL 0.8.8 src/main/native/pairhmm/Context.h Context<float> ph2pr; avx-pairhmm-template.h stripeINITIALIZATION VEC_DIV",
    );
    kv(
        "rust_formula",
        "match_mismatch_prior_table: err=10f64.powf(-(q as f64)/10.0); match=1.0-err (unchanged); mismatch=err/3.0 (f64)",
    );
    kv(
        "rust_source",
        "gatk-haplotypecaller/src/pairhmm_logless.rs match_mismatch_prior_table / logless_match_mismatch_prior",
    );
    kv(
        "counterfactual_rule",
        "CF1 injects f64::from(gkl_ph2pr/3.f) as mismatch only; CF2 injects f64(gkl_ph2pr)/3.0; match stays Rust; no other substitution",
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
    assert_eq!(TARGET_REF, "T");
    assert_eq!(TARGET_ALT, "TGTTTG");

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

    let mut unique_q = BTreeSet::new();
    let mut n_bq_total = 0usize;
    let mut rust_pooled = vec![0.0f64; 123];
    let mut cf1_pooled = vec![0.0f64; 123];
    let mut cf2_pooled = vec![0.0f64; 123];

    for (ri, &read_idx) in snap.ad_row_read_index.iter().enumerate() {
        let key = (
            snap.ad_row_qname[ri].as_bytes().to_vec(),
            snap.ad_row_flags[ri],
        );
        let rec = java_map
            .get(&key)
            .unwrap_or_else(|| panic!("missing Java-clip read row={ri}"));
        let p = planes(rec, &cfg);
        n_bq_total += p.bq.len();
        unique_q.extend(p.bq.iter().copied());

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

        let mut cf1_raw = [0.0f64; 5];
        let mut cf2_raw = [0.0f64; 5];
        for k in 0..5 {
            cf1_raw[k] = rust_diag_ll(
                &p.bases,
                &p.bq,
                &java_haps[k],
                &p.iq,
                &p.dq,
                &p.gcp,
                cf1_gkl_mismatch_only,
            );
            cf2_raw[k] = rust_diag_ll(
                &p.bases,
                &p.bq,
                &java_haps[k],
                &p.iq,
                &p.dq,
                &p.gcp,
                cf2_gkl_ph2pr_f64_div3,
            );
            let scalar =
                logless_pairhmm_likelihood(&p.bases, &p.bq, &java_haps[k], &p.iq, &p.dq, &p.gcp)
                    .expect("scalar");
            assert_eq!(
                scalar.to_bits(),
                raw[k].to_bits(),
                "NEON must remain bit-identical to scalar (6R.258)"
            );
        }
        cf1_pooled[ri] = floor5(cf1_raw, other_best)
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        cf2_pooled[ri] = floor5(cf2_raw, other_best)
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
    }

    kv(
        "unique_bq_values",
        unique_q
            .iter()
            .map(|q| q.to_string())
            .collect::<Vec<_>>()
            .join(","),
    );
    kv("n_unique_bq", unique_q.len().to_string());
    kv("n_bq_total", n_bq_total.to_string());

    let mut n_exact = 0usize;
    let mut n_inherited = 0usize;
    let mut n_new_div = 0usize;
    let mut first_mismatch_unequal: Option<u8> = None;
    let mut first_new_div: Option<u8> = None;
    for &q in &unique_q {
        dump_mismatch_row(q);
        let class = classify_mismatch(q);
        match class {
            "exact-match" => n_exact += 1,
            "inherited-only" => n_inherited += 1,
            "new-f32-division-rounding" => n_new_div += 1,
            _ => {}
        }
        if class != "exact-match" && first_mismatch_unequal.is_none() {
            first_mismatch_unequal = Some(q);
        }
        if class == "new-f32-division-rounding" && first_new_div.is_none() {
            first_new_div = Some(q);
        }
    }
    kv("n_unique_bq_exact_match", n_exact.to_string());
    kv("n_unique_bq_inherited_only", n_inherited.to_string());
    kv(
        "n_unique_bq_new_f32_division_rounding",
        n_new_div.to_string(),
    );
    kv(
        "first_differing_mismatch_q",
        first_mismatch_unequal
            .map(|q| q.to_string())
            .unwrap_or_else(|| "none".into()),
    );
    kv(
        "first_new_f32_division_rounding_q",
        first_new_div
            .map(|q| q.to_string())
            .unwrap_or_else(|| "none".into()),
    );

    assert!(unique_q.contains(&20), "frozen object must include Q=20");
    assert!(unique_q.contains(&30), "frozen object must include Q=30");
    kv("q20_gkl_ph2pr", fmt_f32(gkl_ph2pr(20)));
    kv("q20_gkl_mismatch", fmt_f32(gkl_mismatch(20)));
    kv("q20_rust_mismatch", fmt_f64(rust_mismatch(20)));
    kv(
        "q20_f64_div_same_gkl_ph2pr",
        fmt_f64(f64::from(gkl_ph2pr(20)) / 3.0),
    );
    kv("q20_class", classify_mismatch(20));
    kv(
        "q20_mismatch_abs_delta",
        format!(
            "{:.17e}",
            (f64::from(gkl_mismatch(20)) - rust_mismatch(20)).abs()
        ),
    );
    kv("q30_gkl_ph2pr", fmt_f32(gkl_ph2pr(30)));
    kv("q30_gkl_mismatch", fmt_f32(gkl_mismatch(30)));
    kv("q30_rust_mismatch", fmt_f64(rust_mismatch(30)));
    kv("q30_class", classify_mismatch(30));
    kv(
        "q30_mismatch_abs_delta",
        format!(
            "{:.17e}",
            (f64::from(gkl_mismatch(30)) - rust_mismatch(30)).abs()
        ),
    );

    kv(
        "counterfactual_note",
        "NEON-equivalent scalar (6R.258 bit-identical); CF1 replaces only mismatch with f64::from(ph2pr/3.f); CF2 replaces only mismatch with f64(ph2pr)/3.0; match/transitions/init/log10 unchanged",
    );

    let rust_gls = biallelic_gls(&l0, &rust_pooled);
    let (rr_rel, rr_cont, rr_pl) = emitted_hom_alt(&rust_gls);
    kv("baseline_emitted_GL", fmt_f64(rr_rel));
    kv("baseline_continuous_PL", fmt_f64(rr_cont));
    kv("baseline_integer_PL", rr_pl.to_string());
    assert!((rr_rel - JOINT_GL).abs() < 1e-9);
    assert_eq!(rr_pl, 3519);

    let cf1_gls = biallelic_gls(&l0, &cf1_pooled);
    let (cf1_rel, cf1_cont, cf1_pl) = emitted_hom_alt(&cf1_gls);
    kv("first_divergent_operation", "GKL_FLOAT_PH2PR_DIV3");
    kv(
        "first_divergent_detail",
        format!(
            "first_mismatch_q={} first_new_div_q={} q20_class={} q20_gkl_mm={} q20_rust_mm={}",
            first_mismatch_unequal
                .map(|q| q.to_string())
                .unwrap_or_else(|| "none".into()),
            first_new_div
                .map(|q| q.to_string())
                .unwrap_or_else(|| "none".into()),
            classify_mismatch(20),
            fmt_f32(gkl_mismatch(20)),
            fmt_f64(rust_mismatch(20))
        ),
    );
    kv("cf1_emitted_GL", fmt_f64(cf1_rel));
    kv("cf1_continuous_PL", fmt_f64(cf1_cont));
    kv("cf1_integer_PL", cf1_pl.to_string());
    kv("cf1_delta_GL", fmt_f64(cf1_rel - rr_rel));
    kv("cf1_delta_continuous_PL", fmt_f64(cf1_cont - rr_cont));
    kv("cf1_delta_integer_PL", (cf1_pl - rr_pl).to_string());
    kv("required_pl_movement", format!("{REQUIRED_MOVE:.17}"));
    kv("cf1_crosses_3517_5", (cf1_cont < BOUNDARY).to_string());

    let cf2_gls = biallelic_gls(&l0, &cf2_pooled);
    let (cf2_rel, cf2_cont, cf2_pl) = emitted_hom_alt(&cf2_gls);
    kv("cf2_emitted_GL", fmt_f64(cf2_rel));
    kv("cf2_continuous_PL", fmt_f64(cf2_cont));
    kv("cf2_integer_PL", cf2_pl.to_string());
    kv("cf2_delta_GL", fmt_f64(cf2_rel - rr_rel));
    kv("cf2_delta_continuous_PL", fmt_f64(cf2_cont - rr_cont));
    kv("cf2_delta_integer_PL", (cf2_pl - rr_pl).to_string());
    kv("cf2_crosses_3517_5", (cf2_cont < BOUNDARY).to_string());

    let n_unequal = n_inherited + n_new_div;
    let classification = if n_unequal == 0 {
        "PAIRHMM_MISMATCH_PRIOR_RESULT_MATCH"
    } else if cf1_cont < BOUNDARY && cf1_pl == 3517 {
        "PAIRHMM_MISMATCH_PRIOR_DIVERGENCE_CAUSAL"
    } else {
        "PAIRHMM_MISMATCH_PRIOR_DIVERGENCE_NOT_CAUSAL"
    };
    kv("classification", classification);
    kv("production_change", "NONE");
    kv(
        "next_arrow",
        if classification == "PAIRHMM_MISMATCH_PRIOR_DIVERGENCE_CAUSAL" {
            "ph2pr/3.f mismatch is PL-causal; production still unchanged this round"
        } else {
            "GKL_AVX_DISTM_BLEND: VEC_BLENDV(distm, _1_distm, mask) after constructing match/mismatch; diagnostic-only"
        },
    );
    assert!(
        classification == "PAIRHMM_MISMATCH_PRIOR_DIVERGENCE_CAUSAL"
            || classification == "PAIRHMM_MISMATCH_PRIOR_DIVERGENCE_NOT_CAUSAL"
            || classification == "PAIRHMM_MISMATCH_PRIOR_RESULT_MATCH"
    );
}
