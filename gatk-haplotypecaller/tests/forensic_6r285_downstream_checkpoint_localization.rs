//! 6R.285: first downstream checkpoint after the closed PairHMM matrix.
//! PRODUCTION CHANGE: NONE.
//!
//! Java values are from one GATK 4.4.0.0 `HaplotypeCallerEngine` pass
//! (`AVX_LOGLESS_CACHING`, native threads 1, double precision off) at
//! `20:29455649`. Continuous PLs are the doubles inside `GLsToPLs`
//! before `Math.round`.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r285_downstream_checkpoint_localization -- --nocapture --test-threads=1
//! ```

use gatk_haplotypecaller::{
    biallelic_genotype_log10_likelihoods_gatk, logless_pairhmm_likelihood, ReadLikelihoodRow,
};
use std::collections::HashMap;

const FLOOR: f64 = -4.5;
const JOINT_GL: f64 = -351.91571812977724676;
const BASELINE_PL_BITS: u64 = 0x40ab7e507a112af4;
const INPUTS: &str = include_str!("6r284_frozen_inputs.tsv");
const CAPTURE: &str = include_str!("6r284_live_gkl_capture.tsv");
const JAVA: &str = include_str!("6r285_java_genotyping.tsv");

struct ReadIn {
    qname: String,
    flags: u16,
    skipped: bool,
    bases: Vec<u8>,
    bq: Vec<u8>,
    iq: Vec<u8>,
    dq: Vec<u8>,
    gcp: Vec<u8>,
    l0_bits: u64,
    other_bits: u64,
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R285\t{key}\t{}", value.as_ref());
}

fn fmt_f64(x: f64) -> String {
    format!("{x:.17} bits=0x{:016x}", x.to_bits())
}

fn parse_u8s(csv: &str) -> Vec<u8> {
    if csv.is_empty() {
        return Vec::new();
    }
    csv.split(',').map(|s| s.parse::<u8>().unwrap()).collect()
}

fn bits_cell(cell: &str) -> f64 {
    let hex = cell.split("bits=").nth(1).unwrap().trim();
    f64::from_bits(u64::from_str_radix(hex.trim_start_matches("0x"), 16).unwrap())
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

fn biallelic_gls(l0: &[f64], l2: &[f64], floor: bool) -> Vec<f64> {
    let rows: Vec<ReadLikelihoodRow> = l0
        .iter()
        .zip(l2.iter())
        .enumerate()
        .map(|(i, (a, b))| {
            let (fa, fb) = if floor { floor_pair(*a, *b) } else { (*a, *b) };
            ReadLikelihoodRow {
                read_index: i,
                read_id: String::new(),
                haplotype_log10_likelihoods: vec![fa, fb],
            }
        })
        .collect();
    biallelic_genotype_log10_likelihoods_gatk(&rows, 0, 1)
}

fn emitted_hom_alt(gls: &[f64]) -> (f64, i32) {
    let best = gls.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let cont = -10.0 * (gls[2] - best);
    let pl = (cont + 0.5).floor() as i32;
    (cont, pl)
}

fn cont_vector(gls: &[f64]) -> Vec<f64> {
    let best = gls.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    gls.iter().map(|g| -10.0 * (g - best)).collect()
}

#[test]
fn forensic_6r285_downstream_checkpoint_localization() {
    let mut haps: Vec<Vec<u8>> = Vec::new();
    let mut reads: Vec<ReadIn> = Vec::new();
    for line in INPUTS.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<_> = line.split('\t').collect();
        if f[0] == "HAP" {
            haps.push(f[3].as_bytes().to_vec());
        } else if f[0] == "READ" {
            reads.push(ReadIn {
                qname: f[2].to_string(),
                flags: f[3].parse().unwrap(),
                skipped: f[4] == "1",
                bases: f[6].as_bytes().to_vec(),
                bq: parse_u8s(f[7]),
                iq: parse_u8s(f[8]),
                dq: parse_u8s(f[9]),
                gcp: parse_u8s(f[10]),
                l0_bits: f[11].parse().unwrap(),
                other_bits: f[12].parse().unwrap(),
            });
        }
    }
    let mut live: HashMap<(usize, usize), u64> = HashMap::new();
    for line in CAPTURE.lines() {
        let f: Vec<_> = line.split('\t').collect();
        if f.first() == Some(&"CELL") {
            live.insert(
                (f[1].parse().unwrap(), f[2].parse().unwrap()),
                u64::from_str_radix(f[7].trim_start_matches("0x"), 16).unwrap(),
            );
        }
    }

    let mut l0 = Vec::new();
    let mut rust_alt = Vec::new();
    let mut cf1_alt = Vec::new();
    for (ri, read) in reads.iter().enumerate() {
        let l0v = f64::from_bits(read.l0_bits);
        l0.push(l0v);
        if read.skipped {
            rust_alt.push(0.0);
            cf1_alt.push(0.0);
            continue;
        }
        let other = f64::from_bits(read.other_bits);
        let mut raw = [0.0; 5];
        let mut gkl = [0.0; 5];
        for k in 0..5 {
            raw[k] = logless_pairhmm_likelihood(
                &read.bases,
                &read.bq,
                &haps[k],
                &read.iq,
                &read.dq,
                &read.gcp,
            )
            .unwrap();
            gkl[k] = f64::from_bits(live[&(ri, k)]);
        }
        rust_alt.push(
            floor5(raw, other)
                .into_iter()
                .fold(f64::NEG_INFINITY, f64::max),
        );
        cf1_alt.push(
            floor5(gkl, other)
                .into_iter()
                .fold(f64::NEG_INFINITY, f64::max),
        );
    }
    let rust_gls = biallelic_gls(&l0, &rust_alt, true);
    let (rust_cont, rust_pl) = emitted_hom_alt(&rust_gls);
    assert!(
        (rust_gls[2] - rust_gls.iter().copied().fold(f64::NEG_INFINITY, f64::max) - JOINT_GL).abs()
            < 1e-9
    );
    assert_eq!(rust_cont.to_bits(), BASELINE_PL_BITS);
    assert_eq!(rust_pl, 3519);
    let cf1_gls = biallelic_gls(&l0, &cf1_alt, true);
    let (cf1_cont, cf1_pl) = emitted_hom_alt(&cf1_gls);
    assert_eq!(cf1_pl, 3519);

    let mut java_index = String::new();
    let mut java_alleles = String::new();
    let mut java_gl = String::new();
    let mut java_cont = String::new();
    let mut java_ipl = String::new();
    let mut java_rows: HashMap<(String, u16), (f64, f64, f64)> = HashMap::new();
    let mut java_ref_sum = 0.0;
    let mut java_alt_sum = 0.0;
    let mut java_tttg_sum = 0.0;
    for line in JAVA.lines() {
        let f: Vec<_> = line.split('\t').collect();
        if f.len() < 3 || f[0] != "6R285" {
            continue;
        }
        match f[1] {
            "frozen_hap_index" => java_index = f[2].to_string(),
            "allele_columns" => java_alleles = f[2].to_string(),
            "java_gl_vector" => java_gl = f[2..].join("\t"),
            "java_continuous_pl_vector" => java_cont = f[2..].join("\t"),
            "java_integer_pl" => java_ipl = f[2].to_string(),
            "allele_ll" => {
                let qname = f[2].to_string();
                let flags: u16 = f[3].parse().unwrap();
                let t = bits_cell(f[4]);
                let tttg = bits_cell(f[5]);
                let tgtttg = bits_cell(f[6]);
                java_ref_sum += t;
                java_tttg_sum += tttg;
                java_alt_sum += tgtttg;
                java_rows.insert((qname, flags), (t, tttg, tgtttg));
            }
            _ => {}
        }
    }
    assert_eq!(java_index, "[19, 20, 21, 22, 23]");
    assert_eq!(java_alleles, "T,TTTG,TGTTTG");
    assert_eq!(java_rows.len(), 122);
    assert_eq!(java_ipl, "570,148,3485,0,2762,3517");

    let mut rust_ref_sum = 0.0;
    let mut rust_alt_sum = 0.0;
    let mut max_ref = 0.0f64;
    let mut max_alt = 0.0f64;
    let mut j_l0 = Vec::new();
    let mut j_alt = Vec::new();
    for (ri, read) in reads.iter().enumerate() {
        if read.skipped {
            continue;
        }
        let (t, _tttg, tg) = java_rows[&(read.qname.clone(), read.flags)];
        let dref = t - l0[ri];
        let dalt = tg - rust_alt[ri];
        max_ref = max_ref.max(dref.abs());
        max_alt = max_alt.max(dalt.abs());
        rust_ref_sum += l0[ri];
        rust_alt_sum += rust_alt[ri];
        j_l0.push(t);
        j_alt.push(tg);
    }
    let cf2_gls = biallelic_gls(&j_l0, &j_alt, false);
    let (cf2_cont, cf2_pl) = emitted_hom_alt(&cf2_gls);
    let java_cont_vals: Vec<f64> = java_cont.split('\t').map(bits_cell).collect();
    let java_hom = java_cont_vals[5];
    let required = (rust_cont - java_hom).abs();
    let cf1_delta = cf1_cont - rust_cont;

    kv(
        "execution_environment",
        "broadinstitute/gatk:4.4.0.0 --platform linux/amd64; HaplotypeCallerEngine; pairHMM AVX_LOGLESS_CACHING; pairHmmNativeThreads=1; useDoublePrecision=false",
    );
    kv("frozen_hap_index", &java_index);
    kv("allele_columns", &java_alleles);
    kv("java_continuous_PL_hom_alt", fmt_f64(java_hom));
    kv("java_integer_pl", &java_ipl);
    kv("java_gl_vector", &java_gl);
    kv("java_continuous_pl_vector", &java_cont);
    kv(
        "rust_gl_vector",
        format!(
            "GL00={} GL01={} GL22={}",
            fmt_f64(rust_gls[0]),
            fmt_f64(rust_gls[1]),
            fmt_f64(rust_gls[2]),
        ),
    );
    kv(
        "rust_continuous_pl_vector",
        cont_vector(&rust_gls)
            .iter()
            .map(|v| fmt_f64(*v))
            .collect::<Vec<_>>()
            .join("\t"),
    );
    kv(
        "aggregated",
        format!(
            "java_REF={} rust_REF={} dREF={} java_TGTTTG={} rust_TGTTTG={} dALT={} java_TTTG={} max_abs_read_ref={} max_abs_read_alt={}",
            fmt_f64(java_ref_sum),
            fmt_f64(rust_ref_sum),
            fmt_f64(java_ref_sum - rust_ref_sum),
            fmt_f64(java_alt_sum),
            fmt_f64(rust_alt_sum),
            fmt_f64(java_alt_sum - rust_alt_sum),
            fmt_f64(java_tttg_sum),
            fmt_f64(max_ref),
            fmt_f64(max_alt),
        ),
    );
    kv("baseline_continuous_PL", fmt_f64(rust_cont));
    kv("cf1_continuous_PL", fmt_f64(cf1_cont));
    kv("cf1_delta", fmt_f64(cf1_delta));
    kv("cf1_integer_PL", cf1_pl.to_string());
    kv(
        "cf2_gl_vector",
        format!(
            "GL00={} GL01={} GL22={}",
            fmt_f64(cf2_gls[0]),
            fmt_f64(cf2_gls[1]),
            fmt_f64(cf2_gls[2]),
        ),
    );
    kv("cf2_continuous_PL_hom_alt", fmt_f64(cf2_cont));
    kv("cf2_integer_PL", cf2_pl.to_string());
    kv("required_movement", fmt_f64(required));
    let allele_material = max_ref > 0.1 || max_alt > 0.1;
    let cf2_closes = (cf2_cont - java_hom).abs() < 1e-3;
    let classification = if allele_material && cf2_closes {
        "DOWNSTREAM_PER_READ_ALLELE_LIKELIHOOD_MATERIAL"
    } else if allele_material {
        "DOWNSTREAM_GENOTYPE_LIKELIHOOD_CALCULATOR_MATERIAL"
    } else {
        "DOWNSTREAM_PL_NORMALIZATION_MATERIAL"
    };
    kv("material_allele_rows", allele_material.to_string());
    kv("cf2_matches_java_continuous", cf2_closes.to_string());
    kv("classification", classification);
    kv("production_change", "NONE");
    assert_eq!(
        classification,
        "DOWNSTREAM_PER_READ_ALLELE_LIKELIHOOD_MATERIAL"
    );
    assert!(required > 1.0);
    assert_eq!(cf2_pl, 3517);
}
