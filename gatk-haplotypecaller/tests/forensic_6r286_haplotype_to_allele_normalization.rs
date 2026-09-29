//! 6R.286: haplotype likelihoods to allele likelihoods at 20:29455649.
//! PRODUCTION CHANGE: NONE.
//!
//! The live GATK 4.4 matrix is already past `normalizeLikelihoods`.
//! This gate replays `createAlleleMapper` grouping and the marginalize
//! max, then feeds those allele columns into the existing genotype calculator.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r286_haplotype_to_allele_normalization -- --nocapture --test-threads=1
//! ```

use gatk_haplotypecaller::{
    biallelic_genotype_log10_likelihoods_gatk, logless_pairhmm_likelihood, ReadLikelihoodRow,
};
use std::collections::HashMap;

const FLOOR: f64 = -4.5;
const JOINT_GL: f64 = -351.91571812977724676;
const BASELINE_PL_BITS: u64 = 0x40ab7e507a112af4;
const JAVA_HOM_BITS: u64 = 0x40ab7aee44bb569a;
const INPUTS: &str = include_str!("6r284_frozen_inputs.tsv");
const JAVA285: &str = include_str!("6r285_java_genotyping.tsv");
const JAVA286: &str = include_str!("6r286_haplotype_allele_map.tsv");

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

struct HapRow {
    values: Vec<f64>,
    best: f64,
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R286\t{key}\t{}", value.as_ref());
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

fn hex_f64(hex: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(hex, 16).unwrap())
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

fn pool_max(values: &[f64], idx: &[usize]) -> f64 {
    idx.iter()
        .map(|&i| values[i])
        .fold(f64::NEG_INFINITY, f64::max)
}

#[test]
fn forensic_6r286_haplotype_to_allele_normalization() {
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

    let mut java_rows: HashMap<(String, u16), (f64, f64, f64)> = HashMap::new();
    let mut java_hom = 0.0;
    for line in JAVA285.lines() {
        let f: Vec<_> = line.split('\t').collect();
        if f.len() < 3 || f[0] != "6R285" {
            continue;
        }
        match f[1] {
            "allele_ll" => {
                java_rows.insert(
                    (f[2].to_string(), f[3].parse().unwrap()),
                    (bits_cell(f[4]), bits_cell(f[5]), bits_cell(f[6])),
                );
            }
            "java_continuous_pl_vector" => {
                java_hom = bits_cell(f[7]);
            }
            _ => {}
        }
    }
    assert_eq!(java_hom.to_bits(), JAVA_HOM_BITS);
    assert_eq!(java_rows.len(), 122);

    let mut allele_of: Vec<String> = Vec::new();
    let mut is_ref: Vec<bool> = Vec::new();
    let mut slot: Vec<i32> = Vec::new();
    let mut hap_rows: HashMap<(String, u16), HapRow> = HashMap::new();
    let mut ref_hap_index = -1i32;
    for line in JAVA286.lines() {
        let f: Vec<_> = line.split('\t').collect();
        if f.len() < 3 || f[0] != "6R286" {
            continue;
        }
        match f[1] {
            "map" => {
                allele_of.push(f[3].to_string());
                is_ref.push(f[4] == "1");
                slot.push(f[5].parse().unwrap());
            }
            "ref_hap_index" => ref_hap_index = f[2].parse().unwrap(),
            "hap_row" => {
                let qname = f[2].to_string();
                let flags: u16 = f[3].parse().unwrap();
                let n = f.len() - 2;
                assert!(n >= 24 + 2);
                let values: Vec<f64> = f[4..28].iter().copied().map(hex_f64).collect();
                let best = hex_f64(f[28]);
                hap_rows.insert((qname, flags), HapRow { values, best });
            }
            _ => {}
        }
    }
    assert_eq!(allele_of.len(), 24);
    assert_eq!(ref_hap_index, 4);
    assert!(is_ref[4]);
    for i in 0..14 {
        assert_eq!(allele_of[i], "T", "hap {i}");
    }
    for i in 14..19 {
        assert_eq!(allele_of[i], "TTTG", "hap {i}");
    }
    for i in 19..24 {
        assert_eq!(allele_of[i], "TGTTTG", "hap {i}");
        assert_eq!(slot[i], (i - 19) as i32);
    }
    let t_idx: Vec<usize> = (0..14).collect();
    let tttg_idx: Vec<usize> = (14..19).collect();
    let tg_idx: Vec<usize> = (19..24).collect();

    let mut l0 = Vec::new();
    let mut rust_alt = Vec::new();
    let mut rust_raw: Vec<[f64; 5]> = Vec::new();
    let mut rust_floor: Vec<[f64; 5]> = Vec::new();
    for read in &reads {
        l0.push(f64::from_bits(read.l0_bits));
        if read.skipped {
            rust_alt.push(0.0);
            rust_raw.push([0.0; 5]);
            rust_floor.push([0.0; 5]);
            continue;
        }
        let other = f64::from_bits(read.other_bits);
        let mut raw = [0.0; 5];
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
        }
        let floored = floor5(raw, other);
        rust_alt.push(floored.into_iter().fold(f64::NEG_INFINITY, f64::max));
        rust_raw.push(raw);
        rust_floor.push(floored);
    }
    let mut base_l0 = Vec::new();
    let mut base_alt = Vec::new();
    for (i, _) in reads.iter().enumerate() {
        base_l0.push(l0[i]);
        base_alt.push(rust_alt[i]);
    }
    let rust_gls = biallelic_gls(&base_l0, &base_alt, true);
    let (rust_cont, rust_pl) = emitted_hom_alt(&rust_gls);
    assert!(
        (rust_gls[2] - rust_gls.iter().copied().fold(f64::NEG_INFINITY, f64::max) - JOINT_GL).abs()
            < 1e-9
    );
    assert_eq!(rust_cont.to_bits(), BASELINE_PL_BITS);
    assert_eq!(rust_pl, 3519);

    let mut max_allele_err = 0.0f64;
    let mut max_below_floor = 0.0f64;
    let mut n_checked = 0usize;
    let mut max_ref_vs_single = 0.0f64;
    let mut n_ref_pool_differs = 0usize;
    let mut first: Option<usize> = None;
    let mut j_l0 = Vec::new();
    let mut j_alt = Vec::new();
    for (ri, read) in reads.iter().enumerate() {
        if read.skipped {
            continue;
        }
        let key = (read.qname.clone(), read.flags);
        let (jt, jtt, jtg) = java_rows[&key];
        let row = &hap_rows[&key];
        let best = row.values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        assert_eq!(best.to_bits(), row.best.to_bits());
        let floor = best + FLOOR;
        for &v in &row.values {
            max_below_floor = max_below_floor.max(floor - v);
        }
        let t = pool_max(&row.values, &t_idx);
        let tttg = pool_max(&row.values, &tttg_idx);
        let tg = pool_max(&row.values, &tg_idx);
        max_allele_err = max_allele_err.max((t - jt).abs());
        max_allele_err = max_allele_err.max((tttg - jtt).abs());
        max_allele_err = max_allele_err.max((tg - jtg).abs());
        let ref_only = row.values[ref_hap_index as usize];
        let dref = (t - ref_only).abs();
        max_ref_vs_single = max_ref_vs_single.max(dref);
        if dref > 1e-12 {
            n_ref_pool_differs += 1;
        }
        if first.is_none() && ((jt - l0[ri]).abs() > 0.1 || (jtg - rust_alt[ri]).abs() > 0.1) {
            first = Some(ri);
        }
        j_l0.push(t);
        j_alt.push(tg);
        n_checked += 1;
    }
    assert_eq!(n_checked, 122);
    assert!(
        max_below_floor <= 1e-12,
        "haplotype value below best-4.5 by {max_below_floor}"
    );
    assert!(
        max_allele_err <= 1e-12,
        "marginalize max mismatch {max_allele_err}"
    );

    let ri = first.expect("a material read");
    let read = &reads[ri];
    let row = &hap_rows[&(read.qname.clone(), read.flags)];
    let (jt, jtt, jtg) = java_rows[&(read.qname.clone(), read.flags)];
    let java_five: Vec<f64> = (0..5).map(|k| row.values[19 + k]).collect();
    let ref_ll = row.values[ref_hap_index as usize];

    let cf1_gls = biallelic_gls(&j_l0, &j_alt, false);
    let (cf1_cont, cf1_pl) = emitted_hom_alt(&cf1_gls);
    let required = (rust_cont - java_hom).abs();
    let toward = rust_cont - cf1_cont;
    let residual = (cf1_cont - java_hom).abs();
    let classification = if residual < 1e-3 && toward > 1.0 {
        "DOWNSTREAM_HAPLOTYPE_TO_ALLELE_NORMALIZATION_CAUSAL"
    } else {
        "DOWNSTREAM_HAPLOTYPE_TO_ALLELE_NORMALIZATION_NOT_CAUSAL"
    };

    kv(
        "call_path",
        "PairHMMLikelihoodCalculationEngine.computeReadLikelihoods:200 normalizeLikelihoods(qualToErrorProbLog10(45), true); HaplotypeCallerGenotypingEngine.assignGenotypeLikelihoods:182 createAlleleMapper; :191 marginalize",
    );
    kv(
        "operation",
        "per read, best = max over all haplotype columns; values below best-4.5 are raised to that floor; then each allele is the max of the haplotypes createAlleleMapper assigned to it",
    );
    kv(
        "mapping",
        "0-13 -> T (4 is reference); 14-18 -> TTTG; 19-23 -> TGTTTG; combine = max",
    );
    kv("first_read_index", ri.to_string());
    kv("first_read_qname", &read.qname);
    kv("first_read_flags", read.flags.to_string());
    kv(
        "first_java_haps",
        java_five
            .iter()
            .map(|v| fmt_f64(*v))
            .collect::<Vec<_>>()
            .join("\t"),
    );
    kv("first_java_ref_hap", fmt_f64(ref_ll));
    kv("first_java_best", fmt_f64(row.best));
    kv("first_java_floor", fmt_f64(row.best + FLOOR));
    kv(
        "first_java_alleles",
        format!(
            "T={} TTTG={} TGTTTG={}",
            fmt_f64(jt),
            fmt_f64(jtt),
            fmt_f64(jtg)
        ),
    );
    kv(
        "first_rust_raw",
        rust_raw[ri]
            .iter()
            .map(|v| fmt_f64(*v))
            .collect::<Vec<_>>()
            .join("\t"),
    );
    kv(
        "first_rust_floor5",
        rust_floor[ri]
            .iter()
            .map(|v| fmt_f64(*v))
            .collect::<Vec<_>>()
            .join("\t"),
    );
    kv("first_rust_l0", fmt_f64(l0[ri]));
    kv("first_rust_alt", fmt_f64(rust_alt[ri]));
    kv("max_allele_recompute_err", fmt_f64(max_allele_err));
    kv("max_below_floor", fmt_f64(max_below_floor));
    kv("n_reads_T_ne_reference_hap", n_ref_pool_differs.to_string());
    kv("max_abs_T_minus_reference_hap", fmt_f64(max_ref_vs_single));
    kv("baseline_continuous_PL", fmt_f64(rust_cont));
    kv("cf1_continuous_PL", fmt_f64(cf1_cont));
    kv("cf1_integer_PL", cf1_pl.to_string());
    kv("java_continuous_PL", fmt_f64(java_hom));
    kv("movement_toward_java", fmt_f64(toward));
    kv("required_movement", fmt_f64(required));
    kv("residual_vs_java", fmt_f64(residual));
    kv("classification", classification);
    kv("production_change", "NONE");

    assert_eq!(
        classification,
        "DOWNSTREAM_HAPLOTYPE_TO_ALLELE_NORMALIZATION_CAUSAL"
    );
    assert_eq!(cf1_pl, 3517);
    assert!(required > 1.0);
}
