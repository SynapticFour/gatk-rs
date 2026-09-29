//! 6R.295: after the Java trim span is forced, where the remaining
//! continuous-PL residual first appears.
//!
//! Java reference is the AVX capture in `6r285_java_genotyping.tsv`
//! (the capture that produced hom-alt PL 3517.46536813194416027).
//! That file stores post-normalization haplotype columns 19–23 and the
//! full 122×3 allele matrix. Production code is unchanged.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r295_residual_allele_likelihood_localization -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::hc_genotyping_engine::diploid_genotype_log10_likelihoods_from_allele_rows;
use gatk_haplotypecaller::{
    begin_likelihood_pipeline_observe, call_disposition, flatten_assembly_regions,
    set_forensic_6r294_padded_span, take_likelihood_pipeline_cells,
    traverse_assembly_region_walker, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, ReadFilterParams, ReadLikelihoodRow, WalkerTraversalConfig,
};
use std::collections::HashMap;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_TSV: &str = "gatk-haplotypecaller/tests/6r285_java_genotyping.tsv";
const LOCUS: u64 = 29_455_649;
const JAVA_SPAN: (u64, u64) = (29_455_560, 29_455_728);
const FLOOR_DELTA: f64 = -4.5;
const POOLS: [(usize, usize); 3] = [(0, 14), (14, 19), (19, 24)];

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R295\t{key}\t{}", value.as_ref());
}

fn fmt_f(x: f64) -> String {
    format!("{x:.17}")
}

fn cell_bits(cell: &str) -> f64 {
    let hex = cell.split("bits=").nth(1).unwrap().trim();
    f64::from_bits(u64::from_str_radix(hex.trim_start_matches("0x"), 16).unwrap())
}

struct JavaRow {
    qname: String,
    flags: u16,
    values: Vec<f64>,
}

fn java_rows(kind: &str) -> Vec<JavaRow> {
    let text = std::fs::read_to_string(repo_root().join(JAVA_TSV)).expect("java tsv");
    let mut out = Vec::new();
    for line in text.lines() {
        let mut p = line.split('\t');
        if p.next() != Some("6R285") || p.next() != Some(kind) {
            continue;
        }
        let qname = p.next().unwrap().to_string();
        let flags: u16 = p.next().unwrap().parse().unwrap();
        let values = p.map(cell_bits).collect();
        out.push(JavaRow {
            qname,
            flags,
            values,
        });
    }
    out
}

fn java_gls() -> Vec<f64> {
    let text = std::fs::read_to_string(repo_root().join(JAVA_TSV)).expect("java tsv");
    for line in text.lines() {
        let mut p = line.split('\t');
        if p.next() == Some("6R285") && p.next() == Some("java_gl_vector") {
            return p.map(cell_bits).collect();
        }
    }
    panic!("java gl vector");
}

fn floor_row(raw: &[f64]) -> Vec<f64> {
    let best = raw.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let floor = best + FLOOR_DELTA;
    raw.iter()
        .map(|v| {
            if v.is_finite() && *v < floor {
                floor
            } else {
                *v
            }
        })
        .collect()
}

fn pool_max(row: &[f64]) -> [f64; 3] {
    std::array::from_fn(|a| {
        let (lo, hi) = POOLS[a];
        row[lo..hi]
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max)
    })
}

struct Delta {
    max_abs: f64,
    n_nonzero: usize,
    n_cells: usize,
    sum: f64,
    where_max: String,
}

fn deltas(pairs: impl IntoIterator<Item = (f64, f64, String)>) -> Delta {
    let mut max_abs = 0.0;
    let mut n_nonzero = 0;
    let mut n_cells = 0;
    let mut sum = 0.0;
    let mut where_max = String::new();
    for (rust, java, label) in pairs {
        let d = rust - java;
        n_cells += 1;
        if d != 0.0 {
            n_nonzero += 1;
        }
        sum += d;
        if d.abs() > max_abs {
            max_abs = d.abs();
            where_max = format!(
                "{label} rust={} java={} delta={}",
                fmt_f(rust),
                fmt_f(java),
                fmt_f(d)
            );
        }
    }
    Delta {
        max_abs,
        n_nonzero,
        n_cells,
        sum,
        where_max,
    }
}

fn continuous_hom(gls: &[f64]) -> f64 {
    let best = gls.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    -10.0 * (gls[5] - best)
}

#[test]
fn forensic_6r295_residual_allele_likelihood_localization() {
    gatk_haplotypecaller::hc_genotyping_engine::set_forensic_6r287_global_hap_floor(None, None);
    gatk_haplotypecaller::hc_genotyping_engine::set_forensic_6r288_haplotype_substitute(
        None, None, 0, None,
    );
    gatk_haplotypecaller::likelihood_engine::set_forensic_6r289_target(None, 0);
    gatk_haplotypecaller::likelihood_engine::set_forensic_6r289_substitute(None, None);
    set_forensic_6r294_padded_span(Some(JAVA_SPAN.0), Some(JAVA_SPAN.1));
    begin_likelihood_pipeline_observe();

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
    let region = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= LOCUS
                && r.end.get() >= LOCUS
        })
        .expect("ActiveFull");
    let _outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    set_forensic_6r294_padded_span(None, None);
    let cells = take_likelihood_pipeline_cells();

    let mut raw: HashMap<(String, u16), Vec<f64>> = HashMap::new();
    for cell in cells.iter().filter(|c| c.stage == "post_kernel") {
        let row = raw
            .entry((cell.qname.clone(), cell.flags))
            .or_insert_with(|| vec![f64::NAN; 24]);
        if cell.hap_index < 24 && row[cell.hap_index].is_nan() {
            row[cell.hap_index] = cell.log10_likelihood;
        }
    }

    let java_haps = java_rows("hap_ll");
    let java_alleles = java_rows("allele_ll");
    assert_eq!(java_alleles.len(), 122);
    assert!(java_haps.iter().all(|r| r.values.len() == 5));

    let mut raw_pairs = Vec::new();
    let mut floor_pairs = Vec::new();
    let mut allele_pairs = Vec::new();
    let mut n_lifted = 0usize;
    let mut rust_allele_rows = Vec::new();
    let mut java_allele_rows = Vec::new();
    let mut missing = 0usize;
    for (i, j) in java_alleles.iter().enumerate() {
        let Some(raw_row) = raw.get(&(j.qname.clone(), j.flags)) else {
            missing += 1;
            continue;
        };
        assert!(
            raw_row.iter().all(|v| v.is_finite()),
            "incomplete {}",
            j.qname
        );
        let floored = floor_row(raw_row);
        if let Some(h) = java_haps
            .iter()
            .find(|h| h.qname == j.qname && h.flags == j.flags)
        {
            for k in 0..5 {
                let label = format!("read={} flags={} hap={}", j.qname, j.flags, 19 + k);
                // Java's stored column is post-floor. A raw value below Rust's
                // own floor was raised in Java too; that gap is the floor, not
                // a kernel difference.
                if raw_row[19 + k] < floored[19 + k] {
                    n_lifted += 1;
                } else {
                    raw_pairs.push((raw_row[19 + k], h.values[k], label.clone()));
                }
                floor_pairs.push((floored[19 + k], h.values[k], label));
            }
        }
        let rust_al = pool_max(&floored);
        let mut argmax = [0usize; 3];
        for a in 0..3 {
            let (lo, hi) = POOLS[a];
            argmax[a] = (lo..hi)
                .max_by(|&i, &j| floored[i].partial_cmp(&floored[j]).unwrap())
                .unwrap();
        }
        for a in 0..3 {
            let name = ["T", "TTTG", "TGTTTG"][a];
            allele_pairs.push((
                rust_al[a],
                j.values[a],
                format!(
                    "read={} flags={} allele={name} hap={}",
                    j.qname, j.flags, argmax[a]
                ),
            ));
        }
        rust_allele_rows.push(ReadLikelihoodRow {
            read_index: i,
            read_id: String::new(),
            haplotype_log10_likelihoods: rust_al.to_vec(),
        });
        java_allele_rows.push(ReadLikelihoodRow {
            read_index: i,
            read_id: String::new(),
            haplotype_log10_likelihoods: j.values.clone(),
        });
    }
    assert_eq!(
        missing, 0,
        "java genotyping reads missing from rust PairHMM"
    );

    let raw_d = deltas(raw_pairs);
    let floor_d = deltas(floor_pairs);
    let allele_d = deltas(allele_pairs);
    let java_gl = java_gls();
    let rust_gl = diploid_genotype_log10_likelihoods_from_allele_rows(&rust_allele_rows, 3);
    let replay_gl = diploid_genotype_log10_likelihoods_from_allele_rows(&java_allele_rows, 3);
    let java_pl = continuous_hom(&java_gl);
    let rust_pl = continuous_hom(&rust_gl);
    let replay_pl = continuous_hom(&replay_gl);
    let gl_delta: Vec<f64> = rust_gl
        .iter()
        .zip(java_gl.iter())
        .map(|(r, j)| r - j)
        .collect();
    let pl_delta = rust_pl - java_pl;

    // Unfloored haplotype cells already differ, and the same deltas are the
    // allele maxima. Their sum is the continuous-PL residual. The floor does
    // not create a larger gap than that PairHMM output.
    let classification = "RESIDUAL_AT_PAIRHMM_OUTPUT";
    assert!(
        raw_d.max_abs < 1e-4,
        "unfloored PairHMM delta {}",
        raw_d.max_abs
    );
    assert!(
        raw_d.max_abs > 1e-6,
        "unfloored PairHMM delta {}",
        raw_d.max_abs
    );
    assert!((pl_delta - 0.00005013171085011).abs() < 1e-16);

    kv("geometry", "29455560-29455728");
    kv("java_geometry", "29455560-29455728");
    kv("pairhmm_lifted_cells", n_lifted.to_string());
    kv(
        "pairhmm_compared_cells",
        format!("{} unfloored of haps 19-23", raw_d.n_cells),
    );
    kv("pairhmm_max_abs_delta", fmt_f(raw_d.max_abs));
    kv("pairhmm_nonzero_cells", raw_d.n_nonzero.to_string());
    kv("pairhmm_largest", &raw_d.where_max);
    kv("post_norm_max_abs_delta", fmt_f(floor_d.max_abs));
    kv("post_norm_nonzero_cells", floor_d.n_nonzero.to_string());
    kv("post_norm_largest", &floor_d.where_max);
    kv("allele_max_abs_delta", fmt_f(allele_d.max_abs));
    kv("allele_nonzero_cells", allele_d.n_nonzero.to_string());
    kv("allele_n_cells", allele_d.n_cells.to_string());
    kv("allele_sum_delta", fmt_f(allele_d.sum));
    kv("allele_largest", &allele_d.where_max);
    kv(
        "gl_deltas",
        gl_delta
            .iter()
            .map(|d| fmt_f(*d))
            .collect::<Vec<_>>()
            .join(","),
    );
    kv("java_gl_replay_hom_pl", fmt_f(replay_pl));
    kv("rust_hom_pl_from_alleles", fmt_f(rust_pl));
    kv("java_hom_pl", fmt_f(java_pl));
    kv("final_continuous_pl_delta", fmt_f(pl_delta));
    kv("classification", classification);
    kv("production_change", "NONE");
}
