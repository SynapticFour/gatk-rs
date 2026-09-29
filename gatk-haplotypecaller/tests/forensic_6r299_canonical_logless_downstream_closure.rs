//! 6R.299: canonical LOGLESS Java downstream versus span-forced Rust.
//! The historical AVX file is not read.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r299_canonical_logless_downstream_closure -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::hc_genotyping_engine::diploid_genotype_log10_likelihoods_from_allele_rows;
use gatk_haplotypecaller::{
    begin_likelihood_pipeline_observe, call_disposition, flatten_assembly_regions,
    set_forensic_6r294_padded_span, take_colocated_merge_numerics, take_likelihood_pipeline_cells,
    traverse_assembly_region_walker, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, ReadFilterParams, ReadLikelihoodRow, WalkerTraversalConfig,
};
use std::collections::HashMap;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_TSV: &str = "gatk-haplotypecaller/tests/6r297_java_hap_ll.tsv";
const LOCUS: u64 = 29_455_649;
const ALT: &str = "TGTTTG";
const JAVA_SPAN: (u64, u64) = (29_455_560, 29_455_728);
const FLOOR_DELTA: f64 = -4.5;
const POOLS: [(usize, usize); 3] = [(0, 14), (14, 19), (19, 24)];
const TARGET_QNAME: &str = "HWI-D00360:5:H814YADXX:2:2207:19511:63503";
const TARGET_FLAGS: u16 = 163;
const LOGLESS_T_BITS: u64 = 0xc0029ba6bfb37500;
/// Live span-forced hom-alt after 6R.303's Jacobian heterozygote combine.
const RUST_CITED_HOM: &str = "3517.46534871093808761";
const AVX_CONTINUOUS_HOM: &str = "3517.46536813194416027";
const AVX_INTEGER: &str = "570,148,3485,0,2762,3517";

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R299\t{key}\t{}", value.as_ref());
}

fn fmt_f(x: f64) -> String {
    format!("{x:.17}")
}

fn fmt_vec(xs: &[f64]) -> String {
    xs.iter().map(|x| fmt_f(*x)).collect::<Vec<_>>().join(",")
}

fn cell_bits(cell: &str) -> f64 {
    let hex = cell.split("bits=").nth(1).unwrap().trim();
    f64::from_bits(u64::from_str_radix(hex.trim_start_matches("0x"), 16).unwrap())
}

fn field_cells(text: &str, key: &str) -> Vec<f64> {
    for line in text.lines() {
        let mut p = line.split('\t');
        if p.next() == Some("6R297") && p.next() == Some(key) {
            return p.map(cell_bits).collect();
        }
    }
    panic!("missing {key}");
}

fn field_raw(text: &str, key: &str) -> String {
    for line in text.lines() {
        let mut p = line.split('\t');
        if p.next() == Some("6R297") && p.next() == Some(key) {
            return p.collect::<Vec<_>>().join("\t");
        }
    }
    panic!("missing {key}");
}

struct HapRow {
    evidence_index: usize,
    values: Vec<f64>,
}

fn load(
    text: &str,
) -> (
    Vec<(String, u16)>,
    HashMap<(String, u16), HapRow>,
    HashMap<(String, u16), [f64; 3]>,
) {
    let mut haps = HashMap::new();
    let mut alleles = HashMap::new();
    let mut order = Vec::new();
    for line in text.lines() {
        let p: Vec<_> = line.split('\t').collect();
        if p.len() < 5 || p[0] != "6R297" {
            continue;
        }
        if p[1] == "hap_ll" {
            let qname = p[2].to_string();
            let flags: u16 = p[3].parse().unwrap();
            let evidence_index: usize = p[4].parse().unwrap();
            haps.insert(
                (qname, flags),
                HapRow {
                    evidence_index,
                    values: p[5..].iter().copied().map(cell_bits).collect(),
                },
            );
        } else if p[1] == "allele_ll" {
            let qname = p[2].to_string();
            let flags: u16 = p[3].parse().unwrap();
            let vals: Vec<f64> = p[4..].iter().copied().map(cell_bits).collect();
            order.push((qname.clone(), flags));
            alleles.insert((qname, flags), [vals[0], vals[1], vals[2]]);
        }
    }
    (order, haps, alleles)
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

fn continuous_pl(gls: &[f64]) -> Vec<f64> {
    let best = gls.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    gls.iter().map(|g| -10.0 * (g - best)).collect()
}

#[test]
fn forensic_6r299_canonical_logless_downstream_closure() {
    let text = std::fs::read_to_string(repo_root().join(JAVA_TSV)).unwrap();
    assert_eq!(field_raw(&text, "pairhmm_impl"), "LOGLESS_CACHING");
    assert_eq!(field_raw(&text, "allele_columns"), "T,TTTG,TGTTTG");
    assert!(!text.contains("6r285_java_genotyping"));
    let (order, haps, alleles) = load(&text);
    assert_eq!(order.len(), 122);
    let mut java_rows = Vec::new();
    let mut marginal_mismatch = 0usize;
    for (i, key) in order.iter().enumerate() {
        let hap = haps
            .get(key)
            .unwrap_or_else(|| panic!("hap missing {}", key.0));
        assert_eq!(hap.values.len(), 24);
        let pooled = pool_max(&hap.values);
        let stored = alleles[key];
        for a in 0..3 {
            if pooled[a].to_bits() != stored[a].to_bits() {
                marginal_mismatch += 1;
            }
        }
        java_rows.push(ReadLikelihoodRow {
            read_index: i,
            read_id: String::new(),
            haplotype_log10_likelihoods: stored.to_vec(),
        });
    }
    let target = &haps[&(TARGET_QNAME.to_string(), TARGET_FLAGS)];
    let target_allele = alleles[&(TARGET_QNAME.to_string(), TARGET_FLAGS)];
    assert_eq!(target_allele[0].to_bits(), LOGLESS_T_BITS);
    let t_max = target.values[..14]
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let t_argmax: Vec<usize> = (0..14)
        .filter(|&h| target.values[h].to_bits() == t_max.to_bits())
        .collect();
    assert_eq!(t_max.to_bits(), LOGLESS_T_BITS);

    let java_gl = field_cells(&text, "java_gl_vector");
    let java_pl = field_cells(&text, "java_continuous_pl_vector");
    let java_int = field_raw(&text, "java_integer_pl");
    assert_eq!(java_gl.len(), 6);
    assert_eq!(java_pl.len(), 6);
    let replay_gl = diploid_genotype_log10_likelihoods_from_allele_rows(&java_rows, 3);
    let replay_pl = continuous_pl(&replay_gl);

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
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing chr20_tiny BAM/REF");
        return;
    }
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
    let numerics = take_colocated_merge_numerics();
    let site = numerics
        .iter()
        .find(|n| n.loc == LOCUS && n.alts.iter().any(|a| a == ALT))
        .expect("site");
    let mut raw: HashMap<(String, u16), Vec<f64>> = HashMap::new();
    for cell in cells.iter().filter(|c| c.stage == "post_kernel") {
        let row = raw
            .entry((cell.qname.clone(), cell.flags))
            .or_insert_with(|| vec![f64::NAN; 24]);
        if cell.hap_index < 24 && row[cell.hap_index].is_nan() {
            row[cell.hap_index] = cell.log10_likelihood;
        }
    }
    let mut rust_rows = Vec::new();
    let mut allele_mismatch = 0usize;
    for (i, key) in order.iter().enumerate() {
        let raw_row = raw
            .get(key)
            .unwrap_or_else(|| panic!("rust missing {}", key.0));
        let pooled = pool_max(&floor_row(raw_row));
        let stored = alleles[key];
        for a in 0..3 {
            if pooled[a].to_bits() != stored[a].to_bits() {
                allele_mismatch += 1;
            }
        }
        rust_rows.push(ReadLikelihoodRow {
            read_index: i,
            read_id: String::new(),
            haplotype_log10_likelihoods: pooled.to_vec(),
        });
    }
    let rust_from_alleles = diploid_genotype_log10_likelihoods_from_allele_rows(&rust_rows, 3);
    let rust_gl = site.merged_gls.clone();
    let rust_pl = continuous_pl(&rust_gl);
    let rust_int = site
        .merged_pl
        .iter()
        .map(|p| p.to_string())
        .collect::<Vec<_>>()
        .join(",");
    let gl_delta: Vec<f64> = rust_gl
        .iter()
        .zip(java_gl.iter())
        .map(|(r, j)| r - j)
        .collect();
    let pl_delta: Vec<f64> = rust_pl
        .iter()
        .zip(java_pl.iter())
        .map(|(r, j)| r - j)
        .collect();
    let gl_identical = gl_delta.iter().all(|d| *d == 0.0);
    let pl_identical = pl_delta.iter().all(|d| *d == 0.0);
    let (checkpoint, classification) = if marginal_mismatch != 0 || allele_mismatch != 0 {
        (
            "allele_marginalization",
            "CANONICAL_LOGLESS_DOWNSTREAM_DIVERGENCE",
        )
    } else if !gl_identical {
        (
            "genotype_likelihood",
            "CANONICAL_LOGLESS_DOWNSTREAM_DIVERGENCE",
        )
    } else if !pl_identical {
        ("pl_conversion", "CANONICAL_LOGLESS_DOWNSTREAM_DIVERGENCE")
    } else {
        ("NONE", "CANONICAL_LOGLESS_JAVA_RUST_PL_IDENTICAL")
    };

    assert_eq!(marginal_mismatch, 0);
    assert_eq!(allele_mismatch, 0);
    assert_eq!(checkpoint, "genotype_likelihood");
    assert_eq!(classification, "CANONICAL_LOGLESS_DOWNSTREAM_DIVERGENCE");
    assert_eq!(fmt_f(rust_pl[5]), RUST_CITED_HOM);
    assert_eq!(java_int, "570,148,3485,0,2762,3517");
    assert_eq!(rust_int, java_int);
    assert_eq!(replay_gl[0].to_bits(), java_gl[0].to_bits());
    assert_eq!(replay_gl[1].to_bits(), java_gl[1].to_bits());

    kv(
        "allele_source",
        "6r297_java_hap_ll.tsv allele_ll from Capture297.dumpEvent",
    );
    kv(
        "allele_dimensions",
        format!("122x3 mismatch_vs_pool_max={marginal_mismatch}"),
    );
    kv("diagnostic_T", fmt_f(target_allele[0]));
    kv(
        "diagnostic_T_argmax",
        t_argmax
            .iter()
            .map(|h| h.to_string())
            .collect::<Vec<_>>()
            .join(","),
    );
    kv("diagnostic_evidence", target.evidence_index.to_string());
    kv("java_gl", fmt_vec(&java_gl));
    kv("rust_live_gl", fmt_vec(&rust_gl));
    kv("rust_allele_replay_gl", fmt_vec(&rust_from_alleles));
    kv("java_allele_replay_gl", fmt_vec(&replay_gl));
    kv("gl_deltas_live_minus_java", fmt_vec(&gl_delta));
    kv("java_continuous_pl", fmt_vec(&java_pl));
    kv("rust_continuous_pl", fmt_vec(&rust_pl));
    kv("java_allele_replay_pl", fmt_vec(&replay_pl));
    kv("pl_deltas_live_minus_java", fmt_vec(&pl_delta));
    kv("rust_cited_hom", RUST_CITED_HOM);
    kv("rust_live_hom", fmt_f(rust_pl[5]));
    kv("java_integer_pl", &java_int);
    kv("rust_integer_pl", &rust_int);
    kv("avx_integer_pl_reference_only", AVX_INTEGER);
    kv("avx_continuous_hom_reference_only", AVX_CONTINUOUS_HOM);
    kv("allele_cell_mismatches", allele_mismatch.to_string());
    kv("first_divergence", checkpoint);
    kv("classification", classification);
    kv("production_change", "NONE");
}
