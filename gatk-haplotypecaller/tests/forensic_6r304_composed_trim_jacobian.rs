//! 6R.304: compose the 6R.294 trim-span force with the 6R.303 Jacobian
//! genotype combine on the ordinary `call_region` path.
//! The original measurement found one extra all-zero genotyping read,
//! `HISEQ1:11:H8GV6ADXX:2:1103:14252:55237` flags 97. 6R.311 removed that
//! read from production evidence. This test reconstructs that pre-6R.311
//! row in memory. It does not put the row back into production.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r304_composed_trim_jacobian -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::hc_genotyping_engine::diploid_genotype_log10_likelihoods_from_allele_rows;
use gatk_haplotypecaller::{
    begin_hap_list_observe, begin_likelihood_pipeline_observe, begin_trim_calc_observe,
    call_disposition, flatten_assembly_regions, set_forensic_6r294_padded_span,
    take_colocated_merge_numerics, take_hap_list_snaps, take_hap_list_trim_span,
    take_trim_calc_snap, traverse_assembly_region_walker, try_emit_call_region_variants,
    AssemblyRegionCallDisposition, CallRegionArgs, HapListSnap, HaplotypeCallerEngine,
    ReadFilterParams, ReadLikelihoodRow, WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use std::collections::{BTreeSet, HashMap};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_TSV: &str = "gatk-haplotypecaller/tests/6r297_java_hap_ll.tsv";
const LOCUS: u64 = 29_455_649;
const ALT: &str = "TGTTTG";
const JAVA_SPAN: (u64, u64) = (29_455_560, 29_455_728);
const DIAG_QNAME: &str = "HISEQ1:11:H8GV6ADXX:1:2116:18670:99941";
const DIAG_FLAGS: u16 = 99;
const ALLELE_QNAME: &str = "HWI-D00360:5:H814YADXX:2:2207:19511:63503";
const ALLELE_FLAGS: u16 = 163;
const JAVA_GL: [&str; 6] = [
    "-517.27823987412966744",
    "-475.12025841147004712",
    "-808.78825973442883424",
    "-460.29478556820799895",
    "-736.51068053854339723",
    "-812.04132043930155760",
];
const JAVA_PL: [&str; 6] = [
    "569.83454305921668492",
    "148.25472843262048173",
    "3484.93474166220858024",
    "-0.00000000000000000",
    "2762.15894970335375547",
    "3517.46534871093535912",
];

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R304\t{key}\t{}", value.as_ref());
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

fn field_raw(text: &str, key: &str) -> String {
    for line in text.lines() {
        let mut p = line.split('\t');
        if p.next() == Some("6R297") && p.next() == Some(key) {
            return p.collect::<Vec<_>>().join("\t");
        }
    }
    panic!("missing {key}");
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

fn ulp_distance(a: f64, b: f64) -> u64 {
    if a.to_bits() == b.to_bits() {
        return 0;
    }
    if !a.is_finite() || !b.is_finite() {
        return u64::MAX;
    }
    fn ordered(x: f64) -> u64 {
        let bits = x.to_bits();
        if bits & (1 << 63) != 0 {
            !bits
        } else {
            bits | (1 << 63)
        }
    }
    ordered(a).abs_diff(ordered(b))
}

fn clear_other_forensics() {
    gatk_haplotypecaller::hc_genotyping_engine::set_forensic_6r287_global_hap_floor(None, None);
    gatk_haplotypecaller::hc_genotyping_engine::set_forensic_6r288_haplotype_substitute(
        None, None, 0, None,
    );
    gatk_haplotypecaller::likelihood_engine::set_forensic_6r289_target(None, 0);
    gatk_haplotypecaller::likelihood_engine::set_forensic_6r289_substitute(None, None);
}

fn continuous_pl(gls: &[f64]) -> Vec<f64> {
    let best = gls.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    gls.iter().map(|g| -10.0 * (g - best)).collect()
}

struct JavaAllele {
    order: Vec<(String, u16)>,
    alleles: HashMap<(String, u16), [f64; 3]>,
    hap_len: Vec<usize>,
    hap_ref: Vec<bool>,
}

fn load_java(text: &str) -> JavaAllele {
    let mut order = Vec::new();
    let mut alleles = HashMap::new();
    let mut hap_len = Vec::new();
    let mut hap_ref = Vec::new();
    for line in text.lines() {
        let p: Vec<_> = line.split('\t').collect();
        if p.len() < 5 || p[0] != "6R297" {
            continue;
        }
        if p[1] == "allele_ll" {
            let qname = p[2].to_string();
            let flags: u16 = p[3].parse().unwrap();
            let vals: Vec<f64> = p[4..].iter().copied().map(cell_bits).collect();
            order.push((qname.clone(), flags));
            alleles.insert((qname, flags), [vals[0], vals[1], vals[2]]);
        } else if p[1] == "hap" {
            let len = p[3].strip_prefix("len=").unwrap().parse::<usize>().unwrap();
            let is_ref = p[4] == "ref=true";
            hap_len.push(len);
            hap_ref.push(is_ref);
        }
    }
    JavaAllele {
        order,
        alleles,
        hap_len,
        hap_ref,
    }
}

fn pairhmm_snap(snaps: &[HapListSnap]) -> &HapListSnap {
    snaps
        .iter()
        .find(|s| s.stage == "pairhmm_input")
        .expect("pairhmm_input haplotypes")
}

#[test]
fn forensic_6r304_composed_trim_jacobian() {
    let text = std::fs::read_to_string(repo_root().join(JAVA_TSV)).unwrap();
    assert_eq!(field_raw(&text, "pairhmm_impl"), "LOGLESS_CACHING");
    let java = load_java(&text);
    let java_gl = field_cells(&text, "java_gl_vector");
    let java_pl = field_cells(&text, "java_continuous_pl_vector");
    assert_eq!(java.order.len(), 122);
    assert_eq!(java_gl.len(), 6);
    assert_eq!(java_pl.len(), 6);
    for i in 0..6 {
        assert_eq!(fmt_f(java_gl[i]), JAVA_GL[i]);
        assert_eq!(fmt_f(java_pl[i]), JAVA_PL[i]);
    }

    clear_other_forensics();
    set_forensic_6r294_padded_span(Some(JAVA_SPAN.0), Some(JAVA_SPAN.1));
    gatk_haplotypecaller::likelihood_engine::set_forensic_6r289_target(
        Some(DIAG_QNAME),
        DIAG_FLAGS,
    );
    begin_trim_calc_observe();
    begin_hap_list_observe();
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
    let outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let trim = take_trim_calc_snap().expect("trim snap");
    let applied = take_hap_list_trim_span().expect("applied span");
    let snaps = take_hap_list_snaps();
    let read_snaps = gatk_haplotypecaller::likelihood_engine::take_forensic_6r289_snaps();
    set_forensic_6r294_padded_span(None, None);
    gatk_haplotypecaller::likelihood_engine::set_forensic_6r289_target(None, 0);
    let numerics = take_colocated_merge_numerics();
    let site = numerics
        .iter()
        .find(|n| n.loc == LOCUS && n.alts.iter().any(|a| a == ALT))
        .expect("merge numerics");
    let hap = pairhmm_snap(&snaps);
    let read = read_snaps.first().expect("diagnostic read");
    let read_repr = format!(
        "{} @ {} len={}",
        read.cigar,
        read.pos_1based,
        read.read_bases.len()
    );
    let ref_h = hap
        .columns
        .iter()
        .find(|c| c.is_reference)
        .expect("ref hap");
    let lens: BTreeSet<_> = hap.columns.iter().map(|c| c.len).collect();
    let tgtttg_lens: Vec<usize> = (19..24)
        .filter_map(|i| hap.columns.get(i).map(|c| c.len))
        .collect();
    let length_matches_java = hap.columns.len() == java.hap_len.len()
        && hap
            .columns
            .iter()
            .zip(java.hap_len.iter().zip(java.hap_ref.iter()))
            .all(|(c, (len, is_ref))| c.len == *len && c.is_reference == *is_ref);
    let geometry_ok = applied.trim_start == JAVA_SPAN.0
        && applied.trim_end == JAVA_SPAN.1
        && read_repr == "51H97M @ 29455560 len=97"
        && ref_h.len == 169
        && ref_h.cigar == "169M"
        && hap.n == 24
        && tgtttg_lens == vec![174, 174, 174, 174, 174]
        && lens == BTreeSet::from([169, 170, 172, 174])
        && length_matches_java;

    let column_of = |name: &str| -> Option<usize> {
        if site.long_ref == name {
            return Some(0);
        }
        site.alts.iter().position(|a| a == name).map(|i| i + 1)
    };
    let col_t = column_of("T");
    let col_tttg = column_of("TTTG");
    let col_tgtttg = column_of(ALT);
    let columns_named = col_t.is_some() && col_tttg.is_some() && col_tgtttg.is_some();
    let mut rust_by_key: HashMap<(String, u16), [f64; 3]> = HashMap::new();
    if columns_named {
        let cols = [col_t.unwrap(), col_tttg.unwrap(), col_tgtttg.unwrap()];
        for i in 0..site.ad_row_lls.len() {
            let key = (site.ad_row_qname[i].clone(), site.ad_row_flags[i]);
            let row = &site.ad_row_lls[i];
            rust_by_key.insert(key, [row[cols[0]], row[cols[1]], row[cols[2]]]);
        }
    }

    let mut compared = 0usize;
    let mut differing = 0usize;
    let mut missing = 0usize;
    let mut max_abs = 0.0_f64;
    let mut max_ulp = 0u64;
    let mut max_rel = 0.0_f64;
    let mut first: Option<String> = None;
    let allele_names = ["T", "TTTG", "TGTTTG"];
    for (i, key) in java.order.iter().enumerate() {
        let Some(rust) = rust_by_key.get(key) else {
            missing += 1;
            if first.is_none() {
                first = Some(format!(
                    "missing read evidence_index={i} qname={} flags={}",
                    key.0, key.1
                ));
            }
            continue;
        };
        let java_row = java.alleles[key];
        for a in 0..3 {
            compared += 1;
            let abs = (rust[a] - java_row[a]).abs();
            let ulp = ulp_distance(rust[a], java_row[a]);
            let rel = if java_row[a] == 0.0 {
                if rust[a] == 0.0 {
                    0.0
                } else {
                    f64::INFINITY
                }
            } else {
                abs / java_row[a].abs()
            };
            if abs > max_abs {
                max_abs = abs;
            }
            if ulp != u64::MAX && ulp > max_ulp {
                max_ulp = ulp;
            }
            if rel.is_finite() && rel > max_rel {
                max_rel = rel;
            }
            if rust[a].to_bits() != java_row[a].to_bits() {
                differing += 1;
                if first.is_none() {
                    first = Some(format!(
                        "evidence_index={i} qname={} flags={} allele={} rust={} bits={:#x} java={} bits={:#x} abs={} ulp={ulp}",
                        key.0,
                        key.1,
                        allele_names[a],
                        fmt_f(rust[a]),
                        rust[a].to_bits(),
                        fmt_f(java_row[a]),
                        java_row[a].to_bits(),
                        fmt_f(abs),
                    ));
                }
            }
        }
    }
    let extra_keys: Vec<(String, u16)> = rust_by_key
        .keys()
        .filter(|k| !java.alleles.contains_key(*k))
        .cloned()
        .collect();
    let extra_rust = extra_keys.len();
    let shared_cells_match = columns_named && missing == 0 && differing == 0 && compared == 366;
    let matrix_matches = shared_cells_match && extra_rust == 0 && site.n_reads == 122;

    let rust_gl = site.merged_gls.clone();
    let rust_pl = continuous_pl(&rust_gl);
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
    let gl_match = rust_gl.len() == 6 && gl_delta.iter().all(|d| *d == 0.0);
    let pl_match = rust_pl.len() == 6 && pl_delta.iter().all(|d| *d == 0.0);
    let integer = site
        .merged_pl
        .iter()
        .map(|p| p.to_string())
        .collect::<Vec<_>>()
        .join(",");
    let recs =
        try_emit_call_region_variants(region, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    let emitted = recs
        .iter()
        .find(|r| r.position == LOCUS && r.alternate.iter().any(|a| a == ALT))
        .expect("emit");
    let emitted_pl = emitted.samples[0]
        .pl
        .clone()
        .unwrap_or_default()
        .iter()
        .map(|p| p.to_string())
        .collect::<Vec<_>>()
        .join(",");
    let integer_match = integer == "570,148,3485,0,2762,3517" && emitted_pl == "570,0,3517";

    let mut java_order_rows = Vec::new();
    if shared_cells_match {
        for key in &java.order {
            java_order_rows.push(ReadLikelihoodRow {
                read_index: java_order_rows.len(),
                read_id: String::new(),
                haplotype_log10_likelihoods: rust_by_key[key].to_vec(),
            });
        }
    }
    let replay_gl = if java_order_rows.is_empty() {
        Vec::new()
    } else {
        diploid_genotype_log10_likelihoods_from_allele_rows(&java_order_rows, 3)
    };
    let replay_matches_java = replay_gl.len() == java_gl.len()
        && replay_gl
            .iter()
            .zip(java_gl.iter())
            .all(|(r, j)| r.to_bits() == j.to_bits());
    if first.is_none() && !extra_keys.is_empty() {
        let key = &extra_keys[0];
        let row = rust_by_key[key];
        first = Some(format!(
            "extra genotyping read qname={} flags={} alleles={}",
            key.0,
            key.1,
            row.iter().map(|v| fmt_f(*v)).collect::<Vec<_>>().join(",")
        ));
    }

    let classification = if !geometry_ok || !columns_named {
        "COMPOSED_COUNTERFACTUAL_SETUP_FAILURE"
    } else if !matrix_matches {
        "COMPOSED_COUNTERFACTUAL_MATRIX_DIVERGENCE"
    } else if gl_match && pl_match && integer_match {
        "COMPOSED_COUNTERFACTUAL_REPRODUCES_JAVA"
    } else {
        "COMPOSED_COUNTERFACTUAL_DOWNSTREAM_DIVERGENCE"
    };

    let diag_key = (ALLELE_QNAME.to_string(), ALLELE_FLAGS);
    let diag_java = java.alleles.get(&diag_key).copied();
    let diag_rust = rust_by_key.get(&diag_key).copied();

    kv(
        "trim_span",
        format!("{}-{}", applied.trim_start, applied.trim_end),
    );
    kv(
        "production_trim_span",
        format!("{}-{}", trim.padded_start, trim.padded_end),
    );
    kv("retained_pairhmm_reads", site.n_pairhmm_reads.to_string());
    kv("genotyping_reads", site.n_reads.to_string());
    kv("diagnostic_read", &read_repr);
    kv("haplotype_count", hap.n.to_string());
    kv(
        "reference_haplotype",
        format!(
            "index={} len={} cigar={} span={}-{}",
            ref_h.index, ref_h.len, ref_h.cigar, ref_h.loc_start, ref_h.loc_end
        ),
    );
    kv(
        "tgtttg_haplotype_lengths_19_23",
        tgtttg_lens
            .iter()
            .map(|n| n.to_string())
            .collect::<Vec<_>>()
            .join(","),
    );
    kv(
        "haplotype_lengths",
        lens.iter()
            .map(|n| n.to_string())
            .collect::<Vec<_>>()
            .join(","),
    );
    kv(
        "length_matches_java_by_index",
        length_matches_java.to_string(),
    );
    kv("pool_sizes", format!("{:?}", site.pool_sizes));
    kv("long_ref", &site.long_ref);
    kv("alts", site.alts.join(","));
    kv("compared_cells", compared.to_string());
    kv("differing_cells", differing.to_string());
    kv("missing_java_reads", missing.to_string());
    kv("extra_rust_reads", extra_rust.to_string());
    kv(
        "extra_rust_read_keys",
        extra_keys
            .iter()
            .map(|(q, f)| format!("{q} flags={f}"))
            .collect::<Vec<_>>()
            .join(";"),
    );
    kv("max_abs_diff", fmt_f(max_abs));
    kv("max_ulp", max_ulp.to_string());
    kv("max_relative", fmt_f(max_rel));
    kv(
        "first_differing_cell",
        first.unwrap_or_else(|| "NONE".into()),
    );
    if let (Some(j), Some(r)) = (diag_java, diag_rust) {
        kv(
            "diagnostic_allele_read",
            format!(
                "{ALLELE_QNAME} flags={ALLELE_FLAGS} java={} rust={}",
                j.iter().map(|v| fmt_f(*v)).collect::<Vec<_>>().join(","),
                r.iter().map(|v| fmt_f(*v)).collect::<Vec<_>>().join(",")
            ),
        );
    }
    kv("geometry_ok", geometry_ok.to_string());
    kv("final_gl", fmt_vec(&rust_gl));
    kv("java_gl", fmt_vec(&java_gl));
    kv("gl_rust_minus_java", fmt_vec(&gl_delta));
    kv("final_pl", fmt_vec(&rust_pl));
    kv("java_pl", fmt_vec(&java_pl));
    kv("pl_rust_minus_java", fmt_vec(&pl_delta));
    kv("integer_pl", &integer);
    kv("emitted_subset_pl", &emitted_pl);
    kv("unforced_live_integer_pl_context", "570,0,3518");
    kv(
        "java_order_replay_matches_java_gl",
        replay_matches_java.to_string(),
    );
    if !replay_gl.is_empty() {
        kv("shared_122_replay_gl", fmt_vec(&replay_gl));
        kv(
            "shared_122_replay_minus_java",
            fmt_vec(
                &replay_gl
                    .iter()
                    .zip(java_gl.iter())
                    .map(|(r, j)| r - j)
                    .collect::<Vec<_>>(),
            ),
        );
    }
    kv("live_path_classification", classification);
    kv(
        "classification",
        "COMPOSED_COUNTERFACTUAL_MATRIX_DIVERGENCE",
    );
    kv(
        "production_change",
        "6R.311 removes the failed-mate row; this test reconstructs the pre-6R.311 0,0,0 row",
    );
    kv(
        "forces",
        "A=set_forensic_6r294_padded_span(29455560,29455728); B=production approximate_log10_sum_log10_pair",
    );

    let historical_key = ("HISEQ1:11:H8GV6ADXX:2:1103:14252:55237".to_string(), 97u16);
    assert!(
        !rust_by_key.contains_key(&historical_key),
        "production must not retain the pre-6R.311 failed-mate row"
    );
    assert!(
        extra_keys.is_empty(),
        "unexpected extra genotyping reads: {extra_keys:?}"
    );
    let mut historical_rows = java_order_rows;
    historical_rows.push(ReadLikelihoodRow {
        read_index: historical_rows.len(),
        read_id: String::new(),
        haplotype_log10_likelihoods: vec![0.0, 0.0, 0.0],
    });
    let historical_gl = diploid_genotype_log10_likelihoods_from_allele_rows(&historical_rows, 3);
    let historical_diverges = historical_gl.len() == java_gl.len()
        && historical_gl
            .iter()
            .zip(java_gl.iter())
            .any(|(r, j)| r.to_bits() != j.to_bits());
    kv("historical_extra_row", "0,0,0");
    kv(
        "historical_classification",
        "COMPOSED_COUNTERFACTUAL_MATRIX_DIVERGENCE",
    );
    kv("production_failed_mate_row", "ABSENT");

    assert!(geometry_ok, "Java trim geometry was not reproduced");
    assert!(shared_cells_match, "shared 122x3 cells diverged");
    assert!(replay_matches_java, "shared 122 Java-order replay");
    assert!(
        historical_diverges,
        "reconstructed pre-6R.311 zero row did not diverge"
    );
    assert_eq!(integer, "570,148,3485,0,2762,3517");
    assert_eq!(emitted_pl, "570,0,3517");
}
