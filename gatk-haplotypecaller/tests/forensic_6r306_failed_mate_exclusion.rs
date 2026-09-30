//! 6R.306: diagnostic exclusion of mate-contig failures from PairHMM evidence,
//! composed with the 6R.294 trim span and the production Jacobian combine.
//! The exclusion is the existing mate-contig predicate. No likelihood, haplotype,
//! or genotype values are injected. The call restores the pre-6R.312 name
//! sort, which is the order that left the 2 ULP / 1 ULP gap.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r306_failed_mate_exclusion -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::assembly_region_finalize::set_forensic_6r308_java_read_coordinate_order;
use gatk_haplotypecaller::hc_genotyping_engine::diploid_genotype_log10_likelihoods_from_allele_rows;
use gatk_haplotypecaller::read_assembly_filter::{passes_assembly_read, AssemblyReadFilterConfig};
use gatk_haplotypecaller::ReadLikelihoodRow;
use gatk_haplotypecaller::{
    begin_hap_list_observe, begin_likelihood_pipeline_observe, begin_trim_calc_observe,
    call_disposition, flatten_assembly_regions, forensic_6r306_excluded_count,
    forensic_6r306_mate_fail_noted, set_forensic_6r294_padded_span,
    set_forensic_6r306_exclude_failed_mate, take_colocated_merge_numerics, take_hap_list_snaps,
    take_hap_list_trim_span, take_likelihood_pipeline_cells, take_trim_calc_snap,
    traverse_assembly_region_walker, try_emit_call_region_variants, AssemblyRegionCallDisposition,
    CallRegionArgs, HapListSnap, HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
    DEFAULT_STAND_EMIT_CONFIDENCE,
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
const MATE_QNAME: &str = "HISEQ1:11:H8GV6ADXX:2:1103:14252:55237";
const MATE_FLAGS: u16 = 97;
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
    println!("6R306\t{key}\t{}", value.as_ref());
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
            hap_len.push(len);
            hap_ref.push(p[4] == "ref=true");
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
fn forensic_6r306_failed_mate_exclusion() {
    let text = std::fs::read_to_string(repo_root().join(JAVA_TSV)).unwrap();
    assert_eq!(field_raw(&text, "pairhmm_impl"), "LOGLESS_CACHING");
    let java = load_java(&text);
    let java_gl = field_cells(&text, "java_gl_vector");
    let java_pl = field_cells(&text, "java_continuous_pl_vector");
    assert_eq!(java.order.len(), 122);
    for i in 0..6 {
        assert_eq!(fmt_f(java_gl[i]), JAVA_GL[i]);
        assert_eq!(fmt_f(java_pl[i]), JAVA_PL[i]);
    }

    clear_other_forensics();
    set_forensic_6r294_padded_span(Some(JAVA_SPAN.0), Some(JAVA_SPAN.1));
    set_forensic_6r306_exclude_failed_mate(true);
    set_forensic_6r308_java_read_coordinate_order(false);
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
    let filter_cfg = AssemblyReadFilterConfig::gatk_defaults();
    let active_before = region.reads.len();
    let passing_filter = region
        .reads
        .iter()
        .filter(|r| passes_assembly_read(r, &filter_cfg))
        .count();
    let outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let excluded = forensic_6r306_excluded_count();
    let mate_fail_noted = forensic_6r306_mate_fail_noted();
    let trim = take_trim_calc_snap().expect("trim snap");
    let applied = take_hap_list_trim_span().expect("applied span");
    let snaps = take_hap_list_snaps();
    let cells = take_likelihood_pipeline_cells();
    let read_snaps = gatk_haplotypecaller::likelihood_engine::take_forensic_6r289_snaps();
    set_forensic_6r294_padded_span(None, None);
    set_forensic_6r306_exclude_failed_mate(false);
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

    let mate_key = (MATE_QNAME.to_string(), MATE_FLAGS);
    let mate_kernel_cells = cells
        .iter()
        .filter(|c| c.qname == MATE_QNAME && c.flags == MATE_FLAGS)
        .count();
    let pairhmm_keys: BTreeSet<(String, u16)> = cells
        .iter()
        .filter(|c| c.stage == "post_kernel")
        .map(|c| (c.qname.clone(), c.flags))
        .collect();

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
    let java_only: Vec<_> = java
        .order
        .iter()
        .filter(|k| !rust_by_key.contains_key(*k))
        .cloned()
        .collect();
    let rust_only: Vec<_> = rust_by_key
        .keys()
        .filter(|k| !java.alleles.contains_key(*k))
        .cloned()
        .collect();
    let shared = java.order.len() - java_only.len();
    let membership_ok = columns_named
        && java.order.len() == 122
        && rust_by_key.len() == 122
        && shared == 122
        && java_only.is_empty()
        && rust_only.is_empty()
        && !rust_by_key.contains_key(&mate_key)
        && mate_kernel_cells == 0
        && excluded >= 1;

    let mut compared = 0usize;
    let mut differing = 0usize;
    let mut max_abs = 0.0_f64;
    let mut max_ulp = 0u64;
    let mut first: Option<String> = None;
    let allele_names = ["T", "TTTG", "TGTTTG"];
    if membership_ok {
        for (i, key) in java.order.iter().enumerate() {
            let rust = rust_by_key[key];
            let java_row = java.alleles[key];
            for a in 0..3 {
                compared += 1;
                let abs = (rust[a] - java_row[a]).abs();
                let ulp = ulp_distance(rust[a], java_row[a]);
                if abs > max_abs {
                    max_abs = abs;
                }
                if ulp != u64::MAX && ulp > max_ulp {
                    max_ulp = ulp;
                }
                if rust[a].to_bits() != java_row[a].to_bits() {
                    differing += 1;
                    if first.is_none() {
                        first = Some(format!(
                            "evidence_index={i} qname={} flags={} allele={} rust={} java={} ulp={ulp}",
                            key.0,
                            key.1,
                            allele_names[a],
                            fmt_f(rust[a]),
                            fmt_f(java_row[a]),
                        ));
                    }
                }
            }
        }
    }
    let matrix_ok = membership_ok && compared == 366 && differing == 0 && max_ulp == 0;

    let mut java_order_rows = Vec::new();
    let mut rust_order_rows = Vec::new();
    if matrix_ok {
        for key in &java.order {
            java_order_rows.push(ReadLikelihoodRow {
                read_index: java_order_rows.len(),
                read_id: String::new(),
                haplotype_log10_likelihoods: rust_by_key[key].to_vec(),
            });
        }
        for i in 0..site.ad_row_qname.len() {
            let key = (site.ad_row_qname[i].clone(), site.ad_row_flags[i]);
            rust_order_rows.push(ReadLikelihoodRow {
                read_index: i,
                read_id: String::new(),
                haplotype_log10_likelihoods: rust_by_key[&key].to_vec(),
            });
        }
    }
    let replay_java_order = if java_order_rows.is_empty() {
        Vec::new()
    } else {
        diploid_genotype_log10_likelihoods_from_allele_rows(&java_order_rows, 3)
    };
    let replay_rust_order = if rust_order_rows.is_empty() {
        Vec::new()
    } else {
        diploid_genotype_log10_likelihoods_from_allele_rows(&rust_order_rows, 3)
    };
    let java_order_matches = replay_java_order.len() == java_gl.len()
        && replay_java_order
            .iter()
            .zip(java_gl.iter())
            .all(|(r, j)| r.to_bits() == j.to_bits());
    let first_order_mismatch = site
        .ad_row_qname
        .iter()
        .zip(site.ad_row_flags.iter())
        .enumerate()
        .find(|(i, (q, f))| {
            let rust_q = q.as_str();
            let rust_f = **f;
            java.order
                .get(*i)
                .map(|(jq, jf)| jq.as_str() != rust_q || *jf != rust_f)
                .unwrap_or(true)
        })
        .map(|(i, (q, f))| {
            let java_at = java
                .order
                .get(i)
                .map(|(jq, jf)| format!("{jq} flags={jf}"))
                .unwrap_or_else(|| "NONE".into());
            format!("index={i} rust={q} flags={f} java={java_at}")
        });

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
    let gl_match = matrix_ok
        && rust_gl.len() == 6
        && rust_gl
            .iter()
            .zip(java_gl.iter())
            .all(|(r, j)| r.to_bits() == j.to_bits());
    let pl_match = matrix_ok
        && rust_pl.len() == 6
        && rust_pl
            .iter()
            .zip(java_pl.iter())
            .all(|(r, j)| r.to_bits() == j.to_bits());
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

    let live_matches_rust_order = replay_rust_order.len() == rust_gl.len()
        && replay_rust_order
            .iter()
            .zip(rust_gl.iter())
            .all(|(r, j)| r.to_bits() == j.to_bits());
    let classification = if !geometry_ok || !columns_named || excluded == 0 {
        "FAILED_MATE_EXCLUSION_SETUP_FAILURE"
    } else if !membership_ok {
        "FAILED_MATE_EXCLUSION_MEMBERSHIP_DIVERGENCE"
    } else if matrix_ok && gl_match && pl_match && integer_match {
        "FAILED_MATE_EXCLUSION_REPRODUCES_JAVA"
    } else if matrix_ok {
        "FAILED_MATE_EXCLUSION_PARTIAL"
    } else {
        "FAILED_MATE_EXCLUSION_NOT_CAUSAL"
    };

    let diag_key = (ALLELE_QNAME.to_string(), ALLELE_FLAGS);
    kv("active_region_reads", active_before.to_string());
    kv("passing_assembly_filter", passing_filter.to_string());
    kv("mate_fail_noted", mate_fail_noted.to_string());
    kv("failed_mate_excluded", excluded.to_string());
    kv("post_kernel_read_keys", pairhmm_keys.len().to_string());
    kv("genotyping_reads", rust_by_key.len().to_string());
    kv("n_reads_field", site.n_reads.to_string());
    kv("java_genotyping_reads", "122");
    kv("shared_reads", shared.to_string());
    kv("java_only", java_only.len().to_string());
    kv("rust_only", rust_only.len().to_string());
    kv(
        "java_only_keys",
        java_only
            .iter()
            .map(|(q, f)| format!("{q} flags={f}"))
            .collect::<Vec<_>>()
            .join(";"),
    );
    kv(
        "rust_only_keys",
        rust_only
            .iter()
            .map(|(q, f)| format!("{q} flags={f}"))
            .collect::<Vec<_>>()
            .join(";"),
    );
    kv(
        "failed_mate_in_genotyping",
        rust_by_key.contains_key(&mate_key).to_string(),
    );
    kv("failed_mate_kernel_cells", mate_kernel_cells.to_string());
    kv(
        "trim_span",
        format!("{}-{}", applied.trim_start, applied.trim_end),
    );
    kv(
        "production_trim_span",
        format!("{}-{}", trim.padded_start, trim.padded_end),
    );
    kv("diagnostic_read", &read_repr);
    kv(
        "reference_haplotype",
        format!("len={} cigar={}", ref_h.len, ref_h.cigar),
    );
    kv(
        "tgtttg_lengths",
        tgtttg_lens
            .iter()
            .map(|n| n.to_string())
            .collect::<Vec<_>>()
            .join(","),
    );
    kv("compared_cells", compared.to_string());
    kv("differing_cells", differing.to_string());
    kv("max_abs_diff", fmt_f(max_abs));
    kv("max_ulp", max_ulp.to_string());
    kv(
        "first_differing_cell",
        first.clone().unwrap_or_else(|| "NONE".into()),
    );
    if let (Some(j), Some(r)) = (java.alleles.get(&diag_key), rust_by_key.get(&diag_key)) {
        kv(
            "diagnostic_allele_read",
            format!(
                "{ALLELE_QNAME} flags={ALLELE_FLAGS} java={} rust={}",
                j.iter().map(|v| fmt_f(*v)).collect::<Vec<_>>().join(","),
                r.iter().map(|v| fmt_f(*v)).collect::<Vec<_>>().join(",")
            ),
        );
    }
    kv(
        "java_order_replay_matches_java_gl",
        java_order_matches.to_string(),
    );
    kv(
        "live_matches_rust_row_order_replay",
        live_matches_rust_order.to_string(),
    );
    kv(
        "first_read_order_mismatch",
        first_order_mismatch.unwrap_or_else(|| "NONE".into()),
    );
    kv("final_gl", fmt_vec(&rust_gl));
    kv("gl_rust_minus_java", fmt_vec(&gl_delta));
    kv("final_pl", fmt_vec(&rust_pl));
    kv("pl_rust_minus_java", fmt_vec(&pl_delta));
    kv("integer_pl", &integer);
    kv("emitted_subset_pl", &emitted_pl);
    kv("baseline_emitted_pl", "570,0,3518");
    kv("6r304_emitted_pl", "570,0,3517");
    kv("classification", classification);
    kv("production_change", "NONE");

    assert!(geometry_ok, "trim geometry");
    assert!(membership_ok, "evidence membership: {classification}");
    assert!(matrix_ok, "allele matrix: {}", first.unwrap_or_default());
    assert!(java_order_matches, "Java-order replay of the Rust matrix");
    assert!(live_matches_rust_order, "live GL is the Rust row-order sum");
    assert!(!gl_match, "live GL unexpectedly bit-identical");
    assert!(!pl_match, "live continuous PL unexpectedly bit-identical");
    assert!(integer_match, "emitted {integer} subset {emitted_pl}");
    assert_eq!(ulp_distance(rust_gl[3], java_gl[3]), 2);
    assert!(rust_gl[3] > java_gl[3]);
    assert_eq!(ulp_distance(rust_gl[5], java_gl[5]), 1);
    assert!(rust_gl[5] < java_gl[5]);
    for i in [0, 1, 2, 4] {
        assert_eq!(rust_gl[i].to_bits(), java_gl[i].to_bits());
    }
    assert_eq!(classification, "FAILED_MATE_EXCLUSION_PARTIAL");
}
