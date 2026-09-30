//! 6R.307: localize why the same 122 genotyping reads accumulate in a different
//! order in Rust than in Java. Ordering is observed, not changed.
//! 6R.312 made Java's comparator the production sort. This test restores
//! the pre-6R.312 `(tid, pos, qname)` sort so the divergence it measured
//! stays on record.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r307_genotyping_order -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::assembly_region_finalize::set_forensic_6r308_java_read_coordinate_order;
use gatk_haplotypecaller::hc_genotyping_engine::diploid_genotype_log10_likelihoods_from_allele_rows;
use gatk_haplotypecaller::{
    begin_forensic_6r307_order_observe, call_disposition, flatten_assembly_regions,
    set_forensic_6r294_padded_span, set_forensic_6r306_exclude_failed_mate,
    take_colocated_merge_numerics, take_forensic_6r307_order, traverse_assembly_region_walker,
    AssemblyRegionCallDisposition, CallRegionArgs, Forensic6r307Read, HaplotypeCallerEngine,
    ReadFilterParams, ReadLikelihoodRow, WalkerTraversalConfig,
};
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_TSV: &str = "gatk-haplotypecaller/tests/6r297_java_hap_ll.tsv";
const LOCUS: u64 = 29_455_649;
const ALT: &str = "TGTTTG";
const JAVA_SPAN: (u64, u64) = (29_455_560, 29_455_728);

type Key = (String, u16);

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R307\t{key}\t{}", value.as_ref());
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

fn load_java_order(text: &str) -> Vec<Key> {
    let mut order = Vec::new();
    for line in text.lines() {
        let p: Vec<_> = line.split('\t').collect();
        if p.len() >= 4 && p[0] == "6R297" && p[1] == "allele_ll" {
            order.push((p[2].to_string(), p[3].parse().unwrap()));
        }
    }
    order
}

fn load_java_alleles(text: &str) -> HashMap<Key, [f64; 3]> {
    let mut alleles = HashMap::new();
    for line in text.lines() {
        let p: Vec<_> = line.split('\t').collect();
        if p.len() >= 7 && p[0] == "6R297" && p[1] == "allele_ll" {
            let vals: Vec<f64> = p[4..].iter().copied().map(cell_bits).collect();
            alleles.insert(
                (p[2].to_string(), p[3].parse().unwrap()),
                [vals[0], vals[1], vals[2]],
            );
        }
    }
    alleles
}

fn clear_hooks() {
    gatk_haplotypecaller::hc_genotyping_engine::set_forensic_6r287_global_hap_floor(None, None);
    gatk_haplotypecaller::hc_genotyping_engine::set_forensic_6r288_haplotype_substitute(
        None, None, 0, None,
    );
    gatk_haplotypecaller::likelihood_engine::set_forensic_6r289_target(None, 0);
    gatk_haplotypecaller::likelihood_engine::set_forensic_6r289_substitute(None, None);
    set_forensic_6r294_padded_span(Some(JAVA_SPAN.0), Some(JAVA_SPAN.1));
    set_forensic_6r306_exclude_failed_mate(true);
    set_forensic_6r308_java_read_coordinate_order(false);
}

fn arm_off() {
    set_forensic_6r294_padded_span(None, None);
    set_forensic_6r306_exclude_failed_mate(false);
}

struct Pass {
    geno: Vec<Key>,
    stages: Vec<(String, Vec<Forensic6r307Read>)>,
    gl: Vec<f64>,
    alleles: HashMap<Key, [f64; 3]>,
}

fn one_pass(root: &Path) -> Pass {
    clear_hooks();
    begin_forensic_6r307_order_observe();
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
    let _ = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let stages = take_forensic_6r307_order();
    arm_off();
    let numerics = take_colocated_merge_numerics();
    let site = numerics
        .iter()
        .find(|n| n.loc == LOCUS && n.alts.iter().any(|a| a == ALT))
        .expect("site");
    let mut geno = Vec::new();
    let mut alleles = HashMap::new();
    let col = |name: &str| -> usize {
        if site.long_ref == name {
            0
        } else {
            site.alts.iter().position(|a| a == name).unwrap() + 1
        }
    };
    let cols = [col("T"), col("TTTG"), col(ALT)];
    for i in 0..site.ad_row_qname.len() {
        let key = (site.ad_row_qname[i].clone(), site.ad_row_flags[i]);
        let row = &site.ad_row_lls[i];
        alleles.insert(key.clone(), [row[cols[0]], row[cols[1]], row[cols[2]]]);
        geno.push(key);
    }
    Pass {
        geno,
        stages,
        gl: site.merged_gls.clone(),
        alleles,
    }
}

fn keys_of(rows: &[Forensic6r307Read]) -> Vec<Key> {
    rows.iter().map(|r| (r.qname.clone(), r.flags)).collect()
}

fn relative(order: &[Key], keep: &HashSet<Key>) -> Vec<Key> {
    order
        .iter()
        .filter(|k| keep.contains(*k))
        .cloned()
        .collect()
}

fn first_mismatch(a: &[Key], b: &[Key]) -> Option<usize> {
    a.iter()
        .zip(b.iter())
        .position(|(x, y)| x != y)
        .or_else(|| {
            if a.len() != b.len() {
                Some(a.len().min(b.len()))
            } else {
                None
            }
        })
}

fn dup_count(keys: &[Key]) -> usize {
    let mut seen = HashSet::new();
    keys.iter().filter(|k| !seen.insert((*k).clone())).count()
}

/// Java `ReadCoordinateComparator` prefix: reference start, forward before reverse, name, flags.
fn java_prefix_cmp(a: &Forensic6r307Read, b: &Forensic6r307Read) -> Ordering {
    a.start_1based
        .cmp(&b.start_1based)
        .then_with(|| a.reverse.cmp(&b.reverse))
        .then_with(|| a.qname.cmp(&b.qname))
        .then_with(|| a.flags.cmp(&b.flags))
}

fn rust_clip_cmp(a: &Forensic6r307Read, b: &Forensic6r307Read) -> Ordering {
    a.start_1based
        .cmp(&b.start_1based)
        .then_with(|| a.qname.cmp(&b.qname))
}

fn rows_gl(order: &[Key], alleles: &HashMap<Key, [f64; 3]>) -> Vec<f64> {
    let rows: Vec<_> = order
        .iter()
        .enumerate()
        .map(|(i, k)| ReadLikelihoodRow {
            read_index: i,
            read_id: String::new(),
            haplotype_log10_likelihoods: alleles[k].to_vec(),
        })
        .collect();
    diploid_genotype_log10_likelihoods_from_allele_rows(&rows, 3)
}

#[test]
fn forensic_6r307_genotyping_order() {
    let root = repo_root();
    let text = std::fs::read_to_string(root.join(JAVA_TSV)).unwrap();
    let java = load_java_order(&text);
    let java_alleles = load_java_alleles(&text);
    let java_gl = field_cells(&text, "java_gl_vector");
    assert_eq!(java.len(), 122);
    let java_set: HashSet<Key> = java.iter().cloned().collect();
    assert_eq!(java_set.len(), 122);
    assert_eq!(dup_count(&java), 0);
    if !root.join(REF_REL).is_file() || !root.join(BAM_REL).is_file() {
        eprintln!("skip: missing chr20_tiny BAM/REF");
        return;
    }

    let a = one_pass(&root);
    let b = one_pass(&root);
    let deterministic = a.geno == b.geno;
    assert_eq!(a.geno.len(), 122);
    assert_eq!(dup_count(&a.geno), 0);
    let rust_set: HashSet<Key> = a.geno.iter().cloned().collect();
    assert_eq!(rust_set, java_set);

    let mut rust_at: HashMap<Key, usize> = HashMap::new();
    for (i, k) in a.geno.iter().enumerate() {
        rust_at.insert(k.clone(), i);
    }
    let mut java_at: HashMap<Key, usize> = HashMap::new();
    for (i, k) in java.iter().enumerate() {
        java_at.insert(k.clone(), i);
    }
    let first = first_mismatch(&java, &a.geno).unwrap();
    let n_pos_diff = java
        .iter()
        .zip(a.geno.iter())
        .filter(|(j, r)| j != r)
        .count();

    let stage = |name: &str| -> Vec<Key> {
        keys_of(
            &a.stages
                .iter()
                .find(|(s, _)| s == name)
                .unwrap_or_else(|| panic!("missing {name}"))
                .1,
        )
    };
    let post_trim = stage("post_trim");
    let post_filter = stage("post_filter");
    let pairhmm = stage("pairhmm_input");
    let rel_trim = relative(&post_trim, &java_set);
    let rel_filter = relative(&post_filter, &java_set);
    let rel_pairhmm = relative(&pairhmm, &java_set);
    assert_eq!(
        rel_filter.len(),
        122,
        "post_filter missing a genotyping read"
    );
    assert_eq!(rel_pairhmm.len(), 122, "pairhmm missing a genotyping read");

    let pairhmm_rows = &a
        .stages
        .iter()
        .find(|(s, _)| s == "pairhmm_input")
        .unwrap()
        .1;
    let mut by_key: HashMap<Key, Forensic6r307Read> = HashMap::new();
    for row in pairhmm_rows {
        by_key.insert((row.qname.clone(), row.flags), row.clone());
    }
    let mut java_sorted: Vec<Key> = a.geno.clone();
    java_sorted.sort_by(|x, y| java_prefix_cmp(&by_key[x], &by_key[y]));
    let mut rust_sorted: Vec<Key> = a.geno.clone();
    rust_sorted.sort_by(|x, y| rust_clip_cmp(&by_key[x], &by_key[y]));
    let java_cmp_matches = java_sorted == java;
    let rust_cmp_matches = rust_sorted == a.geno;
    let filter_matches_java = rel_filter == java;
    let pairhmm_matches_java = rel_pairhmm == java;
    let geno_matches_pairhmm = a.geno == rel_pairhmm;
    let pairhmm_changes_filter = rel_pairhmm != rel_filter;

    let earliest = if rel_trim != java {
        "post_trim"
    } else if rel_filter != java {
        "post_filter"
    } else if rel_pairhmm != java {
        "pairhmm_input"
    } else if a.geno != java {
        "genotype_accumulation"
    } else {
        "NONE"
    };

    let j0 = &java[first];
    let r0 = &a.geno[first];
    let jr = &by_key[j0];
    let rr = &by_key[r0];

    for i in 0..20 {
        kv(
            "pos",
            format!(
                "{i}\t{}\t{}\t{}\t{}",
                java[i].0, java[i].1, a.geno[i].0, a.geno[i].1
            ),
        );
    }
    let perm: Vec<String> = java
        .iter()
        .enumerate()
        .map(|(i, k)| format!("{i}:{}", rust_at[k]))
        .collect();
    kv("permutation_java_to_rust", perm.join(","));
    kv("n_positional_differences", n_pos_diff.to_string());
    kv("first_index", first.to_string());
    kv(
        "first_java",
        format!(
            "{} flags={} start={} cigar={} reverse={} java_index={} rust_index={}",
            j0.0, j0.1, jr.start_1based, jr.cigar, jr.reverse, java_at[j0], rust_at[j0]
        ),
    );
    kv(
        "first_rust",
        format!(
            "{} flags={} start={} cigar={} reverse={} java_index={} rust_index={}",
            r0.0, r0.1, rr.start_1based, rr.cigar, rr.reverse, java_at[r0], rust_at[r0]
        ),
    );
    kv("java_duplicate_identities", dup_count(&java).to_string());
    kv("rust_duplicate_identities", dup_count(&a.geno).to_string());
    kv("post_trim_n", post_trim.len().to_string());
    kv("post_filter_n", post_filter.len().to_string());
    kv("pairhmm_n", pairhmm.len().to_string());
    kv("post_trim_122_matches_java", (rel_trim == java).to_string());
    kv(
        "post_filter_122_matches_java",
        filter_matches_java.to_string(),
    );
    kv("pairhmm_122_matches_java", pairhmm_matches_java.to_string());
    kv(
        "genotyping_matches_pairhmm_122",
        geno_matches_pairhmm.to_string(),
    );
    kv(
        "pairhmm_reorders_post_filter",
        pairhmm_changes_filter.to_string(),
    );
    kv(
        "java_comparator_on_pairhmm_records_matches_java",
        java_cmp_matches.to_string(),
    );
    kv(
        "rust_pos_qname_on_pairhmm_records_matches_rust",
        rust_cmp_matches.to_string(),
    );
    kv("earliest_122_order_divergence", earliest);
    kv("deterministic", deterministic.to_string());

    let java_order_gl = rows_gl(&java, &a.alleles);
    let rust_order_gl = rows_gl(&a.geno, &a.alleles);
    let java_order_bits = java_order_gl
        .iter()
        .zip(java_gl.iter())
        .all(|(x, y)| x.to_bits() == y.to_bits());
    let rust_order_is_live = rust_order_gl
        .iter()
        .zip(a.gl.iter())
        .all(|(x, y)| x.to_bits() == y.to_bits());
    kv(
        "java_order_replay_matches_java_gl",
        java_order_bits.to_string(),
    );
    kv(
        "rust_order_replay_matches_live_gl",
        rust_order_is_live.to_string(),
    );
    let _ = java_alleles;

    assert!(deterministic);
    assert_eq!(first, 0);
    assert_eq!(j0.0, "HISEQ1:11:H8GV6ADXX:1:1116:4033:65919");
    assert_eq!(j0.1, 99);
    assert_eq!(r0.0, "HISEQ1:11:H8GV6ADXX:1:1109:9994:7054");
    assert_eq!(r0.1, 83);
    assert!(java_order_bits);
    assert!(rust_order_is_live);
    assert!(geno_matches_pairhmm);
    assert_eq!(earliest, "post_trim");
    assert!(
        java_cmp_matches,
        "Java comparator did not reproduce Java order"
    );
    assert!(
        rust_cmp_matches,
        "pos,qname sort did not reproduce Rust order"
    );
    assert_eq!(jr.start_1based, rr.start_1based);
    assert!(!jr.reverse);
    assert!(rr.reverse);
    assert!(pairhmm_changes_filter);
    assert!(!filter_matches_java);
    kv("classification", "GENOTYPING_ORDER_EARLY_DIVERGENCE");
}
