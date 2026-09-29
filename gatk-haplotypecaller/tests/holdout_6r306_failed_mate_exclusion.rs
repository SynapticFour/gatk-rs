//! 6R.306 holdout: independent re-run of the mate-exclusion counterfactual.
//! Skipped unless `HOLDOUT_6R306=1`. Does not run HOLDOUT_6R243.
//!
//! ```text
//! HOLDOUT_6R306=1 cargo test -p gatk-haplotypecaller --test holdout_6r306_failed_mate_exclusion -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, forensic_6r306_excluded_count,
    set_forensic_6r294_padded_span, set_forensic_6r306_exclude_failed_mate,
    take_colocated_merge_numerics, traverse_assembly_region_walker, try_emit_call_region_variants,
    AssemblyRegionCallDisposition, CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams,
    WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
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
const MATE_QNAME: &str = "HISEQ1:11:H8GV6ADXX:2:1103:14252:55237";
const MATE_FLAGS: u16 = 97;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
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

fn continuous_pl(gls: &[f64]) -> Vec<f64> {
    let best = gls.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    gls.iter().map(|g| -10.0 * (g - best)).collect()
}

#[test]
fn holdout_6r306_failed_mate_exclusion() {
    if std::env::var("HOLDOUT_6R306").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R306=1");
        return;
    }
    let text = std::fs::read_to_string(repo_root().join(JAVA_TSV)).unwrap();
    let java_gl = field_cells(&text, "java_gl_vector");
    let java_pl = field_cells(&text, "java_continuous_pl_vector");
    let mut java_keys = std::collections::BTreeSet::new();
    for line in text.lines() {
        let p: Vec<_> = line.split('\t').collect();
        if p.len() >= 4 && p[0] == "6R297" && p[1] == "allele_ll" {
            java_keys.insert((p[2].to_string(), p[3].parse::<u16>().unwrap()));
        }
    }
    assert_eq!(java_keys.len(), 122);

    gatk_haplotypecaller::hc_genotyping_engine::set_forensic_6r287_global_hap_floor(None, None);
    gatk_haplotypecaller::hc_genotyping_engine::set_forensic_6r288_haplotype_substitute(
        None, None, 0, None,
    );
    gatk_haplotypecaller::likelihood_engine::set_forensic_6r289_target(None, 0);
    gatk_haplotypecaller::likelihood_engine::set_forensic_6r289_substitute(None, None);
    set_forensic_6r294_padded_span(Some(JAVA_SPAN.0), Some(JAVA_SPAN.1));
    set_forensic_6r306_exclude_failed_mate(true);

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
    let outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let excluded = forensic_6r306_excluded_count();
    set_forensic_6r294_padded_span(None, None);
    set_forensic_6r306_exclude_failed_mate(false);

    let numerics = take_colocated_merge_numerics();
    let site = numerics
        .iter()
        .find(|n| n.loc == LOCUS && n.alts.iter().any(|a| a == ALT))
        .expect("site");
    let mut rust_keys = HashMap::new();
    for i in 0..site.ad_row_qname.len() {
        rust_keys.insert(
            (site.ad_row_qname[i].clone(), site.ad_row_flags[i]),
            site.ad_row_lls[i].clone(),
        );
    }
    assert!(!rust_keys.contains_key(&(MATE_QNAME.to_string(), MATE_FLAGS)));
    assert_eq!(rust_keys.len(), 122);
    assert!(java_keys.iter().all(|k| rust_keys.contains_key(k)));
    assert!(excluded >= 1);
    let rust_pl = continuous_pl(&site.merged_gls);
    for i in [0, 1, 2, 4] {
        assert_eq!(site.merged_gls[i].to_bits(), java_gl[i].to_bits());
    }
    assert_ne!(site.merged_gls[3].to_bits(), java_gl[3].to_bits());
    assert_ne!(site.merged_gls[5].to_bits(), java_gl[5].to_bits());
    assert_ne!(rust_pl[5].to_bits(), java_pl[5].to_bits());
    let recs =
        try_emit_call_region_variants(region, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    let emitted = recs
        .iter()
        .find(|r| r.position == LOCUS && r.alternate.iter().any(|a| a == ALT))
        .expect("emit");
    let subset = emitted.samples[0].pl.clone().unwrap_or_default();
    assert_eq!(subset.as_slice(), &[570, 0, 3517]);
    let integer = site
        .merged_pl
        .iter()
        .map(|p| p.to_string())
        .collect::<Vec<_>>()
        .join(",");
    assert_eq!(integer, "570,148,3485,0,2762,3517");
    println!("6R306\tclassification\tFAILED_MATE_EXCLUSION_PARTIAL");
    println!("6R306\tholdout\tPASS");
}
