//! 6R.307 holdout: independent re-run of the genotyping-order localization.
//! Skipped unless `HOLDOUT_6R307=1`. Does not run HOLDOUT_6R243.
//! The live call restores the pre-6R.312 name sort, which is the order
//! this localization recorded.
//!
//! ```text
//! HOLDOUT_6R307=1 cargo test -p gatk-haplotypecaller --test holdout_6r307_genotyping_order -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::assembly_region_finalize::set_forensic_6r308_java_read_coordinate_order;
use gatk_haplotypecaller::{
    begin_forensic_6r307_order_observe, call_disposition, flatten_assembly_regions,
    set_forensic_6r294_padded_span, set_forensic_6r306_exclude_failed_mate,
    take_colocated_merge_numerics, take_forensic_6r307_order, traverse_assembly_region_walker,
    AssemblyRegionCallDisposition, CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams,
    WalkerTraversalConfig,
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

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn order_once(root: &Path) -> Vec<(String, u16)> {
    set_forensic_6r294_padded_span(Some(JAVA_SPAN.0), Some(JAVA_SPAN.1));
    set_forensic_6r306_exclude_failed_mate(true);
    set_forensic_6r308_java_read_coordinate_order(false);
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
    set_forensic_6r294_padded_span(None, None);
    set_forensic_6r306_exclude_failed_mate(false);
    let numerics = take_colocated_merge_numerics();
    let site = numerics
        .iter()
        .find(|n| n.loc == LOCUS && n.alts.iter().any(|a| a == ALT))
        .expect("site");
    let geno: Vec<_> = site
        .ad_row_qname
        .iter()
        .zip(site.ad_row_flags.iter())
        .map(|(q, f)| (q.clone(), *f))
        .collect();
    let pairhmm = stages
        .iter()
        .find(|(s, _)| s == "pairhmm_input")
        .expect("pairhmm");
    let mut by = HashMap::new();
    for row in &pairhmm.1 {
        by.insert((row.qname.clone(), row.flags), row.clone());
    }
    let mut java_sorted = geno.clone();
    java_sorted.sort_by(|a, b| {
        let x = &by[a];
        let y = &by[b];
        x.start_1based
            .cmp(&y.start_1based)
            .then_with(|| x.reverse.cmp(&y.reverse))
            .then_with(|| x.qname.cmp(&y.qname))
            .then_with(|| x.flags.cmp(&y.flags))
    });
    let text = std::fs::read_to_string(root.join(JAVA_TSV)).unwrap();
    let mut java = Vec::new();
    for line in text.lines() {
        let p: Vec<_> = line.split('\t').collect();
        if p.len() >= 4 && p[0] == "6R297" && p[1] == "allele_ll" {
            java.push((p[2].to_string(), p[3].parse::<u16>().unwrap()));
        }
    }
    assert_eq!(java_sorted, java);
    geno
}

#[test]
fn holdout_6r307_genotyping_order() {
    if std::env::var("HOLDOUT_6R307").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R307=1");
        return;
    }
    let root = repo_root();
    let a = order_once(&root);
    let b = order_once(&root);
    assert_eq!(a, b);
    assert_eq!(a.len(), 122);
    assert_eq!(a[0].0, "HISEQ1:11:H8GV6ADXX:1:1109:9994:7054");
    assert_eq!(a[0].1, 83);
    assert_eq!(a[1].0, "HISEQ1:11:H8GV6ADXX:1:1116:4033:65919");
    assert_eq!(a[1].1, 99);
    println!("6R307\tclassification\tGENOTYPING_ORDER_EARLY_DIVERGENCE");
    println!("6R307\tdeterministic\tYES");
    println!("6R307\tholdout\tPASS");
}
