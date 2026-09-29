//! 6R.308 holdout: independent re-run of the Java coordinate-comparator counterfactual.
//! Skipped unless `HOLDOUT_6R308=1`. Does not run HOLDOUT_6R243.
//!
//! ```text
//! HOLDOUT_6R308=1 cargo test -p gatk-haplotypecaller --test holdout_6r308_java_read_coordinate_comparator -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::assembly_region_finalize::set_forensic_6r308_java_read_coordinate_order;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, set_forensic_6r294_padded_span,
    set_forensic_6r306_exclude_failed_mate, take_colocated_merge_numerics,
    traverse_assembly_region_walker, try_emit_call_region_variants, AssemblyRegionCallDisposition,
    CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
    DEFAULT_STAND_EMIT_CONFIDENCE,
};
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

fn once(root: &Path) -> (Vec<(String, u16)>, Vec<f64>, String) {
    set_forensic_6r294_padded_span(Some(JAVA_SPAN.0), Some(JAVA_SPAN.1));
    set_forensic_6r306_exclude_failed_mate(true);
    set_forensic_6r308_java_read_coordinate_order(true);
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
    set_forensic_6r294_padded_span(None, None);
    set_forensic_6r306_exclude_failed_mate(false);
    set_forensic_6r308_java_read_coordinate_order(false);
    let numerics = take_colocated_merge_numerics();
    let site = numerics
        .iter()
        .find(|n| n.loc == LOCUS && n.alts.iter().any(|a| a == ALT))
        .expect("site");
    let geno = site
        .ad_row_qname
        .iter()
        .zip(site.ad_row_flags.iter())
        .map(|(q, f)| (q.clone(), *f))
        .collect();
    let recs =
        try_emit_call_region_variants(region, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    let emitted = recs
        .iter()
        .find(|r| r.position == LOCUS && r.alternate.iter().any(|a| a == ALT))
        .expect("emit");
    let subset = emitted.samples[0]
        .pl
        .clone()
        .unwrap_or_default()
        .iter()
        .map(|p| p.to_string())
        .collect::<Vec<_>>()
        .join(",");
    (geno, site.merged_gls.clone(), subset)
}

#[test]
fn holdout_6r308_java_read_coordinate_comparator() {
    if std::env::var("HOLDOUT_6R308").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R308=1");
        return;
    }
    let root = repo_root();
    let text = std::fs::read_to_string(root.join(JAVA_TSV)).unwrap();
    let java_gl = field_cells(&text, "java_gl_vector");
    let mut java = Vec::new();
    for line in text.lines() {
        let p: Vec<_> = line.split('\t').collect();
        if p.len() >= 4 && p[0] == "6R297" && p[1] == "allele_ll" {
            java.push((p[2].to_string(), p[3].parse::<u16>().unwrap()));
        }
    }
    let (a, gl, subset) = once(&root);
    let (b, _, _) = once(&root);
    assert_eq!(a, b);
    assert_eq!(a, java);
    assert!(gl
        .iter()
        .zip(java_gl.iter())
        .all(|(r, j)| r.to_bits() == j.to_bits()));
    assert_eq!(subset, "570,0,3517");
    println!("6R308\tclassification\tJAVA_READ_COORDINATE_COMPARATOR_REPRODUCES_JAVA");
    println!("6R308\tdeterministic\tYES");
    println!("6R308\tholdout\tPASS");
}
