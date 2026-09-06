//! 6R.70 forensic: Java vs Rust PCR repeat length on the canonical first read.
//!
//! Skipped unless `HOLDOUT_6R70=1`. Coordinate-free proof is
//! `forensic_6r70_pcr_repeat_contract`.
//!
//! ```text
//! HOLDOUT_6R70=1 cargo test -p gatk-haplotypecaller --test holdout_6r70_pcr_repeat -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::pcr_error_model::{
    error_model_adjusted_qual, find_tandem_repeat_units, PcrErrorModel,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, take_colocated_merge_numerics,
    traverse_assembly_region_walker, try_emit_call_region_variants, AssemblyRegionCallDisposition,
    CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
    DEFAULT_STAND_EMIT_CONFIDENCE,
};
use rust_htslib::bam::record::Aux;
use serde_json::json;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const POS_SNP: u64 = 29_456_344;
const JAVA_MAX_STR_UNIT: usize = 8;
const JAVA_MAX_REPEAT: usize = 20;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn aux_phred_z(rec: &rust_htslib::bam::Record, tag: &[u8]) -> Option<Vec<u8>> {
    match rec.aux(tag) {
        Ok(Aux::String(s)) => Some(s.bytes().map(|b| b.saturating_sub(33)).collect()),
        _ => None,
    }
}

fn java_fast_round(d: f64) -> i32 {
    if d > 0.0 {
        (d + 0.5) as i32
    } else {
        (d - 0.5) as i32
    }
}

fn java_44_cache(r: usize) -> u8 {
    let q = 40.0 - (r as f64 / (3.0 * std::f64::consts::PI)).exp() + 1.0;
    java_fast_round(q).max(10) as u8
}

fn equal_range(a: &[u8], a_off: usize, b: &[u8], b_off: usize, len: usize) -> bool {
    a.get(a_off..a_off + len) == b.get(b_off..b_off + len)
}

fn find_number_of_repetitions(
    unit: &[u8],
    unit_off: usize,
    unit_len: usize,
    test: &[u8],
    test_off: usize,
    test_len: usize,
    leading: bool,
) -> usize {
    if unit_len == 0 || test_len == 0 {
        return 0;
    }
    let length_difference = test_len as isize - unit_len as isize;
    if leading {
        let mut n = 0usize;
        let mut start = 0isize;
        while start <= length_difference {
            if equal_range(test, start as usize + test_off, unit, unit_off, unit_len) {
                n += 1;
                start += unit_len as isize;
            } else {
                return n;
            }
        }
        n
    } else {
        let mut n = 0usize;
        let mut start = length_difference;
        while start >= 0 {
            if equal_range(test, start as usize + test_off, unit, unit_off, unit_len) {
                n += 1;
                start -= unit_len as isize;
            } else {
                return n;
            }
        }
        n
    }
}

fn java_find_tandem_repeat_units(read: &[u8], offset: usize) -> (Vec<u8>, usize) {
    let mut max_bw = 0usize;
    let mut best_bw: Vec<u8> = vec![read[offset]];
    for str_len in 1..=JAVA_MAX_STR_UNIT {
        if (offset + 1).checked_sub(str_len).is_none() {
            break;
        }
        max_bw = find_number_of_repetitions(
            read,
            offset + 1 - str_len,
            str_len,
            read,
            0,
            offset + 1,
            false,
        );
        if max_bw > 1 {
            best_bw = read[offset + 1 - str_len..=offset].to_vec();
            break;
        }
    }
    let mut best_unit = best_bw.clone();
    let mut max_rl = max_bw;
    if offset < read.len() - 1 {
        let mut best_fw: Vec<u8> = vec![read[offset + 1]];
        let mut max_fw = 0usize;
        for str_len in 1..=JAVA_MAX_STR_UNIT {
            if offset + str_len + 1 > read.len() {
                break;
            }
            max_fw = find_number_of_repetitions(
                read,
                offset + 1,
                str_len,
                read,
                offset + 1,
                read.len() - offset - 1,
                true,
            );
            if max_fw > 1 {
                best_fw = read[offset + 1..offset + 1 + str_len].to_vec();
                break;
            }
        }
        if best_fw == best_bw {
            max_rl = max_bw + max_fw;
            best_unit = best_fw;
        } else {
            let test = &read[..=offset];
            max_bw =
                find_number_of_repetitions(&best_fw, 0, best_fw.len(), test, 0, test.len(), false);
            max_rl = max_fw + max_bw;
            best_unit = best_fw;
        }
    }
    if max_rl > JAVA_MAX_REPEAT {
        max_rl = JAVA_MAX_REPEAT;
    }
    (best_unit, max_rl)
}

#[test]
fn holdout_6r70_canonical_read0_base0_repeat_and_cache() {
    if std::env::var("HOLDOUT_6R70").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R70=1");
        return;
    }
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
    let covering: Vec<_> = regions
        .iter()
        .filter(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= POS_SNP
                && r.end.get() >= POS_SNP
        })
        .collect();
    assert_eq!(covering.len(), 1);
    let region = covering[0];
    let outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let emitted =
        try_emit_call_region_variants(region, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    let vcf = emitted
        .iter()
        .find(|r| {
            r.position == POS_SNP && r.reference == "T" && r.alternate.iter().any(|a| a == "C")
        })
        .expect("T/C");
    assert_eq!(vcf.samples[0].pl.clone().unwrap(), vec![542, 0, 1353]);
    let live = take_colocated_merge_numerics();
    // 6R.84: Java createAlleleMapper leaves spanning-del haplotypes out of REF.
    assert_eq!(
        live.iter().find(|n| n.loc == POS_SNP).unwrap().pool_sizes,
        vec![35, 6, 21, 6]
    );

    let rec = outcome
        .genotyping_reads
        .iter()
        .find(|r| aux_phred_z(r, b"BI").is_some() && aux_phred_z(r, b"BD").is_some())
        .expect("BI/BD genotyping read");
    let bases = rec.seq().as_bytes();
    let bi = aux_phred_z(rec, b"BI").expect("BI");
    let bd = aux_phred_z(rec, b"BD").expect("BD");
    let (j_unit, j_len) = java_find_tandem_repeat_units(&bases, 0);
    let (r_unit, r_len) = find_tandem_repeat_units(&bases, 0);
    let rust_cache = error_model_adjusted_qual(
        r_len.min(20),
        PcrErrorModel::Conservative.rate_factor().unwrap(),
    );
    let java_cap = java_44_cache(j_len.min(20));
    let prefix: String = bases.iter().take(8).map(|&b| b as char).collect();

    let doc = json!({
        "prefix8": prefix,
        "bi0": bi[0],
        "bd0": bd[0],
        "java_unit": String::from_utf8_lossy(&j_unit),
        "java_repeat_len": j_len,
        "java_cache": java_cap,
        "java_gop": bi[0].min(java_cap),
        "rust_unit": String::from_utf8_lossy(&r_unit),
        "rust_repeat_len": r_len,
        "rust_cache": rust_cache,
        "rust_gop_if_bi_sub": bi[0].min(rust_cache),
        "rust_production_gop": bi[0].min(rust_cache),
    });
    println!("{}", serde_json::to_string_pretty(&doc).unwrap());

    assert_eq!(r_unit, j_unit, "production tandem unit matches Java");
    assert_eq!(r_len, j_len, "production tandem count matches Java");
    assert_eq!(rust_cache, java_cap, "CONSERVATIVE cache matches Java");
    assert_eq!(bi[0].min(java_cap), bi[0].min(rust_cache));
    assert_ne!(bi[0], 45u8, "pre-PCR GOP source is BI, not Q45");
    assert_eq!(java_44_cache(1), 40);
    assert_eq!(java_44_cache(2), 40);
}
