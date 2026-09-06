//! 6R.103 holdout: emitted VCF QUAL column at the canonical T/C site.
//!
//! Skipped unless `HOLDOUT_6R103=1`. Coordinate-free contract lives in
//! `forensic_6r103_qual_representation_contract`.
//!
//! 6R.102 closed the QUAL **formula**. This holdout reads the QUAL **text**
//! written by `VcfWriter`, not the internal `VcfRecord.quality` diagnostic.
//!
//! ```text
//! HOLDOUT_6R103=1 cargo test -p gatk-haplotypecaller --test holdout_6r103_qual_representation -- --nocapture
//! ```

use gatk_core::io::vcf::{VcfHeader, VcfWriter};
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, take_colocated_merge_numerics,
    traverse_assembly_region_walker, try_emit_call_region_variants, AssemblyRegionCallDisposition,
    CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
    DEFAULT_STAND_EMIT_CONFIDENCE,
};
use serde_json::json;
use std::io::Read;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const POS_SNP: u64 = 29_456_344;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn write_and_read_body_line(rec: &gatk_core::io::vcf::VcfRecord) -> String {
    let dir = tempfile::tempdir().expect("tempdir");
    let path = dir.path().join("live.vcf");
    let mut header = VcfHeader::default();
    header.samples.push("SAMPLE".to_string());
    gatk_haplotypecaller::region_vcf_emit::populate_hc_vcf_header_schema(&mut header);
    let mut writer = VcfWriter::new(&path, header).expect("writer");
    writer.write_header().expect("header");
    writer.write_record(rec).expect("record");
    drop(writer);
    let mut text = String::new();
    std::fs::File::open(&path)
        .expect("open")
        .read_to_string(&mut text)
        .expect("read");
    text.lines()
        .find(|l| !l.starts_with('#'))
        .expect("body")
        .to_string()
}

fn format_field<'a>(format: &[&'a str], sample: &'a str, key: &str) -> &'a str {
    let idx = format
        .iter()
        .position(|k| *k == key)
        .unwrap_or_else(|| panic!("FORMAT {key}"));
    sample
        .split(':')
        .nth(idx)
        .unwrap_or_else(|| panic!("sample {key}"))
}

#[test]
fn holdout_6r103_emitted_vcf_qual_column() {
    if std::env::var("HOLDOUT_6R103").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R103=1");
        return;
    }
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    assert!(ref_fasta.is_file() && bam.is_file());

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
    let live = take_colocated_merge_numerics();
    let snap = live
        .iter()
        .find(|n| n.loc == POS_SNP)
        .expect("colocated merge numerics");

    let raw_qual = vcf.quality.expect("QUAL");
    let log10_p_no_variant = raw_qual / -10.0;
    let line = write_and_read_body_line(vcf);
    let cols: Vec<&str> = line.split('\t').collect();
    assert!(cols.len() >= 10, "VCF columns: {line}");
    let chrom = cols[0];
    let pos = cols[1];
    let rec_id = cols[2];
    let ref_al = cols[3];
    let alt = cols[4];
    let qual_col = cols[5];
    let filter = cols[6];
    let info = cols[7];
    let format = cols[8];
    let sample = cols[9];
    let format_keys: Vec<&str> = format.split(':').collect();
    let gt = format_field(&format_keys, sample, "GT");
    let ad = format_field(&format_keys, sample, "AD");
    let pl = format_field(&format_keys, sample, "PL");

    let doc = json!({
        "locus": "20:29456344 T/C",
        "internal": {
            "qual_decimal": raw_qual,
            "qual_bits": format!("0x{:016x}", raw_qual.to_bits()),
            "log10PNoVariant": log10_p_no_variant,
            "log10PNoVariant_bits": format!("0x{:016x}", log10_p_no_variant.to_bits()),
        },
        "emitted_line": line,
        "emitted": {
            "chrom": chrom,
            "pos": pos,
            "id": rec_id,
            "ref": ref_al,
            "alt": alt,
            "qual": qual_col,
            "filter": filter,
            "info": info,
            "format": format,
            "gt": gt,
            "ad": ad,
            "pl": pl,
        },
        "java_oracle": {"gt": "0/1", "ad": "36,19", "pl": "542,0,1353", "qual": "510.06"},
        "merged_has_span_del": snap.alts.iter().any(|a| a == "*"),
        "classification": "QUAL_CLOSED",
        "production_change": "NONE",
    });
    println!("{}", serde_json::to_string_pretty(&doc).unwrap());

    assert_eq!(chrom, vcf.chromosome);
    assert_eq!(pos, "29456344");
    assert_eq!(ref_al, "T");
    assert_eq!(alt, "C");
    assert_eq!(gt, "0/1");
    assert_eq!(ad, "36,19");
    assert_eq!(pl, "542,0,1353");
    assert_eq!(qual_col, "510.06");
    assert_eq!(format!("{:.2}", raw_qual), "510.06");
}
