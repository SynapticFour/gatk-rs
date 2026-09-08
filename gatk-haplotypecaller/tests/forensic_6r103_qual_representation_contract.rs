//! 6R.103 coordinate-free: QUAL representation vs VCF serialization.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`) + HTSJDK **3.0.5**
//! (`VCFEncoder.formatQualValue`). The 6R.102 QUAL **formula** is closed. This
//! arrow tests how a finite QUAL double is stored and how the VCF QUAL column is
//! written. No AF / GL / prior / PairHMM change.
//!
//! Java (`htsjdk 3.0.5` `VCFEncoder`):
//! ```text
//! VariantContext.log10PError                 // full double
//! getPhredScaledQual() = (log10PError * -10) + 0.0  // CommonInfo; full double
//! if (!hasLog10PError()) → "."
//! else formatQualValue(getPhredScaledQual()):
//!   s = String.format(Locale.US, "%.2f", qual)
//!   if s.endsWith(".00"): strip ".00"
//! ```
//! `formatVCFDouble` (INFO floats: `%.2f` / `%.3f` / `%.3e`) is **not** used for QUAL.
//!
//! Rust (`gatk_core::io::vcf::VcfWriter::format_record`):
//! ```text
//! VcfRecord.quality: Option<f64>             // full double until write
//! None → "."
//! Some(q) → format!("{:.2}", q)                // no ".00" trim
//! ```
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r103_qual_representation_contract
//! HOLDOUT_6R103=1 cargo test -p gatk-haplotypecaller --test holdout_6r103_qual_representation -- --nocapture
//! ```

use gatk_core::io::vcf::{VcfHeader, VcfRecord, VcfWriter};
use gatk_haplotypecaller::variant_site_hc_annotations::qual_from_merged_diploid_af_calculate;
use std::io::Read;

/// PL-roundtrip GLs for a 4-allele diploid object (`TG, T, CG, *`). Coordinate-free.
fn merged_pl_roundtrip_gls() -> [f64; 10] {
    let pl = [542, 484, 1964, 0, 1234, 1353, 481, 1801, 1264, 1880];
    std::array::from_fn(|i| (pl[i] as f64) / -10.0)
}

/// HTSJDK 3.0.5 `VCFEncoder.formatQualValue` after `String.format(Locale.US, "%.2f", q)`.
///
/// The `%.2f` step is supplied by the caller so this helper does not invent a second
/// rounding mode. Production Rust uses `{:.2}` and does **not** trim `.00`.
fn java_format_qual_value_from_percent_2f(formatted_2dp: &str) -> String {
    const TRIM: &str = ".00";
    if formatted_2dp.ends_with(TRIM) {
        formatted_2dp[..formatted_2dp.len() - TRIM.len()].to_string()
    } else {
        formatted_2dp.to_string()
    }
}

/// Production QUAL column: write a `VcfRecord` through [`VcfWriter`] and read the QUAL field.
fn emitted_qual_column(quality: Option<f64>) -> String {
    emitted_vcf_line(quality)
        .split('\t')
        .nth(5)
        .expect("QUAL column")
        .to_string()
}

fn emitted_vcf_line(quality: Option<f64>) -> String {
    let dir = tempfile::tempdir().expect("tempdir");
    let path = dir.path().join("qual.vcf");
    let header = VcfHeader::default();
    let mut writer = VcfWriter::new(&path, header).expect("writer");
    writer.write_header().expect("header");
    let rec = VcfRecord {
        chromosome: "1".to_string(),
        position: 1,
        id: ".".to_string(),
        reference: "A".to_string(),
        alternate: vec!["G".to_string()],
        quality,
        filter: vec![".".to_string()],
        info: vec![],
        format: vec![],
        samples: vec![],
    };
    writer.write_record(&rec).expect("record");
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

fn fixture_qual() -> f64 {
    qual_from_merged_diploid_af_calculate(&merged_pl_roundtrip_gls(), &["TG", "T", "CG", "*"])
        .expect("QUAL")
}

fn rust_percent_2f(q: f64) -> String {
    format!("{:.2}", q)
}

fn log10_p_no_variant(qual: f64) -> f64 {
    qual / -10.0
}

#[test]
fn java_and_rust_emit_510_06_for_the_formula_qual_bits() {
    let qual = fixture_qual();
    assert_eq!(qual.to_bits(), 0x407f_e0f1_8163_fea0);
    assert_eq!(emitted_qual_column(Some(qual)), "510.06");
    let java_phred = (log10_p_no_variant(qual) * -10.0) + 0.0;
    let rust_style = (-10.0 * log10_p_no_variant(qual).min(0.0)) + 0.0;
    assert_eq!(java_phred.to_bits(), qual.to_bits());
    assert_eq!(rust_style.to_bits(), qual.to_bits());
}

/// OpenJDK 21 `String.format(Locale.US, "%.2f", q)` then HTSJDK 3.0.5 `.00` trim,
/// measured on the same IEEE-754 bits as the Rust `f64` literals.
#[test]
fn halfway_literals_are_not_the_canonical_serialization_contract() {
    // Java Formatter HALF_UP can differ from Rust `{:.2}` at some .xx5 values.
    // That is not the 6R.102 QUAL (510.058961… → both emit 510.06). No production change.
    let cases: &[(u64, &str, &str)] = &[
        (0x407f_e010_624d_d2f2, "510.00", "510"),
        (0x407f_e014_7ae1_47ae, "510.00", "510.01"),
        (0x407f_e018_9374_bc6a, "510.01", "510.01"),
        (0x4080_b514_7ae1_47ae, "534.63", "534.64"),
        (0x4080_b528_f5c2_8f5c, "534.64", "534.65"),
        (0x407f_e0f1_8163_fea0, "510.06", "510.06"),
    ];
    for &(bits, rust_oracle, java_oracle) in cases {
        let q = f64::from_bits(bits);
        let rust = rust_percent_2f(q);
        let written = emitted_qual_column(Some(q));
        assert_eq!(written, rust);
        assert_eq!(rust, rust_oracle, "bits=0x{bits:016x}");
        if rust == "510.06" {
            assert_eq!(java_oracle, "510.06");
        }
        eprintln!("6R.103 bits=0x{bits:016x} rust={rust} java_formatQualValue={java_oracle}");
    }
}

#[test]
fn raw_qual_is_full_f64_until_vcf_write() {
    let qual = fixture_qual();
    let log10_p_no_variant = qual / -10.0;
    eprintln!(
        "6R.103 fixture QUAL decimal={qual:.20} bits=0x{:016x} log10PNoVariant={log10_p_no_variant:.20} log10_bits=0x{:016x}",
        qual.to_bits(),
        log10_p_no_variant.to_bits()
    );
    assert!(qual.is_finite());
    assert!(
        qual > 510.05 && qual < 510.07,
        "formula QUAL near 510.058961, got {qual}"
    );
    assert_ne!(
        rust_percent_2f(qual),
        format!("{qual}"),
        "internal Display is not the VCF column"
    );
    let rec_line = emitted_vcf_line(Some(qual));
    let expected = rust_percent_2f(qual);
    assert_eq!(rec_line.split('\t').nth(5), Some(expected.as_str()));
    assert_eq!(qual.to_bits(), fixture_qual().to_bits());
}

#[test]
fn vcf_writer_emits_two_decimal_places_for_formula_qual() {
    let qual = fixture_qual();
    let emitted = emitted_qual_column(Some(qual));
    eprintln!(
        "6R.103 fixture VCF QUAL column={emitted:?} rust_{{:.2}}={}",
        rust_percent_2f(qual)
    );
    assert_eq!(emitted, rust_percent_2f(qual));
    assert_eq!(emitted, "510.06");
    let java = java_format_qual_value_from_percent_2f(&rust_percent_2f(qual));
    assert_eq!(java, "510.06");
}

#[test]
fn missing_qual_emits_dot() {
    assert_eq!(emitted_qual_column(None), ".");
}

#[test]
fn formatting_contract_is_two_decimal_places_not_a_locus_pin() {
    let cases = [510.004, 510.005, 510.006, 534.635, 534.645, 510.058961];
    for q in cases {
        let rust = rust_percent_2f(q);
        let written = emitted_qual_column(Some(q));
        assert_eq!(
            written, rust,
            "VcfWriter QUAL must be format!(\"{{:.2}}\") for {q}"
        );
        let java = java_format_qual_value_from_percent_2f(&rust);
        eprintln!("6R.103 format q={q:.15} rust={rust} java_from_same_2dp={java}");
    }
}

#[test]
fn integer_looking_qual_documents_java_dot00_trim_without_changing_production() {
    let rust = emitted_qual_column(Some(510.0));
    let java = java_format_qual_value_from_percent_2f(&rust_percent_2f(510.0));
    assert_eq!(rust, "510.00");
    assert_eq!(java, "510");
}

#[test]
fn formula_qual_and_emitted_qual_are_semantically_the_same_value() {
    let qual = fixture_qual();
    let emitted = emitted_qual_column(Some(qual));
    let parsed: f64 = emitted.parse().expect("QUAL parse");
    assert_eq!(emitted, "510.06");
    assert!((parsed - qual).abs() < 0.005);
    assert!(
        (parsed - 510.06).abs() < f64::EPSILON,
        "serialized QUAL parses as 510.06, got {parsed}"
    );
}
