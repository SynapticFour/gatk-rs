//! GATK-compatible `VariantFiltration` hard-filtering (INFO/QUAL JEXL subset).
//! Observable Java contract (GATK 4.4 `VariantFiltration`):
//! Each `--filter-expression` / `--filter-name` pair is evaluated independently.
//! If the referenced annotation is **missing**, that expression does **not** fail
//! (site stays unfiltered for that name) — same as GATK JEXL missing-attribute behavior
//! for single-attribute expressions.
//! Failing expressions append their filter name to FILTER (`;`-joined, application order).
//! Sites that fail nothing are emitted with `PASS`.
//! This is hard-filtering, **not** VQSR. Official GATK guidance recommends VQSR when the
//! cohort is large enough to train; hard filters are the pragmatic fallback for smaller
//! callsets (see GATK “Filter variants either with VQSR or by hard-filtering”).

use crate::io::vcf::{FilterField, InfoValue, VcfReader, VcfRecord, VcfWriter};
use gatk_common::{GatkError, GatkResult};
use std::path::{Path, PathBuf};

/// Official GATK Best Practices hard-filter expressions for **SNPs**
/// ([How to: Filter variants either with VQSR or by hard-filtering](https://gatk.broadinstitute.org/hc/en-us/articles/360035531112)).
pub const GATK_HARD_FILTER_SNP: &[(&str, &str)] = &[
    ("QD < 2.0", "QD2"),
    ("QUAL < 30.0", "QUAL30"),
    ("SOR > 3.0", "SOR3"),
    ("FS > 60.0", "FS60"),
    ("MQ < 40.0", "MQ40"),
    ("MQRankSum < -12.5", "MQRankSum-12.5"),
    ("ReadPosRankSum < -8.0", "ReadPosRankSum-8"),
];

/// Official GATK Best Practices hard-filter expressions for **indels**
/// (same article; `InbreedingCoeff` omitted — only defined for ≥10 samples).
pub const GATK_HARD_FILTER_INDEL: &[(&str, &str)] = &[
    ("QD < 2.0", "QD2"),
    ("QUAL < 30.0", "QUAL30"),
    ("FS > 200.0", "FS200"),
    ("ReadPosRankSum < -20.0", "ReadPosRankSum-20"),
    ("SOR > 10.0", "SOR10"),
];

/// One named filter expression (GATK `--filter-expression` / `--filter-name`).
#[derive(Debug, Clone, PartialEq)]
pub struct FilterSpec {
    pub expression: String,
    pub name: String,
}

/// CLI / library args for VariantFiltration.
#[derive(Debug, Clone)]
pub struct VariantFiltrationArgs {
    pub variant: PathBuf,
    pub output: PathBuf,
    pub filters: Vec<FilterSpec>,
    /// Optional reference (accepted for GATK CLI familiarity; unused for hard filters).
    pub reference: Option<PathBuf>,
}

/// Comparison operator in a simple JEXL-like filter expression.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CmpOp {
    Lt,
    Le,
    Gt,
    Ge,
    Eq,
    Ne,
}

/// Parsed `ANNOTATION OP NUMBER` expression.
#[derive(Debug, Clone, PartialEq)]
pub struct ParsedFilterExpr {
    pub annotation: String,
    pub op: CmpOp,
    pub threshold: f64,
}

/// Parse a single-attribute comparison, e.g. `QD < 2.0` or `MQRankSum < -12.5`.
pub fn parse_filter_expression(expr: &str) -> GatkResult<ParsedFilterExpr> {
    let s = expr.trim();
    if s.is_empty() {
        return Err(GatkError::argument("empty filter expression"));
    }
    // Longest operators first.
    let ops: &[(&str, CmpOp)] = &[
        ("<=", CmpOp::Le),
        (">=", CmpOp::Ge),
        ("==", CmpOp::Eq),
        ("!=", CmpOp::Ne),
        ("<", CmpOp::Lt),
        (">", CmpOp::Gt),
        ("=", CmpOp::Eq),
    ];
    for &(tok, op) in ops {
        if let Some(idx) = s.find(tok) {
            let left = s[..idx].trim();
            let right = s[idx + tok.len()..].trim();
            if left.is_empty() || right.is_empty() {
                continue;
            }
            // Reject compound expressions in this slice (&& / ||).
            if left.contains("&&")
                || left.contains("||")
                || right.contains("&&")
                || right.contains("||")
            {
                return Err(GatkError::argument(
                    "compound filter expressions (&&/||) are not supported; \
                     pass each criterion as its own --filter-expression (GATK recommendation)",
                ));
            }
            let threshold: f64 = right.parse().map_err(|_| {
                GatkError::argument(format!(
                    "filter expression threshold is not a number: '{right}' in '{expr}'"
                ))
            })?;
            if !is_ident(left) {
                return Err(GatkError::argument(format!(
                    "filter annotation must be an identifier (got '{left}')"
                )));
            }
            return Ok(ParsedFilterExpr {
                annotation: left.to_string(),
                op,
                threshold,
            });
        }
    }
    Err(GatkError::argument(format!(
        "unsupported filter expression (need 'ANNOTATION OP NUMBER'): '{expr}'"
    )))
}

fn is_ident(s: &str) -> bool {
    let mut chars = s.chars();
    match chars.next() {
        Some(c) if c.is_ascii_alphabetic() || c == '_' => {}
        _ => return false,
    }
    chars.all(|c| c.is_ascii_alphanumeric() || c == '_')
}

/// Build filter specs from GATK SNP hard-filter recommendations.
pub fn gatk_snp_hard_filters() -> Vec<FilterSpec> {
    GATK_HARD_FILTER_SNP
        .iter()
        .map(|&(e, n)| FilterSpec {
            expression: e.to_string(),
            name: n.to_string(),
        })
        .collect()
}

/// Build filter specs from GATK indel hard-filter recommendations.
pub fn gatk_indel_hard_filters() -> Vec<FilterSpec> {
    GATK_HARD_FILTER_INDEL
        .iter()
        .map(|&(e, n)| FilterSpec {
            expression: e.to_string(),
            name: n.to_string(),
        })
        .collect()
}

/// Pair expression/name vectors (must be equal length).
pub fn zip_filter_pairs(expressions: &[String], names: &[String]) -> GatkResult<Vec<FilterSpec>> {
    if expressions.len() != names.len() {
        return Err(GatkError::argument(format!(
            "VariantFiltration requires matching --filter-expression / --filter-name counts \
             (got {} expressions, {} names)",
            expressions.len(),
            names.len()
        )));
    }
    if expressions.is_empty() {
        return Err(GatkError::argument(
            "VariantFiltration requires at least one --filter-expression / --filter-name pair \
             (or --preset snp|indel)",
        ));
    }
    Ok(expressions
        .iter()
        .zip(names.iter())
        .map(|(e, n)| FilterSpec {
            expression: e.clone(),
            name: n.clone(),
        })
        .collect())
}

/// Evaluate whether a site **fails** a single filter (true → apply filter name).
/// Missing annotation → `Ok(false)` (do not fail), matching GATK JEXL missing behavior.
pub fn expression_fails(rec: &VcfRecord, parsed: &ParsedFilterExpr) -> bool {
    let Some(value) = annotation_value(rec, &parsed.annotation) else {
        return false;
    };
    compare(value, parsed.op, parsed.threshold)
}

fn compare(value: f64, op: CmpOp, threshold: f64) -> bool {
    match op {
        CmpOp::Lt => value < threshold,
        CmpOp::Le => value <= threshold,
        CmpOp::Gt => value > threshold,
        CmpOp::Ge => value >= threshold,
        CmpOp::Eq => (value - threshold).abs() < 1e-12,
        CmpOp::Ne => (value - threshold).abs() >= 1e-12,
    }
}

/// htsjdk / GATK sentinel for missing QUAL (`VCFConstants.MISSING_VALUE_v4` → −10.0).
/// So `QUAL < 30.0` fails when QUAL is `.` in the VCF — matching Java VariantFiltration.
pub const MISSING_QUAL_SENTINEL: f64 = -10.0;

fn annotation_value(rec: &VcfRecord, name: &str) -> Option<f64> {
    if name.eq_ignore_ascii_case("QUAL") {
        return Some(rec.quality.unwrap_or(MISSING_QUAL_SENTINEL));
    }
    for info in &rec.info {
        match info {
            InfoValue::Float(id, vals) if id == name => return vals.first().copied(),
            InfoValue::Integer(id, vals) if id == name => {
                return vals.first().map(|&v| v as f64);
            }
            _ => {}
        }
    }
    None
}

/// Apply filters to one record; returns FILTER tokens (`PASS` or failing names).
pub fn apply_filters_to_record(rec: &VcfRecord, filters: &[FilterSpec]) -> GatkResult<Vec<String>> {
    let mut failed = Vec::new();
    for spec in filters {
        let parsed = parse_filter_expression(&spec.expression)?;
        if expression_fails(rec, &parsed) {
            // CLONE: needed because owned element into collection.
            failed.push(spec.name.clone());
        }
    }
    if failed.is_empty() {
        Ok(vec!["PASS".to_string()])
    } else {
        Ok(failed)
    }
}

/// Soft-filter all records in memory.
pub fn filter_records(records: &[VcfRecord], filters: &[FilterSpec]) -> GatkResult<Vec<VcfRecord>> {
    // Pre-parse expressions once.
    let parsed: Vec<(String, ParsedFilterExpr)> = filters
        .iter()
        .map(|f| parse_filter_expression(&f.expression).map(|p| (f.name.clone(), p)))
        .collect::<GatkResult<_>>()?;

    let mut out = Vec::with_capacity(records.len());
    for rec in records {
        let mut failed = Vec::new();
        for (name, expr) in &parsed {
            if expression_fails(rec, expr) {
                // CLONE: needed because owned element into collection.
                failed.push(name.clone());
            }
        }
        let mut cloned = rec.clone();
        cloned.filter = if failed.is_empty() {
            vec!["PASS".to_string()]
        } else {
            failed
        };
        out.push(cloned);
    }
    Ok(out)
}

/// Run VariantFiltration from filesystem paths.
pub fn run_variant_filtration(args: &VariantFiltrationArgs) -> GatkResult<()> {
    if args.filters.is_empty() {
        return Err(GatkError::argument(
            "VariantFiltration requires at least one filter",
        ));
    }
    let mut reader = VcfReader::from_file(&args.variant)?;
    let mut header = reader.header().clone();
    for f in &args.filters {
        // Validate expression early.
        let _ = parse_filter_expression(&f.expression)?;
        if !header.filter_fields.iter().any(|h| h.id == f.name) {
            header.filter_fields.push(FilterField {
                id: f.name.clone(),
                description: format!("Hard filter: {}", f.expression),
            });
        }
    }
    if !header.filter_fields.iter().any(|h| h.id == "PASS") {
        header.filter_fields.insert(
            0,
            FilterField {
                id: "PASS".to_string(),
                description: "All filters passed".to_string(),
            },
        );
    }
    header.source = Some("gatk-rs VariantFiltration".to_string());
    let parsed = args
        .filters
        .iter()
        .map(|f| parse_filter_expression(&f.expression).map(|p| (f.name.clone(), p)))
        .collect::<GatkResult<Vec<_>>>()?;
    let mut writer = VcfWriter::new(&args.output, header)?;
    writer.write_header()?;
    while let Some(rec) = reader.read_next_record()? {
        let mut failed = Vec::new();
        for (name, expr) in &parsed {
            if expression_fails(&rec, expr) {
                failed.push(name.clone());
            }
        }
        let mut cloned = rec;
        cloned.filter = if failed.is_empty() {
            vec!["PASS".to_string()]
        } else {
            failed
        };
        writer.write_record(&cloned)?;
    }
    let _ = args.reference.as_ref(); // accepted for CLI familiarity
    Ok(())
}

/// Convenience: write filtered VCF using a preset.
pub fn run_variant_filtration_preset(
    variant: &Path,
    output: &Path,
    preset: &str,
) -> GatkResult<()> {
    let filters = match preset.to_ascii_lowercase().as_str() {
        "snp" | "snps" => gatk_snp_hard_filters(),
        "indel" | "indels" => gatk_indel_hard_filters(),
        other => {
            return Err(GatkError::argument(format!(
                "unknown --preset '{other}' (use snp or indel)"
            )));
        }
    };
    run_variant_filtration(&VariantFiltrationArgs {
        variant: variant.to_path_buf(),
        output: output.to_path_buf(),
        filters,
        reference: None,
    })
}

// --------------------------------------------------------------------------
// Tests — boundary decisions matching GATK JEXL / VariantFiltration
// --------------------------------------------------------------------------

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;
    use crate::io::vcf::{Genotype, SampleData};

    fn rec_with(info: Vec<InfoValue>, qual: Option<f64>) -> VcfRecord {
        VcfRecord {
            chromosome: "chr1".to_string(),
            position: 100,
            id: ".".to_string(),
            reference: "A".to_string(),
            alternate: vec!["G".to_string()],
            quality: qual,
            filter: vec![".".to_string()],
            info,
            format: vec!["GT".to_string()],
            samples: vec![SampleData {
                gt: Some(Genotype {
                    alleles: vec![0, 1],
                    phased: false,
                }),
                gq: None,
                dp: None,
                ad: None,
                pl: None,
                other: Vec::new(),
            }],
        }
    }

    fn f(id: &str, v: f64) -> InfoValue {
        InfoValue::Float(id.to_string(), vec![v])
    }

    fn fails(expr: &str, rec: &VcfRecord) -> bool {
        let p = parse_filter_expression(expr).unwrap();
        expression_fails(rec, &p)
    }

    #[test]
    fn parse_qd_lt() {
        let p = parse_filter_expression("QD < 2.0").unwrap();
        assert_eq!(p.annotation, "QD");
        assert_eq!(p.op, CmpOp::Lt);
        assert!((p.threshold - 2.0).abs() < 1e-12);
    }

    #[test]
    fn parse_negative_threshold() {
        let p = parse_filter_expression("MQRankSum < -12.5").unwrap();
        assert!((p.threshold - (-12.5)).abs() < 1e-12);
    }

    #[test]
    fn qd_boundary_strict_lt() {
        // Java: QD < 2.0 → fail only when strictly below 2.0
        let at = rec_with(vec![f("QD", 2.0)], Some(40.0));
        let below = rec_with(vec![f("QD", 1.999)], Some(40.0));
        let above = rec_with(vec![f("QD", 2.001)], Some(40.0));
        assert!(!fails("QD < 2.0", &at));
        assert!(fails("QD < 2.0", &below));
        assert!(!fails("QD < 2.0", &above));
    }

    #[test]
    fn fs_boundary_strict_gt_snp() {
        let at = rec_with(vec![f("FS", 60.0)], Some(40.0));
        let above = rec_with(vec![f("FS", 60.001)], Some(40.0));
        let below = rec_with(vec![f("FS", 59.999)], Some(40.0));
        assert!(!fails("FS > 60.0", &at));
        assert!(fails("FS > 60.0", &above));
        assert!(!fails("FS > 60.0", &below));
    }

    #[test]
    fn fs_boundary_indel_200() {
        let at = rec_with(vec![f("FS", 200.0)], Some(40.0));
        let above = rec_with(vec![f("FS", 200.1)], Some(40.0));
        assert!(!fails("FS > 200.0", &at));
        assert!(fails("FS > 200.0", &above));
    }

    #[test]
    fn mq_boundary_strict_lt() {
        let at = rec_with(vec![f("MQ", 40.0)], Some(40.0));
        let below = rec_with(vec![f("MQ", 39.999)], Some(40.0));
        assert!(!fails("MQ < 40.0", &at));
        assert!(fails("MQ < 40.0", &below));
    }

    #[test]
    fn mqranksum_boundary() {
        let at = rec_with(vec![f("MQRankSum", -12.5)], Some(40.0));
        let below = rec_with(vec![f("MQRankSum", -12.5001)], Some(40.0));
        assert!(!fails("MQRankSum < -12.5", &at));
        assert!(fails("MQRankSum < -12.5", &below));
    }

    #[test]
    fn readposranksum_snp_and_indel_thresholds() {
        let mid = rec_with(vec![f("ReadPosRankSum", -10.0)], Some(40.0));
        assert!(fails("ReadPosRankSum < -8.0", &mid)); // SNP threshold
        assert!(!fails("ReadPosRankSum < -20.0", &mid)); // indel threshold
        let deep = rec_with(vec![f("ReadPosRankSum", -20.1)], Some(40.0));
        assert!(fails("ReadPosRankSum < -20.0", &deep));
    }

    #[test]
    fn sor_snp_and_indel_thresholds() {
        let v = rec_with(vec![f("SOR", 5.0)], Some(40.0));
        assert!(fails("SOR > 3.0", &v));
        assert!(!fails("SOR > 10.0", &v));
        let high = rec_with(vec![f("SOR", 10.1)], Some(40.0));
        assert!(fails("SOR > 10.0", &high));
    }

    #[test]
    fn qual_boundary() {
        let at = rec_with(vec![], Some(30.0));
        let below = rec_with(vec![], Some(29.999));
        assert!(!fails("QUAL < 30.0", &at));
        assert!(fails("QUAL < 30.0", &below));
    }

    #[test]
    fn missing_annotation_does_not_fail() {
        let rec = rec_with(vec![f("FS", 100.0)], Some(40.0));
        // QD missing → QD2 must not apply
        assert!(!fails("QD < 2.0", &rec));
        assert!(fails("FS > 60.0", &rec));
    }

    #[test]
    fn missing_qual_uses_htsjdk_sentinel_and_fails_qual_lt_30() {
        // Java/htsjdk: missing QUAL ≡ -10.0 → fails `QUAL < 30.0`.
        let rec = rec_with(vec![f("QD", 10.0)], None);
        assert!(fails("QUAL < 30.0", &rec));
        assert!(!fails("QUAL < -10.0", &rec)); // -10 < -10 is false
    }

    #[test]
    fn multiple_filters_accumulate_names() {
        let rec = rec_with(
            vec![f("QD", 1.0), f("FS", 100.0), f("MQ", 50.0)],
            Some(40.0),
        );
        let filters = vec![
            FilterSpec {
                expression: "QD < 2.0".into(),
                name: "QD2".into(),
            },
            FilterSpec {
                expression: "FS > 60.0".into(),
                name: "FS60".into(),
            },
            FilterSpec {
                expression: "MQ < 40.0".into(),
                name: "MQ40".into(),
            },
        ];
        let names = apply_filters_to_record(&rec, &filters).unwrap();
        assert_eq!(names, vec!["QD2", "FS60"]);
    }

    #[test]
    fn all_pass_emits_pass() {
        let rec = rec_with(
            vec![
                f("QD", 10.0),
                f("FS", 1.0),
                f("SOR", 0.5),
                f("MQ", 60.0),
                f("MQRankSum", 0.0),
                f("ReadPosRankSum", 0.0),
            ],
            Some(99.0),
        );
        let names = apply_filters_to_record(&rec, &gatk_snp_hard_filters()).unwrap();
        assert_eq!(names, vec!["PASS"]);
    }

    #[test]
    fn snp_preset_matches_official_table() {
        assert_eq!(GATK_HARD_FILTER_SNP.len(), 7);
        assert!(GATK_HARD_FILTER_SNP
            .iter()
            .any(|&(e, _)| e.contains("FS > 60")));
        assert!(GATK_HARD_FILTER_INDEL
            .iter()
            .any(|&(e, _)| e.contains("FS > 200")));
    }

    #[test]
    fn zip_pairs_length_mismatch_errors() {
        let err = zip_filter_pairs(&["QD < 2.0".into()], &[]).unwrap_err();
        assert!(format!("{err}").contains("matching"));
    }

    #[test]
    fn integer_info_compared_as_float() {
        let rec = rec_with(
            vec![InfoValue::Integer("DP".to_string(), vec![5])],
            Some(40.0),
        );
        assert!(fails("DP < 10", &rec));
        assert!(!fails("DP < 5", &rec));
    }

    #[test]
    fn filter_records_sets_pass_and_fail() {
        let good = rec_with(vec![f("QD", 10.0)], Some(40.0));
        let bad = rec_with(vec![f("QD", 1.0)], Some(40.0));
        let filters = gatk_snp_hard_filters()
            .into_iter()
            .filter(|f| f.name == "QD2")
            .collect::<Vec<_>>();
        let out = filter_records(&[good, bad], &filters).unwrap();
        assert_eq!(out[0].filter, vec!["PASS"]);
        assert_eq!(out[1].filter, vec!["QD2"]);
    }
}
