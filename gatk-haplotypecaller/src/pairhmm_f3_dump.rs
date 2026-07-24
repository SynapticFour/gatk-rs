//! Parity dumps: BQ caps and haplotype support filtering.

use crate::pairhmm_log10::log10_pairhmm_likelihood_parity_defaults;
use crate::pairhmm_qual::cap_read_base_qualities;
use gatk_common::{GatkError, GatkResult};
use std::io::Write;
use std::path::Path;

fn parse_quals_csv(s: &str) -> GatkResult<Vec<u8>> {
    s.split(',')
        .map(|p| {
            p.trim()
                .parse::<u8>()
                .map_err(|_| GatkError::argument(format!("invalid qual: {p}")))
        })
        .collect()
}

fn format_ll(v: f64) -> String {
    if v.is_infinite() && v.is_sign_negative() {
        "-inf".to_string()
    } else {
        format!("{v:.17}")
    }
}

/// `pairhmm-bq-cap <cases.tsv>` — capped base qualities per row.
pub fn dump_pairhmm_bq_cap_tsv(cases_path: &Path, out: &mut impl Write) -> GatkResult<()> {
    let text = std::fs::read_to_string(cases_path)
        .map_err(|e| GatkError::generic(format!("read {}: {e}", cases_path.display())))?;
    writeln!(out, "case_id\tcapped_quals")?;
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 4 {
            return Err(GatkError::argument("bq-cap cases need 4 cols"));
        }
        let case_id = cols[0];
        let mut quals = parse_quals_csv(cols[1])?;
        let threshold: u8 = cols[2]
            .parse()
            .map_err(|_| GatkError::argument("invalid threshold"))?;
        let mapq: u8 = cols[3]
            .parse()
            .map_err(|_| GatkError::argument("invalid mapq"))?;
        cap_read_base_qualities(&mut quals, mapq, threshold, true);
        let capped = quals
            .iter()
            .map(|q| q.to_string())
            .collect::<Vec<_>>()
            .join(",");
        writeln!(out, "{case_id}\t{capped}")?;
    }
    Ok(())
}

/// `pairhmm-haplotype-filter <cases.tsv>` — keep haplotypes with max read LL above threshold.
pub fn dump_pairhmm_haplotype_filter_tsv(
    cases_path: &Path,
    out: &mut impl Write,
) -> GatkResult<()> {
    let text = std::fs::read_to_string(cases_path)
        .map_err(|e| GatkError::generic(format!("read {}: {e}", cases_path.display())))?;
    writeln!(out, "case_id\thaplotype_index\tkept\tmax_log10_likelihood")?;
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 5 {
            return Err(GatkError::argument("hap-filter cases need 5 cols"));
        }
        let case_id = cols[0];
        let read = cols[1].as_bytes();
        let quals = parse_quals_csv(cols[2])?;
        let threshold: f64 = cols[3]
            .parse()
            .map_err(|_| GatkError::argument("invalid ll threshold"))?;
        let haps: Vec<&[u8]> = cols[4..].iter().map(|s| s.as_bytes()).collect();
        for (hi, hap) in haps.iter().enumerate() {
            let ll = log10_pairhmm_likelihood_parity_defaults(read, &quals, hap)?;
            let kept = if ll > threshold { "true" } else { "false" };
            writeln!(out, "{case_id}\t{hi}\t{kept}\t{}", format_ll(ll))?;
        }
    }
    Ok(())
}
