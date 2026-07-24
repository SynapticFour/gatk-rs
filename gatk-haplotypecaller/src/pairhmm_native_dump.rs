//! GATK `Log10PairHMM` (GAP-F-02).

use crate::pairhmm_dump::load_pairhmm_cases_tsv;
use crate::pairhmm_log10::log10_pairhmm_likelihood_parity_defaults;
use gatk_common::{GatkError, GatkResult};
use std::io::Write;
use std::path::Path;

pub fn dump_pairhmm_native_likelihoods_tsv(
    cases_path: &Path,
    out: &mut impl Write,
) -> GatkResult<()> {
    let cases = load_pairhmm_cases_tsv(cases_path)?;
    writeln!(out, "case_id\tlog10_likelihood")
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    for row in &cases {
        let read = row.read_bases.as_bytes();
        let hap = row.haplotype.as_bytes();
        let ll = log10_pairhmm_likelihood_parity_defaults(read, &row.read_base_quals, hap)?;
        writeln!(out, "{}\t{ll:.17}", row.case_id)
            .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    }
    Ok(())
}
