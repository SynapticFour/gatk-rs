//! Parity dumps for likelihood engine modes (F.4–F.6).

use crate::likelihood_engine::{HcLikelihoodEngineConfig, HcLikelihoodImplementation};
use crate::pairhmm_log10::log10_pairhmm_likelihood_parity_defaults;
use crate::pcr_error_model::{error_model_adjusted_qual, PcrErrorModel};
use gatk_common::{GatkError, GatkResult};
use std::io::Write;
use std::path::Path;

/// `likelihood-engine-config` — HC default engine flags.
pub fn dump_likelihood_engine_config_tsv(out: &mut impl Write) -> GatkResult<()> {
    let cfg = HcLikelihoodEngineConfig::default();
    writeln!(out, "key\tvalue").map_err(|e| GatkError::generic(format!("write: {e}")))?;
    writeln!(out, "primary_engine\t{}", cfg.primary_engine_label())
        .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    writeln!(
        out,
        "stepwise_filtering\t{}",
        cfg.filter_step_engine_active()
    )
    .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    writeln!(out, "dragstr_pair_hmm\t{}", cfg.uses_dragstr_pair_hmm())
        .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    writeln!(
        out,
        "flow_based\t{}",
        cfg.implementation == HcLikelihoodImplementation::FlowBased
    )
    .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    Ok(())
}

/// `pcr-error-model` — adjusted ins qual at repeat lengths (CONSERVATIVE).
pub fn dump_pcr_error_model_tsv(out: &mut impl Write) -> GatkResult<()> {
    writeln!(out, "repeat_length\tadjusted_ins_qual")
        .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    let rate = PcrErrorModel::Conservative.rate_factor().unwrap();
    for repeat in 0..=10 {
        writeln!(out, "{repeat}\t{}", error_model_adjusted_qual(repeat, rate))
            .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    }
    Ok(())
}

/// `likelihood-pcr-read` — one read with PCR-adjusted likelihood vs unadjusted.
pub fn dump_likelihood_pcr_read_tsv(cases_path: &Path, out: &mut impl Write) -> GatkResult<()> {
    let text = std::fs::read_to_string(cases_path)
        .map_err(|e| GatkError::generic(format!("read {}: {e}", cases_path.display())))?;
    writeln!(out, "case_id\tmode\tlog10_likelihood")
        .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 4 {
            return Err(GatkError::argument("pcr-read cases need 4 cols"));
        }
        let case_id = cols[0];
        let read = cols[1].as_bytes();
        let quals: Vec<u8> = cols[2]
            .split(',')
            .map(|p| p.parse::<u8>().map_err(|_| GatkError::argument("bad qual")))
            .collect::<Result<_, _>>()?;
        let hap = cols[3].as_bytes();
        let ll_none = log10_pairhmm_likelihood_parity_defaults(read, &quals, hap)?;
        let mut cfg = HcLikelihoodEngineConfig::default();
        cfg.pcr_error_model = PcrErrorModel::None;
        let mapq: u8 = cols.get(4).and_then(|s| s.parse().ok()).unwrap_or(60);
        let ll_none2 = crate::likelihood_engine::log10_read_haplotype_likelihood(
            &cfg, read, &quals, mapq, hap,
        )?;
        let ll_pcr = crate::likelihood_engine::log10_read_haplotype_likelihood(
            &HcLikelihoodEngineConfig::default(),
            read,
            &quals,
            mapq,
            hap,
        )?;
        let _ = ll_none2;
        writeln!(out, "{case_id}\tnone\t{ll_none:.17}")?;
        writeln!(out, "{case_id}\tconservative\t{ll_pcr:.17}")?;
    }
    Ok(())
}
