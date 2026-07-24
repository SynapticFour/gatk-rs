//! PRE-D01 — DRAGSTR calibration parity scaffold.

use crate::likelihood_engine::HcLikelihoodEngineConfig;
use gatk_common::GatkResult;
use std::io::Write;

pub fn dump_dragstr_calibration_tsv(params_loaded: bool, out: &mut impl Write) -> GatkResult<()> {
    let mut cfg = HcLikelihoodEngineConfig::default();
    cfg.dragstr_params_loaded = params_loaded;
    writeln!(out, "dragstr_params_loaded\t{params_loaded}")?;
    writeln!(
        out,
        "dragstr_pair_hmm_active\t{}",
        cfg.uses_dragstr_pair_hmm()
    )?;
    writeln!(out, "calibration_ready\t{}", params_loaded)?;
    Ok(())
}
