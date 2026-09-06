//! HC read×haplotype likelihood scoring (, **F.4**, **F.5**, **F.6**).
//! Free functions + [`HcLikelihoodImplementation`] enum. PairHMM kernel selection is
//! [`crate::pairhmm_simd::PairHmmImpl`] (`FastestAvailable` → host SIMD when present).

use crate::pairhmm_log10::{GATK_PARITY_DEFAULT_GCP, GATK_PARITY_DEFAULT_INS_QUAL};
use crate::pairhmm_qual::cap_read_base_qualities;
use crate::pairhmm_simd::{
    resolve_pair_hmm_impl, score_read_haps_logless, PairHmmBackend, PairHmmImpl,
};
use crate::pcr_error_model::PcrErrorModel;
use gatk_common::GatkResult;
use std::cell::RefCell;

/// Per-read PairHMM input planes (BQ-capped quals + GOP/GCP). TLS avoids alloc-per-read
/// on dense regions (observable contract: same fills as fresh `Vec`s each call).
struct PairHmmReadScratch {
    capped: Vec<u8>,
    ins: Vec<u8>,
    del: Vec<u8>,
    gcp: Vec<u8>,
}

impl PairHmmReadScratch {
    fn empty() -> Self {
        Self {
            capped: Vec::new(),
            ins: Vec::new(),
            del: Vec::new(),
            gcp: Vec::new(),
        }
    }

    fn ensure(&mut self, n: usize) {
        if self.capped.len() < n {
            self.capped.resize(n, 0);
            self.ins.resize(n, 0);
            self.del.resize(n, 0);
            self.gcp.resize(n, 0);
        }
    }
}

thread_local! {
    static PAIRHMM_READ_SCRATCH: RefCell<PairHmmReadScratch> =
        RefCell::new(PairHmmReadScratch::empty());
}

/// Diagnostic dump of kernel-boundary arrays. Scoring is unchanged.
/// Set `GATK_RS_PAIRHMM_INPUT_DUMP` to a writable path (Java `--pair-hmm-results-file` layout
/// without the likelihood column).
fn dump_pairhmm_kernel_inputs_if_enabled(
    read_bases: &[u8],
    read_quals: &[u8],
    haplotype_bases: &[&[u8]],
    ins: &[u8],
    del: &[u8],
    gcp: &[u8],
) {
    use std::io::Write;
    use std::sync::{Mutex, OnceLock};
    static FILE: OnceLock<Option<Mutex<std::fs::File>>> = OnceLock::new();
    let Some(slot) = FILE.get_or_init(|| {
        let path = crate::runtime_config::pairhmm_input_dump_path()?;
        let mut f = std::fs::File::create(&path).ok()?;
        let _ = writeln!(
            f,
            "# hap-bases read-bases read-qual read-ins-qual read-del-qual gcp"
        );
        Some(Mutex::new(f))
    }) else {
        return;
    };
    let Ok(mut f) = slot.lock() else {
        return;
    };
    fn fastq(q: &[u8]) -> String {
        q.iter()
            .map(|&b| char::from(b.saturating_add(33)))
            .collect()
    }
    let read = std::str::from_utf8(read_bases).unwrap_or("");
    let bq = fastq(read_quals);
    let iq = fastq(ins);
    let dq = fastq(del);
    let gcp_s = fastq(gcp);
    let mut block = String::new();
    for hap in haplotype_bases {
        let hap_s = std::str::from_utf8(hap).unwrap_or("");
        block.push_str(&format!("{hap_s} {read} {bq} {iq} {dq} {gcp_s}\n"));
    }
    let _ = f.write_all(block.as_bytes());
}

/// GATK `PairHMMLikelihoodCalculationEngine` default (`--base-quality-score-threshold`).
pub const HC_DEFAULT_BASE_QUALITY_SCORE_THRESHOLD: u8 = 18;

/// GATK `ReadLikelihoodCalculationEngine.Implementation` slice for parity / production defaults.
/// # Invariants
/// Production `FastestAvailable` still scores via scalar `Log10PairHMM` until SIMD gates close
/// (`pair_hmm_impl` defaults to [`PairHmmImpl::FastestAvailable`] — host SIMD when present).
/// `FlowBased` is unused unless explicitly enabled.
/// # Ownership
/// [`Copy`] implementation discriminant.
/// # Mutation
/// Selected via [`HcLikelihoodEngineConfig`]; immutable during scoring.
/// # Biological assumptions
/// None beyond which PairHMM engine family scores reads.
/// # Java equivalence
/// GATK `ReadLikelihoodCalculationEngine.Implementation`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HcLikelihoodImplementation {
    /// Standard HC likelihood path (PairHMM kernel from [`HcLikelihoodEngineConfig::pair_hmm_impl`]).
    FastestAvailable,
    /// Flow-based engine — not used unless explicitly enabled.
    FlowBased,
}

/// DRAGSTR + stepwise filtering configuration (HC defaults: all off).
/// # Invariants
/// Production defaults: `FastestAvailable`, conservative PCR, BQ threshold 18, DRAGSTR off.
/// Stepwise filtering only activates with a non-flow primary engine.
/// # Ownership
/// Cloneable config nested in [`crate::engine::CallRegionArgs`].
/// # Mutation
/// Snapshot for a region; prepare-quals helpers may mutate quality buffers, not this config.
/// # Biological assumptions
/// Controls PairHMM quality capping, PCR indel model, and optional DRAGSTR/stepwise paths.
/// # Java equivalence
/// GATK `PairHMMLikelihoodCalculationEngine` argument collection slice.
#[derive(Debug, Clone)]
pub struct HcLikelihoodEngineConfig {
    pub implementation: HcLikelihoodImplementation,
    /// PairHMM kernel (`--pair-hmm`). Default: Log10 for F.1/F.2 / P12 stability.
    pub pair_hmm_impl: PairHmmImpl,
    pub pcr_error_model: PcrErrorModel,
    /// When true, a second FlowBased engine would score reads for allele filtering (Java stepwise).
    pub stepwise_filtering: bool,
    /// Parsed DRAGSTR params present (parity: path set in Java).
    pub dragstr_params_loaded: bool,
    pub dont_use_dragstr_pair_hmm: bool,
    /// GATK `base-quality-score-threshold` before PairHMM.
    pub base_quality_score_threshold: u8,
    /// GATK `disable-cap-base-qualities-to-map-quality`.
    pub disable_cap_read_qualities_to_mapq: bool,
}

impl HcLikelihoodEngineConfig {
    /// Production HC: scalar `Log10PairHMM` until SIMD GIAB gates promote the default.
    pub fn gatk_haplotype_caller_production() -> Self {
        Self::default()
    }

    pub fn with_pair_hmm_impl(mut self, pair_hmm_impl: PairHmmImpl) -> Self {
        self.pair_hmm_impl = pair_hmm_impl;
        self
    }

    pub fn resolved_pair_hmm_backend(&self) -> PairHmmBackend {
        resolve_pair_hmm_impl(self.pair_hmm_impl)
    }
}

impl Default for HcLikelihoodEngineConfig {
    fn default() -> Self {
        Self {
            implementation: HcLikelihoodImplementation::FastestAvailable,
            // Keep Log10 until unit + GIAB SIMD gates pass (plan step 4/6).
            pair_hmm_impl: PairHmmImpl::FastestAvailable,
            pcr_error_model: PcrErrorModel::Conservative,
            stepwise_filtering: false,
            dragstr_params_loaded: false,
            dont_use_dragstr_pair_hmm: false,
            base_quality_score_threshold: HC_DEFAULT_BASE_QUALITY_SCORE_THRESHOLD,
            disable_cap_read_qualities_to_mapq: false,
        }
    }
}

/// Apply GATK production BQ/MQ capping before PairHMM (matches `PairHMMLikelihoodCalculationEngine`).
pub fn prepare_read_quals_for_pairhmm(
    read_quals: &[u8],
    read_mapq: u8,
    config: &HcLikelihoodEngineConfig,
) -> Vec<u8> {
    let mut quals = read_quals.to_vec();
    prepare_read_quals_for_pairhmm_inplace(&mut quals, read_mapq, config);
    quals
}

/// In-place BQ/MQ cap (avoids an extra allocation when the caller already owns a buffer).
pub fn prepare_read_quals_for_pairhmm_inplace(
    read_quals: &mut [u8],
    read_mapq: u8,
    config: &HcLikelihoodEngineConfig,
) {
    cap_read_base_qualities(
        read_quals,
        read_mapq,
        config.base_quality_score_threshold,
        config.disable_cap_read_qualities_to_mapq,
    );
}

impl HcLikelihoodEngineConfig {
    pub fn uses_dragstr_pair_hmm(&self) -> bool {
        self.dragstr_params_loaded && !self.dont_use_dragstr_pair_hmm
    }

    pub fn filter_step_engine_active(&self) -> bool {
        self.stepwise_filtering && self.implementation != HcLikelihoodImplementation::FlowBased
    }

    /// Stable label for parity dumps / logging.
    ///
    /// Reports the **configured** PairHMM selection (host-independent). Resolved
    /// backends (AVX2 / NEON / scalar) vary by runner and must not appear in L2
    /// frozen dumps.
    pub fn primary_engine_label(&self) -> &'static str {
        if self.implementation == HcLikelihoodImplementation::FlowBased {
            return "FlowBased";
        }
        match self.pair_hmm_impl {
            PairHmmImpl::Log10PairHmm => "Log10PairHMM",
            PairHmmImpl::LoglessPairHmm => "LoglessPairHMM",
            PairHmmImpl::Simd => "Simd",
            PairHmmImpl::SimdF32 => "SimdF32",
            PairHmmImpl::Wavefront => "Wavefront",
            PairHmmImpl::FastestAvailable => "FastestAvailable",
        }
    }
}

/// GATK 4.4 `ReadUtils.getBaseInsertionQualities` / `getBaseDeletionQualities`.
///
/// `tag` is already FastQ-33→Phred. Missing tag → fill
/// [`GATK_PARITY_DEFAULT_INS_QUAL`] (Q45). Length mismatch is a Java
/// `createQualityModifiedRead` failure (no silent fill).
pub fn fill_indel_gop_from_optional_tag(dst: &mut [u8], tag: Option<&[u8]>) -> GatkResult<()> {
    match tag {
        None => {
            dst.fill(GATK_PARITY_DEFAULT_INS_QUAL);
            Ok(())
        }
        Some(q) if q.len() == dst.len() => {
            dst.copy_from_slice(q);
            Ok(())
        }
        Some(q) => Err(gatk_common::GatkError::algorithm(format!(
            "indel GOP tag length {} != read length {}",
            q.len(),
            dst.len()
        ))),
    }
}

/// Convenience wrapper around [`fill_indel_gop_from_optional_tag`].
pub fn indel_gop_from_optional_tag(tag: Option<&[u8]>, read_len: usize) -> GatkResult<Vec<u8>> {
    let mut v = vec![0u8; read_len];
    fill_indel_gop_from_optional_tag(&mut v, tag)?;
    Ok(v)
}

/// Score one read×haplotype with production `Log10PairHMM` (BQ cap + PCR + default indel/GCP).
pub fn log10_read_haplotype_likelihood(
    config: &HcLikelihoodEngineConfig,
    read_bases: &[u8],
    read_quals: &[u8],
    read_mapq: u8,
    haplotype_bases: &[u8],
) -> GatkResult<f64> {
    let scores = score_read_against_haplotypes(
        config,
        read_bases,
        read_quals,
        read_mapq,
        &[haplotype_bases],
        None,
        None,
    )?;
    Ok(scores[0])
}

/// Score one read against many haplotypes (BQ/PCR prepared once; haplotypes in parallel).
/// Result order matches `haplotype_bases` (algorithm-identical to sequential scoring).
/// Under `GATK_RS_HC_SEQUENTIAL=1`, scalar Log10/Logless use a plain iterator (no rayon) so
/// haplotype score order still matches sequential collect; SIMD path unchanged.
///
/// When already inside a Rayon worker (engine read-parallel), haplotype scoring stays
/// sequential — one parallel axis only (Peak TLS + oversubscription).
/// Dispatches on [`HcLikelihoodEngineConfig::implementation`] (FlowBased errors until enabled).
pub fn score_read_against_haplotypes(
    config: &HcLikelihoodEngineConfig,
    read_bases: &[u8],
    read_quals: &[u8],
    read_mapq: u8,
    haplotype_bases: &[&[u8]],
    insertion_gop: Option<&[u8]>,
    deletion_gop: Option<&[u8]>,
) -> GatkResult<Vec<f64>> {
    if config.implementation == HcLikelihoodImplementation::FlowBased {
        return Err(gatk_common::GatkError::algorithm(
            "FlowBased likelihood engine is not enabled in this Rust build (use standard HC)",
        ));
    }
    if haplotype_bases.is_empty() {
        return Ok(Vec::new());
    }
    let n = read_bases.len();
    if n == 0 {
        return Ok(vec![0.0; haplotype_bases.len()]);
    }
    let backend = config.resolved_pair_hmm_backend();

    PAIRHMM_READ_SCRATCH.with(|cell| {
        let mut scratch = cell.borrow_mut();
        scratch.ensure(n);
        scratch.capped[..n].copy_from_slice(read_quals);
        prepare_read_quals_for_pairhmm_inplace(&mut scratch.capped[..n], read_mapq, config);
        fill_indel_gop_from_optional_tag(&mut scratch.ins[..n], insertion_gop)?;
        fill_indel_gop_from_optional_tag(&mut scratch.del[..n], deletion_gop)?;
        scratch.gcp[..n].fill(GATK_PARITY_DEFAULT_GCP);
        if !config.uses_dragstr_pair_hmm() {
            let PairHmmReadScratch {
                ins,
                del,
                capped: _,
                gcp: _,
            } = &mut *scratch;
            crate::pcr_error_model::apply_pcr_error_model(
                read_bases,
                &mut ins[..n],
                &mut del[..n],
                config.pcr_error_model,
            );
        }
        let capped = &scratch.capped[..n];
        let ins = &scratch.ins[..n];
        let del = &scratch.del[..n];
        let gcp = &scratch.gcp[..n];
        dump_pairhmm_kernel_inputs_if_enabled(read_bases, capped, haplotype_bases, ins, del, gcp);
        match backend {
            PairHmmBackend::Log10Scalar => haplotype_bases
                .iter()
                .map(|hap| {
                    crate::pairhmm_log10::log10_pairhmm_likelihood(
                        read_bases, capped, hap, ins, del, gcp,
                    )
                })
                .collect(),
            PairHmmBackend::LoglessScalar => haplotype_bases
                .iter()
                .map(|hap| {
                    crate::pairhmm_logless::logless_pairhmm_likelihood(
                        read_bases, capped, hap, ins, del, gcp,
                    )
                })
                .collect(),
            _ => {
                score_read_haps_logless(backend, read_bases, capped, haplotype_bases, ins, del, gcp)
            }
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn production_pairhmm_caps_low_base_qualities() {
        let cfg = HcLikelihoodEngineConfig::default();
        let mut quals = vec![5u8, 30, 30, 30];
        let capped = prepare_read_quals_for_pairhmm(&quals, 60, &cfg);
        assert_eq!(capped[0], crate::pairhmm_qual::MIN_USABLE_Q_SCORE);
        cap_read_base_qualities(&mut quals, 60, cfg.base_quality_score_threshold, true);
        assert_eq!(capped, quals);
    }

    #[test]
    fn missing_indel_tag_fills_q45_and_present_tag_is_copied() {
        let mut dst = vec![0u8; 4];
        fill_indel_gop_from_optional_tag(&mut dst, None).unwrap();
        assert_eq!(dst, vec![45, 45, 45, 45]);
        fill_indel_gop_from_optional_tag(&mut dst, Some(&[30, 31, 32, 33])).unwrap();
        assert_eq!(dst, vec![30, 31, 32, 33]);
    }

    #[test]
    fn mismatched_indel_tag_length_is_an_error() {
        let mut dst = vec![0u8; 4];
        assert!(fill_indel_gop_from_optional_tag(&mut dst, Some(&[30, 31])).is_err());
    }
}
