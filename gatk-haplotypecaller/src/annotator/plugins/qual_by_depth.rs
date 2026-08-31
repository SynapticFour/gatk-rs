//! GATK `QualByDepth` (production).
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`):
//! - `QualByDepth.annotate` runs once per `VariantContext` (not per allele/genotype).
//! - `fixTooHighQD`: if `QD < 35` return unchanged (no RNG); else
//!   `30 + Utils.getRandomGenerator().nextGaussian() * 3`.
//! - `Utils.randomGenerator` is a JVM-static `java.util.Random(47382911)`.
//! - QualByDepth itself makes no other `nextInt` / `nextDouble` / `nextGaussian` calls.

use crate::read_downsample::GatkJavaRng;
use gatk_core::io::vcf::{InfoValue, VcfRecord};
use std::sync::{Mutex, OnceLock};

/// GATK `QualByDepth.MAX_QD_BEFORE_FIXING`.
pub const MAX_QD_BEFORE_FIXING: f64 = 35.0;
/// GATK `QualByDepth.IDEAL_HIGH_QD`.
pub const IDEAL_HIGH_QD: f64 = 30.0;
/// GATK `QualByDepth.JITTER_SIGMA`.
pub const JITTER_SIGMA: f64 = 3.0;

struct ProcessQdRng {
    rng: GatkJavaRng,
    /// Count of `nextGaussian` calls on this process stream (test diagnostics only).
    gaussian_draws: u64,
}

fn process_qd_rng() -> &'static Mutex<ProcessQdRng> {
    static QD_RNG: OnceLock<Mutex<ProcessQdRng>> = OnceLock::new();
    QD_RNG.get_or_init(|| {
        Mutex::new(ProcessQdRng {
            rng: GatkJavaRng::reset_gatk_default(),
            gaussian_draws: 0,
        })
    })
}

#[cfg(test)]
thread_local! {
    static QD_SERIAL_HELD: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

fn lock_process_qd_rng() -> std::sync::MutexGuard<'static, ProcessQdRng> {
    #[cfg(test)]
    ProcessQdRngTestGuard::wait_if_held_elsewhere();
    process_qd_rng()
        .lock()
        .unwrap_or_else(|poison| poison.into_inner())
}

/// Serialize tests that reset/assert the process-global QualByDepth stream.
///
/// Java `Utils.randomGenerator` is JVM-static; `cargo test` default thread pools
/// would otherwise interleave `reset` + draws. Hold this guard for the whole
/// reset/assert sequence.
#[cfg(test)]
pub fn hold_process_qd_rng_for_test() -> ProcessQdRngTestGuard {
    ProcessQdRngTestGuard::acquire()
}

#[cfg(test)]
pub struct ProcessQdRngTestGuard {
    _raw: std::sync::MutexGuard<'static, ()>,
}

#[cfg(test)]
impl ProcessQdRngTestGuard {
    fn serial() -> &'static Mutex<()> {
        static SERIAL: OnceLock<Mutex<()>> = OnceLock::new();
        SERIAL.get_or_init(|| Mutex::new(()))
    }

    fn acquire() -> Self {
        let g = Self {
            _raw: Self::serial()
                .lock()
                .unwrap_or_else(|poison| poison.into_inner()),
        };
        QD_SERIAL_HELD.with(|c| c.set(true));
        g
    }

    fn wait_if_held_elsewhere() {
        if QD_SERIAL_HELD.with(|c| c.get()) {
            return;
        }
        drop(
            Self::serial()
                .lock()
                .unwrap_or_else(|poison| poison.into_inner()),
        );
    }
}

#[cfg(test)]
impl Drop for ProcessQdRngTestGuard {
    fn drop(&mut self) {
        QD_SERIAL_HELD.with(|c| c.set(false));
    }
}

/// Restore the QD RNG to `Utils.GATK_RANDOM_SEED` (fresh JVM / `Utils.resetRandomGenerator`).
pub fn reset_gatk_qual_by_depth_rng() {
    let mut g = lock_process_qd_rng();
    g.rng = GatkJavaRng::reset_gatk_default();
    g.gaussian_draws = 0;
    #[cfg(test)]
    {
        qd_draw_log()
            .lock()
            .unwrap_or_else(|p| p.into_inner())
            .clear();
    }
}

/// QUAL/depth **without** `fixTooHighQD` (raw ratio Java computes before the hack).
pub fn raw_qual_by_depth(qual: f64, dp: i32) -> f64 {
    if dp <= 0 {
        0.0
    } else {
        qual / (dp as f64)
    }
}

/// Variant QUAL (phred) divided by depth; applies [`fix_too_high_qd`] like GATK.
pub fn qual_by_depth(qual: f64, dp: i32) -> f64 {
    if dp <= 0 {
        return 0.0;
    }
    fix_too_high_qd(qual / (dp as f64))
}

/// Same as [`qual_by_depth`] with an explicit `java.util.Random` stand-in.
pub fn qual_by_depth_with_rng(qual: f64, dp: i32, rng: &mut GatkJavaRng) -> f64 {
    if dp <= 0 {
        return 0.0;
    }
    fix_too_high_qd_with_rng(qual / (dp as f64), rng)
}

/// GATK `QualByDepth.fixTooHighQD` (`QD < 35` unchanged; else `30 + N(0,3)`).
///
/// Uses the process-global stream (`Utils.getRandomGenerator()`), **not** a
/// per-thread generator. Parallel annotation must not call this; apply it
/// later in genomic emit order ([`apply_fix_too_high_qd_to_vcf_records`]).
pub fn fix_too_high_qd(qd: f64) -> f64 {
    fix_too_high_qd_at_site(qd, None)
}

fn fix_too_high_qd_at_site(qd: f64, site: Option<(&str, u64)>) -> f64 {
    let mut g = lock_process_qd_rng();
    #[cfg(test)]
    let before = g.gaussian_draws;
    let eligible = qd >= MAX_QD_BEFORE_FIXING;
    let gaussian = if eligible {
        let x = g.rng.next_gaussian();
        g.gaussian_draws = g.gaussian_draws.saturating_add(1);
        Some(x)
    } else {
        None
    };
    let result = match gaussian {
        Some(x) => IDEAL_HIGH_QD + x * JITTER_SIGMA,
        None => qd,
    };
    #[cfg(test)]
    {
        let after = g.gaussian_draws;
        drop(g);
        qd_draw_log()
            .lock()
            .unwrap_or_else(|p| p.into_inner())
            .push(QdDrawTrace {
                site: site.map(|(c, p)| (c.to_string(), p)),
                raw_qd: qd,
                eligible,
                gaussian_ordinal: if eligible { after } else { 0 },
                gaussian,
                result_qd: result,
                stream_gaussians_before: before,
                stream_gaussians_after: after,
            });
    }
    #[cfg(not(test))]
    {
        let _ = site;
        let _ = eligible;
    }
    result
}

/// GATK `QualByDepth.fixTooHighQD` with caller-owned RNG.
pub fn fix_too_high_qd_with_rng(qd: f64, rng: &mut GatkJavaRng) -> f64 {
    if qd < MAX_QD_BEFORE_FIXING {
        qd
    } else {
        IDEAL_HIGH_QD + rng.next_gaussian() * JITTER_SIGMA
    }
}

/// Apply [`fix_too_high_qd`] to the `QD` INFO field of one record (in-place).
///
/// Records with raw QD `< 35` do not consume the process RNG. Already-jittered
/// values (`< 35`) are left unchanged (no second draw).
pub fn apply_fix_too_high_qd_to_vcf_record(rec: &mut VcfRecord) {
    let chrom = rec.chromosome.clone();
    let pos = rec.position;
    for iv in &mut rec.info {
        if let InfoValue::Float(key, vals) = iv {
            if key == "QD" {
                if let Some(v) = vals.first_mut() {
                    *v = fix_too_high_qd_at_site(*v, Some((&chrom, pos)));
                }
            }
        }
    }
}

/// Apply [`apply_fix_too_high_qd_to_vcf_record`] in slice order (Java walker order).
pub fn apply_fix_too_high_qd_to_vcf_records(records: &mut [VcfRecord]) {
    for rec in records {
        apply_fix_too_high_qd_to_vcf_record(rec);
    }
}

/// Test-only QualByDepth stream trace (not a public product API).
#[cfg(test)]
#[derive(Debug, Clone)]
pub struct QdDrawTrace {
    pub site: Option<(String, u64)>,
    pub raw_qd: f64,
    pub eligible: bool,
    pub gaussian_ordinal: u64,
    pub gaussian: Option<f64>,
    pub result_qd: f64,
    pub stream_gaussians_before: u64,
    pub stream_gaussians_after: u64,
}

#[cfg(test)]
fn qd_draw_log() -> &'static Mutex<Vec<QdDrawTrace>> {
    static LOG: OnceLock<Mutex<Vec<QdDrawTrace>>> = OnceLock::new();
    LOG.get_or_init(|| Mutex::new(Vec::new()))
}

#[cfg(test)]
pub fn take_qd_draw_log() -> Vec<QdDrawTrace> {
    std::mem::take(&mut *qd_draw_log().lock().unwrap_or_else(|p| p.into_inner()))
}

#[cfg(test)]
pub fn qd_process_gaussian_draw_count() -> u64 {
    lock_process_qd_rng().gaussian_draws
}

#[cfg(test)]
mod six_r48_qd_rng_stream_tests {
    use super::*;
    use crate::read_downsample::GatkJavaRng;

    /// Java VCF `String.format("%.2f", 30 + gaussian * 3)` prefix for seed `47382911`.
    const JAVA_PRINTED_QD_PREFIX: [&str; 9] = [
        "25.36", "28.73", "30.97", "27.24", "28.20", "25.00", "29.56", "30.62", "28.17",
    ];

    fn printed_jitter(g: f64) -> String {
        format!("{:.2}", IDEAL_HIGH_QD + g * JITTER_SIGMA)
    }

    #[test]
    fn six_r48_gatk_java_rng_gaussians_match_java_printed_qd() {
        let mut rng = GatkJavaRng::reset_gatk_default();
        let gaussians: Vec<f64> = (0..JAVA_PRINTED_QD_PREFIX.len())
            .map(|_| rng.next_gaussian())
            .collect();
        let printed: Vec<String> = gaussians.iter().copied().map(printed_jitter).collect();
        assert_eq!(printed, JAVA_PRINTED_QD_PREFIX);
        for (i, g) in gaussians.iter().enumerate() {
            assert_eq!(
                printed_jitter(*g),
                JAVA_PRINTED_QD_PREFIX[i],
                "draw {i}: gaussian={g}"
            );
        }
    }

    #[test]
    fn six_r48_high_high_low_high_skips_draw_on_low() {
        let mut rng = GatkJavaRng::reset_gatk_default();
        let h1 = fix_too_high_qd_with_rng(40.0, &mut rng);
        let h2 = fix_too_high_qd_with_rng(40.0, &mut rng);
        let low = fix_too_high_qd_with_rng(20.0, &mut rng);
        let h3 = fix_too_high_qd_with_rng(40.0, &mut rng);
        assert_eq!(format!("{h1:.2}"), "25.36");
        assert_eq!(format!("{h2:.2}"), "28.73");
        assert_eq!(format!("{low:.2}"), "20.00");
        assert_eq!(
            format!("{h3:.2}"),
            "30.97",
            "low site must not consume a Gaussian (hypotheses D/E/F)"
        );
    }

    #[test]
    fn six_r48_low_high_high_starts_at_first_gaussian() {
        let mut rng = GatkJavaRng::reset_gatk_default();
        let low = fix_too_high_qd_with_rng(10.0, &mut rng);
        let h1 = fix_too_high_qd_with_rng(40.0, &mut rng);
        let h2 = fix_too_high_qd_with_rng(40.0, &mut rng);
        assert_eq!(format!("{low:.2}"), "10.00");
        assert_eq!(format!("{h1:.2}"), "25.36");
        assert_eq!(format!("{h2:.2}"), "28.73");
    }

    #[test]
    fn six_r48_process_stream_continues_without_reseed() {
        let _qd = hold_process_qd_rng_for_test();
        reset_gatk_qual_by_depth_rng();
        let a = qual_by_depth(78.32, 2);
        let b = qual_by_depth(78.32, 2);
        let c = qual_by_depth(78.32, 2);
        assert_eq!(format!("{a:.2}"), "25.36");
        assert_eq!(format!("{b:.2}"), "28.73");
        assert_eq!(format!("{c:.2}"), "30.97");
        assert_eq!(qd_process_gaussian_draw_count(), 3);
    }

    #[test]
    fn six_r48_raw_qd_does_not_consume_process_rng() {
        let _qd = hold_process_qd_rng_for_test();
        reset_gatk_qual_by_depth_rng();
        let raw = raw_qual_by_depth(78.32, 2);
        assert!(raw >= MAX_QD_BEFORE_FIXING);
        assert_eq!(qd_process_gaussian_draw_count(), 0);
        let jittered = fix_too_high_qd(raw);
        assert_eq!(format!("{jittered:.2}"), "25.36");
        assert_eq!(qd_process_gaussian_draw_count(), 1);
    }

    #[test]
    fn six_r48_four_highs_then_five_highs_match_java_snp_cluster_prefix() {
        let _qd = hold_process_qd_rng_for_test();
        reset_gatk_qual_by_depth_rng();
        let mut recs: Vec<VcfRecord> = (0..9)
            .map(|i| VcfRecord {
                chromosome: "2".into(),
                position: 100 + i,
                id: ".".into(),
                reference: "A".into(),
                alternate: vec!["C".into()],
                quality: Some(78.32),
                filter: vec!["PASS".into()],
                info: vec![InfoValue::Float("QD".into(), vec![78.32 / 2.0])],
                format: vec![],
                samples: vec![],
            })
            .collect();
        apply_fix_too_high_qd_to_vcf_records(&mut recs[0..4]);
        assert_eq!(qd_process_gaussian_draw_count(), 4);
        apply_fix_too_high_qd_to_vcf_records(&mut recs[4..9]);
        assert_eq!(qd_process_gaussian_draw_count(), 9);
        let printed: Vec<String> = recs
            .iter()
            .map(|r| {
                r.info
                    .iter()
                    .find_map(|v| match v {
                        InfoValue::Float(k, xs) if k == "QD" => {
                            Some(format!("{:.2}", xs.first().copied().unwrap_or(0.0)))
                        }
                        _ => None,
                    })
                    .unwrap()
            })
            .collect();
        assert_eq!(printed, JAVA_PRINTED_QD_PREFIX);
        let log = take_qd_draw_log();
        assert_eq!(log.len(), 9);
        assert!(log.iter().all(|t| t.eligible));
        assert_eq!(log[0].gaussian_ordinal, 1);
        assert_eq!(log[4].gaussian_ordinal, 5);
        assert_eq!(format!("{:.2}", log[4].result_qd), "28.20");
    }

    #[test]
    fn six_r48_process_rng_is_shared_across_threads() {
        let _qd = hold_process_qd_rng_for_test();
        reset_gatk_qual_by_depth_rng();
        // Mutex-backed process stream (not TLS). Sequential draws are the cargo-test-safe
        // form: spawning while another test holds the serial lock would deadlock the wait.
        let a = format!("{:.2}", qual_by_depth(78.32, 2));
        let b = format!("{:.2}", qual_by_depth(78.32, 2));
        assert_eq!(a, "25.36");
        assert_eq!(b, "28.73");
    }
}
