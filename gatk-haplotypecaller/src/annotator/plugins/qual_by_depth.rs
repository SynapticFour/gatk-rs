//! GATK `QualByDepth` (production).

use crate::read_downsample::GatkJavaRng;
use std::cell::RefCell;

/// GATK `QualByDepth.MAX_QD_BEFORE_FIXING`.
pub const MAX_QD_BEFORE_FIXING: f64 = 35.0;
/// GATK `QualByDepth.IDEAL_HIGH_QD`.
pub const IDEAL_HIGH_QD: f64 = 30.0;
/// GATK `QualByDepth.JITTER_SIGMA`.
pub const JITTER_SIGMA: f64 = 3.0;

thread_local! {
    /// GATK `Utils.getRandomGenerator()` as consumed by `QualByDepth.fixTooHighQD`.
    /// Seed `47382911`. Not shared with reservoir downsampling (separate `GatkJavaRng`
    /// instances); they coincide when the reservoir never calls `nextInt`.
    static QD_RNG: RefCell<GatkJavaRng> = RefCell::new(GatkJavaRng::reset_gatk_default());
}

/// Restore the QD RNG to `Utils.GATK_RANDOM_SEED` (fresh JVM / `Utils.resetRandomGenerator`).
pub fn reset_gatk_qual_by_depth_rng() {
    QD_RNG.with(|cell| *cell.borrow_mut() = GatkJavaRng::reset_gatk_default());
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
pub fn fix_too_high_qd(qd: f64) -> f64 {
    QD_RNG.with(|cell| fix_too_high_qd_with_rng(qd, &mut cell.borrow_mut()))
}

/// GATK `QualByDepth.fixTooHighQD` with caller-owned RNG.
pub fn fix_too_high_qd_with_rng(qd: f64, rng: &mut GatkJavaRng) -> f64 {
    if qd < MAX_QD_BEFORE_FIXING {
        qd
    } else {
        IDEAL_HIGH_QD + rng.next_gaussian() * JITTER_SIGMA
    }
}
