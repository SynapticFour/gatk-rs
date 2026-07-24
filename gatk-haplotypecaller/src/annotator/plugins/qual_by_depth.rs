//! GATK `QualByDepth` (production).

/// GATK `QualByDepth.MAX_QD_BEFORE_FIXING`.
pub const MAX_QD_BEFORE_FIXING: f64 = 35.0;
/// GATK `QualByDepth.IDEAL_HIGH_QD`.
pub const IDEAL_HIGH_QD: f64 = 30.0;

/// Variant QUAL (phred) divided by depth; applies [`fix_too_high_qd`] like GATK.
pub fn qual_by_depth(qual: f64, dp: i32) -> f64 {
    if dp <= 0 {
        return 0.0;
    }
    fix_too_high_qd(qual / (dp as f64))
}

/// GATK `QualByDepth.fixTooHighQD` without RNG (parity: zero jitter for reproducible L2).
pub fn fix_too_high_qd(qd: f64) -> f64 {
    if qd < MAX_QD_BEFORE_FIXING {
        qd
    } else {
        IDEAL_HIGH_QD
    }
}
