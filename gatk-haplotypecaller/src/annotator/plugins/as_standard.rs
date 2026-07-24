//! Allele-specific standard annotations (parity scaffold, I-D02).

/// Per-alt `AS_AF` from site-level AF and alt index (biallelic scaffold).
pub fn as_af(site_af: f64, alt_index: usize) -> f64 {
    if alt_index == 0 {
        site_af
    } else {
        (1.0 - site_af).max(0.0)
    }
}

/// Per-alt `AS_QUAL` scaffold: scales site QUAL by alt AF share.
pub fn as_qual(site_qual: f64, site_af: f64, alt_index: usize) -> f64 {
    let share = as_af(site_af, alt_index);
    site_qual * share
}
