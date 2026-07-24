//! Fixed-width traversal tiles over closed genomic intervals.
//! GATK’s adaptive active regions come later; this only splits each closed interval into
//! non-overlapping tiles for engine wiring and logging.

/// Default tile size (bp) for scaffolding traversal; real assembly regions will be adaptive.
pub const DEFAULT_TRAVERSAL_TILE_BP: u64 = 100;

/// One contiguous 1-based inclusive tile on a contig.
/// # Invariants
/// `start` / `end` are **1-based inclusive**; `end >= start` for non-empty tiles.
/// Tiles from [`tile_closed_interval`] are non-overlapping and cover the closed interval.
/// # Ownership
/// Owns contig name; coordinates are plain `u64`.
/// # Mutation
/// Immutable scaffold record; engine traversal does not mutate tiles in place.
/// # Biological assumptions
/// Fixed-width scaffold only — not adaptive GATK active regions.
/// # Java equivalence
/// Rust-native Phase 2 scaffold; adaptive `AssemblyRegion` discovery comes later.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TraversalTile {
    pub contig: String,
    pub start: u64,
    pub end: u64,
}

/// Split `[start, end]` into non-overlapping tiles of at most `tile_bp` bases.
pub fn tile_closed_interval(
    contig: &str,
    start: u64,
    end: u64,
    tile_bp: u64,
) -> Vec<TraversalTile> {
    if tile_bp == 0 {
        return Vec::new();
    }
    let mut out = Vec::new();
    let mut pos = start;
    while pos <= end {
        let region_end = (pos + tile_bp - 1).min(end);
        out.push(TraversalTile {
            contig: contig.to_string(),
            start: pos,
            end: region_end,
        });
        pos = region_end.saturating_add(1);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tiles_short_interval() {
        let t = tile_closed_interval("chr1", 1, 32, 100);
        assert_eq!(t.len(), 1);
        assert_eq!(t[0].start, 1);
        assert_eq!(t[0].end, 32);
    }

    #[test]
    fn tiles_span_multiple() {
        let t = tile_closed_interval("chr1", 1, 250, 100);
        assert_eq!(t.len(), 3);
        assert_eq!(t[0].end, 100);
        assert_eq!(t[1].start, 101);
        assert_eq!(t[2].end, 250);
    }
}
