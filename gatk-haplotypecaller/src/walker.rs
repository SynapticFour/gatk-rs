//! GATK4 `AssemblyRegionWalker`–style read sharding. See `docs/ARCHITECTURE.md`.
//! Mirrors the first step of `AssemblyRegionWalker.makeReadShards`:
//! group user intervals **by contig**, merge overlaps on the input grid, then expand each span by
//! `assemblyRegionPadding` (default **100** bp per side, matching
//! `AssemblyRegionArgumentCollection.DEFAULT_ASSEMBLY_REGION_PADDING` in GATK 4.x).
//! **B.2:** region stream over each padded span lives in [`crate::assembly_region_iterator`] (`AssemblyRegionIterator`).
//! **B.3:** `apply` vs `callRegion` disposition — [`crate::walker_apply`] and [`crate::HaplotypeCallerEngine::walker_apply_stats`].

use gatk_common::GatkResult;
use gatk_core::reference::{IntervalSpec, SequenceDictionary};
use std::collections::HashMap;

/// GATK 4 `AssemblyRegionArgumentCollection.DEFAULT_ASSEMBLY_REGION_PADDING`.
pub const GATK_DEFAULT_ASSEMBLY_REGION_PADDING: u64 = 100;

/// One contig’s traversal shard: disjoint padded spans used to bound read fetches (B.1).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReadShard {
    pub contig: String,
    /// User `-L` intervals on this contig (merged, **unpadded**) — GATK `readShard.getIntervals`.
    pub user_spans: Vec<(u64, u64)>,
    /// Disjoint 1-based inclusive intervals after padding + merge (read query bounds).
    pub padded_spans: Vec<(u64, u64)>,
}

fn merge_closed_intervals(mut v: Vec<(u64, u64)>) -> Vec<(u64, u64)> {
    if v.is_empty() {
        return v;
    }
    v.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
    let mut out: Vec<(u64, u64)> = Vec::new();
    let mut cur = v[0];
    for &(s, e) in v.iter().skip(1) {
        if s <= cur.1.saturating_add(1) {
            cur.1 = cur.1.max(e);
        } else {
            out.push(cur);
            cur = (s, e);
        }
    }
    out.push(cur);
    out
}

fn pad_interval(contig_len: u64, start: u64, end: u64, pad: u64) -> (u64, u64) {
    let lo = start.saturating_sub(pad).max(1);
    let hi = end.saturating_add(pad).min(contig_len);
    (lo, hi)
}

/// Build read shards matching GATK’s **per-contig** shard list (one shard per contig that appears
/// in `specs`), in **reference dictionary order** among those contigs.
pub fn make_read_shards(
    dictionary: &SequenceDictionary,
    specs: &[IntervalSpec],
    assembly_region_padding: u64,
) -> GatkResult<Vec<ReadShard>> {
    let mut by_contig: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
    for spec in specs {
        dictionary.validate_interval(spec)?;
        let (c, s, e) = spec.resolve_closed_ends(dictionary)?;
        by_contig.entry(c).or_default().push((s, e));
    }

    let mut shards: Vec<ReadShard> = Vec::new();
    for rec in dictionary.contig_records() {
        // Lifetime: each contig's segment list is consumed once into merge; remove
        // moves it out of the map so segs is not cloned.
        let Some(segs) = by_contig.remove(&rec.name) else {
            continue;
        };
        let user_spans = merge_closed_intervals(segs);
        let mut padded: Vec<(u64, u64)> = user_spans
            .iter()
            .copied()
            .map(|(s, e)| pad_interval(rec.length, s, e, assembly_region_padding))
            .collect();
        padded = merge_closed_intervals(padded);
        shards.push(ReadShard {
            // CLONE: needed because owned contig id for output record.
            contig: rec.name.clone(),
            user_spans,
            padded_spans: padded,
        });
    }
    Ok(shards)
}

/// Same as [`make_read_shards`] with GATK default padding.
pub fn make_read_shards_default_padding(
    dictionary: &SequenceDictionary,
    specs: &[IntervalSpec],
) -> GatkResult<Vec<ReadShard>> {
    make_read_shards(dictionary, specs, GATK_DEFAULT_ASSEMBLY_REGION_PADDING)
}

#[cfg(test)]
mod tests {
    use super::*;
    use gatk_core::reference::IntervalSpec;

    fn tiny_dict_two_contigs() -> SequenceDictionary {
        let mut d = SequenceDictionary::new();
        d.add_contig("1".to_string(), 1000);
        d.add_contig("2".to_string(), 500);
        d
    }

    /// Two intervals on contig `1` overlap after merge; padding 100 → single span clamped to [1,1000].
    #[test]
    fn fixture_merge_same_contig_then_pad() {
        let d = tiny_dict_two_contigs();
        let specs = vec![
            IntervalSpec {
                contig: "1".into(),
                start: Some(10),
                end: Some(20),
            },
            IntervalSpec {
                contig: "1".into(),
                start: Some(15),
                end: Some(25),
            },
        ];
        let shards = make_read_shards(&d, &specs, 100).unwrap();
        assert_eq!(shards.len(), 1);
        assert_eq!(shards[0].contig, "1");
        // merged input 10-25 → pad → 1..125 (len 1000)
        assert_eq!(shards[0].padded_spans, vec![(1, 125)]);
    }

    /// Disjoint intervals on one contig → two padded blocks (no merge across gap).
    #[test]
    fn fixture_disjoint_same_contig() {
        let d = tiny_dict_two_contigs();
        let specs = vec![
            IntervalSpec {
                contig: "1".into(),
                start: Some(1),
                end: Some(10),
            },
            IntervalSpec {
                contig: "1".into(),
                start: Some(50),
                end: Some(60),
            },
        ];
        let shards = make_read_shards(&d, &specs, 5).unwrap();
        assert_eq!(shards.len(), 1);
        assert_eq!(shards[0].padded_spans, vec![(1, 15), (45, 65)]);
    }

    /// Two contigs → two shards in dictionary order (`1` then `2`).
    #[test]
    fn fixture_two_contigs_two_shards() {
        let d = tiny_dict_two_contigs();
        let specs = vec![
            IntervalSpec {
                contig: "2".into(),
                start: Some(100),
                end: Some(110),
            },
            IntervalSpec {
                contig: "1".into(),
                start: Some(200),
                end: Some(210),
            },
        ];
        let shards = make_read_shards(&d, &specs, 10).unwrap();
        assert_eq!(shards.len(), 2);
        assert_eq!(shards[0].contig, "1");
        assert_eq!(shards[0].padded_spans, vec![(190, 220)]);
        assert_eq!(shards[1].contig, "2");
        assert_eq!(shards[1].padded_spans, vec![(90, 120)]);
    }

    #[test]
    fn default_padding_constant_matches_gatk4() {
        assert_eq!(GATK_DEFAULT_ASSEMBLY_REGION_PADDING, 100);
    }
}
