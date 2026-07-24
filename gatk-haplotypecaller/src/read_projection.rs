//! Read ↔ reference projection and clipping helpers.
//! Coordinates follow BAM conventions: `alignment_start` is 0-based inclusive first aligned
//! reference position; `query_index` is 0-based into the read sequence as stored in the BAM
//! record (soft-clipped bases are present; hard-clipped bases are absent).

use rust_htslib::bam::record::{Cigar, CigarString};

/// Total hard-clip bases in the CIGAR (SAM: hard clips are not stored in SEQ).
pub fn cigar_hard_clip_length(cigar: &CigarString) -> u32 {
    cigar
        .iter()
        .filter_map(|c| match c {
            Cigar::HardClip(n) => Some(*n),
            _ => None,
        })
        .sum()
}

/// Leading and trailing soft-clip lengths (typical 5'/3' clips; internal `S` is rare).
pub fn cigar_soft_clip_ends(cigar: &CigarString) -> (u32, u32) {
    let mut lead = 0u32;
    for c in cigar.iter() {
        match c {
            Cigar::SoftClip(n) => lead += n,
            Cigar::HardClip(_) => {}
            _ => break,
        }
    }
    let mut trail = 0u32;
    for c in cigar.iter().rev() {
        match c {
            Cigar::SoftClip(n) => trail += n,
            Cigar::HardClip(_) => {}
            _ => break,
        }
    }
    (lead, trail)
}

fn cigar_consumes_reference(op: &Cigar) -> u32 {
    match op {
        Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) | Cigar::Del(n) | Cigar::RefSkip(n) => {
            *n
        }
        Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) => 0,
        Cigar::Pad(n) => *n,
    }
}

fn cigar_consumes_query(op: &Cigar) -> u32 {
    match op {
        Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) | Cigar::Ins(n) | Cigar::SoftClip(n) => {
            *n
        }
        Cigar::HardClip(_) | Cigar::Del(_) | Cigar::RefSkip(_) => 0,
        Cigar::Pad(n) => *n,
    }
}

/// Map a **0-based** query index into the aligned reference position, or `None` if the index
/// falls in a soft clip / insertion-only block (no unique reference base).
pub fn reference_position_at_query_index(
    alignment_start: i64,
    cigar: &CigarString,
    query_index: usize,
) -> Option<i64> {
    let mut r = alignment_start;
    let mut q: usize = 0;
    for op in cigar.iter() {
        let dq = cigar_consumes_query(op) as usize;
        let dr = cigar_consumes_reference(op) as i64;
        if dq == 0 && dr > 0 {
            r += dr;
            continue;
        }
        if dq > 0 && dr == 0 {
            if query_index < q + dq {
                return None;
            }
            q += dq;
            continue;
        }
        if dq > 0 && dr > 0 {
            if query_index < q + dq {
                let off = (query_index - q) as i64;
                return Some(r + off);
            }
            r += dr;
            q += dq;
            continue;
        }
    }
    None
}

/// Map a **0-based inclusive** reference position to the read query index, or `None` if the
/// read does not align that reference base (deletion gap / before alignment / past end).
pub fn query_index_at_reference_position(
    alignment_start: i64,
    cigar: &CigarString,
    ref_pos: i64,
) -> Option<usize> {
    let mut r = alignment_start;
    let mut q: usize = 0;
    for op in cigar.iter() {
        let dq = cigar_consumes_query(op) as usize;
        let dr = cigar_consumes_reference(op) as i64;
        if dq == 0 && dr > 0 {
            if ref_pos >= r && ref_pos < r + dr {
                return None;
            }
            r += dr;
            continue;
        }
        if dq > 0 && dr == 0 {
            q += dq;
            continue;
        }
        if dq > 0 && dr > 0 {
            if ref_pos >= r && ref_pos < r + dr {
                let off = (ref_pos - r) as usize;
                return Some(q + off);
            }
            r += dr;
            q += dq;
            continue;
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::Cigar;

    fn cig(v: Vec<Cigar>) -> CigarString {
        CigarString(v)
    }

    #[test]
    fn soft_clip_ends_2s8m() {
        let c = cig(vec![Cigar::SoftClip(2), Cigar::Match(8)]);
        assert_eq!(cigar_soft_clip_ends(&c), (2, 0));
        assert_eq!(cigar_hard_clip_length(&c), 0);
    }

    #[test]
    fn projection_roundtrip_8m() {
        let c = cig(vec![Cigar::Match(8)]);
        let pos = 100i64;
        for qi in 0..8usize {
            let rp = reference_position_at_query_index(pos, &c, qi).unwrap();
            assert_eq!(query_index_at_reference_position(pos, &c, rp).unwrap(), qi);
        }
    }

    #[test]
    fn projection_soft_clip_offset() {
        let c = cig(vec![Cigar::SoftClip(2), Cigar::Match(8)]);
        let pos = 100i64;
        assert!(reference_position_at_query_index(pos, &c, 0).is_none());
        assert_eq!(reference_position_at_query_index(pos, &c, 2).unwrap(), 100);
        assert_eq!(query_index_at_reference_position(pos, &c, 105).unwrap(), 7);
    }

    #[test]
    fn deletion_skips_reference_bases() {
        let c = cig(vec![Cigar::Match(3), Cigar::Del(2), Cigar::Match(5)]);
        let pos = 10i64;
        assert_eq!(query_index_at_reference_position(pos, &c, 12).unwrap(), 2);
        assert!(query_index_at_reference_position(pos, &c, 13).is_none());
        assert_eq!(query_index_at_reference_position(pos, &c, 10).unwrap(), 0);
        assert_eq!(query_index_at_reference_position(pos, &c, 15).unwrap(), 3);
    }

    #[test]
    fn insertion_query_has_no_unique_ref_base() {
        let c = cig(vec![Cigar::Match(4), Cigar::Ins(2), Cigar::Match(4)]);
        let pos = 100i64;
        assert_eq!(reference_position_at_query_index(pos, &c, 3).unwrap(), 103);
        assert!(reference_position_at_query_index(pos, &c, 4).is_none());
        assert!(reference_position_at_query_index(pos, &c, 5).is_none());
        assert_eq!(reference_position_at_query_index(pos, &c, 6).unwrap(), 104);
    }

    #[test]
    fn hard_clip_length_is_counted() {
        let c = cig(vec![
            Cigar::HardClip(3),
            Cigar::SoftClip(2),
            Cigar::Match(5),
            Cigar::SoftClip(1),
            Cigar::HardClip(4),
        ]);
        assert_eq!(cigar_hard_clip_length(&c), 7);
        assert_eq!(cigar_soft_clip_ends(&c), (2, 1));
    }

    #[test]
    fn ref_pos_outside_alignment_returns_none() {
        let c = cig(vec![Cigar::Match(6)]);
        let pos = 50i64;
        assert!(query_index_at_reference_position(pos, &c, 49).is_none());
        assert!(query_index_at_reference_position(pos, &c, 56).is_none());
    }
}
