//! GATK `PositionalDownsampler` — cap reads per alignment start.
//! Reservoir mode matches GATK `ReservoirDownsampler` + `Utils.getRandomGenerator` (seed `47382911`).

use rust_htslib::bam;
use std::cmp::Ordering;
#[cfg(any(feature = "dev-dumps", test))]
use std::io::Write;

/// GATK `Utils.GATK_RANDOM_SEED` — `java.util.Random` used by `ReservoirDownsampler`.
pub const GATK_JAVA_RANDOM_SEED: i64 = 47_382_911;

/// Deprecated alias; use [`GATK_JAVA_RANDOM_SEED`].
pub const GATK_DOWNSAMPLE_RNG_SEED: u64 = GATK_JAVA_RANDOM_SEED as u64;

/// GATK `AssemblyRegionArgumentCollection` typical default when enabled.
pub const GATK_DEFAULT_MAX_READS_PER_ALIGNMENT_START: u32 = 50;

/// `java.util.Random` (48-bit LCG) for parity with GATK `ReservoirDownsampler`.
/// # Invariants
/// Internal seed uses OpenJDK 48-bit LCG (`0x5DEECE66D` multiplier); matches [`GATK_JAVA_RANDOM_SEED`] after reset.
/// [`Self::next_int`] rejects out-of-range draws like OpenJDK 21.
/// # Ownership
/// Owns RNG state; callers hold `&mut Self` across downsampling passes.
/// # Mutation
/// Each `next_*` call advances internal seed; [`Self::reset_gatk_default`] restores GATK seed.
/// # Biological assumptions
/// None documented (deterministic reservoir sampling for read-cap parity).
/// # Java equivalence
/// GATK `Utils.getRandomGenerator` / `java.util.Random` used by `ReservoirDownsampler`.
#[derive(Debug, Clone)]
pub struct GatkJavaRng {
    seed: u64,
}

impl GatkJavaRng {
    pub fn reset_gatk_default() -> Self {
        Self::new(GATK_JAVA_RANDOM_SEED)
    }

    fn new(seed: i64) -> Self {
        const MULTIPLIER: u64 = 0x5DEECE66D;
        const MASK: u64 = (1u64 << 48) - 1;
        Self {
            seed: ((seed as u64) ^ MULTIPLIER) & MASK,
        }
    }

    fn next(&mut self, bits: u32) -> u32 {
        const MULTIPLIER: u64 = 0x5DEECE66D;
        const ADDEND: u64 = 0xB;
        const MASK: u64 = (1u64 << 48) - 1;
        self.seed = self.seed.wrapping_mul(MULTIPLIER).wrapping_add(ADDEND) & MASK;
        (self.seed >> (48 - bits)) as u32
    }

    /// `Random.nextInt` (32 bits).
    pub fn next_int_bounded_i32(&mut self) -> i32 {
        self.next(32) as i32
    }

    /// `Random.nextInt(bound)` — OpenJDK 21 (`u - (r = u % bound) + m < 0` rejection).
    pub fn next_int(&mut self, bound: u32) -> u32 {
        assert!(bound > 0, "bound must be positive");
        let mut r = self.next(31);
        let m = bound - 1;
        if (bound & m) == 0 {
            return ((bound as u64 * r as u64) >> 31) as u32;
        }
        loop {
            let u = r;
            r = u % bound;
            if (u as i32).wrapping_sub(r as i32).wrapping_add(m as i32) >= 0 {
                break;
            }
            r = self.next(31);
        }
        r
    }
}

/// Cap on reads sharing the same alignment start (GATK positional downsampler).
/// # Invariants
/// `max_reads_per_alignment_start == 0` disables downsampling ([`Self::disabled`]).
/// `rng_seed` matches GATK Java RNG when using reservoir mode.
/// # Ownership
/// [`Copy`] config passed into downsampling functions.
/// # Mutation
/// Immutable per downsample pass; RNG state lives outside this struct.
/// # Biological assumptions
/// Limits pileup depth spikes at identical alignment starts without biasing toward higher MAPQ.
/// # Java equivalence
/// GATK `PositionalDownsampler` / `ReservoirDownsampler` + `AssemblyRegionArgumentCollection`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PositionalDownsamplerConfig {
    pub max_reads_per_alignment_start: u32,
    pub non_random_downsampling_mode: bool,
    pub rng_seed: u64,
}

impl PositionalDownsamplerConfig {
    pub fn disabled() -> Self {
        Self {
            max_reads_per_alignment_start: 0,
            non_random_downsampling_mode: true,
            rng_seed: GATK_DOWNSAMPLE_RNG_SEED,
        }
    }

    pub fn gatk_haplotype_caller_defaults() -> Self {
        Self {
            max_reads_per_alignment_start: GATK_DEFAULT_MAX_READS_PER_ALIGNMENT_START,
            non_random_downsampling_mode: false,
            rng_seed: GATK_DOWNSAMPLE_RNG_SEED,
        }
    }
}

impl Default for PositionalDownsamplerConfig {
    fn default() -> Self {
        Self::disabled()
    }
}

/// `SAMRecord.getAlignmentStart` (1-based inclusive) as used by `HcFullParityGateDump` TreeMap keys.
#[inline]
pub fn sam_alignment_start_1based_tree_key(rec: &bam::Record) -> i32 {
    if rec.tid() < 0 {
        return 0;
    }
    (rec.pos() + 1) as i32
}

/// `ReadUtils.readHasNoAssignedPosition` for rust_htslib records (SAM / `SAMRecordToGATKReadAdapter` semantics).
pub fn read_has_no_assigned_position(rec: &bam::Record, header: &bam::HeaderView) -> bool {
    if rec.tid() < 0 {
        return true;
    }
    if (rec.tid() as u32) >= header.target_count() {
        return true;
    }
    if header.tid2name(rec.tid() as u32) == b"*" {
        return true;
    }
    sam_alignment_start_1based_tree_key(rec) <= 0
}

/// `ReadUtils.getAssignedReferenceIndex` + `ReadCoordinateComparator.compareCoordinates` (coordinates only).
pub fn compare_read_coordinates_java(
    a: &bam::Record,
    b: &bam::Record,
    header: &bam::HeaderView,
) -> i32 {
    let ai = assigned_reference_index(a, header);
    let bi = assigned_reference_index(b, header);
    if ai == -1 {
        return if bi == -1 { 0 } else { 1 };
    }
    if bi == -1 {
        return -1;
    }
    let rd = ai - bi;
    if rd != 0 {
        return rd.signum();
    }
    match sam_alignment_start_1based_tree_key(a).cmp(&sam_alignment_start_1based_tree_key(b)) {
        Ordering::Less => -1,
        Ordering::Equal => 0,
        Ordering::Greater => 1,
    }
}

fn assigned_reference_index(rec: &bam::Record, header: &bam::HeaderView) -> i32 {
    if rec.tid() < 0 {
        return -1;
    }
    if (rec.tid() as u32) >= header.target_count() {
        return -1;
    }
    if header.tid2name(rec.tid() as u32) == b"*" {
        return -1;
    }
    rec.tid()
}

fn hash_qname_java(qname: &[u8]) -> i32 {
    let mut h: i32 = 0;
    for &b in qname {
        h = h.wrapping_mul(31).wrapping_add(b as i32);
    }
    h
}

/// One `ReservoirDownsampler.submit` step (Java `ReservoirDownsampler` + non-random branch).
fn reservoir_submit_one(
    buf: &mut Vec<bam::Record>,
    cap: usize,
    cfg: &PositionalDownsamplerConfig,
    rng: &mut GatkJavaRng,
    total_reads_seen: &mut u32,
    rec: bam::Record,
) {
    *total_reads_seen = total_reads_seen.saturating_add(1);
    let t = *total_reads_seen;
    if buf.len() < cap {
        buf.push(rec);
        return;
    }
    let random_slot = if cfg.non_random_downsampling_mode {
        let h = hash_qname_java(rec.qname());
        (h % t as i32).unsigned_abs()
    } else {
        rng.next_int(t)
    };
    if (random_slot as usize) < cap {
        buf[random_slot as usize] = rec;
    }
}

/// `HcFullParityGateDump.downsamplePositional` semantics for one `TreeMap` bucket (file order preserved).
/// Consumes owned records from `group` slots (taken via `Option::take`) — no `bam::Record` clones.
fn downsample_positional_bucket_java_gate(
    group: &mut [Option<bam::Record>],
    header: &bam::HeaderView,
    cap: usize,
    cfg: &PositionalDownsamplerConfig,
    rng: &mut GatkJavaRng,
) -> Vec<bam::Record> {
    let mut out = Vec::new();
    let mut reservoir: Vec<bam::Record> = Vec::new();
    let mut total_seen = 0u32;
    for slot in group.iter_mut() {
        let rec = slot.take().expect("downsample bucket record");
        if read_has_no_assigned_position(&rec, header) {
            out.push(rec);
            continue;
        }
        reservoir_submit_one(&mut reservoir, cap, cfg, rng, &mut total_seen, rec);
    }
    out.extend(reservoir);
    out
}

/// Reservoir indices (into `group`) matching GATK `ReservoirDownsampler.submit` when **all** reads are reservoir-eligible.
fn reservoir_keep_indices(
    group: &[Option<bam::Record>],
    cap: usize,
    cfg: &PositionalDownsamplerConfig,
    rng: &mut GatkJavaRng,
) -> Vec<usize> {
    let mut slots: Vec<usize> = Vec::with_capacity(cap.min(group.len()));
    let mut total_seen = 0u32;
    for (idx, slot) in group.iter().enumerate() {
        let rec = slot.as_ref().expect("reservoir group record");
        total_seen += 1;
        if slots.len() < cap {
            slots.push(idx);
            continue;
        }
        let random_slot = if cfg.non_random_downsampling_mode {
            let h = hash_qname_java(rec.qname());
            (h % total_seen as i32).unsigned_abs()
        } else {
            rng.next_int(total_seen)
        };
        if (random_slot as usize) < cap {
            slots[random_slot as usize] = idx;
        }
    }
    slots
}

/// Filter `records` in coordinate order like GATK `PositionalDownsampler` on a sorted stream.
/// With `header: Some`, enforces `ReadCoordinateComparator.compareCoordinates` monotonicity (no `cmp == 1`)
/// and uses coordinate-based run boundaries. With `None`, falls back to `pos`-only grouping (unit tests).
pub fn apply_positional_downsampler(
    records: &mut Vec<bam::Record>,
    header: Option<&bam::HeaderView>,
    cfg: &PositionalDownsamplerConfig,
    rng: &mut GatkJavaRng,
) -> Result<(), String> {
    let cap = cfg.max_reads_per_alignment_start as usize;
    if cap == 0 {
        return Ok(());
    }
    if let Some(h) = header {
        for w in records.windows(2) {
            if compare_read_coordinates_java(&w[0], &w[1], h) == 1 {
                return Err(format!(
                    "Reads must be coordinate sorted (earlier read {:?} later read {:?})",
                    String::from_utf8_lossy(w[0].qname()),
                    String::from_utf8_lossy(w[1].qname())
                ));
            }
        }
    }
    // Take ownership so keep/reservoir paths move records instead of cloning BAM payloads.
    let mut owned: Vec<Option<bam::Record>> =
        std::mem::take(records).into_iter().map(Some).collect();
    let mut out = Vec::new();
    let mut i = 0;
    while i < owned.len() {
        let mut j = i + 1;
        if let Some(h) = header {
            while j < owned.len()
                && compare_read_coordinates_java(
                    owned[i].as_ref().expect("coord record"),
                    owned[j].as_ref().expect("coord record"),
                    h,
                ) == 0
            {
                j += 1;
            }
            let group = &mut owned[i..j];
            if group
                .iter()
                .all(|r| !read_has_no_assigned_position(r.as_ref().expect("group record"), h))
            {
                if group.len() <= cap {
                    for slot in group.iter_mut() {
                        out.push(slot.take().expect("keep record"));
                    }
                } else {
                    let keep = reservoir_keep_indices(group, cap, cfg, rng);
                    for idx in keep {
                        out.push(group[idx].take().expect("reservoir keep record"));
                    }
                }
            } else {
                out.extend(downsample_positional_bucket_java_gate(
                    group, h, cap, cfg, rng,
                ));
            }
        } else {
            while j < owned.len()
                && owned[j].as_ref().expect("pos record").pos()
                    == owned[i].as_ref().expect("pos record").pos()
            {
                j += 1;
            }
            let group = &mut owned[i..j];
            if group.len() <= cap {
                for slot in group.iter_mut() {
                    out.push(slot.take().expect("keep record"));
                }
            } else {
                let keep = reservoir_keep_indices(group, cap, cfg, rng);
                for idx in keep {
                    out.push(group[idx].take().expect("reservoir keep record"));
                }
            }
        }
        i = j;
    }
    *records = out;
    Ok(())
}

/// Summary after positional downsampling: one row per distinct alignment start (0-based `alignment_start` column).
#[cfg(any(feature = "dev-dumps", test))]
pub fn dump_positional_downsample_summary_tsv(
    alignment_path: &std::path::Path,
    cap: u32,
    non_random_downsampling_mode: bool,
    out: &mut impl Write,
) -> Result<(), String> {
    use rust_htslib::bam::Read as _;
    use std::collections::BTreeMap;
    writeln!(out, "alignment_start\tkept_count\tkept_qnames").map_err(|e| format!("write: {e}"))?;
    let mut reader = bam::Reader::from_path(alignment_path).map_err(|e| format!("open: {e}"))?;
    let header = reader.header().clone();
    let mut by_key: BTreeMap<i32, Vec<bam::Record>> = BTreeMap::new();
    for res in reader.records() {
        let rec = res.map_err(|e| format!("read: {e}"))?;
        let k = sam_alignment_start_1based_tree_key(&rec);
        by_key.entry(k).or_default().push(rec);
    }
    let cfg = PositionalDownsamplerConfig {
        max_reads_per_alignment_start: cap,
        non_random_downsampling_mode,
        rng_seed: GATK_DOWNSAMPLE_RNG_SEED,
    };
    let mut rng = GatkJavaRng::reset_gatk_default();
    for (start1, group) in by_key {
        let alignment_start0 = if start1 <= 0 {
            -1
        } else {
            i64::from(start1 - 1)
        };
        let mut slots: Vec<Option<bam::Record>> = group.into_iter().map(Some).collect();
        let kept = downsample_positional_bucket_java_gate(
            &mut slots,
            &header,
            cap as usize,
            &cfg,
            &mut rng,
        );
        let mut names: Vec<String> = kept
            .iter()
            .map(|r| String::from_utf8_lossy(r.qname()).into_owned())
            .collect();
        names.sort();
        writeln!(
            out,
            "{}\t{}\t{}",
            alignment_start0,
            names.len(),
            names.join(",")
        )
        .map_err(|e| format!("write: {e}"))?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::header::HeaderRecord;
    use rust_htslib::bam::{Header, HeaderView, Record};

    fn test_header() -> HeaderView {
        let mut header = Header::new();
        header.push_record(HeaderRecord::new(b"HD").push_tag(b"VN", &"1.6"));
        header.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", &"chr1")
                .push_tag(b"LN", &100),
        );
        HeaderView::from_header(&header)
    }

    fn make_record(header: &HeaderView, qname: &[u8], pos: i64) -> Record {
        let seq = b"AAAAAAAAAA";
        let qual = [30u8; 10];
        let mut rec = Record::new();
        rec.set(qname, None, seq, &qual);
        rec.set_tid(header.tid(b"chr1").unwrap() as i32);
        rec.set_pos(pos);
        rec.set_mapq(60);
        rec
    }

    #[test]
    fn cap_zero_is_noop() {
        let hv = test_header();
        let mut recs = vec![make_record(&hv, b"a", 5), make_record(&hv, b"b", 5)];
        let before = recs.len();
        let mut rng = GatkJavaRng::reset_gatk_default();
        apply_positional_downsampler(
            &mut recs,
            None,
            &PositionalDownsamplerConfig::disabled(),
            &mut rng,
        )
        .unwrap();
        assert_eq!(recs.len(), before);
    }

    #[test]
    fn cap_one_per_alignment_start_keeps_one_read_per_pos() {
        let hv = test_header();
        let mut recs = vec![
            make_record(&hv, b"first", 5),
            make_record(&hv, b"second", 5),
            make_record(&hv, b"other", 10),
        ];
        let mut rng = GatkJavaRng::reset_gatk_default();
        apply_positional_downsampler(
            &mut recs,
            Some(&hv),
            &PositionalDownsamplerConfig {
                max_reads_per_alignment_start: 1,
                non_random_downsampling_mode: true,
                rng_seed: GATK_DOWNSAMPLE_RNG_SEED,
            },
            &mut rng,
        )
        .unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].pos(), 5);
        assert_eq!(recs[1].pos(), 10);
    }

    #[test]
    fn java_rng_next_int_300_sequence_matches_openjdk_21() {
        let mut rng = GatkJavaRng::reset_gatk_default();
        let expected = [279u32, 185, 297, 53, 215];
        for &e in &expected {
            assert_eq!(rng.next_int(300), e, "next_int(300) mismatch");
        }
    }

    #[test]
    fn reservoir_mode_keeps_cap_per_start() {
        let hv = test_header();
        let mut recs = vec![
            make_record(&hv, b"a", 5),
            make_record(&hv, b"b", 5),
            make_record(&hv, b"c", 5),
        ];
        let mut rng = GatkJavaRng::reset_gatk_default();
        apply_positional_downsampler(
            &mut recs,
            Some(&hv),
            &PositionalDownsamplerConfig {
                max_reads_per_alignment_start: 2,
                non_random_downsampling_mode: false,
                rng_seed: GATK_DOWNSAMPLE_RNG_SEED,
            },
            &mut rng,
        )
        .unwrap();
        assert_eq!(recs.len(), 2);
    }

    #[test]
    fn unmapped_star_reads_skip_reservoir_cap() {
        let hv = test_header();
        let seq = b"A";
        let qual = [30u8];
        let mut recs: Vec<Record> = (0..5)
            .map(|i| {
                let name = format!("u{i}");
                let mut rec = Record::new();
                rec.set(name.as_bytes(), None, seq, &qual);
                rec.set_tid(-1);
                rec.set_pos(-1);
                rec.set_unmapped();
                rec.set_mapq(0);
                rec
            })
            .collect();
        let mut rng = GatkJavaRng::reset_gatk_default();
        apply_positional_downsampler(
            &mut recs,
            Some(&hv),
            &PositionalDownsamplerConfig {
                max_reads_per_alignment_start: 2,
                non_random_downsampling_mode: true,
                rng_seed: GATK_DOWNSAMPLE_RNG_SEED,
            },
            &mut rng,
        )
        .unwrap();
        assert_eq!(
            recs.len(),
            5,
            "reads with no assigned position must not be capped"
        );
    }
}
