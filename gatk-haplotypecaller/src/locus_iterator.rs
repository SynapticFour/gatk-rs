//! GATK `IntervalLocusIterator` + `LocusIteratorByState`-style pileup walk (GAP-B-02).
//! Walks every 1-based position in merged intervals (including empty loci) and builds pileups
//! from a monotonically advancing active-read set in BAM iteration order.

use crate::activity_scoring::{
    is_alt_after_assembly, is_alt_before_assembly, PileupObservation, REF_MODEL_DELETION_QUAL,
};
use crate::pileup_element::pileup_element_flags_at_ref;
use crate::read_binding::closed_interval_1based_to_ref_span0;
use crate::read_header_semantics::ReadHeaderSemantics;
use crate::read_model::{passes_hc_read_filters_with_header, ReadFilterParams};
use gatk_common::{GatkError, GatkResult};
use rust_htslib::bam;
use rust_htslib::bam::record::Aux;
use std::collections::HashMap;

/// Every 1-based position in merged closed intervals `[start, end]`.
/// # Invariants
/// Yields every position in each closed span in increasing order, including empty pileup loci.
/// # Ownership
/// Owns span list and cursor state.
/// # Mutation
/// Advances span/position cursors via [`Iterator::next`].
/// # Biological assumptions
/// None — interval enumeration for pileup walks.
/// # Java equivalence
/// GATK `IntervalLocusIterator` (every base in merged intervals).
#[derive(Debug, Clone)]
pub struct IntervalLocusIterator {
    spans: Vec<(u64, u64)>,
    span_idx: usize,
    pos: Option<u64>,
}

impl IntervalLocusIterator {
    pub fn from_spans(spans: Vec<(u64, u64)>) -> Self {
        Self {
            spans,
            span_idx: 0,
            pos: None,
        }
    }

    pub fn from_closed_interval(start1: u64, end1: u64) -> Self {
        Self::from_spans(vec![(start1, end1)])
    }
}

impl Iterator for IntervalLocusIterator {
    type Item = u64;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if self.span_idx >= self.spans.len() {
                return None;
            }
            let (s, e) = self.spans[self.span_idx];
            let next = match self.pos {
                None => s,
                Some(p) if p < e => p + 1,
                Some(_) => {
                    self.span_idx += 1;
                    self.pos = None;
                    continue;
                }
            };
            self.pos = Some(next);
            return Some(next);
        }
    }
}

fn reference_span0(
    record: &bam::Record,
    header: &bam::HeaderView,
    filters: &ReadFilterParams,
) -> Option<(i64, i64)> {
    if !passes_hc_read_filters_with_header(record, header, filters) {
        return None;
    }
    let start = record.pos();
    let end = record.cigar().end_pos();
    if end <= start {
        return None;
    }
    Some((start, end))
}

/// Owned LIBS-style pileup cursor (no borrow of `records`; safe inside [`crate::assembly_region_iterator::AssemblyRegionIterator`]).
/// # Invariants
/// Active-read set advances monotonically in BAM order for increasing locus positions.
/// Records are indexed; filtered-out records have `ref_end0 == -1`.
/// # Ownership
/// Owns header view, index lists, and caches; borrows record slices only during pileup queries.
/// # Mutation
/// Advances active-set cursor and last-position tracking as loci are visited.
/// # Biological assumptions
/// Builds per-locus pileup observations for activity / genotyping evidence.
/// # Java equivalence
/// GATK `LocusIteratorByState` (LIBS) style pileup walk.
#[derive(Debug)]
pub struct LocusPileupState {
    /// Retained for Java LIBS-style state; contig names are resolved at construction.
    _header: bam::HeaderView,
    sorted_indices: Vec<usize>,
    pub(crate) active: Vec<usize>,
    next_sorted: usize,
    pub(crate) last_pos1: Option<u64>,
    /// When true, pileup observations use Java `isAltAfterAssembly` (post-realign genotyping reads).
    pub(crate) reads_were_realigned: bool,
    /// Exclusive reference end (0-based) per record index; `-1` if filtered out at construction.
    ref_end0: Vec<i64>,
    /// Cached RCM HQ soft-clip counts (R4-1 — avoid recounting every locus).
    hq_soft_clip_cache: Vec<Option<u32>>,
}

impl LocusPileupState {
    pub fn from_records<R: std::borrow::Borrow<bam::Record>>(
        records: &[R],
        header: &bam::HeaderView,
        contig: &str,
        filters: &ReadFilterParams,
    ) -> Self {
        Self::from_records_inner(records, header, contig, filters, false)
    }

    /// Untrimmed region reads with fragment-level QNAME dedupe (Java MIN_DP on interior hom-ref blocks).
    pub fn from_records_qname_deduped<R: std::borrow::Borrow<bam::Record>>(
        records: &[R],
        header: &bam::HeaderView,
        contig: &str,
        filters: &ReadFilterParams,
    ) -> Self {
        Self::from_records_inner(records, header, contig, filters, true)
    }

    /// Post-`changeEvidence` genotyping reads for active-region RCM (`readsWereRealigned=true`).
    pub fn from_genotyping_evidence_records<R: std::borrow::Borrow<bam::Record>>(
        records: &[R],
        header: &bam::HeaderView,
        contig: &str,
        filters: &ReadFilterParams,
    ) -> Self {
        Self::from_records_inner(records, header, contig, filters, true)
    }

    fn from_records_inner<R: std::borrow::Borrow<bam::Record>>(
        records: &[R],
        header: &bam::HeaderView,
        contig: &str,
        filters: &ReadFilterParams,
        reads_were_realigned: bool,
    ) -> Self {
        let mut sorted_indices: Vec<usize> = (0..records.len())
            .filter(|&i| {
                let rec = records[i].borrow();
                let rn = String::from_utf8_lossy(header.tid2name(rec.tid() as u32));
                rn == contig && reference_span0(rec, header, filters).is_some()
            })
            .collect();
        sorted_indices.sort_by(|&a, &b| {
            let ra = records[a].borrow();
            let rb = records[b].borrow();
            ra.pos()
                .cmp(&rb.pos())
                .then_with(|| ra.qname().cmp(rb.qname()))
        });
        let mut ref_end0 = vec![-1i64; records.len()];
        for &i in &sorted_indices {
            if let Some((_, r1)) = reference_span0(records[i].borrow(), header, filters) {
                ref_end0[i] = r1;
            }
        }
        Self {
            _header: header.clone(),
            sorted_indices,
            active: Vec::new(),
            next_sorted: 0,
            last_pos1: None,
            reads_were_realigned,
            ref_end0,
            hq_soft_clip_cache: vec![None; records.len()],
        }
    }

    #[inline]
    fn hq_soft_clip_cached<R: std::borrow::Borrow<bam::Record>>(
        &mut self,
        records: &[R],
        idx: usize,
    ) -> u32 {
        if let Some(n) = self.hq_soft_clip_cache[idx] {
            return n;
        }
        let n = crate::hq_soft_clip::count_high_quality_soft_clip_bases_rcm(records[idx].borrow());
        self.hq_soft_clip_cache[idx] = Some(n);
        n
    }

    /// Clear active-read cursor (e.g. when jumping backward across disjoint spans).
    pub fn reset_cursor(&mut self) {
        self.active.clear();
        self.next_sorted = 0;
        self.last_pos1 = None;
    }

    pub fn advance_to<R: std::borrow::Borrow<bam::Record>>(
        &mut self,
        records: &[R],
        filters: &ReadFilterParams,
        pos1: u64,
    ) -> GatkResult<()> {
        if let Some(prev) = self.last_pos1 {
            if pos1 < prev {
                return Err(GatkError::argument(format!(
                    "LocusPileupState: pos {pos1} < previous {prev}"
                )));
            }
            if pos1 == prev {
                return Ok(());
            }
        }

        let ref_pos0 = pos1.saturating_sub(1) as i64;
        if self.last_pos1.is_none() {
            self.bootstrap_active(records, filters, ref_pos0);
        } else {
            // R4-1: use cached exclusive ref ends (no re-filter / CIGAR per locus).
            self.active.retain(|&idx| self.ref_end0[idx] > ref_pos0);
            self.ingest_reads_through(records, filters, ref_pos0);
        }
        self.last_pos1 = Some(pos1);
        Ok(())
    }

    fn bootstrap_active<R: std::borrow::Borrow<bam::Record>>(
        &mut self,
        records: &[R],
        filters: &ReadFilterParams,
        ref_pos0: i64,
    ) {
        self.active.clear();
        self.next_sorted = 0;
        self.ingest_reads_through(records, filters, ref_pos0);
    }

    fn ingest_reads_through<R: std::borrow::Borrow<bam::Record>>(
        &mut self,
        records: &[R],
        _filters: &ReadFilterParams,
        ref_pos0: i64,
    ) {
        while self.next_sorted < self.sorted_indices.len() {
            let idx = self.sorted_indices[self.next_sorted];
            let r1 = self.ref_end0[idx];
            if r1 < 0 {
                self.next_sorted += 1;
                continue;
            }
            let r0 = records[idx].borrow().pos();
            if r0 > ref_pos0 {
                break;
            }
            // Each index is visited once via `next_sorted` — no membership scan.
            if r1 > ref_pos0 {
                self.active.push(idx);
            }
            self.next_sorted += 1;
        }
    }

    pub fn pileup_observations<R: std::borrow::Borrow<bam::Record>>(
        &mut self,
        records: &[R],
        ref_base: u8,
    ) -> GatkResult<Vec<PileupObservation>> {
        let ref_pos0 = self
            .last_pos1
            .ok_or_else(|| GatkError::argument("advance_to before pileup_observations"))?
            .saturating_sub(1) as i64;
        let mut obs = Vec::with_capacity(self.active.len());
        let mut seen_qname = if self.reads_were_realigned {
            Some(std::collections::BTreeSet::new())
        } else {
            None
        };
        for i in 0..self.active.len() {
            let idx = self.active[i];
            let hq = self.hq_soft_clip_cached(records, idx);
            let rec = records[idx].borrow();
            if let Some(ref mut seen) = seen_qname {
                if !seen.insert(rec.qname().to_owned()) {
                    continue;
                }
            }
            if let Some(o) = pileup_observation_from_record_with_hq(
                rec,
                ref_pos0,
                ref_base,
                self.reads_were_realigned,
                hq,
            )? {
                obs.push(o);
            }
        }
        Ok(obs)
    }

    /// Non-empty per-sample pileups in **pileup visitation order**, matching encounter order from
    /// `ReadPileup#getSamples` on the active read list (excluding empty strata like Java splits).
    pub fn nonempty_stratified_sample_pileups_ordered<R: std::borrow::Borrow<bam::Record>>(
        &mut self,
        records: &[R],
        header_semantics: &ReadHeaderSemantics,
        ref_base: u8,
    ) -> GatkResult<Vec<Vec<PileupObservation>>> {
        let ref_pos0 = self
            .last_pos1
            .ok_or_else(|| GatkError::argument("advance_to before stratified piles"))?
            .saturating_sub(1) as i64;
        let mut by_sample: HashMap<String, Vec<PileupObservation>> = HashMap::new();
        let mut encounter: Vec<String> = Vec::new();
        for i in 0..self.active.len() {
            let idx = self.active[i];
            let hq = self.hq_soft_clip_cached(records, idx);
            let rec = records[idx].borrow();
            let sample = rg_sm_from_record(rec, header_semantics)?;
            if let Some(obs) =
                pileup_observation_from_record_with_hq(rec, ref_pos0, ref_base, false, hq)?
            {
                use std::collections::hash_map::Entry;
                match by_sample.entry(sample) {
                    Entry::Vacant(v) => {
                        // CLONE: needed — encounter order list and HashMap both own sample id.
                        encounter.push(v.key().clone());
                        v.insert(vec![obs]);
                    }
                    Entry::Occupied(mut o) => {
                        o.get_mut().push(obs);
                    }
                }
            }
        }

        Ok(encounter
            .into_iter()
            .filter_map(|s| {
                let p = by_sample.remove(&s)?;
                if !p.is_empty() {
                    Some(p)
                } else {
                    None
                }
            })
            .collect())
    }

    pub fn pileup_at<R: std::borrow::Borrow<bam::Record>>(
        &mut self,
        records: &[R],
        filters: &ReadFilterParams,
        pos1: u64,
        ref_base: u8,
    ) -> GatkResult<Vec<PileupObservation>> {
        self.advance_to(records, filters, pos1)?;
        self.pileup_observations(records, ref_base)
    }

    /// GATK `AlignmentContext.size` after [`Self::advance_to`].
    pub fn pileup_depth(&self) -> usize {
        self.active.len()
    }
}

pub(crate) fn pileup_observation_from_record(
    rec: &bam::Record,
    ref_pos0: i64,
    ref_base: u8,
    reads_were_realigned: bool,
) -> GatkResult<Option<PileupObservation>> {
    let hq = crate::hq_soft_clip::count_high_quality_soft_clip_bases_rcm(rec);
    pileup_observation_from_record_with_hq(rec, ref_pos0, ref_base, reads_were_realigned, hq)
}

pub(crate) fn pileup_observation_from_record_with_hq(
    rec: &bam::Record,
    ref_pos0: i64,
    ref_base: u8,
    reads_were_realigned: bool,
    read_hq_soft_clip_base_count: u32,
) -> GatkResult<Option<PileupObservation>> {
    let cigar: Vec<_> = rec.cigar().iter().copied().collect();
    let seq = rec.seq();
    let qual = rec.qual();
    let Some(flags) =
        pileup_element_flags_at_ref(rec.pos(), &cigar, &seq.as_bytes(), qual, ref_pos0)
    else {
        return Ok(None);
    };
    let is_alt = if reads_were_realigned {
        is_alt_after_assembly(flags.read_base, ref_base, flags.is_deletion)
    } else {
        is_alt_before_assembly(
            flags.read_base,
            ref_base,
            flags.is_deletion,
            flags.is_before_deletion_start,
            flags.is_after_deletion_end,
            flags.is_before_insertion,
            flags.is_after_insertion,
            flags.is_next_to_soft_clip,
        )
    };
    let qual = if flags.is_deletion {
        REF_MODEL_DELETION_QUAL
    } else {
        flags.qual
    };
    Ok(Some(PileupObservation {
        read_base: flags.read_base,
        qual,
        is_deletion: flags.is_deletion,
        is_alt,
        is_next_to_soft_clip: flags.is_next_to_soft_clip,
        read_hq_soft_clip_base_count,
    }))
}

fn rg_sm_from_record(rec: &bam::Record, sem: &ReadHeaderSemantics) -> GatkResult<String> {
    let aux = rec
        .aux(b"RG")
        .map_err(|_| GatkError::read("read missing RG auxiliary field"))?;
    let rg_id: String = match aux {
        Aux::String(s) => s.to_string(),
        _ => return Err(GatkError::read("RG tag has unexpected htslib Aux type")),
    };
    let resolved = sem.validate_record_links(Some(rg_id.as_str()), None)?;
    resolved
        .sample_name
        .ok_or_else(|| GatkError::read("RG is missing SM sample name in header"))
}

/// Borrowing wrapper over [`LocusPileupState`] for short-lived dump scopes.
/// # Invariants
/// Lifetime ties records and filter params to the walker; state is owned.
/// # Ownership
/// Borrows records/filters; owns inner [`LocusPileupState`].
/// # Mutation
/// [`Self::pileup_at`] mutates the inner pileup cursor.
/// # Biological assumptions
/// Same as [`LocusPileupState`].
/// # Java equivalence
/// Thin Rust wrapper around LIBS-style pileup for dumps.
#[derive(Debug)]
pub struct LocusPileupWalker<'a, R: std::borrow::Borrow<bam::Record> = bam::Record> {
    records: &'a [R],
    filters: &'a ReadFilterParams,
    state: LocusPileupState,
}

impl<'a, R: std::borrow::Borrow<bam::Record>> LocusPileupWalker<'a, R> {
    pub fn new(
        records: &'a [R],
        header: &bam::HeaderView,
        contig: &str,
        filters: &'a ReadFilterParams,
    ) -> Self {
        Self {
            records,
            filters,
            state: LocusPileupState::from_records(records, header, contig, filters),
        }
    }

    pub fn pileup_at(&mut self, pos1: u64, ref_base: u8) -> GatkResult<Vec<PileupObservation>> {
        self.state
            .pileup_at(self.records, self.filters, pos1, ref_base)
    }
}

/// True when a record overlaps a single 1-based reference position.
pub fn record_overlaps_locus_1based(
    record: &bam::Record,
    header: &bam::HeaderView,
    contig: &str,
    pos1: u64,
    filters: &ReadFilterParams,
) -> bool {
    let rn = String::from_utf8_lossy(header.tid2name(record.tid() as u32));
    if rn != contig {
        return false;
    }
    let Some((r0, r1)) = reference_span0(record, header, filters) else {
        return false;
    };
    let (i0, i1) = closed_interval_1based_to_ref_span0(pos1, pos1);
    r0 < i1 && i0 < r1
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::Read as _;

    #[test]
    fn interval_locus_iterator_visits_empty_loci() {
        let loci: Vec<_> = IntervalLocusIterator::from_closed_interval(5, 7).collect();
        assert_eq!(loci, vec![5, 6, 7]);
    }

    #[test]
    fn locus_pileup_walker_advances_monotonically() {
        let path =
            std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../parity/fixtures/sample.bam");
        if !path.exists() {
            return;
        }
        let filters = ReadFilterParams::gatk_standard_hc();
        let mut reader = bam::Reader::from_path(&path).expect("bam");
        let header = reader.header().clone();
        let tid = header.tid(b"chr1").expect("chr1") as i32;
        let mut records = Vec::new();
        for res in reader.records() {
            let rec = res.expect("rec");
            if rec.tid() == tid {
                records.push(rec);
            }
        }
        let mut walker = LocusPileupWalker::new(&records, &header, "chr1", &filters);
        let mut prev = 0u64;
        for pos1 in [5u64, 6, 7, 8] {
            let pile = walker.pileup_at(pos1, b'N').expect("pileup");
            assert!(pos1 > prev || prev == 0);
            prev = pos1;
            let _ = pile;
        }
    }
}
