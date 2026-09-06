//! Test-only poorly-modeled / post-kernel pipeline capture (split from `engine.rs` for N-3).
//! Idle unless [`begin_poorly_modeled_observe`] / [`begin_likelihood_pipeline_observe`]
//! is called on this thread. Does not change PairHMM or QUAL arithmetic.

use crate::haplotype::Haplotype;
use crate::region_read_likelihood::RegionReadLikelihood;

/// Capture of poorly-modeled filter inputs (6R.93/6R.94 TEST-ONLY observe).
#[derive(Clone, Debug)]
pub struct PoorlyModeledObserveRow {
    pub pass: u32,
    pub qname: String,
    pub flags: u16,
    pub read_len: usize,
    pub qual_len: usize,
    pub threshold: f64,
    pub max_ll: f64,
    pub java_equiv_keep: bool,
    pub rust_keep: bool,
    pub extra_retain: bool,
    pub n_hap_cells: usize,
    pub row_index: usize,
    pub start_1based: i64,
    pub end_1based: i64,
    pub cigar: String,
    pub argmax_col: usize,
    pub argmax_fnv: u64,
    pub argmax_is_ref: bool,
    pub argmax_len: usize,
    pub n_columns: usize,
}

/// Haplotype column identity at the exact poorly-modeled filter call (6R.94).
#[derive(Clone, Debug)]
pub struct PoorlyModeledHapColumn {
    pub pass: u32,
    pub index: usize,
    pub is_reference: bool,
    pub len: usize,
    pub fnv1a: u64,
}

/// One read×haplotype cell at the poorly-modeled filter (6R.95 TEST-ONLY).
#[derive(Clone, Debug)]
pub struct PoorlyModeledObserveCell {
    pub pass: u32,
    pub qname: String,
    pub flags: u16,
    pub hap_index: usize,
    pub hap_fnv: u64,
    pub log10_likelihood: f64,
}

thread_local! {
    static POORLY_MODELED_OBSERVE_ON: std::cell::Cell<bool> = std::cell::Cell::new(false);
    static POORLY_MODELED_OBSERVE: std::cell::RefCell<Vec<PoorlyModeledObserveRow>> =
        std::cell::RefCell::new(Vec::new());
    static POORLY_MODELED_PASS: std::cell::Cell<u32> = std::cell::Cell::new(0);
    static POORLY_MODELED_PENDING_HAPS: std::cell::RefCell<Vec<(bool, usize, u64)>> =
        std::cell::RefCell::new(Vec::new());
    static POORLY_MODELED_HAPS: std::cell::RefCell<Vec<PoorlyModeledHapColumn>> =
        std::cell::RefCell::new(Vec::new());
    static POORLY_MODELED_CELLS: std::cell::RefCell<Vec<PoorlyModeledObserveCell>> =
        std::cell::RefCell::new(Vec::new());
}

fn fnv1a64_hap_bases(bases: &[u8]) -> u64 {
    let mut h = 0xcbf29ce484222325u64;
    for &b in bases {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

fn observe_cigar_string(rec: &rust_htslib::bam::Record) -> String {
    use rust_htslib::bam::record::Cigar;
    let mut out = String::new();
    for c in rec.cigar().iter() {
        let (n, op) = match c {
            Cigar::Match(n) => (*n, 'M'),
            Cigar::Ins(n) => (*n, 'I'),
            Cigar::Del(n) => (*n, 'D'),
            Cigar::SoftClip(n) => (*n, 'S'),
            Cigar::HardClip(n) => (*n, 'H'),
            Cigar::Equal(n) => (*n, '='),
            Cigar::Diff(n) => (*n, 'X'),
            Cigar::RefSkip(n) => (*n, 'N'),
            Cigar::Pad(n) => (*n, 'P'),
        };
        out.push_str(&n.to_string());
        out.push(op);
    }
    out
}

/// Enable per-read poorly-modeled input capture for this thread (TEST-ONLY dump).
pub fn begin_poorly_modeled_observe() {
    POORLY_MODELED_OBSERVE_ON.with(|c| c.set(true));
    POORLY_MODELED_OBSERVE.with(|v| v.borrow_mut().clear());
    POORLY_MODELED_PASS.with(|c| c.set(0));
    POORLY_MODELED_PENDING_HAPS.with(|v| v.borrow_mut().clear());
    POORLY_MODELED_HAPS.with(|v| v.borrow_mut().clear());
    POORLY_MODELED_CELLS.with(|v| v.borrow_mut().clear());
}

/// Snapshot haplotype columns used by the next poorly-modeled filter (TEST-ONLY).
pub fn observe_poorly_modeled_haplotypes(haps: &[Haplotype]) {
    if !POORLY_MODELED_OBSERVE_ON.with(|c| c.get()) {
        return;
    }
    POORLY_MODELED_PENDING_HAPS.with(|v| {
        *v.borrow_mut() = haps
            .iter()
            .map(|h| (h.is_reference, h.bases.len(), fnv1a64_hap_bases(&h.bases)))
            .collect();
    });
}

/// Disable capture and take accumulated rows.
pub fn take_poorly_modeled_observe() -> Vec<PoorlyModeledObserveRow> {
    POORLY_MODELED_OBSERVE_ON.with(|c| c.set(false));
    POORLY_MODELED_OBSERVE.with(|v| std::mem::take(&mut *v.borrow_mut()))
}

/// Take haplotype-column snapshots captured during observe.
pub fn take_poorly_modeled_haplotypes() -> Vec<PoorlyModeledHapColumn> {
    POORLY_MODELED_HAPS.with(|v| std::mem::take(&mut *v.borrow_mut()))
}

/// Take per-cell likelihoods captured during observe (6R.95).
pub fn take_poorly_modeled_cells() -> Vec<PoorlyModeledObserveCell> {
    POORLY_MODELED_CELLS.with(|v| std::mem::take(&mut *v.borrow_mut()))
}

pub(super) fn start_poorly_modeled_filter_pass() -> u32 {
    if !POORLY_MODELED_OBSERVE_ON.with(|c| c.get()) {
        return 0;
    }
    POORLY_MODELED_PASS.with(|c| {
        let n = c.get() + 1;
        c.set(n);
        n
    })
}

pub(super) fn record_poorly_modeled_hap_columns(observe_pass: u32) {
    if observe_pass == 0 {
        return;
    }
    let pending = POORLY_MODELED_PENDING_HAPS.with(|v| v.borrow().clone());
    POORLY_MODELED_HAPS.with(|cols| {
        let mut cols = cols.borrow_mut();
        for (index, &(is_reference, len, fnv1a)) in pending.iter().enumerate() {
            cols.push(PoorlyModeledHapColumn {
                pass: observe_pass,
                index,
                is_reference,
                len,
                fnv1a,
            });
        }
    });
}

pub(super) fn record_poorly_modeled_filter_read(
    observe_pass: u32,
    ll: &[RegionReadLikelihood],
    rec: &rust_htslib::bam::Record,
    read_idx: usize,
    best_ll: f64,
    threshold: f64,
    java_equiv_keep: bool,
    retain: bool,
) {
    if observe_pass == 0 {
        return;
    }
    let pending = POORLY_MODELED_PENDING_HAPS.with(|v| v.borrow().clone());
    let n_columns = pending.len();
    let mut n_hap_cells = 0usize;
    let mut argmax_col = 0usize;
    let mut argmax_seen = f64::NEG_INFINITY;
    for e in ll.iter().filter(|e| e.read_index.get() == read_idx) {
        n_hap_cells += 1;
        if e.log10_likelihood > argmax_seen {
            argmax_seen = e.log10_likelihood;
            argmax_col = e.haplotype_index.get();
        }
    }
    let (argmax_is_ref, argmax_len, argmax_fnv) =
        pending.get(argmax_col).copied().unwrap_or((false, 0, 0));
    let qname = String::from_utf8_lossy(rec.qname()).into_owned();
    let flags = rec.flags();
    POORLY_MODELED_CELLS.with(|cells| {
        let mut cells = cells.borrow_mut();
        for e in ll.iter().filter(|e| e.read_index.get() == read_idx) {
            let hap_index = e.haplotype_index.get();
            let hap_fnv = pending.get(hap_index).map(|p| p.2).unwrap_or(0);
            cells.push(PoorlyModeledObserveCell {
                pass: observe_pass,
                qname: qname.clone(),
                flags,
                hap_index,
                hap_fnv,
                log10_likelihood: e.log10_likelihood,
            });
        }
    });
    POORLY_MODELED_OBSERVE.with(|rows| {
        rows.borrow_mut().push(PoorlyModeledObserveRow {
            pass: observe_pass,
            qname,
            flags,
            read_len: rec.seq_len(),
            qual_len: rec.qual().len().max(1),
            threshold,
            max_ll: best_ll,
            java_equiv_keep,
            rust_keep: retain,
            extra_retain: retain && !java_equiv_keep,
            n_hap_cells,
            row_index: read_idx,
            start_1based: rec.pos() + 1,
            end_1based: i64::from(crate::read_unclip::alignment_end_1based(rec)),
            cigar: observe_cigar_string(rec),
            argmax_col,
            argmax_fnv,
            argmax_is_ref,
            argmax_len,
            n_columns,
        });
    });
}

/// One read×haplotype cell at a named post-kernel pipeline stage (6R.96 TEST-ONLY).
#[derive(Clone, Debug)]
pub struct LikelihoodPipelineCell {
    pub seq: u32,
    pub stage: &'static str,
    pub qname: String,
    pub flags: u16,
    pub read_index: usize,
    pub hap_index: usize,
    pub hap_fnv: u64,
    pub log10_likelihood: f64,
    pub n_reads: usize,
    pub n_haps: usize,
}

/// Dimensions of one captured post-kernel stage (6R.96 TEST-ONLY).
#[derive(Clone, Debug)]
pub struct LikelihoodPipelineSnap {
    pub seq: u32,
    pub stage: &'static str,
    pub n_reads: usize,
    pub n_haps: usize,
    pub n_ll_entries: usize,
}

thread_local! {
    static PIPELINE_OBSERVE_ON: std::cell::Cell<bool> = std::cell::Cell::new(false);
    static PIPELINE_SEQ: std::cell::Cell<u32> = std::cell::Cell::new(0);
    static PIPELINE_KERNEL_N: std::cell::Cell<u32> = std::cell::Cell::new(0);
    static PIPELINE_READ_ID: std::cell::RefCell<Vec<(String, u16)>> =
        std::cell::RefCell::new(Vec::new());
    static PIPELINE_CELLS: std::cell::RefCell<Vec<LikelihoodPipelineCell>> =
        std::cell::RefCell::new(Vec::new());
    static PIPELINE_SNAPS: std::cell::RefCell<Vec<LikelihoodPipelineSnap>> =
        std::cell::RefCell::new(Vec::new());
}

/// Enable post-kernel likelihood-pipeline capture for this thread (TEST-ONLY dump).
pub fn begin_likelihood_pipeline_observe() {
    PIPELINE_OBSERVE_ON.with(|c| c.set(true));
    PIPELINE_SEQ.with(|c| c.set(0));
    PIPELINE_KERNEL_N.with(|c| c.set(0));
    PIPELINE_READ_ID.with(|v| v.borrow_mut().clear());
    PIPELINE_CELLS.with(|v| v.borrow_mut().clear());
    PIPELINE_SNAPS.with(|v| v.borrow_mut().clear());
}

/// Disable capture and take accumulated post-kernel pipeline cells.
pub fn take_likelihood_pipeline_cells() -> Vec<LikelihoodPipelineCell> {
    PIPELINE_OBSERVE_ON.with(|c| c.set(false));
    PIPELINE_CELLS.with(|v| std::mem::take(&mut *v.borrow_mut()))
}

/// Take stage-dimension snapshots captured during pipeline observe.
pub fn take_likelihood_pipeline_snaps() -> Vec<LikelihoodPipelineSnap> {
    PIPELINE_SNAPS.with(|v| std::mem::take(&mut *v.borrow_mut()))
}

fn pipeline_observe_on() -> bool {
    PIPELINE_OBSERVE_ON.with(|c| c.get())
}

fn next_pipeline_seq() -> u32 {
    PIPELINE_SEQ.with(|c| {
        let n = c.get() + 1;
        c.set(n);
        n
    })
}

fn hap_fnv_list(haplotypes: &[Haplotype]) -> Vec<u64> {
    haplotypes
        .iter()
        .map(|h| fnv1a64_hap_bases(&h.bases))
        .collect()
}

fn capture_likelihood_pipeline_inner(
    stage: &'static str,
    ll: &[RegionReadLikelihood],
    haplotypes: &[Haplotype],
) {
    if !pipeline_observe_on() {
        return;
    }
    let seq = next_pipeline_seq();
    let hap_fnv = hap_fnv_list(haplotypes);
    let n_haps = haplotypes.len();
    let n_reads = PIPELINE_READ_ID.with(|ids| ids.borrow().len());
    PIPELINE_SNAPS.with(|s| {
        s.borrow_mut().push(LikelihoodPipelineSnap {
            seq,
            stage,
            n_reads,
            n_haps,
            n_ll_entries: ll.len(),
        });
    });
    PIPELINE_CELLS.with(|cells| {
        let mut cells = cells.borrow_mut();
        PIPELINE_READ_ID.with(|ids| {
            let ids = ids.borrow();
            for e in ll {
                let ri = e.read_index.get();
                let hi = e.haplotype_index.get();
                let Some((qname, flags)) = ids.get(ri) else {
                    continue;
                };
                cells.push(LikelihoodPipelineCell {
                    seq,
                    stage,
                    qname: qname.clone(),
                    flags: *flags,
                    read_index: ri,
                    hap_index: hi,
                    hap_fnv: hap_fnv.get(hi).copied().unwrap_or(0),
                    log10_likelihood: e.log10_likelihood,
                    n_reads,
                    n_haps,
                });
            }
        });
    });
}

/// Capture a scored PairHMM matrix using the exact reads that were scored (6R.96).
/// First kernel call is `post_kernel`; later kernels (refresh) are `refresh`.
pub(super) fn capture_scored_likelihood_pipeline<
    R: std::borrow::Borrow<rust_htslib::bam::Record>,
>(
    ll: &[RegionReadLikelihood],
    reads: &[R],
    haplotypes: &[Haplotype],
) {
    if !pipeline_observe_on() {
        return;
    }
    let n = PIPELINE_KERNEL_N.with(|c| {
        let n = c.get() + 1;
        c.set(n);
        n
    });
    let stage = if n == 1 { "post_kernel" } else { "refresh" };
    PIPELINE_READ_ID.with(|ids| {
        *ids.borrow_mut() = reads
            .iter()
            .map(|r| {
                let rec = r.borrow();
                (
                    String::from_utf8_lossy(rec.qname()).into_owned(),
                    rec.flags(),
                )
            })
            .collect();
    });
    capture_likelihood_pipeline_inner(stage, ll, haplotypes);
}

/// Capture a later stage that still uses the last scored read-index list (6R.96).
pub(super) fn capture_likelihood_pipeline_stage(
    stage: &'static str,
    ll: &[RegionReadLikelihood],
    haplotypes: &[Haplotype],
) {
    capture_likelihood_pipeline_inner(stage, ll, haplotypes);
}
