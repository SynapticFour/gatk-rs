//! Production HaplotypeCaller profiling (observe-only).
//!
//! Enable with `GATK_RS_HC_PROFILE=<path>` where `<path>` is a `.json` file or a
//! directory (writes `hc_profile.json` + `hc_profile.md`). When unset / `0` /
//! `off`, every helper is a no-op after a single atomic check.
//!
//! **Must never change genotype / emit results.** Timing and counters only.

mod genotype;
mod pairhmm;
mod report;
mod stages;

pub use genotype::GenotypeSiteSample;
pub use pairhmm::PairHmmCallSample;
pub use stages::Stage;

use std::cell::Cell;
use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Mutex, OnceLock};
use std::time::{Duration, Instant};

use genotype::GenotypeAgg;
use pairhmm::PairHmmAgg;
use stages::StageAgg;

static ENABLED: AtomicBool = AtomicBool::new(false);
static INIT: OnceLock<()> = OnceLock::new();
static OUT_PATH: OnceLock<PathBuf> = OnceLock::new();
static STATE: OnceLock<Mutex<ProfileState>> = OnceLock::new();

thread_local! {
    /// Nested AD pileup wall accumulated while genotype assign runs (observe-only).
    static AD_WALL_NS: Cell<u64> = const { Cell::new(0) };
    /// EventMap / hap-event rebuild wall (observe-only).
    static EVENT_REBUILD_WALL_NS: Cell<u64> = const { Cell::new(0) };
    static ALLELE_MAP_WALL_NS: Cell<u64> = const { Cell::new(0) };
    static MARGINALIZE_WALL_NS: Cell<u64> = const { Cell::new(0) };
    static GENOTYPE_ENUM_WALL_NS: Cell<u64> = const { Cell::new(0) };
    /// Active-region locus scoring wall (TLS; flushed per region to avoid mutex/locus).
    static ACTIVE_REGION_WALL_NS: Cell<u64> = const { Cell::new(0) };
    static ACTIVE_REGION_CALLS: Cell<u64> = const { Cell::new(0) };
}

struct ProfileState {
    stages: StageAgg,
    pairhmm: PairHmmAgg,
    genotype: GenotypeAgg,
    run_wall_start: Instant,
    cpu_start_ns: Option<u64>,
    rayon_threads: usize,
    regions: u64,
    /// Wall spent outside named stages (approx = run wall − Σ stage wall).
    unclassified_hint: bool,
}

fn state() -> &'static Mutex<ProfileState> {
    STATE.get_or_init(|| {
        Mutex::new(ProfileState {
            stages: StageAgg::default(),
            pairhmm: PairHmmAgg::default(),
            genotype: GenotypeAgg::default(),
            run_wall_start: Instant::now(),
            cpu_start_ns: process_cpu_time_ns(),
            rayon_threads: rayon::current_num_threads().max(1),
            regions: 0,
            unclassified_hint: true,
        })
    })
}

/// Parse env and arm profiling. Safe to call multiple times; first wins for path.
/// Prefer [`init`] with a path from [`crate::runtime_config`] (N-2 env allowlist).
pub fn init_from_runtime(cfg: &crate::runtime_config::RuntimeConfig) {
    if let Some(path) = &cfg.debug.hc_profile_path {
        init(path);
    }
}

/// Arm profiling to write JSON/MD at `path`.
pub fn init(path: &std::path::Path) {
    INIT.get_or_init(|| {
        let _ = OUT_PATH.set(path.to_path_buf());
        ENABLED.store(true, Ordering::Relaxed);
        let _ = state();
        eprintln!(
            "HC_PROFILE enabled out={}",
            OUT_PATH
                .get()
                .map(|p| p.display().to_string())
                .unwrap_or_default()
        );
        let _ = std::io::Write::flush(&mut std::io::stderr());
    });
}

#[inline]
pub fn enabled() -> bool {
    ENABLED.load(Ordering::Relaxed)
}

/// RAII wall+CPU timer for a named production stage.
pub struct StageGuard {
    stage: Stage,
    wall0: Instant,
    cpu0: Option<u64>,
}

impl StageGuard {
    pub fn new(stage: Stage) -> Option<Self> {
        if !enabled() {
            return None;
        }
        Some(Self {
            stage,
            wall0: Instant::now(),
            cpu0: process_cpu_time_ns(),
        })
    }
}

impl Drop for StageGuard {
    fn drop(&mut self) {
        if !enabled() {
            return;
        }
        let wall = self.wall0.elapsed();
        let cpu = match (self.cpu0, process_cpu_time_ns()) {
            (Some(a), Some(b)) if b >= a => Some(Duration::from_nanos(b - a)),
            _ => None,
        };
        if let Ok(mut st) = state().lock() {
            st.stages.add(self.stage, wall, cpu);
        }
    }
}

/// Begin a timed stage (no-op when profiling is off).
#[inline]
pub fn begin(stage: Stage) -> Option<StageGuard> {
    StageGuard::new(stage)
}

/// Record wall/cpu for a completed span (when RAII is awkward).
pub fn record_span(stage: Stage, wall: Duration, cpu: Option<Duration>) {
    if !enabled() {
        return;
    }
    if let Ok(mut st) = state().lock() {
        st.stages.add(stage, wall, cpu);
    }
}

/// Dual-write from [`crate::runtime_config::rss_trace_checkpoint`] deltas.
pub fn observe_trace_phase(phase: &str, delta_ms: Option<f64>) {
    if !enabled() {
        return;
    }
    let Some(ms) = delta_ms else {
        return;
    };
    if ms < 0.0 || !ms.is_finite() {
        return;
    }
    let Some(stage) = Stage::from_trace_phase(phase) else {
        return;
    };
    let wall = Duration::from_secs_f64(ms / 1000.0);
    if let Ok(mut st) = state().lock() {
        // TRACE deltas are wall-only; CPU unknown for inter-checkpoint gaps.
        st.stages.add_wall_only(stage, wall);
    }
}

pub fn record_pairhmm_call(sample: PairHmmCallSample) {
    if !enabled() {
        return;
    }
    if let Ok(mut st) = state().lock() {
        st.pairhmm.add(sample);
    }
}

/// Observe-only: record one PairHMM region score call (lengths + pack occupancy).
pub fn note_pairhmm_region<R: std::borrow::Borrow<rust_htslib::bam::Record>>(
    reads: &[R],
    hap_refs: &[&[u8]],
    wall: Duration,
) {
    if !enabled() {
        return;
    }
    let (simd_packs, prefix_reuse_haps, leftover_haps) = take_pairhmm_pack_occupancy();
    let (dp_eval_pack, dp_avoid_pack) = crate::pairhmm_simd::take_pack_dp_cell_stats();
    let mut read_lens = Vec::with_capacity(reads.len());
    let mut read_len_sum = 0u64;
    for rec in reads {
        let len = rec.borrow().seq().len() as u32;
        read_len_sum += u64::from(len);
        read_lens.push(len);
    }
    let mut hap_lens = Vec::with_capacity(hap_refs.len());
    let mut hap_len_sum = 0u64;
    for h in hap_refs {
        let len = h.len() as u32;
        hap_len_sum += u64::from(len);
        hap_lens.push(len);
    }
    let hap_len_total: u64 = hap_refs.iter().map(|h| h.len() as u64).sum();
    let dp_cells_full: u64 = read_lens
        .iter()
        .map(|&r| u64::from(r).saturating_mul(hap_len_total))
        .sum();
    let dp_cells_evaluated = if dp_eval_pack > 0 {
        dp_eval_pack
    } else {
        dp_cells_full
    };
    record_pairhmm_call(PairHmmCallSample {
        reads: reads.len() as u64,
        haplotypes: hap_refs.len() as u64,
        read_len_sum,
        hap_len_sum,
        read_lens,
        hap_lens,
        simd_packs,
        prefix_reuse_haps,
        leftover_haps,
        dp_cells_evaluated,
        dp_cells_avoided_prefix: dp_avoid_pack,
        wall_ns: wall.as_nanos() as u64,
    });
}

fn take_pairhmm_pack_occupancy() -> (u64, u64, u64) {
    #[cfg(target_arch = "aarch64")]
    {
        crate::pairhmm_simd::take_neon_pack_stats()
    }
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        crate::pairhmm_simd::take_avx2_pack_stats()
    }
    #[cfg(not(any(target_arch = "aarch64", target_arch = "x86", target_arch = "x86_64")))]
    {
        (0, 0, 0)
    }
}

pub fn record_genotype_site(sample: GenotypeSiteSample) {
    if !enabled() {
        return;
    }
    if let Ok(mut st) = state().lock() {
        st.genotype.add(sample);
    }
}

pub fn note_region_complete() {
    if !enabled() {
        return;
    }
    flush_active_region_tls();
    if let Ok(mut st) = state().lock() {
        st.regions += 1;
    }
}

/// Accumulate active-region construction wall without locking (flushed per region).
#[inline]
pub fn note_active_region_wall(d: Duration) {
    if !enabled() {
        return;
    }
    let ns = d.as_nanos() as u64;
    ACTIVE_REGION_WALL_NS.with(|c| c.set(c.get().saturating_add(ns)));
    ACTIVE_REGION_CALLS.with(|c| c.set(c.get().saturating_add(1)));
}

fn flush_active_region_tls() {
    let ns = ACTIVE_REGION_WALL_NS.with(|c| c.replace(0));
    let calls = ACTIVE_REGION_CALLS.with(|c| c.replace(0));
    if ns == 0 && calls == 0 {
        return;
    }
    if let Ok(mut st) = state().lock() {
        let s = &mut st.stages.stats_mut(Stage::ActiveRegionConstruction);
        s.calls = s.calls.saturating_add(calls.max(1));
        s.wall_ns = s.wall_ns.saturating_add(ns);
    }
}

/// Estimate a large allocation (bytes) for bandwidth accounting — observe-only.
pub fn note_alloc_bytes(stage: Stage, bytes: u64) {
    if !enabled() || bytes == 0 {
        return;
    }
    if let Ok(mut st) = state().lock() {
        st.stages.add_alloc(stage, bytes, 1);
    }
}

/// Accumulate AD/annotation wall (thread-local; drained into genotype samples).
#[inline]
pub fn note_ad_wall(d: Duration) {
    if !enabled() {
        return;
    }
    let ns = d.as_nanos() as u64;
    AD_WALL_NS.with(|c| c.set(c.get().saturating_add(ns)));
}

/// Accumulate event-rebuild wall (thread-local).
#[inline]
pub fn note_event_rebuild_wall(d: Duration) {
    if !enabled() {
        return;
    }
    let ns = d.as_nanos() as u64;
    EVENT_REBUILD_WALL_NS.with(|c| c.set(c.get().saturating_add(ns)));
}

#[inline]
pub fn note_allele_map_wall(d: Duration) {
    if !enabled() {
        return;
    }
    let ns = d.as_nanos() as u64;
    ALLELE_MAP_WALL_NS.with(|c| c.set(c.get().saturating_add(ns)));
}

#[inline]
pub fn note_marginalize_wall(d: Duration) {
    if !enabled() {
        return;
    }
    let ns = d.as_nanos() as u64;
    MARGINALIZE_WALL_NS.with(|c| c.set(c.get().saturating_add(ns)));
}

#[inline]
pub fn note_genotype_enum_wall(d: Duration) {
    if !enabled() {
        return;
    }
    let ns = d.as_nanos() as u64;
    GENOTYPE_ENUM_WALL_NS.with(|c| c.set(c.get().saturating_add(ns)));
}

/// Take-and-reset AD wall ns accumulated on this thread.
pub fn take_ad_wall_ns() -> u64 {
    let ns = AD_WALL_NS.with(|c| c.replace(0));
    if ns > 0 {
        if let Ok(mut st) = state().lock() {
            st.stages
                .add_wall_only(Stage::AdAnnotation, Duration::from_nanos(ns));
        }
    }
    ns
}

/// Take-and-reset event-rebuild wall ns on this thread.
pub fn take_event_rebuild_wall_ns() -> u64 {
    EVENT_REBUILD_WALL_NS.with(|c| c.replace(0))
}

pub fn take_allele_map_wall_ns() -> u64 {
    ALLELE_MAP_WALL_NS.with(|c| c.replace(0))
}

pub fn take_marginalize_wall_ns() -> u64 {
    MARGINALIZE_WALL_NS.with(|c| c.replace(0))
}

pub fn take_genotype_enum_wall_ns() -> u64 {
    GENOTYPE_ENUM_WALL_NS.with(|c| c.replace(0))
}

/// Reset per-thread genotype nested counters at the start of assignGenotypeLikelihoods.
pub fn reset_genotype_nested_walls() {
    if !enabled() {
        return;
    }
    let _ = take_ad_wall_ns();
    let _ = take_event_rebuild_wall_ns();
    let _ = take_allele_map_wall_ns();
    let _ = take_marginalize_wall_ns();
    let _ = take_genotype_enum_wall_ns();
}

/// Flush JSON + Markdown reports. Idempotent; safe to call from Drop/atexit paths.
pub fn flush() {
    if !enabled() {
        return;
    }
    flush_active_region_tls();
    let Some(json_path) = OUT_PATH.get() else {
        return;
    };
    let Ok(st) = state().lock() else {
        return;
    };
    let run_wall = st.run_wall_start.elapsed();
    let run_cpu = match (st.cpu_start_ns, process_cpu_time_ns()) {
        (Some(a), Some(b)) if b >= a => Some(Duration::from_nanos(b - a)),
        _ => None,
    };
    if let Err(e) = report::write_reports(
        json_path,
        &st.stages,
        &st.pairhmm,
        &st.genotype,
        run_wall,
        run_cpu,
        st.rayon_threads,
        st.regions,
    ) {
        eprintln!("HC_PROFILE flush failed: {e}");
        let _ = std::io::Write::flush(&mut std::io::stderr());
    } else {
        let md = json_path.with_extension("md");
        eprintln!(
            "HC_PROFILE wrote {} and {}",
            json_path.display(),
            md.display()
        );
        let _ = std::io::Write::flush(&mut std::io::stderr());
    }
}

/// Best-effort process CPU time (user+sys) in nanoseconds.
pub fn process_cpu_time_ns() -> Option<u64> {
    #[cfg(unix)]
    {
        #[repr(C)]
        struct TimeVal {
            tv_sec: i64,
            tv_usec: i64,
        }
        #[repr(C)]
        struct Rusage {
            ru_utime: TimeVal,
            ru_stime: TimeVal,
            _rest: [i64; 14],
        }
        extern "C" {
            fn getrusage(who: i32, usage: *mut Rusage) -> i32;
        }
        const RUSAGE_SELF: i32 = 0;
        // SAFETY: getrusage(RUSAGE_SELF) is well-defined; Rusage layout matches common Unix ABI.
        unsafe {
            let mut usage: Rusage = std::mem::zeroed();
            if getrusage(RUSAGE_SELF, &mut usage) != 0 {
                return None;
            }
            let sec = usage.ru_utime.tv_sec as i128 + usage.ru_stime.tv_sec as i128;
            let usec = usage.ru_utime.tv_usec as i128 + usage.ru_stime.tv_usec as i128;
            let ns = sec * 1_000_000_000 + usec * 1_000;
            if ns < 0 {
                None
            } else {
                Some(ns as u64)
            }
        }
    }
    #[cfg(not(unix))]
    {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Duration;

    #[test]
    fn disabled_by_default() {
        // Do not flip process env in unit tests that share the process — just
        // ensure StageGuard::new returns None when ENABLED is false.
        if !ENABLED.load(Ordering::Relaxed) {
            assert!(StageGuard::new(Stage::PairHmm).is_none());
        }
    }

    #[test]
    fn report_json_schema_keys() {
        let stages = StageAgg::default();
        let mut pairhmm = PairHmmAgg::default();
        pairhmm.add(PairHmmCallSample {
            reads: 10,
            haplotypes: 4,
            read_len_sum: 1500,
            hap_len_sum: 800,
            read_lens: vec![150; 10],
            hap_lens: vec![200; 4],
            simd_packs: 2,
            prefix_reuse_haps: 1,
            leftover_haps: 1,
            dp_cells_evaluated: 1000,
            dp_cells_avoided_prefix: 100,
            wall_ns: 1_000_000,
        });
        let mut genotype = GenotypeAgg::default();
        genotype.add(GenotypeSiteSample {
            candidate_alleles: 2,
            genotype_states: 3,
            pl_vector_len: 3,
            samples: 1,
            wall_ns: 50_000,
            ad_wall_ns: 10_000,
            event_rebuild_wall_ns: 5_000,
            allele_map_wall_ns: 2_000,
            marginalize_wall_ns: 3_000,
            genotype_enum_wall_ns: 1_000,
        });
        let json = report::build_json_for_test(
            &stages,
            &pairhmm,
            &genotype,
            Duration::from_secs(1),
            Some(Duration::from_secs(2)),
            2,
            3,
        );
        assert!(json.contains("\"schema\": \"gatk-rs.hc_profile.v1\""));
        assert!(json.contains("\"pairhmm\""));
        assert!(json.contains("\"genotype\""));
        assert!(json.contains("simd_pack_occupancy_pct"));
        assert!(json.contains("time_per_genotype_state_s"));
    }
}
