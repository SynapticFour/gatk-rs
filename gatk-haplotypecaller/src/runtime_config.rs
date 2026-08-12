//! Central runtime configuration parsed from environment variables (Sprint I-4).
//! Production algorithm modules must not call `std::env::var` directly — use this
//! module or [`crate::parity_harness`] (harness-gated flags only).
//! Inventory: `docs/ARCHITECTURE.md`.

/// Execution / CLI surface flags (always readable in production builds).
/// # Invariants
/// `legacy_provisional` is rejected at runtime when set (Sprint B removal).
/// `scaffold_output` and `activate_output_opt_out` are independent emit toggles.
/// # Ownership
/// [`Clone`] sub-config owned by [`RuntimeConfig`]; parsed once at process edges.
/// # Mutation
/// Immutable after [`RuntimeConfig::from_env`]; must not be mutated mid-run for parity.
/// # Biological assumptions
/// None documented (process control flags only).
/// # Java equivalence
/// Rust-native env surface; no direct Java argument collection (see `ENV_INVENTORY.md`).
#[derive(Debug, Clone, Default)]
pub struct ExecutionConfig {
    /// `GATK_RS_HC_SCAFFOLD_OUTPUT=1` — header-only VCF.
    pub scaffold_output: bool,
    /// `GATK_RS_HC_ACTIVATE_OUTPUT=0` — explicit opt-out of variant emit.
    pub activate_output_opt_out: bool,
    /// `GATK_RS_HC_LEGACY_PROVISIONAL=1` — rejected (Sprint B removal).
    pub legacy_provisional: bool,
    /// `GATK_RS_HC_SEQUENTIAL=1` — force assembly-region apply batch size 1 (Peak-RSS).
    pub sequential_regions: bool,
    /// `GATK_RS_HC_LARGE_REGION_READS` — flush a region alone at/above this read count.
    /// `None` → caller default (`run.rs` `LARGE_REGION_READS_SEQUENTIAL_DEFAULT`).
    pub large_region_reads: Option<usize>,
}

/// Diagnostics that must never change emit/genotype results.
/// # Invariants
/// Debug flags affect logging/tracing only; documented as non-behavior-changing.
/// `debug_tile_overlaps` may be O(tiles×BAM) and is for scaffold diagnostics only.
/// # Ownership
/// [`Clone`] sub-config owned by [`RuntimeConfig`].
/// # Mutation
/// Immutable after env parse; helpers may re-read env on each call for dump-only paths.
/// # Biological assumptions
/// None documented.
/// # Java equivalence
/// Rust-native diagnostics env vars; not mirrored in GATK CLI.
#[derive(Debug, Clone, Default)]
pub struct DebugConfig {
    /// `GATK_RS_STRICT_CLUSTER_DEBUG=1` — stderr tracing for cluster materialize.
    pub strict_cluster_debug: bool,
    /// `GATK_RS_HC_DEBUG_TILE_OVERLAPS=1` — O(tiles×BAM) scaffold overlap count (R4-1).
    pub debug_tile_overlaps: bool,
    /// `GATK_RS_SEMANTIC_TRACE=<path>` — observe-only NDJSON semantic checkpoints.
    /// Never changes genotype/emit results; see `semantic_trace` module.
    pub semantic_trace_path: Option<String>,
    /// `GATK_RS_HC_RSS_TRACE=1` — per-region Peak-RSS diagnostics (logging only).
    pub rss_trace: bool,
}

/// Full runtime config. Prefer [`RuntimeConfig::from_env`] at process edges.
/// # Invariants
/// Production algorithm modules must not call `std::env::var` directly (use this module).
/// Harness-gated flags remain in [`crate::parity_harness`], not here.
/// # Ownership
/// Owns [`ExecutionConfig`] and [`DebugConfig`]; cheap to clone at startup.
/// # Mutation
/// Immutable after construction except dump helpers that re-parse env for debug booleans.
/// # Biological assumptions
/// None documented.
/// # Java equivalence
/// Rust-native centralized env config (Sprint I-4); GATK uses CLI args instead.
#[derive(Debug, Clone, Default)]
pub struct RuntimeConfig {
    pub execution: ExecutionConfig,
    pub debug: DebugConfig,
}

impl RuntimeConfig {
    /// Parse production-facing env vars once. Harness flags stay in [`crate::parity_harness`].
    pub fn from_env() -> Self {
        Self {
            execution: ExecutionConfig {
                scaffold_output: env_truthy("GATK_RS_HC_SCAFFOLD_OUTPUT"),
                activate_output_opt_out: env_is("GATK_RS_HC_ACTIVATE_OUTPUT", &["0", "false"]),
                legacy_provisional: env_truthy("GATK_RS_HC_LEGACY_PROVISIONAL"),
                sequential_regions: env_truthy("GATK_RS_HC_SEQUENTIAL"),
                large_region_reads: std::env::var("GATK_RS_HC_LARGE_REGION_READS")
                    .ok()
                    .and_then(|s| s.parse().ok()),
            },
            debug: DebugConfig {
                strict_cluster_debug: env_truthy("GATK_RS_STRICT_CLUSTER_DEBUG"),
                debug_tile_overlaps: env_truthy("GATK_RS_HC_DEBUG_TILE_OVERLAPS"),
                semantic_trace_path: std::env::var("GATK_RS_SEMANTIC_TRACE")
                    .ok()
                    .filter(|p| !p.is_empty() && p != "0" && !p.eq_ignore_ascii_case("off")),
                rss_trace: env_truthy("GATK_RS_HC_RSS_TRACE"),
            },
        }
    }
}

fn env_truthy(name: &str) -> bool {
    std::env::var(name)
        .map(|v| matches!(v.as_str(), "1" | "true" | "TRUE" | "yes" | "YES"))
        .unwrap_or(false)
}

fn env_is(name: &str, values: &[&str]) -> bool {
    std::env::var(name).ok().is_some_and(|v| {
        values
            .iter()
            .any(|want| v == *want || v.eq_ignore_ascii_case(want))
    })
}

/// Process-wide debug helpers (cheap; reads env each call — dump/debug only).
pub fn strict_cluster_debug_enabled() -> bool {
    RuntimeConfig::from_env().debug.strict_cluster_debug
}

/// `GATK_RS_HC_SEQUENTIAL=1` — force one assembly region at a time (16 GiB hosts).
pub fn hc_force_sequential_regions() -> bool {
    RuntimeConfig::from_env().execution.sequential_regions
}

/// Region read-count at which apply flushes alone. Override via `GATK_RS_HC_LARGE_REGION_READS`.
pub fn large_region_reads_sequential(default: usize) -> usize {
    RuntimeConfig::from_env()
        .execution
        .large_region_reads
        .unwrap_or(default)
}

/// `GATK_RS_HC_RSS_TRACE=1` — log per-region RSS (observe-only; no genotype effect).
pub fn hc_rss_trace_enabled() -> bool {
    RuntimeConfig::from_env().debug.rss_trace
}

use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Mutex, Once, OnceLock};

/// Optional soft Peak abort during k-best: `GATK_RS_HC_RSS_ABORT_MIB=<MiB>`.
///
/// When set and current RSS ≥ threshold, k-best returns partial paths instead of
/// expanding into a hard runner OOM. Unset = off (algorithm default). GIAB CI sets
/// this below hosted-runner RAM so dense shards soft-land.
pub fn hc_rss_abort_mib() -> Option<f64> {
    static CACHED: OnceLock<Option<f64>> = OnceLock::new();
    let parsed = *CACHED.get_or_init(|| {
        let parsed = std::env::var("GATK_RS_HC_RSS_ABORT_MIB")
            .ok()
            .and_then(|s| s.parse::<f64>().ok())
            .filter(|&v| v.is_finite() && v > 0.0);
        // Always announce once so CI logs prove the in-process limit (vs env-only in the job).
        // Do not start the sampler here — that would deadlock if the thread re-entered
        // this OnceLock while init is still running.
        match parsed {
            Some(limit) => {
                let rss = current_rss_mib()
                    .map(|v| format!("{v:.1}"))
                    .unwrap_or_else(|| "?".into());
                eprintln!("HC_RSS_ABORT_CONFIG limit_MiB={limit:.0} rss_MiB={rss}");
                let _ = std::io::Write::flush(&mut std::io::stderr());
            }
            None => {
                if std::env::var_os("GATK_RS_HC_RSS_ABORT_MIB").is_some() {
                    eprintln!("HC_RSS_ABORT_CONFIG ignored (unparseable or non-positive)");
                    let _ = std::io::Write::flush(&mut std::io::stderr());
                }
            }
        }
        parsed
    });
    if parsed.is_some() {
        // Watchdog samples even when TRACE is off — proves locus/phase if k-best
        // abort checks are not on the allocating path.
        ensure_rss_sampler();
    }
    parsed
}

/// Peak / phase diagnostics active when TRACE is on **or** an RSS abort limit is configured.
pub fn hc_rss_diagnostics_enabled() -> bool {
    hc_rss_trace_enabled() || hc_rss_abort_mib().is_some()
}

/// True when RSS abort is configured and current RSS is at/above the threshold.
pub fn hc_rss_abort_triggered() -> bool {
    let Some(limit) = hc_rss_abort_mib() else {
        return false;
    };
    let Some(rss) = current_rss_mib() else {
        return false;
    };
    if rss >= limit {
        static LOGGED: AtomicBool = AtomicBool::new(false);
        if !LOGGED.swap(true, Ordering::Relaxed) {
            let locus = rss_trace_locus_lock()
                .lock()
                .map(|s| s.clone())
                .unwrap_or_default();
            // Always log once (not TRACE-only) so GIAB CI shows why a shard soft-landed.
            eprintln!("HC_RSS_ABORT rss_MiB={rss:.1} limit_MiB={limit:.0} locus={locus}");
            let _ = std::io::Write::flush(&mut std::io::stderr());
        }
        true
    } else {
        false
    }
}

static RSS_TRACE_LOCUS: OnceLock<Mutex<String>> = OnceLock::new();
static RSS_SAMPLER_STARTED: Once = Once::new();
static RSS_SAMPLER_STOP: AtomicBool = AtomicBool::new(false);
/// High-water RSS (MiB × 10) seen by the background sampler while TRACE is on.
static RSS_TRACE_PEAK_X10: AtomicU64 = AtomicU64::new(0);

fn rss_trace_locus_lock() -> &'static Mutex<String> {
    RSS_TRACE_LOCUS.get_or_init(|| Mutex::new(String::new()))
}

/// Set the in-flight locus label for mid-phase RSS samples (observe-only).
pub fn rss_trace_set_locus(contig: &str, start: u64, end: u64, detail: &str) {
    if !hc_rss_diagnostics_enabled() {
        return;
    }
    if let Ok(mut s) = rss_trace_locus_lock().lock() {
        s.clear();
        use std::fmt::Write as _;
        let _ = write!(s, "{contig}:{start}-{end}");
        if !detail.is_empty() {
            let _ = write!(s, " {detail}");
        }
    }
    ensure_rss_sampler();
}

/// Clear the in-flight locus label after a region finishes.
pub fn rss_trace_clear_locus() {
    if !hc_rss_diagnostics_enabled() {
        return;
    }
    if let Ok(mut s) = rss_trace_locus_lock().lock() {
        s.clear();
    }
}

/// Log a named phase RSS sample when TRACE is on or an RSS abort limit is set.
pub fn rss_trace_checkpoint(phase: &str, detail: &str) {
    if !hc_rss_diagnostics_enabled() {
        return;
    }
    ensure_rss_sampler();
    let rss = current_rss_mib()
        .map(|v| format!("{v:.1}"))
        .unwrap_or_else(|| "?".into());
    let locus = rss_trace_locus_lock()
        .lock()
        .map(|s| s.clone())
        .unwrap_or_default();
    if detail.is_empty() {
        eprintln!("HC_RSS_TRACE phase={phase} locus={locus} rss_MiB={rss}");
    } else {
        eprintln!("HC_RSS_TRACE phase={phase} locus={locus} {detail} rss_MiB={rss}");
    }
    let _ = std::io::Write::flush(&mut std::io::stderr());
}

fn ensure_rss_sampler() {
    RSS_SAMPLER_STARTED.call_once(|| {
        RSS_SAMPLER_STOP.store(false, Ordering::Relaxed);
        std::thread::Builder::new()
            .name("hc-rss-trace".into())
            .spawn(|| {
                let mut last_logged_mib = 0.0f64;
                let mut abort_watchdog_logged = false;
                while !RSS_SAMPLER_STOP.load(Ordering::Relaxed) {
                    std::thread::sleep(std::time::Duration::from_millis(100));
                    let Some(rss) = current_rss_mib() else {
                        continue;
                    };
                    let x10 = (rss * 10.0) as u64;
                    RSS_TRACE_PEAK_X10.fetch_max(x10, Ordering::Relaxed);
                    let locus = rss_trace_locus_lock()
                        .lock()
                        .map(|s| s.clone())
                        .unwrap_or_default();
                    // Watchdog: prove over-limit even if no abort check is on the hot path.
                    if let Some(limit) = hc_rss_abort_mib() {
                        if rss >= limit && !abort_watchdog_logged {
                            abort_watchdog_logged = true;
                            eprintln!(
                                "HC_RSS_ABORT_WATCHDOG rss_MiB={rss:.1} limit_MiB={limit:.0} locus={locus}"
                            );
                            let _ = std::io::Write::flush(&mut std::io::stderr());
                        }
                    }
                    // Log when crossing ~800 MiB or every ~200 MiB thereafter.
                    if rss >= 800.0 && (last_logged_mib < 800.0 || rss >= last_logged_mib + 200.0) {
                        eprintln!(
                            "HC_RSS_TRACE sample locus={locus} rss_MiB={rss:.1} peak_MiB={:.1}",
                            x10 as f64 / 10.0
                        );
                        let _ = std::io::Write::flush(&mut std::io::stderr());
                        last_logged_mib = rss;
                    }
                }
            })
            .ok();
    });
}

/// Best-effort current process RSS in MiB (macOS/Linux). `None` if unavailable.
pub fn current_rss_mib() -> Option<f64> {
    #[cfg(target_os = "macos")]
    {
        // Prefer mach task info (no fork/`ps`); fall back to `ps`.
        if let Some(mib) = macos_task_rss_mib() {
            return Some(mib);
        }
        let out = std::process::Command::new("ps")
            .args(["-o", "rss=", "-p"])
            .arg(std::process::id().to_string())
            .output()
            .ok()?;
        let s = String::from_utf8_lossy(&out.stdout);
        let kb: f64 = s.trim().replace(',', "").parse().ok()?;
        return Some(kb / 1024.0);
    }
    #[cfg(target_os = "linux")]
    {
        let status = std::fs::read_to_string("/proc/self/status").ok()?;
        for line in status.lines() {
            if let Some(rest) = line.strip_prefix("VmRSS:") {
                let kb: f64 = rest.split_whitespace().next()?.parse().ok()?;
                return Some(kb / 1024.0);
            }
        }
        return None;
    }
    #[cfg(not(any(target_os = "macos", target_os = "linux")))]
    {
        None
    }
}

#[cfg(target_os = "macos")]
fn macos_task_rss_mib() -> Option<f64> {
    // mach_task_basic_info.resident_size is bytes.
    #[repr(C)]
    struct MachTaskBasicInfo {
        virtual_size: u64,
        resident_size: u64,
        resident_size_max: u64,
        user_time: [u32; 2],
        system_time: [u32; 2],
        policy: i32,
        suspend_count: i32,
    }
    const MACH_TASK_BASIC_INFO: u32 = 20;
    const MACH_TASK_BASIC_INFO_COUNT: u32 =
        (std::mem::size_of::<MachTaskBasicInfo>() / std::mem::size_of::<u32>()) as u32;
    extern "C" {
        fn mach_task_self() -> u32;
        fn task_info(
            target_task: u32,
            flavor: u32,
            task_info_out: *mut u32,
            task_info_outCnt: *mut u32,
        ) -> i32;
    }
    let mut info = std::mem::MaybeUninit::<MachTaskBasicInfo>::uninit();
    let mut count = MACH_TASK_BASIC_INFO_COUNT;
    let kr = unsafe {
        task_info(
            mach_task_self(),
            MACH_TASK_BASIC_INFO,
            info.as_mut_ptr() as *mut u32,
            &mut count,
        )
    };
    if kr != 0 {
        return None;
    }
    let info = unsafe { info.assume_init() };
    Some(info.resident_size as f64 / (1024.0 * 1024.0))
}
