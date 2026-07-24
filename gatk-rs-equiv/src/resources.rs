//! Host-resource caps for constrained machines (e.g. MacBook Air M4 16GB).

use anyhow::{bail, Result};
use std::path::Path;
use std::process::Command;

/// Default max threads for callers / hap.py / RTG on 16GB hosts.
pub const DEFAULT_THREADS: u32 = 2;
/// Hard ceiling — refuse runaway `--threads` on laptop-class hosts.
pub const MAX_THREADS: u32 = 4;
/// Default Java HC heap (GATK is memory-hungry; leave headroom for hap.py + OS).
pub const DEFAULT_JAVA_XMX: &str = "4g";
/// Minimum free disk (GiB) before starting a full equiv run.
pub const DEFAULT_MIN_FREE_GB: u64 = 8;

pub fn clamp_threads(requested: u32) -> u32 {
    requested.clamp(1, MAX_THREADS)
}

pub fn apply_process_env(threads: u32) {
    let t = clamp_threads(threads).to_string();
    // Only set if unset — allow explicit operator override.
    if std::env::var_os("RAYON_NUM_THREADS").is_none() {
        std::env::set_var("RAYON_NUM_THREADS", &t);
    }
    if std::env::var_os("JAVA_TOOL_OPTIONS").is_none() {
        std::env::set_var(
            "JAVA_TOOL_OPTIONS",
            format!("-Xmx{DEFAULT_JAVA_XMX} -XX:+UseParallelGC"),
        );
    }
    if std::env::var_os("OMP_NUM_THREADS").is_none() {
        std::env::set_var("OMP_NUM_THREADS", &t);
    }
}

/// Best-effort free space check on the filesystem that holds `path`.
pub fn require_free_gb(path: &Path, need_gb: u64) -> Result<()> {
    let avail = free_gb(path)?;
    if avail < need_gb {
        bail!(
            "only ~{avail} GiB free near {} (need ≥{need_gb} GiB). \
             Free disk or set --min-free-gb lower (not recommended on 16GB hosts).",
            path.display()
        );
    }
    eprintln!("[gatk-rs-equiv] disk ok: ~{avail} GiB free (need ≥{need_gb} GiB)");
    Ok(())
}

fn free_gb(path: &Path) -> Result<u64> {
    // Prefer `df -g` (macOS) / `df -BG` (GNU); fall back to df -k.
    let probe = if path.exists() {
        path.to_path_buf()
    } else {
        path.parent()
            .unwrap_or_else(|| Path::new("."))
            .to_path_buf()
    };

    if let Ok(out) = Command::new("df")
        .args(["-g", &probe.to_string_lossy()])
        .output()
    {
        if out.status.success() {
            if let Some(v) = parse_df_avail_field(&String::from_utf8_lossy(&out.stdout), 1.0) {
                return Ok(v);
            }
        }
    }
    if let Ok(out) = Command::new("df")
        .args(["-k", &probe.to_string_lossy()])
        .output()
    {
        if out.status.success() {
            if let Some(v) =
                parse_df_avail_field(&String::from_utf8_lossy(&out.stdout), 1024.0 * 1024.0)
            {
                return Ok(v);
            }
        }
    }
    // Non-fatal if df parsing fails — warn and continue.
    eprintln!("[gatk-rs-equiv] warning: could not determine free disk; continuing");
    Ok(need_passthrough_ok())
}

fn need_passthrough_ok() -> u64 {
    u64::MAX / 4
}

fn parse_df_avail_field(stdout: &str, unit_to_gib: f64) -> Option<u64> {
    let line = stdout.lines().nth(1)?;
    let cols: Vec<&str> = line.split_whitespace().collect();
    // macOS df -g: Filesystem Size Used Avail Capacity...
    // avail is typically column index 3 (0-based).
    let avail_s = cols.get(3)?;
    let avail: f64 = avail_s.trim_end_matches('G').parse().ok()?;
    Some((avail / unit_to_gib).floor() as u64)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn clamps_threads() {
        assert_eq!(clamp_threads(0), 1);
        assert_eq!(clamp_threads(2), 2);
        assert_eq!(clamp_threads(64), MAX_THREADS);
    }
}
