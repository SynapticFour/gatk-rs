//! Rayon global pool / `RAYON_NUM_THREADS` initialization for HaplotypeCaller.

/// Install the Rayon global pool with `n` threads.
/// When `override_env` is true (CLI `--threads` was set), also writes
/// `RAYON_NUM_THREADS` so nested tools see the same size.
pub fn init_worker_threads(n: usize, override_env: bool) {
    let n = n.max(1);
    if override_env {
        // SAFETY: called once at CLI/process startup before parallel work.
        unsafe {
            std::env::set_var("RAYON_NUM_THREADS", n.to_string());
        }
    }
    if let Err(e) = rayon::ThreadPoolBuilder::new()
        .num_threads(n)
        .build_global()
    {
        // Pool may already exist (tests / re-entry); keep installed pool.
        tracing::debug!("rayon global pool already initialized: {e}");
    }
}
