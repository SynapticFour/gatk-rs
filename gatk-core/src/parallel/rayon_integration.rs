use crate::parallel::ParallelConfig;
use rayon::prelude::*;

/// Dedicated Rayon thread pool wrapper for data-parallel genomic batches.
/// # Invariants
/// Pool size from `rayon_pool_size` or `worker_threads` in config (minimum 1).
/// # Ownership
/// Owns `rayon::ThreadPool`; share processor via `&self` across threads.
/// # Mutation
/// Pool immutable after construction; work installed via `process_items_parallel`.
/// # Biological assumptions
/// None (generic parallel executor).
/// # Java equivalence
/// None / Rust-native (replaces Java `ForkJoinPool` usage patterns).
pub struct RayonProcessor {
    pool: rayon::ThreadPool,
}

impl RayonProcessor {
    pub fn new(config: ParallelConfig) -> gatk_common::GatkResult<Self> {
        let threads = config
            .rayon_pool_size
            .unwrap_or(config.worker_threads.max(1));
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .map_err(|e| {
                gatk_common::GatkError::generic(format!("Failed to build rayon pool: {e}"))
            })?;
        Ok(Self { pool })
    }

    pub fn process_items_parallel<T, R, F>(
        &self,
        items: Vec<T>,
        processor: F,
    ) -> gatk_common::GatkResult<Vec<R>>
    where
        T: Send,
        R: Send,
        F: Fn(T) -> R + Send + Sync,
    {
        let out = self
            .pool
            .install(|| items.into_par_iter().map(processor).collect());
        Ok(out)
    }

    pub fn reverse_complement_parallel(
        &self,
        sequences: Vec<Vec<u8>>,
    ) -> gatk_common::GatkResult<Vec<Vec<u8>>> {
        self.process_items_parallel(sequences, |seq| {
            seq.into_iter()
                .rev()
                .map(|b| match b {
                    b'A' => b'T',
                    b'T' => b'A',
                    b'G' => b'C',
                    b'C' => b'G',
                    b'a' => b't',
                    b't' => b'a',
                    b'g' => b'c',
                    b'c' => b'g',
                    _ => b,
                })
                .collect()
        })
    }

    pub fn gc_content_parallel(
        &self,
        sequences: Vec<Vec<u8>>,
    ) -> gatk_common::GatkResult<Vec<f64>> {
        self.process_items_parallel(sequences, |seq| {
            if seq.is_empty() {
                return 0.0;
            }
            let gc = seq
                .iter()
                .filter(|&&b| matches!(b, b'G' | b'C' | b'g' | b'c'))
                .count();
            gc as f64 / seq.len() as f64
        })
    }
}
