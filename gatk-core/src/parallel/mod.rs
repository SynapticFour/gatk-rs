//! Distributed Computing and Parallel Processing for GATK-RS
//! This module provides high-performance parallel processing capabilities
//! using Rayon for data parallelism, Tokio for async I/O, and Polars for
//! efficient data processing of genomic datasets.

pub mod async_io;
pub mod distributed;
pub mod memory_mapped;
pub mod parallel_algorithms;
pub mod polars_integration;
pub mod rayon_integration;
pub mod tokio_integration;

use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::sync::Arc;
use std::time::Duration;

/// Parallel processing configuration (Rayon, Tokio, chunking, memory limits).
/// # Invariants
/// `worker_threads == 0` means auto-detect at use sites that honor it.
/// `chunk_size` and `memory_limit_mb` guide schedulers, not hard OS caps.
/// # Ownership
/// Owns optional pool sizes; clone for worker contexts.
/// # Mutation
/// Typically immutable per pipeline run; public fields for tuning.
/// # Biological assumptions
/// None (infrastructure).
/// # Java equivalence
/// Loosely similar to GATK `-nt` / Spark partition settings (Rust-native struct).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ParallelConfig {
    /// Number of worker threads (0 = auto-detect)
    pub worker_threads: usize,
    /// Thread pool size for Rayon
    pub rayon_pool_size: Option<usize>,
    /// Tokio runtime configuration
    pub tokio_threads: Option<usize>,
    /// Chunk size for parallel processing
    pub chunk_size: usize,
    /// Memory limit per worker (MB)
    pub memory_limit_mb: usize,
    /// Enable memory-mapped file processing
    pub enable_memory_mapping: bool,
    /// Enable async I/O operations
    pub enable_async_io: bool,
    /// Enable distributed processing
    pub enable_distributed: bool,
}

impl Default for ParallelConfig {
    fn default() -> Self {
        let num_cpus = num_cpus::get();
        Self {
            worker_threads: num_cpus,
            rayon_pool_size: Some(num_cpus),
            tokio_threads: Some(num_cpus),
            chunk_size: 1000,
            memory_limit_mb: 1024, // 1GB per worker
            enable_memory_mapping: true,
            enable_async_io: true,
            enable_distributed: false,
        }
    }
}

/// Timing and efficiency metrics from a parallel execution.
/// # Invariants
/// `parallel_efficiency` updated via [`ParallelStats::calculate_efficiency`] when workers > 0.
/// # Ownership
/// Plain scalars and `Duration`; clone for reporting.
/// # Mutation
/// Fields may be filled incrementally by profilers.
/// # Biological assumptions
/// None (infrastructure).
/// # Java equivalence
/// None / Rust-native.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ParallelStats {
    /// Total processing time
    pub total_time: Duration,
    /// Parallel processing time
    pub parallel_time: Duration,
    /// Sequential processing time
    pub sequential_time: Duration,
    /// Number of worker threads used
    pub workers_used: usize,
    /// Memory usage per worker (MB)
    pub memory_per_worker: f64,
    /// Parallel efficiency (speedup / workers)
    pub parallel_efficiency: f64,
    /// Load balancing score (0-1, higher is better)
    pub load_balance_score: f64,
    /// Cache hit rate
    pub cache_hit_rate: f64,
}

impl ParallelStats {
    /// Calculate parallel efficiency
    pub fn calculate_efficiency(&mut self) {
        if self.workers_used > 0 {
            let theoretical_speedup = self.workers_used as f64;
            let actual_speedup =
                self.sequential_time.as_secs_f64() / self.parallel_time.as_secs_f64();
            self.parallel_efficiency = actual_speedup / theoretical_speedup;
        }
    }
}

/// Shared Rayon thread pool and Tokio runtime for parallel genomic workloads.
/// # Invariants
/// Constructed once per config; thread pools outlive individual tasks installed on them.
/// # Ownership
/// Owns `Arc` handles to Rayon pool and Tokio runtime; cheap to clone not implemented—share via `Arc<ParallelContext>` externally.
/// # Mutation
/// Runtimes are thread-safe; config is immutable after construction.
/// # Biological assumptions
/// None (infrastructure).
/// # Java equivalence
/// None / Rust-native (replaces ad-hoc Java thread pools in tooling).
pub struct ParallelContext {
    config: ParallelConfig,
    rayon_pool: Arc<rayon::ThreadPool>,
    tokio_runtime: Arc<tokio::runtime::Runtime>,
}

impl ParallelContext {
    /// Create a new parallel processing context
    pub fn new(config: ParallelConfig) -> gatk_common::GatkResult<Self> {
        // Create Rayon thread pool
        let rayon_pool = if let Some(pool_size) = config.rayon_pool_size {
            Arc::new(
                rayon::ThreadPoolBuilder::new()
                    .num_threads(pool_size)
                    .thread_name(|index| format!("gatk-rayon-{}", index))
                    .build()
                    .map_err(|_e| {
                        gatk_common::GatkError::generic("Failed to create Rayon thread pool")
                    })?,
            )
        } else {
            Arc::new(
                rayon::ThreadPoolBuilder::new()
                    .thread_name(|index| format!("gatk-rayon-{}", index))
                    .build()
                    .map_err(|_e| {
                        gatk_common::GatkError::generic("Failed to create Rayon thread pool")
                    })?,
            )
        };

        // Create Tokio runtime
        let _tokio_threads = config.tokio_threads.unwrap_or_else(num_cpus::get);
        let tokio_runtime = Arc::new(
            tokio::runtime::Runtime::new()
                .map_err(|_e| gatk_common::GatkError::generic("Failed to create Tokio runtime"))?,
        );

        Ok(Self {
            config,
            rayon_pool,
            tokio_runtime,
        })
    }

    /// Get the parallel configuration
    pub fn config(&self) -> &ParallelConfig {
        &self.config
    }

    /// Get the Rayon thread pool
    pub fn rayon_pool(&self) -> &rayon::ThreadPool {
        &self.rayon_pool
    }

    /// Get the Tokio runtime
    pub fn tokio_runtime(&self) -> &tokio::runtime::Runtime {
        &self.tokio_runtime
    }

    /// Execute a parallel operation using Rayon
    pub fn execute_parallel<F, R>(&self, operation: F) -> gatk_common::GatkResult<R>
    where
        F: FnOnce() -> R + Send,
        R: Send,
    {
        let result = self.rayon_pool.install(operation);
        Ok(result)
    }

    /// Execute an async operation using Tokio
    pub async fn execute_async<F, R>(&self, operation: F) -> gatk_common::GatkResult<R>
    where
        F: std::future::Future<Output = gatk_common::GatkResult<R>> + Send,
        R: Send,
    {
        operation.await
    }

    /// Block on an async operation
    pub fn block_on<F, R>(&self, operation: F) -> gatk_common::GatkResult<R>
    where
        F: std::future::Future<Output = gatk_common::GatkResult<R>> + Send,
        R: Send,
    {
        self.tokio_runtime.block_on(operation)
    }
}

/// Splits owned data into chunks and processes them with [`ParallelConfig`].
/// # Invariants
/// Chunk count equals input vector length (one item per chunk entry).
/// # Ownership
/// Owns input `Vec<T>`; output `Vec<R>` owned by caller.
/// # Mutation
/// Parallel path clones chunk items; sequential path clones on iterate.
/// # Biological assumptions
/// None (generic batching).
/// # Java equivalence
/// None / Rust-native.
pub struct ChunkedProcessor<T, R> {
    chunks: Vec<T>,
    config: ParallelConfig,
    phantom: std::marker::PhantomData<R>,
}

impl<T, R: Sync> ChunkedProcessor<T, R> {
    /// Create a new chunked processor
    pub fn new(data: Vec<T>, config: ParallelConfig) -> Self {
        Self {
            chunks: data,
            config,
            phantom: std::marker::PhantomData,
        }
    }

    /// Process chunks in parallel
    pub fn process_chunks_parallel<F>(&self, processor: F) -> gatk_common::GatkResult<Vec<R>>
    where
        F: Fn(T) -> R + Send + Sync,
        T: Send + Clone + Sync,
        R: Send,
    {
        let context = ParallelContext::new(self.config.clone())?;
        context.execute_parallel(|| {
            self.chunks
                .par_iter()
                .map(|chunk| processor(chunk.clone()))
                .collect()
        })
    }

    /// Process chunks sequentially
    pub fn process_chunks_sequential<F>(&self, processor: F) -> Vec<R>
    where
        F: Fn(T) -> R,
        T: Clone,
    {
        self.chunks
            .iter()
            .map(|chunk| processor(chunk.clone()))
            .collect()
    }
}

/// Async task gate that blocks scheduling when estimated memory would exceed limit.
/// # Invariants
/// Tracks reserved MB atomically; waits with 10ms sleeps when over budget.
/// # Ownership
/// Shared via `&self`; tasks run on Tokio blocking pool.
/// # Mutation
/// Atomic counter updated around each scheduled task.
/// # Biological assumptions
/// None (infrastructure).
/// # Java equivalence
/// None / Rust-native.
pub struct MemoryAwareScheduler {
    memory_limit_mb: usize,
    current_memory_usage: std::sync::atomic::AtomicUsize,
}

impl MemoryAwareScheduler {
    /// Create a new memory-aware scheduler
    pub fn new(memory_limit_mb: usize) -> Self {
        Self {
            memory_limit_mb,
            current_memory_usage: std::sync::atomic::AtomicUsize::new(0),
        }
    }

    /// Schedule a task with memory requirements
    pub async fn schedule_task<F, R>(&self, memory_mb: usize, task: F) -> gatk_common::GatkResult<R>
    where
        F: FnOnce() -> R + Send + 'static,
        R: Send + 'static,
    {
        // Wait if memory limit would be exceeded
        while self
            .current_memory_usage
            .load(std::sync::atomic::Ordering::Relaxed)
            + memory_mb
            > self.memory_limit_mb
        {
            tokio::time::sleep(Duration::from_millis(10)).await;
        }

        // Reserve memory
        self.current_memory_usage
            .fetch_add(memory_mb, std::sync::atomic::Ordering::Relaxed);

        // Execute task
        let result = tokio::task::spawn_blocking(task)
            .await
            .map_err(|_e| gatk_common::GatkError::generic("Task execution failed"))?;

        // Release memory
        self.current_memory_usage
            .fetch_sub(memory_mb, std::sync::atomic::Ordering::Relaxed);

        Ok(result)
    }

    /// Get current memory usage
    pub fn current_memory_usage(&self) -> usize {
        self.current_memory_usage
            .load(std::sync::atomic::Ordering::Relaxed)
    }

    /// Get memory limit
    pub fn memory_limit(&self) -> usize {
        self.memory_limit_mb
    }
}

/// Worker load tracker with least-loaded and round-robin selection.
/// # Invariants
/// Worker ids are `0..num_workers`; load counters saturate at usize bounds.
/// # Ownership
/// Shared via `&self` with atomic loads.
/// # Mutation
/// `increment_load` / `decrement_load` must be paired by callers.
/// # Biological assumptions
/// None (infrastructure).
/// # Java equivalence
/// None / Rust-native.
pub struct LoadBalancer {
    worker_loads: Vec<std::sync::atomic::AtomicUsize>,
    round_robin_index: std::sync::atomic::AtomicUsize,
}

impl LoadBalancer {
    /// Create a new load balancer
    pub fn new(num_workers: usize) -> Self {
        Self {
            worker_loads: (0..num_workers)
                .map(|_| std::sync::atomic::AtomicUsize::new(0))
                .collect(),
            round_robin_index: std::sync::atomic::AtomicUsize::new(0),
        }
    }

    /// Get the least loaded worker
    pub fn get_least_loaded_worker(&self) -> usize {
        let mut min_load = usize::MAX;
        let mut min_worker = 0;

        for (i, load) in self.worker_loads.iter().enumerate() {
            let current_load = load.load(std::sync::atomic::Ordering::Relaxed);
            if current_load < min_load {
                min_load = current_load;
                min_worker = i;
            }
        }

        min_worker
    }

    /// Get worker using round-robin
    pub fn get_round_robin_worker(&self) -> usize {
        let index = self
            .round_robin_index
            .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        index % self.worker_loads.len()
    }

    /// Increment worker load
    pub fn increment_load(&self, worker_id: usize) {
        if worker_id < self.worker_loads.len() {
            self.worker_loads[worker_id].fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        }
    }

    /// Decrement worker load
    pub fn decrement_load(&self, worker_id: usize) {
        if worker_id < self.worker_loads.len() {
            self.worker_loads[worker_id].fetch_sub(1, std::sync::atomic::Ordering::Relaxed);
        }
    }

    /// Get load balance score (0-1, higher is better)
    pub fn load_balance_score(&self) -> f64 {
        if self.worker_loads.is_empty() {
            return 1.0;
        }

        let loads: Vec<usize> = self
            .worker_loads
            .iter()
            .map(|load| load.load(std::sync::atomic::Ordering::Relaxed))
            .collect();

        let total_load: usize = loads.iter().sum();
        if total_load == 0 {
            return 1.0;
        }

        let expected_load = total_load as f64 / loads.len() as f64;
        let variance: f64 = loads
            .iter()
            .map(|&load| {
                let diff = load as f64 - expected_load;
                diff * diff
            })
            .sum();

        let std_dev = variance.sqrt();
        let coefficient_of_variation = std_dev / expected_load;

        // Convert to 0-1 scale (lower CV = higher score)
        (1.0 / (1.0 + coefficient_of_variation)).min(1.0)
    }
}

/// Collects per-worker timings and memory for [`ParallelStats`] generation.
/// # Invariants
/// Worker ids must be `< worker_metrics.len` when recording.
/// # Ownership
/// Owns per-worker metric vectors; not `Sync`—use from one thread or wrap externally.
/// # Mutation
/// `record_worker_activity` mutates internal aggregates.
/// # Biological assumptions
/// None (infrastructure).
/// # Java equivalence
/// None / Rust-native.
pub struct ParallelProfiler {
    start_time: std::time::Instant,
    worker_metrics: Vec<WorkerMetrics>,
}

impl ParallelProfiler {
    /// Create a new profiler
    pub fn new(num_workers: usize) -> Self {
        Self {
            start_time: std::time::Instant::now(),
            worker_metrics: (0..num_workers).map(|_| WorkerMetrics::new()).collect(),
        }
    }

    /// Record worker activity
    pub fn record_worker_activity(
        &mut self,
        worker_id: usize,
        duration: Duration,
        memory_mb: usize,
    ) {
        if worker_id < self.worker_metrics.len() {
            self.worker_metrics[worker_id].record_activity(duration, memory_mb);
        }
    }

    /// Generate performance statistics
    pub fn generate_stats(&self) -> ParallelStats {
        let total_time = self.start_time.elapsed();

        let total_parallel_time: Duration = self
            .worker_metrics
            .iter()
            .map(|m| m.total_processing_time)
            .sum();

        let avg_memory_per_worker: f64 = self
            .worker_metrics
            .iter()
            .map(|m| m.peak_memory_mb as f64)
            .sum::<f64>()
            / self.worker_metrics.len() as f64;

        ParallelStats {
            total_time,
            parallel_time: total_parallel_time,
            sequential_time: Duration::from_secs(0), // Would need separate measurement
            workers_used: self.worker_metrics.len(),
            memory_per_worker: avg_memory_per_worker,
            parallel_efficiency: 0.0, // Calculated separately
            load_balance_score: 0.0,  // Would need load balancer reference
            cache_hit_rate: 0.0,      // Would need cache metrics
        }
    }
}

/// Worker performance metrics
#[derive(Debug, Clone)]
struct WorkerMetrics {
    total_processing_time: Duration,
    peak_memory_mb: usize,
    tasks_completed: usize,
}

impl WorkerMetrics {
    fn new() -> Self {
        Self {
            total_processing_time: Duration::from_secs(0),
            peak_memory_mb: 0,
            tasks_completed: 0,
        }
    }

    fn record_activity(&mut self, duration: Duration, memory_mb: usize) {
        self.total_processing_time += duration;
        self.peak_memory_mb = self.peak_memory_mb.max(memory_mb);
        self.tasks_completed += 1;
    }
}

/// Utility functions for parallel processing
pub mod utils {
    use super::*;

    /// Split data into chunks for parallel processing
    pub fn chunk_data<T: Clone>(data: &[T], chunk_size: usize) -> Vec<Vec<T>> {
        data.chunks(chunk_size)
            .map(|chunk| chunk.to_vec())
            .collect()
    }

    /// Calculate optimal chunk size based on data size and worker count
    pub fn calculate_optimal_chunk_size(data_size: usize, worker_count: usize) -> usize {
        let base_chunk_size = data_size / worker_count;
        // Ensure chunks aren't too small (overhead) or too large (load imbalance)
        base_chunk_size.clamp(100, 10000)
    }

    /// Estimate memory usage for data processing
    pub fn estimate_memory_usage<T>(data: &[T]) -> usize {
        std::mem::size_of_val(data) / (1024 * 1024) // Convert to MB
    }

    /// Check if parallel processing is beneficial
    pub fn should_use_parallel(data_size: usize, sequential_time: Duration) -> bool {
        // Use parallel if data is large enough to justify overhead
        let min_data_size = 1000;
        let min_sequential_time = Duration::from_millis(100);

        data_size >= min_data_size && sequential_time >= min_sequential_time
    }
}
