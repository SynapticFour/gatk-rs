//! Memory management utilities for GATK-RS
//! This module provides memory-efficient data structures and utilities for handling
//! large genomic datasets with zero-copy operations where possible.

use bytes::Bytes;
use dashmap::DashMap;
use lru::LruCache;
use parking_lot::RwLock;
use rayon::prelude::*;
use rayon::slice::Iter as ParSliceIter;
use std::sync::Arc;

/// Memory pool for reusable byte buffer allocations.
/// # Invariants
/// Buffers are keyed by power-of-two capacity classes; per-class pool capped at 10 buffers.
/// Individual pooled buffers larger than 1 MiB are not retained.
/// # Ownership
/// Shared via `&self`; callers own returned `Vec<u8>` until returned with [`MemoryPool::return_buffer`].
/// # Mutation
/// Interior mutability via `DashMap`; safe to share across threads.
/// # Biological assumptions
/// None (infrastructure).
/// # Java equivalence
/// None / Rust-native (general allocation pooling).
pub struct MemoryPool {
    pools: DashMap<usize, Vec<Vec<u8>>>,
    max_pool_size: usize,
}

impl MemoryPool {
    /// Create a new memory pool
    pub fn new(max_memory_mb: usize) -> Self {
        let max_pool_size = max_memory_mb * 1024 * 1024;
        Self {
            pools: DashMap::new(),
            max_pool_size,
        }
    }

    /// Get a buffer from the pool or allocate a new one
    pub fn get_buffer(&self, size: usize) -> Vec<u8> {
        // Find the appropriate pool size (round up to next power of 2)
        let pool_size = size.next_power_of_two();

        if let Some(mut pool) = self.pools.get_mut(&pool_size) {
            if let Some(buffer) = pool.pop() {
                drop(pool);
                let mut buffer = buffer;
                buffer.resize(size, 0);
                return buffer;
            }
        }

        // No buffer available, allocate new one
        vec![0u8; size]
    }

    /// Return a buffer to the pool
    pub fn return_buffer(&self, mut buffer: Vec<u8>) {
        let pool_size = buffer.capacity().next_power_of_two();

        // Don't pool buffers that are too large
        if pool_size > 1024 * 1024 {
            // 1MB limit per buffer
            return;
        }

        let mut pool = self.pools.entry(pool_size).or_default();
        if pool.len() < 10 {
            // Limit pool size per size class
            buffer.clear();
            pool.push(buffer);
        }
    }

    /// Get current memory usage statistics
    pub fn memory_usage(&self) -> MemoryUsage {
        let mut total_pooled = 0;
        let mut pool_count = 0;

        for pool in self.pools.iter() {
            total_pooled += pool.value().len() * pool.key();
            pool_count += pool.value().len();
        }

        MemoryUsage {
            total_pooled_bytes: total_pooled,
            pool_count,
            max_pool_size: self.max_pool_size,
        }
    }
}

/// Snapshot of pooled buffer memory usage.
/// # Invariants
/// `total_pooled_bytes <= max_pool_size` in normal operation (approximate accounting).
/// # Ownership
/// `Copy`-less small struct; clone for reporting.
/// # Mutation
/// Immutable snapshot from [`MemoryPool::memory_usage`].
/// # Biological assumptions
/// None (infrastructure).
/// # Java equivalence
/// None / Rust-native.
#[derive(Debug, Clone)]
pub struct MemoryUsage {
    pub total_pooled_bytes: usize,
    pub pool_count: usize,
    pub max_pool_size: usize,
}

impl MemoryUsage {
    /// Get memory usage in MB
    pub fn pooled_mb(&self) -> f64 {
        self.total_pooled_bytes as f64 / (1024.0 * 1024.0)
    }
}

/// Zero-copy view into shared [`bytes::Bytes`] storage.
/// # Invariants
/// `len` and sub-slices stay within the backing `Bytes` bounds.
/// # Ownership
/// Holds `Arc`-like shared storage via `Bytes`; clone is shallow.
/// # Mutation
/// Immutable view; create new slices via [`ByteSlice::slice`].
/// # Biological assumptions
/// None (infrastructure); used for file/sequence byte views.
/// # Java equivalence
/// None / Rust-native.
#[derive(Debug, Clone)]
pub struct ByteSlice {
    data: Bytes,
    offset: usize,
    len: usize,
}

impl ByteSlice {
    /// Create a new byte slice
    pub fn new(data: Bytes) -> Self {
        let len = data.len();
        Self {
            data,
            offset: 0,
            len,
        }
    }

    /// Create a sub-slice
    pub fn slice(&self, offset: usize, len: usize) -> Option<Self> {
        if offset + len <= self.len {
            Some(Self {
                data: self
                    .data
                    .slice((self.offset + offset)..(self.offset + offset + len)),
                offset: 0,
                len,
            })
        } else {
            None
        }
    }

    /// Get the underlying data
    pub fn as_bytes(&self) -> &[u8] {
        &self.data
    }

    /// Get the length
    pub fn len(&self) -> usize {
        self.len
    }

    /// Check if empty
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }
}

/// Thread-safe LRU cache for genomic key/value pairs.
/// # Invariants
/// Capacity is at least 1; LRU order maintained by inner `LruCache`.
/// # Ownership
/// Keys and values are cloned on `get`/`put`; cache shared via `Arc<RwLock<...>>`.
/// # Mutation
/// `get`/`put` take write lock (LRU promotion); concurrent access safe.
/// # Biological assumptions
/// Generic cache; typical keys are contig/window identifiers.
/// # Java equivalence
/// None / Rust-native.
pub struct GenomicCache<K, V> {
    cache: Arc<RwLock<LruCache<K, V>>>,
    max_size: usize,
}

impl<K: Clone + Eq + std::hash::Hash, V: Clone> GenomicCache<K, V> {
    /// Create a new cache
    pub fn new(max_size: usize) -> Self {
        let effective_size = max_size.max(1);
        Self {
            cache: Arc::new(RwLock::new(LruCache::new(
                std::num::NonZeroUsize::new(effective_size).expect("non-zero by construction"),
            ))),
            max_size: effective_size,
        }
    }

    /// Get a value from cache
    pub fn get(&self, key: &K) -> Option<V> {
        self.cache.write().get(key).cloned()
    }

    /// Put a value in cache
    pub fn put(&self, key: K, value: V) {
        self.cache.write().put(key, value);
    }

    /// Get cache statistics
    pub fn stats(&self) -> CacheStats {
        let cache = self.cache.read();
        CacheStats {
            len: cache.len(),
            max_size: self.max_size,
            is_full: cache.len() >= self.max_size,
        }
    }
}

/// LRU cache utilization snapshot.
/// # Invariants
/// `len <= max_size`; `utilization` is percentage in `[0, 100]`.
/// # Ownership
/// Small owned snapshot; clone freely.
/// # Mutation
/// Read-only stats from [`GenomicCache::stats`].
/// # Biological assumptions
/// None (infrastructure).
/// # Java equivalence
/// None / Rust-native.
#[derive(Debug, Clone)]
pub struct CacheStats {
    pub len: usize,
    pub max_size: usize,
    pub is_full: bool,
}

impl CacheStats {
    /// Get cache utilization as percentage
    pub fn utilization(&self) -> f64 {
        if self.max_size == 0 {
            0.0
        } else {
            (self.len as f64 / self.max_size as f64) * 100.0
        }
    }
}

/// Linear interval collection with parallel overlap queries.
/// # Invariants
/// Intervals are stored unsorted until [`IntervalTree::sort`]; overlap tests use inclusive coordinates.
/// # Ownership
/// Owns `Vec<GenomicInterval<T>>`; returns borrowed `&T` from queries.
/// # Mutation
/// `insert` and `sort` require `&mut self`; queries are `&self`.
/// # Biological assumptions
/// Intervals indexed by chromosome name string and 1-based-style numeric range (caller convention).
/// # Java equivalence
/// Loosely similar to GATK `IntervalTree` / `GenomeLocTree` (simpler, linear scan + sort).
pub struct IntervalTree<T> {
    intervals: Vec<GenomicInterval<T>>,
}

/// Genomic interval with attached payload of type `T`.
/// # Invariants
/// `start <= end` expected; same chromosome string compared literally in queries.
/// # Ownership
/// Owns `chromosome` and `data`; clone to duplicate interval + payload.
/// # Mutation
/// Public fields; typically built then inserted into [`IntervalTree`].
/// # Biological assumptions
/// Closed interval on a named contig; strand not modeled.
/// # Java equivalence
/// Approximates `SimpleInterval` + optional feature payload (Rust-native pairing).
#[derive(Debug, Clone)]
pub struct GenomicInterval<T> {
    pub chromosome: String,
    pub start: u64,
    pub end: u64,
    pub data: T,
}

impl<T: Send + Sync> IntervalTree<T> {
    /// Create a new interval tree
    pub fn new() -> Self {
        Self {
            intervals: Vec::new(),
        }
    }

    /// Add an interval
    pub fn insert(&mut self, interval: GenomicInterval<T>) {
        self.intervals.push(interval);
    }

    /// Find intervals overlapping a position
    pub fn find_overlapping(&self, chromosome: &str, position: u64) -> Vec<&T> {
        self.intervals
            .par_iter()
            .filter(|interval| {
                interval.chromosome == chromosome
                    && interval.start <= position
                    && interval.end >= position
            })
            .map(|interval| &interval.data)
            .collect()
    }

    /// Find intervals overlapping a range
    pub fn find_overlapping_range(&self, chromosome: &str, start: u64, end: u64) -> Vec<&T> {
        self.intervals
            .par_iter()
            .filter(|interval| {
                interval.chromosome == chromosome && interval.start <= end && interval.end >= start
            })
            .map(|interval| &interval.data)
            .collect()
    }

    /// Sort intervals by position for faster queries
    pub fn sort(&mut self) {
        self.intervals.par_sort_by(|a, b| {
            a.chromosome
                .cmp(&b.chromosome)
                .then(a.start.cmp(&b.start))
                .then(a.end.cmp(&b.end))
        });
    }

    /// Get the number of intervals
    pub fn len(&self) -> usize {
        self.intervals.len()
    }

    /// Check if empty
    pub fn is_empty(&self) -> bool {
        self.intervals.is_empty()
    }
}

impl<'a, T: Sync + 'a> rayon::iter::IntoParallelRefIterator<'a> for IntervalTree<T> {
    type Item = &'a GenomicInterval<T>;
    type Iter = ParSliceIter<'a, GenomicInterval<T>>;

    fn par_iter(&'a self) -> Self::Iter {
        self.intervals.par_iter()
    }
}

impl<T> Default for IntervalTree<T> {
    fn default() -> Self {
        Self {
            intervals: Vec::new(),
        }
    }
}

/// Chunked stream reader using a shared [`MemoryPool`].
/// # Invariants
/// Reads up to `chunk_size` bytes per iteration until EOF.
/// # Ownership
/// Borrows [`MemoryPool`] via `Arc`; does not own the input reader beyond each call.
/// # Mutation
/// `process_chunks` mutates caller-supplied reader and callback state.
/// # Biological assumptions
/// None (infrastructure I/O).
/// # Java equivalence
/// None / Rust-native.
pub struct StreamProcessor {
    chunk_size: usize,
    buffer_pool: Arc<MemoryPool>,
}

impl StreamProcessor {
    /// Create a new stream processor
    pub fn new(chunk_size: usize, memory_pool: Arc<MemoryPool>) -> Self {
        Self {
            chunk_size,
            buffer_pool: memory_pool,
        }
    }

    /// Process data in chunks
    pub fn process_chunks<F, R>(&self, reader: R, mut processor: F) -> gatk_common::GatkResult<()>
    where
        F: FnMut(&[u8]) -> gatk_common::GatkResult<()>,
        R: std::io::Read,
    {
        let mut buffer = self.buffer_pool.get_buffer(self.chunk_size);
        let mut reader = reader;

        loop {
            let bytes_read = reader
                .read(&mut buffer)
                .map_err(|e| gatk_common::GatkError::io("Read error", e))?;
            if bytes_read == 0 {
                break;
            }

            processor(&buffer[..bytes_read])?;
        }

        self.buffer_pool.return_buffer(buffer);
        Ok(())
    }
}

/// Memory-mapped read-only file wrapper.
/// # Invariants
/// Mapping covers entire file; slices must stay within `file_size`.
/// # Ownership
/// Owns `memmap2::Mmap`; exposes borrowed `&[u8]` slices.
/// # Mutation
/// Immutable after map; no write mapping.
/// # Biological assumptions
/// None (infrastructure); used for large FASTA/FASTQ inputs.
/// # Java equivalence
/// Similar to memory-mapped reference/read access patterns in htsjdk (Rust-native API).
pub struct MemoryMappedFile {
    mmap: memmap2::Mmap,
    file_size: usize,
}

impl MemoryMappedFile {
    /// Open a file for memory-mapped access
    pub fn open(path: &str) -> gatk_common::GatkResult<Self> {
        let file = std::fs::File::open(path)
            .map_err(|e| gatk_common::GatkError::io("File open error", e))?;
        let mmap = unsafe {
            memmap2::Mmap::map(&file).map_err(|e| gatk_common::GatkError::io("Memory map error", e))
        }?;
        let file_size = mmap.len();

        Ok(Self { mmap, file_size })
    }

    /// Get a slice of the file
    pub fn slice(&self, offset: usize, len: usize) -> Option<&[u8]> {
        if offset + len <= self.file_size {
            Some(&self.mmap[offset..offset + len])
        } else {
            None
        }
    }

    /// Get the entire file as bytes
    pub fn as_bytes(&self) -> &[u8] {
        &self.mmap
    }

    /// Get file size
    pub fn len(&self) -> usize {
        self.file_size
    }

    /// Check if empty
    pub fn is_empty(&self) -> bool {
        self.file_size == 0
    }
}

/// Rayon-based parallel map over cloned items with coarse memory chunking.
/// # Invariants
/// Items must be `Send + Sync + Clone`; processor invoked once per item.
/// # Ownership
/// Takes ownership of input `Vec<T>`; returns owned `Vec<R>`.
/// # Mutation
/// Stateless unit struct; no instance fields.
/// # Biological assumptions
/// None (infrastructure).
/// # Java equivalence
/// None / Rust-native.
pub struct ParallelProcessor;

impl ParallelProcessor {
    /// Process items in parallel with controlled memory usage
    pub fn process_parallel<T, R, F>(items: Vec<T>, processor: F, max_memory_mb: usize) -> Vec<R>
    where
        T: Send + Sync + Clone,
        R: Send,
        F: Fn(T) -> R + Send + Sync,
    {
        // Calculate optimal chunk size based on memory constraints
        let chunk_size = (items.len() / (max_memory_mb * 2)).max(1);

        items
            .par_chunks(chunk_size)
            .flat_map(|chunk| {
                chunk
                    .iter()
                    .map(|item| processor(item.clone()))
                    .collect::<Vec<_>>()
            })
            .collect()
    }
}

/// Process memory usage tracker (placeholder platform hook).
/// # Invariants
/// `initial_memory` captured at construction; delta relative to that baseline.
/// # Ownership
/// Stateless aside from baseline snapshot; cheap to clone not implemented—create per scope.
/// # Mutation
/// `current_memory_usage` is a static placeholder returning 0 until platform integration.
/// # Biological assumptions
/// None (infrastructure).
/// # Java equivalence
/// None / Rust-native.
pub struct MemoryMonitor {
    initial_memory: usize,
}

impl MemoryMonitor {
    /// Create a new memory monitor
    pub fn new() -> Self {
        Self {
            initial_memory: Self::current_memory_usage(),
        }
    }

    /// Get current memory usage in bytes
    pub fn current_memory_usage() -> usize {
        // This is a simplified implementation
        // In production, you'd use platform-specific APIs
        0 // Placeholder
    }

    /// Get memory usage since monitoring started
    pub fn memory_delta(&self) -> isize {
        Self::current_memory_usage() as isize - self.initial_memory as isize
    }

    /// Check if memory usage is within limits
    pub fn within_limit(&self, max_mb: usize) -> bool {
        Self::current_memory_usage() <= max_mb * 1024 * 1024
    }
}

impl Default for MemoryMonitor {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_memory_pool() {
        let pool = MemoryPool::new(100); // 100MB

        let buffer1 = pool.get_buffer(1024);
        let buffer2 = pool.get_buffer(2048);

        assert_eq!(buffer1.len(), 1024);
        assert_eq!(buffer2.len(), 2048);

        pool.return_buffer(buffer1);
        pool.return_buffer(buffer2);

        let stats = pool.memory_usage();
        assert!(stats.pooled_mb() > 0.0);
    }

    #[test]
    fn test_byte_slice() {
        let data = Bytes::from("Hello, World!");
        let slice = ByteSlice::new(data);

        assert_eq!(slice.len(), 13);
        assert!(!slice.is_empty());

        let sub_slice = slice.slice(7, 5).unwrap();
        assert_eq!(sub_slice.as_bytes(), b"World");
    }

    #[test]
    fn test_genomic_cache() {
        let cache = GenomicCache::new(10);

        cache.put("chr1".to_string(), "sequence1".to_string());
        cache.put("chr2".to_string(), "sequence2".to_string());

        assert_eq!(
            cache.get(&"chr1".to_string()),
            Some("sequence1".to_string())
        );
        assert_eq!(
            cache.get(&"chr2".to_string()),
            Some("sequence2".to_string())
        );

        let stats = cache.stats();
        assert_eq!(stats.len, 2);
        assert_eq!(stats.max_size, 10);
    }

    #[test]
    fn test_interval_tree() {
        let mut tree = IntervalTree::new();

        tree.insert(GenomicInterval {
            chromosome: "chr1".to_string(),
            start: 100,
            end: 200,
            data: "interval1",
        });

        tree.insert(GenomicInterval {
            chromosome: "chr1".to_string(),
            start: 150,
            end: 250,
            data: "interval2",
        });

        let results = tree.find_overlapping("chr1", 180);
        assert_eq!(results.len(), 2);
    }
}
