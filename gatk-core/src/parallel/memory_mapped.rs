use crate::parallel::ParallelConfig;
use std::collections::HashMap;

/// Helper for memory-mapped genomic file access and line indexing.
/// # Invariants
/// Holds config for future mmap tuning; current reads may load full file into RAM.
/// # Ownership
/// Owns [`ParallelConfig`]; returns owned [`ReadOnlyMmapFile`] handles.
/// # Mutation
/// `create_index` writes sidecar index files.
/// # Biological assumptions
/// FASTA/text genomic formats with newline-delimited records.
/// # Java equivalence
/// None / Rust-native (htsjdk uses different indexing models).
pub struct MemoryMappedProcessor {
    _config: ParallelConfig,
}

impl MemoryMappedProcessor {
    pub fn new(config: ParallelConfig) -> Self {
        Self { _config: config }
    }

    pub fn open_read_only<P: AsRef<std::path::Path>>(
        &self,
        path: P,
    ) -> gatk_common::GatkResult<ReadOnlyMmapFile> {
        ReadOnlyMmapFile::open(path)
    }

    pub fn process_fasta_mapped<P: AsRef<std::path::Path>>(
        &self,
        path: P,
    ) -> gatk_common::GatkResult<ReadOnlyMmapFile> {
        ReadOnlyMmapFile::open(path)
    }

    pub fn create_index<P: AsRef<std::path::Path>>(
        &self,
        input: P,
        output: P,
    ) -> gatk_common::GatkResult<SimpleLineIndex> {
        let content = std::fs::read(input.as_ref())
            .map_err(|e| gatk_common::GatkError::io("Failed to read file for index", e))?;
        let mut offsets = vec![0usize];
        for (idx, b) in content.iter().enumerate() {
            if *b == b'\n' {
                offsets.push(idx + 1);
            }
        }
        let body = offsets
            .iter()
            .map(|v| v.to_string())
            .collect::<Vec<_>>()
            .join("\n");
        std::fs::write(output.as_ref(), body)
            .map_err(|e| gatk_common::GatkError::io("Failed to write index", e))?;
        Ok(SimpleLineIndex { offsets })
    }
}

/// Byte offsets of each line start in a text genomic file.
/// # Invariants
/// `offsets[0] == 0`; monotonically increasing byte positions.
/// # Ownership
/// Owns offset vector; cheap clone for sharing index in memory.
/// # Mutation
/// Immutable after construction.
/// # Biological assumptions
/// None (text indexing infrastructure).
/// # Java equivalence
/// None / Rust-native.
pub struct SimpleLineIndex {
    offsets: Vec<usize>,
}

impl SimpleLineIndex {
    pub fn line_count(&self) -> usize {
        self.offsets.len()
    }
}

/// Read-only in-memory file view (full read today; mmap-ready API).
/// # Invariants
/// `data` holds entire file contents after `open`.
/// # Ownership
/// Owns byte buffer; exposes borrowed slices via `data`.
/// # Mutation
/// Immutable after load.
/// # Biological assumptions
/// Used for FASTA/FASTQ/text parsing stubs.
/// # Java equivalence
/// None / Rust-native.
pub struct ReadOnlyMmapFile {
    data: Vec<u8>,
}

impl ReadOnlyMmapFile {
    pub fn open<P: AsRef<std::path::Path>>(path: P) -> gatk_common::GatkResult<Self> {
        let data = std::fs::read(path.as_ref())
            .map_err(|e| gatk_common::GatkError::io("Failed to open mapped file", e))?;
        Ok(Self { data })
    }

    pub fn data(&self) -> &[u8] {
        &self.data
    }

    pub fn get_sequence(&self, name: &str) -> Option<Vec<u8>> {
        let text = std::str::from_utf8(&self.data).ok()?;
        let mut current_name: Option<&str> = None;
        let mut seqs: HashMap<String, Vec<u8>> = HashMap::new();
        for line in text.lines() {
            if let Some(rest) = line.strip_prefix('>') {
                current_name = Some(rest.trim());
                seqs.entry(rest.trim().to_string()).or_default();
            } else if let Some(n) = current_name {
                seqs.entry(n.to_string())
                    .or_default()
                    .extend_from_slice(line.trim().as_bytes());
            }
        }
        seqs.remove(name)
    }
}
