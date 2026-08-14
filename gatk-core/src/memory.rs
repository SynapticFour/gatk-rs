//! Memory-mapped read-only file access for FASTA/FASTQ backends.

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
    /// Open a file for memory-mapped access.
    pub fn open(path: &str) -> gatk_common::GatkResult<Self> {
        let file = std::fs::File::open(path)
            .map_err(|e| gatk_common::GatkError::io("File open error", e))?;
        // SAFETY: `Mmap` keeps the mapping alive after `file` is dropped; the OS
        // holds the inode. We never write through this mapping (`Mmap` is read-only).
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
