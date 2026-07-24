//! FASTA file parser and writer for GATK-RS
//! This module provides efficient parsing and writing of FASTA files
//! with support for large files, memory-mapped access, and
//! streaming operations.

#![warn(clippy::unwrap_used, clippy::expect_used)]

use crate::memory::MemoryMappedFile;
use indexmap::IndexMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom, Write};
use std::path::Path;

/// FASTA sequence representation.
/// # Invariants
/// `length` equals `sequence.len` after construction.
/// # Ownership
/// Owns name, optional description, and sequence bytes.
/// # Mutation
/// Immutable helpers return new sequences (e.g. reverse complement).
/// # Biological assumptions
/// Reference or assembly contig sequence; IUPAC bytes allowed in `sequence`.
/// # Java equivalence
/// htsjdk `ReferenceSequence` / FASTA record.
#[derive(Debug, Clone)]
pub struct FastaSequence {
    pub name: String,
    pub description: Option<String>,
    pub sequence: Vec<u8>,
    pub length: usize,
}

impl FastaSequence {
    pub fn new(name: String, sequence: Vec<u8>) -> Self {
        let length = sequence.len();
        Self {
            name,
            description: None,
            sequence,
            length,
        }
    }

    pub fn with_description(name: String, description: Option<String>, sequence: Vec<u8>) -> Self {
        let length = sequence.len();
        Self {
            name,
            description,
            sequence,
            length,
        }
    }

    pub fn as_bytes(&self) -> &[u8] {
        &self.sequence
    }

    pub fn len(&self) -> usize {
        self.length
    }

    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    pub fn gc_content(&self) -> f64 {
        if self.sequence.is_empty() {
            return 0.0;
        }
        let gc = self
            .sequence
            .iter()
            .filter(|&&b| matches!(b, b'G' | b'C' | b'g' | b'c'))
            .count();
        gc as f64 / self.sequence.len() as f64
    }

    pub fn reverse_complement(&self) -> Self {
        let rev = self
            .sequence
            .iter()
            .rev()
            .map(|&b| match b {
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
            .collect::<Vec<_>>();
        Self {
            name: self.name.clone(),
            description: self.description.clone(),
            length: rev.len(),
            sequence: rev,
        }
    }

    pub fn subsequence(&self, start: usize, len: usize) -> Vec<u8> {
        if start >= self.sequence.len() {
            return Vec::new();
        }
        let end = (start + len).min(self.sequence.len());
        self.sequence[start..end].to_vec()
    }

    pub fn count_pattern(&self, pattern: &[u8]) -> usize {
        if pattern.is_empty() || pattern.len() > self.sequence.len() {
            return 0;
        }
        self.sequence
            .windows(pattern.len())
            .filter(|w| *w == pattern)
            .count()
    }

    pub fn is_valid_dna(&self) -> bool {
        self.sequence.iter().all(|&b| {
            matches!(
                b,
                b'A' | b'C' | b'G' | b'T' | b'N' | b'a' | b'c' | b'g' | b't' | b'n'
            )
        })
    }
}

/// FASTA-level metadata collected during reads.
/// # Invariants
/// `total_length` is sum of parsed sequence lengths; `sequence_count` increments per record.
/// # Ownership
/// Plain scalars; clone with parent metadata snapshots.
/// # Mutation
/// Updated by readers while streaming sequences.
/// # Biological assumptions
/// Summarizes reference FASTA composition.
/// # Java equivalence
/// Similar to htsjdk `ReferenceSequenceFile` statistics (Rust-native struct).
#[derive(Debug, Clone, Default)]
pub struct FastaMetadata {
    pub sequence_count: usize,
    pub total_length: usize,
}

/// Single-sequence entry in a samtools-style FASTA index (`.fai`).
/// # Invariants
/// `offset` points to first base byte in the FASTA file; line width fields match wrapped FASTA layout.
/// `length` is total bases in the sequence (excluding newlines).
/// # Ownership
/// Owns sequence `name`; clone for index maps.
/// # Mutation
/// Immutable index row after load/build.
/// # Biological assumptions
/// One reference contig per entry; coordinates elsewhere are 1-based against `length`.
/// # Java equivalence
/// Same role as htsjdk `FastaSequenceIndex` / samtools `.fai` record.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct FastaIndexEntry {
    pub name: String,
    pub offset: u64,
    pub length: usize,
    pub line_length: usize,
    pub bytes_per_line: usize,
}

/// In-memory map of FASTA sequence names to [`FastaIndexEntry`] offsets.
/// # Invariants
/// `file_path` is the indexed FASTA; `sequences` keys match contig names in the file headers.
/// [`Self::sequences`] is an [`IndexMap`] preserving **FASTA encounter order** so
/// [`Self::sequence_names`] is deterministic (not hash-iteration order).
/// # Ownership
/// Owns path and index map; clone duplicates index in memory.
/// # Mutation
/// Built once; random access uses immutable lookups.
/// # Biological assumptions
/// Reference index for O(1) contig lookup and byte-range seeks.
/// # Java equivalence
/// htsjdk `ReferenceSequenceFile` + `FastaSequenceIndex`.
#[derive(Debug, Clone)]
pub struct FastaIndex {
    pub sequences: IndexMap<String, FastaIndexEntry>,
    pub file_path: std::path::PathBuf,
}

/// FASTA file reader with streaming and memory-mapped support.
/// # Invariants
/// Backend chosen by file size (>100 MiB → mmap); exposes sequential `read_next_sequence`.
/// # Ownership
/// Owns boxed reader backend; borrows nothing from caller paths.
/// # Mutation
/// `read_next_sequence` mutates reader cursor.
/// # Biological assumptions
/// Standard FASTA with wrapped sequence lines.
/// # Java equivalence
/// htsjdk `ReferenceSequenceFile`.
pub struct FastaReader {
    reader: Box<dyn FastaReaderBackend>,
}

/// FASTA file writer.
/// # Invariants
/// Wraps buffered file handle; caller supplies valid FASTA records.
/// # Ownership
/// Owns writer and underlying file.
/// # Mutation
/// Write methods append to file.
/// # Biological assumptions
/// Emits reference sequences for downstream indexing.
/// # Java equivalence
/// htsjdk FASTA writer utilities (conceptual).
pub struct FastaWriter {
    writer: std::io::BufWriter<File>,
}

/// Backend trait for FASTA reading implementations
trait FastaReaderBackend {
    fn read_next_sequence(&mut self) -> gatk_common::GatkResult<Option<FastaSequence>>;
    fn get_metadata(&self) -> &FastaMetadata;
}

/// Standard buffered FASTA reader
struct BufferedFastaReader {
    reader: BufReader<File>,
    metadata: FastaMetadata,
    line_buffer: String,
    pending_header: Option<String>,
}

/// Memory-mapped FASTA reader
struct MemoryMappedFastaReader {
    mmap_file: MemoryMappedFile,
    metadata: FastaMetadata,
    position: usize,
}

impl FastaReader {
    /// Create a new FASTA reader from file path
    pub fn from_file<P: AsRef<Path>>(path: P) -> gatk_common::GatkResult<Self> {
        let file = File::open(&path)
            .map_err(|e| gatk_common::GatkError::io("Failed to open FASTA file", e))?;

        // Determine if we should use memory mapping
        let file_size = file
            .metadata()
            .map_err(|e| gatk_common::GatkError::io("Failed to get FASTA file metadata", e))?
            .len();

        let reader: Box<dyn FastaReaderBackend> = if file_size > 100 * 1024 * 1024 {
            // > 100MB
            // Use memory mapping for large files
            let mmap_file = MemoryMappedFile::open(path.as_ref().to_string_lossy().as_ref())?;
            Box::new(MemoryMappedFastaReader::new(mmap_file)?)
        } else {
            // Use buffered reading for smaller files
            Box::new(BufferedFastaReader::new(file)?)
        };

        Ok(Self { reader })
    }

    /// Create a new FASTA reader with explicit buffering
    pub fn from_file_buffered<P: AsRef<Path>>(path: P) -> gatk_common::GatkResult<Self> {
        let file = File::open(&path)
            .map_err(|e| gatk_common::GatkError::io("Failed to open FASTA file", e))?;

        let reader = Box::new(BufferedFastaReader::new(file)?);

        Ok(Self { reader })
    }

    /// Create a new FASTA reader with memory mapping
    pub fn from_file_memory_mapped<P: AsRef<Path>>(path: P) -> gatk_common::GatkResult<Self> {
        let mmap_file = MemoryMappedFile::open(path.as_ref().to_string_lossy().as_ref())?;
        let reader = Box::new(MemoryMappedFastaReader::new(mmap_file)?);

        Ok(Self { reader })
    }

    /// Read next sequence from FASTA file
    pub fn read_next_sequence(&mut self) -> gatk_common::GatkResult<Option<FastaSequence>> {
        self.reader.read_next_sequence()
    }

    /// Read all sequences from FASTA file
    pub fn read_all_sequences(&mut self) -> gatk_common::GatkResult<Vec<FastaSequence>> {
        let mut sequences = Vec::new();

        while let Some(sequence) = self.read_next_sequence()? {
            sequences.push(sequence);
        }

        Ok(sequences)
    }

    /// Get sequence by name (requires index)
    pub fn get_sequence_by_name(
        &mut self,
        name: &str,
    ) -> gatk_common::GatkResult<Option<FastaSequence>> {
        // For now, iterate through sequences
        // In practice, would use an index for efficient lookup
        while let Some(sequence) = self.read_next_sequence()? {
            if sequence.name == name {
                return Ok(Some(sequence));
            }
        }
        Ok(None)
    }

    /// Get FASTA metadata
    pub fn metadata(&self) -> &FastaMetadata {
        self.reader.get_metadata()
    }

    /// Create an iterator over sequences
    pub fn iter(&mut self) -> FastaIterator<'_> {
        FastaIterator { reader: self }
    }

    /// Get sequence by position (requires indexed FASTA)
    pub fn get_sequence_by_position(
        &mut self,
        position: usize,
    ) -> gatk_common::GatkResult<Option<FastaSequence>> {
        // For now, iterate to position
        // In practice, would use index for efficient lookup
        for i in 0..=position {
            if let Some(sequence) = self.read_next_sequence()? {
                if i == position {
                    return Ok(Some(sequence));
                }
            } else {
                return Ok(None);
            }
        }
        Ok(None)
    }

    /// Count sequences in file
    pub fn count_sequences(&mut self) -> gatk_common::GatkResult<usize> {
        let mut count = 0;
        while self.read_next_sequence()?.is_some() {
            count += 1;
        }
        Ok(count)
    }

    /// Get total sequence length
    pub fn total_sequence_length(&mut self) -> gatk_common::GatkResult<usize> {
        let mut total_length = 0;
        while let Some(sequence) = self.read_next_sequence()? {
            total_length += sequence.sequence.len();
        }
        Ok(total_length)
    }
}

impl FastaWriter {
    pub fn new<P: AsRef<Path>>(path: P) -> gatk_common::GatkResult<Self> {
        let file = File::create(path)
            .map_err(|e| gatk_common::GatkError::io("Failed to create FASTA file", e))?;
        Ok(Self {
            writer: std::io::BufWriter::new(file),
        })
    }

    pub fn write_sequence(&mut self, sequence: &FastaSequence) -> gatk_common::GatkResult<()> {
        self.writer.write_all(b">")?;
        self.writer.write_all(sequence.name.as_bytes())?;
        if let Some(desc) = &sequence.description {
            self.writer.write_all(b"\t")?;
            self.writer.write_all(desc.as_bytes())?;
        }
        self.writer.write_all(b"\n")?;
        self.writer.write_all(&sequence.sequence)?;
        self.writer.write_all(b"\n")?;
        Ok(())
    }

    pub fn write_sequences(&mut self, sequences: &[FastaSequence]) -> gatk_common::GatkResult<()> {
        for sequence in sequences {
            self.write_sequence(sequence)?;
        }
        Ok(())
    }

    pub fn finish(&mut self) -> gatk_common::GatkResult<()> {
        self.writer
            .flush()
            .map_err(|e| gatk_common::GatkError::io("Failed to flush FASTA writer", e))
    }
}

/// Iterator over FASTA sequences borrowing [`FastaReader`].
/// # Invariants
/// Each `next` yields one complete FASTA record.
/// # Ownership
/// Borrows reader mutably.
/// # Mutation
/// Advances parse cursor.
/// # Biological assumptions
/// None (I/O adapter).
/// # Java equivalence
/// htsjdk reference sequence iteration (conceptual).
pub struct FastaIterator<'a> {
    reader: &'a mut FastaReader,
}

impl<'a> Iterator for FastaIterator<'a> {
    type Item = gatk_common::GatkResult<FastaSequence>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.read_next_sequence() {
            Ok(Some(sequence)) => Some(Ok(sequence)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

impl BufferedFastaReader {
    /// Create new buffered FASTA reader
    fn new(file: File) -> gatk_common::GatkResult<Self> {
        let reader = BufReader::new(file);

        Ok(Self {
            reader,
            metadata: FastaMetadata::default(),
            line_buffer: String::new(),
            pending_header: None,
        })
    }
}

fn parse_fasta_header_name_and_description(header_line: &str) -> (String, Option<String>) {
    let trimmed = header_line.trim();
    if trimmed.is_empty() {
        return (String::new(), None);
    }
    let mut parts = trimmed.splitn(2, char::is_whitespace);
    let name = parts.next().unwrap_or_default().to_string();
    let description = parts
        .next()
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(ToString::to_string);
    (name, description)
}

impl FastaReaderBackend for BufferedFastaReader {
    fn read_next_sequence(&mut self) -> gatk_common::GatkResult<Option<FastaSequence>> {
        let mut header_line = self.pending_header.take().unwrap_or_default();
        let mut sequence_lines = Vec::new();
        let mut found_header = !header_line.is_empty();

        // Read until we find a header line
        while !found_header {
            self.line_buffer.clear();
            let bytes_read = self
                .reader
                .read_line(&mut self.line_buffer)
                .map_err(|e| gatk_common::GatkError::io("Failed to read FASTA line", e))?;

            if bytes_read == 0 {
                return Ok(None); // EOF
            }

            let line = self.line_buffer.trim();
            if line.is_empty() {
                continue;
            }

            if let Some(stripped) = line.strip_prefix('>') {
                header_line = stripped.to_string();
                found_header = true;
                self.metadata.sequence_count += 1;
            }
        }

        // Read sequence lines until next header or EOF
        loop {
            self.line_buffer.clear();
            let bytes_read = self
                .reader
                .read_line(&mut self.line_buffer)
                .map_err(|e| gatk_common::GatkError::io("Failed to read FASTA line", e))?;

            if bytes_read == 0 {
                break; // EOF
            }

            let line = self.line_buffer.trim();
            if line.is_empty() {
                continue;
            }

            if let Some(stripped) = line.strip_prefix('>') {
                // Preserve next header for the next call.
                self.pending_header = Some(stripped.to_string());
                break;
            }

            sequence_lines.push(line.to_string());
        }

        if sequence_lines.is_empty() {
            return Ok(None);
        }

        // Combine sequence lines and remove whitespace
        let sequence: String = sequence_lines.join("");
        let sequence_bytes = sequence.into_bytes();

        let (name, description) = parse_fasta_header_name_and_description(&header_line);

        let fasta_sequence = FastaSequence {
            name,
            description,
            sequence: sequence_bytes.clone(),
            length: sequence_bytes.len(),
        };

        self.metadata.total_length += fasta_sequence.length;

        Ok(Some(fasta_sequence))
    }

    fn get_metadata(&self) -> &FastaMetadata {
        &self.metadata
    }
}

impl MemoryMappedFastaReader {
    /// Create new memory-mapped FASTA reader
    fn new(mmap_file: MemoryMappedFile) -> gatk_common::GatkResult<Self> {
        let metadata = FastaMetadata::default();

        Ok(Self {
            mmap_file,
            metadata,
            position: 0,
        })
    }

    /// Find next header position
    fn find_next_header(&self, start_pos: usize) -> Option<usize> {
        let data = self.mmap_file.as_bytes();
        (start_pos..data.len()).find(|&i| data[i] == b'>')
    }

    /// Find end of current sequence (next header or EOF)
    fn find_sequence_end(&self, start_pos: usize) -> usize {
        let data = self.mmap_file.as_bytes();
        (start_pos..data.len())
            .find(|&i| data[i] == b'>')
            .unwrap_or(data.len())
    }

    /// Extract header and description from header line
    fn extract_header_info(&self, start_pos: usize) -> (String, Option<String>, usize) {
        let data = self.mmap_file.as_bytes();
        let mut end_pos = start_pos;

        // Find end of header line
        while end_pos < data.len() && data[end_pos] != b'\n' && data[end_pos] != b'\r' {
            end_pos += 1;
        }

        if start_pos >= end_pos {
            return (String::new(), None, end_pos);
        }

        let header_bytes = &data[start_pos..end_pos];
        let header_str = String::from_utf8_lossy(header_bytes);
        let (name, description) = parse_fasta_header_name_and_description(&header_str);
        while end_pos < data.len() && (data[end_pos] == b'\n' || data[end_pos] == b'\r') {
            end_pos += 1;
        }
        (name, description, end_pos)
    }

    /// Extract sequence from data range
    fn extract_sequence(&self, start_pos: usize, end_pos: usize) -> Vec<u8> {
        let data = self.mmap_file.as_bytes();
        let sequence_bytes = &data[start_pos..end_pos];

        // Remove whitespace and newlines
        sequence_bytes
            .iter()
            .filter(|&&b| b != b'\n' && b != b'\r' && b != b' ' && b != b'\t')
            .copied()
            .collect()
    }
}

impl FastaReaderBackend for MemoryMappedFastaReader {
    fn read_next_sequence(&mut self) -> gatk_common::GatkResult<Option<FastaSequence>> {
        // Find next header
        let header_pos = match self.find_next_header(self.position) {
            Some(pos) => pos,
            None => return Ok(None), // No more headers
        };

        // Extract header info
        let (name, description, sequence_start) = self.extract_header_info(header_pos + 1);

        // Find end of sequence
        let sequence_end = self.find_sequence_end(sequence_start);

        // Extract sequence
        let sequence = self.extract_sequence(sequence_start, sequence_end);

        // Update position
        self.position = sequence_end;

        let fasta_sequence = FastaSequence {
            name: name.clone(),
            description,
            length: sequence.len(),
            sequence,
        };

        self.metadata.sequence_count += 1;
        self.metadata.total_length += fasta_sequence.length;

        Ok(Some(fasta_sequence))
    }

    fn get_metadata(&self) -> &FastaMetadata {
        &self.metadata
    }
}

impl FastaIndex {
    /// Create new FASTA index
    pub fn new() -> Self {
        Self {
            sequences: IndexMap::new(),
            file_path: std::path::PathBuf::new(),
        }
    }

    /// Build index from FASTA file
    pub fn build_from_file<P: AsRef<Path>>(path: P) -> gatk_common::GatkResult<Self> {
        let file_path = path.as_ref().to_path_buf();
        let mut reader = FastaReader::from_file(&file_path)?;
        let mut index = Self {
            sequences: IndexMap::new(),
            file_path,
        };

        let mut current_offset = 0u64;

        while let Some(sequence) = reader.read_next_sequence()? {
            let entry = FastaIndexEntry {
                name: sequence.name.clone(),
                offset: current_offset,
                length: sequence.length,
                line_length: 80,    // Default, would need to detect
                bytes_per_line: 81, // 80 + newline
            };

            // CLONE: needed because owned HashMap/BTree/HashSet key or value.
            index.sequences.insert(sequence.name.clone(), entry);

            current_offset += (sequence.name.len() + 1) as u64; // +1 for '>'
            current_offset += sequence.length.div_ceil(80) as u64; // Lines + newlines
        }

        Ok(index)
    }

    /// Get sequence by name using index
    pub fn get_sequence(
        &self,
        name: &str,
        start: usize,
        length: usize,
    ) -> gatk_common::GatkResult<Option<Vec<u8>>> {
        if let Some(entry) = self.sequences.get(name) {
            let file = File::open(&self.file_path).map_err(|e| {
                gatk_common::GatkError::io("Failed to open FASTA file for indexed access", e)
            })?;

            let mut reader = std::io::BufReader::new(file);

            // Seek to sequence offset
            reader
                .seek(SeekFrom::Start(entry.offset))
                .map_err(|e| gatk_common::GatkError::io("Failed to seek in FASTA file", e))?;

            // Read sequence portion
            let mut buffer = vec![0u8; length];
            let bytes_read = reader
                .read(&mut buffer)
                .map_err(|e| gatk_common::GatkError::io("Failed to read from FASTA file", e))?;

            if bytes_read < length {
                buffer.truncate(bytes_read);
            }

            // Remove newlines and adjust for line wrapping
            let cleaned_sequence =
                self.remove_newlines(&buffer, start, entry.line_length, entry.bytes_per_line);

            Ok(Some(cleaned_sequence))
        } else {
            Ok(None)
        }
    }

    /// Remove newlines from sequence data
    fn remove_newlines(
        &self,
        data: &[u8],
        start: usize,
        line_length: usize,
        _bytes_per_line: usize,
    ) -> Vec<u8> {
        let mut result = Vec::new();
        let mut data_pos = 0;
        let mut line_pos = start % line_length;

        while data_pos < data.len() && result.len() < data.len() - start {
            if line_pos < line_length {
                result.push(data[data_pos]);
                line_pos += 1;
            } else {
                // Skip newline
                line_pos = 0;
            }
            data_pos += 1;
        }

        result
    }

    /// Save index to file
    pub fn save_to_file<P: AsRef<Path>>(&self, path: P) -> gatk_common::GatkResult<()> {
        let serialized = serde_json::to_string(&self.sequences)
            .map_err(|_e| gatk_common::GatkError::generic("Failed to serialize FASTA index"))?;

        std::fs::write(path, serialized)
            .map_err(|e| gatk_common::GatkError::io("Failed to write FASTA index", e))?;

        Ok(())
    }

    /// Load index from file
    pub fn load_from_file<P: AsRef<Path>>(path: P) -> gatk_common::GatkResult<Self> {
        let index_path = path.as_ref().to_path_buf();
        let content = std::fs::read_to_string(&index_path)
            .map_err(|e| gatk_common::GatkError::io("Failed to read FASTA index", e))?;

        let sequences: IndexMap<String, FastaIndexEntry> = serde_json::from_str(&content)
            .map_err(|_e| gatk_common::GatkError::generic("Failed to deserialize FASTA index"))?;

        Ok(Self {
            sequences,
            file_path: index_path,
        })
    }

    /// Get all sequence names in FASTA / index insertion order (deterministic).
    pub fn sequence_names(&self) -> Vec<String> {
        self.sequences.keys().cloned().collect()
    }

    /// Get sequence entry
    pub fn get_entry(&self, name: &str) -> Option<&FastaIndexEntry> {
        self.sequences.get(name)
    }
}

impl Default for FastaIndex {
    fn default() -> Self {
        Self::new()
    }
}
