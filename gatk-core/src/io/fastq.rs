//! FASTQ file parser and writer for GATK-RS
//! This module provides efficient parsing and writing of FASTQ files
//! with support for large files, memory-mapped access, streaming
//! operations, and quality score handling.

use crate::memory::MemoryMappedFile;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

/// FASTQ file reader with streaming and memory-mapped support.
/// # Invariants
/// Four-line FASTQ records; quality encoding detected or assumed Phred+33.
/// # Ownership
/// Owns reader backend and metadata.
/// # Mutation
/// Sequential reads mutate cursor.
/// # Biological assumptions
/// Unaligned sequencing reads with per-base qualities.
/// # Java equivalence
/// htsjdk `FastqReader` (conceptual).
pub struct FastqReader {
    reader: Box<dyn FastqReaderBackend>,
}

/// Backend trait for FASTQ reading implementations
trait FastqReaderBackend {
    fn read_next_read(&mut self) -> gatk_common::GatkResult<Option<FastqRead>>;
    fn get_metadata(&self) -> &FastqMetadata;
}

/// Standard buffered FASTQ reader
struct BufferedFastqReader {
    reader: BufReader<File>,
    metadata: FastqMetadata,
    line_buffer: String,
}

/// Memory-mapped FASTQ reader
struct MemoryMappedFastqReader {
    mmap_file: MemoryMappedFile,
    metadata: FastqMetadata,
    position: usize,
}

impl FastqReader {
    /// Create a new FASTQ reader from file path
    pub fn from_file<P: AsRef<Path>>(path: P) -> gatk_common::GatkResult<Self> {
        let file = File::open(&path)
            .map_err(|e| gatk_common::GatkError::io("Failed to open FASTQ file", e))?;

        // Determine if we should use memory mapping
        let file_size = file
            .metadata()
            .map_err(|e| gatk_common::GatkError::io("Failed to get FASTQ file metadata", e))?
            .len();

        let reader: Box<dyn FastqReaderBackend> = if file_size > 100 * 1024 * 1024 {
            // > 100MB
            // Use memory mapping for large files
            let mmap_file = MemoryMappedFile::open(path.as_ref().to_string_lossy().as_ref())?;
            Box::new(MemoryMappedFastqReader::new(mmap_file)?)
        } else {
            // Use buffered reading for smaller files
            Box::new(BufferedFastqReader::new(file)?)
        };

        Ok(Self { reader })
    }

    /// Create a new FASTQ reader with explicit buffering
    pub fn from_file_buffered<P: AsRef<Path>>(path: P) -> gatk_common::GatkResult<Self> {
        let file = File::open(&path)
            .map_err(|e| gatk_common::GatkError::io("Failed to open FASTQ file", e))?;

        let reader = Box::new(BufferedFastqReader::new(file)?);

        Ok(Self { reader })
    }

    /// Create a new FASTQ reader with memory mapping
    pub fn from_file_memory_mapped<P: AsRef<Path>>(path: P) -> gatk_common::GatkResult<Self> {
        let mmap_file = MemoryMappedFile::open(path.as_ref().to_string_lossy().as_ref())?;
        let reader = Box::new(MemoryMappedFastqReader::new(mmap_file)?);

        Ok(Self { reader })
    }

    /// Read next read from FASTQ file
    pub fn read_next_read(&mut self) -> gatk_common::GatkResult<Option<FastqRead>> {
        self.reader.read_next_read()
    }

    /// Read all reads from FASTQ file
    pub fn read_all_reads(&mut self) -> gatk_common::GatkResult<Vec<FastqRead>> {
        let mut reads = Vec::new();

        while let Some(read) = self.read_next_read()? {
            reads.push(read);
        }

        Ok(reads)
    }

    /// Get FASTQ metadata
    pub fn metadata(&self) -> &FastqMetadata {
        self.reader.get_metadata()
    }

    /// Create an iterator over reads
    pub fn iter(&mut self) -> FastqIterator<'_> {
        FastqIterator { reader: self }
    }

    /// Count reads in file
    pub fn count_reads(&mut self) -> gatk_common::GatkResult<usize> {
        let mut count = 0;
        while self.read_next_read()?.is_some() {
            count += 1;
        }
        Ok(count)
    }

    /// Get total bases in file
    pub fn total_bases(&mut self) -> gatk_common::GatkResult<usize> {
        let mut total_bases = 0;
        while let Some(read) = self.read_next_read()? {
            total_bases += read.sequence.len();
        }
        Ok(total_bases)
    }

    /// Filter reads by quality
    pub fn filter_by_quality(
        &mut self,
        min_quality: u8,
    ) -> gatk_common::GatkResult<Vec<FastqRead>> {
        let mut filtered_reads = Vec::new();

        while let Some(read) = self.read_next_read()? {
            if read.average_quality() >= f64::from(min_quality) {
                filtered_reads.push(read);
            }
        }

        Ok(filtered_reads)
    }

    /// Filter reads by length
    pub fn filter_by_length(
        &mut self,
        min_length: usize,
        max_length: usize,
    ) -> gatk_common::GatkResult<Vec<FastqRead>> {
        let mut filtered_reads = Vec::new();

        while let Some(read) = self.read_next_read()? {
            let read_len = read.sequence.len();
            if read_len >= min_length && read_len <= max_length {
                filtered_reads.push(read);
            }
        }

        Ok(filtered_reads)
    }

    /// Sample reads from file
    pub fn sample_reads(&mut self, sample_size: usize) -> gatk_common::GatkResult<Vec<FastqRead>> {
        let mut sampled_reads = Vec::new();
        let mut total_reads = 0;

        // First, count total reads
        let original_position = self.get_position();
        while self.read_next_read()?.is_some() {
            total_reads += 1;
        }
        self.seek_to_position(original_position)?;

        // Calculate sampling interval
        let interval = if total_reads > 0 {
            total_reads / sample_size
        } else {
            1
        };

        // Sample reads
        let mut count = 0;
        while let Some(read) = self.read_next_read()? {
            if count % interval == 0 && sampled_reads.len() < sample_size {
                sampled_reads.push(read);
            }
            count += 1;
        }

        Ok(sampled_reads)
    }

    /// Get current position (simplified)
    fn get_position(&self) -> usize {
        0 // Would need to track actual position
    }

    /// Seek to position (simplified)
    fn seek_to_position(&mut self, _position: usize) -> gatk_common::GatkResult<()> {
        // Would need to implement proper seeking
        Ok(())
    }
}

/// Iterator over FASTQ reads borrowing [`FastqReader`].
/// # Invariants
/// Each item is one four-line FASTQ record.
/// # Ownership
/// Borrows reader mutably.
/// # Mutation
/// Advances reader on `next`.
/// # Biological assumptions
/// None (I/O adapter).
/// # Java equivalence
/// htsjdk FASTQ iterators (conceptual).
pub struct FastqIterator<'a> {
    reader: &'a mut FastqReader,
}

impl<'a> Iterator for FastqIterator<'a> {
    type Item = gatk_common::GatkResult<FastqRead>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.read_next_read() {
            Ok(Some(read)) => Some(Ok(read)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

impl BufferedFastqReader {
    /// Create new buffered FASTQ reader
    fn new(file: File) -> gatk_common::GatkResult<Self> {
        let reader = BufReader::new(file);

        Ok(Self {
            reader,
            metadata: FastqMetadata::default(),
            line_buffer: String::new(),
        })
    }
}

impl FastqReaderBackend for BufferedFastqReader {
    fn read_next_read(&mut self) -> gatk_common::GatkResult<Option<FastqRead>> {
        // Read header line
        self.line_buffer.clear();
        let bytes_read = self
            .reader
            .read_line(&mut self.line_buffer)
            .map_err(|e| gatk_common::GatkError::io("Failed to read FASTQ header line", e))?;

        if bytes_read == 0 {
            return Ok(None); // EOF
        }

        if !self.line_buffer.starts_with('@') {
            return Err(gatk_common::GatkError::generic(
                "Invalid FASTQ format: expected header line starting with '@'",
            ));
        }

        let header = self.line_buffer[1..].trim_end().to_string(); // Remove '@' and trim

        // Read sequence line
        self.line_buffer.clear();
        self.reader
            .read_line(&mut self.line_buffer)
            .map_err(|e| gatk_common::GatkError::io("Failed to read FASTQ sequence line", e))?;
        let sequence = self.line_buffer.trim_end().as_bytes().to_vec();

        // Read plus line
        self.line_buffer.clear();
        self.reader
            .read_line(&mut self.line_buffer)
            .map_err(|e| gatk_common::GatkError::io("Failed to read FASTQ plus line", e))?;

        if !self.line_buffer.starts_with('+') {
            return Err(gatk_common::GatkError::generic(
                "Invalid FASTQ format: expected plus line starting with '+'",
            ));
        }

        // Read quality line
        self.line_buffer.clear();
        self.reader
            .read_line(&mut self.line_buffer)
            .map_err(|e| gatk_common::GatkError::io("Failed to read FASTQ quality line", e))?;
        let quality = self.line_buffer.trim_end().as_bytes().to_vec();

        // Validate sequence and quality lengths match
        if sequence.len() != quality.len() {
            return Err(gatk_common::GatkError::generic(format!(
                "FASTQ sequence length ({}) does not match quality length ({})",
                sequence.len(),
                quality.len()
            )));
        }

        let fastq_read = FastqRead {
            header,
            sequence,
            quality,
        };

        self.metadata.read_count += 1;
        self.metadata.total_bases += fastq_read.sequence.len();

        Ok(Some(fastq_read))
    }

    fn get_metadata(&self) -> &FastqMetadata {
        &self.metadata
    }
}

impl MemoryMappedFastqReader {
    /// Create new memory-mapped FASTQ reader
    fn new(mmap_file: MemoryMappedFile) -> gatk_common::GatkResult<Self> {
        let metadata = FastqMetadata::default();

        Ok(Self {
            mmap_file,
            metadata,
            position: 0,
        })
    }

    /// Find next header position
    fn find_next_header(&self, start_pos: usize) -> Option<usize> {
        let data = self.mmap_file.as_bytes();
        for (i, byte) in data.iter().enumerate().skip(start_pos) {
            if *byte == b'@' {
                // Check if it's a valid header (not in quality line)
                if i == 0 || data[i - 1] == b'\n' || data[i - 1] == b'\r' {
                    return Some(i);
                }
            }
        }
        None
    }

    /// Find end of current line
    fn find_line_end(&self, start_pos: usize) -> usize {
        let data = self.mmap_file.as_bytes();
        for (i, byte) in data.iter().enumerate().skip(start_pos) {
            if *byte == b'\n' || *byte == b'\r' {
                return i;
            }
        }
        data.len()
    }

    /// Extract line as string
    fn extract_line(&self, start_pos: usize, end_pos: usize) -> String {
        let data = self.mmap_file.as_bytes();
        String::from_utf8_lossy(&data[start_pos..end_pos]).to_string()
    }

    fn extract_line_bytes(&self, start_pos: usize, end_pos: usize) -> Vec<u8> {
        let data = self.mmap_file.as_bytes();
        data[start_pos..end_pos].to_vec()
    }
}

impl FastqReaderBackend for MemoryMappedFastqReader {
    fn read_next_read(&mut self) -> gatk_common::GatkResult<Option<FastqRead>> {
        // Find next header
        let header_pos = match self.find_next_header(self.position) {
            Some(pos) => pos,
            None => return Ok(None), // No more headers
        };

        // Extract header
        let header_end = self.find_line_end(header_pos);
        let header = self.extract_line(header_pos + 1, header_end);

        // Extract sequence
        let sequence_start = header_end + 1;
        let sequence_end = self.find_line_end(sequence_start);
        let sequence = self.extract_line_bytes(sequence_start, sequence_end);

        // Extract plus line (and validate)
        let plus_start = sequence_end + 1;
        let plus_end = self.find_line_end(plus_start);
        let plus_line = self.extract_line(plus_start, plus_end);

        if !plus_line.starts_with('+') {
            return Err(gatk_common::GatkError::generic(
                "Invalid FASTQ format: expected plus line",
            ));
        }

        // Extract quality
        let quality_start = plus_end + 1;
        let quality_end = self.find_line_end(quality_start);
        let quality = self.extract_line_bytes(quality_start, quality_end);

        // Validate sequence and quality lengths match
        if sequence.len() != quality.len() {
            return Err(gatk_common::GatkError::generic(format!(
                "FASTQ sequence length ({}) does not match quality length ({})",
                sequence.len(),
                quality.len()
            )));
        }

        let fastq_read = FastqRead {
            header,
            sequence,
            quality,
        };

        // Update position
        self.position = quality_end + 1;
        self.metadata.read_count += 1;
        self.metadata.total_bases += fastq_read.sequence.len();

        Ok(Some(fastq_read))
    }

    fn get_metadata(&self) -> &FastqMetadata {
        &self.metadata
    }
}

/// FASTQ read representation (header, sequence, qualities).
/// # Invariants
/// `sequence.len` should equal `quality.len` for valid FASTQ.
/// # Ownership
/// Owns header string and byte vectors.
/// # Mutation
/// Immutable helpers compute derived metrics (GC, averages).
/// # Biological assumptions
/// One sequencing read with Phred+33 qualities by default.
/// # Java equivalence
/// htsjdk `FastqRecord`.
#[derive(Debug, Clone)]
pub struct FastqRead {
    /// Header line (without '@')
    pub header: String,
    /// Sequence data
    pub sequence: Vec<u8>,
    /// Quality scores (Phred+33)
    pub quality: Vec<u8>,
}

impl FastqRead {
    /// Create a new FASTQ read
    pub fn new(header: String, sequence: Vec<u8>, quality: Vec<u8>) -> Self {
        Self {
            header,
            sequence,
            quality,
        }
    }

    /// Get read name (first part of header)
    pub fn name(&self) -> &str {
        self.header
            .split_whitespace()
            .next()
            .unwrap_or(&self.header)
    }

    /// Get read description (rest of header after name)
    pub fn description(&self) -> Option<&str> {
        self.header.split_whitespace().nth(1)
    }

    /// Get read length
    pub fn length(&self) -> usize {
        self.sequence.len()
    }

    /// Get GC content
    pub fn gc_content(&self) -> f64 {
        if self.sequence.is_empty() {
            return 0.0;
        }

        let gc_count = self
            .sequence
            .iter()
            .filter(|&&base| matches!(base, b'G' | b'C' | b'g' | b'c'))
            .count();

        gc_count as f64 / self.sequence.len() as f64
    }

    /// Get average quality score
    pub fn average_quality(&self) -> f64 {
        if self.quality.is_empty() {
            return 0.0;
        }

        let sum: f64 = self.quality.iter().map(|&q| q as f64).sum();
        sum / self.quality.len() as f64
    }

    /// Get minimum quality score
    pub fn min_quality(&self) -> u8 {
        self.quality.iter().min().copied().unwrap_or(0)
    }

    /// Get maximum quality score
    pub fn max_quality(&self) -> u8 {
        self.quality.iter().max().copied().unwrap_or(0)
    }

    /// Get median quality score
    pub fn median_quality(&self) -> f64 {
        if self.quality.is_empty() {
            return 0.0;
        }

        let mut sorted_quality = self.quality.clone();
        sorted_quality.sort_unstable();

        let len = sorted_quality.len();
        if len % 2 == 0 {
            (sorted_quality[len / 2 - 1] as f64 + sorted_quality[len / 2] as f64) / 2.0
        } else {
            sorted_quality[len / 2] as f64
        }
    }

    /// Get quality score at position
    pub fn quality_at(&self, position: usize) -> Option<u8> {
        self.quality.get(position).copied()
    }

    /// Get base at position
    pub fn base_at(&self, position: usize) -> Option<u8> {
        self.sequence.get(position).copied()
    }

    /// Trim low-quality bases from both ends
    pub fn trim_quality(&self, min_quality: u8) -> FastqRead {
        let start = self
            .quality
            .iter()
            .position(|&q| q >= min_quality)
            .unwrap_or(0);

        let end = self
            .quality
            .iter()
            .rposition(|&q| q >= min_quality)
            .map(|pos| pos + 1)
            .unwrap_or(self.quality.len());
        // Clamp to both streams so mismatched synthetic/bench inputs cannot panic.
        let end = end.min(self.sequence.len()).min(self.quality.len());
        let start = start.min(end);

        FastqRead {
            header: self.header.clone(),
            sequence: self.sequence[start..end].to_vec(),
            quality: self.quality[start..end].to_vec(),
        }
    }

    /// Trim bases from start and end
    pub fn trim_bases(&self, trim_start: usize, trim_end: usize) -> FastqRead {
        let start = trim_start;
        let end = self.sequence.len().saturating_sub(trim_end);

        FastqRead {
            header: self.header.clone(),
            sequence: self.sequence[start..end].to_vec(),
            quality: self.quality[start..end].to_vec(),
        }
    }

    /// Reverse complement the sequence
    pub fn reverse_complement(&self) -> FastqRead {
        let mut rev_comp_sequence = self.sequence.clone();
        rev_comp_sequence.reverse();

        for base in rev_comp_sequence.iter_mut() {
            *base = match base {
                b'A' => b'T',
                b'T' => b'A',
                b'G' => b'C',
                b'C' => b'G',
                b'a' => b't',
                b't' => b'a',
                b'g' => b'c',
                b'c' => b'g',
                _ => *base,
            };
        }

        let mut rev_quality = self.quality.clone();
        rev_quality.reverse();

        FastqRead {
            header: self.header.clone(),
            sequence: rev_comp_sequence,
            quality: rev_quality,
        }
    }

    /// Check if read is valid
    pub fn is_valid(&self) -> bool {
        // Check sequence and quality lengths match
        if self.sequence.len() != self.quality.len() {
            return false;
        }

        // Check for valid DNA bases
        for &base in &self.sequence {
            if !matches!(
                base,
                b'A' | b'T' | b'G' | b'C' | b'N' | b'a' | b't' | b'g' | b'c' | b'n'
            ) {
                return false;
            }
        }

        // Check for valid quality scores
        for &quality in &self.quality {
            if !(33..=126).contains(&quality) {
                // Phred+33 range
                return false;
            }
        }

        true
    }

    /// Get quality histogram
    pub fn quality_histogram(&self) -> Vec<(u8, usize)> {
        let mut histogram = std::collections::HashMap::new();

        for &quality in &self.quality {
            *histogram.entry(quality).or_insert(0) += 1;
        }

        let mut sorted_histogram: Vec<_> = histogram.into_iter().collect();
        sorted_histogram.sort_by_key(|&(quality, _)| quality);
        sorted_histogram
    }

    /// Convert quality scores to Phred+33
    pub fn quality_to_phred33(&self) -> Vec<u8> {
        self.quality.iter().map(|&q| q.saturating_sub(33)).collect()
    }

    /// Convert quality scores from Phred+33
    pub fn quality_from_phred33(&self) -> Vec<u8> {
        self.quality.iter().map(|&q| q.saturating_add(33)).collect()
    }
}

/// FASTQ file metadata accumulated while streaming.
/// # Invariants
/// `average_length` derived from total bases / read count when updated by reader.
/// # Ownership
/// Owns [`QualityFormat`]; clone with reader state.
/// # Mutation
/// Updated incrementally during FASTQ parse.
/// # Biological assumptions
/// Run-level QC summary for FASTQ inputs.
/// # Java equivalence
/// None / Rust-native.
#[derive(Debug, Clone, Default)]
pub struct FastqMetadata {
    /// Number of reads in file
    pub read_count: usize,
    /// Total number of bases
    pub total_bases: usize,
    /// Average read length
    pub average_length: f64,
    /// Quality encoding format
    pub quality_format: QualityFormat,
}

/// Quality score encoding scheme for FASTQ files.
/// # Invariants
/// Default is Phred+33 per modern Illumina conventions.
/// # Ownership
/// `Copy`/`Clone` enum.
/// # Mutation
/// N/A.
/// # Biological assumptions
/// Maps ASCII quality characters to Phred scores.
/// # Java equivalence
/// htsjdk quality encoding utilities (conceptual).
#[derive(Debug, Clone, Default)]
pub enum QualityFormat {
    #[default]
    Phred33,
    Phred64,
    Solexa,
}

impl FastqMetadata {
    /// Calculate average read length
    pub fn calculate_average_length(&mut self) {
        if self.read_count > 0 {
            self.average_length = self.total_bases as f64 / self.read_count as f64;
        }
    }

    /// Get file size estimate
    pub fn estimated_file_size(&self) -> usize {
        // Rough estimate: 4 lines per read + overhead
        self.total_bases + (self.read_count * 100) // 100 bytes per read overhead
    }
}

/// FASTQ writer with efficient output.
/// # Invariants
/// Emits four-line records per write call.
/// # Ownership
/// Owns buffered writer handle.
/// # Mutation
/// Append-only file writes.
/// # Biological assumptions
/// Writes unaligned reads for aligner input.
/// # Java equivalence
/// htsjdk FASTQ writers (conceptual).
pub struct FastqWriter {
    writer: Box<dyn Write>,
}

impl FastqWriter {
    /// Create new FASTQ writer
    pub fn new<P: AsRef<Path>>(path: P) -> gatk_common::GatkResult<Self> {
        let file = File::create(&path)
            .map_err(|e| gatk_common::GatkError::io("Failed to create FASTQ file", e))?;

        Ok(Self {
            writer: Box::new(file),
        })
    }

    /// Write a single read
    pub fn write_read(&mut self, read: &FastqRead) -> gatk_common::GatkResult<()> {
        // Write header
        writeln!(self.writer, "@{}", read.header)
            .map_err(|e| gatk_common::GatkError::io("Failed to write FASTQ header", e))?;

        // Write sequence
        let sequence_str = String::from_utf8_lossy(&read.sequence);
        writeln!(self.writer, "{}", sequence_str)
            .map_err(|e| gatk_common::GatkError::io("Failed to write FASTQ sequence", e))?;

        // Write plus line
        writeln!(self.writer, "+")
            .map_err(|e| gatk_common::GatkError::io("Failed to write FASTQ plus line", e))?;

        // Write quality
        let quality_str = String::from_utf8_lossy(&read.quality);
        writeln!(self.writer, "{}", quality_str)
            .map_err(|e| gatk_common::GatkError::io("Failed to write FASTQ quality", e))?;

        Ok(())
    }

    /// Write multiple reads
    pub fn write_reads(&mut self, reads: &[FastqRead]) -> gatk_common::GatkResult<()> {
        for read in reads {
            self.write_read(read)?;
        }
        Ok(())
    }

    /// Flush writer
    pub fn flush(&mut self) -> gatk_common::GatkResult<()> {
        self.writer
            .flush()
            .map_err(|e| gatk_common::GatkError::io("Failed to flush FASTQ writer", e))
    }

    /// Finish writing and close file
    pub fn finish(mut self) -> gatk_common::GatkResult<()> {
        self.flush()?;
        Ok(())
    }
}

/// Aggregated quality metrics computed from FASTQ reads.
/// # Invariants
/// Means/min/max derived from scanned reads in collector.
/// # Ownership
/// Plain numeric snapshot; clone for reports.
/// # Mutation
/// Immutable stats product.
/// # Biological assumptions
/// QC summary for sequencing run (quality distribution).
/// # Java equivalence
/// None / Rust-native.
pub struct FastqQualityStats {
    pub count: usize,
    pub total_bases: usize,
    pub average_length: f64,
    pub min_length: usize,
    pub max_length: usize,
    pub average_quality: f64,
    pub min_quality: u8,
    pub max_quality: u8,
    pub median_quality: f64,
    pub quality_histogram: Vec<(u8, usize)>,
}

impl FastqQualityStats {
    /// Calculate statistics from reads
    pub fn from_reads(reads: &[FastqRead]) -> Self {
        if reads.is_empty() {
            return Self {
                count: 0,
                total_bases: 0,
                average_length: 0.0,
                min_length: 0,
                max_length: 0,
                average_quality: 0.0,
                min_quality: 0,
                max_quality: 0,
                median_quality: 0.0,
                quality_histogram: Vec::new(),
            };
        }

        let count = reads.len();
        let total_bases: usize = reads.iter().map(|r| r.length()).sum();
        let lengths: Vec<usize> = reads.iter().map(|r| r.length()).collect();
        let qualities: Vec<u8> = reads.iter().flat_map(|r| r.quality.clone()).collect();
        let mut sorted_qualities = qualities.clone();
        let mut quality_histogram = std::collections::HashMap::new();
        for &quality in &qualities {
            *quality_histogram.entry(quality).or_insert(0) += 1;
        }

        let mut sorted_histogram: Vec<_> = quality_histogram.into_iter().collect();
        sorted_histogram.sort_by_key(|&(quality, _)| quality);

        sorted_qualities.sort_unstable();

        Self {
            count,
            total_bases,
            average_length: total_bases as f64 / count as f64,
            min_length: *lengths.iter().min().unwrap_or(&0),
            max_length: *lengths.iter().max().unwrap_or(&0),
            average_quality: qualities.iter().map(|&q| q as f64).sum::<f64>()
                / qualities.len() as f64,
            min_quality: *sorted_qualities.iter().min().unwrap_or(&0),
            max_quality: *sorted_qualities.iter().max().unwrap_or(&0),
            median_quality: if sorted_qualities.len() % 2 == 0 {
                (sorted_qualities[sorted_qualities.len() / 2 - 1] as f64
                    + sorted_qualities[sorted_qualities.len() / 2] as f64)
                    / 2.0
            } else {
                sorted_qualities[sorted_qualities.len() / 2] as f64
            },
            quality_histogram: sorted_histogram,
        }
    }

    /// Get quality distribution
    pub fn quality_distribution(&self) -> Vec<(u8, f64)> {
        let total: usize = self.quality_histogram.iter().map(|(_, count)| *count).sum();
        self.quality_histogram
            .iter()
            .map(|(quality, count)| (*quality, *count as f64 / total as f64))
            .collect()
    }
}
