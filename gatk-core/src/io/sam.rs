//! SAM file parser and writer for GATK-RS
//! This module provides efficient parsing and writing of SAM files
//! with support for large files, streaming operations, and
//! full SAM format specification compliance.

use super::*;
use std::io::{BufRead, Write};

/// SAM file reader with streaming support.
/// # Invariants
/// Text SAM with header `@` lines then alignment rows.
/// # Ownership
/// Owns backend and [`SamHeader`].
/// # Mutation
/// Record reads advance cursor.
/// # Biological assumptions
/// Aligned reads identical semantics to BAM but text encoded.
/// # Java equivalence
/// htsjdk `SamReader`.
pub struct SamReader {
    reader: Box<dyn SamReaderBackend>,
    header: SamHeader,
}

/// SAM file writer with efficient output.
/// # Invariants
/// Header precedes alignment lines in output.
/// # Ownership
/// Owns writer and header.
/// # Mutation
/// Append-only.
/// # Biological assumptions
/// Text SAM interchange format.
/// # Java equivalence
/// htsjdk `SAMFileWriter` (SAM mode).
pub struct SamWriter {
    writer: Box<dyn Write>,
    header: SamHeader,
}

/// SAM header containing reference sequences and read groups.
/// # Invariants
/// Includes optional HD sort/grouping directives when parsed.
/// # Ownership
/// Owns nested BAM-shared header structs.
/// # Mutation
/// Vectors mutable while assembling header.
/// # Biological assumptions
/// SAM metadata block for alignment files.
/// # Java equivalence
/// htsjdk `SAMFileHeader`.
#[derive(Debug, Clone, Default)]
pub struct SamHeader {
    pub reference_sequences: Vec<ReferenceSequence>,
    pub read_groups: Vec<ReadGroup>,
    pub programs: Vec<Program>,
    pub comments: Vec<String>,
    pub sort_order: Option<String>,
    pub grouping: Option<String>,
}

/// SAM alignment record (text form).
/// # Invariants
/// CIGAR stored as string; optional fields parsed to [`OptionalField`].
/// # Ownership
/// Owns all column strings; clone per record.
/// # Mutation
/// Public fields for normalization pipelines.
/// # Biological assumptions
/// One aligned read line equivalent to [`BamRecord`].
/// # Java equivalence
/// htsjdk `SAMRecord`.
#[derive(Debug, Clone)]
pub struct SamRecord {
    pub qname: String,
    pub flag: u16,
    pub rname: String,
    pub pos: i64,
    pub mapq: u8,
    pub cigar: String,
    pub rnext: String,
    pub pnext: i64,
    pub tlen: i64,
    pub seq: String,
    pub qual: String,
    pub optional: Vec<OptionalField>,
}

impl SamReader {
    /// Create new SAM reader from file
    pub fn from_file<P: AsRef<std::path::Path>>(path: P) -> gatk_common::GatkResult<Self> {
        let file = std::fs::File::open(path)
            .map_err(|e| gatk_common::GatkError::io("Failed to open SAM file", e))?;

        let mut reader = Box::new(BufferedSamReader::new(file)?);
        let header = reader.read_header()?;

        Ok(Self { reader, header })
    }

    /// Read next alignment record
    pub fn read_next_record(&mut self) -> gatk_common::GatkResult<Option<SamRecord>> {
        self.reader.read_next_record()
    }

    /// Read all records
    pub fn read_all_records(&mut self) -> gatk_common::GatkResult<Vec<SamRecord>> {
        let mut records = Vec::new();
        while let Some(record) = self.read_next_record()? {
            records.push(record);
        }
        Ok(records)
    }

    /// Get header
    pub fn header(&self) -> &SamHeader {
        &self.header
    }

    /// Create iterator
    pub fn iter(&mut self) -> SamIterator<'_> {
        SamIterator { reader: self }
    }
}

/// Iterator over SAM records borrowing [`SamReader`].
/// # Invariants
/// Yields parsed records until EOF.
/// # Ownership
/// Borrows reader mutably.
/// # Mutation
/// Advances underlying reader.
/// # Biological assumptions
/// None (I/O adapter).
/// # Java equivalence
/// htsjdk `SAMRecordIterator`.
pub struct SamIterator<'a> {
    reader: &'a mut SamReader,
}

impl<'a> Iterator for SamIterator<'a> {
    type Item = gatk_common::GatkResult<SamRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.read_next_record() {
            Ok(Some(record)) => Some(Ok(record)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

impl SamWriter {
    /// Create new SAM writer
    pub fn new<P: AsRef<std::path::Path>>(
        path: P,
        header: SamHeader,
    ) -> gatk_common::GatkResult<Self> {
        let file = std::fs::File::create(path)
            .map_err(|e| gatk_common::GatkError::io("Failed to create SAM file", e))?;

        Ok(Self {
            writer: Box::new(file),
            header,
        })
    }

    /// Write header
    pub fn write_header(&mut self) -> gatk_common::GatkResult<()> {
        // Write header lines
        self.writer.write_all(b"@HD\tVN:1.6\tSO:coordinate\n")?;

        // Write reference sequences
        for ref_seq in &self.header.reference_sequences {
            let line = format!("@SQ\tSN:{}\tLN:{}", ref_seq.name, ref_seq.length);
            if let Some(ref md5) = ref_seq.md5 {
                let line_with_md5 = format!("{}\tM5:{}", line, md5);
                self.writer.write_all(line_with_md5.as_bytes())?;
            } else {
                self.writer.write_all(line.as_bytes())?;
            }
            self.writer.write_all(b"\n")?;
        }

        // Write read groups
        for rg in &self.header.read_groups {
            let mut line = format!("@RG\tID:{}", rg.id);
            if let Some(ref desc) = rg.description {
                line.push_str(&format!("\tDS:{}", desc));
            }
            if let Some(ref platform) = rg.platform {
                line.push_str(&format!("\tPL:{}", platform));
            }
            if let Some(ref sample) = rg.sample {
                line.push_str(&format!("\tSM:{}", sample));
            }
            self.writer.write_all(line.as_bytes())?;
            self.writer.write_all(b"\n")?;
        }

        // Write programs
        for prog in &self.header.programs {
            let mut line = format!("@PG\tID:{}", prog.id);
            if let Some(ref name) = prog.name {
                line.push_str(&format!("\tPN:{}", name));
            }
            if let Some(ref version) = prog.version {
                line.push_str(&format!("\tVN:{}", version));
            }
            self.writer.write_all(line.as_bytes())?;
            self.writer.write_all(b"\n")?;
        }

        // Write comments
        for comment in &self.header.comments {
            self.writer
                .write_all(format!("@CO\t{}\n", comment).as_bytes())?;
        }

        Ok(())
    }

    /// Write record
    pub fn write_record(&mut self, record: &SamRecord) -> gatk_common::GatkResult<()> {
        let line = self.format_record(record);
        self.writer.write_all(line.as_bytes())?;
        self.writer.write_all(b"\n")?;
        Ok(())
    }

    pub fn write_records(&mut self, records: &[SamRecord]) -> gatk_common::GatkResult<()> {
        for record in records {
            self.write_record(record)?;
        }
        Ok(())
    }

    /// Format record as SAM line
    fn format_record(&self, record: &SamRecord) -> String {
        let mut line = String::new();

        // Required fields
        line.push_str(&record.qname);
        line.push('\t');
        line.push_str(&record.flag.to_string());
        line.push('\t');
        line.push_str(&record.rname);
        line.push('\t');
        line.push_str(&record.pos.to_string());
        line.push('\t');
        line.push_str(&record.mapq.to_string());
        line.push('\t');
        line.push_str(&record.cigar);
        line.push('\t');
        line.push_str(&record.rnext);
        line.push('\t');
        line.push_str(&record.pnext.to_string());
        line.push('\t');
        line.push_str(&record.tlen.to_string());
        line.push('\t');
        line.push_str(&record.seq);
        line.push('\t');
        line.push_str(&record.qual);

        // Optional fields
        for field in &record.optional {
            line.push('\t');
            line.push_str(&match field {
                OptionalField::Char(tag, value) => format!("{}:A:{}", tag, value),
                OptionalField::Int(tag, value) => format!("{}:i:{}", tag, value),
                OptionalField::Float(tag, value) => format!("{}:f:{}", tag, value),
                OptionalField::String(tag, value) => format!("{}:Z:{}", tag, value),
                OptionalField::Hex(tag, value) => {
                    let hex_str: String = value.iter().map(|b| format!("{:02X}", b)).collect();
                    format!("{}:H:{}", tag, hex_str)
                }
                OptionalField::Array(tag, values) => {
                    let values_str: String = values.iter().map(|v| v.to_string()).collect();
                    format!("{}:B:i,{}", tag, values_str)
                }
            });
        }

        line
    }
}

/// Buffered SAM reader implementation
struct BufferedSamReader {
    reader: std::io::BufReader<std::fs::File>,
    line_buffer: String,
}

impl BufferedSamReader {
    fn new(file: std::fs::File) -> gatk_common::GatkResult<Self> {
        Ok(Self {
            reader: std::io::BufReader::new(file),
            line_buffer: String::new(),
        })
    }

    fn read_header(&mut self) -> gatk_common::GatkResult<SamHeader> {
        let mut header = SamHeader::default();

        // Read header lines
        loop {
            self.line_buffer.clear();
            let bytes_read = self
                .reader
                .read_line(&mut self.line_buffer)
                .map_err(|e| gatk_common::GatkError::io("Failed to read SAM header line", e))?;

            if bytes_read == 0 {
                break;
            }

            let line = self.line_buffer.trim();
            if line.is_empty() {
                continue;
            }

            if !line.starts_with('@') {
                break; // End of header
            }

            // Parse header line
            if line.starts_with("@HD") {
                header.sort_order = Self::parse_sort_order(line);
                header.grouping = Self::parse_grouping(line);
            } else if line.starts_with("@SQ") {
                header
                    .reference_sequences
                    .push(Self::parse_reference_sequence(line));
            } else if line.starts_with("@RG") {
                header.read_groups.push(Self::parse_read_group(line));
            } else if line.starts_with("@PG") {
                header.programs.push(Self::parse_program(line));
            } else if let Some(stripped) = line.strip_prefix("@CO") {
                header.comments.push(stripped.trim().to_string());
            }
        }

        Ok(header)
    }

    fn parse_sort_order(line: &str) -> Option<String> {
        for field in line.split('\t').skip(1) {
            if let Some((key, value)) = field.split_once(':') {
                if key == "SO" {
                    return Some(value.to_string());
                }
            }
        }
        None
    }

    fn parse_grouping(line: &str) -> Option<String> {
        for field in line.split('\t').skip(1) {
            if let Some((key, value)) = field.split_once(':') {
                if key == "GO" {
                    return Some(value.to_string());
                }
            }
        }
        None
    }

    fn parse_reference_sequence(line: &str) -> ReferenceSequence {
        let mut ref_seq = ReferenceSequence {
            name: String::new(),
            length: 0,
            md5: None,
            assembly: None,
            uri: None,
            species: None,
        };

        for field in line.split('\t').skip(1) {
            if let Some((key, value)) = field.split_once(':') {
                match key {
                    "SN" => ref_seq.name = value.to_string(),
                    "LN" => ref_seq.length = value.parse().unwrap_or(0),
                    "M5" => ref_seq.md5 = Some(value.to_string()),
                    "AS" => ref_seq.assembly = Some(value.to_string()),
                    "UR" => ref_seq.uri = Some(value.to_string()),
                    "SP" => ref_seq.species = Some(value.to_string()),
                    _ => {}
                }
            }
        }

        ref_seq
    }

    fn parse_read_group(line: &str) -> ReadGroup {
        let mut rg = ReadGroup {
            id: String::new(),
            description: None,
            flow_order: None,
            key_sequence: None,
            library: None,
            platform_unit: None,
            platform: None,
            sample: None,
        };

        for field in line.split('\t').skip(1) {
            if let Some((key, value)) = field.split_once(':') {
                match key {
                    "ID" => rg.id = value.to_string(),
                    "DS" => rg.description = Some(value.to_string()),
                    "FO" => rg.flow_order = Some(value.to_string()),
                    "KS" => rg.key_sequence = Some(value.to_string()),
                    "LB" => rg.library = Some(value.to_string()),
                    "PU" => rg.platform_unit = Some(value.to_string()),
                    "PL" => rg.platform = Some(value.to_string()),
                    "SM" => rg.sample = Some(value.to_string()),
                    _ => {}
                }
            }
        }

        rg
    }

    fn parse_program(line: &str) -> Program {
        let mut prog = Program {
            id: String::new(),
            name: None,
            command_line: None,
            previous_id: None,
            version: None,
        };

        for field in line.split('\t').skip(1) {
            if let Some((key, value)) = field.split_once(':') {
                match key {
                    "ID" => prog.id = value.to_string(),
                    "PN" => prog.name = Some(value.to_string()),
                    "CL" => prog.command_line = Some(value.to_string()),
                    "PP" => prog.previous_id = Some(value.to_string()),
                    "VN" => prog.version = Some(value.to_string()),
                    _ => {}
                }
            }
        }

        prog
    }

    fn read_next_record(&mut self) -> gatk_common::GatkResult<Option<SamRecord>> {
        // Check if we have a line from header parsing
        if !self.line_buffer.is_empty() && !self.line_buffer.starts_with('@') {
            let record = Self::parse_record(self.line_buffer.trim_end())?;
            self.line_buffer.clear();
            return Ok(Some(record));
        }

        // Read next line
        self.line_buffer.clear();
        let bytes_read = self
            .reader
            .read_line(&mut self.line_buffer)
            .map_err(|e| gatk_common::GatkError::io("Failed to read SAM record line", e))?;

        if bytes_read == 0 {
            return Ok(None);
        }

        let line = self.line_buffer.trim();
        if line.is_empty() {
            return self.read_next_record(); // Skip empty lines
        }

        let record = Self::parse_record(line)?;
        Ok(Some(record))
    }

    fn parse_record(line: &str) -> gatk_common::GatkResult<SamRecord> {
        let fields: Vec<&str> = line.split('\t').collect();

        if fields.len() < 11 {
            return Err(gatk_common::GatkError::generic(
                "SAM record must have at least 11 fields",
            ));
        }

        // Parse required fields
        let qname = fields[0].to_string();
        let flag = fields[1].parse().unwrap_or(0);
        let rname = fields[2].to_string();
        let pos = fields[3].parse().unwrap_or(-1);
        let mapq = fields[4].parse().unwrap_or(255);
        let cigar = fields[5].to_string();
        let rnext = fields[6].to_string();
        let pnext = fields[7].parse().unwrap_or(0);
        let tlen = fields[8].parse().unwrap_or(0);
        let seq = fields[9].to_string();
        let qual = fields[10].to_string();

        // Parse optional fields
        let mut optional = Vec::new();
        for field in fields.iter().skip(11) {
            if let Some(parsed_field) = Self::parse_optional_field(field) {
                optional.push(parsed_field);
            }
        }

        Ok(SamRecord {
            qname,
            flag,
            rname,
            pos,
            mapq,
            cigar,
            rnext,
            pnext,
            tlen,
            seq,
            qual,
            optional,
        })
    }

    fn parse_optional_field(field: &str) -> Option<OptionalField> {
        let parts: Vec<&str> = field.split(':').collect();
        if parts.len() < 3 {
            return None;
        }

        let tag = parts[0].to_string();
        let type_char = parts[1].chars().next()?;

        match type_char {
            'A' => {
                if parts.len() >= 3 {
                    Some(OptionalField::Char(tag, parts[2].chars().next()?))
                } else {
                    None
                }
            }
            'i' => {
                if parts.len() >= 3 {
                    parts[2]
                        .parse()
                        .ok()
                        .map(|value| OptionalField::Int(tag, value))
                } else {
                    None
                }
            }
            'f' => {
                if parts.len() >= 3 {
                    parts[2]
                        .parse()
                        .ok()
                        .map(|value| OptionalField::Float(tag, value))
                } else {
                    None
                }
            }
            'Z' => {
                if parts.len() >= 3 {
                    Some(OptionalField::String(tag, parts[2..].join(":")))
                } else {
                    None
                }
            }
            'H' => {
                if parts.len() >= 3 {
                    let hex_str = parts[2..].join("");
                    let mut hex_bytes = Vec::new();
                    for chunk in hex_str.as_bytes().chunks(2) {
                        if chunk.len() == 2 {
                            if let Ok(s) = std::str::from_utf8(chunk) {
                                if let Ok(byte) = u8::from_str_radix(s, 16) {
                                    hex_bytes.push(byte);
                                }
                            }
                        }
                    }
                    Some(OptionalField::Hex(tag, hex_bytes))
                } else {
                    None
                }
            }
            'B' => {
                if parts.len() >= 4 {
                    let values: Result<Vec<i32>, _> =
                        parts[3..].iter().map(|s| s.parse()).collect();
                    values.ok().map(|vals| OptionalField::Array(tag, vals))
                } else {
                    None
                }
            }
            _ => None,
        }
    }
}

trait SamReaderBackend {
    #[allow(dead_code)]
    fn read_header(&mut self) -> gatk_common::GatkResult<SamHeader>;
    fn read_next_record(&mut self) -> gatk_common::GatkResult<Option<SamRecord>>;
}

impl SamReaderBackend for BufferedSamReader {
    fn read_header(&mut self) -> gatk_common::GatkResult<SamHeader> {
        self.read_header()
    }

    fn read_next_record(&mut self) -> gatk_common::GatkResult<Option<SamRecord>> {
        self.read_next_record()
    }
}

impl SamRecord {
    /// Check if read is properly paired
    pub fn is_paired(&self) -> bool {
        (self.flag & 0x1) != 0
    }

    /// Check if read is properly mapped
    pub fn is_proper_pair(&self) -> bool {
        (self.flag & 0x2) != 0
    }

    /// Check if read is unmapped
    pub fn is_unmapped(&self) -> bool {
        (self.flag & 0x4) != 0
    }

    /// Check if mate is unmapped
    pub fn is_mate_unmapped(&self) -> bool {
        (self.flag & 0x8) != 0
    }

    /// Check if read is reverse strand
    pub fn is_reverse_strand(&self) -> bool {
        (self.flag & 0x10) != 0
    }

    /// Check if mate is reverse strand
    pub fn is_mate_reverse_strand(&self) -> bool {
        (self.flag & 0x20) != 0
    }

    /// Check if read is first in pair
    pub fn is_first_in_pair(&self) -> bool {
        (self.flag & 0x40) != 0
    }

    /// Check if read is second in pair
    pub fn is_second_in_pair(&self) -> bool {
        (self.flag & 0x80) != 0
    }

    /// Check if read is secondary alignment
    pub fn is_secondary(&self) -> bool {
        (self.flag & 0x100) != 0
    }

    /// Check if read fails QC
    pub fn is_qc_fail(&self) -> bool {
        (self.flag & 0x200) != 0
    }

    /// Check if read is duplicate
    pub fn is_duplicate(&self) -> bool {
        (self.flag & 0x400) != 0
    }

    /// Check if read is supplementary
    pub fn is_supplementary(&self) -> bool {
        (self.flag & 0x800) != 0
    }

    /// Parse CIGAR string to operations
    pub fn parse_cigar(&self) -> Vec<CigarOp> {
        parse_cigar_str(&self.cigar)
    }

    /// Get read length from CIGAR
    pub fn read_length(&self) -> u32 {
        self.parse_cigar()
            .iter()
            .map(|op| match op {
                CigarOp::Match(len)
                | CigarOp::Insertion(len)
                | CigarOp::SoftClip(len)
                | CigarOp::HardClip(len)
                | CigarOp::Pad(len)
                | CigarOp::Equal(len)
                | CigarOp::Diff(len) => *len,
                CigarOp::Deletion(_) | CigarOp::RefSkip(_) => 0,
            })
            .sum()
    }

    /// Get reference span from CIGAR
    pub fn reference_span(&self) -> u32 {
        self.parse_cigar()
            .iter()
            .map(|op| match op {
                CigarOp::Match(len)
                | CigarOp::Deletion(len)
                | CigarOp::RefSkip(len)
                | CigarOp::Equal(len)
                | CigarOp::Diff(len) => *len,
                CigarOp::Insertion(_)
                | CigarOp::SoftClip(_)
                | CigarOp::HardClip(_)
                | CigarOp::Pad(_) => 0,
            })
            .sum()
    }

    /// Get average quality score
    pub fn average_quality(&self) -> f64 {
        if self.qual.is_empty() || self.qual == "*" {
            return 0.0;
        }

        self.qual
            .chars()
            .filter_map(|c| c.to_digit(10))
            .map(|d| d as f64)
            .sum::<f64>()
            / self.qual.len() as f64
    }

    /// Get GC content
    pub fn gc_content(&self) -> f64 {
        if self.seq.is_empty() || self.seq == "*" {
            return 0.0;
        }

        let gc_count = self
            .seq
            .chars()
            .filter(|&c| matches!(c, 'G' | 'C' | 'g' | 'c'))
            .count();
        gc_count as f64 / self.seq.len() as f64
    }

    /// Get optional field by tag
    pub fn get_optional_field(&self, tag: &str) -> Option<&OptionalField> {
        self.optional.iter().find(|field| match field {
            OptionalField::Char(t, _) => t == tag,
            OptionalField::Int(t, _) => t == tag,
            OptionalField::Float(t, _) => t == tag,
            OptionalField::String(t, _) => t == tag,
            OptionalField::Hex(t, _) => t == tag,
            OptionalField::Array(t, _) => t == tag,
        })
    }

    /// Get optional field value by tag
    pub fn get_optional_value(&self, tag: &str) -> Option<String> {
        self.get_optional_field(tag).map(|field| match field {
            OptionalField::Char(_, value) => value.to_string(),
            OptionalField::Int(_, value) => value.to_string(),
            OptionalField::Float(_, value) => value.to_string(),
            OptionalField::String(_, value) => value.clone(),
            OptionalField::Hex(_, value) => value.iter().map(|b| format!("{:02X}", b)).collect(),
            OptionalField::Array(_, values) => values.iter().map(|v| v.to_string()).collect(),
        })
    }
}
