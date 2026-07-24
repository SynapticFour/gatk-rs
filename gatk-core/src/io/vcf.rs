//! VCF file parser and writer for GATK-RS
//! This module provides efficient parsing and writing of VCF files
//! with support for large files, streaming operations, and
//! full VCF specification compliance.

#![warn(clippy::unwrap_used, clippy::expect_used)]

use std::io::{BufRead, Read, Write};

/// VCF file reader with streaming support.
/// # Invariants
/// Header parsed at open; body records follow VCF 4.x column layout.
/// # Ownership
/// Owns backend reader and [`VcfHeader`].
/// # Mutation
/// Sequential record reads mutate cursor.
/// # Biological assumptions
/// Variant calls encoded as VCF rows with INFO/FORMAT fields.
/// # Java equivalence
/// htsjdk `VCFCodec` / `VariantContextBuilder` input side.
pub struct VcfReader {
    reader: Box<dyn VcfReaderBackend>,
    header: VcfHeader,
}

/// VCF file writer with efficient output.
/// # Invariants
/// Header written before variant records.
/// # Ownership
/// Owns writer trait object and header snapshot.
/// # Mutation
/// Append-only writes.
/// # Biological assumptions
/// Emits variant records for downstream GATK tools.
/// # Java equivalence
/// htsjdk `VariantContextWriter`.
pub struct VcfWriter {
    writer: Box<dyn Write>,
    header: VcfHeader,
}

/// VCF header containing metadata, samples, and format fields.
/// # Invariants
/// Declares contigs, INFO/FORMAT/FILTER dictionaries used by body rows.
/// # Ownership
/// Owns sample list and field metadata vectors.
/// # Mutation
/// Built during parse or writer initialization.
/// # Biological assumptions
/// Describes variant file schema and cohort sample names.
/// # Java equivalence
/// htsjdk `VCFHeader`.
#[derive(Debug, Clone)]
pub struct VcfHeader {
    pub file_format: String,
    pub source: Option<String>,
    pub reference: Option<String>,
    pub contigs: Vec<Contig>,
    pub samples: Vec<String>,
    pub info_fields: Vec<InfoField>,
    pub format_fields: Vec<FormatField>,
    pub filter_fields: Vec<FilterField>,
    pub other_headers: Vec<(String, String)>,
}

/// Contig metadata from VCF `##contig` header lines.
/// # Invariants
/// `id` matches CHROM values in body records.
/// # Ownership
/// Owns contig id and optional INFO strings.
/// # Mutation
/// Immutable after header parse.
/// # Biological assumptions
/// Reference sequence dictionary entry for variants.
/// # Java equivalence
/// htsjdk `VCFContigLine`.
#[derive(Debug, Clone)]
pub struct Contig {
    pub id: String,
    pub length: Option<u64>,
    pub md5: Option<String>,
    pub assembly: Option<String>,
    pub species: Option<String>,
    pub uri: Option<String>,
}

/// VCF INFO field definition from header.
/// # Invariants
/// `number` and `type_field` follow VCF header grammar.
/// # Ownership
/// Owns id/description strings.
/// # Mutation
/// Schema metadata; immutable in parsed headers.
/// # Biological assumptions
/// Describes per-variant annotations (AF, DP, etc.).
/// # Java equivalence
/// htsjdk `VCFInfoHeaderLine`.
#[derive(Debug, Clone)]
pub struct InfoField {
    pub id: String,
    pub number: String,
    pub type_field: String,
    pub description: String,
    pub source: Option<String>,
    pub version: Option<String>,
}

/// VCF FORMAT field definition from header.
/// # Invariants
/// Applies to genotype columns in body rows.
/// # Ownership
/// Owns field metadata strings.
/// # Mutation
/// Immutable header schema.
/// # Biological assumptions
/// Per-sample genotype annotations (GT, GQ, PL, AD).
/// # Java equivalence
/// htsjdk `VCFFormatHeaderLine`.
#[derive(Debug, Clone)]
pub struct FormatField {
    pub id: String,
    pub number: String,
    pub type_field: String,
    pub description: String,
}

/// VCF FILTER field definition from header.
/// # Invariants
/// FILTER ids referenced in record FILTER column must appear here or as `.`.
/// # Ownership
/// Owns id and description.
/// # Mutation
/// Immutable header metadata.
/// # Biological assumptions
/// Variant QC failure codes.
/// # Java equivalence
/// htsjdk `VCFFilterHeaderLine`.
#[derive(Debug, Clone)]
pub struct FilterField {
    pub id: String,
    pub description: String,
}

/// VCF variant record (single site or symbolic allele row).
/// # Invariants
/// `position` is 1-based; `alternate` may be empty for monomorphic sites depending on caller.
/// # Ownership
/// Owns allele strings, INFO values, and per-sample [`SampleData`].
/// # Mutation
/// Public fields for transform/filter stages.
/// # Biological assumptions
/// One VCF row representing variant hypotheses and sample genotypes.
/// # Java equivalence
/// htsjdk `VariantContext`.
#[derive(Debug, Clone)]
pub struct VcfRecord {
    pub chromosome: String,
    pub position: u64,
    pub id: String,
    pub reference: String,
    pub alternate: Vec<String>,
    pub quality: Option<f64>,
    pub filter: Vec<String>,
    pub info: Vec<InfoValue>,
    pub format: Vec<String>,
    pub samples: Vec<SampleData>,
}

/// Parsed INFO column value with SAM-style typing.
/// # Invariants
/// Variant discriminant matches VCF INFO type declaration.
/// # Ownership
/// Owns tag id and typed payload vectors.
/// # Mutation
/// Immutable parsed value.
/// # Biological assumptions
/// Variant-level annotation attached to VCF row.
/// # Java equivalence
/// htsjdk INFO attribute representations on `VariantContext`.
#[derive(Debug, Clone)]
pub enum InfoValue {
    Flag(String),
    Integer(String, Vec<i32>),
    Float(String, Vec<f64>),
    String(String, Vec<String>),
    Character(String, Vec<char>),
}

/// Per-sample FORMAT column values for one VCF row.
/// # Invariants
/// Fields optional depending on FORMAT keys present in the row.
/// # Ownership
/// Owns genotype and numeric vectors; clone per sample column.
/// # Mutation
/// Public optional fields for parsing pipelines.
/// # Biological assumptions
/// One sample's genotypes and supporting evidence at a locus.
/// # Java equivalence
/// htsjdk `Genotype` / `GenotypeBuilder` fields.
#[derive(Debug, Clone)]
pub struct SampleData {
    pub gt: Option<Genotype>,
    pub gq: Option<f64>,
    pub dp: Option<u32>,
    pub ad: Option<Vec<u32>>,
    pub pl: Option<Vec<u32>>,
    pub other: Vec<(String, String)>,
}

/// VCF genotype allele indices with phasing flag.
/// # Invariants
/// Allele indices follow VCF GT encoding (0 = REF, 1+ = ALT); `-1` or missing handled by caller.
/// # Ownership
/// Owns allele index vector.
/// # Mutation
/// Immutable once parsed from GT string.
/// # Biological assumptions
/// Unphased or phased genotype call at a locus.
/// # Java equivalence
/// htsjdk `Genotype`.
#[derive(Debug, Clone)]
pub struct Genotype {
    pub alleles: Vec<i32>,
    pub phased: bool,
}

impl VcfReader {
    /// Create new VCF reader from file
    pub fn from_file<P: AsRef<std::path::Path>>(path: P) -> gatk_common::GatkResult<Self> {
        let path_ref = path.as_ref();
        let is_bgzip = path_ref
            .extension()
            .and_then(|e| e.to_str())
            .is_some_and(|e| e.eq_ignore_ascii_case("gz"));

        let file = std::fs::File::open(path_ref)
            .map_err(|e| gatk_common::GatkError::io("Failed to open VCF file", e))?;

        let data: Box<dyn Read> = if is_bgzip {
            Box::new(flate2::read::MultiGzDecoder::new(file))
        } else {
            Box::new(file)
        };

        let mut reader = Box::new(BufferedVcfReader::new(data));
        let header = reader.read_header()?;

        Ok(Self { reader, header })
    }

    /// Read next variant record
    pub fn read_next_record(&mut self) -> gatk_common::GatkResult<Option<VcfRecord>> {
        self.reader.read_next_record()
    }

    /// Read all records
    pub fn read_all_records(&mut self) -> gatk_common::GatkResult<Vec<VcfRecord>> {
        let mut records = Vec::new();
        while let Some(record) = self.read_next_record()? {
            records.push(record);
        }
        Ok(records)
    }

    /// Get header
    pub fn header(&self) -> &VcfHeader {
        &self.header
    }

    /// Create iterator
    pub fn iter(&mut self) -> VcfIterator<'_> {
        VcfIterator { reader: self }
    }
}

/// Iterator borrowing a [`VcfReader`] for streaming variants.
/// # Invariants
/// Propagates parse errors as `GatkResult` items.
/// # Ownership
/// Holds `&mut VcfReader`.
/// # Mutation
/// Advances reader on each `next`.
/// # Biological assumptions
/// None (I/O adapter).
/// # Java equivalence
/// htsjdk `CloseableIterator<VariantContext>`.
pub struct VcfIterator<'a> {
    reader: &'a mut VcfReader,
}

impl<'a> Iterator for VcfIterator<'a> {
    type Item = gatk_common::GatkResult<VcfRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.read_next_record() {
            Ok(Some(record)) => Some(Ok(record)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

impl VcfWriter {
    /// Create new VCF writer
    pub fn new<P: AsRef<std::path::Path>>(
        path: P,
        header: VcfHeader,
    ) -> gatk_common::GatkResult<Self> {
        let file = std::fs::File::create(path)
            .map_err(|e| gatk_common::GatkError::io("Failed to create VCF file", e))?;

        Ok(Self {
            writer: Box::new(file),
            header,
        })
    }

    /// Write header
    pub fn write_header(&mut self) -> gatk_common::GatkResult<()> {
        // Write file format
        self.writer
            .write_all(format!("##fileformat={}\n", self.header.file_format).as_bytes())?;

        // Write source
        if let Some(ref source) = self.header.source {
            self.writer
                .write_all(format!("##source={}\n", source).as_bytes())?;
        }

        // Write reference
        if let Some(ref reference) = self.header.reference {
            self.writer
                .write_all(format!("##reference={}\n", reference).as_bytes())?;
        }

        // Write contigs
        for contig in &self.header.contigs {
            let mut contig_line = format!("##contig=<ID={}", contig.id);
            if let Some(ref length) = contig.length {
                contig_line.push_str(&format!(",length={}", length));
            }
            if let Some(ref md5) = contig.md5 {
                contig_line.push_str(&format!(",MD5={}", md5));
            }
            if let Some(ref assembly) = contig.assembly {
                contig_line.push_str(&format!(",assembly={}", assembly));
            }
            contig_line.push_str(">\n");
            self.writer.write_all(contig_line.as_bytes())?;
        }

        // Write INFO fields
        for info in &self.header.info_fields {
            self.writer.write_all(
                format!(
                    "##INFO=<ID={},Number={},Type={},Description=\"{}\">\n",
                    info.id, info.number, info.type_field, info.description
                )
                .as_bytes(),
            )?;
        }

        // Write FORMAT fields
        for format in &self.header.format_fields {
            self.writer.write_all(
                format!(
                    "##FORMAT=<ID={},Number={},Type={},Description=\"{}\">\n",
                    format.id, format.number, format.type_field, format.description
                )
                .as_bytes(),
            )?;
        }

        // Write FILTER fields
        for filter in &self.header.filter_fields {
            self.writer.write_all(
                format!(
                    "##FILTER=<ID={},Description=\"{}\">\n",
                    filter.id, filter.description
                )
                .as_bytes(),
            )?;
        }

        // Write other headers
        for (key, value) in &self.header.other_headers {
            self.writer
                .write_all(format!("##{}={}\n", key, value).as_bytes())?;
        }

        // Write column header
        self.writer
            .write_all(b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")?;
        if !self.header.samples.is_empty() {
            self.writer.write_all(b"\tFORMAT")?;
            for sample in &self.header.samples {
                self.writer.write_all(b"\t")?;
                self.writer.write_all(sample.as_bytes())?;
            }
        }
        self.writer.write_all(b"\n")?;

        Ok(())
    }

    /// Write record
    pub fn write_record(&mut self, record: &VcfRecord) -> gatk_common::GatkResult<()> {
        let line = self.format_record(record);
        self.writer.write_all(line.as_bytes())?;
        self.writer.write_all(b"\n")?;
        Ok(())
    }

    pub fn write_records(&mut self, records: &[VcfRecord]) -> gatk_common::GatkResult<()> {
        for record in records {
            self.write_record(record)?;
        }
        Ok(())
    }

    /// Format record as VCF line
    fn format_record(&self, record: &VcfRecord) -> String {
        let mut line = String::new();

        // Required fields
        line.push_str(&record.chromosome);
        line.push('\t');
        line.push_str(&record.position.to_string());
        line.push('\t');
        line.push_str(&record.id);
        line.push('\t');
        line.push_str(&record.reference);
        line.push('\t');
        line.push_str(&record.alternate.join(","));
        line.push('\t');
        line.push_str(
            &record
                .quality
                .map_or(".".to_string(), |q| format!("{:.2}", q)),
        );
        line.push('\t');
        line.push_str(&record.filter.join(";"));
        line.push('\t');
        line.push_str(&self.format_info(&record.info));

        // Sample data
        if !record.samples.is_empty() {
            line.push('\t');
            line.push_str(&record.format.join(":"));
            for sample in &record.samples {
                line.push('\t');
                line.push_str(&self.format_sample(sample, &record.format));
            }
        }

        line
    }

    /// Format INFO fields
    fn format_info(&self, info: &[InfoValue]) -> String {
        if info.is_empty() {
            return ".".to_string();
        }

        info.iter()
            .map(|info_val| match info_val {
                InfoValue::Flag(id) => id.clone(),
                InfoValue::Integer(id, values) => {
                    if values.len() == 1 {
                        format!("{}={}", id, values[0])
                    } else {
                        format!(
                            "{}={}",
                            id,
                            values
                                .iter()
                                .map(|v| v.to_string())
                                .collect::<Vec<_>>()
                                .join(",")
                        )
                    }
                }
                InfoValue::Float(id, values) => {
                    if values.len() == 1 {
                        format!("{}={}", id, values[0])
                    } else {
                        format!(
                            "{}={}",
                            id,
                            values
                                .iter()
                                .map(|v| format!("{:.3}", v))
                                .collect::<Vec<_>>()
                                .join(",")
                        )
                    }
                }
                InfoValue::String(id, values) => {
                    if values.len() == 1 {
                        format!("{}={}", id, values[0])
                    } else {
                        format!("{}={}", id, values.join(","))
                    }
                }
                InfoValue::Character(id, values) => {
                    if values.len() == 1 {
                        format!("{}={}", id, values[0])
                    } else {
                        format!("{}={}", id, values.iter().collect::<String>())
                    }
                }
            })
            .collect::<Vec<_>>()
            .join(";")
    }

    /// Format sample data
    fn format_sample(&self, sample: &SampleData, format: &[String]) -> String {
        format
            .iter()
            .map(|field| match field.as_str() {
                "GT" => {
                    if let Some(ref gt) = sample.gt {
                        gt.to_string()
                    } else {
                        ".".to_string()
                    }
                }
                "GQ" => sample
                    .gq
                    .map_or_else(|| ".".to_string(), |gq| format!("{:.0}", gq)),
                "DP" => sample
                    .dp
                    .map_or_else(|| ".".to_string(), |dp| dp.to_string()),
                "AD" => sample.ad.as_ref().map_or_else(
                    || ".".to_string(),
                    |ad| {
                        ad.iter()
                            .map(|v| v.to_string())
                            .collect::<Vec<_>>()
                            .join(",")
                    },
                ),
                "PL" => sample.pl.as_ref().map_or_else(
                    || ".".to_string(),
                    |pl| {
                        pl.iter()
                            .map(|v| v.to_string())
                            .collect::<Vec<_>>()
                            .join(",")
                    },
                ),
                _ => {
                    // Look for other fields
                    sample
                        .other
                        .iter()
                        .find(|(key, _)| key == field)
                        .map(|(_, value)| value.clone())
                        .unwrap_or_else(|| ".".to_string())
                }
            })
            .collect::<Vec<_>>()
            .join(":")
    }
}

/// Buffered VCF reader implementation for streaming parse.
/// # Invariants
/// Maintains line buffer for partial records; header parsed before body.
/// # Ownership
/// Owns buffered reader state; used behind [`VcfReader`] trait object.
/// # Mutation
/// Internal buffer mutated on each line read.
/// # Biological assumptions
/// None (parser infrastructure).
/// # Java equivalence
/// htsjdk `VCFCodec` buffered decode path.
pub struct BufferedVcfReader {
    reader: std::io::BufReader<Box<dyn Read>>,
    line_buffer: String,
}

impl BufferedVcfReader {
    fn new(data: Box<dyn Read>) -> Self {
        Self {
            reader: std::io::BufReader::new(data),
            line_buffer: String::new(),
        }
    }

    fn read_header(&mut self) -> gatk_common::GatkResult<VcfHeader> {
        let mut header = VcfHeader::default();
        let mut saw_column_header = false;

        // Read header lines
        loop {
            self.line_buffer.clear();
            let bytes_read = self
                .reader
                .read_line(&mut self.line_buffer)
                .map_err(|e| gatk_common::GatkError::io("Failed to read VCF header line", e))?;

            if bytes_read == 0 {
                break;
            }

            let line = self.line_buffer.trim();
            if line.is_empty() {
                continue;
            }

            if line.starts_with("##") {
                // Parse header line
                Self::parse_header_line(line, &mut header);
            } else if line.starts_with('#') {
                // Column header line
                Self::parse_column_header(line, &mut header);
                saw_column_header = true;
                break;
            } else {
                // No header found
                break;
            }
        }

        if !saw_column_header {
            return Err(gatk_common::GatkError::file_format(
                "Invalid VCF: missing #CHROM column header (empty or non-VCF input)",
            ));
        }

        Ok(header)
    }

    fn parse_header_line(line: &str, header: &mut VcfHeader) {
        if let Some(stripped) = line.strip_prefix("##fileformat=") {
            header.file_format = stripped.to_string();
        } else if let Some(stripped) = line.strip_prefix("##source=") {
            header.source = Some(stripped.to_string());
        } else if let Some(stripped) = line.strip_prefix("##reference=") {
            header.reference = Some(stripped.to_string());
        } else if line.starts_with("##contig=") {
            header.contigs.push(Self::parse_contig(line));
        } else if line.starts_with("##INFO=") {
            header.info_fields.push(Self::parse_info_field(line));
        } else if line.starts_with("##FORMAT=") {
            header.format_fields.push(Self::parse_format_field(line));
        } else if line.starts_with("##FILTER=") {
            header.filter_fields.push(Self::parse_filter_field(line));
        } else {
            // Other header line
            if let Some((key, value)) = line[2..].split_once('=') {
                header
                    .other_headers
                    .push((key.to_string(), value.to_string()));
            }
        }
    }

    fn parse_contig(line: &str) -> Contig {
        let mut contig = Contig {
            id: String::new(),
            length: None,
            md5: None,
            assembly: None,
            species: None,
            uri: None,
        };

        // `##contig=<` is 10 bytes; drop trailing `>`.
        let content = line
            .strip_prefix("##contig=<")
            .and_then(|s| s.strip_suffix('>'))
            .unwrap_or("");
        for field in content.split(',') {
            if let Some((key, value)) = field.split_once('=') {
                match key {
                    "ID" => contig.id = value.to_string(),
                    "length" => contig.length = value.parse().ok(),
                    "MD5" => contig.md5 = Some(value.to_string()),
                    "assembly" => contig.assembly = Some(value.to_string()),
                    "species" => contig.species = Some(value.to_string()),
                    "URI" => contig.uri = Some(value.to_string()),
                    _ => {}
                }
            }
        }

        contig
    }

    fn parse_info_field(line: &str) -> InfoField {
        let mut info = InfoField {
            id: String::new(),
            number: String::new(),
            type_field: String::new(),
            description: String::new(),
            source: None,
            version: None,
        };

        let content = &line[8..line.len() - 1]; // Remove ##INFO=< and >
        for field in content.split(',') {
            if let Some((key, value)) = field.split_once('=') {
                let value = value.trim_matches('"');
                match key {
                    "ID" => info.id = value.to_string(),
                    "Number" => info.number = value.to_string(),
                    "Type" => info.type_field = value.to_string(),
                    "Description" => info.description = value.to_string(),
                    "Source" => info.source = Some(value.to_string()),
                    "Version" => info.version = Some(value.to_string()),
                    _ => {}
                }
            }
        }

        info
    }

    fn parse_format_field(line: &str) -> FormatField {
        let mut format = FormatField {
            id: String::new(),
            number: String::new(),
            type_field: String::new(),
            description: String::new(),
        };

        let content = &line[10..line.len() - 1]; // Remove ##FORMAT=< and >
        for field in content.split(',') {
            if let Some((key, value)) = field.split_once('=') {
                let value = value.trim_matches('"');
                match key {
                    "ID" => format.id = value.to_string(),
                    "Number" => format.number = value.to_string(),
                    "Type" => format.type_field = value.to_string(),
                    "Description" => format.description = value.to_string(),
                    _ => {}
                }
            }
        }

        format
    }

    fn parse_filter_field(line: &str) -> FilterField {
        let mut filter = FilterField {
            id: String::new(),
            description: String::new(),
        };

        let content = &line[11..line.len() - 1]; // Remove ##FILTER=< and >
        for field in content.split(',') {
            if let Some((key, value)) = field.split_once('=') {
                let value = value.trim_matches('"');
                match key {
                    "ID" => filter.id = value.to_string(),
                    "Description" => filter.description = value.to_string(),
                    _ => {}
                }
            }
        }

        filter
    }

    fn parse_column_header(line: &str, header: &mut VcfHeader) {
        let fields: Vec<&str> = line[1..].split('\t').collect();
        if fields.len() > 8 {
            header.samples = fields[9..].iter().map(|s| s.to_string()).collect();
        }
    }

    fn read_next_record(&mut self) -> gatk_common::GatkResult<Option<VcfRecord>> {
        // Check if we have a line from header parsing
        if !self.line_buffer.is_empty() && !self.line_buffer.starts_with('#') {
            let record = Self::parse_record(&self.line_buffer)?;
            self.line_buffer.clear();
            return Ok(Some(record));
        }

        // Read next line
        self.line_buffer.clear();
        let bytes_read = self
            .reader
            .read_line(&mut self.line_buffer)
            .map_err(|e| gatk_common::GatkError::io("Failed to read VCF record line", e))?;

        if bytes_read == 0 {
            return Ok(None);
        }

        let line = self.line_buffer.trim();
        if line.is_empty() {
            return self.read_next_record(); // Skip empty lines
        }
        if line.starts_with('#') {
            self.line_buffer.clear();
            return self.read_next_record();
        }

        let record = Self::parse_record(line)?;
        self.line_buffer.clear();
        Ok(Some(record))
    }

    fn parse_record(line: &str) -> gatk_common::GatkResult<VcfRecord> {
        let fields: Vec<&str> = line.split('\t').collect();

        if fields.len() < 8 {
            return Err(gatk_common::GatkError::generic(
                "VCF record must have at least 8 fields",
            ));
        }

        // Parse required fields
        let chromosome = fields[0].to_string();
        let position = fields[1].parse::<u64>().map_err(|_| {
            gatk_common::GatkError::generic(format!(
                "{} is not a valid start position in the VCF format",
                fields[1]
            ))
        })?;
        let id = fields[2].to_string();
        let reference = fields[3].to_string();
        let alternate = if fields[4] == "." {
            Vec::new()
        } else {
            fields[4].split(',').map(|s| s.to_string()).collect()
        };
        let quality = if fields[5] == "." {
            None
        } else {
            fields[5].parse().ok()
        };
        let filter = if fields[6] == "." {
            Vec::new()
        } else {
            fields[6].split(';').map(|s| s.to_string()).collect()
        };
        let info = if fields[7] == "." {
            Vec::new()
        } else {
            Self::parse_info(fields[7])?
        };

        // Parse format and samples
        let (format, samples) = if fields.len() > 8 {
            let format = if fields[8] == "." {
                Vec::new()
            } else {
                fields[8].split(':').map(|s| s.to_string()).collect()
            };

            let samples = if fields.len() > 9 {
                fields[9..]
                    .iter()
                    .map(|sample_str| Self::parse_sample(sample_str, &format))
                    .collect()
            } else {
                Vec::new()
            };

            (format, samples)
        } else {
            (Vec::new(), Vec::new())
        };

        Ok(VcfRecord {
            chromosome,
            position,
            id,
            reference,
            alternate,
            quality,
            filter,
            info,
            format,
            samples,
        })
    }

    pub fn parse_info(info_str: &str) -> gatk_common::GatkResult<Vec<InfoValue>> {
        if info_str.is_empty() || info_str == "." {
            return Ok(Vec::new());
        }

        let mut info_values = Vec::new();
        for field in info_str.split(';') {
            if let Some(info_val) = Self::parse_info_field_value(field) {
                info_values.push(info_val);
            }
        }

        Ok(info_values)
    }

    fn parse_info_field_value(field: &str) -> Option<InfoValue> {
        if let Some((id, value)) = field.split_once('=') {
            // Determine type based on value format
            if value.contains('.') {
                // Float
                let float_values: Vec<f64> =
                    value.split(',').filter_map(|s| s.parse().ok()).collect();
                if !float_values.is_empty() {
                    return Some(InfoValue::Float(id.to_string(), float_values));
                }
            } else if value
                .chars()
                .all(|c| c.is_ascii_digit() || c == '-' || c == ',')
            {
                // Integer
                let int_values: Vec<i32> =
                    value.split(',').filter_map(|s| s.parse().ok()).collect();
                if !int_values.is_empty() {
                    return Some(InfoValue::Integer(id.to_string(), int_values));
                }
            }

            // String
            let string_values: Vec<String> = value.split(',').map(|s| s.to_string()).collect();
            Some(InfoValue::String(id.to_string(), string_values))
        } else {
            // Flag
            Some(InfoValue::Flag(field.to_string()))
        }
    }

    pub fn parse_sample(sample_str: &str, format: &[String]) -> SampleData {
        let fields: Vec<&str> = sample_str.split(':').collect();
        let mut sample = SampleData {
            gt: None,
            gq: None,
            dp: None,
            ad: None,
            pl: None,
            other: Vec::new(),
        };

        for (i, field) in fields.iter().enumerate() {
            if i >= format.len() {
                break;
            }

            let format_field = &format[i];
            match format_field.as_str() {
                "GT" => {
                    sample.gt = Self::parse_genotype(field);
                }
                "GQ" => {
                    sample.gq = field.parse().ok();
                }
                "DP" => {
                    sample.dp = field.parse().ok();
                }
                "AD" => {
                    if *field != "." {
                        let values = field
                            .split(',')
                            .filter_map(|s| s.parse::<u32>().ok())
                            .collect::<Vec<_>>();
                        if !values.is_empty() {
                            sample.ad = Some(values);
                        }
                    }
                }
                "PL" => {
                    if *field != "." {
                        let values = field
                            .split(',')
                            .filter_map(|s| s.parse::<u32>().ok())
                            .collect::<Vec<_>>();
                        if !values.is_empty() {
                            sample.pl = Some(values);
                        }
                    }
                }
                _ => {
                    // CLONE: needed because owned element into collection.
                    sample.other.push((format_field.clone(), field.to_string()));
                }
            }
        }

        sample
    }

    fn parse_genotype(gt_str: &str) -> Option<Genotype> {
        if gt_str == "." || gt_str == "./." {
            return None;
        }

        let phased = gt_str.contains('|');
        let alleles: Vec<i32> = gt_str
            .replace('|', "/")
            .split('/')
            .filter_map(|s| s.parse().ok())
            .collect();

        if alleles.is_empty() {
            return None;
        }

        Some(Genotype { alleles, phased })
    }
}

trait VcfReaderBackend {
    #[allow(dead_code)]
    fn read_header(&mut self) -> gatk_common::GatkResult<VcfHeader>;
    fn read_next_record(&mut self) -> gatk_common::GatkResult<Option<VcfRecord>>;
}

impl VcfReaderBackend for BufferedVcfReader {
    fn read_header(&mut self) -> gatk_common::GatkResult<VcfHeader> {
        self.read_header()
    }

    fn read_next_record(&mut self) -> gatk_common::GatkResult<Option<VcfRecord>> {
        self.read_next_record()
    }
}

impl Default for VcfHeader {
    fn default() -> Self {
        Self {
            file_format: "VCFv4.2".to_string(),
            source: None,
            reference: None,
            contigs: Vec::new(),
            samples: Vec::new(),
            info_fields: Vec::new(),
            format_fields: Vec::new(),
            filter_fields: Vec::new(),
            other_headers: Vec::new(),
        }
    }
}

impl std::fmt::Display for Genotype {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let separator = if self.phased { '|' } else { '/' };
        write!(
            f,
            "{}",
            self.alleles
                .iter()
                .map(|a| a.to_string())
                .collect::<Vec<_>>()
                .join(&separator.to_string())
        )
    }
}

impl VcfRecord {
    /// Check if variant is a SNP
    pub fn is_snp(&self) -> bool {
        self.reference.len() == 1
            && self.alternate.iter().all(|alt| alt.len() == 1)
            && !self.alternate.is_empty()
    }

    /// Check if variant is an insertion
    pub fn is_insertion(&self) -> bool {
        self.alternate
            .iter()
            .any(|alt| alt.len() > self.reference.len())
    }

    /// Check if variant is a deletion
    pub fn is_deletion(&self) -> bool {
        self.alternate
            .iter()
            .any(|alt| alt.len() < self.reference.len())
    }

    /// Check if variant is an indel
    pub fn is_indel(&self) -> bool {
        self.is_insertion() || self.is_deletion()
    }

    /// Get variant type as string
    pub fn variant_type(&self) -> &'static str {
        if self.alternate.is_empty() {
            "reference"
        } else if self.is_snp() {
            "SNP"
        } else if self.is_insertion() {
            "insertion"
        } else if self.is_deletion() {
            "deletion"
        } else {
            "complex"
        }
    }

    /// Get allele count
    pub fn allele_count(&self) -> usize {
        self.alternate.len() + 1
    }

    /// Check if variant is filtered
    pub fn is_filtered(&self) -> bool {
        !self.filter.is_empty() && (self.filter.len() != 1 || self.filter[0] != "PASS")
    }

    /// Get reference allele frequency (if available)
    pub fn get_af(&self) -> Option<f64> {
        self.info.iter().find_map(|info| match info {
            InfoValue::Float(id, values) if id == "AF" => values.first().copied(),
            _ => None,
        })
    }

    /// Get allele depths (if available)
    pub fn get_ad(&self, sample_idx: usize) -> Option<&Vec<u32>> {
        self.samples.get(sample_idx).and_then(|s| s.ad.as_ref())
    }

    /// Get genotype quality (if available)
    pub fn get_gq(&self, sample_idx: usize) -> Option<f64> {
        self.samples.get(sample_idx).and_then(|s| s.gq)
    }

    /// Get depth (if available)
    pub fn get_dp(&self, sample_idx: usize) -> Option<u32> {
        self.samples.get(sample_idx).and_then(|s| s.dp)
    }

    /// Get genotype (if available)
    pub fn get_gt(&self, sample_idx: usize) -> Option<&Genotype> {
        self.samples.get(sample_idx).and_then(|s| s.gt.as_ref())
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]
mod reader_tests {
    use super::*;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::io::Write;

    #[test]
    fn from_file_reads_gzipped_vcf() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("sites.vcf.gz");
        let f = std::fs::File::create(&path).unwrap();
        let mut enc = GzEncoder::new(f, Compression::default());
        writeln!(enc, "##fileformat=VCFv4.2").unwrap();
        writeln!(enc, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
        writeln!(enc, "chr1\t9\t.\tC\tG\t60\tPASS\t.").unwrap();
        enc.finish().unwrap();

        let mut r = VcfReader::from_file(&path).unwrap();
        let recs = r.read_all_records().unwrap();
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].chromosome, "chr1");
        assert_eq!(recs[0].position, 9);
    }
}
