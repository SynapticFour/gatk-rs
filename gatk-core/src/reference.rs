//! Reference dictionary and interval parsing utilities.

#![warn(clippy::unwrap_used, clippy::expect_used)]

use crate::io::fasta::FastaReader;
use crate::io::{BamHeader, SamHeader, VcfHeader};
use gatk_common::{GatkError, GatkResult};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};
use std::fs;
use std::io::{Read, Seek, SeekFrom};
use std::path::Path;
use std::sync::{Arc, Mutex, OnceLock};

/// A contig entry in a sequence dictionary.
/// # Invariants
/// `index` is the zero-based insertion order in the dictionary.
/// `length` is the declared sequence length (bases), typically from FASTA or `@SQ` LN.
/// # Ownership
/// Owns contig `name`; clone for independent copies.
/// # Mutation
/// Public fields; normally populated only via [`SequenceDictionary::add_contig`].
/// # Biological assumptions
/// One reference sequence (chromosome/contig) with stable name and length.
/// # Java equivalence
/// Approximates htsjdk `SAMSequenceRecord` / GATK `SAMSequenceDictionary` entry.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ContigRecord {
    pub name: String,
    pub length: u64,
    pub index: usize,
}

/// Sequence dictionary derived from reference metadata.
/// # Invariants
/// Contig names are unique; `by_name` stays consistent with `contigs` order.
/// Lookup tolerates optional `chr` prefix mismatch via [`SequenceDictionary::contig`].
/// # Ownership
/// Owns contig records; cheap to [`Clone`] for read-only sharing across threads.
/// # Mutation
/// Mutate via `add_contig` and loaders; internal maps updated together.
/// # Biological assumptions
/// Reference assembly metadata; coordinates elsewhere are 1-based against these lengths.
/// # Java equivalence
/// Approximates htsjdk `SAMSequenceDictionary` / GATK `ReferenceDictionary`.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct SequenceDictionary {
    contigs: Vec<ContigRecord>,
    by_name: HashMap<String, usize>,
}

/// Parsed interval specification (1-based, inclusive).
/// # Invariants
/// `start`/`end`, when present, are 1-based inclusive and should satisfy `start <= end`.
/// Omitted bounds mean open-ended contig span at parse time.
/// # Ownership
/// Owns contig name string; clone for interval lists.
/// # Mutation
/// Public fields for pipeline staging; validate before engine use.
/// # Biological assumptions
/// GATK-style `-L` interval: one contig, closed interval on reference.
/// # Java equivalence
/// Approximates GATK interval parsing / `Interval` CLI specs (not a single Java type).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct IntervalSpec {
    pub contig: String,
    pub start: Option<u64>,
    pub end: Option<u64>,
}

impl SequenceDictionary {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_contig(&mut self, name: String, length: u64) {
        let index = self.contigs.len();
        // CLONE: needed because owned HashMap/BTree/HashSet key or value.
        self.by_name.insert(name.clone(), index);
        self.contigs.push(ContigRecord {
            name,
            length,
            index,
        });
    }

    /// Build a sequence dictionary from a FASTA path.
    ///
    /// When a samtools `.fai` sits beside the FASTA, contig names/lengths are read from the
    /// index only — we must **not** mmap/scan multi-gigabase FASTA bodies just for LN values
    /// (that path previously drove ~2–3 GiB Peak-RSS before any reads were processed).
    pub fn from_fasta_path<P: AsRef<Path>>(path: P) -> GatkResult<Self> {
        let path = path.as_ref();
        if let Some(dict) = Self::try_from_fai_beside(path) {
            return Ok(dict);
        }
        let mut reader = FastaReader::from_file_buffered(path)?;
        let mut dict = Self::new();

        while let Some(sequence) = reader.read_next_sequence()? {
            dict.add_contig(sequence.name, sequence.length as u64);
        }

        Ok(dict)
    }

    /// Dictionary from `path.fai` in file order. `None` if the index is missing or unusable.
    fn try_from_fai_beside(fasta_path: &Path) -> Option<Self> {
        let mut fai_os = fasta_path.as_os_str().to_owned();
        fai_os.push(".fai");
        let text = fs::read_to_string(Path::new(&fai_os)).ok()?;
        let mut dict = Self::new();
        for line in text.lines() {
            let line = line.trim();
            if line.is_empty() {
                continue;
            }
            let mut parts = line.split('\t');
            let name = parts.next()?.to_string();
            let length: u64 = parts.next()?.parse().ok()?;
            let _offset = parts.next()?;
            let line_bases: u64 = parts.next()?.parse().ok()?;
            let line_width: u64 = parts.next()?.parse().ok()?;
            if line_bases == 0 || line_width < line_bases {
                continue;
            }
            dict.add_contig(name, length);
        }
        if dict.contig_count() == 0 {
            None
        } else {
            Some(dict)
        }
    }

    pub fn contig_count(&self) -> usize {
        self.contigs.len()
    }

    pub fn contig_names(&self) -> impl Iterator<Item = &str> {
        self.contigs.iter().map(|c| c.name.as_str())
    }

    /// Ordered contig records (FASTA order).
    pub fn contig_records(&self) -> &[ContigRecord] {
        &self.contigs
    }

    /// One closed interval per contig covering the full sequence ([1, LN]) — used when `-L` is omitted.
    pub fn whole_genome_interval_specs(&self) -> Vec<IntervalSpec> {
        self.contigs
            .iter()
            .map(|c| IntervalSpec {
                // CLONE: needed because owned contig id for output record.
                contig: c.name.clone(),
                start: Some(1),
                end: Some(c.length),
            })
            .collect()
    }

    pub fn contig(&self, name: &str) -> Option<&ContigRecord> {
        if let Some(idx) = self.by_name.get(name) {
            return self.contigs.get(*idx);
        }
        if let Some(without_chr) = name.strip_prefix("chr") {
            if let Some(idx) = self.by_name.get(without_chr) {
                return self.contigs.get(*idx);
            }
        }
        let with_chr = format!("chr{name}");
        self.by_name
            .get(&with_chr)
            .and_then(|idx| self.contigs.get(*idx))
    }

    pub fn validate_interval(&self, interval: &IntervalSpec) -> GatkResult<()> {
        let contig = self.contig(&interval.contig).ok_or_else(|| {
            GatkError::argument(format!("Unknown contig in interval: {}", interval.contig))
        })?;

        match (interval.start, interval.end) {
            (None, None) => Ok(()),
            (Some(start), Some(end)) => {
                if start == 0 {
                    return Err(GatkError::argument("Interval start must be >= 1"));
                }
                if start > end {
                    return Err(GatkError::argument("Interval start must be <= end"));
                }
                if end > contig.length {
                    return Err(GatkError::argument(format!(
                        "Interval end {} exceeds contig {} length {}",
                        end, contig.name, contig.length
                    )));
                }
                Ok(())
            }
            _ => Err(GatkError::argument(
                "Interval must specify both start and end or neither",
            )),
        }
    }

    pub fn validate_vcf_header(&self, header: &VcfHeader) -> GatkResult<()> {
        for contig in &header.contigs {
            let dict_contig = self.contig(&contig.id).ok_or_else(|| {
                GatkError::argument(format!(
                    "VCF header contig not found in reference dictionary: {}",
                    contig.id
                ))
            })?;
            if let Some(vcf_len) = contig.length {
                if vcf_len != dict_contig.length {
                    return Err(GatkError::argument(format!(
                        "VCF contig length mismatch for {}: header={}, reference={}",
                        contig.id, vcf_len, dict_contig.length
                    )));
                }
            }
        }
        Ok(())
    }

    pub fn validate_bam_header(&self, header: &BamHeader) -> GatkResult<()> {
        for contig in &header.reference_sequences {
            let dict_contig = self.contig(&contig.name).ok_or_else(|| {
                GatkError::argument(format!(
                    "BAM header contig not found in reference dictionary: {}",
                    contig.name
                ))
            })?;
            if contig.length != dict_contig.length {
                return Err(GatkError::argument(format!(
                    "BAM contig length mismatch for {}: header={}, reference={}",
                    contig.name, contig.length, dict_contig.length
                )));
            }
        }
        Ok(())
    }

    pub fn validate_sam_header(&self, header: &SamHeader) -> GatkResult<()> {
        for contig in &header.reference_sequences {
            let dict_contig = self.contig(&contig.name).ok_or_else(|| {
                GatkError::argument(format!(
                    "SAM header contig not found in reference dictionary: {}",
                    contig.name
                ))
            })?;
            if contig.length != dict_contig.length {
                return Err(GatkError::argument(format!(
                    "SAM contig length mismatch for {}: header={}, reference={}",
                    contig.name, contig.length, dict_contig.length
                )));
            }
        }
        Ok(())
    }
}

impl IntervalSpec {
    pub fn parse(input: &str) -> GatkResult<Self> {
        let text = input.trim();
        if text.is_empty() {
            return Err(GatkError::argument("Interval cannot be empty"));
        }

        if let Some((contig, range_part)) = text.split_once(':') {
            let cleaned = range_part.replace(',', "");
            let (start_s, end_s) = cleaned.split_once('-').ok_or_else(|| {
                GatkError::argument(format!("Invalid interval range syntax: {text}"))
            })?;
            let start = start_s
                .parse::<u64>()
                .map_err(|_| GatkError::argument(format!("Invalid interval start: {start_s}")))?;
            let end = end_s
                .parse::<u64>()
                .map_err(|_| GatkError::argument(format!("Invalid interval end: {end_s}")))?;
            if start == 0 {
                return Err(GatkError::argument("Interval start must be >= 1"));
            }
            if start > end {
                return Err(GatkError::argument(format!(
                    "Interval start must be <= end ({start} > {end})"
                )));
            }

            Ok(Self {
                contig: contig.to_string(),
                start: Some(start),
                end: Some(end),
            })
        } else {
            Ok(Self {
                contig: text.to_string(),
                start: None,
                end: None,
            })
        }
    }

    pub fn parse_list(input: &str) -> GatkResult<Vec<Self>> {
        let (includes, excludes) = Self::parse_include_exclude_list(input)?;
        if !excludes.is_empty() {
            return Err(GatkError::argument(
                "Exclusion intervals (^prefix) require a sequence dictionary; use parse_intervals_cli_string instead of parse_list",
            ));
        }
        if includes.is_empty() {
            return Err(GatkError::argument("No intervals specified"));
        }
        Ok(includes)
    }

    /// Semicolon-separated tokens; tokens starting with `^` are **exclusions** (GATK-style complement within the union of includes).
    pub(crate) fn parse_include_exclude_list(input: &str) -> GatkResult<(Vec<Self>, Vec<Self>)> {
        let mut includes = Vec::new();
        let mut excludes = Vec::new();
        for token in input.split(';') {
            let trimmed = token.trim();
            if trimmed.is_empty() {
                continue;
            }
            if let Some(rest) = trimmed.strip_prefix('^') {
                let rest = rest.trim();
                if rest.is_empty() {
                    return Err(GatkError::argument("Empty exclusion interval"));
                }
                excludes.push(Self::parse(rest)?);
            } else {
                includes.push(Self::parse(trimmed)?);
            }
        }
        if includes.is_empty() {
            return Err(GatkError::argument(
                "At least one include interval is required (exclusions alone are invalid)",
            ));
        }
        Ok((includes, excludes))
    }

    pub fn parse_list_file<P: AsRef<Path>>(path: P) -> GatkResult<Vec<Self>> {
        let content = fs::read_to_string(path)
            .map_err(|e| GatkError::io("Failed to read interval list file", e))?;
        let mut intervals = Vec::new();
        for line in content.lines() {
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') || trimmed.starts_with('@') {
                continue;
            }
            if trimmed.contains('\t') {
                // Picard/GATK interval_list row: contig<TAB>start<TAB>end...
                let cols: Vec<&str> = trimmed.split('\t').collect();
                if cols.len() < 3 {
                    return Err(GatkError::argument(format!(
                        "Invalid interval_list row (expected >=3 columns): {trimmed}"
                    )));
                }
                let contig = cols[0];
                let start = cols[1].parse::<u64>().map_err(|_| {
                    GatkError::argument(format!("Invalid interval_list start: {}", cols[1]))
                })?;
                let end = cols[2].parse::<u64>().map_err(|_| {
                    GatkError::argument(format!("Invalid interval_list end: {}", cols[2]))
                })?;
                intervals.push(Self {
                    contig: contig.to_string(),
                    start: Some(start),
                    end: Some(end),
                });
            } else {
                intervals.push(Self::parse(trimmed)?);
            }
        }
        if intervals.is_empty() {
            return Err(GatkError::argument(
                "No intervals found in interval list file",
            ));
        }
        Ok(intervals)
    }

    pub fn parse_list_file_with_dictionary<P: AsRef<Path>>(
        path: P,
        dictionary: &SequenceDictionary,
    ) -> GatkResult<Vec<Self>> {
        let content = fs::read_to_string(path)
            .map_err(|e| GatkError::io("Failed to read interval list file", e))?;
        let mut includes = Vec::new();
        let mut excludes = Vec::new();

        for line in content.lines() {
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }

            let line_is_exclude = trimmed.starts_with('^');
            let body = if line_is_exclude {
                trimmed[1..].trim_start()
            } else {
                trimmed
            };

            if body.starts_with("@SQ") {
                let mut name: Option<&str> = None;
                let mut length: Option<u64> = None;
                for field in body.split('\t').skip(1) {
                    if !field.contains(':') {
                        return Err(GatkError::argument(format!(
                            "Invalid @SQ field in interval list (expected KEY:VALUE): {field}"
                        )));
                    }
                    if let Some(sn) = field.strip_prefix("SN:") {
                        if name.is_some() {
                            return Err(GatkError::argument(format!(
                                "Invalid @SQ line in interval list (duplicate SN): {body}"
                            )));
                        }
                        name = Some(sn);
                    } else if let Some(ln) = field.strip_prefix("LN:") {
                        if length.is_some() {
                            return Err(GatkError::argument(format!(
                                "Invalid @SQ line in interval list (duplicate LN): {body}"
                            )));
                        }
                        length = Some(ln.parse::<u64>().map_err(|_| {
                            GatkError::argument(format!(
                                "Invalid @SQ LN value in interval list: {ln}"
                            ))
                        })?);
                    }
                }
                let sn = name.ok_or_else(|| {
                    GatkError::argument(format!(
                        "Invalid @SQ line in interval list (missing SN): {body}"
                    ))
                })?;
                let ln = length.ok_or_else(|| {
                    GatkError::argument(format!(
                        "Invalid @SQ line in interval list (missing LN): {body}"
                    ))
                })?;
                if ln == 0 {
                    return Err(GatkError::argument(format!(
                        "Invalid @SQ LN value in interval list (must be >= 1): {ln}"
                    )));
                }
                let dict_contig = dictionary.contig(sn).ok_or_else(|| {
                    GatkError::argument(format!(
                        "Interval list @SQ contig not found in reference dictionary: {sn}"
                    ))
                })?;
                if dict_contig.length != ln {
                    return Err(GatkError::argument(format!(
                        "Interval list @SQ length mismatch for {}: interval_list={}, reference={}",
                        sn, ln, dict_contig.length
                    )));
                }
                continue;
            }

            if body.starts_with('@') {
                continue;
            }

            if body.contains('\t') {
                let cols: Vec<&str> = body.split('\t').collect();
                if cols.len() < 3 {
                    return Err(GatkError::argument(format!(
                        "Invalid interval_list row (expected >=3 columns): {body}"
                    )));
                }
                let contig = cols[0];
                let start = cols[1].parse::<u64>().map_err(|_| {
                    GatkError::argument(format!("Invalid interval_list start: {}", cols[1]))
                })?;
                let end = cols[2].parse::<u64>().map_err(|_| {
                    GatkError::argument(format!("Invalid interval_list end: {}", cols[2]))
                })?;
                if start == 0 {
                    return Err(GatkError::argument("Interval start must be >= 1"));
                }
                if start > end {
                    return Err(GatkError::argument("Interval start must be <= end"));
                }
                if cols.len() >= 4 && !cols[3].is_empty() && cols[3] != "+" && cols[3] != "-" {
                    return Err(GatkError::argument(format!(
                        "Invalid interval_list strand value: {}",
                        cols[3]
                    )));
                }
                let interval = Self {
                    contig: contig.to_string(),
                    start: Some(start),
                    end: Some(end),
                };
                dictionary.validate_interval(&interval)?;
                if line_is_exclude {
                    excludes.push(interval);
                } else {
                    includes.push(interval);
                }
            } else {
                let interval = Self::parse(body)?;
                dictionary.validate_interval(&interval)?;
                if line_is_exclude {
                    excludes.push(interval);
                } else {
                    includes.push(interval);
                }
            }
        }

        if includes.is_empty() {
            return Err(GatkError::argument(
                "No include intervals found in interval list file (only exclusions are invalid)",
            ));
        }
        if excludes.is_empty() {
            Ok(includes)
        } else {
            resolve_interval_specs_includes_excludes(dictionary, &includes, &excludes)
        }
    }

    /// Resolve to 1-based inclusive `[start, end]` (whole contig when both are [`None`]).
    pub fn resolve_closed_ends(
        &self,
        dictionary: &SequenceDictionary,
    ) -> GatkResult<(String, u64, u64)> {
        let c = dictionary.contig(&self.contig).ok_or_else(|| {
            GatkError::argument(format!("Unknown contig in interval: {}", self.contig))
        })?;
        match (self.start, self.end) {
            (Some(start), Some(end)) => {
                if start == 0 {
                    return Err(GatkError::argument("Interval start must be >= 1"));
                }
                if start > end {
                    return Err(GatkError::argument("Interval start must be <= end"));
                }
                Ok((c.name.clone(), start, end))
            }
            (None, None) => Ok((c.name.clone(), 1, c.length)),
            _ => Err(GatkError::argument(
                "Interval must specify both start and end, or neither for whole contig",
            )),
        }
    }
}

fn merge_closed_intervals(mut v: Vec<(u64, u64)>) -> Vec<(u64, u64)> {
    if v.is_empty() {
        return v;
    }
    v.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
    let mut out: Vec<(u64, u64)> = Vec::new();
    let mut cur = v[0];
    for &(s, e) in v.iter().skip(1) {
        if s <= cur.1.saturating_add(1) {
            cur.1 = cur.1.max(e);
        } else {
            out.push(cur);
            cur = (s, e);
        }
    }
    out.push(cur);
    out
}

fn subtract_one_closed_segment(lo: u64, hi: u64, cuts: &[(u64, u64)]) -> Vec<(u64, u64)> {
    let mut lo = lo;
    let mut res: Vec<(u64, u64)> = Vec::new();
    for &(cs, ce) in cuts {
        if ce < lo {
            continue;
        }
        if cs > hi {
            break;
        }
        let c0 = cs.max(lo);
        let c1 = ce.min(hi);
        if lo < c0 {
            res.push((lo, c0 - 1));
        }
        lo = c1.saturating_add(1);
        if lo > hi {
            return res;
        }
    }
    if lo <= hi {
        res.push((lo, hi));
    }
    res
}

/// Merge include intervals per contig, subtract excludes, emit disjoint closed intervals in deterministic contig order.
pub fn resolve_interval_specs_includes_excludes(
    dictionary: &SequenceDictionary,
    includes: &[IntervalSpec],
    excludes: &[IntervalSpec],
) -> GatkResult<Vec<IntervalSpec>> {
    let mut inc_map: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
    for spec in includes {
        let (c, s, e) = spec.resolve_closed_ends(dictionary)?;
        inc_map.entry(c).or_default().push((s, e));
    }
    let mut exc_map: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
    for spec in excludes {
        let (c, s, e) = spec.resolve_closed_ends(dictionary)?;
        exc_map.entry(c).or_default().push((s, e));
    }

    let mut contig_keys: Vec<String> = inc_map.keys().cloned().collect();
    contig_keys.sort();

    let mut out: Vec<IntervalSpec> = Vec::new();
    for cname in contig_keys {
        // Lifetime: `cname` is owned for this iteration; `remove` moves interval vectors
        // into `merge_closed_intervals` so the per-contig working sets are not cloned.
        let merged_inc = merge_closed_intervals(inc_map.remove(&cname).unwrap_or_default());
        let merged_exc = exc_map
            .remove(&cname)
            .map(merge_closed_intervals)
            .unwrap_or_default();
        let mut after = if merged_exc.is_empty() {
            merged_inc
        } else {
            let mut tmp = Vec::new();
            for seg in merged_inc {
                tmp.extend(subtract_one_closed_segment(seg.0, seg.1, &merged_exc));
            }
            merge_closed_intervals(tmp)
        };
        // Move contig into the final IntervalSpec; clone only for earlier segments.
        if let Some((last_s, last_e)) = after.pop() {
            for (s, e) in after {
                out.push(IntervalSpec {
                    // CLONE: needed because owned contig id for output record.
                    contig: cname.clone(),
                    start: Some(s),
                    end: Some(e),
                });
            }
            out.push(IntervalSpec {
                contig: cname,
                start: Some(last_s),
                end: Some(last_e),
            });
        }
    }

    if out.is_empty() {
        return Err(GatkError::argument(
            "No positions remain after applying inclusion and exclusion intervals",
        ));
    }
    Ok(out)
}

fn uppercase_acgtn(b: u8) -> u8 {
    match b {
        b'a' => b'A',
        b'c' => b'C',
        b'g' => b'G',
        b't' => b'T',
        b'n' => b'N',
        _ => b,
    }
}

/// Samtools-compatible `.fai` entry (`name length offset linebases linewidth`).
#[derive(Debug, Clone)]
struct SamtoolsFaiEntry {
    length: u64,
    offset: u64,
    line_bases: u64,
    line_width: u64,
}

#[derive(Debug, Clone, Default)]
struct SamtoolsFai {
    /// Keyed by FASTA contig name as written in the `.fai`.
    by_name: HashMap<String, SamtoolsFaiEntry>,
}

impl SamtoolsFai {
    fn load_beside_fasta(fasta_path: &Path) -> Option<Self> {
        let mut fai_os = fasta_path.as_os_str().to_owned();
        fai_os.push(".fai");
        let fai_path = Path::new(&fai_os);
        let text = fs::read_to_string(fai_path).ok()?;
        let mut by_name = HashMap::new();
        for line in text.lines() {
            let line = line.trim();
            if line.is_empty() {
                continue;
            }
            let mut parts = line.split('\t');
            let name = parts.next()?.to_string();
            let length: u64 = parts.next()?.parse().ok()?;
            let offset: u64 = parts.next()?.parse().ok()?;
            let line_bases: u64 = parts.next()?.parse().ok()?;
            let line_width: u64 = parts.next()?.parse().ok()?;
            if line_bases == 0 || line_width < line_bases {
                continue;
            }
            by_name.insert(
                name,
                SamtoolsFaiEntry {
                    length,
                    offset,
                    line_bases,
                    line_width,
                },
            );
        }
        if by_name.is_empty() {
            None
        } else {
            Some(Self { by_name })
        }
    }

    fn entry_for_query(&self, contig: &str) -> Option<&SamtoolsFaiEntry> {
        if let Some(e) = self.by_name.get(contig) {
            return Some(e);
        }
        self.by_name
            .iter()
            .find(|(name, _)| contig_name_matches_fasta_entry(name, contig))
            .map(|(_, e)| e)
    }
}

fn read_interval_via_samtools_fai(
    fasta_path: &Path,
    entry: &SamtoolsFaiEntry,
    start_1based: u64,
    end_1based_inclusive: u64,
) -> GatkResult<Vec<u8>> {
    if end_1based_inclusive > entry.length {
        return Err(GatkError::argument(
            "Interval extends past end of contig sequence in FASTA",
        ));
    }
    let start0 = start_1based - 1;
    let end0_excl = end_1based_inclusive;
    let want = (end0_excl - start0) as usize;
    if want == 0 {
        return Ok(Vec::new());
    }
    let end_pos = end0_excl - 1;
    let start_row = start0 / entry.line_bases;
    let end_row = end_pos / entry.line_bases;
    let start_col = start0 % entry.line_bases;
    let end_col = end_pos % entry.line_bases;
    let file_start = entry.offset + start_row * entry.line_width + start_col;
    let file_end = entry.offset + end_row * entry.line_width + end_col + 1;
    let raw_len = (file_end - file_start) as usize;
    let mut raw = vec![0u8; raw_len];
    let mut file = fs::File::open(fasta_path)
        .map_err(|e| GatkError::io("Failed to open FASTA for indexed read", e))?;
    file.seek(SeekFrom::Start(file_start))
        .map_err(|e| GatkError::io("Failed to seek in indexed FASTA", e))?;
    file.read_exact(&mut raw)
        .map_err(|e| GatkError::io("Failed to read indexed FASTA bases", e))?;
    let mut out = Vec::with_capacity(want);
    for b in raw {
        if b == b'\n' || b == b'\r' {
            continue;
        }
        out.push(uppercase_acgtn(b));
        if out.len() == want {
            break;
        }
    }
    if out.len() != want {
        return Err(GatkError::argument(format!(
            "Indexed FASTA read for {} yielded {} bases, expected {want}",
            fasta_path.display(),
            out.len()
        )));
    }
    Ok(out)
}

fn read_contig_bases_sequential(fasta_path: &Path, contig: &str) -> GatkResult<Vec<u8>> {
    let mut reader = FastaReader::from_file_buffered(fasta_path)?;
    while let Some(seq) = reader.read_next_sequence()? {
        if contig_name_matches_fasta_entry(&seq.name, contig) {
            return Ok(seq.sequence.into_iter().map(uppercase_acgtn).collect());
        }
    }
    Err(GatkError::argument(format!(
        "Contig {contig} not found in reference FASTA"
    )))
}

fn load_contig_bases(
    fasta_path: &Path,
    contig: &str,
    fai: Option<&SamtoolsFai>,
) -> GatkResult<Vec<u8>> {
    if let Some(fai) = fai {
        if let Some(entry) = fai.entry_for_query(contig) {
            return read_interval_via_samtools_fai(fasta_path, entry, 1, entry.length);
        }
    }
    read_contig_bases_sequential(fasta_path, contig)
}

type ContigArcCache = Mutex<HashMap<(String, String), Arc<Vec<u8>>>>;

fn global_contig_arc_cache() -> &'static ContigArcCache {
    static CACHE: OnceLock<ContigArcCache> = OnceLock::new();
    CACHE.get_or_init(|| Mutex::new(HashMap::new()))
}

fn shared_contig_bases(
    fasta_path: &Path,
    contig: &str,
    fai: Option<&SamtoolsFai>,
) -> GatkResult<Arc<Vec<u8>>> {
    let key = (
        fasta_path.to_string_lossy().into_owned(),
        contig.to_string(),
    );
    if let Ok(guard) = global_contig_arc_cache().lock() {
        if let Some(hit) = guard.get(&key) {
            return Ok(Arc::clone(hit));
        }
    }
    let loaded = Arc::new(load_contig_bases(fasta_path, contig, fai)?);
    let mut guard = global_contig_arc_cache()
        .lock()
        .map_err(|_| GatkError::generic("reference contig cache lock poisoned"))?;
    Ok(Arc::clone(
        guard.entry(key).or_insert_with(|| Arc::clone(&loaded)),
    ))
}

fn read_fasta_interval_bytes<P: AsRef<Path>>(
    fasta_path: P,
    contig: &str,
    start_1based: u64,
    end_1based_inclusive: u64,
) -> GatkResult<Vec<u8>> {
    if start_1based == 0 {
        return Err(GatkError::argument("Genomic position must be >= 1"));
    }
    if start_1based > end_1based_inclusive {
        return Err(GatkError::argument("Interval start must be <= end"));
    }
    let path = fasta_path.as_ref();
    if let Some(fai) = SamtoolsFai::load_beside_fasta(path) {
        if let Some(entry) = fai.entry_for_query(contig) {
            return read_interval_via_samtools_fai(path, entry, start_1based, end_1based_inclusive);
        }
    }
    let s = (start_1based - 1) as usize;
    let e = (end_1based_inclusive - 1) as usize;
    let contig_bases = read_contig_bases_sequential(path, contig)?;
    let slice = contig_bases.get(s..=e).ok_or_else(|| {
        GatkError::argument("Interval extends past end of contig sequence in FASTA")
    })?;
    Ok(slice.to_vec())
}

/// Count A/C/G/T/N over merged intervals (double-count avoided by merging overlaps per contig).
pub fn count_acgtn_histogram_for_interval_specs<P: AsRef<Path>>(
    fasta_path: P,
    dictionary: &SequenceDictionary,
    interval_specs: &[IntervalSpec],
) -> GatkResult<[u64; 5]> {
    let mut counts = [0u64; 5];
    let mut by_contig: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
    for spec in interval_specs {
        let (c, s, e) = spec.resolve_closed_ends(dictionary)?;
        by_contig.entry(c).or_default().push((s, e));
    }
    for (cname, segs) in by_contig {
        let merged = merge_closed_intervals(segs);
        for (s, e) in merged {
            let bytes = read_fasta_interval_bytes(&fasta_path, &cname, s, e)?;
            for b in bytes {
                let idx = match b {
                    b'A' => 0,
                    b'C' => 1,
                    b'G' => 2,
                    b'T' => 3,
                    b'N' => 4,
                    _ => 4,
                };
                counts[idx] += 1;
            }
        }
    }
    Ok(counts)
}

/// Contig-level reference cache (R4-1): load each contig once (via `.fai` when present), serve windows by slice.
/// Contig bytes are retained in a process-wide `Arc` map so short-lived caches (per `callRegion`)
/// do not re-read the FASTA. Local eviction drops the **smallest** contig name (BTreeMap order).
/// # Invariants
/// Cached sequences are uppercase; window requests are 1-based inclusive and clipped to contig length.
/// At most `capacity` contigs retained locally; eviction is deterministic by name order.
/// # Ownership
/// Owns `fasta_path` and `Arc<Vec<u8>>` contig blobs; borrows slices from cache on hit.
/// # Mutation
/// `get_interval_bytes` mutates cache (load/evict); not internally synchronized—one instance per thread or external lock.
/// # Biological assumptions
/// Reference bases for variant calling windows; N handling is raw FASTA bytes.
/// # Java equivalence
/// Similar role to GATK `ReferenceDataSource` / htsjdk indexed FASTA access (Rust-native layout).
pub struct ReferenceWindowCache {
    fasta_path: std::path::PathBuf,
    capacity: usize,
    fai: Option<SamtoolsFai>,
    /// Uppercase contig sequences keyed by dictionary-canonical contig name.
    contig_cache: BTreeMap<String, Arc<Vec<u8>>>,
}

impl ReferenceWindowCache {
    pub fn new(fasta_path: impl Into<std::path::PathBuf>, capacity: usize) -> Self {
        let fasta_path = fasta_path.into();
        let fai = SamtoolsFai::load_beside_fasta(&fasta_path);
        Self {
            fasta_path,
            capacity: capacity.max(1),
            fai,
            contig_cache: BTreeMap::new(),
        }
    }

    pub fn get_interval_bytes(
        &mut self,
        dictionary: &SequenceDictionary,
        contig_query: &str,
        start_1based: u64,
        end_1based_inclusive: u64,
    ) -> GatkResult<&[u8]> {
        if start_1based == 0 {
            return Err(GatkError::argument("Genomic position must be >= 1"));
        }
        if start_1based > end_1based_inclusive {
            return Err(GatkError::argument("Interval start must be <= end"));
        }
        let canon = dictionary
            .contig(contig_query)
            .map(|c| c.name.clone())
            .ok_or_else(|| GatkError::argument(format!("Unknown contig: {contig_query}")))?;
        let contig_len = dictionary.contig(&canon).map(|c| c.length).unwrap_or(0);
        // Contigs larger than this are fetched window-wise via `.fai` instead of retaining
        // the whole chromosome in the process-wide Arc map (chr1 alone is ~250 MiB).
        const WHOLE_CONTIG_CACHE_MAX_BP: u64 = 2_000_000;
        if contig_len > WHOLE_CONTIG_CACHE_MAX_BP {
            if let Some(fai) = self.fai.as_ref() {
                if let Some(entry) = fai.entry_for_query(&canon) {
                    let key = format!("{canon}:{start_1based}-{end_1based_inclusive}");
                    if !self.contig_cache.contains_key(&key) {
                        while self.contig_cache.len() >= self.capacity {
                            self.contig_cache.pop_first();
                        }
                        let bytes = read_interval_via_samtools_fai(
                            &self.fasta_path,
                            entry,
                            start_1based,
                            end_1based_inclusive,
                        )?;
                        self.contig_cache.insert(key.clone(), Arc::new(bytes));
                    }
                    #[allow(clippy::expect_used)]
                    let seq = self.contig_cache.get(&key).expect("window cache");
                    return Ok(seq.as_slice());
                }
            }
            // Missing/unusable .fai: whole-contig cache can be multi-GiB (Peak-RSS).
            eprintln!(
                "gatk-core: ReferenceWindowCache loading whole contig {canon} ({contig_len} bp) \
                 without usable .fai — Peak-RSS risk; provide a samtools .fai beside the FASTA"
            );
        }
        if !self.contig_cache.contains_key(&canon) {
            while self.contig_cache.len() >= self.capacity {
                self.contig_cache.pop_first();
            }
            let bytes = shared_contig_bases(&self.fasta_path, &canon, self.fai.as_ref())?;
            // CLONE: needed because owned HashMap/BTree/HashSet key or value.
            self.contig_cache.insert(canon.clone(), bytes);
        }
        // INVARIANT: `canon` was just inserted (or already present) in `contig_cache` above.
        #[allow(clippy::expect_used)]
        let seq = self.contig_cache.get(&canon).expect("contig cache");
        let s = (start_1based - 1) as usize;
        let e = (end_1based_inclusive - 1) as usize;
        seq.get(s..=e).ok_or_else(|| {
            GatkError::argument("Interval extends past end of contig sequence in FASTA")
        })
    }

    /// Number of locally retained contigs (historical name: used to mean exact windows).
    pub fn cached_windows(&self) -> usize {
        self.contig_cache.len()
    }
}

/// Parse `-L` / `--intervals` the same way as `gatk-cli` HaplotypeCaller (file, `;`-separated list, or single token).
pub fn parse_intervals_cli_string(
    dictionary: &SequenceDictionary,
    intervals_val: &str,
) -> GatkResult<Vec<IntervalSpec>> {
    let parsed = if Path::new(intervals_val).exists() {
        IntervalSpec::parse_list_file_with_dictionary(intervals_val, dictionary)?
    } else if intervals_val.contains(';') {
        let (includes, excludes) = IntervalSpec::parse_include_exclude_list(intervals_val)?;
        for i in &includes {
            dictionary.validate_interval(i)?;
        }
        for e in &excludes {
            dictionary.validate_interval(e)?;
        }
        if excludes.is_empty() {
            includes
        } else {
            resolve_interval_specs_includes_excludes(dictionary, &includes, &excludes)?
        }
    } else {
        let t = intervals_val.trim();
        if t.starts_with('^') {
            return Err(GatkError::argument(
                "A leading exclusion (^) is not valid as the only -L token; provide at least one include interval (e.g. 'chr1:1-10;^chr1:5-5').",
            ));
        }
        let one = IntervalSpec::parse(intervals_val)?;
        dictionary.validate_interval(&one)?;
        vec![one]
    };
    Ok(parsed)
}

/// Intervals for traversal: one closed range per contig over the full FASTA when `maybe_intervals` is [`None`].
pub fn intervals_for_haplotype_caller(
    dictionary: &SequenceDictionary,
    maybe_intervals: Option<&str>,
) -> GatkResult<Vec<IntervalSpec>> {
    match maybe_intervals {
        None => Ok(dictionary.whole_genome_interval_specs()),
        Some(s) => parse_intervals_cli_string(dictionary, s),
    }
}

fn contig_name_matches_fasta_entry(fasta_contig: &str, query: &str) -> bool {
    if fasta_contig == query {
        return true;
    }
    let fa = fasta_contig.strip_prefix("chr").unwrap_or(fasta_contig);
    let q = query.strip_prefix("chr").unwrap_or(query);
    fa == q
}

/// Return the reference base (uppercase ASCII) at **1-based inclusive** `pos` on `contig`.
/// Phase 1 / Step 20: small FASTA inputs (linear scan per call). Pair with [`ReferenceWindowCache`] for repeated windows.
pub fn reference_base_at_1based<P: AsRef<Path>>(
    fasta_path: P,
    contig: &str,
    pos_1based: u64,
) -> GatkResult<u8> {
    let v = read_fasta_interval_bytes(fasta_path, contig, pos_1based, pos_1based)?;
    v.first()
        .copied()
        .ok_or_else(|| GatkError::argument("Empty interval"))
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]
mod tests {
    use super::*;

    #[test]
    fn parses_interval_with_range() {
        let iv = IntervalSpec::parse("chr1:1,000-2,000").unwrap();
        assert_eq!(iv.contig, "chr1");
        assert_eq!(iv.start, Some(1000));
        assert_eq!(iv.end, Some(2000));
    }

    #[test]
    fn parses_whole_contig_interval() {
        let iv = IntervalSpec::parse("chr2").unwrap();
        assert_eq!(iv.contig, "chr2");
        assert_eq!(iv.start, None);
        assert_eq!(iv.end, None);
    }

    #[test]
    fn validates_alias_contig_names() {
        let mut dict = SequenceDictionary::new();
        dict.add_contig("1".to_string(), 1000);

        let iv = IntervalSpec::parse("chr1:10-20").unwrap();
        assert!(dict.validate_interval(&iv).is_ok());
    }

    #[test]
    fn parses_interval_list() {
        let intervals = IntervalSpec::parse_list("chr1:1-10; chr2:20-30").unwrap();
        assert_eq!(intervals.len(), 2);
        assert_eq!(intervals[0].contig, "chr1");
        assert_eq!(intervals[1].contig, "chr2");
    }

    #[test]
    fn whole_genome_interval_specs_one_contig() {
        let mut d = SequenceDictionary::new();
        d.add_contig("chr1".to_string(), 32);
        let w = d.whole_genome_interval_specs();
        assert_eq!(w.len(), 1);
        assert_eq!(w[0].contig, "chr1");
        assert_eq!(w[0].start, Some(1));
        assert_eq!(w[0].end, Some(32));
    }

    #[test]
    fn intervals_for_hc_none_is_whole_genome() {
        let mut d = SequenceDictionary::new();
        d.add_contig("chr1".to_string(), 100);
        let v = intervals_for_haplotype_caller(&d, None).unwrap();
        assert_eq!(v.len(), 1);
        let (c, s, e) = v[0].resolve_closed_ends(&d).unwrap();
        assert_eq!(c, "chr1");
        assert_eq!(s, 1);
        assert_eq!(e, 100);
    }

    #[test]
    fn reference_base_at_1based_chr_alias_and_uppercase() {
        let fa = Path::new(env!("CARGO_MANIFEST_DIR")).join("src/tests/test_data/reference.fa");
        assert_eq!(reference_base_at_1based(&fa, "chr1", 1).unwrap(), b'A');
        assert_eq!(reference_base_at_1based(&fa, "1", 4).unwrap(), b'G');
        assert_eq!(reference_base_at_1based(&fa, "chr2", 1).unwrap(), b'G');
    }

    #[test]
    fn reference_base_at_1based_rejects_zero() {
        let fa = Path::new(env!("CARGO_MANIFEST_DIR")).join("src/tests/test_data/reference.fa");
        assert!(reference_base_at_1based(&fa, "chr1", 0).is_err());
    }

    #[test]
    fn parse_list_rejects_exclusion_without_dictionary() {
        assert!(IntervalSpec::parse_list("chr1:1-10;^chr1:5-5").is_err());
    }

    #[test]
    fn resolve_excludes_mid_segment() {
        let mut d = SequenceDictionary::new();
        d.add_contig("chr1".to_string(), 100);
        let inc = vec![IntervalSpec::parse("chr1:1-20").unwrap()];
        let exc = vec![IntervalSpec::parse("chr1:10-15").unwrap()];
        let v = resolve_interval_specs_includes_excludes(&d, &inc, &exc).unwrap();
        assert_eq!(v.len(), 2);
        assert_eq!(v[0].start, Some(1));
        assert_eq!(v[0].end, Some(9));
        assert_eq!(v[1].start, Some(16));
        assert_eq!(v[1].end, Some(20));
    }

    #[test]
    fn parse_intervals_cli_string_union_excludes() {
        let mut d = SequenceDictionary::new();
        d.add_contig("chr1".to_string(), 32);
        let v = parse_intervals_cli_string(&d, "chr1:1-20;^chr1:10-15").unwrap();
        assert_eq!(v.len(), 2);
    }

    #[test]
    fn count_acgtn_histogram_chr1_16() {
        let fa = Path::new(env!("CARGO_MANIFEST_DIR")).join("src/tests/test_data/reference.fa");
        let dict = SequenceDictionary::from_fasta_path(&fa).unwrap();
        let specs = vec![IntervalSpec::parse("chr1:1-16").unwrap()];
        let c = count_acgtn_histogram_for_interval_specs(&fa, &dict, &specs).unwrap();
        assert_eq!(c[0], 4);
        assert_eq!(c[1], 4);
        assert_eq!(c[2], 4);
        assert_eq!(c[3], 4);
        assert_eq!(c[4], 0);
    }

    #[test]
    fn count_histogram_buckets_non_acgtn_as_n() {
        let td = tempfile::tempdir().unwrap();
        let p = td.path().join("mask.fa");
        std::fs::write(&p, ">chr1\nACGTN\n").unwrap();
        let dict = SequenceDictionary::from_fasta_path(&p).unwrap();
        let specs = vec![IntervalSpec::parse("chr1:1-5").unwrap()];
        let c = count_acgtn_histogram_for_interval_specs(&p, &dict, &specs).unwrap();
        assert_eq!(c[4], 1);
    }

    #[test]
    fn reference_window_cache_evicts_deterministically() {
        let fa = Path::new(env!("CARGO_MANIFEST_DIR")).join("src/tests/test_data/reference.fa");
        let dict = SequenceDictionary::from_fasta_path(&fa).unwrap();
        let mut cache = ReferenceWindowCache::new(&fa, 1);
        let _ = cache.get_interval_bytes(&dict, "chr1", 1, 2).unwrap();
        assert_eq!(cache.cached_windows(), 1);
        // Same contig: no extra cache entry (contig-level cache on small fixtures).
        let _ = cache.get_interval_bytes(&dict, "chr1", 5, 6).unwrap();
        assert_eq!(cache.cached_windows(), 1);
    }

    #[test]
    fn dictionary_prefers_fai_without_reading_fasta_body() {
        let td = tempfile::tempdir().unwrap();
        let fa = td.path().join("huge.fa");
        // Header only — a full FASTA body parse cannot recover LN=63025520.
        std::fs::write(&fa, ">chr20\n").unwrap();
        std::fs::write(
            td.path().join("huge.fa.fai"),
            "chr20\t63025520\t8\t60\t61\n",
        )
        .unwrap();
        let dict = SequenceDictionary::from_fasta_path(&fa).expect("dict from fai");
        assert_eq!(dict.contig_count(), 1);
        assert_eq!(dict.contig("chr20").unwrap().length, 63_025_520);
    }

    #[test]
    fn samtools_fai_interval_matches_sequential() {
        let td = tempfile::tempdir().unwrap();
        let fa = td.path().join("tiny.fa");
        // 10 bases, 5 per line → linewidth 6 (5 + newline).
        std::fs::write(&fa, ">chr1\nACGTA\nCGTAN\n").unwrap();
        let fai = td.path().join("tiny.fa.fai");
        // offset: ">chr1\n" = 6 bytes.
        std::fs::write(&fai, "chr1\t10\t6\t5\t6\n").unwrap();
        let via_fai = read_fasta_interval_bytes(&fa, "chr1", 3, 8).unwrap();
        assert_eq!(via_fai, b"GTACGT");
    }
}
