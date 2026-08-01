//! GATK4 `AssemblyRegionIterator`–style region stream. See `docs/ARCHITECTURE.md`.
//! Consumes a [`ReadShard`](`crate::walker::ReadShard`) padded span, walks 1-based loci with the same
//! pileup → activity → band-pass path as the smoothed-activity parity exporters (feature
//! `dev-dumps`), and drains [`BandPassActivityProfile::try_pop_next_ready_region`] like Java’s
//! iterator (non-forcing pops after each locus, forcing flush at span boundaries).
//! Shard load uses [`load_records_for_shard_raw`] + [`crate::read_transformer::apply_shard_read_pipeline`]
//! in [`crate::walker_traversal::drain_assembly_regions_for_shard`]. [`load_all_records_for_contig`] still
//! filters in one step for parity exporters that do not use the B.4 pipeline.
//! **B.3:** after draining regions, use [`crate::walker_apply::WalkerApplyStats`] (or [`crate::HaplotypeCallerEngine::walker_apply_stats`]) to mirror Java `apply` count and inactive `callRegion` branching.

use crate::activity_profile::{
    ActivityProfileRegion, BandPassActivityProfile, BandPassActivityProfileParams,
};
use crate::activity_scoring::HaplotypeCallerActivityScoringParams;
use crate::assembly_region_evaluator::add_locus_for_smoothed_activity;
use crate::feature_context::{FeatureContext, FeatureDataSources};
use crate::genome_loc::GenomePosition;
use crate::locus_iterator::LocusPileupState;
use crate::read_binding::{
    locus_is_strictly_after_closed_interval_1based, record_overlaps_closed_interval_1based,
};
#[cfg(any(feature = "dev-dumps", test))]
use crate::read_binding::{record_alignment_end_1based, record_alignment_start_1based};
use crate::read_header_semantics::ReadHeaderSemantics;
use crate::read_model::{passes_hc_read_filters_with_header, ReadFilterParams};
use crate::reference_context::ReferenceContext;
use crate::region_pileup::RegionPileupLocus;
use crate::walker::ReadShard;
use gatk_common::{AssemblyRegionConfig, GatkError, GatkResult};
use gatk_core::reference::{ReferenceWindowCache, SequenceDictionary};
use rust_htslib::bam;
use rust_htslib::bam::Read as _;
use std::collections::VecDeque;
#[cfg(any(feature = "dev-dumps", test))]
use std::io::Write;
use std::path::Path;

/// Operational safety ceiling on reads retained in one assembly region after positional DS.
///
/// GATK 4.4 caps **per alignment start** (default 50), not total depth. Staggered starts
/// (e.g. centromere / P12 full-30×) can leave hundreds of thousands of reads per region;
/// PairHMM / cloning then OOMs a 16 GiB host. This refuse matches PairHMM oversized-DP
/// fail-closed behavior — it is **not** a genotype-contract downsampler.
pub const MAX_READS_PER_ASSEMBLY_REGION: usize = 100_000;

/// Fail closed when an assembly region retains more reads than
/// [`MAX_READS_PER_ASSEMBLY_REGION`] (after GATK positional DS).
pub fn refuse_oversized_assembly_region_reads(
    contig: &str,
    start_1based: u64,
    end_1based: u64,
    read_count: usize,
) -> GatkResult<()> {
    if read_count <= MAX_READS_PER_ASSEMBLY_REGION {
        return Ok(());
    }
    Err(GatkError::algorithm(format!(
        "assembly region refused oversized read set ({contig}:{start_1based}-{end_1based}, \
         reads={read_count}, max={MAX_READS_PER_ASSEMBLY_REGION}); \
         GATK positional DS (max-reads-per-alignment-start) does not bound total depth when \
         starts are staggered — exclude ultra-deep intervals or use a smaller -L"
    )))
}

/// One assembly region handed to downstream `apply` / `callRegion` (GATK `AssemblyRegion` core fields).
/// # Invariants
/// Unpadded and extended spans are **1-based inclusive**; extended = start/end ± extension clipped to contig.
/// `reads` overlap the extended span; `read_qnames` are sorted unique names derived from those reads.
/// # Ownership
/// Owns reads, reference context, features, and optional pileup loci for the region.
/// # Mutation
/// Built by the iterator; `call_region` / finalize may transform reads but treat span fields as fixed.
/// # Biological assumptions
/// Active regions are locally reassembled; inactive may emit reference-confidence output.
/// # Java equivalence
/// GATK `AssemblyRegion` core fields for HC walker apply.
#[derive(Debug, Clone)]
pub struct AssemblyRegion {
    pub contig: String,
    /// Unpadded region span (1-based inclusive).
    pub start: GenomePosition,
    pub end: GenomePosition,
    pub is_active: bool,
    /// Reference span including extension (`assemblyRegionExtension`), clipped to contig.
    pub extended_start: GenomePosition,
    pub extended_end: GenomePosition,
    pub extension: u32,
    /// Full alignments overlapping the extended span (GAP-B-01 / B.5.1).
    pub reads: Vec<bam::Record>,
    /// Sorted unique qnames derived from [`Self::reads`] (parity / legacy dumps).
    pub read_qnames: Vec<String>,
    /// Padded-span reference bases for `apply` (B.5.3).
    pub reference: ReferenceContext,
    /// Optional resource features on padded span (B.5.4).
    pub features: FeatureContext,
    /// Per-locus pileups on active span when tracking is enabled (B.5.8).
    pub pileup_loci: Vec<RegionPileupLocus>,
}

impl PartialEq for AssemblyRegion {
    fn eq(&self, other: &Self) -> bool {
        self.contig == other.contig
            && self.start == other.start
            && self.end == other.end
            && self.is_active == other.is_active
            && self.extended_start == other.extended_start
            && self.extended_end == other.extended_end
            && self.extension == other.extension
            && self.read_qnames == other.read_qnames
    }
}

impl From<ActivityProfileRegion> for AssemblyRegion {
    fn from(r: ActivityProfileRegion) -> Self {
        Self {
            contig: r.contig,
            start: GenomePosition::new_1based(r.start),
            end: GenomePosition::new_1based(r.end),
            is_active: r.is_active,
            extended_start: GenomePosition::new_1based(r.padded_start),
            extended_end: GenomePosition::new_1based(r.padded_end),
            extension: r.extension,
            reads: Vec::new(),
            read_qnames: Vec::new(),
            reference: ReferenceContext::empty(),
            features: FeatureContext::empty(),
            pileup_loci: Vec::new(),
        }
    }
}

/// Rebuild [`AssemblyRegion::read_qnames`] from [`AssemblyRegion::reads`].
pub fn sync_read_qnames(region: &mut AssemblyRegion) {
    let mut names: Vec<String> = region
        .reads
        .iter()
        .map(|r| String::from_utf8_lossy(r.qname()).into_owned())
        .collect();
    names.sort();
    names.dedup();
    region.read_qnames = names;
}

fn sort_record_indices(records: &[bam::Record]) -> Vec<usize> {
    let mut indices: Vec<usize> = (0..records.len()).collect();
    indices.sort_by(|&a, &b| {
        records[a]
            .pos()
            .cmp(&records[b].pos())
            .then_with(|| records[a].qname().cmp(records[b].qname()))
    });
    indices
}

/// Write region rows + per-read lines for (`b5-reads` gate).
#[cfg(any(feature = "dev-dumps", test))]
pub fn dump_assembly_region_reads_tsv(
    regions: &[AssemblyRegion],
    out: &mut impl Write,
    header: &bam::HeaderView,
    filters: &ReadFilterParams,
) -> GatkResult<()> {
    writeln!(out, "contig\tstart\tend\tis_active\tread_count")
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    for r in regions {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}",
            r.contig,
            r.start.get(),
            r.end.get(),
            r.is_active,
            r.reads.len()
        )
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
        let mut read_rows: Vec<(String, u64, u64)> = r
            .reads
            .iter()
            .filter_map(|rec| {
                let qname = String::from_utf8_lossy(rec.qname()).into_owned();
                let s = record_alignment_start_1based(rec, header, filters)?;
                let e = record_alignment_end_1based(rec, header, filters)?;
                Some((qname, s, e))
            })
            .collect();
        read_rows.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        for (qname, s, e) in read_rows {
            writeln!(out, "read\t{qname}\t{s}\t{e}")
                .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
        }
    }
    Ok(())
}

/// Gate: padded-span [`ReferenceContext`] per region.
#[cfg(any(feature = "dev-dumps", test))]
pub fn dump_assembly_region_reference_tsv(
    regions: &[AssemblyRegion],
    out: &mut impl Write,
) -> GatkResult<()> {
    writeln!(
        out,
        "contig\tstart\tend\tis_active\textended_start\textended_end\tref_window_start\tref_window_end\tref_len"
    )
    .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    for r in regions {
        let refc = &r.reference;
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.contig,
            r.start.get(),
            r.end.get(),
            r.is_active,
            r.extended_start.get(),
            r.extended_end.get(),
            refc.window_start,
            refc.window_end,
            refc.len()
        )
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
        writeln!(out, "ref_bases\t{}", refc.bases_ascii())
            .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    }
    Ok(())
}

/// Gate: [`FeatureContext`] per region (padded span).
#[cfg(any(feature = "dev-dumps", test))]
pub fn dump_assembly_region_features_tsv(
    regions: &[AssemblyRegion],
    out: &mut impl Write,
) -> GatkResult<()> {
    writeln!(
        out,
        "contig\tstart\tend\tis_active\textended_start\textended_end\tfeat_has_backing\tfeat_count"
    )
    .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    for r in regions {
        let fc = &r.features;
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.contig,
            r.start.get(),
            r.end.get(),
            r.is_active,
            r.extended_start.get(),
            r.extended_end.get(),
            fc.has_backing_data_source,
            fc.features.len()
        )
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
        for f in &fc.features {
            let alts = f.alternates.join(",");
            writeln!(
                out,
                "feature\t{}\t{}\t{}\t{}\t{}\t{alts}",
                f.source, f.start, f.end, f.contig, f.reference
            )
            .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
        }
    }
    Ok(())
}

/// Gate: per-region attached pileup loci (`AlignmentContext` spans).
#[cfg(any(feature = "dev-dumps", test))]
pub fn dump_assembly_region_pileup_track_tsv(
    regions: &[AssemblyRegion],
    track_enabled: bool,
    out: &mut impl Write,
) -> GatkResult<()> {
    writeln!(out, "track_enabled\t{}", track_enabled)
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    writeln!(out, "contig\tstart\tend\tis_active\tpileup_count")
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    for r in regions {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}",
            r.contig,
            r.start.get(),
            r.end.get(),
            r.is_active,
            r.pileup_loci.len()
        )
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
        for p in &r.pileup_loci {
            writeln!(out, "pileup\t{}\t{}\t{}", p.contig, p.pos, p.depth)
                .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
        }
    }
    Ok(())
}

/// Arguments mirroring GATK `AssemblyRegionArgumentCollection` + band-pass wiring passed into `AssemblyRegionIterator`.
/// # Invariants
/// `min_region_size` ≤ `max_region_size`; extension pads active/inactive cuts.
/// HC defaults: padding 100, min/max region 50/300, adaptive band-pass.
/// # Ownership
/// Owns nested band-pass/scoring configs and optional feature sources.
/// # Mutation
/// Snapshot for iterator construction; not mutated during drain.
/// # Biological assumptions
/// Controls how activity profiles are cut into assembly regions for HC.
/// # Java equivalence
/// GATK `AssemblyRegionArgumentCollection` + iterator band-pass wiring.
#[derive(Debug, Clone)]
pub struct AssemblyRegionIteratorConfig {
    pub assembly_region_extension: u32,
    pub min_region_size: u32,
    pub max_region_size: u32,
    pub band_pass: BandPassActivityProfileParams,
    pub scoring: HaplotypeCallerActivityScoringParams,
    /// GATK `assemblyRegionArgs.forceActive` — mark every emitted region active.
    pub force_active: bool,
    /// Optional VCF / resource features (B.5.4); `None` = HC default (no `-dbsnp` etc.).
    pub feature_sources: Option<FeatureDataSources>,
    /// GATK `shouldTrackPileupsForAssemblyRegions` (HC: `usePileupDetection`, default false).
    pub track_pileups: bool,
}

impl AssemblyRegionIteratorConfig {
    /// Defaults from [`AssemblyRegionConfig::default`] (padding 100, min/max assembly region 50/300, HC band-pass).
    /// Region bounds are statically valid (`min ≤ max`); no panic path.
    pub fn gatk_haplotype_caller_defaults() -> Self {
        let ar = AssemblyRegionConfig::default();
        debug_assert!(ar.min_assembly_region_size <= ar.max_assembly_region_size);
        Self {
            assembly_region_extension: ar.assembly_region_padding,
            min_region_size: ar.min_assembly_region_size,
            max_region_size: ar.max_assembly_region_size,
            band_pass: BandPassActivityProfileParams::gatk_haplotype_caller_defaults(),
            scoring: HaplotypeCallerActivityScoringParams::default(),
            force_active: false,
            feature_sources: None,
            track_pileups: false,
        }
    }

    /// Reject inverted region size bounds (`min > max`) for user-built configs.
    pub fn try_validate(self) -> Result<Self, gatk_common::GatkError> {
        if self.min_region_size > self.max_region_size {
            return Err(gatk_common::GatkError::invalid_argument(
                "min_region_size",
                format!(
                    "min_region_size ({}) must be ≤ max_region_size ({})",
                    self.min_region_size, self.max_region_size
                ),
            ));
        }
        Ok(self)
    }
}

impl Default for AssemblyRegionIteratorConfig {
    fn default() -> Self {
        Self::gatk_haplotype_caller_defaults()
    }
}

/// Linear BAM/SAM scan: all records on `contig` without HC filtering (B.4 shard pipeline).
pub fn load_all_records_for_contig_raw(
    bam_path: &Path,
    contig: &str,
) -> GatkResult<(bam::HeaderView, Vec<bam::Record>)> {
    let mut reader = bam::Reader::from_path(bam_path)
        .map_err(|e| GatkError::io_message(format!("open BAM {}: {e}", bam_path.display())))?;
    let header = reader.header().clone();
    let tid = header.tid(contig.as_bytes()).ok_or_else(|| {
        GatkError::invalid_argument("contig", format!("BAM header missing contig {contig}"))
    })? as i32;
    let mut out = Vec::new();
    for res in reader.records() {
        let rec = res
            .map_err(|e| GatkError::io_message(format!("read BAM {}: {e}", bam_path.display())))?;
        if rec.tid() == tid {
            out.push(rec);
        }
    }
    Ok((header, out))
}

/// Linear BAM/SAM scan: all primary-filtered records on `contig` (for small fixtures).
pub fn load_all_records_for_contig(
    bam_path: &Path,
    contig: &str,
    filters: &ReadFilterParams,
) -> GatkResult<(bam::HeaderView, Vec<bam::Record>)> {
    let (header, mut out) = load_all_records_for_contig_raw(bam_path, contig)?;
    out.retain(|rec| passes_hc_read_filters_with_header(rec, &header, filters));
    Ok((header, out))
}

/// Records for a shard via indexed interval query (Java `SamReaderQueryingIterator` / `queryOverlapping`).
pub fn load_records_for_shard_raw(
    bam_path: &Path,
    shard: &ReadShard,
) -> GatkResult<(bam::HeaderView, Vec<bam::Record>)> {
    let mut reader = bam::IndexedReader::from_path(bam_path).map_err(|e| {
        GatkError::io_message(format!("open indexed BAM {}: {e}", bam_path.display()))
    })?;
    let header = reader.header().clone();
    let tid = header.tid(shard.contig.as_bytes()).ok_or_else(|| {
        GatkError::invalid_argument(
            "contig",
            format!("BAM header missing contig {}", shard.contig),
        )
    })? as i32;
    let mut out = Vec::new();
    for &(start1, end1) in &shard.padded_spans {
        let start0 = start1.saturating_sub(1);
        reader
            .fetch((tid, start0 as i64, end1 as i64))
            .map_err(|e| {
                GatkError::io_message(format!(
                    "fetch {}:{}-{} from {}: {e}",
                    shard.contig,
                    start1,
                    end1,
                    bam_path.display()
                ))
            })?;
        for res in reader.records() {
            let rec = res.map_err(|e| {
                GatkError::io_message(format!("read BAM {}: {e}", bam_path.display()))
            })?;
            out.push(rec);
        }
    }
    dedupe_records_by_alignment_start(&mut out);
    Ok((header, out))
}

/// One alignment per `(qname, pos, flags)` when padded shard spans overlap.
/// Preserves first-seen order; keeps both mates (same POS, different `flags`).
fn dedupe_records_by_alignment_start(records: &mut Vec<bam::Record>) {
    let mut seen = std::collections::HashSet::new();
    records.retain(|rec| seen.insert((rec.qname().to_vec(), rec.pos(), rec.flags())));
}

/// Iterator over [`AssemblyRegion`] values for one [`ReadShard`].
/// # Invariants
/// Walks loci in increasing 1-based order; non-forcing pops after each locus, force flush at span end.
/// Emitted regions respect min/max size and extension from config.
/// # Ownership
/// Owns contig state, reference bytes, BAM records, header, and pending region queue.
/// # Mutation
/// Advances locus/span state in place; popping drains ready regions from the band-pass profile.
/// # Biological assumptions
/// Converts per-locus activity into contiguous active/inactive assembly regions.
/// # Java equivalence
/// GATK `AssemblyRegionIterator` drain semantics.
pub struct AssemblyRegionIterator {
    contig: String,
    contig_len: u64,
    /// Locus walk spans (user `-L`, merged).
    spans: Vec<(u64, u64)>,
    /// Read overlap bounds (padded).
    _read_spans: Vec<(u64, u64)>,
    span_idx: usize,
    next_pos: Option<u64>,
    span_start: u64,
    span_end: u64,
    ref_bytes: Vec<u8>,
    all_records: Vec<bam::Record>,
    header: bam::HeaderView,
    header_semantics: ReadHeaderSemantics,
    profile: BandPassActivityProfile,
    ref_cache: ReferenceWindowCache,
    dict: SequenceDictionary,
    read_filters: ReadFilterParams,
    cfg: AssemblyRegionIteratorConfig,
    /// GATK `pendingRegions` — popped only after locus iterator passes `paddedSpan.end`.
    pending: VecDeque<AssemblyRegion>,
    finished: bool,
    pileup_state: LocusPileupState,
    /// Coordinate-sorted indices into [`Self::all_records`] (GATK read-shard order).
    sorted_record_indices: Vec<usize>,
    /// GATK `readCache` — indices of reads not yet assigned to a region.
    read_cache: VecDeque<usize>,
    ingest_cursor: usize,
    /// GATK `previousRegionReads` as indices into [`Self::all_records`] (R4-1: no BAM clones).
    previous_region_indices: Vec<usize>,
    /// Cached half-open reference spans `(r0, r1)` per `all_records` index (`r1 < 0` = filtered out).
    read_ref_span0: Vec<(i64, i64)>,
    /// GATK `pendingAlignmentData` when [`AssemblyRegionIteratorConfig::track_pileups`].
    pending_pileups: Option<VecDeque<RegionPileupLocus>>,
}

impl AssemblyRegionIterator {
    /// Build an iterator for `shard` using pre-loaded alignments on that contig.
    pub fn try_new(
        shard: &ReadShard,
        dictionary: &SequenceDictionary,
        reference_fasta: &Path,
        all_records: Vec<bam::Record>,
        header: bam::HeaderView,
        read_filters: ReadFilterParams,
        cfg: AssemblyRegionIteratorConfig,
    ) -> GatkResult<Self> {
        let contig = shard.contig.clone();
        let contig_len = dictionary
            .contig(&contig)
            .ok_or_else(|| GatkError::argument(format!("unknown contig {contig} in dictionary")))?
            .length;
        // GATK walks user intervals for loci (`getIntervals`), not the padded read bounds.
        let locus_spans = if shard.user_spans.is_empty() {
            shard.padded_spans.clone()
        } else {
            shard.user_spans.clone()
        };
        let read_spans = shard.padded_spans.clone();
        let finished = locus_spans.is_empty();
        let mut filtered_records = all_records;
        filtered_records.retain(|r| {
            read_spans.iter().any(|&(rs, re)| {
                record_overlaps_closed_interval_1based(r, &header, &contig, rs, re, &read_filters)
            })
        });
        let sorted_record_indices = sort_record_indices(&filtered_records);
        // Spans for already-filtered records (geometry only; filter already applied above).
        let mut read_ref_span0 = vec![(-1i64, -1i64); filtered_records.len()];
        for (i, rec) in filtered_records.iter().enumerate() {
            let r0 = rec.pos();
            let r1 = rec.cigar().end_pos();
            if r1 > r0 {
                read_ref_span0[i] = (r0, r1);
            }
        }
        let pileup_state =
            LocusPileupState::from_records(&filtered_records, &header, &contig, &read_filters);
        let header_semantics = ReadHeaderSemantics::from_bam_header_view(&header)?;
        let track_pileups = cfg.track_pileups;
        let mut it = Self {
            // CLONE: needed because owned contig id for output record.
            contig: contig.clone(),
            contig_len,
            spans: locus_spans,
            _read_spans: read_spans,
            span_idx: 0,
            next_pos: None,
            span_start: 0,
            span_end: 0,
            ref_bytes: Vec::new(),
            all_records: filtered_records,
            header,
            header_semantics,
            profile: BandPassActivityProfile::new(
                contig.as_str(),
                contig_len,
                cfg.band_pass.clone(),
            ),
            ref_cache: ReferenceWindowCache::new(reference_fasta.to_path_buf(), 4),
            dict: dictionary.clone(),
            read_filters,
            cfg,
            pending: VecDeque::new(),
            finished,
            pileup_state,
            sorted_record_indices,
            read_cache: VecDeque::new(),
            ingest_cursor: 0,
            previous_region_indices: Vec::new(),
            read_ref_span0,
            pending_pileups: if track_pileups {
                Some(VecDeque::new())
            } else {
                None
            },
        };
        if !finished {
            it.begin_span()?;
        }
        Ok(it)
    }

    fn begin_span(&mut self) -> GatkResult<()> {
        let (s, e) = self.spans[self.span_idx];
        if self.pileup_state.last_pos1.is_some_and(|prev| s < prev) {
            self.pileup_state.reset_cursor();
        }
        self.span_start = s;
        self.span_end = e;
        self.ref_bytes = self
            .ref_cache
            .get_interval_bytes(&self.dict, &self.contig, s, e)
            .map_err(|e| GatkError::generic(e.to_string()))?
            .to_vec();
        // Note: do not clone shard reads into a span-local vec — `all_records` already
        // holds them; a prior `span_records` field was written and never read (pure Peak-RSS tax).
        self.profile = BandPassActivityProfile::new(
            self.contig.as_str(),
            self.contig_len,
            self.cfg.band_pass.clone(),
        );
        self.next_pos = Some(s);
        Ok(())
    }

    fn flush_profile_to_pending(&mut self, force_conversion: bool) -> GatkResult<()> {
        let ext = self.cfg.assembly_region_extension;
        let min = self.cfg.min_region_size;
        let max = self.cfg.max_region_size;
        for reg in self
            .profile
            .pop_ready_regions(ext, min, max, force_conversion)?
        {
            // Observe-only semantic checkpoint (no algorithm effect).
            crate::semantic_trace::emit_activity_profile_cut(
                &reg.contig,
                reg.start,
                reg.end,
                reg.is_active,
                reg.padded_start,
                reg.padded_end,
                reg.extension,
            );
            self.pending.push_back(AssemblyRegion::from(reg));
        }
        Ok(())
    }

    fn finalize_profile_at_shard_end(&mut self) -> GatkResult<()> {
        if !self.profile.is_empty() {
            self.flush_profile_to_pending(true)?;
        }
        Ok(())
    }

    /// GATK `ReadCachingIterator` + LIBS: ingest alignments with start ≤ current locus into `read_cache`.
    fn ingest_reads_through_locus(&mut self, pos1: u64) {
        let ref_pos0 = pos1.saturating_sub(1) as i64;
        while self.ingest_cursor < self.sorted_record_indices.len() {
            let idx = self.sorted_record_indices[self.ingest_cursor];
            if self.all_records[idx].pos() > ref_pos0 {
                break;
            }
            self.read_cache.push_back(idx);
            self.ingest_cursor += 1;
        }
    }

    #[inline]
    fn idx_overlaps_closed(&self, idx: usize, start1: u64, end1: u64) -> bool {
        let (r0, r1) = self.read_ref_span0[idx];
        if r1 < 0 {
            return false;
        }
        let (i0, i1) = crate::read_binding::closed_interval_1based_to_ref_span0(start1, end1);
        crate::read_binding::half_open_overlaps(r0, r1, i0, i1)
    }

    #[inline]
    fn idx_strictly_after_closed(&self, idx: usize, end1: u64) -> bool {
        let (r0, r1) = self.read_ref_span0[idx];
        if r1 < 0 {
            return false;
        }
        let (_i0, i1) = crate::read_binding::closed_interval_1based_to_ref_span0(1, end1);
        r0 >= i1
    }

    /// GATK `fillNextAssemblyRegionWithReads`.
    fn fill_region_with_reads(&mut self, region: &mut AssemblyRegion) {
        let mut indices = Vec::new();
        for &idx in &self.previous_region_indices {
            if self.idx_overlaps_closed(idx, region.extended_start.get(), region.extended_end.get())
            {
                indices.push(idx);
            }
        }
        while let Some(&idx) = self.read_cache.front() {
            if self.idx_strictly_after_closed(idx, region.extended_end.get()) {
                break;
            }
            self.read_cache.pop_front();
            if self.idx_overlaps_closed(idx, region.extended_start.get(), region.extended_end.get())
            {
                indices.push(idx);
            }
        }
        indices.sort_by(|&a, &b| {
            let ra = &self.all_records[a];
            let rb = &self.all_records[b];
            ra.pos()
                .cmp(&rb.pos())
                .then_with(|| ra.qname().cmp(rb.qname()))
                .then_with(|| ra.flags().cmp(&rb.flags()))
        });
        // GATK keeps both mates when they share POS; dedup only identical alignments.
        indices.dedup_by(|&mut a, &mut b| {
            let ra = &self.all_records[a];
            let rb = &self.all_records[b];
            ra.qname() == rb.qname() && ra.pos() == rb.pos() && ra.flags() == rb.flags()
        });
        // CLONE: needed — `AssemblyRegion` owns `bam::Record`s while `all_records` remains
        // the shard-wide source for overlapping / previous-region reuse. Follow-up: `Arc<Record>`.
        region.reads = indices
            .iter()
            .map(|&i| self.all_records[i].clone())
            .collect();
        sync_read_qnames(region);
        self.previous_region_indices = indices;
    }

    fn fill_region_with_pileup_data(&mut self, region: &mut AssemblyRegion) {
        let Some(pending) = self.pending_pileups.as_mut() else {
            return;
        };
        while let Some(front) = pending.front() {
            if front.contig != region.contig || front.pos < region.start.get() {
                pending.pop_front();
            } else {
                break;
            }
        }
        let mut overlapping = Vec::new();
        while let Some(front) = pending.front() {
            if front.contig != region.contig {
                break;
            }
            if front.pos <= region.end.get() {
                if let Some(locus) = pending.pop_front() {
                    overlapping.push(locus);
                }
            } else {
                break;
            }
        }
        let mut new_pending: VecDeque<RegionPileupLocus> = overlapping.iter().cloned().collect();
        new_pending.extend(pending.drain(..));
        region.pileup_loci = overlapping;
        *pending = new_pending;
    }

    fn finish_pending_region(&mut self, mut region: AssemblyRegion) -> GatkResult<AssemblyRegion> {
        self.fill_region_with_reads(&mut region);
        refuse_oversized_assembly_region_reads(
            &region.contig,
            region.start.get(),
            region.end.get(),
            region.reads.len(),
        )?;
        self.fill_region_with_pileup_data(&mut region);
        if self.cfg.force_active {
            region.is_active = true;
        }
        let mut out = region;
        out.reference = ReferenceContext::from_interval(
            &self.dict,
            &mut self.ref_cache,
            &out.contig,
            out.extended_start.get(),
            out.extended_end.get(),
        )?;
        out.features = FeatureContext::from_padded_span(
            &out.contig,
            out.extended_start.get(),
            out.extended_end.get(),
            self.cfg.feature_sources.as_ref(),
        );
        // `previous_region_indices` already updated in `fill_region_with_reads`.
        crate::semantic_trace::emit_active_region(&out);
        Ok(out)
    }

    /// Pop the front pending region when `pos1` is strictly after its padded span (GAP-B-03).
    fn try_poll_pending_after_locus(&mut self, pos1: u64) -> GatkResult<Option<AssemblyRegion>> {
        if let Some(front) = self.pending.front() {
            if locus_is_strictly_after_closed_interval_1based(pos1, front.extended_end.get()) {
                if let Some(region) = self.pending.pop_front() {
                    return Ok(Some(self.finish_pending_region(region)?));
                }
            }
        }
        Ok(None)
    }

    fn take_pending_without_after_check(&mut self) -> GatkResult<Option<AssemblyRegion>> {
        match self.pending.pop_front() {
            Some(r) => Ok(Some(self.finish_pending_region(r)?)),
            None => Ok(None),
        }
    }

    /// Next region in traversal order, or `None` when the shard is exhausted.
    pub fn next_region(&mut self) -> GatkResult<Option<AssemblyRegion>> {
        if self.finished {
            return Ok(None);
        }

        loop {
            if self.span_idx >= self.spans.len() {
                self.finalize_profile_at_shard_end()?;
                if let Some(region) = self.take_pending_without_after_check()? {
                    return Ok(Some(region));
                }
                self.finished = true;
                return Ok(None);
            }

            let Some(pos) = self.next_pos else {
                if !self.profile.is_empty() {
                    self.flush_profile_to_pending(true)?;
                }
                if let Some(region) = self.try_poll_pending_after_locus(self.span_end)? {
                    return Ok(Some(region));
                }
                self.span_idx += 1;
                if self.span_idx >= self.spans.len() {
                    self.finalize_profile_at_shard_end()?;
                    if let Some(region) = self.take_pending_without_after_check()? {
                        return Ok(Some(region));
                    }
                    self.finished = true;
                    return Ok(None);
                }
                self.begin_span()?;
                continue;
            };

            if pos > self.span_end {
                return Err(GatkError::argument(format!(
                    "assembly iterator pos {pos} > span_end {}",
                    self.span_end
                )));
            }

            let ref_base = *self
                .ref_bytes
                .get((pos - self.span_start) as usize)
                .ok_or_else(|| {
                    GatkError::argument("reference window index out of range for assembly span")
                })?;

            if !self.profile.is_empty() {
                let force_conversion = self
                    .profile
                    .last_input_pos()
                    .is_some_and(|end| pos != end + 1);
                self.flush_profile_to_pending(force_conversion)?;
            }

            add_locus_for_smoothed_activity(
                &mut self.profile,
                &self.all_records,
                &self.header,
                &self.header_semantics,
                &self.contig,
                pos,
                &self.read_filters,
                ref_base,
                &self.cfg.scoring,
                Some(&mut self.pileup_state),
                false,
            )?;
            if let Some(pending) = self.pending_pileups.as_mut() {
                let depth = self.pileup_state.pileup_depth();
                pending.push_back(RegionPileupLocus {
                    // CLONE: needed because owned contig id for output record.
                    contig: self.contig.clone(),
                    pos,
                    depth,
                });
            }
            self.ingest_reads_through_locus(pos);

            if let Some(region) = self.try_poll_pending_after_locus(pos)? {
                self.next_pos = if pos < self.span_end {
                    Some(pos + 1)
                } else {
                    None
                };
                return Ok(Some(region));
            }

            self.next_pos = if pos < self.span_end {
                Some(pos + 1)
            } else {
                None
            };
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::activity_profile::{
        ActivityProfileRegion, ActivityProfileState, BandPassActivityProfile,
    };
    use crate::walker::make_read_shards;

    #[test]
    fn oversized_assembly_region_read_count_is_refused() {
        assert!(
            refuse_oversized_assembly_region_reads("2", 1, 100, MAX_READS_PER_ASSEMBLY_REGION)
                .is_ok()
        );
        let err = refuse_oversized_assembly_region_reads(
            "2",
            92_300_000,
            92_350_000,
            MAX_READS_PER_ASSEMBLY_REGION + 1,
        )
        .expect_err("must refuse above ceiling");
        let msg = err.to_string();
        assert!(
            msg.contains("refused oversized read set"),
            "unexpected message: {msg}"
        );
        assert!(
            msg.contains(&format!("max={MAX_READS_PER_ASSEMBLY_REGION}")),
            "unexpected message: {msg}"
        );
    }

    #[test]
    fn assembly_region_from_activity_profile_region() {
        let r = ActivityProfileRegion {
            contig: "1".into(),
            start: 10,
            end: 20,
            is_active: true,
            padded_start: 1,
            padded_end: 30,
            extension: 100,
        };
        let a = AssemblyRegion::from(r.clone());
        assert_eq!(a.contig, "1");
        assert_eq!((a.start.get(), a.end.get(), a.is_active), (10, 20, true));
        assert_eq!((a.extended_start.get(), a.extended_end.get()), (1, 30));
        assert_eq!(a.extension, 100);
    }

    #[test]
    fn empty_shard_yields_no_regions() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../parity/fixtures");
        let dict = SequenceDictionary::from_fasta_path(root.join("reference.fa")).unwrap();
        let (header, _recs) = load_all_records_for_contig(
            &root.join("sample.bam"),
            "chr1",
            &ReadFilterParams::default(),
        )
        .unwrap();
        let shard = ReadShard {
            contig: "chr1".into(),
            user_spans: vec![],
            padded_spans: vec![],
        };
        let mut it = AssemblyRegionIterator::try_new(
            &shard,
            &dict,
            root.join("reference.fa").as_path(),
            Vec::new(),
            header,
            ReadFilterParams::default(),
            AssemblyRegionIteratorConfig::default(),
        )
        .unwrap();
        assert!(it.next_region().unwrap().is_none());
    }

    #[test]
    fn profile_only_emits_inactive_tail_matches_iterator_contract() {
        let mut p = BandPassActivityProfile::new(
            "chr1",
            32,
            BandPassActivityProfileParams {
                max_prob_propagation_distance: 0,
                active_prob_threshold: 0.5,
                max_filter_size: 0,
                sigma: crate::activity_profile::PositiveSigma::try_new(1.0).unwrap(),
                adaptive_filter_size: false,
            },
        );
        for pos in 1..=5 {
            p.add(ActivityProfileState::new("chr1", pos, 0.0)).unwrap();
        }
        let mut got = Vec::new();
        while let Some(r) = p.try_pop_next_ready_region(0, 1, 100, true).unwrap() {
            got.push(r);
        }
        assert_eq!(got.len(), 1);
        assert!(!got[0].is_active);
        let a = AssemblyRegion::from(got[0].clone());
        assert_eq!((a.start.get(), a.end.get()), (1, 5));
    }

    #[test]
    fn pending_region_waits_until_locus_after_padded_end() {
        let mut pending = VecDeque::new();
        pending.push_back(AssemblyRegion {
            contig: "chr1".into(),
            start: GenomePosition::new_1based(1),
            end: GenomePosition::new_1based(5),
            is_active: false,
            extended_start: GenomePosition::new_1based(1),
            extended_end: GenomePosition::new_1based(10),
            extension: 5,
            reads: Vec::new(),
            read_qnames: Vec::new(),
            reference: ReferenceContext::empty(),
            features: FeatureContext::empty(),
            pileup_loci: Vec::new(),
        });
        let front = pending.front().unwrap();
        assert!(!locus_is_strictly_after_closed_interval_1based(
            5,
            front.extended_end.get()
        ));
        assert!(locus_is_strictly_after_closed_interval_1based(
            11,
            front.extended_end.get()
        ));
    }

    #[test]
    fn load_shard_dedupes_read_across_overlapping_padded_spans() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../parity/fixtures");
        let bam = root.join("sample.bam");
        let ref_fa = root.join("reference.fa");
        let dict = SequenceDictionary::from_fasta_path(&ref_fa).unwrap();
        let specs =
            gatk_core::reference::parse_intervals_cli_string(&dict, "chr1:1-5;chr1:20-25").unwrap();
        let shards = make_read_shards(&dict, &specs, 5).unwrap();
        let (_, recs) = load_records_for_shard_raw(&bam, &shards[0]).unwrap();
        assert_eq!(
            recs.len(),
            1,
            "single read must not be duplicated across padded fetches"
        );
    }

    #[test]
    fn track_pileups_disjoint_padded_spans_depth_matches_java() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../parity/fixtures");
        let ref_fa = root.join("reference.fa");
        let bam = root.join("sample.bam");
        let dict = SequenceDictionary::from_fasta_path(&ref_fa).unwrap();
        let specs =
            gatk_core::reference::parse_intervals_cli_string(&dict, "chr1:1-5;chr1:20-25").unwrap();
        let shards = make_read_shards(&dict, &specs, 5).unwrap();
        let filters = ReadFilterParams::gatk_standard_hc();
        let (header, recs) = load_records_for_shard_raw(&bam, &shards[0]).unwrap();
        let mut cfg = AssemblyRegionIteratorConfig::gatk_haplotype_caller_defaults();
        cfg.track_pileups = true;
        let mut it =
            AssemblyRegionIterator::try_new(&shards[0], &dict, &ref_fa, recs, header, filters, cfg)
                .unwrap();
        let mut regions = Vec::new();
        while let Some(r) = it.next_region().unwrap() {
            regions.push(r);
        }
        let depths: Vec<_> = regions
            .iter()
            .flat_map(|r| r.pileup_loci.iter().map(|p| p.depth))
            .collect();
        assert!(
            depths.iter().all(|&d| d == 1),
            "expected depth 1 at every pileup locus, got {depths:?}"
        );
    }

    #[test]
    fn track_pileups_attaches_loci_on_active_span() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../parity/fixtures");
        let ref_fa = root.join("reference.fa");
        let bam = root.join("sample.bam");
        let dict = SequenceDictionary::from_fasta_path(&ref_fa).unwrap();
        let specs = gatk_core::reference::parse_intervals_cli_string(&dict, "chr1:5-15").unwrap();
        let shards = make_read_shards(&dict, &specs, 100).unwrap();
        let filters = ReadFilterParams::gatk_standard_hc();
        let (header, recs) = load_all_records_for_contig(&bam, "chr1", &filters).unwrap();
        let mut cfg = AssemblyRegionIteratorConfig::gatk_haplotype_caller_defaults();
        cfg.track_pileups = true;
        let mut it =
            AssemblyRegionIterator::try_new(&shards[0], &dict, &ref_fa, recs, header, filters, cfg)
                .unwrap();
        let region = it.next_region().unwrap().expect("one region");
        assert_eq!(region.start.get(), 5);
        assert_eq!(region.end.get(), 15);
        assert_eq!(region.pileup_loci.len(), 11);
        assert!(region.pileup_loci.iter().all(|p| p.pos >= 5 && p.pos <= 15));
    }

    #[test]
    fn regions_attach_features_from_vcf_on_chr1_fixture() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../parity/fixtures");
        let ref_fa = root.join("reference.fa");
        let bam = root.join("sample.bam");
        let vcf = root.join("hc-full-parity/b5-feature/forced_chr1_10.vcf");
        let dict = SequenceDictionary::from_fasta_path(&ref_fa).unwrap();
        let specs = gatk_core::reference::parse_intervals_cli_string(&dict, "chr1:5-15").unwrap();
        let shards = make_read_shards(&dict, &specs, 100).unwrap();
        let filters = ReadFilterParams::gatk_standard_hc();
        let (header, recs) = load_all_records_for_contig(&bam, "chr1", &filters).unwrap();
        let mut sources = FeatureDataSources::default();
        sources.load_vcf_source("alleles", &vcf).unwrap();
        let mut cfg = AssemblyRegionIteratorConfig::gatk_haplotype_caller_defaults();
        cfg.feature_sources = Some(sources);
        let mut it =
            AssemblyRegionIterator::try_new(&shards[0], &dict, &ref_fa, recs, header, filters, cfg)
                .unwrap();
        let region = it.next_region().unwrap().expect("one region");
        assert!(region.features.has_backing_data_source);
        assert_eq!(region.features.features.len(), 1);
        assert_eq!(region.features.features[0].start, 10);
    }

    #[test]
    fn regions_carry_reference_bases_on_chr1_fixture() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../parity/fixtures");
        let ref_fa = root.join("reference.fa");
        let bam = root.join("sample.bam");
        let dict = SequenceDictionary::from_fasta_path(&ref_fa).unwrap();
        let specs = gatk_core::reference::parse_intervals_cli_string(&dict, "chr1:5-15").unwrap();
        let shards = make_read_shards(&dict, &specs, 100).unwrap();
        let filters = ReadFilterParams::gatk_standard_hc();
        let (header, recs) = load_all_records_for_contig(&bam, "chr1", &filters).unwrap();
        let mut it = AssemblyRegionIterator::try_new(
            &shards[0],
            &dict,
            &ref_fa,
            recs,
            header,
            filters,
            AssemblyRegionIteratorConfig::default(),
        )
        .unwrap();
        let region = it.next_region().unwrap().expect("one region");
        assert_eq!(region.reference.window_start, region.extended_start.get());
        assert_eq!(region.reference.window_end, region.extended_end.get());
        assert_eq!(
            region.reference.len(),
            (region.extended_end.get() - region.extended_start.get() + 1) as usize
        );
        assert!(!region.reference.bases.is_empty());
    }

    #[test]
    fn regions_carry_full_reads_on_chr1_fixture() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../parity/fixtures");
        let ref_fa = root.join("reference.fa");
        let bam = root.join("sample.bam");
        let dict = SequenceDictionary::from_fasta_path(&ref_fa).unwrap();
        let specs = gatk_core::reference::parse_intervals_cli_string(&dict, "chr1:5-15").unwrap();
        let shards = make_read_shards(&dict, &specs, 100).unwrap();
        let filters = ReadFilterParams::gatk_standard_hc();
        let (header, recs) = load_all_records_for_contig(&bam, "chr1", &filters).unwrap();
        let mut it = AssemblyRegionIterator::try_new(
            &shards[0],
            &dict,
            &ref_fa,
            recs,
            header,
            filters,
            AssemblyRegionIteratorConfig::default(),
        )
        .unwrap();
        let mut total_reads = 0usize;
        while let Some(r) = it.next_region().unwrap() {
            assert!(!r.reads.is_empty(), "expected reads on extended span");
            assert_eq!(r.read_qnames.len(), r.reads.len(), "unique qnames vs reads");
            total_reads += r.reads.len();
        }
        assert!(total_reads > 0);
    }

    #[test]
    fn iterator_runs_on_fixture_shard_and_bam() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../parity/fixtures");
        let ref_fa = root.join("reference.fa");
        let bam = root.join("sample.bam");
        let dict = SequenceDictionary::from_fasta_path(&ref_fa).unwrap();
        let specs = gatk_core::reference::parse_intervals_cli_string(&dict, "chr1:5-15").unwrap();
        let shards = make_read_shards(&dict, &specs, 100).unwrap();
        assert_eq!(shards.len(), 1);
        let (header, recs) =
            load_all_records_for_contig(&bam, "chr1", &ReadFilterParams::default()).unwrap();
        let mut it = AssemblyRegionIterator::try_new(
            &shards[0],
            &dict,
            &ref_fa,
            recs,
            header,
            ReadFilterParams::default(),
            AssemblyRegionIteratorConfig::default(),
        )
        .unwrap();
        let mut regions = Vec::new();
        while let Some(r) = it.next_region().unwrap() {
            regions.push(r);
        }
        assert!(
            !regions.is_empty(),
            "expected at least one inactive/active block over chr1:1-32"
        );
        for r in &regions {
            assert_eq!(r.contig, "chr1");
            assert!(r.start <= r.end);
            assert!(r.extended_start <= r.start);
            assert!(r.end <= r.extended_end);
        }
    }
}
