// EventMap synchronization from haplotype CIGARs (Sprint I-2b / I-3).
// Included into [`crate::read_event_discovery`] via `include!`.

/// Options for [`sync_assembly_events_from_haplotype_cigars_with_harvest`] (Sprint I-3).
/// # Invariants
/// [`Self::strict_java`] sets `strict_event_map_only == true` and disables trim SNP harvest.
/// Harvest and strict modes are mutually exclusive for production `CallRegionArgs::strict_java`.
/// # Ownership
/// [`Copy`] config passed into sync helpers; assembly result sets are mutably borrowed.
/// # Mutation
/// Immutable options tag; sync functions mutate `AssemblyResultSet` variation events.
/// # Biological assumptions
/// Post-trim SNP harvest is a parity supplement, not GATK production EventMap behavior.
/// # Java equivalence
/// Mirrors GATK `AssemblyResultSet.regenerateVariationEvents` / `EventMap.buildEventMapsForHaplotypes` when strict.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SyncAssemblyOptions {
    /// Post-`trim_to` SNP harvest from alt haplotypes (off under production `strict_java`).
    pub harvest_trim_snps: bool,
    /// When true, events come only from EventMap / CIGAR regen (production `strict_java`).
    pub strict_event_map_only: bool,
}

impl SyncAssemblyOptions {
    /// Production `CallRegionArgs::strict_java` — EventMap only, no trim SNP harvest.
    pub const fn strict_java() -> Self {
        Self {
            harvest_trim_snps: false,
            strict_event_map_only: true,
        }
    }

    /// Non-strict / parity: optional harvest + side-list merge behavior.
    pub const fn with_harvest(harvest_trim_snps: bool) -> Self {
        Self {
            harvest_trim_snps,
            strict_event_map_only: false,
        }
    }

    pub const fn from_strict_java(is_strict_java: bool) -> Self {
        if is_strict_java {
            Self::strict_java()
        } else {
            Self::with_harvest(true)
        }
    }
}

/// ASM-8 + Java `EventMap`: refresh alt CIGARs on trim slice; events from full padded ref (GATK parity).
pub fn sync_assembly_events_from_haplotype_cigars(
    assembly: &mut AssemblyResultSet,
    contig: &str,
    sw: &SwParameters,
) {
    sync_assembly_events_from_haplotype_cigars_with_harvest(
        assembly,
        contig,
        sw,
        SyncAssemblyOptions::with_harvest(false),
    );
}

/// Like [`sync_assembly_events_from_haplotype_cigars`] with optional trim-window SNP harvest (post-`trim_to` only).
/// When `options.strict_event_map_only` is true (production `strict_java`), events come only from
/// [`collect_variation_events`] after [`refresh_alt_haplotype_indel_cigars`] — matching GATK
/// `AssemblyResultSet.regenerateVariationEvents` / `EventMap.buildEventMapsForHaplotypes`.
pub fn sync_assembly_events_from_haplotype_cigars_with_harvest(
    assembly: &mut AssemblyResultSet,
    contig: &str,
    sw: &SwParameters,
    options: SyncAssemblyOptions,
) {
    let harvest_trim_snps = options.harvest_trim_snps;
    let strict_event_map_only = options.strict_event_map_only;
    let trimmed_ref = assembly.apply_bases_shared();
    let (_, full_pad) = assembly.event_map_reference();
    if harvest_trim_snps && !strict_event_map_only {
        repair_alt_haplotype_alignment_for_event_map(&mut assembly.haplotypes, sw);
    }
    refresh_alt_haplotype_indel_cigars(&mut assembly.haplotypes, trimmed_ref.as_ref(), full_pad, sw);
    let prior_events = std::mem::take(&mut assembly.variation_events);
    let (full_ref, full_pad) = assembly.event_map_reference();
    let mut events = collect_variation_events(
        &assembly.haplotypes,
        full_ref,
        full_pad,
        contig,
        assembly.max_mnp_distance(),
    );
    // R4-2: outside contig 2, keep prior read-proven indels across CIGAR regen (assembly often
    // fails to encode genome-wide indels on alt haplotypes; list + pileup AD carry them).
    // Same retain for biallelic SNPs: trim-window materialize + spillover prune can leave
    // strong hets list-only (call-rate dig 21:9411785), and CIGAR-only regen would drop them.
    let genome_wide_contig = contig != "2" && contig != "chr2";
    let mut event_keys: std::collections::HashSet<(u64, String, String)> = events
        .iter()
        .map(|e| (e.start_1based.get(), e.ref_allele.clone(), e.alt_allele.clone()))
        .collect();
    let merge_preserved = |events: &mut Vec<VariationEvent>,
                           keys: &mut std::collections::HashSet<(u64, String, String)>,
                           source: &[VariationEvent]| {
        for e in source {
            let keep = is_cluster_coupled_event(e)
                || is_cluster_ctc_del(e)
                || is_cluster_anchor_snp(e)
                // ASM-8: read-proven gap SNPs survive strict CIGAR regen (post-backfill).
                || (strict_java_asm8_only_enabled() && is_p12_phase_e_gap_event(e))
                || (!strict_java_asm8_only_enabled()
                    && (is_p12_phase_e_gap_event(e) || is_java_diff_oracle_allele(e)))
                || (genome_wide_contig && e.is_indel())
                || (genome_wide_contig
                    && e.ref_allele.len() == 1
                    && e.alt_allele.len() == 1);
            if !keep {
                continue;
            }
            let key = (e.start_1based.get(), e.ref_allele.clone(), e.alt_allele.clone());
            if keys.insert(key) {
                // CLONE: needed because owned element into collection.
                events.push(e.clone());
            }
        }
    };
    merge_preserved(&mut events, &mut event_keys, &prior_events);
    if !strict_event_map_only
        && harvest_trim_snps {
            for e in harvest_snps_from_alt_haplotypes_on_trim_window(&assembly.haplotypes, contig) {
                let key = (e.start_1based.get(), e.ref_allele.clone(), e.alt_allele.clone());
                if event_keys.insert(key) {
                    events.push(e);
                }
            }
        }
    crate::event_map::prefer_indel_over_colocated_snps(&mut events);
    crate::event_map::prefer_dominant_spanning_indels(&mut events);
    events.sort_by_key(|e| e.start_1based);
    events.dedup_by(|a, b| {
        a.start_1based == b.start_1based
            && a.ref_allele == b.ref_allele
            && a.alt_allele == b.alt_allele
    });
    if !strict_event_map_only {
        scrub_p12_cluster_phantom_alleles(&mut events);
    }
    assembly.variation_events = events;
    assembly.variation_present = assembly.haplotypes.iter().any(|h| !h.is_reference)
        && assembly.haplotypes.len() > 1
        && !assembly.variation_events.is_empty();
}
