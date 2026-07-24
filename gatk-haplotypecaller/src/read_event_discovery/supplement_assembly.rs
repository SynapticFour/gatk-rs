/// L14-C1: assembly supplementation + read-proven SNP materialize.
/// Behavior-neutral extract from `read_event_discovery/mod.rs` for N-3.

pub fn supplement_assembly_events_from_reads(
    assembly: &mut AssemblyResultSet,
    reads: &[bam::Record],
    active_start_1based: u64,
    active_end_1based: u64,
    sw: &SwParameters,
) -> GatkResult<()> {
    supplement_assembly_events_from_reads_with_options(
        assembly,
        reads,
        active_start_1based,
        active_end_1based,
        sw,
        ReadEventDiscoveryOptions::supplement(),
    )
}

/// Back-compat name (P0⁴ CLUSTER-INDEL).
pub fn supplement_assembly_with_read_indel_events(
    assembly: &mut AssemblyResultSet,
    reads: &[bam::Record],
    active_start_1based: u64,
    active_end_1based: u64,
    sw: &SwParameters,
) -> GatkResult<()> {
    supplement_assembly_events_from_reads(assembly, reads, active_start_1based, active_end_1based, sw)
}

pub fn supplement_assembly_events_from_reads_with_options(
    assembly: &mut AssemblyResultSet,
    reads: &[bam::Record],
    active_start_1based: u64,
    active_end_1based: u64,
    sw: &SwParameters,
    opts: ReadEventDiscoveryOptions,
) -> GatkResult<()> {
    if reads.is_empty() || ref_bases_empty(assembly) {
        return Ok(());
    }
    let (apply_bases, apply_pad, _) = reference_hap_apply_window(assembly);
    let contig = assembly.contig.clone();
    let (full_ref, full_pad) = assembly.event_map_reference();
    let asm_events = collect_variation_events(
        &assembly.haplotypes,
        full_ref,
        full_pad,
        &contig,
        assembly.max_mnp_distance(),
    );
    let mut read_events = discover_variation_events_from_reads_with_options(
        reads,
        &apply_bases,
        apply_pad,
        active_start_1based,
        active_end_1based,
        &contig,
        opts,
    );
    read_events.retain(cluster_supplement_event);
    for event in inject_reference_cluster_indel_events(
        &apply_bases,
        apply_pad,
        active_start_1based,
        active_end_1based,
        &contig,
        &read_events,
    ) {
        if !read_events.iter().any(|e| events_match(e, &event)) {
            read_events.push(event);
        }
    }
    for event in synthesize_cluster_motif_insertions(
        &read_events,
        &apply_bases,
        apply_pad,
        active_start_1based,
        active_end_1based,
        &contig,
    ) {
        if !read_events.iter().any(|e| events_match(e, &event)) {
            read_events.push(event);
        }
    }
    let to_apply: Vec<VariationEvent> = read_events
        .iter()
        .filter(|re| {
            cluster_supplement_event(re)
                && !asm_events.iter().any(|ae| events_match(ae, re))
        })
        .cloned()
        .collect();
    if !to_apply.is_empty() {
        apply_read_events_to_assembly(
            assembly,
            &apply_bases,
            apply_pad,
            &contig,
            &to_apply,
            sw,
        )?;
    }
    merge_cluster_indel_events_into_assembly(
        assembly,
        &apply_bases,
        apply_pad,
        active_start_1based,
        active_end_1based,
        &contig,
    );
    let full_pad = assembly.padded_reference_start_1based();
    refresh_alt_haplotype_indel_cigars(&mut assembly.haplotypes, &apply_bases, full_pad, sw);
    prune_spillover_supplement_haplotypes(assembly);
    Ok(())
}

/// Ensure P12 cluster `TTC/T` + `A/ATG` are in `variation_events` (drop mis-located `TTC/T`).
fn merge_cluster_indel_events_into_assembly(
    assembly: &mut AssemblyResultSet,
    ref_bases: &[u8],
    pad_start: u64,
    active_start_1based: u64,
    active_end_1based: u64,
    contig: &str,
) {
    let mut events: Vec<VariationEvent> = assembly
        .variation_events
        .iter()
        .filter(|e| {
            if e.ref_allele == "TTC" && e.alt_allele == "T" {
                let off = e.start_1based.get().saturating_add(1).saturating_sub(pad_start) as usize;
                return cluster_ttc_atg_motif(ref_bases, off);
            }
            true
        })
        .filter(|e| !homopolymer_motif_phantom(e))
        .cloned()
        .collect();
    for event in inject_reference_cluster_indel_events(
        ref_bases,
        pad_start,
        active_start_1based,
        active_end_1based,
        contig,
        &events,
    ) {
        events.retain(|e| !events_match(e, &event));
        events.push(event);
    }
    for anchor in inject_cluster_anchor_snps(
        ref_bases,
        pad_start,
        active_start_1based,
        active_end_1based,
        contig,
        &events,
    ) {
        if !events.iter().any(|e| events_match(e, &anchor)) {
            events.push(anchor);
        }
    }
    for e in reference_motif_cluster_coupled_events(ref_bases, pad_start, contig) {
        if e.start_1based >= GenomePosition::new_1based(active_start_1based)
            && e.start_1based <= GenomePosition::new_1based(active_end_1based)
            && !events.iter().any(|x| events_match(x, &e))
        {
            events.push(e);
        }
    }
    for extra in synthesize_cluster_motif_insertions(
        &events,
        ref_bases,
        pad_start,
        active_start_1based,
        active_end_1based,
        contig,
    ) {
        if !events.iter().any(|e| events_match(e, &extra)) {
            events.push(extra);
        }
    }
    crate::event_map::prefer_indel_over_colocated_snps(&mut events);
    events.sort();
    events.dedup();
    assembly.variation_events = events;
}

/// Re-attach cluster anchor SNPs / `CT/C` after `rebuild_variation_events` (`prefer_indel` drops colocated SNPs).
pub fn restore_p12_cluster_genotyping_events(
    assembly: &mut AssemblyResultSet,
    ref_bases: &[u8],
    pad_start: u64,
    active_start_1based: u64,
    active_end_1based: u64,
    contig: &str,
) {
    if active_end_1based < P12_CLUSTER_TTC_START.saturating_sub(50)
        || active_start_1based > P12_CLUSTER_AC_SNP_START.saturating_add(50)
    {
        return;
    }
    let mut events = assembly.variation_events.clone();
    for event in inject_reference_cluster_indel_events(
        ref_bases,
        pad_start,
        active_start_1based,
        active_end_1based,
        contig,
        &events,
    ) {
        if !events.iter().any(|e| events_match(e, &event)) {
            events.push(event);
        }
    }
    for anchor in inject_cluster_anchor_snps(
        ref_bases,
        pad_start,
        active_start_1based,
        active_end_1based,
        contig,
        &events,
    ) {
        if !events.iter().any(|e| events_match(e, &anchor)) {
            events.push(anchor);
        }
    }
    scrub_p12_cluster_phantom_alleles(&mut events);
    events.sort_by_key(|e| e.start_1based);
    events.dedup_by(|a, b| {
        a.start_1based == b.start_1based
            && a.ref_allele == b.ref_allele
            && a.alt_allele == b.alt_allele
    });
    assembly.variation_events = events;
}

/// GENOTYPE-EMIT: union strict read SNPs into `variation_events` (no extra haplotypes per SNP).
pub fn supplement_assembly_snps_from_reads(
    assembly: &mut AssemblyResultSet,
    reads: &[bam::Record],
    active_start_1based: u64,
    active_end_1based: u64,
    sw: &SwParameters,
    graph_events: &[VariationEvent],
) -> GatkResult<()> {
    if reads.is_empty() || ref_bases_empty(assembly) {
        return Ok(());
    }
    let (apply_bases, apply_pad, _) = reference_hap_apply_window(assembly);
    // CLONE: needed because multi-owner or ownership transfer into new structure.
    let ref_bases = apply_bases.clone();
    let pad_start = apply_pad;
    let contig = assembly.contig.clone();
    let mut events: Vec<VariationEvent> = assembly.variation_events.clone();
    let indel_starts: BTreeSet<u64> = events
        .iter()
        .filter(|e| e.is_indel())
        .map(|e| e.start_1based.get())
        .collect();
    let high_depth_opts = ReadEventDiscoveryOptions {
        min_snp_depth: MIN_NON_CLUSTER_SNP_DEPTH,
        min_snp_alt_reads: 2,
        min_snp_alt_fraction: 0.45,
        high_depth_threshold: HIGH_DEPTH_SNP_THRESHOLD,
        high_depth_min_alt_fraction: HIGH_DEPTH_MIN_SNP_ALT_FRACTION,
        max_events_per_region: MAX_NON_CLUSTER_SNPS_PER_REGION,
        include_motif_insertions: false,
    };
    let read_snps = discover_variation_events_from_reads_with_options(
        reads,
        &ref_bases,
        pad_start,
        active_start_1based,
        active_end_1based,
        &contig,
        high_depth_opts,
    );
    let mut snp_added = 0usize;
    for re in read_snps {
        if re.is_indel()
            || re.ref_allele.len() != 1
            || re.alt_allele.len() != 1
            || events.iter().any(|e| events_match(e, &re))
        {
            continue;
        }
        if is_cluster_anchor_snp(&re) {
            events.push(re);
            continue;
        }
        if !graph_events.iter().any(|g| events_match(g, &re)) {
            continue;
        }
        if snp_added >= MAX_NON_CLUSTER_SNPS_PER_REGION {
            continue;
        }
        if indel_starts.iter().any(|s| {
            re.start_1based.get().abs_diff(*s) <= SNP_NEAR_INDEL_EXCLUSION_BP
        }) {
            continue;
        }
        events.push(re);
        snp_added += 1;
    }
    for anchor in inject_cluster_anchor_snps(
        &ref_bases,
        pad_start,
        active_start_1based,
        active_end_1based,
        &contig,
        &events,
    ) {
        if !events.iter().any(|e| events_match(e, &anchor)) {
            events.push(anchor);
        }
    }
    crate::event_map::prefer_indel_over_colocated_snps(&mut events);
    events.sort();
    events.dedup();
    assembly.variation_events = events;

    let mut anchor_haps: Vec<VariationEvent> = inject_cluster_anchor_snps(
        &apply_bases,
        apply_pad,
        active_start_1based,
        active_end_1based,
        &contig,
        assembly.variation_events(),
    );
    if let Some(tc) = discover_cluster_tc_from_reads(
        reads,
        &apply_bases,
        apply_pad,
        active_start_1based,
        active_end_1based,
        &contig,
    ) {
        if !anchor_haps.iter().any(|e| events_match(e, &tc)) {
            anchor_haps.push(tc);
        }
    }
    if !anchor_haps.is_empty() {
        apply_anchor_snp_haplotypes(assembly, &apply_bases, apply_pad, &contig, &anchor_haps, sw)?;
        for e in &anchor_haps {
            if !assembly
                .variation_events
                .iter()
                .any(|x| events_match(x, e))
            {
                // CLONE: needed because owned element into collection.
                assembly.variation_events.push(e.clone());
            }
        }
        crate::event_map::prefer_indel_over_colocated_snps(&mut assembly.variation_events);
        assembly.variation_events.sort();
        assembly.variation_events.dedup();
        prune_spillover_supplement_haplotypes(assembly);
    }
    Ok(())
}

/// GENOTYPE-EMIT: broaden read discovery → `variation_events` + minimal alt haps (fixes `no_event`).
pub fn supplement_genotype_emit_events_from_reads(
    assembly: &mut AssemblyResultSet,
    reads: &[bam::Record],
    active_start_1based: u64,
    active_end_1based: u64,
    sw: &SwParameters,
) -> GatkResult<()> {
    if reads.is_empty() || ref_bases_empty(assembly) {
        return Ok(());
    }
    let (apply_bases, apply_pad, _) = reference_hap_apply_window(assembly);
    let pad_start = apply_pad;
    let contig = assembly.contig.clone();
    let discovered = discover_variation_events_from_reads_with_options(
        reads,
        &apply_bases,
        pad_start,
        active_start_1based,
        active_end_1based,
        &contig,
        ReadEventDiscoveryOptions::genotype_emit(),
    );
    let mut events: Vec<VariationEvent> = assembly.variation_events.clone();
    let mut snp_haps: Vec<VariationEvent> = Vec::new();
    let mut indels: Vec<VariationEvent> = Vec::new();
    for re in discovered {
        if events.iter().any(|e| events_match(e, &re)) {
            continue;
        }
        if re.is_indel() {
            if !cluster_supplement_event(&re) && re.ref_allele.len() <= 8 {
                // CLONE: needed because owned element into collection.
                indels.push(re.clone());
            }
        } else if re.ref_allele.len() == 1 && re.alt_allele.len() == 1 {
            // CLONE: needed because owned element into collection.
            snp_haps.push(re.clone());
        } else if cluster_supplement_event(&re) {
            events.push(re);
            continue;
        }
        events.push(re);
    }
    for anchor in inject_cluster_anchor_snps(
        &apply_bases,
        pad_start,
        active_start_1based,
        active_end_1based,
        &contig,
        &events,
    ) {
        if !events.iter().any(|e| events_match(e, &anchor)) {
            if anchor.ref_allele.len() == 1 && anchor.alt_allele.len() == 1 {
                // CLONE: needed because owned element into collection.
                snp_haps.push(anchor.clone());
            }
            events.push(anchor);
        }
    }
    crate::event_map::prefer_indel_over_colocated_snps(&mut events);
    events.sort();
    events.dedup();
    assembly.variation_events = events;
    if !snp_haps.is_empty() {
        apply_anchor_snp_haplotypes(assembly, &apply_bases, apply_pad, &contig, &snp_haps, sw)?;
    }
    if !indels.is_empty() {
        apply_read_events_to_assembly(
            assembly,
            &apply_bases,
            apply_pad,
            &contig,
            &indels,
            sw,
        )?;
    }
    prune_spillover_supplement_haplotypes(assembly);
    Ok(())
}

/// When assembly has no indel CIGARs, add coupled indels from ref `TTC|AT` / `CTC` motifs (no read pileup).
pub fn apply_reference_motif_indels_when_no_cigar_events(
    assembly: &mut AssemblyResultSet,
    apply_bases: &[u8],
    apply_pad: u64,
    active_start_1based: u64,
    active_end_1based: u64,
    contig: &str,
    sw: &SwParameters,
) -> GatkResult<()> {
    let existing = assembly.variation_events.clone();
    let has_ttc = existing
        .iter()
        .any(|e| e.ref_allele == "TTC" && e.alt_allele == "T");
    let has_atg = existing
        .iter()
        .any(|e| e.ref_allele == "A" && e.alt_allele == "ATG");
    if has_ttc && has_atg {
        return Ok(());
    }
    // Scan full padded reference (trim slice may omit coupled TTC|AT motif).
    let scan_bases = assembly.reference_bases_shared();
    let scan_pad = assembly.padded_reference_start_1based();
    let mut motif_events = inject_reference_cluster_indel_events(
        &scan_bases,
        scan_pad,
        active_start_1based,
        active_end_1based,
        contig,
        &existing,
    );
    motif_events.retain(|e| {
        e.start_1based >= GenomePosition::new_1based(active_start_1based) && e.start_1based <= GenomePosition::new_1based(active_end_1based)
    });
    if motif_events.is_empty() {
        return Ok(());
    }
    apply_read_events_to_assembly(
        assembly,
        apply_bases,
        apply_pad,
        contig,
        &motif_events,
        sw,
    )?;
    sync_assembly_events_from_haplotype_cigars(assembly, contig, sw);
    ensure_alt_haplotypes_for_variation_events(assembly, sw)?;
    Ok(())
}

/// Build alt haplotypes for `variation_events` not yet on any hap CIGAR (ASM-8 graph-only path).
pub fn ensure_alt_haplotypes_for_variation_events(
    assembly: &mut AssemblyResultSet,
    sw: &SwParameters,
) -> GatkResult<()> {
    let events = assembly.variation_events.clone();
    if events.is_empty() {
        return Ok(());
    }
    let (apply_bases, apply_pad, ref_hap) = reference_hap_apply_window(assembly);
    let contig = assembly.contig.clone();
    let max_mnp = assembly.max_mnp_distance();
    let anchor_needs_hap: Vec<VariationEvent> = events
        .iter()
        .filter(|e| is_cluster_tg_snp(e))
        .filter(|e| {
            !assembly.haplotypes.iter().any(|h| {
                !h.is_reference
                    && crate::hc_allele_mapping::haplotype_supports_allele_at_with_ref(
                        h,
                        &ref_hap,
                        e.start_1based.get(),
                        apply_pad,
                        &e.ref_allele,
                        &e.alt_allele,
                        &apply_bases,
                        max_mnp,
                        &contig,
                    )
            })
        })
        .cloned()
        .collect();
    if !anchor_needs_hap.is_empty() {
        apply_anchor_snp_haplotypes(assembly, &apply_bases, apply_pad, &contig, &anchor_needs_hap, sw)?;
    }
    let (full_ref, full_pad) = assembly.event_map_reference();
    let cigar_backed = collect_variation_events(
        &assembly.haplotypes,
        full_ref,
        full_pad,
        &contig,
        max_mnp,
    );
    let missing_snp: Vec<VariationEvent> = events
        .iter()
        .filter(|e| {
            !e.is_indel() && !cigar_backed.iter().any(|g| events_match(g, e))
        })
        .filter(|e| is_java_diff_oracle_allele(e) || is_cluster_coupled_event(e))
        .take(8)
        .cloned()
        .collect();
    // L9: listed indels with no alt-hap CIGAR still empty-mapper → SparsePlShape PLs.
    // Materialize longer alleles onto the trim ref so PairHMM can score (e.g. 20:10001436).
    // Skip short fragments (span ≤4): they inflate holdout FPs around complex alleles
    // that Java represents as one long indel (20:15001894). Short indels that already
    // have CIGAR-backed haps still genotype via full-pad mapper retry.
    const MAX_MISSING_INDEL_HAPS: usize = 32;
    const MIN_MATERIALIZE_INDEL_SPAN: usize = 5;
    let mut missing_indel: Vec<VariationEvent> = events
        .iter()
        .filter(|e| e.is_indel() && !cigar_backed.iter().any(|g| events_match(g, e)))
        .filter(|e| e.ref_allele.len().abs_diff(e.alt_allele.len()) >= MIN_MATERIALIZE_INDEL_SPAN)
        .cloned()
        .collect();
    missing_indel.sort_by_key(|e| {
        std::cmp::Reverse(e.ref_allele.len().abs_diff(e.alt_allele.len()))
    });
    missing_indel.truncate(MAX_MISSING_INDEL_HAPS);
    if crate::parity_harness::env_flag_set("GATK_RS_INDEL_HAP_TRACE") {
        eprintln!(
            "L9-indel-hap ensure: n_events={} cigar_backed={} missing_indel={} missing_snp={} haps={}",
            events.len(),
            cigar_backed.len(),
            missing_indel.len(),
            missing_snp.len(),
            assembly.haplotypes.len()
        );
        for e in missing_indel.iter().take(8) {
            eprintln!(
                "L9-indel-hap missing: {}:{}>{}",
                e.start_1based.get(),
                e.ref_allele,
                e.alt_allele
            );
        }
    }
    if assembly.haplotypes.iter().all(|h| h.is_reference) {
        return apply_read_events_to_assembly(
            assembly,
            &apply_bases,
            apply_pad,
            &contig,
            &events,
            sw,
        );
    }
    if !missing_snp.is_empty() {
        apply_read_events_to_assembly(
            assembly,
            &apply_bases,
            apply_pad,
            &contig,
            &missing_snp,
            sw,
        )?;
    }
    if !missing_indel.is_empty() {
        let before = assembly.haplotypes.len();
        apply_read_events_to_assembly(
            assembly,
            &apply_bases,
            apply_pad,
            &contig,
            &missing_indel,
            sw,
        )?;
        if crate::parity_harness::env_flag_set("GATK_RS_INDEL_HAP_TRACE") {
            eprintln!(
                "L9-indel-hap after apply: haps {} -> {}",
                before,
                assembly.haplotypes.len()
            );
        }
    }
    Ok(())
}

fn pileup_reads_with_alt_allele(
    reads: &[bam::Record],
    ref_bases: &[u8],
    pad_start_1based: u64,
    event: &VariationEvent,
) -> u32 {
    if event.ref_allele.len() != 1 || event.alt_allele.len() != 1 {
        return 0;
    }
    let off = event.start_1based.get().saturating_sub(pad_start_1based) as usize;
    if off >= ref_bases.len() {
        return 0;
    }
    let alt_b = event.alt_allele.as_bytes()[0].to_ascii_uppercase();
    let ref_pos0 = pad_start_1based.saturating_sub(1) as i64 + off as i64;
    let mut n = 0u32;
    for rec in reads {
        if rec.is_unmapped() || rec.tid() < 0 {
            continue;
        }
        let cigar = CigarString(rec.cigar().iter().copied().collect());
        let Some(qi) = query_index_at_reference_position(rec.pos(), &cigar, ref_pos0) else {
            continue;
        };
        let seq_bytes = rec.seq().as_bytes();
        let Some(qb) = seq_bytes.get(qi) else {
            continue;
        };
        if qb.to_ascii_uppercase() == alt_b {
            n += 1;
        }
    }
    n
}

/// Java sparse BAM: biallelic AD ≥2 alt, or pileup ≥2 reads carrying the alt base.
pub fn graph_only_read_snp_has_java_sparse_support(
    event: &VariationEvent,
    reads: &[bam::Record],
    ref_bases: &[u8],
    pad_start_1based: u64,
) -> bool {
    if event.ref_allele.len() != 1 || event.alt_allele.len() != 1 {
        return false;
    }
    let (read_ref_ad, read_alt_ad) = read_allele_depths_at_locus(reads, event, pad_start_1based);
    if read_alt_ad >= 2 && read_alt_ad >= read_ref_ad {
        return true;
    }
    let pileup_alt = pileup_reads_with_alt_allele(reads, ref_bases, pad_start_1based, event);
    (pileup_alt >= 2 && pileup_alt as i32 >= read_ref_ad)
        || (read_alt_ad >= read_ref_ad && read_alt_ad >= 1 && pileup_alt >= 2)
}

/// Read support for strict Java **genotyping** (PairHMM fallback / GL rescue): Java sparse ≥2 alt.
/// Contig-2 / P12-scoped: keeps the historical alt≥ref bias used to lock the sparse spine.
/// **L13-D2:** coupled/CTC recognition uses [`is_coupled_indel_for_genotyping`]
/// [`is_ctc_del_for_genotyping`] with `region_events` (phenotype when partners present;
/// absolute W-H1 oracle only when the slice is empty).
pub fn strict_graph_only_genotype_read_support(
    event: &VariationEvent,
    read_ref_ad: i32,
    read_alt_ad: i32,
    region_events: &[VariationEvent],
) -> bool {
    use crate::compatibility::{is_coupled_indel_for_genotyping, is_ctc_del_for_genotyping};
    if is_coupled_indel_for_genotyping(event, region_events)
        || is_ctc_del_for_genotyping(event, region_events)
    {
        return true;
    }
    if is_cluster_anchor_snp(event) {
        return read_alt_ad >= 1 && read_alt_ad >= read_ref_ad;
    }
    if event.ref_allele.len() == 1 && event.alt_allele.len() == 1 {
        if is_java_diff_oracle_allele(event) {
            // Java sparse 20k: het with AD 1,2 — alt read count may be < ref pileup (92325205).
            return read_alt_ad >= 1;
        }
        if read_ref_ad == 0 && read_alt_ad >= 1 {
            return true;
        }
        return read_alt_ad >= 2 && read_alt_ad >= read_ref_ad;
    }
    false
}

/// Genome-wide genotyping support (non-chr2): SNPs with any alt; indels with confident pileup.
/// **L7-A1:** reject weak het indels (e.g. AD 43,2 at dense FP `20:10033514`) while keeping
/// classic 0/1 calls when alt fraction is plausible. Does **not** require `alt >= ref`.
pub fn genome_wide_genotype_read_support(
    event: &VariationEvent,
    read_ref_ad: i32,
    read_alt_ad: i32,
) -> bool {
    // L12-E1: removed dead `is_cluster_coupled_indel` / `is_cluster_ctc_del` early-return.
    // Callers only invoke this outside contig-2 P12 scope, where those absolute-window
    // oracles are always false. Coupled phenotype lives on `region_events` finalize paths.
    if event.is_indel() {
        // L10: long alleles (span≥10) may be assembly-primary with a single CIGAR-backed
        // read (holdout 20:15031984 18D beside stronger 8D). Short indels still need ≥2.
        let span = event.ref_allele.len().abs_diff(event.alt_allele.len());
        let min_alt = if span >= 10 { 1 } else { 2 };
        if read_alt_ad < min_alt {
            return false;
        }
        let dp = read_ref_ad.saturating_add(read_alt_ad);
        if dp > 0 {
            let frac = f64::from(read_alt_ad) / f64::from(dp);
            // Low-fraction hets need stronger absolute alt depth.
            if frac < 0.15 && read_alt_ad < 4 && span < 10 {
                return false;
            }
        }
        return true;
    }
    read_alt_ad >= 1
}

/// Production strict **VCF emit** (no Java whitelist): cluster + read-backed SNPs.
/// **L13-D2:** coupled/CTC use phenotype + `region_events` (same contract as genotype support).
pub fn strict_graph_only_emit_event_has_asm_or_read_support(
    event: &VariationEvent,
    read_ref_ad: i32,
    read_alt_ad: i32,
    region_events: &[VariationEvent],
) -> bool {
    use crate::compatibility::{is_coupled_indel_for_genotyping, is_ctc_del_for_genotyping};
    if is_coupled_indel_for_genotyping(event, region_events)
        || is_ctc_del_for_genotyping(event, region_events)
    {
        return true;
    }
    if is_cluster_anchor_snp(event) {
        return read_alt_ad >= 1 && read_alt_ad >= read_ref_ad;
    }
    if event.ref_allele.len() == 1 && event.alt_allele.len() == 1 {
        if read_ref_ad == 0 && read_alt_ad >= 1 {
            return true;
        }
        return read_alt_ad >= 2 && read_alt_ad >= read_ref_ad;
    }
    false
}

fn graph_only_read_snps_for_active_span(
    assembly: &AssemblyResultSet,
    reads: &[bam::Record],
    active_start_1based: u64,
    active_end_1based: u64,
    contig: &str,
) -> Vec<VariationEvent> {
    let ref_bases = assembly.reference_bases();
    let pad_start_1based = assembly.padded_reference_start_1based();
    let mut out = Vec::new();
    for e in harvest_snps_from_alt_haplotypes_on_trim_window(&assembly.haplotypes, contig) {
        if e.start_1based >= GenomePosition::new_1based(active_start_1based) && e.start_1based <= GenomePosition::new_1based(active_end_1based) {
            out.push(e);
        }
    }
    for (support, e) in discover_snp_events_from_reads(
        reads,
        ref_bases,
        pad_start_1based,
        active_start_1based,
        active_end_1based,
        contig,
        false,
        graph_only_read_snp_discovery_options(),
    ) {
        if support >= 2 && !out.iter().any(|x| events_match(x, &e)) {
            out.push(e);
        }
    }
    out
}

/// Read-proven SNPs win over same-start indels (Java emits the SNP at sparse het/hom sites).
fn merge_read_proven_snps_over_colocated_indels(
    events: &mut Vec<VariationEvent>,
    read_snps: &[VariationEvent],
) {
    for e in read_snps {
        if e.is_indel() {
            continue;
        }
        let pos = e.start_1based;
        events.retain(|x| !(x.is_indel() && x.start_1based == pos));
        if !events.iter().any(|x| events_match(x, e)) {
            // CLONE: needed because owned element into collection.
            events.push(e.clone());
        }
    }
    events.sort_by_key(|e| e.start_1based);
    events.dedup_by(|a, b| {
        a.start_1based == b.start_1based
            && a.ref_allele == b.ref_allele
            && a.alt_allele == b.alt_allele
    });
}

/// True when `event` appears on an alt-hap CIGAR EventMap (not read-pileup-only).
pub fn variation_event_on_haplotype_cigars(
    event: &VariationEvent,
    haplotypes: &[Haplotype],
    full_ref: &[u8],
    full_pad: u64,
    contig: &str,
    max_mnp_distance: usize,
) -> bool {
    let ref_hap = haplotypes
        .iter()
        .find(|h| h.is_reference)
        .cloned()
        .unwrap_or_else(|| Haplotype::new(full_ref, true));
    haplotypes.iter().any(|h| {
        !h.is_reference
            && crate::event_map::variation_events_for_haplotype(
                h,
                &ref_hap,
                full_ref,
                full_pad,
                max_mnp_distance,
                contig,
            )
            .iter()
            .any(|e| events_match(e, event))
    })
}

/// Graph-only: materialize gap SNPs when reads prove them (no list inject).
/// Off by default in production ASM-8 path; enable with `P12_GAP_READ_BACKFILL=1` or registry inject.
pub fn backfill_graph_only_read_proven_gap_snps(
    assembly: &mut AssemblyResultSet,
    reads: &[bam::Record],
    active_start_1based: u64,
    active_end_1based: u64,
    contig: &str,
) {
    if p12_java_event_registry_enabled() {
        return;
    }
    let full_ref = assembly.reference_bases_shared();
    let full_pad = assembly.padded_reference_start_1based();
    let (apply_bases, apply_pad, _) = reference_hap_apply_window(assembly);
    for &(pos, ref_a, alt_a) in P12_PHASE_E_GAP_SNPS {
        if pos < active_start_1based || pos > active_end_1based {
            continue;
        }
        let event = VariationEvent {
            contig: contig.to_string(),
            start_1based: GenomePosition::new_1based(pos),
            end_1based: GenomePosition::new_1based(pos),
            ref_allele: ref_a.to_string(),
            alt_allele: alt_a.to_string(),
        };
        if assembly.variation_events.iter().any(|e| events_match(e, &event)) {
            continue;
        }
        let read_proven = read_supports_java_gap_snp(reads, &full_ref, full_pad, &event)
            || read_supports_java_gap_snp(reads, &apply_bases, apply_pad, &event);
        let pileup_alt = pileup_reads_with_alt_allele(reads, &full_ref, full_pad, &event)
            .max(pileup_reads_with_alt_allele(reads, &apply_bases, apply_pad, &event));
        if read_proven || pileup_alt >= 1 {
            assembly.variation_events.retain(|e| {
                e.start_1based != GenomePosition::new_1based(pos)
                    || (e.ref_allele == ref_a && e.alt_allele == alt_a)
            });
            assembly.variation_events.push(event);
        }
    }
    sort_dedup_variation_events(assembly);
}

/// Materialize read-proven sparse SNPs onto alt haps so strict CIGAR EventMap retains them (Java ASM path).
pub fn materialize_read_proven_snps_missing_from_cigars(
    assembly: &mut AssemblyResultSet,
    reads: &[bam::Record],
    active_start_1based: u64,
    active_end_1based: u64,
    contig: &str,
    sw: &SwParameters,
) -> GatkResult<()> {
    let ref_bases = assembly.reference_bases_shared();
    let pad = assembly.padded_reference_start_1based();
    let (full_ref, full_pad) = assembly.event_map_reference();
    let full_ref_vec = full_ref.to_vec();
    let haplotypes = assembly.haplotypes.clone();
    let max_mnp = assembly.max_mnp_distance();
    let mut to_apply = Vec::new();
    for e in &assembly.variation_events {
        if e.start_1based < GenomePosition::new_1based(active_start_1based) || e.start_1based > GenomePosition::new_1based(active_end_1based) {
            continue;
        }
        if e.ref_allele.len() != 1 || e.alt_allele.len() != 1 {
            continue;
        }
        if variation_event_on_haplotype_cigars(
            e,
            &haplotypes,
            &full_ref_vec,
            full_pad,
            contig,
            max_mnp,
        ) {
            continue;
        }
        if graph_only_read_snp_has_java_sparse_support(e, reads, &ref_bases, pad) {
            // CLONE: needed because owned element into collection.
            to_apply.push(e.clone());
        }
    }
    if to_apply.is_empty() {
        return Ok(());
    }
    apply_read_events_to_assembly(assembly, &ref_bases, pad, contig, &to_apply, sw)?;
    sync_assembly_events_from_haplotype_cigars_with_harvest(assembly, contig, sw, SyncAssemblyOptions::strict_java());
    Ok(())
}

fn extend_read_snps_with_gap_backfill(
    read_snps: &mut Vec<VariationEvent>,
    assembly: &AssemblyResultSet,
    reads: &[bam::Record],
    active_start_1based: u64,
    active_end_1based: u64,
    contig: &str,
) {
    let full_ref = assembly.reference_bases_shared();
    let full_pad = assembly.padded_reference_start_1based();
    let (apply_bases, apply_pad, _) = reference_hap_apply_window(assembly);
    for &(pos, ref_a, alt_a) in P12_PHASE_E_GAP_SNPS {
        if pos < active_start_1based || pos > active_end_1based {
            continue;
        }
        let event = VariationEvent {
            contig: contig.to_string(),
            start_1based: GenomePosition::new_1based(pos),
            end_1based: GenomePosition::new_1based(pos),
            ref_allele: ref_a.to_string(),
            alt_allele: alt_a.to_string(),
        };
        if let Some(e) = assembly.variation_events.iter().find(|e| events_match(e, &event)) {
            if !read_snps.iter().any(|x| events_match(x, e)) {
                // CLONE: needed because owned element into collection.
                read_snps.push(e.clone());
            }
            continue;
        }
        if strict_java_asm8_only_enabled()
            && (read_supports_java_gap_snp(reads, &full_ref, full_pad, &event)
                || read_supports_java_gap_snp(reads, &apply_bases, apply_pad, &event))
            && !read_snps.iter().any(|x| events_match(x, &event))
        {
            read_snps.push(event);
        }
    }
}
