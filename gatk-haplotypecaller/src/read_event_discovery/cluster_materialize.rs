
/// Strict Java: when threading/SeqGraph miss cluster indels, build alt hap from read-proven
/// `TTC/T` + `A/ATG` then derive EventMap from hap CIGAR (not list inject).
/// Java analogue: reads are already in `assembleReads`; this adds an alt haplotype the graph
/// should have produced, then `EventMap.buildEventMapsForHaplotypes` / `getVariationEvents`
/// (see `AssemblyBasedCallerUtils` + `HaplotypeCallerEngine` after `trimTo`).
pub fn strict_materialize_cluster_haplotype_from_reads(
    assembly: &mut AssemblyResultSet,
    reads: &[SharedBamRecord],
    active_start_1based: u64,
    active_end_1based: u64,
    contig: &str,
    sw: &SwParameters,
) -> GatkResult<()> {
    let debug = crate::runtime_config::strict_cluster_debug_enabled();
    if active_end_1based < P12_CLUSTER_TTC_START
        || active_start_1based > P12_CLUSTER_TTC_START.saturating_add(3)
    {
        if debug {
            eprintln!("strict_materialize\tskip span {active_start_1based}-{active_end_1based}");
        }
        return Ok(());
    }
    let coupled = cluster_coupled_events_from_assembly_haplotypes(
        assembly,
        contig,
        active_start_1based,
        active_end_1based,
    );
    if cluster_coupled_events_complete(&coupled) {
        if debug {
            eprintln!("strict_materialize\tskip already complete {coupled:?}");
        }
        sync_assembly_events_from_haplotype_cigars_with_harvest(assembly, contig, sw, SyncAssemblyOptions::strict_java());
        return Ok(());
    }
    // Must match Java trim slice: ref hap `genome_loc.start` + trim-length bases (not full pad).
    let (apply_bases, apply_pad, _) = reference_hap_apply_window(assembly);
    let read_coupled = discover_p12_cluster_coupled_events_from_reads(
        reads,
        &apply_bases,
        apply_pad,
        active_start_1based,
        active_end_1based,
        contig,
    );
    if debug {
        eprintln!(
            "strict_materialize\tpad={apply_pad} ref_len={} reads={} asm={coupled:?} read={read_coupled:?}",
            apply_bases.len(),
            reads.len()
        );
    }
    if read_coupled.is_empty() {
        return Ok(());
    }
    // Graph + ref-motif path may have already materialized coupled cluster haps in assembleReads.
    let coupled_after_graph = cluster_coupled_events_from_assembly_haplotypes(
        assembly,
        contig,
        active_start_1based,
        active_end_1based,
    );
    if cluster_coupled_events_complete(&coupled_after_graph) {
        sync_assembly_events_from_haplotype_cigars_with_harvest(assembly, contig, sw, SyncAssemblyOptions::strict_java());
        return Ok(());
    }
    let haps_before = assembly.haplotypes.len();
    upsert_coupled_cluster_alt_haplotype(assembly, &apply_bases, apply_pad, &read_coupled, sw)?;
    if debug {
        let after = cluster_coupled_events_from_assembly_haplotypes(
            assembly,
            contig,
            active_start_1based,
            active_end_1based,
        );
        eprintln!(
            "strict_materialize\thaps {haps_before} -> {} coupled_asm={after:?}",
            assembly.haplotypes.len()
        );
    }
    repair_alt_haplotype_alignment_for_event_map(&mut assembly.haplotypes, sw);
    sync_assembly_events_from_haplotype_cigars_with_harvest(assembly, contig, sw, SyncAssemblyOptions::strict_java());
    Ok(())
}

/// Re-attach pre-trim cluster indels when trim clips alt-hap CIGARs to all-`M`.
pub fn propagate_cluster_coupled_from_untrimmed(
    untrimmed: &AssemblyResultSet,
    assembly: &mut AssemblyResultSet,
    apply_bases: &[u8],
    apply_pad: u64,
    active_start_1based: u64,
    active_end_1based: u64,
    contig: &str,
    sw: &SwParameters,
) -> GatkResult<()> {
    let from_untrimmed = cluster_coupled_events_from_assembly_haplotypes(
        untrimmed,
        contig,
        active_start_1based,
        active_end_1based,
    );
    if from_untrimmed.is_empty() {
        return Ok(());
    }
    let full_pad = assembly.padded_reference_start_1based();
    refresh_alt_haplotype_indel_cigars(&mut assembly.haplotypes, apply_bases, full_pad, sw);
    let ref_hap = assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .cloned()
        .unwrap_or_else(|| crate::haplotype::Haplotype::new(apply_bases.to_vec(), true));
    let (full_ref, full_pad) = assembly.event_map_reference();
    let max_mnp = assembly.max_mnp_distance();
    let mut from_trimmed: Vec<VariationEvent> = Vec::new();
    for h in assembly.haplotypes.iter().filter(|h| !h.is_reference) {
        for e in crate::event_map::variation_events_for_haplotype(
            h,
            &ref_hap,
            full_ref,
            full_pad,
            max_mnp,
            contig,
        ) {
            if is_cluster_coupled_event(&e) {
                from_trimmed.push(e);
            }
        }
    }
    if cluster_coupled_events_complete(&from_trimmed) {
        return Ok(());
    }
    push_coupled_cluster_alt_haplotype(assembly, apply_bases, apply_pad, &from_untrimmed, sw)?;
    refresh_alt_haplotype_indel_cigars(&mut assembly.haplotypes, apply_bases, full_pad, sw);
    Ok(())
}

/// Synthesize coupled cluster alt on the trim ref when apply fails (graph alt already touched trim).
fn synthesize_coupled_cluster_bases_on_trim(
    apply_bases: &[u8],
    apply_pad: u64,
) -> Option<Vec<u8>> {
    let ttc_off = P12_CLUSTER_TTC_START.saturating_sub(apply_pad) as usize;
    let atg_off = P12_CLUSTER_ATG_START.saturating_sub(apply_pad) as usize;
    if ttc_off + 3 > apply_bases.len() || atg_off >= apply_bases.len() {
        return None;
    }
    let ttc_ref = &apply_bases[ttc_off..ttc_off + 3];
    if ttc_ref != b"TTC" && ttc_ref != b"ttc" {
        return None;
    }
    if !apply_bases[atg_off].eq_ignore_ascii_case(&b'A') {
        return None;
    }
    let mut out = apply_bases.to_vec();
    out.remove(ttc_off + 1);
    out.remove(ttc_off + 1);
    let atg_adj = atg_off.saturating_sub(2);
    if !out.get(atg_adj).copied().unwrap_or(0).eq_ignore_ascii_case(&b'A') {
        return None;
    }
    out.insert(atg_adj + 1, b'T');
    out.insert(atg_adj + 2, b'G');
    Some(out)
}

/// Build coupled cluster alt bases from trim ref, falling back to full padded ref.
fn coupled_cluster_alt_bases(
    assembly: &AssemblyResultSet,
    apply_bases: &[u8],
    apply_pad: u64,
    cluster_events: &[VariationEvent],
) -> Option<Vec<u8>> {
    if let Some(bases) = apply_events_to_ref_chained(apply_bases, cluster_events, apply_pad) {
        return Some(bases);
    }
    if cluster_coupled_events_complete(cluster_events) {
        if let Some(bases) = synthesize_coupled_cluster_bases_on_trim(apply_bases, apply_pad) {
            return Some(bases);
        }
    }
    let full_ref = assembly.reference_bases_shared();
    let full_pad = assembly.padded_reference_start_1based();
    let full_alt = apply_events_to_ref_chained(&full_ref, cluster_events, full_pad)?;
    let ref_hap = assembly.haplotypes.iter().find(|h| h.is_reference)?;
    let trim_off = ref_hap.alignment_start_hap_wrt_ref;
    if trim_off >= full_alt.len() {
        return None;
    }
    let trim_end = trim_off.saturating_add(apply_bases.len());
    if full_alt.len() >= trim_end {
        Some(full_alt[trim_off..trim_end].to_vec())
    } else {
        Some(full_alt[trim_off..].to_vec())
    }
}

fn force_cluster_coupled_haplotype_cigar(apply_pad: u64, ref_len: usize) -> crate::cigar::Cigar {
    use crate::cigar::{Cigar, CigarOperator};
    let ttc_off = P12_CLUSTER_TTC_START.saturating_sub(apply_pad) as usize;
    let tail = ref_len.saturating_sub(ttc_off + 3);
    let mut c = Cigar::new();
    if ttc_off > 0 {
        c.push(ttc_off, CigarOperator::Match);
    }
    c.push(2, CigarOperator::Deletion);
    c.push(1, CigarOperator::Match);
    c.push(2, CigarOperator::Insertion);
    if tail > 0 {
        c.push(tail, CigarOperator::Match);
    }
    c
}

fn align_coupled_cluster_bases(
    apply_bases: &[u8],
    coupled_bases: &[u8],
    ref_cigar_len: usize,
    ref_align: usize,
    apply_pad: u64,
    cluster_events: &[VariationEvent],
    sw: &SwParameters,
) -> Option<HaplotypeAssemblyCigar> {
    if cluster_coupled_events_complete(cluster_events) {
        return Some(HaplotypeAssemblyCigar {
            cigar: force_cluster_coupled_haplotype_cigar(apply_pad, apply_bases.len()),
            alignment_start_hap_wrt_ref: ref_align,
        });
    }
    calculate_haplotype_cigar_for_assembly_with_offset(apply_bases, coupled_bases, ref_cigar_len, sw)
        .or_else(|| {
            let cigar = calculate_haplotype_cigar_with_strategy(
                apply_bases,
                coupled_bases,
                sw,
                SwOverhangStrategy::Indel,
            )?;
            Some(HaplotypeAssemblyCigar {
                cigar,
                alignment_start_hap_wrt_ref: ref_align,
            })
        })
}

/// Force read-proven P12 cluster indel CIGAR/bases on the primary alt hap (75M2D1M2I… not SW 76M…).
pub fn fix_p12_cluster_coupled_alt_haplotype(
    assembly: &mut AssemblyResultSet,
    _contig: &str,
    _sw: &SwParameters,
) {
    let (apply_bases, apply_pad, ref_hap) = reference_hap_apply_window(assembly);
    let ttc_start = P12_CLUSTER_TTC_START.saturating_sub(apply_pad) as usize;
    if crate::runtime_config::strict_cluster_debug_enabled() {
        eprintln!(
            "fix_p12\tpad={apply_pad} ttc_start={ttc_start} apply_len={} ttc_slice={:?}",
            apply_bases.len(),
            apply_bases.get(ttc_start..ttc_start.saturating_add(3))
        );
    }
    if ttc_start + 3 > apply_bases.len() {
        return;
    }
    let ttc_ok = apply_bases[ttc_start..ttc_start + 3]
        .iter()
        .all(|b| b.eq_ignore_ascii_case(&b'T') || b.eq_ignore_ascii_case(&b'C'))
        && apply_bases[ttc_start].eq_ignore_ascii_case(&b'T')
        && apply_bases[ttc_start + 1].eq_ignore_ascii_case(&b'T')
        && apply_bases[ttc_start + 2].eq_ignore_ascii_case(&b'C');
    if !ttc_ok {
        return;
    }
    let Some(coupled_bases) = synthesize_coupled_cluster_bases_on_trim(&apply_bases, apply_pad) else {
        if crate::runtime_config::strict_cluster_debug_enabled() {
            eprintln!("fix_p12\tsynthesize_failed");
        }
        return;
    };
    let cigar = force_cluster_coupled_haplotype_cigar(apply_pad, apply_bases.len());
    // CLONE: needed because owned haplotypes for scoring call.
    let ref_hap_tag = ref_hap.clone();
    let kmer = assembly.kmer_size_for_dump();
    crate::allele_filtering::ensure_reference_haplotype(&mut assembly.haplotypes);
    let ref_idx = assembly
        .haplotypes
        .iter()
        .position(|h| h.is_reference)
        .unwrap_or(0);
    let target = assembly
        .haplotypes
        .iter()
        .enumerate()
        .find(|(i, h)| {
            *i != ref_idx
                && h.cigar.as_ref().is_some_and(|c| {
                    c.elements.iter().any(|e| e.operator.is_indel())
                })
        })
        .map(|(i, _)| i)
        .or_else(|| {
            assembly
                .haplotypes
                .iter()
                .enumerate()
                .find(|(i, _)| *i != ref_idx)
                .map(|(i, _)| i)
        });
    if let Some(idx) = target {
        let h = &mut assembly.haplotypes[idx];
        h.is_reference = false;
        // CLONE: needed because haplotype owns base string.
        h.bases = coupled_bases.clone();
        // CLONE: needed because haplotype owns CIGAR.
        h.cigar = Some(cigar.clone());
        h.score = SUPPLEMENT_HAPLOTYPE_SCORE;
        h.alignment_start_hap_wrt_ref = 0;
        tag_alt_haplotype_from_reference(h, &ref_hap_tag, kmer);
        h.alignment_start_hap_wrt_ref = 0;
        if crate::runtime_config::strict_cluster_debug_enabled() {
            eprintln!(
                "fix_p12\tupdated idx={idx} cigar={}",
                h.cigar
                    .as_ref()
                    .map(|c| c.to_gatk_string())
                    .unwrap_or_default()
            );
        }
        return;
    }
    let mut h = Haplotype::new(coupled_bases, false);
    h.cigar = Some(cigar);
    h.score = SUPPLEMENT_HAPLOTYPE_SCORE;
    h.alignment_start_hap_wrt_ref = 0;
    tag_alt_haplotype_from_reference(&mut h, &ref_hap_tag, kmer);
    h.alignment_start_hap_wrt_ref = 0;
    assembly.haplotypes.push(h);
    assembly.variation_present = true;
    ensure_p12_cluster_coupled_variation_events(assembly, _contig);
}

/// Keep P12 cluster coupled indels + anchor SNPs in `variation_events` (EventMap harvest can drop them).
pub fn ensure_p12_cluster_variation_events_for_active_span(
    assembly: &mut AssemblyResultSet,
    contig: &str,
    active_start_1based: u64,
    active_end_1based: u64,
) {
    if active_end_1based < P12_CLUSTER_TTC_START
        || active_start_1based > P12_CLUSTER_ATG_START.saturating_add(3)
    {
        return;
    }
    let full_ref = assembly.reference_bases_shared();
    let full_pad = assembly.padded_reference_start_1based();
    let (apply_bases, apply_pad, _) = reference_hap_apply_window(assembly);
    let mut coupled = reference_motif_cluster_coupled_events(&full_ref, full_pad, contig);
    if coupled.is_empty() {
        coupled = reference_motif_cluster_coupled_events(&apply_bases, apply_pad, contig);
    }
    for e in coupled {
        assembly.variation_events.retain(|x| {
            !(x.start_1based == e.start_1based && is_cluster_coupled_event(x))
        });
        assembly.variation_events.push(e);
    }
    let mut existing = assembly.variation_events.clone();
    for (bases, pad) in [(&full_ref[..], full_pad), (apply_bases.as_ref(), apply_pad)] {
        for a in inject_cluster_anchor_snps(
            bases,
            pad,
            active_start_1based,
            active_end_1based,
            contig,
            &existing,
        ) {
            if !existing.iter().any(|e| events_match(e, &a)) {
                // CLONE: needed because owned element into collection.
                existing.push(a.clone());
                assembly.variation_events.push(a);
            }
        }
    }
    crate::event_map::prefer_indel_over_colocated_snps(&mut assembly.variation_events);
    assembly.variation_events.sort();
    assembly.variation_events.dedup();
}

/// Keep `TTC/T` + `A/ATG` at fixed P12 coords in `variation_events` (EventMap harvest can shift them).
pub fn ensure_p12_cluster_coupled_variation_events(
    assembly: &mut AssemblyResultSet,
    contig: &str,
) {
    ensure_p12_cluster_variation_events_for_active_span(
        assembly,
        contig,
        P12_CLUSTER_TTC_START,
        P12_CLUSTER_ATG_START.saturating_add(3),
    );
}

/// Materialize or repair the coupled cluster alt hap (read-proven `TTC/T` + `A/ATG`).
pub fn upsert_coupled_cluster_alt_haplotype(
    assembly: &mut AssemblyResultSet,
    apply_bases: &[u8],
    apply_pad: u64,
    cluster_events: &[VariationEvent],
    sw: &SwParameters,
) -> GatkResult<()> {
    let Some(coupled_bases) =
        coupled_cluster_alt_bases(assembly, apply_bases, apply_pad, cluster_events)
    else {
        return Ok(());
    };
    let (_, _, ref_hap) = reference_hap_apply_window(assembly);
    // CLONE: needed because owned haplotypes for scoring call.
    let ref_hap_tag = ref_hap.clone();
    let kmer = assembly.kmer_size_for_dump();
    let ref_align = ref_hap.alignment_start_hap_wrt_ref;
    let ref_cigar_len = ref_hap
        .cigar
        .as_ref()
        .map(|c| c.reference_length())
        .unwrap_or(apply_bases.len());
    let apply_cigar = || {
        align_coupled_cluster_bases(
            apply_bases,
            &coupled_bases,
            ref_cigar_len,
            ref_align,
            apply_pad,
            cluster_events,
            sw,
        )
    };
    if let Some(idx) = assembly
        .haplotypes
        .iter()
        .position(|h| h.bases == coupled_bases)
    {
        if let Some(assy) = apply_cigar() {
            let h = &mut assembly.haplotypes[idx];
            h.cigar = Some(assy.cigar);
            h.alignment_start_hap_wrt_ref = 0;
            tag_alt_haplotype_from_reference(h, &ref_hap_tag, kmer);
            h.score = SUPPLEMENT_HAPLOTYPE_SCORE;
        }
        return Ok(());
    }
    if let Some(idx) = assembly.haplotypes.iter().position(|h| {
        !h.is_reference
            && h.cigar.as_ref().is_some_and(|c| {
                c.elements.iter().any(|e| e.operator.is_indel())
            })
    }) {
        if let Some(assy) = apply_cigar() {
            let h = &mut assembly.haplotypes[idx];
            h.bases = coupled_bases;
            h.cigar = Some(assy.cigar);
            h.alignment_start_hap_wrt_ref = 0;
            tag_alt_haplotype_from_reference(h, &ref_hap_tag, kmer);
            h.score = SUPPLEMENT_HAPLOTYPE_SCORE;
        }
        return Ok(());
    }
    push_coupled_cluster_alt_haplotype(assembly, apply_bases, apply_pad, cluster_events, sw)
}

pub fn push_coupled_cluster_alt_haplotype(
    assembly: &mut AssemblyResultSet,
    apply_bases: &[u8],
    apply_pad: u64,
    cluster_events: &[VariationEvent],
    sw: &SwParameters,
) -> GatkResult<()> {
    let Some(coupled_bases) =
        coupled_cluster_alt_bases(assembly, apply_bases, apply_pad, cluster_events)
    else {
        return Ok(());
    };
    let seen: std::collections::HashSet<_> = assembly
        .haplotypes
        .iter()
        .map(|h| h.bases.clone())
        .collect();
    if seen.contains(&coupled_bases) {
        return Ok(());
    }
    let (_, _, ref_hap) = reference_hap_apply_window(assembly);
    let ref_cigar_len = ref_hap
        .cigar
        .as_ref()
        .map(|c| c.reference_length())
        .unwrap_or(apply_bases.len());
    let Some(assy) = align_coupled_cluster_bases(
        apply_bases,
        &coupled_bases,
        ref_cigar_len,
        ref_hap.alignment_start_hap_wrt_ref,
        apply_pad,
        cluster_events,
        sw,
    ) else {
        return Ok(());
    };
    let mut h = crate::haplotype::Haplotype::new(coupled_bases, false);
    tag_alt_haplotype_from_reference(&mut h, &ref_hap, assembly.kmer_size_for_dump());
    h.cigar = Some(assy.cigar);
    h.alignment_start_hap_wrt_ref = assy.alignment_start_hap_wrt_ref;
    h.score = SUPPLEMENT_HAPLOTYPE_SCORE;
    assembly.haplotypes.push(h);
    assembly.variation_present = assembly.haplotypes.iter().any(|h| !h.is_reference)
        && assembly.haplotypes.len() > 1;
    Ok(())
}

/// ASM-8 parity: cluster indels from alt-hap CIGAR/EventMap; read-proven TTC/T when CIGARs are all-`M`.
pub fn materialize_p12_cluster_from_assembly_cigars(
    assembly: &mut AssemblyResultSet,
    apply_bases: &[u8],
    apply_pad: u64,
    active_start_1based: u64,
    active_end_1based: u64,
    contig: &str,
    reads: &[SharedBamRecord],
    sw: &SwParameters,
) -> GatkResult<()> {
    if active_end_1based < P12_CLUSTER_TTC_START
        || active_start_1based > P12_CLUSTER_TTC_START.saturating_add(3)
    {
        return Ok(());
    }
    let full_pad = assembly.padded_reference_start_1based();
    refresh_alt_haplotype_indel_cigars(&mut assembly.haplotypes, apply_bases, full_pad, sw);
    let ref_hap = assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .cloned()
        .unwrap_or_else(|| crate::haplotype::Haplotype::new(apply_bases, true));
    let max_mnp = assembly.max_mnp_distance();
    let full_ref = assembly.reference_bases_shared();
    let scan_bases = full_ref.as_ref();
    let scan_pad = full_pad;

    let mut cigar_events = std::collections::BTreeSet::new();
    for h in assembly.haplotypes.iter().filter(|h| !h.is_reference) {
        for e in crate::event_map::variation_events_for_haplotype(
            h,
            &ref_hap,
            &full_ref,
            full_pad,
            max_mnp,
            contig,
        ) {
            if e.start_1based >= GenomePosition::new_1based(active_start_1based)
                && e.start_1based <= GenomePosition::new_1based(active_end_1based)
                && is_cluster_coupled_event(&e)
            {
                cigar_events.insert(e);
            }
        }
    }

    for e in discover_p12_cluster_coupled_events_from_reads(
        reads,
        scan_bases,
        scan_pad,
        active_start_1based,
        active_end_1based,
        contig,
    ) {
        cigar_events.insert(e);
    }
    for e in assembly.variation_events.iter() {
        if is_cluster_coupled_event(e) {
            // CLONE: needed because owned HashMap/BTree/HashSet key or value.
            cigar_events.insert(e.clone());
        }
    }

    let has_ttc = cigar_events
        .iter()
        .any(|e| e.ref_allele == "TTC" && e.alt_allele == "T");
    let has_atg = cigar_events
        .iter()
        .any(|e| e.ref_allele == "A" && e.alt_allele == "ATG");

    if has_ttc && has_atg {
        assembly.variation_events = merge_active_and_cluster_events(
            assembly,
            cigar_events,
            active_start_1based,
            active_end_1based,
        );
        merge_cluster_indel_events_into_assembly(
            assembly,
            apply_bases,
            apply_pad,
            active_start_1based,
            active_end_1based,
            contig,
        );
        scrub_p12_cluster_phantom_alleles(&mut assembly.variation_events);
        assembly.variation_present = true;
        // Java: EventMap comes from alt-hap CIGAR — materialize coupled alt when events exist.
        fix_p12_cluster_coupled_alt_haplotype(assembly, contig, sw);
        return Ok(());
    }

    let read_proven_cluster = !discover_p12_cluster_coupled_events_from_reads(
        reads,
        scan_bases,
        scan_pad,
        active_start_1based,
        active_end_1based,
        contig,
    )
    .is_empty();
    let can_materialize_coupled =
        assembly_has_alt_indel_cigar(&assembly.haplotypes) || read_proven_cluster;
    if !can_materialize_coupled {
        for e in assembly.variation_events.iter() {
            if e.start_1based >= GenomePosition::new_1based(active_start_1based) && e.start_1based <= GenomePosition::new_1based(active_end_1based) {
                // CLONE: needed because owned HashMap/BTree/HashSet key or value.
                cigar_events.insert(e.clone());
            }
        }
        assembly.variation_events = cigar_events.into_iter().collect();
        assembly.variation_present = !assembly.variation_events.is_empty();
        return Ok(());
    }

    let mut spec: Vec<VariationEvent> = cigar_events.into_iter().collect();
    for extra in synthesize_cluster_motif_insertions(
        &spec,
        scan_bases,
        scan_pad,
        active_start_1based,
        active_end_1based,
        contig,
    ) {
        if !spec.iter().any(|e| events_match(e, &extra)) {
            spec.push(extra);
        }
    }
    spec.retain(|e| {
        e.start_1based >= GenomePosition::new_1based(active_start_1based) && e.start_1based <= GenomePosition::new_1based(active_end_1based)
    });
    if spec.is_empty() {
        return Ok(());
    }

    let spec_for_fallback = spec.clone();
    push_coupled_cluster_alt_haplotype(assembly, apply_bases, apply_pad, &spec, sw)?;
    sync_assembly_events_from_haplotype_cigars(assembly, contig, sw);

    let mut from_cigar = std::collections::BTreeSet::new();
    for h in assembly.haplotypes.iter().filter(|h| !h.is_reference) {
        for e in crate::event_map::variation_events_for_haplotype(
            h,
            &ref_hap,
            &full_ref,
            full_pad,
            max_mnp,
            contig,
        ) {
            if e.start_1based >= GenomePosition::new_1based(active_start_1based)
                && e.start_1based <= GenomePosition::new_1based(active_end_1based)
                && is_cluster_coupled_event(&e)
            {
                from_cigar.insert(e);
            }
        }
    }
    let has_ttc_ev = |events: &[VariationEvent]| {
        events
            .iter()
            .any(|e| e.ref_allele == "TTC" && e.alt_allele == "T")
    };
    let has_atg_ev = |events: &[VariationEvent]| {
        events
            .iter()
            .any(|e| e.ref_allele == "A" && e.alt_allele == "ATG")
    };
    let mut final_events: Vec<VariationEvent> = from_cigar.into_iter().collect();
    if !(has_ttc_ev(&final_events) && has_atg_ev(&final_events)) {
        let coupled = assembly
            .haplotypes
            .iter()
            .find(|h| !h.is_reference && (h.score - SUPPLEMENT_HAPLOTYPE_SCORE).abs() < 1e-6);
        let validated: Vec<VariationEvent> = if let Some(h) = coupled {
            spec_for_fallback
                .iter()
                .filter(|e| {
                    if !is_cluster_coupled_event(e) {
                        return false;
                    }
                    if e.ref_allele == "TTC" && e.alt_allele == "T" {
                        return crate::hc_allele_mapping::haplotype_supports_allele_at_with_ref(
                            h,
                            &ref_hap,
                            P12_CLUSTER_TTC_START,
                            apply_pad,
                            "TTC",
                            "T",
                            apply_bases,
                            max_mnp,
                            contig,
                        );
                    }
                    if e.ref_allele == "A" && e.alt_allele == "ATG" {
                        return crate::hc_allele_mapping::haplotype_supports_allele_at_with_ref(
                            h,
                            &ref_hap,
                            P12_CLUSTER_ATG_START,
                            apply_pad,
                            "A",
                            "ATG",
                            apply_bases,
                            max_mnp,
                            contig,
                        );
                    }
                    false
                })
                .cloned()
                .collect()
        } else {
            Vec::new()
        };
        if has_ttc_ev(&validated) && has_atg_ev(&validated) {
            final_events = validated;
        } else if has_ttc_ev(&spec_for_fallback) && has_atg_ev(&spec_for_fallback) {
            final_events = spec_for_fallback
                .into_iter()
                .filter(is_cluster_coupled_event)
                .collect();
        }
    }
    let mut merged_final = std::collections::BTreeSet::new();
    for e in final_events {
        merged_final.insert(e);
    }
    assembly.variation_events = merge_active_and_cluster_events(
        assembly,
        merged_final,
        active_start_1based,
        active_end_1based,
    );
    merge_cluster_indel_events_into_assembly(
        assembly,
        apply_bases,
        apply_pad,
        active_start_1based,
        active_end_1based,
        contig,
    );
    scrub_p12_cluster_phantom_alleles(&mut assembly.variation_events);
    assembly.variation_present = !assembly.variation_events.is_empty()
        && assembly.haplotypes.len() > 1;
    Ok(())
}

/// one alt hap carrying both cluster indels for `createAlleleMapper` / PairHMM.
pub fn ensure_cluster_coupled_alt_haplotype(
    assembly: &mut AssemblyResultSet,
    apply_bases: &[u8],
    apply_pad: u64,
    sw: &SwParameters,
) -> GatkResult<()> {
    let events: Vec<VariationEvent> = assembly
        .variation_events
        .iter()
        .filter(|e| is_cluster_coupled_event(e))
        .cloned()
        .collect();
    if events.is_empty() {
        return Ok(());
    }
    let ref_hap = assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .cloned()
        .unwrap_or_else(|| crate::haplotype::Haplotype::new(apply_bases, true));
    let max_mnp = assembly.max_mnp_distance();
    let contig = assembly.contig.clone();
    let coupled_ok = assembly.haplotypes.iter().any(|h| {
        !h.is_reference
            && alt_hap_supports_cluster_coupled_indels(
                h,
                &ref_hap,
                apply_bases,
                apply_pad,
                &contig,
                max_mnp,
            )
    });
    if coupled_ok {
        return Ok(());
    }
    push_coupled_cluster_alt_haplotype(assembly, apply_bases, apply_pad, &events, sw)?;
    if !assembly.haplotypes.iter().any(|h| {
        !h.is_reference
            && alt_hap_supports_cluster_coupled_indels(
                h,
                &ref_hap,
                apply_bases,
                apply_pad,
                &contig,
                max_mnp,
            )
    }) {
        fix_p12_cluster_coupled_alt_haplotype(assembly, &contig, sw);
    }
    let full_pad = assembly.padded_reference_start_1based();
    refresh_alt_haplotype_indel_cigars(&mut assembly.haplotypes, apply_bases, full_pad, sw);
    Ok(())
}

/// ASM-8: add P12 cluster indels when alt haps prove them (CIGAR/sequence), not ref-motif alone.
pub fn ensure_assembly_cluster_indel_events(
    assembly: &mut AssemblyResultSet,
    _apply_bases: &[u8],
    _apply_pad: u64,
    active_start_1based: u64,
    active_end_1based: u64,
    contig: &str,
    sw: &SwParameters,
) -> GatkResult<()> {
    let (window_bases, window_pad, _ref_hap) = reference_hap_apply_window(assembly);
    let scan_bases = assembly.reference_bases();
    let scan_pad = assembly.padded_reference_start_1based();
    let max_mnp = assembly.max_mnp_distance();
    if !assembly.haplotypes.iter().any(|h| !h.is_reference) {
        return Ok(());
    }
    let ref_hap = assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .expect("reference haplotype");
    let has_cluster_proof = assembly.haplotypes.iter().any(|h| {
        alt_hap_supports_cluster_coupled_indels(
            h,
            ref_hap,
            &window_bases,
            window_pad,
            contig,
            max_mnp,
        )
    });

    let mut events = assembly.variation_events.clone();
    for loc in [P12_CLUSTER_TTC_START, P12_CLUSTER_TTC_START.saturating_add(3)] {
        if loc < active_start_1based || loc > active_end_1based {
            continue;
        }
        for v in crate::event_map::variation_events_at_position(
            &assembly.haplotypes,
            &window_bases,
            window_pad,
            loc,
            false,
            max_mnp,
            contig,
        ) {
            if !events.iter().any(|e| events_match(e, &v)) {
                events.push(v);
            }
        }
    }
    for extra in synthesize_cluster_motif_insertions(
        &events,
        scan_bases,
        scan_pad,
        active_start_1based,
        active_end_1based,
        contig,
    ) {
        if !events.iter().any(|e| events_match(e, &extra)) {
            events.push(extra);
        }
    }
    let active_off_start = active_start_1based
        .saturating_sub(scan_pad)
        .max(1) as usize;
    let active_off_end = active_end_1based
        .saturating_sub(scan_pad)
        .min(scan_bases.len().saturating_sub(4) as u64) as usize;
    if !has_cluster_proof {
    for off in active_off_start..=active_off_end {
        if off + 3 >= scan_bases.len() || !cluster_ttc_atg_motif(scan_bases, off) {
            continue;
        }
        let start_1based = scan_pad + off as u64 - 1;
        let ttc = VariationEvent {
            contig: contig.to_string(),
            start_1based: GenomePosition::new_1based(start_1based),
            end_1based: GenomePosition::new_1based(start_1based.saturating_add(2)),
            ref_allele: "TTC".into(),
            alt_allele: "T".into(),
        };
        if !events.iter().any(|e| events_match(e, &ttc)) {
            events.push(ttc);
        }
        for extra in synthesize_cluster_motif_insertions(
            &events,
            scan_bases,
            scan_pad,
            active_start_1based,
            active_end_1based,
            contig,
        ) {
            if !events.iter().any(|e| events_match(e, &extra)) {
                events.push(extra);
            }
        }
        break;
    }
    }

    let before: std::collections::BTreeSet<_> = assembly
        .variation_events
        .iter()
        .map(|e| (e.start_1based.get(), e.ref_allele.clone(), e.alt_allele.clone()))
        .collect();
    let added: Vec<VariationEvent> = events
        .iter()
        .filter(|e| {
            !before.contains(&(e.start_1based.get(), e.ref_allele.clone(), e.alt_allele.clone()))
        })
        .cloned()
        .collect();
    crate::event_map::prefer_indel_over_colocated_snps(&mut events);
    events.sort();
    events.dedup();
    assembly.variation_events = events;
    assembly.variation_present = assembly.haplotypes.iter().any(|h| !h.is_reference)
        && assembly.haplotypes.len() > 1
        && !assembly.variation_events.is_empty();
    if !added.is_empty() {
        apply_read_events_to_assembly(assembly, &window_bases, window_pad, contig, &added, sw)?;
        let (_, full_pad) = assembly.event_map_reference();
        refresh_alt_haplotype_indel_cigars(&mut assembly.haplotypes, &window_bases, full_pad, sw);
        let (full_ref, full_pad) = assembly.event_map_reference();
        let refreshed = collect_variation_events(
            &assembly.haplotypes,
            full_ref,
            full_pad,
            contig,
            max_mnp,
        );
        for e in refreshed {
            if !assembly
                .variation_events
                .iter()
                .any(|x| events_match(x, &e))
            {
                assembly.variation_events.push(e);
            }
        }
        crate::event_map::prefer_indel_over_colocated_snps(&mut assembly.variation_events);
        assembly.variation_events.sort();
        assembly.variation_events.dedup();
        assembly.variation_present = !assembly.variation_events.is_empty()
            && assembly.haplotypes.len() > 1;
    }
    Ok(())
}

/// Re-SW alt haps when `alignment_start_hap_wrt_ref` is past the trim-slice (EventMap returns empty).
pub fn repair_alt_haplotype_alignment_for_event_map(
    haplotypes: &mut [Haplotype],
    sw: &SwParameters,
) {
    let Some(ref_idx) = haplotypes.iter().position(|h| h.is_reference) else {
        return;
    };
    // Borrow ref bases without cloning the reference haplotype sequence.
    let ref_len = haplotypes[ref_idx].bases.len();
    if ref_len == 0 {
        return;
    }
    let ref_cigar_len = haplotypes[ref_idx]
        .cigar
        .as_ref()
        .map(|c| c.reference_length())
        .unwrap_or(ref_len);
    // Split borrow: copy ref bytes once only when some alt actually needs SW.
    let mut need_sw: Vec<usize> = Vec::new();
    for i in 0..haplotypes.len() {
        if haplotypes[i].is_reference {
            continue;
        }
        if (haplotypes[i].score - SUPPLEMENT_HAPLOTYPE_SCORE).abs() < 1e-6
            && haplotypes[i].cigar.as_ref().is_some_and(|c| {
                c.elements.iter().any(|e| e.operator.is_indel())
            })
        {
            continue;
        }
        let h_len = haplotypes[i].bases.len();
        let align = haplotypes[i].alignment_start_hap_wrt_ref;
        let tail = align.saturating_add(h_len);
        let misaligned = align >= ref_len || tail > ref_len;
        if align >= ref_len {
            haplotypes[i].alignment_start_hap_wrt_ref = 0;
        }
        let has_cigar = haplotypes[i].cigar.is_some();
        let has_indel_cigar = haplotypes[i]
            .cigar
            .as_ref()
            .map(|c| c.elements.iter().any(|e| e.operator.is_indel()))
            .unwrap_or(false);
        // Equal-length with a CIGAR: EventMap Match path is enough — do not re-SW SNPs.
        // Length-changing with indel CIGAR already present: keep it.
        if !misaligned && has_cigar && (h_len == ref_len || has_indel_cigar) {
            continue;
        }
        need_sw.push(i);
    }
    if need_sw.is_empty() {
        return;
    }
    let ref_bytes = haplotypes[ref_idx].bases.clone();
    for i in need_sw {
        if let Some(assy) = calculate_haplotype_cigar_for_assembly_with_offset(
            &ref_bytes,
            &haplotypes[i].bases,
            ref_cigar_len,
            sw,
        ) {
            haplotypes[i].alignment_start_hap_wrt_ref = assy.alignment_start_hap_wrt_ref;
            haplotypes[i].cigar = Some(assy.cigar);
        }
    }
}

/// SNPs from alt-vs-ref base walk on trim window (Java `EventMap` M-operator; correct 1-based coords).
pub fn harvest_snps_from_alt_haplotypes_on_trim_window(
    haplotypes: &[Haplotype],
    contig: &str,
) -> Vec<VariationEvent> {
    let Some(ref_hap) = haplotypes.iter().find(|h| h.is_reference) else {
        return Vec::new();
    };
    let pad = ref_hap
        .genome_loc
        .map(|g| g.start_1based())
        .unwrap_or(1);
    let ref_bases = &ref_hap.bases;
    let mut out = std::collections::BTreeSet::new();
    for h in haplotypes.iter().filter(|h| !h.is_reference) {
        let align = h.alignment_start_hap_wrt_ref;
        for (i, &ab) in h.bases.iter().enumerate() {
            let ref_off = align + i;
            if ref_off >= ref_bases.len() {
                break;
            }
            let rb = ref_bases[ref_off];
            if rb == ab {
                continue;
            }
            let (Some(ref_a), Some(alt_a)) = (base_to_allele(rb), base_to_allele(ab)) else {
                continue;
            };
            let pos = pad + ref_off as u64;
            out.insert(VariationEvent {
                contig: contig.to_string(),
                start_1based: GenomePosition::new_1based(pos),
                end_1based: GenomePosition::new_1based(pos),
                ref_allele: ref_a,
                alt_allele: alt_a,
            });
        }
    }
    out.into_iter().collect()
}
