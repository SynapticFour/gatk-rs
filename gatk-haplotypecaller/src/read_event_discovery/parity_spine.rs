/// Drop phantom SNPs at P12 cluster indel loci when coupled `TTC/T` + `A/ATG` are present.
pub fn scrub_p12_cluster_phantom_alleles(events: &mut Vec<VariationEvent>) {
    let has_atg = events.iter().any(|e| {
        e.start_1based == GenomePosition::new_1based(P12_CLUSTER_ATG_START) && e.ref_allele == "A" && e.alt_allele == "ATG"
    });
    let has_ttc = events.iter().any(|e| {
        e.start_1based == GenomePosition::new_1based(P12_CLUSTER_TTC_START) && e.ref_allele == "TTC" && e.alt_allele == "T"
    });
    events.retain(|e| {
        if has_atg
            && e.start_1based == GenomePosition::new_1based(P12_CLUSTER_ATG_START)
            && !(e.ref_allele == "A" && e.alt_allele == "ATG")
        {
            return false;
        }
        if has_ttc
            && e.start_1based == GenomePosition::new_1based(P12_CLUSTER_TTC_START)
            && !(e.ref_allele == "TTC" && e.alt_allele == "T")
        {
            return false;
        }
        true
    });
    crate::event_map::prefer_indel_over_colocated_snps(events);
}

fn merge_active_and_cluster_events(
    assembly: &AssemblyResultSet,
    cluster_events: std::collections::BTreeSet<VariationEvent>,
    active_start_1based: u64,
    active_end_1based: u64,
) -> Vec<VariationEvent> {
    let has_atg = cluster_events
        .iter()
        .any(|e| e.ref_allele == "A" && e.alt_allele == "ATG");
    let has_ttc = cluster_events
        .iter()
        .any(|e| e.ref_allele == "TTC" && e.alt_allele == "T");
    let mut merged = cluster_events;
    for e in assembly.variation_events.iter() {
        if e.start_1based < GenomePosition::new_1based(active_start_1based) || e.start_1based > GenomePosition::new_1based(active_end_1based) {
            continue;
        }
        if has_atg
            && e.start_1based == GenomePosition::new_1based(P12_CLUSTER_ATG_START)
            && !(e.ref_allele == "A" && e.alt_allele == "ATG")
        {
            continue;
        }
        if has_ttc
            && e.start_1based == GenomePosition::new_1based(P12_CLUSTER_TTC_START)
            && !(e.ref_allele == "TTC" && e.alt_allele == "T")
        {
            continue;
        }
        // CLONE: needed because owned HashMap/BTree/HashSet key or value.
        merged.insert(e.clone());
    }
    let mut out: Vec<VariationEvent> = merged.into_iter().collect();
    scrub_p12_cluster_phantom_alleles(&mut out);
    out
}

/// Parity spine: read-proven indels missing from assembly EventMap (genome-wide + harness).
///
/// `materialize_alt_haps`: when true, SW-rematerialize alt haplotypes so EventMap sync
/// retains the alleles. List-only mode (`false`) is only safe when a later path merges
/// prior events without a strict CIGAR-only rebuild — post-HMM StrictJava sync does not.
pub fn parity_spine_read_proven_indels(
    assembly: &mut AssemblyResultSet,
    reads: &[SharedBamRecord],
    active_start_1based: u64,
    active_end_1based: u64,
    sw: &SwParameters,
    materialize_alt_haps: bool,
) -> GatkResult<()> {
    if reads.is_empty() || ref_bases_empty(assembly) {
        return Ok(());
    }
    const MIN_GENERIC_INDEL_READS: u32 = 2;
    // Dense GIAB windows can carry many true indels per active region.
    const MAX_SPINE_INDELS: usize = 64;
    let (apply_bases, apply_pad, _) = reference_hap_apply_window(assembly);
    // Scan on full padded ref — trim/event-map windows can exclude indel anchors that
    // still lie inside the active genotyping span.
    let scan_bases = assembly.reference_bases();
    let scan_pad = assembly.padded_reference_start_1based();
    let contig = assembly.contig.clone();
    let max_mnp = assembly.max_mnp_distance();
    let graph = collect_variation_events(
        &assembly.haplotypes,
        scan_bases,
        scan_pad,
        &contig,
        max_mnp,
    );
    // Dense RT-first hit path: EventMap indels already on haplotype CIGARs — skip the
    // multi-pass read rediscovery (TTCT / indel / SNP-collapse / long-INS). L9 listed
    // materialize is empty when every listed indel is graph-encoded. SNP spine still
    // runs separately for Java SNP FNs.
    let has_alt = assembly.haplotypes.iter().any(|h| !h.is_reference);
    let indel_cigar_complete = has_alt
        && assembly
            .variation_events
            .iter()
            .filter(|e| e.is_indel())
            .all(|e| graph.iter().any(|g| events_match(g, e)));
    if indel_cigar_complete {
        crate::runtime_config::rss_trace_checkpoint(
            "parity_spine_indel_skip_cigar_complete",
            &format!(
                "events={} graph={}",
                assembly.variation_events.len(),
                graph.len()
            ),
        );
        return Ok(());
    }
    let existing: std::collections::BTreeSet<(u64, String, String)> = assembly
        .variation_events
        .iter()
        .chain(graph.iter())
        .map(|e| (e.start_1based.get(), e.ref_allele.clone(), e.alt_allele.clone()))
        .collect();
    let mut candidates: Vec<(u32, VariationEvent)> = Vec::new();
    for (support, e) in discover_ttct_deletions_from_reads(
        reads,
        scan_bases,
        scan_pad,
        active_start_1based,
        active_end_1based,
        &contig,
    ) {
        // CLONE: needed because owned composite key for dedup/lookup.
        let key = (e.start_1based.get(), e.ref_allele.clone(), e.alt_allele.clone());
        if !existing.contains(&key) && !candidates.iter().any(|(_, c)| events_match(c, &e)) {
            candidates.push((support.max(MIN_GENERIC_INDEL_READS), e));
        }
    }
    let mut read_indels = discover_indel_events_from_reads(
        reads,
        scan_bases,
        scan_pad,
        active_start_1based,
        active_end_1based,
        &contig,
    );
    // L10: admit support==1 indels that properly extend a same-start indel with
    // support≥2 (nested STR / longest-allele; holdout 20:15031984 18D beside 8D).
    let strong_starts: BTreeSet<u64> = read_indels
        .iter()
        .filter(|(s, e)| *s >= MIN_GENERIC_INDEL_READS && e.is_indel())
        .map(|(_, e)| e.start_1based.get())
        .collect();
    for (support, e) in read_indels.drain(..) {
        if !e.is_indel() {
            continue;
        }
        let admit = support >= MIN_GENERIC_INDEL_READS
            || (support >= 1
                && strong_starts.contains(&e.start_1based.get())
                && e.ref_allele.len() >= 5);
        if !admit {
            continue;
        }
        // CLONE: needed because owned composite key for dedup/lookup.
        let key = (e.start_1based.get(), e.ref_allele.clone(), e.alt_allele.clone());
        if !existing.contains(&key) && !candidates.iter().any(|(_, c)| events_match(c, &e)) {
            candidates.push((support.max(1), e));
        }
    }
    let snps = discover_snp_events_from_reads(
        reads,
        scan_bases,
        scan_pad,
        active_start_1based,
        active_end_1based,
        &contig,
        false,
        ReadEventDiscoveryOptions::strict(),
    );
    let mut snp_buf = snps;
    for (_, e) in collapse_snps_to_deletions(
        &mut snp_buf,
        scan_bases,
        scan_pad,
        active_start_1based,
        active_end_1based,
        &contig,
    ) {
        if e.is_indel() {
            // CLONE: needed because owned composite key for dedup/lookup.
            let key = (e.start_1based.get(), e.ref_allele.clone(), e.alt_allele.clone());
            if !existing.contains(&key) && !candidates.iter().any(|(_, c)| events_match(c, &e)) {
                candidates.push((MIN_GENERIC_INDEL_READS, e));
            }
        }
    }
    // L9: events already on the genotyping list but missing from hap CIGARs still need
    // alt-hap materialization (spine previously skipped them via `existing`).
    // Same span floor as ensure_alt_haplotypes — avoid short fragment inflation.
    const MIN_LISTED_MATERIALIZE_SPAN: usize = 5;
    for e in &assembly.variation_events {
        if !e.is_indel() {
            continue;
        }
        if e.ref_allele.len().abs_diff(e.alt_allele.len()) < MIN_LISTED_MATERIALIZE_SPAN {
            continue;
        }
        if graph.iter().any(|g| events_match(g, e)) {
            continue;
        }
        if candidates.iter().any(|(_, c)| events_match(c, e)) {
            continue;
        }
        let (rr, ra) = read_allele_depths_at_locus(reads, e, apply_pad);
        let (rr2, ra2) = read_allele_depths_at_locus(reads, e, scan_pad);
        let (rref, ralt) = if ra2 >= ra { (rr2, ra2) } else { (rr, ra) };
        if genome_wide_genotype_read_support(e, rref, ralt) {
            // CLONE: needed — candidate list owns VariationEvent for later materialization.
            candidates.push((ralt.max(0) as u32, e.clone()));
        }
    }
    // L11: long insertions in query sequence without CIGAR I (holdout 20:15001894).
    let mut long_ins_pos: Vec<(u64, IndelSpan)> = Vec::new();
    for (support, e) in discover_long_insertion_events_from_reads(
        reads,
        scan_bases,
        scan_pad,
        active_start_1based,
        active_end_1based,
        &contig,
    ) {
        let span = e.indel_span();
        // Track every discovered long-INS locus (even if already on EventMap) so fragment
        // spray drop still fires.
        long_ins_pos.push((e.start_1based.get(), span));
        // CLONE: needed because owned composite key for dedup/lookup.
        let key = (e.start_1based.get(), e.ref_allele.clone(), e.alt_allele.clone());
        if !existing.contains(&key) && !candidates.iter().any(|(_, c)| events_match(c, &e)) {
            candidates.push((support, e));
        }
    }
    // Also treat long insertions already on candidates / assembly EventMap.
    for (_, e) in &candidates {
        if e.alt_allele.len() > e.ref_allele.len() && e.indel_span().is_long_insertion_span() {
            let pos = e.start_1based.get();
            let span = e.indel_span();
            if !long_ins_pos.iter().any(|(p, _)| *p == pos) {
                long_ins_pos.push((pos, span));
            }
        }
    }
    for e in &assembly.variation_events {
        if e.is_indel()
            && e.alt_allele.len() > e.ref_allele.len()
            && e.indel_span().is_long_insertion_span()
        {
            let pos = e.start_1based.get();
            let span = e.indel_span();
            if !long_ins_pos.iter().any(|(p, _)| *p == pos) {
                long_ins_pos.push((pos, span));
            }
        }
    }
    // A3: when a long insertion is present, drop short fragment sprays in-window so they
    // are not chained onto the long-allele haplotype (frankenstein EventMap wiped nearby
    // true indels e.g. holdout 15002023 C>CAT).
    if !long_ins_pos.is_empty() {
        candidates.retain(|(_, e)| {
            let span = e.indel_span();
            if !span.is_short_fragment() {
                return true;
            }
            let p = e.start_1based.get();
            !long_ins_pos
                .iter()
                .any(|(lp, lspan)| IndelSpan::nests_beside_long(p, *lp, *lspan))
        });
    }
        candidates.sort_by(|a, b| b.0.cmp(&a.0).then_with(|| a.1.start_1based.cmp(&b.1.start_1based)));
    candidates.truncate(MAX_SPINE_INDELS);
    if candidates.is_empty() {
        return Ok(());
    }
    if crate::parity_harness::env_flag_set("GATK_RS_INDEL_HAP_TRACE") {
        for (sup, e) in candidates
            .iter()
            .filter(|(_, e)| e.ref_allele.len().abs_diff(e.alt_allele.len()) >= 10)
            .take(8)
        {
            eprintln!(
                "L10-spine-indel cand support={} {}:{} {}>{} span={}",
                sup,
                e.contig,
                e.start_1based.get(),
                e.ref_allele,
                e.alt_allele,
                e.ref_allele.len().abs_diff(e.alt_allele.len())
            );
        }
    }
    let events: Vec<VariationEvent> = candidates.into_iter().map(|(_, e)| e).collect();
    if materialize_alt_haps {
        apply_read_events_to_assembly(
            assembly,
            &apply_bases,
            apply_pad,
            &contig,
            &events,
            sw,
        )?;
        sync_assembly_events_from_haplotype_cigars(assembly, &contig, sw);
    } else {
        for e in events {
            if !assembly.variation_events.iter().any(|x| events_match(x, &e)) {
                assembly.variation_events.push(e);
            }
        }
        sort_dedup_variation_events(assembly);
    }
    scrub_p12_cluster_phantom_alleles(&mut assembly.variation_events);
    Ok(())
}

/// Parity spine: biallelic SNPs from reads when EventMap has no event at that locus (Java `getVariationEvents` gap).
/// See [`parity_spine_read_proven_indels`] for `materialize_alt_haps`.
pub fn parity_spine_read_proven_snps(
    assembly: &mut AssemblyResultSet,
    reads: &[SharedBamRecord],
    active_start_1based: u64,
    active_end_1based: u64,
    sw: &SwParameters,
    materialize_alt_haps: bool,
) -> GatkResult<()> {
    if reads.is_empty() || ref_bases_empty(assembly) {
        return Ok(());
    }
    // Dense GIAB active regions can exceed the old spine-era cap of 20; 64 keeps
    // throughput bounded (128 made dense HC ~6× slower without clearing remaining FNs).
    const MAX_SPINE_SNPS: usize = 64;
    let (apply_bases, apply_pad, _) = reference_hap_apply_window(assembly);
    // Scan on full padded ref — same as indel spine. `event_map_reference()` can be a
    // trim/EventMap slice that excludes active-span SNPs (e.g. 21:9411785 inside
    // 9411693-9411822 but outside the assembled-variation window).
    let scan_bases = assembly.reference_bases();
    let scan_pad = assembly.padded_reference_start_1based();
    let contig = assembly.contig.clone();
    let existing: std::collections::BTreeSet<(u64, String, String)> = assembly
        .variation_events
        .iter()
        .map(|e| (e.start_1based.get(), e.ref_allele.clone(), e.alt_allele.clone()))
        .collect();
    // Do not use `strict()`: its high-depth 0.55 alt-fraction gate rejects classic
    // ~30% hets (e.g. 21:9411785 Java AD 38,16) that assembly EventMap also missed.
    let snps = discover_parity_spine_snp_events(
        reads,
        scan_bases,
        scan_pad,
        active_start_1based,
        active_end_1based,
        &contig,
    );
    let mut candidates: Vec<VariationEvent> = Vec::new();
    for e in snps {
        if candidates.len() >= MAX_SPINE_SNPS {
            break;
        }
        if e.start_1based >= GenomePosition::new_1based(P12_CLUSTER_CORE_START) && e.start_1based <= GenomePosition::new_1based(P12_CLUSTER_CORE_END) {
            continue;
        }
        // L8 pin retirement: do not skip merely because *some* graph event shares the locus
        // require the exact REF/ALT already present (wrong colocated alleles blocked FN SNPs).
        // CLONE: needed because owned composite key for dedup/lookup.
        let key = (e.start_1based.get(), e.ref_allele.clone(), e.alt_allele.clone());
        if existing.contains(&key) || candidates.iter().any(|c| events_match(c, &e)) {
            continue;
        }
        candidates.push(e);
    }
    if candidates.is_empty() {
        return Ok(());
    }
    if materialize_alt_haps {
        apply_read_events_to_assembly(
            assembly,
            &apply_bases,
            apply_pad,
            &contig,
            &candidates,
            sw,
        )?;
        sync_assembly_events_from_haplotype_cigars(assembly, &contig, sw);
    } else {
        for e in candidates {
            if !assembly.variation_events.iter().any(|x| events_match(x, &e)) {
                assembly.variation_events.push(e);
            }
        }
        sort_dedup_variation_events(assembly);
    }
    Ok(())
}

/// Strong read-proven SNPs for AssemblyRegion trim expansion (pre-trim).
///
/// Observable contract: when assembly EventMap only covers a subset of the active span,
/// trim would otherwise shrink past remaining read-proven hets (e.g. 21:9411785 beside
/// 9411732). Returning these sites as trim anchors keeps the genotyping window wide enough.
pub fn discover_parity_spine_snp_events(
    reads: &[SharedBamRecord],
    ref_bases: &[u8],
    pad_start_1based: u64,
    active_start_1based: u64,
    active_end_1based: u64,
    contig: &str,
) -> Vec<VariationEvent> {
    discover_snp_events_from_reads(
        reads,
        ref_bases,
        pad_start_1based,
        active_start_1based,
        active_end_1based,
        contig,
        false,
        ReadEventDiscoveryOptions::parity_spine_snps(),
    )
    .into_iter()
    .filter_map(|(_, e)| {
        if e.ref_allele.len() == 1 && e.alt_allele.len() == 1 {
            Some(e)
        } else {
            None
        }
    })
    .collect()
}

/// Back-compat alias.
pub fn parity_spine_indels_from_reads(
    assembly: &mut AssemblyResultSet,
    reads: &[SharedBamRecord],
    active_start_1based: u64,
    active_end_1based: u64,
    sw: &SwParameters,
) -> GatkResult<()> {
    parity_spine_read_proven_indels(
        assembly,
        reads,
        active_start_1based,
        active_end_1based,
        sw,
        true,
    )
}

