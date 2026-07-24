pub fn prune_spillover_supplement_haplotypes(assembly: &mut AssemblyResultSet) {
    let ref_len = assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .map(|h| h.bases.len())
        .unwrap_or(0);
    if ref_len == 0 {
        return;
    }
    // L11: allow long insertions up to EventMap allele cap (was ref+8, which pruned
    // +36 INS haps → empty-mapper rescue + prefer_dominant dropped nearby true indels).
    let max_alt_len = ref_len.saturating_add(MAX_VARIATION_EVENT_ALLELE_LENGTH);
    assembly
        .haplotypes
        .retain(|h| h.is_reference || h.bases.len() <= max_alt_len);
}

fn apply_read_events_to_assembly(
    assembly: &mut AssemblyResultSet,
    _ref_bases: &[u8],
    _pad_start: u64,
    contig: &str,
    read_events: &[VariationEvent],
    sw: &SwParameters,
) -> GatkResult<()> {
    if read_events.is_empty() {
        return Ok(());
    }
    let (apply_bases, apply_pad, ref_hap) = reference_hap_apply_window(assembly);
    let ref_cigar_len = ref_hap
        .cigar
        .as_ref()
        .map(|c| c.reference_length())
        .unwrap_or(apply_bases.len());

    let (full_ref, full_pad) = assembly.event_map_reference();
    let full_ref_vec = assembly.reference_bases_shared();
    let full_pad_genomic = assembly.padded_reference_start_1based();
    let asm_events = collect_variation_events(
        &assembly.haplotypes,
        full_ref,
        full_pad,
        contig,
        assembly.max_mnp_distance(),
    );

    let mut seen: std::collections::HashSet<Vec<u8>> = assembly
        .haplotypes
        .iter()
        .map(|h| h.bases.clone())
        .collect();

    let kmer = assembly.kmer_size_for_dump();
    let mut push_alt =
        |ref_for_cigar: &[u8],
         pad_for_cigar: u64,
         event: Option<&VariationEvent>,
         alt_bases: Vec<u8>|
         -> GatkResult<()> {
            // CLONE: needed because owned HashMap/BTree/HashSet key or value.
            if !seen.insert(alt_bases.clone()) {
                return Ok(());
            }
            let cigar_len = if ref_for_cigar.len() == apply_bases.len() {
                ref_cigar_len
            } else {
                ref_for_cigar.len()
            };
            // L11: long insertions in STR must use the deterministic single-event CIGAR.
            // SW often emits a frankenstein CIGAR → EventMap floods (~189 events) and
            // prefer_dominant_spanning_indels drops nearby true indels (15002023 C>CAT).
            let force_single = event.is_some_and(|ev| {
                ev.alt_allele.len() > ev.ref_allele.len()
                    && ev.indel_span().is_long_insertion_span()
            });
            let mut cigar = if force_single {
                event.and_then(|ev| {
                    cigar_for_single_indel_event(ref_for_cigar, pad_for_cigar, ev)
                })
            } else {
                None
            };
            if cigar.is_none() {
                cigar = calculate_haplotype_cigar_with_strategy(
                    ref_for_cigar,
                    &alt_bases,
                    sw,
                    SwOverhangStrategy::Indel,
                )
                .or_else(|| {
                    calculate_haplotype_cigar_for_assembly(
                        ref_for_cigar,
                        &alt_bases,
                        cigar_len,
                        sw,
                    )
                });
            }
            let cigar_has_indel = cigar
                .as_ref()
                .is_some_and(|c| c.elements.iter().any(|e| e.operator.is_indel()));
            if !cigar_has_indel {
                if let Some(ev) = event {
                    cigar = cigar_for_single_indel_event(ref_for_cigar, pad_for_cigar, ev);
                }
            }
            let Some(cigar) = cigar else {
                return Ok(());
            };
            let mut h = Haplotype::new(alt_bases, false);
            tag_alt_haplotype_from_reference(&mut h, &ref_hap, kmer);
            h.cigar = Some(cigar);
            h.score = SUPPLEMENT_HAPLOTYPE_SCORE;
            assembly.haplotypes.push(h);
            Ok(())
        };

    // L11: never chain a long insertion with other spine events — SW on the
    // frankenstein haplotype explodes EventMap and drops nearby true indels.
    let has_long_ins = read_events.iter().any(|e| {
        e.alt_allele.len() > e.ref_allele.len() && e.indel_span().is_long_insertion_span()
    });
    let chained = if has_long_ins {
        None
    } else {
        apply_events_to_ref_chained(&apply_bases, read_events, apply_pad)
    };
    if let Some(cluster_hap) = chained {
        push_alt(&apply_bases, apply_pad, None, cluster_hap)?;
    } else {
        for event in read_events {
            if let Some(alt_bases) = apply_event_to_ref(&apply_bases, event, apply_pad) {
                push_alt(&apply_bases, apply_pad, Some(event), alt_bases)?;
            } else if let Some(alt_bases) =
                apply_event_to_ref(&full_ref_vec, event, full_pad_genomic)
            {
                push_alt(&full_ref_vec, full_pad_genomic, Some(event), alt_bases)?;
            }
        }
    }

    let asm_indels_for_restore: Vec<VariationEvent> = asm_events
        .iter()
        .filter(|e| e.is_indel() && !homopolymer_motif_phantom(e))
        .cloned()
        .collect();
    let mut events: Vec<VariationEvent> = asm_events
        .into_iter()
        .filter(|e| !homopolymer_motif_phantom(e))
        .collect();
    for event in read_events {
        if !events.iter().any(|e| events_match(e, event)) {
            // CLONE: needed because owned element into collection.
            events.push(event.clone());
        }
    }
    crate::event_map::prefer_indel_over_colocated_snps(&mut events);
    crate::event_map::prefer_dominant_spanning_indels(&mut events);
    // L11: prefer_dominant can drop distant true indels when the EventMap flood around a
    // newly injected long-INS includes long DELs. Restore asm indels far from new long-INS.
    const RESTORE_FAR_FROM_LONG_INS_BP: u64 = IndelSpan::FRAGMENT_WINDOW_BP;
    let new_long_ins: Vec<(u64, IndelSpan)> = read_events
        .iter()
        .filter(|e| e.alt_allele.len() > e.ref_allele.len() && e.indel_span().is_long_insertion_span())
        .map(|e| (e.start_1based.get(), e.indel_span()))
        .collect();
    if !new_long_ins.is_empty() {
        for e in &asm_indels_for_restore {
            if events.iter().any(|x| events_match(x, e)) {
                continue;
            }
            let p = e.start_1based.get();
            let near = new_long_ins.iter().any(|(lp, lspan)| {
                p.abs_diff(*lp) <= RESTORE_FAR_FROM_LONG_INS_BP.max(lspan.get() as u64)
            });
            if !near {
                // CLONE: needed because owned element into collection.
                events.push(e.clone());
            }
        }
    }
    events.sort();
    events.dedup();
    assembly.variation_events = events;
    assembly.variation_present = assembly.haplotypes.iter().any(|h| !h.is_reference)
        && assembly.haplotypes.len() > 1;
    prune_spillover_supplement_haplotypes(assembly);
    Ok(())
}
