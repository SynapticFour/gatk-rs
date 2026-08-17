/// L13-C2: read-backed SNP/indel/motif discovery helpers.
/// Behavior-neutral extract from `read_event_discovery/mod.rs` for N-3.

/// P12 `92307364 T/C` from reads on the trimmed apply window.
fn discover_cluster_tc_from_reads(
    reads: &[SharedBamRecord],
    apply_bases: &[u8],
    apply_pad: u64,
    active_start: u64,
    active_end: u64,
    contig: &str,
) -> Option<VariationEvent> {
    if P12_CLUSTER_TC_SNP_START < active_start || P12_CLUSTER_TC_SNP_START > active_end {
        return None;
    }
    let off = P12_CLUSTER_TC_SNP_START.saturating_sub(apply_pad) as usize;
    if off >= apply_bases.len() {
        return None;
    }
    let ref_b = apply_bases[off].to_ascii_uppercase();
    if ref_b != b'T' {
        return None;
    }
    let mut alt_count = 0u32;
    for rec in reads {
        if rec.is_unmapped() {
            continue;
        }
        let cigar = CigarString(rec.cigar().iter().copied().collect());
        let start0 = rec.pos();
        let seq = rec.seq().as_bytes();
        let ref_pos0 = apply_pad.saturating_sub(1) as i64 + off as i64;
        let Some(qi) = query_index_at_reference_position(start0, &cigar, ref_pos0) else {
            continue;
        };
        let Some(qb) = seq.get(qi) else {
            continue;
        };
        if qb.eq_ignore_ascii_case(&b'C') {
            alt_count += 1;
        }
    }
    if alt_count == 0 {
        return None;
    }
    Some(VariationEvent {
        contig: contig.to_string(),
        start_1based: GenomePosition::new_1based(P12_CLUSTER_TC_SNP_START),
        end_1based: GenomePosition::new_1based(P12_CLUSTER_TC_SNP_START),
        ref_allele: "T".to_string(),
        alt_allele: "C".to_string(),
    })
}

const MIN_NON_CLUSTER_SNP_DEPTH: u32 = 4;

/// Homopolymer motif insertions from assembly paths (e.g. `T/TGAT`) — not Java cluster calls.
fn homopolymer_motif_phantom(e: &VariationEvent) -> bool {
    e.ref_allele.len() == 1
        && e.alt_allele.len() > e.ref_allele.len()
        && e.alt_allele.starts_with(&e.ref_allele)
        && !(e.ref_allele == "A" && e.alt_allele == "ATG")
}

/// Safe query subsequence; returns `None` if `[start, start+len)` is out of bounds.
fn query_subseq(seq: &[u8], start: usize, len: usize) -> Option<&[u8]> {
    seq.get(start..start.checked_add(len)?)
}

fn allele_len_ok(event: &VariationEvent) -> bool {
    event.ref_allele.len() <= MAX_VARIATION_EVENT_ALLELE_LENGTH
        && event.alt_allele.len() <= MAX_VARIATION_EVENT_ALLELE_LENGTH
}

/// Fast check: reads in the active span carry indels or SNP mismatches (assembly retry trigger).
pub fn reads_support_variation_in_active_span(
    reads: &[SharedBamRecord],
    ref_bases: &[u8],
    pad_start_1based: u64,
    active_start_1based: u64,
    active_end_1based: u64,
) -> bool {
    let pad_start0 = pad_start_1based.saturating_sub(1) as i64;
    for rec in reads {
        if rec.is_unmapped() || rec.tid() < 0 {
            continue;
        }
        let cigar = CigarString(rec.cigar().iter().copied().collect());
        let start0 = rec.pos();
        let end0 = alignment_end0(rec);
        let active0_start = active_start_1based.saturating_sub(1) as i64;
        let active0_end = active_end_1based.saturating_sub(1) as i64;
        if end0 < active0_start || start0 > active0_end {
            continue;
        }
        for op in cigar.iter() {
            match op {
                Cigar::Ins(n) | Cigar::Del(n) if *n > 0 && *n <= 4 => return true,
                _ => {}
            }
        }
        let seq = rec.seq().as_bytes();
        for off in 0..ref_bases.len() {
            let ref_pos0 = pad_start0 + off as i64;
            if ref_pos0 < start0 || ref_pos0 > end0 {
                continue;
            }
            let pos_1 = ref_pos0 as u64 + 1;
            if !in_active_span(pos_1, active_start_1based, active_end_1based) {
                continue;
            }
            let Some(qi) = query_index_at_reference_position(start0, &cigar, ref_pos0) else {
                continue;
            };
            let Some(qb) = seq.get(qi) else {
                continue;
            };
            if !ref_bases[off].eq_ignore_ascii_case(qb) {
                return true;
            }
        }
    }
    false
}

fn discover_snp_events_from_reads(
    reads: &[SharedBamRecord],
    ref_bases: &[u8],
    pad_start_1based: u64,
    active_start_1based: u64,
    active_end_1based: u64,
    contig: &str,
    ref_only_fallback: bool,
    opts: ReadEventDiscoveryOptions,
) -> Vec<(u32, VariationEvent)> {
    let pad_start0 = pad_start_1based.saturating_sub(1) as i64;
    let mut alt_counts: BTreeMap<usize, [u32; 4]> = BTreeMap::new();
    let mut depth: BTreeMap<usize, u32> = BTreeMap::new();
    let idx = |b: u8| -> usize {
        match b.to_ascii_uppercase() {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => 4,
        }
    };

    for rec in reads {
        if rec.is_unmapped() || rec.tid() < 0 {
            continue;
        }
        let cigar = CigarString(rec.cigar().iter().copied().collect());
        let start0 = rec.pos();
        let seq = rec.seq().as_bytes();
        for off in 0..ref_bases.len() {
            let ref_pos0 = pad_start0 + off as i64;
            if ref_pos0 < start0 {
                continue;
            }
            let pos_1 = ref_pos0 as u64 + 1;
            if !in_active_span(pos_1, active_start_1based, active_end_1based) {
                continue;
            }
            let Some(qi) = query_index_at_reference_position(start0, &cigar, ref_pos0) else {
                continue;
            };
            let Some(qb) = seq.get(qi) else {
                continue;
            };
            if idx(*qb) >= 4 {
                continue;
            }
            *depth.entry(off).or_default() += 1;
            let rb = ref_bases[off];
            if !rb.eq_ignore_ascii_case(qb) {
                alt_counts.entry(off).or_default()[idx(*qb)] += 1;
            }
        }
    }

    let mut scored = Vec::new();
    for (off, alt) in alt_counts {
        let d = *depth.get(&off).unwrap_or(&0);
        let min_depth = if ref_only_fallback {
            REF_ONLY_MIN_SNP_DEPTH.max(opts.min_snp_depth)
        } else {
            opts.min_snp_depth
        };
        if d < min_depth {
            continue;
        }
        let rb = ref_bases[off];
        let Some(ref_allele) = base_to_allele(rb) else {
            continue;
        };
        let (best_alt_idx, best_alt_count) = alt
            .iter()
            .enumerate()
            .max_by_key(|(_, c)| *c)
            .map(|(i, c)| (i, *c))
            .unwrap_or((0, 0));
        if best_alt_count < opts.min_snp_alt_reads {
            continue;
        }
        let min_frac = if ref_only_fallback {
            REF_ONLY_MIN_SNP_ALT_FRACTION
        } else if d >= opts.high_depth_threshold {
            opts.high_depth_min_alt_fraction
        } else {
            opts.min_snp_alt_fraction
        };
        if (best_alt_count as f64) / (d as f64) < min_frac {
            continue;
        }
        let alt_byte = match best_alt_idx {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        };
        let Some(alt_allele) = base_to_allele(alt_byte) else {
            continue;
        };
        if ref_allele == alt_allele {
            continue;
        }
        let pos = pad_start_1based + off as u64;
        let event = VariationEvent {
            contig: contig.to_string(),
            start_1based: GenomePosition::new_1based(pos),
            end_1based: GenomePosition::new_1based(pos),
            ref_allele,
            alt_allele,
        };
        if allele_len_ok(&event) {
            scored.push((best_alt_count, event));
        }
    }
    scored
}

/// Indels from read CIGARs (GATK `EventMap` padding semantics).
fn discover_indel_events_from_reads(
    reads: &[SharedBamRecord],
    ref_bases: &[u8],
    pad_start_1based: u64,
    active_start_1based: u64,
    active_end_1based: u64,
    contig: &str,
) -> Vec<(u32, VariationEvent)> {
    let pad_start0 = pad_start_1based.saturating_sub(1) as i64;
    let mut support: std::collections::HashMap<(u64, String, String), u32> =
        std::collections::HashMap::new();

    for rec in reads {
        if rec.is_unmapped() || rec.tid() < 0 {
            continue;
        }
        let cigar = CigarString(rec.cigar().iter().copied().collect());
        let start0 = rec.pos();
        let seq = rec.seq().as_bytes();
        let mut ref_pos0 = start0;
        let mut query_pos: usize = 0;

        for op in cigar.iter() {
            match op {
                Cigar::Del(n) => {
                    let len = *n as usize;
                    let off = (ref_pos0 - pad_start0) as usize;
                    if len > 0
                        && len <= MAX_VARIATION_EVENT_ALLELE_LENGTH
                        && off > 0
                        && off + len <= ref_bases.len()
                    {
                        let anchor = ref_bases[off - 1];
                        if is_regular_base(anchor) {
                            let mut ref_allele = vec![anchor];
                            ref_allele.extend_from_slice(&ref_bases[off..off + len]);
                            if ref_allele.iter().all(|&b| is_regular_base(b)) {
                                let pos = pad_start_1based + off as u64 - 1;
                                if in_active_span(pos, active_start_1based, active_end_1based) {
                                    let ref_allele = allele_bytes_to_string(ref_allele);
                                    let alt_allele =
                                        allele_bytes_to_string(vec![anchor]);
                                    *support
                                        // CLONE: needed because owned HashMap entry key.
                                        .entry((pos, ref_allele.clone(), alt_allele.clone()))
                                        .or_default() += 1;
                                }
                            }
                        }
                    }
                    ref_pos0 += len as i64;
                }
                Cigar::Ins(n) => {
                    let len = *n as usize;
                    // After preceding ref-consuming ops, ref_pos0 is past the anchor; anchor is ref_pos0-1.
                    let ref_index = (ref_pos0 - pad_start0 - 1) as usize;
                    if len > 0
                        && len < MAX_VARIATION_EVENT_ALLELE_LENGTH
                        && ref_index < ref_bases.len()
                    {
                        let anchor = ref_bases[ref_index];
                        if is_regular_base(anchor) {
                            let Some(inserted) = query_subseq(&seq, query_pos, len) else {
                                query_pos += len;
                                continue;
                            };
                            let mut alt_allele = vec![anchor];
                            alt_allele.extend_from_slice(inserted);
                            if alt_allele.iter().all(|&b| is_regular_base(b)) {
                                let pos = pad_start_1based + ref_index as u64;
                                if in_active_span(pos, active_start_1based, active_end_1based) {
                                    let ref_allele = allele_bytes_to_string(vec![anchor]);
                                    let alt_allele = allele_bytes_to_string(alt_allele);
                                    *support
                                        // CLONE: needed because owned HashMap entry key.
                                        .entry((pos, ref_allele.clone(), alt_allele.clone()))
                                        .or_default() += 1;
                                }
                            }
                        }
                    }
                    query_pos += len;
                }
                Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                    ref_pos0 += *n as i64;
                    query_pos += *n as usize;
                }
                Cigar::SoftClip(n) => {
                    query_pos += *n as usize;
                }
                Cigar::HardClip(_) | Cigar::RefSkip(_) | Cigar::Pad(_) => {}
            }
        }
    }

    let mut scored = Vec::new();
    for ((pos, ref_allele, alt_allele), count) in support {
        if count < MIN_INDEL_READ_SUPPORT {
            continue;
        }
        let event = VariationEvent {
            contig: contig.to_string(),
            start_1based: GenomePosition::new_1based(pos),
            end_1based: GenomePosition::new_1based(pos + ref_allele.len().saturating_sub(1) as u64),
            ref_allele,
            alt_allele,
        };
        if allele_len_ok(&event) && event.ref_allele != event.alt_allele {
            scored.push((count, event));
        }
    }
    scored
}

/// Insertions between ref[off] and ref[off+1] when query carries extra bases (no `I` in CIGAR at locus).
/// Test-only: production uses CIGAR indels + [`synthesize_cluster_motif_insertions`] (GATK `EventMap` path).
#[cfg(test)]
fn discover_plug_insertion_events_from_reads(
    reads: &[SharedBamRecord],
    ref_bases: &[u8],
    pad_start_1based: u64,
    active_start_1based: u64,
    active_end_1based: u64,
    contig: &str,
) -> Vec<(u32, VariationEvent)> {
    let pad_start0 = pad_start_1based.saturating_sub(1) as i64;
    let active_off_start = active_start_1based
        .saturating_sub(pad_start_1based) as usize;
    let active_off_end = active_end_1based
        .saturating_sub(pad_start_1based)
        .min(ref_bases.len().saturating_sub(2) as u64) as usize;
    const MIN_PLUG_INSERTION_READ_SUPPORT: u32 = 1;
    let mut support: std::collections::HashMap<(u64, String, String), u32> = std::collections::HashMap::new();

    for ins_len in 1..=MAX_VARIATION_EVENT_ALLELE_LENGTH.saturating_sub(1) {
        for off in active_off_start..=active_off_end {
            if off + 1 + ins_len > ref_bases.len() {
                continue;
            }
            let anchor = ref_bases[off];
            let after = ref_bases[off + 1];
            if !is_regular_base(anchor) || !is_regular_base(after) {
                continue;
            }
            let pos = pad_start_1based + off as u64;
            if !in_active_span(pos, active_start_1based, active_end_1based) {
                continue;
            }
            let left0 = pad_start0 + off as i64;
            let right0 = pad_start0 + off as i64 + 1;

            for rec in reads {
                if rec.is_unmapped() || rec.tid() < 0 {
                    continue;
                }
                let cigar = CigarString(rec.cigar().iter().copied().collect());
                let start0 = rec.pos();
                let end0 = alignment_end0(rec);
                if end0 < left0 || start0 > right0 + ins_len as i64 {
                    continue;
                }
                let seq = rec.seq().as_bytes();
                let Some(ql) = query_index_at_reference_position(start0, &cigar, left0) else {
                    continue;
                };
                let Some(_qr) = query_index_at_reference_position(start0, &cigar, right0) else {
                    continue;
                };
                if seq.get(ql).copied().unwrap_or(0).to_ascii_uppercase()
                    != anchor.to_ascii_uppercase()
                {
                    continue;
                }
                let flank_qi = ql + 1 + ins_len;
                let right_ok = query_index_at_reference_position(start0, &cigar, right0)
                    .and_then(|qr| seq.get(qr).copied())
                    .is_some_and(|b| b.to_ascii_uppercase() == after.to_ascii_uppercase())
                    || seq
                        .get(flank_qi)
                        .copied()
                        .is_some_and(|b| b.to_ascii_uppercase() == after.to_ascii_uppercase());
                if !right_ok {
                    continue;
                }
                let Some(inserted) = query_subseq(&seq, ql + 1, ins_len) else {
                    continue;
                };
                if !inserted.iter().all(|&b| is_regular_base(b)) {
                    continue;
                }
                let mut alt_allele = vec![anchor];
                alt_allele.extend_from_slice(inserted);
                let ref_allele_s = allele_bytes_to_string(vec![anchor]);
                let alt_allele_s = allele_bytes_to_string(alt_allele);
                *support
                    // CLONE: needed because owned HashMap entry key.
                    .entry((pos, ref_allele_s.clone(), alt_allele_s.clone()))
                    .or_default() += 1;
            }
        }
    }

    let mut scored = Vec::new();
    for ((pos, ref_allele, alt_allele), count) in support {
        if count < MIN_PLUG_INSERTION_READ_SUPPORT {
            continue;
        }
        let event = VariationEvent {
            contig: contig.to_string(),
            start_1based: GenomePosition::new_1based(pos),
            end_1based: GenomePosition::new_1based(pos),
            ref_allele,
            alt_allele,
        };
        if allele_len_ok(&event) && event.ref_allele != event.alt_allele {
            scored.push((count, event));
        }
    }
    scored
}

/// Anchor + inserted motif in read (e.g. `A` + `TG` → `ATG` when ref is `AT`).
fn discover_motif_insertion_events_from_reads(
    reads: &[SharedBamRecord],
    ref_bases: &[u8],
    pad_start_1based: u64,
    active_start_1based: u64,
    active_end_1based: u64,
    contig: &str,
) -> Vec<(u32, VariationEvent)> {
    let pad_start0 = pad_start_1based.saturating_sub(1) as i64;
    let active_off_start = active_start_1based
        .saturating_sub(pad_start_1based) as usize;
    let active_off_end = active_end_1based
        .saturating_sub(pad_start_1based)
        .min(ref_bases.len().saturating_sub(3) as u64) as usize;
    const MIN_MOTIF_INSERTION_READ_SUPPORT: u32 = 1;
    let mut support: std::collections::HashMap<(u64, String, String), u32> = std::collections::HashMap::new();

    for ins_len in 1..=MAX_VARIATION_EVENT_ALLELE_LENGTH.saturating_sub(1) {
        for off in active_off_start..=active_off_end {
            if off + 1 + ins_len >= ref_bases.len() {
                continue;
            }
            let anchor = ref_bases[off];
            let after = ref_bases[off + 1];
            if !is_regular_base(anchor) || !is_regular_base(after) {
                continue;
            }
            let pos = pad_start_1based + off as u64;
            if !in_active_span(pos, active_start_1based, active_end_1based) {
                continue;
            }
            let left0 = pad_start0 + off as i64;

            for rec in reads {
                if rec.is_unmapped() || rec.tid() < 0 {
                    continue;
                }
                let cigar = CigarString(rec.cigar().iter().copied().collect());
                let start0 = rec.pos();
                let seq = rec.seq().as_bytes();
                let Some(ql) = query_index_at_reference_position(start0, &cigar, left0) else {
                    continue;
                };
                if !seq.get(ql).copied().unwrap_or(0).eq_ignore_ascii_case(&anchor)
                {
                    continue;
                }
                let Some(inserted) = query_subseq(&seq, ql + 1, ins_len) else {
                    continue;
                };
                if !seq
                    .get(ql + 1 + ins_len)
                    .copied()
                    .unwrap_or(0).eq_ignore_ascii_case(&after)
                {
                    continue;
                }
                if !inserted.iter().all(|&b| is_regular_base(b)) {
                    continue;
                }
                let mut alt_allele = vec![anchor];
                alt_allele.extend_from_slice(inserted);
                let ref_allele_s = allele_bytes_to_string(vec![anchor]);
                let alt_allele_s = allele_bytes_to_string(alt_allele);
                *support
                    // CLONE: needed because owned HashMap entry key.
                    .entry((pos, ref_allele_s.clone(), alt_allele_s.clone()))
                    .or_default() += 1;
            }
        }
    }

    let mut scored = Vec::new();
    for ((pos, ref_allele, alt_allele), count) in support {
        if count < MIN_MOTIF_INSERTION_READ_SUPPORT {
            continue;
        }
        let event = VariationEvent {
            contig: contig.to_string(),
            start_1based: GenomePosition::new_1based(pos),
            end_1based: GenomePosition::new_1based(pos),
            ref_allele,
            alt_allele,
        };
        if allele_len_ok(&event) && event.ref_allele != event.alt_allele {
            scored.push((count, event));
        }
    }
    scored
}

