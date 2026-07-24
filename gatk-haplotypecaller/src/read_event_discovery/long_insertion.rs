/// True when one ACGT base is ≥70% of the insertion (poly-A/T plugs in homopolymer runs).
fn insertion_is_low_complexity(inserted: &[u8]) -> bool {
    if inserted.is_empty() {
        return true;
    }
    let mut counts = [0u32; 4];
    for &b in inserted {
        match b.to_ascii_uppercase() {
            b'A' => counts[0] += 1,
            b'C' => counts[1] += 1,
            b'G' => counts[2] += 1,
            b'T' => counts[3] += 1,
            _ => {}
        }
    }
    let n = inserted.len() as u32;
    let max = counts.into_iter().max().unwrap_or(0);
    max.saturating_mul(10) >= n.saturating_mul(7)
}

/// Long insertions present in read sequence between ref anchors but not as CIGAR `I`
/// (L11 holdout `20:15001894`: Java +36 INS; BAM shows 2D/4D while query carries the motif).
/// Narrower than [`discover_motif_insertion_events_from_reads`] (which is off on `strict`):
/// only `ins_len ∈ [10, MAX)`, support ≥2 — cheaper on M4 and fewer short-motif FPs.
/// STR self-similarity makes naive plug matching explode (same read matches many offsets).
/// Production filters (phenotype, not locus pins):
/// 1. inserted ≠ reference continuation (true alt vs ref-matching read-ahead);
/// 2. `LONG_INS_REMATCH_BP` bases after the plug rematch reference (back on haplotype);
/// 3. left-normalized: last inserted base ≠ anchor (reject tandem-shiftable forms);
/// 4. reject low-complexity (homopolymer-dominated) insertions.
fn discover_long_insertion_events_from_reads(
    reads: &[bam::Record],
    ref_bases: &[u8],
    pad_start_1based: u64,
    active_start_1based: u64,
    active_end_1based: u64,
    contig: &str,
) -> Vec<(u32, VariationEvent)> {
    use crate::read_projection::query_index_at_reference_position;

    const MIN_LONG_INS_LEN: usize = IndelSpan::LONG_INSERTION_MIN.get();
    const MIN_LONG_INS_READ_SUPPORT: u32 = 2;
    const LONG_INS_REMATCH_BP: usize = 10;
    let pad_start0 = pad_start_1based.saturating_sub(1) as i64;
    let active_off_start = active_start_1based
        .saturating_sub(pad_start_1based) as usize;
    let active_off_end = active_end_1based
        .saturating_sub(pad_start_1based)
        .min(ref_bases.len().saturating_sub(3) as u64) as usize;
    let mut support: BTreeMap<(u64, String, String), u32> = BTreeMap::new();

    let max_ins = MAX_VARIATION_EVENT_ALLELE_LENGTH.saturating_sub(1);
    if max_ins < MIN_LONG_INS_LEN {
        return Vec::new();
    }
    for ins_len in MIN_LONG_INS_LEN..=max_ins {
        for off in active_off_start..=active_off_end {
            // Need anchor, after, and rematch bases on reference.
            if off + 2 + LONG_INS_REMATCH_BP > ref_bases.len() {
                continue;
            }
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
            let ref_cont = &ref_bases[off + 1..off + 1 + ins_len];
            let ref_rematch = &ref_bases[off + 2..off + 2 + LONG_INS_REMATCH_BP];

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
                if !seq
                    .get(ql)
                    .copied()
                    .unwrap_or(0).eq_ignore_ascii_case(&anchor)
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
                // Reject low-complexity plugs (poly-A/T etc.) — holdout 15040030 false +20
                // over true CA>C deletion in an A-run.
                if insertion_is_low_complexity(inserted) {
                    continue;
                }
                // Reject ref-identical read-ahead (common in M-aligned STR).
                if inserted.eq_ignore_ascii_case(ref_cont) {
                    continue;
                }
                // Require rematch onto reference after the insertion plug.
                let rematch_start = ql + 1 + ins_len + 1;
                let Some(q_rematch) = query_subseq(&seq, rematch_start, LONG_INS_REMATCH_BP) else {
                    continue;
                };
                if !q_rematch.eq_ignore_ascii_case(ref_rematch) {
                    continue;
                }
                // Left-normalized: last inserted base ≠ anchor (tandem-shiftable otherwise).
                if inserted
                    .last()
                    .copied()
                    .unwrap_or(0).eq_ignore_ascii_case(&anchor)
                {
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
        if count < MIN_LONG_INS_READ_SUPPORT {
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
    // Prefer longest ALT at the same POS (nested STR lengths collapse to Java-scale allele).
    scored.sort_by(|a, b| {
        a.1.start_1based
            .cmp(&b.1.start_1based)
            .then_with(|| b.1.alt_allele.len().cmp(&a.1.alt_allele.len()))
            .then_with(|| b.0.cmp(&a.0))
    });
    let mut best_at_pos: BTreeMap<u64, (u32, VariationEvent)> = BTreeMap::new();
    for (count, event) in scored {
        let pos = event.start_1based.get();
        best_at_pos.entry(pos).or_insert((count, event));
    }
    // Drop satellite long-INS that start inside a longer upstream insertion's span
    // (same misaligned STR haplotype; holdout 15001894 +36 vs 15001902/924 plugs).
    let mut ordered: Vec<(u32, VariationEvent)> = best_at_pos.into_values().collect();
    ordered.sort_by(|a, b| {
        b.1.alt_allele
            .len()
            .cmp(&a.1.alt_allele.len())
            .then_with(|| b.0.cmp(&a.0))
            .then_with(|| a.1.start_1based.cmp(&b.1.start_1based))
    });
    let mut kept: Vec<(u32, VariationEvent)> = Vec::new();
    for (count, event) in ordered {
        let pos = event.start_1based.get();
        let span = event.alt_allele.len().saturating_sub(event.ref_allele.len());
        let satellite = kept.iter().any(|(_, k)| {
            let kpos = k.start_1based.get();
            let kspan = k.alt_allele.len().saturating_sub(k.ref_allele.len());
            kspan >= span
                && kpos < pos
                && pos <= kpos.saturating_add(kspan as u64)
        });
        if !satellite {
            kept.push((count, event));
        }
    }
    kept
}

