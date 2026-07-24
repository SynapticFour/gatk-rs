/// L14-C1: motif / TTCT / cluster indel discovery helpers.
/// Behavior-neutral extract from `read_event_discovery/mod.rs` for N-3.

fn synthesize_cluster_motif_insertions(
    discovered: &[VariationEvent],
    ref_bases: &[u8],
    pad_start_1based: u64,
    active_start_1based: u64,
    active_end_1based: u64,
    contig: &str,
) -> Vec<VariationEvent> {
    let existing: BTreeSet<(u64, String, String)> = discovered
        .iter()
        .map(|e| (e.start_1based.get(), e.ref_allele.clone(), e.alt_allele.clone()))
        .collect();

    let mut out = Vec::new();
    for ttc in discovered.iter().filter(|e| e.ref_allele == "TTC" && e.alt_allele == "T") {
        let off = ttc.start_1based.get().saturating_add(1).saturating_sub(pad_start_1based) as usize;
        if !cluster_ttc_atg_motif(ref_bases, off) {
            continue;
        }
        // Coupled site: A/ATG starts 3 bp after TTC/T start (92307324 → 92307327).
        let pos = ttc.start_1based.get().saturating_add(3);
        if pos < active_start_1based
            || pos > active_end_1based
            || ttc.start_1based.get().abs_diff(pos) > CLUSTER_MOTIF_MAX_DISTANCE_BP
        {
            continue;
        }
        // Avoid allocating HashSet lookup keys on the miss-common path.
        if existing
            .iter()
            .any(|(p, r, a)| *p == pos && r == "A" && a == "ATG")
        {
            continue;
        }
        let event = VariationEvent {
            contig: contig.to_string(),
            start_1based: GenomePosition::new_1based(pos),
            end_1based: GenomePosition::new_1based(pos),
            ref_allele: "A".to_string(),
            alt_allele: "ATG".to_string(),
        };
        if allele_len_ok(&event) {
            out.push(event);
        }
    }
    out
}

/// Deletions implied by reads that span a ref segment without aligning to it (no `D` in CIGAR).
fn discover_gap_deletion_events_from_reads(
    reads: &[bam::Record],
    ref_bases: &[u8],
    pad_start_1based: u64,
    active_start_1based: u64,
    active_end_1based: u64,
    contig: &str,
) -> Vec<(u32, VariationEvent)> {
    let pad_start0 = pad_start_1based.saturating_sub(1) as i64;
    let mut support: BTreeMap<(u64, String, String), u32> = BTreeMap::new();

    // Scan only the active span (not full assembly padding) — O(padding×reads) was ~80× slower.
    let active_off_start = active_start_1based
        .saturating_sub(pad_start_1based)
        .max(1) as usize;
    let active_off_end = active_end_1based
        .saturating_sub(pad_start_1based)
        .min(ref_bases.len().saturating_sub(1) as u64) as usize;

    for del_len in 1..=MAX_VARIATION_EVENT_ALLELE_LENGTH.saturating_sub(1) {
        if active_off_end + del_len >= ref_bases.len() {
            continue;
        }
        for off in active_off_start..=active_off_end.min(ref_bases.len().saturating_sub(del_len)) {
            let anchor = ref_bases[off - 1];
            if !is_regular_base(anchor) {
                continue;
            }
            let deleted = &ref_bases[off..off + del_len];
            if !deleted.iter().all(|&b| is_regular_base(b)) {
                continue;
            }
            let mut ref_allele = vec![anchor];
            ref_allele.extend_from_slice(deleted);
            if ref_allele.len() > MAX_VARIATION_EVENT_ALLELE_LENGTH {
                continue;
            }
            let pos = pad_start_1based + off as u64 - 1;
            if !in_active_span(pos, active_start_1based, active_end_1based) {
                continue;
            }
            let ref_allele_s = allele_bytes_to_string(ref_allele);
            let alt_allele_s = allele_bytes_to_string(vec![anchor]);
            let left0 = pad_start0 + off as i64 - 1;
            let right0 = pad_start0 + off as i64 + del_len as i64;

            for rec in reads {
                if rec.is_unmapped() || rec.tid() < 0 {
                    continue;
                }
                let cigar = CigarString(rec.cigar().iter().copied().collect());
                let start0 = rec.pos();
                let end0 = alignment_end0(rec);
                if end0 < left0 || start0 > right0 {
                    continue;
                }
                let seq = rec.seq().as_bytes();
                let Some(ql) = query_index_at_reference_position(start0, &cigar, left0) else {
                    continue;
                };
                let Some(qr) = query_index_at_reference_position(start0, &cigar, right0) else {
                    continue;
                };
                if !seq.get(ql).copied().unwrap_or(0).eq_ignore_ascii_case(&anchor)
                {
                    continue;
                }
                if !seq.get(qr).copied().unwrap_or(0).eq_ignore_ascii_case(&ref_bases[off + del_len])
                {
                    continue;
                }
                let mut gap = true;
                for k in 0..del_len {
                    let ref0 = pad_start0 + off as i64 + k as i64;
                    if query_index_at_reference_position(start0, &cigar, ref0).is_some() {
                        gap = false;
                        break;
                    }
                }
                if gap {
                    *support
                        // CLONE: needed because owned HashMap entry key.
                        .entry((pos, ref_allele_s.clone(), alt_allele_s.clone()))
                        .or_default() += 1;
                }
            }
        }
    }

    let mut scored = Vec::new();
    for ((pos, ref_allele, alt_allele), count) in support {
        if count < MIN_GAP_DELETION_READ_SUPPORT {
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

fn alignment_end0(rec: &bam::Record) -> i64 {
    let mut end = rec.pos();
    for c in rec.cigar().iter() {
        match c {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) | Cigar::Del(n) | Cigar::RefSkip(n) => {
                end += *n as i64;
            }
            _ => {}
        }
    }
    end
}

fn is_regular_base(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')
}

/// GATK `EventMap` merges adjacent mismatches into MNPs/indels from haplotype CIGAR; when assembly
/// is ref-only we approximate deletions from adjacent read-pileup SNPs (e.g. TTC/T at 92307324).
fn collapse_snps_to_deletions(
    snps: &mut Vec<(u32, VariationEvent)>,
    ref_bases: &[u8],
    pad_start_1based: u64,
    active_start_1based: u64,
    active_end_1based: u64,
    contig: &str,
) -> Vec<(u32, VariationEvent)> {
    let mut by_off: BTreeMap<usize, (u32, VariationEvent)> = BTreeMap::new();
    for (count, event) in snps.drain(..) {
        let off = event.start_1based.get().saturating_sub(pad_start_1based) as usize;
        by_off.insert(off, (count, event));
    }
    let mut collapsed = Vec::new();
    let mut consumed: BTreeSet<usize> = BTreeSet::new();

    // Only 2-bp deletions from exactly two adjacent pileup SNPs (P12: TTC/T at 92307324).
    // Longer spans created spurious TTCA/T-style calls when many SNPs cluster.
    const COLLAPSE_DEL_LEN: usize = 2;
    let del_len = COLLAPSE_DEL_LEN;
    {
        for off in 1..ref_bases.len().saturating_sub(del_len) {
            if !(off..off + del_len).all(|o| by_off.contains_key(&o)) {
                continue;
            }
            if (off..off + del_len).any(|o| consumed.contains(&o)) {
                continue;
            }
            let anchor = ref_bases[off - 1];
            if !is_regular_base(anchor) {
                continue;
            }
            let deleted = &ref_bases[off..off + del_len];
            if !deleted.iter().all(|&b| is_regular_base(b)) {
                continue;
            }
            let mut ref_allele = vec![anchor];
            ref_allele.extend_from_slice(deleted);
            if ref_allele.len() > MAX_VARIATION_EVENT_ALLELE_LENGTH {
                continue;
            }
            let pos = pad_start_1based + off as u64 - 1;
            if !in_active_span(pos, active_start_1based, active_end_1based) {
                continue;
            }
            let support: u32 = (off..off + del_len).map(|o| by_off[&o].0).sum();
            let ref_allele_s = allele_bytes_to_string(ref_allele);
            let alt_allele_s = allele_bytes_to_string(vec![anchor]);
            let event = VariationEvent {
                contig: contig.to_string(),
                start_1based: GenomePosition::new_1based(pos),
                end_1based: GenomePosition::new_1based(pos + ref_allele_s.len().saturating_sub(1) as u64),
                ref_allele: ref_allele_s,
                alt_allele: alt_allele_s,
            };
            if allele_len_ok(&event) && event.ref_allele != event.alt_allele {
                collapsed.push((support, event));
                consumed.extend(off..off + del_len);
            }
        }
    }

    snps.extend(
        by_off
            .into_iter()
            .filter(|(o, _)| !consumed.contains(o))
            .map(|(_, v)| v),
    );
    collapsed
}

/// P12 `TTC/T`: anchor `T` + deleted `TC` (2 bp) when ref slice is `TTC`.
fn discover_ttct_deletions_from_reads(
    reads: &[bam::Record],
    ref_bases: &[u8],
    pad_start_1based: u64,
    active_start_1based: u64,
    active_end_1based: u64,
    contig: &str,
) -> Vec<(u32, VariationEvent)> {
    let mut out = Vec::new();
    let active_off_start = active_start_1based
        .saturating_sub(pad_start_1based)
        .max(1) as usize;
    let active_off_end = active_end_1based
        .saturating_sub(pad_start_1based)
        .min(ref_bases.len().saturating_sub(3) as u64) as usize;
    for off in active_off_start..=active_off_end {
        if off + 1 >= ref_bases.len() {
            continue;
        }
        if !cluster_ttc_atg_motif(ref_bases, off) {
            continue;
        }
        let pos = pad_start_1based + off as u64 - 1;
        if !in_active_span(pos, active_start_1based, active_end_1based) {
            continue;
        }
        let event = VariationEvent {
            contig: contig.to_string(),
            start_1based: GenomePosition::new_1based(pos),
            end_1based: GenomePosition::new_1based(pos + 2),
            ref_allele: "TTC".to_string(),
            alt_allele: "T".to_string(),
        };
        if !allele_len_ok(&event) {
            continue;
        }
        if ttct_deletion_read_support(reads, ref_bases, pad_start_1based, off) {
            out.push((1, event));
        }
    }
    out
}

/// P12 NA12878 cluster anchor (hs37d5 20k).
pub const P12_CLUSTER_CORE_START: u64 = 92307229;
pub const P12_CLUSTER_CORE_END: u64 = 92307422;
pub const P12_CLUSTER_TTC_START: u64 = 92307324;
pub const P12_CLUSTER_TG_SNP_START: u64 = 92307333;
pub const P12_CLUSTER_CTC_START: u64 = 92307359;
pub const P12_CLUSTER_TC_SNP_START: u64 = 92307364;
pub const P12_CLUSTER_AC_SNP_START: u64 = 92307383;
/// Upstream of cluster core: hom-alt `PL=130,9,0` block (92305716–92305728).
pub const P12_CLUSTER_UPSTREAM_START: u64 = 92305716;
pub const P12_CLUSTER_UPSTREAM_END: u64 = 92305728;
/// Mid-B sparse soft-clip PairHMM band — aliases [`crate::java_hc_site_semantics`] (Sprint I-1).
pub const P12_SPARSE_SOFTCLIP_PAIRHMM_START: u64 =
    crate::java_hc_site_semantics::SPARSE_SOFTCLIP_PAIRHMM_START;
pub const P12_SPARSE_SOFTCLIP_PAIRHMM_END: u64 =
    crate::java_hc_site_semantics::SPARSE_SOFTCLIP_PAIRHMM_END;

/// P12 cluster: `TTC` followed by `AT` (coupled `TTC/T` + `A/ATG`).
fn cluster_ttc_atg_motif(ref_bases: &[u8], off: usize) -> bool {
    if off < 1 || off + 3 >= ref_bases.len() {
        return false;
    }
    ref_bases[off - 1].eq_ignore_ascii_case(&b'T')
        && ref_bases[off].eq_ignore_ascii_case(&b'T')
        && ref_bases[off + 1].eq_ignore_ascii_case(&b'C')
        && ref_bases[off + 2].eq_ignore_ascii_case(&b'A')
        && ref_bases[off + 3].eq_ignore_ascii_case(&b'T')
}

/// P12 cluster: `CT` → `C` (92307359).
fn cluster_ctc_deletion_motif(ref_bases: &[u8], off: usize) -> bool {
    if off < 1 || off + 1 >= ref_bases.len() {
        return false;
    }
    ref_bases[off - 1].eq_ignore_ascii_case(&b'C')
        && ref_bases[off].eq_ignore_ascii_case(&b'T')
}

/// At least one read shows a gap across `ref[off..off+2]` with anchor `ref[off-1]`.
pub(crate) fn ttct_deletion_read_support(
    reads: &[bam::Record],
    ref_bases: &[u8],
    pad_start_1based: u64,
    off: usize,
) -> bool {
    reads
        .iter()
        .any(|rec| read_supports_ttct_deletion_at_off(rec, ref_bases, pad_start_1based, off))
}

fn read_supports_ttct_deletion_at_off(
    rec: &bam::Record,
    ref_bases: &[u8],
    pad_start_1based: u64,
    off: usize,
) -> bool {
    let pad_start0 = pad_start_1based.saturating_sub(1) as i64;
    let anchor = ref_bases[off - 1];
    let del_len = 2usize;
    if off + del_len >= ref_bases.len() {
        return false;
    }
    let left0 = pad_start0 + off as i64 - 1;
    let right0 = pad_start0 + off as i64 + del_len as i64;
    if rec.is_unmapped() || rec.tid() < 0 {
        return false;
    }
    let cigar = CigarString(rec.cigar().iter().copied().collect());
    let start0 = rec.pos();
    let end0 = alignment_end0(rec);
    if end0 < left0 || start0 > right0 {
        return false;
    }
    let seq = rec.seq().as_bytes();
    let Some(ql) = query_index_at_reference_position(start0, &cigar, left0) else {
        return false;
    };
    let Some(qr) = query_index_at_reference_position(start0, &cigar, right0) else {
        return false;
    };
    if !seq.get(ql).copied().unwrap_or(0).eq_ignore_ascii_case(&anchor) {
        return false;
    }
    if !seq.get(qr).copied().unwrap_or(0).eq_ignore_ascii_case(&ref_bases[off + del_len])
    {
        return false;
    }
    for k in 0..del_len {
        let ref0 = pad_start0 + off as i64 + k as i64;
        if query_index_at_reference_position(start0, &cigar, ref0).is_some() {
            return false;
        }
    }
    true
}

fn read_supports_motif_insertion_at_off(
    rec: &bam::Record,
    ref_bases: &[u8],
    pad_start_1based: u64,
    off: usize,
    ins_len: usize,
) -> bool {
    if off + 1 + ins_len >= ref_bases.len() {
        return false;
    }
    let pad_start0 = pad_start_1based.saturating_sub(1) as i64;
    let anchor = ref_bases[off];
    let after = ref_bases[off + 1];
    if !is_regular_base(anchor) || !is_regular_base(after) {
        return false;
    }
    let left0 = pad_start0 + off as i64;
    if rec.is_unmapped() || rec.tid() < 0 {
        return false;
    }
    let cigar = CigarString(rec.cigar().iter().copied().collect());
    let start0 = rec.pos();
    let seq = rec.seq().as_bytes();
    let Some(ql) = query_index_at_reference_position(start0, &cigar, left0) else {
        return false;
    };
    if !seq.get(ql).copied().unwrap_or(0).eq_ignore_ascii_case(&anchor) {
        return false;
    }
    let Some(inserted) = query_subseq(&seq, ql + 1, ins_len) else {
        return false;
    };
    if !seq
        .get(ql + 1 + ins_len)
        .copied()
        .unwrap_or(0).eq_ignore_ascii_case(&after)
    {
        return false;
    }
    inserted.iter().all(|&b| is_regular_base(b))
}

/// QNAMEs of untrimmed pileup reads with BAM support for P12 cluster coupled indels.
pub fn p12_cluster_coupled_indel_supporting_read_qnames(
    reads: &[bam::Record],
    event: &VariationEvent,
    ref_bases: &[u8],
    pad_start_1based: u64,
) -> std::collections::BTreeSet<Vec<u8>> {
    let mut out = std::collections::BTreeSet::new();
    if is_cluster_coupled_indel(event) && event.ref_allele == "TTC" && event.alt_allele == "T" {
        let off = event
            .start_1based
            .get()
            .saturating_add(1)
            .saturating_sub(pad_start_1based) as usize;
        if off >= 1 && off + 2 < ref_bases.len() {
            for rec in reads {
                if read_supports_ttct_deletion_at_off(rec, ref_bases, pad_start_1based, off) {
                    out.insert(rec.qname().to_owned());
                }
            }
        }
    } else if is_cluster_coupled_indel(event) && event.ref_allele == "A" && event.alt_allele == "ATG"
    {
        let off = event.start_1based.get().saturating_sub(pad_start_1based) as usize;
        let ins_len = event.alt_allele.len().saturating_sub(event.ref_allele.len());
        if ins_len > 0 && off + 1 + ins_len < ref_bases.len() {
            for rec in reads {
                if read_supports_motif_insertion_at_off(rec, ref_bases, pad_start_1based, off, ins_len)
                {
                    out.insert(rec.qname().to_owned());
                }
            }
        }
    } else if is_cluster_ctc_del(event) {
        let off = event
            .start_1based
            .get()
            .saturating_add(1)
            .saturating_sub(pad_start_1based) as usize;
        if off >= 1 && off + 2 < ref_bases.len() {
            for rec in reads {
                if read_supports_ttct_deletion_at_off(rec, ref_bases, pad_start_1based, off) {
                    out.insert(rec.qname().to_owned());
                }
            }
        }
    }
    out
}

/// P12 cluster downstream SNPs (92307403–92307422) when ref anchor matches.
pub const P12_CLUSTER_DOWNSTREAM_SNPS: &[(u64, &str, &str)] = &[
    (92307403, "C", "A"),
    (92307418, "T", "A"),
    (92307420, "T", "G"),
    (92307421, "C", "G"),
    (92307422, "T", "C"),
];

/// P12 cluster SNPs (`92307364 T/C`, `92307383 A/C`, downstream tail) when ref anchor matches.
pub fn inject_cluster_anchor_snps(
    ref_bases: &[u8],
    pad_start_1based: u64,
    active_start_1based: u64,
    active_end_1based: u64,
    contig: &str,
    existing: &[VariationEvent],
) -> Vec<VariationEvent> {
    let mut anchors: Vec<(u64, &str, &str)> = vec![
        (P12_CLUSTER_TG_SNP_START, "T", "G"),
        (P12_CLUSTER_TC_SNP_START, "T", "C"),
        (P12_CLUSTER_AC_SNP_START, "A", "C"),
    ];
    anchors.extend_from_slice(P12_CLUSTER_DOWNSTREAM_SNPS);
    let mut out = Vec::new();
    for (pos, ref_a, alt_a) in anchors {
        if pos < active_start_1based || pos > active_end_1based {
            continue;
        }
        let off = pos.saturating_sub(pad_start_1based) as usize;
        if off >= ref_bases.len() {
            continue;
        }
        let Some(rb) = base_to_allele(ref_bases[off]) else {
            continue;
        };
        if rb != ref_a {
            continue;
        }
        let event = VariationEvent {
            contig: contig.to_string(),
            start_1based: GenomePosition::new_1based(pos),
            end_1based: GenomePosition::new_1based(pos),
            ref_allele: ref_a.to_string(),
            alt_allele: alt_a.to_string(),
        };
        if allele_len_ok(&event) && !existing.iter().any(|e| events_match(e, &event)) {
            out.push(event);
        }
    }
    out
}

/// Inject P12 cluster indels from ref motif when assembly CIGARs are all-`M` (no read pileup required).
pub fn inject_reference_cluster_indel_events(
    ref_bases: &[u8],
    pad_start_1based: u64,
    active_start_1based: u64,
    active_end_1based: u64,
    contig: &str,
    existing: &[VariationEvent],
) -> Vec<VariationEvent> {
    let mut out = Vec::new();
    let active_off_start = active_start_1based
        .saturating_sub(pad_start_1based)
        .max(1) as usize;
    let active_off_end = active_end_1based
        .saturating_sub(pad_start_1based)
        .min(ref_bases.len().saturating_sub(3) as u64) as usize;
    let mut best_ttc: Option<(u64, VariationEvent)> = None;
    for off in active_off_start..=active_off_end {
        if off + 1 >= ref_bases.len() || !cluster_ttc_atg_motif(ref_bases, off) {
            continue;
        }
        let event = VariationEvent {
            contig: contig.to_string(),
            start_1based: GenomePosition::new_1based(pad_start_1based + off as u64 - 1),
            end_1based: GenomePosition::new_1based(pad_start_1based + off as u64 + 1),
            ref_allele: "TTC".to_string(),
            alt_allele: "T".to_string(),
        };
        if !allele_len_ok(&event) {
            continue;
        }
        let dist = event.start_1based.get().abs_diff(P12_CLUSTER_TTC_START);
        if best_ttc.as_ref().is_none_or(|(d, _)| dist < *d) {
            best_ttc = Some((dist, event));
        }
    }
    if let Some((_, event)) = best_ttc {
        if !existing.iter().any(|e| events_match(e, &event)) {
            out.push(event);
        }
    }
    let mut best_ctc: Option<(u64, VariationEvent)> = None;
    for off in active_off_start..=active_off_end {
        if off + 1 >= ref_bases.len() || !cluster_ctc_deletion_motif(ref_bases, off) {
            continue;
        }
        let event = VariationEvent {
            contig: contig.to_string(),
            start_1based: GenomePosition::new_1based(pad_start_1based + off as u64 - 1),
            end_1based: GenomePosition::new_1based(pad_start_1based + off as u64),
            ref_allele: "CT".to_string(),
            alt_allele: "C".to_string(),
        };
        if !allele_len_ok(&event) {
            continue;
        }
        let dist = event.start_1based.get().abs_diff(P12_CLUSTER_CTC_START);
        if best_ctc.as_ref().is_none_or(|(d, _)| dist < *d) {
            best_ctc = Some((dist, event));
        }
    }
    if let Some((_, event)) = best_ctc {
        if !existing.iter().any(|e| events_match(e, &event)) {
            out.push(event);
        }
    }
    for event in inject_cluster_anchor_snps(
        ref_bases,
        pad_start_1based,
        active_start_1based,
        active_end_1based,
        contig,
        existing,
    ) {
        if !out.iter().any(|e| events_match(e, &event)) {
            out.push(event);
        }
    }
    let mut with_ttc: Vec<VariationEvent> = existing.to_vec();
    with_ttc.extend(out.iter().cloned());
    for event in synthesize_cluster_motif_insertions(
        &with_ttc,
        ref_bases,
        pad_start_1based,
        active_start_1based,
        active_end_1based,
        contig,
    ) {
        if !existing.iter().any(|e| events_match(e, &event))
            && !out.iter().any(|e| events_match(e, &event))
        {
            out.push(event);
        }
    }
    out
}

fn merge_read_events(
    mut snps: Vec<(u32, VariationEvent)>,
    mut indels: Vec<(u32, VariationEvent)>,
    gap_dels: Vec<(u32, VariationEvent)>,
    plug_ins: Vec<(u32, VariationEvent)>,
    ttct_dels: Vec<(u32, VariationEvent)>,
    ref_bases: &[u8],
    pad_start_1based: u64,
    active_start_1based: u64,
    active_end_1based: u64,
    contig: &str,
    max_events: usize,
) -> Vec<(u32, VariationEvent)> {
    let collapsed_dels = collapse_snps_to_deletions(
        &mut snps,
        ref_bases,
        pad_start_1based,
        active_start_1based,
        active_end_1based,
        contig,
    );
    indels.extend(collapsed_dels);
    indels.extend(gap_dels);
    indels.extend(plug_ins);
    indels.extend(ttct_dels);
    let indel_spans: Vec<(u64, u64)> = indels
        .iter()
        .map(|(_, e)| (e.start_1based.get(), e.end_1based.get()))
        .collect();
    snps.retain(|(_, e)| {
        !indel_spans.iter().any(|(s, end)| {
            e.start_1based.get() + SNP_NEAR_INDEL_EXCLUSION_BP >= *s
                && e.start_1based <= GenomePosition::new_1based(end.saturating_add(SNP_NEAR_INDEL_EXCLUSION_BP))
        })
    });
    indels.sort_by(|a, b| b.0.cmp(&a.0).then_with(|| a.1.start_1based.cmp(&b.1.start_1based)));
    snps.sort_by(|a, b| b.0.cmp(&a.0).then_with(|| a.1.start_1based.cmp(&b.1.start_1based)));
    let mut out = indels;
    if out.len() > max_events {
        out.truncate(max_events);
    } else {
        let snp_cap = max_events.saturating_sub(out.len());
        out.extend(snps.into_iter().take(snp_cap));
    }
    out
}
pub fn discover_variation_events_from_reads(
    reads: &[bam::Record],
    ref_bases: &[u8],
    pad_start_1based: u64,
    active_start_1based: u64,
    active_end_1based: u64,
    contig: &str,
) -> Vec<VariationEvent> {
    discover_variation_events_from_reads_with_options(
        reads,
        ref_bases,
        pad_start_1based,
        active_start_1based,
        active_end_1based,
        contig,
        ReadEventDiscoveryOptions::strict(),
    )
}
pub fn discover_variation_events_from_reads_with_options(
    reads: &[bam::Record],
    ref_bases: &[u8],
    pad_start_1based: u64,
    active_start_1based: u64,
    active_end_1based: u64,
    contig: &str,
    opts: ReadEventDiscoveryOptions,
) -> Vec<VariationEvent> {
    let snps = discover_snp_events_from_reads(
        reads,
        ref_bases,
        pad_start_1based,
        active_start_1based,
        active_end_1based,
        contig,
        false,
        opts,
    );
    let indels = discover_indel_events_from_reads(
        reads,
        ref_bases,
        pad_start_1based,
        active_start_1based,
        active_end_1based,
        contig,
    );
    let gap_dels = discover_gap_deletion_events_from_reads(
        reads,
        ref_bases,
        pad_start_1based,
        active_start_1based,
        active_end_1based,
        contig,
    );
    let plug_ins = if opts.include_motif_insertions {
        discover_motif_insertion_events_from_reads(
            reads,
            ref_bases,
            pad_start_1based,
            active_start_1based,
            active_end_1based,
            contig,
        )
    } else {
        Vec::new()
    };
    let ttct = if opts.include_motif_insertions {
        discover_ttct_deletions_from_reads(
            reads,
            ref_bases,
            pad_start_1based,
            active_start_1based,
            active_end_1based,
            contig,
        )
    } else {
        Vec::new()
    };
    merge_read_events(
        snps,
        indels,
        gap_dels,
        plug_ins,
        ttct,
        ref_bases,
        pad_start_1based,
        active_start_1based,
        active_end_1based,
        contig,
        opts.max_events_per_region,
    )
    .into_iter()
    .map(|(_, e)| e)
    .collect()
}
