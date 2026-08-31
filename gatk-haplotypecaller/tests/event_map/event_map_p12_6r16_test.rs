//! 6R.16 TEST-ONLY: production `trim_to` + post-trim `fix_p12` / `ensure_p12` /
//! `ensure_alt_haplotypes_for_variation_events` on the canonical control hap.
//! Does not call W-H1. Does not change production algorithms.

#[cfg(test)]
mod traces {
    use crate::alignment::SwParameters;
    use crate::assembly_region_iterator::AssemblyRegion;
    use crate::assembly_region_trimmer::{
        AssemblyRegionTrimmer, AssemblyRegionTrimmerConfig, TrimVariant,
    };
    use crate::assembly_result_set::{AssemblyResultSet, DEFAULT_MAX_MNP_DISTANCE};
    use crate::cigar::{Cigar, CigarOperator};
    use crate::engine::preserve_untrimmed_indel_haplotypes;
    use crate::event_map::{variation_events_for_haplotype, VariationEvent};
    use crate::event_map_rebuild::{rebuild_variation_events, RebuildVariationEventsOpts};
    use crate::genome_loc::{GenomeLoc, GenomePosition};
    use crate::haplotype::Haplotype;
    use crate::hc_allele_mapping::create_allele_mapper;
    use crate::read_event_discovery::{
        ensure_alt_haplotypes_for_variation_events, ensure_cluster_coupled_alt_haplotype,
        ensure_p12_cluster_variation_events_for_active_span, fix_p12_cluster_coupled_alt_haplotype,
        materialize_p12_cluster_from_assembly_cigars, prune_spillover_supplement_haplotypes,
        P12_CLUSTER_ATG_START, P12_CLUSTER_TTC_START, SUPPLEMENT_HAPLOTYPE_SCORE,
    };
    use crate::read_threading_assembler::AssemblyStatus;
    use gatk_core::reference::SequenceDictionary;
    use std::path::Path;

    const ALT_WIN: &[u8] = b"CTTTTATGTGATGTAT";
    const REAL_P12_ACTIVE_START: u64 = P12_CLUSTER_TTC_START - 96;
    const REAL_P12_ACTIVE_END: u64 = P12_CLUSTER_TTC_START + 76;

    fn fixture_paths() -> Option<(std::path::PathBuf, std::path::PathBuf)> {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
        let ref_path = root.join("parity/realworld/assets/hs37d5.simple.fa");
        let bam = root.join("parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam");
        if !ref_path.is_file() || !bam.is_file() {
            return None;
        }
        Some((ref_path, bam))
    }

    fn load_real_p12() -> Option<(Vec<u8>, u64, AssemblyRegion, SequenceDictionary)> {
        use crate::assembly_region_finalize::{
            assembly_reference_read, finalize_region_reads_for_assembly,
            gatk_min_tail_quality_for_assembly, records_to_assembly_reads,
        };
        use crate::read_model::ReadFilterParams;
        use crate::walker_apply::{call_disposition, AssemblyRegionCallDisposition};
        use crate::walker_traversal::{
            flatten_assembly_regions, traverse_assembly_region_walker, WalkerTraversalConfig,
        };
        use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache};

        let (ref_path, bam) = fixture_paths()?;
        let dict = SequenceDictionary::from_fasta_path(&ref_path).ok()?;
        let interval = format!("2:{REAL_P12_ACTIVE_START}-{REAL_P12_ACTIVE_END}");
        let specs = parse_intervals_cli_string(&dict, &interval).ok()?;
        let walk = traverse_assembly_region_walker(
            &dict,
            &specs,
            &ref_path,
            &bam,
            &ReadFilterParams::gatk_standard_hc(),
            &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
        )
        .ok()?;
        let regions = flatten_assembly_regions(&walk);
        let region = regions
            .iter()
            .find(|r| {
                matches!(
                    call_disposition(r),
                    AssemblyRegionCallDisposition::ActiveFull
                ) && r.start.get() <= P12_CLUSTER_TTC_START
                    && r.end.get() >= P12_CLUSTER_ATG_START
            })?
            .clone();
        let mut ref_cache = ReferenceWindowCache::new(ref_path.clone(), 4);
        let reference = assembly_reference_read(&dict, &mut ref_cache, &region).ok()?;
        let _finalized = finalize_region_reads_for_assembly(
            &region.reads,
            &region,
            true,
            gatk_min_tail_quality_for_assembly(10),
            false,
        );
        let _reads = records_to_assembly_reads(&_finalized);
        let pad = region
            .extended_start
            .get()
            .saturating_sub(crate::assembly_region_finalize::GATK_REFERENCE_PADDING_FOR_ASSEMBLY)
            .max(1);
        Some((reference.bases, pad, region, dict))
    }

    fn control_alt_bases(ref_bytes: &[u8], pad: u64) -> Option<Vec<u8>> {
        let ttc_off = P12_CLUSTER_TTC_START.saturating_sub(pad) as usize;
        let atg_off = P12_CLUSTER_ATG_START.saturating_sub(pad) as usize;
        if ttc_off + 3 > ref_bytes.len() || atg_off >= ref_bytes.len() {
            return None;
        }
        if &ref_bytes[ttc_off..ttc_off + 3] != b"TTC" {
            return None;
        }
        if ref_bytes[atg_off] != b'A' && ref_bytes[atg_off] != b'a' {
            return None;
        }
        let mut out = ref_bytes.to_vec();
        out.remove(ttc_off + 1);
        out.remove(ttc_off + 1);
        let atg_adj = atg_off.saturating_sub(2);
        if !out
            .get(atg_adj)
            .copied()
            .unwrap_or(0)
            .eq_ignore_ascii_case(&b'A')
        {
            return None;
        }
        out.insert(atg_adj + 1, b'T');
        out.insert(atg_adj + 2, b'G');
        Some(out)
    }

    fn forced_cigar(pad: u64, ref_len: usize) -> Cigar {
        let ttc_off = P12_CLUSTER_TTC_START.saturating_sub(pad) as usize;
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

    fn control_haplotype(ref_bytes: &[u8], pad: u64) -> Option<Haplotype> {
        let bases = control_alt_bases(ref_bytes, pad)?;
        let cigar = forced_cigar(pad, ref_bytes.len());
        let mut h = Haplotype::new(bases, false);
        h.cigar = Some(cigar);
        h.score = SUPPLEMENT_HAPLOTYPE_SCORE;
        h.alignment_start_hap_wrt_ref = 0;
        h.genome_loc = Some(GenomeLoc::new(
            pad,
            pad.saturating_add(ref_bytes.len() as u64).saturating_sub(1),
        ));
        Some(h)
    }

    fn is_regular(b: u8) -> bool {
        matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')
    }

    #[derive(Clone, Debug, PartialEq, Eq)]
    struct JavaVc {
        start: u64,
        end: u64,
        ref_al: String,
        alt_al: String,
    }

    fn java_source_eventmap_proposed(hap: &Haplotype, ref_bytes: &[u8], pad: u64) -> Vec<JavaVc> {
        let Some(cigar) = hap.cigar.as_ref() else {
            return Vec::new();
        };
        let alignment = hap.bases.as_slice();
        let mut ref_pos = hap.alignment_start_hap_wrt_ref;
        let mut alignment_pos = 0usize;
        let mut proposed = Vec::new();
        for (cigar_index, el) in cigar.elements.iter().enumerate() {
            let n = el.length;
            match el.operator {
                CigarOperator::Insertion => {
                    if ref_pos > 0 && ref_pos <= ref_bytes.len() {
                        let insertion_start = pad + ref_pos as u64 - 1;
                        let ref_byte = ref_bytes[ref_pos - 1];
                        let is_edge = cigar_index == 0 || cigar_index == cigar.elements.len() - 1;
                        if is_regular(ref_byte) && !is_edge {
                            let mut ins = vec![ref_byte];
                            let end = alignment_pos.saturating_add(n).min(alignment.len());
                            ins.extend_from_slice(&alignment[alignment_pos..end]);
                            if ins.len() >= 2 && ins.iter().copied().all(is_regular) {
                                proposed.push(JavaVc {
                                    start: insertion_start,
                                    end: insertion_start,
                                    ref_al: String::from_utf8(vec![ref_byte]).unwrap(),
                                    alt_al: String::from_utf8(ins).unwrap(),
                                });
                            }
                        }
                    }
                    alignment_pos += n;
                }
                CigarOperator::Deletion => {
                    if ref_pos > 0 && ref_pos + n <= ref_bytes.len() {
                        let mut del = vec![ref_bytes[ref_pos - 1]];
                        del.extend_from_slice(&ref_bytes[ref_pos..ref_pos + n]);
                        let deletion_start = pad + ref_pos as u64 - 1;
                        if is_regular(ref_bytes[ref_pos - 1]) && del.iter().copied().all(is_regular)
                        {
                            proposed.push(JavaVc {
                                start: deletion_start,
                                end: deletion_start + n as u64,
                                ref_al: String::from_utf8(del).unwrap(),
                                alt_al: String::from_utf8(vec![ref_bytes[ref_pos - 1]]).unwrap(),
                            });
                        }
                    }
                    ref_pos += n;
                }
                CigarOperator::Match => {
                    for offset in 0..n {
                        if ref_pos + offset >= ref_bytes.len() {
                            break;
                        }
                        let rb = ref_bytes[ref_pos + offset];
                        let ab = alignment
                            .get(alignment_pos + offset)
                            .copied()
                            .unwrap_or(b'N');
                        if rb != ab && is_regular(rb) && is_regular(ab) {
                            let start = pad + ref_pos as u64 + offset as u64;
                            proposed.push(JavaVc {
                                start,
                                end: start,
                                ref_al: String::from_utf8(vec![rb]).unwrap(),
                                alt_al: String::from_utf8(vec![ab]).unwrap(),
                            });
                        }
                    }
                    ref_pos += n;
                    alignment_pos += n;
                }
                CigarOperator::SoftClip => {
                    alignment_pos += n;
                }
                _ => {
                    ref_pos += n;
                    alignment_pos += n;
                }
            }
        }
        proposed
    }

    fn java_make_block(vc1: JavaVc, vc2: JavaVc) -> JavaVc {
        let vc1_snp = vc1.ref_al.len() == 1 && vc1.alt_al.len() == 1;
        if vc1_snp {
            if vc1.ref_al == vc2.ref_al {
                JavaVc {
                    start: vc1.start,
                    end: vc1.end,
                    ref_al: vc1.ref_al,
                    alt_al: format!("{}{}", vc1.alt_al, &vc2.alt_al[1..]),
                }
            } else {
                JavaVc {
                    start: vc1.start,
                    end: vc2.end,
                    ref_al: vc2.ref_al,
                    alt_al: vc1.alt_al,
                }
            }
        } else {
            let (ins, del) = if vc1.alt_al.len() > vc1.ref_al.len() {
                (vc1, vc2)
            } else {
                (vc2, vc1)
            };
            JavaVc {
                start: del.start,
                end: del.end,
                ref_al: del.ref_al,
                alt_al: ins.alt_al,
            }
        }
    }

    fn java_source_eventmap_merged(proposed: &[JavaVc]) -> Vec<JavaVc> {
        let mut by_start: std::collections::BTreeMap<u64, JavaVc> =
            std::collections::BTreeMap::new();
        for vc in proposed {
            if let Some(prev) = by_start.remove(&vc.start) {
                by_start.insert(vc.start, java_make_block(prev, vc.clone()));
            } else {
                by_start.insert(vc.start, vc.clone());
            }
        }
        by_start.into_values().collect()
    }

    const P12_FOCUS_LO: u64 = P12_CLUSTER_TTC_START - 8;
    const P12_FOCUS_HI: u64 = P12_CLUSTER_ATG_START + 60;

    fn dump_events(label: &str, events: &[VariationEvent]) {
        eprintln!("EVENTS {label} n={}", events.len());
        for e in events {
            eprintln!(
                "  {}-{} REF={} ALT={} indel={}",
                e.start_1based.get(),
                e.end_1based.get(),
                e.ref_allele,
                e.alt_allele,
                e.is_indel()
            );
        }
    }

    fn dump_events_compact(label: &str, events: &[VariationEvent]) {
        let focus: Vec<_> = events
            .iter()
            .filter(|e| {
                let s = e.start_1based.get();
                s >= P12_FOCUS_LO && s <= P12_FOCUS_HI
            })
            .collect();
        eprintln!(
            "EVENTS {label} n={} p12_focus={}",
            events.len(),
            focus.len()
        );
        for e in &focus {
            eprintln!(
                "  {}-{} REF={} ALT={} indel={}",
                e.start_1based.get(),
                e.end_1based.get(),
                e.ref_allele,
                e.alt_allele,
                e.is_indel()
            );
        }
        if events.len() > focus.len() {
            eprintln!(
                "  ... {} events outside P12 focus window omitted",
                events.len() - focus.len()
            );
        }
    }

    fn dump_java_compact(label: &str, vcs: &[JavaVc]) {
        let focus: Vec<_> = vcs
            .iter()
            .filter(|v| v.start >= P12_FOCUS_LO && v.start <= P12_FOCUS_HI)
            .collect();
        eprintln!("{label} n={} p12_focus={}", vcs.len(), focus.len());
        for v in &focus {
            eprintln!(
                "    {}-{} REF={} ALT={}",
                v.start, v.end, v.ref_al, v.alt_al
            );
        }
        if vcs.len() > focus.len() {
            eprintln!(
                "    ... {} events outside P12 focus window omitted",
                vcs.len() - focus.len()
            );
        }
    }

    fn dump_hap(label: &str, i: usize, h: &Haplotype, full_pad: u64) {
        let ttc_off = P12_CLUSTER_TTC_START
            .saturating_sub(h.genome_loc.map(|g| g.start_1based()).unwrap_or(full_pad))
            as usize;
        let win = h
            .bases
            .get(ttc_off.saturating_sub(4)..ttc_off.saturating_add(12))
            .map(|b| String::from_utf8_lossy(b).into_owned());
        eprintln!(
            "HAP {label} idx={i} len={} cigar={} align={} score={} is_ref={} kmer={} genome={:?} alt_win={} cluster_win={win:?}",
            h.bases.len(),
            h.cigar
                .as_ref()
                .map(|c| c.to_gatk_string())
                .unwrap_or_default(),
            h.alignment_start_hap_wrt_ref,
            h.score,
            h.is_reference,
            h.kmer_size,
            h.genome_loc,
            h.bases.windows(ALT_WIN.len()).any(|w| w == ALT_WIN)
        );
        if h.bases.len() <= 320 {
            eprintln!("  BASES_FULL idx={i} {}", String::from_utf8_lossy(&h.bases));
        }
    }

    fn has_2d1m2i(h: &Haplotype) -> bool {
        let Some(c) = h.cigar.as_ref() else {
            return false;
        };
        let ops: Vec<_> = c.elements.iter().map(|e| (e.operator, e.length)).collect();
        ops.windows(3).any(|w| {
            w[0].0 == CigarOperator::Deletion
                && w[0].1 == 2
                && w[1].0 == CigarOperator::Match
                && w[1].1 == 1
                && w[2].0 == CigarOperator::Insertion
                && w[2].1 == 2
        })
    }

    fn hap_has_ttc(events: &[VariationEvent]) -> bool {
        events
            .iter()
            .any(|e| e.ref_allele == "TTC" && e.alt_allele == "T")
    }
    fn hap_has_atg(events: &[VariationEvent]) -> bool {
        events
            .iter()
            .any(|e| e.ref_allele == "A" && e.alt_allele == "ATG")
    }

    fn stage_table(
        stage: &str,
        assembly: &AssemblyResultSet,
        ref_hap: &Haplotype,
        full_ref: &[u8],
        full_pad: u64,
    ) {
        let mut ttc = false;
        let mut atg = false;
        let mut coupled = false;
        let mut canon = false;
        for h in &assembly.haplotypes {
            if h.is_reference {
                continue;
            }
            let ev = variation_events_for_haplotype(h, ref_hap, full_ref, full_pad, 0, "2");
            let ht = hap_has_ttc(&ev);
            let ha = hap_has_atg(&ev);
            ttc |= ht;
            atg |= ha;
            coupled |= ht && ha;
            canon |= h.bases.windows(ALT_WIN.len()).any(|w| w == ALT_WIN) && has_2d1m2i(h);
        }
        eprintln!(
            "STAGE {stage} hap_count={} idxs=0..{} canon_p12={canon} ttc_hap={ttc} atg_hap={atg} coupled={coupled}",
            assembly.haplotypes.len(),
            assembly.haplotypes.len().saturating_sub(1)
        );
    }

    #[test]
    fn six_r16_trimmed_p12_post_supplement_audit() {
        let Some((ref_bytes, pad, orig_region, dict)) = load_real_p12() else {
            eprintln!("Real-data P12 comparison unavailable");
            return;
        };
        let Some(alt) = control_haplotype(&ref_bytes, pad) else {
            panic!("control haplotype construction failed");
        };
        let mut ref_hap = Haplotype::new(ref_bytes.as_slice(), true);
        let mut rc = Cigar::new();
        rc.push(ref_bytes.len(), CigarOperator::Match);
        ref_hap.cigar = Some(rc);
        ref_hap.genome_loc = alt.genome_loc;

        assert_eq!(
            alt.cigar.as_ref().map(|c| c.to_gatk_string()).as_deref(),
            Some("696M2D1M2I674M")
        );
        assert!(alt.bases.windows(ALT_WIN.len()).any(|w| w == ALT_WIN));

        let untrimmed = AssemblyResultSet::from_assembly_for_calling_owned(
            AssemblyStatus::AssembledSomeVariation,
            85,
            vec![ref_hap.clone(), alt.clone()],
            ref_bytes.clone(),
            pad,
            "2",
            DEFAULT_MAX_MNP_DISTANCE,
        );
        eprintln!(
            "UNTRIMMED orig_region active={}..{} extended={}..{} pad={pad} ref_len={} events={}",
            orig_region.start.get(),
            orig_region.end.get(),
            orig_region.extended_start.get(),
            orig_region.extended_end.get(),
            ref_bytes.len(),
            untrimmed.variation_events().len()
        );
        dump_events("untrimmed_variation_events", untrimmed.variation_events());

        let trim_variants: Vec<TrimVariant> = untrimmed
            .variation_events()
            .iter()
            .map(|e| TrimVariant {
                contig: e.contig.clone(),
                start: e.start_1based.get(),
                end: e.end_1based.get(),
                is_indel: e.is_indel(),
            })
            .collect();
        let trimmer =
            AssemblyRegionTrimmer::new(AssemblyRegionTrimmerConfig::gatk_defaults(), &dict, "2");
        let trim_result = trimmer.trim(&orig_region, &trim_variants, Some(&orig_region.reference));
        eprintln!(
            "TRIM_RESULT present={} variant={:?}..{:?} padded={:?}..{:?}",
            trim_result.variation_present,
            trim_result.variant_start,
            trim_result.variant_end,
            trim_result.padded_variant_start,
            trim_result.padded_variant_end
        );
        let region_for_genotyping = AssemblyRegionTrimmer::apply_trim(&orig_region, &trim_result);
        eprintln!(
            "TRIM_TO_SPAN active={}..{} extended={}..{} (this is the span trim_to uses)",
            region_for_genotyping.start.get(),
            region_for_genotyping.end.get(),
            region_for_genotyping.extended_start.get(),
            region_for_genotyping.extended_end.get()
        );

        let mut assembly = untrimmed.trim_to(&region_for_genotyping).expect("trim_to");
        eprintln!(
            "AFTER_trim_to n_haps={} events={}",
            assembly.haplotypes.len(),
            assembly.variation_events().len()
        );
        dump_events("after_trim_to", assembly.variation_events());
        for (i, h) in assembly.haplotypes.iter().enumerate() {
            dump_hap("after_trim_to", i, h, pad);
        }
        let trimmed_canon = assembly
            .haplotypes
            .iter()
            .find(|h| !h.is_reference)
            .expect("trimmed canonical hap");
        eprintln!(
            "PROD_TRIM_STATE len={} cigar={} align={} (6R.14 active-span expectation was 173/96M2D1M2I74M/600)",
            trimmed_canon.bases.len(),
            trimmed_canon
                .cigar
                .as_ref()
                .map(|c| c.to_gatk_string())
                .unwrap_or_default(),
            trimmed_canon.alignment_start_hap_wrt_ref
        );

        let active_trim = GenomeLoc::new(REAL_P12_ACTIVE_START, REAL_P12_ACTIVE_END);
        if let Some(t) = alt.trim(&active_trim, true) {
            eprintln!(
                "COMPARE_6R14_active_trim len={} cigar={} align={}",
                t.bases.len(),
                t.cigar
                    .as_ref()
                    .map(|c| c.to_gatk_string())
                    .unwrap_or_default(),
                t.alignment_start_hap_wrt_ref
            );
        }

        let sw = SwParameters::gatk_haplotype_to_reference();
        let n_before_preserve = assembly.haplotypes.len();
        preserve_untrimmed_indel_haplotypes(&untrimmed, &mut assembly, &region_for_genotyping, &sw);
        eprintln!(
            "AFTER_preserve n_haps {}→{} (score=1000 indel hap skipped by refresh)",
            n_before_preserve,
            assembly.haplotypes.len()
        );
        for (i, h) in assembly.haplotypes.iter().enumerate() {
            dump_hap("after_preserve", i, h, pad);
        }

        let ref_after_trim = assembly
            .haplotypes
            .iter()
            .find(|h| h.is_reference)
            .cloned()
            .expect("ref after trim");
        stage_table(
            "after_trim_preserve",
            &assembly,
            &ref_after_trim,
            &ref_bytes,
            pad,
        );

        let before_fix: Vec<Haplotype> = assembly.haplotypes.clone();
        let events_before_fix = assembly.variation_events().to_vec();
        dump_events("BEFORE_fix_p12", &events_before_fix);
        for (i, h) in before_fix.iter().enumerate() {
            dump_hap("BEFORE_fix_p12", i, h, pad);
        }

        fix_p12_cluster_coupled_alt_haplotype(&mut assembly, "2", &sw);

        dump_events("AFTER_fix_p12", assembly.variation_events());
        eprintln!(
            "FIX_P12 events_changed={}",
            events_before_fix != assembly.variation_events()
        );
        for (i, h) in assembly.haplotypes.iter().enumerate() {
            dump_hap("AFTER_fix_p12", i, h, pad);
            if i < before_fix.len() {
                let b = &before_fix[i];
                eprintln!(
                    "  DELTA idx={i} bases_eq={} cigar_eq={} align {}→{} genome_eq={} score_eq={} kmer {}→{} is_ref_eq={}",
                    b.bases == h.bases,
                    b.cigar == h.cigar,
                    b.alignment_start_hap_wrt_ref,
                    h.alignment_start_hap_wrt_ref,
                    b.genome_loc == h.genome_loc,
                    (b.score - h.score).abs() < 1e-12,
                    b.kmer_size,
                    h.kmer_size,
                    b.is_reference == h.is_reference
                );
            }
        }
        stage_table("after_fix_p12", &assembly, &ref_after_trim, &ref_bytes, pad);

        let hap_events_before_ensure: Vec<(usize, Vec<VariationEvent>)> = assembly
            .haplotypes
            .iter()
            .enumerate()
            .map(|(i, h)| {
                (
                    i,
                    variation_events_for_haplotype(h, &ref_after_trim, &ref_bytes, pad, 0, "2"),
                )
            })
            .collect();
        for (i, ev) in &hap_events_before_ensure {
            dump_events(&format!("hap{i}_EventMap_BEFORE_ensure"), ev);
        }
        let events_before_ensure = assembly.variation_events().to_vec();
        dump_events("BEFORE_ensure_p12", &events_before_ensure);

        ensure_p12_cluster_variation_events_for_active_span(
            &mut assembly,
            "2",
            orig_region.start.get(),
            orig_region.end.get(),
        );
        dump_events("AFTER_ensure_p12", assembly.variation_events());
        let hap_events_after_ensure: Vec<(usize, Vec<VariationEvent>)> = assembly
            .haplotypes
            .iter()
            .enumerate()
            .map(|(i, h)| {
                (
                    i,
                    variation_events_for_haplotype(h, &ref_after_trim, &ref_bytes, pad, 0, "2"),
                )
            })
            .collect();
        for (i, ev) in &hap_events_after_ensure {
            dump_events(&format!("hap{i}_EventMap_AFTER_ensure"), ev);
            if let Some((_, before)) = hap_events_before_ensure.iter().find(|(j, _)| *j == *i) {
                eprintln!("  hap{i}_EventMap_changed={}", before != ev);
            }
        }
        eprintln!(
            "ENSURE has_TTC={} has_ATG={} has_TTT={} has_TAT={} has_CAT={} has_AG={}",
            assembly
                .variation_events()
                .iter()
                .any(|e| e.start_1based.get() == P12_CLUSTER_TTC_START
                    && e.ref_allele == "TTC"
                    && e.alt_allele == "T"),
            assembly
                .variation_events()
                .iter()
                .any(|e| e.start_1based.get() == P12_CLUSTER_ATG_START
                    && e.ref_allele == "A"
                    && e.alt_allele == "ATG"),
            assembly
                .variation_events()
                .iter()
                .any(|e| e.start_1based.get() == P12_CLUSTER_TTC_START - 1
                    && e.ref_allele == "TTT"
                    && e.alt_allele == "T"),
            assembly
                .variation_events()
                .iter()
                .any(|e| e.start_1based.get() == P12_CLUSTER_TTC_START + 2
                    && e.ref_allele == "C"
                    && e.alt_allele == "TAT"),
            assembly
                .variation_events()
                .iter()
                .any(|e| e.start_1based.get() == P12_CLUSTER_TTC_START + 2
                    && e.ref_allele == "C"
                    && e.alt_allele == "CAT"),
            assembly
                .variation_events()
                .iter()
                .any(|e| e.start_1based.get() == P12_CLUSTER_ATG_START
                    && e.ref_allele == "A"
                    && e.alt_allele == "G"),
        );
        let events_after_ensure_p12 = assembly.variation_events().to_vec();
        let ensure_p12_has_ttc = events_after_ensure_p12
            .iter()
            .any(|e| e.ref_allele == "TTC" && e.alt_allele == "T");
        let ensure_p12_has_atg = events_after_ensure_p12
            .iter()
            .any(|e| e.ref_allele == "A" && e.alt_allele == "ATG");
        stage_table(
            "after_ensure_p12",
            &assembly,
            &ref_after_trim,
            &ref_bytes,
            pad,
        );

        let apply_bases = ref_after_trim.bases.as_slice();
        let apply_pad = ref_after_trim
            .genome_loc
            .map(|g| g.start_1based())
            .unwrap_or(pad);

        let haps_before_alt: Vec<Haplotype> = assembly.haplotypes.clone();
        let n_before_alt = haps_before_alt.len();
        ensure_alt_haplotypes_for_variation_events(&mut assembly, &sw).expect("ensure_alt");
        eprintln!(
            "AFTER_ensure_alt n_haps {}→{}",
            n_before_alt,
            assembly.haplotypes.len()
        );
        dump_events(
            "AFTER_ensure_alt_variation_events",
            assembly.variation_events(),
        );
        eprintln!(
            "ENSURE_ALT clobbered_TTC={} clobbered_ATG={}",
            ensure_p12_has_ttc
                && !assembly
                    .variation_events()
                    .iter()
                    .any(|e| e.ref_allele == "TTC" && e.alt_allele == "T"),
            ensure_p12_has_atg
                && !assembly
                    .variation_events()
                    .iter()
                    .any(|e| e.ref_allele == "A" && e.alt_allele == "ATG"),
        );
        for (i, h) in assembly.haplotypes.iter().enumerate() {
            dump_hap("AFTER_ensure_alt", i, h, pad);
            let ev = variation_events_for_haplotype(h, &ref_after_trim, &ref_bytes, pad, 0, "2");
            dump_events_compact(&format!("RUST_EVENTMAP hap{i}"), &ev);
            let java_full_u = java_source_eventmap_proposed(h, &ref_bytes, pad);
            let java_full_m = java_source_eventmap_merged(&java_full_u);
            dump_java_compact("  JAVA_SOURCE_FULLREF_UNMERGED", &java_full_u);
            dump_java_compact("  JAVA_SOURCE_FULLREF_MERGED", &java_full_m);
            let mut hap_on_apply = h.clone();
            if hap_on_apply.alignment_start_hap_wrt_ref
                >= ref_after_trim.alignment_start_hap_wrt_ref
                && ref_after_trim.alignment_start_hap_wrt_ref > 0
            {
                hap_on_apply.alignment_start_hap_wrt_ref -=
                    ref_after_trim.alignment_start_hap_wrt_ref;
            }
            let java_apply_u = java_source_eventmap_proposed(&hap_on_apply, apply_bases, apply_pad);
            let java_apply_m = java_source_eventmap_merged(&java_apply_u);
            dump_java_compact("  JAVA_SOURCE_APPLYWIN_UNMERGED", &java_apply_u);
            dump_java_compact("  JAVA_SOURCE_APPLYWIN_MERGED", &java_apply_m);
            if i >= n_before_alt {
                let created = if h.bases.windows(16).any(|w| w == b"CTTTTTCATGATGGAT")
                    && h.cigar.as_ref().is_some_and(|c| {
                        c.elements.len() == 1 && c.elements[0].operator == CigarOperator::Match
                    }) {
                    "apply_anchor_snp_haplotypes(T→G @92307333)"
                } else if hap_has_ttc(&ev) || hap_has_atg(&ev) {
                    "unexpected TTC/ATG synthetic hap"
                } else {
                    "apply_read_events_to_assembly(missing_snp chained)"
                };
                eprintln!(
                    "  NEW_HAP idx={i} created_by={created} alt_win={} 2d1m2i={} ttc={} atg={} rust_em_n={} java_apply_u_n={}",
                    h.bases.windows(ALT_WIN.len()).any(|w| w == ALT_WIN),
                    has_2d1m2i(h),
                    hap_has_ttc(&ev),
                    hap_has_atg(&ev),
                    ev.len(),
                    java_apply_u.len()
                );
            }
        }
        stage_table(
            "after_ensure_alt",
            &assembly,
            &ref_after_trim,
            &ref_bytes,
            pad,
        );

        let alt_after = assembly
            .haplotypes
            .iter()
            .find(|h| !h.is_reference && h.bases.windows(ALT_WIN.len()).any(|w| w == ALT_WIN));
        assert!(
            alt_after.is_some(),
            "canonical ALT_WIN must remain after trim/fix/ensure_alt"
        );
        let canon = alt_after.unwrap();
        // 6R.39: trim_modern uses Math.max(end+padding). Pre-fix this haplotype was
        // 256 bp (`79M2D1M2I174M`) because padding accumulated until clipped to the
        // extended-region end. Java-contract span is 159 bp.
        assert_eq!(canon.bases.len(), 159);
        assert_eq!(
            canon.cigar.as_ref().map(|c| c.to_gatk_string()).as_deref(),
            Some("79M2D1M2I77M")
        );
        assert_eq!(canon.alignment_start_hap_wrt_ref, 0);
        assert_eq!(canon.kmer_size, 85);
        assert!(
            ensure_p12_has_ttc && ensure_p12_has_atg,
            "ensure_p12 must inject TTC→T and A→ATG on trimmed state"
        );
        assert_eq!(n_before_alt, 2);
        assert_eq!(assembly.haplotypes.len(), 4);
        assert!(
            !assembly
                .variation_events()
                .iter()
                .any(|e| e.ref_allele == "TTC" && e.alt_allele == "T"),
            "document: ensure_alt clobbers injected TTC→T via apply_read_events_to_assembly"
        );
        assert!(
            !assembly
                .variation_events()
                .iter()
                .any(|e| e.ref_allele == "A" && e.alt_allele == "ATG"),
            "document: ensure_alt clobbers injected A→ATG via apply_read_events_to_assembly"
        );
        let any_ttc_hap = assembly.haplotypes.iter().any(|h| {
            hap_has_ttc(&variation_events_for_haplotype(
                h,
                &ref_after_trim,
                &ref_bytes,
                pad,
                0,
                "2",
            ))
        });
        let any_atg_hap = assembly.haplotypes.iter().any(|h| {
            hap_has_atg(&variation_events_for_haplotype(
                h,
                &ref_after_trim,
                &ref_bytes,
                pad,
                0,
                "2",
            ))
        });
        assert!(
            !any_ttc_hap && !any_atg_hap,
            "ensure_alt must not create TTC→T or A→ATG synthetic haplotypes (indel span < 5)"
        );
    }

    /// 6R.16 closure: production pre-first-PairHMM path, then post-HMM ensure_alt,
    /// then the hap/event set that would enter a second PairHMM if the hap count grew.
    #[test]
    fn six_r16_trimmed_post_supplement_haplotype_closure() {
        let Some((ref_bytes, pad, orig_region, dict)) = load_real_p12() else {
            eprintln!("Real-data P12 comparison unavailable");
            return;
        };
        let Some(alt) = control_haplotype(&ref_bytes, pad) else {
            panic!("control haplotype construction failed");
        };
        let mut ref_hap = Haplotype::new(ref_bytes.as_slice(), true);
        let mut rc = Cigar::new();
        rc.push(ref_bytes.len(), CigarOperator::Match);
        ref_hap.cigar = Some(rc);
        ref_hap.genome_loc = alt.genome_loc;

        let untrimmed = AssemblyResultSet::from_assembly_for_calling_owned(
            AssemblyStatus::AssembledSomeVariation,
            85,
            vec![ref_hap.clone(), alt.clone()],
            ref_bytes.clone(),
            pad,
            "2",
            DEFAULT_MAX_MNP_DISTANCE,
        );
        let trim_variants: Vec<TrimVariant> = untrimmed
            .variation_events()
            .iter()
            .map(|e| TrimVariant {
                contig: e.contig.clone(),
                start: e.start_1based.get(),
                end: e.end_1based.get(),
                is_indel: e.is_indel(),
            })
            .collect();
        let trimmer =
            AssemblyRegionTrimmer::new(AssemblyRegionTrimmerConfig::gatk_defaults(), &dict, "2");
        let trim_result = trimmer.trim(&orig_region, &trim_variants, Some(&orig_region.reference));
        let region_for_genotyping = AssemblyRegionTrimmer::apply_trim(&orig_region, &trim_result);
        let mut assembly = untrimmed.trim_to(&region_for_genotyping).expect("trim_to");
        let sw = SwParameters::gatk_haplotype_to_reference();
        preserve_untrimmed_indel_haplotypes(&untrimmed, &mut assembly, &region_for_genotyping, &sw);

        let n_after_trim = assembly.haplotypes.len();
        let apply_bases = assembly.apply_bases_shared();
        let apply_pad = assembly
            .haplotypes
            .iter()
            .find(|h| h.is_reference)
            .and_then(|h| h.genome_loc.map(|g| g.start_1based()))
            .unwrap_or(pad);

        materialize_p12_cluster_from_assembly_cigars(
            &mut assembly,
            apply_bases.as_ref(),
            apply_pad,
            orig_region.start.get(),
            orig_region.end.get(),
            "2",
            &orig_region.reads,
            &sw,
        )
        .expect("materialize_p12");
        eprintln!(
            "AFTER_materialize_p12 n_haps {}→{} events={}",
            n_after_trim,
            assembly.haplotypes.len(),
            assembly.variation_events().len()
        );
        dump_events("AFTER_materialize_p12", assembly.variation_events());
        for (i, h) in assembly.haplotypes.iter().enumerate() {
            dump_hap("AFTER_materialize_p12", i, h, pad);
        }

        let n_after_mat = assembly.haplotypes.len();
        ensure_cluster_coupled_alt_haplotype(&mut assembly, apply_bases.as_ref(), apply_pad, &sw)
            .expect("ensure_cluster_coupled");
        eprintln!(
            "AFTER_ensure_cluster_coupled n_haps {}→{} events={}",
            n_after_mat,
            assembly.haplotypes.len(),
            assembly.variation_events().len()
        );
        dump_events("AFTER_ensure_cluster_coupled", assembly.variation_events());
        for (i, h) in assembly.haplotypes.iter().enumerate() {
            dump_hap("AFTER_ensure_cluster_coupled", i, h, pad);
        }

        let (full_ref, full_pad) = assembly.event_map_reference();
        let prior = assembly.variation_events.clone();
        assembly.variation_events = rebuild_variation_events(
            &assembly.haplotypes,
            full_ref,
            full_pad,
            "2",
            assembly.max_mnp_distance(),
            &prior,
            &[],
            RebuildVariationEventsOpts {
                event_map_only: false,
                merge_read_supplements: false,
            },
        );
        prune_spillover_supplement_haplotypes(&mut assembly);

        eprintln!(
            "PRE_FIRST_PAIRHMM n_haps={} events={} (engine.rs compute_region_read_likelihoods first call)",
            assembly.haplotypes.len(),
            assembly.variation_events().len()
        );
        dump_events("PRE_FIRST_PAIRHMM", assembly.variation_events());
        for (i, h) in assembly.haplotypes.iter().enumerate() {
            dump_hap("PRE_FIRST_PAIRHMM", i, h, pad);
        }
        let n_first_hmm = assembly.haplotypes.len();
        let first_hmm_has_canon = assembly.haplotypes.iter().any(|h| {
            !h.is_reference && h.bases.windows(ALT_WIN.len()).any(|w| w == ALT_WIN) && has_2d1m2i(h)
        });
        eprintln!("PRE_FIRST_PAIRHMM canon_p12={first_hmm_has_canon}");

        let n_before_fix = assembly.haplotypes.len();
        fix_p12_cluster_coupled_alt_haplotype(&mut assembly, "2", &sw);
        eprintln!(
            "CLOSURE_fix_p12 n_haps {}→{}",
            n_before_fix,
            assembly.haplotypes.len()
        );
        ensure_p12_cluster_variation_events_for_active_span(
            &mut assembly,
            "2",
            orig_region.start.get(),
            orig_region.end.get(),
        );
        let events_after_ensure_p12 = assembly.variation_events().to_vec();
        dump_events("CLOSURE_AFTER_ensure_p12", &events_after_ensure_p12);
        let n_before_alt = assembly.haplotypes.len();
        ensure_alt_haplotypes_for_variation_events(&mut assembly, &sw).expect("ensure_alt");
        prune_spillover_supplement_haplotypes(&mut assembly);
        eprintln!(
            "PRE_SECOND_PAIRHMM n_haps {}→{} events={} (engine.rs second compute_region_read_likelihoods if hap count grew)",
            n_before_alt,
            assembly.haplotypes.len(),
            assembly.variation_events().len()
        );
        dump_events("PRE_SECOND_PAIRHMM", assembly.variation_events());
        for (i, h) in assembly.haplotypes.iter().enumerate() {
            dump_hap("PRE_SECOND_PAIRHMM", i, h, pad);
        }

        let mut bases_seen: std::collections::BTreeMap<Vec<u8>, Vec<usize>> =
            std::collections::BTreeMap::new();
        for (i, h) in assembly.haplotypes.iter().enumerate() {
            bases_seen.entry(h.bases.clone()).or_default().push(i);
        }
        let dup_bases: Vec<_> = bases_seen
            .iter()
            .filter(|(_, idxs)| idxs.len() > 1)
            .collect();
        eprintln!("DUP_BASES groups={}", dup_bases.len());
        for (b, idxs) in &dup_bases {
            eprintln!(
                "  duplicate_bases idxs={idxs:?} len={} prefix={}",
                b.len(),
                String::from_utf8_lossy(&b[..b.len().min(16)])
            );
        }

        let ref_after = assembly
            .haplotypes
            .iter()
            .find(|h| h.is_reference)
            .cloned()
            .expect("ref");
        let mut ttc_on_hap = false;
        let mut atg_on_hap = false;
        for h in &assembly.haplotypes {
            let ev = variation_events_for_haplotype(h, &ref_after, &ref_bytes, pad, 0, "2");
            ttc_on_hap |= hap_has_ttc(&ev);
            atg_on_hap |= hap_has_atg(&ev);
        }
        let ttc_in_list = assembly
            .variation_events()
            .iter()
            .any(|e| e.ref_allele == "TTC" && e.alt_allele == "T");
        let atg_in_list = assembly
            .variation_events()
            .iter()
            .any(|e| e.ref_allele == "A" && e.alt_allele == "ATG");
        eprintln!(
            "CLOSURE ttc_in_variation_events={ttc_in_list} atg_in_variation_events={atg_in_list} ttc_on_any_hap_EventMap={ttc_on_hap} atg_on_any_hap_EventMap={atg_on_hap} first_hmm_n={n_first_hmm} second_hmm_n={}",
            assembly.haplotypes.len()
        );

        assert!(
            first_hmm_has_canon,
            "canonical coupled hap must enter first PairHMM"
        );
        assert!(
            assembly.haplotypes.iter().any(|h| !h.is_reference
                && h.bases.windows(ALT_WIN.len()).any(|w| w == ALT_WIN)
                && has_2d1m2i(h)),
            "canonical coupled hap must survive ensure_alt into second-PairHMM hap set"
        );
        assert!(
            events_after_ensure_p12
                .iter()
                .any(|e| e.ref_allele == "TTC" && e.alt_allele == "T"),
            "ensure_p12 injects TTC→T before ensure_alt"
        );
        assert!(
            !ttc_in_list && !atg_in_list,
            "document: variation_events not closed under injected TTC/ATG after ensure_alt"
        );
        assert!(
            !ttc_on_hap && !atg_on_hap,
            "document: no hap EventMap carries TTC→T or A→ATG (CIGAR EventMap is TTT/TAT/A→G)"
        );
        assert!(
            dup_bases.is_empty(),
            "no byte-identical haplotype duplicates after ensure_alt"
        );
        assert!(
            assembly.haplotypes.len() >= n_first_hmm,
            "post-HMM bridges must not shrink the first-PairHMM hap set"
        );
    }

    fn ve(start: u64, end: u64, r: &str, a: &str) -> VariationEvent {
        VariationEvent {
            contig: "2".into(),
            start_1based: GenomePosition::new_1based(start),
            end_1based: GenomePosition::new_1based(end),
            ref_allele: r.into(),
            alt_allele: a.into(),
        }
    }

    fn ref_only_assembly(ref_bytes: &[u8], pad: u64) -> AssemblyResultSet {
        let mut ref_hap = Haplotype::new(ref_bytes, true);
        let mut rc = Cigar::new();
        rc.push(ref_bytes.len(), CigarOperator::Match);
        ref_hap.cigar = Some(rc);
        ref_hap.genome_loc = Some(GenomeLoc::new(
            pad,
            pad.saturating_add(ref_bytes.len() as u64).saturating_sub(1),
        ));
        AssemblyResultSet::from_assembly_for_calling_owned(
            AssemblyStatus::AssembledSomeVariation,
            85,
            vec![ref_hap],
            ref_bytes.to_vec(),
            pad,
            "2",
            DEFAULT_MAX_MNP_DISTANCE,
        )
    }

    fn dump_ensure_alt_result(label: &str, assembly: &AssemblyResultSet, pad: u64) {
        eprintln!(
            "ENSURE_ALT {label} n_haps={} events={}",
            assembly.haplotypes.len(),
            assembly.variation_events().len()
        );
        dump_events(
            &format!("{label}_variation_events"),
            assembly.variation_events(),
        );
        for (i, h) in assembly.haplotypes.iter().enumerate() {
            dump_hap(label, i, h, pad);
        }
    }

    /// 6R.16: event-to-hap materialization — REF-only coupled events, REF-only EventMap
    /// events, and canonical hap + injected TTC/ATG (duplicate check + mapper snapshot).
    #[test]
    fn six_r16_ensure_alt_event_to_hap_materialization() {
        let Some((ref_bytes, pad, _, _)) = load_real_p12() else {
            eprintln!("Real-data P12 comparison unavailable");
            return;
        };
        let sw = SwParameters::gatk_haplotype_to_reference();
        let ttc = ve(P12_CLUSTER_TTC_START, P12_CLUSTER_TTC_START + 2, "TTC", "T");
        let atg = ve(P12_CLUSTER_ATG_START, P12_CLUSTER_ATG_START, "A", "ATG");
        let ttt = ve(
            P12_CLUSTER_TTC_START - 1,
            P12_CLUSTER_TTC_START + 1,
            "TTT",
            "T",
        );
        let cat = ve(
            P12_CLUSTER_TTC_START + 2,
            P12_CLUSTER_TTC_START + 4,
            "C",
            "CAT",
        );
        let ag = ve(P12_CLUSTER_ATG_START, P12_CLUSTER_ATG_START, "A", "G");

        let mut coupled_only = ref_only_assembly(&ref_bytes, pad);
        coupled_only.variation_events = vec![ttc.clone(), atg.clone()];
        coupled_only.variation_present = true;
        let n0 = coupled_only.haplotypes.len();
        ensure_alt_haplotypes_for_variation_events(&mut coupled_only, &sw)
            .expect("ensure_alt coupled");
        dump_ensure_alt_result("REF_ONLY_TTC_ATG", &coupled_only, pad);
        let coupled_alt: Vec<_> = coupled_only
            .haplotypes
            .iter()
            .filter(|h| !h.is_reference)
            .collect();
        eprintln!(
            "REF_ONLY_TTC_ATG n_haps {}→{} n_alt={} any_ALT_WIN={} any_2d1m2i={}",
            n0,
            coupled_only.haplotypes.len(),
            coupled_alt.len(),
            coupled_alt
                .iter()
                .any(|h| h.bases.windows(ALT_WIN.len()).any(|w| w == ALT_WIN)),
            coupled_alt.iter().any(|h| has_2d1m2i(h))
        );
        for h in &coupled_alt {
            eprintln!(
                "  coupled_alt_cigar={} align={} score={} kmer={}",
                h.cigar
                    .as_ref()
                    .map(|c| c.to_gatk_string())
                    .unwrap_or_default(),
                h.alignment_start_hap_wrt_ref,
                h.score,
                h.kmer_size
            );
        }

        let mut emap_only = ref_only_assembly(&ref_bytes, pad);
        emap_only.variation_events = vec![ttt.clone(), cat.clone(), ag.clone()];
        emap_only.variation_present = true;
        ensure_alt_haplotypes_for_variation_events(&mut emap_only, &sw).expect("ensure_alt emap");
        dump_ensure_alt_result("REF_ONLY_TTT_CAT_AG", &emap_only, pad);
        let emap_alt: Vec<_> = emap_only
            .haplotypes
            .iter()
            .filter(|h| !h.is_reference)
            .collect();
        eprintln!(
            "REF_ONLY_TTT_CAT_AG n_alt={} any_ALT_WIN={} any_2d1m2i={}",
            emap_alt.len(),
            emap_alt
                .iter()
                .any(|h| h.bases.windows(ALT_WIN.len()).any(|w| w == ALT_WIN)),
            emap_alt.iter().any(|h| has_2d1m2i(h))
        );

        let Some(alt) = control_haplotype(&ref_bytes, pad) else {
            panic!("control haplotype construction failed");
        };
        let mut ref_hap = Haplotype::new(ref_bytes.as_slice(), true);
        let mut rc = Cigar::new();
        rc.push(ref_bytes.len(), CigarOperator::Match);
        ref_hap.cigar = Some(rc);
        ref_hap.genome_loc = alt.genome_loc;
        let mut with_canon = AssemblyResultSet::from_assembly_for_calling_owned(
            AssemblyStatus::AssembledSomeVariation,
            85,
            vec![ref_hap.clone(), alt.clone()],
            ref_bytes.clone(),
            pad,
            "2",
            DEFAULT_MAX_MNP_DISTANCE,
        );
        let before_bases: Vec<Vec<u8>> = with_canon
            .haplotypes
            .iter()
            .map(|h| h.bases.clone())
            .collect();
        with_canon.variation_events.push(ttc.clone());
        with_canon.variation_events.push(atg.clone());
        let n_before = with_canon.haplotypes.len();
        ensure_alt_haplotypes_for_variation_events(&mut with_canon, &sw)
            .expect("ensure_alt canonical+inject");
        dump_ensure_alt_result("CANON_PLUS_TTC_ATG", &with_canon, pad);
        let new_haps: Vec<_> = with_canon
            .haplotypes
            .iter()
            .enumerate()
            .filter(|(_, h)| !before_bases.iter().any(|b| b == &h.bases))
            .collect();
        eprintln!(
            "CANON_PLUS_TTC_ATG n_haps {}→{} n_new_bases={} canon_survives={} duplicate_ALT_WIN_count={}",
            n_before,
            with_canon.haplotypes.len(),
            new_haps.len(),
            with_canon
                .haplotypes
                .iter()
                .any(|h| h.bases.windows(ALT_WIN.len()).any(|w| w == ALT_WIN)),
            with_canon
                .haplotypes
                .iter()
                .filter(|h| h.bases.windows(ALT_WIN.len()).any(|w| w == ALT_WIN))
                .count()
        );
        for (i, h) in &new_haps {
            eprintln!(
                "  NEW_BASES idx={i} len={} cigar={} alt_win={}",
                h.bases.len(),
                h.cigar
                    .as_ref()
                    .map(|c| c.to_gatk_string())
                    .unwrap_or_default(),
                h.bases.windows(ALT_WIN.len()).any(|w| w == ALT_WIN)
            );
        }

        let mapper_ttc = create_allele_mapper(
            &ttc,
            P12_CLUSTER_TTC_START,
            &with_canon.haplotypes,
            pad,
            &ref_bytes,
            DEFAULT_MAX_MNP_DISTANCE,
            false,
        );
        let mapper_atg = create_allele_mapper(
            &atg,
            P12_CLUSTER_ATG_START,
            &with_canon.haplotypes,
            pad,
            &ref_bytes,
            DEFAULT_MAX_MNP_DISTANCE,
            false,
        );
        eprintln!(
            "MAPPER TTC alt_haps={:?} REF_haps={:?} ATG alt_haps={:?} REF_haps={:?}",
            mapper_ttc.alt_haplotype_indices,
            mapper_ttc.ref_haplotype_indices,
            mapper_atg.alt_haplotype_indices,
            mapper_atg.ref_haplotype_indices
        );

        assert_eq!(
            coupled_alt.len(),
            1,
            "REF-only TTC+ATG: all-reference path should synthesize one alt hap"
        );
        assert!(
            coupled_alt[0]
                .bases
                .windows(ALT_WIN.len())
                .any(|w| w == ALT_WIN),
            "REF-only TTC+ATG chained splice should contain canonical coupled window"
        );
        assert_eq!(
            with_canon
                .haplotypes
                .iter()
                .filter(|h| h.bases.windows(ALT_WIN.len()).any(|w| w == ALT_WIN))
                .count(),
            1,
            "canonical hap must not be duplicated by ensure_alt when already present"
        );
    }
}
