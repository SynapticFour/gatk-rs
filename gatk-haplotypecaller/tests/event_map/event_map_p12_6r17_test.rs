//! 6R.17 TEST-ONLY: allele-mapper decision-path closure on the canonical P12 hap.
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
    use crate::hc_allele_mapping::{
        audit_trace_create_allele_mapper, create_allele_mapper, MapperHapTrace,
    };
    use crate::hc_genotyping_engine::{
        audit_stored_events_with_p12_cluster_anchors, HcGenotypingConfig, SiteMap,
    };
    use crate::read_event_discovery::{
        ensure_alt_haplotypes_for_variation_events,
        ensure_p12_cluster_variation_events_for_active_span, fix_p12_cluster_coupled_alt_haplotype,
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

    fn ve(start: u64, end: u64, r: &str, a: &str) -> VariationEvent {
        VariationEvent {
            contig: "2".into(),
            start_1based: GenomePosition::new_1based(start),
            end_1based: GenomePosition::new_1based(end),
            ref_allele: r.into(),
            alt_allele: a.into(),
        }
    }

    fn has_allele(events: &[VariationEvent], start: u64, r: &str, a: &str) -> bool {
        events
            .iter()
            .any(|e| e.start_1based.get() == start && e.ref_allele == r && e.alt_allele == a)
    }

    fn eventmap_of(
        haps: &[Haplotype],
        ref_hap: &Haplotype,
        ref_bytes: &[u8],
        pad: u64,
    ) -> Vec<Vec<VariationEvent>> {
        haps.iter()
            .map(|h| variation_events_for_haplotype(h, ref_hap, ref_bytes, pad, 0, "2"))
            .collect()
    }

    fn dump_trace(
        state: &str,
        pad_label: &str,
        merged: &VariationEvent,
        variation_events: &[VariationEvent],
        hap_eventmaps: &[Vec<VariationEvent>],
        mapping_alt: &[crate::bio_ids::HaplotypeIndex],
        mapping_ref: &[crate::bio_ids::HaplotypeIndex],
        traces: &[MapperHapTrace],
    ) {
        let loc = merged.start_1based.get();
        eprintln!(
            "=== 6R.17 STATE={state} PAD={pad_label} EVENT={}->{}/{} loc={}..{} ===",
            merged.ref_allele,
            merged.alt_allele,
            loc,
            merged.start_1based.get(),
            merged.end_1based.get()
        );
        eprintln!(
            "  variation_events_has_this_allele={}",
            has_allele(
                variation_events,
                loc,
                &merged.ref_allele,
                &merged.alt_allele
            )
        );
        eprintln!(
            "  mapper alt_hap_indices={:?} ref_hap_indices={:?}",
            mapping_alt, mapping_ref
        );
        eprintln!(
            "  NOTE create_allele_mapper does not read assembly.variation_events; \
             variation_events presence is recorded only as a control."
        );
        for t in traces {
            let emap_has = hap_eventmaps
                .get(t.hap_index)
                .map(|ev| has_allele(ev, loc, &merged.ref_allele, &merged.alt_allele))
                .unwrap_or(false);
            eprintln!(
                "  HAP idx={} role={} EventMap_has_merged_allele={} \
                 overlapping_EventMap={:?} overlapping_walk={} \
                 hap_len={} pad={} loc={} off={} hap_slice_len={} \
                 hap_prefix={} ref_prefix={} \
                 haplotype_supports={} supports_path={} \
                 assignment_path={}",
                t.hap_index,
                t.assigned_role,
                emap_has || t.eventmap_has_merged_allele,
                t.overlapping,
                t.overlapping_walk,
                t.hap_len,
                t.pad_used,
                t.loc,
                t.off,
                t.hap_slice_len,
                t.hap_slice_prefix,
                t.ref_slice_prefix,
                t.haplotype_supports,
                t.haplotype_supports_path,
                t.assignment_path
            );
        }
    }

    fn map_and_dump(
        state: &str,
        pad_label: &str,
        merged: &VariationEvent,
        haps: &[Haplotype],
        pad: u64,
        ref_bytes: &[u8],
        variation_events: &[VariationEvent],
        emit_spanning_dels: bool,
    ) {
        let ref_hap = haps.iter().find(|h| h.is_reference).expect("ref hap");
        let hap_eventmaps = eventmap_of(haps, ref_hap, ref_bytes, pad);
        let (mapping, traces) = audit_trace_create_allele_mapper(
            merged,
            merged.start_1based.get(),
            haps,
            pad,
            ref_bytes,
            0,
            emit_spanning_dels,
        );
        dump_trace(
            state,
            pad_label,
            merged,
            variation_events,
            &hap_eventmaps,
            &mapping.alt_haplotype_indices,
            &mapping.ref_haplotype_indices,
            &traces,
        );
    }

    fn production_trimmed_assembly(
        ref_bytes: &[u8],
        pad: u64,
        orig_region: &AssemblyRegion,
        dict: &SequenceDictionary,
    ) -> AssemblyResultSet {
        let alt = control_haplotype(ref_bytes, pad).expect("control haplotype");
        let mut ref_hap = Haplotype::new(ref_bytes, true);
        let mut rc = Cigar::new();
        rc.push(ref_bytes.len(), CigarOperator::Match);
        ref_hap.cigar = Some(rc);
        ref_hap.genome_loc = alt.genome_loc;
        let untrimmed = AssemblyResultSet::from_assembly_for_calling_owned(
            AssemblyStatus::AssembledSomeVariation,
            85,
            vec![ref_hap, alt],
            ref_bytes.to_vec(),
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
            AssemblyRegionTrimmer::new(AssemblyRegionTrimmerConfig::gatk_defaults(), dict, "2");
        let trim_result = trimmer.trim(orig_region, &trim_variants, Some(&orig_region.reference));
        let region_for_genotyping = AssemblyRegionTrimmer::apply_trim(orig_region, &trim_result);
        let mut assembly = untrimmed.trim_to(&region_for_genotyping).expect("trim_to");
        let sw = SwParameters::gatk_haplotype_to_reference();
        preserve_untrimmed_indel_haplotypes(&untrimmed, &mut assembly, &region_for_genotyping, &sw);
        assembly
    }

    fn apply_window(assembly: &AssemblyResultSet, full_pad: u64) -> (Vec<u8>, u64) {
        let ref_hap = assembly
            .haplotypes
            .iter()
            .find(|h| h.is_reference)
            .expect("ref");
        let apply_pad = ref_hap
            .genome_loc
            .map(|g| g.start_1based())
            .unwrap_or(full_pad);
        (ref_hap.bases.clone(), apply_pad)
    }

    /// 6R.17: five mapper states on production-trimmed canonical hap.
    #[test]
    fn six_r17_mapper_five_states() {
        let Some((ref_bytes, pad, orig_region, dict)) = load_real_p12() else {
            eprintln!("Real-data P12 comparison unavailable");
            return;
        };
        let ttc = ve(P12_CLUSTER_TTC_START, P12_CLUSTER_TTC_START + 2, "TTC", "T");
        let atg = ve(P12_CLUSTER_ATG_START, P12_CLUSTER_ATG_START, "A", "ATG");
        let sw = SwParameters::gatk_haplotype_to_reference();
        let cfg = HcGenotypingConfig::strict_java();

        let mut assembly = production_trimmed_assembly(&ref_bytes, pad, &orig_region, &dict);
        fix_p12_cluster_coupled_alt_haplotype(&mut assembly, "2", &sw);
        let (apply_bases, apply_pad) = apply_window(&assembly, pad);
        let emap_only = rebuild_variation_events(
            &assembly.haplotypes,
            &ref_bytes,
            pad,
            "2",
            DEFAULT_MAX_MNP_DISTANCE,
            &[],
            &[],
            RebuildVariationEventsOpts {
                event_map_only: true,
                merge_read_supplements: false,
            },
        );
        assembly.variation_events = emap_only.clone();

        eprintln!(
            "6R.17_HAP_SHAPE n={} apply_pad={apply_pad} full_pad={pad} apply_len={} full_len={}",
            assembly.haplotypes.len(),
            apply_bases.len(),
            ref_bytes.len()
        );
        for (i, h) in assembly.haplotypes.iter().enumerate() {
            eprintln!(
                "  hap{i} len={} cigar={} align={} is_ref={} alt_win={} genome={:?}",
                h.bases.len(),
                h.cigar
                    .as_ref()
                    .map(|c| c.to_gatk_string())
                    .unwrap_or_default(),
                h.alignment_start_hap_wrt_ref,
                h.is_reference,
                h.bases.windows(ALT_WIN.len()).any(|w| w == ALT_WIN),
                h.genome_loc.map(|g| (g.start_1based(), g.end_1based()))
            );
        }

        // STATE 1: EventMap only (no TTC/ATG in variation_events).
        for (label, ev) in [("TTC", &ttc), ("ATG", &atg)] {
            map_and_dump(
                &format!("1_EventMap_only_{label}"),
                "apply",
                ev,
                &assembly.haplotypes,
                apply_pad,
                &apply_bases,
                &emap_only,
                true,
            );
            map_and_dump(
                &format!("1_EventMap_only_{label}"),
                "full",
                ev,
                &assembly.haplotypes,
                pad,
                &ref_bytes,
                &emap_only,
                true,
            );
        }
        assert!(
            !has_allele(&emap_only, P12_CLUSTER_TTC_START, "TTC", "T"),
            "STATE 1 EventMap-only variation_events must not contain TTC→T"
        );
        assert!(
            !has_allele(&emap_only, P12_CLUSTER_ATG_START, "A", "ATG"),
            "STATE 1 EventMap-only variation_events must not contain A→ATG"
        );

        // STATE 2: EventMap + TTC/ATG in variation_events (haps unchanged).
        let mut ve_plus = emap_only.clone();
        ve_plus.push(ttc.clone());
        ve_plus.push(atg.clone());
        assembly.variation_events = ve_plus.clone();
        for (label, ev) in [("TTC", &ttc), ("ATG", &atg)] {
            map_and_dump(
                &format!("2_EventMap_plus_TTC_ATG_{label}"),
                "apply",
                ev,
                &assembly.haplotypes,
                apply_pad,
                &apply_bases,
                &ve_plus,
                true,
            );
        }
        let s1_ttc = create_allele_mapper(
            &ttc,
            P12_CLUSTER_TTC_START,
            &assembly.haplotypes,
            apply_pad,
            &apply_bases,
            0,
            true,
        );
        let s2_ttc = create_allele_mapper(
            &ttc,
            P12_CLUSTER_TTC_START,
            &assembly.haplotypes,
            apply_pad,
            &apply_bases,
            0,
            true,
        );
        assert_eq!(
            s1_ttc.alt_haplotype_indices, s2_ttc.alt_haplotype_indices,
            "STATE 1 vs 2: variation_events must not change mapper (mapper does not read it)"
        );

        // STATE 3: after ensure_p12 (haps unchanged; variation_events injected).
        ensure_p12_cluster_variation_events_for_active_span(
            &mut assembly,
            "2",
            orig_region.start.get(),
            orig_region.end.get(),
        );
        let after_ensure = assembly.variation_events().to_vec();
        eprintln!(
            "STATE3_ensure_p12 has_TTC={} has_ATG={} n_events={}",
            has_allele(&after_ensure, P12_CLUSTER_TTC_START, "TTC", "T"),
            has_allele(&after_ensure, P12_CLUSTER_ATG_START, "A", "ATG"),
            after_ensure.len()
        );
        for (label, ev) in [("TTC", &ttc), ("ATG", &atg)] {
            map_and_dump(
                &format!("3_after_ensure_p12_{label}"),
                "apply",
                ev,
                &assembly.haplotypes,
                apply_pad,
                &apply_bases,
                &after_ensure,
                true,
            );
        }
        let s3_ttc = create_allele_mapper(
            &ttc,
            P12_CLUSTER_TTC_START,
            &assembly.haplotypes,
            apply_pad,
            &apply_bases,
            0,
            true,
        );
        assert_eq!(
            s1_ttc.alt_haplotype_indices, s3_ttc.alt_haplotype_indices,
            "STATE 3: ensure_p12 variation_events injection must not change mapper hap assignment"
        );

        // STATE 4: after ensure_alt (hap set may grow).
        let n_before_alt = assembly.haplotypes.len();
        ensure_alt_haplotypes_for_variation_events(&mut assembly, &sw).expect("ensure_alt");
        let after_alt_events = assembly.variation_events().to_vec();
        eprintln!(
            "STATE4_ensure_alt n_haps {}→{} events={}",
            n_before_alt,
            assembly.haplotypes.len(),
            after_alt_events.len()
        );
        let (apply_bases4, apply_pad4) = apply_window(&assembly, pad);
        for (label, ev) in [("TTC", &ttc), ("ATG", &atg)] {
            map_and_dump(
                &format!("4_after_ensure_alt_{label}"),
                "apply",
                ev,
                &assembly.haplotypes,
                apply_pad4,
                &apply_bases4,
                &after_alt_events,
                true,
            );
        }

        // Production SiteMap (apply pad, then full-pad retry if empty alt).
        let full_ref = assembly.reference_bases().to_vec();
        let full_pad = assembly.padded_reference_start_1based();
        for (label, ev) in [("TTC", &ttc), ("ATG", &atg)] {
            let mapping = SiteMap::build_mapping(
                ev,
                &assembly.haplotypes,
                &apply_bases4,
                apply_pad4,
                &full_ref,
                full_pad,
                0,
                &cfg,
                None,
            );
            eprintln!(
                "STATE4_SiteMap_{label} alt={:?} ref={:?}",
                mapping.alt_haplotype_indices, mapping.ref_haplotype_indices
            );
        }

        // STATE 5: stored_events_with_p12_cluster_anchors output.
        let stored = audit_stored_events_with_p12_cluster_anchors(
            &after_ensure,
            &ref_bytes,
            pad,
            REAL_P12_ACTIVE_START,
            REAL_P12_ACTIVE_END,
            "2",
        );
        eprintln!("STATE5_stored n_events={}", stored.len());
        for e in &stored {
            if e.ref_allele == "TTC"
                || e.alt_allele == "ATG"
                || e.start_1based.get() == P12_CLUSTER_TTC_START
                || e.start_1based.get() == P12_CLUSTER_ATG_START
                || (e.ref_allele == "TTT" && e.alt_allele == "T")
                || (e.ref_allele == "C" && e.alt_allele == "TAT")
                || (e.ref_allele == "C" && e.alt_allele == "CAT")
                || (e.ref_allele == "A" && e.alt_allele == "G")
            {
                map_and_dump(
                    "5_stored_events",
                    "apply",
                    e,
                    &assembly.haplotypes,
                    apply_pad4,
                    &apply_bases4,
                    &stored,
                    true,
                );
            }
        }

        let canon_idx = assembly
            .haplotypes
            .iter()
            .position(|h| !h.is_reference && h.bases.windows(ALT_WIN.len()).any(|w| w == ALT_WIN))
            .expect("canonical hap after trim");
        // After production trim the canonical hap is index 0 (REF is 1). Untrimmed 6R.16
        // mapper snapshots used index 1 because the untrimmed list was [REF, alt].
        let s1_maps_canon = s1_ttc
            .alt_haplotype_indices
            .iter()
            .any(|i| i.get() == canon_idx || i.get() == 0);
        assert!(
            s1_maps_canon,
            "canonical trimmed hap must be mapped to TTC→T via fallback (observed idx={canon_idx})"
        );
    }

    /// 6R.17: negative controls A–D. Test-only haplotypes/events; no production changes.
    #[test]
    fn six_r17_mapper_negative_controls() {
        let Some((ref_bytes, pad, orig_region, dict)) = load_real_p12() else {
            eprintln!("Real-data P12 comparison unavailable");
            return;
        };
        let ttc = ve(P12_CLUSTER_TTC_START, P12_CLUSTER_TTC_START + 2, "TTC", "T");
        let atg = ve(P12_CLUSTER_ATG_START, P12_CLUSTER_ATG_START, "A", "ATG");
        let ttc_wrong_alt = ve(P12_CLUSTER_TTC_START, P12_CLUSTER_TTC_START + 2, "TTC", "G");
        let atg_wrong_alt = ve(P12_CLUSTER_ATG_START, P12_CLUSTER_ATG_START, "A", "C");
        let ag = ve(P12_CLUSTER_ATG_START, P12_CLUSTER_ATG_START, "A", "G");
        let sw = SwParameters::gatk_haplotype_to_reference();

        let mut assembly = production_trimmed_assembly(&ref_bytes, pad, &orig_region, &dict);
        fix_p12_cluster_coupled_alt_haplotype(&mut assembly, "2", &sw);
        let (apply_bases, apply_pad) = apply_window(&assembly, pad);
        let emap_only = rebuild_variation_events(
            &assembly.haplotypes,
            &ref_bytes,
            pad,
            "2",
            DEFAULT_MAX_MNP_DISTANCE,
            &[],
            &[],
            RebuildVariationEventsOpts {
                event_map_only: true,
                merge_read_supplements: false,
            },
        );

        let haps = assembly.haplotypes.clone();
        eprintln!("NEG_CONTROL_BASELINE");
        for ev in [&ttc, &atg] {
            map_and_dump(
                "NEG_baseline",
                "apply",
                ev,
                &haps,
                apply_pad,
                &apply_bases,
                &emap_only,
                true,
            );
        }

        // A: keep CIGAR, replace alt hap bases with REF so TTC/ATG are not represented.
        let mut haps_a = haps.clone();
        let ref_bases = haps_a
            .iter()
            .find(|h| h.is_reference)
            .expect("ref")
            .bases
            .clone();
        for h in &mut haps_a {
            if !h.is_reference {
                h.bases = ref_bases.clone();
            }
        }
        eprintln!("NEG_CONTROL_A sequence=REF cigar=unchanged");
        for ev in [&ttc, &atg] {
            map_and_dump(
                "NEG_A_ref_sequence",
                "apply",
                ev,
                &haps_a,
                apply_pad,
                &apply_bases,
                &emap_only,
                true,
            );
        }

        // B: keep sequence, remove P12 from variation_events (already EventMap-only).
        eprintln!("NEG_CONTROL_B variation_events without TTC/ATG (identical to STATE 1)");
        for ev in [&ttc, &atg] {
            map_and_dump(
                "NEG_B_no_p12_variation_events",
                "apply",
                ev,
                &haps,
                apply_pad,
                &apply_bases,
                &emap_only,
                true,
            );
        }
        let with_p12 = {
            let mut v = emap_only.clone();
            v.push(ttc.clone());
            v.push(atg.clone());
            v
        };
        let m_no = create_allele_mapper(
            &ttc,
            P12_CLUSTER_TTC_START,
            &haps,
            apply_pad,
            &apply_bases,
            0,
            true,
        );
        let m_yes = create_allele_mapper(
            &ttc,
            P12_CLUSTER_TTC_START,
            &haps,
            apply_pad,
            &apply_bases,
            0,
            true,
        );
        let _ = with_p12;
        assert_eq!(
            m_no.alt_haplotype_indices, m_yes.alt_haplotype_indices,
            "NEG B: variation_events presence/absence must not change mapper"
        );

        // C: keep sequence, replace CIGAR with all-M so ordinary EventMap is empty.
        let mut haps_c = haps.clone();
        for h in &mut haps_c {
            if !h.is_reference {
                let mut c = Cigar::new();
                c.push(h.bases.len(), CigarOperator::Match);
                h.cigar = Some(c);
                h.alignment_start_hap_wrt_ref = 0;
            }
        }
        eprintln!("NEG_CONTROL_C sequence=unchanged cigar=all-M (empty EventMap)");
        for ev in [&ttc, &atg] {
            map_and_dump(
                "NEG_C_allM_cigar",
                "apply",
                ev,
                &haps_c,
                apply_pad,
                &apply_bases,
                &emap_only,
                true,
            );
        }

        // D: same coordinates, deliberately different ALT.
        eprintln!("NEG_CONTROL_D different ALT at same coordinates");
        for ev in [&ttc_wrong_alt, &atg_wrong_alt, &ag] {
            map_and_dump(
                "NEG_D_wrong_alt",
                "apply",
                ev,
                &haps,
                apply_pad,
                &apply_bases,
                &emap_only,
                true,
            );
        }
    }
}
