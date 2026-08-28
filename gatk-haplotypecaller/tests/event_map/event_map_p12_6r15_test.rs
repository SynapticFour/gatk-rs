//! 6R.15 TEST-ONLY: post-W-H1 `fix_p12` / `ensure_p12` / stored-events provenance.
//! Does not call W-H1. Does not change production algorithms.

#[cfg(test)]
mod traces {
    use crate::alignment::SwParameters;
    use crate::assembly_region_iterator::AssemblyRegion;
    use crate::assembly_result_set::{AssemblyResultSet, DEFAULT_MAX_MNP_DISTANCE};
    use crate::cigar::{Cigar, CigarOperator};
    use crate::event_map::{variation_events_for_haplotype, VariationEvent};
    use crate::feature_context::FeatureContext;
    use crate::genome_loc::{GenomeLoc, GenomePosition};
    use crate::haplotype::Haplotype;
    use crate::hc_allele_mapping::create_allele_mapper;
    use crate::hc_genotyping_engine::audit_stored_events_with_p12_cluster_anchors;
    use crate::read_event_discovery::{
        ensure_p12_cluster_variation_events_for_active_span, fix_p12_cluster_coupled_alt_haplotype,
        P12_CLUSTER_ATG_START, P12_CLUSTER_TTC_START, SUPPLEMENT_HAPLOTYPE_SCORE,
    };
    use crate::read_threading_assembler::AssemblyStatus;
    use crate::reference_context::ReferenceContext;
    use std::path::Path;

    const ALT_WIN: &[u8] = b"CTTTTATGTGATGTAT";
    const REAL_P12_ACTIVE_START: u64 = P12_CLUSTER_TTC_START - 96;
    const REAL_P12_ACTIVE_END: u64 = P12_CLUSTER_TTC_START + 76;

    fn load_real_p12_ref() -> Option<(Vec<u8>, u64)> {
        use crate::assembly_region_finalize::{
            assembly_reference_read, finalize_region_reads_for_assembly,
            gatk_min_tail_quality_for_assembly, records_to_assembly_reads,
        };
        use crate::read_model::ReadFilterParams;
        use crate::walker_apply::{call_disposition, AssemblyRegionCallDisposition};
        use crate::walker_traversal::{
            flatten_assembly_regions, traverse_assembly_region_walker, WalkerTraversalConfig,
        };
        use gatk_core::reference::{
            parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary,
        };

        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
        let ref_path = root.join("parity/realworld/assets/hs37d5.simple.fa");
        let bam = root.join("parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam");
        if !ref_path.is_file() || !bam.is_file() {
            return None;
        }
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
        let region = regions.iter().find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= P12_CLUSTER_TTC_START
                && r.end.get() >= P12_CLUSTER_ATG_START
        })?;
        let mut ref_cache = ReferenceWindowCache::new(ref_path.clone(), 4);
        let reference = assembly_reference_read(&dict, &mut ref_cache, region).ok()?;
        let _finalized = finalize_region_reads_for_assembly(
            &region.reads,
            region,
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
        Some((reference.bases, pad))
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

    fn dummy_region(start: u64, end: u64, ext_start: u64, ext_end: u64) -> AssemblyRegion {
        AssemblyRegion {
            contig: "2".into(),
            start: GenomePosition::new_1based(start),
            end: GenomePosition::new_1based(end),
            is_active: true,
            extended_start: GenomePosition::new_1based(ext_start),
            extended_end: GenomePosition::new_1based(ext_end),
            extension: 100,
            reads: Vec::new(),
            read_qnames: Vec::new(),
            reference: ReferenceContext::empty(),
            features: FeatureContext::empty(),
            pileup_loci: Vec::new(),
        }
    }

    #[derive(Clone)]
    struct HapSnap {
        bases: Vec<u8>,
        cigar: String,
        align: usize,
        score: f64,
        is_ref: bool,
        kmer: usize,
        genome: Option<GenomeLoc>,
    }

    fn snap_alt(assembly: &AssemblyResultSet) -> HapSnap {
        let h = assembly
            .haplotypes
            .iter()
            .find(|h| !h.is_reference)
            .expect("control alt hap");
        HapSnap {
            bases: h.bases.clone(),
            cigar: h
                .cigar
                .as_ref()
                .map(|c| c.to_gatk_string())
                .unwrap_or_default(),
            align: h.alignment_start_hap_wrt_ref,
            score: h.score,
            is_ref: h.is_reference,
            kmer: h.kmer_size,
            genome: h.genome_loc,
        }
    }

    fn dump_hap(label: &str, s: &HapSnap) {
        eprintln!(
            "HAP {label} len={} cigar={} align={} score={} is_ref={} kmer={} genome={:?} alt_win={}",
            s.bases.len(),
            s.cigar,
            s.align,
            s.score,
            s.is_ref,
            s.kmer,
            s.genome,
            s.bases.windows(ALT_WIN.len()).any(|w| w == ALT_WIN)
        );
    }

    fn in_active(e: &VariationEvent) -> bool {
        e.start_1based.get() >= REAL_P12_ACTIVE_START && e.start_1based.get() <= REAL_P12_ACTIVE_END
    }

    fn dump_events(label: &str, events: &[VariationEvent]) {
        let span: Vec<_> = events.iter().filter(|e| in_active(e)).collect();
        eprintln!(
            "EVENTS {label} n_all={} n_active={}",
            events.len(),
            span.len()
        );
        for e in span {
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

    fn has_allele(events: &[VariationEvent], start: u64, r: &str, a: &str) -> bool {
        events
            .iter()
            .any(|e| e.start_1based.get() == start && e.ref_allele == r && e.alt_allele == a)
    }

    fn classify_fix(before: &HapSnap, after: &HapSnap, events_changed: bool) -> &'static str {
        let bases = before.bases != after.bases;
        let cigar = before.cigar != after.cigar;
        let meta = before.align != after.align
            || (before.score - after.score).abs() > 1e-12
            || before.is_ref != after.is_ref
            || before.kmer != after.kmer
            || before.genome != after.genome;
        if events_changed && !bases && !cigar {
            return "E";
        }
        if bases {
            return "D";
        }
        if cigar && !bases {
            return "C";
        }
        if meta {
            return "B";
        }
        if events_changed {
            return "E";
        }
        "A"
    }

    #[test]
    fn six_r15_p12_post_supplement_audit() {
        let Some((ref_bytes, pad)) = load_real_p12_ref() else {
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

        assert!(
            alt.bases.windows(ALT_WIN.len()).any(|w| w == ALT_WIN),
            "control hap must contain canonical ALT_WIN"
        );
        assert_eq!(
            alt.cigar.as_ref().map(|c| c.to_gatk_string()).as_deref(),
            Some("696M2D1M2I674M")
        );

        let mut assembly = AssemblyResultSet::from_assembly_for_calling_owned(
            AssemblyStatus::AssembledSomeVariation,
            85,
            vec![ref_hap.clone(), alt.clone()],
            ref_bytes.clone(),
            pad,
            "2",
            DEFAULT_MAX_MNP_DISTANCE,
        );

        let hap_events_before =
            variation_events_for_haplotype(&alt, &ref_hap, &ref_bytes, pad, 0, "2");
        dump_events("ordinary_eventmap_hap", &hap_events_before);
        dump_events(
            "assembly_variation_events_before_fix_p12",
            assembly.variation_events(),
        );

        let before = snap_alt(&assembly);
        dump_hap("BEFORE_fix_p12", &before);
        let events_before_fix = assembly.variation_events().to_vec();

        let sw = SwParameters::gatk_haplotype_to_reference();
        fix_p12_cluster_coupled_alt_haplotype(&mut assembly, "2", &sw);

        let after = snap_alt(&assembly);
        dump_hap("AFTER_fix_p12", &after);
        dump_events(
            "assembly_variation_events_after_fix_p12",
            assembly.variation_events(),
        );
        let events_changed_by_fix = events_before_fix != assembly.variation_events();
        let fix_class = classify_fix(&before, &after, events_changed_by_fix);
        eprintln!(
            "FIX_P12_CLASS={fix_class} bases_eq={} cigar_eq={} align_eq={} score_eq={} is_ref_eq={} kmer {}→{} events_changed={events_changed_by_fix} n_haps={}",
            before.bases == after.bases,
            before.cigar == after.cigar,
            before.align == after.align,
            (before.score - after.score).abs() < 1e-12,
            before.is_ref == after.is_ref,
            before.kmer,
            after.kmer,
            assembly.haplotypes.len()
        );
        eprintln!("FIX_P12_SOURCE=fix_p12_cluster_coupled_alt_haplotype");

        assert_eq!(
            before.bases, after.bases,
            "canonical bases must survive fix_p12"
        );
        assert_eq!(
            before.cigar, after.cigar,
            "canonical CIGAR must survive fix_p12"
        );
        assert_eq!(before.align, after.align);
        assert!((before.score - after.score).abs() < 1e-12);
        assert!(!after.is_ref);
        assert_eq!(after.cigar, "696M2D1M2I674M");
        assert!(
            after.bases.windows(ALT_WIN.len()).any(|w| w == ALT_WIN),
            "ALT_WIN must remain after fix_p12"
        );

        let before_ensure = assembly.variation_events().to_vec();
        dump_events("BEFORE_ensure_p12", &before_ensure);

        let region = dummy_region(
            REAL_P12_ACTIVE_START,
            REAL_P12_ACTIVE_END,
            REAL_P12_ACTIVE_START.saturating_sub(100),
            REAL_P12_ACTIVE_END.saturating_add(100),
        );
        ensure_p12_cluster_variation_events_for_active_span(
            &mut assembly,
            "2",
            region.start.get(),
            region.end.get(),
        );
        let after_ensure = assembly.variation_events().to_vec();
        dump_events("AFTER_ensure_p12", &after_ensure);
        eprintln!("ENSURE_P12_SOURCE=ensure_p12_cluster_variation_events_for_active_span");
        eprintln!(
            "ENSURE_has_TTC_T={} ENSURE_has_A_ATG={} ENSURE_has_TTT_T={} ENSURE_has_C_TAT={} ENSURE_has_C_CAT={} ENSURE_has_A_G={}",
            has_allele(&after_ensure, P12_CLUSTER_TTC_START, "TTC", "T"),
            has_allele(&after_ensure, P12_CLUSTER_ATG_START, "A", "ATG"),
            has_allele(&after_ensure, P12_CLUSTER_TTC_START - 1, "TTT", "T"),
            has_allele(&after_ensure, P12_CLUSTER_TTC_START + 2, "C", "TAT"),
            has_allele(&after_ensure, P12_CLUSTER_TTC_START + 2, "C", "CAT"),
            has_allele(&after_ensure, P12_CLUSTER_ATG_START, "A", "G"),
        );

        let hap_after = assembly
            .haplotypes
            .iter()
            .find(|h| !h.is_reference)
            .expect("alt");
        let hap_events_after_ensure =
            variation_events_for_haplotype(hap_after, &ref_hap, &ref_bytes, pad, 0, "2");
        dump_events(
            "ordinary_eventmap_hap_after_ensure",
            &hap_events_after_ensure,
        );

        let stored = audit_stored_events_with_p12_cluster_anchors(
            &after_ensure,
            &ref_bytes,
            pad,
            REAL_P12_ACTIVE_START,
            REAL_P12_ACTIVE_END,
            "2",
        );
        dump_events("stored_events_with_p12_cluster_anchors", &stored);
        eprintln!(
            "STORED_has_TTC_T={} STORED_has_A_ATG={} STORED_has_TTT_T={} STORED_has_C_TAT={} STORED_has_C_CAT={} STORED_has_A_G={}",
            has_allele(&stored, P12_CLUSTER_TTC_START, "TTC", "T"),
            has_allele(&stored, P12_CLUSTER_ATG_START, "A", "ATG"),
            has_allele(&stored, P12_CLUSTER_TTC_START - 1, "TTT", "T"),
            has_allele(&stored, P12_CLUSTER_TTC_START + 2, "C", "TAT"),
            has_allele(&stored, P12_CLUSTER_TTC_START + 2, "C", "CAT"),
            has_allele(&stored, P12_CLUSTER_ATG_START, "A", "G"),
        );

        let stored_from_eventmap = audit_stored_events_with_p12_cluster_anchors(
            &hap_events_before,
            &ref_bytes,
            pad,
            REAL_P12_ACTIVE_START,
            REAL_P12_ACTIVE_END,
            "2",
        );
        dump_events(
            "stored_events_from_ordinary_eventmap_only",
            &stored_from_eventmap,
        );

        let haps = assembly.haplotypes.clone();
        eprintln!("=== MAPPER after ensure/stored biological alleles ===");
        for (r, a, start) in [
            ("TTC", "T", P12_CLUSTER_TTC_START),
            ("A", "ATG", P12_CLUSTER_ATG_START),
            ("TTT", "T", P12_CLUSTER_TTC_START - 1),
            ("C", "TAT", P12_CLUSTER_TTC_START + 2),
            ("C", "CAT", P12_CLUSTER_TTC_START + 2),
            ("A", "G", P12_CLUSTER_ATG_START),
        ] {
            let Some(ev) = stored
                .iter()
                .chain(hap_events_after_ensure.iter())
                .find(|e| e.start_1based.get() == start && e.ref_allele == r && e.alt_allele == a)
                .cloned()
                .or_else(|| {
                    after_ensure
                        .iter()
                        .find(|e| {
                            e.start_1based.get() == start && e.ref_allele == r && e.alt_allele == a
                        })
                        .cloned()
                })
            else {
                eprintln!("  MAPPER missing {r}:{a} @{start}");
                continue;
            };
            let mapper = create_allele_mapper(&ev, start, &haps, pad, &ref_bytes, 0, true);
            eprintln!(
                "  loc={start} {r}:{a} ref_haps={:?} alt_haps={:?}",
                mapper.ref_haplotype_indices, mapper.alt_haplotype_indices
            );
        }

        assert!(
            has_allele(&after_ensure, P12_CLUSTER_TTC_START, "TTC", "T"),
            "ensure_p12 must inject TTC→T"
        );
        assert!(
            has_allele(&after_ensure, P12_CLUSTER_ATG_START, "A", "ATG"),
            "ensure_p12 must inject A→ATG"
        );
        assert!(
            has_allele(&after_ensure, P12_CLUSTER_TTC_START - 1, "TTT", "T"),
            "EventMap TTT→T must remain after ensure_p12"
        );
        assert!(
            has_allele(&after_ensure, P12_CLUSTER_TTC_START + 2, "C", "TAT"),
            "EventMap C→TAT (Java 4.4 makeBlock of C/T + C/CAT) must remain after ensure_p12"
        );
        assert!(
            !has_allele(&after_ensure, P12_CLUSTER_TTC_START + 2, "C", "CAT"),
            "production EventMap must not keep unmerged C→CAT after makeBlock"
        );
        assert!(
            !has_allele(&after_ensure, P12_CLUSTER_ATG_START, "A", "G"),
            "colocated A→G SNP must be dropped once A→ATG is injected"
        );
        assert!(
            !has_allele(&hap_events_before, P12_CLUSTER_TTC_START, "TTC", "T"),
            "ordinary EventMap must not emit TTC→T"
        );
        assert!(
            !has_allele(&hap_events_before, P12_CLUSTER_ATG_START, "A", "ATG"),
            "ordinary EventMap must not emit A→ATG"
        );
        assert!(
            has_allele(&stored, P12_CLUSTER_TTC_START, "TTC", "T")
                && has_allele(&stored, P12_CLUSTER_ATG_START, "A", "ATG"),
            "genotyping stored events must contain biological cluster alleles"
        );
        assert!(
            has_allele(&stored_from_eventmap, P12_CLUSTER_TTC_START, "TTC", "T")
                && has_allele(&stored_from_eventmap, P12_CLUSTER_ATG_START, "A", "ATG"),
            "stored_events injection must add TTC/ATG even from ordinary EventMap input"
        );
    }
}
