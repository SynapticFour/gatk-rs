//! 6R.19 production EventMap Java 4.4 parity regressions.
//!
//! Pins [`crate::event_map::make_block`] / [`crate::event_map::add_vc_merge`] and
//! insertion END (`start + REF.len() − 1`) after promotion from the 6R.18 test helper.

#[cfg(test)]
mod six_r19_tests {
    use crate::cigar::{Cigar, CigarOperator};
    use crate::event_map::{
        add_vc_merge, make_block, overlapping_events, variation_events_for_haplotype, EventMap,
        VariationEvent,
    };
    use crate::genome_loc::GenomeLoc;
    use crate::haplotype::Haplotype;
    use crate::read_event_discovery::{
        P12_CLUSTER_ATG_START, P12_CLUSTER_TTC_START, SUPPLEMENT_HAPLOTYPE_SCORE,
    };
    use std::path::Path;

    const CHR: &str = "20";
    const START: u64 = 10;
    const ALT_WIN: &[u8] = b"CTTTTATGTGATGTAT";

    fn vc(contig: &str, start: u64, r: &str, a: &str) -> VariationEvent {
        VariationEvent::from_alleles(contig, start, r, a)
    }

    fn assert_alleles(got: &VariationEvent, r: &str, a: &str, start: u64, end: u64, case: &str) {
        assert_eq!(got.ref_allele, r, "{case}: REF");
        assert_eq!(got.alt_allele, a, "{case}: ALT");
        assert_eq!(got.start_1based.get(), start, "{case}: start");
        assert_eq!(got.end_1based.get(), end, "{case}: end");
    }

    fn hap_with_cigar(bases: &[u8], elements: &[(usize, CigarOperator)]) -> Haplotype {
        let mut cigar = Cigar::new();
        for (len, op) in elements {
            cigar.push(*len, *op);
        }
        Haplotype {
            bases: bases.to_vec(),
            is_reference: false,
            score: 0.0,
            kmer_size: 10,
            cigar: Some(cigar),
            genome_loc: None,
            alignment_start_hap_wrt_ref: 0,
        }
    }

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

        let active_start = P12_CLUSTER_TTC_START - 96;
        let active_end = P12_CLUSTER_TTC_START + 76;
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
        let ref_path = root.join("parity/realworld/assets/hs37d5.simple.fa");
        let bam = root.join("parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam");
        if !ref_path.is_file() || !bam.is_file() {
            return None;
        }
        let dict = SequenceDictionary::from_fasta_path(&ref_path).ok()?;
        let interval = format!("2:{active_start}-{active_end}");
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

    /// Test 1 — Java `EventMapUnitTest.MakeBlockData` (GATK 4.4.0.0 SHA 2dbc025).
    #[test]
    fn six_r19_java_makeblock_unit_vectors() {
        let vectors: &[((&str, &str), (&str, &str), (&str, &str))] = &[
            (("A", "G"), ("AGT", "A"), ("AGT", "G")),
            (("A", "G"), ("A", "AGT"), ("A", "GGT")),
            (("AC", "A"), ("A", "AGT"), ("AC", "AGT")),
            (("ACGTA", "A"), ("A", "AG"), ("ACGTA", "AG")),
            (("AC", "A"), ("A", "AGCGT"), ("AC", "AGCGT")),
            (("A", "ACGTA"), ("AG", "A"), ("AG", "ACGTA")),
            (("A", "AC"), ("AGCGT", "A"), ("AGCGT", "AC")),
        ];
        for (i, &((r1, a1), (r2, a2), (er, ea))) in vectors.iter().enumerate() {
            let vc1 = vc(CHR, START, r1, a1);
            let vc2 = vc(CHR, START, r2, a2);
            let expected = vc(CHR, START, er, ea);
            let got = make_block(&vc1, &vc2).unwrap_or_else(|e| panic!("vector {i}: {e}"));
            assert_alleles(
                &got,
                er,
                ea,
                START,
                expected.end_1based.get(),
                &format!("MakeBlockData[{i}]"),
            );
        }
    }

    /// Test 2 — P12 same-start SNP + insertion → C/TAT at insertion start.
    #[test]
    fn six_r19_p12_snp_plus_insertion_ct_cat_to_tat() {
        let loc = P12_CLUSTER_TTC_START + 2;
        let got =
            make_block(&vc("2", loc, "C", "T"), &vc("2", loc, "C", "CAT")).expect("C/T + C/CAT");
        assert_alleles(&got, "C", "TAT", loc, loc, "P12 C/T + C/CAT");
    }

    /// Test 3 — A/G at ATG start stays a separate event after addVC.
    #[test]
    fn six_r19_p12_separate_snp_ag() {
        let loc_ins = P12_CLUSTER_TTC_START + 2;
        let loc_snp = P12_CLUSTER_ATG_START;
        let merged = add_vc_merge(vec![
            vc("2", loc_ins, "C", "T"),
            vc("2", loc_ins, "C", "CAT"),
            vc("2", loc_snp, "A", "G"),
        ]);
        assert_eq!(merged.len(), 2);
        assert_alleles(&merged[0], "C", "TAT", loc_ins, loc_ins, "TAT block");
        assert_alleles(&merged[1], "A", "G", loc_snp, loc_snp, "A/G separate");
    }

    /// Test 4 — pure insertion END is start, not start + |ALT| − 1.
    #[test]
    fn six_r19_insertion_end_start_eq_end() {
        let ins = vc("chr", 100, "A", "AGT");
        assert_eq!(ins.start_1based.get(), 100);
        assert_eq!(ins.end_1based.get(), 100, "Java insertion start == end");
        assert_ne!(ins.end_1based.get(), 102);

        let ref_bytes = b"ACGTAAAA";
        let ref_hap = Haplotype::new(ref_bytes, true);
        let alt = hap_with_cigar(
            b"AGTCGTAAAA",
            &[
                (1, CigarOperator::Match),
                (2, CigarOperator::Insertion),
                (7, CigarOperator::Match),
            ],
        );
        let events = EventMap::from_haplotype_and_reference(&alt, &ref_hap, ref_bytes, 100, 0)
            .variation_events("chr", 100);
        let ins_ev = events
            .iter()
            .find(|e| e.ref_allele == "A" && e.alt_allele == "AGT")
            .expect("production CIGAR path must emit A/AGT");
        assert_eq!(ins_ev.start_1based.get(), 100);
        assert_eq!(ins_ev.end_1based.get(), 100);
    }

    /// Test 5 — insertion / deletion / SNP / mixed-block overlap.
    #[test]
    fn six_r19_insertion_deletion_snp_mixed_overlap() {
        let ins = vc("chr", 100, "A", "AGT");
        assert_eq!(overlapping_events(&[ins.clone()], 100).len(), 1);
        assert_eq!(
            overlapping_events(&[ins], 101).len(),
            0,
            "pure insertion must not overlap start+1"
        );

        let del = vc("chr", 100, "AGT", "A");
        assert_eq!(del.end_1based.get(), 102);
        assert_eq!(overlapping_events(&[del.clone()], 100).len(), 1);
        assert_eq!(overlapping_events(&[del.clone()], 101).len(), 1);
        assert_eq!(overlapping_events(&[del.clone()], 102).len(), 1);
        assert_eq!(overlapping_events(&[del], 103).len(), 0);

        let snp = vc("chr", 100, "A", "G");
        assert_eq!(overlapping_events(&[snp.clone()], 100).len(), 1);
        assert_eq!(overlapping_events(&[snp], 101).len(), 0);

        let block = make_block(&vc(CHR, START, "AC", "A"), &vc(CHR, START, "A", "AGT"))
            .expect("del+ins MakeBlockData");
        assert_alleles(&block, "AC", "AGT", START, START + 1, "mixed AC/AGT");
        assert_eq!(overlapping_events(&[block.clone()], START).len(), 1);
        assert_eq!(overlapping_events(&[block.clone()], START + 1).len(), 1);
        assert_eq!(overlapping_events(&[block], START + 2).len(), 0);

        let p12_ins = vc("2", P12_CLUSTER_TTC_START + 2, "C", "CAT");
        assert_eq!(p12_ins.start_1based.get(), p12_ins.end_1based.get());
        assert_eq!(
            overlapping_events(&[p12_ins.clone()], P12_CLUSTER_TTC_START + 2).len(),
            1
        );
        assert_eq!(
            overlapping_events(&[p12_ins], P12_CLUSTER_ATG_START).len(),
            0,
            "C→CAT must not overlap the following A→G start"
        );
    }

    /// Always-on 2D1M2I composition: TTT→T, C→TAT, A→G (same merge as P12).
    #[test]
    fn six_r19_synthetic_2d1m2i_makeblock() {
        let ref_bytes = b"TTTCAGGGGG";
        let ref_hap = Haplotype::new(ref_bytes, true);
        let alt = hap_with_cigar(
            b"TTATGGGGG",
            &[
                (1, CigarOperator::Match),
                (2, CigarOperator::Deletion),
                (1, CigarOperator::Match),
                (2, CigarOperator::Insertion),
                (5, CigarOperator::Match),
            ],
        );
        let events = variation_events_for_haplotype(&alt, &ref_hap, ref_bytes, 1000, 0, "chr");
        eprintln!(
            "SYNTHETIC_2D1M2I {:?}",
            events
                .iter()
                .map(|e| format!(
                    "{}-{} {}→{}",
                    e.start_1based.get(),
                    e.end_1based.get(),
                    e.ref_allele,
                    e.alt_allele
                ))
                .collect::<Vec<_>>()
        );
        assert_eq!(events.len(), 3);
        assert_alleles(&events[0], "TTT", "T", 1000, 1002, "del");
        assert_alleles(&events[1], "C", "TAT", 1003, 1003, "makeBlock SNP+INS");
        assert_alleles(&events[2], "A", "G", 1004, 1004, "separate SNP");
    }

    /// Test 6 — production CIGAR path on canonical `696M2D1M2I674M`.
    #[test]
    fn six_r19_canonical_p12_cigar_production_eventmap() {
        let Some((ref_bytes, pad)) = load_real_p12_ref() else {
            eprintln!("Real-data P12 comparison unavailable");
            return;
        };
        let alt = control_haplotype(&ref_bytes, pad).expect("control haplotype");
        assert_eq!(
            alt.cigar.as_ref().map(|c| c.to_gatk_string()).as_deref(),
            Some("696M2D1M2I674M")
        );
        assert!(alt.bases.windows(ALT_WIN.len()).any(|w| w == ALT_WIN));

        let mut ref_hap = Haplotype::new(ref_bytes.as_slice(), true);
        let mut rc = Cigar::new();
        rc.push(ref_bytes.len(), CigarOperator::Match);
        ref_hap.cigar = Some(rc);
        ref_hap.genome_loc = alt.genome_loc;

        let events = variation_events_for_haplotype(&alt, &ref_hap, &ref_bytes, pad, 0, "2");
        eprintln!(
            "P12_PROD_EVENTMAP pad={pad} {:?}",
            events
                .iter()
                .map(|e| format!(
                    "{}-{} {}→{}",
                    e.start_1based.get(),
                    e.end_1based.get(),
                    e.ref_allele,
                    e.alt_allele
                ))
                .collect::<Vec<_>>()
        );
        assert_eq!(
            events
                .iter()
                .map(|e| (
                    e.start_1based.get(),
                    e.end_1based.get(),
                    e.ref_allele.as_str(),
                    e.alt_allele.as_str()
                ))
                .collect::<Vec<_>>(),
            vec![
                (
                    P12_CLUSTER_TTC_START - 1,
                    P12_CLUSTER_TTC_START + 1,
                    "TTT",
                    "T"
                ),
                (
                    P12_CLUSTER_TTC_START + 2,
                    P12_CLUSTER_TTC_START + 2,
                    "C",
                    "TAT"
                ),
                (P12_CLUSTER_ATG_START, P12_CLUSTER_ATG_START, "A", "G"),
            ]
        );
        assert!(
            !events
                .iter()
                .any(|e| e.ref_allele == "C" && e.alt_allele == "CAT"),
            "merged EventMap must not retain unmerged C→CAT"
        );
        assert!(
            !events
                .iter()
                .any(|e| e.ref_allele == "TTC" && e.alt_allele == "T"),
            "EventMap must not emit injected TTC→T"
        );
        assert!(
            !events
                .iter()
                .any(|e| e.ref_allele == "A" && e.alt_allele == "ATG"),
            "EventMap must not emit injected A→ATG"
        );
    }
}
