//! 6R.18 tests for GATK **4.4.0.0** `EventMap.makeBlock` / `addVC`.
//!
//! Production implementation lives in [`super`] (`make_block`, `add_vc_merge`,
//! `overlapping_events`). This file keeps the 6R.18 names as aliases so existing
//! tests continue to pin the Java 4.4 vectors.

use super::{add_vc_merge, make_block, overlapping_events, MakeBlockError, VariationEvent};

pub type Java44MakeBlockError = MakeBlockError;

pub fn java_44_end_from_ref_allele(start_1based: u64, ref_allele: &str) -> u64 {
    VariationEvent::vcf_end_1based(start_1based, ref_allele)
}

pub fn java_44_vc(
    contig: &str,
    start_1based: u64,
    ref_allele: &str,
    alt_allele: &str,
) -> VariationEvent {
    VariationEvent::from_alleles(contig, start_1based, ref_allele, alt_allele)
}

pub fn make_block_java_44(
    vc1: &VariationEvent,
    vc2: &VariationEvent,
) -> Result<VariationEvent, MakeBlockError> {
    make_block(vc1, vc2)
}

pub fn add_vc_merge_java_44(
    proposed: &[VariationEvent],
) -> Result<Vec<VariationEvent>, MakeBlockError> {
    Ok(add_vc_merge(proposed.to_vec()))
}

pub fn get_overlapping_events_java_44(
    events: &[VariationEvent],
    loc_1based: u64,
) -> Vec<VariationEvent> {
    overlapping_events(events, loc_1based)
}

#[cfg(test)]
mod six_r18_tests {
    use super::*;
    use crate::event_map::prefer_indel_over_colocated_snps;
    use crate::genome_loc::GenomePosition;

    const CHR: &str = "20";
    const START: u64 = 10;

    fn dump(label: &str, events: &[VariationEvent]) {
        eprintln!("EVENTS {label} n={}", events.len());
        for e in events {
            eprintln!(
                "  {}-{} {}→{}",
                e.start_1based.get(),
                e.end_1based.get(),
                e.ref_allele,
                e.alt_allele
            );
        }
    }

    fn assert_vc_eq(got: &VariationEvent, expected: &VariationEvent, case: &str) {
        assert_eq!(
            got.start_1based,
            expected.start_1based,
            "{case}: start {} vs {}",
            got.start_1based.get(),
            expected.start_1based.get()
        );
        assert_eq!(
            got.end_1based,
            expected.end_1based,
            "{case}: end {} vs {}",
            got.end_1based.get(),
            expected.end_1based.get()
        );
        assert_eq!(
            got.ref_allele, expected.ref_allele,
            "{case}: REF {} vs {}",
            got.ref_allele, expected.ref_allele
        );
        assert_eq!(
            got.alt_allele, expected.alt_allele,
            "{case}: ALT {} vs {}",
            got.alt_allele, expected.alt_allele
        );
    }

    /// OBSERVED JAVA TEST VECTOR: `EventMapUnitTest.MakeBlockData`.
    #[test]
    fn six_r18_java_makeblock_unit_vectors() {
        // Each row: (vc1 REF/ALT, vc2 REF/ALT, expected REF/ALT)
        // Source: EventMapUnitTest.java MakeBlockData @ SHA 2dbc025.
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
            let vc1 = java_44_vc(CHR, START, r1, a1);
            let vc2 = java_44_vc(CHR, START, r2, a2);
            let expected = java_44_vc(CHR, START, er, ea);
            let got = make_block_java_44(&vc1, &vc2).unwrap_or_else(|e| panic!("vector {i}: {e}"));
            eprintln!(
                "VECTOR {i}: {r1}/{a1} + {r2}/{a2} → {}/{} (expected {er}/{ea}) start={} end={}",
                got.ref_allele,
                got.alt_allele,
                got.start_1based.get(),
                got.end_1based.get()
            );
            assert_vc_eq(&got, &expected, &format!("MakeBlockData[{i}]"));
            assert_eq!(
                got.start_1based.get(),
                START,
                "Java test asserts start equality"
            );
        }
    }

    /// Case A — SNP + insertion from Java MakeBlockData.
    #[test]
    fn six_r18_case_a_snp_plus_insertion() {
        let vc1 = java_44_vc(CHR, START, "A", "G");
        let vc2 = java_44_vc(CHR, START, "A", "AGT");
        let expected = java_44_vc(CHR, START, "A", "GGT");
        let got = make_block_java_44(&vc1, &vc2).expect("Case A");
        dump("Case A input", &[vc1.clone(), vc2.clone()]);
        dump("Case A Java expected", &[expected.clone()]);
        dump("Case A Rust make_block_java_44", &[got.clone()]);
        assert_vc_eq(&got, &expected, "Case A A/G + A/AGT → A/GGT");
        assert_eq!(
            got.start_1based, got.end_1based,
            "SNP+INS keeps SNP stop (start==end)"
        );
    }

    /// Case B — P12-style SNP + insertion (same composition as Java A/G + A/AGT).
    #[test]
    fn six_r18_case_b_p12_snp_plus_insertion() {
        let vc1 = java_44_vc("2", 92_307_326, "C", "T");
        let vc2 = java_44_vc("2", 92_307_326, "C", "CAT");
        let expected = java_44_vc("2", 92_307_326, "C", "TAT");
        let got = make_block_java_44(&vc1, &vc2).expect("Case B");
        dump("Case B input", &[vc1.clone(), vc2.clone()]);
        dump("Case B Java expected (composition)", &[expected.clone()]);
        dump("Case B Rust make_block_java_44", &[got.clone()]);
        assert_vc_eq(&got, &expected, "Case B C/T + C/CAT → C/TAT");
        assert_eq!(got.end_1based.get(), 92_307_326);
    }

    /// Case C — SNP + deletion from Java MakeBlockData (first vector).
    #[test]
    fn six_r18_case_c_snp_plus_deletion() {
        let vc1 = java_44_vc(CHR, START, "A", "G");
        let vc2 = java_44_vc(CHR, START, "AGT", "A");
        let expected = java_44_vc(CHR, START, "AGT", "G");
        let got = make_block_java_44(&vc1, &vc2).expect("Case C");
        dump("Case C input", &[vc1.clone(), vc2.clone()]);
        dump("Case C Java expected", &[expected.clone()]);
        dump("Case C Rust make_block_java_44", &[got.clone()]);
        assert_vc_eq(&got, &expected, "Case C A/G + AGT/A → AGT/G");
        assert_eq!(
            got.end_1based.get(),
            java_44_end_from_ref_allele(START, "AGT"),
            "SNP+DEL stop is the deletion end"
        );
    }

    /// Case D — non-colocated events remain independent under addVC.
    #[test]
    fn six_r18_case_d_non_colocated() {
        let a = java_44_vc("2", 92_307_323, "TTT", "T");
        let b = java_44_vc("2", 92_307_326, "C", "CAT");
        let c = java_44_vc("2", 92_307_327, "A", "G");
        let got = add_vc_merge_java_44(&[a.clone(), b.clone(), c.clone()]).expect("Case D");
        dump("Case D input", &[a.clone(), b.clone(), c.clone()]);
        dump("Case D add_vc_merge_java_44", &got);
        assert_eq!(got.len(), 3);
        assert_vc_eq(&got[0], &a, "Case D TTT");
        assert_vc_eq(&got[1], &b, "Case D CAT");
        assert_vc_eq(&got[2], &c, "Case D AG");
    }

    /// Case E — Java MakeBlockData indel+indel (del+ins and ins+del) composites.
    #[test]
    fn six_r18_case_e_indel_plus_opposite_indel() {
        let rows: &[((&str, &str), (&str, &str), (&str, &str))] = &[
            (("AC", "A"), ("A", "AGT"), ("AC", "AGT")),
            (("A", "ACGTA"), ("AG", "A"), ("AG", "ACGTA")),
        ];
        for (i, &((r1, a1), (r2, a2), (er, ea))) in rows.iter().enumerate() {
            let vc1 = java_44_vc(CHR, START, r1, a1);
            let vc2 = java_44_vc(CHR, START, r2, a2);
            let expected = java_44_vc(CHR, START, er, ea);
            let got = make_block_java_44(&vc1, &vc2).unwrap_or_else(|e| panic!("Case E {i}: {e}"));
            assert_vc_eq(&got, &expected, &format!("Case E[{i}]"));
        }
    }

    /// Case F — Java insertion start == end; overlap at start but not start+1.
    #[test]
    fn six_r18_case_f_insertion_start_eq_end_overlap() {
        let ins = java_44_vc("2", 92_307_326, "C", "CAT");
        assert_eq!(
            ins.start_1based, ins.end_1based,
            "Java insertion VariantContextBuilder(start, start)"
        );
        let at_start = get_overlapping_events_java_44(&[ins.clone()], 92_307_326);
        let at_plus_one = get_overlapping_events_java_44(&[ins.clone()], 92_307_327);
        dump("Case F insertion", &[ins]);
        eprintln!(
            "Case F overlap @326 n={} @327 n={}",
            at_start.len(),
            at_plus_one.len()
        );
        assert_eq!(at_start.len(), 1, "insertion overlaps its start");
        assert_eq!(
            at_plus_one.len(),
            0,
            "Java insertion does not overlap start+1"
        );
    }

    /// Pre-6R.19 alt-span END (`start + max(|REF|,|ALT|) − 1`) vs Java/production REF-span END.
    #[test]
    fn six_r18_insertion_end_rust_span_vs_java() {
        let java = java_44_vc("2", 92_307_326, "C", "CAT");
        let legacy_span = VariationEvent {
            contig: "2".into(),
            start_1based: GenomePosition::new_1based(92_307_326),
            end_1based: GenomePosition::new_1based(
                92_307_326 + "CAT".len().max("C".len()).saturating_sub(1) as u64,
            ),
            ref_allele: "C".into(),
            alt_allele: "CAT".into(),
        };
        assert_eq!(java.end_1based.get(), 92_307_326);
        assert_eq!(
            legacy_span.end_1based.get(),
            92_307_328,
            "pre-6R.19 EventMap used max(|REF|,|ALT|)-1"
        );
        assert_eq!(
            get_overlapping_events_java_44(&[legacy_span.clone()], 92_307_327).len(),
            1,
            "legacy alt-span CAT overlaps 327; production/Java insertion does not"
        );
        assert_eq!(get_overlapping_events_java_44(&[java], 92_307_327).len(), 0);
        let _ = legacy_span;
    }

    /// `prefer_indel_over_colocated_snps` (union lists) still drops the SNP; EventMap uses makeBlock.
    #[test]
    fn six_r18_current_rust_prefer_indel_vs_java_makeblock() {
        let snp = java_44_vc("2", 92_307_326, "C", "T");
        let ins = java_44_vc("2", 92_307_326, "C", "CAT");
        let ag = java_44_vc("2", 92_307_327, "A", "G");
        let mut rust = vec![snp.clone(), ins.clone(), ag.clone()];
        prefer_indel_over_colocated_snps(&mut rust);
        dump("CURRENT RUST after prefer_indel", &rust);
        assert_eq!(rust.len(), 2, "SNP dropped, CAT and A→G kept");
        assert!(rust
            .iter()
            .any(|e| e.ref_allele == "C" && e.alt_allele == "CAT"));
        assert!(rust
            .iter()
            .any(|e| e.ref_allele == "A" && e.alt_allele == "G"));
        assert!(!rust
            .iter()
            .any(|e| e.ref_allele == "C" && e.alt_allele == "T"));
        assert!(
            !rust
                .iter()
                .any(|e| e.ref_allele == "C" && e.alt_allele == "TAT"),
            "prefer_indel (union lists) does not compose C→TAT"
        );

        let java = add_vc_merge_java_44(&[snp, ins, ag]).expect("java merge");
        dump("JAVA-STYLE addVC", &java);
        assert!(java
            .iter()
            .any(|e| e.ref_allele == "C" && e.alt_allele == "TAT"));
        assert!(!java.iter().any(|e| e.alt_allele == "CAT"));
        assert!(!java
            .iter()
            .any(|e| e.alt_allele == "T" && e.ref_allele == "C"));
    }

    /// P12 proposed EventMap (CIGAR walk order) vs Java addVC vs current prefer_indel.
    ///
    /// Same-start pair is C→T + C→CAT @ 92307326, **not** C→CAT + A→G (different starts).
    #[test]
    fn six_r18_p12_control_addvc_vs_prefer_indel() {
        let ttt = java_44_vc("2", 92_307_323, "TTT", "T");
        let ct = java_44_vc("2", 92_307_326, "C", "T");
        let cat = java_44_vc("2", 92_307_326, "C", "CAT");
        let ag = java_44_vc("2", 92_307_327, "A", "G");
        let proposed = [ttt.clone(), ct.clone(), cat.clone(), ag.clone()];
        dump("P12 INPUT (CIGAR-proposed, Java coords)", &proposed);

        let java = add_vc_merge_java_44(&proposed).expect("P12 java merge");
        dump("P12 JAVA-STYLE RESULT", &java);

        let mut rust = proposed.to_vec();
        prefer_indel_over_colocated_snps(&mut rust);
        dump("P12 CURRENT RUST RESULT (prefer_indel)", &rust);

        eprintln!("P12 DIFFERENCE:");
        eprintln!(
            "  Java n={} alleles={:?}",
            java.len(),
            java.iter()
                .map(|e| format!(
                    "{}-{} {}→{}",
                    e.start_1based.get(),
                    e.end_1based.get(),
                    e.ref_allele,
                    e.alt_allele
                ))
                .collect::<Vec<_>>()
        );
        eprintln!(
            "  Rust n={} alleles={:?}",
            rust.len(),
            rust.iter()
                .map(|e| format!(
                    "{}-{} {}→{}",
                    e.start_1based.get(),
                    e.end_1based.get(),
                    e.ref_allele,
                    e.alt_allele
                ))
                .collect::<Vec<_>>()
        );
        eprintln!(
            "  prefer_indel still drops C→T; production EventMap add_vc_merge matches Java C→TAT"
        );

        assert_eq!(java.len(), 3);
        assert_vc_eq(&java[0], &ttt, "P12 TTT unchanged");
        assert_eq!(java[1].ref_allele, "C");
        assert_eq!(java[1].alt_allele, "TAT");
        assert_eq!(java[1].start_1based.get(), 92_307_326);
        assert_eq!(java[1].end_1based.get(), 92_307_326);
        assert_vc_eq(&java[2], &ag, "P12 A→G unchanged (different start)");

        // Feeding only post-prefer_indel events (no C→T) does NOT invent C→TAT.
        let post_prefer = [ttt.clone(), cat.clone(), ag.clone()];
        let no_snp = add_vc_merge_java_44(&post_prefer).expect("no snp to merge");
        dump("P12 addVC on post-prefer_indel set (no C→T)", &no_snp);
        assert_eq!(no_snp.len(), 3);
        assert!(
            no_snp.iter().any(|e| e.alt_allele == "CAT"),
            "without the colocated SNP, Java addVC cannot compose TAT"
        );
        assert!(
            !no_snp.iter().any(|e| e.alt_allele == "TAT"),
            "C→CAT + A→G are not same-start; makeBlock does not apply"
        );
    }

    /// Java getOverlappingEvents: deletion ending at loc + insertion at loc keeps insertion.
    #[test]
    fn six_r18_overlapping_events_del_end_plus_insertion() {
        // EventMapUnitTest overlappingEvents hap1 at query loc 13.
        let del = java_44_vc(CHR, 10, "ACGG", "A");
        let ins = java_44_vc(CHR, 13, "G", "GTT");
        assert_eq!(del.end_1based.get(), 13);
        assert_eq!(ins.start_1based.get(), ins.end_1based.get());
        let at_13 = get_overlapping_events_java_44(&[del.clone(), ins.clone()], 13);
        dump("overlap loc=13", &at_13);
        assert_eq!(at_13.len(), 1);
        assert_eq!(at_13[0].ref_allele, "G");
        assert_eq!(at_13[0].alt_allele, "GTT");
        let at_12 = get_overlapping_events_java_44(&[del, ins], 12);
        assert_eq!(at_12.len(), 1);
        assert_eq!(at_12[0].ref_allele, "ACGG");
    }
}
