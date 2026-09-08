//! 6R.61 diagnostic: colocated alleles merge **before** genotype calculation.
//!
//! Coordinate-free. Extends the 6R.59/6R.60 merge-order tests: genotype input must be
//! `TG/T,CG` (or another REF/ALT pair) and GLs must be calculated in that 3-allele space
//! — not a SNP-only PL later remapped at emit.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test java_merge_order_colocated_snp_and_deletion -- --nocapture
//! ```

use gatk_core::io::vcf::{Genotype, SampleData, VcfRecord};
use gatk_haplotypecaller::event_map::{
    is_colocated_snp_indel_merged_site, merged_alleles_for_genotyping,
    merged_biallelic_sites_at_position, remap_alt_onto_longer_ref, VariationEvent,
};
use gatk_haplotypecaller::genotyping::ReadLikelihoodRow;
use gatk_haplotypecaller::hc_genotyping_engine::diploid_genotype_log10_likelihoods_from_allele_rows;
use gatk_haplotypecaller::multiallelic_emit::merge_emitted_multiallelic_records;

/// Synthetic POS only — not a genomic pin.
const POS: u64 = 1000;

fn biallelic_vcf(pos: u64, r: &str, a: &str) -> VcfRecord {
    biallelic_vcf_ad(pos, r, a, vec![30, 10], vec![298, 0, 1169])
}

fn biallelic_vcf_ad(pos: u64, r: &str, a: &str, ad: Vec<u32>, pl: Vec<u32>) -> VcfRecord {
    VcfRecord {
        chromosome: "20".into(),
        position: pos,
        id: ".".into(),
        reference: r.into(),
        alternate: vec![a.into()],
        quality: Some(100.0),
        filter: vec![".".into()],
        info: vec![],
        format: vec![
            "GT".into(),
            "AD".into(),
            "DP".into(),
            "GQ".into(),
            "PL".into(),
        ],
        samples: vec![SampleData {
            gt: Some(Genotype {
                alleles: vec![0, 1],
                phased: false,
            }),
            gq: Some(99.0),
            dp: Some(ad.iter().sum::<u32>()),
            ad: Some(ad),
            pl: Some(pl),
            other: Vec::new(),
        }],
    }
}

fn strip_common_suffix(ref_a: &str, alt_a: &str) -> (String, String) {
    let mut r: Vec<char> = ref_a.chars().collect();
    let mut a: Vec<char> = alt_a.chars().collect();
    while r.len() > 1 && a.len() > 1 && r.last() == a.last() {
        r.pop();
        a.pop();
    }
    (r.into_iter().collect(), a.into_iter().collect())
}

fn emitted_alts(recs: &[VcfRecord], pos: u64) -> (Option<String>, Vec<String>) {
    let rec = recs.iter().find(|r| r.position == pos);
    match rec {
        Some(r) => (Some(r.reference.clone()), r.alternate.clone()),
        None => (None, Vec::new()),
    }
}

#[test]
fn java_merge_order_colocated_snp_and_deletion() {
    let snp = VariationEvent::from_alleles("20", POS, "T", "C");
    let del = VariationEvent::from_alleles("20", POS, "TG", "T");

    // Java `GATKVariantContextUtils.createAlleleMapping` / `simpleMerge`:
    // longest REF is TG; SNP alt C is extended by extra base G → CG.
    let remapped = remap_alt_onto_longer_ref("T", "C", "TG");
    let merged = merged_biallelic_sites_at_position(&[snp.clone(), del.clone()], POS);
    let mut merged_keys: Vec<(String, String)> = merged
        .iter()
        .map(|e| (e.ref_allele.clone(), e.alt_allele.clone()))
        .collect();
    merged_keys.sort();

    // Java `reverseTrimAlleles` when unused deletion alt is dropped: TG/CG → T/C.
    let (trim_ref, trim_alt) = match remapped.as_deref() {
        Some(alt) => strip_common_suffix("TG", alt),
        None => (String::new(), String::new()),
    };

    let emitted = merge_emitted_multiallelic_records(
        "20",
        vec![
            biallelic_vcf_ad(POS, "T", "C", vec![30, 10], vec![298, 0, 1169]),
            biallelic_vcf_ad(POS, "TG", "T", vec![59, 5], vec![81, 0, 36]),
        ],
    )
    .expect("merge");

    eprintln!("6R.60 intermediates:");
    eprintln!("  per-haplotype (input): T/C and TG/T at POS={POS}");
    eprintln!("  Java createAlleleMapping analogue: T/C onto TG → {remapped:?}");
    eprintln!("  merged_biallelic_sites_at_position (pre-genotype analogue): {merged_keys:?}");
    eprintln!("  reverseTrim analogue of TG/{{remap}}: {trim_ref}/{trim_alt}");
    eprintln!(
        "  merge_emitted_multiallelic_records (post-genotype, remapped): {:?}",
        emitted
            .iter()
            .map(|r| format!(
                "{} {}/{:?} AD={:?}",
                r.position,
                r.reference,
                r.alternate,
                r.samples.first().and_then(|s| s.ad.clone())
            ))
            .collect::<Vec<_>>()
    );

    assert_eq!(
        remapped.as_deref(),
        Some("CG"),
        "Java createAlleleMapping: extra bases of TG beyond T pad SNP alt C → CG"
    );
    assert!(
        merged_keys.iter().any(|(r, a)| r == "TG" && a == "CG"),
        "Java simpleMerge analogue keeps remapped SNP TG/CG"
    );
    assert!(
        merged_keys.iter().any(|(r, a)| r == "TG" && a == "T"),
        "Java simpleMerge analogue keeps deletion TG/T"
    );
    assert_eq!(
        (trim_ref.as_str(), trim_alt.as_str()),
        ("T", "C"),
        "Java reverseTrimAlleles of TG/CG (after unused-alt drop) is T/C"
    );

    let (eref, ealts) = emitted_alts(&emitted, POS);
    assert_eq!(eref.as_deref(), Some("TG"));
    assert!(
        ealts.contains(&"CG".to_string()),
        "SNP must remap to TG/CG, not disappear: {ealts:?}"
    );
    assert!(
        ealts.contains(&"T".to_string()),
        "deletion TG/T must remain: {ealts:?}"
    );
    let ad = emitted[0].samples[0].ad.as_ref().expect("AD");
    assert!(ad.contains(&10), "SNP AD must remain observable: {ad:?}");
}

#[test]
fn case_b_remap_not_hardcoded_to_t_g() {
    let emitted = merge_emitted_multiallelic_records(
        "20",
        vec![
            biallelic_vcf(2000, "A", "C"),
            biallelic_vcf(2000, "AC", "A"),
        ],
    )
    .expect("merge");
    let (eref, ealts) = emitted_alts(&emitted, 2000);
    assert_eq!(eref.as_deref(), Some("AC"));
    assert!(ealts.contains(&"CC".to_string()), "{ealts:?}");
    assert!(ealts.contains(&"A".to_string()), "{ealts:?}");
}

#[test]
fn case_c_snp_extends_onto_longer_ref() {
    let emitted = merge_emitted_multiallelic_records(
        "20",
        vec![
            biallelic_vcf(3000, "A", "G"),
            biallelic_vcf(3000, "ACGT", "A"),
        ],
    )
    .expect("merge");
    let (eref, ealts) = emitted_alts(&emitted, 3000);
    assert_eq!(eref.as_deref(), Some("ACGT"));
    assert!(ealts.contains(&"GCGT".to_string()), "{ealts:?}");
    assert!(ealts.contains(&"A".to_string()), "{ealts:?}");
}

#[test]
fn case_d_non_colocated_records_stay_independent() {
    let emitted = merge_emitted_multiallelic_records(
        "20",
        vec![
            biallelic_vcf(1000, "T", "C"),
            biallelic_vcf(2000, "TG", "T"),
        ],
    )
    .expect("merge");
    assert_eq!(emitted.len(), 2);
    let (r1, a1) = emitted_alts(&emitted, 1000);
    let (r2, a2) = emitted_alts(&emitted, 2000);
    assert_eq!(r1.as_deref(), Some("T"));
    assert_eq!(a1, vec!["C".to_string()]);
    assert_eq!(r2.as_deref(), Some("TG"));
    assert_eq!(a2, vec!["T".to_string()]);
}

#[test]
fn case_e_already_compatible_alleles_unchanged() {
    let emitted = merge_emitted_multiallelic_records(
        "20",
        vec![
            biallelic_vcf(100, "G", "GTT"),
            biallelic_vcf(100, "G", "GTTT"),
        ],
    )
    .expect("merge");
    let (eref, ealts) = emitted_alts(&emitted, 100);
    assert_eq!(eref.as_deref(), Some("G"));
    assert_eq!(ealts, vec!["GTT".to_string(), "GTTT".to_string()]);
}

#[test]
fn merge_colocated_alleles_before_genotyping() {
    let snp = VariationEvent::from_alleles("20", POS, "T", "C");
    let del = VariationEvent::from_alleles("20", POS, "TG", "T");
    let (long_ref, alts) =
        merged_alleles_for_genotyping(&[snp, del], POS).expect("pre-genotype merge");
    assert_eq!(long_ref, "TG");
    assert_eq!(
        alts,
        vec!["T".to_string(), "CG".to_string()],
        "genotype input must be TG/T,CG — not independent T/C"
    );
    assert!(is_colocated_snp_indel_merged_site(&long_ref, &alts));

    // Joint GLs from 3-allele marginalized rows (REF, T, CG). Evidence supports both ALTs
    // so the best GT is 1/2 from the 6-state calculator — not a SNP-only 3-PL remapped
    // to fabricated emit-merge PL 90,30,60,30,0,60.
    let rows = vec![
        ReadLikelihoodRow {
            read_index: 0,
            read_id: "r0".into(),
            haplotype_log10_likelihoods: vec![-10.0, -0.1, -8.0],
        },
        ReadLikelihoodRow {
            read_index: 1,
            read_id: "r1".into(),
            haplotype_log10_likelihoods: vec![-10.0, -8.0, -0.1],
        },
        ReadLikelihoodRow {
            read_index: 2,
            read_id: "r2".into(),
            haplotype_log10_likelihoods: vec![-9.0, -0.2, -0.2],
        },
        ReadLikelihoodRow {
            read_index: 3,
            read_id: "r3".into(),
            haplotype_log10_likelihoods: vec![-9.0, -0.15, -7.0],
        },
        ReadLikelihoodRow {
            read_index: 4,
            read_id: "r4".into(),
            haplotype_log10_likelihoods: vec![-9.0, -7.0, -0.15],
        },
    ];
    let gls = diploid_genotype_log10_likelihoods_from_allele_rows(&rows, 3);
    assert_eq!(gls.len(), 6);
    // Independent T/C genotyping sees only REF + SNP-like ALT (columns 0 and 2) → 3 GLs.
    let snp_only: Vec<ReadLikelihoodRow> = rows
        .iter()
        .map(|r| ReadLikelihoodRow {
            read_index: r.read_index,
            read_id: r.read_id.clone(),
            haplotype_log10_likelihoods: vec![
                r.haplotype_log10_likelihoods[0],
                r.haplotype_log10_likelihoods[2],
            ],
        })
        .collect();
    let snp_gls = diploid_genotype_log10_likelihoods_from_allele_rows(&snp_only, 2);
    assert_eq!(snp_gls.len(), 3);
    assert_ne!(
        &gls[..3],
        snp_gls.as_slice(),
        "merged 6-GL must not be a padded copy of SNP-only 3-GL"
    );
    let fields = gatk_haplotypecaller::emit_genotype_format_fields(&gls, &[10, 2, 8]).expect("PL");
    let pl: Vec<i32> = fields.pl_as_i32();
    assert_ne!(pl, vec![90, 30, 60, 30, 0, 60]);
    assert_eq!(pl.iter().min().copied().unwrap_or(-1), 0);
    let best = gatk_haplotypecaller::best_pl_index(&fields.pl);
    let gt = gatk_haplotypecaller::diploid_genotype_alleles_from_pl_index(3, best);
    assert_eq!(
        gt,
        vec![1, 2],
        "GT must come from min PL in the 6-state merged list, got {gt:?} PL={pl:?}"
    );
}

#[test]
fn full_padded_eventmap_union_enters_six_state_gls_before_emit() {
    use gatk_haplotypecaller::cigar::{Cigar, CigarOperator};
    use gatk_haplotypecaller::event_map::collect_variation_events;
    use gatk_haplotypecaller::haplotype::Haplotype;

    fn hap_with_cigar(
        bases: &[u8],
        elements: &[(usize, CigarOperator)],
        is_ref: bool,
    ) -> Haplotype {
        let mut cigar = Cigar::new();
        for (len, op) in elements {
            cigar.push(*len, *op);
        }
        Haplotype {
            bases: bases.to_vec(),
            is_reference: is_ref,
            score: 0.0,
            kmer_size: 10,
            cigar: Some(cigar),
            genome_loc: None,
            alignment_start_hap_wrt_ref: 0,
        }
    }

    let mut full_ref = vec![b'A'; 10];
    full_ref.extend_from_slice(b"TG");
    full_ref.extend(std::iter::repeat_n(b'A', 12));
    let full_pad = 100u64;
    let loc = 110u64;
    let mut ref_hap = hap_with_cigar(&full_ref, &[(full_ref.len(), CigarOperator::Match)], true);
    ref_hap.is_reference = true;
    let mut snp_bases = full_ref.clone();
    snp_bases[10] = b'C';
    let snp_hap = hap_with_cigar(
        &snp_bases,
        &[(snp_bases.len(), CigarOperator::Match)],
        false,
    );
    let mut del_bases = full_ref.clone();
    del_bases.remove(11);
    let del_hap = hap_with_cigar(
        &del_bases,
        &[
            (11, CigarOperator::Match),
            (1, CigarOperator::Deletion),
            (12, CigarOperator::Match),
        ],
        false,
    );
    let events =
        collect_variation_events(&[ref_hap, snp_hap, del_hap], &full_ref, full_pad, "20", 0);
    let (long_ref, alts) =
        merged_alleles_for_genotyping(&events, loc).expect("EventMap union merge");
    assert_eq!(long_ref, "TG");
    assert_eq!(
        alts,
        vec!["T".to_string(), "CG".to_string()],
        "EventMap union must enter genotyping as TG/T,CG"
    );
    assert!(is_colocated_snp_indel_merged_site(&long_ref, &alts));

    let rows = vec![
        ReadLikelihoodRow {
            read_index: 0,
            read_id: "r0".into(),
            haplotype_log10_likelihoods: vec![-10.0, -0.1, -8.0],
        },
        ReadLikelihoodRow {
            read_index: 1,
            read_id: "r1".into(),
            haplotype_log10_likelihoods: vec![-10.0, -8.0, -0.1],
        },
    ];
    let gls = diploid_genotype_log10_likelihoods_from_allele_rows(&rows, 1 + alts.len());
    assert_eq!(
        gls.len(),
        6,
        "GLs must be calculated in merged 3-allele space before emit"
    );
}

#[test]
fn unused_alt_subset_remaps_genotype_after_merged_genotyping() {
    // Lifecycle: merge → genotype merged space → unused-ALT subset. Not: genotype T/C then remap.
    let snp = VariationEvent::from_alleles("chrX", POS, "T", "C");
    let del = VariationEvent::from_alleles("chrX", POS, "TG", "T");
    let (long_ref, alts) =
        merged_alleles_for_genotyping(&[snp, del], POS).expect("pre-genotype merge");
    assert_eq!(long_ref, "TG");
    assert_eq!(alts, vec!["T".to_string(), "CG".to_string()]);

    let gls = vec![-29.8, -33.7, -162.0, 0.0, -105.8, -110.3];
    let ad = vec![28, 2, 10];
    let fields = gatk_haplotypecaller::emit_genotype_format_fields(&gls, &ad).expect("PL");
    let best = gatk_haplotypecaller::best_pl_index(&fields.pl);
    let gt = gatk_haplotypecaller::diploid_genotype_alleles_from_pl_index(3, best);
    assert_eq!(gt, vec![0, 2], "merged GT before subset is 0/2");

    let subset =
        gatk_haplotypecaller::subset_unused_alts_after_merged_genotyping(&alts, &gt, &gls, &ad)
            .expect("subset");
    assert_eq!(subset.alt_alleles, vec!["CG".to_string()]);
    assert_eq!(subset.ad, vec![28, 10]);
    let out = gatk_haplotypecaller::emit_genotype_format_fields(&subset.log10_gls, &subset.ad)
        .expect("subset PL");
    assert_eq!(out.pl_as_i32().len(), 3);
    assert_ne!(out.pl_as_i32(), vec![90, 30, 60, 30, 0, 60]);
    let out_gt = gatk_haplotypecaller::diploid_genotype_alleles_from_pl_index(
        1 + subset.alt_alleles.len(),
        gatk_haplotypecaller::best_pl_index(&out.pl),
    );
    assert_eq!(
        out_gt,
        vec![0, 1],
        "0/2 on [TG,T,CG] remaps to 0/1 on [TG,CG]"
    );

    // Forbidden lifecycle: independent T/C genotyping is a 3-GL space, not the subset source.
    assert_ne!(
        subset.log10_gls.len(),
        gls.len(),
        "subset is applied after 6-state merged genotyping"
    );
}

#[test]
fn reverse_trim_common_suffix_after_unused_alt_subset() {
    // Lifecycle: merge → genotype merged space → unused-ALT subset → reverseTrimAlleles.
    // Forbidden: reverse-trim TG/T,CG first (length-1 T blocks suffix clip).
    let snp = VariationEvent::from_alleles("chrX", POS, "T", "C");
    let del = VariationEvent::from_alleles("chrX", POS, "TG", "T");
    let (long_ref, alts) =
        merged_alleles_for_genotyping(&[snp, del], POS).expect("pre-genotype merge");
    assert_eq!(long_ref, "TG");
    assert_eq!(alts, vec!["T".to_string(), "CG".to_string()]);

    let before_subset = gatk_haplotypecaller::reverse_trim_alleles(&long_ref, &alts);
    assert_eq!(
        before_subset,
        (long_ref.clone(), alts.clone()),
        "reverse-trim before unused-ALT subset is a no-op while length-1 T remains"
    );

    let gls = vec![-29.8, -33.7, -162.0, 0.0, -105.8, -110.3];
    let ad = vec![28, 2, 10];
    let subset =
        gatk_haplotypecaller::subset_unused_alts_after_merged_genotyping(&alts, &[0, 2], &gls, &ad)
            .expect("subset");
    assert_eq!(subset.alt_alleles, vec!["CG".to_string()]);
    assert_eq!(subset.ad, vec![28, 10]);
    let pl = gatk_haplotypecaller::emit_genotype_format_fields(&subset.log10_gls, &subset.ad)
        .expect("subset PL");
    assert_eq!(pl.pl_as_i32(), vec![298, 0, 1103]);
    let gt = gatk_haplotypecaller::diploid_genotype_alleles_from_pl_index(
        1 + subset.alt_alleles.len(),
        gatk_haplotypecaller::best_pl_index(&pl.pl),
    );
    assert_eq!(gt, vec![0, 1]);

    let (trim_ref, trim_alts) =
        gatk_haplotypecaller::reverse_trim_alleles(&long_ref, &subset.alt_alleles);
    assert_eq!(trim_ref, "T");
    assert_eq!(trim_alts, vec!["C".to_string()]);
    let after = gatk_haplotypecaller::emit_genotype_format_fields(&subset.log10_gls, &subset.ad)
        .expect("post-trim PL");
    assert_eq!(
        after.pl_as_i32(),
        vec![298, 0, 1103],
        "reverse-trim must not recalculate PL"
    );
    assert_eq!(subset.ad, vec![28, 10], "reverse-trim must not change AD");
    let after_gt = gatk_haplotypecaller::diploid_genotype_alleles_from_pl_index(
        1 + trim_alts.len(),
        gatk_haplotypecaller::best_pl_index(&after.pl),
    );
    assert_eq!(
        after_gt,
        vec![0, 1],
        "GT indices unchanged after reverse-trim"
    );
}

#[test]
fn unused_alt_subset_does_not_drop_used_second_alt() {
    let alts = vec!["T".to_string(), "CG".to_string()];
    let gls = vec![0.0, -1.0, -20.0, -30.0, -40.0, -60.0];
    let ad = vec![10, 8, 1];
    let r =
        gatk_haplotypecaller::subset_unused_alts_after_merged_genotyping(&alts, &[0, 1], &gls, &ad)
            .expect("subset");
    assert_eq!(r.alt_alleles, vec!["T".to_string()]);
}

#[test]
fn unused_alt_subset_keeps_both_alts_when_gt_is_1_2() {
    let alts = vec!["T".to_string(), "CG".to_string()];
    let gls = vec![-9.0, -3.0, -6.0, -3.0, 0.0, -6.0];
    let ad = vec![5, 10, 12];
    let r =
        gatk_haplotypecaller::subset_unused_alts_after_merged_genotyping(&alts, &[1, 2], &gls, &ad)
            .expect("subset");
    assert_eq!(r.alt_alleles, vec!["T".to_string(), "CG".to_string()]);
    assert_eq!(r.log10_gls, gls);
}

#[test]
fn merge_colocated_alleles_before_genotyping_other_bases() {
    let snp = VariationEvent::from_alleles("20", 2000, "A", "C");
    let del = VariationEvent::from_alleles("20", 2000, "AC", "A");
    let (long_ref, alts) =
        merged_alleles_for_genotyping(&[snp, del], 2000).expect("pre-genotype merge");
    assert_eq!(long_ref, "AC");
    assert_eq!(alts, vec!["A".to_string(), "CC".to_string()]);
}
