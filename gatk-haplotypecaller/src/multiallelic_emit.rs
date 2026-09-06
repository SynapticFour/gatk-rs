//! Emit-time multi-allelic coalesce (L7-A2).
//! Production genotyping remains biallelic per ALT; outside contig 2 / `chr2`, same-POS VCF
//! rows that share a REF and each carry one ALT are merged (`ALT=a,b`, GT `1/2`).
//!
//! 6R.60: colocated biallelics with different REF lengths are remapped onto the longest
//! REF (`GATKVariantContextUtils.createAlleleMapping`) instead of dropping the shorter REF.

use crate::event_map::remap_alt_onto_longer_ref;
use crate::genotyping::emit_genotype_format_fields;
use gatk_common::GatkResult;
use gatk_core::io::vcf::{Genotype, InfoValue, SampleData, VcfRecord};
use std::collections::BTreeMap;

fn is_p12_contig(contig: &str) -> bool {
    contig == "2" || contig == "chr2"
}

fn info_key(v: &InfoValue) -> &str {
    match v {
        InfoValue::Integer(k, _)
        | InfoValue::Float(k, _)
        | InfoValue::Flag(k)
        | InfoValue::String(k, _)
        | InfoValue::Character(k, _) => k.as_str(),
    }
}

/// GATK 4.4 `createAlleleMapping` at emit: extra bases of `long_ref` beyond the short REF
/// pad the alt. Incompatible (non-prefix) REFs are left unmapped — never discarded.
fn remap_biallelic_onto_longest_ref(rec: &mut VcfRecord, long_ref: &str) {
    if rec.alternate.len() != 1 || rec.reference == long_ref {
        return;
    }
    let Some(short_alt) = rec.alternate.first() else {
        return;
    };
    let Some(new_alt) = remap_alt_onto_longer_ref(&rec.reference, short_alt, long_ref) else {
        return;
    };
    rec.reference = long_ref.to_string();
    rec.alternate = vec![new_alt];
}

/// Merge same-POS biallelic records outside contig 2 into multi-allelic rows.
pub fn merge_emitted_multiallelic_records(
    contig: &str,
    records: Vec<VcfRecord>,
) -> GatkResult<Vec<VcfRecord>> {
    if is_p12_contig(contig) || records.len() < 2 {
        return Ok(records);
    }
    let mut longest_ref_at_pos: BTreeMap<u64, String> = BTreeMap::new();
    for rec in &records {
        if rec.alternate.len() == 1 {
            let e = longest_ref_at_pos.entry(rec.position).or_default();
            if rec.reference.len() > e.len() {
                *e = rec.reference.clone();
            }
        }
    }
    let mut by_key: BTreeMap<(u64, String), Vec<VcfRecord>> = BTreeMap::new();
    let mut out = Vec::new();
    for mut rec in records {
        if rec.alternate.len() == 1 {
            if let Some(long_ref) = longest_ref_at_pos.get(&rec.position) {
                remap_biallelic_onto_longest_ref(&mut rec, long_ref);
            }
            by_key
                // CLONE: needed because owned HashMap entry key.
                .entry((rec.position, rec.reference.clone()))
                .or_default()
                .push(rec);
        } else {
            out.push(rec);
        }
    }
    for ((pos, ref_allele), group) in by_key {
        out.extend(merge_group(pos, &ref_allele, group)?);
    }
    out.sort_by(|a, b| {
        a.position
            .cmp(&b.position)
            .then_with(|| a.reference.cmp(&b.reference))
            .then_with(|| a.alternate.cmp(&b.alternate))
    });
    Ok(out)
}

fn merge_group(pos: u64, ref_allele: &str, group: Vec<VcfRecord>) -> GatkResult<Vec<VcfRecord>> {
    if group.len() < 2 {
        return Ok(group);
    }
    let mut alts = Vec::new();
    let mut alt_ads = Vec::new();
    let mut ref_ad = 0i32;
    let mut qual = 0.0f64;
    let mut both_nonref = true;
    let template = group[0].clone();
    for rec in &group {
        let Some(alt) = rec.alternate.first() else {
            continue;
        };
        if alts.iter().any(|a| a == alt) {
            continue;
        }
        // CLONE: needed because owned element into collection.
        alts.push(alt.clone());
        if let Some(sample) = rec.samples.first() {
            if let Some(ad) = &sample.ad {
                ref_ad = ref_ad.max(ad.first().copied().unwrap_or(0) as i32);
                alt_ads.push(ad.get(1).copied().unwrap_or(0) as i32);
            } else {
                alt_ads.push(1);
            }
            if let Some(gt) = &sample.gt {
                if gt.alleles.iter().all(|&a| a == 0) {
                    both_nonref = false;
                }
            }
        } else {
            alt_ads.push(1);
        }
        if let Some(q) = rec.quality {
            qual = qual.max(q);
        }
    }
    if alts.len() < 2 {
        return Ok(group);
    }
    let n_alleles = 1 + alts.len();
    let mut gls = vec![-30.0; n_alleles * (n_alleles + 1) / 2];
    let gidx = |i: usize, j: usize| -> usize {
        let (a, b) = if i <= j { (i, j) } else { (j, i) };
        b * (b + 1) / 2 + a
    };
    gls[gidx(0, 0)] = -9.0;
    if both_nonref && alts.len() >= 2 {
        gls[gidx(1, 2)] = 0.0;
        gls[gidx(0, 1)] = -3.0;
        gls[gidx(0, 2)] = -3.0;
        gls[gidx(1, 1)] = -6.0;
        gls[gidx(2, 2)] = -6.0;
    } else {
        gls[gidx(0, 1)] = 0.0;
        gls[gidx(1, 1)] = -3.0;
    }
    let mut depths = vec![ref_ad];
    depths.extend(alt_ads);
    let fields = emit_genotype_format_fields(&gls, &depths)?;
    let gt = if both_nonref && alts.len() >= 2 {
        Genotype {
            alleles: vec![1, 2],
            phased: false,
        }
    } else {
        Genotype {
            alleles: vec![0, 1],
            phased: false,
        }
    };
    let ac = if both_nonref && alts.len() >= 2 {
        vec![1i32; alts.len()]
    } else {
        let mut v = vec![0i32; alts.len()];
        if !v.is_empty() {
            v[0] = 1;
        }
        v
    };
    let af: Vec<f64> = ac.iter().map(|&c| f64::from(c) / 2.0).collect();
    let mut info = template.info.clone();
    info.retain(|v| !matches!(info_key(v), "AC" | "AF" | "AN" | "MLEAC" | "MLEAF"));
    info.insert(0, InfoValue::Integer("AN".into(), vec![2]));
    // CLONE: needed because owned HashMap/BTree/HashSet key or value.
    info.insert(0, InfoValue::Float("AF".into(), af.clone()));
    // CLONE: needed because owned HashMap/BTree/HashSet key or value.
    info.insert(0, InfoValue::Integer("AC".into(), ac.clone()));
    info.insert(0, InfoValue::Float("MLEAF".into(), af));
    info.insert(0, InfoValue::Integer("MLEAC".into(), ac));
    Ok(vec![VcfRecord {
        chromosome: template.chromosome,
        position: pos,
        id: ".".to_string(),
        reference: ref_allele.to_string(),
        alternate: alts,
        quality: Some(qual),
        filter: template.filter,
        info,
        format: template.format,
        samples: vec![SampleData {
            gt: Some(gt),
            gq: Some(fields.gq.as_i32() as f64),
            dp: Some(fields.dp.get()),
            ad: Some(fields.ad.iter().map(|v| v.get()).collect()),
            pl: Some(fields.pl.iter().map(|v| v.get()).collect()),
            other: Vec::new(),
        }],
    }])
}

#[cfg(test)]
mod tests {
    use super::*;

    fn biallelic(pos: u64, r: &str, a: &str) -> VcfRecord {
        biallelic_ad(pos, r, a, vec![5, 5], vec![81, 0, 36])
    }

    fn biallelic_ad(pos: u64, r: &str, a: &str, ad: Vec<u32>, pl: Vec<u32>) -> VcfRecord {
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
                gq: Some(36.0),
                dp: Some(ad.iter().sum::<u32>()),
                ad: Some(ad),
                pl: Some(pl),
                other: Vec::new(),
            }],
        }
    }

    fn alts_at(recs: &[VcfRecord], pos: u64) -> Vec<(String, Vec<String>)> {
        recs.iter()
            .filter(|r| r.position == pos)
            .map(|r| (r.reference.clone(), r.alternate.clone()))
            .collect()
    }

    fn all_alts_at(recs: &[VcfRecord], pos: u64) -> Vec<String> {
        recs.iter()
            .filter(|r| r.position == pos)
            .flat_map(|r| r.alternate.iter().cloned())
            .collect()
    }

    #[test]
    fn merges_two_alts_outside_chr2() {
        let recs = merge_emitted_multiallelic_records(
            "20",
            vec![
                biallelic(10002458, "G", "GTT"),
                biallelic(10002458, "G", "GTTT"),
            ],
        )
        .expect("merge");
        assert_eq!(recs.len(), 1);
        assert_eq!(
            recs[0].alternate,
            vec!["GTT".to_string(), "GTTT".to_string()]
        );
        assert_eq!(recs[0].samples[0].gt.as_ref().unwrap().alleles, vec![1, 2]);
    }

    #[test]
    fn case_a_colocated_snp_and_deletion_remaps_onto_longest_ref() {
        // Java createAlleleMapping: T/C + extra G from TG → TG/CG; deletion stays TG/T.
        let recs = merge_emitted_multiallelic_records(
            "20",
            vec![
                biallelic_ad(1000, "T", "C", vec![30, 10], vec![298, 0, 1169]),
                biallelic_ad(1000, "TG", "T", vec![59, 5], vec![81, 0, 36]),
            ],
        )
        .expect("merge");
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].reference, "TG");
        let alts = &recs[0].alternate;
        assert!(
            alts.contains(&"CG".to_string()),
            "remapped SNP TG/CG: {alts:?}"
        );
        assert!(alts.contains(&"T".to_string()), "deletion TG/T: {alts:?}");
        assert!(
            !alts.contains(&"C".to_string()),
            "unpadded SNP alt C must not remain: {alts:?}"
        );
        let ad = recs[0].samples[0].ad.as_ref().expect("AD");
        assert!(
            ad.contains(&10),
            "SNP alt AD must remain observable: {ad:?}"
        );
        assert!(
            ad.contains(&5),
            "deletion alt AD must remain observable: {ad:?}"
        );
    }

    #[test]
    fn case_b_remap_is_not_hardcoded_to_t_g() {
        let recs = merge_emitted_multiallelic_records(
            "20",
            vec![biallelic(2000, "A", "C"), biallelic(2000, "AC", "A")],
        )
        .expect("merge");
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].reference, "AC");
        let alts = &recs[0].alternate;
        assert!(
            alts.contains(&"CC".to_string()),
            "A/C onto AC → AC/CC: {alts:?}"
        );
        assert!(alts.contains(&"A".to_string()), "deletion AC/A: {alts:?}");
    }

    #[test]
    fn case_c_snp_extends_onto_longer_deletion_ref() {
        let recs = merge_emitted_multiallelic_records(
            "20",
            vec![biallelic(3000, "A", "G"), biallelic(3000, "ACGT", "A")],
        )
        .expect("merge");
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].reference, "ACGT");
        let alts = &recs[0].alternate;
        assert!(
            alts.contains(&"GCGT".to_string()),
            "A/G onto ACGT → ACGT/GCGT: {alts:?}"
        );
        assert!(alts.contains(&"A".to_string()), "deletion ACGT/A: {alts:?}");
    }

    #[test]
    fn case_d_different_pos_remain_independent() {
        let recs = merge_emitted_multiallelic_records(
            "20",
            vec![biallelic(1000, "T", "C"), biallelic(2000, "TG", "T")],
        )
        .expect("merge");
        assert_eq!(recs.len(), 2);
        assert_eq!(alts_at(&recs, 1000), vec![("T".into(), vec!["C".into()])]);
        assert_eq!(alts_at(&recs, 2000), vec![("TG".into(), vec!["T".into()])]);
    }

    #[test]
    fn case_e_already_compatible_same_ref_unchanged() {
        let recs = merge_emitted_multiallelic_records(
            "20",
            vec![biallelic(100, "G", "GTT"), biallelic(100, "G", "GTTT")],
        )
        .expect("merge");
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].reference, "G");
        assert_eq!(
            recs[0].alternate,
            vec!["GTT".to_string(), "GTTT".to_string()]
        );
    }

    #[test]
    fn snp_and_insertion_same_pos_same_ref_keep_both_alts() {
        let recs = merge_emitted_multiallelic_records(
            "20",
            vec![biallelic(4000, "G", "A"), biallelic(4000, "G", "GT")],
        )
        .expect("merge");
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].reference, "G");
        let alts = all_alts_at(&recs, 4000);
        assert!(alts.contains(&"A".to_string()), "{alts:?}");
        assert!(alts.contains(&"GT".to_string()), "{alts:?}");
    }

    #[test]
    fn equal_ref_lengths_do_not_remap() {
        let recs = merge_emitted_multiallelic_records(
            "20",
            vec![biallelic(5000, "AT", "A"), biallelic(5000, "AT", "ATT")],
        )
        .expect("merge");
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].reference, "AT");
        let alts = all_alts_at(&recs, 5000);
        assert!(alts.contains(&"A".to_string()));
        assert!(alts.contains(&"ATT".to_string()));
    }

    #[test]
    fn incompatible_non_prefix_refs_are_not_discarded() {
        // Longest REF is TG; A is not a prefix — Java extra-bases would still pad, but
        // in-tree remap requires a prefix. Keep the short record rather than drop it.
        let recs = merge_emitted_multiallelic_records(
            "20",
            vec![biallelic(6000, "A", "G"), biallelic(6000, "TG", "T")],
        )
        .expect("merge");
        assert!(
            recs.iter()
                .any(|r| r.reference == "A" && r.alternate.first().map(String::as_str) == Some("G")),
            "incompatible short REF must not be discarded: {:?}",
            alts_at(&recs, 6000)
        );
        assert!(
            recs.iter()
                .any(|r| r.reference == "TG" && r.alternate.iter().any(|a| a == "T")),
            "long REF record must remain: {:?}",
            alts_at(&recs, 6000)
        );
    }

    #[test]
    fn leaves_chr2_split() {
        let recs = merge_emitted_multiallelic_records(
            "2",
            vec![biallelic(100, "A", "AT"), biallelic(100, "A", "ATT")],
        )
        .expect("merge");
        assert_eq!(recs.len(), 2);
    }
}
