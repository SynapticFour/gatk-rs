//! Emit-time multi-allelic coalesce (L7-A2).
//! Production genotyping remains biallelic per ALT; outside contig 2 / `chr2`, same-POS VCF
//! rows that share a REF and each carry one ALT are merged (`ALT=a,b`, GT `1/2`).

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

/// Merge same-POS biallelic records outside contig 2 into multi-allelic rows.
pub fn merge_emitted_multiallelic_records(
    contig: &str,
    records: Vec<VcfRecord>,
) -> GatkResult<Vec<VcfRecord>> {
    if is_p12_contig(contig) || records.len() < 2 {
        return Ok(records);
    }
    // L10: when the same POS has multiple REF lengths (nested STR), keep the longest REF
    // group only — shorter representations are fragments of the Java multi-allelic site.
    let mut longest_ref_at_pos: BTreeMap<u64, usize> = BTreeMap::new();
    for rec in &records {
        if rec.alternate.len() == 1 {
            let e = longest_ref_at_pos.entry(rec.position).or_insert(0);
            *e = (*e).max(rec.reference.len());
        }
    }
    let mut by_key: BTreeMap<(u64, String), Vec<VcfRecord>> = BTreeMap::new();
    let mut out = Vec::new();
    for rec in records {
        if rec.alternate.len() == 1 {
            if let Some(&max_len) = longest_ref_at_pos.get(&rec.position) {
                if rec.reference.len() < max_len {
                    continue;
                }
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
                dp: Some(10),
                ad: Some(vec![5, 5]),
                pl: Some(vec![81, 0, 36]),
                other: Vec::new(),
            }],
        }
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
    fn leaves_chr2_split() {
        let recs = merge_emitted_multiallelic_records(
            "2",
            vec![biallelic(100, "A", "AT"), biallelic(100, "A", "ATT")],
        )
        .expect("merge");
        assert_eq!(recs.len(), 2);
    }
}
