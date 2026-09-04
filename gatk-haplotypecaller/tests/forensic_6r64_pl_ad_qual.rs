//! 6R.64 coordinate-free diagnostics: QUAL is a PL consequence; Java AD remarginalizes
//! after unused-ALT drop. No genomic coordinates. No production algorithm change.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r64_pl_ad_qual
//! ```

use gatk_haplotypecaller::genotyping::ReadLikelihoodRow;
use gatk_haplotypecaller::read_realignment::LOG_10_INFORMATIVE_THRESHOLD;
use gatk_haplotypecaller::variant_site_hc_annotations::qual_from_af_calculation;
use gatk_haplotypecaller::InformativeAd;

fn row(lls: Vec<f64>) -> ReadLikelihoodRow {
    ReadLikelihoodRow {
        read_index: 0,
        read_id: String::new(),
        haplotype_log10_likelihoods: lls,
    }
}

fn informative_best(lls: &[f64]) -> Option<usize> {
    let mut best_i = 0usize;
    let mut best = f64::NEG_INFINITY;
    let mut second = f64::NEG_INFINITY;
    for (i, &ll) in lls.iter().enumerate() {
        if ll > best {
            second = best;
            best = ll;
            best_i = i;
        } else if ll > second {
            second = ll;
        }
    }
    if best.is_finite() && (best - second).abs() > LOG_10_INFORMATIVE_THRESHOLD {
        Some(best_i)
    } else {
        None
    }
}

#[test]
fn qual_from_java_emitted_pl_is_not_java_vcf_qual() {
    let gl = [-54.2, 0.0, -135.3];
    let qual = qual_from_af_calculation(&gl).expect("qual");
    assert!(
        (qual - 534.64).abs() < 1.0,
        "Rust AF on Java emitted PL 542,0,1353 → ~534.64, got {qual}"
    );
    assert!(
        (qual - 510.06).abs() > 10.0,
        "Java site QUAL 510.06 is not AF(emitted PL); got {qual}"
    );
}

#[test]
fn qual_from_rust_pl_298_0_1103_matches_emitted() {
    let gl = [-29.8, 0.0, -110.3];
    let qual = qual_from_af_calculation(&gl).expect("qual");
    assert!(
        (qual - 290.64).abs() < 1.0,
        "Rust PL 298,0,1103 must produce QUAL near 290.64, got {qual}"
    );
}

#[test]
fn remarginalize_after_dropping_unused_alt_can_increase_ad() {
    // Alleles: 0=REF, 1=deletion T, 2=CG. Java DepthPerAlleleBySample remarginalizes
    // onto remaining alleles after unused-ALT subset; Rust permutes 3-way counts.
    let rows = [
        row(vec![0.0, -10.0, -10.0]), // informative REF in 3-way
        row(vec![-10.0, -10.0, 0.0]), // informative CG in 3-way
        row(vec![-10.0, 0.0, -10.0]), // informative deletion in 3-way
        row(vec![0.0, -0.05, -10.0]), // 3-way uninformative (REF vs T); 2-way REF vs CG
    ];
    let mut ad3 = [0i32; 3];
    for r in &rows {
        if let Some(i) = informative_best(&r.haplotype_log10_likelihoods) {
            ad3[i] += 1;
        }
    }
    assert_eq!(ad3, [1, 1, 1], "3-way informative AD");
    let permuted = vec![ad3[0], ad3[2]];
    assert_eq!(permuted, vec![1, 1]);

    let two_way: Vec<ReadLikelihoodRow> = rows
        .iter()
        .map(|r| {
            row(vec![
                r.haplotype_log10_likelihoods[0],
                r.haplotype_log10_likelihoods[2],
            ])
        })
        .collect();
    let remarg = InformativeAd::from_marginalized_rows(&two_way, 0, 1, None);
    assert_eq!(
        remarg.as_vec(),
        vec![2, 1],
        "near-tie REF-vs-T becomes informative REF vs CG; deletion-only read is 2-way uninformative"
    );
    assert_ne!(
        remarg.as_vec(),
        permuted,
        "Java annotation AD (remarginalize) ≠ Rust permute-after-3-way"
    );
}
