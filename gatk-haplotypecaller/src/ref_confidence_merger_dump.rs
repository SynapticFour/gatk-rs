//! H-D02 — `ReferenceConfidenceVariantContextMerger` parity dumps.

use crate::ref_confidence_merger::{
    merge_reference_confidence, MergeAllele, MergeGenotype, MergeVcInput, RefConfidenceMergeResult,
    NON_REF_ALLELE,
};
use gatk_common::GatkResult;
use std::io::Write;

/// Per-sample genotype row in H-D02 reference-confidence merge parity dumps.
/// # Invariants
/// Mirrors [`MergeGenotypeOut`](crate::ref_confidence_merger::MergeGenotypeOut) plus optional GQ for dump completeness.
/// # Ownership
/// Owns name and optional PL/AD/GQ vectors for TSV serialization.
/// # Mutation
/// Immutable dump snapshot.
/// # Biological assumptions
/// None — parity fixture/dump shape.
/// # Java equivalence
/// Rust-native dump of GATK `ReferenceConfidenceVariantContextMerger` unit-test outputs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DumpGenotype {
    pub name: String,
    pub pl: Option<Vec<i32>>,
    pub ad: Option<Vec<i32>>,
    pub gq: Option<i32>,
}

/// Full merged-site dump for one H-D02 parity case.
/// # Invariants
/// When `merged_null` is true, allele/genotype vectors are empty (Java null merge).
/// # Ownership
/// Owns case metadata, alleles, and genotype dump rows.
/// # Mutation
/// Immutable after fixture merge or null-merge construction.
/// # Biological assumptions
/// None — test/dump artifact.
/// # Java equivalence
/// GATK `ReferenceConfidenceVariantContextMergerUnitTest` expected merge snapshots.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RefConfidenceMergeDump {
    pub case_id: String,
    pub contig: String,
    pub pos: u64,
    pub alleles: Vec<String>,
    pub genotypes: Vec<DumpGenotype>,
    pub has_non_ref: bool,
    pub merged_null: bool,
}

fn a_ref(base: &str) -> MergeAllele {
    MergeAllele {
        bases: base.to_string(),
        is_reference: true,
    }
}

fn a_alt(base: &str) -> MergeAllele {
    MergeAllele {
        bases: base.to_string(),
        is_reference: false,
    }
}

fn a_non_ref() -> MergeAllele {
    MergeAllele {
        bases: NON_REF_ALLELE.to_string(),
        is_reference: false,
    }
}

fn g_pl(sample: &str, pl: &[i32]) -> MergeGenotype {
    MergeGenotype {
        sample: sample.to_string(),
        pl: Some(pl.to_vec()),
        ad: None,
    }
}

fn g_ad(sample: &str, ad: &[i32]) -> MergeGenotype {
    MergeGenotype {
        sample: sample.to_string(),
        pl: None,
        ad: Some(ad.to_vec()),
    }
}

fn g_pl_ad(sample: &str, pl: &[i32], ad: &[i32]) -> MergeGenotype {
    MergeGenotype {
        sample: sample.to_string(),
        pl: Some(pl.to_vec()),
        ad: Some(ad.to_vec()),
    }
}

fn vc(
    source: &str,
    start: u64,
    alleles: Vec<MergeAllele>,
    genotypes: Vec<MergeGenotype>,
) -> MergeVcInput {
    MergeVcInput {
        source: source.to_string(),
        start,
        alleles,
        genotypes,
    }
}

/// GATK `ReferenceConfidenceVariantContextMergerUnitTest` fixtures.
fn unit_test_vcs(case_id: &str) -> GatkResult<(String, u64, Vec<MergeVcInput>)> {
    let start = 10u64;
    let contig = "20".to_string();
    let standard_pls = vec![30, 20, 10, 71, 72, 73];
    let vcs = match case_id {
        "merge_single_vc" | "test00" => vec![vc(
            "test",
            start,
            vec![a_ref("A"), a_alt("C"), a_non_ref()],
            vec![g_pl("A_C", &standard_pls)],
        )],
        "merge_two_snps" | "test01" => vec![
            vc(
                "test",
                start,
                vec![a_ref("A"), a_alt("C"), a_non_ref()],
                vec![g_pl("A_C", &standard_pls)],
            ),
            vc(
                "test",
                start,
                vec![a_ref("A"), a_alt("G"), a_non_ref()],
                vec![g_pl("A_G", &standard_pls)],
            ),
        ],
        "merge_snp_indel" | "test02" => vec![
            vc(
                "test",
                start,
                vec![a_ref("A"), a_alt("C"), a_non_ref()],
                vec![g_pl("A_C", &standard_pls)],
            ),
            vc(
                "test",
                start,
                vec![a_ref("A"), a_alt("ATC"), a_non_ref()],
                vec![g_pl("A_ATC", &standard_pls)],
            ),
        ],
        "merge_snp_three_alleles" | "test03" => vec![
            vc(
                "test",
                start,
                vec![a_ref("A"), a_alt("C"), a_non_ref()],
                vec![g_pl("A_C", &standard_pls)],
            ),
            vc(
                "test",
                start,
                vec![a_ref("A"), a_alt("C"), a_alt("G"), a_non_ref()],
                vec![g_pl("A_C_G", &[40, 20, 30, 20, 10, 30, 71, 72, 73, 74])],
            ),
        ],
        "merge_snps_ref" | "test04" => vec![
            vc(
                "test",
                start,
                vec![a_ref("A"), a_alt("C"), a_non_ref()],
                vec![g_pl("A_C", &standard_pls)],
            ),
            vc(
                "test",
                start,
                vec![a_ref("A"), a_non_ref()],
                vec![g_pl("A", &[0, 100, 1000])],
            ),
        ],
        "merge_spanning_del" | "test06" => vec![
            vc(
                "test",
                start,
                vec![a_ref("A"), a_alt("C"), a_non_ref()],
                vec![g_pl("A_C", &standard_pls)],
            ),
            vc(
                "test2",
                start - 1,
                vec![a_ref("AA"), a_alt("A"), a_non_ref()],
                vec![g_pl("AA_A", &standard_pls)],
            ),
        ],
        "merge_all_combined" | "test07" => vec![
            vc(
                "test",
                start,
                vec![a_ref("A"), a_alt("C"), a_non_ref()],
                vec![g_pl("A_C", &standard_pls)],
            ),
            vc(
                "test",
                start,
                vec![a_ref("A"), a_alt("G"), a_non_ref()],
                vec![g_pl("A_G", &standard_pls)],
            ),
            vc(
                "test",
                start,
                vec![a_ref("A"), a_alt("ATC"), a_non_ref()],
                vec![g_pl("A_ATC", &standard_pls)],
            ),
            vc(
                "test",
                start,
                vec![a_ref("A"), a_alt("C"), a_alt("G"), a_non_ref()],
                vec![g_pl("A_C_G", &[40, 20, 30, 20, 10, 30, 71, 72, 73, 74])],
            ),
            vc(
                "test",
                start,
                vec![a_ref("A"), a_non_ref()],
                vec![g_pl("A", &[0, 100, 1000])],
            ),
            vc(
                "test",
                start - 1,
                vec![a_ref("AA"), a_non_ref()],
                vec![g_pl("AA", &[0, 80, 800])],
            ),
            vc(
                "test2",
                start - 1,
                vec![a_ref("AA"), a_alt("A"), a_non_ref()],
                vec![g_pl("AA_A", &standard_pls)],
            ),
        ],
        "merge_spanning_ref_only" | "test08" => vec![vc(
            "test",
            start - 1,
            vec![a_ref("AA"), a_non_ref()],
            vec![g_pl("AA", &[0, 80, 800])],
        )],
        "merge_ad_pl_mix" | "test12" => vec![
            vc(
                "test",
                start,
                vec![a_ref("A"), a_alt("ATC"), a_non_ref()],
                vec![g_pl_ad("A_ATC", &[30, 20, 10, 71, 72, 73], &[20, 10])],
            ),
            vc(
                "test",
                start,
                vec![a_ref("A"), a_alt("C"), a_alt("G"), a_non_ref()],
                vec![g_pl_ad(
                    "A_C_G",
                    &[40, 20, 30, 20, 10, 30, 71, 72, 73, 74],
                    &[30, 0, 8],
                )],
            ),
        ],
        "merge_ad_only_overlap" | "test13" => vec![
            vc(
                "test",
                start,
                vec![a_ref("A"), a_alt("C"), a_alt("G"), a_non_ref()],
                vec![g_ad("A_C_G", &[60, 9, 20])],
            ),
            vc(
                "test",
                start,
                vec![a_ref("A"), a_alt("C"), a_alt("G"), a_non_ref()],
                vec![g_ad("A_C_G", &[60, 9, 20])],
            ),
        ],
        "merge_ad_only_distinct" | "test14" => vec![
            vc(
                "test",
                start,
                vec![a_ref("A"), a_alt("C"), a_alt("G"), a_non_ref()],
                vec![g_ad("A_C_G", &[60, 9, 20])],
            ),
            vc(
                "test",
                start,
                vec![a_ref("A"), a_alt("ATC"), a_alt("AA"), a_non_ref()],
                vec![g_ad("A_ATC_AA", &[30, 8, 40])],
            ),
        ],
        _ => {
            return Err(gatk_common::GatkError::argument(format!(
                "unknown merge fixture: {case_id}"
            )))
        }
    };
    Ok((contig, start, vcs))
}

fn result_to_dump(case_id: &str, merged: RefConfidenceMergeResult) -> RefConfidenceMergeDump {
    RefConfidenceMergeDump {
        case_id: case_id.to_string(),
        contig: merged.contig,
        pos: merged.pos,
        alleles: merged.alleles,
        genotypes: merged
            .genotypes
            .into_iter()
            .map(|g| DumpGenotype {
                name: g.name,
                pl: g.pl,
                ad: g.ad,
                gq: None,
            })
            .collect(),
        has_non_ref: merged.has_non_ref,
        merged_null: false,
    }
}

fn null_merge_dump(case_id: &str) -> RefConfidenceMergeDump {
    RefConfidenceMergeDump {
        case_id: case_id.to_string(),
        contig: "20".to_string(),
        pos: 10,
        alleles: Vec::new(),
        genotypes: Vec::new(),
        has_non_ref: false,
        merged_null: true,
    }
}

pub fn merge_reference_confidence_fixture(case_id: &str) -> GatkResult<RefConfidenceMergeDump> {
    if case_id == "merge_spanning_ref_only" || case_id == "test08" {
        let (contig, start, vcs) = unit_test_vcs(case_id)?;
        if merge_reference_confidence(&contig, start, &vcs, None, true, false)?.is_none() {
            return Ok(null_merge_dump(case_id));
        }
    }
    let (contig, start, vcs) = unit_test_vcs(case_id)?;
    let merged = merge_reference_confidence(&contig, start, &vcs, None, true, false)?
        .ok_or_else(|| gatk_common::GatkError::argument("merge returned null"))?;
    Ok(result_to_dump(case_id, merged))
}

pub fn dump_ref_confidence_merge_tsv(
    merge: &RefConfidenceMergeDump,
    out: &mut impl Write,
) -> GatkResult<()> {
    writeln!(out, "case_id\t{}", merge.case_id)?;
    writeln!(out, "contig\t{}", merge.contig)?;
    writeln!(out, "pos\t{}", merge.pos)?;
    writeln!(out, "merged_null\t{}", merge.merged_null)?;
    if merge.merged_null {
        return Ok(());
    }
    writeln!(out, "allele_count\t{}", merge.alleles.len())?;
    for (i, a) in merge.alleles.iter().enumerate() {
        writeln!(out, "allele_{i}\t{a}")?;
    }
    writeln!(out, "sample_count\t{}", merge.genotypes.len())?;
    for (gi, g) in merge.genotypes.iter().enumerate() {
        writeln!(out, "genotype_{gi}_name\t{}", g.name)?;
        if let Some(pl) = &g.pl {
            writeln!(
                out,
                "genotype_{gi}_pl\t{}",
                pl.iter()
                    .map(|v| v.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            )?;
        }
        if let Some(ad) = &g.ad {
            writeln!(
                out,
                "genotype_{gi}_ad\t{}",
                ad.iter()
                    .map(|v| v.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            )?;
        }
        if let Some(gq) = g.gq {
            writeln!(out, "genotype_{gi}_gq\t{gq}")?;
        }
    }
    writeln!(out, "has_non_ref\t{}", merge.has_non_ref)?;
    Ok(())
}

pub fn dump_ref_confidence_merge_case(case_id: &str, out: &mut impl Write) -> GatkResult<()> {
    dump_ref_confidence_merge_tsv(&merge_reference_confidence_fixture(case_id)?, out)
}
