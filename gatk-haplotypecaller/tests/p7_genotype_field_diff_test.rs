use gatk_haplotypecaller::emit_genotype_format_fields;
use std::collections::BTreeMap;
use std::fs;
use std::path::PathBuf;

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("workspace root")
        .to_path_buf()
}

#[derive(Debug, Clone)]
struct InputCase {
    genotype_log10_likelihoods: Vec<f64>,
    allele_depths: Vec<i32>,
}

#[derive(Debug, Clone)]
struct ExpectedCase {
    pl: Vec<i32>,
    gq: i32,
    ad: Vec<i32>,
    dp: i32,
}

fn parse_f64_csv(raw: &str) -> Vec<f64> {
    raw.split(',')
        .map(|v| v.parse::<f64>().expect("f64 csv value"))
        .collect()
}

fn parse_i32_csv(raw: &str) -> Vec<i32> {
    raw.split(',')
        .map(|v| v.parse::<i32>().expect("i32 csv value"))
        .collect()
}

fn load_inputs(path: &PathBuf) -> BTreeMap<String, InputCase> {
    fs::read_to_string(path)
        .expect("input fixture")
        .lines()
        .filter(|l| !l.trim().is_empty() && !l.starts_with('#'))
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            (
                c[0].to_string(),
                InputCase {
                    genotype_log10_likelihoods: parse_f64_csv(c[1]),
                    allele_depths: parse_i32_csv(c[2]),
                },
            )
        })
        .collect()
}

fn load_expected(path: &PathBuf) -> BTreeMap<String, ExpectedCase> {
    fs::read_to_string(path)
        .expect("expected fixture")
        .lines()
        .filter(|l| !l.trim().is_empty() && !l.starts_with('#'))
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            (
                c[0].to_string(),
                ExpectedCase {
                    pl: parse_i32_csv(c[1]),
                    gq: c[2].parse::<i32>().expect("gq"),
                    ad: parse_i32_csv(c[3]),
                    dp: c[4].parse::<i32>().expect("dp"),
                },
            )
        })
        .collect()
}

#[test]
fn genotype_field_emission_matches_frozen_java_smoke_fixture() {
    let root = repo_root();
    let inputs = load_inputs(&root.join("parity/fixtures/p7_genotype_fields_smoke.tsv"));
    let expected = load_expected(&root.join("parity/expected/p7_genotype_fields_smoke.java.tsv"));

    assert_eq!(inputs.len(), expected.len());
    assert!(!inputs.is_empty());

    for (case_id, input) in &inputs {
        let exp = expected
            .get(case_id)
            .unwrap_or_else(|| panic!("missing expected case: {case_id}"));
        let got =
            emit_genotype_format_fields(&input.genotype_log10_likelihoods, &input.allele_depths)
                .unwrap_or_else(|e| panic!("case={case_id} emit failed: {e}"));
        assert_eq!(got.pl_as_i32(), exp.pl, "case={case_id} PL mismatch");
        assert_eq!(got.gq.as_i32(), exp.gq, "case={case_id} GQ mismatch");
        assert_eq!(got.ad_as_i32(), exp.ad, "case={case_id} AD mismatch");
        assert_eq!(got.dp.as_i32(), exp.dp, "case={case_id} DP mismatch");
    }
}
