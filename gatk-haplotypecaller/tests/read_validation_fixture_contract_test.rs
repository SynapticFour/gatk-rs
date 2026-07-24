use gatk_common::GatkError;
use gatk_haplotypecaller::validate_mapped_read_sanity;
use serde::Deserialize;
use std::path::PathBuf;

#[derive(Debug, Deserialize)]
struct MalformedReadCase {
    label: String,
    read_len: usize,
    qual_len: usize,
    reference_start: i64,
    reference_end: i64,
    #[serde(default)]
    expect_ok: bool,
    #[serde(default)]
    expect_error_contains: String,
}

fn cases_fixture_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../parity/fixtures/malformed_read_cases.json")
}

#[test]
fn malformed_read_fixture_corpus_matches_error_contracts() {
    let text = std::fs::read_to_string(cases_fixture_path()).unwrap();
    let cases: Vec<MalformedReadCase> = serde_json::from_str(&text).unwrap();
    assert!(!cases.is_empty());

    for case in cases {
        let res = validate_mapped_read_sanity(
            case.read_len,
            case.qual_len,
            case.reference_start,
            case.reference_end,
        );

        if case.expect_ok {
            assert!(
                res.is_ok(),
                "fixture '{}' expected Ok, got {res:?}",
                case.label
            );
            continue;
        }

        let err = res.unwrap_err();
        match err {
            GatkError::Read { message, .. } => assert!(
                message.contains(&case.expect_error_contains),
                "fixture '{}' message mismatch: '{}'",
                case.label,
                message
            ),
            other => panic!(
                "fixture '{}' expected Read error, got {other:?}",
                case.label
            ),
        }
    }
}
