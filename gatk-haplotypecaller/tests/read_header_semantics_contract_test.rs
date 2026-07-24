use gatk_common::GatkError;
use gatk_haplotypecaller::ReadHeaderSemantics;

#[test]
fn rg_tag_maps_to_sample_contract() {
    let header = "@HD\tVN:1.6\n@RG\tID:rgA\tSM:tumor\n@RG\tID:rgB\tSM:normal\n@PG\tID:pgA\n";
    let semantics = ReadHeaderSemantics::from_sam_header_text(header).unwrap();

    let a = semantics
        .validate_record_links(Some("rgA"), Some("pgA"))
        .unwrap();
    let b = semantics.validate_record_links(Some("rgB"), None).unwrap();

    assert_eq!(a.sample_name.as_deref(), Some("tumor"));
    assert_eq!(b.sample_name.as_deref(), Some("normal"));
}

#[test]
fn missing_rg_in_header_is_validation_error() {
    let header = "@HD\tVN:1.6\n@RG\tID:rgA\tSM:tumor\n";
    let semantics = ReadHeaderSemantics::from_sam_header_text(header).unwrap();
    let err = semantics
        .validate_record_links(Some("rgZ"), None)
        .unwrap_err();

    match err {
        GatkError::Validation { message, .. } => {
            assert!(message.contains("record RG"));
            assert!(message.contains("not found"));
        }
        other => panic!("expected Validation error, got {other:?}"),
    }
}
