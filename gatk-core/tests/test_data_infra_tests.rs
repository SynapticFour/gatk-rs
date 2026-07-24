use gatk_core::tests::fixture_path;

#[test]
fn shared_fixture_files_exist_and_are_readable() {
    let reference = fixture_path("reference.fa");
    let vcf = fixture_path("sample.vcf");
    let fastq = fixture_path("sample.fastq");
    let bam = fixture_path("sample.bam");

    for path in [&reference, &vcf, &fastq, &bam] {
        assert!(path.exists(), "fixture should exist: {}", path.display());
        let bytes = std::fs::read(path).expect("fixture should be readable");
        assert!(
            !bytes.is_empty(),
            "fixture should not be empty: {}",
            path.display()
        );
    }
}
