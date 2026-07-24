use gatk_core::io::{
    fasta::{FastaReader, FastaSequence, FastaWriter},
    fastq::{FastqRead, FastqReader, FastqWriter},
};
use gatk_core::tests::TestData;

#[test]
fn fasta_reader_parses_basic_contract() {
    let test_data = TestData::new();
    let fasta_path =
        test_data.create_file("contract.fa", ">chr1 description\nATCGATCG\n>chr2\nGGTT\n");

    let mut reader = FastaReader::from_file_buffered(&fasta_path).unwrap();
    let records = reader.read_all_sequences().unwrap();

    assert_eq!(records.len(), 2);
    assert_eq!(records[0].name, "chr1");
    assert_eq!(records[0].description.as_deref(), Some("description"));
    assert_eq!(records[0].as_bytes(), b"ATCGATCG");
    assert_eq!(records[1].name, "chr2");
    assert_eq!(records[1].as_bytes(), b"GGTT");
}

#[test]
fn fasta_reader_memory_mapped_header_with_spaces_contract() {
    let test_data = TestData::new();
    let fasta_path = test_data.create_file(
        "contract_mmap.fa",
        ">20 dna:chromosome chromosome:GRCh37:20:1:63025520:1\nATCG\n",
    );
    let mut reader = FastaReader::from_file_memory_mapped(&fasta_path).unwrap();
    let records = reader.read_all_sequences().unwrap();
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].name, "20");
    assert_eq!(
        records[0].description.as_deref(),
        Some("dna:chromosome chromosome:GRCh37:20:1:63025520:1")
    );
    assert_eq!(records[0].as_bytes(), b"ATCG");
}

#[test]
fn fasta_sequence_contract_for_reverse_complement_and_validation() {
    let seq = FastaSequence::new("s1".to_string(), b"ATCGNatcgn".to_vec());
    assert!(seq.is_valid_dna());
    assert_eq!(seq.gc_content(), 0.4);

    let rc = seq.reverse_complement();
    assert_eq!(rc.as_bytes(), b"ncgatNCGAT");
}

#[test]
fn fasta_writer_roundtrip_contract() {
    let test_data = TestData::new();
    let out_path = test_data.path().join("out.fa");
    let records = vec![
        FastaSequence::new("ref1".to_string(), b"ATCG".to_vec()),
        FastaSequence::new("ref2".to_string(), b"GGTT".to_vec()),
    ];

    let mut writer = FastaWriter::new(&out_path).unwrap();
    writer.write_sequences(&records).unwrap();
    writer.finish().unwrap();

    let mut reader = FastaReader::from_file_buffered(&out_path).unwrap();
    let parsed = reader.read_all_sequences().unwrap();
    assert_eq!(parsed.len(), 2);
    assert_eq!(parsed[0].name, "ref1");
    assert_eq!(parsed[0].as_bytes(), b"ATCG");
    assert_eq!(parsed[1].name, "ref2");
    assert_eq!(parsed[1].as_bytes(), b"GGTT");
}

#[test]
fn fastq_read_validity_and_reverse_complement_contract() {
    let read = FastqRead::new(
        "r1 sample".to_string(),
        b"ATCGN".to_vec(),
        b"IIIII".to_vec(), // ASCII Phred+33 values in valid range
    );
    assert_eq!(read.name(), "r1");
    assert_eq!(read.description(), Some("sample"));
    assert!(read.is_valid());

    let rc = read.reverse_complement();
    assert_eq!(rc.sequence, b"NCGAT");
    assert_eq!(rc.quality, b"IIIII");
}

#[test]
fn fastq_reader_writer_roundtrip_contract() {
    let test_data = TestData::new();
    let out_path = test_data.path().join("reads.fastq");
    let reads = vec![
        FastqRead::new("r1".to_string(), b"ATCG".to_vec(), b"IIII".to_vec()),
        FastqRead::new("r2".to_string(), b"GGTT".to_vec(), b"HHHH".to_vec()),
    ];

    let mut writer = FastqWriter::new(&out_path).unwrap();
    writer.write_reads(&reads).unwrap();
    writer.finish().unwrap();

    let mut reader = FastqReader::from_file_buffered(&out_path).unwrap();
    let parsed = reader.read_all_reads().unwrap();
    assert_eq!(parsed.len(), 2);
    assert_eq!(parsed[0].header, "r1");
    assert_eq!(parsed[0].sequence, b"ATCG");
    assert_eq!(parsed[1].header, "r2");
    assert_eq!(parsed[1].sequence, b"GGTT");
}
