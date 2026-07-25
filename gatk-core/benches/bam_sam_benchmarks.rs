//! BAM/SAM benchmarks for GATK-RS
//! Performance benchmarks to ensure BAM/SAM parsers meet GATK standards.

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use gatk_core::io::*;
use gatk_core::tests::*;
use std::time::Duration;

fn benchmark_sam_parsing(c: &mut Criterion) {
    let test_data = TestData::new();
    let sizes = vec![100, 1000, 10000];

    let mut group = c.benchmark_group("sam_parsing");
    group.measurement_time(Duration::from_secs(10));

    for size in sizes {
        let mut content = String::new();
        content.push_str("@HD\tVN:1.6\tSO:coordinate\n");
        content.push_str("@SQ\tSN:chr1\tLN:1000000\n");

        for i in 0..size {
            content.push_str(&format!(
                "read{}\t0\tchr1\t{}\t60\t100M\t*\t0\t0\tACGTACGTACGTACGT\tIIIIIIIIIIIIIIII\n",
                i,
                i * 100
            ));
        }

        let sam_path = test_data.create_file(format!("test_{}.sam", size), &content);

        group.bench_with_input(BenchmarkId::new("parsing", size), &sam_path, |b, path| {
            b.iter(|| {
                let mut reader = SamReader::from_file(black_box(path)).unwrap();
                let records = reader.read_all_records().unwrap();
                black_box(records)
            })
        });
    }

    group.finish();
}

fn benchmark_sam_iterator(c: &mut Criterion) {
    let test_data = TestData::new();
    let mut content = String::new();
    content.push_str("@HD\tVN:1.6\n");
    content.push_str("@SQ\tSN:chr1\tLN:1000000\n");

    for i in 0..1000 {
        content.push_str(&format!(
            "read{}\t0\tchr1\t{}\t60\t100M\t*\t0\t0\tACGTACGTACGTACGT\tIIIIIIIIIIIIIIII\n",
            i,
            i * 100
        ));
    }

    let sam_path = test_data.create_file("iterator.sam", &content);

    let mut group = c.benchmark_group("sam_iterator");
    group.bench_function("iterator", |b| {
        b.iter(|| {
            let mut reader = SamReader::from_file(black_box(&sam_path)).unwrap();
            let count: usize = reader.iter().count();
            black_box(count)
        })
    });

    group.finish();
}

fn benchmark_sam_writing(c: &mut Criterion) {
    let test_data = TestData::new();

    let records: Vec<SamRecord> = (0..1000)
        .map(|i| SamRecord {
            qname: format!("read{}", i),
            flag: 0,
            rname: "chr1".to_string(),
            pos: (i * 100) as i64,
            mapq: 60,
            cigar: "100M".to_string(),
            rnext: "*".to_string(),
            pnext: 0,
            tlen: 0,
            seq: "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_string(),
            qual: "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII".to_string(),
            optional: vec![
                OptionalField::Int("NM".to_string(), 1),
                OptionalField::String("RG".to_string(), "RG1".to_string()),
            ],
        })
        .collect();

    let header = SamHeader {
        reference_sequences: vec![ReferenceSequence {
            name: "chr1".to_string(),
            length: 1000000,
            md5: Some("test_md5".to_string()),
            assembly: None,
            uri: None,
            species: None,
        }],
        read_groups: vec![ReadGroup {
            id: "RG1".to_string(),
            description: Some("Test RG".to_string()),
            flow_order: None,
            key_sequence: None,
            library: None,
            platform_unit: None,
            platform: Some("ILLUMINA".to_string()),
            sample: Some("sample1".to_string()),
        }],
        programs: Vec::new(),
        comments: Vec::new(),
        sort_order: Some("coordinate".to_string()),
        grouping: None,
    };

    let mut group = c.benchmark_group("sam_writing");
    group.bench_function("writing", |b| {
        let sam_path = test_data.path().join("bench.sam");
        b.iter(|| {
            let mut writer =
                SamWriter::new(black_box(&sam_path), black_box(header.clone())).unwrap();
            writer.write_header().unwrap();
            writer.write_records(black_box(&records)).unwrap();
        })
    });

    group.finish();
}

fn benchmark_cigar_parsing(c: &mut Criterion) {
    let cigar_strings = [
        "100M",
        "10M5I10D2M2S",
        "50M10N50M",
        "20M10S20M10H",
        "5=10X5=10X5=",
    ];

    let mut group = c.benchmark_group("cigar_parsing");

    for (i, cigar) in cigar_strings.iter().enumerate() {
        group.bench_with_input(BenchmarkId::new("parse", i), cigar, |b, cigar_str| {
            b.iter(|| {
                let record = SamRecord {
                    qname: "test".to_string(),
                    flag: 0,
                    rname: "chr1".to_string(),
                    pos: 100,
                    mapq: 60,
                    cigar: cigar_str.to_string(),
                    rnext: "*".to_string(),
                    pnext: 0,
                    tlen: 0,
                    seq: "ACGTACGTACGTACGT".to_string(),
                    qual: "IIIIIIIIIIIIIIII".to_string(),
                    optional: Vec::new(),
                };
                let ops = record.parse_cigar();
                black_box(ops)
            })
        });
    }

    group.finish();
}

fn benchmark_optional_fields(c: &mut Criterion) {
    let test_data = TestData::new();

    let sam_content = r#"@HD	VN:1.6
@SQ	SN:chr1	LN:1000000
read1	0	chr1	100	60	100M	*	0	0	ACGTACGTACGTACGT	IIIIIIIIIIIIIIII	NM:i:1	MD:Z:100	AS:i:50	XS:i:+5	BC:B:10,20,30	RG:Z:RG1	H0:H:1A2B3C"#;

    let sam_path = test_data.create_file("optional.sam", sam_content);

    let mut group = c.benchmark_group("optional_fields");

    group.bench_function("parsing", |b| {
        b.iter(|| {
            let mut reader = SamReader::from_file(black_box(&sam_path)).unwrap();
            let record = reader.read_next_record().unwrap().unwrap();
            let nm_field = record.get_optional_field("NM").is_some();
            let md_field = record.get_optional_field("MD").is_some();
            let rg_field = record.get_optional_field("RG").is_some();
            black_box((nm_field, md_field, rg_field))
        })
    });

    group.finish();
}

fn benchmark_bam_parsing(c: &mut Criterion) {
    // `BamReader` is a pure-Rust uncompressed-BAM parser (`BAM\x01` magic). Tracked
    // `parity/fixtures/sample.bam` is BGZF — samtools can read it, BamReader cannot.
    // Skip unless we have an uncompressed BAM the reader can open.
    let candidates = [
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../parity/fixtures/sample.bam"),
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("src/tests/test_data/sample.bam"),
    ];
    let Some(bam_path) = candidates
        .into_iter()
        .find(|p| p.is_file() && BamReader::from_file(p).is_ok())
    else {
        eprintln!("benchmark_bam_parsing: skip (no uncompressed BAM fixture for BamReader)");
        return;
    };

    let mut group = c.benchmark_group("bam_parsing");
    group.measurement_time(Duration::from_secs(5));
    group.bench_with_input(
        BenchmarkId::new("parsing", "sample.bam"),
        &bam_path,
        |b, path| {
            b.iter(|| {
                let mut reader = BamReader::from_file(black_box(path)).unwrap();
                let records = reader.read_all_records().unwrap();
                black_box(records)
            })
        },
    );
    group.finish();
}

criterion_group!(
    benches,
    benchmark_sam_parsing,
    benchmark_sam_iterator,
    benchmark_sam_writing,
    benchmark_cigar_parsing,
    benchmark_optional_fields,
    benchmark_bam_parsing
);

criterion_main!(benches);
