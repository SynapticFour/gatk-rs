//! VCF benchmarks for GATK-RS
//! Performance benchmarks to ensure VCF parsers meet GATK standards.

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use gatk_core::io::*;
use gatk_core::tests::*;
use std::time::Duration;

fn benchmark_vcf_parsing(c: &mut Criterion) {
    let test_data = TestData::new();
    let sizes = vec![100, 1000, 10000];

    let mut group = c.benchmark_group("vcf_parsing");
    group.measurement_time(Duration::from_secs(10));

    for size in sizes {
        let mut content = String::new();
        content.push_str("##fileformat=VCFv4.2\n");
        content.push_str("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
        content.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n");

        for i in 0..size {
            content.push_str(&format!(
                "chr1\t{}\t.\tA\tT\t60.00\tPASS\t.\tGT\t0/1\n",
                i * 100
            ));
        }

        let vcf_path = test_data.create_file(format!("test_{}.vcf", size), &content);

        group.bench_with_input(BenchmarkId::new("parsing", size), &vcf_path, |b, path| {
            b.iter(|| {
                let mut reader = VcfReader::from_file(black_box(path)).unwrap();
                let records = reader.read_all_records().unwrap();
                black_box(records)
            })
        });
    }

    group.finish();
}

fn benchmark_vcf_iterator(c: &mut Criterion) {
    let test_data = TestData::new();
    let mut content = String::new();
    content.push_str("##fileformat=VCFv4.2\n");
    content.push_str("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    content.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n");

    for i in 0..1000 {
        content.push_str(&format!(
            "chr1\t{}\t.\tA\tT\t60.00\tPASS\t.\tGT\t0/1\n",
            i * 100
        ));
    }

    let vcf_path = test_data.create_file("iterator.vcf", &content);

    let mut group = c.benchmark_group("vcf_iterator");
    group.bench_function("iterator", |b| {
        b.iter(|| {
            let mut reader = VcfReader::from_file(black_box(&vcf_path)).unwrap();
            let count: usize = reader.iter().count();
            black_box(count)
        })
    });

    group.finish();
}

fn benchmark_vcf_writing(c: &mut Criterion) {
    let test_data = TestData::new();

    let records: Vec<VcfRecord> = (0..1000)
        .map(|i| VcfRecord {
            chromosome: "chr1".to_string(),
            position: (i * 100) as u64,
            id: ".".to_string(),
            reference: "A".to_string(),
            alternate: vec!["T".to_string()],
            quality: Some(60.0),
            filter: vec!["PASS".to_string()],
            info: vec![
                InfoValue::Float("AF".to_string(), vec![0.5]),
                InfoValue::Integer("DP".to_string(), vec![100]),
            ],
            format: vec!["GT".to_string(), "GQ".to_string(), "DP".to_string()],
            samples: vec![SampleData {
                gt: Some(Genotype {
                    alleles: vec![0, 1],
                    phased: false,
                }),
                gq: Some(60.0),
                dp: Some(50),
                ad: Some(vec![25, 25]),
                pl: Some(vec![100, 0, 200]),
                other: Vec::new(),
            }],
        })
        .collect();

    let header = VcfHeader {
        file_format: "VCFv4.2".to_string(),
        source: Some("GATK".to_string()),
        reference: Some("GRCh38".to_string()),
        contigs: vec![Contig {
            id: "chr1".to_string(),
            length: Some(248956422),
            md5: Some("test_md5".to_string()),
            assembly: None,
            species: None,
            uri: None,
        }],
        samples: vec!["SAMPLE".to_string()],
        info_fields: vec![InfoField {
            id: "AF".to_string(),
            number: "A".to_string(),
            type_field: "Float".to_string(),
            description: "Allele Frequency".to_string(),
            source: None,
            version: None,
        }],
        format_fields: vec![FormatField {
            id: "GT".to_string(),
            number: "1".to_string(),
            type_field: "String".to_string(),
            description: "Genotype".to_string(),
        }],
        filter_fields: Vec::new(),
        other_headers: Vec::new(),
    };

    let mut group = c.benchmark_group("vcf_writing");
    group.bench_function("writing", |b| {
        let vcf_path = test_data.path().join("bench.vcf");
        b.iter(|| {
            let mut writer =
                VcfWriter::new(black_box(&vcf_path), black_box(header.clone())).unwrap();
            writer.write_header().unwrap();
            writer.write_records(black_box(&records)).unwrap();
        })
    });

    group.finish();
}

fn benchmark_vcf_info_parsing(c: &mut Criterion) {
    let info_strings = [
        "AF=0.5",
        "DP=100",
        "AC=50",
        "AF=0.25,0.25",
        "AC=25,25",
        "FLAG",
        "STR=test_string",
        "AF=0.5;DP=100;AC=50;FLAG;STR=test",
    ];

    let mut group = c.benchmark_group("vcf_info_parsing");

    for (i, info_str) in info_strings.iter().enumerate() {
        group.bench_with_input(BenchmarkId::new("parsing", i), info_str, |b, info| {
            b.iter(|| {
                let info = BufferedVcfReader::parse_info(black_box(info)).unwrap();
                black_box(info)
            })
        });
    }

    group.finish();
}

fn benchmark_vcf_genotype_parsing(c: &mut Criterion) {
    let genotype_strings = ["0/1", "1/1", "./.", "0|1", "1|1", "0/2", "1/2", "2/2"];

    let format = vec!["GT".to_string(), "GQ".to_string(), "DP".to_string()];

    let mut group = c.benchmark_group("vcf_genotype_parsing");

    for (i, gt_str) in genotype_strings.iter().enumerate() {
        group.bench_with_input(BenchmarkId::new("parsing", i), gt_str, |b, gt| {
            b.iter(|| {
                let sample = BufferedVcfReader::parse_sample(black_box(gt), black_box(&format));
                black_box(sample)
            })
        });
    }

    group.finish();
}

fn benchmark_vcf_operations(c: &mut Criterion) {
    let records: Vec<VcfRecord> = (0..1000)
        .map(|i| VcfRecord {
            chromosome: "chr1".to_string(),
            position: (i * 100) as u64,
            id: ".".to_string(),
            reference: "A".to_string(),
            alternate: vec!["T".to_string()],
            quality: Some(60.0),
            filter: vec!["PASS".to_string()],
            info: vec![
                InfoValue::Float("AF".to_string(), vec![0.5]),
                InfoValue::Integer("DP".to_string(), vec![100]),
            ],
            format: vec!["GT".to_string(), "GQ".to_string(), "DP".to_string()],
            samples: vec![SampleData {
                gt: Some(Genotype {
                    alleles: vec![0, 1],
                    phased: false,
                }),
                gq: Some(60.0),
                dp: Some(50),
                ad: Some(vec![25, 25]),
                pl: Some(vec![100, 0, 200]),
                other: Vec::new(),
            }],
        })
        .collect();

    let mut group = c.benchmark_group("vcf_operations");

    group.bench_function("variant_type_detection", |b| {
        b.iter(|| {
            let count: usize = black_box(&records)
                .iter()
                .filter(|r| r.is_snp() || r.is_insertion() || r.is_deletion())
                .count();
            black_box(count)
        })
    });

    group.bench_function("allele_frequency_access", |b| {
        b.iter(|| {
            let sum: f64 = black_box(&records).iter().filter_map(|r| r.get_af()).sum();
            black_box(sum)
        })
    });

    group.bench_function("sample_data_access", |b| {
        b.iter(|| {
            let sum: u32 = black_box(&records).iter().filter_map(|r| r.get_dp(0)).sum();
            black_box(sum)
        })
    });

    group.finish();
}

criterion_group!(
    benches,
    benchmark_vcf_parsing,
    benchmark_vcf_iterator,
    benchmark_vcf_writing,
    benchmark_vcf_info_parsing,
    benchmark_vcf_genotype_parsing,
    benchmark_vcf_operations
);

criterion_main!(benches);
