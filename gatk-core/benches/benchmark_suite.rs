//! Comprehensive benchmark suite for GATK-RS

use criterion::{criterion_group, criterion_main, Criterion};
use gatk_core::benchmarking::*;
use std::time::Duration;

use std::hint::black_box;
/// Benchmark the benchmarking framework itself
fn benchmark_framework_overhead(c: &mut Criterion) {
    let mut group = c.benchmark_group("framework_overhead");

    // Benchmark result creation
    group.bench_function("create_result", |b| {
        b.iter(|| {
            let result = BenchmarkResult::success(
                "test".to_string(),
                "test command".to_string(),
                Duration::from_millis(100),
                1024 * 1024 * 100, // 100MB
                50.0,
                1000,
                1024 * 1024, // 1MB
            );
            black_box(result);
        });
    });

    // Benchmark comparison creation
    group.bench_function("create_comparison", |b| {
        let gatk_result = BenchmarkResult::success(
            "GATK".to_string(),
            "test command".to_string(),
            Duration::from_millis(200),
            1024 * 1024 * 200, // 200MB
            80.0,
            1000,
            1024 * 1024, // 1MB
        );

        let gatk_rs_result = BenchmarkResult::success(
            "GATK-RS".to_string(),
            "test command".to_string(),
            Duration::from_millis(100),
            1024 * 1024 * 100, // 100MB
            50.0,
            1000,
            1024 * 1024, // 1MB
        );

        b.iter(|| {
            let comparison = ComparisonResult::new(&gatk_result, &gatk_rs_result);
            black_box(comparison);
        });
    });

    group.finish();
}

/// Benchmark metrics collection
fn benchmark_metrics_collection(c: &mut Criterion) {
    let mut group = c.benchmark_group("metrics_collection");

    group.bench_function("metrics_collector", |b| {
        b.iter(|| {
            let mut collector = MetricsCollector::new();
            collector.start_timing();
            std::thread::sleep(Duration::from_millis(1));
            let _duration = collector.stop_timing("test_metric");
            collector.add_metric("memory", MetricValue::Memory(1024 * 1024));
            collector.add_metric("count", MetricValue::Count(1000));
            black_box(collector);
        });
    });

    group.finish();
}

/// Benchmark suite operations
fn benchmark_suite_operations(c: &mut Criterion) {
    let mut group = c.benchmark_group("suite_operations");

    // Benchmark suite creation
    group.bench_function("create_suite", |b| {
        let config = BenchmarkConfig::default();
        b.iter(|| {
            let suite = BenchmarkSuite::new(config.clone());
            black_box(suite);
        });
    });

    // Benchmark suite with results
    group.bench_function("suite_with_results", |b| {
        let config = BenchmarkConfig::default();
        b.iter(|| {
            let mut suite = BenchmarkSuite::new(config.clone());

            // Add some results
            for i in 0..10 {
                let gatk_result = BenchmarkResult::success(
                    "GATK".to_string(),
                    format!("command {}", i),
                    Duration::from_millis(200 + i * 10),
                    1024 * 1024 * (200 + i * 10),
                    80.0 + i as f64 * 2.0,
                    1000 + i * 10,
                    1024 * 1024,
                );
                suite.add_gatk_result(gatk_result);

                let gatk_rs_result = BenchmarkResult::success(
                    "GATK-RS".to_string(),
                    format!("command {}", i),
                    Duration::from_millis(100 + i * 5),
                    1024 * 1024 * (100 + i * 5),
                    50.0 + i as f64 * 1.0,
                    1000 + i * 10,
                    1024 * 1024,
                );
                suite.add_gatk_rs_result(gatk_rs_result);
            }

            suite.generate_comparisons();
            black_box(suite);
        });
    });

    // Benchmark suite serialization
    group.bench_function("suite_serialization", |b| {
        let config = BenchmarkConfig::default();
        let mut suite = BenchmarkSuite::new(config);

        // Add some results
        for i in 0..100 {
            let gatk_result = BenchmarkResult::success(
                "GATK".to_string(),
                format!("command {}", i),
                Duration::from_millis(200 + i * 10),
                1024 * 1024 * (200 + i * 10),
                80.0 + i as f64 * 2.0,
                1000 + i * 10,
                1024 * 1024,
            );
            suite.add_gatk_result(gatk_result);

            let gatk_rs_result = BenchmarkResult::success(
                "GATK-RS".to_string(),
                format!("command {}", i),
                Duration::from_millis(100 + i * 5),
                1024 * 1024 * (100 + i * 5),
                50.0 + i as f64 * 1.0,
                1000 + i * 10,
                1024 * 1024,
            );
            suite.add_gatk_rs_result(gatk_rs_result);
        }

        suite.generate_comparisons();

        b.iter(|| {
            let json = serde_json::to_string(&suite).unwrap();
            black_box(json);
        });
    });

    group.finish();
}

/// Benchmark performance analysis
fn benchmark_performance_analysis(c: &mut Criterion) {
    let mut group = c.benchmark_group("performance_analysis");

    group.bench_function("analyze_small_dataset", |b| {
        let mut analyzer = PerformanceAnalyzer::new();

        // Add small dataset
        for i in 0..10 {
            analyzer.add_result(BenchmarkResult::success(
                "GATK-RS".to_string(),
                format!("command {}", i),
                Duration::from_millis(100 + i * 5),
                1024 * 1024 * (100 + i * 5),
                50.0 + i as f64 * 1.0,
                1000 + i * 10,
                1024 * 1024,
            ));
        }

        b.iter(|| {
            let stats = analyzer.analyze();
            black_box(stats);
        });
    });

    group.bench_function("analyze_large_dataset", |b| {
        let mut analyzer = PerformanceAnalyzer::new();

        // Add large dataset
        for i in 0..1000 {
            analyzer.add_result(BenchmarkResult::success(
                "GATK-RS".to_string(),
                format!("command {}", i),
                Duration::from_millis(100 + i * 5),
                1024 * 1024 * (100 + i * 5),
                50.0 + i as f64 * 1.0,
                1000 + i * 10,
                1024 * 1024,
            ));
        }

        b.iter(|| {
            let stats = analyzer.analyze();
            black_box(stats);
        });
    });

    group.finish();
}

/// Benchmark regression detection
fn benchmark_regression_detection(c: &mut Criterion) {
    let mut group = c.benchmark_group("regression_detection");

    // Create baseline stats
    let mut baseline_analyzer = PerformanceAnalyzer::new();
    for i in 0..100 {
        baseline_analyzer.add_result(BenchmarkResult::success(
            "GATK-RS".to_string(),
            format!("command {}", i),
            Duration::from_millis(100),
            1024 * 1024 * 100,
            50.0,
            1000,
            1024 * 1024,
        ));
    }
    let baseline_stats = baseline_analyzer.analyze();

    group.bench_function("detect_no_regression", |b| {
        let mut current_analyzer = PerformanceAnalyzer::new();
        for i in 0..100 {
            current_analyzer.add_result(BenchmarkResult::success(
                "GATK-RS".to_string(),
                format!("command {}", i),
                Duration::from_millis(95), // Slightly better
                1024 * 1024 * 95,          // Slightly better
                48.0,                      // Slightly better
                1000,
                1024 * 1024,
            ));
        }
        let current_stats = current_analyzer.analyze();

        let detector = RegressionDetector::new(baseline_stats.clone(), 10.0);
        b.iter(|| {
            let report = detector.detect_regression(&current_stats);
            black_box(report);
        });
    });

    group.bench_function("detect_regression", |b| {
        let mut current_analyzer = PerformanceAnalyzer::new();
        for i in 0..100 {
            current_analyzer.add_result(BenchmarkResult::success(
                "GATK-RS".to_string(),
                format!("command {}", i),
                Duration::from_millis(150), // Worse
                1024 * 1024 * 150,          // Worse
                80.0,                       // Worse
                1000,
                1024 * 1024,
            ));
        }
        let current_stats = current_analyzer.analyze();

        let detector = RegressionDetector::new(baseline_stats.clone(), 10.0);
        b.iter(|| {
            let report = detector.detect_regression(&current_stats);
            black_box(report);
        });
    });

    group.finish();
}

/// Benchmark report generation
fn benchmark_report_generation(c: &mut Criterion) {
    let mut group = c.benchmark_group("report_generation");

    // Create test suite
    let config = BenchmarkConfig::default();
    let mut suite = BenchmarkSuite::new(config);

    for i in 0..50 {
        let gatk_result = BenchmarkResult::success(
            "GATK".to_string(),
            format!("command {}", i),
            Duration::from_millis(200 + i * 10),
            1024 * 1024 * (200 + i * 10),
            80.0 + i as f64 * 2.0,
            1000 + i * 10,
            1024 * 1024,
        );
        suite.add_gatk_result(gatk_result);

        let gatk_rs_result = BenchmarkResult::success(
            "GATK-RS".to_string(),
            format!("command {}", i),
            Duration::from_millis(100 + i * 5),
            1024 * 1024 * (100 + i * 5),
            50.0 + i as f64 * 1.0,
            1000 + i * 10,
            1024 * 1024,
        );
        suite.add_gatk_rs_result(gatk_rs_result);
    }

    suite.generate_comparisons();

    group.bench_function("generate_json_report", |b| {
        let reporter = BenchmarkReporter::new("test_output".to_string(), ReportFormat::Json);
        b.iter(|| {
            let result = reporter.generate_json_report(&suite);
            let _ = black_box(result);
        });
    });

    group.bench_function("generate_markdown_report", |b| {
        let reporter = BenchmarkReporter::new("test_output".to_string(), ReportFormat::Markdown);
        b.iter(|| {
            let result = reporter.generate_markdown_report(&suite);
            let _ = black_box(result);
        });
    });

    group.bench_function("generate_csv_report", |b| {
        let reporter = BenchmarkReporter::new("test_output".to_string(), ReportFormat::Csv);
        b.iter(|| {
            let result = reporter.generate_csv_report(&suite);
            let _ = black_box(result);
        });
    });

    group.finish();
}

/// Benchmark dataset operations
fn benchmark_dataset_operations(c: &mut Criterion) {
    let mut group = c.benchmark_group("dataset_operations");

    group.bench_function("create_dataset_info", |b| {
        b.iter(|| {
            let info = DatasetInfo {
                name: "test".to_string(),
                size: 1024 * 1024 * 1024,
                read_count: 1_000_000,
                reference_size: 3 * 1024 * 1024 * 1024,
                description: "Test dataset".to_string(),
                url: None,
                checksum: None,
                created_at: chrono::Utc::now(),
                category: DatasetCategory::Synthetic,
                complexity: DatasetComplexity::Medium,
            };
            black_box(info);
        });
    });

    group.bench_function("validate_dataset_info", |b| {
        let info = DatasetInfo {
            name: "test".to_string(),
            size: 1024 * 1024 * 1024,
            read_count: 1_000_000,
            reference_size: 3 * 1024 * 1024 * 1024,
            description: "Test dataset".to_string(),
            url: None,
            checksum: None,
            created_at: chrono::Utc::now(),
            category: DatasetCategory::Synthetic,
            complexity: DatasetComplexity::Medium,
        };

        b.iter(|| {
            let validation = info.size > 0 && info.read_count > 0 && info.reference_size > 0;
            black_box(validation);
        });
    });

    group.finish();
}

/// Benchmark comparison operations
fn benchmark_comparison_operations(c: &mut Criterion) {
    let mut group = c.benchmark_group("comparison_operations");

    // Create test results
    let gatk_result = BenchmarkResult::success(
        "GATK".to_string(),
        "test command".to_string(),
        Duration::from_millis(200),
        1024 * 1024 * 200,
        80.0,
        1000,
        1024 * 1024,
    );

    let gatk_rs_result = BenchmarkResult::success(
        "GATK-RS".to_string(),
        "test command".to_string(),
        Duration::from_millis(100),
        1024 * 1024 * 100,
        50.0,
        1000,
        1024 * 1024,
    );

    group.bench_function("create_comparison", |b| {
        b.iter(|| {
            let comparison = ComparisonResult::new(&gatk_result, &gatk_rs_result);
            black_box(comparison);
        });
    });

    group.bench_function("check_targets", |b| {
        let comparison = ComparisonResult::new(&gatk_result, &gatk_rs_result);
        b.iter(|| {
            let meets_targets = comparison.meets_targets();
            black_box(meets_targets);
        });
    });

    group.finish();
}

criterion_group!(
    benches,
    benchmark_framework_overhead,
    benchmark_metrics_collection,
    benchmark_suite_operations,
    benchmark_performance_analysis,
    benchmark_regression_detection,
    benchmark_report_generation,
    benchmark_dataset_operations,
    benchmark_comparison_operations
);

criterion_main!(benches);
