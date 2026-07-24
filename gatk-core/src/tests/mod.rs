//! Lightweight test/benchmark support helpers.
//! Note: legacy in-crate test suites are intentionally not re-exported here.

use std::path::{Path, PathBuf};
use tempfile::TempDir;

type VcfVariantRecord<'a> = (
    &'a str,
    u64,
    &'a str,
    &'a str,
    &'a str,
    f64,
    &'a str,
    &'a str,
);

/// Temporary fixture manager for integration tests and local benchmarks.
/// # Invariants
/// Owns a [`tempfile::TempDir`] destroyed on drop unless persisted externally.
/// # Ownership
/// Owns temp directory; paths borrowed via helper methods.
/// # Mutation
/// File creation methods write into the temp dir; not thread-safe.
/// # Biological assumptions
/// Generates toy FASTA/FASTQ/VCF/BAM fixtures for tests (not production data).
/// # Java equivalence
/// None / Rust-native test harness.
pub struct TestData {
    temp_dir: TempDir,
}

impl TestData {
    /// Create a new test data directory
    pub fn new() -> Self {
        Self {
            temp_dir: TempDir::new().expect("Failed to create temp directory"),
        }
    }

    /// Get the path to the temporary directory
    pub fn path(&self) -> &Path {
        self.temp_dir.path()
    }

    /// Create a test file with the given content
    pub fn create_file<P: AsRef<Path>>(&self, name: P, content: &str) -> PathBuf {
        let file_path = self.temp_dir.path().join(name.as_ref());
        std::fs::write(&file_path, content).expect("Failed to write test file");
        file_path
    }

    /// Create a test FASTA file
    pub fn create_fasta<P: AsRef<Path>>(&self, name: P, sequences: &[(&str, &str)]) -> PathBuf {
        let mut content = String::new();
        for (header, sequence) in sequences {
            content.push('>');
            content.push_str(header);
            content.push('\n');
            content.push_str(sequence);
            content.push('\n');
        }
        self.create_file(name, &content)
    }

    /// Create a test VCF file
    pub fn create_vcf<P: AsRef<Path>>(
        &self,
        name: P,
        variants: &[VcfVariantRecord<'_>],
    ) -> PathBuf {
        let mut content = String::new();
        content.push_str("##fileformat=VCFv4.2\n");
        content.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");

        for (chrom, pos, id, ref_allele, alt_allele, qual, filter, info) in variants {
            content.push_str(&format!(
                "{chrom}\t{pos}\t{id}\t{ref_allele}\t{alt_allele}\t{qual}\t{filter}\t{info}\n"
            ));
        }

        self.create_file(name, &content)
    }

    /// Create a test BAM file (simplified - just creates a placeholder)
    pub fn create_bam<P: AsRef<Path>>(&self, name: P) -> PathBuf {
        // In a real implementation, this would create a proper BAM file
        // For now, we create a placeholder file
        self.create_file(name, "BAM_PLACEHOLDER")
    }
}

impl Default for TestData {
    fn default() -> Self {
        Self::new()
    }
}

/// Get path to checked-in test fixture file.
pub fn fixture_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("src")
        .join("tests")
        .join("test_data")
        .join(name)
}

/// Macro for creating property-based tests
#[macro_export]
macro_rules! prop_test {
    ($name:ident in |$param:ident : $ty:ty| $body:block) => {
        ::proptest::proptest! {
            #[test]
            fn $name($param in ::proptest::arbitrary::any::<$ty>()) $body
        }
    };
    ($name:ident in |$($param:ident : $ty:ty),+| $body:block) => {
        ::proptest::proptest! {
            #[test]
            fn $name($($param in ::proptest::arbitrary::any::<$ty>()),+) $body
        }
    };
}

/// Macro for creating benchmark tests
#[macro_export]
macro_rules! bench_test {
    ($name:ident, $bench:expr) => {
        pub fn $name(c: &mut ::criterion::Criterion) {
            c.bench_function(stringify!($name), |b| {
                ($bench)(b);
            });
        }
    };
}

/// Helper function to create test genomic positions
pub fn test_position(chrom: &str, pos: u64) -> crate::types::GenomicPosition {
    let contig = chrom.trim_start_matches("chr").parse::<u32>().unwrap_or(0);
    crate::types::GenomicPosition {
        contig,
        position: pos,
    }
}

/// Helper function to create test alleles
pub fn test_allele(allele: &str) -> crate::types::Allele {
    let bases: Vec<crate::types::Base> = allele
        .as_bytes()
        .iter()
        .map(|&b| match b {
            b'A' | b'a' => crate::types::Base::A,
            b'T' | b't' => crate::types::Base::T,
            b'G' | b'g' => crate::types::Base::G,
            b'C' | b'c' => crate::types::Base::C,
            _ => crate::types::Base::N,
        })
        .collect();
    crate::types::Allele::new(bases)
}

/// Helper function to create test sequences
pub fn test_sequence(seq: &str) -> Vec<crate::types::Base> {
    seq.as_bytes()
        .iter()
        .map(|&b| match b {
            b'A' | b'a' => crate::types::Base::A,
            b'T' | b't' => crate::types::Base::T,
            b'G' | b'g' => crate::types::Base::G,
            b'C' | b'c' => crate::types::Base::C,
            _ => crate::types::Base::N,
        })
        .collect()
}

/// Helper function to create test variants
pub fn test_variant(
    chrom: &str,
    pos: u64,
    ref_allele: &str,
    alt_allele: &str,
) -> crate::types::VariantContext {
    crate::types::VariantContext::new(
        test_position(chrom, pos),
        test_allele(ref_allele),
        vec![test_allele(alt_allele)],
    )
}

/// Helper function to create test reads
pub fn test_read(name: &str, sequence: &str, qualities: &[u8]) -> crate::types::SequenceRead {
    let bases = sequence
        .bytes()
        .map(|b| match b {
            b'A' | b'a' => crate::types::Base::A,
            b'C' | b'c' => crate::types::Base::C,
            b'G' | b'g' => crate::types::Base::G,
            b'T' | b't' => crate::types::Base::T,
            _ => crate::types::Base::N,
        })
        .collect();
    let read_quality = crate::types::ReadQuality::from_vec(qualities.to_vec());
    crate::types::SequenceRead::new(
        name.to_string(),
        bases,
        read_quality,
        test_position("chr1", 100),
        false, // is_reverse_strand
        false, // is_paired
    )
}

/// Test constants
pub mod constants {
    pub const TEST_CHROMOSOME: &str = "chr1";
    pub const TEST_POSITION: u64 = 1000;
    pub const TEST_SEQUENCE: &str = "ATCGATCGATCG";
    pub const TEST_REF_ALLELE: &str = "A";
    pub const TEST_ALT_ALLELE: &str = "T";
    pub const TEST_QUALITY: u8 = 30;
    pub const TEST_READ_NAME: &str = "read1";
}

/// Integration test helpers module
pub mod integration_test_helpers {
    use super::TestData;

    /// Create a complete test dataset for integration testing
    pub fn create_test_dataset() -> TestData {
        let data = TestData::new();

        // Create reference FASTA
        data.create_fasta(
            "reference.fa",
            &[
                ("chr1", "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"),
                ("chr2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"),
                ("chr3", "TTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAA"),
            ],
        );

        // Create test VCF
        data.create_vcf(
            "variants.vcf",
            &[
                ("chr1", 100, ".", "A", "T", 50.0, "PASS", "."),
                ("chr1", 200, ".", "C", "G", 60.0, "PASS", "."),
                ("chr2", 150, ".", "G", "A", 70.0, "PASS", "."),
            ],
        );

        // Create test BAM placeholder
        data.create_bam("reads.bam");

        data
    }

    /// Run a complete integration test
    pub fn run_integration_test<F>(test_fn: F) -> Result<(), Box<dyn std::error::Error>>
    where
        F: FnOnce(&TestData) -> Result<(), Box<dyn std::error::Error>>,
    {
        let test_data = create_test_dataset();
        test_fn(&test_data)
    }
}

/// Performance testing utilities
pub mod performance {
    use std::time::{Duration, Instant};

    /// Measure execution time of a function
    pub fn measure_time<F, R>(f: F) -> (R, Duration)
    where
        F: FnOnce() -> R,
    {
        let start = Instant::now();
        let result = f();
        let duration = start.elapsed();
        (result, duration)
    }

    /// Benchmark memory usage
    pub fn benchmark_memory<F, R>(f: F) -> (R, usize)
    where
        F: FnOnce() -> R,
    {
        let initial_memory = crate::memory::MemoryMonitor::current_memory_usage();
        let result = f();
        let final_memory = crate::memory::MemoryMonitor::current_memory_usage();
        let memory_used = final_memory.saturating_sub(initial_memory);
        (result, memory_used)
    }

    /// Assert performance requirements
    pub fn assert_performance<F, R>(
        f: F,
        max_time: Duration,
        max_memory: usize,
    ) -> Result<R, String>
    where
        F: FnOnce() -> R + Clone,
    {
        let (result, duration) = measure_time(f.clone());
        let (_, memory_used) = benchmark_memory(f);

        if duration > max_time {
            return Err(format!(
                "Time limit exceeded: {:?} > {:?}",
                duration, max_time
            ));
        }

        if memory_used > max_memory {
            return Err(format!(
                "Memory limit exceeded: {} > {}",
                memory_used, max_memory
            ));
        }

        Ok(result)
    }
}

/// Mock implementations for testing
pub mod mocks {
    /// Mock file reader for testing
    pub struct MockFileReader {
        content: Vec<u8>,
        position: usize,
    }

    impl MockFileReader {
        pub fn new(content: &str) -> Self {
            Self {
                content: content.as_bytes().to_vec(),
                position: 0,
            }
        }
    }

    impl std::io::Read for MockFileReader {
        fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
            if self.position >= self.content.len() {
                return Ok(0);
            }

            let remaining = self.content.len() - self.position;
            let to_read = std::cmp::min(buf.len(), remaining);

            buf[..to_read].copy_from_slice(&self.content[self.position..self.position + to_read]);
            self.position += to_read;

            Ok(to_read)
        }
    }

    /// Mock error injector for testing error handling
    pub struct MockErrorInjector {
        should_fail: bool,
        error_message: String,
    }

    impl MockErrorInjector {
        pub fn new(should_fail: bool, error_message: &str) -> Self {
            Self {
                should_fail,
                error_message: error_message.to_string(),
            }
        }

        pub fn check_error(&self) -> gatk_common::GatkResult<()> {
            if self.should_fail {
                Err(gatk_common::GatkError::generic(&self.error_message))
            } else {
                Ok(())
            }
        }
    }
}
