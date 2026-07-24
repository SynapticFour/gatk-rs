use super::DatasetInfo;
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};

/// Filesystem-backed registry of benchmark datasets under a root directory.
/// # Invariants
/// Each immediate subdirectory of `root` is treated as one dataset name.
/// `list_datasets` best-effort; missing metadata defaults to unknown/low complexity.
/// # Ownership
/// Owns root `PathBuf`; `&self` methods borrow filesystem.
/// # Mutation
/// `cleanup_dataset` removes directories; other methods read-only on manager state.
/// # Biological assumptions
/// Dataset folders expected to contain reference/read fixtures for HC benchmarks.
/// # Java equivalence
/// None / Rust-native.
#[derive(Debug, Clone)]
pub struct DatasetManager {
    root: PathBuf,
}

impl DatasetManager {
    pub fn new<P: AsRef<Path>>(root: P) -> gatk_common::GatkResult<Self> {
        std::fs::create_dir_all(root.as_ref())
            .map_err(|e| gatk_common::GatkError::io("Failed to create dataset root", e))?;
        Ok(Self {
            root: root.as_ref().to_path_buf(),
        })
    }

    pub fn list_datasets(&self) -> Vec<DatasetInfo> {
        let mut out = Vec::new();
        if let Ok(entries) = std::fs::read_dir(&self.root) {
            for entry in entries.flatten() {
                let p = entry.path();
                if p.is_dir() {
                    let name = p
                        .file_name()
                        .map(|s| s.to_string_lossy().to_string())
                        .unwrap_or_default();
                    out.push(DatasetInfo {
                        name: name.clone(),
                        size: dir_size(&p),
                        read_count: 0,
                        reference_size: 0,
                        description: format!("Dataset {name}"),
                        url: None,
                        checksum: None,
                        created_at: chrono::Utc::now(),
                        category: super::DatasetCategory::Unknown,
                        complexity: super::DatasetComplexity::Low,
                    });
                }
            }
        }
        out
    }

    pub fn validate_dataset(&self, name: &str) -> gatk_common::GatkResult<bool> {
        Ok(self.root.join(name).exists())
    }

    pub fn cleanup_dataset(&self, name: &str) -> gatk_common::GatkResult<()> {
        let path = self.root.join(name);
        if path.exists() {
            std::fs::remove_dir_all(&path)
                .map_err(|e| gatk_common::GatkError::io("Failed to remove dataset", e))?;
        }
        Ok(())
    }

    pub fn get_dataset_stats(&self, name: &str) -> gatk_common::GatkResult<DatasetStats> {
        let path = self.root.join(name);
        let total_size = dir_size(&path);
        let modified = std::fs::metadata(&path)
            .ok()
            .and_then(|m| m.modified().ok())
            .map(chrono::DateTime::<chrono::Utc>::from);
        Ok(DatasetStats {
            name: name.to_string(),
            reference_size: total_size / 2,
            reads_size: total_size / 2,
            total_size,
            is_valid: path.exists(),
            last_modified: modified,
        })
    }
}

/// Parameters for synthesizing a small custom benchmark dataset on disk.
/// # Invariants
/// Generated FASTQ capped at 1000 reads regardless of `read_count`.
/// Rates in `[0.0, 1.0]` expected but not enforced here.
/// # Ownership
/// Owns dataset name string; clone to reuse configs.
/// # Mutation
/// Immutable input to generator; public fields for serde/tuning.
/// # Biological assumptions
/// Toy DNA (all-A reference/reads) for smoke benchmarks, not realistic variation.
/// # Java equivalence
/// None / Rust-native.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CustomDatasetConfig {
    pub name: String,
    pub chromosome_count: usize,
    pub chromosome_length: u64,
    pub read_count: usize,
    pub read_length: usize,
    pub include_repetitive_regions: bool,
    pub variant_rate: f64,
    pub error_rate: f64,
}

/// Writes minimal FASTA/FASTQ fixtures from [`CustomDatasetConfig`].
/// # Invariants
/// Output paths: `{root}/{name}/reference.fa` and `reads.fq`.
/// # Ownership
/// Owns generator root path.
/// # Mutation
/// `generate_custom_dataset` writes filesystem state.
/// # Biological assumptions
/// Placeholder sequences for pipeline wiring tests only.
/// # Java equivalence
/// None / Rust-native.
pub struct DatasetGenerator {
    root: PathBuf,
}

impl DatasetGenerator {
    pub fn new<P: AsRef<Path>>(root: P) -> Self {
        Self {
            root: root.as_ref().to_path_buf(),
        }
    }

    pub fn generate_custom_dataset(
        &self,
        config: &CustomDatasetConfig,
    ) -> gatk_common::GatkResult<()> {
        let dir = self.root.join(&config.name);
        std::fs::create_dir_all(&dir)
            .map_err(|e| gatk_common::GatkError::io("Failed to create dataset directory", e))?;

        let reference = dir.join("reference.fa");
        let reads = dir.join("reads.fq");
        std::fs::write(&reference, format!(">chr1\n{}\n", "A".repeat(10_000)))
            .map_err(|e| gatk_common::GatkError::io("Failed to write reference", e))?;

        let mut fastq = String::new();
        let n = config.read_count.min(1000);
        for i in 0..n {
            fastq.push_str(&format!(
                "@read_{i}\n{}\n+\n{}\n",
                "A".repeat(config.read_length),
                "I".repeat(config.read_length)
            ));
        }
        std::fs::write(&reads, fastq)
            .map_err(|e| gatk_common::GatkError::io("Failed to write reads", e))?;
        Ok(())
    }
}

/// On-disk size summary for a named dataset directory.
/// # Invariants
/// `reference_size` and `reads_size` are heuristic half splits of `total_size`.
/// # Ownership
/// Owns name string; clone for reporting.
/// # Mutation
/// Snapshot from [`DatasetManager::get_dataset_stats`].
/// # Biological assumptions
/// None (filesystem metrics).
/// # Java equivalence
/// None / Rust-native.
#[derive(Debug, Clone)]
pub struct DatasetStats {
    pub name: String,
    pub reference_size: u64,
    pub reads_size: u64,
    pub total_size: u64,
    pub is_valid: bool,
    pub last_modified: Option<chrono::DateTime<chrono::Utc>>,
}

fn dir_size(path: &Path) -> u64 {
    let mut total = 0u64;
    if let Ok(entries) = std::fs::read_dir(path) {
        for entry in entries.flatten() {
            let p = entry.path();
            if p.is_dir() {
                total += dir_size(&p);
            } else if let Ok(meta) = entry.metadata() {
                total += meta.len();
            }
        }
    }
    total
}
