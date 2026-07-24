//! `HaplotypeBAMWriter` parity stub.
//! **Deferred product feature (Sprint G):** not wired to production HC CLI.
//! See `docs/ARCHITECTURE.md` (T5-1).

use gatk_common::GatkResult;
#[cfg(any(feature = "dev-dumps", test))]
use std::io::Write;

/// Enable flag for assembled-read BAM output (product stub).
/// # Invariants
/// When `enabled` is false, [`BamoutWriter`] performs no I/O.
/// # Ownership
/// [`Copy`] / cloneable config passed into [`BamoutWriter::new`].
/// # Mutation
/// Callers set `enabled` before constructing a writer; writer does not mutate config.
/// # Biological assumptions
/// None — output formatting stub, not genomics logic.
/// # Java equivalence
/// GATK `HaplotypeBAMWriter` configuration slice (; deferred product feature).
#[derive(Debug, Clone, Default)]
pub struct BamoutWriterConfig {
    pub enabled: bool,
}

/// Placeholder BAM writer counting assembled-read emissions for parity gates.
/// # Invariants
/// `records_written` monotonically increases on each successful placeholder write.
/// # Ownership
/// Caller-owned writer instance; no underlying file handle yet.
/// # Mutation
/// [`Self::write_assembled_read_placeholder`] mutates `records_written` in place.
/// # Biological assumptions
/// None — stub until real BAM serialization is wired.
/// # Java equivalence
/// GATK `HaplotypeBAMWriter` parity stub (not production HC CLI).
#[derive(Debug, Clone, Default)]
pub struct BamoutWriter {
    pub records_written: usize,
}

impl BamoutWriter {
    pub fn new(_config: &BamoutWriterConfig) -> Self {
        Self { records_written: 0 }
    }

    pub fn write_assembled_read_placeholder(&mut self) -> GatkResult<()> {
        self.records_written += 1;
        Ok(())
    }
}

/// J-D04 — stable bamout gate dump.
#[cfg(any(feature = "dev-dumps", test))]
pub fn dump_bamout_stub_tsv(
    enabled: bool,
    write_count: usize,
    out: &mut impl Write,
) -> GatkResult<()> {
    let mut writer = BamoutWriter::new(&BamoutWriterConfig { enabled });
    for _ in 0..write_count {
        writer.write_assembled_read_placeholder()?;
    }
    writeln!(out, "enabled\t{enabled}")?;
    writeln!(out, "records_written\t{}", writer.records_written)?;
    writeln!(out, "writer_ready\ttrue")?;
    Ok(())
}
