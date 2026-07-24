//! E-D01 — assembly failure / graph debug dumps (scaffold).

use gatk_common::GatkResult;
use std::io::Write;

pub fn dump_assembly_debug_stub_tsv(
    assembly_failure_bam: bool,
    graph_dot: bool,
    out: &mut impl Write,
) -> GatkResult<()> {
    writeln!(out, "assembly_failure_bam_enabled\t{assembly_failure_bam}")?;
    writeln!(out, "graph_dot_enabled\t{graph_dot}")?;
    writeln!(out, "dump_ready\tfalse")?;
    Ok(())
}
