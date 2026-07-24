//! `ReadThreadingAssemblerArgumentCollection` defaults dump.

use crate::read_threading_assembler::ReadThreadingAssemblerArgs;
use gatk_common::{GatkError, GatkResult};
use std::io::Write;

pub fn dump_assembler_args_tsv(out: &mut impl Write) -> GatkResult<()> {
    let a = ReadThreadingAssemblerArgs::default();
    writeln!(out, "arg\tvalue").map_err(|e| GatkError::generic(format!("write: {e}")))?;
    macro_rules! row {
        ($k:expr, $v:expr) => {
            writeln!(out, "{}\t{}", $k, $v)
                .map_err(|e| GatkError::generic(format!("write: {e}")))?;
        };
    }
    row!("kmer_sizes", format!("{:?}", a.kmer_sizes));
    row!("min_base_quality", a.min_base_quality);
    row!("min_prune_factor", a.min_prune_factor);
    row!("use_adaptive_pruning", a.use_adaptive_pruning);
    row!("prune_before_cycle_counting", a.prune_before_cycle_counting);
    row!("recover_dangling_branches", a.recover_dangling_branches);
    row!("recover_dangling_heads", a.recover_dangling_heads);
    row!(
        "recover_all_dangling_branches",
        a.recover_all_dangling_branches
    );
    row!("min_dangling_branch_length", a.min_dangling_branch_length);
    row!(
        "allow_non_unique_kmers_in_ref",
        a.allow_non_unique_kmers_in_ref
    );
    row!(
        "dont_increase_kmer_sizes_for_cycles",
        a.dont_increase_kmer_sizes_for_cycles
    );
    row!("allow_low_complexity_graphs", a.allow_low_complexity_graphs);
    row!(
        "remove_paths_not_connected_to_ref",
        a.remove_paths_not_connected_to_ref
    );
    row!("use_seq_graph", a.use_seq_graph);
    row!("ensure_reference_in_result", a.ensure_reference_in_result);
    row!(
        "num_best_haplotypes_per_graph",
        a.num_best_haplotypes_per_graph
    );
    Ok(())
}
