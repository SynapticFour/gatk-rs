# Phase B.4 — AssemblyRegionWalker traversal

L1 gate: `scripts/parity/run_hc_full_parity_b4_walker_traversal.sh`

Runs `hc_full_parity_gate_dump walker-traversal-summary`, which uses `traverse_assembly_region_walker` (all read shards → per-shard `AssemblyRegionIterator` → aggregated `WalkerApplyStats`).

Goldens are byte-identical to B.3 apply-summary: traversal must not change apply/inactive/active counts vs the single-shard drain path when downsampling is disabled (default).
