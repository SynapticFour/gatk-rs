# HC full parity — fixture subtree

Artifacts used by **Phase B** (walker / assembly regions) and future `docs/ARCHITECTURE.md` gates.

| Path | Role |
|------|------|
| [`b1/`](./b1/README.md) | Read-shard padding + merge semantics (`make_read_shards`). |
| [`b2/`](./b2/README.md) | `AssemblyRegionIterator` TSV output on tiny REF+BAM. |
| [`b3/`](./b3/README.md) | B.3: `apply` count + inactive `callRegion` fast-path stats (`WalkerApplyStats`). |
| [`b4/`](./b4/README.md) | B.4: full walker traversal (`traverse_assembly_region_walker`) — same stats schema as B.3. |
| [`c1/`](./c1/README.md) | C.1: per-locus raw `ActivityProfileState` (pre–band-pass). |
| [`c2/`](./c2/README.md) | C.2: band-pass smoothed activity. |
| [`c3/`](./c3/README.md) | C.3: per-locus binary `is_active`. |

**Pin:** Java comparisons (when added) must use the GATK build recorded in [`docs/GATK_PINNED.env`](../../../docs/GATK_PINNED.env).

**Phase B (required before Phase C):** `./scripts/parity/run_hc_full_parity_phase_b.sh` — runs `cargo test -p gatk-haplotypecaller` + B.1–B.4 gates.

**Phase C:** `./scripts/parity/run_hc_full_parity_phase_c.sh` (runs Phase B, then C.1–C.5 incl. genotype-likelihoods and multisample raw activity).

**Phase D:** `./scripts/parity/run_hc_full_parity_phase_d.sh` (runs Phase C, then D.1–D.3).

**Harness (individual):** `run_hc_full_parity_b{1,2,3,4}_*.sh`, `run_hc_full_parity_c{1,2,3}_*.sh`.
