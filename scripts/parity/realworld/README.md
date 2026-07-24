# Real-world staged scripts

Numbered helpers for `docs/ARCHITECTURE.md`. Run from repository root:

```bash
./scripts/parity/realworld/01_check_environment.sh
# …
```

`03_stage_reference_and_truth.sh` sets `REALWORLD_STOP_AFTER_ASSETS=1` and calls `run_p12_p13_realworld_full.sh`, which exits after staging `parity/realworld/assets/` (see that script for details).

**Evidence (foundation + report):**

```bash
./scripts/parity/realworld/run_foundation_evidence.sh
# optional: + P13 on current VCFs
export P12_INTERVAL='20:413419-463418'
./scripts/parity/realworld/run_full_evidence.sh
```

Writes `parity/reports/realworld_parity_evidence.{md,json}`.
