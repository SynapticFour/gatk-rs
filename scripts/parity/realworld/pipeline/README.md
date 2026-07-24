# Real-World HC pipeline runner

**Definition (read this first):** [`docs/ARCHITECTURE.md`](../../../../docs/ARCHITECTURE.md)

**Runner:**

```bash
export RW_REF="$PWD/parity/realworld/assets/hs37d5.simple.fa"
export RW_BAM="$PWD/parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam"
export RW_INTERVAL="20:413419-463418"

# Fast path (ingress + counts only):
RW_SKIP_STEP06=1 RW_SKIP_STEP07=1 ./scripts/parity/realworld/pipeline/run_paired_realworld_pipeline.sh

# Full path (includes long Docker HC — can take hours):
./scripts/parity/realworld/pipeline/run_paired_realworld_pipeline.sh
```

Output: under `parity/reports/realworld_pipeline_run/` (`summary.md`, `equivalence_report.{md,json}`, append-only `pipeline_footer.txt`).

**Strict smoothed float compare (debug only):** `RW_SMOOTHED_ACTIVITY_STRICT=1` adds `--require-continuous-max-diff` to step 06 (usually **fails** vs tri-state IGV; default contract is **binary** activity parity — see `docs/CLAIM_MATRIX.md`).

**Fast regression on existing `OUT_DIR` (no Docker):** regenerates the machine report and asserts JSON gates:

```bash
chmod +x scripts/parity/realworld/pipeline/run_realworld_equivalence_selfcheck.sh
./scripts/parity/realworld/pipeline/run_realworld_equivalence_selfcheck.sh
```

**Single step:** `RW_ONLY_STEP=3 ./scripts/parity/realworld/pipeline/run_paired_realworld_pipeline.sh`
