# Air M4 16 GB GIAB recipe (memory + disk)

Local gold-standard path for `GIAB_MODE=ci-subset` on a constrained laptop.
Does **not** invent non-Java downsampling; Peak-RSS wins come from shared reads, sequential regions, and TLS scratch release.

## Environment

```bash
export GATK_RS_HC_SEQUENTIAL=1          # one assembly region at a time
export RAYON_NUM_THREADS=1
export GATK_RS_HC_LARGE_REGION_READS=4096
export GIAB_KEEP_SHARD_VCFS=0           # delete shard VCFs after concat
export GIAB_KEEP_EQUIV_INTERMEDIATES=0  # delete union BAM + strat beds after equiv
export GIAB_HC_WINDOW_BP=1000000        # start at 1 Mb; raise only after RSS proof
```

## Disk hygiene

- Stage reference + truth **once** under `parity/realworld/assets/` / GIAB cache.
- Prefer window BAM slices over full WGS downloads.
- Finalize hardlinks `hc/*.vcf` → `equiv/` when possible (no duplex).
- Wipe `parity/giab/runs/ci/**` between modes when free space is tight.

## Honesty

- Full HG001 30× P12 on this host remains a dedicated-RAM gate (see `HC_MEMORY_PROFILE.md`).
- Smoke uses NA12878_20k for P12; do not cite smoke as `ci-subset`.
- Do not claim genome-wide equivalence from fitting in swap.
