# Dense loser TRACE recipe (phase5 wall campaign)

Product-shaped wall measurement for NA12878 CI loser windows. **Do not** use
`GATK_RS_HC_SEQUENTIAL=1` for wall pies (that is Peak-only). Keep production
PairHMM default as configured; for wall campaigns pass `--pair-hmm FASTEST_AVAILABLE`.

## Env

```bash
export RAYON_NUM_THREADS=2
export GATK_RS_HC_RSS_TRACE=1
# unset GATK_RS_HC_SEQUENTIAL
BIN=target/release/gatk-rs
REF=parity/realworld/assets/hs37d5.simple.fa
OUT=parity/giab/runs/assemble-wall-campaign/phase5-rematch
mkdir -p "$OUT"
```

## Windows

```bash
# w09 hot 200kb
"$BIN" HaplotypeCaller -R "$REF" \
  -I parity/realworld/na12878_ci_loser_windows/01_chr21_w09.bam \
  -O "$OUT/w09_200kb.vcf" -L 21:9500000-9700000 --threads 2 \
  --pair-hmm FASTEST_AVAILABLE \
  >"$OUT/w09_200kb.stdout" 2>"$OUT/w09_200kb.trace"

# w11 head ~200kb
"$BIN" HaplotypeCaller -R "$REF" \
  -I parity/realworld/na12878_ci_loser_windows/01_chr21_w11.bam \
  -O "$OUT/w11_200kb.vcf" -L 21:11000000-11200000 --threads 2 \
  --pair-hmm FASTEST_AVAILABLE \
  >"$OUT/w11_200kb.stdout" 2>"$OUT/w11_200kb.trace"

# w26 head ~200kb
"$BIN" HaplotypeCaller -R "$REF" \
  -I parity/realworld/na12878_ci_loser_windows/00_chr20_w26.bam \
  -O "$OUT/w26_200kb.vcf" -L 20:26000000-26200000 --threads 2 \
  --pair-hmm FASTEST_AVAILABLE \
  >"$OUT/w26_200kb.stdout" 2>"$OUT/w26_200kb.trace"
```

## Summarize

```bash
python3 scripts/parity/giab/summarize_hc_rss_trace_wall.py "$OUT/w09_200kb.trace"
```

Compare to [`PHASE5_WALL_BASELINE.md`](PHASE5_WALL_BASELINE.md).
Use **hc-mem-probe** for Peak RSS in CI; ignore `/usr/bin/time` Java ~30 MiB.
