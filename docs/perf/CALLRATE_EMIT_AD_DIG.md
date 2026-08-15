# Call-rate dig — ci-subset ΔF1 + L8 holdout (post hang-fix)

## Inputs

- Failed ci-subset [31880000273](https://github.com/SynapticFour/gatk-rs/actions/runs/31880000273):
  `max_|ΔF1|=0.0249` (SNP); INDEL within 0.02.
- Artifacts under `parity/giab/runs/ci-subset-31880000273/`.

## Site-diff (concat VCF keys)

| | Count |
|--|------:|
| Java sites | 200890 |
| Rust sites | 137793 |
| Shared | 123916 |
| Java-only | 76974 (**70%** true misses, not allele swaps) |
| rust/java | **0.686** |

Java-only is **SNP-heavy** (66734 SNP / 10240 indel), concentrated on full chr20/21
(dense clusters e.g. `21:10698011-10726651` n=1774, `20:29.5–29.6Mb`).

Holdout window inside the same concat: java=55 rust=41 rate **0.745** (matches
`docs/perf/L8_HOLDOUT_F1_TRACK.md`).

## Root cause (holdout TRACE)

Miss SNP `20:15009054` **is genotyped** (`after_genotype calls=1`) then dropped in
`try_emit_call_region_variants` by `strict_java_non_p12_region_supports_emit`:

| Context | reads | FORMAT AD | Result |
|---------|------:|----------:|--------|
| Isolated `-L` ~1 kb | 84 | 12/16 (Java-like) | **emitted** |
| Full 50 kb holdout | 68 | **2/1** GQ=99 | **emit_skip_non_p12_support** |

Gate required het `alt_ad≥4 && dp≥10 && GQ≥30` using **FORMAT AD only**. PairHMM
“informative” AD undercounts vs Java `DepthPerAlleleBySample` / region pileup, so
confident hets are suppressed. Not PairHMM SIMD; not cycle-strip memo.

## Fix

In `region_vcf_emit.rs`: for the non-P12 support gate, take
`max(FORMAT AD, region-read pileup AD)`; if pileup alone clears the depth bar,
do not also demand GQ≥30.

## Proof (local holdout)

| | Before | After |
|--|-------:|------:|
| rust/java sites | 0.745 | **0.927** |
| L8 p13 gate | fail F1≈0.844 | **pass** rust F1≈0.904 |

Remaining java-only (6): sparse/weak pileup or indel class — follow-up, not blocking
L8 threshold.

## Next

1. Land emit-gate fix; re-dispatch **ci-subset** (expect |ΔF1| closer to 0.02).
2. **wall-losers** rematch already dispatched: [31884065981](https://github.com/SynapticFour/gatk-rs/actions/runs/31884065981).
3. Wall leaf (PairHMM/EventMap) from TRACE after product-wall numbers land.
