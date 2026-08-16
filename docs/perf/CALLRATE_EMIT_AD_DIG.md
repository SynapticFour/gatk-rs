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

1. ~~Land emit-gate fix; re-dispatch **ci-subset**~~ — #116 on main; ci-subset
   `max_|ΔF1|=0.0099` PASS; wall-losers `max_|ΔF1|=0.0004` PASS
   ([31940745824](https://github.com/SynapticFour/gatk-rs/actions/runs/31940745824)).
2. **Stage-classify remaining wall-losers Java-only** (post-#116): of ~9.6k
   Java-only sites, **~6.9k** are strong SNP hets (`alt≥4`, `dp≥10`) with **no Rust
   site at the same POS** — not emit-depth failures. Highest leverage is
   discovery / EventMap / genotyping-empty (`supplement_assembly_events_from_reads`,
   EventMap retain), not lowering emit thresholds.
3. Holdout residual 6 (sparse pileup / indel / anomalous GQ) remain follow-ups;
   do not widen P12 bands or blindly loosen `alt_ad≥4` / GQ gates (Rust-only already
   ~5k on wall-losers).

### Wall-losers stage hypothesis (31940745824)

| Class | Approx | Implication |
|-------|-------:|-------------|
| Strong SNP het, no Rust POS | ~6917 | Discovery/EventMap undercall |
| Other Java-only | ~2726 | Mixed indel / weak / allele-diff |
| Rust-only | ~5165 | Overcall / allele reshape — watch when raising recall |

Artifacts: `parity/giab/runs/assemble-wall-campaign/wall-losers-31940745824/`,
local TRACE dig under `callrate-stage-dig/` when BAMs available.

### Local stage TRACE (`21:9411500-9414600`, product thr=2)

12 strong Java-only SNP hets from wall-losers (Java AD clear hets):

| Bucket | n | Example |
|--------|--:|---------|
| `emit_skip_non_p12_support` | 2 | `9411732` AD=11/3 (Java 31,13); `9412808` read_AD=14/3 GQ=3 (Java 41,12) |
| Absent from TRACE at POS | 10 | Never discovered/genotyped in this window |
| Emitted | 0 | — |

So the bulk miss class is **discovery/EventMap empty**, with a smaller residual
**pileup/FORMAT AD undercount** that fails `alt_ad≥4` even when Java DepthPerAllele
is strong. Next code bet: evidence-class EventMap/SNP retain in
`supplement_assembly.rs` — not lowering emit thresholds.