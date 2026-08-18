# Absolute Rust leaf campaign (ci-subset HC VCF)

## Scope

Optimize **every** Rust function on the GIAB ci-subset **production VCF** path for
absolute wall on the worst product shard — not Java ratio. Keep algorithm parity
gates; no P12 band widening.

| Item | Value |
|------|------:|
| Workload | `01_chr21_w09` · `21:9000001-10000000` BAM |
| Microbench | mega `21:9825000-9828000` thr=2 |
| Branch | `perf/absolute-rust-leaves` |
| All fns (HC crates) | **3096** → `hc_path_all_functions.tsv` |
| VCF-path filtered | **2786** → `hc_vcf_path_functions.tsv` |

## Baseline TRACE (mega, this campaign)

Pre-dedupe-opt:

| Phase | Σ seconds |
|-------|----------:|
| assign_genotype | **123.1** |
| pairhmm | 8.4 |
| realign | 2.6 |

**sample (gatk-rs PID) top-of-stack:** `dedupe_likelihood_subset_by_qname` ≫ PairHMM
`score_one_f64` ≫ SW. That leaf was O(|subset|²) via repeated max-LL scans.

### Cut 1 — `dedupe_likelihood_subset_by_qname` ✅

Precompute max LL per read_index; HashMap best-per-QNAME; same tie-break
(higher LL, then lower index). Mega rematch: **assign 123.1 → ~0.5–1.1 s** (~100×+);
sites unchanged (62).

### Cut 2 — PairHMM `score_one_f64` ✅ (partial)

- Inline priors (no prior plane writes)
- Qual→match/mismatch LUT (no per-cell `powf`)
- Word-at-a-time `first_hap_divergence`
- Rolling 2-row DP for leftover singles
- Skip full-matrix memset on fresh haps (overwrite-only)
- Kept hapStartIndex prefix reuse for equal-length chains (rolling-all regressed)

Laptop TRACE variance on mega pairhmm is high (~6.5–14 s); use repeated runs /
CI wall-losers for signed claims. Sample still tops at `score_one_f64` then SW.

### Cut 3 — SW `calculate_matrix` ✅

Unchecked indexing in the DP hot loop (bounds proven by refuse-oversized gate).

### Cut 4 — SNP discovery span clamp ✅

`discover_snp_events_from_reads` no longer walks full pad×reads; clamps to
read∩active offsets, HashMap counts, AdDecodeCache for CIGAR/seq.

## Method

1. Inventory (done) — static fn census + pipeline map.
2. Profile deepest leaves on mega (`sample` / Instruments on **gatk-rs PID**, not shell).
3. Isolate one leaf → Rust-native opt → rematch mega TRACE + unit gates.
4. Walk inventory by measured time, then by remaining unoptimized modules.

## Leaf cohort order

1. Genotype site pipeline / AD pileup / softclip decode
2. af_calc / PL / allele mapping
3. PairHMM SIMD
4. SW realign + EventMap
5. Assembly / kbest / dangling
6. Activity / BAM / emit

## Full isolation audit (plan)

Protocol: every VCF-path fn → `FN_VERDICTS.tsv` (`opt`/`tight`/`defer`/`n/a`);
`AUDIT_BOARD.tsv` file status → `audited`. Micro-opts only when same math.

| Artifact | Status |
|----------|--------|
| `FN_VERDICTS.tsv` | **2786** rows (0 pending) |
| `AUDIT_BOARD.tsv` | **156/156** audited |
| Verdict mix | opt **24** · tight **2711** · n/a **51** |

### Stage rematch — mega `21:9825000-9828000` thr=2

After stages 01–04 (+05–09 audit) safe opts (SNP CIGAR-walk candidate **reverted** —
unsigned off wrap changed site set; SW border-only fill **reverted** pending deeper
proof):

| Metric | Value |
|--------|------:|
| sites | **62** (exact allele match vs baseline-mega) |
| assign_genotype Σ | ~0.9 s |
| pairhmm Σ | ~6–10 s (laptop variance) |

### Additional opts this pass (parity-safe)

- Genotype: BTreeSet→HashSet for QNAME/index membership; O(n²) softclip subset → HashSet
- Allele map: `hap_base_at_ref_locus` CIGAR span bulk-advance
- PairHMM scalar Logless: skip full memset + inline prior LUT (mirror pack)
- PCR: `OnceLock` Conservative/Aggressive qual caches
- EventMap: `prefer_indel_over_colocated_snps` HashSet; indel support HashMap
- Rejected after rematch: dense SNP CIGAR walk; SW skip-fill

Helpers: `stage_audit.py`, `mark_file_audited.py`.

## Follow-on pass (NEON/AVX2 + SW + SNP attempt)

| Item | Result |
|------|--------|
| NEON `score_pack2` | skip full SoA memset + inline priors; SIMD↔scalar green |
| AVX2 `score_pack4` | equal-length kernel (no lane masks); same |
| SW `calculate_cigar` | unchecked max-score / backtrack loads |
| `haplotype_cigar` pad | TLS reuse for padded ref/alt |
| SNP CigarWalk v2 | **deferred** — synthetic golden OK, mega allele set drifted (69–76); production stays HEAD query_index full-pad |
| Mega after follow-on | **sites=62** exact vs baseline |
| w09 full shard | exit 0; sites=5207; TRACE pairhmm **37.9s** ≫ assign **5.4s**; wall ~197s; Peak ~472 MiB |

### Cut 5 — post-PairHMM prep membership ✅

Replace O(n²) `candidates.iter().any(events_match)` / merge scans with
`HashSet` keys; `mem::take` prior events (no clone); indel discovery uses
`AdDecodeCache`; `apply_cigar_to_cigar` run-length emit.

| Metric (mega thr=2) | Value |
|---------------------|------:|
| sites | **62** exact vs baseline-mega |
| post_pairhmm_prep Σ | **3.17 s** (`prep_realign` 2.94 · `prep_parity_spine` **0.13**) |
| pairhmm Σ | ~6.1 s |
| genotype Σ | ~3.7 s |

