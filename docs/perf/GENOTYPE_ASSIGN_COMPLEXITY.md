# Genotype assignment complexity & dense-region wall

## Evidence (why this leaf)

Local TRACE (`docs/perf/BEAT_JAVA_WALL_NEXT.md`):

| Window | `assign_genotype` | PairHMM | Realign |
|--------|------------------:|--------:|--------:|
| mega `21:9825–9828k` | **~130–172 s** | ~9–20 s | ~4–7 s |
| hot locus `21:9826233` | **~37 s** alone | — | — |

Assign dominates PairHMM by ~10× on mega dense windows. Product wall-losers
still lose on genotype/AD structure, not PairHMM packs.

## Exact complexity (before optimization)

Symbols:

| Symbol | Meaning |
|--------|---------|
| \(S\) | genotyped start positions in the active window |
| \(A_s\) | merged biallelic alleles at position \(s\) (strict Java: **not** truncated) |
| \(R\) | pileup / likelihood reads |
| \(H\) | haplotypes |
| \(P\) | full pileup AD / softclip passes per allele (**≈5–10** in `try_genotype`) |
| \(C\) | CIGAR walk cost per read/hap |
| \(L\) | PairHMM DP length product (paid **once** per region) |

### Per allele (one `try_genotype_variation_event`)

| Stage | Complexity | Notes |
|-------|------------|-------|
| Allele map (`SiteMap` / `createAlleleMapper`) | \(O(H\cdot C)\) | EventMap cache helps; SNP still CIGAR base walks |
| Pileup AD (multi-pass) | \(O(P\cdot R\cdot C)\) | **Primary hotspot** |
| Softclip alt counts | \(O(R\cdot C)\) × (often 2: early + main) | Duplicate scans |
| Likelihood reshape + marginalize | \(O(R\cdot H)\) | `region_likelihoods_to_rows` + max over pools |
| Diploid PL / genotype states | \(O(R)\) → **3** values | **Not** the bottleneck |

### Region total

\[
T_{\text{geno}} = O(H\cdot C)_{\text{EventMap once}}
  + \sum_{s=1}^{S}\Big(
      O(A_s\log A_s\cdot H\cdot C)_{\text{allele sort}}
      + \sum_{a=1}^{A_s} O(P\cdot R\cdot C + H\cdot C + R\cdot H)
    \Big)
\]

PairHMM (once): \(O(R\cdot H\cdot L)\).

**Superlinearity vs PairHMM:** density raises \(S\), often \(A_s\), and coverage \(R\)
together, while PairHMM does not pay \(P\cdot S\cdot A\). Multi-allelic sites become
**\(A\) separate biallelic trials** (each PL length 3), not
\(\binom{A+\text{ploidy}-1}{\text{ploidy}}\) genotype states.

### What is *not* the bottleneck

| Question | Answer |
|----------|--------|
| Diploid genotype enumeration | Always **3** states (0/0, 0/1, 1/1) per biallelic trial |
| Sparse genotype-state reps | Not beneficial — state space is tiny |
| Caching genotype-index maps | Already trivial (fixed diploid formula) |
| SIMD / dense PL kernels | Marginal vs AD; PL is \(O(R)\) |

## Per-site profile fields (`GATK_RS_HC_PROFILE`)

| Field | Meaning |
|-------|---------|
| `candidate_alleles` | REF+ALT for the call (biallelic → 2) |
| `genotype_states` / `pl_vector_len` | Diploid PL length (typically 3) |
| `samples` | HC single-sample → 1 |
| `ad_wall_s` | Pileup AD / softclip |
| `allele_map_wall_s` | `SiteMap` / mapper |
| `marginalize_wall_s` | AlleleLikelihoods max pools |
| `genotype_enum_wall_s` | PL / posterior from marginalized rows |
| `event_rebuild_wall_s` | Region EventMap build |

## Optimizations landed (semantics preserved)

1. **AD result memo** (`ad_result_memo.rs`) — identical `(reads, loc, pad, alleles, mode)`
   returns cached `(ref,alt)` / softclip counts. Collapses \(P\) rescans to 1 when
   early-template + main path repeat the same scan.
2. **Per-locus hap SNP base precompute** before allele sort — \(O(H)\) once vs
   \(O(A\log A\cdot H)\) CIGAR walks.
3. Nested profile counters for map / marginalize / PL-enum / AD.

Decode cache (`AdDecodeCache`) remains: avoids CIGAR/seq realloc; memo avoids
re-projection.

## Adversarial benches

```bash
cargo bench -p gatk-haplotypecaller --bench genotype_assign --features parity_harness --locked
```

Groups:

- `genotype_dense_ad/multipass_ad/R*_S*_A*_P*` — mega-shaped multi-pass AD
- `genotype_marginalize_pl/marg_then_pl/*` — show PL/marg ≪ AD at scale

## Prove

```bash
cargo test -p gatk-haplotypecaller --lib ad_result_memo --locked
cargo test -p gatk-haplotypecaller --lib hc_profile --locked
cargo test -p gatk-haplotypecaller --test hc_genotyping_marginalization_fixture_test --locked
```

## Deferred

- Collapse `try_genotype` control flow to a **single AD authority** (structural;
  higher parity risk than memo)
- Share `region_likelihoods_to_rows` across alleles at one locus
- Full-pad indel mapper EventMap reuse
