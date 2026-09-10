# Experimental SeqGraph k-best resource policy (6R.140)

**Experimental instrumentation only. No k-best memory contract has been selected** (6R.139). **`docs/CLAIM_MATRIX.md` is unchanged.** The production default remains the historical **1024 `BinaryHeap.len()` Peak-RSS proxy**.

This is not a product claim, CLI feature, or approved byte budget.

## What is separated

| axis | meaning | production default |
|------|---------|--------------------|
| **K** | requested completed paths / per-vertex visit budget | unchanged (`num_best_haplotypes_per_graph`, typically 128) |
| **Resource policy** | how much **frontier** may be live | `legacy_1024` |

Do not equate K with the 1024 heap-object proxy.

## Policies

| name | frontier refusal | status |
|------|------------------|--------|
| `legacy_1024` | `BinaryHeap.len() >= 1024` (same timing as pre-6R.140 SeqGraph k-best) | **production default** |
| `unbounded_diagnostic` | none (path-edge 4096 and expansion 12_000 Peak-RSS caps **remain**) | diagnostic only |
| `byte_budget` | experimental live-frontier **bytes** ≥ caller ceiling | diagnostic; **no ceiling is a product budget** |

## Env (internal / experimental)

Unset = default `legacy_1024`. Invalid values **fail closed** (configuration error). Not a public CLI.

| variable | role |
|----------|------|
| `GATK_RS_EXPERIMENTAL_KBEST_POLICY` | `legacy_1024` \| `unbounded_diagnostic` \| `byte_budget` |
| `GATK_RS_EXPERIMENTAL_KBEST_MEMORY_BUDGET` | positive integer **bytes**; required iff `byte_budget`; forbidden otherwise |
| `GATK_RS_KBEST_DIAGNOSTICS` | `1` → stderr line `KBEST_RESOURCE_DIAG ...` |

## Accounting (experimental, copied `PathState`)

Counted in `peak_frontier_bytes` (and the byte-budget check):

- `heap.len() × size_of::<HeapItem>()` (includes `PathState` / `Vec` header)
- `Σ edges.capacity() × 16` (owned edge payload)

**Not** counted: malloc rounding, `BinaryHeap` capacity slack (`peak_heap_vec_capacity_bytes` is diagnostic only), completed result paths, SeqGraph, PairHMM, BAM, fragmentation.

Process RSS (`rss_before` / `rss_peak` / `rss_after`) is sampled around the run. **It is not k-best-attributable** merely because it was recorded (6R.138).

## Logical vs heap vs bytes

A future compact representation must **not** redefine the resource contract as `BinaryHeap.len()`.

| representation | logical frontier | heap entries |
|----------------|------------------|--------------|
| copied / compact eager | live unpopped partials | same |
| compact lazy heads | live unpopped **+ deferred siblings** | heap heads only |

6R.136: logical peak 1812 with lazy heap peak 1529.

## API

- Production: `find_best_haplotypes_seq_graph` → `SeqKbestResourcePolicy::from_runtime()` (default legacy).
- Tests/diagnostics: `find_best_haplotypes_seq_graph_with_policy`.
- `DIAGNOSTIC_GENEROUS_BYTE_BUDGET` (8 MiB) is **only** for 6R.140 canonical experiments. It is not a selected contract.

## Next

Curves across budgets/workloads are **6R.141**, only if this instrumentation stays trustworthy. Do not treat 8 MiB, 839 KiB, or 1812 states as policy.
