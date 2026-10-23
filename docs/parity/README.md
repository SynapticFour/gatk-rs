# Parity investigation notes

These files are the engineering archive. They are not the acceptance policy
and they do not add a **Yes** row by existing.

Start with [`../PARITY.md`](../PARITY.md) for how gates are defined, and
[`../CLAIM_MATRIX.md`](../CLAIM_MATRIX.md) for what the project asserts.
This directory explains *why a particular site diverged* and what, if
anything, was changed afterward.

## Names

`<series>.<number>_<TOPIC>.md`

- **6R** — HaplotypeCaller investigations that led to the mid-B milestone and
  the chr20 / P12 site work through the frozen witness in
  [`6R.312_FINAL_CLOSURE.md`](6R.312_FINAL_CLOSURE.md) and the holdouts in
  [`6R.314_HOLDOUT.md`](6R.314_HOLDOUT.md). The long narrative also sits in
  the engineering-record half of [`../PARITY.md`](../PARITY.md).
- **7R** — later investigations (k-best frontier, annotation membership,
  post-PairHMM, realignment). A note that says the production code was left
  unchanged is a finding, not a passing gate.
- Numbering is investigation order. A missing number is a missing note, not
  a skipped release. Some 7R notes are local working files and may be
  untracked; treat them as evidence in progress, not as policy.

TSV files beside the notes (`6R.141_KBEST_ORACLE.tsv` and similar) are
frozen comparison tables for those investigations.

## Current policy note inside this directory

[`KBEST_RESOURCE_POLICY.md`](KBEST_RESOURCE_POLICY.md) and
[`7R.22_JAVA_EQUIVALENT_KBEST_DEFAULT.md`](7R.22_JAVA_EQUIVALENT_KBEST_DEFAULT.md)
are the current description of the SeqGraph k-best resource policy.
When `GATK_RS_EXPERIMENTAL_KBEST_POLICY` is unset, the policy is
`unbounded_diagnostic`. `legacy_1024` is selected only by setting that
value explicitly. The claim matrix was not given a new Yes row for this
default.

These notes still say, in the present tense, that an unset variable selects
`legacy_1024`. That was the observation when they were written. It is not
the current policy:

- [`7R.18_KBEST_RESOURCE_POLICY_COMPOSITION.md`](7R.18_KBEST_RESOURCE_POLICY_COMPOSITION.md)
- [`7R.19_KBEST_RESOURCE_ABORT_CONTRACT.md`](7R.19_KBEST_RESOURCE_ABORT_CONTRACT.md)
- [`7R.20_RSS_ABORT_ERROR_CONTRACT.md`](7R.20_RSS_ABORT_ERROR_CONTRACT.md)
- [`7R.21_UNBOUNDED_KBEST_PRODUCTION_COUNTERFACTUAL.md`](7R.21_UNBOUNDED_KBEST_PRODUCTION_COUNTERFACTUAL.md)

Older sentences in the 6R log that call `legacy_1024` the default are
historical in the same way. Do not edit those notes to match the later
default.

## How to record a new divergence

Follow [`../PARITY.md`](../PARITY.md): name the Java contract, the observed
outputs, and whether production code changed. Do not reclassify an
exploratory note as a formal acceptance result, and do not widen a P12 band
to make a comparison pass.
