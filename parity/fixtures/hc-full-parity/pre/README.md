# Phase PRE — read preparation before assembly

Producer: `hc_full_parity_gate_dump read-pre-softclip <sam|bam> [dont_use_soft_clipped_bases] [override_softclip_fragment_check]`

Default args match HC `AssemblyBasedCallerArgumentCollection`: `0` / `0`.

Implements the first soft-clip branch in `AssemblyBasedCallerUtils.finalizeRegion`:

- `hard_clip` when `dontUseSoftClippedBases` **or** fragment size is not well-defined
- `revert` otherwise (`ReadClipper.revertSoftClippedBases` + `os`/`oe` tags on revert path)

Schema: `qname`, `flags`, `fragment_length`, `cigar_in`, `cigar_out`, `seq_len_in`, `seq_len_out`, `action`, `os`, `oe`.
