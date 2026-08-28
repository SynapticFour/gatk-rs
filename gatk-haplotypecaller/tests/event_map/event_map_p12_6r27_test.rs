//! 6R.27 TEST-ONLY: diagnose `call_region(None)` after a Java-equivalent EventMap.
//! Does not change production algorithms, EventMap, P12, W-H1, or k-mer policy.

#[cfg(test)]
mod traces {
    use crate::assembly_based_caller::{assemble_reads_with_finalized, AssembleReadsArgs};
    use crate::assembly_region_finalize::{
        assembly_reference_read, create_graph_reference_read, padded_reference_loc,
    };
    use crate::assembly_region_trimmer::{
        AssemblyRegionTrimmer, AssemblyRegionTrimmerConfig, TrimVariant,
    };
    use crate::engine::{take_call_region_audit, CallRegionArgs, HaplotypeCallerEngine};
    use crate::genome_loc::GenomeLoc;
    use crate::haplotype::Haplotype;
    use crate::haplotype_cigar::get_bases_covering_ref_interval;
    use crate::read_event_discovery::P12_CLUSTER_TTC_START;
    use crate::read_threading_assembler::assemble_from_ref_and_reads;
    use crate::walker_apply::{call_disposition, AssemblyRegionCallDisposition};
    use gatk_core::reference::{
        parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary,
    };
    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};
    use std::path::Path;

    const SITE_CT: u64 = 92_317_361;
    const SITE_CG: u64 = 92_317_371;
    const SITE_CA: u64 = 92_317_399;
    const SITE_TC: u64 = 92_317_407;
    const SITE_GC: u64 = 92_317_412;
    const JAVA_ACTIVE: (u64, u64) = (92_317_262, 92_317_491);
    const JAVA_EXTENDED: (u64, u64) = (92_317_162, 92_317_591);
    const SNP_PAD: u64 = 20;

    fn fixture_paths() -> Option<(std::path::PathBuf, std::path::PathBuf)> {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
        let ref_path = root.join("parity/realworld/assets/hs37d5.simple.fa");
        let bam = root.join("parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam");
        if !ref_path.is_file() || !bam.is_file() {
            return None;
        }
        Some((ref_path, bam))
    }

    fn digest(bytes: &[u8]) -> u64 {
        let mut h = DefaultHasher::new();
        bytes.hash(&mut h);
        h.finish()
    }

    fn cigar_str(h: &Haplotype) -> String {
        h.cigar
            .as_ref()
            .map(|c| c.to_gatk_string())
            .unwrap_or_else(|| "NA".into())
    }

    fn loc_str(h: &Haplotype) -> String {
        h.genome_loc
            .map(|g| format!("{}-{}", g.start_1based(), g.end_1based()))
            .unwrap_or_else(|| "None".into())
    }

    fn event_has(
        events: &[crate::event_map::VariationEvent],
        start: u64,
        r: &str,
        a: &str,
    ) -> bool {
        events
            .iter()
            .any(|e| e.start_1based.get() == start && e.ref_allele == r && e.alt_allele == a)
    }

    fn trim_fail_reason(h: &Haplotype, span: &GenomeLoc) -> &'static str {
        let Some(genome) = h.genome_loc else {
            return "GENOME_LOC_NONE";
        };
        if !genome.contains(span) {
            return "GENOME_LOC_DOES_NOT_CONTAIN_SPAN";
        }
        let Some(cigar) = h.cigar.as_ref() else {
            return "CIGAR_NONE";
        };
        let read_len: usize = cigar
            .elements
            .iter()
            .map(|e| {
                if e.operator.consumes_read_bases() {
                    e.length
                } else {
                    0
                }
            })
            .sum();
        if h.bases.len() != read_len {
            return "CIGAR_READ_LEN_NE_BASES";
        }
        let new_start = span.start_1based().saturating_sub(genome.start_1based()) as usize;
        let new_stop = new_start + span.reference_span_length().saturating_sub(1) as usize;
        match get_bases_covering_ref_interval(new_start, new_stop, &h.bases, 0, cigar) {
            None => "GET_BASES_COVERING_NONE",
            Some(b) if b.is_empty() => "GET_BASES_COVERING_EMPTY",
            Some(_) => {
                if h.trim(span, false).is_none() {
                    "TRIM_LEADING_TRAILING_INSERTION_OR_EMPTY_CIGAR"
                } else {
                    "OK"
                }
            }
        }
    }

    fn java_style_padded_span(
        events: &[crate::event_map::VariationEvent],
        region: &crate::assembly_region_iterator::AssemblyRegion,
    ) -> (u64, u64, u64, u64) {
        let overlapping: Vec<_> = events
            .iter()
            .filter(|e| {
                e.start_1based.get() <= region.end.get() && e.end_1based.get() >= region.start.get()
            })
            .collect();
        let mut min_start = overlapping
            .iter()
            .map(|e| e.start_1based.get())
            .min()
            .unwrap();
        let mut max_end = overlapping
            .iter()
            .map(|e| e.end_1based.get())
            .max()
            .unwrap();
        let variant_start = region.start.get().max(min_start);
        let variant_end = region.end.get().min(max_end);
        for e in &overlapping {
            let padding = if e.is_indel() { 75 } else { SNP_PAD };
            min_start = min_start.min(e.start_1based.get().saturating_sub(padding).max(1));
            max_end = max_end.max(e.end_1based.get().saturating_add(padding));
        }
        let padded_start = region.extended_start.get().max(min_start);
        let padded_end = region.extended_end.get().min(max_end);
        (variant_start, variant_end, padded_start, padded_end)
    }

    fn dump_hap_stage(label: &str, haps: &[Haplotype], ref_bytes: &[u8]) {
        let n_alt = haps.iter().filter(|h| !h.is_reference).count();
        let n_base_diff = haps
            .iter()
            .filter(|h| h.bases.as_slice() != ref_bytes)
            .count();
        eprintln!(
            "{label} n_haps={} n_alt={} n_bases_ne_ref={} ref_len={}",
            haps.len(),
            n_alt,
            n_base_diff,
            ref_bytes.len()
        );
        for (i, h) in haps.iter().enumerate() {
            eprintln!(
                "  HAP[{i}] ref={} cigar={} len={} align={} loc={} digest={:016x} bases_eq_ref={}",
                h.is_reference,
                cigar_str(h),
                h.bases.len(),
                h.alignment_start_hap_wrt_ref,
                loc_str(h),
                digest(&h.bases),
                h.bases.as_slice() == ref_bytes
            );
        }
    }

    #[test]
    fn six_r27_call_region_none_after_java_eventmap() {
        assert!(
            !AssembleReadsArgs::default()
                .assembler
                .allow_non_unique_kmers_in_ref,
            "allowNonUniqueKmersInRef must remain false"
        );
        let assemble_src = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/src/assembly_based_caller.rs"
        ));
        assert!(
            !assemble_src.contains("92317399"),
            "no mid-B coordinate special case in assembleReads"
        );

        let Some((ref_fasta, bam_path)) = fixture_paths() else {
            eprintln!("Real-data mid-B comparison unavailable");
            return;
        };

        let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
        let walk_iv = parse_intervals_cli_string(&dict, "2:92317000-92319000").expect("iv");
        let filters = crate::read_model::ReadFilterParams::gatk_standard_hc();
        let cfg =
            crate::walker_traversal::WalkerTraversalConfig::gatk_haplotype_caller_production(100);
        let walk = crate::walker_traversal::traverse_assembly_region_walker(
            &dict, &walk_iv, &ref_fasta, &bam_path, &filters, &cfg,
        )
        .expect("walk");
        let regions = crate::walker_traversal::flatten_assembly_regions(&walk);
        let region = regions
            .iter()
            .find(|r| {
                matches!(
                    call_disposition(r),
                    AssemblyRegionCallDisposition::ActiveFull
                ) && r.start.get() <= SITE_CA
                    && r.end.get() >= SITE_CA
            })
            .expect("ActiveFull mid-B")
            .clone();

        assert_eq!(region.start.get(), JAVA_ACTIVE.0);
        assert_eq!(region.end.get(), JAVA_ACTIVE.1);
        assert_eq!(region.extended_start.get(), JAVA_EXTENDED.0);
        assert_eq!(region.extended_end.get(), JAVA_EXTENDED.1);
        assert_eq!(region.reads.len(), 2);
        assert!(
            !crate::read_threading_assembler::region_overlaps_p12_cluster(
                region.start.get(),
                region.end.get(),
            ),
            "canonical mid-B must not overlap the P12 cluster"
        );
        assert!(
            region.end.get() < P12_CLUSTER_TTC_START
                || region.start.get() > P12_CLUSTER_TTC_START.saturating_add(3),
            "mid-B is not the TTC/ATG cluster"
        );

        eprintln!("PINNED_JAVA GATK 4.4.0.0 SHA=2dbc025821bc5f686c423ff332a41e6cef892a77");
        eprintln!(
            "REGION active={}..{} extended={}..{} n_reads={}",
            region.start.get(),
            region.end.get(),
            region.extended_start.get(),
            region.extended_end.get(),
            region.reads.len()
        );

        let args = CallRegionArgs::strict_java();
        eprintln!(
            "CALL_REGION_FLAGS is_strict_java={} is_java_compatible={} disable_optimizations={} enable_read_event_supplement={}",
            args.is_strict_java(),
            args.is_java_compatible(),
            args.disable_optimizations,
            args.enable_read_event_supplement
        );
        let outcome =
            HaplotypeCallerEngine::call_region(&region, &dict, &ref_fasta, &args).expect("call");
        let audit = take_call_region_audit();
        let none_reason = audit
            .none_at
            .as_deref()
            .unwrap_or("CALL_REGION_RETURNED_SOME");
        eprintln!(
            "CALL_REGION={}",
            if outcome.is_some() { "Some" } else { "None" }
        );
        eprintln!("NONE_REASON={none_reason}");
        eprintln!("NONE_AT_ENGINE_BRANCH={:?}", audit.none_at);
        eprintln!(
            "PRE_TRIM_AUDIT n_haps={} n_alt={} trim_variation_present={} overlapping={} variant={:?}..{:?} padded={:?}..{:?}",
            audit.n_haps,
            audit.n_alt_haps,
            audit.trim_variation_present,
            audit.n_trim_overlapping,
            audit.trim_variant_start,
            audit.trim_variant_end,
            audit.trim_padded_start,
            audit.trim_padded_end
        );
        for (i, h) in audit.hap_cigars.iter().enumerate() {
            eprintln!("  PRE_TRIM_HAP[{i}] {h}");
        }
        eprintln!("ENGINE_POST_TRIM_STAGES n={}", audit.post_trim_stages.len());
        eprintln!("NEEDS_POST_TRIM_RESYNC={:?}", audit.needs_post_trim_resync);
        for s in &audit.post_trim_stages {
            eprintln!("  STAGE {s}");
        }
        eprintln!("PRE_TRIM_EVENTMAP n={}", audit.events.len());
        for e in &audit.events {
            eprintln!(
                "  EVENT {}-{} {}/{} indel={}",
                e.start, e.end, e.ref_al, e.alt_al, e.is_indel
            );
        }
        eprintln!(
            "AFTER_READ_FILTER n_haps={:?} n_events={:?} has_variation_for_calling={:?} n_reads={:?}",
            audit.n_haps_after_trim,
            audit.n_events_after_trim,
            audit.has_variation_for_calling,
            audit.n_reads_after_filter
        );
        eprintln!("AFTER_TRIM_EVENTMAP n={}", audit.events_after_trim.len());
        for e in &audit.events_after_trim {
            eprintln!(
                "  EVENT {}-{} {}/{} indel={}",
                e.start, e.end, e.ref_al, e.alt_al, e.is_indel
            );
        }

        let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
        let mut owned = region.clone();
        let mut assemble_args = args.assemble.clone();
        assemble_args.strict_java_assembly = true;
        let padded = assembly_reference_read(&dict, &mut ref_cache, &region).expect("pad");
        let (pad_start, pad_end) = padded_reference_loc(&region, &dict);
        let graph_ref = create_graph_reference_read(&padded, &region, &dict);
        eprintln!(
            "REF_WINDOWS graph_len={} graph_digest={:016x} pad={pad_start}-{pad_end} pad_len={}",
            graph_ref.bases.len(),
            digest(&graph_ref.bases),
            padded.bases.len()
        );

        let raw = assemble_from_ref_and_reads(
            &graph_ref,
            &crate::assembly_region_finalize::records_to_assembly_reads(
                &crate::assembly_region_finalize::finalize_region_reads_for_assembly(
                    &region.reads,
                    &region,
                    assemble_args.correct_overlapping_base_qualities,
                    crate::assembly_region_finalize::gatk_min_tail_quality_for_assembly(
                        assemble_args.assembler.min_base_quality,
                    ),
                    false,
                ),
            ),
            &assemble_args.assembler,
        )
        .expect("raw assemble");
        dump_hap_stage(
            "STAGE_RAW_SEQGRAPH_KBEST",
            &raw.haplotypes,
            graph_ref.bases.as_slice(),
        );

        let assembled =
            assemble_reads_with_finalized(&mut owned, &dict, &mut ref_cache, &assemble_args)
                .expect("production assemble");
        let mut untrimmed = assembled.assembly;
        dump_hap_stage(
            "STAGE_AFTER_NORMALIZE_AND_EVENTMAP",
            &untrimmed.haplotypes,
            untrimmed.reference_bases(),
        );
        eprintln!(
            "STAGE_AFTER_NORMALIZE n_events={} has_variation_for_calling={} is_variation_present={} variation_present={} pad_start={} ref_len={}",
            untrimmed.variation_events.len(),
            untrimmed.has_variation_for_calling(),
            untrimmed.is_variation_present(),
            untrimmed.variation_present,
            untrimmed.padded_reference_start_1based(),
            untrimmed.reference_bases().len()
        );
        for e in untrimmed.variation_events() {
            eprintln!(
                "  EVENTMAP {}-{} {}/{}",
                e.start_1based.get(),
                e.end_1based.get(),
                e.ref_allele,
                e.alt_allele
            );
        }
        assert!(
            event_has(untrimmed.variation_events(), SITE_CA, "C", "A"),
            "untrimmed EventMap must contain 92317399 C/A"
        );
        assert!(
            event_has(untrimmed.variation_events(), SITE_TC, "T", "C"),
            "untrimmed EventMap must contain 92317407 T/C"
        );
        assert!(
            event_has(untrimmed.variation_events(), SITE_GC, "G", "C"),
            "untrimmed EventMap must contain 92317412 G/C"
        );

        let apply_pad_u = untrimmed
            .haplotypes
            .iter()
            .find(|h| h.is_reference)
            .and_then(|h| h.genome_loc.map(|g| g.start_1based()))
            .unwrap_or_else(|| untrimmed.padded_reference_start_1based());
        eprintln!(
            "APPLY_PAD_U={apply_pad_u} padded_reference_start={} extended_start={}",
            untrimmed.padded_reference_start_1based(),
            region.extended_start.get()
        );
        Haplotype::tag_padded_reference_span(&mut untrimmed.haplotypes, apply_pad_u);
        dump_hap_stage(
            "STAGE_AFTER_TAG_PADDED_SPAN",
            &untrimmed.haplotypes,
            untrimmed.reference_bases(),
        );

        let trim_variants: Vec<TrimVariant> = untrimmed
            .variation_events()
            .iter()
            .map(|e| TrimVariant {
                contig: e.contig.clone(),
                start: e.start_1based.get(),
                end: e.end_1based.get(),
                is_indel: e.is_indel(),
            })
            .collect();
        let trimmer =
            AssemblyRegionTrimmer::new(AssemblyRegionTrimmerConfig::gatk_defaults(), &dict, "2");
        let trim_result = trimmer.trim(&region, &trim_variants, Some(&region.reference));
        eprintln!(
            "TRIM_RESULT variation_present={} variant={:?}..{:?} padded={:?}..{:?}",
            trim_result.variation_present,
            trim_result.variant_start,
            trim_result.variant_end,
            trim_result.padded_variant_start,
            trim_result.padded_variant_end
        );
        let (jv_s, jv_e, jp_s, jp_e) =
            java_style_padded_span(untrimmed.variation_events(), &region);
        eprintln!(
            "JAVA_STYLE_TRIM variant={jv_s}-{jv_e} padded={jp_s}-{jp_e} rust_padded={:?}-{:?}",
            trim_result.padded_variant_start, trim_result.padded_variant_end
        );
        eprintln!(
            "TRIM_PAD_ACCUMULATION_DIVERGES={}",
            trim_result.padded_variant_end != Some(jp_e)
                || trim_result.padded_variant_start != Some(jp_s)
        );

        let oracle_only: Vec<TrimVariant> = untrimmed
            .variation_events()
            .iter()
            .filter(|e| {
                let s = e.start_1based.get();
                s == SITE_CA || s == SITE_TC || s == SITE_GC
            })
            .map(|e| TrimVariant {
                contig: e.contig.clone(),
                start: e.start_1based.get(),
                end: e.end_1based.get(),
                is_indel: e.is_indel(),
            })
            .collect();
        let oracle_trim = trimmer.trim(&region, &oracle_only, Some(&region.reference));
        eprintln!(
            "ORACLE_ONLY_TRIM variation_present={} variant={:?}..{:?} padded={:?}..{:?}",
            oracle_trim.variation_present,
            oracle_trim.variant_start,
            oracle_trim.variant_end,
            oracle_trim.padded_variant_start,
            oracle_trim.padded_variant_end
        );
        eprintln!(
            "EXTRA_SNPS_CHANGE_TRIM_PRESENT={}",
            oracle_trim.variation_present != trim_result.variation_present
        );

        let rust_span = GenomeLoc::new(
            trim_result.padded_variant_start.expect("padded start"),
            trim_result.padded_variant_end.expect("padded end"),
        );
        let java_span = GenomeLoc::new(jp_s, jp_e);
        let ref_bytes = untrimmed.reference_bases().to_vec();
        eprintln!(
            "PER_HAP_TRIM rust_span={}..{}",
            rust_span.start_1based(),
            rust_span.end_1based()
        );
        let mut first_drop: Option<&'static str> = None;
        let mut n_trim_ok_false = 0usize;
        let mut n_trim_ok_true = 0usize;
        for (i, h) in untrimmed.haplotypes.iter().enumerate() {
            let reason = trim_fail_reason(h, &rust_span);
            let t_false = h.trim(&rust_span, false);
            let t_true = h.trim(&rust_span, true);
            if t_false.is_some() {
                n_trim_ok_false += 1;
            }
            if t_true.is_some() {
                n_trim_ok_true += 1;
            }
            let trimmed_eq_ref = t_false
                .as_ref()
                .map(|t| {
                    let off = rust_span
                        .start_1based()
                        .saturating_sub(untrimmed.padded_reference_start_1based())
                        as usize;
                    let len = t.bases.len();
                    ref_bytes.get(off..off.saturating_add(len)) == Some(t.bases.as_slice())
                })
                .unwrap_or(false);
            eprintln!(
                "  TRIM_HAP[{i}] ref={} cigar={} loc={} contains={} reason={reason} trim(false)={} trim(true)={} trimmed_eq_sliced_ref={trimmed_eq_ref} bases_ne_full_ref={}",
                h.is_reference,
                cigar_str(h),
                loc_str(h),
                h.genome_loc
                    .map(|g| g.contains(&rust_span))
                    .unwrap_or(false),
                t_false.is_some(),
                t_true.is_some(),
                h.bases.as_slice() != ref_bytes.as_slice()
            );
            if !h.is_reference && t_false.is_none() && first_drop.is_none() {
                first_drop = Some("Haplotype::trim returned None");
            }
            let mut retagged = h.clone();
            let ext = GenomeLoc::new(JAVA_EXTENDED.0, JAVA_EXTENDED.1);
            retagged.genome_loc = Some(ext);
            let retag_ok = retagged.trim(&rust_span, false).is_some();
            eprintln!(
                "    RETAG_EXTENDED_LOC trim(false)={retag_ok} reason={}",
                trim_fail_reason(&retagged, &rust_span)
            );
            let java_ok = h.trim(&java_span, true);
            eprintln!(
                "    JAVA_SPAN_ignoreRefState=true trim={} reason={}",
                java_ok.is_some(),
                trim_fail_reason(h, &java_span)
            );
        }
        eprintln!(
            "TRIM_SUCCESS n_ok_ignoreRefState_false={n_trim_ok_false} n_ok_true={n_trim_ok_true} / {}",
            untrimmed.haplotypes.len()
        );

        let trimmed_region = AssemblyRegionTrimmer::apply_trim(&region, &trim_result);
        eprintln!(
            "APPLY_TRIM active={}..{} extended={}..{} n_reads={}",
            trimmed_region.start.get(),
            trimmed_region.end.get(),
            trimmed_region.extended_start.get(),
            trimmed_region.extended_end.get(),
            trimmed_region.reads.len()
        );
        let after_trim_to = untrimmed.trim_to(&trimmed_region).expect("trim_to");
        dump_hap_stage(
            "STAGE_AFTER_TRIM_TO",
            &after_trim_to.haplotypes,
            after_trim_to.reference_bases(),
        );
        eprintln!(
            "STAGE_AFTER_TRIM_TO n_events={} has_variation_for_calling={} is_variation_present={} variation_present={}",
            after_trim_to.variation_events.len(),
            after_trim_to.has_variation_for_calling(),
            after_trim_to.is_variation_present(),
            after_trim_to.variation_present
        );
        for e in after_trim_to.variation_events() {
            eprintln!(
                "  POST_TRIM_EVENT {}-{} {}/{}",
                e.start_1based.get(),
                e.end_1based.get(),
                e.ref_allele,
                e.alt_allele
            );
        }
        if after_trim_to.haplotypes.iter().all(|h| h.is_reference)
            && untrimmed.haplotypes.iter().any(|h| !h.is_reference)
        {
            first_drop =
                Some(first_drop.unwrap_or("AssemblyResultSet::trim_to dropped all ALT haplotypes"));
        }
        eprintln!(
            "FIRST_ALT_DROP={}",
            first_drop.unwrap_or("ALTs survived trim_to")
        );

        let pad_start = untrimmed.padded_reference_start_1based();
        let ref_hap_t = after_trim_to
            .haplotypes
            .iter()
            .find(|h| h.is_reference)
            .expect("trimmed REF");
        for (site, r, a) in [
            (SITE_CT, b'C', b'T'),
            (SITE_CG, b'C', b'G'),
            (SITE_CA, b'C', b'A'),
            (SITE_TC, b'T', b'C'),
            (SITE_GC, b'G', b'C'),
        ] {
            let mut n_exact = 0usize;
            for (i, h) in after_trim_to.haplotypes.iter().enumerate() {
                if h.is_reference {
                    continue;
                }
                let b = crate::hc_allele_mapping::hap_base_at_ref_locus(h, pad_start, site);
                let apply_pad = ref_hap_t
                    .genome_loc
                    .map(|g| g.start_1based())
                    .unwrap_or(pad_start);
                let apply_off = site.saturating_sub(apply_pad) as usize;
                let slice = h.bases.get(apply_off).copied();
                let hit = b == Some(a) || slice == Some(a);
                if hit {
                    n_exact += 1;
                }
                eprintln!(
                    "  SUPPORT_HAP[{i}] {site} hap_base={:?} slice={:?} want={} hit={hit}",
                    b.map(|c| c as char),
                    slice.map(|c| c as char),
                    a as char
                );
                let _ = r;
            }
            eprintln!("SUPPORT {site} n_alt_exact={n_exact}");
        }
        for e in after_trim_to.variation_events() {
            if e.start_1based.get() > SITE_GC {
                break;
            }
            eprintln!(
                "EMIT_CANDIDATE {} {}/{} = {}",
                e.start_1based.get(),
                e.ref_allele,
                e.alt_allele,
                crate::java_hc_site_semantics::is_strict_java_production_emit_candidate(e)
            );
        }
        let mut filtered = after_trim_to.clone();
        let _ = crate::allele_filtering::filter_assembly_and_likelihoods(
            &mut filtered,
            Vec::new(),
            crate::allele_filter_options::AlleleFilterOptions::from_strict_java(
                true,
                Some(region.start.get()),
                Some(region.end.get()),
            ),
        )
        .expect("allele filter");
        dump_hap_stage(
            "STAGE_AFTER_ALLELE_FILTER_EMPTY_LL",
            &filtered.haplotypes,
            filtered.reference_bases(),
        );
        eprintln!(
            "STAGE_AFTER_ALLELE_FILTER_EMPTY_LL n_events={} has_variation_for_calling={} is_variation_present={}",
            filtered.variation_events.len(),
            filtered.has_variation_for_calling(),
            filtered.is_variation_present()
        );
        eprintln!(
            "HAS_C_A untrimmed={} after_trim_to={}",
            event_has(untrimmed.variation_events(), SITE_CA, "C", "A"),
            event_has(after_trim_to.variation_events(), SITE_CA, "C", "A")
        );
        eprintln!(
            "EXTRA_SNPS_PRESENT 92317361 C/T and 92317371 C/G = {}",
            event_has(untrimmed.variation_events(), SITE_CT, "C", "T")
                && event_has(untrimmed.variation_events(), SITE_CG, "C", "G")
        );

        if outcome.is_none() {
            assert!(
                audit.none_at.as_deref().is_some_and(|s| {
                    s.contains("has_variation_for_calling")
                        || s.contains("trim_result.variation_present")
                }),
                "Ok(None) must name the engine branch, got {:?}",
                audit.none_at
            );
            assert!(
                audit.trim_variation_present,
                "trimmer variation_present must be true so the None is post-trim"
            );
            assert_eq!(
                audit.has_variation_for_calling,
                Some(false),
                "the recorded None path is has_variation_for_calling=false"
            );
            assert!(
                audit
                    .post_trim_stages
                    .iter()
                    .any(|s| s.contains("after_trim_to") && s.contains("n_alt=3")),
                "ALTs must still be present immediately after trim_to; got {:?}",
                audit.post_trim_stages
            );
            assert!(
                audit
                    .post_trim_stages
                    .iter()
                    .any(|s| s.contains("after_early_allele_filter") && s.contains("n_alt=0")),
                "early allele filter must drop all ALTs; got {:?}",
                audit.post_trim_stages
            );
        }
        let _ = pad_end;
        let _ = n_trim_ok_true;
    }
}
