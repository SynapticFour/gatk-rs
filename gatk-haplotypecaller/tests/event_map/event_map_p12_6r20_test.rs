//! 6R.20 TEST-ONLY: diagnose `call_region(None)` at mid-B 92317399 after Java EventMap
//! geometry. Does not change production algorithms.

#[cfg(test)]
mod traces {
    use crate::assembly_region_trimmer::{
        AssemblyRegionTrimmer, AssemblyRegionTrimmerConfig, TrimVariant,
    };
    use crate::engine::{
        take_call_region_audit, AuditEvent, CallRegionArgs, HaplotypeCallerEngine,
    };
    use crate::walker_apply::{call_disposition, AssemblyRegionCallDisposition};
    use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
    use std::path::Path;

    const SITE: u64 = 92_317_399;
    const EXPECT_REF: &str = "C";
    const EXPECT_ALT: &str = "A";

    fn fixture_paths() -> Option<(std::path::PathBuf, std::path::PathBuf)> {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
        let ref_path = root.join("parity/realworld/assets/hs37d5.simple.fa");
        let bam = root.join("parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam");
        if !ref_path.is_file() || !bam.is_file() {
            return None;
        }
        Some((ref_path, bam))
    }

    fn old_max_span_end(start: u64, ref_len: usize, alt_len: usize) -> u64 {
        start.saturating_add(ref_len.max(alt_len).saturating_sub(1) as u64)
    }

    fn java_ref_span_end(start: u64, ref_len: usize) -> u64 {
        start.saturating_add(ref_len.saturating_sub(1) as u64)
    }

    fn dump_events(label: &str, events: &[AuditEvent]) {
        eprintln!("EVENTS {label} n={}", events.len());
        for e in events {
            let kind = if e.is_indel {
                if e.ref_al.len() < e.alt_al.len() {
                    "INS"
                } else {
                    "DEL"
                }
            } else {
                "SNP"
            };
            eprintln!("  {}-{} {}→{} {kind}", e.start, e.end, e.ref_al, e.alt_al);
        }
    }

    fn events_near<'a>(events: &'a [AuditEvent], loc: u64, pad: u64) -> Vec<&'a AuditEvent> {
        events
            .iter()
            .filter(|e| e.start <= loc.saturating_add(pad) && e.end + pad >= loc)
            .collect()
    }

    fn trim_vars_from_events(
        events: &[AuditEvent],
        contig: &str,
        end_of: impl Fn(&AuditEvent) -> u64,
    ) -> Vec<TrimVariant> {
        events
            .iter()
            .map(|e| TrimVariant {
                contig: contig.to_string(),
                start: e.start,
                end: end_of(e),
                is_indel: e.is_indel,
            })
            .collect()
    }

    fn rewrite_cat_tat(events: &[AuditEvent], alt: &str) -> Vec<AuditEvent> {
        events
            .iter()
            .map(|e| {
                let mut o = e.clone();
                if e.ref_al == "C" && (e.alt_al == "CAT" || e.alt_al == "TAT") {
                    o.alt_al = alt.to_string();
                    o.is_indel = o.ref_al.len() != o.alt_al.len();
                }
                o
            })
            .collect()
    }

    #[test]
    fn six_r20_mid_b_java_eventmap_geometry() {
        let Some((ref_fasta, bam)) = fixture_paths() else {
            eprintln!("Real-data mid-B comparison unavailable");
            return;
        };
        let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
        let specs = parse_intervals_cli_string(&dict, "2:92317000-92319000").expect("interval");
        let walk = crate::walker_traversal::traverse_assembly_region_walker(
            &dict,
            &specs,
            &ref_fasta,
            &bam,
            &crate::read_model::ReadFilterParams::gatk_standard_hc(),
            &crate::walker_traversal::WalkerTraversalConfig::gatk_haplotype_caller_production(100),
        )
        .expect("walk");
        let regions = crate::walker_traversal::flatten_assembly_regions(&walk);
        let region = regions
            .iter()
            .find(|r| {
                matches!(
                    call_disposition(r),
                    AssemblyRegionCallDisposition::ActiveFull
                ) && r.start.get() <= SITE
                    && r.end.get() >= SITE
            })
            .expect("ActiveFull region covering 92317399");
        eprintln!(
            "REGION active={}..{} extended={}..{} n_reads={}",
            region.start.get(),
            region.end.get(),
            region.extended_start.get(),
            region.extended_end.get(),
            region.reads.len()
        );

        let args = CallRegionArgs::strict_java();
        let outcome =
            HaplotypeCallerEngine::call_region(region, &dict, &ref_fasta, &args).expect("call");
        let audit = take_call_region_audit();
        let java_call = if outcome.is_some() { "Some" } else { "None" };
        eprintln!("JAVA_44_CALL_REGION={java_call}");
        eprintln!("NONE_AT={:?}", audit.none_at);
        eprintln!(
            "TRIM present={} overlapping={} variant={:?}..{:?} padded={:?}..{:?}",
            audit.trim_variation_present,
            audit.n_trim_overlapping,
            audit.trim_variant_start,
            audit.trim_variant_end,
            audit.trim_padded_start,
            audit.trim_padded_end
        );
        eprintln!(
            "RESCUE cluster_reads_support={} read_variation_in_active={} disable_opt={} n_haps={} n_alt={}",
            audit.cluster_reads_support,
            audit.read_variation_in_active,
            audit.disable_optimizations,
            audit.n_haps,
            audit.n_alt_haps
        );
        for (i, h) in audit.hap_cigars.iter().enumerate() {
            eprintln!("  HAP[{i}] {h}");
        }
        eprintln!(
            "AFTER_TRIM n_haps={:?} n_events={:?} has_variation={:?} n_reads_filter={:?}",
            audit.n_haps_after_trim,
            audit.n_events_after_trim,
            audit.has_variation_for_calling,
            audit.n_reads_after_filter
        );
        dump_events("untrimmed_java44_eventmap", &audit.events);
        dump_events("after_trim_eventmap", &audit.events_after_trim);
        let near = events_near(&audit.events, SITE, 50);
        eprintln!("NEAR_{SITE} n={}", near.len());
        for e in &near {
            eprintln!(
                "  {}-{} {}→{} indel={}",
                e.start, e.end, e.ref_al, e.alt_al, e.is_indel
            );
        }
        let has_ca_untrimmed = audit
            .events
            .iter()
            .any(|e| e.start == SITE && e.ref_al == EXPECT_REF && e.alt_al == EXPECT_ALT);
        let has_ca_after = audit
            .events_after_trim
            .iter()
            .any(|e| e.start == SITE && e.ref_al == EXPECT_REF && e.alt_al == EXPECT_ALT);
        eprintln!("HAS_C_A untrimmed={has_ca_untrimmed} after_trim={has_ca_after}");

        eprintln!("TRIM_VARIANTS n={}", audit.trim_variants.len());
        for v in &audit.trim_variants {
            eprintln!(
                "  {}-{} indel={} overlaps_active={}",
                v.start, v.end, v.is_indel, v.overlaps_active
            );
        }

        let java_ends: Vec<(u64, u64, &str, &str)> = audit
            .events
            .iter()
            .map(|e| (e.start, e.end, e.ref_al.as_str(), e.alt_al.as_str()))
            .collect();
        let old_ends: Vec<(u64, u64, String, String)> = audit
            .events
            .iter()
            .map(|e| {
                (
                    e.start,
                    old_max_span_end(e.start, e.ref_al.len(), e.alt_al.len()),
                    e.ref_al.clone(),
                    e.alt_al.clone(),
                )
            })
            .collect();
        let end_geometry_differs = java_ends
            .iter()
            .zip(old_ends.iter())
            .any(|(j, o)| j.0 == o.0 && j.1 != o.1);
        eprintln!("END_GEOMETRY_DIFFERS={end_geometry_differs}");
        for (j, o) in java_ends.iter().zip(old_ends.iter()) {
            if j.1 != o.1 {
                eprintln!(
                    "  END_DIFF {} {}→{} java_end={} old_end={}",
                    j.0, j.2, j.3, j.1, o.1
                );
            }
        }

        let trimmer =
            AssemblyRegionTrimmer::new(AssemblyRegionTrimmerConfig::gatk_defaults(), &dict, "2");
        let java_tv = trim_vars_from_events(&audit.events, "2", |e| e.end);
        let old_tv = trim_vars_from_events(&audit.events, "2", |e| {
            old_max_span_end(e.start, e.ref_al.len(), e.alt_al.len())
        });
        let java_trim = trimmer.trim(region, &java_tv, Some(&region.reference));
        let old_trim = trimmer.trim(region, &old_tv, Some(&region.reference));
        eprintln!(
            "CAUSAL_END java_trim.variation_present={} old_trim.variation_present={}",
            java_trim.variation_present, old_trim.variation_present
        );
        eprintln!(
            "CAUSAL_END java_span={:?}..{:?} old_span={:?}..{:?}",
            java_trim.variant_start,
            java_trim.variant_end,
            old_trim.variant_start,
            old_trim.variant_end
        );

        let has_cat = audit
            .events
            .iter()
            .any(|e| e.ref_al == "C" && e.alt_al == "CAT");
        let has_tat = audit
            .events
            .iter()
            .any(|e| e.ref_al == "C" && e.alt_al == "TAT");
        eprintln!("HAS_C_CAT={has_cat} HAS_C_TAT={has_tat}");
        let as_cat = rewrite_cat_tat(&audit.events, "CAT");
        let as_tat = rewrite_cat_tat(&audit.events, "TAT");
        let cat_tv =
            trim_vars_from_events(&as_cat, "2", |e| java_ref_span_end(e.start, e.ref_al.len()));
        let tat_tv =
            trim_vars_from_events(&as_tat, "2", |e| java_ref_span_end(e.start, e.ref_al.len()));
        let cat_trim = trimmer.trim(region, &cat_tv, Some(&region.reference));
        let tat_trim = trimmer.trim(region, &tat_tv, Some(&region.reference));
        eprintln!(
            "CAUSAL_MAKEBLOCK CAT_trim.variation_present={} TAT_trim.variation_present={}",
            cat_trim.variation_present, tat_trim.variation_present
        );

        let insertion_end_causes_trim_flip =
            java_trim.variation_present != old_trim.variation_present;
        let makeblock_alt_causes_trim_flip =
            (has_cat || has_tat) && cat_trim.variation_present != tat_trim.variation_present;
        eprintln!("INSERTION_END_CAUSES_TRIM_FLIP={insertion_end_causes_trim_flip}");
        eprintln!("MAKEBLOCK_ALT_CAUSES_TRIM_FLIP={makeblock_alt_causes_trim_flip}");
        eprintln!(
            "OLD_CALL_REGION=UNKNOWN (cannot replay production EventMap); OLD_TRIM variation_present={}",
            old_trim.variation_present
        );
        eprintln!("JAVA_CALL_REGION={java_call}");

        if outcome.is_none() {
            assert!(
                audit.none_at.is_some(),
                "Ok(None) must record the exact engine.rs branch"
            );
        }
        let _ = end_geometry_differs;
        let _ = java_ends;
    }
}
