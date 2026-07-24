//! production `strict_java` EventMap path (CIGAR-only, no P12 bridges by default).

use gatk_haplotypecaller::{
    strict_java_asm8_only_enabled, strict_java_p12_ensure_bridges_enabled, CallRegionArgs,
};

#[test]
fn strict_java_assembly_uses_dangling_java_exact() {
    let args = CallRegionArgs::strict_java();
    assert!(args.assemble.strict_java_assembly);
    assert!(args.assemble.assembler.dangling_java_exact);
}

#[test]
fn production_strict_default_disables_p12_bridges() {
    if std::env::var("GATK_RS_P12_ENSURE_BRIDGES").is_err()
        && std::env::var("GATK_RS_ASM8_ONLY").is_err()
    {
        assert!(!strict_java_p12_ensure_bridges_enabled());
        assert!(strict_java_asm8_only_enabled());
    }
}

#[test]
fn graph_only_production_implies_bridges_off() {
    if strict_java_asm8_only_enabled() {
        assert!(!strict_java_p12_ensure_bridges_enabled());
    }
}
