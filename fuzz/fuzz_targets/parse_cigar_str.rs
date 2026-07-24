#![no_main]
use libfuzzer_sys::fuzz_target;
use gatk_core::parse_cigar_str;

fuzz_target!(|data: &[u8]| {
    let Ok(s) = std::str::from_utf8(data) else {
        return;
    };
    let _ = parse_cigar_str(s);
});
