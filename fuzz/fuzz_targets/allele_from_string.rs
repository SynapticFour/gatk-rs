#![no_main]
use libfuzzer_sys::fuzz_target;
use gatk_core::Allele;

fuzz_target!(|data: &[u8]| {
    let Ok(s) = std::str::from_utf8(data) else {
        return;
    };
    let _ = Allele::from_string(s);
});
