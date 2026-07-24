use gatk_core::prop_test;

prop_test! {
    ascii_uppercase_is_still_ascii in |value: u8| {
        let c = (value % 26 + b'a') as char;
        assert!(c.to_ascii_uppercase().is_ascii_uppercase());
    }
}
