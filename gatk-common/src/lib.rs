//! GATK-RS Common Library
//! This crate contains shared utilities and error handling for GATK-RS.

pub mod config;
pub mod error;
pub mod logging;

// Re-export commonly used items
pub use config::*;
pub use error::*;
pub use logging::*;

#[cfg(test)]
mod tests {

    #[test]
    fn test_basic_functionality() {
        // Basic test to ensure the library compiles and links correctly
        let sum = 1 + 1;
        assert_eq!(sum, 2);
    }
}
