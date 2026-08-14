//! Error handling for GATK-RS.
//! Small thiserror enum (`Io` / `Parse` / `Config` / `Algorithm`) plus named constructors
//! that preserve Display prefixes and CLI classification.

use std::backtrace::Backtrace;
use thiserror::Error;

/// Error severity levels (mirroring GATK's logging levels)
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum ErrorSeverity {
    Warning,
    Error,
    Fatal,
}

/// Error context information for better debugging
#[derive(Debug)]
pub struct ErrorContext {
    /// Tool name where error occurred
    pub tool: Option<String>,
    /// File or operation context
    pub operation: Option<String>,
    /// Genomic coordinates if applicable
    pub position: Option<(String, u64)>, // (chromosome, position)
    /// Stack trace information
    pub backtrace: Option<Backtrace>,
    /// Error severity
    pub severity: ErrorSeverity,
    /// Additional context key-value pairs
    pub context: std::collections::HashMap<String, String>,
}

impl Clone for ErrorContext {
    fn clone(&self) -> Self {
        Self {
            tool: self.tool.clone(),
            operation: self.operation.clone(),
            position: self.position.clone(),
            backtrace: None, // Don't clone backtrace for performance
            severity: self.severity,
            context: self.context.clone(),
        }
    }
}

impl ErrorContext {
    /// Create a new error context
    pub fn new() -> Self {
        Self {
            tool: None,
            operation: None,
            position: None,
            backtrace: if std::env::var_os("RUST_BACKTRACE").is_some() {
                Backtrace::capture().into()
            } else {
                None
            },
            severity: ErrorSeverity::Error,
            context: std::collections::HashMap::new(),
        }
    }

    /// Set the tool context
    pub fn with_tool(mut self, tool: impl Into<String>) -> Self {
        self.tool = Some(tool.into());
        self
    }

    /// Set the operation context
    pub fn with_operation(mut self, operation: impl Into<String>) -> Self {
        self.operation = Some(operation.into());
        self
    }

    /// Set the genomic position context
    pub fn with_position(mut self, position: (String, u64)) -> Self {
        self.position = Some(position);
        self
    }

    /// Set the error severity
    pub fn with_severity(mut self, severity: ErrorSeverity) -> Self {
        self.severity = severity;
        self
    }

    /// Add context key-value pair
    pub fn add_context(mut self, key: impl Into<String>, value: impl Into<String>) -> Self {
        self.context.insert(key.into(), value.into());
        self
    }
}

impl Default for ErrorContext {
    fn default() -> Self {
        Self::new()
    }
}

/// High-level error class for CLI exit codes and user vs internal reporting.
/// [`User`](GatkErrorClass::User) — argument / configuration / validation mistakes (Picard `USER_ERROR` analog)
/// [`Io`](GatkErrorClass::Io) — filesystem / stream failures
/// [`Internal`](GatkErrorClass::Internal) — algorithm / unexpected failures
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GatkErrorClass {
    User,
    Io,
    Internal,
}

fn boxed_ctx(context: ErrorContext) -> Box<ErrorContext> {
    Box::new(context)
}

/// Main error type for GATK-RS.
/// Named constructors (`argument`, `read`, `assembly`, …) keep Display prefixes and
/// [`GatkError::classification`] the same as the former Java-shaped variants.
#[derive(Error, Debug)]
pub enum GatkError {
    /// I/O related errors
    #[error("I/O error: {message}")]
    Io {
        message: String,
        #[source]
        source: std::io::Error,
        context: Box<ErrorContext>,
    },

    /// Parsing / file-format errors
    #[error("Parse error: {message} at line {line}")]
    Parse {
        message: String,
        line: usize,
        context: Box<ErrorContext>,
    },

    /// User-facing argument / configuration / validation / interval / reference errors.
    /// `message` already includes the historical Display prefix (e.g. `Argument error: …`).
    #[error("{message}")]
    Config {
        message: String,
        context: Box<ErrorContext>,
    },

    /// Algorithm / assembly / generic internal failures.
    /// `message` already includes the historical Display prefix.
    #[error("{message}")]
    Algorithm {
        message: String,
        context: Box<ErrorContext>,
    },
}

impl GatkError {
    fn config_msg(prefix: &str, message: impl Into<String>, context: ErrorContext) -> Self {
        Self::Config {
            message: format!("{prefix}: {}", message.into()),
            context: boxed_ctx(context),
        }
    }

    fn algorithm_msg(prefix: &str, message: impl Into<String>, context: ErrorContext) -> Self {
        Self::Algorithm {
            message: format!("{prefix}: {}", message.into()),
            context: boxed_ctx(context),
        }
    }

    /// Create a file format error
    pub fn file_format<S: Into<String>>(message: S) -> Self {
        Self::algorithm_msg("File format error", message, ErrorContext::new())
    }

    /// Create a file format error with context
    pub fn file_format_with_context<S: Into<String>>(message: S, context: ErrorContext) -> Self {
        Self::algorithm_msg("File format error", message, context)
    }

    /// Create a parse error
    pub fn parse<S: Into<String>>(message: S, line: usize) -> Self {
        Self::Parse {
            message: message.into(),
            line,
            context: boxed_ctx(ErrorContext::new()),
        }
    }

    /// Create a parse error with context
    pub fn parse_with_context<S: Into<String>>(
        message: S,
        line: usize,
        context: ErrorContext,
    ) -> Self {
        Self::Parse {
            message: message.into(),
            line,
            context: boxed_ctx(context),
        }
    }

    /// Create a validation error
    pub fn validation<S: Into<String>>(message: S) -> Self {
        Self::config_msg("Validation error", message, ErrorContext::new())
    }

    /// Create a validation error with context
    pub fn validation_with_context<S: Into<String>>(message: S, context: ErrorContext) -> Self {
        Self::config_msg("Validation error", message, context)
    }

    /// Create an algorithm error
    pub fn algorithm<S: Into<String>>(message: S) -> Self {
        Self::algorithm_msg("Algorithm error", message, ErrorContext::new())
    }

    /// Create an algorithm error with context
    pub fn algorithm_with_context<S: Into<String>>(message: S, context: ErrorContext) -> Self {
        Self::algorithm_msg("Algorithm error", message, context)
    }

    /// Create a memory error
    pub fn memory<S: Into<String>>(message: S) -> Self {
        Self::algorithm_msg("Memory error", message, ErrorContext::new())
    }

    /// Create a memory error with context
    pub fn memory_with_context<S: Into<String>>(message: S, context: ErrorContext) -> Self {
        Self::algorithm_msg("Memory error", message, context)
    }

    /// Create a configuration error
    pub fn configuration<S: Into<String>>(message: S) -> Self {
        Self::config_msg("Configuration error", message, ErrorContext::new())
    }

    /// Create a configuration error with context
    pub fn configuration_with_context<S: Into<String>>(message: S, context: ErrorContext) -> Self {
        Self::config_msg("Configuration error", message, context)
    }

    /// Create an argument error
    pub fn argument<S: Into<String>>(message: S) -> Self {
        Self::config_msg("Argument error", message, ErrorContext::new())
    }

    /// Create an argument error with context
    pub fn argument_with_context<S: Into<String>>(message: S, context: ErrorContext) -> Self {
        Self::config_msg("Argument error", message, context)
    }

    /// Create a reference error
    pub fn reference<S: Into<String>>(message: S) -> Self {
        Self::config_msg("Reference error", message, ErrorContext::new())
    }

    /// Create a reference error with context
    pub fn reference_with_context<S: Into<String>>(message: S, context: ErrorContext) -> Self {
        Self::config_msg("Reference error", message, context)
    }

    /// Create a read error
    pub fn read<S: Into<String>>(message: S) -> Self {
        Self::algorithm_msg("Read error", message, ErrorContext::new())
    }

    /// Create a read error with context
    pub fn read_with_context<S: Into<String>>(message: S, context: ErrorContext) -> Self {
        Self::algorithm_msg("Read error", message, context)
    }

    /// Create a variant error
    pub fn variant<S: Into<String>>(message: S) -> Self {
        Self::algorithm_msg("Variant calling error", message, ErrorContext::new())
    }

    /// Create a variant error with context
    pub fn variant_with_context<S: Into<String>>(message: S, context: ErrorContext) -> Self {
        Self::algorithm_msg("Variant calling error", message, context)
    }

    /// Create a statistical error
    pub fn statistical<S: Into<String>>(message: S) -> Self {
        Self::algorithm_msg("Statistical error", message, ErrorContext::new())
    }

    /// Create a statistical error with context
    pub fn statistical_with_context<S: Into<String>>(message: S, context: ErrorContext) -> Self {
        Self::algorithm_msg("Statistical error", message, context)
    }

    /// Create a quality error
    pub fn quality<S: Into<String>>(message: S) -> Self {
        Self::algorithm_msg("Quality score error", message, ErrorContext::new())
    }

    /// Create a quality error with context
    pub fn quality_with_context<S: Into<String>>(message: S, context: ErrorContext) -> Self {
        Self::algorithm_msg("Quality score error", message, context)
    }

    /// Create an interval error
    pub fn interval<S: Into<String>>(message: S) -> Self {
        Self::config_msg("Interval error", message, ErrorContext::new())
    }

    /// Create an interval error with context
    pub fn interval_with_context<S: Into<String>>(message: S, context: ErrorContext) -> Self {
        Self::config_msg("Interval error", message, context)
    }

    /// Create an assembly error
    pub fn assembly<S: Into<String>>(message: S) -> Self {
        Self::algorithm_msg("Assembly error", message, ErrorContext::new())
    }

    /// Create an assembly error with context
    pub fn assembly_with_context<S: Into<String>>(message: S, context: ErrorContext) -> Self {
        Self::algorithm_msg("Assembly error", message, context)
    }

    /// Create an alignment error
    pub fn alignment<S: Into<String>>(message: S) -> Self {
        Self::algorithm_msg("Alignment error", message, ErrorContext::new())
    }

    /// Create an alignment error with context
    pub fn alignment_with_context<S: Into<String>>(message: S, context: ErrorContext) -> Self {
        Self::algorithm_msg("Alignment error", message, context)
    }

    /// Create a genotype error
    pub fn genotype<S: Into<String>>(message: S) -> Self {
        Self::algorithm_msg("Genotype error", message, ErrorContext::new())
    }

    /// Create a genotype error with context
    pub fn genotype_with_context<S: Into<String>>(message: S, context: ErrorContext) -> Self {
        Self::algorithm_msg("Genotype error", message, context)
    }

    /// Create a haplotype error
    pub fn haplotype<S: Into<String>>(message: S) -> Self {
        Self::algorithm_msg("Haplotype error", message, ErrorContext::new())
    }

    /// Create a haplotype error with context
    pub fn haplotype_with_context<S: Into<String>>(message: S, context: ErrorContext) -> Self {
        Self::algorithm_msg("Haplotype error", message, context)
    }

    /// Create a generic error
    pub fn generic<S: Into<String>>(message: S) -> Self {
        Self::algorithm_msg("Error", message, ErrorContext::new())
    }

    /// Create a generic error with context
    pub fn generic_with_context<S: Into<String>>(message: S, context: ErrorContext) -> Self {
        Self::algorithm_msg("Error", message, context)
    }

    /// Create an I/O error from std::io::Error
    pub fn io<S: Into<String>>(message: S, source: std::io::Error) -> Self {
        Self::Io {
            message: message.into(),
            source,
            context: boxed_ctx(ErrorContext::new()),
        }
    }

    /// Create an I/O error from std::io::Error with context
    pub fn io_with_context<S: Into<String>>(
        message: S,
        source: std::io::Error,
        context: ErrorContext,
    ) -> Self {
        Self::Io {
            message: message.into(),
            source,
            context: boxed_ctx(context),
        }
    }

    /// Map a non-`std::io::Error` failure (e.g. htslib) into [`GatkError::Io`].
    pub fn io_message<S: Into<String>>(message: S) -> Self {
        Self::io(message, std::io::Error::other("I/O failure"))
    }

    /// User-facing argument error with `parameter` / `reason` context keys.
    pub fn invalid_argument(param: impl Into<String>, message: impl Into<String>) -> Self {
        let param = param.into();
        let message = message.into();
        Self::argument_with_context(
            message.clone(),
            ErrorContext::new()
                .add_context("parameter", param)
                .add_context("reason", message),
        )
    }

    /// User-facing configuration error with `parameter` / `reason` context keys.
    pub fn invalid_configuration(param: impl Into<String>, message: impl Into<String>) -> Self {
        let param = param.into();
        let message = message.into();
        Self::configuration_with_context(
            message.clone(),
            ErrorContext::new()
                .add_context("parameter", param)
                .add_context("reason", message),
        )
    }

    /// Classify this error for CLI exit codes and reporting.
    /// User: argument / configuration / validation / interval / parse / reference.
    /// Io: `Io`.
    /// Internal: algorithm / file-format / generic / …
    pub fn classification(&self) -> GatkErrorClass {
        match self {
            Self::Config { .. } | Self::Parse { .. } => GatkErrorClass::User,
            Self::Io { .. } => GatkErrorClass::Io,
            Self::Algorithm { .. } => GatkErrorClass::Internal,
        }
    }

    /// True for user-facing argument / configuration / validation failures.
    #[inline]
    pub fn is_user_facing(&self) -> bool {
        self.classification() == GatkErrorClass::User
    }

    /// Display message plus tool / operation / position when present (CLI diagnostics).
    /// Does not alter [`std::fmt::Display`] / thiserror templates used by unit tests.
    pub fn display_with_context(&self) -> String {
        let mut out = self.to_string();
        let ctx = self.context();
        if let Some(tool) = &ctx.tool {
            out.push_str(&format!(" [tool={tool}]"));
        }
        if let Some(op) = &ctx.operation {
            out.push_str(&format!(" [operation={op}]"));
        }
        if let Some((contig, pos)) = &ctx.position {
            out.push_str(&format!(" [at {contig}:{pos}]"));
        }
        if let Some(param) = ctx.context.get("parameter") {
            out.push_str(&format!(" [parameter={param}]"));
        }
        out
    }

    /// Get the error context
    pub fn context(&self) -> &ErrorContext {
        match self {
            Self::Io { context, .. }
            | Self::Parse { context, .. }
            | Self::Config { context, .. }
            | Self::Algorithm { context, .. } => context,
        }
    }

    /// Get mutable reference to error context
    pub fn context_mut(&mut self) -> &mut ErrorContext {
        match self {
            Self::Io { context, .. }
            | Self::Parse { context, .. }
            | Self::Config { context, .. }
            | Self::Algorithm { context, .. } => context,
        }
    }
}

/// Process exit code for `gatk-rs` aligned with common GATK/Picard CLI conventions.
/// `2` — [`GatkErrorClass::User`] (Picard `USER_ERROR` analog)
/// `3` — [`GatkErrorClass::Io`]
/// `1` — [`GatkErrorClass::Internal`]
pub fn gatk_cli_exit_code(err: &GatkError) -> i32 {
    match err.classification() {
        GatkErrorClass::User => 2,
        GatkErrorClass::Io => 3,
        GatkErrorClass::Internal => 1,
    }
}

/// Result type for GATK-RS operations
pub type GatkResult<T> = Result<T, GatkError>;

impl From<std::io::Error> for GatkError {
    fn from(source: std::io::Error) -> Self {
        GatkError::io("I/O error", source)
    }
}

/// Error context trait for adding additional information to errors.
/// Preserves the original [`GatkError`] variant (and thus [`GatkError::classification`]).
pub trait ErrorContextExt<T> {
    /// Attach an operation label without changing the error variant.
    fn with_context(self, context: &str) -> GatkResult<T>;

    /// Attach tool + operation labels without changing the error variant.
    fn with_tool_context(self, tool: &str, context: &str) -> GatkResult<T>;

    /// Attach genomic position + operation labels without changing the error variant.
    fn with_position_context(self, position: (String, u64), context: &str) -> GatkResult<T>;
}

impl<T, E> ErrorContextExt<T> for Result<T, E>
where
    E: Into<GatkError>,
{
    fn with_context(self, context: &str) -> GatkResult<T> {
        self.map_err(|e| {
            let mut err = e.into();
            let ctx = err.context_mut();
            if ctx.operation.is_none() {
                ctx.operation = Some(context.to_string());
            }
            ctx.context
                .insert("additional_context".to_string(), context.to_string());
            err
        })
    }

    fn with_tool_context(self, tool: &str, context: &str) -> GatkResult<T> {
        self.map_err(|e| {
            let mut err = e.into();
            let ctx = err.context_mut();
            ctx.tool = Some(tool.to_string());
            if ctx.operation.is_none() {
                ctx.operation = Some(context.to_string());
            }
            ctx.context
                .insert("additional_context".to_string(), context.to_string());
            err
        })
    }

    fn with_position_context(self, position: (String, u64), context: &str) -> GatkResult<T> {
        self.map_err(|e| {
            let mut err = e.into();
            let ctx = err.context_mut();
            ctx.position = Some(position);
            if ctx.operation.is_none() {
                ctx.operation = Some(context.to_string());
            }
            ctx.context
                .insert("additional_context".to_string(), context.to_string());
            err
        })
    }
}

/// Macro for creating errors with formatted strings
#[macro_export]
macro_rules! gatk_error {
    (Validation, $($arg:tt)*) => {
        $crate::error::GatkError::validation(format!($($arg)*))
    };
    (Algorithm, $($arg:tt)*) => {
        $crate::error::GatkError::algorithm(format!($($arg)*))
    };
    ($variant:ident, $($arg:tt)*) => {
        $crate::error::GatkError::algorithm(format!($($arg)*))
    };
}

/// Macro for creating validation errors
#[macro_export]
macro_rules! validation_error {
    ($($arg:tt)*) => {
        $crate::error::GatkError::validation(format!($($arg)*))
    };
}

/// Macro for creating algorithm errors
#[macro_export]
macro_rules! algorithm_error {
    ($($arg:tt)*) => {
        $crate::error::GatkError::algorithm(format!($($arg)*))
    };
}

/// Macro for creating errors with context
#[macro_export]
macro_rules! gatk_error_with_context {
    ($context:expr, $($arg:tt)*) => {
        $crate::error::GatkError::algorithm_with_context(format!($($arg)*), $context)
    };
}

/// Macro for creating validation errors with tool context
#[macro_export]
macro_rules! validation_error_with_tool {
    ($tool:expr, $($arg:tt)*) => {
        $crate::error::GatkError::validation_with_context(
            format!($($arg)*),
            $crate::error::ErrorContext::new().with_tool($tool)
        )
    };
}

/// Macro for creating errors with genomic position context
#[macro_export]
macro_rules! gatk_error_at_position {
    ($pos:expr, $($arg:tt)*) => {
        $crate::error::GatkError::algorithm_with_context(
            format!($($arg)*),
            $crate::error::ErrorContext::new().with_position($pos)
        )
    };
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_error_creation() {
        let err = GatkError::file_format("Invalid VCF format");
        assert!(matches!(err, GatkError::Algorithm { .. }));
        assert!(err.to_string().contains("Invalid VCF format"));

        let err = GatkError::parse("Expected integer", 42);
        assert!(matches!(err, GatkError::Parse { line: 42, .. }));
    }

    #[test]
    fn with_context_preserves_variant_and_classification() {
        let result: Result<i32, GatkError> = Err(GatkError::argument("bad ploidy"));
        let err = result.with_context("activity scoring").unwrap_err();
        assert!(matches!(err, GatkError::Config { .. }));
        assert_eq!(err.classification(), GatkErrorClass::User);
        assert_eq!(err.context().operation.as_deref(), Some("activity scoring"));
        assert!(err
            .display_with_context()
            .contains("operation=activity scoring"));
    }

    #[test]
    fn with_tool_context_preserves_argument_variant() {
        let result: Result<i32, GatkError> = Err(GatkError::argument("x"));
        let err = result
            .with_tool_context("HaplotypeCaller", "parse args")
            .unwrap_err();
        assert!(matches!(err, GatkError::Config { .. }));
        assert_eq!(err.context().tool.as_deref(), Some("HaplotypeCaller"));
        assert_eq!(gatk_cli_exit_code(&err), 2);
    }

    #[test]
    fn test_error_macros() {
        let err = gatk_error!(Validation, "Invalid value: {}", 42);
        match err {
            GatkError::Config { message, .. } => {
                assert!(message.contains("Invalid value: 42"));
            }
            _ => panic!("Expected Config error"),
        }
    }

    #[test]
    fn gatk_cli_exit_code_mappings() {
        assert_eq!(gatk_cli_exit_code(&GatkError::argument("x")), 2);
        assert_eq!(gatk_cli_exit_code(&GatkError::configuration("x")), 2);
        assert_eq!(
            gatk_cli_exit_code(&GatkError::io(
                "x",
                std::io::Error::new(std::io::ErrorKind::Other, "e")
            )),
            3
        );
        assert_eq!(gatk_cli_exit_code(&GatkError::algorithm("x")), 1);
        assert!(GatkError::argument("x").is_user_facing());
        assert!(!GatkError::algorithm("x").is_user_facing());
    }

    #[test]
    fn invalid_argument_attaches_parameter_metadata() {
        let err = GatkError::invalid_argument("kmer_size", "kmer_size must be ≥ 2, got 1");
        assert!(matches!(err, GatkError::Config { .. }));
        assert_eq!(
            err.context().context.get("parameter").map(String::as_str),
            Some("kmer_size")
        );
        assert!(err.display_with_context().contains("parameter=kmer_size"));
    }
}
