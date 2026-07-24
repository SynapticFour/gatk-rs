//! Logging utilities for GATK-RS
//! This module provides logging functionality that mimics GATK's log4j-based logging
//! using Rust's tracing ecosystem.

use tracing::{info, Level};
use tracing_appender::{non_blocking, rolling};
use tracing_subscriber::{fmt, layer::SubscriberExt, util::SubscriberInitExt, EnvFilter, Layer};

/// Logging configuration
#[derive(Debug, Clone)]
pub struct LoggingConfig {
    /// Log level
    pub level: Level,
    /// Whether to log to file
    pub log_to_file: bool,
    /// Log file directory
    pub log_dir: Option<String>,
    /// Log file name prefix
    pub log_file_prefix: String,
    /// Whether to include timestamps
    pub include_timestamps: bool,
    /// Whether to include thread names
    pub include_threads: bool,
    /// Whether to use colors in console output
    pub use_colors: bool,
}

impl Default for LoggingConfig {
    fn default() -> Self {
        Self {
            level: Level::INFO,
            log_to_file: false,
            log_dir: None,
            log_file_prefix: "gatk-rs".to_string(),
            include_timestamps: true,
            include_threads: false,
            use_colors: true,
        }
    }
}

/// Initialize logging with the given configuration
pub fn init_logging(config: LoggingConfig) -> Result<(), Box<dyn std::error::Error>> {
    let env_filter = EnvFilter::try_from_default_env()
        .unwrap_or_else(|_| EnvFilter::new(config.level.to_string()));

    let mut layers = Vec::new();

    // Console layer
    let console_layer = fmt::layer()
        .with_target(false)
        .with_thread_ids(config.include_threads)
        .with_thread_names(config.include_threads)
        .with_ansi(config.use_colors);

    if config.include_timestamps {
        layers.push(
            console_layer
                .with_timer(fmt::time::ChronoUtc::rfc_3339())
                .boxed(),
        );
    } else {
        layers.push(console_layer.boxed());
    }

    // File layer (if configured)
    if config.log_to_file {
        if let Some(log_dir) = &config.log_dir {
            std::fs::create_dir_all(log_dir)?;

            let file_appender = rolling::daily(log_dir, &config.log_file_prefix);
            let (non_blocking, _guard) = non_blocking(file_appender);

            let file_layer = fmt::layer()
                .with_writer(non_blocking)
                .with_target(true)
                .with_thread_ids(config.include_threads)
                .with_thread_names(config.include_threads)
                .with_ansi(false)
                .with_timer(fmt::time::ChronoUtc::rfc_3339());

            layers.push(file_layer.boxed());
        }
    }

    tracing_subscriber::registry()
        .with(env_filter)
        .with(layers)
        .init();

    info!("GATK-RS logging initialized with level: {}", config.level);

    Ok(())
}

/// Initialize logging with default configuration
pub fn init_default_logging() -> Result<(), Box<dyn std::error::Error>> {
    init_logging(LoggingConfig::default())
}

/// Initialize logging from environment variables
pub fn init_logging_from_env() -> Result<(), Box<dyn std::error::Error>> {
    let config = LoggingConfig {
        level: std::env::var("GATK_LOG_LEVEL")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(Level::INFO),
        log_to_file: std::env::var("GATK_LOG_TO_FILE")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(false),
        log_dir: std::env::var("GATK_LOG_DIR").ok(),
        log_file_prefix: std::env::var("GATK_LOG_FILE_PREFIX")
            .ok()
            .unwrap_or_else(|| "gatk-rs".to_string()),
        include_timestamps: std::env::var("GATK_LOG_TIMESTAMPS")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(true),
        include_threads: std::env::var("GATK_LOG_THREADS")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(false),
        use_colors: std::env::var("GATK_LOG_COLORS")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(true),
    };

    init_logging(config)
}

/// Initialize a lightweight `tracing` subscriber from **`GATK_RS_HC_TRACE`**.
/// No-op if **`RUST_LOG`** is already set (standard `tracing-subscriber` env wins) or if the
/// variable is unset, empty, **`0`**, or **`off`**. Returns **`true`** if a subscriber was installed.
/// Numeric shorthand (when the value does not contain `=`): **`1`** → `error`, **`2`** → `warn`,
/// **`3`** → `info`, **`4`** → `debug`, **`5`** / **`trace`** → `trace`. Any other token is applied
/// as `gatk_haplotypecaller=<token>` (e.g. `debug`). See `docs/ARCHITECTURE.md`.
pub fn try_init_from_gatk_rs_hc_trace() -> bool {
    if std::env::var_os("RUST_LOG").is_some() {
        return false;
    }
    let Ok(raw) = std::env::var("GATK_RS_HC_TRACE") else {
        return false;
    };
    let raw = raw.trim();
    if raw.is_empty() || raw == "0" || raw.eq_ignore_ascii_case("off") {
        return false;
    }
    let directive: String = if raw.contains('=') {
        raw.to_string()
    } else {
        let level = match raw {
            "1" => "error",
            "2" => "warn",
            "3" => "info",
            "4" => "debug",
            "5" | "trace" => "trace",
            other => other,
        };
        format!("gatk_haplotypecaller={level},gatk_core=warn,gatk_common=warn")
    };
    let Ok(filter) = EnvFilter::try_new(&directive) else {
        return false;
    };
    tracing_subscriber::fmt()
        .with_env_filter(filter)
        .with_target(true)
        .try_init()
        .is_ok()
}

/// Macro for logging tool execution start
#[macro_export]
macro_rules! log_tool_start {
    ($tool_name:expr) => {
        info!("Starting {}", $tool_name);
        info!("GATK-RS version: {}", env!("CARGO_PKG_VERSION"));
        info!("Rust version: {}", env!("RUSTC_VERSION"));
    };
}

/// Macro for logging tool execution end
#[macro_export]
macro_rules! log_tool_end {
    ($tool_name:expr) => {
        info!("Completed {}", $tool_name);
    };
}

/// Macro for logging progress
#[macro_export]
macro_rules! log_progress {
    ($message:expr) => {
        info!("{}", $message);
    };
}

/// Macro for logging warnings
#[macro_export]
macro_rules! log_warning {
    ($message:expr) => {
        warn!("{}", $message);
    };
}

/// Macro for logging errors
#[macro_export]
macro_rules! log_error {
    ($message:expr) => {
        error!("{}", $message);
    };
}

/// Macro for logging debug information
#[macro_export]
macro_rules! log_debug {
    ($message:expr) => {
        debug!("{}", $message);
    };
}

/// Macro for logging trace information
#[macro_export]
macro_rules! log_trace {
    ($message:expr) => {
        trace!("{}", $message);
    };
}

/// Progress logging utility
pub struct ProgressLogger {
    total: u64,
    current: u64,
    log_interval: u64,
    tool_name: String,
}

impl ProgressLogger {
    /// Create a new progress logger
    pub fn new(tool_name: String, total: u64) -> Self {
        Self {
            total,
            current: 0,
            log_interval: std::cmp::max(1, total / 100), // Log every 1%
            tool_name,
        }
    }

    /// Set custom log interval
    pub fn with_interval(mut self, interval: u64) -> Self {
        self.log_interval = interval;
        self
    }

    /// Increment progress and log if necessary
    pub fn increment(&mut self) {
        self.current += 1;

        if self.current % self.log_interval == 0 {
            let percentage = (self.current as f64 / self.total as f64) * 100.0;
            info!(
                "{} progress: {:.1}% ({}/{})",
                self.tool_name, percentage, self.current, self.total
            );
        }
    }

    /// Add multiple to progress
    pub fn add(&mut self, count: u64) {
        for _ in 0..count {
            self.increment();
        }
    }

    /// Finish progress logging
    pub fn finish(self) {
        info!(
            "{} completed: {}/{} (100.0%)",
            self.tool_name, self.total, self.total
        );
    }

    /// Get current progress as percentage
    pub fn progress_percentage(&self) -> f64 {
        (self.current as f64 / self.total as f64) * 100.0
    }

    /// Get current count
    pub fn current(&self) -> u64 {
        self.current
    }

    /// Get total count
    pub fn total(&self) -> u64 {
        self.total
    }
}

/// Performance logging utility
pub struct PerformanceLogger {
    tool_name: String,
    start_time: std::time::Instant,
}

impl PerformanceLogger {
    /// Create a new performance logger
    pub fn new(tool_name: String) -> Self {
        let start_time = std::time::Instant::now();
        info!("Starting performance logging for {}", tool_name);

        Self {
            tool_name,
            start_time,
        }
    }

    /// Log elapsed time
    pub fn log_elapsed(&self, checkpoint: &str) {
        let elapsed = self.start_time.elapsed();
        info!("{} - {}: {:?}", self.tool_name, checkpoint, elapsed);
    }

    /// Finish performance logging
    pub fn finish(self) {
        let elapsed = self.start_time.elapsed();
        info!("{} total execution time: {:?}", self.tool_name, elapsed);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_progress_logger() {
        let logger = ProgressLogger::new("test".to_string(), 100);

        for _i in 0..100 {
            let message = format!("Test message {}", _i);
            info!("{}", message);
        }
        assert_eq!(logger.current(), 0);
        assert_eq!(logger.total(), 100);
        assert_eq!(logger.progress_percentage(), 0.0);
    }

    #[test]
    fn test_performance_logger() {
        let logger = PerformanceLogger::new("test".to_string());
        std::thread::sleep(std::time::Duration::from_millis(10));
        logger.log_elapsed("checkpoint");
        logger.finish();
    }

    #[test]
    fn test_logging_config_default() {
        let config = LoggingConfig::default();
        assert_eq!(config.level, Level::INFO);
        assert!(!config.log_to_file);
        assert!(config.include_timestamps);
        assert!(!config.include_threads);
        assert!(config.use_colors);
    }
}
