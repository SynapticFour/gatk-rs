//! Process-global optional NDJSON sink for semantic checkpoints.

use super::schema::{RegionKey, SemanticStage, SemanticTraceEvent, SCHEMA_ID, TRACE_IMPL_RUST};
use crate::runtime_config::RuntimeConfig;
use serde_json::Value;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::Mutex;

static ENABLED: AtomicBool = AtomicBool::new(false);
static SEQ: AtomicU64 = AtomicU64::new(0);
static SINK: Mutex<Option<NdjsonFileSink>> = Mutex::new(None);

/// Opaque handle reserved for future non-global sinks.
#[derive(Debug)]
pub struct TraceSinkHandle;

struct NdjsonFileSink {
    writer: BufWriter<File>,
}

impl NdjsonFileSink {
    fn open(path: PathBuf) -> std::io::Result<Self> {
        let file = File::create(path)?;
        Ok(Self {
            writer: BufWriter::new(file),
        })
    }

    fn write_event(&mut self, event: &SemanticTraceEvent) -> std::io::Result<()> {
        serde_json::to_writer(&mut self.writer, event)?;
        self.writer.write_all(b"\n")?;
        Ok(())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.writer.flush()
    }
}

/// Whether semantic tracing is active for this process.
#[inline]
pub fn is_enabled() -> bool {
    ENABLED.load(Ordering::Relaxed)
}

/// Initialize from [`RuntimeConfig`]. No-op if already enabled or path unset.
pub fn try_init_from_runtime(cfg: &RuntimeConfig) {
    if is_enabled() {
        return;
    }
    let Some(path) = cfg.debug.semantic_trace_path.as_ref() else {
        return;
    };
    if path.is_empty() {
        return;
    }
    let Ok(mut guard) = SINK.lock() else {
        return;
    };
    if guard.is_some() {
        ENABLED.store(true, Ordering::Relaxed);
        return;
    }
    match NdjsonFileSink::open(PathBuf::from(path)) {
        Ok(sink) => {
            *guard = Some(sink);
            SEQ.store(0, Ordering::Relaxed);
            ENABLED.store(true, Ordering::Relaxed);
        }
        Err(e) => {
            eprintln!("GATK_RS_SEMANTIC_TRACE: failed to open {path}: {e}");
            ENABLED.store(false, Ordering::Relaxed);
        }
    }
}

/// Emit one event to the global sink (no-op when disabled).
pub fn emit(stage: SemanticStage, region: Option<RegionKey>, payload: Value) {
    if !is_enabled() {
        return;
    }
    let Ok(mut guard) = SINK.lock() else {
        return;
    };
    let Some(sink) = guard.as_mut() else {
        return;
    };
    let seq = SEQ.fetch_add(1, Ordering::Relaxed) + 1;
    let event = SemanticTraceEvent {
        schema: SCHEMA_ID.to_string(),
        seq,
        impl_name: TRACE_IMPL_RUST.to_string(),
        stage,
        region,
        payload,
    };
    if let Err(e) = sink.write_event(&event) {
        eprintln!("GATK_RS_SEMANTIC_TRACE: write failed: {e}");
    }
}

#[cfg(test)]
pub(super) fn reset_for_tests() {
    ENABLED.store(false, Ordering::Relaxed);
    SEQ.store(0, Ordering::Relaxed);
    if let Ok(mut guard) = SINK.lock() {
        *guard = None;
    }
}

#[cfg(test)]
pub(super) fn flush_for_tests() {
    if let Ok(mut guard) = SINK.lock() {
        if let Some(sink) = guard.as_mut() {
            let _ = sink.flush();
        }
    }
}
