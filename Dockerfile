# Multi-stage build for GATK-RS (release CLI binary).
FROM rust:1.85-bookworm AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    pkg-config \
    libssl-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    curl \
    cmake \
    clang \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Full context (see .dockerignore). Avoid the old Cargo.toml-only cache trick —
# [[bench]]/[[example]] entries require their source files to exist for cargo parse.
COPY . .

ENV CARGO_TERM_COLOR=always \
    CARGO_REGISTRIES_CRATES_IO_PROTOCOL=sparse

RUN cargo build --release -p gatk-cli --bin gatk-rs \
    && strip target/release/gatk-rs

FROM debian:bookworm-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    libssl3 \
    zlib1g \
    libbz2-1.0 \
    liblzma5 \
    curl \
    && rm -rf /var/lib/apt/lists/*

RUN useradd -m -u 1000 gatkrs
WORKDIR /app

COPY --from=builder /app/target/release/gatk-rs /usr/local/bin/gatk-rs
RUN chmod +x /usr/local/bin/gatk-rs

USER gatkrs

LABEL org.opencontainers.image.title="GATK-RS"
LABEL org.opencontainers.image.description="A native Rust implementation of the Genome Analysis Toolkit (GATK)"
LABEL org.opencontainers.image.source="https://github.com/SynapticFour/gatk-rs"
LABEL org.opencontainers.image.licenses="Apache-2.0"
LABEL org.opencontainers.image.version="0.1.0"

HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD gatk-rs --version || exit 1

ENTRYPOINT ["gatk-rs"]
CMD ["--help"]
