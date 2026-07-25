# Multi-stage build for GATK-RS
FROM rust:1.75-slim AS builder

# Install build dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    pkg-config \
    libssl-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy Cargo files
COPY Cargo.toml Cargo.lock ./
COPY gatk-core/Cargo.toml ./gatk-core/
COPY gatk-cli/Cargo.toml ./gatk-cli/
COPY gatk-common/Cargo.toml ./gatk-common/
COPY gatk-haplotypecaller/Cargo.toml ./gatk-haplotypecaller/
COPY gatk-rs-equiv/Cargo.toml ./gatk-rs-equiv/

# Create dummy source files to cache dependencies
RUN mkdir -p src gatk-core/src gatk-cli/src gatk-common/src gatk-haplotypecaller/src gatk-rs-equiv/src
RUN echo "fn main() {}" > src/main.rs
RUN echo "pub fn dummy() {}" > gatk-core/src/lib.rs
RUN echo "pub fn dummy() {}" > gatk-cli/src/lib.rs
RUN echo "pub fn dummy() {}" > gatk-common/src/lib.rs
RUN echo "pub fn dummy() {}" > gatk-haplotypecaller/src/lib.rs
RUN echo "fn main() {}" > gatk-rs-equiv/src/main.rs

# Build dependencies
RUN cargo build --workspace --release
RUN rm -rf src gatk-core/src gatk-cli/src gatk-common/src gatk-haplotypecaller/src gatk-rs-equiv/src

# Copy source code
COPY . .

# Build the application
RUN cargo build --workspace --release

# Runtime stage
FROM debian:bullseye-slim

# Install runtime dependencies
RUN apt-get update && apt-get install -y \
    ca-certificates \
    libssl1.1 \
    zlib1g \
    libbz2-1.0 \
    liblzma5 \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Create non-root user
RUN useradd -m -u 1000 gatkrs

# Set working directory
WORKDIR /app

# Copy binary from builder stage
COPY --from=builder /app/target/release/gatk-rs /usr/local/bin/gatk-rs

# Set permissions
RUN chmod +x /usr/local/bin/gatk-rs

# Switch to non-root user
USER gatkrs

# Add metadata
LABEL org.opencontainers.image.title="GATK-RS"
LABEL org.opencontainers.image.description="A native Rust implementation of the Genome Analysis Toolkit (GATK)"
LABEL org.opencontainers.image.source="https://github.com/SynapticFour/gatk-rs"
LABEL org.opencontainers.image.licenses="Apache-2.0"
LABEL org.opencontainers.image.version="0.1.0"

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD gatk-rs --version || exit 1

# Set entrypoint
ENTRYPOINT ["gatk-rs"]

# Default command shows help
CMD ["--help"]
