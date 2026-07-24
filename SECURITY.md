# Security Policy

## Reporting a Vulnerability

Please do **not** open public GitHub issues for security vulnerabilities.

Report vulnerabilities privately to **contact@synapticfour.com** with:
- affected repository and version/commit
- reproduction steps or proof-of-concept
- impact assessment

We will acknowledge receipt as quickly as possible, triage severity, and coordinate a responsible disclosure timeline.

Alternatively, use GitHub Security Advisories if available: **Security** tab → **Report a vulnerability**.

## Scope and Guarantees

This project is maintained on a best-effort basis (maturity: Alpha). Security documentation and test coverage improve over time, but no absolute security guarantee is provided.

In scope: the gatk-rs workspace crates, CLI, and CI infrastructure for this repository.

Out of scope: pinned third-party oracles (e.g. Java GATK jars), GIAB / reference data, and unmodified upstream tools invoked by equivalence harnesses.
