# Rust artifact preparation

The Rust CLI runs without Python at runtime. This contribution proposes it as
an optional implementation alongside the existing Python package.

## Validation and artifact workflow

Run the manual **Rust release artifacts** workflow to check formatting, Clippy,
all targets and features, tests, a locked release build, and a scoring smoke
example on Linux x86_64, macOS x86_64, and Windows x86_64. The workflow generates
dependency license notices and uploads an archive for each platform.

Each archive contains:

- The platform's DockQ executable.
- The original MIT LICENSE and THIRD_PARTY_LICENSES.txt.
- The Rust README, build revision and compiler metadata, and SHA256SUMS.txt.

The workflow only uploads Actions artifacts. Upstream release names, publication,
and long-term distribution remain decisions for the DockQ maintainers.

## Existing validation evidence

See TEST_REPORT.md for hosted validation of the same Rust implementation and
the independently maintained fork release. Those binaries are separate from
official DockQ releases. A new contribution commit needs its own CI results.

## Reproduce checks locally

From the repository root, with stable Rust and Python 3 available:

    cargo fmt --manifest-path rust/Cargo.toml --all -- --check
    cargo clippy --locked --manifest-path rust/Cargo.toml --all-targets --all-features -- -D warnings
    cargo test --locked --manifest-path rust/Cargo.toml --all-targets --all-features
    cargo build --locked --release --manifest-path rust/Cargo.toml
    python scripts/generate_third_party_licenses.py --output THIRD_PARTY_LICENSES.txt
    rust/target/release/DockQ examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb --mapping AB:AB --short

Read README.md and MIGRATION_STATUS.md for compatibility limits. Python bindings,
exhaustive Biopython parser compatibility, and bitwise numerical identity are
outside the current contribution.
