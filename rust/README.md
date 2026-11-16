# DockQ Rust

DockQ Rust is a command-line program and Rust library for assessing protein,
nucleic-acid, peptide, and small-molecule docking models. This contribution
proposes it as an optional implementation alongside the Python package. It runs
without Python at runtime and is not a replacement for the Python API.

## Install a release binary

Independently maintained fork binaries are available from the
[DockQ Rust v2.1.3 GitHub release](https://github.com/enriquead1/DockQ/releases/tag/rust-v2.1.3).
Extract it, read BUILD-INFO.txt, and verify SHA256SUMS.txt before running the
executable.

Release bundles target Linux x86_64, macOS x86_64, and Windows x86_64. CI checks
formatting, Clippy, tests, release builds, and a representative scoring example
on all three operating systems. See TEST_REPORT.md for the latest results.

## Build from source

Use a stable Rust toolchain and the committed lockfile:

    cargo build --locked --release --manifest-path rust/Cargo.toml

Run formatting, lint, and tests from the repository root:

    cargo fmt --manifest-path rust/Cargo.toml --all -- --check
    cargo clippy --locked --manifest-path rust/Cargo.toml --all-targets --all-features -- -D warnings
    cargo test --locked --manifest-path rust/Cargo.toml --all-targets --all-features

The binary is rust/target/release/DockQ on Linux and macOS and
rust/target/release/DockQ.exe on Windows.

## Example

    rust/target/release/DockQ \
      examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb \
      --mapping AB:AB --short

The mapping syntax is MODELCHAINS:NATIVECHAINS. Mapping currently accepts
single-character chain identifiers in the CLI. Use --small_molecule to score
hetero-residue chains and --json results.json to write structured output. Run
DockQ --help for all options.

## Scope and compatibility

The CLI reads PDB and mmCIF coordinate files, including gzip-compressed files.
It scores protein and nucleic-acid interfaces, CAPRI peptide interfaces, and
small-molecule interfaces; searches chain mappings; and writes text or JSON
results.

The reference corpus checks integer counts exactly, RMSD values within 1e-4 Å,
and other floating-point metrics within 1e-5. This is numerical agreement on the
included cases, not bitwise identity or a guarantee for every structure accepted
by Biopython. The recorded CLI comparison has three strict text differences:
Python prints LRMSD 0.000 while Rust prints 0.001 for one interface whose
values lie on opposite sides of the three-decimal rounding boundary. The
underlying values differ by about 4.9e-6 Å. See MIGRATION_STATUS.md and
../compat/CLI_ACCEPTANCE.md for details.

The Rust library uses native Rust data structures. Existing Python callers that
pass Biopython objects are not supported by this CLI and are not preserved as a
Rust API. Parser recovery and detailed diagnostics for malformed inputs can also
differ from Biopython.

## Documentation

- Repository overview and quick start: ../README.md
- Migration status and compatibility details: MIGRATION_STATUS.md
- Cross-platform validation: TEST_REPORT.md
- Release contents and workflow: RELEASE.md
- Reproducible all-against-all benchmark: ../benchmarks/all_vs_all/BENCHMARK_FINDINGS.md

If you use DockQ, cite the DockQ v2 paper at
https://doi.org/10.1093/bioinformatics/btae586.
