# Rust implementation and acceptance status

The Rust CLI and library are proposed as an optional implementation alongside
the Python package. The existing Python source, packaging, and README remain
in place. The Rust API uses native Rust types and does not replace the Python
API.

## Implemented workflows

The CLI loads PDB and mmCIF coordinate files, including gzip-compressed files.
It scores protein and nucleic-acid interfaces, CAPRI peptide interfaces, and
small-molecule interfaces; searches mappings; supports parallel candidate
search; and writes text or JSON output. The executable does not invoke Python
or load Python libraries.

## Validation

The included Rust suite covers unit behavior, Python reference fixtures, and
deterministic mapping search. CI runs formatting, Clippy with warnings denied,
all targets and features, a locked release build, and a CLI scoring smoke test
on Linux x86_64, macOS x86_64, and Windows x86_64. See TEST_REPORT.md for the
latest commit and platform-by-platform outcomes.

The differential CLI corpus has 12 of 15 strict text comparisons passing.
All numeric JSON comparisons pass the documented tolerances. Three strict text
comparisons differ only in LRMSD rounding at a 0.0005 Å boundary: Python prints
0.000 and Rust prints 0.001. For the affected interface, the values differ by
about 4.9e-6 Å. This is documented behavior and is not hidden by widening the
text comparator.

## Compatibility boundaries

- Contact counts, residue and atom correspondence, and selected top-level
  mappings are checked on the supplied corpus. RMSDs use a 1e-4 Å absolute
  tolerance; other floating metrics use 1e-5. These checks do not promise
  bitwise equality or exhaustive agreement on every structure.
- The legacy CLI mapping syntax supports single-character chain identifiers.
  The Rust library keeps longer mmCIF identifiers.
- Typed errors and recovery from malformed records can differ from Biopython.
  Full mmCIF dictionary validation and unrelated Biopython metadata are outside
  this implementation's scope.
- The public Rust API uses native Rust structures. Python bindings and the
  Biopython object API are not included.
- The benchmark data manifest is checked in, but the 6,600-structure CIF corpus
  is not redistributed with this repository.

## Reproduce validation

From the repository root:

    cargo fmt --manifest-path rust/Cargo.toml --all -- --check
    cargo clippy --locked --manifest-path rust/Cargo.toml --all-targets --all-features -- -D warnings
    cargo test --locked --manifest-path rust/Cargo.toml --all-targets --all-features
    cargo build --locked --release --manifest-path rust/Cargo.toml
    rust/target/release/DockQ examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb --mapping AB:AB --short

The cross-platform jobs run these checks on hosted Linux, macOS, and Windows
runners. Regenerating Python comparison fixtures requires the reference
environment described in ../compat/README.md.
