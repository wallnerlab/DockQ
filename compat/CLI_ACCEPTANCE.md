# CLI differential acceptance

The complete captured run is in [cli_comparison.txt](cli_comparison.txt).
It compares the release Rust executable with the pinned Python reference described
in [README.md](README.md), using one worker and the same host. Python's progress
executor is adapted to a serial map in the harness because the sandbox restricts
process-manager sockets; scientific source functions are unchanged.

## Outcome

- All numeric JSON comparisons in the 12 scoring scenarios pass: contact counts,
  native/model assignments, scores, and RMSDs. Top-level best mappings agree.
- Two additional Python-fallback invocations match normalized text.
- A generated no-interface example fails appropriately in both executables.
- **12 of 15 full CLI checks pass.** Three 1EXB cases deliberately remain marked
  FAIL by the strict text comparator: short automatic, partial, and complete mapping.
  Their only remaining failing comparison is one printed LRMSD rounding across
  the 0.0005 Å boundary: Python prints `0.000`, Rust prints `0.001`.
  These are known documented differences, not a claim of exact text parity.

For the affected 1EXB FH:EG interface, Python LRMSD is 0.000497884823376 Å
and Rust LRMSD is 0.000502755378778 Å: a difference of 4.8705554e-6 Å across
the three-decimal rounding boundary.

The largest observed LRMSD JSON difference is **4.44393463e-5 Å**. The largest
per-interface DockQ difference is **8.95861977e-7**. Counts and contact-derived
fractions agree exactly. The numerical gates remain absolute 1e-4 Å for RMSDs
and 1e-5 for other floating metrics. They have not been relaxed to hide text
rounding at three decimal places.

## Python diagnostic cache behavior

Python's cached `run_on_chains` returns mutable result dictionaries.
`run_on_all_native_interfaces` adds `chain_map` to those dictionaries for every
candidate; subsequent candidates can overwrite a dictionary still referenced by
the stored best result. This changes nested diagnostic `chain_map` fields without
changing the actual top-level best mapping or scientific metrics.

Rust keeps per-interface diagnostics consistent with its selected mapping. The
harness explicitly labels this known baseline quirk and excludes only those
nested diagnostic fields when the top-level mappings agree. It still compares
top-level mappings, each interface's actual chain assignments, ligand atom
mappings, counts, and metrics. It does not normalize away the three text failures.

## Illustrative performance

Single observations including process startup, parsing, mapping, scoring, and
output; these are not repeated statistical benchmarks or portability guarantees:

| Scenario | Python seconds | Rust seconds | Python peak RSS | Rust peak RSS |
|---|---:|---:|---:|---:|
| 1A2K automatic mapping | 0.37 | 0.07 | 57.8 MB | 14.6 MB |
| 1EXB automatic mapping | 6.76 | 0.87 | 93.2 MB | 14.6 MB |

Reproduce with a suitable reference environment:

```sh
cargo build --release --locked --manifest-path rust/Cargo.toml
python compat/compare_cli.py --python /path/to/reference/python --memory
```

The strict comparator currently exits 1 for the documented text differences.
Rust's own test suite passes independently of this compatibility report.
