# Python reference fixtures

Generate with a compiled installation of the baseline checkout in an isolated
Python environment, using `python compat/generate_reference.py`. The script records
the git revision, Python/dependency versions, and compiled operations path. It does
not change Python code. `reference.json` contains normalized parser summaries,
unrounded per-interface metrics for explicit mappings, alignment tie cases, and a
controlled compiled/fallback kernel discrepancy. Automatic mapping is checked
separately by end-to-end regression tests.

Counts, residue/atom order, and mappings must agree exactly. Calibrated numerical
targets are absolute 1e-4 Å for RMSD and 1e-5 for scores and other metrics.
The Python baseline performs centroid, covariance, and SVD operations on float32
arrays; Rust performs rigid geometry in float64. A measured LRMSD difference of
2.97e-5 Å on 1A2K motivates the RMSD tolerance. Coordinates and contact decisions
still use float32, and integer counts must match exactly.
Threshold and mapping differences must be investigated rather than hidden by a
larger tolerance. Runtime paths in provenance are informational.

Pinned generation environment: Python 3.12, NumPy 1.26.4, Biopython 1.85,
NetworkX 3.4.2, parallelbar 2.5, Cython 3.3.1, DockQ 2.1.3 from the recorded revision.
Pip may report unrelated inherited-package conflicts in an environment created with
`--system-site-packages`; a clean virtual environment is preferred for reproduction.

## Compatibility decisions

See [CLI_ACCEPTANCE.md](CLI_ACCEPTANCE.md) for the full 15-scenario CLI comparison,
including three strict text-rounding failures and the labeled diagnostic cache
differences. Numeric JSON comparisons pass the stated tolerances; exact text
parity is not claimed.

- The compiled implementation is the initial distance-kernel reference. Its
  early-exit condition can return a non-minimum distance: two atoms at 100 Å and
  0 Å against one at 0 Å yield 10000 Å² compiled versus 0 Å² exhaustive. Preserve
  compiled behavior during parity work and expose an exhaustive helper only if
  separately named; do not silently substitute it.
- Internal mappings preserve native-to-model insertion order. Equal-length
  receptor/ligand assignment and best-score ties depend on order.
- Parser coordinates and compiled contact-distance arithmetic use float32.
- Rust uses float64 transforms, whereas the baseline superimposer retains float32
  arithmetic when its input arrays are float32. See the calibrated tolerances above.
- Known baseline quirks and any accepted deviations must be documented here or in
  the final compatibility report. A baseline error is not evidence Rust should
  fabricate a successful result.
