# Rust all-against-all benchmark: direct vs cached contact matrices

## Workload and method

The benchmark runner uses the 88 complex/interface rows in `benchmark_all_88_results_v0.csv` as a workload manifest. It loads 100 CIF structures per complex, sorts file paths, fuses the listed chains, and attempts all unordered pairs with a fixed model/native direction. It uses one shared Rayon pool with 48 workers. Each of the 88 rows runs 5,050 pair attempts, for 444,400 attempts total. No Python scorer was run or used as a timing reference.

The direct baseline was measured before changing the scorer. After the cache implementation, the same runner compared the direct path with the cached path on the same pairs, alternating which path ran first by interface. Cache construction time is included in the cached total. Output files are written only after scoring.

The manifest CSV is included in this directory. The 6,600 mmCIF files (about 2.81 GB uncompressed) are not included in this repository. To rerun the workload, put one subdirectory per complex (for example, `benchmarks/all_vs_all/1AHW/`) here, containing that complex's CIF files.

## Results

| Measurement | Direct, before change | Direct control, after change | Cached, after change |
| --- | ---: | ---: | ---: |
| CIF loading | 59.870 s | 58.628 s | 58.628 s |
| Chain fusion | included below | 3.773 s | 3.773 s |
| Contact-matrix precomputation | — | — | 2.414 s |
| Pair scoring | included below | 104.736 s | 11.806 s |
| Fusion + precompute + scoring | 109.364 s | 108.510 s | 17.994 s |
| End to end, including loading | 169.234 s | 167.138 s | 76.622 s |

The paired direct control and cached paths ran in the same process. Including precomputation, caching was **6.03× faster** than direct scoring and reduced compute time by **83.4%**. Including CIF loading, the whole run was **2.21× faster**, from 169.234 s to 76.622 s. The direct control after the change was within 0.8% of the pre-change baseline.

All 444,400 attempts ran in both paths. There were 375,392 scored results and 69,008 `None` results in the chosen direction. The direct and cached results matched exactly for every attempted pair: **zero mismatches**, maximum DockQ difference **0.0**. The saved before/after score TSVs also had 375,392 rows each and no differing pair keys or scores.

For the formerly slow 1S78 AB vs C interface, the paired run took 2.494 s for direct pair scoring; cached precomputation took 0.045 s and cached pair scoring took 0.348 s. Including chain fusion, that is about 5.4× faster for this interface.

## Follow-up optimizations

I profiled the cached path and removed additional repeated work: contact masks/counts and interface residue lists are now built once per structure, backbone coordinates are cached in scoring order, identical sequences use a prebuilt identity mapping, and ligand RMSD applies its transform without a temporary coordinate array. Rigid-fit coordinate validation is folded into the centroid pass.

| Cached-path stage | Pair scoring | Contact/backbone precompute | Chain fusion | Fusion + precompute + scoring |
| --- | ---: | ---: | ---: | ---: |
| Initial contact-matrix reuse | 11.806 s | 2.414 s | 3.773 s | 17.994 s |
| Contact summaries | 5.552 s | 2.571 s | 3.812 s | 11.935 s |
| Final optimized path | 1.067 s | 2.693 s | 3.738 s | 7.498 s |

The final path is **2.40× faster** than the initial cached path for compute time, including all precomputation, and reduces that time by **58.3%**. Pair scoring alone fell by **91.0%** from 11.806 s to 1.067 s. Including structure loading, the final run took 66.696 s (59.198 s load + 7.498 s compute), versus 76.622 s for the initial cached run. The final SVD run produced 375,392 scored pairs matching the prior cached scores; the earlier full paired run established exact direct/cached parity across all 444,400 attempts. The current runner also records raw DockQ bits for future comparisons.

The follow-up profile used the 1S78 AB vs C interface with 48 workers and aggregated CPU time across worker threads. Geometry/RMSD used 81.1% of the instrumented scorer phase time; backbone selection used 14.6%, contact statistics 3.8%, and sequence alignment 0.5%. The next major target is the pair-specific Kabsch fitting and ligand RMSD work. A faster 3×3 fit needs careful numerical comparison because it can change DockQ values slightly.

Reproduce the final all-interface run from the repository root with:

```sh
cargo run --locked --release --offline --no-default-features \
  --manifest-path rust/Cargo.toml --example benchmark_all -- benchmarks/all_vs_all cached 48 /tmp/final-cached-scores.tsv
```

## QCP trial for interface RMSD

I tried the quaternion characteristic polynomial (QCP) for the interface fit only; the receptor fit remained on SVD. The full workload produced the same 375,392 pair keys, but **312,730 raw DockQ `f64` values differed** from the SVD baseline. At 12-decimal display precision, 3,334 values differed, with a maximum displayed difference of about `1e-12`. Pair scoring was 1.023 s with QCP versus 1.028 s with SVD, while total cached compute was 7.481 s versus 7.449 s. That is no meaningful speedup, and it does not preserve exact scores, so the QCP path was removed and the default remains SVD.

The benchmark TSV now includes a `DockQ_bits` column containing each score's raw `f64` representation, so future solver trials can check exact equality rather than relying on the displayed decimal precision.

## Implementation

[`contact_matrix`](../../rust/src/scoring/polymer.rs) builds an all-residue distance matrix once for each structure and interface. [`score_interface_with_cached_contacts`](../../rust/src/scoring/polymer.rs) borrows those matrices when the alignment maps every residue to the same index. It falls back to pair-specific distance calculation when the mapping is not the identity, preserving the existing path for gaps or reordered residues. The score formula and distance arithmetic are unchanged.

The full benchmark and direct/cached comparison are in [`benchmark_all.rs`](../../rust/examples/benchmark_all.rs). For example, run the direct version with:

```sh
cargo run --locked --release --offline --no-default-features \
  --manifest-path rust/Cargo.toml --example benchmark_all -- benchmarks/all_vs_all direct 48 /tmp/direct-scores.tsv
```

Run the paired comparison with:

```sh
cargo run --locked --release --offline --no-default-features \
  --manifest-path rust/Cargo.toml --example benchmark_all -- benchmarks/all_vs_all compare 48 /tmp/cached-scores.tsv
```

The Rust checkout used was branch `rust-implementation` at revision `a8d0bc2`. Release builds used the committed lockfile and cached dependencies. I built the benchmark with and without the optional `profile-phases` feature; I did not run the test suite.
