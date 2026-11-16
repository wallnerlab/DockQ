# Package interfaces

P01 contract, version 1. Types are in `src/types.rs`, errors in `src/error.rs`.
Owners may add local helpers. Shared contract changes go through the lead.

- P04 `io::pdb::parse_pdb(text: &str, options: &LoadOptions) -> Result<Structure>`.
  Own `io::normalize` shared atom/residue normalization; agree with P05 directly.
- P05 `io::mmcif::parse_mmcif(text: &str, options: &LoadOptions) -> Result<Structure>`.
- P06 `io::load::load_structure(path: impl AsRef<Path>, options: &LoadOptions) -> Result<Structure>`.
  gzip via flate2, try PDB then mmCIF; set returned id to path. No shell dependencies.
- P07 `alignment::align_chains(model: &Chain, native: &Chain, mode: AlignmentMode) -> Result<Alignment>`.
- P08 `selection::backbone_pairs(model: &Chain, native: &Chain, pairs: &[(usize,usize)], subset: Option<&[usize]>) -> Result<AtomPairs>`.
  Subset indices refer to positions within pairs, deduplicated in ascending order.
  `selection::residue_coordinates(chain: &Chain, indices: &[usize], cb_only: bool) -> Result<(Vec<Point>, Vec<usize>)>`.
- P09 `geometry::fit(source: &[Point], target: &[Point]) -> Result<RigidTransform>`;
  `geometry::apply(point: Point, transform: &RigidTransform) -> Point`;
  `geometry::rmsd(source: &[Point], target: &[Point]) -> Result<f64>`.
- P10 `contacts::residue_distances(a: &[Point], b: &[Point], groups_a: &[usize], groups_b: &[usize]) -> Result<SquaredDistances>`;
  `contacts::contact_stats(model: &SquaredDistances, native: &SquaredDistances, threshold: f32) -> Result<ContactStats>`;
  `contacts::count_contacts(distances: &SquaredDistances, threshold: f32) -> usize`.
  Public thresholds are linear Å, matrices squared Å², strict comparison.
- P11 `scoring::polymer::score_interface(model: [&Chain;2], native: [&Chain;2], options: &ScoringOptions) -> Result<Option<InterfaceScore>>`.
- P12 `scoring::ligand::score_interface(model: [&Chain;2], native: [&Chain;2], options: &ScoringOptions) -> Result<Option<InterfaceScore>>`.
- P13 `mapping::candidate_mappings(model: &Structure, native: &Structure, constraints: &MappingConstraints) -> Result<MappingIter>`;
  `MappingIter: Iterator<Item = ChainMapping>`. Preserve order, generate lazily.
  `mapping::count_mappings` returns the exact count by consuming this stream;
  counting can take factorial time, and the evaluation engine never calls it.
- P14 `engine::evaluate_mapping(model: &Structure, native: &Structure, mapping: &ChainMapping, options: &ScoringOptions) -> Result<Evaluation>`;
  `engine::evaluate(model: &Structure, native: &Structure, constraints: &MappingConstraints, options: &ScoringOptions, workers: usize, max_chunk: usize) -> Result<Evaluation>`.
- P15 `cli::run() -> Result<()>`; main reports error then exits 1. Output helpers
  consume structures and Evaluation; `--mapping` parsed here to MappingConstraints.

Shared semantics: input order matters; internal mapping is native→model; equal-score
search keeps first candidate; parser coordinates round through f32; geometry is f64;
distance operations use f32 to match compiled baseline. Contact counts exact; floating
metrics tolerance calibrated by P03/P16: absolute 1e-4 Å for RMSD, 1e-5 for scores
and other floating metrics (float32 baseline fits versus float64 Rust geometry).
No interface is `Ok(None)`, invalid/empty required atom selections are errors.
No library printing/exiting. Do not change Python source behavior to ease comparisons.

Minimal shared dependencies: serde/serde_json, flate2, nalgebra, rayon. P02 owns manifest.
Owners add module declarations for files they own; P02 coordinates `lib.rs` declarations.
