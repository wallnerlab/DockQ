//! Shared immutable-by-convention data contracts. Input order is scientific state.
use serde::{Deserialize, Serialize};

pub type Point = [f64; 3];

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct Atom {
    pub name: String,
    pub element: String,
    /// Coordinates are rounded to f32 by parsers, then stored as f64.
    pub coord: Point,
    pub occupancy: f64,
    pub altloc: String,
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct Residue {
    pub name: String,
    pub number: i32,
    pub insertion_code: String,
    pub hetero: bool,
    pub atoms: Vec<Atom>,
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct Chain {
    pub id: String,
    pub residues: Vec<Residue>,
    pub sequence: String,
    /// None denotes polymer; Some(residue name) denotes a ligand chain.
    pub hetero: Option<String>,
}
#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub struct Structure {
    pub id: String,
    pub chains: Vec<Chain>,
}
#[derive(Debug, Clone)]
pub struct LoadOptions {
    pub chains: Option<Vec<String>>,
    pub model_index: usize,
    pub include_hetero: bool,
    pub auth_chains: bool,
}
impl Default for LoadOptions {
    fn default() -> Self {
        Self {
            chains: None,
            model_index: 0,
            include_hetero: false,
            auth_chains: true,
        }
    }
}
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum AlignmentMode {
    #[default]
    Sequence,
    Numbering,
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct AlignmentColumn {
    pub model: Option<usize>,
    pub native: Option<usize>,
    pub matched: bool,
}
#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub struct Alignment {
    pub columns: Vec<AlignmentColumn>,
    pub mismatches: usize,
}
impl Alignment {
    pub fn matched_pairs(&self) -> Vec<(usize, usize)> {
        self.columns
            .iter()
            .filter(|c| c.matched)
            .filter_map(|c| Some((c.model?, c.native?)))
            .collect()
    }
}
#[derive(Debug, Clone, Default)]
pub struct AtomPairs {
    pub model: Vec<Point>,
    pub native: Vec<Point>,
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct RigidTransform {
    /// Column convention: target = rotation * source + translation.
    pub rotation: [[f64; 3]; 3],
    pub translation: Point,
    pub rmsd: f64,
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct SquaredDistances {
    pub rows: usize,
    pub cols: usize,
    pub values: Vec<f32>,
}
#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq, Eq)]
pub struct ContactStats {
    pub shared: usize,
    pub non_native: usize,
    pub native: usize,
    pub model: usize,
}
/// Ordered (native chain index, model chain index) pairs.
pub type ChainMapping = Vec<(usize, usize)>;
#[derive(Debug, Clone, Default)]
pub struct MappingConstraints {
    pub fixed: ChainMapping,
    pub model_chains: Option<Vec<usize>>,
    pub native_chains: Option<Vec<usize>>,
    pub allowed_mismatches: usize,
}
#[derive(Debug, Clone, Copy, Default)]
pub struct ScoringOptions {
    pub peptide: bool,
    pub no_align: bool,
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct PolymerMetrics {
    pub irmsd: f64,
    pub lrmsd: f64,
    pub fnat: f64,
    pub fnonnat: f64,
    pub f1: f64,
    pub clashes: usize,
    pub nat_correct: usize,
    pub nat_total: usize,
    pub nonnat_count: usize,
    pub model_total: usize,
    pub len1: usize,
    pub len2: usize,
    pub class1: String,
    pub class2: String,
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct LigandMetrics {
    pub lrmsd: f64,
    pub hetero: String,
    pub atom_mapping: Vec<(usize, usize)>,
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum Metrics {
    Polymer(PolymerMetrics),
    Ligand(LigandMetrics),
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct InterfaceScore {
    pub dockq: f64,
    pub metrics: Metrics,
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct InterfaceResult {
    pub native_chains: [usize; 2],
    pub model_chains: [usize; 2],
    pub score: InterfaceScore,
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct Evaluation {
    pub mapping: ChainMapping,
    pub interfaces: Vec<InterfaceResult>,
    pub total_dockq: f64,
    pub global_dockq: f64,
}
