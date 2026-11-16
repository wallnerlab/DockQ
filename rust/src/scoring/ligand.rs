//! Small-molecule interface scoring with graph-isomorphism symmetry correction.

use std::collections::HashSet;

use crate::alignment::align_chains;
use crate::contacts::{count_contacts, residue_distances};
use crate::error::{DockError, Result};
use crate::geometry::{apply, fit};
use crate::selection::{backbone_pairs, residue_coordinates};
use crate::types::{Chain, InterfaceScore, LigandMetrics, Metrics, ScoringOptions};

const BOND_TOLERANCE: f64 = 0.4;
const INTERFACE_THRESHOLD: f32 = 10.0;

// Covalent radii from DockQ's constants.py (Cordero et al., 2008).
const COVALENT_RADII: &[(&str, f64)] = &[
    ("H", 0.31),
    ("HE", 0.28),
    ("LI", 1.28),
    ("BE", 0.96),
    ("B", 0.84),
    ("C", 0.76),
    ("N", 0.71),
    ("O", 0.66),
    ("F", 0.57),
    ("NE", 0.58),
    ("NA", 1.66),
    ("MG", 1.41),
    ("AL", 1.21),
    ("SI", 1.11),
    ("P", 1.07),
    ("S", 1.05),
    ("CL", 1.02),
    ("AR", 1.06),
    ("K", 2.03),
    ("CA", 1.76),
    ("SC", 1.70),
    ("TI", 1.60),
    ("V", 1.53),
    ("CR", 1.39),
    ("MN", 1.50),
    ("FE", 1.42),
    ("CO", 1.38),
    ("NI", 1.24),
    ("CU", 1.32),
    ("ZN", 1.22),
    ("GA", 1.22),
    ("GE", 1.20),
    ("AS", 1.19),
    ("SE", 1.20),
    ("BR", 1.20),
    ("KR", 1.16),
    ("RB", 2.20),
    ("SR", 1.95),
    ("Y", 1.90),
    ("ZR", 1.75),
    ("NB", 1.64),
    ("MO", 1.54),
    ("TC", 1.47),
    ("RU", 1.46),
    ("RH", 1.42),
    ("PD", 1.39),
    ("AG", 1.45),
    ("CD", 1.44),
    ("IN", 1.42),
    ("SN", 1.39),
    ("SB", 1.39),
    ("TE", 1.38),
    ("I", 1.39),
    ("XE", 1.40),
    ("CS", 2.44),
    ("BA", 2.15),
    ("LA", 2.07),
    ("CE", 2.04),
    ("PR", 2.03),
    ("ND", 2.01),
    ("PM", 1.99),
    ("SM", 1.98),
    ("EU", 1.98),
    ("GD", 1.96),
    ("TB", 1.94),
    ("DY", 1.92),
    ("HO", 1.92),
    ("ER", 1.89),
    ("TM", 1.90),
    ("YB", 1.87),
    ("LU", 1.87),
    ("HF", 1.75),
    ("TA", 1.70),
    ("W", 1.62),
    ("RE", 1.51),
    ("OS", 1.44),
    ("IR", 1.41),
    ("PT", 1.36),
    ("AU", 1.36),
    ("HG", 1.32),
    ("TL", 1.45),
    ("PB", 1.46),
    ("BI", 1.48),
    ("PO", 1.40),
    ("AT", 1.50),
    ("RN", 1.50),
    ("FR", 2.60),
    ("RA", 2.21),
    ("AC", 2.15),
    ("TH", 2.06),
    ("PA", 2.00),
    ("U", 1.96),
    ("NP", 1.90),
    ("PU", 1.87),
    ("AM", 1.80),
    ("CM", 1.69),
    ("BK", 2.00),
    ("CF", 2.00),
    ("ES", 2.00),
    ("FM", 2.00),
    ("MD", 2.00),
    ("NO", 2.00),
    ("LR", 2.00),
    ("RF", 2.00),
    ("DB", 2.00),
    ("SG", 2.00),
    ("BH", 2.00),
    ("HS", 2.00),
    ("MT", 2.00),
    ("DS", 2.00),
    ("RG", 2.00),
    ("CN", 2.00),
    ("UUT", 2.00),
    ("UUQ", 2.00),
    ("UUP", 2.00),
    ("UUH", 2.00),
    ("UUS", 2.00),
    ("UUO", 2.00),
];

#[derive(Clone, Copy)]
struct LigandAtom<'a> {
    element: &'a str,
    coord: [f64; 3],
}

/// Score an interface with exactly one ligand chain using symmetry-corrected LRMSD.
///
/// A pair of ligand chains has no receptor frame and returns `Ok(None)`. If the
/// atom graphs have no isomorphism, this returns an error with both ligand IDs.
pub fn score_interface(
    model: [&Chain; 2],
    native: [&Chain; 2],
    options: &ScoringOptions,
) -> Result<Option<InterfaceScore>> {
    let model_ligand_index = match (model[0].hetero.is_some(), model[1].hetero.is_some()) {
        (true, false) => 0,
        (false, true) => 1,
        (true, true) => return Ok(None),
        (false, false) => {
            return Err(DockError(
                "ligand scoring requires one model ligand chain".into(),
            ))
        }
    };
    let receptor_index = 1 - model_ligand_index;
    let model_ligand = model[model_ligand_index];
    let model_receptor = model[receptor_index];
    let native_ligand = native[model_ligand_index];
    let native_receptor = native[receptor_index];

    let all_native_receptor: Vec<_> = (0..native_receptor.residues.len()).collect();
    let all_native_ligand: Vec<_> = (0..native_ligand.residues.len()).collect();
    let (receptor_coords, receptor_groups) =
        residue_coordinates(native_receptor, &all_native_receptor, false)?;
    let (ligand_coords, ligand_groups) =
        residue_coordinates(native_ligand, &all_native_ligand, false)?;
    let interface_distances = residue_distances(
        &receptor_coords,
        &ligand_coords,
        &receptor_groups,
        &ligand_groups,
    )?;
    if count_contacts(&interface_distances, INTERFACE_THRESHOLD) == 0 {
        return Ok(None);
    }

    let interface_residues: HashSet<usize> = (0..interface_distances.rows)
        .filter(|&row| {
            (0..interface_distances.cols).any(|col| {
                interface_distances.values[row * interface_distances.cols + col]
                    < INTERFACE_THRESHOLD * INTERFACE_THRESHOLD
            })
        })
        .collect();
    let mode = if options.no_align {
        crate::types::AlignmentMode::Numbering
    } else {
        crate::types::AlignmentMode::Sequence
    };
    let receptor_alignment = align_chains(model_receptor, native_receptor, mode)?;
    let matched_pairs = receptor_alignment.matched_pairs();
    let interface_pair_positions: Vec<_> = matched_pairs
        .iter()
        .enumerate()
        // The Python baseline indexes the aligned residue lists directly with
        // full native receptor indices, even when alignment gaps shift them.
        .filter_map(|(position, _)| interface_residues.contains(&position).then_some(position))
        .collect();
    if interface_pair_positions.is_empty() {
        return Ok(None);
    }
    let receptor_pairs = backbone_pairs(
        model_receptor,
        native_receptor,
        &matched_pairs,
        Some(&interface_pair_positions),
    )?;
    let transform = fit(&receptor_pairs.model, &receptor_pairs.native)?;

    let model_atoms = common_ligand_atoms(model_ligand, native_ligand);
    let native_atoms = common_ligand_atoms(native_ligand, model_ligand);
    if model_atoms.is_empty() || native_atoms.is_empty() || model_atoms.len() != native_atoms.len()
    {
        return Err(DockError(format!(
            "ligands {} and {} have no compatible atom-name set",
            model_ligand.id, native_ligand.id
        )));
    }
    let model_graph = create_graph(&model_atoms)?;
    let native_graph = create_graph(&native_atoms)?;
    let rotated_model: Vec<_> = model_atoms
        .iter()
        .map(|atom| apply(atom.coord, &transform))
        .collect();

    let mut best: Option<(f64, Vec<usize>)> = None;
    enumerate_isomorphisms(
        &model_graph,
        &native_graph,
        &rotated_model,
        &native_atoms,
        &mut best,
    );
    let (sum_squared, mapping) = best.ok_or_else(|| {
        DockError(format!(
            "ligand atom graphs for {} and {} are not isomorphic",
            model_ligand.id, native_ligand.id
        ))
    })?;
    let lrmsd = (sum_squared / model_atoms.len() as f64).sqrt();
    let atom_mapping = mapping
        .iter()
        .enumerate()
        .map(|(model_node, &native_node)| (model_node, native_node))
        .collect();
    let hetero = model_ligand.hetero.clone().unwrap_or_default();
    Ok(Some(InterfaceScore {
        dockq: (1.0 + 1.0 / (1.0 + (lrmsd / 8.5).powi(2))) / 3.0,
        metrics: Metrics::Ligand(LigandMetrics {
            lrmsd,
            hetero,
            atom_mapping,
        }),
    }))
}

/// Keep common atom names in each chain's input order.
///
/// The Python baseline filters coordinates but passes the unfiltered element
/// array to graph construction. This implementation filters elements together
/// with coordinates, so ligands with missing atom names remain well-defined.
fn common_ligand_atoms<'a>(ligand: &'a Chain, other: &Chain) -> Vec<LigandAtom<'a>> {
    let other_names: HashSet<&str> = other
        .residues
        .iter()
        .flat_map(|residue| residue.atoms.iter().map(|atom| atom.name.as_str()))
        .collect();
    let mut atoms = Vec::new();
    for residue in &ligand.residues {
        for atom in &residue.atoms {
            if other_names.contains(atom.name.as_str()) {
                atoms.push(LigandAtom {
                    element: atom.element.as_str(),
                    coord: atom.coord,
                });
            }
        }
    }
    atoms
}

fn create_graph(atoms: &[LigandAtom<'_>]) -> Result<Vec<Vec<bool>>> {
    let radii: Vec<_> = atoms
        .iter()
        .map(|atom| covalent_radius(atom.element))
        .collect::<Result<_>>()?;
    let float_coords: Vec<_> = atoms
        .iter()
        .map(|atom| atom.coord.map(|coordinate| coordinate as f32))
        .collect();
    if float_coords
        .iter()
        .flatten()
        .any(|value| !value.is_finite())
    {
        return Err(DockError("ligand coordinates must be finite in f32".into()));
    }

    let mut graph = vec![vec![false; atoms.len()]; atoms.len()];
    for i in 0..atoms.len() {
        for j in 0..atoms.len() {
            let distance = if i == j {
                0.0
            } else {
                let dx = float_coords[i][0] - float_coords[j][0];
                let dy = float_coords[i][1] - float_coords[j][1];
                let dz = float_coords[i][2] - float_coords[j][2];
                (dx * dx + dy * dy + dz * dz).sqrt()
            };
            let threshold = if i == j {
                1.0
            } else {
                radii[i] + radii[j] + BOND_TOLERANCE
            };
            if !distance.is_finite() {
                return Err(DockError("ligand atom distance is non-finite".into()));
            }
            if (distance as f64) < threshold {
                graph[i][j] = true;
            }
        }
    }
    Ok(graph)
}

fn covalent_radius(element: &str) -> Result<f64> {
    let normalized = element.trim().to_ascii_uppercase();
    COVALENT_RADII
        .iter()
        .find_map(|(symbol, radius)| (*symbol == normalized).then_some(*radius))
        .ok_or_else(|| {
            DockError(format!(
                "no covalent radius is defined for element {element:?}"
            ))
        })
}

fn enumerate_isomorphisms(
    model_graph: &[Vec<bool>],
    native_graph: &[Vec<bool>],
    model_coords: &[[f64; 3]],
    native_atoms: &[LigandAtom<'_>],
    best: &mut Option<(f64, Vec<usize>)>,
) {
    let count = model_graph.len();
    if count == 0
        || native_graph.len() != count
        || model_coords.len() != count
        || native_atoms.len() != count
    {
        return;
    }
    let model_degrees: Vec<_> = model_graph
        .iter()
        .map(|row| row.iter().filter(|&&x| x).count())
        .collect();
    let native_degrees: Vec<_> = native_graph
        .iter()
        .map(|row| row.iter().filter(|&&x| x).count())
        .collect();
    if model_degrees
        .iter()
        .any(|degree| !native_degrees.contains(degree))
    {
        return;
    }
    let mut mapping = vec![usize::MAX; count];
    let mut used = vec![false; count];
    enumerate_from(
        0,
        0.0,
        model_graph,
        native_graph,
        &model_degrees,
        &native_degrees,
        model_coords,
        native_atoms,
        &mut mapping,
        &mut used,
        best,
    );
}

#[allow(clippy::too_many_arguments)]
fn enumerate_from(
    model_node: usize,
    squared_error: f64,
    model_graph: &[Vec<bool>],
    native_graph: &[Vec<bool>],
    model_degrees: &[usize],
    native_degrees: &[usize],
    model_coords: &[[f64; 3]],
    native_atoms: &[LigandAtom<'_>],
    mapping: &mut [usize],
    used: &mut [bool],
    best: &mut Option<(f64, Vec<usize>)>,
) {
    if model_node == mapping.len() {
        if best
            .as_ref()
            .is_none_or(|(best_error, _)| squared_error < *best_error)
        {
            *best = Some((squared_error, mapping.to_vec()));
        }
        return;
    }

    for native_node in 0..mapping.len() {
        if used[native_node] || model_degrees[model_node] != native_degrees[native_node] {
            continue;
        }
        if (0..model_node).any(|previous| {
            model_graph[model_node][previous] != native_graph[native_node][mapping[previous]]
        }) {
            continue;
        }

        let native_coord = native_atoms[native_node].coord;
        let model_coord = model_coords[model_node];
        let dx = model_coord[0] - native_coord[0];
        let dy = model_coord[1] - native_coord[1];
        let dz = model_coord[2] - native_coord[2];
        let next_error = squared_error + dx * dx + dy * dy + dz * dz;
        if !next_error.is_finite()
            || best
                .as_ref()
                .is_some_and(|(best_error, _)| next_error > *best_error)
        {
            continue;
        }

        mapping[model_node] = native_node;
        used[native_node] = true;
        enumerate_from(
            model_node + 1,
            next_error,
            model_graph,
            native_graph,
            model_degrees,
            native_degrees,
            model_coords,
            native_atoms,
            mapping,
            used,
            best,
        );
        used[native_node] = false;
    }
}

#[cfg(test)]
mod tests {
    use super::{covalent_radius, score_interface, LigandAtom};
    use crate::error::Result;
    use crate::types::{Atom, Chain, Metrics, Residue, ScoringOptions};

    fn atom(name: &str, element: &str, coord: [f64; 3]) -> Atom {
        Atom {
            name: name.into(),
            element: element.into(),
            coord,
            occupancy: 1.0,
            altloc: String::new(),
        }
    }

    fn chain(id: &str, hetero: Option<&str>, residues: Vec<Residue>) -> Chain {
        Chain {
            id: id.into(),
            residues,
            sequence: "A".into(),
            hetero: hetero.map(str::to_string),
        }
    }

    fn residue(name: &str, atoms: Vec<Atom>, hetero: bool) -> Residue {
        Residue {
            name: name.into(),
            number: 1,
            insertion_code: String::new(),
            hetero,
            atoms,
        }
    }

    fn receptor(id: &str) -> Chain {
        chain(
            id,
            None,
            vec![residue(
                "ALA",
                vec![
                    atom("N", "N", [-1.0, 0.0, 0.0]),
                    atom("CA", "C", [0.0, 0.0, 0.0]),
                    atom("C", "C", [0.0, 1.0, 0.0]),
                    atom("O", "O", [-1.0, 1.0, 0.0]),
                ],
                false,
            )],
        )
    }

    #[test]
    fn both_ligand_chains_have_no_receptor_frame() {
        let ligand_a = chain("A", Some("LIG"), vec![]);
        let ligand_b = chain("B", Some("LIG"), vec![]);
        assert!(score_interface(
            [&ligand_a, &ligand_b],
            [&ligand_a, &ligand_b],
            &ScoringOptions::default()
        )
        .unwrap()
        .is_none());
    }

    #[test]
    fn unknown_element_has_no_radius() {
        assert!(covalent_radius("XX")
            .unwrap_err()
            .to_string()
            .contains("XX"));
    }

    #[test]
    fn symmetry_permutation_recovers_zero_lrmsd() -> Result<()> {
        let ligand_model = chain(
            "LM",
            Some("LIG"),
            vec![residue(
                "LIG",
                vec![
                    atom("C1", "C", [2.0, 2.0, 0.0]),
                    atom("EXTRA", "O", [50.0, 50.0, 50.0]),
                    atom("C2", "C", [3.0, 2.0, 0.0]),
                    atom("C3", "C", [2.0, 3.0, 0.0]),
                ],
                true,
            )],
        );
        let ligand_native = chain(
            "LN",
            Some("LIG"),
            vec![residue(
                "LIG",
                vec![
                    atom("C3", "C", [2.0, 3.0, 0.0]),
                    atom("C1", "C", [2.0, 2.0, 0.0]),
                    atom("C2", "C", [3.0, 2.0, 0.0]),
                ],
                true,
            )],
        );
        let model_receptor = receptor("RM");
        let native_receptor = receptor("RN");

        let score = score_interface(
            [&ligand_model, &model_receptor],
            [&ligand_native, &native_receptor],
            &ScoringOptions::default(),
        )?
        .unwrap();
        let Metrics::Ligand(metrics) = score.metrics else {
            panic!("expected ligand metrics");
        };
        assert!(metrics.lrmsd < 1e-10);
        assert!((score.dockq - 2.0 / 3.0).abs() < 1e-12);
        // EXTRA was filtered out, so indices are in the common-atom arrays.
        assert_eq!(metrics.atom_mapping, vec![(0, 1), (1, 2), (2, 0)]);
        Ok(())
    }

    #[test]
    fn graph_atom_fixture_keeps_element_for_each_selected_atom() {
        let atom = LigandAtom {
            element: "c",
            coord: [0.0; 3],
        };
        assert_eq!(covalent_radius(atom.element).unwrap(), 0.76);
    }

    #[test]
    fn common_atom_filter_preserves_each_input_order_and_filtered_indices() {
        use super::common_ligand_atoms;

        let model = chain(
            "M",
            Some("LIG"),
            vec![residue(
                "LIG",
                vec![
                    atom("C1", "C", [1.0, 0.0, 0.0]),
                    atom("EXTRA", "O", [99.0, 0.0, 0.0]),
                    atom("C2", "C", [2.0, 0.0, 0.0]),
                ],
                true,
            )],
        );
        let native = chain(
            "N",
            Some("LIG"),
            vec![residue(
                "LIG",
                vec![
                    atom("C2", "C", [20.0, 0.0, 0.0]),
                    atom("C1", "C", [10.0, 0.0, 0.0]),
                ],
                true,
            )],
        );
        let model_common = common_ligand_atoms(&model, &native);
        let native_common = common_ligand_atoms(&native, &model);
        assert_eq!(
            model_common.iter().map(|a| a.coord[0]).collect::<Vec<_>>(),
            [1.0, 2.0]
        );
        assert_eq!(
            native_common.iter().map(|a| a.coord[0]).collect::<Vec<_>>(),
            [20.0, 10.0]
        );
    }
}
