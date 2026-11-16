//! Atom correspondence and coordinate selection.

use crate::error::{DockError, Result};
use crate::types::{AtomPairs, Chain, Point};

/// Backbone atom order used by the Python implementation (`constants.py`).
/// Keeping this order preserves the coordinate ordering expected downstream.
const BACKBONE_ATOMS: &[&str] = &[
    "CA", "C", "N", "O", "P", "OP1", "OP2", "O2'", "O3'", "O4'", "O5'", "C1'", "C2'", "C3'", "C4'",
    "C5'",
];

/// Compact, ordered backbone coordinates cached once per chain.
#[derive(Debug, Clone)]
pub(crate) struct CachedBackboneResidue {
    coordinates: [[f64; 3]; 16],
    present: u16,
}

pub(crate) fn cache_backbone(chain: &Chain) -> Vec<CachedBackboneResidue> {
    chain
        .residues
        .iter()
        .map(|residue| {
            let mut cached = CachedBackboneResidue {
                coordinates: [[0.0; 3]; 16],
                present: 0,
            };
            for (index, name) in BACKBONE_ATOMS.iter().enumerate() {
                if let Some(atom) = residue.atoms.iter().find(|atom| atom.name == **name) {
                    cached.coordinates[index] = atom.coord;
                    cached.present |= 1 << index;
                }
            }
            cached
        })
        .collect()
}

/// Select corresponding backbone coordinates for matched residue index pairs.
///
/// `subset` contains positions within `pairs`; duplicates are removed and
/// positions are visited in ascending order. `None` selects every pair.
/// Missing atom names are skipped independently for each residue pair.
pub fn backbone_pairs(
    model: &Chain,
    native: &Chain,
    pairs: &[(usize, usize)],
    subset: Option<&[usize]>,
) -> Result<AtomPairs> {
    let positions = selected_positions(pairs.len(), subset, "matched residue pair")?;
    let mut selected = AtomPairs::default();

    for position in positions {
        let (model_index, native_index) = pairs[position];
        let model_residue = model.residues.get(model_index).ok_or_else(|| {
            DockError(format!(
                "model residue index {model_index} is out of range for chain {}",
                model.id
            ))
        })?;
        let native_residue = native.residues.get(native_index).ok_or_else(|| {
            DockError(format!(
                "native residue index {native_index} is out of range for chain {}",
                native.id
            ))
        })?;

        for name in BACKBONE_ATOMS {
            if let (Some(model_atom), Some(native_atom)) = (
                find_atom(&model_residue.atoms, name),
                find_atom(&native_residue.atoms, name),
            ) {
                selected.model.push(model_atom.coord);
                selected.native.push(native_atom.coord);
            }
        }
    }

    if selected.model.is_empty() {
        return Err(DockError("backbone atom selection is empty".into()));
    }
    Ok(selected)
}

pub(crate) fn backbone_pairs_cached(
    model: &[CachedBackboneResidue],
    native: &[CachedBackboneResidue],
    model_id: &str,
    native_id: &str,
    pairs: &[(usize, usize)],
    subset: Option<&[usize]>,
) -> Result<AtomPairs> {
    let positions = selected_positions(pairs.len(), subset, "matched residue pair")?;
    let mut selected = AtomPairs::default();

    for position in positions {
        let (model_index, native_index) = pairs[position];
        let model_residue = model.get(model_index).ok_or_else(|| {
            DockError(format!(
                "model residue index {model_index} is out of range for chain {model_id}"
            ))
        })?;
        let native_residue = native.get(native_index).ok_or_else(|| {
            DockError(format!(
                "native residue index {native_index} is out of range for chain {native_id}"
            ))
        })?;
        let shared = model_residue.present & native_residue.present;
        for atom_index in 0..BACKBONE_ATOMS.len() {
            if shared & (1 << atom_index) != 0 {
                selected.model.push(model_residue.coordinates[atom_index]);
                selected.native.push(native_residue.coordinates[atom_index]);
            }
        }
    }

    if selected.model.is_empty() {
        return Err(DockError("backbone atom selection is empty".into()));
    }
    Ok(selected)
}

/// Flatten coordinates from selected residues and return each residue's atom count.
///
/// In all-atom mode, atoms retain their input order within each residue. In
/// CB-only mode, each residue contributes CB, falling back to CA when CB is
/// absent.
pub fn residue_coordinates(
    chain: &Chain,
    indices: &[usize],
    cb_only: bool,
) -> Result<(Vec<Point>, Vec<usize>)> {
    let mut coordinates = Vec::new();
    let mut groups = Vec::new();

    for &residue_index in indices {
        let residue = chain.residues.get(residue_index).ok_or_else(|| {
            DockError(format!(
                "residue index {residue_index} is out of range for chain {}",
                chain.id
            ))
        })?;

        if cb_only {
            let atom = find_atom(&residue.atoms, "CB")
                .or_else(|| find_atom(&residue.atoms, "CA"))
                .ok_or_else(|| {
                    DockError(format!(
                        "residue {} ({}) in chain {} has neither CB nor CA",
                        residue_index, residue.name, chain.id
                    ))
                })?;
            coordinates.push(atom.coord);
            groups.push(1);
        } else {
            for atom in &residue.atoms {
                coordinates.push(atom.coord);
            }
            groups.push(residue.atoms.len());
        }
    }

    if coordinates.is_empty() {
        return Err(DockError("residue coordinate selection is empty".into()));
    }
    Ok((coordinates, groups))
}

fn selected_positions(len: usize, subset: Option<&[usize]>, what: &str) -> Result<Vec<usize>> {
    let mut positions: Vec<usize> = match subset {
        None => (0..len).collect(),
        Some(positions) => positions.to_vec(),
    };
    positions.sort_unstable();
    positions.dedup();
    if let Some(&invalid) = positions.iter().find(|&&position| position >= len) {
        return Err(DockError(format!(
            "{what} subset position {invalid} is out of range (length {len})"
        )));
    }
    Ok(positions)
}

fn find_atom<'a>(atoms: &'a [crate::types::Atom], name: &str) -> Option<&'a crate::types::Atom> {
    atoms.iter().find(|atom| atom.name == name)
}

#[cfg(test)]
mod tests {
    use super::{backbone_pairs, residue_coordinates};
    use crate::types::{Atom, Chain, Residue};

    fn atom(name: &str, x: f64) -> Atom {
        Atom {
            name: name.into(),
            element: name.chars().next().unwrap_or('C').to_string(),
            coord: [x, 0.0, 0.0],
            occupancy: 1.0,
            altloc: String::new(),
        }
    }

    fn residue(name: &str, atoms: Vec<Atom>) -> Residue {
        Residue {
            name: name.into(),
            number: 1,
            insertion_code: String::new(),
            hetero: false,
            atoms,
        }
    }

    fn chain(id: &str, residues: Vec<Residue>) -> Chain {
        Chain {
            id: id.into(),
            residues,
            sequence: String::new(),
            hetero: None,
        }
    }

    #[test]
    fn backbone_pairs_skip_missing_atoms_and_sort_deduplicated_subset() {
        let model = chain(
            "M",
            vec![
                residue("ALA", vec![atom("N", 10.0), atom("CA", 11.0)]),
                residue("GLY", vec![atom("O", 20.0), atom("C", 21.0)]),
            ],
        );
        let native = chain(
            "N",
            vec![
                residue("ALA", vec![atom("CA", 1.0), atom("N", 2.0)]),
                residue("GLY", vec![atom("C", 3.0)]),
            ],
        );

        // The subset addresses pair positions, not residue indices.
        let pairs = [(0, 0), (1, 1)];
        let result = backbone_pairs(&model, &native, &pairs, Some(&[1, 0, 1])).unwrap();
        assert_eq!(
            result.model,
            vec![[11.0, 0.0, 0.0], [10.0, 0.0, 0.0], [21.0, 0.0, 0.0]]
        );
        assert_eq!(
            result.native,
            vec![[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]]
        );
    }

    #[test]
    fn backbone_pairs_include_nucleic_backbone_names_in_contract_order() {
        let atoms = vec![
            atom("C5'", 5.0),
            atom("OP2", 2.0),
            atom("P", 1.0),
            atom("O3'", 3.0),
        ];
        let model = chain("M", vec![residue("A", atoms.clone())]);
        let native = chain("N", vec![residue("A", atoms)]);
        let result = backbone_pairs(&model, &native, &[(0, 0)], None).unwrap();
        assert_eq!(
            result.model,
            vec![
                [1.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [3.0, 0.0, 0.0],
                [5.0, 0.0, 0.0]
            ]
        );
    }

    #[test]
    fn residue_coordinates_preserve_atom_order_and_return_per_residue_counts() {
        let chain = chain(
            "A",
            vec![
                residue("ALA", vec![atom("C", 1.0), atom("N", 2.0)]),
                residue("GLY", vec![atom("CA", 3.0)]),
            ],
        );
        let (coords, groups) = residue_coordinates(&chain, &[1, 0], false).unwrap();
        assert_eq!(
            coords,
            vec![[3.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]
        );
        assert_eq!(groups, vec![1, 2]);
    }

    #[test]
    fn cb_only_prefers_cb_and_falls_back_to_ca() {
        let chain = chain(
            "A",
            vec![
                residue("ALA", vec![atom("CA", 1.0), atom("CB", 2.0)]),
                residue("GLY", vec![atom("CA", 3.0)]),
            ],
        );
        let (coords, groups) = residue_coordinates(&chain, &[0, 1], true).unwrap();
        assert_eq!(coords, vec![[2.0, 0.0, 0.0], [3.0, 0.0, 0.0]]);
        assert_eq!(groups, vec![1, 1]);
    }

    #[test]
    fn cb_only_errors_when_neither_cb_nor_ca_exists() {
        let chain = chain("A", vec![residue("LIG", vec![atom("FE", 1.0)])]);
        let error = residue_coordinates(&chain, &[0], true).unwrap_err();
        assert!(error.to_string().contains("neither CB nor CA"));
    }
}
