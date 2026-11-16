//! Shared scoring-relevant normalization for PDB and mmCIF readers.
use crate::error::{DockError, Result};
use crate::types::{Atom, Chain, LoadOptions, Residue, Structure};

pub fn residue_symbol(name: &str) -> char {
    if name.len() <= 2 {
        return name.chars().last().unwrap_or('X');
    }
    match name {
        "ALA" => 'A',
        "ARG" => 'R',
        "ASN" => 'N',
        "ASP" => 'D',
        "CYS" | "CME" => 'C',
        "GLN" => 'Q',
        "GLU" => 'E',
        "GLY" => 'G',
        "HIS" => 'H',
        "ILE" => 'I',
        "LEU" => 'L',
        "LYS" => 'K',
        "MET" | "MSE" => 'M',
        "PHE" => 'F',
        "PRO" => 'P',
        "SER" => 'S',
        "THR" => 'T',
        "TRP" => 'W',
        "TYR" => 'Y',
        "VAL" => 'V',
        "SEC" => 'U',
        "PYL" => 'O',
        "ASX" => 'B',
        "GLX" => 'Z',
        "XLE" => 'J',
        _ => 'X',
    }
}

pub struct Record<'a> {
    pub chain: &'a str,
    pub residue: &'a str,
    pub number: i32,
    pub insertion: &'a str,
    pub hetero: bool,
    pub atom: Atom,
}

pub fn push_record(s: &mut Structure, record: Record<'_>, opts: &LoadOptions) -> Result<()> {
    // Readers filter explicitly declared H before element inference, as Python does.
    if record.hetero && !opts.include_hetero {
        return Ok(());
    }
    if opts
        .chains
        .as_ref()
        .is_some_and(|cs| !cs.iter().any(|c| c == record.chain))
    {
        return Ok(());
    }
    if !record.atom.coord.iter().all(|x| x.is_finite()) {
        return Err(DockError("non-finite atom coordinate".into()));
    }
    let ci = match s.chains.iter().position(|c| c.id == record.chain) {
        Some(i) => i,
        None => {
            s.chains.push(Chain {
                id: record.chain.into(),
                residues: vec![],
                sequence: String::new(),
                hetero: None,
            });
            s.chains.len() - 1
        }
    };
    let chain = &mut s.chains[ci];
    let ri = match chain.residues.iter().position(|r| {
        r.number == record.number
            && r.insertion_code == record.insertion
            && r.hetero == record.hetero
    }) {
        Some(i) => i,
        None => {
            chain.residues.push(Residue {
                name: record.residue.into(),
                number: record.number,
                insertion_code: record.insertion.into(),
                hetero: record.hetero,
                atoms: vec![],
            });
            if record.hetero {
                chain.sequence = record.residue.into();
                chain.hetero = Some(record.residue.into());
            } else {
                chain.sequence.push(residue_symbol(record.residue));
            }
            chain.residues.len() - 1
        }
    };
    let r = &mut chain.residues[ri];
    // Biopython selects the first atom on equal occupancy and ignores duplicate
    // blank-altloc atoms. Distinct alternate labels select maximum occupancy.
    if let Some(i) = r.atoms.iter().position(|a| a.name == record.atom.name) {
        if r.atoms[i].altloc != record.atom.altloc && record.atom.occupancy > r.atoms[i].occupancy {
            r.atoms[i] = record.atom;
        }
    } else {
        r.atoms.push(record.atom);
    }
    Ok(())
}

pub fn inferred_element(atom_field: &str) -> String {
    let trimmed = atom_field.trim();
    if atom_field.starts_with(' ') || trimmed.starts_with(|c: char| c.is_ascii_digit()) {
        return trimmed
            .chars()
            .find(|c| c.is_ascii_alphabetic())
            .map(|c| c.to_ascii_uppercase().to_string())
            .unwrap_or_default();
    }
    let letters: String = trimmed
        .chars()
        .filter(|c| c.is_ascii_alphabetic())
        .take(2)
        .collect::<String>()
        .to_ascii_uppercase();
    match letters.as_str() {
        "FE" | "ZN" | "MG" | "MN" | "CL" | "BR" | "CA" | "NA" | "CU" | "NI" | "CO" | "SE" => {
            letters
        }
        _ => letters
            .chars()
            .next()
            .map(|c| c.to_string())
            .unwrap_or_default(),
    }
}
