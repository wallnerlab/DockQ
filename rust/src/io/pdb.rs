use super::normalize::{inferred_element, push_record, Record};
use crate::error::{DockError, Result};
use crate::types::{Atom, LoadOptions, Structure};

fn field(line: &str, start: usize, end: usize) -> &str {
    line.get(start..end.min(line.len())).unwrap_or("")
}
fn number<T: std::str::FromStr>(
    line: &str,
    start: usize,
    end: usize,
    what: &str,
    ln: usize,
) -> Result<T> {
    field(line, start, end)
        .trim()
        .parse()
        .map_err(|_| DockError(format!("invalid {what} at PDB line {ln}")))
}

pub fn parse_pdb(text: &str, options: &LoadOptions) -> Result<Structure> {
    let mut structure = Structure::default();
    let mut model = 0usize;
    let mut have_model = false;
    let mut in_model = false;
    let mut selected_model_seen = false;
    for (li, line) in text.lines().enumerate() {
        match field(line, 0, 6) {
            "MODEL " => {
                if have_model {
                    model += 1;
                }
                have_model = true;
                in_model = true;
                if model == options.model_index {
                    selected_model_seen = true;
                }
            }
            "ENDMDL" => {
                in_model = false;
            }
            "END   " | "CONECT" => break,
            "ATOM  " | "HETATM" => {
                if !have_model {
                    have_model = true;
                    in_model = true;
                }
                if !in_model || model != options.model_index {
                    continue;
                }
                selected_model_seen = true;
                let hetero = line.starts_with("HETATM");
                if hetero && !options.include_hetero {
                    continue;
                }
                let chain = field(line, 21, 22);
                if options
                    .chains
                    .as_ref()
                    .is_some_and(|c| !c.iter().any(|x| x == chain))
                {
                    continue;
                }
                let atom_field = field(line, 12, 16);
                let element_field = field(line, 76, 78).trim().to_ascii_uppercase();
                if element_field == "H" {
                    continue;
                }
                let element = if element_field.is_empty() {
                    inferred_element(atom_field)
                } else {
                    element_field
                };
                let coord = [
                    number::<f32>(line, 30, 38, "x", li + 1)? as f64,
                    number::<f32>(line, 38, 46, "y", li + 1)? as f64,
                    number::<f32>(line, 46, 54, "z", li + 1)? as f64,
                ];
                push_record(
                    &mut structure,
                    Record {
                        chain,
                        residue: field(line, 17, 20).trim(),
                        number: number(line, 22, 26, "residue number", li + 1)?,
                        insertion: field(line, 26, 27).trim(),
                        hetero,
                        atom: Atom {
                            name: atom_field.trim().into(),
                            element,
                            coord,
                            occupancy: field(line, 54, 60).trim().parse().unwrap_or(0.0),
                            altloc: field(line, 16, 17).trim().into(),
                        },
                    },
                    options,
                )?;
            }
            _ => {}
        }
    }
    if !selected_model_seen || structure.chains.is_empty() {
        return Err(DockError(format!(
            "no selected atoms in PDB model {}",
            options.model_index
        )));
    }
    Ok(structure)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn atom_line(
        serial: usize,
        alt: &str,
        number: i32,
        insert: &str,
        x: f64,
        occupancy: f64,
    ) -> String {
        format!("ATOM  {serial:5}   CA{alt}ALA A{number:4}{insert}   {x:8.3}{:8.3}{:8.3}{occupancy:6.2}{:6.2}           C",0.0,0.0,20.0)
    }
    #[test]
    fn alternates_duplicates_insertion_codes_and_models() {
        let text = format!(
            "MODEL        1\n{}\n{}\n{}\n{}\nENDMDL\nMODEL        2\n{}\nENDMDL\n",
            atom_line(1, "A", -1, " ", 1.0, 0.4),
            atom_line(2, "B", -1, " ", 2.0, 0.6),
            atom_line(3, "C", -1, " ", 3.0, 0.6),
            atom_line(4, " ", -1, "A", 4.0, 1.0),
            atom_line(5, " ", 1, " ", 8.0, 1.0)
        );
        let s = parse_pdb(&text, &LoadOptions::default()).unwrap();
        assert_eq!(s.chains[0].sequence, "AA");
        assert_eq!(s.chains[0].residues.len(), 2);
        assert_eq!(s.chains[0].residues[0].atoms.len(), 1);
        assert_eq!(s.chains[0].residues[0].atoms[0].coord[0], 2.0);
        assert_eq!(s.chains[0].residues[1].insertion_code, "A");
        let opts = LoadOptions {
            model_index: 1,
            ..Default::default()
        };
        let s = parse_pdb(&text, &opts).unwrap();
        assert_eq!(s.chains[0].residues.len(), 1);
        assert_eq!(s.chains[0].residues[0].atoms[0].coord[0], 8.0);
    }
    #[test]
    fn example_chain_order_and_precision() {
        let s = parse_pdb(
            include_str!("../../../examples/1A2K_r_l_b.pdb"),
            &LoadOptions::default(),
        )
        .unwrap();
        assert_eq!(
            s.chains.iter().map(|c| c.id.as_str()).collect::<Vec<_>>(),
            vec!["A", "B", "C"]
        );
        assert!(s.chains.iter().all(|c| !c.residues.is_empty()));
        for c in s.chains {
            for r in c.residues {
                for a in r.atoms {
                    assert_eq!(a.coord[0], a.coord[0] as f32 as f64);
                }
            }
        }
    }
    #[test]
    fn missing_model_and_invalid_text_fail() {
        assert!(parse_pdb("data_foo\n", &LoadOptions::default()).is_err());
        let opt = LoadOptions {
            model_index: 9,
            ..Default::default()
        };
        assert!(parse_pdb(include_str!("../../../examples/1A2K_r_l_b.pdb"), &opt).is_err());
    }
}
