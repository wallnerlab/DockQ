//! Scoring-relevant mmCIF reader, including quoted and semicolon text values.
use super::normalize::{inferred_element, push_record, Record};
use crate::error::{DockError, Result};
use crate::types::{Atom, LoadOptions, Structure};
use std::collections::HashMap;

#[derive(Debug)]
struct Token {
    text: String,
    quoted: bool,
}
fn tokens(text: &str) -> Result<Vec<Token>> {
    let b = text.as_bytes();
    let mut i = 0;
    let mut out = Vec::new();
    while i < b.len() {
        if b[i].is_ascii_whitespace() {
            i += 1;
            continue;
        }
        if b[i] == b'#' {
            while i < b.len() && b[i] != b'\n' {
                i += 1;
            }
            continue;
        }
        if b[i] == b';' && (i == 0 || b[i - 1] == b'\n') {
            i += 1;
            let start = i;
            while i < b.len() && !(b[i] == b';' && b[i - 1] == b'\n') {
                i += 1;
            }
            if i == b.len() {
                return Err(DockError("unterminated mmCIF semicolon text".into()));
            }
            out.push(Token {
                text: text[start..i].trim_end_matches(['\r', '\n']).into(),
                quoted: true,
            });
            i += 1;
            continue;
        }
        if b[i] == b'\'' || b[i] == b'"' {
            let quote = b[i];
            i += 1;
            let start = i;
            while i < b.len()
                && !(b[i] == quote && (i + 1 == b.len() || b[i + 1].is_ascii_whitespace()))
            {
                i += 1;
            }
            if i == b.len() {
                return Err(DockError("unterminated mmCIF quoted value".into()));
            }
            out.push(Token {
                text: text[start..i].into(),
                quoted: true,
            });
            i += 1;
            continue;
        }
        let start = i;
        while i < b.len() && !b[i].is_ascii_whitespace() {
            i += 1;
        }
        out.push(Token {
            text: text[start..i].into(),
            quoted: false,
        });
    }
    Ok(out)
}
fn reserved(t: &Token) -> bool {
    !t.quoted
        && (t.text.starts_with('_')
            || t.text == "loop_"
            || t.text == "stop_"
            || t.text.starts_with("data_")
            || t.text.starts_with("save_"))
}
fn missing(s: &str) -> &str {
    if s == "." || s == "?" {
        ""
    } else {
        s
    }
}

pub fn parse_mmcif(text: &str, options: &LoadOptions) -> Result<Structure> {
    let ts = tokens(text)?;
    let mut i = 0;
    let mut structure = Structure::default();
    let mut models: Vec<String> = Vec::new();
    let mut found = false;
    while i < ts.len() {
        if ts[i].text != "loop_" || ts[i].quoted {
            i += 1;
            continue;
        }
        i += 1;
        let start = i;
        while i < ts.len() && !ts[i].quoted && ts[i].text.starts_with('_') {
            i += 1;
        }
        let tags = &ts[start..i];
        if tags.is_empty() {
            return Err(DockError("mmCIF loop has no tags".into()));
        }
        let atom_loop = tags.iter().any(|t| t.text.starts_with("_atom_site."));
        let indices: HashMap<&str, usize> = tags
            .iter()
            .enumerate()
            .map(|(n, t)| (t.text.as_str(), n))
            .collect();
        while i < ts.len() && !reserved(&ts[i]) {
            if i + tags.len() > ts.len() || ts[i..i + tags.len()].iter().any(reserved) {
                return Err(DockError("incomplete mmCIF loop row".into()));
            }
            let row = &ts[i..i + tags.len()];
            i += tags.len();
            if !atom_loop {
                continue;
            }
            found = true;
            let get =
                |key: &str| -> Option<&str> { indices.get(key).map(|&n| row[n].text.as_str()) };
            let req = |key: &str| -> Result<&str> {
                get(key).ok_or_else(|| DockError(format!("missing mmCIF column {key}")))
            };
            let model = get("_atom_site.pdbx_PDB_model_num").unwrap_or("1");
            let mi = match models.iter().position(|m| m == model) {
                Some(n) => n,
                None => {
                    models.push(model.into());
                    models.len() - 1
                }
            };
            if mi != options.model_index {
                continue;
            }
            let hetero = req("_atom_site.group_PDB")? == "HETATM";
            if hetero && !options.include_hetero {
                continue;
            }
            let chain = req(if options.auth_chains {
                "_atom_site.auth_asym_id"
            } else {
                "_atom_site.label_asym_id"
            })?;
            if options
                .chains
                .as_ref()
                .is_some_and(|cs| !cs.iter().any(|c| c == chain))
            {
                continue;
            }
            let name = req("_atom_site.label_atom_id")?;
            let element = get("_atom_site.type_symbol")
                .map(str::to_ascii_uppercase)
                .unwrap_or_else(|| inferred_element(name));
            if element == "H" {
                continue;
            }
            let resseq = get("_atom_site.auth_seq_id")
                .or_else(|| get("_atom_site.label_seq_id"))
                .ok_or_else(|| DockError("missing mmCIF residue numbering".into()))?;
            if resseq == "." {
                continue;
            }
            let number = resseq
                .parse()
                .map_err(|_| DockError(format!("invalid mmCIF residue number {resseq}")))?;
            let float = |key: &str| -> Result<f64> {
                req(key)?
                    .parse()
                    .map_err(|_| DockError(format!("invalid mmCIF numeric column {key}")))
            };
            let coord = [
                float("_atom_site.Cartn_x")? as f32 as f64,
                float("_atom_site.Cartn_y")? as f32 as f64,
                float("_atom_site.Cartn_z")? as f32 as f64,
            ];
            let occupancy = float("_atom_site.occupancy")?;
            let insertion = missing(get("_atom_site.pdbx_PDB_ins_code").unwrap_or("."));
            let altloc = missing(get("_atom_site.label_alt_id").unwrap_or("."));
            push_record(
                &mut structure,
                Record {
                    chain,
                    residue: req("_atom_site.label_comp_id")?,
                    number,
                    insertion,
                    hetero,
                    atom: Atom {
                        name: name.into(),
                        element,
                        coord,
                        occupancy,
                        altloc: altloc.into(),
                    },
                },
                options,
            )?;
        }
    }
    if !found || structure.chains.is_empty() {
        return Err(DockError(format!(
            "no selected atoms in mmCIF model {}",
            options.model_index
        )));
    }
    Ok(structure)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn quoted_values_and_semicolon_text() {
        let ts = tokens("data_x\n_note\n;hello\nworld\n;\nloop_\n_tag\n'CA' \"O3'\" # comment\n")
            .unwrap();
        assert!(ts.iter().any(|t| t.text == "hello\nworld"));
        assert!(ts.iter().any(|t| t.text == "O3'"));
        assert!(tokens("'unclosed").is_err());
    }
    #[test]
    fn ligand_chain_namespace() {
        let text = include_str!("../../../examples/1HHO_hem.cif");
        let s = parse_mmcif(
            text,
            &LoadOptions {
                include_hetero: true,
                auth_chains: false,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(s.chains.iter().any(|c| c.hetero.as_deref() == Some("HEM")));
        assert!(s.chains.iter().any(|c| c.hetero.is_none()));
    }
}
