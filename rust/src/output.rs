//! Text and JSON formatting for the DockQ command-line interface.
use crate::error::{DockError, Result};
use crate::types::{Evaluation, InterfaceResult, Metrics, Structure};
use serde_json::{Map, Value};

#[derive(Debug, Clone, Copy, Default)]
pub struct RenderOptions {
    pub short: bool,
    pub verbose: bool,
    pub capri_peptide: bool,
    pub small_molecule: bool,
}

pub fn render_text(
    model: &Structure,
    native: &Structure,
    e: &Evaluation,
    o: &RenderOptions,
) -> Result<String> {
    let mapping = mapping_string(model, native, &e.mapping)?;
    let score_name = if o.small_molecule {
        "DockQ-small_molecules"
    } else if o.capri_peptide {
        "DockQ-capri_peptide"
    } else {
        "DockQ"
    };
    let mut out = String::new();
    if o.short {
        out.push_str(&format!("Total {score_name} over {} native interfaces: {:.3} with {mapping} model:native mapping\n", e.interfaces.len(), e.global_dockq));
        for i in &e.interfaces {
            short_interface(&mut out, model, native, i)?;
        }
        return Ok(out);
    }
    out.push_str(&header(o.verbose, o.capri_peptide));
    out.push_str(&format!("Model  : {}\nNative : {}\n", model.id, native.id));
    out.push_str(&format!(
        "Total {score_name} over {} native interfaces: {:.3} with {mapping} model:native mapping\n",
        e.interfaces.len(),
        e.global_dockq
    ));
    for i in &e.interfaces {
        let n = pair_ids(native, i.native_chains, "native")?;
        let m = pair_ids(model, i.model_chains, "model")?;
        let het = match &i.score.metrics {
            Metrics::Ligand(x) => format!(" ({})", x.hetero),
            _ => String::new(),
        };
        out.push_str(&format!(
            "Native chains: {}, {}{het}\n\tModel chains: {}, {}\n",
            n[0], n[1], m[0], m[1]
        ));
        match &i.score.metrics {
            Metrics::Polymer(x) => out.push_str(&format!("\tDockQ: {:.3}\n\tiRMSD: {:.3}\n\tLRMSD: {:.3}\n\tfnat: {:.3}\n\tfnonnat: {:.3}\n\tF1: {:.3}\n\tclashes: {}\n",i.score.dockq,x.irmsd,x.lrmsd,x.fnat,x.fnonnat,x.f1,x.clashes)),
            Metrics::Ligand(x) => out.push_str(&format!("\tLRMSD: {:.3}\n",x.lrmsd)),
        }
    }
    Ok(out)
}

fn short_interface(
    out: &mut String,
    model: &Structure,
    native: &Structure,
    i: &InterfaceResult,
) -> Result<()> {
    let n = pair_ids(native, i.native_chains, "native")?;
    let m = pair_ids(model, i.model_chains, "model")?;
    match &i.score.metrics {
        Metrics::Polymer(x)=>out.push_str(&format!("DockQ {:.3} iRMSD {:.3} LRMSD {:.3} fnat {:.3} fnonnat {:.3} F1 {:.3} clashes {} mapping {}{}:{}{} {} {} {} -> {} {} {}\n",i.score.dockq,x.irmsd,x.lrmsd,x.fnat,x.fnonnat,x.f1,x.clashes,m[0],m[1],n[0],n[1],model.id,m[0],m[1],native.id,n[0],n[1])),
        Metrics::Ligand(x)=>out.push_str(&format!("LRMSD {:.3} mapping {}{}:{}{} ({}) {} {} {} -> {} {} {}\n",x.lrmsd,m[0],m[1],n[0],n[1],x.hetero,model.id,m[0],m[1],native.id,n[0],n[1])),
    }
    Ok(())
}

pub fn render_json(model: &Structure, native: &Structure, e: &Evaluation) -> Result<String> {
    serde_json::to_string(&json_value(model, native, e)?)
        .map_err(|x| DockError(format!("could not serialize results as JSON: {x}")))
}

pub fn json_value(model: &Structure, native: &Structure, e: &Evaluation) -> Result<Value> {
    let mut mapping = Map::new();
    for &(n, m) in &e.mapping {
        mapping.insert(
            chain_id(native, n, "native")?.into(),
            Value::String(chain_id(model, m, "model")?.into()),
        );
    }
    let mapping_value = Value::Object(mapping.clone());
    let mut results = Map::new();
    for i in &e.interfaces {
        let n = pair_ids(native, i.native_chains, "native")?;
        let m = pair_ids(model, i.model_chains, "model")?;
        let mut r = Map::new();
        r.insert("chain1".into(), m[0].clone().into());
        r.insert("chain2".into(), m[1].clone().into());
        r.insert("chain_map".into(), mapping_value.clone());
        match &i.score.metrics {
            Metrics::Polymer(x) => {
                r.insert("DockQ".into(), num(i.score.dockq));
                r.insert("iRMSD".into(), num(x.irmsd));
                r.insert("LRMSD".into(), num(x.lrmsd));
                r.insert("fnat".into(), num(x.fnat));
                r.insert("fnonnat".into(), num(x.fnonnat));
                r.insert("F1".into(), num(x.f1));
                r.insert("clashes".into(), (x.clashes as u64).into());
                r.insert("nat_correct".into(), (x.nat_correct as u64).into());
                r.insert("nat_total".into(), (x.nat_total as u64).into());
                r.insert("nonnat_count".into(), (x.nonnat_count as u64).into());
                r.insert("model_total".into(), (x.model_total as u64).into());
                r.insert("len1".into(), (x.len1 as u64).into());
                r.insert("len2".into(), (x.len2 as u64).into());
                r.insert("class1".into(), x.class1.clone().into());
                r.insert("class2".into(), x.class2.clone().into());
                r.insert("is_het".into(), false.into());
            }
            Metrics::Ligand(x) => {
                r.insert("DockQ".into(), num(i.score.dockq));
                r.insert("LRMSD".into(), num(x.lrmsd));
                let mut atom_mapping = Map::new();
                for &(model_atom, native_atom) in &x.atom_mapping {
                    atom_mapping.insert(model_atom.to_string(), Value::from(native_atom as u64));
                }
                r.insert("mapping".into(), Value::Object(atom_mapping));
                r.insert("is_het".into(), x.hetero.clone().into());
            }
        }
        results.insert(format!("{}{}", n[0], n[1]), Value::Object(r));
    }
    let mut root = Map::new();
    root.insert("model".into(), model.id.clone().into());
    root.insert("native".into(), native.id.clone().into());
    root.insert("best_dockq".into(), num(e.total_dockq));
    root.insert("best_result".into(), Value::Object(results));
    root.insert("GlobalDockQ".into(), num(e.global_dockq));
    root.insert("best_mapping".into(), mapping_value);
    root.insert(
        "best_mapping_str".into(),
        mapping_string(model, native, &e.mapping)?.into(),
    );
    Ok(Value::Object(root))
}

fn num(x: f64) -> Value {
    serde_json::Number::from_f64(x)
        .map(Value::Number)
        .unwrap_or(Value::Null)
}
fn mapping_string(
    model: &Structure,
    native: &Structure,
    mapping: &[(usize, usize)],
) -> Result<String> {
    let (mut m, mut n) = (String::new(), String::new());
    for &(ni, mi) in mapping {
        n.push_str(chain_id(native, ni, "native")?);
        m.push_str(chain_id(model, mi, "model")?);
    }
    Ok(format!("{m}:{n}"))
}
fn pair_ids(s: &Structure, pair: [usize; 2], label: &str) -> Result<[String; 2]> {
    Ok([
        chain_id(s, pair[0], label)?.into(),
        chain_id(s, pair[1], label)?.into(),
    ])
}
fn chain_id<'a>(s: &'a Structure, i: usize, label: &str) -> Result<&'a str> {
    s.chains
        .get(i)
        .map(|c| c.id.as_str())
        .ok_or_else(|| DockError(format!("{label} chain index {i} is out of range")))
}

fn header(verbose: bool, peptide: bool) -> String {
    let base="****************************************************************\n*                       DockQ                                  *\n*   Docking scoring for biomolecular models                    *\n*   DockQ score legend:                                        *\n*    0.00 <= DockQ <  0.23 - Incorrect                         *\n*    0.23 <= DockQ <  0.49 - Acceptable quality                *\n*    0.49 <= DockQ <  0.80 - Medium quality                    *\n*            DockQ >= 0.80 - High quality                      *";
    let reference="*   Ref: Mirabello and Wallner, 'DockQ v2: Improved automatic  *\n*   quality measure for protein multimers, nucleic acids       *\n*   and small molecules'                                       *\n*                                                              *\n*   For comments, please email: bjorn.wallner@.liu.se          *";
    let notice = if verbose {
        format!("*   For the record:                                            *\n*   Definition of contact <{} (Fnat)                           *\n*   Definition of interface <{} *\n****************************************************************",if peptide{"4A"}else{"5A"},if peptide{"8A CB (iRMSD)                    "}else{"10A all heavy atoms (iRMSD)      "})
    } else {
        "****************************************************************".into()
    };
    format!("{base}\n{reference}\n{notice}\n")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{Chain, InterfaceScore, PolymerMetrics};
    fn fixture() -> (Structure, Structure, Evaluation) {
        let s = |id: &str, ids: &[&str]| Structure {
            id: id.into(),
            chains: ids
                .iter()
                .map(|x| Chain {
                    id: (*x).into(),
                    residues: vec![],
                    sequence: String::new(),
                    hetero: None,
                })
                .collect(),
        };
        let x = PolymerMetrics {
            irmsd: 1.25,
            lrmsd: 2.5,
            fnat: 0.5,
            fnonnat: 0.1,
            f1: 0.6,
            clashes: 2,
            nat_correct: 3,
            nat_total: 6,
            nonnat_count: 1,
            model_total: 10,
            len1: 11,
            len2: 12,
            class1: "ligand".into(),
            class2: "receptor".into(),
        };
        (
            s("model.pdb", &["A", "B"]),
            s("native.pdb", &["X", "Y"]),
            Evaluation {
                mapping: vec![(0, 1), (1, 0)],
                interfaces: vec![InterfaceResult {
                    native_chains: [0, 1],
                    model_chains: [1, 0],
                    score: InterfaceScore {
                        dockq: 0.75,
                        metrics: Metrics::Polymer(x),
                    },
                }],
                total_dockq: 0.75,
                global_dockq: 0.75,
            },
        )
    }
    #[test]
    fn text_keeps_mapping_direction() {
        let (m, n, e) = fixture();
        let t = render_text(
            &m,
            &n,
            &e,
            &RenderOptions {
                short: true,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(t.contains("with BA:XY model:native mapping"));
        assert!(t.contains("mapping BA:XY"));
    }
    #[test]
    fn json_uses_legacy_keys() {
        let (m, n, e) = fixture();
        let v = json_value(&m, &n, &e).unwrap();
        assert_eq!(v["best_mapping_str"], "BA:XY");
        assert_eq!(v["best_mapping"]["X"], "B");
        assert_eq!(v["best_result"]["XY"]["iRMSD"], 1.25);
    }
}
