//! Cross-module scientific acceptance against the pinned Python/Cython baseline.
use dockq::io::load::load_structure;
use dockq::types::LoadOptions;
use serde_json::Value;
use std::path::PathBuf;

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .to_path_buf()
}
fn reference() -> Value {
    serde_json::from_str(&std::fs::read_to_string(root().join("compat/reference.json")).unwrap())
        .unwrap()
}

#[test]
fn alignment_columns_match_python_first_optimum() {
    use dockq::alignment::align_chains;
    use dockq::types::{AlignmentMode, Chain, Residue};
    let fixture = reference();
    for case in fixture["alignments"].as_array().unwrap() {
        let chain = |key: &str| {
            let sequence = case[key].as_str().unwrap().to_string();
            Chain {
                id: key.into(),
                residues: sequence
                    .chars()
                    .enumerate()
                    .map(|(i, _)| Residue {
                        name: "ALA".into(),
                        number: i as i32 - 2,
                        insertion_code: String::new(),
                        hetero: false,
                        atoms: vec![],
                    })
                    .collect(),
                sequence,
                hetero: None,
            }
        };
        let m = chain("model");
        let n = chain("native");
        let numbering = case["numbering"].as_bool().unwrap();
        let mode = if numbering {
            AlignmentMode::Numbering
        } else {
            AlignmentMode::Sequence
        };
        let a = align_chains(&m, &n, mode).unwrap();
        let mut model = String::new();
        let mut native = String::new();
        let symbol = |c: &Chain, i: Option<usize>| {
            i.map(|i| {
                if numbering {
                    char::from_u32((c.residues[i].number + 45) as u32).unwrap()
                } else {
                    c.sequence.chars().nth(i).unwrap()
                }
            })
            .unwrap_or('-')
        };
        for col in &a.columns {
            model.push(symbol(&m, col.model));
            native.push(symbol(&n, col.native));
        }
        assert_eq!(model, case["seqA"].as_str().unwrap(), "{case}");
        assert_eq!(native, case["seqB"].as_str().unwrap(), "{case}");
    }
}

#[test]
fn normalized_parser_results_match_python() {
    let reference = reference();
    for (file, expected) in reference["structures"].as_object().unwrap() {
        let ligand = file.contains("hem");
        let options = LoadOptions {
            include_hetero: ligand,
            auth_chains: !ligand,
            ..Default::default()
        };
        let s = load_structure(root().join("examples").join(file), &options).unwrap();
        let chains = expected.as_array().unwrap();
        assert_eq!(s.chains.len(), chains.len(), "{file} chain count");
        for (c, e) in s.chains.iter().zip(chains) {
            assert_eq!(c.id, e["id"].as_str().unwrap(), "{file} chain id");
            assert_eq!(
                c.sequence,
                e["sequence"].as_str().unwrap(),
                "{file} {} sequence",
                c.id
            );
            assert_eq!(
                serde_json::to_value(&c.hetero).unwrap(),
                e["hetero"],
                "{file} hetero"
            );
            let residues = e["residues"].as_array().unwrap();
            assert_eq!(
                c.residues.len(),
                residues.len(),
                "{file} {} residue count",
                c.id
            );
            for (r, er) in c.residues.iter().zip(residues) {
                assert_eq!(r.name, er["name"].as_str().unwrap(), "{file} residue name");
                assert_eq!(r.number, er["number"].as_i64().unwrap() as i32);
                assert_eq!(r.insertion_code, er["insertion_code"].as_str().unwrap());
                assert_eq!(r.hetero, er["hetero"].as_bool().unwrap());
                assert_eq!(
                    serde_json::to_value(r.atoms.iter().map(|a| &a.name).collect::<Vec<_>>())
                        .unwrap(),
                    er["atom_names"],
                    "{file} {} residue {} atom order",
                    c.id,
                    r.number
                );
                for (v, ev) in r.atoms[0]
                    .coord
                    .iter()
                    .zip(er["first_coord"].as_array().unwrap())
                {
                    assert!(
                        (v - ev.as_f64().unwrap()).abs() < 1e-12,
                        "{file} coordinate"
                    );
                }
            }
        }
    }
}

#[test]
fn explicit_mapping_scores_match_python() {
    use dockq::engine::evaluate_mapping;
    use dockq::types::{Metrics, ScoringOptions};
    let r = reference();
    for case in r["cases"].as_array().unwrap() {
        let name = case["name"].as_str().unwrap();
        let small = case["options"]["small_molecule"].as_bool().unwrap_or(false);
        let opt = LoadOptions {
            include_hetero: small,
            auth_chains: !small,
            ..Default::default()
        };
        let m = load_structure(
            root()
                .join("examples")
                .join(case["model"].as_str().unwrap()),
            &opt,
        )
        .unwrap();
        let n = load_structure(
            root()
                .join("examples")
                .join(case["native"].as_str().unwrap()),
            &opt,
        )
        .unwrap();
        let mapping = case["mapping"]
            .as_array()
            .unwrap()
            .iter()
            .map(|p| {
                (
                    n.chains
                        .iter()
                        .position(|c| c.id == p[0].as_str().unwrap())
                        .unwrap(),
                    m.chains
                        .iter()
                        .position(|c| c.id == p[1].as_str().unwrap())
                        .unwrap(),
                )
            })
            .collect();
        let options = ScoringOptions {
            peptide: case["options"]["capri_peptide"].as_bool().unwrap_or(false),
            no_align: case["options"]["no_align"].as_bool().unwrap_or(false),
        };
        let actual =
            evaluate_mapping(&m, &n, &mapping, &options).unwrap_or_else(|e| panic!("{name}: {e}"));
        let expected = case["results"].as_object().unwrap();
        assert_eq!(
            actual.interfaces.len(),
            expected.len(),
            "{name} interface count"
        );
        let close = |field: &str, a: f64, e: f64| {
            let tolerance = if field == "LRMSD" || field == "iRMSD" {
                1e-4
            } else {
                1e-5
            };
            assert!(
                (a - e).abs() <= tolerance,
                "{name} {field}: {a} vs {e}, difference {}",
                (a - e).abs()
            );
        };
        close(
            "total",
            actual.total_dockq,
            case["total_dockq"].as_f64().unwrap(),
        );
        for interface in &actual.interfaces {
            let key = format!(
                "{}{}",
                n.chains[interface.native_chains[0]].id, n.chains[interface.native_chains[1]].id
            );
            let er = &expected[&key];
            close(
                "DockQ",
                interface.score.dockq,
                er["DockQ"].as_f64().unwrap(),
            );
            match &interface.score.metrics {
                Metrics::Polymer(p) => {
                    for (field, v) in [
                        ("iRMSD", p.irmsd),
                        ("LRMSD", p.lrmsd),
                        ("fnat", p.fnat),
                        ("fnonnat", p.fnonnat),
                        ("F1", p.f1),
                    ] {
                        close(field, v, er[field].as_f64().unwrap());
                    }
                    for (field, v) in [
                        ("nat_correct", p.nat_correct),
                        ("nat_total", p.nat_total),
                        ("nonnat_count", p.nonnat_count),
                        ("model_total", p.model_total),
                        ("clashes", p.clashes),
                        ("len1", p.len1),
                        ("len2", p.len2),
                    ] {
                        assert_eq!(
                            v,
                            er[field].as_u64().unwrap() as usize,
                            "{name} {key} {field}"
                        );
                    }
                    assert_eq!(p.class1, er["class1"].as_str().unwrap());
                    assert_eq!(p.class2, er["class2"].as_str().unwrap());
                }
                Metrics::Ligand(l) => {
                    close("LRMSD", l.lrmsd, er["LRMSD"].as_f64().unwrap());
                    assert_eq!(l.hetero, er["is_het"].as_str().unwrap());
                }
            }
        }
    }
}
