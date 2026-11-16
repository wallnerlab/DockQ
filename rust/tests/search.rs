//! Mapping enumeration and deterministic engine acceptance tests.

use dockq::engine::{evaluate, evaluate_mapping};
use dockq::io::load::load_structure;
use dockq::mapping::candidate_mappings;
use dockq::types::{Chain, LoadOptions, MappingConstraints, Residue, ScoringOptions, Structure};
use serde_json::Value;
use std::path::PathBuf;

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .to_path_buf()
}

fn example(file: &str) -> Structure {
    load_structure(root().join("examples").join(file), &LoadOptions::default()).unwrap()
}

fn from_case_chains(field: &Value, prefix: &str) -> Vec<Chain> {
    field
        .as_array()
        .unwrap()
        .iter()
        .enumerate()
        .map(|(index, value)| {
            let (sequence, hetero) = if let Some(sequence) = value.as_str() {
                (sequence, None)
            } else {
                (
                    value["sequence"].as_str().unwrap(),
                    value["hetero"].as_str(),
                )
            };
            Chain {
                id: format!("{prefix}{index}"),
                sequence: sequence.into(),
                hetero: hetero.map(str::to_owned),
                residues: sequence
                    .chars()
                    .enumerate()
                    .map(|(i, _)| Residue {
                        name: "ALA".into(),
                        number: i as i32 + 1,
                        insertion_code: String::new(),
                        hetero: false,
                        atoms: Vec::new(),
                    })
                    .collect(),
            }
        })
        .collect()
}

#[test]
fn mapping_order_matches_python_cases() {
    let fixture: Value =
        serde_json::from_str(include_str!("../../compat/mapping_cases.json")).unwrap();
    for (case_index, case) in fixture["cases"].as_array().unwrap().iter().enumerate() {
        let model_chains = from_case_chains(&case["model"], &format!("M{case_index}"));
        let native_chains = from_case_chains(&case["native"], &format!("N{case_index}"));
        let model = Structure {
            id: "model".into(),
            chains: model_chains,
        };
        let native = Structure {
            id: "native".into(),
            chains: native_chains,
        };
        let all_model: Vec<_> = (0..model.chains.len()).collect();
        let all_native: Vec<_> = (0..native.chains.len()).collect();
        let constraints = MappingConstraints {
            fixed: case["fixed"]
                .as_array()
                .into_iter()
                .flatten()
                .map(|pair| {
                    (
                        pair[0].as_u64().unwrap() as usize,
                        pair[1].as_u64().unwrap() as usize,
                    )
                })
                .collect(),
            model_chains: case["model_selected"]
                .as_array()
                .map(|values| {
                    values
                        .iter()
                        .map(|value| value.as_u64().unwrap() as usize)
                        .collect()
                })
                .or(Some(all_model)),
            native_chains: case["native_selected"]
                .as_array()
                .map(|values| {
                    values
                        .iter()
                        .map(|value| value.as_u64().unwrap() as usize)
                        .collect()
                })
                .or(Some(all_native)),
            allowed_mismatches: case["allowed_mismatches"].as_u64().unwrap_or(0) as usize,
        };
        let result = candidate_mappings(&model, &native, &constraints);
        if !case["expected"]["ok"].as_bool().unwrap() {
            assert!(
                result.is_err(),
                "{} should have no compatible chain",
                case["name"]
            );
            continue;
        }
        let actual: Vec<_> = result.unwrap().collect();
        let expected: Vec<Vec<(usize, usize)>> = case["expected"]["mappings"]
            .as_array()
            .unwrap()
            .iter()
            .map(|mapping| {
                mapping
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|pair| {
                        (
                            pair[0].as_u64().unwrap() as usize,
                            pair[1].as_u64().unwrap() as usize,
                        )
                    })
                    .collect()
            })
            .collect();
        assert_eq!(actual, expected, "{}", case["name"]);
    }
}

#[test]
fn large_mapping_search_yields_prefix_without_materializing_all_candidates() {
    let make_chain = |index: usize| Chain {
        id: format!("{index}"),
        residues: Vec::new(),
        sequence: "A".into(),
        hetero: Some("LIG".into()),
    };
    let model = Structure {
        id: String::new(),
        chains: (0..12).map(make_chain).collect(),
    };
    let native = Structure {
        id: String::new(),
        chains: (0..8).map(make_chain).collect(),
    };
    let first_five: Vec<_> = candidate_mappings(&model, &native, &MappingConstraints::default())
        .unwrap()
        .take(5)
        .collect();
    assert_eq!(first_five.len(), 5);
    assert_eq!(
        first_five[0],
        (0..8).map(|index| (index, index)).collect::<Vec<_>>()
    );
    assert_eq!(
        first_five[1],
        vec![
            (0, 0),
            (1, 1),
            (2, 2),
            (3, 3),
            (4, 4),
            (5, 5),
            (6, 6),
            (7, 8)
        ]
    );
}

fn assert_serial_parallel_equal(model: &Structure, native: &Structure) {
    let constraints = MappingConstraints::default();
    let options = ScoringOptions::default();
    let serial = evaluate(model, native, &constraints, &options, 1, 1).unwrap();
    let parallel = evaluate(model, native, &constraints, &options, 4, 2).unwrap();
    assert_eq!(serial, parallel);
}

#[test]
fn serial_and_parallel_results_match_on_1a2k_and_dimer() {
    let model = example("1A2K_r_l_b.model.pdb");
    let native = example("1A2K_r_l_b.pdb");
    assert_serial_parallel_equal(&model, &native);

    let model = example("dimer_dimer.model.pdb");
    let native = example("dimer_dimer.pdb");
    assert_serial_parallel_equal(&model, &native);
}

#[test]
fn serial_and_parallel_results_match_for_reverse_mapping() {
    let full_model = example("1A2K_r_l_b.model.pdb");
    let native = example("1A2K_r_l_b.pdb");
    let model = Structure {
        id: full_model.id,
        chains: vec![full_model.chains[0].clone(), full_model.chains[2].clone()],
    };
    assert_eq!(
        candidate_mappings(&model, &native, &MappingConstraints::default())
            .unwrap()
            .next()
            .unwrap(),
        vec![(0, 0), (2, 1)]
    );
    assert_serial_parallel_equal(&model, &native);
}

#[test]
fn equal_score_mapping_tie_keeps_first_candidate_across_worker_counts() {
    let model = example("1A2K_r_l_b.model.pdb");
    let native = example("1A2K_r_l_b.pdb");
    let model_chain = model.chains[0].clone();
    let native_chain = native.chains[0].clone();
    let tied_model = Structure {
        id: "tied-model".into(),
        chains: vec![
            model_chain.clone(),
            Chain {
                id: "B".into(),
                ..model_chain
            },
        ],
    };
    let tied_native = Structure {
        id: "tied-native".into(),
        chains: vec![
            native_chain.clone(),
            Chain {
                id: "B".into(),
                ..native_chain
            },
        ],
    };
    let first = vec![(0, 0), (1, 1)];
    let swapped = vec![(0, 1), (1, 0)];
    let score_first = evaluate_mapping(
        &tied_model,
        &tied_native,
        &first,
        &ScoringOptions::default(),
    )
    .unwrap();
    let score_swapped = evaluate_mapping(
        &tied_model,
        &tied_native,
        &swapped,
        &ScoringOptions::default(),
    )
    .unwrap();
    assert_eq!(
        score_first.interfaces[0].score,
        score_swapped.interfaces[0].score
    );
    let serial = evaluate(
        &tied_model,
        &tied_native,
        &MappingConstraints::default(),
        &ScoringOptions::default(),
        1,
        1,
    )
    .unwrap();
    let parallel = evaluate(
        &tied_model,
        &tied_native,
        &MappingConstraints::default(),
        &ScoringOptions::default(),
        4,
        1,
    )
    .unwrap();
    assert_eq!(serial.mapping, first);
    assert_eq!(serial, parallel);
}
