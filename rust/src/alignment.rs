//! Sequence and residue-number alignment.

use crate::error::{DockError, Result};
use crate::types::{Alignment, AlignmentColumn, AlignmentMode, Chain};

const MATCH_SCORE: i64 = 10;
const MISMATCH_SCORE: i64 = 0;
const GAP_OPEN: i64 = -8;
const GAP_EXTEND: i64 = -1;
const NEGATIVE_INFINITY: i64 = i64::MIN / 4;
const MATCH: u8 = 0;
const GAP_IN_NATIVE: u8 = 1;
const GAP_IN_MODEL: u8 = 2;

#[derive(Clone, Copy)]
struct Cell {
    score: i64,
    previous: u8,
}

impl Cell {
    const UNREACHABLE: Self = Self {
        score: NEGATIVE_INFINITY,
        previous: MATCH,
    };
}

/// Align model and native chains with the baseline affine-gap scoring.
///
/// Scores are doubled to represent the half-point extension exactly. Ties use
/// the first-optimum order observed from Bio.Align.PairwiseAligner: match,
/// gap in native, then gap in model.
pub fn align_chains(model: &Chain, native: &Chain, mode: AlignmentMode) -> Result<Alignment> {
    let (model_sequence, native_sequence) = match mode {
        AlignmentMode::Sequence => (
            model.sequence.chars().collect::<Vec<_>>(),
            native.sequence.chars().collect::<Vec<_>>(),
        ),
        AlignmentMode::Numbering => numbering_sequences(model, native)?,
    };
    if model_sequence.is_empty() || native_sequence.is_empty() {
        return Err(DockError(
            "cannot align a chain with an empty sequence".into(),
        ));
    }
    if model_sequence.len() != model.residues.len()
        || native_sequence.len() != native.residues.len()
    {
        return Err(DockError(format!(
            "alignment sequence lengths must match residue counts (model {} vs {}, native {} vs {})",
            model_sequence.len(), model.residues.len(), native_sequence.len(), native.residues.len()
        )));
    }
    if model_sequence == native_sequence {
        return Ok(Alignment {
            columns: (0..model_sequence.len())
                .map(|i| AlignmentColumn {
                    model: Some(i),
                    native: Some(i),
                    matched: true,
                })
                .collect(),
            mismatches: 0,
        });
    }
    affine_global_alignment(&model_sequence, &native_sequence)
}

fn numbering_sequences(model: &Chain, native: &Chain) -> Result<(Vec<char>, Vec<char>)> {
    if model.residues.is_empty() || native.residues.is_empty() {
        return Err(DockError(
            "cannot align a chain with no residues by numbering".into(),
        ));
    }
    let minimum = model
        .residues
        .iter()
        .chain(&native.residues)
        .map(|r| i64::from(r.number))
        .min()
        .unwrap();
    let shift = 45i64.max(-minimum);
    let encode = |chain: &Chain| -> Result<Vec<char>> {
        chain
            .residues
            .iter()
            .map(|residue| {
                let codepoint = i64::from(residue.number) + shift;
                u32::try_from(codepoint)
                    .ok()
                    .and_then(char::from_u32)
                    .ok_or_else(|| {
                        DockError(format!(
                            "residue number {} cannot be encoded for numbering alignment",
                            residue.number
                        ))
                    })
            })
            .collect()
    };
    Ok((encode(model)?, encode(native)?))
}

fn affine_global_alignment(model: &[char], native: &[char]) -> Result<Alignment> {
    let rows = model
        .len()
        .checked_add(1)
        .ok_or_else(|| DockError("alignment dimensions overflow".into()))?;
    let cols = native
        .len()
        .checked_add(1)
        .ok_or_else(|| DockError("alignment dimensions overflow".into()))?;
    let size = rows
        .checked_mul(cols)
        .ok_or_else(|| DockError("alignment dimensions overflow".into()))?;
    let mut matched = vec![Cell::UNREACHABLE; size];
    let mut gap_native = vec![Cell::UNREACHABLE; size];
    let mut gap_model = vec![Cell::UNREACHABLE; size];
    let at = |i: usize, j: usize| i * cols + j;
    matched[0] = Cell {
        score: 0,
        previous: MATCH,
    };
    for i in 1..rows {
        gap_native[at(i, 0)] = Cell {
            score: GAP_OPEN + (i as i64 - 1) * GAP_EXTEND,
            previous: if i == 1 { MATCH } else { GAP_IN_NATIVE },
        };
    }
    for j in 1..cols {
        gap_model[at(0, j)] = Cell {
            score: GAP_OPEN + (j as i64 - 1) * GAP_EXTEND,
            previous: if j == 1 { MATCH } else { GAP_IN_MODEL },
        };
    }

    for i in 1..rows {
        for j in 1..cols {
            let here = at(i, j);
            let diag = at(i - 1, j - 1);
            let substitution = if model[i - 1] == native[j - 1] {
                MATCH_SCORE
            } else {
                MISMATCH_SCORE
            };
            matched[here] = choose([
                (matched[diag], substitution, MATCH),
                (gap_native[diag], substitution, GAP_IN_NATIVE),
                (gap_model[diag], substitution, GAP_IN_MODEL),
            ]);
            let above = at(i - 1, j);
            gap_native[here] = choose([
                (matched[above], GAP_OPEN, MATCH),
                (gap_native[above], GAP_EXTEND, GAP_IN_NATIVE),
                (gap_model[above], GAP_OPEN, GAP_IN_MODEL),
            ]);
            let left = at(i, j - 1);
            gap_model[here] = choose([
                (matched[left], GAP_OPEN, MATCH),
                (gap_native[left], GAP_OPEN, GAP_IN_NATIVE),
                (gap_model[left], GAP_EXTEND, GAP_IN_MODEL),
            ]);
        }
    }
    let last = at(rows - 1, cols - 1);
    let terminal = choose([
        (matched[last], 0, MATCH),
        (gap_native[last], 0, GAP_IN_NATIVE),
        (gap_model[last], 0, GAP_IN_MODEL),
    ]);
    let mut state = terminal.previous;
    let (mut i, mut j) = (model.len(), native.len());
    let mut columns = Vec::with_capacity(model.len().max(native.len()));
    let mut mismatches = 0;
    while i > 0 || j > 0 {
        let cell = at(i, j);
        match state {
            MATCH if i > 0 && j > 0 => {
                let (mi, ni) = (i - 1, j - 1);
                let is_match = model[mi] == native[ni];
                mismatches += usize::from(!is_match);
                columns.push(AlignmentColumn {
                    model: Some(mi),
                    native: Some(ni),
                    matched: is_match,
                });
                state = matched[cell].previous;
                i -= 1;
                j -= 1;
            }
            GAP_IN_NATIVE if i > 0 => {
                columns.push(AlignmentColumn {
                    model: Some(i - 1),
                    native: None,
                    matched: false,
                });
                state = gap_native[cell].previous;
                i -= 1;
            }
            GAP_IN_MODEL if j > 0 => {
                columns.push(AlignmentColumn {
                    model: None,
                    native: Some(j - 1),
                    matched: false,
                });
                state = gap_model[cell].previous;
                j -= 1;
            }
            _ => return Err(DockError("invalid alignment traceback state".into())),
        }
    }
    columns.reverse();
    Ok(Alignment {
        columns,
        mismatches,
    })
}

/// Keep the first max, in match / gap-native / gap-model predecessor order.
fn choose(candidates: [(Cell, i64, u8); 3]) -> Cell {
    let mut result = Cell::UNREACHABLE;
    for (previous, score_delta, previous_state) in candidates {
        let score = previous.score + score_delta;
        if score > result.score {
            result = Cell {
                score,
                previous: previous_state,
            };
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::align_chains;
    use crate::types::{AlignmentColumn, AlignmentMode, Chain, Residue};

    fn chain(id: &str, sequence: &str, numbers: &[i32]) -> Chain {
        Chain {
            id: id.into(),
            sequence: sequence.into(),
            hetero: None,
            residues: numbers
                .iter()
                .map(|&number| Residue {
                    name: "ALA".into(),
                    number,
                    insertion_code: String::new(),
                    hetero: false,
                    atoms: Vec::new(),
                })
                .collect(),
        }
    }

    fn aligned(model: Option<usize>, native: Option<usize>, matched: bool) -> AlignmentColumn {
        AlignmentColumn {
            model,
            native,
            matched,
        }
    }

    fn formatted(
        alignment: &crate::types::Alignment,
        model: &str,
        native: &str,
    ) -> (String, String, String) {
        let model: Vec<_> = model.chars().collect();
        let native: Vec<_> = native.chars().collect();
        let mut seq_a = String::new();
        let mut seq_b = String::new();
        let mut matches = String::new();
        for column in &alignment.columns {
            seq_a.push(column.model.map_or('-', |i| model[i]));
            seq_b.push(column.native.map_or('-', |i| native[i]));
            matches.push(if column.matched {
                '|'
            } else if column.model.is_none() || column.native.is_none() {
                ' '
            } else {
                '.'
            });
        }
        (seq_a, matches, seq_b)
    }

    #[test]
    fn repeated_sequence_tie_matches_reference_first_optimum() {
        let model = chain("M", "AAAA", &[1, 2, 3, 4]);
        let native = chain("N", "AAA", &[1, 2, 3]);
        let alignment = align_chains(&model, &native, AlignmentMode::Sequence).unwrap();
        assert_eq!(
            alignment.columns,
            [
                aligned(Some(0), None, false),
                aligned(Some(1), Some(0), true),
                aligned(Some(2), Some(1), true),
                aligned(Some(3), Some(2), true),
            ]
        );
        assert_eq!(alignment.matched_pairs(), [(1, 0), (2, 1), (3, 2)]);
    }

    #[test]
    fn fixture_tie_matches_biopython_formatted_alignment() {
        let fixture: serde_json::Value =
            serde_json::from_str(include_str!("../../compat/reference.json")).unwrap();
        let record = &fixture["alignments"][0];
        let model_sequence = record["model"].as_str().unwrap();
        let native_sequence = record["native"].as_str().unwrap();
        let model = chain("M", model_sequence, &[0, 1, 2, 3]);
        let native = chain("N", native_sequence, &[0, 1, 2]);
        let alignment = align_chains(&model, &native, AlignmentMode::Sequence).unwrap();
        assert_eq!(
            formatted(&alignment, model_sequence, native_sequence).0,
            record["seqA"]
        );
        assert_eq!(
            formatted(&alignment, model_sequence, native_sequence).1,
            record["matches"]
        );
        assert_eq!(
            formatted(&alignment, model_sequence, native_sequence).2,
            record["seqB"]
        );
    }

    #[test]
    fn exhaustive_small_binary_sequences_match_biopython_ties() {
        let fixture: serde_json::Value =
            serde_json::from_str(include_str!("../../compat/alignment_ties.json")).unwrap();
        for record in fixture["alignments"].as_array().unwrap() {
            let model_sequence = record["model"].as_str().unwrap();
            let native_sequence = record["native"].as_str().unwrap();
            let model_numbers: Vec<_> = (0..model_sequence.chars().count())
                .map(|i| i as i32)
                .collect();
            let native_numbers: Vec<_> = (0..native_sequence.chars().count())
                .map(|i| i as i32)
                .collect();
            let model = chain("M", model_sequence, &model_numbers);
            let native = chain("N", native_sequence, &native_numbers);
            let alignment = align_chains(&model, &native, AlignmentMode::Sequence).unwrap();
            assert_eq!(
                formatted(&alignment, model_sequence, native_sequence),
                (
                    record["model_aligned"].as_str().unwrap().to_string(),
                    record["matches"].as_str().unwrap().to_string(),
                    record["native_aligned"].as_str().unwrap().to_string(),
                ),
                "model={model_sequence}, native={native_sequence}"
            );
        }
    }

    #[test]
    fn affine_gap_and_mismatch_columns_match_pairwise_aligner() {
        let model = chain("M", "ACDE", &[1, 2, 3, 4]);
        let native = chain("N", "ACE", &[1, 2, 3]);
        let alignment = align_chains(&model, &native, AlignmentMode::Sequence).unwrap();
        assert_eq!(
            alignment.columns,
            [
                aligned(Some(0), Some(0), true),
                aligned(Some(1), Some(1), true),
                aligned(Some(2), None, false),
                aligned(Some(3), Some(2), true),
            ]
        );
        assert_eq!(alignment.mismatches, 0);

        let model = chain("M", "ABC", &[1, 2, 3]);
        let native = chain("N", "ADC", &[1, 2, 3]);
        let alignment = align_chains(&model, &native, AlignmentMode::Sequence).unwrap();
        assert_eq!(
            alignment.columns,
            [
                aligned(Some(0), Some(0), true),
                aligned(Some(1), Some(1), false),
                aligned(Some(2), Some(2), true),
            ]
        );
        assert_eq!(alignment.mismatches, 1);
        assert_eq!(alignment.matched_pairs(), [(0, 0), (2, 2)]);
    }

    #[test]
    fn numbering_aligns_integer_residue_numbers_and_ignores_sequence() {
        let model = chain("M", "ABCD", &[-2, -1, 0, 1]);
        let native = chain("N", "WXY", &[-2, -1, 0]);
        let alignment = align_chains(&model, &native, AlignmentMode::Numbering).unwrap();
        assert_eq!(
            alignment.columns,
            [
                aligned(Some(0), Some(0), true),
                aligned(Some(1), Some(1), true),
                aligned(Some(2), Some(2), true),
                aligned(Some(3), None, false),
            ]
        );
        assert_eq!(alignment.mismatches, 0);
    }

    #[test]
    fn empty_or_inconsistent_chains_return_errors() {
        let empty = chain("M", "", &[]);
        let one = chain("N", "A", &[1]);
        assert!(align_chains(&empty, &one, AlignmentMode::Sequence).is_err());
        let inconsistent = chain("M", "AB", &[1]);
        assert!(align_chains(&inconsistent, &one, AlignmentMode::Sequence).is_err());
    }
}
