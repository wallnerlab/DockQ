//! Ordered, lazy injective chain mapping enumeration.
use crate::alignment::align_chains;
use crate::error::{DockError, Result};
use crate::types::{AlignmentMode, Chain, ChainMapping, MappingConstraints, Structure};

pub struct MappingIter {
    fixed: ChainMapping,
    keys: Vec<usize>,
    candidates: Vec<Vec<usize>>,
    positions: Vec<usize>,
    chosen: Vec<usize>,
    depth: usize,
    reverse: bool,
    done: bool,
}

impl Iterator for MappingIter {
    type Item = ChainMapping;
    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        if self.keys.is_empty() {
            self.done = true;
            return Some(self.fixed.clone());
        }
        loop {
            if self.depth == self.keys.len() {
                let mut out = self.fixed.clone();
                for (key, chosen) in self.keys.iter().zip(&self.chosen) {
                    out.push(if self.reverse {
                        (*chosen, *key)
                    } else {
                        (*key, *chosen)
                    });
                }
                self.depth -= 1;
                self.chosen.pop();
                return Some(out);
            }
            let d = self.depth;
            if self.positions[d] == self.candidates[d].len() {
                self.positions[d] = 0;
                if d == 0 {
                    self.done = true;
                    return None;
                }
                self.depth -= 1;
                self.chosen.pop();
                continue;
            }
            let c = self.candidates[d][self.positions[d]];
            self.positions[d] += 1;
            if self.chosen.contains(&c) {
                continue;
            }
            self.chosen.push(c);
            self.depth += 1;
            if self.depth < self.positions.len() {
                self.positions[self.depth] = 0;
            }
        }
    }
}

fn compatible(model: &Chain, native: &Chain, mismatches: usize) -> Result<bool> {
    match (&model.hetero, &native.hetero) {
        (None, None) => {
            Ok(align_chains(model, native, AlignmentMode::Sequence)?.mismatches <= mismatches)
        }
        (Some(a), Some(b)) => Ok(a == b),
        _ => Ok(false),
    }
}

pub fn candidate_mappings(
    model: &Structure,
    native: &Structure,
    c: &MappingConstraints,
) -> Result<MappingIter> {
    let ms = c
        .model_chains
        .clone()
        .unwrap_or_else(|| (0..model.chains.len()).collect());
    let ns = c
        .native_chains
        .clone()
        .unwrap_or_else(|| (0..native.chains.len()).collect());
    if ms.iter().any(|&i| i >= model.chains.len()) || ns.iter().any(|&i| i >= native.chains.len()) {
        return Err(DockError("selected chain index out of range".into()));
    }
    if ms
        .iter()
        .enumerate()
        .any(|(i, chain)| ms[..i].contains(chain))
        || ns
            .iter()
            .enumerate()
            .any(|(i, chain)| ns[..i].contains(chain))
    {
        return Err(DockError("selected chain indices must be unique".into()));
    }
    let mut fixed_n = Vec::new();
    let mut fixed_m = Vec::new();
    for &(n, m) in &c.fixed {
        if !ns.contains(&n) || !ms.contains(&m) || fixed_n.contains(&n) || fixed_m.contains(&m) {
            return Err(DockError("invalid or duplicate fixed chain mapping".into()));
        }
        fixed_n.push(n);
        fixed_m.push(m);
    }
    let ms: Vec<_> = ms.into_iter().filter(|m| !fixed_m.contains(m)).collect();
    let ns: Vec<_> = ns.into_iter().filter(|n| !fixed_n.contains(n)).collect();
    let reverse = ms.len() < ns.len();
    let keys = if reverse { ms.clone() } else { ns.clone() };
    let pool = if reverse { ns } else { ms };
    let mut candidates = Vec::new();
    for &key in &keys {
        let mut choices = Vec::new();
        for &p in &pool {
            let (m, n) = if reverse { (key, p) } else { (p, key) };
            if compatible(&model.chains[m], &native.chains[n], c.allowed_mismatches)? {
                choices.push(p);
            }
        }
        if choices.is_empty() {
            return Err(DockError(format!(
                "no compatible chain for {}",
                if reverse {
                    &model.chains[key].id
                } else {
                    &native.chains[key].id
                }
            )));
        }
        candidates.push(choices);
    }
    Ok(MappingIter {
        fixed: c.fixed.clone(),
        positions: vec![0; keys.len()],
        keys,
        candidates,
        chosen: Vec::new(),
        depth: 0,
        reverse,
        done: false,
    })
}

/// Count the exact injective mappings accepted by the constraints.
///
/// This consumes the lazy candidate stream without storing mappings, but may
/// take factorial time for broad homomer searches. The evaluation engine does
/// not call this function so it can start scoring immediately.
pub fn count_mappings(
    model: &Structure,
    native: &Structure,
    constraints: &MappingConstraints,
) -> Result<u128> {
    candidate_mappings(model, native, constraints)?.try_fold(0u128, |count, _| {
        count
            .checked_add(1)
            .ok_or_else(|| DockError("mapping count exceeds u128".into()))
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    fn s(ids: &[&str]) -> Structure {
        Structure {
            id: String::new(),
            chains: ids
                .iter()
                .map(|id| Chain {
                    id: (*id).into(),
                    residues: vec![],
                    sequence: "A".into(),
                    hetero: Some("LIG".into()),
                })
                .collect(),
        }
    }
    #[test]
    fn permutations_order_and_fixed_constraints() {
        let m = s(&["A", "B", "C"]);
        let n = s(&["X", "Y"]);
        let v: Vec<_> = candidate_mappings(&m, &n, &MappingConstraints::default())
            .unwrap()
            .collect();
        assert_eq!(v.len(), 6);
        assert_eq!(v[0], vec![(0, 0), (1, 1)]);
        assert_eq!(v[1], vec![(0, 0), (1, 2)]);
        let c = MappingConstraints {
            fixed: vec![(1, 2)],
            ..Default::default()
        };
        let v: Vec<_> = candidate_mappings(&m, &n, &c).unwrap().collect();
        assert_eq!(v, vec![vec![(1, 2), (0, 0)], vec![(1, 2), (0, 1)]]);
    }
    #[test]
    fn reverse_maps_larger_native_to_smaller_model() {
        let v: Vec<_> = candidate_mappings(
            &s(&["A", "B"]),
            &s(&["X", "Y", "Z"]),
            &MappingConstraints::default(),
        )
        .unwrap()
        .collect();
        assert_eq!(v.len(), 6);
        assert_eq!(v[0], vec![(0, 0), (1, 1)]);
        assert_eq!(v[1], vec![(0, 0), (2, 1)]);
    }

    #[test]
    fn duplicate_selected_chain_indices_are_invalid() {
        let model = s(&["A", "B"]);
        let native = s(&["X", "Y"]);
        let constraints = MappingConstraints {
            model_chains: Some(vec![0, 0]),
            ..Default::default()
        };
        assert!(candidate_mappings(&model, &native, &constraints).is_err());

        let constraints = MappingConstraints {
            native_chains: Some(vec![1, 1]),
            ..Default::default()
        };
        assert!(candidate_mappings(&model, &native, &constraints).is_err());
    }

    #[test]
    fn exact_count_covers_forward_and_reverse_search() {
        let model = s(&["A", "B", "C"]);
        let native = s(&["X", "Y"]);
        assert_eq!(
            count_mappings(&model, &native, &MappingConstraints::default()).unwrap(),
            6
        );

        let model = s(&["A", "B"]);
        let native = s(&["X", "Y", "Z"]);
        assert_eq!(
            count_mappings(&model, &native, &MappingConstraints::default()).unwrap(),
            6
        );
    }
}
