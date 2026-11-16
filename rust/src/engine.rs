//! Ordered search with bounded batches and deterministic parallel reduction.
use crate::error::{DockError, Result};
use crate::mapping::candidate_mappings;
use crate::scoring::{ligand, polymer};
use crate::types::*;
use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::Mutex;

type Key = (usize, usize, usize, usize);
type Cache = Mutex<HashMap<Key, Option<InterfaceScore>>>;
const CACHE_LIMIT: usize = 65_536;

fn mapping_result(
    model: &Structure,
    native: &Structure,
    mapping: &ChainMapping,
    options: &ScoringOptions,
    cache: Option<&Cache>,
) -> Result<Evaluation> {
    let mut interfaces = Vec::new();
    for i in 0..mapping.len() {
        for j in i + 1..mapping.len() {
            let (n1, m1) = mapping[i];
            let (n2, m2) = mapping[j];
            if m1 == m2 {
                continue;
            }
            let nc = [native.chains.get(n1), native.chains.get(n2)];
            let mc = [model.chains.get(m1), model.chains.get(m2)];
            let (Some(nc1), Some(nc2), Some(mc1), Some(mc2)) = (nc[0], nc[1], mc[0], mc[1]) else {
                return Err(DockError("mapping chain index out of range".into()));
            };
            let key = (n1, n2, m1, m2);
            let cached = cache.and_then(|c| c.lock().ok()?.get(&key).cloned());
            let score = if let Some(s) = cached {
                s
            } else {
                let value = if nc1.hetero.is_some() || nc2.hetero.is_some() {
                    ligand::score_interface([mc1, mc2], [nc1, nc2], options)?
                } else {
                    polymer::score_interface([mc1, mc2], [nc1, nc2], options)?
                };
                if let Some(c) = cache {
                    let mut c = c
                        .lock()
                        .map_err(|_| DockError("interface cache lock poisoned".into()))?;
                    // Bound memory without eviction changing deterministic results.
                    if c.len() < CACHE_LIMIT {
                        c.insert(key, value.clone());
                    }
                }
                value
            };
            if let Some(score) = score {
                interfaces.push(InterfaceResult {
                    native_chains: [n1, n2],
                    model_chains: [m1, m2],
                    score,
                });
            }
        }
    }
    let total_dockq = interfaces.iter().map(|i| i.score.dockq).sum::<f64>();
    let global_dockq = if interfaces.is_empty() {
        0.0
    } else {
        total_dockq / interfaces.len() as f64
    };
    Ok(Evaluation {
        mapping: mapping.clone(),
        interfaces,
        total_dockq,
        global_dockq,
    })
}

pub fn evaluate_mapping(
    model: &Structure,
    native: &Structure,
    mapping: &ChainMapping,
    options: &ScoringOptions,
) -> Result<Evaluation> {
    mapping_result(model, native, mapping, options, None)
}

pub fn evaluate(
    model: &Structure,
    native: &Structure,
    constraints: &MappingConstraints,
    options: &ScoringOptions,
    workers: usize,
    max_chunk: usize,
) -> Result<Evaluation> {
    if workers == 0 || max_chunk == 0 {
        return Err(DockError("workers and max_chunk must be positive".into()));
    }
    let mut mappings = candidate_mappings(model, native, constraints)?;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(workers)
        .build()
        .map_err(|e| DockError(e.to_string()))?;
    let batch_size = workers.saturating_mul(max_chunk).min(65_536);
    let cache = Cache::default();
    let mut best: Option<Evaluation> = None;
    loop {
        let batch: Vec<_> = mappings.by_ref().take(batch_size).collect();
        if batch.is_empty() {
            break;
        }
        let results: Vec<Result<Evaluation>> = if workers == 1 {
            batch
                .iter()
                .map(|m| mapping_result(model, native, m, options, Some(&cache)))
                .collect()
        } else {
            pool.install(|| {
                batch
                    .par_iter()
                    .map(|m| mapping_result(model, native, m, options, Some(&cache)))
                    .collect()
            })
        };
        for result in results {
            let result = result?;
            if best
                .as_ref()
                .is_none_or(|b| result.total_dockq > b.total_dockq)
            {
                best = Some(result);
            }
        }
    }
    let best = best
        .ok_or_else(|| DockError("no injective chain mappings satisfy the constraints".into()))?;
    if best.interfaces.is_empty() {
        return Err(DockError(
            "could not find interfaces in the native model; check inputs and --mapping".into(),
        ));
    }
    Ok(best)
}
