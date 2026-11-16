//! Benchmark direct versus per-structure cached contacts for all CSV interfaces.
//!
//! Usage: cargo run --release --example benchmark_all -- DATA_ROOT MODE WORKERS SCORES_TSV [COMPLEX]
//! MODE is direct, cached, or compare. Compare runs both paths on the same pairs.

use dockq::error::Result as DockResult;
use dockq::io::load::load_structure;
use dockq::scoring::polymer::{
    contact_matrix, score_interface, score_interface_with_cached_contacts, CachedContacts,
};
use dockq::types::{Chain, InterfaceScore, LoadOptions, ScoringOptions, Structure};
use rayon::prelude::*;
use rayon::ThreadPool;
use rayon::ThreadPoolBuilder;
use std::collections::BTreeMap;
use std::error::Error;
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

#[derive(Clone)]
struct InterfaceSpec {
    interface: String,
    side1: String,
    side2: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Mode {
    Direct,
    Cached,
    Compare,
}

fn csv_fields(line: &str) -> Vec<String> {
    let mut fields = Vec::new();
    let mut field = String::new();
    let mut quoted = false;
    let mut chars = line.chars().peekable();
    while let Some(ch) = chars.next() {
        match ch {
            '"' if quoted && chars.peek() == Some(&'"') => {
                field.push('"');
                chars.next();
            }
            '"' => quoted = !quoted,
            ',' if !quoted => fields.push(std::mem::take(&mut field)),
            _ => field.push(ch),
        }
    }
    fields.push(field);
    fields
}

fn read_manifest(path: &Path) -> Result<BTreeMap<String, Vec<InterfaceSpec>>, Box<dyn Error>> {
    let text = std::fs::read_to_string(path)?;
    let mut grouped: BTreeMap<String, Vec<InterfaceSpec>> = BTreeMap::new();
    for (line_number, line) in text.lines().enumerate().skip(1) {
        let columns = csv_fields(line);
        if columns.len() != 11 {
            return Err(
                format!("CSV line {} has {} fields", line_number + 1, columns.len()).into(),
            );
        }
        let (side1, side2) = columns[2]
            .split_once(" vs ")
            .ok_or_else(|| format!("invalid interface on CSV line {}", line_number + 1))?;
        grouped
            .entry(columns[1].clone())
            .or_default()
            .push(InterfaceSpec {
                interface: columns[2].clone(),
                side1: side1.to_owned(),
                side2: side2.to_owned(),
            });
    }
    Ok(grouped)
}

fn cif_paths(directory: &Path) -> Result<Vec<PathBuf>, Box<dyn Error>> {
    let mut paths = std::fs::read_dir(directory)?
        .map(|entry| entry.map(|e| e.path()))
        .collect::<Result<Vec<_>, _>>()?;
    paths.retain(|path| path.extension().is_some_and(|ext| ext == "cif"));
    paths.sort();
    Ok(paths)
}

fn merge_chains(structure: &Structure, ids: &str) -> Result<Chain, Box<dyn Error>> {
    let mut merged = Chain {
        id: ids.to_owned(),
        residues: Vec::new(),
        sequence: String::new(),
        hetero: None,
    };
    for id in ids.chars() {
        let id_string = id.to_string();
        let chain = structure
            .chains
            .iter()
            .find(|chain| chain.id == id_string)
            .ok_or_else(|| format!("chain {id} is missing from {}", structure.id))?;
        merged.residues.extend(chain.residues.iter().cloned());
        merged.sequence.push_str(&chain.sequence);
    }
    Ok(merged)
}

fn score_pairs(
    interfaces: &[[Chain; 2]],
    pairs: &[(usize, usize)],
    cache: Option<&[CachedContacts]>,
    pool: &ThreadPool,
) -> (Duration, Vec<DockResult<Option<InterfaceScore>>>) {
    let start = Instant::now();
    let results = pool.install(|| {
        pairs
            .par_iter()
            .map(|&(model, native)| match cache {
                Some(matrices) => score_interface_with_cached_contacts(
                    [&interfaces[model][0], &interfaces[model][1]],
                    [&interfaces[native][0], &interfaces[native][1]],
                    &ScoringOptions::default(),
                    &matrices[model],
                    &matrices[native],
                ),
                None => score_interface(
                    [&interfaces[model][0], &interfaces[model][1]],
                    [&interfaces[native][0], &interfaces[native][1]],
                    &ScoringOptions::default(),
                ),
            })
            .collect()
    });
    (start.elapsed(), results)
}

fn collect_scores(
    results: Vec<DockResult<Option<InterfaceScore>>>,
) -> Result<Vec<Option<InterfaceScore>>, Box<dyn Error>> {
    results
        .into_iter()
        .collect::<DockResult<Vec<_>>>()
        .map_err(Into::into)
}

fn precompute_contacts(
    interfaces: &[[Chain; 2]],
    pool: &ThreadPool,
) -> Result<Vec<CachedContacts>, Box<dyn Error>> {
    let results = pool.install(|| {
        interfaces
            .par_iter()
            .map(|chains| contact_matrix([&chains[0], &chains[1]]))
            .collect::<Vec<_>>()
    });
    results
        .into_iter()
        .collect::<DockResult<Vec<_>>>()
        .map_err(Into::into)
}

fn append_scores(
    output: &mut String,
    complex: &str,
    spec: &InterfaceSpec,
    paths: &[PathBuf],
    pairs: &[(usize, usize)],
    scores: &[Option<InterfaceScore>],
) -> (usize, f64) {
    let mut scored = 0usize;
    let mut checksum = 0.0;
    for (&(model, native), result) in pairs.iter().zip(scores) {
        if let Some(result) = result {
            scored += 1;
            checksum += result.dockq;
            output.push_str(&format!(
                "{}\t{}\t{}\t{}\t{:.12}\t{:016x}\n",
                complex,
                spec.interface,
                paths[model].file_name().unwrap().to_string_lossy(),
                paths[native].file_name().unwrap().to_string_lossy(),
                result.dockq,
                result.dockq.to_bits()
            ));
        }
    }
    (scored, checksum)
}

fn main() -> Result<(), Box<dyn Error>> {
    let mut args = std::env::args().skip(1);
    let data_root = PathBuf::from(args.next().ok_or("missing data root")?);
    let mode = match args.next().as_deref() {
        Some("direct") => Mode::Direct,
        Some("cached") => Mode::Cached,
        Some("compare") => Mode::Compare,
        _ => return Err("mode must be direct, cached, or compare".into()),
    };
    let workers: usize = args.next().ok_or("missing worker count")?.parse()?;
    let scores_path = PathBuf::from(args.next().ok_or("missing scores output path")?);
    let complex_filter = args.next();
    if workers == 0 || args.next().is_some() {
        return Err("usage: benchmark_all DATA_ROOT MODE WORKERS SCORES_TSV [COMPLEX]".into());
    }

    let manifest: BTreeMap<_, _> =
        read_manifest(&data_root.join("benchmark_all_88_results_v0.csv"))?
            .into_iter()
            .filter(|(complex, _)| {
                complex_filter
                    .as_ref()
                    .is_none_or(|filter| filter == complex)
            })
            .collect();
    let total_complexes = manifest.len();
    let total_interfaces: usize = manifest.values().map(Vec::len).sum();
    let pool = ThreadPoolBuilder::new().num_threads(workers).build()?;
    #[cfg(feature = "profile-phases")]
    dockq::profile::reset();
    let mut score_output = String::from("complex\tinterface\tmodel\tnative\tDockQ\tDockQ_bits\n");
    let mut total_load = Duration::ZERO;
    let mut total_fusion = Duration::ZERO;
    let mut total_precompute = Duration::ZERO;
    let mut total_direct = Duration::ZERO;
    let mut total_cached = Duration::ZERO;
    let mut total_attempted = 0usize;
    let mut total_scored = 0usize;
    let mut total_mismatches = 0usize;
    let mut max_score_diff = 0.0f64;
    let mut interface_ordinal = 0usize;

    for (complex, specs) in manifest {
        let paths = cif_paths(&data_root.join(&complex))?;
        let load_start = Instant::now();
        let structures = paths
            .iter()
            .map(|path| load_structure(path, &LoadOptions::default()).map_err(Into::into))
            .collect::<Result<Vec<Structure>, Box<dyn Error>>>()?;
        let load_elapsed = load_start.elapsed();
        total_load += load_elapsed;

        for spec in specs {
            let fusion_start = Instant::now();
            let interfaces = structures
                .iter()
                .map(|structure| {
                    Ok::<_, Box<dyn Error>>([
                        merge_chains(structure, &spec.side1)?,
                        merge_chains(structure, &spec.side2)?,
                    ])
                })
                .collect::<Result<Vec<_>, _>>()?;
            let fusion_elapsed = fusion_start.elapsed();
            total_fusion += fusion_elapsed;

            let pairs: Vec<_> = (0..interfaces.len())
                // Deterministic direction: later file is model, earlier file is native.
                .flat_map(|native| (native..interfaces.len()).map(move |model| (model, native)))
                .collect();
            total_attempted += pairs.len();

            let mut precompute_elapsed = Duration::ZERO;
            let matrices = if mode != Mode::Direct {
                let start = Instant::now();
                let matrices = precompute_contacts(&interfaces, &pool)?;
                precompute_elapsed = start.elapsed();
                total_precompute += precompute_elapsed;
                Some(matrices)
            } else {
                None
            };

            let (direct_scores, cached_scores, direct_elapsed, cached_elapsed) = match mode {
                Mode::Direct => {
                    let (elapsed, results) = score_pairs(&interfaces, &pairs, None, &pool);
                    (
                        Some(collect_scores(results)?),
                        None,
                        elapsed,
                        Duration::ZERO,
                    )
                }
                Mode::Cached => {
                    let (elapsed, results) =
                        score_pairs(&interfaces, &pairs, matrices.as_deref(), &pool);
                    (
                        None,
                        Some(collect_scores(results)?),
                        Duration::ZERO,
                        elapsed,
                    )
                }
                Mode::Compare => {
                    // Alternate order by interface to reduce warm-cache/order bias.
                    let cached_first = interface_ordinal % 2 == 1;
                    if cached_first {
                        let (cached_elapsed, cached_results) =
                            score_pairs(&interfaces, &pairs, matrices.as_deref(), &pool);
                        let (direct_elapsed, direct_results) =
                            score_pairs(&interfaces, &pairs, None, &pool);
                        (
                            Some(collect_scores(direct_results)?),
                            Some(collect_scores(cached_results)?),
                            direct_elapsed,
                            cached_elapsed,
                        )
                    } else {
                        let (direct_elapsed, direct_results) =
                            score_pairs(&interfaces, &pairs, None, &pool);
                        let (cached_elapsed, cached_results) =
                            score_pairs(&interfaces, &pairs, matrices.as_deref(), &pool);
                        (
                            Some(collect_scores(direct_results)?),
                            Some(collect_scores(cached_results)?),
                            direct_elapsed,
                            cached_elapsed,
                        )
                    }
                }
            };
            total_direct += direct_elapsed;
            total_cached += cached_elapsed;

            if let (Some(direct), Some(cached)) = (&direct_scores, &cached_scores) {
                for (direct, cached) in direct.iter().zip(cached) {
                    match (direct, cached) {
                        (Some(direct), Some(cached)) => {
                            let diff = (direct.dockq - cached.dockq).abs();
                            max_score_diff = max_score_diff.max(diff);
                            if direct != cached {
                                total_mismatches += 1;
                            }
                        }
                        (None, None) => {}
                        _ => total_mismatches += 1,
                    }
                }
            }

            let output_scores = match (&cached_scores, &direct_scores) {
                (Some(cached), _) => cached,
                (_, Some(direct)) => direct,
                _ => unreachable!(),
            };
            let (scored, checksum) = append_scores(
                &mut score_output,
                &complex,
                &spec,
                &paths,
                &pairs,
                output_scores,
            );
            total_scored += scored;
            interface_ordinal += 1;
            println!(
                "complex={} interface=\"{}\" attempted={} scored={} fusion_s={:.3} precompute_s={:.3} direct_s={:.3} cached_pair_s={:.3} output_checksum={:.9} checksums_direct={:.9} checksums_cached={:.9}",
                complex,
                spec.interface,
                pairs.len(),
                scored,
                fusion_elapsed.as_secs_f64(),
                precompute_elapsed.as_secs_f64(),
                direct_elapsed.as_secs_f64(),
                cached_elapsed.as_secs_f64(),
                checksum,
                direct_scores.as_ref().map_or(0.0, |scores| scores.iter().flatten().map(|s| s.dockq).sum()),
                cached_scores.as_ref().map_or(0.0, |scores| scores.iter().flatten().map(|s| s.dockq).sum()),
            );
        }
    }

    std::fs::write(scores_path, score_output)?;
    let direct_total = total_fusion + total_direct;
    let cached_total = total_fusion + total_precompute + total_cached;
    println!(
        "summary mode={:?} complexes={} interfaces={} attempted={} scored={} load_s={:.3} fusion_s={:.3} precompute_s={:.3} direct_s={:.3} cached_pair_s={:.3} direct_total_s={:.3} cached_total_s={:.3} speedup={:.3}x mismatched_pairs={} max_dockq_diff={:.12} workers={}",
        mode,
        total_complexes,
        total_interfaces,
        total_attempted,
        total_scored,
        total_load.as_secs_f64(),
        total_fusion.as_secs_f64(),
        total_precompute.as_secs_f64(),
        total_direct.as_secs_f64(),
        total_cached.as_secs_f64(),
        direct_total.as_secs_f64(),
        cached_total.as_secs_f64(),
        direct_total.as_secs_f64() / cached_total.as_secs_f64(),
        total_mismatches,
        max_score_diff,
        workers,
    );
    #[cfg(feature = "profile-phases")]
    {
        let phase_totals = dockq::profile::totals();
        let measured_phase_nanos: u64 = phase_totals[..6].iter().sum();
        let profile_wall = total_direct + total_cached;
        for (index, phase) in dockq::profile::PHASES.iter().enumerate() {
            let seconds = phase_totals[index] as f64 / 1e9;
            let share = if index == 6 {
                100.0 * seconds / profile_wall.as_secs_f64().max(f64::EPSILON)
            } else {
                100.0 * phase_totals[index] as f64 / measured_phase_nanos.max(1) as f64
            };
            println!(
                "profile_phase={} aggregate_s={seconds:.6} share={share:.1}%",
                phase
            );
        }
    }
    Ok(())
}
