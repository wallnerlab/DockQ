//! Time all unordered pairs for one interface in the extracted CIF dataset.
//!
//! Usage: cargo run --release --example benchmark_interface -- DIR SIDE1 SIDE2 [WORKERS] [SCORES_TSV]
//! For example: ... 1S78 AB C 8

use dockq::io::load::load_structure;
use dockq::scoring::polymer::score_interface;
use dockq::types::{Chain, LoadOptions, ScoringOptions, Structure};
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::error::Error;
use std::path::{Path, PathBuf};
use std::time::Instant;

fn merge_chains(structure: &Structure, ids: &str) -> Result<Chain, Box<dyn Error>> {
    let mut merged = Chain {
        id: ids.to_owned(),
        residues: Vec::new(),
        sequence: String::new(),
        hetero: None,
    };
    for id in ids.chars() {
        let chain = structure
            .chains
            .iter()
            .find(|chain| chain.id == id.to_string())
            .ok_or_else(|| format!("chain {id} is missing from {}", structure.id))?;
        merged.residues.extend(chain.residues.iter().cloned());
        merged.sequence.push_str(&chain.sequence);
    }
    Ok(merged)
}

fn cif_paths(directory: &Path) -> Result<Vec<PathBuf>, Box<dyn Error>> {
    let mut paths = std::fs::read_dir(directory)?
        .map(|entry| entry.map(|e| e.path()))
        .collect::<Result<Vec<_>, _>>()?;
    paths.retain(|path| path.extension().is_some_and(|ext| ext == "cif"));
    paths.sort();
    Ok(paths)
}

fn main() -> Result<(), Box<dyn Error>> {
    let mut args = std::env::args().skip(1);
    let directory = PathBuf::from(args.next().ok_or("missing complex directory")?);
    let side1 = args.next().ok_or("missing first interface side")?;
    let side2 = args.next().ok_or("missing second interface side")?;
    let workers = args
        .next()
        .map(|value| value.parse::<usize>())
        .transpose()?
        .unwrap_or(1);
    let scores_path = args.next().map(PathBuf::from);
    if workers == 0 || args.next().is_some() {
        return Err("usage: benchmark_interface DIR SIDE1 SIDE2 [WORKERS] [SCORES_TSV]".into());
    }

    let paths = cif_paths(&directory)?;
    let load_start = Instant::now();
    let structures = paths
        .iter()
        .map(|path| load_structure(path, &LoadOptions::default()).map_err(Into::into))
        .collect::<Result<Vec<Structure>, Box<dyn Error>>>()?;
    let interfaces = structures
        .iter()
        .map(|structure| {
            Ok::<_, Box<dyn Error>>([
                merge_chains(structure, &side1)?,
                merge_chains(structure, &side2)?,
            ])
        })
        .collect::<Result<Vec<_>, _>>()?;
    let load_elapsed = load_start.elapsed();

    let score_start = Instant::now();
    #[cfg(feature = "profile-phases")]
    dockq::profile::reset();
    let pairs: Vec<_> = (0..interfaces.len())
        // Keep a fixed orientation because DockQ scores are directional.
        .flat_map(|native| (native..interfaces.len()).map(move |model| (model, native)))
        .collect();
    let pool = ThreadPoolBuilder::new().num_threads(workers).build()?;
    let results: Vec<_> = pool.install(|| {
        pairs
            .par_iter()
            .map(|&(model, native)| {
                score_interface(
                    [&interfaces[model][0], &interfaces[model][1]],
                    [&interfaces[native][0], &interfaces[native][1]],
                    &ScoringOptions::default(),
                )
            })
            .collect()
    });
    if let Some(Ok(Some(sample))) = results.get(1) {
        println!(
            "sample_pair={} vs {} DockQ={:.12}",
            paths[1].display(),
            paths[0].display(),
            sample.dockq
        );
    }
    let mut scored = 0usize;
    let mut checksum = 0.0;
    let mut score_rows = String::from("model\tnative\tDockQ\n");
    for (&(model, native), result) in pairs.iter().zip(results) {
        if let Some(result) = result? {
            scored += 1;
            checksum += result.dockq;
            score_rows.push_str(&format!(
                "{}\t{}\t{:.12}\n",
                paths[model].file_name().unwrap().to_string_lossy(),
                paths[native].file_name().unwrap().to_string_lossy(),
                result.dockq
            ));
        }
    }
    if let Some(path) = scores_path {
        std::fs::write(path, score_rows)?;
    }
    let score_elapsed = score_start.elapsed();
    let pairs = interfaces.len() * (interfaces.len() + 1) / 2;
    println!(
        "structures={} pairs={} workers={} scored={} load_s={:.3} score_s={:.3} pairs_per_s={:.1} checksum={:.9}",
        interfaces.len(),
        pairs,
        workers,
        scored,
        load_elapsed.as_secs_f64(),
        score_elapsed.as_secs_f64(),
        pairs as f64 / score_elapsed.as_secs_f64(),
        checksum,
    );
    #[cfg(feature = "profile-phases")]
    {
        let totals = dockq::profile::totals();
        let aggregate_nanos: u64 = totals[..6].iter().sum();
        for (index, label) in dockq::profile::PHASES.iter().enumerate() {
            let seconds = totals[index] as f64 / 1e9;
            let share = if index == 6 {
                100.0 * seconds / score_elapsed.as_secs_f64()
            } else {
                100.0 * totals[index] as f64 / aggregate_nanos as f64
            };
            println!("phase={label} aggregate_s={seconds:.3} share={share:.1}%");
        }
    }
    Ok(())
}
