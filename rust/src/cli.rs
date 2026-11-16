//! Command-line argument parsing and executable orchestration.

use std::path::PathBuf;

use crate::error::{DockError, Result};
use crate::output::{self, RenderOptions};
use crate::types::{LoadOptions, MappingConstraints, ScoringOptions, Structure};

const VERSION: &str = env!("CARGO_PKG_VERSION");

#[derive(Debug, Clone, PartialEq)]
pub enum Command {
    Help,
    Version,
    Score(Arguments),
}

#[derive(Debug, Clone, PartialEq)]
pub struct Arguments {
    pub model: PathBuf,
    pub native: PathBuf,
    pub mapping: Option<String>,
    pub json: Option<PathBuf>,
    pub short: bool,
    pub verbose: bool,
    pub capri_peptide: bool,
    pub small_molecule: bool,
    pub no_align: bool,
    pub n_cpu: usize,
    pub max_chunk: usize,
    pub allowed_mismatches: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ChainMappingSpec {
    /// Optional model and native filters for complete or native-subset mappings.
    pub model_chains: Option<Vec<char>>,
    pub native_chains: Option<Vec<char>>,
    /// Explicit fixed pairs in native:model order.
    pub fixed: Vec<(char, char)>,
}

#[derive(Debug, Clone)]
pub struct CliResponse {
    pub stdout: String,
    pub json_file: Option<(PathBuf, String)>,
}

/// Parse arguments without printing or exiting, which also makes help behavior
/// available to embedding callers.
pub fn parse_args<I, S>(args: I) -> Result<Command>
where
    I: IntoIterator<Item = S>,
    S: Into<String>,
{
    let args: Vec<String> = args.into_iter().map(Into::into).collect();
    let mut model = None;
    let mut native = None;
    let mut mapping = None;
    let mut json = None;
    let mut short = false;
    let mut verbose = false;
    let mut capri_peptide = false;
    let mut small_molecule = false;
    let mut no_align = false;
    let mut n_cpu = 8usize;
    let mut max_chunk = 512usize;
    let mut allowed_mismatches = 0usize;
    let mut help = false;
    let mut version = false;
    let mut positional_only = false;
    let mut index = 0;

    while index < args.len() {
        let arg = &args[index];
        if !positional_only && arg == "--" {
            positional_only = true;
            index += 1;
            continue;
        }
        if !positional_only && (arg == "-h" || arg == "--help") {
            help = true;
            index += 1;
            continue;
        }
        if !positional_only && arg == "--version" {
            version = true;
            index += 1;
            continue;
        }
        if !positional_only && (arg == "-v" || arg == "--verbose") {
            verbose = true;
            index += 1;
            continue;
        }
        if !positional_only {
            match arg.as_str() {
                "--short" => short = true,
                "--capri_peptide" => capri_peptide = true,
                "--small_molecule" => small_molecule = true,
                "--no_align" => no_align = true,
                _ => {
                    let (option, inline_value) = split_option(arg);
                    if matches!(
                        option,
                        "--mapping" | "--json" | "--n_cpu" | "--max_chunk" | "--allowed_mismatches"
                    ) {
                        let value = match inline_value {
                            Some(value) => value.to_owned(),
                            None => {
                                index += 1;
                                args.get(index).cloned().ok_or_else(|| {
                                    DockError(format!("{option} requires a value"))
                                })?
                            }
                        };
                        match option {
                            "--mapping" => mapping = Some(value),
                            "--json" => json = Some(PathBuf::from(value)),
                            "--n_cpu" => n_cpu = parse_positive(option, &value)?,
                            "--max_chunk" => max_chunk = parse_positive(option, &value)?,
                            "--allowed_mismatches" => {
                                allowed_mismatches = value.parse().map_err(|_| {
                                    DockError(format!("{option} must be a non-negative integer"))
                                })?;
                            }
                            _ => unreachable!(),
                        }
                        index += 1;
                        continue;
                    } else if arg.starts_with('-') {
                        return Err(DockError(format!("unrecognized option: {arg}")));
                    }
                }
            }
            if matches!(
                arg.as_str(),
                "--short" | "--capri_peptide" | "--small_molecule" | "--no_align"
            ) {
                index += 1;
                continue;
            }
        }

        if model.is_none() {
            model = Some(PathBuf::from(arg));
        } else if native.is_none() {
            native = Some(PathBuf::from(arg));
        } else {
            return Err(DockError(format!("unexpected positional argument: {arg}")));
        }
        index += 1;
    }

    if help && version {
        return Err(DockError(
            "--help and --version cannot be requested together".into(),
        ));
    }
    if help {
        return Ok(Command::Help);
    }
    if version {
        return Ok(Command::Version);
    }

    let arguments = Arguments {
        model: model.ok_or_else(|| DockError("missing required model path".into()))?,
        native: native.ok_or_else(|| DockError("missing required native path".into()))?,
        mapping,
        json,
        short,
        verbose,
        capri_peptide,
        small_molecule,
        no_align,
        n_cpu,
        max_chunk,
        allowed_mismatches,
    };
    if let Some(mapping) = arguments.mapping.as_deref() {
        parse_mapping_spec(mapping)?;
    }
    Ok(Command::Score(arguments))
}

fn split_option(arg: &str) -> (&str, Option<&str>) {
    match arg.split_once('=') {
        Some((option, value)) => (option, Some(value)),
        None => (arg, None),
    }
}

fn parse_positive(option: &str, value: &str) -> Result<usize> {
    let parsed: usize = value
        .parse()
        .map_err(|_| DockError(format!("{option} must be a positive integer")))?;
    if parsed == 0 {
        return Err(DockError(format!("{option} must be greater than zero")));
    }
    Ok(parsed)
}

/// Parse the legacy one-character `MODELCHAINS:NATIVECHAINS` syntax.
pub fn parse_mapping_spec(spec: &str) -> Result<ChainMappingSpec> {
    let (model, native) = spec
        .split_once(':')
        .ok_or_else(|| DockError("--mapping must use MODELCHAINS:NATIVECHAINS syntax".into()))?;
    if native.contains(':') {
        return Err(DockError(
            "--mapping may contain only one ':' separator".into(),
        ));
    }
    if native.is_empty() {
        return Err(DockError(
            "--mapping requires native chains (for example ABC:ABC or :ABC)".into(),
        ));
    }
    let native_chars = mapping_chars(native)?;
    if model.is_empty() || model == "*" {
        if native_chars.contains(&'*') {
            return Err(DockError(
                "wildcards are not valid in the native chain subset".into(),
            ));
        }
        return Ok(ChainMappingSpec {
            model_chains: None,
            native_chains: Some(native_chars),
            fixed: Vec::new(),
        });
    }
    let model_chars = mapping_chars(model)?;
    if model_chars.len() != native_chars.len() {
        return Err(DockError(format!(
            "--mapping sides must have equal length (got {} model and {} native chains)",
            model_chars.len(),
            native_chars.len()
        )));
    }
    if has_duplicates(&model_chars) || has_duplicates(&native_chars) {
        return Err(DockError(
            "--mapping chain identifiers must not repeat".into(),
        ));
    }

    let fixed = model_chars
        .iter()
        .zip(native_chars.iter())
        .filter_map(|(&model_id, &native_id)| {
            (model_id != '*' && native_id != '*').then_some((native_id, model_id))
        })
        .collect();
    let fully_fixed = !model_chars.contains(&'*') && !native_chars.contains(&'*');
    Ok(ChainMappingSpec {
        model_chains: fully_fixed.then_some(model_chars),
        native_chains: fully_fixed.then_some(native_chars),
        fixed,
    })
}

fn mapping_chars(value: &str) -> Result<Vec<char>> {
    if !value.is_ascii() {
        return Err(DockError(
            "--mapping supports single-character ASCII chain IDs; use an unambiguous selection for multi-character IDs".into(),
        ));
    }
    Ok(value.chars().collect())
}

fn has_duplicates(values: &[char]) -> bool {
    values
        .iter()
        .enumerate()
        .any(|(i, value)| values[..i].contains(value))
}

fn resolve_indices(ids: &[char], structure: &Structure, label: &str) -> Result<Vec<usize>> {
    ids.iter()
        .map(|&wanted| {
            structure
                .chains
                .iter()
                .position(|chain| chain.id == wanted.to_string())
                .ok_or_else(|| DockError(format!("{label} chain '{wanted}' was not found")))
        })
        .collect()
}

fn make_constraints(
    spec: Option<&str>,
    model: &Structure,
    native: &Structure,
    allowed_mismatches: usize,
) -> Result<MappingConstraints> {
    let mut constraints = MappingConstraints {
        allowed_mismatches,
        ..MappingConstraints::default()
    };
    let Some(spec) = spec else {
        return Ok(constraints);
    };
    let parsed = parse_mapping_spec(spec)?;
    if let Some(ids) = parsed.model_chains {
        constraints.model_chains = Some(resolve_indices(&ids, model, "model")?);
    }
    if let Some(ids) = parsed.native_chains {
        constraints.native_chains = Some(resolve_indices(&ids, native, "native")?);
    }
    constraints.fixed = parsed
        .fixed
        .into_iter()
        .map(|(native_id, model_id)| {
            let native_index = resolve_indices(&[native_id], native, "native")?[0];
            let model_index = resolve_indices(&[model_id], model, "model")?[0];
            Ok((native_index, model_index))
        })
        .collect::<Result<Vec<_>>>()?;
    Ok(constraints)
}

/// Run the CLI without printing. The executable owns terminal and file output.
pub fn run() -> Result<CliResponse> {
    run_from(std::env::args().skip(1))
}

pub fn run_from<I, S>(args: I) -> Result<CliResponse>
where
    I: IntoIterator<Item = S>,
    S: Into<String>,
{
    match parse_args(args)? {
        Command::Help => Ok(CliResponse {
            stdout: help_text(),
            json_file: None,
        }),
        Command::Version => Ok(CliResponse {
            stdout: format!("DockQ {VERSION}\n"),
            json_file: None,
        }),
        Command::Score(args) => score(args),
    }
}

fn score(args: Arguments) -> Result<CliResponse> {
    let load_options = LoadOptions {
        include_hetero: args.small_molecule,
        auth_chains: !args.small_molecule,
        ..LoadOptions::default()
    };
    let model = crate::io::load::load_structure(&args.model, &load_options)?;
    let native = crate::io::load::load_structure(&args.native, &load_options)?;
    let constraints = make_constraints(
        args.mapping.as_deref(),
        &model,
        &native,
        args.allowed_mismatches,
    )?;
    let scoring = ScoringOptions {
        peptide: args.capri_peptide,
        no_align: args.no_align,
    };
    let evaluation = crate::engine::evaluate(
        &model,
        &native,
        &constraints,
        &scoring,
        args.n_cpu,
        args.max_chunk,
    )?;
    if evaluation.interfaces.is_empty() {
        return Err(DockError(
            "could not find interfaces in the native structure; check the input files or --mapping"
                .into(),
        ));
    }
    let render_options = RenderOptions {
        short: args.short,
        verbose: args.verbose,
        capri_peptide: args.capri_peptide,
        small_molecule: args.small_molecule,
    };
    let stdout = output::render_text(&model, &native, &evaluation, &render_options)?;
    let json_file = args
        .json
        .map(|path| output::render_json(&model, &native, &evaluation).map(|json| (path, json)))
        .transpose()?;
    Ok(CliResponse { stdout, json_file })
}

pub fn help_text() -> String {
    "DockQ - Quality measure for protein-protein docking models\n\nUSAGE:\n    DockQ <model> <native> [OPTIONS]\n\nARGS:\n    <model>    Path to model file\n    <native>   Path to native file\n\nOPTIONS:\n    --capri_peptide                 Use CAPRI peptide definitions\n    --small_molecule                Include and score small-molecule chains\n    --short                         Print compact output\n    --json <out.json>               Write structured results to a JSON file\n    --verbose, -v                   Print scoring definitions in the header\n    --no_align                      Align residues by numbering instead of sequence\n    --n_cpu <CPU>                   Number of workers [default: 8]\n    --max_chunk <CHUNK>             Maximum mapping chunk size [default: 512]\n    --allowed_mismatches <COUNT>    Allowed sequence mismatches [default: 0]\n    --mapping <MODEL:NATIVE>        Constrain chain mapping; IDs are single characters\n    --help, -h                      Print help\n    --version                       Print version\n"
        .to_owned()
}

#[cfg(test)]
mod tests {
    use super::{parse_args, parse_mapping_spec, Command};

    #[test]
    fn help_and_version_do_not_require_paths() {
        assert_eq!(parse_args(["--help"]).unwrap(), Command::Help);
        assert_eq!(parse_args(["--version"]).unwrap(), Command::Version);
    }

    #[test]
    fn parses_legacy_full_mapping_in_native_model_order() {
        let parsed = parse_mapping_spec("AB:YX").unwrap();
        assert_eq!(parsed.model_chains, Some(vec!['A', 'B']));
        assert_eq!(parsed.native_chains, Some(vec!['Y', 'X']));
        assert_eq!(parsed.fixed, vec![('Y', 'A'), ('X', 'B')]);
    }

    #[test]
    fn parses_native_subset_and_partial_wildcards() {
        let subset = parse_mapping_spec(":WX").unwrap();
        assert_eq!(subset.model_chains, None);
        assert_eq!(subset.native_chains, Some(vec!['W', 'X']));
        let partial = parse_mapping_spec("A*:WY").unwrap();
        assert_eq!(partial.model_chains, None);
        assert_eq!(partial.native_chains, None);
        assert_eq!(partial.fixed, vec![('W', 'A')]);
    }

    #[test]
    fn rejects_ambiguous_mapping_input() {
        assert!(parse_mapping_spec("AB").is_err());
        assert!(parse_mapping_spec("AB:").is_err());
        assert!(parse_mapping_spec("AB:W").is_err());
        assert!(parse_mapping_spec("AA:WX").is_err());
        assert!(parse_mapping_spec("AA:W*").is_err());
        assert!(parse_mapping_spec("α:β").is_err());
    }

    #[test]
    fn parses_all_baseline_flags_and_equals_values() {
        let parsed = parse_args([
            "model.pdb",
            "native.pdb",
            "--capri_peptide",
            "--small_molecule",
            "--short",
            "--json=out.json",
            "-v",
            "--no_align",
            "--n_cpu=2",
            "--max_chunk",
            "10",
            "--allowed_mismatches",
            "3",
            "--mapping=AB:AB",
        ])
        .unwrap();
        let Command::Score(parsed) = parsed else {
            panic!("expected score args")
        };
        assert!(parsed.capri_peptide && parsed.small_molecule && parsed.short);
        assert!(parsed.verbose && parsed.no_align);
        assert_eq!(parsed.n_cpu, 2);
        assert_eq!(parsed.max_chunk, 10);
        assert_eq!(parsed.allowed_mismatches, 3);
        assert_eq!(parsed.mapping.as_deref(), Some("AB:AB"));
    }
}
