//! Protein/nucleic acid and CAPRI-peptide interface scoring.
use crate::alignment::align_chains;
use crate::contacts::{contact_stats, count_contacts, residue_distances};
use crate::error::Result;
use crate::geometry::{fit, fit_transform, rmsd_after_transform};
use crate::selection::{
    backbone_pairs, backbone_pairs_cached, cache_backbone, residue_coordinates,
    CachedBackboneResidue,
};
use crate::types::*;
use std::borrow::Cow;
#[cfg(feature = "profile-phases")]
use std::time::Instant;

pub fn dockq_formula(fnat: f64, irms: f64, lrms: f64) -> f64 {
    (fnat + 1.0 / (1.0 + (irms / 1.5).powi(2)) + 1.0 / (1.0 + (lrms / 8.5).powi(2))) / 3.0
}

/// Per-structure distances and summaries reused by all-against-all scoring.
/// The summary masks use the same strict, squared-distance cutoffs as scoring.
#[derive(Debug, Clone)]
pub struct CachedContacts {
    distances: SquaredDistances,
    contacts_4: Vec<u64>,
    contacts_5: Vec<u64>,
    count_2: usize,
    count_4: usize,
    count_5: usize,
    interface_rows: Vec<usize>,
    interface_cols: Vec<usize>,
    backbone: [Vec<CachedBackboneResidue>; 2],
    identity_pairs: [Vec<(usize, usize)>; 2],
}

impl CachedContacts {
    fn has_shape(&self, rows: usize, cols: usize) -> bool {
        rows.checked_mul(cols) == Some(self.distances.values.len())
            && self.distances.rows == rows
            && self.distances.cols == cols
    }

    fn contacts_below(&self, threshold: f32) -> Option<(&[u64], usize)> {
        match threshold {
            4.0 => Some((&self.contacts_4, self.count_4)),
            5.0 => Some((&self.contacts_5, self.count_5)),
            _ => None,
        }
    }

    fn count_below(&self, threshold: f32) -> Option<usize> {
        match threshold {
            2.0 => Some(self.count_2),
            4.0 => Some(self.count_4),
            5.0 => Some(self.count_5),
            _ => None,
        }
    }
}

fn distances(
    a: &Chain,
    b: &Chain,
    ia: &[usize],
    ib: &[usize],
    cb: bool,
) -> Result<SquaredDistances> {
    let (ca, ga) = residue_coordinates(a, ia, cb)?;
    let (cc, gb) = residue_coordinates(b, ib, cb)?;
    residue_distances(&ca, &cc, &ga, &gb)
}

/// Build the all-residue atom-distance matrix for a two-chain interface.
///
/// All-against-all comparisons can retain this matrix for each input structure
/// and reuse it as either the model or native contact matrix.
pub fn contact_matrix(chains: [&Chain; 2]) -> Result<CachedContacts> {
    let all: [Vec<usize>; 2] = std::array::from_fn(|i| (0..chains[i].residues.len()).collect());
    let distances = distances(chains[0], chains[1], &all[0], &all[1], false)?;
    let backbone = [cache_backbone(chains[0]), cache_backbone(chains[1])];
    let identity_pairs: [Vec<(usize, usize)>; 2] = std::array::from_fn(|i| {
        (0..chains[i].residues.len())
            .map(|index| (index, index))
            .collect()
    });
    let word_count = distances.values.len() / 64 + usize::from(distances.values.len() % 64 != 0);
    let mut contacts_4 = vec![0u64; word_count];
    let mut contacts_5 = vec![0u64; word_count];
    let mut count_2 = 0;
    let mut count_4 = 0;
    let mut count_5 = 0;
    let mut rows_within_interface = vec![false; distances.rows];
    let mut cols_within_interface = vec![false; distances.cols];
    for (index, &distance) in distances.values.iter().enumerate() {
        let bit = 1u64 << (index % 64);
        if distance < 4.0 {
            count_2 += 1;
        }
        if distance < 16.0 {
            contacts_4[index / 64] |= bit;
            count_4 += 1;
        }
        if distance < 25.0 {
            contacts_5[index / 64] |= bit;
            count_5 += 1;
        }
        if distance < 100.0 {
            rows_within_interface[index / distances.cols] = true;
            cols_within_interface[index % distances.cols] = true;
        }
    }
    let interface_rows = rows_within_interface
        .iter()
        .enumerate()
        .filter_map(|(index, &selected)| selected.then_some(index))
        .collect();
    let interface_cols = cols_within_interface
        .iter()
        .enumerate()
        .filter_map(|(index, &selected)| selected.then_some(index))
        .collect();
    Ok(CachedContacts {
        distances,
        contacts_4,
        contacts_5,
        count_2,
        count_4,
        count_5,
        interface_rows,
        interface_cols,
        backbone,
        identity_pairs,
    })
}

pub fn score_interface(
    model: [&Chain; 2],
    native: [&Chain; 2],
    options: &ScoringOptions,
) -> Result<Option<InterfaceScore>> {
    score_interface_inner(model, native, options, None)
}

/// Score an interface using cached all-residue distance matrices when the
/// residue mapping is the identity. Non-identity alignments fall back to the
/// normal pair-specific distance calculations.
pub fn score_interface_with_cached_contacts(
    model: [&Chain; 2],
    native: [&Chain; 2],
    options: &ScoringOptions,
    model_contacts: &CachedContacts,
    native_contacts: &CachedContacts,
) -> Result<Option<InterfaceScore>> {
    score_interface_inner(
        model,
        native,
        options,
        Some((model_contacts, native_contacts)),
    )
}

fn score_interface_inner(
    model: [&Chain; 2],
    native: [&Chain; 2],
    options: &ScoringOptions,
    cached_contacts: Option<(&CachedContacts, &CachedContacts)>,
) -> Result<Option<InterfaceScore>> {
    #[cfg(feature = "profile-phases")]
    let total_start = Instant::now();
    let threshold = if options.peptide { 4.0 } else { 5.0 };
    let interface_threshold: f32 = if options.peptide { 8.0 } else { 10.0 };
    #[cfg(feature = "profile-phases")]
    let phase_start = Instant::now();
    let cached_native = cached_contacts
        .map(|(_, native)| native)
        .filter(|cached_native| {
            cached_native.has_shape(native[0].residues.len(), native[1].residues.len())
        });
    let mut nd = match cached_native {
        Some(cached_native) => Cow::Borrowed(&cached_native.distances),
        _ => {
            let all_n: [Vec<usize>; 2] =
                std::array::from_fn(|i| (0..native[i].residues.len()).collect());
            Cow::Owned(distances(
                native[0], native[1], &all_n[0], &all_n[1], false,
            )?)
        }
    };
    let nat_total = cached_native
        .and_then(|cached| cached.count_below(threshold))
        .unwrap_or_else(|| count_contacts(&nd, threshold));
    #[cfg(feature = "profile-phases")]
    crate::profile::record(0, phase_start);
    if nat_total == 0 {
        #[cfg(feature = "profile-phases")]
        crate::profile::record(6, total_start);
        return Ok(None);
    }
    let mode = if options.no_align {
        AlignmentMode::Numbering
    } else {
        AlignmentMode::Sequence
    };
    #[cfg(feature = "profile-phases")]
    let phase_start = Instant::now();
    let cached_model_candidate = cached_contacts
        .map(|(model, _)| model)
        .filter(|cached_model| {
            cached_model.has_shape(model[0].residues.len(), model[1].residues.len())
        });
    let identical_sequences = !options.no_align
        && (0..2).all(|chain| {
            model[chain].sequence == native[chain].sequence
                && model[chain].residues.len() == native[chain].residues.len()
                && model[chain].sequence.chars().count() == model[chain].residues.len()
                && native[chain].sequence.chars().count() == native[chain].residues.len()
        });
    let pairs: [Cow<'_, [(usize, usize)]>; 2] = if identical_sequences {
        match cached_model_candidate {
            Some(cached_model) => [
                Cow::Borrowed(&cached_model.identity_pairs[0]),
                Cow::Borrowed(&cached_model.identity_pairs[1]),
            ],
            None => [
                Cow::Owned(identity_pairs(model[0].residues.len())),
                Cow::Owned(identity_pairs(model[1].residues.len())),
            ],
        }
    } else {
        [
            Cow::Owned(align_chains(model[0], native[0], mode)?.matched_pairs()),
            Cow::Owned(align_chains(model[1], native[1], mode)?.matched_pairs()),
        ]
    };
    #[cfg(feature = "profile-phases")]
    crate::profile::record(1, phase_start);
    #[cfg(feature = "profile-phases")]
    let phase_start = Instant::now();
    let identity_mapping = identical_sequences
        || (0..2).all(|chain| {
            mapping_is_identity(
                pairs[chain].as_ref(),
                model[chain].residues.len(),
                native[chain].residues.len(),
            )
        });
    let cached_model = cached_model_candidate.filter(|_| identity_mapping);
    let cached_backbones = cached_model
        .zip(cached_native)
        .filter(|(model_cache, native_cache)| {
            model_cache.distances.rows == native_cache.distances.rows
                && model_cache.distances.cols == native_cache.distances.cols
        });
    let md = match cached_model {
        Some(cached_model) => Cow::Borrowed(&cached_model.distances),
        _ => {
            let mi: [Vec<usize>; 2] =
                std::array::from_fn(|i| pairs[i].iter().map(|pair| pair.0).collect());
            Cow::Owned(distances(model[0], model[1], &mi[0], &mi[1], false)?)
        }
    };
    #[cfg(feature = "profile-phases")]
    crate::profile::record(2, phase_start);
    #[cfg(feature = "profile-phases")]
    let phase_start = Instant::now();
    if nd.rows != md.rows || nd.cols != md.cols {
        let ni: [Vec<usize>; 2] =
            std::array::from_fn(|i| pairs[i].iter().map(|pair| pair.1).collect());
        nd = Cow::Owned(distances(native[0], native[1], &ni[0], &ni[1], false)?);
    }
    let cached_stats = if options.peptide {
        None
    } else {
        cached_model
            .zip(cached_native)
            .filter(|(model_cache, native_cache)| {
                model_cache.distances.rows == native_cache.distances.rows
                    && model_cache.distances.cols == native_cache.distances.cols
            })
            .and_then(|(model_cache, native_cache)| {
                Some((
                    model_cache.contacts_below(threshold)?,
                    native_cache.contacts_below(threshold)?,
                ))
            })
    };
    let stats = if let Some(((model_mask, model_count), (native_mask, native_count))) = cached_stats
    {
        let shared = model_mask
            .iter()
            .zip(native_mask)
            .map(|(model_bits, native_bits)| (model_bits & native_bits).count_ones() as usize)
            .sum();
        ContactStats {
            shared,
            non_native: model_count - shared,
            native: native_count,
            model: model_count,
        }
    } else {
        contact_stats(md.as_ref(), nd.as_ref(), threshold)?
    };
    let fnat = stats.shared as f64 / nat_total as f64;
    let fnonnat = if stats.model == 0 {
        0.0
    } else {
        stats.non_native as f64 / stats.model as f64
    };
    if options.peptide {
        let ni: [Vec<usize>; 2] =
            std::array::from_fn(|i| pairs[i].iter().map(|pair| pair.1).collect());
        nd = Cow::Owned(distances(native[0], native[1], &ni[0], &ni[1], true)?);
    }
    #[cfg(feature = "profile-phases")]
    crate::profile::record(3, phase_start);
    #[cfg(feature = "profile-phases")]
    let phase_start = Instant::now();
    let mut subsets = [Vec::new(), Vec::new()];
    if let Some((model_cache, native_cache)) =
        cached_model.zip(cached_native).filter(|_| !options.peptide)
    {
        if model_cache.distances.rows == native_cache.distances.rows
            && model_cache.distances.cols == native_cache.distances.cols
        {
            subsets[0].extend_from_slice(&native_cache.interface_rows);
            subsets[1].extend_from_slice(&native_cache.interface_cols);
        }
    }
    if subsets[0].is_empty() && subsets[1].is_empty() {
        for row in 0..nd.rows {
            for col in 0..nd.cols {
                if nd.values[row * nd.cols + col] < interface_threshold * interface_threshold {
                    subsets[0].push(row);
                    subsets[1].push(col);
                }
            }
        }
    }
    let mut inter = AtomPairs::default();
    for i in 0..2 {
        // Python subset_atoms treats an empty collection like None (all pairs).
        let subset = if subsets[i].is_empty() {
            None
        } else {
            Some(subsets[i].as_slice())
        };
        let ap = match cached_backbones {
            Some((model_cache, native_cache)) => backbone_pairs_cached(
                &model_cache.backbone[i],
                &native_cache.backbone[i],
                &model[i].id,
                &native[i].id,
                pairs[i].as_ref(),
                subset,
            )?,
            None => backbone_pairs(model[i], native[i], pairs[i].as_ref(), subset)?,
        };
        inter.model.extend(ap.model);
        inter.native.extend(ap.native);
    }
    #[cfg(feature = "profile-phases")]
    crate::profile::record(4, phase_start);
    #[cfg(feature = "profile-phases")]
    let phase_start = Instant::now();
    let irms = fit(&inter.model, &inter.native)?.rmsd;
    // Equal length uses group 2 as receptor. Never sort chain pairs here.
    let receptor = if native[0].residues.len() > native[1].residues.len() {
        0
    } else {
        1
    };
    let ligand = 1 - receptor;
    let (rec, lig) = match cached_backbones {
        Some((model_cache, native_cache)) => (
            backbone_pairs_cached(
                &model_cache.backbone[receptor],
                &native_cache.backbone[receptor],
                &model[receptor].id,
                &native[receptor].id,
                pairs[receptor].as_ref(),
                None,
            )?,
            backbone_pairs_cached(
                &model_cache.backbone[ligand],
                &native_cache.backbone[ligand],
                &model[ligand].id,
                &native[ligand].id,
                pairs[ligand].as_ref(),
                None,
            )?,
        ),
        None => (
            backbone_pairs(
                model[receptor],
                native[receptor],
                pairs[receptor].as_ref(),
                None,
            )?,
            backbone_pairs(model[ligand], native[ligand], pairs[ligand].as_ref(), None)?,
        ),
    };
    let transform = fit_transform(&rec.model, &rec.native)?;
    let lrms = rmsd_after_transform(&lig.model, &lig.native, &transform)?;
    #[cfg(feature = "profile-phases")]
    crate::profile::record(5, phase_start);
    let denominator = stats.shared + stats.non_native + nat_total;
    let f1 = 2.0 * stats.shared as f64 / denominator as f64;
    let metrics = PolymerMetrics {
        irmsd: irms,
        lrmsd: lrms,
        fnat,
        fnonnat,
        f1,
        clashes: cached_model
            .and_then(|cached| cached.count_below(2.0))
            .unwrap_or_else(|| count_contacts(md.as_ref(), 2.0)),
        nat_correct: stats.shared,
        nat_total,
        nonnat_count: stats.non_native,
        model_total: stats.model,
        len1: native[0].residues.len(),
        len2: native[1].residues.len(),
        class1: if receptor == 0 { "receptor" } else { "ligand" }.into(),
        class2: if receptor == 1 { "receptor" } else { "ligand" }.into(),
    };
    let result = Some(InterfaceScore {
        dockq: dockq_formula(fnat, irms, lrms),
        metrics: Metrics::Polymer(metrics),
    });
    #[cfg(feature = "profile-phases")]
    crate::profile::record(6, total_start);
    Ok(result)
}

fn mapping_is_identity(pairs: &[(usize, usize)], model_len: usize, native_len: usize) -> bool {
    model_len == native_len
        && pairs.len() == model_len
        && pairs
            .iter()
            .enumerate()
            .all(|(index, &(model, native))| model == index && native == index)
}

fn identity_pairs(length: usize) -> Vec<(usize, usize)> {
    (0..length).map(|index| (index, index)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn formula_limits() {
        assert_eq!(dockq_formula(1.0, 0.0, 0.0), 1.0);
        assert!((dockq_formula(0.0, 1.5, 8.5) - 1.0 / 3.0).abs() < 1e-15);
    }
}
