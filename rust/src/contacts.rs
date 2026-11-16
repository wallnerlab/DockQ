//! Residue distance and contact-count operations.

use crate::error::{DockError, Result};
use crate::types::{ContactStats, Point, SquaredDistances};

const INITIAL_MIN_DISTANCE: f32 = 100_000.0;
const EARLY_EXIT_DISTANCE: f32 = 1_000.0;

/// Calculate the minimum squared distance between every pair of atom groups.
///
/// Coordinates are rounded to `f32` before subtraction, as in the compiled
/// Python implementation. The final three-term sum is accumulated in `f64`
/// (the Cython `pow` expressions promote to double) and rounded to `f32` when
/// stored. For compatibility, the first newly observed minimum above 1000 Å²
/// ends the atom scan for that residue pair.
pub fn residue_distances(
    a: &[Point],
    b: &[Point],
    groups_a: &[usize],
    groups_b: &[usize],
) -> Result<SquaredDistances> {
    let ends_a = validate_groups(a, groups_a, "first")?;
    let ends_b = validate_groups(b, groups_b, "second")?;
    let capacity = groups_a
        .len()
        .checked_mul(groups_b.len())
        .ok_or_else(|| DockError("residue distance matrix dimensions overflow".into()))?;
    let mut values = Vec::with_capacity(capacity);

    let mut start_a = 0;
    for &end_a in &ends_a {
        let mut start_b = 0;
        for &end_b in &ends_b {
            let mut min_distance = INITIAL_MIN_DISTANCE;
            'atoms_a: for atom_a in &a[start_a..end_a] {
                for atom_b in &b[start_b..end_b] {
                    let dx = atom_a[0] as f32 - atom_b[0] as f32;
                    let dy = atom_a[1] as f32 - atom_b[1] as f32;
                    let dz = atom_a[2] as f32 - atom_b[2] as f32;
                    // Cython's `** 2` calls pow (double), then the result is
                    // converted back to float on assignment to `this_d`.
                    let distance =
                        ((dx as f64).powi(2) + (dy as f64).powi(2) + (dz as f64).powi(2)) as f32;
                    if distance < min_distance {
                        min_distance = distance;
                        if min_distance > EARLY_EXIT_DISTANCE {
                            break 'atoms_a;
                        }
                    }
                }
                // Cython also checks after each outer iteration, including when
                // every inner distance exceeded the initial sentinel.
                if min_distance > EARLY_EXIT_DISTANCE {
                    break 'atoms_a;
                }
            }
            values.push(min_distance);
            start_b = end_b;
        }
        start_a = end_a;
    }

    Ok(SquaredDistances {
        rows: groups_a.len(),
        cols: groups_b.len(),
        values,
    })
}

/// Count native and model contacts and their overlap at a linear Å threshold.
pub fn contact_stats(
    model: &SquaredDistances,
    native: &SquaredDistances,
    threshold: f32,
) -> Result<ContactStats> {
    validate_matrix(model, "model")?;
    validate_matrix(native, "native")?;
    if model.rows != native.rows || model.cols != native.cols {
        return Err(DockError(format!(
            "model and native distance matrices have different dimensions ({}x{} vs {}x{})",
            model.rows, model.cols, native.rows, native.cols
        )));
    }
    let threshold_squared = validate_threshold(threshold)?;
    let mut stats = ContactStats::default();
    for (&model_distance, &native_distance) in model.values.iter().zip(&native.values) {
        let model_contact = model_distance < threshold_squared;
        let native_contact = native_distance < threshold_squared;
        if model_contact {
            stats.model += 1;
            if !native_contact {
                stats.non_native += 1;
            }
        }
        if native_contact {
            stats.native += 1;
            if model_contact {
                stats.shared += 1;
            }
        }
    }
    Ok(stats)
}

/// Count entries strictly below a linear Å threshold squared.
pub fn count_contacts(distances: &SquaredDistances, threshold: f32) -> usize {
    let Ok(threshold_squared) = validate_threshold(threshold) else {
        return 0;
    };
    let Some(expected_len) = distances.rows.checked_mul(distances.cols) else {
        return 0;
    };
    if expected_len != distances.values.len() {
        return 0;
    }
    distances
        .values
        .iter()
        .filter(|&&distance| distance < threshold_squared)
        .count()
}

fn validate_groups(points: &[Point], groups: &[usize], label: &str) -> Result<Vec<usize>> {
    let mut end = 0usize;
    let mut ends = Vec::with_capacity(groups.len());
    for (index, &group_size) in groups.iter().enumerate() {
        if group_size == 0 {
            return Err(DockError(format!(
                "{label} residue group {index} contains no atoms"
            )));
        }
        end = end
            .checked_add(group_size)
            .ok_or_else(|| DockError(format!("{label} residue atom-group sizes overflow")))?;
        if end > points.len() {
            return Err(DockError(format!(
                "{label} residue atom groups require {end} coordinates but only {} were supplied",
                points.len()
            )));
        }
        ends.push(end);
    }
    if end != points.len() {
        return Err(DockError(format!(
            "{label} residue atom groups cover {end} coordinates but {} were supplied",
            points.len()
        )));
    }
    if points
        .iter()
        .flatten()
        .any(|coordinate| !coordinate.is_finite() || !(*coordinate as f32).is_finite())
    {
        return Err(DockError(format!(
            "{label} atom coordinates must be finite when represented as f32"
        )));
    }
    Ok(ends)
}

fn validate_matrix(matrix: &SquaredDistances, label: &str) -> Result<()> {
    let expected_len = matrix
        .rows
        .checked_mul(matrix.cols)
        .ok_or_else(|| DockError(format!("{label} distance matrix dimensions overflow")))?;
    if expected_len != matrix.values.len() {
        return Err(DockError(format!(
            "{label} distance matrix has dimensions {}x{} but {} values",
            matrix.rows,
            matrix.cols,
            matrix.values.len()
        )));
    }
    if matrix.values.iter().any(|value| !value.is_finite()) {
        return Err(DockError(format!(
            "{label} squared distances must be finite"
        )));
    }
    Ok(())
}

fn validate_threshold(threshold: f32) -> Result<f32> {
    if !threshold.is_finite() || threshold < 0.0 {
        return Err(DockError(
            "contact threshold must be finite and nonnegative".into(),
        ));
    }
    Ok(threshold * threshold)
}

#[cfg(test)]
mod tests {
    use super::{contact_stats, count_contacts, residue_distances};
    use crate::error::DockError;
    use crate::types::{Point, SquaredDistances};

    fn matrix(rows: usize, cols: usize, values: &[f32]) -> SquaredDistances {
        SquaredDistances {
            rows,
            cols,
            values: values.to_vec(),
        }
    }

    #[test]
    fn residue_distances_groups_atoms_in_input_order() {
        let a: [Point; 3] = [[0.0, 0.0, 0.0], [9.0, 0.0, 0.0], [5.0, 0.0, 0.0]];
        let b: [Point; 2] = [[2.0, 0.0, 0.0], [20.0, 0.0, 0.0]];
        let result = residue_distances(&a, &b, &[2, 1], &[1, 1]).unwrap();
        assert_eq!((result.rows, result.cols), (2, 2));
        assert_eq!(result.values, [4.0, 121.0, 9.0, 225.0]);
    }

    #[test]
    fn reproduces_compiled_early_exit_discrepancy() {
        let a: [Point; 2] = [[0.0, 0.0, 0.0], [100.0, 0.0, 0.0]];
        let b: [Point; 1] = [[100.0, 0.0, 0.0]];
        let result = residue_distances(&a, &b, &[2], &[1]).unwrap();
        // The first pair is 10000 Å², so Cython stops before finding zero.
        assert_eq!(result.values, [10_000.0]);
    }

    #[test]
    fn threshold_is_linear_strict_and_counts_stats_match_reference() {
        let model = matrix(2, 2, &[24.999, 25.0, 4.0, 100.0]);
        let native = matrix(2, 2, &[25.0, 24.999, 9.0, 100.0]);
        assert_eq!(count_contacts(&model, 5.0), 2);
        assert_eq!(
            contact_stats(&model, &native, 5.0).unwrap(),
            crate::types::ContactStats {
                shared: 1,
                non_native: 1,
                native: 2,
                model: 2,
            }
        );
    }

    #[test]
    fn invalid_group_layout_and_nonfinite_inputs_error() {
        let a: [Point; 1] = [[0.0, 0.0, 0.0]];
        assert!(residue_distances(&a, &a, &[0], &[1]).is_err());
        assert!(residue_distances(&a, &a, &[2], &[1]).is_err());
        let bad: [Point; 1] = [[f64::INFINITY, 0.0, 0.0]];
        assert!(residue_distances(&bad, &a, &[1], &[1]).is_err());
    }

    #[test]
    fn compiled_outer_exit_keeps_sentinel_when_first_atom_exceeds_it() {
        let a = [[1000.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        let b = [[0.0, 0.0, 0.0]];
        let d = residue_distances(&a, &b, &[2], &[1]).unwrap();
        assert_eq!(d.values, vec![100_000.0]);
    }

    #[test]
    fn inconsistent_matrices_and_invalid_threshold_error() {
        let model = matrix(1, 2, &[1.0]);
        let native = matrix(1, 1, &[1.0]);
        assert!(contact_stats(&model, &native, 5.0).is_err());
        let valid = matrix(1, 1, &[1.0]);
        assert_eq!(
            contact_stats(&valid, &valid, f32::NAN).unwrap_err(),
            DockError("contact threshold must be finite and nonnegative".into())
        );
    }
}
