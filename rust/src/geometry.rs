//! Rigid fitting and coordinate error calculations.

use nalgebra::{Matrix3, Vector3};

use crate::error::{DockError, Result};
use crate::types::{Point, RigidTransform};

/// Fit a proper rigid transform taking `source` coordinates to `target`.
///
/// Uses the Kabsch algorithm. Rank deficient point sets are accepted when the
/// SVD yields a finite result; reflections are corrected to keep the rotation
/// proper. The convention is `target = rotation * source + translation`.
pub fn fit(source: &[Point], target: &[Point]) -> Result<RigidTransform> {
    let mut transform = fit_transform(source, target)?;
    transform.rmsd = transformed_rmsd(source, target, &transform)?;
    Ok(transform)
}

/// Fit a rigid transform without calculating RMSD for the fitted coordinates.
/// Use this when only the transform is needed, such as aligning the receptor
/// before measuring ligand RMSD.
pub fn fit_transform(source: &[Point], target: &[Point]) -> Result<RigidTransform> {
    validate_pair_lengths(source, target, "rigid fit")?;

    let source_center = centroid(source, "rigid fit")?;
    let target_center = centroid(target, "rigid fit")?;
    if !source_center.iter().all(|value| value.is_finite())
        || !target_center.iter().all(|value| value.is_finite())
    {
        return Err(DockError("rigid fit centroid is non-finite".into()));
    }
    let mut covariance = Matrix3::<f64>::zeros();
    for (source_point, target_point) in source.iter().zip(target) {
        let source_delta = to_vector(*source_point) - source_center;
        let target_delta = to_vector(*target_point) - target_center;
        covariance += source_delta * target_delta.transpose();
    }
    if !covariance.iter().all(|value| value.is_finite()) {
        return Err(DockError("rigid fit covariance is non-finite".into()));
    }

    let svd = covariance.svd(true, true);
    let u = svd
        .u
        .ok_or_else(|| DockError("SVD did not return U".into()))?;
    let v_t = svd
        .v_t
        .ok_or_else(|| DockError("SVD did not return V transpose".into()))?;
    let v = v_t.transpose();

    let mut correction = Matrix3::<f64>::identity();
    if (v * u.transpose()).determinant() < 0.0 {
        correction[(2, 2)] = -1.0;
    }
    let rotation = v * correction * u.transpose();
    let translation = target_center - rotation * source_center;

    if !rotation.iter().all(|value| value.is_finite())
        || !translation.iter().all(|value| value.is_finite())
    {
        return Err(DockError(
            "rigid fit produced a non-finite transform".into(),
        ));
    }

    Ok(RigidTransform {
        rotation: matrix_to_array(&rotation),
        translation: [translation.x, translation.y, translation.z],
        rmsd: 0.0,
    })
}

/// Apply a source-to-target transform to one point.
pub fn apply(point: Point, transform: &RigidTransform) -> Point {
    apply_with_components(
        point,
        &rotation_matrix(transform),
        Vector3::from(transform.translation),
    )
}

/// Apply one transform to a coordinate set, constructing its matrix once.
pub fn apply_all(points: &[Point], transform: &RigidTransform) -> Vec<Point> {
    let rotation = rotation_matrix(transform);
    let translation = Vector3::from(transform.translation);
    points
        .iter()
        .map(|&point| apply_with_components(point, &rotation, translation))
        .collect()
}

fn rotation_matrix(transform: &RigidTransform) -> Matrix3<f64> {
    Matrix3::from_row_slice(&[
        transform.rotation[0][0],
        transform.rotation[0][1],
        transform.rotation[0][2],
        transform.rotation[1][0],
        transform.rotation[1][1],
        transform.rotation[1][2],
        transform.rotation[2][0],
        transform.rotation[2][1],
        transform.rotation[2][2],
    ])
}

fn apply_with_components(
    point: Point,
    rotation: &Matrix3<f64>,
    translation: Vector3<f64>,
) -> Point {
    let result = rotation * to_vector(point) + translation;
    [result.x, result.y, result.z]
}

fn transformed_rmsd(source: &[Point], target: &[Point], transform: &RigidTransform) -> Result<f64> {
    let rotation = rotation_matrix(transform);
    let translation = Vector3::from(transform.translation);
    let sum_squared = source
        .iter()
        .zip(target)
        .map(|(&source_point, target_point)| {
            let fitted = rotation * to_vector(source_point) + translation;
            let dx = fitted.x - target_point[0];
            let dy = fitted.y - target_point[1];
            let dz = fitted.z - target_point[2];
            dx * dx + dy * dy + dz * dz
        })
        .sum::<f64>();
    let value = (sum_squared / source.len() as f64).sqrt();
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DockError("rigid fit RMSD is non-finite".into()))
    }
}

/// Calculate RMSD after applying a previously fitted transform, without
/// allocating a temporary transformed coordinate set.
pub(crate) fn rmsd_after_transform(
    source: &[Point],
    target: &[Point],
    transform: &RigidTransform,
) -> Result<f64> {
    validate_pair_sets(source, target, "RMSD")?;
    transformed_rmsd(source, target, transform)
}

/// Return RMSD between two coordinate sets without refitting them.
pub fn rmsd(source: &[Point], target: &[Point]) -> Result<f64> {
    validate_pair_sets(source, target, "RMSD")?;
    let sum_squared = source
        .iter()
        .zip(target)
        .map(|(a, b)| {
            let dx = a[0] - b[0];
            let dy = a[1] - b[1];
            let dz = a[2] - b[2];
            dx * dx + dy * dy + dz * dz
        })
        .sum::<f64>();
    let value = (sum_squared / source.len() as f64).sqrt();
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DockError(
            "RMSD calculation produced a non-finite value".into(),
        ))
    }
}

fn validate_pair_sets(source: &[Point], target: &[Point], operation: &str) -> Result<()> {
    validate_pair_lengths(source, target, operation)?;
    if source
        .iter()
        .chain(target)
        .flat_map(|point| point.iter())
        .any(|coordinate| !coordinate.is_finite())
    {
        return Err(DockError(format!("{operation} coordinates must be finite")));
    }
    Ok(())
}

fn validate_pair_lengths(source: &[Point], target: &[Point], operation: &str) -> Result<()> {
    if source.is_empty() || target.is_empty() {
        return Err(DockError(format!(
            "{operation} requires non-empty coordinate sets"
        )));
    }
    if source.len() != target.len() {
        return Err(DockError(format!(
            "{operation} coordinate set lengths differ ({} vs {})",
            source.len(),
            target.len()
        )));
    }
    Ok(())
}

fn centroid(points: &[Point], operation: &str) -> Result<Vector3<f64>> {
    let mut sum = Vector3::zeros();
    for point in points {
        if !point.iter().all(|coordinate| coordinate.is_finite()) {
            return Err(DockError(format!("{operation} coordinates must be finite")));
        }
        sum += to_vector(*point);
    }
    Ok(sum / points.len() as f64)
}

fn to_vector(point: Point) -> Vector3<f64> {
    Vector3::new(point[0], point[1], point[2])
}

fn matrix_to_array(matrix: &Matrix3<f64>) -> [[f64; 3]; 3] {
    [
        [matrix[(0, 0)], matrix[(0, 1)], matrix[(0, 2)]],
        [matrix[(1, 0)], matrix[(1, 1)], matrix[(1, 2)]],
        [matrix[(2, 0)], matrix[(2, 1)], matrix[(2, 2)]],
    ]
}

#[cfg(test)]
mod tests {
    use super::{apply, fit, rmsd};
    use crate::types::Point;

    const TOL: f64 = 1e-10;

    fn close(a: f64, b: f64) {
        assert!((a - b).abs() < TOL, "{a} is not close to {b}");
    }

    fn assert_points_close(actual: &[Point], expected: &[Point]) {
        assert_eq!(actual.len(), expected.len());
        for (actual, expected) in actual.iter().zip(expected) {
            for axis in 0..3 {
                close(actual[axis], expected[axis]);
            }
        }
    }

    #[test]
    fn fits_known_rotation_and_translation_source_to_target() {
        let source = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 3.0],
        ];
        // A +90-degree rotation around z, followed by translation [4, -2, 7].
        let target = [
            [4.0, -2.0, 7.0],
            [4.0, -1.0, 7.0],
            [2.0, -2.0, 7.0],
            [4.0, -2.0, 10.0],
        ];
        let transform = fit(&source, &target).unwrap();
        let transformed: Vec<_> = source
            .iter()
            .map(|&point| apply(point, &transform))
            .collect();
        assert_points_close(&transformed, &target);
        close(transform.rmsd, 0.0);
        close(transform.rotation[0][1], -1.0);
        close(transform.rotation[1][0], 1.0);
        for (actual, expected) in transform.translation.iter().zip([4.0, -2.0, 7.0]) {
            close(*actual, expected);
        }
    }

    #[test]
    fn reflection_is_corrected_to_a_proper_rotation() {
        let source = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let reflected: Vec<_> = source.iter().map(|p| [-p[0], p[1], p[2]]).collect();
        let transform = fit(&source, &reflected).unwrap();
        let r = transform.rotation;
        let determinant = r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1])
            - r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0])
            + r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0]);
        close(determinant, 1.0);
        assert!(transform.rmsd > 0.0);
    }

    #[test]
    fn rank_deficient_collinear_fit_is_finite() {
        let source = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let target = [[2.0, 3.0, 4.0], [2.0, 4.0, 4.0], [2.0, 6.0, 4.0]];
        let transform = fit(&source, &target).unwrap();
        let transformed: Vec<_> = source
            .iter()
            .map(|&point| apply(point, &transform))
            .collect();
        assert_points_close(&transformed, &target);
        assert!(transform.rmsd.is_finite());
    }

    #[test]
    fn direct_rmsd_does_not_refit() {
        let source = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let translated = [[0.0, 2.0, 0.0], [1.0, 2.0, 0.0]];
        close(rmsd(&source, &translated).unwrap(), 2.0);
        close(fit(&source, &translated).unwrap().rmsd, 0.0);
    }

    #[test]
    fn empty_and_mismatched_coordinate_sets_are_errors() {
        assert!(fit(&[], &[]).is_err());
        assert!(rmsd(&[[0.0, 0.0, 0.0]], &[]).is_err());
        assert!(fit(&[[0.0, 0.0, 0.0]], &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]).is_err());
    }
}
