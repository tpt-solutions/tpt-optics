//! Zernike wavefront-error merit.
//!
//! Samples the optical path difference (OPD) of a pupil grid of rays through a
//! [`SequentialSystem`], referenced to the chief ray's image point, and fits
//! the result with OSA/ANSI Zernike polynomials over the unit pupil. The fit
//! gives a first-class wavefront merit (RMS error and per-term coefficients)
//! instead of the spot-size proxy in [`merit`](crate::merit).

use alloc::{vec, vec::Vec};

use tpt_optics_core::Vec3;
use tpt_optics_geometric::ray::Ray;
use tpt_optics_geometric::sequential::{SequentialSystem, TraceContext};
use tpt_optics_geometric::surface::zernike_eval;

/// One pupil sample: normalized pupil coordinates `(x, y)` in the unit disk
/// and the measured wavefront error `w` (meters of OPD).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PupilSample {
    /// Normalized pupil x (|(x, y)| ≤ 1).
    pub x: f64,
    /// Normalized pupil y.
    pub y: f64,
    /// Wavefront error at the sample (meters).
    pub w: f64,
}

/// Sample a hexapolar pupil grid: the centre plus `rings` rings with
/// `6·k` points on ring `k`, returning normalized `(x, y)` coordinates.
#[must_use]
pub fn hexapolar_grid(rings: usize) -> Vec<(f64, f64)> {
    let mut pts = vec![(0.0, 0.0)];
    for k in 1..=rings {
        let rho = k as f64 / rings as f64;
        let n = 6 * k;
        for i in 0..n {
            let theta = core::f64::consts::TAU * i as f64 / n as f64;
            pts.push((rho * libm::cos(theta), rho * libm::sin(theta)));
        }
    }
    pts
}

/// Wavefront error of every grid ray at one field, in meters.
///
/// Rays start on a flat plane (so their incident OPD is uniform) and are
/// referenced to a plane through the chief-ray image point, perpendicular to
/// the chief exit direction: `W = OPD_i + n·s_i − OPD_chief`, where `s_i` is
/// the distance from the ray's image-plane hit to that reference plane.
/// Blocked rays are dropped.
///
/// `aperture` is the physical pupil semi-diameter (meters) mapping normalized
/// `(x, y)` to launch positions.
#[must_use]
pub fn wavefront_samples(
    system: &SequentialSystem,
    ctx: TraceContext,
    field_angle_deg: f64,
    grid: &[(f64, f64)],
    aperture: f64,
) -> Vec<PupilSample> {
    let theta = field_angle_deg.to_radians();
    let dir = Vec3::new(0.0, libm::sin(theta), libm::cos(theta));
    let start_z = -0.1;

    let Some(last) = system.len().checked_sub(1) else {
        return Vec::new();
    };
    let n_image = system.index_after(last, ctx);

    let chief = system.trace(&Ray::new(Vec3::new(0.0, 0.0, start_z), dir), ctx);
    if chief.blocked {
        return Vec::new();
    }
    let (Some(&p_chief), exit_c) = (chief.points.last(), chief.exit_direction) else {
        return Vec::new();
    };
    let opd_chief = chief.opd;

    grid.iter()
        .filter_map(|&(x, y)| {
            let origin = Vec3::new(x * aperture, y * aperture, start_z);
            let res = system.trace(&Ray::new(origin, dir), ctx);
            if res.blocked {
                return None;
            }
            let p = *res.points.last()?;
            let d = res.exit_direction;
            // Distance along the ray to the chief-centred reference plane.
            let denom = d.dot(exit_c);
            if denom.abs() < 1e-15 {
                return None;
            }
            let s = (p_chief - p).dot(exit_c) / denom;
            let w = res.opd + n_image * s - opd_chief;
            Some(PupilSample { x, y, w })
        })
        .collect()
}

/// Least-squares Zernike fit `w(x, y) ≈ Σ_j c_j Z_j(x, y)` for `j = 0..=max_j`.
///
/// Returns the coefficient vector (index = OSA `j`, meters). Solves the normal
/// equations with Gaussian elimination and partial pivoting. Returns an empty
/// vector if there are no samples or the basis is degenerate.
#[must_use]
pub fn fit_zernike(samples: &[PupilSample], max_j: u32) -> Vec<f64> {
    let n = max_j as usize + 1;
    if samples.is_empty() {
        return Vec::new();
    }
    // Normal equations AᵀA c = Aᵀw, with A_ij = Z_j(x_i, y_i).
    let mut ata = vec![vec![0.0f64; n]; n];
    let mut atw = vec![0.0f64; n];
    let mut row = vec![0.0f64; n];
    for s in samples {
        let rho = libm::sqrt(s.x * s.x + s.y * s.y);
        let theta = libm::atan2(s.y, s.x);
        for (j, r) in row.iter_mut().enumerate() {
            *r = zernike_eval(j as u32, rho, theta);
        }
        for a in 0..n {
            atw[a] += row[a] * s.w;
            for b in 0..n {
                ata[a][b] += row[a] * row[b];
            }
        }
    }
    solve_linear(ata, atw).unwrap_or_default()
}

/// RMS of the wavefront after removing the best-fit Zernike terms `j ≥ 0`
/// up to `max_j` — i.e. the residual error the basis cannot express.
///
/// Use `max_j = 0` for the plain RMS about the mean (piston removed).
#[must_use]
pub fn wavefront_rms(samples: &[PupilSample], max_j: u32) -> f64 {
    if samples.is_empty() {
        return f64::NAN;
    }
    let coeffs = fit_zernike(samples, max_j);
    let mut sum = 0.0;
    for s in samples {
        let rho = libm::sqrt(s.x * s.x + s.y * s.y);
        let theta = libm::atan2(s.y, s.x);
        let model: f64 = coeffs
            .iter()
            .enumerate()
            .map(|(j, c)| c * zernike_eval(j as u32, rho, theta))
            .sum();
        let r = s.w - model;
        sum += r * r;
    }
    libm::sqrt(sum / samples.len() as f64)
}

/// Gaussian elimination with partial pivoting for a small dense system.
fn solve_linear(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = b.len();
    for col in 0..n {
        let pivot = (col..n).max_by(|&i, &j| a[i][col].abs().total_cmp(&a[j][col].abs()))?;
        if a[pivot][col].abs() < 1e-300 {
            return None;
        }
        a.swap(col, pivot);
        b.swap(col, pivot);
        let (top, rest) = a.split_at_mut(col + 1);
        let pivot_row = &top[col];
        let pivot_rhs = b[col];
        for (row, rhs) in rest.iter_mut().zip(b[col + 1..].iter_mut()) {
            let f = row[col] / pivot_row[col];
            for (x, p) in row[col..].iter_mut().zip(&pivot_row[col..]) {
                *x -= f * p;
            }
            *rhs -= f * pivot_rhs;
        }
    }
    let mut x = vec![0.0f64; n];
    for r in (0..n).rev() {
        let tail: f64 = (r + 1..n).map(|c| a[r][c] * x[c]).sum();
        x[r] = (b[r] - tail) / a[r][r];
    }
    Some(x)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_optics_core::{Kelvin, Wavelength};
    use tpt_optics_geometric::sequential::{Medium, SeqSurface};
    use tpt_optics_geometric::surface::Asphere;

    fn ctx() -> TraceContext {
        TraceContext {
            wavelength: Wavelength::from_nm(587.6),
            temperature: Kelvin::from_celsius(20.0),
        }
    }

    #[test]
    fn fit_recovers_known_coefficients() {
        // Synthetic wavefront: 2 µm of Z4 (defocus) plus 0.5 µm of Z7 (coma).
        let grid = hexapolar_grid(8);
        let samples: Vec<PupilSample> = grid
            .iter()
            .map(|&(x, y)| {
                let rho = libm::sqrt(x * x + y * y);
                let th = libm::atan2(y, x);
                let w = 2e-6 * zernike_eval(4, rho, th) + 0.5e-6 * zernike_eval(7, rho, th);
                PupilSample { x, y, w }
            })
            .collect();
        let c = fit_zernike(&samples, 8);
        assert!((c[4] - 2e-6).abs() < 1e-12, "defocus: {}", c[4]);
        assert!((c[7] - 0.5e-6).abs() < 1e-12, "coma: {}", c[7]);
        assert!(c[1].abs() < 1e-12 && c[2].abs() < 1e-12);
        // Once the true terms are in the basis, nothing is left over.
        assert!(wavefront_rms(&samples, 8) < 1e-12);
    }

    #[test]
    fn flat_plate_has_no_wavefront_error() {
        // A plane-parallel plate at normal incidence adds no aberration.
        let mut sys = SequentialSystem::new();
        sys.push(SeqSurface::flat(10e-3), 5e-3, Medium::Const(1.5));
        sys.push(SeqSurface::flat(10e-3), 0.0, Medium::Vacuum);
        let grid = hexapolar_grid(4);
        let samples = wavefront_samples(&sys, ctx(), 0.0, &grid, 5e-3);
        assert_eq!(samples.len(), grid.len());
        assert!(wavefront_rms(&samples, 0) < 1e-12);
    }

    #[test]
    fn spherical_lens_shows_wavefront_error() {
        // A strongly curved singlet has spherical aberration, so the
        // wavefront RMS must be clearly non-zero.
        let mut sys = SequentialSystem::new();
        sys.push(
            SeqSurface::aspheric(Asphere::sphere(10.0), 12e-3),
            5e-3,
            Medium::Const(1.5168),
        );
        sys.push(
            SeqSurface::aspheric(Asphere::sphere(-10.0), 12e-3),
            40e-3,
            Medium::Vacuum,
        );
        sys.push(SeqSurface::flat(12e-3), 0.0, Medium::Vacuum);
        let grid = hexapolar_grid(6);
        let samples = wavefront_samples(&sys, ctx(), 0.0, &grid, 5e-3);
        assert!(!samples.is_empty());
        assert!(wavefront_rms(&samples, 0) > 1e-9);
    }

    #[test]
    fn grid_has_expected_size() {
        // 1 + 6·(1 + 2 + … + rings) points.
        assert_eq!(hexapolar_grid(1).len(), 7);
        assert_eq!(hexapolar_grid(3).len(), 1 + 6 * (1 + 2 + 3));
    }

    #[test]
    fn solver_handles_centre_only_samples() {
        // Samples all at the pupil centre: only piston is determined, and the
        // fit must return that constant instead of failing on the singular basis.
        // Here Z_j(0, 0) is non-zero only for j = 0 (ρ = 0 ⇒ Z_j = 0 for j > 0).
        let samples = vec![
            PupilSample {
                x: 0.0,
                y: 0.0,
                w: 1.0
            };
            4
        ];
        let c = fit_zernike(&samples, 0);
        assert_eq!(c.len(), 1);
        assert!((c[0] - 1.0).abs() < 1e-12);
    }
}
