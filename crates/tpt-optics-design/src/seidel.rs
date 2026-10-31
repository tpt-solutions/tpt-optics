//! Spherical-aberration Seidel sum from the exact on-axis wavefront.
//!
//! Zemax's Seidel coefficients are tied to the wavefront coefficients by
//! `S1 = 8·W040` (its Seidel-coefficient page). `W040` is the coefficient of
//! `r⁴` in the on-axis wavefront. This module measures it by a least-squares
//! fit of the exact traced wavefront, so it is checked against the trace.
//! Defocus (the `r²` term) is removed by the fit.
//!
//! This is spherical aberration only. Coma, astigmatism, Petzval and distortion
//! need the off-axis fields, and are not computed here.

use crate::wavefront::{hexapolar_grid, wavefront_samples};
use tpt_optics_geometric::sequential::{SequentialSystem, TraceContext};

/// On-axis spherical wavefront coefficient `K` (in m⁻³), with the wavefront
/// modelled as `w(r) = c₀ + c₂·r² + K·r⁴` for a pupil of radius `aperture`
/// (m). The Seidel sum is `S1 = 8·K` in the same units, following Zemax.
///
/// Returns `None` when too few rays reach the image or the fit is singular.
#[must_use]
pub fn spherical_coefficient(
    system: &SequentialSystem,
    ctx: TraceContext,
    aperture: f64,
    rings: usize,
) -> Option<f64> {
    if !(aperture > 0.0 && aperture.is_finite()) {
        return None;
    }
    let grid = hexapolar_grid(rings);
    let samples = wavefront_samples(system, ctx, 0.0, &grid, aperture);
    if samples.len() < 3 {
        return None;
    }
    // Basis in normalised ρ² = (x² + y²): [1, ρ², ρ⁴]. The wavefront is in
    // metres, so c₄ (for ρ) is K·a⁴.
    let mut ata = [[0.0f64; 3]; 3];
    let mut atw = [0.0f64; 3];
    for s in &samples {
        let r2 = s.x * s.x + s.y * s.y;
        let row = [1.0, r2, r2 * r2];
        for a in 0..3 {
            atw[a] += row[a] * s.w;
            for b in 0..3 {
                ata[a][b] += row[a] * row[b];
            }
        }
    }
    let c = solve3(ata, atw)?;
    let a2 = aperture * aperture;
    Some(c[2] / (a2 * a2))
}

/// Solve a 3 × 3 system by Cramer's rule.
fn solve3(a: [[f64; 3]; 3], b: [f64; 3]) -> Option<[f64; 3]> {
    let det = |m: [[f64; 3]; 3]| {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    };
    let d = det(a);
    if d == 0.0 || !d.is_finite() {
        return None;
    }
    let mut out = [0.0; 3];
    for (k, o) in out.iter_mut().enumerate() {
        let mut m = a;
        for r in 0..3 {
            m[r][k] = b[r];
        }
        *o = det(m) / d;
    }
    Some(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_optics_core::material::MaterialDatabase;
    use tpt_optics_core::{Kelvin, Wavelength};
    use tpt_optics_geometric::presets::singlet;

    fn ctx() -> TraceContext {
        TraceContext {
            wavelength: Wavelength::from_nm(587.6),
            temperature: Kelvin::from_celsius(20.0),
        }
    }

    #[test]
    fn spherical_coefficient_is_independent_of_aperture_for_small_apertures() {
        // K is the r⁴ coefficient. At small apertures the r⁶ terms are
        // negligible, so K measured at 1 mm and at 2 mm must agree.
        let db = MaterialDatabase::presets();
        let sys = singlet(&db, 100e-3, -100e-3, 12.5e-3, 5e-3, "bk7", 95e-3).unwrap();
        let k1 = spherical_coefficient(&sys, ctx(), 1e-3, 6).expect("fit");
        let k2 = spherical_coefficient(&sys, ctx(), 2e-3, 6).expect("fit");
        assert!(k1.is_finite() && k1.abs() > 0.0, "{k1}");
        assert!(
            ((k1 - k2) / k1).abs() < 0.02,
            "K(1 mm) = {k1}, K(2 mm) = {k2}"
        );
    }

    #[test]
    fn a_plane_window_has_no_spherical_aberration() {
        // A plane-parallel plate in air is aberration-free on axis.
        let db = MaterialDatabase::presets();
        let flat = singlet(&db, 0.0, 0.0, 12.5e-3, 5e-3, "bk7", 95e-3).unwrap();
        let lens = singlet(&db, 100e-3, -100e-3, 12.5e-3, 5e-3, "bk7", 95e-3).unwrap();
        let k_flat = spherical_coefficient(&flat, ctx(), 2e-3, 6).expect("fit");
        let k_lens = spherical_coefficient(&lens, ctx(), 2e-3, 6).expect("fit");
        assert!(
            k_flat.abs() < 1e-3 * k_lens.abs(),
            "flat {k_flat} vs lens {k_lens}"
        );
    }

    #[test]
    fn non_positive_aperture_is_refused() {
        let db = MaterialDatabase::presets();
        let sys = singlet(&db, 100e-3, -100e-3, 12.5e-3, 5e-3, "bk7", 95e-3).unwrap();
        assert!(spherical_coefficient(&sys, ctx(), 0.0, 6).is_none());
    }
}
