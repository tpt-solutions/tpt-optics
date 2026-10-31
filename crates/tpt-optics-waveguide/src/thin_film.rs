//! Thin-film stacks with complex indices at oblique incidence, for s and p
//! polarisation, by the characteristic matrix method (Born & Wolf, Macleod).
//!
//! The stack sits between a lossless incident medium `n₀` and a substrate
//! `n_s`, both of which may be absorbing in the substrate's case. For each
//! layer the tilted admittance is `η = n·cos θ` (s) or `n/cos θ` (p), and the
//! phase thickness is `δ = (2π/λ)·n·d·cos θ`. Snell's law fixes `θ` in each
//! layer as `n sin θ = n₀ sin θ₀`, which is complex in an absorbing layer.

use libm::{cos, cosh, sin, sinh};
use tpt_optics_core::{Complex, Wavelength};

/// One layer of a thin-film stack.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct FilmLayer {
    /// Complex refractive index `n + ik` (`k ≥ 0`).
    pub n: Complex<f64>,
    /// Physical thickness (m).
    pub thickness: f64,
}

/// Power reflectances and transmittances of a stack for both polarisations.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct FilmPower {
    /// s-polarised reflectance.
    pub r_s: f64,
    /// p-polarised reflectance.
    pub r_p: f64,
    /// s-polarised transmittance.
    pub t_s: f64,
    /// p-polarised transmittance.
    pub t_p: f64,
}

/// Power response of `layers` between `n0` (real) and the substrate `ns` at
/// incidence angle `theta0` (rad) and wavelength `wl`.
///
/// For a lossless stack `R + T = 1`. With an absorbing layer, `R + T < 1` and
/// the difference is the absorbed power.
#[must_use]
pub fn film_response(
    layers: &[FilmLayer],
    n0: f64,
    ns: Complex<f64>,
    theta0: f64,
    wl: Wavelength,
) -> FilmPower {
    let (r_s, t_s) = polarisation(layers, n0, ns, theta0, wl, true);
    let (r_p, t_p) = polarisation(layers, n0, ns, theta0, wl, false);
    FilmPower { r_s, r_p, t_s, t_p }
}

/// Reflectance and transmittance power for one polarisation.
fn polarisation(
    layers: &[FilmLayer],
    n0: f64,
    ns: Complex<f64>,
    theta0: f64,
    wl: Wavelength,
    s_pol: bool,
) -> (f64, f64) {
    let one = Complex::one();
    let i = Complex::new(0.0, 1.0);
    // Snell invariant n sin θ, real because the incident medium is lossless.
    let invariant = Complex::real(n0 * sin(theta0));
    // cos θ in a medium of index n, on the principal branch.
    let cos_in = |n: Complex<f64>| (one - (invariant / n) * (invariant / n)).sqrt();
    // Tilted admittance.
    let eta = |n: Complex<f64>| {
        let c = cos_in(n);
        if s_pol {
            n * c
        } else {
            n / c
        }
    };
    let k0 = core::f64::consts::TAU / wl.m();

    let mut m = [[one, Complex::zero()], [Complex::zero(), one]];
    for layer in layers {
        // The matrix uses the e^{iωt} convention, where absorption is n − ik.
        let n = layer.n.conj();
        let eta_j = eta(n);
        let delta = Complex::real(k0 * layer.thickness) * n * cos_in(n);
        let (c, s) = (ccos(delta), csin(delta));
        let mj = [[c, i * s / eta_j], [i * eta_j * s, c]];
        m = mat_mul(&m, &mj);
    }

    let eta_0 = eta(Complex::real(n0));
    let eta_s = eta(ns.conj());
    let b = m[0][0] + m[0][1] * eta_s;
    let c = m[1][0] + m[1][1] * eta_s;
    let denom = eta_0 * b + c;
    let r = (eta_0 * b - c) / denom;
    let t = Complex::real(2.0) * eta_0 / denom;
    // Power transmittance: Re(η_s)/Re(η₀)·|t|².
    let t_power = eta_s.re / eta_0.re * t.norm_sqr();
    (r.norm_sqr(), t_power)
}

fn mat_mul(a: &[[Complex<f64>; 2]; 2], b: &[[Complex<f64>; 2]; 2]) -> [[Complex<f64>; 2]; 2] {
    let mut out = [[Complex::zero(); 2]; 2];
    for (i, row) in out.iter_mut().enumerate() {
        for (j, v) in row.iter_mut().enumerate() {
            *v = a[i][0] * b[0][j] + a[i][1] * b[1][j];
        }
    }
    out
}

/// cos(a + ib) = cos a·cosh b − i·sin a·sinh b.
fn ccos(z: Complex<f64>) -> Complex<f64> {
    Complex::new(cos(z.re) * cosh(z.im), -sin(z.re) * sinh(z.im))
}

/// sin(a + ib) = sin a·cosh b + i·cos a·sinh b.
fn csin(z: Complex<f64>) -> Complex<f64> {
    Complex::new(sin(z.re) * cosh(z.im), cos(z.re) * sinh(z.im))
}

/// A quarter-wave layer of index `n` at wavelength `wl`: thickness `λ/(4n)`.
#[must_use]
pub fn quarter_wave_layer(n: f64, wl: Wavelength) -> FilmLayer {
    FilmLayer {
        n: Complex::real(n),
        thickness: wl.m() / (4.0 * n),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_optics_core::fresnel::fresnel_energy;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn bare_interface_matches_fresnel_at_oblique_incidence() {
        let wl = Wavelength::from_nm(550.0);
        let theta = 0.5;
        let p = film_response(&[], 1.0, Complex::real(1.5), theta, wl);
        let (rs, rp, ts, tp) = fresnel_energy(1.0, 1.5, theta);
        assert!(close(p.r_s, rs, 1e-12), "{} vs {rs}", p.r_s);
        assert!(close(p.r_p, rp, 1e-12), "{} vs {rp}", p.r_p);
        assert!(close(p.t_s, ts, 1e-12));
        assert!(close(p.t_p, tp, 1e-12));
    }

    #[test]
    fn p_reflectance_vanishes_at_brewster_angle() {
        let wl = Wavelength::from_nm(550.0);
        let brewster = libm::atan(1.5);
        let p = film_response(&[], 1.0, Complex::real(1.5), brewster, wl);
        assert!(p.r_p < 1e-12, "{}", p.r_p);
        assert!(p.r_s > 0.1);
    }

    #[test]
    fn quarter_wave_anti_reflection_coating_is_perfect() {
        // Index sqrt(ns) with a quarter-wave thickness cancels the reflection.
        let wl = Wavelength::from_nm(550.0);
        let ns = 1.5_f64;
        let layer = quarter_wave_layer(ns.sqrt(), wl);
        let p = film_response(&[layer], 1.0, Complex::real(ns), 0.0, wl);
        assert!(p.r_s < 1e-12 && p.r_p < 1e-12, "{p:?}");
    }

    #[test]
    fn lossless_stack_conserves_power_for_both_polarisations() {
        let wl = Wavelength::from_nm(633.0);
        let layers = [
            FilmLayer {
                n: Complex::real(2.3),
                thickness: 40e-9,
            },
            FilmLayer {
                n: Complex::real(1.38),
                thickness: 90e-9,
            },
            FilmLayer {
                n: Complex::real(2.3),
                thickness: 70e-9,
            },
        ];
        let p = film_response(&layers, 1.0, Complex::real(1.52), 0.7, wl);
        assert!(close(p.r_s + p.t_s, 1.0, 1e-10), "s: {}", p.r_s + p.t_s);
        assert!(close(p.r_p + p.t_p, 1.0, 1e-10), "p: {}", p.r_p + p.t_p);
    }

    #[test]
    fn absorbing_layer_loses_power() {
        let wl = Wavelength::from_nm(550.0);
        let absorber = [FilmLayer {
            n: Complex::new(2.0, 0.5),
            thickness: 100e-9,
        }];
        let p = film_response(&absorber, 1.0, Complex::real(1.5), 0.3, wl);
        assert!(p.r_s + p.t_s < 1.0 - 1e-3, "{}", p.r_s + p.t_s);
        assert!(p.r_s + p.t_s >= 0.0);
    }
}
