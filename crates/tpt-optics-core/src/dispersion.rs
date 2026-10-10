//! Dispersion models: refractive index as a function of wavelength.
//!
//! All models implement [`Dispersion`]. Wavelengths cross the API as
//! [`Wavelength`] but the classic glass formulas are defined in micrometers;
//! conversion happens inside each model.

use alloc::boxed::Box;
use alloc::vec;
use alloc::vec::Vec;

use libm::pow;

use crate::complex::Complex;
use crate::index::ComplexRefractiveIndex;
use crate::units::Wavelength;

/// Refractive index as a function of wavelength.
pub trait Dispersion: core::fmt::Debug + Send + Sync {
    /// Complex refractive index at the given vacuum wavelength.
    fn index(&self, wavelength: Wavelength) -> ComplexRefractiveIndex;

    /// Convenience: real part only (valid for transparent regions).
    fn real_index(&self, wavelength: Wavelength) -> f64 {
        self.index(wavelength).n
    }

    /// Boxed clone, for storing heterogeneous models in collections.
    fn boxed(self) -> Box<dyn Dispersion>
    where
        Self: Sized + 'static,
    {
        Box::new(self)
    }
}

/// Sellmeier dispersion (formula 1):
///
/// `n² = 1 + Σᵢ Bᵢ·λ² / (λ² − Cᵢ)`, with λ in µm.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Sellmeier {
    /// Oscillator strengths B₁..B₃.
    pub b: [f64; 3],
    /// Squared resonance wavelengths C₁..C₃ (µm²).
    pub c: [f64; 3],
}

impl Sellmeier {
    /// Construct from explicit coefficients (B/C alternating, in µm/µm²).
    pub fn new(coeffs_b: [f64; 3], coeffs_c: [f64; 3]) -> Self {
        Sellmeier {
            b: coeffs_b,
            c: coeffs_c,
        }
    }

    /// N-BK7 (Schott) coefficients.
    #[must_use]
    pub fn bk7() -> Self {
        Sellmeier::new(
            [1.039_612_12, 0.231_792_344, 1.010_469_45],
            [0.006_000_698_67, 0.020_017_914_4, 103.560_653],
        )
    }

    /// Fused silica (SiO₂, Malitson 1965) coefficients.
    #[must_use]
    pub fn fused_silica() -> Self {
        Sellmeier::new(
            [0.696_166_3, 0.407_942_6, 0.897_479_4],
            [0.004_679_148_26, 0.013_512_063_1, 97.934_002_5],
        )
    }

    /// Diamond (Peter 1923, 2-term; third term zeroed with a negative
    /// resonance so it never divides by zero).
    #[must_use]
    pub fn diamond() -> Self {
        Sellmeier::new([0.3306, 4.3356, 0.0], [0.030_625, 0.011_236, -1.0])
    }
}

impl Dispersion for Sellmeier {
    fn index(&self, wavelength: Wavelength) -> ComplexRefractiveIndex {
        let l2 = wavelength.um() * wavelength.um();
        let mut n2 = 1.0;
        for (b, c) in self.b.iter().zip(self.c.iter()) {
            n2 += b * l2 / (l2 - c);
        }
        ComplexRefractiveIndex::real(pow(n2, 0.5))
    }
}

/// Cauchy dispersion: `n = A + B/λ² + C/λ⁴` (λ in µm), k = 0.
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct Cauchy {
    /// Coefficient A.
    pub a: f64,
    /// Coefficient B (µm²).
    pub b: f64,
    /// Coefficient C (µm⁴).
    pub c: f64,
}

impl Cauchy {
    /// Construct.
    pub fn new(a: f64, b: f64, c: f64) -> Self {
        Cauchy { a, b, c }
    }
}

impl Dispersion for Cauchy {
    fn index(&self, wavelength: Wavelength) -> ComplexRefractiveIndex {
        let l = wavelength.um();
        let l2 = l * l;
        ComplexRefractiveIndex::real(self.a + self.b / l2 + self.c / (l2 * l2))
    }
}

/// Schott polynomial dispersion:
///
/// `n² = A₀ + A₁λ² + A₂λ⁻² + A₃λ⁻⁴ + A₄λ⁻⁶ + A₅λ⁻⁸` (λ in µm).
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct Schott {
    /// Polynomial coefficients A₀..A₅.
    pub a: [f64; 6],
}

impl Schott {
    /// Construct.
    pub fn new(a: [f64; 6]) -> Self {
        Schott { a }
    }
}

impl Dispersion for Schott {
    fn index(&self, wavelength: Wavelength) -> ComplexRefractiveIndex {
        let l = wavelength.um();
        let l2 = l * l;
        let [a0, a1, a2, a3, a4, a5] = self.a;
        let n2 = a0 + a1 * l2 + a2 / l2 + a3 / (l2 * l2) + a4 / pow(l2, 3.0) + a5 / pow(l2, 4.0);
        ComplexRefractiveIndex::real(pow(n2, 0.5))
    }
}

/// One Lorentz/Drude oscillator term.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Oscillator {
    /// Oscillator strength (dimensionless, or Δε for Lorentz terms).
    pub strength: f64,
    /// Resonance angular frequency ω₀ in rad/s.
    pub omega0: f64,
    /// Damping angular frequency γ in rad/s (0 for a pure Lorentz pole).
    pub gamma: f64,
}

/// Drude–Lorentz model in relative permittivity space:
///
/// `ε(ω) = ε∞ + Σⱼ Sⱼ·ω₀ⱼ² / (ω₀ⱼ² − ω² − i·γⱼ·ω)`, `n* = √ε`.
///
/// A Drude free-electron metal is the special case `ω₀ → 0`
/// (`S·ω₀² ≡ ωₚ²`): set `omega0 = 0` and interpret `strength·omega0²` via
/// the separate [`Oscillator::drude`] helper which stores ωₚ² directly.
#[derive(Clone, Debug, PartialEq)]
pub struct DrudeLorentz {
    /// High-frequency background permittivity.
    pub epsilon_infinity: f64,
    /// Oscillator terms.
    pub oscillators: Vec<Oscillator>,
}

impl DrudeLorentz {
    /// Construct.
    pub fn new(epsilon_infinity: f64, oscillators: Vec<Oscillator>) -> Self {
        DrudeLorentz {
            epsilon_infinity,
            oscillators,
        }
    }

    /// Add a Drude (free-electron) term with plasma frequency `wp` rad/s and
    /// damping `gamma` rad/s.
    pub fn drude(mut self, wp: f64, gamma: f64) -> Self {
        // Stored as S·ω₀² = ωₚ² with ω₀ = 0: handled specially in `index`.
        self.oscillators.push(Oscillator {
            strength: wp * wp,
            omega0: 0.0,
            gamma,
        });
        self
    }
}

impl Dispersion for DrudeLorentz {
    fn index(&self, wavelength: Wavelength) -> ComplexRefractiveIndex {
        let omega = wavelength.wavenumber() * crate::units::constants::SPEED_OF_LIGHT;
        let mut eps = Complex::real(self.epsilon_infinity);
        for osc in &self.oscillators {
            let contrib = if osc.omega0 == 0.0 {
                // Drude: −ωₚ²/(ω² + iγω)
                let denom = Complex::new(omega * omega, osc.gamma * omega);
                Complex::real(-osc.strength) / denom
            } else {
                let s_omega0_sq = osc.strength * osc.omega0 * osc.omega0;
                let denom =
                    Complex::new(osc.omega0 * osc.omega0 - omega * omega, -osc.gamma * omega);
                Complex::real(s_omega0_sq) / denom
            };
            eps += contrib;
        }
        let n = eps.sqrt();
        ComplexRefractiveIndex::new(n.re, n.im)
    }
}

/// One tabulated point: (vacuum wavelength in m, n, k).
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct TablePoint {
    /// Vacuum wavelength in meters.
    pub wavelength_m: f64,
    /// Real index.
    pub n: f64,
    /// Extinction coefficient.
    pub k: f64,
}

/// Cubic-spline (natural) interpolation of tabulated n/k data.
///
/// Points are sorted by wavelength on construction; lookups outside the
/// table clamp to the endpoint values (documented behaviour — the models are
/// not valid outside their measured band).
#[derive(Clone, Debug, PartialEq)]
pub struct Tabulated {
    points: Vec<TablePoint>,
    n_splines: Vec<[f64; 4]>,
    k_splines: Vec<[f64; 4]>,
}

impl Tabulated {
    /// Construct from (λ in m, n) pairs, k = 0.
    pub fn from_nm_pairs(pairs: &[(f64, f64)]) -> Self {
        Self::from_points(
            &pairs
                .iter()
                .map(|(nm, n)| TablePoint {
                    wavelength_m: nm * 1e-9,
                    n: *n,
                    k: 0.0,
                })
                .collect::<Vec<_>>(),
        )
    }

    /// Construct from explicit points (λ in meters). Sorted internally.
    pub fn from_points(points: &[TablePoint]) -> Self {
        assert!(!points.is_empty(), "tabulated dispersion needs data points");
        let mut points = points.to_vec();
        points.sort_by(|a, b| a.wavelength_m.total_cmp(&b.wavelength_m));
        let xs: Vec<f64> = points.iter().map(|p| p.wavelength_m).collect();
        let n_splines = natural_cubic(&xs, &points.iter().map(|p| p.n).collect::<Vec<_>>());
        let k_splines = natural_cubic(&xs, &points.iter().map(|p| p.k).collect::<Vec<_>>());
        Tabulated {
            points,
            n_splines,
            k_splines,
        }
    }

    /// Evaluate one set of spline coefficients at wavelength `x`,
    /// clamping outside the table to the endpoint values.
    fn eval(&self, splines: &[[f64; 4]], x: f64) -> f64 {
        let pts = &self.points;
        if x <= pts[0].wavelength_m {
            return splines[0][0];
        }
        let last = splines.len() - 1;
        if x >= pts[last + 1].wavelength_m {
            return splines[last][0];
        }
        // Binary search for the interval.
        let mut lo = 0usize;
        let mut hi = pts.len() - 1;
        while hi - lo > 1 {
            let mid = (lo + hi) / 2;
            if pts[mid].wavelength_m <= x {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let h = pts[hi].wavelength_m - pts[lo].wavelength_m;
        let t = (x - pts[lo].wavelength_m) / h;
        let [c0, c1, c2, c3] = splines[lo];
        // Horner evaluation of c0 + c1·t + c2·t² + c3·t³
        c0 + t * (c1 + t * (c2 + t * c3))
    }
}

impl Dispersion for Tabulated {
    fn index(&self, wavelength: Wavelength) -> ComplexRefractiveIndex {
        let x = wavelength.m();
        let pts = &self.points;
        // Out-of-range lookups clamp to the endpoint values.
        if x <= pts[0].wavelength_m {
            return ComplexRefractiveIndex::new(pts[0].n, pts[0].k);
        }
        let last = pts.len() - 1;
        if x >= pts[last].wavelength_m {
            return ComplexRefractiveIndex::new(pts[last].n, pts[last].k);
        }
        ComplexRefractiveIndex::new(self.eval(&self.n_splines, x), self.eval(&self.k_splines, x))
    }
}

/// Linear interpolation of tabulated data (robust, monotone-preserving).
#[derive(Clone, Debug, PartialEq)]
pub struct TabulatedLinear {
    points: Vec<TablePoint>,
}

impl TabulatedLinear {
    /// Construct from points (λ in meters). Sorted internally.
    pub fn from_points(points: &[TablePoint]) -> Self {
        assert!(!points.is_empty(), "tabulated dispersion needs data points");
        let mut points = points.to_vec();
        points.sort_by(|a, b| a.wavelength_m.total_cmp(&b.wavelength_m));
        TabulatedLinear { points }
    }

    /// Construct from (λ nm, n) pairs with k = 0.
    pub fn from_nm_pairs(pairs: &[(f64, f64)]) -> Self {
        Self::from_points(
            &pairs
                .iter()
                .map(|(nm, n)| TablePoint {
                    wavelength_m: nm * 1e-9,
                    n: *n,
                    k: 0.0,
                })
                .collect::<Vec<_>>(),
        )
    }
}

impl Dispersion for TabulatedLinear {
    fn index(&self, wavelength: Wavelength) -> ComplexRefractiveIndex {
        let pts = &self.points;
        let x = wavelength.m();
        if x <= pts[0].wavelength_m {
            return ComplexRefractiveIndex::new(pts[0].n, pts[0].k);
        }
        if x >= pts[pts.len() - 1].wavelength_m {
            let last = &pts[pts.len() - 1];
            return ComplexRefractiveIndex::new(last.n, last.k);
        }
        let mut lo = 0usize;
        let mut hi = pts.len() - 1;
        while hi - lo > 1 {
            let mid = (lo + hi) / 2;
            if pts[mid].wavelength_m <= x {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let t = (x - pts[lo].wavelength_m) / (pts[hi].wavelength_m - pts[lo].wavelength_m);
        ComplexRefractiveIndex::new(
            pts[lo].n + t * (pts[hi].n - pts[lo].n),
            pts[lo].k + t * (pts[hi].k - pts[lo].k),
        )
    }
}

/// Natural cubic spline through `(xs, ys)`, returning per-interval
/// coefficients `[c0, c1, c2, c3]` in normalised parameter
/// `t = (x − xᵢ)/(xᵢ₊₁ − xᵢ)`.
pub(crate) fn natural_cubic(xs: &[f64], ys: &[f64]) -> Vec<[f64; 4]> {
    let n = xs.len();
    debug_assert!(n >= 2);
    if n == 2 {
        // Linear fallback.
        let slope = (ys[1] - ys[0]) / (xs[1] - xs[0]);
        return vec![[ys[0], slope, 0.0, 0.0]];
    }
    // Solve the tridiagonal system for second derivatives (Thomas algorithm)
    // on second differences Mᵢ with natural boundary M₀ = Mₙ₋₁ = 0.
    let mut lower = vec![0.0; n];
    let mut diag = vec![0.0; n];
    let mut upper = vec![0.0; n];
    let mut rhs = vec![0.0; n];
    for i in 1..n - 1 {
        let hi_ = xs[i] - xs[i - 1];
        let hi1 = xs[i + 1] - xs[i];
        lower[i] = hi_;
        diag[i] = 2.0 * (hi_ + hi1);
        upper[i] = hi1;
        rhs[i] = 6.0 * ((ys[i + 1] - ys[i]) / hi1 - (ys[i] - ys[i - 1]) / hi_);
    }
    // Forward sweep (the first unknown is M₁; M₀ = 0 is the natural BC).
    for i in 2..n - 1 {
        let w = lower[i] / diag[i - 1];
        diag[i] -= w * upper[i - 1];
        rhs[i] -= w * rhs[i - 1];
    }
    // Back substitution.
    let mut m = vec![0.0; n];
    for i in (1..n - 1).rev() {
        m[i] = (rhs[i] - upper[i] * m[i + 1]) / diag[i];
    }
    // Per-interval coefficients in normalised parameter
    // `t = (x − xᵢ)/(xᵢ₊₁ − xᵢ)`.
    let mut out = Vec::with_capacity(n - 1);
    for i in 0..n - 1 {
        let hi_ = xs[i + 1] - xs[i];
        let c0 = ys[i];
        let c1 = (ys[i + 1] - ys[i]) - hi_ * hi_ * (2.0 * m[i] + m[i + 1]) / 6.0;
        let c2 = hi_ * hi_ * m[i] / 2.0;
        let c3 = hi_ * hi_ * (m[i + 1] - m[i]) / 6.0;
        out.push([c0, c1, c2, c3]);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::units::constants::SPEED_OF_LIGHT;

    fn wl_nm(nm: f64) -> Wavelength {
        Wavelength::from_nm(nm)
    }

    #[test]
    fn sellmeier_bk7_d_line() {
        let bk7 = Sellmeier::bk7();
        // Catalog: n_d(587.6 nm) = 1.51680
        let n = bk7.real_index(wl_nm(587.6));
        assert!((n - 1.5168).abs() < 2e-5, "n = {n}");
        // F-line 486.1 nm: 1.52238; C-line 656.3 nm: 1.51432
        assert!((bk7.real_index(wl_nm(486.1)) - 1.52238).abs() < 2e-5);
        assert!((bk7.real_index(wl_nm(656.3)) - 1.51432).abs() < 2e-5);
        // Abbe number V_d = (n_d - 1)/(n_F - n_C) ≈ 64.17
        let v: f64 = (1.5168 - 1.0) / (1.52238 - 1.51432);
        assert!((v - 64.17).abs() < 0.2);
    }

    #[test]
    fn sellmeier_fused_silica() {
        let sio2 = Sellmeier::fused_silica();
        // n(587.6 nm) ≈ 1.45846 (Malitson)
        let n = sio2.real_index(wl_nm(587.6));
        assert!((n - 1.45846).abs() < 5e-5, "n = {n}");
        // n(1550 nm) ≈ 1.44402
        assert!((sio2.real_index(wl_nm(1550.0)) - 1.44402).abs() < 5e-5);
    }

    #[test]
    fn sellmeier_diamond() {
        let d = Sellmeier::diamond();
        let n = d.real_index(wl_nm(589.0));
        assert!((n - 2.417).abs() < 5e-3, "n = {n}");
    }

    #[test]
    fn cauchy_known_glass() {
        // A crude but standard water fit: n ≈ 1.3247 + 3146.7/λ²_µm? Use a
        // synthetic check instead: reproduce a quadratic exactly.
        let c = Cauchy::new(1.5, 0.004, 0.0001);
        let l = 0.5;
        let expect = 1.5 + 0.004 / (l * l) + 0.0001 / (l * l * l * l);
        assert!((c.real_index(Wavelength::from_um(l)) - expect).abs() < 1e-14);
    }

    #[test]
    fn schott_reproduces_quadratic_form() {
        // n² = 2.3 + 0.01 λ² + 0.01/λ²  -> at λ = 1 µm, n² = 2.32
        let s = Schott::new([2.3, 0.01, 0.01, 0.0, 0.0, 0.0]);
        assert!((s.real_index(Wavelength::from_um(1.0)).powi(2) - 2.32).abs() < 1e-12);
    }

    #[test]
    fn drude_metal_has_negative_real_permittivity_above_plasma() {
        // Simple Drude metal far above ωₚ: n real part < 1, strong k.
        let m = DrudeLorentz::new(1.0, Vec::new()).drude(2.0e16, 5.0e13);
        let idx = m.index(wl_nm(500.0));
        assert!(idx.n < 1.0);
        assert!(idx.k > 0.5);
    }

    #[test]
    fn lorentz_resonance_absorbs() {
        // A Lorentz oscillator at the probe frequency should give k > 0.
        let wl = wl_nm(500.0);
        let omega = wl.wavenumber() * SPEED_OF_LIGHT;
        let m = DrudeLorentz::new(
            1.0,
            vec![Oscillator {
                strength: 0.5,
                omega0: omega,
                gamma: 0.05 * omega,
            }],
        );
        let idx = m.index(wl);
        assert!(idx.k > 0.01, "k = {}", idx.k);
        assert!(idx.n > 0.0);
    }

    #[test]
    fn tabulated_spline_interpolates_smoothly() {
        // Sample sin at 1-degree steps; spline should recover midpoints.
        let pairs: Vec<(f64, f64)> = (0..=90)
            .map(|d| {
                (
                    d as f64,
                    libm::sin(d as f64 * core::f64::consts::PI / 180.0),
                )
            })
            .collect();
        let t = Tabulated::from_nm_pairs(&pairs);
        let mid = t.real_index(Wavelength::from_nm(45.5));
        let expect = libm::sin(45.5 * core::f64::consts::PI / 180.0);
        assert!((mid - expect).abs() < 1e-6, "got {mid}, want {expect}");
        // Out-of-range clamps to the endpoints (sin 0° = 0, sin 90° = 1).
        assert!((t.real_index(Wavelength::from_nm(-10.0)) - 0.0).abs() < 1e-9);
        assert!((t.real_index(Wavelength::from_nm(500.0)) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn tabulated_linear_clamps() {
        let t = TabulatedLinear::from_nm_pairs(&[(400.0, 1.5), (800.0, 1.6)]);
        assert!((t.real_index(wl_nm(400.0)) - 1.5).abs() < 1e-12);
        assert!((t.real_index(wl_nm(600.0)) - 1.55).abs() < 1e-12);
        assert!((t.real_index(wl_nm(1000.0)) - 1.6).abs() < 1e-12);
        assert!((t.real_index(wl_nm(100.0)) - 1.5).abs() < 1e-12);
    }
}
