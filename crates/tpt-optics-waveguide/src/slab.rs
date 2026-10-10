//! Slab waveguide mode solving and planar-stack transfer matrices.
//!
//! # Symmetric slab
//!
//! Exact transcendental mode equations for a symmetric high-index core of
//! thickness `d` in a cladding `n_clad` (weak-guidance uses the same TE/TM
//! equations with the polarization factor). Modes are found by bracketing
//! and bisection on `u ∈ (0, π/2 + mπ)` where
//! `u = (d/2)·√(k₀²n_core² − β²)`, `w = (d/2)·√(β² − k₀²n_clad²)`, with the
//! invariance `u² + w² = V²`, `V = (d/2)·k₀·√(n_core² − n_clad²)`.
//!
//! # Planar stacks
//!
//! Characteristic-matrix (transfer matrix) reflectance/transmittance at
//! normal incidence for multilayer dielectric stacks.

use alloc::vec::Vec;

use libm::sqrt;
use tpt_optics_core::Wavelength;

/// Polarization of a slab mode.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum SlabPolarization {
    /// Transverse electric (s).
    Te,
    /// Transverse magnetic (p).
    Tm,
}

/// One solved slab mode.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct SlabMode {
    /// Effective index β/k₀.
    pub n_eff: f64,
    /// Mode order (0 = fundamental, even).
    pub order: usize,
    /// Polarization.
    pub polarization: SlabPolarization,
}

/// Symmetric slab waveguide parameters.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct SlabLayer {
    /// Core index.
    pub n_core: f64,
    /// Cladding index (both sides).
    pub n_clad: f64,
    /// Core thickness (m).
    pub thickness: f64,
}

/// All guided modes of the symmetric slab at `wl`, ordered by mode number.
pub fn slab_modes(slab: &SlabLayer, wl: Wavelength, pol: SlabPolarization) -> Vec<SlabMode> {
    let k0 = wl.wavenumber();
    let v = 0.5 * slab.thickness * k0 * sqrt(slab.n_core * slab.n_core - slab.n_clad * slab.n_clad);
    if v <= 0.0 {
        return Vec::new();
    }
    // Mode m exists when V > m·π/2.
    let n_modes = (libm::floor(v / core::f64::consts::FRAC_PI_2) as usize) + 1;
    let idx2 = slab.n_core * slab.n_core / (slab.n_clad * slab.n_clad);
    let mut modes = Vec::with_capacity(n_modes);
    for m in 0..n_modes {
        // u ∈ (mπ/2, (m+1)π/2]; root of F(u) with w = √(V²−u²).
        let f = |u: f64| -> f64 {
            let w = sqrt((v * v - u * u).max(0.0));
            let tan = libm::tan(u);
            match (m % 2 == 0, pol) {
                (true, SlabPolarization::Te) => u * tan - w,
                (false, SlabPolarization::Te) => u / tan + w,
                (true, SlabPolarization::Tm) => u * tan - idx2 * w,
                (false, SlabPolarization::Tm) => u / tan + idx2 * w,
            }
        };
        let a = m as f64 * core::f64::consts::FRAC_PI_2 + 1e-9;
        let b = ((m + 1) as f64 * core::f64::consts::FRAC_PI_2).min(v) - 1e-9;
        if b <= a {
            break;
        }
        let (fa, fb) = (f(a), f(b));
        if fa * fb > 0.0 {
            continue; // no root bracketed (mode at cutoff edge)
        }
        // Bisection, 80 iterations → full double precision.
        let (mut lo, mut hi) = (a, b);
        for _ in 0..80 {
            let mid = 0.5 * (lo + hi);
            if f(lo) * f(mid) <= 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        let u = 0.5 * (lo + hi);
        let ku = 2.0 * u / slab.thickness;
        let beta = sqrt(k0 * k0 * slab.n_core * slab.n_core - ku * ku);
        if beta > k0 * slab.n_clad && beta.is_finite() {
            modes.push(SlabMode {
                n_eff: beta / k0,
                order: m,
                polarization: pol,
            });
        }
    }
    modes
}

/// One layer of a planar stack (normal incidence).
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct StackLayer {
    /// Refractive index.
    pub index: f64,
    /// Physical thickness (m).
    pub thickness: f64,
}

/// Normal-incidence reflectance and transmittance of a dielectric stack
/// between a semi-infinite incident medium `n0` and substrate `ns`, via the
/// characteristic transfer matrix (MacLeod convention):
/// `M = M₁···M_N`, `[B; C] = M [1; ηs]`,
/// `r = (η₀B − C)/(η₀B + C)`, `t = 2η₀/(η₀B + C)`.
pub fn stack_r_t(layers: &[StackLayer], wl: Wavelength, n0: f64, ns: f64) -> (f64, f64) {
    // M = M_1 · M_2 · ··· M_N (matrix product in entry order).
    let mut m = [
        [Complex2::new(1.0, 0.0), Complex2::new(0.0, 0.0)],
        [Complex2::new(0.0, 0.0), Complex2::new(1.0, 0.0)],
    ];
    for layer in layers {
        let eta = layer.index;
        let delta = core::f64::consts::TAU * eta * layer.thickness / wl.m();
        let (cd, sd) = (libm::cos(delta), libm::sin(delta));
        let mj = [
            [Complex2::new(cd, 0.0), Complex2::new(0.0, sd / eta)],
            [Complex2::new(0.0, sd * eta), Complex2::new(cd, 0.0)],
        ];
        m = mat_mul(&mj, &m);
    }
    let b = m[0][0] + m[0][1].scale(ns);
    let c = m[1][0] + m[1][1].scale(ns);
    let denom = b.scale(n0) + c;
    let r = (b.scale(n0) - c) / denom;
    let t = Complex2::new(2.0 * n0, 0.0) / denom;
    let big_r = r.re * r.re + r.im * r.im;
    let big_t = (t.re * t.re + t.im * t.im) * ns / n0;
    (big_r, big_t)
}

fn mat_mul(a: &[[Complex2; 2]; 2], b: &[[Complex2; 2]; 2]) -> [[Complex2; 2]; 2] {
    let mut out = [[Complex2::new(0.0, 0.0); 2]; 2];
    for (i, row) in out.iter_mut().enumerate() {
        for (j, v) in row.iter_mut().enumerate() {
            *v = a[i][0] * b[0][j] + a[i][1] * b[1][j];
        }
    }
    out
}

/// Minimal complex pair (re, im) sufficient for the stack math.
#[derive(Copy, Clone, Debug, PartialEq)]
pub(crate) struct Complex2 {
    pub re: f64,
    pub im: f64,
}

impl Complex2 {
    pub const fn new(re: f64, im: f64) -> Self {
        Complex2 { re, im }
    }
    pub fn mul(self, o: Self) -> Self {
        Complex2::new(
            self.re * o.re - self.im * o.im,
            self.re * o.im + self.im * o.re,
        )
    }
    pub fn add(self, o: Self) -> Self {
        Complex2::new(self.re + o.re, self.im + o.im)
    }
    pub fn scale(self, s: f64) -> Self {
        Complex2::new(self.re * s, self.im * s)
    }
    pub fn cdiv(self, o: Self) -> Self {
        let d = o.re * o.re + o.im * o.im;
        Complex2::new(
            (self.re * o.re + self.im * o.im) / d,
            (self.im * o.re - self.re * o.im) / d,
        )
    }
}

impl core::ops::Add for Complex2 {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        self.add(o)
    }
}
impl core::ops::Sub for Complex2 {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        Complex2::new(self.re - o.re, self.im - o.im)
    }
}
impl core::ops::Div for Complex2 {
    type Output = Self;
    fn div(self, o: Self) -> Self {
        self.cdiv(o)
    }
}
impl core::ops::Mul for Complex2 {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        self.mul(o)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    fn wl(nm: f64) -> Wavelength {
        Wavelength::from_nm(nm)
    }

    #[test]
    fn single_mode_slab_has_one_te_mode() {
        // Silicon slab 220 nm in SiO2 at 1550 nm: V ≈ 1.86 < π/2 → 1 mode.
        let slab = SlabLayer {
            n_core: 3.48,
            n_clad: 1.444,
            thickness: 220e-9,
        };
        let modes = slab_modes(&slab, wl(1550.0), SlabPolarization::Te);
        assert_eq!(modes.len(), 1, "{modes:?}");
        let neff = modes[0].n_eff;
        // SOI slab TE0 effective index ≈ 2.85 (literature ballpark).
        assert!((neff - 2.85).abs() < 0.15, "n_eff = {neff}");
        assert!(neff > slab.n_clad && neff < slab.n_core);
    }

    #[test]
    fn mode_count_increases_with_thickness() {
        let mut slab = SlabLayer {
            n_core: 3.48,
            n_clad: 1.444,
            thickness: 220e-9,
        };
        let m0 = slab_modes(&slab, wl(1550.0), SlabPolarization::Te).len();
        slab.thickness = 3.0e-6; // thick slab → many modes
        let m1 = slab_modes(&slab, wl(1550.0), SlabPolarization::Te).len();
        assert!(m1 > m0);
        assert!(m1 >= 5);
    }

    #[test]
    fn te_tm_neff_ordered() {
        // For the same slab, TM n_eff ≤ TE n_eff (weak-guidance factor).
        let slab = SlabLayer {
            n_core: 3.48,
            n_clad: 1.444,
            thickness: 500e-9,
        };
        let te = slab_modes(&slab, wl(1550.0), SlabPolarization::Te);
        let tm = slab_modes(&slab, wl(1550.0), SlabPolarization::Tm);
        assert_eq!(te.len(), tm.len());
        assert!(te[0].n_eff >= tm[0].n_eff);
    }

    #[test]
    fn quarter_wave_stack_is_high_reflectivity() {
        // 6 pairs of SiO2/TiO2 quarter-wave layers at 1064 nm → R > 0.95.
        let wl0 = wl(1064.0);
        let d = |n: f64| wl0.m() / (4.0 * n);
        let mut layers = Vec::new();
        for _ in 0..6 {
            layers.push(StackLayer {
                index: 1.45,
                thickness: d(1.45),
            });
            layers.push(StackLayer {
                index: 2.35,
                thickness: d(2.35),
            });
        }
        let (r, t) = stack_r_t(&layers, wl0, 1.0, 1.0);
        assert!(r > 0.95, "R = {r}");
        assert!(close(r + t, 1.0, 1e-9), "R+T = {}", r + t);
        // And it is transmissive away from the band.
        let (r2, _) = stack_r_t(&layers, wl(500.0), 1.0, 1.0);
        assert!(r2 < 0.9);
    }

    #[test]
    fn single_interface_stack_matches_fresnel() {
        // No layers: bare air-glass interface. R = 0.04 for n = 1.5.
        let (r, t) = stack_r_t(&[], wl(550.0), 1.0, 1.5);
        assert!(close(r, 0.04, 1e-12), "R = {r}");
        assert!(close(t, 0.96, 1e-12));
    }
}
