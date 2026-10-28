//! Point defects in 2-D photonic crystals via supercell plane-wave expansion.
//!
//! A defect is modelled by tiling an `M × M` supercell (M odd) of the square
//! rod lattice and replacing the rod at the cell centre with a rod of
//! different radius and permittivity (a missing rod is radius 0). Bloch
//! states of the supercell are eigenstates of a plane-wave problem evaluated
//! on the supercell reciprocal lattice `G = (2π/M)(gx, gy)`. Frequencies are
//! in units of `c/a`, with `a` the lattice constant of the underlying crystal,
//! so they can be compared directly with the bulk band structure.
//!
//! Two polarisations are provided:
//! - [`supercell_tm_freqs`] — E-polarisation (E_z). Solves
//!   `|k+G|² E(G) = (ω/c)² Σ ε(G−G′) E(G′)`, a generalised problem with a
//!   positive-definite ε matrix. The ε Fourier coefficients converge quickly
//!   for rods, so this is the default for defect work.
//! - [`supercell_te_freqs`] — H-polarisation (H_z), with 1/ε coefficients.
//!   These converge slowly at high contrast; kept for comparison with
//!   [`crate::phc::pwem_te_bands`].
//!
//! Defect states appear as supercell bands that have no counterpart in the
//! bulk spectrum at the same `k`.
//!
//! Accuracy limits:
//! - Plane-wave truncation (`n_g` per axis). A supercell of size `M` has
//!   reciprocal spacing `2π/M`, so it needs about `M` times more plane waves
//!   per axis than the bulk solver to reach the same `G` range.
//! - The supercell is periodic, so defects interact with their images at
//!   distance `M·a`. Use `M ≥ 5` when the defect band must be isolated from
//!   its images.

use alloc::vec;
use alloc::vec::Vec;
use core::f64::consts::{PI, TAU};

use libm::{cos, sqrt};
use tpt_optics_core::bessel_j;

use crate::eigen::symmetric_eigenvalues;

/// A dielectric rod placed at an integer cell position of the supercell.
/// Positions are centred: for `M = 5` the cell indices run over `-2..=2`.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Rod {
    /// Cell index along x (centred).
    pub x: i32,
    /// Cell index along y (centred).
    pub y: i32,
    /// Rod radius over lattice constant, `r/a`. Zero means no rod.
    pub radius: f64,
    /// Rod relative permittivity.
    pub epsilon: f64,
}

/// Rods of an `M × M` supercell with every site filled by a rod of radius
/// `r_over_a` and permittivity `epsilon`, except the centre site, which holds
/// a defect rod. Pass `defect_radius = 0` for a missing rod.
pub fn point_defect_rods(
    m: usize,
    epsilon: f64,
    r_over_a: f64,
    defect_epsilon: f64,
    defect_radius: f64,
) -> Vec<Rod> {
    let h = (m as i32 - 1) / 2;
    let mut rods = Vec::with_capacity(m * m);
    for x in -h..=h {
        for y in -h..=h {
            rods.push(if x == 0 && y == 0 {
                Rod {
                    x,
                    y,
                    radius: defect_radius,
                    epsilon: defect_epsilon,
                }
            } else {
                Rod {
                    x,
                    y,
                    radius: r_over_a,
                    epsilon,
                }
            });
        }
    }
    rods
}

/// Fourier table of a supercell material function, indexed by integer
/// differences `(dx, dy)` with `|dx|, |dy| ≤ d_max`.
///
/// `base` is the background value and `weight(rod)` the contrast of each
/// rod relative to the background. The mean is `base + Σ weight·πρ²/M²`,
/// and the off-diagonal term is `Σ weight·e^{−iG·R}·(2πρ/G)J₁(Gρ)/M²`, which
/// is real because the rod set is inversion-symmetric.
struct FourierTable {
    d_max: i32,
    span: usize,
    values: Vec<f64>,
}

impl FourierTable {
    fn new(m: usize, rods: &[Rod], d_max: i32, base: f64, weight: impl Fn(&Rod) -> f64) -> Self {
        let mf = m as f64;
        let span = (2 * d_max + 1) as usize;
        let mut values = vec![0.0; span * span];
        for dx in -d_max..=d_max {
            for dy in -d_max..=d_max {
                let idx = ((dx + d_max) as usize) * span + (dy + d_max) as usize;
                values[idx] = if dx == 0 && dy == 0 {
                    let total: f64 = rods.iter().map(|rod| weight(rod) * PI * rod.radius * rod.radius).sum();
                    base + total / (mf * mf)
                } else {
                    let g = TAU * sqrt((dx * dx + dy * dy) as f64) / mf;
                    let mut sum = 0.0;
                    for rod in rods {
                        if rod.radius == 0.0 {
                            continue;
                        }
                        let phase = TAU * (dx * rod.x + dy * rod.y) as f64 / mf;
                        let gr = g * rod.radius;
                        // (1/A)∫_disk e^{-iG·r} dA = (2πρ/G)·J₁(Gρ).
                        let ft = 2.0 * PI * rod.radius * bessel_j(1, gr) / g;
                        sum += weight(rod) * cos(phase) * ft;
                    }
                    sum / (mf * mf)
                };
            }
        }
        Self { d_max, span, values }
    }

    fn get(&self, dx: i32, dy: i32) -> f64 {
        self.values[((dx + self.d_max) as usize) * self.span + (dy + self.d_max) as usize]
    }
}

/// Supercell reciprocal vectors `(gx, gy)` in units of `2π/M`, for a basis
/// of `n_g` per axis, and the largest index difference they can produce.
fn basis(n_g: usize) -> (Vec<(i32, i32)>, i32) {
    let ng = n_g as i32;
    let g_max = ng / 2;
    let gs = (0..ng)
        .flat_map(|ix| (0..ng).map(move |iy| (ix - g_max, iy - g_max)))
        .collect();
    (gs, 2 * g_max)
}

/// TE (H_z) supercell frequencies at Bloch vector `k` (radians per `a`), in
/// units of `c/a` (`ωa/2πc`), ascending. Returns `n_g²` values.
///
/// * `m` — supercell size in lattice constants (odd),
/// * `rods` — rod set, typically from [`point_defect_rods`]; background is air,
/// * `n_g` — plane waves per axis (odd).
pub fn supercell_te_freqs(m: usize, rods: &[Rod], n_g: usize, k: (f64, f64)) -> Vec<f64> {
    let (gs, d_max) = basis(n_g);
    let eta = FourierTable::new(m, rods, d_max, 1.0, |rod| 1.0 / rod.epsilon - 1.0);
    let n = gs.len();
    let mf = m as f64;
    let (kx, ky) = k;
    let mut h = vec![0.0; n * n];
    for (i, &(gxi, gyi)) in gs.iter().enumerate() {
        let (px, py) = (kx + TAU * gxi as f64 / mf, ky + TAU * gyi as f64 / mf);
        for (j, &(gxj, gyj)) in gs.iter().enumerate() {
            let (qx, qy) = (kx + TAU * gxj as f64 / mf, ky + TAU * gyj as f64 / mf);
            h[i * n + j] = (px * qx + py * qy) * eta.get(gxi - gxj, gyi - gyj);
        }
    }
    let w = symmetric_eigenvalues(&mut h, n);
    w.iter().map(|lam| sqrt(lam.max(0.0)) / TAU).collect()
}

/// TM (E_z) supercell frequencies at Bloch vector `k`, in units of `c/a`,
/// ascending. Returns `n_g²` values. See the module docs for the formulation.
///
/// Solves `A x = λ B x` with `A = diag(|k+G|²)` and `B = ε(G−G′)`. Reducing
/// with the Cholesky factor `B = L Lᵀ` gives the standard problem
/// `(L⁻¹ A L⁻ᵀ) y = λ y`.
pub fn supercell_tm_freqs(m: usize, rods: &[Rod], n_g: usize, k: (f64, f64)) -> Vec<f64> {
    let (gs, d_max) = basis(n_g);
    let eps = FourierTable::new(m, rods, d_max, 1.0, |rod| rod.epsilon - 1.0);
    let n = gs.len();
    let mf = m as f64;
    let (kx, ky) = k;
    let a_diag: Vec<f64> = gs
        .iter()
        .map(|&(gx, gy)| {
            let (px, py) = (kx + TAU * gx as f64 / mf, ky + TAU * gy as f64 / mf);
            px * px + py * py
        })
        .collect();
    let mut b = vec![0.0; n * n];
    for (i, &(gxi, gyi)) in gs.iter().enumerate() {
        for (j, &(gxj, gyj)) in gs.iter().enumerate() {
            b[i * n + j] = eps.get(gxi - gxj, gyi - gyj);
        }
    }
    let lower_inv = cholesky_lower_inverse(&b, n);
    // C = X·A·Xᵀ with X = L⁻¹ lower triangular.
    let mut c = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..=i {
            let mut s = 0.0;
            for kk in 0..=j {
                s += lower_inv[i * n + kk] * a_diag[kk] * lower_inv[j * n + kk];
            }
            c[i * n + j] = s;
            c[j * n + i] = s;
        }
    }
    let w = symmetric_eigenvalues(&mut c, n);
    w.iter().map(|lam| sqrt(lam.max(0.0)) / TAU).collect()
}

/// Cholesky factor `L` of the symmetric positive-definite `b` (row-major,
/// overwritten), returned as its inverse `X = L⁻¹` (lower triangular, dense
/// storage).
fn cholesky_lower_inverse(b: &[f64], n: usize) -> Vec<f64> {
    let mut l = vec![0.0; n * n];
    for j in 0..n {
        let mut s = b[j * n + j];
        for k in 0..j {
            s -= l[j * n + k] * l[j * n + k];
        }
        let d = sqrt(s.max(f64::MIN_POSITIVE));
        l[j * n + j] = d;
        for i in (j + 1)..n {
            let mut t = b[i * n + j];
            for k in 0..j {
                t -= l[i * n + k] * l[j * n + k];
            }
            l[i * n + j] = t / d;
        }
    }
    // Forward substitution per column: L·X = I.
    let mut x = vec![0.0; n * n];
    for col in 0..n {
        x[col * n + col] = 1.0 / l[col * n + col];
        for i in (col + 1)..n {
            let mut s = 0.0;
            for k in col..i {
                s += l[i * n + k] * x[k * n + col];
            }
            x[i * n + col] = -s / l[i * n + i];
        }
    }
    x
}

/// Widest frequency window inside `[f_lo, f_hi]` that contains no bulk state
/// of the uniform crystal, found by sampling the supercell spectrum over the
/// full supercell Brillouin zone on an `(n_k + 1)²` grid.
///
/// This is a sampling estimate. A coarse grid can miss band extrema and
/// report a gap that is too wide, and the result also depends on the
/// plane-wave cutoff. Check it with a denser `n_k` and a larger `n_g` before
/// relying on it. Returns `(lower, upper)` in units of `c/a`.
pub fn bulk_gap_from_supercell(
    m: usize,
    epsilon: f64,
    r_over_a: f64,
    n_g: usize,
    n_k: usize,
    f_lo: f64,
    f_hi: f64,
    tm: bool,
) -> Option<(f64, f64)> {
    let rods = point_defect_rods(m, epsilon, r_over_a, epsilon, r_over_a);
    let kmax = PI / m as f64;
    let mut all: Vec<f64> = Vec::new();
    for i in 0..=n_k {
        for j in 0..=n_k {
            let k = (kmax * i as f64 / n_k as f64, kmax * j as f64 / n_k as f64);
            let freqs = if tm {
                supercell_tm_freqs(m, &rods, n_g, k)
            } else {
                supercell_te_freqs(m, &rods, n_g, k)
            };
            all.extend(freqs.into_iter().filter(|f| *f >= f_lo && *f <= f_hi));
        }
    }
    all.sort_by(f64::total_cmp);
    let mut best: Option<(f64, f64)> = None;
    let mut prev = f_lo;
    for &f in all.iter().chain(core::iter::once(&f_hi)) {
        if f - prev > best.map_or(0.0, |(lo, hi)| hi - lo) {
            best = Some((prev, f));
        }
        prev = f;
    }
    best.filter(|(lo, hi)| hi > lo)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_supercell_gives_folded_plane_waves() {
        // ε = 1 everywhere: the spectrum at k = 0 is |n|/M for n ∈ Z².
        let m = 3;
        let rods = point_defect_rods(m, 1.0, 0.0, 1.0, 0.0);
        for solve in [supercell_te_freqs, supercell_tm_freqs] {
            let f = solve(m, &rods, 7, (0.0, 0.0));
            assert!(f[0].abs() < 1e-9, "{f:?}");
            for &v in &f[1..5] {
                assert!((v - 1.0 / m as f64).abs() < 1e-6, "{f:?}");
            }
        }
    }

    #[test]
    fn missing_rod_defect_differs_from_uniform() {
        // Removing the centre rod must change the spectrum: the uniform and
        // defect supercells share the same bulk, so any difference is the defect.
        let (eps, r) = (11.56, 0.2);
        let uniform = point_defect_rods(5, eps, r, eps, r);
        let missing = point_defect_rods(5, eps, r, 1.0, 0.0);
        let fu = supercell_tm_freqs(5, &uniform, 5, (0.0, 0.0));
        let fm = supercell_tm_freqs(5, &missing, 5, (0.0, 0.0));
        let diff: f64 = fu.iter().zip(&fm).map(|(a, b)| (a - b).abs()).sum();
        assert!(diff > 1e-3, "defect had no effect: {diff}");
    }

    #[test]
    fn single_cell_tm_matches_single_cell_eigenproblem() {
        // M = 1 with uniform rods: the TM problem must reproduce the bulk
        // E_z spectrum at Γ, which the TE and TM solvers share at ε = 1.
        let rods = point_defect_rods(1, 1.0, 0.0, 1.0, 0.0);
        let f = supercell_tm_freqs(1, &rods, 5, (0.0, 0.0));
        assert!(f[0].abs() < 1e-9 && (f[1] - 1.0).abs() < 1e-9, "{f:?}");
    }
}
