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
//! Convergence (checked for ε = 8.9, r/a = 0.2, M = 5):
//! - TM: the bulk Γ/X/M levels agree to ~1e-3 at `n_g = 13`. The full-zone gap
//!   is (0.324, 0.447) c/a at `n_g = 17–21`. A missing-rod defect gives a
//!   level at 0.393 c/a, converged to ~0.3% for `n_g = 17–29`, and its
//!   dispersion across the supercell zone is ~0.004 c/a.
//! - TE: the 1/ε formulation converges slowly at high contrast. The bulk gap
//!   moved by more than its own width between cutoffs, so TE defect levels are
//!   not quantitative.
//! - High-permittivity defects (ε = 40) are not converged at `n_g ≤ 29`.
//!   Their levels move by several hundredths between cutoffs.
//! - The supercell is periodic, so defects interact with their images at
//!   distance `M·a`. For the missing-rod level, the 7×7 cell gives 0.3996 c/a
//!   at `n_g = 21` against 0.3943 for the 5×5 cell, a 1.3% shift. Quote the
//!   isolated-defect level as 0.39–0.40 c/a until larger cells are converged.

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
                    let total: f64 = rods
                        .iter()
                        .map(|rod| weight(rod) * PI * rod.radius * rod.radius)
                        .sum();
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
        Self {
            d_max,
            span,
            values,
        }
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
///
/// Returns an error if the supercell description is invalid (see
/// [`validate`]).
pub fn supercell_te_freqs(
    m: usize,
    rods: &[Rod],
    n_g: usize,
    k: (f64, f64),
) -> Result<Vec<f64>, SupercellError> {
    validate(m, rods, n_g)?;
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
    Ok(w.iter().map(|lam| sqrt(lam.max(0.0)) / TAU).collect())
}

/// TM (E_z) supercell frequencies at Bloch vector `k`, in units of `c/a`,
/// ascending. Returns `n_g²` values. See the module docs for the formulation.
///
/// Solves `A x = λ B x` with `A = diag(|k+G|²)` and `B = ε(G−G′)`. Reducing
/// with the Cholesky factor `B = L Lᵀ` gives the standard problem
/// `(L⁻¹ A L⁻ᵀ) y = λ y`.
///
/// Returns an error if the supercell description is invalid, or if `B` is
/// not positive definite.
pub fn supercell_tm_freqs(
    m: usize,
    rods: &[Rod],
    n_g: usize,
    k: (f64, f64),
) -> Result<Vec<f64>, SupercellError> {
    validate(m, rods, n_g)?;
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
    let lower_inv = cholesky_lower_inverse(&b, n)?;
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
    Ok(w.iter().map(|lam| sqrt(lam.max(0.0)) / TAU).collect())
}

/// Cholesky factor `L` of the symmetric positive-definite `b` (row-major),
/// returned as its inverse `X = L⁻¹` (lower triangular, dense storage).
///
/// A non-positive pivot means `b` is not positive definite. That is reported
/// as an error instead of being clamped to a tiny positive value.
fn cholesky_lower_inverse(b: &[f64], n: usize) -> Result<Vec<f64>, SupercellError> {
    let mut l = vec![0.0; n * n];
    for j in 0..n {
        let mut s = b[j * n + j];
        for k in 0..j {
            s -= l[j * n + k] * l[j * n + k];
        }
        if !(s > 0.0 && s.is_finite()) {
            return Err(SupercellError::NotPositiveDefinite);
        }
        let d = sqrt(s);
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
    Ok(x)
}

/// Errors from the supercell sampling helpers.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum SupercellError {
    /// An argument is out of range or not finite. The message names it.
    InvalidParameter(&'static str),
    /// The eigensolver returned a non-finite frequency. Such values are never
    /// dropped silently, because that would hide a broken spectrum.
    NonFiniteFrequency,
    /// The permittivity matrix lost positive-definiteness in the Cholesky
    /// factorisation. The problem is then not a valid generalised eigenproblem.
    NotPositiveDefinite,
}

impl core::fmt::Display for SupercellError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidParameter(msg) => write!(f, "invalid supercell parameter: {msg}"),
            Self::NonFiniteFrequency => {
                write!(f, "supercell eigensolver gave a non-finite frequency")
            }
            Self::NotPositiveDefinite => write!(
                f,
                "permittivity matrix is not positive definite (check rod epsilon and overlap)"
            ),
        }
    }
}

/// Checks the supercell description before any eigenproblem is built.
///
/// Requires an odd `M` and odd `n_g`. Rods must have finite, non-negative
/// radii and positive finite permittivities, lie inside the cell, and not
/// overlap. The rod set must be inversion-symmetric, since the Fourier
/// coefficients are only real for such a set.
fn validate(m: usize, rods: &[Rod], n_g: usize) -> Result<(), SupercellError> {
    let bad = SupercellError::InvalidParameter;
    if m == 0 || m % 2 == 0 {
        return Err(bad("supercell size m must be odd and at least 1"));
    }
    if n_g == 0 || n_g % 2 == 0 {
        return Err(bad("plane-wave count n_g must be odd and at least 1"));
    }
    let h = ((m - 1) / 2) as i32;
    for r in rods {
        if !(r.radius.is_finite() && r.radius >= 0.0) {
            return Err(bad("rod radius must be finite and non-negative"));
        }
        if !(r.epsilon.is_finite() && r.epsilon > 0.0) {
            return Err(bad("rod permittivity must be finite and positive"));
        }
        if r.x.abs() > h || r.y.abs() > h {
            return Err(bad("rod lies outside the supercell"));
        }
    }
    for (i, a) in rods.iter().enumerate() {
        let mirrored = rods
            .iter()
            .any(|b| b.x == -a.x && b.y == -a.y && b.radius == a.radius && b.epsilon == a.epsilon);
        if !mirrored {
            return Err(bad("rod set must be inversion-symmetric"));
        }
        for b in &rods[i + 1..] {
            if a.x == b.x && a.y == b.y {
                return Err(bad("two rods share a cell"));
            }
            let dx = (a.x - b.x) as f64;
            let dy = (a.y - b.y) as f64;
            if a.radius + b.radius > sqrt(dx * dx + dy * dy) {
                return Err(bad("rods overlap"));
            }
        }
    }
    Ok(())
}

impl core::error::Error for SupercellError {}

/// Widest frequency window inside `window = (f_lo, f_hi)` that contains no
/// bulk state of the uniform crystal, found by sampling the supercell
/// spectrum over the full supercell Brillouin zone on an `(n_k + 1)²` grid.
///
/// This is a sampling estimate. A coarse grid can miss band extrema and
/// report a gap that is too wide, and the result also depends on the
/// plane-wave cutoff. Check it with a denser `n_k` and a larger `n_g` before
/// relying on it. `tm` selects E-polarisation (`true`) or H-polarisation.
///
/// Returns `Ok(None)` when no gap exists inside the window, and
/// `Ok(Some((lower, upper)))` in units of `c/a` otherwise. `n_k` must be at
/// least 1, and the window must be finite with `f_lo < f_hi`.
pub fn bulk_gap_from_supercell(
    m: usize,
    epsilon: f64,
    r_over_a: f64,
    n_g: usize,
    n_k: usize,
    window: (f64, f64),
    tm: bool,
) -> Result<Option<(f64, f64)>, SupercellError> {
    let (f_lo, f_hi) = window;
    if n_k == 0 {
        return Err(SupercellError::InvalidParameter("n_k must be >= 1"));
    }
    if !(f_lo.is_finite() && f_hi.is_finite() && f_lo < f_hi) {
        return Err(SupercellError::InvalidParameter(
            "window must be finite with f_lo < f_hi",
        ));
    }
    let rods = point_defect_rods(m, epsilon, r_over_a, epsilon, r_over_a);
    let kmax = PI / m as f64;
    let mut all: Vec<f64> = Vec::new();
    for i in 0..=n_k {
        for j in 0..=n_k {
            let k = (kmax * i as f64 / n_k as f64, kmax * j as f64 / n_k as f64);
            let freqs = if tm {
                supercell_tm_freqs(m, &rods, n_g, k)?
            } else {
                supercell_te_freqs(m, &rods, n_g, k)?
            };
            for f in freqs {
                if !f.is_finite() {
                    return Err(SupercellError::NonFiniteFrequency);
                }
                if f >= f_lo && f <= f_hi {
                    all.push(f);
                }
            }
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
    Ok(best.filter(|(lo, hi)| hi > lo))
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
            let f = solve(m, &rods, 7, (0.0, 0.0)).unwrap();
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
        let fu = supercell_tm_freqs(5, &uniform, 5, (0.0, 0.0)).unwrap();
        let fm = supercell_tm_freqs(5, &missing, 5, (0.0, 0.0)).unwrap();
        let diff: f64 = fu.iter().zip(&fm).map(|(a, b)| (a - b).abs()).sum();
        assert!(diff > 1e-3, "defect had no effect: {diff}");
    }

    #[test]
    fn empty_single_cell_tm_gives_unit_lightline() {
        // M = 1, ε = 1: the first nonzero level at Γ is |G|/2π = 1 (c/a).
        let rods = point_defect_rods(1, 1.0, 0.0, 1.0, 0.0);
        let f = supercell_tm_freqs(1, &rods, 5, (0.0, 0.0)).unwrap();
        assert!(f[0].abs() < 1e-9 && (f[1] - 1.0).abs() < 1e-9, "{f:?}");
    }

    #[test]
    fn missing_rod_level_sits_in_tm_gap() {
        // Square rods, ε = 8.9, r = 0.2a (E_z). The full-zone bulk gap is
        // about (0.324, 0.447) c/a, and a missing-rod defect puts a level at
        // ≈ 0.393 c/a at Γ (converged to ~0.3% for ng ≥ 21). The uniform
        // crystal has no level in (0.36, 0.42) at Γ.
        let m = 5;
        let ng = 17;
        let uniform = point_defect_rods(m, 8.9, 0.2, 8.9, 0.2);
        let missing = point_defect_rods(m, 8.9, 0.2, 1.0, 0.0);
        let window = |f: &[f64]| -> Vec<f64> {
            f.iter()
                .copied()
                .filter(|v| (0.36..0.42).contains(v))
                .collect()
        };
        let bulk = window(&supercell_tm_freqs(m, &uniform, ng, (0.0, 0.0)).unwrap());
        let defect = window(&supercell_tm_freqs(m, &missing, ng, (0.0, 0.0)).unwrap());
        assert!(
            bulk.is_empty(),
            "uniform crystal has levels in gap: {bulk:?}"
        );
        assert_eq!(defect.len(), 1, "expected one defect level: {defect:?}");
        assert!((defect[0] - 0.393).abs() < 0.01, "defect at {}", defect[0]);
    }

    fn is_invalid(r: &Result<Vec<f64>, SupercellError>) -> bool {
        matches!(r, Err(SupercellError::InvalidParameter(_)))
    }

    #[test]
    fn even_supercell_size_is_rejected() {
        // An even M has no centre cell, so the defect cannot sit at the middle.
        let rods = point_defect_rods(3, 8.9, 0.2, 8.9, 0.2);
        assert!(is_invalid(&supercell_tm_freqs(4, &rods, 5, (0.0, 0.0))));
        assert!(is_invalid(&supercell_te_freqs(0, &rods, 5, (0.0, 0.0))));
    }

    #[test]
    fn even_plane_wave_count_is_rejected() {
        let rods = point_defect_rods(3, 8.9, 0.2, 8.9, 0.2);
        assert!(is_invalid(&supercell_tm_freqs(3, &rods, 6, (0.0, 0.0))));
    }

    #[test]
    fn overlapping_rods_are_rejected() {
        // Two rods of radius 0.6a whose centres are a apart overlap.
        let rods = [
            Rod {
                x: 0,
                y: 0,
                radius: 0.6,
                epsilon: 8.9,
            },
            Rod {
                x: 1,
                y: 0,
                radius: 0.6,
                epsilon: 8.9,
            },
            Rod {
                x: -1,
                y: 0,
                radius: 0.6,
                epsilon: 8.9,
            },
        ];
        assert!(is_invalid(&supercell_tm_freqs(3, &rods, 5, (0.0, 0.0))));
    }

    #[test]
    fn asymmetric_rod_set_is_rejected() {
        // The Fourier formulation assumes inversion symmetry. Without it the
        // coefficients are complex and the code would silently drop them.
        let rods = [
            Rod {
                x: 1,
                y: 0,
                radius: 0.2,
                epsilon: 8.9,
            },
            Rod {
                x: 0,
                y: 0,
                radius: 0.0,
                epsilon: 1.0,
            },
        ];
        assert!(is_invalid(&supercell_tm_freqs(3, &rods, 5, (0.0, 0.0))));
    }

    #[test]
    fn non_physical_rod_values_are_rejected() {
        for (radius, eps) in [
            (0.2, f64::NAN),
            (0.2, 0.0),
            (0.2, -4.0),
            (f64::NAN, 8.9),
            (-0.1, 8.9),
        ] {
            let rods = point_defect_rods(3, eps, radius, eps, radius);
            assert!(
                is_invalid(&supercell_tm_freqs(3, &rods, 5, (0.0, 0.0))),
                "{radius} {eps}"
            );
        }
    }

    #[test]
    fn rod_outside_the_supercell_is_rejected() {
        let rods = [
            Rod {
                x: 2,
                y: 0,
                radius: 0.2,
                epsilon: 8.9,
            },
            Rod {
                x: -2,
                y: 0,
                radius: 0.2,
                epsilon: 8.9,
            },
        ];
        assert!(is_invalid(&supercell_tm_freqs(3, &rods, 5, (0.0, 0.0))));
    }

    #[test]
    fn bulk_gap_rejects_zero_k_grid() {
        // With n_k = 0 the sampling grid has only Γ, which used to be reported
        // as a full-window gap. It must be an error instead.
        let r = bulk_gap_from_supercell(5, 8.9, 0.2, 9, 0, (0.3, 0.5), true);
        assert_eq!(r, Err(SupercellError::InvalidParameter("n_k must be >= 1")));
    }

    #[test]
    fn bulk_gap_rejects_non_finite_window() {
        let r = bulk_gap_from_supercell(5, 8.9, 0.2, 9, 2, (f64::NAN, 0.5), true);
        assert!(
            matches!(r, Err(SupercellError::InvalidParameter(_))),
            "{r:?}"
        );
        let r = bulk_gap_from_supercell(5, 8.9, 0.2, 9, 2, (0.3, f64::INFINITY), true);
        assert!(
            matches!(r, Err(SupercellError::InvalidParameter(_))),
            "{r:?}"
        );
    }

    #[test]
    fn bulk_gap_rejects_inverted_window() {
        let r = bulk_gap_from_supercell(5, 8.9, 0.2, 9, 2, (0.5, 0.3), true);
        assert!(
            matches!(r, Err(SupercellError::InvalidParameter(_))),
            "{r:?}"
        );
    }

    #[test]
    fn bulk_gap_is_found_for_square_rods() {
        // The full-zone TM gap is (0.324, 0.447) c/a. A coarse grid must
        // report a gap that lies inside it, not one that is too wide.
        let (lo, hi) = bulk_gap_from_supercell(5, 8.9, 0.2, 13, 4, (0.2, 0.6), true)
            .unwrap()
            .expect("gap");
        assert!(lo >= 0.30 && hi <= 0.46 && hi > lo, "({lo}, {hi})");
    }
}
