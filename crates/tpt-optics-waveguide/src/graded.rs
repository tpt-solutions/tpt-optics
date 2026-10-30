//! Weak-guidance LP modes of a graded-index circular fiber with a power-law
//! core profile.
//!
//! The core index follows `n²(r) = n₁²·(1 − 2Δ·(r/a)^α)` for `r < a` and is
//! constant at `n₂` in the cladding, with `Δ = (n₁² − n₂²)/(2n₁²)`. Writing
//! `E = k₀²n₁² − β²`, the scalar mode equation in the radial coordinate is
//!
//! `−(1/r)(r u′)′ + (l²/r²) u + V(r) u = E u`, with `V(r) = k₀²(n₁² − n²(r))`.
//!
//! The problem is discretised on a cell-centred radial grid with Dirichlet
//! walls at `r = R` (the truncation radius) and solved as a dense symmetric
//! eigenproblem. Guided modes are the eigenvalues below the cladding barrier
//! `E < k₀²(n₁² − n₂²)`.
//!
//! Accuracy limits: scalar (weak-guidance) approximation, second-order radial
//! discretisation, and a truncation at `R = 4a`. The cladding decay of a mode
//! is `exp(−w·(r − a))` with `w·a = V√(1−b)`; at `R = 4a` the wall-induced
//! error is below 1 % for modes with `V ≳ 2`, and it grows for modes near
//! cutoff, which leak further than `R` allows.
//!
//! Limit checks: for `α → ∞` the profile becomes a step and the result
//! approaches [`crate::fiber::fiber_lp_modes`]; for `α = 2` the ground state
//! matches the exact 2-D harmonic-oscillator eigenvalue `E = 2ω(2m+|l|+1)`
//! with `ω = k₀n₁√(2Δ)/a`.

use alloc::vec;
use alloc::vec::Vec;

use libm::sqrt;
use tpt_optics_core::Wavelength;

use crate::eigen::symmetric_eigenvalues;
use crate::fiber::LpMode;

/// Radial grid size used by [`graded_lp_modes`].
pub const DEFAULT_RADIAL_POINTS: usize = 200;

/// Guided LP modes of a power-law graded-index fiber with `l ≤ l_max`.
///
/// `alpha` is the profile exponent (`2` = parabolic, large values approach a
/// step index). `n_radial` is the number of radial grid cells; use
/// [`DEFAULT_RADIAL_POINTS`] unless convergence is being checked.
pub fn graded_lp_modes(
    core_radius: f64,
    n_core: f64,
    n_clad: f64,
    alpha: f64,
    wl: Wavelength,
    l_max: u32,
    n_radial: usize,
) -> Vec<LpMode> {
    let k0 = wl.wavenumber();
    let (n1s, n2s) = (n_core * n_core, n_clad * n_clad);
    let barrier = k0 * k0 * (n1s - n2s);
    if barrier <= 0.0 || n_radial < 4 || alpha <= 0.0 {
        return Vec::new();
    }
    let a = core_radius;
    let big_r = 4.0 * a;
    let h = big_r / n_radial as f64;
    // Cell centres and faces: r_i = (i + ½)h, faces at i·h.
    let r_cell: Vec<f64> = (0..n_radial).map(|i| (i as f64 + 0.5) * h).collect();
    let r_face: Vec<f64> = (0..=n_radial).map(|i| i as f64 * h).collect();
    // Potential V(r) = k₀²(n₁² − n²(r)) = k₀²(n₁² − n₂²)·min((r/a)^α, 1).
    let v_at = |r: f64| barrier * powf_clamped(r / a, alpha);

    let mut modes: Vec<LpMode> = Vec::new();
    for l in 0..=l_max {
        let ll = (l * l) as f64;
        // Symmetric form of the r-weighted operator: K u = E·diag(r) u.
        let mut k = vec![0.0; n_radial * n_radial];
        for i in 0..n_radial {
            let r_in = r_face[i];
            let r_out = r_face[i + 1];
            // At i = 0 the inner face r = 0 contributes nothing.
            k[i * n_radial + i] =
                (r_in + r_out) / (h * h) + ll / r_cell[i] + r_cell[i] * v_at(r_cell[i]);
            if i + 1 < n_radial {
                let off = -r_out / (h * h);
                k[i * n_radial + i + 1] = off;
                k[(i + 1) * n_radial + i] = off;
            }
        }
        // Diagonal rescaling by √r_i gives a standard symmetric problem.
        let s: Vec<f64> = r_cell.iter().map(|&r| sqrt(r)).collect();
        for i in 0..n_radial {
            for j in 0..n_radial {
                k[i * n_radial + j] /= s[i] * s[j];
            }
        }

        let eig = symmetric_eigenvalues(&mut k, n_radial);

        let mut m = 0_u32;
        for &e in &eig {
            // Guided: β² > k₀²n₂² ⇔ E < barrier. Small margin rejects the
            // box-state continuum just above the barrier.
            if e >= barrier * (1.0 - 1e-9) {
                break;
            }
            if e <= 0.0 {
                continue;
            }
            m += 1;
            modes.push(LpMode {
                l,
                m,
                n_eff: sqrt(n1s - e / (k0 * k0)),
            });
        }
    }
    modes
}

/// `x^p` for `x ≥ 0`, saturating at 1 for `x ≥ 1` (the cladding plateau).
fn powf_clamped(x: f64, p: f64) -> f64 {
    if x >= 1.0 {
        1.0
    } else if x <= 0.0 {
        0.0
    } else {
        libm::pow(x, p)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fiber::fiber_lp_modes;

    fn wl(nm: f64) -> Wavelength {
        Wavelength::from_nm(nm)
    }

    #[test]
    fn parabolic_ground_state_matches_oscillator() {
        // Δ = 1 %, a = 13 µm, λ = 850 nm → ω·a² ≈ 20, so the mode sits well
        // inside the core and the harmonic-oscillator form is accurate.
        let (n1, delta, a) = (1.456, 0.01, 13.0e-6);
        let n2 = n1 * sqrt(1.0 - 2.0 * delta);
        let w = wl(850.0);
        let k0 = w.wavenumber();
        let omega = k0 * n1 * sqrt(2.0 * delta) / a;
        // 2-D oscillator: E = 2ω(2m + |l| + 1); ground state m = 0, l = 0.
        let e_exact = 2.0 * omega;
        let n_exact = sqrt(n1 * n1 - e_exact / (k0 * k0));

        let modes = graded_lp_modes(a, n1, n2, 2.0, w, 0, DEFAULT_RADIAL_POINTS);
        assert_ne!(modes.len(), 0);
        let ground = modes[0];
        assert_eq!((ground.l, ground.m), (0, 1));
        assert!(
            (ground.n_eff - n_exact).abs() < 2e-5,
            "numeric {} vs analytic {}",
            ground.n_eff,
            n_exact
        );
    }

    #[test]
    fn large_alpha_approaches_step_index() {
        // Single-mode step fiber: compare LP01 against the Bessel solver.
        let (a, n1, n2) = (4.1e-6, 1.4515, 1.4447);
        let w = wl(1550.0);
        let step = fiber_lp_modes(a, n1, n2, w, 0);
        let graded = graded_lp_modes(a, n1, n2, 40.0, w, 0, DEFAULT_RADIAL_POINTS);
        assert_eq!(step.len(), 1);
        assert_eq!(graded.len(), 1, "{graded:?}");
        assert!(
            (graded[0].n_eff - step[0].n_eff).abs() < 2e-4,
            "graded {} vs step {}",
            graded[0].n_eff,
            step[0].n_eff
        );
    }

    #[test]
    fn graded_index_supports_more_modes_than_step() {
        // Same core/cladding: a parabolic profile guides modes whose β
        // decreases with order; ordering within each l must be monotone.
        let (a, n1, n2) = (25e-6, 1.456, 1.44);
        let modes = graded_lp_modes(a, n1, n2, 2.0, wl(850.0), 3, DEFAULT_RADIAL_POINTS);
        assert!(modes.len() > 4, "{} modes", modes.len());
        for l in 0..=3 {
            let mut prev = f64::INFINITY;
            for m in modes.iter().filter(|m| m.l == l) {
                assert!(m.n_eff < prev && m.n_eff > n2 && m.n_eff < n1);
                prev = m.n_eff;
            }
        }
    }

    #[test]
    fn zero_contrast_has_no_modes() {
        assert_eq!(
            graded_lp_modes(5e-6, 1.444, 1.444, 2.0, wl(1550.0), 1, 64).len(),
            0
        );
    }
}
