//! Weak-guidance LP modes of a step-index circular fiber.
//!
//! Eigenvalue equations (Gloge):
//! - `l = 0`: `u·J₁(u)/J₀(u) = w·K₁(w)/K₀(w)`
//! - `l ≥ 1`: `u·J_{l−1}(u)/J_l(u) = −w·K_{l−1}(w)/K_l(w)`
//!
//! with `u = a·√(k₀²n₁² − β²)`, `w = a·√(β² − k₀²n₂²)`, `u² + w² = V²`,
//! `V = k₀·a·√(n₁² − n₂²)`. LPₗₘ modes have cutoff at the m-th zero of the
//! appropriate Bessel function (LP₀₁ is never cut off; LP₁₁ cuts off at
//! V = 2.4048).

use alloc::vec::Vec;

use libm::sqrt;
use tpt_optics_core::{bessel_j, bessel_k, Wavelength};

/// One solved LP mode.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct LpMode {
    /// Azimuthal order l.
    pub l: u32,
    /// Radial order m (1-based).
    pub m: u32,
    /// Effective index β/k₀.
    pub n_eff: f64,
}

/// All guided LP modes with `l ≤ l_max` at wavelength `wl`.
pub fn fiber_lp_modes(
    core_radius: f64,
    n_core: f64,
    n_clad: f64,
    wl: Wavelength,
    l_max: u32,
) -> Vec<LpMode> {
    let k0 = wl.wavenumber();
    let v = k0 * core_radius * sqrt(n_core * n_core - n_clad * n_clad);
    if v <= 0.0 {
        return Vec::new();
    }
    let mut modes: Vec<LpMode> = Vec::new();
    for l in 0..=l_max {
        // Radial order m runs while the bracket contains a root.
        // Pole-free forms: multiply through by J_l(u)·K_l(w) so Bessel
        // zeros of the denominator don't create spurious sign flips.
        let f = |u: f64| -> f64 {
            let w = sqrt((v * v - u * u).max(1e-12));
            let li = l as i32;
            if l == 0 {
                u * bessel_j(1, u) * bessel_k(0, w) - w * bessel_k(1, w) * bessel_j(0, u)
            } else {
                u * bessel_j(li - 1, u) * bessel_k(li, w)
                    + w * bessel_k(li - 1, w) * bessel_j(li, u)
            }
        };
        // Scan u ∈ (max(ε, previous root), min(V, ...)) for sign changes.
        let n_scan = 400;
        let mut prev_u = 1e-6_f64;
        let mut prev_f = f(prev_u);
        for i in 1..=n_scan {
            let u = (v - 1e-6) * i as f64 / n_scan as f64;
            let fu = f(u);
            if prev_f.is_finite() && fu.is_finite() && prev_f * fu <= 0.0 {
                // Bisect.
                let (mut lo, mut hi) = (prev_u, u);
                for _ in 0..70 {
                    let mid = 0.5 * (lo + hi);
                    if f(lo) * f(mid) <= 0.0 {
                        hi = mid;
                    } else {
                        lo = mid;
                    }
                }
                let u_root = 0.5 * (lo + hi);
                // Reject roots pinned at the wall u ≈ V (w → 0, K diverges).
                if v - u_root < 0.02 * v {
                    continue;
                }
                let kn = k0 * n_core;
                let ku = u_root / core_radius;
                let beta2 = kn * kn - ku * ku;
                if beta2 > (k0 * n_clad) * (k0 * n_clad) {
                    modes.push(LpMode {
                        l,
                        m: (modes.iter().filter(|m| m.l == l).count() + 1) as u32,
                        n_eff: sqrt(beta2) / k0,
                    });
                }
            }
            prev_u = u;
            prev_f = fu;
        }
    }
    modes
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
    fn single_mode_fiber_at_1550() {
        // SMF-28-like: a = 4.1 µm, Δ ≈ 0.36% → V ≈ 2.33 < 2.405 → LP01 only.
        let a = 4.1e-6;
        let (n1, n2) = (1.4515, 1.4447);
        let modes = fiber_lp_modes(a, n1, n2, wl(1550.0), 2);
        assert_eq!(modes.len(), 1, "{modes:?}");
        assert_eq!(modes[0].l, 0);
        let v = wl(1550.0).wavenumber() * a * sqrt(n1 * n1 - n2 * n2);
        assert!(close(v, 2.33, 0.08), "V = {v}");
        assert!(modes[0].n_eff > n2 && modes[0].n_eff < n1);
        // Weak guidance: n_eff ≈ n2 + (n1−n2)·b with b ∈ (0,1).
        assert!(modes[0].n_eff - n2 < n1 - n2);
    }

    #[test]
    fn lp11_cutoff_near_2p405() {
        // Scan V by varying wavelength; LP11 (l=1,m=1) appears just above
        // V = 2.4048.
        let a = 5.0e-6;
        let (n1, n2) = (1.4436, 1.44);
        let v_at = |nm: f64| wl(nm).wavenumber() * a * sqrt(n1 * n1 - n2 * n2);
        let count_at = |nm: f64| fiber_lp_modes(a, n1, n2, wl(nm), 1).len();
        let v_below = v_at(1700.0); // larger λ → smaller V
        let v_above = v_at(1100.0);
        assert!(v_below < 2.4048 && v_above > 2.4048, "{v_below} {v_above}");
        assert_eq!(count_at(1700.0), 1);
        assert_eq!(count_at(1100.0), 2, "LP11 should exist above cutoff");
    }

    #[test]
    fn v_number_zero_no_modes() {
        assert!(fiber_lp_modes(1e-6, 1.444, 1.444, wl(1550.0), 1).is_empty());
    }

    #[test]
    fn many_modes_in_multimode_fiber() {
        // 50 µm core, NA ≈ 0.2 at 850 nm: V ≈ 37 → many mode groups.
        let modes = fiber_lp_modes(25e-6, 1.456, 1.44, wl(850.0), 5);
        assert!(modes.len() > 8, "{} modes", modes.len());
        // All n_eff bounded and decreasing in β with order.
        for m in &modes {
            assert!(m.n_eff > 1.44 && m.n_eff < 1.456);
        }
    }
}
