//! Waveguide presets: a single-mode fibre, a silicon ring resonator, and a
//! photonic-crystal line-defect waveguide.
//!
//! The parameters are typical published values, not vendor data. Check them
//! against your datasheet or process before using a result for design.

use alloc::vec::Vec;

use tpt_optics_core::Wavelength;

use crate::coupled::{ring_fsr, ring_resonances, CoupledError};
use crate::fiber::{fiber_lp_modes, LpMode};
use crate::supercell::{point_defect_rods, Rod};

/// Step-index fibre parameters.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct StepFiber {
    /// Core radius (m).
    pub core_radius: f64,
    /// Core index.
    pub n_core: f64,
    /// Cladding index.
    pub n_clad: f64,
}

/// Typical standard single-mode fibre (SMF-28 class): 8.2 µm core diameter,
/// Δn ≈ 0.0057 at 1550 nm. This is illustrative, not a vendor specification.
pub const SMF28: StepFiber = StepFiber {
    core_radius: 4.1e-6,
    n_core: 1.4504,
    n_clad: 1.4447,
};

/// Guided LP modes of a fibre at a wavelength.
#[must_use]
pub fn fiber_modes(fiber: StepFiber, wl: Wavelength, l_max: u32) -> Vec<LpMode> {
    fiber_lp_modes(fiber.core_radius, fiber.n_core, fiber.n_clad, wl, l_max)
}

/// Silicon-wire ring resonator (500 nm × 220 nm wire class) in silica.
///
/// `resonances` uses the phase index, so it ignores dispersion. `fsr` uses
/// the group index. The two differ by about `n_g/n_eff`.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct SiliconRing {
    /// Ring radius (m).
    pub radius: f64,
    /// Effective index of the wire mode near 1550 nm.
    pub n_eff: f64,
    /// Group index near 1550 nm.
    pub n_group: f64,
}

impl SiliconRing {
    /// A ring of radius `radius` (m) with typical silicon-wire indices.
    #[must_use]
    pub fn new(radius: f64) -> Self {
        SiliconRing {
            radius,
            n_eff: 2.4,
            n_group: 4.2,
        }
    }

    /// Ring circumference (m).
    #[must_use]
    pub fn circumference(&self) -> f64 {
        core::f64::consts::TAU * self.radius
    }

    /// Free spectral range at `wl`, or `None` for a degenerate ring.
    #[must_use]
    pub fn fsr(&self, wl: Wavelength) -> Option<f64> {
        ring_fsr(self.n_group, self.circumference(), wl)
    }

    /// Resonance wavelengths (m) within `[wl_min, wl_max]`.
    ///
    /// # Errors
    /// Returns the [`CoupledError`] from [`ring_resonances`] for a bad window.
    pub fn resonances(&self, wl_min: f64, wl_max: f64) -> Result<Vec<f64>, CoupledError> {
        ring_resonances(self.n_eff, self.circumference(), wl_min, wl_max)
    }
}

/// Rods of a square photonic-crystal supercell with one full row of rods
/// removed through the centre. The result is a line-defect waveguide along
/// x, for use with the supercell solvers.
///
/// `m` is the supercell size (odd). The rods have radius `r_over_a` and
/// permittivity `epsilon`.
#[must_use]
pub fn pc_line_defect_rods(m: usize, epsilon: f64, r_over_a: f64) -> Vec<Rod> {
    point_defect_rods(m, epsilon, r_over_a, epsilon, r_over_a)
        .into_iter()
        .filter(|rod| rod.y != 0)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::supercell::supercell_tm_freqs;

    #[test]
    fn smf28_guides_only_the_fundamental_mode_at_1550() {
        // V ≈ 2.13 at 1550 nm, below the LP11 cutoff of 2.405.
        let modes = fiber_modes(SMF28, Wavelength::from_nm(1550.0), 2);
        assert!(modes.iter().any(|m| m.l == 0 && m.m == 1));
        assert!(modes.iter().all(|m| !(m.l == 1 && m.m == 1)), "{modes:?}");
    }

    #[test]
    fn smf28_guides_lp11_at_1200_nm() {
        // V ≈ 2.7 at 1200 nm, above the LP11 cutoff.
        let modes = fiber_modes(SMF28, Wavelength::from_nm(1200.0), 2);
        assert!(modes.iter().any(|m| m.l == 1 && m.m == 1), "{modes:?}");
    }

    #[test]
    fn ring_fsr_is_about_two_nanometres_for_fifty_micron_radius() {
        let ring = SiliconRing::new(50e-6);
        let fsr = ring.fsr(Wavelength::from_nm(1550.0)).unwrap();
        assert!(fsr > 1.0e-9 && fsr < 3.0e-9, "{fsr}");
    }

    #[test]
    fn ring_resonances_are_spaced_by_lambda_squared_over_n_eff_l() {
        // Positions come from m·λ = n_eff·L with no dispersion, so adjacent
        // resonances are λ²/(n_eff·L) apart. The group-index FSR is the
        // dispersive value and is not what this model predicts.
        let ring = SiliconRing::new(50e-6);
        let res = ring.resonances(1540e-9, 1560e-9).unwrap();
        assert!(res.len() >= 2);
        let spacing = res[1] - res[0];
        let expected = (1550e-9_f64).powi(2) / (ring.n_eff * ring.circumference());
        assert!(
            (spacing - expected).abs() < 0.02 * expected,
            "{spacing} vs {expected}"
        );
    }

    #[test]
    fn line_defect_removes_one_row_of_rods() {
        let rods = pc_line_defect_rods(5, 8.9, 0.2);
        assert_eq!(rods.len(), 5 * 5 - 5);
        assert!(rods.iter().all(|r| r.y != 0));
    }

    #[test]
    fn line_defect_supercell_solves_and_has_a_level_in_the_gap() {
        // Uniform crystal: no level in (0.36, 0.42) c/a at Γ. With the row
        // removed, a guided band appears in the gap.
        let (m, eps, r) = (5, 8.9, 0.2);
        let rods = pc_line_defect_rods(m, eps, r);
        let freqs = supercell_tm_freqs(m, &rods, 13, (0.0, 0.0)).expect("valid supercell");
        let in_gap: Vec<f64> = freqs
            .into_iter()
            .filter(|f| (0.36..0.42).contains(f))
            .collect();
        assert!(!in_gap.is_empty(), "no guided level in the gap");
    }
}
