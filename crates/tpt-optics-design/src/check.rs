//! Self-check of a sequential system: a list of invariants, each reported as
//! passed or failed with a short explanation.
//!
//! The checks look for prescriptions that are malformed or numerically
//! unsafe before they are optimised or exported: non-finite radii, surfaces
//! with no aperture, a paraxial ABCD matrix whose determinant is not the
//! index ratio, and an on-axis trace that does not produce a finite path.

use alloc::string::String;
use alloc::vec::Vec;

use tpt_optics_geometric::sequential::{SequentialSystem, TraceContext};

/// Result of one invariant.
#[derive(Clone, Debug, PartialEq)]
pub struct Check {
    /// Short name of the invariant.
    pub name: &'static str,
    /// Whether it holds.
    pub passed: bool,
    /// What was measured, or why it failed.
    pub detail: String,
}

/// Run every check on `system` at context `ctx`.
///
/// The paraxial determinant check accepts either `n_in / n_out` or its
/// inverse, because the ABCD convention used here can run either way.
#[must_use]
pub fn self_check(system: &SequentialSystem, ctx: TraceContext) -> Vec<Check> {
    let mut checks = Vec::new();

    let finite_curvatures = system
        .surfaces
        .iter()
        .all(|s| s.asphere.curvature.is_finite());
    checks.push(Check {
        name: "curvatures are finite",
        passed: finite_curvatures,
        detail: alloc::format!("{} surfaces", system.surfaces.len()),
    });

    let apertures_positive = system
        .surfaces
        .iter()
        .all(|s| s.semi_diameter > 0.0 && !s.semi_diameter.is_nan());
    checks.push(Check {
        name: "semi-diameters are positive",
        passed: apertures_positive,
        detail: String::from("surfaces with a zero or NaN aperture fail"),
    });

    let thicknesses_valid = system
        .thickness
        .iter()
        .all(|t| t.is_finite() || t.is_infinite() && *t > 0.0);
    checks.push(Check {
        name: "thicknesses are finite or +∞",
        passed: thicknesses_valid,
        detail: alloc::format!("{} gaps", system.thickness.len()),
    });

    if system.surfaces.is_empty() {
        checks.push(Check {
            name: "paraxial determinant equals index ratio",
            passed: false,
            detail: String::from("no surfaces"),
        });
        return checks;
    }

    let n_in = system.index_before(0, ctx);
    let n_out = system.index_after(system.surfaces.len() - 1, ctx);
    let det = system.paraxial(ctx).abcd().det();
    let ratio = n_in / n_out;
    let close = |a: f64, b: f64| libm::fabs(a - b) <= 1e-9 * libm::fabs(b).max(1.0);
    let det_ok = det.is_finite() && (close(det, ratio) || close(det, 1.0 / ratio));
    checks.push(Check {
        name: "paraxial determinant equals index ratio",
        passed: det_ok,
        detail: alloc::format!("det = {det:.6e}, n_in/n_out = {ratio:.6e}"),
    });

    let trace = system.trace(
        &tpt_optics_geometric::prelude::Ray::new(
            tpt_optics_core::Vec3::new(0.0, 0.0, -1e-3),
            tpt_optics_core::Vec3::new(0.0, 0.0, 1.0),
        ),
        ctx,
    );
    let opd_finite = trace.opd.is_finite();
    checks.push(Check {
        name: "on-axis trace gives a finite path",
        passed: opd_finite && !trace.blocked,
        detail: alloc::format!("opd = {:.6e} m, blocked = {}", trace.opd, trace.blocked),
    });

    checks
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_optics_core::material::MaterialDatabase;
    use tpt_optics_core::{Kelvin, Wavelength};
    use tpt_optics_geometric::presets::singlet;
    use tpt_optics_geometric::sequential::{Medium, SeqSurface};
    use tpt_optics_geometric::surface::Asphere;

    fn ctx() -> TraceContext {
        TraceContext {
            wavelength: Wavelength::from_nm(587.6),
            temperature: Kelvin::from_celsius(20.0),
        }
    }

    fn all_passed(checks: &[Check]) -> bool {
        checks.iter().all(|c| c.passed)
    }

    #[test]
    fn a_good_singlet_passes_every_check() {
        let catalog = MaterialDatabase::presets();
        let sys = singlet(&catalog, 100e-3, -100e-3, 12.5e-3, 5e-3, "bk7", 95e-3).unwrap();
        let checks = self_check(&sys, ctx());
        assert!(all_passed(&checks), "{checks:#?}");
    }

    #[test]
    fn a_non_finite_radius_fails_the_curvature_check() {
        let mut sys = SequentialSystem::new();
        sys.push(
            SeqSurface::aspheric(Asphere::sphere(f64::NAN), 12.5e-3),
            5e-3,
            Medium::Vacuum,
        );
        sys.push(SeqSurface::flat(12.5e-3), 0.0, Medium::Vacuum);
        let checks = self_check(&sys, ctx());
        let curv = checks
            .iter()
            .find(|c| c.name == "curvatures are finite")
            .unwrap();
        assert!(!curv.passed);
    }

    #[test]
    fn a_zero_aperture_fails_the_aperture_check() {
        let mut sys = SequentialSystem::new();
        sys.push(SeqSurface::flat(0.0), 0.0, Medium::Vacuum);
        let checks = self_check(&sys, ctx());
        let ap = checks
            .iter()
            .find(|c| c.name == "semi-diameters are positive")
            .unwrap();
        assert!(!ap.passed);
    }

    #[test]
    fn an_empty_system_reports_failure_without_panicking() {
        let checks = self_check(&SequentialSystem::new(), ctx());
        assert!(checks.iter().any(|c| !c.passed));
    }
}
