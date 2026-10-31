//! System-level invariant checks: the verify crate's trace invariants run on a
//! fixed set of rays through a sequential system.
//!
//! Each check returns a name, pass or fail, and a short detail. The checks are
//! the same invariants the property tests use: a trace retraces the mirrored
//! system, and the optical path is unchanged by rotation about the axis.

use alloc::string::String;
use alloc::vec::Vec;

use tpt_optics_core::Vec3;
use tpt_optics_geometric::prelude::Ray;
use tpt_optics_geometric::sequential::{SequentialSystem, TraceContext};

use crate::invariants::{opl_rotation_invariance, trace_reversibility};

/// Result of one system invariant.
#[derive(Clone, Debug, PartialEq)]
pub struct InvariantCheck {
    /// Name of the invariant and the ray it was run on.
    pub name: &'static str,
    /// Whether it held.
    pub passed: bool,
    /// The failure message, or an empty string when it passed.
    pub detail: String,
}

/// Rays used by [`verify_system`]: on axis, off axis, and tilted, all starting
/// 1 mm before the first surface.
fn rays() -> [(&'static str, Ray<f64>); 3] {
    let tilt = 2f64.to_radians();
    [
        (
            "on-axis ray",
            Ray::new(Vec3::new(0.0, 0.0, -1e-3), Vec3::new(0.0, 0.0, 1.0)),
        ),
        (
            "2 mm off-axis ray",
            Ray::new(Vec3::new(2e-3, 0.0, -1e-3), Vec3::new(0.0, 0.0, 1.0)),
        ),
        (
            "2° tilted ray",
            Ray::new(
                Vec3::new(1e-3, 0.0, -1e-3),
                Vec3::new(libm::sin(tilt), 0.0, libm::cos(tilt)),
            ),
        ),
    ]
}

/// Run the trace invariants through `system`: reversibility on the axis
/// (tolerance 1e-6 m), and optical-path invariance under rotation for each
/// ray (tolerance 1e-9 m). A blocked ray gives a pass, because no path is
/// defined past the block.
#[must_use]
pub fn verify_system(system: &SequentialSystem, ctx: TraceContext) -> Vec<InvariantCheck> {
    let mut out = Vec::new();
    for (label, ray) in rays() {
        let reversal_name = match label {
            "on-axis ray" => "trace reversibility (on-axis ray)",
            "2 mm off-axis ray" => "trace reversibility (2 mm off-axis ray)",
            _ => "trace reversibility (2° tilted ray)",
        };
        out.push(to_check(
            reversal_name,
            trace_reversibility(system, &ray, ctx, 1e-6),
        ));
        let rotation_name = match label {
            "on-axis ray" => "optical path invariant under rotation (on-axis ray)",
            "2 mm off-axis ray" => "optical path invariant under rotation (2 mm off-axis ray)",
            _ => "optical path invariant under rotation (2° tilted ray)",
        };
        out.push(to_check(
            rotation_name,
            opl_rotation_invariance(system, &ray, ctx, 0.7, 1e-9),
        ));
    }
    out
}

fn to_check(name: &'static str, result: Result<(), String>) -> InvariantCheck {
    match result {
        Ok(()) => InvariantCheck {
            name,
            passed: true,
            detail: String::new(),
        },
        Err(detail) => InvariantCheck {
            name,
            passed: false,
            detail,
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_optics_core::material::MaterialDatabase;
    use tpt_optics_core::{Kelvin, Wavelength};
    use tpt_optics_geometric::sequential::{Medium, SeqSurface};
    use tpt_optics_geometric::surface::Asphere;

    fn ctx() -> TraceContext {
        TraceContext {
            wavelength: Wavelength::from_nm(587.6),
            temperature: Kelvin::from_celsius(20.0),
        }
    }

    fn singlet() -> SequentialSystem {
        let bk7 = MaterialDatabase::presets().get("bk7").unwrap();
        let mut sys = SequentialSystem::new();
        sys.push(
            SeqSurface::aspheric(Asphere::sphere(10.0), 12.5e-3),
            5e-3,
            Medium::Material(bk7),
        );
        sys.push(
            SeqSurface::aspheric(Asphere::sphere(-10.0), 12.5e-3),
            95e-3,
            Medium::Vacuum,
        );
        sys
    }

    #[test]
    fn a_singlet_passes_every_system_invariant() {
        let checks = verify_system(&singlet(), ctx());
        assert_eq!(checks.len(), 6);
        assert!(checks.iter().all(|c| c.passed), "{checks:#?}");
    }
}
