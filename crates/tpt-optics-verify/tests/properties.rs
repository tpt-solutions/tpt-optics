//! Property-based verification of the optical laws (Phase 3 exit criteria),
//! exercising the shared strategies and invariants.

use proptest::prelude::*;
use tpt_optics_core::{Kelvin, Matrix2, Vec3};
use tpt_optics_geometric::abcd;
use tpt_optics_geometric::sequential::{Medium, SeqSurface, SequentialSystem};
use tpt_optics_geometric::surface::Asphere;
use tpt_optics_verify::invariants as inv;
use tpt_optics_verify::strategies as strat;

proptest! {
    #![proptest_config(ProptestConfig::with_cases(1_000))]

    #[test]
    fn snell_reversibility(
        n1 in strat::refractive_index(),
        n2 in strat::refractive_index(),
        d in strat::incident_direction(),
    ) {
        let n = Vec3::new(0.0, 0.0, 1.0);
        inv::snell_reversible(d, n, n1, n2, 1e-10).map_err(proptest::test_runner::TestCaseError::fail)?;
    }

    #[test]
    fn snell_momentum(
        n1 in strat::refractive_index(),
        n2 in strat::refractive_index(),
        d in strat::incident_direction(),
    ) {
        let n = Vec3::new(0.0, 0.0, 1.0);
        inv::snell_tangential_momentum(d, n, n1, n2, 1e-9).map_err(proptest::test_runner::TestCaseError::fail)?;
    }

    #[test]
    fn fresnel_energy(
        n1 in strat::refractive_index(),
        n2 in strat::refractive_index(),
        sin_t in 0.0_f64..0.95,
    ) {
        inv::fresnel_energy_conserved(n1, n2, sin_t.asin(), 1e-9).map_err(proptest::test_runner::TestCaseError::fail)?;
    }

    #[test]
    fn arbitrary_unit_index_compositions_symplectic(
        f1 in strat::curvature(),
        f2 in strat::curvature(),
        d in 0.0_f64..1.0,
    ) {
        let focal1 = 1.0 / f1.abs().max(0.05) * f1.signum();
        let focal2 = 1.0 / f2.abs().max(0.05) * f2.signum();
        let m = thin(focal1).mul_matrix(abcd::free_space(d)).mul_matrix(thin(focal2));
        inv::abcd_symplectic(m, 1e-9).map_err(proptest::test_runner::TestCaseError::fail)?;
    }

    #[test]
    fn refraction_det_is_index_ratio(
        n1 in strat::refractive_index(),
        n2 in strat::refractive_index(),
        c in strat::curvature(),
    ) {
        inv::abcd_refraction_det(n1, n2, c, 1e-10).map_err(proptest::test_runner::TestCaseError::fail)?;
    }

    #[test]
    fn singlet_trace_is_reversible(
        thickness in 20.0_f64..200.0,
    ) {
        let mut sys = singlet();
        sys.thickness[1] = thickness * 1e-3;
        let ctx = tpt_optics_geometric::sequential::TraceContext {
            wavelength: inv::d_line(),
            temperature: Kelvin::from_celsius(20.0),
        };
        let ray = tpt_optics_geometric::prelude::Ray::new(
            Vec3::new(0.5e-3, 0.0, -1e-6),
            Vec3::new(0.0, 0.0, 1.0),
        );
        inv::trace_reversibility(&sys, &ray, ctx, 1e-6).map_err(proptest::test_runner::TestCaseError::fail)?;
    }

    #[test]
    fn singlet_opl_rotation_invariant(
        azimuth in 0.0_f64..core::f64::consts::TAU,
        height in 0.0_f64..5e-3,
    ) {
        let sys = singlet();
        let ctx = tpt_optics_geometric::sequential::TraceContext {
            wavelength: inv::d_line(),
            temperature: Kelvin::from_celsius(20.0),
        };
        let ray = tpt_optics_geometric::prelude::Ray::new(
            Vec3::new(height, 0.0, -0.01),
            Vec3::new(0.0, 0.0, 1.0),
        );
        inv::opl_rotation_invariance(&sys, &ray, ctx, azimuth, 1e-9).map_err(proptest::test_runner::TestCaseError::fail)?;
    }
}

fn thin(f: f64) -> Matrix2<f64> {
    abcd::thin_lens(f)
}

/// Biconvex BK7-like singlet used by the trace-level properties.
fn singlet() -> SequentialSystem {
    let mut sys = SequentialSystem::new();
    sys.push(
        SeqSurface::aspheric(Asphere::sphere(10.0), 12e-3),
        5e-3,
        Medium::Const(1.5168),
    );
    sys.push(
        SeqSurface::aspheric(Asphere::sphere(-10.0), 12e-3),
        0.0,
        Medium::Const(1.0),
    );
    sys.push(SeqSurface::flat(12e-3), 0.0, Medium::Vacuum);
    sys
}
