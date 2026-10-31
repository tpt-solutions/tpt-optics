//! Property-based verification of the optical laws (Phase 3 exit criteria),
//! exercising the shared strategies and invariants.

use proptest::prelude::*;
use tpt_optics_core::fresnel;
use tpt_optics_core::{Kelvin, Matrix2, Vec3};
use tpt_optics_geometric::abcd;
use tpt_optics_geometric::prelude::{Ray, Sphere};
use tpt_optics_geometric::sequential::{Medium, SeqSurface, SequentialSystem};
use tpt_optics_geometric::surface::{Asphere, OpticalSurface};
use tpt_optics_verify::bounded::{self, Axis};
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

// The intersection and energy checks run on fixed grids, so every point is
// checked on every run and a failure names its exact input.
#[test]
fn sphere_hit_is_ahead_and_on_surface() {
    bounded::assert_holds(
        [
            Axis::new(1.0, 10.0, 10),
            Axis::new(-3.5, 3.5, 8),
            Axis::new(-3.5, 3.5, 8),
            Axis::new(-100.0, -1e-3, 6),
            Axis::new(0.5, 1.0, 4),
        ],
        |[curvature, ox, oy, oz, dz]| {
            let ray = Ray::new(Vec3::new(ox, oy, oz), Vec3::new(0.0, 0.0, dz));
            let Some(hit) = Sphere::new(curvature).intersect(&ray, 1e-9) else {
                return Ok(());
            };
            if hit.t <= 0.0 {
                return Err(format!("hit parameter {} not positive", hit.t));
            }
            // |p − C|² = R² with the centre C = (0, 0, R).
            let r = 1.0 / curvature;
            let (px, py, pz) = (hit.point.x, hit.point.y, hit.point.z - r);
            let residual = px * px + py * py + pz * pz - r * r;
            if residual.abs() >= 1e-6 * r * r {
                return Err(format!("off the sphere by {residual}"));
            }
            Ok(())
        },
    );
}

#[test]
fn sphere_returns_the_nearest_positive_root() {
    bounded::assert_holds(
        [Axis::new(0.5, 5.0, 10), Axis::new(0.0, 1.0, 12)],
        |[r, u]| {
            // z0 in (−2r, −1e-3): in front of the vertex, inside the sphere's span.
            let z0 = -(1e-3 + u * (2.0 * r - 2e-3));
            let ray = Ray::new(Vec3::new(0.0, 0.0, z0), Vec3::new(0.0, 0.0, 1.0));
            match Sphere::new(1.0 / r).intersect(&ray, 1e-9) {
                Some(hit) if hit.point.z < -1e-9 || hit.point.z >= r => Err(format!(
                    "hit at z = {} is not the front surface",
                    hit.point.z
                )),
                _ => Ok(()),
            }
        },
    );
}

#[test]
fn lossless_split_conserves_energy_below_tir() {
    bounded::assert_holds(
        [
            Axis::new(1.0, 2.0, 11),
            Axis::new(0.0, 0.999, 8),
            Axis::new(0.0, 0.999, 12),
        ],
        |[n1, f, sin_t]| {
            let n2 = 1.0 + (n1 - 1.0) * f;
            let theta = sin_t.asin();
            match fresnel::critical_angle(n1, n2) {
                Some(crit) if theta < crit => {
                    let (rs, rp, ts, tp) = fresnel::fresnel_energy(n1, n2, theta);
                    if (rs + ts - 1.0).abs() >= 1e-9 {
                        return Err(format!("s: R + T = {}", rs + ts));
                    }
                    if (rp + tp - 1.0).abs() >= 1e-9 {
                        return Err(format!("p: R + T = {}", rp + tp));
                    }
                    Ok(())
                }
                _ => Ok(()),
            }
        },
    );
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
