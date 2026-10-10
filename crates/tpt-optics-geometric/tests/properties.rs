//! Property-based verification of optical laws (Phase 1 exit criteria).
//!
//! - Snell's law is reversible: refracting A→B then B→A along the same
//!   path returns the original direction.
//! - ABCD matrices: a refraction element has determinant n₁/n₂ and a full
//!   system in a unit-index environment is symplectic (det = 1).
//! - The vector form of Snell's law conserves tangential optical momentum.

use proptest::prelude::*;
use tpt_optics_core::Vec3;
use tpt_optics_geometric::abcd;
use tpt_optics_geometric::snell;

proptest! {
    #![proptest_config(ProptestConfig::with_cases(2_000))]

    #[test]
    fn snell_is_reversible(
        n1 in 1.0_f64..3.0,
        n2 in 1.0_f64..3.0,
        sin_theta in 0.0_f64..0.95,
        az in 0.0_f64..core::f64::consts::TAU,
    ) {
        // Incident direction on the z=0 plane at angle asin(sin_theta).
        let cos_t = (1.0_f64 - sin_theta * sin_theta).sqrt();
        let d = Vec3::new(sin_theta * az.cos(), sin_theta * az.sin(), -cos_t);
        let n = Vec3::new(0.0, 0.0, 1.0);
        // Skip TIR configurations (refract returns None).
        let Some(t1) = snell::refract(d, n, n1, n2) else { return Ok(()) };
        // Reverse pass: approach from the far side, so the normal flips.
        let Some(t2) = snell::refract(-t1, -n, n2, n1) else {
            panic!("asymmetric TIR: forward refracts, reverse does not");
        };
        prop_assert!((-t2).approx_eq(d, 1e-10), "forward {d:?} vs back {:?}'", -t2);
    }

    #[test]
    fn snell_conserves_tangential_momentum(
        n1 in 1.0_f64..3.0,
        n2 in 1.0_f64..3.0,
        sin_theta in 0.0_f64..0.95,
    ) {
        let cos_t = (1.0_f64 - sin_theta * sin_theta).sqrt();
        let d = Vec3::new(sin_theta, 0.0, -cos_t);
        let n = Vec3::new(0.0, 0.0, 1.0);
        let Some(t) = snell::refract(d, n, n1, n2) else { return Ok(()) };
        prop_assert!((t.x * n2 - d.x * n1).abs() < 1e-12);
    }

    #[test]
    fn refraction_abcd_det_is_index_ratio(
        n1 in 1.0_f64..3.0,
        n2 in 1.0_f64..3.0,
        c in -10.0_f64..10.0,
    ) {
        let m = abcd::refracting_surface::<f64>(n1, n2, c);
        prop_assert!((m.det() - n1 / n2).abs() < 1e-12);
    }

    #[test]
    fn unit_index_systems_are_symplectic(
        f1 in 0.05_f64..2.0,
        f2 in 0.05_f64..2.0,
        sep in 0.0_f64..3.0,
    ) {
        // Arbitrary compositions of lenses and free space in air: det = 1.
        let m = abcd::thin_lens::<f64>(f1)
            .mul_matrix(abcd::free_space(sep))
            .mul_matrix(abcd::thin_lens(f2))
            .mul_matrix(abcd::free_space(sep * 0.5));
        prop_assert!((m.det() - 1.0).abs() < 1e-9, "det = {}", m.det());
    }

    #[test]
    fn reflection_preserves_direction_norm(
        dx in -1.0_f64..1.0,
        dy in -1.0_f64..1.0,
        dz in -1.0_f64..0.0,
    ) {
        let d = Vec3::new(dx, dy, dz);
        let Some(d) = d.normalized() else { return Ok(()) };
        let n = Vec3::new(0.0, 0.0, 1.0);
        let r = snell::reflect(d, n);
        prop_assert!((r.norm() - 1.0).abs() < 1e-12);
        // Incidence = reflection angle.
        prop_assert!((d.dot(n) + r.dot(n)).abs() < 1e-12);
    }

    #[test]
    fn mirror_abcd_is_symplectic(c in -10.0_f64..10.0) {
        let m = abcd::mirror::<f64>(c);
        prop_assert!((m.det() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn free_space_q_transform_adds_distance(d in 0.0_f64..5.0, q_re in -2.0_f64..2.0) {
        use tpt_optics_core::Complex;
        let m = abcd::free_space::<f64>(d);
        let q = Complex::real(q_re);
        prop_assert!(m.apply_q(q).approx_eq(Complex::real(q_re + d), 1e-10));
    }
}
