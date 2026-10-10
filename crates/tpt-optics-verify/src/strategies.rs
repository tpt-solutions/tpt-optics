//! proptest strategies for optical quantities.

use proptest::prelude::*;

use tpt_optics_core::{Complex, Vec3};
use tpt_optics_geometric::prelude::Ray;

/// Finite, reasonable f64 in an inclusive range.
pub fn f64_in(lo: f64, hi: f64) -> impl Strategy<Value = f64> {
    (lo..=hi).prop_filter("finite", |v| v.is_finite())
}

/// Unit 3-direction with z component negative (propagating toward +z
/// surfaces from the −z side): `sinθ ∈ [0, 0.95]`, uniform azimuth.
pub fn incident_direction() -> impl Strategy<Value = Vec3<f64>> {
    (0.0_f64..0.95, 0.0_f64..core::f64::consts::TAU).prop_map(|(sin_t, az)| {
        let cos_t = (1.0 - sin_t * sin_t).sqrt();
        Vec3::new(sin_t * az.cos(), sin_t * az.sin(), -cos_t)
    })
}

/// Ray starting slightly before the origin plane, unit direction.
pub fn incident_ray() -> impl Strategy<Value = Ray<f64>> {
    incident_direction().prop_map(|d| Ray::new(Vec3::new(0.0, 0.0, -0.01), d))
}

/// Refractive index in [1.0, 3.5] (covers gases to semiconductors).
pub fn refractive_index() -> impl Strategy<Value = f64> {
    f64_in(1.0, 3.5)
}

/// Pair (n1, n2) of indices, n1 > n2 for TIR-capable geometries.
pub fn index_pair_descending() -> impl Strategy<Value = (f64, f64)> {
    (f64_in(1.0, 3.5), f64_in(0.0_f64, 1.0)).prop_map(|(n1, frac)| {
        let n2 = n1 * (0.2 + 0.8 * frac); // n2 ∈ (0.2, 1.0]·n1
        (n1, n2)
    })
}

/// Curvature in ±20 1/m (radii ≥ 50 mm).
pub fn curvature() -> impl Strategy<Value = f64> {
    f64_in(-20.0, 20.0)
}

/// Conic constant in [−10, 2].
pub fn conic_constant() -> impl Strategy<Value = f64> {
    f64_in(-10.0, 2.0)
}

/// Unit complex number (fully polarized state).
pub fn unit_complex() -> impl Strategy<Value = Complex<f64>> {
    (0.0_f64..core::f64::consts::TAU).prop_map(Complex::cis)
}

/// Sellmeier-style coefficients producing a physical index (B ≥ 0,
/// distinct resonances).
pub fn sellmeier_coeffs() -> impl Strategy<Value = ([f64; 3], [f64; 3])> {
    (
        f64_in(0.0, 2.0),
        f64_in(0.0, 2.0),
        f64_in(0.0, 2.0),
        f64_in(0.001, 0.05),
        f64_in(0.05, 0.3),
        f64_in(30.0, 200.0),
    )
        .prop_map(|(b1, b2, b3, c1, c2, c3)| ([b1, b2, b3], [c1, c2, c3]))
}
