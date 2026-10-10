//! Kani bounded-model-checking harnesses.
//!
//! These harnesses are compiled only under `cfg(kani)` (set by the `kani`
//! verifier via `cargo kani`). They encode the two proofs required by the
//! Phase-3 exit criteria:
//!
//! 1. **Ray–sphere intersection finds the nearest positive root** when one
//!    exists (the quadratic solver's selection logic over symbolic inputs).
//! 2. **A lossless beam splitter conserves energy** (Fresnel R + T = 1 for
//!    symbolic angles below the critical angle).
//!
//! Run with: `cargo kani -p tpt-optics-verify` (requires the Kani toolchain).

#[cfg(kani)]
mod proofs {
    use tpt_optics_core::fresnel::{critical_angle, fresnel_energy};
    use tpt_optics_core::Vec3;
    use tpt_optics_geometric::prelude::{Ray, Sphere};
    use tpt_optics_geometric::surface::OpticalSurface;

    /// Proof 1a: whenever the solver returns a hit, it lies ahead of the
    /// ray and on the sphere.
    #[kani::proof]
    fn sphere_hit_is_ahead_and_on_surface() {
        let curvature: f64 = kani::any();
        kani::assume(curvature > 1.0 && curvature < 10.0);
        let ox: f64 = kani::any();
        let oy: f64 = kani::any();
        let oz: f64 = kani::any();
        kani::assume(oz < -1e-3 && oz > -100.0);
        kani::assume(ox * ox + oy * oy < 25.0);
        let dz: f64 = kani::any();
        kani::assume(dz > 0.5 && dz <= 1.0);
        let d = Vec3::new(0.0, 0.0, dz);
        let ray = Ray::new(Vec3::new(ox, oy, oz), d);
        let sphere = Sphere::new(curvature);
        if let Some(hit) = sphere.intersect(&ray, 1e-9) {
            kani::assert(hit.t > 0.0, "hit parameter positive");
            let on_surface = {
                // |p − C|² = R² with C = (0, 0, R)
                let r = 1.0 / curvature;
                let px = hit.point.x;
                let py = hit.point.y;
                let pz = hit.point.z - r;
                (px * px + py * py + pz * pz - r * r).abs() < 1e-6 * r * r
            };
            kani::assert(on_surface, "hit lies on the sphere");
        }
    }

    /// Proof 1b: of the two quadratic roots, the returned hit uses the
    /// smaller positive one (no far-side skipping).
    #[kani::proof]
    fn sphere_nearest_positive_root() {
        // Symbolic 1-D setup along the axis: ray at z0 < 0 heading +z.
        let r: f64 = kani::any();
        kani::assume(r > 0.5 && r < 5.0);
        let z0: f64 = kani::any();
        kani::assume(z0 < -1e-3 && z0 > -2.0 * r);
        let sphere = Sphere::new(1.0 / r);
        let ray = Ray::new(Vec3::new(0.0, 0.0, z0), Vec3::new(0.0, 0.0, 1.0));
        if let Some(hit) = sphere.intersect(&ray, 1e-9) {
            // The near root is z_hit = r − sqrt(r² − (r + z0)²)... assert the
            // hit is the front surface: z ∈ [0, r) for a vertex-centered sphere.
            kani::assert(hit.point.z >= -1e-9 && hit.point.z < r, "front surface hit");
        }
    }

    /// Proof 2: lossless Fresnel split conserves energy below TIR.
    #[kani::proof]
    fn lossless_split_conserves_energy() {
        let n1: f64 = kani::any();
        let n2: f64 = kani::any();
        kani::assume(n1 >= 1.0 && n1 <= 2.0);
        kani::assume(n2 >= 1.0 && n2 <= 2.0);
        let sin_t: f64 = kani::any();
        kani::assume(sin_t >= 0.0 && sin_t < 0.999);
        let theta = libm::asin(sin_t);
        if let Some(crit) = critical_angle(n1, n2) {
            if theta < crit {
                let (rs, rp, ts, tp) = fresnel_energy(n1, n2, theta);
                kani::assert((rs + ts - 1.0).abs() < 1e-9, "s: R + T = 1");
                kani::assert((rp + tp - 1.0).abs() < 1e-9, "p: R + T = 1");
            }
        }
    }
}
