//! Differentiable ray tracing: exact gradients of traced quantities.
//!
//! tpt-optics traces are generic over [`Scalar`], so seeding a single design
//! variable as a [`Dual`] number turns the whole trace into the Jacobian
//! column for that variable. [`gradient_of`] orchestrates this: give it the
//! parameters to differentiate and a projection from the trace result to the
//! quantity of interest (e.g. the final intersection x, or the OPD), and it
//! returns exact first derivatives.
//!
//! This is the hook used by `tpt-optics-design` for damped-least-squares
//! optimisation and by downstream differentiable-rendering integrations.

use alloc::vec::Vec;

use tpt_optics_core::Dual;

use crate::ray::Ray;
use crate::sequential::{ParamId, SequentialSystem, TraceContext, TraceResult, TunedParams};

/// Exact gradient of a scalar traced quantity with respect to the listed
/// parameters, via one dual-number trace per parameter.
///
/// `extract` projects the (dual) trace result to the quantity of interest.
/// It must be finite for the provided ray; NaN derivatives indicate the
/// projection was evaluated on a blocked trace.
pub fn gradient_of<F>(
    system: &SequentialSystem,
    ray: &Ray<f64>,
    ctx: TraceContext,
    params: &[ParamId],
    extract: F,
) -> Vec<f64>
where
    F: Fn(&TraceResult<Dual>) -> Dual,
{
    let mut grads = Vec::with_capacity(params.len());
    for &id in params {
        let tuned = TunedParams::<Dual>::seeded(system, id);
        // The ray itself enters as constants (derivative 0); only the
        // geometry is seeded.
        let dual_ray = Ray::new(ray.origin.map(Dual::new), ray.direction.map(Dual::new));
        let result = system.trace_with(&tuned, &dual_ray, ctx);
        let q = extract(&result);
        grads.push(q.d);
    }
    grads
}

/// Gradient of the x-coordinate of the final intersection.
pub fn gradient_final_x(
    system: &SequentialSystem,
    ray: &Ray<f64>,
    ctx: TraceContext,
    params: &[ParamId],
) -> Vec<f64> {
    gradient_of(system, ray, ctx, params, |r| {
        r.points.last().map_or(Dual::new(f64::NAN), |p| p.x)
    })
}

/// Gradient of the accumulated optical path length.
pub fn gradient_opd(
    system: &SequentialSystem,
    ray: &Ray<f64>,
    ctx: TraceContext,
    params: &[ParamId],
) -> Vec<f64> {
    gradient_of(system, ray, ctx, params, |r| r.opd)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequential::{Medium, SeqSurface, SequentialSystem};
    use crate::surface::Asphere;
    use tpt_optics_core::Vec3;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    fn singlet() -> SequentialSystem {
        let mut sys = SequentialSystem::new();
        sys.push(
            SeqSurface::aspheric(Asphere::sphere(10.0), 12e-3),
            5e-3,
            Medium::Const(1.5168),
        );
        sys.push(
            SeqSurface::aspheric(Asphere::sphere(-10.0), 12e-3),
            30e-3,
            Medium::Const(1.0),
        );
        sys.push(SeqSurface::flat(12e-3), 0.0, Medium::Vacuum);
        sys
    }

    #[test]
    fn focal_shift_vs_curvature_matches_finite_difference() {
        let sys = singlet();
        let ctx = TraceContext::default();
        let ray = Ray::new(Vec3::new(1e-3, 0.0, -0.01), Vec3::new(0.0, 0.0, 1.0));

        // Quantity: x of the final intersection on the (fixed) image plane.
        // d(x_image)/d(c1): compare AD against central differences.
        let ad = gradient_final_x(&sys, &ray, ctx, &[ParamId::Curvature(0)])[0];

        let h = 1e-7;
        let mut sys_p = sys.clone();
        sys_p.surfaces[0].asphere.curvature += h;
        let mut sys_m = sys.clone();
        sys_m.surfaces[0].asphere.curvature -= h;
        let xp = sys_p.trace(&ray, ctx).points.last().unwrap().x;
        let xm = sys_m.trace(&ray, ctx).points.last().unwrap().x;
        let fd = (xp - xm) / (2.0 * h);

        assert!(
            close(ad, fd, 1e-5 * fd.abs().max(1e-12)),
            "ad {ad} vs fd {fd}"
        );
    }

    #[test]
    fn opd_gradient_vs_thickness_fd() {
        let sys = singlet();
        let ctx = TraceContext::default();
        let ray = Ray::new(Vec3::new(0.5e-3, 0.0, -0.01), Vec3::new(0.0, 0.0, 1.0));
        let ad = gradient_opd(&sys, &ray, ctx, &[ParamId::Thickness(0)])[0];

        let h = 1e-9;
        let mut sys_p = sys.clone();
        sys_p.thickness[0] += h;
        let mut sys_m = sys.clone();
        sys_m.thickness[0] -= h;
        let fd = (sys_p.trace(&ray, ctx).opd - sys_m.trace(&ray, ctx).opd) / (2.0 * h);
        assert!(close(ad, fd, 1e-6), "ad {ad} vs fd {fd}");
    }

    #[test]
    fn opd_gradient_of_thickness_is_glass_index() {
        // thickness[0] is the glass gap after surface 0; stretching it moves
        // every later surface rigidly with it, so dOPD/dt0 = n_glass.
        let sys = singlet();
        let ctx = TraceContext::default();
        let ray = Ray::new(Vec3::new(0.0, 0.0, -0.01), Vec3::new(0.0, 0.0, 1.0));
        let ad = gradient_opd(&sys, &ray, ctx, &[ParamId::Thickness(0)])[0];
        assert!(close(ad, 1.5168, 1e-9), "ad = {ad}");
    }
}
