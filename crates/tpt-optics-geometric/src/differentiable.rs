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
use crate::sequential::{
    ParamId, SequentialSystem, SurfaceAction, TraceContext, TraceResult, TunedParams,
};

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
        let Some(tuned) = TunedParams::<Dual>::seeded(system, id) else {
            grads.push(f64::NAN);
            continue;
        };
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

/// Paraxial ray `(y, u)` at the last surface for a ray that starts at height
/// 1 with zero slope, computed with dual numbers. The curvatures and gaps come
/// from `tuned`, and the indices from the system's media.
fn paraxial_exit(
    system: &SequentialSystem,
    ctx: TraceContext,
    tuned: &TunedParams<Dual>,
) -> (Dual, Dual) {
    let n = system.surfaces.len();
    let mut y = Dual::new(1.0);
    let mut u = Dual::new(0.0);
    for i in 0..n {
        let c = tuned.curvature[i];
        if system.surfaces[i].action == SurfaceAction::Reflect {
            u -= c * y * Dual::new(2.0);
        } else {
            let n_before = system.index_before(i, ctx);
            let n_after = system.index_after(i, ctx);
            let phi = c * Dual::new(n_after - n_before);
            u = (u * Dual::new(n_before) - y * phi) * Dual::new(1.0 / n_after);
        }
        if i + 1 < n {
            y += u * tuned.thickness[i];
        }
    }
    (y, u)
}

/// Exact derivatives of the paraxial effective focal length and back focal
/// length with respect to the parameters `params`.
///
/// Returns `(d EFL/dp, d BFL/dp)`, one entry per parameter, or `None` for a
/// system without a finite focus. The derivatives come from one dual-number
/// trace per parameter, so they are exact up to rounding.
pub fn paraxial_focus_gradients(
    system: &SequentialSystem,
    ctx: TraceContext,
    params: &[ParamId],
) -> Option<(Vec<f64>, Vec<f64>)> {
    let mut efl = Vec::with_capacity(params.len());
    let mut bfl = Vec::with_capacity(params.len());
    for &id in params {
        let tuned = TunedParams::<Dual>::seeded(system, id)?;
        let (y, u) = paraxial_exit(system, ctx, &tuned);
        if libm::fabs(u.v) < 1e-14 {
            return None;
        }
        // EFL = −1/u, so dEFL = du/u².
        efl.push(u.d / (u.v * u.v));
        // BFL = −y/u, so dBFL = −(dy·u − y·du)/u².
        bfl.push(-(y.d * u.v - y.v * u.d) / (u.v * u.v));
    }
    Some((efl, bfl))
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
    fn parameter_outside_the_system_gives_nan_not_a_panic() {
        // A ParamId that names no surface or gap has no derivative to give.
        // NaN marks it, so a caller cannot take it for a zero slope.
        let sys = singlet();
        let ray = Ray::new(Vec3::new(2e-3, 0.0, -0.1), Vec3::new(0.0, 0.0, 1.0));
        let ctx = TraceContext::default();
        let grads = gradient_of(&sys, &ray, ctx, &[ParamId::Curvature(99)], |res| {
            res.points.last().map_or(Dual::new(0.0), |p| p.x)
        });
        assert!(grads[0].is_nan(), "{grads:?}");
        assert!(TunedParams::<Dual>::seeded(&sys, ParamId::Thickness(99)).is_none());
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

#[cfg(test)]
mod paraxial_gradient_tests {
    use super::*;
    use tpt_optics_core::material::MaterialDatabase;
    use tpt_optics_core::{Kelvin, Wavelength};

    fn ctx() -> TraceContext {
        TraceContext {
            wavelength: Wavelength::from_nm(587.6),
            temperature: Kelvin::from_celsius(20.0),
        }
    }

    #[test]
    fn focus_gradients_match_central_differences() {
        let catalog = MaterialDatabase::presets();
        let mut sys = SequentialSystem::new();
        sys.push(
            crate::sequential::SeqSurface::aspheric(crate::surface::Asphere::sphere(10.0), 12.5e-3),
            5e-3,
            crate::sequential::Medium::Material(catalog.get("bk7").unwrap()),
        );
        sys.push(
            crate::sequential::SeqSurface::aspheric(
                crate::surface::Asphere::sphere(-10.0),
                12.5e-3,
            ),
            95e-3,
            crate::sequential::Medium::Vacuum,
        );
        let params = [
            ParamId::Curvature(0),
            ParamId::Curvature(1),
            ParamId::Thickness(0),
        ];
        let (efl_g, bfl_g) = paraxial_focus_gradients(&sys, ctx(), &params).unwrap();
        for (k, &id) in params.iter().enumerate() {
            let h = match id {
                ParamId::Curvature(_) => 1e-4,
                _ => 1e-8,
            };
            let shifted = |d: f64| {
                let mut s = sys.clone();
                match id {
                    ParamId::Curvature(i) => s.surfaces[i].asphere.curvature += d,
                    ParamId::Thickness(i) => s.thickness[i] += d,
                    _ => {}
                }
                s
            };
            let efl = |s: &SequentialSystem| s.paraxial(ctx()).effective_focal_length().unwrap();
            let bfl = |s: &SequentialSystem| s.paraxial(ctx()).back_focal_length().unwrap();
            let fd_e = (efl(&shifted(h)) - efl(&shifted(-h))) / (2.0 * h);
            let fd_b = (bfl(&shifted(h)) - bfl(&shifted(-h))) / (2.0 * h);
            assert!(
                (efl_g[k] - fd_e).abs() < 1e-5 * fd_e.abs().max(1e-6),
                "EFL {k}: {} vs {fd_e}",
                efl_g[k]
            );
            assert!(
                (bfl_g[k] - fd_b).abs() < 1e-5 * fd_b.abs().max(1e-6),
                "BFL {k}: {} vs {fd_b}",
                bfl_g[k]
            );
        }
    }
}
