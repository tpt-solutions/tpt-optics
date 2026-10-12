//! Lens merit functions over sequential ray tracing.
//!
//! Each function takes a [`SequentialSystem`](tpt_optics_geometric::sequential::SequentialSystem)
//! plus a trace context and reduces it to a scalar (or a small vector for
//! least-squares). Variables are handled by the caller: rebuild the system
//! from a parameter vector inside the residual closure.

use alloc::vec::Vec;

use tpt_optics_geometric::sequential::{rms_spot_radius, SequentialSystem, TraceContext};

/// RMS transverse spot radius (meters) over a pupil fan of fractional
/// heights (0..1).
pub fn rms_spot(
    system: &SequentialSystem,
    ctx: TraceContext,
    field_angle_deg: f64,
    fractions: &[f64],
    pupil_semi_diameter: f64,
) -> f64 {
    let heights: Vec<f64> = fractions.iter().map(|f| f * pupil_semi_diameter).collect();
    rms_spot_radius(system, ctx, field_angle_deg, &heights)
}

/// Chromatic focal shift: axial image-plane z that minimizes the spot for
/// each wavelength, returned as (max − min) shift across the band.
///
/// Uses the paraxial back focal length per wavelength as a fast proxy.
pub fn chromatic_focal_shift(system: &SequentialSystem, contexts: &[TraceContext]) -> f64 {
    let mut bfls: Vec<f64> = Vec::with_capacity(contexts.len());
    for &ctx in contexts {
        if let Some(b) = system.paraxial(ctx).back_focal_length() {
            bfls.push(b);
        }
    }
    if bfls.len() < 2 {
        return 0.0;
    }
    let max = bfls.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let min = bfls.iter().copied().fold(f64::INFINITY, f64::min);
    max - min
}

/// Relative distortion at a field angle: the real chief-ray image height
/// versus the paraxial chief-ray height, `(y_img − y_ideal)/y_ideal`.
///
/// The paraxial reference is `B·tan(θ)` from the system ABCD matrix, for a
/// chief ray launched from the first surface (the stop in these designs).
/// This does not assume the image plane sits at the focus, and it accounts
/// for the principal-plane offset of thick lenses. The real chief ray is
/// traced with [`field_ray`](tpt_optics_geometric::sequential::field_ray).
pub fn distortion(system: &SequentialSystem, ctx: TraceContext, field_angle_deg: f64) -> f64 {
    let theta = field_angle_deg.to_radians();
    let abcd = system.paraxial(ctx).abcd();
    // Chief ray at the field angle, crossing the first vertex plane on axis.
    let ray = tpt_optics_geometric::sequential::field_ray(0.0, field_angle_deg);
    let res = system.trace(&ray, ctx);
    if res.blocked {
        return f64::NAN;
    }
    let Some(point) = res.points.last() else {
        return f64::NAN;
    };
    let y_ideal = abcd.b * libm::tan(theta);
    if y_ideal.abs() < 1e-15 {
        return 0.0;
    }
    (point.y - y_ideal) / y_ideal
}

/// Aggregate spot merit: squared RMS spot radii over the given fields,
/// ready for least-squares (one residual per field).
pub fn spot_residuals(
    system: &SequentialSystem,
    ctx: TraceContext,
    fields_deg: &[f64],
    fractions: &[f64],
    pupil_semi_diameter: f64,
) -> Vec<f64> {
    fields_deg
        .iter()
        .map(|&f| rms_spot(system, ctx, f, fractions, pupil_semi_diameter))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::{boxed::Box, sync::Arc};
    use tpt_optics_core::{dispersion::Sellmeier, material::Material, Kelvin, Wavelength};
    use tpt_optics_geometric::sequential::{Medium, SeqSurface, SequentialSystem};
    use tpt_optics_geometric::surface::Asphere;

    fn ctx(nm: f64) -> TraceContext {
        TraceContext {
            wavelength: Wavelength::from_nm(nm),
            temperature: Kelvin::from_celsius(20.0),
        }
    }

    /// Simple singlet; `dispersive` selects real BK7 instead of a constant
    /// index so chromatic effects exist.
    fn singlet(dispersive: bool) -> SequentialSystem {
        let mut sys = SequentialSystem::new();
        let glass = if dispersive {
            Medium::Material(Arc::new(Material::new(
                "bk7",
                Box::new(Sellmeier::bk7()),
                0.0,
                Kelvin::from_celsius(20.0),
            )))
        } else {
            Medium::Const(1.5168)
        };
        sys.push(
            SeqSurface::aspheric(Asphere::sphere(10.0), 12e-3),
            5e-3,
            glass,
        );
        sys.push(
            SeqSurface::aspheric(Asphere::sphere(-10.0), 40e-3),
            40e-3,
            Medium::Const(1.0),
        );
        sys.push(SeqSurface::flat(12e-3), 0.0, Medium::Vacuum);
        sys
    }

    #[test]
    fn spot_minimized_at_best_focus_via_dls() {
        // Optimize the image-plane distance for minimum RMS spot with the
        // damped-least-squares optimizer: the integration of optimizer +
        // merit on a real lens.
        let build = |image_dist: f64| -> SequentialSystem {
            let mut sys = singlet(false);
            sys.thickness[1] = image_dist;
            sys
        };
        // Paraxial seed: back focal length measured from the last lens
        // surface (zero the image-plane offset for the reference).
        let mut seed_sys = singlet(false);
        seed_sys.thickness[1] = 0.0;
        let bfl = seed_sys.paraxial(ctx(587.6)).back_focal_length().unwrap();
        let r = |x: &[f64], out: &mut Vec<f64>| {
            let sys = build(x[0]);
            out[0] = rms_spot(&sys, ctx(587.6), 0.0, &[0.0, 0.3, 0.6, 0.8], 2e-3);
        };
        let res = crate::optimize::damped_least_squares(&r, &[bfl], 1, None, 25);
        assert!(
            res.f < 100e-6,
            "best-focus rms = {} at x = {:?}",
            res.f,
            res.x
        );
        // Best focus sits close to the paraxial back focus (within the
        // few-mm scale of primary spherical aberration).
        assert!(
            (res.x[0] - bfl).abs() < 5e-3,
            "best focus {} vs bfl {bfl}",
            res.x[0]
        );
    }

    #[test]
    fn chromatic_shift_positive_for_bk7_singlet() {
        let sys = singlet(true);
        let shift = chromatic_focal_shift(&sys, &[ctx(486.1), ctx(587.6), ctx(656.3)]);
        // BK7 singlet: focus shift of order 0.5–2 mm.
        assert!(shift > 100e-6 && shift < 5e-3, "shift = {shift}");
        // And a constant-index singlet has none by construction.
        assert_eq!(
            chromatic_focal_shift(&singlet(false), &[ctx(486.1), ctx(587.6), ctx(656.3)]),
            0.0
        );
    }

    #[test]
    fn distortion_finite_and_small_on_axis() {
        let sys = singlet(false);
        let d0 = distortion(&sys, ctx(587.6), 0.0);
        assert!(d0.abs() < 1e-9, "on-axis distortion {d0}");
        let d5 = distortion(&sys, ctx(587.6), 2.0);
        assert!(d5.is_finite());
        assert!(d5.abs() < 0.5, "distortion {d5}");
    }

    #[test]
    fn spot_residuals_shape() {
        let sys = singlet(false);
        let res = spot_residuals(&sys, ctx(587.6), &[0.0, 1.0], &[0.0, 0.5], 6e-3);
        assert_eq!(res.len(), 2);
        assert!(res.iter().all(|v| v.is_finite()));
    }
}

#[cfg(test)]
mod distortion_focus_tests {
    use super::*;
    use tpt_optics_core::{Kelvin, Wavelength};
    use tpt_optics_geometric::sequential::{Medium, SeqSurface, SequentialSystem, TraceContext};
    use tpt_optics_geometric::surface::Asphere;

    #[test]
    fn distortion_is_small_when_the_lens_is_focused() {
        // Distortion is measured against f·tanθ, which assumes the image plane
        // sits at the focus. Place the image at the paraxial back focal length.
        let ctx = TraceContext {
            wavelength: Wavelength::from_nm(587.6),
            temperature: Kelvin::from_celsius(20.0),
        };
        let mut sys = SequentialSystem::new();
        sys.push(
            SeqSurface::aspheric(Asphere::sphere(10.0), 12e-3),
            5e-3,
            Medium::Const(1.5168),
        );
        sys.push(
            SeqSurface::aspheric(Asphere::sphere(-10.0), 40e-3),
            40e-3,
            Medium::Const(1.0),
        );
        sys.push(SeqSurface::flat(12e-3), 0.0, Medium::Vacuum);
        let bfl = sys.paraxial(ctx).back_focal_length().expect("BFL");
        sys.thickness[1] = bfl;
        let d = distortion(&sys, ctx, 2.0);
        assert!(d.is_finite());
        let efl = sys.paraxial(ctx).effective_focal_length();
        assert!(
            d.abs() < 0.05,
            "focused distortion at 2°: {d}, bfl {bfl}, efl {efl:?}"
        );
    }
}
