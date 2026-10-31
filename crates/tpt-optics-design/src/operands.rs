//! Operand-based merit functions: each operand is one weighted residual, and
//! the residuals of a list of operands form a least-squares problem.
//!
//! An operand reads one quantity of a sequential system and compares it with
//! a target or a range: focal length, back focal length, spot size, a
//! thickness, the height of a ray at a surface, or the edge thickness of a
//! gap. Build the residuals from a parameter vector with [`optimize_system`],
//! which runs Levenberg–Marquardt on them.

use alloc::vec::Vec;

use tpt_optics_geometric::sequential::{
    field_ray, rms_spot_radius, SequentialSystem, TraceContext,
};
use tpt_optics_geometric::surface::OpticalSurface;

use crate::optimize::{damped_least_squares, OptimizeResult};

/// One requirement on a sequential system.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Operand {
    /// Effective focal length (m) should equal `target`.
    Efl {
        /// Target focal length (m).
        target: f64,
        /// Weight of this residual.
        weight: f64,
    },
    /// Back focal length (m) should equal `target`.
    Bfl {
        /// Target back focal length (m).
        target: f64,
        /// Weight of this residual.
        weight: f64,
    },
    /// RMS spot radius (m) at `field_deg` should equal `target`.
    RmsSpot {
        /// Field angle (degrees).
        field_deg: f64,
        /// Pupil ray heights (m).
        heights: Vec<f64>,
        /// Target spot radius (m).
        target: f64,
        /// Weight of this residual.
        weight: f64,
    },
    /// Thickness `system.thickness[index]` (m) should lie in `[min, max]`.
    /// Inside the range the residual is zero.
    ThicknessRange {
        /// Index of the gap.
        index: usize,
        /// Lower bound (m).
        min: f64,
        /// Upper bound (m).
        max: f64,
        /// Weight of this residual.
        weight: f64,
    },
    /// Height (m) of a ray at surface `surface` should equal `target`. The
    /// ray starts at `pupil_height` on the field `field_deg`. `axis` picks
    /// the transverse coordinate read from the ray.
    RayHeight {
        /// Index of the surface.
        surface: usize,
        /// Field angle (degrees).
        field_deg: f64,
        /// Ray height at the first surface (m), along the `x` pupil axis.
        pupil_height: f64,
        /// Target height at the surface (m).
        target: f64,
        /// Coordinate read from the ray: [`RayAxis::X`] (sagittal, the
        /// default) or [`RayAxis::Y`] (meridional, where the field lies).
        #[cfg_attr(feature = "serde", serde(default))]
        axis: RayAxis,
        /// Weight of this residual.
        weight: f64,
    },
    /// Edge thickness of gap `index` should lie in `[min, max]` (m). The edge
    /// is taken at the smaller of the two clear semi-diameters, so it is the
    /// thinnest point of the element.
    EdgeThickness {
        /// Index of the gap.
        index: usize,
        /// Lower bound (m).
        min: f64,
        /// Upper bound (m).
        max: f64,
        /// Weight of this residual.
        weight: f64,
    },
}

/// Edge thickness of gap `index` (m): the centre thickness plus the sag of
/// the rear surface minus the sag of the front surface, both at the smaller
/// of their clear semi-diameters. `None` for a gap with no following surface.
#[must_use]
pub fn edge_thickness(system: &SequentialSystem, index: usize) -> Option<f64> {
    let front = system.surfaces.get(index)?;
    let rear = system.surfaces.get(index + 1)?;
    let t = *system.thickness.get(index)?;
    let r = front.semi_diameter.min(rear.semi_diameter);
    Some(t + rear.asphere.sag(r, 0.0) - front.asphere.sag(r, 0.0))
}

/// Residual of one operand. A quantity that does not exist for the system
/// (for example a focal length for an afocal system, or a ray that is
/// blocked) gives a penalty of `weight · max(|target|, 1)`, so the optimiser
/// moves away from it instead of seeing a NaN.
#[must_use]
pub fn residual(op: &Operand, system: &SequentialSystem, ctx: TraceContext) -> f64 {
    let penalty = |weight: f64, target: f64| weight * libm::fabs(target).max(1.0);
    match op {
        Operand::Efl { target, weight } => match system.paraxial(ctx).effective_focal_length() {
            Some(v) => weight * (v - target),
            None => penalty(*weight, *target),
        },
        Operand::Bfl { target, weight } => match system.paraxial(ctx).back_focal_length() {
            Some(v) => weight * (v - target),
            None => penalty(*weight, *target),
        },
        Operand::RmsSpot {
            field_deg,
            heights,
            target,
            weight,
        } => {
            let spot = rms_spot_radius(system, ctx, *field_deg, heights);
            if spot.is_finite() {
                weight * (spot - target)
            } else {
                penalty(*weight, *target)
            }
        }
        Operand::ThicknessRange {
            index,
            min,
            max,
            weight,
        } => range_residual(
            system.thickness.get(*index).copied(),
            *min,
            *max,
            *weight,
            penalty,
        ),
        Operand::RayHeight {
            surface,
            field_deg,
            pupil_height,
            target,
            axis,
            weight,
        } => {
            let trace = system.trace(&field_ray(*pupil_height, *field_deg), ctx);
            match trace.points.get(*surface) {
                Some(p) if !trace.blocked => weight * (axis.pick(p.x, p.y) - target),
                _ => penalty(*weight, *target),
            }
        }
        Operand::EdgeThickness {
            index,
            min,
            max,
            weight,
        } => range_residual(edge_thickness(system, *index), *min, *max, *weight, penalty),
    }
}

/// Transverse coordinate of a ray read by [`Operand::RayHeight`].
///
/// Field angles tilt the ray in the y–z plane, so `Y` is the meridional
/// height (the chief ray has a nonzero `y` off axis) and `X` is the sagittal
/// one. A pupil offset of zero leaves `X` at zero for every field.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum RayAxis {
    /// Sagittal (x) height. This is the default, matching earlier releases.
    #[default]
    X,
    /// Meridional (y) height, where the field angle acts.
    Y,
}

impl RayAxis {
    /// The coordinate of a point along this axis.
    #[must_use]
    pub fn pick<T>(self, x: T, y: T) -> T {
        match self {
            Self::X => x,
            Self::Y => y,
        }
    }
}

/// Residual of a value that should lie in `[min, max]`. A missing value gets
/// the penalty for `max`.
fn range_residual(
    value: Option<f64>,
    min: f64,
    max: f64,
    weight: f64,
    penalty: impl Fn(f64, f64) -> f64,
) -> f64 {
    match value {
        Some(t) if t < min => weight * (min - t),
        Some(t) if t > max => weight * (t - max),
        Some(_) => 0.0,
        None => penalty(weight, max),
    }
}

/// Residual vector of a list of operands, one entry per operand.
#[must_use]
pub fn residuals(ops: &[Operand], system: &SequentialSystem, ctx: TraceContext) -> Vec<f64> {
    ops.iter().map(|op| residual(op, system, ctx)).collect()
}

/// Penalty of an operand when the system cannot be built at all.
fn build_penalty(op: &Operand) -> f64 {
    match op {
        Operand::Efl { target, weight }
        | Operand::Bfl { target, weight }
        | Operand::RmsSpot { target, weight, .. }
        | Operand::RayHeight { target, weight, .. } => weight * libm::fabs(*target).max(1.0),
        Operand::ThicknessRange { max, weight, .. }
        | Operand::EdgeThickness { max, weight, .. } => weight * libm::fabs(*max).max(1.0),
    }
}

/// Optimise a parameter vector so that the built system meets the operands.
///
/// `build` maps parameters to a system. Returning `None` (for example for an
/// invalid prescription) gives every residual the penalty of its operand.
pub fn optimize_system<B>(
    build: &B,
    x0: &[f64],
    ops: &[Operand],
    ctx: TraceContext,
    iterations: usize,
) -> OptimizeResult
where
    B: Fn(&[f64]) -> Option<SequentialSystem>,
{
    let r = |x: &[f64], out: &mut Vec<f64>| {
        out.clear();
        match build(x) {
            Some(sys) => out.extend(residuals(ops, &sys, ctx)),
            None => out.extend(ops.iter().map(build_penalty)),
        }
    };
    damped_least_squares(&r, x0, ops.len(), None, iterations)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_optics_core::material::MaterialDatabase;
    use tpt_optics_core::{Kelvin, Wavelength};
    use tpt_optics_geometric::presets::singlet;

    fn ctx() -> TraceContext {
        TraceContext {
            wavelength: Wavelength::from_nm(587.6),
            temperature: Kelvin::from_celsius(20.0),
        }
    }

    fn build(catalog: &MaterialDatabase, x: &[f64]) -> Option<SequentialSystem> {
        singlet(catalog, x[0], x[1], 12.5e-3, 5e-3, "bk7", 95e-3).ok()
    }

    #[test]
    fn efl_residual_is_zero_at_the_exact_focal_length() {
        let catalog = MaterialDatabase::presets();
        let sys = build(&catalog, &[100e-3, -100e-3]).unwrap();
        let efl = sys.paraxial(ctx()).effective_focal_length().unwrap();
        let r = residual(
            &Operand::Efl {
                target: efl,
                weight: 1.0,
            },
            &sys,
            ctx(),
        );
        assert!(libm::fabs(r) < 1e-12, "{r}");
    }

    #[test]
    fn thickness_inside_its_range_has_no_residual() {
        let catalog = MaterialDatabase::presets();
        let sys = build(&catalog, &[100e-3, -100e-3]).unwrap();
        let op = Operand::ThicknessRange {
            index: 0,
            min: 2e-3,
            max: 8e-3,
            weight: 1.0,
        };
        assert_eq!(residual(&op, &sys, ctx()), 0.0);
        let tight = Operand::ThicknessRange {
            index: 0,
            min: 6e-3,
            max: 8e-3,
            weight: 2.0,
        };
        assert!((residual(&tight, &sys, ctx()) - 2.0 * 1e-3).abs() < 1e-12);
    }

    #[test]
    fn edge_thickness_matches_the_sag_formula() {
        // Front R = 100 mm (c = 10/m), rear R = −100 mm, edge at 12.5 mm:
        // sag = c·r²/(1 + √(1 − c²r²)) ≈ 0.7855 mm each way.
        let catalog = MaterialDatabase::presets();
        let sys = build(&catalog, &[100e-3, -100e-3]).unwrap();
        let r = 12.5e-3;
        let sag = |c: f64| c * r * r / (1.0 + libm::sqrt(1.0 - c * c * r * r));
        let expected = 5e-3 + sag(-10.0) - sag(10.0);
        let got = edge_thickness(&sys, 0).unwrap();
        assert!((got - expected).abs() < 1e-12, "{got} vs {expected}");
        assert!(got < 5e-3, "a biconvex edge is thinner than the centre");
    }

    #[test]
    fn edge_thickness_range_penalises_a_thin_edge() {
        let catalog = MaterialDatabase::presets();
        let sys = build(&catalog, &[100e-3, -100e-3]).unwrap();
        let edge = edge_thickness(&sys, 0).unwrap();
        let op = Operand::EdgeThickness {
            index: 0,
            min: edge + 1e-3,
            max: 10e-3,
            weight: 1.0,
        };
        assert!((residual(&op, &sys, ctx()) - 1e-3).abs() < 1e-12);
    }

    #[test]
    fn ray_height_residual_is_the_offset_from_the_target() {
        let catalog = MaterialDatabase::presets();
        let sys = build(&catalog, &[100e-3, -100e-3]).unwrap();
        let trace = sys.trace(&field_ray(2e-3, 0.0), ctx());
        let x_front = trace.points[0].x;
        let op = Operand::RayHeight {
            surface: 0,
            field_deg: 0.0,
            pupil_height: 2e-3,
            target: 0.0,
            axis: RayAxis::X,
            weight: 3.0,
        };
        assert!((residual(&op, &sys, ctx()) - 3.0 * x_front).abs() < 1e-12);
    }

    #[test]
    fn chief_ray_height_needs_the_meridional_axis() {
        // The chief ray has no pupil offset, so its sagittal (x) height is
        // zero for every field. The meridional (y) height moves with the
        // field angle.
        let catalog = MaterialDatabase::presets();
        let sys = build(&catalog, &[100e-3, -100e-3]).unwrap();
        let trace = sys.trace(&field_ray(0.0, 3.0), ctx());
        assert!(!trace.blocked, "the test chief ray must reach the image");
        let y_image = trace.points[1].y;
        assert!(y_image.abs() > 1e-4, "chief ray should be off axis");
        let mut op = Operand::RayHeight {
            surface: 1,
            field_deg: 3.0,
            pupil_height: 0.0,
            target: 0.0,
            axis: RayAxis::X,
            weight: 2.0,
        };
        assert!(residual(&op, &sys, ctx()).abs() < 1e-15);
        if let Operand::RayHeight { axis, .. } = &mut op {
            *axis = RayAxis::Y;
        }
        assert!((residual(&op, &sys, ctx()) - 2.0 * y_image).abs() < 1e-12);
    }

    #[test]
    fn missing_focal_length_is_penalised_not_nan() {
        // Two flat surfaces in air: afocal, so there is no focal length.
        let mut sys = SequentialSystem::new();
        sys.push(
            tpt_optics_geometric::sequential::SeqSurface::flat(12e-3),
            5e-3,
            tpt_optics_geometric::sequential::Medium::Const(1.0),
        );
        sys.push(
            tpt_optics_geometric::sequential::SeqSurface::flat(12e-3),
            0.0,
            tpt_optics_geometric::sequential::Medium::Vacuum,
        );
        let r = residual(
            &Operand::Efl {
                target: 100e-3,
                weight: 1.0,
            },
            &sys,
            ctx(),
        );
        assert!(r.is_finite() && r > 0.0, "{r}");
    }

    #[test]
    fn optimizer_moves_the_efl_to_its_target() {
        let catalog = MaterialDatabase::presets();
        let target = 100e-3;
        let ops = [Operand::Efl {
            target,
            weight: 1.0,
        }];
        let result = optimize_system(
            &|x: &[f64]| build(&catalog, x),
            &[120e-3, -120e-3],
            &ops,
            ctx(),
            80,
        );
        let sys = build(&catalog, &result.x).unwrap();
        let efl = sys.paraxial(ctx()).effective_focal_length().unwrap();
        assert!((efl - target).abs() < 1e-3 * target, "EFL {efl}");
    }
}

#[cfg(all(test, feature = "serde"))]
mod serde_tests {
    use super::*;

    #[test]
    fn operands_serialise_with_their_kind() {
        let op = Operand::RayHeight {
            surface: 2,
            field_deg: 1.5,
            pupil_height: 5e-3,
            target: 4e-3,
            axis: RayAxis::Y,
            weight: 2.0,
        };
        let text = serde_json::to_string(&op).unwrap();
        assert!(
            text.contains("RayHeight") && text.contains("\"surface\":2"),
            "{text}"
        );
        let back: Operand = serde_json::from_str(&text).unwrap();
        assert_eq!(back, op);
    }
}
