//! Exact gradients of operands, and an optimiser that uses them.
//!
//! Operands that are read from a traced ray (`RayHeight`) have exact
//! derivatives with respect to the surface curvatures and gaps. Those come
//! from one dual-number trace per parameter, through
//! [`tpt_optics_geometric::differentiable::gradient_of`]. The other operands
//! depend on paraxial rays or on several rays, so their derivatives are
//! central differences in parameter space. The chain from the optimiser's
//! variables to the parameters is a central difference of the builder, which
//! is cheap.
//!
//! The Jacobian comes out as the derivative of each residual with respect to
//! each variable, so [`optimize_exact`] uses Levenberg–Marquardt with the
//! analytic structure rather than a full set of finite differences in the
//! variables.

use alloc::vec::Vec;

use tpt_optics_core::dual::Dual;
use tpt_optics_geometric::differentiable::gradient_of;
use tpt_optics_geometric::sequential::{field_ray, ParamId, SequentialSystem, TraceContext};

use crate::operands::{residual, residuals, Operand};
use crate::optimize::{damped_least_squares_with_jacobian, Control, OptimizeResult};

/// Step for central differences in parameter space (curvature in 1/m,
/// thickness in m).
const PARAM_STEP: f64 = 1e-7;

/// Every curvature and gap of a system, as parameter identifiers and values.
#[must_use]
pub fn parameter_space(sys: &SequentialSystem) -> (Vec<ParamId>, Vec<f64>) {
    let n = sys.surfaces.len();
    let mut ids = Vec::with_capacity(2 * n);
    let mut values = Vec::with_capacity(2 * n);
    for (i, s) in sys.surfaces.iter().enumerate() {
        ids.push(ParamId::Curvature(i));
        values.push(s.asphere.curvature);
    }
    for i in 0..n.saturating_sub(1) {
        ids.push(ParamId::Thickness(i));
        values.push(sys.thickness[i]);
    }
    (ids, values)
}

/// A copy of `sys` with parameter `id` moved by `delta`.
fn shifted(sys: &SequentialSystem, id: ParamId, delta: f64) -> SequentialSystem {
    let mut s = sys.clone();
    match id {
        ParamId::Curvature(i) => s.surfaces[i].asphere.curvature += delta,
        ParamId::Thickness(i) => s.thickness[i] += delta,
        // Conic and asphere terms are not variables of this optimiser.
        _ => {}
    }
    s
}

/// Derivative of one residual with respect to each parameter in `ids`.
#[must_use]
pub fn operand_gradient(
    op: &Operand,
    sys: &SequentialSystem,
    ctx: TraceContext,
    ids: &[ParamId],
) -> Vec<f64> {
    if let Operand::RayHeight {
        surface,
        field_deg,
        pupil_height,
        weight,
        ..
    } = op
    {
        let ray = field_ray(*pupil_height, *field_deg);
        let trace = sys.trace(&ray, ctx);
        if !trace.blocked && trace.points.get(*surface).is_some() {
            let grads = gradient_of(sys, &ray, ctx, ids, |res| {
                res.points.get(*surface).map_or(Dual::new(0.0), |p| p.x)
            });
            return grads.into_iter().map(|g| weight * g).collect();
        }
        return alloc::vec![0.0; ids.len()];
    }
    ids.iter()
        .map(|&id| {
            let h = PARAM_STEP;
            let plus = residual(op, &shifted(sys, id, h), ctx);
            let minus = residual(op, &shifted(sys, id, -h), ctx);
            (plus - minus) / (2.0 * h)
        })
        .collect()
}

/// Jacobian of the residuals with respect to the variables `x`, row-major
/// `ops.len() × x.len()`. `build` maps the variables to a system.
pub fn system_jacobian<B>(build: &B, x: &[f64], ops: &[Operand], ctx: TraceContext) -> Vec<f64>
where
    B: Fn(&[f64]) -> Option<SequentialSystem>,
{
    let n = x.len();
    let m = ops.len();
    let mut out = alloc::vec![0.0; m * n];
    let Some(sys) = build(x) else {
        return out;
    };
    let (ids, _) = parameter_space(&sys);
    // dr/dp: m × np.
    let dr_dp: Vec<Vec<f64>> = ops
        .iter()
        .map(|op| operand_gradient(op, &sys, ctx, &ids))
        .collect();
    // dp/dx: np × n, by central differences of the (cheap) builder.
    let np = ids.len();
    let mut dp_dx = alloc::vec![0.0; np * n];
    let mut xs = x.to_vec();
    for j in 0..n {
        let h = PARAM_STEP.max(1e-9 * libm::fabs(x[j]));
        xs[j] = x[j] + h;
        let plus = build(&xs).map(|s| parameter_space(&s).1);
        xs[j] = x[j] - h;
        let minus = build(&xs).map(|s| parameter_space(&s).1);
        xs[j] = x[j];
        if let (Some(p), Some(q)) = (plus, minus) {
            for k in 0..np {
                dp_dx[k * n + j] = (p[k] - q[k]) / (2.0 * h);
            }
        }
    }
    for i in 0..m {
        for j in 0..n {
            out[i * n + j] = (0..np).map(|k| dr_dp[i][k] * dp_dx[k * n + j]).sum();
        }
    }
    out
}

/// First-order tolerance budget: the standard deviation of each residual
/// when the variables `x` are independently perturbed with standard
/// deviations `sigmas`. It is the square root of the sum over variables of
/// (dr/dx_j times sigma_j) squared, from one Jacobian, with no Monte-Carlo
/// sampling. It is accurate while the errors are small enough for the
/// response to be linear.
pub fn tolerance_budget<B>(
    build: &B,
    x: &[f64],
    sigmas: &[f64],
    ops: &[Operand],
    ctx: TraceContext,
) -> Vec<f64>
where
    B: Fn(&[f64]) -> Option<SequentialSystem>,
{
    let n = x.len();
    let jac = system_jacobian(build, x, ops, ctx);
    (0..ops.len())
        .map(|i| {
            libm::sqrt(
                (0..n)
                    .map(|j| {
                        let s = jac[i * n + j] * sigmas.get(j).copied().unwrap_or(0.0);
                        s * s
                    })
                    .sum::<f64>(),
            )
        })
        .collect()
}

/// Levenberg–Marquardt on the operands, using [`system_jacobian`] for the
/// derivatives. Returns the same kind of result as the other optimisers.
pub fn optimize_exact<B>(
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
            None => out.extend(ops.iter().map(|_| 1.0)),
        }
    };
    let jacobian = |x: &[f64]| system_jacobian(build, x, ops, ctx);
    damped_least_squares_with_jacobian(&r, &jacobian, x0, ops.len(), None, iterations, &mut |_| {
        Control::Continue
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::optimize::numeric_jacobian;
    use tpt_optics_core::material::MaterialDatabase;
    use tpt_optics_core::{Kelvin, Wavelength};
    use tpt_optics_geometric::presets::singlet;

    fn ctx() -> TraceContext {
        TraceContext {
            wavelength: Wavelength::from_nm(587.6),
            temperature: Kelvin::from_celsius(20.0),
        }
    }

    fn build(x: &[f64]) -> Option<SequentialSystem> {
        let catalog = MaterialDatabase::presets();
        singlet(&catalog, x[0], x[1], 12.5e-3, 5e-3, "bk7", 95e-3).ok()
    }

    fn operands() -> Vec<Operand> {
        alloc::vec![
            Operand::Efl {
                target: 100e-3,
                weight: 1.0
            },
            Operand::RayHeight {
                surface: 1,
                field_deg: 0.0,
                pupil_height: 5e-3,
                target: 4.8e-3,
                weight: 1.0,
            },
        ]
    }

    #[test]
    fn ray_height_gradient_is_exact_against_central_differences() {
        // Build the operand residuals by hand, using the same geometry.
        let sys = build(&[120e-3, -120e-3]).unwrap();
        let (ids, _) = parameter_space(&sys);
        let op = operands()[1].clone();
        let exact = operand_gradient(&op, &sys, ctx(), &ids);
        for (k, &id) in ids.iter().enumerate() {
            let h = 1e-7;
            let fd = (residual(&op, &shifted(&sys, id, h), ctx())
                - residual(&op, &shifted(&sys, id, -h), ctx()))
                / (2.0 * h);
            let scale = fd.abs().max(1e-9);
            assert!(
                (exact[k] - fd).abs() < 1e-4 * scale + 1e-12,
                "param {k}: exact {} vs fd {fd}",
                exact[k]
            );
        }
    }

    #[test]
    fn system_jacobian_matches_numeric_jacobian_of_the_residuals() {
        let ops = operands();
        let x = [120e-3, -120e-3];
        let exact = system_jacobian(&build, &x, &ops, ctx());
        let r = |xv: &[f64], out: &mut Vec<f64>| {
            out.clear();
            let sys = build(xv).unwrap();
            out.extend(residuals(&ops, &sys, ctx()));
        };
        let numeric = numeric_jacobian(&r, &x, ops.len(), 1e-7);
        for (a, b) in exact.iter().zip(&numeric) {
            let scale = b.abs().max(1.0);
            assert!((a - b).abs() < 1e-3 * scale, "{exact:?} vs {numeric:?}");
        }
    }

    #[test]
    fn exact_optimiser_matches_the_numerical_one() {
        let ops = operands();
        let x0 = [120e-3, -120e-3];
        let exact = optimize_exact(&build, &x0, &ops, ctx(), 60);
        let numeric = crate::optimize::damped_least_squares(
            &|xv: &[f64], out: &mut Vec<f64>| {
                out.clear();
                let sys = build(xv).unwrap();
                out.extend(residuals(&ops, &sys, ctx()));
            },
            &x0,
            ops.len(),
            None,
            60,
        );
        let sys = build(&exact.x).unwrap();
        let efl = sys.paraxial(ctx()).effective_focal_length().unwrap();
        assert!((efl - 100e-3).abs() < 1e-3 * 100e-3, "EFL {efl}");
        assert!(
            exact.f <= 2.0 * numeric.f + 1e-12,
            "exact {} vs numeric {}",
            exact.f,
            numeric.f
        );
    }
}

#[cfg(test)]
mod budget_tests {
    use super::*;
    use tpt_optics_core::material::MaterialDatabase;
    use tpt_optics_core::{Kelvin, Rng, Wavelength};
    use tpt_optics_geometric::presets::singlet;

    fn ctx() -> TraceContext {
        TraceContext {
            wavelength: Wavelength::from_nm(587.6),
            temperature: Kelvin::from_celsius(20.0),
        }
    }

    fn build(x: &[f64]) -> Option<SequentialSystem> {
        let catalog = MaterialDatabase::presets();
        singlet(&catalog, x[0], x[1], 12.5e-3, 5e-3, "bk7", 95e-3).ok()
    }

    #[test]
    fn linear_budget_agrees_with_monte_carlo_for_small_errors() {
        let ops = alloc::vec![Operand::Efl {
            target: 100e-3,
            weight: 1.0
        }];
        let x = [120e-3, -120e-3];
        let sigmas = [0.05e-3, 0.05e-3];
        let linear = tolerance_budget(&build, &x, &sigmas, &ops, ctx())[0];

        let mut rng = Rng::new(3);
        let mut samples = Vec::new();
        for _ in 0..400 {
            let xs = [
                x[0] + rng.normal(0.0, sigmas[0]),
                x[1] + rng.normal(0.0, sigmas[1]),
            ];
            let sys = build(&xs).unwrap();
            samples.push(sys.paraxial(ctx()).effective_focal_length().unwrap());
        }
        let mean = samples.iter().sum::<f64>() / samples.len() as f64;
        let var =
            samples.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / samples.len() as f64;
        let mc = libm::sqrt(var);
        assert!(
            (linear - mc).abs() < 0.15 * mc,
            "linear {linear} vs Monte Carlo {mc}"
        );
    }
}
