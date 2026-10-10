//! Optical invariants as checkable predicates.
//!
//! Each function returns `Ok(())` when the invariant holds within the given
//! tolerance and `Err(message)` otherwise, so property tests, debug
//! assertions and CI gates share one definition of "the law holds".

use alloc::format;

use tpt_optics_core::{Vec3, Wavelength};
use tpt_optics_geometric::prelude::Ray;
use tpt_optics_geometric::sequential::{SequentialSystem, TraceContext};
use tpt_optics_geometric::snell;

/// Snell's law is reversible: refracting down and back returns −d.
pub fn snell_reversible(
    d: Vec3<f64>,
    n: Vec3<f64>,
    n1: f64,
    n2: f64,
    tol: f64,
) -> Result<(), alloc::string::String> {
    let Some(t1) = snell::refract(d, n, n1, n2) else {
        return Ok(()); // TIR forward — nothing to check
    };
    let Some(t2) = snell::refract(-t1, -n, n2, n1) else {
        return Err(format!(
            "asymmetric TIR: {n1}→{n2} refracts, reverse does not"
        ));
    };
    if (-t2).approx_eq(d, tol) {
        Ok(())
    } else {
        Err(format!("reversed direction {:?} ≠ original {:?}", -t2, d))
    }
}

/// Tangential optical momentum conservation: n₂ sinθ₂ = n₁ sinθ₁.
pub fn snell_tangential_momentum(
    d: Vec3<f64>,
    n: Vec3<f64>,
    n1: f64,
    n2: f64,
    tol: f64,
) -> Result<(), alloc::string::String> {
    let Some(t) = snell::refract(d, n, n1, n2) else {
        return Ok(());
    };
    // In-plane component relative to the tangent direction of d.
    let cos_i = d.dot(n);
    let tangent = (d - n * cos_i).normalized().ok_or("degenerate direction")?;
    let lhs = t.dot(tangent) * n2;
    let rhs = d.dot(tangent) * n1;
    if (lhs - rhs).abs() < tol {
        Ok(())
    } else {
        Err(format!("momentum: n2·s2t = {lhs} vs n1·s1t = {rhs}"))
    }
}

/// ABCD matrix of a unit-index system is symplectic (det = 1).
pub fn abcd_symplectic(
    m: tpt_optics_core::Matrix2<f64>,
    tol: f64,
) -> Result<(), alloc::string::String> {
    let det = m.det();
    if (det - 1.0).abs() < tol {
        Ok(())
    } else {
        Err(format!("det = {det} ≠ 1"))
    }
}

/// Refraction ABCD determinant equals the index ratio.
pub fn abcd_refraction_det(
    n1: f64,
    n2: f64,
    c: f64,
    tol: f64,
) -> Result<(), alloc::string::String> {
    let det = tpt_optics_geometric::abcd::refracting_surface::<f64>(n1, n2, c).det();
    if (det - n1 / n2).abs() < tol {
        Ok(())
    } else {
        Err(format!("det = {det} ≠ n1/n2 = {}", n1 / n2))
    }
}

/// Optical reversibility: tracing the system backwards through the same
/// physical path accumulates the same optical path length.
///
/// The reverse pass uses a mirrored copy of the system (surfaces, media and
/// thicknesses reversed), since a sequential tracer consumes surfaces in
/// propagation order.
pub fn trace_reversibility(
    system: &SequentialSystem,
    ray: &Ray<f64>,
    ctx: TraceContext,
    tol: f64,
) -> Result<(), alloc::string::String> {
    let fwd = system.trace(ray, ctx);
    if fwd.blocked {
        return Ok(());
    }
    let Some(last) = fwd.points.last() else {
        return Ok(());
    };
    // Mirrored system.
    let mut rev = SequentialSystem::new();
    rev.media = {
        let mut m: alloc::vec::Vec<_> = system.media.iter().rev().cloned().collect();
        m.truncate(system.surfaces.len() + 1);
        m
    };
    let n = system.surfaces.len();
    for i in (0..n).rev() {
        // Gap between mirrored surfaces j and j+1 equals the original gap
        // between surfaces i−1 and i (i is the mirrored source).
        let thickness_after_prev = if i > 0 { system.thickness[i - 1] } else { 0.0 };
        // Mirror the surface about its vertex plane: z-flip negates the
        // curvature and polynomial coefficients (conic unchanged).
        let mut mirrored = system.surfaces[i].clone();
        mirrored.asphere.curvature = -mirrored.asphere.curvature;
        for cf in &mut mirrored.asphere.coeffs {
            *cf = -*cf;
        }
        if let Some(phase) = mirrored.phase.as_mut() {
            for c in phase.iter_mut() {
                *c = -*c;
            }
        }
        rev.push(mirrored, thickness_after_prev, system.media[i].clone());
    }
    // Note: `push` appended medium_after again; rebuild media precisely.
    rev.media = {
        let mut m: alloc::vec::Vec<_> = system.media.iter().rev().cloned().collect();
        m.truncate(n + 1);
        m
    };
    // In the mirrored frame (z' = L − z) the reversed ray starts at the
    // image surface (z' = 0) heading +z' with the transverse velocity
    // negated.
    let e = fwd.exit_direction;
    let Some(dir) = tpt_optics_core::Vec3::new(-e.x, -e.y, e.z).normalized() else {
        return Ok(());
    };
    let origin = tpt_optics_core::Vec3::new(last.x, last.y, 0.0) - dir * 1e-6;
    let back_ray = Ray::new(origin, dir);
    let back = rev.trace(&back_ray, ctx);
    if back.blocked {
        return Err("reverse trace blocked where forward traced".into());
    }
    // Equalize spans: the forward OPD includes the lead-in from the ray
    // origin to the first surface; the backward OPD includes the epsilon
    // offset before the mirrored image. Subtract both lead-ins and compare
    // the shared surface-to-surface span.
    let n0 = system.media[0].index_at(ctx.wavelength, ctx.temperature);
    let lead = (fwd.points[0] - ray.origin).norm() * n0;
    let fwd_core = fwd.opd - lead;
    let back_core = back.opd - 1e-6;
    let d = (fwd_core - back_core).abs();
    if d < tol {
        Ok(())
    } else {
        Err(format!(
            "OPD mismatch: forward {fwd_core} vs back {back_core} (span diff {d})"
        ))
    }
}

/// Fresnel energy conservation for lossless media: R + T = 1 per
/// polarization.
pub fn fresnel_energy_conserved(
    n1: f64,
    n2: f64,
    theta_i: f64,
    tol: f64,
) -> Result<(), alloc::string::String> {
    if let Some(crit) = tpt_optics_core::fresnel::critical_angle(n1, n2) {
        if theta_i >= crit {
            return Ok(()); // TIR: R = 1 trivially
        }
    }
    let (rs, rp, ts, tp) = tpt_optics_core::fresnel::fresnel_energy(n1, n2, theta_i);
    for (label, r, t) in [("s", rs, ts), ("p", rp, tp)] {
        if (r + t - 1.0).abs() < tol {
            continue;
        }
        return Err(format!("{label}: R + T = {}", r + t));
    }
    Ok(())
}

/// Optical path length is invariant under rotation about the optical axis
/// for rotationally symmetric (on-axis aspheric) systems.
pub fn opl_rotation_invariance(
    system: &SequentialSystem,
    ray: &Ray<f64>,
    ctx: TraceContext,
    azimuth: f64,
    tol: f64,
) -> Result<(), alloc::string::String> {
    let fwd = system.trace(ray, ctx);
    // Rotate origin and direction about z.
    let (c, s) = (libm::cos(azimuth), libm::sin(azimuth));
    let rot = |v: &tpt_optics_core::Vec3<f64>| {
        tpt_optics_core::Vec3::new(c * v.x - s * v.y, s * v.x + c * v.y, v.z)
    };
    let rotated = Ray::new(rot(&ray.origin), rot(&ray.direction));
    let moved = system.trace(&rotated, ctx);
    if fwd.blocked || moved.blocked {
        return Ok(());
    }
    let d = (fwd.opd - moved.opd).abs();
    if d < tol {
        Ok(())
    } else {
        Err(format!(
            "OPL differs under rotation: {} vs {}",
            fwd.opd, moved.opd
        ))
    }
}

/// Wavelength of the sodium d line for default traces.
#[must_use]
pub fn d_line() -> Wavelength {
    Wavelength::from_nm(587.6)
}
