//! Sequential ray tracing through lens systems.
//!
//! A [`SequentialSystem`] is a list of surfaces along a common optical axis
//! (z), each with an action (refract / reflect / diffract), an aspheric
//! geometry, an aperture semi-diameter, and a medium on each side. Media
//! resolve to refractive indices at the trace wavelength/temperature, so a
//! single trace is a single-wavelength monochromatic event.
//!
//! # Differentiability
//!
//! The trace kernels are generic over [`Scalar`]. Geometry that the
//! optimiser varies lives in a [`TunedParams`] value: build
//! [`TunedParams::<Dual>::seeded`] with one [`ParamId`] and the same trace
//! returns exact derivatives of any traced quantity with respect to that
//! parameter (see [`crate::differentiable`]).

use alloc::sync::Arc;
use alloc::vec;
use alloc::vec::Vec;

use libm::sqrt;
use tpt_optics_core::index::IndexAt;
use tpt_optics_core::material::Material;
use tpt_optics_core::{Dual, Kelvin, Scalar, Vec2, Vec3, Wavelength};

use crate::ray::Ray;
use crate::snell;
use crate::surface::Asphere;

/// Generic asphere sag from scalar parameters (keeps dual derivatives of
/// curvature/conic/polynomial coefficients alive through the trace).
fn asphere_sag<T: Scalar>(c: T, k: T, a: &[T; 8], x: T, y: T) -> T {
    let r2 = x * x + y * y;
    let denom = T::ONE + (T::ONE - (T::ONE + k) * c * c * r2).sqrt();
    let mut z = c * r2 / denom;
    let mut rpow = r2;
    for ai in a {
        z += *ai * rpow;
        rpow *= r2;
    }
    z
}

/// Generic asphere sag gradient.
fn asphere_grad<T: Scalar>(c: T, k: T, a: &[T; 8], x: T, y: T) -> Vec2<T> {
    let r2 = x * x + y * y;
    let denom = (T::ONE - (T::ONE + k) * c * c * r2).sqrt();
    let mut gx = c * x / denom;
    let mut gy = c * y / denom;
    let mut r2pow = T::ONE;
    for (i, ai) in a.iter().enumerate() {
        let coef = T::from_f64(2.0 * (i + 1) as f64) * *ai;
        gx += coef * r2pow * x;
        gy += coef * r2pow * y;
        r2pow *= r2;
    }
    Vec2::new(gx, gy)
}

/// Newton-Raphson intersection of a parametric asphere (quadric seed),
/// fully generic over the scalar type. Returns `(t, point, normal)` with
/// the normal oriented against the ray.
fn intersect_asphere<T: Scalar>(
    c: T,
    k: T,
    a: &[T; 8],
    ray: &Ray<T>,
    t_min: T,
) -> Option<(T, Vec3<T>, Vec3<T>)> {
    let d = ray.direction;
    let o = ray.origin;
    let mut t;
    if c == T::ZERO {
        if d.z == T::ZERO {
            return None;
        }
        t = -o.z / d.z;
        if t <= t_min {
            return None;
        }
    } else {
        // Quadric seed: x^2 + y^2 + (1+k)z^2 - 2z/c = 0.
        let ok = T::ONE + k;
        let aa = d.x * d.x + d.y * d.y + ok * d.z * d.z;
        let inv_c = T::ONE / c;
        let b = T::from_f64(2.0) * (o.x * d.x + o.y * d.y + ok * o.z * d.z - d.z * inv_c);
        let cc = o.x * o.x + o.y * o.y + ok * o.z * o.z - T::from_f64(2.0) * o.z * inv_c;
        let disc = b * b - T::from_f64(4.0) * aa * cc;
        if disc < T::ZERO || aa == T::ZERO {
            return None;
        }
        let sq = disc.sqrt();
        let inv = T::ONE / (T::from_f64(2.0) * aa);
        let t1 = (-b - sq) * inv;
        let t2 = (-b + sq) * inv;
        t = if t1 > t_min {
            t1
        } else if t2 > t_min {
            t2
        } else {
            return None;
        };
    }
    // Newton refinement on F(t) = ray.z(t) - sag(...).
    for _ in 0..16 {
        let pt = ray.at(t);
        let f = pt.z - asphere_sag(c, k, a, pt.x, pt.y);
        let g = asphere_grad(c, k, a, pt.x, pt.y);
        let df = d.z - (g.x * d.x + g.y * d.y);
        if df == T::ZERO {
            break;
        }
        let step = f / df;
        t -= step;
        if !t.is_finite() {
            return None;
        }
        if T::to_f64(step.abs()) < 1e-15 {
            break;
        }
    }
    if t <= t_min {
        return None;
    }
    let point = ray.at(t);
    let residual = point.z - asphere_sag(c, k, a, point.x, point.y);
    if T::to_f64(residual.abs()) > 1e-7 {
        return None;
    }
    let g = asphere_grad(c, k, a, point.x, point.y);
    let up = Vec3::new(-g.x, -g.y, T::ONE)
        .normalized()
        .unwrap_or_else(|| Vec3::new(T::ZERO, T::ZERO, T::ONE));
    let normal = if d.dot(up) > T::ZERO { -up } else { up };
    Some((t, point, normal))
}

/// What happens to a ray hitting a surface.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum SurfaceAction {
    /// Refract into the next medium.
    Refract,
    /// Specular reflection (mirror).
    Reflect,
    /// Refraction plus a diffractive phase kick (HOE / DOE / Fresnel lens).
    Diffract,
}

/// Geometric description of one sequential surface.
#[derive(Clone, Debug, PartialEq)]
pub struct SeqSurface {
    /// Refraction / reflection / diffraction behaviour.
    pub action: SurfaceAction,
    /// Base aspheric geometry (plane when curvature and all terms vanish).
    pub asphere: Asphere,
    /// Aperture semi-diameter (meters). Rays outside are blocked.
    pub semi_diameter: f64,
    /// Diffractive phase polynomial (row-major, side 3: c₀₀ c₀₁ c₀₂ c₁₀ …
    /// c₂₂) — used only when `action == Diffract`.
    pub phase: Option<[f64; 9]>,
}

impl SeqSurface {
    /// A flat refracting surface (e.g. a stop or plain interface).
    pub fn flat(semi_diameter: f64) -> Self {
        SeqSurface {
            action: SurfaceAction::Refract,
            asphere: Asphere::flat(),
            semi_diameter,
            phase: None,
        }
    }

    /// A spherical/conic/aspheric refracting surface.
    pub fn aspheric(asphere: Asphere, semi_diameter: f64) -> Self {
        SeqSurface {
            action: SurfaceAction::Refract,
            asphere,
            semi_diameter,
            phase: None,
        }
    }

    /// A spherical mirror (positive curvature = concave toward +z).
    pub fn mirror(curvature: f64, semi_diameter: f64) -> Self {
        SeqSurface {
            action: SurfaceAction::Reflect,
            asphere: Asphere::sphere(curvature),
            semi_diameter,
            phase: None,
        }
    }

    /// A diffractive surface on a (usually flat) base.
    pub fn diffractive(asphere: Asphere, phase: [f64; 9], semi_diameter: f64) -> Self {
        SeqSurface {
            action: SurfaceAction::Diffract,
            asphere,
            semi_diameter,
            phase: Some(phase),
        }
    }

    fn phase_gradient(&self, x: f64, y: f64) -> Vec2<f64> {
        let mut gx = 0.0;
        let mut gy = 0.0;
        if let Some(p) = &self.phase {
            for i in 0..3usize {
                for j in 0..3usize {
                    let cij = p[i * 3 + j];
                    if cij == 0.0 {
                        continue;
                    }
                    if i > 0 {
                        gx +=
                            cij * i as f64 * libm::pow(x, (i - 1) as f64) * libm::pow(y, j as f64);
                    }
                    if j > 0 {
                        gy +=
                            cij * j as f64 * libm::pow(x, i as f64) * libm::pow(y, (j - 1) as f64);
                    }
                }
            }
        }
        Vec2::new(gx, gy)
    }
}

/// A medium between surfaces: either a constant index or a material model.
#[derive(Clone, Debug)]
pub enum Medium {
    /// Constant real refractive index (air approximations, interstice).
    Const(f64),
    /// Full material model (dispersion + temperature); carries the
    /// [`Material`](tpt_optics_core::material::Material) so names survive
    /// file-format round trips.
    Material(Arc<Material>),
    /// Vacuum.
    Vacuum,
}

impl Medium {
    /// Refractive index at (wavelength, temperature).
    pub fn index_at(&self, wl: Wavelength, temp: Kelvin) -> f64 {
        match self {
            Medium::Const(n) => *n,
            Medium::Vacuum => 1.0,
            Medium::Material(m) => m.index_at(wl, temp).n,
        }
    }
}

/// Monochromatic trace context.
#[derive(Copy, Clone, Debug)]
pub struct TraceContext {
    /// Vacuum wavelength of the trace.
    pub wavelength: Wavelength,
    /// System temperature.
    pub temperature: Kelvin,
}

impl Default for TraceContext {
    fn default() -> Self {
        TraceContext {
            wavelength: Wavelength::from_nm(587.6),
            temperature: Kelvin::from_celsius(20.0),
        }
    }
}

/// Variable geometry handed to the generic tracer.
///
/// The `f64` instantiation is the nominal system; a `Dual` instantiation
/// seeds one parameter for automatic differentiation.
#[derive(Clone, Debug)]
pub struct TunedParams<T: Scalar> {
    /// Curvature (1/m) per surface.
    pub curvature: Vec<T>,
    /// Conic constant per surface.
    pub conic: Vec<T>,
    /// Even asphere coefficients per surface.
    pub asphere: Vec<[T; 8]>,
    /// Axial separation from surface i to surface i+1 (meters).
    pub thickness: Vec<T>,
}

impl TunedParams<f64> {
    /// Nominal parameters of a system.
    pub fn from_system(sys: &SequentialSystem) -> Self {
        TunedParams {
            curvature: sys.surfaces.iter().map(|s| s.asphere.curvature).collect(),
            conic: sys.surfaces.iter().map(|s| s.asphere.conic).collect(),
            asphere: sys.surfaces.iter().map(|s| s.asphere.coeffs).collect(),
            thickness: sys.thickness.clone(),
        }
    }
}

impl TunedParams<Dual> {
    /// Parameters with exactly one variable seeded (derivative 1).
    pub fn seeded(sys: &SequentialSystem, id: ParamId) -> Self {
        let mut p = TunedParams {
            curvature: sys
                .surfaces
                .iter()
                .map(|s| Dual::new(s.asphere.curvature))
                .collect(),
            conic: sys
                .surfaces
                .iter()
                .map(|s| Dual::new(s.asphere.conic))
                .collect(),
            asphere: sys
                .surfaces
                .iter()
                .map(|s| {
                    let mut a = [Dual::new(0.0); 8];
                    for (i, v) in s.asphere.coeffs.iter().enumerate() {
                        a[i] = Dual::new(*v);
                    }
                    a
                })
                .collect(),
            thickness: sys.thickness.iter().map(|t| Dual::new(*t)).collect(),
        };
        match id {
            ParamId::Curvature(i) => p.curvature[i] = Dual::var(sys.surfaces[i].asphere.curvature),
            ParamId::Conic(i) => p.conic[i] = Dual::var(sys.surfaces[i].asphere.conic),
            ParamId::AsphereCoeff(i, j) => {
                p.asphere[i][j] = Dual::var(sys.surfaces[i].asphere.coeffs[j]);
            }
            ParamId::Thickness(i) => p.thickness[i] = Dual::var(sys.thickness[i]),
        }
        p
    }
}

/// A design variable of the sequential system.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum ParamId {
    /// Curvature of surface `i`.
    Curvature(usize),
    /// Conic constant of surface `i`.
    Conic(usize),
    /// Even asphere coefficient `j` of surface `i`.
    AsphereCoeff(usize, usize),
    /// Axial separation after surface `i`.
    Thickness(usize),
}

/// Outcome of a sequential trace.
#[derive(Clone, Debug)]
pub struct TraceResult<T: Scalar> {
    /// Intersection point on each surface (blocked traces end early).
    pub points: Vec<Vec3<T>>,
    /// Optical path length from the start point to the final surface,
    /// accumulated as n·|Δs| in each medium (meters).
    pub opd: T,
    /// Direction after the final surface.
    pub exit_direction: Vec3<T>,
    /// True when the ray was blocked by an aperture before the last surface.
    pub blocked: bool,
}

/// A sequential optical system.
#[derive(Clone, Debug)]
pub struct SequentialSystem {
    /// Surfaces in propagation order; the last one is the image surface.
    pub surfaces: Vec<SeqSurface>,
    /// Axial separation from surface i to i+1 (meters); the final entry is
    /// unused (kept so indices line up with `surfaces`).
    pub thickness: Vec<f64>,
    /// Media: `media[i]` is the medium between surface i−1 and surface i.
    /// Length must be `surfaces.len() + 1` (index 0 = object space).
    pub media: Vec<Medium>,
}

impl SequentialSystem {
    /// An empty system (object and image space both vacuum).
    pub fn new() -> Self {
        SequentialSystem {
            surfaces: Vec::new(),
            thickness: Vec::new(),
            // media[i] = medium before surface i; starts with object space.
            media: vec![Medium::Vacuum],
        }
    }

    /// Append a surface.
    pub fn push(&mut self, surface: SeqSurface, thickness: f64, medium_after: Medium) {
        self.surfaces.push(surface);
        self.thickness.push(thickness);
        self.media.push(medium_after);
    }

    /// Number of surfaces.
    pub fn len(&self) -> usize {
        self.surfaces.len()
    }

    /// True when the system has no surfaces.
    pub fn is_empty(&self) -> bool {
        self.surfaces.is_empty()
    }

    /// Refractive index before surface `i` at the context.
    pub fn index_before(&self, i: usize, ctx: TraceContext) -> f64 {
        self.media[i].index_at(ctx.wavelength, ctx.temperature)
    }

    /// Refractive index after surface `i`.
    pub fn index_after(&self, i: usize, ctx: TraceContext) -> f64 {
        self.media[i + 1].index_at(ctx.wavelength, ctx.temperature)
    }

    /// Vertex z of surface `i` under the given thickness parameters.
    fn vertex_z<T: Scalar>(thickness: &[T], i: usize) -> T {
        let mut z = T::ZERO;
        for t in thickness.iter().take(i) {
            z += *t;
        }
        z
    }

    /// Monochromatic exact trace, generic over the scalar type.
    ///
    /// `ray` is expressed in the system frame with surface 0's vertex at the
    /// origin. Rays starting before the first surface propagate to it; the
    /// trace ends at the last surface (the image surface).
    pub fn trace_with<T: Scalar>(
        &self,
        params: &TunedParams<T>,
        ray: &Ray<T>,
        ctx: TraceContext,
    ) -> TraceResult<T> {
        let n_surfaces = self.surfaces.len();
        let mut points = Vec::with_capacity(n_surfaces);
        let mut opd = T::ZERO;
        // Refractive index of the current medium (kept in f64: material
        // derivatives wrt geometry parameters are not modelled).
        let mut n_cur = self.index_before(0, ctx);
        let mut current = Ray::new(
            ray.origin,
            ray.direction.normalized().unwrap_or(ray.direction),
        );
        let mut blocked = false;
        let mut exit_dir = current.direction;

        for i in 0..n_surfaces {
            let surf = &self.surfaces[i];
            let z_i = Self::vertex_z(&params.thickness, i);
            // Ray in the surface frame.
            let local = Ray::new(
                current.origin - Vec3::new(T::ZERO, T::ZERO, z_i),
                current.direction,
            );
            // Fully generic intersection keeps Dual derivatives alive.
            let Some((hit_t, hit_point, hit_normal)) = intersect_asphere(
                params.curvature[i],
                params.conic[i],
                &params.asphere[i],
                &local,
                T::from_f64(1e-9),
            ) else {
                blocked = true;
                break;
            };
            let point_global = hit_point + Vec3::new(T::ZERO, T::ZERO, z_i);
            // Aperture check (branch on values — stable for AD).
            let r2 = T::to_f64(hit_point.x * hit_point.x + hit_point.y * hit_point.y);
            if r2 > surf.semi_diameter * surf.semi_diameter {
                blocked = true;
                break;
            }
            // Optical path: n·|t·d| (direction is unit).
            opd += T::from_f64(n_cur) * hit_t;

            let n_next = self.index_after(i, ctx);
            let d = current.direction;
            let new_dir = match surf.action {
                SurfaceAction::Refract => {
                    // TIR mirrors the ray and keeps the current medium.
                    if let Some(t) = snell::refract(d, hit_normal, n_cur, n_next) {
                        t
                    } else {
                        let mirrored = snell::reflect(d, hit_normal);
                        points.push(point_global);
                        exit_dir = mirrored;
                        current = Ray::new(point_global, mirrored);
                        continue;
                    }
                }
                SurfaceAction::Reflect => {
                    let mirrored = snell::reflect(d, hit_normal);
                    points.push(point_global);
                    exit_dir = mirrored;
                    current = Ray::new(point_global, mirrored);
                    continue;
                }
                SurfaceAction::Diffract => {
                    // Refraction + phase kick of the polynomial in surface
                    // coordinates. The phase gradient is evaluated on the
                    // nominal geometry (phase coefficients are not variables
                    // for the AD path).
                    let xg = T::to_f64(hit_point.x);
                    let yg = T::to_f64(hit_point.y);
                    let grad = surf.phase_gradient(xg, yg);
                    let refracted = snell::refract(d, hit_normal, n_cur, n_next)
                        .unwrap_or_else(|| snell::reflect(d, hit_normal));
                    snell::diffract(
                        refracted,
                        hit_normal,
                        n_next,
                        n_next,
                        ctx.wavelength.m(),
                        Vec2::new(T::from_f64(grad.x), T::from_f64(grad.y)),
                    )
                    .unwrap_or(refracted)
                }
            };
            n_cur = n_next;
            points.push(point_global);
            exit_dir = new_dir;
            current = Ray::new(point_global, new_dir);
        }

        TraceResult {
            points,
            opd,
            exit_direction: exit_dir,
            blocked,
        }
    }

    /// Trace with the nominal (f64) parameters.
    pub fn trace(&self, ray: &Ray<f64>, ctx: TraceContext) -> TraceResult<f64> {
        let params = TunedParams::from_system(self);
        self.trace_with(&params, ray, ctx)
    }

    /// Paraxial (first-order) description at the context wavelength.
    pub fn paraxial(&self, ctx: TraceContext) -> ParaxialSystem {
        let mut surfaces = Vec::with_capacity(self.surfaces.len());
        for i in 0..self.surfaces.len() {
            surfaces.push(ParaxialSurface {
                curvature: self.surfaces[i].asphere.curvature,
                thickness: self.thickness[i],
                n_before: self.index_before(i, ctx),
                n_after: self.index_after(i, ctx),
                reflective: self.surfaces[i].action == SurfaceAction::Reflect,
                semi_diameter: self.surfaces[i].semi_diameter,
            });
        }
        ParaxialSystem { surfaces }
    }
}

impl Default for SequentialSystem {
    fn default() -> Self {
        Self::new()
    }
}

/// One surface in the paraxial (first-order) model.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct ParaxialSurface {
    /// Curvature (1/m).
    pub curvature: f64,
    /// Separation to the next surface (m).
    pub thickness: f64,
    /// Index before.
    pub n_before: f64,
    /// Index after.
    pub n_after: f64,
    /// Mirror flag.
    pub reflective: bool,
    /// Aperture semi-diameter (m).
    pub semi_diameter: f64,
}

/// First-order (paraxial) model of a system: `y`-`u` ray tracing, ABCD
/// matrix, effective and back focal lengths, and entrance pupil.
#[derive(Clone, Debug, Default)]
pub struct ParaxialSystem {
    /// Surfaces in order.
    pub surfaces: Vec<ParaxialSurface>,
}

impl ParaxialSystem {
    /// Refraction/reflection transfer of (y, u) at surface i.
    fn refract_pair(&self, i: usize, y: f64, u: f64) -> (f64, f64) {
        let s = &self.surfaces[i];
        if s.reflective {
            // Mirror: u' = u − 2 c y (index sign convention folded in);
            // direction sense flips but paraxial formalism keeps going.
            (y, u - 2.0 * s.curvature * y)
        } else {
            let phi = (s.n_after - s.n_before) * s.curvature;
            let u2 = (s.n_before * u - y * phi) / s.n_after;
            (y, u2)
        }
    }

    /// Trace a paraxial ray `(y, u)` from before the first surface to the
    /// last surface; returns the state at every surface.
    pub fn trace_ray(&self, y0: f64, u0: f64) -> Vec<(f64, f64)> {
        let mut states = Vec::with_capacity(self.surfaces.len());
        let mut y = y0;
        let mut u = u0;
        for i in 0..self.surfaces.len() {
            let (y2, u2) = self.refract_pair(i, y, u);
            states.push((y2, u2));
            // Transfer to next surface.
            if i + 1 < self.surfaces.len() {
                y = y2 + u2 * self.surfaces[i].thickness;
                u = u2;
            }
        }
        states
    }

    /// Effective focal length for an object at infinity in the object-space
    /// medium (EFL = −y₀/u_exit, y₀ = 1). Returns `None` for afocal systems.
    pub fn effective_focal_length(&self) -> Option<f64> {
        let states = self.trace_ray(1.0, 0.0);
        let (_, u_exit) = states.last()?;
        if u_exit.abs() < 1e-14 {
            return None;
        }
        Some(-1.0 / u_exit)
    }

    /// Back focal distance from the last surface vertex to the axial focus.
    pub fn back_focal_length(&self) -> Option<f64> {
        let states = self.trace_ray(1.0, 0.0);
        let (y_last, u_exit) = states.last()?;
        if u_exit.abs() < 1e-14 {
            return None;
        }
        Some(-y_last / u_exit)
    }

    /// ABCD matrix from just before the first surface to just after the
    /// last (composed element matrices; det = n₀/n_image).
    pub fn abcd(&self) -> tpt_optics_core::Matrix2<f64> {
        use tpt_optics_core::Matrix2;
        let mut m = Matrix2::identity();
        for i in 0..self.surfaces.len() {
            let s = &self.surfaces[i];
            // Element matrix at surface i.
            let e = if s.reflective {
                crate::abcd::mirror(s.curvature)
            } else {
                crate::abcd::refracting_surface(s.n_before, s.n_after, s.curvature)
            };
            m = e.mul_matrix(m);
            // Translation (unused after the final surface).
            if i + 1 < self.surfaces.len() {
                m = crate::abcd::free_space(s.thickness).mul_matrix(m);
            }
        }
        m
    }

    /// Entrance pupil semi-diameter given a stop (aperture-limiting)
    /// surface index, and the paraxial F/# using the EFL.
    ///
    /// Traces a unit-height parallel ray and scales it so it just clears
    /// the stop.
    pub fn entrance_pupil(&self, stop_index: usize) -> Option<(f64, f64)> {
        let states = self.trace_ray(1.0, 0.0);
        let (y_stop, _) = states.get(stop_index)?;
        if y_stop.abs() < 1e-14 {
            return None;
        }
        let epd = self.surfaces[stop_index].semi_diameter / y_stop.abs();
        // `epd` is a semi-diameter, so the diameter (and the F/#
        // denominator) is 2·epd.
        let fnum = self.effective_focal_length().map(|f| f / (2.0 * epd));
        Some((epd, fnum.unwrap_or(f64::NAN)))
    }
}

/// Ray at field angle `field_angle_deg` (degrees) that crosses the first
/// surface plane (`z = 0`) at transverse position `(x, 0)`.
///
/// It starts 100 mm before that plane, back along the tilted direction, so a
/// field tilt does not shift the bundle off the pupil point.
#[must_use]
pub fn field_ray(x: f64, field_angle_deg: f64) -> Ray<f64> {
    let theta = field_angle_deg.to_radians();
    let (s, c) = (libm::sin(theta), libm::cos(theta));
    // Back along `dir` to z = −0.1 is a distance 0.1/c; the transverse
    // offset over that distance is −0.1·s/c.
    Ray::new(Vec3::new(x, -0.1 * s / c, -0.1), Vec3::new(0.0, s, c))
}

/// RMS transverse spot radius at the image surface for a fan of rays,
/// measured about the spot centroid.
///
/// `heights` are pupil offsets in meters; each ray is launched by
/// [`field_ray`] at the given field angle (degrees).
pub fn rms_spot_radius(
    system: &SequentialSystem,
    ctx: TraceContext,
    field_angle_deg: f64,
    heights: &[f64],
) -> f64 {
    let mut spot: Vec<(f64, f64)> = Vec::with_capacity(heights.len());
    for &h in heights {
        let res = system.trace(&field_ray(h, field_angle_deg), ctx);
        if res.blocked {
            continue;
        }
        if let Some(p) = res.points.last() {
            spot.push((p.x, p.y));
        }
    }
    if spot.is_empty() {
        return f64::NAN;
    }
    let n = spot.len() as f64;
    let cx = spot.iter().map(|p| p.0).sum::<f64>() / n;
    let cy = spot.iter().map(|p| p.1).sum::<f64>() / n;
    let mean_sq = spot
        .iter()
        .map(|p| (p.0 - cx) * (p.0 - cx) + (p.1 - cy) * (p.1 - cy))
        .sum::<f64>()
        / n;
    sqrt(mean_sq)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::surface::{OpticalSurface, Sphere};

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    /// Symmetric biconvex BK7-like singlet, R = 100 mm, t = 5 mm.
    fn singlet() -> SequentialSystem {
        let mut sys = SequentialSystem::new();
        // Object space air; surfaces: R1 +100mm, R2 −100mm.
        sys.media[0] = Medium::Const(1.0);
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
        // Image plane.
        sys.push(SeqSurface::flat(12e-3), 0.0, Medium::Vacuum);
        sys
    }

    #[test]
    fn thick_lens_paraxial_efl_matches_formula() {
        let sys = singlet();
        let px = sys.paraxial(TraceContext::default());
        let efl = px.effective_focal_length().unwrap();
        // 1/f = (n−1)[c1 − c2] + (n−1)²/n·t·c1c2
        let n = 1.5168;
        let c1 = 10.0;
        let c2 = -10.0;
        let t = 5e-3;
        let inv_f = (n - 1.0) * (c1 - c2) + (n - 1.0) * (n - 1.0) / n * t * c1 * c2;
        let expect = 1.0 / inv_f;
        assert!(close(efl, expect, 1e-6), "efl {efl} vs {expect}");
    }

    #[test]
    fn exact_trace_near_paraxial_focus() {
        let sys = singlet();
        let px = sys.paraxial(TraceContext::default());
        let bfl = px.back_focal_length().unwrap();
        // Paraxial-height ray: hits image plane near axis when image plane
        // is at bfl. Our system's image surface is at vertex 5mm + 0:
        // instead evaluate the crossing after the last lens surface.
        let ray = Ray::new(Vec3::new(1e-3, 0.0, -0.01), Vec3::new(0.0, 0.0, 1.0));
        let res = sys.trace(&ray, TraceContext::default());
        assert!(!res.blocked);
        // Exit direction must converge toward the axis.
        assert!(res.exit_direction.x < 0.0, "converging expected");
        // Absolute focus position = second vertex + bfl (bfl is measured
        // from the last surface).
        let x2 = res.points[1].x;
        let z2 = 5e-3; // second surface vertex
        let cross = z2 - x2 / res.exit_direction.x * res.exit_direction.z;
        assert!(
            close(cross, z2 + bfl, 2e-4),
            "cross {cross} vs z2+bfl {bfl}"
        );
    }

    #[test]
    fn mirror_folds_and_reflects() {
        let mut sys = SequentialSystem::new();
        sys.push(SeqSurface::mirror(1.0 / 0.5, 20e-3), 0.0, Medium::Vacuum);
        // Axial ray reflects back on itself.
        let ray = Ray::new(Vec3::new(0.0, 0.0, -0.2), Vec3::new(0.0, 0.0, 1.0));
        let res = sys.trace(&ray, TraceContext::default());
        assert!(close(res.exit_direction.z, -1.0, 1e-12));
        // OPD = n·0.2
        assert!(close(res.opd, 0.2, 1e-12));
    }

    #[test]
    fn aperture_blocks() {
        let mut sys = SequentialSystem::new();
        sys.push(SeqSurface::flat(1e-3), 0.0, Medium::Vacuum);
        sys.push(SeqSurface::flat(1e-3), 0.0, Medium::Vacuum);
        let ray = Ray::new(Vec3::new(5e-3, 0.0, -0.1), Vec3::new(0.0, 0.0, 1.0));
        let res = sys.trace(&ray, TraceContext::default());
        assert!(res.blocked);
    }

    #[test]
    fn diffractive_lens_focuses() {
        // DOE thin-lens-equivalent: phase −π r²/(λ f), f = 25 mm.
        let f = 25e-3;
        let wl = Wavelength::from_nm(550.0);
        let phase = {
            let c = -core::f64::consts::PI / (wl.m() * f);
            // Row-major [i*3 + j] multiplies x^i*y^j: the x^2 term is slot 6.
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, c, 0.0, 0.0]
        };
        let mut sys = SequentialSystem::new();
        sys.push(
            SeqSurface::diffractive(Asphere::flat(), phase, 5e-3),
            60e-3,
            Medium::Const(1.0),
        );
        sys.push(SeqSurface::flat(5e-3), 0.0, Medium::Vacuum);
        // Collimated rays at several heights focus at ≈ f.
        for h in [1e-3, 2e-3, 3e-3] {
            let ray = Ray::new(Vec3::new(h, 0.0, -0.01), Vec3::new(0.0, 0.0, 1.0));
            let res = sys.trace(&ray, wl_ctx(wl));
            assert!(!res.blocked);
            let x2 = res.points[0].x;
            let cross = -x2 / res.exit_direction.x * res.exit_direction.z;
            assert!((cross - f).abs() < 0.02 * f, "h={h} cross {cross} vs f {f}");
        }
    }

    fn wl_ctx(wl: Wavelength) -> TraceContext {
        TraceContext {
            wavelength: wl,
            temperature: Kelvin::from_celsius(20.0),
        }
    }

    #[test]
    fn opd_accumulates_index_times_length() {
        // Two flat surfaces, 10 mm of glass n = 1.5.
        let mut sys = SequentialSystem::new();
        sys.media[0] = Medium::Const(1.0);
        sys.push(SeqSurface::flat(10e-3), 10e-3, Medium::Const(1.5));
        sys.push(SeqSurface::flat(10e-3), 0.0, Medium::Vacuum);
        let ray = Ray::new(Vec3::new(0.0, 0.0, -5e-3), Vec3::new(0.0, 0.0, 1.0));
        let res = sys.trace(&ray, TraceContext::default());
        // 5 mm air + 10 mm glass
        assert!(close(res.opd, 5e-3 + 10e-3 * 1.5, 1e-12));
    }

    #[test]
    fn sphere_surface_direct_use() {
        // Sanity: Sphere trait object intersection agrees with Asphere.
        let s = Sphere::from_radius(50e-3);
        let a = Asphere::sphere(1.0 / 50e-3);
        let ray = Ray::new(Vec3::new(5e-3, 0.0, -10e-3), Vec3::new(0.0, 0.0, 1.0));
        let h1 = s.intersect(&ray, 1e-12).unwrap();
        let h2 = a.intersect(&ray, 1e-12).unwrap();
        assert!(h1.point.approx_eq(h2.point, 1e-12));
    }

    #[test]
    fn rms_spot_on_axis_is_tiny() {
        let sys = singlet();
        // Move the image plane to the paraxial focus for a small spot.
        let mut sys = sys;
        let px = sys.paraxial(TraceContext::default());
        let bfl = px.back_focal_length().unwrap();
        sys.thickness[1] = bfl;
        let rms = rms_spot_radius(&sys, TraceContext::default(), 0.0, &[0.0, 0.3e-3, 0.6e-3]);
        assert!(rms < 200e-6, "rms = {rms}");
    }
}

#[cfg(test)]
mod seed_and_fnum_tests {
    use super::*;
    use crate::ray::Ray;
    use tpt_optics_core::Vec3;

    #[test]
    fn quadric_seed_finds_off_axis_hit() {
        // Sphere R = 50 mm (c = 20 m⁻¹). A ray at height 30 mm hits the
        // near cap at z = 10 mm, so from z = −10 mm the parameter is 0.02.
        // A seed with radius 1/(2c) finds no root here.
        let ray = Ray::new(Vec3::new(0.03, 0.0, -0.01), Vec3::new(0.0, 0.0, 1.0));
        let hit = intersect_asphere::<f64>(20.0, 0.0, &[0.0; 8], &ray, 0.0);
        let (t, _, _) = hit.expect("valid off-axis intersection");
        assert!((t - 0.02).abs() < 1e-9, "t = {t}");
    }

    #[test]
    fn fnum_uses_diameter_not_semi_diameter() {
        // Same BK7-like singlet as the paraxial tests: R = ±100 mm, t = 5 mm.
        let mut sys = SequentialSystem::new();
        sys.media[0] = Medium::Const(1.0);
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
        let px = sys.paraxial(TraceContext::default());
        let (epd, fnum) = px.entrance_pupil(0).expect("stop at surface 0");
        let efl = px.effective_focal_length().unwrap();
        assert!((epd - 12e-3).abs() < 1e-12, "semi-diameter {epd}");
        assert!((fnum - efl / (2.0 * epd)).abs() < 1e-9, "F/# {fnum}");
    }
}

#[cfg(test)]
mod field_launch_tests {
    use super::*;

    #[test]
    fn field_ray_crosses_first_plane_at_pupil_point() {
        // The ray must pass through (x, 0, 0) for any field angle.
        for &(x, deg) in &[(1e-3, 0.0), (2e-3, 10.0), (-3e-3, 25.0)] {
            let r = field_ray(x, deg);
            let t = -r.origin.z / r.direction.z;
            let p = r.at(t);
            assert!((p.x - x).abs() < 1e-15, "x {} vs {x}", p.x);
            assert!(p.y.abs() < 1e-15, "y {} at field {deg}", p.y);
        }
    }

    #[test]
    fn flat_window_spot_is_independent_of_field_angle() {
        // A plane window is imaged onto itself: its spot must not depend on
        // the field angle. The old launch offset shifted the bundle by
        // 0.1·tanθ, which broke this.
        let mut sys = SequentialSystem::new();
        sys.push(SeqSurface::flat(12e-3), 0.0, Medium::Vacuum);
        let heights = [0.0, 1e-3, 2e-3, 3e-3];
        let on = rms_spot_radius(&sys, TraceContext::default(), 0.0, &heights);
        let off = rms_spot_radius(&sys, TraceContext::default(), 10.0, &heights);
        assert!((on - off).abs() < 1e-12, "on-axis {on} vs 10° {off}");
    }
}
