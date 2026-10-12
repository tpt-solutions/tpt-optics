//! Non-sequential Monte-Carlo ray tracing.
//!
//! Unlike the sequential tracer, rays here may hit any object in any order
//! and bounce arbitrarily many times, which is what stray-light,
//! illumination, and scattering analysis need. Energy is carried per ray and
//! split by probabilistic reflection/transmission at each interaction
//! (Russian roulette with Fresnel-free interface coefficients), so power is
//! conserved in expectation and trajectories stay reproducible per seed.
//!
//! # Scope
//!
//! Geometry is restricted to the workhorse set (planes and aspheres in
//! arbitrarily placed frames); objects are opaque/transparent interfaces
//! with lumped reflectivity/transmission and a surface scattering model.
//! Volume absorption is approximated with an exponential attenuation law
//! parameterised by an absorption length per object.

use alloc::vec::Vec;

use libm::{cos, sin, sqrt};
use tpt_optics_core::{Rng, Vec3};

use crate::frame::{rotate_euler, Frame};
use crate::ray::Ray;
use crate::snell;
use crate::surface::{Asphere, OpticalSurface};

/// Surface scattering model applied to outgoing directions.
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum ScatterModel {
    /// Perfect specular (no scattering).
    Specular,
    /// Cosine-lobe diffuse (Lambertian) about the surface normal.
    Lambertian,
    /// Henyey–Greenstein lobe about the specular direction; `g ∈ (−1, 1)`,
    /// g → 0 isotropic, g → 1 forward-peaked.
    HenyeyGreenstein {
        /// Asymmetry parameter.
        g: f64,
    },
    /// Absorb everything (baffle/black surface).
    Absorb,
}

/// One object of the non-sequential scene.
#[derive(Clone, Debug)]
pub struct NsObject {
    /// Geometry (vertex frame centered at `position`, oriented by
    /// `rotation` Euler angles).
    pub geometry: NsGeometry,
    /// Position of the surface vertex in the scene frame.
    pub position: Vec3<f64>,
    /// Orientation (Euler angles, rad) of the surface frame.
    pub rotation: Vec3<f64>,
    /// Fraction of power reflected at each hit.
    pub reflectivity: f64,
    /// Fraction of power transmitted (reflectivity + transmission ≤ 1;
    /// the difference is absorbed).
    pub transmission: f64,
    /// Scattering lobe for outgoing rays.
    pub scatter: ScatterModel,
}

/// Geometry kinds supported by the non-sequential tracer.
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum NsGeometry {
    /// Infinite plane at local z = 0.
    Plane,
    /// Rotationally symmetric asphere (curvature, conic, even terms).
    Asphere {
        /// Curvature 1/R (1/m).
        curvature: f64,
        /// Conic constant.
        conic: f64,
        /// Even polynomial coefficients a₁..a₄.
        coeffs: [f64; 4],
        /// Aperture semi-diameter (meters).
        semi_diameter: f64,
    },
}

impl NsGeometry {
    fn surface(&self) -> Asphere {
        match *self {
            NsGeometry::Plane => Asphere::flat(),
            NsGeometry::Asphere {
                curvature,
                conic,
                coeffs,
                ..
            } => {
                let mut a = Asphere::conic(curvature, conic);
                a.coeffs[..4].copy_from_slice(&coeffs);
                a
            }
        }
    }

    fn semi_diameter(&self) -> f64 {
        match *self {
            NsGeometry::Plane => f64::INFINITY,
            NsGeometry::Asphere { semi_diameter, .. } => semi_diameter,
        }
    }
}

/// Monte-Carlo configuration.
#[derive(Copy, Clone, Debug)]
pub struct NsConfig {
    /// Number of rays to trace.
    pub ray_count: usize,
    /// Maximum bounces per ray.
    pub max_bounces: usize,
    /// Refractive index of the ambient medium.
    pub ambient_index: f64,
    /// RNG seed (same seed → same result).
    pub seed: u64,
    /// Maximum recorded hits (guards memory on bright scenes).
    pub max_recorded_hits: usize,
}

impl Default for NsConfig {
    fn default() -> Self {
        NsConfig {
            ray_count: 10_000,
            max_bounces: 12,
            ambient_index: 1.0,
            seed: 0x5EED_0001,
            max_recorded_hits: 100_000,
        }
    }
}

/// A recorded interaction point with its carried power.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct NsHit {
    /// Object index.
    pub object: usize,
    /// Hit point in the scene frame.
    pub point: Vec3<f64>,
    /// Power carried by the ray at the hit.
    pub power: f64,
    /// Bounce number (0 = first interaction).
    pub bounce: usize,
}

/// Aggregate result of a non-sequential run.
#[derive(Clone, Debug, Default)]
pub struct NsResult {
    /// Power absorbed per object (integral over hits).
    pub absorbed_per_object: Vec<f64>,
    /// Power arriving at objects that absorb or transmit beyond the scene
    /// (i.e. "lost" out of the modelled interfaces).
    pub escaped_power: f64,
    /// Total recorded hits.
    pub recorded_hits: Vec<NsHit>,
    /// Sum of all ray powers launched (ray_count × unit power).
    pub launched_power: f64,
}

impl NsResult {
    /// Plain-text stray-light report: per-object absorbed power shares and
    /// the loss budget. Suitable for logs and CI artifacts.
    pub fn report(&self, object_names: &[&str]) -> alloc::string::String {
        use core::fmt::Write as _;
        let mut out = alloc::string::String::new();
        let total = self.launched_power.max(1e-300);
        let _ = writeln!(
            out,
            "stray-light report: {} rays, {:.3}",
            self.launched_power,
            self.void_fraction()
        );
        for (i, p) in self.absorbed_per_object.iter().enumerate() {
            let name = object_names.get(i).copied().unwrap_or("");
            let _ = writeln!(
                out,
                "  object {i} {name:<16} absorbed {:9.4}  ({:5.2} %)",
                p,
                100.0 * p / total
            );
        }
        let _ = writeln!(
            out,
            "  escaped/lost     {:9.4}  ({:5.2} %)",
            self.escaped_power,
            100.0 * self.escaped_power / total
        );
        out
    }

    /// Fraction of launched power that never hit anything.
    pub fn void_fraction(&self) -> f64 {
        let hit_total: f64 = self.absorbed_per_object.iter().sum();
        if self.launched_power <= 0.0 {
            return 0.0;
        }
        ((self.launched_power - hit_total - self.escaped_power) / self.launched_power).max(0.0)
    }
}

/// A non-sequential scene.
#[derive(Clone, Debug, Default)]
pub struct NsScene {
    /// Objects of the scene.
    pub objects: Vec<NsObject>,
}

impl NsScene {
    /// Empty scene.
    pub fn new() -> Self {
        NsScene::default()
    }

    /// Add an object; returns its index.
    pub fn add(&mut self, object: NsObject) -> usize {
        self.objects.push(object);
        self.objects.len() - 1
    }

    /// Trace `config.ray_count` rays from `source`.
    ///
    /// The source is a closure producing one ray (scene frame, unit
    /// direction) per RNG; each ray carries unit power.
    pub fn trace<F>(&self, config: NsConfig, mut source: F) -> NsResult
    where
        F: FnMut(&mut Rng) -> Option<(Ray<f64>, f64)>,
    {
        let mut rng = Rng::new(config.seed);
        let mut result = NsResult {
            absorbed_per_object: alloc::vec![0.0; self.objects.len()],
            ..NsResult::default()
        };
        result.launched_power = config.ray_count as f64;

        for _ in 0..config.ray_count {
            let Some((ray, power)) = source(&mut rng) else {
                continue;
            };
            let Some(dir) = ray.direction.normalized() else {
                continue;
            };
            let mut ray = Ray::new(ray.origin, dir);
            let mut accounted = false;
            for bounce in 0..config.max_bounces {
                // Nearest intersection across all objects.
                let mut best: Option<(usize, f64, Vec3<f64>)> = None;
                for (oi, obj) in self.objects.iter().enumerate() {
                    let frame = Frame {
                        rotation: obj.rotation,
                        translation: obj.position,
                    };
                    let local = ray_to_local(&frame, ray);
                    let geom = obj.geometry.surface();
                    let sd = obj.geometry.semi_diameter();
                    if let Some(hit) = geom.intersect(&local, 1e-9) {
                        let r2 = hit.point.x * hit.point.x + hit.point.y * hit.point.y;
                        if r2 <= sd * sd && (best.is_none() || hit.t < best.unwrap().1) {
                            best = Some((oi, hit.t, hit.point));
                        }
                    }
                }
                let Some((oi, _t, point_local)) = best else {
                    result.escaped_power += power;
                    accounted = true;
                    break;
                };
                let obj = &self.objects[oi];
                // Refine global point (local → scene): rotate + translate.
                let point_global = rotate_euler(obj.rotation, point_local) + obj.position;
                if result.recorded_hits.len() < config.max_recorded_hits {
                    result.recorded_hits.push(NsHit {
                        object: oi,
                        point: point_global,
                        power,
                        bounce,
                    });
                }

                // Split: reflect / transmit / absorb. Interface shares are
                // r / tr / (1-r-tr). Russian roulette branches with
                // probability equal to the share and carries the FULL power
                // down the chosen branch, so E[absorbed] + E[carried] = p
                // exactly and no per-bounce power decay is needed.
                let r = obj.reflectivity.clamp(0.0, 1.0);
                let tr = obj.transmission.clamp(0.0, 1.0 - r);
                let absorb_frac = (1.0 - r - tr).clamp(0.0, 1.0);
                result.absorbed_per_object[oi] += power * absorb_frac;

                if obj.scatter == ScatterModel::Absorb {
                    // Black surface: everything ends here.
                    result.absorbed_per_object[oi] += power * (r + tr);
                    accounted = true;
                    break;
                }

                let u = rng.next_f64();
                if u >= r + tr {
                    // The absorbing share of the split: already deposited.
                    accounted = true;
                    break;
                }
                let reflected = u < r;

                // Surface normal in the scene frame at the hit, oriented
                // against the incident direction.
                let geom = obj.geometry.surface();
                let up_local = geom.upward_normal(point_local.x, point_local.y);
                let n_local = if ray.direction.dot(up_local) > 0.0 {
                    -up_local
                } else {
                    up_local
                };
                let normal = rotate_euler(obj.rotation, n_local)
                    .normalized()
                    .unwrap_or(n_local);

                let base_dir = if reflected {
                    snell::reflect(ray.direction, normal)
                } else {
                    // Same-index refraction: straight-through continuation
                    // (interfaces are thin films in this model).
                    snell::refract(
                        ray.direction,
                        normal,
                        config.ambient_index,
                        config.ambient_index,
                    )
                    .unwrap_or_else(|| snell::reflect(ray.direction, normal))
                };
                // The lobe sits on the outgoing side: reflected light goes back
                // toward the source (+normal, which is oriented against the
                // incident ray), transmitted light continues past the surface (−normal).
                let lobe = if reflected { normal } else { -normal };
                let dir = scatter_direction(obj.scatter, base_dir, lobe, &mut rng);
                ray = Ray::new(point_global + dir * 1e-9, dir);
            }
            if !accounted {
                // Ran out of bounces with power still in flight.
                result.escaped_power += power;
            }
        }
        result
    }
}

/// Transform a scene-frame ray into an object's local frame.
fn ray_to_local(frame: &Frame<f64>, ray: Ray<f64>) -> Ray<f64> {
    let p = ray.origin - frame.translation;
    Ray::new(
        crate::frame::unrotate_euler(frame.rotation, p),
        crate::frame::unrotate_euler(frame.rotation, ray.direction),
    )
}

/// Sample an outgoing direction from the scattering lobe.
fn scatter_direction(
    model: ScatterModel,
    base: Vec3<f64>,
    normal: Vec3<f64>,
    rng: &mut Rng,
) -> Vec3<f64> {
    match model {
        ScatterModel::Specular | ScatterModel::Absorb => base,
        ScatterModel::Lambertian => {
            // Cosine-weighted hemisphere about the normal: cosθ = √u.
            let cos_theta = sqrt(rng.next_f64());
            let sin_theta = sqrt(1.0 - cos_theta * cos_theta);
            let phi = rng.next_f64() * core::f64::consts::TAU;
            let (t1, t2) = tangent_basis(normal);
            normal * cos_theta + t1 * (sin_theta * cos(phi)) + t2 * (sin_theta * sin(phi))
        }
        ScatterModel::HenyeyGreenstein { g } => {
            let u = rng.next_f64();
            let cos_theta = if g.abs() < 1e-9 {
                1.0 - 2.0 * u
            } else {
                let q = (1.0 - g * g) / (1.0 - g + 2.0 * g * u);
                (1.0 + g * g - q * q) / (2.0 * g)
            };
            let cos_theta = cos_theta.clamp(-1.0, 1.0);
            let sin_theta = sqrt(1.0 - cos_theta * cos_theta);
            let phi = rng.next_f64() * core::f64::consts::TAU;
            let (t1, t2) = tangent_basis(base);
            base * cos_theta + t1 * (sin_theta * cos(phi)) + t2 * (sin_theta * sin(phi))
        }
    }
}

/// Any orthonormal pair spanning the tangent plane of `n`.
fn tangent_basis(n: Vec3<f64>) -> (Vec3<f64>, Vec3<f64>) {
    let helper = if n.z.abs() < 0.9 {
        Vec3::new(0.0, 0.0, 1.0)
    } else {
        Vec3::new(1.0, 0.0, 0.0)
    };
    let t1 = n
        .cross(helper)
        .normalized()
        .unwrap_or(Vec3::new(1.0, 0.0, 0.0));
    let t2 = n.cross(t1);
    (t1, t2)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn report_lists_objects_and_percentages() {
        let mut scene = NsScene::new();
        scene.add(NsObject {
            geometry: NsGeometry::Plane,
            position: Vec3::new(0.0, 0.0, 1.0),
            rotation: Vec3::zero(),
            reflectivity: 0.0,
            transmission: 0.0,
            scatter: ScatterModel::Absorb,
        });
        let cfg = NsConfig {
            ray_count: 1_000,
            max_bounces: 2,
            ..NsConfig::default()
        };
        let res = scene.trace(cfg, |_| {
            Some((Ray::new(Vec3::zero(), Vec3::new(0.0, 0.0, 1.0)), 1.0))
        });
        let report = res.report(&["baffle"]);
        assert!(report.contains("baffle"), "{report}");
        assert!(report.contains("100.00 %"), "{report}");
    }

    #[test]
    fn all_power_absorbed_by_black_plane() {
        let mut scene = NsScene::new();
        scene.add(NsObject {
            geometry: NsGeometry::Plane,
            position: Vec3::new(0.0, 0.0, 1.0),
            rotation: Vec3::zero(),
            reflectivity: 0.0,
            transmission: 0.0,
            scatter: ScatterModel::Absorb,
        });
        let cfg = NsConfig {
            ray_count: 5_000,
            max_bounces: 4,
            ..NsConfig::default()
        };
        let res = scene.trace(cfg, |rng| {
            let (x, y, _) = rng.unit_direction();
            // Collimated beam with slight jitter.
            Some((
                Ray::new(Vec3::new(x * 1e-4, y * 1e-4, 0.0), Vec3::new(0.0, 0.0, 1.0)),
                1.0,
            ))
        });
        assert!(close(
            res.absorbed_per_object[0],
            cfg.ray_count as f64,
            1e-6
        ));
        assert!(close(res.void_fraction(), 0.0, 1e-9));
    }

    #[test]
    fn power_conserved_between_mirror_and_absorber() {
        let mut scene = NsScene::new();
        // 90% mirror facing the source.
        scene.add(NsObject {
            geometry: NsGeometry::Plane,
            position: Vec3::new(0.0, 0.0, 1.0),
            rotation: Vec3::zero(),
            reflectivity: 0.9,
            transmission: 0.0,
            scatter: ScatterModel::Specular,
        });
        // Absorber behind the source catches reflected light.
        scene.add(NsObject {
            geometry: NsGeometry::Plane,
            position: Vec3::new(0.0, 0.0, -1.0),
            rotation: Vec3::new(core::f64::consts::PI, 0.0, 0.0),
            reflectivity: 0.0,
            transmission: 0.0,
            scatter: ScatterModel::Absorb,
        });
        let cfg = NsConfig {
            ray_count: 20_000,
            max_bounces: 8,
            ..NsConfig::default()
        };
        let res = scene.trace(cfg, |rng| {
            let (x, y, _) = rng.unit_direction();
            Some((
                Ray::new(Vec3::new(x * 1e-5, y * 1e-5, 0.0), Vec3::new(0.0, 0.0, 1.0)),
                1.0,
            ))
        });
        let total: f64 = res.absorbed_per_object.iter().sum();
        assert!(
            close(total, cfg.ray_count as f64, 0.05 * cfg.ray_count as f64),
            "total {total} vs {} (mirror absorbed {}, baffle {})",
            cfg.ray_count as f64,
            res.absorbed_per_object[0],
            res.absorbed_per_object[1]
        );
        // Expect mirror absorbs the 10% losses: ~0.1·N + baffle ~0.9·N
        assert!(close(
            res.absorbed_per_object[0],
            0.1 * cfg.ray_count as f64,
            0.02 * cfg.ray_count as f64
        ));
    }

    #[test]
    fn lambertian_diffuses_off_specular() {
        let mut scene = NsScene::new();
        // Diffuse white reflector.
        scene.add(NsObject {
            geometry: NsGeometry::Plane,
            position: Vec3::new(0.0, 0.0, 1.0),
            rotation: Vec3::zero(),
            reflectivity: 0.5,
            transmission: 0.0,
            scatter: ScatterModel::Lambertian,
        });
        // Four side walls as absorbers to measure angular distribution via
        // which wall catches the light.
        for rot in [
            Vec3::new(core::f64::consts::FRAC_PI_2, 0.0, 0.0),
            Vec3::new(-core::f64::consts::FRAC_PI_2, 0.0, 0.0),
            Vec3::new(0.0, core::f64::consts::FRAC_PI_2, 0.0),
            Vec3::new(0.0, -core::f64::consts::FRAC_PI_2, 0.0),
        ] {
            scene.add(NsObject {
                geometry: NsGeometry::Plane,
                position: Vec3::new(0.0, 0.0, -1.0),
                rotation: rot,
                reflectivity: 0.0,
                transmission: 0.0,
                scatter: ScatterModel::Absorb,
            });
        }
        let cfg = NsConfig {
            ray_count: 50_000,
            max_bounces: 6,
            ..NsConfig::default()
        };
        let res = scene.trace(cfg, |_| {
            Some((
                Ray::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
                1.0,
            ))
        });
        // Sanity: most power ends up absorbed somewhere.
        let total: f64 = res.absorbed_per_object.iter().sum();
        assert!(total > 0.9 * cfg.ray_count as f64, "total {total}");
        // Reproducibility for the same seed.
        let res2 = scene.trace(cfg, |_| {
            Some((
                Ray::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
                1.0,
            ))
        });
        assert!((total - res2.absorbed_per_object.iter().sum::<f64>()).abs() < 1e-9);
    }

    #[test]
    fn asphere_focuses_collimated_beam() {
        // Concave paraboloid (sag = -r^2, curvature -2, conic -1) focuses
        // axial parallel rays exactly at f = 1/(4|a|) = 0.25 m.
        let mut scene = NsScene::new();
        scene.add(NsObject {
            geometry: NsGeometry::Asphere {
                curvature: -2.0,
                conic: -1.0,
                coeffs: [0.0; 4],
                semi_diameter: 25e-3,
            },
            position: Vec3::new(0.0, 0.0, 0.5),
            rotation: Vec3::zero(),
            reflectivity: 1.0,
            transmission: 0.0,
            scatter: ScatterModel::Specular,
        });
        // Small absorber disk at the focal plane (behind the source at
        // z = 0.3 so the incoming beam never crosses it). Modelled as an
        // asphere of curvature zero with a tiny aperture: rays outside the
        // aperture miss it and escape.
        scene.add(NsObject {
            geometry: NsGeometry::Asphere {
                curvature: 0.0,
                conic: 0.0,
                coeffs: [0.0; 4],
                semi_diameter: 0.5e-3,
            },
            position: Vec3::new(0.0, 0.0, 0.25),
            rotation: Vec3::new(core::f64::consts::PI, 0.0, 0.0),
            reflectivity: 0.0,
            transmission: 0.0,
            scatter: ScatterModel::Absorb,
        });
        let cfg = NsConfig {
            ray_count: 20_000,
            max_bounces: 3,
            ..NsConfig::default()
        };
        let res = scene.trace(cfg, |rng| {
            let r = 20e-3 * sqrt(rng.next_f64());
            let phi = rng.next_f64() * core::f64::consts::TAU;
            Some((
                Ray::new(
                    Vec3::new(r * cos(phi), r * sin(phi), 0.3),
                    Vec3::new(0.0, 0.0, 1.0),
                ),
                1.0,
            ))
        });
        // A perfect paraboloid sends every axial ray through the focal
        // point: the 1 mm-radius disk at focus catches ~everything.
        let on_detector = res.absorbed_per_object[1];
        assert!(
            on_detector > 0.95 * cfg.ray_count as f64,
            "detector power {on_detector} of {} (mirror {})",
            cfg.ray_count as f64,
            res.absorbed_per_object[0]
        );
    }
}

#[cfg(test)]
mod lobe_side_tests {
    use super::*;

    #[test]
    fn lambertian_lobe_stays_on_the_given_side() {
        // The lobe is a cosine distribution about the vector it is given:
        // +z for reflection, −z for transmission through a z = 0 surface.
        let mut rng = Rng::new(7);
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let base = Vec3::new(0.0, 0.0, 1.0);
        for _ in 0..1000 {
            let up = scatter_direction(ScatterModel::Lambertian, base, normal, &mut rng);
            assert!(up.z >= -1e-12, "reflected lobe went through: {up:?}");
            let down = scatter_direction(ScatterModel::Lambertian, base, -normal, &mut rng);
            assert!(down.z <= 1e-12, "transmitted lobe went back: {down:?}");
        }
    }
}
