//! One-call analysis of a sequential system: paraxial properties, and per
//! field the spot, ray fan, OPD, distortion, wavefront error and best-focus
//! shift. The report exports to CSV, JSON, and SVG spot diagrams.
//!
//! The functions it calls are all available on their own. This module gives
//! them one entry point and one report type. The image plane is the last
//! surface, so the system should end with a flat image plane (as the presets
//! do).

use alloc::string::String;
use alloc::vec::Vec;
use core::fmt::Write as _;

use tpt_optics_core::Vec3;
use tpt_optics_geometric::ray::Ray;
use tpt_optics_geometric::sequential::{
    field_ray, rms_spot_radius, SequentialSystem, TraceContext,
};

use crate::merit::{chromatic_focal_shift, distortion};
use crate::mtf::{system_mtf, MtfCurves};
use crate::wavefront::{hexapolar_grid, wavefront_rms, wavefront_samples};

/// Zernike order used for the wavefront fit (OSA index up to 8).
const WAVEFRONT_MAX_J: u32 = 8;

/// Measurements at one field angle.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct FieldReport {
    /// Field angle (degrees).
    pub angle_deg: f64,
    /// RMS transverse spot radius over the pupil fan (m).
    pub rms_spot: f64,
    /// Distortion at this field, as a fraction of the paraxial image height.
    pub distortion: f64,
    /// RMS wavefront error after a Zernike fit (m). `NaN` if no ray reached
    /// the image.
    pub wavefront_rms: f64,
    /// Best-focus position relative to the on-axis best focus (m). Positive
    /// means the image plane moves away from the last surface. This is a
    /// through-focus search with about 0.1 mm resolution for typical lenses.
    pub best_focus_shift: f64,
    /// Image-plane intersection of each unblocked ray (x, y in m).
    pub spot: Vec<(f64, f64)>,
    /// Sagittal ray fan: (pupil height along x, image-plane x) for each
    /// unblocked ray (m).
    pub ray_fan: Vec<(f64, f64)>,
    /// Tangential (meridional) ray fan: (pupil height along y, image-plane y)
    /// for each unblocked ray (m). The field angle acts in this plane.
    pub ray_fan_tangential: Vec<(f64, f64)>,
    /// OPD of each unblocked ray relative to the on-axis chief ray (m),
    /// one entry per entry of `ray_fan`.
    pub opd: Vec<f64>,
    /// System MTF along the tangential and sagittal axes, or `None` when no focus
    /// or no wavefront is available. Frequencies are normalised to the cutoff.
    pub mtf: Option<MtfCurves>,
}

/// Summary of a sequential system. Quantities that do not exist for the
/// system (for example a focal length for an afocal system) are `None`.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct AnalysisReport {
    /// Effective focal length (m).
    pub efl: Option<f64>,
    /// Back focal length (m).
    pub bfl: Option<f64>,
    /// Working F-number, from the entrance pupil of the stop surface.
    pub f_number: Option<f64>,
    /// Spread of the paraxial focus across the wavelength set (m).
    pub chromatic_shift: Option<f64>,
    /// One entry per requested field angle, in the order given.
    pub fields: Vec<FieldReport>,
}

fn json_number(v: f64) -> String {
    if v.is_finite() {
        alloc::format!("{v}")
    } else {
        String::from("null")
    }
}

fn json_option(v: Option<f64>) -> String {
    v.map_or_else(|| String::from("null"), json_number)
}

impl AnalysisReport {
    /// CSV of the per-field scalars, with a header line.
    #[must_use]
    pub fn to_csv(&self) -> String {
        let mut out =
            String::from("field_deg,rms_spot_m,distortion,wavefront_rms_m,best_focus_shift_m\n");
        for f in &self.fields {
            let _ = writeln!(
                out,
                "{},{},{},{},{}",
                f.angle_deg, f.rms_spot, f.distortion, f.wavefront_rms, f.best_focus_shift
            );
        }
        out
    }

    /// JSON of the whole report. Non-finite numbers are written as `null`.
    #[must_use]
    pub fn to_json(&self) -> String {
        let mut out = String::from("{");
        let _ = write!(
            out,
            "\"efl\":{},\"bfl\":{},\"f_number\":{},\"chromatic_shift\":{},\"fields\":[",
            json_option(self.efl),
            json_option(self.bfl),
            json_option(self.f_number),
            json_option(self.chromatic_shift)
        );
        for (i, f) in self.fields.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            let _ = write!(
                out,
                "{{\"angle_deg\":{},\"rms_spot\":{},\"distortion\":{},\"wavefront_rms\":{},\"best_focus_shift\":{},\"ray_fan\":[",
                json_number(f.angle_deg),
                json_number(f.rms_spot),
                json_number(f.distortion),
                json_number(f.wavefront_rms),
                json_number(f.best_focus_shift)
            );
            for (k, (h, y)) in f.ray_fan.iter().enumerate() {
                if k > 0 {
                    out.push(',');
                }
                let _ = write!(out, "[{},{}]", json_number(*h), json_number(*y));
            }
            out.push_str("],\"opd\":[");
            for (k, v) in f.opd.iter().enumerate() {
                if k > 0 {
                    out.push(',');
                }
                out.push_str(&json_number(*v));
            }
            out.push_str("],\"mtf\":");
            match &f.mtf {
                Some(m) => {
                    out.push_str("{\"tangential\":[");
                    for (k, (nu, v)) in m.tangential.iter().enumerate() {
                        if k > 0 {
                            out.push(',');
                        }
                        let _ = write!(out, "[{},{}]", json_number(*nu), json_number(*v));
                    }
                    out.push_str("],\"sagittal\":[");
                    for (k, (nu, v)) in m.sagittal.iter().enumerate() {
                        if k > 0 {
                            out.push(',');
                        }
                        let _ = write!(out, "[{},{}]", json_number(*nu), json_number(*v));
                    }
                    out.push_str("]}");
                }
                None => out.push_str("null"),
            }
            out.push('}');
        }
        out.push_str("]}");
        out
    }

    /// SVG ray-fan plot for field `index`: image-plane height (horizontal) against
    /// pupil height (vertical), with a polyline through the unblocked rays.
    /// Returns `None` for a bad index or an empty fan.
    #[must_use]
    pub fn ray_fan_svg(&self, index: usize) -> Option<String> {
        let field = self.fields.get(index)?;
        if field.ray_fan.len() < 2 {
            return None;
        }
        let size = 200.0;
        let pad = 12.0;
        let span = |values: &mut dyn Iterator<Item = f64>| {
            let (lo, hi) = values.fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), v| {
                (l.min(v), h.max(v))
            });
            if hi - lo > 0.0 {
                (lo, hi)
            } else {
                (lo - 1.0, hi + 1.0)
            }
        };
        let (x_lo, x_hi) = span(&mut field.ray_fan.iter().map(|p| p.1));
        let (y_lo, y_hi) = span(&mut field.ray_fan.iter().map(|p| p.0));
        let to_px = |x: f64, y: f64| {
            let px = pad + (x - x_lo) / (x_hi - x_lo) * (size - 2.0 * pad);
            let py = size - pad - (y - y_lo) / (y_hi - y_lo) * (size - 2.0 * pad);
            (px, py)
        };
        let mut out = alloc::format!(
            "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{size}\" height=\"{size}\" viewBox=\"0 0 {size} {size}\">\n"
        );
        out.push_str("<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n<polyline fill=\"none\" stroke=\"black\" points=\"");
        for (h, x) in &field.ray_fan {
            let (px, py) = to_px(*x, *h);
            let _ = write!(out, "{px:.2},{py:.2} ");
        }
        out.push_str("\"/>\n</svg>\n");
        Some(out)
    }

    /// SVG spot diagram for field `index`, centred on the spot centroid and
    /// scaled to the largest spot radius. Returns `None` for a bad index or
    /// an empty spot.
    #[must_use]
    pub fn spot_svg(&self, index: usize) -> Option<String> {
        let field = self.fields.get(index)?;
        if field.spot.is_empty() {
            return None;
        }
        let n = field.spot.len() as f64;
        let cx = field.spot.iter().map(|p| p.0).sum::<f64>() / n;
        let cy = field.spot.iter().map(|p| p.1).sum::<f64>() / n;
        let radius = field
            .spot
            .iter()
            .map(|p| libm::hypot(p.0 - cx, p.1 - cy))
            .fold(0.0_f64, f64::max)
            .max(1e-12);
        let size = 200.0;
        let scale = (size / 2.0 - 10.0) / radius;
        let mut out = alloc::format!(
            "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{size}\" height=\"{size}\" viewBox=\"0 0 {size} {size}\">\n"
        );
        out.push_str("<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n");
        for (x, y) in &field.spot {
            let px = size / 2.0 + (x - cx) * scale;
            let py = size / 2.0 - (y - cy) * scale;
            let _ = writeln!(
                out,
                "<circle cx=\"{px:.2}\" cy=\"{py:.2}\" r=\"2\" fill=\"black\"/>"
            );
        }
        out.push_str("</svg>\n");
        Some(out)
    }
}

/// Best-focus image offset (m) for a field, from a two-stage through-focus
/// scan of the RMS spot. The offset is relative to the nominal image plane.
fn best_focus(
    system: &SequentialSystem,
    ctx: TraceContext,
    angle_deg: f64,
    heights: &[f64],
    nominal_efl: Option<f64>,
) -> f64 {
    let n = system.thickness.len();
    if n < 2 {
        return f64::NAN;
    }
    // The paraxial focus lies `bfl` beyond the image plane, so the scan is
    // centred there. A scan centred on the nominal plane misses the focus
    // whenever the design is far from focus.
    let center = system.paraxial(ctx).back_focal_length().unwrap_or(0.0);
    let range = nominal_efl.map_or(5e-3, |f| 0.1 * libm::fabs(f)).max(1e-6);
    let spot_at = |dz: f64| {
        let mut trial = system.clone();
        trial.thickness[n - 2] += dz;
        rms_spot_radius(&trial, ctx, angle_deg, heights)
    };
    let argmin = |center: f64, half: f64, steps: usize| {
        let mut best = (f64::INFINITY, center);
        for k in 0..=steps {
            let dz = center - half + 2.0 * half * k as f64 / steps as f64;
            let s = spot_at(dz);
            if s.is_finite() && s < best.0 {
                best = (s, dz);
            }
        }
        best.1
    };
    let coarse = argmin(center, range, 40);
    argmin(coarse, range / 20.0, 20)
}

/// Analyse `system` at the reference context `ctx`.
///
/// * `stop_index` — surface that limits the aperture, used for the F-number,
/// * `pupil_heights` — ray heights (m) of the pupil fan for each field,
///   including `0.0` if the on-axis chief ray is wanted for the OPD,
/// * `fields_deg` — field angles (degrees) to evaluate,
/// * `wavelength_contexts` — contexts for the chromatic shift. Pass an empty
///   slice to skip it.
pub fn analyze(
    system: &SequentialSystem,
    ctx: TraceContext,
    stop_index: usize,
    pupil_heights: &[f64],
    fields_deg: &[f64],
    wavelength_contexts: &[TraceContext],
) -> AnalysisReport {
    let paraxial = system.paraxial(ctx);
    let efl = paraxial.effective_focal_length();
    let bfl = paraxial.back_focal_length();
    let f_number = paraxial
        .entrance_pupil(stop_index)
        .and_then(|(_, fnum)| fnum);
    let chromatic_shift = if wavelength_contexts.is_empty() {
        None
    } else {
        // NaN means fewer than two wavelengths have a focus; report that as
        // "not computed" rather than as a zero shift.
        Some(chromatic_focal_shift(system, wavelength_contexts)).filter(|s| s.is_finite())
    };

    let aperture = pupil_heights
        .iter()
        .fold(0.0_f64, |m, h| m.max(libm::fabs(*h)));
    let grid = hexapolar_grid(4);
    let on_axis_focus = best_focus(system, ctx, 0.0, pupil_heights, efl);

    let fields = fields_deg
        .iter()
        .map(|&angle_deg| {
            let mut spot = Vec::new();
            let mut ray_fan = Vec::new();
            let mut ray_fan_tangential = Vec::new();
            let mut opd = Vec::new();
            // OPD is relative to this field's own chief ray (pupil height 0).
            let chief_opd = system.trace(&field_ray(0.0, angle_deg), ctx).opd;
            let (sin_t, cos_t) = (
                libm::sin(angle_deg.to_radians()),
                libm::cos(angle_deg.to_radians()),
            );
            for &h in pupil_heights {
                let res = system.trace(&field_ray(h, angle_deg), ctx);
                if res.blocked {
                    continue;
                }
                if let Some(p) = res.points.last() {
                    spot.push((p.x, p.y));
                    ray_fan.push((h, p.x));
                    opd.push(res.opd - chief_opd);
                }
                // The tangential ray has its pupil offset along y, in the plane
                // the field tilts in. Launched like `field_ray`, it crosses the
                // first vertex plane at (0, h, 0).
                let origin = Vec3::new(0.0, h - 0.1 * sin_t / cos_t, -0.1);
                let tangential = system.trace(&Ray::new(origin, Vec3::new(0.0, sin_t, cos_t)), ctx);
                if !tangential.blocked {
                    if let Some(p) = tangential.points.last() {
                        ray_fan_tangential.push((h, p.y));
                    }
                }
            }
            let wavefront = if aperture > 0.0 {
                let samples = wavefront_samples(system, ctx, angle_deg, &grid, aperture);
                wavefront_rms(&samples, WAVEFRONT_MAX_J)
            } else {
                f64::NAN
            };
            let focus = best_focus(system, ctx, angle_deg, pupil_heights, efl);
            FieldReport {
                angle_deg,
                rms_spot: rms_spot_radius(system, ctx, angle_deg, pupil_heights),
                distortion: distortion(system, ctx, angle_deg),
                wavefront_rms: wavefront,
                best_focus_shift: focus - on_axis_focus,
                spot,
                ray_fan,
                ray_fan_tangential,
                opd,
                mtf: system_mtf(system, ctx, angle_deg, aperture, 64),
            }
        })
        .collect();
    AnalysisReport {
        efl,
        bfl,
        f_number,
        chromatic_shift,
        fields,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_optics_core::material::MaterialDatabase;
    use tpt_optics_core::{Kelvin, Wavelength};
    use tpt_optics_geometric::presets::singlet;

    fn ctx(nm: f64) -> TraceContext {
        TraceContext {
            wavelength: Wavelength::from_nm(nm),
            temperature: Kelvin::from_celsius(20.0),
        }
    }

    fn heights() -> Vec<f64> {
        (0..=10).map(|i| (i as f64 - 5.0) * 1e-3).collect()
    }

    fn lens() -> SequentialSystem {
        let catalog = MaterialDatabase::presets();
        singlet(&catalog, 100e-3, -100e-3, 12.5e-3, 5e-3, "bk7", 95e-3).unwrap()
    }

    #[test]
    fn report_matches_the_individual_measurements() {
        let sys = lens();
        let report = analyze(&sys, ctx(587.6), 0, &heights(), &[0.0, 3.0], &[]);

        let efl = sys.paraxial(ctx(587.6)).effective_focal_length().unwrap();
        assert_eq!(report.efl, Some(efl));
        assert_eq!(report.fields.len(), 2);
        assert_eq!(report.fields[0].angle_deg, 0.0);
        let direct = rms_spot_radius(&sys, ctx(587.6), 3.0, &heights());
        assert_eq!(report.fields[1].rms_spot, direct);
        assert!(report.chromatic_shift.is_none());
    }

    #[test]
    fn chromatic_shift_is_reported_when_contexts_are_given() {
        let sys = lens();
        let contexts = [ctx(486.1), ctx(656.3)];
        let report = analyze(&sys, ctx(587.6), 0, &heights(), &[0.0], &contexts);
        let shift = report.chromatic_shift.expect("shift requested");
        // Normal dispersion: the blue focus is shorter, so the spread is positive.
        assert!(shift > 0.0 && shift < 5e-3, "{shift}");
    }

    #[test]
    fn ray_fan_and_opd_have_one_entry_per_unblocked_ray() {
        let report = analyze(&lens(), ctx(587.6), 0, &heights(), &[0.0], &[]);
        let f = &report.fields[0];
        assert_eq!(f.ray_fan.len(), heights().len());
        assert_eq!(f.opd.len(), f.ray_fan.len());
        // The chief ray (height 0, index 5) is its own reference.
        assert!(f.opd[5].abs() < 1e-12, "{}", f.opd[5]);
        // The on-axis fan is symmetric in height.
        assert!((f.ray_fan[0].1 + f.ray_fan[10].1).abs() < 1e-12);
        // Guard against a fan that is zero everywhere, which would pass the symmetry check.
        assert!(
            f.ray_fan.iter().any(|(_, x)| x.abs() > 1e-6),
            "{:?}",
            f.ray_fan
        );
    }

    #[test]
    fn wavefront_and_focus_are_finite_on_axis() {
        let report = analyze(&lens(), ctx(587.6), 0, &heights(), &[0.0, 3.0], &[]);
        let on_axis = &report.fields[0];
        assert!(on_axis.wavefront_rms.is_finite() && on_axis.wavefront_rms >= 0.0);
        assert!(on_axis.best_focus_shift.abs() < 1e-9);
        assert!(report.fields[1].best_focus_shift.is_finite());
    }

    #[test]
    fn spot_svg_has_one_circle_per_ray() {
        let report = analyze(&lens(), ctx(587.6), 0, &heights(), &[0.0], &[]);
        let svg = report.spot_svg(0).expect("spot present");
        assert!(svg.starts_with("<svg"));
        assert_eq!(svg.matches("<circle").count(), heights().len());
        assert!(report.spot_svg(9).is_none());
        let fan = report.ray_fan_svg(0).expect("fan present");
        assert!(fan.contains("<polyline"));
    }

    #[test]
    fn json_and_csv_carry_the_same_fields() {
        let report = analyze(&lens(), ctx(587.6), 0, &heights(), &[0.0, 2.0], &[]);
        let json = report.to_json();
        assert!(json.starts_with('{') && json.ends_with('}'));
        assert!(json.contains("\"efl\":"));
        assert!(json.contains("\"ray_fan\":[["));
        assert!(json.contains("\"mtf\":{\"tangential\":[["), "{json}");
        let csv = report.to_csv();
        assert_eq!(csv.lines().count(), 3);
    }
}

#[cfg(test)]
mod field_opd_tests {
    use super::*;
    use tpt_optics_core::material::MaterialDatabase;
    use tpt_optics_core::{Kelvin, Wavelength};
    use tpt_optics_geometric::presets::singlet;

    #[test]
    fn chief_ray_opd_is_zero_at_every_field() {
        // OPD is measured against each field's own chief ray. The chief ray
        // at 10° has a different OPD from the on-axis one, so measuring
        // against the axis gave a nonzero value for the chief ray itself.
        let catalog = MaterialDatabase::presets();
        let sys = singlet(&catalog, 100e-3, -100e-3, 12.5e-3, 5e-3, "bk7", 95e-3).unwrap();
        let ctx = TraceContext {
            wavelength: Wavelength::from_nm(587.6),
            temperature: Kelvin::from_celsius(20.0),
        };
        let heights: Vec<f64> = (0..=10).map(|i| (i as f64 - 5.0) * 1e-3).collect();
        let report = analyze(&sys, ctx, 0, &heights, &[2.0], &[]);
        let field = &report.fields[0];
        assert_eq!(field.opd.len(), heights.len(), "some rays were blocked");
        assert!(field.opd[5].abs() < 1e-12, "chief OPD {}", field.opd[5]);
    }
}

#[cfg(test)]
mod tangential_fan_tests {
    use super::*;
    use tpt_optics_core::material::MaterialDatabase;
    use tpt_optics_core::{Kelvin, Wavelength};
    use tpt_optics_geometric::presets::singlet;

    fn report(angle: f64) -> FieldReport {
        let catalog = MaterialDatabase::presets();
        let sys = singlet(&catalog, 100e-3, -100e-3, 12.5e-3, 5e-3, "bk7", 95e-3).unwrap();
        let ctx = TraceContext {
            wavelength: Wavelength::from_nm(587.6),
            temperature: Kelvin::from_celsius(20.0),
        };
        let heights: Vec<f64> = (0..=10).map(|i| (i as f64 - 5.0) * 1e-3).collect();
        analyze(&sys, ctx, 0, &heights, &[angle], &[])
            .fields
            .remove(0)
    }

    #[test]
    fn on_axis_tangential_fan_matches_the_sagittal_one() {
        // A rotationally symmetric lens on axis: the two fans are the same
        // curve, one in x and one in y.
        let f = report(0.0);
        assert_eq!(f.ray_fan_tangential.len(), f.ray_fan.len());
        for (t, s) in f.ray_fan_tangential.iter().zip(&f.ray_fan) {
            assert!((t.1 - s.1).abs() < 1e-12, "{t:?} vs {s:?}");
        }
    }

    #[test]
    fn off_axis_tangential_fan_differs_from_the_sagittal_one() {
        // Off axis the meridional and sagittal aberrations differ, so the
        // tangential fan is not a copy of the sagittal one.
        let f = report(2.0);
        let max_diff = f
            .ray_fan_tangential
            .iter()
            .zip(&f.ray_fan)
            .map(|(t, s)| (t.1 - s.1).abs())
            .fold(0.0_f64, f64::max);
        assert!(max_diff > 1e-6, "{max_diff}");
    }
}
