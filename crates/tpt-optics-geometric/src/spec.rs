//! A serialisable form of a sequential system.
//!
//! Glasses are stored by name and resolved against a [`MaterialDatabase`] on
//! load, so a file does not carry dispersion data. Values that JSON cannot
//! represent (infinite apertures and gaps) are stored as `None`.

use alloc::string::String;
use alloc::vec::Vec;
use core::fmt;

use tpt_optics_core::material::MaterialDatabase;

use crate::sequential::{Medium, SeqSurface, SequentialSystem, SurfaceAction};
use crate::surface::Asphere;

/// Why a [`SystemSpec`] could not become a system.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SpecError {
    /// A glass name is not in the catalogue.
    UnknownMaterial(String),
    /// The surfaces, gaps and media do not fit together.
    Shape(&'static str),
}

impl fmt::Display for SpecError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnknownMaterial(name) => write!(f, "material '{name}' is not in the catalogue"),
            Self::Shape(msg) => write!(f, "system shape: {msg}"),
        }
    }
}

impl core::error::Error for SpecError {}

/// Surface behaviour, as in [`SurfaceAction`].
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum ActionSpec {
    /// Refract.
    Refract,
    /// Reflect.
    Reflect,
    /// Refract with a diffractive phase.
    Diffract,
}

/// Medium after a surface or before the first one.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum MediumSpec {
    /// Vacuum (index 1).
    Vacuum,
    /// A constant index.
    Const(f64),
    /// A catalogue material, by name.
    Material(String),
}

/// One surface.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SurfaceSpec {
    /// What happens to a ray at the surface.
    pub action: ActionSpec,
    /// Curvature (1/m).
    pub curvature: f64,
    /// Conic constant.
    pub conic: f64,
    /// Even asphere coefficients (eight values).
    pub coeffs: Vec<f64>,
    /// Clear semi-diameter (m); `None` for an unlimited aperture.
    pub semi_diameter: Option<f64>,
    /// Diffractive phase coefficients (nine values), if any.
    pub phase: Option<Vec<f64>>,
}

/// A whole system.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SystemSpec {
    /// Surfaces, object side first.
    pub surfaces: Vec<SurfaceSpec>,
    /// Gap after each surface (m); `None` for an infinite gap.
    pub thickness: Vec<Option<f64>>,
    /// Media: the medium before the first surface, then the one after each.
    pub media: Vec<MediumSpec>,
}

fn finite(value: f64) -> Option<f64> {
    value.is_finite().then_some(value)
}

impl SystemSpec {
    /// Describe a system. Materials are stored by name.
    #[must_use]
    pub fn from_system(sys: &SequentialSystem) -> Self {
        SystemSpec {
            surfaces: sys
                .surfaces
                .iter()
                .map(|s| SurfaceSpec {
                    action: match s.action {
                        SurfaceAction::Refract => ActionSpec::Refract,
                        SurfaceAction::Reflect => ActionSpec::Reflect,
                        SurfaceAction::Diffract => ActionSpec::Diffract,
                    },
                    curvature: s.asphere.curvature,
                    conic: s.asphere.conic,
                    coeffs: s.asphere.coeffs.to_vec(),
                    semi_diameter: finite(s.semi_diameter),
                    phase: s.phase.map(|p| p.to_vec()),
                })
                .collect(),
            thickness: sys.thickness.iter().map(|t| finite(*t)).collect(),
            media: sys
                .media
                .iter()
                .map(|m| match m {
                    Medium::Vacuum => MediumSpec::Vacuum,
                    Medium::Const(v) => MediumSpec::Const(*v),
                    Medium::Material(mat) => MediumSpec::Material(String::from(mat.name())),
                })
                .collect(),
        }
    }

    /// Build the system, resolving material names in `catalog`.
    ///
    /// # Errors
    /// [`SpecError::UnknownMaterial`] for a name the catalogue does not have,
    /// and [`SpecError::Shape`] when the counts do not match: there must be one
    /// gap per surface and one more medium than surfaces, and each surface
    /// needs eight asphere coefficients (and nine phase values, if given).
    pub fn to_system(&self, catalog: &MaterialDatabase) -> Result<SequentialSystem, SpecError> {
        if self.thickness.len() != self.surfaces.len() {
            return Err(SpecError::Shape("one gap per surface is required"));
        }
        if self.media.len() != self.surfaces.len() + 1 {
            return Err(SpecError::Shape(
                "one medium before and one after each surface",
            ));
        }
        let mut sys = SequentialSystem::new();
        sys.surfaces.clear();
        for s in &self.surfaces {
            let coeffs: [f64; 8] = s
                .coeffs
                .as_slice()
                .try_into()
                .map_err(|_| SpecError::Shape("eight asphere coefficients per surface"))?;
            let phase = match &s.phase {
                Some(p) => Some(
                    <[f64; 9]>::try_from(p.as_slice())
                        .map_err(|_| SpecError::Shape("nine phase coefficients"))?,
                ),
                None => None,
            };
            sys.surfaces.push(SeqSurface {
                action: match s.action {
                    ActionSpec::Refract => SurfaceAction::Refract,
                    ActionSpec::Reflect => SurfaceAction::Reflect,
                    ActionSpec::Diffract => SurfaceAction::Diffract,
                },
                asphere: Asphere {
                    curvature: s.curvature,
                    conic: s.conic,
                    coeffs,
                },
                semi_diameter: s.semi_diameter.unwrap_or(f64::INFINITY),
                phase,
            });
        }
        sys.thickness = self
            .thickness
            .iter()
            .map(|t| t.unwrap_or(f64::INFINITY))
            .collect();
        sys.media = Vec::with_capacity(self.media.len());
        for m in &self.media {
            sys.media.push(match m {
                MediumSpec::Vacuum => Medium::Vacuum,
                MediumSpec::Const(v) => Medium::Const(*v),
                MediumSpec::Material(name) => Medium::Material(
                    catalog
                        .get(name)
                        .ok_or_else(|| SpecError::UnknownMaterial(name.clone()))?,
                ),
            });
        }
        Ok(sys)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequential::TraceContext;
    use tpt_optics_core::{Kelvin, Wavelength};

    fn ctx() -> TraceContext {
        TraceContext {
            wavelength: Wavelength::from_nm(587.6),
            temperature: Kelvin::from_celsius(20.0),
        }
    }

    fn singlet() -> SequentialSystem {
        let catalog = MaterialDatabase::presets();
        let mut sys = SequentialSystem::new();
        sys.push(
            SeqSurface::aspheric(Asphere::sphere(10.0), 12.5e-3),
            5e-3,
            Medium::Material(catalog.get("bk7").unwrap()),
        );
        sys.push(
            SeqSurface::aspheric(Asphere::sphere(-10.0), 12.5e-3),
            95e-3,
            Medium::Vacuum,
        );
        sys
    }

    #[test]
    fn round_trip_keeps_the_focal_length() {
        let catalog = MaterialDatabase::presets();
        let sys = singlet();
        let spec = SystemSpec::from_system(&sys);
        let back = spec.to_system(&catalog).unwrap();
        let a = sys.paraxial(ctx()).effective_focal_length().unwrap();
        let b = back.paraxial(ctx()).effective_focal_length().unwrap();
        assert!((a - b).abs() < 1e-15, "{a} vs {b}");
    }

    #[test]
    fn materials_are_stored_by_name() {
        let spec = SystemSpec::from_system(&singlet());
        assert_eq!(spec.media[1], MediumSpec::Material(String::from("bk7")));
    }

    #[test]
    fn an_unknown_material_is_reported_by_name() {
        let mut spec = SystemSpec::from_system(&singlet());
        spec.media[1] = MediumSpec::Material(String::from("unobtainium"));
        let err = spec.to_system(&MaterialDatabase::presets()).unwrap_err();
        assert_eq!(err, SpecError::UnknownMaterial(String::from("unobtainium")));
    }

    #[test]
    fn a_mismatched_gap_count_is_a_shape_error() {
        let mut spec = SystemSpec::from_system(&singlet());
        spec.thickness.pop();
        assert!(matches!(
            spec.to_system(&MaterialDatabase::presets()),
            Err(SpecError::Shape(_))
        ));
    }

    #[test]
    fn infinite_aperture_survives_as_none() {
        let mut sys = singlet();
        sys.surfaces[0].semi_diameter = f64::INFINITY;
        let spec = SystemSpec::from_system(&sys);
        assert_eq!(spec.surfaces[0].semi_diameter, None);
        let back = spec.to_system(&MaterialDatabase::presets()).unwrap();
        assert!(back.surfaces[0].semi_diameter.is_infinite());
    }

    #[cfg(feature = "serde")]
    #[test]
    fn json_round_trip() {
        let spec = SystemSpec::from_system(&singlet());
        let text = serde_json::to_string(&spec).unwrap();
        let back: SystemSpec = serde_json::from_str(&text).unwrap();
        assert_eq!(back, spec);
    }
}
