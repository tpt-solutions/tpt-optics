//! The design file: a TOML description of a lens, the quantities to report,
//! and the operands and variables for an optimisation.
//!
//! ```toml
//! name = "singlet"
//! stop = 0
//! fields_deg = [0.0, 2.0]
//! heights_mm = [-5.0, -2.5, 0.0, 2.5, 5.0]
//!
//! [context]
//! wavelength_nm = 587.6
//! temperature_c = 20.0
//!
//! [[surfaces]]
//! radius_mm = 100.0
//! semi_diameter_mm = 12.5
//! thickness_mm = 5.0
//! glass = "N-BK7"
//!
//! [[surfaces]]
//! radius_mm = -100.0
//! semi_diameter_mm = 12.5
//! thickness_mm = 95.0
//!
//! [[surfaces]]            # flat image plane
//! radius_mm = 0.0
//! semi_diameter_mm = 12.5
//! thickness_mm = 0.0
//!
//! [[operands]]
//! kind = "efl"
//! target_mm = 100.0
//!
//! [optimize]
//! variables = [0, 1]
//! iterations = 50
//! ```

use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use tpt_optics::core::material::MaterialDatabase;
use tpt_optics::core::{Kelvin, Wavelength};
use tpt_optics::design::analysis::{analyze, AnalysisReport};
use tpt_optics::design::operands::{optimize_system, Operand};
use tpt_optics::design::optimize::{OptimizeResult, StopReason};
use tpt_optics::geometric::builder::SystemBuilder;
use tpt_optics::geometric::sequential::{SequentialSystem, TraceContext};
use tpt_optics::io::glass::reference_catalog;

/// Wavelength and temperature of the analysis.
#[derive(Clone, Debug, Serialize, Deserialize, JsonSchema, PartialEq)]
pub struct Context {
    /// Vacuum wavelength (nm).
    #[serde(default = "default_wavelength")]
    pub wavelength_nm: f64,
    /// Temperature (°C).
    #[serde(default = "default_temperature")]
    pub temperature_c: f64,
}

fn default_wavelength() -> f64 {
    587.6
}

fn default_temperature() -> f64 {
    20.0
}

impl Default for Context {
    fn default() -> Self {
        Context {
            wavelength_nm: default_wavelength(),
            temperature_c: default_temperature(),
        }
    }
}

/// One surface row of the prescription.
#[derive(Clone, Debug, Serialize, Deserialize, JsonSchema, PartialEq)]
pub struct SurfaceSpec {
    /// Radius of curvature (mm); zero is flat.
    pub radius_mm: f64,
    /// Clear semi-diameter (mm).
    pub semi_diameter_mm: f64,
    /// Thickness to the next surface (mm).
    pub thickness_mm: f64,
    /// Glass after the surface, or empty for air.
    #[serde(default)]
    pub glass: String,
}

/// One requirement. `kind` selects the operand; the other fields apply to
/// the kinds that use them.
#[derive(Clone, Debug, Serialize, Deserialize, JsonSchema, PartialEq)]
pub struct OperandSpec {
    /// `efl`, `bfl`, `rms_spot`, `thickness`, `ray_height`, or `edge_thickness`.
    pub kind: String,
    /// Target (mm), or for spots the target spot radius (µm).
    #[serde(default)]
    pub target: f64,
    /// Lower bound (mm) for range kinds.
    #[serde(default)]
    pub min_mm: f64,
    /// Upper bound (mm) for range kinds.
    #[serde(default)]
    pub max_mm: f64,
    /// Weight of the residual.
    #[serde(default = "default_weight")]
    pub weight: f64,
    /// Gap index for `thickness` and `edge_thickness`.
    #[serde(default)]
    pub index: usize,
    /// Surface index for `ray_height`.
    #[serde(default)]
    pub surface: usize,
    /// Field angle (degrees) for `rms_spot` and `ray_height`.
    #[serde(default)]
    pub field_deg: f64,
    /// Pupil height (mm) for `ray_height`.
    #[serde(default)]
    pub pupil_mm: f64,
}

fn default_weight() -> f64 {
    1.0
}

/// Optimisation settings: which surface radii vary, and for how long.
#[derive(Clone, Debug, Serialize, Deserialize, JsonSchema, PartialEq)]
pub struct OptimizeSpec {
    /// Indices of the surfaces whose radius is a variable.
    pub variables: Vec<usize>,
    /// Levenberg–Marquardt iterations.
    #[serde(default = "default_iterations")]
    pub iterations: usize,
}

fn default_iterations() -> usize {
    50
}

/// A whole design file.
#[derive(Clone, Debug, Serialize, Deserialize, JsonSchema, PartialEq)]
pub struct Design {
    /// Name of the design.
    #[serde(default)]
    pub name: String,
    /// Analysis wavelength and temperature.
    #[serde(default)]
    pub context: Context,
    /// Surfaces, object side first. The last row is the image plane.
    pub surfaces: Vec<SurfaceSpec>,
    /// Index of the aperture stop surface.
    #[serde(default)]
    pub stop: usize,
    /// Field angles (degrees) to analyse.
    #[serde(default)]
    pub fields_deg: Vec<f64>,
    /// Pupil heights (mm) for the fans and spots.
    #[serde(default)]
    pub heights_mm: Vec<f64>,
    /// Requirements for optimisation.
    #[serde(default)]
    pub operands: Vec<OperandSpec>,
    /// Optimisation settings, if the design is to be optimised.
    #[serde(default)]
    pub optimize: Option<OptimizeSpec>,
}

/// Why a design could not be used.
#[derive(Clone, Debug, PartialEq)]
pub struct DesignError(pub String);

impl core::fmt::Display for DesignError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(&self.0)
    }
}

impl Design {
    /// Parse a design from TOML text.
    ///
    /// # Errors
    /// Returns the TOML error text when the file does not match the schema.
    pub fn from_toml(text: &str) -> Result<Self, DesignError> {
        toml::from_str(text).map_err(|e| DesignError(e.to_string()))
    }

    /// Serialise to TOML text.
    ///
    /// # Errors
    /// Returns the serialisation error text.
    pub fn to_toml(&self) -> Result<String, DesignError> {
        toml::to_string_pretty(self).map_err(|e| DesignError(e.to_string()))
    }

    /// Trace context for the analysis wavelength.
    #[must_use]
    pub fn context(&self) -> TraceContext {
        TraceContext {
            wavelength: Wavelength::from_nm(self.context.wavelength_nm),
            temperature: Kelvin::from_celsius(self.context.temperature_c),
        }
    }

    /// Glass catalogue used to resolve names.
    #[must_use]
    pub fn catalog() -> MaterialDatabase {
        reference_catalog()
    }

    /// Build the system with the radii `radii_mm` (one per surface).
    ///
    /// # Errors
    /// Returns the builder's error text for an unknown glass or a bad row.
    pub fn build_with_radii(&self, radii_mm: &[f64]) -> Result<SequentialSystem, DesignError> {
        let catalog = Self::catalog();
        let mut builder = SystemBuilder::new(&catalog);
        for (spec, &r) in self.surfaces.iter().zip(radii_mm) {
            builder = builder.surface(
                r * 1e-3,
                spec.semi_diameter_mm * 1e-3,
                spec.thickness_mm * 1e-3,
                &spec.glass,
            );
        }
        builder.build().map_err(|e| DesignError(e.to_string()))
    }

    /// Build with the radii in the file.
    ///
    /// # Errors
    /// As [`Design::build_with_radii`].
    pub fn build(&self) -> Result<SequentialSystem, DesignError> {
        let radii: Vec<f64> = self.surfaces.iter().map(|s| s.radius_mm).collect();
        self.build_with_radii(&radii)
    }

    /// Heights (m) of the pupil fan.
    #[must_use]
    pub fn heights_m(&self) -> Vec<f64> {
        self.heights_mm.iter().map(|h| h * 1e-3).collect()
    }

    /// The operands, converted to the design crate's form.
    ///
    /// # Errors
    /// Returns an error for an unknown `kind`.
    pub fn operands(&self) -> Result<Vec<Operand>, DesignError> {
        self.operands
            .iter()
            .map(|spec| operand_of(spec, &self.heights_m()))
            .collect()
    }

    /// Analyse the design as it is in the file.
    ///
    /// # Errors
    /// As [`Design::build`].
    pub fn analyze(&self) -> Result<AnalysisReport, DesignError> {
        let sys = self.build()?;
        Ok(analyze(
            &sys,
            self.context(),
            self.stop,
            &self.heights_m(),
            &self.fields_deg,
            &[],
        ))
    }

    /// Optimise the variables named in `optimize`, starting from the file's
    /// radii. Returns the radii found (mm) and the optimiser's result.
    ///
    /// # Errors
    /// Returns an error when there is no `[optimize]` table, or for an unknown
    /// operand kind.
    pub fn optimize(&self) -> Result<(Vec<f64>, OptimizeResult), DesignError> {
        let spec = self
            .optimize
            .as_ref()
            .ok_or_else(|| DesignError(String::from("the design has no [optimize] table")))?;
        let ops = self.operands()?;
        let start: Vec<f64> = spec
            .variables
            .iter()
            .map(|&i| self.surfaces[i].radius_mm)
            .collect();
        let base: Vec<f64> = self.surfaces.iter().map(|s| s.radius_mm).collect();
        let build = |x: &[f64]| {
            let mut radii = base.clone();
            for (&i, &v) in spec.variables.iter().zip(x) {
                radii[i] = v;
            }
            self.build_with_radii(&radii).ok()
        };
        let result = optimize_system(&build, &start, &ops, self.context(), spec.iterations);
        let mut radii = base;
        for (&i, &v) in spec.variables.iter().zip(&result.x) {
            radii[i] = v;
        }
        Ok((radii, result))
    }
}

/// Text form of a stop reason, for reports.
#[must_use]
pub fn stop_text(stop: StopReason) -> &'static str {
    match stop {
        StopReason::Converged => "converged",
        StopReason::MaxIterations => "max iterations",
        StopReason::Stalled => "stalled",
        StopReason::Cancelled => "cancelled",
    }
}

fn operand_of(spec: &OperandSpec, heights_m: &[f64]) -> Result<Operand, DesignError> {
    let mm = |v: f64| v * 1e-3;
    Ok(match spec.kind.as_str() {
        "efl" => Operand::Efl {
            target: mm(spec.target),
            weight: spec.weight,
        },
        "bfl" => Operand::Bfl {
            target: mm(spec.target),
            weight: spec.weight,
        },
        "rms_spot" => Operand::RmsSpot {
            field_deg: spec.field_deg,
            heights: heights_m.to_vec(),
            target: spec.target * 1e-6,
            weight: spec.weight,
        },
        "thickness" => Operand::ThicknessRange {
            index: spec.index,
            min: mm(spec.min_mm),
            max: mm(spec.max_mm),
            weight: spec.weight,
        },
        "edge_thickness" => Operand::EdgeThickness {
            index: spec.index,
            min: mm(spec.min_mm),
            max: mm(spec.max_mm),
            weight: spec.weight,
        },
        "ray_height" => Operand::RayHeight {
            surface: spec.surface,
            field_deg: spec.field_deg,
            pupil_height: mm(spec.pupil_mm),
            target: mm(spec.target),
            weight: spec.weight,
        },
        other => return Err(DesignError(format!("unknown operand kind '{other}'"))),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const SINGLET: &str = r#"
name = "singlet"
fields_deg = [0.0]
heights_mm = [-5.0, 0.0, 5.0]

[[surfaces]]
radius_mm = 100.0
semi_diameter_mm = 12.5
thickness_mm = 5.0
glass = "N-BK7"

[[surfaces]]
radius_mm = -100.0
semi_diameter_mm = 12.5
thickness_mm = 95.0

[[surfaces]]
radius_mm = 0.0
semi_diameter_mm = 12.5
thickness_mm = 0.0

[[operands]]
kind = "efl"
target = 100.0

[optimize]
variables = [0]
iterations = 40
"#;

    #[test]
    fn parses_and_builds_the_singlet() {
        let d = Design::from_toml(SINGLET).unwrap();
        assert_eq!(d.name, "singlet");
        assert_eq!(d.context, Context::default());
        let sys = d.build().unwrap();
        assert_eq!(sys.surfaces.len(), 3);
    }

    #[test]
    fn round_trips_through_toml() {
        let d = Design::from_toml(SINGLET).unwrap();
        let text = d.to_toml().unwrap();
        assert_eq!(Design::from_toml(&text).unwrap(), d);
    }

    #[test]
    fn unknown_glass_is_a_build_error() {
        let bad = SINGLET.replace("N-BK7", "N-BK9");
        let d = Design::from_toml(&bad).unwrap();
        assert!(d.build().is_err());
    }

    #[test]
    fn unknown_operand_kind_is_rejected() {
        let bad = SINGLET.replace("kind = \"efl\"", "kind = \"magic\"");
        let d = Design::from_toml(&bad).unwrap();
        assert!(d.operands().is_err());
    }

    #[test]
    fn optimisation_moves_the_focal_length_to_the_target() {
        let d = Design::from_toml(SINGLET).unwrap();
        let (radii, result) = d.optimize().unwrap();
        assert!(matches!(
            result.stop,
            StopReason::Converged | StopReason::MaxIterations
        ));
        let sys = d.build_with_radii(&radii).unwrap();
        let efl = sys.paraxial(d.context()).effective_focal_length().unwrap();
        assert!((efl - 0.1).abs() < 1e-3, "EFL {efl}");
    }

    #[test]
    fn a_design_without_optimize_table_says_so() {
        let text = SINGLET.split("[optimize]").next().unwrap();
        let d = Design::from_toml(text).unwrap();
        assert!(d.optimize().is_err());
    }
}
