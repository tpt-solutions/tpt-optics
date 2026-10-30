//! Build a [`SequentialSystem`] from a lens prescription.
//!
//! Each row describes one surface: its radius, clear semi-diameter, the
//! thickness after it, and the glass after it. Glasses are looked up by name
//! in a [`MaterialDatabase`]. A radius of zero makes a flat surface, so a
//! flat last row can serve as the image plane.

use alloc::string::String;
use alloc::vec::Vec;
use core::fmt;

use tpt_optics_core::material::MaterialDatabase;

use crate::sequential::{Medium, SeqSurface, SequentialSystem};
use crate::surface::Asphere;

/// Why a prescription could not be built.
#[derive(Clone, Debug, PartialEq)]
pub enum BuildError {
    /// The builder has no rows.
    NoSurfaces,
    /// A glass name is not in the catalog.
    UnknownGlass {
        /// The name as given.
        name: String,
        /// Closest catalog name, when one is near enough.
        suggestion: Option<String>,
    },
    /// A row has a non-finite value or a non-positive semi-diameter.
    InvalidRow {
        /// Zero-based row index.
        index: usize,
        /// What is wrong with the row.
        reason: &'static str,
    },
}

impl fmt::Display for BuildError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NoSurfaces => write!(f, "prescription has no surfaces"),
            Self::UnknownGlass { name, suggestion } => {
                write!(f, "glass '{name}' is not in the catalog")?;
                if let Some(near) = suggestion {
                    write!(f, "; did you mean '{near}'?")?;
                }
                Ok(())
            }
            Self::InvalidRow { index, reason } => write!(f, "row {index}: {reason}"),
        }
    }
}

impl core::error::Error for BuildError {}

#[derive(Clone, Debug)]
struct Row {
    radius: f64,
    semi_diameter: f64,
    thickness: f64,
    glass: String,
}

/// Builder for a sequential system, one surface at a time.
#[derive(Clone, Debug)]
pub struct SystemBuilder<'a> {
    catalog: &'a MaterialDatabase,
    rows: Vec<Row>,
}

impl<'a> SystemBuilder<'a> {
    /// Empty prescription that resolves glass names against `catalog`.
    #[must_use]
    pub fn new(catalog: &'a MaterialDatabase) -> Self {
        SystemBuilder {
            catalog,
            rows: Vec::new(),
        }
    }

    /// Appends a surface.
    ///
    /// * `radius` — radius of curvature in metres; `0.0` for a flat surface.
    /// * `semi_diameter` — clear semi-diameter in metres.
    /// * `thickness` — distance in metres to the next surface.
    /// * `glass` — medium after this surface: a catalog name, or `""` or
    ///   `"AIR"` for air.
    #[must_use]
    pub fn surface(mut self, radius: f64, semi_diameter: f64, thickness: f64, glass: &str) -> Self {
        self.rows.push(Row {
            radius,
            semi_diameter,
            thickness,
            glass: String::from(glass),
        });
        self
    }

    /// Builds the system. Glass names are checked here, so a typo is reported
    /// with a suggestion instead of being silently treated as air.
    ///
    /// # Errors
    /// [`BuildError::NoSurfaces`] for an empty prescription,
    /// [`BuildError::InvalidRow`] for a non-finite value or a non-positive
    /// semi-diameter, and [`BuildError::UnknownGlass`] for a glass that is not
    /// in the catalog.
    pub fn build(&self) -> Result<SequentialSystem, BuildError> {
        if self.rows.is_empty() {
            return Err(BuildError::NoSurfaces);
        }
        let mut sys = SequentialSystem::new();
        for (index, row) in self.rows.iter().enumerate() {
            let invalid = |reason| BuildError::InvalidRow { index, reason };
            if !(row.radius.is_finite() && row.thickness.is_finite()) {
                return Err(invalid("radius and thickness must be finite"));
            }
            if !(row.semi_diameter.is_finite() && row.semi_diameter > 0.0) {
                return Err(invalid("semi-diameter must be positive and finite"));
            }
            let medium = self.medium_for(&row.glass)?;
            let surface = if row.radius == 0.0 {
                SeqSurface::flat(row.semi_diameter)
            } else {
                SeqSurface::aspheric(Asphere::sphere(1.0 / row.radius), row.semi_diameter)
            };
            sys.push(surface, row.thickness, medium);
        }
        Ok(sys)
    }

    fn medium_for(&self, glass: &str) -> Result<Medium, BuildError> {
        if glass.is_empty() || glass.eq_ignore_ascii_case("AIR") {
            return Ok(Medium::Vacuum);
        }
        self.catalog
            .get(glass)
            .map(Medium::Material)
            .ok_or_else(|| BuildError::UnknownGlass {
                name: String::from(glass),
                suggestion: self.catalog.closest_name(glass, 2).map(String::from),
            })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn builds_a_two_surface_system() {
        let catalog = MaterialDatabase::presets();
        let sys = SystemBuilder::new(&catalog)
            .surface(100e-3, 12.5e-3, 5e-3, "bk7")
            .surface(-100e-3, 12.5e-3, 95e-3, "")
            .build()
            .unwrap();
        assert_eq!(sys.surfaces.len(), 2);
    }

    #[test]
    fn unknown_glass_is_reported_with_a_suggestion() {
        let catalog = MaterialDatabase::presets();
        let err = SystemBuilder::new(&catalog)
            .surface(100e-3, 12.5e-3, 5e-3, "bk8")
            .build()
            .unwrap_err();
        assert_eq!(
            err,
            BuildError::UnknownGlass {
                name: String::from("bk8"),
                suggestion: Some(String::from("bk7")),
            }
        );
    }

    #[test]
    fn empty_prescription_is_an_error() {
        let catalog = MaterialDatabase::presets();
        assert_eq!(
            SystemBuilder::new(&catalog).build().unwrap_err(),
            BuildError::NoSurfaces
        );
    }

    #[test]
    fn bad_rows_are_rejected_with_their_index() {
        let catalog = MaterialDatabase::presets();
        let err = SystemBuilder::new(&catalog)
            .surface(100e-3, 12.5e-3, 5e-3, "bk7")
            .surface(100e-3, 0.0, 5e-3, "")
            .build()
            .unwrap_err();
        assert!(
            matches!(err, BuildError::InvalidRow { index: 1, .. }),
            "{err:?}"
        );
    }
}
