//! Python bindings for tpt-optics (PyO3).
//!
//! The module is `tpt_optics`. Lengths are in millimetres at this boundary, and
//! analysis results come back as JSON text, which `json.loads` turns into
//! dictionaries. Build a wheel with `maturin build` in this crate, or a
//! development build with `cargo build -p tpt-optics-py`.

// The PyO3 macros expand to conversions that clippy reports as no-ops.
#![allow(clippy::useless_conversion)]

use pyo3::exceptions::{PyKeyError, PyValueError};
use pyo3::prelude::*;
use tpt_optics::core::index::IndexAt;
use tpt_optics::core::{Kelvin, Wavelength};
use tpt_optics::design::analysis::analyze;
use tpt_optics::design::glassmap::{rank_achromats, MapSpec};
use tpt_optics::geometric::presets::singlet;
use tpt_optics::geometric::sequential::TraceContext;

fn ctx(wavelength_nm: f64) -> TraceContext {
    TraceContext {
        wavelength: Wavelength::from_nm(wavelength_nm),
        temperature: Kelvin::from_celsius(20.0),
    }
}

/// Refractive index of a catalogue glass at `wavelength_nm`, at 20 °C.
#[pyfunction]
fn glass_index(name: &str, wavelength_nm: f64) -> PyResult<f64> {
    let catalog = tpt_optics::io::glass::reference_catalog();
    let material = catalog
        .get(name)
        .ok_or_else(|| PyKeyError::new_err(format!("glass '{name}' is not in the catalogue")))?;
    Ok(material
        .index_at(
            Wavelength::from_nm(wavelength_nm),
            Kelvin::from_celsius(20.0),
        )
        .n)
}

/// Singlet with radii `r1_mm`, `r2_mm` (zero is flat), centre thickness
/// `thickness_mm`, and semi-diameter `semi_diameter_mm`. The image plane is at
/// the paraxial focus. Returns the analysis as a JSON string.
#[pyfunction]
fn singlet_report(
    r1_mm: f64,
    r2_mm: f64,
    thickness_mm: f64,
    semi_diameter_mm: f64,
    glass: &str,
) -> PyResult<String> {
    let catalog = tpt_optics::io::glass::reference_catalog();
    let sd = semi_diameter_mm * 1e-3;
    let t = thickness_mm * 1e-3;
    let build = |image: f64| {
        singlet(&catalog, r1_mm * 1e-3, r2_mm * 1e-3, sd, t, glass, image)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    };
    let provisional = build(95e-3)?;
    let bfl = provisional
        .paraxial(ctx(587.6))
        .back_focal_length()
        .unwrap_or(95e-3);
    let sys = build(95e-3 + bfl)?;
    let heights: Vec<f64> = (0..=10).map(|i| (i as f64 - 5.0) / 5.0 * sd).collect();
    Ok(analyze(&sys, ctx(587.6), 0, &heights, &[0.0], &[]).to_json())
}

/// Rank crown–flint pairs for an achromat of focal length `focal_mm`. Returns a
/// list of `(crown, flint, rms_spot_um, axial_colour_um)` tuples, best first.
// PyO3 extracts owned argument types, so the lint's by-reference advice does not apply.
#[allow(clippy::needless_pass_by_value)]
#[pyfunction]
fn rank_achromat_pairs(
    focal_mm: f64,
    crowns: Vec<String>,
    flints: Vec<String>,
) -> Vec<(String, String, f64, f64)> {
    let catalog = tpt_optics::io::glass::reference_catalog();
    let spec = MapSpec {
        focal_length: focal_mm * 1e-3,
        semi_diameter: 6e-3,
        crown_thickness: 4e-3,
        flint_thickness: 2e-3,
        image_distance: 95e-3,
    };
    let crown_refs: Vec<&str> = crowns.iter().map(String::as_str).collect();
    let flint_refs: Vec<&str> = flints.iter().map(String::as_str).collect();
    rank_achromats(&catalog, spec, &crown_refs, &flint_refs)
        .into_iter()
        .map(|p| (p.crown, p.flint, p.rms_spot * 1e6, p.axial_colour * 1e6))
        .collect()
}

/// The `tpt_optics` Python module.
#[pymodule]
#[pyo3(name = "tpt_optics")]
fn tpt_optics_module(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(glass_index, m)?)?;
    m.add_function(wrap_pyfunction!(singlet_report, m)?)?;
    m.add_function(wrap_pyfunction!(rank_achromat_pairs, m)?)?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
