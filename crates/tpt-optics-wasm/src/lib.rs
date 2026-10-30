//! WebAssembly bindings for tpt-optics: glass lookup and a singlet report,
//! for the browser demo in `www/` and for Node.
//!
//! Build with `wasm-pack build crates/tpt-optics-wasm --target web` (or
//! `--target nodejs`). Lengths cross the boundary in millimetres and the
//! report comes back as JSON text.

use tpt_optics::core::index::IndexAt;
use tpt_optics::core::{Kelvin, Wavelength};
use tpt_optics::design::analysis::analyze;
use tpt_optics::geometric::presets::singlet;
use tpt_optics::geometric::sequential::{SequentialSystem, TraceContext};
use wasm_bindgen::prelude::*;

/// Refractive index of a catalogue glass at a wavelength (nm), at 20 °C.
///
/// # Errors
/// Returns an error for a glass name that is not in the catalogue.
#[wasm_bindgen]
pub fn glass_index(name: &str, wavelength_nm: f64) -> Result<f64, JsValue> {
    let catalog = tpt_optics::io::glass::reference_catalog();
    let material = catalog
        .get(name)
        .ok_or_else(|| JsValue::from_str(&format!("glass '{name}' is not in the catalogue")))?;
    Ok(material
        .index_at(
            Wavelength::from_nm(wavelength_nm),
            Kelvin::from_celsius(20.0),
        )
        .n)
}

/// Biconvex or meniscus singlet in `glass`, with its image plane at the
/// paraxial focus. Returns the EFL, BFL and an on-axis analysis as JSON.
///
/// Radii are in millimetres, with zero for a flat surface. The semi-diameter
/// is in millimetres.
///
/// # Errors
/// Returns an error for an unknown glass or an invalid prescription.
#[wasm_bindgen]
pub fn singlet_report(
    r1_mm: f64,
    r2_mm: f64,
    thickness_mm: f64,
    semi_diameter_mm: f64,
    glass: &str,
) -> Result<String, JsValue> {
    let catalog = tpt_optics::io::glass::reference_catalog();
    let sd = semi_diameter_mm * 1e-3;
    let t = thickness_mm * 1e-3;
    // Provisional image distance, then move the image plane to the focus.
    let provisional = singlet(&catalog, r1_mm * 1e-3, r2_mm * 1e-3, sd, t, glass, 95e-3)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    let ctx = TraceContext {
        wavelength: Wavelength::from_nm(587.6),
        temperature: Kelvin::from_celsius(20.0),
    };
    let bfl = provisional
        .paraxial(ctx)
        .back_focal_length()
        .unwrap_or(95e-3);
    // The BFL is measured from the provisional image plane, so the new plane is
    // that far beyond it.
    let image_distance = 95e-3 + bfl;
    let sys: SequentialSystem = singlet(
        &catalog,
        r1_mm * 1e-3,
        r2_mm * 1e-3,
        sd,
        t,
        glass,
        image_distance,
    )
    .map_err(|e| JsValue::from_str(&e.to_string()))?;
    let heights: Vec<f64> = (0..=10).map(|i| (i as f64 - 5.0) / 5.0 * sd).collect();
    let report = analyze(&sys, ctx, 0, &heights, &[0.0], &[]);
    Ok(report.to_json())
}
