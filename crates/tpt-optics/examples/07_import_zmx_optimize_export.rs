//! Import a Zemax `.zmx` lens, optimise its front radius to a target focal
//! length, and export the result back to `.zmx` text.
//!
//! Run with: `cargo run -p tpt-optics --example 07_import_zmx_optimize_export --features full`

use tpt_optics::core::{Kelvin, Wavelength};
use tpt_optics::design::optimize::damped_least_squares;
use tpt_optics::geometric::sequential::TraceContext;
use tpt_optics::io::glass::reference_catalog;
use tpt_optics::io::zmx;

const SOURCE: &str = "\
MODE SEQ
NAME singlet example
SURF
CURV 0 0 0 0 0 0 0 0 0
DISZ INFINITY
SURF
STOP
DIAM 12.5
CURV 0.01 0 0 0 0 0 0 0 0
GLAS N-BK7 0 0 0 0
CONI 0.0
DISZ 5.0
SURF
DIAM 12.0
CURV -0.01 0 0 0 0 0 0 0 0
GLAS AIR
DISZ 95.0
SURF
CURV 0 0 0 0 0 0 0 0 0
DISZ 0
";

/// Target effective focal length (m).
const TARGET_EFL: f64 = 100e-3;

fn main() {
    let doc = zmx::parse(SOURCE).expect("valid zmx text");
    let catalog = reference_catalog();
    let wl = Wavelength::from_nm(587.6);
    let temp = Kelvin::from_celsius(20.0);
    let ctx = TraceContext {
        wavelength: wl,
        temperature: temp,
    };

    // Residual: focal length error, with the front radius's curvature (1/m)
    // as the single free parameter.
    let residual = |x: &[f64], out: &mut Vec<f64>| {
        let mut trial = doc.clone();
        trial.surfaces[1].curvature = x[0];
        let efl = zmx::to_sequential(&trial, &catalog, wl, temp)
            .ok()
            .and_then(|sys| sys.paraxial(ctx).effective_focal_length())
            .unwrap_or(0.0);
        out.clear();
        out.push(efl - TARGET_EFL);
    };

    let start = doc.surfaces[1].curvature;
    println!("start:  front curvature {start:.3} 1/m");

    let result = damped_least_squares(&residual, &[start], 1, None, 50);
    println!(
        "result: front curvature {:.3} 1/m (R = {:.2} mm), |error| = {:.2e} m after {} evaluations",
        result.x[0],
        1e3 / result.x[0],
        result.f.abs(),
        result.evaluations
    );

    let mut optimised = doc.clone();
    optimised.surfaces[1].curvature = result.x[0];
    println!("\nexported .zmx:\n{}", zmx::write(&optimised));
}
