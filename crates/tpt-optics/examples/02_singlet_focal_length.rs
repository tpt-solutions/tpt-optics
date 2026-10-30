//! Effective focal length of a biconvex N-BK7 singlet, traced and compared
//! with the thick-lens formula.
//!
//! Run with: `cargo run -p tpt-optics --example 02_singlet_focal_length --features full`

use tpt_optics::core::index::IndexAt;
use tpt_optics::core::material::MaterialDatabase;
use tpt_optics::core::{Kelvin, Wavelength};
use tpt_optics::geometric::sequential::{Medium, SeqSurface, SequentialSystem, TraceContext};
use tpt_optics::geometric::surface::Asphere;

/// Thickness of the lens centre (m).
const THICKNESS: f64 = 5e-3;
/// Clear semi-diameter of both surfaces (m).
const SEMI_DIAMETER: f64 = 12.5e-3;

/// Biconvex singlet with radii `r1`, `r2` (m) in N-BK7, air on both sides.
fn singlet(r1: f64, r2: f64) -> SequentialSystem {
    let bk7 = MaterialDatabase::presets().get("bk7").expect("bk7 preset");
    let mut sys = SequentialSystem::new();
    sys.push(
        SeqSurface::aspheric(Asphere::sphere(1.0 / r1), SEMI_DIAMETER),
        THICKNESS,
        Medium::Material(bk7),
    );
    sys.push(
        SeqSurface::aspheric(Asphere::sphere(1.0 / r2), SEMI_DIAMETER),
        95e-3,
        Medium::Vacuum,
    );
    sys
}

fn main() {
    let (r1, r2) = (100e-3, -100e-3);
    let sys = singlet(r1, r2);
    let wl = Wavelength::from_nm(587.6);
    let temp = Kelvin::from_celsius(20.0);
    let ctx = TraceContext {
        wavelength: wl,
        temperature: temp,
    };

    let n = MaterialDatabase::presets()
        .get("bk7")
        .expect("bk7 preset")
        .index_at(wl, temp)
        .n;
    // Thick-lens formula: 1/f = (n−1)[1/R1 − 1/R2 + (n−1)t/(n R1 R2)].
    let power = (n - 1.0) * (1.0 / r1 - 1.0 / r2 + (n - 1.0) * THICKNESS / (n * r1 * r2));
    let f_formula = 1.0 / power;

    let f_traced = sys
        .paraxial(ctx)
        .effective_focal_length()
        .expect("singlet has finite focal length");

    println!(
        "Biconvex N-BK7 singlet, R1 = {:.0} mm, R2 = {:.0} mm, t = {:.0} mm",
        r1 * 1e3,
        r2 * 1e3,
        THICKNESS * 1e3
    );
    println!("  index at 587.6 nm        : {n:.5}");
    println!("  thick-lens focal length  : {:.3} mm", f_formula * 1e3);
    println!("  traced focal length (EFL): {:.3} mm", f_traced * 1e3);
}
