//! RMS spot radius of the singlet from example 02 at on-axis and off-axis
//! field angles. This is a singlet, not a Cooke triplet: its aberrations are
//! larger, which is what the spot sizes show.
//!
//! Run with: `cargo run -p tpt-optics --example 03_lens_spot_diagram --features full`

use tpt_optics::core::material::MaterialDatabase;
use tpt_optics::core::{Kelvin, Wavelength};
use tpt_optics::geometric::sequential::{
    rms_spot_radius, Medium, SeqSurface, SequentialSystem, TraceContext,
};
use tpt_optics::geometric::surface::Asphere;

fn singlet(r1: f64, r2: f64) -> SequentialSystem {
    let bk7 = MaterialDatabase::presets().get("bk7").expect("bk7 preset");
    let mut sys = SequentialSystem::new();
    sys.push(
        SeqSurface::aspheric(Asphere::sphere(1.0 / r1), 12.5e-3),
        5e-3,
        Medium::Material(bk7),
    );
    sys.push(
        SeqSurface::aspheric(Asphere::sphere(1.0 / r2), 12.5e-3),
        95e-3,
        Medium::Vacuum,
    );
    sys.push(SeqSurface::flat(12.5e-3), 0.0, Medium::Vacuum);
    sys
}

fn main() {
    let sys = singlet(100e-3, -100e-3);
    let ctx = TraceContext {
        wavelength: Wavelength::from_nm(587.6),
        temperature: Kelvin::from_celsius(20.0),
    };

    // Pupil fan: 21 rays across the 10 mm aperture radius.
    let heights: Vec<f64> = (0..=20).map(|i| (i as f64 - 10.0) * 0.5e-3).collect();

    println!("RMS spot radius of the singlet (21-ray fan):");
    for angle in [0.0, 2.0, 5.0] {
        let rms = rms_spot_radius(&sys, ctx, angle, &heights);
        println!("  field {angle:>3.1}°  RMS spot = {:>8.2} µm", rms * 1e6);
    }
}
