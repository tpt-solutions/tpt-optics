//! Monte-Carlo tolerance analysis of the singlet from example 02: random
//! radius errors, then the spread of the on-axis RMS spot.
//!
//! Run with: `cargo run -p tpt-optics --example 10_tolerance_monte_carlo --features full`

use tpt_optics::core::material::MaterialDatabase;
use tpt_optics::core::{Kelvin, Wavelength};
use tpt_optics::design::tolerance::{analyze, Tolerance};
use tpt_optics::geometric::sequential::{
    rms_spot_radius, Medium, SeqSurface, SequentialSystem, TraceContext,
};
use tpt_optics::geometric::surface::Asphere;

/// Singlet from the parameter vector `[R1, R2]` (m), 5 mm thick, BK7.
fn singlet(params: &[f64]) -> SequentialSystem {
    let bk7 = MaterialDatabase::presets().get("bk7").expect("bk7 preset");
    let mut sys = SequentialSystem::new();
    sys.push(
        SeqSurface::aspheric(Asphere::sphere(1.0 / params[0]), 12.5e-3),
        5e-3,
        Medium::Material(bk7),
    );
    sys.push(
        SeqSurface::aspheric(Asphere::sphere(1.0 / params[1]), 12.5e-3),
        95e-3,
        Medium::Vacuum,
    );
    sys.push(SeqSurface::flat(12.5e-3), 0.0, Medium::Vacuum);
    sys
}

fn main() {
    let ctx = TraceContext {
        wavelength: Wavelength::from_nm(587.6),
        temperature: Kelvin::from_celsius(20.0),
    };
    let heights: Vec<f64> = (0..=20).map(|i| (i as f64 - 10.0) * 0.5e-3).collect();

    // Merit: on-axis RMS spot radius (m). Smaller is better.
    let merit = |p: &[f64]| rms_spot_radius(&singlet(p), ctx, 0.0, &heights);

    let nominal = [100e-3, -100e-3];
    println!("nominal RMS spot: {:.2} µm", merit(&nominal) * 1e6);

    // Each radius is made with a 0.1 mm standard deviation.
    let tolerances = [
        Tolerance::Normal {
            index: 0,
            sigma: 0.1e-3,
        },
        Tolerance::Normal {
            index: 1,
            sigma: 0.1e-3,
        },
    ];
    let stats = analyze(&nominal, &tolerances, None, &merit, 500, 25e-6, 2024)
        .expect("the example's tolerances name valid parameters");

    println!("over 500 random builds (radius error σ = 0.1 mm):");
    println!("  mean RMS spot      : {:.2} µm", stats.mean * 1e6);
    println!("  std deviation      : {:.2} µm", stats.sigma * 1e6);
    println!("  90th percentile    : {:.2} µm", stats.p90 * 1e6);
    println!(
        "  fails 25 µm limit  : {:.1} %",
        stats.failure_fraction * 100.0
    );
}
