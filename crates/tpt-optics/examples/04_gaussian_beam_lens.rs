//! A Gaussian beam focused by a thin lens, using the complex-beam (q) and
//! ABCD formalism.
//!
//! Run with: `cargo run -p tpt-optics --example 04_gaussian_beam_lens --features full`

use tpt_optics::core::linalg::Matrix2;
use tpt_optics::core::Wavelength;
use tpt_optics::wave::gaussian::GaussianBeam;

fn main() {
    let wl = Wavelength::from_nm(633.0);
    // A 0.5 mm waist radius in air.
    let input = GaussianBeam::at_waist(0.5e-3, wl, 1.0);
    println!(
        "Input beam: waist w0 = 0.500 mm, z_R = {:.1} mm",
        input.rayleigh_range() * 1e3
    );

    // Thin lens of focal length 100 mm, placed directly at the waist.
    let f = 100e-3;
    let lens = Matrix2 {
        a: 1.0,
        b: 0.0,
        c: -1.0 / f,
        d: 1.0,
    };
    let after_lens = input.through(&lens);

    // Propagate the focused beam and report the spot size along the axis.
    println!("\nDistance after lens (mm)   spot radius w (µm)");
    for z_mm in [0.0, 50.0, 90.0, 100.0, 110.0, 150.0] {
        let beam = after_lens.propagated(z_mm * 1e-3);
        println!("  {z_mm:>8.1}                {:>10.2}", beam.w() * 1e6);
    }
    println!(
        "\nWaist behind the lens: z_R of the focused beam = {:.2} mm",
        after_lens.rayleigh_range() * 1e3
    );
}
