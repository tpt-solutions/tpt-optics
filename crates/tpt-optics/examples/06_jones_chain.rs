//! Polarisation through a chain of optical elements, with Jones calculus.
//!
//! Light polarised along x passes a polariser at 45°, then a quarter-wave
//! plate with its fast axis at 45°. The example prints the intensity after
//! each element and the final polarisation state.
//!
//! Run with: `cargo run -p tpt-optics --example 06_jones_chain --features full`

use core::f64::consts::{FRAC_PI_2, FRAC_PI_4};

use tpt_optics::polarization::jones::{JonesMatrix, JonesVector};

fn main() {
    let input = JonesVector::linear_x();
    println!(
        "input                       intensity {:.4}",
        input.intensity()
    );

    let after_polariser = input.apply(&JonesMatrix::polarizer(FRAC_PI_4));
    println!(
        "after polariser at 45°      intensity {:.4}",
        after_polariser.intensity()
    );

    let after_waveplate = after_polariser.apply(&JonesMatrix::waveplate(FRAC_PI_2, FRAC_PI_4));
    println!(
        "after λ/4 plate at 45°      intensity {:.4}",
        after_waveplate.intensity()
    );
    println!("\nfinal Jones vector: {after_waveplate:?}");
}
