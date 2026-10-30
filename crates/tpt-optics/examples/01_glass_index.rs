//! Refractive index of N-BK7 from the Sellmeier model, and its Abbe number.
//!
//! Run with: `cargo run -p tpt-optics --example 01_glass_index --features full`

use tpt_optics::core::dispersion::{Dispersion, Sellmeier};
use tpt_optics::core::Wavelength;

fn main() {
    let bk7 = Sellmeier::bk7();

    let n_d = bk7.real_index(Wavelength::from_nm(587.6));
    let n_f = bk7.real_index(Wavelength::from_nm(486.1));
    let n_c = bk7.real_index(Wavelength::from_nm(656.3));
    let n_he_ne = bk7.real_index(Wavelength::from_nm(632.8));

    println!("N-BK7 refractive index (Sellmeier):");
    println!("  n_d  (587.6 nm) = {n_d:.5}   published 1.5168");
    println!("  n_F  (486.1 nm) = {n_f:.5}   published 1.52238");
    println!("  n_C  (656.3 nm) = {n_c:.5}   published 1.51432");
    println!("  He-Ne (632.8 nm) = {n_he_ne:.5}");

    let abbe = (n_d - 1.0) / (n_f - n_c);
    println!("Abbe number V_d = {abbe:.2}   published 64.17");
}
