//! Guided modes of a silicon slab and of a step-index single-mode fiber.
//!
//! The fiber parameters are illustrative values for a standard single-mode
//! fiber at 1550 nm. Check them against your fiber's datasheet.
//!
//! Run with: `cargo run -p tpt-optics --example 05_slab_fiber_modes --features full`

use tpt_optics::core::Wavelength;
use tpt_optics::waveguide::fiber::fiber_lp_modes;
use tpt_optics::waveguide::slab::{slab_modes, SlabLayer, SlabPolarization};

fn main() {
    let wl = Wavelength::from_nm(1550.0);

    // 220 nm silicon slab in silica.
    let slab = SlabLayer {
        n_core: 3.48,
        n_clad: 1.444,
        thickness: 220e-9,
    };
    println!("Silicon slab (220 nm) at 1550 nm, TE modes:");
    for mode in slab_modes(&slab, wl, SlabPolarization::Te) {
        println!("  order {}: n_eff = {:.4}", mode.order, mode.n_eff);
    }

    // Step-index fiber: 8.2 µm core diameter, n_core − n_clad ≈ 0.0057.
    let modes = fiber_lp_modes(4.1e-6, 1.4504, 1.4447, wl, 2);
    println!("\nStep-index fiber (core radius 4.1 µm) at 1550 nm:");
    for m in &modes {
        println!("  LP{}{}: n_eff = {:.5}", m.l, m.m, m.n_eff);
    }
}
