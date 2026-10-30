//! Band gap and point-defect level of a square lattice of dielectric rods
//! (E-polarisation), from the supercell plane-wave expansion.
//!
//! Run with: `cargo run -p tpt-optics --example 09_photonic_crystal_gap --features full`

use tpt_optics::waveguide::supercell::{
    bulk_gap_from_supercell, point_defect_rods, supercell_tm_freqs,
};

fn main() {
    // Square rods, ε = 8.9, r/a = 0.2, in a 5×5 supercell with a rod at the centre.
    let (m, eps, r_over_a, n_g) = (5, 8.9, 0.2, 17);

    let gap = bulk_gap_from_supercell(m, eps, r_over_a, 17, 8, (0.2, 0.6), true)
        .expect("valid supercell")
        .expect("TM gap inside the window");
    println!(
        "E-polarised full-zone gap: ({:.3}, {:.3}) c/a   (earlier estimate 0.324–0.447)",
        gap.0, gap.1
    );

    // Remove the centre rod and look for a level inside the gap.
    let missing = point_defect_rods(m, eps, r_over_a, 1.0, 0.0);
    let levels = supercell_tm_freqs(m, &missing, n_g, (0.0, 0.0)).expect("valid supercell");
    let defect: Vec<f64> = levels
        .into_iter()
        .filter(|f| (0.36..0.42).contains(f))
        .collect();
    println!("missing-rod defect levels in (0.36, 0.42) c/a at Γ: {defect:.4?}");
}
