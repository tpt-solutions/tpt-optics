//! `verify_system` on sequential systems beyond the singlet: a mirror and a
//! glass plate. Every invariant must hold for each.

use tpt_optics_core::{Kelvin, Wavelength};
use tpt_optics_geometric::sequential::{Medium, SeqSurface, SequentialSystem, TraceContext};
use tpt_optics_verify::system_check::verify_system;

fn ctx() -> TraceContext {
    TraceContext {
        wavelength: Wavelength::from_nm(587.6),
        temperature: Kelvin::from_celsius(20.0),
    }
}

fn assert_all_pass(sys: &SequentialSystem) {
    let checks = verify_system(sys, ctx());
    assert_ne!(checks.len(), 0);
    for c in &checks {
        assert!(c.passed, "{}: {}", c.name, c.detail);
    }
}

#[test]
fn concave_mirror_passes_every_invariant() {
    // A spherical mirror reflects, so the reversed trace must retrace the ray.
    let mut sys = SequentialSystem::new();
    sys.push(SeqSurface::mirror(1.0 / 0.5, 20e-3), 0.0, Medium::Vacuum);
    assert_all_pass(&sys);
}

#[test]
fn glass_plate_passes_every_invariant() {
    // A plane-parallel plate of constant index, air on both sides.
    let mut sys = SequentialSystem::new();
    sys.push(SeqSurface::flat(10e-3), 5e-3, Medium::Const(1.5));
    sys.push(SeqSurface::flat(10e-3), 0.0, Medium::Vacuum);
    assert_all_pass(&sys);
}
