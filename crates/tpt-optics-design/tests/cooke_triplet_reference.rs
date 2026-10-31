//! Cooke triplet checked against a published optimised prescription.
//!
//! Source: the Optiland documentation, "Tutorial 3d: Optimization Case Study:
//! Cooke Triplet" (MIT-licensed project). Its final prescription, optimised to
//! a 50 mm focal length, is reproduced below in millimetres. The tutorial
//! reports a final focal length of 50.002 mm.
//!
//! Glasses: SK16 and F2 are registered at their published d-line indices
//! (Schott datasheets: N-SK16 `nd = 1.62041`, F2 `nd = 1.62004`). The check is
//! therefore the d-line focal length, which is what the chromatic-independent
//! part of the design fixes.

extern crate alloc;

use tpt_optics_core::material::{ConstIndex, Material, MaterialDatabase};
use tpt_optics_core::units::REFERENCE_TEMPERATURE;
use tpt_optics_core::{Kelvin, Wavelength};
use tpt_optics_geometric::presets;
use tpt_optics_geometric::sequential::TraceContext;

fn catalog() -> MaterialDatabase {
    let mut db = MaterialDatabase::presets();
    for (name, nd) in [("sk16", 1.62041), ("f2", 1.62004)] {
        db.insert(Material::new(
            name,
            alloc::boxed::Box::new(ConstIndex(nd)),
            0.0,
            REFERENCE_TEMPERATURE,
        ));
    }
    db
}

#[test]
fn published_cooke_triplet_focal_length_is_reproduced() {
    let db = catalog();
    let sys = presets::cooke_triplet_50mm(&db).expect("published prescription builds");
    let ctx = TraceContext {
        wavelength: Wavelength::from_nm(587.6),
        temperature: Kelvin::from_celsius(20.0),
    };
    let efl_mm = sys
        .paraxial(ctx)
        .effective_focal_length()
        .expect("has a focus")
        * 1e3;
    assert!(
        (efl_mm - 50.002).abs() < 0.01,
        "EFL = {efl_mm} mm, published 50.002 mm"
    );
}
