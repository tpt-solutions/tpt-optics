//! Double Gauss lens checked against a published prescription.
//!
//! Source: the COMSOL "Double Gauss Lens" model (Ray Optics, model 1095781),
//! Table 1, which reproduces the f/1.7, 100.2 mm lens of Lautenbacher and
//! Brendel (Agfa, U.S. Patent 2,784,643), with glasses from the Schott and
//! OHARA catalogues. Blank radius cells are flat surfaces. Cemented pairs have
//! no air gap between them. The table gives the focal length and the back focal
//! distance (61.418 mm), which are checked here.
//!
//! Glass constants come from the datasheets: Schott N-SF5 (2007), and OHARA
//! S-LAM3 and S-BAH11 (OHARA's datasheets, decrypted for this test).

extern crate alloc;

use tpt_optics_core::dispersion::Sellmeier;
use tpt_optics_core::index::IndexAt;
use tpt_optics_core::material::{Material, MaterialDatabase};
use tpt_optics_core::units::REFERENCE_TEMPERATURE;
use tpt_optics_core::{Kelvin, Wavelength};

use tpt_optics_geometric::presets;
use tpt_optics_geometric::sequential::TraceContext;

fn catalog() -> MaterialDatabase {
    let mut db = MaterialDatabase::presets();
    let glasses = [
        // Schott N-SF5 (datasheet 2007): B and C of the Sellmeier formula.
        (
            "n-sf5",
            [1.524_818_89, 0.187_085_527, 1.427_290_15],
            [0.011_254_756, 0.058_899_539_2, 129.141_675],
        ),
        // OHARA S-LAM3: A_i and B_i (µm²) of the Sellmeier formula.
        (
            "s-lam3",
            [1.642_587_13, 0.239_634_610, 1.224_830_26],
            [0.008_682_460_20, 0.035_122_624_2, 116.604_369],
        ),
        // OHARA S-BAH11.
        (
            "s-bah11",
            [1.571_388_60, 0.147_869_313, 1.280_928_46],
            [0.009_108_079_36, 0.040_240_168_4, 130.399_367],
        ),
    ];
    for (name, a, b) in glasses {
        db.insert(Material::new(
            name,
            alloc::boxed::Box::new(Sellmeier::new(a, b)),
            0.0,
            REFERENCE_TEMPERATURE,
        ));
    }
    db
}

fn ctx() -> TraceContext {
    TraceContext {
        wavelength: Wavelength::from_nm(587.6),
        temperature: Kelvin::from_celsius(20.0),
    }
}

#[test]
fn glass_indices_match_the_datasheets() {
    let db = catalog();
    // Datasheets: N-SF5 nd 1.67271 vd 32.25; S-LAM3 nd 1.71700 vd 47.92;
    // S-BAH11 nd 1.66672 vd 48.32.
    for (name, nd, vd) in [
        ("n-sf5", 1.672_71, 32.25),
        ("s-lam3", 1.717_00, 47.92),
        ("s-bah11", 1.666_72, 48.32),
    ] {
        let m = db.get(name).unwrap();
        let n = m
            .index_at(Wavelength::from_nm(587.56), Kelvin::from_celsius(20.0))
            .n;
        let nf = m
            .index_at(Wavelength::from_nm(486.13), Kelvin::from_celsius(20.0))
            .n;
        let nc = m
            .index_at(Wavelength::from_nm(656.27), Kelvin::from_celsius(20.0))
            .n;
        assert!((n - nd).abs() < 2e-5, "{name} nd {n}");
        assert!(((n - 1.0) / (nf - nc) - vd).abs() < 0.05, "{name} vd");
    }
}

#[test]
fn double_gauss_focal_length_and_back_focal_distance() {
    let db = catalog();
    let sys = presets::double_gauss_100mm(&db).expect("prescription builds");
    let p = sys.paraxial(ctx());
    let efl = p.effective_focal_length().expect("focus") * 1e3;
    let bfl = p.back_focal_length().expect("focus") * 1e3;
    assert!(
        (efl - 100.2).abs() < 0.02,
        "EFL {efl} mm, published 100.2 mm"
    );
    assert!(
        (bfl - 61.418).abs() < 0.01,
        "BFL {bfl} mm, published 61.418 mm"
    );
}
