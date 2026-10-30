//! Glass map: rank crown–flint pairs for an achromatic doublet.
//!
//! For each crown and flint in the given lists, the search builds the
//! achromat with [`achromatic_doublet`], measures its on-axis RMS spot and
//! its axial colour, and sorts the pairs by spot size. The search is
//! exhaustive over the pairs given. It does not use a genetic algorithm.

use alloc::string::String;
use alloc::vec::Vec;

use tpt_optics_core::material::MaterialDatabase;
use tpt_optics_core::{Kelvin, Wavelength};
use tpt_optics_geometric::presets::{achromatic_doublet, AchromatSpec};
use tpt_optics_geometric::sequential::{rms_spot_radius, SequentialSystem, TraceContext};

use crate::optimize::{genetic_minimize, ParamBound};

/// One crown–flint pair and how well its achromat performs.
#[derive(Clone, Debug, PartialEq)]
pub struct GlassPair {
    /// Crown glass name.
    pub crown: String,
    /// Flint glass name.
    pub flint: String,
    /// On-axis RMS spot radius (m) over a pupil fan.
    pub rms_spot: f64,
    /// Axial colour: difference of the F and C focal lengths (m).
    pub axial_colour: f64,
}

/// Parameters shared by every pair in the map.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct MapSpec {
    /// Target focal length (m).
    pub focal_length: f64,
    /// Clear semi-diameter (m).
    pub semi_diameter: f64,
    /// Crown centre thickness (m).
    pub crown_thickness: f64,
    /// Flint centre thickness (m).
    pub flint_thickness: f64,
    /// Back-surface to image distance (m).
    pub image_distance: f64,
}

fn ctx(nm: f64) -> TraceContext {
    TraceContext {
        wavelength: Wavelength::from_nm(nm),
        temperature: Kelvin::from_celsius(20.0),
    }
}

/// Rank every crown–flint pair by its on-axis spot size, best first.
///
/// Pairs that cannot be built (an unknown glass, or a design with no focal
/// length) are left out. The list is sorted by `rms_spot`, and ties keep the
/// input order.
#[must_use]
pub fn rank_achromats(
    catalog: &MaterialDatabase,
    spec: MapSpec,
    crowns: &[&str],
    flints: &[&str],
) -> Vec<GlassPair> {
    let heights: Vec<f64> = (0..=10)
        .map(|i| (i as f64 - 5.0) / 5.0 * spec.semi_diameter)
        .collect();
    let mut out = Vec::new();
    for &crown in crowns {
        for &flint in flints {
            let spec = AchromatSpec {
                focal_length: spec.focal_length,
                crown,
                flint,
                semi_diameter: spec.semi_diameter,
                crown_thickness: spec.crown_thickness,
                flint_thickness: spec.flint_thickness,
                image_distance: spec.image_distance,
            };
            let Ok(sys) = achromatic_doublet(catalog, &spec) else {
                continue;
            };
            let Some(colour) = axial_colour(&sys) else {
                continue;
            };
            out.push(GlassPair {
                crown: crown.into(),
                flint: flint.into(),
                rms_spot: rms_spot_radius(&sys, ctx(587.6), 0.0, &heights),
                axial_colour: colour,
            });
        }
    }
    out.sort_by(|a, b| a.rms_spot.total_cmp(&b.rms_spot));
    out
}

/// One achromat found by [`genetic_achromat`].
#[derive(Clone, Debug, PartialEq)]
pub struct GeneticAchromat {
    /// Crown glass name.
    pub crown: String,
    /// Flint glass name.
    pub flint: String,
    /// Crown centre thickness (m).
    pub crown_thickness: f64,
    /// Flint centre thickness (m).
    pub flint_thickness: f64,
    /// On-axis RMS spot radius (m).
    pub rms_spot: f64,
    /// Axial colour (m).
    pub axial_colour: f64,
    /// Objective value at the optimum (lower is better).
    pub cost: f64,
}

/// Genetic search over crown and flint choices and element thicknesses.
///
/// The genome is `[crown index, flint index, crown thickness (mm), flint
/// thickness (mm)]`. For each genome the achromat is built with
/// [`achromatic_doublet`], and the cost combines the focal-length error, the
/// axial colour, and the RMS spot, all relative to the focal length. The
/// search is seeded and so reproducible. It is a global search over the
/// discrete glass choice and the thicknesses together.
///
/// Returns `None` when the glass lists are empty or no genome builds.
#[must_use]
pub fn genetic_achromat(
    catalog: &MaterialDatabase,
    spec: MapSpec,
    crowns: &[&str],
    flints: &[&str],
    population: usize,
    generations: usize,
    seed: u64,
) -> Option<GeneticAchromat> {
    if crowns.is_empty() || flints.is_empty() {
        return None;
    }
    let heights: Vec<f64> = (0..=10)
        .map(|i| (i as f64 - 5.0) / 5.0 * spec.semi_diameter)
        .collect();
    let decode = |g: &[f64]| {
        let pick = |v: f64, n: usize| ((v.max(0.0) as usize).min(n - 1)).min(n - 1);
        (
            crowns[pick(g[0], crowns.len())],
            flints[pick(g[1], flints.len())],
            g[2] * 1e-3,
            g[3] * 1e-3,
        )
    };
    let cost = |g: &[f64]| -> f64 {
        let (crown, flint, tc, tf) = decode(g);
        let aspec = AchromatSpec {
            focal_length: spec.focal_length,
            crown,
            flint,
            semi_diameter: spec.semi_diameter,
            crown_thickness: tc,
            flint_thickness: tf,
            image_distance: spec.image_distance,
        };
        let Ok(sys) = achromatic_doublet(catalog, &aspec) else {
            return 1e6;
        };
        let Some(efl) = sys.paraxial(ctx(587.6)).effective_focal_length() else {
            return 1e6;
        };
        let Some(colour) = axial_colour(&sys) else {
            return 1e6;
        };
        let spot = rms_spot_radius(&sys, ctx(587.6), 0.0, &heights);
        if !spot.is_finite() {
            return 1e6;
        }
        let f = spec.focal_length;
        let e = (efl - f) / f;
        let c = colour / f;
        let s = spot / f;
        e * e + c * c + s * s
    };
    let bounds = [
        ParamBound {
            lo: 0.0,
            hi: crowns.len() as f64 - 1e-9,
        },
        ParamBound {
            lo: 0.0,
            hi: flints.len() as f64 - 1e-9,
        },
        ParamBound { lo: 1.0, hi: 8.0 },
        ParamBound { lo: 1.0, hi: 6.0 },
    ];
    let result = genetic_minimize(&cost, &bounds, population, generations, seed);
    let (crown, flint, tc, tf) = decode(&result.x);
    let aspec = AchromatSpec {
        focal_length: spec.focal_length,
        crown,
        flint,
        semi_diameter: spec.semi_diameter,
        crown_thickness: tc,
        flint_thickness: tf,
        image_distance: spec.image_distance,
    };
    let sys = achromatic_doublet(catalog, &aspec).ok()?;
    Some(GeneticAchromat {
        crown: crown.into(),
        flint: flint.into(),
        crown_thickness: tc,
        flint_thickness: tf,
        rms_spot: rms_spot_radius(&sys, ctx(587.6), 0.0, &heights),
        axial_colour: axial_colour(&sys)?,
        cost: result.f,
    })
}

/// F-line minus C-line focal length of a system (m).
fn axial_colour(sys: &SequentialSystem) -> Option<f64> {
    let f = sys.paraxial(ctx(486.1)).effective_focal_length()?;
    let c = sys.paraxial(ctx(656.3)).effective_focal_length()?;
    Some(f - c)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_optics_core::dispersion::Sellmeier;
    use tpt_optics_core::material::Material;
    use tpt_optics_core::units::REFERENCE_TEMPERATURE;

    fn catalog() -> MaterialDatabase {
        let mut db = MaterialDatabase::presets();
        db.insert(Material::new(
            "F2",
            alloc::boxed::Box::new(Sellmeier::new(
                [1.345_333_59, 0.209_073_176, 0.937_357_162],
                [0.009_977_438_71, 0.047_045_076_7, 111.886_764],
            )),
            6.0e-6,
            REFERENCE_TEMPERATURE,
        ));
        db.insert(Material::new(
            "N-SF11",
            alloc::boxed::Box::new(Sellmeier::new(
                [1.737_596_95, 0.313_747_346, 1.898_781_01],
                [0.013_188_707, 0.062_306_814_2, 155.236_29],
            )),
            7.0e-6,
            REFERENCE_TEMPERATURE,
        ));
        db
    }

    fn spec() -> MapSpec {
        MapSpec {
            focal_length: 100e-3,
            semi_diameter: 6e-3,
            crown_thickness: 4e-3,
            flint_thickness: 2e-3,
            image_distance: 95e-3,
        }
    }

    #[test]
    fn every_pair_is_ranked_and_sorted_by_spot() {
        let cat = catalog();
        let pairs = rank_achromats(&cat, spec(), &["bk7"], &["f2", "n-sf11"]);
        assert_eq!(pairs.len(), 2);
        assert!(pairs[0].rms_spot <= pairs[1].rms_spot);
        assert_eq!(pairs[0].crown, "bk7");
    }

    #[test]
    fn unknown_glasses_are_skipped() {
        let cat = catalog();
        let pairs = rank_achromats(&cat, spec(), &["bk7", "nope"], &["f2"]);
        assert_eq!(pairs.len(), 1);
        assert_eq!(pairs[0].crown, "bk7");
    }

    #[test]
    fn the_map_reports_small_axial_colour_for_every_pair() {
        let cat = catalog();
        for p in rank_achromats(&cat, spec(), &["bk7"], &["f2", "n-sf11"]) {
            assert!(p.axial_colour.abs() < 0.02e-3, "{p:?}");
        }
    }
}

#[cfg(test)]
mod genetic_tests {
    use super::*;
    use tpt_optics_core::dispersion::Sellmeier;
    use tpt_optics_core::material::Material;
    use tpt_optics_core::units::REFERENCE_TEMPERATURE;

    fn catalog() -> MaterialDatabase {
        let mut db = MaterialDatabase::presets();
        db.insert(Material::new(
            "F2",
            alloc::boxed::Box::new(Sellmeier::new(
                [1.345_333_59, 0.209_073_176, 0.937_357_162],
                [0.009_977_438_71, 0.047_045_076_7, 111.886_764],
            )),
            6.0e-6,
            REFERENCE_TEMPERATURE,
        ));
        db
    }

    fn spec() -> MapSpec {
        MapSpec {
            focal_length: 100e-3,
            semi_diameter: 6e-3,
            crown_thickness: 4e-3,
            flint_thickness: 2e-3,
            image_distance: 95e-3,
        }
    }

    #[test]
    fn genetic_search_finds_an_achromat_near_the_target_focal_length() {
        let cat = catalog();
        let best = genetic_achromat(&cat, spec(), &["bk7"], &["f2"], 16, 12, 5).unwrap();
        assert_eq!(best.crown, "bk7");
        assert_eq!(best.flint, "f2");
        assert!(best.axial_colour.abs() < 0.5e-3, "{best:?}");
    }

    #[test]
    fn the_search_is_reproducible_for_a_seed() {
        let cat = catalog();
        let a = genetic_achromat(&cat, spec(), &["bk7", "f2"], &["f2"], 12, 8, 9).unwrap();
        let b = genetic_achromat(&cat, spec(), &["bk7", "f2"], &["f2"], 12, 8, 9).unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn empty_glass_lists_give_no_result() {
        assert!(genetic_achromat(&catalog(), spec(), &[], &["f2"], 8, 4, 1).is_none());
    }
}
