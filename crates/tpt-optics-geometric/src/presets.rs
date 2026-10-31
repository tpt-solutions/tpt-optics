//! Common lens prescriptions, built with [`SystemBuilder`].
//!
//! Every preset takes its parameters explicitly and none is optimised. Each
//! one ends with a flat image plane, so the system can be traced directly
//! (for example with `rms_spot_radius`).

use tpt_optics_core::material::MaterialDatabase;

use crate::builder::{BuildError, SystemBuilder};
use crate::sequential::SequentialSystem;

/// Plano-convex lens: a spherical front surface of `radius` (m) and a flat
/// back surface. `thickness` is the centre thickness, and `image_distance` is
/// the distance from the back surface to the image plane.
///
/// # Errors
/// Returns the [`BuildError`] of the underlying prescription, for example an
/// unknown `glass`.
pub fn plano_convex(
    catalog: &MaterialDatabase,
    radius: f64,
    semi_diameter: f64,
    thickness: f64,
    glass: &str,
    image_distance: f64,
) -> Result<SequentialSystem, BuildError> {
    SystemBuilder::new(catalog)
        .surface(radius, semi_diameter, thickness, glass)
        .surface(0.0, semi_diameter, image_distance, "")
        .surface(0.0, semi_diameter, 0.0, "")
        .build()
}

/// Biconvex or meniscus singlet with radii `r1` and `r2` (m), centre
/// thickness `thickness` (m), and image plane `image_distance` after the back
/// surface.
///
/// # Errors
/// Returns the [`BuildError`] of the underlying prescription, for example an
/// unknown `glass`.
pub fn singlet(
    catalog: &MaterialDatabase,
    r1: f64,
    r2: f64,
    semi_diameter: f64,
    thickness: f64,
    glass: &str,
    image_distance: f64,
) -> Result<SequentialSystem, BuildError> {
    SystemBuilder::new(catalog)
        .surface(r1, semi_diameter, thickness, glass)
        .surface(r2, semi_diameter, image_distance, "")
        .surface(0.0, semi_diameter, 0.0, "")
        .build()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequential::TraceContext;
    use tpt_optics_core::index::IndexAt;
    use tpt_optics_core::{Kelvin, Wavelength};

    /// Thick-lens focal length of a lens in air (m), from the lensmaker's
    /// formula with the centre thickness.
    fn thick_lens_efl(n: f64, r1: f64, r2: f64, t: f64) -> f64 {
        let power = (n - 1.0) * (1.0 / r1 - 1.0 / r2 + (n - 1.0) * t / (n * r1 * r2));
        1.0 / power
    }

    fn ctx() -> TraceContext {
        TraceContext {
            wavelength: Wavelength::from_nm(587.6),
            temperature: Kelvin::from_celsius(20.0),
        }
    }

    #[test]
    fn singlet_efl_matches_thick_lens_formula() {
        let catalog = MaterialDatabase::presets();
        let (r1, r2, t) = (100e-3, -100e-3, 5e-3);
        let sys = singlet(&catalog, r1, r2, 12.5e-3, t, "bk7", 95e-3).unwrap();
        let n = catalog
            .get("bk7")
            .unwrap()
            .index_at(ctx().wavelength, ctx().temperature)
            .n;
        let expected = thick_lens_efl(n, r1, r2, t);
        let traced = sys.paraxial(ctx()).effective_focal_length().unwrap();
        assert!(
            (traced - expected).abs() < 1e-9 * expected,
            "{traced} vs {expected}"
        );
    }

    #[test]
    fn plano_convex_efl_is_r_over_n_minus_one() {
        // With a flat back surface and a thin lens, f = R/(n − 1).
        let catalog = MaterialDatabase::presets();
        let r = 100e-3;
        let sys = plano_convex(&catalog, r, 12.5e-3, 1e-6, "bk7", 95e-3).unwrap();
        let n = catalog
            .get("bk7")
            .unwrap()
            .index_at(ctx().wavelength, ctx().temperature)
            .n;
        let traced = sys.paraxial(ctx()).effective_focal_length().unwrap();
        assert!((traced - r / (n - 1.0)).abs() < 1e-4 * r, "{traced}");
    }

    #[test]
    fn unknown_glass_propagates_from_presets() {
        let catalog = MaterialDatabase::presets();
        let err = plano_convex(&catalog, 100e-3, 12.5e-3, 5e-3, "bk8", 95e-3).unwrap_err();
        assert!(matches!(err, BuildError::UnknownGlass { .. }), "{err:?}");
    }
}

/// Index and Abbe number of a glass at the d, F and C lines.
fn glass_constants(catalog: &MaterialDatabase, glass: &str) -> Option<(f64, f64)> {
    use tpt_optics_core::index::IndexAt;
    use tpt_optics_core::{Kelvin, Wavelength};
    let m = catalog.get(glass)?;
    let t = Kelvin::from_celsius(20.0);
    let n = |nm: f64| m.index_at(Wavelength::from_nm(nm), t).n;
    let nd = n(587.6);
    let abbe = (nd - 1.0) / (n(486.1) - n(656.3));
    Some((nd, abbe))
}

/// Geometry and glasses of a cemented achromat.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct AchromatSpec<'a> {
    /// Target focal length (m).
    pub focal_length: f64,
    /// Crown glass name.
    pub crown: &'a str,
    /// Flint glass name.
    pub flint: &'a str,
    /// Clear semi-diameter of every surface (m).
    pub semi_diameter: f64,
    /// Crown centre thickness (m).
    pub crown_thickness: f64,
    /// Flint centre thickness (m).
    pub flint_thickness: f64,
    /// Distance from the last surface to the image plane (m).
    pub image_distance: f64,
}

/// Cemented achromatic doublet described by an [`AchromatSpec`].
///
/// The thin-lens achromat condition sets the element powers in proportion to
/// the Abbe numbers, `φ_crown = φ·V_c/(V_c − V_f)`, and the cemented surface
/// is flat. Element thicknesses shift those powers, so a 2×2 Newton
/// iteration on the traced d-line focal length and the F–C focus difference
/// corrects them. The result has the requested focal length and equal F and
/// C foci. It does not correct spherical aberration, so optimise the radii
/// with the operands module if the residual aberrations matter.
///
///
/// # Errors
/// [`BuildError::UnknownGlass`] for a glass that is not in the catalog, and
/// [`BuildError::InvalidRow`] (via the builder) when the power is zero.
pub fn achromatic_doublet(
    catalog: &MaterialDatabase,
    spec: &AchromatSpec<'_>,
) -> Result<SequentialSystem, BuildError> {
    let AchromatSpec {
        focal_length,
        crown,
        flint,
        semi_diameter,
        crown_thickness,
        flint_thickness,
        image_distance,
    } = *spec;
    let (nc, vc) = glass_constants(catalog, crown).ok_or_else(|| BuildError::UnknownGlass {
        name: crown.into(),
        suggestion: catalog.closest_name(crown, 2).map(Into::into),
    })?;
    let (nf, vf) = glass_constants(catalog, flint).ok_or_else(|| BuildError::UnknownGlass {
        name: flint.into(),
        suggestion: catalog.closest_name(flint, 2).map(Into::into),
    })?;
    // Unknowns: `x[0]` scales both element powers, and `x[1]` scales the
    // flint power relative to the crown. The thin-lens powers ignore the
    // element thicknesses, so Newton's method corrects them against the traced
    // focal lengths: the d-line focal length must equal the request, and the
    // F and C foci must coincide.
    let build = |x: [f64; 2]| {
        let power = x[0] / focal_length;
        let phi_crown = power * vc / (vc - vf);
        let phi_flint = -power * vf / (vc - vf) * x[1];
        // Flat cemented surface: φ_crown = (n_c − 1)/R₁ and φ_flint = −(n_f − 1)/R₃.
        let r1 = (nc - 1.0) / phi_crown;
        let r3 = -(nf - 1.0) / phi_flint;
        SystemBuilder::new(catalog)
            .surface(r1, semi_diameter, crown_thickness, crown)
            .surface(0.0, semi_diameter, flint_thickness, flint)
            .surface(r3, semi_diameter, image_distance, "")
            .surface(0.0, semi_diameter, 0.0, "")
            .build()
    };
    let efl_at = |sys: &SequentialSystem, nm: f64| {
        let ctx = crate::sequential::TraceContext {
            wavelength: tpt_optics_core::Wavelength::from_nm(nm),
            temperature: tpt_optics_core::Kelvin::from_celsius(20.0),
        };
        sys.paraxial(ctx).effective_focal_length()
    };
    let residual = |x: [f64; 2]| -> Result<[f64; 2], BuildError> {
        let sys = build(x)?;
        let (Some(d), Some(f), Some(c)) = (
            efl_at(&sys, 587.6),
            efl_at(&sys, 486.1),
            efl_at(&sys, 656.3),
        ) else {
            return Ok([f64::NAN, f64::NAN]);
        };
        Ok([d - focal_length, f - c])
    };
    let mut x = [1.0, 1.0];
    for _ in 0..30 {
        let r = residual(x)?;
        if !(r[0].is_finite() && r[1].is_finite()) {
            break;
        }
        if libm::fabs(r[0]) < 1e-9 * focal_length && libm::fabs(r[1]) < 1e-9 * focal_length {
            break;
        }
        // Finite-difference Jacobian.
        let h = 1e-6;
        let rx = residual([x[0] + h, x[1]])?;
        let ry = residual([x[0], x[1] + h])?;
        let j = [
            [(rx[0] - r[0]) / h, (ry[0] - r[0]) / h],
            [(rx[1] - r[1]) / h, (ry[1] - r[1]) / h],
        ];
        let det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
        if det == 0.0 || !det.is_finite() {
            break;
        }
        let dx = [
            (j[1][1] * r[0] - j[0][1] * r[1]) / det,
            (-j[1][0] * r[0] + j[0][0] * r[1]) / det,
        ];
        x = [x[0] - dx[0], x[1] - dx[1]];
    }
    build(x)
}

/// Three-element crown–flint–crown lens, 100 mm focal length, about f/5.
///
/// This is an optimised design, not the classical Cooke prescription. Its radii
/// and gaps were found with the operand optimiser (EFL 100 mm, focus on the
/// image plane, on-axis and 5° spots), with the gaps bounded to physical
/// values. It needs the glasses `bk7` (crown) and `n-sf11` (flint) in the
/// catalogue. Traced with the built-in values, it has a focal length of about
/// 100.0 mm and an on-axis RMS spot of about 10 µm at a 10 mm semi-diameter.
/// Check the performance with your own catalogue before relying on it.
///
/// # Errors
/// Returns the [`BuildError`] of the underlying prescription, for example an
/// unknown glass.
pub fn optimised_triplet_100mm(catalog: &MaterialDatabase) -> Result<SequentialSystem, BuildError> {
    // Radii and gaps in millimetres, from the optimiser.
    let sd = 10.0e-3;
    let r = [42.22, -260.0, -55.85, 153.51, 43.12, -290.11];
    let t = [6.13, 4.47, 1.0, 19.05, 8.8];
    SystemBuilder::new(catalog)
        .surface(r[0] * 1e-3, sd, t[0] * 1e-3, "bk7")
        .surface(r[1] * 1e-3, sd, t[1] * 1e-3, "")
        .surface(r[2] * 1e-3, sd, t[2] * 1e-3, "n-sf11")
        .surface(r[3] * 1e-3, sd, t[3] * 1e-3, "")
        .surface(r[4] * 1e-3, sd, t[4] * 1e-3, "bk7")
        .surface(r[5] * 1e-3, sd, 89.22e-3, "")
        .surface(0.0, sd, 0.0, "")
        .build()
}

#[cfg(test)]
mod achromat_tests {
    use super::*;
    use crate::sequential::TraceContext;
    use tpt_optics_core::dispersion::Sellmeier;
    use tpt_optics_core::material::Material;
    use tpt_optics_core::units::REFERENCE_TEMPERATURE;
    use tpt_optics_core::{Kelvin, Wavelength};

    /// Catalogue with a crown (N-BK7) and a flint (F2, from its Sellmeier
    /// coefficients), so the test does not depend on the io crate.
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

    fn efl_at(sys: &SequentialSystem, nm: f64) -> f64 {
        let ctx = TraceContext {
            wavelength: Wavelength::from_nm(nm),
            temperature: Kelvin::from_celsius(20.0),
        };
        sys.paraxial(ctx).effective_focal_length().unwrap()
    }

    #[test]
    fn achromat_has_the_requested_focal_length_and_equal_f_and_c_foci() {
        let cat = catalog();
        let f = 100e-3;
        let spec = AchromatSpec {
            focal_length: f,
            crown: "bk7",
            flint: "f2",
            semi_diameter: 12.5e-3,
            crown_thickness: 4e-3,
            flint_thickness: 2e-3,
            image_distance: 95e-3,
        };
        let sys = achromatic_doublet(&cat, &spec).unwrap();
        let efl_d = efl_at(&sys, 587.6);
        assert!((efl_d - f).abs() < 0.03 * f, "d-line EFL {efl_d}");
        let efl_f = efl_at(&sys, 486.1);
        let efl_c = efl_at(&sys, 656.3);
        // A single crown would differ by about f/V ≈ 1.6 mm between F and C.
        assert!((efl_f - efl_c).abs() < 0.02e-3, "F {efl_f} vs C {efl_c}");
    }

    #[test]
    fn unknown_flint_is_reported() {
        let cat = catalog();
        let spec = AchromatSpec {
            focal_length: 100e-3,
            crown: "bk7",
            flint: "f9",
            semi_diameter: 12.5e-3,
            crown_thickness: 4e-3,
            flint_thickness: 2e-3,
            image_distance: 95e-3,
        };
        let err = achromatic_doublet(&cat, &spec).unwrap_err();
        assert!(
            matches!(err, BuildError::UnknownGlass { ref name, .. } if name == "f9"),
            "{err:?}"
        );
    }
}

#[cfg(test)]
mod triplet_tests {
    use super::*;
    use crate::sequential::rms_spot_radius;
    use crate::sequential::TraceContext;
    use tpt_optics_core::dispersion::Sellmeier;
    use tpt_optics_core::material::Material;
    use tpt_optics_core::units::REFERENCE_TEMPERATURE;
    use tpt_optics_core::{Kelvin, Wavelength};

    fn catalog() -> MaterialDatabase {
        let mut db = MaterialDatabase::presets();
        db.insert(Material::new(
            "n-sf11",
            alloc::boxed::Box::new(Sellmeier::new(
                [1.737_596_95, 0.313_747_346, 1.898_781_01],
                [0.013_188_707, 0.062_306_814_2, 155.236_29],
            )),
            7.0e-6,
            REFERENCE_TEMPERATURE,
        ));
        db
    }

    #[test]
    fn optimised_triplet_has_its_focal_length_and_a_small_axial_spot() {
        let cat = catalog();
        let sys = optimised_triplet_100mm(&cat).unwrap();
        let ctx = TraceContext {
            wavelength: Wavelength::from_nm(587.6),
            temperature: Kelvin::from_celsius(20.0),
        };
        let efl = sys.paraxial(ctx).effective_focal_length().unwrap();
        assert!((efl - 0.100).abs() < 0.5e-3, "EFL {efl}");
        let heights: alloc::vec::Vec<f64> = (0..=8).map(|i| (i as f64 - 4.0) * 2e-3).collect();
        let spot = rms_spot_radius(&sys, ctx, 0.0, &heights);
        assert!(spot < 20e-6, "on-axis spot {spot}");
    }
}
