//! Materials and the material database.

use alloc::boxed::Box;
use alloc::collections::BTreeMap;
use alloc::string::String;
use alloc::sync::Arc;
use alloc::vec::Vec;

use crate::dispersion::{Cauchy, Dispersion, DrudeLorentz, Sellmeier, TablePoint, Tabulated};
use crate::index::{ComplexRefractiveIndex, IndexAt};
use crate::units::{Kelvin, Wavelength};

/// A named optical material: a dispersion model plus thermo-optic behaviour.
///
/// The temperature dependence is the linear thermo-optic correction
/// `n(T) = n(T₀) + (dn/dT)·(T − T₀)` applied to the real part, which is the
/// standard engineering model around the catalog reference temperature.
#[derive(Debug)]
pub struct Material {
    name: String,
    dispersion: Box<dyn Dispersion>,
    dn_dtemp: f64,
    reference_temp: Kelvin,
    validity_nm: Option<(f64, f64)>,
    thermo: Option<SchottThermo>,
}

/// Schott's temperature model of the refractive index (the form in Schott's
/// TIE-19 technical note), with `ΔT = T − T₀` in kelvin and `λ` in µm:
///
/// `Δn_abs(λ, T) = (n² − 1)/(2n) · [D₀ΔT + D₁ΔT² + D₂ΔT³ + (E₀ΔT + E₁ΔT²)/(λ² − λ_TK²)]`
///
/// The coefficients are printed on each Schott datasheet. A material with this
/// model is evaluated at its own index `n`, so the change follows the glass's
/// dispersion. Without it, [`Material`] uses the linear `dn/dT` term.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct SchottThermo {
    /// `D₀` (1/K).
    pub d0: f64,
    /// `D₁` (1/K²).
    pub d1: f64,
    /// `D₂` (1/K³).
    pub d2: f64,
    /// `E₀` (1/K).
    pub e0: f64,
    /// `E₁` (1/K²).
    pub e1: f64,
    /// Wavelength constant `λ_TK` (µm).
    pub lambda_tk_um: f64,
}

impl SchottThermo {
    /// Change of the absolute index at wavelength `wl_um` (µm) for a glass of
    /// index `n` at that wavelength, when the temperature differs by `dt` (K).
    #[must_use]
    pub fn delta_n_abs(&self, n: f64, wl_um: f64, dt: f64) -> f64 {
        let factor = (n * n - 1.0) / (2.0 * n);
        let polynomial = self.d0 * dt + self.d1 * dt * dt + self.d2 * dt * dt * dt;
        let resonant = (self.e0 * dt + self.e1 * dt * dt)
            / (wl_um * wl_um - self.lambda_tk_um * self.lambda_tk_um);
        factor * (polynomial + resonant)
    }
}

impl Material {
    /// Construct.
    pub fn new(
        name: impl Into<String>,
        dispersion: Box<dyn Dispersion>,
        dn_dtemp: f64,
        reference_temp: Kelvin,
    ) -> Self {
        Material {
            name: name.into(),
            dispersion,
            dn_dtemp,
            reference_temp,
            validity_nm: None,
            thermo: None,
        }
    }

    /// Use Schott's temperature model instead of the linear `dn/dT` term.
    #[must_use]
    pub fn with_schott_thermo(mut self, thermo: SchottThermo) -> Self {
        self.thermo = Some(thermo);
        self
    }

    /// Declare the wavelength range (nm) over which the dispersion model was
    /// fitted. Wavelengths outside it are flagged by
    /// [`Material::outside_validity`].
    #[must_use]
    pub fn with_validity(mut self, lo_nm: f64, hi_nm: f64) -> Self {
        self.validity_nm = Some((lo_nm, hi_nm));
        self
    }

    /// Declared validity range (nm), if any.
    pub fn validity_nm(&self) -> Option<(f64, f64)> {
        self.validity_nm
    }

    /// True when `wl` lies outside the declared validity range. A material
    /// without a declared range never reports outside.
    pub fn outside_validity(&self, wl: Wavelength) -> bool {
        match self.validity_nm {
            Some((lo, hi)) => {
                let nm = wl.nm();
                nm < lo || nm > hi
            }
            None => false,
        }
    }

    /// Abbe number `V_d = (n_d − 1)/(n_F − n_C)` from the dispersion model.
    /// It is infinite for a material with no dispersion (`n_F = n_C`).
    pub fn abbe_number(&self) -> f64 {
        let t = self.reference_temp;
        let n = |nm: f64| self.index_at(Wavelength::from_nm(nm), t).n;
        (n(587.6) - 1.0) / (n(486.1) - n(656.3))
    }

    /// Six-digit glass code: `nd` as three digits after the point (`n_d =
    /// 1.517` gives `517`), then `V_d` times ten (`V_d = 64.2` gives `642`).
    /// The value comes from the dispersion model, so the last digit can differ
    /// from a vendor code. Each part is clamped to three digits: an infinite
    /// `V_d` (no dispersion) reads as `999`, and a `V_d` that is not a number
    /// reads as `000`.
    pub fn glass_code(&self) -> String {
        let t = self.reference_temp;
        let nd = self.index_at(Wavelength::from_nm(587.6), t).n;
        let nd_digits = libm::round((nd - 1.0) * 1000.0).clamp(0.0, 999.0) as i64;
        let vd = self.abbe_number();
        let vd_digits = if vd.is_finite() {
            libm::round(vd * 10.0).clamp(0.0, 999.0) as i64
        } else if vd > 0.0 {
            999
        } else {
            0
        };
        alloc::format!("{nd_digits:03}{vd_digits:03}")
    }

    /// Material name.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Thermo-optic coefficient dn/dT in 1/K.
    pub fn dn_dtemp(&self) -> f64 {
        self.dn_dtemp
    }

    /// Catalog reference temperature.
    pub fn reference_temp(&self) -> Kelvin {
        self.reference_temp
    }
}

impl IndexAt for Material {
    fn index_at(&self, wavelength: Wavelength, temperature: Kelvin) -> ComplexRefractiveIndex {
        let mut idx = self.dispersion.index(wavelength);
        let dt = temperature.k() - self.reference_temp.k();
        match self.thermo {
            Some(thermo) => idx.n += thermo.delta_n_abs(idx.n, wavelength.um(), dt),
            None => idx.n += self.dn_dtemp * dt,
        }
        idx
    }
}

/// Vacuum: n = 1 at every wavelength and temperature.
#[derive(Copy, Clone, Debug, Default)]
pub struct Vacuum;

impl IndexAt for Vacuum {
    fn index_at(&self, _wavelength: Wavelength, _temperature: Kelvin) -> ComplexRefractiveIndex {
        ComplexRefractiveIndex::real(1.0)
    }
}

/// A wavelength-independent real index (homogeneous approximate medium).
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct ConstIndex(pub f64);

impl IndexAt for ConstIndex {
    fn index_at(&self, _wavelength: Wavelength, _temperature: Kelvin) -> ComplexRefractiveIndex {
        ComplexRefractiveIndex::real(self.0)
    }
}

impl Dispersion for ConstIndex {
    fn index(&self, _wavelength: Wavelength) -> ComplexRefractiveIndex {
        ComplexRefractiveIndex::real(self.0)
    }
}

/// Dry air at 15 °C, 101.325 kPa via the Edlén equation:
///
/// `(n − 1)·10⁸ = 6432.8 + 2949810/(146 − σ²) + 25540/(41 − σ²)`,
/// with σ = 1/λ in µm⁻¹. Valid for 0.2–2 µm. Outside that range the index is
/// NaN.
#[derive(Copy, Clone, Debug, Default)]
pub struct EdlenAir;

impl Dispersion for EdlenAir {
    fn index(&self, wavelength: Wavelength) -> ComplexRefractiveIndex {
        let um = wavelength.um();
        // Outside 0.2–2 µm the formula is not valid (it has poles near
        // 0.156 µm and 0.083 µm), so return NaN rather than a plausible number.
        if !(0.2..=2.0).contains(&um) {
            return ComplexRefractiveIndex::real(f64::NAN);
        }
        let sigma2 = 1.0 / (um * um);
        let n_minus_1 =
            1e-8 * (6432.8 + 2_949_810.0 / (146.0 - sigma2) + 25_540.0 / (41.0 - sigma2));
        ComplexRefractiveIndex::real(1.0 + n_minus_1)
    }
}

impl IndexAt for EdlenAir {
    fn index_at(&self, wavelength: Wavelength, _temperature: Kelvin) -> ComplexRefractiveIndex {
        Dispersion::index(self, wavelength)
    }
}

/// Refractive index of air at `wl` under the Edlén reference condition
/// (15 °C, 101 325 Pa). Catalogue indices are quoted for air, so this is the
/// factor between them and vacuum values.
#[must_use]
pub fn air_index(wl: Wavelength) -> f64 {
    Dispersion::index(&EdlenAir, wl).n
}

/// Convert a catalogue index quoted for air to vacuum: `n_vac = n_air · n_air(λ)`.
#[must_use]
pub fn vacuum_from_air_index(n_air: f64, wl: Wavelength) -> f64 {
    n_air * air_index(wl)
}

/// Convert a vacuum index to the value quoted for air: `n_air = n_vac / n_air(λ)`.
#[must_use]
pub fn air_from_vacuum_index(n_vac: f64, wl: Wavelength) -> f64 {
    n_vac / air_index(wl)
}

/// Table for a built-in material. The data is fixed in the source, so a
/// failure here is a bug in this crate and panics.
fn builtin_table(points: &[TablePoint]) -> Tabulated {
    Tabulated::from_points(points).expect("built-in dispersion table is valid")
}

/// A registry of named materials.
#[derive(Debug, Default)]
pub struct MaterialDatabase {
    /// Materials keyed by lower-cased name: lookup is logarithmic and ignores case.
    entries: BTreeMap<String, Arc<Material>>,
}

impl MaterialDatabase {
    /// Empty database.
    pub fn new() -> Self {
        MaterialDatabase::default()
    }

    /// Insert (replaces an existing entry with the same name).
    pub fn insert(&mut self, material: Material) {
        let key = material.name().to_ascii_lowercase();
        self.entries.insert(key, Arc::new(material));
    }

    /// Look up by name. Case is ignored, so `"n-bk7"` finds `"N-BK7"`.
    pub fn get(&self, name: &str) -> Option<Arc<Material>> {
        self.entries.get(&name.to_ascii_lowercase()).map(Arc::clone)
    }

    /// Number of materials.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// True when empty.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Registered name closest to `query`, ignoring case, if it is within
    /// `max_distance` edits. Use it for "did you mean" messages.
    pub fn closest_name(&self, query: &str, max_distance: usize) -> Option<&str> {
        let query = query.to_ascii_lowercase();
        self.entries
            .values()
            .map(|m| {
                (
                    edit_distance(&query, &m.name().to_ascii_lowercase()),
                    m.name(),
                )
            })
            .filter(|(d, _)| *d <= max_distance)
            .min_by_key(|(d, _)| *d)
            .map(|(_, name)| name)
    }

    /// Names of all materials, sorted.
    pub fn names(&self) -> Vec<&str> {
        let mut names: Vec<&str> = self.entries.values().map(|m| m.name()).collect();
        names.sort_unstable();
        names
    }

    /// Database with the built-in reference materials.
    ///
    /// The embedded catalog data is an *illustrative reference subset*
    /// sufficient for testing and examples; production work should import
    /// vendor catalogs via `tpt-optics-io`.
    pub fn presets() -> Self {
        let mut db = MaterialDatabase::new();
        let mut add = |name: &str, disp: Box<dyn Dispersion>, dn: f64| {
            let material = Material::new(name, disp, dn, crate::units::REFERENCE_TEMPERATURE);
            let material = match name {
                // Schott's fit for BK7 covers 0.3–2.5 µm.
                "bk7" | "n_bk7" => material.with_validity(300.0, 2500.0),
                // Malitson's fit for fused silica covers 0.21–3.71 µm.
                "fused_silica" => material.with_validity(210.0, 3710.0),
                _ => material,
            };
            db.insert(material);
        };

        add("vacuum", ConstIndex(1.0).boxed(), 0.0);
        add("air", Box::new(EdlenAir), 0.0);
        add(
            "fused_silica",
            Sellmeier::fused_silica().boxed(),
            1.0e-5, // dn/dT ≈ +1.0e-5 /K near 589 nm
        );
        add("bk7", Sellmeier::bk7().boxed(), 3.0e-6);
        add("n_bk7", Sellmeier::bk7().boxed(), 3.0e-6);
        add(
            "water",
            builtin_table(&[
                TablePoint {
                    wavelength_m: 0.30e-6,
                    n: 1.349,
                    k: 0.0,
                },
                TablePoint {
                    wavelength_m: 0.40e-6,
                    n: 1.344,
                    k: 0.0,
                },
                TablePoint {
                    wavelength_m: 0.50e-6,
                    n: 1.336,
                    k: 0.0,
                },
                TablePoint {
                    wavelength_m: 0.589e-6,
                    n: 1.3330,
                    k: 0.0,
                },
                TablePoint {
                    wavelength_m: 0.70e-6,
                    n: 1.331,
                    k: 0.0,
                },
                TablePoint {
                    wavelength_m: 1.00e-6,
                    n: 1.328,
                    k: 0.0,
                },
                TablePoint {
                    wavelength_m: 1.50e-6,
                    n: 1.321,
                    k: 0.0,
                },
            ])
            .boxed(),
            -1.0e-4,
        );
        add("diamond", Sellmeier::diamond().boxed(), 1.0e-5);
        add(
            "silicon",
            builtin_table(&[
                TablePoint {
                    wavelength_m: 1.20e-6,
                    n: 3.528,
                    k: 1.0e-6,
                },
                TablePoint {
                    wavelength_m: 1.50e-6,
                    n: 3.481,
                    k: 0.0,
                },
                TablePoint {
                    wavelength_m: 2.00e-6,
                    n: 3.453,
                    k: 0.0,
                },
                TablePoint {
                    wavelength_m: 2.50e-6,
                    n: 3.439,
                    k: 0.0,
                },
                TablePoint {
                    wavelength_m: 3.00e-6,
                    n: 3.433,
                    k: 0.0,
                },
                TablePoint {
                    wavelength_m: 4.00e-6,
                    n: 3.425,
                    k: 0.0,
                },
            ])
            .boxed(),
            1.6e-4,
        );
        add(
            "gold",
            DrudeLorentz::new(1.0, Vec::new())
                .drude(1.37e16, 1.2e14)
                .boxed(),
            0.0,
        );
        add(
            "soda_lime",
            Cauchy::new(1.513, 0.00359, 0.0).boxed(),
            2.5e-6,
        );
        db
    }
}

/// Levenshtein distance between two strings, counted in characters.
fn edit_distance(a: &str, b: &str) -> usize {
    let b: Vec<char> = b.chars().collect();
    let mut prev: Vec<usize> = (0..=b.len()).collect();
    for (i, ca) in a.chars().enumerate() {
        let mut cur = Vec::with_capacity(b.len() + 1);
        cur.push(i + 1);
        for (j, cb) in b.iter().enumerate() {
            let subst = prev[j] + usize::from(ca != *cb);
            cur.push(subst.min(prev[j + 1] + 1).min(cur[j] + 1));
        }
        prev = cur;
    }
    prev[b.len()]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn edlen_air_index() {
        // n_air(589 nm) ≈ 1.000277 (15 °C, 1 atm)
        let n = EdlenAir.real_index(Wavelength::from_nm(589.0));
        assert!((n - 1.000277).abs() < 1e-6, "n = {n}");
        // Dispersion: blue bends more.
        assert!(EdlenAir.real_index(Wavelength::from_nm(450.0)) > n);
    }

    #[test]
    fn presets_lookup_and_thermal() {
        let db = MaterialDatabase::presets();
        assert!(db.len() >= 9);
        let bk7 = db.get("bk7").expect("bk7 preset");
        let n20 = bk7.index_at(Wavelength::from_nm(587.6), Kelvin::from_celsius(20.0));
        let n40 = bk7.index_at(Wavelength::from_nm(587.6), Kelvin::from_celsius(40.0));
        assert!((n20.n - 1.5168).abs() < 2e-4);
        // dn/dT = 3e-6/K over 20 K
        assert!((n40.n - n20.n - 6.0e-5).abs() < 1e-7);

        let water = db.get("water").unwrap();
        let nw = water.index_at(Wavelength::from_nm(589.0), Kelvin::from_celsius(20.0));
        assert!((nw.n - 1.333).abs() < 5e-4);
        assert!(db.get("unobtainium").is_none());
    }

    #[test]
    fn vacuum_and_const() {
        let wl = Wavelength::from_nm(500.0);
        assert_eq!(Vacuum.index_at(wl, Kelvin::from_k(300.0)).n, 1.0);
        assert_eq!(ConstIndex(1.77).index_at(wl, Kelvin::from_k(300.0)).n, 1.77);
    }

    #[test]
    fn insert_replaces() {
        let mut db = MaterialDatabase::new();
        db.insert(Material::new(
            "test",
            ConstIndex(1.5).boxed(),
            0.0,
            Kelvin::from_k(293.0),
        ));
        db.insert(Material::new(
            "test",
            ConstIndex(2.5).boxed(),
            0.0,
            Kelvin::from_k(293.0),
        ));
        assert_eq!(db.len(), 1);
        assert!(
            (db.get("test")
                .unwrap()
                .index_at(Wavelength::from_nm(500.0), Kelvin::from_k(293.0))
                .n
                - 2.5)
                .abs()
                < 1e-15
        );
    }
}

#[cfg(test)]
mod lookup_tests {
    use super::*;
    use crate::units::Wavelength;

    #[test]
    fn lookup_ignores_case() {
        // Zemax and catalog files do not agree on capitalisation.
        let mut db = MaterialDatabase::new();
        db.insert(Material::new(
            "N-SF11",
            Box::new(ConstIndex(1.78)),
            0.0,
            crate::units::REFERENCE_TEMPERATURE,
        ));
        assert!(db.get("n-sf11").is_some());
        assert!(db.get("N-SF11").is_some());
        assert!(db.get("N-SF12").is_none());
    }

    #[test]
    fn edlen_outside_its_range_is_nan() {
        // The formula has poles at about 0.156 µm and 0.083 µm, and it is
        // only valid for 0.2–2 µm. Outside that range it must not return a
        // number that looks valid.
        assert!(EdlenAir.real_index(Wavelength::from_um(0.1)).is_nan());
        assert!(EdlenAir.real_index(Wavelength::from_um(3.0)).is_nan());
        assert!(EdlenAir.real_index(Wavelength::from_um(0.5)).is_finite());
    }
}

#[cfg(test)]
mod suggestion_tests {
    use super::*;

    #[test]
    fn suggests_a_near_miss_name() {
        let db = MaterialDatabase::presets();
        assert_eq!(db.closest_name("n-bk7", 1), Some("n_bk7"));
        assert_eq!(db.closest_name("bk8", 1), Some("bk7"));
    }

    #[test]
    fn no_suggestion_for_an_unrelated_name() {
        let db = MaterialDatabase::presets();
        assert_eq!(db.closest_name("unobtainium", 2), None);
    }
}

#[cfg(test)]
mod metadata_tests {
    use super::*;

    #[test]
    fn bk7_abbe_number_and_glass_code_match_the_datasheet() {
        let db = MaterialDatabase::presets();
        let bk7 = db.get("bk7").unwrap();
        assert!(
            (bk7.abbe_number() - 64.17).abs() < 0.1,
            "{}",
            bk7.abbe_number()
        );
        // Schott code 517642; the model's V_d may land one unit off.
        let code = bk7.glass_code();
        assert!(code.starts_with("517"), "{code}");
        let vd: i64 = code[3..].parse().unwrap();
        assert!((vd - 642).abs() <= 1, "{code}");
    }

    #[test]
    fn validity_range_is_reported_outside_the_fit() {
        let db = MaterialDatabase::presets();
        let bk7 = db.get("bk7").unwrap();
        assert!(!bk7.outside_validity(Wavelength::from_nm(587.6)));
        assert!(bk7.outside_validity(Wavelength::from_nm(5000.0)));
        assert!(bk7.outside_validity(Wavelength::from_nm(200.0)));
    }

    #[test]
    fn material_without_a_declared_range_never_warns() {
        let db = MaterialDatabase::presets();
        let water = db.get("water").unwrap();
        assert_eq!(water.validity_nm(), None);
        assert!(!water.outside_validity(Wavelength::from_nm(50_000.0)));
    }
}

#[cfg(test)]
mod glass_code_tests {
    use super::*;
    #[test]
    fn glass_code_stays_six_digits_without_dispersion() {
        // n_F = n_C, so V_d is infinite. The code used to saturate and print
        // a 19-digit number.
        let flat = Material::new(
            "flat",
            alloc::boxed::Box::new(ConstIndex(1.5)),
            0.0,
            Kelvin::from_celsius(20.0),
        );
        assert!(flat.abbe_number().is_infinite());
        assert_eq!(flat.glass_code(), "500999");
    }
}

#[cfg(test)]
mod air_vacuum_tests {
    use super::*;

    #[test]
    fn glass_index_in_vacuum_is_about_4e4_above_its_index_in_air() {
        // Catalogue indices are for air. At 587.6 nm, n_air ≈ 1.000277, so a
        // glass with n = 1.5168 in air is about 4.2e-4 higher in vacuum.
        let wl = Wavelength::from_nm(587.6);
        let vac = vacuum_from_air_index(1.5168, wl);
        let diff = vac - 1.5168;
        assert!(diff > 4.0e-4 && diff < 4.4e-4, "difference {diff}");
        assert!((air_from_vacuum_index(vac, wl) - 1.5168).abs() < 1e-12);
    }
}

#[cfg(test)]
mod schott_n_sf5_tests {
    use super::*;
    use crate::dispersion::Sellmeier;
    use crate::units::REFERENCE_TEMPERATURE;

    /// N-SF5 from the Schott datasheet (2007): Sellmeier constants, indices
    /// at the standard lines, and the dn/dT coefficients. Expected values are
    /// the datasheet's, so the tests check this crate against that source.
    fn n_sf5() -> Material {
        Material::new(
            "N-SF5",
            alloc::boxed::Box::new(Sellmeier::new(
                [1.524_818_89, 0.187_085_527, 1.427_290_15],
                [0.011_254_756, 0.058_899_539_2, 129.141_675],
            )),
            0.0,
            REFERENCE_TEMPERATURE,
        )
        .with_schott_thermo(SchottThermo {
            d0: -2.51e-7,
            d1: 1.07e-8,
            d2: -2.40e-11,
            e0: 7.85e-7,
            e1: 1.15e-9,
            lambda_tk_um: 0.278,
        })
    }

    fn n_at(m: &Material, nm: f64, celsius: f64) -> f64 {
        m.index_at(Wavelength::from_nm(nm), Kelvin::from_celsius(celsius))
            .n
    }

    #[test]
    fn sellmeier_reproduces_the_datasheet_indices() {
        let m = n_sf5();
        // Datasheet: n_d 1.67271, n_F 1.68750, n_C 1.66664, n_e 1.67763.
        assert!((n_at(&m, 587.6, 20.0) - 1.672_71).abs() < 2e-5);
        assert!((n_at(&m, 486.1, 20.0) - 1.687_50).abs() < 2e-5);
        assert!((n_at(&m, 656.3, 20.0) - 1.666_64).abs() < 2e-5);
        assert!((n_at(&m, 546.1, 20.0) - 1.677_63).abs() < 2e-5);
    }

    #[test]
    fn thermo_optic_slope_matches_the_datasheet_table() {
        let m = n_sf5();
        // Table "Temperature Coefficients of Refractive Index", +20/+40 °C,
        // absolute dn/dT in units of 1e-6 /K: 0.4 at 1060 nm, 2.0 at the e-line.
        let slope = |nm: f64, lo: f64, hi: f64| (n_at(&m, nm, hi) - n_at(&m, nm, lo)) / (hi - lo);
        let s1060 = slope(1060.0, 20.0, 40.0) * 1e6;
        let se = slope(546.1, 20.0, 40.0) * 1e6;
        assert!((s1060 - 0.4).abs() < 0.2, "1060 nm slope {s1060} vs 0.4");
        assert!((se - 2.0).abs() < 0.3, "e-line slope {se} vs 2.0");
        // And −40/−20 °C at 1060 nm: −0.5e-6 /K.
        let cold = slope(1060.0, -40.0, -20.0) * 1e6;
        assert!(
            (cold + 0.5).abs() < 0.3,
            "cold 1060 nm slope {cold} vs -0.5"
        );
    }
}
