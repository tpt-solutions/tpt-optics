//! Materials and the material database.

use alloc::boxed::Box;
use alloc::string::{String, ToString};
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
        }
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
        idx.n += self.dn_dtemp * (temperature.k() - self.reference_temp.k());
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
/// with σ = 1/λ in µm⁻¹. Valid ≈ 0.2–2 µm.
#[derive(Copy, Clone, Debug, Default)]
pub struct EdlenAir;

impl Dispersion for EdlenAir {
    fn index(&self, wavelength: Wavelength) -> ComplexRefractiveIndex {
        let sigma2 = 1.0 / (wavelength.um() * wavelength.um());
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

/// A registry of named materials.
#[derive(Debug, Default)]
pub struct MaterialDatabase {
    entries: Vec<(String, Arc<Material>)>,
}

impl MaterialDatabase {
    /// Empty database.
    pub fn new() -> Self {
        MaterialDatabase::default()
    }

    /// Insert (replaces an existing entry with the same name).
    pub fn insert(&mut self, material: Material) {
        let name = material.name().to_string();
        self.entries.retain(|(n, _)| n != &name);
        self.entries.push((name, Arc::new(material)));
    }

    /// Look up by name (case-sensitive).
    pub fn get(&self, name: &str) -> Option<Arc<Material>> {
        self.entries
            .iter()
            .find(|(n, _)| n == name)
            .map(|(_, m)| Arc::clone(m))
    }

    /// Number of materials.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// True when empty.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Names of all materials, sorted.
    pub fn names(&self) -> Vec<&str> {
        let mut names: Vec<&str> = self.entries.iter().map(|(n, _)| n.as_str()).collect();
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
            db.insert(Material::new(
                name,
                disp,
                dn,
                crate::units::REFERENCE_TEMPERATURE,
            ));
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
            Tabulated::from_points(&[
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
            Tabulated::from_points(&[
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
