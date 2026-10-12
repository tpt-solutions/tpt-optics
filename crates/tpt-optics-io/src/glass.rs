//! Glass catalogs: Zemax AGF-style parsing and an embedded reference set.
//!
//! The AGF subset covers `NM` (name), `ED` (extended data), `CD`
//! (Sellmeier `K1 L1 K2 L2 K3 L3` on one line), `TD` (thermal dn/dT) and
//! `CC` (concatenation) lines. The embedded [`reference_catalog`] carries a
//! small set of well-known glasses for tests and examples — it is an
//! *illustrative subset*, not a substitute for vendor catalogs.

use std::string::{String, ToString};
use std::sync::Arc;
use std::vec::Vec;

use tpt_optics_core::dispersion::Sellmeier;
use tpt_optics_core::material::{Material, MaterialDatabase};
use tpt_optics_core::units::REFERENCE_TEMPERATURE;

/// A parse error with a line number.
#[derive(Debug, Clone, PartialEq)]
pub struct AgfError {
    /// 1-based line number.
    pub line: usize,
    /// Description.
    pub message: String,
}

impl core::fmt::Display for AgfError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "AGF line {}: {}", self.line, self.message)
    }
}

#[cfg(feature = "std")]
impl std::error::Error for AgfError {}

/// One glass record under construction.
#[derive(Default, Clone, Debug)]
struct AgfGlass {
    name: String,
    b: [f64; 3],
    c: [f64; 3],
    next_cd_slot: usize,
    dn_dtemp: f64,
}

/// Parse AGF text into a material database.
///
/// A `CD` line carries the Sellmeier coefficients `K1 L1 K2 L2 K3 L3` in
/// that order, with `L` in µm².
pub fn parse_agf(text: &str) -> Result<MaterialDatabase, AgfError> {
    let mut db = MaterialDatabase::new();
    let mut current: Option<AgfGlass> = None;
    for (idx, raw) in text.lines().enumerate() {
        let line_no = idx + 1;
        let line = raw.trim();
        if line.is_empty() || line.starts_with('!') {
            continue;
        }
        let mut tokens = line.split_whitespace();
        let key = tokens.next().unwrap_or("").to_ascii_uppercase();
        let args: Vec<&str> = tokens.collect();
        let num = |i: usize| -> Result<f64, AgfError> {
            args.get(i)
                .and_then(|s| s.parse::<f64>().ok())
                .ok_or(AgfError {
                    line: line_no,
                    message: std::format!("'{key}' expected a number at position {i}"),
                })
        };
        match key.as_str() {
            "NM" => {
                if let Some(g) = current.take() {
                    finish_glass(&mut db, &g, line_no)?;
                }
                let name = args.first().copied().unwrap_or("unnamed").to_string();
                current = Some(AgfGlass {
                    name,
                    ..AgfGlass::default()
                });
            }
            "ED" => {
                // Extended data: position 6 (index 5) is dn/dT in some
                // catalogs; tolerated if absent.
                if let Ok(v) = num(5) {
                    if let Some(g) = current.as_mut() {
                        g.dn_dtemp = v;
                    }
                }
            }
            "CD" => {
                let Some(g) = current.as_mut() else {
                    return Err(AgfError {
                        line: line_no,
                        message: "CD before NM".into(),
                    });
                };
                // Sellmeier K1 L1 K2 L2 K3 L3: each (K, L) pair is one term.
                let vals: Vec<f64> = args.iter().filter_map(|s| s.parse::<f64>().ok()).collect();
                for pair in vals.chunks_exact(2).take(3) {
                    g.b[g.next_cd_slot] = pair[0];
                    g.c[g.next_cd_slot] = pair[1];
                    g.next_cd_slot += 1;
                }
            }
            _ => {} // GC (comments), OD, etc. tolerated
        }
    }
    if let Some(g) = current.take() {
        finish_glass(&mut db, &g, text.lines().count())?;
    }
    Ok(db)
}

fn finish_glass(db: &mut MaterialDatabase, g: &AgfGlass, line: usize) -> Result<(), AgfError> {
    if g.name.is_empty() {
        return Err(AgfError {
            line,
            message: "glass without a name".into(),
        });
    }
    let material = if g.b.iter().any(|v| *v != 0.0) {
        Material::new(
            g.name.clone(),
            std::boxed::Box::new(Sellmeier::new(g.b, g.c)),
            g.dn_dtemp,
            REFERENCE_TEMPERATURE,
        )
    } else {
        return Err(AgfError {
            line,
            message: std::format!("glass '{}' has no dispersion data", g.name),
        });
    };
    db.insert(material);
    Ok(())
}

/// Embedded reference glass catalog: a handful of well-known glasses with
/// published Sellmeier coefficients (BK7, fused silica, SF11, F2, diamond,
/// water, silicon, gold — see `tpt-optics-core` presets), plus the
/// hyphenated vendor spellings used by Zemax files.
#[must_use]
pub fn reference_catalog() -> MaterialDatabase {
    let mut db = MaterialDatabase::presets();
    // Hyphenated aliases of the presets for .zmx GLAS names.
    for alias in ["N-BK7", "BK7"] {
        db.insert(Material::new(
            alias,
            std::boxed::Box::new(Sellmeier::bk7()),
            3.0e-6,
            REFERENCE_TEMPERATURE,
        ));
    }
    // Extra lens glasses commonly referenced in tests.
    db.insert(Material::new(
        "N-SF11",
        std::boxed::Box::new(Sellmeier::new(
            [1.737_596_95, 0.313_747_346, 1.898_781_01],
            [0.013_188_707, 0.062_306_814_2, 155.236_29],
        )),
        7.0e-6,
        REFERENCE_TEMPERATURE,
    ));
    db.insert(Material::new(
        "F2",
        std::boxed::Box::new(Sellmeier::new(
            [1.345_333_59, 0.209_073_176, 0.937_357_162],
            [0.009_977_438_71, 0.047_045_076_7, 111.886_764],
        )),
        6.0e-6,
        REFERENCE_TEMPERATURE,
    ));
    let _ = Arc::new(()); // silence unused-import lint when features drift
    db
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_optics_core::index::IndexAt;
    use tpt_optics_core::{Kelvin, Wavelength};

    fn wl(nm: f64) -> Wavelength {
        Wavelength::from_nm(nm)
    }

    const SAMPLE_AGF: &str = "\
! reference glasses
CC
NM N-BK7 1 2 3.0 2.0 0.0 0.0 0
ED 1.51680 64.17 0.0 0.0 0.0 3.0E-6 0.0 0
CD 1.03961212 0.00600069867 0.231792344 0.0200179144 1.01046945 103.560653
TD 0.0 0.0 0.0 3.0e-6
";

    #[test]
    fn parses_sellmeier_glass() {
        let db = parse_agf(SAMPLE_AGF).unwrap();
        let bk7 = db.get("N-BK7").expect("N-BK7");
        let n = bk7.index_at(wl(587.6), Kelvin::from_celsius(20.0)).n;
        assert!((n - 1.5168).abs() < 1e-3, "n = {n}");
    }

    #[test]
    fn bk7_single_line_cd_matches_published_indices() {
        // Schott N-BK7: n_F (486.1 nm) = 1.52238, n_C (656.3 nm) = 1.51432.
        let db = parse_agf(SAMPLE_AGF).unwrap();
        let bk7 = db.get("N-BK7").expect("N-BK7");
        let t = Kelvin::from_celsius(20.0);
        let nf = bk7.index_at(wl(486.1), t).n;
        let nc = bk7.index_at(wl(656.3), t).n;
        assert!((nf - 1.52238).abs() < 2e-4, "n_F = {nf}");
        assert!((nc - 1.51432).abs() < 2e-4, "n_C = {nc}");
    }

    #[test]
    fn reference_catalog_bk7_and_sf11() {
        let db = reference_catalog();
        let bk7 = db.get("bk7").expect("bk7");
        assert!((bk7.index_at(wl(587.6), Kelvin::from_celsius(20.0)).n - 1.5168).abs() < 2e-4);
        let sf11 = db.get("N-SF11").expect("N-SF11");
        // n_d ≈ 1.78472
        assert!(
            (sf11.index_at(wl(587.6), Kelvin::from_celsius(20.0)).n - 1.7847).abs() < 2e-3,
            "n = {}",
            sf11.index_at(wl(587.6), Kelvin::from_celsius(20.0)).n
        );
    }

    #[test]
    fn glass_without_data_rejected() {
        let err = parse_agf("NM EMPTY 1 1 1 1 0 0 0\n").unwrap_err();
        assert!(err.message.contains("no dispersion"));
    }
}
