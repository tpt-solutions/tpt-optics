//! Glass catalogs: Zemax AGF-style parsing and an embedded reference set.
//!
//! The AGF subset covers `NM` (name), `ED` (extended data), `CD`
//! (Sellmeier `K1 L1 K2 L2 K3 L3` on one line), `TD` (thermal dn/dT) and
//! `CC` (concatenation) lines. The embedded [`reference_catalog`] carries a
//! small set of well-known glasses for tests and examples — it is an
//! *illustrative subset*, not a substitute for vendor catalogs.

use std::string::{String, ToString};
use std::vec::Vec;

use tpt_optics_core::dispersion::{Dispersion, Schott, Sellmeier};
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

impl core::error::Error for AgfError {}

/// Number of coefficients in a formula-1 (Schott) or formula-2 (Sellmeier) record.
const FORMULA_COEFFICIENTS: usize = 6;

/// One glass record under construction.
#[derive(Default, Clone, Debug)]
struct AgfGlass {
    name: String,
    /// Dispersion formula number from the `NM` line.
    formula: Option<u32>,
    /// `CD` coefficients in file order, collected across all `CD` lines.
    cd: Vec<f64>,
    dn_dtemp: f64,
}

/// Parse AGF text into a material database.
///
/// The `NM` line gives the dispersion formula number as its second field.
/// Formula 1 (Schott) and formula 2 (Sellmeier) are supported. Other formula
/// numbers are an error rather than being read as Sellmeier.
///
/// For formula 2 the `CD` coefficients are `K1 L1 K2 L2 K3 L3` (with `L` in
/// µm²). For formula 1 they are `A0 … A5`. A `CD` record may wrap across
/// several lines. Non-numeric or non-finite tokens are errors.
///
/// `TD` (thermal) and `CC` (concatenation) lines are recognised and skipped.
/// The temperature coefficient comes from the `ED` line.
pub fn parse_agf(text: &str) -> Result<MaterialDatabase, AgfError> {
    // A UTF-8 byte-order mark is not part of the first keyword.
    let text = text.trim_start_matches('\u{feff}');
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
                .filter(|v| v.is_finite())
                .ok_or(AgfError {
                    line: line_no,
                    message: std::format!("'{key}' expected a finite number at position {i}"),
                })
        };
        match key.as_str() {
            "NM" => {
                if let Some(g) = current.take() {
                    finish_glass(&mut db, &g, line_no)?;
                }
                let name = args.first().copied().unwrap_or("unnamed").to_string();
                let formula = match args.get(1) {
                    None => None,
                    Some(s) => Some(s.parse::<u32>().map_err(|_| AgfError {
                        line: line_no,
                        message: std::format!("glass '{name}' has a bad dispersion formula '{s}'"),
                    })?),
                };
                current = Some(AgfGlass {
                    name,
                    formula,
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
                for s in &args {
                    let v = s
                        .parse::<f64>()
                        .ok()
                        .filter(|v| v.is_finite())
                        .ok_or_else(|| AgfError {
                            line: line_no,
                            message: std::format!("'CD' value '{s}' is not a finite number"),
                        })?;
                    g.cd.push(v);
                }
                if g.cd.len() > FORMULA_COEFFICIENTS {
                    return Err(AgfError {
                        line: line_no,
                        message: std::format!(
                            "glass '{}' has more than {FORMULA_COEFFICIENTS} CD coefficients",
                            g.name
                        ),
                    });
                }
            }
            _ => {} // GC (comments), TD, CC, OD, etc. are recognised and skipped
        }
    }
    if let Some(g) = current.take() {
        finish_glass(&mut db, &g, text.lines().count())?;
    }
    Ok(db)
}

fn finish_glass(db: &mut MaterialDatabase, g: &AgfGlass, line: usize) -> Result<(), AgfError> {
    let fail = |message: std::string::String| AgfError { line, message };
    if g.name.is_empty() {
        return Err(fail("glass without a name".into()));
    }
    let formula = g
        .formula
        .ok_or_else(|| fail(std::format!("glass '{}' has no dispersion formula", g.name)))?;
    if formula != 1 && formula != 2 {
        return Err(fail(std::format!(
            "glass '{}' uses dispersion formula {formula}; only 1 (Schott) and 2 (Sellmeier) are supported",
            g.name
        )));
    }
    if g.cd.iter().all(|v| *v == 0.0) {
        return Err(fail(std::format!(
            "glass '{}' has no dispersion data",
            g.name
        )));
    }
    if g.cd.len() != FORMULA_COEFFICIENTS {
        return Err(fail(std::format!(
            "glass '{}' has {} CD coefficients; formula {formula} takes {FORMULA_COEFFICIENTS}",
            g.name,
            g.cd.len()
        )));
    }
    let c = &g.cd;
    let dispersion: std::boxed::Box<dyn Dispersion> = if formula == 1 {
        std::boxed::Box::new(Schott::new([c[0], c[1], c[2], c[3], c[4], c[5]]))
    } else {
        std::boxed::Box::new(Sellmeier::new([c[0], c[2], c[4]], [c[1], c[3], c[5]]))
    };
    db.insert(Material::new(
        g.name.clone(),
        dispersion,
        g.dn_dtemp,
        REFERENCE_TEMPERATURE,
    ));
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
NM N-BK7 2 0 1.51680 64.17 0 3 3 0
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

    #[test]
    fn split_cd_lines_match_single_line() {
        // CD coefficients may wrap onto a second line; the result must equal
        // the single-line form, not drop or shift values.
        let single = parse_agf("NM X 2 0 1.5 60\nCD 1.0 0.01 0.2 0.02 1.0 100.0\n").unwrap();
        let split = parse_agf("NM X 2 0 1.5 60\nCD 1.0 0.01 0.2 0.02\nCD 1.0 100.0\n").unwrap();
        let t = Kelvin::from_celsius(20.0);
        let a = single.get("X").unwrap().index_at(wl(587.6), t).n;
        let b = split.get("X").unwrap().index_at(wl(587.6), t).n;
        assert!((a - b).abs() < 1e-12, "single {a} vs split {b}");
    }

    #[test]
    fn extra_cd_coefficients_are_an_error_not_a_panic() {
        // A second CD line after six coefficients used to index past the
        // three Sellmeier slots and panic.
        let text = "NM X 2 0 1.5 60\nCD 1 0.01 0.2 0.02 1 100\nCD 1 0.01 0.2 0.02 1 100\n";
        assert!(parse_agf(text).is_err());
    }

    #[test]
    fn non_numeric_cd_token_is_an_error() {
        // A bad token used to be skipped, which moved every later value into
        // the wrong slot.
        let text = "NM X 2 0 1.5 60\nCD 1 0.01 abc 0.02 1 100\n";
        let err = parse_agf(text).unwrap_err();
        assert_eq!(err.line, 2, "{err}");
    }

    #[test]
    fn non_finite_coefficients_are_rejected() {
        for bad in ["nan", "inf", "-inf"] {
            let text = std::format!("NM X 2 0 1.5 60\nCD {bad} 0.01 0.2 0.02 1 100\n");
            assert!(parse_agf(&text).is_err(), "accepted {bad}");
        }
    }

    #[test]
    fn unsupported_dispersion_formula_is_an_error() {
        // Formula 4 is not implemented. It must not be read as Sellmeier.
        let text = "NM X 4 0 1.5 60\nCD 1 0.01 0.2 0.02 1 100\n";
        let err = parse_agf(text).unwrap_err();
        assert!(err.message.contains("formula"), "{err}");
    }

    #[test]
    fn schott_formula_one_uses_polynomial_form() {
        // Formula 1 stores A0..A5 of n² = A0 + A1 λ² + A2 λ⁻² + ... (λ in µm).
        let db = parse_agf("NM S 1 0 1.5 60\nCD 2.3 0.01 0 0 0 0\n").unwrap();
        let t = Kelvin::from_celsius(20.0);
        let n = db.get("S").unwrap().index_at(wl(500.0), t).n;
        let expected = (2.3 + 0.01 * 0.25_f64).sqrt();
        assert!((n - expected).abs() < 1e-12, "n = {n}, expected {expected}");
    }
}
