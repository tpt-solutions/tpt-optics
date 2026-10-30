//! Zemax `.zmx` text-format sequential files (pragmatic subset).
//!
//! Supported keywords per surface: `SURF`, `CURV` (curvature + asphere
//! terms), `DISZ` (thickness; `INFINITY` supported), `GLAS` (name + d/Abbe
//! + index hint), `CONI`, `DIAM`, `STOP`, `NAME`, `TYPE STANDARD`. Header: `MODE SEQ`.
//!
//! Keywords that would change the geometry but are not modelled are errors,
//! not silent drops: `TYPE` other than `STANDARD`, and `PARM` with a non-zero
//! value. Other unknown keywords (`FNUM`, `WAVL`, ...) are ignored. `STOP`
//! is kept in the document but does not change the trace. The object-space
//! `DISZ` is kept on surface 0 but is not used, because rays are launched at
//! the first traced surface.

use std::string::{String, ToString};
use std::vec::Vec;

use tpt_optics_core::{index::IndexAt, material::MaterialDatabase, Kelvin, Wavelength};
use tpt_optics_geometric::sequential::{Medium, SeqSurface, SequentialSystem, SurfaceAction};
use tpt_optics_geometric::surface::Asphere;

/// A parse or write error with a line number.
#[derive(Debug, Clone, PartialEq)]
pub struct ZmxError {
    /// 1-based line number (0 for structural errors).
    pub line: usize,
    /// Description.
    pub message: String,
}

impl core::fmt::Display for ZmxError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, ".zmx line {}: {}", self.line, self.message)
    }
}

impl core::error::Error for ZmxError {}

/// One parsed surface record.
#[derive(Clone, Debug, PartialEq)]
pub struct ZmxSurface {
    /// Surface name (may be empty).
    pub name: String,
    /// Curvature (1/m). The file stores CURV in 1/mm (lens units).
    pub curvature: f64,
    /// Even asphere terms following CURV, stored in lens units (mm) as in
    /// the file; [`to_sequential`] converts them to metres.
    pub asphere: [f64; 8],
    /// Conic constant.
    pub conic: f64,
    /// Thickness to the next surface (meters); the last surface's value is
    /// ignored on export.
    pub thickness: f64,
    /// Glass name (empty = air/mirror per action).
    pub glass: String,
    /// Semi-diameter (meters).
    pub semi_diameter: f64,
    /// Stop flag.
    pub stop: bool,
}

impl Default for ZmxSurface {
    fn default() -> Self {
        ZmxSurface {
            name: String::new(),
            curvature: 0.0,
            asphere: [0.0; 8],
            conic: 0.0,
            thickness: 0.0,
            glass: String::new(),
            semi_diameter: 0.0,
            stop: false,
        }
    }
}

/// A parsed sequential .zmx document.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ZmxDocument {
    /// System name.
    pub name: String,
    /// Surfaces in order (object surface first, image surface last).
    pub surfaces: Vec<ZmxSurface>,
}

/// Parse a text `.zmx` document.
pub fn parse(text: &str) -> Result<ZmxDocument, ZmxError> {
    // A UTF-8 byte-order mark is not part of the first keyword.
    let text = text.trim_start_matches('\u{feff}');
    let mut doc = ZmxDocument::default();
    let mut current: Option<ZmxSurface> = None;
    for (idx, raw) in text.lines().enumerate() {
        let line_no = idx + 1;
        let line = raw.trim();
        if line.is_empty() || line.starts_with('!') || line.starts_with('#') {
            continue;
        }
        let mut tokens = line.split_whitespace();
        let key = tokens.next().unwrap_or("").to_ascii_uppercase();
        let args: Vec<&str> = tokens.collect();
        match key.as_str() {
            "MODE" => {
                if args
                    .first()
                    .copied()
                    .map(str::to_ascii_uppercase)
                    .as_deref()
                    != Some("SEQ")
                {
                    return Err(ZmxError {
                        line: line_no,
                        message: "only MODE SEQ files are supported".into(),
                    });
                }
            }
            "NAME" => {
                // Inside a SURF block the name belongs to that surface.
                let name = args.join(" ");
                match current.as_mut() {
                    Some(c) => c.name = name,
                    None => doc.name = name,
                }
            }
            "SURF" => {
                if let Some(s) = current.take() {
                    doc.surfaces.push(s);
                }
                current = Some(ZmxSurface::default());
            }
            "CURV" => {
                let c = current.as_mut().ok_or(ZmxError {
                    line: line_no,
                    message: "CURV before SURF".into(),
                })?;
                // CURV is in 1/mm; store 1/m.
                c.curvature = parse_f64(args.first(), line_no)? * 1e3;
                for (i, v) in args.iter().skip(1).take(8).enumerate() {
                    c.asphere[i] = parse_f64(Some(v), line_no)?;
                }
            }
            "CONI" => {
                let c = current.as_mut().ok_or(ZmxError {
                    line: line_no,
                    message: "CONI before SURF".into(),
                })?;
                c.conic = parse_f64(args.first(), line_no)?;
            }
            "DISZ" => {
                let c = current.as_mut().ok_or(ZmxError {
                    line: line_no,
                    message: "DISZ before SURF".into(),
                })?;
                c.thickness = if args
                    .first()
                    .is_some_and(|s| s.eq_ignore_ascii_case("INFINITY"))
                {
                    f64::INFINITY
                } else {
                    parse_f64(args.first(), line_no)? * 1e-3 // mm → m
                };
            }
            "GLAS" => {
                let c = current.as_mut().ok_or(ZmxError {
                    line: line_no,
                    message: "GLAS before SURF".into(),
                })?;
                if let Some(name) = args.first() {
                    if !name.eq_ignore_ascii_case("AIR") {
                        c.glass = (*name).to_string();
                    }
                }
            }
            "DIAM" => {
                let c = current.as_mut().ok_or(ZmxError {
                    line: line_no,
                    message: "DIAM before SURF".into(),
                })?;
                c.semi_diameter = parse_f64(args.first(), line_no)? * 1e-3; // mm → m
            }
            "STOP" => {
                let c = current.as_mut().ok_or(ZmxError {
                    line: line_no,
                    message: "STOP before SURF".into(),
                })?;
                c.stop = true;
            }
            "TYPE" => {
                // Only standard surfaces are modelled. Other types (EVENASPH,
                // ...) keep their shape in PARM lines, which are not read.
                let kind = args.first().copied().unwrap_or("");
                if !kind.eq_ignore_ascii_case("STANDARD") {
                    return Err(ZmxError {
                        line: line_no,
                        message: std::format!(
                            "surface type '{kind}' is not supported (only STANDARD)"
                        ),
                    });
                }
            }
            "PARM" => {
                // A non-zero parameter belongs to a non-standard surface type
                // and would be dropped, so it is an error.
                let value = parse_f64(args.get(1), line_no)?;
                if value != 0.0 {
                    return Err(ZmxError {
                        line: line_no,
                        message: "non-zero surface parameters (PARM) are not supported".into(),
                    });
                }
            }
            _ => {} // tolerated unknown keywords (FNUM, WAVL, ...)
        }
    }
    if let Some(s) = current.take() {
        doc.surfaces.push(s);
    }
    Ok(doc)
}

/// Parse a finite number. `nan` and `inf` parse as `f64` but would feed
/// non-finite curvature or thickness into the tracer, so they are rejected.
fn parse_f64(v: Option<&&str>, line: usize) -> Result<f64, ZmxError> {
    v.and_then(|s| s.parse::<f64>().ok())
        .filter(|x| x.is_finite())
        .ok_or(ZmxError {
            line,
            message: "expected a number".into(),
        })
}

/// Serialize a document back to text form.
pub fn write(doc: &ZmxDocument) -> String {
    use std::fmt::Write as _;
    let mut out = String::new();
    out.push_str("MODE SEQ\n");
    if !doc.name.is_empty() {
        let _ = writeln!(out, "NAME {}", doc.name);
    }
    for s in &doc.surfaces {
        out.push_str("SURF\n");
        if !s.name.is_empty() {
            let _ = writeln!(out, "NAME {}", s.name);
        }
        if s.stop {
            out.push_str("STOP\n");
        }
        // CURV: curvature in 1/mm, then all eight asphere terms in lens units.
        // Shortest exact form (`{:e}`), so coefficients of any size round-trip.
        let mut curv = std::format!("CURV {:.12}", s.curvature / 1e3);
        for a in &s.asphere {
            let _ = write!(curv, " {a:e}");
        }
        let _ = writeln!(out, "{curv}");
        let _ = writeln!(out, "CONI {:.12}", s.conic);
        if !s.glass.is_empty() {
            let _ = writeln!(out, "GLAS {} 0 0 0 0", s.glass);
        }
        if s.semi_diameter > 0.0 {
            let _ = writeln!(out, "DIAM {:.6}", s.semi_diameter * 1e3);
        }
        let th = if s.thickness.is_infinite() {
            "INFINITY".to_string()
        } else {
            std::format!("{:.9}", s.thickness * 1e3)
        };
        let _ = writeln!(out, "DISZ {th}");
    }
    out
}

/// Factor converting the `i`-th (0-based) even-asphere coefficient from lens
/// units to metres. A term `a·r^{2(i+1)}` with `z` and `r` in mm becomes
/// `a·10^{6(i+1)−3}·r_m^{2(i+1)}` in metres, so the factor is `10^{6(i+1)−3}`.
fn asphere_mm_to_m(i: usize) -> f64 {
    10f64.powi(6 * (i as i32 + 1) - 3)
}

/// Convert a parsed document into a sequential system, resolving glasses
/// against `catalog` at the reference wavelength.
///
/// Surfaces carrying the `MIRROR` name convention (glass name starting
/// with `MIRROR`) are exported as reflecting.
pub fn to_sequential(
    doc: &ZmxDocument,
    catalog: &MaterialDatabase,
    wl: Wavelength,
    temperature: Kelvin,
) -> Result<SequentialSystem, ZmxError> {
    let mut sys = SequentialSystem::new();
    if doc.surfaces.is_empty() {
        return Err(ZmxError {
            line: 0,
            message: "document has no surfaces".into(),
        });
    }
    // Object space: air. Zemax files always begin with SURF 0 (the object
    // surface, thickness = distance to the first surface); it is not a
    // traced surface and is skipped here.
    sys.media[0] = Medium::Const(1.0);
    if doc.surfaces.is_empty() {
        return Err(ZmxError {
            line: 0,
            message: "document has no surfaces".into(),
        });
    }
    // The medium the ray is travelling in, starting in object-space air.
    let mut travelling = Medium::Const(1.0);
    for (i, s) in doc.surfaces.iter().enumerate().skip(1) {
        let is_last = i + 1 == doc.surfaces.len();
        let is_mirror = s.glass.to_ascii_uppercase().starts_with("MIRROR");
        let medium_after = if is_mirror {
            // A mirror reflects: the ray stays in the medium it was in.
            travelling.clone()
        } else if is_last || s.glass.is_empty() {
            Medium::Const(1.0)
        } else if let Some(mat) = catalog.get(&s.glass) {
            Medium::Material(mat)
        } else {
            // Unknown glass: try to parse a plain index ("1.5168") or fail.
            if let Ok(v) = s.glass.parse::<f64>() {
                Medium::Const(v)
            } else {
                let hint = match catalog.closest_name(&s.glass, 2) {
                    Some(near) => std::format!("; did you mean '{near}'?"),
                    None => String::new(),
                };
                return Err(ZmxError {
                    line: 0,
                    message: std::format!("glass '{}' not in catalog{hint}", s.glass),
                });
            }
        };
        let action = if is_mirror {
            SurfaceAction::Reflect
        } else {
            SurfaceAction::Refract
        };
        // Asphere terms are in lens units (mm); convert each to metres.
        let mut coeffs = [0.0; 8];
        for (ci, z) in s.asphere.iter().enumerate() {
            coeffs[ci] = z * asphere_mm_to_m(ci);
        }
        let surface = SeqSurface {
            action,
            asphere: Asphere {
                curvature: s.curvature,
                conic: s.conic,
                coeffs,
            },
            semi_diameter: if s.semi_diameter > 0.0 {
                s.semi_diameter
            } else {
                f64::INFINITY
            },
            phase: None,
        };
        let thickness = if is_last { 0.0 } else { s.thickness };
        travelling = medium_after.clone();
        sys.push(surface, thickness, medium_after);
    }
    let _ = (wl, temperature);
    Ok(sys)
}

/// Build a document from a sequential system (round-trip support),
/// tagging constant-index media with their index and material media with
/// their names.
pub fn from_sequential(system: &SequentialSystem) -> ZmxDocument {
    let mut doc = ZmxDocument::default();
    // Zemax files start with the object surface, which `to_sequential`
    // skips; emit one so a round trip keeps every traced surface.
    doc.surfaces.push(ZmxSurface {
        thickness: f64::INFINITY,
        ..ZmxSurface::default()
    });
    let n = system.surfaces.len();
    for (i, s) in system.surfaces.iter().enumerate() {
        let mut z = ZmxSurface {
            curvature: s.asphere.curvature,
            conic: s.asphere.conic,
            ..ZmxSurface::default()
        };
        for (ci, a) in s.asphere.coeffs.iter().enumerate() {
            z.asphere[ci] = a / asphere_mm_to_m(ci);
        }
        z.semi_diameter = if s.semi_diameter.is_finite() {
            s.semi_diameter
        } else {
            0.0
        };
        if i + 1 < n {
            z.thickness = system.thickness[i];
        }
        match (&system.media[i + 1], s.action) {
            // A mirror is a mirror whatever medium it sits in, so this arm
            // comes before the vacuum check.
            (_, SurfaceAction::Reflect) => {
                z.glass = "MIRROR".into();
            }
            (_, SurfaceAction::Diffract) => {
                z.glass = String::new();
            }
            (Medium::Vacuum, _) => {}
            (Medium::Const(v), SurfaceAction::Refract) => {
                // `{v}` is the shortest exact form, so the index survives a round trip.
                if *v != 1.0 {
                    z.glass = std::format!("{v}");
                }
            }
            (Medium::Material(m), SurfaceAction::Refract) => {
                z.glass = m.name().to_string();
            }
        }
        doc.surfaces.push(z);
    }
    doc
}

/// Reference: index of a `GLAS` name via the embedded catalog, falling back
/// to a plain-number parse.
pub fn glass_index(name: &str, wl: Wavelength) -> Option<f64> {
    if let Some(mat) = crate::glass::reference_catalog().get(name) {
        return Some(mat.index_at(wl, Kelvin::from_celsius(20.0)).n);
    }
    name.parse().ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = "\
MODE SEQ
NAME singlet test
! a comment
SURF
CURV 0 0 0 0 0 0 0 0 0
DISZ INFINITY
SURF
STOP
DIAM 12.5
CURV 0.01 0 0 0 0 0 0 0 0
GLAS N-BK7 0 0 0 0
CONI 0.0
DISZ 5.0
SURF
DIAM 12.0
CURV -0.01 0 0 0 0 0 0 0 0
GLAS AIR
DISZ 95.0
SURF
CURV 0 0 0 0 0 0 0 0 0
DISZ 0
";

    #[test]
    fn parses_sample() {
        let doc = parse(SAMPLE).unwrap();
        assert_eq!(doc.name, "singlet test");
        assert_eq!(doc.surfaces.len(), 4);
        assert_eq!(doc.surfaces[1].curvature, 10.0);
        assert_eq!(doc.surfaces[1].thickness, 5e-3);
        assert_eq!(doc.surfaces[1].glass, "N-BK7");
        assert!(doc.surfaces[1].stop);
        assert!(doc.surfaces[0].thickness.is_infinite());
        assert_eq!(doc.surfaces[2].glass, "");
    }

    #[test]
    fn round_trip_preserves() {
        let doc = parse(SAMPLE).unwrap();
        let text = write(&doc);
        let doc2 = parse(&text).unwrap();
        assert_eq!(doc, doc2);
    }

    #[test]
    fn rejects_garbage_numbers() {
        let bad = "MODE SEQ\nSURF\nCURV abc\n";
        let err = parse(bad).unwrap_err();
        assert_eq!(err.line, 3);
    }

    #[test]
    fn rejects_non_sequential_mode() {
        let err = parse("MODE NSC\n").unwrap_err();
        assert!(err.message.contains("SEQ"));
    }

    #[test]
    fn converts_to_sequential_and_traces() {
        let doc = parse(SAMPLE).unwrap();
        let catalog = crate::glass::reference_catalog();
        let sys = to_sequential(
            &doc,
            &catalog,
            Wavelength::from_nm(587.6),
            Kelvin::from_celsius(20.0),
        )
        .unwrap();
        // Object surface skipped: three traced surfaces.
        assert_eq!(sys.surfaces.len(), 3);
        // BK7 resolved from the catalog (medium after the first traced
        // surface R1).
        let n = match &sys.media[1] {
            Medium::Material(m) => {
                m.index_at(Wavelength::from_nm(587.6), Kelvin::from_celsius(20.0))
                    .n
            }
            _ => panic!("expected material"),
        };
        assert!((n - 1.5168).abs() < 1e-3);
        // Traces to a sensible focus.
        let ctx = tpt_optics_geometric::sequential::TraceContext {
            wavelength: Wavelength::from_nm(587.6),
            temperature: Kelvin::from_celsius(20.0),
        };
        let ray = tpt_optics_geometric::prelude::Ray::new(
            tpt_optics_core::Vec3::new(1e-3, 0.0, -1e-6),
            tpt_optics_core::Vec3::new(0.0, 0.0, 1.0),
        );
        let res = sys.trace(&ray, ctx);
        assert!(!res.blocked);
        assert!(res.exit_direction.x < 0.0);
    }

    #[test]
    fn unknown_glass_is_an_error() {
        let bad = "MODE SEQ\nSURF\nCURV 0\nDISZ 5\nSURF\nGLAS UNOBTAINIUM 0 0 0 0\nDISZ 5\nSURF\n";
        let doc = parse(bad).unwrap();
        let catalog = crate::glass::reference_catalog();
        let err = to_sequential(
            &doc,
            &catalog,
            Wavelength::from_nm(550.0),
            Kelvin::from_celsius(20.0),
        );
        assert!(err.is_err());
    }

    #[test]
    fn from_sequential_round_trip() {
        let doc = parse(SAMPLE).unwrap();
        let catalog = crate::glass::reference_catalog();
        let sys = to_sequential(
            &doc,
            &catalog,
            Wavelength::from_nm(587.6),
            Kelvin::from_celsius(20.0),
        )
        .unwrap();
        let doc2 = from_sequential(&sys);
        // Object surface first, then the three traced surfaces.
        assert_eq!(doc2.surfaces.len(), 4);
        assert!(doc2.surfaces[0].thickness.is_infinite());
        assert!((doc2.surfaces[1].curvature - 10.0).abs() < 1e-15);
        assert_eq!(doc2.surfaces[1].glass, "N-BK7");
        // Text round trip keeps every traced surface and its curvature.
        let sys2 = to_sequential(
            &parse(&write(&doc2)).unwrap(),
            &catalog,
            Wavelength::from_nm(587.6),
            Kelvin::from_celsius(20.0),
        )
        .unwrap();
        assert_eq!(sys2.surfaces.len(), sys.surfaces.len());
        assert!((sys2.surfaces[0].asphere.curvature - 10.0).abs() < 1e-12);
    }

    #[test]
    fn curv_is_read_in_inverse_millimetres() {
        // CURV 0.02 mm⁻¹ is R = 50 mm, i.e. 20 m⁻¹.
        let doc = parse("MODE SEQ\nSURF\nCURV 0.02\n").unwrap();
        assert!((doc.surfaces[0].curvature - 20.0).abs() < 1e-12);
    }

    #[test]
    fn asphere_factor_matches_mm_to_m_expansion() {
        // z = a·r⁴ with z, r in mm: z_m = 1e-3·a·(1e3·r_m)⁴ ⇒ a_m = a·1e9.
        assert_eq!(asphere_mm_to_m(0), 1e3);
        assert_eq!(asphere_mm_to_m(1), 1e9);
        let a_mm = 2.0e-6_f64;
        let r_m = 1.0e-3_f64;
        let z_m = 1e-3 * a_mm * (1e3 * r_m).powi(4);
        let z_via_metres = a_mm * asphere_mm_to_m(1) * r_m.powi(4);
        assert!(
            (z_m - z_via_metres).abs() < 1e-18,
            "{z_m} vs {z_via_metres}"
        );
    }

    #[test]
    fn surface_name_inside_surf_block() {
        let doc = parse("MODE SEQ\nNAME lens\nSURF\nNAME front\nCURV 0\nSURF\n").unwrap();
        assert_eq!(doc.name, "lens");
        assert_eq!(doc.surfaces[0].name, "front");
    }
}

#[cfg(test)]
mod mirror_medium_tests {
    use super::*;

    #[test]
    fn mirror_keeps_the_ray_in_its_medium() {
        // Object air → BK7 lens surface → mirror (reflects inside glass)
        // → image. The medium after the mirror must still be BK7.
        let text = "MODE SEQ\nSURF\nDISZ INFINITY\nSURF\nGLAS N-BK7 0 0 0 0\nDISZ 5\nSURF\nGLAS MIRROR\nDISZ 5\nSURF\nDISZ 0\n";
        let doc = parse(text).unwrap();
        let catalog = crate::glass::reference_catalog();
        let sys = to_sequential(
            &doc,
            &catalog,
            Wavelength::from_nm(587.6),
            Kelvin::from_celsius(20.0),
        )
        .unwrap();
        assert!(matches!(&sys.media[1], Medium::Material(m) if m.name() == "N-BK7"));
        assert!(matches!(&sys.media[2], Medium::Material(m) if m.name() == "N-BK7"));
    }
}

#[cfg(test)]
mod unsupported_and_precision_tests {
    use super::*;

    #[test]
    fn aspheric_surface_type_is_an_error() {
        // An EVENASPH surface keeps its shape in PARM lines. Loading it as a
        // sphere would give the wrong geometry, so it must be rejected.
        let text = "MODE SEQ\nSURF\nTYPE EVENASPH\nCURV 0.01\nPARM 1 0.5\n";
        let err = parse(text).unwrap_err();
        assert_eq!(err.line, 3, "{err}");
    }

    #[test]
    fn standard_surface_type_is_accepted() {
        let doc = parse("MODE SEQ\nSURF\nTYPE STANDARD\nCURV 0.01\n").unwrap();
        assert_eq!(doc.surfaces.len(), 1);
    }

    #[test]
    fn nonzero_parameter_is_an_error() {
        let err = parse("MODE SEQ\nSURF\nPARM 2 0.25\n").unwrap_err();
        assert_eq!(err.line, 3, "{err}");
    }

    #[test]
    fn asphere_coefficients_survive_write_and_parse_exactly() {
        // Writing with {:.12} rounded these to 12 decimals.
        let mut doc = ZmxDocument::default();
        let mut s = ZmxSurface {
            curvature: 10.0,
            ..ZmxSurface::default()
        };
        s.asphere = [
            1.234_567_890_123_45e-9,
            -3.141_592_653_589_79e-6,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            2.5e-14,
        ];
        doc.surfaces.push(s);
        let doc2 = parse(&write(&doc)).unwrap();
        assert_eq!(doc2.surfaces[0].asphere, doc.surfaces[0].asphere);
    }
}

#[cfg(test)]
mod export_and_bom_tests {
    use super::*;

    #[test]
    fn written_text_has_no_blank_lines() {
        // Each record used to end with an extra newline, which gave a blank
        // line after every keyword.
        let doc = parse("MODE SEQ\nNAME t\nSURF\nSTOP\nCURV 0.01\nDISZ 5\n").unwrap();
        let text = write(&doc);
        assert!(!text.contains("\n\n"), "{text:?}");
        assert_eq!(parse(&text).unwrap(), doc);
    }

    #[test]
    fn mirror_in_vacuum_is_exported_as_mirror() {
        // A mirror in air used to lose its MIRROR tag on export.
        let mut sys = SequentialSystem::new();
        sys.push(SeqSurface::mirror(1.0 / 0.5, 20e-3), 0.0, Medium::Vacuum);
        let doc = from_sequential(&sys);
        assert_eq!(doc.surfaces[1].glass, "MIRROR");
    }

    #[test]
    fn index_hint_survives_export_exactly() {
        // The index was written with six decimals, which changed the value.
        let n = 1.516_812_345_678_9;
        let mut sys = SequentialSystem::new();
        sys.push(SeqSurface::flat(5e-3), 0.0, Medium::Const(n));
        let doc = from_sequential(&sys);
        assert_eq!(doc.surfaces[1].glass.parse::<f64>().unwrap(), n);
    }

    #[test]
    fn utf8_bom_is_ignored() {
        let doc = parse("\u{feff}MODE SEQ\nSURF\nCURV 0.01\n").unwrap();
        assert_eq!(doc.surfaces.len(), 1);
    }
}

#[cfg(test)]
mod finite_number_tests {
    use super::*;

    #[test]
    fn non_finite_curvature_is_rejected() {
        for bad in ["nan", "inf", "-inf", "1e999"] {
            let text = std::format!("MODE SEQ\nSURF\nCURV {bad}\n");
            assert!(parse(&text).is_err(), "accepted CURV {bad}");
        }
    }
}

#[cfg(test)]
mod glass_hint_tests {
    use super::*;

    #[test]
    fn unknown_glass_error_suggests_a_near_miss() {
        let text = "MODE SEQ
SURF
CURV 0
DISZ 5
SURF
GLAS N-BK8 0 0 0 0
DISZ 5
SURF
";
        let doc = parse(text).unwrap();
        let err = to_sequential(
            &doc,
            &crate::glass::reference_catalog(),
            Wavelength::from_nm(550.0),
            Kelvin::from_celsius(20.0),
        )
        .unwrap_err();
        assert!(
            err.message.contains("did you mean 'N-BK7'"),
            "{}",
            err.message
        );
    }
}
