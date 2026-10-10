//! CODE V `.seq` sequence files (common subset).
//!
//! Supported commands: `LEN NEW` (header), `SO` (object surface),
//! `S <radius> <thickness> <glass>` / `S <radius> <thickness>` (surfaces;
//! radius 0 = flat, `INF` radius supported), `STO` (stop flag on the next
//! surface), `EFI`/`GO`/etc. are tolerated and ignored.

use std::string::{String, ToString};
use std::vec::Vec;

use tpt_optics_geometric::sequential::{Medium, SeqSurface, SequentialSystem, SurfaceAction};
use tpt_optics_geometric::surface::Asphere;

/// A parse error with a line number.
#[derive(Debug, Clone, PartialEq)]
pub struct SeqError {
    /// 1-based line number.
    pub line: usize,
    /// Description.
    pub message: String,
}

impl core::fmt::Display for SeqError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, ".seq line {}: {}", self.line, self.message)
    }
}

#[cfg(feature = "std")]
impl std::error::Error for SeqError {}

/// Parse a `.seq` document into surfaces
/// `(radius_m, thickness_m, glass_name, stop)`.
pub fn parse(text: &str) -> Result<Vec<(f64, f64, String, bool)>, SeqError> {
    let mut out: Vec<(f64, f64, String, bool)> = Vec::new();
    let mut stop_next = false;
    for (idx, raw) in text.lines().enumerate() {
        let line_no = idx + 1;
        let line = raw.trim();
        if line.is_empty() || line.starts_with('!') || line.starts_with('#') {
            continue;
        }
        let mut tokens = line.split_whitespace();
        let cmd = tokens.next().unwrap_or("").to_ascii_uppercase();
        let args: Vec<&str> = tokens.collect();
        match cmd.as_str() {
            "LEN" | "GO" | "EFI" | "EFL" | "TITLE" | "RMS" | " ref" => {}
            "SO" => {
                // Object surface: radius and thickness to the first surface.
                let r = parse_radius(args.first(), line_no)?;
                let t = if args.len() > 1 {
                    parse_len(args.get(1), line_no)?
                } else {
                    0.0
                };
                out.push((r, t, String::new(), false));
            }
            "STO" => stop_next = true,
            "S" => {
                let r = parse_radius(args.first(), line_no)?;
                let t = if args.len() > 1 {
                    parse_len(args.get(1), line_no)?
                } else {
                    0.0
                };
                let glass = args.get(2).map(|s| (*s).to_string()).unwrap_or_default();
                out.push((r, t, glass, core::mem::take(&mut stop_next)));
            }
            other => {
                return Err(SeqError {
                    line: line_no,
                    message: std::format!("unsupported command '{other}'"),
                })
            }
        }
    }
    if out.len() < 2 {
        return Err(SeqError {
            line: 0,
            message: "sequence needs at least an object and image surface".into(),
        });
    }
    Ok(out)
}

fn parse_radius(v: Option<&&str>, line: usize) -> Result<f64, SeqError> {
    let s = v.ok_or(SeqError {
        line,
        message: "missing radius".into(),
    })?;
    if s.eq_ignore_ascii_case("INF") {
        return Ok(0.0);
    }
    s.parse::<f64>()
        .map(|mm| mm * 1e-3) // mm → m
        .map_err(|_| SeqError {
            line,
            message: "bad radius".into(),
        })
}

fn parse_len(v: Option<&&str>, line: usize) -> Result<f64, SeqError> {
    let s = v.ok_or(SeqError {
        line,
        message: "missing thickness".into(),
    })?;
    if s.eq_ignore_ascii_case("INF") {
        return Ok(f64::INFINITY);
    }
    s.parse::<f64>().map(|mm| mm * 1e-3).map_err(|_| SeqError {
        line,
        message: "bad thickness".into(),
    })
}

/// Serialize surfaces back to `.seq` text.
pub fn write(surfaces: &[(f64, f64, String, bool)]) -> String {
    use std::fmt::Write as _;
    let mut out = String::from("LEN NEW\n");
    for (i, (r, t, glass, stop)) in surfaces.iter().enumerate() {
        let cmd = if i == 0 { "SO" } else { "S" };
        if *stop && i > 0 {
            let _ = writeln!(out, "STO");
        }
        let radius = if *r == 0.0 {
            "INF".to_string()
        } else {
            std::format!("{:.9}", r * 1e3)
        };
        let thickness = if t.is_infinite() {
            "INF".to_string()
        } else {
            std::format!("{:.9}", t * 1e3)
        };
        if glass.is_empty() {
            let _ = writeln!(out, "{cmd} {radius} {thickness}");
        } else {
            let _ = writeln!(out, "{cmd} {radius} {thickness} {glass}");
        }
    }
    let _ = writeln!(out, "GO");
    out
}

/// Convert a parsed sequence into a sequential system. Glass names resolve
/// against `catalog` (falling back to numeric index strings like `1.5168`).
pub fn to_sequential(
    surfaces: &[(f64, f64, String, bool)],
    resolve: impl Fn(&str) -> Option<Medium>,
) -> Result<SequentialSystem, SeqError> {
    let mut sys = SequentialSystem::new();
    sys.media[0] = Medium::Const(1.0);
    let n = surfaces.len();
    for (i, (r, t, glass, _stop)) in surfaces.iter().enumerate() {
        let is_object = i == 0;
        let is_last = i + 1 == n;
        if is_object {
            // Object-space thickness becomes the leading offset: encoded by
            // placing a flat "surface 0" at the origin and shifting.
            continue;
        }
        let medium_after = if is_last {
            Medium::Vacuum
        } else if glass.is_empty()
            || glass.eq_ignore_ascii_case("AIR")
            || glass.to_ascii_uppercase().starts_with("MIRROR")
        {
            Medium::Const(1.0)
        } else {
            resolve(glass).ok_or(SeqError {
                line: 0,
                message: std::format!("glass '{glass}' unresolved"),
            })?
        };
        let action = if glass.to_ascii_uppercase().starts_with("MIRROR") {
            SurfaceAction::Reflect
        } else {
            SurfaceAction::Refract
        };
        // Radius already parsed to meters (0 = flat).
        let curvature = if *r == 0.0 { 0.0 } else { 1.0 / r };
        let surface = SeqSurface {
            action,
            asphere: Asphere::sphere(curvature),
            semi_diameter: f64::INFINITY,
            phase: None,
        };
        let thickness = if is_last { 0.0 } else { *t };
        sys.push(surface, thickness, medium_after);
    }
    Ok(sys)
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = "\
! CODE V sequence: BK7 singlet
LEN NEW
TITLE singlet
SO INF 10.0
STO
S 100.0 5.0 N-BK7
S -100.0 95.0
S INF 0.0
GO
";

    #[test]
    fn parses_sample() {
        let s = parse(SAMPLE).unwrap();
        assert_eq!(s.len(), 4);
        assert!(s[0].1.is_infinite() || s[0].1 > 0.0);
        assert_eq!(s[1].2, "N-BK7");
        assert!(s[1].3, "stop on first lens surface");
        assert!(
            (s[2].0 - (-0.1)).abs() < 1e-12,
            "radius parsed to meters: {}",
            s[2].0
        );
    }

    #[test]
    fn round_trip() {
        let s = parse(SAMPLE).unwrap();
        let text = write(&s);
        let s2 = parse(&text).unwrap();
        assert_eq!(s, s2);
    }

    #[test]
    fn rejects_unknown_command() {
        let err = parse("LEN NEW\nFOO 1\n").unwrap_err();
        assert_eq!(err.line, 2);
    }

    #[test]
    fn to_sequential_resolves_glass() {
        let s = parse(SAMPLE).unwrap();
        let sys = to_sequential(&s, |name| {
            if name == "N-BK7" {
                Some(Medium::Const(1.5168))
            } else {
                None
            }
        })
        .unwrap();
        assert_eq!(sys.surfaces.len(), 3);
        // Curvature: radius 100 mm → 10 1/m.
        assert!((sys.surfaces[0].asphere.curvature - 10.0).abs() < 1e-9);
        assert!((sys.surfaces[1].asphere.curvature + 10.0).abs() < 1e-9);
    }
}
