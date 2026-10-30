//! Conversion between Zemax `.zmx` and CODE V `.seq` text.
//!
//! Both formats describe the same sequence of surfaces, so the conversion
//! maps each surface's radius, thickness and glass. Asphere terms and
//! semi-diameters are not part of `.seq`, so a `.seq` → `.zmx` conversion
//! sets them to zero, and a `.zmx` → `.seq` conversion drops them.

use tpt_optics::io::{seq, zmx};

/// Convert `.zmx` text to `.seq` text.
///
/// # Errors
/// Returns the parse error text when the input is not a valid `.zmx` file.
pub fn zmx_to_seq(text: &str) -> Result<String, String> {
    let doc = zmx::parse(text).map_err(|e| e.to_string())?;
    let surfaces: Vec<(f64, f64, String, bool)> = doc
        .surfaces
        .iter()
        .map(|s| {
            let radius = if s.curvature == 0.0 {
                0.0
            } else {
                1.0 / s.curvature
            };
            (radius, s.thickness, s.glass.clone(), s.stop)
        })
        .collect();
    Ok(seq::write(&surfaces))
}

/// Convert `.seq` text to `.zmx` text.
///
/// # Errors
/// Returns the parse error text when the input is not a valid `.seq` file.
pub fn seq_to_zmx(text: &str) -> Result<String, String> {
    let surfaces = seq::parse(text).map_err(|e| e.to_string())?;
    let doc = zmx::ZmxDocument {
        name: String::new(),
        surfaces: surfaces
            .into_iter()
            .map(|(radius, thickness, glass, stop)| zmx::ZmxSurface {
                curvature: if radius == 0.0 { 0.0 } else { 1.0 / radius },
                thickness,
                glass,
                stop,
                ..zmx::ZmxSurface::default()
            })
            .collect(),
    };
    Ok(zmx::write(&doc))
}

#[cfg(test)]
mod tests {
    use super::*;

    const ZMX: &str = "\
MODE SEQ
SURF
CURV 0
DISZ INFINITY
SURF
STOP
CURV 0.01
GLAS N-BK7 0 0 0 0
DISZ 5.0
SURF
CURV -0.01
GLAS AIR
DISZ 95.0
SURF
CURV 0
DISZ 0
";

    #[test]
    fn zmx_to_seq_keeps_radii_thicknesses_and_glasses() {
        let seq_text = zmx_to_seq(ZMX).unwrap();
        let back = seq::parse(&seq_text).unwrap();
        assert_eq!(back.len(), 4);
        assert!((back[1].0 - 0.1).abs() < 1e-12, "R1 = {}", back[1].0);
        assert!((back[1].1 - 5e-3).abs() < 1e-12);
        assert_eq!(back[1].2, "N-BK7");
        assert!(back[1].3, "stop flag survives");
    }

    #[test]
    fn seq_to_zmx_round_trips_the_prescription() {
        let seq_text = zmx_to_seq(ZMX).unwrap();
        let zmx_text = seq_to_zmx(&seq_text).unwrap();
        let doc = zmx::parse(&zmx_text).unwrap();
        assert_eq!(doc.surfaces.len(), 4);
        assert!((doc.surfaces[1].curvature - 10.0).abs() < 1e-9);
        assert_eq!(doc.surfaces[1].glass, "N-BK7");
    }

    #[test]
    fn invalid_input_is_an_error_text() {
        assert!(zmx_to_seq("MODE NSC\n").is_err());
        assert!(seq_to_zmx("FOO 1\n").is_err());
    }
}
