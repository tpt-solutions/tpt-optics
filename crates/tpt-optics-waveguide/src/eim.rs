//! Effective index method (EIM) for channel (rib/strip) waveguides.
//!
//! The 2-D index profile is decomposed into two cascaded 1-D slab problems:
//! first solve the vertical slab to get an effective index, then use it as
//! the core index of the lateral slab. Accuracy is a few percent for
//! high-contrast wires and degrades near cutoff.

use alloc::vec::Vec;

use tpt_optics_core::Wavelength;

use crate::slab::{slab_modes, SlabLayer, SlabPolarization};

/// Channel waveguide: core `n_core`, background `n_clad`, rectangular core
/// of `width` (x) × `height` (y).
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct ChannelGuide {
    /// Core index.
    pub n_core: f64,
    /// Cladding/background index.
    pub n_clad: f64,
    /// Core width (x, m).
    pub width: f64,
    /// Core height (y, m).
    pub height: f64,
}

/// Fundamental quasi-TE mode effective index of the channel guide via EIM.
///
/// The quasi-TE mode's dominant field is across the lateral interfaces, so
/// the vertical slab is TE and the lateral slab is TM (see
/// [`channel_wguide_neff_pol`]).
///
/// Returns `None` when one of the two slab problems has no guided mode
/// (below cutoff).
pub fn channel_wguide_neff(guide: &ChannelGuide, wl: Wavelength) -> Option<f64> {
    channel_wguide_neff_pol(guide, wl, SlabPolarization::Te)
}

/// Fundamental effective index for a chosen vertical polarization.
///
/// `vertical = Te` gives the quasi-TE mode (lateral slab TM); `vertical = Tm`
/// gives the quasi-TM mode (lateral slab TE). The lateral slab uses the
/// polarization whose dominant field lies across the lateral interfaces.
pub fn channel_wguide_neff_pol(
    guide: &ChannelGuide,
    wl: Wavelength,
    vertical: SlabPolarization,
) -> Option<f64> {
    let lateral = match vertical {
        SlabPolarization::Te => SlabPolarization::Tm,
        SlabPolarization::Tm => SlabPolarization::Te,
    };
    // Vertical slab (infinite in x).
    let vert = SlabLayer {
        n_core: guide.n_core,
        n_clad: guide.n_clad,
        thickness: guide.height,
    };
    let mut vmodes = slab_modes(&vert, wl, vertical);
    vmodes.sort_by(|a, b| b.n_eff.total_cmp(&a.n_eff));
    let n_vert = vmodes.first().map(|m| m.n_eff)?;
    // Lateral slab with the effective core index.
    let lat = SlabLayer {
        n_core: n_vert,
        n_clad: guide.n_clad,
        thickness: guide.width,
    };
    let mut lmodes: Vec<_> = slab_modes(&lat, wl, lateral);
    lmodes.sort_by(|a, b| b.n_eff.total_cmp(&a.n_eff));
    lmodes.first().map(|m| m.n_eff)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn wl(nm: f64) -> Wavelength {
        Wavelength::from_nm(nm)
    }

    #[test]
    fn soi_wire_neff_in_bounds() {
        // 450×220 nm SOI wire at 1550 nm: n_eff ≈ 2.3–2.5 (EIM ballpark).
        let guide = ChannelGuide {
            n_core: 3.48,
            n_clad: 1.444,
            width: 450e-9,
            height: 220e-9,
        };
        let neff = channel_wguide_neff(&guide, wl(1550.0)).expect("guided");
        assert!((neff - 2.4).abs() < 0.3, "n_eff = {neff}");
        assert!(neff > guide.n_clad && neff < guide.n_core);
    }

    #[test]
    fn wider_is_better_guided() {
        let base = ChannelGuide {
            n_core: 3.48,
            n_clad: 1.444,
            width: 400e-9,
            height: 220e-9,
        };
        let wide = ChannelGuide {
            width: 800e-9,
            ..base
        };
        let n1 = channel_wguide_neff(&base, wl(1550.0)).unwrap();
        let n2 = channel_wguide_neff(&wide, wl(1550.0)).unwrap();
        assert!(
            n2 > n1,
            "wider guide should have higher n_eff: {n1} vs {n2}"
        );
    }

    #[test]
    fn narrow_guide_near_cutoff() {
        // 120 nm wire at 4.5 µm: vertical slab may be below cutoff → None or
        // n_eff barely above cladding.
        let guide = ChannelGuide {
            n_core: 3.48,
            n_clad: 1.444,
            width: 120e-9,
            height: 60e-9,
        };
        if let Some(neff) = channel_wguide_neff(&guide, wl(4500.0)) {
            assert!(neff < 1.7, "near-cutoff n_eff should be low: {neff}");
        }
    }
}

#[cfg(test)]
mod quasi_polarization_tests {
    use super::*;

    #[test]
    fn quasi_te_index_exceeds_quasi_tm_for_soi_wire() {
        // A 450 × 220 nm Si wire in SiO₂ at 1550 nm: the quasi-TE mode (field
        // mostly across the wide, lateral interfaces) has the higher index.
        let guide = ChannelGuide {
            n_core: 3.476,
            n_clad: 1.444,
            width: 450e-9,
            height: 220e-9,
        };
        let wl = Wavelength::from_nm(1550.0);
        let te = channel_wguide_neff(&guide, wl).expect("quasi-TE guided");
        let tm =
            channel_wguide_neff_pol(&guide, wl, SlabPolarization::Tm).expect("quasi-TM guided");
        assert!(te > tm, "quasi-TE {te} vs quasi-TM {tm}");
        assert!(te < guide.n_core && tm > guide.n_clad);
    }
}
