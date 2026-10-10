//! Polarization ray tracing: tracking the Jones state through a sequence of
//! Fresnel interfaces with rotating s/p frames.
//!
//! Each element applies, in its local (s, p) frame, the diagonal Fresnel
//! amplitude matrix `diag(rs, rp)`; the frames are connected by the rotation
//! angle between the plane of incidence of consecutive surfaces. This is the
//! standard reduced (2×2) polarization ray trace for isotropic coatings.

use crate::jones::{JonesMatrix, JonesVector};
use tpt_optics_core::fresnel::FresnelCoefficients;

/// One optical surface seen by a polarized ray.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct PolarizedElement {
    /// Fresnel amplitude coefficients at this incidence.
    pub fresnel: FresnelCoefficients,
    /// Angle (rad) from this surface's s-axis to the next element's s-axis
    /// (0 for a coplanar system).
    pub frame_rotation: f64,
}

impl PolarizedElement {
    /// Local (s,p) Jones operator: `diag(rs, rp)`.
    #[must_use]
    pub fn local_operator(&self) -> JonesMatrix {
        JonesMatrix::new(
            self.fresnel.rs,
            tpt_optics_core::Complex::zero(),
            tpt_optics_core::Complex::zero(),
            self.fresnel.rp,
        )
    }

    /// Local transmission operator: `diag(ts, tp)`.
    #[must_use]
    pub fn local_transmission(&self) -> JonesMatrix {
        JonesMatrix::new(
            self.fresnel.ts,
            tpt_optics_core::Complex::zero(),
            tpt_optics_core::Complex::zero(),
            self.fresnel.tp,
        )
    }
}

/// A sequence of polarizing interfaces.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct PolarizedSystem {
    /// Entrance azimuth of the input Jones reference frame relative to the
    /// first surface's s-axis (rad).
    pub entry_azimuth: f64,
    /// Interfaces in order.
    pub elements: [PolarizedElement; 2],
}

impl PolarizedSystem {
    /// Two-element (double-interface) system — window / plate traversal.
    #[must_use]
    pub fn window(fresnel_front: FresnelCoefficients, fresnel_back: FresnelCoefficients) -> Self {
        PolarizedSystem {
            entry_azimuth: 0.0,
            elements: [
                PolarizedElement {
                    fresnel: fresnel_front,
                    frame_rotation: 0.0,
                },
                PolarizedElement {
                    fresnel: fresnel_back,
                    frame_rotation: 0.0,
                },
            ],
        }
    }

    /// Reflected Jones vector through the system (reflection at each
    /// element, accounting for frame rotations).
    #[must_use]
    pub fn trace_reflection(&self, input: &JonesVector) -> JonesVector {
        let mut state = rotate_to_frame(input, self.entry_azimuth);
        for el in &self.elements {
            state = state.apply(&el.local_operator());
            if el.frame_rotation != 0.0 {
                state = rotate_to_frame(&state, el.frame_rotation);
            }
        }
        state
    }

    /// Transmitted Jones vector through the system.
    #[must_use]
    pub fn trace_transmission(&self, input: &JonesVector) -> JonesVector {
        let mut state = rotate_to_frame(input, self.entry_azimuth);
        for el in &self.elements {
            state = state.apply(&el.local_transmission());
            if el.frame_rotation != 0.0 {
                state = rotate_to_frame(&state, el.frame_rotation);
            }
        }
        state
    }
}

/// Project a Jones vector onto a frame rotated by `theta` (the components
/// are re-expressed in the rotated basis).
fn rotate_to_frame(j: &JonesVector, theta: f64) -> JonesVector {
    let c = libm::cos(theta);
    let s = libm::sin(theta);
    JonesVector::new(
        j.ex * tpt_optics_core::Complex::real(c) + j.ey * tpt_optics_core::Complex::real(s),
        -j.ex * tpt_optics_core::Complex::real(s) + j.ey * tpt_optics_core::Complex::real(c),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_optics_core::fresnel::fresnel_amplitudes;
    use tpt_optics_core::index::ComplexRefractiveIndex;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn brewster_window_transmits_p_fully() {
        // Air→glass at Brewster: rp = 0 on reflection; the transmitted p is
        // strong and the window's two surfaces leave p intact.
        let theta_b = libm::atan(1.5);
        let theta_t = libm::asin(libm::sin(theta_b) / 1.5); // inside the glass
        let fr = fresnel_amplitudes(
            ComplexRefractiveIndex::real(1.0),
            ComplexRefractiveIndex::real(1.5),
            theta_b,
        );
        let ft = fresnel_amplitudes(
            ComplexRefractiveIndex::real(1.5),
            ComplexRefractiveIndex::real(1.0),
            theta_t,
        );
        let sys = PolarizedSystem::window(fr, ft);
        // p-polarized input: (0, 1) in the s/p frame.
        let p_in = JonesVector::new(
            tpt_optics_core::Complex::zero(),
            tpt_optics_core::Complex::one(),
        );
        let t = sys.trace_transmission(&p_in);
        let t_p = t.ey.norm_sqr();
        // T = tp(front)·tp(back); at Brewster this is high (~0.92).
        assert!(t_p > 0.85, "T_p = {t_p}");
        // s-polarized: noticeable reflection loss.
        let s_in = JonesVector::linear_x();
        let ts = sys.trace_transmission(&s_in);
        assert!(
            ts.ex.norm_sqr() < t_p,
            "s should transmit less than p at Brewster"
        );
        // Reflected p vanishes.
        let r = sys.trace_reflection(&p_in);
        assert!(r.ey.abs() < 1e-12, "Brewster reflection rp must vanish");
    }

    #[test]
    fn normal_incidence_s_p_degenerate() {
        let fr = fresnel_amplitudes(
            ComplexRefractiveIndex::real(1.0),
            ComplexRefractiveIndex::real(1.5),
            0.0,
        );
        let ft = fresnel_amplitudes(
            ComplexRefractiveIndex::real(1.5),
            ComplexRefractiveIndex::real(1.0),
            0.0,
        );
        let sys = PolarizedSystem::window(fr, ft);
        for az in [0.0, 0.7, 1.3] {
            let input = JonesVector::linear_at(az);
            let t = sys.trace_transmission(&input);
            // Intensity is polarization-independent at normal incidence:
            // 0.96 per surface, 0.9216 through both.
            assert!(
                close(t.intensity(), 0.9216_f64, 1e-9),
                "az {az}: {}",
                t.intensity()
            );
        }
    }

    #[test]
    fn frame_rotation_mixes_states() {
        // Two crossed-plane elements: rotation must rotate the output
        // azimuth accordingly.
        let fr = fresnel_amplitudes(
            ComplexRefractiveIndex::real(1.0),
            ComplexRefractiveIndex::real(1.5),
            0.5,
        );
        let mut sys = PolarizedSystem::window(fr, fr);
        sys.elements[0].frame_rotation = 0.9;
        let input = JonesVector::linear_x();
        let t_rotated = sys.trace_transmission(&input);
        // With a frame rotation, x-light develops a y-component (beyond
        // what the unrotated trace would give).
        let mut plain = PolarizedSystem::window(fr, fr);
        plain.elements[0].frame_rotation = 0.0;
        let t_plain = plain.trace_transmission(&input);
        assert!(
            t_rotated.ey.abs() > t_plain.ey.abs() + 1e-3,
            "rotation should mix components"
        );
    }
}
