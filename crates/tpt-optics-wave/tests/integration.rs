//! Phase 2 exit criterion: physical wave propagation must agree with the
//! Gaussian q-parameter ABCD law.
//!
//! A Gaussian is propagated by the angular-spectrum method through a
//! thin-lens phase and free space; the resulting beam radii are compared
//! with the q-parameter prediction `q → (Aq+B)/(Cq+D)`.
//!
//! Geometry uses CO₂-laser scales (λ = 10.6 µm, mm beams) so the focused
//! spots are well resolved on the sampling grid.

use tpt_optics_core::{Complex, Wavelength};
use tpt_optics_geometric::abcd::{free_space, thin_lens};
use tpt_optics_wave::diffraction::angular_spectrum;
use tpt_optics_wave::gaussian::GaussianBeam;

const N: usize = 512;
const DX: f64 = 20e-6; // 10.24 mm window
const WL_M: f64 = 10.6e-6;

/// 2-D Gaussian field (1/e amplitude radius `w`) on an n×n grid.
fn gaussian_field(n: usize, dx: f64, w: f64, curvature_r: Option<f64>) -> Vec<Complex<f64>> {
    let k = core::f64::consts::TAU / WL_M;
    (0..n * n)
        .map(|i| {
            let x = ((i / n) as f64 - n as f64 / 2.0) * dx;
            let y = ((i % n) as f64 - n as f64 / 2.0) * dx;
            let r2 = x * x + y * y;
            let phase = curvature_r.map_or(0.0, |r| -k * r2 / (2.0 * r));
            Complex::cis(phase) * Complex::real((-r2 / (w * w)).exp())
        })
        .collect()
}

/// Apply a thin-lens phase `exp(−i·k·r²/(2f))` to a field.
fn apply_lens_phase(field: &mut [Complex<f64>], n: usize, dx: f64, f: f64) {
    let k = core::f64::consts::TAU / WL_M;
    for (i, v) in field.iter_mut().enumerate() {
        let x = ((i / n) as f64 - n as f64 / 2.0) * dx;
        let y = ((i % n) as f64 - n as f64 / 2.0) * dx;
        let r2 = x * x + y * y;
        *v *= Complex::cis(-k * r2 / (2.0 * f));
    }
}

/// 1/e (amplitude) radius along +x from the grid center.
fn measure_radius(field: &[Complex<f64>], n: usize, dx: f64) -> f64 {
    let c = n / 2;
    let peak = field[c * n + c].abs();
    for (i, v) in field.iter().enumerate().skip(c * n + c + 1).take(n - c - 1) {
        if v.abs() < peak / core::f64::consts::E {
            return ((i % n) - c) as f64 * dx;
        }
    }
    f64::NAN
}

fn wl() -> Wavelength {
    Wavelength::from_m(WL_M)
}

#[test]
fn gaussian_through_lens_matches_q_parameter_abcd_law() {
    let w0 = 2.0e-3;
    let f = 80.0e-3;
    let z_planes = [60.0e-3, 80.0e-3, 100.0e-3];

    let beam_in = GaussianBeam::at_waist(w0, wl(), 1.0);
    let mut field = gaussian_field(N, DX, w0, None);
    apply_lens_phase(&mut field, N, DX, f);

    for &z in &z_planes {
        let physical = angular_spectrum(&field, N, DX, WL_M, z, 1.0);
        let w_phys = measure_radius(&physical, N, DX);

        // Analytic: lens then free space — the ABCD law on q.
        let m = free_space::<f64>(z).mul_matrix(thin_lens(f));
        let q_out = m.apply_q(beam_in.q);
        let w_abcd = GaussianBeam {
            q: q_out,
            wavelength: wl(),
            index: 1.0,
        }
        .w();

        let rel = (w_phys - w_abcd).abs() / w_abcd;
        assert!(
            rel < 0.06,
            "z = {z}: physical w = {w_phys:.4e}, ABCD w = {w_abcd:.4e} (rel {rel:.3})"
        );
    }
}

#[test]
fn collimated_beam_focuses_where_abcd_predicts() {
    // A near-collimated beam (waist 1 m upstream, Rayleigh 1.19 m) focused
    // by f = 100 mm: the q-law predicts both the focus position (≈ f) and
    // the focal spot size.
    let w0_upstream = 2.0e-3;
    let f = 100.0e-3;

    let beam0 = GaussianBeam::at_waist(w0_upstream, wl(), 1.0);
    let beam_lens = beam0.propagated(1.0); // 1 m downstream of the waist
    let w_lens = beam_lens.w();
    let r_lens = beam_lens.r().unwrap();

    // Physical: analytic field (amplitude + curvature) at the lens, lens
    // phase, then propagate to the nominal focus.
    let mut field = gaussian_field(N, DX, w_lens, Some(r_lens));
    apply_lens_phase(&mut field, N, DX, f);
    let out = angular_spectrum(&field, N, DX, WL_M, f, 1.0);
    let w_phys = measure_radius(&out, N, DX);

    // Analytic: ABCD through lens then f.
    let m = free_space::<f64>(f).mul_matrix(thin_lens(f));
    let w_abcd = GaussianBeam {
        q: m.apply_q(beam_lens.q),
        wavelength: wl(),
        index: 1.0,
    }
    .w();

    // Sanity: the confocal-ish spot must be much smaller than the beam.
    assert!(w_abcd < 0.1 * w_lens);
    assert!(
        (w_phys - w_abcd).abs() < 0.1 * w_abcd,
        "physical {w_phys:.4e} vs ABCD {w_abcd:.4e}"
    );
}
