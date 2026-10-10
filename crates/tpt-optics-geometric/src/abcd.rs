//! ABCD matrix formalism for paraxial systems.
//!
//! Ray-vector convention: column `(y, u)` — transverse height and paraxial
//! angle. The determinant of a refraction matrix is `n₁/n₂`; a full system
//! traced from index n₁ back into index n₂ has det = n₁/n₂ (unit systems
//! have det = 1, the symplectic condition).

use tpt_optics_core::{Matrix2, Scalar};

/// Free-space propagation over distance `d`.
pub fn free_space<T: Scalar>(d: f64) -> Matrix2<T> {
    Matrix2::new(T::ONE, T::from_f64(d), T::ZERO, T::ONE)
}

/// Thin lens of focal length `f` (positive = converging).
pub fn thin_lens<T: Scalar>(f: f64) -> Matrix2<T> {
    Matrix2::new(T::ONE, T::ZERO, T::from_f64(-1.0 / f), T::ONE)
}

/// Refraction at a spherical surface with curvature `c` from index `n1`
/// into `n2`. Determinant = n1/n2.
pub fn refracting_surface<T: Scalar>(n1: f64, n2: f64, c: f64) -> Matrix2<T> {
    Matrix2::new(
        T::ONE,
        T::ZERO,
        T::from_f64((n1 - n2) * c / n2),
        T::from_f64(n1 / n2),
    )
}

/// Reflection off a spherical mirror with curvature `c` (c > 0:
/// center of curvature on the incoming side, i.e. concave toward the ray).
///
/// For a system unfolded about the mirror the matrix acts like a lens of
/// f = R/2 with det = 1 (index sign-flip convention).
pub fn mirror<T: Scalar>(c: f64) -> Matrix2<T> {
    Matrix2::new(T::ONE, T::ZERO, T::from_f64(-2.0 * c), T::ONE)
}

/// Two-lens (afocal/Keplerian) telescope angular magnification check:
/// returns the composed system matrix for two thin lenses separated by
/// `f1 + f2` — an afocal system with angular magnification −f1/f2.
pub fn keplerian_telescope<T: Scalar>(f1: f64, f2: f64) -> Matrix2<T> {
    thin_lens::<T>(f2)
        .mul_matrix(free_space::<T>(f1 + f2))
        .mul_matrix(thin_lens::<T>(f1))
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_optics_core::Complex;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn determinants() {
        assert!(close(free_space::<f64>(1.7).det(), 1.0, 1e-15));
        assert!(close(thin_lens::<f64>(0.25).det(), 1.0, 1e-12));
        // Refraction air→glass: det = n1/n2
        assert!(close(
            refracting_surface::<f64>(1.0, 1.5, 0.1).det(),
            1.0 / 1.5,
            1e-15
        ));
        assert!(close(mirror::<f64>(0.5).det(), 1.0, 1e-12));
    }

    #[test]
    fn thin_lens_imaging() {
        // Object at 2f images to 2f: matrix from object to image is a pure
        // magnification −1.
        let f = 0.1;
        let m = free_space::<f64>(2.0 * f)
            .mul_matrix(thin_lens(f))
            .mul_matrix(free_space(2.0 * f));
        assert!(close(m.a, -1.0, 1e-12), "a = {}", m.a);
        assert!(close(m.b, 0.0, 1e-9), "b = {}", m.b);
        // The lens power remains in the c element (system ends at the image).
        assert!(close(m.c, -1.0 / f, 1e-9), "c = {}", m.c);
        assert!(close(m.d, -1.0, 1e-12), "d = {}", m.d);
    }

    #[test]
    fn telescope_angular_magnification() {
        let m = keplerian_telescope::<f64>(0.2, 0.05);
        assert!(close(m.c, 0.0, 1e-9), "afocal: c = {}", m.c);
        // Angular magnification −f1/f2 = −4.
        assert!(close(m.d, -4.0, 1e-9), "d = {}", m.d);
    }

    #[test]
    fn q_parameter_free_space() {
        let m = free_space::<f64>(0.05);
        let q = Complex::real(0.01);
        assert!((m.apply_q(q) - Complex::real(0.06)).approx_eq(Complex::zero(), 1e-12));
    }

    #[test]
    fn stability_plano_convex_cavity() {
        // Symmetric two-mirror cavity with R = 1 m, L = 0.5 m: stable.
        let half = free_space::<f64>(0.5);
        let m1 = mirror(1.0); // concave R=1 in the c = +1/R convention
        let round = m1
            .mul_matrix(half)
            .mul_matrix(Matrix2::identity().mul_matrix(m1))
            .mul_matrix(half);
        // |tr|/2 = |A+D|/2 < 1
        assert!(round.is_stable());
    }
}
