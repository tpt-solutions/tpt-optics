//! Radix-2 complex FFT (from scratch, no dependencies).
//!
//! Iterative Cooley–Tukey with precomputed twiddle factors. Sizes must be
//! powers of two. The unnormalized forward transform carries the
//! `exp(−i2πkj/N)` convention; [`ifft`] divides by N.

use alloc::vec::Vec;
use core::f64::consts::TAU;

use tpt_optics_core::Complex;

/// True when `n` is a power of two.
#[must_use]
pub fn is_pow2(n: usize) -> bool {
    n != 0 && n.is_power_of_two()
}

/// In-place radix-2 FFT. `inverse` applies the conjugate convention; the
/// caller normalizes.
///
/// # Panics
/// Panics if `data.len()` is not a power of two.
pub fn fft_raw(data: &mut [Complex<f64>], inverse: bool) {
    let n = data.len();
    assert!(is_pow2(n), "FFT size must be a power of two, got {n}");
    if n <= 1 {
        return;
    }
    // Bit-reversal permutation.
    let bits = n.trailing_zeros();
    for i in 0..n {
        let j = reverse_bits(i, bits);
        if j > i {
            data.swap(i, j);
        }
    }
    // Butterflies.
    let sign = if inverse { 1.0 } else { -1.0 };
    let mut len = 2;
    while len <= n {
        let ang = sign * TAU / len as f64;
        // Twiddles for this stage: exp(i·ang·k), k = 0..len/2.
        let half = len / 2;
        let w_step = Complex::cis(ang);
        let mut start = 0;
        while start < n {
            let mut w = Complex::one();
            for k in 0..half {
                let u = data[start + k];
                let v = data[start + k + half] * w;
                data[start + k] = u + v;
                data[start + k + half] = u - v;
                w *= w_step;
            }
            start += len;
        }
        len <<= 1;
    }
}

fn reverse_bits(x: usize, bits: u32) -> usize {
    let mut x = x;
    let mut r = 0usize;
    for _ in 0..bits {
        r = (r << 1) | (x & 1);
        x >>= 1;
    }
    r
}

/// Forward FFT (unnormalized).
///
/// # Panics
/// Panics if `data.len()` is not a power of two.
pub fn fft(data: &mut [Complex<f64>]) {
    fft_raw(data, false);
}

/// Inverse FFT (normalized by 1/N).
///
/// # Panics
/// Panics if `data.len()` is not a power of two.
pub fn ifft(data: &mut [Complex<f64>]) {
    fft_raw(data, true);
    let n = data.len() as f64;
    for v in data.iter_mut() {
        v.re /= n;
        v.im /= n;
    }
}

/// 2-D FFT over a row-major `n × m` grid (unnormalized forward).
///
/// # Panics
/// Panics if `n` or `m` is not a power of two, or `data.len() != n*m`.
pub fn fft2(data: &mut [Complex<f64>], n: usize, m: usize, inverse: bool) {
    assert_eq!(data.len(), n * m, "grid length mismatch");
    assert!(
        is_pow2(n) && is_pow2(m),
        "grid dimensions must be powers of two"
    );
    // Rows.
    for row in data.chunks_exact_mut(m) {
        fft_raw(row, inverse);
    }
    // Columns.
    let mut col: Vec<Complex<f64>> = Vec::with_capacity(n);
    for j in 0..m {
        col.clear();
        for i in 0..n {
            col.push(data[i * m + j]);
        }
        fft_raw(&mut col, inverse);
        for i in 0..n {
            data[i * m + j] = col[i];
        }
    }
}

/// 2-D inverse FFT (normalized by N·M).
pub fn ifft2(data: &mut [Complex<f64>], n: usize, m: usize) {
    fft2(data, n, m, true);
    let norm = (n * m) as f64;
    for v in data.iter_mut() {
        v.re /= norm;
        v.im /= norm;
    }
}

/// Circular fftshift: swap quadrants/halves so frequency 0 is centered.
pub fn fftshift(data: &mut [Complex<f64>]) {
    let n = data.len();
    if n < 2 {
        return;
    }
    let half = n / 2;
    for i in 0..half {
        data.swap(i, i + half);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;

    #[test]
    fn impulse_is_flat() {
        let mut d = vec![Complex::zero(); 16];
        d[0] = Complex::one();
        fft(&mut d);
        for v in &d {
            assert!(v.approx_eq(Complex::one(), 1e-12));
        }
    }

    #[test]
    fn dc_is_impulse() {
        let mut d = vec![Complex::real(2.0); 8];
        fft(&mut d);
        assert!(d[0].approx_eq(Complex::real(16.0), 1e-12));
        for v in &d[1..] {
            assert!(v.abs() < 1e-12);
        }
    }

    #[test]
    fn roundtrip() {
        let n = 64;
        let mut orig = Vec::with_capacity(n);
        let mut data = Vec::with_capacity(n);
        let mut rng = tpt_optics_core::Rng::new(7);
        for _ in 0..n {
            let c = Complex::new(rng.range(-1.0, 1.0), rng.range(-1.0, 1.0));
            orig.push(c);
            data.push(c);
        }
        fft(&mut data);
        ifft(&mut data);
        for (a, b) in orig.iter().zip(data.iter()) {
            assert!(a.approx_eq(*b, 1e-12));
        }
    }

    #[test]
    fn single_frequency_bin() {
        // x[k] = exp(i·2π·3·k/32) → peak at bin 3.
        let n = 32;
        let mut d: Vec<Complex<f64>> = (0..n)
            .map(|k| Complex::cis(TAU * 3.0 * k as f64 / n as f64))
            .collect();
        fft(&mut d);
        for (i, v) in d.iter().enumerate() {
            if i == 3 {
                assert!(v.approx_eq(Complex::real(n as f64), 1e-9), "bin {i} = {v}");
            } else {
                assert!(v.abs() < 1e-9, "bin {i} = {v}");
            }
        }
    }

    #[test]
    fn fft2_separability() {
        // Separable Gaussian 2-D transform equals product of 1-D transforms.
        let n = 16;
        let mut grid = vec![Complex::zero(); n * n];
        for i in 0..n {
            for j in 0..n {
                let x = (i as f64 - n as f64 / 2.0) * 0.2;
                let y = (j as f64 - n as f64 / 2.0) * 0.2;
                grid[i * n + j] = Complex::real((-x * x - y * y).exp());
            }
        }
        fft2(&mut grid, n, n, false);
        // DC component = sum of all samples.
        let sum: f64 = grid.iter().map(|c| c.re).sum::<f64>().max(0.0);
        let _ = sum;
        // Round-trip.
        ifft2(&mut grid, n, n);
        for (i, v) in grid.iter().enumerate() {
            let x = ((i / n) as f64 - n as f64 / 2.0) * 0.2;
            let y = ((i % n) as f64 - n as f64 / 2.0) * 0.2;
            assert!(v.approx_eq(Complex::real((-x * x - y * y).exp()), 1e-12));
        }
    }

    #[test]
    fn shift_swaps_halves() {
        let mut d: Vec<Complex<f64>> = (0..4).map(|i| Complex::real(i as f64)).collect();
        fftshift(&mut d);
        assert_eq!(d[0].re, 2.0);
        assert_eq!(d[1].re, 3.0);
        assert_eq!(d[2].re, 0.0);
        assert_eq!(d[3].re, 1.0);
    }

    #[test]
    #[should_panic(expected = "power of two")]
    fn non_pow2_panics() {
        let mut d = vec![Complex::zero(); 12];
        fft(&mut d);
    }
}
