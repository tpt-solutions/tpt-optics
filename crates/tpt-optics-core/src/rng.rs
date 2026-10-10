//! Deterministic pseudo-random number generation for Monte-Carlo optics.
//!
//! A small xoshiro256** generator (seeded through SplitMix64) plus the
//! sampling helpers the non-sequential tracer and tolerancing need. Fully
//! reproducible: same seed → same trajectory, which keeps Monte-Carlo tests
//! deterministic.

/// xoshiro256** PRNG.
#[derive(Clone, Debug, PartialEq)]
pub struct Rng {
    s: [u64; 4],
}

impl Rng {
    /// Seed from any `u64` via SplitMix64.
    #[must_use]
    pub fn new(seed: u64) -> Self {
        let mut sm = seed;
        let mut next = || {
            sm = sm.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = sm;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^ (z >> 31)
        };
        Rng {
            s: [next(), next(), next(), next()],
        }
    }

    /// Next raw 64-bit value.
    pub fn next_u64(&mut self) -> u64 {
        let result = self.s[1].wrapping_mul(5).rotate_left(7).wrapping_mul(9);
        let t = self.s[1] << 17;
        self.s[2] ^= self.s[0];
        self.s[3] ^= self.s[1];
        self.s[1] ^= self.s[2];
        self.s[0] ^= self.s[3];
        self.s[2] ^= t;
        self.s[3] = self.s[3].rotate_left(45);
        result
    }

    /// Uniform on [0, 1).
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }

    /// Uniform on [a, b).
    pub fn range(&mut self, a: f64, b: f64) -> f64 {
        a + (b - a) * self.next_f64()
    }

    /// Standard normal via Box–Muller.
    pub fn normal(&mut self, mean: f64, std_dev: f64) -> f64 {
        let u1 = self.next_f64().max(1e-15);
        let u2 = self.next_f64();
        let r = libm::sqrt(-2.0 * libm::log(u1));
        let theta = core::f64::consts::TAU * u2;
        mean + std_dev * r * libm::sin(theta)
    }

    /// Uniform direction on the unit sphere (for scattering models).
    pub fn unit_direction(&mut self) -> (f64, f64, f64) {
        let mu = self.range(-1.0, 1.0);
        let phi = self.range(0.0, core::f64::consts::TAU);
        let s = libm::sqrt(1.0 - mu * mu);
        (s * libm::cos(phi), s * libm::sin(phi), mu)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reproducible() {
        let mut a = Rng::new(42);
        let mut b = Rng::new(42);
        for _ in 0..100 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
        let mut c = Rng::new(43);
        assert_ne!(a.next_u64(), c.next_u64());
    }

    #[test]
    fn uniform_statistics() {
        let mut rng = Rng::new(7);
        let n = 200_000;
        let mut sum = 0.0;
        let mut sum2 = 0.0;
        let mut in_bounds = true;
        for _ in 0..n {
            let v = rng.next_f64();
            in_bounds &= (0.0..1.0).contains(&v);
            sum += v;
            sum2 += v * v;
        }
        assert!(in_bounds);
        let mean = sum / n as f64;
        let var = sum2 / n as f64 - mean * mean;
        assert!((mean - 0.5).abs() < 0.005, "mean {mean}");
        assert!((var - 1.0 / 12.0).abs() < 0.005, "var {var}");
    }

    #[test]
    fn normal_statistics() {
        let mut rng = Rng::new(11);
        let n = 200_000;
        let mut sum = 0.0;
        let mut sum2 = 0.0;
        for _ in 0..n {
            let v = rng.normal(0.0, 1.0);
            sum += v;
            sum2 += v * v;
        }
        let mean = sum / n as f64;
        let var = sum2 / n as f64 - mean * mean;
        assert!(mean.abs() < 0.01, "mean {mean}");
        assert!((var - 1.0).abs() < 0.02, "var {var}");
    }

    #[test]
    fn unit_directions_on_sphere() {
        let mut rng = Rng::new(3);
        let mut mean_z = 0.0;
        for _ in 0..10_000 {
            let (x, y, z) = rng.unit_direction();
            assert!((x * x + y * y + z * z - 1.0).abs() < 1e-12);
            mean_z += z;
        }
        assert!((mean_z / 10_000.0).abs() < 0.03);
    }

    #[test]
    fn range_within_bounds() {
        let mut rng = Rng::new(5);
        for _ in 0..1000 {
            let v = rng.range(2.0, 3.0);
            assert!((2.0..3.0).contains(&v));
        }
    }
}
