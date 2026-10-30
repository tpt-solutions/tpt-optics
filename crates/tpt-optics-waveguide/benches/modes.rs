//! Mode solvers: scalar and full-vector finite-element fundamental modes.
#![allow(missing_docs)]

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use tpt_optics_core::Wavelength;
use tpt_optics_waveguide::fem::{fem_fundamental_mode, RectDomain};
use tpt_optics_waveguide::vector_fem::vector_mode_neff;

fn channel(x: f64, y: f64) -> f64 {
    if (x - 1.5e-6).abs() < 0.25e-6 && (y - 1.5e-6).abs() < 0.25e-6 {
        3.48
    } else {
        1.444
    }
}

fn modes(c: &mut Criterion) {
    let wl = Wavelength::from_nm(1550.0);
    let domain = RectDomain {
        width: 3e-6,
        height: 3e-6,
        nx: 16,
        ny: 16,
    };
    let mut group = c.benchmark_group("modes");
    group.sample_size(10);
    group.bench_function("scalar_fem_16x16", |b| {
        b.iter(|| fem_fundamental_mode(black_box(&domain), &channel, wl, 30));
    });
    group.bench_function("vector_fem_16x16", |b| {
        b.iter(|| vector_mode_neff(black_box(&domain), &channel, wl, 2.4, 30));
    });
    group.finish();
}

criterion_group!(benches, modes);
criterion_main!(benches);
