//! Wave-optics kernels: FFT and split-step BPM.
#![allow(missing_docs)]

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use tpt_optics_core::Complex;
use tpt_optics_wave::bpm::{bpm_2d, Bpm2dConfig};
use tpt_optics_wave::fft::fft;

fn fft_4096(c: &mut Criterion) {
    let template: Vec<Complex<f64>> = (0..4096)
        .map(|i| Complex::new((i as f64 * 0.01).sin(), 0.0))
        .collect();
    c.bench_function("fft_4096", |b| {
        b.iter(|| {
            let mut data = template.clone();
            fft(black_box(&mut data));
            data
        });
    });
}

fn bpm_256x200(c: &mut Criterion) {
    let cfg = Bpm2dConfig {
        nx: 256,
        nz: 200,
        dx: 0.1e-6,
        dz: 0.5e-6,
        wavelength_m: 1.0e-6,
        n0: 1.5,
    };
    let field: Vec<Complex<f64>> = (0..cfg.nx)
        .map(|i| {
            let x = (i as f64 - 128.0) * cfg.dx;
            Complex::new((-(x * x) / 1e-12).exp(), 0.0)
        })
        .collect();
    let profile = vec![1.5; cfg.nx * cfg.nz];
    c.bench_function("bpm_2d_256x200", |b| {
        b.iter(|| bpm_2d(black_box(&field), black_box(&profile), black_box(&cfg)));
    });
}

criterion_group!(benches, fft_4096, bpm_256x200);
criterion_main!(benches);
