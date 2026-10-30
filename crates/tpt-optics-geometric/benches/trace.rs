//! Sequential ray tracing through a two-surface glass slab.
#![allow(missing_docs)]

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use tpt_optics_core::{Kelvin, Vec3, Wavelength};
use tpt_optics_geometric::prelude::Ray;
use tpt_optics_geometric::sequential::{Medium, SeqSurface, SequentialSystem, TraceContext};

fn trace_slab(c: &mut Criterion) {
    let mut sys = SequentialSystem::new();
    sys.push(SeqSurface::flat(12e-3), 5e-3, Medium::Const(1.5));
    sys.push(SeqSurface::flat(12e-3), 0.0, Medium::Vacuum);
    let ctx = TraceContext {
        wavelength: Wavelength::from_nm(587.6),
        temperature: Kelvin::from_celsius(20.0),
    };
    let ray = Ray::new(Vec3::new(1e-3, 0.0, -5e-3), Vec3::new(0.0, 0.0, 1.0));
    c.bench_function("trace_slab_two_surfaces", |b| {
        b.iter(|| sys.trace(black_box(&ray), black_box(ctx)));
    });
}

criterion_group!(benches, trace_slab);
criterion_main!(benches);
