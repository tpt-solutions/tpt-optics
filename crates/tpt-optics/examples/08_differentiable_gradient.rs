//! Exact gradient of a traced quantity with respect to a lens curvature,
//! from dual-number ray tracing, checked against a finite difference.
//!
//! Run with: `cargo run -p tpt-optics --example 08_differentiable_gradient --features full`

use tpt_optics::core::material::MaterialDatabase;
use tpt_optics::core::{Kelvin, Vec3, Wavelength};
use tpt_optics::geometric::differentiable::gradient_of;
use tpt_optics::geometric::prelude::Ray;
use tpt_optics::geometric::sequential::{
    Medium, ParamId, SeqSurface, SequentialSystem, TraceContext,
};
use tpt_optics::geometric::surface::Asphere;

/// Singlet whose front surface has curvature `c1` (1/m).
fn singlet(c1: f64) -> SequentialSystem {
    let bk7 = MaterialDatabase::presets().get("bk7").expect("bk7 preset");
    let mut sys = SequentialSystem::new();
    sys.push(
        SeqSurface::aspheric(Asphere::sphere(c1), 12.5e-3),
        5e-3,
        Medium::Material(bk7),
    );
    sys.push(
        SeqSurface::aspheric(Asphere::sphere(-10.0), 12.5e-3),
        95e-3,
        Medium::Vacuum,
    );
    sys
}

fn main() {
    let ctx = TraceContext {
        wavelength: Wavelength::from_nm(587.6),
        temperature: Kelvin::from_celsius(20.0),
    };
    // A ray 4 mm off axis, travelling parallel to the axis.
    let ray = Ray::new(Vec3::new(4e-3, 0.0, -20e-3), Vec3::new(0.0, 0.0, 1.0));

    let c1 = 10.0;
    let sys = singlet(c1);

    // Quantity: the x component of the exit direction (a ray-angle proxy).
    let analytic = gradient_of(&sys, &ray, ctx, &[ParamId::Curvature(0)], |res| {
        res.exit_direction.x
    })[0];

    let h = 1e-6;
    let exit_x = |c: f64| singlet(c).trace(&ray, ctx).exit_direction.x;
    let finite_difference = (exit_x(c1 + h) - exit_x(c1 - h)) / (2.0 * h);

    println!("d(exit_x)/d(front curvature) at c1 = {c1} 1/m:");
    println!("  dual-number gradient : {analytic:.6e}");
    println!("  central difference   : {finite_difference:.6e}");
    println!(
        "  relative difference  : {:.2e}",
        ((analytic - finite_difference) / finite_difference).abs()
    );
}
