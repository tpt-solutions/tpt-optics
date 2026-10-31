//! Tests of the bounded grid checker itself.

use tpt_optics_verify::bounded::{check, Axis};

#[test]
fn axis_samples_include_both_endpoints_exactly() {
    let a = Axis::new(-1.0, 1.0, 5);
    assert_eq!(a.sample(0), -1.0);
    assert_eq!(a.sample(2), 0.0);
    assert_eq!(a.sample(4), 1.0);
}

#[test]
fn single_sample_axis_is_its_lower_end() {
    assert_eq!(Axis::new(3.0, 7.0, 1).sample(0), 3.0);
}

#[test]
fn checks_the_full_product_grid() {
    let cases = check([Axis::new(0.0, 1.0, 4), Axis::new(0.0, 1.0, 3)], |_| Ok(()));
    assert_eq!(cases, Ok(12));
}

#[test]
fn a_point_with_no_axes_is_checked_once() {
    assert_eq!(check([], |_: [f64; 0]| Ok(())), Ok(1));
}

#[test]
fn reports_the_first_failing_point_in_enumeration_order() {
    // Fails only when x is above 0.5. The first axis varies fastest, so the
    // first failing point is x = 1 with y = 0.
    let err = check([Axis::new(0.0, 1.0, 3), Axis::new(0.0, 1.0, 2)], |[x, _]| {
        if x > 0.5 {
            Err(format!("x = {x} too large"))
        } else {
            Ok(())
        }
    })
    .unwrap_err();
    assert_eq!(err.point, [1.0, 0.0]);
    assert_eq!(err.message, "x = 1 too large");
}

#[test]
fn an_edge_only_failure_is_found_because_endpoints_are_sampled() {
    // A bug that only shows at the top of the domain.
    let err = check([Axis::new(0.0, 10.0, 6)], |[x]| {
        if x == 10.0 {
            Err(String::from("edge"))
        } else {
            Ok(())
        }
    });
    assert_eq!(err.unwrap_err().point, [10.0]);
}

#[test]
#[should_panic(expected = "each axis needs at least one sample")]
fn a_zero_step_axis_is_rejected() {
    let _ = check([Axis::new(0.0, 1.0, 0)], |_| Ok(()));
}
