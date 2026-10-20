//! Data in one variable: a cubic spline through five points, evaluated at one
//! point, outside its data, and at a batch.
//!
//! `cargo run --example one_variable`

use interpx::{BoundaryCondition, CubicSpline, Extrapolate, InterpolateError, Interpolator};

fn main() -> Result<(), InterpolateError> {
    let x = [0.0, 1.0, 2.0, 3.5, 5.0];
    let y = [1.0, 2.7, 5.8, 4.0, 2.2];

    // The not-a-knot condition at both ends.
    let spline = CubicSpline::new(&x, &y, BoundaryCondition::NotAKnot)?;
    println!("s(2.5) = {:?}", spline.eval(&2.5));

    // Outside the data the end pieces go on. Without that, a point
    // outside has no value, and `eval_or` gives the fill value instead.
    println!("s(6.0) = {:?}", spline.eval(&6.0));
    let inside_only = spline.with_extrapolate(Extrapolate::Forbid);
    println!(
        "s(6.0) = {:?} without extrapolation",
        inside_only.eval(&6.0)
    );
    println!(
        "s(6.0) = {} with a fill of NaN",
        inside_only.eval_or(&6.0, f64::NAN)
    );

    // A batch into a buffer of the caller's: nothing is allocated.
    let queries = [0.5, 1.5, 4.0, 9.0];
    let mut out = [0.0; 4];
    inside_only.eval_batch_into(&queries, &mut out, f64::NAN)?;
    println!("batch = {out:?}");
    Ok(())
}
