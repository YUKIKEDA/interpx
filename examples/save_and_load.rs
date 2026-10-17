//! An interpolator written as JSON and read back without being built again
//! (the `serde` feature).
//!
//! `cargo run --example save_and_load --features serde`

use interpx::{BSpline1D, BoundaryCondition, Interpolator};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let x = [0.0, 0.5, 1.2, 2.0, 3.1, 4.0];
    let y = [0.0, 0.4, 0.9, 0.9, 0.1, -0.7];
    let spline = BSpline1D::new(&x, &y, 3, BoundaryCondition::NotAKnot)?;

    // What is stored is what it was built into (knots and control points),
    // with a version of the form. Reading it runs the checks of the
    // constructor.
    let text = serde_json::to_string(&spline)?;
    let back: BSpline1D<f64> = serde_json::from_str(&text)?;
    println!("{} bytes of JSON", text.len());
    println!(
        "before: {:?}, after: {:?}",
        spline.eval(&1.7),
        back.eval(&1.7)
    );
    Ok(())
}
