//! Scattered data in two dimensions: the nearest data point, linear and
//! cubic interpolation on the Delaunay triangulation, and a radial basis
//! function interpolant over all points and over the nearest ones.
//!
//! `cargo run --example scattered`

use interpx::{
    CloughTocher2DInterpolator, InterpolateError, Interpolator, LinearNDInterpolator,
    NearestNDInterpolator, RBFInterpolator, RbfKernel, RbfOptions,
};

fn main() -> Result<(), InterpolateError> {
    // 500 points spread evenly over the unit square (a low-discrepancy
    // sequence), and a smooth function of them.
    let f = |[x, y]: [f64; 2]| (3.0 * x).sin() * (2.0 * y).cos();
    let points: Vec<[f64; 2]> = (1..=500)
        .map(|i| {
            let i = f64::from(i);
            [
                (0.5 + i * 0.754_877_666).fract(),
                (0.5 + i * 0.569_840_291).fract(),
            ]
        })
        .collect();
    let values: Vec<f64> = points.iter().copied().map(f).collect();
    let query = [0.4, 0.6];
    println!("exact: {}", f(query));

    let nearest = NearestNDInterpolator::new(&points, &values)?;
    println!("nearest: {:?}", nearest.eval(&query));

    // On the triangles of the points: no value outside their convex hull.
    let linear = LinearNDInterpolator::new(&points, &values)?;
    println!("linear: {:?}", linear.eval(&query));
    let cubic = CloughTocher2DInterpolator::new(&points, &values)?;
    println!("Clough-Tocher: {:?}", cubic.eval(&query));

    // One system over all points: O(N^2) memory, and an estimate of how well
    // it was solved.
    let rbf = RBFInterpolator::new(&points, &values, RbfKernel::ThinPlateSpline)?;
    println!("thin plate spline: {:?}", rbf.eval(&query));
    println!("reciprocal condition: {:?}", rbf.reciprocal_condition());

    // A small system for each query, over its 30 nearest points.
    let options = RbfOptions::new().with_neighbors(30);
    let local =
        RBFInterpolator::with_options(&points, &values, RbfKernel::ThinPlateSpline, &options)?;
    println!("over 30 neighbors: {:?}", local.eval(&query));
    Ok(())
}
