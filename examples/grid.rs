//! Values on a regular grid in two dimensions: the same data by four methods.
//!
//! `cargo run --example grid`

use interpx::{GridMethod, InterpolateError, Interpolator, RegularGridInterpolator};

fn main() -> Result<(), InterpolateError> {
    let x: [f64; 6] = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
    let y: [f64; 6] = [0.0, 0.5, 1.5, 3.0, 5.0, 7.5];
    // In C order, the last axis fastest, as NumPy lays out `values[i, j]`.
    let values: Vec<f64> = x
        .iter()
        .flat_map(|xi| y.iter().map(move |yj| (xi * 0.7).sin() * (yj * 0.3).cos()))
        .collect();

    let point = [2.3, 4.1];
    for method in [
        GridMethod::Nearest,
        GridMethod::Linear,
        GridMethod::Cubic,
        GridMethod::Pchip,
    ] {
        let grid = RegularGridInterpolator::new([&x, &y], &values, method)?;
        println!("{method:?}: {:?}", grid.eval(&point));
    }
    // Outside the grid there is no value until an extrapolation is chosen.
    let grid = RegularGridInterpolator::new([&x, &y], &values, GridMethod::Linear)?;
    println!("outside: {:?}", grid.eval(&[6.0, 1.0]));
    Ok(())
}
