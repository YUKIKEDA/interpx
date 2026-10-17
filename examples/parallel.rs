//! A large batch evaluated on several threads (the `rayon` feature). The
//! results are those of the serial batch, bit for bit.
//!
//! `cargo run --release --example parallel --features rayon`

use interpx::{
    GridMethod, InterpolateError, Interpolator, ParInterpolator, RegularGridInterpolator,
};

fn main() -> Result<(), InterpolateError> {
    let n = 200;
    let axis: Vec<f64> = (0..n).map(f64::from).collect();
    let values: Vec<f64> = (0..n * n).map(|k| f64::from(k % 97)).collect();
    let grid = RegularGridInterpolator::new([&axis, &axis], &values, GridMethod::Cubic)?;

    let queries: Vec<[f64; 2]> = (0..1_000_000)
        .map(|k| {
            let t = f64::from(k) * 1.0e-6;
            [t * 198.0, (1.0 - t) * 198.0]
        })
        .collect();
    let mut serial = vec![0.0; queries.len()];
    let mut parallel = vec![0.0; queries.len()];
    grid.eval_batch_into(&queries, &mut serial, f64::NAN)?;
    grid.eval_batch_into_par(&queries, &mut parallel, f64::NAN)?;
    let same = serial
        .iter()
        .zip(&parallel)
        .all(|(a, b)| a.to_bits() == b.to_bits());
    println!(
        "{} points, the same bits on several threads: {same}",
        queries.len()
    );
    Ok(())
}
