//! Analytic test functions of docs/design.ja.md section 8.4, for the tests that
//! compare with a closed form (convergence order, exact reproduction).
//!
//! These are Rust copies of `compare/common/datasets.py`. `tests/goldens.rs`
//! checks them against values that the Python side wrote to
//! `compare/goldens/functions.json`.

use std::f64::consts::PI;

/// Franke's test function `k` (1 to 6) on the unit square. `None` for another `k`.
///
/// F6 is `sqrt(64 - 81 r^2) / 9 - 1/2`: a sphere of radius 8/9.
pub fn franke(k: u64, x: f64, y: f64) -> Option<f64> {
    let r2 = (x - 0.5).powi(2) + (y - 0.5).powi(2);
    let value = match k {
        1 => {
            0.75 * (-((9.0 * x - 2.0).powi(2) + (9.0 * y - 2.0).powi(2)) / 4.0).exp()
                + 0.75 * (-(9.0 * x + 1.0).powi(2) / 49.0 - (9.0 * y + 1.0).powi(2) / 10.0).exp()
                + 0.5 * (-((9.0 * x - 7.0).powi(2) + (9.0 * y - 3.0).powi(2)) / 4.0).exp()
                - 0.2 * (-(9.0 * x - 4.0).powi(2) - (9.0 * y - 7.0).powi(2)).exp()
        }
        2 => ((9.0 * (y - x)).tanh() + 1.0) / 9.0,
        3 => (1.25 + (5.4 * y).cos()) / (6.0 + 6.0 * (3.0 * x - 1.0).powi(2)),
        4 => (-81.0 / 16.0 * r2).exp() / 3.0,
        5 => (-81.0 / 4.0 * r2).exp() / 3.0,
        6 => (64.0 - 81.0 * r2).sqrt() / 9.0 - 0.5,
        _ => return None,
    };
    Some(value)
}

/// The Genz parameters `a_i = 2 + i / d`. An implementation choice shared with
/// the Python side, not Genz's difficulty scaling.
pub fn genz_a(d: usize) -> Vec<f64> {
    (0..d).map(|i| 2.0 + i as f64 / d as f64).collect()
}

/// A Genz (1984) family at `x`, with `u_i = 1/2`. `None` for an unknown `kind`.
pub fn genz(kind: &str, x: &[f64]) -> Option<f64> {
    let d = x.len();
    let a = genz_a(d);
    let u = 0.5;
    let dot: f64 = a.iter().zip(x).map(|(a, x)| a * x).sum();
    let value = match kind {
        "oscillatory" => (2.0 * PI * u + dot).cos(),
        "product_peak" => a
            .iter()
            .zip(x)
            .map(|(a, x)| 1.0 / (a.powi(-2) + (x - u).powi(2)))
            .product(),
        "corner_peak" => (1.0 + dot).powi(-(d as i32 + 1)),
        "gaussian_peak" => (-a
            .iter()
            .zip(x)
            .map(|(a, x)| a * a * (x - u).powi(2))
            .sum::<f64>())
        .exp(),
        "continuous" => (-a.iter().zip(x).map(|(a, x)| a * (x - u).abs()).sum::<f64>()).exp(),
        _ => return None,
    };
    Some(value)
}

/// The Rosenbrock function: the sum of `100 (x[i+1] - x[i]^2)^2 + (1 - x[i])^2`.
pub fn rosenbrock(x: &[f64]) -> f64 {
    x.windows(2)
        .map(|w| 100.0 * (w[1] - w[0] * w[0]).powi(2) + (1.0 - w[0]).powi(2))
        .sum()
}
