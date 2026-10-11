//! Checks the golden pipeline end to end: the file the generator wrote, the
//! loader, and the tolerance helper.
//!
//! No interpx code exists yet, so the reference here is a small independent
//! linear interpolation written in this test. It stands for `numpy.interp`
//! (`compare/generate_pipeline_check.py`). Real families replace it with the
//! interpolator under test.

mod common;

use common::{assert_close, close, load, reals, Real, DIRECT_F64};
use serde::Deserialize;

#[derive(Debug, Deserialize)]
struct Options {
    left: Option<f64>,
    right: Option<f64>,
}

#[derive(Debug, Deserialize)]
struct Case {
    name: String,
    reference: String,
    options: Options,
    x: Vec<Real>,
    y: Vec<Real>,
    queries: Vec<Real>,
    expected: Vec<Real>,
}

/// `numpy.interp` semantics: `left` below the first knot, `right` above the
/// last, NaN in gives NaN out.
fn linear(x: &[f64], y: &[f64], q: f64, left: f64, right: f64) -> f64 {
    let last = x.len() - 1;
    if q.is_nan() {
        return f64::NAN;
    }
    if q < x[0] {
        return left;
    }
    if q > x[last] {
        return right;
    }
    if q >= x[last] {
        return y[last];
    }
    let i = x.partition_point(|&v| v <= q) - 1;
    let slope = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
    slope * (q - x[i]) + y[i]
}

#[test]
fn file_records_the_baseline() {
    let file = load::<Case>("pipeline_check").expect("golden file pipeline_check");
    assert!(
        file.meta.scipy.starts_with("1.18."),
        "SciPy baseline is 1.18, file says {}",
        file.meta.scipy
    );
    assert_eq!(file.meta.generator, "generate_pipeline_check.py");
    assert!(!file.cases.is_empty());
}

#[test]
fn golden_matches_the_independent_reference() {
    let file = load::<Case>("pipeline_check").expect("golden file pipeline_check");
    for case in &file.cases {
        assert_eq!(case.reference, "numpy.interp", "{}", case.name);
        let (x, y) = (reals(&case.x), reals(&case.y));
        let left = case.options.left.unwrap_or(y[0]);
        let right = case.options.right.unwrap_or(y[y.len() - 1]);
        assert_eq!(case.queries.len(), case.expected.len(), "{}", case.name);
        for (k, (q, want)) in case.queries.iter().zip(&case.expected).enumerate() {
            let got = linear(&x, &y, q.0, left, right);
            assert_close(got, want.0, DIRECT_F64, &format!("{} query {k}", case.name));
        }
    }
}

#[test]
fn file_has_a_nan_query_and_both_out_of_range_sides() {
    let file = load::<Case>("pipeline_check").expect("golden file pipeline_check");
    let case = file
        .cases
        .iter()
        .find(|c| c.name == "small_fill")
        .expect("case small_fill");
    let x = reals(&case.x);
    let q = reals(&case.queries);
    assert!(q.iter().any(|v| v.is_nan()));
    assert!(q.iter().any(|&v| v < x[0]));
    assert!(q.iter().any(|&v| v > x[x.len() - 1]));
}

#[test]
fn tolerance_rejects_a_tampered_value() {
    let tol = DIRECT_F64;
    assert!(close(1.0, 1.0, tol));
    assert!(close(1.0, 1.0 + 1e-13, tol));
    assert!(!close(1.0, 1.0 + 1e-9, tol));
}

#[test]
fn tolerance_handles_non_finite_values() {
    let tol = DIRECT_F64;
    assert!(close(f64::NAN, f64::NAN, tol));
    assert!(!close(f64::NAN, 0.0, tol));
    assert!(!close(0.0, f64::NAN, tol));
    assert!(close(f64::INFINITY, f64::INFINITY, tol));
    assert!(!close(f64::INFINITY, f64::NEG_INFINITY, tol));
    assert!(!close(f64::INFINITY, 1.0e308, tol));
}
