//! `CloughTocher2DInterpolator`: SciPy's values, the reproduction of linear
//! functions, the gradient estimate, the input checks, and the properties.

#![cfg(feature = "scattered")]

mod common;

use common::{load, reals, Real};
use interpx::{CloughTocher2DInterpolator as Cubic, Interpolator};
use serde_json::Value as Json;

fn nums(v: &Json) -> Vec<f64> {
    reals(&serde_json::from_value::<Vec<Real>>(v.clone()).unwrap_or_default())
}

fn rows(v: &Json) -> Vec<Vec<f64>> {
    serde_json::from_value::<Vec<Vec<Real>>>(v.clone())
        .unwrap_or_default()
        .iter()
        .map(|r| reals(r))
        .collect()
}

fn pairs(v: &Json) -> Vec<[f64; 2]> {
    rows(v).iter().map(|r| [r[0], r[1]]).collect()
}

/// The golden sets whose points include four or more on one circle: one
/// diagonal of `franke_33` and the grid with holes. Qhull and `convx`
/// triangulate them otherwise, and since the gradients are the solution of
/// one system over all edges, an edge that differs moves the values
/// everywhere (6e-3 of the largest value on `franke_33`, where one edge of 88
/// differs).
const COCIRCULAR: &[&str] = &[
    "clough_tocher/franke_33/f1",
    "clough_tocher/grid_holes_10/f1",
];

/// How far a value may be from SciPy's, relative to the largest magnitude of
/// the values: when SciPy has iterated to convergence (the `converged`
/// cases; the largest seen is 4.5e-13), and when it stopped at its default
/// `tol` of 1e-6 (2.8e-8). The test prints them.
const CONVERGED: f64 = 1e-12;
const DEFAULT_TOL: f64 = 1e-7;

// Marked as test code so that `expect` is allowed here as in a `#[test]` function
// (clippy.toml allows it in tests only).
#[cfg(test)]
fn check(case: &Json) -> (bool, f64) {
    let name = case["name"].as_str().expect("name");
    let points = pairs(&case["points"]);
    let queries = pairs(&case["queries"]);
    let values = nums(&case["values"]);
    let expected = nums(&case["expected"]);
    let rescale = case["options"]["rescale"].as_bool().unwrap_or(false);
    let cubic = if rescale {
        Cubic::rescaled(&points, &values)
    } else {
        Cubic::new(&points, &values)
    }
    .expect("valid points");
    let cubic = match case["options"]["simplex_tolerance"].as_f64() {
        Some(t) => cubic.with_simplex_tolerance(t).expect("above zero"),
        None => cubic,
    };
    let converged = name.starts_with("clough_tocher/converged/");
    let limit = if converged { CONVERGED } else { DEFAULT_TOL };
    let scale = values.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    let mut batch = vec![0.0; queries.len()];
    assert_eq!(
        cubic.eval_batch_into(&queries, &mut batch, f64::NAN),
        Ok(())
    );
    let mut largest = 0.0_f64;
    for (k, (query, want)) in queries.iter().zip(&expected).enumerate() {
        let got = cubic.eval(query);
        assert_eq!(
            batch[k].to_bits(),
            got.unwrap_or(f64::NAN).to_bits(),
            "{name} {query:?}"
        );
        if want.is_nan() {
            assert_eq!(got, None, "{name} {query:?}");
            continue;
        }
        let got = got.unwrap_or_else(|| panic!("{name} {query:?}: no value"));
        let error = (got - want).abs() / scale;
        largest = largest.max(error);
        if !COCIRCULAR.contains(&name) {
            assert!(error <= limit, "{name} {query:?}: {got} against {want}");
        }
    }
    if let Some(columns) = case.get("values_columns") {
        let both: Vec<[f64; 2]> = pairs(columns);
        let two = Cubic::new(&points, &both).expect("valid points");
        let second = both.iter().map(|v| v[1].abs()).fold(0.0_f64, f64::max);
        for (query, want) in queries.iter().zip(rows(&case["expected_columns"])) {
            match two.eval(query) {
                None => assert!(want[0].is_nan() && want[1].is_nan(), "{name}"),
                Some([a, b]) => {
                    // The first column is the interpolant of the first values.
                    assert_eq!(Some(a), cubic.eval(query), "{name}");
                    let error = (b - want[1]).abs() / second;
                    assert!(error <= limit, "{name}: {b} against {}", want[1]);
                }
            }
        }
    }
    println!("{name}: largest {largest:.1e}");
    (converged, largest)
}

// --- Golden data -------------------------------------------------------------

#[test]
fn scipy_clough_tocher_interpolants_evaluate_to_the_golden_values() {
    let file = load::<Json>("scattered").expect("golden file scattered");
    let (mut converged, mut stopped) = (0.0_f64, 0.0_f64);
    let mut cases = [0, 0];
    for case in &file.cases {
        let name = case["name"].as_str().expect("name");
        if !name.starts_with("clough_tocher/") {
            continue;
        }
        let (is_converged, largest) = check(case);
        cases[usize::from(is_converged)] += 1;
        if COCIRCULAR.contains(&name) {
            continue;
        }
        if is_converged {
            converged = converged.max(largest);
        } else {
            stopped = stopped.max(largest);
        }
    }
    assert_eq!(cases, [11, 10]);
    println!(
        "worst: against converged SciPy {converged:.1e}, against SciPy at tol 1e-6 {stopped:.1e}"
    );
}
