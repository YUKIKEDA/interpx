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

// --- Input checks ---------------------------------------------------------------

#[test]
fn the_constructors_report_each_fault_in_the_order_of_the_design() {
    use interpx::{DegeneracyKind, InterpolateError as E};
    let points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
    type Build = fn(&[[f64; 2]], &[f64]) -> Result<Cubic, E>;
    let constructors: [Build; 2] = [Cubic::new, Cubic::rescaled];
    for build in constructors {
        let err = |points: &[[f64; 2]], values: &[f64]| build(points, values).err();
        assert_eq!(err(&points, &[1.0, 2.0, 3.0]), None);
        assert_eq!(
            err(&points, &[1.0, 2.0]),
            Some(E::LengthMismatch {
                expected: 3,
                actual: 2
            })
        );
        assert_eq!(
            err(&points[..2], &[1.0, 2.0]),
            Some(E::InsufficientPoints {
                required: 3,
                actual: 2
            })
        );
        let bad = [[0.0, 0.0], [f64::NAN, 0.0], [0.0, f64::INFINITY]];
        assert_eq!(
            err(&bad, &[f64::NAN, 2.0, 3.0]),
            Some(E::NonFinite { index: 1 })
        );
        assert_eq!(
            err(&points, &[1.0, f64::NEG_INFINITY, f64::NAN]),
            Some(E::NonFinite { index: 1 })
        );
        let line = [[0.0, 0.0], [1.0, 1.0], [2.0, 2.0], [0.5, 0.5]];
        assert_eq!(
            err(&line, &[1.0; 4]),
            Some(E::DegenerateGeometry(DegeneracyKind::LowerDimensional {
                affine_dim: 1
            }))
        );
    }
}

// --- The rules ------------------------------------------------------------------

/// A small generator, so that the points are the same on every run.
fn generator(seed: u64) -> impl FnMut() -> f64 {
    let mut state = seed;
    move || {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (state >> 11) as f64 / (1_u64 << 53) as f64
    }
}

#[test]
fn a_linear_function_comes_back_with_its_gradient_and_data_values_at_the_points() {
    let mut next = generator(5);
    let points: Vec<[f64; 2]> = (0..200).map(|_| [next(), next()]).collect();
    let f = |p: &[f64; 2]| 0.5 - 2.0 * p[0] + 3.0 * p[1];
    let values: Vec<f64> = points.iter().map(f).collect();
    let cubic = Cubic::new(&points, &values).expect("valid");
    // The estimate gives the gradient of a plane: its energy is then zero.
    for gradient in cubic.gradients() {
        assert!((gradient[0] + 2.0).abs() < 1e-10 && (gradient[1] - 3.0).abs() < 1e-10);
    }
    for _ in 0..500 {
        let q = [next(), next()];
        if let Some(got) = cubic.eval(&q) {
            assert!((got - f(&q)).abs() < 1e-12, "{q:?}");
        }
    }
    // Any data at the points comes back there.
    let wavy: Vec<f64> = points.iter().map(|p| (5.0 * p[0]).sin() * p[1]).collect();
    let cubic = Cubic::new(&points, &wavy).expect("valid");
    for (p, v) in points.iter().zip(&wavy) {
        assert!((cubic.eval(p).expect("a data point") - v).abs() < 1e-13);
    }
}

#[test]
fn each_column_is_solved_as_it_would_be_alone() {
    // Columns that the conjugate gradients solve in different numbers of
    // steps (a plane at once, a smooth surface, noise), and one that is
    // constant, whose system has nothing to solve: every column comes out bit
    // for bit as it does alone, and the constant one has no gradient.
    let mut next = generator(11);
    let points: Vec<[f64; 2]> = (0..300).map(|_| [next(), next()]).collect();
    let columns: [Vec<f64>; 4] = [
        points.iter().map(|p| 1.0 + p[0] - 2.0 * p[1]).collect(),
        points
            .iter()
            .map(|p| (6.0 * p[0]).sin() * (4.0 * p[1]).cos())
            .collect(),
        (0..points.len()).map(|_| next()).collect(),
        vec![7.0; points.len()],
    ];
    let values: Vec<[f64; 4]> = (0..points.len())
        .map(|i| core::array::from_fn(|c| columns[c][i]))
        .collect();
    let together = Cubic::new(&points, &values).expect("valid");
    for (c, column) in columns.iter().enumerate() {
        let alone = Cubic::new(&points, column).expect("valid");
        for (both, one) in together.gradients().iter().zip(alone.gradients()) {
            for k in 0..2 {
                assert_eq!(both[k][c].to_bits(), one[k].to_bits(), "column {c}");
            }
        }
    }
    assert!(together
        .gradients()
        .iter()
        .all(|g| g[0][3] == 0.0 && g[1][3] == 0.0));
}

#[test]
fn values_scaled_by_any_power_of_ten_give_the_same_interpolant_scaled() {
    // The interpolant is linear in the values, so values times `scale` give
    // the interpolant times `scale`, whatever the scale, as long as it is a
    // finite number. The gradient estimate compares sums of squares with a
    // threshold: before the right-hand side was brought near one, values
    // near 1e±160 overflowed or underflowed those squares, and the
    // gradients came back zero with `Ok` (2.4e-2 off here) or the system
    // was not solved (the review of #187).
    let mut next = generator(13);
    let points: Vec<[f64; 2]> = (0..60).map(|_| [next(), next()]).collect();
    let f = |p: &[f64; 2]| (3.0 * p[0]).sin() * p[1] + 0.5 * p[0] * p[0];
    let values: Vec<f64> = points.iter().map(f).collect();
    let plain = Cubic::new(&points, &values).expect("valid");
    let queries: Vec<[f64; 2]> = (0..200).map(|_| [next(), next()]).collect();
    for scale in [1e-300, 1e-200, 1e-160, 1e-150, 1e150, 1e160, 1e200, 1e300] {
        let scaled: Vec<f64> = values.iter().map(|v| v * scale).collect();
        let cubic = Cubic::new(&points, &scaled).expect("solved at every scale");
        for q in &queries {
            if let (Some(a), Some(b)) = (plain.eval(q), cubic.eval(q)) {
                assert!(
                    (b / scale - a).abs() < 1e-12,
                    "{scale:e} {q:?}: {a} {}",
                    b / scale
                );
            }
        }
    }
    // Coordinates times `scale`: the same interpolant at the scaled queries
    // (the gradients go as one over the scale). Coordinates near 1e±100 gave
    // zero gradients or no solution before the system was brought near one.
    for scale in [1e-300, 1e-160, 1e-100, 1e100, 1e160, 1e300] {
        let moved: Vec<[f64; 2]> = points
            .iter()
            .map(|p| [p[0] * scale, p[1] * scale])
            .collect();
        let cubic = Cubic::new(&moved, &values).expect("solved at every scale");
        for q in &queries {
            if let (Some(a), Some(b)) = (plain.eval(q), cubic.eval(&[q[0] * scale, q[1] * scale])) {
                assert!((b - a).abs() < 1e-12, "{scale:e} {q:?}: {a} {b}");
            }
        }
    }
}

#[test]
fn a_system_that_overflows_is_an_error_and_not_zero_gradients() {
    // Two points 1e-200 apart in a unit square: the weight of their edge,
    // which goes as one over its length, overflows the sums of the
    // iteration. That is an error, not gradients left at zero with `Ok`.
    let mut next = generator(17);
    let mut points: Vec<[f64; 2]> = (0..30).map(|_| [next(), next()]).collect();
    points.extend([[0.0, 0.0], [1e-200, 1e-200]]);
    let values: Vec<f64> = points.iter().map(|p| p[0] + 2.0 * p[1]).collect();
    assert!(matches!(
        Cubic::new(&points, &values),
        Err(interpx::InterpolateError::SolverFailed(
            interpx::SolverError::NotConverged
        ))
    ));
}

#[test]
fn the_interpolant_and_its_gradient_are_continuous_across_edges() {
    // C1: across an edge the difference of the one-sided slopes goes to zero
    // with the step (it is the second derivative times the step), where the
    // linear interpolant's stays (its slope jumps there). At the midpoint of
    // each inner edge, along the normal, steps of 1e-5 and 1e-6. Near the
    // hull the second derivative of a thin triangle's patch reached 1e4: the
    // difference there went 1.02, 0.10, 0.010, 0.001 for steps of 1e-4 down
    // to 1e-7.
    let mut next = generator(9);
    let points: Vec<[f64; 2]> = (0..150).map(|_| [next(), next()]).collect();
    let values: Vec<f64> = points
        .iter()
        .map(|p| (4.0 * p[0]).sin() * (3.0 * p[1]).cos())
        .collect();
    let cubic = Cubic::new(&points, &values).expect("valid");
    let linear = interpx::LinearNDInterpolator::new(&points, &values).expect("valid");
    let mut edges = 0;
    for i in 0..points.len() {
        for j in i + 1..points.len() {
            let (a, b) = (points[i], points[j]);
            let m = [(a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0];
            let e = [b[0] - a[0], b[1] - a[1]];
            let l = (e[0] * e[0] + e[1] * e[1]).sqrt();
            let n = [-e[1] / l, e[0] / l];
            let at = |s: f64| [m[0] + s * n[0], m[1] + s * n[1]];
            // The slopes on the two sides differ by this, at a step.
            let jump = |f: &dyn Fn([f64; 2]) -> Option<f64>, step: f64| {
                let (c, p, q) = (f(at(0.0))?, f(at(step))?, f(at(-step))?);
                Some(((p - c) - (c - q)).abs() / step)
            };
            let lin = |x: [f64; 2]| linear.eval(&x);
            let (Some(kink), Some(kink_fine)) = (jump(&lin, 1e-5), jump(&lin, 1e-6)) else {
                continue;
            };
            // Only on an edge: where the linear interpolant's slope jumps by
            // the same at both steps (a midpoint that is not on the edge is
            // straddled by one step and not the other).
            if kink < 1e-2 || (kink_fine - kink).abs() > 0.1 * kink {
                continue;
            }
            let cub = |x: [f64; 2]| cubic.eval(&x);
            let coarse = jump(&cub, 1e-5).expect("inside");
            let fine = jump(&cub, 1e-6).expect("inside");
            assert!(fine <= 0.2 * coarse + 1e-6, "{m:?}: {coarse} then {fine}");
            edges += 1;
        }
    }
    assert!(edges > 50, "{edges}");
}

#[test]
fn points_given_more_than_once_count_once_and_queries_outside_have_no_value() {
    let points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 0.0], [1.0, 1.0]];
    let cubic = Cubic::new(&points, &[1.0, 2.0, 3.0, 20.0, 4.0]).expect("valid");
    assert_eq!(cubic.len(), 5);
    assert!((cubic.eval(&[1.0, 0.0]).expect("a point") - 2.0).abs() < 1e-14);
    // The point given again has no gradient of its own.
    assert_eq!(cubic.gradients()[3], [0.0, 0.0]);
    for query in [
        [1.5, 0.5],
        [-0.1, 0.5],
        [f64::NAN, 0.5],
        [0.5, f64::INFINITY],
    ] {
        assert_eq!(cubic.eval(&query), None, "{query:?}");
    }
    let queries = [[0.5, 0.5], [2.0, 2.0]];
    let mut out = [0.0; 2];
    assert_eq!(cubic.eval_batch_into(&queries, &mut out, -1.0), Ok(()));
    assert_eq!(out[1], -1.0);
    fn assert_send_sync<T: Send + Sync + 'static>() {}
    assert_send_sync::<Cubic>();
    assert_send_sync::<Cubic<[f64; 3]>>();
}

#[test]
fn the_order_of_convergence_is_that_of_the_estimated_gradients() {
    // With the exact gradients of the function in place of the estimate,
    // SciPy's patch converged at orders 3.25, 3.66, 3.69, 3.05 for h = 1/20
    // to 1/160 on such points: the estimate is what sets the order. On Franke's
    // first function, a lattice of spacing h with its inner points moved by
    // up to 0.15 h, the largest error inside the square went 2.06e-2,
    // 3.17e-3, 5.21e-4, 1.55e-4, 2.79e-5 for h = 1/10 to 1/160: orders 2.70,
    // 2.61, 1.75, 2.47, 2.38 on average. SciPy iterated to convergence gives
    // the same errors within 15% (2.06e-2 to 3.02e-5).
    use common::functions::franke;
    let f = |p: [f64; 2]| franke(1, p[0], p[1]).expect("F1");
    let halton = |i: usize, base: usize| {
        let (mut x, mut share, mut i) = (0.0, 1.0, i);
        while i > 0 {
            share /= base as f64;
            x += share * (i % base) as f64;
            i /= base;
        }
        x
    };
    let error = |m: u32| {
        let h = 1.0 / f64::from(m);
        let mut points = Vec::new();
        for i in 0..=m {
            for j in 0..=m {
                let k = (i * (m + 1) + j + 1) as usize;
                let s = [halton(k, 2), halton(k, 3)].map(|u| (u - 0.5) * 0.3 * h);
                let moved = |c: u32, s: f64| if c == 0 || c == m { 0.0 } else { s };
                points.push([
                    f64::from(i) * h + moved(i, s[0]),
                    f64::from(j) * h + moved(j, s[1]),
                ]);
            }
        }
        let values: Vec<f64> = points.iter().map(|p| f(*p)).collect();
        let cubic = Cubic::from_owned(points, values).expect("valid");
        let mut worst = 0.0_f64;
        for i in 0..=40 {
            for j in 0..=40 {
                let p = [i, j].map(|q| 0.1 + 0.8 * f64::from(q) / 40.0);
                worst = worst.max((cubic.eval(&p).expect("inside") - f(p)).abs());
            }
        }
        worst
    };
    let (coarse, fine) = (error(10), error(160));
    let order = (coarse / fine).log2() / 4.0;
    assert!(
        (2.2..=2.6).contains(&order),
        "{coarse:.2e} {fine:.2e}: {order}"
    );
}
