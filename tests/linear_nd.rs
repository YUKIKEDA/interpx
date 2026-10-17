//! `LinearNDInterpolator`: SciPy's values, the exact reproduction of affine
//! functions, points given more than once, the hull, rescaled coordinates, the
//! input checks and their order, and the properties.

#![cfg(feature = "scattered")]

mod common;

use common::{load, reals, Real};
use interpx::{
    DegeneracyKind, InterpolateError as E, Interpolator, LinearNDInterpolator as Linear,
    ParamReason,
};
use proptest::prelude::*;
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

/// The golden sets whose points include four or more on one circle: a grid
/// with points taken out. There Qhull and `convx` may split a cell otherwise,
/// and the values inside it differ (docs/design.ja.md section 6.2).
const COCIRCULAR: &[&str] = &["linear/grid_holes_10/f1"];

/// How far a value may be from SciPy's, relative to the largest magnitude of
/// the values of the case: both sum the same three or four products, from
/// transforms inverted by different eliminations. The largest seen over the
/// golden cases is below 1e-14 (see the test).
const TOLERANCE: f64 = 1e-13;

// Marked as test code so that `expect` is allowed here as in a `#[test]` function
// (clippy.toml allows it in tests only).
#[cfg(test)]
fn check<const D: usize>(case: &Json) -> (usize, usize, f64) {
    let name = case["name"].as_str().expect("name");
    let points: Vec<[f64; D]> = rows(&case["points"])
        .iter()
        .map(|r| r.as_slice().try_into().expect("D coordinates"))
        .collect();
    let queries: Vec<[f64; D]> = rows(&case["queries"])
        .iter()
        .map(|r| r.as_slice().try_into().expect("D coordinates"))
        .collect();
    let values = nums(&case["values"]);
    let rescale = case["options"]["rescale"].as_bool().unwrap_or(false);
    let linear = if rescale {
        Linear::rescaled(&points, &values)
    } else {
        Linear::new(&points, &values)
    }
    .expect("valid points");
    assert_eq!(linear.len(), points.len(), "{name}");
    let scale = values.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    let expected = nums(&case["expected"]);
    let mut batch = vec![0.0; queries.len()];
    assert_eq!(
        linear.eval_batch_into(&queries, &mut batch, f64::NAN),
        Ok(())
    );
    let (mut apart, mut outside, mut worst) = (0, 0, 0.0_f64);
    for (k, (query, want)) in queries.iter().zip(&expected).enumerate() {
        let what = format!("{name} query {k} {query:?}");
        let got = linear.eval(query);
        assert_eq!(
            batch[k].to_bits(),
            got.unwrap_or(f64::NAN).to_bits(),
            "{what}"
        );
        if want.is_nan() {
            assert_eq!(got, None, "{what}: outside the hull for SciPy");
            outside += 1;
            continue;
        }
        let got = got.unwrap_or_else(|| panic!("{what}: no value, SciPy has {want}"));
        let error = (got - want).abs() / scale;
        if error > TOLERANCE {
            assert!(
                COCIRCULAR.contains(&name),
                "{what}: {got} against SciPy's {want}"
            );
            apart += 1;
        } else {
            worst = worst.max(error);
        }
    }

    if let Some(columns) = case.get("values_columns") {
        let pairs: Vec<[f64; 2]> = rows(columns).iter().map(|r| [r[0], r[1]]).collect();
        let pair = Linear::new(&points, &pairs).expect("valid points");
        let plain = Linear::new(&points, &values).expect("valid points");
        for (query, want) in queries.iter().zip(rows(&case["expected_columns"])) {
            match pair.eval(query) {
                None => assert!(want[0].is_nan() && want[1].is_nan(), "{name}"),
                Some([a, b]) => {
                    // The first column is the interpolant of the first values.
                    assert_eq!(Some(a), plain.eval(query), "{name}");
                    assert!((a - want[0]).abs() <= TOLERANCE, "{name}");
                    assert!((b - want[1]).abs() <= TOLERANCE, "{name}");
                }
            }
        }
    }
    (apart, outside, worst)
}

// --- Golden data -------------------------------------------------------------

#[test]
fn scipy_linear_interpolants_evaluate_to_the_golden_values() {
    let file = load::<Json>("scattered").expect("golden file scattered");
    let mut cases = [0; 5];
    let (mut apart, mut outside, mut worst) = (0, 0, 0.0_f64);
    for case in &file.cases {
        let name = case["name"].as_str().expect("name");
        if !name.starts_with("linear/") {
            continue;
        }
        let dimension = rows(&case["points"])[0].len();
        let (a, o, w) = match dimension {
            2 => check::<2>(case),
            3 => check::<3>(case),
            4 => check::<4>(case),
            other => panic!("points in {other} dimensions"),
        };
        println!("{name}: {a} apart, {o} outside, worst {w:.1e}");
        cases[dimension] += 1;
        (apart, outside, worst) = (apart + a, outside + o, worst.max(w));
    }
    assert_eq!(cases, [0, 0, 16, 4, 1]);
    println!("apart {apart}, outside {outside}, worst {worst:.1e}");
    assert!(worst < 1e-14, "{worst}");
}

#[test]
fn the_cocircular_set_differs_from_scipy_only_inside_its_cells_and_reproduces_affine_values() {
    // The grid with holes of the golden data: where SciPy's values differ,
    // both are still linear on their own cells, so an affine function comes
    // back exactly from both.
    let file = load::<Json>("scattered").expect("golden file scattered");
    let case = file
        .cases
        .iter()
        .find(|c| c["name"] == COCIRCULAR[0])
        .expect("the grid case");
    let points: Vec<[f64; 2]> = rows(&case["points"]).iter().map(|r| [r[0], r[1]]).collect();
    let queries: Vec<[f64; 2]> = rows(&case["queries"])
        .iter()
        .map(|r| [r[0], r[1]])
        .collect();
    let affine = |p: &[f64; 2]| 0.5 - 2.0 * p[0] + 3.0 * p[1];
    let values: Vec<f64> = points.iter().map(affine).collect();
    let linear = Linear::new(&points, &values).expect("valid points");
    let mut inside = 0;
    for query in &queries {
        if let Some(got) = linear.eval(query) {
            assert!((got - affine(query)).abs() < 1e-14, "{query:?}");
            inside += 1;
        }
    }
    assert_eq!(inside, 54);
}

// --- Input checks ---------------------------------------------------------------

#[test]
fn the_constructors_report_each_fault_in_the_order_of_the_design() {
    let points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
    type Build = fn(&[[f64; 2]], &[f64]) -> Result<Linear<2>, E>;
    let constructors: [Build; 2] = [Linear::new, Linear::rescaled];
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
        // The length before the count.
        assert_eq!(
            err(&points[..2], &[1.0]),
            Some(E::LengthMismatch {
                expected: 2,
                actual: 1
            })
        );
        // A coordinate before a value, and the first of each.
        let bad = [[0.0, 0.0], [f64::NAN, 0.0], [0.0, f64::INFINITY]];
        assert_eq!(
            err(&bad, &[f64::NAN, 2.0, 3.0]),
            Some(E::NonFinite { index: 1 })
        );
        assert_eq!(
            err(&points, &[1.0, f64::NEG_INFINITY, f64::NAN]),
            Some(E::NonFinite { index: 1 })
        );
        // Enough points, too few of them distinct: the count of the distinct.
        let repeated = [[0.0, 0.0], [1.0, 0.0], [0.0, 0.0], [1.0, 0.0]];
        assert_eq!(
            err(&repeated, &[1.0; 4]),
            Some(E::InsufficientPoints {
                required: 3,
                actual: 2
            })
        );
        // On one line: the dimension they span.
        let line = [[0.0, 0.0], [1.0, 1.0], [2.0, 2.0], [0.5, 0.5]];
        assert_eq!(
            err(&line, &[1.0; 4]),
            Some(E::DegenerateGeometry(DegeneracyKind::LowerDimensional {
                affine_dim: 1
            }))
        );
    }
    // In space, on one plane.
    let plane = [
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 1.0],
        [0.0, 1.0, 1.0],
        [1.0, 1.0, 1.0],
    ];
    assert_eq!(
        Linear::new(&plane, &[1.0; 4]).err(),
        Some(E::DegenerateGeometry(DegeneracyKind::LowerDimensional {
            affine_dim: 2
        }))
    );
    // Rescaled: a mean that overflows is refused.
    let huge = 1.7e308;
    let same = [[huge, 0.0], [huge, 1.0], [huge * 0.5, 2.0]];
    assert_eq!(
        Linear::rescaled(&same, &[1.0, 2.0, 3.0]).err(),
        Some(E::NonFinite { index: 1 })
    );

    // More points than the triangulation can number. Points without
    // coordinates and values without columns take no memory.
    #[cfg(target_pointer_width = "64")]
    {
        let count = u32::MAX as usize + 1;
        let (points, values) = (vec![[0.0_f64; 0]; count], vec![[0.0_f64; 0]; count]);
        assert_eq!(
            Linear::new(&points, &values).err(),
            Some(E::InvalidParameter {
                name: "points",
                reason: ParamReason::TooLarge {
                    max: u32::MAX as usize
                },
            })
        );
    }
}

// --- The rules ------------------------------------------------------------------

#[test]
fn points_given_more_than_once_count_once_with_the_value_given_first() {
    let points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 0.0], [0.0, 0.0]];
    let linear = Linear::new(&points, &[1.0, 2.0, 3.0, 20.0, 10.0]).expect("valid");
    assert_eq!(linear.len(), 5);
    assert_eq!(linear.simplex_count(), 1);
    assert_eq!(linear.eval(&[0.0, 0.0]), Some(1.0));
    assert_eq!(linear.eval(&[1.0, 0.0]), Some(2.0));
    let middle = linear.eval(&[0.5, 0.0]).expect("on an edge");
    assert!((middle - 1.5).abs() < 1e-15);
}

#[test]
fn queries_outside_the_hull_or_not_finite_have_no_value_and_a_batch_fills_them() {
    let points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]];
    let linear = Linear::new(&points, &[1.0, 2.0, 3.0, 4.0]).expect("valid");
    // On the edges and corners of the hull: a value.
    for query in [[0.5, 0.0], [1.0, 0.5], [0.0, 0.0], [1.0, 1.0]] {
        assert!(linear.eval(&query).is_some(), "{query:?}");
    }
    // Outside by rounding: still a value. Outside by more: none.
    assert!(linear.eval(&[1.0 + 1e-15, 0.5]).is_some());
    for query in [
        [1.0 + 1e-9, 0.5],
        [-0.1, 0.5],
        [0.5, 2.0],
        [1e300, 0.0],
        [f64::NAN, 0.5],
        [0.5, f64::INFINITY],
    ] {
        assert_eq!(linear.eval(&query), None, "{query:?}");
    }
    let queries = [[0.5, 0.5], [2.0, 2.0], [f64::NAN, 0.0]];
    let mut out = [0.0; 3];
    assert_eq!(linear.eval_batch_into(&queries, &mut out, -1.0), Ok(()));
    assert!((out[0] - 2.5).abs() < 1e-15);
    assert_eq!(&out[1..], &[-1.0, -1.0]);
    assert_eq!(
        linear.eval_batch_into(&queries, &mut out[..2], -1.0),
        Err(E::LengthMismatch {
            expected: 3,
            actual: 2
        })
    );
}

#[test]
fn handed_over_and_rescaled_give_the_interpolants_they_say() {
    let points: Vec<[f64; 2]> = (0..40)
        .map(|i| {
            let t = f64::from(i);
            [(t * 0.618_034).fract(), 1000.0 * (t * 0.414_214).fract()]
        })
        .collect();
    let values: Vec<f64> = points
        .iter()
        .map(|p| (p[0] * 3.0).sin() + p[1] / 500.0)
        .collect();
    let queries: Vec<[f64; 2]> = (0..50)
        .map(|i| {
            let t = f64::from(i) + 0.5;
            [(t * 0.3819).fract(), 1000.0 * (t * 0.7317).fract()]
        })
        .collect();
    let copied = Linear::new(&points, &values).expect("valid");
    let owned = Linear::from_owned(points.clone(), values.clone()).expect("valid");
    let rescaled = Linear::rescaled(&points, &values).expect("valid");
    let rescaled_owned =
        Linear::rescaled_from_owned(points.clone(), values.clone()).expect("valid");
    // Rescaled is new on the points brought to one range.
    let low = [0, 1].map(|d| points.iter().map(|p| p[d]).fold(f64::INFINITY, f64::min));
    let high = [0, 1].map(|d| {
        points
            .iter()
            .map(|p| p[d])
            .fold(f64::NEG_INFINITY, f64::max)
    });
    let mean = [0, 1].map(|d| points.iter().map(|p| p[d]).sum::<f64>() / 40.0);
    let change = |p: &[f64; 2]| [0, 1].map(|d| (p[d] - mean[d]) / (high[d] - low[d]));
    let moved: Vec<[f64; 2]> = points.iter().map(change).collect();
    let on_moved = Linear::new(&moved, &values).expect("valid");
    let mut inside = 0;
    for query in &queries {
        assert_eq!(copied.eval(query), owned.eval(query));
        assert_eq!(rescaled.eval(query), rescaled_owned.eval(query));
        assert_eq!(rescaled.eval(query), on_moved.eval(&change(query)));
        inside += usize::from(rescaled.eval(query).is_some());
    }
    assert!(inside > 30, "{inside}");
}

#[test]
fn no_dimensions_one_dimension_columns_and_thread_safety() {
    // Without dimensions: the first value, everywhere.
    let none = Linear::<0>::new(&[[], [], []], &[4.0, 5.0, 6.0]).expect("valid");
    assert_eq!(none.eval(&[]), Some(4.0));
    assert_eq!(none.simplex_count(), 0);
    // One dimension: linear interpolation between neighboring points, in any
    // order given.
    let line = Linear::new(&[[3.0], [0.0], [1.0]], &[5.0, 0.0, 1.0]).expect("valid");
    assert_eq!(line.simplex_count(), 2);
    let between = |x: f64| line.eval(&[x]).expect("inside");
    assert!((between(0.5) - 0.5).abs() < 1e-15);
    assert!((between(2.0) - 3.0).abs() < 1e-15);
    assert_eq!(line.eval(&[3.5]), None);
    // Two values at every point: each column as alone.
    let points = [[0.0, 0.0], [2.0, 0.0], [0.0, 2.0]];
    let pairs = Linear::new(&points, &[[0.0, 1.0], [2.0, 1.0], [4.0, 1.0]]).expect("valid");
    let [a, b] = pairs.eval(&[0.5, 0.5]).expect("inside");
    assert!((a - 1.5).abs() < 1e-15 && (b - 1.0).abs() < 1e-15);

    fn assert_send_sync<T: Send + Sync + 'static>() {}
    assert_send_sync::<Linear<3>>();
    assert_send_sync::<Linear<2, [f64; 4]>>();
}

#[test]
fn a_thin_simplex_is_looked_through_to_its_neighbors() {
    // The triangle along the bottom edge is 1e-13 high: the reciprocal of the
    // condition number of its matrix is 2e-13, below SciPy's limit of
    // 1000 * EPSILON, so it is degenerate. A query inside it, 2.5e-14 below
    // the edge of the triangle to its left, is outside that triangle by more
    // than EPS and less than EPS_BROAD: found there by the wider leeway toward
    // the degenerate one.
    let points = [[0.0, 0.0], [1.0, 0.0], [0.5, 1e-13], [0.5, 1.0]];
    // f(x, y) = x + 1.5 y at each point.
    let values = [0.0, 1.0, 0.5 + 1.5e-13, 2.0];
    let linear = Linear::new(&points, &values).expect("valid");
    assert_eq!(linear.simplex_count(), 3);
    let inside = linear
        .eval(&[0.25, 2.5e-14])
        .expect("found through the sliver");
    assert!((inside - 0.25).abs() < 1e-12, "{inside}");
    for query in [[0.5, 5e-14], [0.25, 0.0], [0.75, 1e-14], [0.5, -1e-12]] {
        if let Some(value) = linear.eval(&query) {
            assert!(
                (value - query[0] - 1.5 * query[1]).abs() < 1e-12,
                "{query:?}"
            );
        }
    }
    let above = linear.eval(&[0.5, 0.5]).expect("inside");
    assert!((above - 1.25).abs() < 1e-12, "{above}");
}

#[test]
fn the_data_values_come_back_at_the_vertices_of_a_thin_tetrahedron() {
    // Found by the property test below on CI: a data point is a vertex of a
    // tetrahedron with a condition number of 10^4. With the inverse of its
    // matrix the value there came out 1.5e-12 off (SciPy, which finds the
    // point in a well-conditioned tetrahedron, 1.8e-15); with its LU
    // decomposition the coordinates of a vertex come out as they are.
    let mut points = vec![
        [0.042693381611414366, 0.3890558422117542, 0.7940582959686024],
        [0.6125880437154237, 0.9472302659307226, 0.8112851235989772],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.7785790445215961, 0.9893039824225027, 0.8887853078116238],
        [0.28942934674679066, 0.7174686133446314, 0.0],
        [0.7553268439579649, 0.6297163868201694, 0.0],
        [0.29863708187555077, 0.0, 0.0],
    ];
    points.extend([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ]);
    let f = |p: &[f64; 3]| 4.920027376003285 * p[0] + 3.343642051541309 * p[1];
    let values: Vec<f64> = points.iter().map(f).collect();
    let linear = Linear::new(&points, &values).expect("valid");
    for p in &points {
        let got = linear.eval(p).expect("a data point is in the hull");
        assert!((got - f(p)).abs() < 1e-14, "{p:?}: {got} against {}", f(p));
    }
}

// --- Properties -----------------------------------------------------------------

proptest! {
    #![proptest_config(ProptestConfig { cases: 64, ..ProptestConfig::default() })]

    /// An affine function comes back exactly (to rounding) inside the hull,
    /// and the data values come back at the data points, for points in
    /// general position and for points on a lattice (four on one circle).
    #[test]
    fn affine_functions_come_back_and_data_values_at_the_points(
        raw in proptest::collection::vec((0.0..1.0_f64, 0.0..1.0_f64, 0.0..1.0_f64), 4..40),
        lattice in any::<bool>(),
        coefficients in (-5.0..5.0_f64, -5.0..5.0_f64, -5.0..5.0_f64, -5.0..5.0_f64),
        query in (0.0..1.0_f64, 0.0..1.0_f64, 0.0..1.0_f64),
    ) {
        let snap = |c: f64| if lattice { (c * 4.0).round() / 4.0 } else { c };
        let mut points: Vec<[f64; 3]> =
            raw.iter().map(|&(x, y, z)| [snap(x), snap(y), snap(z)]).collect();
        // A tetrahedron that spans space, so that the points are never flat.
        points.extend([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
        let (a, b, c, d) = coefficients;
        let f = |p: &[f64; 3]| a + b * p[0] + c * p[1] + d * p[2];
        let values: Vec<f64> = points.iter().map(f).collect();
        let linear = Linear::new(&points, &values).expect("valid");
        let q = [query.0, query.1, query.2];
        if let Some(got) = linear.eval(&q) {
            prop_assert!((got - f(&q)).abs() < 1e-12, "{got} against {}", f(&q));
        }
        for p in &points {
            let got = linear.eval(p).expect("a data point is in the hull");
            prop_assert!((got - f(p)).abs() < 1e-12);
        }
    }
}
