//! `CubicHermiteSpline`: the golden data, the input checks with the slopes
//! through the constructor, what is reproduced of the values and slopes given,
//! the slopes of a `CubicSpline`, and the properties.

mod common;

use common::{assert_close, load, reals, Real, DIRECT_F64};
use interpx::{
    BoundaryCondition, CubicHermiteSpline, CubicSpline, Extrapolate, Hint, InterpolateError as E,
    Interpolator,
};
use proptest::prelude::*;
use serde_json::Value as Json;

fn bits(v: Option<f64>) -> Option<u64> {
    v.map(f64::to_bits)
}

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

/// The value and the slope of piece `i` at distance `dx` from its left end,
/// written apart from the library.
fn piece(c: &[f64], i: usize, dx: f64) -> [f64; 2] {
    let c = &c[4 * i..4 * i + 4];
    [
        ((c[0] * dx + c[1]) * dx + c[2]) * dx + c[3],
        (3.0 * c[0] * dx + 2.0 * c[1]) * dx + c[2],
    ]
}

// --- Golden data -------------------------------------------------------------

/// A golden column as a function of the query.
type Column<'a> = &'a dyn Fn(f64) -> Option<f64>;

#[test]
fn scipy_hermite_splines_are_rebuilt_from_their_data() {
    let file = load::<Json>("oned").expect("golden file oned");
    let (mut cases, mut small) = (0, 0);
    let mut modes = std::collections::BTreeMap::new();
    for case in &file.cases {
        let name = case["name"].as_str().expect("name");
        if case["reference"] != "scipy.interpolate.CubicHermiteSpline"
            || case.get("coefficients").is_none()
        {
            continue;
        }
        let (x, y, dydx) = (nums(&case["x"]), nums(&case["y"]), nums(&case["dydx"]));
        let built = CubicHermiteSpline::new(&x, &y, &dydx).expect("valid data");
        assert_eq!(
            built.extrapolate(),
            Extrapolate::Extend,
            "{name}: SciPy's default"
        );
        let mode = match &case["options"]["extrapolate"] {
            Json::Bool(true) => Extrapolate::Extend,
            Json::Bool(false) => Extrapolate::Forbid,
            Json::String(s) if s == "periodic" => Extrapolate::Periodic,
            other => panic!("{name}: extrapolate = {other}"),
        };
        *modes.entry(format!("{mode:?}")).or_insert(0) += 1;
        let spline = built.with_extrapolate(mode);

        // A coefficient is compared by what it adds to the value at the far end
        // of its interval (`tests/cubic_spline.rs` says why): SciPy stores
        // `c[power][interval]`.
        let scipy = rows(&case["coefficients"]);
        let poly = spline.as_ppoly();
        assert_eq!(poly.breakpoints().len(), x.len(), "{name}");
        for (i, piece) in poly.coefficients().as_chunks::<4>().0.iter().enumerate() {
            let width = x[i + 1] - x[i];
            for (power, got) in piece.iter().enumerate() {
                let reach = width.powi(3 - i32::try_from(power).expect("0 to 3"));
                let what = format!("{name} interval {i} coefficient {power}");
                assert_close(got * reach, scipy[power][i] * reach, DIRECT_F64, &what);
            }
        }

        let queries = nums(&case["queries"]);
        let columns: [(&str, Column); 3] = [
            ("expected", &|q| spline.eval(&q)),
            ("derivative_1", &|q| poly.eval_derivative(q, 1)),
            ("derivative_2", &|q| poly.eval_derivative(q, 2)),
        ];
        for (field, f) in columns {
            for (k, (q, want)) in queries.iter().zip(nums(&case[field])).enumerate() {
                let what = format!("{name} {field} query {k} (x = {q})");
                match f(*q) {
                    None => assert!(want.is_nan(), "{what}: None, SciPy has {want}"),
                    Some(got) => assert_close(got, want, DIRECT_F64, &what),
                }
            }
        }
        for triple in rows(&case["integrals"]) {
            let got = poly
                .integrate(triple[0], triple[1])
                .expect("inside the data");
            assert_close(got, triple[2], DIRECT_F64, &format!("{name} integral"));
        }
        cases += 1;
        small += usize::from(x.len() <= 3);
    }
    // 4 on the standard data, one of them without extrapolation, then 8 chosen
    // for the type: two of them with three points or fewer, two periodic.
    assert_eq!((cases, small), (12, 2));
    let modes: Vec<(String, i32)> = modes.into_iter().collect();
    let want = [("Extend", 9), ("Forbid", 1), ("Periodic", 2)].map(|(m, n)| (m.to_owned(), n));
    assert_eq!(modes, want);
}

#[test]
fn scipy_hermite_of_two_columns_is_rebuilt_column_by_column() {
    let file = load::<Json>("oned").expect("golden file oned");
    let mut cases = 0;
    for case in &file.cases {
        let name = case["name"].as_str().expect("name");
        if case["reference"] != "scipy.interpolate.CubicHermiteSpline"
            || case.get("coefficients").is_some()
        {
            continue;
        }
        let (y, dydx) = (pairs(&case["y"]), pairs(&case["dydx"]));
        let spline = CubicHermiteSpline::new(&nums(&case["x"]), &y, &dydx).expect("valid data");
        let poly = spline.as_ppoly();
        let queries = nums(&case["queries"]);
        for (field, nu) in [("expected", 0), ("derivative_1", 1), ("derivative_2", 2)] {
            for (k, (q, want)) in queries.iter().zip(rows(&case[field])).enumerate() {
                let got = poly.eval_derivative(*q, nu).expect("the case extrapolates");
                for column in 0..2 {
                    let what = format!("{name} {field} query {k} column {column}");
                    assert_close(got[column], want[column], DIRECT_F64, &what);
                }
            }
        }
        cases += 1;
    }
    assert_eq!(cases, 1);
}

// --- Input checks through the constructor --------------------------------------

/// Coordinates, values, and slopes with one fault, and its error.
type Fault<'a> = (&'a [f64], &'a [f64], &'a [f64], E);

#[test]
fn the_constructor_reports_each_fault() {
    assert!(CubicHermiteSpline::new(&[0.0, 1.0], &[1.0, 2.0], &[0.0, 0.0]).is_ok());
    let ok: &[f64] = &[1.0, 2.0, 3.0];
    let cases: [Fault; 9] = [
        (
            &[0.0, 1.0, 2.0],
            &[1.0, 2.0],
            ok,
            E::LengthMismatch {
                expected: 3,
                actual: 2,
            },
        ),
        (
            &[0.0, 1.0, 2.0],
            ok,
            &[1.0, 2.0, 3.0, 4.0],
            E::LengthMismatch {
                expected: 3,
                actual: 4,
            },
        ),
        (
            &[0.0],
            &[1.0],
            &[1.0],
            E::InsufficientPoints {
                required: 2,
                actual: 1,
            },
        ),
        (&[0.0, f64::NAN, 2.0], ok, ok, E::NonFinite { index: 1 }),
        (
            &[0.0, 1.0, 2.0],
            &[1.0, 2.0, f64::INFINITY],
            ok,
            E::NonFinite { index: 2 },
        ),
        (
            &[0.0, 1.0, 2.0],
            ok,
            &[f64::NAN, 2.0, 3.0],
            E::NonFinite { index: 0 },
        ),
        (
            &[0.0, 1.0, 2.0],
            ok,
            &[1.0, f64::NEG_INFINITY, 3.0],
            E::NonFinite { index: 1 },
        ),
        (&[0.0, 2.0, 1.0], ok, ok, E::NotSorted { index: 2 }),
        (
            &[0.0, 1.0, 1.0],
            ok,
            ok,
            E::DuplicatePoints {
                first: 1,
                second: 2,
            },
        ),
    ];
    for (x, y, dydx, want) in cases {
        assert_eq!(CubicHermiteSpline::new(x, y, dydx).err(), Some(want));
    }
    // Finite values whose difference is not: the secant of interval 1. The
    // coefficients of interval 0 are finite (the largest is three times its
    // secant), and those of interval 2 are not, but it comes later.
    let y = [0.0, -5e307, 1.5e308, 0.0];
    let err = CubicHermiteSpline::new(&[0.0, 1.0, 2.0, 3.0], &y, &[0.0; 4]);
    assert_eq!(err.err(), Some(E::NonFinite { index: 1 }));
}

#[test]
fn the_slopes_are_checked_after_the_values_within_each_step() {
    let new = |x: &[f64], y: &[f64], dydx: &[f64]| CubicHermiteSpline::new(x, y, dydx).err();
    // Both lengths are wrong: that of the values is reported.
    assert_eq!(
        new(&[0.0, 1.0], &[1.0], &[1.0, 2.0, 3.0]),
        Some(E::LengthMismatch {
            expected: 2,
            actual: 1
        })
    );
    // The length of the slopes comes before the count of points.
    assert_eq!(
        new(&[0.0], &[1.0], &[]),
        Some(E::LengthMismatch {
            expected: 1,
            actual: 0
        })
    );
    // A value and a slope are NaN, the slope's first by index: the value.
    assert_eq!(
        new(
            &[0.0, 1.0, 2.0],
            &[1.0, 2.0, f64::NAN],
            &[f64::NAN, 0.0, 0.0]
        ),
        Some(E::NonFinite { index: 2 })
    );
    // A slope that is NaN comes before the order of the coordinates.
    assert_eq!(
        new(&[1.0, 0.0, 2.0], &[1.0, 2.0, 3.0], &[0.0, 0.0, f64::NAN]),
        Some(E::NonFinite { index: 2 })
    );
    // The order of the coordinates comes before an overflow.
    assert_eq!(
        new(&[1.0, 0.0], &[-1e308, 1e308], &[0.0, 0.0]),
        Some(E::NotSorted { index: 1 })
    );
    // An array value: one column of one slope.
    let err = CubicHermiteSpline::new(
        &[0.0, 1.0],
        &[[1.0, 2.0]; 2],
        &[[0.0, 0.0], [0.0, f64::INFINITY]],
    );
    assert_eq!(err.err(), Some(E::NonFinite { index: 1 }));
}

#[test]
fn an_overflowing_coefficient_names_the_first_such_interval() {
    // Finite slopes of 1e300 on an interval 1e-10 wide: the cubic coefficient
    // is their sum over the width squared. Intervals 1 and 3 overflow.
    let x = [0.0, 1.0, 1.0 + 1e-10, 2.0, 2.0 + 1e-10];
    let y = [0.0; 5];
    let dydx = [0.0, 1e300, 1e300, 1e300, 1e300];
    assert_eq!(
        CubicHermiteSpline::new(&x, &y, &dydx).err(),
        Some(E::NonFinite { index: 1 })
    );
}

// --- What is reproduced ----------------------------------------------------------

#[test]
fn two_points_give_one_cubic() {
    // 3 x^2 - 2 x^3 on [0, 1], written from 0.
    let step = CubicHermiteSpline::new(&[0.0, 1.0], &[0.0, 1.0], &[0.0, 0.0]).expect("valid");
    let c = step.as_ppoly().coefficients();
    assert_eq!(
        [c[0], c[1], c[2], c[3]].map(f64::to_bits),
        [-2.0_f64, 3.0, 0.0, 0.0].map(f64::to_bits)
    );
    assert_eq!(bits(step.eval(&0.5)), bits(Some(0.5)));
    assert_eq!(bits(step.eval(&2.0)), bits(Some(-4.0)), "extended");

    // Widths other than 1: secant 3, t = (1 + 2 - 6) / 2, and the piece is
    // -0.75 dx^3 + 2.5 dx^2 + dx + 2.
    let wide = CubicHermiteSpline::new(&[1.0, 3.0], &[2.0, 8.0], &[1.0, 2.0]).expect("valid");
    let c = wide.as_ppoly().coefficients();
    assert_eq!(
        [c[0], c[1], c[2], c[3]].map(f64::to_bits),
        [-0.75_f64, 2.5, 1.0, 2.0].map(f64::to_bits)
    );
    assert_eq!(bits(wide.eval(&3.0)), bits(Some(8.0)));
    assert_eq!(
        bits(wide.as_ppoly().eval_derivative(3.0, 1)),
        bits(Some(2.0))
    );
}

#[test]
fn the_slope_of_each_point_is_used_at_that_point() {
    // Slopes that are all different: a slope taken from a neighbor shows.
    let x = [0.0, 1.0, 3.0, 3.5, 6.0];
    let y = [1.0, -2.0, 0.5, 4.0, 4.0];
    let dydx = [0.25, -3.0, 7.0, 0.0, -1.5];
    let spline = CubicHermiteSpline::new(&x, &y, &dydx).expect("valid");
    let poly = spline.as_ppoly();
    for i in 0..5 {
        let value = spline.eval(&x[i]).expect("inside");
        let slope = poly.eval_derivative(x[i], 1).expect("inside");
        assert_close(value, y[i], DIRECT_F64, &format!("value at point {i}"));
        assert_close(slope, dydx[i], DIRECT_F64, &format!("slope at point {i}"));
    }
    // From the left as well: the piece before a point ends with its value and
    // its slope.
    let c = poly.coefficients();
    for i in 0..4 {
        let at_right = piece(c, i, x[i + 1] - x[i]);
        assert_close(
            at_right[0],
            y[i + 1],
            DIRECT_F64,
            &format!("piece {i} value"),
        );
        assert_close(
            at_right[1],
            dydx[i + 1],
            DIRECT_F64,
            &format!("piece {i} slope"),
        );
    }
}

#[test]
fn each_column_has_its_own_slopes() {
    let x = [0.0, 1.0, 2.5, 3.0, 4.0];
    let (first, second) = ([0.0, 1.0, 1.5, 4.0, 4.5], [2.0, -1.0, 3.0, 3.0, 0.0]);
    let (slopes_1, slopes_2) = ([1.0, 0.0, -2.0, 0.5, 3.0], [-4.0, 2.0, 0.0, 1.0, -0.5]);
    let zip = |a: [f64; 5], b: [f64; 5]| -> Vec<[f64; 2]> {
        a.iter().zip(b).map(|(a, b)| [*a, b]).collect()
    };
    let pair =
        CubicHermiteSpline::new(&x, &zip(first, second), &zip(slopes_1, slopes_2)).expect("valid");
    let one = CubicHermiteSpline::new(&x, &first, &slopes_1).expect("valid");
    let two = CubicHermiteSpline::new(&x, &second, &slopes_2).expect("valid");
    for q in [-1.0, 0.0, 0.3, 1.0, 2.0, 2.75, 3.5, 4.0, 5.0] {
        let got = pair.eval(&q).expect("extends");
        assert_eq!(Some(got[0].to_bits()), bits(one.eval(&q)), "x = {q}");
        assert_eq!(Some(got[1].to_bits()), bits(two.eval(&q)), "x = {q}");
    }
}

#[test]
fn periodic_data_is_not_asked_for() {
    // The ends differ in value and in slope. The spline is built, and repeats
    // as it is: one period on, the value at the left end is the one given
    // there, not the one given at the right end.
    let x = [0.0, 1.0, 2.0];
    let spline = CubicHermiteSpline::new(&x, &[0.0, 1.0, 5.0], &[1.0, 0.0, -2.0])
        .expect("no check of periodicity")
        .with_extrapolate(Extrapolate::Periodic);
    assert_eq!(spline.extrapolate(), Extrapolate::Periodic);
    assert_eq!(
        bits(spline.eval(&2.0)),
        bits(Some(0.0)),
        "wraps to the left"
    );
    assert_eq!(bits(spline.eval(&4.0)), bits(Some(0.0)));
    assert_eq!(bits(spline.eval(&2.5)), bits(spline.eval(&0.5)));
    // Just before the end of a period the value is near 5: a jump.
    let before = spline.eval(&1.999_999).expect("periodic");
    assert!((before - 5.0).abs() < 1e-5, "{before}");
}

#[test]
fn the_slopes_of_a_cubic_spline_give_its_coefficients_back() {
    let x = [0.0, 0.4, 1.0, 1.5, 2.75, 3.0, 4.25, 5.0];
    let y = [1.0, -0.5, 2.0, 2.25, -1.0, 0.0, 3.0, 1.5];
    let n = x.len();
    for bc in [BoundaryCondition::NotAKnot, BoundaryCondition::natural()] {
        let spline = CubicSpline::new(&x, &y, bc).expect("valid");
        let source = spline.as_ppoly().coefficients();
        // The slope at each point: the linear coefficient of the piece that
        // starts there, and for the last point the slope where the last piece
        // ends.
        let mut dydx: Vec<f64> = (0..n - 1).map(|i| source[4 * i + 2]).collect();
        dydx.push(piece(source, n - 2, x[n - 1] - x[n - 2])[1]);
        let hermite = CubicHermiteSpline::new(&x, &y, &dydx).expect("valid");
        let rebuilt = hermite.as_ppoly().coefficients();
        for (i, (got, want)) in rebuilt.iter().zip(source).enumerate() {
            let width = x[i / 4 + 1] - x[i / 4];
            let reach = width.powi(3 - i32::try_from(i % 4).expect("0 to 3"));
            let what = format!("{bc:?} coefficient {i}");
            assert_close(got * reach, want * reach, DIRECT_F64, &what);
        }
    }
}

#[test]
fn modes_f32_and_thread_safety() {
    let spline: CubicHermiteSpline<f64> =
        CubicHermiteSpline::new(&[0.0, 1.0, 2.0], &[0.0, 1.0, 4.0], &[0.0, 2.0, 4.0])
            .expect("valid");
    assert_eq!(spline.extrapolate(), Extrapolate::Extend);
    assert_eq!(spline.domain().0.to_bits(), 0.0_f64.to_bits());
    assert_eq!(spline.domain().1.to_bits(), 2.0_f64.to_bits());
    // x^2 from its values and slopes, extended.
    assert_eq!(bits(spline.eval(&3.0)), bits(Some(9.0)));
    let forbid = spline.clone().with_extrapolate(Extrapolate::Forbid);
    assert_eq!(forbid.extrapolate(), Extrapolate::Forbid);
    assert_eq!(forbid.as_ppoly().extrapolate(), Extrapolate::Forbid);
    assert_eq!(forbid.eval(&3.0), None);
    assert_eq!(forbid.eval(&f64::NAN), None);
    let clamp = spline.clone().with_extrapolate(Extrapolate::Clamp);
    assert_eq!(bits(clamp.eval(&3.0)), bits(Some(4.0)));
    let periodic = spline.with_extrapolate(Extrapolate::Periodic);
    assert_eq!(bits(periodic.eval(&2.5)), bits(Some(0.25)));

    let single: CubicHermiteSpline<f32> =
        CubicHermiteSpline::new(&[0.0, 1.0], &[0.0, 1.0], &[0.0, 0.0]).expect("valid");
    assert_eq!(single.eval(&0.5).map(f32::to_bits), Some(0.5_f32.to_bits()));
    assert_eq!(
        single.eval(&0.25).map(f32::to_bits),
        Some(0.15625_f32.to_bits())
    );

    fn assert_send_sync<T: Send + Sync>() {}
    assert_send_sync::<CubicHermiteSpline<f64>>();
    assert_send_sync::<CubicHermiteSpline<f32, [f32; 3]>>();
}

// --- Properties ----------------------------------------------------------------

/// Coordinates that increase strictly, at least two.
fn coordinates() -> impl Strategy<Value = Vec<f64>> {
    (-5.0_f64..5.0, prop::collection::vec(0.05_f64..3.0, 1..12)).prop_map(|(start, widths)| {
        let mut x = vec![start];
        for w in widths {
            x.push(x[x.len() - 1] + w);
        }
        x
    })
}

/// Values, or slopes, for those coordinates.
fn values() -> impl Strategy<Value = Vec<f64>> {
    prop::collection::vec(prop_oneof![-10.0_f64..10.0, Just(0.0)], 13)
}

proptest! {
    #[test]
    fn the_values_and_the_slopes_given_are_reproduced(
        x in coordinates(),
        values in values(),
        slopes in values(),
    ) {
        let n = x.len();
        let (y, dydx) = (&values[..n], &slopes[..n]);
        let spline = CubicHermiteSpline::new(&x, y, dydx).expect("valid data");
        let c = spline.as_ppoly().coefficients();
        let value_scale = y.iter().fold(1.0_f64, |m, v| m.max(v.abs()));
        let slope_scale = dydx.iter().fold(1.0_f64, |m, v| m.max(v.abs()));
        for i in 0..n - 1 {
            // At the left end of a piece: the coefficients themselves.
            prop_assert_eq!(c[4 * i + 3].to_bits(), y[i].to_bits());
            prop_assert_eq!(c[4 * i + 2].to_bits(), dydx[i].to_bits());
            // At its right end, to rounding.
            let h = x[i + 1] - x[i];
            let at_right = piece(c, i, h);
            let value_tol = 1e-12 * (value_scale + slope_scale * h);
            prop_assert!((at_right[0] - y[i + 1]).abs() <= value_tol, "piece {} ends at {} for {}", i, at_right[0], y[i + 1]);
            let slope_tol = 1e-12 * (slope_scale + value_scale / h);
            prop_assert!((at_right[1] - dydx[i + 1]).abs() <= slope_tol, "piece {} ends with slope {} for {}", i, at_right[1], dydx[i + 1]);
        }
    }

    #[test]
    fn a_cubic_polynomial_is_reproduced(
        x in coordinates(),
        a in prop::array::uniform4(-3.0_f64..3.0),
        ts in prop::collection::vec(-0.2_f64..1.2, 1..30),
    ) {
        let f = |x: f64| ((a[0] * x + a[1]) * x + a[2]) * x + a[3];
        let df = |x: f64| (3.0 * a[0] * x + 2.0 * a[1]) * x + a[2];
        let y: Vec<f64> = x.iter().map(|x| f(*x)).collect();
        let dydx: Vec<f64> = x.iter().map(|x| df(*x)).collect();
        let spline = CubicHermiteSpline::new(&x, &y, &dydx).expect("valid data");
        let (lo, hi) = spline.domain();
        // The size of the terms of the polynomial over the queries, which reach
        // past both ends.
        let reach = lo.abs().max(hi.abs()) + 0.2 * (hi - lo);
        let scale = 1.0 + 3.0 * (reach + reach.powi(2) + reach.powi(3));
        for t in &ts {
            let q = lo + t * (hi - lo);
            let got = spline.eval(&q).expect("extends");
            prop_assert!((got - f(q)).abs() <= 1e-10 * scale, "x = {}: {} for {}", q, got, f(q));
        }
    }

    #[test]
    fn columns_are_interpolated_on_their_own(
        x in coordinates(),
        first in values(),
        second in values(),
        slopes_1 in values(),
        slopes_2 in values(),
        magnitude in -12_i32..=12,
        qs in prop::collection::vec(-12.0_f64..40.0, 1..40),
    ) {
        let n = x.len();
        let factor = 10.0_f64.powi(magnitude);
        let second: Vec<f64> = second[..n].iter().map(|v| v * factor).collect();
        let slopes_2: Vec<f64> = slopes_2[..n].iter().map(|v| v * factor).collect();
        let y: Vec<[f64; 2]> = first[..n].iter().zip(&second).map(|(a, b)| [*a, *b]).collect();
        let dydx: Vec<[f64; 2]> = slopes_1[..n].iter().zip(&slopes_2).map(|(a, b)| [*a, *b]).collect();
        let pair = CubicHermiteSpline::new(&x, &y, &dydx).expect("valid data");
        let one = CubicHermiteSpline::new(&x, &first[..n], &slopes_1[..n]).expect("valid data");
        let two = CubicHermiteSpline::new(&x, &second, &slopes_2).expect("valid data");
        for q in &qs {
            let got = pair.eval(q).expect("extends");
            prop_assert_eq!(Some(got[0].to_bits()), bits(one.eval(q)), "x = {}", q);
            prop_assert_eq!(Some(got[1].to_bits()), bits(two.eval(q)), "x = {}", q);
        }
    }

    #[test]
    fn evaluation_is_that_of_the_piecewise_polynomial(
        x in coordinates(),
        values in values(),
        slopes in values(),
        qs in prop::collection::vec(-12.0_f64..40.0, 1..100),
    ) {
        let (y, dydx) = (&values[..x.len()], &slopes[..x.len()]);
        for mode in [Extrapolate::Forbid, Extrapolate::Clamp, Extrapolate::Extend, Extrapolate::Periodic] {
            let spline = CubicHermiteSpline::new(&x, y, dydx).expect("valid data").with_extrapolate(mode);
            let poly = spline.as_ppoly();
            prop_assert_eq!(poly.extrapolate(), mode);
            let mut hint = Hint::new();
            for q in &qs {
                prop_assert_eq!(bits(spline.eval(q)), bits(poly.eval(q)), "{:?} x = {}", mode, q);
                prop_assert_eq!(bits(spline.eval_hinted(*q, &mut hint)), bits(poly.eval(q)));
            }
            let (mut out, mut want) = (vec![0.0; qs.len()], vec![0.0; qs.len()]);
            prop_assert_eq!(spline.eval_batch_into(&qs, &mut out, -7.0), Ok(()));
            prop_assert_eq!(poly.eval_batch_into(&qs, &mut want, -7.0), Ok(()));
            for (got, want) in out.iter().zip(&want) {
                prop_assert_eq!(got.to_bits(), want.to_bits());
            }
        }
    }
}
