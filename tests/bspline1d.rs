//! `BSpline1D`: the golden data, the input checks and their order through the
//! constructor, the knots of each case, the other interpolators it agrees with,
//! and the properties.

mod common;

use common::{assert_close, load, reals, Real, Tol, DIRECT_F64};
use interpx::{
    BSpline1D, BoundaryCondition as Bc, CubicSpline, EndCondition as End, Extrapolate, Hint,
    InterpolateError as E, Interpolator, Linear1D, ParamReason, BSPLINE_MAX_DEGREE,
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

/// One end of SciPy's `bc_type`: a name, or `[order, value]` with one value per
/// column.
fn end_of<const M: usize>(v: &Json) -> Option<End<[f64; M]>> {
    if v == "not-a-knot" {
        return Some(End::NotAKnot);
    }
    let value: [f64; M] = match &v[1] {
        Json::Array(_) => nums(&v[1]).try_into().ok()?,
        number => [number.as_f64()?; M],
    };
    match v[0].as_u64()? {
        1 => Some(End::FirstDerivative(value)),
        2 => Some(End::SecondDerivative(value)),
        _ => None,
    }
}

/// SciPy's `bc_type` for values of `M` columns. No condition is not-a-knot.
fn bc_of<const M: usize>(v: &Json) -> Option<Bc<[f64; M]>> {
    if v.is_null() {
        return Some(Bc::NotAKnot);
    }
    match v.as_str() {
        Some("not-a-knot") => Some(Bc::NotAKnot),
        Some("periodic") => Some(Bc::Periodic),
        Some("natural") => Some(Bc::natural()),
        Some("clamped") => Some(Bc::clamped()),
        Some(_) => None,
        None => Some(Bc::Ends(end_of(&v[0])?, end_of(&v[1])?)),
    }
}

/// The same condition for scalar values.
fn scalar(bc: Bc<[f64; 1]>) -> Bc<f64> {
    let end = |e: End<[f64; 1]>| match e {
        End::NotAKnot => End::NotAKnot,
        End::FirstDerivative([v]) => End::FirstDerivative(v),
        End::SecondDerivative([v]) => End::SecondDerivative(v),
    };
    match bc {
        Bc::NotAKnot => Bc::NotAKnot,
        Bc::Periodic => Bc::Periodic,
        Bc::Ends(a, b) => Bc::Ends(end(a), end(b)),
    }
}

// Marked as test code so that `expect` is allowed here as in a `#[test]` function
// (clippy.toml allows it in tests only). The same for the two helpers below.
#[cfg(test)]
fn degree_of(case: &Json) -> usize {
    usize::try_from(case["options"]["k"].as_u64().expect("k")).expect("a small degree")
}

/// The spline of a golden case with scalar values, in the mode of the case.
#[cfg(test)]
fn built(case: &Json) -> BSpline1D<f64> {
    let name = case["name"].as_str().expect("name");
    let bc = scalar(bc_of::<1>(&case["options"]["bc_type"]).expect("a boundary condition"));
    let spline = BSpline1D::new(&nums(&case["x"]), &nums(&case["y"]), degree_of(case), bc)
        .expect("valid data");
    let default = if bc == Bc::Periodic {
        Extrapolate::Periodic
    } else {
        Extrapolate::Extend
    };
    assert_eq!(spline.extrapolate(), default, "{name}: default mode");
    if case["options"]["extrapolate"] == Json::Bool(false) {
        spline.with_extrapolate(Extrapolate::Forbid)
    } else {
        spline
    }
}

/// The tolerance for a second derivative of a spline with a narrow piece at an
/// end (see [`narrow_end`]), inside the data.
const SECOND_DERIVATIVE_NARROW_END: Tol = Tol {
    rtol: 1e-11,
    atol: 1e-11,
};

/// Whether the first or the last piece of the spline is under a hundredth of
/// the data: `random_20` ends with two points 0.018 apart in a span of 9.2,
/// and a spline with a derivative given at that end has a knot at the last
/// point but one.
///
/// On such a piece the basis functions are steep. One rounding of the control
/// points (which agree with SciPy's to 2.4e-16) moves the second derivative
/// by 6e-12 inside the data, and the value by 2e-11 where the piece is
/// continued 50 widths past the end: more than the tolerance for a direct
/// method. So such a case is compared inside the data only, and its second
/// derivative with a tolerance of its own. (Measured on
/// `bspline/k3/first_and_second/random_20`, the one case that needs it;
/// `tests/cubic_spline.rs` treats its `uneven` cases the same way.)
fn narrow_end(spline: &BSpline1D<f64>) -> bool {
    let b = spline.as_bspline();
    let (knots, count) = (b.knots(), b.control_points().len());
    let (lo, hi) = spline.domain();
    let first = knots[b.degree() + 1] - lo;
    let last = hi - knots[count - 1];
    first.min(last) < (hi - lo) / 100.0
}

/// Compares the spline with the columns of values the case has. Returns
/// whether the case was compared as one with a narrow end.
#[cfg(test)]
fn check_values(case: &Json, spline: &BSpline1D<f64>) -> bool {
    let name = case["name"].as_str().expect("name");
    let queries = nums(&case["queries"]);
    let narrow = narrow_end(spline) && spline.extrapolate() != Extrapolate::Periodic;
    let (lo, hi) = spline.domain();
    for (field, nu) in [("expected", 0), ("derivative_1", 1), ("derivative_2", 2)] {
        let Some(want) = case.get(field) else {
            assert_eq!(field, "derivative_2", "{name}: only this is ever left out");
            continue;
        };
        let tol = if narrow && nu == 2 {
            SECOND_DERIVATIVE_NARROW_END
        } else {
            DIRECT_F64
        };
        for (k, (q, want)) in queries.iter().zip(nums(want)).enumerate() {
            if narrow && !(lo..=hi).contains(q) {
                continue;
            }
            let what = format!("{name} {field} query {k} (x = {q})");
            match spline.as_bspline().eval_derivative(*q, nu) {
                None => assert!(want.is_nan(), "{what}: None, SciPy has {want}"),
                Some(got) => assert_close(got, want, tol, &what),
            }
            if nu == 0 {
                assert_eq!(
                    bits(spline.eval(q)),
                    bits(spline.as_bspline().eval(q)),
                    "{what}"
                );
            }
        }
    }
    narrow
}

// --- Golden data -------------------------------------------------------------

#[test]
fn scipy_interpolating_splines_are_rebuilt_from_their_data() {
    let file = load::<Json>("bspline1d").expect("golden file bspline1d");
    let mut degrees = std::collections::BTreeMap::new();
    let (mut cases, mut periodic, mut with_derivatives, mut forbid) = (0, 0, 0, 0);
    let mut narrow = 0;
    for case in &file.cases {
        let name = case["name"].as_str().expect("name");
        if case["reference"] != "scipy.interpolate.make_interp_spline"
            || case.get("knots").is_none()
        {
            continue;
        }
        let spline = built(case);
        let b = spline.as_bspline();
        assert_eq!(b.degree(), degree_of(case), "{name}");
        assert_eq!(spline.degree(), degree_of(case), "{name}");

        let knots = nums(&case["knots"]);
        assert_eq!(b.knots().len(), knots.len(), "{name}: knots");
        for (i, (got, want)) in b.knots().iter().zip(&knots).enumerate() {
            assert_close(*got, *want, DIRECT_F64, &format!("{name} knot {i}"));
        }
        let control_points = nums(&case["control_points"]);
        assert_eq!(b.control_points().len(), control_points.len(), "{name}");
        for (i, (got, want)) in b.control_points().iter().zip(&control_points).enumerate() {
            let what = format!("{name} control point {i}");
            assert_close(*got, *want, DIRECT_F64, &what);
        }

        narrow += usize::from(check_values(case, &spline));
        if let Some(integrals) = case.get("integrals") {
            for triple in rows(integrals) {
                let got = b.integrate(triple[0], triple[1]).expect("inside the data");
                assert_close(got, triple[2], DIRECT_F64, &format!("{name} integral"));
            }
        }

        cases += 1;
        *degrees.entry(degree_of(case)).or_insert(0) += 1;
        periodic += usize::from(spline.extrapolate() == Extrapolate::Periodic);
        with_derivatives += usize::from(case["options"]["bc_type"].is_array());
        forbid += usize::from(spline.extrapolate() == Extrapolate::Forbid);
        if matches!(
            case["options"]["bc_type"].as_str(),
            Some("natural" | "clamped")
        ) {
            with_derivatives += 1;
        }
    }
    // 33 from before `BSpline1D` (degrees 1 to 5 on five data sets, the cubic
    // conditions, one without extrapolation, one periodic) and 28 chosen for
    // it.
    assert_eq!((cases, periodic, with_derivatives, forbid), (61, 14, 8, 1));
    // The cases compared with the rule for a narrow end, all on `random_20`:
    // degree 1 (a knot at every point) and the three cubic conditions with a
    // derivative at the right end.
    assert_eq!(narrow, 4);
    let degrees: Vec<(usize, i32)> = degrees.into_iter().collect();
    assert_eq!(
        degrees,
        [(1, 5), (2, 12), (3, 22), (4, 8), (5, 10), (6, 2), (7, 2)]
    );
}

#[test]
fn scipy_splines_of_two_columns_are_rebuilt_column_by_column() {
    let file = load::<Json>("bspline1d").expect("golden file bspline1d");
    let mut cases = 0;
    for case in &file.cases {
        let name = case["name"].as_str().expect("name");
        if !name.ends_with("/columns") {
            continue;
        }
        let bc = bc_of::<2>(&case["options"]["bc_type"]).expect("a boundary condition");
        let spline = BSpline1D::new(&nums(&case["x"]), &pairs(&case["y"]), degree_of(case), bc)
            .expect("valid data");
        let queries = nums(&case["queries"]);
        for (field, nu) in [("expected", 0), ("derivative_1", 1), ("derivative_2", 2)] {
            for (k, (q, want)) in queries.iter().zip(rows(&case[field])).enumerate() {
                let got = spline
                    .as_bspline()
                    .eval_derivative(*q, nu)
                    .expect("the case extrapolates");
                for column in 0..2 {
                    let what = format!("{name} {field} query {k} column {column}");
                    assert_close(got[column], want[column], DIRECT_F64, &what);
                }
            }
        }
        cases += 1;
    }
    assert_eq!(cases, 3);
}

#[test]
fn a_cubic_with_not_a_knot_at_one_end_is_scipys_cubic_spline() {
    // `make_interp_spline` does not take these conditions; `CubicSpline` does,
    // and the curve is the same.
    let file = load::<Json>("bspline1d").expect("golden file bspline1d");
    let mut cases = 0;
    for case in &file.cases {
        if case["reference"] != "scipy.interpolate.CubicSpline" {
            continue;
        }
        let spline = built(case);
        assert_eq!(spline.degree(), 3);
        let _ = check_values(case, &spline);
        cases += 1;
    }
    assert_eq!(cases, 3);
}

// --- Input checks through the constructor --------------------------------------

/// Coordinates, values, a degree, a boundary condition, and the error they give.
type Fault<'a> = (&'a [f64], &'a [f64], usize, Bc<f64>, E);

fn invalid(name: &'static str, reason: ParamReason) -> E {
    E::InvalidParameter { name, reason }
}

#[test]
fn the_constructor_reports_each_fault() {
    let x5: &[f64] = &[0.0, 1.0, 2.0, 3.0, 4.0];
    let y5: &[f64] = &[1.0, 2.0, 0.0, 3.0, 1.0];
    let first = Bc::Ends(End::FirstDerivative(0.5), End::NotAKnot);
    let cases: [Fault; 14] = [
        (
            x5,
            &[1.0, 2.0],
            3,
            Bc::NotAKnot,
            E::LengthMismatch {
                expected: 5,
                actual: 2,
            },
        ),
        (
            x5,
            y5,
            0,
            Bc::NotAKnot,
            invalid("degree", ParamReason::TooSmall { min: 1 }),
        ),
        (
            x5,
            y5,
            16,
            Bc::NotAKnot,
            invalid("degree", ParamReason::TooLarge { max: 15 }),
        ),
        (
            x5,
            y5,
            2,
            Bc::natural(),
            invalid("bc", ParamReason::Unsupported),
        ),
        (x5, y5, 4, first, invalid("bc", ParamReason::Unsupported)),
        (x5, y5, 1, first, invalid("bc", ParamReason::Unsupported)),
        (
            x5,
            y5,
            5,
            Bc::NotAKnot,
            E::InsufficientPoints {
                required: 6,
                actual: 5,
            },
        ),
        (
            &[0.0, f64::NAN, 2.0, 3.0, 4.0],
            y5,
            3,
            Bc::NotAKnot,
            E::NonFinite { index: 1 },
        ),
        (
            x5,
            &[1.0, 2.0, 0.0, f64::INFINITY, 1.0],
            3,
            Bc::NotAKnot,
            E::NonFinite { index: 3 },
        ),
        (
            &[0.0, 1.0, 3.0, 2.0, 4.0],
            y5,
            3,
            Bc::NotAKnot,
            E::NotSorted { index: 3 },
        ),
        (
            &[0.0, 1.0, 1.0, 2.0, 4.0],
            y5,
            3,
            Bc::NotAKnot,
            E::DuplicatePoints {
                first: 1,
                second: 2,
            },
        ),
        (
            x5,
            y5,
            3,
            Bc::Ends(End::SecondDerivative(f64::NAN), End::NotAKnot),
            invalid("bc", ParamReason::NonFinite),
        ),
        (
            x5,
            y5,
            3,
            Bc::Ends(End::NotAKnot, End::FirstDerivative(f64::INFINITY)),
            invalid("bc", ParamReason::NonFinite),
        ),
        (
            x5,
            &[1.0, 2.0, 0.0, 3.0, 2.0],
            4,
            Bc::Periodic,
            E::NotPeriodic,
        ),
    ];
    for (x, y, degree, bc, want) in cases {
        assert_eq!(
            BSpline1D::new(x, y, degree, bc).err(),
            Some(want),
            "degree {degree} {bc:?}"
        );
    }
    assert_eq!(BSPLINE_MAX_DEGREE, 15);
    assert!(BSpline1D::new(&[0.0; 16].map(f64::from), &[0.0; 16], 15, Bc::NotAKnot).is_err());
}

#[test]
fn the_checks_come_in_the_order_of_the_design() {
    let new = |x: &[f64], y: &[f64], degree: usize, bc: Bc<f64>| BSpline1D::new(x, y, degree, bc);
    let nan = Bc::Ends(End::FirstDerivative(f64::NAN), End::NotAKnot);
    // The length before the degree.
    assert_eq!(
        new(&[0.0, 1.0], &[1.0], 0, Bc::NotAKnot).err(),
        Some(E::LengthMismatch {
            expected: 2,
            actual: 1
        })
    );
    // The degree before the boundary condition: degree 0 with a derivative.
    assert_eq!(
        new(&[0.0, 1.0], &[1.0, 2.0], 0, Bc::natural()).err(),
        Some(invalid("degree", ParamReason::TooSmall { min: 1 }))
    );
    assert_eq!(
        new(&[0.0, 1.0], &[1.0, 2.0], 16, Bc::natural()).err(),
        Some(invalid("degree", ParamReason::TooLarge { max: 15 }))
    );
    // The boundary condition and the degree before the count of points: one
    // point, degree 2, a derivative.
    assert_eq!(
        new(&[0.0], &[1.0], 2, Bc::clamped()).err(),
        Some(invalid("bc", ParamReason::Unsupported))
    );
    // The count before a NaN in the data.
    assert_eq!(
        new(&[f64::NAN, 1.0, 2.0], &[1.0, 2.0, 3.0], 3, Bc::NotAKnot).err(),
        Some(E::InsufficientPoints {
            required: 4,
            actual: 3
        })
    );
    // The data before the value of the boundary condition: unsorted
    // coordinates and a derivative that is NaN.
    assert_eq!(
        new(&[0.0, 2.0, 1.0, 3.0], &[1.0; 4], 3, nan).err(),
        Some(E::NotSorted { index: 2 })
    );
    // The value of the boundary condition before anything is solved.
    assert_eq!(
        new(&[0.0, 1.0, 2.0, 3.0], &[1.0; 4], 3, nan).err(),
        Some(invalid("bc", ParamReason::NonFinite))
    );
    // Periodicity is checked for every degree, 1 included.
    assert_eq!(
        new(&[0.0, 1.0, 2.0], &[1.0, 2.0, 3.0], 1, Bc::Periodic).err(),
        Some(E::NotPeriodic)
    );
    // Its tolerance is that of `CubicSpline`: a rounding apart passes.
    let near = [0.0_f64.sin(), 1.0, (2.0 * std::f64::consts::PI).sin()];
    assert_ne!(near[0].to_bits(), near[2].to_bits());
    assert!(new(&[0.0, 1.0, 2.0], &near, 2, Bc::Periodic).is_ok());
}

#[test]
fn the_fewest_points_follow_the_degree_and_the_derivatives_given() {
    let x: Vec<f64> = (0..20).map(f64::from).collect();
    let y: Vec<f64> = x.iter().map(|v| (0.7 * v).sin()).collect();
    let fewest = |degree: usize, bc: Bc<f64>| {
        (1..=20)
            .find(|n| BSpline1D::new(&x[..*n], &y[..*n], degree, bc).is_ok())
            .expect("twenty points are enough")
    };
    // Not-a-knot: one polynomial of the degree, and two points at least.
    for degree in 1..=BSPLINE_MAX_DEGREE {
        assert_eq!(
            fewest(degree, Bc::NotAKnot),
            (degree + 1).max(2),
            "{degree}"
        );
        assert_eq!(
            fewest(degree, Bc::Ends(End::NotAKnot, End::NotAKnot)),
            (degree + 1).max(2)
        );
        let short = BSpline1D::new(&x[..degree], &y[..degree], degree, Bc::NotAKnot).err();
        let required = (degree + 1).max(2);
        assert_eq!(
            short,
            Some(E::InsufficientPoints {
                required,
                actual: degree
            })
        );
    }
    // Degree 3: each derivative stands for one point.
    let one = Bc::Ends(End::NotAKnot, End::SecondDerivative(0.0));
    assert_eq!(fewest(3, one), 3);
    assert_eq!(fewest(3, Bc::natural()), 2);
    assert_eq!(
        BSpline1D::new(&x[..2], &y[..2], 3, one).err(),
        Some(E::InsufficientPoints {
            required: 3,
            actual: 2
        })
    );
    assert_eq!(
        BSpline1D::new(&x[..1], &y[..1], 3, Bc::natural()).err(),
        Some(E::InsufficientPoints {
            required: 2,
            actual: 1
        })
    );
    // Periodic: two points for any degree.
    for degree in 1..=BSPLINE_MAX_DEGREE {
        let ends = [1.5, 1.5];
        let spline = BSpline1D::new(&x[..2], &ends, degree, Bc::Periodic).expect("two points");
        // Two equal values: a constant.
        for q in [0.0, 0.25, 0.5, 1.0, 7.3, -2.9] {
            let got = spline.eval(&q).expect("periodic");
            assert!((got - 1.5).abs() <= 1e-13, "degree {degree} at {q}: {got}");
        }
        assert_eq!(
            BSpline1D::new(&x[..1], &y[..1], degree, Bc::Periodic).err(),
            Some(E::InsufficientPoints {
                required: 2,
                actual: 1
            })
        );
    }
}

#[test]
fn control_points_that_overflow_are_an_error_and_not_a_panic() {
    // Values near the largest number, alternating in sign: the control points
    // of the spline through them are larger.
    let x: Vec<f64> = (0..8).map(f64::from).collect();
    let y: Vec<f64> = (0..8)
        .map(|i| if i % 2 == 0 { 1.7e308 } else { -1.7e308 })
        .collect();
    for degree in [2, 3, 5] {
        let err = BSpline1D::new(&x, &y, degree, Bc::NotAKnot).err();
        assert!(
            matches!(err, Some(E::NonFinite { .. })),
            "{degree}: {err:?}"
        );
    }
    let mut closed = y.clone();
    closed[7] = closed[0];
    let err = BSpline1D::new(&x, &closed, 3, Bc::Periodic).err();
    assert!(matches!(err, Some(E::NonFinite { .. })), "{err:?}");
}

#[test]
fn an_even_degree_refuses_points_with_no_number_between_them() {
    // Neighboring floating-point numbers from the third point on: the number
    // halfway between two of them is one of them.
    let eps = f64::EPSILON;
    let x = [
        0.0,
        0.5,
        1.0,
        1.0 + eps,
        1.0 + 2.0 * eps,
        2.0,
        3.0,
        4.0,
        5.0,
    ];
    let y = [0.0, 1.0, 0.0, 2.0, 1.0, 0.0, 1.0, 2.0, 0.0];
    let too_close = E::DuplicatePoints {
        first: 2,
        second: 3,
    };
    for degree in [2, 4, 6, 8] {
        let open = BSpline1D::new(&x, &y, degree, Bc::NotAKnot).err();
        assert_eq!(open, Some(too_close), "degree {degree}");
        let periodic = BSpline1D::new(&x, &y, degree, Bc::Periodic).err();
        assert_eq!(periodic, Some(too_close), "degree {degree}, periodic");
    }
    // An odd degree has its knots at the data points, and takes the data. So
    // does degree 1.
    for degree in [1, 3, 5, 7] {
        assert!(
            BSpline1D::new(&x, &y, degree, Bc::NotAKnot).is_ok(),
            "{degree}"
        );
        assert!(
            BSpline1D::new(&x, &y, degree, Bc::Periodic).is_ok(),
            "{degree}"
        );
    }
    // The pair is refused even where its midpoint is not a knot (the knots of
    // degree 8 on 9 points are the two ends only): one rule for every pair.
    // In the last interval a periodic spline has no knot of its own, and takes
    // the pair.
    let mut at_the_end = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 7.0];
    at_the_end[7] = 7.0 - 4.0 * eps;
    assert!(at_the_end[7] < at_the_end[8]);
    let values = [1.0, 0.0, 2.0, 1.0, 0.0, 1.0, 2.0, 0.0, 1.0];
    assert_eq!(
        BSpline1D::new(&at_the_end, &values, 2, Bc::NotAKnot).err(),
        Some(E::DuplicatePoints {
            first: 7,
            second: 8
        })
    );
    assert!(BSpline1D::new(&at_the_end, &values, 2, Bc::Periodic).is_ok());
    // The smallest numbers there are: halving loses the last bit.
    let tiny: Vec<f64> = (0..5).map(|i| f64::from(i) * 5e-324).collect();
    assert_eq!(
        BSpline1D::new(&tiny, &[0.0; 5], 2, Bc::NotAKnot).err(),
        Some(E::DuplicatePoints {
            first: 0,
            second: 1
        })
    );
}

#[test]
fn coordinates_near_the_largest_number_are_an_error_or_a_spline_and_not_a_panic() {
    // The midpoints of an even degree do not overflow: 1.5e308 lies between.
    let x = [1.0e308_f64, 1.2e308, 1.4e308, 1.6e308, 1.7e308];
    let y = [0.0_f64, 1.0, 0.0, 1.0, 0.0];
    let spline = BSpline1D::new(&x, &y, 2, Bc::NotAKnot).expect("no overflow");
    let at = spline.eval(&1.4e308).expect("inside");
    assert!(at.abs() < 1e-9, "{at}");
    // Coordinates so far apart that their distance overflows: the basis
    // functions cannot be evaluated, and the solver reports it.
    let wide = [-1.7e308, -1.0e308, 0.0, 1.0e308, 1.7e308];
    for degree in [2, 3, 4] {
        let err = BSpline1D::new(&wide, &y, degree, Bc::NotAKnot).err();
        assert!(
            matches!(err, Some(E::NonFinite { .. })),
            "{degree}: {err:?}"
        );
        let err = BSpline1D::new(&wide, &y, degree, Bc::Periodic).err();
        assert!(
            matches!(err, Some(E::NonFinite { .. })),
            "{degree}: {err:?}"
        );
    }
    // Two neighbors whose distance overflows: the knot an even periodic degree
    // puts between them is not finite. The later of the two is named.
    let apart = [-1.7e308, -1.6e308, 1.6e308, 1.7e308];
    assert_eq!(
        BSpline1D::new(&apart, &[0.0, 1.0, 2.0, 0.0], 2, Bc::Periodic).err(),
        Some(E::NonFinite { index: 2 })
    );
    // The knots of a periodic spline reach past the data, and overflow there
    // although the data does not: the first knot, then the last.
    let right = [1.0e308, 1.2e308, 1.4e308, 1.6e308, 1.7e308];
    assert_eq!(
        BSpline1D::new(&right, &y, 5, Bc::Periodic).err(),
        Some(E::NonFinite { index: 4 })
    );
    let left = right.map(|v: f64| -v);
    let mut ascending = left;
    ascending.reverse();
    assert_eq!(
        BSpline1D::new(&ascending, &y, 5, Bc::Periodic).err(),
        Some(E::NonFinite { index: 0 })
    );
}

// --- The knots and the curve of each case ----------------------------------------

#[test]
fn degree_one_is_the_polygon_and_solves_nothing() {
    let x = [0.0_f64, 1.0, 2.5, 4.0];
    let y = [1.0_f64, -1.0, 4.0, 0.5];
    let spline = BSpline1D::new(&x, &y, 1, Bc::NotAKnot).expect("valid");
    let b = spline.as_bspline();
    assert_eq!(
        b.knots().iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
        [0.0_f64, 0.0, 1.0, 2.5, 4.0, 4.0].map(f64::to_bits)
    );
    // The control points are the values themselves.
    for (c, v) in b.control_points().iter().zip(&y) {
        assert_eq!(c.to_bits(), v.to_bits());
    }
    let linear = Linear1D::new(&x, &y)
        .expect("valid")
        .with_extrapolate(Extrapolate::Extend);
    for q in [-1.0, 0.0, 0.3, 1.0, 1.75, 2.5, 3.9, 4.0, 5.0] {
        let (got, want) = (
            spline.eval(&q).expect("extends"),
            linear.eval(&q).expect("extends"),
        );
        assert_close(got, want, DIRECT_F64, &format!("x = {q}"));
    }
}

#[test]
fn a_periodic_spline_of_degree_one_is_periodic() {
    // SciPy leaves it extrapolating; here every degree repeats.
    let x = [0.0, 1.0, 2.0, 4.0];
    let y = [1.0, 3.0, -1.0, 1.0];
    let spline = BSpline1D::new(&x, &y, 1, Bc::Periodic).expect("valid");
    assert_eq!(spline.extrapolate(), Extrapolate::Periodic);
    assert_eq!(bits(spline.eval(&4.5)), bits(spline.eval(&0.5)));
    assert_eq!(bits(spline.eval(&-3.0)), bits(spline.eval(&1.0)));
    assert_eq!(bits(spline.eval(&0.5)), bits(Some(2.0)));
}

#[test]
fn the_not_a_knot_knots_leave_out_the_points_next_to_the_ends() {
    let x: Vec<f64> = (0..9).map(|i| f64::from(i * i)).collect();
    let y = vec![0.0; 9];
    let inner = |degree: usize| {
        let spline = BSpline1D::new(&x, &y, degree, Bc::NotAKnot).expect("valid");
        let knots = spline.as_bspline().knots().to_vec();
        // Both ends with multiplicity `degree + 1`.
        assert!(knots[..=degree]
            .iter()
            .all(|t| t.to_bits() == x[0].to_bits()));
        let tail = &knots[knots.len() - degree - 1..];
        assert!(tail.iter().all(|t| t.to_bits() == x[8].to_bits()));
        knots[degree + 1..knots.len() - degree - 1].to_vec()
    };
    // Odd degrees: the data points but (degree + 1) / 2 at each end.
    assert_eq!(inner(3), x[2..7]);
    assert_eq!(inner(5), x[3..6]);
    assert_eq!(inner(7), x[4..5]);
    // Even degrees: the midpoints but degree / 2 at each end.
    let mids: Vec<f64> = x.windows(2).map(|w| (w[0] + w[1]) / 2.0).collect();
    assert_eq!(inner(2), mids[1..7]);
    assert_eq!(inner(4), mids[2..6]);
    assert_eq!(inner(8), Vec::<f64>::new());
    // Degree 3 with a derivative at an end keeps the point next to that end.
    let knots = |bc: Bc<f64>| {
        let spline = BSpline1D::new(&x, &y, 3, bc).expect("valid");
        let knots = spline.as_bspline().knots().to_vec();
        knots[4..knots.len() - 4].to_vec()
    };
    assert_eq!(knots(Bc::natural()), x[1..8]);
    assert_eq!(
        knots(Bc::Ends(End::FirstDerivative(1.0), End::NotAKnot)),
        x[1..7]
    );
    assert_eq!(
        knots(Bc::Ends(End::NotAKnot, End::SecondDerivative(1.0))),
        x[2..8]
    );
}

#[test]
fn the_periodic_knots_repeat_the_spacings_of_the_data() {
    let x = [0.0_f64, 1.0, 3.0, 3.5, 6.0];
    let y = [1.0_f64, 2.0, 0.0, -1.0, 1.0];
    // Odd: the data points, then the spacings 1, 2, 0.5, 2.5 again on each side.
    let odd = BSpline1D::new(&x, &y, 3, Bc::Periodic).expect("valid");
    assert_eq!(
        odd.as_bspline().knots(),
        [-5.0, -3.0, -2.5, 0.0, 1.0, 3.0, 3.5, 6.0, 7.0, 9.0, 9.5]
    );
    // Even: the inner points move half an interval to the left (0.5, 2, 3.25),
    // and the spacings of those knots repeat (0.5, 1.5, 1.25, 2.75).
    let even = BSpline1D::new(&x, &y, 2, Bc::Periodic).expect("valid");
    assert_eq!(
        even.as_bspline().knots(),
        [-4.0, -2.75, 0.0, 0.5, 2.0, 3.25, 6.0, 6.5, 8.0]
    );
    // More knots outside than the data has spacings (1 and 2): they go around
    // again. To the left: 0 - 2, - 1, - 2, - 1, - 2.
    let short = BSpline1D::new(&[0.0, 1.0, 3.0], &[1.0, 2.0, 1.0], 5, Bc::Periodic).expect("valid");
    assert_eq!(
        short.as_bspline().knots(),
        [-8.0, -6.0, -5.0, -3.0, -2.0, 0.0, 1.0, 3.0, 4.0, 6.0, 7.0, 9.0, 10.0]
    );
    // The control points repeat with the period: n - 1 of them are different.
    let c = odd.as_bspline().control_points();
    assert_eq!(c.len(), 7);
    for j in 0..3 {
        assert_eq!(c[j].to_bits(), c[j + 4].to_bits(), "control point {j}");
    }
}

#[test]
fn as_many_points_as_the_degree_takes_give_one_polynomial() {
    for degree in 2..=7_usize {
        let x: Vec<f64> = (0..=degree)
            .map(|i| 0.5 * f64::from(u8::try_from(i).expect("small")).powf(1.3))
            .collect();
        // (x / 3)^degree - x.
        let power = i32::try_from(degree).expect("small");
        let f = |v: f64| (v / 3.0).powi(power) - v;
        let y: Vec<f64> = x.iter().map(|v| f(*v)).collect();
        let spline = BSpline1D::new(&x, &y, degree, Bc::NotAKnot).expect("valid");
        assert_eq!(spline.as_bspline().knots().len(), 2 * (degree + 1));
        let hi = x[degree];
        for i in 0..=40 {
            let q = -0.2 * hi + 1.4 * hi * f64::from(i) / 40.0;
            let got = spline.eval(&q).expect("extends");
            let tol = Tol {
                rtol: 1e-9,
                atol: 1e-9,
            };
            assert_close(got, f(q), tol, &format!("degree {degree} at {q}"));
        }
    }
}

#[test]
fn the_derivatives_given_at_the_ends_are_met() {
    let x = [0.0, 0.5, 2.0, 3.0, 4.5];
    let y = [1.0, -1.0, 0.5, 2.0, 0.0];
    let ends = [
        (End::FirstDerivative(0.75), End::SecondDerivative(-2.0)),
        (End::SecondDerivative(1.5), End::FirstDerivative(-0.25)),
        (End::NotAKnot, End::FirstDerivative(3.0)),
        (End::SecondDerivative(-4.0), End::NotAKnot),
    ];
    for (left, right) in ends {
        let spline = BSpline1D::new(&x, &y, 3, Bc::Ends(left, right)).expect("valid");
        let b = spline.as_bspline();
        for (end, at) in [(left, 0.0), (right, 4.5)] {
            let (order, want) = match end {
                End::FirstDerivative(v) => (1, v),
                End::SecondDerivative(v) => (2, v),
                End::NotAKnot => continue,
            };
            let got = b.eval_derivative(at, order).expect("inside");
            assert_close(
                got,
                want,
                DIRECT_F64,
                &format!("{left:?} {right:?} at {at}"),
            );
        }
        for (xi, yi) in x.iter().zip(&y) {
            let got = spline.eval(xi).expect("inside");
            assert_close(got, *yi, DIRECT_F64, &format!("{left:?} {right:?} at {xi}"));
        }
    }
}

#[test]
fn modes_f32_and_thread_safety() {
    let x = [0.0, 1.0, 2.0, 3.0];
    let y = [0.0, 1.0, 4.0, 9.0];
    let spline: BSpline1D<f64> = BSpline1D::new(&x, &y, 2, Bc::NotAKnot).expect("valid");
    assert_eq!(spline.extrapolate(), Extrapolate::Extend);
    assert_eq!(spline.domain().0.to_bits(), 0.0_f64.to_bits());
    assert_eq!(spline.domain().1.to_bits(), 3.0_f64.to_bits());
    assert_eq!(spline.degree(), 2);
    assert_close(spline.eval(&4.0).expect("extends"), 16.0, DIRECT_F64, "x^2");
    let forbid = spline.clone().with_extrapolate(Extrapolate::Forbid);
    assert_eq!(forbid.extrapolate(), Extrapolate::Forbid);
    assert_eq!(forbid.as_bspline().extrapolate(), Extrapolate::Forbid);
    assert_eq!(forbid.eval(&4.0), None);
    assert_eq!(forbid.eval(&f64::NAN), None);
    let clamp = spline.clone().with_extrapolate(Extrapolate::Clamp);
    assert_eq!(bits(clamp.eval(&4.0)), bits(spline.eval(&3.0)));
    let periodic = spline.with_extrapolate(Extrapolate::Periodic);
    assert_eq!(bits(periodic.eval(&3.5)), bits(periodic.eval(&0.5)));

    let single: BSpline1D<f32> = BSpline1D::new(
        &[0.0, 1.0, 2.0, 3.0],
        &[0.0, 1.0, 4.0, 9.0],
        3,
        Bc::NotAKnot,
    )
    .expect("valid");
    let at = single.eval(&1.5).expect("inside");
    assert!((at - 2.25).abs() < 1e-5, "{at}");
    let closed: BSpline1D<f32, [f32; 2]> = BSpline1D::new(
        &[0.0, 1.0, 2.0, 3.0],
        &[[0.0, 1.0], [1.0, 0.0], [0.0, -1.0], [0.0, 1.0]],
        4,
        Bc::Periodic,
    )
    .expect("valid");
    assert_eq!(closed.extrapolate(), Extrapolate::Periodic);

    fn assert_send_sync<T: Send + Sync>() {}
    assert_send_sync::<BSpline1D<f64>>();
    assert_send_sync::<BSpline1D<f32, [f32; 3]>>();
}

// --- Properties ----------------------------------------------------------------

/// Coordinates that increase strictly, `least` at least.
fn coordinates(least: usize) -> impl Strategy<Value = Vec<f64>> {
    (
        -5.0_f64..5.0,
        prop::collection::vec(0.2_f64..2.0, least - 1..14),
    )
        .prop_map(|(start, widths)| {
            let mut x = vec![start];
            for w in widths {
                x.push(x[x.len() - 1] + w);
            }
            x
        })
}

fn values() -> impl Strategy<Value = Vec<f64>> {
    prop::collection::vec(-10.0_f64..10.0, 15)
}

/// The conditions a cubic takes, with the fewest points each needs.
fn cubic_conditions() -> impl Strategy<Value = (Bc<f64>, usize)> {
    let end = prop_oneof![
        Just(End::NotAKnot),
        (-3.0_f64..3.0).prop_map(End::FirstDerivative),
        (-3.0_f64..3.0).prop_map(End::SecondDerivative),
    ];
    prop_oneof![
        Just((Bc::NotAKnot, 4)),
        Just((Bc::Periodic, 2)),
        (end.clone(), end).prop_map(|(a, b)| {
            let free = usize::from(a == End::NotAKnot) + usize::from(b == End::NotAKnot);
            (Bc::Ends(a, b), 2 + free)
        }),
    ]
}

proptest! {
    #[test]
    fn the_spline_goes_through_the_data(
        x in coordinates(9),
        values in values(),
        degree in 1..=8_usize,
        periodic in any::<bool>(),
    ) {
        let n = x.len();
        let mut y = values[..n].to_vec();
        if periodic {
            y[n - 1] = y[0];
        }
        let bc = if periodic { Bc::Periodic } else { Bc::NotAKnot };
        let spline = BSpline1D::new(&x, &y, degree, bc).expect("valid data");
        // A spline of a high degree through uneven points has large control
        // points: the error at a data point is a rounding of those.
        let size = spline.as_bspline().control_points().iter().fold(10.0_f64, |m, c| m.max(c.abs()));
        for (xi, yi) in x.iter().zip(&y) {
            let got = spline.eval(xi).expect("inside");
            prop_assert!((got - yi).abs() <= 1e-12 * size, "degree {} at {}: {} for {}", degree, xi, got, yi);
        }
        if periodic {
            // The derivatives up to `degree - 1` agree at the two ends.
            let b = spline.as_bspline().clone().with_extrapolate(Extrapolate::Extend);
            let period = x[n - 1] - x[0];
            for nu in 1..degree.min(4) {
                let (start, end) = (b.eval_derivative(x[0], nu).expect("inside"), b.eval_derivative(x[n - 1], nu).expect("inside"));
                let scale = size / (0.2 / period).min(1.0).powi(i32::try_from(nu).expect("small")) * 1e3;
                prop_assert!((start - end).abs() <= 1e-9 * scale, "degree {} derivative {}: {} and {}", degree, nu, start, end);
            }
        }
    }

    #[test]
    fn a_polynomial_of_the_degree_is_reproduced(
        x in coordinates(7),
        coefficients in prop::array::uniform6(-2.0_f64..2.0),
        degree in 1..=5_usize,
        ts in prop::collection::vec(-0.1_f64..1.1, 1..30),
    ) {
        let (lo, hi) = (x[0], x[x.len() - 1]);
        // A polynomial of the degree in the variable that runs from 0 to 1 over
        // the data.
        let f = |v: f64| {
            let u = (v - lo) / (hi - lo);
            coefficients[..=degree].iter().rev().fold(0.0, |sum, c| sum * u + c)
        };
        let y: Vec<f64> = x.iter().map(|v| f(*v)).collect();
        let spline = BSpline1D::new(&x, &y, degree, Bc::NotAKnot).expect("valid data");
        for t in &ts {
            let q = lo + t * (hi - lo);
            let got = spline.eval(&q).expect("extends");
            prop_assert!((got - f(q)).abs() <= 1e-9, "degree {} at {}: {} for {}", degree, q, got, f(q));
        }
    }

    #[test]
    fn degree_three_is_the_curve_of_the_cubic_spline(
        x in coordinates(4),
        values in values(),
        (bc, least) in cubic_conditions(),
        ts in prop::collection::vec(0.0_f64..=1.0, 1..30),
    ) {
        prop_assume!(x.len() >= least);
        // `CubicSpline` takes fewer points for not-a-knot (a parabola, a
        // line); here the two types are compared where both build a cubic
        // spline.
        let n = x.len();
        let mut y = values[..n].to_vec();
        if bc == Bc::Periodic {
            y[n - 1] = y[0];
        }
        let spline = BSpline1D::new(&x, &y, 3, bc).expect("valid data");
        let reference = CubicSpline::new(&x, &y, bc).expect("valid data");
        prop_assert_eq!(spline.extrapolate(), reference.extrapolate());
        let (lo, hi) = (x[0], x[n - 1]);
        for t in &ts {
            let q = lo + t * (hi - lo);
            let (got, want) = (spline.eval(&q).expect("inside"), reference.eval(&q).expect("inside"));
            prop_assert!((got - want).abs() <= 1e-9 * want.abs().max(10.0), "{:?} at {}: {} for {}", bc, q, got, want);
        }
    }

    #[test]
    fn degree_one_is_linear_interpolation(
        x in coordinates(2),
        values in values(),
        ts in prop::collection::vec(0.0_f64..=1.0, 1..30),
    ) {
        let n = x.len();
        let y = &values[..n];
        let spline = BSpline1D::new(&x, y, 1, Bc::NotAKnot).expect("valid data");
        let linear = Linear1D::new(&x, y).expect("valid data");
        // The control points are the values, bit for bit. The value at a data
        // point is a rounding away: a basis function there is `w * (1 / w)`.
        // (`Linear1D` returns the data exactly.)
        for (c, yi) in spline.as_bspline().control_points().iter().zip(y) {
            prop_assert_eq!(c.to_bits(), yi.to_bits());
        }
        for (xi, yi) in x.iter().zip(y) {
            let got = spline.eval(xi).expect("inside");
            prop_assert!((got - yi).abs() <= 1e-14 * yi.abs().max(1.0), "at {}: {} for {}", xi, got, yi);
        }
        for t in &ts {
            let q = x[0] + t * (x[n - 1] - x[0]);
            let (got, want) = (spline.eval(&q).expect("inside"), linear.eval(&q).expect("inside"));
            prop_assert!((got - want).abs() <= 1e-12 * want.abs().max(10.0), "at {}: {} for {}", q, got, want);
        }
    }

    #[test]
    fn columns_are_interpolated_on_their_own(
        x in coordinates(9),
        first in values(),
        second in values(),
        degree in 1..=8_usize,
        periodic in any::<bool>(),
        magnitude in -12_i32..=12,
        qs in prop::collection::vec(-12.0_f64..40.0, 1..30),
    ) {
        let n = x.len();
        let factor = 10.0_f64.powi(magnitude);
        let mut first = first[..n].to_vec();
        let mut second: Vec<f64> = second[..n].iter().map(|v| v * factor).collect();
        if periodic {
            first[n - 1] = first[0];
            second[n - 1] = second[0];
        }
        let both: Vec<[f64; 2]> = first.iter().zip(&second).map(|(a, b)| [*a, *b]).collect();
        let pair = BSpline1D::new(&x, &both, degree, if periodic { Bc::Periodic } else { Bc::NotAKnot }).expect("valid data");
        let bc = if periodic { Bc::Periodic } else { Bc::NotAKnot };
        let one = BSpline1D::new(&x, &first, degree, bc).expect("valid data");
        let two = BSpline1D::new(&x, &second, degree, bc).expect("valid data");
        for q in &qs {
            let got = pair.eval(q).expect("extends, or periodic");
            prop_assert_eq!(Some(got[0].to_bits()), bits(one.eval(q)), "x = {}", q);
            prop_assert_eq!(Some(got[1].to_bits()), bits(two.eval(q)), "x = {}", q);
        }
    }

    #[test]
    fn evaluation_is_that_of_the_b_spline(
        x in coordinates(6),
        values in values(),
        degree in 1..=5_usize,
        qs in prop::collection::vec(-12.0_f64..40.0, 1..100),
    ) {
        let y = &values[..x.len()];
        for mode in [Extrapolate::Forbid, Extrapolate::Clamp, Extrapolate::Extend, Extrapolate::Periodic] {
            let spline = BSpline1D::new(&x, y, degree, Bc::NotAKnot).expect("valid data").with_extrapolate(mode);
            let b = spline.as_bspline();
            prop_assert_eq!(b.extrapolate(), mode);
            let mut hint = Hint::new();
            for q in &qs {
                prop_assert_eq!(bits(spline.eval(q)), bits(b.eval(q)), "{:?} x = {}", mode, q);
                prop_assert_eq!(bits(spline.eval_hinted(*q, &mut hint)), bits(b.eval(q)));
            }
            let (mut out, mut want) = (vec![0.0; qs.len()], vec![0.0; qs.len()]);
            prop_assert_eq!(spline.eval_batch_into(&qs, &mut out, -7.0), Ok(()));
            prop_assert_eq!(b.eval_batch_into(&qs, &mut want, -7.0), Ok(()));
            for (got, want) in out.iter().zip(&want) {
                prop_assert_eq!(got.to_bits(), want.to_bits());
            }
        }
    }
}
