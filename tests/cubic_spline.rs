//! `CubicSpline`: the golden data, the input checks through the constructor, few
//! points, the periodic check, the end conditions, and the properties that hold
//! for any cubic spline.

mod common;

use common::{assert_close, load, reals, Real, Tol, DIRECT_F64};
use interpx::{
    BoundaryCondition as Bc, CubicSpline, EndCondition as End, Extrapolate, Hint,
    InterpolateError as E, Interpolator, ParamReason,
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

// --- Golden data -------------------------------------------------------------

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

/// SciPy's `bc_type` for values of `M` columns.
fn bc_of<const M: usize>(v: &Json) -> Option<Bc<[f64; M]>> {
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

/// The mode SciPy gave the case, and whether it is the default for the boundary.
fn mode_of(case: &Json) -> (Extrapolate, bool) {
    match &case["options"]["extrapolate"] {
        Json::Bool(false) => (Extrapolate::Forbid, false),
        Json::String(s) if s == "periodic" => (Extrapolate::Periodic, true),
        _ => (Extrapolate::Extend, true),
    }
}

/// A golden column as a function of the query.
type Column<'a> = &'a dyn Fn(f64) -> Option<f64>;

/// Coordinates, values, a boundary condition, and the error they give.
type Fault<'a> = (&'a [f64], &'a [f64], Bc<f64>, E);

/// The tolerance for the second derivatives of a case, the largest of which is
/// `largest`.
///
/// A second derivative is a difference of slopes divided by a width, so on a
/// short interval it carries the rounding of the slopes many times over, and it
/// is compared with a looser tolerance than a value. The error is also measured
/// against the largest one, not against each on its own: on
/// `square_with_near_point`, where two points are 1e-5 apart, the second
/// derivative runs from 1.7e4 to -1.7e4 across that interval and passes 7 in the
/// middle. There SciPy's and this one differ by 2e-8, which is 1e-12 of the
/// swing and 3e-9 of the 7.
fn second_derivative_tolerance(largest: f64) -> Tol {
    Tol {
        rtol: 1e-10,
        atol: 1e-10 * largest.max(1.0),
    }
}

#[test]
fn scipy_cubic_splines_are_rebuilt_from_their_data() {
    let file = load::<Json>("cubic_spline").expect("golden file cubic_spline");
    let (mut cases, mut small, mut periodic, mut mixed) = (0, 0, 0, 0);
    for case in &file.cases {
        let name = case["name"].as_str().expect("name");
        if case.get("coefficients").is_none() {
            continue;
        }
        let (x, y) = (nums(&case["x"]), nums(&case["y"]));
        let bc = scalar(bc_of::<1>(&case["options"]["bc_type"]).expect("a boundary condition"));
        let (mode, is_default) = mode_of(case);
        let built = CubicSpline::new(&x, &y, bc).expect("valid data");
        if is_default {
            assert_eq!(built.extrapolate(), mode, "{name}: default mode");
        }
        let spline = built.with_extrapolate(mode);

        // The coefficients: SciPy stores `c[power][interval]`.
        let scipy = rows(&case["coefficients"]);
        let poly = spline.as_ppoly();
        assert_eq!(poly.degree(), 3, "{name}");
        assert_eq!(poly.breakpoints().len(), x.len(), "{name}");
        // A coefficient is compared by what it adds to the value at the far end
        // of its interval: the coefficient of `(x - x[i])^m` times `h^m`. Compared
        // as it is, the cubic coefficient of an interval 1e-3 wide would differ
        // from SciPy's by 1e-9, relatively: it is a difference of slopes over
        // `h^2`, and both solvers round the slopes.
        for (i, piece) in poly.coefficients().as_chunks::<4>().0.iter().enumerate() {
            let width = x[i + 1] - x[i];
            for (power, got) in piece.iter().enumerate() {
                let reach = width.powi(3 - i32::try_from(power).expect("0 to 3"));
                let what = format!("{name} interval {i} coefficient {power}");
                assert_close(got * reach, scipy[power][i] * reach, DIRECT_F64, &what);
            }
        }

        // The values, and the two derivatives through the piecewise polynomial.
        let queries = nums(&case["queries"]);
        let largest = nums(&case["derivative_2"])
            .iter()
            .filter(|v| v.is_finite())
            .fold(0.0_f64, |m, v| m.max(v.abs()));
        let columns: [(&str, Column, Tol); 3] = [
            ("expected", &|q| spline.eval(&q), DIRECT_F64),
            ("derivative_1", &|q| poly.eval_derivative(q, 1), DIRECT_F64),
            (
                "derivative_2",
                &|q| poly.eval_derivative(q, 2),
                second_derivative_tolerance(largest),
            ),
        ];
        // The `uneven` data ends with a piece 1e-3 wide. Continued 400 widths
        // past it, the rounding of its cubic coefficient (1e-9, see above) shows
        // in the value, in SciPy's as in this one: those cases are compared
        // inside the data.
        let inside_only = name.contains("uneven");
        for (field, f, tol) in columns {
            for (k, (q, want)) in queries.iter().zip(nums(&case[field])).enumerate() {
                if inside_only && !(x[0]..=x[x.len() - 1]).contains(q) {
                    continue;
                }
                let what = format!("{name} {field} query {k} (x = {q})");
                match f(*q) {
                    None => assert!(want.is_nan(), "{what}: None, SciPy has {want}"),
                    Some(got) => assert_close(got, want, tol, &what),
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
        small += usize::from(x.len() <= 4);
        periodic += usize::from(bc == Bc::Periodic);
        mixed += usize::from(name.contains("_and_"));
    }
    // 32 on the standard data, 12 on two to four points, 5 more.
    assert_eq!((cases, small, periodic, mixed), (49, 12, 5, 6));
}

/// One case with `M` values per point.
// Marked as test code so that `expect` is allowed here as in a `#[test]` function.
#[cfg(test)]
fn check_columns<const M: usize>(case: &Json) {
    let name = case["name"].as_str().expect("name");
    let to_array = |r: &Vec<f64>| -> [f64; M] { r.as_slice().try_into().expect("M columns") };
    let y: Vec<[f64; M]> = rows(&case["y"]).iter().map(to_array).collect();
    let bc = bc_of::<M>(&case["options"]["bc_type"]).expect("a boundary condition");
    let spline = CubicSpline::new(&nums(&case["x"]), &y, bc).expect("valid data");
    assert_eq!(spline.extrapolate(), mode_of(case).0, "{name}");
    let poly = spline.as_ppoly();
    let queries = nums(&case["queries"]);
    for (field, nu) in [("expected", 0), ("derivative_1", 1), ("derivative_2", 2)] {
        for (k, (q, want)) in queries.iter().zip(rows(&case[field])).enumerate() {
            let got = poly.eval_derivative(*q, nu).expect("the case extrapolates");
            for column in 0..M {
                let what = format!("{name} {field} query {k} column {column}");
                assert_close(got[column], want[column], DIRECT_F64, &what);
            }
            if nu == 0 {
                let direct = spline.eval(q).expect("the case extrapolates");
                assert_eq!(direct.map(f64::to_bits), got.map(f64::to_bits), "{name}");
            }
        }
    }
}

#[test]
fn scipy_cubic_splines_of_two_columns_are_rebuilt_column_by_column() {
    let file = load::<Json>("cubic_spline").expect("golden file cubic_spline");
    let mut cases = 0;
    for case in &file.cases {
        if case.get("coefficients").is_none() {
            check_columns::<2>(case);
            cases += 1;
        }
    }
    assert_eq!(cases, 4);
}

// --- Input checks through the constructor --------------------------------------

#[test]
fn the_constructor_reports_each_fault() {
    let ok = CubicSpline::new(&[0.0, 1.0], &[1.0, 2.0], Bc::default());
    assert!(ok.is_ok(), "two points are enough");

    let nan_end = Bc::Ends(End::FirstDerivative(f64::NAN), End::NotAKnot);
    let infinite_end = Bc::Ends(End::NotAKnot, End::SecondDerivative(f64::INFINITY));
    let not_finite = E::InvalidParameter {
        name: "bc",
        reason: ParamReason::NonFinite,
    };
    let cases: [Fault; 10] = [
        (
            &[0.0, 1.0, 2.0],
            &[1.0, 2.0],
            Bc::NotAKnot,
            E::LengthMismatch {
                expected: 3,
                actual: 2,
            },
        ),
        (
            &[0.0],
            &[1.0],
            Bc::NotAKnot,
            E::InsufficientPoints {
                required: 2,
                actual: 1,
            },
        ),
        (
            &[0.0, f64::NAN, 2.0],
            &[1.0, 2.0, 3.0],
            Bc::NotAKnot,
            E::NonFinite { index: 1 },
        ),
        (
            &[0.0, 1.0, 2.0],
            &[1.0, 2.0, f64::INFINITY],
            Bc::NotAKnot,
            E::NonFinite { index: 2 },
        ),
        (
            &[0.0, 2.0, 1.0],
            &[1.0, 2.0, 3.0],
            Bc::NotAKnot,
            E::NotSorted { index: 2 },
        ),
        (
            &[0.0, 1.0, 1.0],
            &[1.0, 2.0, 3.0],
            Bc::NotAKnot,
            E::DuplicatePoints {
                first: 1,
                second: 2,
            },
        ),
        (&[0.0, 1.0, 2.0], &[1.0, 2.0, 3.0], nan_end, not_finite),
        (&[0.0, 1.0, 2.0], &[1.0, 2.0, 3.0], infinite_end, not_finite),
        (
            &[0.0, 1.0, 2.0],
            &[1.0, 2.0, 3.0],
            Bc::Periodic,
            E::NotPeriodic,
        ),
        // Finite values whose difference is not: the secant to point 2.
        (
            &[0.0, 1.0, 2.0],
            &[0.0, -1e308, 1e308],
            Bc::NotAKnot,
            E::NonFinite { index: 2 },
        ),
    ];
    for (x, y, bc, want) in cases {
        assert_eq!(CubicSpline::new(x, y, bc).err(), Some(want));
    }
}

#[test]
fn the_constructor_keeps_the_order_of_the_checks() {
    let nan_end = Bc::Ends(End::FirstDerivative(f64::NAN), End::NotAKnot);
    // Unsorted coordinates, and a NaN in the boundary condition.
    assert_eq!(
        CubicSpline::new(&[1.0, 0.0], &[1.0, 2.0], nan_end).err(),
        Some(E::NotSorted { index: 1 })
    );
    // A NaN value, and end values that differ.
    assert_eq!(
        CubicSpline::new(&[0.0, 1.0, 2.0], &[1.0, f64::NAN, 3.0], Bc::Periodic).err(),
        Some(E::NonFinite { index: 1 })
    );
    // End values that differ, and a secant that overflows.
    assert_eq!(
        CubicSpline::new(&[0.0, 1.0, 2.0], &[1e308, -1e308, 0.0], Bc::Periodic).err(),
        Some(E::NotPeriodic)
    );
}

#[test]
fn an_overflow_while_solving_is_an_error_and_not_a_panic() {
    // Points 1e-160 apart: the secants are finite, the coefficients are not.
    let x = [0.0, 1e-160, 2e-160, 3e-160, 4e-160];
    let y = [0.0, 1.0, 0.0, 1.0, 0.0];
    for bc in [Bc::NotAKnot, Bc::natural(), Bc::clamped()] {
        let err = CubicSpline::new(&x, &y, bc).err();
        assert!(
            matches!(err, Some(E::NonFinite { index }) if index < 5),
            "{bc:?}: {err:?}"
        );
    }
    let err = CubicSpline::new(&x, &y, Bc::Periodic).err();
    assert!(
        matches!(err, Some(E::NonFinite { index }) if index < 5),
        "{err:?}"
    );
    // Widths near the largest number: the entries of the system overflow.
    let wide = [-1.5e308, -0.5e308, 0.5e308, 1.5e308];
    for bc in [Bc::NotAKnot, Bc::natural(), Bc::Periodic] {
        let err = CubicSpline::new(&wide, &[0.0, 1.0, 1.0, 0.0], bc).err();
        assert!(
            matches!(err, Some(E::NonFinite { index }) if index < 4),
            "{bc:?}: {err:?}"
        );
    }
}

#[test]
fn an_overflowing_coefficient_names_its_interval() {
    // Values that swing by 1e300. On the two intervals of width 1 every
    // coefficient is finite. On the third, 1e-3 wide, the secant is 1e303 and the
    // cubic coefficient, of the order of a secant over the width squared, is not.
    let x = [0.0, 1.0, 2.0, 2.001, 2.002];
    let y = [0.0, 1e300, 0.0, 1e300, 0.0];
    assert!(x[3] > x[2], "the two coordinates differ");
    for bc in [Bc::NotAKnot, Bc::natural(), Bc::clamped()] {
        assert_eq!(
            CubicSpline::new(&x, &y, bc).err(),
            Some(E::NonFinite { index: 2 }),
            "{bc:?}"
        );
    }
}

// --- The periodic check ----------------------------------------------------------

#[test]
fn a_periodic_spline_needs_equal_end_values_up_to_rounding() {
    let x = [0.0, 1.0, 2.0, 3.0];
    let ends =
        |first: f64, last: f64| CubicSpline::new(&x, &[first, 2.0, -1.0, last], Bc::Periodic);
    assert!(ends(1.0, 1.0).is_ok());
    // SciPy's tolerance is 1e-15 * (1 + |last|): 2e-15 around 1.
    assert!(ends(1.0, 1.0 + 4.0 * f64::EPSILON).is_ok());
    assert_eq!(
        ends(1.0, 1.0 + 16.0 * f64::EPSILON).err(),
        Some(E::NotPeriodic)
    );
    assert_eq!(
        ends(1.0, 1.0 - 16.0 * f64::EPSILON).err(),
        Some(E::NotPeriodic)
    );
    // Around zero the absolute part decides.
    assert!(ends(0.0, 5e-16).is_ok());
    assert_eq!(ends(0.0, 2e-15).err(), Some(E::NotPeriodic));
    // The sine of a whole period is not zero in floating point.
    let period = [0.0, 0.5, 1.0, 1.5, 2.0].map(|t: f64| t * std::f64::consts::PI);
    assert!(period[4].sin().abs() > 0.0);
    assert!(CubicSpline::new(&period, &period.map(f64::sin), Bc::Periodic).is_ok());

    // Every column has to close.
    let closed = [[1.0, 5.0], [2.0, 6.0], [0.0, 7.0], [1.0, 5.0]];
    assert!(CubicSpline::new(&x, &closed, Bc::Periodic).is_ok());
    let open = [[1.0, 5.0], [2.0, 6.0], [0.0, 7.0], [1.0, 5.5]];
    assert_eq!(
        CubicSpline::new(&x, &open, Bc::Periodic).err(),
        Some(E::NotPeriodic)
    );

    // In single precision the same multiple of its own rounding unit.
    let x32 = [0.0_f32, 1.0, 2.0, 3.0];
    let ends = |last: f32| CubicSpline::new(&x32, &[1.0, 2.0, -1.0, last], Bc::Periodic);
    assert!(ends(1.0 + 4.0 * f32::EPSILON).is_ok());
    assert_eq!(ends(1.0 + 16.0 * f32::EPSILON).err(), Some(E::NotPeriodic));
}

// --- Few points --------------------------------------------------------------------

#[test]
fn two_points_with_not_a_knot_give_the_line_through_them() {
    let spline: CubicSpline<f64> =
        CubicSpline::new(&[1.0, 3.0], &[2.0, 8.0], Bc::NotAKnot).expect("valid");
    let c = spline.as_ppoly().coefficients();
    assert_eq!(c.len(), 4);
    assert_eq!(c.map_bits(), [0.0, 0.0, 3.0, 2.0].map(f64::to_bits));
    assert_eq!(bits(spline.eval(&2.0)), bits(Some(5.0)));
    assert_eq!(bits(spline.eval(&5.0)), bits(Some(14.0)), "extended");

    // One not-a-knot end takes the slope of the line; the other end keeps its own.
    let half: CubicSpline<f64> = CubicSpline::new(
        &[1.0, 3.0],
        &[2.0, 8.0],
        Bc::Ends(End::NotAKnot, End::FirstDerivative(0.0)),
    )
    .expect("valid");
    let poly = half.as_ppoly();
    assert_close(
        poly.eval_derivative(1.0, 1).expect("inside"),
        3.0,
        DIRECT_F64,
        "left",
    );
    assert_close(
        poly.eval_derivative(3.0, 1).expect("inside"),
        0.0,
        DIRECT_F64,
        "right",
    );

    // Periodic through two equal values: a constant.
    let flat: CubicSpline<f64> =
        CubicSpline::new(&[1.0, 3.0], &[4.0, 4.0], Bc::Periodic).expect("valid");
    assert_eq!(flat.extrapolate(), Extrapolate::Periodic);
    for x in [1.0, 1.7, 3.0, 9.5, -4.0] {
        assert_eq!(bits(flat.eval(&x)), bits(Some(4.0)), "x = {x}");
    }
    // End values that differ by rounding pass the check, and are joined by the
    // line through them: the slope at both ends is their secant, sign included.
    let rise = 8.0 * f64::EPSILON;
    let nearly: CubicSpline<f64> =
        CubicSpline::new(&[1.0, 3.0], &[4.0, 4.0 + rise], Bc::Periodic).expect("within rounding");
    let c = nearly.as_ppoly().coefficients();
    assert_eq!(c[2].to_bits(), (rise / 2.0).to_bits(), "the slope");
    assert_eq!(
        c.map_bits()[..2],
        [0.0_f64, 0.0].map(f64::to_bits),
        "a line"
    );
}

/// The coefficients as bits, for an exact comparison.
trait MapBits {
    fn map_bits(&self) -> [u64; 4];
}

impl MapBits for [f64] {
    fn map_bits(&self) -> [u64; 4] {
        [self[0], self[1], self[2], self[3]].map(f64::to_bits)
    }
}

#[test]
fn three_points_with_not_a_knot_give_the_parabola_through_them() {
    // x^2 - 2x + 3 at 0, 1, 4.
    let f = |x: f64| x * x - 2.0 * x + 3.0;
    let x = [0.0, 1.0, 4.0];
    let spline = CubicSpline::new(&x, &x.map(f), Bc::NotAKnot).expect("valid");
    for piece in spline.as_ppoly().coefficients().as_chunks::<4>().0 {
        assert_close(piece[0], 0.0, DIRECT_F64, "no cubic term");
        assert_close(piece[1], 1.0, DIRECT_F64, "the square term");
    }
    for q in [-1.0, 0.5, 2.5, 4.0, 6.0] {
        assert_close(
            spline.eval(&q).expect("extends"),
            f(q),
            DIRECT_F64,
            "the parabola",
        );
    }
}

#[test]
fn three_points_periodic_have_one_slope_and_close_smoothly() {
    let spline =
        CubicSpline::new(&[0.0, 1.0, 3.0], &[1.0, -1.0, 1.0], Bc::Periodic).expect("valid");
    let poly = spline.as_ppoly();
    // SciPy's formula gives every point the same slope.
    let slopes = [0.0, 1.0, 3.0].map(|x| {
        poly.clone()
            .with_extrapolate(Extrapolate::Extend)
            .eval_derivative(x, 1)
            .expect("inside")
    });
    assert_close(slopes[1], slopes[0], DIRECT_F64, "slope at the middle");
    assert_close(slopes[2], slopes[0], DIRECT_F64, "slope at the end");
    // (d0 / h0 + d1 / h1) / (1 / h0 + 1 / h1) with secants -2 and 1, widths 1 and 2.
    assert_close(slopes[0], (-2.0 + 0.5) / 1.5, DIRECT_F64, "the formula");
    assert_eq!(
        bits(spline.eval(&3.5)),
        bits(spline.eval(&0.5)),
        "one period on"
    );
}

// --- End conditions and modes ----------------------------------------------------

#[test]
fn each_end_condition_is_met() {
    let x = [0.0, 0.5, 1.5, 2.0, 3.5, 4.0];
    let y = [1.0, -0.5, 2.0, 0.25, 0.0, 1.5];
    let at = |bc: Bc<f64>, q: f64, nu: usize| {
        let spline = CubicSpline::new(&x, &y, bc).expect("valid");
        spline
            .as_ppoly()
            .clone()
            .with_extrapolate(Extrapolate::Extend)
            .eval_derivative(q, nu)
            .expect("extends")
    };
    let close = |got: f64, want: f64, what: &str| {
        assert!((got - want).abs() <= 1e-10, "{what}: {got} for {want}");
    };

    let natural = Bc::natural();
    close(at(natural, 0.0, 2), 0.0, "natural, left");
    close(at(natural, 4.0, 2), 0.0, "natural, right");
    let clamped = Bc::clamped();
    close(at(clamped, 0.0, 1), 0.0, "clamped, left");
    close(at(clamped, 4.0, 1), 0.0, "clamped, right");
    let mixed = Bc::Ends(End::FirstDerivative(0.75), End::SecondDerivative(-1.25));
    close(at(mixed, 0.0, 1), 0.75, "first derivative on the left");
    close(at(mixed, 4.0, 2), -1.25, "second derivative on the right");
    let mixed = Bc::Ends(End::SecondDerivative(2.5), End::FirstDerivative(-0.5));
    close(at(mixed, 0.0, 2), 2.5, "second derivative on the left");
    close(at(mixed, 4.0, 1), -0.5, "first derivative on the right");

    // Not-a-knot: the third derivative does not jump at the second point, nor
    // at the one before the last. It does at the points between.
    let third = |bc, i: usize| {
        let spline = CubicSpline::new(&x, &y, bc).expect("valid");
        let c = spline.as_ppoly().coefficients();
        (c[4 * (i - 1)], c[4 * i])
    };
    let (before, after) = third(Bc::NotAKnot, 1);
    close(before, after, "not-a-knot at the second point");
    let (before, after) = third(Bc::NotAKnot, 4);
    close(before, after, "not-a-knot at the point before the last");
    let (before, after) = third(Bc::NotAKnot, 2);
    assert!((before - after).abs() > 1e-3, "an inner point is a knot");
    // One end only.
    let one = Bc::Ends(End::NotAKnot, End::FirstDerivative(0.0));
    let (before, after) = third(one, 1);
    close(before, after, "not-a-knot on the left only");
    let (before, after) = third(one, 4);
    assert!((before - after).abs() > 1e-3, "the right end is clamped");
    close(at(one, 4.0, 1), 0.0, "and its slope is zero");

    // Periodic: value, slope, and curvature agree at the two ends.
    let closed = [1.0, -0.5, 2.0, 0.25, 0.0, 1.0];
    let spline = CubicSpline::new(&x, &closed, Bc::Periodic).expect("valid");
    let poly = spline
        .as_ppoly()
        .clone()
        .with_extrapolate(Extrapolate::Extend);
    for nu in 0..3 {
        let (left, right) = (poly.eval_derivative(0.0, nu), poly.eval_derivative(4.0, nu));
        close(right.expect("inside"), left.expect("inside"), "periodic");
    }
}

#[test]
fn the_default_mode_follows_the_boundary_and_can_be_changed() {
    let x = [0.0, 1.0, 2.0, 3.0];
    let y = [1.0, 2.0, 0.0, 1.0];
    for bc in [Bc::NotAKnot, Bc::natural(), Bc::clamped()] {
        let spline = CubicSpline::new(&x, &y, bc).expect("valid");
        assert_eq!(spline.extrapolate(), Extrapolate::Extend, "{bc:?}");
        assert!(spline.eval(&4.0).is_some());
    }
    let spline = CubicSpline::new(&x, &y, Bc::Periodic).expect("valid");
    assert_eq!(spline.extrapolate(), Extrapolate::Periodic);
    assert_eq!(spline.as_ppoly().extrapolate(), Extrapolate::Periodic);
    assert_eq!(bits(spline.eval(&3.5)), bits(spline.eval(&0.5)));
    assert_eq!(spline.domain().0.to_bits(), 0.0_f64.to_bits());
    assert_eq!(spline.domain().1.to_bits(), 3.0_f64.to_bits());

    let forbid = spline.clone().with_extrapolate(Extrapolate::Forbid);
    assert_eq!(forbid.extrapolate(), Extrapolate::Forbid);
    assert_eq!(forbid.eval(&3.5), None);
    assert_eq!(forbid.eval(&f64::NAN), None);
    let clamp = spline.with_extrapolate(Extrapolate::Clamp);
    assert_eq!(bits(clamp.eval(&-5.0)), bits(clamp.eval(&0.0)));
}

#[test]
fn array_values_f32_and_thread_safety() {
    // Two columns in single precision: x^2 and 1 - x, both reproduced by the
    // parabola through three points.
    let x = [0.0_f32, 1.0, 3.0];
    let y = x.map(|x| [x * x, 1.0 - x]);
    let spline = CubicSpline::new(&x, &y, Bc::default()).expect("valid");
    let got = spline.eval(&2.0).expect("inside");
    assert!(
        (got[0] - 4.0).abs() <= 1e-5 && (got[1] + 1.0).abs() <= 1e-5,
        "{got:?}"
    );

    // A slope per column at each end.
    let ends = Bc::Ends(
        End::FirstDerivative([0.0, -1.0]),
        End::FirstDerivative([6.0, -1.0]),
    );
    let spline = CubicSpline::new(&x, &y, ends).expect("valid");
    let got = spline.eval(&2.0).expect("inside");
    assert!(
        (got[0] - 4.0).abs() <= 1e-5 && (got[1] + 1.0).abs() <= 1e-5,
        "{got:?}"
    );
    let nan = Bc::Ends(End::FirstDerivative([0.0, f32::NAN]), End::NotAKnot);
    assert_eq!(
        CubicSpline::new(&x, &y, nan).err(),
        Some(E::InvalidParameter {
            name: "bc",
            reason: ParamReason::NonFinite,
        })
    );

    fn assert_send_sync<T: Send + Sync>() {}
    assert_send_sync::<CubicSpline<f64>>();
    assert_send_sync::<CubicSpline<f32, [f32; 3]>>();
}

// --- Properties ----------------------------------------------------------------

/// Coordinates that increase strictly, at least `least` of them.
fn coordinates(least: usize) -> impl Strategy<Value = Vec<f64>> {
    (
        -5.0_f64..5.0,
        prop::collection::vec(0.1_f64..3.0, least - 1..12),
    )
        .prop_map(|(start, widths)| {
            let mut x = vec![start];
            for w in widths {
                x.push(x[x.len() - 1] + w);
            }
            x
        })
}

/// A boundary condition of every kind but the periodic one.
fn conditions() -> impl Strategy<Value = Bc<f64>> {
    let end = prop_oneof![
        Just(End::NotAKnot),
        (-3.0_f64..3.0).prop_map(End::FirstDerivative),
        (-3.0_f64..3.0).prop_map(End::SecondDerivative),
    ];
    (end.clone(), end).prop_map(|(a, b)| Bc::Ends(a, b))
}

/// The value and the first two derivatives of piece `i` at distance `dx` from
/// its left end, written apart from the library.
fn piece(c: &[f64], i: usize, dx: f64) -> [f64; 3] {
    let c = &c[4 * i..4 * i + 4];
    [
        ((c[0] * dx + c[1]) * dx + c[2]) * dx + c[3],
        (3.0 * c[0] * dx + 2.0 * c[1]) * dx + c[2],
        6.0 * c[0] * dx + 2.0 * c[1],
    ]
}

proptest! {
    #[test]
    fn the_spline_goes_through_the_data_and_is_twice_differentiable(
        x in coordinates(2),
        values in prop::collection::vec(-10.0_f64..10.0, 12),
        bc in conditions(),
        periodic in any::<bool>(),
    ) {
        let n = x.len();
        let mut y = values[..n].to_vec();
        let bc = if periodic {
            y[n - 1] = y[0];
            Bc::Periodic
        } else {
            bc
        };
        let spline = CubicSpline::new(&x, &y, bc).expect("valid data");
        let c = spline.as_ppoly().coefficients();
        // Widths are 0.1 at least and values 10 at most: a second derivative is
        // of the order of 10 / 0.1^2, times the growth of the solve.
        let tolerance = 1e-8;
        for i in 0..n - 1 {
            // The left end of a piece has the data value as its constant term.
            prop_assert_eq!(c[4 * i + 3].to_bits(), y[i].to_bits());
            let at_right = piece(c, i, x[i + 1] - x[i]);
            prop_assert!((at_right[0] - y[i + 1]).abs() <= tolerance, "piece {} ends at {} for {}", i, at_right[0], y[i + 1]);
            if i + 2 < n {
                // The next piece starts with the same slope and curvature.
                let next = piece(c, i + 1, 0.0);
                prop_assert!((at_right[1] - next[1]).abs() <= tolerance, "slope at point {}: {} and {}", i + 1, at_right[1], next[1]);
                prop_assert!((at_right[2] - next[2]).abs() <= tolerance * 100.0, "curvature at point {}: {} and {}", i + 1, at_right[2], next[2]);
            }
        }
        // The ends.
        let (first, last) = (piece(c, 0, 0.0), piece(c, n - 2, x[n - 1] - x[n - 2]));
        match bc {
            Bc::Periodic => {
                prop_assert!((first[1] - last[1]).abs() <= tolerance, "periodic slope: {} and {}", first[1], last[1]);
                prop_assert!((first[2] - last[2]).abs() <= tolerance * 100.0, "periodic curvature: {} and {}", first[2], last[2]);
            }
            Bc::Ends(start, end) => {
                for (condition, got) in [(start, first), (end, last)] {
                    match condition {
                        End::FirstDerivative(v) => prop_assert!((got[1] - v).abs() <= tolerance, "slope {} for {}", got[1], v),
                        End::SecondDerivative(v) => prop_assert!((got[2] - v).abs() <= tolerance * 100.0, "curvature {} for {}", got[2], v),
                        End::NotAKnot => {}
                    }
                }
            }
            Bc::NotAKnot => {}
        }
    }

    #[test]
    fn a_cubic_polynomial_is_reproduced(
        x in coordinates(4),
        a in prop::array::uniform4(-2.0_f64..2.0),
        qs in prop::collection::vec(0.0_f64..=1.0, 1..20),
    ) {
        let f = |t: f64| ((a[0] * t + a[1]) * t + a[2]) * t + a[3];
        let df = |t: f64| (3.0 * a[0] * t + 2.0 * a[1]) * t + a[2];
        let d2f = |t: f64| 6.0 * a[0] * t + 2.0 * a[1];
        let y: Vec<f64> = x.iter().map(|t| f(*t)).collect();
        let (lo, hi) = (x[0], x[x.len() - 1]);
        // Not-a-knot reproduces any cubic; so do the exact end derivatives.
        let exact = [
            Bc::NotAKnot,
            Bc::Ends(End::FirstDerivative(df(lo)), End::FirstDerivative(df(hi))),
            Bc::Ends(End::SecondDerivative(d2f(lo)), End::SecondDerivative(d2f(hi))),
            Bc::Ends(End::NotAKnot, End::FirstDerivative(df(hi))),
            Bc::Ends(End::SecondDerivative(d2f(lo)), End::NotAKnot),
        ];
        let scale = y.iter().fold(1.0_f64, |m, v| m.max(v.abs()));
        for bc in exact {
            let spline = CubicSpline::new(&x, &y, bc).expect("valid data");
            for t in &qs {
                let q = lo + t * (hi - lo);
                let got = spline.eval(&q).expect("inside");
                prop_assert!((got - f(q)).abs() <= 1e-9 * scale, "{:?} at {}: {} for {}", bc, q, got, f(q));
            }
        }
    }

    #[test]
    fn evaluation_is_that_of_the_piecewise_polynomial(
        x in coordinates(2),
        values in prop::collection::vec(-10.0_f64..10.0, 12),
        bc in conditions(),
        qs in prop::collection::vec(-12.0_f64..40.0, 1..100),
    ) {
        let y = &values[..x.len()];
        for mode in [Extrapolate::Forbid, Extrapolate::Clamp, Extrapolate::Extend, Extrapolate::Periodic] {
            let spline = CubicSpline::new(&x, y, bc).expect("valid data").with_extrapolate(mode);
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
