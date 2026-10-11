//! Derivatives and integrals of `PPoly` and `BSpline`: the golden data, the
//! cases where the two types are made to agree and SciPy's are not, the four
//! extrapolation modes, and the properties that tie the four functions together.

mod common;

use common::{assert_close, load, reals, Real, Tol, DIRECT_F64};
use interpx::{
    BSpline, Extrapolate, InterpolateError as E, Interpolator, PPoly, ParamReason,
    BSPLINE_MAX_DEGREE,
};
use proptest::prelude::*;
use serde_json::Value as Json;

const MODES: [Extrapolate; 4] = [
    Extrapolate::Forbid,
    Extrapolate::Clamp,
    Extrapolate::Extend,
    Extrapolate::Periodic,
];

/// The tolerance for a value that went through a change of basis (see
/// `tests/bspline.rs`).
const CONVERTED: Tol = Tol {
    rtol: 1e-9,
    atol: 1e-9,
};

fn bits(v: Option<f64>) -> Option<u64> {
    v.map(f64::to_bits)
}

/// The bits of `0 - v`: what reversed limits make of an integral `v`. (For a
/// zero this is `+0`, where `-v` would be `-0`.)
fn reversed(v: Option<f64>) -> Option<u64> {
    bits(v.map(|v| 0.0 - v))
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

/// SciPy's `c[power][interval]` as `[interval][power]`, and the degree.
fn interval_major(coefficients: &Json) -> (Vec<f64>, usize) {
    let rows = rows(coefficients);
    let intervals = rows.first().map_or(0, Vec::len);
    let flat = (0..intervals)
        .flat_map(|i| rows.iter().map(move |row| row[i]))
        .collect();
    (flat, rows.len().saturating_sub(1))
}

fn ppoly_of(case: &Json) -> Result<PPoly<f64>, E> {
    let (flat, degree) = interval_major(&case["coefficients"]);
    let mode = match &case["options"]["extrapolate"] {
        Json::Bool(true) => Extrapolate::Extend,
        Json::String(s) if s == "periodic" => Extrapolate::Periodic,
        // `false`, and `None` for Akima, whose default is no extrapolation.
        _ => Extrapolate::Forbid,
    };
    Ok(PPoly::from_owned(nums(&case["breakpoints"]), flat, degree)?.with_extrapolate(mode))
}

fn degree_of(case: &Json) -> usize {
    usize::try_from(case["k"].as_u64().unwrap_or(u64::MAX)).unwrap_or(usize::MAX)
}

fn bspline_of(case: &Json) -> Result<BSpline<f64>, E> {
    let mode = match &case["options"]["extrapolate"] {
        Json::Bool(false) => Extrapolate::Forbid,
        Json::String(s) if s == "periodic" => Extrapolate::Periodic,
        _ => Extrapolate::Extend,
    };
    Ok(BSpline::from_owned(
        nums(&case["knots"]),
        nums(&case["control_points"]),
        degree_of(case),
    )?
    .with_extrapolate(mode))
}

/// Checks one column of golden values against a function of the query. SciPy
/// has NaN where it does not extrapolate; the function has `None` there.
fn check_column(
    queries: &[f64],
    expected: &[f64],
    tol: Tol,
    what: &str,
    f: impl Fn(f64) -> Option<f64>,
) {
    assert_eq!(queries.len(), expected.len(), "{what}");
    for (k, (x, want)) in queries.iter().zip(expected).enumerate() {
        let what = format!("{what}, query {k} (x = {x})");
        match f(*x) {
            None => assert!(want.is_nan(), "{what}: None, SciPy has {want}"),
            Some(got) => assert_close(got, *want, tol, &what),
        }
    }
}

/// `[[a, b, integral], ...]` against a function of the limits.
// Marked as test code so that `panic!` is allowed here as in a `#[test]` function.
#[cfg(test)]
fn check_integrals(v: &Json, what: &str, f: impl Fn(f64, f64) -> Option<f64>) -> usize {
    let triples = rows(v);
    for t in &triples {
        let (a, b, want) = (t[0], t[1], t[2]);
        let what = format!("{what}, integral from {a} to {b}");
        match f(a, b) {
            None => panic!("{what}: None, SciPy has {want}"),
            Some(got) => assert_close(got, want, DIRECT_F64, &what),
        }
    }
    triples.len()
}

// --- Golden data: the fields every 1D case has ----------------------------------

#[test]
fn scipy_piecewise_polynomials_differentiate_and_integrate_to_the_golden_values() {
    let file = load::<Json>("oned").expect("golden file oned");
    let (mut cases, mut integrals, mut modes) = (0, 0, std::collections::BTreeSet::new());
    for case in &file.cases {
        let name = case["name"].as_str().expect("name");
        if case.get("coefficients").is_none() || case.get("expected").is_none() {
            continue;
        }
        let p = ppoly_of(case).expect("SciPy's coefficients are valid");
        let queries = nums(&case["queries"]);

        for nu in [1, 2] {
            let expected = nums(&case[format!("derivative_{nu}").as_str()]);
            check_column(
                &queries,
                &expected,
                DIRECT_F64,
                &format!("{name} eval_derivative {nu}"),
                |x| p.eval_derivative(x, nu),
            );
            let d = p.derivative(nu).expect("finite");
            check_column(
                &queries,
                &expected,
                DIRECT_F64,
                &format!("{name} derivative({nu})"),
                |x| d.eval(&x),
            );
        }
        let primitive = p.antiderivative(1).expect("finite");
        check_column(
            &queries,
            &nums(&case["antiderivative"]),
            DIRECT_F64,
            &format!("{name} antiderivative"),
            |x| primitive.eval(&x),
        );
        integrals += check_integrals(&case["integrals"], name, |a, b| p.integrate(a, b));

        cases += 1;
        modes.insert(format!("{:?}", p.extrapolate()));
    }
    assert_eq!(cases, 58, "the cases tests/ppoly.rs evaluates");
    assert_eq!(integrals, 58 * 4);
    assert_eq!(modes.len(), 3, "extend, forbid, and periodic: {modes:?}");
}

#[test]
fn scipy_b_splines_differentiate_and_integrate_to_the_golden_values() {
    let file = load::<Json>("oned").expect("golden file oned");
    let (mut cases, mut second, mut integrals) = (0, 0, 0);
    for case in &file.cases {
        let name = case["name"].as_str().expect("name");
        if case.get("knots").is_none() || case.get("expected").is_none() {
            continue;
        }
        let spline = bspline_of(case).expect("SciPy's B-spline is valid");
        let queries = nums(&case["queries"]);

        for nu in [1, 2] {
            let Some(expected) = case.get(format!("derivative_{nu}").as_str()) else {
                continue;
            };
            let expected = nums(expected);
            check_column(
                &queries,
                &expected,
                DIRECT_F64,
                &format!("{name} eval_derivative {nu}"),
                |x| spline.eval_derivative(x, nu),
            );
            // The control points of a derivative are differences of control
            // points over knot spacings: a change of basis.
            let d = spline.derivative(nu).expect("finite");
            check_column(
                &queries,
                &expected,
                CONVERTED,
                &format!("{name} derivative({nu})"),
                |x| d.eval(&x),
            );
            second += usize::from(nu == 2);
        }
        if let Some(expected) = case.get("antiderivative") {
            let primitive = spline.antiderivative(1).expect("finite");
            check_column(
                &queries,
                &nums(expected),
                DIRECT_F64,
                &format!("{name} antiderivative"),
                |x| primitive.eval(&x),
            );
            integrals += check_integrals(&case["integrals"], name, |a, b| spline.integrate(a, b));
        }
        cases += 1;
    }
    assert_eq!(cases, 33, "the cases tests/bspline.rs evaluates");
    // Degree 1 has no second derivative in the file, nor has the case without
    // extrapolation; the periodic case and that one have no integrals.
    assert_eq!(second, 27);
    assert_eq!(integrals, 31 * 4);
}

/// The derivatives of one curve case with `D` outputs per point.
// Marked as test code so that `expect` is allowed here as in a `#[test]` function.
#[cfg(test)]
fn check_curve<const D: usize>(case: &Json) -> usize {
    let name = case["name"].as_str().expect("name");
    let to_array = |r: &Vec<f64>| -> [f64; D] { r.as_slice().try_into().expect("D columns") };
    let control_points: Vec<[f64; D]> =
        rows(&case["control_points"]).iter().map(to_array).collect();
    let spline = BSpline::from_owned(nums(&case["knots"]), control_points, degree_of(case))
        .expect("SciPy's B-spline is valid")
        .with_extrapolate(match &case["options"]["extrapolate"] {
            Json::String(s) if s == "periodic" => Extrapolate::Periodic,
            _ => Extrapolate::Extend,
        });
    let queries = nums(&case["queries"]);
    let mut checked = 0;
    for nu in [1, 2] {
        let Some(expected) = case.get(format!("derivative_{nu}").as_str()) else {
            continue;
        };
        let d = spline.derivative(nu).expect("finite");
        for (k, (x, want)) in queries.iter().zip(rows(expected)).enumerate() {
            let direct = spline
                .eval_derivative(*x, nu)
                .expect("the curves extrapolate");
            let through = d.eval(x).expect("the curves extrapolate");
            for column in 0..D {
                let what = format!("{name} derivative {nu} query {k} column {column}");
                assert_close(direct[column], want[column], DIRECT_F64, &what);
                assert_close(through[column], want[column], DIRECT_F64, &what);
            }
        }
        checked += 1;
    }
    checked
}

#[test]
fn scipy_curves_differentiate_column_by_column() {
    let file = load::<Json>("curve").expect("golden file curve");
    let mut checked = 0;
    for case in &file.cases {
        checked += match rows(&case["control_points"]).first().map_or(0, Vec::len) {
            2 => check_curve::<2>(case),
            3 => check_curve::<3>(case),
            other => panic!("a curve with {other} columns"),
        };
    }
    assert!(checked >= 20, "{checked} derivative columns of curves");
}

// --- Golden data: derivatives and antiderivatives as objects ---------------------

#[test]
fn derivatives_and_antiderivatives_of_piecewise_polynomials_match_scipys_coefficients() {
    let file = load::<Json>("calculus").expect("golden file calculus");
    let (mut cases, mut periodic) = (0, 0);
    for case in &file.cases {
        let name = case["name"].as_str().expect("name");
        if !name.starts_with("calculus/ppoly/") {
            continue;
        }
        let p = ppoly_of(case).expect("valid");
        let queries = nums(&case["queries"]);

        let derivatives = case["derivatives"].as_array().expect("derivatives");
        assert_eq!(derivatives.len(), p.degree() + 1, "{name}");
        for entry in derivatives {
            let nu = usize::try_from(entry["nu"].as_u64().expect("nu")).expect("usize");
            let d = p.derivative(nu).expect("finite");
            let (want, degree) = interval_major(&entry["coefficients"]);
            assert_eq!(d.degree(), degree, "{name} derivative {nu}");
            assert_eq!(d.extrapolate(), p.extrapolate(), "{name} derivative {nu}");
            assert_eq!(d.coefficients().len(), want.len(), "{name} derivative {nu}");
            for (i, (got, want)) in d.coefficients().iter().zip(&want).enumerate() {
                let what = format!("{name} derivative {nu} coefficient {i}");
                assert_close(*got, *want, DIRECT_F64, &what);
            }
            let values = nums(&entry["values"]);
            let what = format!("{name} derivative {nu}");
            check_column(&queries, &values, DIRECT_F64, &what, |x| d.eval(&x));
            check_column(&queries, &values, DIRECT_F64, &what, |x| {
                p.eval_derivative(x, nu)
            });
        }

        if let Some(antiderivatives) = case["antiderivatives"].as_array() {
            for entry in antiderivatives {
                let nu = usize::try_from(entry["nu"].as_u64().expect("nu")).expect("usize");
                let primitive = p.antiderivative(nu).expect("finite");
                let (want, degree) = interval_major(&entry["coefficients"]);
                assert_eq!(primitive.degree(), degree, "{name} antiderivative {nu}");
                for (i, (got, want)) in primitive.coefficients().iter().zip(&want).enumerate() {
                    let what = format!("{name} antiderivative {nu} coefficient {i}");
                    assert_close(*got, *want, DIRECT_F64, &what);
                }
                check_column(
                    &queries,
                    &nums(&entry["values"]),
                    DIRECT_F64,
                    &format!("{name} antiderivative {nu}"),
                    |x| primitive.eval(&x),
                );
            }
        } else {
            periodic += 1;
        }
        assert!(check_integrals(&case["integrals"], name, |a, b| p.integrate(a, b)) >= 4);
        cases += 1;
    }
    assert_eq!((cases, periodic), (4, 1));
}

#[test]
fn derivatives_and_antiderivatives_of_b_splines_match_scipys_knots_and_control_points() {
    let file = load::<Json>("calculus").expect("golden file calculus");
    let (mut cases, mut periodic) = (0, 0);
    for case in &file.cases {
        let name = case["name"].as_str().expect("name");
        if !name.starts_with("calculus/bspline/") {
            continue;
        }
        let spline = bspline_of(case).expect("valid");
        let queries = nums(&case["queries"]);
        // Checks that `got` is SciPy's object: the same knots, bit for bit, and
        // control points within the tolerance.
        let same = |got: &BSpline<f64>, entry: &Json, what: &str| {
            assert_eq!(got.degree(), degree_of(entry), "{what}");
            let knots = nums(&entry["knots"]);
            assert_eq!(got.knots().len(), knots.len(), "{what}");
            for (g, w) in got.knots().iter().zip(&knots) {
                assert_eq!(g.to_bits(), w.to_bits(), "{what}: knot");
            }
            let control_points = nums(&entry["control_points"]);
            assert_eq!(got.control_points().len(), control_points.len(), "{what}");
            for (i, (g, w)) in got.control_points().iter().zip(&control_points).enumerate() {
                assert_close(*g, *w, DIRECT_F64, &format!("{what} control point {i}"));
            }
        };

        let derivatives = case["derivatives"].as_array().expect("derivatives");
        assert_eq!(derivatives.len(), spline.degree() + 1, "{name}");
        for entry in derivatives {
            let nu = usize::try_from(entry["nu"].as_u64().expect("nu")).expect("usize");
            let what = format!("{name} derivative {nu}");
            let d = spline.derivative(nu).expect("finite");
            assert_eq!(d.extrapolate(), spline.extrapolate(), "{what}");
            if nu <= spline.degree() {
                same(&d, entry, &what);
            } else {
                // SciPy has no object for this order, only the zeros.
                assert_eq!(d.degree(), 0, "{what}");
            }
            let values = nums(&entry["values"]);
            check_column(&queries, &values, DIRECT_F64, &what, |x| d.eval(&x));
            check_column(&queries, &values, DIRECT_F64, &what, |x| {
                spline.eval_derivative(x, nu)
            });
        }

        if let Some(antiderivatives) = case["antiderivatives"].as_array() {
            for entry in antiderivatives {
                let nu = usize::try_from(entry["nu"].as_u64().expect("nu")).expect("usize");
                let what = format!("{name} antiderivative {nu}");
                let primitive = spline.antiderivative(nu).expect("finite");
                same(&primitive, entry, &what);
                check_column(&queries, &nums(&entry["values"]), DIRECT_F64, &what, |x| {
                    primitive.eval(&x)
                });
            }
        } else {
            periodic += 1;
        }
        assert!(check_integrals(&case["integrals"], name, |a, b| spline.integrate(a, b)) >= 4);
        cases += 1;
    }
    assert_eq!((cases, periodic), (7, 1));
}

// --- Where the two types agree and SciPy's do not --------------------------------

/// x^2 on [0, 1], then 1 + 3(x - 1) on [1, 3]: degree 2, a kink at 1.
fn bent() -> Result<PPoly<f64>, E> {
    PPoly::new(&[0.0, 1.0, 3.0], &[1.0, 0.0, 0.0, 0.0, 3.0, 1.0], 2)
}

#[test]
fn an_order_above_the_degree_gives_zero_on_both_types() {
    let poly = bent().expect("valid");
    let spline = BSpline::from_ppoly(&poly).expect("degree 2");

    let zero = poly.derivative(3).expect("zeros are finite");
    assert_eq!(zero.degree(), 0);
    assert_eq!(zero.breakpoints().len(), 3, "the break points stay");
    assert!(zero
        .coefficients()
        .iter()
        .all(|c| c.to_bits() == 0.0_f64.to_bits()));

    // SciPy raises for this one.
    let zero = spline.derivative(3).expect("zeros are finite");
    assert_eq!(zero.degree(), 0);
    assert_eq!(zero.knots().len(), 2);
    assert_eq!(zero.knots()[0].to_bits(), 0.0_f64.to_bits());
    assert_eq!(zero.knots()[1].to_bits(), 3.0_f64.to_bits());
    assert_eq!(zero.control_points().len(), 1);

    for nu in [3, 4, 100, usize::MAX] {
        for x in [0.0, 0.5, 1.0, 3.0, 7.0] {
            assert_eq!(
                bits(poly.eval_derivative(x, nu)),
                bits(Some(0.0)),
                "{nu} {x}"
            );
            assert_eq!(
                bits(spline.eval_derivative(x, nu)),
                bits(Some(0.0)),
                "{nu} {x}"
            );
            assert_eq!(
                bits(poly.derivative(nu).expect("zeros").eval(&x)),
                bits(Some(0.0))
            );
            assert_eq!(
                bits(spline.derivative(nu).expect("zeros").eval(&x)),
                bits(Some(0.0))
            );
        }
        // Zero where there is a value, and no value where there is none.
        let forbid = poly.clone().with_extrapolate(Extrapolate::Forbid);
        assert_eq!(forbid.eval_derivative(7.0, nu), None);
        assert_eq!(forbid.eval_derivative(f64::NAN, nu), None);
        let forbid = spline.clone().with_extrapolate(Extrapolate::Forbid);
        assert_eq!(forbid.eval_derivative(7.0, nu), None);
        assert_eq!(forbid.eval_derivative(f64::NAN, nu), None);
    }
}

#[test]
fn order_zero_gives_a_copy_with_the_same_mode() {
    for mode in MODES {
        let poly = bent().expect("valid").with_extrapolate(mode);
        let spline = BSpline::from_ppoly(&poly).expect("degree 2");
        for copy in [
            poly.derivative(0).expect("copy"),
            poly.antiderivative(0).expect("copy"),
        ] {
            assert_eq!(copy.extrapolate(), mode);
            assert_eq!(copy.degree(), 2);
            assert_eq!(bits(copy.eval(&7.5)), bits(poly.eval(&7.5)));
        }
        for copy in [
            spline.derivative(0).expect("copy"),
            spline.antiderivative(0).expect("copy"),
        ] {
            assert_eq!(copy.extrapolate(), mode);
            assert_eq!(copy.knots().len(), spline.knots().len());
            assert_eq!(bits(copy.eval(&7.5)), bits(spline.eval(&7.5)));
        }
    }
}

#[test]
fn a_b_spline_with_full_multiplicity_inside_is_differentiated_piece_by_piece() {
    // `from_ppoly` repeats every break point `degree + 1` times. SciPy refuses to
    // differentiate such a spline.
    let poly = bent().expect("valid");
    let spline = BSpline::from_ppoly(&poly).expect("degree 2");
    assert_eq!(spline.knots().len(), 9);

    let slope = spline.derivative(1).expect("finite");
    assert_eq!(slope.degree(), 1);
    assert_eq!(slope.knots().len(), 7);
    // 2x on [0, 1), then 3 on [1, 3].
    for (x, want) in [(0.0, 0.0), (0.5, 1.0), (1.0, 3.0), (2.0, 3.0), (3.0, 3.0)] {
        assert_close(slope.eval(&x).expect("inside"), want, DIRECT_F64, "slope");
        assert_close(
            spline.eval_derivative(x, 1).expect("inside"),
            want,
            DIRECT_F64,
            "eval_derivative",
        );
    }
    // 2 on [0, 1), then 0.
    let curvature = spline.derivative(2).expect("finite");
    for (x, want) in [(0.0, 2.0), (0.5, 2.0), (1.0, 0.0), (3.0, 0.0)] {
        assert_close(
            curvature.eval(&x).expect("inside"),
            want,
            DIRECT_F64,
            "curvature",
        );
    }
}

#[test]
fn the_antiderivative_is_zero_at_the_left_end_of_the_domain_whatever_the_knots() {
    // Degree 3 on the knots 0..=9 with every control point 1: the spline is 1 on
    // its domain [3, 6]. SciPy's antiderivative is 1.5 at 3.
    let knots: Vec<f64> = (0..10).map(f64::from).collect();
    let spline = BSpline::new(&knots, &[1.0; 6], 3).expect("valid");
    let primitive = spline.antiderivative(1).expect("finite");
    assert_eq!(primitive.degree(), 4);
    assert_eq!(primitive.knots().len(), 12);
    assert_close(
        primitive.eval(&3.0).expect("inside"),
        0.0,
        DIRECT_F64,
        "left end",
    );
    assert_close(
        primitive.eval(&4.5).expect("inside"),
        1.5,
        DIRECT_F64,
        "middle",
    );
    assert_close(
        primitive.eval(&6.0).expect("inside"),
        3.0,
        DIRECT_F64,
        "right end",
    );
    assert_close(
        spline.integrate(3.0, 6.0).expect("inside"),
        3.0,
        DIRECT_F64,
        "integral",
    );

    // Order 2: x^2 / 2 from the left end, zero there with its first derivative.
    let second = spline.antiderivative(2).expect("finite");
    for x in [3.0, 3.5, 5.0, 6.0] {
        let want = (x - 3.0) * (x - 3.0) / 2.0;
        assert_close(
            second.eval(&x).expect("inside"),
            want,
            DIRECT_F64,
            "order 2",
        );
    }

    // The same for a piecewise polynomial that does not start at zero.
    let poly = PPoly::new(&[2.0, 3.0, 5.0], &[4.0, 6.0], 0).expect("valid");
    let primitive = poly.antiderivative(1).expect("finite");
    assert_eq!(bits(primitive.eval(&2.0)), bits(Some(0.0)));
    assert_eq!(bits(primitive.eval(&3.0)), bits(Some(4.0)));
    assert_eq!(bits(primitive.eval(&5.0)), bits(Some(16.0)));
}

#[test]
fn derivatives_and_antiderivatives_take_the_documented_mode() {
    use Extrapolate::{Clamp, Extend, Forbid, Periodic};
    let table = [
        (Extend, Extend, Extend),
        (Forbid, Forbid, Forbid),
        (Periodic, Periodic, Forbid),
        (Clamp, Forbid, Forbid),
    ];
    for (mode, of_derivative, of_antiderivative) in table {
        let poly = bent().expect("valid").with_extrapolate(mode);
        let spline = BSpline::from_ppoly(&poly).expect("degree 2");
        // An order above the degree takes the same path of the table.
        for nu in [1, 2, 3] {
            assert_eq!(
                poly.derivative(nu).expect("finite").extrapolate(),
                of_derivative
            );
            assert_eq!(
                spline.derivative(nu).expect("finite").extrapolate(),
                of_derivative
            );
        }
        for nu in [1, 2] {
            assert_eq!(
                poly.antiderivative(nu).expect("finite").extrapolate(),
                of_antiderivative
            );
            assert_eq!(
                spline.antiderivative(nu).expect("finite").extrapolate(),
                of_antiderivative
            );
        }
    }
}

#[test]
fn an_antiderivative_too_large_to_hold_is_an_error() {
    let line = BSpline::new(&[0.0, 0.0, 1.0, 1.0], &[0.0, 2.0], 1).expect("valid");
    assert_eq!(BSPLINE_MAX_DEGREE, 15);
    assert_eq!(line.antiderivative(14).expect("degree 15").degree(), 15);
    for nu in [15, 16, usize::MAX] {
        assert_eq!(
            line.antiderivative(nu).err(),
            Some(E::InvalidParameter {
                name: "nu",
                reason: ParamReason::TooLarge { max: 14 },
            })
        );
    }
    // At the limit already: order 0 still gives the copy.
    let top = line.antiderivative(14).expect("degree 15");
    assert_eq!(top.antiderivative(0).expect("copy").degree(), 15);
    assert_eq!(
        top.antiderivative(1).err(),
        Some(E::InvalidParameter {
            name: "nu",
            reason: ParamReason::TooLarge { max: 0 },
        })
    );

    // A piecewise polynomial has no degree limit, only the size of a slice.
    let poly = PPoly::new(&[0.0, 1.0, 2.0], &[1.0, 0.0, 2.0, 1.0], 1).expect("valid");
    assert_eq!(poly.antiderivative(20).expect("finite").degree(), 21);
    // Two intervals of `f64`: a `Vec` holds `isize::MAX / 8` of them, so an
    // order of half that at most, and this polynomial has order 2 already.
    let max = isize::MAX.unsigned_abs() / 8 / 2 - 2;
    let too_large = E::InvalidParameter {
        name: "nu",
        reason: ParamReason::TooLarge { max },
    };
    // Too many to count, too many to count in bytes, and one over the limit.
    for nu in [
        usize::MAX,
        usize::MAX - 1,
        usize::MAX / 2,
        usize::MAX / 4,
        usize::MAX / 8,
        usize::MAX / 16,
        max + 1,
    ] {
        assert_eq!(poly.antiderivative(nu).err(), Some(too_large), "nu = {nu}");
    }
    // The limit itself is not tried: a `Vec` could hold that many, so the call
    // would go on to ask for the memory, and running out of memory is left to
    // the allocator (`.cursor/rules/rust.mdc`).

    // The limit counts bytes: with three columns it is a third.
    let wide: PPoly<f64, [f64; 3]> = PPoly::new(&[0.0, 1.0], &[[1.0; 3]], 0).expect("valid");
    assert_eq!(
        wide.antiderivative(usize::MAX / 16).err(),
        Some(E::InvalidParameter {
            name: "nu",
            reason: ParamReason::TooLarge {
                max: isize::MAX.unsigned_abs() / 24 - 1
            },
        })
    );
    // A value of no columns takes no memory; only the count limits it.
    let empty: PPoly<f64, [f64; 0]> = PPoly::new(&[0.0, 1.0], &[[]], 0).expect("valid");
    assert_eq!(empty.antiderivative(3).expect("no columns").degree(), 3);
    assert!(empty.antiderivative(usize::MAX).is_err());
}

#[test]
fn a_coefficient_that_overflows_is_an_error() {
    // 5e307 x^3: the first derivative has 1.5e308, the second 3e308.
    let poly = PPoly::new(&[0.0, 1.0], &[5e307, 0.0, 0.0, 0.0], 3).expect("valid");
    assert!(poly.derivative(1).is_ok());
    assert_eq!(poly.derivative(2).err(), Some(E::NonFinite { index: 0 }));
    // 1e308 on an interval 100 wide: the integral over it is the constant of the
    // next piece.
    let poly = PPoly::new(&[0.0, 100.0, 101.0], &[1e308, 0.0], 0).expect("valid");
    assert_eq!(
        poly.antiderivative(1).err(),
        Some(E::NonFinite { index: 3 })
    );

    // Control points 2e308 apart.
    let spline = BSpline::new(&[0.0, 0.0, 1.0, 1.0], &[-1e308, 1e308], 1).expect("valid");
    assert_eq!(spline.derivative(1).err(), Some(E::NonFinite { index: 0 }));
    let spline = BSpline::new(&[0.0, 0.0, 100.0, 100.0], &[1e308, 1e308], 1).expect("valid");
    assert_eq!(
        spline.antiderivative(1).err(),
        Some(E::NonFinite { index: 1 })
    );
}

// --- eval_derivative and integrate outside the domain ----------------------------

#[test]
fn eval_derivative_follows_each_mode_outside_the_domain() {
    // x^2 on [0, 1], then 1 + 3(x - 1) on [1, 3].
    let poly = bent().expect("valid");
    let spline = BSpline::from_ppoly(&poly).expect("degree 2");
    for mode in MODES {
        let poly = poly.clone().with_extrapolate(mode);
        let spline = spline.clone().with_extrapolate(mode);
        let both = |x: f64, nu: usize| {
            let (p, s) = (poly.eval_derivative(x, nu), spline.eval_derivative(x, nu));
            match (p, s) {
                (Some(p), Some(s)) => assert_close(s, p, DIRECT_F64, &format!("{mode:?} {x}")),
                (None, None) => {}
                other => panic!("{mode:?} at {x}: the two types disagree: {other:?}"),
            }
            p
        };
        // Inside, the mode does not matter. At the kink the piece on the right
        // decides, and at the right end the piece on the left. A period holds
        // its left end only, so there the right end is the left end again.
        assert_eq!(bits(both(0.5, 1)), bits(Some(1.0)), "{mode:?}");
        assert_eq!(bits(both(1.0, 1)), bits(Some(3.0)), "{mode:?}");
        let at_right_end = if mode == Extrapolate::Periodic {
            0.0
        } else {
            3.0
        };
        assert_eq!(bits(both(3.0, 1)), bits(Some(at_right_end)), "{mode:?}");
        assert_eq!(bits(both(0.5, 2)), bits(Some(2.0)), "{mode:?}");
        assert_eq!(bits(both(1.0, 2)), bits(Some(0.0)), "{mode:?}");
        assert_eq!(both(f64::NAN, 1), None, "{mode:?}");

        let (below, above) = (both(-1.0, 1), both(4.5, 1));
        match mode {
            Extrapolate::Forbid => assert_eq!((below, above), (None, None)),
            // Constant outside: the slope is zero, not the slope at the end.
            Extrapolate::Clamp => {
                assert_eq!(bits(below), bits(Some(0.0)));
                assert_eq!(bits(above), bits(Some(0.0)));
                assert_eq!(bits(both(-1.0, 2)), bits(Some(0.0)));
            }
            // 2x continued to -1, and the line continued.
            Extrapolate::Extend => {
                assert_eq!(bits(below), bits(Some(-2.0)));
                assert_eq!(bits(above), bits(Some(3.0)));
            }
            // The period is 3: -1 is 2, and 4.5 is 1.5.
            _ => {
                assert_eq!(bits(below), bits(Some(3.0)));
                assert_eq!(bits(above), bits(Some(3.0)));
                assert_eq!(bits(both(3.5, 1)), bits(Some(1.0)));
            }
        }
    }
}

#[test]
fn integrate_follows_each_mode_outside_the_domain() {
    // x^2 on [0, 1], then 1 + 3(x - 1) on [1, 3]. Its integral over [0, 1] is
    // 1/3, over [1, 3] it is 8, and the end values are 0 and 7.
    let poly = bent().expect("valid");
    let spline = BSpline::from_ppoly(&poly).expect("degree 2");
    let whole = 1.0 / 3.0 + 8.0;
    for mode in MODES {
        let poly = poly.clone().with_extrapolate(mode);
        let spline = spline.clone().with_extrapolate(mode);
        let both = |a: f64, b: f64| {
            let (p, s) = (poly.integrate(a, b), spline.integrate(a, b));
            match (p, s) {
                (Some(p), Some(s)) => assert_close(s, p, DIRECT_F64, &format!("{mode:?} {a} {b}")),
                (None, None) => {}
                other => panic!("{mode:?} from {a} to {b}: the two types disagree: {other:?}"),
            }
            // Reversed limits change the sign and nothing else.
            if a < b {
                assert_eq!(bits(poly.integrate(b, a)), reversed(p), "{mode:?} {a} {b}");
                assert_eq!(
                    bits(spline.integrate(b, a)),
                    reversed(s),
                    "{mode:?} {a} {b}"
                );
            }
            p
        };
        let close = |got: Option<f64>, want: f64, what: &str| {
            assert_close(
                got.expect("a value"),
                want,
                DIRECT_F64,
                &format!("{mode:?} {what}"),
            );
        };

        // Inside, the mode does not matter.
        close(both(0.0, 3.0), whole, "whole");
        close(both(0.5, 2.0), 7.0 / 24.0 + 2.5, "across the kink");
        close(both(1.25, 1.25), 0.0, "equal limits");
        close(both(0.0, 0.0), 0.0, "equal limits at the end");
        for (a, b) in [
            (f64::NAN, 1.0),
            (1.0, f64::NAN),
            (f64::INFINITY, 1.0),
            (0.0, f64::NEG_INFINITY),
        ] {
            assert_eq!(poly.integrate(a, b), None, "{mode:?} {a} {b}");
            assert_eq!(spline.integrate(a, b), None, "{mode:?} {a} {b}");
        }

        match mode {
            Extrapolate::Forbid => {
                assert_eq!(both(-1.0, 2.0), None);
                assert_eq!(both(1.0, 3.5), None);
                assert_eq!(both(-2.0, -1.0), None);
                // The ends themselves are inside.
                close(both(0.0, 3.0), whole, "the ends");
            }
            // 0 on the left of the domain, 7 on the right.
            Extrapolate::Clamp => {
                close(both(-2.0, 3.0), whole, "the left adds nothing");
                close(both(0.0, 5.0), whole + 14.0, "two wide on the right");
                close(both(-2.0, 5.0), whole + 14.0, "both sides");
                close(both(4.0, 6.0), 14.0, "all on the right");
                close(both(-5.0, -1.0), 0.0, "all on the left");
                close(both(2.0, 4.0), 5.5 + 7.0, "from inside");
            }
            // x^2 on the left, the line on the right.
            Extrapolate::Extend => {
                close(both(-1.0, 0.0), 1.0 / 3.0, "x^2 continued");
                close(both(-1.0, 3.0), whole + 1.0 / 3.0, "from the left");
                // The line is 3x - 2: from 3 to 5 it gives 20.
                close(both(3.0, 5.0), 20.0, "the line continued");
                close(both(4.0, 5.0), 11.5, "all on the right");
                close(both(-2.0, -1.0), 7.0 / 3.0, "all on the left");
            }
            // The period is 3.
            _ => {
                close(both(0.0, 9.0), 3.0 * whole, "three periods");
                close(both(3.0, 6.0), whole, "the second period");
                close(both(-3.0, 0.0), whole, "the period before");
                close(both(3.5, 5.0), 7.0 / 24.0 + 2.5, "0.5 to 2, one period on");
                // From 2 to 3, then from 0 to 1.
                close(both(2.0, 4.0), 5.5 + 1.0 / 3.0, "across the end");
                close(both(-1.0, 1.0), 5.5 + 1.0 / 3.0, "across the start");
                close(
                    both(2.0, 10.0),
                    5.5 + 2.0 * whole + 1.0 / 3.0,
                    "two periods and a part",
                );
            }
        }
    }
}

#[test]
fn array_values_and_f32_go_through_the_same_code() {
    // Two columns: x^2 and 2 - x on [0, 2].
    let poly: PPoly<f32, [f32; 2]> =
        PPoly::new(&[0.0, 2.0], &[[1.0, 0.0], [0.0, -1.0], [0.0, 2.0]], 2).expect("valid");
    let spline = BSpline::from_ppoly(&poly).expect("degree 2");
    let same = |got: Option<[f32; 2]>, want: [f32; 2], what: &str| {
        let got = got.expect("a value");
        for (g, w) in got.iter().zip(want) {
            assert!((g - w).abs() <= 1e-5, "{what}: {got:?} for {want:?}");
        }
    };
    same(poly.eval_derivative(1.0, 1), [2.0, -1.0], "PPoly slope");
    same(spline.eval_derivative(1.0, 1), [2.0, -1.0], "BSpline slope");
    same(poly.integrate(0.0, 2.0), [8.0 / 3.0, 2.0], "PPoly integral");
    same(
        spline.integrate(0.0, 2.0),
        [8.0 / 3.0, 2.0],
        "BSpline integral",
    );
    same(
        poly.derivative(1).expect("finite").eval(&1.0),
        [2.0, -1.0],
        "PPoly derivative",
    );
    same(
        spline.derivative(1).expect("finite").eval(&1.0),
        [2.0, -1.0],
        "BSpline derivative",
    );
    same(
        poly.antiderivative(1).expect("finite").eval(&2.0),
        [8.0 / 3.0, 2.0],
        "PPoly antiderivative",
    );
    same(
        spline.antiderivative(1).expect("finite").eval(&2.0),
        [8.0 / 3.0, 2.0],
        "BSpline antiderivative",
    );
}

#[test]
fn a_spline_of_the_largest_degree_can_be_integrated() {
    // One piece of degree 15 on [0, 2] in Bernstein form. With every control
    // point 1 it is the constant 1; with the last one only it is (x / 2)^15.
    let k = BSPLINE_MAX_DEGREE;
    let knots: Vec<f64> = [vec![0.0; k + 1], vec![2.0; k + 1]].concat();
    let one = BSpline::new(&knots, &vec![1.0; k + 1], k).expect("valid");
    assert_close(
        one.integrate(0.5, 1.75).expect("inside"),
        1.25,
        DIRECT_F64,
        "constant",
    );
    let mut last = vec![0.0; k + 1];
    last[k] = 1.0;
    let power = BSpline::new(&knots, &last, k).expect("valid");
    assert_close(
        power.integrate(0.0, 2.0).expect("inside"),
        2.0 / 16.0,
        DIRECT_F64,
        "(x / 2)^15",
    );
    assert_close(
        power.integrate(1.0, 2.0).expect("inside"),
        (2.0 - 2.0 * 0.5_f64.powi(16)) / 16.0,
        DIRECT_F64,
        "(x / 2)^15 from 1",
    );
    // Its derivatives of every order exist, and order 15 is the constant 15! / 2^15.
    let top = power.derivative(k).expect("finite");
    assert_eq!(top.degree(), 0);
    let want = (1..=15).map(f64::from).product::<f64>() / 2.0_f64.powi(15);
    assert_close(
        top.eval(&1.0).expect("inside"),
        want,
        DIRECT_F64,
        "order 15",
    );
    assert_close(
        power.eval_derivative(1.0, k).expect("inside"),
        want,
        DIRECT_F64,
        "eval_derivative 15",
    );
}

// --- Properties ----------------------------------------------------------------

/// Break points, coefficients, and a degree of a random piecewise polynomial.
fn pieces() -> impl Strategy<Value = (Vec<f64>, Vec<f64>, usize)> {
    (1_usize..8, 0_usize..5)
        .prop_flat_map(|(intervals, degree)| {
            (
                -5.0_f64..5.0,
                prop::collection::vec(0.25_f64..3.0, intervals),
                prop::collection::vec(-10.0_f64..10.0, intervals * (degree + 1)),
                Just(degree),
            )
        })
        .prop_map(|(start, widths, coefficients, degree)| {
            let mut breakpoints = vec![start];
            for w in widths {
                breakpoints.push(breakpoints[breakpoints.len() - 1] + w);
            }
            (breakpoints, coefficients, degree)
        })
}

/// Knots, control points, and a degree of a random B-spline on [0, 3]. The inner
/// knots come from a coarse grid, so some repeat, up to `degree + 1` times. With
/// `clamped` false the end knots are spread out instead of repeated.
fn splines() -> impl Strategy<Value = (Vec<f64>, Vec<f64>, usize)> {
    (0_usize..5, 0_usize..7, any::<bool>())
        .prop_flat_map(|(degree, inner, clamped)| {
            (
                prop::collection::vec(1_u8..12, inner),
                prop::collection::vec(-10.0_f64..10.0, degree + 1 + inner),
                Just(degree),
                Just(clamped),
            )
        })
        .prop_map(|(mut inner, control_points, degree, clamped)| {
            inner.sort_unstable();
            let mut kept: Vec<u8> = Vec::new();
            for knot in inner {
                if kept.iter().filter(|k| **k == knot).count() <= degree {
                    kept.push(knot);
                }
            }
            let n = degree + 1 + kept.len();
            let end = |i: usize| if clamped { 0.0 } else { 0.5 * i as f64 };
            let mut knots: Vec<f64> = (0..=degree).rev().map(|i| -end(i)).collect();
            knots.extend(kept.iter().map(|k| f64::from(*k) * 0.25));
            knots.extend((0..=degree).map(|i| 3.0 + end(i)));
            (knots, control_points[..n].to_vec(), degree)
        })
}

/// The derivative of order `nu` of one piece, written apart from the library:
/// term by term, with the factor each power gets.
fn piece_derivative(piece: &[f64], dx: f64, nu: usize) -> f64 {
    let degree = piece.len() - 1;
    piece
        .iter()
        .enumerate()
        .filter(|(j, _)| degree - j >= nu)
        .map(|(j, c)| {
            let power = degree - j;
            let factor: f64 = (0..nu).map(|r| (power - r) as f64).product();
            c * factor * dx.powi(i32::try_from(power - nu).unwrap_or(i32::MAX))
        })
        .sum()
}

/// The integral of one piece from `p` to `q`, both measured from its left break
/// point, term by term.
fn piece_integral(piece: &[f64], p: f64, q: f64) -> f64 {
    let degree = piece.len() - 1;
    piece
        .iter()
        .enumerate()
        .map(|(j, c)| {
            let n = i32::try_from(degree - j + 1).unwrap_or(i32::MAX);
            c * (q.powi(n) - p.powi(n)) / f64::from(n)
        })
        .sum()
}

/// A bound on the size of what the tests compare, for a tolerance.
fn scale(numbers: &[f64]) -> f64 {
    numbers.iter().fold(1.0_f64, |m, c| m.max(c.abs()))
}

/// `integrate`, `integrate` under a given mode, and `eval`, of either type: the
/// properties below run on both through these.
type Integral<'a> = &'a dyn Fn(f64, f64) -> Option<f64>;
type IntegralIn<'a> = &'a dyn Fn(Extrapolate, f64, f64) -> Option<f64>;
type Eval<'a> = &'a dyn Fn(f64) -> Option<f64>;

/// Whether two results agree: both absent, or both within `tolerance`.
fn agree(a: Option<f64>, b: Option<f64>, tolerance: f64) -> bool {
    match (a, b) {
        (Some(a), Some(b)) => (a - b).abs() <= tolerance * a.abs().max(b.abs()).max(1.0),
        (None, None) => true,
        _ => false,
    }
}

proptest! {
    #[test]
    fn order_zero_is_eval_bit_for_bit(
        (breakpoints, coefficients, degree) in pieces(),
        (knots, control_points, k) in splines(),
        xs in prop::collection::vec(-12.0_f64..12.0, 1..30),
    ) {
        for mode in MODES {
            let poly = PPoly::new(&breakpoints, &coefficients, degree).expect("valid").with_extrapolate(mode);
            let spline = BSpline::new(&knots, &control_points, k).expect("valid").with_extrapolate(mode);
            for x in xs.iter().copied().chain([f64::NAN, f64::INFINITY]) {
                prop_assert_eq!(bits(poly.eval_derivative(x, 0)), bits(poly.eval(&x)), "{:?} {}", mode, x);
                prop_assert_eq!(bits(spline.eval_derivative(x, 0)), bits(spline.eval(&x)), "{:?} {}", mode, x);
            }
        }
    }

    #[test]
    fn a_piecewise_polynomial_differentiates_and_integrates_term_by_term(
        (breakpoints, coefficients, degree) in pieces(),
        t in 0.0_f64..1.0,
        u in 0.0_f64..1.0,
        nu in 1_usize..7,
    ) {
        let poly = PPoly::new(&breakpoints, &coefficients, degree).expect("valid");
        let order = degree + 1;
        let (lo, hi) = (breakpoints[0], breakpoints[breakpoints.len() - 1]);
        let width = (hi - lo).max(1.0);
        let tolerance = 1e-10 * scale(&coefficients) * width.powi(6);

        // The derivative at a point of each interval, and at its left end.
        for (i, piece) in coefficients.chunks_exact(order).enumerate() {
            let (left, right) = (breakpoints[i], breakpoints[i + 1]);
            for x in [left, left + t * (right - left)] {
                if x < right {
                    let want = piece_derivative(piece, x - left, nu);
                    let got = poly.eval_derivative(x, nu).expect("inside");
                    prop_assert!((got - want).abs() <= tolerance, "interval {} x = {}: {} for {}", i, x, got, want);
                }
            }
        }

        // The integral between two points of the domain: the pieces between them.
        let (a, b) = (lo + t.min(u) * (hi - lo), lo + t.max(u) * (hi - lo));
        let mut want = 0.0;
        for (i, piece) in coefficients.chunks_exact(order).enumerate() {
            let (left, right) = (breakpoints[i], breakpoints[i + 1]);
            let (p, q) = (a.max(left), b.min(right));
            if p < q {
                want += piece_integral(piece, p - left, q - left);
            }
        }
        let got = poly.integrate(a, b).expect("inside");
        prop_assert!((got - want).abs() <= tolerance, "from {} to {}: {} for {}", a, b, got, want);
    }

    #[test]
    fn the_derivative_object_and_eval_derivative_agree(
        (breakpoints, coefficients, degree) in pieces(),
        (knots, control_points, k) in splines(),
        xs in prop::collection::vec(-12.0_f64..12.0, 1..30),
        nu in 1_usize..7,
    ) {
        // Extend, Forbid, and Periodic are inherited, so the two agree
        // everywhere. Clamp is not: there they agree inside the domain.
        for mode in MODES {
            let poly = PPoly::new(&breakpoints, &coefficients, degree).expect("valid").with_extrapolate(mode);
            let d = poly.derivative(nu).expect("finite");
            let (lo, hi) = (breakpoints[0], breakpoints[breakpoints.len() - 1]);
            for x in &xs {
                let direct = poly.eval_derivative(*x, nu);
                if mode == Extrapolate::Clamp && !(lo..=hi).contains(x) {
                    prop_assert_eq!(bits(direct), bits(Some(0.0)));
                    prop_assert_eq!(d.eval(x), None);
                } else {
                    prop_assert!(agree(direct, d.eval(x), 1e-9), "PPoly {:?} x = {}: {:?} and {:?}", mode, x, direct, d.eval(x));
                }
            }

            let spline = BSpline::new(&knots, &control_points, k).expect("valid").with_extrapolate(mode);
            let d = spline.derivative(nu).expect("finite");
            prop_assert_eq!(d.degree(), k.saturating_sub(nu));
            for x in &xs {
                let direct = spline.eval_derivative(*x, nu);
                if mode == Extrapolate::Clamp && !(0.0..=3.0).contains(x) {
                    prop_assert_eq!(bits(direct), bits(Some(0.0)));
                    prop_assert_eq!(d.eval(x), None);
                } else {
                    prop_assert!(agree(direct, d.eval(x), 1e-7), "BSpline {:?} x = {}: {:?} and {:?}", mode, x, direct, d.eval(x));
                }
            }
        }
    }

    #[test]
    fn the_derivative_of_the_antiderivative_is_the_function(
        (breakpoints, coefficients, degree) in pieces(),
        (knots, control_points, k) in splines(),
        xs in prop::collection::vec(-6.0_f64..9.0, 1..30),
        nu in 1_usize..4,
    ) {
        let poly = PPoly::new(&breakpoints, &coefficients, degree).expect("valid");
        let primitive = poly.antiderivative(nu).expect("finite");
        prop_assert_eq!(primitive.degree(), degree + nu);
        let back = primitive.derivative(nu).expect("finite");
        prop_assert_eq!(back.degree(), degree);
        for x in &xs {
            prop_assert!(agree(back.eval(x), poly.eval(x), 1e-9), "PPoly x = {}", x);
            prop_assert!(agree(primitive.eval_derivative(*x, nu), poly.eval(x), 1e-9), "PPoly x = {}", x);
        }
        // Zero at the left end, with the derivatives below the order.
        for lower in 0..nu {
            let at_left = primitive.eval_derivative(breakpoints[0], lower).expect("inside");
            prop_assert!(at_left.abs() <= 1e-9, "PPoly order {} at the left end: {}", lower, at_left);
        }

        let spline = BSpline::new(&knots, &control_points, k).expect("valid");
        let primitive = spline.antiderivative(nu).expect("finite");
        prop_assert_eq!(primitive.degree(), k + nu);
        prop_assert_eq!(primitive.knots().len(), knots.len() + 2 * nu);
        prop_assert_eq!(primitive.domain(), spline.domain());
        let back = primitive.derivative(nu).expect("finite");
        for x in &xs {
            prop_assert!(agree(back.eval(x), spline.eval(x), 1e-7), "BSpline x = {}: {:?} and {:?}", x, back.eval(x), spline.eval(x));
        }
        for lower in 0..nu {
            let at_left = primitive.eval_derivative(0.0, lower).expect("inside");
            prop_assert!(at_left.abs() <= 1e-9 * scale(&control_points), "BSpline order {} at the left end: {}", lower, at_left);
        }
    }

    #[test]
    fn the_integral_is_the_difference_of_the_antiderivative(
        (breakpoints, coefficients, degree) in pieces(),
        (knots, control_points, k) in splines(),
        a in -6.0_f64..9.0,
        b in -6.0_f64..9.0,
    ) {
        // With Extend the antiderivative extends too, so this holds outside the
        // domain as well.
        let poly = PPoly::new(&breakpoints, &coefficients, degree).expect("valid");
        let primitive = poly.antiderivative(1).expect("finite");
        let want = primitive.eval(&b).expect("extends") - primitive.eval(&a).expect("extends");
        let got = poly.integrate(a, b).expect("extends");
        let bound = primitive.eval(&a).expect("extends").abs().max(primitive.eval(&b).expect("extends").abs()).max(1.0);
        prop_assert!((got - want).abs() <= 1e-9 * bound, "PPoly from {} to {}: {} for {}", a, b, got, want);

        let spline = BSpline::new(&knots, &control_points, k).expect("valid");
        let primitive = spline.antiderivative(1).expect("finite");
        let want = primitive.eval(&b).expect("extends") - primitive.eval(&a).expect("extends");
        let got = spline.integrate(a, b).expect("extends");
        let bound = primitive.eval(&a).expect("extends").abs().max(primitive.eval(&b).expect("extends").abs()).max(1.0);
        prop_assert!((got - want).abs() <= 1e-9 * bound, "BSpline from {} to {}: {} for {}", a, b, got, want);
    }

    #[test]
    fn the_integral_adds_over_intervals_and_changes_sign_when_reversed(
        (breakpoints, coefficients, degree) in pieces(),
        (knots, control_points, k) in splines(),
        a in -6.0_f64..9.0,
        b in -6.0_f64..9.0,
        c in -6.0_f64..9.0,
    ) {
        for mode in MODES {
            let poly = PPoly::new(&breakpoints, &coefficients, degree).expect("valid").with_extrapolate(mode);
            let spline = BSpline::new(&knots, &control_points, k).expect("valid").with_extrapolate(mode);
            let both: [(&str, Integral); 2] = [
                ("PPoly", &|p, q| poly.integrate(p, q)),
                ("BSpline", &|p, q| spline.integrate(p, q)),
            ];
            for (name, integrate) in both {
                if a < b {
                    prop_assert_eq!(bits(integrate(b, a)), reversed(integrate(a, b)), "{} {:?}", name, mode);
                }
                prop_assert!(integrate(a, a).is_none_or(|v| v.abs() <= 1e-12), "{} {:?}: equal limits", name, mode);
                if let (Some(ab), Some(bc), Some(ac)) = (integrate(a, b), integrate(b, c), integrate(a, c)) {
                    let bound = ab.abs().max(bc.abs()).max(1.0);
                    prop_assert!((ab + bc - ac).abs() <= 1e-8 * bound, "{} {:?}: {} + {} and {}", name, mode, ab, bc, ac);
                } else {
                    prop_assert_eq!(mode, Extrapolate::Forbid, "{}: only Forbid has no value", name);
                }
            }
        }
    }

    #[test]
    fn clamp_and_periodic_integrals_are_built_from_the_integral_inside(
        (breakpoints, coefficients, degree) in pieces(),
        (knots, control_points, k) in splines(),
        a in -6.0_f64..9.0,
        length in 0.0_f64..12.0,
        periods in 1_u8..4,
    ) {
        let b = a + length;
        let poly = PPoly::new(&breakpoints, &coefficients, degree).expect("valid");
        let spline = BSpline::new(&knots, &control_points, k).expect("valid");
        let both: [(&str, (f64, f64), IntegralIn, Eval); 2] = [
            (
                "PPoly",
                (breakpoints[0], breakpoints[breakpoints.len() - 1]),
                &|mode, p, q| poly.clone().with_extrapolate(mode).integrate(p, q),
                &|x| poly.eval(&x),
            ),
            (
                "BSpline",
                (0.0, 3.0),
                &|mode, p, q| spline.clone().with_extrapolate(mode).integrate(p, q),
                &|x| spline.eval(&x),
            ),
        ];
        for (name, (lo, hi), integrate, eval) in both {
            let inside = |p: f64, q: f64| integrate(Extrapolate::Forbid, p, q).expect("inside the domain");
            let whole = inside(lo, hi);
            let bound = whole.abs().max(1.0) * f64::from(periods) * 8.0;

            // Clamp: the end values over the parts outside, and the part inside.
            let mut want = 0.0;
            if a < lo {
                want += eval(lo).expect("the end") * (b.min(lo) - a);
            }
            if b > hi {
                want += eval(hi).expect("the end") * (b - a.max(hi));
            }
            if a.max(lo) < b.min(hi) {
                want += inside(a.max(lo), b.min(hi));
            }
            let got = integrate(Extrapolate::Clamp, a, b).expect("clamped");
            let clamp_bound = bound.max(want.abs());
            prop_assert!((got - want).abs() <= 1e-9 * clamp_bound, "{} clamp from {} to {}: {} for {}", name, a, b, got, want);

            // Periodic: a whole number of periods from anywhere is that many times
            // the integral over the domain, and a shift by a period changes nothing.
            let period = hi - lo;
            let n = f64::from(periods);
            let got = integrate(Extrapolate::Periodic, a, a + n * period).expect("periodic");
            prop_assert!((got - n * whole).abs() <= 1e-9 * bound, "{} {} periods from {}: {} for {}", name, periods, a, got, n * whole);
            let here = integrate(Extrapolate::Periodic, a, b).expect("periodic");
            let there = integrate(Extrapolate::Periodic, a + n * period, b + n * period).expect("periodic");
            prop_assert!((here - there).abs() <= 1e-8 * bound, "{} from {} to {}: {} here, {} a period on", name, a, b, here, there);
            // Inside the domain it is the plain integral.
            let (p, q) = (lo + 0.25 * period, lo + 0.75 * period);
            prop_assert_eq!(bits(integrate(Extrapolate::Periodic, p, q)), bits(Some(inside(p, q))), "{}", name);
        }
    }

    #[test]
    fn both_types_give_the_same_derivatives_and_integrals(
        (knots, control_points, k) in splines(),
        xs in prop::collection::vec(-2.0_f64..5.0, 1..30),
        a in -2.0_f64..5.0,
        b in -2.0_f64..5.0,
        nu in 1_usize..6,
    ) {
        let tolerance = 1e-7 * scale(&control_points);
        for mode in MODES {
            let spline = BSpline::new(&knots, &control_points, k).expect("valid").with_extrapolate(mode);
            let poly = PPoly::from_bspline(&spline).expect("finite coefficients");
            for x in &xs {
                let (s, p) = (spline.eval_derivative(*x, nu), poly.eval_derivative(*x, nu));
                match (s, p) {
                    (Some(s), Some(p)) => prop_assert!((s - p).abs() <= tolerance * s.abs().max(1.0), "{:?} order {} at {}: {} and {}", mode, nu, x, s, p),
                    (None, None) => {}
                    other => prop_assert!(false, "{:?} at {}: {:?}", mode, x, other),
                }
            }
            let (s, p) = (spline.integrate(a, b), poly.integrate(a, b));
            match (s, p) {
                (Some(s), Some(p)) => prop_assert!((s - p).abs() <= tolerance * s.abs().max(1.0), "{:?} from {} to {}: {} and {}", mode, a, b, s, p),
                (None, None) => {}
                other => prop_assert!(false, "{:?} from {} to {}: {:?}", mode, a, b, other),
            }
        }
    }
}
