//! `BSpline` and the conversions with `PPoly`: the golden data, the input checks
//! through the constructor, repeated knots, the extrapolation modes, and the
//! properties that hold for any B-spline.

mod common;

use common::{assert_close, load, reals, Real, Tol, DIRECT_F64};
use interpx::{
    BSpline, Extrapolate, Hint, InterpolateError as E, Interpolator, PPoly, ParamReason,
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

/// The tolerance for a value that went through a change of basis. The Taylor
/// coefficients of a spline are differences of control points divided by knot
/// spacings, so a conversion loses a few digits that an evaluation does not.
const CONVERTED: Tol = Tol {
    rtol: 1e-9,
    atol: 1e-9,
};

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

fn degree_of(case: &Json) -> usize {
    usize::try_from(case["k"].as_u64().unwrap_or(u64::MAX)).unwrap_or(usize::MAX)
}

/// The mode a golden case was generated with (`extrapolate` in its options).
fn mode_of(case: &Json) -> Extrapolate {
    match &case["options"]["extrapolate"] {
        Json::Bool(false) => Extrapolate::Forbid,
        Json::String(s) if s == "periodic" => Extrapolate::Periodic,
        _ => Extrapolate::Extend,
    }
}

// --- Golden data -------------------------------------------------------------

#[test]
fn scipy_b_splines_evaluate_to_the_golden_values() {
    let file = load::<Json>("oned").expect("golden file oned");
    let (mut cases, mut degrees, mut modes) = (0, std::collections::BTreeSet::new(), 0_u8);
    for case in &file.cases {
        let name = case["name"].as_str().expect("name");
        let (Some(knots), Some(control_points), Some(queries), Some(expected)) = (
            case.get("knots"),
            case.get("control_points"),
            case.get("queries"),
            case.get("expected"),
        ) else {
            continue;
        };
        let mode = mode_of(case);
        let spline = BSpline::from_owned(nums(knots), nums(control_points), degree_of(case))
            .expect("SciPy's B-spline is valid")
            .with_extrapolate(mode);
        let (queries, expected) = (nums(queries), nums(expected));
        assert_eq!(queries.len(), expected.len(), "{name}");

        for (k, (x, want)) in queries.iter().zip(&expected).enumerate() {
            let what = format!("{name} query {k} (x = {x})");
            match spline.eval(x) {
                // SciPy returns NaN where it does not extrapolate.
                None => assert!(want.is_nan(), "{what}: None, SciPy has {want}"),
                Some(got) => assert_close(got, *want, DIRECT_F64, &what),
            }
        }

        let mut out = vec![0.0; queries.len()];
        assert_eq!(spline.eval_batch_into(&queries, &mut out, f64::NAN), Ok(()));
        for (k, (got, want)) in out.iter().zip(&expected).enumerate() {
            assert_close(*got, *want, DIRECT_F64, &format!("{name} batch {k}"));
        }

        // The piecewise polynomial of the same spline gives the same values.
        let poly = PPoly::from_bspline(&spline).expect("finite coefficients");
        assert_eq!(poly.extrapolate(), mode, "{name}");
        for (k, (x, want)) in queries.iter().zip(&expected).enumerate() {
            let what = format!("{name} as PPoly, query {k} (x = {x})");
            match poly.eval(x) {
                None => assert!(want.is_nan(), "{what}"),
                Some(got) => assert_close(got, *want, CONVERTED, &what),
            }
        }

        cases += 1;
        degrees.insert(spline.degree());
        modes |= match mode {
            Extrapolate::Forbid => 1,
            Extrapolate::Periodic => 2,
            _ => 4,
        };
    }
    // Degrees 1 to 5 on five data sets, three boundary conditions on two, one
    // without extrapolation, one periodic.
    assert_eq!(cases, 33, "cases with a B-spline and expected values");
    assert_eq!(degrees.into_iter().collect::<Vec<_>>(), [1, 2, 3, 4, 5]);
    assert_eq!(modes, 7, "extend, forbid, and periodic");
}

/// One curve case with `D` outputs per point.
// Marked as test code so that `expect` is allowed here as in a `#[test]` function
// (clippy.toml allows it in tests only; this helper is generic, so it cannot be one).
#[cfg(test)]
fn check_curve<const D: usize>(case: &Json) {
    let name = case["name"].as_str().expect("name");
    let to_array = |r: &Vec<f64>| -> [f64; D] { r.as_slice().try_into().expect("D columns") };
    let control_points: Vec<[f64; D]> =
        rows(&case["control_points"]).iter().map(to_array).collect();
    let expected: Vec<[f64; D]> = rows(&case["expected"]).iter().map(to_array).collect();
    let queries = nums(&case["queries"]);
    let spline = BSpline::from_owned(nums(&case["knots"]), control_points, degree_of(case))
        .expect("SciPy's B-spline is valid")
        .with_extrapolate(mode_of(case));

    let mut out = vec![[0.0; D]; queries.len()];
    assert_eq!(
        spline.eval_batch_into(&queries, &mut out, [f64::NAN; D]),
        Ok(())
    );
    for (k, (x, want)) in queries.iter().zip(&expected).enumerate() {
        let got = spline.eval(x).expect("the curve cases extrapolate");
        for (column, (g, w)) in got.iter().zip(want).enumerate() {
            assert_close(
                *g,
                *w,
                DIRECT_F64,
                &format!("{name} query {k} column {column}"),
            );
            assert_eq!(g.to_bits(), out[k][column].to_bits(), "{name}: batch");
        }
    }
}

#[test]
fn scipy_curves_evaluate_to_the_golden_values_column_by_column() {
    let file = load::<Json>("curve").expect("golden file curve");
    let (mut two, mut three) = (0, 0);
    for case in &file.cases {
        match rows(&case["control_points"]).first().map_or(0, Vec::len) {
            2 => {
                check_curve::<2>(case);
                two += 1;
            }
            3 => {
                check_curve::<3>(case);
                three += 1;
            }
            other => panic!("a curve with {other} columns"),
        }
    }
    assert!(two >= 4 && three >= 4, "{two} curves in 2D, {three} in 3D");
}

#[test]
fn from_bspline_matches_scipy_on_the_intervals_of_positive_width() {
    let file = load::<Json>("oned").expect("golden file oned");
    let mut cases = 0;
    for case in &file.cases {
        let name = case["name"].as_str().expect("name");
        if case["reference"] != "scipy.interpolate.PPoly.from_spline" {
            continue;
        }
        let degree = usize::try_from(case["options"]["k"].as_u64().expect("k")).expect("usize");
        let spline =
            BSpline::from_owned(nums(&case["knots"]), nums(&case["control_points"]), degree)
                .expect("valid");
        let poly = PPoly::from_bspline(&spline).expect("finite coefficients");

        // SciPy keeps every knot as a break point: `c[power][interval]`, with an
        // interval of zero width for each repeated knot.
        let scipy_breaks = nums(&case["breakpoints"]);
        let scipy_rows = rows(&case["coefficients"]);
        let kept: Vec<usize> = (0..scipy_breaks.len() - 1)
            .filter(|i| scipy_breaks[*i] < scipy_breaks[*i + 1])
            .collect();
        assert!(
            kept.len() < scipy_breaks.len() - 1,
            "{name}: SciPy has no empty interval"
        );

        assert_eq!(poly.breakpoints().len(), kept.len() + 1, "{name}");
        for (ours, i) in poly.breakpoints().iter().zip(&kept) {
            assert_eq!(
                ours.to_bits(),
                scipy_breaks[*i].to_bits(),
                "{name}: break point"
            );
        }
        let order = degree + 1;
        for (piece, i) in poly.coefficients().chunks_exact(order).zip(&kept) {
            for (power, got) in piece.iter().enumerate() {
                let want = scipy_rows[power][*i];
                assert_close(
                    *got,
                    want,
                    CONVERTED,
                    &format!("{name} interval {i} row {power}"),
                );
            }
        }
        cases += 1;
    }
    assert_eq!(cases, 2);
}

#[test]
fn from_ppoly_keeps_the_values_of_every_scipy_piecewise_polynomial() {
    let all_cases = common::load_cases(&common::PIECEWISE).expect("golden files");
    let mut cases = 0;
    for case in &all_cases {
        let name = case["name"].as_str().expect("name");
        let (Some(breakpoints), Some(coefficients), Some(queries), Some(expected)) = (
            case.get("breakpoints"),
            case.get("coefficients"),
            case.get("queries"),
            case.get("expected"),
        ) else {
            continue;
        };
        let scipy_rows = rows(coefficients);
        let intervals = scipy_rows.first().map_or(0, Vec::len);
        let flat: Vec<f64> = (0..intervals)
            .flat_map(|i| scipy_rows.iter().map(move |row| row[i]))
            .collect();
        let mode = match &case["options"]["extrapolate"] {
            Json::Bool(true) => Extrapolate::Extend,
            Json::String(s) if s == "periodic" => Extrapolate::Periodic,
            _ => Extrapolate::Forbid,
        };
        let poly = PPoly::from_owned(nums(breakpoints), flat, scipy_rows.len() - 1)
            .expect("valid")
            .with_extrapolate(mode);
        let spline = BSpline::from_ppoly(&poly).expect("degree 3");
        assert_eq!(spline.extrapolate(), mode, "{name}");
        assert_eq!(spline.degree(), poly.degree(), "{name}");
        assert_eq!(spline.knots().len(), poly.breakpoints().len() * 4, "{name}");

        // Where an end piece is very short, its Bézier form loses digits far
        // outside it: the `uneven` data has an end piece of width 0.001, and a
        // query 700 widths away is off by 4e-8. Those cases are compared inside
        // the data only (`BSpline::from_ppoly` documents this).
        let (lo, hi) = (
            poly.breakpoints()[0],
            poly.breakpoints()[poly.breakpoints().len() - 1],
        );
        let inside_only = name.contains("uneven");
        for (k, (x, want)) in nums(queries).iter().zip(nums(expected)).enumerate() {
            if inside_only && !(lo..=hi).contains(x) {
                continue;
            }
            let what = format!("{name} as BSpline, query {k} (x = {x})");
            match spline.eval(x) {
                None => assert!(want.is_nan(), "{what}"),
                Some(got) => assert_close(got, want, CONVERTED, &what),
            }
        }
        cases += 1;
    }
    assert_eq!(cases, 84, "the cases tests/ppoly.rs evaluates");
}

// --- Input checks through the constructor --------------------------------------

#[test]
fn the_constructor_reports_each_fault() {
    assert!(BSpline::new(&[0.0, 0.0, 1.0, 2.0, 2.0], &[0.0, 1.0, 0.0], 1).is_ok());

    let cases: [(&[f64], &[f64], usize, E); 8] = [
        (
            &[0.0; 34],
            &[0.0; 17],
            16,
            E::InvalidParameter {
                name: "degree",
                reason: ParamReason::TooLarge { max: 15 },
            },
        ),
        (
            &[0.0, 0.0, 1.0, 1.0, 1.0],
            &[0.0, 1.0],
            1,
            E::LengthMismatch {
                expected: 4,
                actual: 5,
            },
        ),
        (
            &[0.0, 1.0, 2.0],
            &[1.0],
            1,
            E::InsufficientPoints {
                required: 2,
                actual: 1,
            },
        ),
        (
            &[0.0, f64::NAN, 1.0, 1.0],
            &[0.0, 1.0],
            1,
            E::NonFinite { index: 1 },
        ),
        (
            &[0.0, 0.0, 1.0, 1.0],
            &[0.0, f64::INFINITY],
            1,
            E::NonFinite { index: 1 },
        ),
        (
            &[0.0, 0.0, 1.0, 0.5],
            &[0.0, 1.0],
            1,
            E::NotSorted { index: 3 },
        ),
        // The domain is from knot 1 to knot 2, and they are equal.
        (
            &[0.0, 1.0, 1.0, 2.0],
            &[0.0, 1.0],
            1,
            E::DuplicatePoints {
                first: 1,
                second: 2,
            },
        ),
        // Every knot equal: in order, but the domain is empty.
        (
            &[3.0; 6],
            &[0.0, 1.0, 2.0],
            2,
            E::DuplicatePoints {
                first: 2,
                second: 3,
            },
        ),
    ];
    for (knots, control_points, degree, want) in cases {
        assert_eq!(
            BSpline::new(knots, control_points, degree).err(),
            Some(want)
        );
        assert_eq!(
            BSpline::from_owned(knots.to_vec(), control_points.to_vec(), degree).err(),
            Some(want)
        );
    }
}

#[test]
fn the_constructor_keeps_the_order_of_the_checks() {
    let too_large = E::InvalidParameter {
        name: "degree",
        reason: ParamReason::TooLarge { max: 15 },
    };
    // The degree is over the limit, and the knots have the wrong length.
    assert_eq!(BSpline::new(&[0.0, 1.0], &[1.0], 16).err(), Some(too_large));
    // The wrong number of knots, and too few control points.
    assert_eq!(
        BSpline::new(&[0.0, 1.0], &[1.0], 1).err(),
        Some(E::LengthMismatch {
            expected: 3,
            actual: 2
        })
    );
    // Too few control points, and a NaN knot.
    assert_eq!(
        BSpline::new(&[f64::NAN, 1.0, 2.0, 3.0], &[1.0], 2).err(),
        Some(E::InsufficientPoints {
            required: 3,
            actual: 1
        })
    );
    // A NaN control point, and decreasing knots: knots are checked first for
    // NaN, then control points, then the order.
    assert_eq!(
        BSpline::new(&[1.0, 0.0, 2.0, 3.0], &[f64::NAN, 1.0], 1).err(),
        Some(E::NonFinite { index: 0 })
    );
    // Decreasing knots, and an empty domain.
    assert_eq!(
        BSpline::new(&[0.0, 1.0, 1.0, 0.5], &[0.0, 1.0], 1).err(),
        Some(E::NotSorted { index: 3 })
    );
    // The largest degree is accepted.
    let k = BSPLINE_MAX_DEGREE;
    assert_eq!(k, 15);
    let knots: Vec<f64> = [vec![0.0; k + 1], vec![1.0; k + 1]].concat();
    assert!(BSpline::new(&knots, &vec![1.0; k + 1], k).is_ok());
}

#[test]
fn from_ppoly_rejects_a_degree_a_b_spline_cannot_have() {
    let poly = PPoly::new(&[0.0, 1.0], &[1.0; 17], 16).expect("valid");
    assert_eq!(
        BSpline::from_ppoly(&poly).err(),
        Some(E::InvalidParameter {
            name: "degree",
            reason: ParamReason::TooLarge { max: 15 },
        })
    );
}

// --- Knots, intervals, and extrapolation -----------------------------------------

/// The polygon through (0, 0), (1, 1), (2, 0).
fn hat() -> Result<BSpline<f64>, E> {
    BSpline::new(&[0.0, 0.0, 1.0, 2.0, 2.0], &[0.0, 1.0, 0.0], 1)
}

#[test]
fn each_mode_treats_the_outside_as_documented() {
    let spline = hat().expect("valid");
    assert_eq!(spline.extrapolate(), Extrapolate::Extend, "SciPy's default");
    assert_eq!(spline.domain().0.to_bits(), 0.0_f64.to_bits());
    assert_eq!(spline.domain().1.to_bits(), 2.0_f64.to_bits());
    let at = |mode, x: f64| bits(spline.clone().with_extrapolate(mode).eval(&x));

    assert_eq!(at(Extrapolate::Forbid, -0.5), None);
    assert_eq!(at(Extrapolate::Forbid, 2.5), None);
    assert_eq!(at(Extrapolate::Clamp, -0.5), bits(Some(0.0)));
    assert_eq!(at(Extrapolate::Clamp, 2.5), bits(Some(0.0)));
    // The end pieces continue: x on the left, 2 - x on the right.
    assert_eq!(at(Extrapolate::Extend, -0.5), bits(Some(-0.5)));
    assert_eq!(at(Extrapolate::Extend, 2.5), bits(Some(-0.5)));
    // The period is 2.
    assert_eq!(at(Extrapolate::Periodic, 2.5), bits(Some(0.5)));
    assert_eq!(at(Extrapolate::Periodic, -0.5), bits(Some(0.5)));

    for mode in MODES {
        assert_eq!(at(mode, 0.0), bits(Some(0.0)), "{mode:?}");
        assert_eq!(at(mode, 1.0), bits(Some(1.0)), "{mode:?}");
        assert_eq!(at(mode, 2.0), bits(Some(0.0)), "{mode:?}");
        assert!(spline
            .clone()
            .with_extrapolate(mode)
            .eval(&f64::NAN)
            .is_none());
    }
}

#[test]
fn a_repeated_inner_knot_belongs_to_the_interval_on_its_right() {
    // Degree 1 with the knot 1 twice: the spline jumps there, from the line
    // 0 -> 1 on [0, 1) to the line 5 -> 6 on [1, 2].
    let spline =
        BSpline::new(&[0.0, 0.0, 1.0, 1.0, 2.0, 2.0], &[0.0, 1.0, 5.0, 6.0], 1).expect("valid");
    assert_eq!(bits(spline.eval(&0.5)), bits(Some(0.5)));
    assert_eq!(bits(spline.eval(&1.0)), bits(Some(5.0)));
    assert_eq!(bits(spline.eval(&1.5)), bits(Some(5.5)));
    assert_eq!(bits(spline.eval(&2.0)), bits(Some(6.0)));
    let mut hint = Hint::new();
    for x in [0.5, 1.0, 1.5, 2.0, 0.0, 1.0] {
        assert_eq!(
            bits(spline.eval_hinted(x, &mut hint)),
            bits(spline.eval(&x))
        );
    }

    // The same function as a piecewise polynomial, and back.
    let poly = PPoly::from_bspline(&spline).expect("finite");
    assert_eq!(poly.breakpoints().len(), 3);
    assert_eq!(bits(poly.eval(&1.0)), bits(Some(5.0)));
    let back = BSpline::from_ppoly(&poly).expect("degree 1");
    assert_eq!(back.knots().len(), 6);
    for x in [0.0, 0.5, 1.0, 1.5, 2.0] {
        assert_eq!(bits(back.eval(&x)), bits(spline.eval(&x)), "x = {x}");
    }
}

#[test]
fn knots_outside_the_domain_do_not_move_it() {
    // Degree 2 without repeated end knots: the domain is [t[2], t[4]] = [2, 4].
    let spline: BSpline<f64> =
        BSpline::new(&[0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0], &[1.0; 4], 2).expect("valid");
    assert_eq!(spline.domain().0.to_bits(), 2.0_f64.to_bits());
    assert_eq!(spline.domain().1.to_bits(), 4.0_f64.to_bits());
    // The basis functions sum to 1 on the domain, so equal control points give
    // that value.
    for x in [2.0, 2.5, 3.0, 3.75, 4.0] {
        assert_close(
            spline.eval(&x).expect("inside"),
            1.0,
            DIRECT_F64,
            "sum of the basis",
        );
    }
    let forbid = spline.with_extrapolate(Extrapolate::Forbid);
    assert!(
        forbid.eval(&1.5).is_none(),
        "between knots, but outside the domain"
    );
    assert!(forbid.eval(&4.5).is_none());
}

#[test]
fn end_knots_repeated_more_than_needed_leave_an_empty_interval_at_the_end() {
    // Degree 1 needs the end knots twice. Three times leaves a knot interval of
    // zero width at that end of the domain, which no point belongs to: the end
    // pieces, and so the values outside, come from the interval next to it.

    // Left end: the domain [0, 1] is the line from 2 to 5. The first control
    // point has no basis function on the domain.
    let left: BSpline<f64> =
        BSpline::new(&[0.0, 0.0, 0.0, 1.0, 1.0], &[100.0, 2.0, 5.0], 1).expect("valid");
    assert_eq!(left.domain().0.to_bits(), 0.0_f64.to_bits());
    assert_eq!(bits(left.eval(&0.0)), bits(Some(2.0)));
    assert_eq!(bits(left.eval(&0.5)), bits(Some(3.5)));
    assert_eq!(bits(left.eval(&1.0)), bits(Some(5.0)));
    assert_eq!(
        bits(left.eval(&-1.0)),
        bits(Some(-1.0)),
        "the line continued"
    );
    let clamp = left.clone().with_extrapolate(Extrapolate::Clamp);
    assert_eq!(bits(clamp.eval(&-1.0)), bits(Some(2.0)));

    // Right end: the line from 2 to 5 again, and the last control point idle.
    let right: BSpline<f64> =
        BSpline::new(&[0.0, 0.0, 1.0, 1.0, 1.0], &[2.0, 5.0, 100.0], 1).expect("valid");
    assert_eq!(right.domain().1.to_bits(), 1.0_f64.to_bits());
    assert_eq!(bits(right.eval(&0.0)), bits(Some(2.0)));
    assert_eq!(
        bits(right.eval(&1.0)),
        bits(Some(5.0)),
        "the right end of the domain"
    );
    assert_eq!(
        bits(right.eval(&2.0)),
        bits(Some(8.0)),
        "the line continued"
    );
    let clamp = right.clone().with_extrapolate(Extrapolate::Clamp);
    assert_eq!(bits(clamp.eval(&2.0)), bits(Some(5.0)));

    // Both as piecewise polynomials: one piece each, no interval of zero width.
    for spline in [&left, &right] {
        let poly = PPoly::from_bspline(spline).expect("finite");
        assert_eq!(poly.breakpoints().len(), 2);
        let mut hint = Hint::new();
        for x in [-1.0, 0.0, 0.25, 1.0, 2.0] {
            assert_eq!(bits(poly.eval(&x)), bits(spline.eval(&x)), "x = {x}");
            assert_eq!(
                bits(spline.eval_hinted(x, &mut hint)),
                bits(spline.eval(&x))
            );
        }
    }
}

#[test]
fn accessors_and_thread_safety() {
    let spline: BSpline<f32, [f32; 2]> =
        BSpline::new(&[0.0, 0.0, 1.0, 1.0], &[[0.0, 1.0], [2.0, 5.0]], 1).expect("valid");
    assert_eq!(spline.degree(), 1);
    assert_eq!(spline.knots().len(), 4);
    assert_eq!(spline.control_points().len(), 2);
    assert_eq!(
        spline.eval(&0.5).map(|v| v.map(f32::to_bits)),
        Some([1.0_f32, 3.0].map(f32::to_bits))
    );

    fn assert_send_sync<T: Send + Sync>() {}
    assert_send_sync::<BSpline<f64>>();
    assert_send_sync::<BSpline<f32, [f32; 3]>>();
    assert_send_sync::<ParamReason>();
}

#[test]
fn the_batch_rejects_a_length_mismatch_and_leaves_the_buffer() {
    let spline = hat().expect("valid");
    for (n_points, n_out) in [(2, 3), (3, 2), (0, 1), (1, 0)] {
        let points = vec![0.5; n_points];
        let mut out = vec![7.0; n_out];
        assert_eq!(
            spline.eval_batch_into(&points, &mut out, 0.0),
            Err(E::LengthMismatch {
                expected: n_points,
                actual: n_out,
            })
        );
        assert!(out.iter().all(|v| v.to_bits() == 7.0_f64.to_bits()));
    }
}

// --- Properties ----------------------------------------------------------------

/// Knots, control points, and a degree of a random B-spline. The end knots are
/// repeated `degree + 1` times; inner knots come from a coarse grid, so some of
/// them repeat.
fn splines() -> impl Strategy<Value = (Vec<f64>, Vec<f64>, usize)> {
    (0_usize..5, 0_usize..7)
        .prop_flat_map(|(degree, inner)| {
            (
                prop::collection::vec(1_u8..12, inner),
                prop::collection::vec(-10.0_f64..10.0, degree + 1 + inner),
                Just(degree),
            )
        })
        .prop_map(|(mut inner, control_points, degree)| {
            inner.sort_unstable();
            // No inner knot more than `degree + 1` times: more would leave a
            // control point with no basis function.
            let mut kept: Vec<u8> = Vec::new();
            for knot in inner {
                if kept.iter().filter(|k| **k == knot).count() <= degree {
                    kept.push(knot);
                }
            }
            let n = degree + 1 + kept.len();
            let mut knots = vec![0.0; degree + 1];
            knots.extend(kept.iter().map(|k| f64::from(*k) * 0.25));
            knots.extend(std::iter::repeat_n(3.0, degree + 1));
            (knots, control_points[..n].to_vec(), degree)
        })
}

fn queries() -> impl Strategy<Value = Vec<f64>> {
    prop::collection::vec(-4.0_f64..8.0, 1..150)
}

proptest! {
    #[test]
    fn a_hint_never_changes_the_result(
        (knots, control_points, degree) in splines(),
        xs in queries(),
    ) {
        for mode in MODES {
            let spline = BSpline::new(&knots, &control_points, degree)
                .expect("valid")
                .with_extrapolate(mode);
            let mut hint = Hint::new();
            for x in &xs {
                prop_assert_eq!(bits(spline.eval_hinted(*x, &mut hint)), bits(spline.eval(x)), "{:?} x = {}", mode, x);
            }
        }
    }

    #[test]
    fn the_batch_equals_the_loop_in_every_order(
        (knots, control_points, degree) in splines(),
        xs in queries(),
    ) {
        let mut ascending = xs.clone();
        ascending.sort_by(f64::total_cmp);
        let descending: Vec<f64> = ascending.iter().rev().copied().collect();
        let dense: Vec<f64> = (0..xs.len()).map(|i| 3.0 * i as f64 / xs.len() as f64).collect();
        for mode in MODES {
            let spline = BSpline::new(&knots, &control_points, degree)
                .expect("valid")
                .with_extrapolate(mode);
            for batch in [&xs, &ascending, &descending, &dense] {
                let mut out = vec![0.0; batch.len()];
                prop_assert_eq!(spline.eval_batch_into(batch, &mut out, -7.0), Ok(()));
                for (x, got) in batch.iter().zip(&out) {
                    prop_assert_eq!(got.to_bits(), spline.eval_or(x, -7.0).to_bits(), "{:?} x = {}", mode, x);
                }
            }
        }
    }

    #[test]
    fn equal_control_points_give_that_value_on_the_domain(
        (knots, control_points, degree) in splines(),
        value in -10.0_f64..10.0,
        t in 0.0_f64..=1.0,
    ) {
        // The basis functions sum to 1 wherever at least one is defined.
        let constant = vec![value; control_points.len()];
        let spline = BSpline::new(&knots, &constant, degree).expect("valid");
        let got = spline.eval(&(3.0 * t)).expect("inside the domain");
        prop_assert!((got - value).abs() <= 1e-12 * value.abs().max(1.0), "{} for {}", got, value);
    }

    #[test]
    fn both_conversions_keep_the_values(
        (knots, control_points, degree) in splines(),
        xs in prop::collection::vec(0.0_f64..=3.0, 1..40),
    ) {
        let spline = BSpline::new(&knots, &control_points, degree).expect("valid");
        let poly = PPoly::from_bspline(&spline).expect("finite coefficients");
        let back = BSpline::from_ppoly(&poly).expect("a degree a B-spline can have");
        prop_assert_eq!(poly.degree(), degree);
        prop_assert_eq!(back.degree(), degree);
        // Inner knots are 0.25 apart at least, and the degree is 4 at most: a
        // derivative of order 4 is at most 10 * (4 / 0.25)^4 times a small factor.
        let scale = control_points.iter().fold(1.0_f64, |m, c| m.max(c.abs()));
        let tolerance = 1e-9 * scale;
        for x in &xs {
            let want = spline.eval(x).expect("inside the domain");
            let as_poly = poly.eval(x).expect("inside the domain");
            let as_spline = back.eval(x).expect("inside the domain");
            prop_assert!((as_poly - want).abs() <= tolerance, "PPoly at {}: {} for {}", x, as_poly, want);
            prop_assert!((as_spline - want).abs() <= tolerance, "BSpline at {}: {} for {}", x, as_spline, want);
        }
    }
}
