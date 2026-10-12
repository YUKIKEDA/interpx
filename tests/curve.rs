//! `BSplineCurve`: the golden data, the parameters each way chooses, the input
//! checks through the two constructors, and the properties of a curve.

mod common;

use common::{assert_close, load, reals, Real, DIRECT_F64};
use interpx::{
    BSpline1D, BSplineCurve, BoundaryCondition as Bc, EndCondition as End, Extrapolate, Hint,
    InterpolateError as E, Interpolator, ParamReason, Parametrization as Way,
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

// Marked as test code so that `expect` is allowed here as in a `#[test]` function
// (clippy.toml allows it in tests only). The same for the helpers below.
#[cfg(test)]
fn points_of<const D: usize>(v: &Json) -> Vec<[f64; D]> {
    rows(v)
        .iter()
        .map(|r| r.as_slice().try_into().expect("D columns"))
        .collect()
}

fn bits<const D: usize>(v: Option<[f64; D]>) -> Option<[u64; D]> {
    v.map(|point| point.map(f64::to_bits))
}

/// SciPy's `bc_type` for points of `D` coordinates: no condition is not-a-knot.
#[cfg(test)]
fn bc_of<const D: usize>(v: &Json) -> Bc<[f64; D]> {
    match v.as_str() {
        None => Bc::NotAKnot,
        Some("periodic") => Bc::Periodic,
        Some("natural") => Bc::natural(),
        Some("clamped") => Bc::clamped(),
        Some(other) => panic!("a boundary condition {other}"),
    }
}

#[cfg(test)]
fn way_of(v: &Json) -> Option<Way> {
    match v.as_str()? {
        "uniform" => Some(Way::Uniform),
        "chord_length" => Some(Way::ChordLength),
        "centripetal" => Some(Way::Centripetal),
        other => panic!("a parametrization {other}"),
    }
}

/// Rebuilds one golden case of `D` coordinates and compares it with SciPy.
/// Returns whether its parameters were chosen from the points.
#[cfg(test)]
fn check<const D: usize>(case: &Json) -> bool {
    let name = case["name"].as_str().expect("name");
    let points = points_of::<D>(&case["points"]);
    let degree = usize::try_from(case["options"]["k"].as_u64().expect("k")).expect("small");
    let bc = bc_of::<D>(&case["options"]["bc_type"]);
    let way = way_of(&case["options"]["parametrization"]);
    let curve = match way {
        Some(way) => BSplineCurve::new(&points, degree, bc, way).expect("valid points"),
        None => {
            let t = nums(&case["t"]);
            let curve = BSplineCurve::with_parameters(&t, &points, degree, bc).expect("valid");
            // The parameters given, bit for bit.
            assert_eq!(curve.parameters().len(), t.len(), "{name}");
            for (got, want) in curve.parameters().iter().zip(&t) {
                assert_eq!(got.to_bits(), want.to_bits(), "{name}");
            }
            curve
        }
    };
    let default = if bc == Bc::Periodic {
        Extrapolate::Periodic
    } else {
        Extrapolate::Extend
    };
    assert_eq!(curve.extrapolate(), default, "{name}: default mode");
    assert_eq!(curve.degree(), degree, "{name}");

    if way.is_some() {
        // The parameters the generator computed by the same formulas.
        let parameters = nums(&case["parameters"]);
        assert_eq!(curve.parameters().len(), parameters.len(), "{name}");
        for (i, (got, want)) in curve.parameters().iter().zip(&parameters).enumerate() {
            assert_close(*got, *want, DIRECT_F64, &format!("{name} parameter {i}"));
        }
        assert_eq!(curve.parameters()[0].to_bits(), 0.0_f64.to_bits(), "{name}");
        let last = curve.parameters()[parameters.len() - 1];
        assert_eq!(last.to_bits(), 1.0_f64.to_bits(), "{name}");
    }

    let b = curve.as_bspline();
    let knots = nums(&case["knots"]);
    assert_eq!(b.knots().len(), knots.len(), "{name}: knots");
    for (i, (got, want)) in b.knots().iter().zip(&knots).enumerate() {
        assert_close(*got, *want, DIRECT_F64, &format!("{name} knot {i}"));
    }
    let control_points = points_of::<D>(&case["control_points"]);
    assert_eq!(b.control_points().len(), control_points.len(), "{name}");
    for (i, (got, want)) in b.control_points().iter().zip(&control_points).enumerate() {
        for column in 0..D {
            let what = format!("{name} control point {i} column {column}");
            assert_close(got[column], want[column], DIRECT_F64, &what);
        }
    }

    let queries = nums(&case["queries"]);
    for (field, nu) in [("expected", 0), ("derivative_1", 1), ("derivative_2", 2)] {
        let Some(want) = case.get(field) else {
            assert_eq!(degree, 1, "{name}: only degree 1 has no second derivative");
            continue;
        };
        for (k, (q, want)) in queries.iter().zip(points_of::<D>(want)).enumerate() {
            // The derivative of the order of the degree is constant between
            // the knots and jumps at them. With chosen parameters a knot is a
            // rounding away from SciPy's (the distances are summed
            // differently), so a query that SciPy has on a knot may lie on
            // the other side of it here: nothing to compare.
            let on_a_knot = knots.iter().any(|knot| (knot - q).abs() <= 1e-12);
            if way.is_some() && nu == degree && on_a_knot {
                continue;
            }
            let got = b.eval_derivative(*q, nu).expect("the cases extrapolate");
            for column in 0..D {
                let what = format!("{name} {field} query {k} column {column}");
                assert_close(got[column], want[column], DIRECT_F64, &what);
            }
            if nu == 0 {
                assert_eq!(bits(curve.eval(q)), bits(Some(got)), "{name} query {k}");
            }
        }
    }
    way.is_some()
}

// --- Golden data -------------------------------------------------------------

#[test]
fn scipy_curves_are_rebuilt_from_their_points() {
    let file = load::<Json>("curve").expect("golden file curve");
    let mut given = [0; 5];
    let mut chosen = [0; 5];
    for case in &file.cases {
        let dimension = rows(&case["points"]).first().map_or(0, Vec::len);
        let from_points = match dimension {
            2 => check::<2>(case),
            3 => check::<3>(case),
            4 => check::<4>(case),
            other => panic!("a curve in {other} dimensions"),
        };
        if from_points {
            chosen[dimension] += 1;
        } else {
            given[dimension] += 1;
        }
    }
    // With parameters given: 5 in 2D, 6 in 3D. With parameters chosen from
    // the points: 10 in 2D, 4 in 3D, 1 in 4D.
    assert_eq!(given, [0, 0, 5, 6, 0]);
    assert_eq!(chosen, [0, 0, 10, 4, 1]);
}

// --- The parameters ---------------------------------------------------------------

#[test]
fn each_way_chooses_its_parameters() {
    // Distances 5 (a 3-4-5 triangle), 1, and 4.
    let points = [[0.0, 0.0], [3.0, 4.0], [3.0, 5.0], [7.0, 5.0]];
    let parameters = |way: Way| {
        BSplineCurve::new(&points, 1, Bc::NotAKnot, way)
            .expect("valid")
            .parameters()
            .to_vec()
    };
    let exact = |got: Vec<f64>, want: [f64; 4]| {
        assert_eq!(
            got.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
            want.map(f64::to_bits)
        );
    };
    exact(parameters(Way::Uniform), [0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0]);
    exact(parameters(Way::ChordLength), [0.0, 0.5, 0.6, 1.0]);
    // Square roots: sqrt(5), 1, 2.
    let root = 5.0_f64.sqrt();
    let total = root + 1.0 + 2.0;
    let centripetal = parameters(Way::Centripetal);
    for (got, want) in centripetal
        .iter()
        .zip([0.0, root / total, (root + 1.0) / total, 1.0])
    {
        assert_close(*got, want, DIRECT_F64, "centripetal");
    }
    // The domain of a curve with chosen parameters is [0, 1].
    let curve = BSplineCurve::new(&points, 3, Bc::NotAKnot, Way::ChordLength).expect("valid");
    assert_eq!(curve.domain().0.to_bits(), 0.0_f64.to_bits());
    assert_eq!(curve.domain().1.to_bits(), 1.0_f64.to_bits());
}

#[test]
fn the_distance_neither_overflows_nor_vanishes() {
    // Coordinates whose squares overflow: distances 3e200, 4e200 (the legs) and
    // then 5e200 (the hypotenuse back).
    let large = [[0.0, 0.0], [3e200, 0.0], [3e200, 4e200], [0.0, 0.0]];
    let curve = BSplineCurve::new(&large, 1, Bc::NotAKnot, Way::ChordLength).expect("finite");
    for (got, want) in curve.parameters().iter().zip([0.0, 0.25, 7.0 / 12.0, 1.0]) {
        assert_close(*got, want, DIRECT_F64, "large");
    }
    // Differences whose squares are zero: the same triangle at 1e-200.
    let small = large.map(|p: [f64; 2]| p.map(|v| v / 1e200 / 1e200));
    assert!(small[1][0] > 0.0 && small[1][0] * small[1][0] == 0.0);
    let curve = BSplineCurve::new(&small, 1, Bc::NotAKnot, Way::ChordLength).expect("distinct");
    for (got, want) in curve.parameters().iter().zip([0.0, 0.25, 7.0 / 12.0, 1.0]) {
        assert_close(*got, want, DIRECT_F64, "small");
    }
    let curve = BSplineCurve::new(&small, 1, Bc::NotAKnot, Way::Centripetal).expect("distinct");
    assert!(curve.parameters().windows(2).all(|w| w[0] < w[1]));
    // Points so far apart that their difference overflows: no finite distance.
    let apart = [[-1.5e308, 0.0], [1.5e308, 0.0], [1.5e308, 1.0]];
    // The later of the two points is named.
    for way in [Way::ChordLength, Way::Centripetal] {
        let err = BSplineCurve::new(&apart, 1, Bc::NotAKnot, way).err();
        assert_eq!(err, Some(E::NonFinite { index: 1 }), "{way:?}");
    }
    // The sum of finite distances overflows.
    let long = [[-1.0e308, 0.0], [0.0, 0.0], [1.0e308, 0.0], [0.0, 1.0e308]];
    // The sum is finite up to the second point, and not from the third on.
    let err = BSplineCurve::new(&long, 1, Bc::NotAKnot, Way::ChordLength).err();
    assert_eq!(err, Some(E::NonFinite { index: 2 }));
    // Uniform parameters do not look at the distances.
    assert!(BSplineCurve::new(&apart, 1, Bc::NotAKnot, Way::Uniform).is_ok());
}

// --- Input checks -------------------------------------------------------------------

fn invalid(name: &'static str, reason: ParamReason) -> E {
    E::InvalidParameter { name, reason }
}

#[test]
fn new_reports_each_fault_in_the_order_of_the_design() {
    let square = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.0, 0.0]];
    let new = |points: &[[f64; 2]], degree: usize, bc: Bc<[f64; 2]>, way: Way| {
        BSplineCurve::new(points, degree, bc, way).err()
    };
    let chord = Way::ChordLength;
    assert_eq!(new(&square, 3, Bc::Periodic, chord), None);
    assert_eq!(
        new(&square, 0, Bc::NotAKnot, chord),
        Some(invalid("degree", ParamReason::TooSmall { min: 1 }))
    );
    assert_eq!(
        new(&square, 16, Bc::NotAKnot, chord),
        Some(invalid("degree", ParamReason::TooLarge { max: 15 }))
    );
    assert_eq!(
        new(&square, 2, Bc::natural(), chord),
        Some(invalid("bc", ParamReason::Unsupported))
    );
    assert_eq!(
        new(&square, 5, Bc::NotAKnot, chord),
        Some(E::InsufficientPoints {
            required: 6,
            actual: 5
        })
    );
    assert_eq!(
        new(&[], 1, Bc::NotAKnot, chord),
        Some(E::InsufficientPoints {
            required: 2,
            actual: 0
        })
    );
    assert_eq!(
        new(&square[..1], 3, Bc::Periodic, Way::Uniform),
        Some(E::InsufficientPoints {
            required: 2,
            actual: 1
        })
    );
    let mut nan = square;
    nan[3][1] = f64::NAN;
    for way in [Way::Uniform, chord, Way::Centripetal] {
        assert_eq!(
            new(&nan, 3, Bc::NotAKnot, way),
            Some(E::NonFinite { index: 3 })
        );
    }
    let mut infinite = square;
    infinite[2][0] = f64::NEG_INFINITY;
    assert_eq!(
        new(&infinite, 1, Bc::NotAKnot, chord),
        Some(E::NonFinite { index: 2 })
    );
    // The same point twice in a row.
    let stalled = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [1.0, 1.0], [0.0, 0.0]];
    let duplicate = E::DuplicatePoints {
        first: 2,
        second: 3,
    };
    assert_eq!(new(&stalled, 1, Bc::NotAKnot, chord), Some(duplicate));
    assert_eq!(
        new(&stalled, 1, Bc::NotAKnot, Way::Centripetal),
        Some(duplicate)
    );
    assert_eq!(new(&stalled, 1, Bc::NotAKnot, Way::Uniform), None);
    // Every point the same: there is no length at all to divide by, and the
    // first pair is named.
    let still = [[2.0, -1.0]; 4];
    let first_pair = E::DuplicatePoints {
        first: 0,
        second: 1,
    };
    assert_eq!(new(&still, 1, Bc::NotAKnot, chord), Some(first_pair));
    assert_eq!(
        new(&still, 1, Bc::NotAKnot, Way::Centripetal),
        Some(first_pair)
    );
    assert_eq!(new(&still, 1, Bc::NotAKnot, Way::Uniform), None);
    // An open curve under the periodic condition.
    let open = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    assert_eq!(new(&open, 3, Bc::Periodic, chord), Some(E::NotPeriodic));
    // A derivative that is NaN, degree 3.
    let bad = Bc::Ends(End::FirstDerivative([f64::NAN, 0.0]), End::NotAKnot);
    assert_eq!(
        new(&square, 3, bad, chord),
        Some(invalid("bc", ParamReason::NonFinite))
    );

    // The order: the degree before the boundary condition, both before the
    // count, the count before a NaN, a NaN before equal points, equal points
    // before the periodic check.
    assert_eq!(
        new(&[], 0, Bc::natural(), chord),
        Some(invalid("degree", ParamReason::TooSmall { min: 1 }))
    );
    assert_eq!(
        new(&[], 2, Bc::natural(), chord),
        Some(invalid("bc", ParamReason::Unsupported))
    );
    assert_eq!(
        new(&[[f64::NAN, 0.0]], 1, Bc::NotAKnot, chord),
        Some(E::InsufficientPoints {
            required: 2,
            actual: 1
        })
    );
    let mut both = stalled;
    both[4][0] = f64::NAN;
    assert_eq!(
        new(&both, 1, Bc::NotAKnot, chord),
        Some(E::NonFinite { index: 4 })
    );
    let mut open_and_stalled = stalled;
    open_and_stalled[4] = [5.0, 5.0];
    assert_eq!(
        new(&open_and_stalled, 3, Bc::Periodic, chord),
        Some(duplicate)
    );
}

#[test]
fn points_too_close_for_the_length_of_the_whole_are_duplicates() {
    // Two points 1e-20 apart on a curve of length 2: their parameters round to
    // the same number, though the points differ.
    let points = [[0.0, 0.0], [1.0, 0.0], [1.0, 1e-20], [2.0, 0.0]];
    let err = BSplineCurve::new(&points, 1, Bc::NotAKnot, Way::ChordLength).err();
    assert_eq!(
        err,
        Some(E::DuplicatePoints {
            first: 1,
            second: 2
        })
    );
    // The square root of the distance is 1e-10: not too close.
    assert!(BSplineCurve::new(&points, 1, Bc::NotAKnot, Way::Centripetal).is_ok());
}

#[test]
fn with_parameters_reports_the_faults_of_the_spline() {
    let points = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let with = |t: &[f64], degree: usize| {
        BSplineCurve::with_parameters(t, &points, degree, Bc::NotAKnot).err()
    };
    assert_eq!(with(&[0.0, 1.0, 2.0, 3.0], 3), None);
    assert_eq!(
        with(&[0.0, 1.0, 2.0], 3),
        Some(E::LengthMismatch {
            expected: 3,
            actual: 4
        })
    );
    assert_eq!(
        with(&[0.0, 2.0, 1.0, 3.0], 3),
        Some(E::NotSorted { index: 2 })
    );
    assert_eq!(
        with(&[0.0, 1.0, 1.0, 3.0], 3),
        Some(E::DuplicatePoints {
            first: 1,
            second: 2
        })
    );
    assert_eq!(
        with(&[0.0, f64::NAN, 2.0, 3.0], 3),
        Some(E::NonFinite { index: 1 })
    );
    assert_eq!(
        with(&[0.0, 1.0, 2.0, 3.0], 4),
        Some(E::InsufficientPoints {
            required: 5,
            actual: 4
        })
    );
}

#[test]
fn modes_f32_dimensions_and_thread_safety() {
    let points = [[0.0, 0.0], [1.0, 2.0], [3.0, 2.0], [4.0, 0.0]];
    let curve = BSplineCurve::new(&points, 3, Bc::NotAKnot, Way::ChordLength).expect("valid");
    assert_eq!(curve.extrapolate(), Extrapolate::Extend);
    assert!(curve.eval(&1.5).is_some());
    let forbid = curve.clone().with_extrapolate(Extrapolate::Forbid);
    assert_eq!(forbid.extrapolate(), Extrapolate::Forbid);
    assert_eq!(forbid.as_bspline().extrapolate(), Extrapolate::Forbid);
    assert_eq!(forbid.eval(&1.5), None);
    assert_eq!(forbid.eval(&f64::NAN), None);
    assert_eq!(forbid.parameters(), curve.parameters());
    let clamp = curve.clone().with_extrapolate(Extrapolate::Clamp);
    assert_eq!(bits(clamp.eval(&1.5)), bits(curve.eval(&1.0)));
    let periodic = curve.with_extrapolate(Extrapolate::Periodic);
    assert_eq!(bits(periodic.eval(&1.25)), bits(periodic.eval(&0.25)));

    let single: BSplineCurve<f32, 3> = BSplineCurve::new(
        &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 2.0, 0.0]],
        1,
        Bc::NotAKnot,
        Way::ChordLength,
    )
    .expect("valid");
    assert_eq!(single.parameters()[1].to_bits(), (1.0_f32 / 3.0).to_bits());
    assert_eq!(
        single.eval(&0.5).map(|p| p.map(f32::to_bits)),
        Some([1.0_f32, 0.5, 0.0].map(f32::to_bits))
    );
    // One dimension: a curve on a line.
    let line: BSplineCurve<f64, 1> =
        BSplineCurve::new(&[[0.0], [1.0], [4.0]], 2, Bc::NotAKnot, Way::ChordLength)
            .expect("valid");
    let at = line.eval(&0.5).expect("inside");
    assert_close(at[0], 2.0, DIRECT_F64, "the parameter is the length");

    fn assert_send_sync<T: Send + Sync>() {}
    assert_send_sync::<BSplineCurve<f64, 2>>();
    assert_send_sync::<BSplineCurve<f32, 7>>();
}

// --- Properties ----------------------------------------------------------------

/// Points in the plane with neighbors at least 0.1 apart in x: no two equal.
fn polyline() -> impl Strategy<Value = Vec<[f64; 2]>> {
    prop::collection::vec((0.1_f64..2.0, -3.0_f64..3.0), 7..14).prop_map(|steps| {
        let mut x = 0.0;
        steps
            .into_iter()
            .map(|(step, y)| {
                x += step;
                [x * if y > 0.0 { 1.0 } else { 0.7 }, y]
            })
            .collect()
    })
}

fn ways() -> impl Strategy<Value = Way> {
    prop::sample::select(vec![Way::Uniform, Way::ChordLength, Way::Centripetal])
}

/// The largest coordinate of the control points, and 1 at least.
fn size(curve: &BSplineCurve<f64, 2>) -> f64 {
    curve
        .as_bspline()
        .control_points()
        .iter()
        .flatten()
        .fold(1.0_f64, |m, c| m.max(c.abs()))
}

proptest! {
    #[test]
    fn the_curve_is_at_each_point_at_its_parameter(
        points in polyline(),
        way in ways(),
        degree in 1..=6_usize,
        closed in any::<bool>(),
    ) {
        let mut points = points;
        if closed {
            points.push(points[0]);
        }
        let bc = if closed { Bc::Periodic } else { Bc::NotAKnot };
        let curve = BSplineCurve::new(&points, degree, bc, way).expect("valid points");
        let t = curve.parameters();
        prop_assert_eq!(t.len(), points.len());
        prop_assert_eq!(t[0].to_bits(), 0.0_f64.to_bits());
        prop_assert_eq!(t[t.len() - 1].to_bits(), 1.0_f64.to_bits());
        prop_assert!(t.windows(2).all(|w| w[0] < w[1]));
        let tol = 1e-11 * size(&curve);
        for (u, point) in t.iter().zip(&points) {
            let at = curve.eval(u).expect("inside");
            prop_assert!((at[0] - point[0]).abs() <= tol && (at[1] - point[1]).abs() <= tol, "at {}: {:?} for {:?}", u, at, point);
        }
    }

    #[test]
    fn given_parameters_make_the_spline_of_the_same_arguments(
        points in polyline(),
        degree in 1..=6_usize,
        qs in prop::collection::vec(-0.5_f64..12.0, 1..40),
    ) {
        // Parameters of the caller's own: about 3 apart, moved by the points.
        let t: Vec<f64> = points
            .iter()
            .zip(0_u8..)
            .map(|(p, i)| 3.0 * f64::from(i) + 0.3 * p[1].abs())
            .collect();
        let curve = BSplineCurve::with_parameters(&t, &points, degree, Bc::NotAKnot).expect("valid");
        let spline = BSpline1D::new(&t, &points, degree, Bc::NotAKnot).expect("valid");
        for q in &qs {
            prop_assert_eq!(bits(curve.eval(q)), bits(spline.eval(q)));
        }
        // And the curve with chosen parameters is the one with those
        // parameters given.
        let chosen = BSplineCurve::new(&points, degree, Bc::NotAKnot, Way::Centripetal).expect("valid");
        let again = BSplineCurve::with_parameters(chosen.parameters(), &points, degree, Bc::NotAKnot).expect("valid");
        for q in &qs {
            prop_assert_eq!(bits(chosen.eval(q)), bits(again.eval(q)));
        }
    }

    #[test]
    fn a_motion_and_a_scaling_carry_the_curve_along(
        points in polyline(),
        way in ways(),
        degree in 1..=5_usize,
        angle in 0.0_f64..6.2,
        scale in 0.25_f64..4.0,
        shift in prop::array::uniform2(-5.0_f64..5.0),
        us in prop::collection::vec(0.0_f64..=1.0, 1..20),
    ) {
        let (sin, cos) = angle.sin_cos();
        let moved = |p: [f64; 2]| [
            scale * (cos * p[0] - sin * p[1]) + shift[0],
            scale * (sin * p[0] + cos * p[1]) + shift[1],
        ];
        let image: Vec<[f64; 2]> = points.iter().map(|p| moved(*p)).collect();
        let curve = BSplineCurve::new(&points, degree, Bc::NotAKnot, way).expect("valid");
        let other = BSplineCurve::new(&image, degree, Bc::NotAKnot, way).expect("valid");
        // The parameters go by ratios of distances, which a motion and a
        // uniform scaling keep.
        for (a, b) in curve.parameters().iter().zip(other.parameters()) {
            prop_assert!((a - b).abs() <= 1e-12, "{} and {}", a, b);
        }
        let tol = 1e-9 * (size(&curve) * scale + 5.0);
        for u in &us {
            let want = moved(curve.eval(u).expect("inside"));
            let got = other.eval(u).expect("inside");
            prop_assert!((got[0] - want[0]).abs() <= tol && (got[1] - want[1]).abs() <= tol, "at {}: {:?} for {:?}", u, got, want);
        }
    }

    #[test]
    fn reversed_points_give_the_curve_backward(
        points in polyline(),
        way in ways(),
        degree in 1..=5_usize,
        us in prop::collection::vec(0.0_f64..=1.0, 1..20),
    ) {
        let reversed: Vec<[f64; 2]> = points.iter().rev().copied().collect();
        let curve = BSplineCurve::new(&points, degree, Bc::NotAKnot, way).expect("valid");
        let back = BSplineCurve::new(&reversed, degree, Bc::NotAKnot, way).expect("valid");
        let n = points.len();
        for (i, u) in curve.parameters().iter().enumerate() {
            let mirrored = 1.0 - back.parameters()[n - 1 - i];
            prop_assert!((u - mirrored).abs() <= 1e-12, "{} and {}", u, mirrored);
        }
        let tol = 1e-8 * size(&curve);
        for u in &us {
            let (there, here) = (curve.eval(u).expect("inside"), back.eval(&(1.0 - u)).expect("inside"));
            prop_assert!((there[0] - here[0]).abs() <= tol && (there[1] - here[1]).abs() <= tol, "at {}: {:?} and {:?}", u, there, here);
        }
    }

    #[test]
    fn evaluation_is_that_of_the_b_spline(
        points in polyline(),
        way in ways(),
        degree in 1..=5_usize,
        qs in prop::collection::vec(-1.0_f64..3.0, 1..60),
    ) {
        for mode in [Extrapolate::Forbid, Extrapolate::Clamp, Extrapolate::Extend, Extrapolate::Periodic] {
            let curve = BSplineCurve::new(&points, degree, Bc::NotAKnot, way).expect("valid").with_extrapolate(mode);
            let b = curve.as_bspline();
            prop_assert_eq!(b.extrapolate(), mode);
            let mut hint = Hint::new();
            for q in &qs {
                prop_assert_eq!(bits(curve.eval(q)), bits(b.eval(q)), "{:?} u = {}", mode, q);
                prop_assert_eq!(bits(curve.eval_hinted(*q, &mut hint)), bits(b.eval(q)));
            }
            let (mut out, mut want) = (vec![[0.0; 2]; qs.len()], vec![[0.0; 2]; qs.len()]);
            prop_assert_eq!(curve.eval_batch_into(&qs, &mut out, [-7.0; 2]), Ok(()));
            prop_assert_eq!(b.eval_batch_into(&qs, &mut want, [-7.0; 2]), Ok(()));
            for (got, want) in out.iter().zip(&want) {
                prop_assert_eq!(got.map(f64::to_bits), want.map(f64::to_bits));
            }
        }
    }
}
