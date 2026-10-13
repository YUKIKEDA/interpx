//! `NearestNDInterpolator`: SciPy's values, a scan of all points as the
//! reference, points given more than once, points at the same distance, the
//! bound on the distance, rescaled coordinates, the input checks and their
//! order, and the properties.

#![cfg(feature = "scattered")]

mod common;

use common::{load, reals, Real};
use interpx::{InterpolateError as E, Interpolator, NearestNDInterpolator as Nearest, ParamReason};
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

fn bits(v: Option<f64>) -> Option<u64> {
    v.map(f64::to_bits)
}

/// The change of coordinates of SciPy's `rescale`: the mean of the points on
/// each axis is taken off, and the range divides (1 where the range is 0).
fn scaled<const D: usize>(
    points: &[[f64; D]],
    queries: &[[f64; D]],
) -> (Vec<[f64; D]>, Vec<[f64; D]>) {
    let mut offset = [0.0; D];
    let mut scale = [1.0; D];
    for d in 0..D {
        let column: Vec<f64> = points.iter().map(|p| p[d]).collect();
        offset[d] = column.iter().sum::<f64>() / column.len() as f64;
        let least = column.iter().copied().fold(f64::INFINITY, f64::min);
        let most = column.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        if most > least {
            scale[d] = most - least;
        }
    }
    let change = |p: &[f64; D]| core::array::from_fn(|d| (p[d] - offset[d]) / scale[d]);
    (
        points.iter().map(change).collect(),
        queries.iter().map(change).collect(),
    )
}

/// The reference: a scan of all points. Returns the distance to the nearest
/// point, and the points at that distance (to the rounding of a square).
fn scan<const D: usize>(points: &[[f64; D]], query: &[f64; D]) -> (f64, Vec<usize>) {
    let squared = |p: &[f64; D]| {
        (0..D)
            .map(|d| (p[d] - query[d]) * (p[d] - query[d]))
            .sum::<f64>()
    };
    let least = points.iter().map(squared).fold(f64::INFINITY, f64::min);
    let at_least = (0..points.len())
        .filter(|i| squared(&points[*i]) <= least * (1.0 + 1e-12))
        .collect();
    (least.sqrt(), at_least)
}

// Marked as test code so that `expect` is allowed here as in a `#[test]` function
// (clippy.toml allows it in tests only).
#[cfg(test)]
fn check<const D: usize>(case: &Json) -> [usize; 3] {
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
    let bound = case["options"]["distance_upper_bound"].as_f64();
    let build = |values: &[f64]| {
        let built = if rescale {
            Nearest::rescaled(&points, values)
        } else {
            Nearest::new(&points, values)
        }
        .expect("valid points");
        assert_eq!(built.distance_bound(), None, "{name}: no bound at first");
        match bound {
            Some(bound) => built.with_distance_bound(bound),
            None => built,
        }
    };
    let nearest = build(&values);
    assert_eq!(nearest.len(), points.len(), "{name}");
    let expected = nums(&case["expected"]);
    let tied: Vec<bool> = serde_json::from_value(case["tied"].clone()).unwrap_or_default();
    // The distances as the interpolant measures them.
    let (measured, asked) = if rescale {
        scaled(&points, &queries)
    } else {
        (points.clone(), queries.clone())
    };
    let mut batch = vec![0.0; queries.len()];
    assert_eq!(
        nearest.eval_batch_into(&queries, &mut batch, f64::NAN),
        Ok(())
    );
    let (mut without, mut several) = (0, 0);
    for (k, (query, want)) in queries.iter().zip(&expected).enumerate() {
        let what = format!("{name} query {k} {query:?}");
        let got = nearest.eval(query);
        assert_eq!(
            bits(Some(batch[k])),
            bits(Some(got.unwrap_or(f64::NAN))),
            "{what}"
        );
        let (distance, at_least) = scan(&measured, &asked[k]);
        if want.is_nan() {
            // Only a bound leaves a query without a value.
            assert_eq!(got, None, "{what}");
            assert!(
                bound.is_some_and(|b| distance >= b),
                "{what}: at {distance}"
            );
            without += 1;
        } else if tied.get(k).copied().unwrap_or(false) {
            // Several points at the same distance: the value of one of them.
            assert!(at_least.len() > 1, "{what}: marked as tied");
            let got = got.expect("a value");
            assert!(
                at_least
                    .iter()
                    .any(|i| values[*i].to_bits() == got.to_bits()),
                "{what}: {got} is the value of no nearest point"
            );
            several += 1;
        } else {
            // A value of the data, bit for bit.
            assert_eq!(bits(got), bits(Some(*want)), "{what}");
        }
    }

    // Two values at every point: the same point is found for both.
    if let Some(columns) = case.get("values_columns") {
        let pairs: Vec<[f64; 2]> = rows(columns).iter().map(|r| [r[0], r[1]]).collect();
        let pair = Nearest::new(&points, &pairs).expect("valid points");
        for (query, want) in queries.iter().zip(rows(&case["expected_columns"])) {
            let got = pair.eval(query).expect("a value");
            assert_eq!(
                got.map(f64::to_bits),
                [want[0], want[1]].map(f64::to_bits),
                "{name}"
            );
        }
    }
    [queries.len(), without, several]
}

// --- Golden data -------------------------------------------------------------

#[test]
fn scipy_nearest_neighbors_evaluate_to_the_golden_values() {
    let file = load::<Json>("scattered").expect("golden file scattered");
    let mut cases = [0; 5];
    let mut totals = [0; 3];
    let (mut rescaled, mut bounded, mut columns) = (0, 0, 0);
    for case in &file.cases {
        if !case["name"].as_str().expect("name").starts_with("nearest/") {
            continue;
        }
        let dimension = rows(&case["points"])[0].len();
        let counts = match dimension {
            2 => check::<2>(case),
            3 => check::<3>(case),
            4 => check::<4>(case),
            other => panic!("points in {other} dimensions"),
        };
        cases[dimension] += 1;
        for (total, count) in totals.iter_mut().zip(counts) {
            *total += count;
        }
        rescaled += usize::from(case["options"]["rescale"].as_bool().unwrap_or(false));
        bounded += usize::from(case["options"]["distance_upper_bound"].is_number());
        columns += usize::from(case.get("values_columns").is_some());
    }
    assert_eq!(cases, [0, 0, 19, 4, 1]);
    assert_eq!((rescaled, bounded, columns), (3, 3, 1));
    // Queries, those the bound leaves without a value, and those with several
    // points at the same distance.
    assert_eq!(totals, [1245, 75, 12]);
}

// --- Input checks ---------------------------------------------------------------

#[test]
fn the_constructors_report_each_fault_in_the_order_of_the_design() {
    let points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
    type Build = fn(&[[f64; 2]], &[f64]) -> Result<Nearest<2>, E>;
    let constructors: [Build; 2] = [Nearest::new, Nearest::rescaled];
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
            err(&[], &[]),
            Some(E::InsufficientPoints {
                required: 1,
                actual: 0
            })
        );
        // The length before the count.
        assert_eq!(
            err(&[], &[1.0]),
            Some(E::LengthMismatch {
                expected: 0,
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
        // An infinite coordinate, of either sign, is refused like a NaN.
        for infinite in [f64::INFINITY, f64::NEG_INFINITY] {
            let far = [[0.0, 0.0], [1.0, 0.0], [0.0, infinite]];
            assert_eq!(err(&far, &[1.0, 2.0, 3.0]), Some(E::NonFinite { index: 2 }));
        }
        // One point is enough.
        let one = build(&points[..1], &[4.0]).expect("one point");
        assert_eq!(bits(one.eval(&[100.0, -100.0])), bits(Some(4.0)));
    }
    // Rescaled: a mean or a range that overflows is refused, and not left to
    // put every point at one place. The sum passes the largest number at the
    // second point; the range of the second axis is too wide, and the point
    // with its largest coordinate is named.
    let huge = 1.7e308;
    let same = [[huge, 0.0], [huge, 1.0], [huge, 2.0]];
    assert_eq!(
        Nearest::rescaled(&same, &[1.0, 2.0, 3.0]).err(),
        Some(E::NonFinite { index: 1 })
    );
    let wide = [[0.0, -huge], [1.0, 0.0], [2.0, huge], [3.0, 1.0]];
    assert_eq!(
        Nearest::rescaled(&wide, &[1.0, 2.0, 3.0, 4.0]).err(),
        Some(E::NonFinite { index: 2 })
    );
    // As the numbers stand both are sets of points like any other.
    let plain = Nearest::new(&same, &[1.0, 2.0, 3.0]).expect("finite points");
    assert_eq!(bits(plain.eval(&[huge, 1.2])), bits(Some(2.0)));
    assert!(Nearest::new(&wide, &[1.0, 2.0, 3.0, 4.0]).is_ok());
    // Large, and still with a mean and a range.
    let large = [[8e307, 0.0], [-8e307, 1.0]];
    let rescaled = Nearest::rescaled(&large, &[1.0, 2.0]).expect("finite mean and range");
    assert_eq!(bits(rescaled.eval(&[7e307, 0.1])), bits(Some(1.0)));

    // More points than the tree can number. Points without coordinates and
    // values without columns take no memory, so the count can be reached.
    #[cfg(target_pointer_width = "64")]
    {
        let count = u32::MAX as usize + 1;
        let (points, values) = (vec![[0.0_f64; 0]; count], vec![[0.0_f64; 0]; count]);
        let too_many = E::InvalidParameter {
            name: "points",
            reason: ParamReason::TooLarge {
                max: u32::MAX as usize,
            },
        };
        assert_eq!(Nearest::new(&points, &values).err(), Some(too_many));
        // The length still comes first.
        assert_eq!(
            Nearest::new(&points, &values[1..]).err(),
            Some(E::LengthMismatch {
                expected: count,
                actual: count - 1
            })
        );
    }
}

// --- The rules ------------------------------------------------------------------

#[test]
fn values_handed_over_give_the_interpolant_of_values_copied() {
    // Points given twice, values of two columns, a query on a point and ones
    // between: the same answers, plain and rescaled.
    let points = [
        [0.0, 0.0],
        [1000.0, 0.0],
        [0.0, 1.0],
        [1000.0, 0.0],
        [400.0, 0.5],
    ];
    let values = [
        [1.0, -1.0],
        [2.0, -2.0],
        [3.0, -3.0],
        [9.0, -9.0],
        [5.0, -5.0],
    ];
    // The first query has another nearest point once the axes are rescaled.
    let queries = [
        [300.0, 0.0],
        [400.0, 0.3],
        [1000.0, 0.0],
        [600.0, 0.9],
        [f64::NAN, 0.0],
    ];
    let pairs = [
        (
            Nearest::new(&points, &values),
            Nearest::from_owned(&points, values.to_vec()),
        ),
        (
            Nearest::rescaled(&points, &values),
            Nearest::rescaled_from_owned(&points, values.to_vec()),
        ),
    ];
    let mut first = Vec::new();
    for (copied, moved) in pairs {
        let (copied, moved) = (copied.expect("valid"), moved.expect("valid"));
        assert_eq!(moved.len(), copied.len());
        assert_eq!(moved.values(), copied.values());
        assert_eq!(moved.eval_batch(&queries), copied.eval_batch(&queries));
        first.push(moved.eval(&queries[0]));
    }
    assert_eq!(first, [Some([5.0, -5.0]), Some([1.0, -1.0])]);
    // The same faults, in the same order.
    let short = Nearest::from_owned(&points, vec![[0.0, 0.0]; 4]);
    assert_eq!(
        short.err(),
        Some(E::LengthMismatch {
            expected: 5,
            actual: 4
        })
    );
    let mut bad = values.to_vec();
    bad[3][1] = f64::INFINITY;
    assert_eq!(
        Nearest::from_owned(&points, bad.clone()).err(),
        Some(E::NonFinite { index: 3 })
    );
    assert_eq!(
        Nearest::rescaled_from_owned(&points, bad).err(),
        Some(E::NonFinite { index: 3 })
    );
    let none: [[f64; 2]; 0] = [];
    let empty = Nearest::from_owned(&none, Vec::<f64>::new());
    assert_eq!(
        empty.err(),
        Some(E::InsufficientPoints {
            required: 1,
            actual: 0
        })
    );
}

#[test]
fn points_given_more_than_once_count_once_with_the_value_given_first() {
    // The same point at places 1, 3, and 4, and one more twice with zeros of
    // both signs.
    let points = [
        [5.0, 5.0],
        [1.0, 2.0],
        [0.0, -0.0],
        [1.0, 2.0],
        [1.0, 2.0],
        [-0.0, 0.0],
    ];
    let values = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0];
    for build in [Nearest::new, Nearest::rescaled] {
        let nearest = build(&points, &values).expect("valid");
        assert_eq!(bits(nearest.eval(&[1.0, 2.0])), bits(Some(10.0)));
        assert_eq!(bits(nearest.eval(&[1.2, 2.1])), bits(Some(10.0)));
        assert_eq!(bits(nearest.eval(&[0.1, -0.1])), bits(Some(20.0)));
        assert_eq!(bits(nearest.eval(&[9.0, 9.0])), bits(Some(0.0)));
        // All the values are kept, also those no query returns.
        assert_eq!(nearest.values(), &values);
        assert_eq!(nearest.len(), 6);
    }
    // The two zeros are one coordinate. In six dimensions the origin can be
    // written in 64 ways, more than a leaf of the tree holds: were they 64
    // points, a query would get one of them by the layout of the tree. They
    // are one point, with the value given first.
    let origins: Vec<[f64; 6]> = (0..64_u32)
        .map(|signs| core::array::from_fn(|d| if signs >> d & 1 == 1 { -0.0 } else { 0.0 }))
        .collect();
    let numbered: Vec<f64> = (0..64).map(f64::from).collect();
    for turn in [0, 1, 17, 40, 63] {
        let mut points = origins.clone();
        points.rotate_left(turn);
        let nearest = Nearest::new(&points, &numbered).expect("valid");
        for query in [
            [0.0; 6],
            [-0.0; 6],
            [0.3, -0.2, 0.0, 1.0, -1.0, 0.5],
            [9.0; 6],
        ] {
            assert_eq!(
                bits(nearest.eval(&query)),
                bits(Some(0.0)),
                "{turn} {query:?}"
            );
        }
    }
    // Enough equal points to fill more than one leaf of the tree (32).
    let many = vec![[0.5, 0.5, 0.5]; 100];
    let numbered: Vec<f64> = (0..100).map(f64::from).collect();
    let nearest = Nearest::new(&many, &numbered).expect("valid");
    assert_eq!(bits(nearest.eval(&[0.0, 1.0, 2.0])), bits(Some(0.0)));
}

#[test]
fn a_bound_leaves_a_query_with_a_value_only_when_its_point_is_nearer() {
    let points = [[0.0, 0.0], [10.0, 0.0]];
    let nearest = Nearest::new(&points, &[1.0, 2.0]).expect("valid");
    let near = nearest.clone().with_distance_bound(5.0);
    // The clone has the bound; the first is as it was.
    assert_eq!(nearest.distance_bound(), None);
    assert_eq!(
        near.distance_bound().map(f64::to_bits),
        Some(5.0_f64.to_bits())
    );
    // 3-4-5: at the bound exactly there is no value, a little nearer there is.
    assert_eq!(near.eval(&[3.0, 4.0]), None);
    assert_eq!(bits(near.eval(&[3.0, 3.999_999])), bits(Some(1.0)));
    assert_eq!(near.eval(&[3.0, 4.000_001]), None);
    assert_eq!(bits(near.eval(&[10.0, 0.0])), bits(Some(2.0)));
    assert_eq!(bits(nearest.eval(&[3.0, 4.0])), bits(Some(1.0)));
    assert_eq!(
        bits(Some(near.eval_or(&[3.0, 4.0], -7.0))),
        bits(Some(-7.0))
    );
    // No distance is under zero, or under NaN; every distance is under infinity.
    for none in [0.0, -1.0, f64::NAN] {
        assert_eq!(
            nearest.clone().with_distance_bound(none).eval(&[0.0, 0.0]),
            None
        );
    }
    let open = nearest.clone().with_distance_bound(f64::INFINITY);
    assert_eq!(
        bits(open.eval(&[1e100, 0.0])),
        bits(nearest.eval(&[1e100, 0.0]))
    );
    assert!(open.eval(&[1e100, 0.0]).is_some());
    // Rescaled, the bound is on the distance after the change: the axes have
    // the ranges 10 and 1, so a step of 4 along the first is 0.4.
    let wide = [[0.0, 0.0], [10.0, 0.0], [0.0, 1.0]];
    let rescaled = Nearest::rescaled(&wide, &[1.0, 2.0, 3.0])
        .expect("valid")
        .with_distance_bound(0.5);
    assert_eq!(bits(rescaled.eval(&[4.0, 0.0])), bits(Some(1.0)));
    assert_eq!(rescaled.eval(&[4.0, -0.4]), None);
}

#[test]
fn a_query_that_is_not_finite_or_too_far_to_measure_has_no_value() {
    let points = [[0.0, 0.0], [1.0, 1.0]];
    for build in [Nearest::new, Nearest::rescaled] {
        let nearest = build(&points, &[1.0, 2.0]).expect("valid");
        for query in [
            [f64::NAN, 0.0],
            [0.0, f64::NAN],
            [f64::INFINITY, 0.0],
            [0.0, f64::NEG_INFINITY],
            // The square of the distance does not fit.
            [1e200, 0.0],
            [-1e160, 1e160],
        ] {
            assert_eq!(nearest.eval(&query), None, "{query:?}");
            assert_eq!(bits(Some(nearest.eval_or(&query, -7.0))), bits(Some(-7.0)));
        }
        // Far, and still measured (both points are as far, to the rounding).
        assert!(nearest.eval(&[1e150, 1e150]).is_some());
    }
}

#[test]
fn a_batch_gives_the_values_of_eval_and_fills_the_rest() {
    let points = [[0.0], [1.0], [3.0]];
    let nearest = Nearest::new(&points, &[10.0, 20.0, 30.0])
        .expect("valid")
        .with_distance_bound(0.75);
    let queries = [[0.2], [0.8], [2.0], [f64::NAN], [3.5], [0.5]];
    let mut out = [0.0; 6];
    assert_eq!(nearest.eval_batch_into(&queries, &mut out, -7.0), Ok(()));
    // 0.5 is as far from 0 as from 1: either value.
    assert!(out[5].to_bits() == 10.0_f64.to_bits() || out[5].to_bits() == 20.0_f64.to_bits());
    assert_eq!(
        out[..5].iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
        [10.0_f64, 20.0, -7.0, -7.0, 30.0].map(f64::to_bits)
    );
    let mut short = [1.5; 5];
    assert_eq!(
        nearest.eval_batch_into(&queries, &mut short, -7.0),
        Err(E::LengthMismatch {
            expected: 6,
            actual: 5
        })
    );
    assert!(short.iter().all(|v| v.to_bits() == 1.5_f64.to_bits()));
    assert_eq!(nearest.eval_batch(&queries[..2]).len(), 2);
}

#[test]
fn no_dimensions_columns_and_thread_safety() {
    // Without coordinates every point is the same point: the first value.
    let point: Nearest<0> = Nearest::new(&[[], []], &[4.5, 5.5]).expect("valid");
    assert_eq!(bits(point.eval(&[])), bits(Some(4.5)));
    // Arrays for values.
    let pairs = [[1.0, -1.0], [2.0, -2.0]];
    let both: Nearest<1, [f64; 2]> = Nearest::new(&[[0.0], [1.0]], &pairs).expect("valid");
    assert_eq!(
        both.eval(&[0.9]).map(|v| v.map(f64::to_bits)),
        Some([2.0_f64, -2.0].map(f64::to_bits))
    );
    let nan = Nearest::new(&[[0.0], [1.0]], &[[1.0, 1.0], [2.0, f64::NAN]]);
    assert_eq!(nan.err(), Some(E::NonFinite { index: 1 }));

    fn assert_send_sync<T: Send + Sync>() {}
    assert_send_sync::<Nearest<3>>();
    assert_send_sync::<Nearest<2, [f64; 4]>>();
    // Clones share the tree and answer alike, from another thread too.
    let nearest = Nearest::new(&[[0.0, 0.0], [1.0, 1.0]], &[1.0, 2.0]).expect("valid");
    let other = nearest.clone();
    let from_thread = std::thread::spawn(move || other.eval(&[0.9, 0.8]))
        .join()
        .expect("joins");
    assert_eq!(bits(from_thread), bits(nearest.eval(&[0.9, 0.8])));
    assert_eq!(bits(from_thread), bits(Some(2.0)));
}

// --- Properties ----------------------------------------------------------------

/// Points in three dimensions with coordinates from a few values, so that
/// points repeat and distances tie; and queries, half of them from the same
/// values (ties) and half anywhere.
type Case = (Vec<[f64; 3]>, Vec<[f64; 3]>);

fn cases() -> impl Strategy<Value = Case> {
    let coordinate = prop_oneof![
        prop::sample::select(vec![-1.0, 0.0, 0.5, 1.0, 2.0]),
        -1.5_f64..2.5,
    ];
    (
        prop::collection::vec(prop::array::uniform3(coordinate.clone()), 1..120),
        prop::collection::vec(prop::array::uniform3(coordinate), 1..40),
    )
}

proptest! {
    #[test]
    fn the_value_is_that_of_a_nearest_point_and_of_the_first_of_equal_points(
        (points, queries) in cases(),
        rescale in any::<bool>(),
        bound in prop::option::of(0.0_f64..3.0),
    ) {
        // The value of a point is its place, so that the point shows.
        let values: Vec<f64> = (0..points.len()).map(|i| i as f64).collect();
        let built = if rescale { Nearest::rescaled(&points, &values) } else { Nearest::new(&points, &values) };
        let mut nearest = built.expect("valid");
        if let Some(bound) = bound {
            nearest = nearest.with_distance_bound(bound);
        }
        let (measured, asked) = if rescale { scaled(&points, &queries) } else { (points.clone(), queries.clone()) };
        let mut out = vec![0.0; queries.len()];
        prop_assert_eq!(nearest.eval_batch_into(&queries, &mut out, -7.0), Ok(()));
        for ((query, asked), slot) in queries.iter().zip(&asked).zip(&out) {
            let got = nearest.eval(query);
            prop_assert_eq!(bits(Some(*slot)), bits(Some(got.unwrap_or(-7.0))));
            let (distance, at_least) = scan(&measured, asked);
            match got {
                None => {
                    let bound = bound.expect("only a bound leaves a query without a value");
                    prop_assert!(distance >= bound * (1.0 - 1e-12), "{:?}: no value at {} under {}", query, distance, bound);
                }
                Some(value) => {
                    if let Some(bound) = bound {
                        prop_assert!(distance <= bound * (1.0 + 1e-12), "{:?}: a value at {} over {}", query, distance, bound);
                    }
                    let found = value as usize;
                    prop_assert!(at_least.contains(&found), "{:?}: point {} is not nearest ({:?})", query, found, at_least);
                    // Of the points with its coordinates, it is the first.
                    let same = |a: &[f64; 3], b: &[f64; 3]| (0..3).all(|d| a[d].total_cmp(&b[d]).is_eq() || (a[d] + 0.0).to_bits() == (b[d] + 0.0).to_bits());
                    let first = points.iter().position(|p| same(p, &points[found])).expect("itself");
                    prop_assert_eq!(found, first, "{:?}", query);
                }
            }
        }
    }

    #[test]
    fn in_one_dimension_the_nearest_point_is_the_nearest_number(
        xs in prop::collection::vec(-10.0_f64..10.0, 1..60),
        queries in prop::collection::vec(-12.0_f64..12.0, 1..40),
    ) {
        let points: Vec<[f64; 1]> = xs.iter().map(|x| [*x]).collect();
        let nearest = Nearest::new(&points, &xs).expect("valid");
        for q in queries {
            let got = nearest.eval(&[q]).expect("a value");
            let least = xs.iter().map(|x| (x - q).abs()).fold(f64::INFINITY, f64::min);
            prop_assert!((got - q).abs() <= least * (1.0 + 1e-12), "{}: {} is {} away, the nearest {}", q, got, (got - q).abs(), least);
        }
    }

    #[test]
    fn rescaled_is_new_on_the_coordinates_brought_to_one_range(
        (points, queries) in cases(),
        stretch in prop::array::uniform3(0.01_f64..1000.0),
    ) {
        // Other units on each axis.
        let stretched = |p: &[f64; 3]| core::array::from_fn(|d| p[d] * stretch[d]);
        let points: Vec<[f64; 3]> = points.iter().map(stretched).collect();
        let queries: Vec<[f64; 3]> = queries.iter().map(stretched).collect();
        let values: Vec<f64> = (0..points.len()).map(|i| i as f64).collect();
        let rescaled = Nearest::rescaled(&points, &values).expect("valid");
        let (measured, asked) = scaled(&points, &queries);
        let plain = Nearest::new(&measured, &values).expect("valid");
        for (query, asked) in queries.iter().zip(&asked) {
            prop_assert_eq!(bits(rescaled.eval(query)), bits(plain.eval(asked)), "{:?}", query);
        }
    }
}
