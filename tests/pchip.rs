//! `PchipInterpolator`: the golden data, the input checks through the
//! constructor, the rules for the slopes, and the properties that make it
//! shape-preserving.

mod common;

use common::{assert_close, load, reals, Real, DIRECT_F64};
use interpx::{Extrapolate, Hint, InterpolateError as E, Interpolator, PchipInterpolator};
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

/// The slope at data point `i`, read from the piecewise polynomial: the linear
/// coefficient of the piece that starts there, or for the last point the slope
/// at the end of the last piece.
fn slope_at(pchip: &PchipInterpolator<f64>, i: usize) -> f64 {
    let poly = pchip.as_ppoly();
    let (x, c) = (poly.breakpoints(), poly.coefficients());
    if i + 1 < x.len() {
        c[4 * i + 2]
    } else {
        let (h, c) = (x[i] - x[i - 1], &c[4 * (i - 1)..4 * i]);
        (3.0 * c[0] * h + 2.0 * c[1]) * h + c[2]
    }
}

// --- Golden data -------------------------------------------------------------

/// A golden column as a function of the query.
type Column<'a> = &'a dyn Fn(f64) -> Option<f64>;

#[test]
fn scipy_pchip_interpolants_are_rebuilt_from_their_data() {
    let file = load::<Json>("oned").expect("golden file oned");
    let (mut cases, mut small, mut forbid) = (0, 0, 0);
    for case in &file.cases {
        let name = case["name"].as_str().expect("name");
        if case["reference"] != "scipy.interpolate.PchipInterpolator"
            || case.get("coefficients").is_none()
        {
            continue;
        }
        let (x, y) = (nums(&case["x"]), nums(&case["y"]));
        let built = PchipInterpolator::new(&x, &y).expect("valid data");
        assert_eq!(
            built.extrapolate(),
            Extrapolate::Extend,
            "{name}: SciPy's default"
        );
        let extrapolates = case["options"]["extrapolate"] == Json::Bool(true);
        let pchip = if extrapolates {
            built
        } else {
            forbid += 1;
            built.with_extrapolate(Extrapolate::Forbid)
        };

        // A coefficient is compared by what it adds to the value at the far end
        // of its interval (`tests/cubic_spline.rs` says why): SciPy stores
        // `c[power][interval]`.
        let scipy = rows(&case["coefficients"]);
        let poly = pchip.as_ppoly();
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
            ("expected", &|q| pchip.eval(&q)),
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
        small += usize::from(x.len() <= 4);
    }
    // 9 on the standard data and one of them without extrapolation, then 9
    // chosen for the rules, 7 of them with four points or fewer.
    assert_eq!((cases, small, forbid), (19, 7, 1));
}

#[test]
fn scipy_pchip_of_two_columns_is_rebuilt_column_by_column() {
    let file = load::<Json>("oned").expect("golden file oned");
    let mut cases = 0;
    for case in &file.cases {
        let name = case["name"].as_str().expect("name");
        if case["reference"] != "scipy.interpolate.PchipInterpolator"
            || case.get("coefficients").is_some()
        {
            continue;
        }
        let y: Vec<[f64; 2]> = rows(&case["y"])
            .iter()
            .map(|r| r.as_slice().try_into().expect("two columns"))
            .collect();
        let pchip = PchipInterpolator::new(&nums(&case["x"]), &y).expect("valid data");
        let poly = pchip.as_ppoly();
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

#[test]
fn the_constructor_reports_each_fault() {
    assert!(PchipInterpolator::new(&[0.0, 1.0], &[1.0, 2.0]).is_ok());
    let cases: [(&[f64], &[f64], E); 8] = [
        (
            &[0.0, 1.0, 2.0],
            &[1.0, 2.0],
            E::LengthMismatch {
                expected: 3,
                actual: 2,
            },
        ),
        (
            &[0.0],
            &[1.0],
            E::InsufficientPoints {
                required: 2,
                actual: 1,
            },
        ),
        (
            &[0.0, f64::NAN, 2.0],
            &[1.0, 2.0, 3.0],
            E::NonFinite { index: 1 },
        ),
        (
            &[0.0, 1.0, 2.0],
            &[1.0, 2.0, f64::INFINITY],
            E::NonFinite { index: 2 },
        ),
        (
            &[0.0, 2.0, 1.0],
            &[1.0, 2.0, 3.0],
            E::NotSorted { index: 2 },
        ),
        (
            &[0.0, 1.0, 1.0],
            &[1.0, 2.0, 3.0],
            E::DuplicatePoints {
                first: 1,
                second: 2,
            },
        ),
        // Finite values whose difference is not: the secant to point 2.
        (
            &[0.0, 1.0, 2.0],
            &[0.0, -1e308, 1e308],
            E::NonFinite { index: 2 },
        ),
        // The same in the first interval.
        (
            &[0.0, 1.0, 2.0],
            &[1e308, -1e308, 0.0],
            E::NonFinite { index: 1 },
        ),
    ];
    for (x, y, want) in cases {
        assert_eq!(PchipInterpolator::new(x, y).err(), Some(want));
    }
    // The order of the checks: lengths, the count, NaN, then the order of `x`.
    assert_eq!(
        PchipInterpolator::new(&[0.0], &[1.0, 2.0]).err(),
        Some(E::LengthMismatch {
            expected: 1,
            actual: 2
        })
    );
    assert_eq!(
        PchipInterpolator::new(&[1.0, 0.0], &[f64::NAN, 1.0]).err(),
        Some(E::NonFinite { index: 0 })
    );
    assert_eq!(
        PchipInterpolator::new(&[1.0, 0.0], &[-1e308, 1e308]).err(),
        Some(E::NotSorted { index: 1 })
    );
}

#[test]
fn an_overflowing_coefficient_names_its_interval() {
    // Values that swing by 1e300 on intervals 1e-3 wide: the secants are finite,
    // the cubic coefficient of the second interval, where the slopes at the two
    // ends are zero (the data turns at both), is the secant over the width
    // squared.
    let x = [0.0, 1.0, 1.001, 1.002, 2.0];
    let y = [0.0, 1.0, 1e300, 0.0, 1.0];
    assert_eq!(
        PchipInterpolator::new(&x, &y).err(),
        Some(E::NonFinite { index: 1 })
    );
}

// --- The rules for the slopes ----------------------------------------------------

#[test]
fn two_points_give_the_line_through_them() {
    let pchip = PchipInterpolator::new(&[1.0, 3.0], &[2.0, 8.0]).expect("valid");
    let c = pchip.as_ppoly().coefficients();
    assert_eq!(
        [c[0], c[1], c[2], c[3]].map(f64::to_bits),
        [0.0_f64, 0.0, 3.0, 2.0].map(f64::to_bits)
    );
    assert_eq!(bits(pchip.eval(&2.0)), bits(Some(5.0)));
    assert_eq!(bits(pchip.eval(&5.0)), bits(Some(14.0)), "extended");
}

#[test]
fn the_slope_is_zero_where_the_data_turns_or_is_flat() {
    // Flat, a rise, a peak, a fall onto a flat run, a rise.
    let x = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
    let y = [0.0, 0.0, 1.0, 3.0, 2.0, 2.0, 2.0, 5.0];
    let pchip = PchipInterpolator::new(&x, &y).expect("valid");
    // Beside a flat interval (1, 4, 5, 6) and at the peak (3).
    for i in [1, 3, 4, 5, 6] {
        assert_eq!(slope_at(&pchip, i).abs().to_bits(), 0, "point {i}");
    }
    // Between two rising secants, 1 and 2, of equal widths: their harmonic mean.
    assert_close(slope_at(&pchip, 2), 4.0 / 3.0, DIRECT_F64, "point 2");
    // A flat run stays flat, exactly.
    for q in [0.25, 0.5, 4.5, 5.25, 5.999] {
        let want = if q < 1.0 { 0.0 } else { 2.0 };
        assert_eq!(bits(pchip.eval(&q)), bits(Some(want)), "x = {q}");
    }
}

#[test]
fn the_inner_slope_is_the_harmonic_mean_weighted_by_the_widths() {
    // Widths 1 and 3, secants 2 and 4: w1 = 2 * 3 + 1 = 7, w2 = 3 + 2 * 1 = 5,
    // and the slope is 12 / (7 / 2 + 5 / 4).
    let pchip = PchipInterpolator::new(&[0.0, 1.0, 4.0], &[0.0, 2.0, 14.0]).expect("valid");
    assert_close(
        slope_at(&pchip, 1),
        12.0 / (3.5 + 1.25),
        DIRECT_F64,
        "the middle",
    );
    // Falling data gives the same with the other sign.
    let falling = PchipInterpolator::new(&[0.0, 1.0, 4.0], &[0.0, -2.0, -14.0]).expect("valid");
    assert_close(
        slope_at(&falling, 1),
        -12.0 / (3.5 + 1.25),
        DIRECT_F64,
        "falling",
    );
}

#[test]
fn the_end_slope_is_the_three_point_estimate_cut_to_keep_the_shape() {
    let ends = |y: [f64; 4]| {
        let pchip = PchipInterpolator::new(&[0.0, 1.0, 2.0, 3.0], &y).expect("valid");
        (slope_at(&pchip, 0), slope_at(&pchip, 3))
    };
    // Secants 1 and 2 at the left: (3 * 1 - 2) / 2 = 0.5, kept as it is. At the
    // right 4 and 2: (3 * 4 - 2) / 2 = 5.
    let (left, right) = ends([0.0, 1.0, 3.0, 7.0]);
    assert_close(left, 0.5, DIRECT_F64, "left, plain");
    assert_close(right, 5.0, DIRECT_F64, "right, plain");
    // Secants 1 and 9: (3 - 9) / 2 = -3 points against the first secant: zero.
    let (left, _) = ends([0.0, 1.0, 10.0, 11.0]);
    assert_eq!(left.abs().to_bits(), 0, "left, cut to zero");
    let (_, right) = ends([11.0, 10.0, 1.0, 0.0]);
    assert_eq!(right.abs().to_bits(), 0, "right, cut to zero");
    // Secants 1 and -11: (3 + 11) / 2 = 7 is over three times the first secant,
    // and the second has the other sign: 3.
    let (left, _) = ends([0.0, 1.0, -10.0, -9.0]);
    assert_eq!(left.to_bits(), 3.0_f64.to_bits(), "left, limited");
    let (_, right) = ends([-9.0, -10.0, 1.0, 0.0]);
    assert_eq!(right.to_bits(), (-3.0_f64).to_bits(), "right, limited");
    // Secants 1 and -1: 2 is under three times the first secant, so it stays,
    // although the second secant has the other sign.
    let (left, _) = ends([0.0, 1.0, 0.0, 1.0]);
    assert_close(left, 2.0, DIRECT_F64, "left, under the limit");
    // A flat first interval: the estimate, -m1 / 2, has a sign and the secant
    // has none, so the slope is zero.
    let (left, _) = ends([1.0, 1.0, 3.0, 4.0]);
    assert_eq!(left.abs().to_bits(), 0, "left, flat");
}

#[test]
fn each_column_keeps_its_own_shape() {
    let x = [0.0, 1.0, 2.5, 3.0, 4.0];
    let rising = [0.0, 1.0, 1.5, 4.0, 4.5];
    let turning = [2.0, -1.0, 3.0, 3.0, 0.0];
    let both: Vec<[f64; 2]> = rising.iter().zip(turning).map(|(a, b)| [*a, b]).collect();
    let pair = PchipInterpolator::new(&x, &both).expect("valid");
    let first = PchipInterpolator::new(&x, &rising).expect("valid");
    let second = PchipInterpolator::new(&x, &turning).expect("valid");
    for q in [-1.0, 0.0, 0.3, 1.0, 2.0, 2.75, 3.5, 4.0, 5.0] {
        let got = pair.eval(&q).expect("extends");
        assert_eq!(Some(got[0].to_bits()), bits(first.eval(&q)), "x = {q}");
        assert_eq!(Some(got[1].to_bits()), bits(second.eval(&q)), "x = {q}");
    }
}

#[test]
fn modes_f32_and_thread_safety() {
    let pchip: PchipInterpolator<f64> =
        PchipInterpolator::new(&[0.0, 1.0, 2.0], &[0.0, 1.0, 4.0]).expect("valid");
    assert_eq!(pchip.extrapolate(), Extrapolate::Extend);
    assert_eq!(pchip.domain().0.to_bits(), 0.0_f64.to_bits());
    assert_eq!(pchip.domain().1.to_bits(), 2.0_f64.to_bits());
    assert!(pchip.eval(&3.0).is_some());
    let forbid = pchip.clone().with_extrapolate(Extrapolate::Forbid);
    assert_eq!(forbid.extrapolate(), Extrapolate::Forbid);
    assert_eq!(forbid.as_ppoly().extrapolate(), Extrapolate::Forbid);
    assert_eq!(forbid.eval(&3.0), None);
    assert_eq!(forbid.eval(&f64::NAN), None);
    let clamp = pchip.clone().with_extrapolate(Extrapolate::Clamp);
    assert_eq!(bits(clamp.eval(&3.0)), bits(pchip.eval(&2.0)));
    let periodic = pchip.with_extrapolate(Extrapolate::Periodic);
    assert_eq!(bits(periodic.eval(&2.5)), bits(periodic.eval(&0.5)));

    let single: PchipInterpolator<f32> =
        PchipInterpolator::new(&[0.0, 1.0, 2.0, 3.0], &[0.0, 0.0, 1.0, 1.0]).expect("valid");
    assert_eq!(single.eval(&1.5).map(f32::to_bits), Some(0.5_f32.to_bits()));
    assert_eq!(single.eval(&0.5).map(f32::to_bits), Some(0.0_f32.to_bits()));

    fn assert_send_sync<T: Send + Sync>() {}
    assert_send_sync::<PchipInterpolator<f64>>();
    assert_send_sync::<PchipInterpolator<f32, [f32; 3]>>();
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

/// Values for those coordinates. Some are drawn from a short list, so that equal
/// neighbors (flat runs) and turns are common.
fn values() -> impl Strategy<Value = Vec<f64>> {
    prop::collection::vec(
        prop_oneof![
            -10.0_f64..10.0,
            prop::sample::select(vec![-1.0, 0.0, 0.0, 2.0])
        ],
        13,
    )
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

proptest! {
    #[test]
    fn the_interpolant_goes_through_the_data_with_a_continuous_slope(
        x in coordinates(),
        values in values(),
    ) {
        let n = x.len();
        let y = &values[..n];
        let pchip = PchipInterpolator::new(&x, y).expect("valid data");
        let c = pchip.as_ppoly().coefficients();
        let scale = y.iter().fold(1.0_f64, |m, v| m.max(v.abs()));
        for i in 0..n - 1 {
            prop_assert_eq!(c[4 * i + 3].to_bits(), y[i].to_bits());
            let h = x[i + 1] - x[i];
            let at_right = piece(c, i, h);
            prop_assert!((at_right[0] - y[i + 1]).abs() <= 1e-12 * scale, "piece {} ends at {} for {}", i, at_right[0], y[i + 1]);
            if i + 2 < n {
                let slope = c[4 * (i + 1) + 2];
                prop_assert!((at_right[1] - slope).abs() <= 1e-10 * scale / h, "slope at point {}: {} and {}", i + 1, at_right[1], slope);
            }
        }
    }

    #[test]
    fn no_piece_leaves_the_range_of_its_two_data_values(
        x in coordinates(),
        values in values(),
        ts in prop::collection::vec(0.0_f64..=1.0, 1..12),
    ) {
        let n = x.len();
        let y = &values[..n];
        let pchip = PchipInterpolator::new(&x, y).expect("valid data");
        let scale = y.iter().fold(1.0_f64, |m, v| m.max(v.abs()));
        for i in 0..n - 1 {
            let (low, high) = (y[i].min(y[i + 1]), y[i].max(y[i + 1]));
            for t in &ts {
                let q = (x[i] + t * (x[i + 1] - x[i])).min(x[i + 1]);
                let got = pchip.eval(&q).expect("inside");
                prop_assert!(low - 1e-12 * scale <= got && got <= high + 1e-12 * scale, "interval {} at {}: {} outside [{}, {}]", i, q, got, low, high);
            }
        }
    }

    #[test]
    fn monotone_data_gives_a_monotone_interpolant(
        x in coordinates(),
        steps in prop::collection::vec(prop_oneof![0.0_f64..5.0, Just(0.0)], 13),
        falling in any::<bool>(),
        ts in prop::collection::vec(0.0_f64..=1.0, 2..40),
    ) {
        let n = x.len();
        let mut y = vec![0.0];
        for step in &steps[..n - 1] {
            y.push(y[y.len() - 1] + if falling { -step } else { *step });
        }
        let pchip = PchipInterpolator::new(&x, &y).expect("valid data");
        let (lo, hi) = pchip.domain();
        let mut qs: Vec<f64> = ts.iter().map(|t| lo + t * (hi - lo)).collect();
        qs.sort_by(f64::total_cmp);
        let scale = y.iter().fold(1.0_f64, |m, v| m.max(v.abs()));
        for pair in qs.windows(2) {
            let (a, b) = (pchip.eval(&pair[0]).expect("inside"), pchip.eval(&pair[1]).expect("inside"));
            let rise = if falling { a - b } else { b - a };
            prop_assert!(rise >= -1e-12 * scale, "from {} to {}: {} then {}", pair[0], pair[1], a, b);
        }
    }

    #[test]
    fn a_slope_is_zero_at_a_turn_and_never_over_three_times_a_secant_beside_it(
        x in coordinates(),
        values in values(),
    ) {
        let n = x.len();
        let y = &values[..n];
        let pchip = PchipInterpolator::new(&x, y).expect("valid data");
        let secant = |i: usize| (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
        for i in 1..n - 1 {
            let (before, after, slope) = (secant(i - 1), secant(i), slope_at(&pchip, i));
            if before * after <= 0.0 {
                // A turn, or flat on one side.
                prop_assert_eq!(slope.abs().to_bits(), 0, "point {}: secants {} and {}", i, before, after);
            } else {
                prop_assert!(slope * before > 0.0, "point {}: the slope has the sign of the secants", i);
                let bound = 3.0 * before.abs().min(after.abs());
                prop_assert!(slope.abs() <= bound * (1.0 + 1e-12), "point {}: {} over {}", i, slope, bound);
            }
        }
        // An end slope never points against its secant, and is three times it at
        // most.
        if n >= 3 {
            for (i, m) in [(0, secant(0)), (n - 1, secant(n - 2))] {
                let slope = slope_at(&pchip, i);
                let scale = m.abs().max(1.0);
                prop_assert!(slope * m >= -1e-9 * scale, "end {}: slope {} against secant {}", i, slope, m);
                prop_assert!(slope.abs() <= 3.0 * m.abs() + 1e-9 * scale, "end {}: {} over three times {}", i, slope, m);
            }
        }
    }

    #[test]
    fn evaluation_is_that_of_the_piecewise_polynomial(
        x in coordinates(),
        values in values(),
        qs in prop::collection::vec(-12.0_f64..40.0, 1..100),
    ) {
        let y = &values[..x.len()];
        for mode in [Extrapolate::Forbid, Extrapolate::Clamp, Extrapolate::Extend, Extrapolate::Periodic] {
            let pchip = PchipInterpolator::new(&x, y).expect("valid data").with_extrapolate(mode);
            let poly = pchip.as_ppoly();
            prop_assert_eq!(poly.extrapolate(), mode);
            let mut hint = Hint::new();
            for q in &qs {
                prop_assert_eq!(bits(pchip.eval(q)), bits(poly.eval(q)), "{:?} x = {}", mode, q);
                prop_assert_eq!(bits(pchip.eval_hinted(*q, &mut hint)), bits(poly.eval(q)));
            }
            let (mut out, mut want) = (vec![0.0; qs.len()], vec![0.0; qs.len()]);
            prop_assert_eq!(pchip.eval_batch_into(&qs, &mut out, -7.0), Ok(()));
            prop_assert_eq!(poly.eval_batch_into(&qs, &mut want, -7.0), Ok(()));
            for (got, want) in out.iter().zip(&want) {
                prop_assert_eq!(got.to_bits(), want.to_bits());
            }
        }
    }
}
