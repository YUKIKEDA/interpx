//! `Akima1D`: the golden data, the input checks through the constructor, the
//! rules for the slopes of both methods, and the properties of the interpolant.

mod common;

use common::{assert_close, load, reals, Real, Tol, DIRECT_F64};
use interpx::{Akima1D, AkimaMethod, Extrapolate, Hint, InterpolateError as E, Interpolator};
use proptest::prelude::*;
use serde_json::Value as Json;

const METHODS: [AkimaMethod; 2] = [AkimaMethod::Akima, AkimaMethod::Makima];

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

/// The method of a golden case, when it is one of `Akima1DInterpolator`.
fn method_of(case: &Json) -> Option<AkimaMethod> {
    if case["reference"] != "scipy.interpolate.Akima1DInterpolator" {
        return None;
    }
    match case["options"]["method"].as_str()? {
        "akima" => Some(AkimaMethod::Akima),
        "makima" => Some(AkimaMethod::Makima),
        _ => None,
    }
}

/// The slope at data point `i`, read from the piecewise polynomial: the linear
/// coefficient of the piece that starts there, or for the last point the slope
/// at the end of the last piece.
fn slope_at(akima: &Akima1D<f64>, i: usize) -> f64 {
    let poly = akima.as_ppoly();
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
fn scipy_akima_interpolants_are_rebuilt_from_their_data() {
    let file = load::<Json>("oned").expect("golden file oned");
    let (mut akima, mut makima, mut extended, mut small) = (0, 0, 0, 0);
    for case in &file.cases {
        let (Some(method), Some(_)) = (method_of(case), case.get("coefficients")) else {
            continue;
        };
        let name = case["name"].as_str().expect("name");
        let (x, y) = (nums(&case["x"]), nums(&case["y"]));
        let built = Akima1D::new(&x, &y, method).expect("valid data");
        assert_eq!(
            built.extrapolate(),
            Extrapolate::Forbid,
            "{name}: SciPy's default"
        );
        assert_eq!(built.method(), method, "{name}");
        // `extrapolate=None` is SciPy's default for Akima: no value outside.
        let interpolant = if case["options"]["extrapolate"] == Json::Bool(true) {
            extended += 1;
            built.with_extrapolate(Extrapolate::Extend)
        } else {
            built
        };

        // A coefficient is compared by what it adds to the value at the far end
        // of its interval (`tests/cubic_spline.rs` says why): SciPy stores
        // `c[power][interval]`.
        let scipy = rows(&case["coefficients"]);
        let poly = interpolant.as_ppoly();
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
            ("expected", &|q| interpolant.eval(&q)),
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
        match method {
            AkimaMethod::Makima => makima += 1,
            _ => akima += 1,
        }
        small += usize::from(x.len() <= 4);
    }
    // Per method: 5 on the standard data, one of them again with extrapolation,
    // and 6 chosen for the rules, 3 of them with four points or fewer.
    assert_eq!((akima, makima, extended, small), (12, 12, 2, 6));
}

#[test]
fn scipy_akima_of_two_columns_is_rebuilt_column_by_column() {
    let file = load::<Json>("oned").expect("golden file oned");
    let mut cases = 0;
    for case in &file.cases {
        let Some(method) = method_of(case) else {
            continue;
        };
        if case.get("derivative_2").is_none() || case.get("coefficients").is_some() {
            continue;
        }
        let name = case["name"].as_str().expect("name");
        let y: Vec<[f64; 2]> = rows(&case["y"])
            .iter()
            .map(|r| r.as_slice().try_into().expect("two columns"))
            .collect();
        let interpolant = Akima1D::new(&nums(&case["x"]), &y, method).expect("valid data");
        let poly = interpolant.as_ppoly();
        let queries = nums(&case["queries"]);
        for (field, nu) in [("expected", 0), ("derivative_1", 1), ("derivative_2", 2)] {
            for (k, (q, want)) in queries.iter().zip(rows(&case[field])).enumerate() {
                let what = format!("{name} {field} query {k}");
                match poly.eval_derivative(*q, nu) {
                    None => assert!(want.iter().all(|w| w.is_nan()), "{what}: None"),
                    Some(got) => {
                        for column in 0..2 {
                            assert_close(got[column], want[column], DIRECT_F64, &what);
                        }
                    }
                }
            }
        }
        cases += 1;
    }
    assert_eq!(cases, 2);
}

#[test]
fn a_column_is_interpolated_on_its_own_where_scipy_lets_another_one_in() {
    // The one difference from SciPy (docs/design.ja.md section 6.1). The first
    // column has values near 1e-10, the second near 1. SciPy takes the cutoff
    // for vanishing weights from both together, so beside the second column the
    // first one gets other slopes than it has alone.
    let file = load::<Json>("oned").expect("golden file oned");
    // The values are near 1e-10: an absolute tolerance would hide everything.
    let relative = Tol {
        rtol: 1e-12,
        atol: 1e-24,
    };
    let mut cases = 0;
    for case in &file.cases {
        let (Some(method), Some(alone)) = (method_of(case), case.get("expected_alone")) else {
            continue;
        };
        let name = case["name"].as_str().expect("name");
        let x = nums(&case["x"]);
        let y: Vec<[f64; 2]> = rows(&case["y"])
            .iter()
            .map(|r| r.as_slice().try_into().expect("two columns"))
            .collect();
        let first: Vec<f64> = y.iter().map(|v| v[0]).collect();
        let paired = Akima1D::new(&x, &y, method).expect("valid data");
        let single = Akima1D::new(&x, &first, method).expect("valid data");

        let (alone, with_the_other) = (nums(alone), rows(&case["expected"]));
        let mut scipy_differs: f64 = 0.0;
        for (k, q) in nums(&case["queries"]).iter().enumerate() {
            let got = paired.eval(q).expect("inside the data");
            // Beside the other column, the same bits as alone.
            assert_eq!(
                Some(got[0].to_bits()),
                bits(single.eval(q)),
                "{name} query {k}"
            );
            // And what SciPy gives the column alone.
            assert_close(got[0], alone[k], relative, &format!("{name} query {k}"));
            scipy_differs = scipy_differs.max((with_the_other[k][0] - alone[k]).abs());
        }
        // SciPy's own two results for this column are apart: the case shows the
        // difference, it does not only claim it.
        assert!(
            scipy_differs > 1e-11,
            "{name}: SciPy differs by {scipy_differs:e} only"
        );
        cases += 1;
    }
    assert_eq!(cases, 2, "one per method");
}

#[test]
fn the_cutoff_for_vanishing_weights_is_scipys() {
    // Small secants beside large ones, on each side of the cutoff: far from it
    // (1e-12 and 1e-8 of the largest sum of weights) and close to it (0.75e-9 and
    // 1.5e-9; the generator says how). The small slopes are compared relatively:
    // an absolute tolerance would pass anything.
    let file = load::<Json>("oned").expect("golden file oned");
    let relative = Tol {
        rtol: 1e-9,
        atol: 0.0,
    };
    let (mut under, mut over) = (0, 0);
    for case in &file.cases {
        let (Some(method), Some(slopes)) = (method_of(case), case.get("slopes")) else {
            continue;
        };
        let name = case["name"].as_str().expect("name");
        let (x, y) = (nums(&case["x"]), nums(&case["y"]));
        let interpolant = Akima1D::new(&x, &y, method).expect("valid data");
        let slopes = nums(slopes);
        for (i, want) in slopes.iter().enumerate() {
            let what = format!("{name} slope at point {i}");
            assert_close(slope_at(&interpolant, i), *want, relative, &what);
        }
        // At point 2 the secants are `scale` times 1, 2 before it and 5, 3 after
        // it. Under the cutoff the slope is the mean of the outer two, 2 `scale`.
        // Over it, it is the weighted mean, which is not.
        let scale = y[1];
        let outer_mean = 2.0 * scale;
        let is_outer_mean = (slopes[2] - outer_mean).abs() <= 1e-6 * scale;
        if name.contains("under_the_cutoff") {
            assert!(is_outer_mean, "{name}: SciPy has {}", slopes[2]);
            under += 1;
        } else {
            assert!(!is_outer_mean, "{name}: SciPy has {}", slopes[2]);
            over += 1;
        }
    }
    assert_eq!(
        (under, over),
        (4, 4),
        "both methods, far from the cutoff and close to it"
    );
}

// --- Input checks through the constructor --------------------------------------

#[test]
fn the_constructor_reports_each_fault() {
    let cases: [(&[f64], &[f64], E); 7] = [
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
    ];
    for method in METHODS {
        assert!(Akima1D::new(&[0.0, 1.0], &[1.0, 2.0], method).is_ok());
        for (x, y, want) in cases {
            assert_eq!(Akima1D::new(x, y, method).err(), Some(want), "{method:?}");
        }
        // The order of the checks: lengths, the count, NaN, then the order of `x`.
        assert_eq!(
            Akima1D::new(&[0.0], &[1.0, 2.0], method).err(),
            Some(E::LengthMismatch {
                expected: 1,
                actual: 2
            })
        );
        assert_eq!(
            Akima1D::new(&[1.0, 0.0], &[f64::NAN, 1.0], method).err(),
            Some(E::NonFinite { index: 0 })
        );
        assert_eq!(
            Akima1D::new(&[1.0, 0.0], &[-1e308, 1e308], method).err(),
            Some(E::NotSorted { index: 1 })
        );
        // Secants near the largest number: the ones made up past the ends, and so
        // the slopes, are not finite. An error, not a panic.
        let err = Akima1D::new(&[0.0, 1.0, 2.0, 3.0], &[0.0, 1.5e308, 0.0, 1.5e308], method).err();
        assert!(
            matches!(err, Some(E::NonFinite { index }) if index < 4),
            "{method:?}: {err:?}"
        );
    }
}

// --- The rules for the slopes ----------------------------------------------------

#[test]
fn two_points_give_the_line_through_them() {
    for method in METHODS {
        let line = Akima1D::new(&[1.0, 3.0], &[2.0, 8.0], method).expect("valid");
        let c = line.as_ppoly().coefficients();
        assert_eq!(
            [c[0], c[1], c[2], c[3]].map(f64::to_bits),
            [0.0_f64, 0.0, 3.0, 2.0].map(f64::to_bits),
            "{method:?}"
        );
        assert_eq!(bits(line.eval(&2.0)), bits(Some(5.0)));
    }
}

#[test]
fn a_slope_is_the_mean_of_the_two_secants_beside_it_by_the_weights_of_the_method() {
    // Secants 1, 2, 3, 4, continued to 0 and -1 on the left, 5 and 6 on the right.
    let (x, y) = ([0.0, 1.0, 2.0, 3.0, 4.0], [0.0, 1.0, 3.0, 6.0, 10.0]);
    let akima = Akima1D::new(&x, &y, AkimaMethod::Akima).expect("valid");
    // Akima: every weight is 1, so every slope is the plain mean.
    for (i, want) in [0.5, 1.5, 2.5, 3.5, 4.5].into_iter().enumerate() {
        assert_close(
            slope_at(&akima, i),
            want,
            DIRECT_F64,
            &format!("akima, point {i}"),
        );
    }
    // Makima at point 2, between the secants 2 and 3: the weight of the side
    // after it is |4 - 3| + |4 + 3| / 2 = 4.5, of the side before it
    // |2 - 1| + |2 + 1| / 2 = 2.5, and the slope is 2 + 2.5 / 7.
    let makima = Akima1D::new(&x, &y, AkimaMethod::Makima).expect("valid");
    assert_close(
        slope_at(&makima, 2),
        2.0 + 2.5 / 7.0,
        DIRECT_F64,
        "makima, point 2",
    );
    // At point 0, between the made-up secant 0 and the secant 1: weights
    // |2 - 1| + 1.5 = 2.5 and |0 + 1| + 0.5 = 1.5, so 0 + 1.5 / 4.
    assert_close(slope_at(&makima, 0), 0.375, DIRECT_F64, "makima, point 0");
    // At the last point, between 4 and the made-up 5: weights |6 - 5| + 5.5 and
    // |4 - 3| + 3.5, so 4 + 4.5 / 11.
    assert_close(
        slope_at(&makima, 4),
        4.0 + 4.5 / 11.0,
        DIRECT_F64,
        "makima, point 4",
    );
}

#[test]
fn where_the_weights_vanish_the_slope_is_the_mean_of_the_outer_secants() {
    // A flat run, then a straight rise: at point 2 the secants before it are
    // 0 and 0, and those after it 1 and 1.
    let (x, y) = (
        [0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
        [1.0, 1.0, 1.0, 2.0, 3.0, 4.0],
    );
    let akima = Akima1D::new(&x, &y, AkimaMethod::Akima).expect("valid");
    // Akima's weights are both zero there: (0 + 1) / 2, and the curve dips under
    // the flat run.
    assert_eq!(slope_at(&akima, 2).to_bits(), 0.5_f64.to_bits());
    assert!(akima.eval(&1.5).expect("inside") < 1.0);
    // Makima's are 0 and 1: the slope is that of the flat side, and the flat run
    // is flat, exactly.
    let makima = Akima1D::new(&x, &y, AkimaMethod::Makima).expect("valid");
    assert_eq!(slope_at(&makima, 2).abs().to_bits(), 0);
    for q in [0.0, 0.25, 1.0, 1.5, 1.999, 2.0] {
        assert_eq!(bits(makima.eval(&q)), bits(Some(1.0)), "x = {q}");
    }

    // Data on a line: every weight of Akima vanishes, the cutoff is zero, and the
    // mean of two equal secants is the secant.
    let line = |t: f64| 1.0 + 2.0 * t;
    let x = [0.0, 1.0, 2.5, 3.0, 4.5, 6.0];
    for method in METHODS {
        let on_a_line = Akima1D::new(&x, &x.map(line), method).expect("valid");
        for i in 0..x.len() {
            assert_close(
                slope_at(&on_a_line, i),
                2.0,
                DIRECT_F64,
                &format!("{method:?} point {i}"),
            );
        }
        for q in [0.5, 2.0, 2.75, 5.0] {
            assert_close(
                on_a_line.eval(&q).expect("inside"),
                line(q),
                DIRECT_F64,
                "the line",
            );
        }
    }
}

#[test]
fn a_weight_far_under_the_largest_one_counts_as_zero() {
    // Secants 0, 0, 1e-12, 1e-12, 5, -5: at point 2 the weights add up to 1.5e-12
    // for makima, which is under 1e-9 of the largest sum (15 at point 4), and to
    // zero for Akima. Both take the mean of the outer secants, 0.5e-12, where
    // makima's weights would give 0.
    let x = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
    let y = [0.0, 0.0, 0.0, 1e-12, 2e-12, 5.0, 0.0];
    for method in METHODS {
        let interpolant = Akima1D::new(&x, &y, method).expect("valid");
        assert_eq!(
            slope_at(&interpolant, 2).to_bits(),
            0.5e-12_f64.to_bits(),
            "{method:?}"
        );
    }
    // Without the large secants the same weights are over the cutoff, and makima
    // keeps the flat run flat.
    let makima = Akima1D::new(&x[..5], &y[..5], AkimaMethod::Makima).expect("valid");
    assert_eq!(slope_at(&makima, 2).abs().to_bits(), 0);
}

#[test]
fn modes_f32_and_thread_safety() {
    let akima: Akima1D<f64> =
        Akima1D::new(&[0.0, 1.0, 2.0], &[0.0, 1.0, 4.0], AkimaMethod::Akima).expect("valid");
    assert_eq!(akima.extrapolate(), Extrapolate::Forbid);
    assert_eq!(akima.domain().0.to_bits(), 0.0_f64.to_bits());
    assert_eq!(akima.domain().1.to_bits(), 2.0_f64.to_bits());
    assert_eq!(akima.eval(&3.0), None);
    assert_eq!(akima.eval(&f64::NAN), None);
    let extend = akima.clone().with_extrapolate(Extrapolate::Extend);
    assert_eq!(extend.extrapolate(), Extrapolate::Extend);
    assert_eq!(extend.as_ppoly().extrapolate(), Extrapolate::Extend);
    assert!(extend.eval(&3.0).is_some());
    let clamp = akima.clone().with_extrapolate(Extrapolate::Clamp);
    assert_eq!(bits(clamp.eval(&3.0)), bits(akima.eval(&2.0)));
    let periodic = akima.with_extrapolate(Extrapolate::Periodic);
    assert_eq!(bits(periodic.eval(&2.5)), bits(periodic.eval(&0.5)));

    // Single precision, on a line.
    let single: Akima1D<f32> = Akima1D::new(
        &[0.0, 1.0, 3.0, 4.0],
        &[1.0, 3.0, 7.0, 9.0],
        AkimaMethod::Makima,
    )
    .expect("valid");
    assert_eq!(single.eval(&2.0).map(f32::to_bits), Some(5.0_f32.to_bits()));

    fn assert_send_sync<T: Send + Sync>() {}
    assert_send_sync::<Akima1D<f64>>();
    assert_send_sync::<Akima1D<f32, [f32; 3]>>();
    assert_send_sync::<AkimaMethod>();
}

// --- Properties ----------------------------------------------------------------

/// Coordinates that increase strictly, at least `least` of them.
fn coordinates(least: usize) -> impl Strategy<Value = Vec<f64>> {
    (
        -5.0_f64..5.0,
        prop::collection::vec(0.05_f64..3.0, least - 1..12),
    )
        .prop_map(|(start, widths)| {
            let mut x = vec![start];
            for w in widths {
                x.push(x[x.len() - 1] + w);
            }
            x
        })
}

fn methods() -> impl Strategy<Value = AkimaMethod> {
    prop::sample::select(METHODS.to_vec())
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
        x in coordinates(2),
        values in prop::collection::vec(-10.0_f64..10.0, 13),
        method in methods(),
    ) {
        let n = x.len();
        let y = &values[..n];
        let interpolant = Akima1D::new(&x, y, method).expect("valid data");
        let c = interpolant.as_ppoly().coefficients();
        for i in 0..n - 1 {
            prop_assert_eq!(c[4 * i + 3].to_bits(), y[i].to_bits());
            let h = x[i + 1] - x[i];
            let at_right = piece(c, i, h);
            // Widths are 0.05 at least and values 10 at most, so a slope is up to
            // a few hundred.
            prop_assert!((at_right[0] - y[i + 1]).abs() <= 1e-10, "piece {} ends at {} for {}", i, at_right[0], y[i + 1]);
            if i + 2 < n {
                let slope = c[4 * (i + 1) + 2];
                prop_assert!((at_right[1] - slope).abs() <= 1e-8, "slope at point {}: {} and {}", i + 1, at_right[1], slope);
            }
        }
    }

    #[test]
    fn data_on_a_line_is_reproduced(
        x in coordinates(2),
        slope in -5.0_f64..5.0,
        offset in -5.0_f64..5.0,
        method in methods(),
        ts in prop::collection::vec(0.0_f64..=1.0, 1..20),
    ) {
        let line = |t: f64| offset + slope * t;
        let y: Vec<f64> = x.iter().map(|t| line(*t)).collect();
        let interpolant = Akima1D::new(&x, &y, method).expect("valid data");
        let (lo, hi) = interpolant.domain();
        for t in &ts {
            let q = lo + t * (hi - lo);
            let got = interpolant.eval(&q).expect("inside");
            // The secants of rounded data differ in their last bits, and a
            // weight is a difference of two of them.
            prop_assert!((got - line(q)).abs() <= 1e-9, "{:?} at {}: {} for {}", method, q, got, line(q));
        }
    }

    #[test]
    fn makima_keeps_a_flat_run_of_three_points_flat(
        x in coordinates(3),
        values in prop::collection::vec(-10.0_f64..10.0, 13),
        start in 0_usize..10,
        ts in prop::collection::vec(0.0_f64..=1.0, 1..10),
    ) {
        let n = x.len();
        let mut y = values[..n].to_vec();
        // Three equal values in a row, somewhere.
        let start = start % (n - 2);
        y[start + 1] = y[start];
        y[start + 2] = y[start];
        let makima = Akima1D::new(&x, &y, AkimaMethod::Makima).expect("valid data");
        for i in [start, start + 1] {
            for t in &ts {
                let q = (x[i] + t * (x[i + 1] - x[i])).min(x[i + 1]);
                prop_assert_eq!(bits(makima.eval(&q)), bits(Some(y[start])), "interval {} at {}", i, q);
            }
        }
    }

    #[test]
    fn each_column_is_what_it_is_alone(
        x in coordinates(2),
        first in prop::collection::vec(-10.0_f64..10.0, 13),
        // The second column is much larger or much smaller than the first.
        second in prop::collection::vec(-1.0_f64..1.0, 13),
        exponent in -12_i32..12,
        method in methods(),
        qs in prop::collection::vec(-12.0_f64..40.0, 1..40),
    ) {
        let n = x.len();
        let scale = 10.0_f64.powi(exponent);
        let pairs: Vec<[f64; 2]> = (0..n).map(|i| [first[i], second[i] * scale]).collect();
        let paired = Akima1D::new(&x, &pairs, method).expect("valid data").with_extrapolate(Extrapolate::Extend);
        for column in 0..2 {
            let alone: Vec<f64> = pairs.iter().map(|v| v[column]).collect();
            let single = Akima1D::new(&x, &alone, method).expect("valid data").with_extrapolate(Extrapolate::Extend);
            for q in &qs {
                let got = paired.eval(q).expect("extends");
                prop_assert_eq!(Some(got[column].to_bits()), bits(single.eval(q)), "column {} at {}", column, q);
            }
        }
    }

    #[test]
    fn evaluation_is_that_of_the_piecewise_polynomial(
        x in coordinates(2),
        values in prop::collection::vec(-10.0_f64..10.0, 13),
        method in methods(),
        qs in prop::collection::vec(-12.0_f64..40.0, 1..100),
    ) {
        let y = &values[..x.len()];
        for mode in [Extrapolate::Forbid, Extrapolate::Clamp, Extrapolate::Extend, Extrapolate::Periodic] {
            let interpolant = Akima1D::new(&x, y, method).expect("valid data").with_extrapolate(mode);
            let poly = interpolant.as_ppoly();
            prop_assert_eq!(poly.extrapolate(), mode);
            let mut hint = Hint::new();
            for q in &qs {
                prop_assert_eq!(bits(interpolant.eval(q)), bits(poly.eval(q)), "{:?} x = {}", mode, q);
                prop_assert_eq!(bits(interpolant.eval_hinted(*q, &mut hint)), bits(poly.eval(q)));
            }
            let (mut out, mut want) = (vec![0.0; qs.len()], vec![0.0; qs.len()]);
            prop_assert_eq!(interpolant.eval_batch_into(&qs, &mut out, -7.0), Ok(()));
            prop_assert_eq!(poly.eval_batch_into(&qs, &mut want, -7.0), Ok(()));
            for (got, want) in out.iter().zip(&want) {
                prop_assert_eq!(got.to_bits(), want.to_bits());
            }
        }
    }
}
