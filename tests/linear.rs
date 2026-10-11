//! `Linear1D`: the golden data, the input checks through the constructor, exact
//! values at the data points, the four extrapolation modes, and the properties
//! that hold for any polygon.

mod common;

use common::{assert_close, load, reals, Real, DIRECT_F64};
use interpx::{Extrapolate, Hint, InterpolateError as E, Interpolator, Linear1D};
use proptest::prelude::*;
use serde_json::Value as Json;

const MODES: [Extrapolate; 4] = [
    Extrapolate::Forbid,
    Extrapolate::Clamp,
    Extrapolate::Extend,
    Extrapolate::Periodic,
];

fn bits(v: Option<f64>) -> Option<u64> {
    v.map(f64::to_bits)
}

fn nums(v: &Json) -> Vec<f64> {
    reals(&serde_json::from_value::<Vec<Real>>(v.clone()).unwrap_or_default())
}

// --- Golden data -------------------------------------------------------------

#[test]
fn numpy_interp_cases_evaluate_to_the_golden_values() {
    let file = load::<Json>("oned").expect("golden file oned");
    let (mut clamped, mut filled, mut infinite) = (0, 0, 0);
    for case in &file.cases {
        if case["reference"] != "numpy.interp" {
            continue;
        }
        let name = case["name"].as_str().expect("name");
        let (x, y) = (nums(&case["x"]), nums(&case["y"]));
        let (queries, expected) = (nums(&case["queries"]), nums(&case["expected"]));
        assert_eq!(queries.len(), expected.len(), "{name}");
        let line = Linear1D::new(&x, &y).expect("valid data");
        assert_eq!(line.extrapolate(), Extrapolate::Clamp, "NumPy's default");

        // `left` and `right` are the values outside: `Forbid` with a fill.
        let fills = case["options"]["left"]
            .as_f64()
            .zip(case["options"]["right"].as_f64());
        let got: Vec<f64> = match fills {
            None => queries
                .iter()
                .map(|q| line.eval(q).expect("clamped"))
                .collect(),
            Some((left, right)) => {
                let forbid = line.clone().with_extrapolate(Extrapolate::Forbid);
                queries
                    .iter()
                    .map(|q| forbid.eval_or(q, if *q < x[0] { left } else { right }))
                    .collect()
            }
        };
        for (k, (g, w)) in got.iter().zip(&expected).enumerate() {
            let what = format!("{name} query {k} (x = {})", queries[k]);
            assert_close(*g, *w, DIRECT_F64, &what);
        }

        // The batch writes the same values.
        if fills.is_none() {
            let mut out = vec![0.0; queries.len()];
            assert_eq!(line.eval_batch_into(&queries, &mut out, f64::NAN), Ok(()));
            for (g, b) in got.iter().zip(&out) {
                assert_eq!(g.to_bits(), b.to_bits(), "{name}: batch");
            }
        }

        match (fills.is_some(), queries.iter().any(|q| q.is_infinite())) {
            (true, _) => filled += 1,
            (false, true) => infinite += 1,
            (false, false) => clamped += 1,
        }
    }
    assert_eq!((clamped, filled, infinite), (9, 1, 2));
}

// --- Input checks through the constructor --------------------------------------

#[test]
fn the_constructor_reports_each_fault() {
    assert!(Linear1D::new(&[0.0, 1.0], &[1.0, 2.0]).is_ok());

    let cases: [(&[f64], &[f64], E); 8] = [
        (
            &[0.0, 1.0, 2.0],
            &[1.0, 2.0],
            E::LengthMismatch {
                expected: 3,
                actual: 2,
            },
        ),
        // NumPy takes one point as a constant.
        (
            &[0.0],
            &[1.0],
            E::InsufficientPoints {
                required: 2,
                actual: 1,
            },
        ),
        (
            &[],
            &[],
            E::InsufficientPoints {
                required: 2,
                actual: 0,
            },
        ),
        (&[0.0, f64::NAN], &[1.0, 2.0], E::NonFinite { index: 1 }),
        // NumPy takes an infinite value.
        (
            &[0.0, 1.0],
            &[f64::INFINITY, 2.0],
            E::NonFinite { index: 0 },
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
        // Finite values whose difference is not: the slope to point 2 overflows.
        (
            &[0.0, 1.0, 2.0],
            &[0.0, -1e308, 1e308],
            E::NonFinite { index: 2 },
        ),
    ];
    for (x, y, want) in cases {
        assert_eq!(Linear1D::new(x, y).err(), Some(want));
    }
    // A slope too steep to hold, between points that are too close.
    assert_eq!(
        Linear1D::new(&[0.0, 1e-300], &[0.0, 1e10]).err(),
        Some(E::NonFinite { index: 1 })
    );
}

#[test]
fn the_constructor_keeps_the_order_of_the_checks() {
    // Different lengths, and too few points.
    assert_eq!(
        Linear1D::new(&[0.0], &[1.0, 2.0]).err(),
        Some(E::LengthMismatch {
            expected: 1,
            actual: 2
        })
    );
    // Too few points, and a NaN.
    assert_eq!(
        Linear1D::new(&[f64::NAN], &[1.0]).err(),
        Some(E::InsufficientPoints {
            required: 2,
            actual: 1
        })
    );
    // A NaN value, and unsorted coordinates.
    assert_eq!(
        Linear1D::new(&[1.0, 0.0], &[f64::NAN, 1.0]).err(),
        Some(E::NonFinite { index: 0 })
    );
    // Unsorted coordinates, and a slope that overflows.
    assert_eq!(
        Linear1D::new(&[1.0, 0.0], &[-1e308, 1e308]).err(),
        Some(E::NotSorted { index: 1 })
    );
}

// --- Values at the data points and outside ---------------------------------------

#[test]
fn each_mode_treats_the_outside_as_documented() {
    let line: Linear1D<f64> = Linear1D::new(&[0.0, 1.0, 3.0], &[0.0, 10.0, 30.5]).expect("valid");
    assert_eq!(line.domain().0.to_bits(), 0.0_f64.to_bits());
    assert_eq!(line.domain().1.to_bits(), 3.0_f64.to_bits());
    let at = |mode, x: f64| bits(line.clone().with_extrapolate(mode).eval(&x));

    assert_eq!(at(Extrapolate::Forbid, -0.5), None);
    assert_eq!(at(Extrapolate::Forbid, 3.5), None);
    // The end values, to infinity.
    assert_eq!(at(Extrapolate::Clamp, -0.5), bits(Some(0.0)));
    assert_eq!(at(Extrapolate::Clamp, 3.5), bits(Some(30.5)));
    assert_eq!(at(Extrapolate::Clamp, f64::NEG_INFINITY), bits(Some(0.0)));
    assert_eq!(at(Extrapolate::Clamp, f64::INFINITY), bits(Some(30.5)));
    // The end segments continued: slope 10 on the left, 10.25 on the right.
    assert_eq!(at(Extrapolate::Extend, -0.5), bits(Some(-5.0)));
    assert_eq!(at(Extrapolate::Extend, 5.0), bits(Some(10.0 + 10.25 * 4.0)));
    // The period is 3: 3.5 is 0.5, -1 is 2, and 3 itself is 0.
    assert_eq!(at(Extrapolate::Periodic, 3.5), bits(Some(5.0)));
    assert_eq!(at(Extrapolate::Periodic, -1.0), bits(Some(20.25)));
    assert_eq!(at(Extrapolate::Periodic, 3.0), bits(Some(0.0)));
    assert_eq!(at(Extrapolate::Periodic, f64::INFINITY), None);

    for mode in MODES {
        assert_eq!(at(mode, 0.5), bits(Some(5.0)), "{mode:?}");
        assert_eq!(at(mode, 2.0), bits(Some(20.25)), "{mode:?}");
        assert_eq!(at(mode, f64::NAN), None, "{mode:?}");
    }
}

#[test]
fn the_last_data_value_is_returned_as_it_is_where_the_formula_is_not_exact() {
    // On [0.1, 0.7] from 0.3 to 0.9 the formula gives 0.9000000000000001 at 0.7.
    let (x, y) = ([0.1_f64, 0.7], [0.3_f64, 0.9]);
    let slope = (y[1] - y[0]) * (x[1] - x[0]).recip();
    let by_formula = y[0] + (x[1] - x[0]) * slope;
    assert_ne!(
        by_formula.to_bits(),
        y[1].to_bits(),
        "the formula is exact here: pick other data"
    );

    let line = Linear1D::new(&x, &y).expect("valid");
    assert_eq!(bits(line.eval(&0.7)), bits(Some(0.9)));
    assert_eq!(bits(line.eval(&5.0)), bits(Some(0.9)), "clamped");
    let mut hint = Hint::new();
    assert_eq!(bits(line.eval_hinted(0.7, &mut hint)), bits(Some(0.9)));
    let mut out = [0.0; 2];
    assert_eq!(
        line.eval_batch_into(&[0.7, 9.0], &mut out, f64::NAN),
        Ok(())
    );
    assert_eq!(out.map(f64::to_bits), [0.9_f64, 0.9].map(f64::to_bits));
    // The piecewise polynomial evaluates the formula.
    assert_eq!(bits(line.as_ppoly().eval(&0.7)), bits(Some(by_formula)));

    // Continued past the end, it is the line again, not the data value.
    let extend = line.with_extrapolate(Extrapolate::Extend);
    assert_eq!(bits(extend.eval(&0.7)), bits(Some(0.9)));
    assert_close(
        extend.eval(&1.3).expect("extends"),
        1.5,
        DIRECT_F64,
        "the line continued",
    );
}

#[test]
fn the_piecewise_polynomial_gives_derivatives_and_integrals() {
    let line = Linear1D::new(&[0.0, 1.0, 3.0], &[0.0, 10.0, 30.0]).expect("valid");
    let poly = line.as_ppoly();
    assert_eq!(poly.degree(), 1);
    assert_eq!(poly.breakpoints().len(), 3);
    assert_eq!(poly.extrapolate(), Extrapolate::Clamp);
    assert_eq!(bits(poly.eval_derivative(0.5, 1)), bits(Some(10.0)));
    assert_eq!(bits(poly.integrate(0.0, 3.0)), bits(Some(45.0)));
    // Clamped: 30 over two more units.
    assert_eq!(bits(poly.integrate(0.0, 5.0)), bits(Some(105.0)));

    // The mode follows the interpolator.
    for mode in MODES {
        let line = line.clone().with_extrapolate(mode);
        assert_eq!(line.extrapolate(), mode);
        assert_eq!(line.as_ppoly().extrapolate(), mode);
    }
}

#[test]
fn array_values_f32_and_thread_safety() {
    let path: Linear1D<f32, [f32; 2]> =
        Linear1D::new(&[0.0, 2.0, 3.0], &[[0.0, 1.0], [4.0, 1.0], [4.0, -1.0]]).expect("valid");
    let same = |got: Option<[f32; 2]>, want: [f32; 2]| {
        assert_eq!(
            got.map(|v| v.map(f32::to_bits)),
            Some(want.map(f32::to_bits))
        );
    };
    same(path.eval(&0.5), [1.0, 1.0]);
    same(path.eval(&2.5), [4.0, 0.0]);
    same(path.eval(&3.0), [4.0, -1.0]);
    same(path.eval(&9.0), [4.0, -1.0]);

    fn assert_send_sync<T: Send + Sync>() {}
    assert_send_sync::<Linear1D<f64>>();
    assert_send_sync::<Linear1D<f32, [f32; 3]>>();
}

#[test]
fn the_batch_rejects_a_length_mismatch_and_leaves_the_buffer() {
    let line: Linear1D<f64> = Linear1D::new(&[0.0, 1.0], &[0.0, 2.0]).expect("valid");
    for (n_points, n_out) in [(2, 3), (3, 2), (0, 1), (1, 0)] {
        let mut out = vec![7.0; n_out];
        assert_eq!(
            line.eval_batch_into(&vec![0.5; n_points], &mut out, 0.0),
            Err(E::LengthMismatch {
                expected: n_points,
                actual: n_out,
            })
        );
        assert!(out.iter().all(|v| v.to_bits() == 7.0_f64.to_bits()));
    }
}

// --- Properties ----------------------------------------------------------------

/// Coordinates that increase strictly, and a value for each.
fn data() -> impl Strategy<Value = (Vec<f64>, Vec<f64>)> {
    (2_usize..12)
        .prop_flat_map(|n| {
            (
                -5.0_f64..5.0,
                prop::collection::vec(0.01_f64..3.0, n - 1),
                prop::collection::vec(-10.0_f64..10.0, n),
            )
        })
        .prop_map(|(start, widths, y)| {
            let mut x = vec![start];
            for w in widths {
                x.push(x[x.len() - 1] + w);
            }
            (x, y)
        })
}

proptest! {
    #[test]
    fn every_data_value_is_returned_as_it_is(
        (x, y) in data(),
        distance in 0.001_f64..50.0,
    ) {
        let last = x.len() - 1;
        for mode in MODES {
            let line = Linear1D::new(&x, &y).expect("valid").with_extrapolate(mode);
            let mut hint = Hint::new();
            for i in 0..=last {
                // A period holds its left end only: the last point is the first.
                let want = if mode == Extrapolate::Periodic && i == last { y[0] } else { y[i] };
                prop_assert_eq!(bits(line.eval(&x[i])), bits(Some(want)), "{:?} point {}", mode, i);
                prop_assert_eq!(bits(line.eval_hinted(x[i], &mut hint)), bits(Some(want)));
            }
        }
        // Clamped, the end values hold outside as well.
        let clamp = Linear1D::new(&x, &y).expect("valid");
        prop_assert_eq!(bits(clamp.eval(&(x[0] - distance))), bits(Some(y[0])));
        prop_assert_eq!(bits(clamp.eval(&(x[last] + distance))), bits(Some(y[last])));
    }

    #[test]
    fn a_value_lies_between_the_two_data_values_beside_it(
        (x, y) in data(),
        t in 0.0_f64..1.0,
    ) {
        let line = Linear1D::new(&x, &y).expect("valid");
        for i in 0..x.len() - 1 {
            let q = x[i] + t * (x[i + 1] - x[i]);
            let got = line.eval(&q).expect("inside");
            let (low, high) = (y[i].min(y[i + 1]), y[i].max(y[i + 1]));
            // One rounding of the product can carry it past a neighbor.
            let slack = 1e-13 * (high - low).max(1.0);
            prop_assert!(low - slack <= got && got <= high + slack, "{} outside [{}, {}]", got, low, high);
            // The formula, written apart from the library.
            let want = y[i] + (q - x[i]) * ((y[i + 1] - y[i]) / (x[i + 1] - x[i]));
            if q < x[i + 1] {
                prop_assert!((got - want).abs() <= 1e-13 * want.abs().max(1.0), "{} for {}", got, want);
            }
        }
    }

    #[test]
    fn a_hint_and_the_batch_never_change_the_result(
        (x, y) in data(),
        qs in prop::collection::vec(-12.0_f64..40.0, 1..150),
    ) {
        let mut ascending = qs.clone();
        ascending.sort_by(f64::total_cmp);
        // Dense and in order, ending on the last point: the batch takes its
        // hinted path there.
        let (lo, hi) = (x[0], x[x.len() - 1]);
        let dense: Vec<f64> = (0..=qs.len()).map(|i| lo + (hi - lo) * i as f64 / qs.len() as f64).chain([hi]).collect();
        for mode in MODES {
            let line = Linear1D::new(&x, &y).expect("valid").with_extrapolate(mode);
            let mut hint = Hint::new();
            for q in qs.iter().chain([&f64::NAN, &hi]) {
                prop_assert_eq!(bits(line.eval_hinted(*q, &mut hint)), bits(line.eval(q)), "{:?} x = {}", mode, q);
            }
            for batch in [&qs, &ascending, &dense] {
                let mut out = vec![0.0; batch.len()];
                prop_assert_eq!(line.eval_batch_into(batch, &mut out, -7.0), Ok(()));
                for (q, got) in batch.iter().zip(&out) {
                    prop_assert_eq!(got.to_bits(), line.eval_or(q, -7.0).to_bits(), "{:?} x = {}", mode, q);
                }
            }
        }
    }

    #[test]
    fn the_piecewise_polynomial_is_the_same_function(
        (x, y) in data(),
        qs in prop::collection::vec(-12.0_f64..40.0, 1..60),
    ) {
        let hi = x[x.len() - 1];
        for mode in MODES {
            let line = Linear1D::new(&x, &y).expect("valid").with_extrapolate(mode);
            let poly = line.as_ppoly();
            for q in &qs {
                let (got, want) = (line.eval(q), poly.eval(q));
                // Left of the last point the two run the same arithmetic. At it,
                // and right of it when clamped, the polynomial uses the formula.
                // A wrapped query can round onto the last point from anywhere.
                let same_arithmetic = mode != Extrapolate::Periodic
                    && (*q < hi || mode != Extrapolate::Clamp);
                if same_arithmetic {
                    prop_assert_eq!(bits(got), bits(want), "{:?} x = {}", mode, q);
                } else {
                    let (got, want) = (got.expect("a value"), want.expect("a value"));
                    prop_assert!((got - want).abs() <= 1e-12 * want.abs().max(1.0), "{:?} x = {}: {} and {}", mode, q, got, want);
                }
            }
        }
    }
}
