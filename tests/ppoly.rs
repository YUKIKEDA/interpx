//! `PPoly`: the golden data, the input checks through the constructor, the
//! meaning of the intervals and of the four extrapolation modes, and the
//! properties that hold for any piecewise polynomial.

mod common;

use common::{assert_close, reals, Real, DIRECT_F64};
use interpx::{Extrapolate, Hint, InterpolateError as E, Interpolator, PPoly};
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

// --- Golden data -------------------------------------------------------------

fn nums(v: &Json) -> Vec<f64> {
    reals(&serde_json::from_value::<Vec<Real>>(v.clone()).unwrap_or_default())
}

/// SciPy stores `c[power][interval]`; `PPoly` takes `[interval][power]`.
fn from_scipy(breakpoints: &Json, coefficients: &Json) -> Result<PPoly<f64>, E> {
    let rows: Vec<Vec<f64>> = serde_json::from_value::<Vec<Vec<Real>>>(coefficients.clone())
        .unwrap_or_default()
        .iter()
        .map(|r| reals(r))
        .collect();
    let intervals = rows.first().map_or(0, Vec::len);
    let flat: Vec<f64> = (0..intervals)
        .flat_map(|i| rows.iter().map(move |row| row[i]))
        .collect();
    PPoly::from_owned(nums(breakpoints), flat, rows.len().saturating_sub(1))
}

/// The mode a golden case was generated with (`extrapolate` in its options).
fn mode_of(case: &Json) -> Extrapolate {
    match &case["options"]["extrapolate"] {
        Json::Bool(true) => Extrapolate::Extend,
        Json::String(s) if s == "periodic" => Extrapolate::Periodic,
        // `false`, and `None` for Akima, whose default is no extrapolation.
        _ => Extrapolate::Forbid,
    }
}

#[test]
fn scipy_piecewise_polynomials_evaluate_to_the_golden_values() {
    let all_cases = common::load_cases(&common::PIECEWISE).expect("golden files");
    let (mut cases, mut modes) = (0, std::collections::BTreeSet::new());
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
        let mode = mode_of(case);
        let p = from_scipy(breakpoints, coefficients)
            .expect("SciPy's coefficients are valid")
            .with_extrapolate(mode);
        let (queries, expected) = (nums(queries), nums(expected));
        assert_eq!(queries.len(), expected.len(), "{name}");

        for (k, (x, want)) in queries.iter().zip(&expected).enumerate() {
            let what = format!("{name} query {k} (x = {x})");
            match p.eval(x) {
                // SciPy returns NaN where it does not extrapolate.
                None => assert!(want.is_nan(), "{what}: None, SciPy has {want}"),
                Some(got) => assert_close(got, *want, DIRECT_F64, &what),
            }
        }

        // The batch form writes the same values, with NaN where there is none.
        let mut out = vec![0.0; queries.len()];
        assert_eq!(p.eval_batch_into(&queries, &mut out, f64::NAN), Ok(()));
        for (k, (got, want)) in out.iter().zip(&expected).enumerate() {
            assert_close(*got, *want, DIRECT_F64, &format!("{name} batch {k}"));
        }

        cases += 1;
        modes.insert(format!("{mode:?}"));
    }
    // CubicSpline 49, PCHIP 19, Akima 24, Hermite 4.
    assert_eq!(cases, 96, "cases with a piecewise polynomial");
    assert_eq!(modes.len(), 3, "extend, forbid, and periodic: {modes:?}");
}

// --- Input checks through the constructor --------------------------------------

#[test]
fn the_constructor_reports_each_fault() {
    let ok = PPoly::new(&[0.0, 1.0, 2.0], &[1.0, 0.0, 2.0, 1.0], 1);
    assert!(ok.is_ok());

    let cases: [(&[f64], &[f64], usize, E); 7] = [
        (
            &[0.0, 1.0, 2.0],
            &[1.0, 0.0, 2.0],
            1,
            E::LengthMismatch {
                expected: 4,
                actual: 3,
            },
        ),
        (
            &[0.0],
            &[],
            1,
            E::InsufficientPoints {
                required: 2,
                actual: 1,
            },
        ),
        (
            &[],
            &[],
            0,
            E::InsufficientPoints {
                required: 2,
                actual: 0,
            },
        ),
        (&[0.0, f64::NAN], &[1.0], 0, E::NonFinite { index: 1 }),
        (&[0.0, 1.0], &[f64::INFINITY], 0, E::NonFinite { index: 0 }),
        (&[0.0, 2.0, 1.0], &[1.0, 2.0], 0, E::NotSorted { index: 2 }),
        (
            &[0.0, 0.0, 1.0],
            &[1.0, 2.0],
            0,
            E::DuplicatePoints {
                first: 0,
                second: 1,
            },
        ),
    ];
    for (breakpoints, coefficients, degree, want) in cases {
        assert_eq!(
            PPoly::new(breakpoints, coefficients, degree).err(),
            Some(want)
        );
        assert_eq!(
            PPoly::from_owned(breakpoints.to_vec(), coefficients.to_vec(), degree).err(),
            Some(want)
        );
    }
}

#[test]
fn the_constructor_keeps_the_order_of_the_checks() {
    // The wrong number of coefficients, and too few break points.
    assert_eq!(
        PPoly::new(&[0.0], &[1.0], 0).err(),
        Some(E::LengthMismatch {
            expected: 0,
            actual: 1
        })
    );
    // A NaN break point and unsorted break points.
    assert_eq!(
        PPoly::new(&[2.0, f64::NAN, 1.0], &[1.0, 2.0], 0).err(),
        Some(E::NonFinite { index: 1 })
    );
    // A degree that cannot be counted asks for more coefficients than any slice has.
    assert_eq!(
        PPoly::new(&[0.0, 1.0], &[1.0], usize::MAX).err(),
        Some(E::LengthMismatch {
            expected: usize::MAX,
            actual: 1
        })
    );
}

// --- Intervals and extrapolation -----------------------------------------------

/// A step function: 1 on [0, 1), 2 on [1, 2], so the two sides of a break point
/// can be told apart.
fn steps() -> Result<PPoly<f64>, E> {
    PPoly::new(&[0.0, 1.0, 2.0], &[1.0, 2.0], 0)
}

#[test]
fn an_interval_owns_its_left_end_and_the_last_one_both_ends() {
    let p = steps().expect("valid");
    assert_eq!(bits(p.eval(&0.0)), bits(Some(1.0)));
    assert_eq!(bits(p.eval(&0.999_999)), bits(Some(1.0)));
    assert_eq!(bits(p.eval(&1.0)), bits(Some(2.0)));
    assert_eq!(bits(p.eval(&2.0)), bits(Some(2.0)));
}

#[test]
fn each_mode_treats_the_outside_as_documented() {
    // x on [0, 1], then 1 + 2(x - 1) on [1, 2].
    let line = PPoly::new(&[0.0, 1.0, 2.0], &[1.0, 0.0, 2.0, 1.0], 1).expect("valid");
    assert_eq!(line.extrapolate(), Extrapolate::Extend, "SciPy's default");
    let at = |mode, x: f64| bits(line.clone().with_extrapolate(mode).eval(&x));

    assert_eq!(at(Extrapolate::Forbid, -0.5), None);
    assert_eq!(at(Extrapolate::Forbid, 2.5), None);

    assert_eq!(at(Extrapolate::Clamp, -0.5), bits(Some(0.0)));
    assert_eq!(at(Extrapolate::Clamp, 2.5), bits(Some(3.0)));

    // The end pieces continue: x on the left, 1 + 2(x - 1) on the right.
    assert_eq!(at(Extrapolate::Extend, -0.5), bits(Some(-0.5)));
    assert_eq!(at(Extrapolate::Extend, 2.5), bits(Some(4.0)));

    // The period is 2: 2.5 is 0.5, -0.5 is 1.5, 4.25 is 0.25.
    assert_eq!(at(Extrapolate::Periodic, 2.5), bits(Some(0.5)));
    assert_eq!(at(Extrapolate::Periodic, -0.5), bits(Some(2.0)));
    assert_eq!(at(Extrapolate::Periodic, 4.25), bits(Some(0.25)));

    // Inside, the mode does not matter.
    for mode in MODES {
        assert_eq!(at(mode, 1.5), bits(Some(2.0)), "{mode:?}");
        assert_eq!(at(mode, 0.0), bits(Some(0.0)), "{mode:?}");
        // The right end belongs to the last interval, except that a period holds
        // its left end only: there it is the left end again (as in SciPy).
        let at_right_end = if mode == Extrapolate::Periodic {
            0.0
        } else {
            3.0
        };
        assert_eq!(at(mode, 2.0), bits(Some(at_right_end)), "{mode:?}");
    }
}

#[test]
fn nan_and_infinite_queries_never_panic() {
    for mode in MODES {
        let p = steps().expect("valid").with_extrapolate(mode);
        let mut hint = Hint::new();
        assert!(p.eval(&f64::NAN).is_none(), "{mode:?}");
        assert!(p.eval_hinted(f64::NAN, &mut hint).is_none(), "{mode:?}");
        for x in [f64::INFINITY, f64::NEG_INFINITY] {
            let got = p.eval(&x);
            assert_eq!(bits(got), bits(p.eval_hinted(x, &mut hint)), "{mode:?}");
            match mode {
                // An infinite point is outside, and it is at no place in a period.
                Extrapolate::Forbid | Extrapolate::Periodic => assert!(got.is_none(), "{mode:?}"),
                // A step function is constant, so its ends reach to infinity.
                _ => assert!(got.is_some(), "{mode:?}"),
            }
        }
    }
}

#[test]
fn the_batch_rejects_a_length_mismatch_and_leaves_the_buffer() {
    // PPoly replaces the provided `eval_batch_into`, so it repeats the check.
    for mode in MODES {
        let p = steps().expect("valid").with_extrapolate(mode);
        for (n_points, n_out) in [(2, 3), (3, 2), (0, 1), (1, 0)] {
            let points = vec![0.5; n_points];
            let mut out = vec![7.0; n_out];
            assert_eq!(
                p.eval_batch_into(&points, &mut out, 0.0),
                Err(E::LengthMismatch {
                    expected: n_points,
                    actual: n_out,
                }),
                "{mode:?}"
            );
            assert!(
                out.iter().all(|v| v.to_bits() == 7.0_f64.to_bits()),
                "{mode:?}: buffer was written"
            );
        }
        let mut empty: [f64; 0] = [];
        assert_eq!(p.eval_batch_into(&[], &mut empty, 0.0), Ok(()));
    }
}

#[test]
fn accessors_return_what_was_given() {
    let p: PPoly<f32, [f32; 2]> =
        PPoly::new(&[0.0, 1.0, 3.0], &[[1.0, 2.0], [3.0, 4.0]], 0).expect("valid");
    assert_eq!(p.degree(), 0);
    assert_eq!(p.breakpoints().len(), 3);
    assert_eq!(p.coefficients().len(), 2);
    assert_eq!(
        p.eval(&2.0).map(|v| v.map(f32::to_bits)),
        Some([3.0_f32, 4.0].map(f32::to_bits))
    );
    assert_eq!(Hint::new(), Hint::default());

    fn assert_send_sync<T: Send + Sync>() {}
    assert_send_sync::<PPoly<f64>>();
    assert_send_sync::<PPoly<f32, [f32; 3]>>();
    assert_send_sync::<Hint>();
}

#[test]
fn a_stale_hint_far_from_the_query_still_finds_the_interval() {
    // Twenty intervals of width 1, piece i is the constant i.
    let breakpoints: Vec<f64> = (0..=20).map(f64::from).collect();
    let coefficients: Vec<f64> = (0..20).map(f64::from).collect();
    let p = PPoly::from_owned(breakpoints, coefficients, 0).expect("valid");
    let mut hint = Hint::new();
    // Forward by one, a jump to the far end, back by one, a jump to the start.
    for x in [0.5, 1.5, 2.5, 19.5, 18.5, 20.0, 0.0, 10.0] {
        assert_eq!(
            bits(p.eval_hinted(x, &mut hint)),
            bits(p.eval(&x)),
            "x = {x}"
        );
    }
}

// --- Properties ----------------------------------------------------------------

/// Break points, coefficients, and a degree of a random piecewise polynomial.
fn pieces() -> impl Strategy<Value = (Vec<f64>, Vec<f64>, usize)> {
    (1_usize..8, 0_usize..5)
        .prop_flat_map(|(intervals, degree)| {
            (
                -5.0_f64..5.0,
                prop::collection::vec(0.01_f64..3.0, intervals),
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

fn queries() -> impl Strategy<Value = Vec<f64>> {
    prop::collection::vec(-30.0_f64..30.0, 1..40)
}

/// Horner's rule on one interval, written apart from the library.
fn horner(coefficients: &[f64], dx: f64) -> f64 {
    coefficients.iter().fold(0.0, |acc, c| acc * dx + c)
}

proptest! {
    #[test]
    fn a_hint_never_changes_the_result(
        (breakpoints, coefficients, degree) in pieces(),
        xs in queries(),
        start in 0_usize..12,
    ) {
        for mode in MODES {
            let p = PPoly::new(&breakpoints, &coefficients, degree)
                .expect("valid")
                .with_extrapolate(mode);
            // A hint left by an earlier, unrelated run of queries.
            let mut hint = Hint::new();
            for _ in 0..start {
                p.eval_hinted(breakpoints[breakpoints.len() - 1], &mut hint);
            }
            for x in &xs {
                prop_assert_eq!(bits(p.eval_hinted(*x, &mut hint)), bits(p.eval(x)), "{:?} x = {}", mode, x);
            }
        }
    }

    #[test]
    fn the_batch_equals_the_loop(
        (breakpoints, coefficients, degree) in pieces(),
        xs in queries(),
    ) {
        for mode in MODES {
            let p = PPoly::new(&breakpoints, &coefficients, degree)
                .expect("valid")
                .with_extrapolate(mode);
            let mut out = vec![0.0; xs.len()];
            prop_assert_eq!(p.eval_batch_into(&xs, &mut out, -7.0), Ok(()));
            for (x, got) in xs.iter().zip(&out) {
                prop_assert_eq!(got.to_bits(), p.eval_or(x, -7.0).to_bits(), "{:?} x = {}", mode, x);
            }
            let batch = p.eval_batch(&xs);
            for (x, got) in xs.iter().zip(batch) {
                prop_assert_eq!(bits(got), bits(p.eval(x)));
            }
        }
    }

    #[test]
    fn the_batch_equals_the_loop_in_every_order_and_past_one_run(
        (breakpoints, coefficients, degree) in pieces(),
        xs in prop::collection::vec(-30.0_f64..30.0, 1..300),
        nan_at in 0_usize..300,
    ) {
        // The batch works by runs of 64 points and uses a hint only on a run that
        // is in order and dense. Ascending, descending, and unordered batches,
        // longer than one run, with a NaN somewhere, take all of its paths.
        let mut ascending = xs.clone();
        ascending.sort_by(f64::total_cmp);
        let descending: Vec<f64> = ascending.iter().rev().copied().collect();
        // Dense and in order: many points inside the domain, a few intervals apart
        // at most.
        let (lo, hi) = (breakpoints[0], breakpoints[breakpoints.len() - 1]);
        let dense: Vec<f64> = (0..xs.len())
            .map(|i| lo + (hi - lo) * i as f64 / xs.len() as f64)
            .collect();
        let mut with_nan = ascending.clone();
        let at = nan_at % with_nan.len();
        with_nan[at] = f64::NAN;

        for mode in MODES {
            let p = PPoly::new(&breakpoints, &coefficients, degree)
                .expect("valid")
                .with_extrapolate(mode);
            for batch in [&xs, &ascending, &descending, &dense, &with_nan] {
                let mut out = vec![0.0; batch.len()];
                prop_assert_eq!(p.eval_batch_into(batch, &mut out, -7.0), Ok(()));
                for (x, got) in batch.iter().zip(&out) {
                    prop_assert_eq!(got.to_bits(), p.eval_or(x, -7.0).to_bits(), "{:?} x = {}", mode, x);
                }
            }
        }
    }

    #[test]
    fn inside_the_domain_the_value_is_horner_on_the_interval_of_the_point(
        (breakpoints, coefficients, degree) in pieces(),
        t in 0.0_f64..1.0,
    ) {
        let p = PPoly::new(&breakpoints, &coefficients, degree).expect("valid");
        let order = degree + 1;
        let last = breakpoints.len() - 2;
        for i in 0..=last {
            let piece = &coefficients[i * order..(i + 1) * order];
            let (left, right) = (breakpoints[i], breakpoints[i + 1]);

            // The left break point belongs to this interval: the constant term.
            prop_assert_eq!(bits(p.eval(&left)), bits(Some(horner(piece, 0.0))));

            // A point strictly inside.
            let x = left + t * (right - left);
            if left < x && x < right {
                prop_assert_eq!(bits(p.eval(&x)), bits(Some(horner(piece, x - left))));
            }
        }
        // The last break point belongs to the last interval.
        let piece = &coefficients[last * order..];
        let (left, right) = (breakpoints[last], breakpoints[last + 1]);
        prop_assert_eq!(bits(p.eval(&right)), bits(Some(horner(piece, right - left))));
    }

    #[test]
    fn outside_the_domain_each_mode_keeps_its_promise(
        (breakpoints, coefficients, degree) in pieces(),
        distance in 0.001_f64..50.0,
    ) {
        let p = PPoly::new(&breakpoints, &coefficients, degree).expect("valid");
        let order = degree + 1;
        let (lo, hi) = (breakpoints[0], breakpoints[breakpoints.len() - 1]);
        let (below, above) = (lo - distance, hi + distance);
        let first = &coefficients[..order];
        let last_left = breakpoints[breakpoints.len() - 2];
        let last = &coefficients[coefficients.len() - order..];

        let forbid = p.clone().with_extrapolate(Extrapolate::Forbid);
        prop_assert!(forbid.eval(&below).is_none());
        prop_assert!(forbid.eval(&above).is_none());

        let clamp = p.clone().with_extrapolate(Extrapolate::Clamp);
        prop_assert_eq!(bits(clamp.eval(&below)), bits(p.eval(&lo)));
        prop_assert_eq!(bits(clamp.eval(&above)), bits(p.eval(&hi)));

        let extend = p.clone().with_extrapolate(Extrapolate::Extend);
        prop_assert_eq!(bits(extend.eval(&below)), bits(Some(horner(first, below - lo))));
        prop_assert_eq!(bits(extend.eval(&above)), bits(Some(horner(last, above - last_left))));

        // A wrapped query is evaluated at some point of the domain.
        let periodic = p.with_extrapolate(Extrapolate::Periodic);
        let bound: f64 = coefficients.iter().map(|c| c.abs()).sum::<f64>()
            * (1.0_f64).max(hi - lo).powi(i32::try_from(degree).expect("small degree"));
        for x in [below, above] {
            let got = periodic.eval(&x).expect("a finite query has a place in the period");
            prop_assert!(got.is_finite() && got.abs() <= bound, "x = {}: {} over {}", x, got, bound);
        }
    }
}
