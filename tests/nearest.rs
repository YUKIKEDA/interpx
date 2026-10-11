//! `Nearest1D`: the golden data, the input checks through the constructor, which
//! point each kind returns, the four extrapolation modes, and the properties
//! that hold for any data.

mod common;

use common::{load, reals, Real};
use interpx::{Extrapolate, Hint, InterpolateError as E, Interpolator, Nearest1D, NearestKind};
use proptest::prelude::*;
use serde_json::Value as Json;

const MODES: [Extrapolate; 4] = [
    Extrapolate::Forbid,
    Extrapolate::Clamp,
    Extrapolate::Extend,
    Extrapolate::Periodic,
];

const KINDS: [NearestKind; 4] = [
    NearestKind::Nearest,
    NearestKind::NearestUp,
    NearestKind::Previous,
    NearestKind::Next,
];

fn bits(v: Option<f64>) -> Option<u64> {
    v.map(f64::to_bits)
}

fn nums(v: &Json) -> Vec<f64> {
    reals(&serde_json::from_value::<Vec<Real>>(v.clone()).unwrap_or_default())
}

/// A list as rows of one column, or a list of rows as it is.
fn rows(v: &Json) -> Vec<Vec<f64>> {
    match serde_json::from_value::<Vec<Vec<Real>>>(v.clone()) {
        Ok(rows) => rows.iter().map(|r| reals(r)).collect(),
        Err(_) => nums(v).into_iter().map(|x| vec![x]).collect(),
    }
}

fn kind_of(case: &Json) -> Option<NearestKind> {
    match case["options"]["kind"].as_str()? {
        "nearest" => Some(NearestKind::Nearest),
        "nearest-up" => Some(NearestKind::NearestUp),
        "previous" => Some(NearestKind::Previous),
        "next" => Some(NearestKind::Next),
        _ => None,
    }
}

// --- Golden data -------------------------------------------------------------

/// One case with `M` values per point. A data value is returned as it is, so
/// the comparison is bit for bit; SciPy has NaN where there is no value.
// Marked as test code so that `expect` is allowed here as in a `#[test]` function.
#[cfg(test)]
fn check<const M: usize>(case: &Json, kind: NearestKind, modes: &[Extrapolate]) {
    let name = case["name"].as_str().expect("name");
    let to_array = |r: &Vec<f64>| -> [f64; M] { r.as_slice().try_into().expect("M columns") };
    let y: Vec<[f64; M]> = rows(&case["y"]).iter().map(to_array).collect();
    let expected: Vec<[f64; M]> = rows(&case["expected"]).iter().map(to_array).collect();
    let queries = nums(&case["queries"]);
    assert_eq!(queries.len(), expected.len(), "{name}");
    for mode in modes {
        let f = Nearest1D::new(&nums(&case["x"]), &y, kind)
            .expect("valid data")
            .with_extrapolate(*mode);
        let mut out = vec![[0.0; M]; queries.len()];
        assert_eq!(f.eval_batch_into(&queries, &mut out, [f64::NAN; M]), Ok(()));
        for (k, (q, want)) in queries.iter().zip(&expected).enumerate() {
            let what = format!("{name} {mode:?} query {k} (x = {q})");
            match f.eval(q) {
                None => assert!(want.iter().all(|w| w.is_nan()), "{what}: None for {want:?}"),
                Some(got) => assert_eq!(got.map(f64::to_bits), want.map(f64::to_bits), "{what}"),
            }
            // NaN has many bit patterns; the batch wrote the one it was given.
            let wrote = out[k].map(|v| if v.is_nan() { f64::NAN } else { v });
            let wants = want.map(|v| if v.is_nan() { f64::NAN } else { v });
            assert_eq!(
                wrote.map(f64::to_bits),
                wants.map(f64::to_bits),
                "{what}: batch"
            );
        }
    }
}

#[test]
fn scipy_nearest_kinds_evaluate_to_the_golden_values() {
    let file = load::<Json>("oned").expect("golden file oned");
    let (mut inside, mut outside, mut wide) = (0, 0, 0);
    let mut kinds = std::collections::BTreeSet::new();
    for case in &file.cases {
        let Some(kind) = kind_of(case) else {
            continue;
        };
        // `fill_value="extrapolate"` gives the end values, as `Clamp` and
        // `Extend` both do. Otherwise SciPy was told to fill with NaN.
        let extrapolates = case["options"]["fill_value"] == "extrapolate";
        let modes: &[Extrapolate] = if extrapolates {
            outside += 1;
            &[Extrapolate::Clamp, Extrapolate::Extend]
        } else {
            &[Extrapolate::Forbid]
        };
        match rows(&case["y"]).first().map_or(0, Vec::len) {
            1 => {
                check::<1>(case, kind, modes);
                inside += usize::from(!extrapolates);
            }
            2 => {
                check::<2>(case, kind, modes);
                wide += 1;
            }
            other => panic!("values with {other} columns"),
        }
        kinds.insert(format!("{kind:?}"));
    }
    assert_eq!((inside, outside, wide), (12, 2, 4));
    assert_eq!(kinds.len(), 4);
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
        // SciPy takes one point.
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
        (
            &[0.0, f64::INFINITY],
            &[1.0, 2.0],
            E::NonFinite { index: 1 },
        ),
        (&[0.0, 1.0], &[1.0, f64::NAN], E::NonFinite { index: 1 }),
        // SciPy sorts, and takes a repeated coordinate.
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
    ];
    for kind in KINDS {
        assert!(Nearest1D::new(&[0.0, 1.0], &[1.0, 2.0], kind).is_ok());
        for (x, y, want) in cases {
            assert_eq!(Nearest1D::new(x, y, kind).err(), Some(want), "{kind:?}");
        }
        // The order of the checks: lengths, then the count, then NaN, then order.
        assert_eq!(
            Nearest1D::new(&[0.0], &[1.0, 2.0], kind).err(),
            Some(E::LengthMismatch {
                expected: 1,
                actual: 2
            })
        );
        assert_eq!(
            Nearest1D::new(&[f64::NAN], &[1.0], kind).err(),
            Some(E::InsufficientPoints {
                required: 2,
                actual: 1
            })
        );
        assert_eq!(
            Nearest1D::new(&[1.0, 0.0], &[f64::NAN, 1.0], kind).err(),
            Some(E::NonFinite { index: 0 })
        );
    }
}

// --- Which point each kind returns -----------------------------------------------

/// Points at 0, 1, 3 with values 0, 10, 30: the midpoints are 0.5 and 2.
fn at(kind: NearestKind, mode: Extrapolate, q: f64) -> Option<u64> {
    let f = Nearest1D::new(&[0.0, 1.0, 3.0], &[0.0, 10.0, 30.0], kind).ok()?;
    bits(f.with_extrapolate(mode).eval(&q))
}

#[test]
fn each_kind_returns_the_point_the_table_names() {
    use NearestKind::{Nearest, NearestUp, Next, Previous};
    // Query, then the value for Nearest, NearestUp, Previous, Next.
    let table = [
        (0.0, [0.0, 0.0, 0.0, 0.0]),
        (0.25, [0.0, 0.0, 0.0, 10.0]),
        // Halfway: Nearest goes left, NearestUp right.
        (0.5, [0.0, 10.0, 0.0, 10.0]),
        (0.75, [10.0, 10.0, 0.0, 10.0]),
        // A data point: every kind returns its value.
        (1.0, [10.0, 10.0, 10.0, 10.0]),
        (1.5, [10.0, 10.0, 10.0, 30.0]),
        (2.0, [10.0, 30.0, 10.0, 30.0]),
        (2.5, [30.0, 30.0, 10.0, 30.0]),
        (3.0, [30.0, 30.0, 30.0, 30.0]),
    ];
    for (q, want) in table {
        for (kind, w) in [Nearest, NearestUp, Previous, Next].into_iter().zip(want) {
            for mode in [Extrapolate::Forbid, Extrapolate::Clamp, Extrapolate::Extend] {
                assert_eq!(at(kind, mode, q), bits(Some(w)), "{kind:?} {mode:?} at {q}");
            }
        }
    }
}

#[test]
fn each_mode_treats_the_outside_as_documented() {
    for kind in KINDS {
        let f: Nearest1D<f64> =
            Nearest1D::new(&[0.0, 1.0, 3.0], &[0.0, 10.0, 30.0], kind).expect("valid");
        assert_eq!(
            f.extrapolate(),
            Extrapolate::Forbid,
            "SciPy raises by default"
        );
        assert_eq!(f.kind(), kind);
        assert_eq!(f.domain().0.to_bits(), 0.0_f64.to_bits());
        assert_eq!(f.domain().1.to_bits(), 3.0_f64.to_bits());

        assert_eq!(at(kind, Extrapolate::Forbid, -0.5), None);
        assert_eq!(at(kind, Extrapolate::Forbid, 3.5), None);
        // The end values on both sides, for every kind. (SciPy's `previous` has
        // NaN on the left and its `next` on the right.)
        for mode in [Extrapolate::Clamp, Extrapolate::Extend] {
            for below in [-0.5, -100.0, f64::NEG_INFINITY] {
                assert_eq!(at(kind, mode, below), bits(Some(0.0)), "{kind:?} {mode:?}");
            }
            for above in [3.5, 100.0, f64::INFINITY] {
                assert_eq!(at(kind, mode, above), bits(Some(30.0)), "{kind:?} {mode:?}");
            }
        }
        // The period is 3: 3.25 is 0.25, -0.25 is 2.75, and 3 itself is 0.
        let near = |q: f64| at(kind, Extrapolate::Periodic, q);
        assert_eq!(near(3.25), at(kind, Extrapolate::Forbid, 0.25), "{kind:?}");
        assert_eq!(near(-0.25), at(kind, Extrapolate::Forbid, 2.75), "{kind:?}");
        assert_eq!(near(3.0), bits(Some(0.0)), "{kind:?}");
        assert_eq!(near(f64::INFINITY), None, "{kind:?}");
        // SciPy returns the last value for NaN.
        for mode in MODES {
            assert_eq!(at(kind, mode, f64::NAN), None, "{kind:?} {mode:?}");
        }
    }
}

#[test]
fn the_midpoint_does_not_overflow() {
    // The sum of these two coordinates is over the largest number.
    let (x, y) = ([1.0e308, 1.7e308], [1.0, 2.0]);
    for (kind, q, want) in [
        (NearestKind::Nearest, 1.3e308, 1.0),
        (NearestKind::Nearest, 1.4e308, 2.0),
        (NearestKind::NearestUp, 1.3e308, 1.0),
        (NearestKind::NearestUp, 1.4e308, 2.0),
    ] {
        let f = Nearest1D::new(&x, &y, kind).expect("valid");
        assert_eq!(bits(f.eval(&q)), bits(Some(want)), "{kind:?} at {q}");
    }
}

#[test]
fn array_values_f32_and_thread_safety() {
    let f: Nearest1D<f32, [f32; 2]> =
        Nearest1D::new(&[0.0, 1.0], &[[0.0, 5.0], [1.0, 6.0]], NearestKind::Next).expect("valid");
    assert_eq!(
        f.eval(&0.5).map(|v| v.map(f32::to_bits)),
        Some([1.0_f32, 6.0].map(f32::to_bits))
    );
    assert_eq!(
        f.eval(&0.0).map(|v| v.map(f32::to_bits)),
        Some([0.0_f32, 5.0].map(f32::to_bits))
    );

    fn assert_send_sync<T: Send + Sync>() {}
    assert_send_sync::<Nearest1D<f64>>();
    assert_send_sync::<Nearest1D<f32, [f32; 3]>>();
    assert_send_sync::<NearestKind>();
}

#[test]
fn the_batch_rejects_a_length_mismatch_and_leaves_the_buffer() {
    let f: Nearest1D<f64> =
        Nearest1D::new(&[0.0, 1.0], &[0.0, 2.0], NearestKind::Nearest).expect("valid");
    for (n_points, n_out) in [(2, 3), (3, 2), (0, 1), (1, 0)] {
        let mut out = vec![7.0; n_out];
        assert_eq!(
            f.eval_batch_into(&vec![0.5; n_points], &mut out, 0.0),
            Err(E::LengthMismatch {
                expected: n_points,
                actual: n_out,
            })
        );
        assert!(out.iter().all(|v| v.to_bits() == 7.0_f64.to_bits()));
    }
}

// --- Properties ----------------------------------------------------------------

/// Coordinates that increase strictly. The value of point `i` is `i`, so that a
/// result says which point it came from.
fn coordinates() -> impl Strategy<Value = Vec<f64>> {
    (-5.0_f64..5.0, prop::collection::vec(0.01_f64..3.0, 1..12)).prop_map(|(start, widths)| {
        let mut x = vec![start];
        for w in widths {
            x.push(x[x.len() - 1] + w);
        }
        x
    })
}

proptest! {
    #[test]
    fn inside_the_data_the_point_returned_is_the_one_the_kind_names(
        x in coordinates(),
        t in prop::collection::vec(0.0_f64..=1.0, 1..40),
    ) {
        let y: Vec<f64> = (0..x.len()).map(|i| i as f64).collect();
        let (lo, hi) = (x[0], x[x.len() - 1]);
        for kind in KINDS {
            let f = Nearest1D::new(&x, &y, kind).expect("valid");
            // Every data point returns its own value.
            for (xi, yi) in x.iter().zip(&y) {
                prop_assert_eq!(bits(f.eval(xi)), bits(Some(*yi)), "{:?} at a point", kind);
            }
            for t in &t {
                let q = (lo + t * (hi - lo)).min(hi);
                let got = f.eval(&q).expect("inside");
                let k = got as usize;
                prop_assert_eq!((k as f64).to_bits(), got.to_bits(), "a data value");
                // No other point does better than the one returned.
                match kind {
                    NearestKind::Previous => {
                        prop_assert!(x[k] <= q, "{} is after {}", x[k], q);
                        prop_assert!(k + 1 == x.len() || x[k + 1] > q, "a later point is not after {}", q);
                    }
                    NearestKind::Next => {
                        prop_assert!(x[k] >= q, "{} is before {}", x[k], q);
                        prop_assert!(k == 0 || x[k - 1] < q, "an earlier point is not before {}", q);
                    }
                    _ => {
                        // SciPy's rule: against the halved sum of the neighbors.
                        let left_of = |i: usize| x[i] / 2.0 + x[i + 1] / 2.0;
                        if kind == NearestKind::Nearest {
                            prop_assert!(k + 1 == x.len() || q <= left_of(k), "{:?}: {} is past the midpoint", kind, q);
                            prop_assert!(k == 0 || q > left_of(k - 1), "{:?}: {} is not past the midpoint", kind, q);
                        } else {
                            prop_assert!(k + 1 == x.len() || q < left_of(k), "{:?}: {} is at or past the midpoint", kind, q);
                            prop_assert!(k == 0 || q >= left_of(k - 1), "{:?}: {} is before the midpoint", kind, q);
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn a_hint_and_the_batch_never_change_the_result(
        x in coordinates(),
        qs in prop::collection::vec(-12.0_f64..40.0, 1..150),
    ) {
        let y: Vec<f64> = (0..x.len()).map(|i| i as f64).collect();
        let mut ascending = qs.clone();
        ascending.sort_by(f64::total_cmp);
        let (lo, hi) = (x[0], x[x.len() - 1]);
        let dense: Vec<f64> = (0..=qs.len()).map(|i| lo + (hi - lo) * i as f64 / qs.len() as f64).chain([hi]).collect();
        for kind in KINDS {
            for mode in MODES {
                let f = Nearest1D::new(&x, &y, kind).expect("valid").with_extrapolate(mode);
                let mut hint = Hint::new();
                for q in qs.iter().chain([&f64::NAN, &hi, &lo]) {
                    prop_assert_eq!(bits(f.eval_hinted(*q, &mut hint)), bits(f.eval(q)), "{:?} {:?} x = {}", kind, mode, q);
                }
                for batch in [&qs, &ascending, &dense] {
                    let mut out = vec![0.0; batch.len()];
                    prop_assert_eq!(f.eval_batch_into(batch, &mut out, -7.0), Ok(()));
                    for (q, got) in batch.iter().zip(&out) {
                        prop_assert_eq!(got.to_bits(), f.eval_or(q, -7.0).to_bits(), "{:?} {:?} x = {}", kind, mode, q);
                    }
                }
            }
        }
    }

    #[test]
    fn outside_the_data_clamp_and_extend_give_the_end_values(
        x in coordinates(),
        distance in 0.001_f64..50.0,
    ) {
        let y: Vec<f64> = (0..x.len()).map(|i| i as f64).collect();
        let last = x.len() - 1;
        for kind in KINDS {
            let f = Nearest1D::new(&x, &y, kind).expect("valid");
            prop_assert_eq!(f.eval(&(x[0] - distance)), None);
            prop_assert_eq!(f.eval(&(x[last] + distance)), None);
            for mode in [Extrapolate::Clamp, Extrapolate::Extend] {
                let f = f.clone().with_extrapolate(mode);
                prop_assert_eq!(bits(f.eval(&(x[0] - distance))), bits(Some(y[0])), "{:?} {:?}", kind, mode);
                prop_assert_eq!(bits(f.eval(&(x[last] + distance))), bits(Some(y[last])), "{:?} {:?}", kind, mode);
            }
        }
    }
}
