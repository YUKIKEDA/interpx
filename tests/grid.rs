//! `RegularGridInterpolator` with `Nearest` and `Linear`: the golden data, the
//! input checks and their order, the layout of the values, the extrapolation
//! modes on a grid, and the properties.

mod common;

use common::{assert_close, load, reals, Real, DIRECT_F64};
use interpx::{
    AxisFault, Extrapolate, GridMethod, Hint, InterpolateError as E, Interpolator,
    RegularGridInterpolator as Grid,
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

fn bits(v: Option<f64>) -> Option<u64> {
    v.map(f64::to_bits)
}

// Marked as test code so that `expect` is allowed here as in a `#[test]` function
// (clippy.toml allows it in tests only). The same for the helpers below.
#[cfg(test)]
fn method_of(case: &Json) -> Option<GridMethod> {
    match case["options"]["method"].as_str().expect("method") {
        "nearest" => Some(GridMethod::Nearest),
        "linear" => Some(GridMethod::Linear),
        _ => None,
    }
}

/// Rebuilds one golden case in `D` dimensions and compares it with SciPy.
#[cfg(test)]
fn check<const D: usize>(case: &Json, method: GridMethod) {
    let name = case["name"].as_str().expect("name");
    let axes = rows(&case["axes"]);
    let slices: [&[f64]; D] = core::array::from_fn(|d| axes[d].as_slice());
    let built = Grid::new(slices, &nums(&case["values"]), method).expect("valid grid");
    assert_eq!(built.extrapolate(), Extrapolate::Forbid, "{name}: default");
    // `fill_value` is NaN (no value outside), or none (continued outside).
    let grid = if case["options"]["fill_value"].is_null() {
        built.with_extrapolate(Extrapolate::Extend)
    } else {
        built
    };
    let queries: Vec<[f64; D]> = rows(&case["queries"])
        .iter()
        .map(|r| r.as_slice().try_into().expect("D coordinates"))
        .collect();
    let expected = nums(&case["expected"]);
    assert_eq!(queries.len(), expected.len(), "{name}");
    let mut hints = [Hint::new(); D];
    let mut batch = vec![0.0; queries.len()];
    assert_eq!(grid.eval_batch_into(&queries, &mut batch, f64::NAN), Ok(()));
    for (k, (q, want)) in queries.iter().zip(&expected).enumerate() {
        let what = format!("{name} query {k} {q:?}");
        let got = grid.eval(q);
        match got {
            None => assert!(want.is_nan(), "{what}: None, SciPy has {want}"),
            // The values of the nearest point are the data, bit for bit.
            Some(got) if method == GridMethod::Nearest => {
                assert_eq!(got.to_bits(), want.to_bits(), "{what}");
            }
            Some(got) => assert_close(got, *want, DIRECT_F64, &what),
        }
        assert_eq!(bits(grid.eval_hinted(q, &mut hints)), bits(got), "{what}");
        assert_eq!(bits(Some(batch[k])), bits(Some(got.unwrap_or(f64::NAN))));
    }

    // Two values at every point: each column is its own interpolant.
    if let Some(columns) = case.get("values_columns") {
        let values: Vec<[f64; 2]> = rows(columns).iter().map(|r| [r[0], r[1]]).collect();
        let pair = Grid::new(slices, &values, method).expect("valid grid");
        for (k, (q, want)) in queries
            .iter()
            .zip(rows(&case["expected_columns"]))
            .enumerate()
        {
            match pair.eval(q) {
                None => assert!(want[0].is_nan() && want[1].is_nan(), "{name} query {k}"),
                Some(got) => {
                    for column in 0..2 {
                        let what = format!("{name} query {k} column {column}");
                        assert_close(got[column], want[column], DIRECT_F64, &what);
                    }
                    assert_eq!(bits(Some(got[0])), bits(grid.eval(q)), "{name} query {k}");
                }
            }
        }
    }
}

// --- Golden data -------------------------------------------------------------

#[test]
fn scipy_grids_evaluate_to_the_golden_values() {
    let file = load::<Json>("grid").expect("golden file grid");
    let mut counts = std::collections::BTreeMap::new();
    let (mut extended, mut columns, mut descending) = (0, 0, 0);
    for case in &file.cases {
        let Some(method) = method_of(case) else {
            continue;
        };
        let dimension = rows(&case["axes"]).len();
        match dimension {
            2 => check::<2>(case, method),
            3 => check::<3>(case, method),
            4 => check::<4>(case, method),
            other => panic!("a grid in {other} dimensions"),
        }
        *counts
            .entry((dimension, format!("{method:?}")))
            .or_insert(0) += 1;
        extended += usize::from(case["options"]["fill_value"].is_null());
        columns += usize::from(case.get("values_columns").is_some());
        let axes = rows(&case["axes"]);
        descending += usize::from(axes.iter().any(|a| a[0] > a[1]));
    }
    let counts: Vec<((usize, String), i32)> = counts.into_iter().collect();
    let want = [
        ((2, "Linear"), 12),
        ((2, "Nearest"), 14),
        ((3, "Linear"), 4),
        ((3, "Nearest"), 4),
        ((4, "Linear"), 1),
        ((4, "Nearest"), 1),
    ]
    .map(|((d, m), n)| ((d, m.to_owned()), n));
    assert_eq!(counts, want);
    assert_eq!((extended, columns, descending), (6, 2, 5));
}

// --- Input checks -------------------------------------------------------------------

fn on_axis(axis: usize, fault: AxisFault) -> E {
    E::GridAxis { axis, fault }
}

#[test]
fn the_constructor_reports_each_fault() {
    let (x, y) = ([0.0, 1.0, 2.0], [0.0, 1.0]);
    let new =
        |x: &[f64], y: &[f64], values: &[f64]| Grid::new([x, y], values, GridMethod::Linear).err();
    assert_eq!(new(&x, &y, &[0.0; 6]), None);
    assert_eq!(
        new(&x, &y, &[0.0; 5]),
        Some(E::LengthMismatch {
            expected: 6,
            actual: 5
        })
    );
    assert_eq!(
        new(&x, &[0.0], &[0.0; 3]),
        Some(on_axis(
            1,
            AxisFault::InsufficientPoints {
                required: 2,
                actual: 1
            }
        ))
    );
    assert_eq!(
        new(&[], &y, &[]),
        Some(on_axis(
            0,
            AxisFault::InsufficientPoints {
                required: 2,
                actual: 0
            }
        ))
    );
    assert_eq!(
        new(&[0.0, f64::NAN, 2.0], &y, &[0.0; 6]),
        Some(on_axis(0, AxisFault::NonFinite { index: 1 }))
    );
    assert_eq!(
        new(&x, &[0.0, f64::INFINITY], &[0.0; 6]),
        Some(on_axis(1, AxisFault::NonFinite { index: 1 }))
    );
    // Of several, the first: in the axis, among the axes, among the values.
    assert_eq!(
        new(&[f64::NAN, 1.0, f64::NAN], &[f64::NAN, 1.0], &[0.0; 6]),
        Some(on_axis(0, AxisFault::NonFinite { index: 0 }))
    );
    assert_eq!(
        new(&x, &y, &[0.0, f64::NAN, 0.0, f64::INFINITY, 0.0, 0.0]),
        Some(E::NonFinite { index: 1 })
    );
    assert_eq!(
        new(&x, &y, &[0.0, 0.0, 0.0, 0.0, f64::NAN, 0.0]),
        Some(E::NonFinite { index: 4 })
    );
    assert_eq!(
        new(&[0.0, 0.0, 2.0], &y, &[0.0; 6]),
        Some(on_axis(
            0,
            AxisFault::DuplicatePoints {
                first: 0,
                second: 1
            }
        ))
    );
    assert_eq!(
        new(&[0.0, 1.0, 1.0], &y, &[0.0; 6]),
        Some(on_axis(
            0,
            AxisFault::DuplicatePoints {
                first: 1,
                second: 2
            }
        ))
    );
    assert_eq!(
        new(&[0.0, 2.0, 1.0], &y, &[0.0; 6]),
        Some(on_axis(0, AxisFault::NotMonotonic { index: 2 }))
    );
    // A descending axis that turns up, and one that stalls.
    assert_eq!(
        new(&[2.0, 0.0, 1.0], &y, &[0.0; 6]),
        Some(on_axis(0, AxisFault::NotMonotonic { index: 2 }))
    );
    assert_eq!(
        new(&[2.0, 1.0, 1.0], &y, &[0.0; 6]),
        Some(on_axis(
            0,
            AxisFault::DuplicatePoints {
                first: 1,
                second: 2
            }
        ))
    );
    // The same checks for the method that needs no weights.
    let nearest = Grid::new([&x[..], &[0.0][..]], &[0.0; 3], GridMethod::Nearest).err();
    assert_eq!(
        nearest,
        Some(on_axis(
            1,
            AxisFault::InsufficientPoints {
                required: 2,
                actual: 1
            }
        ))
    );
}

#[test]
fn the_checks_come_in_the_order_of_the_design() {
    let new =
        |x: &[f64], y: &[f64], values: &[f64]| Grid::new([x, y], values, GridMethod::Linear).err();
    // The length before the count of an axis.
    assert_eq!(
        new(&[0.0], &[0.0, 1.0], &[0.0; 3]),
        Some(E::LengthMismatch {
            expected: 2,
            actual: 3
        })
    );
    // The count of a later axis before a NaN in an earlier one.
    assert_eq!(
        new(&[f64::NAN, 1.0], &[0.0], &[0.0; 2]),
        Some(on_axis(
            1,
            AxisFault::InsufficientPoints {
                required: 2,
                actual: 1
            }
        ))
    );
    // A NaN in a later axis before a NaN among the values.
    assert_eq!(
        new(&[0.0, 1.0], &[0.0, f64::NAN], &[f64::NAN; 4]),
        Some(on_axis(1, AxisFault::NonFinite { index: 1 }))
    );
    // A NaN among the values before the order of an axis.
    assert_eq!(
        new(&[1.0, 1.0], &[0.0, 1.0], &[0.0, 0.0, 0.0, f64::NAN]),
        Some(E::NonFinite { index: 3 })
    );
    // Of two axes out of order, the first.
    assert_eq!(
        new(&[0.0, 2.0, 1.0], &[5.0, 5.0], &[0.0; 6]),
        Some(on_axis(0, AxisFault::NotMonotonic { index: 2 }))
    );
    // A product of lengths that does not fit: the length of no slice.
    let long = vec![0.0; 1 << 17];
    let axes: [&[f64]; 4] = [&long, &long, &long, &long];
    assert_eq!(
        Grid::new(axes, &[0.0], GridMethod::Linear).err(),
        Some(E::LengthMismatch {
            expected: usize::MAX,
            actual: 1
        })
    );
}

// --- The layout and the two methods --------------------------------------------------

#[test]
fn values_are_in_c_order_with_the_last_axis_fastest() {
    // The value at (i, j, k) says where it is: 100 i + 10 j + k.
    let (a, b, c) = ([0.0, 1.0], [0.0, 1.0, 2.0], [0.0, 1.0, 2.0, 3.0]);
    let values: Vec<f64> = (0..24)
        .map(|n| f64::from(100 * (n / 12) + 10 * (n / 4 % 3) + n % 4))
        .collect();
    let grid = Grid::new([&a[..], &b[..], &c[..]], &values, GridMethod::Nearest).expect("valid");
    assert_eq!(grid.shape(), [2, 3, 4]);
    assert_eq!(bits(grid.eval(&[1.0, 2.0, 3.0])), bits(Some(123.0)));
    assert_eq!(bits(grid.eval(&[0.0, 1.0, 2.0])), bits(Some(12.0)));
    assert_eq!(bits(grid.eval(&[0.9, 0.2, 2.6])), bits(Some(103.0)));
    let linear = Grid::new([&a[..], &b[..], &c[..]], &values, GridMethod::Linear).expect("valid");
    // The values are linear in each index, so the interpolant is exact.
    assert_close(
        linear.eval(&[0.5, 1.25, 2.5]).expect("inside"),
        50.0 + 12.5 + 2.5,
        DIRECT_F64,
        "inside a cell",
    );
    assert_eq!(linear.method(), GridMethod::Linear);
    assert_eq!(grid.method(), GridMethod::Nearest);
    let (low, high) = grid.domain();
    assert_eq!(low.map(f64::to_bits), [0.0_f64; 3].map(f64::to_bits));
    assert_eq!(
        high.map(f64::to_bits),
        [1.0_f64, 2.0, 3.0].map(f64::to_bits)
    );
}

#[test]
fn the_corners_are_added_with_the_first_axis_outermost() {
    // At the center of a cell every corner weighs a quarter. The values are
    // chosen so that the order of the sum shows: (0, 0), (0, 1), (1, 0), (1, 1)
    // gives 2.5e15 + 0.25 (lost in the rounding) - 2.5e15 + 0.25 = 0.25, while
    // the two large values first would give 0.5.
    let side = [0.0, 1.0];
    let values = [1e16, 1.0, -1e16, 1.0];
    let grid = Grid::new([&side[..]; 2], &values, GridMethod::Linear).expect("valid");
    assert_eq!(bits(grid.eval(&[0.5, 0.5])), bits(Some(0.25)));
}

#[test]
fn nearest_takes_the_lower_point_halfway() {
    let x = [0.0, 2.0, 3.0];
    let y = [0.0, 1.0];
    let values = [0.0, 1.0, 10.0, 11.0, 20.0, 21.0];
    let grid = Grid::new([&x[..], &y[..]], &values, GridMethod::Nearest).expect("valid");
    // Halfway in x, then in y, then in both; and a rounding past halfway.
    assert_eq!(bits(grid.eval(&[1.0, 0.0])), bits(Some(0.0)));
    assert_eq!(bits(grid.eval(&[2.0, 0.5])), bits(Some(10.0)));
    assert_eq!(bits(grid.eval(&[2.5, 0.5])), bits(Some(10.0)));
    assert_eq!(
        bits(grid.eval(&[1.0 + 1e-12, 0.5 + 1e-12])),
        bits(Some(11.0))
    );
    // The upper ends belong to the grid.
    assert_eq!(bits(grid.eval(&[3.0, 1.0])), bits(Some(21.0)));
    // A descending axis is put in order first: halfway between its points 2
    // and 0 the lower one, 0, is taken, which was given second.
    let falling = [3.0, 2.0, 0.0];
    let flipped = [20.0, 21.0, 10.0, 11.0, 0.0, 1.0];
    let grid = Grid::new([&falling[..], &y[..]], &flipped, GridMethod::Nearest).expect("valid");
    assert_eq!(bits(grid.eval(&[1.0, 0.0])), bits(Some(0.0)));
    assert_eq!(bits(grid.eval(&[2.5, 1.0])), bits(Some(11.0)));
}

#[test]
fn a_descending_axis_is_kept_ascending_with_its_values() {
    // Three axes of 2, 3, and 4 points; the middle one and the last descend.
    let (a, b, c) = ([0.0, 1.0], [5.0, 3.0, 2.0], [9.0, 8.0, 6.0, 5.0]);
    let values: Vec<f64> = (0..24).map(f64::from).collect();
    let grid = Grid::new([&a[..], &b[..], &c[..]], &values, GridMethod::Linear).expect("valid");
    assert_eq!(grid.axes()[1], [2.0, 3.0, 5.0]);
    assert_eq!(grid.axes()[2], [5.0, 6.0, 8.0, 9.0]);
    // Every grid point keeps its value: the point (i, j, k) as given has the
    // value 12 i + 4 j + k.
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            for (k, z) in c.iter().enumerate() {
                let want = f64::from(u8::try_from(12 * i + 4 * j + k).expect("small"));
                assert_eq!(
                    bits(grid.eval(&[*x, *y, *z])),
                    bits(Some(want)),
                    "{i} {j} {k}"
                );
            }
        }
    }
    // What is kept: the values of the grid in ascending order.
    let kept = grid.values();
    assert_eq!(kept.len(), 24);
    assert_eq!(kept[0].to_bits(), 11.0_f64.to_bits(), "(0, 2, 3) as given");
    assert_eq!(kept[23].to_bits(), 12.0_f64.to_bits(), "(1, 0, 0) as given");
    // Taken over, the same.
    let owned = Grid::from_owned(
        [a.to_vec(), b.to_vec(), c.to_vec()],
        values,
        GridMethod::Linear,
    )
    .expect("valid");
    assert_eq!(owned.values(), kept);
    assert_eq!(owned.axes(), grid.axes());
}

#[test]
fn each_mode_does_its_part_outside_the_grid() {
    // 1 + x + 10 y on [0, 2] x [0, 1], which the linear method reproduces.
    let (x, y) = ([0.0, 1.0, 2.0], [0.0, 1.0]);
    let values: Vec<f64> = x
        .iter()
        .flat_map(|x| y.iter().map(move |y| 1.0 + x + 10.0 * y))
        .collect();
    let plane = |x: f64, y: f64| 1.0 + x + 10.0 * y;
    for method in [GridMethod::Linear, GridMethod::Nearest] {
        let grid = Grid::new([&x[..], &y[..]], &values, method).expect("valid");
        let with = |mode| grid.clone().with_extrapolate(mode);
        // One coordinate outside is enough to be outside.
        let outside = [[3.0, 0.5], [1.0, -0.5], [-1.0, 2.0]];
        for q in outside {
            assert_eq!(grid.eval(&q), None, "{method:?} {q:?}");
            assert_eq!(bits(Some(grid.eval_or(&q, -7.0))), bits(Some(-7.0)));
        }
        assert_eq!(with(Extrapolate::Forbid).eval(&[f64::NAN, 0.5]), None);
        for mode in [
            Extrapolate::Clamp,
            Extrapolate::Extend,
            Extrapolate::Periodic,
        ] {
            assert_eq!(with(mode).eval(&[1.0, f64::NAN]), None, "{mode:?}");
            assert_eq!(with(mode).extrapolate(), mode);
        }
        // Clamp: each coordinate moves to the nearest end of its axis.
        let clamp = with(Extrapolate::Clamp);
        assert_eq!(bits(clamp.eval(&[3.0, 0.5])), bits(grid.eval(&[2.0, 0.5])));
        assert_eq!(bits(clamp.eval(&[-1.0, 2.0])), bits(grid.eval(&[0.0, 1.0])));
        assert_eq!(
            bits(clamp.eval(&[f64::INFINITY, f64::NEG_INFINITY])),
            bits(grid.eval(&[2.0, 0.0]))
        );
        // Periodic: each coordinate wraps into its axis; the upper end is the
        // lower end again.
        let periodic = with(Extrapolate::Periodic);
        assert_eq!(
            bits(periodic.eval(&[2.5, 1.25])),
            bits(grid.eval(&[0.5, 0.25]))
        );
        assert_eq!(
            bits(periodic.eval(&[-0.5, 0.5])),
            bits(grid.eval(&[1.5, 0.5]))
        );
        assert_eq!(
            bits(periodic.eval(&[2.0, 1.0])),
            bits(grid.eval(&[0.0, 0.0]))
        );
        assert_eq!(periodic.eval(&[f64::INFINITY, 0.5]), None);
        // Extend: the linear method continues its plane, the nearest method
        // stays at the outermost point.
        let extend = with(Extrapolate::Extend);
        let far = [3.5, -0.25];
        let want = match method {
            GridMethod::Linear => plane(3.5, -0.25),
            _ => plane(2.0, 0.0),
        };
        assert_close(
            extend.eval(&far).expect("extends"),
            want,
            DIRECT_F64,
            "extend",
        );
        if method == GridMethod::Nearest {
            assert_eq!(bits(extend.eval(&far)), bits(clamp.eval(&far)));
        }
    }
}

#[test]
fn a_batch_gives_the_bits_of_eval_along_a_path_and_off_it() {
    // 30 by 20 points, unevenly spaced.
    let x: Vec<f64> = (0..30)
        .map(|i| f64::from(i) + 0.3 * f64::from(i % 3))
        .collect();
    let y: Vec<f64> = (0..20).map(|j| 0.5 * f64::from(j * j)).collect();
    let values: Vec<f64> = (0..600).map(|n| f64::from(n % 17) - 8.0).collect();
    let (x_end, y_end) = (x[29], y[19]);
    // A path from corner to corner and a little past it, the same backward, a
    // path that turns around halfway, and points in no order. 200 points each:
    // several runs of the batch.
    let forward: Vec<[f64; 2]> = (0..200)
        .map(|k| {
            let t = f64::from(k) / 190.0;
            [t * x_end, t * y_end]
        })
        .collect();
    let backward: Vec<[f64; 2]> = forward.iter().rev().copied().collect();
    let turning: Vec<[f64; 2]> = forward
        .iter()
        .enumerate()
        .map(|(k, p)| {
            if k < 100 {
                *p
            } else {
                [p[0], forward[199 - k][1]]
            }
        })
        .collect();
    let scattered: Vec<[f64; 2]> = (0..200_u32)
        .map(|k| {
            [
                f64::from(k * 37 % 200) / 200.0 * x_end,
                f64::from(k * 91 % 200) / 200.0 * y_end,
            ]
        })
        .collect();
    let mut with_nan = forward.clone();
    with_nan[70][1] = f64::NAN;
    for method in [GridMethod::Linear, GridMethod::Nearest] {
        for mode in [
            Extrapolate::Forbid,
            Extrapolate::Extend,
            Extrapolate::Periodic,
        ] {
            let grid = Grid::new([&x[..], &y[..]], &values, method)
                .expect("valid")
                .with_extrapolate(mode);
            for points in [&forward, &backward, &turning, &scattered, &with_nan] {
                let mut out = vec![0.0; points.len()];
                assert_eq!(grid.eval_batch_into(points, &mut out, -7.0), Ok(()));
                for (point, got) in points.iter().zip(&out) {
                    let want = grid.eval(point).unwrap_or(-7.0);
                    assert_eq!(
                        got.to_bits(),
                        want.to_bits(),
                        "{method:?} {mode:?} {point:?}"
                    );
                }
            }
            // A buffer of another length is refused, and left alone.
            let mut short = vec![1.5; 199];
            assert_eq!(
                grid.eval_batch_into(&forward, &mut short, -7.0),
                Err(E::LengthMismatch {
                    expected: 200,
                    actual: 199
                })
            );
            assert!(short.iter().all(|v| v.to_bits() == 1.5_f64.to_bits()));
        }
    }
}

#[test]
fn f32_one_dimension_columns_and_thread_safety() {
    let x = [0.0_f32, 1.0, 3.0];
    let line: Grid<f32, 1> =
        Grid::new([&x[..]], &[0.0, 2.0, 6.0], GridMethod::Linear).expect("valid");
    assert_eq!(line.eval(&[2.0]).map(f32::to_bits), Some(4.0_f32.to_bits()));
    let pairs = [[0.0_f32, 1.0], [2.0, 1.0], [6.0, 5.0]];
    let both: Grid<f32, 1, [f32; 2]> =
        Grid::new([&x[..]], &pairs, GridMethod::Linear).expect("valid");
    assert_eq!(
        both.eval(&[2.0]).map(|v| v.map(f32::to_bits)),
        Some([4.0_f32, 3.0].map(f32::to_bits))
    );
    let nan = Grid::new(
        [&x[..]],
        &[[0.0, 0.0], [0.0, f32::NAN], [0.0, 0.0]],
        GridMethod::Nearest,
    );
    assert_eq!(nan.err(), Some(E::NonFinite { index: 1 }));

    fn assert_send_sync<T: Send + Sync>() {}
    assert_send_sync::<Grid<f64, 3>>();
    assert_send_sync::<Grid<f32, 2, [f32; 4]>>();
}

// --- Properties ----------------------------------------------------------------

/// An axis of 2 to 6 points, increasing.
fn axis() -> impl Strategy<Value = Vec<f64>> {
    (-3.0_f64..3.0, prop::collection::vec(0.1_f64..2.0, 1..6)).prop_map(|(start, widths)| {
        let mut x = vec![start];
        for w in widths {
            x.push(x[x.len() - 1] + w);
        }
        x
    })
}

/// Three axes, values for their grid, and points given as fractions of the box
/// stretched by half its size on each side.
type Case = ([Vec<f64>; 3], Vec<f64>, Vec<[f64; 3]>);

fn cases() -> impl Strategy<Value = Case> {
    (
        [axis(), axis(), axis()],
        prop::collection::vec(-10.0_f64..10.0, 216),
        prop::collection::vec(prop::array::uniform3(-0.5_f64..1.5), 1..30),
    )
        .prop_map(|(axes, values, fractions)| {
            let count = axes.iter().map(Vec::len).product();
            let points = fractions
                .iter()
                .map(|f| {
                    core::array::from_fn(|d| {
                        let (lo, hi) = (axes[d][0], axes[d][axes[d].len() - 1]);
                        lo + f[d] * (hi - lo)
                    })
                })
                .collect();
            (axes, values[..count].to_vec(), points)
        })
}

fn slices(axes: &[Vec<f64>; 3]) -> [&[f64]; 3] {
    [&axes[0], &axes[1], &axes[2]]
}

fn methods() -> impl Strategy<Value = GridMethod> {
    prop::sample::select(vec![GridMethod::Nearest, GridMethod::Linear])
}

fn modes() -> impl Strategy<Value = Extrapolate> {
    prop::sample::select(vec![
        Extrapolate::Forbid,
        Extrapolate::Clamp,
        Extrapolate::Extend,
        Extrapolate::Periodic,
    ])
}

proptest! {
    #[test]
    fn the_values_are_reproduced_at_the_grid_points(
        (axes, values, _) in cases(),
        method in methods(),
        mode in modes(),
    ) {
        let grid = Grid::new(slices(&axes), &values, method).expect("valid").with_extrapolate(mode);
        let shape = grid.shape();
        for i in 0..shape[0] {
            for j in 0..shape[1] {
                for k in 0..shape[2] {
                    let want = values[(i * shape[1] + j) * shape[2] + k];
                    let at = [axes[0][i], axes[1][j], axes[2][k]];
                    // A period holds its lower end and not its upper one.
                    let wrapped = mode == Extrapolate::Periodic
                        && (i + 1 == shape[0] || j + 1 == shape[1] || k + 1 == shape[2]);
                    if !wrapped {
                        prop_assert_eq!(bits(grid.eval(&at)), bits(Some(want)), "{:?} at {:?}", method, at);
                    }
                }
            }
        }
    }

    #[test]
    fn a_function_linear_in_each_coordinate_is_reproduced(
        (axes, _, points) in cases(),
        c in prop::array::uniform8(-2.0_f64..2.0),
    ) {
        // Every product of the coordinates, each to the power 0 or 1.
        let f = |p: &[f64; 3]| {
            let (x, y, z) = (p[0], p[1], p[2]);
            c[0] + c[1] * x + c[2] * y + c[3] * z + c[4] * x * y + c[5] * y * z + c[6] * x * z + c[7] * x * y * z
        };
        let mut values = Vec::new();
        for x in &axes[0] {
            for y in &axes[1] {
                for z in &axes[2] {
                    values.push(f(&[*x, *y, *z]));
                }
            }
        }
        let grid = Grid::new(slices(&axes), &values, GridMethod::Linear).expect("valid").with_extrapolate(Extrapolate::Extend);
        for p in &points {
            let (got, want) = (grid.eval(p).expect("extends"), f(p));
            prop_assert!((got - want).abs() <= 1e-9 * want.abs().max(100.0), "at {:?}: {} for {}", p, got, want);
        }
    }

    #[test]
    fn linear_stays_between_the_corners_and_nearest_returns_a_value(
        (axes, values, points) in cases(),
    ) {
        let (lo, hi) = values.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), v| (lo.min(*v), hi.max(*v)));
        let linear = Grid::new(slices(&axes), &values, GridMethod::Linear).expect("valid").with_extrapolate(Extrapolate::Clamp);
        let nearest = Grid::new(slices(&axes), &values, GridMethod::Nearest).expect("valid").with_extrapolate(Extrapolate::Extend);
        for p in &points {
            let got = linear.eval(p).expect("clamped");
            prop_assert!(lo - 1e-12 <= got && got <= hi + 1e-12, "at {:?}: {} outside [{}, {}]", p, got, lo, hi);
            let got = nearest.eval(p).expect("extends");
            prop_assert!(values.iter().any(|v| v.to_bits() == got.to_bits()), "at {:?}: {} is no value of the grid", p, got);
        }
    }

    #[test]
    fn a_descending_axis_gives_the_grid_of_the_ascending_one(
        (axes, values, points) in cases(),
        method in methods(),
        mode in modes(),
        flip in prop::array::uniform3(any::<bool>()),
    ) {
        let grid = Grid::new(slices(&axes), &values, method).expect("valid").with_extrapolate(mode);
        // The same grid with some axes given backward, and the values to match.
        let shape = grid.shape();
        let given: [Vec<f64>; 3] = core::array::from_fn(|d| {
            let mut axis = axes[d].clone();
            if flip[d] {
                axis.reverse();
            }
            axis
        });
        let mut moved = vec![0.0; values.len()];
        for i in 0..shape[0] {
            for j in 0..shape[1] {
                for k in 0..shape[2] {
                    let from = [i, j, k];
                    let to: [usize; 3] = core::array::from_fn(|d| if flip[d] { shape[d] - 1 - from[d] } else { from[d] });
                    moved[(to[0] * shape[1] + to[1]) * shape[2] + to[2]] = values[(i * shape[1] + j) * shape[2] + k];
                }
            }
        }
        let other = Grid::new(slices(&given), &moved, method).expect("valid").with_extrapolate(mode);
        prop_assert_eq!(other.axes(), grid.axes());
        for (a, b) in other.values().iter().zip(grid.values()) {
            prop_assert_eq!(a.to_bits(), b.to_bits());
        }
        for p in &points {
            prop_assert_eq!(bits(other.eval(p)), bits(grid.eval(p)));
        }
    }

    #[test]
    fn clamp_and_periodic_evaluate_at_the_moved_point(
        (axes, values, points) in cases(),
        method in methods(),
    ) {
        let grid = Grid::new(slices(&axes), &values, method).expect("valid");
        let clamp = grid.clone().with_extrapolate(Extrapolate::Clamp);
        let periodic = grid.clone().with_extrapolate(Extrapolate::Periodic);
        let (lo, hi) = grid.domain();
        for p in &points {
            let clamped: [f64; 3] = core::array::from_fn(|d| p[d].clamp(lo[d], hi[d]));
            prop_assert_eq!(bits(clamp.eval(p)), bits(grid.eval(&clamped)), "clamp at {:?}", p);
            // One period on, the same value.
            let shifted: [f64; 3] = core::array::from_fn(|d| p[d] + (hi[d] - lo[d]));
            let (here, there) = (periodic.eval(p).expect("periodic"), periodic.eval(&shifted).expect("periodic"));
            if method == GridMethod::Linear {
                prop_assert!((here - there).abs() <= 1e-9 * here.abs().max(10.0), "periodic at {:?}: {} and {}", p, here, there);
            }
            // Inside the grid, away from the upper faces, the mode changes nothing.
            if (0..3).all(|d| lo[d] <= p[d] && p[d] < hi[d]) {
                prop_assert_eq!(bits(periodic.eval(p)), bits(grid.eval(p)));
                prop_assert_eq!(bits(clamp.eval(p)), bits(grid.eval(p)));
            }
        }
    }

    #[test]
    fn columns_hints_and_batches_agree_with_eval(
        (axes, values, points) in cases(),
        method in methods(),
        mode in modes(),
        scale in -3.0_f64..3.0,
    ) {
        let grid = Grid::new(slices(&axes), &values, method).expect("valid").with_extrapolate(mode);
        let second: Vec<f64> = values.iter().rev().map(|v| v * scale + 1.0).collect();
        let both: Vec<[f64; 2]> = values.iter().zip(&second).map(|(a, b)| [*a, *b]).collect();
        let pair = Grid::new(slices(&axes), &both, method).expect("valid").with_extrapolate(mode);
        let other = Grid::new(slices(&axes), &second, method).expect("valid").with_extrapolate(mode);
        let mut hints = [Hint::new(); 3];
        let mut out = vec![0.0; points.len()];
        prop_assert_eq!(grid.eval_batch_into(&points, &mut out, -7.0), Ok(()));
        for (p, slot) in points.iter().zip(&out) {
            let got = grid.eval(p);
            prop_assert_eq!(bits(grid.eval_hinted(p, &mut hints)), bits(got));
            prop_assert_eq!(bits(Some(*slot)), bits(Some(got.unwrap_or(-7.0))));
            let columns = pair.eval(p);
            prop_assert_eq!(bits(columns.map(|c| c[0])), bits(got));
            prop_assert_eq!(bits(columns.map(|c| c[1])), bits(other.eval(p)));
        }
    }
}
