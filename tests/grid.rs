//! `RegularGridInterpolator`: the golden data of every method, the input checks
//! and their order, the layout of the values, the extrapolation modes on a
//! grid, the spline methods against the splines in one variable, the
//! shape-preserving method against the interpolant in one variable, and the
//! properties.

mod common;

use common::{assert_close, load, reals, Real, DIRECT_F64};
use interpx::{
    AxisFault, BSpline1D, BoundaryCondition, Extrapolate, GridMethod, Hint, InterpolateError as E,
    Interpolator, PchipInterpolator, RegularGridInterpolator as Grid,
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
fn method_of(case: &Json) -> GridMethod {
    match case["options"]["method"].as_str().expect("method") {
        "nearest" => GridMethod::Nearest,
        // SciPy's two names for one function.
        "linear" | "slinear" => GridMethod::Linear,
        "cubic" => GridMethod::Cubic,
        "quintic" => GridMethod::Quintic,
        "pchip" => GridMethod::Pchip,
        other => panic!("a method of no grid: {other}"),
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

/// How many cases a golden file has of each dimension and method.
type Counts = Vec<((usize, String), i32)>;

/// Compares every case of a golden file with SciPy. Returns the counts of the
/// cases, and how many are continued outside the grid, have two values a
/// point, and have an axis that descends.
#[cfg(test)]
fn check_file(name: &str) -> (Counts, [usize; 3]) {
    let file = load::<Json>(name).expect("golden file");
    let mut counts = std::collections::BTreeMap::new();
    let (mut extended, mut columns, mut descending) = (0, 0, 0);
    for case in &file.cases {
        let method = method_of(case);
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
    (
        counts.into_iter().collect(),
        [extended, columns, descending],
    )
}

#[test]
fn scipy_grids_evaluate_to_the_golden_values() {
    let (counts, [extended, columns, descending]) = check_file("grid");
    // `Linear` takes SciPy's `linear` and `slinear` cases.
    let want = [
        ((2, "Cubic"), 13),
        ((2, "Linear"), 21),
        ((2, "Nearest"), 14),
        ((2, "Pchip"), 9),
        ((2, "Quintic"), 13),
        ((3, "Cubic"), 4),
        ((3, "Linear"), 6),
        ((3, "Nearest"), 4),
        ((3, "Pchip"), 2),
        ((3, "Quintic"), 4),
        ((4, "Cubic"), 1),
        ((4, "Linear"), 2),
        ((4, "Nearest"), 1),
        ((4, "Pchip"), 1),
        ((4, "Quintic"), 1),
    ]
    .map(|((d, m), n)| ((d, m.to_owned()), n));
    assert_eq!(counts, want);
    assert_eq!((extended, columns, descending), (18, 4, 7));
}

#[test]
fn scipy_pchip_grids_evaluate_to_the_golden_values() {
    // The cases chosen for the method: values in no order, steps, grids
    // continued outside, two and three points on an axis (where the reference
    // is SciPy's interpolant in one variable applied by hand), a descending
    // axis, two values a point.
    let (counts, [extended, columns, descending]) = check_file("grid_pchip");
    let want = [((2, "Pchip"), 11), ((3, "Pchip"), 7)].map(|((d, m), n)| ((d, m.to_owned()), n));
    assert_eq!(counts, want);
    assert_eq!((extended, columns, descending), (12, 1, 1));
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
    for method in [GridMethod::Linear, GridMethod::Nearest, GridMethod::Pchip] {
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

// --- The spline methods ----------------------------------------------------------

fn splines() -> [(GridMethod, usize); 2] {
    [(GridMethod::Cubic, 3), (GridMethod::Quintic, 5)]
}

#[test]
fn a_spline_method_takes_one_point_more_than_its_degree_on_every_axis() {
    let x: Vec<f64> = (0..8).map(f64::from).collect();
    for (method, degree) in splines() {
        let fewest = degree + 1;
        let values = vec![0.0; 8 * fewest];
        assert!(Grid::new([&x[..], &x[..fewest]], &values, method).is_ok());
        assert_eq!(
            Grid::new([&x[..], &x[..fewest]], &values, method).map(|g| g.method()),
            Ok(method)
        );
        let short = Grid::new([&x[..], &x[..degree]], &values[..8 * degree], method).err();
        assert_eq!(
            short,
            Some(on_axis(
                1,
                AxisFault::InsufficientPoints {
                    required: fewest,
                    actual: degree
                }
            ))
        );
        // The count of an axis still comes after the length of the values.
        assert_eq!(
            Grid::new([&x[..], &x[..degree]], &values, method).err(),
            Some(E::LengthMismatch {
                expected: 8 * degree,
                actual: 8 * fewest
            })
        );
    }
}

#[test]
fn what_cannot_be_solved_is_an_error_and_not_a_panic() {
    // The second axis reaches from one end of the numbers to the other: the
    // distances between its points overflow, and so do the splines on it.
    let x: Vec<f64> = (0..6).map(f64::from).collect();
    let wide = [-1.7e308, -1.0e308, -1.0, 1.0, 1.0e308, 1.7e308];
    for (method, _) in splines() {
        let err = Grid::new([&x[..], &wide[..]], &[0.0; 36], method).err();
        assert!(
            matches!(
                err,
                Some(E::GridAxis {
                    axis: 1,
                    fault: AxisFault::NonFinite { .. }
                })
            ),
            "{method:?}: {err:?}"
        );
        // Values near the largest number, alternating in sign: the
        // coefficients of the spline through them are larger.
        let values: Vec<f64> = (0..36)
            .map(|n| {
                if (n / 6 + n % 6) % 2 == 0 {
                    1.7e308
                } else {
                    -1.7e308
                }
            })
            .collect();
        let err = Grid::new([&x[..], &x[..]], &values, method).err();
        assert!(
            matches!(err, Some(E::NonFinite { .. })),
            "{method:?}: {err:?}"
        );
        // The place of the coefficient is that in the values as given: with
        // both axes given backward and the values to match, point (i, j) is
        // given as (5 - i, 5 - j), which is place 35 - n for place n.
        let Some(E::NonFinite { index }) = err else {
            panic!("{method:?}: {err:?}");
        };
        let backward: Vec<f64> = x.iter().rev().copied().collect();
        let moved: Vec<f64> = values.iter().rev().copied().collect();
        assert_eq!(
            Grid::new([&backward[..], &backward[..]], &moved, method).err(),
            Some(E::NonFinite { index: 35 - index }),
            "{method:?}"
        );
        // One axis backward: the row of the matrix that could not be set up
        // is named as the point was given.
        let Some(E::GridAxis {
            axis: 1,
            fault: AxisFault::NonFinite { index: row },
        }) = Grid::new([&x[..], &wide[..]], &[0.0; 36], method).err()
        else {
            panic!("{method:?}");
        };
        let turned: Vec<f64> = wide.iter().rev().copied().collect();
        assert_eq!(
            Grid::new([&x[..], &turned[..]], &[0.0; 36], method).err(),
            Some(on_axis(1, AxisFault::NonFinite { index: 5 - row })),
            "{method:?}"
        );
        // The methods that solve nothing take both.
        assert!(Grid::new([&x[..], &x[..]], &values, GridMethod::Linear).is_ok());
    }
    assert!(Grid::new([&x[..], &wide[..]], &[0.0; 36], GridMethod::Nearest).is_ok());
}

#[test]
fn in_one_dimension_a_spline_method_is_the_spline_through_the_points() {
    let x = [0.0, 0.5, 1.5, 2.0, 3.5, 4.0, 5.5, 6.0, 7.0];
    let y = [1.0, -1.0, 0.5, 2.0, 0.0, -2.0, 1.5, 0.5, 3.0];
    for (method, degree) in splines() {
        let grid = Grid::new([&x[..]], &y, method)
            .expect("valid")
            .with_extrapolate(Extrapolate::Extend);
        let spline = BSpline1D::new(&x, &y, degree, BoundaryCondition::NotAKnot).expect("valid");
        // The same knots, the same system, the same sum: the same bits.
        for k in 0..=160 {
            let q = -0.5 + 0.05 * f64::from(k);
            assert_eq!(
                bits(grid.eval(&[q])),
                bits(spline.eval(&q)),
                "{method:?} at {q}"
            );
        }
    }
}

#[test]
fn a_tensor_product_polynomial_of_the_degree_is_reproduced() {
    // Unevenly spaced axes of 6, 7, and 8 points.
    let a = [0.0, 0.3, 1.0, 1.4, 2.5, 3.0];
    let b = [-1.0, -0.5, 0.25, 1.0, 1.5, 2.75, 3.0];
    let c = [0.0, 0.5, 1.0, 2.0, 2.25, 3.0, 4.0, 4.5];
    for (method, degree) in splines() {
        let power = i32::try_from(degree).expect("small");
        // A product of polynomials of the full degree in each coordinate.
        let f = |p: &[f64; 3]| {
            (1.0 + p[0] - 0.2 * p[0].powi(power))
                * (2.0 - 0.5 * p[1] * p[1] + 0.05 * p[1].powi(power))
                * (0.5 + 0.1 * p[2].powi(power - 1) - 0.01 * p[2].powi(power))
        };
        let mut values = Vec::new();
        for x in a {
            for y in b {
                for z in c {
                    values.push(f(&[x, y, z]));
                }
            }
        }
        let grid = Grid::new([&a[..], &b[..], &c[..]], &values, method)
            .expect("valid")
            .with_extrapolate(Extrapolate::Extend);
        for i in 0..=8 {
            for j in 0..=8 {
                for k in 0..=8 {
                    // From a little outside the grid to a little outside.
                    let p = [
                        -0.2 + 0.425 * f64::from(i),
                        -1.2 + 0.55 * f64::from(j),
                        -0.2 + 0.6125 * f64::from(k),
                    ];
                    let (got, want) = (grid.eval(&p).expect("extends"), f(&p));
                    let tol = 1e-9 * want.abs().max(100.0);
                    assert!(
                        (got - want).abs() <= tol,
                        "{method:?} at {p:?}: {got} for {want}"
                    );
                }
            }
        }
    }
}

#[test]
fn exchanging_the_axes_exchanges_the_coordinates() {
    let a = [0.0, 0.3, 1.0, 1.4, 2.5, 3.0, 3.25];
    let b = [-1.0, -0.5, 0.25, 1.0, 1.5, 2.75];
    let values: Vec<f64> = (0..42_u32).map(|n| f64::from(n * 29 % 17) - 8.0).collect();
    // The transposed grid: the value at (j, i) is the one at (i, j).
    let mut transposed = vec![0.0; 42];
    for i in 0..7 {
        for j in 0..6 {
            transposed[j * 7 + i] = values[i * 6 + j];
        }
    }
    for method in [GridMethod::Linear, GridMethod::Cubic, GridMethod::Quintic] {
        let grid = Grid::new([&a[..], &b[..]], &values, method).expect("valid");
        let other = Grid::new([&b[..], &a[..]], &transposed, method).expect("valid");
        for i in 0..=20 {
            for j in 0..=20 {
                let (x, y) = (
                    3.25 * f64::from(i) / 20.0,
                    -1.0 + 3.75 * f64::from(j) / 20.0,
                );
                let (here, there) = (
                    grid.eval(&[x, y]).expect("inside"),
                    other.eval(&[y, x]).expect("inside"),
                );
                let tol = 1e-11 * here.abs().max(10.0);
                assert!(
                    (here - there).abs() <= tol,
                    "{method:?} at ({x}, {y}): {here} and {there}"
                );
            }
        }
    }
}

#[test]
fn a_grid_of_no_dimensions_is_its_one_value() {
    for method in [
        GridMethod::Nearest,
        GridMethod::Linear,
        GridMethod::Cubic,
        GridMethod::Quintic,
        GridMethod::Pchip,
    ] {
        let point: Grid<f64, 0> = Grid::new([], &[4.5], method).expect("one value");
        assert_eq!(bits(point.eval(&[])), bits(Some(4.5)), "{method:?}");
        assert_eq!(point.shape(), [0_usize; 0]);
        assert_eq!(point.method(), method);
        assert_eq!(
            Grid::<f64, 0>::new([], &[], method).err(),
            Some(E::LengthMismatch {
                expected: 1,
                actual: 0
            })
        );
    }
}

// --- The shape-preserving method ---------------------------------------------------

/// [`PchipInterpolator`] along the last axis, then along each axis before it,
/// with the mode on every line: what `GridMethod::Pchip` is defined as. `None`
/// where a line has no value at its coordinate.
fn folded(axes: &[Vec<f64>], values: &[f64], point: &[f64], mode: Extrapolate) -> Option<f64> {
    let mut values = values.to_vec();
    for (axis, x) in axes.iter().zip(point).rev() {
        values = values
            .chunks_exact(axis.len())
            .map(|line| {
                let along = PchipInterpolator::new(axis, line).ok()?;
                along.with_extrapolate(mode).eval(x)
            })
            .collect::<Option<Vec<f64>>>()?;
    }
    Some(values[0])
}

#[test]
fn pchip_takes_the_axes_in_scipys_order_the_last_one_first() {
    // One grid of values in no order, and the same grid with its two axes
    // exchanged, asked at the same points.
    let file = load::<Json>("grid_pchip").expect("golden file grid_pchip");
    let case = |name: &str| {
        let found = file.cases.iter().find(|c| c["name"] == name);
        found.unwrap_or_else(|| panic!("no case {name}"))
    };
    let (given, exchanged) = (
        case("grid2d/random/pchip"),
        case("grid2d/random/exchanged/pchip"),
    );
    let build = |case: &Json| {
        let axes = rows(&case["axes"]);
        Grid::new(
            [&axes[0][..], &axes[1][..]],
            &nums(&case["values"]),
            GridMethod::Pchip,
        )
        .expect("valid grid")
    };
    let (grid, other) = (build(given), build(exchanged));
    let (axes, values) = (rows(&given["axes"]), nums(&given["values"]));
    let mut apart = (0.0_f64, 0.0_f64);
    let queries = rows(&given["queries"]);
    let scipy = nums(&given["expected"])
        .into_iter()
        .zip(nums(&exchanged["expected"]));
    for (q, (want, want_exchanged)) in queries.iter().zip(scipy) {
        let here = grid.eval(&[q[0], q[1]]).expect("inside");
        let there = other.eval(&[q[1], q[0]]).expect("inside");
        // The last axis first: the bits of the interpolant in one variable
        // applied in that order.
        assert_eq!(
            bits(Some(here)),
            bits(folded(&axes, &values, q, Extrapolate::Extend)),
            "{q:?}"
        );
        assert_close(here, want, DIRECT_F64, "as given");
        assert_close(there, want_exchanged, DIRECT_F64, "exchanged");
        apart = (
            apart.0.max((here - there).abs()),
            apart.1.max((want - want_exchanged).abs()),
        );
    }
    // The other order is another interpolant, here and in SciPy: the order is
    // part of the method. (The linear and the spline methods give one
    // interpolant in either order, see
    // `exchanging_the_axes_exchanges_the_coordinates`.)
    assert!(apart.0 > 1e-3 && apart.1 > 1e-3, "{apart:?}");
}

#[test]
fn pchip_is_flat_where_the_data_is_flat_and_stays_inside_its_cell() {
    // Steps and plateaus along both axes.
    let side: Vec<f64> = (0..8).map(f64::from).collect();
    let level = [0.0, 0.0, 0.0, 1.0, 1.0, 3.0, 3.0, 3.0];
    let at = |i: usize, j: usize| level[i] + 2.0 * level[7 - j];
    let values: Vec<f64> = (0..64).map(|n| at(n / 8, n % 8)).collect();
    let grid = Grid::new([&side[..], &side[..]], &values, GridMethod::Pchip).expect("valid");
    let cubic = Grid::new([&side[..], &side[..]], &values, GridMethod::Cubic).expect("valid");
    let (mut flat, mut below) = (0, f64::INFINITY);
    for a in 0..=70 {
        for b in 0..=70 {
            let p = [0.1 * f64::from(a), 0.1 * f64::from(b)];
            let got = grid.eval(&p).expect("inside");
            // The cell of the point (the last cell at the upper ends).
            let (i, j) = ((a as usize / 10).min(6), (b as usize / 10).min(6));
            let corners = [at(i, j), at(i, j + 1), at(i + 1, j), at(i + 1, j + 1)];
            let lo = corners.iter().copied().fold(f64::INFINITY, f64::min);
            let hi = corners.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            assert!(
                lo - 1e-13 <= got && got <= hi + 1e-13,
                "{p:?}: {got} outside [{lo}, {hi}]"
            );
            if lo.to_bits() == hi.to_bits() {
                // A flat cell is flat, to the bit.
                assert_eq!(got.to_bits(), lo.to_bits(), "{p:?}");
                flat += 1;
            }
            below = below.min(cubic.eval(&p).expect("inside"));
        }
    }
    assert!(flat > 1000, "only {flat} points in flat cells");
    // The property is the method's: the cubic spline through the same values
    // goes below every one of them.
    assert!(below < -0.1, "the cubic spline stays above {below}");
}

#[test]
fn pchip_takes_two_points_an_axis_and_refuses_slopes_that_overflow() {
    let new = |x: &[f64], y: &[f64], values: &[f64]| Grid::new([x, y], values, GridMethod::Pchip);
    // Two points on each axis: a line along each, which is the linear method.
    let (x, y) = ([0.0, 2.0], [1.0, 1.5]);
    let values = [1.0, -2.0, 0.5, 4.0];
    let grid = new(&x, &y, &values).expect("two points an axis are enough");
    assert_eq!(grid.method(), GridMethod::Pchip);
    let linear = Grid::new([&x[..], &y[..]], &values, GridMethod::Linear).expect("valid");
    for p in [[0.0, 1.0], [0.5, 1.1], [1.25, 1.45], [2.0, 1.5]] {
        let (got, want) = (
            grid.eval(&p).expect("inside"),
            linear.eval(&p).expect("inside"),
        );
        assert_close(got, want, DIRECT_F64, "two points an axis");
    }
    assert_eq!(
        new(&x, &[1.0], &[0.0; 2]).err(),
        Some(on_axis(
            1,
            AxisFault::InsufficientPoints {
                required: 2,
                actual: 1
            }
        ))
    );
    // A secant along the last axis that overflows: between the first two
    // points of the second line. The later of the two is named, by its place
    // in the values.
    let (x, y) = ([0.0, 1.0], [0.0, 1.0, 2.0]);
    let huge = 1.7e308;
    assert_eq!(
        new(&x, &y, &[0.0, 0.0, 0.0, huge, -huge, 0.0]).err(),
        Some(E::NonFinite { index: 4 })
    );
    // Secants that are finite, and a slope at the end of the line that is
    // not: three times 8e307 is past the largest number.
    assert_eq!(
        new(&x, &y, &[0.0, 0.0, 0.0, 0.0, 8e307, -8e307]).err(),
        Some(E::NonFinite { index: 3 })
    );
    // A place is named as it was given, also along an axis that descends
    // and is kept the other way round: here the overflow is between the last
    // two values given, and the earlier of them is the later one as kept.
    let one = |x: &[f64], values: &[f64]| Grid::new([x], values, GridMethod::Pchip).err();
    assert_eq!(
        one(&[3.0, 2.0, 1.0, 0.0], &[0.0, 0.0, 1e308, -1e308]),
        Some(E::NonFinite { index: 2 })
    );
    assert_eq!(
        one(&[0.0, 1.0, 2.0, 3.0], &[-1e308, 1e308, 0.0, 0.0]),
        Some(E::NonFinite { index: 1 })
    );
    // The same values along the first axis are not looked at when the grid is
    // built: the lines along the last axis are flat. An evaluation between
    // them overflows, and says so with a value that is not finite.
    let (x, y) = ([0.0, 1.0, 2.0], [0.0, 1.0]);
    let values = [huge, huge, -huge, -huge, huge, huge];
    let grid = new(&x, &y, &values).expect("flat along the last axis");
    assert!(grid.eval(&[0.5, 0.5]).is_some_and(|v| !v.is_finite()));
    let linear = Grid::new([&x[..], &y[..]], &values, GridMethod::Linear).expect("valid");
    assert!(linear.eval(&[0.5, 0.5]).is_some_and(f64::is_finite));
}

#[test]
fn in_one_dimension_pchip_is_the_interpolant_in_one_variable() {
    let x = [0.0, 0.5, 1.5, 2.0, 3.5, 4.0, 5.5, 6.0, 7.0];
    let y = [1.0, -1.0, 0.5, 2.0, 2.0, -2.0, 1.5, 0.5, 3.0];
    for n in 2..=x.len() {
        let grid = Grid::new([&x[..n]], &y[..n], GridMethod::Pchip)
            .expect("valid")
            .with_extrapolate(Extrapolate::Extend);
        let line = PchipInterpolator::new(&x[..n], &y[..n]).expect("valid");
        let mut hints = [Hint::new(); 1];
        for k in 0..=160 {
            let q = -0.5 + 0.05 * f64::from(k);
            assert_eq!(bits(grid.eval(&[q])), bits(line.eval(&q)), "{n} at {q}");
            assert_eq!(
                bits(grid.eval_hinted(&[q], &mut hints)),
                bits(line.eval(&q))
            );
        }
    }
    // f32, two columns: each column is its own interpolant.
    let x32 = [0.0_f32, 1.0, 1.5, 4.0];
    let pairs = [[0.0_f32, 1.0], [2.0, 1.0], [1.0, -5.0], [6.0, 5.0]];
    let both: Grid<f32, 1, [f32; 2]> =
        Grid::new([&x32[..]], &pairs, GridMethod::Pchip).expect("valid");
    let line = PchipInterpolator::new(&x32, &pairs).expect("valid");
    for k in 0..=40_u8 {
        let q = 0.1 * f32::from(k);
        assert_eq!(
            both.eval(&[q]).map(|v| v.map(f32::to_bits)),
            line.eval(&q).map(|v| v.map(f32::to_bits)),
            "at {q}"
        );
    }
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

/// An axis of 6 or 7 points: enough for the quintic.
fn long_axis() -> impl Strategy<Value = Vec<f64>> {
    (-3.0_f64..3.0, prop::collection::vec(0.2_f64..1.5, 5..7)).prop_map(|(start, widths)| {
        let mut x = vec![start];
        for w in widths {
            x.push(x[x.len() - 1] + w);
        }
        x
    })
}

/// [`cases`] for the spline methods, with points up to a tenth of the box
/// outside it.
fn spline_cases() -> impl Strategy<Value = Case> {
    (
        [long_axis(), long_axis(), long_axis()],
        prop::collection::vec(-10.0_f64..10.0, 343),
        prop::collection::vec(prop::array::uniform3(-0.1_f64..1.1), 1..20),
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

    #[test]
    fn a_spline_method_has_the_properties_of_the_grid(
        (axes, values, points) in spline_cases(),
        cubic in any::<bool>(),
        mode in modes(),
        flip in prop::array::uniform3(any::<bool>()),
        scale in -3.0_f64..3.0,
    ) {
        let method = if cubic { GridMethod::Cubic } else { GridMethod::Quintic };
        let grid = Grid::new(slices(&axes), &values, method).expect("valid").with_extrapolate(mode);
        let shape = grid.shape();
        // The values at the grid points. The coefficients of a spline through
        // random values are larger than the values, and round accordingly.
        let size = values.iter().fold(10.0_f64, |m, v| m.max(v.abs()));
        for i in 0..shape[0] {
            for j in 0..shape[1] {
                for k in 0..shape[2] {
                    let wrapped = mode == Extrapolate::Periodic
                        && (i + 1 == shape[0] || j + 1 == shape[1] || k + 1 == shape[2]);
                    if !wrapped {
                        let at = [axes[0][i], axes[1][j], axes[2][k]];
                        let (got, want) = (grid.eval(&at).expect("on the grid"), values[(i * shape[1] + j) * shape[2] + k]);
                        prop_assert!((got - want).abs() <= 1e-9 * size, "{:?} at {:?}: {} for {}", method, at, got, want);
                    }
                }
            }
        }
        // Axes given backward, with the values to match: the same grid.
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
        let backward = Grid::new(slices(&given), &moved, method).expect("valid").with_extrapolate(mode);
        // Two values at every point, and each alone.
        let second: Vec<f64> = values.iter().rev().map(|v| v * scale + 1.0).collect();
        let both: Vec<[f64; 2]> = values.iter().zip(&second).map(|(a, b)| [*a, *b]).collect();
        let pair = Grid::new(slices(&axes), &both, method).expect("valid").with_extrapolate(mode);
        let other = Grid::new(slices(&axes), &second, method).expect("valid").with_extrapolate(mode);
        let clamp = grid.clone().with_extrapolate(Extrapolate::Clamp);
        let plain = grid.clone().with_extrapolate(Extrapolate::Forbid);
        let (lo, hi) = grid.domain();
        let mut hints = [Hint::new(); 3];
        let mut out = vec![0.0; points.len()];
        prop_assert_eq!(grid.eval_batch_into(&points, &mut out, -7.0), Ok(()));
        for (p, slot) in points.iter().zip(&out) {
            let got = grid.eval(p);
            prop_assert_eq!(bits(backward.eval(p)), bits(got), "backward at {:?}", p);
            prop_assert_eq!(bits(grid.eval_hinted(p, &mut hints)), bits(got));
            prop_assert_eq!(bits(Some(*slot)), bits(Some(got.unwrap_or(-7.0))));
            let columns = pair.eval(p);
            prop_assert_eq!(bits(columns.map(|c| c[0])), bits(got));
            prop_assert_eq!(bits(columns.map(|c| c[1])), bits(other.eval(p)));
            let clamped: [f64; 3] = core::array::from_fn(|d| p[d].clamp(lo[d], hi[d]));
            prop_assert_eq!(bits(clamp.eval(p)), bits(plain.eval(&clamped)), "clamp at {:?}", p);
        }
    }

    #[test]
    fn pchip_has_the_properties_of_the_grid(
        (axes, values, points) in cases(),
        mode in modes(),
        flip in prop::array::uniform3(any::<bool>()),
        scale in -3.0_f64..3.0,
    ) {
        let method = GridMethod::Pchip;
        let grid = Grid::new(slices(&axes), &values, method).expect("valid").with_extrapolate(mode);
        let shape = grid.shape();
        // The values at the grid points. At the upper end of an axis the
        // point is the far end of the last interval, where the cubic meets
        // the value to the rounding only.
        for i in 0..shape[0] {
            for j in 0..shape[1] {
                for k in 0..shape[2] {
                    let wrapped = mode == Extrapolate::Periodic
                        && (i + 1 == shape[0] || j + 1 == shape[1] || k + 1 == shape[2]);
                    if !wrapped {
                        let at = [axes[0][i], axes[1][j], axes[2][k]];
                        let (got, want) = (grid.eval(&at).expect("on the grid"), values[(i * shape[1] + j) * shape[2] + k]);
                        prop_assert!((got - want).abs() <= 1e-11, "at {:?}: {} for {}", at, got, want);
                    }
                }
            }
        }
        // Axes given backward, with the values to match: the same grid.
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
        let backward = Grid::new(slices(&given), &moved, method).expect("valid").with_extrapolate(mode);
        // Two values at every point, and each alone.
        let second: Vec<f64> = values.iter().rev().map(|v| v * scale + 1.0).collect();
        let both: Vec<[f64; 2]> = values.iter().zip(&second).map(|(a, b)| [*a, *b]).collect();
        let pair = Grid::new(slices(&axes), &both, method).expect("valid").with_extrapolate(mode);
        let other = Grid::new(slices(&axes), &second, method).expect("valid").with_extrapolate(mode);
        let clamp = grid.clone().with_extrapolate(Extrapolate::Clamp);
        let plain = grid.clone().with_extrapolate(Extrapolate::Forbid);
        let (lo, hi) = grid.domain();
        let mut hints = [Hint::new(); 3];
        let mut out = vec![0.0; points.len()];
        prop_assert_eq!(grid.eval_batch_into(&points, &mut out, -7.0), Ok(()));
        for (p, slot) in points.iter().zip(&out) {
            let got = grid.eval(p);
            // The interpolant in one variable, applied from the last axis to
            // the first, with the mode on every line: to the bit.
            prop_assert_eq!(bits(got), bits(folded(&axes, &values, p, mode)), "folded at {:?}", p);
            prop_assert_eq!(bits(backward.eval(p)), bits(got), "backward at {:?}", p);
            prop_assert_eq!(bits(grid.eval_hinted(p, &mut hints)), bits(got));
            prop_assert_eq!(bits(Some(*slot)), bits(Some(got.unwrap_or(-7.0))));
            let columns = pair.eval(p);
            prop_assert_eq!(bits(columns.map(|c| c[0])), bits(got));
            prop_assert_eq!(bits(columns.map(|c| c[1])), bits(other.eval(p)));
            let clamped: [f64; 3] = core::array::from_fn(|d| p[d].clamp(lo[d], hi[d]));
            prop_assert_eq!(bits(clamp.eval(p)), bits(plain.eval(&clamped)), "clamp at {:?}", p);
            // Inside the grid, no value past the corners of the cell.
            if let Some(inside) = plain.eval(p) {
                let cell: [usize; 3] = core::array::from_fn(|d| {
                    axes[d].partition_point(|a| *a <= p[d]).clamp(1, shape[d] - 1) - 1
                });
                let (mut least, mut most) = (f64::INFINITY, f64::NEG_INFINITY);
                for corner in 0..8_usize {
                    let at: [usize; 3] = core::array::from_fn(|d| cell[d] + (corner >> d & 1));
                    let value = values[(at[0] * shape[1] + at[1]) * shape[2] + at[2]];
                    least = least.min(value);
                    most = most.max(value);
                }
                prop_assert!(least - 1e-11 <= inside && inside <= most + 1e-11, "at {:?}: {} outside [{}, {}]", p, inside, least, most);
            }
        }
    }
}
