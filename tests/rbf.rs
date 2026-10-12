//! `RBFInterpolator` solved over all its points: SciPy's values, the input
//! checks and their order, the systems without a solution, and the properties.

#![cfg(feature = "scattered")]

mod common;

use common::{load, reals, Real};
use interpx::{
    InterpolateError as E, Interpolator, ParamReason, RBFInterpolator as Rbf, RbfDegree, RbfKernel,
    RbfOptions, SolverError,
};
use proptest::prelude::*;
use serde_json::Value as Json;

/// How far a value may be from SciPy's, as a share of the largest value of
/// the case.
///
/// The values come out of a dense system, which SciPy solves another way (a
/// symmetric indefinite factorization, where this crate has an LU): the two
/// differ by the rounding of the solve, grown by the condition of the system.
/// Measured over the cases of `rbf.json` that are solved globally, the
/// largest difference is 9.0e-13 (the Gaussian on 64 points; 5.4e-13 for the
/// quintic on 33). A symmetric factorization of `faer` is as far from SciPy
/// (5.5e-13), so the difference is the system's and not the method's.
const SOLVED: f64 = 1e-11;

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

/// Every kernel, with `epsilon` where it takes one.
fn kernels(epsilon: f64) -> [(&'static str, RbfKernel); 8] {
    [
        ("linear", RbfKernel::Linear),
        ("thin_plate_spline", RbfKernel::ThinPlateSpline),
        ("cubic", RbfKernel::Cubic),
        ("quintic", RbfKernel::Quintic),
        ("multiquadric", RbfKernel::Multiquadric { epsilon }),
        (
            "inverse_multiquadric",
            RbfKernel::InverseMultiquadric { epsilon },
        ),
        ("inverse_quadratic", RbfKernel::InverseQuadratic { epsilon }),
        ("gaussian", RbfKernel::Gaussian { epsilon }),
    ]
}

/// The kernel and the options a golden case was made with.
// Marked as test code so that `expect` is allowed here as in a `#[test]` function
// (clippy.toml allows it in tests only). The same for the helpers below.
#[cfg(test)]
fn asked(case: &Json) -> (RbfKernel, RbfOptions) {
    let options = &case["options"];
    // A kernel that takes no `epsilon` has none in the case.
    let epsilon = options["epsilon"].as_f64().unwrap_or(1.0);
    let (_, kernel) = kernels(epsilon)
        .into_iter()
        .find(|(named, _)| options["kernel"] == *named)
        .expect("a kernel");
    let degree = match options["degree"].as_i64() {
        None => RbfDegree::Default,
        // SciPy's way to write no polynomial terms.
        Some(-1) => RbfDegree::None,
        Some(degree) => RbfDegree::Degree(usize::try_from(degree).expect("a degree")),
    };
    let smoothing = options["smoothing"].as_f64().expect("smoothing");
    (
        kernel,
        RbfOptions::new()
            .with_smoothing(smoothing)
            .with_degree(degree),
    )
}

/// Rebuilds one golden case in `D` dimensions and compares it with SciPy.
/// Returns the largest difference, as a share of the largest value.
#[cfg(test)]
fn check<const D: usize>(case: &Json) -> f64 {
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
    let (kernel, options) = asked(case);
    let built = Rbf::with_options(&points, &values, kernel, &options);
    if case.get("expected_error").is_some() {
        // What SciPy calls a singular matrix.
        assert_eq!(
            built.err(),
            Some(E::SolverFailed(SolverError::Singular)),
            "{name}"
        );
        return 0.0;
    }
    let rbf = built.expect("a system with a solution");
    assert_eq!(rbf.len(), points.len(), "{name}");
    assert_eq!(rbf.kernel(), kernel, "{name}");
    assert_eq!(
        rbf.smoothing().to_bits(),
        options.smoothing().to_bits(),
        "{name}"
    );
    let expected = nums(&case["expected"]);
    let largest = expected.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    let mut batch = vec![0.0; queries.len()];
    assert_eq!(rbf.eval_batch_into(&queries, &mut batch, f64::NAN), Ok(()));
    let mut worst = 0.0_f64;
    for (k, (query, want)) in queries.iter().zip(&expected).enumerate() {
        let got = rbf.eval(query).expect("a value");
        let apart = (got - want).abs() / largest;
        assert!(
            apart <= SOLVED,
            "{name} query {k}: {got} for {want} ({apart:e} of {largest})"
        );
        assert_eq!(batch[k].to_bits(), got.to_bits(), "{name} query {k}");
        worst = worst.max(apart);
    }

    // Two values at every point: each column is its own interpolant, and is
    // solved with the same steps as when it is alone.
    if let Some(columns) = case.get("values_columns") {
        let pairs: Vec<[f64; 2]> = rows(columns).iter().map(|r| [r[0], r[1]]).collect();
        let pair = Rbf::with_options(&points, &pairs, kernel, &options).expect("a solution");
        let wanted = rows(&case["expected_columns"]);
        let size = |column: usize| wanted.iter().fold(0.0_f64, |m, w| m.max(w[column].abs()));
        for (query, want) in queries.iter().zip(&wanted) {
            let got = pair.eval(query).expect("a value");
            for column in 0..2 {
                let apart = (got[column] - want[column]).abs() / size(column);
                assert!(apart <= SOLVED, "{name} column {column}: {apart:e}");
            }
            assert_eq!(bits(Some(got[0])), bits(rbf.eval(query)), "{name}");
        }
    }
    worst
}

// --- Golden data -------------------------------------------------------------

#[test]
fn scipy_rbf_interpolants_solved_globally_evaluate_to_the_golden_values() {
    let file = load::<Json>("rbf").expect("golden file rbf");
    let mut by_kernel = std::collections::BTreeMap::new();
    let (mut local, mut refused, mut worst) = (0, 0, (0.0_f64, String::new()));
    for case in &file.cases {
        // The `neighbors` mode is another row (M3-4).
        if !case["options"]["neighbors"].is_null() {
            local += 1;
            continue;
        }
        let apart = match rows(&case["points"])[0].len() {
            1 => check::<1>(case),
            2 => check::<2>(case),
            3 => check::<3>(case),
            4 => check::<4>(case),
            other => panic!("points in {other} dimensions"),
        };
        refused += usize::from(case.get("expected_error").is_some());
        let kernel = case["options"]["kernel"]
            .as_str()
            .expect("kernel")
            .to_owned();
        *by_kernel.entry(kernel).or_insert(0) += 1;
        if apart > worst.0 {
            worst = (apart, case["name"].as_str().expect("name").to_owned());
        }
    }
    let counts: Vec<(String, i32)> = by_kernel.into_iter().collect();
    let want = [
        ("cubic", 11),
        ("gaussian", 10),
        ("inverse_multiquadric", 9),
        ("inverse_quadratic", 9),
        ("linear", 10),
        ("multiquadric", 10),
        ("quintic", 10),
        ("thin_plate_spline", 14),
    ]
    .map(|(kernel, n)| (kernel.to_owned(), n));
    assert_eq!(counts, want);
    assert_eq!((local, refused), (24, 1));
    // The tolerance has room, and not more than it needs: the worst case is
    // within a tenth of it, and over a thousandth.
    assert!(
        worst.0 < SOLVED / 10.0 && worst.0 > SOLVED / 1000.0,
        "the worst case is {} at {:e}",
        worst.1,
        worst.0
    );
}

// --- Input checks ---------------------------------------------------------------

fn invalid(name: &'static str, reason: ParamReason) -> E {
    E::InvalidParameter { name, reason }
}

#[test]
fn the_constructor_reports_each_fault() {
    let points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0], [0.4, 0.7]];
    let values = [1.0, 2.0, 0.0, 5.0, 2.0];
    let tps = RbfKernel::ThinPlateSpline;
    let err = |points: &[[f64; 2]], values: &[f64], kernel, options: RbfOptions| {
        Rbf::with_options(points, values, kernel, &options).err()
    };
    let plain = RbfOptions::new();
    assert_eq!(err(&points, &values, tps, plain), None);
    assert_eq!(
        err(&points, &values[..4], tps, plain),
        Some(E::LengthMismatch {
            expected: 5,
            actual: 4
        })
    );
    assert_eq!(
        err(&[], &[], tps, plain),
        Some(E::InsufficientPoints {
            required: 1,
            actual: 0
        })
    );
    let mut bad = points;
    bad[3][1] = f64::NAN;
    assert_eq!(
        err(&bad, &values, tps, plain),
        Some(E::NonFinite { index: 3 })
    );
    bad[3][1] = f64::NEG_INFINITY;
    assert_eq!(
        err(&bad, &values, tps, plain),
        Some(E::NonFinite { index: 3 })
    );
    assert_eq!(
        err(
            &points,
            &[1.0, f64::INFINITY, 0.0, 5.0, f64::NAN],
            tps,
            plain
        ),
        Some(E::NonFinite { index: 1 })
    );
    // `epsilon`: finite and above zero, for each kernel that has one.
    for (name, kernel) in kernels(f64::NAN).into_iter().skip(4) {
        let want = Some(invalid("epsilon", ParamReason::NonFinite));
        assert_eq!(err(&points, &values, kernel, plain), want, "{name}");
    }
    for (name, kernel) in kernels(f64::INFINITY).into_iter().skip(4) {
        let want = Some(invalid("epsilon", ParamReason::NonFinite));
        assert_eq!(err(&points, &values, kernel, plain), want, "{name}");
    }
    for epsilon in [0.0, -0.0, -2.0] {
        for (name, kernel) in kernels(epsilon).into_iter().skip(4) {
            let want = Some(invalid("epsilon", ParamReason::NotPositive));
            assert_eq!(
                err(&points, &values, kernel, plain),
                want,
                "{name} {epsilon}"
            );
        }
        // The other four have none to be wrong.
        for (name, kernel) in kernels(epsilon).into_iter().take(3) {
            assert_eq!(err(&points, &values, kernel, plain), None, "{name}");
        }
    }
    // Smoothing: finite and not under zero.
    for smoothing in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            err(&points, &values, tps, plain.with_smoothing(smoothing)),
            Some(invalid("smoothing", ParamReason::NonFinite))
        );
    }
    assert_eq!(
        err(&points, &values, tps, plain.with_smoothing(-1e-300)),
        Some(invalid("smoothing", ParamReason::Negative))
    );
    assert_eq!(err(&points, &values, tps, plain.with_smoothing(0.0)), None);
    // The degree: at the smallest the kernel takes, and not under it.
    let degree = |k| plain.with_degree(RbfDegree::Degree(k));
    assert_eq!(
        err(&points, &values, tps, degree(0)),
        Some(invalid("degree", ParamReason::TooSmall { min: 1 }))
    );
    assert_eq!(err(&points, &values, tps, degree(1)), None);
    assert_eq!(
        err(&points, &values, tps, degree(usize::MAX)),
        // A power is an `i32`.
        Some(invalid(
            "degree",
            ParamReason::TooLarge {
                max: i32::MAX as usize
            }
        ))
    );
    // Degree 2 has 6 monomials in two variables, and there are 5 points.
    assert_eq!(
        err(&points, &values, tps, degree(2)),
        Some(E::InsufficientPoints {
            required: 6,
            actual: 5
        })
    );
    // The same point twice: the earlier one first; the zeros are one number.
    let mut twice = points;
    twice[4] = twice[1];
    assert_eq!(
        err(&twice, &values, tps, plain),
        Some(E::DuplicatePoints {
            first: 1,
            second: 4
        })
    );
    twice[4] = [-0.0, 0.0];
    assert_eq!(
        err(&twice, &values, tps, plain),
        Some(E::DuplicatePoints {
            first: 0,
            second: 4
        })
    );
    // With smoothing the system could be solved; the rule stays.
    assert_eq!(
        err(&twice, &values, tps, plain.with_smoothing(0.5)),
        Some(E::DuplicatePoints {
            first: 0,
            second: 4
        })
    );
}

#[test]
fn the_checks_come_in_the_order_of_the_design() {
    let nan = f64::NAN;
    let gaussian = |epsilon| RbfKernel::Gaussian { epsilon };
    let err = |points: &[[f64; 1]], values: &[f64], kernel, options: RbfOptions| {
        Rbf::with_options(points, values, kernel, &options).err()
    };
    let plain = RbfOptions::new();
    let wrong = plain.with_smoothing(-1.0).with_degree(RbfDegree::Degree(9));
    // Everything wrong at once: the length.
    assert_eq!(
        err(&[], &[nan], gaussian(nan), wrong),
        Some(E::LengthMismatch {
            expected: 0,
            actual: 1
        })
    );
    // Then the count, a coordinate, a value.
    assert_eq!(
        err(&[], &[], gaussian(nan), wrong),
        Some(E::InsufficientPoints {
            required: 1,
            actual: 0
        })
    );
    assert_eq!(
        err(&[[0.0], [nan]], &[nan, 0.0], gaussian(nan), wrong),
        Some(E::NonFinite { index: 1 })
    );
    assert_eq!(
        err(&[[0.0], [0.0]], &[nan, 0.0], gaussian(nan), wrong),
        Some(E::NonFinite { index: 0 })
    );
    // Then `epsilon`, the smoothing, the degree.
    assert_eq!(
        err(&[[0.0], [0.0]], &[1.0, 0.0], gaussian(-1.0), wrong),
        Some(invalid("epsilon", ParamReason::NotPositive))
    );
    assert_eq!(
        err(&[[0.0], [0.0]], &[1.0, 0.0], gaussian(1.0), wrong),
        Some(invalid("smoothing", ParamReason::Negative))
    );
    let low = plain.with_degree(RbfDegree::Degree(1));
    assert_eq!(
        err(&[[0.0], [0.0]], &[1.0, 0.0], RbfKernel::Quintic, low),
        Some(invalid("degree", ParamReason::TooSmall { min: 2 }))
    );
    // Then too few points for the monomials, before the same point twice.
    let high = plain.with_degree(RbfDegree::Degree(9));
    assert_eq!(
        err(&[[0.0], [0.0]], &[1.0, 0.0], gaussian(1.0), high),
        Some(E::InsufficientPoints {
            required: 10,
            actual: 2
        })
    );
    assert_eq!(
        err(&[[0.0], [0.0]], &[1.0, 0.0], gaussian(1.0), plain),
        Some(E::DuplicatePoints {
            first: 0,
            second: 1
        })
    );
}

// --- Systems without a solution ------------------------------------------------------

#[test]
fn a_system_without_one_solution_is_an_error_and_not_a_value() {
    let singular = Some(E::SolverFailed(SolverError::Singular));
    // Points along an axis, with a polynomial of degree 1 in two variables:
    // the second variable is the same at every point, and its coefficient is
    // not determined.
    let along: Vec<[f64; 2]> = (0..8).map(|k| [f64::from(k), 3.0]).collect();
    let values: Vec<f64> = (0..8).map(|k| f64::from(k * k)).collect();
    for kernel in [RbfKernel::ThinPlateSpline, RbfKernel::Cubic] {
        assert_eq!(
            Rbf::new(&along, &values, kernel).err(),
            singular,
            "{kernel:?}"
        );
    }
    // The linear kernel takes a constant only, which the points determine.
    assert!(Rbf::new(&along, &values, RbfKernel::Linear).is_ok());
    // In one dimension the same points are no special case.
    let line: Vec<[f64; 1]> = along.iter().map(|p| [p[0]]).collect();
    assert!(Rbf::new(&line, &values, RbfKernel::ThinPlateSpline).is_ok());
    // One point and a kernel that is 0 at the distance 0, without a
    // polynomial: the system is the equation `0 w = value`.
    let none = RbfOptions::new().with_degree(RbfDegree::None);
    for kernel in [
        RbfKernel::Linear,
        RbfKernel::ThinPlateSpline,
        RbfKernel::Cubic,
    ] {
        let built = Rbf::with_options(&[[1.0, 2.0]], &[5.0], kernel, &none);
        assert_eq!(built.err(), singular, "{kernel:?}");
    }
    // Smoothing puts a number on that diagonal.
    let smooth = none.with_smoothing(2.0);
    let built = Rbf::with_options(&[[1.0, 2.0]], &[5.0], RbfKernel::Linear, &smooth);
    assert_eq!(
        bits(built.expect("solvable").eval(&[1.0, 2.0])),
        bits(Some(0.0))
    );
    // Values near the largest number: the solution does not fit.
    let far = [[0.0], [1.0], [2.0], [3.0]];
    let huge = [1.7e308, -1.7e308, 1.7e308, -1.7e308];
    let built = Rbf::new(&far, &huge, RbfKernel::Gaussian { epsilon: 0.1 });
    assert_eq!(built.err(), singular);
}

// --- Evaluation ------------------------------------------------------------------

#[test]
fn a_query_or_a_value_that_is_not_finite_is_no_value() {
    let points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0], [0.4, 0.7]];
    let values = [1.0, 2.0, 0.0, 5.0, 2.0];
    for (name, kernel) in kernels(1.5) {
        // The quintic takes degree 2, which needs a sixth point.
        let options = RbfOptions::new().with_degree(RbfDegree::None);
        let rbf = Rbf::with_options(&points, &values, kernel, &options).expect("valid");
        for query in [
            [f64::NAN, 0.0],
            [0.0, f64::NAN],
            [f64::INFINITY, 0.0],
            [0.0, f64::NEG_INFINITY],
        ] {
            assert_eq!(rbf.eval(&query), None, "{name} {query:?}");
            assert_eq!(bits(Some(rbf.eval_or(&query, -7.0))), bits(Some(-7.0)));
        }
        // There is no domain: a value far from the data.
        assert!(rbf.eval(&[50.0, -80.0]).is_some(), "{name}");
        assert_eq!(rbf.degree(), None, "{name}");
    }
    // So far that the distance does not fit.
    let cubic = Rbf::new(&points, &values, RbfKernel::Cubic).expect("valid");
    assert_eq!(cubic.eval(&[1e200, 0.0]), None);
    // A batch fills those, and refuses a buffer of another length.
    let queries = [[0.5, 0.5], [f64::NAN, 0.0], [1e200, 0.0], [0.0, 0.0]];
    let mut out = [0.0; 4];
    assert_eq!(cubic.eval_batch_into(&queries, &mut out, -7.0), Ok(()));
    for (query, got) in queries.iter().zip(out) {
        assert_eq!(got.to_bits(), cubic.eval_or(query, -7.0).to_bits());
    }
    let mut short = [1.5; 3];
    assert_eq!(
        cubic.eval_batch_into(&queries, &mut short, -7.0),
        Err(E::LengthMismatch {
            expected: 4,
            actual: 3
        })
    );
    assert!(short.iter().all(|v| v.to_bits() == 1.5_f64.to_bits()));
}

#[test]
fn the_degree_asked_for_is_the_degree_used() {
    let points: Vec<[f64; 2]> = (0..12_u32)
        .map(|k| [f64::from(k * 5 % 12) * 0.3, f64::from(k * k % 7) * 0.4])
        .collect();
    let values: Vec<f64> = (0..12).map(|k| f64::from(k % 5)).collect();
    // The options keep what they are given, and start from the defaults.
    let options = RbfOptions::new();
    assert_eq!(options, RbfOptions::default());
    assert_eq!(
        (options.smoothing().to_bits(), options.degree()),
        (0.0_f64.to_bits(), RbfDegree::Default)
    );
    let changed = options
        .with_smoothing(0.25)
        .with_degree(RbfDegree::Degree(3));
    assert_eq!(
        (changed.smoothing().to_bits(), changed.degree()),
        (0.25_f64.to_bits(), RbfDegree::Degree(3))
    );
    assert_eq!(
        changed.with_degree(RbfDegree::None).degree(),
        RbfDegree::None
    );
    let defaults = [0, 1, 1, 2, 0, 0, 0, 0];
    for ((name, kernel), default) in kernels(1.2).into_iter().zip(defaults) {
        let with = |degree| {
            let options = RbfOptions::new().with_degree(degree);
            Rbf::with_options(&points, &values, kernel, &options).map(|rbf| rbf.degree())
        };
        assert_eq!(with(RbfDegree::Default), Ok(Some(default)), "{name}");
        assert_eq!(
            Rbf::new(&points, &values, kernel).map(|rbf| rbf.degree()),
            Ok(Some(default)),
            "{name}"
        );
        assert_eq!(with(RbfDegree::Degree(2)), Ok(Some(2)), "{name}");
        assert_eq!(with(RbfDegree::None), Ok(None), "{name}");
    }
}

#[test]
fn no_dimensions_columns_and_thread_safety() {
    // Without coordinates there is one point at most, and its value
    // everywhere.
    let point: Rbf<0> = Rbf::new(&[[]], &[4.5], RbfKernel::Linear).expect("valid");
    assert_eq!(bits(point.eval(&[])), bits(Some(4.5)));
    assert_eq!(
        Rbf::<0>::new(&[[], []], &[1.0, 2.0], RbfKernel::Linear).err(),
        Some(E::DuplicatePoints {
            first: 0,
            second: 1
        })
    );
    // Arrays for values: a NaN in any column is found.
    let nan = Rbf::new(
        &[[0.0], [1.0]],
        &[[1.0, 1.0], [2.0, f64::NAN]],
        RbfKernel::Linear,
    );
    assert_eq!(nan.err(), Some(E::NonFinite { index: 1 }));

    fn assert_send_sync<T: Send + Sync>() {}
    assert_send_sync::<Rbf<3>>();
    assert_send_sync::<Rbf<2, [f64; 4]>>();
    let rbf = Rbf::new(&[[0.0], [1.0], [3.0]], &[1.0, 2.0, 0.0], RbfKernel::Linear).expect("valid");
    let other = rbf.clone();
    let from_thread = std::thread::spawn(move || other.eval(&[2.0]))
        .join()
        .expect("joins");
    assert_eq!(bits(from_thread), bits(rbf.eval(&[2.0])));
}

// --- Properties ----------------------------------------------------------------

/// Points in the unit square that are no closer than a twentieth to one
/// another (close points make a system that is hard to solve, which is not
/// what these properties are about), and values for them.
fn data() -> impl Strategy<Value = (Vec<[f64; 2]>, Vec<f64>)> {
    prop::collection::vec(prop::array::uniform2(0.0_f64..1.0), 12..40).prop_flat_map(|candidates| {
        let mut points: Vec<[f64; 2]> = Vec::new();
        for c in candidates {
            let near = |p: &[f64; 2]| (p[0] - c[0]).hypot(p[1] - c[1]) < 0.05;
            if !points.iter().any(near) {
                points.push(c);
            }
        }
        let count = points.len();
        (Just(points), prop::collection::vec(-1.0_f64..1.0, count))
    })
}

/// Any kernel. `epsilon` is large enough for the bumps to be narrower than
/// the square: with wide bumps the points are told apart by little, and the
/// system is solved to few digits (which is the kernel's doing).
fn kernel() -> impl Strategy<Value = RbfKernel> {
    (0..8_usize, 4.0_f64..10.0).prop_map(|(which, epsilon)| kernels(epsilon)[which].1)
}

proptest! {
    #[test]
    fn without_smoothing_the_interpolant_passes_through_the_data(
        (points, values) in data(),
        kernel in kernel(),
    ) {
        prop_assume!(points.len() >= 6);
        let rbf = Rbf::new(&points, &values, kernel).expect("valid");
        for (point, value) in points.iter().zip(&values) {
            let got = rbf.eval(point).expect("a value");
            prop_assert!((got - value).abs() <= 1e-6, "{:?} at {:?}: {} for {}", kernel, point, got, value);
        }
    }

    #[test]
    fn a_polynomial_of_the_degree_is_reproduced_everywhere(
        (points, _) in data(),
        kernel in kernel(),
        c in prop::array::uniform6(-2.0_f64..2.0),
        queries in prop::collection::vec(prop::array::uniform2(-1.0_f64..2.0), 1..20),
    ) {
        prop_assume!(points.len() >= 8);
        // Degree 2 in two variables: six monomials. Every kernel takes it.
        let f = |p: &[f64; 2]| c[0] + c[1] * p[0] + c[2] * p[1] + c[3] * p[0] * p[0] + c[4] * p[0] * p[1] + c[5] * p[1] * p[1];
        let values: Vec<f64> = points.iter().map(f).collect();
        let options = RbfOptions::new().with_degree(RbfDegree::Degree(2));
        let rbf = Rbf::with_options(&points, &values, kernel, &options).expect("valid");
        // The polynomial alone solves the system: the weights are zero, and
        // the interpolant is the polynomial also away from the data.
        for query in &queries {
            let (got, want) = (rbf.eval(query).expect("a value"), f(query));
            prop_assert!((got - want).abs() <= 1e-5 * want.abs().max(1.0), "{:?} at {:?}: {} for {}", kernel, query, got, want);
        }
    }

    #[test]
    fn a_kernel_without_epsilon_does_not_depend_on_the_scale(
        (points, values) in data(),
        which in 0..4_usize,
        stretch in 0.01_f64..100.0,
        queries in prop::collection::vec(prop::array::uniform2(0.0_f64..1.0), 1..10),
    ) {
        prop_assume!(points.len() >= 8);
        let kernel = kernels(1.0)[which].1;
        let rbf = Rbf::new(&points, &values, kernel).expect("valid");
        let stretched: Vec<[f64; 2]> = points.iter().map(|p| p.map(|c| c * stretch)).collect();
        let other = Rbf::new(&stretched, &values, kernel).expect("valid");
        for query in &queries {
            let here = rbf.eval(query).expect("a value");
            let there = other.eval(&query.map(|c| c * stretch)).expect("a value");
            prop_assert!((here - there).abs() <= 1e-6 * here.abs().max(1.0), "{:?} stretched by {}: {} and {}", kernel, stretch, here, there);
        }
    }

    #[test]
    fn columns_and_batches_agree_with_eval(
        (points, values) in data(),
        kernel in kernel(),
        smoothing in prop::option::of(0.0_f64..2.0),
        scale in -3.0_f64..3.0,
        queries in prop::collection::vec(prop::array::uniform2(-0.5_f64..1.5), 1..20),
    ) {
        prop_assume!(points.len() >= 6);
        let options = RbfOptions::new().with_smoothing(smoothing.unwrap_or(0.0));
        let second: Vec<f64> = values.iter().rev().map(|v| v * scale + 1.0).collect();
        let both: Vec<[f64; 2]> = values.iter().zip(&second).map(|(a, b)| [*a, *b]).collect();
        let rbf = Rbf::with_options(&points, &values, kernel, &options).expect("valid");
        let other = Rbf::with_options(&points, &second, kernel, &options).expect("valid");
        let pair = Rbf::with_options(&points, &both, kernel, &options).expect("valid");
        let mut out = vec![0.0; queries.len()];
        prop_assert_eq!(rbf.eval_batch_into(&queries, &mut out, -7.0), Ok(()));
        for (query, slot) in queries.iter().zip(&out) {
            let got = rbf.eval(query);
            prop_assert_eq!(bits(Some(*slot)), bits(Some(got.unwrap_or(-7.0))));
            // Each column goes through the same steps as when it is alone.
            let columns = pair.eval(query);
            prop_assert_eq!(bits(columns.map(|c| c[0])), bits(got));
            prop_assert_eq!(bits(columns.map(|c| c[1])), bits(other.eval(query)));
        }
    }

    #[test]
    fn more_smoothing_leaves_the_data_farther_and_the_weights_smaller(
        (points, values) in data(),
        which in 0..4_usize,
    ) {
        prop_assume!(points.len() >= 8);
        let kernel = kernels(1.0)[which].1;
        // How far the interpolant is from the data, summed over the points.
        let miss = |smoothing: f64| {
            let options = RbfOptions::new().with_smoothing(smoothing);
            let rbf = Rbf::with_options(&points, &values, kernel, &options).expect("valid");
            points.iter().zip(&values).map(|(p, v)| (rbf.eval(p).expect("a value") - v).powi(2)).sum::<f64>()
        };
        let (none, some, much) = (miss(0.0), miss(0.05), miss(5.0));
        prop_assert!(none <= 1e-10, "{:?}: {} without smoothing", kernel, none);
        prop_assert!(some > none && much > some, "{:?}: {} then {} then {}", kernel, none, some, much);
    }
}
