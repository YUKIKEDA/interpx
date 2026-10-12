//! Allocation counts of the evaluation APIs (docs/design.ja.md section 8.1).
//!
//! `eval`, `eval_or`, and `eval_batch_into` allocate nothing, nor do
//! `eval_derivative` and `integrate`. `eval_batch` returns a `Vec`, so it
//! allocates, as do `derivative` and `antiderivative`, which build an object.
//!
//! The counter is global to the process, so nothing else may allocate while a
//! region is measured. This file is therefore built without the test harness
//! (`harness = false` in `Cargo.toml`): the harness runs a test on a thread of
//! its own and keeps its records on the main thread, and an allocation there is
//! counted too. Here `main` is the only thread. `cargo test` runs it like any
//! test; a failed assertion ends it with a non-zero status.
//!
//! A row that adds an interpolator adds it to `main`.

use interpx::{
    Akima1D, AkimaMethod, BSpline, BSpline1D, BSplineCurve, BoundaryCondition, CubicHermiteSpline,
    CubicSpline, Extrapolate, GridMethod, Hint, InterpolateError, Interpolator, Linear1D,
    Nearest1D, NearestKind, PPoly, Parametrization, PchipInterpolator, RegularGridInterpolator,
};
use stats_alloc::{Region, StatsAlloc, INSTRUMENTED_SYSTEM};
use std::alloc::System;

#[global_allocator]
static GLOBAL: &StatsAlloc<System> = &INSTRUMENTED_SYSTEM;

/// The number of allocations and reallocations `f` makes.
fn allocations(f: impl FnOnce()) -> usize {
    let region = Region::new(GLOBAL);
    f();
    let stats = region.change();
    stats.allocations + stats.reallocations
}

/// Checks the allocation contract of one interpolator on the given points.
fn check<I>(name: &str, interp: &I, points: &[I::Point], fill: I::Output)
where
    I: Interpolator,
    I::Output: Clone,
{
    assert!(!points.is_empty(), "{name}: no points to evaluate");
    let mut out = vec![fill.clone(); points.len()];

    let n = allocations(|| {
        for p in points {
            std::hint::black_box(interp.eval(p));
        }
    });
    assert_eq!(n, 0, "{name}: eval allocates");

    let n = allocations(|| {
        for p in points {
            std::hint::black_box(interp.eval_or(p, fill.clone()));
        }
    });
    assert_eq!(n, 0, "{name}: eval_or allocates");

    let n = allocations(|| {
        std::hint::black_box(interp.eval_batch_into(points, &mut out, fill.clone())).ok();
    });
    assert_eq!(n, 0, "{name}: eval_batch_into allocates");

    // The same call with a wrong length returns early; it must not allocate either.
    let n = allocations(|| {
        std::hint::black_box(interp.eval_batch_into(&points[1..], &mut out, fill.clone())).ok();
    });
    assert_eq!(
        n, 0,
        "{name}: eval_batch_into allocates on a length mismatch"
    );

    // `eval_batch` is the one that allocates: this shows the counter is live.
    let n = allocations(|| {
        std::hint::black_box(interp.eval_batch(points));
    });
    assert!(
        n > 0,
        "{name}: eval_batch did not allocate, is the counter on?"
    );
}

/// `2x` on `[0, 1]`: uses the provided methods of the trait as they are.
struct Double;

impl Interpolator for Double {
    type Point = f64;
    type Output = f64;

    fn eval(&self, x: &f64) -> Option<f64> {
        (0.0..=1.0).contains(x).then(|| 2.0 * x)
    }
}

/// Three outputs per point, to cover an array value.
struct Triple;

impl Interpolator for Triple {
    type Point = [f64; 2];
    type Output = [f64; 3];

    fn eval(&self, p: &[f64; 2]) -> Option<[f64; 3]> {
        p.iter()
            .all(|c| c.is_finite())
            .then(|| [p[0], p[1], p[0] + p[1]])
    }
}

/// The interpolators of scattered data.
#[cfg(feature = "scattered")]
fn scattered() -> Result<(), InterpolateError> {
    use interpx::NearestNDInterpolator;

    // NearestNDInterpolator: enough points for a tree with several leaves,
    // some of them given twice; plain and rescaled, with and without a bound;
    // scalar and array values. The queries are near, far, past the bound, too
    // far to measure, and not finite.
    let points: Vec<[f64; 3]> = (0..200_u32)
        .map(|n| {
            [
                f64::from(n * 37 % 100),
                f64::from(n * 61 % 50),
                f64::from(n % 7),
            ]
        })
        .collect();
    let values: Vec<f64> = (0..200).map(f64::from).collect();
    let pairs: Vec<[f64; 2]> = values.iter().map(|v| [*v, -*v]).collect();
    let queries = [
        [10.2, 20.1, 3.3],
        [99.0, 0.0, 6.0],
        [500.0, 500.0, 500.0],
        [1e200, 0.0, 0.0],
        [f64::NAN, 1.0, 1.0],
        [0.0, f64::INFINITY, 1.0],
    ];
    type Build<V> = fn(&[[f64; 3]], &[V]) -> Result<NearestNDInterpolator<3, V>, InterpolateError>;
    let scalars: [(&str, Build<f64>); 2] = [
        ("new", NearestNDInterpolator::new),
        ("rescaled", NearestNDInterpolator::rescaled),
    ];
    for (label, build) in scalars {
        let nearest = build(&points, &values)?;
        check(
            &format!("NearestNDInterpolator<3> {label}"),
            &nearest,
            &queries,
            f64::NAN,
        );
        check(
            &format!("NearestNDInterpolator<3> {label}, bounded"),
            &nearest.with_distance_bound(2.5),
            &queries,
            f64::NAN,
        );
    }
    let arrays = NearestNDInterpolator::new(&points, &pairs)?;
    check(
        "NearestNDInterpolator<3, [f64; 2]>",
        &arrays,
        &queries,
        [f64::NAN; 2],
    );
    // From a thread that has made no query before.
    let nearest = NearestNDInterpolator::new(&points, &values)?;
    let fresh = std::thread::spawn(move || {
        allocations(|| {
            for query in &queries {
                std::hint::black_box(nearest.eval(query));
            }
        })
    })
    .join()
    .unwrap_or(usize::MAX);
    assert_eq!(
        fresh, 0,
        "NearestNDInterpolator: the first queries of a thread allocate"
    );

    // RBFInterpolator: each kind of kernel (with and without `epsilon`), with
    // and without polynomial terms, scalar and array values. The queries are
    // near, far, too far for a finite value, and not finite.
    use interpx::{RBFInterpolator, RbfDegree, RbfKernel, RbfOptions};
    // Every ninth point: every seventh would have the same third coordinate.
    let centers: Vec<[f64; 3]> = points.iter().step_by(9).copied().collect();
    let heights: Vec<f64> = (0..centers.len() as u32)
        .map(|n| f64::from(n % 5))
        .collect();
    let pairs: Vec<[f64; 2]> = heights.iter().map(|v| [*v, -*v]).collect();
    for kernel in [
        RbfKernel::ThinPlateSpline,
        RbfKernel::Quintic,
        RbfKernel::Gaussian { epsilon: 0.05 },
        RbfKernel::InverseMultiquadric { epsilon: 0.05 },
    ] {
        for degree in [RbfDegree::Default, RbfDegree::None, RbfDegree::Degree(2)] {
            let options = RbfOptions::new().with_smoothing(0.5).with_degree(degree);
            let rbf = RBFInterpolator::with_options(&centers, &heights, kernel, &options)?;
            check(
                &format!("RBFInterpolator<3> {kernel:?} {degree:?}"),
                &rbf,
                &queries,
                f64::NAN,
            );
        }
        let arrays = RBFInterpolator::new(&centers, &pairs, kernel)?;
        check(
            &format!("RBFInterpolator<3, [f64; 2]> {kernel:?}"),
            &arrays,
            &queries,
            [f64::NAN; 2],
        );
    }
    Ok(())
}

fn main() -> Result<(), InterpolateError> {
    // Inside, outside, and NaN: each takes a different path through `eval`.
    check("Double", &Double, &[0.25, 0.5, 2.0, f64::NAN], f64::NAN);
    check(
        "Triple",
        &Triple,
        &[[0.0, 1.0], [f64::NAN, 0.0], [2.0, 3.0]],
        [f64::NAN; 3],
    );

    // PPoly: every extrapolation mode, since each takes its own path outside the
    // domain. The queries are inside, at a break point, outside on both sides,
    // and NaN. Building the polynomials allocates; that happens before `check`.
    let queries = [0.5, 1.0, 1.75, -3.0, 9.0, f64::NAN];
    for mode in [
        Extrapolate::Forbid,
        Extrapolate::Clamp,
        Extrapolate::Extend,
        Extrapolate::Periodic,
    ] {
        let scalar = PPoly::new(&[0.0, 1.0, 2.0], &[1.0, 0.0, 2.0, 1.0], 1)?.with_extrapolate(mode);
        check(&format!("PPoly<f64> {mode:?}"), &scalar, &queries, f64::NAN);

        let array =
            PPoly::new(&[0.0, 1.0, 2.0], &[[1.0_f64, 2.0], [3.0, 4.0]], 0)?.with_extrapolate(mode);
        check(
            &format!("PPoly<f64, [f64; 2]> {mode:?}"),
            &array,
            &queries,
            [f64::NAN; 2],
        );

        // The hinted form, which the batch uses.
        let n = allocations(|| {
            let mut hint = Hint::new();
            for x in queries {
                std::hint::black_box(scalar.eval_hinted(x, &mut hint));
            }
        });
        assert_eq!(n, 0, "PPoly {mode:?}: eval_hinted allocates");

        // Derivatives of every order up to past the degree, and integrals with
        // each query as a limit: inside, across the domain, and reversed.
        let n = allocations(|| {
            for x in queries {
                for nu in 0..4 {
                    std::hint::black_box(scalar.eval_derivative(x, nu));
                    std::hint::black_box(array.eval_derivative(x, nu));
                }
                for other in [1.25, -7.0, 30.0] {
                    std::hint::black_box(scalar.integrate(x, other));
                    std::hint::black_box(array.integrate(other, x));
                }
            }
        });
        assert_eq!(
            n, 0,
            "PPoly {mode:?}: eval_derivative or integrate allocates"
        );

        // BSpline: the same queries, a scalar and an array value. Degree 3 puts
        // four basis values on the stack; they must stay there.
        let knots = [0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 2.0];
        let scalar = BSpline::new(&knots, &[1.0, -2.0, 0.5, 3.0, 1.5], 3)?.with_extrapolate(mode);
        check(
            &format!("BSpline<f64> {mode:?}"),
            &scalar,
            &queries,
            f64::NAN,
        );

        let array = BSpline::new(&[0.0, 0.0, 1.0, 2.0, 2.0], &[[0.0_f64, 1.0]; 3], 1)?
            .with_extrapolate(mode);
        check(
            &format!("BSpline<f64, [f64; 2]> {mode:?}"),
            &array,
            &queries,
            [f64::NAN; 2],
        );

        let n = allocations(|| {
            let mut hint = Hint::new();
            for x in queries {
                std::hint::black_box(scalar.eval_hinted(x, &mut hint));
            }
        });
        assert_eq!(n, 0, "BSpline {mode:?}: eval_hinted allocates");

        // The integral keeps the basis values of the antiderivative, one more
        // than the spline has, on the stack as well.
        let n = allocations(|| {
            for x in queries {
                for nu in 0..5 {
                    std::hint::black_box(scalar.eval_derivative(x, nu));
                    std::hint::black_box(array.eval_derivative(x, nu));
                }
                for other in [1.25, -7.0, 30.0] {
                    std::hint::black_box(scalar.integrate(x, other));
                    std::hint::black_box(array.integrate(other, x));
                }
            }
        });
        assert_eq!(
            n, 0,
            "BSpline {mode:?}: eval_derivative or integrate allocates"
        );
    }

    // Linear1D and Nearest1D: the same queries and modes, a scalar and an array
    // value, every kind of nearest. The last data point (2.0) is among them:
    // Linear1D returns its value by another path than the formula.
    let at_points = [0.5, 1.0, 1.75, 2.0, -3.0, 9.0, f64::NAN];
    for mode in [
        Extrapolate::Forbid,
        Extrapolate::Clamp,
        Extrapolate::Extend,
        Extrapolate::Periodic,
    ] {
        let scalar = Linear1D::new(&[0.0, 1.0, 2.0], &[1.0, -1.0, 4.0])?.with_extrapolate(mode);
        check(
            &format!("Linear1D<f64> {mode:?}"),
            &scalar,
            &at_points,
            f64::NAN,
        );
        let array = Linear1D::new(&[0.0, 1.0, 2.0], &[[1.0_f64, 2.0], [3.0, 4.0], [0.0, 0.5]])?
            .with_extrapolate(mode);
        check(
            &format!("Linear1D<f64, [f64; 2]> {mode:?}"),
            &array,
            &at_points,
            [f64::NAN; 2],
        );
        let n = allocations(|| {
            let mut hint = Hint::new();
            for x in at_points {
                std::hint::black_box(scalar.eval_hinted(x, &mut hint));
            }
            std::hint::black_box(scalar.as_ppoly());
        });
        assert_eq!(n, 0, "Linear1D {mode:?}: eval_hinted or as_ppoly allocates");

        // CubicSpline: evaluation is that of its PPoly, through the spline's
        // own methods. Building it allocates; that happens before `check`.
        for bc in [BoundaryCondition::NotAKnot, BoundaryCondition::Periodic] {
            let x = [0.0, 0.5, 1.0, 1.5, 2.0];
            let scalar =
                CubicSpline::new(&x, &[1.0, -1.0, 4.0, 0.0, 1.0], bc)?.with_extrapolate(mode);
            check(
                &format!("CubicSpline<f64> {bc:?} {mode:?}"),
                &scalar,
                &at_points,
                f64::NAN,
            );
            let n = allocations(|| {
                let mut hint = Hint::new();
                for x in at_points {
                    std::hint::black_box(scalar.eval_hinted(x, &mut hint));
                }
                std::hint::black_box(scalar.as_ppoly());
            });
            assert_eq!(
                n, 0,
                "CubicSpline {mode:?}: eval_hinted or as_ppoly allocates"
            );
        }
        let columns = [[1.0_f64, 2.0], [3.0, 4.0], [0.0, 0.5], [1.0, 2.0]];
        let array = CubicSpline::new(&[0.0, 1.0, 1.5, 2.0], &columns, BoundaryCondition::Periodic)?
            .with_extrapolate(mode);
        check(
            &format!("CubicSpline<f64, [f64; 2]> {mode:?}"),
            &array,
            &at_points,
            [f64::NAN; 2],
        );

        // Akima1D, both methods: evaluation is that of its PPoly as well.
        for method in [AkimaMethod::Akima, AkimaMethod::Makima] {
            let x = [0.0, 0.5, 1.0, 1.5, 2.0];
            let scalar =
                Akima1D::new(&x, &[1.0, -1.0, 4.0, 4.0, 1.0], method)?.with_extrapolate(mode);
            check(
                &format!("Akima1D<f64> {method:?} {mode:?}"),
                &scalar,
                &at_points,
                f64::NAN,
            );
            let array =
                Akima1D::new(&[0.0, 1.0, 1.5, 2.0], &columns, method)?.with_extrapolate(mode);
            check(
                &format!("Akima1D<f64, [f64; 2]> {method:?} {mode:?}"),
                &array,
                &at_points,
                [f64::NAN; 2],
            );
            let n = allocations(|| {
                let mut hint = Hint::new();
                for x in at_points {
                    std::hint::black_box(scalar.eval_hinted(x, &mut hint));
                }
                std::hint::black_box(scalar.as_ppoly());
            });
            assert_eq!(n, 0, "Akima1D {mode:?}: eval_hinted or as_ppoly allocates");
        }

        // RegularGridInterpolator: the methods that take two points an
        // axis, scalar and array values, 2 and 3 dimensions.
        let (side, other) = ([0.0_f64, 0.5, 1.0, 2.0], [0.0_f64, 1.0, 1.5]);
        let on_grid: Vec<f64> = (0..36).map(|n| f64::from(n % 7) - 2.0).collect();
        let in_grid = [
            [0.5, 1.0, 0.25],
            [1.75, 0.2, 1.5],
            [2.0, 1.5, 0.0],
            [-3.0, 0.5, 0.5],
            [9.0, 9.0, 9.0],
            [f64::NAN, 0.5, 0.5],
        ];
        for method in [GridMethod::Linear, GridMethod::Nearest, GridMethod::Pchip] {
            let cube = RegularGridInterpolator::new(
                [&side[..], &other[..], &other[..]],
                &on_grid,
                method,
            )?
            .with_extrapolate(mode);
            check(
                &format!("RegularGridInterpolator<f64, 3> {method:?} {mode:?}"),
                &cube,
                &in_grid,
                f64::NAN,
            );
            let plane = RegularGridInterpolator::new(
                [&side[..], &other[..]],
                &[[1.0_f64, -1.0]; 12],
                method,
            )?
            .with_extrapolate(mode);
            let in_plane = in_grid.map(|p| [p[0], p[1]]);
            check(
                &format!("RegularGridInterpolator<f64, 2, [f64; 2]> {method:?} {mode:?}"),
                &plane,
                &in_plane,
                [f64::NAN; 2],
            );
            let n = allocations(|| {
                let mut hints = [Hint::new(); 3];
                for point in &in_grid {
                    std::hint::black_box(cube.eval_hinted(point, &mut hints));
                }
                std::hint::black_box((cube.axes(), cube.values(), cube.shape(), cube.domain()));
            });
            assert_eq!(
                n, 0,
                "RegularGridInterpolator {method:?} {mode:?}: eval_hinted or an accessor allocates"
            );
        }

        // The spline methods of the grid: six points an axis for the quintic.
        let long = [0.0_f64, 0.5, 1.0, 1.25, 1.5, 2.0];
        let on_long: Vec<f64> = (0..216).map(|n| f64::from(n % 11) - 5.0).collect();
        for method in [GridMethod::Cubic, GridMethod::Quintic] {
            let cube = RegularGridInterpolator::new([&long[..]; 3], &on_long, method)?
                .with_extrapolate(mode);
            check(
                &format!("RegularGridInterpolator<f64, 3> {method:?} {mode:?}"),
                &cube,
                &in_grid,
                f64::NAN,
            );
            let plane =
                RegularGridInterpolator::new([&long[..]; 2], &[[1.0_f64, -1.0]; 36], method)?
                    .with_extrapolate(mode);
            let in_plane = in_grid.map(|p| [p[0], p[1]]);
            check(
                &format!("RegularGridInterpolator<f64, 2, [f64; 2]> {method:?} {mode:?}"),
                &plane,
                &in_plane,
                [f64::NAN; 2],
            );
            let n = allocations(|| {
                let mut hints = [Hint::new(); 3];
                for point in &in_grid {
                    std::hint::black_box(cube.eval_hinted(point, &mut hints));
                }
            });
            assert_eq!(
                n, 0,
                "RegularGridInterpolator {method:?} {mode:?}: eval_hinted allocates"
            );
        }

        // BSplineCurve: evaluation is that of its BSpline. Parameters chosen
        // from the points, and given.
        let corners = [
            [0.0_f64, 0.0],
            [1.0, 0.0],
            [1.0, 2.0],
            [0.0, 2.0],
            [0.0, 0.0],
        ];
        let chosen = BSplineCurve::new(
            &corners,
            3,
            BoundaryCondition::Periodic,
            Parametrization::ChordLength,
        )?
        .with_extrapolate(mode);
        check(
            &format!("BSplineCurve<f64, 2> {mode:?}"),
            &chosen,
            &at_points,
            [f64::NAN; 2],
        );
        let given = BSplineCurve::with_parameters(
            &[0.0, 0.5, 1.0, 1.5, 2.0],
            &corners,
            2,
            BoundaryCondition::NotAKnot,
        )?
        .with_extrapolate(mode);
        check(
            &format!("BSplineCurve<f64, 2> with parameters {mode:?}"),
            &given,
            &at_points,
            [f64::NAN; 2],
        );
        let n = allocations(|| {
            let mut hint = Hint::new();
            for u in at_points {
                std::hint::black_box(chosen.eval_hinted(u, &mut hint));
            }
            std::hint::black_box((chosen.as_bspline(), chosen.parameters()));
        });
        assert_eq!(
            n, 0,
            "BSplineCurve {mode:?}: eval_hinted, as_bspline, or parameters allocates"
        );

        // BSpline1D: evaluation is that of its BSpline. An even and an odd
        // degree, and a periodic one.
        for degree in [2, 3] {
            let scalar = BSpline1D::new(
                &[0.0, 0.5, 1.0, 1.5, 2.0],
                &[1.0, -1.0, 4.0, 4.0, 1.0],
                degree,
                BoundaryCondition::NotAKnot,
            )?
            .with_extrapolate(mode);
            check(
                &format!("BSpline1D<f64> degree {degree} {mode:?}"),
                &scalar,
                &at_points,
                f64::NAN,
            );
            let array = BSpline1D::new(
                &[0.0, 1.0, 1.5, 2.0],
                &columns,
                degree,
                BoundaryCondition::Periodic,
            )?
            .with_extrapolate(mode);
            check(
                &format!("BSpline1D<f64, [f64; 2]> degree {degree} {mode:?}"),
                &array,
                &at_points,
                [f64::NAN; 2],
            );
            let n = allocations(|| {
                let mut hint = Hint::new();
                for x in at_points {
                    std::hint::black_box(scalar.eval_hinted(x, &mut hint));
                }
                std::hint::black_box(scalar.as_bspline());
            });
            assert_eq!(
                n, 0,
                "BSpline1D degree {degree} {mode:?}: eval_hinted or as_bspline allocates"
            );
        }

        // CubicHermiteSpline: evaluation is that of its PPoly as well.
        let scalar = CubicHermiteSpline::new(
            &[0.0, 0.5, 1.0, 1.5, 2.0],
            &[1.0, -1.0, 4.0, 4.0, 1.0],
            &[0.0, 2.0, -3.0, 0.5, 1.0],
        )?
        .with_extrapolate(mode);
        check(
            &format!("CubicHermiteSpline<f64> {mode:?}"),
            &scalar,
            &at_points,
            f64::NAN,
        );
        let array = CubicHermiteSpline::new(&[0.0, 1.0, 1.5, 2.0], &columns, &columns)?
            .with_extrapolate(mode);
        check(
            &format!("CubicHermiteSpline<f64, [f64; 2]> {mode:?}"),
            &array,
            &at_points,
            [f64::NAN; 2],
        );
        let n = allocations(|| {
            let mut hint = Hint::new();
            for x in at_points {
                std::hint::black_box(scalar.eval_hinted(x, &mut hint));
            }
            std::hint::black_box(scalar.as_ppoly());
        });
        assert_eq!(
            n, 0,
            "CubicHermiteSpline {mode:?}: eval_hinted or as_ppoly allocates"
        );

        // PchipInterpolator: evaluation is that of its PPoly as well.
        let scalar =
            PchipInterpolator::new(&[0.0, 0.5, 1.0, 1.5, 2.0], &[1.0, -1.0, 4.0, 4.0, 1.0])?
                .with_extrapolate(mode);
        check(
            &format!("PchipInterpolator<f64> {mode:?}"),
            &scalar,
            &at_points,
            f64::NAN,
        );
        let array = PchipInterpolator::new(&[0.0, 1.0, 1.5, 2.0], &columns)?.with_extrapolate(mode);
        check(
            &format!("PchipInterpolator<f64, [f64; 2]> {mode:?}"),
            &array,
            &at_points,
            [f64::NAN; 2],
        );
        let n = allocations(|| {
            let mut hint = Hint::new();
            for x in at_points {
                std::hint::black_box(scalar.eval_hinted(x, &mut hint));
            }
            std::hint::black_box(scalar.as_ppoly());
        });
        assert_eq!(
            n, 0,
            "PchipInterpolator {mode:?}: eval_hinted or as_ppoly allocates"
        );

        for kind in [
            NearestKind::Nearest,
            NearestKind::NearestUp,
            NearestKind::Previous,
            NearestKind::Next,
        ] {
            let scalar =
                Nearest1D::new(&[0.0, 1.0, 2.0], &[1.0, -1.0, 4.0], kind)?.with_extrapolate(mode);
            check(
                &format!("Nearest1D<f64> {kind:?} {mode:?}"),
                &scalar,
                &at_points,
                f64::NAN,
            );
            let array = Nearest1D::new(&[0.0, 1.0], &[[1.0_f64, 2.0], [3.0, 4.0]], kind)?
                .with_extrapolate(mode);
            check(
                &format!("Nearest1D<f64, [f64; 2]> {kind:?} {mode:?}"),
                &array,
                &at_points,
                [f64::NAN; 2],
            );
            let n = allocations(|| {
                let mut hint = Hint::new();
                for x in at_points {
                    std::hint::black_box(scalar.eval_hinted(x, &mut hint));
                }
            });
            assert_eq!(n, 0, "Nearest1D {kind:?} {mode:?}: eval_hinted allocates");
        }
    }

    // A batch long and ordered enough for `eval_batch_into` to take its hinted
    // path (the batches above are too short for it), and one in no order.
    let ordered: Vec<f64> = (0..200).map(|i| f64::from(i) * 0.01).collect();
    let unordered: Vec<f64> = (0..200).map(|i| f64::from((i * 37) % 200) * 0.01).collect();
    let poly = PPoly::new(&[0.0, 1.0, 2.0], &[1.0, 0.0, 2.0, 1.0], 1)?;
    let spline = BSpline::new(&[0.0, 0.0, 1.0, 2.0, 2.0], &[0.0, 1.0, 0.0], 1)?;
    let line = Linear1D::new(&[0.0, 1.0, 2.0], &[1.0, -1.0, 4.0])?;
    let nearest = Nearest1D::new(&[0.0, 1.0, 2.0], &[1.0, -1.0, 4.0], NearestKind::Nearest)?;
    let mut out = vec![0.0; 200];
    for points in [&ordered, &unordered] {
        let n = allocations(|| {
            std::hint::black_box(poly.eval_batch_into(points, &mut out, f64::NAN)).ok();
        });
        assert_eq!(n, 0, "PPoly: a long batch allocates");
        let n = allocations(|| {
            std::hint::black_box(spline.eval_batch_into(points, &mut out, f64::NAN)).ok();
        });
        assert_eq!(n, 0, "BSpline: a long batch allocates");
        let n = allocations(|| {
            std::hint::black_box(line.eval_batch_into(points, &mut out, f64::NAN)).ok();
            std::hint::black_box(nearest.eval_batch_into(points, &mut out, f64::NAN)).ok();
        });
        assert_eq!(n, 0, "Linear1D or Nearest1D: a long batch allocates");
    }

    #[cfg(feature = "scattered")]
    scattered()?;

    // The two that build an object do allocate: the counter sees these calls.
    let n = allocations(|| {
        std::hint::black_box(poly.derivative(1)).ok();
        std::hint::black_box(spline.antiderivative(1)).ok();
    });
    assert!(n >= 2, "derivative and antiderivative build an object");

    println!("tests/alloc.rs: the evaluation APIs do not allocate");
    Ok(())
}
