//! Allocation counts of the evaluation APIs (docs/design.ja.md section 8.1).
//!
//! `eval`, `eval_or`, and `eval_batch_into` allocate nothing, nor do
//! `eval_derivative` and `integrate`. `eval_batch` returns a `Vec`, so it
//! allocates, as do `derivative` and `antiderivative`, which build an object.
//!
//! The counter is global to the process. Everything is checked from the one test
//! function below, so no other test thread allocates while a region is measured.
//! A row that adds an interpolator adds it to `evaluation_apis_do_not_allocate`.

use interpx::{BSpline, Extrapolate, Hint, Interpolator, PPoly};
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

#[test]
fn evaluation_apis_do_not_allocate() {
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
        let scalar = PPoly::new(&[0.0, 1.0, 2.0], &[1.0, 0.0, 2.0, 1.0], 1)
            .expect("valid")
            .with_extrapolate(mode);
        check(&format!("PPoly<f64> {mode:?}"), &scalar, &queries, f64::NAN);

        let array = PPoly::new(&[0.0, 1.0, 2.0], &[[1.0_f64, 2.0], [3.0, 4.0]], 0)
            .expect("valid")
            .with_extrapolate(mode);
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
        let scalar = BSpline::new(&knots, &[1.0, -2.0, 0.5, 3.0, 1.5], 3)
            .expect("valid")
            .with_extrapolate(mode);
        check(
            &format!("BSpline<f64> {mode:?}"),
            &scalar,
            &queries,
            f64::NAN,
        );

        let array = BSpline::new(&[0.0, 0.0, 1.0, 2.0, 2.0], &[[0.0_f64, 1.0]; 3], 1)
            .expect("valid")
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

    // A batch long and ordered enough for `eval_batch_into` to take its hinted
    // path (the batches above are too short for it), and one in no order.
    let ordered: Vec<f64> = (0..200).map(|i| f64::from(i) * 0.01).collect();
    let unordered: Vec<f64> = (0..200).map(|i| f64::from((i * 37) % 200) * 0.01).collect();
    let poly = PPoly::new(&[0.0, 1.0, 2.0], &[1.0, 0.0, 2.0, 1.0], 1).expect("valid");
    let spline = BSpline::new(&[0.0, 0.0, 1.0, 2.0, 2.0], &[0.0, 1.0, 0.0], 1).expect("valid");
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
    }

    // The two that build an object do allocate: the counter sees these calls.
    let n = allocations(|| {
        std::hint::black_box(poly.derivative(1)).ok();
        std::hint::black_box(spline.antiderivative(1)).ok();
    });
    assert!(n >= 2, "derivative and antiderivative build an object");
}
