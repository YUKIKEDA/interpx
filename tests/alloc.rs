//! Allocation counts of the evaluation APIs (docs/design.ja.md section 8.1).
//!
//! `eval`, `eval_or`, and `eval_batch_into` allocate nothing. `eval_batch`
//! returns a `Vec`, so it allocates.
//!
//! The counter is global to the process. Everything is checked from the one test
//! function below, so no other test thread allocates while a region is measured.
//! A row that adds an interpolator adds it to `evaluation_apis_do_not_allocate`.

use interpx::Interpolator;
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
}
