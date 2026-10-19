//! The same evaluation called from two places of one program (P-9, #132;
//! P-10, #138).
//!
//! The time of an evaluation at points in no order once moved by a factor of
//! up to two with where the search by halves was built in: `Nearest1D` on 100
//! points not evenly spaced, 40 000 queries in no order, took 328 µs from one
//! place of a program and 707 µs from another (#137). This program evaluates
//! one interpolator from two places, a function on its concrete type and a
//! function generic over `Interpolator`, and prints both times, so that the
//! check can be made again after a change of the search, of the compiler, or
//! of the machine.
//!
//! Run with `cargo bench --bench search_places`. Not a gate; CI does not run
//! it, and it uses no harness: each time is the median of repeated batches.

// The crate denies `missing_docs`; that lint is about the library's API, not a
// benchmark binary (docs/design.ja.md section 8.3).
#![allow(missing_docs)]

use interpx::{BSpline1D, BoundaryCondition, Interpolator, Linear1D, Nearest1D, NearestKind};
use std::hint::black_box;
use std::time::{Duration, Instant};

/// A small generator, so that the inputs are the same on every run.
struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 11) as f64 / (1_u64 << 53) as f64
    }
}

/// Place one: the batch of the concrete type.
#[inline(never)]
fn concrete(interpolator: &Linear1D<f64>, points: &[f64], out: &mut [f64]) {
    let done = interpolator.eval_batch_into(points, out, f64::NAN);
    black_box(done.is_ok());
}

/// Place two: the same batch through a function generic over the trait.
#[inline(never)]
fn generic<I: Interpolator<Point = f64, Output = f64>>(
    interpolator: &I,
    points: &[f64],
    out: &mut [f64],
) {
    let done = interpolator.eval_batch_into(points, out, f64::NAN);
    black_box(done.is_ok());
}

/// The median time of one call of `f`, after a warm-up.
fn median(mut f: impl FnMut()) -> Duration {
    let start = Instant::now();
    while start.elapsed() < Duration::from_millis(300) {
        f();
    }
    let mut samples: Vec<Duration> = (0..31)
        .map(|_| {
            let t = Instant::now();
            f();
            t.elapsed()
        })
        .collect();
    samples.sort();
    samples[15]
}

fn main() {
    let mut rng = Lcg(1);
    // 100 points, not evenly spaced: increasing gaps of random width.
    let mut x = vec![0.0];
    for _ in 1..100 {
        let gap = 0.5 + rng.next();
        x.push(x.last().copied().unwrap_or(0.0) + gap);
    }
    let y: Vec<f64> = x.iter().map(|t| (0.1 * t).sin()).collect();
    let high = x[x.len() - 1];
    let points: Vec<f64> = (0..40_000).map(|_| rng.next() * high).collect();
    let mut out = vec![0.0; points.len()];

    let (Ok(linear), Ok(nearest), Ok(bspline)) = (
        Linear1D::new(&x, &y),
        Nearest1D::new(&x, &y, NearestKind::Nearest),
        BSpline1D::new(&x, &y, 1, BoundaryCondition::NotAKnot),
    ) else {
        return;
    };

    let a = median(|| concrete(&linear, black_box(&points), &mut out));
    let b = median(|| generic(&linear, black_box(&points), &mut out));
    println!("Linear1D     concrete {a:>10.2?}  generic {b:>10.2?}");
    let b = median(|| generic(&nearest, black_box(&points), &mut out));
    let c = median(|| {
        let done = nearest.eval_batch_into(black_box(&points), &mut out, f64::NAN);
        black_box(done.is_ok());
    });
    println!("Nearest1D    inline   {c:>10.2?}  generic {b:>10.2?}");
    let b = median(|| generic(&bspline, black_box(&points), &mut out));
    let c = median(|| {
        let done = bspline.eval_batch_into(black_box(&points), &mut out, f64::NAN);
        black_box(done.is_ok());
    });
    println!("BSpline1D(1) inline   {c:>10.2?}  generic {b:>10.2?}");
}
