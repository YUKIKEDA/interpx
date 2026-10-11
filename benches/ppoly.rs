//! Time to evaluate a `PPoly`.
//!
//! The shape follows SciPy's `BenchPPoly` (docs/design.ja.md section 8.4.3): 55
//! cubic intervals, evaluated at 10 to 10000 points. Three ways to evaluate are
//! timed, on points in increasing order and on the same points shuffled:
//!
//! - `eval`: a loop over single evaluations (binary search each time).
//! - `eval_hinted`: the same loop with a caller-held hint.
//! - `eval_batch_into`: the batch API, which holds its own hint.
//!
//! Run with `just bench`. Not a gate; CI does not run it.

// `criterion_group!` defines a public function without a doc comment, and the
// crate denies `missing_docs`. The lint is about the library's API, not about a
// benchmark binary (docs/design.ja.md section 8.3).
#![allow(missing_docs)]

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use interpx::{Hint, Interpolator, PPoly};
use std::hint::black_box;

const INTERVALS: usize = 55;
const DEGREE: usize = 3;
const SIZES: [usize; 4] = [10, 100, 1000, 10000];

/// A small generator, so that the inputs are the same on every run and machine.
struct Lcg(u64);

impl Lcg {
    /// The next number in [0, 1).
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 11) as f64 / (1_u64 << 53) as f64
    }
}

fn polynomial() -> Option<PPoly<f64>> {
    let mut rng = Lcg(1);
    let mut breakpoints = vec![0.0];
    for _ in 0..INTERVALS {
        let next = breakpoints[breakpoints.len() - 1] + 0.5 + rng.next();
        breakpoints.push(next);
    }
    let coefficients = (0..INTERVALS * (DEGREE + 1))
        .map(|_| 2.0 * rng.next() - 1.0)
        .collect();
    PPoly::from_owned(breakpoints, coefficients, DEGREE).ok()
}

/// `n` points across the domain in increasing order, and the same points shuffled.
fn points(p: &PPoly<f64>, n: usize) -> (Vec<f64>, Vec<f64>) {
    let (lo, hi) = (p.breakpoints()[0], p.breakpoints()[INTERVALS]);
    let sorted: Vec<f64> = (0..n)
        .map(|i| lo + (hi - lo) * (i as f64 + 0.5) / n as f64)
        .collect();
    let mut shuffled = sorted.clone();
    let mut rng = Lcg(2);
    for i in (1..n).rev() {
        let j = (rng.next() * (i + 1) as f64) as usize;
        shuffled.swap(i, j.min(i));
    }
    (sorted, shuffled)
}

fn bench(c: &mut Criterion) {
    let Some(p) = polynomial() else {
        return;
    };
    let mut group = c.benchmark_group("ppoly");
    for n in SIZES {
        let (sorted, shuffled) = points(&p, n);
        group.throughput(Throughput::Elements(n as u64));
        for (order, xs) in [("sorted", &sorted), ("shuffled", &shuffled)] {
            let mut out = vec![0.0; n];

            group.bench_with_input(BenchmarkId::new(format!("eval/{order}"), n), xs, |b, xs| {
                b.iter(|| {
                    for (x, slot) in xs.iter().zip(out.iter_mut()) {
                        *slot = p.eval_or(black_box(x), f64::NAN);
                    }
                    black_box(&out);
                });
            });

            group.bench_with_input(
                BenchmarkId::new(format!("eval_hinted/{order}"), n),
                xs,
                |b, xs| {
                    b.iter(|| {
                        let mut hint = Hint::new();
                        for (x, slot) in xs.iter().zip(out.iter_mut()) {
                            *slot = p.eval_hinted(black_box(*x), &mut hint).unwrap_or(f64::NAN);
                        }
                        black_box(&out);
                    });
                },
            );

            group.bench_with_input(
                BenchmarkId::new(format!("eval_batch_into/{order}"), n),
                xs,
                |b, xs| {
                    b.iter(|| {
                        let done = p.eval_batch_into(black_box(xs), &mut out, f64::NAN);
                        black_box((&out, done.is_ok()));
                    });
                },
            );
        }
    }
    group.finish();
}

/// Derivatives and integrals.
///
/// The first derivative at 1000 points in increasing order, three ways: a loop
/// over `eval_derivative`, a loop over `eval` on the derivative built ahead, and
/// `eval_batch_into` on it. Then the cost of building the derivative and the
/// antiderivative, and of one definite integral over the whole domain and over
/// a hundredth of it.
fn calculus(c: &mut Criterion) {
    let Some(p) = polynomial() else {
        return;
    };
    let Ok(slope) = p.derivative(1) else {
        return;
    };
    let n = 1000;
    let (xs, _) = points(&p, n);
    let mut out = vec![0.0; n];
    let mut group = c.benchmark_group("ppoly_calculus");

    group.throughput(Throughput::Elements(n as u64));
    group.bench_function("eval_derivative", |b| {
        b.iter(|| {
            for (x, slot) in xs.iter().zip(out.iter_mut()) {
                *slot = p.eval_derivative(black_box(*x), 1).unwrap_or(f64::NAN);
            }
            black_box(&out);
        });
    });
    group.bench_function("derivative_then_eval", |b| {
        b.iter(|| {
            for (x, slot) in xs.iter().zip(out.iter_mut()) {
                *slot = slope.eval_or(black_box(x), f64::NAN);
            }
            black_box(&out);
        });
    });
    group.bench_function("derivative_then_eval_batch_into", |b| {
        b.iter(|| {
            let done = slope.eval_batch_into(black_box(&xs), &mut out, f64::NAN);
            black_box((&out, done.is_ok()));
        });
    });

    group.throughput(Throughput::Elements(1));
    group.bench_function("derivative", |b| {
        b.iter(|| black_box(black_box(&p).derivative(1).is_ok()));
    });
    group.bench_function("antiderivative", |b| {
        b.iter(|| black_box(black_box(&p).antiderivative(1).is_ok()));
    });
    let (lo, hi) = (xs[0], xs[n - 1]);
    group.bench_function("integrate/whole", |b| {
        b.iter(|| black_box(p.integrate(black_box(lo), black_box(hi))));
    });
    let near = lo + (hi - lo) / 100.0;
    group.bench_function("integrate/hundredth", |b| {
        b.iter(|| black_box(p.integrate(black_box(lo), black_box(near))));
    });
    group.finish();
}

criterion_group!(benches, bench, calculus);
criterion_main!(benches);
