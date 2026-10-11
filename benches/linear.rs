//! Time to evaluate a `Linear1D`.
//!
//! The shape follows `benches/ppoly.rs`: 55 intervals, evaluated at 10 to 10000
//! points in increasing order and shuffled, by a loop over `eval`, the same loop
//! with a hint, and the batch API.
//!
//! The group `linear_vs_ppoly` times the same polygon through `as_ppoly()`. The
//! two differ in one comparison per evaluation, which `Linear1D` makes to return
//! the last data value as it is (docs/design.ja.md section 6.1).
//!
//! Run with `just bench`. Not a gate; CI does not run it.

// `criterion_group!` defines a public function without a doc comment, and the
// crate denies `missing_docs`. The lint is about the library's API, not about a
// benchmark binary (docs/design.ja.md section 8.3).
#![allow(missing_docs)]

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use interpx::{Hint, Interpolator, Linear1D};
use std::hint::black_box;

const INTERVALS: usize = 55;
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

fn polygon() -> Option<Linear1D<f64>> {
    let mut rng = Lcg(1);
    let mut x = vec![0.0];
    for _ in 0..INTERVALS {
        let next = x[x.len() - 1] + 0.5 + rng.next();
        x.push(next);
    }
    let y: Vec<f64> = (0..=INTERVALS).map(|_| 2.0 * rng.next() - 1.0).collect();
    Linear1D::new(&x, &y).ok()
}

/// `n` points across the data in increasing order, and the same points shuffled.
fn points(f: &Linear1D<f64>, n: usize) -> (Vec<f64>, Vec<f64>) {
    let (lo, hi) = f.domain();
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
    let Some(f) = polygon() else {
        return;
    };
    let mut group = c.benchmark_group("linear");
    for n in SIZES {
        let (sorted, shuffled) = points(&f, n);
        group.throughput(Throughput::Elements(n as u64));
        for (order, xs) in [("sorted", &sorted), ("shuffled", &shuffled)] {
            let mut out = vec![0.0; n];

            group.bench_with_input(BenchmarkId::new(format!("eval/{order}"), n), xs, |b, xs| {
                b.iter(|| {
                    for (x, slot) in xs.iter().zip(out.iter_mut()) {
                        *slot = f.eval_or(black_box(x), f64::NAN);
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
                            *slot = f.eval_hinted(black_box(*x), &mut hint).unwrap_or(f64::NAN);
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
                        let done = f.eval_batch_into(black_box(xs), &mut out, f64::NAN);
                        black_box((&out, done.is_ok()));
                    });
                },
            );
        }
    }
    group.finish();
}

/// The same polygon through `Linear1D` and through its `PPoly`, at 1000 points.
fn against_ppoly(c: &mut Criterion) {
    let Some(f) = polygon() else {
        return;
    };
    let poly = f.as_ppoly();
    let n = 1000;
    let (sorted, shuffled) = points(&f, n);
    let mut out = vec![0.0; n];
    let mut group = c.benchmark_group("linear_vs_ppoly");
    group.throughput(Throughput::Elements(n as u64));
    for (order, xs) in [("sorted", &sorted), ("shuffled", &shuffled)] {
        group.bench_function(format!("linear/eval/{order}"), |b| {
            b.iter(|| {
                for (x, slot) in xs.iter().zip(out.iter_mut()) {
                    *slot = f.eval_or(black_box(x), f64::NAN);
                }
                black_box(&out);
            });
        });
        group.bench_function(format!("ppoly/eval/{order}"), |b| {
            b.iter(|| {
                for (x, slot) in xs.iter().zip(out.iter_mut()) {
                    *slot = poly.eval_or(black_box(x), f64::NAN);
                }
                black_box(&out);
            });
        });
        group.bench_function(format!("linear/eval_batch_into/{order}"), |b| {
            b.iter(|| {
                let done = f.eval_batch_into(black_box(xs), &mut out, f64::NAN);
                black_box((&out, done.is_ok()));
            });
        });
        group.bench_function(format!("ppoly/eval_batch_into/{order}"), |b| {
            b.iter(|| {
                let done = poly.eval_batch_into(black_box(xs), &mut out, f64::NAN);
                black_box((&out, done.is_ok()));
            });
        });
    }
    group.finish();
}

criterion_group!(benches, bench, against_ppoly);
criterion_main!(benches);
