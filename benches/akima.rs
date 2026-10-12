//! Time to build and to evaluate an `Akima1D`.
//!
//! Construction is timed for 10 to 10000 points and both methods, beside
//! `PchipInterpolator` and `CubicSpline` on the same data. Akima goes over the
//! points twice: once for the cutoff of its weights, once for the slopes.
//!
//! Evaluation follows `benches/ppoly.rs`: 55 intervals, 1000 points in
//! increasing order and shuffled.
//!
//! Run with `just bench`. Not a gate; CI does not run it.

// `criterion_group!` defines a public function without a doc comment, and the
// crate denies `missing_docs`. The lint is about the library's API, not about a
// benchmark binary (docs/design.ja.md section 8.3).
#![allow(missing_docs)]

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use interpx::{
    Akima1D, AkimaMethod, BoundaryCondition, CubicSpline, Hint, Interpolator, PchipInterpolator,
};
use std::hint::black_box;

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

/// `n` points at random spacings, with values that rise and fall.
fn data(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut rng = Lcg(1);
    let mut x = vec![0.0];
    for _ in 1..n {
        let next = x[x.len() - 1] + 0.5 + rng.next();
        x.push(next);
    }
    let y = (0..n).map(|_| 2.0 * rng.next() - 1.0).collect();
    (x, y)
}

fn construction(c: &mut Criterion) {
    let mut group = c.benchmark_group("akima_new");
    for n in SIZES {
        let (x, y) = data(n);
        let input = (&x, &y);
        group.throughput(Throughput::Elements(n as u64));
        for (name, method) in [
            ("akima", AkimaMethod::Akima),
            ("makima", AkimaMethod::Makima),
        ] {
            group.bench_with_input(BenchmarkId::new(name, n), &input, |b, (x, y)| {
                b.iter(|| black_box(Akima1D::new(black_box(x), black_box(y), method).is_ok()));
            });
        }
        group.bench_with_input(BenchmarkId::new("pchip", n), &input, |b, (x, y)| {
            b.iter(|| black_box(PchipInterpolator::new(black_box(x), black_box(y)).is_ok()));
        });
        group.bench_with_input(BenchmarkId::new("cubic_spline", n), &input, |b, (x, y)| {
            let bc = BoundaryCondition::NotAKnot;
            b.iter(|| black_box(CubicSpline::new(black_box(x), black_box(y), bc).is_ok()));
        });
    }
    group.finish();
}

fn evaluation(c: &mut Criterion) {
    let (x, y) = data(56);
    let Ok(akima) = Akima1D::new(&x, &y, AkimaMethod::Makima) else {
        return;
    };
    let n = 1000;
    let (lo, hi) = akima.domain();
    let sorted: Vec<f64> = (0..n)
        .map(|i| lo + (hi - lo) * (i as f64 + 0.5) / n as f64)
        .collect();
    let mut shuffled = sorted.clone();
    let mut rng = Lcg(2);
    for i in (1..n).rev() {
        let j = (rng.next() * (i + 1) as f64) as usize;
        shuffled.swap(i, j.min(i));
    }
    let mut out = vec![0.0; n];
    let mut group = c.benchmark_group("akima");
    group.throughput(Throughput::Elements(n as u64));
    for (order, xs) in [("sorted", &sorted), ("shuffled", &shuffled)] {
        group.bench_function(format!("eval/{order}"), |b| {
            b.iter(|| {
                for (x, slot) in xs.iter().zip(out.iter_mut()) {
                    *slot = akima.eval_or(black_box(x), f64::NAN);
                }
                black_box(&out);
            });
        });
        group.bench_function(format!("eval_hinted/{order}"), |b| {
            b.iter(|| {
                let mut hint = Hint::new();
                for (x, slot) in xs.iter().zip(out.iter_mut()) {
                    *slot = akima
                        .eval_hinted(black_box(*x), &mut hint)
                        .unwrap_or(f64::NAN);
                }
                black_box(&out);
            });
        });
        group.bench_function(format!("eval_batch_into/{order}"), |b| {
            b.iter(|| {
                let done = akima.eval_batch_into(black_box(xs), &mut out, f64::NAN);
                black_box((&out, done.is_ok()));
            });
        });
    }
    group.finish();
}

criterion_group!(benches, construction, evaluation);
criterion_main!(benches);
