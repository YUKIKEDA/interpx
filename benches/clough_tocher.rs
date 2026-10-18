//! Time to build and to evaluate a `CloughTocher2DInterpolator`.
//!
//! Points at random in the unit square: 100, 1000, and 100 000 of them.
//! Construction triangulates the points (exact predicates), inverts the matrix
//! of every triangle, and solves for the gradients at the points (conjugate
//! gradients until the residual is 1e-14 of the right-hand side). Evaluation is
//! 1000 queries at random in the square (some outside the hull), one by one and
//! as a batch (which is a loop over the one).
//!
//! Run with `just bench`. Not a gate; CI does not run it.

// `criterion_group!` defines a public function without a doc comment, and the
// crate denies `missing_docs`. The lint is about the library's API, not about a
// benchmark binary (docs/design.ja.md section 8.3).
#![allow(missing_docs)]

use criterion::{criterion_group, criterion_main, Criterion, Throughput};
use interpx::{CloughTocher2DInterpolator, Interpolator};
use std::hint::black_box;

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

fn all(c: &mut Criterion) {
    let mut group = c.benchmark_group("clough_tocher_2d");
    for n in [100, 1000, 100_000] {
        let mut rng = Lcg(1);
        let points: Vec<[f64; 2]> = (0..n).map(|_| [rng.next(), rng.next()]).collect();
        // A smooth function, so that the gradients are those of a surface and
        // not of noise.
        let values: Vec<f64> = points
            .iter()
            .map(|[x, y]| (3.0 * x).sin() * (2.0 * y).cos())
            .collect();
        let queries: Vec<[f64; 2]> = (0..1000).map(|_| [rng.next(), rng.next()]).collect();
        group.throughput(Throughput::Elements(n as u64));
        // A construction of 100 000 points takes long enough that a few
        // samples say what it takes.
        group.sample_size(if n <= 1000 { 100 } else { 10 });
        group.bench_function(format!("new/{n}"), |b| {
            b.iter(|| {
                black_box(CloughTocher2DInterpolator::new(black_box(&points), &values).is_ok())
            });
        });
        group.sample_size(100);
        let Ok(clough) = CloughTocher2DInterpolator::new(&points, &values) else {
            return;
        };
        let mut out = vec![0.0; queries.len()];
        group.throughput(Throughput::Elements(queries.len() as u64));
        group.bench_function(format!("eval/{n}"), |b| {
            b.iter(|| {
                for (query, slot) in queries.iter().zip(out.iter_mut()) {
                    *slot = clough.eval_or(black_box(query), f64::NAN);
                }
                black_box(&out);
            });
        });
        group.bench_function(format!("eval_batch_into/{n}"), |b| {
            b.iter(|| {
                let done = clough.eval_batch_into(black_box(&queries), &mut out, f64::NAN);
                black_box((&out, done.is_ok()));
            });
        });
    }
    group.finish();
}

criterion_group!(benches, all);
criterion_main!(benches);
