//! Time to build and to evaluate a `NearestNDInterpolator`.
//!
//! Points at random in the unit cube, in 2 and 3 dimensions: 100, 1000, and
//! 100 000 of them. Construction looks for points with the same coordinates, to
//! take them as one, and builds the k-d tree; `new_lattice` times it on a
//! square lattice, where every coordinate is shared by many points. Evaluation
//! is 1000 queries at random, one by one and as a batch (which is a loop over
//! the one).
//!
//! Run with `just bench`. Not a gate; CI does not run it.

// `criterion_group!` defines a public function without a doc comment, and the
// crate denies `missing_docs`. The lint is about the library's API, not about a
// benchmark binary (docs/design.ja.md section 8.3).
#![allow(missing_docs)]

use criterion::{criterion_group, criterion_main, Criterion, Throughput};
use interpx::{Interpolator, NearestNDInterpolator};
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

fn dimension<const D: usize>(c: &mut Criterion) {
    let mut group = c.benchmark_group(format!("nearest_nd_{D}d"));
    for n in [100, 1000, 100_000] {
        let mut rng = Lcg(1);
        let points: Vec<[f64; D]> = (0..n)
            .map(|_| core::array::from_fn(|_| rng.next()))
            .collect();
        let values: Vec<f64> = (0..n).map(|_| rng.next()).collect();
        let queries: Vec<[f64; D]> = (0..1000)
            .map(|_| core::array::from_fn(|_| rng.next()))
            .collect();
        group.throughput(Throughput::Elements(n as u64));
        group.bench_function(format!("new/{n}"), |b| {
            b.iter(|| black_box(NearestNDInterpolator::new(black_box(&points), &values).is_ok()));
        });
        group.bench_function(format!("rescaled/{n}"), |b| {
            b.iter(|| {
                black_box(NearestNDInterpolator::rescaled(black_box(&points), &values).is_ok())
            });
        });
        let Ok(nearest) = NearestNDInterpolator::new(&points, &values) else {
            return;
        };
        let bounded = nearest.clone().with_distance_bound(0.05);
        let mut out = vec![0.0; queries.len()];
        group.throughput(Throughput::Elements(queries.len() as u64));
        group.bench_function(format!("eval/{n}"), |b| {
            b.iter(|| {
                for (query, slot) in queries.iter().zip(out.iter_mut()) {
                    *slot = nearest.eval_or(black_box(query), f64::NAN);
                }
                black_box(&out);
            });
        });
        group.bench_function(format!("eval_bounded/{n}"), |b| {
            b.iter(|| {
                for (query, slot) in queries.iter().zip(out.iter_mut()) {
                    *slot = bounded.eval_or(black_box(query), f64::NAN);
                }
                black_box(&out);
            });
        });
        group.bench_function(format!("eval_batch_into/{n}"), |b| {
            b.iter(|| {
                let done = nearest.eval_batch_into(black_box(&queries), &mut out, f64::NAN);
                black_box((&out, done.is_ok()));
            });
        });
    }
    group.finish();
}

/// Construction on a square lattice of 100 000 points, each given twice.
fn lattice(c: &mut Criterion) {
    let side = 224;
    let mut points: Vec<[f64; 2]> = Vec::new();
    for i in 0..side {
        for j in 0..side {
            points.push([f64::from(i), f64::from(j)]);
        }
    }
    let once = points.clone();
    points.extend(once.iter().copied());
    let values = vec![0.0; points.len()];
    let mut group = c.benchmark_group("nearest_nd_2d");
    group.throughput(Throughput::Elements(points.len() as u64));
    group.bench_function("new_lattice_twice/100352", |b| {
        b.iter(|| black_box(NearestNDInterpolator::new(black_box(&points), &values).is_ok()));
    });
    group.finish();
}

fn all(c: &mut Criterion) {
    dimension::<2>(c);
    dimension::<3>(c);
    lattice(c);
}

criterion_group!(benches, all);
criterion_main!(benches);
