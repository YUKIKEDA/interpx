//! Time to build and to evaluate an `RBFInterpolator` solved over all its
//! points.
//!
//! Points at random in the unit square: 100, 300, and 1000 of them, for three
//! kernels: the thin plate spline (a logarithm an entry), the Gaussian (an
//! exponential), and the linear one (a square root). Construction puts
//! together a matrix with a row for every point and factors it; its time grows
//! with the cube of the points. Evaluation is 1000 queries; each takes one
//! kernel for every point.
//!
//! Run with `just bench`. Not a gate; CI does not run it.

// `criterion_group!` defines a public function without a doc comment, and the
// crate denies `missing_docs`. The lint is about the library's API, not about a
// benchmark binary (docs/design.ja.md section 8.3).
#![allow(missing_docs)]

use criterion::{criterion_group, criterion_main, Criterion, Throughput};
use interpx::{Interpolator, RBFInterpolator, RbfKernel};
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

fn rbf(c: &mut Criterion) {
    let kernels = [
        ("thin_plate_spline", RbfKernel::ThinPlateSpline),
        ("gaussian", RbfKernel::Gaussian { epsilon: 10.0 }),
        ("linear", RbfKernel::Linear),
    ];
    let mut group = c.benchmark_group("rbf");
    group.sample_size(20);
    for n in [100, 300, 1000] {
        let mut rng = Lcg(1);
        let points: Vec<[f64; 2]> = (0..n).map(|_| [rng.next(), rng.next()]).collect();
        let values: Vec<f64> = (0..n).map(|_| 2.0 * rng.next() - 1.0).collect();
        let queries: Vec<[f64; 2]> = (0..1000).map(|_| [rng.next(), rng.next()]).collect();
        let mut out = vec![0.0; queries.len()];
        for (label, kernel) in kernels {
            group.throughput(Throughput::Elements(n as u64));
            group.bench_function(format!("new/{label}/{n}"), |b| {
                b.iter(|| {
                    black_box(RBFInterpolator::new(black_box(&points), &values, kernel).is_ok())
                });
            });
            let Ok(built) = RBFInterpolator::new(&points, &values, kernel) else {
                return;
            };
            group.throughput(Throughput::Elements(queries.len() as u64));
            group.bench_function(format!("eval/{label}/{n}"), |b| {
                b.iter(|| {
                    for (query, slot) in queries.iter().zip(out.iter_mut()) {
                        *slot = built.eval_or(black_box(query), f64::NAN);
                    }
                    black_box(&out);
                });
            });
        }
    }
    group.finish();
}

criterion_group!(benches, rbf);
criterion_main!(benches);
