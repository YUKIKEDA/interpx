//! The parallel batch against the serial one (the `rayon` feature).
//!
//! Three evaluations whose cost a point is far apart, from 100 to 1 000 000
//! points: a piecewise linear function at points in no order (about 5 ns a
//! point), a grid of three dimensions (about 25 ns), and an RBF interpolant
//! with 20 neighbors (about 5 µs; to 10 000 points). A parallel batch times
//! its first points and shares the rest out only when that is worth it, so
//! the small batches show what that decision costs and the large ones what
//! the threads give.
//!
//! Run with `just bench --bench par` (needs `--features rayon`, which the
//! manifest asks for). Not a gate; CI does not run it.

// `criterion_group!` defines a public function without a doc comment, and the
// crate denies `missing_docs`. The lint is about the library's API, not about a
// benchmark binary (docs/design.ja.md section 8.3).
#![allow(missing_docs)]

use criterion::{criterion_group, criterion_main, Criterion, Throughput};
use interpx::{
    GridMethod, Interpolator, Linear1D, ParInterpolator, RBFInterpolator, RbfKernel, RbfOptions,
    RegularGridInterpolator,
};
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

/// The serial and the parallel batch of `interp` on the first `n` of `points`,
/// for each `n` of `sizes`.
fn pair<I>(c: &mut Criterion, name: &str, interp: &I, points: &[I::Point], sizes: &[usize])
where
    I: Interpolator<Output = f64> + Sync,
    I::Point: Sync,
{
    let mut group = c.benchmark_group(format!("par/{name}"));
    group.sample_size(20);
    for n in sizes.iter().copied().filter(|n| *n <= points.len()) {
        let points = &points[..n];
        let mut out = vec![0.0; n];
        group.throughput(Throughput::Elements(n as u64));
        group.bench_function(format!("serial/{n}"), |b| {
            b.iter(|| {
                let done = interp.eval_batch_into(black_box(points), &mut out, f64::NAN);
                black_box((&out, done.is_ok()));
            });
        });
        group.bench_function(format!("parallel/{n}"), |b| {
            b.iter(|| {
                let done = interp.eval_batch_into_par(black_box(points), &mut out, f64::NAN);
                black_box((&out, done.is_ok()));
            });
        });
    }
    group.finish();
}

fn batches(c: &mut Criterion) {
    let mut rng = Lcg(1);
    let sizes = [100, 1000, 10_000, 100_000, 1_000_000];

    let x: Vec<f64> = (0..10_000).map(f64::from).collect();
    let y: Vec<f64> = x.iter().map(|x| (0.01 * x).sin()).collect();
    let queries: Vec<f64> = (0..1_000_000).map(|_| 9999.0 * rng.next()).collect();
    if let Ok(line) = Linear1D::new(&x, &y) {
        pair(c, "linear", &line, &queries, &sizes);
    }

    let axis: Vec<f64> = (0..40).map(|i| f64::from(i) + 0.3 * rng.next()).collect();
    let values: Vec<f64> = (0..64_000).map(|_| rng.next()).collect();
    let span = axis[39] - axis[0];
    let queries: Vec<[f64; 3]> = (0..1_000_000)
        .map(|_| core::array::from_fn(|_| axis[0] + span * rng.next()))
        .collect();
    if let Ok(grid) = RegularGridInterpolator::new([&axis[..]; 3], &values, GridMethod::Linear) {
        pair(c, "grid_3d", &grid, &queries, &sizes);
    }

    let data: Vec<[f64; 2]> = (0..10_000).map(|_| [rng.next(), rng.next()]).collect();
    let values: Vec<f64> = data.iter().map(|p| (3.0 * p[0]).sin() * p[1]).collect();
    let queries: Vec<[f64; 2]> = (0..10_000).map(|_| [rng.next(), rng.next()]).collect();
    let options = RbfOptions::new().with_neighbors(20);
    let built = RBFInterpolator::with_options(&data, &values, RbfKernel::ThinPlateSpline, &options);
    if let Ok(rbf) = built {
        pair(c, "rbf_neighbors", &rbf, &queries, &sizes[..3]);
    }
}

criterion_group!(benches, batches);
criterion_main!(benches);
