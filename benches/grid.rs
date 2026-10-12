//! Time to evaluate a `RegularGridInterpolator` (`Linear` and `Nearest`).
//!
//! Grids in 2, 3, and 4 dimensions (200, 40, and 10 points an axis, unevenly
//! spaced), 1000 points each: along a path that crosses the grid from corner
//! to corner, where a point is in the cell of the one before or next to it,
//! and the same points shuffled. `eval` searches every axis for every point;
//! `eval_hinted` starts each search where the last one ended.
//!
//! Construction copies the axes and the values and solves nothing; it is timed
//! once per dimension, with one axis given in descending order (its values are
//! reordered).
//!
//! Run with `just bench`. Not a gate; CI does not run it.

// `criterion_group!` defines a public function without a doc comment, and the
// crate denies `missing_docs`. The lint is about the library's API, not about a
// benchmark binary (docs/design.ja.md section 8.3).
#![allow(missing_docs)]

use criterion::{criterion_group, criterion_main, Criterion, Throughput};
use interpx::{GridMethod, Hint, Interpolator, RegularGridInterpolator};
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

/// An axis of `n` points from 0 on, at random spacings.
fn axis(rng: &mut Lcg, n: usize) -> Vec<f64> {
    let mut x = vec![0.0];
    for _ in 1..n {
        let next = x[x.len() - 1] + 0.5 + rng.next();
        x.push(next);
    }
    x
}

fn dimension<const D: usize>(c: &mut Criterion, points_per_axis: usize) {
    let mut rng = Lcg(1);
    let axes: [Vec<f64>; D] = core::array::from_fn(|_| axis(&mut rng, points_per_axis));
    let count = points_per_axis.pow(D as u32);
    let values: Vec<f64> = (0..count).map(|_| 2.0 * rng.next() - 1.0).collect();
    let slices: [&[f64]; D] = core::array::from_fn(|d| axes[d].as_slice());

    let n = 1000;
    let path: Vec<[f64; D]> = (0..n)
        .map(|i| {
            let t = (i as f64 + 0.5) / n as f64;
            core::array::from_fn(|d| t * axes[d][points_per_axis - 1])
        })
        .collect();
    let mut shuffled = path.clone();
    for i in (1..n).rev() {
        let j = (rng.next() * (i + 1) as f64) as usize;
        shuffled.swap(i, j.min(i));
    }

    let mut group = c.benchmark_group(format!("grid_{D}d"));
    group.throughput(Throughput::Elements(n as u64));
    let mut out = vec![0.0; n];
    for (label, method) in [
        ("linear", GridMethod::Linear),
        ("nearest", GridMethod::Nearest),
    ] {
        let Ok(grid) = RegularGridInterpolator::new(slices, &values, method) else {
            return;
        };
        for (order, points) in [("path", &path), ("shuffled", &shuffled)] {
            group.bench_function(format!("{label}/eval/{order}"), |b| {
                b.iter(|| {
                    for (point, slot) in points.iter().zip(out.iter_mut()) {
                        *slot = grid.eval_or(black_box(point), f64::NAN);
                    }
                    black_box(&out);
                });
            });
            group.bench_function(format!("{label}/eval_hinted/{order}"), |b| {
                b.iter(|| {
                    let mut hints = [Hint::new(); D];
                    for (point, slot) in points.iter().zip(out.iter_mut()) {
                        *slot = grid
                            .eval_hinted(black_box(point), &mut hints)
                            .unwrap_or(f64::NAN);
                    }
                    black_box(&out);
                });
            });
            group.bench_function(format!("{label}/eval_batch_into/{order}"), |b| {
                b.iter(|| {
                    let done = grid.eval_batch_into(black_box(points), &mut out, f64::NAN);
                    black_box((&out, done.is_ok()));
                });
            });
        }
    }
    // Construction, with the last axis given backward.
    let mut backward = axes.clone();
    backward[D - 1].reverse();
    let given: [&[f64]; D] = core::array::from_fn(|d| backward[d].as_slice());
    group.throughput(Throughput::Elements(count as u64));
    group.bench_function("new/ascending", |b| {
        b.iter(|| {
            let built =
                RegularGridInterpolator::new(slices, black_box(&values), GridMethod::Linear);
            black_box(built.is_ok())
        });
    });
    group.bench_function("new/one_axis_descending", |b| {
        b.iter(|| {
            let built = RegularGridInterpolator::new(given, black_box(&values), GridMethod::Linear);
            black_box(built.is_ok())
        });
    });
    group.finish();
}

fn grids(c: &mut Criterion) {
    dimension::<2>(c, 200);
    dimension::<3>(c, 40);
    dimension::<4>(c, 10);
}

criterion_group!(benches, grids);
criterion_main!(benches);
