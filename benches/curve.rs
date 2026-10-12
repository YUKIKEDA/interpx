//! Time to build and to evaluate a `BSplineCurve`.
//!
//! Construction is timed for 10 to 10000 points in 2 and in 3 dimensions, of
//! degree 3, for each way of choosing the parameters. `BSpline1D` with the
//! same points and the chord-length parameters given stands beside them: what
//! a curve adds is the parameters (a square root per point) and their copy.
//!
//! Evaluation follows `benches/bspline.rs`: 55 intervals, 1000 parameters in
//! increasing order and shuffled, in 3 dimensions.
//!
//! Run with `just bench`. Not a gate; CI does not run it.

// `criterion_group!` defines a public function without a doc comment, and the
// crate denies `missing_docs`. The lint is about the library's API, not about a
// benchmark binary (docs/design.ja.md section 8.3).
#![allow(missing_docs)]

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use interpx::{BSpline1D, BSplineCurve, BoundaryCondition, Hint, Interpolator, Parametrization};
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

/// `n` points of a random walk: every step moves each coordinate by up to 1,
/// and the first by 0.1 at least, so that no two neighbors are equal.
fn walk<const D: usize>(n: usize) -> Vec<[f64; D]> {
    let mut rng = Lcg(1);
    let mut at = [0.0; D];
    (0..n)
        .map(|_| {
            for coordinate in &mut at {
                *coordinate += 2.0 * rng.next() - 1.0;
            }
            at[0] += 1.1;
            at
        })
        .collect()
}

fn construct<const D: usize>(c: &mut Criterion) {
    let mut group = c.benchmark_group(format!("curve_new_{D}d"));
    let bc = BoundaryCondition::NotAKnot;
    for n in SIZES {
        let points = walk::<D>(n);
        group.throughput(Throughput::Elements(n as u64));
        for (label, way) in [
            ("uniform", Parametrization::Uniform),
            ("chord_length", Parametrization::ChordLength),
            ("centripetal", Parametrization::Centripetal),
        ] {
            group.bench_with_input(BenchmarkId::new(label, n), &points, |b, points| {
                b.iter(|| black_box(BSplineCurve::new(black_box(points), 3, bc, way).is_ok()));
            });
        }
        let Ok(curve) = BSplineCurve::new(&points, 3, bc, Parametrization::ChordLength) else {
            continue;
        };
        let t = curve.parameters().to_vec();
        group.bench_with_input(
            BenchmarkId::new("with_parameters", n),
            &(&t, &points),
            |b, (t, points)| {
                b.iter(|| {
                    let built =
                        BSplineCurve::with_parameters(black_box(t), black_box(points), 3, bc);
                    black_box(built.is_ok())
                });
            },
        );
        group.bench_with_input(
            BenchmarkId::new("bspline1d", n),
            &(&t, &points),
            |b, (t, points)| {
                b.iter(|| {
                    black_box(BSpline1D::new(black_box(t), black_box(points), 3, bc).is_ok())
                });
            },
        );
    }
    group.finish();
}

fn construction(c: &mut Criterion) {
    construct::<2>(c);
    construct::<3>(c);
}

fn evaluation(c: &mut Criterion) {
    let points = walk::<3>(56);
    let bc = BoundaryCondition::NotAKnot;
    let Ok(curve) = BSplineCurve::new(&points, 3, bc, Parametrization::ChordLength) else {
        return;
    };
    let n = 1000;
    let sorted: Vec<f64> = (0..n).map(|i| (i as f64 + 0.5) / n as f64).collect();
    let mut shuffled = sorted.clone();
    let mut rng = Lcg(2);
    for i in (1..n).rev() {
        let j = (rng.next() * (i + 1) as f64) as usize;
        shuffled.swap(i, j.min(i));
    }
    let mut out = vec![[0.0; 3]; n];
    let mut group = c.benchmark_group("curve_3d");
    group.throughput(Throughput::Elements(n as u64));
    for (order, us) in [("sorted", &sorted), ("shuffled", &shuffled)] {
        group.bench_function(format!("eval/{order}"), |b| {
            b.iter(|| {
                for (u, slot) in us.iter().zip(out.iter_mut()) {
                    *slot = curve.eval_or(black_box(u), [f64::NAN; 3]);
                }
                black_box(&out);
            });
        });
        group.bench_function(format!("eval_hinted/{order}"), |b| {
            b.iter(|| {
                let mut hint = Hint::new();
                for (u, slot) in us.iter().zip(out.iter_mut()) {
                    *slot = curve
                        .eval_hinted(black_box(*u), &mut hint)
                        .unwrap_or([f64::NAN; 3]);
                }
                black_box(&out);
            });
        });
        group.bench_function(format!("eval_batch_into/{order}"), |b| {
            b.iter(|| {
                let done = curve.eval_batch_into(black_box(us), &mut out, [f64::NAN; 3]);
                black_box((&out, done.is_ok()));
            });
        });
    }
    group.finish();
}

criterion_group!(benches, construction, evaluation);
criterion_main!(benches);
