//! `ParInterpolator`: the parallel batches against the serial ones, bit for
//! bit, for every type, with one thread and with several.

#![cfg(feature = "rayon")]

use interpx::{
    Akima1D, AkimaMethod, BSpline1D, BSplineCurve, BoundaryCondition, CubicHermiteSpline,
    CubicSpline, Extrapolate, GridMethod, InterpolateError as E, Interpolator, Linear1D, Nearest1D,
    NearestKind, ParInterpolator, Parametrization, PchipInterpolator, RegularGridInterpolator,
};
use std::fmt::Debug;

/// A small generator, so that the inputs are the same on every run.
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

/// Runs `work` in a pool of `threads` threads.
// Marked as test code so that `expect` is allowed here as in a `#[test]`
// function (clippy.toml allows it in tests only).
#[cfg(test)]
fn in_pool<R: Send>(threads: usize, work: impl FnOnce() -> R + Send) -> R {
    let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build();
    pool.expect("a pool").install(work)
}

/// Holds the two parallel batches against the serial ones, for a short batch
/// (never shared out), a long one, and lengths around a chunk; with one
/// thread and with five. The values are compared as `Debug` text of their
/// bits would be: through `PartialEq` on `Option`, after NaN is ruled out by
/// the fill and the `None`.
fn check<I>(name: &str, interp: &I, points: &[I::Point], fill: I::Output)
where
    I: Interpolator + Sync,
    I::Point: Sync,
    I::Output: Send + Sync + Clone + PartialEq + Debug,
{
    for n in [0, 1, 17, 64, 65, 1000, points.len()] {
        let points = &points[..n.min(points.len())];
        let serial = interp.eval_batch(points);
        let mut serial_into = vec![fill.clone(); points.len()];
        let done = interp.eval_batch_into(points, &mut serial_into, fill.clone());
        assert_eq!(done, Ok(()), "{name}");
        for threads in [1, 5] {
            let (options, into) = in_pool(threads, || {
                let mut into = vec![fill.clone(); points.len()];
                let done = interp.eval_batch_into_par(points, &mut into, fill.clone());
                assert_eq!(done, Ok(()), "{name}");
                (interp.eval_batch_par(points), into)
            });
            assert_eq!(
                options, serial,
                "{name}: eval_batch_par, {n} points, {threads} threads"
            );
            assert_eq!(
                into, serial_into,
                "{name}: eval_batch_into_par, {n} points, {threads} threads"
            );
        }
    }
    // The lengths are held against each other before anything is evaluated.
    let mut short = vec![fill.clone(); points.len() - 1];
    let mismatch = E::LengthMismatch {
        expected: points.len(),
        actual: points.len() - 1,
    };
    assert_eq!(
        interp.eval_batch_into_par(points, &mut short, fill.clone()),
        Err(mismatch),
        "{name}"
    );
    assert_eq!(
        interp.eval_batch_into(points, &mut short, fill),
        Err(mismatch),
        "{name}"
    );
}

/// Numbers with a bit pattern of their own as a fill, so that a filled slot
/// is not taken for a value.
const FILL: f64 = -12345.678;

#[test]
fn every_type_of_one_variable_gives_the_serial_results() {
    let mut rng = Lcg(11);
    let x: Vec<f64> = (0..200)
        .scan(0.0, |at, _| {
            Some(std::mem::replace(at, *at + 0.5 + rng.next()))
        })
        .collect();
    let y: Vec<f64> = x.iter().map(|x| (0.3 * x).sin()).collect();
    let slopes: Vec<f64> = x.iter().map(|x| 0.3 * (0.3 * x).cos()).collect();
    let (low, high) = (x[0], x[x.len() - 1]);
    // Inside, outside both ends (no value where the type gives none), in no
    // order; then the same in order, where a batch takes its hints.
    let mut points: Vec<f64> = (0..20_000)
        .map(|_| low - 5.0 + rng.next() * (high - low + 10.0))
        .collect();
    let mut sorted = points.clone();
    sorted.sort_by(f64::total_cmp);
    points.extend(sorted);
    points.push(f64::NAN);

    check(
        "Linear1D",
        &Linear1D::new(&x, &y)
            .expect("valid")
            .with_extrapolate(Extrapolate::Forbid),
        &points,
        FILL,
    );
    check(
        "Nearest1D",
        &Nearest1D::new(&x, &y, NearestKind::Nearest).expect("valid"),
        &points,
        FILL,
    );
    check(
        "CubicSpline",
        &CubicSpline::new(&x, &y, BoundaryCondition::NotAKnot).expect("valid"),
        &points,
        FILL,
    );
    let periodic: Vec<f64> = y
        .iter()
        .enumerate()
        .map(|(i, v)| if i + 1 == y.len() { y[0] } else { *v })
        .collect();
    check(
        "CubicSpline, periodic",
        &CubicSpline::new(&x, &periodic, BoundaryCondition::Periodic).expect("valid"),
        &points,
        FILL,
    );
    check(
        "PchipInterpolator",
        &PchipInterpolator::new(&x, &y)
            .expect("valid")
            .with_extrapolate(Extrapolate::Forbid),
        &points,
        FILL,
    );
    check(
        "Akima1D",
        &Akima1D::new(&x, &y, AkimaMethod::Makima).expect("valid"),
        &points,
        FILL,
    );
    check(
        "CubicHermiteSpline",
        &CubicHermiteSpline::new(&x, &y, &slopes).expect("valid"),
        &points,
        FILL,
    );
    let spline = BSpline1D::new(&x, &y, 5, BoundaryCondition::NotAKnot).expect("valid");
    check(
        "BSpline1D",
        &spline.clone().with_extrapolate(Extrapolate::Forbid),
        &points,
        FILL,
    );
    check(
        "BSpline",
        &spline
            .as_bspline()
            .clone()
            .with_extrapolate(Extrapolate::Clamp),
        &points,
        FILL,
    );
    let cubic = CubicSpline::new(&x, &y, BoundaryCondition::NotAKnot).expect("valid");
    check(
        "PPoly",
        &cubic
            .as_ppoly()
            .clone()
            .with_extrapolate(Extrapolate::Forbid),
        &points,
        FILL,
    );

    // Values of two columns, and a curve.
    let pairs: Vec<[f64; 2]> = y.iter().map(|v| [*v, 1.0 - v]).collect();
    check(
        "Linear1D, two columns",
        &Linear1D::new(&x, &pairs)
            .expect("valid")
            .with_extrapolate(Extrapolate::Forbid),
        &points,
        [FILL; 2],
    );
    let curve = BSplineCurve::new(
        &pairs,
        3,
        BoundaryCondition::NotAKnot,
        Parametrization::ChordLength,
    )
    .expect("valid");
    let (start, end) = curve.domain();
    let along: Vec<f64> = (0..5000)
        .map(|_| start - 0.1 + rng.next() * (end - start + 0.2))
        .collect();
    check(
        "BSplineCurve",
        &curve.with_extrapolate(Extrapolate::Forbid),
        &along,
        [FILL; 2],
    );
}

#[test]
fn grids_give_the_serial_results() {
    let mut rng = Lcg(12);
    let mut axis = |n: usize| -> Vec<f64> {
        (0..n)
            .scan(0.0, |at, _| {
                Some(std::mem::replace(at, *at + 0.5 + rng.next()))
            })
            .collect()
    };
    let axes = [axis(9), axis(7), axis(8)];
    let values: Vec<f64> = (0..9 * 7 * 8_u32)
        .map(|n| f64::from(n * 37 % 101) / 17.0)
        .collect();
    let mut rng = Lcg(13);
    // Inside and a little outside; then a path, where a batch takes hints.
    let mut points: Vec<[f64; 3]> = (0..12_000)
        .map(|_| core::array::from_fn(|d| -0.3 + rng.next() * (axes[d][axes[d].len() - 1] + 0.6)))
        .collect();
    points.extend((0..8000).map(|i| {
        let share = f64::from(i) / 8000.0;
        core::array::from_fn(|d| share * axes[d][axes[d].len() - 1])
    }));
    points.push([f64::NAN, 1.0, 1.0]);
    let slices: [&[f64]; 3] = core::array::from_fn(|d| axes[d].as_slice());
    let methods = [
        GridMethod::Nearest,
        GridMethod::Linear,
        GridMethod::Cubic,
        GridMethod::Quintic,
        GridMethod::Pchip,
    ];
    for method in methods {
        let grid = RegularGridInterpolator::new(slices, &values, method).expect("valid");
        check(&format!("grid, {method:?}"), &grid, &points, FILL);
        let outside = grid.with_extrapolate(Extrapolate::Extend);
        check(
            &format!("grid, {method:?}, continued"),
            &outside,
            &points[..3000],
            FILL,
        );
    }
}

#[cfg(feature = "scattered")]
#[test]
fn scattered_types_give_the_serial_results() {
    use interpx::{
        CloughTocher2DInterpolator, LinearNDInterpolator, NearestNDInterpolator, RBFInterpolator,
        RbfKernel, RbfOptions,
    };

    let mut rng = Lcg(14);
    let data: Vec<[f64; 2]> = (0..300).map(|_| [rng.next(), rng.next()]).collect();
    let values: Vec<f64> = data.iter().map(|p| (3.0 * p[0]).sin() * p[1]).collect();
    let mut points: Vec<[f64; 2]> = (0..3000).map(|_| [rng.next(), rng.next()]).collect();
    // The pixels of a row: with neighbors, a batch keeps a factored system
    // from one point to the next.
    points.extend((0..2000).map(|i| [f64::from(i) / 2000.0, 0.4]));
    points.push([f64::INFINITY, 0.5]);

    let nearest = NearestNDInterpolator::new(&data, &values).expect("valid");
    check(
        "NearestND",
        &nearest.clone().with_distance_bound(0.05),
        &points,
        FILL,
    );
    check("NearestND, no bound", &nearest, &points, FILL);
    let whole = RBFInterpolator::new(&data, &values, RbfKernel::ThinPlateSpline).expect("valid");
    check("RBF, all points", &whole, &points, FILL);
    let options = RbfOptions::new().with_neighbors(12);
    let near = RBFInterpolator::with_options(&data, &values, RbfKernel::Cubic, &options);
    check("RBF, neighbors", &near.expect("valid"), &points, FILL);
    let linear = LinearNDInterpolator::new(&data, &values).expect("valid");
    check("LinearND", &linear, &points, FILL);
    let cubic = CloughTocher2DInterpolator::new(&data, &values).expect("valid");
    check("CloughTocher2D", &cubic, &points, FILL);
}

/// An interpolator of the caller's own, slow enough that the first points
/// say the rest is worth sharing out whatever the machine: every batch of
/// more than the first points goes through the chunks.
struct Slow;

impl Interpolator for Slow {
    type Point = u32;
    type Output = u64;

    fn eval(&self, point: &u32) -> Option<u64> {
        // About 20 µs of work that the compiler cannot take away.
        let mut state = u64::from(*point) + 1;
        for _ in 0..20_000 {
            state = std::hint::black_box(state.wrapping_mul(6_364_136_223_846_793_005) >> 1 | 1);
        }
        (*point % 7 != 3).then_some(state)
    }
}

#[test]
fn an_interpolator_of_the_callers_own_has_the_parallel_batches() {
    let points: Vec<u32> = (0..400).collect();
    check("Slow", &Slow, &points, 0);
    // No value where `eval` gives none, in the same places.
    let options = in_pool(4, || Slow.eval_batch_par(&points));
    assert_eq!(options.iter().filter(|v| v.is_none()).count(), 57);
    // As a trait object as well.
    let object: &(dyn Interpolator<Point = u32, Output = u64> + Sync) = &Slow;
    assert_eq!(
        object.eval_batch_par(&points[..70]),
        Slow.eval_batch(&points[..70])
    );
}

#[test]
fn the_work_is_done_on_the_threads_of_the_pool() {
    use std::collections::HashSet;
    use std::sync::Mutex;

    /// Notes the thread every evaluation runs on.
    struct Noting(Mutex<HashSet<std::thread::ThreadId>>);

    impl Interpolator for Noting {
        type Point = u32;
        type Output = u64;

        fn eval(&self, point: &u32) -> Option<u64> {
            if let Ok(mut seen) = self.0.lock() {
                seen.insert(std::thread::current().id());
            }
            Slow.eval(point)
        }
    }

    let points: Vec<u32> = (0..2000).collect();
    let threads = |pool: usize| {
        let noting = Noting(Mutex::new(HashSet::new()));
        in_pool(pool, || noting.eval_batch_par(&points));
        noting.0.lock().map(|seen| seen.len()).unwrap_or(0)
    };
    assert_eq!(threads(1), 1);
    // 40 ms of work in chunks of 16 points: more than one of four threads
    // takes some of it.
    assert!(threads(4) > 1);
}
