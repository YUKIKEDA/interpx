//! Batches evaluated on several threads (the `rayon` feature).

use crate::{InterpolateError, Interpolator};
use rayon::prelude::*;
use std::convert::Infallible;
use std::ops::Range;
use std::time::{Duration, Instant};

/// The points timed first. Few, so that an interpolator that takes
/// microseconds a point is not held on one thread for long: with 64, a batch
/// of 100 points with 20 neighbors (0.5 ms) was 1.2 times as fast as the
/// serial one, with 16 it is 2.0 times (`benches/par.rs`, 16 threads).
const FIRST: usize = 16;

/// How many times as many points are timed when the ones before took too
/// short a time to be told from the step of the clock.
const MORE: usize = 16;

/// The shortest time that is taken for a measurement: twenty steps of a clock
/// that goes by 100 ns.
const TRUSTED: Duration = Duration::from_micros(2);

/// The work that is worth waking the other threads for. Under it, a batch is
/// finished on the thread that called.
///
/// The first points of a batch are slower than the rest (nothing is in the
/// caches yet), so the work ahead is put too high. At 30 µs a grid batch of
/// 1000 points, 26 µs of work in all, was shared out and took 37 µs; at
/// 100 µs it is left alone (`benches/par.rs`).
const WORTH: Duration = Duration::from_micros(100);

/// The work of one chunk.
const CHUNK: Duration = Duration::from_micros(10);

/// The fewest points of a chunk, and what a longer chunk is a multiple of.
const SHORTEST: usize = 16;

/// How the points after the first ones are to be evaluated.
struct Rest {
    /// The first point not evaluated yet.
    from: usize,
    /// The length of the chunks to share the rest out in, or `None` to go on
    /// alone.
    chunk: Option<usize>,
}

/// Evaluates the first of `n` points with `serial`, on this thread, and
/// decides from the time they took what to do with the rest. `serial`
/// evaluates the points of a range and says how long that took ([`timed`]).
///
/// The time a point takes differs a thousandfold between interpolators (5 ns
/// for a piecewise polynomial, 5 µs for a system of neighbors to solve), so
/// no fixed number of points says when sharing a batch out pays, or how long
/// a chunk should be. The first points say it. They are evaluated either
/// way, so timing them wastes nothing but two readings of the clock.
fn lead<E>(
    n: usize,
    mut serial: impl FnMut(Range<usize>) -> Result<Duration, E>,
) -> Result<Rest, E> {
    let mut from = 0;
    let mut timed = FIRST;
    loop {
        let to = n.min(from + timed);
        let took = serial(from..to)?;
        let points = to - from;
        from = to;
        if took >= TRUSTED {
            // With no point left, no work is ahead, and nothing is shared out.
            let per_point = took.as_secs_f64() / points as f64;
            let ahead = per_point * (n - from) as f64;
            let chunk = (ahead >= WORTH.as_secs_f64()).then(|| {
                // Truncating is the intent: a whole number of points.
                let fitting = (CHUNK.as_secs_f64() / per_point) as usize;
                (fitting / SHORTEST).max(1) * SHORTEST
            });
            return Ok(Rest { from, chunk });
        }
        // Too short to measure. With few points left there is little work
        // left either; with many, more of them are timed.
        timed *= MORE;
        if n - from <= timed {
            return Ok(Rest { from, chunk: None });
        }
    }
}

/// Runs `work` and says how long it took.
fn timed<E>(work: impl FnOnce() -> Result<(), E>) -> Result<Duration, E> {
    let start = Instant::now();
    work()?;
    Ok(start.elapsed())
}

/// Batches evaluated on several threads.
///
/// Every interpolator that can be shared between threads has these methods:
/// the types of this crate, and a type of yours that implements
/// [`Interpolator`] and is `Sync`. They are there with the `rayon` feature.
///
/// The results are those of [`Interpolator::eval_batch`] and
/// [`Interpolator::eval_batch_into`], bit for bit, whatever the number of
/// threads: the points are cut into chunks, and each chunk is evaluated by
/// the serial batch, which gives every point the value `eval` gives it.
///
/// # When to use them
///
/// A batch times its own first points and is shared out only when the rest
/// is work worth waking other threads for (a tenth of a millisecond). A batch
/// that is not shared out costs more than the serial form: nothing more up to
/// 16 points, and above that the two readings of the clock and, more than
/// those, the batch being evaluated as two serial batches (the points timed,
/// then the rest). Over the 277 cases of this crate's comparison with SciPy
/// (8 threads), that was up to 0.14 µs on batches under a microsecond (up to
/// two thirds more time), and up to a third more time on batches from a
/// microsecond to a tenth of a millisecond (7 µs at the most). For batches
/// you know to be small, call the serial form.
///
/// The threads are those of `rayon`'s global pool. To use another pool, call
/// from inside its `install`.
///
/// # Examples
///
/// ```
/// use interpx::{InterpolateError, Interpolator, Linear1D, ParInterpolator};
///
/// # fn main() -> Result<(), InterpolateError> {
/// let line = Linear1D::new(&[0.0, 1.0, 3.0], &[0.0, 10.0, 30.0])?;
/// let points: Vec<f64> = (0..10_000).map(|i| f64::from(i) * 3e-4).collect();
///
/// let mut serial = vec![0.0; points.len()];
/// let mut parallel = vec![0.0; points.len()];
/// line.eval_batch_into(&points, &mut serial, f64::NAN)?;
/// line.eval_batch_into_par(&points, &mut parallel, f64::NAN)?;
/// assert_eq!(serial, parallel);
///
/// assert_eq!(line.eval_batch_par(&[0.5, 2.0]), [Some(5.0), Some(20.0)]);
/// # Ok(())
/// # }
/// ```
pub trait ParInterpolator: Interpolator + Sync
where
    Self::Point: Sync,
    Self::Output: Send + Sync + Clone,
{
    /// Evaluates at every point, as [`Interpolator::eval_batch`] does, on
    /// several threads.
    ///
    /// Besides the result, one `Vec` is allocated for each chunk (what the
    /// serial batch of the chunk returns).
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{InterpolateError, Interpolator, Linear1D, ParInterpolator};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let line = Linear1D::new(&[0.0, 1.0], &[0.0, 10.0])?;
    /// let points: Vec<f64> = (0..1000).map(|i| f64::from(i) / 1000.0).collect();
    /// assert_eq!(line.eval_batch_par(&points), line.eval_batch(&points));
    /// # Ok(())
    /// # }
    /// ```
    fn eval_batch_par(&self, points: &[Self::Point]) -> Vec<Option<Self::Output>> {
        if points.len() <= FIRST {
            // Nothing to decide: these would all be evaluated to be timed.
            return self.eval_batch(points);
        }
        let mut out = Vec::with_capacity(points.len());
        let first = lead(points.len(), |range| {
            timed(|| {
                out.extend(self.eval_batch(&points[range]));
                Ok::<(), Infallible>(())
            })
        });
        let Ok(Rest { from, chunk }) = first;
        let rest = &points[from..];
        match chunk {
            Some(chunk) => {
                let chunks = rest.par_chunks(chunk);
                out.par_extend(chunks.flat_map_iter(|chunk| self.eval_batch(chunk)));
            }
            None => out.extend(self.eval_batch(rest)),
        }
        out
    }

    /// Evaluates at every point and writes the results to `out`, as
    /// [`Interpolator::eval_batch_into`] does, on several threads.
    ///
    /// Nothing is allocated here beyond what the serial batch of a chunk
    /// allocates (and what `rayon` does to share the chunks out).
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{InterpolateError, Linear1D, ParInterpolator};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let line = Linear1D::new(&[0.0, 1.0], &[0.0, 10.0])?;
    /// let mut out = [0.0; 3];
    /// line.eval_batch_into_par(&[0.25, 0.5, 2.0], &mut out, -1.0)?;
    /// assert_eq!(out, [2.5, 5.0, 10.0]);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// Returns [`InterpolateError::LengthMismatch`] when `points` and `out`
    /// differ in length, before anything is evaluated.
    fn eval_batch_into_par(
        &self,
        points: &[Self::Point],
        out: &mut [Self::Output],
        fill: Self::Output,
    ) -> Result<(), InterpolateError> {
        if points.len() != out.len() {
            return Err(InterpolateError::LengthMismatch {
                expected: points.len(),
                actual: out.len(),
            });
        }
        if points.len() <= FIRST {
            // Nothing to decide: these would all be evaluated to be timed.
            return self.eval_batch_into(points, out, fill);
        }
        let Rest { from, chunk } = lead(points.len(), |range| {
            timed(|| self.eval_batch_into(&points[range.clone()], &mut out[range], fill.clone()))
        })?;
        let (points, out) = (&points[from..], &mut out[from..]);
        match chunk {
            Some(chunk) => points
                .par_chunks(chunk)
                .zip(out.par_chunks_mut(chunk))
                .try_for_each(|(points, out)| self.eval_batch_into(points, out, fill.clone())),
            None => self.eval_batch_into(points, out, fill),
        }
    }
}

impl<I> ParInterpolator for I
where
    I: Interpolator + Sync + ?Sized,
    I::Point: Sync,
    I::Output: Send + Sync + Clone,
{
}

#[cfg(test)]
mod tests {
    use super::{lead, Rest, CHUNK, FIRST, MORE, SHORTEST, WORTH};
    use std::time::Duration;

    /// Goes through `lead` with a serial evaluation that says it took
    /// `per_point` for each point (no clock is read: the decision is tested,
    /// not the machine).
    fn decide(n: usize, per_point: Duration) -> (Rest, Vec<(usize, usize)>) {
        let mut asked = Vec::new();
        let rest = lead(n, |range| {
            asked.push((range.start, range.end));
            Ok::<_, ()>(per_point * u32::try_from(range.len()).unwrap_or(u32::MAX))
        });
        (
            rest.unwrap_or(Rest {
                from: 0,
                chunk: None,
            }),
            asked,
        )
    }

    #[test]
    fn the_first_points_are_evaluated_once_and_in_order() {
        for n in [0, 1, FIRST - 1, FIRST, FIRST + 1, 3 * FIRST] {
            let (rest, asked) = decide(n, Duration::ZERO);
            // Points that take no time: more and more are timed, to the end.
            let covered: usize = asked.iter().map(|(start, end)| end - start).sum();
            assert_eq!(covered, rest.from, "{n}");
            assert!(asked.windows(2).all(|pair| pair[0].1 == pair[1].0), "{n}");
            assert_eq!(asked.first().map_or(0, |first| first.0), 0, "{n}");
            assert!(rest.from <= n);
            assert_eq!(rest.chunk, None, "{n}: nothing worth sharing out");
        }
    }

    #[test]
    fn few_points_that_take_no_time_are_left_to_the_caller() {
        // More than the first points, and no more than would be timed next.
        let n = FIRST + FIRST * MORE;
        let (rest, asked) = decide(n, Duration::ZERO);
        assert_eq!(asked, [(0, FIRST)]);
        assert_eq!((rest.from, rest.chunk), (FIRST, None));
        // One more, and the next ones are timed as well.
        let (rest, asked) = decide(n + 1, Duration::ZERO);
        assert_eq!(asked, [(0, FIRST), (FIRST, n)]);
        assert_eq!((rest.from, rest.chunk), (n, None));
    }

    #[test]
    fn work_worth_sharing_is_shared_in_chunks_of_about_the_work_of_one() {
        // 50 µs a point: after the first points, each point alone is more
        // than the work of a chunk, so a chunk is as short as chunks get.
        let slow = Duration::from_micros(50);
        let (rest, asked) = decide(FIRST + 10, slow);
        assert_eq!(asked, [(0, FIRST)]);
        assert_eq!((rest.from, rest.chunk), (FIRST, Some(SHORTEST)));
        assert!(slow * 10 >= WORTH && slow > CHUNK);
        // The last points evaluated leave nothing to decide.
        let (rest, _) = decide(FIRST, slow);
        assert_eq!((rest.from, rest.chunk), (FIRST, None));
        // 200 ns a point (the first points take long enough to be believed):
        // a chunk is the points that make up its work, in whole sixteens; and
        // little work ahead is not shared out.
        let quick = Duration::from_nanos(200);
        let (rest, asked) = decide(100_000, quick);
        assert_eq!(asked, [(0, FIRST)]);
        let fitting = CHUNK.as_nanos() / quick.as_nanos();
        assert_eq!(
            rest.chunk,
            usize::try_from(fitting)
                .ok()
                .map(|c| c / SHORTEST * SHORTEST)
        );
        let few = FIRST + usize::try_from(WORTH.as_nanos() / quick.as_nanos()).unwrap_or(0);
        // (One point either side of the line: the division rounds.)
        assert_eq!(decide(few - 1, quick).0.chunk, None);
        assert!(decide(few + 1, quick).0.chunk.is_some());
    }

    #[test]
    fn an_error_of_the_first_points_is_the_error() {
        let failed = lead(1000, |_| Err::<Duration, u8>(7));
        assert_eq!(failed.err(), Some(7));
    }
}
