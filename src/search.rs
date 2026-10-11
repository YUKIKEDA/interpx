//! Finding the interval of a query on a grid of break points or knots, shared by
//! [`PPoly`](crate::PPoly) and [`BSpline`](crate::BSpline).
//!
//! A grid here is a slice in non-decreasing order with at least two entries,
//! whose first and last intervals have positive width. Interval `i` runs from
//! `grid[i]` to `grid[i + 1]`: it contains its left end and not its right end,
//! and the last interval contains both. An interval of zero width (two equal
//! knots) contains no point.

use crate::{Extrapolate, InterpolateError, Scalar};

/// The interval the last evaluation fell into, kept by the caller to speed up
/// the next one.
///
/// Queries that come in order mostly hit the same interval or the next. A hint
/// lets [`PPoly::eval_hinted`](crate::PPoly::eval_hinted) and
/// [`BSpline::eval_hinted`](crate::BSpline::eval_hinted) look there first and
/// skip the binary search. The result never depends on the hint, only the time
/// does.
///
/// # Examples
///
/// ```
/// use interpx::{Hint, InterpolateError, Interpolator, PPoly};
///
/// # fn main() -> Result<(), InterpolateError> {
/// // x on [0, 1], 1 + 2(x - 1) on [1, 2].
/// let p = PPoly::new(&[0.0, 1.0, 2.0], &[1.0, 0.0, 2.0, 1.0], 1)?;
/// let mut hint = Hint::new();
/// for x in [0.25, 0.5, 1.5] {
///     assert_eq!(p.eval_hinted(x, &mut hint), p.eval(&x));
/// }
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Hint {
    interval: usize,
}

impl Hint {
    /// Returns a hint that points at the first interval.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::Hint;
    ///
    /// assert_eq!(Hint::new(), Hint::default());
    /// ```
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }
}

/// The number of points `eval_batch_into` decides on at a time: long enough that
/// the two extra searches of the decision are a small share, short enough that
/// a batch which is in order only in places still uses its hint there.
const BATCH_RUN: usize = 64;

/// The shortest run worth the decision. On a run shorter than this the decision
/// costs more than a hint can save (about 20 ns on a 10-point batch, measured in
/// `benches/ppoly.rs`), so the points are searched plainly.
const MIN_HINT_RUN: usize = 16;

/// The index of the last interval of `grid`.
fn last<T>(grid: &[T]) -> usize {
    // Two entries at least, by the contract of this module.
    grid.len() - 2
}

/// The interval of `x`, for `x` between the first and the last entry of `grid`.
fn locate<T: Scalar>(grid: &[T], x: T) -> usize {
    let after = grid.partition_point(|b| *b <= x);
    after.saturating_sub(1).min(last(grid))
}

/// Whether interval `i` contains `x`.
fn contains<T: Scalar>(grid: &[T], i: usize, x: T) -> bool {
    let (left, right) = (grid[i], grid[i + 1]);
    left <= x && (x < right || (i == last(grid) && x <= right))
}

/// The same interval as [`locate`], trying `start` and its two neighbors before
/// the binary search.
fn locate_near<T: Scalar>(grid: &[T], x: T, start: usize) -> usize {
    let last = last(grid);
    let i = start.min(last);
    if contains(grid, i, x) {
        i
    } else if i < last && contains(grid, i + 1, x) {
        i + 1
    } else if i > 0 && contains(grid, i - 1, x) {
        i - 1
    } else {
        locate(grid, x)
    }
}

/// Whether a hint pays for this run of points: they go one way, and on average
/// no more than one interval apart.
///
/// A hinted search starts from the interval of the point before, so each point
/// waits for the one before it. For points that are close together the wait is
/// short and the search is skipped. For points in no order the wait is a whole
/// binary search, and plain searches, which do not wait for each other, are
/// about twice as fast (measured in `benches/ppoly.rs`).
fn hint_pays<T: Scalar>(grid: &[T], points: &[T]) -> bool {
    // A NaN fails every comparison, so a run with one is not in order.
    let ascending = points.windows(2).all(|w| w[0] <= w[1]);
    let ordered = ascending || points.windows(2).all(|w| w[0] >= w[1]);
    // An empty run has no ends and nothing to speed up.
    ordered
        && points
            .first()
            .zip(points.last())
            .is_some_and(|(first, last)| {
                let (low, high) = if ascending {
                    (*first, *last)
                } else {
                    (*last, *first)
                };
                let crossed =
                    grid.partition_point(|b| *b <= high) - grid.partition_point(|b| *b <= low);
                crossed <= points.len()
            })
}

/// Maps a query to the interval to evaluate and the point to evaluate it at, or
/// `None` when the query has no value.
fn resolve<T: Scalar>(
    grid: &[T],
    extrapolate: Extrapolate,
    x: T,
    locate: impl FnOnce(T) -> usize,
) -> Option<(usize, T)> {
    if x.is_nan() {
        return None;
    }
    let last = last(grid);
    let (lo, hi) = (grid[0], grid[last + 1]);
    if lo <= x && x <= hi {
        return Some((locate(x), x));
    }
    match extrapolate {
        Extrapolate::Forbid => None,
        Extrapolate::Clamp => Some(if x < lo { (0, lo) } else { (last, hi) }),
        Extrapolate::Extend => Some((if x < lo { 0 } else { last }, x)),
        Extrapolate::Periodic => {
            let period = hi - lo;
            let mut offset = (x - lo) % period;
            if offset < T::zero() {
                offset = offset + period;
            }
            // An infinite query has no place in the period.
            if !offset.is_finite() {
                return None;
            }
            // `lo + offset` can round up to just past `hi`.
            let wrapped = (lo + offset).min(hi);
            Some((locate(wrapped), wrapped))
        }
    }
}

/// Evaluates at `x`: finds the interval by binary search and calls `at` with the
/// interval and the point to evaluate it at.
pub(crate) fn eval<T: Scalar, V>(
    grid: &[T],
    extrapolate: Extrapolate,
    x: T,
    at: impl FnOnce(usize, T) -> V,
) -> Option<V> {
    let (interval, x) = resolve(grid, extrapolate, x, |x| locate(grid, x))?;
    Some(at(interval, x))
}

/// The same as [`eval`], looking first at the interval in `hint` and its
/// neighbors, and leaving the interval used in `hint`.
pub(crate) fn eval_hinted<T: Scalar, V>(
    grid: &[T],
    extrapolate: Extrapolate,
    x: T,
    hint: &mut Hint,
    at: impl FnOnce(usize, T) -> V,
) -> Option<V> {
    let start = hint.interval;
    let (interval, x) = resolve(grid, extrapolate, x, |x| locate_near(grid, x, start))?;
    hint.interval = interval;
    Some(at(interval, x))
}

/// Evaluates at every point and writes the results to `out`, with `fill` where
/// a point has no value.
///
/// The points are taken by runs of [`BATCH_RUN`]: with a hint where it pays,
/// plainly where it does not. Either way each point gets the interval [`eval`]
/// finds, so the results are bit for bit those of a loop over [`eval`].
pub(crate) fn eval_batch_into<T: Scalar, V: Copy>(
    grid: &[T],
    extrapolate: Extrapolate,
    points: &[T],
    out: &mut [V],
    fill: V,
    at: impl Fn(usize, T) -> V,
) -> Result<(), InterpolateError> {
    if points.len() != out.len() {
        return Err(InterpolateError::LengthMismatch {
            expected: points.len(),
            actual: out.len(),
        });
    }
    let mut hint = Hint::new();
    for (xs, slots) in points.chunks(BATCH_RUN).zip(out.chunks_mut(BATCH_RUN)) {
        if xs.len() >= MIN_HINT_RUN && hint_pays(grid, xs) {
            for (x, slot) in xs.iter().zip(slots) {
                *slot = eval_hinted(grid, extrapolate, *x, &mut hint, &at).unwrap_or(fill);
            }
        } else {
            for (x, slot) in xs.iter().zip(slots) {
                *slot = eval(grid, extrapolate, *x, &at).unwrap_or(fill);
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{contains, locate, locate_near};

    #[test]
    fn an_interval_of_zero_width_holds_no_point() {
        // Knots with a repeated inner value: intervals 1 and 2 are empty.
        let grid = [0.0, 1.0, 1.0, 1.0, 2.0];
        for x in [0.0, 0.5, 1.0, 1.5, 2.0] {
            assert!(!contains(&grid, 1, x), "x = {x}");
            assert!(!contains(&grid, 2, x), "x = {x}");
        }
        // The repeated value belongs to the interval on its right.
        assert_eq!(locate(&grid, 1.0), 3);
        assert_eq!(locate(&grid, 0.999), 0);
        assert_eq!(locate(&grid, 2.0), 3);
        for start in 0..6 {
            for x in [0.0, 0.5, 1.0, 1.5, 2.0] {
                assert_eq!(
                    locate_near(&grid, x, start),
                    locate(&grid, x),
                    "{start} {x}"
                );
            }
        }
    }
}
