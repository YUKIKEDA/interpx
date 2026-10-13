//! Finding the interval of a query on a grid of break points or knots, shared by
//! [`PPoly`](crate::PPoly) and [`BSpline`](crate::BSpline).
//!
//! A grid here is a slice in non-decreasing order with at least two entries,
//! whose first and last intervals have positive width. Interval `i` runs from
//! `grid[i]` to `grid[i + 1]`: it contains its left end and not its right end,
//! and the last interval contains both. An interval of zero width (two equal
//! knots) contains no point.

use crate::{Extrapolate, InterpolateError, Scalar, Value};

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
pub(crate) const BATCH_RUN: usize = 64;

/// The shortest run worth the decision. On a run shorter than this the decision
/// costs more than a hint can save (about 20 ns on a 10-point batch, measured in
/// `benches/ppoly.rs`), so the points are searched plainly.
pub(crate) const MIN_HINT_RUN: usize = 16;

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

/// [`locate`] for a grid that is [`evenly_spaced`]: the interval `x` would be
/// in were the spacing exactly even is looked at first. When it holds `x` it
/// is the answer, in one step and whatever the rounding of the guess, because
/// one interval alone holds a point. When it does not, the entries are
/// searched by halves.
///
/// On 10 000 evenly spaced points this is 3.6 times as fast for points in no
/// order (16 ns a point to 4.5 ns; `docs/performance.md`, #123). On a grid
/// that is not evenly spaced the guess is a division and a look that find
/// nothing: tried on every grid, it made the evaluation on grids of 10 to 40
/// points an axis up to twice as slow. So the caller says which grid it has.
fn locate_guessing<T: Scalar>(grid: &[T], x: T) -> usize {
    let last = last(grid);
    let (low, high) = (grid[0], grid[last + 1]);
    // No guess (an `x` that is not a number) is a miss.
    let guess = T::from(last + 1)
        .and_then(|intervals| ((x - low) / (high - low) * intervals).to_usize())
        .map_or(0, |guess| guess.min(last));
    if contains(grid, guess, x) {
        guess
    } else {
        locate(grid, x)
    }
}

/// Whether the entries of `grid` are evenly spaced closely enough for
/// [`locate_guessing`] to pay: every spacing within a sixteenth of the mean
/// spacing. Found once, at construction. It decides only how the search
/// starts, never what it finds.
///
/// The spacings are compared, not the places of the entries: one subtraction
/// and one comparison an entry, several entries at once. (Each entry against
/// the place it would have, which takes its index as a number, cost 2 ns an
/// entry: as much as the rest of the construction of `Linear1D`.) Spacings
/// that are all a little long and then all a little short pass, and the
/// guess then misses more often than it would; it is still only a guess.
pub(crate) fn evenly_spaced<T: Scalar>(grid: &[T]) -> bool {
    let intervals = grid.len() - 1;
    let numbers = T::from(intervals).zip(T::from(16));
    numbers.is_some_and(|(count, sixteen)| {
        let step = (grid[intervals] - grid[0]) / count;
        let slack = step / sixteen;
        // A grid of no width has no step, and no spacing is within no slack
        // of it but its own, which are all zero: it is refused by name.
        let near = |pair: &[T]| ((pair[1] - pair[0]) - step).abs() <= slack;
        step > T::zero() && grid.windows(2).fold(true, |all, pair| all & near(pair))
    })
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
    hint_pays_along(grid, points.iter().copied())
}

/// [`hint_pays`] for coordinates that are not a slice of their own: one
/// coordinate of each point of a run of points in several dimensions.
pub(crate) fn hint_pays_along<T: Scalar>(
    grid: &[T],
    coordinates: impl ExactSizeIterator<Item = T> + Clone,
) -> bool {
    let pairs = || coordinates.clone().zip(coordinates.clone().skip(1));
    // A NaN fails every comparison, so a run with one is not in order.
    let ascending = pairs().all(|(a, b)| a <= b);
    let ordered = ascending || pairs().all(|(a, b)| a >= b);
    // An empty run has no ends and nothing to speed up.
    ordered
        && coordinates
            .clone()
            .next()
            .zip(coordinates.clone().last())
            .is_some_and(|(first, last)| {
                let (low, high) = if ascending {
                    (first, last)
                } else {
                    (last, first)
                };
                let crossed =
                    grid.partition_point(|b| *b <= high) - grid.partition_point(|b| *b <= low);
                crossed <= coordinates.len()
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
    // Every query inside the grid passes through this test, so it is kept to
    // two comparisons: a NaN fails them, and the right end is left to
    // `resolve_at_or_past_an_end`. (With the mode tested here as well, a hinted
    // evaluation took about a tenth longer, measured in `benches/ppoly.rs`.)
    if grid[0] <= x && x < grid[last(grid) + 1] {
        return Some((locate(x), x));
    }
    resolve_at_or_past_an_end(grid, extrapolate, x, locate)
}

/// [`resolve`] for a query that is not inside the grid: at its right end,
/// outside, or NaN.
fn resolve_at_or_past_an_end<T: Scalar>(
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
    // The right end of the grid belongs to the last interval, except that a
    // period holds its left end and not its right one: there the right end is
    // the left end again, as in SciPy.
    if lo <= x && x <= hi && extrapolate != Extrapolate::Periodic {
        return Some((last, x));
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

/// Evaluates at `x`: finds the interval and calls `at` with the interval and
/// the point to evaluate it at. `even` is what [`evenly_spaced`] found for
/// `grid`: the search then starts from a guess.
pub(crate) fn eval<T: Scalar, V>(
    grid: &[T],
    even: bool,
    extrapolate: Extrapolate,
    x: T,
    at: impl FnOnce(usize, T) -> V,
) -> Option<V> {
    let find = |x| {
        if even {
            locate_guessing(grid, x)
        } else {
            locate(grid, x)
        }
    };
    let (interval, x) = resolve(grid, extrapolate, x, find)?;
    Some(at(interval, x))
}

/// Evaluates a derivative of order one or more at `x`.
///
/// The same as [`eval`], except under [`Extrapolate::Clamp`] outside the grid:
/// there the function is constant, so the derivative is `zero`.
pub(crate) fn eval_derivative<T: Scalar, V>(
    grid: &[T],
    extrapolate: Extrapolate,
    x: T,
    zero: V,
    at: impl FnOnce(usize, T) -> V,
) -> Option<V> {
    // A NaN is on neither side, and `eval` turns it away.
    if extrapolate == Extrapolate::Clamp && (x < grid[0] || x > grid[last(grid) + 1]) {
        return Some(zero);
    }
    eval(grid, false, extrapolate, x, at)
}

/// The integral from `a` to `b` of the function [`eval`] defines, or `None` when
/// a limit is not finite or the function has no value somewhere between them.
///
/// `at(i, x)` is the value of the piece of interval `i` at `x`. `between(i, p, j,
/// q)` is the integral from `p` to `q`, `p <= q`, where `p` is taken to be in
/// interval `i` and `q` in interval `j`, each on its piece continued when it is
/// outside the grid.
pub(crate) fn integrate<T: Scalar, V: Value<T>>(
    grid: &[T],
    extrapolate: Extrapolate,
    a: T,
    b: T,
    at: impl Fn(usize, T) -> V,
    between: impl Fn(usize, T, usize, T) -> V,
) -> Option<V> {
    if !a.is_finite() || !b.is_finite() {
        return None;
    }
    let reversed = a > b;
    let (a, b) = if reversed { (b, a) } else { (a, b) };
    let last = last(grid);
    let (lo, hi) = (grid[0], grid[last + 1]);
    // `locate` gives the first interval to a point below the grid and the last
    // one to a point above it.
    let pieces = |p: T, q: T| between(locate(grid, p), p, locate(grid, q), q);
    let inside = lo <= a && b <= hi;
    let total = match extrapolate {
        // Inside the grid the mode does not matter.
        _ if inside => pieces(a, b),
        Extrapolate::Forbid => return None,
        Extrapolate::Extend => pieces(a, b),
        Extrapolate::Clamp => {
            let mut sum = V::zero();
            // Outside, the function is the value at the end.
            if a < lo {
                sum = sum.axpy(b.min(lo) - a, at(0, lo));
            }
            if b > hi {
                sum = sum.axpy(b - a.max(hi), at(last, hi));
            }
            let (p, q) = (a.max(lo), b.min(hi));
            if p < q {
                sum = sum.axpy(T::one(), pieces(p, q));
            }
            sum
        }
        Extrapolate::Periodic => {
            let period = hi - lo;
            // Whole periods first, then what is left, from the place of `a` in
            // the period.
            let periods = ((b - a) / period).floor();
            let rest = ((b - a) - periods * period).max(T::zero());
            let mut sum = V::zero();
            if periods > T::zero() {
                sum = sum.axpy(periods, pieces(lo, hi));
            }
            let mut offset = (a - lo) % period;
            if offset < T::zero() {
                offset = offset + period;
            }
            // `lo + offset` can round up to just past `hi`.
            let start = (lo + offset).min(hi);
            let end = start + rest;
            if end <= hi {
                sum.axpy(T::one(), pieces(start, end))
            } else {
                // Over the right end: up to it, then on from the left end.
                sum.axpy(T::one(), pieces(start, hi))
                    .axpy(T::one(), pieces(lo, lo + (end - hi)))
            }
        }
    };
    Some(if reversed {
        V::zero().axpy(-T::one(), total)
    } else {
        total
    })
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
    even: bool,
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
                *slot = eval(grid, even, extrapolate, *x, &at).unwrap_or(fill);
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{contains, evenly_spaced, locate, locate_guessing, locate_near};

    /// The interval by halving alone: what `locate` was before it guessed.
    fn by_halves(grid: &[f64], x: f64) -> usize {
        let after = grid.partition_point(|b| *b <= x);
        after.saturating_sub(1).min(grid.len() - 2)
    }

    #[test]
    fn the_guess_changes_no_interval() {
        // Even spacing (the guess is right, up to its rounding), spacing that
        // grows (the guess is wrong nearly everywhere), entries that repeat,
        // two entries, and a grid of no width.
        let even: Vec<f64> = (0..97).map(|i| 0.1 * f64::from(i) - 3.0).collect();
        let growing: Vec<f64> = (0..60).map(|i| 1.07_f64.powi(i)).collect();
        let stepped: Vec<f64> = (0..40).map(|i| f64::from(i / 3)).collect();
        let grids = [even, growing, stepped, vec![2.0, 5.0], vec![1.0, 1.0]];
        for grid in &grids {
            let (low, high) = (grid[0], grid[grid.len() - 1]);
            let mut queries = vec![low, high];
            for pair in grid.windows(2) {
                // Each entry, the numbers next to it on both sides, and a
                // point inside the interval.
                queries.extend([pair[0].next_down(), pair[0], pair[0].next_up()]);
                queries.push(pair[0] + 0.37 * (pair[1] - pair[0]));
            }
            for x in queries.into_iter().filter(|x| (low..=high).contains(x)) {
                let found = locate_guessing(grid, x);
                assert_eq!(found, by_halves(grid, x), "{x} in {grid:?}");
                assert_eq!(locate(grid, x), found);
                // The interval holds the point, unless it is the repeated
                // entry of a grid of no width.
                assert!(contains(grid, found, x), "{x} in {grid:?}");
            }
        }
        // A number that is not one has no guess, and the search by halves
        // answers as it did.
        let nan = f64::NAN;
        assert_eq!(locate_guessing(&grids[0], nan), by_halves(&grids[0], nan));
        // Single precision goes the same way.
        let single: Vec<f32> = (0..33_u8).map(|i| 0.3 * f32::from(i)).collect();
        for (i, entry) in single.iter().enumerate() {
            assert_eq!(locate_guessing(&single, *entry), i.min(31));
            let below = entry.next_down().max(0.0);
            assert_eq!(locate_guessing(&single, below), i.saturating_sub(1));
        }
    }

    #[test]
    fn a_grid_is_evenly_spaced_when_every_entry_is_near_its_even_place() {
        let even: Vec<f64> = (0..97).map(|i| 0.1 * f64::from(i) - 3.0).collect();
        assert!(evenly_spaced(&even));
        assert!(evenly_spaced(&[2.0, 5.0]));
        assert!(evenly_spaced(&[0.0_f32, 0.5, 1.0]));
        // A spacing a sixteenth of the mean off is taken, more is not:
        // wherever it is, and whether it is long or short.
        for at in 1..96 {
            let mut near = even.clone();
            near[at] += 0.006;
            assert!(evenly_spaced(&near), "{at}");
            near[at] += 0.001;
            assert!(!evenly_spaced(&near), "{at}");
            near[at] -= 0.014;
            assert!(!evenly_spaced(&near), "{at}");
        }
        // One spacing short and the rest a little long: the short one is
        // found as a long one is.
        let mut short: Vec<f64> = (0..18).map(f64::from).collect();
        short[9..].iter_mut().for_each(|entry| *entry -= 0.2);
        assert!(!evenly_spaced(&short));
        // No width, entries that repeat, spacing that grows.
        assert!(!evenly_spaced(&[1.0, 1.0]));
        assert!(!evenly_spaced(&[0.0, 1.0, 1.0, 1.0, 2.0]));
        let growing: Vec<f64> = (0..60).map(|i| 1.07_f64.powi(i)).collect();
        assert!(!evenly_spaced(&growing));
    }

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
