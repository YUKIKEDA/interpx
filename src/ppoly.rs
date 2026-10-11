//! Piecewise polynomials: the representation the 1D splines evaluate through.

use crate::{validate, Extrapolate, InterpolateError, Interpolator, Scalar, Value};

/// The interval the last evaluation fell into, kept by the caller to speed up
/// the next one.
///
/// Queries that come in order mostly hit the same interval or the next. A hint
/// lets [`PPoly::eval_hinted`] look there first and skip the binary search. The
/// result never depends on the hint, only the time does.
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

/// A piecewise polynomial in one variable.
///
/// SciPy's `PPoly`. Between the break points `x[i]` and `x[i + 1]` the value is
///
/// ```text
/// c[i][0] (x - x[i])^k + c[i][1] (x - x[i])^(k - 1) + ... + c[i][k]
/// ```
///
/// with `k` the degree. The coefficients are stored interval by interval, the
/// highest power first, so that one interval's coefficients are adjacent and
/// are read front to back. SciPy stores the same numbers power by power
/// (`c[power][interval]`).
///
/// An interval contains its left break point and not its right one; the last
/// interval contains both. Outside the break points the result follows the
/// [`Extrapolate`] mode, which is [`Extrapolate::Extend`] unless
/// [`with_extrapolate`](PPoly::with_extrapolate) changes it.
///
/// Cloning copies the break points and the coefficients.
///
/// # Examples
///
/// ```
/// use interpx::{Extrapolate, InterpolateError, Interpolator, PPoly};
///
/// # fn main() -> Result<(), InterpolateError> {
/// // Two linear pieces: x on [0, 1], then 1 + 2(x - 1) on [1, 2].
/// let p = PPoly::new(&[0.0, 1.0, 2.0], &[1.0, 0.0, 2.0, 1.0], 1)?;
/// assert_eq!(p.eval(&0.5), Some(0.5));
/// assert_eq!(p.eval(&1.5), Some(2.0));
///
/// // By default the end pieces continue past the ends.
/// assert_eq!(p.eval(&3.0), Some(5.0));
///
/// let p = p.with_extrapolate(Extrapolate::Forbid);
/// assert_eq!(p.eval(&3.0), None);
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone)]
pub struct PPoly<T: Scalar, V: Value<T> = T> {
    breakpoints: Vec<T>,
    coefficients: Vec<V>,
    /// Coefficients per interval: the degree plus one.
    order: usize,
    extrapolate: Extrapolate,
}

impl<T: Scalar, V: Value<T>> PPoly<T, V> {
    /// Builds a piecewise polynomial of the given degree from copies of the break
    /// points and the coefficients.
    ///
    /// `coefficients` holds `degree + 1` numbers for each of the
    /// `breakpoints.len() - 1` intervals, interval by interval, the highest power
    /// first.
    ///
    /// # Errors
    ///
    /// The first of these that applies:
    ///
    /// - [`InterpolateError::LengthMismatch`] when `coefficients` does not hold
    ///   `degree + 1` numbers per interval.
    /// - [`InterpolateError::InsufficientPoints`] when there are fewer than two
    ///   break points.
    /// - [`InterpolateError::NonFinite`] when a break point, or else a coefficient,
    ///   is NaN or infinite.
    /// - [`InterpolateError::DuplicatePoints`] or [`InterpolateError::NotSorted`]
    ///   when the break points do not increase strictly.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{InterpolateError, Interpolator, PPoly};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// // (x - 1)^2 on [1, 3]: one interval, degree 2.
    /// let p = PPoly::new(&[1.0, 3.0], &[1.0, 0.0, 0.0], 2)?;
    /// assert_eq!(p.eval(&2.0), Some(1.0));
    ///
    /// // Two outputs per point: the value type is an array.
    /// let p = PPoly::new(&[0.0, 1.0], &[[1.0, -1.0], [0.0, 5.0]], 1)?;
    /// assert_eq!(p.eval(&0.5), Some([0.5, 4.5]));
    ///
    /// let err = PPoly::new(&[0.0, 1.0], &[1.0, 0.0, 0.0], 1);
    /// assert!(matches!(err, Err(InterpolateError::LengthMismatch { expected: 2, actual: 3 })));
    /// # Ok(())
    /// # }
    /// ```
    pub fn new(
        breakpoints: &[T],
        coefficients: &[V],
        degree: usize,
    ) -> Result<Self, InterpolateError> {
        let order = check(breakpoints, coefficients, degree)?;
        Ok(Self {
            breakpoints: breakpoints.to_vec(),
            coefficients: coefficients.to_vec(),
            order,
            extrapolate: Extrapolate::Extend,
        })
    }

    /// Builds a piecewise polynomial of the given degree, taking the break points
    /// and the coefficients without copying them.
    ///
    /// The layout and the checks are those of [`PPoly::new`].
    ///
    /// # Errors
    ///
    /// The same as [`PPoly::new`].
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{InterpolateError, Interpolator, PPoly};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let p = PPoly::from_owned(vec![0.0, 2.0], vec![3.0, 1.0], 1)?;
    /// assert_eq!(p.eval(&1.0), Some(4.0));
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_owned(
        breakpoints: Vec<T>,
        coefficients: Vec<V>,
        degree: usize,
    ) -> Result<Self, InterpolateError> {
        let order = check(&breakpoints, &coefficients, degree)?;
        Ok(Self {
            breakpoints,
            coefficients,
            order,
            extrapolate: Extrapolate::Extend,
        })
    }

    /// Returns the same polynomial with another behavior outside its domain.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{Extrapolate, InterpolateError, Interpolator, PPoly};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let p = PPoly::new(&[0.0, 1.0], &[2.0, 0.0], 1)?; // 2x
    /// assert_eq!(p.clone().with_extrapolate(Extrapolate::Clamp).eval(&5.0), Some(2.0));
    /// assert_eq!(p.clone().with_extrapolate(Extrapolate::Periodic).eval(&1.25), Some(0.5));
    /// assert_eq!(p.with_extrapolate(Extrapolate::Forbid).eval(&5.0), None);
    /// # Ok(())
    /// # }
    /// ```
    #[must_use]
    pub fn with_extrapolate(mut self, extrapolate: Extrapolate) -> Self {
        self.extrapolate = extrapolate;
        self
    }

    /// Returns the break points, in increasing order.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{InterpolateError, PPoly};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let p = PPoly::new(&[0.0, 1.0, 4.0], &[1.0, 2.0], 0)?;
    /// assert_eq!(p.breakpoints(), [0.0, 1.0, 4.0]);
    /// assert_eq!(p.coefficients(), [1.0, 2.0]);
    /// assert_eq!(p.degree(), 0);
    /// # Ok(())
    /// # }
    /// ```
    pub fn breakpoints(&self) -> &[T] {
        &self.breakpoints
    }

    /// Returns the coefficients: `degree + 1` per interval, interval by interval,
    /// the highest power first. See [`breakpoints`](PPoly::breakpoints) for an
    /// example.
    pub fn coefficients(&self) -> &[V] {
        &self.coefficients
    }

    /// Returns the degree of the pieces. See [`breakpoints`](PPoly::breakpoints)
    /// for an example.
    pub fn degree(&self) -> usize {
        self.order - 1
    }

    /// Returns the behavior outside the domain.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{Extrapolate, InterpolateError, PPoly};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let p = PPoly::new(&[0.0, 1.0], &[1.0], 0)?;
    /// assert_eq!(p.extrapolate(), Extrapolate::Extend);
    /// # Ok(())
    /// # }
    /// ```
    pub fn extrapolate(&self) -> Extrapolate {
        self.extrapolate
    }

    /// Evaluates the polynomial at `x`, looking first at the interval in `hint`
    /// and its neighbors.
    ///
    /// The result is bit for bit that of [`eval`](Interpolator::eval). `hint` is
    /// updated to the interval used, so a run of nearby queries skips the binary
    /// search.
    ///
    /// Use this for queries that come in order. For queries in no order the hint
    /// misses, and a miss costs more than the binary search alone: call `eval`
    /// instead, or [`eval_batch_into`](Interpolator::eval_batch_into), which
    /// stops using its hint while it keeps missing.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{Hint, InterpolateError, PPoly};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let p = PPoly::new(&[0.0, 1.0, 2.0], &[1.0, 0.0, 2.0, 1.0], 1)?;
    /// let mut hint = Hint::new();
    /// assert_eq!(p.eval_hinted(1.5, &mut hint), Some(2.0));
    /// assert_eq!(p.eval_hinted(f64::NAN, &mut hint), None);
    /// # Ok(())
    /// # }
    /// ```
    pub fn eval_hinted(&self, x: T, hint: &mut Hint) -> Option<V> {
        let start = hint.interval;
        let (interval, x) = self.resolve(x, |p, x| p.locate_near(x, start))?;
        hint.interval = interval;
        Some(self.horner(interval, x))
    }

    /// The index of the last interval.
    fn last(&self) -> usize {
        // Two break points at least, by construction.
        self.breakpoints.len() - 2
    }

    /// The interval of `x`, for `x` between the first and the last break point.
    fn locate(&self, x: T) -> usize {
        let after = self.breakpoints.partition_point(|b| *b <= x);
        after.saturating_sub(1).min(self.last())
    }

    /// Whether interval `i` contains `x`: closed on the left, open on the right,
    /// and closed on the right too for the last interval.
    fn contains(&self, i: usize, x: T) -> bool {
        let (left, right) = (self.breakpoints[i], self.breakpoints[i + 1]);
        left <= x && (x < right || (i == self.last() && x <= right))
    }

    /// The same interval as [`locate`](Self::locate), trying `start` and its two
    /// neighbors before the binary search.
    fn locate_near(&self, x: T, start: usize) -> usize {
        let last = self.last();
        let i = start.min(last);
        if self.contains(i, x) {
            i
        } else if i < last && self.contains(i + 1, x) {
            i + 1
        } else if i > 0 && self.contains(i - 1, x) {
            i - 1
        } else {
            self.locate(x)
        }
    }

    /// Whether a hint pays for this run of points: they go one way, and on
    /// average no more than one interval apart.
    ///
    /// A hinted search starts from the interval of the point before, so each
    /// point waits for the one before it. For points that are close together the
    /// wait is short and the search is skipped. For points in no order the wait
    /// is a whole binary search, and plain searches, which do not wait for each
    /// other, are about twice as fast (measured in `benches/ppoly.rs`).
    fn hint_pays(&self, points: &[T]) -> bool {
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
                    let crossed = self.breakpoints.partition_point(|b| *b <= high)
                        - self.breakpoints.partition_point(|b| *b <= low);
                    crossed <= points.len()
                })
    }

    /// Maps a query to the interval to evaluate and the point to evaluate it at,
    /// or `None` when the query has no value.
    fn resolve(&self, x: T, locate: impl FnOnce(&Self, T) -> usize) -> Option<(usize, T)> {
        if x.is_nan() {
            return None;
        }
        let last = self.last();
        let (lo, hi) = (self.breakpoints[0], self.breakpoints[last + 1]);
        if lo <= x && x <= hi {
            return Some((locate(self, x), x));
        }
        match self.extrapolate {
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
                Some((locate(self, wrapped), wrapped))
            }
        }
    }

    /// Evaluates the polynomial of interval `i` at `x` by Horner's rule.
    fn horner(&self, i: usize, x: T) -> V {
        let dx = x - self.breakpoints[i];
        let piece = &self.coefficients[i * self.order..(i + 1) * self.order];
        piece.iter().fold(V::zero(), |acc, c| c.axpy(dx, acc))
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

/// Runs the input checks of [`PPoly::new`] and returns the number of
/// coefficients per interval.
fn check<T: Scalar, V: Value<T>>(
    breakpoints: &[T],
    coefficients: &[V],
    degree: usize,
) -> Result<usize, InterpolateError> {
    let intervals = breakpoints.len().saturating_sub(1);
    // A degree too large to count cannot match any slice.
    let order = degree.saturating_add(1);
    let expected = intervals.saturating_mul(order);
    validate::coordinates_and_values(breakpoints, coefficients, expected, 2)?;
    Ok(order)
}

impl<T: Scalar, V: Value<T>> Interpolator for PPoly<T, V> {
    type Point = T;
    type Output = V;

    fn eval(&self, x: &T) -> Option<V> {
        let (interval, x) = self.resolve(*x, Self::locate)?;
        Some(self.horner(interval, x))
    }

    fn eval_batch_into(
        &self,
        points: &[T],
        out: &mut [V],
        fill: V,
    ) -> Result<(), InterpolateError> {
        if points.len() != out.len() {
            return Err(InterpolateError::LengthMismatch {
                expected: points.len(),
                actual: out.len(),
            });
        }
        // By runs of `BATCH_RUN` points: with a hint where it pays, plainly where
        // it does not. Either way each point gets the interval `eval` finds.
        let mut hint = Hint::new();
        for (xs, slots) in points.chunks(BATCH_RUN).zip(out.chunks_mut(BATCH_RUN)) {
            if xs.len() >= MIN_HINT_RUN && self.hint_pays(xs) {
                for (x, slot) in xs.iter().zip(slots) {
                    *slot = self.eval_hinted(*x, &mut hint).unwrap_or(fill);
                }
            } else {
                for (x, slot) in xs.iter().zip(slots) {
                    *slot = self.eval(x).unwrap_or(fill);
                }
            }
        }
        Ok(())
    }
}
