//! Shape-preserving cubic interpolation of data in one variable.

use super::hermite::{self, finite_secant};
use crate::{Extrapolate, Hint, InterpolateError, Interpolator, PPoly, Scalar, Value};

/// The shape-preserving cubic interpolant through data points: piecewise cubic
/// with a continuous first derivative, monotone wherever the data is.
///
/// PCHIP: piecewise cubic Hermite interpolating polynomial. The slope at each
/// data point is chosen so that the interpolant does not overshoot: between two neighboring points it stays between their
/// values, and where the data turns (a local maximum or minimum, or a flat run)
/// the slope is zero. The price is smoothness: the second derivative jumps at
/// the data points, where a [`CubicSpline`](crate::CubicSpline) has it
/// continuous.
///
/// The slopes: at an inner point the weighted harmonic mean of the
/// two secants beside it (zero when they differ in sign or one is zero), at an
/// end a three-point estimate limited to keep the shape, and for two points the
/// slope of the line through them. For an array value the rule applies to each
/// column on its own.
///
/// Outside the data the result follows the [`Extrapolate`] mode. It starts as
/// [`Extrapolate::Extend`]. The end pieces continued are cubics
/// like any other: what holds between the data points, monotonicity and no
/// overshoot, does not hold out there.
///
/// The interpolant is a piecewise polynomial:
/// [`as_ppoly`](PchipInterpolator::as_ppoly) lends it, and its derivatives and
/// integrals come from there.
///
/// # Examples
///
/// ```
/// use interpx::{InterpolateError, Interpolator, PchipInterpolator};
///
/// # fn main() -> Result<(), InterpolateError> {
/// // A step: flat, a rise, flat again.
/// let x = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
/// let y = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
/// let pchip = PchipInterpolator::new(&x, &y)?;
/// assert_eq!(pchip.eval(&2.5), Some(0.5));
/// // No overshoot: flat where the data is flat.
/// assert_eq!(pchip.eval(&1.5), Some(0.0));
/// assert_eq!(pchip.eval(&3.5), Some(1.0));
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone)]
pub struct PchipInterpolator<T: Scalar, V: Value<T> = T> {
    poly: PPoly<T, V>,
}

impl<T: Scalar, V: Value<T>> PchipInterpolator<T, V> {
    /// Builds the interpolant through `(x[i], y[i])`.
    ///
    /// # Errors
    ///
    /// The first of these that applies:
    ///
    /// - [`InterpolateError::LengthMismatch`] when `x` and `y` differ in length.
    /// - [`InterpolateError::InsufficientPoints`] when there are fewer than two
    ///   points.
    /// - [`InterpolateError::NonFinite`] when a coordinate, or else a value, is
    ///   NaN or infinite. Its `index` is the position in `x` when a coordinate
    ///   is at fault; only when every coordinate is finite is it a position in
    ///   `y`.
    /// - [`InterpolateError::DuplicatePoints`] or [`InterpolateError::NotSorted`]
    ///   when `x` does not increase strictly.
    /// - [`InterpolateError::NonFinite`] when a secant or a coefficient
    ///   overflows. Its `index` is that of a data point near the overflow.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{InterpolateError, Interpolator, PchipInterpolator};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// // Two columns, each kept in shape on its own: the first rises, the
    /// // second has its maximum at the middle point.
    /// let x = [0.0, 1.0, 2.0];
    /// let y = [[0.0, 0.0], [1.0, 4.0], [5.0, 0.0]];
    /// let pchip = PchipInterpolator::new(&x, &y)?;
    /// let slope_at_the_middle = pchip.as_ppoly().eval_derivative(1.0, 1);
    /// assert!(slope_at_the_middle.is_some_and(|s| s[0] > 0.0 && s[1] == 0.0));
    ///
    /// let err = PchipInterpolator::new(&[0.0], &[1.0]);
    /// assert!(matches!(err, Err(InterpolateError::InsufficientPoints { required: 2, actual: 1 })));
    /// # Ok(())
    /// # }
    /// ```
    pub fn new(x: &[T], y: &[V]) -> Result<Self, InterpolateError> {
        crate::validate::coordinates_and_values(x, y, x.len(), 2)?;
        let last = x.len() - 1;
        // The slopes go where the coefficients will be (`hermite::into_ppoly`).
        let mut buffer = vec![V::zero(); 4 * last];
        let per_width = hermite::per_widths(x);
        slopes(x, &per_width, y, &mut buffer)
            .map_err(|index| InterpolateError::NonFinite { index })?;
        Ok(Self {
            poly: hermite::into_ppoly(x, y, buffer, &per_width)?,
        })
    }

    /// Returns the same interpolant with another behavior outside the data.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{Extrapolate, InterpolateError, Interpolator, PchipInterpolator};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let pchip = PchipInterpolator::new(&[0.0, 1.0, 2.0], &[0.0, 1.0, 4.0])?;
    /// assert_eq!(pchip.extrapolate(), Extrapolate::Extend);
    /// assert!(pchip.eval(&3.0).is_some());
    ///
    /// let pchip = pchip.with_extrapolate(Extrapolate::Clamp);
    /// assert_eq!(pchip.eval(&3.0), Some(4.0));
    /// assert_eq!(pchip.domain(), (0.0, 2.0));
    /// # Ok(())
    /// # }
    /// ```
    #[must_use]
    pub fn with_extrapolate(mut self, extrapolate: Extrapolate) -> Self {
        self.poly = self.poly.with_extrapolate(extrapolate);
        self
    }

    /// Returns the behavior outside the data. See
    /// [`with_extrapolate`](PchipInterpolator::with_extrapolate) for an example.
    pub fn extrapolate(&self) -> Extrapolate {
        self.poly.extrapolate()
    }

    /// Returns the first and the last coordinate of the data. See
    /// [`with_extrapolate`](PchipInterpolator::with_extrapolate) for an example.
    pub fn domain(&self) -> (T, T) {
        let x = self.poly.breakpoints();
        (x[0], x[x.len() - 1])
    }

    /// Returns the interpolant as a piecewise polynomial of degree 3, with the
    /// same extrapolation mode.
    ///
    /// This is the way to its derivatives and integrals: see
    /// [`PPoly::derivative`], [`PPoly::antiderivative`], [`PPoly::integrate`],
    /// and [`PPoly::eval_derivative`]. Its break points are the coordinates of
    /// the data, and it evaluates to the same bits as this interpolant.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{InterpolateError, PchipInterpolator};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// // A peak at 1: the slope there is zero.
    /// let pchip = PchipInterpolator::new(&[0.0, 1.0, 2.0], &[0.0, 3.0, 1.0])?;
    /// let poly = pchip.as_ppoly();
    /// assert_eq!(poly.degree(), 3);
    /// assert_eq!(poly.eval_derivative(1.0, 1), Some(0.0));
    /// # Ok(())
    /// # }
    /// ```
    pub fn as_ppoly(&self) -> &PPoly<T, V> {
        &self.poly
    }

    /// Evaluates at `x`, looking first at the interval in `hint` and its
    /// neighbors.
    ///
    /// The result is bit for bit that of [`eval`](Interpolator::eval). Use this
    /// for queries that come in order; for queries in no order `eval` is faster.
    /// [`PPoly::eval_hinted`] explains why, and
    /// [`eval_batch_into`](Interpolator::eval_batch_into) makes the choice itself.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{Hint, InterpolateError, Interpolator, PchipInterpolator};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let pchip = PchipInterpolator::new(&[0.0, 1.0, 2.0], &[0.0, 1.0, 4.0])?;
    /// let mut hint = Hint::new();
    /// for x in [0.25, 0.5, 1.75] {
    ///     assert_eq!(pchip.eval_hinted(x, &mut hint), pchip.eval(&x));
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn eval_hinted(&self, x: T, hint: &mut Hint) -> Option<V> {
        self.poly.eval_hinted(x, hint)
    }
}

/// Writes the slope at every point of `(x[i], y[i])` to the first `x.len()`
/// entries of `slopes`. `per_width` has the reciprocal of the width of every
/// interval ([`hermite::per_widths`]).
///
/// `x` and `y` have passed the input checks, with two points at least.
///
/// # Errors
///
/// Returns the later point of the first interval whose secant overflows.
pub(crate) fn slopes<T: Scalar, V: Value<T>>(
    x: &[T],
    per_width: &[T],
    y: &[V],
    slopes: &mut [V],
) -> Result<(), usize> {
    let n = x.len();
    let last = n - 1;
    let width = |i: usize| x[i + 1] - x[i];

    // The inner points. The secants of the first two and the last two
    // intervals are kept for the ends.
    let mut before = finite_secant(y, 0, per_width)?;
    let mut first = [before; 2];
    let mut tail = [before; 2];
    for (i, slope) in slopes.iter_mut().enumerate().take(last).skip(1) {
        let after = finite_secant(y, i, per_width)?;
        *slope = inner_slope(width(i - 1), width(i), before, after);
        if i == 1 {
            first[1] = after;
        }
        tail = [before, after];
        before = after;
    }
    if n == 2 {
        // The line through the two points.
        slopes[0] = first[0];
        slopes[1] = first[0];
    } else {
        slopes[0] = end_slope(width(0), width(1), first[0], first[1]);
        slopes[last] = end_slope(width(last - 1), width(last - 2), tail[1], tail[0]);
    }
    Ok(())
}

/// What an evaluation in one interval takes of the coordinates of a line.
///
/// The interpolant in an interval depends on the two points of the interval
/// and on one more on each side, where there is one: the slope at a point
/// looks no farther than its neighbors. So the values at a window of four
/// points at most give what a [`PchipInterpolator`] through the whole line
/// would, to the same bits, without building it. A grid folds the values of
/// one axis this way at each evaluation.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Local<T> {
    /// The first point of the window.
    pub(crate) first: usize,
    /// How many points the window has: 2 to 4.
    pub(crate) len: usize,
    /// Where in the window the lower point of the interval is: 0 or 1.
    at: usize,
    /// The widths of the intervals of the window, and their reciprocals.
    widths: [T; 3],
    per_width: [T; 3],
    /// Whether the upper point of the interval is the last of the line.
    ends: bool,
    /// How far the coordinate is from the lower point of the interval.
    dx: T,
}

impl<T: Scalar> Local<T> {
    /// The window for the coordinate `at` in interval `i` of `x`, which has
    /// two points at least.
    pub(crate) fn new(x: &[T], i: usize, at: T) -> Self {
        let first = i.saturating_sub(1);
        let len = (i + 3).min(x.len()) - first;
        let mut widths = [T::one(); 3];
        let mut per_width = [T::one(); 3];
        for (j, pair) in x[first..first + len].windows(2).enumerate() {
            widths[j] = pair[1] - pair[0];
            per_width[j] = widths[j].recip();
        }
        Self {
            first,
            len,
            at: i - first,
            widths,
            per_width,
            ends: i + 2 == x.len(),
            dx: at - x[i],
        }
    }

    /// The value at the coordinate, from the values at the points of the
    /// window. The entries past its length are not read.
    pub(crate) fn eval<V: Value<T>>(&self, y: &[V; 4]) -> V {
        let secant = |j: usize| V::zero().axpy(self.per_width[j], y[j + 1].axpy(-T::one(), y[j]));
        let (at, h) = (self.at, &self.widths);
        let here = secant(at);
        let slopes = if self.len == 2 {
            // The line through the two points.
            [here, here]
        } else {
            let left = if at == 0 {
                end_slope(h[0], h[1], here, secant(1))
            } else {
                inner_slope(h[at - 1], h[at], secant(at - 1), here)
            };
            let right = if self.ends {
                end_slope(h[at], h[at - 1], here, secant(at - 1))
            } else {
                inner_slope(h[at], h[at + 1], here, secant(at + 1))
            };
            [left, right]
        };
        self.eval_between([y[at], y[at + 1]], slopes)
    }

    /// The value at the coordinate, from the values and the slopes at the two
    /// points of the interval: for a line whose slopes are known already.
    pub(crate) fn eval_between<V: Value<T>>(&self, y: [V; 2], slopes: [V; 2]) -> V {
        let piece = hermite::piece(self.per_width[self.at], y, slopes);
        hermite::at(&piece, self.dx)
    }
}

/// The sign of a number as SciPy's `sign` has it: positive, negative, or zero.
fn sign<T: Scalar>(v: T) -> i8 {
    if v > T::zero() {
        1
    } else if v < T::zero() {
        -1
    } else {
        0
    }
}

/// The slope at an inner point, column by column, from the widths of the two
/// intervals beside it and their secants.
///
/// Zero when the secants differ in sign or one of them is zero: the data turns
/// there, or is flat on one side. Otherwise their harmonic mean, weighted by
/// the widths (Fritsch and Butland), which is never steeper than three times
/// the flatter secant: that keeps each piece monotone.
fn inner_slope<T: Scalar, V: Value<T>>(h_before: T, h_after: T, before: V, after: V) -> V {
    let two = T::one() + T::one();
    let (w1, w2) = (two * h_after + h_before, h_after + two * h_before);
    before.zip_with(after, |m0, m1| {
        // The test for a zero secant is SciPy's, and is kept for that. It does
        // not change a result: one zero secant has another sign than its
        // neighbor, and for two the mean below is `1 / infinity`, zero as well.
        if sign(m0) != sign(m1) || sign(m0) == 0 {
            T::zero()
        } else {
            ((w1 / m0 + w2 / m1) / (w1 + w2)).recip()
        }
    })
}

/// The slope at an end, column by column: `h0` and `m0` are the width and the
/// secant of the interval at the end, `h1` and `m1` those of the next one.
///
/// A three-point estimate, one-sided, limited so that the shape is kept
/// (Moler, Numerical Computing with MATLAB, `pchiptx`): zero when it points
/// against the first secant, and three times that secant at most when the
/// second secant has the other sign.
fn end_slope<T: Scalar, V: Value<T>>(h0: T, h1: T, m0: V, m1: V) -> V {
    let two = T::one() + T::one();
    let three = two + T::one();
    m0.zip_with(m1, |m0, m1| {
        let d = ((two * h0 + h1) * m0 - h0 * m1) / (h0 + h1);
        // The test of the two signs in the second branch is SciPy's, and is kept
        // for that. It does not change a result: with secants of one sign the
        // estimate is under twice the first, so it is never over three times it.
        if sign(d) != sign(m0) {
            T::zero()
        } else if sign(m0) != sign(m1) && d.abs() > three * m0.abs() {
            three * m0
        } else {
            d
        }
    })
}

impl<T: Scalar, V: Value<T>> Interpolator for PchipInterpolator<T, V> {
    type Point = T;
    type Output = V;

    fn eval(&self, x: &T) -> Option<V> {
        self.poly.eval(x)
    }

    fn eval_batch_into(
        &self,
        points: &[T],
        out: &mut [V],
        fill: V,
    ) -> Result<(), InterpolateError> {
        self.poly.eval_batch_into(points, out, fill)
    }
}

/// The stored form of the `serde` feature: the piecewise cubic.
#[cfg(feature = "serde")]
mod stored {
    use super::PchipInterpolator;
    use crate::ppoly::stored::Cubic;
    use crate::{Scalar, Value};
    use serde::{Deserialize, Deserializer, Serialize, Serializer};

    impl<T: Scalar, V: Value<T>> Serialize for PchipInterpolator<T, V> {
        fn serialize<S: Serializer>(&self, to: S) -> Result<S::Ok, S::Error> {
            Cubic::of(&self.poly).serialize(to)
        }
    }

    impl<'de, T: Scalar, V: Value<T>> Deserialize<'de> for PchipInterpolator<T, V> {
        fn deserialize<D: Deserializer<'de>>(from: D) -> Result<Self, D::Error> {
            let poly = Cubic::deserialize(from)?.cubic()?;
            Ok(Self { poly })
        }
    }
}
