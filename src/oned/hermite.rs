//! From values and slopes at the data points to a cubic piecewise polynomial:
//! the cubic Hermite form.
//!
//! The spline interpolators of this module differ in how they choose the slopes.
//! Once the slopes are there, they all end the same way, here.

use crate::{InterpolateError, PPoly, Scalar, Value};

/// The secant of interval `i`: the rise of the data over its width. The
/// reciprocal of the width is left in `per_width[i]`, for the slopes that need
/// it and for [`into_ppoly`].
///
/// # Errors
///
/// Returns [`InterpolateError::NonFinite`] with the later of the two points when
/// the secant overflows.
pub(crate) fn secant<T: Scalar, V: Value<T>>(
    x: &[T],
    y: &[V],
    i: usize,
    per_width: &mut [T],
) -> Result<V, InterpolateError> {
    per_width[i] = (x[i + 1] - x[i]).recip();
    finite_secant(y, i, per_width).map_err(|index| InterpolateError::NonFinite { index })
}

/// The reciprocal of the width of every interval of `x`, as [`secant`] leaves
/// them: for a caller with many lines of values on the same coordinates, which
/// finds them once.
pub(crate) fn per_widths<T: Scalar>(x: &[T]) -> Vec<T> {
    x.windows(2)
        .map(|pair| (pair[1] - pair[0]).recip())
        .collect()
}

/// [`secant`] from the reciprocals of the widths, for a caller that has them
/// ([`per_widths`]) and its own way to name the place: a line of a grid is at
/// an offset in the values of the grid.
///
/// # Errors
///
/// Returns the later of the two points when the secant overflows.
pub(crate) fn finite_secant<T: Scalar, V: Value<T>>(
    y: &[V],
    i: usize,
    per_width: &[T],
) -> Result<V, usize> {
    let rise = y[i + 1].axpy(-T::one(), y[i]);
    let secant = V::zero().axpy(per_width[i], rise);
    if secant.all_finite() {
        Ok(secant)
    } else {
        Err(i + 1)
    }
}

/// The cubic of one interval from the values and the slopes at its two ends,
/// highest power first, in the distance from the lower end. `per_width` is the
/// reciprocal of the width of the interval.
///
/// With the secant `d` of an interval of width `h` and the slopes `s0`, `s1` at
/// its ends, and `t = (s0 + s1 - 2 d) / h`, the piece is
/// `(t / h) dx^3 + ((d - s0) / h - t) dx^2 + s0 dx + y0`: SciPy's
/// `CubicHermiteSpline`.
pub(crate) fn piece<T: Scalar, V: Value<T>>(per_width: T, y: [V; 2], slopes: [V; 2]) -> [V; 4] {
    let two = T::one() + T::one();
    let [left, right] = slopes;
    let secant = V::zero().axpy(per_width, y[1].axpy(-T::one(), y[0]));
    // (s0 + s1 - 2 secant) / h
    let t = V::zero().axpy(per_width, left.axpy(T::one(), right).axpy(-two, secant));
    let cubic = V::zero().axpy(per_width, t);
    let square = V::zero()
        .axpy(per_width, secant.axpy(-T::one(), left))
        .axpy(-T::one(), t);
    [cubic, square, left, y[0]]
}

/// The value of a [`piece`] at the distance `dx` from the lower end of its
/// interval, by Horner's rule: the steps of [`PPoly`], to the same bits.
pub(crate) fn at<T: Scalar, V: Value<T>>(piece: &[V; 4], dx: T) -> V {
    piece.iter().fold(V::zero(), |acc, c| c.axpy(dx, acc))
}

/// Builds the cubic piecewise polynomial through `(x[i], y[i])` with the slopes
/// in the first `x.len()` entries of `buffer`.
///
/// `x` and `y` have passed the input checks, with two points at least.
/// `buffer` has four entries per interval, and `per_width[i]` is the reciprocal
/// of the width of interval `i`, as [`secant`] leaves it.
///
/// The coefficients are written over the slopes, from the last interval to the
/// first: interval `i` needs the slopes at `i` and `i + 1`, and writes at `4 i`
/// and above, where no slope it or an earlier interval needs is kept. So the
/// buffer of the slopes becomes the coefficients, without another allocation.
///
/// Only the first `x.len()` entries of `buffer` are read. The rest is the
/// caller's until this call, and is overwritten here: a caller may keep its
/// own working values there while it finds the slopes (`Akima1D` keeps the
/// secants there). The secants this function needs it computes itself, from
/// `y` and `per_width`; it must not come to read them from `buffer`.
///
/// Each interval gets its [`piece`].
///
/// # Errors
///
/// Returns [`InterpolateError::NonFinite`] with the first interval that has a
/// coefficient that is not finite.
pub(crate) fn into_ppoly<T: Scalar, V: Value<T>>(
    x: &[T],
    y: &[V],
    mut buffer: Vec<V>,
    per_width: &[T],
) -> Result<PPoly<T, V>, InterpolateError> {
    let n = x.len();
    debug_assert_eq!(buffer.len(), 4 * (n - 1));
    // The first interval with a coefficient that is not finite: the loop goes
    // down, so the last one it finds.
    let mut overflowed = None;
    for i in (0..n - 1).rev() {
        let piece = piece(per_width[i], [y[i], y[i + 1]], [buffer[i], buffer[i + 1]]);
        if !piece.iter().all(V::all_finite) {
            overflowed = Some(i);
        }
        buffer[4 * i..4 * i + 4].copy_from_slice(&piece);
    }
    if let Some(index) = overflowed {
        return Err(InterpolateError::NonFinite { index });
    }
    // The caller checked the coordinates, and every coefficient is finite.
    Ok(PPoly::from_checked(x.to_vec(), buffer, 3))
}
