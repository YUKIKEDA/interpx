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
    let rise = y[i + 1].axpy(-T::one(), y[i]);
    let secant = V::zero().axpy(per_width[i], rise);
    if secant.all_finite() {
        Ok(secant)
    } else {
        Err(InterpolateError::NonFinite { index: i + 1 })
    }
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
/// With the secant `d` of an interval of width `h` and the slopes `s0`, `s1` at
/// its ends, and `t = (s0 + s1 - 2 d) / h`, the piece is
/// `(t / h) dx^3 + ((d - s0) / h - t) dx^2 + s0 dx + y0`: SciPy's
/// `CubicHermiteSpline`.
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
    let two = T::one() + T::one();
    // The first interval with a coefficient that is not finite: the loop goes
    // down, so the last one it finds.
    let mut overflowed = None;
    for i in (0..n - 1).rev() {
        let per_width = per_width[i];
        let secant = V::zero().axpy(per_width, y[i + 1].axpy(-T::one(), y[i]));
        let (left, right) = (buffer[i], buffer[i + 1]);
        // (s[i] + s[i + 1] - 2 secant) / h
        let t = V::zero().axpy(per_width, left.axpy(T::one(), right).axpy(-two, secant));
        let cubic = V::zero().axpy(per_width, t);
        let square = V::zero()
            .axpy(per_width, secant.axpy(-T::one(), left))
            .axpy(-T::one(), t);
        let piece = [cubic, square, left, y[i]];
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
