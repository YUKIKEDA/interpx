//! Cubic interpolation of data in one variable with given slopes.

use super::hermite;
use crate::{Extrapolate, Hint, InterpolateError, Interpolator, PPoly, Scalar, Value};

/// The cubic interpolant through data points with a given slope at each:
/// piecewise cubic with a continuous first derivative.
///
/// SciPy's `CubicHermiteSpline`. Where [`CubicSpline`](crate::CubicSpline),
/// [`PchipInterpolator`](crate::PchipInterpolator), and
/// [`Akima1D`](crate::Akima1D) each choose the slopes by a rule of their own,
/// this type takes them as they are given. Use it when the derivative of the
/// data is known (the solution of a differential equation comes with one), or
/// to interpolate with slopes chosen by another rule.
///
/// The slopes have the type of the values: for an array value there is one
/// slope per column, and each column is interpolated on its own.
///
/// On each interval the constant and the linear coefficient are the value and
/// the slope given at its left end, bit for bit.
///
/// Outside the data the result follows the [`Extrapolate`] mode. It starts as
/// [`Extrapolate::Extend`], as in SciPy. [`Extrapolate::Periodic`] repeats the
/// interpolant as it is: nothing checks that the values and the slopes at the
/// two ends agree, and where they do not, the result jumps from one period to
/// the next.
///
/// The interpolant is a piecewise polynomial:
/// [`as_ppoly`](CubicHermiteSpline::as_ppoly) lends it, and its derivatives and
/// integrals come from there.
///
/// # Examples
///
/// ```
/// use interpx::{CubicHermiteSpline, InterpolateError, Interpolator};
///
/// # fn main() -> Result<(), InterpolateError> {
/// // From 0 to 1 with no slope at either end: 3 x^2 - 2 x^3.
/// let step = CubicHermiteSpline::new(&[0.0, 1.0], &[0.0, 1.0], &[0.0, 0.0])?;
/// assert_eq!(step.eval(&0.5), Some(0.5));
/// assert_eq!(step.eval(&0.25), Some(0.15625));
/// assert_eq!(step.as_ppoly().eval_derivative(1.0, 1), Some(0.0));
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone)]
pub struct CubicHermiteSpline<T: Scalar, V: Value<T> = T> {
    poly: PPoly<T, V>,
}

impl<T: Scalar, V: Value<T>> CubicHermiteSpline<T, V> {
    /// Builds the interpolant through `(x[i], y[i])` with the slope `dydx[i]`
    /// there.
    ///
    /// # Errors
    ///
    /// The first of these that applies:
    ///
    /// - [`InterpolateError::LengthMismatch`] when `y`, or else `dydx`, differs
    ///   from `x` in length.
    /// - [`InterpolateError::InsufficientPoints`] when there are fewer than two
    ///   points.
    /// - [`InterpolateError::NonFinite`] when a coordinate, or else a value, or
    ///   else a slope, is NaN or infinite. Its `index` is the position in `x`
    ///   when a coordinate is at fault; when every coordinate is finite it is a
    ///   position in `y`; only when every value is finite as well is it a
    ///   position in `dydx`.
    /// - [`InterpolateError::DuplicatePoints`] or [`InterpolateError::NotSorted`]
    ///   when `x` does not increase strictly.
    /// - [`InterpolateError::NonFinite`] when a coefficient overflows. Its
    ///   `index` is that of the first interval with such a coefficient.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{CubicHermiteSpline, InterpolateError, Interpolator};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// // x^3 from its values and its derivative 3 x^2 at three points: a cubic
    /// // is reproduced.
    /// let x = [-1.0, 0.0, 2.0];
    /// let cube = CubicHermiteSpline::new(&x, &[-1.0, 0.0, 8.0], &[3.0, 0.0, 12.0])?;
    /// assert_eq!(cube.eval(&1.0), Some(1.0));
    /// assert_eq!(cube.eval(&-0.5), Some(-0.125));
    ///
    /// // Two columns, each with its own slopes.
    /// let y = [[0.0, 1.0], [1.0, 1.0]];
    /// let dydx = [[1.0, 0.0], [1.0, 0.0]];
    /// let pair = CubicHermiteSpline::new(&[0.0, 1.0], &y, &dydx)?;
    /// assert_eq!(pair.eval(&0.5), Some([0.5, 1.0]));
    ///
    /// let err = CubicHermiteSpline::new(&[0.0, 1.0], &[0.0, 1.0], &[0.0]);
    /// assert!(matches!(err, Err(InterpolateError::LengthMismatch { expected: 2, actual: 1 })));
    /// # Ok(())
    /// # }
    /// ```
    pub fn new(x: &[T], y: &[V], dydx: &[V]) -> Result<Self, InterpolateError> {
        crate::validate::coordinates_values_and_slopes(x, y, dydx, 2)?;
        let intervals = x.len() - 1;
        // The slopes go where the coefficients will be (`hermite::into_ppoly`).
        let mut buffer = Vec::with_capacity(4 * intervals);
        buffer.extend_from_slice(dydx);
        buffer.resize(4 * intervals, V::zero());
        let per_width: Vec<T> = x
            .windows(2)
            .map(|pair| (pair[1] - pair[0]).recip())
            .collect();
        Ok(Self {
            poly: hermite::into_ppoly(x, y, buffer, &per_width)?,
        })
    }

    /// Returns the same interpolant with another behavior outside the data.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{CubicHermiteSpline, Extrapolate, InterpolateError, Interpolator};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let x = [0.0, 1.0, 2.0];
    /// let spline = CubicHermiteSpline::new(&x, &[0.0, 1.0, 4.0], &[0.0, 2.0, 4.0])?;
    /// assert_eq!(spline.extrapolate(), Extrapolate::Extend);
    /// assert!(spline.eval(&3.0).is_some());
    ///
    /// let spline = spline.with_extrapolate(Extrapolate::Clamp);
    /// assert_eq!(spline.eval(&3.0), Some(4.0));
    /// assert_eq!(spline.domain(), (0.0, 2.0));
    /// # Ok(())
    /// # }
    /// ```
    #[must_use]
    pub fn with_extrapolate(mut self, extrapolate: Extrapolate) -> Self {
        self.poly = self.poly.with_extrapolate(extrapolate);
        self
    }

    /// Returns the behavior outside the data. See
    /// [`with_extrapolate`](CubicHermiteSpline::with_extrapolate) for an example.
    pub fn extrapolate(&self) -> Extrapolate {
        self.poly.extrapolate()
    }

    /// Returns the first and the last coordinate of the data. See
    /// [`with_extrapolate`](CubicHermiteSpline::with_extrapolate) for an example.
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
    /// use interpx::{CubicHermiteSpline, InterpolateError};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let x = [0.0, 1.0, 3.0];
    /// let spline = CubicHermiteSpline::new(&x, &[1.0, 2.0, 0.0], &[0.5, -1.0, 4.0])?;
    /// let poly = spline.as_ppoly();
    /// assert_eq!(poly.degree(), 3);
    /// // The slopes given, at the left end of each interval.
    /// assert_eq!(poly.eval_derivative(0.0, 1), Some(0.5));
    /// assert_eq!(poly.eval_derivative(1.0, 1), Some(-1.0));
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
    /// use interpx::{CubicHermiteSpline, Hint, InterpolateError, Interpolator};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let x = [0.0, 1.0, 2.0];
    /// let spline = CubicHermiteSpline::new(&x, &[0.0, 1.0, 4.0], &[0.0, 2.0, 4.0])?;
    /// let mut hint = Hint::new();
    /// for x in [0.25, 0.5, 1.75] {
    ///     assert_eq!(spline.eval_hinted(x, &mut hint), spline.eval(&x));
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn eval_hinted(&self, x: T, hint: &mut Hint) -> Option<V> {
        self.poly.eval_hinted(x, hint)
    }
}

impl<T: Scalar, V: Value<T>> Interpolator for CubicHermiteSpline<T, V> {
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
