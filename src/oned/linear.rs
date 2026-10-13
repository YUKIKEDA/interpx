//! Linear interpolation of data in one variable.

use crate::{
    search, validate, Extrapolate, Hint, InterpolateError, Interpolator, PPoly, Scalar, Value,
};

/// The polygon through data points: linear interpolation in one variable.
///
/// NumPy's `numpy.interp`. Between two neighboring points `(x[i], y[i])` and
/// `(x[i + 1], y[i + 1])` the value is `y[i] + s (x - x[i])` with the slope
/// `s = (y[i + 1] - y[i]) / (x[i + 1] - x[i])`. At every `x[i]`, the last one
/// included, the result is `y[i]` exactly.
///
/// Outside the data the result follows the [`Extrapolate`] mode. It starts as
/// [`Extrapolate::Clamp`], which gives the first or the last value, as
/// `numpy.interp` does. `numpy.interp` has more arguments, which map as follows:
///
/// - `left` and `right`, the values to return outside: use
///   [`Extrapolate::Forbid`] and [`eval_or`](Interpolator::eval_or).
/// - `period`: not provided. [`Extrapolate::Periodic`] wraps a query into the
///   range of the data, which is another thing (`numpy.interp` wraps the data
///   points as well).
///
/// The data is checked more strictly than NumPy checks it: two points at least
/// (NumPy takes one, as a constant), `x` strictly increasing (NumPy does not
/// look), and finite values (NumPy takes infinities).
///
/// # Examples
///
/// ```
/// use interpx::{Extrapolate, InterpolateError, Interpolator, Linear1D};
///
/// # fn main() -> Result<(), InterpolateError> {
/// let line = Linear1D::new(&[0.0, 1.0, 3.0], &[0.0, 10.0, 30.0])?;
/// assert_eq!(line.eval(&0.5), Some(5.0));
/// assert_eq!(line.eval(&3.0), Some(30.0));
///
/// // Outside, the end values.
/// assert_eq!(line.eval(&-1.0), Some(0.0));
/// assert_eq!(line.eval(&9.0), Some(30.0));
///
/// // `numpy.interp(x, xp, fp, left=-1, right=-1)`.
/// let line = line.with_extrapolate(Extrapolate::Forbid);
/// assert_eq!(line.eval_or(&9.0, -1.0), -1.0);
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone)]
pub struct Linear1D<T: Scalar, V: Value<T> = T> {
    /// Degree 1: for each interval the slope, then the value at its left end.
    poly: PPoly<T, V>,
    /// The value at the last point, which no interval has at its left end.
    last: V,
}

impl<T: Scalar, V: Value<T>> Linear1D<T, V> {
    /// Builds the interpolator from copies of the coordinates and the values.
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
    /// - [`InterpolateError::NonFinite`] when the slope between two points
    ///   overflows. Its `index` is that of the later point.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{InterpolateError, Interpolator, Linear1D};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// // Two outputs per point: the value type is an array.
    /// let path = Linear1D::new(&[0.0, 2.0], &[[0.0, 1.0], [4.0, 1.0]])?;
    /// assert_eq!(path.eval(&0.5), Some([1.0, 1.0]));
    ///
    /// let err = Linear1D::new(&[0.0, 0.0, 1.0], &[1.0, 2.0, 3.0]);
    /// assert!(matches!(err, Err(InterpolateError::DuplicatePoints { first: 0, second: 1 })));
    /// # Ok(())
    /// # }
    /// ```
    pub fn new(x: &[T], y: &[V]) -> Result<Self, InterpolateError> {
        validate::coordinates_and_values(x, y, x.len(), 2)?;
        // Every slope first, then one look at them all: a loop that may leave
        // at every slope is compiled one slope at a time.
        let mut coefficients = vec![V::zero(); 2 * (x.len() - 1)];
        let pairs = x.windows(2).zip(y.windows(2));
        for (piece, (xs, ys)) in coefficients.as_chunks_mut::<2>().0.iter_mut().zip(pairs) {
            let rise = ys[1].axpy(-T::one(), ys[0]);
            piece[0] = V::zero().axpy((xs[1] - xs[0]).recip(), rise);
            piece[1] = ys[0];
        }
        // The values are finite, so a coefficient that is not is a slope.
        if let Some(at) = validate::first_failing(&coefficients, V::all_finite) {
            return Err(InterpolateError::NonFinite { index: at / 2 + 1 });
        }
        // The coordinates passed the checks above, and every slope is finite.
        let poly =
            PPoly::from_checked(x.to_vec(), coefficients, 1).with_extrapolate(Extrapolate::Clamp);
        Ok(Self {
            poly,
            last: y[y.len() - 1],
        })
    }

    /// Returns the same interpolator with another behavior outside the data.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{Extrapolate, InterpolateError, Interpolator, Linear1D};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let line = Linear1D::new(&[0.0, 1.0, 3.0], &[0.0, 10.0, 30.0])?;
    /// // The last segment continued.
    /// assert_eq!(line.clone().with_extrapolate(Extrapolate::Extend).eval(&4.0), Some(40.0));
    /// // The range of the data repeated: 4 is 1.
    /// assert_eq!(line.clone().with_extrapolate(Extrapolate::Periodic).eval(&4.0), Some(10.0));
    /// assert_eq!(line.with_extrapolate(Extrapolate::Forbid).eval(&4.0), None);
    /// # Ok(())
    /// # }
    /// ```
    #[must_use]
    pub fn with_extrapolate(mut self, extrapolate: Extrapolate) -> Self {
        self.poly = self.poly.with_extrapolate(extrapolate);
        self
    }

    /// Returns the behavior outside the data.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{Extrapolate, InterpolateError, Linear1D};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let line = Linear1D::new(&[0.0, 1.0], &[0.0, 2.0])?;
    /// assert_eq!(line.extrapolate(), Extrapolate::Clamp);
    /// # Ok(())
    /// # }
    /// ```
    pub fn extrapolate(&self) -> Extrapolate {
        self.poly.extrapolate()
    }

    /// Returns the first and the last coordinate of the data.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{InterpolateError, Linear1D};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let line = Linear1D::new(&[0.0, 1.0, 3.0], &[0.0, 10.0, 30.0])?;
    /// assert_eq!(line.domain(), (0.0, 3.0));
    /// # Ok(())
    /// # }
    /// ```
    pub fn domain(&self) -> (T, T) {
        let x = self.poly.breakpoints();
        (x[0], x[x.len() - 1])
    }

    /// Returns the interpolant as a piecewise polynomial of degree 1, with the
    /// same extrapolation mode.
    ///
    /// This is the way to its derivatives and integrals: see
    /// [`PPoly::derivative`], [`PPoly::antiderivative`], [`PPoly::integrate`],
    /// and [`PPoly::eval_derivative`]. Its break points are the coordinates of
    /// the data.
    ///
    /// [`PPoly::eval`](Interpolator::eval) computes the last coordinate from the
    /// segment that ends there, so at that one point (and to its right under
    /// [`Extrapolate::Clamp`]) it can differ from this interpolator by a
    /// rounding error.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{InterpolateError, Linear1D};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let line = Linear1D::new(&[0.0, 1.0, 3.0], &[0.0, 10.0, 30.0])?;
    /// let poly = line.as_ppoly();
    /// assert_eq!(poly.breakpoints(), [0.0, 1.0, 3.0]);
    /// // The slope on [1, 3], and the area under the polygon.
    /// assert_eq!(poly.eval_derivative(2.0, 1), Some(10.0));
    /// assert_eq!(poly.integrate(0.0, 3.0), Some(45.0));
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
    /// use interpx::{Hint, InterpolateError, Linear1D};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let line = Linear1D::new(&[0.0, 1.0, 3.0], &[0.0, 10.0, 30.0])?;
    /// let mut hint = Hint::new();
    /// assert_eq!(line.eval_hinted(0.5, &mut hint), Some(5.0));
    /// assert_eq!(line.eval_hinted(2.0, &mut hint), Some(20.0));
    /// assert_eq!(line.eval_hinted(f64::NAN, &mut hint), None);
    /// # Ok(())
    /// # }
    /// ```
    pub fn eval_hinted(&self, x: T, hint: &mut Hint) -> Option<V> {
        search::eval_hinted(
            self.poly.breakpoints(),
            self.poly.extrapolate(),
            x,
            hint,
            |i, x| self.segment(i, x),
        )
    }

    /// The value at `x` on the segment of interval `i`, or on its line continued.
    fn segment(&self, i: usize, x: T) -> V {
        let ends = &self.poly.breakpoints()[i..i + 2];
        // The search gives the right end of an interval for the last coordinate
        // only: every other one is the left end of the next interval, where the
        // formula returns the data value as it is. At the last one the formula
        // can be a rounding error off, so the data value is returned instead
        // (docs/design.ja.md section 6.1 names this comparison).
        if x == ends[1] {
            return self.last;
        }
        let piece = &self.poly.coefficients()[2 * i..2 * i + 2];
        piece[1].axpy(x - ends[0], piece[0])
    }
}

impl<T: Scalar, V: Value<T>> Interpolator for Linear1D<T, V> {
    type Point = T;
    type Output = V;

    fn eval(&self, x: &T) -> Option<V> {
        search::eval(
            self.poly.breakpoints(),
            self.poly.extrapolate(),
            *x,
            |i, x| self.segment(i, x),
        )
    }

    fn eval_batch_into(
        &self,
        points: &[T],
        out: &mut [V],
        fill: V,
    ) -> Result<(), InterpolateError> {
        search::eval_batch_into(
            self.poly.breakpoints(),
            self.poly.extrapolate(),
            points,
            out,
            fill,
            |i, x| self.segment(i, x),
        )
    }
}
