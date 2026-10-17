//! The trait every interpolator implements.

use crate::InterpolateError;

/// An interpolator: a function built from data that can be evaluated at a point.
///
/// Only [`eval`](Interpolator::eval) must be written. The other methods are
/// provided, and an implementation may replace them with a faster version that
/// gives the same results.
///
/// `eval`, `eval_or`, and `eval_batch_into` of the interpolators in this crate do
/// not allocate, with one exception: an `RBFInterpolator` over neighbors
/// allocates a working space at its first evaluation on a thread, and widens it
/// when an interpolator with more neighbors is evaluated on that thread.
///
/// # Examples
///
/// ```
/// use interpx::Interpolator;
///
/// /// The step function that is 0 below 1/2 and 1 from there on, on [0, 1].
/// struct Step;
///
/// impl Interpolator for Step {
///     type Point = f64;
///     type Output = f64;
///
///     fn eval(&self, x: &f64) -> Option<f64> {
///         (0.0..=1.0)
///             .contains(x)
///             .then(|| if *x < 0.5 { 0.0 } else { 1.0 })
///     }
/// }
///
/// assert_eq!(Step.eval(&0.75), Some(1.0));
/// assert_eq!(Step.eval(&2.0), None);
/// assert_eq!(Step.eval_or(&2.0, -1.0), -1.0);
/// assert_eq!(Step.eval_batch(&[0.25, 3.0]), [Some(0.0), None]);
/// ```
pub trait Interpolator {
    /// The type of a query point: a scalar in one dimension, `[T; D]` in `D`.
    type Point;

    /// The type of a result: a scalar, or `[T; M]` for `M` outputs.
    type Output;

    /// Evaluates the interpolator at `point`.
    ///
    /// Returns `None` when `point` is outside the domain and the interpolator does
    /// not extrapolate, and when a coordinate of `point` is NaN. The scattered
    /// types also return `None` for a coordinate that is infinite, and
    /// `RBFInterpolator` for a result that is not finite.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{Extrapolate, InterpolateError, Interpolator, Linear1D};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let line = Linear1D::new(&[0.0, 1.0], &[0.0, 10.0])?.with_extrapolate(Extrapolate::Forbid);
    /// assert_eq!(line.eval(&0.5), Some(5.0));
    /// assert_eq!(line.eval(&2.0), None);
    /// assert_eq!(line.eval(&f64::NAN), None);
    /// # Ok(())
    /// # }
    /// ```
    fn eval(&self, point: &Self::Point) -> Option<Self::Output>;

    /// Evaluates the interpolator at `point`, or returns `fill` where
    /// [`eval`](Interpolator::eval) returns `None`.
    ///
    /// This is SciPy's `fill_value`, given at the call. Pass NaN to mark the
    /// points that have no value.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{Extrapolate, InterpolateError, Interpolator, Linear1D};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let line = Linear1D::new(&[0.0, 1.0], &[0.0, 10.0])?.with_extrapolate(Extrapolate::Forbid);
    /// assert_eq!(line.eval_or(&0.5, -1.0), 5.0);
    /// assert_eq!(line.eval_or(&2.0, -1.0), -1.0);
    /// # Ok(())
    /// # }
    /// ```
    fn eval_or(&self, point: &Self::Point, fill: Self::Output) -> Self::Output {
        self.eval(point).unwrap_or(fill)
    }

    /// Evaluates the interpolator at every point and returns the results in a new
    /// vector.
    ///
    /// This allocates the vector. Use
    /// [`eval_batch_into`](Interpolator::eval_batch_into) to write into a buffer
    /// that already exists.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{Extrapolate, InterpolateError, Interpolator, Linear1D};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let line = Linear1D::new(&[0.0, 1.0], &[0.0, 10.0])?.with_extrapolate(Extrapolate::Forbid);
    /// assert_eq!(line.eval_batch(&[0.5, 2.0]), [Some(5.0), None]);
    /// # Ok(())
    /// # }
    /// ```
    fn eval_batch(&self, points: &[Self::Point]) -> Vec<Option<Self::Output>> {
        points.iter().map(|p| self.eval(p)).collect()
    }

    /// Evaluates the interpolator at every point and writes the results to `out`,
    /// with `fill` where [`eval`](Interpolator::eval) returns `None`.
    ///
    /// # Errors
    ///
    /// Returns [`InterpolateError::LengthMismatch`] when `points` and `out` differ
    /// in length. `out` is then left as it was.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{InterpolateError, Interpolator};
    ///
    /// struct Double;
    ///
    /// impl Interpolator for Double {
    ///     type Point = f64;
    ///     type Output = f64;
    ///
    ///     fn eval(&self, x: &f64) -> Option<f64> {
    ///         (0.0..=1.0).contains(x).then(|| 2.0 * x)
    ///     }
    /// }
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let mut out = [0.0; 2];
    /// Double.eval_batch_into(&[0.25, 9.0], &mut out, f64::INFINITY)?;
    /// assert_eq!(out, [0.5, f64::INFINITY]);
    ///
    /// let err = Double.eval_batch_into(&[0.25], &mut out, 0.0);
    /// assert_eq!(
    ///     err,
    ///     Err(InterpolateError::LengthMismatch { expected: 1, actual: 2 }),
    /// );
    /// # Ok(())
    /// # }
    /// ```
    fn eval_batch_into(
        &self,
        points: &[Self::Point],
        out: &mut [Self::Output],
        fill: Self::Output,
    ) -> Result<(), InterpolateError>
    where
        Self::Output: Clone,
    {
        if points.len() != out.len() {
            return Err(InterpolateError::LengthMismatch {
                expected: points.len(),
                actual: out.len(),
            });
        }
        for (point, slot) in points.iter().zip(out.iter_mut()) {
            *slot = self.eval(point).unwrap_or_else(|| fill.clone());
        }
        Ok(())
    }
}
