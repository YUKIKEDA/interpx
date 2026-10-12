//! The cubic spline through data in one variable.

use super::hermite::{self, secant};
use crate::tridiagonal::{self, Breakdown};
use crate::{
    validate, BoundaryCondition, EndCondition, Extrapolate, Hint, InterpolateError, Interpolator,
    PPoly, ParamReason, Scalar, Value,
};

/// The cubic spline through data points: piecewise cubic, with a continuous
/// first and second derivative.
///
/// SciPy's `CubicSpline`. The [`BoundaryCondition`] fixes the two degrees of
/// freedom the data leaves: not-a-knot (SciPy's default), periodic, or a first
/// or second derivative at each end.
///
/// With few points the spline is what SciPy makes of them:
///
/// - Two points: a not-a-knot end takes the slope of the line through them, so
///   not-a-knot at both ends gives that line.
/// - Three points with not-a-knot at both ends: the parabola through them.
///
/// Outside the data the result follows the [`Extrapolate`] mode. It starts as
/// [`Extrapolate::Extend`], and as [`Extrapolate::Periodic`] for the periodic
/// boundary, as in SciPy.
///
/// The interpolant is a piecewise polynomial: [`as_ppoly`](CubicSpline::as_ppoly)
/// lends it, and its derivatives and integrals come from there.
///
/// The data is checked more strictly than SciPy checks it: `x` must increase
/// strictly as given (SciPy raises for a repeated `x` too, but this crate does
/// not sort), and the values must be finite.
///
/// # Examples
///
/// ```
/// use interpx::{BoundaryCondition, CubicSpline, InterpolateError, Interpolator};
///
/// # fn main() -> Result<(), InterpolateError> {
/// let x = [0.0, 1.0, 2.0, 3.0, 4.0];
/// let y = [0.0, 1.0, 0.0, 1.0, 0.0];
/// let spline = CubicSpline::new(&x, &y, BoundaryCondition::default())?;
/// assert_eq!(spline.eval(&1.0), Some(1.0));
/// assert!(spline.eval(&0.5).is_some_and(|v| v > 0.5));
///
/// // The slope at a data point, from the piecewise polynomial.
/// let natural = CubicSpline::new(&x, &y, BoundaryCondition::natural())?;
/// let curvature = natural.as_ppoly().eval_derivative(0.0, 2);
/// assert!(curvature.is_some_and(|c: f64| c.abs() < 1e-12));
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone)]
pub struct CubicSpline<T: Scalar, V: Value<T> = T> {
    poly: PPoly<T, V>,
}

impl<T: Scalar, V: Value<T>> CubicSpline<T, V> {
    /// Builds the spline through `(x[i], y[i])` with the given boundary
    /// condition.
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
    /// - [`InterpolateError::InvalidParameter`] with the name `"bc"` and
    ///   [`ParamReason::NonFinite`] when a derivative of the boundary condition
    ///   is NaN or infinite.
    /// - [`InterpolateError::NotPeriodic`] when the boundary is periodic and the
    ///   first and the last value differ by more than rounding (SciPy's
    ///   tolerance: `1e-15`, relative and absolute, for `f64`).
    /// - [`InterpolateError::NonFinite`] when a slope or a coefficient
    ///   overflows. Its `index` is that of a data point near the overflow.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{
    ///     BoundaryCondition, CubicSpline, EndCondition, InterpolateError, Interpolator,
    /// };
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// // x^3 on [0, 1] is a cubic spline: its slopes at the ends are 0 and 3.
    /// let ends = BoundaryCondition::Ends(
    ///     EndCondition::FirstDerivative(0.0),
    ///     EndCondition::FirstDerivative(3.0),
    /// );
    /// let cube = CubicSpline::new(&[0.0, 1.0], &[0.0, 1.0], ends)?;
    /// assert!(cube.eval(&0.5).is_some_and(|v: f64| (v - 0.125).abs() < 1e-15));
    ///
    /// // A closed curve: two columns, periodic.
    /// let t = [0.0, 1.0, 2.0, 3.0, 4.0];
    /// let points = [[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -1.0], [1.0, 0.0]];
    /// let curve = CubicSpline::new(&t, &points, BoundaryCondition::Periodic)?;
    /// assert_eq!(curve.eval(&5.0), curve.eval(&1.0));
    ///
    /// let open = [[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -1.0], [2.0, 0.0]];
    /// let err = CubicSpline::new(&t, &open, BoundaryCondition::Periodic);
    /// assert!(matches!(err, Err(InterpolateError::NotPeriodic)));
    /// # Ok(())
    /// # }
    /// ```
    pub fn new(x: &[T], y: &[V], bc: BoundaryCondition<V>) -> Result<Self, InterpolateError> {
        validate::coordinates_and_values(x, y, x.len(), 2)?;
        let n = x.len();
        if let BoundaryCondition::Ends(start, end) = bc {
            for end in [start, end] {
                let (EndCondition::FirstDerivative(v) | EndCondition::SecondDerivative(v)) = end
                else {
                    continue;
                };
                if !v.all_finite() {
                    return Err(InterpolateError::InvalidParameter {
                        name: "bc",
                        reason: ParamReason::NonFinite,
                    });
                }
            }
        }
        let periodic = matches!(bc, BoundaryCondition::Periodic);
        if periodic && !y[0].within(&y[n - 1], periodic_tolerance()) {
            return Err(InterpolateError::NotPeriodic);
        }

        // The slopes are solved into the buffer the coefficients will be in
        // (`hermite::into_ppoly`). The secants are computed where they are used,
        // so the system takes one more allocation: its three diagonals, the
        // reciprocal of every width (a division made once, for the secants and
        // again for the coefficients), and the work space of the cyclic solver.
        let mut buffer = vec![V::zero(); 4 * (n - 1)];
        let mut scratch = vec![T::zero(); if periodic { 5 * n } else { 4 * n }];
        match bc {
            BoundaryCondition::Periodic => periodic_slopes(x, y, &mut buffer[..n], &mut scratch)?,
            BoundaryCondition::NotAKnot => {
                let end = EndCondition::NotAKnot;
                end_slopes(x, y, end, end, &mut buffer[..n], &mut scratch)?;
            }
            BoundaryCondition::Ends(start, end) => {
                end_slopes(x, y, start, end, &mut buffer[..n], &mut scratch)?;
            }
        }

        let per_width = &scratch[PER_WIDTH * n..PER_WIDTH * n + n - 1];
        let poly = hermite::into_ppoly(x, y, buffer, per_width)?;
        Ok(Self {
            poly: if periodic {
                poly.with_extrapolate(Extrapolate::Periodic)
            } else {
                poly
            },
        })
    }

    /// Returns the same spline with another behavior outside the data.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{BoundaryCondition, CubicSpline, Extrapolate, InterpolateError, Interpolator};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let spline = CubicSpline::new(&[0.0, 1.0, 2.0], &[0.0, 1.0, 4.0], BoundaryCondition::default())?;
    /// assert_eq!(spline.extrapolate(), Extrapolate::Extend);
    /// // The parabola x^2 continued.
    /// assert!(spline.eval(&3.0).is_some_and(|v: f64| (v - 9.0).abs() < 1e-12));
    ///
    /// let spline = spline.with_extrapolate(Extrapolate::Forbid);
    /// assert_eq!(spline.eval(&3.0), None);
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
    /// [`with_extrapolate`](CubicSpline::with_extrapolate) for an example.
    pub fn extrapolate(&self) -> Extrapolate {
        self.poly.extrapolate()
    }

    /// Returns the first and the last coordinate of the data. See
    /// [`with_extrapolate`](CubicSpline::with_extrapolate) for an example.
    pub fn domain(&self) -> (T, T) {
        let x = self.poly.breakpoints();
        (x[0], x[x.len() - 1])
    }

    /// Returns the spline as a piecewise polynomial of degree 3, with the same
    /// extrapolation mode.
    ///
    /// This is the way to its derivatives and integrals: see
    /// [`PPoly::derivative`], [`PPoly::antiderivative`], [`PPoly::integrate`],
    /// and [`PPoly::eval_derivative`]. Its break points are the coordinates of
    /// the data, and it evaluates to the same bits as this spline.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{BoundaryCondition, CubicSpline, InterpolateError};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// // The parabola x^2 through three points.
    /// let spline = CubicSpline::new(&[0.0, 1.0, 2.0], &[0.0, 1.0, 4.0], BoundaryCondition::default())?;
    /// let poly = spline.as_ppoly();
    /// assert_eq!(poly.degree(), 3);
    /// assert!(poly.eval_derivative(1.5, 1).is_some_and(|d: f64| (d - 3.0).abs() < 1e-12));
    /// assert!(poly.integrate(0.0, 2.0).is_some_and(|a: f64| (a - 8.0 / 3.0).abs() < 1e-12));
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
    /// use interpx::{BoundaryCondition, CubicSpline, Hint, InterpolateError, Interpolator};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let spline = CubicSpline::new(&[0.0, 1.0, 2.0], &[0.0, 1.0, 4.0], BoundaryCondition::default())?;
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

/// How far apart the two end values of a periodic spline may be: within
/// `tolerance * (1 + |last|)`, column by column. SciPy uses `1e-15` for `f64`,
/// which is about 4.5 times the rounding unit; the same multiple serves `f32`.
pub(super) fn periodic_tolerance<T: Scalar>() -> T {
    let two = T::one() + T::one();
    (two * two + two.recip()) * T::epsilon()
}

/// `a u + b v`.
fn combine<T: Scalar, V: Value<T>>(a: T, u: V, b: T, v: V) -> V {
    V::zero().axpy(a, u).axpy(b, v)
}

/// The error for a solve that broke down: with data that passed the checks, a
/// pivot is zero or not finite only when an entry overflowed.
fn overflow(breakdown: Breakdown, points: usize) -> InterpolateError {
    InterpolateError::NonFinite {
        index: breakdown.row.min(points - 1),
    }
}

/// Which part of the scratch space holds the reciprocal widths: after the three
/// diagonals, each as long as the data.
const PER_WIDTH: usize = 3;

/// Writes the slopes at the data points of a periodic spline to `slopes`. The
/// last one repeats the first.
///
/// The unknowns are the slopes at all points but the last, one per interval;
/// row `i` is the inner row of a spline, with its indices taken around the
/// period. `scratch` has room for the three diagonals, the reciprocal widths,
/// and the work space of the cyclic solver: five times the number of points.
fn periodic_slopes<T: Scalar, V: Value<T>>(
    x: &[T],
    y: &[V],
    slopes: &mut [V],
    scratch: &mut [T],
) -> Result<(), InterpolateError> {
    let n = x.len();
    let m = n - 1;
    let three = T::one() + T::one() + T::one();
    let width = |i: usize| x[i + 1] - x[i];
    let (lower, rest) = scratch.split_at_mut(n);
    let (diag, rest) = rest.split_at_mut(n);
    let (upper, rest) = rest.split_at_mut(n);
    let (per_width, work) = rest.split_at_mut(n);
    // Two points: the end values are equal, so the secant is zero (or a rounding
    // error), and SciPy takes it for the slope at both ends.
    if m == 1 {
        slopes.fill(secant(x, y, 0, per_width)?);
        return Ok(());
    }
    // Three points leave two unknowns, fewer than a cyclic system has rows.
    // SciPy's formula for them gives every point the same slope.
    if m == 2 {
        let (first, second) = (secant(x, y, 0, per_width)?, secant(x, y, 1, per_width)?);
        let (a, b) = (per_width[0], per_width[1]);
        let sum = combine(a, first, b, second);
        slopes.fill(V::zero().axpy((a + b).recip(), sum));
        return Ok(());
    }
    // Row `i` joins interval `i` with the one before it; for row 0 that is the
    // last interval.
    let (mut width_before, mut secant_before) = (width(m - 1), secant(x, y, m - 1, per_width)?);
    for i in 0..m {
        let (h, d) = (width(i), secant(x, y, i, per_width)?);
        lower[i] = h;
        diag[i] = (width_before + h) * (T::one() + T::one());
        upper[i] = width_before;
        slopes[i] = V::zero().axpy(three, combine(h, secant_before, width_before, d));
        (width_before, secant_before) = (h, d);
    }
    tridiagonal::solve_cyclic(
        &lower[..m],
        &mut diag[..m],
        &upper[..m],
        &mut slopes[..m],
        &mut work[..m],
    )
    .map_err(|b| overflow(b, n))?;
    slopes[m] = slopes[0];
    Ok(())
}

/// Writes the slopes at the data points of a spline with a condition at each
/// end to `slopes`. `scratch` has room for the three diagonals and the
/// reciprocal widths: four times the number of points.
fn end_slopes<T: Scalar, V: Value<T>>(
    x: &[T],
    y: &[V],
    start: EndCondition<V>,
    end: EndCondition<V>,
    slopes: &mut [V],
    scratch: &mut [T],
) -> Result<(), InterpolateError> {
    let n = x.len();
    let last = n - 1;
    let two = T::one() + T::one();
    let three = two + T::one();
    let half = two.recip();
    let width = |i: usize| x[i + 1] - x[i];
    let (lower, rest) = scratch.split_at_mut(n);
    let (diag, rest) = rest.split_at_mut(n);
    let (upper, per_width) = rest.split_at_mut(n);

    // The inner rows. The secants of the first two and the last two intervals
    // are kept for the end rows.
    let mut before = secant(x, y, 0, per_width)?;
    let mut first = [before; 2];
    let mut tail = [before; 2];
    for i in 1..last {
        let d = secant(x, y, i, per_width)?;
        let (h_before, h) = (width(i - 1), width(i));
        lower[i] = h;
        diag[i] = two * (h_before + h);
        upper[i] = h_before;
        slopes[i] = V::zero().axpy(three, combine(h, before, h_before, d));
        if i == 1 {
            first[1] = d;
        }
        tail = [before, d];
        before = d;
    }

    // Two points have no second interval for not-a-knot to join the first
    // with: SciPy gives such an end the slope of the line through them.
    let line = |end| match end {
        EndCondition::NotAKnot if n == 2 => EndCondition::FirstDerivative(first[0]),
        other => other,
    };
    let (start, end) = (line(start), line(end));
    let not_a_knot = |end| matches!(end, EndCondition::NotAKnot);
    if n == 3 && not_a_knot(start) && not_a_knot(end) {
        // The two conditions say the same thing, that the two pieces are one
        // polynomial: the parabola through the three points. Its middle row is
        // the inner row above.
        diag[0] = T::one();
        upper[0] = T::one();
        slopes[0] = V::zero().axpy(two, first[0]);
        lower[2] = T::one();
        diag[2] = T::one();
        slopes[2] = V::zero().axpy(two, first[1]);
    } else {
        match start {
            EndCondition::NotAKnot => {
                let (h0, h1, d) = (width(0), width(1), x[2] - x[0]);
                diag[0] = h1;
                upper[0] = d;
                let sum = combine((h0 + two * d) * h1, first[0], h0 * h0, first[1]);
                slopes[0] = V::zero().axpy(d.recip(), sum);
            }
            EndCondition::FirstDerivative(v) => {
                diag[0] = T::one();
                slopes[0] = v;
            }
            EndCondition::SecondDerivative(v) => {
                let h = width(0);
                diag[0] = two * h;
                upper[0] = h;
                slopes[0] = combine(-half * h * h, v, three, y[1].axpy(-T::one(), y[0]));
            }
        }
        match end {
            EndCondition::NotAKnot => {
                let (h0, h1, d) = (width(last - 2), width(last - 1), x[last] - x[last - 2]);
                diag[last] = h0;
                lower[last] = d;
                let sum = combine(h1 * h1, tail[0], (two * d + h1) * h0, tail[1]);
                slopes[last] = V::zero().axpy(d.recip(), sum);
            }
            EndCondition::FirstDerivative(v) => {
                diag[last] = T::one();
                slopes[last] = v;
            }
            EndCondition::SecondDerivative(v) => {
                let h = width(last - 1);
                diag[last] = two * h;
                lower[last] = h;
                let rise = y[last].axpy(-T::one(), y[last - 1]);
                slopes[last] = combine(half * h * h, v, three, rise);
            }
        }
    }
    tridiagonal::solve(lower, diag, upper, slopes).map_err(|b| overflow(b, n))
}

impl<T: Scalar, V: Value<T>> Interpolator for CubicSpline<T, V> {
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
