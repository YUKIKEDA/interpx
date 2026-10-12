//! A parametric curve through points in several dimensions.

use crate::oned::fewest_points;
use crate::{
    validate, BSpline, BSpline1D, BoundaryCondition, Extrapolate, Hint, InterpolateError,
    Interpolator, Scalar, Value,
};

/// How [`BSplineCurve::new`] chooses the parameter of each point.
///
/// With `d` the distance between two neighboring points, the parameter grows
/// from one point to the next by the amount below; the parameters are then
/// divided by the last one, so that they run from 0 to 1.
///
/// The enum is `#[non_exhaustive]`: later versions may add ways.
///
/// # Examples
///
/// ```
/// use interpx::{BSplineCurve, BoundaryCondition, InterpolateError, Parametrization};
///
/// # fn main() -> Result<(), InterpolateError> {
/// // Three points on a line, 1 and then 4 apart.
/// let points = [[0.0_f64, 0.0], [1.0, 0.0], [5.0, 0.0]];
/// let parameters = |way| {
///     BSplineCurve::new(&points, 1, BoundaryCondition::NotAKnot, way)
///         .map(|curve| curve.parameters().to_vec())
/// };
/// assert_eq!(parameters(Parametrization::Uniform)?, [0.0, 0.5, 1.0]);
/// assert_eq!(parameters(Parametrization::ChordLength)?, [0.0, 0.2, 1.0]);
/// // The square roots of the distances: 1 and 2.
/// let centripetal = parameters(Parametrization::Centripetal)?;
/// assert!((centripetal[1] - 1.0 / 3.0).abs() < 1e-15);
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum Parametrization {
    /// By 1: the points are equally spaced in the parameter, however far apart
    /// they lie. A curve through unevenly spaced points then wiggles.
    Uniform,
    /// By `d`: the parameter follows the length of the polygon through the
    /// points, which approximates the length along the curve. The choice of
    /// FITPACK's `splprep`.
    ChordLength,
    /// By the square root of `d`: between the two above, with less overshoot
    /// at sharp turns.
    Centripetal,
}

/// The spline curve through points in `D` dimensions: `u` to `[T; D]`,
/// piecewise polynomial of a chosen degree in every coordinate.
///
/// A curve is a [`BSpline1D`] whose values are points, and every point needs a
/// parameter. [`new`](BSplineCurve::new) chooses the parameters from the points
/// by a [`Parametrization`], from 0 at the first point to 1 at the last;
/// [`with_parameters`](BSplineCurve::with_parameters) takes them from the
/// caller, as SciPy's `make_interp_spline` does with a two-dimensional `y`.
/// [`parameters`](BSplineCurve::parameters) returns them.
///
/// The degree and the [`BoundaryCondition`] are those of [`BSpline1D`]: degrees
/// 1 to [`BSPLINE_MAX_DEGREE`](crate::BSPLINE_MAX_DEGREE), not-a-knot or
/// periodic for every degree, derivatives at the ends for degree 3. A closed
/// curve is periodic, and repeats its first point as its last.
///
/// Outside the parameters of the data the result follows the [`Extrapolate`]
/// mode: [`Extrapolate::Extend`] at first, and [`Extrapolate::Periodic`] for
/// the periodic condition.
///
/// The curve is a [`BSpline`]: [`as_bspline`](BSplineCurve::as_bspline) lends
/// it. Its first derivative is the tangent.
///
/// # Examples
///
/// ```
/// use interpx::{
///     BSplineCurve, BoundaryCondition, InterpolateError, Interpolator, Parametrization,
/// };
///
/// # fn main() -> Result<(), InterpolateError> {
/// // A closed curve through four points of the unit circle.
/// let points = [[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -1.0], [1.0, 0.0]];
/// let curve = BSplineCurve::new(
///     &points,
///     3,
///     BoundaryCondition::Periodic,
///     Parametrization::ChordLength,
/// )?;
/// assert_eq!(curve.parameters(), &[0.0, 0.25, 0.5, 0.75, 1.0]);
/// // Halfway between the first two points, near the circle.
/// let [x, y] = curve.eval(&0.125).unwrap_or([f64::NAN; 2]);
/// assert!((x - y).abs() < 1e-12 && (x.hypot(y) - 1.0).abs() < 0.05);
/// // It goes around again.
/// assert_eq!(curve.eval(&1.25), curve.eval(&0.25));
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone)]
pub struct BSplineCurve<T: Scalar, const D: usize> {
    spline: BSpline1D<T, [T; D]>,
    parameters: Vec<T>,
}

impl<T: Scalar, const D: usize> BSplineCurve<T, D> {
    /// Builds the curve of the given degree through `points`, with parameters
    /// chosen by `parametrization`.
    ///
    /// # Errors
    ///
    /// The first of these that applies:
    ///
    /// - [`InterpolateError::InvalidParameter`] for a degree or a boundary
    ///   condition [`BSpline1D::new`] does not take.
    /// - [`InterpolateError::InsufficientPoints`] when there are fewer points
    ///   than [`BSpline1D::new`] needs.
    /// - [`InterpolateError::NonFinite`] when a coordinate of a point is NaN or
    ///   infinite. Its `index` is that of the point.
    /// - [`InterpolateError::DuplicatePoints`] when two neighboring points are
    ///   equal and the parametrization goes by distance
    ///   ([`Parametrization::ChordLength`], [`Parametrization::Centripetal`]),
    ///   or when they are so close, against the length of the whole, that
    ///   their parameters round to the same number.
    /// - [`InterpolateError::NonFinite`] when the points are so far apart that
    ///   a distance, or the sum of the distances, overflows.
    /// - The remaining errors of [`BSpline1D::new`]:
    ///   [`InterpolateError::NotPeriodic`] when the condition is periodic and
    ///   the first and the last point differ, and the ones after it.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{
    ///     BSplineCurve, BoundaryCondition, InterpolateError, Interpolator, Parametrization,
    /// };
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let points = [[0.0, 0.0, 0.0], [1.0, 0.0, 1.0], [1.0, 1.0, 2.0], [0.0, 1.0, 3.0]];
    /// let way = Parametrization::Centripetal;
    /// let curve = BSplineCurve::new(&points, 3, BoundaryCondition::NotAKnot, way)?;
    /// // The curve is at each point at the parameter of that point.
    /// for (u, point) in curve.parameters().iter().zip(&points) {
    ///     let at = curve.eval(u).unwrap_or([f64::NAN; 3]);
    ///     assert!(at.iter().zip(point).all(|(a, p)| (a - p).abs() < 1e-12));
    /// }
    ///
    /// // The same point twice in a row has no distance to go by.
    /// let stalled = [[0.0, 0.0], [1.0, 0.0], [1.0, 0.0], [2.0, 1.0]];
    /// let err = BSplineCurve::new(&stalled, 1, BoundaryCondition::NotAKnot, way);
    /// let duplicate = InterpolateError::DuplicatePoints { first: 1, second: 2 };
    /// assert_eq!(err.err(), Some(duplicate));
    /// # Ok(())
    /// # }
    /// ```
    pub fn new(
        points: &[[T; D]],
        degree: usize,
        bc: BoundaryCondition<[T; D]>,
        parametrization: Parametrization,
    ) -> Result<Self, InterpolateError> {
        let required = fewest_points(degree, &bc)?;
        validate::enough(points, required)?;
        if let Some(index) = points.iter().position(|point| !point.all_finite()) {
            return Err(InterpolateError::NonFinite { index });
        }
        let parameters = parameters_of(points, parametrization)?;
        let spline = BSpline1D::new(&parameters, points, degree, bc)?;
        Ok(Self { spline, parameters })
    }

    /// Builds the curve of the given degree that is at `points[i]` at the
    /// parameter `t[i]`.
    ///
    /// SciPy's `make_interp_spline(t, points)`. The result is, bit for bit, the
    /// [`BSpline1D`] of the same arguments.
    ///
    /// # Errors
    ///
    /// Those of [`BSpline1D::new`], with `t` for its coordinates and `points`
    /// for its values.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{BSplineCurve, BoundaryCondition, InterpolateError, Interpolator};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// // Positions at given times.
    /// let times = [0.0, 1.0, 3.0];
    /// let positions = [[0.0, 0.0], [2.0, 1.0], [6.0, 5.0]];
    /// let bc = BoundaryCondition::NotAKnot;
    /// let path = BSplineCurve::with_parameters(&times, &positions, 1, bc)?;
    /// assert_eq!(path.eval(&2.0), Some([4.0, 3.0]));
    /// assert_eq!(path.parameters(), &times);
    /// assert_eq!(path.domain(), (0.0, 3.0));
    /// # Ok(())
    /// # }
    /// ```
    pub fn with_parameters(
        t: &[T],
        points: &[[T; D]],
        degree: usize,
        bc: BoundaryCondition<[T; D]>,
    ) -> Result<Self, InterpolateError> {
        let spline = BSpline1D::new(t, points, degree, bc)?;
        Ok(Self {
            spline,
            parameters: t.to_vec(),
        })
    }

    /// Returns the same curve with another behavior outside the parameters of
    /// the data.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{
    ///     BSplineCurve, BoundaryCondition, Extrapolate, InterpolateError, Interpolator,
    ///     Parametrization,
    /// };
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let points = [[0.0, 0.0], [1.0, 2.0], [3.0, 2.0]];
    /// let way = Parametrization::Uniform;
    /// let curve = BSplineCurve::new(&points, 2, BoundaryCondition::NotAKnot, way)?;
    /// assert_eq!(curve.extrapolate(), Extrapolate::Extend);
    /// assert_eq!(curve.degree(), 2);
    /// assert_eq!(curve.domain(), (0.0, 1.0));
    /// assert!(curve.eval(&1.5).is_some());
    ///
    /// let curve = curve.with_extrapolate(Extrapolate::Forbid);
    /// assert_eq!(curve.eval(&1.5), None);
    /// # Ok(())
    /// # }
    /// ```
    #[must_use]
    pub fn with_extrapolate(mut self, extrapolate: Extrapolate) -> Self {
        self.spline = self.spline.with_extrapolate(extrapolate);
        self
    }

    /// Returns the behavior outside the parameters of the data. See
    /// [`with_extrapolate`](BSplineCurve::with_extrapolate) for an example.
    pub fn extrapolate(&self) -> Extrapolate {
        self.spline.extrapolate()
    }

    /// Returns the degree of the curve. See
    /// [`with_extrapolate`](BSplineCurve::with_extrapolate) for an example.
    pub fn degree(&self) -> usize {
        self.spline.degree()
    }

    /// Returns the parameters of the first and of the last point. See
    /// [`with_extrapolate`](BSplineCurve::with_extrapolate) for an example.
    pub fn domain(&self) -> (T, T) {
        self.spline.domain()
    }

    /// Returns the parameter of each point: the curve is at point `i` at
    /// `parameters()[i]`. See [`Parametrization`] for an example.
    pub fn parameters(&self) -> &[T] {
        &self.parameters
    }

    /// Returns the curve as a B-spline with points for control points, with the
    /// same extrapolation mode.
    ///
    /// This is the way to its knots and control points, and to its derivatives:
    /// [`BSpline::eval_derivative`] of order 1 is the tangent (its length the
    /// speed at which the parameter moves the point), of order 2 what the
    /// curvature is computed from. The B-spline evaluates to the same bits as
    /// the curve.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{BSplineCurve, BoundaryCondition, InterpolateError, Parametrization};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// // A straight path, twice as fast in x as in y.
    /// let points = [[0.0, 1.0], [2.0, 2.0], [4.0, 3.0]];
    /// let way = Parametrization::ChordLength;
    /// let curve = BSplineCurve::new(&points, 2, BoundaryCondition::NotAKnot, way)?;
    /// let [dx, dy] = curve.as_bspline().eval_derivative(0.3, 1).unwrap_or([f64::NAN; 2]);
    /// assert!((dx - 4.0).abs() < 1e-12 && (dy - 2.0).abs() < 1e-12);
    /// # Ok(())
    /// # }
    /// ```
    pub fn as_bspline(&self) -> &BSpline<T, [T; D]> {
        self.spline.as_bspline()
    }

    /// Evaluates at `u`, looking first at the interval in `hint` and its
    /// neighbors.
    ///
    /// The result is bit for bit that of [`eval`](Interpolator::eval). Use this
    /// for parameters that come in order, as when a curve is drawn; for
    /// parameters in no order `eval` is faster. [`BSpline::eval_hinted`]
    /// explains why, and [`eval_batch_into`](Interpolator::eval_batch_into)
    /// makes the choice itself.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{
    ///     BSplineCurve, BoundaryCondition, Hint, InterpolateError, Interpolator,
    ///     Parametrization,
    /// };
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let points = [[0.0, 0.0], [1.0, 2.0], [3.0, 2.0], [4.0, 0.0]];
    /// let way = Parametrization::Centripetal;
    /// let curve = BSplineCurve::new(&points, 3, BoundaryCondition::NotAKnot, way)?;
    /// let mut hint = Hint::new();
    /// for u in [0.1, 0.2, 0.9] {
    ///     assert_eq!(curve.eval_hinted(u, &mut hint), curve.eval(&u));
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn eval_hinted(&self, u: T, hint: &mut Hint) -> Option<[T; D]> {
        self.spline.eval_hinted(u, hint)
    }
}

/// The distance between two points.
///
/// The differences are divided by the largest of them before they are squared:
/// the squares of coordinates around `1e200` would overflow, and those of
/// differences around `1e-170` would vanish, which would make two different
/// points look equal. So the distance is zero for equal points only. It is
/// infinite when a difference itself overflows.
fn distance<T: Scalar, const D: usize>(from: &[T; D], to: &[T; D]) -> T {
    let mut differences = [T::zero(); D];
    let mut largest = T::zero();
    for ((difference, a), b) in differences.iter_mut().zip(from).zip(to) {
        *difference = *b - *a;
        largest = largest.max(difference.abs());
    }
    if largest == T::zero() || largest.is_infinite() {
        return largest;
    }
    let squares = differences.iter().fold(T::zero(), |sum, difference| {
        let scaled = *difference / largest;
        sum + scaled * scaled
    });
    largest * squares.sqrt()
}

/// The parameter of each point, from 0 to 1. The points are finite, and two at
/// least.
///
/// # Errors
///
/// Returns [`InterpolateError::DuplicatePoints`] when two neighboring points are
/// equal and the parametrization goes by distance. A sum that overflows leaves
/// parameters that are not finite, for the checks of the spline to report.
fn parameters_of<T: Scalar, const D: usize>(
    points: &[[T; D]],
    parametrization: Parametrization,
) -> Result<Vec<T>, InterpolateError> {
    let mut parameters = Vec::with_capacity(points.len());
    let mut sum = T::zero();
    parameters.push(sum);
    for (i, pair) in points.windows(2).enumerate() {
        let step = match parametrization {
            Parametrization::Uniform => T::one(),
            Parametrization::ChordLength | Parametrization::Centripetal => {
                let between = distance(&pair[0], &pair[1]);
                if between == T::zero() {
                    return Err(InterpolateError::DuplicatePoints {
                        first: i,
                        second: i + 1,
                    });
                }
                if parametrization == Parametrization::Centripetal {
                    between.sqrt()
                } else {
                    between
                }
            }
        };
        sum = sum + step;
        parameters.push(sum);
    }
    for parameter in &mut parameters {
        *parameter = *parameter / sum;
    }
    Ok(parameters)
}

impl<T: Scalar, const D: usize> Interpolator for BSplineCurve<T, D> {
    type Point = T;
    type Output = [T; D];

    fn eval(&self, u: &T) -> Option<[T; D]> {
        self.spline.eval(u)
    }

    fn eval_batch_into(
        &self,
        points: &[T],
        out: &mut [[T; D]],
        fill: [T; D],
    ) -> Result<(), InterpolateError> {
        self.spline.eval_batch_into(points, out, fill)
    }
}
