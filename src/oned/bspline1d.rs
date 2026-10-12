//! Interpolation of data in one variable by a B-spline of a chosen degree.

use super::cubic_spline::periodic_tolerance;
use crate::banded::{self, Band, Cyclic};
use crate::bspline::{basis_derivatives, basis_values, MAX_ORDER};
use crate::tridiagonal::Breakdown;
use crate::{
    validate, BSpline, BoundaryCondition, EndCondition, Extrapolate, Hint, InterpolateError,
    Interpolator, ParamReason, Scalar, Value, BSPLINE_MAX_DEGREE,
};

/// The spline of a chosen degree through data points: piecewise polynomial of
/// degree `k`, with `k - 1` continuous derivatives.
///
/// SciPy's `make_interp_spline`. The degree is 1 to
/// [`BSPLINE_MAX_DEGREE`]: degree 1 is the polygon through the points, degree
/// 3 the cubic spline, and the higher degrees are smoother. (The step
/// function, which SciPy has as degree 0, is
/// [`Nearest1D`](crate::Nearest1D) with
/// [`NearestKind::Previous`](crate::NearestKind::Previous).)
///
/// What the data leaves open is settled by the [`BoundaryCondition`]:
///
/// - [`BoundaryCondition::NotAKnot`], for every degree: the knots are chosen
///   so that the spline is one polynomial across the data points next to each
///   end. `degree + 1` points give one polynomial of that degree.
/// - [`BoundaryCondition::Periodic`], for every degree: the spline and its
///   first `degree - 1` derivatives repeat from the last point on. The first
///   and the last value must agree.
/// - [`BoundaryCondition::Ends`] with a first or a second derivative at an
///   end, for degree 3 only (as in SciPy, where `"natural"` and `"clamped"`
///   fit no other degree). The other end may be a derivative or
///   [`EndCondition::NotAKnot`].
///
/// A spline of degree 3 is the curve of [`CubicSpline`](crate::CubicSpline)
/// for the same data and condition, to rounding; `CubicSpline` is the faster
/// to build, and this type is for the other degrees.
///
/// Outside the data the result follows the [`Extrapolate`] mode. It starts as
/// [`Extrapolate::Extend`], and as [`Extrapolate::Periodic`] for the periodic
/// condition.
///
/// The interpolant is a [`BSpline`]: [`as_bspline`](BSpline1D::as_bspline)
/// lends it, and its derivatives and integrals come from there.
///
/// # Examples
///
/// ```
/// use interpx::{BSpline1D, BoundaryCondition, InterpolateError, Interpolator};
///
/// # fn main() -> Result<(), InterpolateError> {
/// let x = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
/// let y = [0.0, 1.0, 4.0, 9.0, 16.0, 25.0];
/// // Degree 2 reproduces a parabola.
/// let spline = BSpline1D::new(&x, &y, 2, BoundaryCondition::NotAKnot)?;
/// let at = spline.eval(&2.5).unwrap_or(f64::NAN);
/// assert!((at - 6.25).abs() < 1e-12);
/// // Degree 1 is the polygon through the points.
/// let polygon = BSpline1D::new(&x, &y, 1, BoundaryCondition::NotAKnot)?;
/// assert_eq!(polygon.eval(&2.5), Some(6.5));
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone)]
pub struct BSpline1D<T: Scalar, V: Value<T> = T> {
    spline: BSpline<T, V>,
}

/// A derivative given at an end of a cubic spline: its order (1 or 2) and its
/// value.
type EndDerivative<V> = Option<(usize, V)>;

/// The derivative an end condition gives, if it gives one.
fn end_derivative<V>(end: EndCondition<V>) -> EndDerivative<V> {
    match end {
        EndCondition::NotAKnot => None,
        EndCondition::FirstDerivative(value) => Some((1, value)),
        EndCondition::SecondDerivative(value) => Some((2, value)),
    }
}

impl<T: Scalar, V: Value<T>> BSpline1D<T, V> {
    /// Builds the spline of the given degree through `(x[i], y[i])`.
    ///
    /// # Errors
    ///
    /// The first of these that applies:
    ///
    /// - [`InterpolateError::LengthMismatch`] when `x` and `y` differ in length.
    /// - [`InterpolateError::InvalidParameter`] with the name `"degree"` when
    ///   the degree is 0 ([`ParamReason::TooSmall`]) or over
    ///   [`BSPLINE_MAX_DEGREE`] ([`ParamReason::TooLarge`]).
    /// - [`InterpolateError::InvalidParameter`] with the name `"bc"` and
    ///   [`ParamReason::Unsupported`] when a derivative is given at an end and
    ///   the degree is not 3.
    /// - [`InterpolateError::InsufficientPoints`] when the points, with the
    ///   derivatives given at the ends, do not fix one polynomial of the
    ///   degree: there must be `degree + 1` points, less one for each end
    ///   with a derivative, and two at least. The periodic condition takes two
    ///   points for any degree.
    /// - [`InterpolateError::NonFinite`] when a coordinate, or else a value, is
    ///   NaN or infinite.
    /// - [`InterpolateError::DuplicatePoints`] or [`InterpolateError::NotSorted`]
    ///   when `x` does not increase strictly.
    /// - [`InterpolateError::InvalidParameter`] with the name `"bc"` and
    ///   [`ParamReason::NonFinite`] when a derivative given at an end is NaN or
    ///   infinite.
    /// - [`InterpolateError::NotPeriodic`] when the condition is periodic and
    ///   the first and the last value differ.
    /// - [`InterpolateError::DuplicatePoints`] when the degree is even and two
    ///   neighboring points are so close that the number halfway between them
    ///   is one of them: an even degree puts its knots there.
    /// - [`InterpolateError::NonFinite`] when a knot overflows, the system for
    ///   the control points cannot be solved, or a control point overflows.
    ///   Its `index` is that of a data point, a row of the system, or a
    ///   control point near the overflow.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{BSpline1D, BoundaryCondition, InterpolateError, Interpolator, ParamReason};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// // A closed curve of degree 5: two columns, periodic.
    /// let t = [0.0, 1.0, 2.0, 3.0, 4.0];
    /// let points = [[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -1.0], [1.0, 0.0]];
    /// let curve = BSpline1D::new(&t, &points, 5, BoundaryCondition::Periodic)?;
    /// assert_eq!(curve.eval(&5.0), curve.eval(&1.0));
    ///
    /// // A natural cubic spline: degree 3 takes derivatives at the ends.
    /// let x = [0.0, 1.0, 2.0, 3.0];
    /// let natural = BSpline1D::new(&x, &[0.0, 1.0, 0.0, 1.0], 3, BoundaryCondition::natural())?;
    /// let bend = natural.as_bspline().eval_derivative(0.0, 2).unwrap_or(f64::NAN);
    /// assert!(bend.abs() < 1e-12);
    ///
    /// // Another degree does not.
    /// let err = BSpline1D::new(&x, &[0.0, 1.0, 0.0, 1.0], 2, BoundaryCondition::natural());
    /// let unsupported = InterpolateError::InvalidParameter {
    ///     name: "bc",
    ///     reason: ParamReason::Unsupported,
    /// };
    /// assert_eq!(err.err(), Some(unsupported));
    /// # Ok(())
    /// # }
    /// ```
    pub fn new(
        x: &[T],
        y: &[V],
        degree: usize,
        bc: BoundaryCondition<V>,
    ) -> Result<Self, InterpolateError> {
        validate::length(y, x.len())?;
        if degree == 0 {
            return Err(InterpolateError::InvalidParameter {
                name: "degree",
                reason: ParamReason::TooSmall { min: 1 },
            });
        }
        if degree > BSPLINE_MAX_DEGREE {
            return Err(InterpolateError::InvalidParameter {
                name: "degree",
                reason: ParamReason::TooLarge {
                    max: BSPLINE_MAX_DEGREE,
                },
            });
        }
        let periodic = matches!(bc, BoundaryCondition::Periodic);
        let (left, right) = match bc {
            BoundaryCondition::Ends(left, right) => (end_derivative(left), end_derivative(right)),
            BoundaryCondition::NotAKnot | BoundaryCondition::Periodic => (None, None),
        };
        let given = usize::from(left.is_some()) + usize::from(right.is_some());
        if given > 0 && degree != 3 {
            return Err(InterpolateError::InvalidParameter {
                name: "bc",
                reason: ParamReason::Unsupported,
            });
        }
        // The conditions fix one polynomial of the degree at least. A periodic
        // spline of two points is a constant.
        let required = if periodic {
            2
        } else {
            (degree + 1 - given).max(2)
        };
        validate::coordinates_and_values(x, y, x.len(), required)?;
        for (_, value) in [left, right].into_iter().flatten() {
            if !value.all_finite() {
                return Err(InterpolateError::InvalidParameter {
                    name: "bc",
                    reason: ParamReason::NonFinite,
                });
            }
        }
        let n = x.len();
        if periodic && !y[0].within(&y[n - 1], periodic_tolerance()) {
            return Err(InterpolateError::NotPeriodic);
        }

        let spline = if degree == 1 {
            polygon(x, y)?
        } else if periodic {
            periodic_spline(x, y, degree)?
        } else {
            open_spline(x, y, degree, left, right)?
        };
        let extrapolate = if periodic {
            Extrapolate::Periodic
        } else {
            Extrapolate::Extend
        };
        Ok(Self {
            spline: spline.with_extrapolate(extrapolate),
        })
    }

    /// Returns the same interpolant with another behavior outside the data.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{BSpline1D, BoundaryCondition, Extrapolate, InterpolateError, Interpolator};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let x = [0.0, 1.0, 2.0];
    /// let spline = BSpline1D::new(&x, &[0.0, 1.0, 4.0], 2, BoundaryCondition::NotAKnot)?;
    /// assert_eq!(spline.extrapolate(), Extrapolate::Extend);
    /// assert!(spline.eval(&3.0).is_some());
    /// assert_eq!(spline.degree(), 2);
    ///
    /// let spline = spline.with_extrapolate(Extrapolate::Clamp);
    /// assert_eq!(spline.eval(&3.0), spline.eval(&2.0));
    /// assert_eq!(spline.domain(), (0.0, 2.0));
    /// # Ok(())
    /// # }
    /// ```
    #[must_use]
    pub fn with_extrapolate(mut self, extrapolate: Extrapolate) -> Self {
        self.spline = self.spline.with_extrapolate(extrapolate);
        self
    }

    /// Returns the behavior outside the data. See
    /// [`with_extrapolate`](BSpline1D::with_extrapolate) for an example.
    pub fn extrapolate(&self) -> Extrapolate {
        self.spline.extrapolate()
    }

    /// Returns the degree of the spline. See
    /// [`with_extrapolate`](BSpline1D::with_extrapolate) for an example.
    pub fn degree(&self) -> usize {
        self.spline.degree()
    }

    /// Returns the first and the last coordinate of the data. See
    /// [`with_extrapolate`](BSpline1D::with_extrapolate) for an example.
    pub fn domain(&self) -> (T, T) {
        self.spline.domain()
    }

    /// Returns the interpolant as a B-spline, with the same extrapolation mode.
    ///
    /// This is the way to its knots and control points, and to its derivatives
    /// and integrals: see [`BSpline::derivative`], [`BSpline::antiderivative`],
    /// [`BSpline::integrate`], and [`BSpline::eval_derivative`]. It evaluates to
    /// the same bits as this interpolant.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{BSpline1D, BoundaryCondition, InterpolateError};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let x = [0.0, 1.0, 2.0, 3.0];
    /// let y = [0.0, 1.0, 8.0, 27.0];
    /// // Four points and degree 3: one cubic, x^3.
    /// let spline = BSpline1D::new(&x, &y, 3, BoundaryCondition::NotAKnot)?;
    /// let b = spline.as_bspline();
    /// assert_eq!(b.degree(), 3);
    /// assert_eq!(b.knots(), &[0.0, 0.0, 0.0, 0.0, 3.0, 3.0, 3.0, 3.0]);
    /// let slope = b.eval_derivative(2.0, 1).unwrap_or(f64::NAN);
    /// assert!((slope - 12.0).abs() < 1e-12);
    /// # Ok(())
    /// # }
    /// ```
    pub fn as_bspline(&self) -> &BSpline<T, V> {
        &self.spline
    }

    /// Evaluates at `x`, looking first at the interval in `hint` and its
    /// neighbors.
    ///
    /// The result is bit for bit that of [`eval`](Interpolator::eval). Use this
    /// for queries that come in order; for queries in no order `eval` is faster.
    /// [`BSpline::eval_hinted`] explains why, and
    /// [`eval_batch_into`](Interpolator::eval_batch_into) makes the choice itself.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{BSpline1D, BoundaryCondition, Hint, InterpolateError, Interpolator};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let x = [0.0, 1.0, 2.0, 3.0];
    /// let spline = BSpline1D::new(&x, &[0.0, 1.0, 0.0, 2.0], 2, BoundaryCondition::NotAKnot)?;
    /// let mut hint = Hint::new();
    /// for x in [0.25, 0.5, 2.75] {
    ///     assert_eq!(spline.eval_hinted(x, &mut hint), spline.eval(&x));
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn eval_hinted(&self, x: T, hint: &mut Hint) -> Option<V> {
        self.spline.eval_hinted(x, hint)
    }
}

/// A breakdown of the solver as the error of the constructor.
fn unsolved(at: Breakdown) -> InterpolateError {
    InterpolateError::NonFinite { index: at.row }
}

/// Degree 1: the control points are the values, and nothing is solved.
fn polygon<T: Scalar, V: Value<T>>(x: &[T], y: &[V]) -> Result<BSpline<T, V>, InterpolateError> {
    let n = x.len();
    let mut knots = Vec::with_capacity(n + 2);
    knots.push(x[0]);
    knots.extend_from_slice(x);
    knots.push(x[n - 1]);
    BSpline::from_owned(knots, y.to_vec(), 1)
}

/// Degree 2 or more, not periodic: `left` and `right` are the derivatives given
/// at the ends (degree 3 only).
///
/// Both ends of the knots have multiplicity `degree + 1`. Between them, the
/// knots are data points (odd degrees) or midpoints of neighboring data points
/// (even degrees), without the outermost ones: that leaves as many control
/// points as there are conditions. An end with a derivative keeps one knot more
/// than a not-a-knot end.
fn open_spline<T: Scalar, V: Value<T>>(
    x: &[T],
    y: &[V],
    degree: usize,
    left: EndDerivative<V>,
    right: EndDerivative<V>,
) -> Result<BSpline<T, V>, InterpolateError> {
    let n = x.len();
    let before = usize::from(left.is_some());
    let unknowns = n + before + usize::from(right.is_some());
    let mut knots = Vec::with_capacity(unknowns + degree + 1);
    knots.resize(degree + 1, x[0]);
    if degree % 2 == 1 {
        // With a derivative at an end (degree 3) the data point next to that
        // end stays a knot.
        let skip = degree.div_ceil(2);
        let (first, last) = (skip - before, n - skip + (unknowns - n - before));
        knots.extend_from_slice(&x[first..last]);
    } else {
        let skip = degree / 2;
        let two = T::one() + T::one();
        for (i, pair) in x.windows(2).enumerate() {
            // Halved before they are added: the sum of two large coordinates
            // would overflow. The result is the same otherwise.
            let middle = between(pair[0] / two + pair[1] / two, i, pair[0], pair[1])?;
            if (skip..n - 1 - skip).contains(&i) {
                knots.push(middle);
            }
        }
    }
    knots.resize(unknowns + degree + 1, x[n - 1]);

    // The rows: a derivative at the left end, the data points, a derivative at
    // the right end. In this order the matrix has `degree` diagonals at most on
    // each side.
    //
    // Without a derivative it has one fewer. A point lies in the interval of
    // `degree + 1` splines, and only two rows reach the outermost diagonal:
    // the first, in the first interval, and the last, in the last one. There
    // the entry is the spline that is zero at that end of the interval, and it
    // is left out. (The cost of the elimination grows with the square of the
    // number of diagonals; `benches/bspline1d.rs` measures the construction.)
    let open_ends = left.is_none() && right.is_none();
    let side = if open_ends { degree - 1 } else { degree };
    let band = Band {
        n: unknowns,
        lower: side,
        upper: side,
    };
    let mut ab = vec![T::zero(); band.len()];
    let mut rhs = Vec::with_capacity(unknowns);
    // The knot interval of the point, which only moves right. The last point is
    // the right end of the last interval.
    let last_span = unknowns - 1;
    let mut span = degree;
    if let Some((order, value)) = left {
        let row = basis_derivatives(&knots, degree, span, x[0], order);
        fill(band, &mut ab, 0, span - degree, &row[..=degree]);
        rhs.push(value);
    }
    for (i, at) in x.iter().enumerate() {
        while span < last_span && knots[span + 1] <= *at {
            span += 1;
        }
        let row = basis_values::<T, MAX_ORDER>(&knots, degree, span, *at);
        if open_ends && i == 0 {
            debug_assert!(
                row[degree] == T::zero(),
                "zero at the left end of its interval"
            );
            fill(band, &mut ab, 0, 0, &row[..degree]);
        } else if open_ends && i == n - 1 {
            debug_assert!(row[0] == T::zero(), "zero at the right end of its interval");
            fill(band, &mut ab, i, span - degree + 1, &row[1..=degree]);
        } else {
            fill(band, &mut ab, before + i, span - degree, &row[..=degree]);
        }
    }
    rhs.extend_from_slice(y);
    if let Some((order, value)) = right {
        let row = basis_derivatives(&knots, degree, last_span, x[n - 1], order);
        fill(
            band,
            &mut ab,
            unknowns - 1,
            last_span - degree,
            &row[..=degree],
        );
        rhs.push(value);
    }

    let mut pivots = vec![0; unknowns];
    banded::factor(band, &mut ab, &mut pivots).map_err(unsolved)?;
    banded::solve(band, &ab, &pivots, &mut rhs);
    BSpline::from_owned(knots, rhs, degree)
}

/// A knot of an even degree, which lies between the data points `index` and
/// `index + 1`, at `low` and `high`.
///
/// # Errors
///
/// Returns [`InterpolateError::DuplicatePoints`] when the knot is not strictly
/// between the two points: they are so close that the number halfway is one of
/// them. A knot on a data point can make the system singular, so such points
/// are refused as two the method cannot tell apart. Returns
/// [`InterpolateError::NonFinite`] with the later point when the knot is not
/// finite (the distance between the two overflows).
fn between<T: Scalar>(knot: T, index: usize, low: T, high: T) -> Result<T, InterpolateError> {
    if !knot.is_finite() {
        return Err(InterpolateError::NonFinite { index: index + 1 });
    }
    if low < knot && knot < high {
        Ok(knot)
    } else {
        Err(InterpolateError::DuplicatePoints {
            first: index,
            second: index + 1,
        })
    }
}

/// Writes the entries of one row of a banded matrix: `entries[j]` is in column
/// `first_column + j`.
fn fill<T: Scalar>(band: Band, ab: &mut [T], row: usize, first_column: usize, entries: &[T]) {
    for (j, entry) in entries.iter().enumerate() {
        ab[band.at(row, first_column + j)] = *entry;
    }
}

/// Degree 2 or more, periodic.
///
/// The knots inside the data are the data points (odd degrees), or the data
/// points moved half an interval to the left, all but the two ends (even
/// degrees). `degree` more knots on each side repeat their spacings with the
/// period. On such knots the control points repeat with the period too: there
/// are `n - 1` unknowns, and the conditions at all points but the last make a
/// cyclic banded matrix.
fn periodic_spline<T: Scalar, V: Value<T>>(
    x: &[T],
    y: &[V],
    degree: usize,
) -> Result<BSpline<T, V>, InterpolateError> {
    let n = x.len();
    let unknowns = n - 1;
    let mut knots = vec![T::zero(); n + 2 * degree];
    let inside = &mut knots[degree..degree + n];
    inside.copy_from_slice(x);
    if degree.is_multiple_of(2) {
        let two = T::one() + T::one();
        for (i, (knot, pair)) in inside[1..n - 1].iter_mut().zip(x.windows(2)).enumerate() {
            *knot = between(pair[1] - (pair[1] - pair[0]) / two, i, pair[0], pair[1])?;
        }
    }
    // Spacing `i` of the knots inside the data, read from the knots themselves.
    for i in 0..degree {
        let from_the_end = degree + n - 1 - i % unknowns;
        knots[degree - 1 - i] = knots[degree - i] - (knots[from_the_end] - knots[from_the_end - 1]);
        let from_the_start = degree + i % unknowns;
        knots[degree + n + i] =
            knots[degree + n - 1 + i] + (knots[from_the_start + 1] - knots[from_the_start]);
    }
    // The knots outside the data move away from it step by step: the outermost
    // two overflow first.
    for (knot, index) in [(knots[0], 0), (knots[n + 2 * degree - 1], n - 1)] {
        if !knot.is_finite() {
            return Err(InterpolateError::NonFinite { index });
        }
    }

    // `x[i]` lies in the interval of the splines `i` to `i + degree`. For an
    // odd degree it is the left end of that interval, where the last of them is
    // zero: it is left out, and the band has `degree / 2` diagonals on each
    // side for every degree. Unknown `r` is every control point `r + side`
    // plus a multiple of the period.
    let side = degree / 2;
    let shape = Cyclic {
        n: unknowns,
        lower: side,
        upper: side,
    };
    let width = 2 * side + 1;
    let mut rows = Vec::with_capacity(unknowns * width);
    for (i, at) in x[..unknowns].iter().enumerate() {
        let row = basis_values::<T, MAX_ORDER>(&knots, degree, degree + i, *at);
        debug_assert!(
            width > degree || row[degree] == T::zero(),
            "zero at the left end of its interval"
        );
        rows.extend_from_slice(&row[..width]);
    }
    let mut solution = y[..unknowns].to_vec();
    let mut scratch = vec![T::zero(); shape.scratch_len()];
    let mut pivots = vec![0; shape.pivots_len()];
    let mut carry = vec![V::zero(); shape.carry_len()];
    shape
        .solve(&rows, &mut solution, &mut scratch, &mut pivots, &mut carry)
        .map_err(unsolved)?;

    // A multiple of the period that is at least `side`.
    let turns = unknowns * side.div_ceil(unknowns);
    let control_points = (0..n + degree - 1)
        .map(|j| solution[(j + turns - side) % unknowns])
        .collect();
    BSpline::from_owned(knots, control_points, degree)
}

impl<T: Scalar, V: Value<T>> Interpolator for BSpline1D<T, V> {
    type Point = T;
    type Output = V;

    fn eval(&self, x: &T) -> Option<V> {
        self.spline.eval(x)
    }

    fn eval_batch_into(
        &self,
        points: &[T],
        out: &mut [V],
        fill: V,
    ) -> Result<(), InterpolateError> {
        self.spline.eval_batch_into(points, out, fill)
    }
}
