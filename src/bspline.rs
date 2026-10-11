//! B-splines, and the conversions between B-splines and piecewise polynomials.

use crate::{
    search, Extrapolate, Hint, InterpolateError, Interpolator, PPoly, ParamReason, Scalar, Value,
};

/// The largest degree a [`BSpline`] can have: 15.
///
/// Evaluation keeps `degree + 1` basis values in an array on the stack, so that
/// it does not allocate. SciPy's `make_interp_spline` is used with degrees up to
/// 5. A larger degree is rejected with [`InterpolateError::InvalidParameter`].
///
/// # Examples
///
/// ```
/// use interpx::{BSpline, BSPLINE_MAX_DEGREE};
///
/// assert_eq!(BSPLINE_MAX_DEGREE, 15);
///
/// // One piece of degree 15 on [0, 1] is accepted, degree 16 is not.
/// let knots = |k: usize| [vec![0.0; k + 1], vec![1.0; k + 1]].concat();
/// assert!(BSpline::new(&knots(15), &[1.0; 16], 15).is_ok());
/// assert!(BSpline::new(&knots(16), &[1.0; 17], 16).is_err());
/// ```
pub const BSPLINE_MAX_DEGREE: usize = 15;

/// The same limit under a short name, for this module.
const MAX_DEGREE: usize = BSPLINE_MAX_DEGREE;

/// The length of the basis arrays: one value per basis function that is not
/// zero on an interval.
const MAX_ORDER: usize = MAX_DEGREE + 1;

/// A spline in the B-spline basis.
///
/// SciPy's `BSpline`. With `n` control points `c` and degree `k`, the knots are
/// `t[0] <= t[1] <= ... <= t[n + k]` and the value is the sum of `c[i] B[i](x)`,
/// where `B[i]` is the B-spline basis function of degree `k` on the knots
/// `t[i]..=t[i + k + 1]`.
///
/// The domain is `[t[k], t[n]]`. Knots may repeat. A point belongs to the knot
/// interval that has it at its left end or inside; the right end of the domain
/// belongs to the last interval of positive width. Outside the domain the result
/// follows the [`Extrapolate`] mode, which is [`Extrapolate::Extend`] unless
/// [`with_extrapolate`](BSpline::with_extrapolate) changes it.
///
/// The degree is at most [`BSPLINE_MAX_DEGREE`].
///
/// Cloning copies the knots and the control points.
///
/// # Examples
///
/// ```
/// use interpx::{BSpline, Extrapolate, InterpolateError, Interpolator};
///
/// # fn main() -> Result<(), InterpolateError> {
/// // Degree 1: the polygon through (0, 0), (1, 1), (2, 0).
/// let hat = BSpline::new(&[0.0, 0.0, 1.0, 2.0, 2.0], &[0.0, 1.0, 0.0], 1)?;
/// assert_eq!(hat.eval(&0.5), Some(0.5));
/// assert_eq!(hat.eval(&1.0), Some(1.0));
/// assert_eq!(hat.eval(&1.5), Some(0.5));
///
/// // By default the end pieces continue past the ends.
/// assert_eq!(hat.eval(&3.0), Some(-1.0));
///
/// let hat = hat.with_extrapolate(Extrapolate::Forbid);
/// assert_eq!(hat.eval(&3.0), None);
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone)]
pub struct BSpline<T: Scalar, V: Value<T> = T> {
    knots: Vec<T>,
    control_points: Vec<V>,
    degree: usize,
    /// The first knot interval of positive width in the domain.
    first: usize,
    /// The last knot interval of positive width in the domain.
    last: usize,
    extrapolate: Extrapolate,
}

impl<T: Scalar, V: Value<T>> BSpline<T, V> {
    /// Builds a B-spline of the given degree from copies of the knots and the
    /// control points.
    ///
    /// # Errors
    ///
    /// The first of these that applies:
    ///
    /// - [`InterpolateError::InvalidParameter`] when `degree` is over
    ///   [`BSPLINE_MAX_DEGREE`].
    /// - [`InterpolateError::LengthMismatch`] when there are not
    ///   `control_points.len() + degree + 1` knots.
    /// - [`InterpolateError::InsufficientPoints`] when there are fewer than
    ///   `degree + 1` control points.
    /// - [`InterpolateError::NonFinite`] when a knot, or else a control point, is
    ///   NaN or infinite. Its `index` is the position in `knots` when a knot is
    ///   at fault; only when every knot is finite is it a position in
    ///   `control_points`.
    /// - [`InterpolateError::NotSorted`] when a knot is smaller than the one
    ///   before it. Equal knots are allowed.
    /// - [`InterpolateError::DuplicatePoints`] when the domain is empty: the knots
    ///   at `degree` and at `control_points.len()` are equal.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{BSpline, InterpolateError, Interpolator};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// // One quadratic piece on [0, 1] in Bernstein form: x^2.
    /// let square = BSpline::new(&[0.0, 0.0, 0.0, 1.0, 1.0, 1.0], &[0.0, 0.0, 1.0], 2)?;
    /// assert_eq!(square.eval(&0.5), Some(0.25));
    ///
    /// // Two outputs per point: the value type is an array.
    /// let line = BSpline::new(&[0.0, 0.0, 1.0, 1.0], &[[0.0, 1.0], [2.0, 5.0]], 1)?;
    /// assert_eq!(line.eval(&0.5), Some([1.0, 3.0]));
    ///
    /// let err = BSpline::new(&[0.0, 0.0, 1.0], &[1.0, 2.0], 1);
    /// assert!(matches!(err, Err(InterpolateError::LengthMismatch { expected: 4, actual: 3 })));
    /// # Ok(())
    /// # }
    /// ```
    pub fn new(knots: &[T], control_points: &[V], degree: usize) -> Result<Self, InterpolateError> {
        let (first, last) = check(knots, control_points, degree)?;
        Ok(Self {
            knots: knots.to_vec(),
            control_points: control_points.to_vec(),
            degree,
            first,
            last,
            extrapolate: Extrapolate::Extend,
        })
    }

    /// Builds a B-spline of the given degree, taking the knots and the control
    /// points without copying them.
    ///
    /// # Errors
    ///
    /// The same as [`BSpline::new`].
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{BSpline, InterpolateError, Interpolator};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let line = BSpline::from_owned(vec![0.0, 0.0, 2.0, 2.0], vec![1.0, 5.0], 1)?;
    /// assert_eq!(line.eval(&1.0), Some(3.0));
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_owned(
        knots: Vec<T>,
        control_points: Vec<V>,
        degree: usize,
    ) -> Result<Self, InterpolateError> {
        let (first, last) = check(&knots, &control_points, degree)?;
        Ok(Self {
            knots,
            control_points,
            degree,
            first,
            last,
            extrapolate: Extrapolate::Extend,
        })
    }

    /// Builds the B-spline that equals a piecewise polynomial.
    ///
    /// Every break point becomes a knot repeated `degree + 1` times, and every
    /// piece gets its own `degree + 1` control points (its Bézier points). No
    /// continuity between the pieces is assumed, so a piecewise polynomial that
    /// jumps at a break point is represented exactly, too. The extrapolation mode
    /// is kept.
    ///
    /// # Errors
    ///
    /// Returns [`InterpolateError::InvalidParameter`] when the degree of `poly` is
    /// over [`BSPLINE_MAX_DEGREE`], and [`InterpolateError::NonFinite`] when a
    /// control point overflows.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{BSpline, InterpolateError, Interpolator, PPoly};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// // x on [0, 1], then 1 + 2(x - 1) on [1, 2].
    /// let poly = PPoly::new(&[0.0, 1.0, 2.0], &[1.0, 0.0, 2.0, 1.0], 1)?;
    /// let spline = BSpline::from_ppoly(&poly)?;
    /// assert_eq!(spline.knots(), [0.0, 0.0, 1.0, 1.0, 2.0, 2.0]);
    /// assert_eq!(spline.control_points(), [0.0, 1.0, 1.0, 3.0]);
    /// assert_eq!(spline.eval(&1.5), poly.eval(&1.5));
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_ppoly(poly: &PPoly<T, V>) -> Result<Self, InterpolateError> {
        let degree = poly.degree();
        if degree > MAX_DEGREE {
            return Err(too_large());
        }
        let order = degree + 1;
        let breakpoints = poly.breakpoints();
        // The integers 0..=degree as scalars, by counting.
        let mut ints = [T::zero(); MAX_ORDER];
        for i in 1..order {
            ints[i] = ints[i - 1] + T::one();
        }

        let knots = breakpoints
            .iter()
            .flat_map(|b| std::iter::repeat_n(*b, order))
            .collect();
        let mut control_points = Vec::with_capacity(poly.coefficients().len());
        for (piece, ends) in poly
            .coefficients()
            .chunks_exact(order)
            .zip(breakpoints.windows(2))
        {
            let width = ends[1] - ends[0];
            // With a[m] the coefficient of (x - left)^m, Bézier point j is the sum
            // over m <= j of C(j, m) / C(degree, m) * width^m * a[m].
            for j in 0..order {
                let mut point = piece[degree];
                let mut weight = T::one();
                for m in 1..=j {
                    weight =
                        weight * (ints[j] - ints[m - 1]) / (ints[degree] - ints[m - 1]) * width;
                    point = point.axpy(weight, piece[degree - m]);
                }
                control_points.push(point);
            }
        }
        Ok(Self::from_owned(knots, control_points, degree)?.with_extrapolate(poly.extrapolate()))
    }

    /// Returns the same spline with another behavior outside its domain.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{BSpline, Extrapolate, InterpolateError, Interpolator};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let line = BSpline::new(&[0.0, 0.0, 1.0, 1.0], &[0.0, 2.0], 1)?; // 2x on [0, 1]
    /// assert_eq!(line.clone().with_extrapolate(Extrapolate::Clamp).eval(&5.0), Some(2.0));
    /// assert_eq!(line.clone().with_extrapolate(Extrapolate::Periodic).eval(&1.25), Some(0.5));
    /// assert_eq!(line.with_extrapolate(Extrapolate::Forbid).eval(&5.0), None);
    /// # Ok(())
    /// # }
    /// ```
    #[must_use]
    pub fn with_extrapolate(mut self, extrapolate: Extrapolate) -> Self {
        self.extrapolate = extrapolate;
        self
    }

    /// Returns the knots, in non-decreasing order.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{BSpline, InterpolateError};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let hat = BSpline::new(&[0.0, 0.0, 1.0, 2.0, 2.0], &[0.0, 1.0, 0.0], 1)?;
    /// assert_eq!(hat.knots(), [0.0, 0.0, 1.0, 2.0, 2.0]);
    /// assert_eq!(hat.control_points(), [0.0, 1.0, 0.0]);
    /// assert_eq!(hat.degree(), 1);
    /// assert_eq!(hat.domain(), (0.0, 2.0));
    /// # Ok(())
    /// # }
    /// ```
    pub fn knots(&self) -> &[T] {
        &self.knots
    }

    /// Returns the control points: the coefficients of the basis functions. See
    /// [`knots`](BSpline::knots) for an example.
    pub fn control_points(&self) -> &[V] {
        &self.control_points
    }

    /// Returns the degree. See [`knots`](BSpline::knots) for an example.
    pub fn degree(&self) -> usize {
        self.degree
    }

    /// Returns the two ends of the domain: the knots at `degree` and at the number
    /// of control points. See [`knots`](BSpline::knots) for an example.
    pub fn domain(&self) -> (T, T) {
        (self.knots[self.first], self.knots[self.last + 1])
    }

    /// Returns the behavior outside the domain.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{BSpline, Extrapolate, InterpolateError};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let line = BSpline::new(&[0.0, 0.0, 1.0, 1.0], &[0.0, 2.0], 1)?;
    /// assert_eq!(line.extrapolate(), Extrapolate::Extend);
    /// # Ok(())
    /// # }
    /// ```
    pub fn extrapolate(&self) -> Extrapolate {
        self.extrapolate
    }

    /// Evaluates the spline at `x`, looking first at the knot interval in `hint`
    /// and its neighbors.
    ///
    /// The result is bit for bit that of [`eval`](Interpolator::eval). Use this
    /// for queries that come in order; for queries in no order `eval` is faster.
    /// [`PPoly::eval_hinted`] explains why, and
    /// [`eval_batch_into`](Interpolator::eval_batch_into) makes the choice itself.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{BSpline, Hint, InterpolateError};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let hat = BSpline::new(&[0.0, 0.0, 1.0, 2.0, 2.0], &[0.0, 1.0, 0.0], 1)?;
    /// let mut hint = Hint::new();
    /// assert_eq!(hat.eval_hinted(0.5, &mut hint), Some(0.5));
    /// assert_eq!(hat.eval_hinted(1.5, &mut hint), Some(0.5));
    /// assert_eq!(hat.eval_hinted(f64::NAN, &mut hint), None);
    /// # Ok(())
    /// # }
    /// ```
    pub fn eval_hinted(&self, x: T, hint: &mut Hint) -> Option<V> {
        search::eval_hinted(self.grid(), self.extrapolate, x, hint, |i, x| {
            self.combine(self.first + i, x, 0)
        })
    }

    /// The knots of the domain, from the first interval of positive width to the
    /// end of the last one: the grid the interval search works on.
    fn grid(&self) -> &[T] {
        &self.knots[self.first..=self.last + 1]
    }

    /// The derivative of order `derivative` (0 for the value) at `x`, taking `x`
    /// to be in knot interval `span`, or on its polynomial continued.
    fn combine(&self, span: usize, x: T, derivative: usize) -> V {
        let basis = if derivative == 0 {
            basis_values(&self.knots, self.degree, span, x)
        } else {
            basis_derivatives(&self.knots, self.degree, span, x, derivative)
        };
        let points = &self.control_points[span - self.degree..=span];
        basis
            .iter()
            .zip(points)
            .fold(V::zero(), |sum, (b, c)| sum.axpy(*b, *c))
    }
}

impl<T: Scalar, V: Value<T>> PPoly<T, V> {
    /// Builds the piecewise polynomial that equals a B-spline.
    ///
    /// SciPy's `PPoly.from_spline`. The break points are the distinct knots of the
    /// domain of `spline`, and the coefficients of a piece are the Taylor
    /// coefficients of the spline at the left end of the piece. SciPy keeps every
    /// knot as a break point, so repeated knots give it intervals of zero width;
    /// here those intervals are left out. The extrapolation mode is kept.
    ///
    /// # Errors
    ///
    /// Returns [`InterpolateError::NonFinite`] when a coefficient overflows.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{BSpline, InterpolateError, Interpolator, PPoly};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// // The polygon through (0, 0), (1, 1), (2, 0).
    /// let hat = BSpline::new(&[0.0, 0.0, 1.0, 2.0, 2.0], &[0.0, 1.0, 0.0], 1)?;
    /// let poly = PPoly::from_bspline(&hat)?;
    /// assert_eq!(poly.breakpoints(), [0.0, 1.0, 2.0]);
    /// // x on [0, 1], then 1 - (x - 1) on [1, 2].
    /// assert_eq!(poly.coefficients(), [1.0, 0.0, -1.0, 1.0]);
    /// assert_eq!(poly.eval(&1.5), hat.eval(&1.5));
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_bspline(spline: &BSpline<T, V>) -> Result<Self, InterpolateError> {
        let knots = &spline.knots;
        let degree = spline.degree;
        let mut breakpoints = Vec::new();
        let mut coefficients = Vec::new();
        for span in spline.first..=spline.last {
            let left = knots[span];
            if left >= knots[span + 1] {
                continue;
            }
            breakpoints.push(left);
            // The coefficient of (x - left)^m is the derivative of order m at
            // `left`, over m!. The highest power comes first.
            for m in (0..=degree).rev() {
                let derivative = spline.combine(span, left, m);
                coefficients.push(V::zero().axpy(factorial::<T>(m).recip(), derivative));
            }
        }
        breakpoints.push(knots[spline.last + 1]);
        Ok(Self::from_owned(breakpoints, coefficients, degree)?
            .with_extrapolate(spline.extrapolate))
    }
}

/// `m!` as a scalar, by counting.
fn factorial<T: Scalar>(m: usize) -> T {
    let mut product = T::one();
    let mut i = T::zero();
    for _ in 0..m {
        i = i + T::one();
        product = product * i;
    }
    product
}

/// The error for a degree over [`MAX_DEGREE`].
fn too_large() -> InterpolateError {
    InterpolateError::InvalidParameter {
        name: "degree",
        reason: ParamReason::TooLarge { max: MAX_DEGREE },
    }
}

/// Runs the input checks of [`BSpline::new`] and returns the first and the last
/// knot interval of positive width in the domain.
fn check<T: Scalar, V: Value<T>>(
    knots: &[T],
    control_points: &[V],
    degree: usize,
) -> Result<(usize, usize), InterpolateError> {
    if degree > MAX_DEGREE {
        return Err(too_large());
    }
    let n = control_points.len();
    // `degree` is small, so this cannot overflow for a slice that exists.
    let expected = n + degree + 1;
    if knots.len() != expected {
        return Err(InterpolateError::LengthMismatch {
            expected,
            actual: knots.len(),
        });
    }
    if n < degree + 1 {
        return Err(InterpolateError::InsufficientPoints {
            required: degree + 1,
            actual: n,
        });
    }
    if let Some(index) = knots.iter().position(|t| !t.is_finite()) {
        return Err(InterpolateError::NonFinite { index });
    }
    if let Some(index) = control_points.iter().position(|c| !c.all_finite()) {
        return Err(InterpolateError::NonFinite { index });
    }
    if let Some(i) = knots.windows(2).position(|w| w[1] < w[0]) {
        return Err(InterpolateError::NotSorted { index: i + 1 });
    }
    let (lo, hi) = (knots[degree], knots[n]);
    if lo >= hi {
        return Err(InterpolateError::DuplicatePoints {
            first: degree,
            second: n,
        });
    }
    // The last knot equal to `lo` starts the first interval of positive width,
    // and the last knot below `hi` starts the last one.
    let first = knots.partition_point(|t| *t <= lo) - 1;
    let last = knots.partition_point(|t| *t < hi) - 1;
    Ok((first, last))
}

/// The derivatives of order `derivative` (1 or more) at `x` of the `degree + 1`
/// basis functions that are not zero on knot interval `span`. The values
/// themselves come from [`basis_values`].
///
/// Entry `j` belongs to basis function `span - degree + j`; the entries past
/// `degree` are zero. This is de Boor's recurrence as SciPy's `_deBoor_D` runs
/// it: it builds the basis functions degree by degree, and takes a derivative
/// step instead of an evaluation step for the last `derivative` degrees.
///
/// `span` is an interval of positive width and `derivative` is at most `degree`
/// (a higher order is zero, and no caller asks for it).
/// Every difference of knots the recurrence divides by then reaches across
/// `span`, so it is positive even when knots repeat: SciPy's test for equal
/// knots is not needed here.
fn basis_derivatives<T: Scalar>(
    knots: &[T],
    degree: usize,
    span: usize,
    x: T,
    derivative: usize,
) -> [T; MAX_ORDER] {
    debug_assert!(
        (1..=degree).contains(&derivative),
        "the callers ask for orders 1..=degree"
    );
    debug_assert!(knots[span] < knots[span + 1], "a span of positive width");
    let mut now = [T::zero(); MAX_ORDER];
    let mut before = [T::zero(); MAX_ORDER];
    now[0] = T::one();
    // `j` as a scalar, by counting.
    let mut j_scalar = T::zero();
    for j in 1..=degree {
        j_scalar = j_scalar + T::one();
        before[..j].copy_from_slice(&now[..j]);
        now[0] = T::zero();
        // The last `derivative` steps differentiate.
        let differentiate = j + derivative > degree;
        for n in 1..=j {
            // `left <= knots[span] < knots[span + 1] <= right`: never equal.
            let (left, right) = (knots[span + n - j], knots[span + n]);
            if differentiate {
                let w = j_scalar * before[n - 1] / (right - left);
                now[n - 1] = now[n - 1] - w;
                now[n] = w;
            } else {
                let w = before[n - 1] / (right - left);
                now[n - 1] = now[n - 1] + w * (right - x);
                now[n] = w * (x - left);
            }
        }
    }
    now
}

/// The values at `x` of the `degree + 1` basis functions that are not zero on
/// knot interval `span`. Entry `j` belongs to basis function `span - degree + j`.
///
/// This is the evaluation path, so it is the triangular form of the recurrence
/// (de Boor's; Piegl and Tiller's algorithm A2.2): it keeps the distances from
/// `x` to the knots on either side and updates the values in place.
/// [`basis_derivatives`] would give the same values for order zero, but it
/// copies its work array at every degree, which made evaluation about twice as
/// slow (measured in `benches/bspline.rs`).
///
/// `span` is an interval of positive width, so every sum of distances the
/// recurrence divides by is a difference of knots that reaches across it.
fn basis_values<T: Scalar>(knots: &[T], degree: usize, span: usize, x: T) -> [T; MAX_ORDER] {
    debug_assert!(knots[span] < knots[span + 1], "a span of positive width");
    let mut values = [T::zero(); MAX_ORDER];
    let mut left = [T::zero(); MAX_ORDER];
    let mut right = [T::zero(); MAX_ORDER];
    values[0] = T::one();
    for j in 1..=degree {
        left[j] = x - knots[span + 1 - j];
        right[j] = knots[span + j] - x;
        let mut saved = T::zero();
        for r in 0..j {
            let temp = values[r] / (right[r + 1] + left[j - r]);
            values[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        values[j] = saved;
    }
    values
}

impl<T: Scalar, V: Value<T>> Interpolator for BSpline<T, V> {
    type Point = T;
    type Output = V;

    fn eval(&self, x: &T) -> Option<V> {
        search::eval(self.grid(), self.extrapolate, *x, |i, x| {
            self.combine(self.first + i, x, 0)
        })
    }

    fn eval_batch_into(
        &self,
        points: &[T],
        out: &mut [V],
        fill: V,
    ) -> Result<(), InterpolateError> {
        search::eval_batch_into(self.grid(), self.extrapolate, points, out, fill, |i, x| {
            self.combine(self.first + i, x, 0)
        })
    }
}
