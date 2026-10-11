//! Piecewise polynomials: the representation the 1D splines evaluate through.

use crate::{search, validate, Extrapolate, Hint, InterpolateError, Interpolator, Scalar, Value};

/// A piecewise polynomial in one variable.
///
/// SciPy's `PPoly`. Between the break points `x[i]` and `x[i + 1]` the value is
///
/// ```text
/// c[i][0] (x - x[i])^k + c[i][1] (x - x[i])^(k - 1) + ... + c[i][k]
/// ```
///
/// with `k` the degree. The coefficients are stored interval by interval, the
/// highest power first, so that one interval's coefficients are adjacent and
/// are read front to back. SciPy stores the same numbers power by power
/// (`c[power][interval]`).
///
/// An interval contains its left break point and not its right one; the last
/// interval contains both. Outside the break points the result follows the
/// [`Extrapolate`] mode, which is [`Extrapolate::Extend`] unless
/// [`with_extrapolate`](PPoly::with_extrapolate) changes it.
///
/// Cloning copies the break points and the coefficients.
///
/// # Examples
///
/// ```
/// use interpx::{Extrapolate, InterpolateError, Interpolator, PPoly};
///
/// # fn main() -> Result<(), InterpolateError> {
/// // Two linear pieces: x on [0, 1], then 1 + 2(x - 1) on [1, 2].
/// let p = PPoly::new(&[0.0, 1.0, 2.0], &[1.0, 0.0, 2.0, 1.0], 1)?;
/// assert_eq!(p.eval(&0.5), Some(0.5));
/// assert_eq!(p.eval(&1.5), Some(2.0));
///
/// // By default the end pieces continue past the ends.
/// assert_eq!(p.eval(&3.0), Some(5.0));
///
/// let p = p.with_extrapolate(Extrapolate::Forbid);
/// assert_eq!(p.eval(&3.0), None);
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone)]
pub struct PPoly<T: Scalar, V: Value<T> = T> {
    breakpoints: Vec<T>,
    coefficients: Vec<V>,
    /// Coefficients per interval: the degree plus one.
    order: usize,
    extrapolate: Extrapolate,
}

impl<T: Scalar, V: Value<T>> PPoly<T, V> {
    /// Builds a piecewise polynomial of the given degree from copies of the break
    /// points and the coefficients.
    ///
    /// `coefficients` holds `degree + 1` numbers for each of the
    /// `breakpoints.len() - 1` intervals, interval by interval, the highest power
    /// first.
    ///
    /// # Errors
    ///
    /// The first of these that applies:
    ///
    /// - [`InterpolateError::LengthMismatch`] when `coefficients` does not hold
    ///   `degree + 1` numbers per interval.
    /// - [`InterpolateError::InsufficientPoints`] when there are fewer than two
    ///   break points.
    /// - [`InterpolateError::NonFinite`] when a break point, or else a coefficient,
    ///   is NaN or infinite. Its `index` is the position in `breakpoints` when a
    ///   break point is at fault; only when every break point is finite is it a
    ///   position in `coefficients`.
    /// - [`InterpolateError::DuplicatePoints`] or [`InterpolateError::NotSorted`]
    ///   when the break points do not increase strictly.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{InterpolateError, Interpolator, PPoly};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// // (x - 1)^2 on [1, 3]: one interval, degree 2.
    /// let p = PPoly::new(&[1.0, 3.0], &[1.0, 0.0, 0.0], 2)?;
    /// assert_eq!(p.eval(&2.0), Some(1.0));
    ///
    /// // Two outputs per point: the value type is an array.
    /// let p = PPoly::new(&[0.0, 1.0], &[[1.0, -1.0], [0.0, 5.0]], 1)?;
    /// assert_eq!(p.eval(&0.5), Some([0.5, 4.5]));
    ///
    /// let err = PPoly::new(&[0.0, 1.0], &[1.0, 0.0, 0.0], 1);
    /// assert!(matches!(err, Err(InterpolateError::LengthMismatch { expected: 2, actual: 3 })));
    /// # Ok(())
    /// # }
    /// ```
    pub fn new(
        breakpoints: &[T],
        coefficients: &[V],
        degree: usize,
    ) -> Result<Self, InterpolateError> {
        let order = check(breakpoints, coefficients, degree)?;
        Ok(Self {
            breakpoints: breakpoints.to_vec(),
            coefficients: coefficients.to_vec(),
            order,
            extrapolate: Extrapolate::Extend,
        })
    }

    /// Builds a piecewise polynomial of the given degree, taking the break points
    /// and the coefficients without copying them.
    ///
    /// The layout and the checks are those of [`PPoly::new`].
    ///
    /// # Errors
    ///
    /// The same as [`PPoly::new`].
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{InterpolateError, Interpolator, PPoly};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let p = PPoly::from_owned(vec![0.0, 2.0], vec![3.0, 1.0], 1)?;
    /// assert_eq!(p.eval(&1.0), Some(4.0));
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_owned(
        breakpoints: Vec<T>,
        coefficients: Vec<V>,
        degree: usize,
    ) -> Result<Self, InterpolateError> {
        let order = check(&breakpoints, &coefficients, degree)?;
        Ok(Self {
            breakpoints,
            coefficients,
            order,
            extrapolate: Extrapolate::Extend,
        })
    }

    /// Returns the same polynomial with another behavior outside its domain.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{Extrapolate, InterpolateError, Interpolator, PPoly};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let p = PPoly::new(&[0.0, 1.0], &[2.0, 0.0], 1)?; // 2x
    /// assert_eq!(p.clone().with_extrapolate(Extrapolate::Clamp).eval(&5.0), Some(2.0));
    /// assert_eq!(p.clone().with_extrapolate(Extrapolate::Periodic).eval(&1.25), Some(0.5));
    /// assert_eq!(p.with_extrapolate(Extrapolate::Forbid).eval(&5.0), None);
    /// # Ok(())
    /// # }
    /// ```
    #[must_use]
    pub fn with_extrapolate(mut self, extrapolate: Extrapolate) -> Self {
        self.extrapolate = extrapolate;
        self
    }

    /// Returns the break points, in increasing order.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{InterpolateError, PPoly};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let p = PPoly::new(&[0.0, 1.0, 4.0], &[1.0, 2.0], 0)?;
    /// assert_eq!(p.breakpoints(), [0.0, 1.0, 4.0]);
    /// assert_eq!(p.coefficients(), [1.0, 2.0]);
    /// assert_eq!(p.degree(), 0);
    /// # Ok(())
    /// # }
    /// ```
    pub fn breakpoints(&self) -> &[T] {
        &self.breakpoints
    }

    /// Returns the coefficients: `degree + 1` per interval, interval by interval,
    /// the highest power first. See [`breakpoints`](PPoly::breakpoints) for an
    /// example.
    pub fn coefficients(&self) -> &[V] {
        &self.coefficients
    }

    /// Returns the degree of the pieces. See [`breakpoints`](PPoly::breakpoints)
    /// for an example.
    pub fn degree(&self) -> usize {
        self.order - 1
    }

    /// Returns the behavior outside the domain.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{Extrapolate, InterpolateError, PPoly};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let p = PPoly::new(&[0.0, 1.0], &[1.0], 0)?;
    /// assert_eq!(p.extrapolate(), Extrapolate::Extend);
    /// # Ok(())
    /// # }
    /// ```
    pub fn extrapolate(&self) -> Extrapolate {
        self.extrapolate
    }

    /// Evaluates the polynomial at `x`, looking first at the interval in `hint`
    /// and its neighbors.
    ///
    /// The result is bit for bit that of [`eval`](Interpolator::eval). `hint` is
    /// updated to the interval used, so a run of nearby queries skips the binary
    /// search.
    ///
    /// Use this for queries that come in order. Each hinted search starts from the
    /// interval of the query before it, so it waits for that one. For queries in
    /// no order the wait is a whole binary search, and plain searches, which do
    /// not wait for each other, are about twice as fast: call `eval` instead.
    ///
    /// [`eval_batch_into`](Interpolator::eval_batch_into) makes this choice
    /// itself. It takes the points 64 at a time and uses a hint only on a run that
    /// is in order (rising or falling) and crosses no more intervals than it has
    /// points.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{Hint, InterpolateError, PPoly};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let p = PPoly::new(&[0.0, 1.0, 2.0], &[1.0, 0.0, 2.0, 1.0], 1)?;
    /// let mut hint = Hint::new();
    /// assert_eq!(p.eval_hinted(1.5, &mut hint), Some(2.0));
    /// assert_eq!(p.eval_hinted(f64::NAN, &mut hint), None);
    /// # Ok(())
    /// # }
    /// ```
    pub fn eval_hinted(&self, x: T, hint: &mut Hint) -> Option<V> {
        search::eval_hinted(&self.breakpoints, self.extrapolate, x, hint, |i, x| {
            self.horner(i, x)
        })
    }

    /// Evaluates the polynomial of interval `i` at `x` by Horner's rule.
    fn horner(&self, i: usize, x: T) -> V {
        let dx = x - self.breakpoints[i];
        let piece = &self.coefficients[i * self.order..(i + 1) * self.order];
        piece.iter().fold(V::zero(), |acc, c| c.axpy(dx, acc))
    }
}

/// Runs the input checks of [`PPoly::new`] and returns the number of
/// coefficients per interval.
fn check<T: Scalar, V: Value<T>>(
    breakpoints: &[T],
    coefficients: &[V],
    degree: usize,
) -> Result<usize, InterpolateError> {
    let intervals = breakpoints.len().saturating_sub(1);
    // A degree too large to count cannot match any slice.
    let order = degree.saturating_add(1);
    let expected = intervals.saturating_mul(order);
    validate::coordinates_and_values(breakpoints, coefficients, expected, 2)?;
    Ok(order)
}

impl<T: Scalar, V: Value<T>> Interpolator for PPoly<T, V> {
    type Point = T;
    type Output = V;

    fn eval(&self, x: &T) -> Option<V> {
        search::eval(&self.breakpoints, self.extrapolate, *x, |i, x| {
            self.horner(i, x)
        })
    }

    fn eval_batch_into(
        &self,
        points: &[T],
        out: &mut [V],
        fill: V,
    ) -> Result<(), InterpolateError> {
        search::eval_batch_into(
            &self.breakpoints,
            self.extrapolate,
            points,
            out,
            fill,
            |i, x| self.horner(i, x),
        )
    }
}
