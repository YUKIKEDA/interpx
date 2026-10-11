//! Piecewise polynomials: the representation the 1D splines evaluate through.

use crate::scalar::{count, falling};
use crate::{
    search, validate, Extrapolate, Hint, InterpolateError, Interpolator, ParamReason, Scalar, Value,
};

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
/// [`with_extrapolate`](PPoly::with_extrapolate) changes it. Under
/// [`Extrapolate::Periodic`] a period holds its left end only: the last break
/// point is the first one again, as in SciPy.
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

    /// Returns the derivative of order `nu` as a piecewise polynomial.
    ///
    /// SciPy's `PPoly.derivative`. The break points stay and the degree drops by
    /// `nu`. An order above the degree gives the polynomial of degree 0 that is
    /// zero everywhere, and order 0 gives a copy.
    ///
    /// The derivative keeps the extrapolation mode, except that
    /// [`Extrapolate::Clamp`] becomes [`Extrapolate::Forbid`]: a clamped function
    /// is constant outside its domain, and no mode gives the zero its derivative
    /// is there. [`eval_derivative`](PPoly::eval_derivative) does return it.
    ///
    /// # Errors
    ///
    /// Returns [`InterpolateError::NonFinite`] when a coefficient overflows.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{InterpolateError, Interpolator, PPoly};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// // x^3 on [0, 2].
    /// let cube = PPoly::new(&[0.0, 2.0], &[1.0, 0.0, 0.0, 0.0], 3)?;
    /// let slope = cube.derivative(1)?; // 3 x^2
    /// assert_eq!(slope.degree(), 2);
    /// assert_eq!(slope.coefficients(), [3.0, 0.0, 0.0]);
    /// assert_eq!(slope.eval(&2.0), Some(12.0));
    ///
    /// assert_eq!(cube.derivative(3)?.coefficients(), [6.0]);
    /// assert_eq!(cube.derivative(4)?.coefficients(), [0.0]);
    /// # Ok(())
    /// # }
    /// ```
    pub fn derivative(&self, nu: usize) -> Result<Self, InterpolateError> {
        if nu == 0 {
            return Ok(self.clone());
        }
        let extrapolate = self.extrapolate.of_derivative();
        let degree = self.degree();
        if nu > degree {
            let zeros = vec![V::zero(); self.breakpoints.len() - 1];
            return Ok(
                Self::from_owned(self.breakpoints.clone(), zeros, 0)?.with_extrapolate(extrapolate)
            );
        }
        let kept = self.order - nu;
        let mut coefficients = Vec::with_capacity(self.coefficients.len() / self.order * kept);
        for piece in self.coefficients.chunks_exact(self.order) {
            // The power of the first coefficient, going down by one each step.
            let mut power = count::<T>(degree);
            for c in &piece[..kept] {
                coefficients.push(V::zero().axpy(falling(power, nu), *c));
                power = power - T::one();
            }
        }
        Ok(
            Self::from_owned(self.breakpoints.clone(), coefficients, degree - nu)?
                .with_extrapolate(extrapolate),
        )
    }

    /// Returns the antiderivative of order `nu` as a piecewise polynomial.
    ///
    /// SciPy's `PPoly.antiderivative`. The break points stay and the degree rises
    /// by `nu`. The antiderivative is continuous, and it is zero at the first
    /// break point, as are its derivatives below order `nu`. Order 0 gives a
    /// copy.
    ///
    /// [`Extrapolate::Extend`] is kept. Every other mode becomes
    /// [`Extrapolate::Forbid`]: the antiderivative of a periodic function is not
    /// periodic in general (SciPy drops the mode as well), and that of a clamped
    /// function is linear outside the domain. [`integrate`](PPoly::integrate)
    /// follows the mode of this polynomial.
    ///
    /// # Errors
    ///
    /// Returns [`InterpolateError::NonFinite`] when a coefficient overflows, and
    /// [`InterpolateError::InvalidParameter`] when `nu` is so large that the
    /// number of coefficients cannot be counted.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{InterpolateError, Interpolator, PPoly};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// // 2 on [0, 1], then 4 on [1, 3].
    /// let steps = PPoly::new(&[0.0, 1.0, 3.0], &[2.0, 4.0], 0)?;
    /// let area = steps.antiderivative(1)?;
    /// // 2x, then 2 + 4(x - 1): continuous at 1.
    /// assert_eq!(area.coefficients(), [2.0, 0.0, 4.0, 2.0]);
    /// assert_eq!(area.eval(&3.0), Some(10.0));
    /// # Ok(())
    /// # }
    /// ```
    pub fn antiderivative(&self, nu: usize) -> Result<Self, InterpolateError> {
        if nu == 0 {
            return Ok(self.clone());
        }
        let intervals = self.breakpoints.len() - 1;
        let Some(new_order) = self
            .order
            .checked_add(nu)
            .filter(|order| intervals.checked_mul(*order).is_some())
        else {
            return Err(InterpolateError::InvalidParameter {
                name: "nu",
                reason: ParamReason::TooLarge {
                    max: usize::MAX / intervals - self.order,
                },
            });
        };
        // Each piece keeps its coefficients at the front, and gets one more
        // behind them at every step.
        let mut coefficients = vec![V::zero(); intervals * new_order];
        for (new, old) in coefficients
            .chunks_exact_mut(new_order)
            .zip(self.coefficients.chunks_exact(self.order))
        {
            new[..self.order].copy_from_slice(old);
        }
        // The highest power after the step.
        let mut top = count::<T>(self.order);
        for step in 0..nu {
            let filled = self.order + step;
            // The integral from the first break point to the left end of the
            // piece: its constant, which makes the result continuous.
            let mut constant = V::zero();
            for (piece, ends) in coefficients
                .chunks_exact_mut(new_order)
                .zip(self.breakpoints.windows(2))
            {
                // The powers after this step: `filled` down to 1.
                let mut power = top;
                for c in &mut piece[..filled] {
                    *c = V::zero().axpy(power.recip(), *c);
                    power = power - T::one();
                }
                piece[filled] = constant;
                let width = ends[1] - ends[0];
                constant = piece[..=filled]
                    .iter()
                    .fold(V::zero(), |acc, c| c.axpy(width, acc));
            }
            top = top + T::one();
        }
        Ok(
            Self::from_owned(self.breakpoints.clone(), coefficients, new_order - 1)?
                .with_extrapolate(self.extrapolate.of_antiderivative()),
        )
    }

    /// Returns the derivative of order `nu` at `x`.
    ///
    /// SciPy's `PPoly.__call__(x, nu)`. This is the derivative of the function
    /// [`eval`](Interpolator::eval) gives, piece by piece:
    ///
    /// - Order 0 is `eval`, bit for bit.
    /// - `x` belongs to the same interval as in `eval`. Where the derivative
    ///   jumps at a break point, the result is that of the piece on the right;
    ///   at the last break point, of the piece on the left.
    /// - Outside the domain, [`Extrapolate::Extend`] differentiates the end piece
    ///   continued, [`Extrapolate::Periodic`] wraps `x`, and
    ///   [`Extrapolate::Forbid`] gives `None`. Under [`Extrapolate::Clamp`] the
    ///   function is constant there, so an order of 1 or more gives zero.
    /// - An order above the degree gives zero wherever there is a value.
    ///
    /// This does not allocate. For many points, build
    /// [`derivative`](PPoly::derivative) once and call
    /// [`eval_batch_into`](Interpolator::eval_batch_into) on it.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{Extrapolate, InterpolateError, PPoly};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// // x on [0, 1], then 1 - (x - 1) on [1, 2]: the slope jumps at 1.
    /// let hat = PPoly::new(&[0.0, 1.0, 2.0], &[1.0, 0.0, -1.0, 1.0], 1)?;
    /// assert_eq!(hat.eval_derivative(0.5, 1), Some(1.0));
    /// assert_eq!(hat.eval_derivative(1.0, 1), Some(-1.0));
    /// assert_eq!(hat.eval_derivative(2.0, 1), Some(-1.0));
    /// assert_eq!(hat.eval_derivative(0.5, 2), Some(0.0));
    ///
    /// let hat = hat.with_extrapolate(Extrapolate::Clamp);
    /// assert_eq!(hat.eval_derivative(5.0, 1), Some(0.0));
    /// # Ok(())
    /// # }
    /// ```
    pub fn eval_derivative(&self, x: T, nu: usize) -> Option<V> {
        if nu == 0 {
            return self.eval(&x);
        }
        search::eval_derivative(&self.breakpoints, self.extrapolate, x, V::zero(), |i, x| {
            self.horner_derivative(i, x, nu)
        })
    }

    /// Returns the integral from `a` to `b`.
    ///
    /// SciPy's `PPoly.integrate`. This is the integral of the function
    /// [`eval`](Interpolator::eval) gives, so outside the domain it follows the
    /// [`Extrapolate`] mode:
    ///
    /// - [`Extrapolate::Extend`] integrates the end pieces continued.
    /// - [`Extrapolate::Clamp`] integrates the end values, which are constant.
    /// - [`Extrapolate::Periodic`] integrates the function repeated.
    /// - [`Extrapolate::Forbid`] gives `None` when `a` or `b` is outside the
    ///   domain. (SciPy returns NaN.)
    ///
    /// With `a > b` the result is minus the integral from `b` to `a`. A limit
    /// that is NaN or infinite gives `None`. (SciPy raises for a NaN.)
    ///
    /// This does not allocate. Its cost grows with the number of intervals
    /// between `a` and `b`.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{Extrapolate, InterpolateError, PPoly};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// // x on [0, 1], then 1 on [1, 2].
    /// let p = PPoly::new(&[0.0, 1.0, 2.0], &[1.0, 0.0, 0.0, 1.0], 1)?;
    /// assert_eq!(p.integrate(0.0, 2.0), Some(1.5));
    /// assert_eq!(p.integrate(2.0, 0.0), Some(-1.5));
    /// // The last piece continues past 2.
    /// assert_eq!(p.integrate(1.0, 4.0), Some(3.0));
    ///
    /// let p = p.with_extrapolate(Extrapolate::Forbid);
    /// assert_eq!(p.integrate(1.0, 4.0), None);
    /// assert_eq!(p.integrate(f64::NAN, 1.0), None);
    /// # Ok(())
    /// # }
    /// ```
    pub fn integrate(&self, a: T, b: T) -> Option<V> {
        search::integrate(
            &self.breakpoints,
            self.extrapolate,
            a,
            b,
            |i, x| self.horner(i, x),
            |i, a, j, b| self.integral_between(i, a, j, b),
        )
    }

    /// The coefficients of interval `i`.
    fn piece(&self, i: usize) -> &[V] {
        &self.coefficients[i * self.order..(i + 1) * self.order]
    }

    /// Evaluates the polynomial of interval `i` at `x` by Horner's rule.
    fn horner(&self, i: usize, x: T) -> V {
        let dx = x - self.breakpoints[i];
        self.piece(i)
            .iter()
            .fold(V::zero(), |acc, c| c.axpy(dx, acc))
    }

    /// Evaluates the derivative of order `nu` of the polynomial of interval `i`
    /// at `x`: Horner's rule on the coefficients the derivative would have.
    fn horner_derivative(&self, i: usize, x: T, nu: usize) -> V {
        let Some(kept) = self.order.checked_sub(nu) else {
            return V::zero();
        };
        let dx = x - self.breakpoints[i];
        let mut power = count::<T>(self.order - 1);
        let mut acc = V::zero();
        for c in &self.piece(i)[..kept] {
            acc = V::zero().axpy(falling(power, nu), *c).axpy(dx, acc);
            power = power - T::one();
        }
        acc
    }

    /// The integral from `a` to `b`, with `a` taken to be in interval `i` and `b`
    /// in interval `j`, `i <= j`.
    fn integral_between(&self, i: usize, a: T, j: usize, b: T) -> V {
        let mut sum = self.primitive(j, b).axpy(-T::one(), self.primitive(i, a));
        // The whole intervals from `i` up to, and not including, `j`.
        for whole in i..j {
            sum = sum.axpy(T::one(), self.primitive(whole, self.breakpoints[whole + 1]));
        }
        sum
    }

    /// The integral of the polynomial of interval `i` from its left break point
    /// to `x`.
    fn primitive(&self, i: usize, x: T) -> V {
        let dx = x - self.breakpoints[i];
        // The powers after integrating: the order down to 1.
        let mut power = count::<T>(self.order);
        let mut acc = V::zero();
        for c in self.piece(i) {
            acc = V::zero().axpy(power.recip(), *c).axpy(dx, acc);
            power = power - T::one();
        }
        V::zero().axpy(dx, acc)
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
