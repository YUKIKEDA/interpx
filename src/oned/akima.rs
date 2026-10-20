//! Akima and modified Akima interpolation of data in one variable.

use super::hermite::{self, secant};
use crate::scalar::count;
use crate::{Extrapolate, Hint, InterpolateError, Interpolator, PPoly, Scalar, Value};

/// How an [`Akima1D`] weighs the secants around a data point.
///
/// # Examples
///
/// ```
/// use interpx::{Akima1D, AkimaMethod, InterpolateError, Interpolator};
///
/// # fn main() -> Result<(), InterpolateError> {
/// // A flat run, then a straight rise.
/// let x = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
/// let y = [1.0, 1.0, 1.0, 2.0, 3.0, 4.0];
/// // Makima keeps the flat run flat up to its last point.
/// let makima = Akima1D::new(&x, &y, AkimaMethod::Makima)?;
/// assert_eq!(makima.eval(&1.5), Some(1.0));
/// // Akima's weights vanish at that point (the secants are equal on each side
/// // of it), its slope there is the mean of the two, and the curve dips under
/// // the flat run.
/// let akima = Akima1D::new(&x, &y, AkimaMethod::Akima)?;
/// assert!(akima.eval(&1.5).is_some_and(|v: f64| v < 1.0));
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[non_exhaustive]
pub enum AkimaMethod {
    /// Akima's weights: how much the secant changes on each side.
    Akima,
    /// The modified weights, which also grow with the size of the secants. A
    /// flat run of three points or more stays flat, where [`AkimaMethod::Akima`]
    /// can overshoot into it.
    Makima,
}

/// The Akima interpolant through data points: piecewise cubic with a continuous
/// first derivative, each slope taken from the five points around it.
///
/// Because a slope depends on nearby points only, an outlier bends the curve next to it and nowhere else, and the curve
/// does not ripple as a spline can. [`AkimaMethod`] chooses between Akima's
/// weights and the modified ones (makima).
///
/// The slope at a point is a weighted mean of the secants before and after it.
/// The weights look one secant further on each side, so past each end two
/// secants are made up by continuing the last two in a straight line. Where
/// the weights vanish (the two secants before the point are equal, and so are
/// the two after it) the slope is the mean of the outer two secants. Two
/// points give the line through them.
///
/// Outside the data the result follows the [`Extrapolate`] mode. It starts as
/// [`Extrapolate::Forbid`].
///
/// The interpolant is a piecewise polynomial: [`as_ppoly`](Akima1D::as_ppoly)
/// lends it, and its derivatives and integrals come from there.
///
/// With array values, each column is interpolated on its own: the result for a
/// column is, bit for bit, what that column gives alone, and a column of large
/// values leaves the others as they are.
///
/// # Examples
///
/// ```
/// use interpx::{Akima1D, AkimaMethod, Extrapolate, InterpolateError, Interpolator};
///
/// # fn main() -> Result<(), InterpolateError> {
/// let x = [0.0, 1.0, 2.0, 3.0, 4.0];
/// let y = [0.0, 1.0, 4.0, 9.0, 16.0];
/// let akima = Akima1D::new(&x, &y, AkimaMethod::Akima)?;
/// assert_eq!(akima.eval(&2.0), Some(4.0));
/// assert!(akima.eval(&2.5).is_some_and(|v: f64| (v - 6.25).abs() < 1e-12));
///
/// // No value outside the data, unless asked for.
/// assert_eq!(akima.eval(&5.0), None);
/// let akima = akima.with_extrapolate(Extrapolate::Extend);
/// assert!(akima.eval(&5.0).is_some());
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone)]
pub struct Akima1D<T: Scalar, V: Value<T> = T> {
    poly: PPoly<T, V>,
    method: AkimaMethod,
}

impl<T: Scalar, V: Value<T>> Akima1D<T, V> {
    /// Builds the interpolant through `(x[i], y[i])` with the weights of
    /// `method`.
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
    /// - [`InterpolateError::NonFinite`] when a secant or a coefficient
    ///   overflows. Its `index` is that of a data point near the overflow.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{Akima1D, AkimaMethod, InterpolateError, Interpolator};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// // Data on a line is reproduced, by either method.
    /// let x = [0.0, 1.0, 3.0, 4.0];
    /// let y = [1.0, 3.0, 7.0, 9.0];
    /// for method in [AkimaMethod::Akima, AkimaMethod::Makima] {
    ///     let line = Akima1D::new(&x, &y, method)?;
    ///     assert_eq!(line.eval(&2.0), Some(5.0));
    ///     assert_eq!(line.method(), method);
    /// }
    ///
    /// let err = Akima1D::new(&[0.0], &[1.0], AkimaMethod::Akima);
    /// assert!(matches!(err, Err(InterpolateError::InsufficientPoints { required: 2, actual: 1 })));
    /// # Ok(())
    /// # }
    /// ```
    pub fn new(x: &[T], y: &[V], method: AkimaMethod) -> Result<Self, InterpolateError> {
        crate::validate::coordinates_and_values(x, y, x.len(), 2)?;
        let n = x.len();
        let intervals = n - 1;
        // The slopes go where the coefficients will be (`hermite::into_ppoly`).
        // Until then the secants are kept behind them in the same buffer, with
        // two made-up ones before the first and two after the last
        // (`continue_past_the_ends`): `n + 3` in all. Four entries per interval
        // are room for both from four points on; three points take one more.
        let mut buffer = vec![V::zero(); (4 * intervals).max(2 * n + 3)];
        let mut per_width = vec![T::zero(); intervals];
        for i in 0..intervals {
            buffer[n + 2 + i] = secant(x, y, i, &mut per_width)?;
        }
        if n == 2 {
            // The line through the two points.
            buffer[0] = buffer[n + 2];
            buffer[1] = buffer[n + 2];
        } else {
            let (slopes, secants) = buffer.split_at_mut(n);
            let secants = &mut secants[..n + 3];
            continue_past_the_ends(secants);
            write_slopes(secants, method, slopes);
        }
        buffer.truncate(4 * intervals);
        Ok(Self {
            poly: hermite::into_ppoly(x, y, buffer, &per_width)?
                .with_extrapolate(Extrapolate::Forbid),
            method,
        })
    }

    /// Returns the same interpolant with another behavior outside the data.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{Akima1D, AkimaMethod, Extrapolate, InterpolateError, Interpolator};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let akima = Akima1D::new(&[0.0, 1.0, 2.0], &[0.0, 1.0, 4.0], AkimaMethod::Makima)?;
    /// assert_eq!(akima.extrapolate(), Extrapolate::Forbid);
    /// assert_eq!(akima.eval(&3.0), None);
    ///
    /// let akima = akima.with_extrapolate(Extrapolate::Clamp);
    /// assert_eq!(akima.eval(&3.0), Some(4.0));
    /// assert_eq!(akima.domain(), (0.0, 2.0));
    /// # Ok(())
    /// # }
    /// ```
    #[must_use]
    pub fn with_extrapolate(mut self, extrapolate: Extrapolate) -> Self {
        self.poly = self.poly.with_extrapolate(extrapolate);
        self
    }

    /// Returns the behavior outside the data. See
    /// [`with_extrapolate`](Akima1D::with_extrapolate) for an example.
    pub fn extrapolate(&self) -> Extrapolate {
        self.poly.extrapolate()
    }

    /// Returns the weights the slopes were chosen with. See [`Akima1D::new`] for
    /// an example.
    pub fn method(&self) -> AkimaMethod {
        self.method
    }

    /// Returns the first and the last coordinate of the data. See
    /// [`with_extrapolate`](Akima1D::with_extrapolate) for an example.
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
    /// use interpx::{Akima1D, AkimaMethod, InterpolateError};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// // The line 1 + 2x: its slope is 2 everywhere, its integral over [0, 4] is 20.
    /// let line = Akima1D::new(&[0.0, 1.0, 3.0, 4.0], &[1.0, 3.0, 7.0, 9.0], AkimaMethod::Akima)?;
    /// let poly = line.as_ppoly();
    /// assert_eq!(poly.degree(), 3);
    /// assert_eq!(poly.eval_derivative(2.0, 1), Some(2.0));
    /// assert_eq!(poly.integrate(0.0, 4.0), Some(20.0));
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
    /// use interpx::{Akima1D, AkimaMethod, Hint, InterpolateError, Interpolator};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let akima = Akima1D::new(&[0.0, 1.0, 2.0], &[0.0, 1.0, 4.0], AkimaMethod::Akima)?;
    /// let mut hint = Hint::new();
    /// for x in [0.25, 0.5, 1.75] {
    ///     assert_eq!(akima.eval_hinted(x, &mut hint), akima.eval(&x));
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn eval_hinted(&self, x: T, hint: &mut Hint) -> Option<V> {
        self.poly.eval_hinted(x, hint)
    }
}

/// `2 a - b`: the secant after `a`, when `b`, `a`, and it lie on a line.
fn continued<T: Scalar, V: Value<T>>(a: V, b: V) -> V {
    V::zero().axpy(T::one() + T::one(), a).axpy(-T::one(), b)
}

/// The weight of two neighboring secants: how much the secant changes there,
/// and for makima half the size of their sum as well.
fn weight<T: Scalar>(method: AkimaMethod, a: T, b: T) -> T {
    let change = (b - a).abs();
    match method {
        AkimaMethod::Akima => change,
        AkimaMethod::Makima => change + (a + b).abs() / (T::one() + T::one()),
    }
}

/// Fills the two secants before the first and the two after the last of
/// `secants`, by continuing the nearest two in a line, as SciPy does. The
/// secants of the intervals are in `secants[2..len - 2]`, two at least. With
/// these there are `n + 3` secants for `n` points, and point `i` has the
/// secants `i..i + 4` around it.
fn continue_past_the_ends<T: Scalar, V: Value<T>>(secants: &mut [V]) {
    let last = secants.len() - 3;
    secants[1] = continued(secants[2], secants[3]);
    secants[0] = continued(secants[1], secants[2]);
    secants[last + 1] = continued(secants[last], secants[last - 1]);
    secants[last + 2] = continued(secants[last + 1], secants[last]);
}

/// Writes the slope at every data point to `slopes`, from the secants around
/// the points ([`continue_past_the_ends`]): those of point `i` are
/// `secants[i..i + 4]`.
///
/// The weight of two neighboring secants says how much the secant changes
/// between them. A point has the weight of the two secants before it and of
/// the two after it. The secants are gone through twice: once to find, column
/// by column, the largest sum of the two weights, which sets the cutoff under
/// which a sum counts as zero, and once to write the slopes. Each time the
/// loop reads four secants next to each other and nothing it wrote before, so
/// that several points are worked on at once. (Carrying the last four secants
/// and three weights from point to point, each weight computed once, took
/// 19 ns a point; this takes a third of it. #122.)
fn write_slopes<T: Scalar, V: Value<T>>(secants: &[V], method: AkimaMethod, slopes: &mut [V]) {
    let weigh = |a: V, b: V| a.zip_with(b, |a, b| weight(method, a, b));
    let mut largest = V::zero();
    for around in secants.windows(4) {
        let (before, after) = (weigh(around[0], around[1]), weigh(around[2], around[3]));
        largest = largest.zip_with(before.axpy(T::one(), after), T::max);
    }
    // SciPy's `break_mult`: a sum of weights under this share of the largest one
    // is taken for zero.
    let cutoff = count::<T>(10).powi(-9);
    let half = (T::one() + T::one()).recip();
    for (slope, around) in slopes.iter_mut().zip(secants.windows(4)) {
        let (before, after) = (weigh(around[0], around[1]), weigh(around[2], around[3]));
        *slope = V::combine(
            [
                around[0], around[1], around[2], around[3], before, after, largest,
            ],
            |[m1, m2, m3, m4, f2, f1, largest]| {
                let f12 = f1 + f2;
                if f12 > cutoff * largest {
                    m2 + f2 / f12 * (m3 - m2)
                } else {
                    // The slope is not defined by the weights: the mean of the
                    // outer two secants.
                    half * (m4 + m1)
                }
            },
        );
    }
}

impl<T: Scalar, V: Value<T>> Interpolator for Akima1D<T, V> {
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

/// The stored form of the `serde` feature.
#[cfg(feature = "serde")]
mod stored {
    use super::{Akima1D, AkimaMethod};
    use crate::{PPoly, Scalar, Value};
    use serde::{Deserialize, Deserializer, Serialize, Serializer};
    use std::borrow::Cow;

    /// The version of the stored form.
    const VERSION: u32 = 1;

    /// What is stored: the piecewise cubic, and the method it was built with.
    #[derive(Serialize, Deserialize)]
    #[serde(bound = "")]
    struct Form<'a, T: Scalar, V: Value<T>> {
        version: u32,
        poly: Cow<'a, PPoly<T, V>>,
        method: AkimaMethod,
    }

    impl<T: Scalar, V: Value<T>> Serialize for Akima1D<T, V> {
        fn serialize<S: Serializer>(&self, to: S) -> Result<S::Ok, S::Error> {
            Form {
                version: VERSION,
                poly: Cow::Borrowed(&self.poly),
                method: self.method,
            }
            .serialize(to)
        }
    }

    impl<'de, T: Scalar, V: Value<T>> Deserialize<'de> for Akima1D<T, V> {
        fn deserialize<D: Deserializer<'de>>(from: D) -> Result<Self, D::Error> {
            let form = Form::<T, V>::deserialize(from)?;
            crate::stored::version(form.version, VERSION)?;
            crate::stored::degree(form.poly.degree(), 3)?;
            Ok(Self {
                poly: form.poly.into_owned(),
                method: form.method,
            })
        }
    }
}
