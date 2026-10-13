//! Nearest, previous, and next neighbor interpolation of data in one variable.

use crate::{search, validate, Extrapolate, Hint, InterpolateError, Interpolator, Scalar, Value};

/// Which data point a [`Nearest1D`] takes its value from.
///
/// The `kind` argument of SciPy's `interp1d`, for its four kinds that return a
/// data value as it is.
///
/// # Examples
///
/// ```
/// use interpx::{InterpolateError, Interpolator, Nearest1D, NearestKind};
///
/// # fn main() -> Result<(), InterpolateError> {
/// let (x, y) = ([0.0, 1.0, 3.0], [0.0, 10.0, 30.0]);
/// let at = |kind, q: f64| Nearest1D::new(&x, &y, kind).map(|f| f.eval(&q));
///
/// // 1.5 is nearer to 1 than to 3.
/// assert_eq!(at(NearestKind::Nearest, 1.5)?, Some(10.0));
/// assert_eq!(at(NearestKind::Previous, 1.5)?, Some(10.0));
/// assert_eq!(at(NearestKind::Next, 1.5)?, Some(30.0));
///
/// // 2 is as near to 1 as to 3.
/// assert_eq!(at(NearestKind::Nearest, 2.0)?, Some(10.0));
/// assert_eq!(at(NearestKind::NearestUp, 2.0)?, Some(30.0));
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[non_exhaustive]
pub enum NearestKind {
    /// The nearest point; halfway between two points, the one on the left.
    /// SciPy's `kind="nearest"`.
    Nearest,
    /// The nearest point; halfway between two points, the one on the right.
    /// SciPy's `kind="nearest-up"`.
    NearestUp,
    /// The last point at or before the query. SciPy's `kind="previous"`.
    Previous,
    /// The first point at or after the query. SciPy's `kind="next"`.
    Next,
}

/// Interpolation that returns the value of one data point: the nearest, the
/// previous, or the next one.
///
/// SciPy's `interp1d` with `kind` set to `"nearest"`, `"nearest-up"`,
/// `"previous"`, or `"next"`; [`NearestKind`] says which. At a data point every
/// kind returns the value of that point.
///
/// Outside the data the result follows the [`Extrapolate`] mode. It starts as
/// [`Extrapolate::Forbid`] (`interp1d` raises an error there by default).
/// [`Extrapolate::Clamp`] and [`Extrapolate::Extend`] both give the first or
/// the last value, and [`Extrapolate::Periodic`] wraps the query.
///
/// Three things differ from SciPy:
///
/// - A NaN query has no value. SciPy returns the last data value.
/// - With `fill_value="extrapolate"`, SciPy's `"previous"` returns NaN left of
///   the data and its `"next"` right of it. Here [`Extrapolate::Clamp`] and
///   [`Extrapolate::Extend`] return the end value on both sides, for every kind.
/// - The data is checked more strictly: two points at least (SciPy takes one)
///   and `x` strictly increasing (SciPy sorts, and takes repeated `x`).
///
/// # Examples
///
/// ```
/// use interpx::{Extrapolate, InterpolateError, Interpolator, Nearest1D, NearestKind};
///
/// # fn main() -> Result<(), InterpolateError> {
/// let nearest = Nearest1D::new(&[0.0, 1.0, 3.0], &[0.0, 10.0, 30.0], NearestKind::Nearest)?;
/// assert_eq!(nearest.eval(&0.4), Some(0.0));
/// assert_eq!(nearest.eval(&0.6), Some(10.0));
/// assert_eq!(nearest.eval(&3.0), Some(30.0));
///
/// assert_eq!(nearest.eval(&3.5), None);
/// let nearest = nearest.with_extrapolate(Extrapolate::Clamp);
/// assert_eq!(nearest.eval(&3.5), Some(30.0));
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone)]
pub struct Nearest1D<T: Scalar, V: Value<T> = T> {
    x: Vec<T>,
    y: Vec<V>,
    kind: NearestKind,
    extrapolate: Extrapolate,
    /// Whether the coordinates are evenly spaced (`search::evenly_spaced`): found at
    /// construction, and only a matter of how a search starts.
    even: bool,
}

impl<T: Scalar, V: Value<T>> Nearest1D<T, V> {
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
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{InterpolateError, Interpolator, Nearest1D, NearestKind};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// // Two outputs per point: the value type is an array.
    /// let steps = Nearest1D::new(&[0.0, 1.0], &[[0.0, 5.0], [1.0, 6.0]], NearestKind::Previous)?;
    /// assert_eq!(steps.eval(&0.9), Some([0.0, 5.0]));
    /// assert_eq!(steps.eval(&1.0), Some([1.0, 6.0]));
    ///
    /// let err = Nearest1D::new(&[0.0], &[1.0], NearestKind::Nearest);
    /// assert!(matches!(err, Err(InterpolateError::InsufficientPoints { required: 2, actual: 1 })));
    /// # Ok(())
    /// # }
    /// ```
    pub fn new(x: &[T], y: &[V], kind: NearestKind) -> Result<Self, InterpolateError> {
        validate::coordinates_and_values(x, y, x.len(), 2)?;
        Ok(Self::from_checked(x.to_vec(), y.to_vec(), kind))
    }

    /// The interpolator of data that passed the checks of [`Nearest1D::new`].
    fn from_checked(x: Vec<T>, y: Vec<V>, kind: NearestKind) -> Self {
        Self {
            even: search::evenly_spaced(&x),
            x,
            y,
            kind,
            extrapolate: Extrapolate::Forbid,
        }
    }

    /// Returns the same interpolator with another behavior outside the data.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{Extrapolate, InterpolateError, Interpolator, Nearest1D, NearestKind};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let next = Nearest1D::new(&[0.0, 1.0, 3.0], &[0.0, 10.0, 30.0], NearestKind::Next)?;
    /// assert_eq!(next.eval(&4.0), None);
    /// assert_eq!(next.clone().with_extrapolate(Extrapolate::Clamp).eval(&4.0), Some(30.0));
    /// assert_eq!(next.clone().with_extrapolate(Extrapolate::Extend).eval(&4.0), Some(30.0));
    /// // The range of the data repeated: 3.5 is 0.5, and the next point is 1.
    /// assert_eq!(next.with_extrapolate(Extrapolate::Periodic).eval(&3.5), Some(10.0));
    /// # Ok(())
    /// # }
    /// ```
    #[must_use]
    pub fn with_extrapolate(mut self, extrapolate: Extrapolate) -> Self {
        self.extrapolate = extrapolate;
        self
    }

    /// Returns the behavior outside the data.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{Extrapolate, InterpolateError, Nearest1D, NearestKind};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let nearest = Nearest1D::new(&[0.0, 1.0], &[0.0, 2.0], NearestKind::NearestUp)?;
    /// assert_eq!(nearest.extrapolate(), Extrapolate::Forbid);
    /// assert_eq!(nearest.kind(), NearestKind::NearestUp);
    /// assert_eq!(nearest.domain(), (0.0, 1.0));
    /// # Ok(())
    /// # }
    /// ```
    pub fn extrapolate(&self) -> Extrapolate {
        self.extrapolate
    }

    /// Returns which data point a value is taken from. See
    /// [`extrapolate`](Nearest1D::extrapolate) for an example.
    pub fn kind(&self) -> NearestKind {
        self.kind
    }

    /// Returns the first and the last coordinate of the data. See
    /// [`extrapolate`](Nearest1D::extrapolate) for an example.
    pub fn domain(&self) -> (T, T) {
        (self.x[0], self.x[self.x.len() - 1])
    }

    /// Evaluates at `x`, looking first at the interval in `hint` and its
    /// neighbors.
    ///
    /// The result is bit for bit that of [`eval`](Interpolator::eval). Use this
    /// for queries that come in order; for queries in no order `eval` is faster.
    /// [`PPoly::eval_hinted`](crate::PPoly::eval_hinted) explains why, and
    /// [`eval_batch_into`](Interpolator::eval_batch_into) makes the choice itself.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{Hint, InterpolateError, Nearest1D, NearestKind};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let nearest = Nearest1D::new(&[0.0, 1.0, 3.0], &[0.0, 10.0, 30.0], NearestKind::Nearest)?;
    /// let mut hint = Hint::new();
    /// assert_eq!(nearest.eval_hinted(0.4, &mut hint), Some(0.0));
    /// assert_eq!(nearest.eval_hinted(2.5, &mut hint), Some(30.0));
    /// assert_eq!(nearest.eval_hinted(f64::NAN, &mut hint), None);
    /// # Ok(())
    /// # }
    /// ```
    pub fn eval_hinted(&self, x: T, hint: &mut Hint) -> Option<V> {
        search::eval_hinted(&self.x, self.extrapolate, x, hint, |i, x| self.pick(i, x))
    }

    /// The value for `x`, taking `x` to be between points `i` and `i + 1`, or
    /// beyond them when it is outside the data.
    fn pick(&self, i: usize, x: T) -> V {
        let (left, right) = (self.x[i], self.x[i + 1]);
        let takes_left = match self.kind {
            NearestKind::Nearest => x <= midpoint(left, right),
            NearestKind::NearestUp => x < midpoint(left, right),
            // The search gives the right end of an interval for the last
            // coordinate only, and beyond it when the last segment is continued.
            NearestKind::Previous => x < right,
            NearestKind::Next => x <= left,
        };
        self.y[if takes_left { i } else { i + 1 }]
    }
}

/// The point halfway between two coordinates, by SciPy's formula: each is halved
/// before they are added, so the sum does not overflow, and the rounding is
/// SciPy's.
fn midpoint<T: Scalar>(left: T, right: T) -> T {
    let half = (T::one() + T::one()).recip();
    left * half + right * half
}

impl<T: Scalar, V: Value<T>> Interpolator for Nearest1D<T, V> {
    type Point = T;
    type Output = V;

    fn eval(&self, x: &T) -> Option<V> {
        search::eval(&self.x, self.even, self.extrapolate, *x, |i, x| {
            self.pick(i, x)
        })
    }

    fn eval_batch_into(
        &self,
        points: &[T],
        out: &mut [V],
        fill: V,
    ) -> Result<(), InterpolateError> {
        let even = self.even;
        search::eval_batch_into(
            &self.x,
            even,
            self.extrapolate,
            points,
            out,
            fill,
            |i, x| self.pick(i, x),
        )
    }
}

/// The stored form of the `serde` feature.
#[cfg(feature = "serde")]
mod stored {
    use super::{Nearest1D, NearestKind};
    use crate::{validate, Extrapolate, Scalar, Value};
    use serde::{de, Deserialize, Deserializer, Serialize, Serializer};
    use std::borrow::Cow;

    /// The version of the stored form.
    const VERSION: u32 = 1;

    /// What is stored: the data and the rules. Whether the coordinates are
    /// evenly spaced is found again when the form is read.
    #[derive(Serialize, Deserialize)]
    #[serde(bound = "")]
    struct Form<'a, T: Scalar, V: Value<T>> {
        version: u32,
        #[serde(with = "crate::stored")]
        x: Cow<'a, [T]>,
        #[serde(with = "crate::stored")]
        y: Cow<'a, [V]>,
        kind: NearestKind,
        extrapolate: Extrapolate,
    }

    impl<T: Scalar, V: Value<T>> Serialize for Nearest1D<T, V> {
        fn serialize<S: Serializer>(&self, to: S) -> Result<S::Ok, S::Error> {
            Form {
                version: VERSION,
                x: Cow::Borrowed(&self.x),
                y: Cow::Borrowed(&self.y),
                kind: self.kind,
                extrapolate: self.extrapolate,
            }
            .serialize(to)
        }
    }

    /// Reads what [`Nearest1D::new`] accepts: the checks are its checks.
    impl<'de, T: Scalar, V: Value<T>> Deserialize<'de> for Nearest1D<T, V> {
        fn deserialize<D: Deserializer<'de>>(from: D) -> Result<Self, D::Error> {
            let form = Form::<T, V>::deserialize(from)?;
            crate::stored::version(form.version, VERSION)?;
            let (x, y) = (form.x.into_owned(), form.y.into_owned());
            validate::coordinates_and_values(&x, &y, x.len(), 2).map_err(de::Error::custom)?;
            Ok(Self::from_checked(x, y, form.kind).with_extrapolate(form.extrapolate))
        }
    }
}
