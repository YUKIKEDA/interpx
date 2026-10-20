//! The value of the nearest data point, in any number of dimensions.

use super::Scaling;
use crate::error::ParamReason;
use crate::{validate, InterpolateError, Interpolator, Value};
use kiddo::{ImmutableKdTree, SquaredEuclidean};
use std::borrow::Cow;
use std::collections::HashSet;
use std::sync::Arc;

/// The interpolant that gives a query point the value of the data point
/// nearest to it, for points anywhere in `D` dimensions.
///
/// The points need no order and no structure.
/// The distance is the Euclidean one. A value is a scalar or a fixed-length
/// array, and the result is always one of the values given: nothing is
/// averaged. Coordinates are `f64`.
///
/// The points go into a k-d tree, built once. It is meant for fewer than
/// about ten dimensions: beyond that, a search of the tree comes near to
/// looking at every point. A query of 100 000 random points took about
/// 0.14 µs in two dimensions, 0.22 µs in three, and 0.77 µs in five.
///
/// Points with the same coordinates are taken as one, with the value of the
/// one given first. Of several different points at the same distance from a
/// query, one is taken; which one is fixed for a given set of points, and
/// follows no rule.
///
/// [`rescaled`](NearestNDInterpolator::rescaled) measures the distance after
/// each coordinate is brought to the same range, for data whose coordinates
/// have different units. [`with_distance_bound`] leaves queries far from every
/// data point without a value.
///
/// A query with a coordinate that is NaN or infinite has no value.
///
/// [`with_distance_bound`]: NearestNDInterpolator::with_distance_bound
///
/// # Examples
///
/// ```
/// use interpx::{InterpolateError, Interpolator, NearestNDInterpolator};
///
/// # fn main() -> Result<(), InterpolateError> {
/// let points = [[0.0, 0.0], [1.0, 0.0], [0.0, 2.0]];
/// let values = [10.0, 20.0, 30.0];
/// let nearest = NearestNDInterpolator::new(&points, &values)?;
/// assert_eq!(nearest.eval(&[0.9, 0.3]), Some(20.0));
/// assert_eq!(nearest.eval(&[-5.0, 9.0]), Some(30.0));
/// assert_eq!(nearest.eval(&[f64::NAN, 0.0]), None);
/// assert_eq!(nearest.len(), 3);
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone)]
pub struct NearestNDInterpolator<const D: usize, V: Value<f64> = f64> {
    /// The points, each set of coordinates once, moved and scaled when the
    /// interpolant is rescaled. The tree is never changed, and clones of the
    /// interpolant share it. There is none without dimensions (`D` is 0):
    /// all points are then one point, which needs no search (and `kiddo`
    /// does not take it).
    tree: Option<Arc<ImmutableKdTree<f64, D>>>,
    /// For each point of the tree, the first of the points given that has its
    /// coordinates.
    first: Vec<u32>,
    /// The values, as given.
    values: Vec<V>,
    /// What is taken off each coordinate and what it is divided by, before a
    /// distance is measured.
    scaling: Option<Scaling<D>>,
    /// A query has a value when its nearest point is nearer than this.
    bound: Option<f64>,
}

impl<const D: usize, V: Value<f64>> NearestNDInterpolator<D, V> {
    /// Builds the interpolant from copies of the points and of their values.
    ///
    /// # Errors
    ///
    /// The first of these that applies:
    ///
    /// - [`InterpolateError::LengthMismatch`] when `points` and `values` differ
    ///   in length.
    /// - [`InterpolateError::InsufficientPoints`] when there is no point.
    /// - [`InterpolateError::InvalidParameter`] for `points` when there are
    ///   more than `u32::MAX` of them.
    /// - [`InterpolateError::NonFinite`] when a coordinate, or else a value, is
    ///   NaN or infinite. Its `index` is that of the point.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{InterpolateError, Interpolator, NearestNDInterpolator};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// // Three dimensions, two values at every point.
    /// let points = [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]];
    /// let values = [[1.0, -1.0], [2.0, -2.0]];
    /// let nearest = NearestNDInterpolator::new(&points, &values)?;
    /// assert_eq!(nearest.eval(&[0.8, 0.9, 0.4]), Some([2.0, -2.0]));
    ///
    /// // The same point twice: the value given first.
    /// let twice = NearestNDInterpolator::new(&[[0.0], [0.0], [1.0]], &[7.0, 8.0, 9.0])?;
    /// assert_eq!(twice.eval(&[0.1]), Some(7.0));
    ///
    /// let err = NearestNDInterpolator::<2>::new(&[], &[]);
    /// assert!(matches!(err, Err(InterpolateError::InsufficientPoints { required: 1, actual: 0 })));
    /// # Ok(())
    /// # }
    /// ```
    pub fn new(points: &[[f64; D]], values: &[V]) -> Result<Self, InterpolateError> {
        Self::build(points, Cow::Borrowed(values), false)
    }

    /// Builds the interpolant from values that are handed over, which are
    /// kept as they are and not copied. The points are borrowed: they go into
    /// a tree, and nothing of the slice is kept.
    ///
    /// The checks and the interpolant are those of
    /// [`new`](NearestNDInterpolator::new).
    ///
    /// # Errors
    ///
    /// As [`new`](NearestNDInterpolator::new).
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{InterpolateError, Interpolator, NearestNDInterpolator};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
    /// let values = vec![10.0, 20.0, 30.0];
    /// let nearest = NearestNDInterpolator::from_owned(&points, values)?;
    /// assert_eq!(nearest.eval(&[0.9, 0.2]), Some(20.0));
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_owned(points: &[[f64; D]], values: Vec<V>) -> Result<Self, InterpolateError> {
        Self::build(points, Cow::Owned(values), false)
    }

    /// Builds the interpolant with every coordinate brought to the same range
    /// before a distance is measured.
    ///
    /// On each axis the mean of the points is taken off and the result is
    /// divided by their range (the largest coordinate less the smallest; an
    /// axis on which all points agree is left as it is). Queries go through
    /// the same change. Use it when the coordinates have different units, so
    /// that the axis with the larger numbers does not decide alone which point
    /// is nearest.
    ///
    /// # Errors
    ///
    /// Those of [`new`](NearestNDInterpolator::new), and then
    /// [`InterpolateError::NonFinite`] when the mean or the range of the
    /// coordinates of an axis overflows. Its `index` is that of a point at
    /// which it does: the one that takes the sum of the axis past the largest
    /// number, or the one with the largest coordinate.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{InterpolateError, Interpolator, NearestNDInterpolator};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// // Positions in kilometres, heights in metres.
    /// let points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1000.0]];
    /// let values = [1.0, 2.0, 3.0];
    /// let query = [0.9, 600.0];
    /// // As the numbers stand, the height decides alone: 600 is nearer 1000
    /// // than 0, and the position does not count.
    /// let plain = NearestNDInterpolator::new(&points, &values)?;
    /// assert_eq!(plain.eval(&query), Some(3.0));
    /// // Rescaled, the query is 0.9 of the way along one axis and 0.6 along
    /// // the other: nearest the second point.
    /// let rescaled = NearestNDInterpolator::rescaled(&points, &values)?;
    /// assert_eq!(rescaled.eval(&query), Some(2.0));
    /// # Ok(())
    /// # }
    /// ```
    pub fn rescaled(points: &[[f64; D]], values: &[V]) -> Result<Self, InterpolateError> {
        Self::build(points, Cow::Borrowed(values), true)
    }

    /// Builds as [`rescaled`](NearestNDInterpolator::rescaled) does, from
    /// values that are handed over, as [`from_owned`](NearestNDInterpolator::from_owned) is
    /// to [`new`](NearestNDInterpolator::new).
    ///
    /// # Errors
    ///
    /// As [`rescaled`](NearestNDInterpolator::rescaled).
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{InterpolateError, Interpolator, NearestNDInterpolator};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// // Metres along the first axis, kilometres along the second.
    /// let points = [[0.0, 0.0], [1000.0, 0.0], [0.0, 1.0]];
    /// let scaled = NearestNDInterpolator::rescaled_from_owned(&points, vec![1.0, 2.0, 3.0])?;
    /// assert_eq!(scaled.eval(&[400.0, 0.3]), Some(1.0));
    /// # Ok(())
    /// # }
    /// ```
    pub fn rescaled_from_owned(
        points: &[[f64; D]],
        values: Vec<V>,
    ) -> Result<Self, InterpolateError> {
        Self::build(points, Cow::Owned(values), true)
    }

    fn build(
        points: &[[f64; D]],
        values: Cow<'_, [V]>,
        rescale: bool,
    ) -> Result<Self, InterpolateError> {
        validate::length(&values, points.len())?;
        validate::enough(points, 1)?;
        if u32::try_from(points.len()).is_err() {
            return Err(too_many());
        }
        let finite = |point: &[f64; D]| point.iter().all(|c| c.is_finite());
        if let Some(index) = validate::first_failing(points, finite) {
            return Err(InterpolateError::NonFinite { index });
        }
        if let Some(index) = validate::first_failing(&values, V::all_finite) {
            return Err(InterpolateError::NonFinite { index });
        }

        let first = first_of_each(points);
        let scaling = if rescale {
            Some(Scaling::of(points)?)
        } else {
            None
        };
        let kept: Vec<[f64; D]> = first
            .iter()
            .map(|i| {
                let point = &points[*i as usize];
                scaling.map_or(*point, |scaling| scaling.apply(point))
            })
            .collect();
        Ok(Self {
            tree: tree_of(&kept)?,
            first,
            values: values.into_owned(),
            scaling,
            bound: None,
        })
    }

    /// Returns the same interpolant, which gives a value only to queries whose
    /// nearest point is nearer than `bound`.
    ///
    /// A query at the distance `bound` exactly, or farther, has no value. For
    /// a [`rescaled`](NearestNDInterpolator::rescaled) interpolant the distance
    /// is the one measured after the change of coordinates. A `bound` of
    /// positive infinity is no bound, and the interpolant has none afterwards
    /// ([`distance_bound`](NearestNDInterpolator::distance_bound) returns
    /// `None`). A `bound` that is NaN, zero, or below zero (negative infinity
    /// among them) leaves no query with a value.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{InterpolateError, Interpolator, NearestNDInterpolator};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let nearest = NearestNDInterpolator::new(&[[0.0, 0.0], [10.0, 0.0]], &[1.0, 2.0])?;
    /// assert_eq!(nearest.distance_bound(), None);
    /// assert_eq!(nearest.eval(&[4.0, 0.0]), Some(1.0));
    ///
    /// let near = nearest.with_distance_bound(3.0);
    /// assert_eq!(near.distance_bound(), Some(3.0));
    /// assert_eq!(near.eval(&[2.5, 0.0]), Some(1.0));
    /// assert_eq!(near.eval(&[3.0, 0.0]), None);
    /// assert_eq!(near.eval(&[4.0, 0.0]), None);
    ///
    /// let open = near.clone().with_distance_bound(f64::INFINITY);
    /// assert_eq!(open.distance_bound(), None);
    /// assert_eq!(open.eval(&[4.0, 0.0]), Some(1.0));
    /// assert_eq!(near.with_distance_bound(f64::NEG_INFINITY).eval(&[0.0, 0.0]), None);
    /// # Ok(())
    /// # }
    /// ```
    #[must_use]
    pub fn with_distance_bound(mut self, bound: f64) -> Self {
        // Every distance is under positive infinity: that is no bound, and
        // is held as none, so that two interpolants that give the same
        // values are the same.
        let open = bound.is_infinite() && bound.is_sign_positive();
        self.bound = (!open).then_some(bound);
        self
    }

    /// Returns the bound on the distance, when there is one. See
    /// [`with_distance_bound`](NearestNDInterpolator::with_distance_bound) for
    /// an example.
    pub fn distance_bound(&self) -> Option<f64> {
        self.bound
    }

    /// Returns the number of points given. See [`NearestNDInterpolator`] for
    /// an example.
    #[expect(
        clippy::len_without_is_empty,
        reason = "a constructor refuses an empty set of points"
    )]
    pub fn len(&self) -> usize {
        self.values.len()
    }

    /// Returns the values, in the order of the points given.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{InterpolateError, NearestNDInterpolator};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let nearest = NearestNDInterpolator::new(&[[0.0], [1.0]], &[5.0, 6.0])?;
    /// assert_eq!(nearest.values(), &[5.0, 6.0]);
    /// # Ok(())
    /// # }
    /// ```
    pub fn values(&self) -> &[V] {
        &self.values
    }
}

/// The points that are the first with their coordinates, by index, in
/// increasing order: every set of coordinates once. `points` has no more than
/// `u32::MAX` entries.
///
/// Coordinates are compared as numbers: `0.0` and `-0.0` are the same point.
///
/// A set of the coordinates seen so far, with the hasher of the standard
/// library. Measured against sorting the points (ns a point, random points in
/// 2D): 33 for 32 at 100 points, 25 for 43 at 1000, 45 for 97 at 100 000, 52
/// for 231 at a million. Sorting wins only on points that come sorted already
/// (a lattice row by row: 12 for 27). A cheaper hasher (multiply and rotate)
/// was twice as fast on random points and five times slower on a lattice,
/// whose coordinates differ in few bits.
fn first_of_each<const D: usize>(points: &[[f64; D]]) -> Vec<u32> {
    let mut seen = HashSet::with_capacity(points.len());
    let mut first = Vec::with_capacity(points.len());
    for (i, point) in points.iter().enumerate() {
        if seen.insert(super::key(point)) {
            first.push(i as u32);
        }
    }
    first
}

/// The error of more points than the tree numbers: it numbers them with `u32`.
fn too_many() -> InterpolateError {
    InterpolateError::InvalidParameter {
        name: "points",
        reason: ParamReason::TooLarge {
            max: u32::MAX as usize,
        },
    }
}

/// The tree of the points kept, or none without dimensions.
fn tree_of<const D: usize>(
    kept: &[[f64; D]],
) -> Result<Option<Arc<ImmutableKdTree<f64, D>>>, InterpolateError> {
    if D == 0 {
        return Ok(None);
    }
    let tree = ImmutableKdTree::new_from_slice(kept).map_err(|_| too_many())?;
    Ok(Some(Arc::new(tree)))
}

impl<const D: usize, V: Value<f64>> Interpolator for NearestNDInterpolator<D, V> {
    type Point = [f64; D];
    type Output = V;

    fn eval(&self, point: &[f64; D]) -> Option<V> {
        if !point.iter().all(|c| c.is_finite()) {
            return None;
        }
        let query = self.scaling.map_or(*point, |scaling| scaling.apply(point));
        // The place in the tree of the nearest point, and the square of its
        // distance.
        let (item, squared) = match &self.tree {
            Some(tree) => {
                let nearest = tree
                    .query(&query)
                    .nearest_one::<SquaredEuclidean<f64>>()
                    .execute();
                (nearest.item as usize, nearest.distance)
            }
            None => (0, 0.0),
        };
        // Not finite when the query is too far from every point for its
        // distance to be measured.
        if !squared.is_finite() {
            return None;
        }
        if let Some(bound) = self.bound {
            // SciPy compares the distance itself, and takes what is nearer.
            let nearer = squared.sqrt() < bound;
            if !nearer {
                return None;
            }
        }
        Some(self.values[self.first[item] as usize])
    }
}

/// The stored form of the `serde` feature.
#[cfg(feature = "serde")]
mod stored {
    use super::{too_many, tree_of, NearestNDInterpolator, Scaling};
    use crate::{validate, InterpolateError, Value};
    use serde::{de, ser, Deserialize, Deserializer, Serialize, Serializer};
    use std::borrow::Cow;

    /// The version of the stored form.
    const VERSION: u32 = 1;

    /// What is stored: the points of the tree (each set of coordinates once,
    /// after the change of coordinates of a rescaled interpolant), which of
    /// the values each of them has, the values, the change of coordinates,
    /// and the bound. The tree is not: it is built again from its points when
    /// the form is read, so that the form does not depend on how a version of
    /// the tree's crate lays it out.
    #[derive(Serialize, Deserialize)]
    #[serde(bound = "")]
    struct Form<'a, const D: usize, V: Value<f64>> {
        version: u32,
        #[serde(with = "crate::stored")]
        points: Cow<'a, [[f64; D]]>,
        first: Cow<'a, [u32]>,
        #[serde(with = "crate::stored")]
        values: Cow<'a, [V]>,
        scaling: Option<Scaling<D>>,
        bound: Option<f64>,
    }

    impl<const D: usize, V: Value<f64>> Serialize for NearestNDInterpolator<D, V> {
        /// # Errors
        ///
        /// Besides the errors of the serializer: a distance bound that is NaN
        /// or negative infinity is not stored. With either, no query has a
        /// value; a text format has neither and would read them back as no
        /// bound, with which every query has one.
        fn serialize<S: Serializer>(&self, to: S) -> Result<S::Ok, S::Error> {
            // (Positive infinity is never held: `with_distance_bound`.)
            if self.bound.is_some_and(|bound| !bound.is_finite()) {
                return Err(ser::Error::custom(
                    "a distance bound that is not finite cannot be stored",
                ));
            }
            // The tree gives its points back with the place each was given at.
            let mut points = vec![[0.0; D]; self.first.len()];
            if let Some(tree) = &self.tree {
                for (place, point) in tree.iter() {
                    points[place as usize] = point;
                }
            }
            Form {
                version: VERSION,
                points: Cow::Owned(points),
                first: Cow::Borrowed(&self.first),
                values: Cow::Borrowed(&self.values),
                scaling: self.scaling,
                bound: self.bound,
            }
            .serialize(to)
        }
    }

    /// The points and the values pass the checks of
    /// [`NearestNDInterpolator::new`]. The points are those of the tree, so
    /// they are also checked to be distinct, and the places of their values to
    /// increase and to be places of the values.
    impl<'de, const D: usize, V: Value<f64>> Deserialize<'de> for NearestNDInterpolator<D, V> {
        fn deserialize<D2: Deserializer<'de>>(from: D2) -> Result<Self, D2::Error> {
            let form = Form::<D, V>::deserialize(from)?;
            crate::stored::version(form.version, VERSION)?;
            let (points, values) = (form.points.into_owned(), form.values.into_owned());
            let first = form.first.into_owned();
            let checked = || {
                validate::length(&first, points.len())?;
                validate::enough(&values, 1)?;
                validate::enough(&points, 1)?;
                if u32::try_from(values.len()).is_err() {
                    return Err(too_many());
                }
                let finite = |point: &[f64; D]| point.iter().all(|c| c.is_finite());
                if let Some(index) = validate::first_failing(&points, finite) {
                    return Err(InterpolateError::NonFinite { index });
                }
                if let Some(index) = validate::first_failing(&values, V::all_finite) {
                    return Err(InterpolateError::NonFinite { index });
                }
                if let Some((first, second)) = crate::scattered::first_repeat(&points) {
                    return Err(InterpolateError::DuplicatePoints { first, second });
                }
                Ok(())
            };
            checked().map_err(de::Error::custom)?;
            let in_order = first.windows(2).all(|pair| pair[0] < pair[1]);
            let in_range = first
                .last()
                .is_some_and(|last| (*last as usize) < values.len());
            if !(in_order && in_range) {
                return Err(de::Error::custom(
                    "the places of the values of the points do not increase within the values",
                ));
            }
            if let Some(scaling) = &form.scaling {
                let usable = scaling.offset.iter().all(|c| c.is_finite())
                    && scaling.scale.iter().all(|c| c.is_finite() && *c > 0.0);
                if !usable {
                    return Err(de::Error::custom(
                        "the change of coordinates is not finite, or its scale is not above zero",
                    ));
                }
            }
            // As none is stored: see `serialize`.
            if form.bound.is_some_and(|bound| !bound.is_finite()) {
                return Err(de::Error::custom("the distance bound is not finite"));
            }
            Ok(Self {
                tree: tree_of(&points).map_err(de::Error::custom)?,
                first,
                values,
                scaling: form.scaling,
                bound: form.bound,
            })
        }
    }
}
