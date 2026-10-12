//! The value of the nearest data point, in any number of dimensions.

use crate::error::ParamReason;
use crate::{validate, InterpolateError, Interpolator, Value};
use kiddo::{ImmutableKdTree, SquaredEuclidean};
use std::collections::HashSet;
use std::sync::Arc;

/// The interpolant that gives a query point the value of the data point
/// nearest to it, for points anywhere in `D` dimensions.
///
/// SciPy's `NearestNDInterpolator`. The points need no order and no structure.
/// The distance is the Euclidean one. A value is a scalar or a fixed-length
/// array, and the result is always one of the values given: nothing is
/// averaged. Coordinates are `f64`.
///
/// The points go into a k-d tree ([`kiddo`](https://docs.rs/kiddo)), built
/// once. It suits few dimensions, up to about ten: beyond that a query looks
/// at most of the points.
///
/// Points with the same coordinates are taken as one, with the value of the
/// one given first. Of several different points at the same distance from a
/// query, one is taken; which one is fixed for a given set of points, and
/// follows no rule (SciPy's choice follows none either, and may differ).
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

/// The change of coordinates of a rescaled interpolant: `(x - offset) / scale`
/// on each axis.
#[derive(Debug, Clone, Copy)]
struct Scaling<const D: usize> {
    offset: [f64; D],
    scale: [f64; D],
}

impl<const D: usize> Scaling<D> {
    /// The mean of the points on each axis, and their range. An axis on which
    /// all points agree keeps its scale.
    ///
    /// # Errors
    ///
    /// Returns [`InterpolateError::NonFinite`] when the sum of the coordinates
    /// of an axis overflows (with the point at which it does), or their range
    /// (with the point that has the largest coordinate). With a mean or a
    /// range that is not finite every point would land on one place, or on
    /// none.
    fn of(points: &[[f64; D]]) -> Result<Self, InterpolateError> {
        let count = points.len() as f64;
        let mut offset = [0.0; D];
        let mut scale = [1.0; D];
        for d in 0..D {
            let (mut sum, mut least, mut most) = (0.0, f64::INFINITY, f64::NEG_INFINITY);
            let mut largest = 0;
            for (index, point) in points.iter().enumerate() {
                sum += point[d];
                if !sum.is_finite() {
                    return Err(InterpolateError::NonFinite { index });
                }
                least = least.min(point[d]);
                if point[d] > most {
                    (most, largest) = (point[d], index);
                }
            }
            offset[d] = sum / count;
            let range = most - least;
            if !range.is_finite() {
                return Err(InterpolateError::NonFinite { index: largest });
            }
            if range > 0.0 {
                scale[d] = range;
            }
        }
        Ok(Self { offset, scale })
    }

    fn apply(&self, point: &[f64; D]) -> [f64; D] {
        core::array::from_fn(|d| (point[d] - self.offset[d]) / self.scale[d])
    }
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
        Self::build(points, values, false)
    }

    /// Builds the interpolant with every coordinate brought to the same range
    /// before a distance is measured: SciPy's `rescale=True`.
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
        Self::build(points, values, true)
    }

    fn build(points: &[[f64; D]], values: &[V], rescale: bool) -> Result<Self, InterpolateError> {
        validate::length(values, points.len())?;
        validate::enough(points, 1)?;
        // The tree numbers its points with `u32`.
        let too_many = || InterpolateError::InvalidParameter {
            name: "points",
            reason: ParamReason::TooLarge {
                max: u32::MAX as usize,
            },
        };
        if u32::try_from(points.len()).is_err() {
            return Err(too_many());
        }
        let not_finite = |point: &[f64; D]| !point.iter().all(|c| c.is_finite());
        if let Some(index) = points.iter().position(not_finite) {
            return Err(InterpolateError::NonFinite { index });
        }
        if let Some(index) = values.iter().position(|v| !v.all_finite()) {
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
        let tree = if D == 0 {
            None
        } else {
            let tree = ImmutableKdTree::new_from_slice(&kept).map_err(|_| too_many())?;
            Some(Arc::new(tree))
        };
        Ok(Self {
            tree,
            first,
            values: values.to_vec(),
            scaling,
            bound: None,
        })
    }

    /// Returns the same interpolant, which gives a value only to queries whose
    /// nearest point is nearer than `bound`: SciPy's `distance_upper_bound`.
    ///
    /// A query at the distance `bound` exactly, or farther, has no value. For
    /// a [`rescaled`](NearestNDInterpolator::rescaled) interpolant the distance
    /// is the one measured after the change of coordinates. A `bound` that is
    /// NaN leaves no query with a value; an infinite one is no bound.
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
    /// # Ok(())
    /// # }
    /// ```
    #[must_use]
    pub fn with_distance_bound(mut self, bound: f64) -> Self {
        self.bound = Some(bound);
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
        // `-0.0 + 0.0` is `0.0`: the two zeros have the same bits after it.
        if seen.insert(point.map(|c| (c + 0.0).to_bits())) {
            first.push(i as u32);
        }
    }
    first
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
