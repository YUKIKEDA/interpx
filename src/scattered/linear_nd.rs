//! Linear interpolation on the Delaunay triangulation of the points, in any
//! number of dimensions.

use super::mesh::Mesh;
use crate::{InterpolateError, Interpolator, Value};
use std::borrow::Cow;

/// The interpolant that is linear on each simplex of the Delaunay
/// triangulation of the points, for points anywhere in `D` dimensions.
///
/// The points need no order and no structure.
/// At a query point the value is the weighted sum of the values at the
/// vertices of the simplex that holds it, the weights being the point's
/// barycentric coordinates in that simplex. The interpolant is continuous,
/// reproduces any affine function exactly, and gives the data values at the
/// data points. A value is a scalar or a fixed-length array. Coordinates are
/// `f64`.
///
/// A query outside the convex hull of the points has no value, nor has a query
/// with a coordinate that is NaN or infinite. A query on a face of the hull,
/// or outside it by no more than rounding (barycentric coordinates down to
/// `-100 * f64::EPSILON`), has a value.
///
/// It is meant for few dimensions. The triangulation of `n` points can have of
/// the order of `n` to the power `D / 2` (rounded up) simplices, and each keeps
/// `2 * D * D + D` numbers, so that in high dimensions the memory
/// can be far more than that of the points.
///
/// Points with the same coordinates are taken as one, with the value of the
/// one given first. [`rescaled`](LinearNDInterpolator::rescaled) triangulates
/// the points after each coordinate is brought to the same range, for data
/// whose coordinates have different units.
///
/// The simplex of a query is found by a walk across the faces toward it, from
/// a simplex that the construction chose for the part of the box of the
/// points the query is in (the box is cut into at most twice as many cells
/// as there are points: `p` along each axis, the largest `p` with `p^D` no
/// more than twice the number of points). Where the walk starts depends on
/// the query alone, so a
/// batch, and a parallel one, gives what `eval` gives, bit for bit. A walk
/// that does not arrive within as many steps as there are simplices, or that
/// meets a degenerate simplex, gives way to a look at every simplex.
///
/// # Ties and small dimensions
///
/// - Where four or more points lie on one circle (five on one sphere in three
///   dimensions, and so on), more than one Delaunay triangulation exists.
///   This crate's triangulation divides the region they bound by a rule of
///   the points' indices, so that the same points in the same order always
///   give the same interpolant. Points on a grid are such points.
/// - One dimension (`D = 1`) is linear interpolation between neighboring
///   points. Without dimensions (`D = 0`) every point is the one point there
///   is, and its value is the first one given.
///
/// # Examples
///
/// ```
/// use interpx::{InterpolateError, Interpolator, LinearNDInterpolator};
///
/// # fn main() -> Result<(), InterpolateError> {
/// // The corners of the unit square, and f(x, y) = 1 + 2x + 3y at each.
/// let points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]];
/// let values = [1.0, 3.0, 4.0, 6.0];
/// let linear = LinearNDInterpolator::new(&points, &values)?;
/// // An affine function comes back exactly (to rounding).
/// let inside = linear.eval(&[0.25, 0.5]).unwrap();
/// assert!((inside - 3.0).abs() < 1e-12);
/// // Outside the square there is no value.
/// assert_eq!(linear.eval(&[1.5, 0.5]), None);
/// assert_eq!(linear.eval(&[f64::NAN, 0.5]), None);
/// assert_eq!(linear.len(), 4);
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone)]
pub struct LinearNDInterpolator<const D: usize, V: Value<f64> = f64> {
    /// The triangulation of the points, and what a query is located by.
    mesh: Mesh<D>,
    /// The values, as given.
    values: Vec<V>,
}

impl<const D: usize, V: Value<f64>> LinearNDInterpolator<D, V> {
    /// Builds the interpolant from copies of the points and of their values.
    ///
    /// # Errors
    ///
    /// The first of these that applies:
    ///
    /// - [`InterpolateError::LengthMismatch`] when `points` and `values` differ
    ///   in length.
    /// - [`InterpolateError::InsufficientPoints`] when there are fewer than
    ///   `D + 1` points.
    /// - [`InterpolateError::InvalidParameter`] for `points` when there are
    ///   more than `u32::MAX` of them.
    /// - [`InterpolateError::NonFinite`] when a coordinate, or else a value, is
    ///   NaN or infinite. Its `index` is that of the point.
    /// - [`InterpolateError::InsufficientPoints`] when fewer than `D + 1`
    ///   points are distinct; `actual` is the number of distinct points.
    /// - [`InterpolateError::DegenerateGeometry`] with
    ///   [`DegeneracyKind::LowerDimensional`](crate::DegeneracyKind::LowerDimensional) when the points lie in an affine
    ///   subspace of fewer than `D` dimensions (on one line in the plane), or
    ///   with [`DegeneracyKind::NumericalLimit`](crate::DegeneracyKind::NumericalLimit) when the working space of the
    ///   exact arithmetic that decides the triangulation cannot be allocated.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{InterpolateError, Interpolator, LinearNDInterpolator};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// // Three dimensions, two values at every point.
    /// let points = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    /// let values = [[0.0, 1.0], [1.0, 1.0], [2.0, 1.0], [3.0, 1.0]];
    /// let linear = LinearNDInterpolator::new(&points, &values)?;
    /// let [a, b] = linear.eval(&[0.25, 0.25, 0.25]).unwrap();
    /// assert!((a - 1.5).abs() < 1e-12 && (b - 1.0).abs() < 1e-12);
    ///
    /// // The same point twice: the value given first.
    /// let twice = LinearNDInterpolator::new(&[[0.0], [1.0], [0.0]], &[7.0, 9.0, 8.0])?;
    /// assert_eq!(twice.eval(&[0.0]), Some(7.0));
    ///
    /// let err = LinearNDInterpolator::<2>::new(&[[0.0, 0.0], [1.0, 0.0]], &[1.0, 2.0]);
    /// assert!(matches!(err, Err(InterpolateError::InsufficientPoints { required: 3, actual: 2 })));
    /// # Ok(())
    /// # }
    /// ```
    pub fn new(points: &[[f64; D]], values: &[V]) -> Result<Self, InterpolateError> {
        Self::build(Cow::Borrowed(points), Cow::Borrowed(values), false)
    }

    /// Builds the interpolant from points and values that are handed over,
    /// which are kept and not copied.
    ///
    /// The checks and the interpolant are those of
    /// [`new`](LinearNDInterpolator::new).
    ///
    /// # Errors
    ///
    /// As [`new`](LinearNDInterpolator::new).
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{InterpolateError, Interpolator, LinearNDInterpolator};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let points = vec![[0.0, 0.0], [2.0, 0.0], [0.0, 2.0]];
    /// let linear = LinearNDInterpolator::from_owned(points, vec![0.0, 2.0, 4.0])?;
    /// assert_eq!(linear.eval(&[1.0, 0.0]), Some(1.0));
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_owned(points: Vec<[f64; D]>, values: Vec<V>) -> Result<Self, InterpolateError> {
        Self::build(Cow::Owned(points), Cow::Owned(values), false)
    }

    /// Builds the interpolant with every coordinate brought to the same range
    /// before the points are triangulated.
    ///
    /// On each axis the mean of the points is taken off and the result is
    /// divided by their range (the largest coordinate less the smallest; an
    /// axis on which all points agree is left as it is). Queries go through
    /// the same change. It is the change of
    /// [`NearestNDInterpolator::rescaled`](crate::NearestNDInterpolator::rescaled).
    /// The triangulation, and so the interpolant, can differ from that of the
    /// points as they are: the Delaunay triangulation depends on the units.
    ///
    /// # Errors
    ///
    /// Those of [`new`](LinearNDInterpolator::new) up to the checks of the
    /// coordinates and values, then [`InterpolateError::NonFinite`] when the
    /// mean or the range of the coordinates of an axis overflows (its `index`
    /// is that of a point at which it does: the one that takes the sum of the
    /// axis past the largest number, or the one with the largest coordinate),
    /// then those of the triangulation.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{InterpolateError, Interpolator, LinearNDInterpolator};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// // Positions in kilometres, heights in metres.
    /// let points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1000.0], [1.0, 1000.0]];
    /// let values = [1.0, 2.0, 3.0, 4.0];
    /// let rescaled = LinearNDInterpolator::rescaled(&points, &values)?;
    /// let value = rescaled.eval(&[0.5, 500.0]).unwrap();
    /// assert!((value - 2.5).abs() < 1e-12);
    /// # Ok(())
    /// # }
    /// ```
    pub fn rescaled(points: &[[f64; D]], values: &[V]) -> Result<Self, InterpolateError> {
        Self::build(Cow::Borrowed(points), Cow::Borrowed(values), true)
    }

    /// Builds as [`rescaled`](LinearNDInterpolator::rescaled) does, from points
    /// and values that are handed over, as
    /// [`from_owned`](LinearNDInterpolator::from_owned) is to
    /// [`new`](LinearNDInterpolator::new). The points are changed in place.
    ///
    /// # Errors
    ///
    /// As [`rescaled`](LinearNDInterpolator::rescaled).
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{InterpolateError, Interpolator, LinearNDInterpolator};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let points = vec![[0.0, 0.0], [1000.0, 0.0], [0.0, 1.0]];
    /// let scaled = LinearNDInterpolator::rescaled_from_owned(points, vec![1.0, 2.0, 3.0])?;
    /// assert_eq!(scaled.eval(&[0.0, 0.0]), Some(1.0));
    /// # Ok(())
    /// # }
    /// ```
    pub fn rescaled_from_owned(
        points: Vec<[f64; D]>,
        values: Vec<V>,
    ) -> Result<Self, InterpolateError> {
        Self::build(Cow::Owned(points), Cow::Owned(values), true)
    }

    fn build(
        points: Cow<'_, [[f64; D]]>,
        values: Cow<'_, [V]>,
        rescale: bool,
    ) -> Result<Self, InterpolateError> {
        Mesh::check(&points, &values)?;
        let mesh = Mesh::build(points, rescale)?;
        Ok(Self {
            mesh,
            values: values.into_owned(),
        })
    }

    /// Returns the number of points given. See [`LinearNDInterpolator`] for
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
    /// use interpx::{InterpolateError, LinearNDInterpolator};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let linear = LinearNDInterpolator::new(&[[0.0], [1.0]], &[5.0, 6.0])?;
    /// assert_eq!(linear.values(), &[5.0, 6.0]);
    /// # Ok(())
    /// # }
    /// ```
    pub fn values(&self) -> &[V] {
        &self.values
    }

    /// Returns the points, as triangulated: as given, or after the change of
    /// coordinates of a [`rescaled`](LinearNDInterpolator::rescaled)
    /// interpolant.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{InterpolateError, LinearNDInterpolator};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let points = [[0.0, 0.0], [2.0, 0.0], [0.0, 4.0]];
    /// let plain = LinearNDInterpolator::new(&points, &[0.0; 3])?;
    /// assert_eq!(plain.points(), &points);
    /// // Rescaled: the mean taken off, divided by the range.
    /// let rescaled = LinearNDInterpolator::rescaled(&points, &[0.0; 3])?;
    /// assert_eq!(rescaled.points()[1], [(2.0 - 2.0 / 3.0) / 2.0, (0.0 - 4.0 / 3.0) / 4.0]);
    /// # Ok(())
    /// # }
    /// ```
    pub fn points(&self) -> &[[f64; D]] {
        &self.mesh.points
    }

    /// Returns the number of simplices of the triangulation (none without
    /// dimensions).
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{InterpolateError, LinearNDInterpolator};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// // A square and its center: four triangles.
    /// let points = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.5, 0.5]];
    /// let linear = LinearNDInterpolator::new(&points, &[0.0; 5])?;
    /// assert_eq!(linear.simplex_count(), 4);
    /// # Ok(())
    /// # }
    /// ```
    pub fn simplex_count(&self) -> usize {
        self.mesh.simplex_count()
    }

    /// Returns the same interpolant, with the tolerance by which a query is
    /// taken to be in a simplex multiplied by `tolerance`.
    ///
    /// A barycentric coordinate down to `-100 * f64::EPSILON * tolerance`
    /// (and up to one more than that) counts as inside. Raise it when points
    /// on the edges of very oblique simplices find none. The wider leeway of
    /// the search toward a degenerate simplex (`sqrt(f64::EPSILON)`) does not
    /// change. The default is 1.
    ///
    /// # Errors
    ///
    /// [`InterpolateError::InvalidParameter`] for `simplex_tolerance` when
    /// `tolerance` is NaN or infinite ([`ParamReason::NonFinite`](crate::ParamReason::NonFinite)), or zero
    /// or below ([`ParamReason::NotPositive`](crate::ParamReason::NotPositive)).
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{InterpolateError, Interpolator, LinearNDInterpolator};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let linear = LinearNDInterpolator::new(&[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], &[0.0; 3])?;
    /// assert_eq!(linear.simplex_tolerance(), 1.0);
    /// // Outside the edge x = 0 by 1e-13: farther than the default tolerance.
    /// let near = [-1e-13, 0.5];
    /// assert_eq!(linear.eval(&near), None);
    /// let wider = linear.with_simplex_tolerance(10.0)?;
    /// assert_eq!(wider.eval(&near), Some(0.0));
    /// assert!(wider.clone().with_simplex_tolerance(0.0).is_err());
    /// # Ok(())
    /// # }
    /// ```
    pub fn with_simplex_tolerance(mut self, tolerance: f64) -> Result<Self, InterpolateError> {
        self.mesh.set_tolerance(tolerance)?;
        Ok(self)
    }

    /// Returns the simplex tolerance. See
    /// [`with_simplex_tolerance`](LinearNDInterpolator::with_simplex_tolerance)
    /// for an example.
    pub fn simplex_tolerance(&self) -> f64 {
        self.mesh.tolerance
    }
}

impl<const D: usize, V: Value<f64>> Interpolator for LinearNDInterpolator<D, V> {
    type Point = [f64; D];
    type Output = V;

    fn eval(&self, point: &[f64; D]) -> Option<V> {
        if !point.iter().all(|c| c.is_finite()) {
            return None;
        }
        if D == 0 {
            // Without dimensions every point is the first one.
            return Some(self.values[0]);
        }
        let (simplex, c, last) = self.mesh.find(point)?;
        // SciPy's order: from zero, the vertices in turn.
        let mut value = V::zero();
        for (k, &v) in self.mesh.vertices(simplex).iter().enumerate() {
            let weight = if k < D { c[k] } else { last };
            value = value.axpy(weight, self.values[v as usize]);
        }
        Some(value)
    }
}

/// The stored form of the `serde` feature.
#[cfg(feature = "serde")]
mod stored {
    use super::{LinearNDInterpolator, Mesh};
    use crate::scattered::Scaling;
    use crate::Value;
    use serde::{de, Deserialize, Deserializer, Serialize, Serializer};
    use std::borrow::Cow;

    /// The version of the stored form.
    const VERSION: u32 = 1;

    /// What is stored: the points (after the change of coordinates of a
    /// rescaled interpolant), the values, the simplices and their neighbors,
    /// and the change of coordinates. The triangulation is not made again
    /// when the form is read, so that the form does not depend on how a
    /// version of the geometry crate orders it; the transforms, the box, and
    /// the starts of its cells are made again from what is stored.
    #[derive(Serialize, Deserialize)]
    #[serde(bound = "")]
    struct Form<'a, const D: usize, V: Value<f64>> {
        version: u32,
        #[serde(with = "crate::stored")]
        points: Cow<'a, [[f64; D]]>,
        #[serde(with = "crate::stored")]
        values: Cow<'a, [V]>,
        simplices: Cow<'a, [u32]>,
        neighbors: Cow<'a, [u32]>,
        scaling: Option<Scaling<D>>,
        simplex_tolerance: f64,
    }

    impl<const D: usize, V: Value<f64>> Serialize for LinearNDInterpolator<D, V> {
        fn serialize<S: Serializer>(&self, to: S) -> Result<S::Ok, S::Error> {
            Form {
                version: VERSION,
                points: Cow::Borrowed(&self.mesh.points),
                values: Cow::Borrowed(&self.values),
                simplices: Cow::Borrowed(&self.mesh.simplices),
                neighbors: Cow::Borrowed(&self.mesh.neighbors),
                scaling: self.mesh.scaling,
                simplex_tolerance: self.mesh.tolerance,
            }
            .serialize(to)
        }
    }

    /// The points and the values pass the checks of
    /// [`LinearNDInterpolator::new`] up to the triangulation, and the
    /// triangulation those of the shared mesh (`Mesh::from_stored`). That it
    /// is a Delaunay triangulation is not checked again: a walk over a form
    /// that passes these and is still wrong ends at its limit, as many steps
    /// as there are simplices, and the search over every simplex follows.
    /// Building the starts of the cells takes such walks too, at most one a
    /// cell, and there are no more cells than twice the points.
    impl<'de, const D: usize, V: Value<f64>> Deserialize<'de> for LinearNDInterpolator<D, V> {
        fn deserialize<D2: Deserializer<'de>>(from: D2) -> Result<Self, D2::Error> {
            let form = Form::<D, V>::deserialize(from)?;
            crate::stored::version(form.version, VERSION)?;
            let points = form.points.into_owned();
            let values = form.values.into_owned();
            Mesh::check(&points, &values).map_err(de::Error::custom)?;
            let mesh = Mesh::from_stored(
                points,
                form.simplices.into_owned(),
                form.neighbors.into_owned(),
                form.scaling,
                form.simplex_tolerance,
            )
            .map_err(de::Error::custom)?;
            Ok(Self { mesh, values })
        }
    }
}
