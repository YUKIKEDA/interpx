//! Linear interpolation on the Delaunay triangulation of the points, in any
//! number of dimensions.

use super::Scaling;
use crate::error::{DegeneracyKind, ParamReason};
use crate::{validate, InterpolateError, Interpolator, Value};
use std::borrow::Cow;

/// The neighbor across a face of the hull: no simplex. `convx` writes it so and
/// promises it in its documentation; it has no constant of its own for it.
const NO_NEIGHBOR: u32 = u32::MAX;

/// How far below zero (or above one) a barycentric coordinate may be for the
/// point to be inside the simplex: SciPy's `eps`.
const EPS: f64 = 100.0 * f64::EPSILON;

/// The wider leeway toward a degenerate simplex, in the search over every
/// simplex: SciPy's `eps_broad`, the square root of `f64::EPSILON` (2^-26).
const EPS_BROAD: f64 = 1.490_116_119_384_765_6e-8;

/// A simplex whose matrix has a reciprocal condition number below this is
/// degenerate: it has no transform, and no point is found inside it by its own
/// coordinates. SciPy's limit, `1000 * eps`.
const RCOND_LIMIT: f64 = 1000.0 * f64::EPSILON;

/// The interpolant that is linear on each simplex of the Delaunay
/// triangulation of the points, for points anywhere in `D` dimensions.
///
/// SciPy's `LinearNDInterpolator`. The points need no order and no structure.
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
/// Points with the same coordinates are taken as one, with the value of the
/// one given first. [`rescaled`](LinearNDInterpolator::rescaled) triangulates
/// the points after each coordinate is brought to the same range, for data
/// whose coordinates have different units.
///
/// The simplex of a query is found by a walk across the faces toward it, from
/// a simplex that the construction chose for the part of the box of the
/// points the query is in (the box is cut into about as many cells as there
/// are simplices). Where the walk starts depends on the query alone, so a
/// batch, and a parallel one, gives what `eval` gives, bit for bit. A walk
/// that does not arrive within as many steps as there are simplices, or that
/// meets a degenerate simplex, gives way to a look at every simplex.
///
/// # Differences from SciPy
///
/// - Where four or more points lie on one circle (five on one sphere in three
///   dimensions, and so on), more than one Delaunay triangulation exists.
///   This crate's triangulation divides the region they bound by a rule of
///   the points' indices; SciPy's (Qhull) divides it another way. Inside
///   such a region the two interpolants differ; elsewhere they agree to
///   rounding. Points on a grid are such points.
/// - Of points with the same coordinates, the value of the one given first
///   is taken. SciPy (Qhull) keeps one of them by its own choice, which is
///   not always the first.
/// - One dimension (`D = 1`) is linear interpolation between neighboring
///   points; SciPy refuses one-dimensional points. Without dimensions
///   (`D = 0`) every point is the one point there is, and its value is the
///   first one given.
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
    /// The points, as given, moved and scaled when the interpolant is
    /// rescaled. Duplicates stay, and are no vertex.
    points: Vec<[f64; D]>,
    /// The values, as given.
    values: Vec<V>,
    /// The `D + 1` vertices of each simplex: places in `points`.
    simplices: Vec<u32>,
    /// The `D + 1` simplices across the faces of each simplex: the one across
    /// the face opposite each vertex, or [`NO_NEIGHBOR`] on the hull.
    neighbors: Vec<u32>,
    /// For each simplex, what its barycentric coordinates are solved with
    /// (see `push_transform`): `L^-1 P` and `U` of the matrix whose columns
    /// are its first `D` vertices less the last one (`D * D` numbers each),
    /// then the last vertex (`D`). All NaN for a degenerate simplex.
    transforms: Vec<f64>,
    /// The smallest and the largest coordinate of the points on each axis.
    lower: [f64; D],
    upper: [f64; D],
    /// The cells of the box, `per_axis` along each axis, the last axis
    /// fastest: for each, the simplex a walk to a query in it starts from.
    /// None without dimensions.
    cells: Vec<u32>,
    per_axis: usize,
    /// The cells per unit of each coordinate.
    cell_scale: [f64; D],
    /// What is taken off each coordinate and what it is divided by, before
    /// the points are triangulated and a query is located.
    scaling: Option<Scaling<D>>,
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
    ///   [`DegeneracyKind::LowerDimensional`] when the points lie in an affine
    ///   subspace of fewer than `D` dimensions (on one line in the plane), or
    ///   with [`DegeneracyKind::NumericalLimit`] when the working space of the
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
    /// before the points are triangulated: SciPy's `rescale=True`.
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
        validate::length(&values, points.len())?;
        validate::enough(&points, D + 1)?;
        if u32::try_from(points.len()).is_err() {
            return Err(too_many());
        }
        let finite = |point: &[f64; D]| point.iter().all(|c| c.is_finite());
        if let Some(index) = validate::first_failing(&points, finite) {
            return Err(InterpolateError::NonFinite { index });
        }
        if let Some(index) = validate::first_failing(&values, V::all_finite) {
            return Err(InterpolateError::NonFinite { index });
        }

        let scaling = if rescale {
            Some(Scaling::of(&points)?)
        } else {
            None
        };
        let mut points = points.into_owned();
        if let Some(scaling) = &scaling {
            for point in &mut points {
                *point = scaling.apply(point);
            }
        }
        let (simplices, neighbors) = triangulate(&points)?;
        Self::assemble(points, values.into_owned(), simplices, neighbors, scaling)
    }

    /// The interpolant of a triangulation: the transforms of the simplices,
    /// the box of the points, and a start for the walk in each cell of it.
    fn assemble(
        points: Vec<[f64; D]>,
        values: Vec<V>,
        simplices: Vec<u32>,
        neighbors: Vec<u32>,
        scaling: Option<Scaling<D>>,
    ) -> Result<Self, InterpolateError> {
        let mut transforms = Vec::with_capacity(simplices.len() / (D + 1) * (2 * D * D + D));
        for simplex in simplices.chunks_exact(D + 1) {
            push_transform(&points, simplex, &mut transforms);
        }
        let mut lower = [f64::INFINITY; D];
        let mut upper = [f64::NEG_INFINITY; D];
        for point in &points {
            for d in 0..D {
                lower[d] = lower[d].min(point[d]);
                upper[d] = upper[d].max(point[d]);
            }
        }
        let mut linear = Self {
            points,
            values,
            simplices,
            neighbors,
            transforms,
            lower,
            upper,
            cells: Vec::new(),
            per_axis: 0,
            cell_scale: [0.0; D],
            scaling,
        };
        linear.place_cells();
        Ok(linear)
    }

    /// Cuts the box into about as many cells as there are simplices, and
    /// walks to the middle of each from the start of the cell before it: the
    /// last simplex of that walk (the one that holds the middle, or one on
    /// the hull toward it) is the start of the cell.
    fn place_cells(&mut self) {
        let count = self.simplex_count();
        if D == 0 || count == 0 {
            return;
        }
        let per_axis = ((self.points.len() as f64).powf(1.0 / D as f64).ceil() as usize).max(1);
        // Every axis has a range above zero: the points span D dimensions.
        let width: [f64; D] =
            core::array::from_fn(|d| (self.upper[d] - self.lower[d]) / per_axis as f64);
        self.cell_scale = core::array::from_fn(|d| 1.0 / width[d]);
        self.per_axis = per_axis;
        let total = per_axis.pow(D as u32);
        self.cells = Vec::with_capacity(total);
        let mut start = 0;
        let mut at = [0_usize; D];
        for _ in 0..total {
            let middle: [f64; D] =
                core::array::from_fn(|d| self.lower[d] + (at[d] as f64 + 0.5) * width[d]);
            start = match self.walk(&middle, start) {
                Walked::Found(simplex, ..)
                | Walked::Outside(simplex, _)
                | Walked::Lost(simplex) => simplex,
            };
            self.cells.push(start);
            // The next cell: the last axis fastest.
            for d in (0..D).rev() {
                at[d] += 1;
                if at[d] < per_axis {
                    break;
                }
                at[d] = 0;
            }
        }
    }

    /// The simplex the walk to `x` starts from: that of the cell `x` is in
    /// (the nearest cell, for a point just outside the box).
    fn start(&self, x: &[f64; D]) -> u32 {
        let mut cell = 0;
        for (d, coordinate) in x.iter().enumerate() {
            let along = ((coordinate - self.lower[d]) * self.cell_scale[d]).floor();
            let place = along.clamp(0.0, (self.per_axis - 1) as f64) as usize;
            cell = cell * self.per_axis + place;
        }
        self.cells[cell]
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
    /// interpolant (SciPy's `points`, which are rescaled too).
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
        &self.points
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
        self.simplices.len() / (D + 1)
    }

    /// The transform of a simplex: `L^-1 P`, `U`, the last vertex.
    fn transform(&self, simplex: u32) -> &[f64] {
        let stride = 2 * D * D + D;
        let start = simplex as usize * stride;
        &self.transforms[start..start + stride]
    }

    /// The barycentric coordinates of `x` in `simplex`: those of the first
    /// `D` vertices, and that of the last. NaN for a degenerate simplex.
    fn barycentric(&self, simplex: u32, x: &[f64; D]) -> ([f64; D], f64) {
        barycentric(self.transform(simplex), x)
    }

    fn neighbor(&self, simplex: u32, face: usize) -> u32 {
        self.neighbors[simplex as usize * (D + 1) + face]
    }

    /// Whether `x` is outside the box of the points by more than [`EPS`] on
    /// some axis: then it is outside the hull. SciPy's
    /// `_is_point_fully_outside`.
    fn fully_outside(&self, x: &[f64; D]) -> bool {
        (0..D).any(|d| x[d] < self.lower[d] - EPS || x[d] > self.upper[d] + EPS)
    }

    /// The walk from `start` toward `x`: across the face opposite the first
    /// vertex whose coordinate is below `-EPS`, as many steps at most as there
    /// are simplices. SciPy's `_find_simplex_directed`.
    fn walk(&self, x: &[f64; D], start: u32) -> Walked<D> {
        let mut simplex = start;
        for _ in 0..self.simplex_count() {
            let (c, last) = self.barycentric(simplex, x);
            let mut inside = true;
            let mut across = None;
            for (k, ck) in c.iter().copied().chain([last]).enumerate() {
                if ck < -EPS {
                    across = Some((k, ck));
                    break;
                }
                // Above one, or NaN: a degenerate simplex is in the way.
                if ck.is_nan() || ck > 1.0 + EPS {
                    inside = false;
                }
            }
            match across {
                Some((face, ck)) => {
                    let next = self.neighbor(simplex, face);
                    if next == NO_NEIGHBOR {
                        return Walked::Outside(simplex, ck);
                    }
                    simplex = next;
                }
                None if inside => return Walked::Found(simplex, c, last),
                None => return Walked::Lost(simplex),
            }
        }
        Walked::Lost(simplex)
    }

    /// The simplex that holds `x`, with the coordinates of `x` in it: by the
    /// walk from `start`, or by the look at every simplex when the walk is
    /// lost.
    ///
    /// Where the walk would cross a face of the hull, `x` is outside it,
    /// unless its coordinate is within [`EPS_BROAD`] of zero: then the
    /// transform's rounding may have put a point of the hull outside (a vertex
    /// of the simplex was measured at `-2.8e-14`), and every simplex is looked
    /// at. SciPy's walk stops at a face of the hull without this look; on
    /// the set where this was found, SciPy (another triangulation, other
    /// transforms) gives every data point a value, from whatever point its
    /// walk starts.
    fn locate(&self, x: &[f64; D], start: u32) -> Option<(u32, [f64; D], f64)> {
        match self.walk(x, start) {
            Walked::Found(simplex, c, last) => Some((simplex, c, last)),
            Walked::Outside(_, ck) if ck < -EPS_BROAD => None,
            Walked::Outside(..) | Walked::Lost(_) => self.search(x),
        }
    }

    /// The simplex that holds `x`, looked for in every simplex. A degenerate
    /// simplex is looked through: `x` is taken in a neighbor of it when it is
    /// inside the neighbor with [`EPS_BROAD`] of leeway toward the degenerate
    /// one. SciPy's `_find_simplex_bruteforce`.
    fn search(&self, x: &[f64; D]) -> Option<(u32, [f64; D], f64)> {
        let within = |c: &[f64; D], last: f64, low: &dyn Fn(usize) -> f64| {
            c.iter()
                .copied()
                .chain([last])
                .enumerate()
                .all(|(k, ck)| -low(k) <= ck && ck <= 1.0 + EPS)
        };
        for simplex in 0..self.simplex_count() as u32 {
            if !self.transform(simplex)[0].is_nan() {
                let (c, last) = self.barycentric(simplex, x);
                if within(&c, last, &|_| EPS) {
                    return Some((simplex, c, last));
                }
                continue;
            }
            for face in 0..=D {
                let next = self.neighbor(simplex, face);
                if next == NO_NEIGHBOR || self.transform(next)[0].is_nan() {
                    continue;
                }
                let (c, last) = self.barycentric(next, x);
                let toward = |k: usize| {
                    if self.neighbor(next, k) == simplex {
                        EPS_BROAD
                    } else {
                        EPS
                    }
                };
                if within(&c, last, &toward) {
                    return Some((next, c, last));
                }
            }
        }
        None
    }
}

/// Where a walk ended: in the simplex that holds the point (with its
/// coordinates there), beyond a face of the hull of a simplex (with the
/// coordinate that sent it there), or neither (its steps ran out, or a
/// degenerate simplex was in the way) in a simplex.
enum Walked<const D: usize> {
    Found(u32, [f64; D], f64),
    Outside(u32, f64),
    Lost(u32),
}

/// The simplices of the Delaunay triangulation of the points, and their
/// neighbors, `D + 1` of each a simplex. Without dimensions there is none.
fn triangulate<const D: usize>(
    points: &[[f64; D]],
) -> Result<(Vec<u32>, Vec<u32>), InterpolateError> {
    if D == 0 {
        return Ok((Vec::new(), Vec::new()));
    }
    let delaunay = convx::DelaunayBuilder::new(D, points.as_flattened())
        .build()
        .map_err(|error| from_convx(error, points.len()))?;
    let simplices = delaunay.simplices();
    let mut vertices = Vec::with_capacity(simplices.len() * (D + 1));
    let mut neighbors = Vec::with_capacity(simplices.len() * (D + 1));
    for simplex in simplices.iter() {
        vertices.extend_from_slice(simplex.vertices());
        neighbors.extend_from_slice(simplex.neighbors());
    }
    Ok((vertices, neighbors))
}

/// The error of this crate for a failure of the triangulation. The checks of
/// the constructor come first, so of the failures of `convx` only the count
/// of distinct points, the dimension they span, and the exact arithmetic can
/// happen; the others are mapped all the same, so that none is lost.
fn from_convx(error: convx::Error, len: usize) -> InterpolateError {
    match error {
        convx::Error::InsufficientPoints { actual, required } => {
            InterpolateError::InsufficientPoints { required, actual }
        }
        convx::Error::DegenerateDimension { actual_dim, .. } => {
            InterpolateError::DegenerateGeometry(DegeneracyKind::LowerDimensional {
                affine_dim: actual_dim,
            })
        }
        convx::Error::NonFiniteCoordinate { index } => InterpolateError::NonFinite { index },
        convx::Error::LengthMismatch { .. } => InterpolateError::LengthMismatch {
            expected: len,
            actual: len,
        },
        convx::Error::TooManyPoints { .. } => too_many(),
        // `ExactEvaluationExhausted`; a dimension of zero, which is not given
        // to it; the failures of the hull's planes and of the Voronoi diagram,
        // which a triangulation does not return; and a failure a later version
        // adds.
        _ => InterpolateError::DegenerateGeometry(DegeneracyKind::NumericalLimit),
    }
}

/// The error of more points than the triangulation numbers: `u32`.
fn too_many() -> InterpolateError {
    InterpolateError::InvalidParameter {
        name: "points",
        reason: ParamReason::TooLarge {
            max: u32::MAX as usize,
        },
    }
}

/// Appends the transform of a simplex, from the LU decomposition with
/// partial pivoting `P A = L U` of the matrix `A` whose column `j` is vertex
/// `j` less the last vertex: `M = L^-1 P` (`D * D` numbers, row by row), then
/// `U` with the reciprocals of its diagonal in place of the diagonal (`D * D`,
/// row by row, what is below the diagonal unused), then the last vertex. The
/// barycentric coordinates of `x` solve `U c = M (x - last)`, and the last is
/// one less their sum.
///
/// A matrix that is singular, or whose reciprocal condition number in the
/// 1-norm is below [`RCOND_LIMIT`], gives NaN throughout: the simplex is
/// degenerate (SciPy's rule; SciPy estimates the condition number where this
/// computes it from the inverse).
///
/// SciPy keeps the inverse of `A` and multiplies by it. That is less
/// accurate: the inverse holds the large numbers of `U^-1`, which the solve
/// with `U` does not form. On a thin tetrahedron (a condition number of
/// 10^4), the inverse put a vertex of it at coordinates off by 1.8e-13 where
/// the solve put it exactly, and the value at a data point came out 1.5e-12
/// off. `M` and the solve with `U` are 1.5 times the multiplications of the
/// inverse, with no division and no reordering at a query.
fn push_transform<const D: usize>(points: &[[f64; D]], simplex: &[u32], out: &mut Vec<f64>) {
    let last = points[simplex[D] as usize];
    let mut a = [[0.0; D]; D];
    for (j, &v) in simplex[..D].iter().enumerate() {
        let vertex = points[v as usize];
        for (k, row) in a.iter_mut().enumerate() {
            row[j] = vertex[k] - last[k];
        }
    }
    let norm = norm1(&a);
    let usable = decompose(a).and_then(|(lu, order)| {
        let m = forward(&lu, &order);
        let inverse = back_all(&lu, &m);
        let finite = inverse
            .iter()
            .flatten()
            .chain(m.iter().flatten())
            .all(|c| c.is_finite());
        (finite && 1.0 / (norm * norm1(&inverse)) >= RCOND_LIMIT).then_some((lu, m))
    });
    match usable {
        Some((mut lu, m)) => {
            for row in &m {
                out.extend_from_slice(row);
            }
            for (i, row) in lu.iter_mut().enumerate() {
                row[i] = 1.0 / row[i];
                out.extend_from_slice(row);
            }
            out.extend_from_slice(&last);
        }
        None => out.extend(core::iter::repeat_n(f64::NAN, 2 * D * D + D)),
    }
}

/// `M = L^-1 P`, column by column: `L z = P e_k` forward.
fn forward<const D: usize>(lu: &[[f64; D]; D], order: &[usize; D]) -> [[f64; D]; D] {
    let mut m = [[0.0; D]; D];
    for k in 0..D {
        let mut z = [0.0; D];
        for i in 0..D {
            let mut sum = f64::from(u8::from(order[i] == k));
            for j in 0..i {
                sum -= lu[i][j] * z[j];
            }
            z[i] = sum;
        }
        for (row, value) in m.iter_mut().zip(z) {
            row[k] = value;
        }
    }
    m
}

/// `U^-1 M`, the inverse of `A`, column by column: `U x = M e_k` back.
fn back_all<const D: usize>(lu: &[[f64; D]; D], m: &[[f64; D]; D]) -> [[f64; D]; D] {
    let mut inverse = [[0.0; D]; D];
    for k in 0..D {
        let mut x = [0.0; D];
        for i in (0..D).rev() {
            let mut sum = m[i][k];
            for j in i + 1..D {
                sum -= lu[i][j] * x[j];
            }
            x[i] = sum / lu[i][i];
        }
        for (row, value) in inverse.iter_mut().zip(x) {
            row[k] = value;
        }
    }
    inverse
}

/// The largest sum of the magnitudes of a column.
fn norm1<const D: usize>(a: &[[f64; D]; D]) -> f64 {
    (0..D)
        .map(|j| a.iter().map(|row| row[j].abs()).sum::<f64>())
        .fold(0.0, f64::max)
}

/// The LU decomposition of `a` with partial pivoting, in place, and the row of
/// `a` each row came from; none when a pivot is zero.
fn decompose<const D: usize>(mut a: [[f64; D]; D]) -> Option<([[f64; D]; D], [usize; D])> {
    let mut order: [usize; D] = core::array::from_fn(|k| k);
    for col in 0..D {
        let pivot = (col..D).max_by(|&i, &j| a[i][col].abs().total_cmp(&a[j][col].abs()))?;
        if a[pivot][col] == 0.0 {
            return None;
        }
        a.swap(col, pivot);
        order.swap(col, pivot);
        let (above, below) = a.split_at_mut(col + 1);
        let top = &above[col];
        for row in below {
            let factor = row[col] / top[col];
            row[col] = factor;
            for (entry, &up) in row[col + 1..].iter_mut().zip(&top[col + 1..]) {
                *entry -= factor * up;
            }
        }
    }
    Some((a, order))
}

/// The barycentric coordinates of `x` from a transform: the first `D` solve
/// `U c = M (x - last)`, and the last is one less their sum (summed in order,
/// as SciPy sums them).
fn barycentric<const D: usize>(transform: &[f64], x: &[f64; D]) -> ([f64; D], f64) {
    let (m, rest) = transform.split_at(D * D);
    let (u, origin) = rest.split_at(D * D);
    let y: [f64; D] = core::array::from_fn(|k| x[k] - origin[k]);
    let mut c: [f64; D] =
        core::array::from_fn(|i| m[D * i..D * i + D].iter().zip(&y).map(|(a, b)| a * b).sum());
    for i in (0..D).rev() {
        let row = &u[D * i..D * i + D];
        let mut sum = c[i];
        for j in i + 1..D {
            sum -= row[j] * c[j];
        }
        c[i] = sum * row[i];
    }
    let mut last = 1.0;
    for value in &c {
        last -= value;
    }
    (c, last)
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
        let x = self.scaling.map_or(*point, |scaling| scaling.apply(point));
        if self.fully_outside(&x) {
            return None;
        }
        let (simplex, c, last) = self.locate(&x, self.start(&x))?;
        let vertices = &self.simplices[simplex as usize * (D + 1)..][..D + 1];
        // SciPy's order: from zero, the vertices in turn.
        let mut value = V::zero();
        for (k, &v) in vertices.iter().enumerate() {
            let weight = if k < D { c[k] } else { last };
            value = value.axpy(weight, self.values[v as usize]);
        }
        Some(value)
    }
}

/// The stored form of the `serde` feature.
#[cfg(feature = "serde")]
mod stored {
    use super::{too_many, LinearNDInterpolator, Scaling, NO_NEIGHBOR};
    use crate::{validate, InterpolateError, Value};
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
    }

    impl<const D: usize, V: Value<f64>> Serialize for LinearNDInterpolator<D, V> {
        fn serialize<S: Serializer>(&self, to: S) -> Result<S::Ok, S::Error> {
            Form {
                version: VERSION,
                points: Cow::Borrowed(&self.points),
                values: Cow::Borrowed(&self.values),
                simplices: Cow::Borrowed(&self.simplices),
                neighbors: Cow::Borrowed(&self.neighbors),
                scaling: self.scaling,
            }
            .serialize(to)
        }
    }

    /// The points and the values pass the checks of
    /// [`LinearNDInterpolator::new`] up to the triangulation. The simplices
    /// and their neighbors are checked to be whole and to name points and
    /// simplices that exist; that they are a Delaunay triangulation is not
    /// checked again (a walk over a broken triangulation ends at its limit,
    /// and the search over every simplex follows).
    impl<'de, const D: usize, V: Value<f64>> Deserialize<'de> for LinearNDInterpolator<D, V> {
        fn deserialize<D2: Deserializer<'de>>(from: D2) -> Result<Self, D2::Error> {
            let form = Form::<D, V>::deserialize(from)?;
            crate::stored::version(form.version, VERSION)?;
            let points = form.points.into_owned();
            let values = form.values.into_owned();
            let checked = || {
                validate::length(&values, points.len())?;
                validate::enough(&points, D + 1)?;
                if u32::try_from(points.len()).is_err() {
                    return Err(too_many());
                }
                let finite = |point: &[f64; D]| point.iter().all(|c| c.is_finite());
                if let Some(index) = validate::first_failing(&points, finite) {
                    return Err(InterpolateError::NonFinite { index });
                }
                if let Some(index) = validate::first_failing(&values, V::all_finite) {
                    return Err(InterpolateError::NonFinite { index });
                }
                Ok(())
            };
            checked().map_err(de::Error::custom)?;
            let (simplices, neighbors) = (form.simplices.into_owned(), form.neighbors.into_owned());
            let count = simplices.len() / (D + 1);
            let whole = simplices.len() % (D + 1) == 0
                && neighbors.len() == simplices.len()
                && (D == 0 || count > 0);
            if !whole {
                return Err(de::Error::custom(
                    "the simplices and their neighbors do not have D + 1 entries each",
                ));
            }
            if simplices.iter().any(|&v| v as usize >= points.len()) {
                return Err(de::Error::custom(
                    "a simplex names a point that does not exist",
                ));
            }
            if neighbors
                .iter()
                .any(|&n| n != NO_NEIGHBOR && n as usize >= count)
            {
                return Err(de::Error::custom(
                    "a neighbor names a simplex that does not exist",
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
            LinearNDInterpolator::assemble(points, values, simplices, neighbors, form.scaling)
                .map_err(de::Error::custom)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{LinearNDInterpolator, Walked, EPS};

    /// Points of the cube, many on its faces, the corners of a tetrahedron
    /// last: from some simplices a walk to a data point leaves the hull by
    /// rounding (found by a search over 20 000 such sets, 539 of which do
    /// this; the first set the property test found stopped doing it when the
    /// coordinates came from an LU decomposition instead of the inverse).
    fn faces() -> Vec<[f64; 3]> {
        let mut points = vec![
            [0.5840273880650761, 0.0, 0.29776683281999905],
            [0.7887697228607375, 0.0, 0.0],
            [0.7876449995383649, 0.8482870597656559, 0.0],
            [0.11013092240014277, 0.5642652099978088, 0.766618708512923],
            [0.0, 0.12887281267446582, 0.0],
            [0.06562029169092043, 0.7760534993474162, 0.7190231567724078],
            [0.0, 0.39189613554591596, 0.43322302906748944],
            [0.0, 0.84377177832942, 0.0],
            [0.015872575120935206, 0.0, 0.13668197412623584],
            [0.8962620857160368, 0.0, 0.2963248552207749],
            [0.0, 0.18161392058810888, 0.3633748190675532],
            [0.3993336928819621, 0.0, 0.0],
            [0.46366610551977316, 0.0, 0.7452603066909818],
            [0.4403637067996585, 0.0, 0.0],
            [0.759357153932266, 0.0, 0.20163702652700322],
            [0.0, 0.1242513710432015, 0.5428371255840139],
            [0.9972511880636947, 0.5069192550161593, 0.0],
            [0.925627308946488, 0.1844727134377805, 0.977486742272565],
            [0.0, 0.0571139320833719, 0.7637884613349173],
            [0.0, 0.0, 0.0],
            [0.0, 0.4044753825230042, 0.0],
            [0.0, 0.2881970878997515, 0.9808982126587361],
            [0.7724729514036233, 0.0, 0.8768541173058809],
            [0.00420881472699941, 0.0, 0.35924430590745415],
            [0.5430749905749186, 0.8800838478470157, 0.1601621397566283],
            [0.3592201432304133, 0.0, 0.0],
            [0.7624082887132098, 0.45579881603055306, 0.0],
        ];
        points.extend([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]);
        points
    }

    #[test]
    fn every_data_point_is_found_from_every_simplex() {
        let points = faces();
        let linear = LinearNDInterpolator::new(&points, &vec![0.0; points.len()]).expect("valid");
        let mut outside = 0;
        for start in 0..linear.simplex_count() as u32 {
            for (i, point) in points.iter().enumerate() {
                if matches!(linear.walk(point, start), Walked::Outside(..)) {
                    outside += 1;
                }
                assert!(
                    linear.locate(point, start).is_some(),
                    "point {i} from {start}"
                );
            }
        }
        // Some of these walks leave the hull by rounding, and are taken back.
        assert!(outside > 0);
    }

    #[test]
    fn the_middle_of_each_cell_is_in_the_simplex_its_walk_starts_from() {
        let mut state = 7_u64;
        let mut next = || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (state >> 11) as f64 / (1_u64 << 53) as f64
        };
        let points: Vec<[f64; 2]> = (0..500).map(|_| [next(), 3.0 * next()]).collect();
        let linear = LinearNDInterpolator::new(&points, &vec![0.0; 500]).expect("valid");
        let per_axis = linear.per_axis;
        assert_eq!(linear.cells.len(), per_axis * per_axis);
        let width = [0, 1].map(|d| (linear.upper[d] - linear.lower[d]) / per_axis as f64);
        let mut inside = 0;
        for i in 0..per_axis {
            for j in 0..per_axis {
                let middle = [
                    linear.lower[0] + (i as f64 + 0.5) * width[0],
                    linear.lower[1] + (j as f64 + 0.5) * width[1],
                ];
                let start = linear.start(&middle);
                assert_eq!(start, linear.cells[i * per_axis + j]);
                if linear.search(&middle).is_some() {
                    let (c, last) = linear.barycentric(start, &middle);
                    let all = [c[0], c[1], last];
                    assert!(all.iter().all(|&ck| ck >= -EPS), "cell {i}, {j}");
                    inside += 1;
                }
            }
        }
        assert!(inside > per_axis * per_axis / 2, "{inside}");
        // A walk across the box arrives by itself: from the start of the
        // first cell to the middle of a triangle at the point farthest from it.
        let corner = points.iter().copied().fold([f64::MIN; 2], |a, p| {
            if p[0] + p[1] > a[0] + a[1] {
                p
            } else {
                a
            }
        });
        let (far, ..) = linear.search(&corner).expect("a data point");
        let vertices = &linear.simplices[far as usize * 3..far as usize * 3 + 3];
        let middle =
            [0, 1].map(|d| vertices.iter().map(|&v| points[v as usize][d]).sum::<f64>() / 3.0);
        assert!(
            matches!(linear.walk(&middle, linear.cells[0]), Walked::Found(found, ..) if found == far)
        );
    }
}
