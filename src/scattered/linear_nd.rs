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
///   This crate's triangulation splits such a cell by a rule of the points'
///   indices; SciPy's (Qhull) splits it another way. Inside such a cell the
///   two interpolants differ; elsewhere they agree to rounding. Points on a
///   grid are such points.
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
    /// For each simplex, the inverse of the matrix whose columns are its first
    /// `D` vertices less the last one (`D * D` numbers, row by row), then the
    /// last vertex (`D` numbers). All NaN for a degenerate simplex.
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
        let mut transforms = Vec::with_capacity(simplices.len() / (D + 1) * (D * D + D));
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

    /// The transform of a simplex: the inverse matrix, then the last vertex.
    fn transform(&self, simplex: u32) -> &[f64] {
        let stride = D * D + D;
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
    /// at. SciPy takes such a point as outside.
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

/// Appends the transform of a simplex: the inverse of the matrix `A` whose
/// column `j` is vertex `j` less the last vertex, row by row, then the last
/// vertex. The barycentric coordinates of `x` are `A^-1 (x - last)`, and one
/// less their sum. A matrix that is singular, or whose reciprocal condition
/// number in the 1-norm is below [`RCOND_LIMIT`], gives NaN throughout: the
/// simplex is degenerate (SciPy's rule, which estimates the condition number
/// where this computes it).
fn push_transform<const D: usize>(points: &[[f64; D]], simplex: &[u32], out: &mut Vec<f64>) {
    let last = points[simplex[D] as usize];
    // The matrix and the identity beside it, reduced together.
    let mut a = [[0.0; D]; D];
    for (j, &v) in simplex[..D].iter().enumerate() {
        let vertex = points[v as usize];
        for k in 0..D {
            a[k][j] = vertex[k] - last[k];
        }
    }
    let norm = norm1(&a);
    match invert(a) {
        Some(inverse) if 1.0 / (norm * norm1(&inverse)) >= RCOND_LIMIT => {
            for row in &inverse {
                out.extend_from_slice(row);
            }
            out.extend_from_slice(&last);
        }
        _ => out.extend(core::iter::repeat_n(f64::NAN, D * D + D)),
    }
}

/// The largest sum of the magnitudes of a column.
fn norm1<const D: usize>(a: &[[f64; D]; D]) -> f64 {
    (0..D)
        .map(|j| (0..D).map(|k| a[k][j].abs()).sum::<f64>())
        .fold(0.0, f64::max)
}

/// The inverse of `a`, by elimination with partial pivoting; none when a pivot
/// is zero or a number is not finite.
fn invert<const D: usize>(mut a: [[f64; D]; D]) -> Option<[[f64; D]; D]> {
    let mut inverse = [[0.0; D]; D];
    for (k, row) in inverse.iter_mut().enumerate() {
        row[k] = 1.0;
    }
    for col in 0..D {
        let pivot = (col..D).max_by(|&i, &j| a[i][col].abs().total_cmp(&a[j][col].abs()))?;
        if a[pivot][col] == 0.0 {
            return None;
        }
        a.swap(col, pivot);
        inverse.swap(col, pivot);
        let p = a[col][col];
        for j in 0..D {
            a[col][j] /= p;
            inverse[col][j] /= p;
        }
        for i in 0..D {
            if i != col {
                let f = a[i][col];
                if f != 0.0 {
                    for j in 0..D {
                        a[i][j] -= f * a[col][j];
                        inverse[i][j] -= f * inverse[col][j];
                    }
                }
            }
        }
    }
    let finite = inverse.iter().flatten().all(|c| c.is_finite());
    finite.then_some(inverse)
}

/// The barycentric coordinates of `x` from a transform, in SciPy's order of
/// operations (`_barycentric_coordinates`): each of the first `D` a sum over
/// the axes, and the last one less their sum.
fn barycentric<const D: usize>(transform: &[f64], x: &[f64; D]) -> ([f64; D], f64) {
    let origin = &transform[D * D..];
    let mut c = [0.0; D];
    let mut last = 1.0;
    for i in 0..D {
        let mut sum = 0.0;
        for j in 0..D {
            sum += transform[D * i + j] * (x[j] - origin[j]);
        }
        c[i] = sum;
        last -= sum;
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
    /// last: from one simplex a vertex of it was measured at -2.8e-14, beyond
    /// a face of the hull (found by the property test of the crate's tests).
    fn faces() -> Vec<[f64; 3]> {
        let mut points = vec![
            [0.01855029732426409, 0.4189100396648029, 0.9878574621282309],
            [0.36193556987925546, 0.0, 0.0],
            [0.7843633632952304, 0.653731011210593, 0.33376661413161524],
            [0.0, 0.9709120991315361, 0.8511377088721644],
            [0.8315244649727277, 0.9290769372512558, 0.5567199593697773],
            [0.7806416091605959, 0.7537200910363946, 0.0],
            [0.8732260078640668, 0.18876222144571253, 0.7324975500244549],
            [0.9563382026194267, 0.24790925355179655, 0.8152721782673668],
            [0.9469698999924424, 0.8620380150821337, 0.478236157411085],
            [0.0, 0.4436120884715504, 0.0],
            [0.898631835153527, 0.10775013641625633, 0.25155815584246183],
            [0.0, 0.0, 0.07626251099751173],
            [0.1836053525164637, 0.12840909092806518, 0.9230920256445325],
            [0.0, 0.0, 0.0],
            [0.8728169810583641, 0.0, 0.6330264227718598],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.171654320720206],
            [0.5073135347364379, 0.651661049074369, 0.3530479578165536],
            [0.054945725420766214, 0.0, 0.38628792332297174],
            [0.15096494212266162, 0.8966487811053174, 0.0],
            [0.7688916685756781, 0.1733312401021074, 0.9846908700392136],
            [0.9068401373632985, 0.3219005406327637, 0.7122078440897608],
            [0.0, 0.3275146390404187, 0.3553889670200437],
            [0.5432056873632348, 0.8867315711200432, 0.20128715609062459],
            [0.3179714535738961, 0.4041479123423767, 0.7199495947101067],
            [0.758088262493234, 0.5563597536656779, 0.0],
            [0.0, 0.8957314652636468, 0.29646767136394586],
            [0.506951858029448, 0.8298281977491287, 0.9652216748206536],
            [0.0, 0.0, 0.8213083640932333],
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
