//! What `LinearNDInterpolator` and `CloughTocher2DInterpolator` share: the
//! Delaunay triangulation of the points, the transform of each simplex, the
//! cells where a walk starts, and the simplex that holds a query (M4-2: one
//! part, so that what is put right in it is put right for both).

use super::Scaling;
use crate::error::{DegeneracyKind, ParamReason};
use crate::{validate, InterpolateError, Value};
use std::borrow::Cow;

/// The neighbor across a face of the hull: no simplex. `convx` writes it so and
/// promises it in its documentation; it has no constant of its own for it.
pub(super) const NO_NEIGHBOR: u32 = u32::MAX;

/// How far below zero (or above one) a barycentric coordinate may be for the
/// point to be inside the simplex: SciPy's `eps`, with a simplex tolerance of
/// one.
const EPS: f64 = 100.0 * f64::EPSILON;

/// The wider leeway toward a degenerate simplex, in the search over every
/// simplex: SciPy's `eps_broad`, the square root of `f64::EPSILON` (2^-26).
const EPS_BROAD: f64 = 1.490_116_119_384_765_6e-8;

/// A simplex whose matrix has a reciprocal condition number below this is
/// degenerate: it has no transform, and no point is found inside it by its own
/// coordinates. SciPy's limit, `1000 * eps`.
const RCOND_LIMIT: f64 = 1000.0 * f64::EPSILON;

/// The Delaunay triangulation of points in `D` dimensions, with what a query
/// is located by.
#[derive(Debug, Clone)]
pub(super) struct Mesh<const D: usize> {
    /// The points, as given, moved and scaled when the interpolant is
    /// rescaled. Duplicates stay, and are no vertex.
    pub(super) points: Vec<[f64; D]>,
    /// The `D + 1` vertices of each simplex: places in `points`.
    pub(super) simplices: Vec<u32>,
    /// The `D + 1` simplices across the faces of each simplex: the one across
    /// the face opposite each vertex, or [`NO_NEIGHBOR`] on the hull.
    pub(super) neighbors: Vec<u32>,
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
    pub(super) scaling: Option<Scaling<D>>,
    /// The multiplier of [`EPS`]: SciPy's `simplex_tolerance`.
    pub(super) tolerance: f64,
    /// [`EPS`] times the tolerance.
    eps: f64,
}

impl<const D: usize> Mesh<D> {
    /// The checks of the constructors before the triangulation, in their
    /// order: the lengths, the number of points (`D + 1`), the limit of their
    /// number, the coordinates, the values.
    pub(super) fn check<V: Value<f64>>(
        points: &[[f64; D]],
        values: &[V],
    ) -> Result<(), InterpolateError> {
        validate::length(values, points.len())?;
        validate::enough(points, D + 1)?;
        if u32::try_from(points.len()).is_err() {
            return Err(too_many());
        }
        let finite = |point: &[f64; D]| point.iter().all(|c| c.is_finite());
        if let Some(index) = validate::first_failing(points, finite) {
            return Err(InterpolateError::NonFinite { index });
        }
        if let Some(index) = validate::first_failing(values, V::all_finite) {
            return Err(InterpolateError::NonFinite { index });
        }
        Ok(())
    }

    /// The triangulation of points that passed [`check`](Mesh::check),
    /// rescaled first when asked.
    pub(super) fn build(
        points: Cow<'_, [[f64; D]]>,
        rescale: bool,
    ) -> Result<Self, InterpolateError> {
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
        Self::assemble(points, simplices, neighbors, scaling)
    }

    /// The mesh of a triangulation: the transforms of the simplices, the box
    /// of the points, and a start for the walk in each cell of it.
    pub(super) fn assemble(
        points: Vec<[f64; D]>,
        simplices: Vec<u32>,
        neighbors: Vec<u32>,
        scaling: Option<Scaling<D>>,
    ) -> Result<Self, InterpolateError> {
        // As many simplices as points or many more: by `reserve`.
        let mut transforms = Vec::new();
        crate::reserve::room(&mut transforms, simplices.len() / (D + 1) * (2 * D * D + D))?;
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
        let mut mesh = Self {
            points,
            simplices,
            neighbors,
            transforms,
            lower,
            upper,
            cells: Vec::new(),
            per_axis: 0,
            cell_scale: [0.0; D],
            scaling,
            tolerance: 1.0,
            eps: EPS,
        };
        mesh.place_cells();
        Ok(mesh)
    }

    /// Cuts the box into at most twice as many cells as there are points,
    /// and walks to the middle of each from the start of the cell before it:
    /// the last simplex of that walk (the one that holds the middle, or one
    /// on the hull toward it) is the start of the cell.
    ///
    /// The number of cells along each axis is the largest `p` with `p^D` no
    /// more than twice the number of points. A rounding up of `n^(1/D)` made
    /// `p` two or more whatever the dimension, and `2^D` cells for a few
    /// points: a million at 20 dimensions for 62 simplices (0.37 s to build;
    /// the review of #183 measured 3.5 s at 22). In high dimensions `p` is 1,
    /// and a walk starts from where the walk to the middle of the box ended.
    /// Twice rather than once the number of points (measured in three trees
    /// in turn, three rounds): evaluation 5 to 15% faster, construction 0 to
    /// 5% slower; at 100 points in 3-D the cells are 125 and not 64.
    fn place_cells(&mut self) {
        let count = self.simplex_count();
        if D == 0 || count == 0 {
            return;
        }
        let per_axis = per_axis(2 * self.points.len(), D);
        // Every axis has a range above zero: the points span D dimensions.
        let width: [f64; D] =
            core::array::from_fn(|d| (self.upper[d] - self.lower[d]) / per_axis as f64);
        self.cell_scale = core::array::from_fn(|d| 1.0 / width[d]);
        self.per_axis = per_axis;
        // No more than twice the number of points (`per_axis`).
        let total = per_axis.pow(D as u32);
        let mut cells = vec![0_u32; total];
        let mut start = 0;
        let mut at = [0_usize; D];
        // The cells in the order of a snake: the last axis fastest, each axis
        // turning back where it ends instead of jumping to its start, so that
        // each walk starts from the cell beside it (a jump back across a row
        // was a walk across the whole triangulation, once a row).
        let mut forward = [true; D];
        for _ in 0..total {
            let middle: [f64; D] =
                core::array::from_fn(|d| self.lower[d] + (at[d] as f64 + 0.5) * width[d]);
            start = match self.walk(&middle, start) {
                Walked::Found(simplex, ..)
                | Walked::Outside(simplex, _)
                | Walked::Lost(simplex) => simplex,
            };
            let index = at.iter().fold(0, |cell, &place| cell * per_axis + place);
            cells[index] = start;
            for d in (0..D).rev() {
                if forward[d] && at[d] + 1 < per_axis {
                    at[d] += 1;
                    break;
                }
                if !forward[d] && at[d] > 0 {
                    at[d] -= 1;
                    break;
                }
                forward[d] = !forward[d];
            }
        }
        self.cells = cells;
    }

    /// The simplex the walk to `x` starts from: that of the cell `x` is in
    /// (the nearest cell, for a point just outside the box).
    fn start(&self, x: &[f64; D]) -> u32 {
        let mut cell = 0;
        for (d, coordinate) in x.iter().enumerate() {
            // A cast saturates (below zero to zero), as `floor` and the clamp
            // would; the coordinate is finite here.
            let along = ((coordinate - self.lower[d]) * self.cell_scale[d]) as usize;
            let place = along.min(self.per_axis - 1);
            cell = cell * self.per_axis + place;
        }
        self.cells[cell]
    }

    /// The number of simplices.
    pub(super) fn simplex_count(&self) -> usize {
        self.simplices.len() / (D + 1)
    }

    /// The `D + 1` vertices of a simplex.
    #[inline]
    pub(super) fn vertices(&self, simplex: u32) -> &[u32] {
        &self.simplices[simplex as usize * (D + 1)..][..D + 1]
    }

    /// Sets SciPy's `simplex_tolerance`; see
    /// `LinearNDInterpolator::with_simplex_tolerance`.
    pub(super) fn set_tolerance(&mut self, tolerance: f64) -> Result<(), InterpolateError> {
        let reason = if !tolerance.is_finite() {
            Some(ParamReason::NonFinite)
        } else if tolerance <= 0.0 {
            Some(ParamReason::NotPositive)
        } else {
            None
        };
        if let Some(reason) = reason {
            return Err(InterpolateError::InvalidParameter {
                name: "simplex_tolerance",
                reason,
            });
        }
        self.tolerance = tolerance;
        self.eps = EPS * tolerance;
        Ok(())
    }

    /// The simplex that holds `point`, with its barycentric coordinates
    /// there: none for a point that is not finite or is outside the hull.
    /// The point is changed as the points were, first.
    #[inline]
    pub(super) fn find(&self, point: &[f64; D]) -> Option<(u32, [f64; D], f64)> {
        if !point.iter().all(|c| c.is_finite()) {
            return None;
        }
        let x = self.scaling.map_or(*point, |scaling| scaling.apply(point));
        if self.fully_outside(&x) {
            return None;
        }
        self.locate(&x, self.start(&x))
    }

    /// The mesh of a stored form whose points passed [`check`](Mesh::check):
    /// the simplices and their neighbors whole, naming points and simplices
    /// that exist, each simplex `D + 1` distinct points, every neighbor
    /// having the simplex among its neighbors (the review of #183), and a
    /// usable change of coordinates and tolerance.
    #[cfg(feature = "serde")]
    pub(super) fn from_stored(
        points: Vec<[f64; D]>,
        simplices: Vec<u32>,
        neighbors: Vec<u32>,
        scaling: Option<Scaling<D>>,
        tolerance: f64,
    ) -> Result<Self, String> {
        let refuse = |why: &str| Err::<Self, String>(why.to_string());
        let count = simplices.len() / (D + 1);
        let whole = simplices.len().is_multiple_of(D + 1)
            && neighbors.len() == simplices.len()
            && (D == 0 || count > 0);
        if !whole {
            return refuse("the simplices and their neighbors do not have D + 1 entries each");
        }
        if simplices.iter().any(|&v| v as usize >= points.len()) {
            return refuse("a simplex names a point that does not exist");
        }
        if neighbors
            .iter()
            .any(|&n| n != NO_NEIGHBOR && n as usize >= count)
        {
            return refuse("a neighbor names a simplex that does not exist");
        }
        fn row(list: &[u32], s: usize, width: usize) -> &[u32] {
            &list[s * width..(s + 1) * width]
        }
        for s in 0..count {
            let vertices = row(&simplices, s, D + 1);
            let repeated = (0..=D).any(|i| vertices[..i].contains(&vertices[i]));
            if repeated {
                return refuse("a simplex names a point twice");
            }
            for &t in row(&neighbors, s, D + 1) {
                if t != NO_NEIGHBOR && !row(&neighbors, t as usize, D + 1).contains(&(s as u32)) {
                    return refuse("a simplex is a neighbor of one that is not its neighbor");
                }
            }
        }
        if let Some(scaling) = &scaling {
            let usable = scaling.offset.iter().all(|c| c.is_finite())
                && scaling.scale.iter().all(|c| c.is_finite() && *c > 0.0);
            if !usable {
                return refuse(
                    "the change of coordinates is not finite, or its scale is not above zero",
                );
            }
        }
        let mut mesh =
            Self::assemble(points, simplices, neighbors, scaling).map_err(|e| e.to_string())?;
        mesh.set_tolerance(tolerance).map_err(|e| e.to_string())?;
        Ok(mesh)
    }

    /// The transform of a simplex: `L^-1 P`, `U`, the last vertex.
    fn transform(&self, simplex: u32) -> &[f64] {
        let stride = 2 * D * D + D;
        let start = simplex as usize * stride;
        &self.transforms[start..start + stride]
    }

    /// The barycentric coordinates of `x` in `simplex`: those of the first
    /// `D` vertices, and that of the last. NaN for a degenerate simplex.
    pub(super) fn barycentric(&self, simplex: u32, x: &[f64; D]) -> ([f64; D], f64) {
        barycentric(self.transform(simplex), x)
    }

    fn neighbor(&self, simplex: u32, face: usize) -> u32 {
        self.neighbors[simplex as usize * (D + 1) + face]
    }

    /// Whether `x` is outside the box of the points by more than [`EPS`] on
    /// some axis: then it is outside the hull. SciPy's
    /// `_is_point_fully_outside`.
    fn fully_outside(&self, x: &[f64; D]) -> bool {
        (0..D).any(|d| x[d] < self.lower[d] - self.eps || x[d] > self.upper[d] + self.eps)
    }

    /// The walk from `start` toward `x`: across the face opposite the first
    /// vertex whose coordinate is below `-EPS` (SciPy's
    /// `_find_simplex_directed`), as many steps at most as there are
    /// simplices. SciPy stops at one plus a quarter of them and looks at
    /// every simplex; here the start is near, and a walk that long is one
    /// that has gone wrong either way.
    fn walk(&self, x: &[f64; D], start: u32) -> Walked<D> {
        let mut simplex = start;
        for _ in 0..self.simplex_count() {
            let (c, last) = self.barycentric(simplex, x);
            let mut inside = true;
            let mut across = None;
            for (k, ck) in c.iter().copied().chain([last]).enumerate() {
                if ck < -self.eps {
                    across = Some((k, ck));
                    break;
                }
                // Above one, or NaN: a degenerate simplex is in the way.
                if ck.is_nan() || ck > 1.0 + self.eps {
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
                .all(|(k, ck)| -low(k) <= ck && ck <= 1.0 + self.eps)
        };
        for simplex in 0..self.simplex_count() as u32 {
            if !self.transform(simplex)[0].is_nan() {
                let (c, last) = self.barycentric(simplex, x);
                if within(&c, last, &|_| self.eps) {
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
                        self.eps
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

/// The largest `p` with `p^dim` no more than `n`, and 1 at least.
fn per_axis(n: usize, dim: usize) -> usize {
    let fits = |p: usize| {
        u32::try_from(dim)
            .ok()
            .and_then(|dim| p.checked_pow(dim))
            .is_some_and(|cells| cells <= n)
    };
    let mut p = ((n as f64).powf(1.0 / dim as f64).floor() as usize).max(1);
    // The root in floating point can be a little off either way.
    while p > 1 && !fits(p) {
        p -= 1;
    }
    while fits(p + 1) {
        p += 1;
    }
    p
}

/// The simplices of the Delaunay triangulation of the points, and their
/// neighbors, `D + 1` of each a simplex. Without dimensions there is none.
///
/// Every distinct point is a vertex, and a point given more than once is
/// represented by the smallest index of the points equal to it (`convx`'s
/// documentation of `DelaunayTriangulation`): so a vertex has the value of
/// the first of equal points, and the others are no vertex.
///
/// The number of simplices can be far more than the number of points in high
/// dimensions, so the copies are made by `reserve`.
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
    let mut vertices = Vec::new();
    let mut neighbors = Vec::new();
    crate::reserve::room(&mut vertices, simplices.len() * (D + 1))?;
    crate::reserve::room(&mut neighbors, simplices.len() * (D + 1))?;
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

#[cfg(test)]
mod tests {
    use super::{Mesh, Walked, EPS};
    use std::borrow::Cow;

    /// The mesh of points as they are.
    fn mesh_of<const D: usize>(points: &[[f64; D]]) -> Mesh<D> {
        Mesh::build(Cow::Borrowed(points), false).expect("valid")
    }

    /// A small generator, so that the sets are the same on every run.
    fn generator(seed: u64) -> impl FnMut() -> f64 {
        let mut state = seed;
        move || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (state >> 11) as f64 / (1_u64 << 53) as f64
        }
    }

    #[test]
    fn every_data_point_is_found_from_every_simplex() {
        // Sets of points of the cube, three in ten coordinates on a face, and
        // the corners of a tetrahedron: from some simplices a walk to a data
        // point leaves the hull by rounding (539 of 20 000 such sets did; the
        // first set the property test found stopped doing it when the
        // coordinates came from an LU decomposition instead of the inverse,
        // which is why the test goes through a series of sets).
        let mut next = generator(11);
        let mut outside = 0;
        for _ in 0..200 {
            let n = 4 + (next() * 30.0) as usize;
            let mut points: Vec<[f64; 3]> = (0..n)
                .map(|_| {
                    core::array::from_fn(|_| {
                        let c = next();
                        if next() < 0.3 {
                            0.0
                        } else {
                            c
                        }
                    })
                })
                .collect();
            points.extend([
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ]);
            let linear = mesh_of(&points);
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
        }
        // Some of these walks leave the hull by rounding, and are taken back.
        assert!(outside > 0);
    }

    #[test]
    fn the_cells_are_no_more_than_twice_the_points_in_any_dimension() {
        // The largest p with p^D no more than n.
        assert_eq!(super::per_axis(19, 16), 1);
        assert_eq!(super::per_axis(7, 3), 1);
        assert_eq!(super::per_axis(8, 3), 2);
        assert_eq!(super::per_axis(26, 3), 2);
        assert_eq!(super::per_axis(27, 3), 3);
        assert_eq!(super::per_axis(100_000, 2), 316);
        assert_eq!(super::per_axis(1000, 64), 1);
        assert_eq!(super::per_axis(usize::MAX, 2), 4_294_967_295);
        // 23 points in 20 dimensions: one cell, 46 being under 2^20 (it
        // was 2^20 cells, and 0.37 s).
        let mut next = generator(5);
        let points: Vec<[f64; 20]> = (0..23).map(|_| core::array::from_fn(|_| next())).collect();
        let linear = mesh_of(&points);
        assert_eq!(linear.cells.len(), 1);
        for point in &points {
            assert!(linear.find(point).is_some());
        }
    }

    #[test]
    fn the_middle_of_each_cell_is_in_the_simplex_its_walk_starts_from() {
        let mut next = generator(7);
        let points: Vec<[f64; 2]> = (0..500).map(|_| [next(), 3.0 * next()]).collect();
        let linear = mesh_of(&points);
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
