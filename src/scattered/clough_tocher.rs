//! Piecewise cubic C1 interpolation on the Delaunay triangulation of points of
//! the plane: the Clough-Tocher scheme, with the gradients at the points from
//! a global estimate.

use super::mesh::{Mesh, NO_NEIGHBOR};
use crate::error::SolverError;
use crate::{InterpolateError, Interpolator, Value};
use std::borrow::Cow;

/// The relative residual at which the gradient estimate is taken as solved,
/// in every column of the values. The preconditioned conjugate gradients
/// bring it down about four times a step, and below 1e-16 within 26 steps on
/// every set tried (random points up to 10^5, coordinates a thousand times
/// apart, clustered points, two points 1e-12 apart, triangles 1e-9 thin):
/// 1e-14 is reached near the 20th. Against SciPy iterated to convergence the
/// values are then within 4.5e-13 of the largest value; with 1e-12 they were
/// within 3e-11.
const SOLVED: f64 = 1e-14;

/// The interpolant that is a cubic on each third of every triangle of the
/// Delaunay triangulation of the points, continuous with its first
/// derivatives (C1), for points of the plane.
///
/// SciPy's `CloughTocher2DInterpolator`. Each triangle is cut at its
/// centroid into three, and the interpolant is a cubic Bernstein-Bezier
/// patch on each part: it takes the data values at the points, and the
/// gradients there. The gradients are estimated from the data by the global
/// method of Nielson (1983) and Renka and Cline (1984): those that make the
/// interpolant bend least along the edges of the triangulation. The direction
/// in which the derivative across an edge is made linear is taken from the
/// centroids of the two triangles beside it, so that the interpolant does not
/// change under an affine change of coordinates (as in SciPy).
///
/// A query outside the convex hull of the points has no value, nor has a
/// query with a coordinate that is NaN or infinite. A value is a scalar or a
/// fixed-length array; each column has its gradients estimated on its own.
/// Coordinates are `f64`.
///
/// Points with the same coordinates are taken as one, with the value of the
/// one given first. The simplex of a query is found as
/// [`LinearNDInterpolator`](crate::LinearNDInterpolator) finds it, with the
/// same [`with_simplex_tolerance`](CloughTocher2DInterpolator::with_simplex_tolerance)
/// and [`rescaled`](CloughTocher2DInterpolator::rescaled).
///
/// # Differences from SciPy
///
/// - The gradients are the solution of the system of the estimate, solved
///   until its residual is below `1e-14` of its right-hand side in every
///   column. SciPy stops its Gauss-Seidel iteration at its `tol` (`1e-6`) or
///   `maxiter` (400), and only warns when it has not converged; its
///   gradients are that far from the solution, and its values with them.
///   There are no `tol` and `maxiter`, and a system that is not solved is an
///   error.
/// - Where four or more points lie on one circle, the triangulation, and so
///   the interpolant, can differ from SciPy's (as for
///   [`LinearNDInterpolator`](crate::LinearNDInterpolator)).
/// - Of points with the same coordinates, the value of the one given first
///   is taken. SciPy (Qhull) keeps one of them by its own choice.
///
/// # Examples
///
/// ```
/// use interpx::{CloughTocher2DInterpolator, InterpolateError, Interpolator};
///
/// # fn main() -> Result<(), InterpolateError> {
/// // f(x, y) = 1 + 2x + 3y at a few points: a plane comes back as it is.
/// let points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0], [0.4, 0.6]];
/// let values: Vec<f64> = points.iter().map(|p| 1.0 + 2.0 * p[0] + 3.0 * p[1]).collect();
/// let cubic = CloughTocher2DInterpolator::new(&points, &values)?;
/// let inside = cubic.eval(&[0.25, 0.5]).unwrap();
/// assert!((inside - 3.0).abs() < 1e-10);
/// assert_eq!(cubic.eval(&[1.5, 0.5]), None);
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone)]
pub struct CloughTocher2DInterpolator<V: Value<f64> = f64> {
    /// The triangulation of the points, and what a query is located by.
    mesh: Mesh<2>,
    /// The values, as given.
    values: Vec<V>,
    /// The gradient at each point, by `x` and by `y` (after the change of
    /// coordinates of a rescaled interpolant); zero at a point that is no
    /// vertex.
    gradients: Vec<[V; 2]>,
    /// For each triangle, SciPy's `g[k]`: how the direction across the edge
    /// opposite vertex `k` leans, from the centroid of the triangle beyond
    /// it (`-1/2` on the hull).
    lean: Vec<[f64; 3]>,
}

impl<V: Value<f64>> CloughTocher2DInterpolator<V> {
    /// Builds the interpolant from copies of the points and of their values.
    ///
    /// # Errors
    ///
    /// Those of [`LinearNDInterpolator::new`](crate::LinearNDInterpolator::new)
    /// in two dimensions, in their order, and then
    /// [`InterpolateError::SolverFailed`] with [`SolverError::NotConverged`]
    /// when the system of the gradient estimate is not solved.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{CloughTocher2DInterpolator, InterpolateError, Interpolator};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// // Two values at every point.
    /// let points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]];
    /// let values = [[0.0, 1.0], [1.0, 1.0], [1.0, 1.0], [2.0, 1.0]];
    /// let cubic = CloughTocher2DInterpolator::new(&points, &values)?;
    /// let [a, b] = cubic.eval(&[0.5, 0.5]).unwrap();
    /// assert!((a - 1.0).abs() < 1e-10 && (b - 1.0).abs() < 1e-10);
    ///
    /// let err = CloughTocher2DInterpolator::new(&[[0.0, 0.0], [1.0, 0.0]], &[1.0, 2.0]);
    /// assert!(matches!(err, Err(InterpolateError::InsufficientPoints { required: 3, actual: 2 })));
    /// # Ok(())
    /// # }
    /// ```
    pub fn new(points: &[[f64; 2]], values: &[V]) -> Result<Self, InterpolateError> {
        Self::build(Cow::Borrowed(points), Cow::Borrowed(values), false)
    }

    /// Builds the interpolant from points and values that are handed over,
    /// which are kept and not copied.
    ///
    /// # Errors
    ///
    /// As [`new`](CloughTocher2DInterpolator::new).
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{CloughTocher2DInterpolator, InterpolateError, Interpolator};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let points = vec![[0.0, 0.0], [2.0, 0.0], [0.0, 2.0]];
    /// let cubic = CloughTocher2DInterpolator::from_owned(points, vec![0.0, 2.0, 4.0])?;
    /// assert!((cubic.eval(&[1.0, 0.0]).unwrap() - 1.0).abs() < 1e-12);
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_owned(points: Vec<[f64; 2]>, values: Vec<V>) -> Result<Self, InterpolateError> {
        Self::build(Cow::Owned(points), Cow::Owned(values), false)
    }

    /// Builds the interpolant with both coordinates brought to the same range
    /// before the points are triangulated: SciPy's `rescale=True`, the change
    /// of [`LinearNDInterpolator::rescaled`](crate::LinearNDInterpolator::rescaled).
    ///
    /// # Errors
    ///
    /// As [`LinearNDInterpolator::rescaled`](crate::LinearNDInterpolator::rescaled),
    /// and then those of the gradient estimate.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{CloughTocher2DInterpolator, InterpolateError, Interpolator};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// // Positions in kilometres, heights in metres; a plane comes back.
    /// let points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1000.0], [1.0, 1000.0], [0.3, 400.0]];
    /// let values: Vec<f64> = points.iter().map(|p| p[0] + p[1] / 1000.0).collect();
    /// let rescaled = CloughTocher2DInterpolator::rescaled(&points, &values)?;
    /// assert!((rescaled.eval(&[0.5, 500.0]).unwrap() - 1.0).abs() < 1e-10);
    /// # Ok(())
    /// # }
    /// ```
    pub fn rescaled(points: &[[f64; 2]], values: &[V]) -> Result<Self, InterpolateError> {
        Self::build(Cow::Borrowed(points), Cow::Borrowed(values), true)
    }

    /// Builds as [`rescaled`](CloughTocher2DInterpolator::rescaled) does, from
    /// points and values that are handed over. The points are changed in
    /// place.
    ///
    /// # Errors
    ///
    /// As [`rescaled`](CloughTocher2DInterpolator::rescaled).
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{CloughTocher2DInterpolator, InterpolateError, Interpolator};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let points = vec![[0.0, 0.0], [1000.0, 0.0], [0.0, 1.0]];
    /// let cubic = CloughTocher2DInterpolator::rescaled_from_owned(points, vec![1.0, 2.0, 3.0])?;
    /// assert!((cubic.eval(&[0.0, 0.0]).unwrap() - 1.0).abs() < 1e-12);
    /// # Ok(())
    /// # }
    /// ```
    pub fn rescaled_from_owned(
        points: Vec<[f64; 2]>,
        values: Vec<V>,
    ) -> Result<Self, InterpolateError> {
        Self::build(Cow::Owned(points), Cow::Owned(values), true)
    }

    fn build(
        points: Cow<'_, [[f64; 2]]>,
        values: Cow<'_, [V]>,
        rescale: bool,
    ) -> Result<Self, InterpolateError> {
        Mesh::check(&points, &values)?;
        let mesh = Mesh::build(points, rescale)?;
        let values = values.into_owned();
        let gradients = estimate_gradients(&mesh, &values)?;
        let lean = leans(&mesh);
        Ok(Self {
            mesh,
            values,
            gradients,
            lean,
        })
    }

    /// Returns the number of points given. See
    /// [`CloughTocher2DInterpolator`] for an example.
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
    /// use interpx::{CloughTocher2DInterpolator, InterpolateError};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let cubic = CloughTocher2DInterpolator::new(&[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], &[5.0, 6.0, 7.0])?;
    /// assert_eq!(cubic.values(), &[5.0, 6.0, 7.0]);
    /// # Ok(())
    /// # }
    /// ```
    pub fn values(&self) -> &[V] {
        &self.values
    }

    /// Returns the points, as triangulated: as given, or after the change of
    /// coordinates of a [`rescaled`](CloughTocher2DInterpolator::rescaled)
    /// interpolant.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{CloughTocher2DInterpolator, InterpolateError};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let points = [[0.0, 0.0], [2.0, 0.0], [0.0, 4.0]];
    /// let cubic = CloughTocher2DInterpolator::new(&points, &[0.0; 3])?;
    /// assert_eq!(cubic.points(), &points);
    /// # Ok(())
    /// # }
    /// ```
    pub fn points(&self) -> &[[f64; 2]] {
        &self.mesh.points
    }

    /// Returns the estimated gradient at each point, by `x` and by `y`, in the
    /// order of the points given (after the change of coordinates of a
    /// rescaled interpolant; zero at a point given again, which is no
    /// vertex). SciPy's `grad`.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{CloughTocher2DInterpolator, InterpolateError};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// // f(x, y) = 2x - y: the estimate gives the gradient of a plane.
    /// let points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]];
    /// let values: Vec<f64> = points.iter().map(|p| 2.0 * p[0] - p[1]).collect();
    /// let cubic = CloughTocher2DInterpolator::new(&points, &values)?;
    /// let [dx, dy] = cubic.gradients()[0];
    /// assert!((dx - 2.0).abs() < 1e-10 && (dy + 1.0).abs() < 1e-10);
    /// # Ok(())
    /// # }
    /// ```
    pub fn gradients(&self) -> &[[V; 2]] {
        &self.gradients
    }

    /// Returns the number of triangles.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{CloughTocher2DInterpolator, InterpolateError};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let points = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.5, 0.5]];
    /// let cubic = CloughTocher2DInterpolator::new(&points, &[0.0; 5])?;
    /// assert_eq!(cubic.simplex_count(), 4);
    /// # Ok(())
    /// # }
    /// ```
    pub fn simplex_count(&self) -> usize {
        self.mesh.simplex_count()
    }

    /// Returns the same interpolant, with the tolerance by which a query is
    /// taken to be in a triangle multiplied by `tolerance`: SciPy's
    /// `simplex_tolerance` (SciPy 1.18), as
    /// [`LinearNDInterpolator::with_simplex_tolerance`](crate::LinearNDInterpolator::with_simplex_tolerance).
    ///
    /// # Errors
    ///
    /// [`InterpolateError::InvalidParameter`] for `simplex_tolerance` when
    /// `tolerance` is not finite or not above zero.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{CloughTocher2DInterpolator, InterpolateError, Interpolator};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let cubic = CloughTocher2DInterpolator::new(&[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], &[0.0; 3])?;
    /// let near = [-1e-13, 0.5];
    /// assert_eq!(cubic.eval(&near), None);
    /// let wider = cubic.with_simplex_tolerance(10.0)?;
    /// assert_eq!(wider.simplex_tolerance(), 10.0);
    /// assert!(wider.eval(&near).is_some());
    /// # Ok(())
    /// # }
    /// ```
    pub fn with_simplex_tolerance(mut self, tolerance: f64) -> Result<Self, InterpolateError> {
        self.mesh.set_tolerance(tolerance)?;
        Ok(self)
    }

    /// Returns the simplex tolerance. See
    /// [`with_simplex_tolerance`](CloughTocher2DInterpolator::with_simplex_tolerance)
    /// for an example.
    pub fn simplex_tolerance(&self) -> f64 {
        self.mesh.tolerance
    }
}

/// The product of two values, column by column.
fn times<V: Value<f64>>(a: V, b: V) -> V {
    a.zip_with(b, |x, y| x * y)
}

/// A value divided by a number, as the number itself would be divided (not
/// multiplied by its reciprocal).
fn over<V: Value<f64>>(a: V, d: f64) -> V {
    a.zip_with(a, |x, _| x / d)
}

/// The derivative of a gradient along `e`.
fn along<V: Value<f64>>(g: &[V; 2], e: [f64; 2]) -> V {
    V::zero().axpy(e[0], g[0]).axpy(e[1], g[1])
}

/// The gradients of the global estimate (Nielson 1983, Renka and Cline 1984):
/// those that minimize the sum over the edges of the integral of the square
/// of the second derivative of the patch along the edge. Setting the
/// derivative of that sum to zero gives, at each vertex `i`, with `E = p_j -
/// p_i` and `L = |E|` over its edges,
///
/// `sum 4 E E^T / L^3 G_i + sum 2 E E^T / L^3 G_j = -sum 6 (f_i - f_j) E / L^3`,
///
/// which is SciPy's update of one vertex (`_estimate_gradients_2d_global`)
/// as a whole: a symmetric positive definite system, solved here by conjugate
/// gradients with the 2 by 2 block of each vertex as the preconditioner,
/// every column of the values at once.
fn estimate_gradients<V: Value<f64>>(
    mesh: &Mesh<2>,
    values: &[V],
) -> Result<Vec<[V; 2]>, InterpolateError> {
    let n = mesh.points.len();
    // The edges, each once: from the triangle with the smaller index of the
    // two beside it, or from its only triangle on the hull.
    let mut edges: Vec<(u32, u32)> = Vec::with_capacity(mesh.simplices.len());
    for (s, t) in mesh.simplices.as_chunks::<3>().0.iter().enumerate() {
        for k in 0..3 {
            let beyond = mesh.neighbors[s * 3 + k];
            if beyond == NO_NEIGHBOR || beyond as usize > s {
                edges.push((t[(k + 1) % 3], t[(k + 2) % 3]));
            }
        }
    }
    // Per edge, twice its weight `E E^T / L^3` (xx, xy, yy): the block it
    // adds between its two ends. Per vertex, the sum of four times the
    // weights of its edges, and the right-hand side.
    let mut couple = Vec::with_capacity(edges.len());
    let mut block = vec![[0.0_f64; 3]; n];
    let mut active = vec![false; n];
    let mut rhs = vec![[V::zero(); 2]; n];
    for &(a, b) in &edges {
        let (pa, pb) = (mesh.points[a as usize], mesh.points[b as usize]);
        let e = [pb[0] - pa[0], pb[1] - pa[1]];
        let l = (e[0] * e[0] + e[1] * e[1]).sqrt();
        let l3 = l * l * l;
        let w = [e[0] * e[0] / l3, e[0] * e[1] / l3, e[1] * e[1] / l3];
        let s = [e[0] / l3, e[1] / l3];
        couple.push((a, b, w.map(|x| 2.0 * x)));
        for (i, j, sign) in [(a, b, 1.0), (b, a, -1.0)] {
            let (i, j) = (i as usize, j as usize);
            active[i] = true;
            for k in 0..3 {
                block[i][k] += 4.0 * w[k];
            }
            // -6 (f_i - f_j) E / L^3, with E from i to j.
            let df = values[i].axpy(-1.0, values[j]);
            rhs[i][0] = rhs[i][0].axpy(-6.0 * sign * s[0], df);
            rhs[i][1] = rhs[i][1].axpy(-6.0 * sign * s[1], df);
        }
    }
    let product = |x: &[[V; 2]], out: &mut [[V; 2]]| {
        for i in 0..n {
            let q = block[i];
            out[i] = [
                V::zero().axpy(q[0], x[i][0]).axpy(q[1], x[i][1]),
                V::zero().axpy(q[1], x[i][0]).axpy(q[2], x[i][1]),
            ];
        }
        for &(a, b, w) in &couple {
            let (a, b) = (a as usize, b as usize);
            let (xa, xb) = (x[a], x[b]);
            out[a][0] = out[a][0].axpy(w[0], xb[0]).axpy(w[1], xb[1]);
            out[a][1] = out[a][1].axpy(w[1], xb[0]).axpy(w[2], xb[1]);
            out[b][0] = out[b][0].axpy(w[0], xa[0]).axpy(w[1], xa[1]);
            out[b][1] = out[b][1].axpy(w[1], xa[0]).axpy(w[2], xa[1]);
        }
    };
    // The inverse of each vertex's block: the preconditioner.
    let inverse: Vec<[f64; 3]> = block
        .iter()
        .zip(&active)
        .map(|(q, &on)| {
            let det = q[0] * q[2] - q[1] * q[1];
            if on {
                [q[2] / det, -q[1] / det, q[0] / det]
            } else {
                [0.0; 3]
            }
        })
        .collect();
    let precondition = |i: usize, r: &[V; 2]| {
        let m = inverse[i];
        [
            V::zero().axpy(m[0], r[0]).axpy(m[1], r[1]),
            V::zero().axpy(m[1], r[0]).axpy(m[2], r[1]),
        ]
    };
    // A sum of products in the order of the vertices, as each dot product
    // of the iteration is taken.
    let add = |sum: V, a: &[V; 2], b: &[V; 2]| {
        sum.axpy(1.0, times(a[0], b[0]))
            .axpy(1.0, times(a[1], b[1]))
    };
    // Conjugate gradients from zero.
    let mut x = vec![[V::zero(); 2]; n];
    let mut r = rhs.clone();
    let z: Vec<[V; 2]> = (0..n).map(|i| precondition(i, &r[i])).collect();
    let mut p = z.clone();
    let (mut rr, mut rz, mut bb) = (V::zero(), V::zero(), V::zero());
    for i in 0..n {
        rr = add(rr, &r[i], &r[i]);
        rz = add(rz, &r[i], &z[i]);
        bb = add(bb, &rhs[i], &rhs[i]);
    }
    let mut z = z;
    let mut ap = vec![[V::zero(); 2]; n];
    // One in a column that is still to be solved, zero in one that is: a
    // column stops changing once it is solved, so that each column comes
    // out as it would alone.
    let open = |rr: V| rr.zip_with(bb, |r, b| f64::from(u8::from(r > SOLVED * SOLVED * b)));
    let limit = 10 * n + 100;
    for _ in 0..limit {
        let still = open(rr);
        if still.largest() == 0.0 {
            return Ok(x);
        }
        product(&p, &mut ap);
        let pap = (0..n).fold(V::zero(), |sum, i| add(sum, &p[i], &ap[i]));
        let alpha = rz.zip_with(pap, |a, b| if b == 0.0 { 0.0 } else { a / b });
        let alpha = times(alpha, still);
        // The new solution, residual, and preconditioned residual, and the
        // two dot products of the residual, in one pass.
        let (mut rr_next, mut rz_next) = (V::zero(), V::zero());
        for i in 0..n {
            for k in 0..2 {
                x[i][k] = x[i][k].axpy(1.0, times(alpha, p[i][k]));
                r[i][k] = r[i][k].axpy(-1.0, times(alpha, ap[i][k]));
            }
            z[i] = precondition(i, &r[i]);
            rr_next = add(rr_next, &r[i], &r[i]);
            rz_next = add(rz_next, &r[i], &z[i]);
        }
        let beta = rz_next.zip_with(rz, |a, b| if b == 0.0 { 0.0 } else { a / b });
        for i in 0..n {
            for k in 0..2 {
                p[i][k] = z[i][k].axpy(1.0, times(beta, p[i][k]));
            }
        }
        rr = rr_next;
        rz = rz_next;
    }
    if open(rr).largest() == 0.0 {
        return Ok(x);
    }
    Err(InterpolateError::SolverFailed(SolverError::NotConverged))
}

/// SciPy's `g[k]` of each triangle: the centroid of the triangle across the
/// edge opposite vertex `k`, in the barycentric coordinates of this one,
/// gives how the direction across that edge leans; `-1/2` on the hull.
fn leans(mesh: &Mesh<2>) -> Vec<[f64; 3]> {
    (0..mesh.simplex_count() as u32)
        .map(|s| {
            let mut g = [-0.5; 3];
            for (k, slot) in g.iter_mut().enumerate() {
                let t = mesh.neighbors[s as usize * 3 + k];
                if t == NO_NEIGHBOR {
                    continue;
                }
                let v = mesh.vertices(t);
                let centroid: [f64; 2] = core::array::from_fn(|d| {
                    (mesh.points[v[0] as usize][d]
                        + mesh.points[v[1] as usize][d]
                        + mesh.points[v[2] as usize][d])
                        / 3.0
                });
                let (c, last) = mesh.barycentric(s, &centroid);
                let c = [c[0], c[1], last];
                *slot = match k {
                    0 => (2.0 * c[2] + c[1] - 1.0) / (2.0 - 3.0 * c[2] - 3.0 * c[1]),
                    1 => (2.0 * c[0] + c[2] - 1.0) / (2.0 - 3.0 * c[0] - 3.0 * c[2]),
                    _ => (2.0 * c[1] + c[0] - 1.0) / (2.0 - 3.0 * c[1] - 3.0 * c[0]),
                };
            }
            g
        })
        .collect()
}

impl<V: Value<f64>> CloughTocher2DInterpolator<V> {
    /// The patch of triangle `simplex` at barycentric coordinates `b`:
    /// SciPy's `_clough_tocher_2d_single`, in its order of operations.
    fn patch(&self, simplex: u32, b: [f64; 3]) -> V {
        let v = self.mesh.vertices(simplex);
        let p = |k: usize| self.mesh.points[v[k] as usize];
        let (p1, p2, p3) = (p(0), p(1), p(2));
        let e12 = [p2[0] - p1[0], p2[1] - p1[1]];
        let e23 = [p3[0] - p2[0], p3[1] - p2[1]];
        let e31 = [p1[0] - p3[0], p1[1] - p3[1]];
        let (f1, f2, f3) = (
            self.values[v[0] as usize],
            self.values[v[1] as usize],
            self.values[v[2] as usize],
        );
        let (g1, g2, g3) = (
            &self.gradients[v[0] as usize],
            &self.gradients[v[1] as usize],
            &self.gradients[v[2] as usize],
        );
        let neg = |x: V| V::zero().axpy(-1.0, x);
        let df12 = along(g1, e12);
        let df21 = neg(along(g2, e12));
        let df23 = along(g2, e23);
        let df32 = neg(along(g3, e23));
        let df31 = along(g3, e31);
        let df13 = neg(along(g1, e31));

        // (a + 3 b) / 3, and the like, as SciPy writes them.
        let third = |a: V, b: V| over(a.axpy(3.0, b), 3.0);
        let mean3 = |a: V, b: V, c: V| over(a.axpy(1.0, b).axpy(1.0, c), 3.0);
        let c3000 = f1;
        let c2100 = third(df12, c3000);
        let c2010 = third(df13, c3000);
        let c0300 = f2;
        let c1200 = third(df21, c0300);
        let c0210 = third(df23, c0300);
        let c0030 = f3;
        let c1020 = third(df31, c0030);
        let c0120 = third(df32, c0030);
        let c2001 = mean3(c2100, c2010, c3000);
        let c0201 = mean3(c1200, c0300, c0210);
        let c0021 = mean3(c1020, c0120, c0030);

        let g = self.lean[simplex as usize];
        // (g * (-a + 3 b - 3 c + d) + (-a + 2 b - c + e + f)) / 2
        let side = |g: f64, a: V, b: V, c: V, d: V, e: V, f: V| {
            let cubic = neg(a).axpy(3.0, b).axpy(-3.0, c).axpy(1.0, d);
            let rest = neg(a).axpy(2.0, b).axpy(-1.0, c).axpy(1.0, e).axpy(1.0, f);
            over(V::zero().axpy(g, cubic).axpy(1.0, rest), 2.0)
        };
        let c0111 = side(g[0], c0300, c0210, c0120, c0030, c0021, c0201);
        let c1011 = side(g[1], c0030, c1020, c2010, c3000, c2001, c0021);
        let c1101 = side(g[2], c3000, c2100, c1200, c0300, c2001, c0201);
        let c1002 = mean3(c1101, c1011, c2001);
        let c0102 = mean3(c1101, c0111, c0201);
        let c0012 = mean3(c1011, c0111, c0021);
        let c0003 = mean3(c1002, c0102, c0012);

        // The extended barycentric coordinates: one of the four is zero.
        let least = b[0].min(b[1]).min(b[2]);
        let (b1, b2, b3, b4) = (b[0] - least, b[1] - least, b[2] - least, 3.0 * least);
        V::zero()
            .axpy(b1 * b1 * b1, c3000)
            .axpy(3.0 * b1 * b1 * b2, c2100)
            .axpy(3.0 * b1 * b1 * b3, c2010)
            .axpy(3.0 * b1 * b1 * b4, c2001)
            .axpy(3.0 * b1 * b2 * b2, c1200)
            .axpy(6.0 * b1 * b2 * b4, c1101)
            .axpy(3.0 * b1 * b3 * b3, c1020)
            .axpy(6.0 * b1 * b3 * b4, c1011)
            .axpy(3.0 * b1 * b4 * b4, c1002)
            .axpy(b2 * b2 * b2, c0300)
            .axpy(3.0 * b2 * b2 * b3, c0210)
            .axpy(3.0 * b2 * b2 * b4, c0201)
            .axpy(3.0 * b2 * b3 * b3, c0120)
            .axpy(6.0 * b2 * b3 * b4, c0111)
            .axpy(3.0 * b2 * b4 * b4, c0102)
            .axpy(b3 * b3 * b3, c0030)
            .axpy(3.0 * b3 * b3 * b4, c0021)
            .axpy(3.0 * b3 * b4 * b4, c0012)
            .axpy(b4 * b4 * b4, c0003)
    }
}

impl<V: Value<f64>> Interpolator for CloughTocher2DInterpolator<V> {
    type Point = [f64; 2];
    type Output = V;

    fn eval(&self, point: &[f64; 2]) -> Option<V> {
        let (simplex, c, last) = self.mesh.find(point)?;
        Some(self.patch(simplex, [c[0], c[1], last]))
    }
}

/// The stored form of the `serde` feature.
#[cfg(feature = "serde")]
mod stored {
    use super::{leans, CloughTocher2DInterpolator, Mesh};
    use crate::scattered::Scaling;
    use crate::{validate, InterpolateError, Value};
    use serde::{de, Deserialize, Deserializer, Serialize, Serializer};
    use std::borrow::Cow;

    /// The version of the stored form.
    const VERSION: u32 = 1;

    /// What is stored: what `LinearNDInterpolator` stores (the points, the
    /// values, the simplices and their neighbors, the change of coordinates,
    /// the simplex tolerance), and the gradients, so that the system of the
    /// estimate is not solved again when the form is read. The leans of the
    /// triangles are made again from the triangulation.
    #[derive(Serialize, Deserialize)]
    #[serde(bound = "")]
    struct Form<'a, V: Value<f64>> {
        version: u32,
        #[serde(with = "crate::stored")]
        points: Cow<'a, [[f64; 2]]>,
        #[serde(with = "crate::stored")]
        values: Cow<'a, [V]>,
        simplices: Cow<'a, [u32]>,
        neighbors: Cow<'a, [u32]>,
        scaling: Option<Scaling<2>>,
        simplex_tolerance: f64,
        #[serde(with = "crate::stored")]
        gradients: Cow<'a, [[V; 2]]>,
    }

    impl<V: Value<f64>> Serialize for CloughTocher2DInterpolator<V> {
        fn serialize<S: Serializer>(&self, to: S) -> Result<S::Ok, S::Error> {
            Form {
                version: VERSION,
                points: Cow::Borrowed(&self.mesh.points),
                values: Cow::Borrowed(&self.values),
                simplices: Cow::Borrowed(&self.mesh.simplices),
                neighbors: Cow::Borrowed(&self.mesh.neighbors),
                scaling: self.mesh.scaling,
                simplex_tolerance: self.mesh.tolerance,
                gradients: Cow::Borrowed(&self.gradients),
            }
            .serialize(to)
        }
    }

    /// The points, the values, and the triangulation pass the checks of the
    /// reading of `LinearNDInterpolator`; the gradients are one a point, and
    /// finite.
    impl<'de, V: Value<f64>> Deserialize<'de> for CloughTocher2DInterpolator<V> {
        fn deserialize<D2: Deserializer<'de>>(from: D2) -> Result<Self, D2::Error> {
            let form = Form::<V>::deserialize(from)?;
            crate::stored::version(form.version, VERSION)?;
            let points = form.points.into_owned();
            let values = form.values.into_owned();
            Mesh::check(&points, &values).map_err(de::Error::custom)?;
            let gradients = form.gradients.into_owned();
            validate::length(&gradients, points.len()).map_err(de::Error::custom)?;
            let finite = |g: &[V; 2]| g[0].all_finite() && g[1].all_finite();
            if let Some(index) = validate::first_failing(&gradients, finite) {
                return Err(de::Error::custom(InterpolateError::NonFinite { index }));
            }
            let mesh = Mesh::from_stored(
                points,
                form.simplices.into_owned(),
                form.neighbors.into_owned(),
                form.scaling,
                form.simplex_tolerance,
            )
            .map_err(de::Error::custom)?;
            let lean = leans(&mesh);
            Ok(Self {
                mesh,
                values,
                gradients,
                lean,
            })
        }
    }
}
