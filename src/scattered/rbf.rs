//! The radial basis function interpolant: its kernels, its polynomial terms,
//! the system they make, the solve, and the evaluation.
//!
//! The interpolant through the points `y_i` is a sum of one kernel centered at
//! every point, plus a polynomial of low degree:
//!
//! ```text
//! f(x) = sum over i of w_i phi(epsilon |x - y_i|) + sum over j of c_j p_j(x^)
//! ```
//!
//! The weights `w` and the coefficients `c` solve one symmetric system
//! ([`system`], [`solve`]); the value at a point is the sum above, term by
//! term. The formulas, the signs, and the order of every array are those of
//! SciPy 1.18's `RBFInterpolator`.

use super::kdtree::{self, Tree};
use crate::error::{ParamReason, SolverError};
use crate::{reserve, validate, InterpolateError, Interpolator, Value};
use faer::dyn_stack::{MemBuffer, MemStack};
use faer::linalg::lu::partial_pivoting::factor::{lu_in_place, lu_in_place_scratch};
use faer::{MatMut, Par};
use std::cell::RefCell;

/// The kernel of an [`RBFInterpolator`]: a function of the distance `r` from a
/// data point.
///
/// The kernels and their signs are SciPy's. Four of them take a shape
/// parameter `epsilon`, which the distance is multiplied by: a larger
/// `epsilon` makes a narrower bump. It must be finite and above zero. The
/// other four do not depend on a scale (stretching all distances leaves the
/// interpolant as it is), and take none.
///
/// A kernel may need a polynomial beside it for its system to have one
/// solution: the degree is given with each. [`RbfDegree::Default`] takes that
/// degree, and a constant for the kernels that need no polynomial.
///
/// SciPy takes an `epsilon` for the scale-free kernels as well. Without
/// smoothing it changes nothing. With smoothing `s` it is the same as the
/// smoothing `s / epsilon^p` here, with `p` 1 for the linear kernel, 2 for the
/// thin plate spline, 3 for the cubic, and 5 for the quintic.
///
/// The enum is `#[non_exhaustive]`: later versions add kernels.
///
/// # Examples
///
/// ```
/// use interpx::{InterpolateError, ParamReason, RBFInterpolator, RbfKernel};
///
/// let points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]];
/// let values = [1.0, 2.0, 3.0, 4.0];
/// let narrow = RbfKernel::Gaussian { epsilon: 3.0 };
/// assert!(RBFInterpolator::new(&points, &values, narrow).is_ok());
///
/// let err = RBFInterpolator::new(&points, &values, RbfKernel::Gaussian { epsilon: 0.0 });
/// let not_positive = InterpolateError::InvalidParameter {
///     name: "epsilon",
///     reason: ParamReason::NotPositive,
/// };
/// assert_eq!(err.err(), Some(not_positive));
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum RbfKernel {
    /// `-r`. Needs a polynomial of degree 0.
    Linear,
    /// `r^2 log r`. Needs a polynomial of degree 1. SciPy's default kernel.
    ThinPlateSpline,
    /// `r^3`. Needs a polynomial of degree 1.
    Cubic,
    /// `-r^5`. Needs a polynomial of degree 2.
    Quintic,
    /// `-sqrt(r^2 + 1)`, with `r` the distance times `epsilon`. Needs a
    /// polynomial of degree 0.
    Multiquadric {
        /// What the distance is multiplied by.
        epsilon: f64,
    },
    /// `1 / sqrt(r^2 + 1)`, with `r` the distance times `epsilon`. Needs no
    /// polynomial.
    InverseMultiquadric {
        /// What the distance is multiplied by.
        epsilon: f64,
    },
    /// `1 / (r^2 + 1)`, with `r` the distance times `epsilon`. Needs no
    /// polynomial.
    InverseQuadratic {
        /// What the distance is multiplied by.
        epsilon: f64,
    },
    /// `exp(-r^2)`, with `r` the distance times `epsilon`. Needs no
    /// polynomial.
    Gaussian {
        /// What the distance is multiplied by.
        epsilon: f64,
    },
}

impl RbfKernel {
    /// What the distances are multiplied by: 1 for a kernel that does not
    /// depend on a scale.
    fn epsilon(self) -> f64 {
        match self {
            Self::Linear | Self::ThinPlateSpline | Self::Cubic | Self::Quintic => 1.0,
            Self::Multiquadric { epsilon }
            | Self::InverseMultiquadric { epsilon }
            | Self::InverseQuadratic { epsilon }
            | Self::Gaussian { epsilon } => epsilon,
        }
    }

    /// The smallest degree of the polynomial the kernel needs for its system
    /// to have one solution, where it needs one: a kernel that is
    /// conditionally positive definite of order `m` needs the degree `m - 1`
    /// (Fasshauer, Meshfree Approximation Methods with MATLAB, chapter 8).
    fn min_degree(self) -> Option<usize> {
        match self {
            Self::Linear | Self::Multiquadric { .. } => Some(0),
            Self::ThinPlateSpline | Self::Cubic => Some(1),
            Self::Quintic => Some(2),
            Self::InverseMultiquadric { .. }
            | Self::InverseQuadratic { .. }
            | Self::Gaussian { .. } => None,
        }
    }

    /// The degree used when none is asked for: the smallest one, and a
    /// constant term for the kernels that need no polynomial (SciPy's default).
    fn default_degree(self) -> usize {
        self.min_degree().unwrap_or(0)
    }

    /// The value at the scaled distance `r`.
    fn at(self, r: f64) -> f64 {
        match self {
            Self::Linear => -r,
            Self::ThinPlateSpline if r > 0.0 => r * r * r.ln(),
            // The limit at 0, where the logarithm is not finite.
            Self::ThinPlateSpline => 0.0,
            Self::Cubic => r * r * r,
            Self::Quintic => -(r * r * r * r * r),
            Self::Multiquadric { .. } => -(r * r + 1.0).sqrt(),
            Self::InverseMultiquadric { .. } => (r * r + 1.0).sqrt().recip(),
            Self::InverseQuadratic { .. } => (r * r + 1.0).recip(),
            Self::Gaussian { .. } => (-(r * r)).exp(),
        }
    }

    /// The value for two points: the coordinates are scaled first and their
    /// distance taken then, as SciPy does.
    fn between<const D: usize>(self, a: &[f64; D], b: &[f64; D]) -> f64 {
        let epsilon = self.epsilon();
        let squared: f64 = a
            .iter()
            .zip(b)
            .map(|(a, b)| {
                let apart = a * epsilon - b * epsilon;
                apart * apart
            })
            .sum();
        self.at(squared.sqrt())
    }
}

/// The number of monomials of degree `degree` or less in `variables`
/// variables, `C(variables + degree, variables)`, when it fits.
fn monomial_count(variables: usize, degree: usize) -> Option<usize> {
    // C(degree + i, i) from C(degree + i - 1, i - 1): times `degree + i`,
    // over `i`. The division is exact at every step.
    let mut count: u128 = 1;
    for i in 1..=variables as u128 {
        count = count.checked_mul(degree as u128 + i)? / i;
    }
    usize::try_from(count).ok()
}

/// The number of monomials of the polynomial terms for `kernel`: of `degree`
/// or less, or none.
///
/// # Errors
///
/// Returns [`InterpolateError::InvalidParameter`] for `degree`: with
/// [`ParamReason::TooSmall`] when it is under the smallest degree the kernel
/// needs (no terms at all is accepted for every kernel, as in SciPy), and with
/// [`ParamReason::TooLarge`] when the monomials are too many to count.
fn count_of<const D: usize>(
    kernel: RbfKernel,
    degree: Option<usize>,
) -> Result<usize, InterpolateError> {
    let invalid = |reason| InterpolateError::InvalidParameter {
        name: "degree",
        reason,
    };
    if let (Some(degree), Some(min)) = (degree, kernel.min_degree()) {
        if degree < min {
            return Err(invalid(ParamReason::TooSmall { min }));
        }
    }
    // A power is an `i32` (what `powi` takes), and the count a `usize`.
    let fits = |degree: usize| i32::try_from(degree).is_ok() && monomial_count(D, degree).is_some();
    match degree {
        None => Ok(0),
        Some(degree) if fits(degree) => Ok(monomial_count(D, degree).unwrap_or(usize::MAX)),
        Some(degree) => {
            // The largest degree that fits is between 0, which does, and
            // this one, which does not.
            let (mut fitting, mut not) = (0, degree);
            while not - fitting > 1 {
                let middle = fitting + (not - fitting) / 2;
                if fits(middle) {
                    fitting = middle;
                } else {
                    not = middle;
                }
            }
            Err(invalid(ParamReason::TooLarge { max: fitting }))
        }
    }
}

/// The power of each variable in each monomial of the polynomial terms, in
/// SciPy's order: by degree, and among those of one degree by the variables
/// they multiply, taken in increasing order (`1, x, y, x^2, x y, y^2`).
///
/// The terms are for `kernel` on `points` points: of `degree` or less, or
/// none.
///
/// # Errors
///
/// Those of [`count_of`], and then [`InterpolateError::InsufficientPoints`]
/// when there are fewer points than monomials: the terms could not be told
/// apart on them.
fn monomials<const D: usize>(
    points: usize,
    kernel: RbfKernel,
    degree: Option<usize>,
) -> Result<Vec<[i32; D]>, InterpolateError> {
    let count = count_of::<D>(kernel, degree)?;
    if points < count {
        return Err(InterpolateError::InsufficientPoints {
            required: count,
            actual: points,
        });
    }

    let mut powers = Vec::with_capacity(count);
    // Without variables the constant is the only monomial.
    let degrees = degree.map(|degree| 0..=if D == 0 { 0 } else { degree });
    for total in degrees.into_iter().flatten() {
        // The variables a monomial multiplies, with repeats, in increasing
        // order: `[0, 1, 1]` is `x y^2`. All such lists of `total` variables,
        // in the order of a dictionary.
        let mut variables = vec![0; total];
        loop {
            let mut power = [0; D];
            for variable in &variables {
                power[*variable] += 1;
            }
            powers.push(power);
            // The next list: the last entry that can still grow grows, and
            // the entries after it start again from its value.
            let Some(at) = variables.iter().rposition(|v| v + 1 < D) else {
                break;
            };
            let next = variables[at] + 1;
            variables[at..].fill(next);
        }
    }
    Ok(powers)
}

/// The polynomial terms on a set of points: the monomials, and the change of
/// coordinates they are evaluated in. The coordinates of the points go to
/// `[-1, 1]` on each axis before a monomial is evaluated, so that the powers
/// stay of a size.
///
/// The monomials are those of the interpolant. The change of coordinates is
/// that of the points the terms are used on: all of them, or the neighbors of
/// a query.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Terms<'a, const D: usize> {
    powers: &'a [[i32; D]],
    /// The middle of the points on each axis.
    shift: [f64; D],
    /// Half the range of the points on each axis; 1 where they all agree.
    scale: [f64; D],
}

impl<'a, const D: usize> Terms<'a, D> {
    /// The terms with the monomials `powers`, on `points`.
    fn on(powers: &'a [[i32; D]], points: impl Iterator<Item = [f64; D]>) -> Self {
        let (mut least, mut most) = ([f64::INFINITY; D], [f64::NEG_INFINITY; D]);
        let mut any = false;
        for point in points {
            any = true;
            for d in 0..D {
                least[d] = least[d].min(point[d]);
                most[d] = most[d].max(point[d]);
            }
        }
        let mut shift = [0.0; D];
        let mut scale = [1.0; D];
        for d in 0..D {
            if any {
                shift[d] = (most[d] + least[d]) / 2.0;
                let half = (most[d] - least[d]) / 2.0;
                if half > 0.0 {
                    scale[d] = half;
                }
            }
        }
        Self {
            powers,
            shift,
            scale,
        }
    }

    /// The number of monomials.
    fn len(&self) -> usize {
        self.powers.len()
    }

    /// The value of every monomial at `point`, in their order.
    fn at(&self, point: &[f64; D]) -> impl Iterator<Item = f64> + 'a {
        let moved: [f64; D] = core::array::from_fn(|d| (point[d] - self.shift[d]) / self.scale[d]);
        self.powers
            .iter()
            .map(move |power| moved.iter().zip(power).map(|(x, p)| x.powi(*p)).product())
    }

    /// Writes the value of every monomial at `point` to `out`, which has one
    /// entry a monomial.
    fn row(&self, point: &[f64; D], out: &mut [f64]) {
        for (value, monomial) in out.iter_mut().zip(self.at(point)) {
            *value = monomial;
        }
    }
}

/// The side of the tiles the kernels of the matrix are worked out by: two
/// tiles of this side fit the first cache.
const TILE: usize = 64;

/// The system of equations of the interpolant through `values` at `points`:
/// the matrix, row by row, and the right-hand side.
///
/// With `n` points and `r` monomials the matrix has `n + r` rows and is
/// symmetric:
///
/// ```text
/// | A    P |      A[i][j] = kernel(points[i], points[j]), plus `smoothing`
/// | P^T  0 |                on the diagonal
///                 P[i][j] = monomial j at point i
/// ```
///
/// The right-hand side is the values, and a zero for every monomial.
///
/// # Errors
///
/// Returns [`InterpolateError::InvalidParameter`] for `points` with
/// [`ParamReason::TooLarge`] when the matrix would be larger than a `Vec`
/// holds, and [`InterpolateError::AllocationFailed`] when it would fit one and
/// there is not the memory for it.
fn system<const D: usize, V: Value<f64>>(
    points: &[[f64; D]],
    values: &[V],
    smoothing: f64,
    kernel: RbfKernel,
    terms: &Terms<D>,
) -> Result<(Vec<f64>, Vec<V>), InterpolateError> {
    let (n, r) = (points.len(), terms.len());
    let side = side_of(n, r).map_err(|max| InterpolateError::InvalidParameter {
        name: "points",
        reason: ParamReason::TooLarge { max },
    })?;
    // The square of the points: far more than the points themselves.
    let mut matrix = reserve::filled(side * side, 0.0)?;
    fill(&mut matrix, n, |i| points[i], smoothing, kernel, terms);
    let mut rhs = vec![V::zero(); side];
    rhs[..n].copy_from_slice(values);
    Ok((matrix, rhs))
}

/// The number of rows of the matrix of a [`system`] with `n` points and `r`
/// monomials.
///
/// # Errors
///
/// Returns the largest number of points whose matrix a `Vec` can hold (it
/// holds `isize::MAX` bytes at most, and asking it for more panics), when
/// this one cannot be held.
fn side_of(n: usize, r: usize) -> Result<usize, usize> {
    let most = isize::MAX.unsigned_abs() / size_of::<f64>();
    n.checked_add(r)
        .filter(|side| {
            side.checked_mul(*side)
                .is_some_and(|entries| entries <= most)
        })
        .ok_or_else(|| most.isqrt().saturating_sub(r))
}

/// Writes the matrix of a [`system`] on the `n` points that `point` gives to
/// `matrix`, which has the square of `n` and the number of monomials entries,
/// all zero.
fn fill<const D: usize>(
    matrix: &mut [f64],
    n: usize,
    point: impl Fn(usize) -> [f64; D],
    smoothing: f64,
    kernel: RbfKernel,
    terms: &Terms<D>,
) {
    let r = terms.len();
    let side = n + r;
    // The kernel of two points does not depend on their order: each is worked
    // out once and written at both places, tile by tile, so that both places
    // are in the cache. Against working out every entry, in ns an entry:
    //
    //   points   thin plate spline   Gaussian     linear, cubic, quintic
    //   100      3.4 for 6.3         3.2 for 6.1  1.0 to 1.5 for 1.6 to 2.1
    //   1000     5.0 for 7.0         4.9 for 6.6  2.6 to 2.9 for 2.2 to 2.3
    //   3000     5.3 for 7.2         5.0 for 6.8  2.7 to 3.1 for 2.2 to 2.3
    //
    // The kernels with a logarithm or an exponential gain at every size. The
    // cheap ones lose on a large matrix, where a second place to write costs
    // more than their arithmetic. A large matrix is where this step counts
    // least: it grows with the square of the points, and the solve that
    // follows with the cube. Filling the lower half in a pass of its own was
    // slower still for them (3.1 to 3.5 at 3000 points).
    for rows in (0..n).step_by(TILE) {
        for columns in (rows..n).step_by(TILE) {
            for i in rows..n.min(rows + TILE) {
                let from = point(i);
                for j in columns.max(i)..n.min(columns + TILE) {
                    let value = kernel.between(&from, &point(j));
                    matrix[i * side + j] = value;
                    matrix[j * side + i] = value;
                }
            }
        }
    }
    for i in 0..n {
        matrix[i * side + i] += smoothing;
        // The monomials at the point: a row of the upper right block, and a
        // column of the lower left one.
        let start = i * side + n;
        terms.row(&point(i), &mut matrix[start..start + r]);
        for j in 0..r {
            matrix[(n + j) * side + i] = matrix[start + j];
        }
    }
}

/// Solves the [`system`] of an interpolant: the weights of the kernels, then
/// the coefficients of the monomials.
///
/// # Errors
///
/// Returns [`InterpolateError::SolverFailed`] when the solution is not finite:
/// a pivot was zero (the system has no solution, or not one alone), or the
/// numbers did not fit.
fn solve<V: Value<f64>>(mut matrix: Vec<f64>, rhs: &[V]) -> Result<Vec<V>, InterpolateError> {
    let side = rhs.len();
    let (mut order, mut inverse) = (vec![0_usize; side], vec![0_usize; side]);
    factor(
        &mut matrix,
        &mut order,
        &mut inverse,
        &mut factor_scratch(side),
    );
    let mut solution: Vec<V> = order.iter().map(|row| rhs[*row]).collect();
    substitute(&matrix, &mut solution);
    if solution.iter().all(V::all_finite) {
        Ok(solution)
    } else {
        Err(InterpolateError::SolverFailed(SolverError::Singular))
    }
}

/// The working space [`factor`] needs for a matrix of `side` rows.
fn factor_scratch(side: usize) -> MemBuffer {
    MemBuffer::new(lu_in_place_scratch::<usize, f64>(
        side,
        side,
        Par::Seq,
        Default::default(),
    ))
}

/// Factors the matrix of a [`system`] in place, by `faer` (LU with partial
/// pivoting). The matrix is symmetric, so its rows are the columns `faer`
/// works on, and no second matrix is made.
///
/// `order` and `inverse` have one entry a row. Afterwards row `i` of the
/// factored matrix is row `order[i]` of the system. The factors are in
/// `matrix`, column after column: `U` on the diagonal and above it, and under
/// it `L`, whose diagonal of ones is not stored.
fn factor(matrix: &mut [f64], order: &mut [usize], inverse: &mut [usize], scratch: &mut MemBuffer) {
    let side = order.len();
    lu_in_place(
        MatMut::from_column_major_slice_mut(matrix, side, side),
        order,
        inverse,
        Par::Seq,
        MemStack::new(scratch),
        Default::default(),
    );
}

/// Solves with a matrix [`factor`] has factored. `solution` holds the
/// right-hand side in the order of the factored rows (`rhs[order[i]]` at `i`),
/// and is left with the solution.
///
/// The substitutions are on the values as they are: a value of several columns
/// is solved in one pass, by linear combinations of whole values. A pivot of
/// zero leaves entries that are infinite or NaN, which the caller looks for.
fn substitute<V: Value<f64>>(matrix: &[f64], solution: &mut [V]) {
    let side = solution.len();
    let column = |j: usize| &matrix[j * side..(j + 1) * side];
    // `L y = P b`: once `y[j]` is final it leaves every entry below it.
    for j in 0..side {
        let (lower, pivot) = (column(j), solution[j]);
        for i in j + 1..side {
            solution[i] = solution[i].axpy(-lower[i], pivot);
        }
    }
    // `U x = y`, from the last column back.
    for j in (0..side).rev() {
        let upper = column(j);
        solution[j] = V::zero().axpy(upper[j].recip(), solution[j]);
        let solved = solution[j];
        for i in 0..j {
            solution[i] = solution[i].axpy(-upper[i], solved);
        }
    }
}

/// What an evaluation in the `neighbors` mode works in: the system of the
/// neighbors, the row of the query, and what the search and the factorization
/// need.
///
/// One a thread ([`WORKING`]). The vectors grow to what the largest evaluation
/// on the thread asked for and stay: after the first evaluation, one of the
/// same size allocates nothing.
#[derive(Default)]
struct Working {
    /// The neighbors of the query.
    found: Vec<kdtree::Found>,
    /// The neighbors `matrix` holds the factored system of, in the order
    /// given. Empty when it holds none that may be used again: the queries of
    /// one batch share a factorization, and no others do.
    factored: Vec<u32>,
    matrix: Vec<f64>,
    /// The row of the query.
    row: Vec<f64>,
    /// The solution of the system with that row.
    solution: Vec<f64>,
    order: Vec<usize>,
    inverse: Vec<usize>,
    /// The working space of the factorization, and the largest side it has
    /// room for.
    factoring: Option<(usize, MemBuffer)>,
}

thread_local! {
    /// The working space of this thread.
    static WORKING: RefCell<Working> = RefCell::new(Working::default());
}

/// The degree of the polynomial an [`RBFInterpolator`] adds to its kernels.
///
/// # Examples
///
/// ```
/// use interpx::{RbfDegree, RbfOptions};
///
/// let options = RbfOptions::new().with_degree(RbfDegree::Degree(2));
/// assert_eq!(options.degree(), RbfDegree::Degree(2));
/// assert_eq!(RbfOptions::new().degree(), RbfDegree::Default);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum RbfDegree {
    /// The degree the kernel needs, and a constant term for the kernels that
    /// need none: SciPy's default.
    #[default]
    Default,
    /// No polynomial: the kernels alone (SciPy's `degree=-1`). With a kernel
    /// that needs a polynomial the system may have no solution.
    None,
    /// Every monomial of this degree or less. A degree under the one the
    /// kernel needs is refused.
    Degree(usize),
}

/// The options of an [`RBFInterpolator`] beside its kernel: the smoothing and
/// the degree of the polynomial.
///
/// # Examples
///
/// ```
/// use interpx::{RbfDegree, RbfOptions};
///
/// let options = RbfOptions::new().with_smoothing(0.5).with_degree(RbfDegree::None);
/// assert_eq!(options.smoothing(), 0.5);
/// assert_eq!(options.degree(), RbfDegree::None);
/// // What `RBFInterpolator::new` uses.
/// assert_eq!(RbfOptions::new().smoothing(), 0.0);
/// assert_eq!(RbfOptions::new(), RbfOptions::default());
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct RbfOptions {
    smoothing: f64,
    degree: RbfDegree,
    neighbors: Option<usize>,
}

impl RbfOptions {
    /// The defaults: no smoothing, the degree of [`RbfDegree::Default`], and
    /// every point for every query. See [`RbfOptions`] for an example.
    pub fn new() -> Self {
        Self::default()
    }

    /// Returns the options with another smoothing.
    ///
    /// With 0 the interpolant passes through the data. A larger value lets it
    /// miss the data to be smoother; as it grows the interpolant tends to the
    /// polynomial that fits the data best. It must be finite and not negative.
    /// See [`RbfOptions`] for an example.
    #[must_use]
    pub fn with_smoothing(mut self, smoothing: f64) -> Self {
        self.smoothing = smoothing;
        self
    }

    /// Returns the options with another degree of the polynomial. See
    /// [`RbfOptions`] for an example.
    #[must_use]
    pub fn with_degree(mut self, degree: RbfDegree) -> Self {
        self.degree = degree;
        self
    }

    /// Returns the options with a number of neighbors: every query is then
    /// answered from the `neighbors` data points nearest to it alone (SciPy's
    /// `neighbors`).
    ///
    /// For many points, where one system over all of them is too large. See
    /// [`RBFInterpolator`] for what changes. With as many neighbors as there
    /// are points, or more, all points are used as without this option.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::RbfOptions;
    ///
    /// assert_eq!(RbfOptions::new().neighbors(), None);
    /// assert_eq!(RbfOptions::new().with_neighbors(20).neighbors(), Some(20));
    /// ```
    #[must_use]
    pub fn with_neighbors(mut self, neighbors: usize) -> Self {
        self.neighbors = Some(neighbors);
        self
    }

    /// Returns the smoothing. See [`RbfOptions`] for an example.
    pub fn smoothing(&self) -> f64 {
        self.smoothing
    }

    /// Returns the number of neighbors, when one was given. See
    /// [`with_neighbors`](RbfOptions::with_neighbors) for an example.
    pub fn neighbors(&self) -> Option<usize> {
        self.neighbors
    }

    /// Returns the degree of the polynomial. See [`RbfOptions`] for an example.
    pub fn degree(&self) -> RbfDegree {
        self.degree
    }
}

/// The interpolant of values at scattered points by radial basis functions: a
/// sum of one kernel centered at every data point, and a polynomial of low
/// degree.
///
/// SciPy's `RBFInterpolator`. The points need no order and no structure, in
/// any number of dimensions. The weights of the kernels solve one system of
/// equations with a row for every point, at construction. Coordinates are
/// `f64`; a value is a scalar or a fixed-length array.
///
/// The interpolant passes through the data, unless
/// [smoothing](RbfOptions::with_smoothing) is asked for. It has a value
/// everywhere: there is no domain to be outside of. Far from the data the
/// value is what the kernel and the polynomial make of it, and is to be read
/// with care.
///
/// # Cost
///
/// Construction solves a dense system: with `n` points, a matrix of `n * n`
/// numbers (800 MB at 10 000 points, 80 GB at 100 000) and a time that grows
/// with the cube of `n`. A matrix there is no memory for is an error
/// ([`InterpolateError::AllocationFailed`]), not the end of the process. An
/// evaluation takes one kernel for every point. This suits up to a few
/// thousand points.
///
/// Two things are out of the reach of that error. `faer`, which factors the
/// matrix, allocates working space of its own (small beside the matrix). And
/// a system that promises more memory than it has (Linux, as it comes) lets
/// the allocation succeed and stops the process when the matrix is written.
///
/// # Many points: neighbors
///
/// With [`RbfOptions::with_neighbors`] a query is answered from the `k` data
/// points nearest to it alone. Construction then only sorts the points into a
/// tree. Every evaluation finds the neighbors and solves their system of `k`
/// rows: it costs far more than an evaluation of the other mode, and does not
/// depend on the number of points beyond the search. The result is the
/// interpolant of the neighbors, so it can jump where the set of neighbors
/// changes from one query to the next. A query whose neighbors have a system
/// without one solution has no value. Nor has any query when there is no
/// memory for a system of `k` rows (the square of `k` numbers: 80 GB at
/// 100 000 neighbors).
///
/// An evaluation in this mode works in a space its thread keeps. The space is
/// allocated at the first evaluation on the thread; evaluations after it do
/// not allocate.
///
/// [`eval_batch_into`](Interpolator::eval_batch_into) uses the system of a
/// point again for the next point when both have the same neighbors, with
/// the same results: put points that are close together next to each other
/// (the rows of an image, a path) and a batch is up to 1.6 times as fast as
/// a loop over `eval`.
///
/// # Points
///
/// Two points with the same coordinates are refused, with or without
/// smoothing. Points that all lie in a line (or, in three dimensions, in a
/// plane) leave the polynomial undetermined: when they lie along an axis the
/// system is found to have no solution, and otherwise the result can be
/// without meaning, as in SciPy.
///
/// # Examples
///
/// ```
/// use interpx::{InterpolateError, Interpolator, RBFInterpolator, RbfKernel};
///
/// # fn main() -> Result<(), InterpolateError> {
/// // x + 2 y at five points.
/// let points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0], [0.5, 0.25]];
/// let values = [0.0, 1.0, 2.0, 3.0, 1.0];
/// let rbf = RBFInterpolator::new(&points, &values, RbfKernel::ThinPlateSpline)?;
/// // The thin plate spline comes with a polynomial of degree 1, which takes
/// // up all of a plane: it is reproduced, between the points and away.
/// let close = |got: Option<f64>, want: f64| got.is_some_and(|v| (v - want).abs() < 1e-9);
/// assert!(close(rbf.eval(&[0.3, 0.6]), 1.5));
/// assert!(close(rbf.eval(&[4.0, -1.0]), 2.0));
/// assert_eq!(rbf.len(), 5);
/// assert_eq!(rbf.degree(), Some(1));
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone)]
pub struct RBFInterpolator<const D: usize, V: Value<f64> = f64> {
    points: Vec<[f64; D]>,
    kernel: RbfKernel,
    /// The monomials of the polynomial terms.
    powers: Vec<[i32; D]>,
    /// As asked for.
    smoothing: f64,
    degree: Option<usize>,
    solved: Solved<D, V>,
}

/// What an interpolant evaluates from.
#[derive(Debug, Clone)]
enum Solved<const D: usize, V> {
    /// One system over all points, solved at construction.
    Globally {
        /// The weight of the kernel at every point, then the coefficient of
        /// every monomial.
        coefficients: Vec<V>,
        /// The middle and half the range of the points on each axis: the
        /// change of coordinates of the polynomial terms.
        shift: [f64; D],
        scale: [f64; D],
    },
    /// A system over the neighbors of each query, solved at evaluation.
    Locally {
        values: Vec<V>,
        tree: Tree<D>,
        /// How many neighbors: fewer than there are points.
        neighbors: usize,
    },
}

impl<const D: usize, V: Value<f64>> RBFInterpolator<D, V> {
    /// Builds the interpolant through the values at the points, without
    /// smoothing and with the polynomial the kernel takes by default.
    ///
    /// # Errors
    ///
    /// Those of [`with_options`](RBFInterpolator::with_options).
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{InterpolateError, Interpolator, RBFInterpolator, RbfKernel};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// // A Gaussian takes the scale of its bump; values of two columns.
    /// let points = [[0.0], [1.0], [2.5], [4.0]];
    /// let values = [[0.0, 1.0], [1.0, 0.0], [0.5, 0.5], [2.0, -1.0]];
    /// let kernel = RbfKernel::Gaussian { epsilon: 0.8 };
    /// let rbf = RBFInterpolator::new(&points, &values, kernel)?;
    /// // Through the data.
    /// let at = rbf.eval(&[2.5]).unwrap_or([f64::NAN; 2]);
    /// assert!((at[0] - 0.5).abs() < 1e-9 && (at[1] - 0.5).abs() < 1e-9);
    /// # Ok(())
    /// # }
    /// ```
    pub fn new(
        points: &[[f64; D]],
        values: &[V],
        kernel: RbfKernel,
    ) -> Result<Self, InterpolateError> {
        Self::with_options(points, values, kernel, &RbfOptions::new())
    }

    /// Builds the interpolant with a smoothing or a degree of its polynomial
    /// other than the defaults.
    ///
    /// # Errors
    ///
    /// The first of these that applies:
    ///
    /// - [`InterpolateError::LengthMismatch`] when `points` and `values` differ
    ///   in length.
    /// - [`InterpolateError::InsufficientPoints`] when there is no point.
    /// - With neighbors, [`InterpolateError::InvalidParameter`] for `points`
    ///   when there are more than `u32::MAX` of them.
    /// - [`InterpolateError::NonFinite`] when a coordinate, or else a value, is
    ///   NaN or infinite. Its `index` is that of the point.
    /// - [`InterpolateError::InvalidParameter`] for `epsilon` when the one of
    ///   the kernel is not finite ([`ParamReason::NonFinite`]) or not above
    ///   zero ([`ParamReason::NotPositive`]).
    /// - [`InterpolateError::InvalidParameter`] for `smoothing` when it is not
    ///   finite ([`ParamReason::NonFinite`]) or under zero
    ///   ([`ParamReason::Negative`]).
    /// - [`InterpolateError::InvalidParameter`] for `degree` when it is under
    ///   the degree the kernel needs ([`ParamReason::TooSmall`]), or has more
    ///   monomials than can be counted ([`ParamReason::TooLarge`]).
    /// - [`InterpolateError::InvalidParameter`] for `neighbors` when it is 0
    ///   ([`ParamReason::TooSmall`]).
    /// - [`InterpolateError::InsufficientPoints`] when there are fewer points,
    ///   or fewer neighbors, than monomials in the polynomial.
    /// - [`InterpolateError::DuplicatePoints`] when two points have the same
    ///   coordinates.
    ///
    /// And, when all points are used:
    ///
    /// - [`InterpolateError::InvalidParameter`] for `points` when the matrix of
    ///   the system would be larger than a `Vec` can hold (`isize::MAX` bytes).
    /// - [`InterpolateError::AllocationFailed`] when the matrix fits a `Vec`
    ///   and there is not the memory for it (see "Cost" on
    ///   [`RBFInterpolator`]).
    /// - [`InterpolateError::SolverFailed`] when the system has no solution, or
    ///   not one alone: the points do not tell the monomials apart (all in a
    ///   line along an axis, with a polynomial of degree 1), or the kernel was
    ///   left without the polynomial it needs.
    ///
    /// With neighbors no system is solved here. In place of the three above:
    ///
    /// - [`InterpolateError::InvalidParameter`] for `neighbors` when the system
    ///   of that many neighbors would be larger than a `Vec` can hold.
    ///
    /// An evaluation whose neighbors have a system without one solution, or
    /// that cannot get the memory for that system, returns `None`.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{
    ///     InterpolateError, Interpolator, ParamReason, RBFInterpolator, RbfDegree, RbfKernel,
    ///     RbfOptions,
    /// };
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0], [0.4, 0.7]];
    /// let values = [1.0, 2.0, 0.0, 5.0, 2.0];
    /// // Smoothed: the interpolant no longer passes through the data.
    /// let options = RbfOptions::new().with_smoothing(0.5);
    /// let smooth = RBFInterpolator::with_options(&points, &values, RbfKernel::Linear, &options)?;
    /// assert_eq!(smooth.smoothing(), 0.5);
    /// assert!(smooth.eval(&[1.0, 1.0]).is_some_and(|v| (v - 5.0).abs() > 0.1));
    ///
    /// // The cubic kernel needs a polynomial of degree 1 at least.
    /// let options = RbfOptions::new().with_degree(RbfDegree::Degree(0));
    /// let err = RBFInterpolator::with_options(&points, &values, RbfKernel::Cubic, &options);
    /// let too_small = InterpolateError::InvalidParameter {
    ///     name: "degree",
    ///     reason: ParamReason::TooSmall { min: 1 },
    /// };
    /// assert_eq!(err.err(), Some(too_small));
    /// # Ok(())
    /// # }
    /// ```
    pub fn with_options(
        points: &[[f64; D]],
        values: &[V],
        kernel: RbfKernel,
        options: &RbfOptions,
    ) -> Result<Self, InterpolateError> {
        validate::length(values, points.len())?;
        validate::enough(points, 1)?;
        // With neighbors the points go into a tree that numbers them with
        // `u32`. (Counted before the points are looked at one by one.)
        let in_a_tree = options.neighbors.is_some_and(|k| k < points.len());
        if in_a_tree && u32::try_from(points.len()).is_err() {
            return Err(InterpolateError::InvalidParameter {
                name: "points",
                reason: ParamReason::TooLarge {
                    max: u32::MAX as usize,
                },
            });
        }
        let not_finite = |point: &[f64; D]| !point.iter().all(|c| c.is_finite());
        if let Some(index) = points.iter().position(not_finite) {
            return Err(InterpolateError::NonFinite { index });
        }
        if let Some(index) = values.iter().position(|v| !v.all_finite()) {
            return Err(InterpolateError::NonFinite { index });
        }
        let invalid = |name, reason| InterpolateError::InvalidParameter { name, reason };
        let epsilon = kernel.epsilon();
        if !epsilon.is_finite() {
            return Err(invalid("epsilon", ParamReason::NonFinite));
        }
        if epsilon <= 0.0 {
            return Err(invalid("epsilon", ParamReason::NotPositive));
        }
        let smoothing = options.smoothing;
        if !smoothing.is_finite() {
            return Err(invalid("smoothing", ParamReason::NonFinite));
        }
        if smoothing < 0.0 {
            return Err(invalid("smoothing", ParamReason::Negative));
        }
        let degree = match options.degree {
            RbfDegree::Default => Some(kernel.default_degree()),
            RbfDegree::None => None,
            RbfDegree::Degree(degree) => Some(degree),
        };
        count_of::<D>(kernel, degree)?;
        if options.neighbors == Some(0) {
            return Err(invalid("neighbors", ParamReason::TooSmall { min: 1 }));
        }
        // With as many neighbors as points, every query has all the points.
        let neighbors = options.neighbors.filter(|k| *k < points.len());
        let powers = monomials(neighbors.unwrap_or(points.len()), kernel, degree)?;
        if let Some((first, second)) = super::first_repeat(points) {
            return Err(InterpolateError::DuplicatePoints { first, second });
        }
        // The system of the neighbors is made at every evaluation, in a
        // matrix that has to be one a `Vec` can hold.
        if let Some(neighbors) = neighbors {
            side_of(neighbors, powers.len())
                .map_err(|max| invalid("neighbors", ParamReason::TooLarge { max }))?;
        }
        let solved = match neighbors {
            Some(neighbors) => Solved::Locally {
                values: values.to_vec(),
                tree: Tree::new(points),
                neighbors,
            },
            None => {
                let terms = Terms::on(&powers, points.iter().copied());
                let (matrix, rhs) = system(points, values, smoothing, kernel, &terms)?;
                Solved::Globally {
                    coefficients: solve(matrix, &rhs)?,
                    shift: terms.shift,
                    scale: terms.scale,
                }
            }
        };
        Ok(Self {
            points: points.to_vec(),
            kernel,
            powers,
            smoothing,
            degree,
            solved,
        })
    }

    /// Returns the number of data points. See [`RBFInterpolator`] for an
    /// example.
    #[expect(
        clippy::len_without_is_empty,
        reason = "a constructor refuses an empty set of points"
    )]
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// Returns the kernel.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{InterpolateError, RBFInterpolator, RbfKernel};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let kernel = RbfKernel::Multiquadric { epsilon: 2.0 };
    /// let rbf = RBFInterpolator::new(&[[0.0], [1.0], [3.0]], &[1.0, 2.0, 0.0], kernel)?;
    /// assert_eq!(rbf.kernel(), kernel);
    /// assert_eq!(rbf.smoothing(), 0.0);
    /// // The polynomial that was used: for this kernel, a constant.
    /// assert_eq!(rbf.degree(), Some(0));
    /// # Ok(())
    /// # }
    /// ```
    pub fn kernel(&self) -> RbfKernel {
        self.kernel
    }

    /// Returns the degree of the polynomial that was used, or `None` when
    /// there is none. See [`kernel`](RBFInterpolator::kernel) for an example.
    pub fn degree(&self) -> Option<usize> {
        self.degree
    }

    /// Returns the smoothing. See [`kernel`](RBFInterpolator::kernel) for an
    /// example.
    pub fn smoothing(&self) -> f64 {
        self.smoothing
    }

    /// Returns how many neighbors a query is answered from, or `None` when it
    /// is answered from all points: also when as many neighbors as there are
    /// points, or more, were asked for.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{InterpolateError, Interpolator, RBFInterpolator, RbfKernel, RbfOptions};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// // 100 points on a line, and 3 x at each.
    /// let points: Vec<[f64; 1]> = (0..100).map(|k| [f64::from(k)]).collect();
    /// let values: Vec<f64> = points.iter().map(|p| 3.0 * p[0]).collect();
    /// let options = RbfOptions::new().with_neighbors(8);
    /// let rbf = RBFInterpolator::with_options(&points, &values, RbfKernel::Linear, &options)?;
    /// assert_eq!(rbf.neighbors(), Some(8));
    /// // From the 8 points around 41.5.
    /// assert!(rbf.eval(&[41.5]).is_some_and(|v| (v - 124.5).abs() < 1e-9));
    ///
    /// let all = RbfOptions::new().with_neighbors(100);
    /// let rbf = RBFInterpolator::with_options(&points, &values, RbfKernel::Linear, &all)?;
    /// assert_eq!(rbf.neighbors(), None);
    /// # Ok(())
    /// # }
    /// ```
    pub fn neighbors(&self) -> Option<usize> {
        match &self.solved {
            Solved::Globally { .. } => None,
            Solved::Locally { neighbors, .. } => Some(*neighbors),
        }
    }

    /// Hands the value at each of `points` to `put`, in their order. With
    /// neighbors, the points share one working space, and with it the
    /// factored system of the point before.
    fn each(&self, points: &[[f64; D]], mut put: impl FnMut(Option<V>)) {
        let Solved::Locally {
            values,
            tree,
            neighbors,
        } = &self.solved
        else {
            for point in points {
                put(self.eval(point));
            }
            return;
        };
        let mut all = |working: &mut Working| {
            // What is factored there may be of another interpolator.
            working.factored.clear();
            for point in points {
                let finite = point.iter().all(|c| c.is_finite());
                put(finite
                    .then(|| self.eval_near(point, values, tree, *neighbors, working))
                    .flatten());
            }
        };
        // The working space of this thread, as `eval` takes it.
        let done = WORKING
            .try_with(|working| working.try_borrow_mut().ok().map(|mut w| all(&mut w)))
            .ok()
            .flatten();
        if done.is_none() {
            all(&mut Working::default());
        }
    }

    /// The value at `point` from its neighbors: their system is put together
    /// and factored in `working`, and solved for the row of the point. The sum
    /// of the neighbors' values, each times its entry of that solution, is the
    /// value of their interpolant at the point (the matrix is symmetric: the
    /// row times the solution for the values is the solution for the row
    /// times the values). Solving for the row takes one right-hand side of
    /// numbers whatever a value is.
    ///
    /// When `working.factored` names the neighbors of this point, their
    /// system is there already, factored for the query before, and is used as
    /// it is: the same numbers, so the same value. The caller empties
    /// `factored` before the first query of a batch.
    fn eval_near(
        &self,
        point: &[f64; D],
        values: &[V],
        tree: &Tree<D>,
        neighbors: usize,
        working: &mut Working,
    ) -> Option<V> {
        let Working {
            found,
            factored,
            matrix,
            row,
            solution,
            order,
            inverse,
            factoring,
        } = working;
        tree.nearest(point, neighbors, found);
        // In the order given, as SciPy has them.
        found.sort_unstable_by_key(|(_, place)| *place);
        let near = |i: usize| self.points[found[i].1 as usize];
        let terms = Terms::on(&self.powers, (0..neighbors).map(near));
        let side = neighbors + terms.len();

        let places = || found.iter().map(|(_, place)| *place);
        if !places().eq(factored.iter().copied()) {
            // Nothing is kept until the new system is factored.
            factored.clear();
            // The square of the neighbors: with many of them, more than there
            // may be memory for. There is no value then.
            matrix.clear();
            reserve::room(matrix, side * side).ok()?;
            matrix.resize(side * side, 0.0);
            fill(matrix, neighbors, near, self.smoothing, self.kernel, &terms);
            order.resize(side, 0);
            inverse.resize(side, 0);
            // Room for a larger matrix is room for this one.
            if factoring.as_ref().is_none_or(|(of, _)| *of < side) {
                *factoring = Some((side, factor_scratch(side)));
            }
            let (_, scratch) = factoring.as_mut()?;
            factor(matrix, order, inverse, scratch);
            factored.extend(places());
        }
        row.clear();
        row.resize(side, 0.0);
        let (kernels, tail) = row.split_at_mut(neighbors);
        for (i, entry) in kernels.iter_mut().enumerate() {
            *entry = self.kernel.between(point, &near(i));
        }
        terms.row(point, tail);
        // The row in the order of the factored rows.
        solution.clear();
        solution.extend(order.iter().map(|at| row[*at]));
        substitute(matrix, solution);
        let mut value = V::zero();
        for (weight, (_, place)) in solution.iter().zip(found.iter()) {
            value = value.axpy(*weight, values[*place as usize]);
        }
        value.all_finite().then_some(value)
    }
}

impl<const D: usize, V: Value<f64>> Interpolator for RBFInterpolator<D, V> {
    type Point = [f64; D];
    type Output = V;

    /// The kernels and the monomials, each times its coefficient, added one
    /// by one: no row of them is built. `None` for a point with a coordinate
    /// that is NaN or infinite, and for a value that is not finite.
    fn eval(&self, point: &[f64; D]) -> Option<V> {
        if !point.iter().all(|c| c.is_finite()) {
            return None;
        }
        let (coefficients, shift, scale) = match &self.solved {
            Solved::Globally {
                coefficients,
                shift,
                scale,
            } => (coefficients, *shift, *scale),
            Solved::Locally {
                values,
                tree,
                neighbors,
            } => {
                // The working space of this thread. Should it not be there
                // (the thread is ending, and its space is gone), one is made
                // for this call.
                let local = |working: &mut Working| {
                    // What is factored there may be of another interpolator.
                    working.factored.clear();
                    self.eval_near(point, values, tree, *neighbors, working)
                };
                return WORKING
                    .try_with(|working| working.try_borrow_mut().ok().map(|mut w| local(&mut w)))
                    .ok()
                    .flatten()
                    .unwrap_or_else(|| local(&mut Working::default()));
            }
        };
        let terms = Terms {
            powers: &self.powers,
            shift,
            scale,
        };
        let (weights, coefficients) = coefficients.split_at(self.points.len());
        // One running sum. Four of them, so that a term would not wait for the
        // one before, changed nothing (0.95 to 1.00 of the time): a term costs
        // what its kernel costs, a square root and for some a logarithm or an
        // exponential.
        let mut value = V::zero();
        for (weight, center) in weights.iter().zip(&self.points) {
            value = value.axpy(self.kernel.between(point, center), *weight);
        }
        for (coefficient, monomial) in coefficients.iter().zip(terms.at(point)) {
            value = value.axpy(monomial, *coefficient);
        }
        value.all_finite().then_some(value)
    }

    /// Evaluates at every point (so does [`eval_batch`](Self::eval_batch), in
    /// the same way). With neighbors, a point that has the same
    /// neighbors as the point before it uses their system as it was factored
    /// for that point: points closer together than the data (the pixels of an
    /// image, a path) come in runs with the same neighbors. The results are
    /// bit for bit those of a loop over `eval`.
    fn eval_batch_into(
        &self,
        points: &[[f64; D]],
        out: &mut [V],
        fill: V,
    ) -> Result<(), InterpolateError> {
        if points.len() != out.len() {
            return Err(InterpolateError::LengthMismatch {
                expected: points.len(),
                actual: out.len(),
            });
        }
        let mut slots = out.iter_mut();
        self.each(points, |value| {
            if let Some(slot) = slots.next() {
                *slot = value.unwrap_or(fill);
            }
        });
        Ok(())
    }

    /// As [`eval_batch_into`](Interpolator::eval_batch_into), with `None` for
    /// a point without a value.
    fn eval_batch(&self, points: &[[f64; D]]) -> Vec<Option<V>> {
        let mut out = Vec::with_capacity(points.len());
        self.each(points, |value| out.push(value));
        out
    }
}

#[cfg(test)]
mod tests {
    use super::{monomial_count, monomials, side_of, system, RbfKernel, Terms};
    use crate::common::{assert_close, load, reals, Real, Tol, DIRECT_F64};
    use crate::{InterpolateError as E, ParamReason};
    use serde_json::Value as Json;

    /// The polynomial terms on all of a set of points, as the constructor
    /// makes them for an interpolant that uses all its points.
    #[derive(Debug, PartialEq)]
    struct Polynomial<const D: usize> {
        powers: Vec<[i32; D]>,
        shift: [f64; D],
        scale: [f64; D],
    }

    impl<const D: usize> Polynomial<D> {
        fn new(points: &[[f64; D]], kernel: RbfKernel, degree: Option<usize>) -> Result<Self, E> {
            let powers = monomials(points.len(), kernel, degree)?;
            let terms = Terms::on(&powers, points.iter().copied());
            let (shift, scale) = (terms.shift, terms.scale);
            Ok(Self {
                powers,
                shift,
                scale,
            })
        }

        fn terms(&self) -> Terms<'_, D> {
            Terms {
                powers: &self.powers,
                shift: self.shift,
                scale: self.scale,
            }
        }

        fn len(&self) -> usize {
            self.powers.len()
        }

        fn row(&self, point: &[f64; D], out: &mut [f64]) {
            self.terms().row(point, out);
        }
    }

    fn nums(v: &Json) -> Vec<f64> {
        reals(&serde_json::from_value::<Vec<Real>>(v.clone()).unwrap_or_default())
    }

    fn rows(v: &Json) -> Vec<Vec<f64>> {
        serde_json::from_value::<Vec<Vec<Real>>>(v.clone())
            .unwrap_or_default()
            .iter()
            .map(|r| reals(r))
            .collect()
    }

    fn points_of<const D: usize>(v: &Json) -> Vec<[f64; D]> {
        rows(v)
            .iter()
            .map(|r| r.as_slice().try_into().expect("D coordinates"))
            .collect()
    }

    /// The row of an evaluation: the kernel between the point and every data
    /// point, then every monomial. Its inner product with the solution of the
    /// system is the value of the interpolant.
    fn evaluation_row<const D: usize>(
        at: &[f64; D],
        points: &[[f64; D]],
        kernel: RbfKernel,
        polynomial: &Polynomial<D>,
        out: &mut [f64],
    ) {
        let (kernels, monomials) = out.split_at_mut(points.len());
        for (entry, point) in kernels.iter_mut().zip(points) {
            *entry = kernel.between(at, point);
        }
        polynomial.row(at, monomials);
    }

    fn bits(v: &[f64]) -> Vec<u64> {
        v.iter().map(|x| x.to_bits()).collect()
    }

    /// Every kernel, with `epsilon` where it takes one.
    fn kernels(epsilon: f64) -> [(&'static str, RbfKernel); 8] {
        [
            ("linear", RbfKernel::Linear),
            ("thin_plate_spline", RbfKernel::ThinPlateSpline),
            ("cubic", RbfKernel::Cubic),
            ("quintic", RbfKernel::Quintic),
            ("multiquadric", RbfKernel::Multiquadric { epsilon }),
            (
                "inverse_multiquadric",
                RbfKernel::InverseMultiquadric { epsilon },
            ),
            ("inverse_quadratic", RbfKernel::InverseQuadratic { epsilon }),
            ("gaussian", RbfKernel::Gaussian { epsilon }),
        ]
    }

    /// One golden case in `D` dimensions: every part against SciPy's.
    fn check<const D: usize>(case: &Json) {
        let name = case["name"].as_str().expect("name");
        let options = &case["options"];
        let epsilon = options["epsilon"].as_f64().expect("epsilon");
        let (_, kernel) = kernels(epsilon)
            .into_iter()
            .find(|(named, _)| options["kernel"] == *named)
            .expect("a kernel");
        // -1 is SciPy's way to write no polynomial terms.
        let degree = usize::try_from(options["degree"].as_i64().expect("degree")).ok();
        let smoothing = options["smoothing"].as_f64().expect("smoothing");
        let points: Vec<[f64; D]> = points_of(&case["points"]);
        let values = nums(&case["values"]);

        let polynomial = Polynomial::new(&points, kernel, degree).expect("valid terms");
        let powers: Vec<Vec<i32>> = serde_json::from_value(case["powers"].clone()).expect("powers");
        assert_eq!(polynomial.len(), powers.len(), "{name}");
        for (got, want) in polynomial.powers.iter().zip(&powers) {
            assert_eq!(got.as_slice(), want.as_slice(), "{name}: the monomials");
        }
        // The middle and half the range: the same two operations, the same bits.
        assert_eq!(
            bits(&polynomial.shift),
            bits(&nums(&case["shift"])),
            "{name}"
        );
        assert_eq!(
            bits(&polynomial.scale),
            bits(&nums(&case["scale"])),
            "{name}"
        );

        let (matrix, rhs) =
            system(&points, &values, smoothing, kernel, &polynomial.terms()).expect("a system");
        let side = points.len() + powers.len();
        let want: Vec<f64> = rows(&case["matrix"]).concat();
        assert_eq!(
            (matrix.len(), want.len()),
            (side * side, side * side),
            "{name}"
        );
        for (k, (got, want)) in matrix.iter().zip(&want).enumerate() {
            let what = format!("{name}: matrix[{}][{}]", k / side, k % side);
            assert_close(*got, *want, DIRECT_F64, &what);
        }
        assert_eq!(
            bits(&rhs),
            bits(&nums(&case["rhs"])),
            "{name}: right-hand side"
        );

        let queries: Vec<[f64; D]> = points_of(&case["queries"]);
        let mut row = vec![0.0; side];
        for (q, (query, want)) in queries.iter().zip(rows(&case["expected"])).enumerate() {
            evaluation_row(query, &points, kernel, &polynomial, &mut row);
            for (k, (got, want)) in row.iter().zip(&want).enumerate() {
                let what = format!("{name}: row {q}, entry {k}");
                assert_close(*got, *want, DIRECT_F64, &what);
            }
        }
    }

    #[test]
    fn the_parts_are_those_of_scipys_internal_functions() {
        let file = load::<Json>("rbf_system").expect("golden file rbf_system");
        let mut by_dimension = [0; 4];
        let mut seen = std::collections::BTreeSet::new();
        for case in &file.cases {
            let dimension = rows(&case["points"])[0].len();
            match dimension {
                1 => check::<1>(case),
                2 => check::<2>(case),
                3 => check::<3>(case),
                other => panic!("points in {other} dimensions"),
            }
            by_dimension[dimension] += 1;
            seen.insert(
                case["options"]["kernel"]
                    .as_str()
                    .expect("kernel")
                    .to_owned(),
            );
        }
        assert_eq!(by_dimension, [0, 1, 15, 1]);
        assert_eq!(seen.len(), 8, "every kernel");
    }

    #[test]
    fn each_kernel_has_scipys_formula_and_sign() {
        // SciPy's `_kernel_matrix` at the distance 0.7.
        let want = [
            -0.7,
            -0.174_770_722_529_978_88,
            0.343,
            -0.168_07,
            -1.220_655_561_573_370_3,
            0.819_231_920_519_040_5,
            0.671_140_939_597_315_5,
            0.612_626_394_184_416_1,
        ];
        let tight = Tol {
            rtol: 1e-15,
            atol: 0.0,
        };
        for ((name, kernel), want) in kernels(1.0).into_iter().zip(want) {
            assert_close(kernel.at(0.7), want, tight, name);
            assert_close(kernel.between(&[0.0, 0.3], &[0.7, 0.3]), want, tight, name);
        }
        // The distance is multiplied by `epsilon` where the kernel has one:
        // half the distance and twice the `epsilon` is the same value.
        for ((name, kernel), want) in kernels(2.0).into_iter().zip(want).skip(4) {
            assert_close(kernel.between(&[0.35], &[0.0]), want, tight, name);
        }
        for (name, kernel) in kernels(2.0).into_iter().take(4) {
            let (got, want) = (kernel.between(&[0.35], &[0.0]), kernel.at(0.35));
            assert_eq!(got.to_bits(), want.to_bits(), "{name} takes no epsilon");
        }
        // At the distance 0: the thin plate spline is 0, not the NaN of
        // `0 * log 0`.
        let at_zero = [0.0, 0.0, 0.0, 0.0, -1.0, 1.0, 1.0, 1.0];
        for ((name, kernel), want) in kernels(3.0).into_iter().zip(at_zero) {
            let got = kernel.between(&[1.5, -2.0], &[1.5, -2.0]);
            assert_eq!((got + 0.0).to_bits(), f64::to_bits(want + 0.0), "{name}");
        }
    }

    #[test]
    fn a_degree_is_checked_against_the_kernel_and_the_points() {
        let points = [
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
            [0.5, 0.2],
            [0.1, 0.9],
        ];
        let smallest = [
            Some(0),
            Some(1),
            Some(1),
            Some(2),
            Some(0),
            None,
            None,
            None,
        ];
        for ((name, kernel), min) in kernels(1.0).into_iter().zip(smallest) {
            assert_eq!(kernel.min_degree(), min, "{name}");
            // The default: the smallest, and a constant where none is needed.
            assert_eq!(kernel.default_degree(), min.unwrap_or(0), "{name}");
            let terms = |degree| Polynomial::new(&points, kernel, degree).map(|p| p.len());
            // No terms at all is taken by every kernel, as in SciPy.
            assert_eq!(terms(None), Ok(0), "{name}");
            for degree in 0..=2 {
                let want = match min {
                    Some(min) if degree < min => Err(E::InvalidParameter {
                        name: "degree",
                        reason: ParamReason::TooSmall { min },
                    }),
                    // 1, 3, and 6 monomials in two variables.
                    _ => Ok((degree + 1) * (degree + 2) / 2),
                };
                assert_eq!(terms(Some(degree)), want, "{name} degree {degree}");
            }
            // Degree 3 has 10 monomials, and there are 6 points.
            assert_eq!(
                terms(Some(3)),
                Err(E::InsufficientPoints {
                    required: 10,
                    actual: 6
                }),
                "{name}"
            );
        }
        // One point fewer than the monomials is too few, as many is enough.
        let gaussian = RbfKernel::Gaussian { epsilon: 1.0 };
        assert_eq!(
            Polynomial::new(&points[..5], gaussian, Some(2)).err(),
            Some(E::InsufficientPoints {
                required: 6,
                actual: 5
            })
        );
        assert!(Polynomial::new(&points[..3], gaussian, Some(1)).is_ok());
        // Without points and without terms there is nothing to shift by: the
        // coordinates are left as they are.
        let none = Polynomial::new(&points[..0], gaussian, None).expect("no terms");
        assert_eq!((none.len(), bits(&none.shift)), (0, bits(&[0.0, 0.0])));
        assert_eq!(bits(&none.scale), bits(&[1.0, 1.0]));
        // A degree under the smallest is reported before too few points.
        assert_eq!(
            Polynomial::new(&points[..0], RbfKernel::Quintic, Some(1)).err(),
            Some(E::InvalidParameter {
                name: "degree",
                reason: ParamReason::TooSmall { min: 2 },
            })
        );
    }

    /// `C(n, k)` by Pascal's triangle, or `None` when it does not fit.
    fn binomial(n: usize, k: usize) -> Option<u128> {
        let mut row = vec![1_u128];
        for _ in 0..n {
            let mut next = vec![1_u128];
            for pair in row.windows(2) {
                next.push(pair[0].checked_add(pair[1])?);
            }
            next.push(1);
            row = next;
        }
        row.get(k).copied()
    }

    fn monomials_in<const D: usize>() {
        // Enough points for every count below.
        let points = vec![[0.5; D]; 300];
        let gaussian = RbfKernel::Gaussian { epsilon: 1.0 };
        for degree in 0..=5_usize {
            let terms = Polynomial::new(&points, gaussian, Some(degree)).expect("valid terms");
            let count = binomial(D + degree, D).expect("small");
            assert_eq!(terms.len() as u128, count, "{D} variables, degree {degree}");
            assert_eq!(monomial_count(D, degree).map(|c| c as u128), Some(count));
            let totals: Vec<i32> = terms.powers.iter().map(|p| p.iter().sum()).collect();
            assert!(
                totals.windows(2).all(|pair| pair[0] <= pair[1]),
                "by degree"
            );
            assert!(totals.iter().all(|total| *total <= degree as i32));
            let distinct: std::collections::BTreeSet<_> = terms.powers.iter().collect();
            assert_eq!(distinct.len(), terms.len(), "each once");
            assert!(terms.powers.iter().flatten().all(|p| *p >= 0));
        }
    }

    #[test]
    fn the_monomials_are_all_of_the_degree_or_less_once_each_by_degree() {
        monomials_in::<0>();
        monomials_in::<1>();
        monomials_in::<2>();
        monomials_in::<3>();
        monomials_in::<4>();
        // Two variables, degree 2, in SciPy's order: 1, x, y, x^2, x y, y^2.
        let points = [[0.0, 0.0]; 6];
        let terms = Polynomial::new(&points, RbfKernel::Linear, Some(2)).expect("valid terms");
        assert_eq!(
            terms.powers,
            [[0, 0], [1, 0], [0, 1], [2, 0], [1, 1], [0, 2]]
        );
    }

    #[test]
    fn a_count_of_monomials_that_does_not_fit_names_the_largest_degree_that_does() {
        let too_large = |max| {
            Some(E::InvalidParameter {
                name: "degree",
                reason: ParamReason::TooLarge { max },
            })
        };
        let gaussian = RbfKernel::Gaussian { epsilon: 1.0 };
        // In one variable the count always fits; the power must fit `i32`.
        let line: [[f64; 1]; 0] = [];
        let largest_power = i32::MAX as usize;
        for degree in [largest_power + 1, largest_power + 2, usize::MAX] {
            assert_eq!(
                Polynomial::new(&line, gaussian, Some(degree)).err(),
                too_large(largest_power)
            );
        }
        // Under that, it is the points that are too few.
        assert_eq!(
            Polynomial::new(&line, gaussian, Some(largest_power)).err(),
            Some(E::InsufficientPoints {
                required: largest_power + 1,
                actual: 0
            })
        );
        // In 40 variables the count is what does not fit: the largest degree
        // whose count does, by Pascal's triangle.
        let fits = |degree| binomial(40 + degree, 40).is_some_and(|c| c <= usize::MAX as u128);
        let largest = (0..200)
            .take_while(|degree| fits(*degree))
            .last()
            .expect("degree 0 fits");
        assert!(largest > 5 && !fits(largest + 1));
        let space: [[f64; 40]; 0] = [];
        for degree in [largest + 1, largest + 2, 1000, usize::MAX] {
            assert_eq!(
                Polynomial::new(&space, gaussian, Some(degree)).err(),
                too_large(largest),
                "degree {degree}"
            );
        }
        assert!(matches!(
            Polynomial::new(&space, gaussian, Some(largest)).err(),
            Some(E::InsufficientPoints { actual: 0, .. })
        ));
        assert_eq!(monomial_count(40, largest + 1), None);
    }

    #[test]
    fn the_system_is_symmetric_and_a_row_at_a_data_point_is_its_row_of_the_matrix() {
        let points = [
            [0.0, 0.5, 1.0],
            [2.0, -1.0, 0.0],
            [1.0, 1.0, 1.0],
            [-1.0, 0.25, 3.0],
            [0.5, 0.5, -2.0],
        ];
        // Two values at every point.
        let values = [
            [1.0, -1.0],
            [2.0, -2.0],
            [3.0, -3.0],
            [4.0, -4.0],
            [5.0, -5.0],
        ];
        let smoothing = 0.375;
        let mut checked = 0;
        for (name, kernel) in kernels(0.6) {
            let degree = kernel.default_degree();
            let Ok(polynomial) = Polynomial::new(&points, kernel, Some(degree)) else {
                // The quintic needs degree 2: 10 monomials for 5 points.
                assert_eq!(kernel, RbfKernel::Quintic);
                continue;
            };
            let (n, side) = (points.len(), points.len() + polynomial.len());
            let (matrix, rhs) =
                system(&points, &values, smoothing, kernel, &polynomial.terms()).expect("a system");
            for i in 0..side {
                for j in 0..side {
                    let (here, there) = (matrix[i * side + j], matrix[j * side + i]);
                    assert_eq!(here.to_bits(), there.to_bits(), "{name}: [{i}][{j}]");
                    if i >= n && j >= n {
                        assert_eq!(here.to_bits(), 0.0_f64.to_bits(), "{name}");
                    }
                }
            }
            // The values, then a zero for every monomial.
            assert_eq!(rhs.len(), side, "{name}");
            for (got, want) in rhs.iter().zip(values.iter().chain(&[[0.0, 0.0]; 4])) {
                assert_eq!(bits(got), bits(want), "{name}");
            }
            // The row of an evaluation at a data point: its row of the
            // matrix, without the smoothing on the diagonal.
            let mut row = vec![0.0; side];
            for i in 0..n {
                evaluation_row(&points[i], &points, kernel, &polynomial, &mut row);
                for j in 0..side {
                    let on_diagonal = if i == j { smoothing } else { 0.0 };
                    let want = matrix[i * side + j] - on_diagonal;
                    let what = format!("{name}: point {i}, entry {j}");
                    assert_close(row[j], want, DIRECT_F64, &what);
                }
            }
            checked += 1;
        }
        assert_eq!(checked, 7);
    }

    #[test]
    fn a_matrix_of_several_tiles_has_every_kernel_at_both_of_its_places() {
        // More points than two tiles are wide, at spacings that differ.
        let points: Vec<[f64; 2]> = (0..150_u32)
            .map(|k| [f64::from(k * 37 % 150) * 0.11, f64::from(k * k % 97) * 0.07])
            .collect();
        let values: Vec<f64> = (0..150).map(f64::from).collect();
        let smoothing = 0.5;
        for (name, kernel) in kernels(0.3) {
            let terms = Polynomial::new(&points, kernel, Some(kernel.default_degree()))
                .expect("valid terms");
            let side = points.len() + terms.len();
            let (matrix, _) =
                system(&points, &values, smoothing, kernel, &terms.terms()).expect("a system");
            let mut monomials = vec![0.0; terms.len()];
            for (i, a) in points.iter().enumerate() {
                for (j, b) in points.iter().enumerate() {
                    let on_diagonal = if i == j { smoothing } else { 0.0 };
                    let want = kernel.between(a, b) + on_diagonal;
                    let got = matrix[i * side + j];
                    assert_eq!(got.to_bits(), want.to_bits(), "{name}: [{i}][{j}]");
                }
                terms.row(a, &mut monomials);
                for (j, want) in monomials.iter().enumerate() {
                    let (right, below) = (matrix[i * side + 150 + j], matrix[(150 + j) * side + i]);
                    assert_eq!(right.to_bits(), want.to_bits(), "{name}: P[{i}][{j}]");
                    assert_eq!(below.to_bits(), want.to_bits(), "{name}: P^T[{j}][{i}]");
                }
            }
        }
    }

    #[test]
    fn a_matrix_there_is_no_memory_for_is_an_error_and_not_the_end_of_the_process() {
        // Points without coordinates and values without columns take no
        // memory, so their number can be anything.
        #[cfg(target_pointer_width = "64")]
        {
            // 2^28 points and a constant: a matrix of (2^28 + 1)^2 numbers, which
            // a `Vec` could hold (under 2^60) and no machine has.
            let count = 1_usize << 28;
            let (points, values) = (vec![[0.0_f64; 0]; count], vec![[0.0_f64; 0]; count]);
            let terms = Polynomial {
                powers: monomials::<0>(count, RbfKernel::Linear, Some(0)).expect("valid terms"),
                shift: [],
                scale: [],
            };
            let side = count + 1;
            assert_eq!(side_of(count, 1), Ok(side));
            assert_eq!(
                system(&points, &values, 0.0, RbfKernel::Linear, &terms.terms()).err(),
                Some(E::AllocationFailed {
                    bytes: side * side * 8
                })
            );
        }
    }

    #[test]
    fn the_side_of_a_system_is_checked_against_what_a_vec_holds() {
        assert_eq!(side_of(5, 3), Ok(8));
        assert_eq!(side_of(0, 0), Ok(0));
        #[cfg(target_pointer_width = "64")]
        {
            // A `Vec` of `f64` holds under 2^60 numbers: a side of 2^30 is one
            // too many.
            let most = (1_usize << 30) - 1;
            assert_eq!(side_of(most, 0), Ok(most));
            assert_eq!(side_of(most - 3, 3), Ok(most));
            assert_eq!(side_of(most + 1, 0), Err(most));
            // The monomials take their rows: fewer points fit.
            assert_eq!(side_of(most - 2, 3), Err(most - 3));
            assert_eq!(side_of(usize::MAX, 10), Err(most - 10));
            assert_eq!(side_of(7, usize::MAX), Err(0));
        }
    }

    #[test]
    fn a_matrix_larger_than_a_vec_holds_is_refused_before_it_is_asked_for() {
        // Points without coordinates and values without columns take no
        // memory, so their number can be anything.
        #[cfg(target_pointer_width = "64")]
        {
            let count = 1_usize << 31;
            let (points, values) = (vec![[0.0_f64; 0]; count], vec![[0.0_f64; 0]; count]);
            // The largest matrix of `f64` a `Vec` holds has under 2^60
            // entries: a side of 2^30 is one too many. One row and one column
            // are for the constant term.
            let powers = monomials::<0>(count, RbfKernel::Linear, Some(0)).expect("valid terms");
            assert_eq!(powers.len(), 1);
            // Without coordinates there is nothing to shift or scale by (and
            // no need to go through the points for it).
            let terms = Polynomial {
                powers,
                shift: [],
                scale: [],
            };
            let most = (isize::MAX.unsigned_abs() / 8).isqrt();
            assert_eq!(most, (1 << 30) - 1);
            assert_eq!(
                system(&points, &values, 0.0, RbfKernel::Linear, &terms.terms()).err(),
                Some(E::InvalidParameter {
                    name: "points",
                    reason: ParamReason::TooLarge { max: most - 1 },
                })
            );
        }
    }
}
