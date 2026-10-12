//! Interpolation of values on a rectilinear grid in several dimensions.

use crate::bspline::basis_values;
use crate::error::AxisFault;
use crate::oned::{hermite, pchip, Collocation};
use crate::{search, Extrapolate, Hint, InterpolateError, Interpolator, Scalar, Value};

/// How a [`RegularGridInterpolator`] interpolates between the points of its
/// grid.
///
/// The method is chosen at construction and stays. The enum is
/// `#[non_exhaustive]`: later versions add methods.
///
/// SciPy's `slinear` is [`GridMethod::Linear`]: the two are one function
/// under two names.
///
/// # Examples
///
/// ```
/// use interpx::{GridMethod, InterpolateError, Interpolator, RegularGridInterpolator};
///
/// # fn main() -> Result<(), InterpolateError> {
/// // x + 10 y on the grid {0, 1} x {0, 1, 2}.
/// let (x, y) = ([0.0, 1.0], [0.0, 1.0, 2.0]);
/// let values = [0.0, 10.0, 20.0, 1.0, 11.0, 21.0];
/// let at = |method| {
///     RegularGridInterpolator::new([&x[..], &y[..]], &values, method)
///         .map(|grid| grid.eval(&[0.25, 1.75]))
/// };
/// assert_eq!(at(GridMethod::Linear)?, Some(17.75));
/// assert_eq!(at(GridMethod::Nearest)?, Some(20.0));
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum GridMethod {
    /// The value of the nearest grid point, axis by axis. Halfway between two
    /// points of an axis the lower one is taken, as in SciPy.
    Nearest,
    /// Linear in each coordinate inside a cell of the grid (multilinear): the
    /// weighted mean of the values at the corners of the cell. Continuous, with
    /// kinks across the faces of the cells.
    Linear,
    /// A cubic spline in each coordinate (a tensor-product spline of degree
    /// 3): two continuous derivatives. Needs four points an axis.
    ///
    /// The spline is the one SciPy's `cubic_legacy` evaluates. SciPy's
    /// `cubic` aims at the same spline with an iterative solver by default,
    /// and stops well short of it: its results are not reproduced.
    Cubic,
    /// A quintic spline in each coordinate (degree 5): four continuous
    /// derivatives. Needs six points an axis. As for
    /// [`Cubic`](GridMethod::Cubic), the spline of SciPy's `quintic_legacy`.
    Quintic,
    /// The shape-preserving cubic of
    /// [`PchipInterpolator`](crate::PchipInterpolator) in each coordinate:
    /// SciPy's `pchip`. One continuous derivative. Inside the grid no value
    /// goes past the values at the corners of its cell, and where the data is
    /// flat the interpolant is; the spline methods overshoot there. Needs two
    /// points an axis, where SciPy asks for four.
    ///
    /// The axes are taken in SciPy's order, the last one first: the lines
    /// along the last axis are interpolated at the last coordinate, then those
    /// results along the axis before it, and so on down to the first. The
    /// order matters, because the slopes of this cubic are not linear in the
    /// values: the same values with two axes exchanged give another
    /// interpolant.
    ///
    /// Only the slopes along the last axis are found, and checked, at
    /// construction. With values so far apart that a slope along an earlier
    /// axis overflows, an evaluation returns a value that is not finite.
    ///
    /// ```
    /// use interpx::{GridMethod, InterpolateError, Interpolator, RegularGridInterpolator};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// // A step along x from 0 to 1, the same for every y.
    /// let (x, y) = ([0.0, 1.0, 2.0, 3.0, 4.0, 5.0], [0.0, 1.0, 2.0, 3.0]);
    /// let step = |x: &f64| if *x < 2.5 { 0.0 } else { 1.0 };
    /// let values: Vec<f64> = x.iter().flat_map(|x| y.map(|_| step(x))).collect();
    /// let at = |method| {
    ///     RegularGridInterpolator::new([&x[..], &y[..]], &values, method)
    ///         .map(|grid| grid.eval(&[1.5, 0.5]))
    /// };
    /// // Flat where the data is flat.
    /// assert_eq!(at(GridMethod::Pchip)?, Some(0.0));
    /// // The cubic spline swings below the data before the step.
    /// assert!(at(GridMethod::Cubic)?.is_some_and(|v| v < -0.05));
    /// # Ok(())
    /// # }
    /// ```
    Pchip,
}

impl GridMethod {
    /// The degree of the spline of a spline method.
    fn degree(self) -> Option<usize> {
        match self {
            Self::Nearest | Self::Linear | Self::Pchip => None,
            Self::Cubic => Some(3),
            Self::Quintic => Some(5),
        }
    }
}

/// The largest order (degree + 1) of a spline method: that of the quintic.
const MAX_ORDER: usize = 6;

/// What a grid evaluates by.
#[derive(Debug, Clone)]
enum Rule<T, const D: usize, V> {
    Nearest,
    Linear,
    Spline(Tensor<T, D, V>),
    /// The slopes of the shape-preserving cubic along the last axis at every
    /// point of the grid, in the order of its values.
    Pchip(Vec<V>),
}

#[cfg(test)]
impl<T, const D: usize, V> Rule<T, D, V> {
    /// The spline of a spline method.
    fn spline(&self) -> Option<&Tensor<T, D, V>> {
        match self {
            Self::Spline(spline) => Some(spline),
            Self::Nearest | Self::Linear | Self::Pchip(_) => None,
        }
    }
}

/// A tensor-product B-spline on a grid: the same degree and the not-a-knot
/// knots of its data points on every axis.
#[derive(Debug, Clone)]
struct Tensor<T, const D: usize, V> {
    /// 3 or 5.
    degree: usize,
    /// The knots of each axis.
    knots: [Vec<T>; D],
    /// The coefficients, as many as the grid has points and in the order of
    /// its values.
    coefficients: Vec<V>,
}

/// The interpolant of values on a rectilinear grid in `D` dimensions: the
/// points are all combinations of one coordinate from each axis.
///
/// SciPy's `RegularGridInterpolator`. The axes need not be evenly spaced. Each
/// increases strictly or decreases strictly, and has two points at least. The
/// values are one flat slice in C order, the last axis varying fastest: with
/// axes of `n0` and `n1` points, the value at point `i0` of the first axis and
/// `i1` of the second is `values[i0 * n1 + i1]` (the layout of a NumPy array
/// and of `ndarray`). A value is a scalar or a fixed-length array.
///
/// An axis given in decreasing order is kept in increasing order, with the
/// values reordered to match: [`axes`](RegularGridInterpolator::axes) and
/// [`values`](RegularGridInterpolator::values) return what is kept.
///
/// Outside the grid the result follows the [`Extrapolate`] mode, applied to
/// each coordinate on its own. It starts as [`Extrapolate::Forbid`]: a point
/// with one coordinate outside has no value. [`Extrapolate::Extend`] continues
/// the outermost cell (for [`GridMethod::Nearest`], the outermost point);
/// [`Extrapolate::Clamp`] moves the coordinate to the nearest end of its axis;
/// [`Extrapolate::Periodic`] wraps it around.
///
/// The spline methods ([`GridMethod::Cubic`], [`GridMethod::Quintic`]) solve
/// for the coefficients of their spline once, at construction, and keep them
/// beside the values: twice the memory of the grid. An evaluation reads the
/// coefficients around the point only. [`GridMethod::Pchip`] keeps its slopes
/// along the last axis beside the values, twice the memory as well, and works
/// out the rest at each evaluation, from four points an axis at most.
///
/// Values that are NaN or infinite are refused at construction, where SciPy
/// takes them and returns NaN near them: fill the gaps of a grid before it is
/// interpolated.
///
/// # Examples
///
/// ```
/// use interpx::{
///     Extrapolate, GridMethod, InterpolateError, Interpolator, RegularGridInterpolator,
/// };
///
/// # fn main() -> Result<(), InterpolateError> {
/// // x * y on the grid {0, 1, 3} x {0, 2}.
/// let (x, y) = ([0.0, 1.0, 3.0], [0.0, 2.0]);
/// let values = [0.0, 0.0, 0.0, 2.0, 0.0, 6.0];
/// let grid = RegularGridInterpolator::new([&x[..], &y[..]], &values, GridMethod::Linear)?;
/// assert_eq!(grid.shape(), [3, 2]);
/// // Multilinear interpolation reproduces x * y.
/// assert_eq!(grid.eval(&[2.0, 1.0]), Some(2.0));
/// assert_eq!(grid.eval(&[0.5, 0.5]), Some(0.25));
/// // Outside the grid there is no value, unless the mode says otherwise.
/// assert_eq!(grid.eval(&[4.0, 1.0]), None);
/// assert_eq!(grid.eval_or(&[4.0, 1.0], -1.0), -1.0);
/// let extended = grid.with_extrapolate(Extrapolate::Extend);
/// assert_eq!(extended.eval(&[4.0, 1.0]), Some(4.0));
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone)]
pub struct RegularGridInterpolator<T: Scalar, const D: usize, V: Value<T> = T> {
    /// The coordinates of each axis, increasing.
    axes: [Vec<T>; D],
    /// The values in C order, for the axes as they are kept.
    values: Vec<V>,
    /// How far apart in `values` two neighbors along each axis are.
    strides: [usize; D],
    rule: Rule<T, D, V>,
    extrapolate: Extrapolate,
}

impl<T: Scalar, const D: usize, V: Value<T>> RegularGridInterpolator<T, D, V> {
    /// Builds the interpolant from copies of the axes and of the values.
    ///
    /// # Errors
    ///
    /// Those of [`from_owned`](RegularGridInterpolator::from_owned).
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{GridMethod, InterpolateError, Interpolator, RegularGridInterpolator};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// // Three dimensions, two points an axis: the corners of the unit cube,
    /// // with x + 2 y + 4 z at each.
    /// let side = [0.0, 1.0];
    /// let values = [0.0, 4.0, 2.0, 6.0, 1.0, 5.0, 3.0, 7.0];
    /// let cube =
    ///     RegularGridInterpolator::new([&side[..]; 3], &values, GridMethod::Linear)?;
    /// assert_eq!(cube.eval(&[0.5, 0.5, 0.5]), Some(3.5));
    ///
    /// // A value is a scalar or an array: two columns here.
    /// let pairs = [[0.0, 1.0], [2.0, 1.0], [4.0, 1.0], [6.0, 1.0]];
    /// let plane =
    ///     RegularGridInterpolator::new([&side[..]; 2], &pairs, GridMethod::Linear)?;
    /// assert_eq!(plane.eval(&[0.5, 0.5]), Some([3.0, 1.0]));
    /// # Ok(())
    /// # }
    /// ```
    pub fn new(
        axes: [&[T]; D],
        values: &[V],
        method: GridMethod,
    ) -> Result<Self, InterpolateError> {
        Self::from_owned(axes.map(<[T]>::to_vec), values.to_vec(), method)
    }

    /// Builds the interpolant from the axes and the values, which it takes
    /// over.
    ///
    /// # Errors
    ///
    /// The first of these that applies:
    ///
    /// - [`InterpolateError::LengthMismatch`] when `values` does not hold one
    ///   value for every point of the grid: the product of the lengths of the
    ///   axes.
    /// - [`InterpolateError::GridAxis`] with [`AxisFault::InsufficientPoints`]
    ///   for the first axis with fewer points than the method needs: four for
    ///   [`GridMethod::Cubic`], six for [`GridMethod::Quintic`], two for the
    ///   others.
    /// - [`InterpolateError::GridAxis`] with [`AxisFault::NonFinite`] for the
    ///   first axis with a coordinate that is NaN or infinite.
    /// - [`InterpolateError::NonFinite`] when a value is NaN or infinite. Its
    ///   `index` is the position in `values`.
    /// - [`InterpolateError::GridAxis`] with [`AxisFault::DuplicatePoints`] or
    ///   [`AxisFault::NotMonotonic`] for the first axis that neither increases
    ///   strictly nor decreases strictly.
    /// - For a spline method, [`InterpolateError::GridAxis`] with
    ///   [`AxisFault::NonFinite`] for an axis whose spline cannot be set up
    ///   (coordinates so far apart that their distance overflows), and
    ///   [`InterpolateError::NonFinite`] with a position in `values` when a
    ///   coefficient of the spline overflows.
    /// - For [`GridMethod::Pchip`], [`InterpolateError::NonFinite`] with a
    ///   position in `values` when a slope along the last axis overflows: a
    ///   point of the interval whose values are too far apart.
    ///
    /// Positions are those in the axes and in `values` as they were given, also
    /// where an axis descends and is kept in the other order.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{
    ///     AxisFault, GridMethod, InterpolateError, Interpolator, RegularGridInterpolator,
    /// };
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// // The first axis decreases: it is kept increasing, with its values.
    /// let axes = [vec![2.0, 1.0, 0.0], vec![0.0, 1.0]];
    /// let values = vec![20.0, 21.0, 10.0, 11.0, 0.0, 1.0];
    /// let grid = RegularGridInterpolator::from_owned(axes, values, GridMethod::Linear)?;
    /// assert_eq!(grid.axes()[0], [0.0, 1.0, 2.0]);
    /// assert_eq!(grid.values(), &[0.0, 1.0, 10.0, 11.0, 20.0, 21.0]);
    /// assert_eq!(grid.eval(&[1.5, 0.5]), Some(15.5));
    ///
    /// // An axis that turns around.
    /// let axes = [vec![0.0, 1.0], vec![0.0, 2.0, 1.0]];
    /// let err = RegularGridInterpolator::from_owned(axes, vec![0.0; 6], GridMethod::Nearest);
    /// let turned = InterpolateError::GridAxis {
    ///     axis: 1,
    ///     fault: AxisFault::NotMonotonic { index: 2 },
    /// };
    /// assert_eq!(err.err(), Some(turned));
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_owned(
        mut axes: [Vec<T>; D],
        mut values: Vec<V>,
        method: GridMethod,
    ) -> Result<Self, InterpolateError> {
        // A product that does not fit `usize` is the length of no slice.
        let points = axes
            .iter()
            .try_fold(1_usize, |product, axis| product.checked_mul(axis.len()));
        if points != Some(values.len()) {
            return Err(InterpolateError::LengthMismatch {
                expected: points.unwrap_or(usize::MAX),
                actual: values.len(),
            });
        }
        let on_axis =
            |axis: usize| move |fault: AxisFault| InterpolateError::GridAxis { axis, fault };
        // A spline of degree `k` takes `k + 1` points to fix one polynomial.
        let required = method.degree().map_or(2, |degree| degree + 1);
        for (d, axis) in axes.iter().enumerate() {
            if axis.len() < required {
                return Err(on_axis(d)(AxisFault::InsufficientPoints {
                    required,
                    actual: axis.len(),
                }));
            }
        }
        for (d, axis) in axes.iter().enumerate() {
            if let Some(index) = axis.iter().position(|c| !c.is_finite()) {
                return Err(on_axis(d)(AxisFault::NonFinite { index }));
            }
        }
        if let Some(index) = values.iter().position(|v| !v.all_finite()) {
            return Err(InterpolateError::NonFinite { index });
        }
        let mut descending = [false; D];
        for (d, axis) in axes.iter().enumerate() {
            descending[d] = decreases(axis).map_err(on_axis(d))?;
        }

        // Neighbors along the last axis are next to each other; along an
        // earlier axis they are a whole block of the later axes apart.
        let mut strides = [1; D];
        for d in (0..D.saturating_sub(1)).rev() {
            strides[d] = strides[d + 1] * axes[d + 1].len();
        }
        for d in 0..D {
            if descending[d] {
                axes[d].reverse();
                reverse_along(&mut values, axes[d].len(), strides[d]);
            }
        }
        let rule = match method {
            GridMethod::Nearest => Ok(Rule::Nearest),
            GridMethod::Linear => Ok(Rule::Linear),
            GridMethod::Cubic => tensor(&axes, &values, &strides, 3).map(Rule::Spline),
            GridMethod::Quintic => tensor(&axes, &values, &strides, 5).map(Rule::Spline),
            GridMethod::Pchip => slopes_along_last(&axes, &values).map(Rule::Pchip),
        }
        .map_err(|overflow| overflow.into_given(&axes, &strides, &descending))?;
        Ok(Self {
            axes,
            values,
            strides,
            rule,
            extrapolate: Extrapolate::Forbid,
        })
    }

    /// Returns the same interpolant with another behavior outside the grid.
    ///
    /// The mode applies to each coordinate on its own.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{
    ///     Extrapolate, GridMethod, InterpolateError, Interpolator, RegularGridInterpolator,
    /// };
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// // x + 10 y on the unit square.
    /// let side = [0.0, 1.0];
    /// let values = [0.0, 10.0, 1.0, 11.0];
    /// let grid = RegularGridInterpolator::new([&side[..]; 2], &values, GridMethod::Linear)?;
    /// assert_eq!(grid.extrapolate(), Extrapolate::Forbid);
    /// assert_eq!(grid.method(), GridMethod::Linear);
    /// assert_eq!(grid.domain(), ([0.0, 0.0], [1.0, 1.0]));
    ///
    /// // 2 is outside in x, 0.5 is inside in y.
    /// let outside = [2.0, 0.5];
    /// assert_eq!(grid.eval(&outside), None);
    /// let extend = grid.clone().with_extrapolate(Extrapolate::Extend);
    /// assert_eq!(extend.eval(&outside), Some(7.0));
    /// let clamp = grid.clone().with_extrapolate(Extrapolate::Clamp);
    /// assert_eq!(clamp.eval(&outside), Some(6.0));
    /// let periodic = grid.with_extrapolate(Extrapolate::Periodic);
    /// assert_eq!(periodic.eval(&[2.25, 0.5]), Some(5.25));
    /// # Ok(())
    /// # }
    /// ```
    #[must_use]
    pub fn with_extrapolate(mut self, extrapolate: Extrapolate) -> Self {
        self.extrapolate = extrapolate;
        self
    }

    /// Returns the behavior outside the grid. See
    /// [`with_extrapolate`](RegularGridInterpolator::with_extrapolate) for an
    /// example.
    pub fn extrapolate(&self) -> Extrapolate {
        self.extrapolate
    }

    /// Returns the method. See
    /// [`with_extrapolate`](RegularGridInterpolator::with_extrapolate) for an
    /// example.
    pub fn method(&self) -> GridMethod {
        match &self.rule {
            Rule::Nearest => GridMethod::Nearest,
            Rule::Linear => GridMethod::Linear,
            Rule::Spline(spline) if spline.degree == 3 => GridMethod::Cubic,
            Rule::Spline(_) => GridMethod::Quintic,
            Rule::Pchip(_) => GridMethod::Pchip,
        }
    }

    /// Returns the corner of the grid with the smallest coordinates and the one
    /// with the largest. See
    /// [`with_extrapolate`](RegularGridInterpolator::with_extrapolate) for an
    /// example.
    pub fn domain(&self) -> ([T; D], [T; D]) {
        (
            core::array::from_fn(|d| self.axes[d][0]),
            core::array::from_fn(|d| self.axes[d][self.axes[d].len() - 1]),
        )
    }

    /// Returns the coordinates of each axis, in increasing order. See
    /// [`from_owned`](RegularGridInterpolator::from_owned) for an example.
    pub fn axes(&self) -> &[Vec<T>; D] {
        &self.axes
    }

    /// Returns the number of points of each axis. See
    /// [`RegularGridInterpolator`] for an example.
    pub fn shape(&self) -> [usize; D] {
        core::array::from_fn(|d| self.axes[d].len())
    }

    /// Returns the values in C order, for the axes as
    /// [`axes`](RegularGridInterpolator::axes) returns them. See
    /// [`from_owned`](RegularGridInterpolator::from_owned) for an example.
    pub fn values(&self) -> &[V] {
        &self.values
    }

    /// Evaluates at `point`, looking on each axis first at the interval in its
    /// hint and the neighbors of that interval.
    ///
    /// The result is bit for bit that of [`eval`](Interpolator::eval). Use this
    /// for points that follow a path through the grid; for points in no order
    /// `eval` is faster, as [`PPoly::eval_hinted`](crate::PPoly::eval_hinted)
    /// explains.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{GridMethod, Hint, InterpolateError, Interpolator, RegularGridInterpolator};
    ///
    /// # fn main() -> Result<(), InterpolateError> {
    /// let (x, y) = ([0.0, 1.0, 2.0, 4.0], [0.0, 1.0, 3.0]);
    /// let values: Vec<f64> = (0..12).map(f64::from).collect();
    /// let grid = RegularGridInterpolator::new([&x[..], &y[..]], &values, GridMethod::Linear)?;
    /// let mut hints = [Hint::new(); 2];
    /// for point in [[0.1, 0.2], [0.4, 0.9], [1.5, 2.5], [3.9, 2.9]] {
    ///     assert_eq!(grid.eval_hinted(&point, &mut hints), grid.eval(&point));
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn eval_hinted(&self, point: &[T; D], hints: &mut [Hint; D]) -> Option<V> {
        let mut cell = [(0, T::zero()); D];
        for (d, place) in cell.iter_mut().enumerate() {
            let axis = &self.axes[d];
            *place =
                search::eval_hinted(axis, self.extrapolate, point[d], &mut hints[d], |i, x| {
                    (i, x)
                })?;
        }
        Some(self.at(&cell))
    }

    /// The value for a point given by the interval of each of its coordinates
    /// (the index of the grid point below it) and the coordinate to evaluate
    /// at.
    fn at(&self, cell: &[(usize, T); D]) -> V {
        // How far each coordinate is from the lower point of its interval
        // toward the upper one.
        let fractions = || {
            let along = cell.iter().zip(&self.axes);
            along.map(|((lower, x), axis)| (*lower, fraction(axis, *lower, *x)))
        };
        match &self.rule {
            Rule::Nearest => {
                let half = (T::one() + T::one()).recip();
                let offset = fractions()
                    .zip(&self.strides)
                    .map(|((lower, fraction), stride)| {
                        (lower + usize::from(fraction > half)) * stride
                    })
                    .sum::<usize>();
                self.values[offset]
            }
            Rule::Linear => {
                // For each axis, the weight and the place in `values` of the
                // lower and of the upper side of the cell.
                let mut weights = [[T::zero(); 2]; D];
                let mut offsets = [[0; 2]; D];
                for (d, ((lower, fraction), stride)) in fractions().zip(&self.strides).enumerate() {
                    weights[d] = [T::one() - fraction, fraction];
                    offsets[d] = [lower * stride, (lower + 1) * stride];
                }
                // The corners of the cell, the first axis outermost (the order
                // SciPy adds them in): bit `d` from the top says whether the
                // corner is on the upper side along axis `d`.
                let mut value = V::zero();
                for corner in 0..1_usize << D {
                    let mut weight = T::one();
                    let mut offset = 0;
                    for d in 0..D {
                        let side = corner >> (D - 1 - d) & 1;
                        weight = weight * weights[d][side];
                        offset += offsets[d][side];
                    }
                    value = value.axpy(weight, self.values[offset]);
                }
                value
            }
            Rule::Spline(spline) => self.spline_at(spline, cell),
            Rule::Pchip(slopes) => self.pchip_at(slopes, cell),
        }
    }

    /// The value of the tensor-product spline, by the degree: the loops over
    /// the splines of an axis then have a length the compiler knows.
    fn spline_at(&self, spline: &Tensor<T, D, V>, cell: &[(usize, T); D]) -> V {
        if spline.degree == 3 {
            self.spline_sum::<4>(spline, cell)
        } else {
            self.spline_sum::<MAX_ORDER>(spline, cell)
        }
    }

    /// The value of a tensor-product spline of degree `N - 1`: the sum over
    /// the `N^D` splines that are not zero at the point, each the product of
    /// one spline of every axis, times its coefficient.
    fn spline_sum<const N: usize>(&self, spline: &Tensor<T, D, V>, cell: &[(usize, T); D]) -> V {
        let Some(last) = D.checked_sub(1) else {
            // No axes: the one coefficient, which is the one value.
            return spline.coefficients[0];
        };
        let degree = N - 1;
        // For each axis, the values of its `N` splines at the coordinate, and
        // the place in the coefficients of the first of them.
        let mut basis = [[T::zero(); N]; D];
        let mut first = [0; D];
        for d in 0..D {
            let (interval, x) = cell[d];
            // The knots inside the axis are its points without `skip` at each
            // end, so the knot interval of the point is ahead of its data
            // interval by `degree + 1 - skip`, and stays inside the knots.
            let skip = degree.div_ceil(2);
            let inner = self.axes[d].len() - degree - 1;
            let span = degree + (interval + 1).saturating_sub(skip).min(inner);
            basis[d] = basis_values::<T, N>(&spline.knots[d], degree, span, x);
            first[d] = (span - degree) * self.strides[d];
        }
        // The sum is nested, one level an axis:
        //
        //   sum over i0 of b0[i0] * (sum over i1 of b1[i1] * (... the sum along
        //   the last axis of b[i] * coefficient))
        //
        // `partial[d]` is the running sum of level `d`: over the splines of
        // axis `d` taken so far, for the splines chosen on the axes before it.
        // The coefficients along the last axis lie side by side. A sum of one
        // level does not wait for the sum of another, where one running sum
        // over all the terms would make every term wait for the one before it.
        let mut partial = [V::zero(); D];
        let mut chosen = [0; D];
        loop {
            let mut offset = first[last];
            for d in 0..last {
                offset += first[d] + chosen[d] * self.strides[d];
            }
            let row = &spline.coefficients[offset..offset + N];
            let mut done = V::zero();
            for (coefficient, spline) in row.iter().zip(&basis[last]) {
                done = done.axpy(*spline, *coefficient);
            }
            // `done` is a finished sum of the level below `d`. It goes into
            // level `d`; when that completes level `d` too, its sum goes up
            // in turn.
            let mut d = last;
            loop {
                if d == 0 {
                    return done;
                }
                d -= 1;
                partial[d] = partial[d].axpy(basis[d][chosen[d]], done);
                chosen[d] += 1;
                if chosen[d] <= degree {
                    break;
                }
                done = partial[d];
                partial[d] = V::zero();
                chosen[d] = 0;
            }
        }
    }

    /// The value of the shape-preserving cubic applied along the last axis,
    /// then along each axis before it.
    ///
    /// The cubic along an axis depends, in an interval, on a window of four
    /// points at most ([`pchip::Local`]). So the lines along the last axis
    /// that matter are those through the windows of the earlier axes:
    /// `4^(D - 1)` at most. Each is evaluated from its values and its slopes,
    /// which are kept. `folded[d]` collects the results along the window of
    /// axis `d`, for the points chosen on the axes before it. A window that is
    /// full is folded into one value, which goes into the window of the axis
    /// before: the levels of [`spline_sum`](Self::spline_sum), with a cubic
    /// through the values of a level where that has a sum.
    fn pchip_at(&self, slopes: &[V], cell: &[(usize, T); D]) -> V {
        let Some(last) = D.checked_sub(1) else {
            // No axes: the one value.
            return self.values[0];
        };
        let windows: [pchip::Local<T>; D] =
            core::array::from_fn(|d| pchip::Local::new(&self.axes[d], cell[d].0, cell[d].1));
        let mut folded = [[V::zero(); 4]; D];
        let mut chosen = [0; D];
        loop {
            // The lower point of the interval on the line chosen.
            let mut offset = cell[last].0;
            for d in 0..last {
                offset += (windows[d].first + chosen[d]) * self.strides[d];
            }
            let mut done = windows[last].eval_between(
                [self.values[offset], self.values[offset + 1]],
                [slopes[offset], slopes[offset + 1]],
            );
            let mut d = last;
            loop {
                if d == 0 {
                    return done;
                }
                d -= 1;
                folded[d][chosen[d]] = done;
                chosen[d] += 1;
                if chosen[d] < windows[d].len {
                    break;
                }
                done = windows[d].eval(&folded[d]);
                chosen[d] = 0;
            }
        }
    }
}

/// Where the construction of a method overflowed, in the grid as it is kept.
enum Overflow {
    /// A number worked out for this place of the values (a coefficient, a
    /// slope) is not finite.
    Value(usize),
    /// The matrix of this axis could not be set up, from this row on.
    Axis { axis: usize, row: usize },
}

impl Overflow {
    /// The error for the caller, with the place in the grid as it was given:
    /// along an axis that was given in descending order, point `i` of `n` was
    /// given as point `n - 1 - i`.
    fn into_given<T, const D: usize>(
        self,
        axes: &[Vec<T>; D],
        strides: &[usize; D],
        descending: &[bool; D],
    ) -> InterpolateError {
        let along = |d: usize, i: usize| {
            if descending[d] {
                axes[d].len() - 1 - i
            } else {
                i
            }
        };
        match self {
            Self::Value(index) => InterpolateError::NonFinite {
                index: (0..D)
                    .map(|d| along(d, index / strides[d] % axes[d].len()) * strides[d])
                    .sum(),
            },
            Self::Axis { axis, row } => InterpolateError::GridAxis {
                axis,
                fault: AxisFault::NonFinite {
                    index: along(axis, row),
                },
            },
        }
    }
}

/// The slopes of the shape-preserving cubic along every line of the last
/// axis, in the order of the values.
///
/// They are those of a [`PchipInterpolator`](crate::PchipInterpolator) through
/// the line. The slopes along the other axes are those of values folded
/// along the later axes, which differ from one point of evaluation to the
/// next: they are found at evaluation.
///
/// # Errors
///
/// Returns the place in the grid where a secant along the last axis
/// overflows (the later of its two points), or a slope.
fn slopes_along_last<T: Scalar, const D: usize, V: Value<T>>(
    axes: &[Vec<T>; D],
    values: &[V],
) -> Result<Vec<V>, Overflow> {
    let Some(axis) = axes.last() else {
        // No axes: one value, and no line.
        return Ok(Vec::new());
    };
    let n = axis.len();
    let mut slopes = vec![V::zero(); values.len()];
    // The same for every line.
    let per_width = hermite::per_widths(axis);
    let lines = values.chunks_exact(n).zip(slopes.chunks_exact_mut(n));
    for (l, (line, slopes)) in lines.enumerate() {
        pchip::slopes(axis, &per_width, line, slopes)
            .map_err(|index| Overflow::Value(l * n + index))?;
    }
    match slopes.iter().position(|slope| !slope.all_finite()) {
        Some(index) => Err(Overflow::Value(index)),
        None => Ok(slopes),
    }
}

/// Solves for the coefficients of the tensor-product spline of the given
/// degree through the values of a grid.
///
/// Interpolation along one axis commutes with interpolation along another, so
/// the coefficients are found axis by axis: along axis 0 every line of the
/// grid is replaced by the control points of its spline, then along axis 1,
/// and so on. All lines along an axis share one collocation matrix, factored
/// once. A line along an axis that is not the last is not contiguous: its
/// entries are a slab of the later axes apart, and the lines of a block are
/// solved together, slab by slab.
///
/// # Errors
///
/// Returns the axis whose matrix cannot be factored, with the row (the matrix
/// of an axis that passed the checks is regular; its entries are not finite
/// when the coordinates are so far apart that their distance overflows), or
/// the place of a coefficient that overflows.
fn tensor<T: Scalar, const D: usize, V: Value<T>>(
    axes: &[Vec<T>; D],
    values: &[V],
    strides: &[usize; D],
    degree: usize,
) -> Result<Tensor<T, D, V>, Overflow> {
    let mut coefficients = values.to_vec();
    let mut knots: [Vec<T>; D] = core::array::from_fn(|_| Vec::new());
    for (d, axis) in axes.iter().enumerate() {
        // An odd degree has its knots at the points of the axis.
        let knots_of_axis = Collocation::knots_at_points(axis, degree, false, false);
        let system =
            Collocation::factor(axis, knots_of_axis, degree, None, None).map_err(|at| {
                Overflow::Axis {
                    axis: d,
                    row: at.row,
                }
            })?;
        if strides[d] == 1 {
            solve_lines(&system, &mut coefficients, axis.len());
        } else {
            for block in coefficients.chunks_exact_mut(axis.len() * strides[d]) {
                system.solve_slabs(block, strides[d]);
            }
        }
        knots[d] = system.knots;
    }
    if let Some(index) = coefficients.iter().position(|c| !c.all_finite()) {
        return Err(Overflow::Value(index));
    }
    Ok(Tensor {
        degree,
        knots,
        coefficients,
    })
}

/// How many lines along the last axis are solved together.
///
/// Enough for the steps of one line not to wait for each other, few enough for
/// the lines to stay in the cache while they are solved. Measured with 8, 32,
/// and 128 on `benches/grid.rs`: more lines are a little faster on a grid of
/// 100 by 100 (52, 48, and 44 microseconds) and a little slower on one of 400
/// by 400 (1.36, 1.43, and 1.53 milliseconds). Any of them is three times as
/// fast as one line at a time.
const LINES_AT_ONCE: usize = 32;

/// Solves the system of the last axis for every line along it. A line is `n`
/// values side by side.
///
/// Solved one by one, the lines would each be a chain of steps that wait for
/// one another. So the lines are taken a few at a time and laid out the other
/// way round, value `i` of every line side by side: then one step works on all
/// the lines at once, as it does along the other axes, where the grid is laid
/// out that way already. The results are those of solving each line alone.
fn solve_lines<T: Scalar, V: Value<T>>(system: &Collocation<T>, values: &mut [V], n: usize) {
    let mut turned = vec![V::zero(); n * LINES_AT_ONCE];
    for group in values.chunks_mut(n * LINES_AT_ONCE) {
        let lines = group.len() / n;
        let turned = &mut turned[..n * lines];
        for (l, line) in group.chunks_exact(n).enumerate() {
            for (slab, value) in turned.chunks_exact_mut(lines).zip(line) {
                slab[l] = *value;
            }
        }
        system.solve_slabs(turned, lines);
        for (l, line) in group.chunks_exact_mut(n).enumerate() {
            for (slab, value) in turned.chunks_exact(lines).zip(line) {
                *value = slab[l];
            }
        }
    }
}

/// How far `x` is from point `i` of the axis toward point `i + 1`, as a
/// fraction of their distance: 0 at the one, 1 at the other, and beyond that
/// range for an `x` the mode lets continue past an end.
fn fraction<T: Scalar>(axis: &[T], i: usize, x: T) -> T {
    (x - axis[i]) / (axis[i + 1] - axis[i])
}

/// Whether an axis of finite coordinates, two at least, decreases strictly.
/// The first two points give the direction.
///
/// # Errors
///
/// Returns the fault of the first point that is equal to the one before it, or
/// goes against the direction.
fn decreases<T: Scalar>(axis: &[T]) -> Result<bool, AxisFault> {
    let descending = axis[1] < axis[0];
    for (i, pair) in axis.windows(2).enumerate() {
        let (before, after) = if descending {
            (pair[1], pair[0])
        } else {
            (pair[0], pair[1])
        };
        if before < after {
            continue;
        }
        // Finite, and not in the direction: equal, or against it.
        return Err(if after < before {
            AxisFault::NotMonotonic { index: i + 1 }
        } else {
            AxisFault::DuplicatePoints {
                first: i,
                second: i + 1,
            }
        });
    }
    Ok(descending)
}

/// Reverses the order of the values along one axis of `n` points, whose
/// neighbors are `stride` apart.
fn reverse_along<V>(values: &mut [V], n: usize, stride: usize) {
    for block in values.chunks_exact_mut(n * stride) {
        for i in 0..n / 2 {
            let (head, tail) = block.split_at_mut((n - 1 - i) * stride);
            head[i * stride..(i + 1) * stride].swap_with_slice(&mut tail[..stride]);
        }
    }
}

impl<T: Scalar, const D: usize, V: Value<T>> Interpolator for RegularGridInterpolator<T, D, V> {
    type Point = [T; D];
    type Output = V;

    fn eval(&self, point: &[T; D]) -> Option<V> {
        let mut cell = [(0, T::zero()); D];
        for (d, place) in cell.iter_mut().enumerate() {
            *place = search::eval(&self.axes[d], self.extrapolate, point[d], |i, x| (i, x))?;
        }
        Some(self.at(&cell))
    }

    /// Evaluates at every point, by runs: with a hint for each axis where the
    /// points of a run move in one direction along every axis, a cell or less
    /// at a time (a path, or the rows of an image); plainly where they do not.
    /// Either way the results are bit for bit those of a loop over `eval`.
    fn eval_batch_into(
        &self,
        points: &[[T; D]],
        out: &mut [V],
        fill: V,
    ) -> Result<(), InterpolateError> {
        if points.len() != out.len() {
            return Err(InterpolateError::LengthMismatch {
                expected: points.len(),
                actual: out.len(),
            });
        }
        let mut hints = [Hint::new(); D];
        let runs = points.chunks(search::BATCH_RUN);
        for (run, slots) in runs.zip(out.chunks_mut(search::BATCH_RUN)) {
            let hinted = run.len() >= search::MIN_HINT_RUN
                && self.axes.iter().enumerate().all(|(d, axis)| {
                    search::hint_pays_along(axis, run.iter().map(|point| point[d]))
                });
            for (point, slot) in run.iter().zip(slots) {
                let value = if hinted {
                    self.eval_hinted(point, &mut hints)
                } else {
                    self.eval(point)
                };
                *slot = value.unwrap_or(fill);
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::{GridMethod, RegularGridInterpolator};
    use crate::common::{assert_close, load, reals, Real, Tol};
    use serde_json::Value as Json;

    /// The tolerance for a coefficient, as a share of the largest one.
    ///
    /// The coefficients are the solution of one linear system an axis, and
    /// carry the condition number of the worst axis. Against the coefficients
    /// in exact rational arithmetic, as a share of the largest:
    ///
    /// | case | worst condition number | SciPy | this crate | apart |
    /// | --- | --- | --- | --- | --- |
    /// | 2D cubic | 320 | 2.6e-13 | 2.5e-13 | 7.6e-15 |
    /// | 2D quintic | 1100 | 8.9e-14 | 2.0e-14 | 7.6e-14 |
    /// | 3D cubic | 200 | 1.1e-13 | 7.1e-14 | 6.4e-14 |
    /// | 3D quintic | 1400 | 2.2e-12 | 2.1e-12 | 1.6e-12 |
    ///
    /// On the last, both are farther from the exact coefficients than the
    /// tolerance for a direct method (1e-12), and as far as the condition
    /// number lets either be. The values of the spline, which are what the
    /// integration tests compare, agree with SciPy's to 1e-12 on every case.
    const COEFFICIENTS: Tol = Tol {
        rtol: 1e-11,
        atol: 1e-11,
    };

    fn nums(v: &Json) -> Vec<f64> {
        reals(&serde_json::from_value::<Vec<Real>>(v.clone()).unwrap_or_default())
    }

    /// The coefficients a grid in `D` dimensions solved for, and SciPy's.
    fn check<const D: usize>(case: &Json, method: GridMethod) {
        let name = case["name"].as_str().expect("name");
        let axes: Vec<Vec<f64>> = case["axes"]
            .as_array()
            .expect("axes")
            .iter()
            .map(nums)
            .collect();
        let slices: [&[f64]; D] = core::array::from_fn(|d| axes[d].as_slice());
        let grid = RegularGridInterpolator::new(slices, &nums(&case["values"]), method)
            .expect("valid grid");
        let spline = grid.rule.spline().expect("a spline method keeps a spline");
        // The methods that solve nothing keep none.
        for plain in [GridMethod::Nearest, GridMethod::Linear] {
            let grid = RegularGridInterpolator::new(slices, &nums(&case["values"]), plain);
            assert!(grid.expect("valid grid").rule.spline().is_none(), "{name}");
        }
        let want = nums(&case["coefficients"]);
        assert_eq!(spline.coefficients.len(), want.len(), "{name}");
        // A coefficient is compared against the size of the largest one, not
        // against its own: the solve mixes them, and what it rounds is the
        // large ones.
        let largest = want.iter().fold(0.0_f64, |m, c| m.max(c.abs()));
        for (i, (got, want)) in spline.coefficients.iter().zip(&want).enumerate() {
            let what = format!("{name} coefficient {i}, largest {largest}");
            assert_close(got / largest, want / largest, COEFFICIENTS, &what);
        }
    }

    #[test]
    fn the_coefficients_are_those_of_scipys_splines_applied_axis_by_axis() {
        let file = load::<Json>("grid").expect("golden file grid");
        let (mut cases, mut flat) = (0, 0);
        for case in &file.cases {
            if case.get("coefficients").is_none() {
                continue;
            }
            let named = case["options"]["method"].as_str();
            let (_, method) = [
                ("cubic", GridMethod::Cubic),
                ("quintic", GridMethod::Quintic),
            ]
            .into_iter()
            .find(|(name, _)| Some(*name) == named)
            .expect("a case with coefficients is cubic or quintic");
            if case["axes"].as_array().map_or(0, Vec::len) == 2 {
                check::<2>(case, method);
                flat += 1;
            } else {
                check::<3>(case, method);
            }
            cases += 1;
        }
        // Cubic and quintic, in 2D and in 3D.
        assert_eq!((cases, flat), (4, 2));
    }
}
