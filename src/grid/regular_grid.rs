//! Interpolation of values on a rectilinear grid in several dimensions.

use crate::error::AxisFault;
use crate::{search, Extrapolate, Hint, InterpolateError, Interpolator, Scalar, Value};

/// How a [`RegularGridInterpolator`] interpolates between the points of its
/// grid.
///
/// The method is chosen at construction and stays. The enum is
/// `#[non_exhaustive]`: later versions add the spline methods.
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
    method: GridMethod,
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
    ///   for the first axis with fewer than two points.
    /// - [`InterpolateError::GridAxis`] with [`AxisFault::NonFinite`] for the
    ///   first axis with a coordinate that is NaN or infinite.
    /// - [`InterpolateError::NonFinite`] when a value is NaN or infinite. Its
    ///   `index` is the position in `values`.
    /// - [`InterpolateError::GridAxis`] with [`AxisFault::DuplicatePoints`] or
    ///   [`AxisFault::NotMonotonic`] for the first axis that neither increases
    ///   strictly nor decreases strictly.
    ///
    /// Positions in an axis are those in the axis as it was given.
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
        for (d, axis) in axes.iter().enumerate() {
            if axis.len() < 2 {
                return Err(on_axis(d)(AxisFault::InsufficientPoints {
                    required: 2,
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
        Ok(Self {
            axes,
            values,
            strides,
            method,
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
        self.method
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
                    (i, fraction(axis, i, x))
                })?;
        }
        Some(self.at(&cell))
    }

    /// The value in the cell whose lower corner is at the given point of each
    /// axis, at the given fraction of the way to the next point.
    fn at(&self, cell: &[(usize, T); D]) -> V {
        match self.method {
            GridMethod::Nearest => {
                let half = (T::one() + T::one()).recip();
                let offset = cell
                    .iter()
                    .zip(&self.strides)
                    .map(|((lower, fraction), stride)| {
                        (lower + usize::from(*fraction > half)) * stride
                    })
                    .sum::<usize>();
                self.values[offset]
            }
            GridMethod::Linear => {
                // For each axis, the weight and the place in `values` of the
                // lower and of the upper side of the cell.
                let mut weights = [[T::zero(); 2]; D];
                let mut offsets = [[0; 2]; D];
                for (d, ((lower, fraction), stride)) in cell.iter().zip(&self.strides).enumerate() {
                    weights[d] = [T::one() - *fraction, *fraction];
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
            let axis = &self.axes[d];
            *place = search::eval(axis, self.extrapolate, point[d], |i, x| {
                (i, fraction(axis, i, x))
            })?;
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
