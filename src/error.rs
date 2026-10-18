//! The error type of this crate.

/// The reason an interpolator could not be built, or a batch could not be written.
///
/// Evaluation at a single point does not use this type: a point outside the
/// domain, or a NaN query, gives `None`.
///
/// The enum is `#[non_exhaustive]`: later versions add variants, so a `match`
/// needs a wildcard arm.
///
/// # Examples
///
/// ```
/// use interpx::InterpolateError;
///
/// let err = InterpolateError::LengthMismatch { expected: 3, actual: 2 };
/// assert_eq!(err.to_string(), "length mismatch: expected 3, got 2");
/// assert!(matches!(err, InterpolateError::LengthMismatch { actual: 2, .. }));
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, thiserror::Error)]
#[non_exhaustive]
pub enum InterpolateError {
    /// There are fewer data points than the method needs.
    #[error("insufficient points: required {required}, got {actual}")]
    InsufficientPoints {
        /// The smallest number of points the method accepts.
        required: usize,
        /// The number of points given.
        actual: usize,
    },

    /// Two inputs that must have the same length do not.
    ///
    /// For data, `expected` is the number of coordinates and `actual` the number
    /// of values. For a batch, `expected` is the number of query points and
    /// `actual` the length of the output buffer.
    #[error("length mismatch: expected {expected}, got {actual}")]
    LengthMismatch {
        /// The length the other input fixes.
        expected: usize,
        /// The length found.
        actual: usize,
    },

    /// A coordinate or a value is NaN or infinite.
    ///
    /// Coordinates are checked before values. `index` is a position among the
    /// coordinates when one of them is at fault, and a position among the values
    /// only when every coordinate is finite. The variant does not say which of
    /// the two it is: check the coordinates to tell.
    #[error("non-finite value at index {index}")]
    NonFinite {
        /// The index of the first non-finite entry, in the coordinates or, when
        /// those are all finite, in the values.
        index: usize,
    },

    /// The coordinates decrease somewhere.
    #[error("coordinates are not strictly increasing at index {index}")]
    NotSorted {
        /// The index of the first coordinate that is smaller than the one before it.
        index: usize,
    },

    /// Two data points have the same coordinates.
    #[error("duplicate points at indices {first} and {second}")]
    DuplicatePoints {
        /// The index of the first of the two points.
        first: usize,
        /// The index of the second of the two points.
        second: usize,
    },

    /// A periodic boundary condition was asked for, but the first and the last
    /// values differ.
    #[error("periodic boundary needs equal values at both ends")]
    NotPeriodic,

    /// An argument that is not data has a value the method cannot use.
    #[error("invalid parameter `{name}`: {reason}")]
    InvalidParameter {
        /// The name of the argument.
        name: &'static str,
        /// What is wrong with its value.
        reason: ParamReason,
    },

    /// A system of equations the method has to solve could not be solved.
    #[error("the solver failed: {0}")]
    SolverFailed(SolverError),

    /// An allocation that can be far larger than the input that asked for it
    /// could not be made: there is not that much memory.
    ///
    /// A matrix with the square of the number of points as its size, or the
    /// coefficients of an antiderivative of a high order. Smaller input, or
    /// another method (the `neighbors` of an RBF interpolant), is the way out.
    #[error("an allocation of {bytes} bytes failed")]
    AllocationFailed {
        /// The size that was asked for.
        bytes: usize,
    },

    /// The points cannot be divided into simplices: they do not span their
    /// dimensions, or the exact arithmetic of the triangulation could not be
    /// carried out.
    #[error("the points cannot be triangulated: {0}")]
    DegenerateGeometry(DegeneracyKind),

    /// One axis of a grid cannot be used.
    #[error("axis {axis} of the grid: {fault}")]
    GridAxis {
        /// The axis, counted from 0.
        axis: usize,
        /// What is wrong with it.
        fault: AxisFault,
    },
}

/// What is wrong with the axis named by [`InterpolateError::GridAxis`].
///
/// Positions are those in the axis as it was given.
///
/// The enum is `#[non_exhaustive]`: later versions add faults.
///
/// # Examples
///
/// ```
/// use interpx::{AxisFault, GridMethod, InterpolateError, RegularGridInterpolator};
///
/// // The second axis has the same coordinate twice.
/// let (x, y) = ([0.0, 1.0], [0.0, 0.5, 0.5]);
/// let err = RegularGridInterpolator::new([&x[..], &y[..]], &[0.0; 6], GridMethod::Linear).err();
/// assert_eq!(
///     err,
///     Some(InterpolateError::GridAxis {
///         axis: 1,
///         fault: AxisFault::DuplicatePoints { first: 1, second: 2 },
///     }),
/// );
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, thiserror::Error)]
#[non_exhaustive]
pub enum AxisFault {
    /// The axis has fewer points than the method needs.
    #[error("insufficient points: required {required}, got {actual}")]
    InsufficientPoints {
        /// The fewest points the method takes on an axis.
        required: usize,
        /// The number of points of the axis.
        actual: usize,
    },

    /// A coordinate is NaN or infinite.
    #[error("non-finite coordinate at index {index}")]
    NonFinite {
        /// The position of the first such coordinate.
        index: usize,
    },

    /// The axis neither increases all the way nor decreases all the way.
    #[error("the coordinates turn around at index {index}")]
    NotMonotonic {
        /// The position of the first coordinate that goes against the
        /// direction of the first two.
        index: usize,
    },

    /// Two neighboring coordinates are equal.
    #[error("equal coordinates at indices {first} and {second}")]
    DuplicatePoints {
        /// The position of the first of the two.
        first: usize,
        /// The position of the second of the two.
        second: usize,
    },
}

/// What is wrong with the argument named by
/// [`InterpolateError::InvalidParameter`].
///
/// The enum is `#[non_exhaustive]`: later versions add reasons.
///
/// # Examples
///
/// ```
/// use interpx::{BSpline, InterpolateError, ParamReason};
///
/// // A B-spline of degree 16: one more than the largest degree it supports.
/// let err = BSpline::new(&[0.0; 34], &[0.0; 17], 16).err();
/// assert_eq!(
///     err,
///     Some(InterpolateError::InvalidParameter {
///         name: "degree",
///         reason: ParamReason::TooLarge { max: 15 },
///     }),
/// );
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, thiserror::Error)]
#[non_exhaustive]
pub enum ParamReason {
    /// The value is larger than the largest one the method supports.
    #[error("larger than {max}")]
    TooLarge {
        /// The largest value the method supports.
        max: usize,
    },

    /// The value is smaller than the smallest one the method supports.
    #[error("smaller than {min}")]
    TooSmall {
        /// The smallest value the method supports.
        min: usize,
    },

    /// The value is NaN or infinite.
    #[error("not finite")]
    NonFinite,

    /// The value is a valid one, but not together with the other arguments.
    #[error("not supported with the other arguments")]
    Unsupported,

    /// The value is zero or negative, where it must be above zero.
    #[error("not above zero")]
    NotPositive,

    /// The value is negative, where it must be zero or above.
    #[error("negative")]
    Negative,
}

/// Why the points of [`InterpolateError::DegenerateGeometry`] cannot be
/// triangulated.
///
/// The enum is `#[non_exhaustive]`: later versions add reasons.
///
/// # Examples
///
/// ```
/// use interpx::{DegeneracyKind, InterpolateError};
///
/// // What `LinearNDInterpolator` returns for points on one line in the plane.
/// let err = InterpolateError::DegenerateGeometry(DegeneracyKind::LowerDimensional {
///     affine_dim: 1,
/// });
/// assert_eq!(
///     err.to_string(),
///     "the points cannot be triangulated: the points lie in an affine subspace of dimension 1",
/// );
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, thiserror::Error)]
#[non_exhaustive]
pub enum DegeneracyKind {
    /// The distinct points lie in an affine subspace of fewer dimensions than
    /// the points have: on a line in the plane, on a plane in space.
    #[error("the points lie in an affine subspace of dimension {affine_dim}")]
    LowerDimensional {
        /// The dimension of the smallest affine subspace that holds them.
        affine_dim: usize,
    },

    /// The exact arithmetic that decides the triangulation could not be
    /// carried out: its working space could not be allocated.
    #[error("the exact arithmetic could not be carried out")]
    NumericalLimit,
}

/// Why a system of equations could not be solved, for
/// [`InterpolateError::SolverFailed`].
///
/// The enum is `#[non_exhaustive]`: later versions add reasons.
///
/// # Examples
///
/// ```
/// use interpx::{InterpolateError, SolverError};
///
/// let err = InterpolateError::SolverFailed(SolverError::Singular);
/// assert_eq!(err.to_string(), "the solver failed: the system is singular");
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, thiserror::Error)]
#[non_exhaustive]
pub enum SolverError {
    /// The system has no solution, or more than one: its matrix is singular.
    #[error("the system is singular")]
    Singular,

    /// An iterative solution did not come within its tolerance in as many
    /// steps as it is given.
    #[error("the iteration did not converge")]
    NotConverged,
}

#[cfg(test)]
mod tests {
    use super::InterpolateError as E;

    #[test]
    fn every_variant_has_its_message() {
        let cases = [
            (
                E::AllocationFailed { bytes: 1 << 40 },
                "an allocation of 1099511627776 bytes failed",
            ),
            (
                E::DegenerateGeometry(super::DegeneracyKind::LowerDimensional { affine_dim: 1 }),
                "the points cannot be triangulated: the points lie in an affine subspace of dimension 1",
            ),
            (
                E::DegenerateGeometry(super::DegeneracyKind::NumericalLimit),
                "the points cannot be triangulated: the exact arithmetic could not be carried out",
            ),
            (
                E::SolverFailed(super::SolverError::Singular),
                "the solver failed: the system is singular",
            ),
            (
                E::InvalidParameter {
                    name: "epsilon",
                    reason: super::ParamReason::NotPositive,
                },
                "invalid parameter `epsilon`: not above zero",
            ),
            (
                E::InvalidParameter {
                    name: "smoothing",
                    reason: super::ParamReason::Negative,
                },
                "invalid parameter `smoothing`: negative",
            ),
            (
                E::InsufficientPoints {
                    required: 4,
                    actual: 2,
                },
                "insufficient points: required 4, got 2",
            ),
            (
                E::LengthMismatch {
                    expected: 3,
                    actual: 5,
                },
                "length mismatch: expected 3, got 5",
            ),
            (E::NonFinite { index: 7 }, "non-finite value at index 7"),
            (
                E::NotSorted { index: 2 },
                "coordinates are not strictly increasing at index 2",
            ),
            (
                E::DuplicatePoints {
                    first: 1,
                    second: 2,
                },
                "duplicate points at indices 1 and 2",
            ),
            (
                E::NotPeriodic,
                "periodic boundary needs equal values at both ends",
            ),
            (
                E::InvalidParameter {
                    name: "degree",
                    reason: super::ParamReason::TooLarge { max: 15 },
                },
                "invalid parameter `degree`: larger than 15",
            ),
            (
                E::InvalidParameter {
                    name: "bc",
                    reason: super::ParamReason::NonFinite,
                },
                "invalid parameter `bc`: not finite",
            ),
            (
                E::InvalidParameter {
                    name: "degree",
                    reason: super::ParamReason::TooSmall { min: 1 },
                },
                "invalid parameter `degree`: smaller than 1",
            ),
            (
                E::InvalidParameter {
                    name: "bc",
                    reason: super::ParamReason::Unsupported,
                },
                "invalid parameter `bc`: not supported with the other arguments",
            ),
            (
                E::GridAxis {
                    axis: 2,
                    fault: super::AxisFault::InsufficientPoints {
                        required: 2,
                        actual: 1,
                    },
                },
                "axis 2 of the grid: insufficient points: required 2, got 1",
            ),
            (
                E::GridAxis {
                    axis: 0,
                    fault: super::AxisFault::NonFinite { index: 3 },
                },
                "axis 0 of the grid: non-finite coordinate at index 3",
            ),
            (
                E::GridAxis {
                    axis: 1,
                    fault: super::AxisFault::NotMonotonic { index: 4 },
                },
                "axis 1 of the grid: the coordinates turn around at index 4",
            ),
            (
                E::GridAxis {
                    axis: 1,
                    fault: super::AxisFault::DuplicatePoints {
                        first: 0,
                        second: 1,
                    },
                },
                "axis 1 of the grid: equal coordinates at indices 0 and 1",
            ),
        ];
        for (err, message) in cases {
            assert_eq!(err.to_string(), message);
        }
    }

    #[test]
    fn it_is_a_std_error_and_compares_by_value() {
        fn as_error(e: &dyn std::error::Error) -> String {
            e.to_string()
        }
        let err = E::NonFinite { index: 0 };
        assert_eq!(as_error(&err), "non-finite value at index 0");
        assert_eq!(err, E::NonFinite { index: 0 });
        assert_ne!(err, E::NonFinite { index: 1 });
        assert_ne!(err, E::NotPeriodic);
    }
}
