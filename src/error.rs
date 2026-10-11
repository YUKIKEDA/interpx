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

    /// The value is NaN or infinite.
    #[error("not finite")]
    NonFinite,
}

#[cfg(test)]
mod tests {
    use super::InterpolateError as E;

    #[test]
    fn every_variant_has_its_message() {
        let cases = [
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
