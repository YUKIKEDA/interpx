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
    #[error("non-finite value at index {index}")]
    NonFinite {
        /// The index of the first non-finite entry.
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
