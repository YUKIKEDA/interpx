//! Input checks shared by the 1D constructors, in the fixed order of the design
//! (docs/design.ja.md section 3.5).
//!
//! An input with several faults returns the error of the earliest step:
//!
//! 1. `LengthMismatch`: the values against the length the coordinates fix.
//! 2. `InsufficientPoints`.
//! 3. `NonFinite`: the coordinates first, then the values.
//! 4. `DuplicatePoints` or `NotSorted`, whichever comes first.
//!
//! A second input per point (the slopes of `CubicHermiteSpline`) is checked
//! after the values within each step. The order of the steps stays.

use crate::{InterpolateError, Scalar, Value};

/// Checks coordinates and the values that belong to them.
///
/// `expected_values` is the number of values the coordinates call for (one per
/// point for data, several per interval for coefficients). `required` is the
/// smallest number of coordinates the caller accepts.
pub(crate) fn coordinates_and_values<T: Scalar, V: Value<T>>(
    coordinates: &[T],
    values: &[V],
    expected_values: usize,
    required: usize,
) -> Result<(), InterpolateError> {
    length(values, expected_values)?;
    enough(coordinates, required)?;
    finite(coordinates.iter().map(|c| c.is_finite()))?;
    finite(values.iter().map(V::all_finite))?;
    strictly_increasing(coordinates)
}

/// Checks coordinates with the values and the slopes at them, one of each per
/// point.
///
/// `required` is the smallest number of coordinates the caller accepts.
pub(crate) fn coordinates_values_and_slopes<T: Scalar, V: Value<T>>(
    coordinates: &[T],
    values: &[V],
    slopes: &[V],
    required: usize,
) -> Result<(), InterpolateError> {
    length(values, coordinates.len())?;
    length(slopes, coordinates.len())?;
    enough(coordinates, required)?;
    finite(coordinates.iter().map(|c| c.is_finite()))?;
    finite(values.iter().map(V::all_finite))?;
    finite(slopes.iter().map(V::all_finite))?;
    strictly_increasing(coordinates)
}

/// Step 1 for one input: its length against the one the coordinates fix.
///
/// On its own for a constructor that checks its other arguments between this
/// step and the next (`BSpline1D`: the fewest points depend on the degree).
pub(crate) fn length<V>(input: &[V], expected: usize) -> Result<(), InterpolateError> {
    if input.len() == expected {
        Ok(())
    } else {
        Err(InterpolateError::LengthMismatch {
            expected,
            actual: input.len(),
        })
    }
}

/// Step 2: the number of coordinates.
///
/// On its own for a constructor that has no coordinates yet (`BSplineCurve`
/// makes them from its points).
pub(crate) fn enough<T>(coordinates: &[T], required: usize) -> Result<(), InterpolateError> {
    if coordinates.len() < required {
        Err(InterpolateError::InsufficientPoints {
            required,
            actual: coordinates.len(),
        })
    } else {
        Ok(())
    }
}

/// Step 3 for one input, given as whether each entry is finite.
fn finite(mut entries: impl Iterator<Item = bool>) -> Result<(), InterpolateError> {
    match entries.position(|is_finite| !is_finite) {
        Some(index) => Err(InterpolateError::NonFinite { index }),
        None => Ok(()),
    }
}

/// Checks that finite coordinates increase strictly.
fn strictly_increasing<T: Scalar>(coordinates: &[T]) -> Result<(), InterpolateError> {
    for (i, pair) in coordinates.windows(2).enumerate() {
        if pair[0] < pair[1] {
            continue;
        }
        // Not increasing here: the coordinates are finite, so the pair is either
        // decreasing or equal.
        return Err(if pair[1] < pair[0] {
            InterpolateError::NotSorted { index: i + 1 }
        } else {
            InterpolateError::DuplicatePoints {
                first: i,
                second: i + 1,
            }
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::coordinates_and_values as check;
    use super::coordinates_values_and_slopes as with_slopes;
    use crate::InterpolateError as E;

    #[test]
    fn good_input_passes() {
        assert_eq!(check(&[0.0, 1.0, 2.5], &[1.0, 2.0, 3.0], 3, 2), Ok(()));
        assert_eq!(check(&[0.0_f32, 1.0], &[[1.0_f32, 2.0]; 2], 2, 2), Ok(()));
        // A single coordinate is in order by definition.
        assert_eq!(check(&[0.0], &[1.0], 1, 1), Ok(()));
    }

    #[test]
    fn each_fault_has_its_error() {
        assert_eq!(
            check(&[0.0, 1.0], &[1.0], 2, 2),
            Err(E::LengthMismatch {
                expected: 2,
                actual: 1
            })
        );
        assert_eq!(
            check(&[0.0], &[1.0], 1, 2),
            Err(E::InsufficientPoints {
                required: 2,
                actual: 1
            })
        );
        assert_eq!(
            check(&[0.0, f64::NAN, 2.0], &[1.0; 3], 3, 2),
            Err(E::NonFinite { index: 1 })
        );
        assert_eq!(
            check(&[0.0, 1.0, f64::INFINITY], &[1.0; 3], 3, 2),
            Err(E::NonFinite { index: 2 })
        );
        assert_eq!(
            check(&[0.0, 1.0, 2.0], &[1.0, 2.0, f64::NEG_INFINITY], 3, 2),
            Err(E::NonFinite { index: 2 })
        );
        assert_eq!(
            check(&[0.0, 1.0], &[[1.0, 2.0], [f64::NAN, 4.0]], 2, 2),
            Err(E::NonFinite { index: 1 })
        );
        assert_eq!(
            check(&[0.0, 2.0, 1.0], &[1.0; 3], 3, 2),
            Err(E::NotSorted { index: 2 })
        );
        assert_eq!(
            check(&[0.0, 1.0, 1.0, 2.0], &[1.0; 4], 4, 2),
            Err(E::DuplicatePoints {
                first: 1,
                second: 2
            })
        );
    }

    #[test]
    fn two_adjacent_faults_return_the_earlier_step() {
        // Step 1 before step 2: the wrong length and too few points.
        assert_eq!(
            check(&[0.0], &[1.0, 2.0], 1, 2),
            Err(E::LengthMismatch {
                expected: 1,
                actual: 2
            })
        );
        // Step 2 before step 3: too few points, and the one point is NaN.
        assert_eq!(
            check(&[f64::NAN], &[1.0], 1, 2),
            Err(E::InsufficientPoints {
                required: 2,
                actual: 1
            })
        );
        // Step 3, coordinates before values: both hold a NaN, the value's comes first
        // by index, and the coordinate is still the one reported.
        assert_eq!(
            check(&[0.0, 1.0, f64::NAN], &[f64::NAN, 1.0, 1.0], 3, 2),
            Err(E::NonFinite { index: 2 })
        );
        // Step 3 before step 4: unsorted coordinates and a NaN value.
        assert_eq!(
            check(&[1.0, 0.0], &[1.0, f64::NAN], 2, 2),
            Err(E::NonFinite { index: 1 })
        );
    }

    #[test]
    fn slopes_are_checked_after_the_values_within_each_step() {
        let ok = [1.0; 3];
        assert_eq!(
            with_slopes(&[0.0, 1.0, 2.5], &ok, &[0.5, -1.0, 0.0], 2),
            Ok(())
        );
        assert_eq!(
            with_slopes(&[0.0_f32, 1.0], &[[1.0_f32, 2.0]; 2], &[[0.0_f32; 2]; 2], 2),
            Ok(())
        );
        // Step 1: the length of the slopes, after that of the values.
        assert_eq!(
            with_slopes(&[0.0, 1.0, 2.0], &ok, &[1.0; 2], 2),
            Err(E::LengthMismatch {
                expected: 3,
                actual: 2
            })
        );
        assert_eq!(
            with_slopes(&[0.0, 1.0, 2.0], &[1.0; 4], &[1.0; 2], 2),
            Err(E::LengthMismatch {
                expected: 3,
                actual: 4
            })
        );
        // Step 1 before step 2, for the slopes as well.
        assert_eq!(
            with_slopes(&[0.0], &[1.0], &[1.0, 2.0], 2),
            Err(E::LengthMismatch {
                expected: 1,
                actual: 2
            })
        );
        assert_eq!(
            with_slopes(&[0.0], &[1.0], &[1.0], 2),
            Err(E::InsufficientPoints {
                required: 2,
                actual: 1
            })
        );
        // Step 3: coordinates, values, slopes, whatever the positions.
        assert_eq!(
            with_slopes(&[0.0, 1.0, 2.0], &ok, &[0.0, f64::NAN, f64::INFINITY], 2),
            Err(E::NonFinite { index: 1 })
        );
        assert_eq!(
            with_slopes(
                &[0.0, 1.0],
                &[[1.0, 2.0]; 2],
                &[[0.0, 0.0], [0.0, f64::NEG_INFINITY]],
                2
            ),
            Err(E::NonFinite { index: 1 })
        );
        assert_eq!(
            with_slopes(
                &[0.0, 1.0, 2.0],
                &[1.0, 1.0, f64::NAN],
                &[f64::NAN, 0.0, 0.0],
                2
            ),
            Err(E::NonFinite { index: 2 })
        );
        assert_eq!(
            with_slopes(
                &[0.0, f64::NAN, 2.0],
                &[1.0, 1.0, f64::NAN],
                &[f64::NAN, 0.0, 0.0],
                2
            ),
            Err(E::NonFinite { index: 1 })
        );
        // Step 3 before step 4: unsorted coordinates and a slope that is NaN.
        assert_eq!(
            with_slopes(&[1.0, 0.0, 2.0], &ok, &[0.0, 0.0, f64::NAN], 2),
            Err(E::NonFinite { index: 2 })
        );
        // Step 4.
        assert_eq!(
            with_slopes(&[1.0, 0.0, 2.0], &ok, &ok, 2),
            Err(E::NotSorted { index: 1 })
        );
        assert_eq!(
            with_slopes(&[0.0, 0.0, 2.0], &ok, &ok, 2),
            Err(E::DuplicatePoints {
                first: 0,
                second: 1
            })
        );
    }

    #[test]
    fn of_a_duplicate_and_a_decrease_the_first_one_wins() {
        assert_eq!(
            check(&[0.0, 0.0, -1.0], &[1.0; 3], 3, 2),
            Err(E::DuplicatePoints {
                first: 0,
                second: 1
            })
        );
        assert_eq!(
            check(&[0.0, -1.0, -1.0], &[1.0; 3], 3, 2),
            Err(E::NotSorted { index: 1 })
        );
    }
}
