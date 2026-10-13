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
    finite(coordinates, |c| c.is_finite())?;
    finite(values, V::all_finite)?;
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
    finite(coordinates, |c| c.is_finite())?;
    finite(values, V::all_finite)?;
    finite(slopes, V::all_finite)?;
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

/// How many entries are looked at without leaving the loop. A loop that may
/// leave at every entry is compiled one entry at a time; over a run of this
/// many the compiler takes several at once, and data without a fault, which
/// is nearly all data, is checked several times as fast (the values of a grid
/// of 40 x 40 x 40: 11 µs before, measured in #122).
const RUN: usize = 64;

/// The first entry of `entries` that `holds` is false for.
pub(crate) fn first_failing<V>(entries: &[V], holds: impl Fn(&V) -> bool) -> Option<usize> {
    for (run, chunk) in entries.chunks(RUN).enumerate() {
        if !chunk.iter().fold(true, |all, entry| all & holds(entry)) {
            let within = chunk.iter().position(|entry| !holds(entry));
            return within.map(|at| run * RUN + at);
        }
    }
    None
}

/// Step 3 for one input.
fn finite<V>(entries: &[V], is_finite: impl Fn(&V) -> bool) -> Result<(), InterpolateError> {
    match first_failing(entries, is_finite) {
        Some(index) => Err(InterpolateError::NonFinite { index }),
        None => Ok(()),
    }
}

/// Checks that finite coordinates increase strictly.
fn strictly_increasing<T: Scalar>(coordinates: &[T]) -> Result<(), InterpolateError> {
    // Runs that share their end with the next run, so that every pair is in
    // one of them.
    let mut start = 0;
    while start + 1 < coordinates.len() {
        let end = (start + RUN).min(coordinates.len() - 1);
        let run = &coordinates[start..=end];
        if !run
            .windows(2)
            .fold(true, |all, pair| all & (pair[0] < pair[1]))
        {
            break;
        }
        start = end;
    }
    for (i, pair) in coordinates.windows(2).enumerate().skip(start) {
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
    fn a_fault_is_found_wherever_it_is_in_a_run() {
        // The checks go by runs of 64. A fault at every place of 200 entries:
        // the first and the last of a run, the last entry of all, and a pair
        // that two runs share.
        let good: Vec<f64> = (0..200).map(f64::from).collect();
        assert_eq!(check(&good, &good, 200, 2), Ok(()));
        assert_eq!(super::first_failing(&good, |v| v.is_finite()), None);
        assert_eq!(super::first_failing(&[] as &[f64], |v| v.is_finite()), None);
        for at in 0..200 {
            let mut bad = good.clone();
            bad[at] = f64::NAN;
            assert_eq!(check(&bad, &good, 200, 2), Err(E::NonFinite { index: at }));
            assert_eq!(check(&good, &bad, 200, 2), Err(E::NonFinite { index: at }));
            // The first of two faults is the one named.
            bad[199] = f64::INFINITY;
            assert_eq!(check(&good, &bad, 200, 2), Err(E::NonFinite { index: at }));
        }
        for at in 1..200 {
            let mut repeated = good.clone();
            repeated[at] = repeated[at - 1];
            let fault = E::DuplicatePoints {
                first: at - 1,
                second: at,
            };
            assert_eq!(check(&repeated, &good, 200, 2), Err(fault));
            let mut back = good.clone();
            back[at] = back[at - 1] - 0.5;
            assert_eq!(check(&back, &good, 200, 2), Err(E::NotSorted { index: at }));
            // The first of two faults is the one named.
            back[199] = -1.0;
            assert_eq!(check(&back, &good, 200, 2), Err(E::NotSorted { index: at }));
        }
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
