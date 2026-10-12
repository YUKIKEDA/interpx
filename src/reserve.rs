//! Allocations that can be far larger than the input that asked for them.
//!
//! A matrix with a row and a column for every point, or one more coefficient a
//! piece for every order of an antiderivative: a modest input can ask for more
//! memory than the machine has. Such an allocation is made here, so that its
//! failure is an error the caller gets, where `vec!` would end the process.
//!
//! Allocations in proportion to their input (a copy of a slice) are made as
//! they are. Two things are out of reach: allocations inside the dependencies,
//! and systems that promise more memory than they have, where an allocation
//! succeeds and the process is stopped when the memory is written to.

use crate::InterpolateError;

/// Makes room in `vec` for `len` entries in all, without changing what it
/// holds.
///
/// # Errors
///
/// Returns [`InterpolateError::AllocationFailed`] with the size asked for when
/// the memory cannot be had.
pub(crate) fn room<T>(vec: &mut Vec<T>, len: usize) -> Result<(), InterpolateError> {
    vec.try_reserve_exact(len.saturating_sub(vec.len()))
        .map_err(|_| InterpolateError::AllocationFailed {
            bytes: len.saturating_mul(size_of::<T>()),
        })
}

/// A vector of `len` copies of `value`.
///
/// # Errors
///
/// Those of [`room`].
pub(crate) fn filled<T: Clone>(len: usize, value: T) -> Result<Vec<T>, InterpolateError> {
    let mut vec = Vec::new();
    room(&mut vec, len)?;
    vec.resize(len, value);
    Ok(vec)
}

#[cfg(test)]
mod tests {
    use super::{filled, room};
    use crate::InterpolateError as E;

    /// More entries than any machine has memory for, and fewer than a `Vec`
    /// can hold: 2^58 numbers are 2^61 bytes.
    const TOO_MANY: usize = 1 << 58;

    #[test]
    fn what_fits_is_made_and_what_does_not_is_an_error() {
        assert_eq!(filled(3, 1.5_f64), Ok(vec![1.5, 1.5, 1.5]));
        assert_eq!(filled(0, 1.5_f64), Ok(Vec::new()));
        assert_eq!(
            filled(TOO_MANY, 0.0_f64),
            Err(E::AllocationFailed { bytes: 1 << 61 })
        );
        // Past what a `Vec` holds as well (a caller checks that limit first,
        // to name the argument): an error, and no panic.
        assert_eq!(
            filled(usize::MAX, [0.0_f64; 2]),
            Err(E::AllocationFailed { bytes: usize::MAX })
        );
    }

    #[test]
    fn room_is_counted_from_what_the_vector_holds_and_leaves_it_alone() {
        let mut vec = vec![7_u8, 8, 9];
        // Room for what is there already asks for nothing.
        let exact = vec.capacity();
        assert_eq!(room(&mut vec, 2), Ok(()));
        assert_eq!(room(&mut vec, exact), Ok(()));
        assert_eq!(vec.capacity(), exact);
        assert_eq!(room(&mut vec, 10), Ok(()));
        assert!(vec.capacity() >= 10);
        assert_eq!(vec, [7, 8, 9]);
        let before = vec.capacity();
        assert_eq!(
            room(&mut vec, TOO_MANY),
            Err(E::AllocationFailed { bytes: TOO_MANY })
        );
        assert_eq!((vec.as_slice(), vec.capacity()), (&[7, 8, 9][..], before));
    }
}
