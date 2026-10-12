//! The types an interpolator can interpolate: a scalar or a fixed-length array.

use crate::Scalar;
use core::fmt::Debug;

pub(crate) mod sealed {
    use crate::Scalar;

    /// Seals [`Value`](super::Value), and holds what the crate needs from a
    /// value and callers do not.
    pub trait Sealed<T> {
        /// Whether every column of `self` is within `tolerance * (1 + |b|)` of
        /// the same column `b` of `other`.
        fn within(&self, other: &Self, tolerance: T) -> bool;

        /// The value whose every column is `f` of the same column of `self`
        /// and of `other`: what a rule that is not linear in the data needs
        /// (a comparison of signs, a harmonic mean).
        #[must_use]
        fn zip_with(self, other: Self, f: impl Fn(T, T) -> T) -> Self;
    }

    impl<T: Scalar> Sealed<T> for T {
        fn within(&self, other: &Self, tolerance: T) -> bool {
            (*self - *other).abs() <= tolerance * (T::one() + other.abs())
        }

        fn zip_with(self, other: Self, f: impl Fn(T, T) -> T) -> Self {
            f(self, other)
        }
    }

    impl<T: Scalar, const M: usize> Sealed<T> for [T; M] {
        fn within(&self, other: &Self, tolerance: T) -> bool {
            self.iter().zip(other).all(|(a, b)| a.within(b, tolerance))
        }

        fn zip_with(mut self, other: Self, f: impl Fn(T, T) -> T) -> Self {
            for (a, b) in self.iter_mut().zip(other) {
                *a = f(*a, b);
            }
            self
        }
    }
}

/// A value an interpolator returns: a scalar `T`, or an array `[T; M]` of them.
///
/// SciPy accepts `values` with several columns. Here the number of columns is
/// part of the type: an interpolator of `[f64; 3]` values returns `[f64; 3]`.
/// The trait is sealed and is implemented for `T` and `[T; M]` only.
///
/// The operations are the ones a linear combination of data needs.
///
/// # Examples
///
/// ```
/// use interpx::Value;
///
/// // 1 + 0.5 * 4
/// assert_eq!(1.0_f64.axpy(0.5, 4.0), 3.0);
///
/// // The same, column by column.
/// let v = [1.0_f64, 10.0].axpy(0.5, [4.0, 40.0]);
/// assert_eq!(v, [3.0, 30.0]);
///
/// assert_eq!(<[f64; 2] as Value<f64>>::zero(), [0.0, 0.0]);
/// assert!(!Value::<f64>::all_finite(&[1.0, f64::NAN]));
/// ```
pub trait Value<T: Scalar>: Copy + Debug + Send + Sync + sealed::Sealed<T> {
    /// Returns the additive identity: zero in every column.
    fn zero() -> Self;

    /// Returns `self + k * x`, column by column.
    #[must_use]
    fn axpy(self, k: T, x: Self) -> Self;

    /// Returns `true` when every column is neither NaN nor infinite.
    fn all_finite(&self) -> bool;
}

impl<T: Scalar> Value<T> for T {
    fn zero() -> Self {
        T::zero()
    }

    fn axpy(self, k: T, x: Self) -> Self {
        self + k * x
    }

    fn all_finite(&self) -> bool {
        self.is_finite()
    }
}

impl<T: Scalar, const M: usize> Value<T> for [T; M] {
    fn zero() -> Self {
        [T::zero(); M]
    }

    fn axpy(mut self, k: T, x: Self) -> Self {
        for (a, b) in self.iter_mut().zip(x) {
            *a = *a + k * b;
        }
        self
    }

    fn all_finite(&self) -> bool {
        self.iter().all(|c| c.is_finite())
    }
}

#[cfg(test)]
mod tests {
    use super::Value;

    #[test]
    fn scalar_values_combine_linearly() {
        assert_eq!(<f64 as Value<f64>>::zero().to_bits(), 0.0_f64.to_bits());
        assert_eq!(2.0_f64.axpy(3.0, 4.0).to_bits(), 14.0_f64.to_bits());
        assert_eq!(2.0_f32.axpy(3.0, 4.0).to_bits(), 14.0_f32.to_bits());
        assert!(1.0_f64.all_finite());
        assert!(!f64::NAN.all_finite());
        assert!(!f32::INFINITY.all_finite());
        assert!(!f64::NEG_INFINITY.all_finite());
    }

    #[test]
    fn array_values_combine_column_by_column() {
        let zero = <[f32; 3] as Value<f32>>::zero();
        assert!(zero.iter().all(|c| c.to_bits() == 0.0_f32.to_bits()));

        let v = [1.0_f64, 2.0, 3.0].axpy(2.0, [10.0, 20.0, 30.0]);
        let want = [21.0_f64, 42.0, 63.0];
        assert!(v.iter().zip(want).all(|(a, b)| a.to_bits() == b.to_bits()));

        assert!(Value::<f64>::all_finite(&[1.0, 2.0]));
        assert!(!Value::<f64>::all_finite(&[1.0, f64::NAN]));
        assert!(!Value::<f64>::all_finite(&[f64::INFINITY, 1.0]));
    }

    #[test]
    fn an_empty_array_is_a_value() {
        // M = 0 is allowed by the type; it has no columns to be non-finite.
        let v = <[f64; 0] as Value<f64>>::zero();
        assert!(Value::<f64>::all_finite(&v.axpy(1.0, [])));
    }
}
