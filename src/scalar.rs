//! The floating-point types interpolators are generic over.

use core::fmt::Debug;

mod sealed {
    pub trait Sealed: super::Kept {}
    impl Sealed for f32 {}
    impl Sealed for f64 {}
}

/// What the `serde` feature needs of a scalar and of a value: that it can be a
/// part of a stored form. Without the feature, nothing.
///
/// It is a supertrait of the sealed traits and not of [`Scalar`] or `Value`,
/// so those read the same with and without the feature.
#[cfg(feature = "serde")]
pub trait Kept: crate::stored::Stored {}
#[cfg(feature = "serde")]
impl<X: crate::stored::Stored> Kept for X {}

/// What the `serde` feature needs of a scalar and of a value. Without the
/// feature, nothing.
#[cfg(not(feature = "serde"))]
pub trait Kept {}
#[cfg(not(feature = "serde"))]
impl<X> Kept for X {}

/// A floating-point scalar: `f32` or `f64`.
///
/// This is the one bound on the scalar type parameter of every interpolator. It
/// is sealed: no other type can implement it. Because it includes `Send + Sync`,
/// an interpolator built from scalars is `Send + Sync` as well.
///
/// # Examples
///
/// ```
/// use interpx::Scalar;
///
/// fn midpoint<T: Scalar>(a: T, b: T) -> T {
///     (a + b) / (T::one() + T::one())
/// }
///
/// assert_eq!(midpoint(1.0_f32, 3.0), 2.0);
/// assert_eq!(midpoint(1.0_f64, 3.0), 2.0);
/// ```
pub trait Scalar: num_traits::Float + Debug + Send + Sync + sealed::Sealed {}

impl Scalar for f32 {}
impl Scalar for f64 {}

/// The integer `n` as a scalar, by counting: exact while `n` is small enough for
/// the type to hold every integer up to it, as a degree or an order is.
pub(crate) fn count<T: Scalar>(n: usize) -> T {
    (0..n).fold(T::zero(), |sum, _| sum + T::one())
}

/// `m (m - 1) ... (m - nu + 1)`: what differentiating `x^m` `nu` times puts in
/// front of `x^(m - nu)`.
pub(crate) fn falling<T: Scalar>(m: T, nu: usize) -> T {
    let mut product = T::one();
    let mut factor = m;
    for _ in 0..nu {
        product = product * factor;
        factor = factor - T::one();
    }
    product
}
