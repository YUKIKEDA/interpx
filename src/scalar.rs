//! The floating-point types interpolators are generic over.

use core::fmt::Debug;

mod sealed {
    pub trait Sealed {}
    impl Sealed for f32 {}
    impl Sealed for f64 {}
}

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
