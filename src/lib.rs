//! Pure Rust interpolation covering the main classes of `scipy.interpolate`.
//!
//! This crate is being built. The core types are in place: [`Scalar`], [`Value`],
//! the [`Interpolator`] trait, [`InterpolateError`], [`Extrapolate`], and
//! [`BoundaryCondition`]. The interpolators themselves follow.
//!
//! # Examples
//!
//! Any type can be an interpolator by implementing [`Interpolator::eval`]. The
//! batch methods come with it.
//!
//! ```
//! use interpx::{InterpolateError, Interpolator};
//!
//! /// Linear interpolation between (0, 0) and (1, 10).
//! struct Ramp;
//!
//! impl Interpolator for Ramp {
//!     type Point = f64;
//!     type Output = f64;
//!
//!     fn eval(&self, x: &f64) -> Option<f64> {
//!         (0.0..=1.0).contains(x).then(|| 10.0 * x)
//!     }
//! }
//!
//! # fn main() -> Result<(), InterpolateError> {
//! let mut out = [0.0; 3];
//! Ramp.eval_batch_into(&[0.5, 2.0, f64::NAN], &mut out, -1.0)?;
//! assert_eq!(out, [5.0, -1.0, -1.0]);
//! # Ok(())
//! # }
//! ```

mod boundary;
mod error;
mod extrapolate;
mod interpolator;
mod scalar;
mod value;

pub use boundary::{BoundaryCondition, EndCondition};
pub use error::InterpolateError;
pub use extrapolate::Extrapolate;
pub use interpolator::Interpolator;
pub use scalar::Scalar;
pub use value::Value;
