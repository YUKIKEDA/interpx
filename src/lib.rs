//! Pure Rust interpolation covering the main classes of `scipy.interpolate`.
//!
//! This crate is being built. So far it has the core types ([`Scalar`],
//! [`Value`], the [`Interpolator`] trait, [`InterpolateError`], [`Extrapolate`],
//! [`BoundaryCondition`]) and the two spline representations the 1D
//! interpolators evaluate through: [`PPoly`], a piecewise polynomial, and
//! [`BSpline`]. Both can be differentiated and integrated. The first
//! interpolators built from data are [`Linear1D`], [`Nearest1D`],
//! [`CubicSpline`], [`PchipInterpolator`], [`Akima1D`],
//! [`CubicHermiteSpline`], and [`BSpline1D`] for data in one variable, and
//! [`BSplineCurve`] for a curve through points, and
//! [`RegularGridInterpolator`] (nearest, linear, cubic, and quintic) for values
//! on a grid in several dimensions; the others follow.
//!
//! # Examples
//!
//! A piecewise polynomial, evaluated at one point and at a batch:
//!
//! ```
//! use interpx::{InterpolateError, Interpolator, PPoly};
//!
//! # fn main() -> Result<(), InterpolateError> {
//! // x on [0, 1], then 1 + 2(x - 1) on [1, 2]. Coefficients are per interval,
//! // the highest power first.
//! let p = PPoly::new(&[0.0, 1.0, 2.0], &[1.0, 0.0, 2.0, 1.0], 1)?;
//! assert_eq!(p.eval(&1.5), Some(2.0));
//!
//! let mut out = [0.0; 3];
//! p.eval_batch_into(&[0.5, 1.0, f64::NAN], &mut out, -1.0)?;
//! assert_eq!(out, [0.5, 1.0, -1.0]);
//! # Ok(())
//! # }
//! ```
//!
//! Any type can be an interpolator by implementing [`Interpolator::eval`]. The
//! batch methods come with it.
//!
//! ```
//! use interpx::Interpolator;
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
//! assert_eq!(Ramp.eval_batch(&[0.5, 2.0]), [Some(5.0), None]);
//! ```

mod banded;
mod boundary;
mod bspline;
// The golden loader of the integration tests, for the unit tests of the
// crate-private modules (an integration test cannot call their functions).
// Declared here once: a file included from two modules is compiled twice.
#[cfg(test)]
#[path = "../tests/common/mod.rs"]
mod common;
mod curve;
mod error;
mod extrapolate;
mod grid;
mod interpolator;
mod oned;
mod ppoly;
mod scalar;
mod search;
mod tridiagonal;
mod validate;
mod value;

pub use boundary::{BoundaryCondition, EndCondition};
pub use bspline::{BSpline, BSPLINE_MAX_DEGREE};
pub use curve::{BSplineCurve, Parametrization};
pub use error::{AxisFault, InterpolateError, ParamReason};
pub use extrapolate::Extrapolate;
pub use grid::{GridMethod, RegularGridInterpolator};
pub use interpolator::Interpolator;
pub use oned::{
    Akima1D, AkimaMethod, BSpline1D, CubicHermiteSpline, CubicSpline, Linear1D, Nearest1D,
    NearestKind, PchipInterpolator,
};
pub use ppoly::PPoly;
pub use scalar::Scalar;
pub use search::Hint;
pub use value::Value;
