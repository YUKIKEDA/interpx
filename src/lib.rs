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
//! # Features
//!
//! - `scattered` (on by default): the interpolators of scattered data, and
//!   the dependencies only they need.
//! - `rayon`: batches evaluated on several threads, through
//!   `ParInterpolator`, which every interpolator that is `Sync` then has. The
//!   results are those of the serial batches bit for bit. A batch times its
//!   first points and is shared out only when the rest is worth it; for
//!   batches you know to be small, call the serial forms: a batch of less
//!   than a tenth of a millisecond of work is not shared out, and deciding
//!   that has made it up to a third slower than the serial form.
//! - `serde`: `Serialize` and `Deserialize` for every interpolator and for
//!   the parameters they are built with. See below.
//!
//! # Storing an interpolator
//!
//! With the `serde` feature an interpolator is written by any format of
//! `serde`, as text or as bytes, and read back without being built again:
//! what is stored is what it was built into (the coefficients of a spline,
//! not the data), with a version of the form.
//!
//! Reading runs the checks of the constructors, so a form that is broken or
//! was written by hand is an error of the format and never an interpolator
//! that panics or gives values of nothing. A version that this build does
//! not know is an error as well.
//!
//! What comes back evaluates bit for bit as what was stored, when the format
//! gives every number back as it was. Binary formats do. `serde_json` does
//! with its `float_roundtrip` feature; without it a number can come back
//! one bit off. A stored form holds finite numbers only.
//!
#![cfg_attr(feature = "serde", doc = "```")]
#![cfg_attr(not(feature = "serde"), doc = "```ignore")]
//! use interpx::{BoundaryCondition, CubicSpline, Interpolator};
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let x = [0.0, 1.0, 2.5, 4.0];
//! let y = [1.0, 3.0, 2.0, 5.0];
//! let spline = CubicSpline::new(&x, &y, BoundaryCondition::NotAKnot)?;
//!
//! let text = serde_json::to_string(&spline)?;
//! let back: CubicSpline<f64> = serde_json::from_str(&text)?;
//! assert_eq!(back.eval(&1.7), spline.eval(&1.7));
//!
//! // Coordinates out of order: not a spline that could have been built.
//! let broken = text.replace("[0.0,1.0,2.5,4.0]", "[1.0,0.0,2.5,4.0]");
//! assert!(serde_json::from_str::<CubicSpline<f64>>(&broken).is_err());
//! # Ok(())
//! # }
//! ```
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
#[cfg(feature = "rayon")]
mod par;
mod ppoly;
mod reserve;
mod scalar;
#[cfg(feature = "scattered")]
mod scattered;
mod search;
#[cfg(feature = "serde")]
mod stored;
mod tridiagonal;
mod validate;
mod value;

pub use boundary::{BoundaryCondition, EndCondition};
pub use bspline::{BSpline, BSPLINE_MAX_DEGREE};
pub use curve::{BSplineCurve, Parametrization};
pub use error::{AxisFault, InterpolateError, ParamReason, SolverError};
pub use extrapolate::Extrapolate;
pub use grid::{GridMethod, RegularGridInterpolator};
pub use interpolator::Interpolator;
pub use oned::{
    Akima1D, AkimaMethod, BSpline1D, CubicHermiteSpline, CubicSpline, Linear1D, Nearest1D,
    NearestKind, PchipInterpolator,
};
#[cfg(feature = "rayon")]
pub use par::ParInterpolator;
pub use ppoly::PPoly;
pub use scalar::Scalar;
#[cfg(feature = "scattered")]
pub use scattered::{NearestNDInterpolator, RBFInterpolator, RbfDegree, RbfKernel, RbfOptions};
pub use search::Hint;
pub use value::Value;
