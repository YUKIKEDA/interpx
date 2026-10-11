//! Interpolators of data in one variable.

mod cubic_spline;
mod linear;
mod nearest;

pub use cubic_spline::CubicSpline;
pub use linear::Linear1D;
pub use nearest::{Nearest1D, NearestKind};
