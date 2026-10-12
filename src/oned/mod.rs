//! Interpolators of data in one variable.

mod cubic_spline;
mod hermite;
mod linear;
mod nearest;
mod pchip;

pub use cubic_spline::CubicSpline;
pub use linear::Linear1D;
pub use nearest::{Nearest1D, NearestKind};
pub use pchip::PchipInterpolator;
