//! Interpolators of data in one variable.

mod akima;
mod cubic_hermite;
mod cubic_spline;
mod hermite;
mod linear;
mod nearest;
mod pchip;

pub use akima::{Akima1D, AkimaMethod};
pub use cubic_hermite::CubicHermiteSpline;
pub use cubic_spline::CubicSpline;
pub use linear::Linear1D;
pub use nearest::{Nearest1D, NearestKind};
pub use pchip::PchipInterpolator;
