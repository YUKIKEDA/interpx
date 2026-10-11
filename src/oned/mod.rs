//! Interpolators of data in one variable.

mod linear;
mod nearest;

pub use linear::Linear1D;
pub use nearest::{Nearest1D, NearestKind};
