//! Interpolators of values at points without structure, in any number of
//! dimensions. Built with the `scattered` feature.

mod nearest_nd;

pub use nearest_nd::NearestNDInterpolator;
