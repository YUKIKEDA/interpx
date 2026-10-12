//! Interpolators of values at points without structure, in any number of
//! dimensions. Built with the `scattered` feature.

mod nearest_nd;
#[cfg_attr(
    not(test),
    expect(dead_code, reason = "first used by RBFInterpolator (M3-3)")
)]
mod rbf;

pub use nearest_nd::NearestNDInterpolator;
