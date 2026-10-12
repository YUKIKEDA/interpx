//! Interpolators of values at points without structure, in any number of
//! dimensions. Built with the `scattered` feature.

mod kdtree;
mod nearest_nd;
mod rbf;

pub use nearest_nd::NearestNDInterpolator;
pub use rbf::{RBFInterpolator, RbfDegree, RbfKernel, RbfOptions};

use std::collections::HashMap;

/// The coordinates of a point as a key: two points have the same key when
/// their coordinates are equal as numbers (`0.0` and `-0.0` are one
/// coordinate; `-0.0 + 0.0` is `0.0`).
fn key<const D: usize>(point: &[f64; D]) -> [u64; D] {
    point.map(|c| (c + 0.0).to_bits())
}

/// The first point that has the coordinates of an earlier one, with that
/// earlier one first.
fn first_repeat<const D: usize>(points: &[[f64; D]]) -> Option<(usize, usize)> {
    let mut seen = HashMap::with_capacity(points.len());
    points
        .iter()
        .enumerate()
        .find_map(|(index, point)| Some((seen.insert(key(point), index)?, index)))
}
