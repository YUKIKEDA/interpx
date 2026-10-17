//! Interpolators of values at points without structure, in any number of
//! dimensions. Built with the `scattered` feature.

mod kdtree;
mod linear_nd;
mod nearest_nd;
mod rbf;

pub use linear_nd::LinearNDInterpolator;
pub use nearest_nd::NearestNDInterpolator;
pub use rbf::{RBFInterpolator, RbfDegree, RbfKernel, RbfOptions};

use crate::InterpolateError;
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

/// The change of coordinates of a rescaled interpolant: `(x - offset) / scale`
/// on each axis.
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(bound = ""))]
pub(super) struct Scaling<const D: usize> {
    #[cfg_attr(feature = "serde", serde(with = "crate::stored"))]
    offset: [f64; D],
    #[cfg_attr(feature = "serde", serde(with = "crate::stored"))]
    scale: [f64; D],
}

impl<const D: usize> Scaling<D> {
    /// The mean of the points on each axis, and their range. An axis on which
    /// all points agree keeps its scale.
    ///
    /// # Errors
    ///
    /// Returns [`InterpolateError::NonFinite`] when the sum of the coordinates
    /// of an axis overflows (with the point at which it does), or their range
    /// (with the point that has the largest coordinate). With a mean or a
    /// range that is not finite every point would land on one place, or on
    /// none.
    pub(super) fn of(points: &[[f64; D]]) -> Result<Self, InterpolateError> {
        let count = points.len() as f64;
        let mut offset = [0.0; D];
        let mut scale = [1.0; D];
        for d in 0..D {
            let (mut sum, mut least, mut most) = (0.0, f64::INFINITY, f64::NEG_INFINITY);
            let mut largest = 0;
            for (index, point) in points.iter().enumerate() {
                sum += point[d];
                if !sum.is_finite() {
                    return Err(InterpolateError::NonFinite { index });
                }
                least = least.min(point[d]);
                if point[d] > most {
                    (most, largest) = (point[d], index);
                }
            }
            offset[d] = sum / count;
            let range = most - least;
            if !range.is_finite() {
                return Err(InterpolateError::NonFinite { index: largest });
            }
            if range > 0.0 {
                scale[d] = range;
            }
        }
        Ok(Self { offset, scale })
    }

    pub(super) fn apply(&self, point: &[f64; D]) -> [f64; D] {
        core::array::from_fn(|d| (point[d] - self.offset[d]) / self.scale[d])
    }
}
