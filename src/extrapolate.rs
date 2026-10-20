//! What an interpolator returns for a query outside its domain.

/// The behavior of a 1D, curve, or grid interpolator outside its domain.
///
/// There is no one default. Each interpolator starts with the mode that suits
/// its method, and its documentation names it: `Linear1D` clamps, `CubicSpline`
/// extends, and `Akima1D` forbids, for example.
///
/// Scattered-data interpolators have no extrapolation mode. `LinearNDInterpolator`
/// and `CloughTocher2DInterpolator` have no value outside the convex hull of the
/// points; `NearestNDInterpolator` has one at every finite query point unless it
/// is given a distance bound; `RBFInterpolator` has one wherever its kernel and
/// polynomial give a finite number.
///
/// # Examples
///
/// ```
/// use interpx::Extrapolate;
///
/// let mode = Extrapolate::Clamp;
/// assert_ne!(mode, Extrapolate::Forbid);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[non_exhaustive]
pub enum Extrapolate {
    /// Gives no value (`None`) outside the domain.
    Forbid,
    /// Gives the value at the nearest end of the domain.
    Clamp,
    /// Continues the first or the last polynomial piece past the end.
    Extend,
    /// Wraps the query into the domain, as for a periodic function. The right
    /// end of the domain wraps to the left end.
    Periodic,
}

impl Extrapolate {
    /// The mode a derivative takes from a function with this mode.
    ///
    /// A clamped function is constant outside its domain, so its derivative is
    /// zero there. No mode says that, so the derivative has no value outside.
    pub(crate) fn of_derivative(self) -> Self {
        match self {
            Self::Clamp => Self::Forbid,
            mode => mode,
        }
    }

    /// The mode an antiderivative takes from a function with this mode.
    ///
    /// The antiderivative of a periodic function is not periodic unless its mean
    /// is zero, and that of a clamped function is linear outside the domain. No
    /// mode says either, so the antiderivative has no value outside.
    pub(crate) fn of_antiderivative(self) -> Self {
        match self {
            Self::Extend => Self::Extend,
            _ => Self::Forbid,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::Extrapolate;
    use std::collections::HashSet;

    #[test]
    fn the_four_modes_are_distinct() {
        let modes = [
            Extrapolate::Forbid,
            Extrapolate::Clamp,
            Extrapolate::Extend,
            Extrapolate::Periodic,
        ];
        let set: HashSet<Extrapolate> = modes.into_iter().collect();
        assert_eq!(set.len(), 4);
        assert_eq!(format!("{:?}", Extrapolate::Forbid), "Forbid");
    }
}
