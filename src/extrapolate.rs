//! What an interpolator returns for a query outside its domain.

/// The behavior of a 1D, curve, or grid interpolator outside its domain.
///
/// There is no one default. Each interpolator starts with the mode SciPy uses for
/// its counterpart: `numpy.interp` clamps, `CubicSpline` extends, and
/// `Akima1DInterpolator` forbids, for example.
///
/// Scattered-data interpolators have no extrapolation mode: outside the convex
/// hull they return `None`.
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
#[non_exhaustive]
pub enum Extrapolate {
    /// Return `None` outside the domain.
    Forbid,
    /// Return the value at the nearest end of the domain.
    Clamp,
    /// Continue the first or the last polynomial piece past the end.
    Extend,
    /// Wrap the query into the domain, as for a periodic function. The right end
    /// of the domain wraps to the left end.
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
