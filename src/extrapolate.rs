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
    /// Wrap the query into the domain, as for a periodic function.
    Periodic,
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
