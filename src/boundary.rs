//! Boundary conditions of a spline.

use crate::{Scalar, Value};

/// The condition at one end of a spline.
///
/// `V` is the type of the values the spline interpolates: a derivative at an
/// end has a value of that type, so one number per column for an array.
///
/// # Examples
///
/// ```
/// use interpx::EndCondition;
///
/// // The slope at this end is 0.5.
/// let end = EndCondition::FirstDerivative(0.5_f64);
/// assert_ne!(end, EndCondition::NotAKnot);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(bound = "V: crate::stored::Stored"))]
pub enum EndCondition<V> {
    /// The first two pieces at this end are one polynomial.
    NotAKnot,
    /// The first derivative at this end has the given value.
    FirstDerivative(#[cfg_attr(feature = "serde", serde(with = "crate::stored"))] V),
    /// The second derivative at this end has the given value.
    SecondDerivative(#[cfg_attr(feature = "serde", serde(with = "crate::stored"))] V),
}

/// The boundary conditions of a spline.
///
/// SciPy's `bc_type`. The default is [`BoundaryCondition::NotAKnot`], as in SciPy.
/// A periodic spline is periodic at both ends, so there is no way to ask for one
/// periodic end.
///
/// `V` is the type of the values the spline interpolates, as in
/// [`EndCondition`].
///
/// # Examples
///
/// ```
/// use interpx::{BoundaryCondition, EndCondition};
///
/// // SciPy's bc_type="natural": zero second derivative at both ends.
/// let natural = BoundaryCondition::<f64>::natural();
/// assert_eq!(
///     natural,
///     BoundaryCondition::Ends(
///         EndCondition::SecondDerivative(0.0),
///         EndCondition::SecondDerivative(0.0),
///     ),
/// );
///
/// // SciPy's bc_type=((1, 0.5), (2, 0.0)): a slope on the left, no curvature on the right.
/// let mixed = BoundaryCondition::Ends(
///     EndCondition::FirstDerivative(0.5),
///     EndCondition::SecondDerivative(0.0),
/// );
/// assert_ne!(mixed, natural);
///
/// assert_eq!(BoundaryCondition::<f64>::default(), BoundaryCondition::NotAKnot);
///
/// // For a spline of two columns, a slope for each.
/// let slopes = BoundaryCondition::Ends(
///     EndCondition::FirstDerivative([1.0, 0.0]),
///     EndCondition::NotAKnot,
/// );
/// assert_ne!(slopes, BoundaryCondition::<[f64; 2]>::clamped());
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(bound = "V: crate::stored::Stored"))]
pub enum BoundaryCondition<V> {
    /// Not-a-knot at both ends.
    #[default]
    NotAKnot,
    /// The spline and its derivatives match at the two ends.
    Periodic,
    /// A condition for the left end, then one for the right end.
    Ends(EndCondition<V>, EndCondition<V>),
}

impl<V> BoundaryCondition<V> {
    /// Returns the natural boundary: zero second derivative at both ends.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{BoundaryCondition, EndCondition};
    ///
    /// let end = EndCondition::SecondDerivative(0.0_f32);
    /// assert_eq!(BoundaryCondition::natural(), BoundaryCondition::Ends(end, end));
    /// ```
    #[must_use]
    pub fn natural<T: Scalar>() -> Self
    where
        V: Value<T>,
    {
        let end = EndCondition::SecondDerivative(V::zero());
        Self::Ends(end, end)
    }

    /// Returns the clamped boundary: zero first derivative at both ends.
    ///
    /// # Examples
    ///
    /// ```
    /// use interpx::{BoundaryCondition, EndCondition};
    ///
    /// let end = EndCondition::FirstDerivative(0.0_f64);
    /// assert_eq!(BoundaryCondition::clamped(), BoundaryCondition::Ends(end, end));
    /// ```
    #[must_use]
    pub fn clamped<T: Scalar>() -> Self
    where
        V: Value<T>,
    {
        let end = EndCondition::FirstDerivative(V::zero());
        Self::Ends(end, end)
    }
}

#[cfg(test)]
mod tests {
    use super::{BoundaryCondition as Bc, EndCondition as End};

    #[test]
    fn the_default_is_not_a_knot_for_both_scalars() {
        assert_eq!(Bc::<f32>::default(), Bc::NotAKnot);
        assert_eq!(Bc::<f64>::default(), Bc::NotAKnot);
    }

    #[test]
    fn natural_and_clamped_set_the_named_derivative_to_zero() {
        assert_eq!(
            Bc::<f64>::natural(),
            Bc::Ends(End::SecondDerivative(0.0), End::SecondDerivative(0.0))
        );
        assert_eq!(
            Bc::<f32>::clamped(),
            Bc::Ends(End::FirstDerivative(0.0), End::FirstDerivative(0.0))
        );
        assert_ne!(Bc::<f64>::natural(), Bc::clamped());
    }

    #[test]
    fn the_two_ends_are_independent() {
        let mixed = Bc::Ends(End::FirstDerivative(0.5_f64), End::NotAKnot);
        let swapped = Bc::Ends(End::NotAKnot, End::FirstDerivative(0.5_f64));
        assert_ne!(mixed, swapped);
        assert_ne!(mixed, Bc::Periodic);
        assert_ne!(
            End::FirstDerivative(1.0_f64),
            End::SecondDerivative(1.0_f64)
        );
    }
}
