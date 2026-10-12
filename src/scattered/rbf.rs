//! The parts of the radial basis function interpolant, up to its system of
//! equations: the kernels, the polynomial terms, and the matrix they make.
//!
//! The interpolant through the points `y_i` is a sum of one kernel centered at
//! every point, plus a polynomial of low degree:
//!
//! ```text
//! f(x) = sum over i of w_i phi(epsilon |x - y_i|) + sum over j of c_j p_j(x^)
//! ```
//!
//! The weights `w` and the coefficients `c` solve one symmetric system
//! ([`system`]); the value at a point is the inner product of the solution
//! with the row of that point ([`evaluation_row`]). The formulas, the signs,
//! and the order of every array are those of SciPy 1.18's `RBFInterpolator`.

use crate::error::ParamReason;
use crate::{InterpolateError, Value};

/// A radial basis function: a function of the distance from a center, the
/// distance multiplied by `epsilon` where the kernel has one.
///
/// The four kernels without `epsilon` do not depend on a scale: stretching the
/// distances changes the weights and not the interpolant.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum Kernel {
    /// `-r`
    Linear,
    /// `r^2 log r`, and 0 at `r = 0`
    ThinPlateSpline,
    /// `r^3`
    Cubic,
    /// `-r^5`
    Quintic,
    /// `-sqrt(r^2 + 1)`
    Multiquadric { epsilon: f64 },
    /// `1 / sqrt(r^2 + 1)`
    InverseMultiquadric { epsilon: f64 },
    /// `1 / (r^2 + 1)`
    InverseQuadratic { epsilon: f64 },
    /// `exp(-r^2)`
    Gaussian { epsilon: f64 },
}

impl Kernel {
    /// What the distances are multiplied by: 1 for a kernel that does not
    /// depend on a scale.
    fn epsilon(self) -> f64 {
        match self {
            Self::Linear | Self::ThinPlateSpline | Self::Cubic | Self::Quintic => 1.0,
            Self::Multiquadric { epsilon }
            | Self::InverseMultiquadric { epsilon }
            | Self::InverseQuadratic { epsilon }
            | Self::Gaussian { epsilon } => epsilon,
        }
    }

    /// The smallest degree of the polynomial the kernel needs for its system
    /// to have one solution, where it needs one: a kernel that is
    /// conditionally positive definite of order `m` needs the degree `m - 1`
    /// (Fasshauer, Meshfree Approximation Methods with MATLAB, chapter 8).
    fn min_degree(self) -> Option<usize> {
        match self {
            Self::Linear | Self::Multiquadric { .. } => Some(0),
            Self::ThinPlateSpline | Self::Cubic => Some(1),
            Self::Quintic => Some(2),
            Self::InverseMultiquadric { .. }
            | Self::InverseQuadratic { .. }
            | Self::Gaussian { .. } => None,
        }
    }

    /// The degree used when none is asked for: the smallest one, and a
    /// constant term for the kernels that need no polynomial (SciPy's default).
    pub(crate) fn default_degree(self) -> usize {
        self.min_degree().unwrap_or(0)
    }

    /// The value at the scaled distance `r`.
    fn at(self, r: f64) -> f64 {
        match self {
            Self::Linear => -r,
            Self::ThinPlateSpline if r > 0.0 => r * r * r.ln(),
            // The limit at 0, where the logarithm is not finite.
            Self::ThinPlateSpline => 0.0,
            Self::Cubic => r * r * r,
            Self::Quintic => -(r * r * r * r * r),
            Self::Multiquadric { .. } => -(r * r + 1.0).sqrt(),
            Self::InverseMultiquadric { .. } => (r * r + 1.0).sqrt().recip(),
            Self::InverseQuadratic { .. } => (r * r + 1.0).recip(),
            Self::Gaussian { .. } => (-(r * r)).exp(),
        }
    }

    /// The value for two points: the coordinates are scaled first and their
    /// distance taken then, as SciPy does.
    fn between<const D: usize>(self, a: &[f64; D], b: &[f64; D]) -> f64 {
        let epsilon = self.epsilon();
        let squared: f64 = a
            .iter()
            .zip(b)
            .map(|(a, b)| {
                let apart = a * epsilon - b * epsilon;
                apart * apart
            })
            .sum();
        self.at(squared.sqrt())
    }
}

/// The number of monomials of degree `degree` or less in `variables`
/// variables, `C(variables + degree, variables)`, when it fits.
fn monomial_count(variables: usize, degree: usize) -> Option<usize> {
    // C(degree + i, i) from C(degree + i - 1, i - 1): times `degree + i`,
    // over `i`. The division is exact at every step.
    let mut count: u128 = 1;
    for i in 1..=variables as u128 {
        count = count.checked_mul(degree as u128 + i)? / i;
    }
    usize::try_from(count).ok()
}

/// The polynomial terms of an interpolant: the monomials, and the change of
/// coordinates they are evaluated in.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Polynomial<const D: usize> {
    /// The power of each variable in each monomial, in SciPy's order: by
    /// degree, and among those of one degree by the variables they multiply,
    /// taken in increasing order (`1, x, y, x^2, x y, y^2`).
    powers: Vec<[i32; D]>,
    /// The middle of the points on each axis.
    shift: [f64; D],
    /// Half the range of the points on each axis; 1 where they all agree.
    scale: [f64; D],
}

impl<const D: usize> Polynomial<D> {
    /// The terms for `kernel` on `points`: of `degree` or less, or none.
    ///
    /// The coordinates of the points go to `[-1, 1]` on each axis before a
    /// monomial is evaluated, so that the powers stay of a size.
    ///
    /// # Errors
    ///
    /// The first of these that applies:
    ///
    /// - [`InterpolateError::InvalidParameter`] for `degree` with
    ///   [`ParamReason::TooSmall`] when it is under the smallest degree the
    ///   kernel needs. No terms at all is accepted for every kernel, as in
    ///   SciPy.
    /// - [`InterpolateError::InvalidParameter`] for `degree` with
    ///   [`ParamReason::TooLarge`] when the monomials are too many to count.
    /// - [`InterpolateError::InsufficientPoints`] when there are fewer points
    ///   than monomials: the terms could not be told apart on them.
    pub(crate) fn new(
        points: &[[f64; D]],
        kernel: Kernel,
        degree: Option<usize>,
    ) -> Result<Self, InterpolateError> {
        let invalid = |reason| InterpolateError::InvalidParameter {
            name: "degree",
            reason,
        };
        if let (Some(degree), Some(min)) = (degree, kernel.min_degree()) {
            if degree < min {
                return Err(invalid(ParamReason::TooSmall { min }));
            }
        }
        // A power is an `i32` (what `powi` takes), and the count a `usize`.
        let fits =
            |degree: usize| i32::try_from(degree).is_ok() && monomial_count(D, degree).is_some();
        let count = match degree {
            None => 0,
            Some(degree) if fits(degree) => monomial_count(D, degree).unwrap_or(usize::MAX),
            Some(degree) => {
                // The largest degree that fits is between 0, which does, and
                // this one, which does not.
                let (mut fitting, mut not) = (0, degree);
                while not - fitting > 1 {
                    let middle = fitting + (not - fitting) / 2;
                    if fits(middle) {
                        fitting = middle;
                    } else {
                        not = middle;
                    }
                }
                return Err(invalid(ParamReason::TooLarge { max: fitting }));
            }
        };
        if points.len() < count {
            return Err(InterpolateError::InsufficientPoints {
                required: count,
                actual: points.len(),
            });
        }

        let mut powers = Vec::with_capacity(count);
        // Without variables the constant is the only monomial.
        let degrees = degree.map(|degree| 0..=if D == 0 { 0 } else { degree });
        for total in degrees.into_iter().flatten() {
            // The variables a monomial multiplies, with repeats, in
            // increasing order: `[0, 1, 1]` is `x y^2`. All such lists of
            // `total` variables, in the order of a dictionary.
            let mut variables = vec![0; total];
            loop {
                let mut power = [0; D];
                for variable in &variables {
                    power[*variable] += 1;
                }
                powers.push(power);
                // The next list: the last entry that can still grow grows,
                // and the entries after it start again from its value.
                let Some(at) = variables.iter().rposition(|v| v + 1 < D) else {
                    break;
                };
                let next = variables[at] + 1;
                variables[at..].fill(next);
            }
        }

        let mut shift = [0.0; D];
        let mut scale = [1.0; D];
        for d in 0..D {
            let (mut least, mut most) = (f64::INFINITY, f64::NEG_INFINITY);
            for point in points {
                least = least.min(point[d]);
                most = most.max(point[d]);
            }
            if !points.is_empty() {
                shift[d] = (most + least) / 2.0;
                let half = (most - least) / 2.0;
                if half > 0.0 {
                    scale[d] = half;
                }
            }
        }
        Ok(Self {
            powers,
            shift,
            scale,
        })
    }

    /// The number of monomials.
    pub(crate) fn len(&self) -> usize {
        self.powers.len()
    }

    /// Writes the value of every monomial at `point` to `out`, which has one
    /// entry a monomial.
    fn row(&self, point: &[f64; D], out: &mut [f64]) {
        let moved: [f64; D] = core::array::from_fn(|d| (point[d] - self.shift[d]) / self.scale[d]);
        for (value, power) in out.iter_mut().zip(&self.powers) {
            *value = moved.iter().zip(power).map(|(x, p)| x.powi(*p)).product();
        }
    }
}

/// The side of the tiles the kernels of the matrix are worked out by: two
/// tiles of this side fit the first cache.
const TILE: usize = 64;

/// The system of equations of the interpolant through `values` at `points`:
/// the matrix, row by row, and the right-hand side.
///
/// With `n` points and `r` monomials the matrix has `n + r` rows and is
/// symmetric:
///
/// ```text
/// | A    P |      A[i][j] = kernel(points[i], points[j]), plus `smoothing`
/// | P^T  0 |                on the diagonal
///                 P[i][j] = monomial j at point i
/// ```
///
/// The right-hand side is the values, and a zero for every monomial.
///
/// # Errors
///
/// Returns [`InterpolateError::InvalidParameter`] for `points` with
/// [`ParamReason::TooLarge`] when the matrix would be larger than a `Vec`
/// holds.
pub(crate) fn system<const D: usize, V: Value<f64>>(
    points: &[[f64; D]],
    values: &[V],
    smoothing: f64,
    kernel: Kernel,
    polynomial: &Polynomial<D>,
) -> Result<(Vec<f64>, Vec<V>), InterpolateError> {
    let (n, r) = (points.len(), polynomial.len());
    // A `Vec` holds `isize::MAX` bytes at most, and asking it for more
    // panics.
    let most = isize::MAX.unsigned_abs() / size_of::<f64>();
    let side = n.checked_add(r).filter(|side| {
        side.checked_mul(*side)
            .is_some_and(|entries| entries <= most)
    });
    let Some(side) = side else {
        return Err(InterpolateError::InvalidParameter {
            name: "points",
            reason: ParamReason::TooLarge {
                max: most.isqrt().saturating_sub(r),
            },
        });
    };

    let mut matrix = vec![0.0; side * side];
    // The kernel of two points does not depend on their order: each is worked
    // out once and written at both places, tile by tile, so that both places
    // are in the cache. Against working out every entry, in ns an entry:
    //
    //   points   thin plate spline   Gaussian     linear, cubic, quintic
    //   100      3.4 for 6.3         3.2 for 6.1  1.0 to 1.5 for 1.6 to 2.1
    //   1000     5.0 for 7.0         4.9 for 6.6  2.6 to 2.9 for 2.2 to 2.3
    //   3000     5.3 for 7.2         5.0 for 6.8  2.7 to 3.1 for 2.2 to 2.3
    //
    // The kernels with a logarithm or an exponential gain at every size. The
    // cheap ones lose on a large matrix, where a second place to write costs
    // more than their arithmetic. A large matrix is where this step counts
    // least: it grows with the square of the points, and the solve that
    // follows with the cube. Filling the lower half in a pass of its own was
    // slower still for them (3.1 to 3.5 at 3000 points).
    for rows in (0..n).step_by(TILE) {
        for columns in (rows..n).step_by(TILE) {
            for i in rows..n.min(rows + TILE) {
                for j in columns.max(i)..n.min(columns + TILE) {
                    let value = kernel.between(&points[i], &points[j]);
                    matrix[i * side + j] = value;
                    matrix[j * side + i] = value;
                }
            }
        }
    }
    for (i, point) in points.iter().enumerate() {
        matrix[i * side + i] += smoothing;
        // The monomials at the point: a row of the upper right block, and a
        // column of the lower left one.
        let start = i * side + n;
        polynomial.row(point, &mut matrix[start..start + r]);
        for j in 0..r {
            matrix[(n + j) * side + i] = matrix[start + j];
        }
    }
    let mut rhs = vec![V::zero(); side];
    rhs[..n].copy_from_slice(values);
    Ok((matrix, rhs))
}

/// Writes the row of the point `at` to `out`: the kernel between it and every
/// point, then every monomial. Its inner product with the solution of the
/// [`system`] is the value of the interpolant there.
///
/// `out` has one entry a point and one a monomial.
pub(crate) fn evaluation_row<const D: usize>(
    at: &[f64; D],
    points: &[[f64; D]],
    kernel: Kernel,
    polynomial: &Polynomial<D>,
    out: &mut [f64],
) {
    let (kernels, monomials) = out.split_at_mut(points.len());
    for (entry, point) in kernels.iter_mut().zip(points) {
        *entry = kernel.between(at, point);
    }
    polynomial.row(at, monomials);
}

#[cfg(test)]
mod tests {
    use super::{evaluation_row, monomial_count, system, Kernel, Polynomial};
    use crate::common::{assert_close, load, reals, Real, Tol, DIRECT_F64};
    use crate::{InterpolateError as E, ParamReason};
    use serde_json::Value as Json;

    fn nums(v: &Json) -> Vec<f64> {
        reals(&serde_json::from_value::<Vec<Real>>(v.clone()).unwrap_or_default())
    }

    fn rows(v: &Json) -> Vec<Vec<f64>> {
        serde_json::from_value::<Vec<Vec<Real>>>(v.clone())
            .unwrap_or_default()
            .iter()
            .map(|r| reals(r))
            .collect()
    }

    fn points_of<const D: usize>(v: &Json) -> Vec<[f64; D]> {
        rows(v)
            .iter()
            .map(|r| r.as_slice().try_into().expect("D coordinates"))
            .collect()
    }

    fn bits(v: &[f64]) -> Vec<u64> {
        v.iter().map(|x| x.to_bits()).collect()
    }

    /// Every kernel, with `epsilon` where it takes one.
    fn kernels(epsilon: f64) -> [(&'static str, Kernel); 8] {
        [
            ("linear", Kernel::Linear),
            ("thin_plate_spline", Kernel::ThinPlateSpline),
            ("cubic", Kernel::Cubic),
            ("quintic", Kernel::Quintic),
            ("multiquadric", Kernel::Multiquadric { epsilon }),
            (
                "inverse_multiquadric",
                Kernel::InverseMultiquadric { epsilon },
            ),
            ("inverse_quadratic", Kernel::InverseQuadratic { epsilon }),
            ("gaussian", Kernel::Gaussian { epsilon }),
        ]
    }

    /// One golden case in `D` dimensions: every part against SciPy's.
    fn check<const D: usize>(case: &Json) {
        let name = case["name"].as_str().expect("name");
        let options = &case["options"];
        let epsilon = options["epsilon"].as_f64().expect("epsilon");
        let (_, kernel) = kernels(epsilon)
            .into_iter()
            .find(|(named, _)| options["kernel"] == *named)
            .expect("a kernel");
        // -1 is SciPy's way to write no polynomial terms.
        let degree = usize::try_from(options["degree"].as_i64().expect("degree")).ok();
        let smoothing = options["smoothing"].as_f64().expect("smoothing");
        let points: Vec<[f64; D]> = points_of(&case["points"]);
        let values = nums(&case["values"]);

        let polynomial = Polynomial::new(&points, kernel, degree).expect("valid terms");
        let powers: Vec<Vec<i32>> = serde_json::from_value(case["powers"].clone()).expect("powers");
        assert_eq!(polynomial.len(), powers.len(), "{name}");
        for (got, want) in polynomial.powers.iter().zip(&powers) {
            assert_eq!(got.as_slice(), want.as_slice(), "{name}: the monomials");
        }
        // The middle and half the range: the same two operations, the same bits.
        assert_eq!(
            bits(&polynomial.shift),
            bits(&nums(&case["shift"])),
            "{name}"
        );
        assert_eq!(
            bits(&polynomial.scale),
            bits(&nums(&case["scale"])),
            "{name}"
        );

        let (matrix, rhs) =
            system(&points, &values, smoothing, kernel, &polynomial).expect("a system");
        let side = points.len() + powers.len();
        let want: Vec<f64> = rows(&case["matrix"]).concat();
        assert_eq!(
            (matrix.len(), want.len()),
            (side * side, side * side),
            "{name}"
        );
        for (k, (got, want)) in matrix.iter().zip(&want).enumerate() {
            let what = format!("{name}: matrix[{}][{}]", k / side, k % side);
            assert_close(*got, *want, DIRECT_F64, &what);
        }
        assert_eq!(
            bits(&rhs),
            bits(&nums(&case["rhs"])),
            "{name}: right-hand side"
        );

        let queries: Vec<[f64; D]> = points_of(&case["queries"]);
        let mut row = vec![0.0; side];
        for (q, (query, want)) in queries.iter().zip(rows(&case["expected"])).enumerate() {
            evaluation_row(query, &points, kernel, &polynomial, &mut row);
            for (k, (got, want)) in row.iter().zip(&want).enumerate() {
                let what = format!("{name}: row {q}, entry {k}");
                assert_close(*got, *want, DIRECT_F64, &what);
            }
        }
    }

    #[test]
    fn the_parts_are_those_of_scipys_internal_functions() {
        let file = load::<Json>("rbf_system").expect("golden file rbf_system");
        let mut by_dimension = [0; 4];
        let mut seen = std::collections::BTreeSet::new();
        for case in &file.cases {
            let dimension = rows(&case["points"])[0].len();
            match dimension {
                1 => check::<1>(case),
                2 => check::<2>(case),
                3 => check::<3>(case),
                other => panic!("points in {other} dimensions"),
            }
            by_dimension[dimension] += 1;
            seen.insert(
                case["options"]["kernel"]
                    .as_str()
                    .expect("kernel")
                    .to_owned(),
            );
        }
        assert_eq!(by_dimension, [0, 1, 15, 1]);
        assert_eq!(seen.len(), 8, "every kernel");
    }

    #[test]
    fn each_kernel_has_scipys_formula_and_sign() {
        // SciPy's `_kernel_matrix` at the distance 0.7.
        let want = [
            -0.7,
            -0.174_770_722_529_978_88,
            0.343,
            -0.168_07,
            -1.220_655_561_573_370_3,
            0.819_231_920_519_040_5,
            0.671_140_939_597_315_5,
            0.612_626_394_184_416_1,
        ];
        let tight = Tol {
            rtol: 1e-15,
            atol: 0.0,
        };
        for ((name, kernel), want) in kernels(1.0).into_iter().zip(want) {
            assert_close(kernel.at(0.7), want, tight, name);
            assert_close(kernel.between(&[0.0, 0.3], &[0.7, 0.3]), want, tight, name);
        }
        // The distance is multiplied by `epsilon` where the kernel has one:
        // half the distance and twice the `epsilon` is the same value.
        for ((name, kernel), want) in kernels(2.0).into_iter().zip(want).skip(4) {
            assert_close(kernel.between(&[0.35], &[0.0]), want, tight, name);
        }
        for (name, kernel) in kernels(2.0).into_iter().take(4) {
            let (got, want) = (kernel.between(&[0.35], &[0.0]), kernel.at(0.35));
            assert_eq!(got.to_bits(), want.to_bits(), "{name} takes no epsilon");
        }
        // At the distance 0: the thin plate spline is 0, not the NaN of
        // `0 * log 0`.
        let at_zero = [0.0, 0.0, 0.0, 0.0, -1.0, 1.0, 1.0, 1.0];
        for ((name, kernel), want) in kernels(3.0).into_iter().zip(at_zero) {
            let got = kernel.between(&[1.5, -2.0], &[1.5, -2.0]);
            assert_eq!((got + 0.0).to_bits(), f64::to_bits(want + 0.0), "{name}");
        }
    }

    #[test]
    fn a_degree_is_checked_against_the_kernel_and_the_points() {
        let points = [
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
            [0.5, 0.2],
            [0.1, 0.9],
        ];
        let smallest = [
            Some(0),
            Some(1),
            Some(1),
            Some(2),
            Some(0),
            None,
            None,
            None,
        ];
        for ((name, kernel), min) in kernels(1.0).into_iter().zip(smallest) {
            assert_eq!(kernel.min_degree(), min, "{name}");
            // The default: the smallest, and a constant where none is needed.
            assert_eq!(kernel.default_degree(), min.unwrap_or(0), "{name}");
            let terms = |degree| Polynomial::new(&points, kernel, degree).map(|p| p.len());
            // No terms at all is taken by every kernel, as in SciPy.
            assert_eq!(terms(None), Ok(0), "{name}");
            for degree in 0..=2 {
                let want = match min {
                    Some(min) if degree < min => Err(E::InvalidParameter {
                        name: "degree",
                        reason: ParamReason::TooSmall { min },
                    }),
                    // 1, 3, and 6 monomials in two variables.
                    _ => Ok((degree + 1) * (degree + 2) / 2),
                };
                assert_eq!(terms(Some(degree)), want, "{name} degree {degree}");
            }
            // Degree 3 has 10 monomials, and there are 6 points.
            assert_eq!(
                terms(Some(3)),
                Err(E::InsufficientPoints {
                    required: 10,
                    actual: 6
                }),
                "{name}"
            );
        }
        // One point fewer than the monomials is too few, as many is enough.
        let gaussian = Kernel::Gaussian { epsilon: 1.0 };
        assert_eq!(
            Polynomial::new(&points[..5], gaussian, Some(2)).err(),
            Some(E::InsufficientPoints {
                required: 6,
                actual: 5
            })
        );
        assert!(Polynomial::new(&points[..3], gaussian, Some(1)).is_ok());
        // Without points and without terms there is nothing to shift by: the
        // coordinates are left as they are.
        let none = Polynomial::new(&points[..0], gaussian, None).expect("no terms");
        assert_eq!((none.len(), bits(&none.shift)), (0, bits(&[0.0, 0.0])));
        assert_eq!(bits(&none.scale), bits(&[1.0, 1.0]));
        // A degree under the smallest is reported before too few points.
        assert_eq!(
            Polynomial::new(&points[..0], Kernel::Quintic, Some(1)).err(),
            Some(E::InvalidParameter {
                name: "degree",
                reason: ParamReason::TooSmall { min: 2 },
            })
        );
    }

    /// `C(n, k)` by Pascal's triangle, or `None` when it does not fit.
    fn binomial(n: usize, k: usize) -> Option<u128> {
        let mut row = vec![1_u128];
        for _ in 0..n {
            let mut next = vec![1_u128];
            for pair in row.windows(2) {
                next.push(pair[0].checked_add(pair[1])?);
            }
            next.push(1);
            row = next;
        }
        row.get(k).copied()
    }

    fn monomials_in<const D: usize>() {
        // Enough points for every count below.
        let points = vec![[0.5; D]; 300];
        let gaussian = Kernel::Gaussian { epsilon: 1.0 };
        for degree in 0..=5_usize {
            let terms = Polynomial::new(&points, gaussian, Some(degree)).expect("valid terms");
            let count = binomial(D + degree, D).expect("small");
            assert_eq!(terms.len() as u128, count, "{D} variables, degree {degree}");
            assert_eq!(monomial_count(D, degree).map(|c| c as u128), Some(count));
            let totals: Vec<i32> = terms.powers.iter().map(|p| p.iter().sum()).collect();
            assert!(
                totals.windows(2).all(|pair| pair[0] <= pair[1]),
                "by degree"
            );
            assert!(totals.iter().all(|total| *total <= degree as i32));
            let distinct: std::collections::BTreeSet<_> = terms.powers.iter().collect();
            assert_eq!(distinct.len(), terms.len(), "each once");
            assert!(terms.powers.iter().flatten().all(|p| *p >= 0));
        }
    }

    #[test]
    fn the_monomials_are_all_of_the_degree_or_less_once_each_by_degree() {
        monomials_in::<0>();
        monomials_in::<1>();
        monomials_in::<2>();
        monomials_in::<3>();
        monomials_in::<4>();
        // Two variables, degree 2, in SciPy's order: 1, x, y, x^2, x y, y^2.
        let points = [[0.0, 0.0]; 6];
        let terms = Polynomial::new(&points, Kernel::Linear, Some(2)).expect("valid terms");
        assert_eq!(
            terms.powers,
            [[0, 0], [1, 0], [0, 1], [2, 0], [1, 1], [0, 2]]
        );
    }

    #[test]
    fn a_count_of_monomials_that_does_not_fit_names_the_largest_degree_that_does() {
        let too_large = |max| {
            Some(E::InvalidParameter {
                name: "degree",
                reason: ParamReason::TooLarge { max },
            })
        };
        let gaussian = Kernel::Gaussian { epsilon: 1.0 };
        // In one variable the count always fits; the power must fit `i32`.
        let line: [[f64; 1]; 0] = [];
        let largest_power = i32::MAX as usize;
        for degree in [largest_power + 1, largest_power + 2, usize::MAX] {
            assert_eq!(
                Polynomial::new(&line, gaussian, Some(degree)).err(),
                too_large(largest_power)
            );
        }
        // Under that, it is the points that are too few.
        assert_eq!(
            Polynomial::new(&line, gaussian, Some(largest_power)).err(),
            Some(E::InsufficientPoints {
                required: largest_power + 1,
                actual: 0
            })
        );
        // In 40 variables the count is what does not fit: the largest degree
        // whose count does, by Pascal's triangle.
        let fits = |degree| binomial(40 + degree, 40).is_some_and(|c| c <= usize::MAX as u128);
        let largest = (0..200)
            .take_while(|degree| fits(*degree))
            .last()
            .expect("degree 0 fits");
        assert!(largest > 5 && !fits(largest + 1));
        let space: [[f64; 40]; 0] = [];
        for degree in [largest + 1, largest + 2, 1000, usize::MAX] {
            assert_eq!(
                Polynomial::new(&space, gaussian, Some(degree)).err(),
                too_large(largest),
                "degree {degree}"
            );
        }
        assert!(matches!(
            Polynomial::new(&space, gaussian, Some(largest)).err(),
            Some(E::InsufficientPoints { actual: 0, .. })
        ));
        assert_eq!(monomial_count(40, largest + 1), None);
    }

    #[test]
    fn the_system_is_symmetric_and_a_row_at_a_data_point_is_its_row_of_the_matrix() {
        let points = [
            [0.0, 0.5, 1.0],
            [2.0, -1.0, 0.0],
            [1.0, 1.0, 1.0],
            [-1.0, 0.25, 3.0],
            [0.5, 0.5, -2.0],
        ];
        // Two values at every point.
        let values = [
            [1.0, -1.0],
            [2.0, -2.0],
            [3.0, -3.0],
            [4.0, -4.0],
            [5.0, -5.0],
        ];
        let smoothing = 0.375;
        let mut checked = 0;
        for (name, kernel) in kernels(0.6) {
            let degree = kernel.default_degree();
            let Ok(polynomial) = Polynomial::new(&points, kernel, Some(degree)) else {
                // The quintic needs degree 2: 10 monomials for 5 points.
                assert_eq!(kernel, Kernel::Quintic);
                continue;
            };
            let (n, side) = (points.len(), points.len() + polynomial.len());
            let (matrix, rhs) =
                system(&points, &values, smoothing, kernel, &polynomial).expect("a system");
            for i in 0..side {
                for j in 0..side {
                    let (here, there) = (matrix[i * side + j], matrix[j * side + i]);
                    assert_eq!(here.to_bits(), there.to_bits(), "{name}: [{i}][{j}]");
                    if i >= n && j >= n {
                        assert_eq!(here.to_bits(), 0.0_f64.to_bits(), "{name}");
                    }
                }
            }
            // The values, then a zero for every monomial.
            assert_eq!(rhs.len(), side, "{name}");
            for (got, want) in rhs.iter().zip(values.iter().chain(&[[0.0, 0.0]; 4])) {
                assert_eq!(bits(got), bits(want), "{name}");
            }
            // The row of an evaluation at a data point: its row of the
            // matrix, without the smoothing on the diagonal.
            let mut row = vec![0.0; side];
            for i in 0..n {
                evaluation_row(&points[i], &points, kernel, &polynomial, &mut row);
                for j in 0..side {
                    let on_diagonal = if i == j { smoothing } else { 0.0 };
                    let want = matrix[i * side + j] - on_diagonal;
                    let what = format!("{name}: point {i}, entry {j}");
                    assert_close(row[j], want, DIRECT_F64, &what);
                }
            }
            checked += 1;
        }
        assert_eq!(checked, 7);
    }

    #[test]
    fn a_matrix_of_several_tiles_has_every_kernel_at_both_of_its_places() {
        // More points than two tiles are wide, at spacings that differ.
        let points: Vec<[f64; 2]> = (0..150_u32)
            .map(|k| [f64::from(k * 37 % 150) * 0.11, f64::from(k * k % 97) * 0.07])
            .collect();
        let values: Vec<f64> = (0..150).map(f64::from).collect();
        let smoothing = 0.5;
        for (name, kernel) in kernels(0.3) {
            let terms = Polynomial::new(&points, kernel, Some(kernel.default_degree()))
                .expect("valid terms");
            let side = points.len() + terms.len();
            let (matrix, _) =
                system(&points, &values, smoothing, kernel, &terms).expect("a system");
            let mut monomials = vec![0.0; terms.len()];
            for (i, a) in points.iter().enumerate() {
                for (j, b) in points.iter().enumerate() {
                    let on_diagonal = if i == j { smoothing } else { 0.0 };
                    let want = kernel.between(a, b) + on_diagonal;
                    let got = matrix[i * side + j];
                    assert_eq!(got.to_bits(), want.to_bits(), "{name}: [{i}][{j}]");
                }
                terms.row(a, &mut monomials);
                for (j, want) in monomials.iter().enumerate() {
                    let (right, below) = (matrix[i * side + 150 + j], matrix[(150 + j) * side + i]);
                    assert_eq!(right.to_bits(), want.to_bits(), "{name}: P[{i}][{j}]");
                    assert_eq!(below.to_bits(), want.to_bits(), "{name}: P^T[{j}][{i}]");
                }
            }
        }
    }

    #[test]
    fn a_matrix_larger_than_a_vec_holds_is_refused_before_it_is_asked_for() {
        // Points without coordinates and values without columns take no
        // memory, so their number can be anything.
        #[cfg(target_pointer_width = "64")]
        {
            let count = 1_usize << 31;
            let (points, values) = (vec![[0.0_f64; 0]; count], vec![[0.0_f64; 0]; count]);
            // The largest matrix of `f64` a `Vec` holds has under 2^60
            // entries: a side of 2^30 is one too many. One row and one column
            // are for the constant term.
            let terms = Polynomial::new(&points, Kernel::Linear, Some(0)).expect("valid terms");
            assert_eq!(terms.len(), 1);
            let most = (isize::MAX.unsigned_abs() / 8).isqrt();
            assert_eq!(most, (1 << 30) - 1);
            assert_eq!(
                system(&points, &values, 0.0, Kernel::Linear, &terms).err(),
                Some(E::InvalidParameter {
                    name: "points",
                    reason: ParamReason::TooLarge { max: most - 1 },
                })
            );
        }
    }
}
