//! Tridiagonal linear systems: what the construction of a cubic spline solves.
//!
//! A matrix of `n` rows is three slices of length `n`: `lower[i]` is the entry
//! left of the diagonal in row `i`, `diag[i]` the one on it, and `upper[i]` the
//! one right of it. `lower[0]` and `upper[n - 1]` are outside the matrix:
//! [`solve`] does not read them, and [`solve_cyclic`] takes them for the two
//! corners of a cyclic matrix.
//!
//! The solvers work in place and do not allocate. The elimination leaves the
//! reciprocals of the pivots in `diag`, and the substitutions multiply by them.
//!
//! A solve is a recurrence: every row waits for the row before it. What limits
//! its speed is the chain of operations between one row's result and the next
//! one's, more than how many operations there are. So each loop here is written
//! with a short chain: a pivot waits for the pivot above through one division
//! and one subtraction, and a substitution through one multiplication and one
//! addition. Everything else in a row is computed beside the chain. (Counting
//! divisions instead, one per row with a longer chain, was no faster; measured
//! in `benches/cubic_spline.rs`.)
//!
//! The elimination is without pivoting (the Thomas algorithm), which is stable
//! for the matrices of the spline boundary conditions (docs/design.ja.md
//! section 6.1). A pivot that is zero or not finite is reported as a
//! [`Breakdown`]; nothing panics. Whether the solution is finite is for the
//! caller to check.

use crate::{Scalar, Value};

/// The elimination met a pivot that is zero, NaN, or infinite: the matrix is
/// singular, or too close to it for an elimination without pivoting.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Breakdown {
    /// The row of the pivot.
    pub(crate) row: usize,
}

/// Solves `A x = rhs` and leaves `x` in `rhs`. `diag` is left holding the
/// reciprocals of the pivots.
///
/// # Errors
///
/// Returns [`Breakdown`] when a pivot is zero or not finite. `diag` and `rhs`
/// are then partly overwritten.
pub(crate) fn solve<T: Scalar, V: Value<T>>(
    lower: &[T],
    diag: &mut [T],
    upper: &[T],
    rhs: &mut [V],
) -> Result<(), Breakdown> {
    // No second right-hand side.
    let mut none: [T; 0] = [];
    sweep(lower, diag, upper, rhs, &mut none)
}

/// Solves `A x = rhs` for a cyclic matrix and leaves `x` in `rhs`.
///
/// The matrix has, besides its three diagonals, `lower[0]` in the top right
/// corner (row `0`, column `n - 1`) and `upper[n - 1]` in the bottom left one.
/// `work` is scratch space of length `n`; `diag` and `work` are overwritten.
///
/// `n` is 3 or more: a smaller cyclic matrix has no one shape. The caller
/// guarantees this; it is not checked in a release build, where 1 or 2 rows give
/// a wrong solution without an error and no rows panic. A periodic spline has
/// one row fewer than it has points, so one through three points has two rows:
/// its caller handles that case apart (docs/design.ja.md section 6.1).
///
/// This is the Sherman-Morrison formula. With `gamma = -diag[0]`, the matrix is
/// `A' + u v^T` for `u = (gamma, 0, ..., 0, bottom_left)`,
/// `v = (1, 0, ..., 0, top_right / gamma)`, and a tridiagonal `A'` that differs
/// from `A` in the two ends of its diagonal. `A' y = rhs` and `A' z = u` are
/// solved with one elimination, and `x = y - z (v . y) / (1 + v . z)`.
///
/// # Errors
///
/// Returns [`Breakdown`] when a pivot of `A'` is zero or not finite (row 0 when
/// `diag[0]` is zero), and with row `n` when `1 + v . z` is: `A` is singular
/// although `A'` is not.
pub(crate) fn solve_cyclic<T: Scalar, V: Value<T>>(
    lower: &[T],
    diag: &mut [T],
    upper: &[T],
    rhs: &mut [V],
    work: &mut [T],
) -> Result<(), Breakdown> {
    let n = diag.len();
    debug_assert!(n >= 3, "a cyclic system has three rows at least");
    debug_assert_eq!(work.len(), n);
    let last = n - 1;
    let (top_right, bottom_left) = (lower[0], upper[last]);
    let gamma = -diag[0];
    // A `gamma` that is zero or not finite leaves a first pivot that
    // `eliminate` turns away, so it needs no test of its own.
    let ratio = top_right / gamma;
    diag[0] = diag[0] - gamma;
    diag[last] = diag[last] - bottom_left * ratio;
    work.fill(T::zero());
    work[0] = gamma;
    work[last] = bottom_left;
    sweep(lower, diag, upper, rhs, work)?;

    let denominator = T::one() + work[0] + ratio * work[last];
    if unusable(denominator) {
        return Err(Breakdown { row: n });
    }
    // `v . y`, over the denominator.
    let factor = V::zero().axpy(denominator.recip(), rhs[0].axpy(ratio, rhs[last]));
    for (x, z) in rhs.iter_mut().zip(work.iter()) {
        *x = x.axpy(-*z, factor);
    }
    Ok(())
}

/// Whether a pivot cannot be divided by.
fn unusable<T: Scalar>(pivot: T) -> bool {
    pivot == T::zero() || !pivot.is_finite()
}

/// The pivot of `row`, given the pivot `above` it. `diag` still holds the entry
/// of the matrix at `row`.
///
/// The product of the two entries beside the diagonal does not wait for
/// `above`, so the pivot follows it by one division and one subtraction.
fn next_pivot<T: Scalar>(lower: &[T], diag: &[T], upper: &[T], row: usize, above: T) -> T {
    if row == 0 {
        diag[0]
    } else {
        diag[row] - lower[row] * upper[row - 1] / above
    }
}

/// Eliminates, and solves for `rhs` and for `second` with the same matrix.
/// `second` is as long as `rhs`, or empty when there is one right-hand side.
/// `diag` is left holding the reciprocals of the pivots.
///
/// One pass down does the elimination and the forward substitutions, one pass up
/// the back substitutions. In each pass the chains of the pivots and of the
/// right-hand sides run side by side: done in passes of their own, they would
/// add up.
fn sweep<T: Scalar, V: Value<T>>(
    lower: &[T],
    diag: &mut [T],
    upper: &[T],
    rhs: &mut [V],
    second: &mut [T],
) -> Result<(), Breakdown> {
    let n = diag.len();
    debug_assert_eq!(lower.len(), n);
    debug_assert_eq!(upper.len(), n);
    debug_assert_eq!(rhs.len(), n);
    debug_assert!(second.is_empty() || second.len() == n);
    let mut above = T::one();
    for row in 0..n {
        let pivot = next_pivot(lower, diag, upper, row, above);
        let inverse = pivot.recip();
        // A pivot too small to invert is as good as zero.
        if unusable(pivot) || !inverse.is_finite() {
            return Err(Breakdown { row });
        }
        if row > 0 {
            let multiplier = lower[row] * diag[row - 1];
            rhs[row] = rhs[row].axpy(-multiplier, rhs[row - 1]);
            if let Some(pair) = second.get_mut(row - 1..=row) {
                pair[1] = pair[1] - multiplier * pair[0];
            }
        }
        diag[row] = inverse;
        above = pivot;
    }
    // A row is scaled by its reciprocal pivot before the row below is taken off
    // it, and not after: the scaling then does not wait for the row below.
    for row in (0..n).rev() {
        let inverse = diag[row];
        let factor = upper[row] * inverse;
        let scaled = V::zero().axpy(inverse, rhs[row]);
        rhs[row] = match rhs.get(row + 1) {
            Some(below) => scaled.axpy(-factor, *below),
            None => scaled,
        };
        if let Some(this) = second.get(row).copied() {
            let below = second.get(row + 1).copied().unwrap_or_else(T::zero);
            second[row] = this * inverse - factor * below;
        }
    }
    Ok(())
}

// The golden loader of the integration tests. An integration test cannot call
// this crate-private module, so its tests are below and read the file the same
// way (docs/conventions.md, "Golden data"). The module is declared here and not
// inside `mod tests`: there the path would start from `src/tridiagonal/tests/`,
// a directory that does not exist, and only Windows resolves `..` through one.
#[cfg(test)]
#[path = "../tests/common/mod.rs"]
mod common;

#[cfg(test)]
mod tests {
    use super::common::{assert_close, load, reals, Real, Tol, DIRECT_F64};
    use super::{solve, solve_cyclic, Breakdown};
    use proptest::prelude::*;
    use serde_json::Value as Json;

    fn nums(v: &Json) -> Vec<f64> {
        reals(&serde_json::from_value::<Vec<Real>>(v.clone()).unwrap_or_default())
    }

    /// The right-hand side or the solution of a case: a list, or a list of rows
    /// of `M` columns.
    fn columns<const M: usize>(v: &Json) -> Option<Vec<[f64; M]>> {
        let rows: Vec<Vec<Real>> = match serde_json::from_value::<Vec<Vec<Real>>>(v.clone()) {
            Ok(rows) => rows,
            // A plain list is one column.
            Err(_) => nums(v).into_iter().map(|x| vec![Real(x)]).collect(),
        };
        rows.iter().map(|r| reals(r).try_into().ok()).collect()
    }

    fn is_cyclic(case: &Json) -> bool {
        case["options"]["cyclic"] == Json::Bool(true)
    }

    /// Solves one golden case with `M` columns, or returns `None` when the case
    /// has another number of columns.
    fn run<const M: usize>(case: &Json) -> Option<Result<Vec<[f64; M]>, Breakdown>> {
        let mut rhs = columns::<M>(&case["rhs"])?;
        let (lower, upper) = (nums(&case["lower"]), nums(&case["upper"]));
        let mut diag = nums(&case["diag"]);
        let done = if is_cyclic(case) {
            let mut work = vec![0.0; diag.len()];
            solve_cyclic(&lower, &mut diag, &upper, &mut rhs, &mut work)
        } else {
            solve(&lower, &mut diag, &upper, &mut rhs)
        };
        Some(done.map(|()| rhs))
    }

    /// Checks one golden case with `M` columns; `false` when it has another
    /// number of columns.
    fn check<const M: usize>(case: &Json, tol: Tol) -> bool {
        let Some(got) = run::<M>(case) else {
            return false;
        };
        let name = case["name"].as_str().expect("name");
        let got = got.expect("SciPy solved this system");
        let want = columns::<M>(&case["expected"]).expect("as many columns as the right-hand side");
        assert_eq!(got.len(), want.len(), "{name}");
        for (row, (g, w)) in got.iter().zip(&want).enumerate() {
            for column in 0..M {
                let what = format!("{name} row {row} column {column}");
                assert_close(g[column], w[column], tol, &what);
            }
        }
        true
    }

    #[test]
    fn solutions_match_scipy() {
        let file = load::<Json>("tridiagonal").expect("golden file tridiagonal");
        let (mut plain, mut cyclic, mut wide, mut uneven) = (0, 0, 0, 0);
        for case in &file.cases {
            let name = case["name"].as_str().expect("name");
            // The same tolerance for every case, the ill-conditioned ones too:
            // spacings of 1 and 1e-6 side by side, which SciPy solves with
            // pivoting and this solver without.
            let tol = DIRECT_F64;
            uneven += usize::from(name.contains("uneven"));
            if check::<1>(case, tol) {
                // One column also goes through the scalar value type, and gives
                // the same bits.
                let mut rhs = nums(&case["rhs"]);
                let (lower, upper) = (nums(&case["lower"]), nums(&case["upper"]));
                let mut diag = nums(&case["diag"]);
                if is_cyclic(case) {
                    let mut work = vec![0.0; diag.len()];
                    solve_cyclic(&lower, &mut diag, &upper, &mut rhs, &mut work).expect("solved");
                    cyclic += 1;
                } else {
                    solve(&lower, &mut diag, &upper, &mut rhs).expect("solved");
                    plain += 1;
                }
                let as_arrays = run::<1>(case).expect("one column").expect("solved");
                for (scalar, array) in rhs.iter().zip(&as_arrays) {
                    assert_eq!(scalar.to_bits(), array[0].to_bits(), "{name}");
                }
            } else {
                assert!(check::<2>(case, tol) || check::<3>(case, tol), "{name}");
                wide += 1;
            }
        }
        assert_eq!((plain, cyclic, wide, uneven), (33, 10, 2, 5));
    }

    /// Gaussian elimination with partial pivoting on a dense matrix: the
    /// reference the properties compare with, written apart from the solver.
    fn dense_solve(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Vec<f64> {
        let n = b.len();
        for k in 0..n {
            let pivot = (k..n)
                .max_by(|i, j| a[*i][k].abs().total_cmp(&a[*j][k].abs()))
                .expect("a row is left");
            a.swap(k, pivot);
            b.swap(k, pivot);
            for i in k + 1..n {
                let m = a[i][k] / a[k][k];
                let pivot_row = a[k].clone();
                for (entry, above) in a[i].iter_mut().zip(&pivot_row).skip(k) {
                    *entry -= m * above;
                }
                b[i] -= m * b[k];
            }
        }
        for k in (0..n).rev() {
            let tail: f64 = (k + 1..n).map(|j| a[k][j] * b[j]).sum();
            b[k] = (b[k] - tail) / a[k][k];
        }
        b
    }

    /// The dense form of three diagonals; with `cyclic`, the two corners as well.
    fn dense(lower: &[f64], diag: &[f64], upper: &[f64], cyclic: bool) -> Vec<Vec<f64>> {
        let n = diag.len();
        let mut a = vec![vec![0.0; n]; n];
        for i in 0..n {
            a[i][i] = diag[i];
            if i > 0 {
                a[i][i - 1] = lower[i];
            }
            if i + 1 < n {
                a[i][i + 1] = upper[i];
            }
        }
        if cyclic {
            a[0][n - 1] = lower[0];
            a[n - 1][0] = upper[n - 1];
        }
        a
    }

    type System = (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>);

    /// Three diagonals of a matrix whose diagonal exceeds the sum of the entries
    /// beside it by 0.1 at least, and a right-hand side. `least` rows or more.
    fn systems(least: usize) -> impl Strategy<Value = System> {
        (least..40_usize)
            .prop_flat_map(|n| {
                (
                    prop::collection::vec(-1.0_f64..1.0, n),
                    prop::collection::vec(-1.0_f64..1.0, n),
                    prop::collection::vec((0.1_f64..2.0, any::<bool>()), n),
                    prop::collection::vec(-10.0_f64..10.0, n),
                )
            })
            .prop_map(|(lower, upper, margins, rhs)| {
                let diag = lower
                    .iter()
                    .zip(&upper)
                    .zip(margins)
                    .map(|((l, u), (margin, negative))| {
                        let size = l.abs() + u.abs() + margin;
                        if negative {
                            -size
                        } else {
                            size
                        }
                    })
                    .collect();
                (lower, diag, upper, rhs)
            })
    }

    proptest! {
        #[test]
        fn a_dominant_system_is_solved_as_the_dense_elimination_solves_it(
            (lower, diag, upper, rhs) in systems(1),
        ) {
            let want = dense_solve(dense(&lower, &diag, &upper, false), rhs.clone());
            let (mut pivots, mut got) = (diag.clone(), rhs.clone());
            prop_assert_eq!(solve(&lower, &mut pivots, &upper, &mut got), Ok(()));
            for (g, w) in got.iter().zip(&want) {
                prop_assert!((g - w).abs() <= 1e-11 * w.abs().max(1.0), "{} for {}", g, w);
            }
            // The residual, row by row.
            for i in 0..rhs.len() {
                let mut sum = diag[i] * got[i];
                if i > 0 {
                    sum += lower[i] * got[i - 1];
                }
                if i + 1 < rhs.len() {
                    sum += upper[i] * got[i + 1];
                }
                prop_assert!((sum - rhs[i]).abs() <= 1e-10, "row {}: {} for {}", i, sum, rhs[i]);
            }
        }

        #[test]
        fn a_dominant_cyclic_system_is_solved_as_the_dense_elimination_solves_it(
            (lower, diag, upper, rhs) in systems(3),
        ) {
            let want = dense_solve(dense(&lower, &diag, &upper, true), rhs.clone());
            // Whatever the scratch space held before must not matter.
            let (mut pivots, mut got, mut work) = (diag.clone(), rhs.clone(), vec![7.0; rhs.len()]);
            prop_assert_eq!(solve_cyclic(&lower, &mut pivots, &upper, &mut got, &mut work), Ok(()));
            for (g, w) in got.iter().zip(&want) {
                prop_assert!((g - w).abs() <= 1e-10 * w.abs().max(1.0), "{} for {}", g, w);
            }
        }

        #[test]
        fn the_entries_outside_the_matrix_are_not_read(
            (mut lower, diag, mut upper, rhs) in systems(1),
        ) {
            let (mut pivots, mut want) = (diag.clone(), rhs.clone());
            prop_assert_eq!(solve(&lower, &mut pivots, &upper, &mut want), Ok(()));
            lower[0] = f64::NAN;
            let last = upper.len() - 1;
            upper[last] = f64::NAN;
            let (mut pivots, mut got) = (diag.clone(), rhs.clone());
            prop_assert_eq!(solve(&lower, &mut pivots, &upper, &mut got), Ok(()));
            for (g, w) in got.iter().zip(&want) {
                prop_assert_eq!(g.to_bits(), w.to_bits());
            }
        }

        #[test]
        fn columns_are_solved_as_if_each_were_alone(
            (lower, diag, upper, first) in systems(3),
            scale in -3.0_f64..3.0,
            cyclic in any::<bool>(),
        ) {
            let second: Vec<f64> = first.iter().rev().map(|v| v * scale + 1.0).collect();
            let run = |rhs: &mut [f64]| {
                let mut pivots = diag.clone();
                if cyclic {
                    solve_cyclic(&lower, &mut pivots, &upper, rhs, &mut vec![0.0; diag.len()])
                } else {
                    solve(&lower, &mut pivots, &upper, rhs)
                }
            };
            let (mut a, mut b) = (first.clone(), second.clone());
            prop_assert_eq!(run(&mut a), Ok(()));
            prop_assert_eq!(run(&mut b), Ok(()));

            let mut both: Vec<[f64; 2]> = first.iter().zip(&second).map(|(x, y)| [*x, *y]).collect();
            let mut pivots = diag.clone();
            let done = if cyclic {
                solve_cyclic(&lower, &mut pivots, &upper, &mut both, &mut vec![0.0; diag.len()])
            } else {
                solve(&lower, &mut pivots, &upper, &mut both)
            };
            prop_assert_eq!(done, Ok(()));
            for (i, pair) in both.iter().enumerate() {
                prop_assert_eq!(pair[0].to_bits(), a[i].to_bits());
                prop_assert_eq!(pair[1].to_bits(), b[i].to_bits());
            }
        }

        #[test]
        fn single_precision_agrees_with_double_to_its_own_accuracy(
            (lower, diag, upper, rhs) in systems(3),
            cyclic in any::<bool>(),
        ) {
            // The same numbers rounded to f32, solved in both precisions.
            let narrow = |v: &[f64]| -> Vec<f32> { v.iter().map(|x| *x as f32).collect() };
            let widen = |v: &[f32]| -> Vec<f64> { v.iter().map(|x| f64::from(*x)).collect() };
            let (l32, d32, u32, mut x32) = (narrow(&lower), narrow(&diag), narrow(&upper), narrow(&rhs));
            let (l64, d64, u64, mut x64) = (widen(&l32), widen(&d32), widen(&u32), widen(&x32));
            let (mut p32, mut p64) = (d32.clone(), d64.clone());
            if cyclic {
                prop_assert_eq!(solve_cyclic(&l32, &mut p32, &u32, &mut x32, &mut vec![0.0; d32.len()]), Ok(()));
                prop_assert_eq!(solve_cyclic(&l64, &mut p64, &u64, &mut x64, &mut vec![0.0; d64.len()]), Ok(()));
            } else {
                prop_assert_eq!(solve(&l32, &mut p32, &u32, &mut x32), Ok(()));
                prop_assert_eq!(solve(&l64, &mut p64, &u64, &mut x64), Ok(()));
            }
            // The margin of dominance is 0.1 at least, so the condition number is
            // a few tens at most.
            for (single, double) in x32.iter().zip(&x64) {
                prop_assert!((f64::from(*single) - double).abs() <= 1e-3 * double.abs().max(1.0), "{} for {}", single, double);
            }
        }
    }

    #[test]
    fn the_reciprocal_pivots_are_left_in_the_diagonal() {
        // [[2, 1], [1, 3]]: the pivots are 2 and 3 - 1/2.
        let mut diag = [2.0_f64, 3.0];
        let mut rhs = [3.0_f64, 4.0];
        assert_eq!(solve(&[0.0, 1.0], &mut diag, &[1.0, 0.0], &mut rhs), Ok(()));
        assert_eq!(diag.map(f64::to_bits), [0.5_f64, 0.4].map(f64::to_bits));
        assert_eq!(rhs.map(f64::to_bits), [1.0_f64, 1.0].map(f64::to_bits));
    }

    #[test]
    fn an_empty_system_is_solved() {
        let mut rhs: [f64; 0] = [];
        assert_eq!(solve::<f64, f64>(&[], &mut [], &[], &mut rhs), Ok(()));
    }

    #[test]
    fn a_pivot_that_cannot_be_divided_by_is_reported_with_its_row() {
        let run = |lower: &[f64], diag: &[f64], upper: &[f64]| {
            let mut rhs = vec![1.0; diag.len()];
            solve(lower, &mut diag.to_vec(), upper, &mut rhs)
        };
        // A zero on the diagonal of the first row.
        assert_eq!(
            run(&[0.0, 1.0], &[0.0, 1.0], &[1.0, 0.0]),
            Err(Breakdown { row: 0 })
        );
        // [[1, 1], [1, 1]] is singular: the second pivot is 1 - 1.
        assert_eq!(
            run(&[0.0, 1.0], &[1.0, 1.0], &[1.0, 0.0]),
            Err(Breakdown { row: 1 })
        );
        // The same in the last of three rows: 3 - (2 / 2) * 3.
        assert_eq!(
            run(&[0.0, 0.0, 2.0], &[1.0, 2.0, 3.0], &[0.0, 3.0, 0.0]),
            Err(Breakdown { row: 2 })
        );
        // NaN and infinity, in the diagonal and beside it.
        assert_eq!(
            run(&[0.0, 1.0], &[f64::NAN, 1.0], &[1.0, 0.0]),
            Err(Breakdown { row: 0 })
        );
        assert_eq!(
            run(&[0.0, 1.0], &[1.0, f64::INFINITY], &[1.0, 0.0]),
            Err(Breakdown { row: 1 })
        );
        assert_eq!(
            run(&[0.0, f64::NAN], &[1.0, 1.0], &[1.0, 0.0]),
            Err(Breakdown { row: 1 })
        );
        // A pivot too small to invert: its reciprocal is infinite.
        assert_eq!(
            run(&[0.0, 1.0], &[1e-320, 1.0], &[1.0, 0.0]),
            Err(Breakdown { row: 0 })
        );
        assert_eq!(
            run(&[0.0, 1.0], &[1.0, 1.0 + 1e-320], &[1.0, 0.0]),
            Err(Breakdown { row: 1 })
        );
        // A matrix that needs a row exchange is not for this solver: [[0, 1], [1, 0]].
        assert_eq!(
            run(&[0.0, 1.0], &[0.0, 0.0], &[1.0, 0.0]),
            Err(Breakdown { row: 0 })
        );
    }

    #[test]
    fn a_cyclic_breakdown_is_reported_with_its_row() {
        let run = |lower: &[f64], diag: &[f64], upper: &[f64]| {
            let n = diag.len();
            solve_cyclic(
                lower,
                &mut diag.to_vec(),
                upper,
                &mut vec![1.0; n],
                &mut vec![0.0; n],
            )
        };
        // `gamma` is minus the first diagonal entry.
        assert_eq!(
            run(&[1.0, 1.0, 1.0], &[0.0, 4.0, 4.0], &[1.0, 1.0, 1.0]),
            Err(Breakdown { row: 0 })
        );
        assert_eq!(
            run(&[1.0, 1.0, 1.0], &[f64::NAN, 4.0, 4.0], &[1.0, 1.0, 1.0]),
            Err(Breakdown { row: 0 })
        );
        // The changed matrix breaks down in its last row, whose diagonal entry
        // becomes 1 - 2 * (2 / -4) - 2 = 0.
        assert_eq!(
            run(&[2.0, 0.0, 0.0], &[4.0, 1.0, -1.0], &[0.0, 0.0, 2.0]),
            Err(Breakdown { row: 2 })
        );
        // [[2, 2, 2], [0, 1, 0], [1, 0, 1]] is singular (row 0 is twice the sum
        // of the others), but the changed matrix [[4, 2, 0], [0, 1, 0], [0, 0, 2]]
        // is not: the last division is by 1 - 1/2 - 1/2. Every step is exact.
        assert_eq!(
            run(&[2.0, 0.0, 0.0], &[2.0, 1.0, 1.0], &[2.0, 0.0, 1.0]),
            Err(Breakdown { row: 3 })
        );
    }

    #[test]
    fn a_small_cyclic_system_has_the_solution_worked_by_hand() {
        // [[4, 1, 2], [1, 4, 1], [3, 1, 4]] times (1, 2, 3) is (12, 12, 17).
        let mut rhs = [12.0_f64, 12.0, 17.0];
        let mut work = [9.0; 3];
        let done = solve_cyclic(
            &[2.0, 1.0, 1.0],
            &mut [4.0, 4.0, 4.0],
            &[1.0, 1.0, 3.0],
            &mut rhs,
            &mut work,
        );
        assert_eq!(done, Ok(()));
        for (got, want) in rhs.iter().zip([1.0, 2.0, 3.0]) {
            assert_close(*got, want, DIRECT_F64, "worked by hand");
        }
    }
}
