//! Banded linear systems: what the construction of a B-spline of any degree
//! solves.
//!
//! A matrix of `n` rows with `lower` diagonals below the main one and `upper`
//! above it is stored column by column, as LAPACK's `gbsv` stores it: a column
//! takes `2 * lower + upper + 1` entries, the main diagonal at offset
//! `lower + upper`. [`Band::at`] gives the place of an entry. The first `lower`
//! entries of a column lie above the band: the row exchanges of the
//! elimination fill them, and a caller leaves them alone. A band of
//! `n - 1` diagonals on each side is a dense matrix ([`Band::dense`]).
//!
//! The elimination is Gaussian with partial pivoting: in each column the entry
//! of largest magnitude on or below the diagonal becomes the pivot, by a row
//! exchange. [`factor`] does it once; [`solve`] then takes any number of
//! right-hand sides, of any value type. A pivot that is zero or not finite is
//! reported as a [`Breakdown`]; nothing panics. Whether a solution is finite is
//! for the caller to check.
//!
//! [`Cyclic`] solves a banded matrix whose band wraps around the corners (a
//! periodic spline), by the Woodbury formula.
//!
//! Nothing here allocates: the matrix, the record of the exchanges, and the
//! scratch space of the cyclic solver are the caller's.

use crate::tridiagonal::{unusable, Breakdown};
use crate::{Scalar, Value};

/// The shape of a banded matrix: `n` rows, `lower` diagonals below the main
/// one, `upper` above it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Band {
    /// The number of rows, and of columns.
    pub(crate) n: usize,
    /// The number of diagonals below the main one.
    pub(crate) lower: usize,
    /// The number of diagonals above the main one.
    pub(crate) upper: usize,
}

impl Band {
    /// The shape of a dense matrix of `n` rows: every diagonal is in the band.
    pub(crate) fn dense(n: usize) -> Self {
        let side = n.saturating_sub(1);
        Self {
            n,
            lower: side,
            upper: side,
        }
    }

    /// The number of entries a column takes in the storage.
    fn stride(self) -> usize {
        2 * self.lower + self.upper + 1
    }

    /// The length of the storage of the matrix.
    pub(crate) fn len(self) -> usize {
        self.n * self.stride()
    }

    /// The place of entry (`row`, `column`) in the storage. The entry is in the
    /// band: `column <= row + upper` and `row <= column + lower`.
    pub(crate) fn at(self, row: usize, column: usize) -> usize {
        column * self.stride() + self.lower + self.upper + row - column
    }
}

/// Factors the matrix in `ab` in place, as `P A = L U`, and records the row
/// exchanges in `pivots`.
///
/// `ab` has [`Band::len`] entries and `pivots` has `n`. What `ab` holds above
/// the band on entry does not matter. On return the diagonal holds the
/// reciprocals of the pivots, and [`solve`] multiplies by them.
///
/// # Errors
///
/// Returns [`Breakdown`] with the column of a pivot that is zero or not finite:
/// the matrix is singular, or holds a value that is not finite. `ab` is then
/// partly overwritten.
pub(crate) fn factor<T: Scalar>(
    band: Band,
    ab: &mut [T],
    pivots: &mut [usize],
) -> Result<(), Breakdown> {
    let Band { n, lower, upper } = band;
    let stride = band.stride();
    // The offset of the main diagonal in a column.
    let main = lower + upper;
    debug_assert_eq!(ab.len(), band.len());
    debug_assert_eq!(pivots.len(), n);
    for column in ab.chunks_exact_mut(stride) {
        column[..lower].fill(T::zero());
    }
    // The last column an exchange has carried an entry to: the rows above the
    // pivot row reach no further.
    let mut last = 0;
    for (j, exchanged_with) in pivots.iter_mut().enumerate() {
        let below = lower.min(n - 1 - j);
        let diag = j * stride + main;
        let mut exchange = 0;
        let mut largest = ab[diag].abs();
        for i in 1..=below {
            let size = ab[diag + i].abs();
            if size > largest {
                largest = size;
                exchange = i;
            }
        }
        let pivot = ab[diag + exchange];
        if unusable(pivot) {
            return Err(Breakdown { row: j });
        }
        *exchanged_with = j + exchange;
        last = last.max((j + upper + exchange).min(n - 1));
        if exchange != 0 {
            for column in j..=last {
                let top = column * stride + main + j - column;
                ab.swap(top, top + exchange);
            }
        }
        let inverse = pivot.recip();
        ab[diag] = inverse;
        let (done, rest) = ab.split_at_mut((j + 1) * stride);
        let multipliers = &mut done[diag + 1..=diag + below];
        for multiplier in multipliers.iter_mut() {
            *multiplier = *multiplier * inverse;
        }
        // Row j, scaled by each multiplier, leaves the rows below it.
        for (ahead, column) in rest.chunks_exact_mut(stride).take(last - j).enumerate() {
            let top = main - 1 - ahead;
            let above = column[top];
            for (entry, multiplier) in column[top + 1..=top + below].iter_mut().zip(&*multipliers) {
                *entry = *entry - *multiplier * above;
            }
        }
    }
    Ok(())
}

/// Solves `A x = rhs` with the factors [`factor`] left in `ab` and `pivots`,
/// and leaves `x` in `rhs`.
pub(crate) fn solve<T: Scalar, V: Value<T>>(band: Band, ab: &[T], pivots: &[usize], rhs: &mut [V]) {
    let Band { n, lower, upper } = band;
    let stride = band.stride();
    let main = lower + upper;
    debug_assert_eq!(ab.len(), band.len());
    debug_assert_eq!(rhs.len(), n);
    // L y = P rhs.
    for j in 0..n {
        rhs.swap(j, pivots[j]);
        let below = lower.min(n - 1 - j);
        let diag = j * stride + main;
        let y = rhs[j];
        for (slot, multiplier) in rhs[j + 1..=j + below]
            .iter_mut()
            .zip(&ab[diag + 1..=diag + below])
        {
            *slot = slot.axpy(-*multiplier, y);
        }
    }
    // U x = y, column by column from the last: the exchanges widen the upper
    // band to `lower + upper` diagonals.
    for j in (0..n).rev() {
        let diag = j * stride + main;
        let x = V::zero().axpy(ab[diag], rhs[j]);
        rhs[j] = x;
        let above = main.min(j);
        for (slot, entry) in rhs[j - above..j].iter_mut().zip(&ab[diag - above..diag]) {
            *slot = slot.axpy(-*entry, x);
        }
    }
}

/// [`solve`] for right-hand sides laid out in rows: row `i` is the slab
/// `rhs[i * width..(i + 1) * width]`, and every position in a slab is a system
/// of its own. Each position gets, bit for bit, what [`solve`] gives it alone.
///
/// The lines of a grid along an axis that is not its last are laid out this
/// way: the entries of one line are a whole slab apart. A step of the
/// substitution then works on two slabs, which are contiguous.
pub(crate) fn solve_slabs<T: Scalar, V: Value<T>>(
    band: Band,
    ab: &[T],
    pivots: &[usize],
    rhs: &mut [V],
    width: usize,
) {
    let Band { n, lower, upper } = band;
    let stride = band.stride();
    let main = lower + upper;
    debug_assert_eq!(ab.len(), band.len());
    debug_assert_eq!(rhs.len(), n * width);
    // L y = P rhs.
    for j in 0..n {
        let exchanged = pivots[j];
        if exchanged != j {
            let (head, tail) = rhs.split_at_mut(exchanged * width);
            head[j * width..(j + 1) * width].swap_with_slice(&mut tail[..width]);
        }
        let below = lower.min(n - 1 - j);
        let diag = j * stride + main;
        let (head, tail) = rhs.split_at_mut((j + 1) * width);
        let y = &head[j * width..];
        for (slab, multiplier) in tail
            .chunks_exact_mut(width)
            .zip(&ab[diag + 1..=diag + below])
        {
            for (slot, y) in slab.iter_mut().zip(y) {
                *slot = slot.axpy(-*multiplier, *y);
            }
        }
    }
    // U x = y, row by row from the last.
    for j in (0..n).rev() {
        let diag = j * stride + main;
        let above = main.min(j);
        let (head, tail) = rhs.split_at_mut(j * width);
        let x = &mut tail[..width];
        for slot in x.iter_mut() {
            *slot = V::zero().axpy(ab[diag], *slot);
        }
        for (slab, entry) in head[(j - above) * width..]
            .chunks_exact_mut(width)
            .zip(&ab[diag - above..diag])
        {
            for (slot, x) in slab.iter_mut().zip(x.iter()) {
                *slot = slot.axpy(-*entry, *x);
            }
        }
    }
}

/// The shape of a cyclic banded matrix: a band of `lower` diagonals below the
/// main one and `upper` above it that wraps around, so that row `i` holds
/// entries in the columns `i - lower` to `i + upper`, taken modulo `n`.
///
/// The matrix is given row by row, `lower + upper + 1` entries a row: entry `e`
/// of row `i` is in column `i + e - lower` modulo `n`. When the band is wider
/// than the matrix it wraps onto itself, and entries that share a column add
/// up.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Cyclic {
    /// The number of rows, and of columns.
    pub(crate) n: usize,
    /// The number of diagonals below the main one.
    pub(crate) lower: usize,
    /// The number of diagonals above the main one.
    pub(crate) upper: usize,
}

impl Cyclic {
    /// The number of entries in the two corners together, per row of them: the
    /// rank of what the corners add to a banded matrix.
    fn rank(self) -> usize {
        self.lower + self.upper
    }

    /// Whether the band wraps onto itself, so that the matrix is solved as a
    /// dense one.
    fn is_dense(self) -> bool {
        self.n <= self.rank()
    }

    /// The banded matrix without its corners.
    fn band(self) -> Band {
        Band {
            n: self.n,
            lower: self.lower,
            upper: self.upper,
        }
    }

    /// The length of the scratch space [`solve`](Cyclic::solve) needs.
    pub(crate) fn scratch_len(self) -> usize {
        if self.is_dense() {
            Band::dense(self.n).len()
        } else {
            let rank = self.rank();
            self.band().len() + self.n * rank + Band::dense(rank).len()
        }
    }

    /// The length of the record of exchanges [`solve`](Cyclic::solve) needs.
    pub(crate) fn pivots_len(self) -> usize {
        if self.is_dense() {
            self.n
        } else {
            self.n + self.rank()
        }
    }

    /// The length of the scratch space of values [`solve`](Cyclic::solve)
    /// needs.
    pub(crate) fn carry_len(self) -> usize {
        if self.is_dense() {
            0
        } else {
            self.rank()
        }
    }

    /// Solves `A x = rhs` and leaves `x` in `rhs`.
    ///
    /// `rows` holds the matrix as the type describes it. `scratch`, `pivots`,
    /// and `carry` have the lengths [`scratch_len`](Cyclic::scratch_len),
    /// [`pivots_len`](Cyclic::pivots_len), and [`carry_len`](Cyclic::carry_len);
    /// what they hold on entry does not matter.
    ///
    /// The matrix is `B + U V^T`: `B` is the band without the corners, the
    /// columns of `U` hold the corners, and `V^T` picks the unknowns they
    /// multiply. With `B Z = U` and `B y = rhs`, the solution is
    /// `y - Z (I + V^T Z)^-1 V^T y` (the Woodbury formula).
    ///
    /// # Errors
    ///
    /// Returns [`Breakdown`] when a pivot is zero or not finite, in the band
    /// or in `I + V^T Z`: its `row` is then that of the unknown the pivot
    /// belongs to.
    pub(crate) fn solve<T: Scalar, V: Value<T>>(
        self,
        rows: &[T],
        rhs: &mut [V],
        scratch: &mut [T],
        pivots: &mut [usize],
        carry: &mut [V],
    ) -> Result<(), Breakdown> {
        let Self { n, lower, upper } = self;
        let width = lower + upper + 1;
        debug_assert_eq!(rows.len(), n * width);
        debug_assert_eq!(rhs.len(), n);
        debug_assert_eq!(scratch.len(), self.scratch_len());
        debug_assert_eq!(pivots.len(), self.pivots_len());
        debug_assert_eq!(carry.len(), self.carry_len());
        scratch.fill(T::zero());
        if self.is_dense() {
            let dense = Band::dense(n);
            // A multiple of `n` that is at least `lower`: it keeps the column
            // below from going under zero.
            let turns = n * lower.div_ceil(n.max(1));
            for (i, row) in rows.chunks_exact(width).enumerate() {
                for (e, entry) in row.iter().enumerate() {
                    let at = dense.at(i, (i + e + turns - lower) % n);
                    scratch[at] = scratch[at] + *entry;
                }
            }
            factor(dense, scratch, pivots)?;
            solve(dense, scratch, pivots, rhs);
            return Ok(());
        }

        let band = self.band();
        let rank = self.rank();
        let small = Band::dense(rank);
        let (ab, rest) = scratch.split_at_mut(band.len());
        let (z, h) = rest.split_at_mut(n * rank);
        let (pivots, small_pivots) = pivots.split_at_mut(n);
        // Column `a` of `U` multiplies the unknown `picked(a)`: the last `lower`
        // unknowns, which the top right corner reaches, then the first `upper`,
        // which the bottom left corner reaches.
        let picked = |a: usize| if a < lower { n - lower + a } else { a - lower };
        for (i, row) in rows.chunks_exact(width).enumerate() {
            for (e, entry) in row.iter().enumerate() {
                // The column is `i + e - lower`, before it wraps.
                if i + e < lower {
                    // It wraps to `n + i + e - lower`, in the top right corner.
                    z[(i + e) * n + i] = *entry;
                } else if i + e - lower >= n {
                    // It wraps to `i + e - lower - n`, in the bottom left corner.
                    z[(i + e - n) * n + i] = *entry;
                } else {
                    ab[band.at(i, i + e - lower)] = *entry;
                }
            }
        }
        factor(band, ab, pivots)?;
        for column in z.chunks_exact_mut(n) {
            solve(band, ab, pivots, column);
        }
        // I + V^T Z: row `b` is row `picked(b)` of `Z`.
        for a in 0..rank {
            for b in 0..rank {
                let identity = if a == b { T::one() } else { T::zero() };
                h[small.at(b, a)] = identity + z[a * n + picked(b)];
            }
        }
        factor(small, h, small_pivots).map_err(|at| Breakdown {
            row: picked(at.row),
        })?;
        solve(band, ab, pivots, rhs);
        for (b, slot) in carry.iter_mut().enumerate() {
            *slot = rhs[picked(b)];
        }
        solve(small, h, small_pivots, carry);
        for (column, weight) in z.chunks_exact(n).zip(&*carry) {
            for (slot, entry) in rhs.iter_mut().zip(column) {
                *slot = slot.axpy(-*entry, *weight);
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::{factor, solve, solve_slabs, Band, Breakdown, Cyclic};
    use crate::common::{assert_close, load, reals, Real, DIRECT_F64};
    use proptest::prelude::*;
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

    fn size(v: &Json) -> usize {
        usize::try_from(v.as_u64().expect("a count")).expect("fits")
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

    /// The storage of a banded matrix from SciPy's: `diagonals[upper + i - j][j]`
    /// is entry (`i`, `j`). The space above the band is filled with `garbage`.
    fn storage(band: Band, diagonals: &[Vec<f64>], garbage: f64) -> Vec<f64> {
        let mut ab = vec![garbage; band.len()];
        for j in 0..band.n {
            for i in j.saturating_sub(band.upper)..(j + band.lower + 1).min(band.n) {
                ab[band.at(i, j)] = diagonals[band.upper + i - j][j];
            }
        }
        ab
    }

    /// Solves a cyclic system with scratch space that holds `garbage`.
    fn solve_cyclic<const M: usize>(
        shape: Cyclic,
        matrix: &[f64],
        rhs: &mut [[f64; M]],
        garbage: f64,
    ) -> Result<(), Breakdown> {
        let mut scratch = vec![garbage; shape.scratch_len()];
        let mut pivots = vec![usize::MAX; shape.pivots_len()];
        let mut carry = vec![[garbage; M]; shape.carry_len()];
        shape.solve(matrix, rhs, &mut scratch, &mut pivots, &mut carry)
    }

    /// Solves one golden case with `M` columns, or returns `None` when the case
    /// has another number of columns.
    fn run<const M: usize>(case: &Json) -> Option<Result<Vec<[f64; M]>, Breakdown>> {
        let mut rhs = columns::<M>(&case["rhs"])?;
        let (n, lower, upper) = (
            rhs.len(),
            size(&case["options"]["lower"]),
            size(&case["options"]["upper"]),
        );
        let done = if is_cyclic(case) {
            let matrix: Vec<f64> = rows(&case["rows"]).concat();
            solve_cyclic(Cyclic { n, lower, upper }, &matrix, &mut rhs, f64::NAN)
        } else {
            let band = Band { n, lower, upper };
            let mut ab = storage(band, &rows(&case["band"]), f64::NAN);
            let mut pivots = vec![usize::MAX; n];
            factor(band, &mut ab, &mut pivots).map(|()| solve(band, &ab, &pivots, &mut rhs))
        };
        Some(done.map(|()| rhs))
    }

    /// Checks one golden case with `M` columns; `false` when it has another
    /// number of columns.
    fn check<const M: usize>(case: &Json) -> bool {
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
                assert_close(g[column], w[column], DIRECT_F64, &what);
            }
        }
        true
    }

    #[test]
    fn solutions_match_scipy() {
        let file = load::<Json>("banded").expect("golden file banded");
        let (mut plain, mut cyclic, mut wide, mut exchanged, mut wrapped) = (0, 0, 0, 0, 0);
        for case in &file.cases {
            let name = case["name"].as_str().expect("name");
            if check::<1>(case) {
                // One column also goes through the scalar value type, and gives
                // the same bits.
                let mut rhs = nums(&case["rhs"]);
                let (n, lower, upper) = (
                    rhs.len(),
                    size(&case["options"]["lower"]),
                    size(&case["options"]["upper"]),
                );
                if is_cyclic(case) {
                    let shape = Cyclic { n, lower, upper };
                    let matrix: Vec<f64> = rows(&case["rows"]).concat();
                    let mut scratch = vec![0.0; shape.scratch_len()];
                    let mut pivots = vec![0; shape.pivots_len()];
                    let mut carry = vec![0.0; shape.carry_len()];
                    shape
                        .solve(&matrix, &mut rhs, &mut scratch, &mut pivots, &mut carry)
                        .expect("solved");
                    cyclic += 1;
                    wrapped += usize::from(n <= lower + upper);
                } else {
                    let band = Band { n, lower, upper };
                    let mut ab = storage(band, &rows(&case["band"]), 0.0);
                    let mut pivots = vec![0; n];
                    factor(band, &mut ab, &mut pivots).expect("solved");
                    solve(band, &ab, &pivots, &mut rhs);
                    plain += 1;
                    // The generator says whether LAPACK exchanged rows, and
                    // this elimination chooses the same pivots.
                    let swaps = pivots.iter().enumerate().filter(|(j, p)| j != *p).count();
                    assert_eq!(
                        swaps,
                        size(&case["exchanges"]),
                        "{name}: the number of row exchanges"
                    );
                    exchanged += usize::from(swaps > 0);
                }
                let as_arrays = run::<1>(case).expect("one column").expect("solved");
                for (scalar, array) in rhs.iter().zip(&as_arrays) {
                    assert_eq!(scalar.to_bits(), array[0].to_bits(), "{name}");
                }
            } else {
                assert!(check::<2>(case) || check::<3>(case), "{name}");
                wide += 1;
            }
        }
        // With one column: 49 banded cases, in 25 of which rows are exchanged,
        // and 29 cyclic ones, 7 of them so small that the band wraps onto
        // itself. With more columns: 3.
        assert_eq!(
            (plain, cyclic, wide, exchanged, wrapped),
            (49, 29, 3, 25, 7)
        );
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

    /// A banded matrix as dense rows, with its shape and a right-hand side.
    type System = (Band, Vec<Vec<f64>>, Vec<f64>);

    /// A random matrix with the given bands. With `dominant`, its diagonal
    /// exceeds the sum of the other entries of its row, so that it is well
    /// conditioned; without, the pivots are wherever the entries are large, and
    /// rows are exchanged.
    fn systems(dominant: bool) -> impl Strategy<Value = System> {
        (1..24_usize, 0..5_usize, 0..5_usize)
            .prop_flat_map(move |(n, lower, upper)| {
                (
                    Just(Band { n, lower, upper }),
                    prop::collection::vec(prop::collection::vec(-1.0_f64..1.0, n), n),
                    prop::collection::vec((1.0_f64..2.0, any::<bool>()), n),
                    prop::collection::vec(-10.0_f64..10.0, n),
                )
            })
            .prop_map(move |(band, mut a, margins, rhs)| {
                for (i, row) in a.iter_mut().enumerate() {
                    for (j, entry) in row.iter_mut().enumerate() {
                        if j + band.lower < i || i + band.upper < j {
                            *entry = 0.0;
                        }
                    }
                    let (margin, negative) = margins[i];
                    let others: f64 = row.iter().map(|v| v.abs()).sum::<f64>() - row[i].abs();
                    if dominant {
                        row[i] = if negative { -1.0 } else { 1.0 } * (others + margin);
                    }
                }
                (band, a, rhs)
            })
    }

    /// The storage of a dense matrix of the given bands.
    fn pack(band: Band, a: &[Vec<f64>], garbage: f64) -> Vec<f64> {
        let mut ab = vec![garbage; band.len()];
        for j in 0..band.n {
            for i in j.saturating_sub(band.upper)..(j + band.lower + 1).min(band.n) {
                ab[band.at(i, j)] = a[i][j];
            }
        }
        ab
    }

    /// The largest residual of `a x = b`, over the largest entry of `b` and 1.
    fn residual(a: &[Vec<f64>], x: &[f64], b: &[f64]) -> f64 {
        let scale = b.iter().fold(1.0_f64, |m, v| m.max(v.abs()));
        let worst = a
            .iter()
            .zip(b)
            .map(|(row, b)| (row.iter().zip(x).map(|(a, x)| a * x).sum::<f64>() - b).abs())
            .fold(0.0, f64::max);
        worst / scale
    }

    proptest! {
        #[test]
        fn a_dominant_system_is_solved_as_the_dense_elimination_solves_it(
            (band, a, rhs) in systems(true),
        ) {
            let want = dense_solve(a.clone(), rhs.clone());
            let mut ab = pack(band, &a, f64::NAN);
            let (mut pivots, mut got) = (vec![usize::MAX; band.n], rhs.clone());
            prop_assert_eq!(factor(band, &mut ab, &mut pivots), Ok(()));
            solve(band, &ab, &pivots, &mut got);
            for (g, w) in got.iter().zip(&want) {
                prop_assert!((g - w).abs() <= 1e-11 * w.abs().max(1.0), "{} for {}", g, w);
            }
            prop_assert!(residual(&a, &got, &rhs) <= 1e-12);
            // A dominant diagonal is the pivot of every column.
            prop_assert!(pivots.iter().enumerate().all(|(j, p)| j == *p));
        }

        #[test]
        fn a_system_that_needs_exchanges_is_solved_with_the_pivots_of_the_dense_elimination(
            (band, a, rhs) in systems(false),
        ) {
            let mut ab = pack(band, &a, f64::NAN);
            let (mut pivots, mut got) = (vec![usize::MAX; band.n], rhs.clone());
            // A random matrix may be singular, or close to it: then there is
            // nothing to compare.
            prop_assume!(factor(band, &mut ab, &mut pivots).is_ok());
            solve(band, &ab, &pivots, &mut got);
            let size = got.iter().fold(1.0_f64, |m, v| m.max(v.abs()));
            prop_assume!(size < 1e6);
            // Backward stability: the residual is small against the solution.
            prop_assert!(residual(&a, &got, &rhs) <= 1e-11 * size, "residual {}", residual(&a, &got, &rhs));
            let want = dense_solve(a.clone(), rhs.clone());
            for (g, w) in got.iter().zip(&want) {
                prop_assert!((g - w).abs() <= 1e-6 * size, "{} for {}", g, w);
            }
            // Each pivot is the largest entry of its column on or below the
            // diagonal: no multiplier exceeds 1.
            for j in 0..band.n {
                let diag = band.at(j, j);
                for i in 1..=band.lower.min(band.n - 1 - j) {
                    prop_assert!(ab[diag + i].abs() <= 1.0, "column {} multiplier {}", j, ab[diag + i]);
                }
            }
        }

        #[test]
        fn one_factorization_solves_every_right_hand_side(
            (band, a, first) in systems(false),
            scale in -3.0_f64..3.0,
        ) {
            let mut ab = pack(band, &a, 0.0);
            let mut pivots = vec![0; band.n];
            prop_assume!(factor(band, &mut ab, &mut pivots).is_ok());
            let second: Vec<f64> = first.iter().rev().map(|v| v * scale + 1.0).collect();
            let (mut x, mut y) = (first.clone(), second.clone());
            solve(band, &ab, &pivots, &mut x);
            solve(band, &ab, &pivots, &mut y);
            // Two columns at once are the two solutions, bit for bit.
            let mut both: Vec<[f64; 2]> = first.iter().zip(&second).map(|(x, y)| [*x, *y]).collect();
            solve(band, &ab, &pivots, &mut both);
            for ((pair, x), y) in both.iter().zip(&x).zip(&y) {
                prop_assert_eq!(pair[0].to_bits(), x.to_bits());
                prop_assert_eq!(pair[1].to_bits(), y.to_bits());
            }
        }

        #[test]
        fn slabs_are_solved_as_if_each_position_were_alone(
            (band, a, first) in systems(false),
            width in 1..5_usize,
            scale in -3.0_f64..3.0,
        ) {
            let mut ab = pack(band, &a, 0.0);
            let mut pivots = vec![0; band.n];
            prop_assume!(factor(band, &mut ab, &mut pivots).is_ok());
            // Position `p` of every slab holds its own right-hand side.
            let column = |p: usize| -> Vec<f64> {
                first.iter().enumerate().map(|(i, v)| v * scale.powi(i32::try_from(p).unwrap_or(0)) + (p * i) as f64).collect()
            };
            let mut slabs = vec![0.0; band.n * width];
            for p in 0..width {
                for (i, v) in column(p).iter().enumerate() {
                    slabs[i * width + p] = *v;
                }
            }
            solve_slabs(band, &ab, &pivots, &mut slabs, width);
            for p in 0..width {
                let mut alone = column(p);
                solve(band, &ab, &pivots, &mut alone);
                for (i, want) in alone.iter().enumerate() {
                    prop_assert_eq!(slabs[i * width + p].to_bits(), want.to_bits(), "row {} position {}", i, p);
                }
            }
        }

        #[test]
        fn what_the_storage_holds_above_the_band_is_not_read(
            (band, a, rhs) in systems(false),
        ) {
            let run = |garbage: f64| {
                let mut ab = pack(band, &a, garbage);
                let (mut pivots, mut x) = (vec![0; band.n], rhs.clone());
                factor(band, &mut ab, &mut pivots).map(|()| {
                    solve(band, &ab, &pivots, &mut x);
                    x
                })
            };
            let (clean, dirty) = (run(0.0), run(f64::NAN));
            prop_assert_eq!(clean.is_ok(), dirty.is_ok());
            if let (Ok(clean), Ok(dirty)) = (clean, dirty) {
                for (c, d) in clean.iter().zip(&dirty) {
                    prop_assert_eq!(c.to_bits(), d.to_bits());
                }
            }
        }

        #[test]
        fn a_cyclic_system_is_solved_as_the_dense_elimination_solves_it(
            n in 1..24_usize,
            lower in 0..5_usize,
            upper in 0..5_usize,
            entries in prop::collection::vec(-1.0_f64..1.0, 24 * 9),
            margins in prop::collection::vec((1.0_f64..2.0, any::<bool>()), 24),
            rhs in prop::collection::vec(-10.0_f64..10.0, 24),
        ) {
            let width = lower + upper + 1;
            let mut matrix = entries[..n * width].to_vec();
            // The dense form: entries that wrap onto the same column add up.
            let mut a = vec![vec![0.0; n]; n];
            for i in 0..n {
                let row = &mut matrix[i * width..(i + 1) * width];
                let others: f64 = row.iter().map(|v| v.abs()).sum::<f64>() - row[lower].abs();
                let (margin, negative) = margins[i];
                row[lower] = if negative { -1.0 } else { 1.0 } * (others + margin);
                for (e, entry) in row.iter().enumerate() {
                    a[i][(i + e + 5 * n - lower) % n] += entry;
                }
            }
            let rhs = &rhs[..n];
            let want = dense_solve(a.clone(), rhs.to_vec());
            let shape = Cyclic { n, lower, upper };
            let mut got: Vec<[f64; 1]> = rhs.iter().map(|v| [*v]).collect();
            prop_assert_eq!(solve_cyclic(shape, &matrix, &mut got, f64::NAN), Ok(()));
            for (g, w) in got.iter().zip(&want) {
                prop_assert!((g[0] - w).abs() <= 1e-10 * w.abs().max(1.0), "{} for {}", g[0], w);
            }
            // Two columns are solved as if each were alone.
            let mut both: Vec<[f64; 2]> = rhs.iter().map(|v| [*v, 1.0 - 2.0 * v]).collect();
            let mut alone: Vec<[f64; 1]> = rhs.iter().map(|v| [1.0 - 2.0 * v]).collect();
            prop_assert_eq!(solve_cyclic(shape, &matrix, &mut both, 0.0), Ok(()));
            prop_assert_eq!(solve_cyclic(shape, &matrix, &mut alone, 0.0), Ok(()));
            for ((pair, x), y) in both.iter().zip(&got).zip(&alone) {
                prop_assert_eq!(pair[0].to_bits(), x[0].to_bits());
                prop_assert_eq!(pair[1].to_bits(), y[0].to_bits());
            }
        }
    }

    #[test]
    fn a_singular_matrix_is_a_breakdown_at_the_column_without_a_pivot() {
        // The second row is twice the first. Every step is exact: the
        // multipliers are 0.5, and the last pivot is a zero, not a rounding.
        let band = Band::dense(3);
        let a = [
            vec![1.0, 2.0, 3.0],
            vec![2.0, 4.0, 6.0],
            vec![1.0, 1.0, 1.0],
        ];
        let mut ab = pack(band, &a, 0.0);
        assert_eq!(
            factor(band, &mut ab, &mut [0; 3]),
            Err(Breakdown { row: 2 })
        );
        // A column of zeros.
        let band = Band {
            n: 3,
            lower: 1,
            upper: 1,
        };
        let a = [
            vec![2.0, 0.0, 0.0],
            vec![1.0, 0.0, 1.0],
            vec![0.0, 0.0, 3.0],
        ];
        let mut ab = pack(band, &a, 0.0);
        assert_eq!(
            factor(band, &mut ab, &mut [0; 3]),
            Err(Breakdown { row: 1 })
        );
        // No rows at all: nothing to do.
        let empty = Band::dense(0);
        assert_eq!(empty.len(), 0);
        assert_eq!(factor::<f64>(empty, &mut [], &mut []), Ok(()));
    }

    #[test]
    fn a_pivot_that_is_not_finite_is_a_breakdown() {
        let band = Band {
            n: 2,
            lower: 1,
            upper: 1,
        };
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut ab = pack(band, &[vec![bad, 1.0], vec![0.5, 2.0]], 0.0);
            assert_eq!(
                factor(band, &mut ab, &mut [0; 2]),
                Err(Breakdown { row: 0 }),
                "{bad}"
            );
        }
        // In the second column, after the elimination: 1e308 - (-1e308).
        let mut ab = pack(band, &[vec![1.0, -1e308], vec![1.0, 1e308]], 0.0);
        assert_eq!(
            factor(band, &mut ab, &mut [0; 2]),
            Err(Breakdown { row: 1 })
        );
    }

    #[test]
    fn a_zero_on_the_diagonal_is_passed_by_an_exchange() {
        // [[0, 1], [1, 0]] x = [3, 5].
        let band = Band {
            n: 2,
            lower: 1,
            upper: 1,
        };
        let mut ab = pack(band, &[vec![0.0, 1.0], vec![1.0, 0.0]], 0.0);
        let mut pivots = [0; 2];
        assert_eq!(factor(band, &mut ab, &mut pivots), Ok(()));
        assert_eq!(pivots, [1, 1]);
        let mut x = [3.0, 5.0];
        solve(band, &ab, &pivots, &mut x);
        assert_eq!(x.map(f64::to_bits), [5.0_f64, 3.0].map(f64::to_bits));
    }

    #[test]
    fn of_entries_equal_in_size_the_first_is_the_pivot() {
        // As LAPACK chooses (`idamax` returns the first of the largest). Both
        // choices solve the system; the solutions differ in rounding.
        let band = Band::dense(3);
        let a = [
            vec![2.0, 1.0, 0.0],
            vec![-2.0, 3.0, 1.0],
            vec![2.0, 0.0, 5.0],
        ];
        let mut ab = pack(band, &a, 0.0);
        let mut pivots = [usize::MAX; 3];
        assert_eq!(factor(band, &mut ab, &mut pivots), Ok(()));
        // Column 0: three entries of size 2, the first stays. Column 1: 4 and
        // -1 are left, the first stays.
        assert_eq!(pivots, [0, 1, 2]);
    }

    #[test]
    fn a_singular_cyclic_matrix_is_a_breakdown() {
        // The identity with 2 in the top right corner and 0.5 in the bottom
        // left one: the first and the last row are multiples of each other.
        // The band alone is the identity, so the breakdown is in the
        // correction, and it is exact: I + V^T Z is [[1, 0.5], [2, 1]]. The
        // unknown without a pivot is the first.
        let n = 6;
        let shape = Cyclic {
            n,
            lower: 1,
            upper: 1,
        };
        let mut matrix: Vec<f64> = (0..n).flat_map(|_| [0.0, 1.0, 0.0]).collect();
        matrix[0] = 2.0;
        matrix[3 * n - 1] = 0.5;
        let mut rhs = vec![[1.0]; n];
        let done = solve_cyclic(shape, &matrix, &mut rhs, 0.0);
        assert_eq!(done, Err(Breakdown { row: 0 }));
        // A breakdown in the band: its second row is zero.
        let mut matrix: Vec<f64> = (0..n).flat_map(|_| [0.0, 1.0, 0.0]).collect();
        matrix[4] = 0.0;
        let mut rhs = vec![[1.0]; n];
        let done = solve_cyclic(shape, &matrix, &mut rhs, 0.0);
        assert_eq!(done, Err(Breakdown { row: 1 }));
        // A matrix so small that the band wraps onto itself: both rows are
        // [2, -2] and [-2, 2].
        let shape = Cyclic {
            n: 2,
            lower: 1,
            upper: 1,
        };
        let mut rhs = vec![[1.0]; 2];
        let done = solve_cyclic(shape, &[-1.0, 2.0, -1.0, -1.0, 2.0, -1.0], &mut rhs, 0.0);
        assert_eq!(done, Err(Breakdown { row: 1 }));
    }

    #[test]
    fn the_lengths_of_the_cyclic_scratch_space_follow_the_shape() {
        // Wider than the matrix: a dense matrix of 3 rows, 3 * (2 * 2 + 2 + 1).
        let small = Cyclic {
            n: 3,
            lower: 2,
            upper: 2,
        };
        assert_eq!(
            (small.scratch_len(), small.pivots_len(), small.carry_len()),
            (21, 3, 0)
        );
        // The band (10 * 4), Z (10 * 2), and a dense matrix of 2 rows (2 * 4).
        let large = Cyclic {
            n: 10,
            lower: 1,
            upper: 1,
        };
        assert_eq!(
            (large.scratch_len(), large.pivots_len(), large.carry_len()),
            (68, 12, 2)
        );
        // Where the two meet. As many rows as diagonals beside the main one:
        // the corners would touch the band, so a dense matrix, 4 * (2 * 3 + 3 + 1).
        let last_dense = Cyclic {
            n: 4,
            lower: 2,
            upper: 2,
        };
        assert_eq!(
            (
                last_dense.scratch_len(),
                last_dense.pivots_len(),
                last_dense.carry_len()
            ),
            (40, 4, 0)
        );
        // One row more: the band (5 * 7), Z (5 * 4), and a dense matrix of 4
        // rows (4 * 10).
        let first_banded = Cyclic {
            n: 5,
            lower: 2,
            upper: 2,
        };
        assert_eq!(
            (
                first_banded.scratch_len(),
                first_banded.pivots_len(),
                first_banded.carry_len()
            ),
            (95, 9, 4)
        );
    }
}
