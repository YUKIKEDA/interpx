"""Generate `goldens/banded.json`: solutions of banded and cyclic banded systems.

Reference: SciPy 1.18. A banded system is solved by LAPACK's `gbsv` (Gaussian
elimination with partial pivoting) through `scipy.linalg.lapack.dgbsv`, which
is what `make_interp_spline` calls. (`scipy.linalg.solve_banded` calls it too,
but not for one diagonal on each side: there it takes `gtsv`.) A cyclic banded
system (a band that wraps around the corners) is solved as a dense matrix by
`scipy.linalg.solve`. Design: docs/design.ja.md section 6.1.

A banded case has `lower` and `upper` in its options (the number of diagonals
below and above the main one) and the matrix in `band`, in the form of
`solve_banded`: `band[upper + i - j][j]` is entry (i, j). `exchanges` is the
number of columns in which LAPACK exchanged two rows.

A cyclic case has the matrix in `rows`, one row of `lower + upper + 1` entries
per row of the matrix: entry `e` of row `i` is in column `i + e - lower` modulo
`n`. When the band is wider than the matrix, entries that share a column add up.

`rhs` and `expected` are lists, or lists of rows when the right-hand side has
several columns.

The spline matrices are the collocation matrices `make_interp_spline` solves, so
that the solver is tested on the systems it is for.

Run: `just gen-goldens`
"""

import numpy as np
import scipy.interpolate as si
from scipy import linalg
from scipy.interpolate._bsplines import _not_a_knot, _periodic_knots

from common import classic_data, golden

GENERATOR = "generate_banded.py"
SEED = 12


def to_band(a: np.ndarray, lower: int, upper: int) -> np.ndarray:
    """The diagonals of a dense matrix, in the form of `solve_banded`."""
    n = a.shape[0]
    band = np.zeros((lower + upper + 1, n))
    for i in range(n):
        for j in range(max(0, i - lower), min(n, i + upper + 1)):
            band[upper + i - j, j] = a[i, j]
    outside = a - np.triu(np.tril(a, upper), -lower)
    assert not outside.any(), "entries outside the band"
    return band


def banded_case(name: str, a: np.ndarray, lower: int, upper: int, rhs, exchanged=None) -> dict:
    a, rhs = np.asarray(a, dtype=np.float64), np.asarray(rhs, dtype=np.float64)
    band = to_band(a, lower, upper)
    # The storage of LAPACK: `lower` rows above the band for the fill.
    storage = np.zeros((2 * lower + upper + 1, a.shape[0]))
    storage[lower:] = band
    # The pivots LAPACK chose: `piv[j]` is the row exchanged with row `j`.
    _, piv, expected, info = linalg.lapack.dgbsv(lower, upper, storage, rhs)
    assert info == 0, name
    assert np.allclose(expected, linalg.solve(a, rhs), rtol=1e-6, atol=1e-9), name
    exchanges = int(np.count_nonzero(piv != np.arange(len(piv))))
    if exchanged is not None:
        assert (exchanges > 0) == exchanged, (name, exchanges)
    return golden.case(
        f"banded/{name}",
        "scipy.linalg.lapack.dgbsv",
        {"cyclic": False, "lower": lower, "upper": upper},
        band=band,
        rhs=rhs,
        exchanges=exchanges,
        expected=expected,
    )


def cyclic_dense(rows: np.ndarray, lower: int) -> np.ndarray:
    n, width = rows.shape
    a = np.zeros((n, n))
    for i in range(n):
        for e in range(width):
            a[i, (i + e - lower) % n] += rows[i, e]
    return a


def cyclic_case(name: str, rows, lower: int, upper: int, rhs) -> dict:
    rows, rhs = np.asarray(rows, dtype=np.float64), np.asarray(rhs, dtype=np.float64)
    assert rows.shape[1] == lower + upper + 1
    return golden.case(
        f"banded/{name}",
        "scipy.linalg.solve",
        {"cyclic": True, "lower": lower, "upper": upper},
        rows=rows,
        rhs=rhs,
        expected=linalg.solve(cyclic_dense(rows, lower), rhs),
    )


def random_band(rng: np.random.Generator, n: int, lower: int, upper: int, margin: float = 1.0) -> np.ndarray:
    """A random banded matrix whose main diagonal exceeds the sum of the other
    entries of its row by about `margin`: it is well conditioned, and its
    elimination exchanges no rows."""
    a = np.triu(np.tril(rng.uniform(-1.0, 1.0, (n, n)), upper), -lower)
    others = np.abs(a).sum(axis=1) - np.abs(np.diag(a))
    a[np.diag_indices(n)] = rng.choice([-1.0, 1.0], n) * (others + margin * rng.uniform(0.5, 1.5, n))
    return a


def shuffled_band(rng: np.random.Generator, n: int, lower: int, upper: int) -> np.ndarray:
    """A dominant matrix of `lower - 1` and `upper - 1` diagonals with about
    half of its rows exchanged with a neighbor. The elimination has to exchange
    them back, which fills the space above the band; the condition number is
    that of the dominant matrix. (A random matrix with a small diagonal needs
    exchanges too, but is ill conditioned: two correct solvers then differ by
    more than the tolerance of the test.)"""
    a = random_band(rng, n, lower - 1, upper - 1)
    i = 0
    while i + 1 < n:
        if rng.random() < 0.6:
            a[[i, i + 1]] = a[[i + 1, i]]
            i += 2
        else:
            i += 1
    assert np.linalg.cond(a) < 100.0
    return a


def random_cyclic(rng: np.random.Generator, n: int, lower: int, upper: int, margin: float = 1.0) -> np.ndarray:
    rows = rng.uniform(-1.0, 1.0, (n, lower + upper + 1))
    others = np.abs(rows).sum(axis=1) - np.abs(rows[:, lower])
    rows[:, lower] = rng.choice([-1.0, 1.0], n) * (others + margin * rng.uniform(0.5, 1.5, n))
    return rows


def basis_values(t: np.ndarray, k: int, x: np.ndarray) -> np.ndarray:
    """The values of every B-spline of the knots `t` at the points `x`."""
    count = len(t) - k - 1
    out = np.zeros((len(x), count))
    for j in range(count):
        c = np.zeros(count)
        c[j] = 1.0
        out[:, j] = si.BSpline(t, c, k, extrapolate=True)(x)
    return out


def collocation(x: np.ndarray, k: int) -> np.ndarray:
    """The matrix `make_interp_spline` solves for the not-a-knot spline of degree `k`."""
    return basis_values(_not_a_knot(x, k), k, x)


def periodic_rows(x: np.ndarray, k: int) -> tuple[np.ndarray, int, int]:
    """The cyclic matrix of the periodic spline of degree `k`: the control
    points repeat with period `n - 1`, and the conditions are those at all
    points but the last."""
    n = len(x)
    m = n - 1
    lower, upper = k // 2, k - k // 2
    values = basis_values(_periodic_knots(x, k), k, x[:-1])
    rows = np.zeros((m, k + 1))
    for i in range(m):
        # The splines that are not zero at `x[i]` are `i` to `i + k`; unknown
        # `r` is the control point `r + lower`, modulo `m`.
        rows[i] = values[i, i : i + k + 1]
        assert np.count_nonzero(values[i]) <= k + 1
        assert not values[i, :i].any() and not values[i, i + k + 1 :].any()
    return rows, lower, upper


def main() -> None:
    rng = np.random.default_rng(SEED)
    cases = []

    # Small systems, where each column is a special case of the elimination.
    cases.append(banded_case("small/n1", [[4.0]], 0, 0, [2.0]))
    cases.append(banded_case("small/n1/wide_band", [[4.0]], 2, 1, [2.0]))
    cases.append(banded_case("small/n2", [[4.0, -2.0], [1.0, 3.0]], 1, 1, [2.0, 5.0]))
    cases.append(banded_case("small/n2/exchange", [[1.0, -2.0], [4.0, 3.0]], 1, 1, [2.0, 5.0], exchanged=True))
    a = [[4.0, -2.0, 0.0], [1.0, 3.0, 0.5], [0.0, -1.0, 5.0]]
    cases.append(banded_case("small/n3", a, 1, 1, [2.0, 5.0, -1.0]))
    a = [[0.0, 2.0, 1.0], [1.0, 0.0, 3.0], [4.0, 1.0, 0.0]]
    cases.append(banded_case("small/n3/zero_diagonal", a, 2, 2, [2.0, 5.0, -1.0], exchanged=True))
    cases.append(banded_case("small/identity", np.eye(3), 1, 2, [3.0, -1.0, 2.0], exchanged=False))
    cases.append(banded_case("small/diagonal", np.diag([2.0, -4.0, 0.5, 8.0]), 0, 0, [3.0, -1.0, 2.0, 4.0]))

    # Random matrices with a dominant diagonal: no row is exchanged.
    for n, lower, upper in ((5, 1, 1), (20, 2, 3), (20, 0, 2), (20, 2, 0), (60, 3, 3), (80, 5, 5), (80, 1, 4)):
        a = random_band(rng, n, lower, upper)
        name = f"dominant/n{n}/l{lower}u{upper}"
        cases.append(banded_case(name, a, lower, upper, rng.uniform(-1.0, 1.0, n), exchanged=False))

    # Dominant matrices with rows exchanged: the pivots are off the diagonal, and
    # the exchanges of the elimination fill the space above the band.
    for n, lower, upper in ((10, 1, 1), (10, 3, 1), (40, 2, 2), (40, 1, 3), (60, 3, 5), (60, 5, 2)):
        a = shuffled_band(rng, n, lower, upper)
        name = f"exchanges/n{n}/l{lower}u{upper}"
        cases.append(banded_case(name, a, lower, upper, rng.uniform(-1.0, 1.0, n), exchanged=True))
    # A zero diagonal: ones beside it.
    n = 8
    a = np.diag(np.ones(n - 1), 1) + np.diag(np.arange(2.0, n + 1.0), -1)
    cases.append(banded_case("exchanges/zero_diagonal/n8", a, 1, 1, np.arange(1.0, n + 1.0), exchanged=True))

    # The whole matrix in the band: a dense solve.
    for n in (2, 6, 12):
        a = shuffled_band(rng, n, n - 1, n - 1) + 0.05 * rng.uniform(-1.0, 1.0, (n, n))
        assert np.linalg.cond(a) < 100.0
        cases.append(banded_case(f"dense/n{n}", a, n - 1, n - 1, rng.uniform(-1.0, 1.0, n)))
    # More diagonals declared than the matrix has.
    a = shuffled_band(rng, 4, 3, 3)
    cases.append(banded_case("dense/n4/band_wider_than_the_matrix", a, 5, 6, rng.uniform(-1.0, 1.0, 4)))

    # Right-hand sides with several columns.
    a = shuffled_band(rng, 12, 2, 2)
    cases.append(banded_case("columns/n12/3", a, 2, 2, rng.uniform(-1.0, 1.0, (12, 3))))
    a = random_band(rng, 12, 3, 1)
    cases.append(banded_case("columns/n12/2", a, 3, 1, rng.uniform(-1.0, 1.0, (12, 2))))

    # The collocation matrices of interpolating splines on the classic data. The
    # band is declared as `make_interp_spline` declares it, `k` on each side.
    data = {}
    for name, pairs in (("fritsch_carlson", classic_data.FRITSCH_CARLSON), ("titanium", classic_data.TITANIUM)):
        pairs = np.array(pairs)
        data[name] = (pairs[:, 0], pairs[:, 1])
    x = np.sort(rng.uniform(-5.0, 5.0, 20))
    data["random_20"] = (x, np.sin(x) + 0.1 * rng.normal(size=20))
    # Spacings that alternate between 1 and 1e-3.
    x = np.cumsum(np.tile([1.0, 1e-3], 10))
    data["uneven_20"] = (x, np.cos(x))
    for name, (x, y) in data.items():
        for k in (2, 3, 4, 5, 7):
            cases.append(banded_case(f"spline/k{k}/{name}", collocation(x, k), k, k, y))
    # The fewest points: one polynomial of degree `k`.
    for k in (2, 3, 5):
        x = np.cumsum(rng.uniform(0.2, 1.0, k + 1))
        cases.append(banded_case(f"spline/k{k}/fewest_points", collocation(x, k), k, k, rng.normal(size=k + 1)))

    # Cyclic matrices. Small ones first: the band wraps onto itself and the
    # matrix is solved as a dense one.
    cases.append(cyclic_case("cyclic/small/n1", [[1.0, 4.0, -2.0]], 1, 1, [2.0]))
    cases.append(cyclic_case("cyclic/small/n2", [[1.0, 4.0, -2.0], [0.5, -3.0, 1.0]], 1, 1, [2.0, 5.0]))
    cases.append(cyclic_case("cyclic/small/n3/l2u2", random_cyclic(rng, 3, 2, 2), 2, 2, rng.uniform(-1.0, 1.0, 3)))
    cases.append(cyclic_case("cyclic/small/n2/l3u2", random_cyclic(rng, 2, 3, 2), 3, 2, rng.uniform(-1.0, 1.0, 2)))
    # As many rows as the band has diagonals beside the main one: the last
    # dense one. One row more: the first with corners.
    cases.append(cyclic_case("cyclic/small/n4/l2u2", random_cyclic(rng, 4, 2, 2), 2, 2, rng.uniform(-1.0, 1.0, 4)))
    cases.append(cyclic_case("cyclic/small/n5/l2u2", random_cyclic(rng, 5, 2, 2), 2, 2, rng.uniform(-1.0, 1.0, 5)))
    cases.append(cyclic_case("cyclic/small/n3/l1u1", random_cyclic(rng, 3, 1, 1), 1, 1, rng.uniform(-1.0, 1.0, 3)))
    # No band at all: a diagonal matrix.
    cases.append(cyclic_case("cyclic/small/n3/l0u0", [[2.0], [-4.0], [0.5]], 0, 0, [3.0, -1.0, 2.0]))

    for n, lower, upper in ((8, 1, 1), (30, 2, 2), (30, 1, 2), (30, 3, 0), (30, 0, 2), (60, 3, 3), (80, 2, 3)):
        rows = random_cyclic(rng, n, lower, upper)
        cases.append(cyclic_case(f"cyclic/dominant/n{n}/l{lower}u{upper}", rows, lower, upper, rng.uniform(-1.0, 1.0, n)))
    rows = random_cyclic(rng, 40, 2, 2, margin=1e-6)
    cases.append(cyclic_case("cyclic/weakly_dominant/n40", rows, 2, 2, rng.uniform(-1.0, 1.0, 40)))
    rows = random_cyclic(rng, 12, 1, 2)
    cases.append(cyclic_case("cyclic/columns/n12/2", rows, 1, 2, rng.uniform(-1.0, 1.0, (12, 2))))

    # The cyclic matrices of periodic splines.
    x = np.linspace(0.0, 2.0 * np.pi, 13)
    y = np.sin(x) + 0.5 * np.cos(2.0 * x)
    for k in (2, 3, 4, 5, 7):
        rows, lower, upper = periodic_rows(x, k)
        cases.append(cyclic_case(f"cyclic/spline/k{k}/even_13", rows, lower, upper, y[:-1]))
    x = np.concatenate([[0.0], np.cumsum(rng.uniform(0.05, 1.0, 30))])
    y = rng.normal(size=31)
    for k in (2, 3, 5):
        rows, lower, upper = periodic_rows(x, k)
        cases.append(cyclic_case(f"cyclic/spline/k{k}/random_31", rows, lower, upper, y[:-1]))
    x = np.concatenate([[0.0], np.cumsum(np.tile([1.0, 1e-3], 8))])
    y = np.cos(2.0 * np.pi * x / x[-1])
    for k in (2, 3, 4):
        rows, lower, upper = periodic_rows(x, k)
        cases.append(cyclic_case(f"cyclic/spline/k{k}/uneven_17", rows, lower, upper, y[:-1]))
    # Fewer points than the degree: the band wraps onto itself.
    x = np.array([0.0, 1.0, 1.5, 4.0])
    for k in (3, 5):
        rows, lower, upper = periodic_rows(x, k)
        cases.append(cyclic_case(f"cyclic/spline/k{k}/four_points", rows, lower, upper, [1.0, -1.0, 2.0]))

    path = golden.write("banded", GENERATOR, SEED, cases)
    print(f"wrote {path.name} ({len(cases)} cases)")


if __name__ == "__main__":
    main()
