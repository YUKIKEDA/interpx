"""Generate `goldens/tridiagonal.json`: solutions of tridiagonal and cyclic systems.

Reference: SciPy 1.18. A tridiagonal system is solved by
`scipy.linalg.solve_banded` (LAPACK `gbsv`, with partial pivoting), which is what
`CubicSpline` calls. A cyclic system (one entry in each far corner) is solved as
a dense matrix by `scipy.linalg.solve`. Design: docs/design.ja.md section 6.1.

A matrix of `n` rows is stored as three lists of length `n`: `lower[i]` is the
entry left of the diagonal in row `i`, `diag[i]` the one on it, `upper[i]` the
one right of it. For a cyclic matrix `lower[0]` is the top right corner and
`upper[n - 1]` the bottom left one; otherwise both are zero. `rhs` and `expected`
are lists, or lists of rows when the right-hand side has several columns.

The spline matrices are built the way `scipy/interpolate/_cubic.py` builds them,
so that the solver is tested on the systems it is for.

Run: `just gen-goldens`
"""

import numpy as np
from scipy import linalg

from common import classic_data, golden

GENERATOR = "generate_tridiagonal.py"
SEED = 6


def solve(lower: np.ndarray, diag: np.ndarray, upper: np.ndarray, rhs: np.ndarray) -> np.ndarray:
    """The banded solve `CubicSpline` uses."""
    n = len(diag)
    ab = np.zeros((3, n))
    ab[0, 1:] = upper[:-1]
    ab[1] = diag
    ab[2, :-1] = lower[1:]
    return linalg.solve_banded((1, 1), ab, rhs)


def solve_cyclic(lower: np.ndarray, diag: np.ndarray, upper: np.ndarray, rhs: np.ndarray) -> np.ndarray:
    """A dense solve of the cyclic matrix."""
    n = len(diag)
    a = np.diag(diag) + np.diag(upper[:-1], 1) + np.diag(lower[1:], -1)
    a[0, n - 1] = lower[0]
    a[n - 1, 0] = upper[n - 1]
    return linalg.solve(a, rhs)


def case(name: str, lower, diag, upper, rhs, cyclic: bool = False) -> dict:
    lower, diag, upper, rhs = (np.asarray(v, dtype=np.float64) for v in (lower, diag, upper, rhs))
    if cyclic:
        reference, expected = "scipy.linalg.solve", solve_cyclic(lower, diag, upper, rhs)
    else:
        reference, expected = "scipy.linalg.solve_banded", solve(lower, diag, upper, rhs)
    return golden.case(
        f"tridiagonal/{name}",
        reference,
        {"cyclic": cyclic},
        lower=lower,
        diag=diag,
        upper=upper,
        rhs=rhs,
        expected=expected,
    )


def dominant(rng: np.random.Generator, n: int, cyclic: bool, margin: float = 1.0):
    """A random matrix whose diagonal exceeds the two entries beside it by `margin`."""
    lower = rng.uniform(-1.0, 1.0, n)
    upper = rng.uniform(-1.0, 1.0, n)
    if not cyclic:
        lower[0] = 0.0
        upper[-1] = 0.0
    sign = rng.choice([-1.0, 1.0], n)
    diag = sign * (np.abs(lower) + np.abs(upper) + margin * rng.uniform(0.5, 1.5, n))
    return lower, diag, upper


def spline_system(x: np.ndarray, y: np.ndarray, start, end):
    """The system for the slopes of a cubic spline, as `_cubic.py` sets it up.

    `start` and `end` are "not-a-knot", or (order, value) for a derivative of
    order 1 or 2.
    """
    n = len(x)
    dx = np.diff(x)
    slope = np.diff(y) / dx
    lower, diag, upper, b = np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n)
    diag[1:-1] = 2.0 * (dx[:-1] + dx[1:])
    upper[1:-1] = dx[:-1]
    lower[1:-1] = dx[1:]
    b[1:-1] = 3.0 * (dx[1:] * slope[:-1] + dx[:-1] * slope[1:])

    if start == "not-a-knot":
        d = x[2] - x[0]
        diag[0], upper[0] = dx[1], d
        b[0] = ((dx[0] + 2.0 * d) * dx[1] * slope[0] + dx[0] ** 2 * slope[1]) / d
    elif start[0] == 1:
        diag[0], upper[0], b[0] = 1.0, 0.0, start[1]
    else:
        diag[0], upper[0] = 2.0 * dx[0], dx[0]
        b[0] = -0.5 * start[1] * dx[0] ** 2 + 3.0 * (y[1] - y[0])

    if end == "not-a-knot":
        d = x[-1] - x[-3]
        diag[-1], lower[-1] = dx[-2], d
        b[-1] = (dx[-1] ** 2 * slope[-2] + (2.0 * d + dx[-1]) * dx[-2] * slope[-1]) / d
    elif end[0] == 1:
        diag[-1], lower[-1], b[-1] = 1.0, 0.0, end[1]
    else:
        diag[-1], lower[-1] = 2.0 * dx[-1], dx[-1]
        b[-1] = 0.5 * end[1] * dx[-1] ** 2 + 3.0 * (y[-1] - y[-2])
    return lower, diag, upper, b


def periodic_system(x: np.ndarray, y: np.ndarray):
    """The cyclic system for the slopes of a periodic cubic spline at all points but the last."""
    dx = np.diff(x)
    slope = np.diff(y) / dx
    before, slope_before = np.roll(dx, 1), np.roll(slope, 1)
    lower, diag, upper = dx, 2.0 * (before + dx), before
    b = 3.0 * (dx * slope_before + before * slope)
    return lower, diag, upper, b


BOUNDARIES = {
    "not-a-knot": ("not-a-knot", "not-a-knot"),
    "natural": ((2, 0.0), (2, 0.0)),
    "clamped": ((1, 0.0), (1, 0.0)),
    "first_and_second": ((1, 0.5), (2, -0.25)),
}


def main() -> None:
    rng = np.random.default_rng(SEED)
    cases = []

    # Small systems, where each row is a special case of the elimination.
    cases.append(case("small/n1", [0.0], [4.0], [0.0], [2.0]))
    cases.append(case("small/n2", [0.0, 1.0], [4.0, 3.0], [-2.0, 0.0], [2.0, 5.0]))
    cases.append(case("small/n3", [0.0, 1.0, -1.0], [4.0, 3.0, 5.0], [-2.0, 0.5, 0.0], [2.0, 5.0, -1.0]))
    cases.append(case("small/n3/cyclic", [0.75, 1.0, -1.0], [4.0, 3.0, 5.0], [-2.0, 0.5, 1.25], [2.0, 5.0, -1.0], cyclic=True))
    cases.append(case("small/identity", [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [0.0, 0.0, 0.0], [3.0, -1.0, 2.0]))

    # Random diagonally dominant matrices.
    for n in (4, 10, 100, 400):
        lower, diag, upper = dominant(rng, n, cyclic=False)
        cases.append(case(f"dominant/n{n}", lower, diag, upper, rng.uniform(-1.0, 1.0, n)))
        lower, diag, upper = dominant(rng, n, cyclic=True)
        cases.append(case(f"dominant/n{n}/cyclic", lower, diag, upper, rng.uniform(-1.0, 1.0, n), cyclic=True))
    # Barely dominant: the diagonal exceeds its neighbors by about 1e-6.
    lower, diag, upper = dominant(rng, 50, cyclic=False, margin=1e-6)
    cases.append(case("weakly_dominant/n50", lower, diag, upper, rng.uniform(-1.0, 1.0, 50)))
    lower, diag, upper = dominant(rng, 50, cyclic=True, margin=1e-6)
    cases.append(case("weakly_dominant/n50/cyclic", lower, diag, upper, rng.uniform(-1.0, 1.0, 50), cyclic=True))

    # Right-hand sides with several columns.
    lower, diag, upper = dominant(rng, 12, cyclic=False)
    cases.append(case("columns/n12/3", lower, diag, upper, rng.uniform(-1.0, 1.0, (12, 3))))
    lower, diag, upper = dominant(rng, 12, cyclic=True)
    cases.append(case("columns/n12/2/cyclic", lower, diag, upper, rng.uniform(-1.0, 1.0, (12, 2)), cyclic=True))

    # The systems of a cubic spline on the classic data.
    data = {}
    for name, pairs in (("fritsch_carlson", classic_data.FRITSCH_CARLSON), ("titanium", classic_data.TITANIUM)):
        a = np.array(pairs)
        data[name] = (a[:, 0], a[:, 1])
    x = np.linspace(-1.0, 1.0, 11)
    data["runge_11"] = (x, 1.0 / (1.0 + 25.0 * x**2))
    x = np.sort(rng.uniform(-5.0, 5.0, 20))
    data["random_20"] = (x, np.sin(x) + 0.1 * rng.normal(size=20))
    # Spacings that alternate between 1 and 1e-6: far from diagonally dominant in
    # the not-a-knot rows.
    x = np.cumsum(np.tile([1.0, 1e-6], 10))
    data["uneven_20"] = (x, np.cos(x))
    for name, (x, y) in data.items():
        for label, (start, end) in BOUNDARIES.items():
            cases.append(case(f"spline/{label}/{name}", *spline_system(x, y, start, end)))
    # The fewest points the general set-up takes: 4 for not-a-knot at both ends.
    x, y = np.array([0.0, 1.0, 1.5, 4.0]), np.array([1.0, -1.0, 2.0, 0.5])
    for label, (start, end) in BOUNDARIES.items():
        cases.append(case(f"spline/{label}/four_points", *spline_system(x, y, start, end)))

    # The cyclic systems of a periodic cubic spline.
    x = np.linspace(0.0, 2.0 * np.pi, 13)
    y = np.sin(x) + 0.5 * np.cos(2.0 * x)
    cases.append(case("spline/periodic/even_13", *periodic_system(x, y), cyclic=True))
    x = np.concatenate([[0.0], np.cumsum(rng.uniform(0.05, 1.0, 30))])
    y = rng.normal(size=31)
    y[-1] = y[0]
    cases.append(case("spline/periodic/random_31", *periodic_system(x, y), cyclic=True))
    x = np.concatenate([[0.0], np.cumsum(np.tile([1.0, 1e-6], 8))])
    y = np.cos(2.0 * np.pi * x / x[-1])
    cases.append(case("spline/periodic/uneven_17", *periodic_system(x, y), cyclic=True))
    x, y = np.array([0.0, 1.0, 1.5, 4.0]), np.array([1.0, -1.0, 2.0, 1.0])
    cases.append(case("spline/periodic/four_points", *periodic_system(x, y), cyclic=True))

    path = golden.write("tridiagonal", GENERATOR, SEED, cases)
    print(f"wrote {path.name} ({len(cases)} cases)")


if __name__ == "__main__":
    main()
