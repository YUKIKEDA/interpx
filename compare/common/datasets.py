"""Data set definitions of docs/design.ja.md section 8.4.

Functions are vectorised NumPy. Samplers take a seed and return the same points on
every machine. The classic data is in `classic_data.py`.
"""

import numpy as np
from scipy.stats import qmc

# --- Franke's six functions on [0, 1]^2 (section 8.4.1) ---------------------------


def franke(k: int, x, y):
    """Franke test function `k` (1 to 6). F6 is sqrt(64 - 81 r^2) / 9 - 1/2."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    r2 = (x - 0.5) ** 2 + (y - 0.5) ** 2
    if k == 1:
        return (
            0.75 * np.exp(-((9 * x - 2) ** 2 + (9 * y - 2) ** 2) / 4)
            + 0.75 * np.exp(-((9 * x + 1) ** 2) / 49 - ((9 * y + 1) ** 2) / 10)
            + 0.5 * np.exp(-((9 * x - 7) ** 2 + (9 * y - 3) ** 2) / 4)
            - 0.2 * np.exp(-((9 * x - 4) ** 2) - (9 * y - 7) ** 2)
        )
    if k == 2:
        return (np.tanh(9 * (y - x)) + 1) / 9
    if k == 3:
        return (1.25 + np.cos(5.4 * y)) / (6 + 6 * (3 * x - 1) ** 2)
    if k == 4:
        return np.exp(-81 / 16 * r2) / 3
    if k == 5:
        return np.exp(-81 / 4 * r2) / 3
    if k == 6:
        return np.sqrt(64 - 81 * r2) / 9 - 0.5
    raise ValueError(f"no Franke function {k}")


# --- Functions for D >= 3 (section 8.4.4) -----------------------------------------

GENZ_KINDS = ("oscillatory", "product_peak", "corner_peak", "gaussian_peak", "continuous")


def genz_a(d: int) -> np.ndarray:
    """The parameters a_i: an implementation choice, 2 + i/d, not Genz's difficulty scaling."""
    return 2.0 + np.arange(d) / d


def genz(kind: str, x) -> np.ndarray:
    """Genz (1984) family `kind` at the rows of `x` (shape (n, d)), with u_i = 1/2."""
    x = np.atleast_2d(np.asarray(x, dtype=np.float64))
    d = x.shape[1]
    a = genz_a(d)
    u = 0.5
    if kind == "oscillatory":
        return np.cos(2 * np.pi * u + x @ a)
    if kind == "product_peak":
        return np.prod(1.0 / (a ** -2.0 + (x - u) ** 2), axis=1)
    if kind == "corner_peak":
        return (1.0 + x @ a) ** (-(d + 1))
    if kind == "gaussian_peak":
        return np.exp(-np.sum(a**2 * (x - u) ** 2, axis=1))
    if kind == "continuous":
        return np.exp(-np.sum(a * np.abs(x - u), axis=1))
    raise ValueError(f"no Genz family {kind}")


def rosenbrock(x) -> np.ndarray:
    """Sum of 100 (x_{i+1} - x_i^2)^2 + (1 - x_i)^2 at the rows of `x`."""
    x = np.atleast_2d(np.asarray(x, dtype=np.float64))
    return np.sum(100.0 * (x[:, 1:] - x[:, :-1] ** 2) ** 2 + (1.0 - x[:, :-1]) ** 2, axis=1)


# --- Samplers (section 8.4.5) -----------------------------------------------------


def halton(n: int, d: int) -> np.ndarray:
    """The first `n` points of the unscrambled Halton sequence, without the origin."""
    h = qmc.Halton(d=d, scramble=False)
    return h.random(n + 1)[1:]


def uniform(n: int, d: int, seed: int) -> np.ndarray:
    return np.random.default_rng(seed).random((n, d))


def grid_with_holes(m: int, d: int, drop: float, seed: int) -> np.ndarray:
    """An m^d grid on [0, 1]^d with a fraction `drop` of the points removed at random."""
    axes = [np.linspace(0.0, 1.0, m)] * d
    pts = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, d)
    keep = np.random.default_rng(seed).random(len(pts)) >= drop
    return pts[keep]


def geometry(kind: str, n: int, d: int, seed: int) -> np.ndarray:
    """A Qhull `rbox`-like distribution: `cube` in [-0.5, 0.5]^d, `sphere` (the surface
    of the ball of radius 0.5), or `cube_surface`. The distribution matches rbox; the
    generator does not match convx's `benches/sets.txt`."""
    rng = np.random.default_rng(seed)
    if kind == "cube":
        return rng.random((n, d)) - 0.5
    if kind == "sphere":
        g = rng.normal(size=(n, d))
        return 0.5 * g / np.linalg.norm(g, axis=1, keepdims=True)
    if kind == "cube_surface":
        p = rng.random((n, d)) - 0.5
        axis = rng.integers(0, d, size=n)
        sign = rng.integers(0, 2, size=n) * 2 - 1
        p[np.arange(n), axis] = 0.5 * sign
        return p
    raise ValueError(f"no geometry {kind}")


def queries_2d(n: int, seed: int, lo: float = -0.05, hi: float = 1.05) -> np.ndarray:
    """Random queries on a window slightly larger than the unit square."""
    return lo + (hi - lo) * np.random.default_rng(seed).random((n, 2))


# --- Adversarial inputs (section 8.4.6) -------------------------------------------


def _circle_lattice() -> np.ndarray:
    """The 12 integer points with x^2 + y^2 = 25: exactly cocircular."""
    pts = {(0, 5), (0, -5), (5, 0), (-5, 0)}
    for a, b in ((3, 4), (4, 3)):
        for sa in (-1, 1):
            for sb in (-1, 1):
                pts.add((sa * a, sb * b))
    return np.array(sorted(pts), dtype=np.float64)


def _sphere_lattice() -> np.ndarray:
    """The 30 integer points with x^2 + y^2 + z^2 = 9: exactly cospherical."""
    pts = set()
    for s in (-3, 3):
        pts |= {(s, 0, 0), (0, s, 0), (0, 0, s)}
    for perm in ((1, 2, 2), (2, 1, 2), (2, 2, 1)):
        for s0 in (-1, 1):
            for s1 in (-1, 1):
                for s2 in (-1, 1):
                    pts.add((s0 * perm[0], s1 * perm[1], s2 * perm[2]))
    return np.array(sorted(pts), dtype=np.float64)


def adversarial() -> list[dict]:
    """Input point sets for the layer-D tests. The expected behavior is in the design
    (section 8.4.6), not in SciPy: only inputs are stored. `values` is an affine
    function of the points so that `LinearND` can be checked for exact reproduction."""
    out = []

    def add(name, tags, points, queries):
        points = np.asarray(points, dtype=np.float64)
        queries = np.asarray(queries, dtype=np.float64)
        d = points.shape[1]
        coef = np.array([2.0, -3.0, 0.5][:d])
        out.append(
            {
                "name": name,
                "tags": tags,
                "dim": d,
                "points": points,
                "values": 1.0 + points @ coef,
                "queries": queries,
                "queries_values": 1.0 + queries @ coef,
            }
        )

    circ = _circle_lattice()
    add("cocircular_2d", ["cocircular"], circ, [[0.0, 0.0], [1.0, 1.0], [-2.5, 3.0], [5.0, 5.0]])
    sph = _sphere_lattice()
    add("cospherical_3d", ["cospherical"], sph, [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5], [3.0, 3.0, 3.0]])
    add("collinear_2d", ["collinear"], [[0.0, 1.0], [1.0, 3.0], [2.0, 5.0], [3.0, 7.0], [4.0, 9.0]], [[1.5, 4.0]])
    add(
        "coplanar_3d",
        ["coplanar"],
        [[0.0, 0.0, 0.0], [1.0, 0.0, 1.0], [0.0, 1.0, 1.0], [1.0, 1.0, 2.0], [2.0, 1.0, 3.0], [0.5, 0.5, 1.0]],
        [[0.5, 0.25, 0.75]],
    )
    add(
        "sliver_2d",
        ["sliver"],
        [[0.0, 0.0], [1.0, 0.0], [0.5, 1e-9], [0.0, 1.0], [1.0, 1.0]],
        [[0.5, 5e-10], [0.5, 0.5], [0.25, 0.25]],
    )
    add(
        "duplicate_2d",
        ["duplicate"],
        [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0], [0.5, 0.5], [0.5, 0.5]],
        [[0.25, 0.25], [0.75, 0.75]],
    )
    rng = np.random.default_rng(0)
    base = rng.random((12, 2))
    near = np.vstack([base, base[:3] + 1e-12])
    add("near_duplicate_2d", ["near_duplicate"], near, [[0.5, 0.5], [0.3, 0.7]])
    sq = grid_with_holes(5, 2, 0.0, 0)
    add("hull_exterior_2d", ["exterior"], sq, [[-0.5, 0.5], [1.5, 0.5], [0.5, 2.0], [2.0, 2.0], [0.5, 0.5]])
    return out
