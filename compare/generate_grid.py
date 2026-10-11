"""Generate `goldens/grid.json`: layer-A goldens for `RegularGridInterpolator`.

Reference: SciPy 1.18 `RegularGridInterpolator` with the six methods `linear`,
`nearest`, `slinear`, `cubic`, `quintic`, `pchip`. For `cubic` and `quintic` the
expected values come from the direct solver (see `run`). `values` is stored flat in C
order (last axis fastest) with its `shape`, as docs/design.ja.md section 6.4 says.

Run: `just gen-goldens`
"""

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.sparse.linalg import spsolve

from common import datasets, golden

GENERATOR = "generate_grid.py"
SEED = 3
METHODS = ("linear", "nearest", "slinear", "cubic", "quintic", "pchip")
EXACT_SOLVER_METHODS = ("cubic", "quintic")


def nonuniform_axis(rng, n: int) -> np.ndarray:
    return np.concatenate([[0.0], np.sort(rng.uniform(0.05, 0.95, n - 2)), [1.0]])


def mesh(axes) -> np.ndarray:
    """Grid points as rows, in C order."""
    return np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, len(axes))


def run(name, axes, values, method, queries, extrapolate: bool, function: str) -> dict:
    shape = [len(a) for a in axes]
    fill = None if extrapolate else np.nan
    options = {"method": method, "bounds_error": False, "fill_value": None if extrapolate else "NaN", "function": function}
    common = {"axes": [np.asarray(a) for a in axes], "shape": shape, "values": np.asarray(values).ravel(), "queries": queries}
    grid = np.asarray(values).reshape(shape)
    kwargs = {"method": method, "bounds_error": False, "fill_value": fill}
    try:
        if method in EXACT_SOLVER_METHODS:
            # SciPy's default for these two methods is an iterative sparse solver
            # (gcrotmk) whose result is about 1e-5 away from the tensor-product
            # spline. `expected` is the spline itself (direct solver); the default's
            # result is kept beside it so that the difference stays measurable.
            options["solver"] = "scipy.sparse.linalg.spsolve"
            common["expected_scipy_default"] = RegularGridInterpolator(tuple(axes), grid, **kwargs)(queries)
            kwargs["solver"] = spsolve
        f = RegularGridInterpolator(tuple(axes), grid, **kwargs)
        return golden.case(name, "scipy.interpolate.RegularGridInterpolator", options, **common, expected=f(queries))
    except Exception as exc:  # recorded so that the difference is visible in review
        print(f"  {name}: {type(exc).__name__}: {exc}")
        return golden.case(name, "scipy.interpolate.RegularGridInterpolator", options, **common, expected_error=type(exc).__name__)


def main() -> None:
    rng = np.random.default_rng(SEED)
    cases = []

    # 2D, Franke functions.
    uniform = (np.linspace(0.0, 1.0, 11), np.linspace(0.0, 1.0, 9))
    nonuniform = (nonuniform_axis(rng, 10), nonuniform_axis(rng, 9))
    descending = (np.linspace(1.0, 0.0, 11), np.linspace(0.0, 1.0, 9))
    q2 = np.vstack(
        [
            -0.1 + 1.2 * rng.random((40, 2)),
            mesh(uniform)[[0, 13, 50, 98]],
            [[0.5, 0.5], [0.05, 0.95], [1.0, 1.0], [0.0, 0.5]],
        ]
    )
    for k in range(1, 7):
        pts = mesh(uniform)
        vals = datasets.franke(k, pts[:, 0], pts[:, 1])
        for method in METHODS:
            cases.append(run(f"grid2d/uniform/f{k}/{method}", uniform, vals, method, q2, False, f"franke_f{k}"))
    pts = mesh(uniform)
    vals = datasets.franke(1, pts[:, 0], pts[:, 1])
    for method in METHODS:
        cases.append(run(f"grid2d/uniform/f1/{method}/extrapolate", uniform, vals, method, q2, True, "franke_f1"))
    for k in (1, 3):
        pts = mesh(nonuniform)
        vals = datasets.franke(k, pts[:, 0], pts[:, 1])
        for method in METHODS:
            cases.append(run(f"grid2d/nonuniform/f{k}/{method}", nonuniform, vals, method, q2, False, f"franke_f{k}"))
    pts = mesh(descending)
    vals = datasets.franke(1, pts[:, 0], pts[:, 1])
    for method in ("linear", "nearest"):
        cases.append(run(f"grid2d/descending/f1/{method}", descending, vals, method, q2, False, "franke_f1"))

    # 3D, Genz families.
    axes3 = (np.linspace(0.0, 1.0, 6), nonuniform_axis(rng, 7), np.linspace(0.0, 1.0, 6))
    q3 = -0.05 + 1.1 * rng.random((30, 3))
    for kind in ("gaussian_peak", "oscillatory"):
        vals = datasets.genz(kind, mesh(axes3))
        for method in METHODS:
            cases.append(run(f"grid3d/{kind}/{method}", axes3, vals, method, q3, False, f"genz_{kind}"))

    # 4D, one Genz family.
    axes4 = tuple(np.linspace(0.0, 1.0, 6) for _ in range(4))
    q4 = -0.05 + 1.1 * rng.random((20, 4))
    vals = datasets.genz("product_peak", mesh(axes4))
    for method in METHODS:
        cases.append(run(f"grid4d/product_peak/{method}", axes4, vals, method, q4, False, "genz_product_peak"))

    path = golden.write("grid", GENERATOR, SEED, cases)
    print(f"wrote {path.name} ({len(cases)} cases)")


if __name__ == "__main__":
    main()
