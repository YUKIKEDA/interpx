"""Generate `goldens/grid.json`: layer-A goldens for `RegularGridInterpolator`.

Reference: SciPy 1.18 `RegularGridInterpolator` with the six methods `linear`,
`nearest`, `slinear`, `cubic`, `quintic`, `pchip`. For `cubic` and `quintic` the
expected values come from `cubic_legacy` and `quintic_legacy` (see `run`). `values`
is stored flat in C order (last axis fastest) with its `shape`, as
docs/design.ja.md section 6.4 says.

Run: `just gen-goldens`
"""

import numpy as np
from scipy.interpolate import RegularGridInterpolator, make_interp_spline
from scipy.sparse.linalg import spsolve

from common import datasets, golden

GENERATOR = "generate_grid.py"
SEED = 3
METHODS = ("linear", "nearest", "slinear", "cubic", "quintic", "pchip")
SPLINE_DEGREES = {"cubic": 3, "quintic": 5}


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
        if method in SPLINE_DEGREES:
            # SciPy's default for these two methods is an iterative sparse solver
            # (gcrotmk) whose result is about 1e-5 away from the tensor-product
            # spline on these smooth functions (and of the size of the values on
            # random ones). `expected` is the spline itself, from the legacy
            # method, which interpolates along one axis after the other; the
            # default's result is kept beside it so that the difference stays
            # measurable.
            #
            # A direct sparse solve (`solver=spsolve`) gives the same spline, but
            # less accurately: it solves all axes in one matrix, whose condition
            # number is the product of those of the axes. It is checked here to
            # be close, and not stored (docs/design.ja.md section 6.4 has the
            # comparison with exact arithmetic).
            options["reference_method"] = f"{method}_legacy"
            common["expected_scipy_default"] = RegularGridInterpolator(tuple(axes), grid, **kwargs)(queries)
            direct = RegularGridInterpolator(tuple(axes), grid, **kwargs, solver=spsolve)(queries)
            kwargs["method"] = f"{method}_legacy"
            legacy = RegularGridInterpolator(tuple(axes), grid, **kwargs)(queries)
            inside = ~np.isnan(legacy)
            assert np.array_equal(np.isnan(direct), ~inside), name
            assert np.allclose(direct[inside], legacy[inside], rtol=1e-9, atol=1e-9), name
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

    # Chosen for `nearest` and `linear` (M2-1).
    # Queries halfway between two points of an axis, of one axis and of both:
    # `nearest` takes the lower point.
    small = (np.array([0.0, 1.0, 3.0, 4.0]), np.array([-1.0, 0.0, 2.0]))
    pts = mesh(small)
    vals = 1.0 + pts[:, 0] + 10.0 * pts[:, 1]
    halfway = np.array([[0.5, -0.5], [2.0, 1.0], [3.5, 0.0], [1.0, -0.5], [2.0, 0.25], [0.25, 1.0], [0.5, 1.0000001]])
    cases.append(run("grid2d/halfway/nearest", small, vals, "nearest", halfway, False, "plane"))
    # A descending axis among ascending ones, in 3D: the middle axis.
    mixed = (axes3[0], axes3[1][::-1].copy(), axes3[2])
    vals = datasets.genz("gaussian_peak", mesh(mixed))
    for method in ("linear", "nearest"):
        cases.append(run(f"grid3d/descending_middle/{method}", mixed, vals, method, q3, False, "genz_gaussian_peak"))
    # Halfway on a descending axis: the lower point after the axis is put in
    # ascending order, which is the later of the two as given.
    falling = (small[0][::-1].copy(), small[1])
    pts = mesh(falling)
    vals = 1.0 + pts[:, 0] + 10.0 * pts[:, 1]
    cases.append(run("grid2d/halfway/descending/nearest", falling, vals, "nearest", halfway, False, "plane"))
    # Outside the grid, continued: a non-uniform grid in 2D, and 3D.
    pts = mesh(nonuniform)
    vals = datasets.franke(3, pts[:, 0], pts[:, 1])
    far2 = np.vstack([q2[:12], [[-0.5, 0.5], [1.5, 1.5], [0.5, -0.25], [-0.2, 1.3]]])
    for method in ("linear", "nearest"):
        cases.append(run(f"grid2d/nonuniform/f3/{method}/extrapolate", nonuniform, vals, method, far2, True, "franke_f3"))
    vals = datasets.genz("oscillatory", mesh(axes3))
    far3 = np.vstack([q3[:10], [[-0.3, 0.5, 0.5], [0.5, 1.4, -0.2], [1.2, 1.2, 1.2]]])
    for method in ("linear", "nearest"):
        cases.append(run(f"grid3d/oscillatory/{method}/extrapolate", axes3, vals, method, far3, True, "genz_oscillatory"))
    # Two values at every point. `values` and `expected` hold the first column,
    # as in every other case; `values_columns` and `expected_columns` hold both.
    pts = mesh(nonuniform)
    both = np.column_stack([datasets.franke(1, pts[:, 0], pts[:, 1]), datasets.franke(4, pts[:, 0], pts[:, 1])])
    shape2 = [len(a) for a in nonuniform]
    for method in ("linear", "nearest"):
        case = run(f"grid2d/nonuniform/columns/{method}", nonuniform, both[:, 0], method, q2[:16], False, "franke_f1_f4")
        f = RegularGridInterpolator(nonuniform, both.reshape(shape2 + [2]), method=method, bounds_error=False)
        case["values_columns"] = golden.encode(both)
        case["expected_columns"] = golden.encode(f(q2[:16]))
        cases.append(case)

    # Chosen for `cubic` and `quintic` (M2-2).
    def coefficients(axes, values, k: int) -> np.ndarray:
        """The coefficients of the tensor-product spline: the 1D spline of every
        line along axis 0, then along axis 1, and so on."""
        c = np.asarray(values, dtype=np.float64).reshape([len(a) for a in axes])
        for d, x in enumerate(axes):
            c = np.moveaxis(make_interp_spline(x, np.moveaxis(c, d, 0), k=k).c, 0, d)
        return c.ravel()

    # The coefficients themselves, on a non-uniform grid in 2D and in 3D.
    pts = mesh(nonuniform)
    vals = datasets.franke(2, pts[:, 0], pts[:, 1])
    for method, k in SPLINE_DEGREES.items():
        case = run(f"grid2d/nonuniform/f2/{method}/coefficients", nonuniform, vals, method, q2[:16], True, "franke_f2")
        case["coefficients"] = golden.encode(coefficients(nonuniform, vals, k))
        cases.append(case)
    vals = datasets.genz("gaussian_peak", mesh(axes3))
    for method, k in SPLINE_DEGREES.items():
        case = run(f"grid3d/gaussian_peak/{method}/coefficients", axes3, vals, method, q3[:10], True, "genz_gaussian_peak")
        case["coefficients"] = golden.encode(coefficients(axes3, vals, k))
        cases.append(case)
    # As many points an axis as the degree takes: one polynomial in each
    # coordinate. And one point more on one axis.
    for method, k in SPLINE_DEGREES.items():
        fewest = (nonuniform_axis(rng, k + 1), nonuniform_axis(rng, k + 1))
        pts = mesh(fewest)
        vals = datasets.franke(1, pts[:, 0], pts[:, 1])
        cases.append(run(f"grid2d/fewest_points/{method}", fewest, vals, method, q2[:20], True, "franke_f1"))
        one_more = (nonuniform_axis(rng, k + 2), nonuniform_axis(rng, k + 1))
        pts = mesh(one_more)
        vals = datasets.franke(1, pts[:, 0], pts[:, 1])
        cases.append(run(f"grid2d/one_point_more/{method}", one_more, vals, method, q2[:20], True, "franke_f1"))
    # A descending axis among ascending ones, in 3D.
    mixed = (axes3[0], axes3[1][::-1].copy(), axes3[2])
    vals = datasets.genz("gaussian_peak", mesh(mixed))
    for method in SPLINE_DEGREES:
        cases.append(run(f"grid3d/descending_middle/{method}", mixed, vals, method, q3[:16], False, "genz_gaussian_peak"))
    # Two values at every point (the layout of the linear cases above).
    pts = mesh(nonuniform)
    both = np.column_stack([datasets.franke(1, pts[:, 0], pts[:, 1]), datasets.franke(4, pts[:, 0], pts[:, 1])])
    for method in SPLINE_DEGREES:
        case = run(f"grid2d/nonuniform/columns/{method}", nonuniform, both[:, 0], method, q2[:16], False, "franke_f1_f4")
        f = RegularGridInterpolator(nonuniform, both.reshape(shape2 + [2]), method=f"{method}_legacy", bounds_error=False)
        case["values_columns"] = golden.encode(both)
        case["expected_columns"] = golden.encode(f(q2[:16]))
        cases.append(case)

    path = golden.write("grid", GENERATOR, SEED, cases)
    print(f"wrote {path.name} ({len(cases)} cases)")


if __name__ == "__main__":
    main()
