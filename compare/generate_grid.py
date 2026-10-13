"""Generate `goldens/grid.json` and `goldens/grid_pchip.json`: layer-A goldens
for `RegularGridInterpolator`.

Reference: SciPy 1.18 `RegularGridInterpolator` with the six methods `linear`,
`nearest`, `slinear`, `cubic`, `quintic`, `pchip`. For `cubic` and `quintic` the
expected values come from `cubic_legacy` and `quintic_legacy` (see `run`). `values`
is stored flat in C order (last axis fastest) with its `shape`, as
docs/design.ja.md section 6.4 says.

`grid.json` has every method on the same grids. `grid_pchip.json` has the cases
chosen for `pchip` (a file has a size budget). Where an axis has fewer than the
four points SciPy's `RegularGridInterpolator` asks of `pchip`, the reference is
`PchipInterpolator` applied by hand, the last axis first (see `fold`).

Run: `just gen-goldens`
"""

import numpy as np
from scipy.interpolate import PchipInterpolator, RegularGridInterpolator, make_interp_spline
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


def fold(axes, values, queries, extrapolate: bool) -> np.ndarray:
    """`PchipInterpolator` along the last axis, then along the one before it,
    down to the first: what `RegularGridInterpolator` does for `pchip`, without
    its lower bound of four points an axis. The axes increase."""
    shape = [len(a) for a in axes]
    out = []
    for q in np.atleast_2d(queries):
        v = np.asarray(values, dtype=np.float64).reshape(shape)
        for d in reversed(range(len(axes))):
            v = PchipInterpolator(axes[d], v, axis=d)(q[d])
        inside = all(a[0] <= c <= a[-1] for a, c in zip(axes, q))
        out.append(float(v) if inside or extrapolate else np.nan)
    return np.array(out)


def run_pchip(name, axes, values, queries, extrapolate: bool, function: str) -> dict:
    """A `pchip` case. With four points on every axis it is SciPy's
    `RegularGridInterpolator`, checked against `fold` where the axes increase;
    with fewer it is `fold`."""
    increasing = all(a[0] < a[1] for a in axes)
    if min(len(a) for a in axes) >= 4:
        case = run(name, axes, values, "pchip", queries, extrapolate, function)
        if increasing:
            by_hand = fold(axes, values, queries, extrapolate)
            assert np.array_equal(np.asarray(case["expected"], dtype=np.float64), by_hand, equal_nan=True), name
        return case
    assert increasing, name
    options = {
        "method": "pchip",
        "bounds_error": False,
        "fill_value": None if extrapolate else "NaN",
        "function": function,
        "reference_method": "PchipInterpolator along each axis, the last axis first",
    }
    return golden.case(
        name,
        "scipy.interpolate.PchipInterpolator",
        options,
        axes=[np.asarray(a) for a in axes],
        shape=[len(a) for a in axes],
        values=np.asarray(values).ravel(),
        queries=queries,
        expected=fold(axes, values, queries, extrapolate),
    )


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

    # An axis of one point (P-4), for the two methods SciPy takes one with.
    # The grid is constant along it: a query has a value at the coordinate of
    # the point, and elsewhere only when the grid is continued. A generator
    # of its own, so that the cases above and the file of `pchip` below stay
    # as they are.
    alone = np.random.default_rng(SEED + 125)
    one_side = nonuniform_axis(alone, 6)
    one_vals = alone.uniform(-1.0, 1.0, (6, 1))
    one_on = np.column_stack([-0.1 + 1.2 * alone.random(16), np.full(16, 0.3)])
    one_off = np.column_stack([-0.1 + 1.2 * alone.random(12), 0.3 + alone.uniform(-0.5, 0.5, 12)])
    one_q = np.vstack([one_on, one_off, [[one_side[2], 0.3], [one_side[0], 0.3], [one_side[5], 0.3]]])
    for method in ("linear", "nearest"):
        one_last = (one_side, np.array([0.3]))
        cases.append(run(f"grid2d/one_point/last/{method}", one_last, one_vals, method, one_q, False, "random"))
        cases.append(run(f"grid2d/one_point/last/{method}/extrapolate", one_last, one_vals, method, one_q, True, "random"))
        one_first = (np.array([0.3]), one_side[::-1].copy())
        cases.append(run(f"grid2d/one_point/first/descending/{method}", one_first, one_vals.T[:, ::-1], method, one_q[:, ::-1], False, "random"))
        one_middle = (nonuniform_axis(alone, 4), np.array([-2.0]), nonuniform_axis(alone, 5))
        one_vals3 = alone.uniform(-1.0, 1.0, (4, 1, 5))
        one_q3 = np.column_stack([alone.random(20), np.where(np.arange(20) % 3 == 0, -1.5, -2.0), alone.random(20)])
        cases.append(run(f"grid3d/one_point/middle/{method}", one_middle, one_vals3, method, one_q3, False, "random"))
        cases.append(run(f"grid3d/one_point/middle/{method}/extrapolate", one_middle, one_vals3, method, one_q3, True, "random"))

    path = golden.write("grid", GENERATOR, SEED, cases)
    print(f"wrote {path.name} ({len(cases)} cases)")

    # Chosen for `pchip` (M2-3). Drawn after everything above, so that
    # `grid.json` stays as it is.
    pchip = []
    # Values in no order: the slopes turn to zero at about half the points, and
    # the order of the axes shows. The same grid with its axes exchanged is
    # another interpolant.
    rough = (nonuniform_axis(rng, 6), nonuniform_axis(rng, 7))
    vals = rng.uniform(-1.0, 1.0, (6, 7))
    qr = np.vstack([rng.random((40, 2)), mesh(rough)[[0, 9, 20, 41]], [[0.0, 0.5], [0.5, 1.0], [1.0, 0.0]]])
    pchip.append(run_pchip("grid2d/random/pchip", rough, vals, qr, False, "random"))
    pchip.append(run_pchip("grid2d/random/exchanged/pchip", rough[::-1], vals.T, qr[:, ::-1], False, "random"))
    rough3 = (nonuniform_axis(rng, 5), nonuniform_axis(rng, 6), nonuniform_axis(rng, 7))
    vals = rng.uniform(-1.0, 1.0, (5, 6, 7))
    qr3 = np.vstack([rng.random((40, 3)), mesh(rough3)[[0, 77, 209]]])
    pchip.append(run_pchip("grid3d/random/pchip", rough3, vals, qr3, False, "random"))
    pchip.append(run_pchip("grid3d/random/pchip/extrapolate", rough3, vals, -0.2 + 1.4 * rng.random((30, 3)), True, "random"))
    # Steps and plateaus: flat where the data is flat, and no value past the
    # data.
    side = np.linspace(0.0, 1.0, 8)
    level = np.array([0.0, 0.0, 0.0, 1.0, 1.0, 3.0, 3.0, 3.0])
    pchip.append(run_pchip("grid2d/steps/pchip", (side, side), np.add.outer(level, 2.0 * level[::-1]), rng.random((60, 2)), False, "steps"))
    # Outside the grid, continued.
    pts = mesh(nonuniform)
    pchip.append(run_pchip("grid2d/nonuniform/f3/pchip/extrapolate", nonuniform, datasets.franke(3, pts[:, 0], pts[:, 1]), far2, True, "franke_f3"))
    pchip.append(run_pchip("grid3d/oscillatory/pchip/extrapolate", axes3, datasets.genz("oscillatory", mesh(axes3)), far3, True, "genz_oscillatory"))
    # Two and three points on an axis: fewer than SciPy's grid takes.
    for shape in ((2, 2), (2, 5), (5, 2), (3, 3), (3, 4), (4, 3), (2, 3, 4), (3, 2, 2), (4, 4, 2)):
        few = tuple(nonuniform_axis(rng, n) for n in shape)
        vals = rng.uniform(-1.0, 1.0, shape)
        queries = np.vstack([-0.15 + 1.3 * rng.random((24, len(shape))), mesh(few)[[0, len(mesh(few)) - 1]]])
        label = "x".join(str(n) for n in shape)
        pchip.append(run_pchip(f"grid{len(shape)}d/few_points/{label}/pchip", few, vals, queries, True, "random"))
    # A descending axis among ascending ones, in 3D.
    mixed = (axes3[0], axes3[1][::-1].copy(), axes3[2])
    pchip.append(run_pchip("grid3d/descending_middle/pchip", mixed, datasets.genz("gaussian_peak", mesh(mixed)), q3[:16], False, "genz_gaussian_peak"))
    # Two values at every point (the layout of the linear cases above).
    pts = mesh(nonuniform)
    both = np.column_stack([datasets.franke(1, pts[:, 0], pts[:, 1]), datasets.franke(4, pts[:, 0], pts[:, 1])])
    case = run_pchip("grid2d/nonuniform/columns/pchip", nonuniform, both[:, 0], q2[:16], False, "franke_f1_f4")
    f = RegularGridInterpolator(nonuniform, both.reshape(shape2 + [2]), method="pchip", bounds_error=False)
    case["values_columns"] = golden.encode(both)
    case["expected_columns"] = golden.encode(f(q2[:16]))
    pchip.append(case)

    path = golden.write("grid_pchip", GENERATOR, SEED, pchip)
    print(f"wrote {path.name} ({len(pchip)} cases)")


if __name__ == "__main__":
    main()
