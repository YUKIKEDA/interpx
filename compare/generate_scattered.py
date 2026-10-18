"""Generate `goldens/scattered.json`: layer-A goldens for scattered data.

Reference: SciPy 1.18 `NearestNDInterpolator`, `LinearNDInterpolator`
(`fill_value` NaN outside the convex hull), and `CloughTocher2DInterpolator`
(`tol = 1e-6`, `maxiter = 400`). Node sets and samplers: docs/design.ja.md
sections 8.4.1 and 8.4.5.

The cases chosen for `NearestNDInterpolator` (`nearest_chosen`) come after the
others: rescaled coordinates, a bound on the distance, points at the same
distance, points given more than once. Those chosen for `LinearNDInterpolator`
(`linear_chosen`, M4-1) come last: rescaled coordinates, points given more than
once, queries on the edges of the hull and at the points, two values at every
point, and larger sets in two and three dimensions.

Run: `just gen-goldens`
"""

import numpy as np
from scipy import interpolate as si
from scipy import spatial
from scipy.spatial import Delaunay

from common import classic_data, datasets, golden

GENERATOR = "generate_scattered.py"
SEED = 4


def node_sets_2d() -> dict[str, np.ndarray]:
    return {
        "franke_100": np.array(classic_data.FRANKE_100),
        "franke_33": np.array(classic_data.FRANKE_33),
        "lawson_25": np.array(classic_data.LAWSON_25),
        "halton_80": datasets.halton(80, 2),
        "uniform_80": datasets.uniform(80, 2, SEED),
        "grid_holes_10": datasets.grid_with_holes(10, 2, 0.4, SEED),
    }


def run(name, reference, options, points, values, queries, build) -> dict:
    common = {"points": points, "values": values, "queries": queries}
    if points.shape[1] <= 4:
        common["n_simplices"] = int(Delaunay(points).simplices.shape[0])
    try:
        return golden.case(name, reference, options, **common, expected=build()(queries))
    except Exception as exc:  # recorded so that the difference is visible in review
        print(f"  {name}: {type(exc).__name__}: {exc}")
        return golden.case(name, reference, options, **common, expected_error=type(exc).__name__)


def three(prefix, points, values, queries, function: str, with_clough_tocher: bool) -> list[dict]:
    out = [
        run(
            f"nearest/{prefix}",
            "scipy.interpolate.NearestNDInterpolator",
            {"function": function},
            points,
            values,
            queries,
            lambda: si.NearestNDInterpolator(points, values),
        ),
        run(
            f"linear/{prefix}",
            "scipy.interpolate.LinearNDInterpolator",
            {"function": function, "fill_value": "NaN", "rescale": False},
            points,
            values,
            queries,
            lambda: si.LinearNDInterpolator(points, values, fill_value=np.nan, rescale=False),
        ),
    ]
    if with_clough_tocher:
        out.append(
            run(
                f"clough_tocher/{prefix}",
                "scipy.interpolate.CloughTocher2DInterpolator",
                {"function": function, "fill_value": "NaN", "tol": 1e-6, "maxiter": 400, "rescale": False},
                points,
                values,
                queries,
                lambda: si.CloughTocher2DInterpolator(points, values, fill_value=np.nan, tol=1e-6, maxiter=400, rescale=False),
            )
        )
    return out


def tied(points, queries, scale=None) -> list[bool]:
    """Whether a query has more than one point at its smallest distance (by a
    scan of all points; the distance is the one the interpolator measures)."""
    pts = np.asarray(points, dtype=np.float64)
    q = np.asarray(queries, dtype=np.float64)
    if scale is not None:
        pts, q = pts / scale, q / scale
    d = np.linalg.norm(pts[None, :, :] - q[:, None, :], axis=2)
    least = d.min(axis=1, keepdims=True)
    return [bool(n > 1) for n in (d <= least * (1.0 + 1e-12) + 1e-300).sum(axis=1)]


def nearest_chosen(rng) -> list[dict]:
    """Cases chosen for `NearestNDInterpolator` (M3-1). Drawn after everything
    else, so that the cases above stay as they are.

    `tied` marks the queries with several points at the same distance: which of
    them SciPy returns follows from how its tree was built, and is not compared.
    """
    reference = "scipy.interpolate.NearestNDInterpolator"
    out = []

    def add(name, options, points, values, queries, build, call=None, **more):
        f = build()
        expected = f(queries, **(call or {}))
        out.append(golden.case(name, reference, options, points=points, values=values, queries=queries, expected=expected, **more))

    # Coordinates of different sizes: a position in [0, 1] and a height in
    # [0, 1000]. As they stand, and rescaled.
    pts = np.column_stack([rng.random(60), 1000.0 * rng.random(60)])
    vals = rng.uniform(-1.0, 1.0, 60)
    q = np.column_stack([-0.1 + 1.2 * rng.random(50), 1000.0 * (-0.1 + 1.2 * rng.random(50))])
    scale = np.ptp(pts, axis=0)
    add("nearest/units/plain", {"rescale": False}, pts, vals, q, lambda: si.NearestNDInterpolator(pts, vals), tied=tied(pts, q))
    add("nearest/units/rescaled", {"rescale": True}, pts, vals, q, lambda: si.NearestNDInterpolator(pts, vals, rescale=True), tied=tied(pts, q, scale))
    # An axis on which all points agree: rescaled, it keeps its scale.
    flat = np.column_stack([10.0 * rng.random(30), np.full(30, 3.0), rng.random(30)])
    vals3 = rng.uniform(-1.0, 1.0, 30)
    q3 = np.column_stack([10.0 * rng.random(30), 3.0 + rng.uniform(-0.5, 0.5, 30), rng.random(30)])
    add("nearest/flat_axis/rescaled", {"rescale": True}, flat, vals3, q3, lambda: si.NearestNDInterpolator(flat, vals3, rescale=True))

    # A bound on the distance: a query has a value when its nearest point is
    # nearer. The queries include some at the bound exactly (3-4-5 triangles
    # from a data point, in numbers that are exact), which have none.
    lattice = np.array([[float(i), float(j)] for i in range(0, 40, 10) for j in range(0, 40, 10)])
    lvals = rng.uniform(-1.0, 1.0, len(lattice))
    exact = np.array([[3.0, 4.0], [13.0, 6.0], [26.0, 23.0], [4.0, 33.0], [5.0, 0.0], [20.0, 25.0]])
    qb = np.vstack([-5.0 + 40.0 * rng.random((40, 2)), exact, lattice[:3]])
    for bound in (5.0, 2.5):
        add(
            f"nearest/bound/{bound}",
            {"rescale": False, "distance_upper_bound": bound},
            lattice, lvals, qb,
            lambda: si.NearestNDInterpolator(lattice, lvals),
            {"distance_upper_bound": bound},
            tied=tied(lattice, qb),
        )
    add(
        "nearest/units/rescaled/bound",
        {"rescale": True, "distance_upper_bound": 0.08},
        pts, vals, q,
        lambda: si.NearestNDInterpolator(pts, vals, rescale=True),
        {"distance_upper_bound": 0.08},
        tied=tied(pts, q, scale),
    )

    # Points at the same distance: the centres of the cells of a lattice, and
    # the middles of their edges.
    centres = np.array([[i + 0.5, j + 0.5] for i in range(0, 30, 10) for j in range(0, 30, 10)]) + 4.5
    edges = np.array([[5.0, 0.0], [10.0, 15.0], [25.0, 30.0]])
    qt = np.vstack([centres, edges, -5.0 + 40.0 * rng.random((20, 2))])
    add("nearest/lattice/ties", {"rescale": False}, lattice, lvals, qt, lambda: si.NearestNDInterpolator(lattice, lvals), tied=tied(lattice, qt))

    # Points given more than once, with other values: the one given first
    # counts. SciPy's own choice among them follows from its tree, so the
    # reference is SciPy on the points without the repeats.
    base = rng.random((25, 3))
    order = rng.permutation(np.concatenate([np.arange(25), np.arange(25), np.arange(10)]))
    dpts = base[order]
    dvals = rng.uniform(-1.0, 1.0, len(dpts))
    _, firsts = np.unique(dpts, axis=0, return_index=True)
    firsts = np.sort(firsts)
    qd = np.vstack([-0.1 + 1.2 * rng.random((40, 3)), base[:5]])
    add(
        "nearest/duplicates",
        {"rescale": False, "reference_points": "the first of each set of equal points"},
        dpts, dvals, qd,
        lambda: si.NearestNDInterpolator(dpts[firsts], dvals[firsts]),
        tied=tied(dpts[firsts], qd),
    )

    # One point; two values at every point. (SciPy takes no points in one
    # dimension; the Rust tests compare that case with a scan.)
    one = np.array([[0.25, -1.0]])
    add("nearest/one_point", {"rescale": False}, one, np.array([4.5]), -2.0 + 4.0 * rng.random((10, 2)), lambda: si.NearestNDInterpolator(one, np.array([4.5])))
    both = np.column_stack([vals, rng.uniform(-1.0, 1.0, 60)])
    add("nearest/units/columns", {"rescale": False}, pts, vals, q, lambda: si.NearestNDInterpolator(pts, vals), tied=tied(pts, q))
    out[-1]["values_columns"] = golden.encode(both)
    out[-1]["expected_columns"] = golden.encode(si.NearestNDInterpolator(pts, both)(q))
    return out


def linear_chosen(rng) -> list[dict]:
    """Cases chosen for `LinearNDInterpolator` (M4-1). Drawn after everything
    else, so that the cases above stay as they are.

    The points are in general position (random), so that the Delaunay
    triangulation is one and the same for Qhull and `convx`: four points on one
    circle would let the two split a cell differently (docs/design.ja.md
    section 6.2).
    """
    reference = "scipy.interpolate.LinearNDInterpolator"
    out = []

    def add(name, options, points, values, queries, build, **more):
        expected = build()(queries)
        out.append(golden.case(name, reference, options, points=points, values=values, queries=queries, expected=expected, **more))

    # Coordinates of different sizes: a position in [0, 1] and a height in
    # [0, 1000]. As they stand, and rescaled (the triangulation differs).
    pts = np.column_stack([rng.random(60), 1000.0 * rng.random(60)])
    vals = rng.uniform(-1.0, 1.0, 60)
    q = np.column_stack([-0.1 + 1.2 * rng.random(50), 1000.0 * (-0.1 + 1.2 * rng.random(50))])
    add("linear/units/plain", {"rescale": False}, pts, vals, q, lambda: si.LinearNDInterpolator(pts, vals, fill_value=np.nan))
    add("linear/units/rescaled", {"rescale": True}, pts, vals, q, lambda: si.LinearNDInterpolator(pts, vals, fill_value=np.nan, rescale=True))

    # Points given more than once, with other values: the one given first
    # counts in interpx. Qhull keeps one of them by its own choice, not always
    # the first, so the reference is SciPy on the points without the repeats.
    base = rng.random((25, 3))
    order = rng.permutation(np.concatenate([np.arange(25), np.arange(25), np.arange(10)]))
    dpts = base[order]
    dvals = rng.uniform(-1.0, 1.0, len(dpts))
    _, firsts = np.unique(dpts, axis=0, return_index=True)
    firsts = np.sort(firsts)
    qd = np.vstack([-0.1 + 1.2 * rng.random((40, 3)), base[:5]])
    add(
        "linear/duplicates",
        {"rescale": False, "reference_points": "the first of each set of equal points"},
        dpts, dvals, qd,
        lambda: si.LinearNDInterpolator(dpts[firsts], dvals[firsts], fill_value=np.nan),
    )

    # Queries on the edges of the hull (midpoints of the edges of the hull, in
    # numbers that are not exact), at the points, and just outside.
    hp = rng.random((40, 2))
    hv = rng.uniform(-1.0, 1.0, 40)
    hull = spatial.ConvexHull(hp)
    mids = np.array([(hp[a] + hp[b]) / 2.0 for a, b in hull.simplices])
    centre = hp.mean(axis=0)
    beyond = centre + 1.001 * (mids - centre)
    qh = np.vstack([mids, hp, beyond])
    add("linear/hull_edges", {"rescale": False}, hp, hv, qh, lambda: si.LinearNDInterpolator(hp, hv, fill_value=np.nan))

    # Two values at every point.
    both = np.column_stack([vals, rng.uniform(-1.0, 1.0, 60)])
    add("linear/units/columns", {"rescale": False}, pts, vals, q, lambda: si.LinearNDInterpolator(pts, vals, fill_value=np.nan))
    out[-1]["values_columns"] = golden.encode(both)
    out[-1]["expected_columns"] = golden.encode(si.LinearNDInterpolator(pts, both, fill_value=np.nan)(q))

    # Larger sets: 2000 points in two dimensions, 500 in three.
    for n, d in ((2000, 2), (500, 3)):
        lp = rng.random((n, d))
        lv = np.sin(3.0 * lp).sum(axis=1)
        lq = -0.05 + 1.1 * rng.random((200, d))
        add(f"linear/random_{n}_d{d}", {"rescale": False}, lp, lv, lq, lambda: si.LinearNDInterpolator(lp, lv, fill_value=np.nan))

    # SciPy 1.18's simplex_tolerance: queries just outside the edges of the
    # hull (moved out from the centre by 1e-14 to 1e-11 of their distance),
    # with the default tolerance and with 100 times it. No random numbers are
    # drawn, so the cases above stay as they are.
    shifts = [1e-14, 1e-13, 1e-12, 1e-11]
    qt = np.vstack([centre + (1.0 + s) * (mids - centre) for s in shifts])
    for tol in (1.0, 100.0):
        name = "linear/hull_edges/tolerance" + ("" if tol == 1.0 else f"_{tol:g}")
        expected = si.LinearNDInterpolator(hp, hv, fill_value=np.nan)(qt, simplex_tolerance=tol)
        out.append(golden.case(name, reference, {"rescale": False, "simplex_tolerance": tol}, points=hp, values=hv, queries=qt, expected=expected))
    return out


def main() -> None:
    cases = []
    for i, (name, pts) in enumerate(node_sets_2d().items()):
        # Random queries on a window larger than the unit square (some fall outside
        # the hull), three data points themselves, and the centroid.
        q = np.vstack([datasets.queries_2d(60, SEED + 10 + i), pts[:3], pts.mean(axis=0, keepdims=True)])
        functions = range(1, 7) if name == "franke_100" else (1,)
        for k in functions:
            vals = datasets.franke(k, pts[:, 0], pts[:, 1])
            cases += three(f"{name}/f{k}", pts, vals, q, f"franke_f{k}", True)

    rng = np.random.default_rng(SEED)
    higher = [
        ("halton_60_d3", datasets.halton(60, 3), "gaussian_peak"),
        ("uniform_60_d3", datasets.uniform(60, 3, SEED + 1), "oscillatory"),
        ("halton_50_d4", datasets.halton(50, 4), "product_peak"),
    ]
    for name, pts, kind in higher:
        d = pts.shape[1]
        q = np.vstack([-0.05 + 1.1 * rng.random((40, d)), pts[:2]])
        cases += three(f"{name}/{kind}", pts, datasets.genz(kind, pts), q, f"genz_{kind}", False)

    cases += nearest_chosen(rng)
    cases += linear_chosen(rng)

    path = golden.write("scattered", GENERATOR, SEED, cases)
    print(f"wrote {path.name} ({len(cases)} cases)")


if __name__ == "__main__":
    main()
