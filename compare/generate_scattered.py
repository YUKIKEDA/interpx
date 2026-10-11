"""Generate `goldens/scattered.json`: layer-A goldens for scattered data.

Reference: SciPy 1.18 `NearestNDInterpolator`, `LinearNDInterpolator`
(`fill_value` NaN outside the convex hull), and `CloughTocher2DInterpolator`
(`tol = 1e-6`, `maxiter = 400`). Node sets and samplers: docs/design.ja.md
sections 8.4.1 and 8.4.5.

Run: `just gen-goldens`
"""

import numpy as np
from scipy import interpolate as si
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

    path = golden.write("scattered", GENERATOR, SEED, cases)
    print(f"wrote {path.name} ({len(cases)} cases)")


if __name__ == "__main__":
    main()
