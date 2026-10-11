"""Generate `goldens/oned.json`: layer-A goldens for the 1D interpolators.

Reference: SciPy 1.18 (`numpy.interp`, `interp1d` for the nearest kinds,
`CubicSpline`, `PchipInterpolator`, `Akima1DInterpolator`, `CubicHermiteSpline`,
`make_interp_spline`). Data: docs/design.ja.md section 8.4.1.

Each case stores the values at the queries and, where the reference has them, the
first and second derivatives, the antiderivative, a few definite integrals, and
the piecewise-polynomial or B-spline representation itself.

Run: `just gen-goldens`
"""

import numpy as np
from scipy import interpolate as si

from common import classic_data, golden

GENERATOR = "generate_oned.py"
SEED = 1


def data_sets() -> dict[str, tuple[np.ndarray, np.ndarray]]:
    rng = np.random.default_rng(SEED)
    ds = {}
    for name, pairs in (
        ("fritsch_carlson", classic_data.FRITSCH_CARLSON),
        ("titanium", classic_data.TITANIUM),
        ("square_with_near_point", classic_data.SQUARE_WITH_NEAR_POINT),
    ):
        a = np.array(pairs)
        ds[name] = (a[:, 0], a[:, 1])
    for n in (11, 21):
        x = np.linspace(-1.0, 1.0, n)
        ds[f"runge_{n}"] = (x, 1.0 / (1.0 + 25.0 * x**2))
    x = np.linspace(0.0, 1.0, 11)
    ds["step_11"] = (x, (x >= 0.5).astype(np.float64))
    x = np.linspace(0.0, 10.0, 41)
    ds["chirp_41"] = (x, np.sin(x**2))
    x = np.sort(rng.uniform(-5.0, 5.0, 20))
    ds["random_20"] = (x, np.sin(x) + 0.1 * rng.normal(size=20))
    x = np.linspace(0.0, 3.0, 12)
    ds["sine_12"] = (x, np.sin(x))
    x = np.linspace(0.0, 2.0 * np.pi, 13)
    y = np.sin(x) + 0.5 * np.cos(2.0 * x)
    y[-1] = y[0]
    ds["periodic_13"] = (x, y)
    return ds


def queries_for(x: np.ndarray, n: int = 31, extend: float = 0.1) -> np.ndarray:
    """Knots, midpoints, and an even grid that reaches past both ends."""
    span = x[-1] - x[0]
    grid = np.linspace(x[0] - extend * span, x[-1] + extend * span, n)
    return np.unique(np.concatenate([x, (x[:-1] + x[1:]) / 2.0, grid]))


def intervals(x: np.ndarray) -> list[tuple[float, float]]:
    """Integration limits inside the domain, one of them reversed."""
    lo, hi = float(x[0]), float(x[-1])
    span = hi - lo
    a, b = lo + 0.3 * span, hi - 0.2 * span
    return [(lo, hi), (a, b), (b, a), (lo, (lo + hi) / 2.0)]


def ppoly_fields(p, x: np.ndarray, q: np.ndarray) -> dict:
    return {
        "queries": q,
        "expected": p(q),
        "derivative_1": p(q, 1),
        "derivative_2": p(q, 2),
        "antiderivative": p.antiderivative()(q),
        "integrals": [[a, b, float(p.integrate(a, b))] for a, b in intervals(x)],
        "breakpoints": p.x,
        "coefficients": p.c,
    }


def bspline_fields(b, x: np.ndarray, q: np.ndarray, periodic: bool = False) -> dict:
    k = b.k
    out = {"queries": q, "expected": b(q), "k": k, "knots": b.t, "control_points": b.c}
    out["derivative_1"] = b(q, 1)
    if k >= 2:
        out["derivative_2"] = b(q, 2)
    if not periodic:
        out["antiderivative"] = b.antiderivative()(q)
        out["integrals"] = [[a, c, float(b.integrate(a, c))] for a, c in intervals(x)]
    return out


def bc_option(bc):
    return bc if isinstance(bc, str) else [list(bc[0]), list(bc[1])]


def main() -> None:
    ds = data_sets()
    cases = []
    non_periodic = [n for n in ds if n != "periodic_13"]

    # numpy.interp: clamps to the end values by default.
    for name in non_periodic:
        x, y = ds[name]
        q = queries_for(x)
        cases.append(golden.case(f"interp/{name}", "numpy.interp", {}, x=x, y=y, queries=q, expected=np.interp(q, x, y)))
    x, y = ds["fritsch_carlson"]
    q = queries_for(x)
    cases.append(
        golden.case(
            "interp/fritsch_carlson/fill",
            "numpy.interp",
            {"left": -1.0, "right": -2.0},
            x=x,
            y=y,
            queries=q,
            expected=np.interp(q, x, y, left=-1.0, right=-2.0),
        )
    )

    # Nearest kinds. Midpoints are among the queries, so `nearest` (rounds half down)
    # and `nearest-up` differ there.
    for name in ("fritsch_carlson", "random_20", "step_11"):
        x, y = ds[name]
        q = queries_for(x)
        for kind in ("nearest", "nearest-up", "previous", "next"):
            f = si.interp1d(x, y, kind=kind, bounds_error=False, fill_value=np.nan)
            cases.append(
                golden.case(
                    f"nearest/{kind}/{name}",
                    "scipy.interpolate.interp1d",
                    {"kind": kind, "bounds_error": False, "fill_value": "NaN"},
                    x=x,
                    y=y,
                    queries=q,
                    expected=f(q),
                )
            )

    # CubicSpline.
    bcs = [
        "natural",
        "clamped",
        ((1, 0.5), (1, -0.5)),
        ((2, 0.25), (2, -0.25)),
        ((1, 0.0), (2, 0.0)),
    ]
    for name in non_periodic:
        x, y = ds[name]
        q = queries_for(x)
        p = si.CubicSpline(x, y, bc_type="not-a-knot")
        cases.append(
            golden.case(
                f"cubic_spline/not-a-knot/{name}",
                "scipy.interpolate.CubicSpline",
                {"bc_type": "not-a-knot", "extrapolate": True},
                x=x,
                y=y,
                **ppoly_fields(p, x, q),
            )
        )
    for name in ("fritsch_carlson", "titanium", "random_20", "sine_12"):
        x, y = ds[name]
        q = queries_for(x)
        for i, bc in enumerate(bcs):
            p = si.CubicSpline(x, y, bc_type=bc)
            label = bc if isinstance(bc, str) else f"derivative_{i}"
            cases.append(
                golden.case(
                    f"cubic_spline/{label}/{name}",
                    "scipy.interpolate.CubicSpline",
                    {"bc_type": bc_option(bc), "extrapolate": True},
                    x=x,
                    y=y,
                    **ppoly_fields(p, x, q),
                )
            )
    for name in ("fritsch_carlson", "sine_12"):
        x, y = ds[name]
        q = queries_for(x)
        p = si.CubicSpline(x, y, bc_type="not-a-knot", extrapolate=False)
        cases.append(
            golden.case(
                f"cubic_spline/not-a-knot/{name}/no_extrapolate",
                "scipy.interpolate.CubicSpline",
                {"bc_type": "not-a-knot", "extrapolate": False},
                x=x,
                y=y,
                **ppoly_fields(p, x, q),
            )
        )
    x, y = ds["periodic_13"]
    q = queries_for(x, extend=0.6)
    p = si.CubicSpline(x, y, bc_type="periodic")
    cases.append(
        golden.case(
            "cubic_spline/periodic/periodic_13",
            "scipy.interpolate.CubicSpline",
            {"bc_type": "periodic", "extrapolate": "periodic"},
            x=x,
            y=y,
            **ppoly_fields(p, x, q),
        )
    )

    # PCHIP.
    for name in non_periodic:
        x, y = ds[name]
        q = queries_for(x)
        p = si.PchipInterpolator(x, y, extrapolate=True)
        cases.append(
            golden.case(
                f"pchip/{name}", "scipy.interpolate.PchipInterpolator", {"extrapolate": True}, x=x, y=y, **ppoly_fields(p, x, q)
            )
        )
    x, y = ds["fritsch_carlson"]
    q = queries_for(x)
    p = si.PchipInterpolator(x, y, extrapolate=False)
    cases.append(
        golden.case(
            "pchip/fritsch_carlson/no_extrapolate",
            "scipy.interpolate.PchipInterpolator",
            {"extrapolate": False},
            x=x,
            y=y,
            **ppoly_fields(p, x, q),
        )
    )

    # Akima and makima.
    for method in ("akima", "makima"):
        for name in ("fritsch_carlson", "titanium", "runge_11", "random_20", "step_11"):
            x, y = ds[name]
            q = queries_for(x)
            p = si.Akima1DInterpolator(x, y, method=method)
            cases.append(
                golden.case(
                    f"akima/{method}/{name}",
                    "scipy.interpolate.Akima1DInterpolator",
                    {"method": method, "extrapolate": None},
                    x=x,
                    y=y,
                    **ppoly_fields(p, x, q),
                )
            )
        x, y = ds["fritsch_carlson"]
        q = queries_for(x)
        p = si.Akima1DInterpolator(x, y, method=method, extrapolate=True)
        cases.append(
            golden.case(
                f"akima/{method}/fritsch_carlson/extrapolate",
                "scipy.interpolate.Akima1DInterpolator",
                {"method": method, "extrapolate": True},
                x=x,
                y=y,
                **ppoly_fields(p, x, q),
            )
        )

    # Cubic Hermite with given slopes.
    rng = np.random.default_rng(SEED + 1)
    hermite = {
        "sine_12": (ds["sine_12"][0], ds["sine_12"][1], np.cos(ds["sine_12"][0])),
        "runge_11": (
            ds["runge_11"][0],
            ds["runge_11"][1],
            -50.0 * ds["runge_11"][0] / (1.0 + 25.0 * ds["runge_11"][0] ** 2) ** 2,
        ),
        "random_20": (ds["random_20"][0], ds["random_20"][1], rng.normal(size=20)),
    }
    for name, (x, y, dydx) in hermite.items():
        q = queries_for(x)
        for extrapolate in (True, False) if name == "sine_12" else (True,):
            p = si.CubicHermiteSpline(x, y, dydx, extrapolate=extrapolate)
            suffix = "" if extrapolate else "/no_extrapolate"
            cases.append(
                golden.case(
                    f"hermite/{name}{suffix}",
                    "scipy.interpolate.CubicHermiteSpline",
                    {"extrapolate": extrapolate},
                    x=x,
                    y=y,
                    dydx=dydx,
                    **ppoly_fields(p, x, q),
                )
            )

    # make_interp_spline, degrees 1 to 5.
    for name in ("fritsch_carlson", "titanium", "random_20", "sine_12", "runge_11"):
        x, y = ds[name]
        q = queries_for(x)
        for k in range(1, 6):
            b = si.make_interp_spline(x, y, k=k)
            cases.append(
                golden.case(
                    f"bspline/k{k}/{name}",
                    "scipy.interpolate.make_interp_spline",
                    {"k": k, "bc_type": None, "extrapolate": True},
                    x=x,
                    y=y,
                    **bspline_fields(b, x, q),
                )
            )
    for name in ("random_20", "sine_12"):
        x, y = ds[name]
        q = queries_for(x)
        for bc in ("not-a-knot", "natural", "clamped"):
            b = si.make_interp_spline(x, y, k=3, bc_type=bc)
            cases.append(
                golden.case(
                    f"bspline/k3/{bc}/{name}",
                    "scipy.interpolate.make_interp_spline",
                    {"k": 3, "bc_type": bc, "extrapolate": True},
                    x=x,
                    y=y,
                    **bspline_fields(b, x, q),
                )
            )
    x, y = ds["fritsch_carlson"]
    q = queries_for(x)
    b = si.make_interp_spline(x, y, k=3)
    cases.append(
        golden.case(
            "bspline/k3/fritsch_carlson/no_extrapolate",
            "scipy.interpolate.make_interp_spline",
            {"k": 3, "bc_type": None, "extrapolate": False},
            x=x,
            y=y,
            queries=q,
            expected=b(q, extrapolate=False),
            derivative_1=b(q, 1, extrapolate=False),
            k=3,
            knots=b.t,
            control_points=b.c,
        )
    )
    x, y = ds["periodic_13"]
    q = queries_for(x, extend=0.6)
    b = si.make_interp_spline(x, y, k=3, bc_type="periodic")
    cases.append(
        golden.case(
            "bspline/k3/periodic/periodic_13",
            "scipy.interpolate.make_interp_spline",
            {"k": 3, "bc_type": "periodic", "extrapolate": "periodic"},
            x=x,
            y=y,
            **bspline_fields(b, x, q, periodic=True),
        )
    )

    # Conversion between the two representations.
    for name in ("sine_12", "random_20"):
        x, y = ds[name]
        b = si.make_interp_spline(x, y, k=3)
        p = si.PPoly.from_spline(b)
        cases.append(
            golden.case(
                f"convert/bspline_to_ppoly/{name}",
                "scipy.interpolate.PPoly.from_spline",
                {"k": 3},
                x=x,
                y=y,
                knots=b.t,
                control_points=b.c,
                breakpoints=p.x,
                coefficients=p.c,
            )
        )

    path = golden.write("oned", GENERATOR, SEED, cases)
    print(f"wrote {path.name} ({len(cases)} cases)")


if __name__ == "__main__":
    main()
