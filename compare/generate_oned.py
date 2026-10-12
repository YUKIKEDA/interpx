"""Generate `goldens/oned.json`, `cubic_spline.json`, `bspline1d.json`, and `calculus.json`: layer-A goldens for the 1D interpolators.

Reference: SciPy 1.18 (`numpy.interp`, `interp1d` for the nearest kinds,
`CubicSpline`, `PchipInterpolator`, `Akima1DInterpolator`, `CubicHermiteSpline`,
`make_interp_spline`). Data: docs/design.ja.md section 8.4.1.

Each case stores the values at the queries and, where the reference has them, the
first and second derivatives, the antiderivative, a few definite integrals, and
the piecewise-polynomial or B-spline representation itself. The `calculus/*`
cases, which go to `calculus.json` (one file has a size budget, and the
`cubic_spline/*` and `bspline/*` cases have files of their own for the same reason), store the derivatives and antiderivatives as objects (their coefficients, or
knots and control points), derivatives of every order, and integrals that leave
the domain or run over several periods.

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


def outside(x: np.ndarray) -> list[tuple[float, float]]:
    """Integration limits that leave the domain, one pair reversed, one all outside."""
    lo, hi = float(x[0]), float(x[-1])
    span = hi - lo
    return [
        (lo - 0.25 * span, hi + 0.1 * span),
        (hi + 0.1 * span, lo + 0.5 * span),
        (lo - 0.2 * span, lo - 0.1 * span),
        (lo + 0.4 * span, hi + 0.3 * span),
    ]


def over_periods(x: np.ndarray) -> list[tuple[float, float]]:
    """Integration limits for a periodic spline: whole periods, parts, and across the end."""
    lo, period = float(x[0]), float(x[-1] - x[0])
    return [
        (lo, lo + 3.0 * period),
        (lo + 0.3 * period, lo + 2.9 * period),
        (lo - 1.7 * period, lo + 0.4 * period),
        (lo + 2.2 * period, lo - 0.3 * period),
        (lo + 0.8 * period, lo + 1.1 * period),
        (lo + 0.2 * period, lo + 0.7 * period),
    ]


def ppoly_calculus(p, q: np.ndarray, limits: list[tuple[float, float]]) -> dict:
    """A piecewise polynomial with its derivatives, antiderivatives, and integrals."""
    k = p.c.shape[0] - 1
    out = {"queries": q, "breakpoints": p.x, "coefficients": p.c}
    # One order past the degree: SciPy returns the zero polynomial.
    out["derivatives"] = [{"nu": nu, "coefficients": p.derivative(nu).c, "values": p(q, nu)} for nu in range(1, k + 2)]
    if p.extrapolate != "periodic":
        # The antiderivative of a periodic polynomial does not extrapolate.
        out["antiderivatives"] = [
            {"nu": nu, "coefficients": p.antiderivative(nu).c, "values": p.antiderivative(nu)(q)} for nu in (1, 2)
        ]
    out["integrals"] = [[a, b, float(p.integrate(a, b))] for a, b in limits]
    return out


def spline_fields(b) -> dict:
    """Knots and control points. SciPy pads `c` to the length of `t`."""
    return {"k": b.k, "knots": b.t, "control_points": b.c[: len(b.t) - b.k - 1]}


def bspline_calculus(b, q: np.ndarray, limits: list[tuple[float, float]]) -> dict:
    """A B-spline with its derivatives, antiderivatives, and integrals."""
    out = {"queries": q, **spline_fields(b)}
    out["derivatives"] = [{"nu": nu, "values": b(q, nu), **spline_fields(b.derivative(nu))} for nu in range(1, b.k + 1)]
    # One order past the degree: SciPy evaluates it to zero, and has no object for it.
    out["derivatives"].append({"nu": b.k + 1, "values": b(q, b.k + 1)})
    if b.extrapolate != "periodic":
        out["antiderivatives"] = [
            {"nu": nu, "values": b.antiderivative(nu)(q), **spline_fields(b.antiderivative(nu))} for nu in (1, 2)
        ]
    out["integrals"] = [[a, c, float(b.integrate(a, c))] for a, c in limits]
    return out


def akima_weight_sums(x: np.ndarray, y: np.ndarray, method: str) -> np.ndarray:
    """The sum of the two weights at every point, as `Akima1DInterpolator` forms it."""
    m = np.empty(len(x) + 3)
    m[2:-2] = np.diff(y) / np.diff(x)
    m[1] = 2.0 * m[2] - m[3]
    m[0] = 2.0 * m[1] - m[2]
    m[-2] = 2.0 * m[-3] - m[-4]
    m[-1] = 2.0 * m[-2] - m[-3]
    dm = np.abs(np.diff(m))
    if method == "makima":
        pm = np.abs(m[1:] + m[:-1])
        return dm[2:] + 0.5 * pm[2:] + dm[:-2] + 0.5 * pm[:-2]
    return dm[2:] + dm[:-2]


def bc_option(bc):
    """A boundary condition as JSON: a name, or one entry per end, each a name or [order, value]."""
    if isinstance(bc, str):
        return bc
    return [end if isinstance(end, str) else [end[0], np.asarray(end[1], dtype=np.float64).tolist()] for end in bc]


def spline_case(name: str, x, y, bc) -> dict:
    """One `CubicSpline` case with everything a piecewise polynomial has."""
    x, y = np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)
    p = si.CubicSpline(x, y, bc_type=bc)
    q = queries_for(x, extend=0.6 if bc == "periodic" else 0.1)
    return golden.case(
        f"cubic_spline/{name}",
        "scipy.interpolate.CubicSpline",
        {"bc_type": bc_option(bc), "extrapolate": "periodic" if bc == "periodic" else True},
        x=x,
        y=y,
        **ppoly_fields(p, x, q),
    )


def spline_columns_case(name: str, x, y, bc) -> dict:
    """One `CubicSpline` case with several values per point: values and derivatives only."""
    x, y = np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)
    p = si.CubicSpline(x, y, bc_type=bc)
    q = queries_for(x, extend=0.6 if bc == "periodic" else 0.1)
    return golden.case(
        f"cubic_spline/{name}",
        "scipy.interpolate.CubicSpline",
        {"bc_type": bc_option(bc), "extrapolate": "periodic" if bc == "periodic" else True},
        x=x,
        y=y,
        queries=q,
        expected=p(q),
        derivative_1=p(q, 1),
        derivative_2=p(q, 2),
    )


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

    # Queries at infinity: `numpy.interp` returns the end values.
    for name in ("fritsch_carlson", "random_20"):
        x, y = ds[name]
        q = np.concatenate([[-np.inf], queries_for(x), [np.inf]])
        cases.append(
            golden.case(f"interp/{name}/infinite", "numpy.interp", {}, x=x, y=y, queries=q, expected=np.interp(q, x, y))
        )

    # The nearest kinds outside the data, where SciPy returns the end values.
    # (`previous` and `next` return NaN on one side each; interpx does not follow
    # that, so they have no case here. docs/design.ja.md section 6.1.)
    for name in ("fritsch_carlson", "random_20"):
        x, y = ds[name]
        q = queries_for(x)
        for kind in ("nearest", "nearest-up"):
            f = si.interp1d(x, y, kind=kind, bounds_error=False, fill_value="extrapolate")
            cases.append(
                golden.case(
                    f"nearest/{kind}/{name}/extrapolate",
                    "scipy.interpolate.interp1d",
                    {"kind": kind, "bounds_error": False, "fill_value": "extrapolate"},
                    x=x,
                    y=y,
                    queries=q,
                    expected=f(q),
                )
            )

    # Values with two columns.
    x, y = ds["random_20"]
    y2 = np.stack([y, np.cos(x)], axis=1)
    q = queries_for(x)
    for kind in ("nearest", "nearest-up", "previous", "next"):
        f = si.interp1d(x, y2, axis=0, kind=kind, bounds_error=False, fill_value=np.nan)
        cases.append(
            golden.case(
                f"nearest/{kind}/random_20/columns",
                "scipy.interpolate.interp1d",
                {"kind": kind, "bounds_error": False, "fill_value": "NaN"},
                x=x,
                y=y2,
                queries=q,
                expected=f(q),
            )
        )

    splines = []
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
        splines.append(
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
            splines.append(
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
        splines.append(
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
    splines.append(
        golden.case(
            "cubic_spline/periodic/periodic_13",
            "scipy.interpolate.CubicSpline",
            {"bc_type": "periodic", "extrapolate": "periodic"},
            x=x,
            y=y,
            **ppoly_fields(p, x, q),
        )
    )

    # Few points, where SciPy leaves the general set-up: two points, three with
    # not-a-knot at both ends (a parabola), three periodic, and four periodic,
    # the smallest cyclic system.
    two = ([0.0, 2.0], [1.0, 4.0])
    for label, bc in (
        ("not-a-knot", "not-a-knot"),
        ("natural", "natural"),
        ("clamped", "clamped"),
        ("first_and_second", ((1, 0.5), (2, -0.25))),
        ("not-a-knot_and_first", ("not-a-knot", (1, 0.5))),
    ):
        splines.append(spline_case(f"{label}/two_points", *two, bc))
    splines.append(spline_case("periodic/two_points", [0.0, 2.0], [1.0, 1.0], "periodic"))
    three = ([0.0, 1.0, 3.0], [1.0, -1.0, 2.0])
    for label, bc in (
        ("not-a-knot", "not-a-knot"),
        ("natural", "natural"),
        ("not-a-knot_and_second", ("not-a-knot", (2, 0.5))),
        ("first_and_not-a-knot", ((1, 0.3), "not-a-knot")),
    ):
        splines.append(spline_case(f"{label}/three_points", *three, bc))
    splines.append(spline_case("periodic/three_points", [0.0, 1.0, 3.0], [1.0, -1.0, 1.0], "periodic"))
    splines.append(spline_case("periodic/four_points", [0.0, 1.0, 1.5, 4.0], [1.0, -1.0, 2.0, 1.0], "periodic"))

    # Not-a-knot at one end and a derivative at the other.
    x, y = ds["random_20"]
    splines.append(spline_case("not-a-knot_and_first/random_20", x, y, ("not-a-knot", (1, 0.5))))
    splines.append(spline_case("second_and_not-a-knot/random_20", x, y, ((2, 0.25), "not-a-knot")))

    # Spacings that alternate between 1 and 1e-3, and a periodic spline on random
    # spacings.
    rng = np.random.default_rng(SEED + 2)
    x = np.cumsum(np.tile([1.0, 1e-3], 8))
    splines.append(spline_case("not-a-knot/uneven_16", x, np.cos(x), "not-a-knot"))
    splines.append(spline_case("natural/uneven_16", x, np.cos(x), "natural"))
    x = np.concatenate([[0.0], np.cumsum(rng.uniform(0.05, 1.0, 14))])
    y = rng.normal(size=15)
    y[-1] = y[0]
    splines.append(spline_case("periodic/random_15", x, y, "periodic"))

    # Two values per point: a derivative at an end has one entry per column.
    x = ds["sine_12"][0]
    y2 = np.stack([np.sin(x), np.cos(2.0 * x)], axis=1)
    splines.append(spline_columns_case("not-a-knot/sine_12/columns", x, y2, "not-a-knot"))
    splines.append(spline_columns_case("natural/sine_12/columns", x, y2, "natural"))
    splines.append(
        spline_columns_case("first_and_second/sine_12/columns", x, y2, ((1, [0.5, -0.5]), (2, [0.0, 0.25])))
    )
    t = np.linspace(0.0, 2.0 * np.pi, 9)
    circle = np.stack([np.cos(t), np.sin(t)], axis=1)
    circle[-1] = circle[0]
    splines.append(spline_columns_case("periodic/circle_9/columns", t, circle, "periodic"))

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

    # PCHIP on data chosen for its rules: few points, flat runs and turns (a zero
    # slope), monotone data, and ends where the three-point estimate is cut to
    # zero or to three times the secant, on the left and on the right.
    rng = np.random.default_rng(SEED + 3)
    shapes = {
        "two_points": ([0.0, 2.0], [1.0, 4.0]),
        "three_points": ([0.0, 1.0, 3.0], [1.0, -1.0, 2.0]),
        "three_points_rising": ([0.0, 1.0, 3.0], [0.0, 1.0, 5.0]),
        "flat_and_turns": ([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0], [0.0, 0.0, 1.0, 3.0, 3.0, 2.0, 2.0, 5.0]),
        "end_zero": ([0.0, 1.0, 2.0, 3.0], [0.0, 1.0, 10.0, 11.0]),
        "end_limited": ([0.0, 1.0, 2.0, 3.0], [0.0, 1.0, -10.0, -9.0]),
        "end_zero_right": ([0.0, 1.0, 2.0, 3.0], [11.0, 10.0, 1.0, 0.0]),
        "end_limited_right": ([0.0, 1.0, 2.0, 3.0], [-9.0, -10.0, 1.0, 0.0]),
    }
    x = np.concatenate([[0.0], np.cumsum(rng.uniform(0.1, 2.0, 11))])
    shapes["monotone_12"] = (x, np.cumsum(np.maximum(rng.normal(0.5, 1.0, 12), 0.0)))
    for name, (x, y) in shapes.items():
        x, y = np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)
        p = si.PchipInterpolator(x, y, extrapolate=True)
        cases.append(
            golden.case(
                f"pchip/{name}",
                "scipy.interpolate.PchipInterpolator",
                {"extrapolate": True},
                x=x,
                y=y,
                **ppoly_fields(p, x, queries_for(x)),
            )
        )
    # Two values per point: each column has its own signs, and so its own slopes.
    x, y = ds["random_20"]
    y2 = np.stack([y, np.cos(x)], axis=1)
    q = queries_for(x)
    p = si.PchipInterpolator(x, y2, extrapolate=True)
    cases.append(
        golden.case(
            "pchip/random_20/columns",
            "scipy.interpolate.PchipInterpolator",
            {"extrapolate": True},
            x=x,
            y=y2,
            queries=q,
            expected=p(q),
            derivative_1=p(q, 1),
            derivative_2=p(q, 2),
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

    # Akima on data chosen for its rules: few points, and data whose secants are
    # equal on each side of a point, where the weights vanish and the slope is the
    # mean of the outer secants (a line, a flat run before a straight rise, steps).
    shapes = {
        "two_points": ([0.0, 2.0], [1.0, 4.0]),
        "three_points": ([0.0, 1.0, 3.0], [1.0, -1.0, 2.0]),
        "four_points": ([0.0, 1.0, 1.5, 4.0], [1.0, -1.0, 2.0, 0.5]),
        "line_6": ([0.0, 1.0, 2.5, 3.0, 4.5, 6.0], [1.0, 3.0, 6.0, 7.0, 10.0, 13.0]),
        "flat_then_rise": ([0.0, 1.0, 2.0, 3.0, 4.0, 5.0], [1.0, 1.0, 1.0, 2.0, 3.0, 4.0]),
        "steps": ([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0], [0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 2.0]),
    }
    for method in ("akima", "makima"):
        for name, (x, y) in shapes.items():
            x, y = np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)
            p = si.Akima1DInterpolator(x, y, method=method)
            cases.append(
                golden.case(
                    f"akima/{method}/{name}",
                    "scipy.interpolate.Akima1DInterpolator",
                    {"method": method, "extrapolate": None},
                    x=x,
                    y=y,
                    **ppoly_fields(p, x, queries_for(x)),
                )
            )
        # Two values per point.
        x, y = ds["random_20"]
        y2 = np.stack([y, np.cos(x)], axis=1)
        q = queries_for(x)
        p = si.Akima1DInterpolator(x, y2, method=method)
        cases.append(
            golden.case(
                f"akima/{method}/random_20/columns",
                "scipy.interpolate.Akima1DInterpolator",
                {"method": method, "extrapolate": None},
                x=x,
                y=y2,
                queries=q,
                expected=p(q),
                derivative_1=p(q, 1),
                derivative_2=p(q, 2),
            )
        )
    # On both sides of the cutoff. Five points carry small secants (1, 2, 5, 3
    # times `scale`), three more carry secants near 5, which set the largest sum
    # of weights. With `scale` 1e-12 the small weights are 1e-12 of the largest:
    # under the cutoff of 1e-9, and the slope there is the mean of the outer
    # secants. With `scale` 1e-7 they are 1e-8 of it: over the cutoff, and the
    # slope is the weighted mean. The values are too small for an absolute
    # tolerance to tell the two apart, so the slopes are stored and compared
    # relatively.
    for label, scale in (("under_the_cutoff", 1e-12), ("over_the_cutoff", 1e-7)):
        x = np.arange(8.0)
        y = np.array([0.0, 1.0, 3.0, 8.0, 11.0, 0.0, 0.0, 0.0]) * scale
        y[5:] = [5.0, 0.0, 5.0]
        q = queries_for(x)
        for method in ("akima", "makima"):
            p = si.Akima1DInterpolator(x, y, method=method)
            cases.append(
                golden.case(
                    f"akima/{method}/{label}",
                    "scipy.interpolate.Akima1DInterpolator",
                    {"method": method, "extrapolate": None},
                    x=x,
                    y=y,
                    queries=q,
                    expected=p(q),
                    slopes=p(x, 1),
                )
            )

    # The same data, with `scale` chosen so that the sum of weights at point 2 is
    # 0.75e-9 of the largest sum, and then 1.5e-9 of it: just under the cutoff and
    # just over it. An implementation that takes the largest weight of one side
    # only, which is about half the largest sum, puts the first of these over
    # its cutoff.
    for label, share in (("just_under_the_cutoff", 0.75e-9), ("just_over_the_cutoff", 1.5e-9)):
        for method in ("akima", "makima"):
            x = np.arange(8.0)
            pattern = np.array([0.0, 1.0, 3.0, 8.0, 11.0, 0.0, 0.0, 0.0])
            tail = np.array([5.0, 0.0, 5.0])
            y = pattern.copy()
            y[5:] = tail
            sums = akima_weight_sums(x, y, method)
            # The weights at point 2 grow with `scale`; the largest sum, among the
            # large secants, does not.
            scale = share * sums.max() / sums[2]
            y = pattern * scale
            y[5:] = tail
            sums = akima_weight_sums(x, y, method)
            assert abs(sums[2] / sums.max() / share - 1.0) < 1e-6, (label, method)
            p = si.Akima1DInterpolator(x, y, method=method)
            q = queries_for(x)
            cases.append(
                golden.case(
                    f"akima/{method}/{label}",
                    "scipy.interpolate.Akima1DInterpolator",
                    {"method": method, "extrapolate": None},
                    x=x,
                    y=y,
                    queries=q,
                    expected=p(q),
                    slopes=p(x, 1),
                )
            )

    # Where SciPy and interpx differ: SciPy takes the cutoff for vanishing weights
    # from the largest weight over all columns. Beside a column of values near 1,
    # a column of values near 1e-10 has all its weights under the cutoff, and its
    # slopes are not those it has alone. `expected` is the first column beside the
    # other, `expected_alone` the first column by itself; interpx gives the second
    # for both (docs/design.ja.md section 6.1).
    x = np.arange(8.0)
    small = np.array([0.0, 0.0, 0.0, 1e-10, 2e-10, 2e-10, 2e-10, 5e-10])
    big = np.array([0.0, 5.0, -3.0, 8.0, 1.0, 9.0, -4.0, 2.0])
    q = np.linspace(0.0, 7.0, 29)
    for method in ("akima", "makima"):
        paired = si.Akima1DInterpolator(x, np.stack([small, big], axis=1), method=method)
        alone = si.Akima1DInterpolator(x, small, method=method)
        cases.append(
            golden.case(
                f"akima/{method}/coupled_columns",
                "scipy.interpolate.Akima1DInterpolator",
                {"method": method, "extrapolate": None},
                x=x,
                y=np.stack([small, big], axis=1),
                queries=q,
                expected=paired(q),
                expected_alone=alone(q),
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

    # Chosen for the type: two and three points, widths that alternate between 1
    # and 1e-2, slopes that are all zero, a line with its slope, slopes that
    # disagree with the data (the interpolant overshoots), and periodic
    # extrapolation of data whose ends agree and of data whose ends do not.
    def hermite_case(name: str, x, y, dydx, extrapolate=True) -> dict:
        x, y, dydx = (np.asarray(a, dtype=np.float64) for a in (x, y, dydx))
        p = si.CubicHermiteSpline(x, y, dydx, extrapolate=extrapolate)
        q = queries_for(x, extend=0.6 if extrapolate == "periodic" else 0.1)
        return golden.case(
            f"hermite/{name}",
            "scipy.interpolate.CubicHermiteSpline",
            {"extrapolate": extrapolate},
            x=x,
            y=y,
            dydx=dydx,
            **ppoly_fields(p, x, q),
        )

    cases.append(hermite_case("two_points", [0.0, 2.0], [1.0, 3.0], [-1.0, 0.5]))
    cases.append(hermite_case("three_points", [0.0, 1.0, 3.0], [1.0, -1.0, 2.0], [0.5, -2.0, 1.0]))
    x = np.cumsum(np.concatenate([[0.0], np.tile([1.0, 1e-2], 8)]))
    cases.append(hermite_case("uneven", x, np.sin(x), np.cos(x)))
    x, y = ds["random_20"]
    cases.append(hermite_case("zero_slopes", x, y, np.zeros_like(x)))
    cases.append(hermite_case("line", x, 2.0 * x - 1.0, np.full_like(x, 2.0)))
    cases.append(
        hermite_case(
            "against_the_data",
            np.arange(6.0),
            [0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            [5.0, -5.0, 5.0, -5.0, 5.0, -5.0],
        )
    )
    x = np.linspace(0.0, 2.0 * np.pi, 9)
    y = np.sin(x)
    y[-1] = y[0]
    dydx = np.cos(x)
    dydx[-1] = dydx[0]
    cases.append(hermite_case("periodic/sine_9", x, y, dydx, extrapolate="periodic"))
    x, y = ds["sine_12"]
    cases.append(hermite_case("periodic/open_ends", x, y, np.cos(x), extrapolate="periodic"))

    # Two columns, each with its own slopes: values and derivatives only.
    x = ds["random_20"][0]
    p = si.CubicHermiteSpline(
        x, np.stack([np.sin(x), np.exp(-x)], axis=1), np.stack([np.cos(x), -np.exp(-x)], axis=1)
    )
    q = queries_for(x)
    cases.append(
        golden.case(
            "hermite/random_20/columns",
            "scipy.interpolate.CubicHermiteSpline",
            {"extrapolate": True},
            x=x,
            y=np.stack([np.sin(x), np.exp(-x)], axis=1),
            dydx=np.stack([np.cos(x), -np.exp(-x)], axis=1),
            queries=q,
            expected=p(q),
            derivative_1=p(q, 1),
            derivative_2=p(q, 2),
        )
    )

    # A file of their own: `bspline1d.json`.
    bsplines = []
    # make_interp_spline, degrees 1 to 5.
    for name in ("fritsch_carlson", "titanium", "random_20", "sine_12", "runge_11"):
        x, y = ds[name]
        q = queries_for(x)
        for k in range(1, 6):
            b = si.make_interp_spline(x, y, k=k)
            bsplines.append(
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
            bsplines.append(
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
    bsplines.append(
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
    bsplines.append(
        golden.case(
            "bspline/k3/periodic/periodic_13",
            "scipy.interpolate.make_interp_spline",
            {"k": 3, "bc_type": "periodic", "extrapolate": "periodic"},
            x=x,
            y=y,
            **bspline_fields(b, x, q, periodic=True),
        )
    )

    # Chosen for `BSpline1D`. Fewer queries than above: the file has a size budget.
    def interp_case(name: str, x, y, k: int, bc=None) -> dict:
        x, y = np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)
        periodic = bc == "periodic"
        # A derivative at an end is (order, value), as for `CubicSpline`;
        # `make_interp_spline` takes a list of such pairs per end.
        lists = bc if bc is None or isinstance(bc, str) else ([bc[0]], [bc[1]])
        b = si.make_interp_spline(x, y, k=k, bc_type=lists)
        q = queries_for(x, n=15, extend=0.6 if periodic else 0.1)
        return golden.case(
            f"bspline/{name}",
            "scipy.interpolate.make_interp_spline",
            {"k": k, "bc_type": None if bc is None else bc_option(bc), "extrapolate": "periodic" if periodic else True},
            x=x,
            y=y,
            **bspline_fields(b, x, q, periodic=periodic),
        )

    def interp_columns_case(name: str, x, y, k: int, bc=None) -> dict:
        """Several values per point: values and derivatives only."""
        x, y = np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)
        periodic = bc == "periodic"
        b = si.make_interp_spline(x, y, k=k, bc_type=bc)
        q = queries_for(x, n=15, extend=0.6 if periodic else 0.1)
        return golden.case(
            f"bspline/{name}",
            "scipy.interpolate.make_interp_spline",
            {"k": k, "bc_type": None if bc is None else bc_option(bc), "extrapolate": "periodic" if periodic else True},
            x=x,
            y=y,
            queries=q,
            expected=b(q),
            derivative_1=b(q, 1),
            derivative_2=b(q, 2),
        )

    def one_sided_case(name: str, x, y, bc) -> dict:
        """Degree 3 with not-a-knot at one end and a derivative at the other.
        `make_interp_spline` does not take it; `CubicSpline` gives the curve."""
        x, y = np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)
        p = si.CubicSpline(x, y, bc_type=bc)
        q = queries_for(x, n=15)
        return golden.case(
            f"bspline/{name}",
            "scipy.interpolate.CubicSpline",
            {"k": 3, "bc_type": bc_option(bc), "extrapolate": True},
            x=x,
            y=y,
            queries=q,
            expected=p(q),
            derivative_1=p(q, 1),
            derivative_2=p(q, 2),
        )

    # Degrees over 5.
    for name in ("sine_12", "random_20"):
        for k in (6, 7):
            bsplines.append(interp_case(f"k{k}/{name}", *ds[name], k))
    # As many points as one polynomial of the degree takes.
    rng = np.random.default_rng(SEED + 4)
    for k in (2, 3, 4, 5):
        x = np.cumsum(rng.uniform(0.3, 1.2, k + 1))
        bsplines.append(interp_case(f"k{k}/fewest_points", x, rng.normal(size=k + 1), k))
    # Widths that alternate between 1 and 1e-2.
    x = np.cumsum(np.concatenate([[0.0], np.tile([1.0, 1e-2], 8)]))
    for k in (2, 3, 5):
        bsplines.append(interp_case(f"k{k}/uneven_17", x, np.sin(x), k))

    # Periodic, the other degrees (degree 1 is left out: SciPy does not make it
    # periodic, see the design).
    x, y = ds["periodic_13"]
    for k in (2, 4, 5):
        bsplines.append(interp_case(f"k{k}/periodic/periodic_13", x, y, k, "periodic"))
    x = np.concatenate([[0.0], np.cumsum(rng.uniform(0.05, 1.0, 14))])
    y = rng.normal(size=15)
    y[-1] = y[0]
    for k in (2, 3, 5):
        bsplines.append(interp_case(f"k{k}/periodic/random_15", x, y, k, "periodic"))
    # Few points: as many unknowns as the degree or fewer, so that the band of
    # the cyclic matrix wraps onto itself; one more; and a constant.
    four = ([0.0, 1.0, 1.5, 4.0], [1.0, -1.0, 2.0, 1.0])
    for k in (2, 3, 4, 5):
        bsplines.append(interp_case(f"k{k}/periodic/four_points", *four, k, "periodic"))
    bsplines.append(interp_case("k2/periodic/three_points", [0.0, 1.0, 3.0], [1.0, -1.0, 1.0], 2, "periodic"))
    for k in (2, 3):
        bsplines.append(interp_case(f"k{k}/periodic/two_points", [0.0, 2.0], [1.5, 1.5], k, "periodic"))

    # Degree 3 with derivatives at the ends that are not zero, and with the
    # fewest points.
    x, y = ds["random_20"]
    bsplines.append(interp_case("k3/first_and_second/random_20", x, y, 3, ((1, 0.5), (2, -0.25))))
    x, y = ds["sine_12"]
    bsplines.append(interp_case("k3/second_and_first/sine_12", x, y, 3, ((2, 1.0), (1, -2.0))))
    bsplines.append(interp_case("k3/first_and_first/two_points", [0.0, 2.0], [1.0, 3.0], 3, ((1, -1.0), (1, 0.5))))
    bsplines.append(
        interp_case("k3/second_and_second/three_points", [0.0, 1.0, 3.0], [1.0, -1.0, 2.0], 3, ((2, 0.5), (2, -1.0)))
    )
    # Not-a-knot at one end.
    x, y = ds["random_20"]
    bsplines.append(one_sided_case("k3/not-a-knot_and_first/random_20", x, y, ("not-a-knot", (1, 0.5))))
    bsplines.append(one_sided_case("k3/second_and_not-a-knot/random_20", x, y, ((2, 0.25), "not-a-knot")))
    bsplines.append(
        one_sided_case("k3/first_and_not-a-knot/three_points", [0.0, 1.0, 3.0], [1.0, -1.0, 2.0], ((1, 0.5), "not-a-knot"))
    )

    # Two columns.
    x = ds["sine_12"][0]
    y2 = np.stack([np.sin(x), np.cos(2.0 * x)], axis=1)
    bsplines.append(interp_columns_case("k3/sine_12/columns", x, y2, 3))
    bsplines.append(interp_columns_case("k4/sine_12/columns", x, y2, 4))
    t = np.linspace(0.0, 2.0 * np.pi, 9)
    circle = np.stack([np.cos(t), np.sin(t)], axis=1)
    circle[-1] = circle[0]
    bsplines.append(interp_columns_case("k5/periodic/circle_9/columns", t, circle, 5, "periodic"))

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

    calculus = []
    # Derivatives and integrals as objects: the coefficients SciPy gives a
    # derivative and an antiderivative, derivatives of every order up to one past
    # the degree, and definite integrals that leave the domain.
    for name in ("sine_12", "random_20"):
        x, y = ds[name]
        q = queries_for(x)
        p = si.CubicSpline(x, y)
        calculus.append(
            golden.case(
                f"calculus/ppoly/cubic_spline/{name}",
                "scipy.interpolate.PPoly",
                {"extrapolate": True},
                **ppoly_calculus(p, q, outside(x)),
            )
        )
        for k in (2, 3, 5):
            b = si.make_interp_spline(x, y, k=k)
            calculus.append(
                golden.case(
                    f"calculus/bspline/k{k}/{name}",
                    "scipy.interpolate.BSpline",
                    {"extrapolate": True},
                    **bspline_calculus(b, q, outside(x)),
                )
            )
    x, y = ds["fritsch_carlson"]
    p = si.PchipInterpolator(x, y, extrapolate=True)
    calculus.append(
        golden.case(
            "calculus/ppoly/pchip/fritsch_carlson",
            "scipy.interpolate.PPoly",
            {"extrapolate": True},
            **ppoly_calculus(p, queries_for(x), outside(x)),
        )
    )
    x, y = ds["periodic_13"]
    q = queries_for(x, extend=0.6)
    p = si.CubicSpline(x, y, bc_type="periodic")
    calculus.append(
        golden.case(
            "calculus/ppoly/cubic_spline/periodic_13",
            "scipy.interpolate.PPoly",
            {"extrapolate": "periodic"},
            **ppoly_calculus(p, q, over_periods(x)),
        )
    )
    b = si.make_interp_spline(x, y, k=3, bc_type="periodic")
    calculus.append(
        golden.case(
            "calculus/bspline/k3/periodic_13",
            "scipy.interpolate.BSpline",
            {"extrapolate": "periodic"},
            **bspline_calculus(b, q, over_periods(x)),
        )
    )

    path = golden.write("oned", GENERATOR, SEED, cases)
    print(f"wrote {path.name} ({len(cases)} cases)")
    path = golden.write("bspline1d", GENERATOR, SEED, bsplines)
    print(f"wrote {path.name} ({len(bsplines)} cases)")
    path = golden.write("cubic_spline", GENERATOR, SEED, splines)
    print(f"wrote {path.name} ({len(splines)} cases)")
    path = golden.write("calculus", GENERATOR, SEED, calculus)
    print(f"wrote {path.name} ({len(calculus)} cases)")


if __name__ == "__main__":
    main()
