"""Generate `goldens/curve.json`: layer-A goldens for N-D parametric curves.

Reference: SciPy 1.18 `make_interp_spline` with a 2-D `y` (one column per output).

SciPy has no function that chooses the parameter of each point for an
interpolating spline. The `parametrized/*` cases compute the parameters here,
by the formulas of docs/design.ja.md section 6.1 (`BSplineCurve`), pass them to
`make_interp_spline`, and store them as `parameters`: the Rust side is given
the points only.

Run: `just gen-goldens`
"""

import numpy as np
from scipy import interpolate as si

from common import golden

GENERATOR = "generate_curve.py"
SEED = 2


def queries_for(t: np.ndarray, n: int = 21, extend: float = 0.05) -> np.ndarray:
    span = t[-1] - t[0]
    grid = np.linspace(t[0] - extend * span, t[-1] + extend * span, n)
    return np.unique(np.concatenate([t, (t[:-1] + t[1:]) / 2.0, grid]))


def chord_parameter(points: np.ndarray) -> np.ndarray:
    """Cumulative chord length scaled to [0, 1]: a non-uniform parameter."""
    step = np.linalg.norm(np.diff(points, axis=0), axis=1)
    t = np.concatenate([[0.0], np.cumsum(step)])
    return t / t[-1]


def parameters(points: np.ndarray, way: str) -> np.ndarray:
    """The parameter of each point, from 0 to 1: `uniform`, `chord_length`, or `centripetal`."""
    step = np.linalg.norm(np.diff(points, axis=0), axis=1)
    if way == "uniform":
        step = np.ones_like(step)
    elif way == "centripetal":
        step = np.sqrt(step)
    else:
        assert way == "chord_length"
    t = np.concatenate([[0.0], np.cumsum(step)])
    return t / t[-1]


def parametrized(name: str, points: np.ndarray, k: int, way: str, bc=None) -> dict:
    """One case whose parameters come from the points."""
    t = parameters(points, way)
    periodic = bc == "periodic"
    q = queries_for(t, extend=0.5 if periodic else 0.05)
    b = si.make_interp_spline(t, points, k=k, bc_type=bc)
    return golden.case(
        f"curve/parametrized/{way}/k{k}/{name}",
        "scipy.interpolate.make_interp_spline",
        {"k": k, "bc_type": bc, "extrapolate": "periodic" if periodic else True, "parametrization": way},
        points=points,
        parameters=t,
        **fields(b, t, q),
    )


def fields(b, t: np.ndarray, q: np.ndarray) -> dict:
    out = {"queries": q, "expected": b(q), "derivative_1": b(q, 1), "k": b.k, "knots": b.t, "control_points": b.c}
    if b.k >= 2:
        out["derivative_2"] = b(q, 2)
    return out


def main() -> None:
    cases = []
    s = np.linspace(0.0, 1.0, 15)
    helix = np.column_stack([np.cos(4 * np.pi * s), np.sin(4 * np.pi * s), s])
    lissajous = np.column_stack([np.sin(3 * np.pi * s), np.cos(2 * np.pi * s)])
    curves = {"helix_3d": (s, helix), "lissajous_2d_chord": (chord_parameter(lissajous), lissajous)}

    for name, (t, pts) in curves.items():
        q = queries_for(t)
        for k in (1, 2, 3, 5):
            b = si.make_interp_spline(t, pts, k=k)
            cases.append(
                golden.case(
                    f"curve/k{k}/{name}",
                    "scipy.interpolate.make_interp_spline",
                    {"k": k, "bc_type": None, "extrapolate": True},
                    t=t,
                    points=pts,
                    **fields(b, t, q),
                )
            )
    t, pts = curves["helix_3d"]
    q = queries_for(t)
    for bc in ("natural", "clamped"):
        b = si.make_interp_spline(t, pts, k=3, bc_type=bc)
        cases.append(
            golden.case(
                f"curve/k3/{bc}/helix_3d",
                "scipy.interpolate.make_interp_spline",
                {"k": 3, "bc_type": bc, "extrapolate": True},
                t=t,
                points=pts,
                **fields(b, t, q),
            )
        )

    # A closed curve: the last point repeats the first, periodic boundary.
    u = np.linspace(0.0, 2.0 * np.pi, 17)
    ellipse = np.column_stack([2.0 * np.cos(u), np.sin(u)])
    ellipse[-1] = ellipse[0]
    q = queries_for(u, extend=0.5)
    b = si.make_interp_spline(u, ellipse, k=3, bc_type="periodic")
    cases.append(
        golden.case(
            "curve/k3/periodic/ellipse_2d",
            "scipy.interpolate.make_interp_spline",
            {"k": 3, "bc_type": "periodic", "extrapolate": "periodic"},
            t=u,
            points=ellipse,
            **fields(b, u, q),
        )
    )

    # Parameters chosen from the points: each way on an open curve in 2D and in
    # 3D. The Lissajous points are unevenly spaced, so the three ways differ.
    for way in ("uniform", "chord_length", "centripetal"):
        cases.append(parametrized("lissajous_2d", lissajous, 3, way))
        cases.append(parametrized("helix_3d", helix, 2, way))
    # Closed curves.
    closed = ellipse[::2].copy()
    closed[-1] = closed[0]
    for k in (2, 3, 5):
        cases.append(parametrized("ellipse_2d", closed, k, "chord_length", "periodic"))
    cases.append(parametrized("ellipse_2d", closed, 3, "centripetal", "periodic"))
    # Four dimensions.
    rng = np.random.default_rng(SEED)
    cases.append(parametrized("random_4d", rng.normal(size=(9, 4)), 3, "chord_length"))
    # As many points as one polynomial of the degree takes.
    cases.append(parametrized("fewest_points_2d", rng.normal(size=(4, 2)), 3, "centripetal"))
    cases.append(parametrized("fewest_points_3d", rng.normal(size=(6, 3)), 5, "chord_length"))
    # Points in clusters: neighbors 1e-3 apart, then 1 apart.
    steps = np.tile([1.0, 1e-3], 6)
    angle = np.linspace(0.0, 2.0, len(steps))
    clustered = np.concatenate([[[0.0, 0.0]], np.cumsum(np.column_stack([steps * np.cos(angle), steps * np.sin(angle)]), axis=0)])
    for way in ("chord_length", "centripetal"):
        cases.append(parametrized("clustered_2d", clustered, 3, way))

    path = golden.write("curve", GENERATOR, SEED, cases)
    print(f"wrote {path.name} ({len(cases)} cases)")


if __name__ == "__main__":
    main()
