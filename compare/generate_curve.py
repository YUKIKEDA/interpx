"""Generate `goldens/curve.json`: layer-A goldens for N-D parametric curves.

Reference: SciPy 1.18 `make_interp_spline` with a 2-D `y` (one column per output).

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

    path = golden.write("curve", GENERATOR, SEED, cases)
    print(f"wrote {path.name} ({len(cases)} cases)")


if __name__ == "__main__":
    main()
