"""The SciPy side of the comparison: one case, the phases every runner goes
through (see `runners/harness/src/lib.rs`).

Run by `run.py`: `python scipy_runner.py CASE_DIR`.
"""

import json
import os
import sys
import time
from pathlib import Path

# Before anything large is imported: the process the driver started may be a
# launcher (the `python.exe` of a `uv` environment on Windows starts the
# interpreter as its child), and what the driver has to watch and to pin to a
# core is this process. Its memory read 4 MB for every case until it said so
# (#166).
print(json.dumps({"mark": "pid", "pid": os.getpid()}), flush=True)

import numpy as np  # noqa: E402
from scipy.interpolate import (  # noqa: E402
    Akima1DInterpolator,
    BSpline,
    CloughTocher2DInterpolator,
    CubicSpline,
    LinearNDInterpolator,
    NearestNDInterpolator,
    PchipInterpolator,
    PPoly,
    RBFInterpolator,
    RegularGridInterpolator,
    interp1d,
    make_interp_spline,
)

LIBRARY = "scipy"
# The pause after each phase, for the driver to sample the memory of the
# phase: asked for by the driver for the cases whose memory it reports.
PAUSE = float(os.environ.get("PERF_PAUSE_MS", "0")) / 1000.0
BUDGET = 0.05
# A call shorter than this is repeated inside a sample (the clock has steps).
SAMPLE = 0.0002
FEWEST = 5
MOST = 2000
# How long constructions are made, untimed, before the construction is timed
# (`WARM` in the Rust harness, #172).
WARM = 1.0
# A call that takes this many seconds is timed once, by the call of its phase,
# and is not made again (`LONG` in the Rust harness, #168).
LONG = 2.0


def mark(name: str) -> None:
    print(json.dumps({"mark": name, "t_ns": time.time_ns()}), flush=True)


def timed(run) -> list[float]:
    """The seconds one call takes, sampled several times (as `time` in the Rust harness)."""
    start = time.perf_counter()
    run()
    first = time.perf_counter() - start
    if first >= LONG:
        # A long call is its one sample (as in the Rust harness).
        return [first]
    calls = 1
    while True:
        start = time.perf_counter()
        for _ in range(calls):
            run()
        if time.perf_counter() - start >= SAMPLE or calls >= 1 << 24:
            break
        calls *= 2
    fewest = 3 if first > BUDGET else FEWEST
    seconds = []
    begun = time.perf_counter()
    while len(seconds) < fewest or (time.perf_counter() - begun < BUDGET and len(seconds) < MOST):
        start = time.perf_counter()
        for _ in range(calls):
            run()
        seconds.append((time.perf_counter() - start) / calls)
    return seconds


def oned(case, a):
    x, y, q = a("x"), a("y"), a("q")
    method = case["method"]
    builders = {
        "linear": lambda: interp1d(x, y, kind="linear"),
        "nearest": lambda: interp1d(x, y, kind="nearest"),
        "previous": lambda: interp1d(x, y, kind="previous"),
        "bspline1": lambda: make_interp_spline(x, y, k=1),
        "bspline2": lambda: make_interp_spline(x, y, k=2),
        "bspline3": lambda: make_interp_spline(x, y, k=3),
        "cubic_spline": lambda: CubicSpline(x, y),
        "pchip": lambda: PchipInterpolator(x, y),
        "akima": lambda: Akima1DInterpolator(x, y, method="akima"),
        "makima": lambda: Akima1DInterpolator(x, y, method="makima"),
    }
    return builders[method], lambda built: built(q)


def grid(case, a):
    axes = tuple(a(f"axis{d}") for d in range(case["dim"]))
    values = a("values").reshape(case["shape"])
    q = a("q").reshape(-1, case["dim"])
    method = case["method"]
    return (lambda: RegularGridInterpolator(axes, values, method=method)), lambda built: built(q)


def scattered(case, a):
    dim = case.get("dim", 2)
    points, values, q = a("points").reshape(-1, dim), a("values"), a("q").reshape(-1, dim)
    if case["family"] == "nearest_nd":
        return (lambda: NearestNDInterpolator(points, values)), lambda built: built(q)
    if case["family"] == "linear_nd":
        return (lambda: LinearNDInterpolator(points, values)), lambda built: built(q)
    if case["family"] == "clough_tocher":
        return (lambda: CloughTocher2DInterpolator(points, values)), lambda built: built(q)
    options = {"kernel": case["kernel"], "neighbors": case["neighbors"]}
    if case["kernel"] in ("multiquadric", "inverse_multiquadric", "inverse_quadratic", "gaussian"):
        options["epsilon"] = case["epsilon"]
    return (lambda: RBFInterpolator(points, values, **options)), lambda built: built(q)


def main() -> None:
    directory = Path(sys.argv[1])
    case = json.loads((directory / "case.json").read_text(encoding="utf-8"))

    def a(name: str) -> np.ndarray:
        return np.fromfile(directory / f"{name}.f64", dtype=np.float64)

    family = case["family"]
    if family == "oned":
        build, evaluate = oned(case, a)
    elif family == "ppoly":
        breakpoints, q = a("breakpoints"), a("q")
        c = a("coefficients").reshape(-1, case["degree"] + 1).T.copy()
        build, evaluate = (lambda: PPoly(c, breakpoints)), lambda built: built(q)
    elif family == "bspline":
        knots, control, q = a("knots"), a("control"), a("q")
        degree = case["degree"]
        build, evaluate = (lambda: BSpline(knots, control, degree)), lambda built: built(q)
    elif family == "curve":
        t, points, q = a("t"), a("points").reshape(-1, 2), a("q")
        degree = case["degree"]
        build, evaluate = (lambda: make_interp_spline(t, points, k=degree)), lambda built: built(q)
    elif family == "grid":
        build, evaluate = grid(case, a)
    else:
        build, evaluate = scattered(case, a)
    phases(directory, LIBRARY, build, evaluate)


def phases(directory: Path, library: str, build, evaluate) -> None:
    """The phases of a runner for `build` and `evaluate`; the values go to
    `values.LIBRARY.f64`. Also what `solvers.py` runs SciPy's solvers with."""
    mark("loaded")
    time.sleep(PAUSE)
    start = time.perf_counter()
    built = build()
    first_build = time.perf_counter() - start
    mark("built")
    time.sleep(PAUSE)
    start = time.perf_counter()
    values = evaluate(built)
    first_eval = time.perf_counter() - start
    mark("evaluated")
    time.sleep(PAUSE)
    np.ascontiguousarray(values, dtype=np.float64).tofile(directory / f"values.{library}.f64")
    del values

    if first_build >= LONG:
        build_s = [first_build]
    else:
        del built
        begun = time.perf_counter()
        while time.perf_counter() - begun < WARM:
            build()
        build_s = timed(build)
        built = build()
    eval_s = [first_eval] if first_eval >= LONG else timed(lambda: evaluate(built))
    mark("timed")
    print(json.dumps({"mark": "result", "build_s": build_s, "eval_s": eval_s}), flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:  # noqa: BLE001  (whatever SciPy raises is the result)
        # Said to the driver, so that the report can say why: an array too
        # large to allocate is a result of the comparison, not a fault of it.
        why = f"{type(error).__name__}: {error}"
        print(json.dumps({"mark": "failed", "why": why[:300]}), flush=True)
        sys.exit(1)
