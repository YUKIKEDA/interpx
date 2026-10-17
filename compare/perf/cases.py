"""The cases of the comparison, and the data of each.

A case is one interpolation problem at one size: its parameters (`case.json`)
and its arrays, written once into a directory that every runner reads, so that
every library is given the same numbers. The sizes are those of SciPy's own
benchmark suite (`common/bench_cases.py`, docs/design.ja.md section 8.4.3),
and a few larger ones.

An array is a file of 64-bit floats in the byte order of the machine
(`NAME.f64`), rows in C order.

The cases on real-world data (`real.py`, design section 8.4.7) are in the
table when the data has been fetched. They have the families and the arrays
of the other cases, so every runner takes them as it takes those, and one
array more: `truth`, the values at the queries, which are points of the data
that were held out.
"""

import functools
import json
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline, make_interp_spline

# `common` is next to this directory (`compare/common`).
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common.bench_cases import RUST_ONLY, SCIPY  # noqa: E402

import real  # noqa: E402  (next to this file)

SEED = 20261013

# The cases whose memory over time is reported and drawn: the large ones. A
# runner pauses after each phase of these, so that the phase can be told
# apart in the samples; for the rest memory is a few megabytes whatever the
# library, and the pauses were most of the time of a run.
MEMORY = (
    "rbf/thin_plate_spline/global/n3000/scattered",
    "rbf/thin_plate_spline/k50/n100000/image",
    "nearest_nd/n100000/image500",
    "linear_nd/n100000/image500",
    "grid/cubic/2d/each400/eval1000/ascending",
    "grid/linear/2d/each800/eval10000/ascending",
    "grid/linear/4d/max200/eval10000/ascending",
    "oned/cubic_spline/n10000/sorted",
)
PAUSE_MS = 60

ONED_METHODS = (
    "linear",
    "nearest",
    "previous",
    "bspline1",
    "bspline2",
    "bspline3",
    "cubic_spline",
    "pchip",
    "akima",
    "makima",
)


def rng_for(name: str) -> np.random.Generator:
    """A generator that depends on the case alone, not on the cases before it."""
    return np.random.default_rng([SEED, *name.encode()])


def oned(n: int, method: str, order: str, even: bool = True) -> dict:
    name = f"oned/{method}/n{n}/{order}" + ("" if even else "-uneven")
    # Evenly spaced, as in SciPy's suite; or spacings of 0.5 to 1.5, for the
    # libraries that find an interval faster on even spacing than otherwise.
    x = np.arange(n, dtype=float) if even else np.cumsum(rng_for(f"oned/uneven/{n}").uniform(0.5, 1.5, n))
    q = np.linspace(x[0], x[-1], 4 * n)
    if order == "shuffled":
        rng_for(name).shuffle(q)
    return {
        "name": name,
        "family": "oned",
        "method": method,
        "n": n,
        "order": order,
        # A decay over the whole range. SciPy's suite has `exp(-x / 3)` at
        # every size, which is below the smallest normal number from x = 2100
        # on: the times would then be those of arithmetic on such numbers.
        "arrays": {"x": x, "y": np.exp(-3.0 * (x - x[0]) / (x[-1] - x[0] + 1.0)), "q": q},
    }


def piecewise(kind: str, n_eval: int) -> dict:
    name = f"{kind}/intervals55/eval{n_eval}"
    rng = rng_for(f"{kind}/intervals55")
    x = np.cumsum(rng.uniform(0.5, 1.5, 56))
    y = rng.uniform(-1.0, 1.0, 56)
    q = rng_for(name).uniform(x[0], x[-1], n_eval)
    case = {"name": name, "family": kind, "degree": 3, "n": 55, "n_eval": n_eval}
    if kind == "ppoly":
        # One interval after another, the highest power first (interpx's
        # layout; SciPy's `c[k, m]` is the transpose).
        c = CubicSpline(x, y).c
        case["arrays"] = {"breakpoints": x, "coefficients": c.T.copy(), "q": q}
    else:
        spline = make_interp_spline(x, y, k=3)
        case["arrays"] = {"knots": spline.t, "control": spline.c, "q": q}
    return case


def curve(n: int) -> dict:
    name = f"curve/degree3/n{n}"
    angle = np.linspace(0.0, 6.0 * np.pi, n)
    points = np.stack([angle * np.cos(angle), angle * np.sin(angle)], axis=1)
    steps = np.linalg.norm(np.diff(points, axis=0), axis=1)
    t = np.concatenate([[0.0], np.cumsum(steps)])
    return {
        "name": name,
        "family": "curve",
        "degree": 3,
        "n": n,
        "arrays": {"t": t, "points": points, "q": np.linspace(t[0], t[-1], 4 * n)},
    }


def grid(shape: list[int], method: str, n_eval: int, flipped: bool, label: str) -> dict:
    dim = len(shape)
    name = f"grid/{method}/{dim}d/{label}/eval{n_eval}/{'flipped' if flipped else 'ascending'}"
    rng = rng_for(f"grid/{dim}d/{label}")
    axes = [np.cumsum(rng.uniform(0.5, 1.5, size)) for size in shape]
    values = rng.uniform(-1.0, 1.0, shape)
    q = np.stack([rng_for(name).uniform(a[0], a[-1], n_eval) for a in axes], axis=1)
    if flipped:
        # Every axis given from its high end, and the values with them.
        axes = [a[::-1].copy() for a in axes]
        values = values[(slice(None, None, -1),) * dim].copy()
    arrays = {f"axis{d}": a for d, a in enumerate(axes)}
    arrays |= {"values": values, "q": q}
    return {
        "name": name,
        "family": "grid",
        "method": method,
        "dim": dim,
        "shape": shape,
        "n": int(np.prod(shape)),
        "n_eval": n_eval,
        "arrays": arrays,
    }


def image(per_axis: int, low: float, high: float) -> np.ndarray:
    """The pixels of a square image, row after row."""
    side = np.linspace(low, high, per_axis)
    return np.stack(np.meshgrid(side, side, indexing="ij"), axis=-1).reshape(-1, 2)


def nearest_nd(n: int, per_axis: int) -> dict:
    name = f"nearest_nd/n{n}/image{per_axis}"
    rng = rng_for(f"nearest_nd/n{n}")
    points = rng.uniform(0.0, 1.0, (n, 2))
    return {
        "name": name,
        "family": "nearest_nd",
        "n": n,
        "n_eval": per_axis * per_axis,
        "arrays": {"points": points, "values": rng.uniform(-1.0, 1.0, n), "q": image(per_axis, 0.0, 1.0)},
    }


def linear_nd(n: int, per_axis: int) -> dict:
    """The points of `nearest_nd`, evaluated by linear interpolation on their
    Delaunay triangulation (SciPy's suite times it through `griddata`)."""
    name = f"linear_nd/n{n}/image{per_axis}"
    rng = rng_for(f"linear_nd/n{n}")
    points = rng.uniform(0.0, 1.0, (n, 2))
    return {
        "name": name,
        "family": "linear_nd",
        "dim": 2,
        "n": n,
        "n_eval": per_axis * per_axis,
        "arrays": {"points": points, "values": rng.uniform(-1.0, 1.0, n), "q": image(per_axis, 0.0, 1.0)},
    }


def linear_nd_3d(n: int) -> dict:
    """`n` points of the unit cube, evaluated at 10 000 points of the cube in
    no order."""
    name = f"linear_nd/3d/n{n}/scattered"
    rng = rng_for(f"linear_nd/3d/n{n}")
    points = rng.uniform(0.0, 1.0, (n, 3))
    values = np.sin(3.0 * points).sum(axis=1)
    return {
        "name": name,
        "family": "linear_nd",
        "dim": 3,
        "n": n,
        "n_eval": 10000,
        "arrays": {"points": points, "values": values, "q": rng.uniform(0.0, 1.0, (10000, 3))},
    }


def rbf(n: int, kernel: str, neighbors: int | None, queries: str) -> dict:
    mode = "global" if neighbors is None else f"k{neighbors}"
    name = f"rbf/{kernel}/{mode}/n{n}/{queries}"
    spec = SCIPY["rbf_interpolator"]
    low, high = spec["domain"]
    rng = rng_for(f"rbf/n{n}")
    points = rng.uniform(low, high, (n, 2))
    values = np.sin(3.0 * points[:, 0]) * np.cos(2.0 * points[:, 1])
    if queries == "scattered":
        q = rng_for(f"rbf/q/{n}").uniform(low, high, (1000, 2))
    else:
        q = image(100, low, high)
    return {
        "name": name,
        "family": "rbf",
        "kernel": kernel,
        "epsilon": spec["epsilon"],
        "neighbors": neighbors,
        "n": n,
        "n_eval": len(q),
        "arrays": {"points": points, "values": values, "q": q},
    }


# The grid methods run on the elevation tile: all of them.
REAL_GRID_METHODS = ("nearest", "linear", "cubic", "quintic", "pchip")
# The kernels run on the stations: those without a shape parameter. A
# parameter would be ours to choose, and the result would follow the choice.
REAL_KERNELS = ("linear", "thin_plate_spline", "cubic", "quintic")
REAL_NEIGHBORS = (None, 50)


@functools.cache
def elevation_arrays() -> dict:
    return real.elevation_split()


@functools.cache
def weather_arrays() -> dict:
    return real.weather_split()


def real_cases() -> list[dict]:
    """The cases on real-world data. Their arrays are read when a case is
    written, not when the table is made: the tile is a quarter of a gigabyte
    of queries."""
    tile, day = real.ELEVATION["name"], real.WEATHER["name"]
    # The sizes, without reading the data: they are those of the sources.
    rows, columns = len(real.kept(7200)), len(real.kept(4800))
    cases = []
    for method in REAL_GRID_METHODS:
        cases.append({
            "name": f"grid/{method}/2d/{tile}/heldout/ascending",
            "family": "grid",
            "method": method,
            "dim": 2,
            "shape": [rows, columns],
            "n": rows * columns,
            "n_eval": 7200 * 4800 - rows * columns,
            "real": tile,
            "arrays": elevation_arrays,
        })
    cases.append({"name": f"nearest_nd/{day}/heldout", "family": "nearest_nd", "real": day, "arrays": weather_arrays})
    cases.append({"name": f"linear_nd/{day}/heldout", "family": "linear_nd", "dim": 2, "real": day, "arrays": weather_arrays})
    for kernel in REAL_KERNELS:
        for neighbors in REAL_NEIGHBORS:
            mode = "global" if neighbors is None else f"k{neighbors}"
            cases.append({
                "name": f"rbf/{kernel}/{mode}/{day}/heldout",
                "family": "rbf",
                "kernel": kernel,
                "epsilon": 1.0,
                "neighbors": neighbors,
                "real": day,
                "arrays": weather_arrays,
            })
    return cases


def all_cases() -> list[dict]:
    cases = []
    for n in SCIPY["interp1d"]["n_samples"]:
        for method in ONED_METHODS:
            for order in ("sorted", "shuffled"):
                cases.append(oned(n, method, order))
    for n in (100, 10000):
        for method in ("linear", "cubic_spline", "pchip"):
            for order in ("sorted", "shuffled"):
                cases.append(oned(n, method, order, even=False))
    for n_eval in SCIPY["ppoly_bspline"]["n_eval"]:
        cases.append(piecewise("ppoly", n_eval))
        cases.append(piecewise("bspline", n_eval))
    for n in (100, 10000):
        cases.append(curve(n))

    spec = SCIPY["regular_grid"]
    for dim in spec["ndim"]:
        for size in spec["max_coord_size"]:
            # The sizes of SciPy's suite: each axis half the one before.
            shape = [size // 2**i for i in range(dim)]
            for n_eval in spec["n_samples"]:
                for flip in (False, True):
                    cases.append(grid(shape, "linear", n_eval, flip, f"max{size}"))
    for method in ("nearest", "pchip", "cubic"):
        for dim in (2, 3):
            cases.append(grid([40] * dim, method, 10000, False, "each40"))
    for per_axis in SCIPY["regular_grid_cubic_2d"]["n_per_axis"]:
        cases.append(grid([per_axis] * 2, "cubic", 1000, False, f"each{per_axis}"))
    for per_axis in SCIPY["regular_grid_quintic_2d"]["n_per_axis"]:
        cases.append(grid([per_axis] * 2, "quintic", 1000, False, f"each{per_axis}"))
    for per_axis in RUST_ONLY["regular_grid"]["n_per_axis"]:
        cases.append(grid([per_axis] * 2, "linear", 10000, False, f"each{per_axis}"))

    spec = SCIPY["nearest_nd"]
    for n in spec["n_samples"] + RUST_ONLY["nearest_nd"]["n_samples"]:
        for per_axis in spec["grid_per_axis"]:
            cases.append(nearest_nd(n, per_axis))
    for n in spec["n_samples"] + RUST_ONLY["linear_nd"]["n_samples"]:
        for per_axis in spec["grid_per_axis"]:
            cases.append(linear_nd(n, per_axis))
    for n in (1000, *RUST_ONLY["linear_nd"]["n_samples"]):
        cases.append(linear_nd_3d(n))

    spec = SCIPY["rbf_interpolator"]
    for n in spec["n_samples"]:
        for kernel in spec["kernels"]:
            for neighbors in spec["neighbors"]:
                cases.append(rbf(n, kernel, neighbors, "scattered"))
    for n in RUST_ONLY["rbf_global"]["n_samples"][:1]:
        cases.append(rbf(n, "thin_plate_spline", None, "scattered"))
    for n in RUST_ONLY["rbf_neighbors"]["n_samples"]:
        for queries in ("scattered", "image"):
            cases.append(rbf(n, "thin_plate_spline", RUST_ONLY["rbf_neighbors"]["neighbors"], queries))
    if real.available():
        cases += real_cases()
    names = [case["name"] for case in cases]
    assert len(names) == len(set(names)), "case names are unique"
    return cases


def write(case: dict, root: Path) -> Path:
    """Writes the case under `root`, and returns its directory."""
    directory = root / case["name"]
    directory.mkdir(parents=True, exist_ok=True)
    # The arrays of a case on real data are read here, when they are needed.
    arrays = case["arrays"]() if callable(case["arrays"]) else case["arrays"]
    for name, array in arrays.items():
        np.ascontiguousarray(array, dtype=np.float64).tofile(directory / f"{name}.f64")
    spec = {key: value for key, value in case.items() if key != "arrays"}
    spec["shapes"] = {name: list(np.shape(array)) for name, array in arrays.items()}
    if "points" in arrays:
        spec |= {"n": len(arrays["points"]), "n_eval": len(arrays["q"])}
    (directory / "case.json").write_text(json.dumps(spec, indent=1), encoding="utf-8")
    return directory
