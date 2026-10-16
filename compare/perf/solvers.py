"""Grid splines: SciPy's ways of finding the coefficients, against each other
and against interpx (P-13, #162).

    uv run --directory compare python -X utf8 perf/solvers.py [--filter REGEX] [--rounds N]

The ways, on the same grid, values, and queries (`cubic` and `quintic`):

- `scipy_iterative`: `RegularGridInterpolator(method=...)`, SciPy's default:
  one sparse system over all the axes, solved by `gcrotmk`.
- `scipy_spsolve`: the same system, solved directly (`solver=spsolve`).
- `scipy_legacy`: `method="..._legacy"`, which solves nothing at construction
  and fits one-dimensional splines along the axes at every evaluation.
- `scipy_axiswise`: a direct banded solve along each axis in turn
  (`make_interp_spline` with `axis=`), then `NdBSpline`. This is how interpx
  poses the problem, done by SciPy, and the reference of the distances: it is
  the tensor-product spline to rounding (docs/design.ja.md section 6.4).
- `interpx`: the comparison's interpx runner.

So the difference between the iterative and the direct solve of one system is
the solver's; between one system and the axes in turn, the posing's; between
SciPy's axes in turn and interpx, the implementation's.

Every run is a process of the comparison (`run.run_one`): one core, one
thread, its memory sampled from outside, stopped at the time limit. A run is
also stopped when the machine has less than `FREE_FLOOR` of memory left (one
system over all the axes of the elevation tile might take all of it), and is
reported as that.

Results: `out/solvers/results.jsonl`; the table: `out/solvers/table.md`,
printed at the end.
"""

import argparse
import json
import re
import sys
import threading
from pathlib import Path

import numpy as np
import psutil

HERE = Path(__file__).resolve().parent
OUT = HERE / "out" / "solvers"
WAYS = ("scipy_iterative", "scipy_spsolve", "scipy_legacy", "scipy_axiswise", "interpx")
REFERENCE = "scipy_axiswise"
METHODS = {"cubic": 3, "quintic": 5}
SIZES = {2: (10, 40, 100, 200, 400, 1000), 3: (10, 20, 40)}
N_EVAL = 10_000
FREE_FLOOR = 4 * 2**30
OUT_OF_MEMORY = f"stopped: less than {FREE_FLOOR // 2**30} GB of memory left on the machine"


def runner(way: str, directory: Path) -> None:
    """One SciPy way on one case, through the phases of the SciPy runner."""
    sys.path.insert(0, str(HERE))
    import scipy_runner  # says its process number first
    from scipy.interpolate import NdBSpline, RegularGridInterpolator, make_interp_spline
    from scipy.sparse.linalg import spsolve

    case = json.loads((directory / "case.json").read_text(encoding="utf-8"))

    def a(name: str) -> np.ndarray:
        return np.fromfile(directory / f"{name}.f64", dtype=np.float64)

    dim, method = case["dim"], case["method"]
    axes = tuple(a(f"axis{d}") for d in range(dim))
    values = a("values").reshape(case["shape"])
    q = a("q").reshape(-1, dim)

    def axiswise():
        coefficients, knots = values, []
        for d, x in enumerate(axes):
            spline = make_interp_spline(x, coefficients, k=METHODS[method], axis=d)
            knots.append(spline.t)
            # The spline keeps the axis it was fitted along first.
            coefficients = np.moveaxis(spline.c, 0, d)
        return NdBSpline(tuple(knots), coefficients, METHODS[method])

    build = {
        "scipy_iterative": lambda: RegularGridInterpolator(axes, values, method=method),
        "scipy_spsolve": lambda: RegularGridInterpolator(axes, values, method=method, solver=spsolve),
        "scipy_legacy": lambda: RegularGridInterpolator(axes, values, method=f"{method}_legacy"),
        "scipy_axiswise": axiswise,
    }[way]
    scipy_runner.phases(directory, way, build, lambda built: built(q))


def synthetic(method: str, dim: int, per_axis: int) -> dict:
    """A smooth function on an unevenly spaced grid, and queries in no order."""
    import cases as case_table

    rng = case_table.rng_for(f"solvers/{dim}d/{per_axis}")
    axes = [np.cumsum(rng.uniform(0.5, 1.5, per_axis)) for _ in range(dim)]
    unit = np.meshgrid(*[(x - x[0]) / (x[-1] - x[0]) for x in axes], indexing="ij")
    values = np.prod([np.cos(3.0 * u) for u in unit], axis=0) + sum(u * u for u in unit)
    q = np.stack([rng.uniform(x[0], x[-1], N_EVAL) for x in axes], axis=1)
    arrays = {f"axis{d}": x for d, x in enumerate(axes)} | {"values": values, "q": q}
    shape = [per_axis] * dim
    return {"name": f"{method}/{dim}d/each{per_axis}", "family": "grid", "method": method, "dim": dim, "shape": shape, "n": per_axis**dim, "n_eval": N_EVAL, "arrays": arrays}


def tile(method: str) -> dict:
    """The elevation tile of the comparison, as it is (docs/design.ja.md section 8.4.7)."""
    import cases as case_table

    real = next(case for case in case_table.real_cases() if case["name"] == f"grid/{method}/2d/{case_table.real.ELEVATION['name']}/heldout/ascending")
    return real | {"name": f"{method}/2d/{case_table.real.ELEVATION['name']}"}


def all_cases() -> list[dict]:
    import cases as case_table

    chosen = [synthetic(method, dim, n) for method in METHODS for dim, sizes in SIZES.items() for n in sizes]
    if case_table.real.available():
        chosen += [tile(method) for method in METHODS]
    return chosen


def apart(directory: Path, way: str) -> float | None:
    """How far the way's values are from the reference's, as a share of the largest of those."""
    ours, theirs = directory / f"values.{way}.f64", directory / f"values.{REFERENCE}.f64"
    if not ours.exists() or not theirs.exists():
        return None
    got, want = np.fromfile(ours), np.fromfile(theirs)
    if got.shape != want.shape or not np.all(np.isfinite(got)):
        return float("inf")
    return float(np.max(np.abs(got - want)) / (np.max(np.abs(want)) or 1.0))


def measure(case: dict, way: str, directory: Path) -> dict:
    """One run, stopped if the machine runs out of memory."""
    import run

    short = threading.Event()
    done = threading.Event()

    def watch() -> None:
        while not done.wait(0.2):
            if psutil.virtual_memory().available < FREE_FLOOR:
                short.set()
                for child in psutil.Process().children(recursive=True):
                    try:
                        child.kill()
                    except psutil.Error:
                        pass
                return

    watcher = threading.Thread(target=watch, daemon=True)
    watcher.start()
    try:
        result = run.run_one(way, directory, pause=False)
    finally:
        done.set()
        watcher.join()
    if short.is_set():
        return {"failed": OUT_OF_MEMORY, "seconds": result["seconds"]}
    return result


def table(rows: list[dict]) -> str:
    """Median over the rounds of each way on each case: construction, evaluation, peak memory, distance."""
    by = {}
    for row in rows:
        by.setdefault((row["case"], row["way"]), []).append(row)

    def seconds(value: float) -> str:
        for unit, scale in (("s", 1.0), ("ms", 1e-3), ("µs", 1e-6)):
            if value >= scale:
                return f"{value / scale:.3g} {unit}"
        return f"{value / 1e-9:.3g} ns"

    lines = ["| case | way | construction | evaluation | memory added (MB) | distance |", "| :--- | :--- | ---: | ---: | ---: | ---: |"]
    for case in dict.fromkeys(row["case"] for row in rows):
        for way in WAYS:
            got = by.get((case, way), [])
            done = [row for row in got if "build" in row]
            if not done:
                why = next((row["failed"] for row in got if "failed" in row), "not run")
                lines.append(f"| `{case}` | {way} | {why} | | | |")
                continue
            build = float(np.median([row["build"]["median"] for row in done]))
            evaluate = float(np.median([row["eval"]["median"] for row in done]))
            added = max((row["rss"]["overall_peak"] or 0) - (row["rss"]["loaded"] or 0) for row in done) / 2**20
            distance = done[0]["apart"]
            shown = "reference" if way == REFERENCE else ("" if distance is None else f"{distance:.1e}")
            lines.append(f"| `{case}` | {way} | {seconds(build)} | {seconds(evaluate)} | {added:.1f} | {shown} |")
    return "\n".join(lines) + "\n"


def main() -> None:
    if len(sys.argv) > 1 and sys.argv[1] == "--run":
        try:
            runner(sys.argv[2], Path(sys.argv[3]))
        except Exception as error:  # noqa: BLE001  (whatever SciPy raises is the result)
            print(json.dumps({"mark": "failed", "why": f"{type(error).__name__}: {error}"[:300]}), flush=True)
            sys.exit(1)
        return

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--filter", default="", help="run the cases whose name matches this regular expression")
    parser.add_argument("--rounds", type=int, default=3)
    args = parser.parse_args()

    sys.path.insert(0, str(HERE))
    import cases as case_table
    import run

    run.command = lambda way, directory: [str(run.executable("interpx")), str(directory)] if way == "interpx" else [sys.executable, "-X", "utf8", __file__, "--run", way, str(directory)]
    run.build_runners(["interpx"])
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "environment.json").write_text(json.dumps(run.environment(), indent=1), encoding="utf-8")
    chosen = [case for case in all_cases() if re.search(args.filter, case["name"])]
    results = OUT / "results.jsonl"
    have = set()
    if results.exists():
        for line in results.read_text(encoding="utf-8").splitlines():
            row = json.loads(line)
            have.add((row["case"], row["way"], row["round"]))
            if run.once(row):
                have |= {(row["case"], row["way"], later) for later in range(args.rounds)}

    with results.open("a", encoding="utf-8") as out:
        for round_ in range(args.rounds):
            for case in chosen:
                directory = OUT / "cases" / case["name"]
                if not (directory / "case.json").exists():
                    directory = case_table.write(case, OUT / "cases")
                # The reference first: the others are held against its values.
                for way in sorted(WAYS, key=lambda name: name != REFERENCE):
                    if (case["name"], way, round_) in have:
                        continue
                    measured = measure(case, way, directory)
                    row = {"case": case["name"], "way": way, "round": round_, "method": case["method"], "dim": case["dim"], "shape": case["shape"], "seconds": measured["seconds"]}
                    if "marks_ns" in measured:
                        row |= {"build": run.summary(measured["build_s"]), "eval": run.summary(measured["eval_s"]), "rss": run.memory(measured), "apart": apart(directory, way)}
                        shown = f"build {row['build']['median']:.3e} s  eval {row['eval']['median']:.3e} s  apart {row['apart']}"
                    else:
                        row["failed"] = measured.get("failed", "no result")
                        shown = row["failed"]
                    out.write(json.dumps(row) + "\n")
                    out.flush()
                    if run.once(row):
                        have |= {(case["name"], way, later) for later in range(args.rounds)}
                    print(f"[{round_ + 1}/{args.rounds}] {case['name']:28} {way:16} {shown}  {row['seconds']:.0f} s", flush=True)

    rows = [json.loads(line) for line in results.read_text(encoding="utf-8").splitlines()]
    text = table(rows)
    (OUT / "table.md").write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
