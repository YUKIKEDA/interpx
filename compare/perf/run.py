"""Runs the comparison: every case with every library, one process each.

    uv run --directory compare python -X utf8 perf/run.py [--filter REGEX]
        [--rounds N] [--libraries a,b] [--list]

For each case the data is written once (`out/cases/NAME/`). Each library then
runs the case in a process of its own, pinned to one core, and reports its
times. While it runs, this script samples the resident set size of that
process, so that memory is a series over time with the phases marked on it
(data read, built, evaluated), measured the same way for every library.

The libraries take turns case by case, and the whole list is gone through
`--rounds` times, so that a slow minute of the machine does not fall on one
library alone.

Results: `out/results.jsonl` (a line for each run), `out/series/` (the memory
samples), `out/environment.json`. `out/` is not committed; `report.py` makes
the summary that is.
"""

import argparse
import json
import os
import platform
import re
import subprocess
import sys
import threading
import time
from pathlib import Path

import numpy as np
import psutil

import cases as case_table

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
# Other processes start on the first cores; the runner gets one further on.
CORE = min(4, (os.cpu_count() or 1) - 1)
SAMPLE_EVERY = 0.001
# One library on one case: longer than this is a fault, not a result.
TIME_LIMIT = 900.0
# How a run that the time limit stopped is recorded (`failed`).
STOPPED = "stopped at the time limit"
# A run whose process takes this long is named in the report as near the
# limit. The same run took 760 s alone and was stopped at 900 s inside a run
# of the whole comparison (#168): the state of the machine moves a run by a
# fifth or more, and by how much more is not known. Half the limit leaves it
# twice its time.
NEAR = TIME_LIMIT / 2

RUST_RUNNERS = ("interpx", "scirs2", "ndarray_interp", "interpn", "spade")
# interpx again, built with `from_owned` (the data handed over, not copied).
# and with its parallel batch, on `PAR_THREADS` threads and every core.
VARIANTS = {"interpx_owned": ("interpx", "owned"), "interpx_par": ("interpx", "par")}
PAR_THREADS = min(8, os.cpu_count() or 1)
LIBRARIES = ("scipy", *RUST_RUNNERS, *VARIANTS)


def executable(runner: str) -> Path:
    suffix = ".exe" if os.name == "nt" else ""
    return HERE / "runners" / runner / "target" / "release" / f"perf-{runner}{suffix}"


def build_runners(libraries: list[str]) -> None:
    for runner in RUST_RUNNERS:
        if runner in libraries or any(VARIANTS.get(name, ("",))[0] == runner for name in libraries):
            manifest = HERE / "runners" / runner / "Cargo.toml"
            subprocess.run(["cargo", "build", "--release", "--manifest-path", str(manifest)], check=True)


def supports(library: str, case: dict) -> bool:
    """Whether the library has the case. A runner says so itself as well
    (`unsupported`); this saves starting it to hear that."""
    family, method = case["family"], case.get("method")
    ascending = not case["name"].endswith("/flipped")
    if library == "interpx_owned":
        return family in ("grid", "nearest_nd", "linear_nd", "rbf")
    if library == "spade":
        return family == "linear_nd" and case.get("dim", 2) == 2
    if family == "linear_nd":
        return library in ("scipy", "interpx", "interpx_par") or (library == "scirs2" and case.get("dim", 2) == 2)
    if library == "interpn":
        return family == "grid" and ascending and method in ("linear", "nearest", "cubic")
    if library == "ndarray_interp":
        if family == "oned":
            return method in ("linear", "cubic_spline")
        return family == "grid" and ascending and method == "linear" and case["dim"] == 2
    if library == "scirs2":
        if family == "oned":
            return method in ("linear", "nearest", "cubic_spline", "pchip", "akima")
        if family == "grid":
            return ascending and method in ("linear", "nearest")
        return family == "rbf" and case["neighbors"] is None and case["kernel"] != "inverse_quadratic"
    return True


def command(library: str, directory: Path) -> list[str]:
    if library == "scipy":
        return [sys.executable, "-X", "utf8", str(HERE / "scipy_runner.py"), str(directory)]
    runner, *arguments = VARIANTS.get(library, (library,))
    return [str(executable(runner)), str(directory), *arguments]


def run_one(library: str, directory: Path, pause: bool) -> dict:
    """One process: its marks, its times, its memory over time, and how long
    the whole process took (`seconds`, whatever its outcome)."""
    begun = time.perf_counter()
    run = measure_one(library, directory, pause)
    return run | {"seconds": time.perf_counter() - begun}


def measure_one(library: str, directory: Path, pause: bool) -> dict:
    several = library == "interpx_par"
    env = os.environ | ({"PERF_PAUSE_MS": str(case_table.PAUSE_MS)} if pause else {}) | {name: "1" for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "RAYON_NUM_THREADS")}
    if several:
        env["RAYON_NUM_THREADS"] = str(PAR_THREADS)
    start_ns = time.time_ns()
    process = subprocess.Popen(command(library, directory), stdout=subprocess.PIPE, text=True, env=env)
    # One core for the serial runs; the parallel one has them all.
    cores = list(range(os.cpu_count() or 1)) if several else [CORE]
    # The process that does the work, whose memory is sampled: the one that
    # was started, until a runner says it is another (`pid`, below).
    watched = [psutil.Process(process.pid)]
    pinned: list[str] = []

    def pin(worker: psutil.Process) -> None:
        """Puts the process on the cores, and reads back where it is: a run
        that is not where the report says it was is not a result."""
        try:
            worker.cpu_affinity(cores)
            at = worker.cpu_affinity()
        except (psutil.Error, AttributeError) as error:
            pinned.append(f"the cores could not be set: {error}")
            return
        if sorted(at) != cores:
            pinned.append(f"on cores {at}, not {cores}")

    pin(watched[0])
    samples: list[tuple[int, int]] = []
    done = threading.Event()

    def sample() -> None:
        while not done.is_set():
            try:
                samples.append((time.time_ns() - start_ns, watched[0].memory_info().rss))
            except psutil.Error:
                return
            time.sleep(SAMPLE_EVERY)

    sampler = threading.Thread(target=sample, daemon=True)
    sampler.start()
    marks, result = {}, {}
    stopped = threading.Event()

    def stop() -> None:
        stopped.set()
        # The worker first: killing a launcher leaves its child running.
        for running in (watched[0], process):
            try:
                running.kill()
            except (psutil.Error, OSError):
                pass

    killer = threading.Timer(TIME_LIMIT, stop)
    killer.start()
    try:
        for line in process.stdout:
            line = line.strip()
            if not line.startswith("{"):
                continue
            message = json.loads(line)
            if message["mark"] == "result":
                result = message
            elif message["mark"] == "unsupported":
                result = {"unsupported": message["why"]}
            elif message["mark"] == "failed":
                result = {"failed": message["why"]}
            elif message["mark"] == "pid":
                if message["pid"] != watched[0].pid:
                    try:
                        worker = psutil.Process(message["pid"])
                    except psutil.Error as error:
                        pinned.append(f"the worker is gone: {error}")
                        continue
                    pin(worker)
                    # What was sampled so far is the launcher's.
                    samples.clear()
                    watched[0] = worker
            else:
                marks[message["mark"]] = message["t_ns"] - start_ns
        code = process.wait()
    finally:
        killer.cancel()
        done.set()
        sampler.join()
    if stopped.is_set():
        return {"failed": f"{STOPPED} of {TIME_LIMIT / 60:.0f} minutes"}
    if "failed" in result:
        return result
    if pinned and "unsupported" not in result:
        return {"failed": f"not measured where it should be: {pinned[0]}"}
    if code != 0 or not result:
        return {"failed": f"exit code {code}"}
    if "unsupported" in result:
        return result
    return {"marks_ns": marks, "samples": samples} | result


def memory(run: dict) -> dict:
    """The resident set size at the phases, in bytes."""
    marks, samples = run["marks_ns"], run["samples"]

    def peak(after: int, before: int) -> int | None:
        inside = [rss for at, rss in samples if after <= at <= before]
        return max(inside) if inside else None

    loaded, built, evaluated, timed = (marks[name] for name in ("loaded", "built", "evaluated", "timed"))
    held = [rss for at, rss in samples if at >= built]
    return {
        # Just before construction starts: the process and the data.
        "loaded": next((rss for at, rss in reversed(samples) if at <= loaded), None),
        "build_peak": peak(loaded, built),
        # With the interpolator built and nothing else going on.
        "built": held[0] if held else None,
        "eval_peak": peak(built, evaluated),
        "overall_peak": peak(0, timed),
    }


def summary(seconds: list[float]) -> dict:
    return {"median": float(np.median(seconds)), "min": min(seconds), "max": max(seconds), "repeats": len(seconds)}


def agreement(directory: Path, library: str) -> float | None:
    """How far the library's values are from SciPy's, as a share of the largest value."""
    ours, theirs = directory / f"values.{library}.f64", directory / "values.scipy.f64"
    if library == "scipy" or not ours.exists() or not theirs.exists():
        return None
    a, b = np.fromfile(ours), np.fromfile(theirs)
    if a.shape != b.shape:
        return float("inf")
    if not np.array_equal(np.isnan(a), np.isnan(b)):
        return float("inf")
    scale = np.nanmax(np.abs(b)) or 1.0
    return float(np.nanmax(np.abs(a - b)) / scale) if a.size else 0.0


def held_out(directory: Path, library: str) -> dict | None:
    """For a case on real data: how far the library's values are from the
    values that were held out, and how many queries it gave no value."""
    ours, truth = directory / f"values.{library}.f64", directory / "truth.f64"
    if not ours.exists() or not truth.exists():
        return None
    got, want = np.fromfile(ours), np.fromfile(truth)
    if got.shape != want.shape:
        return {"shape": list(got.shape)}
    has = np.isfinite(got)
    off = got[has] - want[has]
    if not off.size:
        return {"without_a_value": int(got.size)}
    return {
        "max": float(np.abs(off).max()),
        "rms": float(np.sqrt(np.mean(off * off))),
        "without_a_value": int((~has).sum()),
    }


def once(row: dict) -> bool:
    """Whether a library that gave this row on a case is not run on it again
    in a later round: it was stopped at the time limit, or one of its calls
    was timed by one call (the runner does that for a call of two seconds or
    more, `LONG` in the runners).

    The rounds are there for the spread between processes, which for a call
    that long is a few percent (three constructions of 107 s read 106 to
    124 s); run three times, the cases on real data took four hours, most of
    it one library solving one system again. The rule goes by the time and
    not by the library. It goes by how the runner timed the call and not by
    the median it reported: a call near two seconds that was timed several
    times read 1.75 to 1.99 s in one round and 2.2 to 3.7 s in the next, and
    a rule on the median left the round after out, so that SciPy's time on
    that case was the median of one usual round and one slow one (#171)."""
    if str(row.get("failed", "")).startswith(STOPPED):
        return True
    return any(row[what]["repeats"] == 1 for what in ("build", "eval") if what in row)


def environment() -> dict:
    def version(*args: str) -> str:
        try:
            return subprocess.run(args, capture_output=True, text=True, check=True).stdout.strip()
        except (OSError, subprocess.CalledProcessError):
            return "unknown"

    import scipy

    return {
        "machine": platform.processor() or platform.machine(),
        "cores": os.cpu_count(),
        "memory_bytes": psutil.virtual_memory().total,
        "os": platform.platform(),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "rustc": version("rustc", "--version"),
        "commit": version("git", "rev-parse", "HEAD"),
        "dirty": bool(version("git", "status", "--porcelain", "--untracked-files=no")),
        "core": CORE,
        "par_threads": PAR_THREADS,
        "date": time.strftime("%Y-%m-%d"),
        "command": " ".join(sys.argv),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--filter", default="", help="run the cases whose name matches this regular expression")
    parser.add_argument("--rounds", type=int, default=3)
    parser.add_argument("--libraries", default=",".join(LIBRARIES))
    parser.add_argument("--list", action="store_true", help="print the names of the cases and stop")
    args = parser.parse_args()

    chosen = [case for case in case_table.all_cases() if re.search(args.filter, case["name"])]
    if not case_table.real.available():
        print("The real-world data is not fetched (`just compare-fetch`): its cases are left out.", flush=True)
    if args.list:
        print("\n".join(case["name"] for case in chosen))
        print(f"{len(chosen)} cases")
        return
    libraries = args.libraries.split(",")
    build_runners(libraries)
    OUT.mkdir(exist_ok=True)
    (OUT / "series").mkdir(exist_ok=True)
    (OUT / "environment.json").write_text(json.dumps(environment(), indent=1), encoding="utf-8")
    directories = {case["name"]: case_table.write(case, OUT / "cases") for case in chosen}

    # A run that was stopped goes on where it stopped: what `results.jsonl`
    # has is not run again. Delete `out/` to start over.
    have, stopped, again = set(), set(), set()
    if (OUT / "results.jsonl").exists():
        for line in (OUT / "results.jsonl").read_text(encoding="utf-8").splitlines():
            row = json.loads(line)
            have.add((row["case"], row["library"], row["round"]))
            if row.get("again"):
                again.add((row["case"], row["library"]))
            elif str(row.get("failed", "")).startswith(STOPPED):
                stopped.add((row["case"], row["library"]))
            if once(row):
                have |= {(row["case"], row["library"], later) for later in range(args.rounds)}

    def record(case: dict, library: str, round_: int, results, where: str, again: bool = False) -> dict:
        """Runs one library on one case, and writes its row (`where` says
        how far the run has gone, for the progress line)."""
        directory = directories[case["name"]]
        run = run_one(library, directory, case["name"] in case_table.MEMORY)
        row = {"case": case["name"], "library": library, "round": round_} | ({"again": True} if again else {})
        row |= {key: value for key, value in case.items() if key not in ("arrays", "name")}
        row["seconds"] = run["seconds"]
        if "marks_ns" in run:
            row |= {
                "build": summary(run["build_s"]),
                "eval": summary(run["eval_s"]),
                "rss": memory(run),
                "marks_ns": run["marks_ns"],
                "apart": agreement(directory, library),
            }
            if "real" in case:
                row["error"] = held_out(directory, library)
                # What the case table could not say before the data was read.
                row |= {key: json.loads((directory / "case.json").read_text(encoding="utf-8"))[key] for key in ("n", "n_eval")}
            series = OUT / "series" / f"{case['name'].replace('/', '__')}.{library}.{round_}.csv"
            series.write_text("".join(f"{at},{rss}\n" for at, rss in run["samples"]), encoding="utf-8")
            shown = f"build {row['build']['median']:.3e} s  eval {row['eval']['median']:.3e} s"
        else:
            row |= {key: value for key, value in run.items() if key != "seconds"}
            shown = str(run)
        results.write(json.dumps(row) + "\n")
        results.flush()
        near = "  (near the time limit)" if row["seconds"] >= NEAR else ""
        print(f"[{where}] {case['name']:58} {library:10} {shown}  {row['seconds']:.0f} s{near}", flush=True)
        return row

    with (OUT / "results.jsonl").open("a", encoding="utf-8") as results:
        for round_ in range(args.rounds):
            for number, case in enumerate(chosen, 1):
                # SciPy first: the others are compared with its values.
                for library in sorted(libraries, key=lambda name: name != "scipy"):
                    if (case["name"], library, round_) in have or not supports(library, case):
                        continue
                    row = record(case, library, round_, results, f"{round_ + 1}/{args.rounds} {number}/{len(chosen)}")
                    if once(row):
                        have |= {(case["name"], library, later) for later in range(args.rounds)}
                    if str(row.get("failed", "")).startswith(STOPPED):
                        stopped.add((case["name"], library))
        # A run that the time limit stopped is run once more, alone, after
        # all the others: the report calls it stopped only if it is stopped
        # again. Whether one run of #166 finished was decided by the state of
        # the machine (#168).
        for case in chosen:
            for library in libraries:
                if (case["name"], library) in stopped and (case["name"], library) not in again:
                    record(case, library, args.rounds, results, "again", again=True)


if __name__ == "__main__":
    main()
