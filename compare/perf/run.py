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

RUST_RUNNERS = ("interpx", "scirs2", "ndarray_interp", "interpn")
LIBRARIES = ("scipy", *RUST_RUNNERS)


def executable(runner: str) -> Path:
    suffix = ".exe" if os.name == "nt" else ""
    return HERE / "runners" / runner / "target" / "release" / f"perf-{runner}{suffix}"


def build_runners(libraries: list[str]) -> None:
    for runner in RUST_RUNNERS:
        if runner in libraries:
            manifest = HERE / "runners" / runner / "Cargo.toml"
            subprocess.run(["cargo", "build", "--release", "--manifest-path", str(manifest)], check=True)


def command(library: str, directory: Path) -> list[str]:
    if library == "scipy":
        return [sys.executable, "-X", "utf8", str(HERE / "scipy_runner.py"), str(directory)]
    return [str(executable(library)), str(directory)]


def run_one(library: str, directory: Path) -> dict:
    """One process: its marks, its times, and its memory over time."""
    env = os.environ | {name: "1" for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "RAYON_NUM_THREADS")}
    start_ns = time.time_ns()
    process = subprocess.Popen(command(library, directory), stdout=subprocess.PIPE, text=True, env=env)
    watched = psutil.Process(process.pid)
    try:
        watched.cpu_affinity([CORE])
    except (psutil.Error, AttributeError):
        pass
    samples: list[tuple[int, int]] = []
    done = threading.Event()

    def sample() -> None:
        while not done.is_set():
            try:
                samples.append((time.time_ns() - start_ns, watched.memory_info().rss))
            except psutil.Error:
                return
            time.sleep(SAMPLE_EVERY)

    sampler = threading.Thread(target=sample, daemon=True)
    sampler.start()
    marks, result = {}, {}
    killer = threading.Timer(TIME_LIMIT, process.kill)
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
            else:
                marks[message["mark"]] = message["t_ns"] - start_ns
        code = process.wait()
    finally:
        killer.cancel()
        done.set()
        sampler.join()
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
    have = set()
    if (OUT / "results.jsonl").exists():
        for line in (OUT / "results.jsonl").read_text(encoding="utf-8").splitlines():
            row = json.loads(line)
            have.add((row["case"], row["library"], row["round"]))

    with (OUT / "results.jsonl").open("a", encoding="utf-8") as results:
        for round_ in range(args.rounds):
            for number, case in enumerate(chosen, 1):
                directory = directories[case["name"]]
                # SciPy first: the others are compared with its values.
                for library in sorted(libraries, key=lambda name: name != "scipy"):
                    if (case["name"], library, round_) in have:
                        continue
                    run = run_one(library, directory)
                    row = {"case": case["name"], "library": library, "round": round_}
                    row |= {key: value for key, value in case.items() if key not in ("arrays", "name")}
                    if "marks_ns" in run:
                        row |= {
                            "build": summary(run["build_s"]),
                            "eval": summary(run["eval_s"]),
                            "rss": memory(run),
                            "marks_ns": run["marks_ns"],
                            "apart": agreement(directory, library),
                        }
                        series = OUT / "series" / f"{case['name'].replace('/', '__')}.{library}.{round_}.csv"
                        series.write_text("".join(f"{at},{rss}\n" for at, rss in run["samples"]), encoding="utf-8")
                        shown = f"build {row['build']['median']:.3e} s  eval {row['eval']['median']:.3e} s"
                    else:
                        row |= run
                        shown = str(run)
                    results.write(json.dumps(row) + "\n")
                    results.flush()
                    print(f"[{round_ + 1}/{args.rounds} {number}/{len(chosen)}] {case['name']:58} {library:10} {shown}", flush=True)


if __name__ == "__main__":
    main()
