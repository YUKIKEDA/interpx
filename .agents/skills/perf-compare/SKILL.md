---
name: perf-compare
description: >-
  Run interpx's manual comparison against SciPy 1.18 and other pure-Rust
  interpolation crates, for time and memory. Use when a speed or memory claim
  needs numbers, when a hot-path PR needs a baseline, or when the comparison
  tables must be refreshed.
---

# Performance comparison

The rules are in `.cursor/rules/bench.mdc` and `.cursor/rules/perf.mdc`. The design is `docs/design.ja.md` §8.3. The comparison is manual and is not a CI gate.

## Layout

- `compare/perf/cases.py`: the cases and their data. The sizes are those of SciPy's own suite (`compare/common/bench_cases.py`, `docs/design.ja.md` §8.4.3) and a few larger ones. A case is written once (`out/cases/NAME/`: `case.json` and one file of 64-bit floats for each array), and every library reads the same files.
- `compare/perf/scipy_runner.py`: the SciPy side.
- `compare/perf/runners/`: one Cargo package for each Rust library (`interpx`, `scirs2`, `ndarray_interp`, `interpn`), each with its own `[workspace]`, so that its dependencies never reach the root crate. `harness/` is the part they share: reading a case, the phases, the timing.
- `compare/perf/run.py`: the driver. `compare/perf/report.py`: writes `docs/performance.md` and `docs/performance/*.svg` from a run.
- Output goes to `compare/perf/out/` (ignored). Raw output is not committed; the summary is.

## Run

```text
just compare-perf                       # every case, three rounds
just compare-perf --filter "^grid/"     # the cases whose name matches
just compare-perf --list                # the names
just compare-report                     # docs/performance.md from what out/ holds
```

How long a run takes is set by three things, not by the work measured: the number of processes (cases x libraries x rounds), the time each timed loop is given (50 ms), and the start of Python and SciPy (0.7 s a process). The first runs gave each loop 250 ms, paused 60 ms after every phase of every case, and started every library on every case to hear that it had none: all cases took two hours. Write the time of the next full run here when it is measured.

A run that was stopped goes on where it stopped (`out/results.jsonl` says what is done). Delete `compare/perf/out/` to start over. Nothing else that takes a core runs meanwhile: no build, no test, no second comparison.

## Method

1. One process for one case and one library, pinned to one core, with the thread counts of BLAS and `rayon` set to one. The libraries take turns case by case, and the list is gone through three times.
2. Every runner goes through the same phases and prints a line at the end of each: the data is read, the interpolator is built once, the queries are evaluated once and the values written, then construction and evaluation are each repeated and timed.
3. A time is that of one call. A call shorter than 200 µs is repeated inside one sample (the clock has steps of 100 ns). The report gives the median of the samples of a process, then the median of the rounds, and counts interpx as slower when its median is more than 1.10 times another's and its fastest sample is slower than the other's fastest.
4. Memory is the resident set size of the runner, sampled every millisecond by the driver from outside, so that it is measured the same way for Python and Rust and is a series over time. The report gives it at the phases (data read, while building, built, while evaluating) and draws it, for the large cases named in `cases.py` (`MEMORY`). A runner pauses after each phase of those cases only, so that the phase shows in the samples. The counting allocator of `tests/alloc.rs` stays the tool for "does this allocate"; it is not used across libraries.
5. The values every library writes are compared with SciPy's. Where they are apart by more than rounding, the library computes another interpolant under the same name (another end condition, a local cubic for a spline): the report says so next to its time, and such a time is not a loss to answer.
6. The report records the machine, the OS, the versions, the commit, and the command.

## A new library, a new type

- A library: its cases in `supports` (`run.py`), so that it is not started for the rest; a package under `runners/` with `perf-harness` as a path dependency and a `main` that calls `run` for the cases it has and `unsupported` for the rest; its name in `RUST_RUNNERS` (`run.py`) and in `LIBRARIES` and `TITLES` (`report.py`). Chosen when maintained, pure Rust, and computing something this crate computes; say in the PR why it and not its neighbors.
- A type (the rows of M4): its cases in `cases.py`, an arm in each runner that has it, and its family in `FAMILIES` (`report.py`).

## What a loss becomes

Each comparison where interpx is slower than a library that computes the same thing is an Issue and a roadmap row before the features that sit on top of the types, the largest ratio first (`.cursor/rules/bench.mdc`). Profile before changing code.

## Named baselines

`criterion` baselines are saved per milestone (`--save-baseline m1`, `m2`, ...). A hot-path PR compares against the latest baseline and pastes the comparison under `## Verification`. The log of baselines is local (`.dev/bench-log.md`, untracked).

## When the result is surprising

Profile before changing code: `.agents/skills/cpu-profile/SKILL.md`. Write on the Issue where the time goes.
