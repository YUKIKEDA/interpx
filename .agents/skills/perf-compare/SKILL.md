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

> The runners (`compare/perf/`, its Python harness, and the isolated Rust runner packages) are created by roadmap row M5-5. Until that row is `Done`, this skill states the contract. The `criterion` benches in `benches/` are created by the rows that need them. Do not create files ahead of their rows (`.cursor/rules/layout.mdc`).

## Layout

- `compare/perf/` holds the Python harness (SciPy side, timing, peak RSS) and one Rust runner package per other crate. Each Rust runner package has its own `[workspace]`, so its dependencies never reach the root crate.
- The other crates are chosen when M5-5 starts: maintained, pure Rust, and one per category where a comparable crate exists (1D, grid, scattered). A category without a comparable crate says "not applicable" in the table.
- Output goes to `compare/perf/out/` (ignored). Do not commit raw output; record the summary where the claim is made (PR, ADR, or docs).

## Method

1. Use the same data and queries for every library. The cases are fixed in JSON (sizes, dimensions, seed). Generate the data with a seeded NumPy `Generator`.
2. Build every side with the release profile and the same flags. Stop leftover processes from earlier runs.
3. Time a warmed-up repeat of the evaluation, and report the median and the min and max. A difference inside the spread is not a result.
4. Time build (construction) and evaluation separately. For the batch API, time `eval_batch_into` against SciPy's vectorized call, and the `_par` form against the serial form with the thread count named.
5. Measure memory as the peak (the Windows peak working set, `VmHWM` on Linux, `ru_maxrss` elsewhere) inside the runner, and for Rust runners also the counting-allocator peak, build and evaluate separately. Include RBF at several `N` to show the O(N²) growth.
6. Record: machine, OS, `rustc` version, crate versions, SciPy and NumPy versions, the command, and the commit.

## Named baselines

`criterion` baselines are saved per milestone (`--save-baseline m1`, `m2`, ...). A hot-path PR compares against the latest baseline and pastes the comparison under `## Verification`. The log of baselines is local (`.dev/bench-log.md`, untracked).

## When the result is surprising

Profile before changing code: `.agents/skills/cpu-profile/SKILL.md`. Write on the Issue where the time goes.
