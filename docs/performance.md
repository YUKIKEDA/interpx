# Performance: interpx against SciPy and other Rust crates

<!-- Written by compare/perf/report.py from a run of compare/perf/run.py. Do not edit by hand: run it again. -->

Measured on 2026-10-13, at commit `f82956ec63` with uncommitted changes.

## How it was measured

- Machine: Intel64 Family 6 Model 191 Stepping 2, GenuineIntel, 16 logical cores, 48 GB. Windows-11-10.0.26300-SP0.
- rustc 1.97.1 (8bab26f4f 2026-07-14); Python 3.13.1, NumPy 2.5.4, SciPy 1.18.1.
- The other Rust crates, at the versions in `compare/perf/runners/*/Cargo.toml`: `scirs2-interpolate` (the widest of the Rust crates: one variable, grids, RBF), `ndarray-interp` (the most used for one variable; also bilinear), `interpn` (grids).
- Every library ran each case in a process of its own, pinned to one core, on one thread, on the same data. The libraries took turns case by case, 3 rounds.
- A time is that of one call: construction, or the evaluation of all the queries of the case (`eval_batch_into` for interpx, the call with an array for the others). It is the median of the repeats of a process, and then the median of the rounds.
- interpx is called slower where its median is more than 1.10 times another library's and its fastest repeat is slower than that library's fastest.
- Memory is the resident set size of the process, sampled every millisecond from outside it.
- Command: `perf/run.py --rounds 3`.

## Where interpx is slower

| family | operation | cases | interpx fastest | slower than SciPy | slower than a Rust crate |
| :--- | :--- | ---: | ---: | ---: | ---: |
| One variable | evaluation | 100 | 96 | 0 | 4 |
| One variable | construction | 100 | 70 | 1 | 30 |
| `PPoly` | evaluation | 4 | 4 | 0 | 0 |
| `PPoly` | construction | 4 | 4 | 0 | 0 |
| `BSpline` | evaluation | 4 | 4 | 0 | 0 |
| `BSpline` | construction | 4 | 4 | 0 | 0 |
| `BSplineCurve` | evaluation | 2 | 2 | 0 | 0 |
| `BSplineCurve` | construction | 2 | 2 | 0 | 0 |
| Regular grid | evaluation | 79 | 53 | 0 | 26 |
| Regular grid | construction | 79 | 19 | 36 | 43 |
| `NearestNDInterpolator` | evaluation | 15 | 15 | 0 | 0 |
| `NearestNDInterpolator` | construction | 15 | 15 | 0 | 0 |
| `RBFInterpolator` | evaluation | 53 | 52 | 1 | 0 |
| `RBFInterpolator` | construction | 53 | 53 | 0 | 0 |

192 comparisons (a case, an operation, and a library), by method. "Same values" says whether that library's values are SciPy's up to rounding: where they are not, it computes another interpolant under the same name, and its time is not that of the same work.

| method | operation | against | same values | cases | slower by | the largest |
| :--- | :--- | :--- | :--- | ---: | ---: | :--- |
| `grid/linear` | construction | interpn | yes | 34 | 181x to 4.4e+07x | `grid/linear/4d/max200/eval100/ascending`: 49 ms against 1.11 ns |
| `grid/cubic` | construction | interpn | no | 7 | 1.63e+03x to 1.61e+06x | `grid/cubic/2d/each400/eval1000/ascending`: 1.44 ms against 0 |
| `grid/nearest` | construction | interpn | yes | 2 | 1.05e+03x to 2.97e+04x | `grid/nearest/3d/each40/eval10000/ascending`: 26.4 µs against 0 |
| `grid/linear` | construction | SciPy | yes | 34 | 1.24x to 5.84e+03x | `grid/linear/4d/max200/eval100/flipped`: 120 ms against 20.5 µs |
| `grid/linear` | construction | ndarray-interp | yes | 14 | 1.39x to 522x | `grid/linear/2d/each800/eval10000/ascending`: 1.17 ms against 2.24 µs |
| `grid/pchip` | construction | SciPy | yes | 1 | 52.6x | `grid/pchip/3d/each40/eval10000/ascending`: 507 µs against 9.65 µs |
| `oned/linear` | construction | ndarray-interp | yes | 10 | 1.59x to 4.97x | `oned/linear/n10000/shuffled`: 86.4 µs against 17.4 µs |
| `oned/akima` | construction | scirs2-interpolate | yes | 10 | 1.57x to 4.62x | `oned/akima/n1000/shuffled`: 21.1 µs against 4.56 µs |
| `grid/linear` | construction | scirs2-interpolate | yes | 27 | 1.12x to 4.35x | `grid/linear/4d/max40/eval10/ascending`: 23.1 µs against 5.32 µs |
| `oned/linear` | construction | scirs2-interpolate | yes | 8 | 1.33x to 3.63x | `oned/linear/n10000/shuffled`: 86.4 µs against 23.8 µs |
| `grid/nearest` | construction | SciPy | yes | 1 | 3.18x | `grid/nearest/3d/each40/eval10000/ascending`: 26.4 µs against 8.31 µs |
| `grid/nearest` | construction | scirs2-interpolate | yes | 2 | 2.42x to 2.78x | `grid/nearest/3d/each40/eval10000/ascending`: 26.4 µs against 9.51 µs |
| `oned/linear` | evaluation | ndarray-interp | yes | 2 | 1.49x to 2.42x | `oned/linear/n10000/shuffled`: 633 µs against 262 µs |
| `oned/cubic_spline` | evaluation | ndarray-interp | yes | 2 | 1.35x to 2.18x | `oned/cubic_spline/n10000/shuffled`: 753 µs against 345 µs |
| `oned/pchip` | construction | scirs2-interpolate | yes | 6 | 1.34x to 2.1x | `oned/pchip/n1000/sorted`: 9.51 µs against 4.52 µs |
| `oned/linear` | construction | SciPy | yes | 1 | 1.77x | `oned/linear/n10000/shuffled`: 86.4 µs against 48.7 µs |
| `oned/nearest` | construction | scirs2-interpolate | yes | 2 | 1.56x to 1.69x | `oned/nearest/n10000/shuffled`: 39.6 µs against 23.4 µs |
| `grid/linear` | evaluation | interpn | yes | 26 | 1.12x to 1.59x | `grid/linear/4d/max40/eval10/ascending`: 428 ns against 269 ns |
| `oned/cubic_spline` | construction | scirs2-interpolate | no | 2 | 1.16x to 1.21x | `oned/cubic_spline/n1000/shuffled`: 19.5 µs against 16.1 µs |
| `rbf/quintic` | evaluation | SciPy | yes | 1 | 1.12x | `rbf/quintic/global/n10/scattered`: 75 µs against 66.8 µs |

<details><summary>Every comparison where interpx is slower</summary>

"Values apart" is the distance of that library's values from SciPy's, as a share of the largest value.

| case | operation | against | interpx | the other | slower by | values apart |
| :--- | :--- | :--- | ---: | ---: | ---: | ---: |
| `grid/linear/4d/max200/eval100/ascending` | construction | interpn | 49 ms | 1.11 ns | 44023538.67x | 3e-16 |
| `grid/linear/4d/max200/eval1000/ascending` | construction | interpn | 50.3 ms | 1.15 ns | 43682162.82x | 3e-16 |
| `grid/linear/4d/max200/eval10000/ascending` | construction | interpn | 48.4 ms | 1.11 ns | 43469894.14x | 4e-16 |
| `grid/linear/4d/max200/eval10/ascending` | construction | interpn | 47 ms | 1.12 ns | 41813524.15x | 2e-16 |
| `grid/linear/3d/max200/eval10/ascending` | construction | interpn | 1.86 ms | 0 | 2091933.81x | 3e-16 |
| `grid/linear/3d/max200/eval100/ascending` | construction | interpn | 1.89 ms | 0 | 2088133.54x | 3e-16 |
| `grid/linear/3d/max200/eval10000/ascending` | construction | interpn | 1.82 ms | 0 | 2044296.77x | 5e-16 |
| `grid/linear/3d/max200/eval1000/ascending` | construction | interpn | 1.8 ms | 0 | 2023648.44x | 3e-16 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | interpn | 1.17 ms | 0 | 1751965.81x | 3e-16 |
| `grid/cubic/2d/each400/eval1000/ascending` | construction | interpn | 1.44 ms | 0 | 1611635.73x | 5e-01 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | interpn | 329 µs | 0 | 480012.84x | 3e-16 |
| `grid/cubic/3d/each40/eval10000/ascending` | construction | interpn | 341 µs | 0 | 352910.06x | 5e-01 |
| `grid/cubic/2d/each200/eval1000/ascending` | construction | interpn | 190 µs | 0 | 213553.10x | 2e-01 |
| `grid/cubic/2d/each100/eval1000/ascending` | construction | interpn | 58.2 µs | 0 | 65311.56x | 2e-01 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | interpn | 26.4 µs | 0 | 29667.95x | 0e+00 |
| `grid/linear/4d/max40/eval10/ascending` | construction | interpn | 23.1 µs | 1.11 ns | 20771.77x | 2e-16 |
| `grid/linear/4d/max40/eval10000/ascending` | construction | interpn | 23.1 µs | 1.11 ns | 20732.49x | 5e-16 |
| `grid/linear/4d/max40/eval100/ascending` | construction | interpn | 23.2 µs | 1.14 ns | 20267.01x | 2e-16 |
| `grid/linear/4d/max40/eval1000/ascending` | construction | interpn | 22.9 µs | 1.16 ns | 19706.57x | 6e-16 |
| `grid/linear/2d/max200/eval100/ascending` | construction | interpn | 8.13 µs | 0 | 12161.75x | 1e-16 |
| `grid/linear/2d/max200/eval10000/ascending` | construction | interpn | 8.2 µs | 0 | 12100.09x | 3e-16 |
| `grid/linear/2d/max200/eval10/ascending` | construction | interpn | 7.7 µs | 0 | 11521.17x | 2e-16 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | interpn | 7.6 µs | 0 | 11330.95x | 2e-16 |
| `grid/cubic/2d/each40/eval1000/ascending` | construction | interpn | 8.87 µs | 0 | 9952.44x | 7e-01 |
| `grid/cubic/2d/each40/eval10000/ascending` | construction | interpn | 9.03 µs | 0 | 9233.57x | 7e-01 |
| `grid/linear/4d/max200/eval100/flipped` | construction | SciPy | 120 ms | 20.5 µs | 5841.49x |  |
| `grid/linear/4d/max200/eval10/flipped` | construction | SciPy | 118 ms | 20.8 µs | 5683.82x |  |
| `grid/linear/4d/max200/eval1000/flipped` | construction | SciPy | 116 ms | 20.5 µs | 5670.61x |  |
| `grid/linear/3d/max40/eval100/ascending` | construction | interpn | 5 µs | 0 | 5566.18x | 2e-16 |
| `grid/linear/4d/max200/eval10000/flipped` | construction | SciPy | 118 ms | 22.2 µs | 5321.17x |  |
| `grid/linear/4d/max200/eval10000/ascending` | construction | SciPy | 48.4 ms | 9.69 µs | 4998.27x |  |
| `grid/linear/4d/max200/eval1000/ascending` | construction | SciPy | 50.3 ms | 10.1 µs | 4987.26x |  |
| `grid/linear/3d/max40/eval1000/ascending` | construction | interpn | 4.39 µs | 0 | 4930.63x | 2e-16 |
| `grid/linear/4d/max200/eval10/ascending` | construction | SciPy | 47 ms | 9.64 µs | 4875.94x |  |
| `grid/linear/4d/max200/eval100/ascending` | construction | SciPy | 49 ms | 10.4 µs | 4737.92x |  |
| `grid/linear/3d/max40/eval10000/ascending` | construction | interpn | 4.01 µs | 0 | 4444.41x | 3e-16 |
| `grid/linear/3d/max40/eval10/ascending` | construction | interpn | 3.7 µs | 0 | 3850.37x | 2e-16 |
| `grid/cubic/2d/each10/eval1000/ascending` | construction | interpn | 1.45 µs | 0 | 1632.44x | 3e-01 |
| `grid/nearest/2d/each40/eval10000/ascending` | construction | interpn | 694 ns | 0 | 1050.47x | 0e+00 |
| `grid/linear/2d/max40/eval10000/ascending` | construction | interpn | 409 ns | 0 | 612.70x | 4e-16 |
| `grid/linear/2d/max40/eval1000/ascending` | construction | interpn | 391 ns | 0 | 584.77x | 3e-16 |
| `grid/linear/2d/max40/eval100/ascending` | construction | interpn | 382 ns | 0 | 571.05x | 2e-16 |
| `grid/linear/2d/max40/eval10/ascending` | construction | interpn | 377 ns | 0 | 563.78x | 2e-16 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | ndarray-interp | 1.17 ms | 2.24 µs | 521.85x | 6e-16 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | ndarray-interp | 329 µs | 1.15 µs | 285.94x | 5e-16 |
| `grid/linear/3d/max200/eval10/flipped` | construction | SciPy | 4.19 ms | 16.3 µs | 256.60x |  |
| `grid/linear/3d/max200/eval100/flipped` | construction | SciPy | 3.99 ms | 16.3 µs | 244.48x |  |
| `grid/linear/3d/max200/eval10000/flipped` | construction | SciPy | 3.93 ms | 16.4 µs | 240.30x |  |
| `grid/linear/3d/max200/eval100/ascending` | construction | SciPy | 1.89 ms | 7.91 µs | 238.58x |  |
| `grid/linear/3d/max200/eval1000/flipped` | construction | SciPy | 3.88 ms | 16.5 µs | 235.62x |  |
| `grid/linear/3d/max200/eval10000/ascending` | construction | SciPy | 1.82 ms | 7.81 µs | 233.27x |  |
| `grid/linear/3d/max200/eval1000/ascending` | construction | SciPy | 1.8 ms | 7.8 µs | 231.28x |  |
| `grid/linear/3d/max200/eval10/ascending` | construction | SciPy | 1.86 ms | 8.49 µs | 219.63x |  |
| `grid/linear/3d/max10/eval100/ascending` | construction | interpn | 189 ns | 0 | 212.16x | 3e-16 |
| `grid/linear/3d/max10/eval10000/ascending` | construction | interpn | 186 ns | 0 | 208.88x | 4e-16 |
| `grid/linear/3d/max10/eval1000/ascending` | construction | interpn | 182 ns | 0 | 204.05x | 5e-16 |
| `grid/linear/3d/max10/eval10/ascending` | construction | interpn | 181 ns | 0 | 200.32x | 2e-16 |
| `grid/linear/2d/max10/eval1000/ascending` | construction | interpn | 127 ns | 0 | 189.90x | 3e-16 |
| `grid/linear/2d/max10/eval10/ascending` | construction | interpn | 125 ns | 0 | 187.10x | 1e-16 |
| `grid/linear/2d/max10/eval100/ascending` | construction | interpn | 127 ns | 0 | 184.34x | 2e-16 |
| `grid/linear/2d/max10/eval10000/ascending` | construction | interpn | 127 ns | 0 | 181.35x | 3e-16 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | SciPy | 1.17 ms | 6.56 µs | 178.42x |  |
| `grid/pchip/3d/each40/eval10000/ascending` | construction | SciPy | 507 µs | 9.65 µs | 52.56x |  |
| `grid/linear/2d/each400/eval10000/ascending` | construction | SciPy | 329 µs | 6.37 µs | 51.55x |  |
| `grid/linear/2d/max200/eval10000/ascending` | construction | ndarray-interp | 8.2 µs | 471 ns | 17.39x | 5e-16 |
| `grid/linear/2d/max200/eval10/ascending` | construction | ndarray-interp | 7.7 µs | 474 ns | 16.23x | 2e-16 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | ndarray-interp | 7.6 µs | 474 ns | 16.03x | 4e-16 |
| `grid/linear/2d/max200/eval100/ascending` | construction | ndarray-interp | 8.13 µs | 513 ns | 15.85x | 3e-16 |
| `oned/linear/n10000/shuffled` | construction | ndarray-interp | 86.4 µs | 17.4 µs | 4.97x | 0e+00 |
| `grid/linear/4d/max40/eval100/flipped` | construction | SciPy | 104 µs | 22 µs | 4.70x |  |
| `grid/linear/4d/max40/eval10/flipped` | construction | SciPy | 95.3 µs | 20.4 µs | 4.68x |  |
| `grid/linear/4d/max40/eval1000/flipped` | construction | SciPy | 94.8 µs | 20.3 µs | 4.66x |  |
| `grid/linear/4d/max40/eval10000/flipped` | construction | SciPy | 95.3 µs | 20.5 µs | 4.66x |  |
| `oned/akima/n1000/shuffled` | construction | scirs2-interpolate | 21.1 µs | 4.56 µs | 4.62x | 3e-17 |
| `oned/linear/n10000/sorted` | construction | ndarray-interp | 81.3 µs | 17.9 µs | 4.54x | 0e+00 |
| `grid/linear/4d/max40/eval10/ascending` | construction | scirs2-interpolate | 23.1 µs | 5.32 µs | 4.35x | 2e-16 |
| `grid/linear/4d/max40/eval100/ascending` | construction | scirs2-interpolate | 23.2 µs | 5.36 µs | 4.33x | 3e-16 |
| `oned/akima/n1000/sorted` | construction | scirs2-interpolate | 21.3 µs | 4.93 µs | 4.32x | 3e-17 |
| `grid/linear/4d/max40/eval1000/ascending` | construction | scirs2-interpolate | 22.9 µs | 5.5 µs | 4.16x | 2e-16 |
| `grid/linear/4d/max40/eval10000/ascending` | construction | scirs2-interpolate | 23.1 µs | 5.91 µs | 3.91x | 4e-16 |
| `grid/linear/2d/max200/eval10000/flipped` | construction | SciPy | 47.4 µs | 12.3 µs | 3.84x |  |
| `oned/akima/n100/shuffled` | construction | scirs2-interpolate | 2.08 µs | 549 ns | 3.79x | 1e-16 |
| `grid/linear/3d/max40/eval100/ascending` | construction | scirs2-interpolate | 5 µs | 1.33 µs | 3.78x | 1e-16 |
| `oned/akima/n100/sorted` | construction | scirs2-interpolate | 2.07 µs | 551 ns | 3.75x | 1e-16 |
| `oned/linear/n10000/shuffled` | construction | scirs2-interpolate | 86.4 µs | 23.8 µs | 3.63x | 0e+00 |
| `grid/linear/3d/max40/eval1000/ascending` | construction | scirs2-interpolate | 4.39 µs | 1.32 µs | 3.33x | 2e-16 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | SciPy | 26.4 µs | 8.31 µs | 3.18x |  |
| `grid/linear/3d/max40/eval10000/ascending` | construction | scirs2-interpolate | 4.01 µs | 1.32 µs | 3.04x | 2e-16 |
| `oned/akima/n50/sorted` | construction | scirs2-interpolate | 1.03 µs | 347 ns | 2.96x | 1e-16 |
| `oned/akima/n50/shuffled` | construction | scirs2-interpolate | 987 ns | 337 ns | 2.93x | 1e-16 |
| `grid/linear/3d/max40/eval10/ascending` | construction | scirs2-interpolate | 3.7 µs | 1.32 µs | 2.80x | 1e-16 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | scirs2-interpolate | 26.4 µs | 9.51 µs | 2.78x | 0e+00 |
| `grid/linear/2d/max200/eval100/ascending` | construction | scirs2-interpolate | 8.13 µs | 3 µs | 2.71x | 2e-16 |
| `grid/linear/2d/max200/eval10/flipped` | construction | SciPy | 33 µs | 12.3 µs | 2.68x |  |
| `grid/linear/2d/max200/eval10000/ascending` | construction | scirs2-interpolate | 8.2 µs | 3.06 µs | 2.68x | 2e-16 |
| `grid/linear/2d/max200/eval10/ascending` | construction | scirs2-interpolate | 7.7 µs | 3 µs | 2.56x | 1e-16 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | scirs2-interpolate | 7.6 µs | 2.98 µs | 2.55x | 2e-16 |
| `grid/linear/2d/max40/eval100/ascending` | construction | ndarray-interp | 382 ns | 151 ns | 2.53x | 2e-16 |
| `grid/linear/2d/max40/eval10/ascending` | construction | ndarray-interp | 377 ns | 151 ns | 2.50x | 2e-16 |
| `grid/linear/2d/max40/eval10000/ascending` | construction | ndarray-interp | 409 ns | 165 ns | 2.49x | 4e-16 |
| `oned/linear/n10000/sorted` | construction | scirs2-interpolate | 81.3 µs | 32.9 µs | 2.47x | 0e+00 |
| `grid/linear/2d/max200/eval100/flipped` | construction | SciPy | 32.3 µs | 13.3 µs | 2.43x |  |
| `grid/linear/2d/max200/eval1000/flipped` | construction | SciPy | 32.1 µs | 13.3 µs | 2.42x |  |
| `oned/linear/n10000/shuffled` | evaluation | ndarray-interp | 633 µs | 262 µs | 2.42x | 0e+00 |
| `grid/nearest/2d/each40/eval10000/ascending` | construction | scirs2-interpolate | 694 ns | 287 ns | 2.42x | 0e+00 |
| `grid/linear/4d/max40/eval100/ascending` | construction | SciPy | 23.2 µs | 9.67 µs | 2.40x |  |
| `grid/linear/4d/max40/eval10/ascending` | construction | SciPy | 23.1 µs | 9.65 µs | 2.40x |  |
| `grid/linear/4d/max40/eval10000/ascending` | construction | SciPy | 23.1 µs | 9.67 µs | 2.39x |  |
| `grid/linear/2d/max40/eval1000/ascending` | construction | ndarray-interp | 391 ns | 164 ns | 2.39x | 4e-16 |
| `oned/linear/n1000/shuffled` | construction | scirs2-interpolate | 2.44 µs | 1.05 µs | 2.31x | 0e+00 |
| `grid/linear/4d/max40/eval1000/ascending` | construction | SciPy | 22.9 µs | 10.4 µs | 2.20x |  |
| `oned/cubic_spline/n10000/shuffled` | evaluation | ndarray-interp | 753 µs | 345 µs | 2.18x | 2e-16 |
| `oned/pchip/n1000/sorted` | construction | scirs2-interpolate | 9.51 µs | 4.52 µs | 2.10x | 2e-16 |
| `oned/linear/n1000/sorted` | construction | scirs2-interpolate | 2.37 µs | 1.14 µs | 2.07x | 0e+00 |
| `oned/linear/n10/sorted` | construction | ndarray-interp | 79.6 ns | 39.2 ns | 2.03x | 0e+00 |
| `oned/pchip/n1000/shuffled` | construction | scirs2-interpolate | 8.89 µs | 4.6 µs | 1.93x | 2e-16 |
| `oned/akima/n10000/shuffled` | construction | scirs2-interpolate | 246 µs | 132 µs | 1.87x | 1e-16 |
| `grid/linear/2d/max40/eval10000/ascending` | construction | scirs2-interpolate | 409 ns | 224 ns | 1.83x | 2e-16 |
| `oned/linear/n10/shuffled` | construction | ndarray-interp | 78.1 ns | 43 ns | 1.81x | 0e+00 |
| `grid/linear/2d/max40/eval1000/ascending` | construction | scirs2-interpolate | 391 ns | 220 ns | 1.78x | 3e-16 |
| `oned/linear/n10000/shuffled` | construction | SciPy | 86.4 µs | 48.7 µs | 1.77x |  |
| `oned/linear/n100/shuffled` | construction | ndarray-interp | 292 ns | 165 ns | 1.77x | 0e+00 |
| `oned/linear/n1000/shuffled` | construction | ndarray-interp | 2.44 µs | 1.39 µs | 1.76x | 0e+00 |
| `grid/linear/2d/max40/eval100/ascending` | construction | scirs2-interpolate | 382 ns | 218 ns | 1.75x | 3e-16 |
| `oned/linear/n100/sorted` | construction | ndarray-interp | 291 ns | 167 ns | 1.75x | 0e+00 |
| `oned/linear/n1000/sorted` | construction | ndarray-interp | 2.37 µs | 1.38 µs | 1.72x | 0e+00 |
| `oned/linear/n50/shuffled` | construction | ndarray-interp | 171 ns | 99.7 ns | 1.71x | 0e+00 |
| `grid/linear/2d/max40/eval10/ascending` | construction | scirs2-interpolate | 377 ns | 221 ns | 1.71x | 4e-17 |
| `oned/nearest/n10000/shuffled` | construction | scirs2-interpolate | 39.6 µs | 23.4 µs | 1.69x | 0e+00 |
| `oned/pchip/n100/shuffled` | construction | scirs2-interpolate | 904 ns | 539 ns | 1.68x | 2e-16 |
| `oned/akima/n10000/sorted` | construction | scirs2-interpolate | 236 µs | 142 µs | 1.66x | 1e-16 |
| `oned/pchip/n100/sorted` | construction | scirs2-interpolate | 902 ns | 544 ns | 1.66x | 2e-16 |
| `oned/linear/n50/sorted` | construction | ndarray-interp | 159 ns | 100 ns | 1.59x | 0e+00 |
| `grid/linear/4d/max40/eval10/ascending` | evaluation | interpn | 428 ns | 269 ns | 1.59x | 2e-16 |
| `oned/linear/n100/sorted` | construction | scirs2-interpolate | 291 ns | 184 ns | 1.58x | 0e+00 |
| `grid/linear/4d/max40/eval100/ascending` | evaluation | interpn | 4.2 µs | 2.66 µs | 1.58x | 2e-16 |
| `oned/akima/n10/shuffled` | construction | scirs2-interpolate | 309 ns | 196 ns | 1.57x | 3e-17 |
| `grid/linear/4d/max40/eval1000/ascending` | evaluation | interpn | 42.1 µs | 26.8 µs | 1.57x | 6e-16 |
| `oned/akima/n10/sorted` | construction | scirs2-interpolate | 288 ns | 184 ns | 1.57x | 3e-17 |
| `grid/linear/4d/max40/eval10000/ascending` | evaluation | interpn | 421 µs | 269 µs | 1.57x | 5e-16 |
| `oned/nearest/n10000/sorted` | construction | scirs2-interpolate | 39.4 µs | 25.2 µs | 1.56x | 0e+00 |
| `grid/linear/2d/max10/eval10000/ascending` | construction | ndarray-interp | 127 ns | 84.8 ns | 1.50x | 4e-16 |
| `grid/linear/3d/max200/eval10000/ascending` | construction | scirs2-interpolate | 1.82 ms | 1.22 ms | 1.49x | 3e-16 |
| `oned/linear/n1000/shuffled` | evaluation | ndarray-interp | 37.2 µs | 25 µs | 1.49x | 0e+00 |
| `grid/linear/3d/max200/eval100/ascending` | construction | scirs2-interpolate | 1.89 ms | 1.29 ms | 1.46x | 3e-16 |
| `oned/linear/n100/shuffled` | construction | scirs2-interpolate | 292 ns | 201 ns | 1.45x | 0e+00 |
| `oned/linear/n50/shuffled` | construction | scirs2-interpolate | 171 ns | 119 ns | 1.43x | 0e+00 |
| `grid/linear/3d/max200/eval10/ascending` | construction | scirs2-interpolate | 1.86 ms | 1.31 ms | 1.43x | 3e-16 |
| `grid/linear/4d/max200/eval100/ascending` | evaluation | interpn | 5.18 µs | 3.63 µs | 1.43x | 3e-16 |
| `grid/linear/2d/max10/eval1000/ascending` | construction | ndarray-interp | 127 ns | 89.8 ns | 1.41x | 3e-16 |
| `grid/linear/2d/max10/eval100/ascending` | construction | ndarray-interp | 127 ns | 89.9 ns | 1.41x | 4e-16 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | scirs2-interpolate | 1.17 ms | 831 µs | 1.41x | 2e-16 |
| `grid/linear/3d/max200/eval1000/ascending` | construction | scirs2-interpolate | 1.8 ms | 1.28 ms | 1.40x | 3e-16 |
| `grid/linear/4d/max200/eval10/ascending` | evaluation | interpn | 504 ns | 360 ns | 1.40x | 2e-16 |
| `grid/linear/2d/max10/eval10/ascending` | construction | ndarray-interp | 125 ns | 89.8 ns | 1.39x | 1e-16 |
| `oned/pchip/n50/shuffled` | construction | scirs2-interpolate | 470 ns | 338 ns | 1.39x | 1e-16 |
| `oned/cubic_spline/n1000/shuffled` | evaluation | ndarray-interp | 44.8 µs | 33.2 µs | 1.35x | 2e-16 |
| `grid/linear/4d/max200/eval10000/ascending` | evaluation | interpn | 498 µs | 372 µs | 1.34x | 4e-16 |
| `oned/pchip/n50/sorted` | construction | scirs2-interpolate | 469 ns | 350 ns | 1.34x | 1e-16 |
| `grid/linear/2d/max200/eval100/ascending` | construction | SciPy | 8.13 µs | 6.08 µs | 1.34x |  |
| `oned/linear/n50/sorted` | construction | scirs2-interpolate | 159 ns | 120 ns | 1.33x | 0e+00 |
| `grid/linear/2d/max10/eval100/ascending` | evaluation | interpn | 1.23 µs | 929 ns | 1.32x | 2e-16 |
| `grid/linear/2d/max10/eval10/ascending` | evaluation | interpn | 123 ns | 95.2 ns | 1.30x | 1e-16 |
| `grid/linear/3d/max10/eval10000/ascending` | evaluation | interpn | 178 µs | 137 µs | 1.29x | 4e-16 |
| `grid/linear/3d/max10/eval1000/ascending` | evaluation | interpn | 17.6 µs | 13.7 µs | 1.28x | 5e-16 |
| `grid/linear/4d/max200/eval1000/ascending` | construction | scirs2-interpolate | 50.3 ms | 39.2 ms | 1.28x | 3e-16 |
| `grid/linear/3d/max10/eval10/ascending` | evaluation | interpn | 178 ns | 139 ns | 1.28x | 2e-16 |
| `grid/linear/4d/max200/eval100/ascending` | construction | scirs2-interpolate | 49 ms | 38.7 ms | 1.27x | 1e-16 |
| `grid/linear/2d/max200/eval10/ascending` | construction | SciPy | 7.7 µs | 6.09 µs | 1.26x |  |
| `grid/linear/3d/max10/eval100/ascending` | evaluation | interpn | 1.72 µs | 1.37 µs | 1.26x | 3e-16 |
| `grid/linear/4d/max200/eval1000/ascending` | evaluation | interpn | 49.8 µs | 39.6 µs | 1.26x | 3e-16 |
| `grid/linear/2d/max200/eval10000/ascending` | construction | SciPy | 8.2 µs | 6.58 µs | 1.25x |  |
| `grid/linear/4d/max200/eval10000/ascending` | construction | scirs2-interpolate | 48.4 ms | 38.8 ms | 1.25x | 4e-16 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | SciPy | 7.6 µs | 6.12 µs | 1.24x |  |
| `grid/linear/2d/max10/eval10000/ascending` | evaluation | interpn | 123 µs | 101 µs | 1.21x | 3e-16 |
| `grid/linear/2d/max10/eval1000/ascending` | evaluation | interpn | 12.3 µs | 10.1 µs | 1.21x | 3e-16 |
| `oned/cubic_spline/n1000/shuffled` | construction | scirs2-interpolate | 19.5 µs | 16.1 µs | 1.21x | 4e-04 |
| `grid/linear/2d/max200/eval1000/ascending` | evaluation | interpn | 20.4 µs | 17.3 µs | 1.18x | 2e-16 |
| `grid/linear/2d/max40/eval100/ascending` | evaluation | interpn | 1.55 µs | 1.32 µs | 1.17x | 2e-16 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | scirs2-interpolate | 329 µs | 281 µs | 1.17x | 2e-16 |
| `grid/linear/2d/max40/eval1000/ascending` | evaluation | interpn | 15.3 µs | 13.1 µs | 1.17x | 3e-16 |
| `grid/linear/2d/max40/eval10000/ascending` | evaluation | interpn | 152 µs | 131 µs | 1.16x | 4e-16 |
| `oned/cubic_spline/n1000/sorted` | construction | scirs2-interpolate | 20.2 µs | 17.4 µs | 1.16x | 4e-04 |
| `grid/linear/4d/max200/eval10/ascending` | construction | scirs2-interpolate | 47 ms | 40.7 ms | 1.16x | 2e-16 |
| `grid/linear/3d/max40/eval10/ascending` | evaluation | interpn | 228 ns | 198 ns | 1.15x | 2e-16 |
| `grid/linear/3d/max200/eval100/ascending` | evaluation | interpn | 3 µs | 2.63 µs | 1.14x | 3e-16 |
| `grid/linear/3d/max40/eval1000/ascending` | evaluation | interpn | 22.3 µs | 19.8 µs | 1.13x | 2e-16 |
| `grid/linear/3d/max40/eval10000/ascending` | evaluation | interpn | 222 µs | 197 µs | 1.13x | 3e-16 |
| `grid/linear/2d/max40/eval10/ascending` | evaluation | interpn | 154 ns | 136 ns | 1.13x | 2e-16 |
| `grid/linear/3d/max40/eval100/ascending` | evaluation | interpn | 2.23 µs | 1.98 µs | 1.12x | 2e-16 |
| `grid/linear/3d/max10/eval100/ascending` | construction | scirs2-interpolate | 189 ns | 168 ns | 1.12x | 2e-16 |
| `rbf/quintic/global/n10/scattered` | evaluation | SciPy | 75 µs | 66.8 µs | 1.12x |  |

</details>

## Runs that failed

| case | library | how |
| :--- | :--- | :--- |
| `grid/linear/4d/max10/eval10/ascending` | interpn | exit code 101 |
| `grid/linear/4d/max10/eval10/ascending` | interpx | exit code 101 |
| `grid/linear/4d/max10/eval10/ascending` | scirs2-interpolate | exit code 101 |
| `grid/linear/4d/max10/eval10/flipped` | interpx | exit code 101 |
| `grid/linear/4d/max10/eval100/ascending` | interpn | exit code 101 |
| `grid/linear/4d/max10/eval100/ascending` | interpx | exit code 101 |
| `grid/linear/4d/max10/eval100/ascending` | scirs2-interpolate | exit code 101 |
| `grid/linear/4d/max10/eval100/flipped` | interpx | exit code 101 |
| `grid/linear/4d/max10/eval1000/ascending` | interpn | exit code 101 |
| `grid/linear/4d/max10/eval1000/ascending` | interpx | exit code 101 |
| `grid/linear/4d/max10/eval1000/ascending` | scirs2-interpolate | exit code 101 |
| `grid/linear/4d/max10/eval1000/flipped` | interpx | exit code 101 |
| `grid/linear/4d/max10/eval10000/ascending` | interpn | exit code 101 |
| `grid/linear/4d/max10/eval10000/ascending` | interpx | exit code 101 |
| `grid/linear/4d/max10/eval10000/ascending` | scirs2-interpolate | exit code 101 |
| `grid/linear/4d/max10/eval10000/flipped` | interpx | exit code 101 |

## One variable

Data of `n` points, evaluated at `4 n` points: in order (`sorted`) and in no order (`shuffled`).

### Evaluation

| case | interpx | SciPy | scirs2-interpolate | ndarray-interp | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: |
| `linear/n10/sorted` | 127 ns | 5.4 µs | 304 ns | 255 ns | 2x faster |
| `linear/n10/shuffled` | 184 ns | 5.41 µs | 300 ns | 272 ns | 1.5x faster |
| `nearest/n10/sorted` | 162 ns | 10.2 µs | 272 ns |  | 1.7x faster |
| `nearest/n10/shuffled` | 204 ns | 10.4 µs | 280 ns |  | 1.4x faster |
| `previous/n10/sorted` | 157 ns | 11.2 µs |  |  | 71x faster |
| `previous/n10/shuffled` | 195 ns | 11.3 µs |  |  | 58x faster |
| `bspline1/n10/sorted` | 303 ns | 2.76 µs |  |  | 9.1x faster |
| `bspline1/n10/shuffled` | 343 ns | 2.9 µs |  |  | 8.5x faster |
| `bspline2/n10/sorted` | 407 ns | 3.16 µs |  |  | 7.8x faster |
| `bspline2/n10/shuffled` | 474 ns | 3.3 µs |  |  | 7x faster |
| `bspline3/n10/sorted` | 538 ns | 3.7 µs |  |  | 6.9x faster |
| `bspline3/n10/shuffled` | 613 ns | 3.71 µs |  |  | 6.1x faster |
| `cubic_spline/n10/sorted` | 135 ns | 3.81 µs | 276 ns | 340 ns | 2x faster |
| `cubic_spline/n10/shuffled` | 227 ns | 3.94 µs | 280 ns | 365 ns | 1.2x faster |
| `pchip/n10/sorted` | 138 ns | 3.68 µs | 292 ns |  | 2.1x faster |
| `pchip/n10/shuffled` | 201 ns | 3.8 µs | 292 ns |  | 1.5x faster |
| `akima/n10/sorted` | 131 ns | 3.78 µs | 251 ns |  | 1.9x faster |
| `akima/n10/shuffled` | 210 ns | 3.8 µs | 285 ns |  | 1.4x faster |
| `makima/n10/sorted` | 138 ns | 3.76 µs |  |  | 27x faster |
| `makima/n10/shuffled` | 209 ns | 3.83 µs |  |  | 18x faster |
| `linear/n50/sorted` | 641 ns | 6.49 µs | 1.66 µs | 1.35 µs | 2.1x faster |
| `linear/n50/shuffled` | 1.2 µs | 8.46 µs | 1.76 µs | 1.28 µs | 1.1x faster |
| `nearest/n50/sorted` | 830 ns | 13 µs | 1.83 µs |  | 2.2x faster |
| `nearest/n50/shuffled` | 1.32 µs | 12.5 µs | 1.89 µs |  | 1.4x faster |
| `previous/n50/sorted` | 821 ns | 14.6 µs |  |  | 18x faster |
| `previous/n50/shuffled` | 1.28 µs | 13.3 µs |  |  | 10x faster |
| `bspline1/n50/sorted` | 1.46 µs | 4.06 µs |  |  | 2.8x faster |
| `bspline1/n50/shuffled` | 1.96 µs | 5.03 µs |  |  | 2.6x faster |
| `bspline2/n50/sorted` | 1.97 µs | 6.37 µs |  |  | 3.2x faster |
| `bspline2/n50/shuffled` | 2.77 µs | 6.98 µs |  |  | 2.5x faster |
| `bspline3/n50/sorted` | 2.66 µs | 8.98 µs |  |  | 3.4x faster |
| `bspline3/n50/shuffled` | 3.75 µs | 9.49 µs |  |  | 2.5x faster |
| `cubic_spline/n50/sorted` | 693 ns | 5.32 µs | 4.23 µs | 1.64 µs | 2.4x faster |
| `cubic_spline/n50/shuffled` | 1.34 µs | 5.98 µs | 4.02 µs | 1.63 µs | 1.2x faster |
| `pchip/n50/sorted` | 751 ns | 5.25 µs | 3.92 µs |  | 5.2x faster |
| `pchip/n50/shuffled` | 1.34 µs | 5.6 µs | 3.65 µs |  | 2.7x faster |
| `akima/n50/sorted` | 695 ns | 5.45 µs | 4.42 µs |  | 6.4x faster |
| `akima/n50/shuffled` | 1.36 µs | 5.6 µs | 3.89 µs |  | 2.9x faster |
| `makima/n50/sorted` | 710 ns | 5.25 µs |  |  | 7.4x faster |
| `makima/n50/shuffled` | 1.33 µs | 5.69 µs |  |  | 4.3x faster |
| `linear/n100/sorted` | 1.31 µs | 7.93 µs | 3.51 µs | 2.54 µs | 1.9x faster |
| `linear/n100/shuffled` | 2.69 µs | 13 µs | 3.81 µs | 2.49 µs | **1.1x slower** |
| `nearest/n100/sorted` | 1.48 µs | 14.4 µs | 3.61 µs |  | 2.4x faster |
| `nearest/n100/shuffled` | 2.93 µs | 14.3 µs | 3.77 µs |  | 1.3x faster |
| `previous/n100/sorted` | 1.49 µs | 16 µs |  |  | 11x faster |
| `previous/n100/shuffled` | 2.87 µs | 16 µs |  |  | 5.6x faster |
| `bspline1/n100/sorted` | 2.97 µs | 5.88 µs |  |  | 2x faster |
| `bspline1/n100/shuffled` | 4.33 µs | 12.8 µs |  |  | 3x faster |
| `bspline2/n100/sorted` | 3.95 µs | 10.1 µs |  |  | 2.6x faster |
| `bspline2/n100/shuffled` | 6.11 µs | 12.7 µs |  |  | 2.1x faster |
| `bspline3/n100/sorted` | 5.24 µs | 15.1 µs |  |  | 2.9x faster |
| `bspline3/n100/shuffled` | 8.02 µs | 18 µs |  |  | 2.2x faster |
| `cubic_spline/n100/sorted` | 1.48 µs | 7.1 µs | 17.1 µs | 3.27 µs | 2.2x faster |
| `cubic_spline/n100/shuffled` | 3.04 µs | 8.05 µs | 17.5 µs | 3.27 µs | 1.1x faster |
| `pchip/n100/sorted` | 1.46 µs | 6.99 µs | 17.1 µs |  | 4.8x faster |
| `pchip/n100/shuffled` | 3.05 µs | 8.57 µs | 15.8 µs |  | 2.8x faster |
| `akima/n100/sorted` | 1.46 µs | 6.97 µs | 16.7 µs |  | 4.8x faster |
| `akima/n100/shuffled` | 3.11 µs | 8.08 µs | 17.6 µs |  | 2.6x faster |
| `makima/n100/sorted` | 1.46 µs | 6.94 µs |  |  | 4.8x faster |
| `makima/n100/shuffled` | 3.35 µs | 8.71 µs |  |  | 2.6x faster |
| `linear/n1000/sorted` | 13.3 µs | 31.4 µs | 108 µs | 24.9 µs | 1.9x faster |
| `linear/n1000/shuffled` | 37.2 µs | 143 µs | 213 µs | 25 µs | **1.5x slower** |
| `nearest/n1000/sorted` | 14.7 µs | 57 µs | 141 µs |  | 3.9x faster |
| `nearest/n1000/shuffled` | 42.2 µs | 205 µs | 230 µs |  | 4.8x faster |
| `previous/n1000/sorted` | 14.7 µs | 61.1 µs |  |  | 4.2x faster |
| `previous/n1000/shuffled` | 40.6 µs | 222 µs |  |  | 5.5x faster |
| `bspline1/n1000/sorted` | 31 µs | 31.4 µs |  |  | 1x faster |
| `bspline1/n1000/shuffled` | 55.3 µs | 553 µs |  |  | 10x faster |
| `bspline2/n1000/sorted` | 41.7 µs | 79.3 µs |  |  | 1.9x faster |
| `bspline2/n1000/shuffled` | 73.2 µs | 549 µs |  |  | 7.5x faster |
| `bspline3/n1000/sorted` | 53.7 µs | 130 µs |  |  | 2.4x faster |
| `bspline3/n1000/shuffled` | 97.1 µs | 631 µs |  |  | 6.5x faster |
| `cubic_spline/n1000/sorted` | 15.3 µs | 40.9 µs | 1.53 ms | 35.3 µs | 2.3x faster |
| `cubic_spline/n1000/shuffled` | 44.8 µs | 221 µs | 1.55 ms | 33.2 µs | **1.3x slower** |
| `pchip/n1000/sorted` | 14.7 µs | 37.8 µs | 1.31 ms |  | 2.6x faster |
| `pchip/n1000/shuffled` | 44.4 µs | 237 µs | 1.32 ms |  | 5.4x faster |
| `akima/n1000/sorted` | 14.8 µs | 38.5 µs | 1.51 ms |  | 2.6x faster |
| `akima/n1000/shuffled` | 48.5 µs | 235 µs | 1.55 ms |  | 4.8x faster |
| `makima/n1000/sorted` | 14.6 µs | 40.2 µs |  |  | 2.8x faster |
| `makima/n1000/shuffled` | 44.9 µs | 224 µs |  |  | 5x faster |
| `linear/n10000/sorted` | 140 µs | 241 µs | 1.61 ms | 248 µs | 1.7x faster |
| `linear/n10000/shuffled` | 633 µs | 2.51 ms | 3.09 ms | 262 µs | **2.4x slower** |
| `nearest/n10000/sorted` | 160 µs | 626 µs | 2 ms |  | 3.9x faster |
| `nearest/n10000/shuffled` | 726 µs | 3.18 ms | 3.29 ms |  | 4.4x faster |
| `previous/n10000/sorted` | 166 µs | 661 µs |  |  | 4x faster |
| `previous/n10000/shuffled` | 709 µs | 3.28 ms |  |  | 4.6x faster |
| `bspline1/n10000/sorted` | 319 µs | 307 µs |  |  | **1x slower** |
| `bspline1/n10000/shuffled` | 850 µs | 53.3 ms |  |  | 63x faster |
| `bspline2/n10000/sorted` | 458 µs | 815 µs |  |  | 1.8x faster |
| `bspline2/n10000/shuffled` | 1.03 ms | 52.7 ms |  |  | 51x faster |
| `bspline3/n10000/sorted` | 549 µs | 1.38 ms |  |  | 2.5x faster |
| `bspline3/n10000/shuffled` | 1.39 ms | 52.5 ms |  |  | 38x faster |
| `cubic_spline/n10000/sorted` | 153 µs | 390 µs | 152 ms | 331 µs | 2.2x faster |
| `cubic_spline/n10000/shuffled` | 753 µs | 3.49 ms | 153 ms | 345 µs | **2.2x slower** |
| `pchip/n10000/sorted` | 154 µs | 387 µs | 130 ms |  | 2.5x faster |
| `pchip/n10000/shuffled` | 722 µs | 3.83 ms | 128 ms |  | 5.3x faster |
| `akima/n10000/sorted` | 154 µs | 397 µs | 152 ms |  | 2.6x faster |
| `akima/n10000/shuffled` | 761 µs | 3.64 ms | 152 ms |  | 4.8x faster |
| `makima/n10000/sorted` | 152 µs | 382 µs |  |  | 2.5x faster |
| `makima/n10000/shuffled` | 739 µs | 3.64 ms |  |  | 4.9x faster |

### Construction

| case | interpx | SciPy | scirs2-interpolate | ndarray-interp | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: |
| `linear/n10/sorted` | 79.6 ns | 9.47 µs | 90.9 ns | 39.2 ns | **2x slower** |
| `linear/n10/shuffled` | 78.1 ns | 9.3 µs | 87.5 ns | 43 ns | **1.8x slower** |
| `nearest/n10/sorted` | 64.5 ns | 10.9 µs | 79.4 ns |  | 1.2x faster |
| `nearest/n10/shuffled` | 61.6 ns | 11 µs | 82.2 ns |  | 1.3x faster |
| `previous/n10/sorted` | 60 ns | 10.7 µs |  |  | 1.8e+02x faster |
| `previous/n10/shuffled` | 59.5 ns | 9.83 µs |  |  | 1.7e+02x faster |
| `bspline1/n10/sorted` | 107 ns | 16.2 µs |  |  | 1.5e+02x faster |
| `bspline1/n10/shuffled` | 101 ns | 16.4 µs |  |  | 1.6e+02x faster |
| `bspline2/n10/sorted` | 462 ns | 38.2 µs |  |  | 83x faster |
| `bspline2/n10/shuffled` | 451 ns | 41.2 µs |  |  | 91x faster |
| `bspline3/n10/sorted` | 518 ns | 37.9 µs |  |  | 73x faster |
| `bspline3/n10/shuffled` | 504 ns | 37.6 µs |  |  | 75x faster |
| `cubic_spline/n10/sorted` | 239 ns | 89.9 µs | 417 ns | 477 ns | 1.7x faster |
| `cubic_spline/n10/shuffled` | 220 ns | 95.3 µs | 380 ns | 519 ns | 1.7x faster |
| `pchip/n10/sorted` | 163 ns | 111 µs | 195 ns |  | 1.2x faster |
| `pchip/n10/shuffled` | 164 ns | 105 µs | 193 ns |  | 1.2x faster |
| `akima/n10/sorted` | 288 ns | 78.2 µs | 184 ns |  | **1.6x slower** |
| `akima/n10/shuffled` | 309 ns | 84.7 µs | 196 ns |  | **1.6x slower** |
| `makima/n10/sorted` | 290 ns | 82.1 µs |  |  | 2.8e+02x faster |
| `makima/n10/shuffled` | 290 ns | 80.7 µs |  |  | 2.8e+02x faster |
| `linear/n50/sorted` | 159 ns | 9.47 µs | 120 ns | 100 ns | **1.6x slower** |
| `linear/n50/shuffled` | 171 ns | 10.1 µs | 119 ns | 99.7 ns | **1.7x slower** |
| `nearest/n50/sorted` | 98.7 ns | 10.2 µs | 126 ns |  | 1.3x faster |
| `nearest/n50/shuffled` | 97.2 ns | 11 µs | 119 ns |  | 1.2x faster |
| `previous/n50/sorted` | 97.3 ns | 10.1 µs |  |  | 1e+02x faster |
| `previous/n50/shuffled` | 99.3 ns | 10.1 µs |  |  | 1e+02x faster |
| `bspline1/n50/sorted` | 191 ns | 16.5 µs |  |  | 86x faster |
| `bspline1/n50/shuffled` | 190 ns | 16.7 µs |  |  | 88x faster |
| `bspline2/n50/sorted` | 1.81 µs | 41.8 µs |  |  | 23x faster |
| `bspline2/n50/shuffled` | 1.82 µs | 41.8 µs |  |  | 23x faster |
| `bspline3/n50/sorted` | 2.01 µs | 40.8 µs |  |  | 20x faster |
| `bspline3/n50/shuffled` | 2.04 µs | 41.5 µs |  |  | 20x faster |
| `cubic_spline/n50/sorted` | 822 ns | 91.5 µs | 962 ns | 1.08 µs | 1.2x faster |
| `cubic_spline/n50/shuffled` | 844 ns | 92.8 µs | 925 ns | 1.08 µs | 1.1x faster |
| `pchip/n50/sorted` | 469 ns | 106 µs | 350 ns |  | **1.3x slower** |
| `pchip/n50/shuffled` | 470 ns | 107 µs | 338 ns |  | **1.4x slower** |
| `akima/n50/sorted` | 1.03 µs | 84 µs | 347 ns |  | **3x slower** |
| `akima/n50/shuffled` | 987 ns | 79.4 µs | 337 ns |  | **2.9x slower** |
| `makima/n50/sorted` | 996 ns | 81.4 µs |  |  | 82x faster |
| `makima/n50/shuffled` | 1.02 µs | 81.5 µs |  |  | 80x faster |
| `linear/n100/sorted` | 291 ns | 10.1 µs | 184 ns | 167 ns | **1.7x slower** |
| `linear/n100/shuffled` | 292 ns | 9.54 µs | 201 ns | 165 ns | **1.8x slower** |
| `nearest/n100/sorted` | 182 ns | 10.5 µs | 186 ns |  | 1x faster |
| `nearest/n100/shuffled` | 179 ns | 10.5 µs | 186 ns |  | 1x faster |
| `previous/n100/sorted` | 179 ns | 10.5 µs |  |  | 58x faster |
| `previous/n100/shuffled` | 179 ns | 10.6 µs |  |  | 59x faster |
| `bspline1/n100/sorted` | 335 ns | 16.6 µs |  |  | 50x faster |
| `bspline1/n100/shuffled` | 332 ns | 16.5 µs |  |  | 50x faster |
| `bspline2/n100/sorted` | 3.81 µs | 48.3 µs |  |  | 13x faster |
| `bspline2/n100/shuffled` | 3.54 µs | 47.2 µs |  |  | 13x faster |
| `bspline3/n100/sorted` | 3.91 µs | 49 µs |  |  | 13x faster |
| `bspline3/n100/shuffled` | 3.93 µs | 48.6 µs |  |  | 12x faster |
| `cubic_spline/n100/sorted` | 1.66 µs | 91.5 µs | 1.64 µs | 2.06 µs | **1x slower** |
| `cubic_spline/n100/shuffled` | 1.75 µs | 91.5 µs | 1.66 µs | 2.05 µs | **1.1x slower** |
| `pchip/n100/sorted` | 902 ns | 105 µs | 544 ns |  | **1.7x slower** |
| `pchip/n100/shuffled` | 904 ns | 108 µs | 539 ns |  | **1.7x slower** |
| `akima/n100/sorted` | 2.07 µs | 79.2 µs | 551 ns |  | **3.8x slower** |
| `akima/n100/shuffled` | 2.08 µs | 79.9 µs | 549 ns |  | **3.8x slower** |
| `makima/n100/sorted` | 2.08 µs | 84.1 µs |  |  | 40x faster |
| `makima/n100/shuffled` | 1.91 µs | 89.9 µs |  |  | 47x faster |
| `linear/n1000/sorted` | 2.37 µs | 14 µs | 1.14 µs | 1.38 µs | **2.1x slower** |
| `linear/n1000/shuffled` | 2.44 µs | 13.8 µs | 1.05 µs | 1.39 µs | **2.3x slower** |
| `nearest/n1000/sorted` | 1.13 µs | 15.4 µs | 1.05 µs |  | **1.1x slower** |
| `nearest/n1000/shuffled` | 1.05 µs | 14.9 µs | 1.05 µs |  | **1x slower** |
| `previous/n1000/sorted` | 1.05 µs | 17 µs |  |  | 16x faster |
| `previous/n1000/shuffled` | 1.06 µs | 17.5 µs |  |  | 16x faster |
| `bspline1/n1000/sorted` | 1.92 µs | 18.9 µs |  |  | 9.8x faster |
| `bspline1/n1000/shuffled` | 1.9 µs | 17.7 µs |  |  | 9.3x faster |
| `bspline2/n1000/sorted` | 33.7 µs | 114 µs |  |  | 3.4x faster |
| `bspline2/n1000/shuffled` | 34.6 µs | 120 µs |  |  | 3.5x faster |
| `bspline3/n1000/sorted` | 44.2 µs | 132 µs |  |  | 3x faster |
| `bspline3/n1000/shuffled` | 45.4 µs | 139 µs |  |  | 3.1x faster |
| `cubic_spline/n1000/sorted` | 20.2 µs | 120 µs | 17.4 µs | 15.8 µs | **1.3x slower** |
| `cubic_spline/n1000/shuffled` | 19.5 µs | 124 µs | 16.1 µs | 17.2 µs | **1.2x slower** |
| `pchip/n1000/sorted` | 9.51 µs | 123 µs | 4.52 µs |  | **2.1x slower** |
| `pchip/n1000/shuffled` | 8.89 µs | 126 µs | 4.6 µs |  | **1.9x slower** |
| `akima/n1000/sorted` | 21.3 µs | 105 µs | 4.93 µs |  | **4.3x slower** |
| `akima/n1000/shuffled` | 21.1 µs | 104 µs | 4.56 µs |  | **4.6x slower** |
| `makima/n1000/sorted` | 20.7 µs | 108 µs |  |  | 5.2x faster |
| `makima/n1000/shuffled` | 21.6 µs | 103 µs |  |  | 4.8x faster |
| `linear/n10000/sorted` | 81.3 µs | 47 µs | 32.9 µs | 17.9 µs | **4.5x slower** |
| `linear/n10000/shuffled` | 86.4 µs | 48.7 µs | 23.8 µs | 17.4 µs | **5x slower** |
| `nearest/n10000/sorted` | 39.4 µs | 54.8 µs | 25.2 µs |  | **1.6x slower** |
| `nearest/n10000/shuffled` | 39.6 µs | 55.3 µs | 23.4 µs |  | **1.7x slower** |
| `previous/n10000/sorted` | 39.5 µs | 80.1 µs |  |  | 2x faster |
| `previous/n10000/shuffled` | 46.6 µs | 79.1 µs |  |  | 1.7x faster |
| `bspline1/n10000/sorted` | 71.1 µs | 25.2 µs |  |  | **2.8x slower** |
| `bspline1/n10000/shuffled` | 74.2 µs | 25 µs |  |  | **3x slower** |
| `bspline2/n10000/sorted` | 522 µs | 793 µs |  |  | 1.5x faster |
| `bspline2/n10000/shuffled` | 513 µs | 820 µs |  |  | 1.6x faster |
| `bspline3/n10000/sorted` | 619 µs | 1.01 ms |  |  | 1.6x faster |
| `bspline3/n10000/shuffled` | 633 µs | 977 µs |  |  | 1.5x faster |
| `cubic_spline/n10000/sorted` | 232 µs | 778 µs | 313 µs | 363 µs | 1.4x faster |
| `cubic_spline/n10000/shuffled` | 246 µs | 761 µs | 298 µs | 360 µs | 1.2x faster |
| `pchip/n10000/sorted` | 125 µs | 246 µs | 169 µs |  | 1.4x faster |
| `pchip/n10000/shuffled` | 124 µs | 271 µs | 163 µs |  | 1.3x faster |
| `akima/n10000/sorted` | 236 µs | 678 µs | 142 µs |  | **1.7x slower** |
| `akima/n10000/shuffled` | 246 µs | 663 µs | 132 µs |  | **1.9x slower** |
| `makima/n10000/sorted` | 237 µs | 744 µs |  |  | 3.1x faster |
| `makima/n10000/shuffled` | 249 µs | 766 µs |  |  | 3.1x faster |

## `PPoly`

A cubic of 55 intervals, evaluated at points in no order.

### Evaluation

| case | interpx | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: |
| `intervals55/eval10` | 73.4 ns | 3.83 µs | 52x faster |
| `intervals55/eval100` | 676 ns | 4.46 µs | 6.6x faster |
| `intervals55/eval1000` | 6.72 µs | 16.8 µs | 2.5x faster |
| `intervals55/eval10000` | 68.4 µs | 354 µs | 5.2x faster |

### Construction

| case | interpx | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: |
| `intervals55/eval10` | 173 ns | 5.62 µs | 32x faster |
| `intervals55/eval100` | 172 ns | 6 µs | 35x faster |
| `intervals55/eval1000` | 184 ns | 5.96 µs | 32x faster |
| `intervals55/eval10000` | 182 ns | 5.52 µs | 30x faster |

## `BSpline`

A cubic B-spline of 55 intervals, evaluated at points in no order.

### Evaluation

| case | interpx | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: |
| `intervals55/eval10` | 196 ns | 2.99 µs | 15x faster |
| `intervals55/eval100` | 1.9 µs | 5.96 µs | 3.1x faster |
| `intervals55/eval1000` | 18.9 µs | 38.6 µs | 2x faster |
| `intervals55/eval10000` | 191 µs | 393 µs | 2.1x faster |

### Construction

| case | interpx | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: |
| `intervals55/eval10` | 135 ns | 10.6 µs | 78x faster |
| `intervals55/eval100` | 131 ns | 9.66 µs | 74x faster |
| `intervals55/eval1000` | 138 ns | 9.81 µs | 71x faster |
| `intervals55/eval10000` | 133 ns | 10 µs | 75x faster |

## `BSplineCurve`

A cubic curve through `n` points of the plane, evaluated at `4 n` parameters in order.

### Evaluation

| case | interpx | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: |
| `degree3/n100` | 5.78 µs | 17.3 µs | 3x faster |
| `degree3/n10000` | 622 µs | 1.45 ms | 2.3x faster |

### Construction

| case | interpx | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: |
| `degree3/n100` | 4.39 µs | 48.5 µs | 11x faster |
| `degree3/n10000` | 581 µs | 1.21 ms | 2.1x faster |

## Regular grid

Points in no order inside the grid. `maxN`: the first axis has `N` points and each axis half the one before (SciPy's suite); `eachN`: every axis has `N` points; `flipped`: every axis given from its high end.

### Evaluation

| case | interpx | SciPy | scirs2-interpolate | ndarray-interp | interpn | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/2d/max10/eval10/ascending` | 123 ns | 22.7 µs | 783 ns | 165 ns | 95.2 ns | **1.3x slower** |
| `linear/2d/max10/eval10/flipped` | 123 ns | 22.7 µs |  |  |  | 1.9e+02x faster |
| `linear/2d/max10/eval100/ascending` | 1.23 µs | 25.1 µs | 8.2 µs | 1.37 µs | 929 ns | **1.3x slower** |
| `linear/2d/max10/eval100/flipped` | 1.22 µs | 25.1 µs |  |  |  | 21x faster |
| `linear/2d/max10/eval1000/ascending` | 12.3 µs | 45.3 µs | 90.3 µs | 15.1 µs | 10.1 µs | **1.2x slower** |
| `linear/2d/max10/eval1000/flipped` | 13.3 µs | 43.7 µs |  |  |  | 3.3x faster |
| `linear/2d/max10/eval10000/ascending` | 123 µs | 369 µs | 864 µs | 220 µs | 101 µs | **1.2x slower** |
| `linear/2d/max10/eval10000/flipped` | 134 µs | 385 µs |  |  |  | 2.9x faster |
| `linear/2d/max40/eval10/ascending` | 154 ns | 21.2 µs | 992 ns | 157 ns | 136 ns | **1.1x slower** |
| `linear/2d/max40/eval10/flipped` | 155 ns | 21.1 µs |  |  |  | 1.4e+02x faster |
| `linear/2d/max40/eval100/ascending` | 1.55 µs | 23.8 µs | 10.6 µs | 1.45 µs | 1.32 µs | **1.2x slower** |
| `linear/2d/max40/eval100/flipped` | 1.52 µs | 23.6 µs |  |  |  | 16x faster |
| `linear/2d/max40/eval1000/ascending` | 15.3 µs | 48.5 µs | 99.3 µs | 15.2 µs | 13.1 µs | **1.2x slower** |
| `linear/2d/max40/eval1000/flipped` | 15.1 µs | 48.6 µs |  |  |  | 3.2x faster |
| `linear/2d/max40/eval10000/ascending` | 152 µs | 611 µs | 1.05 ms | 226 µs | 131 µs | **1.2x slower** |
| `linear/2d/max40/eval10000/flipped` | 153 µs | 610 µs |  |  |  | 4x faster |
| `linear/2d/max200/eval10/ascending` | 192 ns | 22.3 µs | 1.95 µs | 189 ns | 175 ns | **1.1x slower** |
| `linear/2d/max200/eval10/flipped` | 192 ns | 21.7 µs |  |  |  | 1.1e+02x faster |
| `linear/2d/max200/eval100/ascending` | 2.04 µs | 26.7 µs | 19.3 µs | 1.99 µs | 1.89 µs | **1.1x slower** |
| `linear/2d/max200/eval100/flipped` | 1.88 µs | 24.3 µs |  |  |  | 13x faster |
| `linear/2d/max200/eval1000/ascending` | 20.4 µs | 75 µs | 207 µs | 23.5 µs | 17.3 µs | **1.2x slower** |
| `linear/2d/max200/eval1000/flipped` | 18.7 µs | 75.6 µs |  |  |  | 4x faster |
| `linear/2d/max200/eval10000/ascending` | 186 µs | 876 µs | 2.28 ms | 452 µs | 172 µs | **1.1x slower** |
| `linear/2d/max200/eval10000/flipped` | 187 µs | 856 µs |  |  |  | 4.6x faster |
| `linear/3d/max10/eval10/ascending` | 178 ns | 61.6 µs | 3.17 µs |  | 139 ns | **1.3x slower** |
| `linear/3d/max10/eval10/flipped` | 173 ns | 59.6 µs |  |  |  | 3.4e+02x faster |
| `linear/3d/max10/eval100/ascending` | 1.72 µs | 65.5 µs | 31.2 µs |  | 1.37 µs | **1.3x slower** |
| `linear/3d/max10/eval100/flipped` | 1.73 µs | 65.9 µs |  |  |  | 38x faster |
| `linear/3d/max10/eval1000/ascending` | 17.6 µs | 121 µs | 332 µs |  | 13.7 µs | **1.3x slower** |
| `linear/3d/max10/eval1000/flipped` | 17.2 µs | 128 µs |  |  |  | 7.5x faster |
| `linear/3d/max10/eval10000/ascending` | 178 µs | 1.13 ms | 3.4 ms |  | 137 µs | **1.3x slower** |
| `linear/3d/max10/eval10000/flipped` | 176 µs | 1.13 ms |  |  |  | 6.4x faster |
| `linear/3d/max40/eval10/ascending` | 228 ns | 59.3 µs | 3.32 µs |  | 198 ns | **1.2x slower** |
| `linear/3d/max40/eval10/flipped` | 227 ns | 60.1 µs |  |  |  | 2.6e+02x faster |
| `linear/3d/max40/eval100/ascending` | 2.23 µs | 66.2 µs | 33.6 µs |  | 1.98 µs | **1.1x slower** |
| `linear/3d/max40/eval100/flipped` | 2.23 µs | 66.5 µs |  |  |  | 30x faster |
| `linear/3d/max40/eval1000/ascending` | 22.3 µs | 149 µs | 342 µs |  | 19.8 µs | **1.1x slower** |
| `linear/3d/max40/eval1000/flipped` | 22.4 µs | 142 µs |  |  |  | 6.3x faster |
| `linear/3d/max40/eval10000/ascending` | 222 µs | 1.43 ms | 3.57 ms |  | 197 µs | **1.1x slower** |
| `linear/3d/max40/eval10000/flipped` | 222 µs | 1.51 ms |  |  |  | 6.8x faster |
| `linear/3d/max200/eval10/ascending` | 281 ns | 61.4 µs | 4.53 µs |  | 265 ns | **1.1x slower** |
| `linear/3d/max200/eval10/flipped` | 302 ns | 60.3 µs |  |  |  | 2e+02x faster |
| `linear/3d/max200/eval100/ascending` | 3 µs | 69.3 µs | 46.8 µs |  | 2.63 µs | **1.1x slower** |
| `linear/3d/max200/eval100/flipped` | 2.75 µs | 69.9 µs |  |  |  | 25x faster |
| `linear/3d/max200/eval1000/ascending` | 27.6 µs | 201 µs | 464 µs |  | 26.2 µs | **1.1x slower** |
| `linear/3d/max200/eval1000/flipped` | 27.6 µs | 188 µs |  |  |  | 6.8x faster |
| `linear/3d/max200/eval10000/ascending` | 277 µs | 1.94 ms | 5.04 ms |  | 262 µs | **1.1x slower** |
| `linear/3d/max200/eval10000/flipped` | 276 µs | 1.87 ms |  |  |  | 6.8x faster |
| `linear/4d/max10/eval10/ascending` |  | 99.5 µs |  |  |  |  |
| `linear/4d/max10/eval10/flipped` |  | 106 µs |  |  |  |  |
| `linear/4d/max10/eval100/ascending` |  | 111 µs |  |  |  |  |
| `linear/4d/max10/eval100/flipped` |  | 113 µs |  |  |  |  |
| `linear/4d/max10/eval1000/ascending` |  | 210 µs |  |  |  |  |
| `linear/4d/max10/eval1000/flipped` |  | 210 µs |  |  |  |  |
| `linear/4d/max10/eval10000/ascending` |  | 1.82 ms |  |  |  |  |
| `linear/4d/max10/eval10000/flipped` |  | 1.92 ms |  |  |  |  |
| `linear/4d/max40/eval10/ascending` | 428 ns | 103 µs | 5.8 µs |  | 269 ns | **1.6x slower** |
| `linear/4d/max40/eval10/flipped` | 427 ns | 100 µs |  |  |  | 2.3e+02x faster |
| `linear/4d/max40/eval100/ascending` | 4.2 µs | 111 µs | 57.6 µs |  | 2.66 µs | **1.6x slower** |
| `linear/4d/max40/eval100/flipped` | 4.6 µs | 111 µs |  |  |  | 24x faster |
| `linear/4d/max40/eval1000/ascending` | 42.1 µs | 256 µs | 580 µs |  | 26.8 µs | **1.6x slower** |
| `linear/4d/max40/eval1000/flipped` | 42.6 µs | 254 µs |  |  |  | 6x faster |
| `linear/4d/max40/eval10000/ascending` | 421 µs | 2.33 ms | 6.31 ms |  | 269 µs | **1.6x slower** |
| `linear/4d/max40/eval10000/flipped` | 428 µs | 2.38 ms |  |  |  | 5.5x faster |
| `linear/4d/max200/eval10/ascending` | 504 ns | 99.8 µs | 7.8 µs |  | 360 ns | **1.4x slower** |
| `linear/4d/max200/eval10/flipped` | 516 ns | 109 µs |  |  |  | 2.1e+02x faster |
| `linear/4d/max200/eval100/ascending` | 5.18 µs | 123 µs | 73.9 µs |  | 3.63 µs | **1.4x slower** |
| `linear/4d/max200/eval100/flipped` | 4.97 µs | 115 µs |  |  |  | 23x faster |
| `linear/4d/max200/eval1000/ascending` | 49.8 µs | 315 µs | 718 µs |  | 39.6 µs | **1.3x slower** |
| `linear/4d/max200/eval1000/flipped` | 49.8 µs | 317 µs |  |  |  | 6.4x faster |
| `linear/4d/max200/eval10000/ascending` | 498 µs | 2.98 ms | 7.77 ms |  | 372 µs | **1.3x slower** |
| `linear/4d/max200/eval10000/flipped` | 501 µs | 2.86 ms |  |  |  | 5.7x faster |
| `nearest/2d/each40/eval10000/ascending` | 152 µs | 654 µs | 1.48 ms |  | 139 µs | **1.1x slower** |
| `nearest/3d/each40/eval10000/ascending` | 212 µs | 993 µs | 1.76 ms |  | 209 µs | **1x slower** |
| `pchip/2d/each40/eval10000/ascending` | 554 µs | 1.29 s |  |  |  | 2.3e+03x faster |
| `pchip/3d/each40/eval10000/ascending` | 1.57 ms | 3.44 s |  |  |  | 2.2e+03x faster |
| `cubic/2d/each40/eval10000/ascending` | 399 µs | 1.28 ms |  |  | 439 µs | 1.1x faster |
| `cubic/3d/each40/eval10000/ascending` | 878 µs | 3.2 ms |  |  | 853 µs | **1x slower** |
| `cubic/2d/each10/eval1000/ascending` | 33.6 µs | 142 µs |  |  | 35.1 µs | 1x faster |
| `cubic/2d/each40/eval1000/ascending` | 36.3 µs | 154 µs |  |  | 39.7 µs | 1.1x faster |
| `cubic/2d/each100/eval1000/ascending` | 38 µs | 167 µs |  |  | 42 µs | 1.1x faster |
| `cubic/2d/each200/eval1000/ascending` | 39.7 µs | 202 µs |  |  | 44.8 µs | 1.1x faster |
| `cubic/2d/each400/eval1000/ascending` | 41.9 µs | 295 µs |  |  | 53 µs | 1.3x faster |
| `quintic/2d/each10/eval1000/ascending` | 60.7 µs | 249 µs |  |  |  | 4.1x faster |
| `quintic/2d/each40/eval1000/ascending` | 65.5 µs | 244 µs |  |  |  | 3.7x faster |
| `linear/2d/each400/eval10000/ascending` | 214 µs | 1.04 ms | 4.3 ms | 548 µs | 219 µs | 1x faster |
| `linear/2d/each800/eval10000/ascending` | 257 µs | 1.19 ms | 8.35 ms | 585 µs | 260 µs | 1x faster |

### Construction

| case | interpx | SciPy | scirs2-interpolate | ndarray-interp | interpn | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/2d/max10/eval10/ascending` | 125 ns | 6.48 µs | 141 ns | 89.8 ns | 0 | **1.9e+02x slower** |
| `linear/2d/max10/eval10/flipped` | 219 ns | 13 µs |  |  |  | 59x faster |
| `linear/2d/max10/eval100/ascending` | 127 ns | 6.51 µs | 148 ns | 89.9 ns | 0 | **1.8e+02x slower** |
| `linear/2d/max10/eval100/flipped` | 214 ns | 12.1 µs |  |  |  | 57x faster |
| `linear/2d/max10/eval1000/ascending` | 127 ns | 6.34 µs | 134 ns | 89.8 ns | 0 | **1.9e+02x slower** |
| `linear/2d/max10/eval1000/flipped` | 205 ns | 12 µs |  |  |  | 59x faster |
| `linear/2d/max10/eval10000/ascending` | 127 ns | 6.09 µs | 134 ns | 84.8 ns | 0 | **1.8e+02x slower** |
| `linear/2d/max10/eval10000/flipped` | 205 ns | 12.9 µs |  |  |  | 63x faster |
| `linear/2d/max40/eval10/ascending` | 377 ns | 6.63 µs | 221 ns | 151 ns | 0 | **5.6e+02x slower** |
| `linear/2d/max40/eval10/flipped` | 1.38 µs | 13 µs |  |  |  | 9.4x faster |
| `linear/2d/max40/eval100/ascending` | 382 ns | 6.13 µs | 218 ns | 151 ns | 0 | **5.7e+02x slower** |
| `linear/2d/max40/eval100/flipped` | 1.38 µs | 13.1 µs |  |  |  | 9.5x faster |
| `linear/2d/max40/eval1000/ascending` | 391 ns | 6.22 µs | 220 ns | 164 ns | 0 | **5.8e+02x slower** |
| `linear/2d/max40/eval1000/flipped` | 1.35 µs | 12 µs |  |  |  | 8.9x faster |
| `linear/2d/max40/eval10000/ascending` | 409 ns | 6.04 µs | 224 ns | 165 ns | 0 | **6.1e+02x slower** |
| `linear/2d/max40/eval10000/flipped` | 1.37 µs | 12.2 µs |  |  |  | 8.9x faster |
| `linear/2d/max200/eval10/ascending` | 7.7 µs | 6.09 µs | 3 µs | 474 ns | 0 | **1.2e+04x slower** |
| `linear/2d/max200/eval10/flipped` | 33 µs | 12.3 µs |  |  |  | **2.7x slower** |
| `linear/2d/max200/eval100/ascending` | 8.13 µs | 6.08 µs | 3 µs | 513 ns | 0 | **1.2e+04x slower** |
| `linear/2d/max200/eval100/flipped` | 32.3 µs | 13.3 µs |  |  |  | **2.4x slower** |
| `linear/2d/max200/eval1000/ascending` | 7.6 µs | 6.12 µs | 2.98 µs | 474 ns | 0 | **1.1e+04x slower** |
| `linear/2d/max200/eval1000/flipped` | 32.1 µs | 13.3 µs |  |  |  | **2.4x slower** |
| `linear/2d/max200/eval10000/ascending` | 8.2 µs | 6.58 µs | 3.06 µs | 471 ns | 0 | **1.2e+04x slower** |
| `linear/2d/max200/eval10000/flipped` | 47.4 µs | 12.3 µs |  |  |  | **3.8x slower** |
| `linear/3d/max10/eval10/ascending` | 181 ns | 7.83 µs | 168 ns |  | 0 | **2e+02x slower** |
| `linear/3d/max10/eval10/flipped` | 446 ns | 16 µs |  |  |  | 36x faster |
| `linear/3d/max10/eval100/ascending` | 189 ns | 7.81 µs | 168 ns |  | 0 | **2.1e+02x slower** |
| `linear/3d/max10/eval100/flipped` | 456 ns | 16.4 µs |  |  |  | 36x faster |
| `linear/3d/max10/eval1000/ascending` | 182 ns | 7.76 µs | 168 ns |  | 0 | **2e+02x slower** |
| `linear/3d/max10/eval1000/flipped` | 440 ns | 16.4 µs |  |  |  | 37x faster |
| `linear/3d/max10/eval10000/ascending` | 186 ns | 8.12 µs | 178 ns |  | 0 | **2.1e+02x slower** |
| `linear/3d/max10/eval10000/flipped` | 438 ns | 17.3 µs |  |  |  | 39x faster |
| `linear/3d/max40/eval10/ascending` | 3.7 µs | 8.18 µs | 1.32 µs |  | 0 | **3.9e+03x slower** |
| `linear/3d/max40/eval10/flipped` | 17.4 µs | 16.1 µs |  |  |  | **1.1x slower** |
| `linear/3d/max40/eval100/ascending` | 5 µs | 7.87 µs | 1.33 µs |  | 0 | **5.6e+03x slower** |
| `linear/3d/max40/eval100/flipped` | 18.7 µs | 16.2 µs |  |  |  | **1.2x slower** |
| `linear/3d/max40/eval1000/ascending` | 4.39 µs | 7.82 µs | 1.32 µs |  | 0 | **4.9e+03x slower** |
| `linear/3d/max40/eval1000/flipped` | 18 µs | 16.1 µs |  |  |  | **1.1x slower** |
| `linear/3d/max40/eval10000/ascending` | 4.01 µs | 7.88 µs | 1.32 µs |  | 0 | **4.4e+03x slower** |
| `linear/3d/max40/eval10000/flipped` | 18.4 µs | 16.1 µs |  |  |  | **1.1x slower** |
| `linear/3d/max200/eval10/ascending` | 1.86 ms | 8.49 µs | 1.31 ms |  | 0 | **2.1e+06x slower** |
| `linear/3d/max200/eval10/flipped` | 4.19 ms | 16.3 µs |  |  |  | **2.6e+02x slower** |
| `linear/3d/max200/eval100/ascending` | 1.89 ms | 7.91 µs | 1.29 ms |  | 0 | **2.1e+06x slower** |
| `linear/3d/max200/eval100/flipped` | 3.99 ms | 16.3 µs |  |  |  | **2.4e+02x slower** |
| `linear/3d/max200/eval1000/ascending` | 1.8 ms | 7.8 µs | 1.28 ms |  | 0 | **2e+06x slower** |
| `linear/3d/max200/eval1000/flipped` | 3.88 ms | 16.5 µs |  |  |  | **2.4e+02x slower** |
| `linear/3d/max200/eval10000/ascending` | 1.82 ms | 7.81 µs | 1.22 ms |  | 0 | **2e+06x slower** |
| `linear/3d/max200/eval10000/flipped` | 3.93 ms | 16.4 µs |  |  |  | **2.4e+02x slower** |
| `linear/4d/max10/eval10/ascending` |  | 10 µs |  |  |  |  |
| `linear/4d/max10/eval10/flipped` |  | 18.4 µs |  |  |  |  |
| `linear/4d/max10/eval100/ascending` |  | 10.8 µs |  |  |  |  |
| `linear/4d/max10/eval100/flipped` |  | 18.3 µs |  |  |  |  |
| `linear/4d/max10/eval1000/ascending` |  | 9.89 µs |  |  |  |  |
| `linear/4d/max10/eval1000/flipped` |  | 18.2 µs |  |  |  |  |
| `linear/4d/max10/eval10000/ascending` |  | 10 µs |  |  |  |  |
| `linear/4d/max10/eval10000/flipped` |  | 18.3 µs |  |  |  |  |
| `linear/4d/max40/eval10/ascending` | 23.1 µs | 9.65 µs | 5.32 µs |  | 1.11 ns | **2.1e+04x slower** |
| `linear/4d/max40/eval10/flipped` | 95.3 µs | 20.4 µs |  |  |  | **4.7x slower** |
| `linear/4d/max40/eval100/ascending` | 23.2 µs | 9.67 µs | 5.36 µs |  | 1.14 ns | **2e+04x slower** |
| `linear/4d/max40/eval100/flipped` | 104 µs | 22 µs |  |  |  | **4.7x slower** |
| `linear/4d/max40/eval1000/ascending` | 22.9 µs | 10.4 µs | 5.5 µs |  | 1.16 ns | **2e+04x slower** |
| `linear/4d/max40/eval1000/flipped` | 94.8 µs | 20.3 µs |  |  |  | **4.7x slower** |
| `linear/4d/max40/eval10000/ascending` | 23.1 µs | 9.67 µs | 5.91 µs |  | 1.11 ns | **2.1e+04x slower** |
| `linear/4d/max40/eval10000/flipped` | 95.3 µs | 20.5 µs |  |  |  | **4.7x slower** |
| `linear/4d/max200/eval10/ascending` | 47 ms | 9.64 µs | 40.7 ms |  | 1.12 ns | **4.2e+07x slower** |
| `linear/4d/max200/eval10/flipped` | 118 ms | 20.8 µs |  |  |  | **5.7e+03x slower** |
| `linear/4d/max200/eval100/ascending` | 49 ms | 10.4 µs | 38.7 ms |  | 1.11 ns | **4.4e+07x slower** |
| `linear/4d/max200/eval100/flipped` | 120 ms | 20.5 µs |  |  |  | **5.8e+03x slower** |
| `linear/4d/max200/eval1000/ascending` | 50.3 ms | 10.1 µs | 39.2 ms |  | 1.15 ns | **4.4e+07x slower** |
| `linear/4d/max200/eval1000/flipped` | 116 ms | 20.5 µs |  |  |  | **5.7e+03x slower** |
| `linear/4d/max200/eval10000/ascending` | 48.4 ms | 9.69 µs | 38.8 ms |  | 1.11 ns | **4.3e+07x slower** |
| `linear/4d/max200/eval10000/flipped` | 118 ms | 22.2 µs |  |  |  | **5.3e+03x slower** |
| `nearest/2d/each40/eval10000/ascending` | 694 ns | 6.12 µs | 287 ns |  | 0 | **1.1e+03x slower** |
| `nearest/3d/each40/eval10000/ascending` | 26.4 µs | 8.31 µs | 9.51 µs |  | 0 | **3e+04x slower** |
| `pchip/2d/each40/eval10000/ascending` | 6.8 µs | 7.05 µs |  |  |  | 1x faster |
| `pchip/3d/each40/eval10000/ascending` | 507 µs | 9.65 µs |  |  |  | **53x slower** |
| `cubic/2d/each40/eval10000/ascending` | 9.03 µs | 1.94 ms |  |  | 0 | **9.2e+03x slower** |
| `cubic/3d/each40/eval10000/ascending` | 341 µs | 152 ms |  |  | 0 | **3.5e+05x slower** |
| `cubic/2d/each10/eval1000/ascending` | 1.45 µs | 1.08 ms |  |  | 0 | **1.6e+03x slower** |
| `cubic/2d/each40/eval1000/ascending` | 8.87 µs | 1.94 ms |  |  | 0 | **1e+04x slower** |
| `cubic/2d/each100/eval1000/ascending` | 58.2 µs | 7.91 ms |  |  | 0 | **6.5e+04x slower** |
| `cubic/2d/each200/eval1000/ascending` | 190 µs | 28.5 ms |  |  | 0 | **2.1e+05x slower** |
| `cubic/2d/each400/eval1000/ascending` | 1.44 ms | 150 ms |  |  | 0 | **1.6e+06x slower** |
| `quintic/2d/each10/eval1000/ascending` | 1.93 µs | 3.75 ms |  |  |  | 1.9e+03x faster |
| `quintic/2d/each40/eval1000/ascending` | 12.4 µs | 8.4 ms |  |  |  | 6.8e+02x faster |
| `linear/2d/each400/eval10000/ascending` | 329 µs | 6.37 µs | 281 µs | 1.15 µs | 0 | **4.8e+05x slower** |
| `linear/2d/each800/eval10000/ascending` | 1.17 ms | 6.56 µs | 831 µs | 2.24 µs | 0 | **1.8e+06x slower** |

## `NearestNDInterpolator`

`n` points of the unit square, evaluated at the pixels of a square image.

### Evaluation

| case | interpx | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: |
| `n10/image50` | 113 µs | 473 µs | 4.2x faster |
| `n10/image100` | 458 µs | 1.86 ms | 4.1x faster |
| `n10/image500` | 12.5 ms | 46.7 ms | 3.7x faster |
| `n100/image50` | 211 µs | 593 µs | 2.8x faster |
| `n100/image100` | 802 µs | 2.37 ms | 3x faster |
| `n100/image500` | 18.8 ms | 58.4 ms | 3.1x faster |
| `n1000/image50` | 340 µs | 859 µs | 2.5x faster |
| `n1000/image100` | 1.26 ms | 3.2 ms | 2.5x faster |
| `n1000/image500` | 25.9 ms | 71.8 ms | 2.8x faster |
| `n10000/image50` | 398 µs | 1.08 ms | 2.7x faster |
| `n10000/image100` | 1.5 ms | 4.13 ms | 2.7x faster |
| `n10000/image500` | 35 ms | 90.1 ms | 2.6x faster |
| `n100000/image50` | 514 µs | 1.48 ms | 2.9x faster |
| `n100000/image100` | 1.99 ms | 5.68 ms | 2.9x faster |
| `n100000/image500` | 47.2 ms | 130 ms | 2.8x faster |

### Construction

| case | interpx | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: |
| `n10/image50` | 552 ns | 10.7 µs | 19x faster |
| `n10/image100` | 565 ns | 10.8 µs | 19x faster |
| `n10/image500` | 537 ns | 10.6 µs | 20x faster |
| `n100/image50` | 3.31 µs | 15.2 µs | 4.6x faster |
| `n100/image100` | 3.33 µs | 15.6 µs | 4.7x faster |
| `n100/image500` | 3.41 µs | 15.7 µs | 4.6x faster |
| `n1000/image50` | 34.5 µs | 127 µs | 3.7x faster |
| `n1000/image100` | 36.6 µs | 125 µs | 3.4x faster |
| `n1000/image500` | 37 µs | 123 µs | 3.3x faster |
| `n10000/image50` | 703 µs | 1.8 ms | 2.6x faster |
| `n10000/image100` | 663 µs | 1.81 ms | 2.7x faster |
| `n10000/image500` | 654 µs | 1.79 ms | 2.7x faster |
| `n100000/image50` | 8.53 ms | 26.7 ms | 3.1x faster |
| `n100000/image100` | 8.75 ms | 26.2 ms | 3x faster |
| `n100000/image500` | 9.23 ms | 26.5 ms | 2.9x faster |

## `RBFInterpolator`

`n` points of a square, evaluated at 1000 points in no order (`scattered`) or at the 10 000 pixels of an image (`image`). `global`: all points; `k50`: the 50 nearest.

### Evaluation

| case | interpx | SciPy | scirs2-interpolate | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `linear/global/n10/scattered` | 24.1 µs | 61.6 µs | 182 µs | 2.6x faster |
| `linear/k50/n10/scattered` | 24.2 µs | 1.38 ms |  | 57x faster |
| `thin_plate_spline/global/n10/scattered` | 92 µs | 116 µs | 303 µs | 1.3x faster |
| `thin_plate_spline/k50/n10/scattered` | 89.8 µs | 1.52 ms |  | 17x faster |
| `cubic/global/n10/scattered` | 40.6 µs | 59.4 µs | 198 µs | 1.5x faster |
| `cubic/k50/n10/scattered` | 40.9 µs | 1.37 ms |  | 33x faster |
| `quintic/global/n10/scattered` | 75 µs | 66.8 µs | 216 µs | **1.1x slower** |
| `quintic/k50/n10/scattered` | 72.5 µs | 1.54 ms |  | 21x faster |
| `multiquadric/global/n10/scattered` | 40.6 µs | 53.3 µs | 236 µs | 1.3x faster |
| `multiquadric/k50/n10/scattered` | 40.6 µs | 1.43 ms |  | 35x faster |
| `inverse_multiquadric/global/n10/scattered` | 47.6 µs | 56.9 µs | 291 µs | 1.2x faster |
| `inverse_multiquadric/k50/n10/scattered` | 51.2 µs | 1.46 ms |  | 28x faster |
| `inverse_quadratic/global/n10/scattered` | 33.1 µs | 47.4 µs |  | 1.4x faster |
| `inverse_quadratic/k50/n10/scattered` | 33.8 µs | 1.5 ms |  | 44x faster |
| `gaussian/global/n10/scattered` | 79.6 µs | 88.8 µs | 358 µs | 1.1x faster |
| `gaussian/k50/n10/scattered` | 85 µs | 1.51 ms |  | 18x faster |
| `linear/global/n100/scattered` | 153 µs | 292 µs | 1.8 ms | 1.9x faster |
| `linear/k50/n100/scattered` | 14.7 ms | 21.4 ms |  | 1.5x faster |
| `thin_plate_spline/global/n100/scattered` | 748 µs | 900 µs | 3.15 ms | 1.2x faster |
| `thin_plate_spline/k50/n100/scattered` | 25.8 ms | 26.5 ms |  | 1x faster |
| `cubic/global/n100/scattered` | 163 µs | 302 µs | 2.06 ms | 1.9x faster |
| `cubic/k50/n100/scattered` | 16.7 ms | 21.4 ms |  | 1.3x faster |
| `quintic/global/n100/scattered` | 232 µs | 314 µs | 2.16 ms | 1.4x faster |
| `quintic/k50/n100/scattered` | 20.4 ms | 23.9 ms |  | 1.2x faster |
| `multiquadric/global/n100/scattered` | 291 µs | 311 µs | 2.38 ms | 1.1x faster |
| `multiquadric/k50/n100/scattered` | 16.8 ms | 21.4 ms |  | 1.3x faster |
| `inverse_multiquadric/global/n100/scattered` | 381 µs | 398 µs | 2.77 ms | 1x faster |
| `inverse_multiquadric/k50/n100/scattered` | 17.8 ms | 21.6 ms |  | 1.2x faster |
| `inverse_quadratic/global/n100/scattered` | 239 µs | 288 µs |  | 1.2x faster |
| `inverse_quadratic/k50/n100/scattered` | 15.7 ms | 20.2 ms |  | 1.3x faster |
| `gaussian/global/n100/scattered` | 634 µs | 732 µs | 3.71 ms | 1.2x faster |
| `gaussian/k50/n100/scattered` | 21.8 ms | 23.9 ms |  | 1.1x faster |
| `linear/global/n1000/scattered` | 1.51 ms | 3.9 ms | 18.7 ms | 2.6x faster |
| `linear/k50/n1000/scattered` | 15.5 ms | 36.7 ms |  | 2.4x faster |
| `thin_plate_spline/global/n1000/scattered` | 7.95 ms | 10.6 ms | 31.9 ms | 1.3x faster |
| `thin_plate_spline/k50/n1000/scattered` | 24.7 ms | 46 ms |  | 1.9x faster |
| `cubic/global/n1000/scattered` | 1.5 ms | 3.99 ms | 20.1 ms | 2.7x faster |
| `cubic/k50/n1000/scattered` | 17.7 ms | 39.6 ms |  | 2.2x faster |
| `quintic/global/n1000/scattered` | 1.7 ms | 3.94 ms | 22.6 ms | 2.3x faster |
| `quintic/k50/n1000/scattered` | 21 ms | 44.3 ms |  | 2.1x faster |
| `multiquadric/global/n1000/scattered` | 2.83 ms | 4.25 ms | 25.2 ms | 1.5x faster |
| `multiquadric/k50/n1000/scattered` | 18.2 ms | 38.2 ms |  | 2.1x faster |
| `inverse_multiquadric/global/n1000/scattered` | 3.9 ms | 5.45 ms | 28.7 ms | 1.4x faster |
| `inverse_multiquadric/k50/n1000/scattered` | 19.7 ms | 40.4 ms |  | 2x faster |
| `inverse_quadratic/global/n1000/scattered` | 2.47 ms | 3.94 ms |  | 1.6x faster |
| `inverse_quadratic/k50/n1000/scattered` | 17.4 ms | 37.8 ms |  | 2.2x faster |
| `gaussian/global/n1000/scattered` | 6.73 ms | 8.58 ms | 37.9 ms | 1.3x faster |
| `gaussian/k50/n1000/scattered` | 24.7 ms | 45.4 ms |  | 1.8x faster |
| `thin_plate_spline/global/n3000/scattered` | 23.1 ms | 31.8 ms | 97.6 ms | 1.4x faster |
| `thin_plate_spline/k50/n10000/scattered` | 25.2 ms | 47.5 ms |  | 1.9x faster |
| `thin_plate_spline/k50/n10000/image` | 250 ms | 491 ms |  | 2x faster |
| `thin_plate_spline/k50/n100000/scattered` | 26.3 ms | 50 ms |  | 1.9x faster |
| `thin_plate_spline/k50/n100000/image` | 262 ms | 484 ms |  | 1.8x faster |

### Construction

| case | interpx | SciPy | scirs2-interpolate | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `linear/global/n10/scattered` | 1.29 µs | 12.6 µs | 5.58 µs | 4.3x faster |
| `linear/k50/n10/scattered` | 1.27 µs | 17 µs |  | 13x faster |
| `thin_plate_spline/global/n10/scattered` | 1.89 µs | 13.7 µs | 6.66 µs | 3.5x faster |
| `thin_plate_spline/k50/n10/scattered` | 1.89 µs | 17.1 µs |  | 9x faster |
| `cubic/global/n10/scattered` | 1.66 µs | 13.6 µs | 6.11 µs | 3.7x faster |
| `cubic/k50/n10/scattered` | 1.66 µs | 18 µs |  | 11x faster |
| `quintic/global/n10/scattered` | 2.38 µs | 15.6 µs | 5.83 µs | 2.5x faster |
| `quintic/k50/n10/scattered` | 2.46 µs | 17.9 µs |  | 7.3x faster |
| `multiquadric/global/n10/scattered` | 1.36 µs | 13.2 µs | 6.64 µs | 4.9x faster |
| `multiquadric/k50/n10/scattered` | 1.35 µs | 17.1 µs |  | 13x faster |
| `inverse_multiquadric/global/n10/scattered` | 1.41 µs | 12.1 µs | 7.02 µs | 5x faster |
| `inverse_multiquadric/k50/n10/scattered` | 1.38 µs | 17.1 µs |  | 12x faster |
| `inverse_quadratic/global/n10/scattered` | 1.39 µs | 12.6 µs |  | 9x faster |
| `inverse_quadratic/k50/n10/scattered` | 1.32 µs | 17.4 µs |  | 13x faster |
| `gaussian/global/n10/scattered` | 1.48 µs | 12.3 µs | 7.86 µs | 5.3x faster |
| `gaussian/k50/n10/scattered` | 1.52 µs | 17 µs |  | 11x faster |
| `linear/global/n100/scattered` | 44.3 µs | 82.7 µs | 802 µs | 1.9x faster |
| `linear/k50/n100/scattered` | 3.33 µs | 23.6 µs |  | 7.1x faster |
| `thin_plate_spline/global/n100/scattered` | 78.2 µs | 117 µs | 916 µs | 1.5x faster |
| `thin_plate_spline/k50/n100/scattered` | 3.77 µs | 24.1 µs |  | 6.4x faster |
| `cubic/global/n100/scattered` | 51.6 µs | 83.2 µs | 764 µs | 1.6x faster |
| `cubic/k50/n100/scattered` | 3.33 µs | 24.4 µs |  | 7.3x faster |
| `quintic/global/n100/scattered` | 55.6 µs | 101 µs | 819 µs | 1.8x faster |
| `quintic/k50/n100/scattered` | 3.29 µs | 25.8 µs |  | 7.8x faster |
| `multiquadric/global/n100/scattered` | 52.4 µs | 84.5 µs | 664 µs | 1.6x faster |
| `multiquadric/k50/n100/scattered` | 3.41 µs | 23 µs |  | 6.7x faster |
| `inverse_multiquadric/global/n100/scattered` | 56.1 µs | 83.3 µs | 888 µs | 1.5x faster |
| `inverse_multiquadric/k50/n100/scattered` | 3.24 µs | 23.3 µs |  | 7.2x faster |
| `inverse_quadratic/global/n100/scattered` | 48.3 µs | 75.6 µs |  | 1.6x faster |
| `inverse_quadratic/k50/n100/scattered` | 3.31 µs | 22.2 µs |  | 6.7x faster |
| `gaussian/global/n100/scattered` | 67.4 µs | 97.9 µs | 828 µs | 1.5x faster |
| `gaussian/k50/n100/scattered` | 3.15 µs | 22.1 µs |  | 7x faster |
| `linear/global/n1000/scattered` | 17.6 ms | 62.9 ms | 197 ms | 3.6x faster |
| `linear/k50/n1000/scattered` | 37.7 µs | 135 µs |  | 3.6x faster |
| `thin_plate_spline/global/n1000/scattered` | 21.1 ms | 66.5 ms | 212 ms | 3.2x faster |
| `thin_plate_spline/k50/n1000/scattered` | 38.3 µs | 141 µs |  | 3.7x faster |
| `cubic/global/n1000/scattered` | 17.9 ms | 62.2 ms | 196 ms | 3.5x faster |
| `cubic/k50/n1000/scattered` | 38.1 µs | 139 µs |  | 3.6x faster |
| `quintic/global/n1000/scattered` | 18.9 ms | 61.7 ms | 198 ms | 3.3x faster |
| `quintic/k50/n1000/scattered` | 36.3 µs | 140 µs |  | 3.9x faster |
| `multiquadric/global/n1000/scattered` | 19.5 ms | 61.7 ms | 201 ms | 3.2x faster |
| `multiquadric/k50/n1000/scattered` | 37.8 µs | 131 µs |  | 3.5x faster |
| `inverse_multiquadric/global/n1000/scattered` | 18.9 ms | 65.4 ms | 210 ms | 3.5x faster |
| `inverse_multiquadric/k50/n1000/scattered` | 38 µs | 139 µs |  | 3.7x faster |
| `inverse_quadratic/global/n1000/scattered` | 18.8 ms | 60.2 ms |  | 3.2x faster |
| `inverse_quadratic/k50/n1000/scattered` | 38.4 µs | 132 µs |  | 3.4x faster |
| `gaussian/global/n1000/scattered` | 20.5 ms | 64 ms | 211 ms | 3.1x faster |
| `gaussian/k50/n1000/scattered` | 38.2 µs | 136 µs |  | 3.6x faster |
| `thin_plate_spline/global/n3000/scattered` | 423 ms | 1.92 s | 6.18 s | 4.5x faster |
| `thin_plate_spline/k50/n10000/scattered` | 720 µs | 1.71 ms |  | 2.4x faster |
| `thin_plate_spline/k50/n10000/image` | 765 µs | 1.73 ms |  | 2.3x faster |
| `thin_plate_spline/k50/n100000/scattered` | 8.99 ms | 26.1 ms |  | 2.9x faster |
| `thin_plate_spline/k50/n100000/image` | 9.44 ms | 26.3 ms |  | 2.8x faster |

## Memory

Megabytes. "Data read" is the process with the data of the case in memory, before anything is built. The columns after it are what was added to that: the highest while building, with the interpolator built, the highest while evaluating. "Highest" is the highest of the whole run, repeats included.

| case | library | data read | while building | built | while evaluating | highest |
| :--- | :--- | ---: | ---: | ---: | ---: | ---: |
| `oned/cubic_spline/n10000/sorted` | interpx | 3.3 | +1.8 | +1.9 | +1.9 | 5.9 |
| `oned/cubic_spline/n10000/sorted` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |
| `oned/cubic_spline/n10000/sorted` | scirs2-interpolate | 4.3 | +0.8 | +0.9 | +1.5 | 5.9 |
| `oned/cubic_spline/n10000/sorted` | ndarray-interp | 3.3 | +1.7 | +1.6 | +1.6 | 5.4 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpx | 4.3 | +3.9 | +4.0 | +4.0 | 8.4 |
| `grid/cubic/2d/each400/eval1000/ascending` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpn | 5.4 | +0.1 | +0.1 | +0.1 | 5.5 |
| `grid/linear/2d/each800/eval10000/ascending` | interpx | 9.4 | +3.6 | +5.5 | +5.5 | 15.0 |
| `grid/linear/2d/each800/eval10000/ascending` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |
| `grid/linear/2d/each800/eval10000/ascending` | scirs2-interpolate | 12.9 | +0.0 | +1.9 | +2.0 | 14.9 |
| `grid/linear/2d/each800/eval10000/ascending` | ndarray-interp | 12.9 | +0.0 | +0.0 | +0.0 | 12.9 |
| `grid/linear/2d/each800/eval10000/ascending` | interpn | 12.7 | +0.0 | +0.0 | +0.0 | 12.7 |
| `nearest_nd/n100000/image500` | interpx | 12.0 | +3.5 | +2.1 | +4.1 | 18.2 |
| `nearest_nd/n100000/image500` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |
| `rbf/thin_plate_spline/global/n3000/scattered` | interpx | 4.3 | +75.9 | +3.8 | +3.8 | 80.3 |
| `rbf/thin_plate_spline/global/n3000/scattered` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |
| `rbf/thin_plate_spline/global/n3000/scattered` | scirs2-interpolate | 2.8 | +289.8 | +1.7 | +1.7 | 292.5 |
| `rbf/thin_plate_spline/k50/n100000/image` | interpx | 7.2 | +3.0 | +4.2 | +4.5 | 11.9 |
| `rbf/thin_plate_spline/k50/n100000/image` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |

### Over time

The resident set size of each library's process over its run, with the ends of the phases. The time after `evaluated` is the timed repeats.

`oned/cubic_spline/n10000/sorted`

![memory over time, oned/cubic_spline/n10000/sorted](performance/oned-cubic_spline-n10000-sorted.svg)

`grid/cubic/2d/each400/eval1000/ascending`

![memory over time, grid/cubic/2d/each400/eval1000/ascending](performance/grid-cubic-2d-each400-eval1000-ascending.svg)

`grid/linear/2d/each800/eval10000/ascending`

![memory over time, grid/linear/2d/each800/eval10000/ascending](performance/grid-linear-2d-each800-eval10000-ascending.svg)

`nearest_nd/n100000/image500`

![memory over time, nearest_nd/n100000/image500](performance/nearest_nd-n100000-image500.svg)

`rbf/thin_plate_spline/global/n3000/scattered`

![memory over time, rbf/thin_plate_spline/global/n3000/scattered](performance/rbf-thin_plate_spline-global-n3000-scattered.svg)

`rbf/thin_plate_spline/k50/n100000/image`

![memory over time, rbf/thin_plate_spline/k50/n100000/image](performance/rbf-thin_plate_spline-k50-n100000-image.svg)

## Do the libraries compute the same thing

The largest distance of a library's values from SciPy's, over the cases of a family, as a share of the largest value. Under 1e-06 the two compute the same interpolant.

| family | interpx | scirs2-interpolate | ndarray-interp | interpn |
| :--- | ---: | ---: | ---: | ---: |
| One variable | 3e-16 | 6e-02 | 2e-16 |  |
| `PPoly` | 9e-16 |  |  |  |
| `BSpline` | 4e-16 |  |  |  |
| `BSplineCurve` | 4e-16 |  |  |  |
| Regular grid | 4e-04 | 4e-16 | 6e-16 | 7e-01 |
| `NearestNDInterpolator` | 0e+00 |  |  |  |
| `RBFInterpolator` | 8e-05 | 3e+01 |  |  |
