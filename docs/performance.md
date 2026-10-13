# Performance: interpx against SciPy and other Rust crates

<!-- Written by compare/perf/report.py from a run of compare/perf/run.py. Do not edit by hand: run it again. -->

Measured on 2026-10-13, at commit `7b20393a17` with uncommitted changes.

## How it was measured

- Machine: Intel64 Family 6 Model 191 Stepping 2, GenuineIntel, 16 logical cores, 48 GB. Windows-11-10.0.26300-SP0.
- rustc 1.99.0 (b940084d7 2026-09-28); Python 3.13.1, NumPy 2.5.4, SciPy 1.18.1.
- The other Rust crates, at the versions in `compare/perf/runners/*/Cargo.toml`: `scirs2-interpolate` (the widest of the Rust crates: one variable, grids, RBF), `ndarray-interp` (the most used for one variable; also bilinear), `interpn` (grids).
- Every library ran each case in a process of its own, pinned to one core, on one thread, on the same data. The libraries took turns case by case, 3 rounds.
- `interpx, from_owned` is interpx built from data that is handed over (`from_owned`) and so not copied. Its column is there for the types that keep their data as it is given (the grid, `NearestNDInterpolator`, `RBFInterpolator`), next to the construction from borrowed data (`new`), which is the one the comparison counts. Each of its samples is given data made just before it, which is not in the nearest cache as the data of the other constructions is: on small grids it reads slower than `new` for that reason.
- A time is that of one call: construction, or the evaluation of all the queries of the case (`eval_batch_into` for interpx, the call with an array for the others). It is the median of the repeats of a process, and then the median of the rounds.
- interpx is called slower where its median is more than 1.10 times another library's and its fastest repeat is slower than that library's fastest.
- Memory is the resident set size of the process, sampled every millisecond from outside it.
- Command: `perf/run.py --rounds 3 --filter ^(oned|grid)/`.

## Where interpx is slower

| family | operation | cases | interpx fastest | slower than SciPy | slower than a Rust crate |
| :--- | :--- | ---: | ---: | ---: | ---: |
| One variable | evaluation | 112 | 110 | 1 | 1 |
| One variable | construction | 112 | 76 | 8 | 32 |
| `PPoly` | evaluation | 4 | 4 | 0 | 0 |
| `PPoly` | construction | 4 | 4 | 0 | 0 |
| `BSpline` | evaluation | 4 | 4 | 0 | 0 |
| `BSpline` | construction | 4 | 4 | 0 | 0 |
| `BSplineCurve` | evaluation | 2 | 2 | 0 | 0 |
| `BSplineCurve` | construction | 2 | 2 | 0 | 0 |
| Regular grid | evaluation | 87 | 67 | 0 | 20 |
| Regular grid | construction | 87 | 25 | 38 | 43 |
| `NearestNDInterpolator` | evaluation | 15 | 15 | 0 | 0 |
| `NearestNDInterpolator` | construction | 15 | 15 | 0 | 0 |
| `RBFInterpolator` | evaluation | 53 | 52 | 1 | 0 |
| `RBFInterpolator` | construction | 53 | 53 | 0 | 0 |

196 comparisons (a case, an operation, and a library), by method. "Same values" says whether that library's values are SciPy's up to rounding: where they are not, it computes another interpolant under the same name, and its time is not that of the same work.

| method | operation | against | same values | cases | slower by | the largest |
| :--- | :--- | :--- | :--- | ---: | ---: | :--- |
| `grid/linear` | construction | interpn | yes | 34 | 173x to 4.68e+07x | `grid/linear/4d/max200/eval10000/ascending`: 59.4 ms against 1.27 ns |
| `grid/cubic` | construction | interpn | no | 7 | 1.68e+03x to 1.76e+06x | `grid/cubic/2d/each400/eval1000/ascending`: 1.54 ms against 0 |
| `grid/nearest` | construction | interpn | yes | 2 | 824x to 2.36e+04x | `grid/nearest/3d/each40/eval10000/ascending`: 23.1 µs against 0 |
| `grid/linear` | construction | SciPy | yes | 36 | 1.18x to 6.1e+03x | `grid/linear/4d/max200/eval100/flipped`: 123 ms against 20.2 µs |
| `grid/linear` | construction | ndarray-interp | yes | 14 | 1.45x to 598x | `grid/linear/2d/each800/eval10000/ascending`: 1.37 ms against 2.29 µs |
| `grid/pchip` | construction | SciPy | yes | 1 | 57.3x | `grid/pchip/3d/each40/eval10000/ascending`: 510 µs against 8.9 µs |
| `oned/linear` | construction | ndarray-interp | yes | 14 | 1.35x to 7.18x | `oned/linear/n10000/sorted`: 107 µs against 14.9 µs |
| `oned/bspline1` | construction | SciPy | yes | 2 | 3.32x to 3.35x | `oned/bspline1/n10000/sorted`: 81.3 µs against 24.3 µs |
| `grid/nearest` | construction | SciPy | yes | 1 | 2.83x | `grid/nearest/3d/each40/eval10000/ascending`: 23.1 µs against 8.17 µs |
| `oned/pchip` | construction | scirs2-interpolate | yes | 8 | 1.46x to 2.74x | `oned/pchip/n1000/shuffled`: 12.5 µs against 4.58 µs |
| `grid/linear` | construction | scirs2-interpolate | yes | 22 | 1.17x to 2.69x | `grid/linear/4d/max40/eval10000/ascending`: 14.8 µs against 5.49 µs |
| `oned/cubic_spline` | construction | scirs2-interpolate | no | 4 | 1.13x to 2.5x | `oned/cubic_spline/n1000/shuffled`: 43.1 µs against 17.2 µs |
| `grid/nearest` | construction | scirs2-interpolate | yes | 2 | 2.23x to 2.47x | `grid/nearest/3d/each40/eval10000/ascending`: 23.1 µs against 9.35 µs |
| `oned/cubic_spline` | construction | ndarray-interp | yes | 1 | 2.42x | `oned/cubic_spline/n1000/sorted`: 42.2 µs against 17.4 µs |
| `oned/akima` | construction | scirs2-interpolate | yes | 6 | 1.39x to 2.16x | `oned/akima/n1000/sorted`: 12.2 µs against 5.64 µs |
| `oned/linear` | construction | SciPy | yes | 4 | 1.86x to 2.16x | `oned/linear/n10000/sorted`: 107 µs against 49.5 µs |
| `oned/linear` | construction | scirs2-interpolate | yes | 13 | 1.12x to 1.67x | `oned/linear/n1000/shuffled`: 2.13 µs against 1.28 µs |
| `grid/linear` | evaluation | interpn | yes | 19 | 1.1x to 1.5x | `grid/linear/4d/max40/eval100/ascending`: 4.14 µs against 2.76 µs |
| `oned/nearest` | construction | SciPy | yes | 2 | 1.24x to 1.3x | `oned/nearest/n10000/sorted`: 70.4 µs against 54.2 µs |
| `rbf/quintic` | evaluation | SciPy | yes | 1 | 1.12x | `rbf/quintic/global/n10/scattered`: 75 µs against 66.8 µs |
| `oned/nearest` | evaluation | scirs2-interpolate | yes | 1 | 1.12x | `oned/nearest/n10/shuffled`: 302 ns against 270 ns |
| `grid/nearest` | evaluation | interpn | yes | 1 | 1.11x | `grid/nearest/2d/each40/eval10000/ascending`: 151 µs against 136 µs |
| `oned/bspline1` | evaluation | SciPy | yes | 1 | 1.11x | `oned/bspline1/n10000/sorted`: 316 µs against 284 µs |

<details><summary>Every comparison where interpx is slower</summary>

"Values apart" is the distance of that library's values from SciPy's, as a share of the largest value.

| case | operation | against | interpx | the other | slower by | values apart |
| :--- | :--- | :--- | ---: | ---: | ---: | ---: |
| `grid/linear/4d/max200/eval10000/ascending` | construction | interpn | 59.4 ms | 1.27 ns | 46794752.00x | 4e-16 |
| `grid/linear/4d/max200/eval1000/ascending` | construction | interpn | 56.2 ms | 1.27 ns | 44266496.00x | 3e-16 |
| `grid/linear/4d/max200/eval10/ascending` | construction | interpn | 55.5 ms | 1.27 ns | 43687778.46x | 2e-16 |
| `grid/linear/4d/max200/eval100/ascending` | construction | interpn | 55.7 ms | 1.37 ns | 40731794.29x | 3e-16 |
| `grid/linear/3d/max200/eval1000/ascending` | construction | interpn | 2.17 ms | 0 | 2222489.60x | 3e-16 |
| `grid/linear/3d/max200/eval10000/ascending` | construction | interpn | 1.99 ms | 0 | 2034483.20x | 5e-16 |
| `grid/linear/3d/max200/eval100/ascending` | construction | interpn | 1.88 ms | 0 | 1928704.00x | 3e-16 |
| `grid/linear/3d/max200/eval10/ascending` | construction | interpn | 1.86 ms | 0 | 1904640.00x | 3e-16 |
| `grid/cubic/2d/each400/eval1000/ascending` | construction | interpn | 1.54 ms | 0 | 1757184.00x | 5e-01 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | interpn | 1.37 ms | 0 | 1748864.00x | 3e-16 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | interpn | 376 µs | 0 | 480832.00x | 3e-16 |
| `grid/cubic/3d/each40/eval10000/ascending` | construction | interpn | 391 µs | 0 | 400179.20x | 5e-01 |
| `grid/cubic/2d/each200/eval1000/ascending` | construction | interpn | 204 µs | 0 | 232106.67x | 2e-01 |
| `grid/cubic/2d/each100/eval1000/ascending` | construction | interpn | 87.3 µs | 0 | 99299.56x | 2e-01 |
| `grid/linear/2d/max200/eval10000/ascending` | construction | interpn | 49.8 µs | 0 | 72832.00x | 3e-16 |
| `grid/linear/3d/max40/eval1000/ascending` | construction | interpn | 23.7 µs | 0 | 24262.40x | 2e-16 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | interpn | 23.1 µs | 0 | 23628.80x | 0e+00 |
| `grid/linear/3d/max40/eval10000/ascending` | construction | interpn | 17.5 µs | 0 | 17881.60x | 3e-16 |
| `grid/cubic/2d/each40/eval1000/ascending` | construction | interpn | 9.58 µs | 0 | 12260.00x | 7e-01 |
| `grid/linear/4d/max40/eval1000/ascending` | construction | interpn | 14.4 µs | 1.27 ns | 11313.23x | 6e-16 |
| `grid/linear/2d/max200/eval100/ascending` | construction | interpn | 7.5 µs | 0 | 10969.14x | 1e-16 |
| `grid/linear/4d/max40/eval10/ascending` | construction | interpn | 14.3 µs | 1.32 ns | 10884.74x | 2e-16 |
| `grid/cubic/2d/each40/eval10000/ascending` | construction | interpn | 9.55 µs | 0 | 10871.11x | 7e-01 |
| `grid/linear/4d/max40/eval10000/ascending` | construction | interpn | 14.8 µs | 1.37 ns | 10825.14x | 5e-16 |
| `grid/linear/4d/max40/eval100/ascending` | construction | interpn | 14.3 µs | 1.37 ns | 10496.00x | 2e-16 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | interpn | 7.55 µs | 0 | 9664.00x | 2e-16 |
| `grid/linear/2d/max200/eval10/ascending` | construction | interpn | 7.47 µs | 0 | 9564.00x | 2e-16 |
| `grid/linear/4d/max200/eval100/flipped` | construction | SciPy | 123 ms | 20.2 µs | 6098.62x |  |
| `grid/linear/4d/max200/eval1000/flipped` | construction | SciPy | 121 ms | 20.1 µs | 6006.90x |  |
| `grid/linear/4d/max200/eval10/ascending` | construction | SciPy | 55.5 ms | 9.49 µs | 5844.94x |  |
| `grid/linear/4d/max200/eval100/ascending` | construction | SciPy | 55.7 ms | 9.6 µs | 5802.72x |  |
| `grid/linear/4d/max200/eval10000/ascending` | construction | SciPy | 59.4 ms | 10.2 µs | 5799.38x |  |
| `grid/linear/4d/max200/eval10000/flipped` | construction | SciPy | 127 ms | 21.9 µs | 5781.11x |  |
| `grid/linear/4d/max200/eval10/flipped` | construction | SciPy | 123 ms | 21.8 µs | 5612.54x |  |
| `grid/linear/4d/max200/eval1000/ascending` | construction | SciPy | 56.2 ms | 10.4 µs | 5424.82x |  |
| `grid/linear/3d/max40/eval100/ascending` | construction | interpn | 2.98 µs | 0 | 3386.67x | 2e-16 |
| `grid/linear/3d/max40/eval10/ascending` | construction | interpn | 2.96 µs | 0 | 3030.40x | 2e-16 |
| `grid/cubic/2d/each10/eval1000/ascending` | construction | interpn | 1.64 µs | 0 | 1677.60x | 3e-01 |
| `grid/nearest/2d/each40/eval10000/ascending` | construction | interpn | 644 ns | 0 | 824.00x | 0e+00 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | ndarray-interp | 1.37 ms | 2.29 µs | 597.70x | 6e-16 |
| `grid/linear/2d/max40/eval1000/ascending` | construction | interpn | 387 ns | 0 | 565.43x | 3e-16 |
| `grid/linear/2d/max40/eval100/ascending` | construction | interpn | 390 ns | 0 | 498.75x | 2e-16 |
| `grid/linear/2d/max40/eval10/ascending` | construction | interpn | 388 ns | 0 | 497.25x | 2e-16 |
| `grid/linear/2d/max40/eval10000/ascending` | construction | interpn | 385 ns | 0 | 492.25x | 4e-16 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | ndarray-interp | 376 µs | 1.16 µs | 325.00x | 5e-16 |
| `grid/linear/3d/max200/eval1000/ascending` | construction | SciPy | 2.17 ms | 7.67 µs | 282.90x |  |
| `grid/linear/3d/max200/eval10000/ascending` | construction | SciPy | 1.99 ms | 7.63 µs | 260.46x |  |
| `grid/linear/3d/max200/eval100/ascending` | construction | SciPy | 1.88 ms | 7.63 µs | 246.71x |  |
| `grid/linear/3d/max200/eval10/ascending` | construction | SciPy | 1.86 ms | 7.6 µs | 244.89x |  |
| `grid/linear/3d/max200/eval1000/flipped` | construction | SciPy | 3.59 ms | 16.2 µs | 221.08x |  |
| `grid/linear/3d/max10/eval100/ascending` | construction | interpn | 194 ns | 0 | 220.89x | 3e-16 |
| `grid/linear/3d/max200/eval10/flipped` | construction | SciPy | 3.7 ms | 17 µs | 217.48x |  |
| `grid/linear/3d/max200/eval10000/flipped` | construction | SciPy | 3.65 ms | 17.5 µs | 208.84x |  |
| `grid/linear/3d/max200/eval100/flipped` | construction | SciPy | 3.53 ms | 17.4 µs | 202.82x |  |
| `grid/linear/3d/max10/eval1000/ascending` | construction | interpn | 195 ns | 0 | 200.10x | 5e-16 |
| `grid/linear/3d/max10/eval10/ascending` | construction | interpn | 195 ns | 0 | 199.30x | 2e-16 |
| `grid/linear/3d/max10/eval10000/ascending` | construction | interpn | 195 ns | 0 | 199.30x | 4e-16 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | SciPy | 1.37 ms | 6.95 µs | 196.68x |  |
| `grid/linear/2d/max10/eval10000/ascending` | construction | interpn | 134 ns | 0 | 196.43x | 3e-16 |
| `grid/linear/2d/max10/eval100/ascending` | construction | interpn | 123 ns | 0 | 180.29x | 2e-16 |
| `grid/linear/2d/max10/eval1000/ascending` | construction | interpn | 139 ns | 0 | 177.88x | 3e-16 |
| `grid/linear/2d/max10/eval10/ascending` | construction | interpn | 135 ns | 0 | 172.88x | 1e-16 |
| `grid/linear/2d/max200/eval10000/ascending` | construction | ndarray-interp | 49.8 µs | 480 ns | 103.66x | 5e-16 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | SciPy | 376 µs | 6.15 µs | 61.08x |  |
| `grid/pchip/3d/each40/eval10000/ascending` | construction | SciPy | 510 µs | 8.9 µs | 57.35x |  |
| `grid/linear/2d/max200/eval1000/ascending` | construction | ndarray-interp | 7.55 µs | 476 ns | 15.87x | 4e-16 |
| `grid/linear/2d/max200/eval10/ascending` | construction | ndarray-interp | 7.47 µs | 477 ns | 15.66x | 2e-16 |
| `grid/linear/2d/max200/eval100/ascending` | construction | ndarray-interp | 7.5 µs | 516 ns | 14.54x | 3e-16 |
| `grid/linear/2d/max200/eval10000/ascending` | construction | SciPy | 49.8 µs | 5.82 µs | 8.55x |  |
| `oned/linear/n10000/sorted` | construction | ndarray-interp | 107 µs | 14.9 µs | 7.18x | 0e+00 |
| `oned/linear/n10000/shuffled` | construction | ndarray-interp | 102 µs | 14.9 µs | 6.86x | 0e+00 |
| `grid/linear/2d/max200/eval10000/flipped` | construction | SciPy | 54.7 µs | 11.9 µs | 4.61x |  |
| `oned/bspline1/n10000/sorted` | construction | SciPy | 81.3 µs | 24.3 µs | 3.35x |  |
| `oned/bspline1/n10000/shuffled` | construction | SciPy | 77.5 µs | 23.4 µs | 3.32x |  |
| `oned/linear/n10000/sorted-uneven` | construction | ndarray-interp | 99 µs | 33.4 µs | 2.97x | 0e+00 |
| `oned/linear/n10000/shuffled-uneven` | construction | ndarray-interp | 97.5 µs | 33.2 µs | 2.94x | 0e+00 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | SciPy | 23.1 µs | 8.17 µs | 2.83x |  |
| `oned/pchip/n1000/shuffled` | construction | scirs2-interpolate | 12.5 µs | 4.58 µs | 2.74x | 2e-16 |
| `grid/linear/4d/max40/eval10000/ascending` | construction | scirs2-interpolate | 14.8 µs | 5.49 µs | 2.69x | 4e-16 |
| `grid/linear/4d/max40/eval10000/flipped` | construction | SciPy | 54.1 µs | 20.3 µs | 2.67x |  |
| `grid/linear/4d/max40/eval100/flipped` | construction | SciPy | 56.4 µs | 21.3 µs | 2.64x |  |
| `grid/linear/4d/max40/eval10/ascending` | construction | scirs2-interpolate | 14.3 µs | 5.44 µs | 2.64x | 2e-16 |
| `grid/linear/4d/max40/eval1000/flipped` | construction | SciPy | 52.3 µs | 19.8 µs | 2.63x |  |
| `grid/linear/4d/max40/eval10/flipped` | construction | SciPy | 56.3 µs | 21.7 µs | 2.59x |  |
| `oned/pchip/n1000/sorted` | construction | scirs2-interpolate | 11.8 µs | 4.55 µs | 2.59x | 2e-16 |
| `grid/linear/2d/max40/eval1000/ascending` | construction | ndarray-interp | 387 ns | 150 ns | 2.57x | 4e-16 |
| `grid/linear/4d/max40/eval100/ascending` | construction | scirs2-interpolate | 14.3 µs | 5.71 µs | 2.51x | 3e-16 |
| `oned/cubic_spline/n1000/shuffled` | construction | scirs2-interpolate | 43.1 µs | 17.2 µs | 2.50x | 4e-04 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | scirs2-interpolate | 7.55 µs | 3.06 µs | 2.47x | 2e-16 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | scirs2-interpolate | 23.1 µs | 9.35 µs | 2.47x | 0e+00 |
| `grid/linear/4d/max40/eval1000/ascending` | construction | scirs2-interpolate | 14.4 µs | 5.89 µs | 2.44x | 2e-16 |
| `oned/cubic_spline/n1000/sorted` | construction | ndarray-interp | 42.2 µs | 17.4 µs | 2.42x | 2e-16 |
| `grid/linear/2d/max40/eval10/ascending` | construction | ndarray-interp | 388 ns | 164 ns | 2.37x | 2e-16 |
| `grid/linear/2d/max40/eval100/ascending` | construction | ndarray-interp | 390 ns | 165 ns | 2.37x | 2e-16 |
| `grid/linear/2d/max200/eval100/ascending` | construction | scirs2-interpolate | 7.5 µs | 3.2 µs | 2.34x | 2e-16 |
| `grid/linear/2d/max40/eval10000/ascending` | construction | ndarray-interp | 385 ns | 164 ns | 2.34x | 4e-16 |
| `grid/linear/2d/max200/eval10/ascending` | construction | scirs2-interpolate | 7.47 µs | 3.25 µs | 2.30x | 1e-16 |
| `oned/cubic_spline/n1000/sorted` | construction | scirs2-interpolate | 42.2 µs | 18.8 µs | 2.24x | 4e-04 |
| `grid/nearest/2d/each40/eval10000/ascending` | construction | scirs2-interpolate | 644 ns | 289 ns | 2.23x | 0e+00 |
| `oned/linear/n10/sorted` | construction | ndarray-interp | 94.6 ns | 43.7 ns | 2.16x | 0e+00 |
| `oned/akima/n1000/sorted` | construction | scirs2-interpolate | 12.2 µs | 5.64 µs | 2.16x | 3e-17 |
| `oned/linear/n10000/sorted` | construction | SciPy | 107 µs | 49.5 µs | 2.16x |  |
| `oned/linear/n10000/sorted-uneven` | construction | SciPy | 99 µs | 46.5 µs | 2.13x |  |
| `oned/linear/n10000/shuffled` | construction | SciPy | 102 µs | 48.9 µs | 2.09x |  |
| `oned/linear/n10/shuffled` | construction | ndarray-interp | 94.9 ns | 47.5 ns | 2.00x | 0e+00 |
| `grid/linear/3d/max40/eval10000/flipped` | construction | SciPy | 31 µs | 16 µs | 1.94x |  |
| `oned/linear/n10000/shuffled-uneven` | construction | SciPy | 97.5 µs | 52.3 µs | 1.86x |  |
| `grid/linear/3d/max40/eval1000/flipped` | construction | SciPy | 29.4 µs | 16 µs | 1.84x |  |
| `grid/linear/2d/max40/eval10/ascending` | construction | scirs2-interpolate | 388 ns | 215 ns | 1.81x | 4e-17 |
| `oned/pchip/n100/sorted-uneven` | construction | scirs2-interpolate | 969 ns | 536 ns | 1.81x | 2e-16 |
| `oned/akima/n100/sorted` | construction | scirs2-interpolate | 918 ns | 536 ns | 1.71x | 1e-16 |
| `grid/linear/2d/max40/eval100/ascending` | construction | scirs2-interpolate | 390 ns | 229 ns | 1.70x | 3e-16 |
| `grid/linear/2d/max40/eval10000/ascending` | construction | scirs2-interpolate | 385 ns | 229 ns | 1.68x | 2e-16 |
| `oned/pchip/n100/shuffled` | construction | scirs2-interpolate | 969 ns | 580 ns | 1.67x | 2e-16 |
| `oned/linear/n1000/shuffled` | construction | scirs2-interpolate | 2.13 µs | 1.28 µs | 1.67x | 0e+00 |
| `oned/linear/n1000/sorted` | construction | scirs2-interpolate | 2.12 µs | 1.28 µs | 1.66x | 0e+00 |
| `oned/pchip/n100/sorted` | construction | scirs2-interpolate | 968 ns | 584 ns | 1.66x | 2e-16 |
| `oned/akima/n100/shuffled` | construction | scirs2-interpolate | 912 ns | 550 ns | 1.66x | 1e-16 |
| `grid/linear/2d/max40/eval1000/ascending` | construction | scirs2-interpolate | 387 ns | 233 ns | 1.66x | 3e-16 |
| `oned/akima/n50/sorted` | construction | scirs2-interpolate | 530 ns | 324 ns | 1.63x | 1e-16 |
| `grid/linear/2d/max10/eval1000/ascending` | construction | ndarray-interp | 139 ns | 85.4 ns | 1.63x | 3e-16 |
| `oned/linear/n10000/sorted` | construction | scirs2-interpolate | 107 µs | 66 µs | 1.62x | 0e+00 |
| `oned/linear/n10000/shuffled` | construction | scirs2-interpolate | 102 µs | 63.3 µs | 1.61x | 0e+00 |
| `oned/pchip/n100/shuffled-uneven` | construction | scirs2-interpolate | 914 ns | 574 ns | 1.59x | 2e-16 |
| `oned/akima/n1000/shuffled` | construction | scirs2-interpolate | 12.9 µs | 8.11 µs | 1.59x | 3e-17 |
| `grid/linear/2d/max10/eval10/ascending` | construction | ndarray-interp | 135 ns | 85.2 ns | 1.59x | 1e-16 |
| `oned/linear/n10000/sorted-uneven` | construction | scirs2-interpolate | 99 µs | 62.7 µs | 1.58x | 6e-17 |
| `oned/linear/n100/sorted-uneven` | construction | ndarray-interp | 266 ns | 169 ns | 1.58x | 0e+00 |
| `grid/linear/2d/max10/eval10000/ascending` | construction | ndarray-interp | 134 ns | 85.3 ns | 1.57x | 4e-16 |
| `grid/linear/4d/max40/eval10000/ascending` | construction | SciPy | 14.8 µs | 9.42 µs | 1.57x |  |
| `grid/linear/4d/max40/eval1000/ascending` | construction | SciPy | 14.4 µs | 9.36 µs | 1.54x |  |
| `oned/linear/n100/shuffled` | construction | ndarray-interp | 265 ns | 173 ns | 1.53x | 0e+00 |
| `oned/pchip/n50/sorted` | construction | scirs2-interpolate | 543 ns | 356 ns | 1.52x | 1e-16 |
| `grid/linear/4d/max40/eval10/ascending` | construction | SciPy | 14.3 µs | 9.43 µs | 1.52x |  |
| `grid/linear/2d/each800/eval10000/ascending` | construction | scirs2-interpolate | 1.37 ms | 899 µs | 1.52x | 2e-16 |
| `oned/linear/n10000/shuffled-uneven` | construction | scirs2-interpolate | 97.5 µs | 64.3 µs | 1.52x | 6e-17 |
| `grid/linear/4d/max40/eval100/ascending` | evaluation | interpn | 4.14 µs | 2.76 µs | 1.50x | 2e-16 |
| `oned/linear/n50/sorted` | construction | ndarray-interp | 165 ns | 110 ns | 1.49x | 0e+00 |
| `grid/linear/3d/max200/eval1000/ascending` | construction | scirs2-interpolate | 2.17 ms | 1.46 ms | 1.49x | 3e-16 |
| `oned/linear/n1000/shuffled` | construction | ndarray-interp | 2.13 µs | 1.45 µs | 1.47x | 0e+00 |
| `oned/pchip/n50/shuffled` | construction | scirs2-interpolate | 508 ns | 347 ns | 1.46x | 1e-16 |
| `oned/linear/n50/shuffled` | construction | ndarray-interp | 161 ns | 111 ns | 1.45x | 0e+00 |
| `grid/linear/2d/max10/eval100/ascending` | construction | ndarray-interp | 123 ns | 85.3 ns | 1.45x | 4e-16 |
| `grid/linear/3d/max10/eval10000/ascending` | evaluation | interpn | 209 µs | 145 µs | 1.44x | 4e-16 |
| `oned/linear/n100/sorted` | construction | ndarray-interp | 266 ns | 185 ns | 1.44x | 0e+00 |
| `grid/linear/3d/max200/eval10000/ascending` | construction | scirs2-interpolate | 1.99 ms | 1.39 ms | 1.42x | 3e-16 |
| `grid/linear/4d/max40/eval10/ascending` | evaluation | interpn | 421 ns | 298 ns | 1.41x | 2e-16 |
| `grid/linear/2d/max10/eval10/ascending` | evaluation | interpn | 146 ns | 104 ns | 1.40x | 1e-16 |
| `grid/linear/4d/max40/eval100/ascending` | construction | SciPy | 14.3 µs | 10.2 µs | 1.40x |  |
| `grid/linear/4d/max40/eval10000/ascending` | evaluation | interpn | 415 µs | 296 µs | 1.40x | 5e-16 |
| `oned/linear/n1000/sorted` | construction | ndarray-interp | 2.12 µs | 1.51 µs | 1.40x | 0e+00 |
| `grid/linear/4d/max40/eval1000/ascending` | evaluation | interpn | 41.5 µs | 29.6 µs | 1.40x | 6e-16 |
| `grid/linear/3d/max10/eval10/ascending` | evaluation | interpn | 216 ns | 155 ns | 1.39x | 2e-16 |
| `oned/akima/n50/shuffled` | construction | scirs2-interpolate | 488 ns | 351 ns | 1.39x | 1e-16 |
| `grid/linear/3d/max200/eval100/ascending` | construction | scirs2-interpolate | 1.88 ms | 1.36 ms | 1.38x | 3e-16 |
| `grid/linear/4d/max200/eval100/ascending` | evaluation | interpn | 5.41 µs | 3.95 µs | 1.37x | 3e-16 |
| `grid/linear/4d/max200/eval1000/ascending` | evaluation | interpn | 54.2 µs | 39.6 µs | 1.37x | 3e-16 |
| `grid/linear/4d/max200/eval10000/ascending` | evaluation | interpn | 543 µs | 398 µs | 1.36x | 4e-16 |
| `grid/linear/2d/max10/eval100/ascending` | evaluation | interpn | 1.26 µs | 931 ns | 1.36x | 2e-16 |
| `oned/linear/n100/shuffled-uneven` | construction | ndarray-interp | 249 ns | 184 ns | 1.35x | 0e+00 |
| `grid/linear/2d/max10/eval1000/ascending` | evaluation | interpn | 13.6 µs | 10.1 µs | 1.34x | 3e-16 |
| `grid/linear/2d/max10/eval10000/ascending` | evaluation | interpn | 127 µs | 94.7 µs | 1.34x | 3e-16 |
| `grid/linear/3d/max10/eval100/ascending` | evaluation | interpn | 1.93 µs | 1.44 µs | 1.33x | 3e-16 |
| `oned/linear/n100/shuffled` | construction | scirs2-interpolate | 265 ns | 199 ns | 1.33x | 0e+00 |
| `grid/linear/4d/max200/eval10000/ascending` | construction | scirs2-interpolate | 59.4 ms | 45.3 ms | 1.31x | 4e-16 |
| `grid/linear/2d/max200/eval100/flipped` | construction | SciPy | 15.8 µs | 12.1 µs | 1.31x |  |
| `oned/nearest/n10000/sorted` | construction | SciPy | 70.4 µs | 54.2 µs | 1.30x |  |
| `grid/linear/4d/max200/eval10/ascending` | evaluation | interpn | 506 ns | 390 ns | 1.30x | 2e-16 |
| `grid/linear/2d/max200/eval10/flipped` | construction | SciPy | 15.9 µs | 12.3 µs | 1.29x |  |
| `oned/linear/n100/sorted-uneven` | construction | scirs2-interpolate | 266 ns | 207 ns | 1.28x | 3e-17 |
| `oned/linear/n100/sorted` | construction | scirs2-interpolate | 266 ns | 207 ns | 1.28x | 0e+00 |
| `grid/linear/2d/max200/eval100/ascending` | construction | SciPy | 7.5 µs | 5.92 µs | 1.27x |  |
| `grid/linear/2d/max40/eval1000/ascending` | evaluation | interpn | 16.6 µs | 13.2 µs | 1.26x | 3e-16 |
| `grid/linear/4d/max200/eval100/ascending` | construction | scirs2-interpolate | 55.7 ms | 44.5 ms | 1.25x | 1e-16 |
| `grid/linear/3d/max10/eval1000/ascending` | evaluation | interpn | 19 µs | 15.3 µs | 1.25x | 5e-16 |
| `grid/linear/4d/max200/eval10/ascending` | construction | scirs2-interpolate | 55.5 ms | 44.6 ms | 1.24x | 2e-16 |
| `grid/linear/4d/max200/eval1000/ascending` | construction | scirs2-interpolate | 56.2 ms | 45.2 ms | 1.24x | 3e-16 |
| `oned/nearest/n10000/shuffled` | construction | SciPy | 70.2 µs | 56.6 µs | 1.24x |  |
| `grid/linear/3d/max200/eval10/ascending` | construction | scirs2-interpolate | 1.86 ms | 1.52 ms | 1.22x | 3e-16 |
| `oned/linear/n50/sorted` | construction | scirs2-interpolate | 165 ns | 135 ns | 1.22x | 0e+00 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | SciPy | 7.55 µs | 6.27 µs | 1.20x |  |
| `grid/linear/2d/max200/eval1000/flipped` | construction | SciPy | 15.7 µs | 13 µs | 1.20x |  |
| `oned/linear/n50/shuffled` | construction | scirs2-interpolate | 161 ns | 134 ns | 1.20x | 0e+00 |
| `oned/linear/n100/shuffled-uneven` | construction | scirs2-interpolate | 249 ns | 208 ns | 1.19x | 3e-17 |
| `grid/linear/2d/max200/eval10/ascending` | construction | SciPy | 7.47 µs | 6.31 µs | 1.18x |  |
| `grid/linear/3d/max10/eval1000/ascending` | construction | scirs2-interpolate | 195 ns | 166 ns | 1.18x | 3e-16 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | scirs2-interpolate | 376 µs | 321 µs | 1.17x | 2e-16 |
| `grid/linear/2d/max40/eval100/ascending` | evaluation | interpn | 1.66 µs | 1.45 µs | 1.15x | 2e-16 |
| `oned/cubic_spline/n100/shuffled-uneven` | construction | scirs2-interpolate | 1.74 µs | 1.53 µs | 1.14x | 6e-03 |
| `oned/cubic_spline/n100/sorted-uneven` | construction | scirs2-interpolate | 1.74 µs | 1.54 µs | 1.13x | 6e-03 |
| `oned/linear/n10/shuffled` | construction | scirs2-interpolate | 94.9 ns | 84.5 ns | 1.12x | 0e+00 |
| `rbf/quintic/global/n10/scattered` | evaluation | SciPy | 75 µs | 66.8 µs | 1.12x |  |
| `oned/nearest/n10/shuffled` | evaluation | scirs2-interpolate | 302 ns | 270 ns | 1.12x | 0e+00 |
| `grid/nearest/2d/each40/eval10000/ascending` | evaluation | interpn | 151 µs | 136 µs | 1.11x | 0e+00 |
| `oned/bspline1/n10000/sorted` | evaluation | SciPy | 316 µs | 284 µs | 1.11x |  |
| `grid/linear/2d/max200/eval10000/ascending` | evaluation | interpn | 190 µs | 172 µs | 1.10x | 3e-16 |

</details>

## Runs that failed

| case | library | how |
| :--- | :--- | :--- |
| `grid/linear/4d/max10/eval10/ascending` | interpn | exit code 101 |
| `grid/linear/4d/max10/eval10/ascending` | scirs2-interpolate | exit code 101 |
| `grid/linear/4d/max10/eval100/ascending` | interpn | exit code 101 |
| `grid/linear/4d/max10/eval100/ascending` | scirs2-interpolate | exit code 101 |
| `grid/linear/4d/max10/eval1000/ascending` | interpn | exit code 101 |
| `grid/linear/4d/max10/eval1000/ascending` | scirs2-interpolate | exit code 101 |
| `grid/linear/4d/max10/eval10000/ascending` | interpn | exit code 101 |
| `grid/linear/4d/max10/eval10000/ascending` | scirs2-interpolate | exit code 101 |

## One variable

Data of `n` points, evaluated at `4 n` points: in order (`sorted`) and in no order (`shuffled`).

### Evaluation

| case | interpx | SciPy | scirs2-interpolate | ndarray-interp | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: |
| `linear/n10/sorted` | 112 ns | 5.39 µs | 260 ns | 255 ns | 2.3x faster |
| `linear/n10/shuffled` | 184 ns | 5.06 µs | 284 ns | 264 ns | 1.4x faster |
| `nearest/n10/sorted` | 167 ns | 9.49 µs | 262 ns |  | 1.6x faster |
| `nearest/n10/shuffled` | 302 ns | 9.43 µs | 270 ns |  | **1.1x slower** |
| `previous/n10/sorted` | 160 ns | 11.6 µs |  |  | 73x faster |
| `previous/n10/shuffled` | 286 ns | 11.1 µs |  |  | 39x faster |
| `bspline1/n10/sorted` | 299 ns | 2.8 µs |  |  | 9.4x faster |
| `bspline1/n10/shuffled` | 654 ns | 2.77 µs |  |  | 4.2x faster |
| `bspline2/n10/sorted` | 402 ns | 3.13 µs |  |  | 7.8x faster |
| `bspline2/n10/shuffled` | 456 ns | 3.13 µs |  |  | 6.9x faster |
| `bspline3/n10/sorted` | 576 ns | 3.87 µs |  |  | 6.7x faster |
| `bspline3/n10/shuffled` | 607 ns | 3.99 µs |  |  | 6.6x faster |
| `cubic_spline/n10/sorted` | 153 ns | 3.98 µs | 283 ns | 354 ns | 1.9x faster |
| `cubic_spline/n10/shuffled` | 223 ns | 3.74 µs | 303 ns | 356 ns | 1.4x faster |
| `pchip/n10/sorted` | 151 ns | 3.73 µs | 301 ns |  | 2x faster |
| `pchip/n10/shuffled` | 225 ns | 3.66 µs | 318 ns |  | 1.4x faster |
| `akima/n10/sorted` | 151 ns | 3.96 µs | 274 ns |  | 1.8x faster |
| `akima/n10/shuffled` | 204 ns | 3.76 µs | 262 ns |  | 1.3x faster |
| `makima/n10/sorted` | 143 ns | 3.65 µs |  |  | 26x faster |
| `makima/n10/shuffled` | 204 ns | 3.64 µs |  |  | 18x faster |
| `linear/n50/sorted` | 703 ns | 6.7 µs | 1.57 µs | 1.35 µs | 1.9x faster |
| `linear/n50/shuffled` | 889 ns | 8.11 µs | 1.53 µs | 1.27 µs | 1.4x faster |
| `nearest/n50/sorted` | 821 ns | 11.7 µs | 1.54 µs |  | 1.9x faster |
| `nearest/n50/shuffled` | 1.58 µs | 11.3 µs | 1.73 µs |  | 1.1x faster |
| `previous/n50/sorted` | 758 ns | 12.5 µs |  |  | 17x faster |
| `previous/n50/shuffled` | 1.36 µs | 13.1 µs |  |  | 9.6x faster |
| `bspline1/n50/sorted` | 1.51 µs | 3.96 µs |  |  | 2.6x faster |
| `bspline1/n50/shuffled` | 3.38 µs | 5.07 µs |  |  | 1.5x faster |
| `bspline2/n50/sorted` | 1.98 µs | 6.1 µs |  |  | 3.1x faster |
| `bspline2/n50/shuffled` | 2.72 µs | 6.87 µs |  |  | 2.5x faster |
| `bspline3/n50/sorted` | 2.64 µs | 8.67 µs |  |  | 3.3x faster |
| `bspline3/n50/shuffled` | 3.49 µs | 9.54 µs |  |  | 2.7x faster |
| `cubic_spline/n50/sorted` | 704 ns | 5.16 µs | 4.34 µs | 1.68 µs | 2.4x faster |
| `cubic_spline/n50/shuffled` | 980 ns | 5.47 µs | 4.47 µs | 1.6 µs | 1.6x faster |
| `pchip/n50/sorted` | 722 ns | 5.15 µs | 4.32 µs |  | 6x faster |
| `pchip/n50/shuffled` | 972 ns | 5.48 µs | 3.84 µs |  | 3.9x faster |
| `akima/n50/sorted` | 715 ns | 5.23 µs | 4.15 µs |  | 5.8x faster |
| `akima/n50/shuffled` | 979 ns | 5.59 µs | 4.02 µs |  | 4.1x faster |
| `makima/n50/sorted` | 713 ns | 5.12 µs |  |  | 7.2x faster |
| `makima/n50/shuffled` | 1.06 µs | 5.61 µs |  |  | 5.3x faster |
| `linear/n100/sorted` | 1.35 µs | 7.45 µs | 3.61 µs | 2.48 µs | 1.8x faster |
| `linear/n100/shuffled` | 1.62 µs | 13 µs | 3.63 µs | 2.72 µs | 1.7x faster |
| `nearest/n100/sorted` | 1.6 µs | 14.8 µs | 3.62 µs |  | 2.3x faster |
| `nearest/n100/shuffled` | 3.14 µs | 13.7 µs | 3.55 µs |  | 1.1x faster |
| `previous/n100/sorted` | 1.69 µs | 15 µs |  |  | 8.9x faster |
| `previous/n100/shuffled` | 2.72 µs | 15.4 µs |  |  | 5.6x faster |
| `bspline1/n100/sorted` | 3.21 µs | 5.49 µs |  |  | 1.7x faster |
| `bspline1/n100/shuffled` | 6.19 µs | 11.9 µs |  |  | 1.9x faster |
| `bspline2/n100/sorted` | 4.01 µs | 9.99 µs |  |  | 2.5x faster |
| `bspline2/n100/shuffled` | 5.61 µs | 12.6 µs |  |  | 2.2x faster |
| `bspline3/n100/sorted` | 5.32 µs | 15.1 µs |  |  | 2.8x faster |
| `bspline3/n100/shuffled` | 7.4 µs | 18 µs |  |  | 2.4x faster |
| `cubic_spline/n100/sorted` | 1.42 µs | 7.18 µs | 17.2 µs | 3.19 µs | 2.3x faster |
| `cubic_spline/n100/shuffled` | 1.98 µs | 8.03 µs | 17.6 µs | 3.21 µs | 1.6x faster |
| `pchip/n100/sorted` | 1.48 µs | 7.41 µs | 17.3 µs |  | 5x faster |
| `pchip/n100/shuffled` | 2.13 µs | 7.83 µs | 16.3 µs |  | 3.7x faster |
| `akima/n100/sorted` | 1.47 µs | 6.89 µs | 16.6 µs |  | 4.7x faster |
| `akima/n100/shuffled` | 1.95 µs | 8.12 µs | 17.2 µs |  | 4.2x faster |
| `makima/n100/sorted` | 1.56 µs | 7.35 µs |  |  | 4.7x faster |
| `makima/n100/shuffled` | 1.94 µs | 8.73 µs |  |  | 4.5x faster |
| `linear/n1000/sorted` | 14.3 µs | 29.4 µs | 93 µs | 24.7 µs | 1.7x faster |
| `linear/n1000/shuffled` | 17.7 µs | 156 µs | 201 µs | 24.7 µs | 1.4x faster |
| `nearest/n1000/sorted` | 14.8 µs | 57.1 µs | 120 µs |  | 3.9x faster |
| `nearest/n1000/shuffled` | 32.1 µs | 200 µs | 232 µs |  | 6.2x faster |
| `previous/n1000/sorted` | 14.9 µs | 61.7 µs |  |  | 4.1x faster |
| `previous/n1000/shuffled` | 27.4 µs | 203 µs |  |  | 7.4x faster |
| `bspline1/n1000/sorted` | 31 µs | 31.5 µs |  |  | 1x faster |
| `bspline1/n1000/shuffled` | 67.6 µs | 578 µs |  |  | 8.5x faster |
| `bspline2/n1000/sorted` | 41.2 µs | 75.9 µs |  |  | 1.8x faster |
| `bspline2/n1000/shuffled` | 72.5 µs | 591 µs |  |  | 8.1x faster |
| `bspline3/n1000/sorted` | 53 µs | 129 µs |  |  | 2.4x faster |
| `bspline3/n1000/shuffled` | 96.7 µs | 580 µs |  |  | 6x faster |
| `cubic_spline/n1000/sorted` | 14.1 µs | 40.4 µs | 1.55 ms | 34.5 µs | 2.4x faster |
| `cubic_spline/n1000/shuffled` | 19.5 µs | 219 µs | 1.43 ms | 31.9 µs | 1.6x faster |
| `pchip/n1000/sorted` | 14 µs | 38.6 µs | 1.24 ms |  | 2.8x faster |
| `pchip/n1000/shuffled` | 19.7 µs | 218 µs | 1.24 ms |  | 11x faster |
| `akima/n1000/sorted` | 14 µs | 38 µs | 1.43 ms |  | 2.7x faster |
| `akima/n1000/shuffled` | 19.6 µs | 223 µs | 1.42 ms |  | 11x faster |
| `makima/n1000/sorted` | 14.4 µs | 37.5 µs |  |  | 2.6x faster |
| `makima/n1000/shuffled` | 19.6 µs | 219 µs |  |  | 11x faster |
| `linear/n10000/sorted` | 141 µs | 235 µs | 1.38 ms | 248 µs | 1.7x faster |
| `linear/n10000/shuffled` | 208 µs | 2.49 ms | 2.85 ms | 284 µs | 1.4x faster |
| `nearest/n10000/sorted` | 155 µs | 583 µs | 1.49 ms |  | 3.8x faster |
| `nearest/n10000/shuffled` | 419 µs | 3.17 ms | 3.08 ms |  | 7.3x faster |
| `previous/n10000/sorted` | 165 µs | 604 µs |  |  | 3.7x faster |
| `previous/n10000/shuffled` | 341 µs | 3 ms |  |  | 8.8x faster |
| `bspline1/n10000/sorted` | 316 µs | 284 µs |  |  | **1.1x slower** |
| `bspline1/n10000/shuffled` | 682 µs | 48.3 ms |  |  | 71x faster |
| `bspline2/n10000/sorted` | 422 µs | 762 µs |  |  | 1.8x faster |
| `bspline2/n10000/shuffled` | 1.05 ms | 48.4 ms |  |  | 46x faster |
| `bspline3/n10000/sorted` | 593 µs | 1.37 ms |  |  | 2.3x faster |
| `bspline3/n10000/shuffled` | 1.28 ms | 49.8 ms |  |  | 39x faster |
| `cubic_spline/n10000/sorted` | 154 µs | 380 µs | 148 ms | 343 µs | 2.2x faster |
| `cubic_spline/n10000/shuffled` | 237 µs | 3.7 ms | 147 ms | 365 µs | 1.5x faster |
| `pchip/n10000/sorted` | 160 µs | 391 µs | 122 ms |  | 2.5x faster |
| `pchip/n10000/shuffled` | 224 µs | 3.53 ms | 124 ms |  | 16x faster |
| `akima/n10000/sorted` | 158 µs | 385 µs | 146 ms |  | 2.4x faster |
| `akima/n10000/shuffled` | 219 µs | 3.46 ms | 150 ms |  | 16x faster |
| `makima/n10000/sorted` | 156 µs | 364 µs |  |  | 2.3x faster |
| `makima/n10000/shuffled` | 216 µs | 3.61 ms |  |  | 17x faster |
| `linear/n100/sorted-uneven` | 1.18 µs | 7.39 µs | 3.31 µs | 3.83 µs | 2.8x faster |
| `linear/n100/shuffled-uneven` | 3.02 µs | 12.8 µs | 3.64 µs | 4.05 µs | 1.2x faster |
| `cubic_spline/n100/sorted-uneven` | 1.44 µs | 8.12 µs | 17.7 µs | 4.41 µs | 3.1x faster |
| `cubic_spline/n100/shuffled-uneven` | 3.07 µs | 8.49 µs | 17.6 µs | 4.43 µs | 1.4x faster |
| `pchip/n100/sorted-uneven` | 1.43 µs | 7.43 µs | 15.9 µs |  | 5.2x faster |
| `pchip/n100/shuffled-uneven` | 3.04 µs | 8 µs | 16.3 µs |  | 2.6x faster |
| `linear/n10000/sorted-uneven` | 209 µs | 273 µs | 1.26 ms | 970 µs | 1.3x faster |
| `linear/n10000/shuffled-uneven` | 677 µs | 2.32 ms | 3.04 ms | 1.77 ms | 2.6x faster |
| `cubic_spline/n10000/sorted-uneven` | 235 µs | 452 µs | 144 ms | 1.05 ms | 1.9x faster |
| `cubic_spline/n10000/shuffled-uneven` | 796 µs | 3.67 ms | 149 ms | 1.97 ms | 2.5x faster |
| `pchip/n10000/sorted-uneven` | 232 µs | 466 µs | 123 ms |  | 2x faster |
| `pchip/n10000/shuffled-uneven` | 741 µs | 3.69 ms | 120 ms |  | 5x faster |

### Construction

| case | interpx | SciPy | scirs2-interpolate | ndarray-interp | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: |
| `linear/n10/sorted` | 94.6 ns | 9.54 µs | 88.8 ns | 43.7 ns | **2.2x slower** |
| `linear/n10/shuffled` | 94.9 ns | 9.03 µs | 84.5 ns | 47.5 ns | **2x slower** |
| `nearest/n10/sorted` | 70.9 ns | 9.47 µs | 88.9 ns |  | 1.3x faster |
| `nearest/n10/shuffled` | 67.2 ns | 9.58 µs | 88.9 ns |  | 1.3x faster |
| `previous/n10/sorted` | 70.3 ns | 9.23 µs |  |  | 1.3e+02x faster |
| `previous/n10/shuffled` | 71.4 ns | 9.49 µs |  |  | 1.3e+02x faster |
| `bspline1/n10/sorted` | 116 ns | 16.8 µs |  |  | 1.5e+02x faster |
| `bspline1/n10/shuffled` | 122 ns | 17.2 µs |  |  | 1.4e+02x faster |
| `bspline2/n10/sorted` | 513 ns | 39 µs |  |  | 76x faster |
| `bspline2/n10/shuffled` | 508 ns | 37.9 µs |  |  | 75x faster |
| `bspline3/n10/sorted` | 566 ns | 38.1 µs |  |  | 67x faster |
| `bspline3/n10/shuffled` | 559 ns | 37.9 µs |  |  | 68x faster |
| `cubic_spline/n10/sorted` | 257 ns | 93.2 µs | 425 ns | 513 ns | 1.7x faster |
| `cubic_spline/n10/shuffled` | 235 ns | 86.8 µs | 389 ns | 512 ns | 1.7x faster |
| `pchip/n10/sorted` | 194 ns | 110 µs | 199 ns |  | 1x faster |
| `pchip/n10/shuffled` | 195 ns | 108 µs | 209 ns |  | 1.1x faster |
| `akima/n10/sorted` | 195 ns | 82.7 µs | 197 ns |  | 1x faster |
| `akima/n10/shuffled` | 196 ns | 75.6 µs | 179 ns |  | **1.1x slower** |
| `makima/n10/sorted` | 190 ns | 84.9 µs |  |  | 4.5e+02x faster |
| `makima/n10/shuffled` | 206 ns | 85.4 µs |  |  | 4.1e+02x faster |
| `linear/n50/sorted` | 165 ns | 8.92 µs | 135 ns | 110 ns | **1.5x slower** |
| `linear/n50/shuffled` | 161 ns | 8.81 µs | 134 ns | 111 ns | **1.5x slower** |
| `nearest/n50/sorted` | 103 ns | 9.63 µs | 131 ns |  | 1.3x faster |
| `nearest/n50/shuffled` | 104 ns | 9.68 µs | 143 ns |  | 1.4x faster |
| `previous/n50/sorted` | 95 ns | 10.3 µs |  |  | 1.1e+02x faster |
| `previous/n50/shuffled` | 95 ns | 10.5 µs |  |  | 1.1e+02x faster |
| `bspline1/n50/sorted` | 219 ns | 17 µs |  |  | 78x faster |
| `bspline1/n50/shuffled` | 218 ns | 17 µs |  |  | 78x faster |
| `bspline2/n50/sorted` | 2.01 µs | 42.8 µs |  |  | 21x faster |
| `bspline2/n50/shuffled` | 1.84 µs | 42.2 µs |  |  | 23x faster |
| `bspline3/n50/sorted` | 2.23 µs | 41.5 µs |  |  | 19x faster |
| `bspline3/n50/shuffled` | 2.23 µs | 39.8 µs |  |  | 18x faster |
| `cubic_spline/n50/sorted` | 838 ns | 94.6 µs | 937 ns | 1.08 µs | 1.1x faster |
| `cubic_spline/n50/shuffled` | 848 ns | 87.5 µs | 874 ns | 1.18 µs | 1x faster |
| `pchip/n50/sorted` | 543 ns | 102 µs | 356 ns |  | **1.5x slower** |
| `pchip/n50/shuffled` | 508 ns | 103 µs | 347 ns |  | **1.5x slower** |
| `akima/n50/sorted` | 530 ns | 82.3 µs | 324 ns |  | **1.6x slower** |
| `akima/n50/shuffled` | 488 ns | 76 µs | 351 ns |  | **1.4x slower** |
| `makima/n50/sorted` | 595 ns | 78.9 µs |  |  | 1.3e+02x faster |
| `makima/n50/shuffled` | 563 ns | 78.2 µs |  |  | 1.4e+02x faster |
| `linear/n100/sorted` | 266 ns | 9.02 µs | 207 ns | 185 ns | **1.4x slower** |
| `linear/n100/shuffled` | 265 ns | 9.15 µs | 199 ns | 173 ns | **1.5x slower** |
| `nearest/n100/sorted` | 167 ns | 11 µs | 207 ns |  | 1.2x faster |
| `nearest/n100/shuffled` | 167 ns | 10 µs | 190 ns |  | 1.1x faster |
| `previous/n100/sorted` | 166 ns | 9.92 µs |  |  | 60x faster |
| `previous/n100/shuffled` | 167 ns | 9.77 µs |  |  | 59x faster |
| `bspline1/n100/sorted` | 319 ns | 16.9 µs |  |  | 53x faster |
| `bspline1/n100/shuffled` | 323 ns | 16.9 µs |  |  | 52x faster |
| `bspline2/n100/sorted` | 3.87 µs | 45.6 µs |  |  | 12x faster |
| `bspline2/n100/shuffled` | 3.87 µs | 46.5 µs |  |  | 12x faster |
| `bspline3/n100/sorted` | 4.32 µs | 43.3 µs |  |  | 10x faster |
| `bspline3/n100/shuffled` | 3.95 µs | 43.3 µs |  |  | 11x faster |
| `cubic_spline/n100/sorted` | 1.58 µs | 95.7 µs | 1.53 µs | 1.88 µs | **1x slower** |
| `cubic_spline/n100/shuffled` | 1.72 µs | 87.6 µs | 1.67 µs | 1.88 µs | **1x slower** |
| `pchip/n100/sorted` | 968 ns | 110 µs | 584 ns |  | **1.7x slower** |
| `pchip/n100/shuffled` | 969 ns | 112 µs | 580 ns |  | **1.7x slower** |
| `akima/n100/sorted` | 918 ns | 77.7 µs | 536 ns |  | **1.7x slower** |
| `akima/n100/shuffled` | 912 ns | 76.4 µs | 550 ns |  | **1.7x slower** |
| `makima/n100/sorted` | 1.07 µs | 87 µs |  |  | 81x faster |
| `makima/n100/shuffled` | 1.08 µs | 87.4 µs |  |  | 81x faster |
| `linear/n1000/sorted` | 2.12 µs | 13.2 µs | 1.28 µs | 1.51 µs | **1.7x slower** |
| `linear/n1000/shuffled` | 2.13 µs | 13.5 µs | 1.28 µs | 1.45 µs | **1.7x slower** |
| `nearest/n1000/sorted` | 1.09 µs | 15.9 µs | 1.23 µs |  | 1.1x faster |
| `nearest/n1000/shuffled` | 997 ns | 16 µs | 1.28 µs |  | 1.3x faster |
| `previous/n1000/sorted` | 1.05 µs | 16.4 µs |  |  | 16x faster |
| `previous/n1000/shuffled` | 994 ns | 17.7 µs |  |  | 18x faster |
| `bspline1/n1000/sorted` | 2.27 µs | 18.3 µs |  |  | 8.1x faster |
| `bspline1/n1000/shuffled` | 2.19 µs | 18.6 µs |  |  | 8.5x faster |
| `bspline2/n1000/sorted` | 35.7 µs | 108 µs |  |  | 3x faster |
| `bspline2/n1000/shuffled` | 37.7 µs | 117 µs |  |  | 3.1x faster |
| `bspline3/n1000/sorted` | 69.4 µs | 125 µs |  |  | 1.8x faster |
| `bspline3/n1000/shuffled` | 75.6 µs | 129 µs |  |  | 1.7x faster |
| `cubic_spline/n1000/sorted` | 42.2 µs | 122 µs | 18.8 µs | 17.4 µs | **2.4x slower** |
| `cubic_spline/n1000/shuffled` | 43.1 µs | 112 µs | 17.2 µs | 16 µs | **2.7x slower** |
| `pchip/n1000/sorted` | 11.8 µs | 124 µs | 4.55 µs |  | **2.6x slower** |
| `pchip/n1000/shuffled` | 12.5 µs | 122 µs | 4.58 µs |  | **2.7x slower** |
| `akima/n1000/sorted` | 12.2 µs | 95.1 µs | 5.64 µs |  | **2.2x slower** |
| `akima/n1000/shuffled` | 12.9 µs | 95.4 µs | 8.11 µs |  | **1.6x slower** |
| `makima/n1000/sorted` | 13.4 µs | 101 µs |  |  | 7.6x faster |
| `makima/n1000/shuffled` | 15.8 µs | 101 µs |  |  | 6.4x faster |
| `linear/n10000/sorted` | 107 µs | 49.5 µs | 66 µs | 14.9 µs | **7.2x slower** |
| `linear/n10000/shuffled` | 102 µs | 48.9 µs | 63.3 µs | 14.9 µs | **6.9x slower** |
| `nearest/n10000/sorted` | 70.4 µs | 54.2 µs | 64.2 µs |  | **1.3x slower** |
| `nearest/n10000/shuffled` | 70.2 µs | 56.6 µs | 63.8 µs |  | **1.2x slower** |
| `previous/n10000/sorted` | 69.2 µs | 75.6 µs |  |  | 1.1x faster |
| `previous/n10000/shuffled` | 68.2 µs | 76.7 µs |  |  | 1.1x faster |
| `bspline1/n10000/sorted` | 81.3 µs | 24.3 µs |  |  | **3.4x slower** |
| `bspline1/n10000/shuffled` | 77.5 µs | 23.4 µs |  |  | **3.3x slower** |
| `bspline2/n10000/sorted` | 525 µs | 807 µs |  |  | 1.5x faster |
| `bspline2/n10000/shuffled` | 538 µs | 746 µs |  |  | 1.4x faster |
| `bspline3/n10000/sorted` | 648 µs | 984 µs |  |  | 1.5x faster |
| `bspline3/n10000/shuffled` | 653 µs | 950 µs |  |  | 1.5x faster |
| `cubic_spline/n10000/sorted` | 258 µs | 745 µs | 295 µs | 400 µs | 1.1x faster |
| `cubic_spline/n10000/shuffled` | 259 µs | 800 µs | 294 µs | 394 µs | 1.1x faster |
| `pchip/n10000/sorted` | 91.7 µs | 258 µs | 178 µs |  | 1.9x faster |
| `pchip/n10000/shuffled` | 96.9 µs | 240 µs | 187 µs |  | 1.9x faster |
| `akima/n10000/sorted` | 91.8 µs | 700 µs | 133 µs |  | 1.4x faster |
| `akima/n10000/shuffled` | 104 µs | 678 µs | 133 µs |  | 1.3x faster |
| `makima/n10000/sorted` | 102 µs | 766 µs |  |  | 7.5x faster |
| `makima/n10000/shuffled` | 102 µs | 771 µs |  |  | 7.6x faster |
| `linear/n100/sorted-uneven` | 266 ns | 9.79 µs | 207 ns | 169 ns | **1.6x slower** |
| `linear/n100/shuffled-uneven` | 249 ns | 9.21 µs | 208 ns | 184 ns | **1.3x slower** |
| `cubic_spline/n100/sorted-uneven` | 1.74 µs | 96.7 µs | 1.54 µs | 1.88 µs | **1.1x slower** |
| `cubic_spline/n100/shuffled-uneven` | 1.74 µs | 91.5 µs | 1.53 µs | 1.89 µs | **1.1x slower** |
| `pchip/n100/sorted-uneven` | 969 ns | 112 µs | 536 ns |  | **1.8x slower** |
| `pchip/n100/shuffled-uneven` | 914 ns | 112 µs | 574 ns |  | **1.6x slower** |
| `linear/n10000/sorted-uneven` | 99 µs | 46.5 µs | 62.7 µs | 33.4 µs | **3x slower** |
| `linear/n10000/shuffled-uneven` | 97.5 µs | 52.3 µs | 64.3 µs | 33.2 µs | **2.9x slower** |
| `cubic_spline/n10000/sorted-uneven` | 260 µs | 772 µs | 303 µs | 377 µs | 1.2x faster |
| `cubic_spline/n10000/shuffled-uneven` | 259 µs | 799 µs | 318 µs | 399 µs | 1.2x faster |
| `pchip/n10000/sorted-uneven` | 93.5 µs | 259 µs | 180 µs |  | 1.9x faster |
| `pchip/n10000/shuffled-uneven` | 99.6 µs | 257 µs | 178 µs |  | 1.8x faster |

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

| case | interpx | interpx, `from_owned` | SciPy | scirs2-interpolate | ndarray-interp | interpn | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/2d/max10/eval10/ascending` | 146 ns | 134 ns | 19.2 µs | 859 ns | 165 ns | 104 ns | **1.4x slower** |
| `linear/2d/max10/eval10/flipped` | 146 ns | 146 ns | 21.2 µs |  |  |  | 1.4e+02x faster |
| `linear/2d/max10/eval100/ascending` | 1.26 µs | 1.25 µs | 23.2 µs | 8.47 µs | 1.36 µs | 931 ns | **1.4x slower** |
| `linear/2d/max10/eval100/flipped` | 1.25 µs | 1.37 µs | 21.5 µs |  |  |  | 17x faster |
| `linear/2d/max10/eval1000/ascending` | 13.6 µs | 13.7 µs | 43.9 µs | 89.2 µs | 15.3 µs | 10.1 µs | **1.3x slower** |
| `linear/2d/max10/eval1000/flipped` | 13.6 µs | 13.6 µs | 43.7 µs |  |  |  | 3.2x faster |
| `linear/2d/max10/eval10000/ascending` | 127 µs | 124 µs | 391 µs | 932 µs | 223 µs | 94.7 µs | **1.3x slower** |
| `linear/2d/max10/eval10000/flipped` | 125 µs | 134 µs | 357 µs |  |  |  | 2.9x faster |
| `linear/2d/max40/eval10/ascending` | 161 ns | 160 ns | 19.3 µs | 968 ns | 170 ns | 147 ns | **1.1x slower** |
| `linear/2d/max40/eval10/flipped` | 166 ns | 171 ns | 19.3 µs |  |  |  | 1.2e+02x faster |
| `linear/2d/max40/eval100/ascending` | 1.66 µs | 1.52 µs | 24.3 µs | 10.9 µs | 1.59 µs | 1.45 µs | **1.1x slower** |
| `linear/2d/max40/eval100/flipped` | 1.54 µs | 1.52 µs | 22.1 µs |  |  |  | 14x faster |
| `linear/2d/max40/eval1000/ascending` | 16.6 µs | 16.2 µs | 50.1 µs | 110 µs | 15.3 µs | 13.2 µs | **1.3x slower** |
| `linear/2d/max40/eval1000/flipped` | 15.2 µs | 15.1 µs | 46 µs |  |  |  | 3x faster |
| `linear/2d/max40/eval10000/ascending` | 153 µs | 164 µs | 583 µs | 1.14 ms | 244 µs | 144 µs | **1.1x slower** |
| `linear/2d/max40/eval10000/flipped` | 152 µs | 166 µs | 582 µs |  |  |  | 3.8x faster |
| `linear/2d/max200/eval10/ascending` | 193 ns | 210 ns | 19.5 µs | 1.94 µs | 205 ns | 193 ns | **1x slower** |
| `linear/2d/max200/eval10/flipped` | 210 ns | 210 ns | 19.7 µs |  |  |  | 94x faster |
| `linear/2d/max200/eval100/ascending` | 1.89 µs | 1.91 µs | 24.8 µs | 18.9 µs | 2.17 µs | 1.89 µs | **1x slower** |
| `linear/2d/max200/eval100/flipped` | 2.06 µs | 2.06 µs | 24.8 µs |  |  |  | 12x faster |
| `linear/2d/max200/eval1000/ascending` | 18.9 µs | 18.9 µs | 72.2 µs | 206 µs | 24.7 µs | 18.8 µs | **1x slower** |
| `linear/2d/max200/eval1000/flipped` | 18.9 µs | 20.7 µs | 70.3 µs |  |  |  | 3.7x faster |
| `linear/2d/max200/eval10000/ascending` | 190 µs | 207 µs | 849 µs | 2.15 ms | 468 µs | 172 µs | **1.1x slower** |
| `linear/2d/max200/eval10000/flipped` | 189 µs | 207 µs | 881 µs |  |  |  | 4.7x faster |
| `linear/3d/max10/eval10/ascending` | 216 ns | 197 ns | 56.8 µs | 3.15 µs |  | 155 ns | **1.4x slower** |
| `linear/3d/max10/eval10/flipped` | 209 ns | 214 ns | 59.4 µs |  |  |  | 2.8e+02x faster |
| `linear/3d/max10/eval100/ascending` | 1.93 µs | 1.91 µs | 62.3 µs | 30.8 µs |  | 1.44 µs | **1.3x slower** |
| `linear/3d/max10/eval100/flipped` | 2.07 µs | 1.92 µs | 62.7 µs |  |  |  | 30x faster |
| `linear/3d/max10/eval1000/ascending` | 19 µs | 19 µs | 115 µs | 335 µs |  | 15.3 µs | **1.2x slower** |
| `linear/3d/max10/eval1000/flipped` | 20.7 µs | 19.2 µs | 116 µs |  |  |  | 5.6x faster |
| `linear/3d/max10/eval10000/ascending` | 209 µs | 207 µs | 1.21 ms | 3.36 ms |  | 145 µs | **1.4x slower** |
| `linear/3d/max10/eval10000/flipped` | 193 µs | 207 µs | 1.26 ms |  |  |  | 6.6x faster |
| `linear/3d/max40/eval10/ascending` | 220 ns | 238 ns | 56.6 µs | 3.7 µs |  | 218 ns | **1x slower** |
| `linear/3d/max40/eval10/flipped` | 218 ns | 239 ns | 61.5 µs |  |  |  | 2.8e+02x faster |
| `linear/3d/max40/eval100/ascending` | 2.13 µs | 2.32 µs | 63.6 µs | 34.5 µs |  | 2.17 µs | 1x faster |
| `linear/3d/max40/eval100/flipped` | 2.33 µs | 2.33 µs | 63.9 µs |  |  |  | 27x faster |
| `linear/3d/max40/eval1000/ascending` | 21.5 µs | 23.2 µs | 142 µs | 342 µs |  | 21.6 µs | 1x faster |
| `linear/3d/max40/eval1000/flipped` | 21.3 µs | 21.8 µs | 130 µs |  |  |  | 6.1x faster |
| `linear/3d/max40/eval10000/ascending` | 212 µs | 212 µs | 1.52 ms | 3.72 ms |  | 219 µs | 1x faster |
| `linear/3d/max40/eval10000/flipped` | 214 µs | 213 µs | 1.52 ms |  |  |  | 7.1x faster |
| `linear/3d/max200/eval10/ascending` | 301 ns | 301 ns | 56.9 µs | 4.95 µs |  | 286 ns | **1.1x slower** |
| `linear/3d/max200/eval10/flipped` | 300 ns | 276 ns | 58.4 µs |  |  |  | 1.9e+02x faster |
| `linear/3d/max200/eval100/ascending` | 2.98 µs | 2.97 µs | 64.3 µs | 51.1 µs |  | 2.86 µs | **1x slower** |
| `linear/3d/max200/eval100/flipped` | 2.8 µs | 2.97 µs | 64.5 µs |  |  |  | 23x faster |
| `linear/3d/max200/eval1000/ascending` | 29.8 µs | 29.8 µs | 185 µs | 503 µs |  | 28.6 µs | **1x slower** |
| `linear/3d/max200/eval1000/flipped` | 29.8 µs | 29.6 µs | 185 µs |  |  |  | 6.2x faster |
| `linear/3d/max200/eval10000/ascending` | 295 µs | 297 µs | 1.87 ms | 5.14 ms |  | 289 µs | **1x slower** |
| `linear/3d/max200/eval10000/flipped` | 298 µs | 295 µs | 1.9 ms |  |  |  | 6.4x faster |
| `linear/4d/max10/eval10/ascending` | 377 ns | 376 ns | 95.5 µs |  |  |  | 2.5e+02x faster |
| `linear/4d/max10/eval10/flipped` | 377 ns | 411 ns | 97.6 µs |  |  |  | 2.6e+02x faster |
| `linear/4d/max10/eval100/ascending` | 4 µs | 4.01 µs | 107 µs |  |  |  | 27x faster |
| `linear/4d/max10/eval100/flipped` | 4.02 µs | 4.01 µs | 115 µs |  |  |  | 29x faster |
| `linear/4d/max10/eval1000/ascending` | 37 µs | 40.5 µs | 206 µs |  |  |  | 5.6x faster |
| `linear/4d/max10/eval1000/flipped` | 40.5 µs | 37 µs | 211 µs |  |  |  | 5.2x faster |
| `linear/4d/max10/eval10000/ascending` | 372 µs | 386 µs | 1.96 ms |  |  |  | 5.3x faster |
| `linear/4d/max10/eval10000/flipped` | 407 µs | 408 µs | 1.95 ms |  |  |  | 4.8x faster |
| `linear/4d/max40/eval10/ascending` | 421 ns | 460 ns | 96 µs | 5.98 µs |  | 298 ns | **1.4x slower** |
| `linear/4d/max40/eval10/flipped` | 421 ns | 462 ns | 97.2 µs |  |  |  | 2.3e+02x faster |
| `linear/4d/max40/eval100/ascending` | 4.14 µs | 4.53 µs | 117 µs | 62.8 µs |  | 2.76 µs | **1.5x slower** |
| `linear/4d/max40/eval100/flipped` | 4.15 µs | 4.54 µs | 108 µs |  |  |  | 26x faster |
| `linear/4d/max40/eval1000/ascending` | 41.5 µs | 45.4 µs | 241 µs | 571 µs |  | 29.6 µs | **1.4x slower** |
| `linear/4d/max40/eval1000/flipped` | 41.4 µs | 45.3 µs | 258 µs |  |  |  | 6.2x faster |
| `linear/4d/max40/eval10000/ascending` | 415 µs | 455 µs | 2.45 ms | 5.86 ms |  | 296 µs | **1.4x slower** |
| `linear/4d/max40/eval10000/flipped` | 414 µs | 454 µs | 2.28 ms |  |  |  | 5.5x faster |
| `linear/4d/max200/eval10/ascending` | 506 ns | 516 ns | 95.5 µs | 8.31 µs |  | 390 ns | **1.3x slower** |
| `linear/4d/max200/eval10/flipped` | 505 ns | 510 ns | 104 µs |  |  |  | 2.1e+02x faster |
| `linear/4d/max200/eval100/ascending` | 5.41 µs | 4.95 µs | 110 µs | 72.3 µs |  | 3.95 µs | **1.4x slower** |
| `linear/4d/max200/eval100/flipped` | 4.96 µs | 4.99 µs | 118 µs |  |  |  | 24x faster |
| `linear/4d/max200/eval1000/ascending` | 54.2 µs | 49.5 µs | 323 µs | 721 µs |  | 39.6 µs | **1.4x slower** |
| `linear/4d/max200/eval1000/flipped` | 49.5 µs | 49.7 µs | 302 µs |  |  |  | 6.1x faster |
| `linear/4d/max200/eval10000/ascending` | 543 µs | 546 µs | 3.06 ms | 7.8 ms |  | 398 µs | **1.4x slower** |
| `linear/4d/max200/eval10000/flipped` | 517 µs | 511 µs | 3.08 ms |  |  |  | 6x faster |
| `nearest/2d/each40/eval10000/ascending` | 151 µs | 165 µs | 701 µs | 1.61 ms |  | 136 µs | **1.1x slower** |
| `nearest/3d/each40/eval10000/ascending` | 212 µs | 233 µs | 1.02 ms | 2 ms |  | 228 µs | 1.1x faster |
| `pchip/2d/each40/eval10000/ascending` | 513 µs | 532 µs | 1.19 s |  |  |  | 2.3e+03x faster |
| `pchip/3d/each40/eval10000/ascending` | 1.5 ms | 1.6 ms | 3.3 s |  |  |  | 2.2e+03x faster |
| `cubic/2d/each40/eval10000/ascending` | 329 µs | 300 µs | 1.32 ms |  |  | 448 µs | 1.4x faster |
| `cubic/3d/each40/eval10000/ascending` | 663 µs | 726 µs | 3.21 ms |  |  | 852 µs | 1.3x faster |
| `cubic/2d/each10/eval1000/ascending` | 26.6 µs | 29.2 µs | 136 µs |  |  | 37.2 µs | 1.4x faster |
| `cubic/2d/each40/eval1000/ascending` | 30.3 µs | 32.9 µs | 146 µs |  |  | 41.1 µs | 1.4x faster |
| `cubic/2d/each100/eval1000/ascending` | 35 µs | 32.9 µs | 181 µs |  |  | 47.1 µs | 1.3x faster |
| `cubic/2d/each200/eval1000/ascending` | 34.6 µs | 34.6 µs | 219 µs |  |  | 50.1 µs | 1.4x faster |
| `cubic/2d/each400/eval1000/ascending` | 38.8 µs | 38.7 µs | 299 µs |  |  | 49.6 µs | 1.3x faster |
| `quintic/2d/each10/eval1000/ascending` | 48.5 µs | 48.4 µs | 229 µs |  |  |  | 4.7x faster |
| `quintic/2d/each40/eval1000/ascending` | 48 µs | 52.6 µs | 248 µs |  |  |  | 5.2x faster |
| `linear/2d/each400/eval10000/ascending` | 239 µs | 241 µs | 970 µs | 4.15 ms | 559 µs | 241 µs | 1x faster |
| `linear/2d/each800/eval10000/ascending` | 263 µs | 263 µs | 1.17 ms | 8.02 ms | 576 µs | 265 µs | 1x faster |

### Construction

| case | interpx | interpx, `from_owned` | SciPy | scirs2-interpolate | ndarray-interp | interpn | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/2d/max10/eval10/ascending` | 135 ns | 133 ns | 5.9 µs | 134 ns | 85.2 ns | 0 | **1.7e+02x slower** |
| `linear/2d/max10/eval10/flipped` | 167 ns | 171 ns | 12.2 µs |  |  |  | 73x faster |
| `linear/2d/max10/eval100/ascending` | 123 ns | 133 ns | 6.35 µs | 138 ns | 85.3 ns | 0 | **1.8e+02x slower** |
| `linear/2d/max10/eval100/flipped` | 168 ns | 168 ns | 12.8 µs |  |  |  | 76x faster |
| `linear/2d/max10/eval1000/ascending` | 139 ns | 140 ns | 6.37 µs | 144 ns | 85.4 ns | 0 | **1.8e+02x slower** |
| `linear/2d/max10/eval1000/flipped` | 168 ns | 162 ns | 12.6 µs |  |  |  | 75x faster |
| `linear/2d/max10/eval10000/ascending` | 134 ns | 133 ns | 6.34 µs | 134 ns | 85.3 ns | 0 | **2e+02x slower** |
| `linear/2d/max10/eval10000/flipped` | 169 ns | 183 ns | 11.7 µs |  |  |  | 69x faster |
| `linear/2d/max40/eval10/ascending` | 388 ns | 825 ns | 5.81 µs | 215 ns | 164 ns | 0 | **5e+02x slower** |
| `linear/2d/max40/eval10/flipped` | 665 ns | 1.01 µs | 12.7 µs |  |  |  | 19x faster |
| `linear/2d/max40/eval100/ascending` | 390 ns | 725 ns | 6.39 µs | 229 ns | 165 ns | 0 | **5e+02x slower** |
| `linear/2d/max40/eval100/flipped` | 693 ns | 1.07 µs | 12.7 µs |  |  |  | 18x faster |
| `linear/2d/max40/eval1000/ascending` | 387 ns | 876 ns | 5.83 µs | 233 ns | 150 ns | 0 | **5.7e+02x slower** |
| `linear/2d/max40/eval1000/flipped` | 703 ns | 1.18 µs | 11.8 µs |  |  |  | 17x faster |
| `linear/2d/max40/eval10000/ascending` | 385 ns | 704 ns | 5.81 µs | 229 ns | 164 ns | 0 | **4.9e+02x slower** |
| `linear/2d/max40/eval10000/flipped` | 700 ns | 1.21 µs | 12.2 µs |  |  |  | 17x faster |
| `linear/2d/max200/eval10/ascending` | 7.47 µs | 18.1 µs | 6.31 µs | 3.25 µs | 477 ns | 0 | **9.6e+03x slower** |
| `linear/2d/max200/eval10/flipped` | 15.9 µs | 26.1 µs | 12.3 µs |  |  |  | **1.3x slower** |
| `linear/2d/max200/eval100/ascending` | 7.5 µs | 19.7 µs | 5.92 µs | 3.2 µs | 516 ns | 0 | **1.1e+04x slower** |
| `linear/2d/max200/eval100/flipped` | 15.8 µs | 24.5 µs | 12.1 µs |  |  |  | **1.3x slower** |
| `linear/2d/max200/eval1000/ascending` | 7.55 µs | 17.8 µs | 6.27 µs | 3.06 µs | 476 ns | 0 | **9.7e+03x slower** |
| `linear/2d/max200/eval1000/flipped` | 15.7 µs | 27.7 µs | 13 µs |  |  |  | **1.2x slower** |
| `linear/2d/max200/eval10000/ascending` | 49.8 µs | 18.6 µs | 5.82 µs | 59.7 µs | 480 ns | 0 | **7.3e+04x slower** |
| `linear/2d/max200/eval10000/flipped` | 54.7 µs | 29.4 µs | 11.9 µs |  |  |  | **4.6x slower** |
| `linear/3d/max10/eval10/ascending` | 195 ns | 235 ns | 7.58 µs | 179 ns |  | 0 | **2e+02x slower** |
| `linear/3d/max10/eval10/flipped` | 322 ns | 346 ns | 17.4 µs |  |  |  | 54x faster |
| `linear/3d/max10/eval100/ascending` | 194 ns | 208 ns | 7.66 µs | 180 ns |  | 0 | **2.2e+02x slower** |
| `linear/3d/max10/eval100/flipped` | 322 ns | 334 ns | 15.9 µs |  |  |  | 49x faster |
| `linear/3d/max10/eval1000/ascending` | 195 ns | 204 ns | 7.59 µs | 166 ns |  | 0 | **2e+02x slower** |
| `linear/3d/max10/eval1000/flipped` | 324 ns | 347 ns | 15.7 µs |  |  |  | 49x faster |
| `linear/3d/max10/eval10000/ascending` | 195 ns | 241 ns | 8.44 µs | 181 ns |  | 0 | **2e+02x slower** |
| `linear/3d/max10/eval10000/flipped` | 299 ns | 337 ns | 17 µs |  |  |  | 57x faster |
| `linear/3d/max40/eval10/ascending` | 2.96 µs | 6.87 µs | 8.35 µs | 3.99 µs |  | 0 | **3e+03x slower** |
| `linear/3d/max40/eval10/flipped` | 8.1 µs | 12.1 µs | 17 µs |  |  |  | 2.1x faster |
| `linear/3d/max40/eval100/ascending` | 2.98 µs | 6.88 µs | 7.67 µs | 9.04 µs |  | 0 | **3.4e+03x slower** |
| `linear/3d/max40/eval100/flipped` | 8.14 µs | 13.5 µs | 17.1 µs |  |  |  | 2.1x faster |
| `linear/3d/max40/eval1000/ascending` | 23.7 µs | 6.29 µs | 8.35 µs | 15.4 µs |  | 0 | **2.4e+04x slower** |
| `linear/3d/max40/eval1000/flipped` | 29.4 µs | 11.6 µs | 16 µs |  |  |  | **1.8x slower** |
| `linear/3d/max40/eval10000/ascending` | 17.5 µs | 6.07 µs | 7.61 µs | 15.4 µs |  | 0 | **1.8e+04x slower** |
| `linear/3d/max40/eval10000/flipped` | 31 µs | 11.7 µs | 16 µs |  |  |  | **1.9x slower** |
| `linear/3d/max200/eval10/ascending` | 1.86 ms | 738 µs | 7.6 µs | 1.52 ms |  | 0 | **1.9e+06x slower** |
| `linear/3d/max200/eval10/flipped` | 3.7 ms | 2.06 ms | 17 µs |  |  |  | **2.2e+02x slower** |
| `linear/3d/max200/eval100/ascending` | 1.88 ms | 766 µs | 7.63 µs | 1.36 ms |  | 0 | **1.9e+06x slower** |
| `linear/3d/max200/eval100/flipped` | 3.53 ms | 2.15 ms | 17.4 µs |  |  |  | **2e+02x slower** |
| `linear/3d/max200/eval1000/ascending` | 2.17 ms | 781 µs | 7.67 µs | 1.46 ms |  | 0 | **2.2e+06x slower** |
| `linear/3d/max200/eval1000/flipped` | 3.59 ms | 2.16 ms | 16.2 µs |  |  |  | **2.2e+02x slower** |
| `linear/3d/max200/eval10000/ascending` | 1.99 ms | 727 µs | 7.63 µs | 1.39 ms |  | 0 | **2e+06x slower** |
| `linear/3d/max200/eval10000/flipped` | 3.65 ms | 2.12 ms | 17.5 µs |  |  |  | **2.1e+02x slower** |
| `linear/4d/max10/eval10/ascending` | 232 ns | 229 ns | 10.4 µs |  |  |  | 45x faster |
| `linear/4d/max10/eval10/flipped` | 373 ns | 384 ns | 18 µs |  |  |  | 48x faster |
| `linear/4d/max10/eval100/ascending` | 211 ns | 232 ns | 9.93 µs |  |  |  | 47x faster |
| `linear/4d/max10/eval100/flipped` | 374 ns | 375 ns | 19.4 µs |  |  |  | 52x faster |
| `linear/4d/max10/eval1000/ascending` | 232 ns | 247 ns | 10 µs |  |  |  | 43x faster |
| `linear/4d/max10/eval1000/flipped` | 376 ns | 403 ns | 19.8 µs |  |  |  | 53x faster |
| `linear/4d/max10/eval10000/ascending` | 229 ns | 266 ns | 9.98 µs |  |  |  | 44x faster |
| `linear/4d/max10/eval10000/flipped` | 364 ns | 394 ns | 19.8 µs |  |  |  | 54x faster |
| `linear/4d/max40/eval10/ascending` | 14.3 µs | 38.4 µs | 9.43 µs | 5.44 µs |  | 1.32 ns | **1.1e+04x slower** |
| `linear/4d/max40/eval10/flipped` | 56.3 µs | 69.6 µs | 21.7 µs |  |  |  | **2.6x slower** |
| `linear/4d/max40/eval100/ascending` | 14.3 µs | 42 µs | 10.2 µs | 5.71 µs |  | 1.37 ns | **1e+04x slower** |
| `linear/4d/max40/eval100/flipped` | 56.4 µs | 70.7 µs | 21.3 µs |  |  |  | **2.6x slower** |
| `linear/4d/max40/eval1000/ascending` | 14.4 µs | 35.2 µs | 9.36 µs | 5.89 µs |  | 1.27 ns | **1.1e+04x slower** |
| `linear/4d/max40/eval1000/flipped` | 52.3 µs | 68.9 µs | 19.8 µs |  |  |  | **2.6x slower** |
| `linear/4d/max40/eval10000/ascending` | 14.8 µs | 36.6 µs | 9.42 µs | 5.49 µs |  | 1.37 ns | **1.1e+04x slower** |
| `linear/4d/max40/eval10000/flipped` | 54.1 µs | 75.6 µs | 20.3 µs |  |  |  | **2.7x slower** |
| `linear/4d/max200/eval10/ascending` | 55.5 ms | 19.5 ms | 9.49 µs | 44.6 ms |  | 1.27 ns | **4.4e+07x slower** |
| `linear/4d/max200/eval10/flipped` | 123 ms | 82.1 ms | 21.8 µs |  |  |  | **5.6e+03x slower** |
| `linear/4d/max200/eval100/ascending` | 55.7 ms | 17.7 ms | 9.6 µs | 44.5 ms |  | 1.37 ns | **4.1e+07x slower** |
| `linear/4d/max200/eval100/flipped` | 123 ms | 90.1 ms | 20.2 µs |  |  |  | **6.1e+03x slower** |
| `linear/4d/max200/eval1000/ascending` | 56.2 ms | 17.7 ms | 10.4 µs | 45.2 ms |  | 1.27 ns | **4.4e+07x slower** |
| `linear/4d/max200/eval1000/flipped` | 121 ms | 83.5 ms | 20.1 µs |  |  |  | **6e+03x slower** |
| `linear/4d/max200/eval10000/ascending` | 59.4 ms | 18 ms | 10.2 µs | 45.3 ms |  | 1.27 ns | **4.7e+07x slower** |
| `linear/4d/max200/eval10000/flipped` | 127 ms | 86.8 ms | 21.9 µs |  |  |  | **5.8e+03x slower** |
| `nearest/2d/each40/eval10000/ascending` | 644 ns | 1.52 µs | 6.33 µs | 289 ns |  | 0 | **8.2e+02x slower** |
| `nearest/3d/each40/eval10000/ascending` | 23.1 µs | 55.4 µs | 8.17 µs | 9.35 µs |  | 0 | **2.4e+04x slower** |
| `pchip/2d/each40/eval10000/ascending` | 6.95 µs | 7.13 µs | 6.75 µs |  |  |  | **1x slower** |
| `pchip/3d/each40/eval10000/ascending` | 510 µs | 502 µs | 8.9 µs |  |  |  | **57x slower** |
| `cubic/2d/each40/eval10000/ascending` | 9.55 µs | 9.32 µs | 1.81 ms |  |  | 0 | **1.1e+04x slower** |
| `cubic/3d/each40/eval10000/ascending` | 391 µs | 341 µs | 171 ms |  |  | 0 | **4e+05x slower** |
| `cubic/2d/each10/eval1000/ascending` | 1.64 µs | 1.73 µs | 984 µs |  |  | 0 | **1.7e+03x slower** |
| `cubic/2d/each40/eval1000/ascending` | 9.58 µs | 9.41 µs | 1.96 ms |  |  | 0 | **1.2e+04x slower** |
| `cubic/2d/each100/eval1000/ascending` | 87.3 µs | 86.2 µs | 7.92 ms |  |  | 0 | **9.9e+04x slower** |
| `cubic/2d/each200/eval1000/ascending` | 204 µs | 200 µs | 30.8 ms |  |  | 0 | **2.3e+05x slower** |
| `cubic/2d/each400/eval1000/ascending` | 1.54 ms | 1.31 ms | 162 ms |  |  | 0 | **1.8e+06x slower** |
| `quintic/2d/each10/eval1000/ascending` | 2.01 µs | 2.18 µs | 3.74 ms |  |  |  | 1.9e+03x faster |
| `quintic/2d/each40/eval1000/ascending` | 13.8 µs | 13.3 µs | 8.22 ms |  |  |  | 5.9e+02x faster |
| `linear/2d/each400/eval10000/ascending` | 376 µs | 129 µs | 6.15 µs | 321 µs | 1.16 µs | 0 | **4.8e+05x slower** |
| `linear/2d/each800/eval10000/ascending` | 1.37 ms | 611 µs | 6.95 µs | 899 µs | 2.29 µs | 0 | **1.7e+06x slower** |

## `NearestNDInterpolator`

`n` points of the unit square, evaluated at the pixels of a square image.

### Evaluation

| case | interpx | interpx, `from_owned` | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `n10/image50` | 113 µs | 125 µs | 473 µs | 4.2x faster |
| `n10/image100` | 458 µs | 495 µs | 1.86 ms | 4.1x faster |
| `n10/image500` | 12.5 ms | 12.1 ms | 46.7 ms | 3.7x faster |
| `n100/image50` | 211 µs | 186 µs | 593 µs | 2.8x faster |
| `n100/image100` | 802 µs | 766 µs | 2.37 ms | 3x faster |
| `n100/image500` | 18.8 ms | 18.1 ms | 58.4 ms | 3.1x faster |
| `n1000/image50` | 340 µs | 331 µs | 859 µs | 2.5x faster |
| `n1000/image100` | 1.26 ms | 1.18 ms | 3.2 ms | 2.5x faster |
| `n1000/image500` | 25.9 ms | 23.7 ms | 71.8 ms | 2.8x faster |
| `n10000/image50` | 398 µs | 325 µs | 1.08 ms | 2.7x faster |
| `n10000/image100` | 1.5 ms | 1.34 ms | 4.13 ms | 2.7x faster |
| `n10000/image500` | 35 ms | 32.6 ms | 90.1 ms | 2.6x faster |
| `n100000/image50` | 514 µs | 472 µs | 1.48 ms | 2.9x faster |
| `n100000/image100` | 1.99 ms | 1.81 ms | 5.68 ms | 2.9x faster |
| `n100000/image500` | 47.2 ms | 45.5 ms | 130 ms | 2.8x faster |

### Construction

| case | interpx | interpx, `from_owned` | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `n10/image50` | 552 ns | 561 ns | 10.7 µs | 19x faster |
| `n10/image100` | 565 ns | 570 ns | 10.8 µs | 19x faster |
| `n10/image500` | 537 ns | 566 ns | 10.6 µs | 20x faster |
| `n100/image50` | 3.31 µs | 3.71 µs | 15.2 µs | 4.6x faster |
| `n100/image100` | 3.33 µs | 3.65 µs | 15.6 µs | 4.7x faster |
| `n100/image500` | 3.41 µs | 3.65 µs | 15.7 µs | 4.6x faster |
| `n1000/image50` | 34.5 µs | 53.2 µs | 127 µs | 3.7x faster |
| `n1000/image100` | 36.6 µs | 43.2 µs | 125 µs | 3.4x faster |
| `n1000/image500` | 37 µs | 53.9 µs | 123 µs | 3.3x faster |
| `n10000/image50` | 703 µs | 653 µs | 1.8 ms | 2.6x faster |
| `n10000/image100` | 663 µs | 727 µs | 1.81 ms | 2.7x faster |
| `n10000/image500` | 654 µs | 663 µs | 1.79 ms | 2.7x faster |
| `n100000/image50` | 8.53 ms | 8.66 ms | 26.7 ms | 3.1x faster |
| `n100000/image100` | 8.75 ms | 8.62 ms | 26.2 ms | 3x faster |
| `n100000/image500` | 9.23 ms | 8.48 ms | 26.5 ms | 2.9x faster |

## `RBFInterpolator`

`n` points of a square, evaluated at 1000 points in no order (`scattered`) or at the 10 000 pixels of an image (`image`). `global`: all points; `k50`: the 50 nearest.

### Evaluation

| case | interpx | interpx, `from_owned` | SciPy | scirs2-interpolate | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: |
| `linear/global/n10/scattered` | 24.1 µs | 26.8 µs | 61.6 µs | 182 µs | 2.6x faster |
| `linear/k50/n10/scattered` | 24.2 µs | 27 µs | 1.38 ms |  | 57x faster |
| `thin_plate_spline/global/n10/scattered` | 92 µs | 97.3 µs | 116 µs | 303 µs | 1.3x faster |
| `thin_plate_spline/k50/n10/scattered` | 89.8 µs | 97.5 µs | 1.52 ms |  | 17x faster |
| `cubic/global/n10/scattered` | 40.6 µs | 40.8 µs | 59.4 µs | 198 µs | 1.5x faster |
| `cubic/k50/n10/scattered` | 40.9 µs | 43.3 µs | 1.37 ms |  | 33x faster |
| `quintic/global/n10/scattered` | 75 µs | 78.5 µs | 66.8 µs | 216 µs | **1.1x slower** |
| `quintic/k50/n10/scattered` | 72.5 µs | 78.3 µs | 1.54 ms |  | 21x faster |
| `multiquadric/global/n10/scattered` | 40.6 µs | 43.4 µs | 53.3 µs | 236 µs | 1.3x faster |
| `multiquadric/k50/n10/scattered` | 40.6 µs | 43.6 µs | 1.43 ms |  | 35x faster |
| `inverse_multiquadric/global/n10/scattered` | 47.6 µs | 50.8 µs | 56.9 µs | 291 µs | 1.2x faster |
| `inverse_multiquadric/k50/n10/scattered` | 51.2 µs | 50.7 µs | 1.46 ms |  | 28x faster |
| `inverse_quadratic/global/n10/scattered` | 33.1 µs | 33.9 µs | 47.4 µs |  | 1.4x faster |
| `inverse_quadratic/k50/n10/scattered` | 33.8 µs | 34.2 µs | 1.5 ms |  | 44x faster |
| `gaussian/global/n10/scattered` | 79.6 µs | 85.4 µs | 88.8 µs | 358 µs | 1.1x faster |
| `gaussian/k50/n10/scattered` | 85 µs | 85.2 µs | 1.51 ms |  | 18x faster |
| `linear/global/n100/scattered` | 153 µs | 152 µs | 292 µs | 1.8 ms | 1.9x faster |
| `linear/k50/n100/scattered` | 14.7 ms | 14.5 ms | 21.4 ms |  | 1.5x faster |
| `thin_plate_spline/global/n100/scattered` | 748 µs | 743 µs | 900 µs | 3.15 ms | 1.2x faster |
| `thin_plate_spline/k50/n100/scattered` | 25.8 ms | 23.9 ms | 26.5 ms |  | 1x faster |
| `cubic/global/n100/scattered` | 163 µs | 178 µs | 302 µs | 2.06 ms | 1.9x faster |
| `cubic/k50/n100/scattered` | 16.7 ms | 16.7 ms | 21.4 ms |  | 1.3x faster |
| `quintic/global/n100/scattered` | 232 µs | 215 µs | 314 µs | 2.16 ms | 1.4x faster |
| `quintic/k50/n100/scattered` | 20.4 ms | 19.5 ms | 23.9 ms |  | 1.2x faster |
| `multiquadric/global/n100/scattered` | 291 µs | 284 µs | 311 µs | 2.38 ms | 1.1x faster |
| `multiquadric/k50/n100/scattered` | 16.8 ms | 16.3 ms | 21.4 ms |  | 1.3x faster |
| `inverse_multiquadric/global/n100/scattered` | 381 µs | 379 µs | 398 µs | 2.77 ms | 1x faster |
| `inverse_multiquadric/k50/n100/scattered` | 17.8 ms | 16.9 ms | 21.6 ms |  | 1.2x faster |
| `inverse_quadratic/global/n100/scattered` | 239 µs | 261 µs | 288 µs |  | 1.2x faster |
| `inverse_quadratic/k50/n100/scattered` | 15.7 ms | 15.5 ms | 20.2 ms |  | 1.3x faster |
| `gaussian/global/n100/scattered` | 634 µs | 636 µs | 732 µs | 3.71 ms | 1.2x faster |
| `gaussian/k50/n100/scattered` | 21.8 ms | 22.4 ms | 23.9 ms |  | 1.1x faster |
| `linear/global/n1000/scattered` | 1.51 ms | 1.39 ms | 3.9 ms | 18.7 ms | 2.6x faster |
| `linear/k50/n1000/scattered` | 15.5 ms | 15.7 ms | 36.7 ms |  | 2.4x faster |
| `thin_plate_spline/global/n1000/scattered` | 7.95 ms | 7.35 ms | 10.6 ms | 31.9 ms | 1.3x faster |
| `thin_plate_spline/k50/n1000/scattered` | 24.7 ms | 25.8 ms | 46 ms |  | 1.9x faster |
| `cubic/global/n1000/scattered` | 1.5 ms | 1.38 ms | 3.99 ms | 20.1 ms | 2.7x faster |
| `cubic/k50/n1000/scattered` | 17.7 ms | 17.9 ms | 39.6 ms |  | 2.2x faster |
| `quintic/global/n1000/scattered` | 1.7 ms | 1.72 ms | 3.94 ms | 22.6 ms | 2.3x faster |
| `quintic/k50/n1000/scattered` | 21 ms | 21 ms | 44.3 ms |  | 2.1x faster |
| `multiquadric/global/n1000/scattered` | 2.83 ms | 2.91 ms | 4.25 ms | 25.2 ms | 1.5x faster |
| `multiquadric/k50/n1000/scattered` | 18.2 ms | 17.5 ms | 38.2 ms |  | 2.1x faster |
| `inverse_multiquadric/global/n1000/scattered` | 3.9 ms | 4.01 ms | 5.45 ms | 28.7 ms | 1.4x faster |
| `inverse_multiquadric/k50/n1000/scattered` | 19.7 ms | 18.8 ms | 40.4 ms |  | 2x faster |
| `inverse_quadratic/global/n1000/scattered` | 2.47 ms | 2.47 ms | 3.94 ms |  | 1.6x faster |
| `inverse_quadratic/k50/n1000/scattered` | 17.4 ms | 16.8 ms | 37.8 ms |  | 2.2x faster |
| `gaussian/global/n1000/scattered` | 6.73 ms | 6.46 ms | 8.58 ms | 37.9 ms | 1.3x faster |
| `gaussian/k50/n1000/scattered` | 24.7 ms | 23.8 ms | 45.4 ms |  | 1.8x faster |
| `thin_plate_spline/global/n3000/scattered` | 23.1 ms | 23.7 ms | 31.8 ms | 97.6 ms | 1.4x faster |
| `thin_plate_spline/k50/n10000/scattered` | 25.2 ms | 25.5 ms | 47.5 ms |  | 1.9x faster |
| `thin_plate_spline/k50/n10000/image` | 250 ms | 245 ms | 491 ms |  | 2x faster |
| `thin_plate_spline/k50/n100000/scattered` | 26.3 ms | 26.5 ms | 50 ms |  | 1.9x faster |
| `thin_plate_spline/k50/n100000/image` | 262 ms | 253 ms | 484 ms |  | 1.8x faster |

### Construction

| case | interpx | interpx, `from_owned` | SciPy | scirs2-interpolate | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: |
| `linear/global/n10/scattered` | 1.29 µs | 1.38 µs | 12.6 µs | 5.58 µs | 4.3x faster |
| `linear/k50/n10/scattered` | 1.27 µs | 1.41 µs | 17 µs |  | 13x faster |
| `thin_plate_spline/global/n10/scattered` | 1.89 µs | 2.04 µs | 13.7 µs | 6.66 µs | 3.5x faster |
| `thin_plate_spline/k50/n10/scattered` | 1.89 µs | 2.03 µs | 17.1 µs |  | 9x faster |
| `cubic/global/n10/scattered` | 1.66 µs | 1.78 µs | 13.6 µs | 6.11 µs | 3.7x faster |
| `cubic/k50/n10/scattered` | 1.66 µs | 1.8 µs | 18 µs |  | 11x faster |
| `quintic/global/n10/scattered` | 2.38 µs | 2.54 µs | 15.6 µs | 5.83 µs | 2.5x faster |
| `quintic/k50/n10/scattered` | 2.46 µs | 2.57 µs | 17.9 µs |  | 7.3x faster |
| `multiquadric/global/n10/scattered` | 1.36 µs | 1.42 µs | 13.2 µs | 6.64 µs | 4.9x faster |
| `multiquadric/k50/n10/scattered` | 1.35 µs | 1.46 µs | 17.1 µs |  | 13x faster |
| `inverse_multiquadric/global/n10/scattered` | 1.41 µs | 1.48 µs | 12.1 µs | 7.02 µs | 5x faster |
| `inverse_multiquadric/k50/n10/scattered` | 1.38 µs | 1.47 µs | 17.1 µs |  | 12x faster |
| `inverse_quadratic/global/n10/scattered` | 1.39 µs | 1.41 µs | 12.6 µs |  | 9x faster |
| `inverse_quadratic/k50/n10/scattered` | 1.32 µs | 1.31 µs | 17.4 µs |  | 13x faster |
| `gaussian/global/n10/scattered` | 1.48 µs | 1.61 µs | 12.3 µs | 7.86 µs | 5.3x faster |
| `gaussian/k50/n10/scattered` | 1.52 µs | 1.61 µs | 17 µs |  | 11x faster |
| `linear/global/n100/scattered` | 44.3 µs | 58.7 µs | 82.7 µs | 802 µs | 1.9x faster |
| `linear/k50/n100/scattered` | 3.33 µs | 3.41 µs | 23.6 µs |  | 7.1x faster |
| `thin_plate_spline/global/n100/scattered` | 78.2 µs | 114 µs | 117 µs | 916 µs | 1.5x faster |
| `thin_plate_spline/k50/n100/scattered` | 3.77 µs | 4.04 µs | 24.1 µs |  | 6.4x faster |
| `cubic/global/n100/scattered` | 51.6 µs | 71.5 µs | 83.2 µs | 764 µs | 1.6x faster |
| `cubic/k50/n100/scattered` | 3.33 µs | 3.46 µs | 24.4 µs |  | 7.3x faster |
| `quintic/global/n100/scattered` | 55.6 µs | 94.7 µs | 101 µs | 819 µs | 1.8x faster |
| `quintic/k50/n100/scattered` | 3.29 µs | 3.56 µs | 25.8 µs |  | 7.8x faster |
| `multiquadric/global/n100/scattered` | 52.4 µs | 67.4 µs | 84.5 µs | 664 µs | 1.6x faster |
| `multiquadric/k50/n100/scattered` | 3.41 µs | 3.19 µs | 23 µs |  | 6.7x faster |
| `inverse_multiquadric/global/n100/scattered` | 56.1 µs | 77.7 µs | 83.3 µs | 888 µs | 1.5x faster |
| `inverse_multiquadric/k50/n100/scattered` | 3.24 µs | 3.39 µs | 23.3 µs |  | 7.2x faster |
| `inverse_quadratic/global/n100/scattered` | 48.3 µs | 56.8 µs | 75.6 µs |  | 1.6x faster |
| `inverse_quadratic/k50/n100/scattered` | 3.31 µs | 3.43 µs | 22.2 µs |  | 6.7x faster |
| `gaussian/global/n100/scattered` | 67.4 µs | 100 µs | 97.9 µs | 828 µs | 1.5x faster |
| `gaussian/k50/n100/scattered` | 3.15 µs | 3.4 µs | 22.1 µs |  | 7x faster |
| `linear/global/n1000/scattered` | 17.6 ms | 17.9 ms | 62.9 ms | 197 ms | 3.6x faster |
| `linear/k50/n1000/scattered` | 37.7 µs | 59.6 µs | 135 µs |  | 3.6x faster |
| `thin_plate_spline/global/n1000/scattered` | 21.1 ms | 21.5 ms | 66.5 ms | 212 ms | 3.2x faster |
| `thin_plate_spline/k50/n1000/scattered` | 38.3 µs | 62.2 µs | 141 µs |  | 3.7x faster |
| `cubic/global/n1000/scattered` | 17.9 ms | 17.5 ms | 62.2 ms | 196 ms | 3.5x faster |
| `cubic/k50/n1000/scattered` | 38.1 µs | 61.4 µs | 139 µs |  | 3.6x faster |
| `quintic/global/n1000/scattered` | 18.9 ms | 18.2 ms | 61.7 ms | 198 ms | 3.3x faster |
| `quintic/k50/n1000/scattered` | 36.3 µs | 60.6 µs | 140 µs |  | 3.9x faster |
| `multiquadric/global/n1000/scattered` | 19.5 ms | 19.5 ms | 61.7 ms | 201 ms | 3.2x faster |
| `multiquadric/k50/n1000/scattered` | 37.8 µs | 60.2 µs | 131 µs |  | 3.5x faster |
| `inverse_multiquadric/global/n1000/scattered` | 18.9 ms | 19.1 ms | 65.4 ms | 210 ms | 3.5x faster |
| `inverse_multiquadric/k50/n1000/scattered` | 38 µs | 59.9 µs | 139 µs |  | 3.7x faster |
| `inverse_quadratic/global/n1000/scattered` | 18.8 ms | 18.4 ms | 60.2 ms |  | 3.2x faster |
| `inverse_quadratic/k50/n1000/scattered` | 38.4 µs | 58.2 µs | 132 µs |  | 3.4x faster |
| `gaussian/global/n1000/scattered` | 20.5 ms | 20.2 ms | 64 ms | 211 ms | 3.1x faster |
| `gaussian/k50/n1000/scattered` | 38.2 µs | 56.7 µs | 136 µs |  | 3.6x faster |
| `thin_plate_spline/global/n3000/scattered` | 423 ms | 461 ms | 1.92 s | 6.18 s | 4.5x faster |
| `thin_plate_spline/k50/n10000/scattered` | 720 µs | 734 µs | 1.71 ms |  | 2.4x faster |
| `thin_plate_spline/k50/n10000/image` | 765 µs | 736 µs | 1.73 ms |  | 2.3x faster |
| `thin_plate_spline/k50/n100000/scattered` | 8.99 ms | 8.99 ms | 26.1 ms |  | 2.9x faster |
| `thin_plate_spline/k50/n100000/image` | 9.44 ms | 8.89 ms | 26.3 ms |  | 2.8x faster |

## Memory

Megabytes. "Data read" is the process with the data of the case in memory, before anything is built. The columns after it are what was added to that: the highest while building, with the interpolator built, the highest while evaluating. "Highest" is the highest of the whole run, repeats included.

| case | library | data read | while building | built | while evaluating | highest |
| :--- | :--- | ---: | ---: | ---: | ---: | ---: |
| `nearest_nd/n100000/image500` | interpx | 12.0 | +3.5 | +2.1 | +4.1 | 18.2 |
| `nearest_nd/n100000/image500` | interpx, `from_owned` | 13.7 | +1.9 | +0.9 | +2.9 | 18.6 |
| `nearest_nd/n100000/image500` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |
| `rbf/thin_plate_spline/global/n3000/scattered` | interpx | 4.3 | +75.9 | +3.8 | +3.8 | 80.3 |
| `rbf/thin_plate_spline/global/n3000/scattered` | interpx, `from_owned` | 4.2 | +76.1 | +3.8 | +3.8 | 80.3 |
| `rbf/thin_plate_spline/global/n3000/scattered` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |
| `rbf/thin_plate_spline/global/n3000/scattered` | scirs2-interpolate | 2.8 | +289.8 | +1.7 | +1.7 | 292.5 |
| `rbf/thin_plate_spline/k50/n100000/image` | interpx | 7.2 | +3.0 | +4.2 | +4.5 | 11.9 |
| `rbf/thin_plate_spline/k50/n100000/image` | interpx, `from_owned` | 7.5 | +5.1 | +3.9 | +4.2 | 13.5 |
| `rbf/thin_plate_spline/k50/n100000/image` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |
| `oned/cubic_spline/n10000/sorted` | interpx | 4.8 | +0.4 | +0.4 | +0.4 | 5.9 |
| `oned/cubic_spline/n10000/sorted` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |
| `oned/cubic_spline/n10000/sorted` | scirs2-interpolate | 4.2 | +0.9 | +1.0 | +1.7 | 5.9 |
| `oned/cubic_spline/n10000/sorted` | ndarray-interp | 3.2 | +1.7 | +1.7 | +1.7 | 5.4 |
| `grid/linear/4d/max200/eval10000/ascending` | interpx | 231.8 | +173.3 | +173.3 | +173.3 | 405.0 |
| `grid/linear/4d/max200/eval10000/ascending` | interpx, `from_owned` | 222.6 | +182.1 | +182.1 | +182.2 | 404.8 |
| `grid/linear/4d/max200/eval10000/ascending` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |
| `grid/linear/4d/max200/eval10000/ascending` | scirs2-interpolate | 225.2 | +172.8 | +179.7 | +179.7 | 404.9 |
| `grid/linear/4d/max200/eval10000/ascending` | interpn | 213.8 | +0.0 | +0.0 | +0.0 | 402.0 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpx | 4.3 | +4.0 | +4.0 | +4.0 | 8.4 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpx, `from_owned` | 6.5 | +1.0 | +1.8 | +1.8 | 8.5 |
| `grid/cubic/2d/each400/eval1000/ascending` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpn | 3.9 | +1.6 | +1.6 | +1.6 | 5.5 |
| `grid/linear/2d/each800/eval10000/ascending` | interpx | 9.4 | +5.5 | +5.5 | +5.5 | 14.9 |
| `grid/linear/2d/each800/eval10000/ascending` | interpx, `from_owned` | 11.1 | +0.0 | +3.7 | +3.7 | 14.9 |
| `grid/linear/2d/each800/eval10000/ascending` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |
| `grid/linear/2d/each800/eval10000/ascending` | scirs2-interpolate | 11.7 | +0.0 | +3.2 | +3.3 | 14.9 |
| `grid/linear/2d/each800/eval10000/ascending` | ndarray-interp | 9.5 | +0.1 | +0.1 | +0.2 | 9.8 |
| `grid/linear/2d/each800/eval10000/ascending` | interpn | 10.3 | +0.0 | +0.0 | +0.0 | 10.6 |

### Over time

The resident set size of each library's process over its run, with the ends of the phases. The time after `evaluated` is the timed repeats.

`nearest_nd/n100000/image500`

![memory over time, nearest_nd/n100000/image500](performance/nearest_nd-n100000-image500.svg)

`rbf/thin_plate_spline/global/n3000/scattered`

![memory over time, rbf/thin_plate_spline/global/n3000/scattered](performance/rbf-thin_plate_spline-global-n3000-scattered.svg)

`rbf/thin_plate_spline/k50/n100000/image`

![memory over time, rbf/thin_plate_spline/k50/n100000/image](performance/rbf-thin_plate_spline-k50-n100000-image.svg)

`oned/cubic_spline/n10000/sorted`

![memory over time, oned/cubic_spline/n10000/sorted](performance/oned-cubic_spline-n10000-sorted.svg)

`grid/linear/4d/max200/eval10000/ascending`

![memory over time, grid/linear/4d/max200/eval10000/ascending](performance/grid-linear-4d-max200-eval10000-ascending.svg)

`grid/cubic/2d/each400/eval1000/ascending`

![memory over time, grid/cubic/2d/each400/eval1000/ascending](performance/grid-cubic-2d-each400-eval1000-ascending.svg)

`grid/linear/2d/each800/eval10000/ascending`

![memory over time, grid/linear/2d/each800/eval10000/ascending](performance/grid-linear-2d-each800-eval10000-ascending.svg)

## Do the libraries compute the same thing

The largest distance of a library's values from SciPy's, over the cases of a family, as a share of the largest value. Under 1e-06 the two compute the same interpolant.

| family | interpx | interpx, `from_owned` | scirs2-interpolate | ndarray-interp | interpn |
| :--- | ---: | ---: | ---: | ---: | ---: |
| One variable | 3e-16 |  | 6e-02 | 1e-04 |  |
| `PPoly` | 9e-16 |  |  |  |  |
| `BSpline` | 4e-16 |  |  |  |  |
| `BSplineCurve` | 4e-16 |  |  |  |  |
| Regular grid | 4e-04 | 4e-04 | 4e-16 | 6e-16 | 7e-01 |
| `NearestNDInterpolator` | 0e+00 | 0e+00 |  |  |  |
| `RBFInterpolator` | 8e-05 | 8e-05 | 3e+01 |  |  |
