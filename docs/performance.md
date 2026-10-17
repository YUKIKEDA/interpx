# Performance: interpx against SciPy and other Rust crates

<!-- Written by compare/perf/report.py from a run of compare/perf/run.py. Do not edit by hand: run it again. -->

Measured on 2026-10-18, at commit `959c7c5d8c`.

## How it was measured

- Machine: Intel64 Family 6 Model 191 Stepping 2, GenuineIntel, 16 logical cores, 48 GB. Windows-11-10.0.26300-SP0.
- rustc 1.99.0 (b940084d7 2026-09-28); Python 3.13.1, NumPy 2.5.4, SciPy 1.18.1.
- The other Rust crates, at the versions in `compare/perf/runners/*/Cargo.toml`: `scirs2-interpolate` (the widest of the Rust crates: one variable, grids, RBF), `ndarray-interp` (the most used for one variable; also bilinear), `interpn` (grids).
- Every library ran each case in a process of its own, pinned to one core, on one thread, on the same data. The libraries took turns case by case, 3 rounds.
- `interpx, from_owned` is interpx built from data that is handed over (`from_owned`) and so not copied. Its column is there for the types that keep their data as it is given (the grid, `NearestNDInterpolator`, `RBFInterpolator`), next to the construction from borrowed data (`new`), which is the one the comparison counts. Each of its samples is given data made just before it, which is not in the nearest cache as the data of the other constructions is: on small grids it reads slower than `new` for that reason.
- `interpx, parallel` is interpx evaluated with `eval_batch_into_par` (the `rayon` feature) on 8 threads, not pinned to a core. Every other column is one thread on one core, and the comparison counts those. A parallel batch times its first points and is shared out only when the rest is worth it, so on small cases the column reads as the serial one does.
- A time is that of one call: construction, or the evaluation of all the queries of the case (`eval_batch_into` for interpx, the call with an array for the others). It is the median of the repeats of a process, and then the median of the rounds.
- A call of two seconds or more is timed once: the construction that made the interpolator, or the evaluation that wrote the values, is its measurement, and it is not made again. A library whose construction or evaluation of a case was timed so is run on that case in the first round only, and its time is then one measurement; a call near two seconds that was timed several times goes on to the later rounds. (The rule goes by the time, whatever the library.)
- interpx is called slower where its median is more than 1.10 times another library's and its fastest repeat is slower than that library's fastest.
- Memory is the resident set size of the process, sampled every millisecond from outside it.
- Command: `perf/run.py`, run in `compare/`.

## Where interpx is slower

| family | operation | cases | interpx fastest | slower than SciPy | slower than a Rust crate |
| :--- | :--- | ---: | ---: | ---: | ---: |
| One variable | evaluation | 112 | 112 | 0 | 0 |
| One variable | construction | 112 | 73 | 0 | 39 |
| `PPoly` | evaluation | 4 | 4 | 0 | 0 |
| `PPoly` | construction | 4 | 4 | 0 | 0 |
| `BSpline` | evaluation | 4 | 4 | 0 | 0 |
| `BSpline` | construction | 4 | 4 | 0 | 0 |
| `BSplineCurve` | evaluation | 2 | 2 | 0 | 0 |
| `BSplineCurve` | construction | 2 | 2 | 0 | 0 |
| Regular grid | evaluation | 92 | 74 | 0 | 18 |
| Regular grid | construction | 92 | 29 | 37 | 46 |
| `NearestNDInterpolator` | evaluation | 16 | 16 | 0 | 0 |
| `NearestNDInterpolator` | construction | 16 | 16 | 0 | 0 |
| `LinearNDInterpolator` | evaluation | 19 | 19 | 0 | 0 |
| `LinearNDInterpolator` | construction | 19 | 3 | 0 | 16 |
| `RBFInterpolator` | evaluation | 61 | 61 | 0 | 0 |
| `RBFInterpolator` | construction | 61 | 61 | 0 | 0 |

209 comparisons (a case, an operation, and a library), by method. "Same values" says whether that library's values are SciPy's up to rounding: where they are not, it computes another interpolant under the same name, and its time is not that of the same work.

| method | operation | against | same values | cases | slower by | the largest |
| :--- | :--- | :--- | :--- | ---: | ---: | :--- |
| `grid/cubic` | construction | interpn | no | 8 | 1.84e+03x to 8.97e+07x | `grid/cubic/2d/auckland-ba32-0401/heldout/ascending`: 78.8 ms against 0.88 ns |
| `grid/linear` | construction | interpn | yes | 35 | 164x to 3.86e+07x | `grid/linear/4d/max200/eval1000/ascending`: 52.8 ms against 1.37 ns |
| `grid/nearest` | construction | interpn | yes | 3 | 834x to 2.15e+07x | `grid/nearest/2d/auckland-ba32-0401/heldout/ascending`: 16.8 ms against 0.78 ns |
| `grid/linear` | construction | SciPy | yes | 34 | 1.15x to 5.2e+03x | `grid/linear/4d/max200/eval10/flipped`: 114 ms against 21.8 µs |
| `grid/nearest` | construction | SciPy | yes | 2 | 2.8x to 2.39e+03x | `grid/nearest/2d/auckland-ba32-0401/heldout/ascending`: 16.8 ms against 7.05 µs |
| `grid/linear` | construction | ndarray-interp | yes | 15 | 1.5x to 2.07e+03x | `grid/linear/2d/auckland-ba32-0401/heldout/ascending`: 17.2 ms against 8.34 µs |
| `grid/pchip` | construction | SciPy | yes | 1 | 54.8x | `grid/pchip/3d/each40/eval10000/ascending`: 518 µs against 9.46 µs |
| `linear_nd` | construction | spade | yes | 16 | 1.57x to 3.06x | `linear_nd/n10/image100`: 5.86 µs against 1.91 µs |
| `oned/pchip` | construction | scirs2-interpolate | yes | 12 | 1.3x to 2.96x | `oned/pchip/n10000/shuffled-uneven`: 138 µs against 46.8 µs |
| `oned/akima` | construction | scirs2-interpolate | yes | 8 | 1.48x to 2.94x | `oned/akima/n10000/sorted`: 143 µs against 48.7 µs |
| `grid/linear` | construction | scirs2-interpolate | yes | 27 | 1.12x to 2.54x | `grid/linear/4d/max40/eval10000/ascending`: 14.4 µs against 5.65 µs |
| `grid/nearest` | construction | scirs2-interpolate | yes | 3 | 1.23x to 2.38x | `grid/nearest/3d/each40/eval10000/ascending`: 23.1 µs against 9.69 µs |
| `oned/linear` | construction | ndarray-interp | yes | 10 | 1.36x to 2.1x | `oned/linear/n10/sorted`: 87.2 ns against 41.5 ns |
| `grid/linear` | evaluation | interpn | yes | 15 | 1.1x to 1.85x | `grid/linear/3d/max200/eval10/ascending`: 532 ns against 287 ns |
| `oned/linear` | construction | scirs2-interpolate | yes | 11 | 1.21x to 1.77x | `oned/linear/n10000/shuffled`: 22.9 µs against 12.9 µs |
| `oned/cubic_spline` | construction | scirs2-interpolate | no | 6 | 1.16x to 1.65x | `oned/cubic_spline/n10000/sorted`: 240 µs against 145 µs |
| `grid/linear` | evaluation | ndarray-interp | yes | 2 | 1.12x to 1.14x | `grid/linear/2d/auckland-ba32-0401/heldout/ascending`: 426 ms against 375 ms |
| `grid/nearest` | evaluation | interpn | yes | 1 | 1.11x | `grid/nearest/2d/each40/eval10000/ascending`: 166 µs against 150 µs |

<details><summary>Every comparison where interpx is slower</summary>

"Values apart" is the distance of that library's values from SciPy's, as a share of the largest value.

| case | operation | against | interpx | the other | slower by | values apart |
| :--- | :--- | :--- | ---: | ---: | ---: | ---: |
| `grid/cubic/2d/auckland-ba32-0401/heldout/ascending` | construction | interpn | 78.8 ms | 0.88 ns | 89651313.78x | 5e-03 |
| `grid/linear/4d/max200/eval1000/ascending` | construction | interpn | 52.8 ms | 1.37 ns | 38623158.86x | 3e-16 |
| `grid/linear/4d/max200/eval10/ascending` | construction | interpn | 49.7 ms | 1.37 ns | 36326765.71x | 2e-16 |
| `grid/linear/4d/max200/eval10000/ascending` | construction | interpn | 49.7 ms | 1.37 ns | 36317184.00x | 4e-16 |
| `grid/linear/4d/max200/eval100/ascending` | construction | interpn | 47.9 ms | 1.37 ns | 35000539.43x | 3e-16 |
| `grid/linear/2d/auckland-ba32-0401/heldout/ascending` | construction | interpn | 17.2 ms | 0.78 ns | 22054528.00x | 0e+00 |
| `grid/nearest/2d/auckland-ba32-0401/heldout/ascending` | construction | interpn | 16.8 ms | 0.78 ns | 21540608.00x | 0e+00 |
| `grid/linear/3d/max200/eval10/ascending` | construction | interpn | 2.36 ms | 0.98 ns | 2414899.20x | 3e-16 |
| `grid/linear/3d/max200/eval1000/ascending` | construction | interpn | 1.82 ms | 0.98 ns | 1862758.40x | 3e-16 |
| `grid/linear/3d/max200/eval100/ascending` | construction | interpn | 1.82 ms | 0.98 ns | 1860403.20x | 3e-16 |
| `grid/linear/3d/max200/eval10000/ascending` | construction | interpn | 1.8 ms | 0.98 ns | 1842073.60x | 5e-16 |
| `grid/cubic/2d/each400/eval1000/ascending` | construction | interpn | 1.53 ms | 0.98 ns | 1564876.80x | 5e-01 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | interpn | 1.21 ms | 0.78 ns | 1543424.00x | 3e-16 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | interpn | 333 µs | 0.78 ns | 426240.00x | 3e-16 |
| `grid/cubic/3d/each40/eval10000/ascending` | construction | interpn | 378 µs | 0.98 ns | 387072.00x | 5e-01 |
| `grid/cubic/2d/each200/eval1000/ascending` | construction | interpn | 197 µs | 0.98 ns | 201932.80x | 2e-01 |
| `grid/cubic/2d/each100/eval1000/ascending` | construction | interpn | 44.6 µs | 0.98 ns | 45708.80x | 2e-01 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | interpn | 23.1 µs | 0.98 ns | 23648.00x | 0e+00 |
| `grid/linear/2d/max200/eval10/ascending` | construction | interpn | 7.41 µs | 0.68 ns | 10838.86x | 2e-16 |
| `grid/cubic/2d/each40/eval1000/ascending` | construction | interpn | 9.42 µs | 0.88 ns | 10716.44x | 7e-01 |
| `grid/linear/4d/max40/eval10/ascending` | construction | interpn | 14.4 µs | 1.37 ns | 10514.29x | 2e-16 |
| `grid/linear/4d/max40/eval10000/ascending` | construction | interpn | 14.4 µs | 1.37 ns | 10507.43x | 5e-16 |
| `grid/linear/4d/max40/eval1000/ascending` | construction | interpn | 14.3 µs | 1.37 ns | 10468.57x | 6e-16 |
| `grid/linear/4d/max40/eval100/ascending` | construction | interpn | 14.3 µs | 1.37 ns | 10464.00x | 2e-16 |
| `grid/cubic/2d/each40/eval10000/ascending` | construction | interpn | 9.65 µs | 0.98 ns | 9881.60x | 7e-01 |
| `grid/linear/2d/max200/eval10000/ascending` | construction | interpn | 7.59 µs | 0.78 ns | 9720.00x | 3e-16 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | interpn | 7.52 µs | 0.78 ns | 9620.00x | 2e-16 |
| `grid/linear/2d/max200/eval100/ascending` | construction | interpn | 6.83 µs | 0.78 ns | 8748.00x | 1e-16 |
| `grid/linear/4d/max200/eval10/flipped` | construction | SciPy | 114 ms | 21.8 µs | 5204.93x |  |
| `grid/linear/4d/max200/eval100/flipped` | construction | SciPy | 114 ms | 22.2 µs | 5149.50x |  |
| `grid/linear/4d/max200/eval1000/ascending` | construction | SciPy | 52.8 ms | 10.4 µs | 5098.86x |  |
| `grid/linear/4d/max200/eval10/ascending` | construction | SciPy | 49.7 ms | 10.3 µs | 4811.67x |  |
| `grid/linear/4d/max200/eval10000/ascending` | construction | SciPy | 49.7 ms | 10.3 µs | 4806.04x |  |
| `grid/linear/4d/max200/eval1000/flipped` | construction | SciPy | 103 ms | 22.2 µs | 4661.41x |  |
| `grid/linear/4d/max200/eval100/ascending` | construction | SciPy | 47.9 ms | 10.3 µs | 4634.61x |  |
| `grid/linear/4d/max200/eval10000/flipped` | construction | SciPy | 100 ms | 21.8 µs | 4609.79x |  |
| `grid/linear/3d/max40/eval1000/ascending` | construction | interpn | 3.27 µs | 0.98 ns | 3351.60x | 2e-16 |
| `grid/linear/3d/max40/eval10/ascending` | construction | interpn | 2.83 µs | 0.88 ns | 3218.67x | 2e-16 |
| `grid/linear/3d/max40/eval10000/ascending` | construction | interpn | 3.05 µs | 0.98 ns | 3120.00x | 3e-16 |
| `grid/linear/3d/max40/eval100/ascending` | construction | interpn | 2.93 µs | 0.98 ns | 3001.60x | 2e-16 |
| `grid/nearest/2d/auckland-ba32-0401/heldout/ascending` | construction | SciPy | 16.8 ms | 7.05 µs | 2387.04x |  |
| `grid/linear/2d/auckland-ba32-0401/heldout/ascending` | construction | SciPy | 17.2 ms | 7.58 µs | 2272.26x |  |
| `grid/linear/2d/auckland-ba32-0401/heldout/ascending` | construction | ndarray-interp | 17.2 ms | 8.34 µs | 2065.03x | 0e+00 |
| `grid/cubic/2d/each10/eval1000/ascending` | construction | interpn | 1.61 µs | 0.88 ns | 1837.33x | 3e-01 |
| `grid/nearest/2d/each40/eval10000/ascending` | construction | interpn | 652 ns | 0.78 ns | 834.00x | 0e+00 |
| `grid/linear/2d/max40/eval1000/ascending` | construction | interpn | 397 ns | 0.68 ns | 580.07x | 3e-16 |
| `grid/linear/2d/max40/eval100/ascending` | construction | interpn | 405 ns | 0.78 ns | 518.62x | 2e-16 |
| `grid/linear/2d/max40/eval10000/ascending` | construction | interpn | 396 ns | 0.78 ns | 507.38x | 4e-16 |
| `grid/linear/2d/max40/eval10/ascending` | construction | interpn | 395 ns | 0.78 ns | 505.00x | 2e-16 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | ndarray-interp | 1.21 ms | 2.45 µs | 492.95x | 6e-16 |
| `grid/linear/3d/max200/eval10/ascending` | construction | SciPy | 2.36 ms | 8.33 µs | 283.07x |  |
| `grid/linear/2d/each400/eval10000/ascending` | construction | ndarray-interp | 333 µs | 1.26 µs | 264.91x | 5e-16 |
| `grid/linear/3d/max10/eval100/ascending` | construction | interpn | 211 ns | 0.88 ns | 240.11x | 3e-16 |
| `grid/linear/3d/max200/eval1000/ascending` | construction | SciPy | 1.82 ms | 8.22 µs | 221.17x |  |
| `grid/linear/3d/max200/eval100/ascending` | construction | SciPy | 1.82 ms | 8.31 µs | 218.73x |  |
| `grid/linear/3d/max200/eval10000/ascending` | construction | SciPy | 1.8 ms | 8.25 µs | 218.13x |  |
| `grid/linear/3d/max200/eval10/flipped` | construction | SciPy | 3.66 ms | 17.8 µs | 205.48x |  |
| `grid/linear/2d/max10/eval1000/ascending` | construction | interpn | 140 ns | 0.68 ns | 205.43x | 3e-16 |
| `grid/linear/3d/max200/eval10000/flipped` | construction | SciPy | 3.49 ms | 17.4 µs | 200.57x |  |
| `grid/linear/3d/max10/eval10000/ascending` | construction | interpn | 195 ns | 0.98 ns | 199.85x | 4e-16 |
| `grid/linear/3d/max10/eval1000/ascending` | construction | interpn | 183 ns | 0.98 ns | 187.70x | 5e-16 |
| `grid/linear/2d/max10/eval10/ascending` | construction | interpn | 140 ns | 0.78 ns | 179.12x | 1e-16 |
| `grid/linear/3d/max200/eval100/flipped` | construction | SciPy | 3.09 ms | 17.7 µs | 174.57x |  |
| `grid/linear/2d/each800/eval10000/ascending` | construction | SciPy | 1.21 ms | 6.93 µs | 174.04x |  |
| `grid/linear/3d/max10/eval10/ascending` | construction | interpn | 203 ns | 1.17 ns | 173.29x | 2e-16 |
| `grid/linear/2d/max10/eval10000/ascending` | construction | interpn | 135 ns | 0.78 ns | 173.13x | 3e-16 |
| `grid/linear/3d/max200/eval1000/flipped` | construction | SciPy | 3.02 ms | 17.5 µs | 172.79x |  |
| `grid/linear/2d/max10/eval100/ascending` | construction | interpn | 128 ns | 0.78 ns | 164.31x | 2e-16 |
| `grid/pchip/3d/each40/eval10000/ascending` | construction | SciPy | 518 µs | 9.46 µs | 54.75x |  |
| `grid/linear/2d/each400/eval10000/ascending` | construction | SciPy | 333 µs | 6.57 µs | 50.72x |  |
| `grid/linear/2d/max200/eval10000/ascending` | construction | ndarray-interp | 7.59 µs | 516 ns | 14.72x | 5e-16 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | ndarray-interp | 7.52 µs | 516 ns | 14.58x | 4e-16 |
| `grid/linear/2d/max200/eval10/ascending` | construction | ndarray-interp | 7.41 µs | 516 ns | 14.36x | 2e-16 |
| `grid/linear/2d/max200/eval100/ascending` | construction | ndarray-interp | 6.83 µs | 518 ns | 13.19x | 3e-16 |
| `linear_nd/n10/image100` | construction | spade | 5.86 µs | 1.91 µs | 3.06x | 8e-16 |
| `linear_nd/n10/image500` | construction | spade | 5.8 µs | 1.9 µs | 3.05x | 9e-16 |
| `linear_nd/n10/image50` | construction | spade | 5.81 µs | 1.91 µs | 3.05x | 6e-16 |
| `oned/pchip/n10000/shuffled-uneven` | construction | scirs2-interpolate | 138 µs | 46.8 µs | 2.96x | 3e-16 |
| `oned/akima/n10000/sorted` | construction | scirs2-interpolate | 143 µs | 48.7 µs | 2.94x | 1e-16 |
| `oned/pchip/n10000/shuffled` | construction | scirs2-interpolate | 134 µs | 46.9 µs | 2.86x | 2e-16 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | SciPy | 23.1 µs | 8.25 µs | 2.80x |  |
| `oned/pchip/n10000/sorted-uneven` | construction | scirs2-interpolate | 124 µs | 47.1 µs | 2.64x | 3e-16 |
| `grid/linear/4d/max40/eval10/flipped` | construction | SciPy | 56.5 µs | 21.5 µs | 2.62x |  |
| `grid/linear/4d/max40/eval100/flipped` | construction | SciPy | 56.5 µs | 21.6 µs | 2.61x |  |
| `grid/linear/4d/max40/eval1000/flipped` | construction | SciPy | 57.2 µs | 21.9 µs | 2.61x |  |
| `grid/linear/4d/max40/eval10000/flipped` | construction | SciPy | 56.8 µs | 21.8 µs | 2.61x |  |
| `grid/linear/4d/max40/eval10000/ascending` | construction | scirs2-interpolate | 14.4 µs | 5.65 µs | 2.54x | 4e-16 |
| `grid/linear/4d/max40/eval10/ascending` | construction | scirs2-interpolate | 14.4 µs | 5.99 µs | 2.40x | 2e-16 |
| `grid/linear/2d/max40/eval100/ascending` | construction | ndarray-interp | 405 ns | 169 ns | 2.40x | 2e-16 |
| `grid/linear/4d/max40/eval1000/ascending` | construction | scirs2-interpolate | 14.3 µs | 6 µs | 2.39x | 2e-16 |
| `grid/linear/2d/max200/eval10/ascending` | construction | scirs2-interpolate | 7.41 µs | 3.11 µs | 2.38x | 1e-16 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | scirs2-interpolate | 23.1 µs | 9.69 µs | 2.38x | 0e+00 |
| `grid/linear/4d/max40/eval100/ascending` | construction | scirs2-interpolate | 14.3 µs | 6.01 µs | 2.38x | 3e-16 |
| `grid/linear/2d/max200/eval10000/ascending` | construction | scirs2-interpolate | 7.59 µs | 3.23 µs | 2.35x | 2e-16 |
| `grid/linear/2d/max40/eval10000/ascending` | construction | ndarray-interp | 396 ns | 169 ns | 2.35x | 4e-16 |
| `grid/linear/2d/max40/eval1000/ascending` | construction | ndarray-interp | 397 ns | 169 ns | 2.35x | 4e-16 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | scirs2-interpolate | 7.52 µs | 3.23 µs | 2.32x | 2e-16 |
| `grid/linear/2d/max40/eval10/ascending` | construction | ndarray-interp | 395 ns | 171 ns | 2.31x | 2e-16 |
| `linear_nd/n10000/image500` | construction | spade | 8.31 ms | 3.72 ms | 2.23x | 3e-15 |
| `oned/akima/n10000/shuffled` | construction | scirs2-interpolate | 101 µs | 45.9 µs | 2.21x | 1e-16 |
| `grid/nearest/2d/each40/eval10000/ascending` | construction | scirs2-interpolate | 652 ns | 296 ns | 2.20x | 0e+00 |
| `grid/linear/2d/max200/eval100/ascending` | construction | scirs2-interpolate | 6.83 µs | 3.2 µs | 2.14x | 2e-16 |
| `oned/linear/n10/sorted` | construction | ndarray-interp | 87.2 ns | 41.5 ns | 2.10x | 0e+00 |
| `oned/pchip/n10000/sorted` | construction | scirs2-interpolate | 98.2 µs | 47.3 µs | 2.08x | 2e-16 |
| `oned/linear/n10/shuffled` | construction | ndarray-interp | 93.9 ns | 45.5 ns | 2.06x | 0e+00 |
| `grid/linear/3d/max40/eval1000/ascending` | construction | scirs2-interpolate | 3.27 µs | 1.6 µs | 2.05x | 2e-16 |
| `linear_nd/n1000/image100` | construction | spade | 618 µs | 302 µs | 2.05x | 1e-15 |
| `linear_nd/n1000/image50` | construction | spade | 614 µs | 302 µs | 2.03x | 1e-15 |
| `grid/linear/3d/max40/eval10/ascending` | construction | scirs2-interpolate | 2.83 µs | 1.4 µs | 2.02x | 1e-16 |
| `oned/pchip/n1000/shuffled` | construction | scirs2-interpolate | 9.17 µs | 4.61 µs | 1.99x | 2e-16 |
| `oned/pchip/n1000/sorted` | construction | scirs2-interpolate | 9.11 µs | 4.59 µs | 1.98x | 2e-16 |
| `linear_nd/n1000/image500` | construction | spade | 605 µs | 308 µs | 1.96x | 3e-15 |
| `grid/linear/3d/max40/eval10000/ascending` | construction | scirs2-interpolate | 3.05 µs | 1.57 µs | 1.94x | 2e-16 |
| `oned/akima/n1000/shuffled` | construction | scirs2-interpolate | 9.17 µs | 4.83 µs | 1.90x | 3e-17 |
| `linear_nd/n10000/image100` | construction | spade | 7.72 ms | 4.08 ms | 1.89x | 8e-16 |
| `oned/akima/n1000/sorted` | construction | scirs2-interpolate | 9.18 µs | 4.85 µs | 1.89x | 3e-17 |
| `grid/linear/3d/max200/eval10/ascending` | evaluation | interpn | 532 ns | 287 ns | 1.85x | 3e-16 |
| `grid/linear/2d/max40/eval100/ascending` | construction | scirs2-interpolate | 405 ns | 220 ns | 1.84x | 3e-16 |
| `linear_nd/n10000/image50` | construction | spade | 7.59 ms | 4.13 ms | 1.84x | 8e-16 |
| `grid/linear/3d/max40/eval100/ascending` | construction | scirs2-interpolate | 2.93 µs | 1.6 µs | 1.83x | 1e-16 |
| `grid/linear/2d/max40/eval10000/ascending` | construction | scirs2-interpolate | 396 ns | 218 ns | 1.82x | 2e-16 |
| `grid/linear/2d/max40/eval1000/ascending` | construction | scirs2-interpolate | 397 ns | 222 ns | 1.78x | 3e-16 |
| `grid/linear/2d/max40/eval10/ascending` | construction | scirs2-interpolate | 395 ns | 222 ns | 1.78x | 4e-17 |
| `oned/linear/n10000/shuffled` | construction | scirs2-interpolate | 22.9 µs | 12.9 µs | 1.77x | 0e+00 |
| `linear_nd/n100/image100` | construction | spade | 35.8 µs | 20.6 µs | 1.74x | 9e-16 |
| `linear_nd/n100/image500` | construction | spade | 35.9 µs | 20.7 µs | 1.74x | 1e-15 |
| `linear_nd/n100/image50` | construction | spade | 35.8 µs | 20.8 µs | 1.72x | 1e-15 |
| `oned/linear/n1000/shuffled` | construction | scirs2-interpolate | 2.17 µs | 1.28 µs | 1.70x | 0e+00 |
| `linear_nd/n100000/image100` | construction | spade | 82.6 ms | 48.8 ms | 1.69x | 7e-16 |
| `linear_nd/ghcn-tmax-20250101/heldout` | construction | spade | 7.76 ms | 4.6 ms | 1.69x | 7e-16 |
| `oned/linear/n10000/sorted` | construction | scirs2-interpolate | 21.7 µs | 12.9 µs | 1.68x | 0e+00 |
| `linear_nd/n100000/image500` | construction | spade | 77.8 ms | 46.7 ms | 1.67x | 6e-15 |
| `oned/linear/n1000/sorted` | construction | scirs2-interpolate | 2.13 µs | 1.28 µs | 1.66x | 0e+00 |
| `oned/linear/n10000/sorted-uneven` | construction | scirs2-interpolate | 21.5 µs | 12.9 µs | 1.66x | 6e-17 |
| `oned/linear/n10000/shuffled-uneven` | construction | scirs2-interpolate | 21.5 µs | 12.9 µs | 1.66x | 6e-17 |
| `oned/akima/n100/shuffled` | construction | scirs2-interpolate | 973 ns | 587 ns | 1.66x | 1e-16 |
| `oned/cubic_spline/n10000/sorted` | construction | scirs2-interpolate | 240 µs | 145 µs | 1.65x | 4e-05 |
| `oned/pchip/n100/shuffled-uneven` | construction | scirs2-interpolate | 980 ns | 594 ns | 1.65x | 2e-16 |
| `oned/pchip/n100/sorted` | construction | scirs2-interpolate | 981 ns | 596 ns | 1.65x | 2e-16 |
| `grid/linear/2d/max10/eval1000/ascending` | construction | ndarray-interp | 140 ns | 85.4 ns | 1.64x | 3e-16 |
| `oned/akima/n100/sorted` | construction | scirs2-interpolate | 967 ns | 589 ns | 1.64x | 1e-16 |
| `grid/linear/2d/max10/eval10/ascending` | construction | ndarray-interp | 140 ns | 85.5 ns | 1.64x | 1e-16 |
| `oned/pchip/n100/shuffled` | construction | scirs2-interpolate | 969 ns | 594 ns | 1.63x | 2e-16 |
| `oned/pchip/n100/sorted-uneven` | construction | scirs2-interpolate | 968 ns | 594 ns | 1.63x | 2e-16 |
| `grid/linear/2d/max10/eval10000/ascending` | construction | ndarray-interp | 135 ns | 85.5 ns | 1.58x | 4e-16 |
| `linear_nd/n100000/image50` | construction | spade | 76.6 ms | 48.7 ms | 1.57x | 8e-16 |
| `oned/linear/n50/shuffled` | construction | ndarray-interp | 164 ns | 104 ns | 1.57x | 0e+00 |
| `oned/akima/n50/shuffled` | construction | scirs2-interpolate | 552 ns | 366 ns | 1.51x | 1e-16 |
| `grid/linear/4d/max40/eval100/ascending` | construction | SciPy | 14.3 µs | 9.51 µs | 1.50x |  |
| `grid/linear/2d/max10/eval100/ascending` | construction | ndarray-interp | 128 ns | 85.6 ns | 1.50x | 4e-16 |
| `oned/akima/n50/sorted` | construction | scirs2-interpolate | 547 ns | 369 ns | 1.48x | 1e-16 |
| `oned/linear/n100/shuffled` | construction | ndarray-interp | 271 ns | 183 ns | 1.48x | 0e+00 |
| `oned/linear/n100/sorted` | construction | ndarray-interp | 270 ns | 183 ns | 1.48x | 0e+00 |
| `oned/linear/n50/sorted` | construction | ndarray-interp | 161 ns | 110 ns | 1.47x | 0e+00 |
| `oned/linear/n100/sorted-uneven` | construction | ndarray-interp | 267 ns | 183 ns | 1.46x | 0e+00 |
| `grid/linear/2d/max10/eval10/ascending` | evaluation | interpn | 150 ns | 104 ns | 1.44x | 1e-16 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | scirs2-interpolate | 1.21 ms | 838 µs | 1.44x | 2e-16 |
| `oned/linear/n1000/shuffled` | construction | ndarray-interp | 2.17 µs | 1.52 µs | 1.43x | 0e+00 |
| `grid/linear/4d/max40/eval10/ascending` | construction | SciPy | 14.4 µs | 10.1 µs | 1.42x |  |
| `oned/linear/n1000/sorted` | construction | ndarray-interp | 2.13 µs | 1.51 µs | 1.40x | 0e+00 |
| `grid/linear/4d/max40/eval1000/ascending` | construction | SciPy | 14.3 µs | 10.2 µs | 1.40x |  |
| `grid/linear/4d/max40/eval10000/ascending` | construction | SciPy | 14.4 µs | 10.3 µs | 1.40x |  |
| `oned/linear/n100/sorted` | construction | scirs2-interpolate | 270 ns | 195 ns | 1.39x | 0e+00 |
| `grid/linear/3d/max10/eval100/ascending` | evaluation | interpn | 2.09 µs | 1.53 µs | 1.37x | 3e-16 |
| `oned/linear/n10000/shuffled` | construction | ndarray-interp | 22.9 µs | 16.8 µs | 1.36x | 0e+00 |
| `oned/pchip/n50/shuffled` | construction | scirs2-interpolate | 532 ns | 390 ns | 1.36x | 1e-16 |
| `grid/linear/2d/max10/eval100/ascending` | evaluation | interpn | 1.38 µs | 1.02 µs | 1.36x | 2e-16 |
| `grid/linear/4d/max200/eval1000/ascending` | construction | scirs2-interpolate | 52.8 ms | 39 ms | 1.35x | 3e-16 |
| `grid/linear/2d/max10/eval10000/ascending` | evaluation | interpn | 137 µs | 101 µs | 1.35x | 3e-16 |
| `grid/linear/3d/max10/eval1000/ascending` | evaluation | interpn | 20.7 µs | 15.3 µs | 1.35x | 5e-16 |
| `grid/linear/2d/max10/eval1000/ascending` | evaluation | interpn | 13.7 µs | 10.1 µs | 1.35x | 3e-16 |
| `grid/linear/3d/max200/eval10000/ascending` | construction | scirs2-interpolate | 1.8 ms | 1.34 ms | 1.35x | 3e-16 |
| `grid/linear/3d/max10/eval10000/ascending` | evaluation | interpn | 209 µs | 155 µs | 1.35x | 4e-16 |
| `grid/linear/3d/max200/eval1000/ascending` | construction | scirs2-interpolate | 1.82 ms | 1.38 ms | 1.32x | 3e-16 |
| `oned/linear/n100/shuffled` | construction | scirs2-interpolate | 271 ns | 209 ns | 1.30x | 0e+00 |
| `oned/pchip/n50/sorted` | construction | scirs2-interpolate | 481 ns | 371 ns | 1.30x | 1e-16 |
| `grid/linear/3d/max10/eval10/ascending` | evaluation | interpn | 200 ns | 156 ns | 1.28x | 2e-16 |
| `oned/linear/n100/sorted-uneven` | construction | scirs2-interpolate | 267 ns | 208 ns | 1.28x | 3e-17 |
| `grid/linear/3d/max200/eval10/ascending` | construction | scirs2-interpolate | 2.36 ms | 1.84 ms | 1.28x | 3e-16 |
| `oned/linear/n50/shuffled` | construction | scirs2-interpolate | 164 ns | 128 ns | 1.28x | 0e+00 |
| `grid/linear/2d/max40/eval1000/ascending` | evaluation | interpn | 16.6 µs | 13.2 µs | 1.26x | 3e-16 |
| `grid/linear/2d/auckland-ba32-0401/heldout/ascending` | construction | scirs2-interpolate | 17.2 ms | 13.8 ms | 1.24x | 0e+00 |
| `grid/linear/4d/max200/eval10/ascending` | construction | scirs2-interpolate | 49.7 ms | 40 ms | 1.24x | 2e-16 |
| `grid/nearest/2d/auckland-ba32-0401/heldout/ascending` | construction | scirs2-interpolate | 16.8 ms | 13.7 ms | 1.23x | 0e+00 |
| `grid/linear/2d/max200/eval10000/flipped` | construction | SciPy | 15.8 µs | 12.9 µs | 1.23x |  |
| `grid/linear/2d/max200/eval100/flipped` | construction | SciPy | 16.1 µs | 13.2 µs | 1.22x |  |
| `grid/linear/4d/max200/eval10000/ascending` | construction | scirs2-interpolate | 49.7 ms | 40.5 ms | 1.22x | 4e-16 |
| `grid/linear/4d/max200/eval100/ascending` | construction | scirs2-interpolate | 47.9 ms | 39.3 ms | 1.22x | 1e-16 |
| `oned/linear/n50/sorted` | construction | scirs2-interpolate | 161 ns | 133 ns | 1.21x | 0e+00 |
| `grid/linear/2d/max40/eval100/ascending` | evaluation | interpn | 1.66 µs | 1.38 µs | 1.21x | 2e-16 |
| `grid/linear/3d/max40/eval10/ascending` | evaluation | interpn | 240 ns | 199 ns | 1.20x | 2e-16 |
| `oned/cubic_spline/n10000/sorted-uneven` | construction | scirs2-interpolate | 174 µs | 146 µs | 1.20x | 7e-06 |
| `grid/linear/2d/max200/eval10/flipped` | construction | SciPy | 14.6 µs | 12.2 µs | 1.19x |  |
| `grid/linear/2d/max200/eval1000/flipped` | construction | SciPy | 15.6 µs | 13.1 µs | 1.19x |  |
| `grid/linear/2d/max200/eval10000/ascending` | construction | SciPy | 7.59 µs | 6.4 µs | 1.19x |  |
| `grid/linear/2d/max200/eval1000/ascending` | construction | SciPy | 7.52 µs | 6.4 µs | 1.17x |  |
| `oned/cubic_spline/n10000/shuffled-uneven` | construction | scirs2-interpolate | 170 µs | 146 µs | 1.17x | 7e-06 |
| `grid/linear/2d/max40/eval10/ascending` | evaluation | interpn | 172 ns | 147 ns | 1.17x | 2e-16 |
| `oned/cubic_spline/n1000/sorted` | construction | scirs2-interpolate | 17 µs | 14.6 µs | 1.17x | 4e-04 |
| `oned/cubic_spline/n1000/shuffled` | construction | scirs2-interpolate | 16.9 µs | 14.6 µs | 1.16x | 4e-04 |
| `oned/cubic_spline/n10000/shuffled` | construction | scirs2-interpolate | 169 µs | 147 µs | 1.16x | 4e-05 |
| `grid/linear/2d/max200/eval10/ascending` | construction | SciPy | 7.41 µs | 6.47 µs | 1.15x |  |
| `grid/linear/4d/max40/eval10/ascending` | evaluation | interpn | 338 ns | 296 ns | 1.14x | 2e-16 |
| `grid/linear/2d/auckland-ba32-0401/heldout/ascending` | evaluation | ndarray-interp | 426 ms | 375 ms | 1.14x | 0e+00 |
| `grid/linear/3d/max10/eval100/ascending` | construction | scirs2-interpolate | 211 ns | 188 ns | 1.12x | 2e-16 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | scirs2-interpolate | 333 µs | 298 µs | 1.12x | 2e-16 |
| `grid/linear/2d/max200/eval10/ascending` | evaluation | ndarray-interp | 211 ns | 189 ns | 1.12x | 2e-16 |
| `grid/nearest/2d/each40/eval10000/ascending` | evaluation | interpn | 166 µs | 150 µs | 1.11x | 0e+00 |
| `grid/linear/4d/max40/eval10000/ascending` | evaluation | interpn | 327 µs | 297 µs | 1.10x | 5e-16 |

</details>

## Runs that failed

| case | library | how |
| :--- | :--- | :--- |
| `grid/linear/4d/max10/eval10/ascending` | interpn | the library refused the case: "All grids must have at least 2 entries" |
| `grid/linear/4d/max10/eval10/ascending` | scirs2-interpolate | the library refused the case: InvalidInput { message: "Dimension 3 has less than 2 points" } |
| `grid/linear/4d/max10/eval100/ascending` | interpn | the library refused the case: "All grids must have at least 2 entries" |
| `grid/linear/4d/max10/eval100/ascending` | scirs2-interpolate | the library refused the case: InvalidInput { message: "Dimension 3 has less than 2 points" } |
| `grid/linear/4d/max10/eval1000/ascending` | interpn | the library refused the case: "All grids must have at least 2 entries" |
| `grid/linear/4d/max10/eval1000/ascending` | scirs2-interpolate | the library refused the case: InvalidInput { message: "Dimension 3 has less than 2 points" } |
| `grid/linear/4d/max10/eval10000/ascending` | interpn | the library refused the case: "All grids must have at least 2 entries" |
| `grid/linear/4d/max10/eval10000/ascending` | scirs2-interpolate | the library refused the case: InvalidInput { message: "Dimension 3 has less than 2 points" } |
| `grid/pchip/2d/auckland-ba32-0401/heldout/ascending` | SciPy | MemoryError: Unable to allocate 695. GiB for an array with shape (25913999, 3601) and data type float64 |

## Runs near the time limit

A run (one library on one case, a process of its own) is stopped at 15 minutes. The longest run that finished took 369 s, 41% of that: scirs2-interpolate on `rbf/cubic/global/ghcn-tmax-20250101/heldout`.

No run took more than half the limit (450 s).

## Real-world data

Two data sets, each taken whole in the unit its source gives it out in; what was taken and how it was divided was fixed before anything was measured (`compare/perf/real.py`, docs/design.ja.md section 8.4.7).

- `auckland-ba32-0401`: Auckland Part 1 LiDAR 1m DEM (2024), tile BA32_10000_0401. 4800 by 7200 cells of 1 m. The grid is every second row and column (and the last ones); every other cell is a query, and its elevation is the value held out. Licensed by Regional Software Holdings Limited for re-use under CC-BY-4.0. Hosted by Toitū Te Whenua Land Information New Zealand (https://registry.opendata.aws/nz-elevation).
- `ghcn-tmax-20250101`: GHCN-Daily, maximum temperature of 2025-01-01, every station; longitude and latitude in degrees, as they are. Of every 10 stations in the order of the file the last is a query, and its temperature is the value held out. The RBF kernels are those without a shape parameter. NOAA Global Historical Climatology Network Daily (GHCN-D), from https://registry.opendata.aws/noaa-ghcn (CC0 1.0).

These cases are also in the tables of their families below and in the counts above. Here they are with the error: the largest and the root mean square distance from the values held out (metres, and degrees Celsius), and the number of queries a library gave no value. A library missing from a case does not have it, or did not finish (see the runs that failed).

`grid/nearest/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 16.8 ms | 376 ms | 15.5 | 0.178 | 0 |
| interpx, parallel | 15.5 ms | 69.1 ms | 15.5 | 0.178 | 0 |
| interpx, `from_owned` | 5.34 ms | 389 ms | 15.5 | 0.178 | 0 |
| SciPy | 7.05 µs | 1.41 s | 15.5 | 0.178 | 0 |
| scirs2-interpolate | 13.7 ms | 172 s | 15.5 | 0.178 | 0 |
| interpn | 0.78 ns | 781 ms | 15.5 | 0.178 | 0 |

`grid/linear/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 17.2 ms | 426 ms | 6.74 | 0.0601 | 0 |
| interpx, parallel | 17.6 ms | 74.6 ms | 6.74 | 0.0601 | 0 |
| interpx, `from_owned` | 6.04 ms | 452 ms | 6.74 | 0.0601 | 0 |
| SciPy | 7.58 µs | 1.21 s | 6.74 | 0.0601 | 0 |
| scirs2-interpolate | 13.8 ms | 81.6 s | 6.74 | 0.0601 | 0 |
| ndarray-interp | 8.34 µs | 375 ms | 6.74 | 0.0601 | 0 |
| interpn | 0.78 ns | 888 ms | 6.74 | 0.0601 | 0 |

`grid/cubic/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 78.8 ms | 803 ms | 6.77 | 0.0507 | 0 |
| interpx, parallel | 79.1 ms | 137 ms | 6.77 | 0.0507 | 0 |
| interpx, `from_owned` | 74.3 ms | 859 ms | 6.77 | 0.0507 | 0 |
| SciPy | 17.8 s | 36.3 s | 6.77 | 0.0507 | 0 |
| interpn | 0.88 ns | 1.83 s | 6.76 | 0.0517 | 0 |

`grid/quintic/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 109 ms | 1.35 s | 6.77 | 0.0517 | 0 |
| interpx, parallel | 97.8 ms | 235 ms | 6.77 | 0.0517 | 0 |
| interpx, `from_owned` | 90.7 ms | 1.32 s | 6.77 | 0.0517 | 0 |
| SciPy | 23.7 s | 39 s | 6.77 | 0.0517 | 0 |

`grid/pchip/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 78.8 ms | 1.1 s | 6.75 | 0.049 | 0 |
| interpx, parallel | 76.6 ms | 184 ms | 6.75 | 0.049 | 0 |
| interpx, `from_owned` | 66.2 ms | 1.09 s | 6.75 | 0.049 | 0 |

`nearest_nd/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 533 µs | 172 µs | 25 | 3.42 | 0 |
| interpx, parallel | 561 µs | 72.1 µs | 25 | 3.42 | 0 |
| interpx, `from_owned` | 559 µs | 173 µs | 25 | 3.42 | 0 |
| SciPy | 1.84 ms | 546 µs | 25 | 3.42 | 0 |

`linear_nd/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 7.76 ms | 99.5 µs | 22.1 | 2.83 | 1 |
| interpx, parallel | 7.55 ms | 105 µs | 22.1 | 2.83 | 1 |
| interpx, `from_owned` | 7.47 ms | 95.7 µs | 22.1 | 2.83 | 1 |
| SciPy | 31.4 ms | 1.19 ms | 22.1 | 2.83 | 1 |
| spade | 4.6 ms | 1.23 ms | 22.1 | 2.83 | 1 |

`rbf/linear/global/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 19.3 s | 20.2 ms | 20.8 | 2.74 | 0 |
| interpx, parallel | 18.8 s | 3.78 ms | 20.8 | 2.74 | 0 |
| interpx, `from_owned` | 19.1 s | 20.3 ms | 20.8 | 2.74 | 0 |
| SciPy | 128 s | 55.5 ms | 20.8 | 2.74 | 0 |
| scirs2-interpolate | 354 s | 252 ms | 20.8 | 2.74 | 0 |

`rbf/linear/k50/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 749 µs | 17.9 ms | 20.8 | 2.75 | 0 |
| interpx, parallel | 829 µs | 3.49 ms | 20.8 | 2.75 | 0 |
| interpx, `from_owned` | 704 µs | 18.2 ms | 20.8 | 2.75 | 0 |
| SciPy | 1.86 ms | 45.4 ms | 20.8 | 2.75 | 0 |

`rbf/thin_plate_spline/global/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 20.5 s | 87.4 ms | 102 | 5.67 | 0 |
| interpx, parallel | 20.2 s | 17.1 ms | 102 | 5.67 | 0 |
| interpx, `from_owned` | 21.2 s | 94.1 ms | 102 | 5.67 | 0 |
| SciPy | 141 s | 144 ms | 102 | 5.67 | 0 |
| scirs2-interpolate | 354 s | 462 ms | 93.7 | 5.41 | 0 |

`rbf/thin_plate_spline/k50/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 724 µs | 29.3 ms | 102 | 5.66 | 0 |
| interpx, parallel | 750 µs | 5.43 ms | 102 | 5.66 | 0 |
| interpx, `from_owned` | 701 µs | 31.7 ms | 102 | 5.66 | 0 |
| SciPy | 1.84 ms | 57.2 ms | 102 | 5.66 | 0 |

`rbf/cubic/global/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 20.5 s | 25.2 ms | 587 | 27.1 | 0 |
| interpx, parallel | 23 s | 4.51 ms | 587 | 27.1 | 0 |
| interpx, `from_owned` | 23.8 s | 25.6 ms | 587 | 27.1 | 0 |
| SciPy | 131 s | 57.5 ms | 597 | 24.4 | 0 |
| scirs2-interpolate | 367 s | 295 ms | 100 | 8.48 | 0 |

`rbf/cubic/k50/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 715 µs | 22.4 ms | 567 | 25.3 | 0 |
| interpx, parallel | 837 µs | 3.94 ms | 567 | 25.3 | 0 |
| interpx, `from_owned` | 724 µs | 21.7 ms | 567 | 25.3 | 0 |
| SciPy | 1.91 ms | 51.3 ms | 567 | 25.3 | 0 |

`rbf/quintic/global/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 18.9 s | 22 ms | 2.55e+05 | 1.36e+04 | 0 |
| interpx, parallel | 21.4 s | 4.27 ms | 2.55e+05 | 1.36e+04 | 0 |
| interpx, `from_owned` | 18.1 s | 22.3 ms | 2.55e+05 | 1.36e+04 | 0 |
| SciPy | 161 s | 63.5 ms | 9.39e+04 | 4.92e+03 | 0 |
| scirs2-interpolate | 310 s | 294 ms | 6.56e+16 | 2.85e+15 | 0 |

`rbf/quintic/k50/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 833 µs | 22.7 ms | 3.74e+03 | 199 | 0 |
| interpx, parallel | 751 µs | 4.13 ms | 3.74e+03 | 199 | 0 |
| interpx, `from_owned` | 706 µs | 21.8 ms | 3.74e+03 | 199 | 0 |
| SciPy | 1.74 ms | 54 ms | 3.74e+03 | 200 | 0 |

## One variable

Data of `n` points, evaluated at `4 n` points: in order (`sorted`) and in no order (`shuffled`).

### Evaluation

| case | interpx | interpx, parallel | SciPy | scirs2-interpolate | ndarray-interp | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/n10/sorted` | 112 ns | 190 ns | 5 µs | 304 ns | 279 ns | 2.5x faster |
| `linear/n10/shuffled` | 188 ns | 247 ns | 5.61 µs | 305 ns | 274 ns | 1.5x faster |
| `nearest/n10/sorted` | 171 ns | 223 ns | 10.3 µs | 296 ns |  | 1.7x faster |
| `nearest/n10/shuffled` | 330 ns | 400 ns | 10.4 µs | 301 ns |  | **1.1x slower** |
| `previous/n10/sorted` | 170 ns | 213 ns | 11.7 µs |  |  | 69x faster |
| `previous/n10/shuffled` | 291 ns | 355 ns | 11.6 µs |  |  | 40x faster |
| `bspline1/n10/sorted` | 229 ns | 317 ns | 2.96 µs |  |  | 13x faster |
| `bspline1/n10/shuffled` | 579 ns | 655 ns | 3 µs |  |  | 5.2x faster |
| `bspline2/n10/sorted` | 271 ns | 374 ns | 3.4 µs |  |  | 13x faster |
| `bspline2/n10/shuffled` | 372 ns | 453 ns | 3.42 µs |  |  | 9.2x faster |
| `bspline3/n10/sorted` | 377 ns | 476 ns | 3.94 µs |  |  | 10x faster |
| `bspline3/n10/shuffled` | 481 ns | 565 ns | 3.96 µs |  |  | 8.2x faster |
| `cubic_spline/n10/sorted` | 153 ns | 210 ns | 3.94 µs | 304 ns | 391 ns | 2x faster |
| `cubic_spline/n10/shuffled` | 225 ns | 292 ns | 3.94 µs | 309 ns | 333 ns | 1.4x faster |
| `pchip/n10/sorted` | 156 ns | 210 ns | 3.83 µs | 320 ns |  | 2.1x faster |
| `pchip/n10/shuffled` | 216 ns | 265 ns | 3.72 µs | 318 ns |  | 1.5x faster |
| `akima/n10/sorted` | 154 ns | 210 ns | 3.76 µs | 268 ns |  | 1.7x faster |
| `akima/n10/shuffled` | 226 ns | 282 ns | 3.65 µs | 281 ns |  | 1.2x faster |
| `makima/n10/sorted` | 140 ns | 210 ns | 3.65 µs |  |  | 26x faster |
| `makima/n10/shuffled` | 225 ns | 256 ns | 3.98 µs |  |  | 18x faster |
| `linear/n50/sorted` | 662 ns | 741 ns | 6.15 µs | 1.66 µs | 1.35 µs | 2x faster |
| `linear/n50/shuffled` | 891 ns | 943 ns | 8.9 µs | 1.84 µs | 1.38 µs | 1.6x faster |
| `nearest/n50/sorted` | 875 ns | 906 ns | 11.4 µs | 1.79 µs |  | 2x faster |
| `nearest/n50/shuffled` | 1.69 µs | 1.76 µs | 12.3 µs | 1.86 µs |  | 1.1x faster |
| `previous/n50/sorted` | 884 ns | 885 ns | 12.6 µs |  |  | 14x faster |
| `previous/n50/shuffled` | 1.39 µs | 1.57 µs | 14 µs |  |  | 10x faster |
| `bspline1/n50/sorted` | 1.16 µs | 1.38 µs | 4.3 µs |  |  | 3.7x faster |
| `bspline1/n50/shuffled` | 3.05 µs | 3.14 µs | 4.96 µs |  |  | 1.6x faster |
| `bspline2/n50/sorted` | 1.48 µs | 1.75 µs | 6.67 µs |  |  | 4.5x faster |
| `bspline2/n50/shuffled` | 2.22 µs | 2.3 µs | 7.37 µs |  |  | 3.3x faster |
| `bspline3/n50/sorted` | 1.86 µs | 2.07 µs | 9.48 µs |  |  | 5.1x faster |
| `bspline3/n50/shuffled` | 3.09 µs | 3.17 µs | 10.4 µs |  |  | 3.4x faster |
| `cubic_spline/n50/sorted` | 791 ns | 815 ns | 5.79 µs | 4.8 µs | 1.76 µs | 2.2x faster |
| `cubic_spline/n50/shuffled` | 1.08 µs | 1.14 µs | 6.1 µs | 4.44 µs | 1.76 µs | 1.6x faster |
| `pchip/n50/sorted` | 721 ns | 828 ns | 5.29 µs | 4.3 µs |  | 6x faster |
| `pchip/n50/shuffled` | 1.08 µs | 1.14 µs | 5.83 µs | 4.04 µs |  | 3.7x faster |
| `akima/n50/sorted` | 786 ns | 851 ns | 5.65 µs | 5.05 µs |  | 6.4x faster |
| `akima/n50/shuffled` | 1.08 µs | 1.15 µs | 5.95 µs | 4.47 µs |  | 4.1x faster |
| `makima/n50/sorted` | 791 ns | 850 ns | 5.7 µs |  |  | 7.2x faster |
| `makima/n50/shuffled` | 1.07 µs | 1.14 µs | 6.06 µs |  |  | 5.7x faster |
| `linear/n100/sorted` | 1.39 µs | 1.49 µs | 8.01 µs | 3.53 µs | 2.7 µs | 1.9x faster |
| `linear/n100/shuffled` | 1.77 µs | 1.89 µs | 14.3 µs | 4.09 µs | 2.73 µs | 1.5x faster |
| `nearest/n100/sorted` | 1.73 µs | 1.85 µs | 14.6 µs | 3.82 µs |  | 2.2x faster |
| `nearest/n100/shuffled` | 3.38 µs | 3.55 µs | 15.7 µs | 4.08 µs |  | 1.2x faster |
| `previous/n100/sorted` | 1.78 µs | 2.79 µs | 27.1 µs |  |  | 15x faster |
| `previous/n100/shuffled` | 3.01 µs | 3.15 µs | 18.1 µs |  |  | 6x faster |
| `bspline1/n100/sorted` | 4.24 µs | 3.02 µs | 6.22 µs |  |  | 1.5x faster |
| `bspline1/n100/shuffled` | 6.13 µs | 6.33 µs | 13.5 µs |  |  | 2.2x faster |
| `bspline2/n100/sorted` | 3.05 µs | 3.51 µs | 11 µs |  |  | 3.6x faster |
| `bspline2/n100/shuffled` | 4.87 µs | 5.02 µs | 14.9 µs |  |  | 3.1x faster |
| `bspline3/n100/sorted` | 3.99 µs | 4.15 µs | 16.5 µs |  |  | 4.1x faster |
| `bspline3/n100/shuffled` | 6.74 µs | 6.89 µs | 19.8 µs |  |  | 2.9x faster |
| `cubic_spline/n100/sorted` | 1.6 µs | 1.69 µs | 7.55 µs | 19.2 µs | 3.51 µs | 2.2x faster |
| `cubic_spline/n100/shuffled` | 2.15 µs | 2.26 µs | 8.71 µs | 19.9 µs | 6.48 µs | 3x faster |
| `pchip/n100/sorted` | 2.66 µs | 1.67 µs | 7.75 µs | 17.4 µs |  | 2.9x faster |
| `pchip/n100/shuffled` | 2.15 µs | 2.27 µs | 8.63 µs | 17.6 µs |  | 4x faster |
| `akima/n100/sorted` | 1.61 µs | 1.68 µs | 7.58 µs | 18.8 µs |  | 4.7x faster |
| `akima/n100/shuffled` | 2.15 µs | 2.29 µs | 8.68 µs | 19.7 µs |  | 4x faster |
| `makima/n100/sorted` | 1.62 µs | 1.71 µs | 7.68 µs |  |  | 4.8x faster |
| `makima/n100/shuffled` | 2.16 µs | 2.26 µs | 8.77 µs |  |  | 4.1x faster |
| `linear/n1000/sorted` | 14.1 µs | 14.3 µs | 32.6 µs | 109 µs | 27.2 µs | 1.9x faster |
| `linear/n1000/shuffled` | 18 µs | 17.8 µs | 186 µs | 233 µs | 27.1 µs | 1.5x faster |
| `nearest/n1000/sorted` | 16.7 µs | 22.7 µs | 62.1 µs | 147 µs |  | 3.7x faster |
| `nearest/n1000/shuffled` | 34 µs | 34.2 µs | 219 µs | 244 µs |  | 6.4x faster |
| `previous/n1000/sorted` | 16.6 µs | 17.2 µs | 65.8 µs |  |  | 4x faster |
| `previous/n1000/shuffled` | 30.2 µs | 30.4 µs | 225 µs |  |  | 7.5x faster |
| `bspline1/n1000/sorted` | 26.6 µs | 26.4 µs | 34 µs |  |  | 1.3x faster |
| `bspline1/n1000/shuffled` | 61.8 µs | 62.1 µs | 583 µs |  |  | 9.4x faster |
| `bspline2/n1000/sorted` | 34.1 µs | 34.2 µs | 85 µs |  |  | 2.5x faster |
| `bspline2/n1000/shuffled` | 63.8 µs | 73.8 µs | 597 µs |  |  | 9.4x faster |
| `bspline3/n1000/sorted` | 39.3 µs | 40 µs | 169 µs |  |  | 4.3x faster |
| `bspline3/n1000/shuffled` | 87 µs | 87 µs | 638 µs |  |  | 7.3x faster |
| `cubic_spline/n1000/sorted` | 15.9 µs | 16.1 µs | 41 µs | 1.56 ms | 35 µs | 2.2x faster |
| `cubic_spline/n1000/shuffled` | 21.9 µs | 22.2 µs | 240 µs | 1.56 ms | 32.7 µs | 1.5x faster |
| `pchip/n1000/sorted` | 15.9 µs | 16.1 µs | 41.3 µs | 1.21 ms |  | 2.6x faster |
| `pchip/n1000/shuffled` | 21.9 µs | 22.1 µs | 239 µs | 1.32 ms |  | 11x faster |
| `akima/n1000/sorted` | 16 µs | 15.9 µs | 41.6 µs | 1.56 ms |  | 2.6x faster |
| `akima/n1000/shuffled` | 21.8 µs | 22 µs | 242 µs | 1.57 ms |  | 11x faster |
| `makima/n1000/sorted` | 16.2 µs | 16 µs | 41.3 µs |  |  | 2.6x faster |
| `makima/n1000/shuffled` | 21.8 µs | 22.1 µs | 241 µs |  |  | 11x faster |
| `linear/n10000/sorted` | 148 µs | 75.6 µs | 242 µs | 1.4 ms | 272 µs | 1.6x faster |
| `linear/n10000/shuffled` | 208 µs | 91.4 µs | 2.42 ms | 3.09 ms | 286 µs | 1.4x faster |
| `nearest/n10000/sorted` | 169 µs | 81.8 µs | 629 µs | 1.59 ms |  | 3.7x faster |
| `nearest/n10000/shuffled` | 417 µs | 105 µs | 3.22 ms | 3.29 ms |  | 7.7x faster |
| `previous/n10000/sorted` | 176 µs | 86 µs | 670 µs |  |  | 3.8x faster |
| `previous/n10000/shuffled` | 377 µs | 97.1 µs | 3.32 ms |  |  | 8.8x faster |
| `bspline1/n10000/sorted` | 274 µs | 93.7 µs | 320 µs |  |  | 1.2x faster |
| `bspline1/n10000/shuffled` | 700 µs | 152 µs | 53.5 ms |  |  | 76x faster |
| `bspline2/n10000/sorted` | 356 µs | 100 µs | 819 µs |  |  | 2.3x faster |
| `bspline2/n10000/shuffled` | 977 µs | 233 µs | 52.2 ms |  |  | 53x faster |
| `bspline3/n10000/sorted` | 402 µs | 125 µs | 1.39 ms |  |  | 3.4x faster |
| `bspline3/n10000/shuffled` | 1.23 ms | 278 µs | 53 ms |  |  | 43x faster |
| `cubic_spline/n10000/sorted` | 168 µs | 80 µs | 402 µs | 154 ms | 322 µs | 1.9x faster |
| `cubic_spline/n10000/shuffled` | 239 µs | 98.8 µs | 3.7 ms | 155 ms | 366 µs | 1.5x faster |
| `pchip/n10000/sorted` | 168 µs | 81.7 µs | 397 µs | 130 ms |  | 2.4x faster |
| `pchip/n10000/shuffled` | 239 µs | 104 µs | 3.69 ms | 132 ms |  | 15x faster |
| `akima/n10000/sorted` | 167 µs | 83.6 µs | 402 µs | 158 ms |  | 2.4x faster |
| `akima/n10000/shuffled` | 239 µs | 98.6 µs | 3.69 ms | 160 ms |  | 15x faster |
| `makima/n10000/sorted` | 163 µs | 79.4 µs | 400 µs |  |  | 2.5x faster |
| `makima/n10000/shuffled` | 239 µs | 96.3 µs | 3.71 ms |  |  | 16x faster |
| `linear/n100/sorted-uneven` | 1.33 µs | 1.4 µs | 8.18 µs | 3.73 µs | 4.13 µs | 2.8x faster |
| `linear/n100/shuffled-uneven` | 3.11 µs | 3.15 µs | 13.8 µs | 6.64 µs | 6.75 µs | 2.1x faster |
| `cubic_spline/n100/sorted-uneven` | 1.53 µs | 1.57 µs | 8.06 µs | 19.1 µs | 6.62 µs | 4.3x faster |
| `cubic_spline/n100/shuffled-uneven` | 3.21 µs | 3.4 µs | 8.77 µs | 19.8 µs | 4.84 µs | 1.5x faster |
| `pchip/n100/sorted-uneven` | 1.46 µs | 1.59 µs | 8.02 µs | 17.2 µs |  | 5.5x faster |
| `pchip/n100/shuffled-uneven` | 3.18 µs | 3.36 µs | 8.71 µs | 17.4 µs |  | 2.7x faster |
| `linear/n10000/sorted-uneven` | 224 µs | 89.1 µs | 281 µs | 1.29 ms | 982 µs | 1.3x faster |
| `linear/n10000/shuffled-uneven` | 686 µs | 179 µs | 2.53 ms | 3.09 ms | 1.93 ms | 2.8x faster |
| `cubic_spline/n10000/sorted-uneven` | 241 µs | 98.3 µs | 481 µs | 168 ms | 1.05 ms | 2x faster |
| `cubic_spline/n10000/shuffled-uneven` | 804 µs | 189 µs | 3.69 ms | 164 ms | 2 ms | 2.5x faster |
| `pchip/n10000/sorted-uneven` | 241 µs | 94.2 µs | 480 µs | 132 ms |  | 2x faster |
| `pchip/n10000/shuffled-uneven` | 800 µs | 189 µs | 3.67 ms | 131 ms |  | 4.6x faster |

### Construction

| case | interpx | interpx, parallel | SciPy | scirs2-interpolate | ndarray-interp | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/n10/sorted` | 87.2 ns | 95.1 ns | 9.64 µs | 89.4 ns | 41.5 ns | **2.1x slower** |
| `linear/n10/shuffled` | 93.9 ns | 101 ns | 9.67 µs | 87.4 ns | 45.5 ns | **2.1x slower** |
| `nearest/n10/sorted` | 82.5 ns | 83 ns | 10.7 µs | 89.5 ns |  | 1.1x faster |
| `nearest/n10/shuffled` | 85 ns | 83 ns | 10.6 µs | 87.6 ns |  | 1x faster |
| `previous/n10/sorted` | 83 ns | 82.6 ns | 10.2 µs |  |  | 1.2e+02x faster |
| `previous/n10/shuffled` | 81.7 ns | 82.2 ns | 10.1 µs |  |  | 1.2e+02x faster |
| `bspline1/n10/sorted` | 106 ns | 114 ns | 17.4 µs |  |  | 1.6e+02x faster |
| `bspline1/n10/shuffled` | 113 ns | 113 ns | 17.4 µs |  |  | 1.5e+02x faster |
| `bspline2/n10/sorted` | 510 ns | 518 ns | 40.2 µs |  |  | 79x faster |
| `bspline2/n10/shuffled` | 518 ns | 510 ns | 39.8 µs |  |  | 77x faster |
| `bspline3/n10/sorted` | 572 ns | 558 ns | 39.1 µs |  |  | 68x faster |
| `bspline3/n10/shuffled` | 560 ns | 562 ns | 38.6 µs |  |  | 69x faster |
| `cubic_spline/n10/sorted` | 254 ns | 253 ns | 93 µs | 380 ns | 488 ns | 1.5x faster |
| `cubic_spline/n10/shuffled` | 255 ns | 275 ns | 94.1 µs | 410 ns | 472 ns | 1.6x faster |
| `pchip/n10/sorted` | 173 ns | 192 ns | 109 µs | 209 ns |  | 1.2x faster |
| `pchip/n10/shuffled` | 189 ns | 191 ns | 109 µs | 208 ns |  | 1.1x faster |
| `akima/n10/sorted` | 194 ns | 196 ns | 81.5 µs | 197 ns |  | 1x faster |
| `akima/n10/shuffled` | 196 ns | 179 ns | 76.6 µs | 199 ns |  | 1x faster |
| `makima/n10/sorted` | 204 ns | 204 ns | 78.3 µs |  |  | 3.8e+02x faster |
| `makima/n10/shuffled` | 205 ns | 207 ns | 85.3 µs |  |  | 4.2e+02x faster |
| `linear/n50/sorted` | 161 ns | 163 ns | 9.42 µs | 133 ns | 110 ns | **1.5x slower** |
| `linear/n50/shuffled` | 164 ns | 166 ns | 9.84 µs | 128 ns | 104 ns | **1.6x slower** |
| `nearest/n50/sorted` | 193 ns | 119 ns | 10.5 µs | 130 ns |  | **1.5x slower** |
| `nearest/n50/shuffled` | 120 ns | 118 ns | 10.6 µs | 130 ns |  | 1.1x faster |
| `previous/n50/sorted` | 124 ns | 128 ns | 9.6 µs |  |  | 77x faster |
| `previous/n50/shuffled` | 116 ns | 128 ns | 10.4 µs |  |  | 90x faster |
| `bspline1/n50/sorted` | 153 ns | 158 ns | 17.5 µs |  |  | 1.1e+02x faster |
| `bspline1/n50/shuffled` | 160 ns | 161 ns | 17.2 µs |  |  | 1.1e+02x faster |
| `bspline2/n50/sorted` | 1.83 µs | 2 µs | 44.1 µs |  |  | 24x faster |
| `bspline2/n50/shuffled` | 1.86 µs | 2 µs | 43.6 µs |  |  | 24x faster |
| `bspline3/n50/sorted` | 2.24 µs | 2.23 µs | 43.5 µs |  |  | 19x faster |
| `bspline3/n50/shuffled` | 2.23 µs | 2.23 µs | 67 µs |  |  | 30x faster |
| `cubic_spline/n50/sorted` | 917 ns | 920 ns | 96.4 µs | 979 ns | 1.22 µs | 1.1x faster |
| `cubic_spline/n50/shuffled` | 918 ns | 917 ns | 95.2 µs | 945 ns | 1.23 µs | 1x faster |
| `pchip/n50/sorted` | 481 ns | 535 ns | 102 µs | 371 ns |  | **1.3x slower** |
| `pchip/n50/shuffled` | 532 ns | 537 ns | 102 µs | 390 ns |  | **1.4x slower** |
| `akima/n50/sorted` | 547 ns | 555 ns | 83.4 µs | 369 ns |  | **1.5x slower** |
| `akima/n50/shuffled` | 552 ns | 559 ns | 83.6 µs | 366 ns |  | **1.5x slower** |
| `makima/n50/sorted` | 601 ns | 602 ns | 87.3 µs |  |  | 1.5e+02x faster |
| `makima/n50/shuffled` | 595 ns | 600 ns | 86.9 µs |  |  | 1.5e+02x faster |
| `linear/n100/sorted` | 270 ns | 269 ns | 9.5 µs | 195 ns | 183 ns | **1.5x slower** |
| `linear/n100/shuffled` | 271 ns | 273 ns | 10 µs | 209 ns | 183 ns | **1.5x slower** |
| `nearest/n100/sorted` | 183 ns | 177 ns | 11.1 µs | 207 ns |  | 1.1x faster |
| `nearest/n100/shuffled` | 184 ns | 179 ns | 11.4 µs | 214 ns |  | 1.2x faster |
| `previous/n100/sorted` | 300 ns | 288 ns | 11.3 µs |  |  | 38x faster |
| `previous/n100/shuffled` | 288 ns | 180 ns | 11.2 µs |  |  | 39x faster |
| `bspline1/n100/sorted` | 231 ns | 222 ns | 28.9 µs |  |  | 1.2e+02x faster |
| `bspline1/n100/shuffled` | 224 ns | 218 ns | 18.2 µs |  |  | 81x faster |
| `bspline2/n100/sorted` | 3.83 µs | 3.82 µs | 47.3 µs |  |  | 12x faster |
| `bspline2/n100/shuffled` | 3.83 µs | 3.82 µs | 48 µs |  |  | 13x faster |
| `bspline3/n100/sorted` | 4.48 µs | 4.38 µs | 49 µs |  |  | 11x faster |
| `bspline3/n100/shuffled` | 4.29 µs | 4.3 µs | 48.6 µs |  |  | 11x faster |
| `cubic_spline/n100/sorted` | 1.73 µs | 1.75 µs | 96.4 µs | 1.66 µs | 2.07 µs | **1x slower** |
| `cubic_spline/n100/shuffled` | 1.74 µs | 1.73 µs | 97 µs | 1.66 µs | 2.79 µs | **1x slower** |
| `pchip/n100/sorted` | 981 ns | 981 ns | 116 µs | 596 ns |  | **1.6x slower** |
| `pchip/n100/shuffled` | 969 ns | 973 ns | 114 µs | 594 ns |  | **1.6x slower** |
| `akima/n100/sorted` | 967 ns | 980 ns | 84 µs | 589 ns |  | **1.6x slower** |
| `akima/n100/shuffled` | 973 ns | 1.73 µs | 84.2 µs | 587 ns |  | **1.7x slower** |
| `makima/n100/sorted` | 1.07 µs | 1.07 µs | 88.1 µs |  |  | 82x faster |
| `makima/n100/shuffled` | 1.07 µs | 1.07 µs | 87 µs |  |  | 81x faster |
| `linear/n1000/sorted` | 2.13 µs | 2.13 µs | 13.7 µs | 1.28 µs | 1.51 µs | **1.7x slower** |
| `linear/n1000/shuffled` | 2.17 µs | 2.13 µs | 14.2 µs | 1.28 µs | 1.52 µs | **1.7x slower** |
| `nearest/n1000/sorted` | 1.09 µs | 1.09 µs | 16 µs | 1.27 µs |  | 1.2x faster |
| `nearest/n1000/shuffled` | 1.07 µs | 1.09 µs | 14.9 µs | 1.28 µs |  | 1.2x faster |
| `previous/n1000/sorted` | 1.09 µs | 1.08 µs | 18 µs |  |  | 17x faster |
| `previous/n1000/shuffled` | 1.09 µs | 1.1 µs | 17.6 µs |  |  | 16x faster |
| `bspline1/n1000/sorted` | 1.15 µs | 1.14 µs | 18.6 µs |  |  | 16x faster |
| `bspline1/n1000/shuffled` | 1.14 µs | 1.13 µs | 18.7 µs |  |  | 16x faster |
| `bspline2/n1000/sorted` | 36.5 µs | 36.5 µs | 119 µs |  |  | 3.2x faster |
| `bspline2/n1000/shuffled` | 36.6 µs | 49.6 µs | 119 µs |  |  | 3.2x faster |
| `bspline3/n1000/sorted` | 43 µs | 42.2 µs | 226 µs |  |  | 5.3x faster |
| `bspline3/n1000/shuffled` | 41.9 µs | 44.1 µs | 138 µs |  |  | 3.3x faster |
| `cubic_spline/n1000/sorted` | 17 µs | 16.9 µs | 124 µs | 14.6 µs | 17.3 µs | **1.2x slower** |
| `cubic_spline/n1000/shuffled` | 16.9 µs | 16.9 µs | 115 µs | 14.6 µs | 17.3 µs | **1.2x slower** |
| `pchip/n1000/sorted` | 9.11 µs | 9.15 µs | 120 µs | 4.59 µs |  | **2x slower** |
| `pchip/n1000/shuffled` | 9.17 µs | 9.08 µs | 130 µs | 4.61 µs |  | **2x slower** |
| `akima/n1000/sorted` | 9.18 µs | 9.07 µs | 101 µs | 4.85 µs |  | **1.9x slower** |
| `akima/n1000/shuffled` | 9.17 µs | 9.13 µs | 105 µs | 4.83 µs |  | **1.9x slower** |
| `makima/n1000/sorted` | 9.96 µs | 9.95 µs | 108 µs |  |  | 11x faster |
| `makima/n1000/shuffled` | 10.2 µs | 9.97 µs | 107 µs |  |  | 11x faster |
| `linear/n10000/sorted` | 21.7 µs | 21.8 µs | 48.3 µs | 12.9 µs | 22.1 µs | **1.7x slower** |
| `linear/n10000/shuffled` | 22.9 µs | 21.5 µs | 47.1 µs | 12.9 µs | 16.8 µs | **1.8x slower** |
| `nearest/n10000/sorted` | 11.9 µs | 11.9 µs | 57.3 µs | 12.9 µs |  | 1.1x faster |
| `nearest/n10000/shuffled` | 11.8 µs | 11.8 µs | 58.9 µs | 12.9 µs |  | 1.1x faster |
| `previous/n10000/sorted` | 11.9 µs | 11.9 µs | 80 µs |  |  | 6.7x faster |
| `previous/n10000/shuffled` | 11.9 µs | 11.8 µs | 81.6 µs |  |  | 6.9x faster |
| `bspline1/n10000/sorted` | 12 µs | 13.1 µs | 26.3 µs |  |  | 2.2x faster |
| `bspline1/n10000/shuffled` | 11.9 µs | 12 µs | 25.6 µs |  |  | 2.2x faster |
| `bspline2/n10000/sorted` | 523 µs | 531 µs | 822 µs |  |  | 1.6x faster |
| `bspline2/n10000/shuffled` | 534 µs | 520 µs | 822 µs |  |  | 1.5x faster |
| `bspline3/n10000/sorted` | 514 µs | 568 µs | 1.03 ms |  |  | 2x faster |
| `bspline3/n10000/shuffled` | 516 µs | 581 µs | 1.02 ms |  |  | 2x faster |
| `cubic_spline/n10000/sorted` | 240 µs | 176 µs | 663 µs | 145 µs | 168 µs | **1.7x slower** |
| `cubic_spline/n10000/shuffled` | 169 µs | 184 µs | 788 µs | 147 µs | 183 µs | **1.2x slower** |
| `pchip/n10000/sorted` | 98.2 µs | 99.7 µs | 257 µs | 47.3 µs |  | **2.1x slower** |
| `pchip/n10000/shuffled` | 134 µs | 92.7 µs | 260 µs | 46.9 µs |  | **2.9x slower** |
| `akima/n10000/sorted` | 143 µs | 108 µs | 253 µs | 48.7 µs |  | **2.9x slower** |
| `akima/n10000/shuffled` | 101 µs | 147 µs | 254 µs | 45.9 µs |  | **2.2x slower** |
| `makima/n10000/sorted` | 101 µs | 114 µs | 716 µs |  |  | 7.1x faster |
| `makima/n10000/shuffled` | 137 µs | 103 µs | 756 µs |  |  | 5.5x faster |
| `linear/n100/sorted-uneven` | 267 ns | 270 ns | 10.2 µs | 208 ns | 183 ns | **1.5x slower** |
| `linear/n100/shuffled-uneven` | 271 ns | 266 ns | 10.1 µs | 339 ns | 257 ns | **1.1x slower** |
| `cubic_spline/n100/sorted-uneven` | 1.73 µs | 1.74 µs | 96.8 µs | 1.69 µs | 2.06 µs | **1x slower** |
| `cubic_spline/n100/shuffled-uneven` | 1.75 µs | 1.74 µs | 96.6 µs | 1.66 µs | 2.08 µs | **1.1x slower** |
| `pchip/n100/sorted-uneven` | 968 ns | 967 ns | 112 µs | 594 ns |  | **1.6x slower** |
| `pchip/n100/shuffled-uneven` | 980 ns | 968 ns | 112 µs | 594 ns |  | **1.6x slower** |
| `linear/n10000/sorted-uneven` | 21.5 µs | 21.7 µs | 51.9 µs | 12.9 µs | 23.1 µs | **1.7x slower** |
| `linear/n10000/shuffled-uneven` | 21.5 µs | 21.6 µs | 48.1 µs | 12.9 µs | 19.6 µs | **1.7x slower** |
| `cubic_spline/n10000/sorted-uneven` | 174 µs | 172 µs | 786 µs | 146 µs | 184 µs | **1.2x slower** |
| `cubic_spline/n10000/shuffled-uneven` | 170 µs | 191 µs | 774 µs | 146 µs | 184 µs | **1.2x slower** |
| `pchip/n10000/sorted-uneven` | 124 µs | 92.6 µs | 255 µs | 47.1 µs |  | **2.6x slower** |
| `pchip/n10000/shuffled-uneven` | 138 µs | 91.8 µs | 257 µs | 46.8 µs |  | **3x slower** |

## `PPoly`

A cubic of 55 intervals, evaluated at points in no order.

### Evaluation

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `intervals55/eval10` | 80.2 ns | 80 ns | 3.68 µs | 46x faster |
| `intervals55/eval100` | 703 ns | 783 ns | 4.79 µs | 6.8x faster |
| `intervals55/eval1000` | 7.12 µs | 7.3 µs | 16.8 µs | 2.4x faster |
| `intervals55/eval10000` | 72.9 µs | 67.8 µs | 362 µs | 5x faster |

### Construction

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `intervals55/eval10` | 174 ns | 179 ns | 5.7 µs | 33x faster |
| `intervals55/eval100` | 179 ns | 173 ns | 5.66 µs | 32x faster |
| `intervals55/eval1000` | 172 ns | 173 ns | 5.65 µs | 33x faster |
| `intervals55/eval10000` | 173 ns | 173 ns | 5.2 µs | 30x faster |

## `BSpline`

A cubic B-spline of 55 intervals, evaluated at points in no order.

### Evaluation

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `intervals55/eval10` | 160 ns | 161 ns | 3 µs | 19x faster |
| `intervals55/eval100` | 1.56 µs | 1.64 µs | 6.42 µs | 4.1x faster |
| `intervals55/eval1000` | 15.6 µs | 15.8 µs | 41.8 µs | 2.7x faster |
| `intervals55/eval10000` | 156 µs | 64.9 µs | 398 µs | 2.6x faster |

### Construction

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `intervals55/eval10` | 175 ns | 179 ns | 9.91 µs | 57x faster |
| `intervals55/eval100` | 175 ns | 162 ns | 10 µs | 57x faster |
| `intervals55/eval1000` | 175 ns | 175 ns | 10.1 µs | 58x faster |
| `intervals55/eval10000` | 160 ns | 175 ns | 10.1 µs | 63x faster |

## `BSplineCurve`

A cubic curve through `n` points of the plane, evaluated at `4 n` parameters in order.

### Evaluation

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `degree3/n100` | 4.56 µs | 4.72 µs | 17.4 µs | 3.8x faster |
| `degree3/n10000` | 497 µs | 157 µs | 1.36 ms | 2.7x faster |

### Construction

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `degree3/n100` | 4.3 µs | 4.66 µs | 50.3 µs | 12x faster |
| `degree3/n10000` | 680 µs | 690 µs | 1.46 ms | 2.1x faster |

## Regular grid

Points in no order inside the grid. `maxN`: the first axis has `N` points and each axis half the one before (SciPy's suite); `eachN`: every axis has `N` points; `flipped`: every axis given from its high end.

### Evaluation

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | ndarray-interp | interpn | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/2d/max10/eval10/ascending` | 150 ns | 153 ns | 150 ns | 20.1 µs | 930 ns | 166 ns | 104 ns | **1.4x slower** |
| `linear/2d/max10/eval10/flipped` | 150 ns | 147 ns | 150 ns | 20.9 µs |  |  |  | 1.4e+02x faster |
| `linear/2d/max10/eval100/ascending` | 1.38 µs | 1.44 µs | 1.39 µs | 23.6 µs | 9.96 µs | 1.75 µs | 1.02 µs | **1.4x slower** |
| `linear/2d/max10/eval100/flipped` | 1.37 µs | 1.44 µs | 1.38 µs | 23.7 µs |  |  |  | 17x faster |
| `linear/2d/max10/eval1000/ascending` | 13.7 µs | 13.8 µs | 13.7 µs | 43.9 µs | 93.3 µs | 15.4 µs | 10.1 µs | **1.3x slower** |
| `linear/2d/max10/eval1000/flipped` | 13.6 µs | 13.8 µs | 13.6 µs | 43.6 µs |  |  |  | 3.2x faster |
| `linear/2d/max10/eval10000/ascending` | 137 µs | 59.4 µs | 137 µs | 395 µs | 977 µs | 239 µs | 101 µs | **1.4x slower** |
| `linear/2d/max10/eval10000/flipped` | 137 µs | 63.2 µs | 137 µs | 393 µs |  |  |  | 2.9x faster |
| `linear/2d/max40/eval10/ascending` | 172 ns | 175 ns | 173 ns | 21.2 µs | 1.19 µs | 172 ns | 147 ns | **1.2x slower** |
| `linear/2d/max40/eval10/flipped` | 172 ns | 176 ns | 173 ns | 20.9 µs |  |  |  | 1.2e+02x faster |
| `linear/2d/max40/eval100/ascending` | 1.66 µs | 1.74 µs | 1.55 µs | 24.3 µs | 11 µs | 1.59 µs | 1.38 µs | **1.2x slower** |
| `linear/2d/max40/eval100/flipped` | 1.66 µs | 1.73 µs | 1.62 µs | 24 µs |  |  |  | 15x faster |
| `linear/2d/max40/eval1000/ascending` | 16.6 µs | 16.8 µs | 16.6 µs | 49.8 µs | 126 µs | 16.8 µs | 13.2 µs | **1.3x slower** |
| `linear/2d/max40/eval1000/flipped` | 15.1 µs | 16.7 µs | 16.5 µs | 46.1 µs |  |  |  | 3.1x faster |
| `linear/2d/max40/eval10000/ascending` | 156 µs | 70.6 µs | 166 µs | 624 µs | 1.29 ms | 246 µs | 145 µs | **1.1x slower** |
| `linear/2d/max40/eval10000/flipped` | 166 µs | 69.6 µs | 166 µs | 623 µs |  |  |  | 3.8x faster |
| `linear/2d/max200/eval10/ascending` | 211 ns | 196 ns | 212 ns | 21.1 µs | 2.04 µs | 189 ns | 192 ns | **1.1x slower** |
| `linear/2d/max200/eval10/flipped` | 211 ns | 215 ns | 212 ns | 19.7 µs |  |  |  | 94x faster |
| `linear/2d/max200/eval100/ascending` | 1.89 µs | 2.14 µs | 2.06 µs | 23 µs | 22.6 µs | 2.18 µs | 1.89 µs | **1x slower** |
| `linear/2d/max200/eval100/flipped` | 2.07 µs | 3.61 µs | 2.06 µs | 25.1 µs |  |  |  | 12x faster |
| `linear/2d/max200/eval1000/ascending` | 20.6 µs | 20.8 µs | 20.7 µs | 76.2 µs | 233 µs | 24.4 µs | 18.8 µs | **1.1x slower** |
| `linear/2d/max200/eval1000/flipped` | 20.8 µs | 20.8 µs | 20.6 µs | 78 µs |  |  |  | 3.8x faster |
| `linear/2d/max200/eval10000/ascending` | 206 µs | 80.7 µs | 195 µs | 883 µs | 2.43 ms | 490 µs | 189 µs | **1.1x slower** |
| `linear/2d/max200/eval10000/flipped` | 207 µs | 78.4 µs | 207 µs | 880 µs |  |  |  | 4.3x faster |
| `linear/3d/max10/eval10/ascending` | 200 ns | 225 ns | 218 ns | 62.5 µs | 4.09 µs |  | 156 ns | **1.3x slower** |
| `linear/3d/max10/eval10/flipped` | 217 ns | 225 ns | 218 ns | 62 µs |  |  |  | 2.9e+02x faster |
| `linear/3d/max10/eval100/ascending` | 2.09 µs | 2.16 µs | 2.08 µs | 68.5 µs | 38.6 µs |  | 1.53 µs | **1.4x slower** |
| `linear/3d/max10/eval100/flipped` | 1.91 µs | 2.16 µs | 2.09 µs | 67.8 µs |  |  |  | 35x faster |
| `linear/3d/max10/eval1000/ascending` | 20.7 µs | 21 µs | 20.9 µs | 126 µs | 401 µs |  | 15.3 µs | **1.4x slower** |
| `linear/3d/max10/eval1000/flipped` | 20.6 µs | 21.1 µs | 20.9 µs | 119 µs |  |  |  | 5.8x faster |
| `linear/3d/max10/eval10000/ascending` | 209 µs | 73.3 µs | 209 µs | 1.24 ms | 3.82 ms |  | 155 µs | **1.3x slower** |
| `linear/3d/max10/eval10000/flipped` | 209 µs | 73 µs | 197 µs | 1.18 ms |  |  |  | 5.6x faster |
| `linear/3d/max40/eval10/ascending` | 240 ns | 247 ns | 221 ns | 57.2 µs | 4.32 µs |  | 199 ns | **1.2x slower** |
| `linear/3d/max40/eval10/flipped` | 243 ns | 246 ns | 241 ns | 60.9 µs |  |  |  | 2.5e+02x faster |
| `linear/3d/max40/eval100/ascending` | 2.34 µs | 2.4 µs | 2.35 µs | 69.9 µs | 40.3 µs |  | 2.17 µs | **1.1x slower** |
| `linear/3d/max40/eval100/flipped` | 2.45 µs | 2.41 µs | 2.33 µs | 70.7 µs |  |  |  | 29x faster |
| `linear/3d/max40/eval1000/ascending` | 23.3 µs | 23.6 µs | 23.5 µs | 144 µs | 406 µs |  | 21.7 µs | **1.1x slower** |
| `linear/3d/max40/eval1000/flipped` | 23.3 µs | 23.5 µs | 21.9 µs | 143 µs |  |  |  | 6.1x faster |
| `linear/3d/max40/eval10000/ascending` | 215 µs | 81.6 µs | 233 µs | 1.61 ms | 4.14 ms |  | 219 µs | 1x faster |
| `linear/3d/max40/eval10000/flipped` | 233 µs | 84.4 µs | 233 µs | 1.64 ms |  |  |  | 7x faster |
| `linear/3d/max200/eval10/ascending` | 532 ns | 309 ns | 303 ns | 63.2 µs | 7.78 µs |  | 287 ns | **1.9x slower** |
| `linear/3d/max200/eval10/flipped` | 318 ns | 324 ns | 313 ns | 66.6 µs |  |  |  | 2.1e+02x faster |
| `linear/3d/max200/eval100/ascending` | 2.97 µs | 3.06 µs | 2.98 µs | 71.1 µs | 58.3 µs |  | 2.87 µs | **1x slower** |
| `linear/3d/max200/eval100/flipped` | 2.97 µs | 3.06 µs | 2.98 µs | 70.3 µs |  |  |  | 24x faster |
| `linear/3d/max200/eval1000/ascending` | 29.7 µs | 30 µs | 29.7 µs | 193 µs | 578 µs |  | 28.9 µs | **1x slower** |
| `linear/3d/max200/eval1000/flipped` | 29.9 µs | 29.9 µs | 29.7 µs | 208 µs |  |  |  | 6.9x faster |
| `linear/3d/max200/eval10000/ascending` | 297 µs | 102 µs | 297 µs | 1.98 ms | 5.83 ms |  | 289 µs | **1x slower** |
| `linear/3d/max200/eval10000/flipped` | 297 µs | 99.8 µs | 297 µs | 2.02 ms |  |  |  | 6.8x faster |
| `linear/4d/max10/eval10/ascending` | 285 ns | 288 ns | 283 ns | 104 µs |  |  |  | 3.6e+02x faster |
| `linear/4d/max10/eval10/flipped` | 286 ns | 503 ns | 284 ns | 105 µs |  |  |  | 3.7e+02x faster |
| `linear/4d/max10/eval100/ascending` | 2.74 µs | 2.72 µs | 2.74 µs | 116 µs |  |  |  | 42x faster |
| `linear/4d/max10/eval100/flipped` | 2.64 µs | 2.71 µs | 2.68 µs | 116 µs |  |  |  | 44x faster |
| `linear/4d/max10/eval1000/ascending` | 27.1 µs | 25.7 µs | 26.6 µs | 228 µs |  |  |  | 8.4x faster |
| `linear/4d/max10/eval1000/flipped` | 26.7 µs | 26.8 µs | 26.7 µs | 227 µs |  |  |  | 8.5x faster |
| `linear/4d/max10/eval10000/ascending` | 265 µs | 94 µs | 270 µs | 2.02 ms |  |  |  | 7.6x faster |
| `linear/4d/max10/eval10000/flipped` | 495 µs | 90 µs | 277 µs | 3.41 ms |  |  |  | 6.9x faster |
| `linear/4d/max40/eval10/ascending` | 338 ns | 345 ns | 333 ns | 104 µs | 6.48 µs |  | 296 ns | **1.1x slower** |
| `linear/4d/max40/eval10/flipped` | 337 ns | 317 ns | 337 ns | 104 µs |  |  |  | 3.1e+02x faster |
| `linear/4d/max40/eval100/ascending` | 3.22 µs | 3.31 µs | 3.26 µs | 116 µs | 59 µs |  | 2.96 µs | **1.1x slower** |
| `linear/4d/max40/eval100/flipped` | 3.28 µs | 3.31 µs | 3.23 µs | 110 µs |  |  |  | 34x faster |
| `linear/4d/max40/eval1000/ascending` | 32.3 µs | 32.7 µs | 32.7 µs | 266 µs | 643 µs |  | 29.6 µs | **1.1x slower** |
| `linear/4d/max40/eval1000/flipped` | 32.3 µs | 32.4 µs | 32.9 µs | 272 µs |  |  |  | 8.4x faster |
| `linear/4d/max40/eval10000/ascending` | 327 µs | 105 µs | 323 µs | 2.47 ms | 6.51 ms |  | 297 µs | **1.1x slower** |
| `linear/4d/max40/eval10000/flipped` | 319 µs | 109 µs | 324 µs | 2.32 ms |  |  |  | 7.3x faster |
| `linear/4d/max200/eval10/ascending` | 414 ns | 422 ns | 420 ns | 105 µs | 8.6 µs |  | 390 ns | **1.1x slower** |
| `linear/4d/max200/eval10/flipped` | 421 ns | 422 ns | 420 ns | 105 µs |  |  |  | 2.5e+02x faster |
| `linear/4d/max200/eval100/ascending` | 4.1 µs | 4.32 µs | 4.14 µs | 122 µs | 93.7 µs |  | 4.07 µs | **1x slower** |
| `linear/4d/max200/eval100/flipped` | 4.17 µs | 4.2 µs | 4.15 µs | 125 µs |  |  |  | 30x faster |
| `linear/4d/max200/eval1000/ascending` | 41 µs | 41.4 µs | 41.2 µs | 343 µs | 801 µs |  | 39.7 µs | **1x slower** |
| `linear/4d/max200/eval1000/flipped` | 41.4 µs | 41.5 µs | 41.2 µs | 332 µs |  |  |  | 8x faster |
| `linear/4d/max200/eval10000/ascending` | 414 µs | 130 µs | 411 µs | 3.6 ms | 8.34 ms |  | 400 µs | **1x slower** |
| `linear/4d/max200/eval10000/flipped` | 409 µs | 127 µs | 413 µs | 3.05 ms |  |  |  | 7.5x faster |
| `nearest/2d/each40/eval10000/ascending` | 166 µs | 68 µs | 166 µs | 706 µs | 1.7 ms |  | 150 µs | **1.1x slower** |
| `nearest/3d/each40/eval10000/ascending` | 233 µs | 89.5 µs | 236 µs | 1.04 ms | 2.04 ms |  | 229 µs | **1x slower** |
| `pchip/2d/each40/eval10000/ascending` | 556 µs | 145 µs | 548 µs | 1.41 s |  |  |  | 2.5e+03x faster |
| `pchip/3d/each40/eval10000/ascending` | 1.62 ms | 619 µs | 1.61 ms | 3.74 s |  |  |  | 2.3e+03x faster |
| `cubic/2d/each40/eval10000/ascending` | 343 µs | 131 µs | 332 µs | 1.4 ms |  |  | 466 µs | 1.4x faster |
| `cubic/3d/each40/eval10000/ascending` | 729 µs | 184 µs | 726 µs | 3.24 ms |  |  | 867 µs | 1.2x faster |
| `cubic/2d/each10/eval1000/ascending` | 29.1 µs | 29.3 µs | 29.1 µs | 145 µs |  |  | 39.4 µs | 1.4x faster |
| `cubic/2d/each40/eval1000/ascending` | 33.1 µs | 33.2 µs | 33.1 µs | 157 µs |  |  | 44.6 µs | 1.3x faster |
| `cubic/2d/each100/eval1000/ascending` | 35.4 µs | 35.7 µs | 35.4 µs | 184 µs |  |  | 47.2 µs | 1.3x faster |
| `cubic/2d/each200/eval1000/ascending` | 38.3 µs | 38.1 µs | 40.5 µs | 220 µs |  |  | 50.2 µs | 1.3x faster |
| `cubic/2d/each400/eval1000/ascending` | 42.3 µs | 42 µs | 42.3 µs | 301 µs |  |  | 54.5 µs | 1.3x faster |
| `quintic/2d/each10/eval1000/ascending` | 48.6 µs | 48.6 µs | 49.1 µs | 251 µs |  |  |  | 5.2x faster |
| `quintic/2d/each40/eval1000/ascending` | 52.8 µs | 52.9 µs | 52.7 µs | 265 µs |  |  |  | 5x faster |
| `linear/2d/each400/eval10000/ascending` | 241 µs | 96.7 µs | 242 µs | 1.07 ms | 4.66 ms | 596 µs | 241 µs | 1x faster |
| `linear/2d/each800/eval10000/ascending` | 262 µs | 106 µs | 262 µs | 1.18 ms | 8.41 ms | 642 µs | 266 µs | 1x faster |
| `nearest/2d/auckland-ba32-0401/heldout/ascending` | 376 ms | 69.1 ms | 389 ms | 1.41 s | 172 s |  | 781 ms | 2.1x faster |
| `linear/2d/auckland-ba32-0401/heldout/ascending` | 426 ms | 74.6 ms | 452 ms | 1.21 s | 81.6 s | 375 ms | 888 ms | **1.1x slower** |
| `cubic/2d/auckland-ba32-0401/heldout/ascending` | 803 ms | 137 ms | 859 ms | 36.3 s |  |  | 1.83 s | 2.3x faster |
| `quintic/2d/auckland-ba32-0401/heldout/ascending` | 1.35 s | 235 ms | 1.32 s | 39 s |  |  |  | 29x faster |
| `pchip/2d/auckland-ba32-0401/heldout/ascending` | 1.1 s | 184 ms | 1.09 s |  |  |  |  |  |

### Construction

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | ndarray-interp | interpn | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/2d/max10/eval10/ascending` | 140 ns | 138 ns | 133 ns | 5.96 µs | 145 ns | 85.5 ns | 0.78 ns | **1.8e+02x slower** |
| `linear/2d/max10/eval10/flipped` | 172 ns | 168 ns | 192 ns | 12.5 µs |  |  |  | 73x faster |
| `linear/2d/max10/eval100/ascending` | 128 ns | 135 ns | 142 ns | 6.33 µs | 144 ns | 85.6 ns | 0.78 ns | **1.6e+02x slower** |
| `linear/2d/max10/eval100/flipped` | 173 ns | 167 ns | 181 ns | 12.7 µs |  |  |  | 74x faster |
| `linear/2d/max10/eval1000/ascending` | 140 ns | 135 ns | 141 ns | 6.32 µs | 145 ns | 85.4 ns | 0.68 ns | **2.1e+02x slower** |
| `linear/2d/max10/eval1000/flipped` | 167 ns | 166 ns | 173 ns | 12.1 µs |  |  |  | 73x faster |
| `linear/2d/max10/eval10000/ascending` | 135 ns | 134 ns | 140 ns | 6.28 µs | 144 ns | 85.5 ns | 0.78 ns | **1.7e+02x slower** |
| `linear/2d/max10/eval10000/flipped` | 175 ns | 156 ns | 190 ns | 12.8 µs |  |  |  | 73x faster |
| `linear/2d/max40/eval10/ascending` | 395 ns | 400 ns | 859 ns | 6.36 µs | 222 ns | 171 ns | 0.78 ns | **5e+02x slower** |
| `linear/2d/max40/eval10/flipped` | 710 ns | 706 ns | 1.23 µs | 12.7 µs |  |  |  | 18x faster |
| `linear/2d/max40/eval100/ascending` | 405 ns | 397 ns | 858 ns | 6.3 µs | 220 ns | 169 ns | 0.78 ns | **5.2e+02x slower** |
| `linear/2d/max40/eval100/flipped` | 714 ns | 656 ns | 1.14 µs | 12.7 µs |  |  |  | 18x faster |
| `linear/2d/max40/eval1000/ascending` | 397 ns | 388 ns | 792 ns | 5.88 µs | 222 ns | 169 ns | 0.68 ns | **5.8e+02x slower** |
| `linear/2d/max40/eval1000/flipped` | 713 ns | 707 ns | 1.15 µs | 11.7 µs |  |  |  | 16x faster |
| `linear/2d/max40/eval10000/ascending` | 396 ns | 401 ns | 796 ns | 6.35 µs | 218 ns | 169 ns | 0.78 ns | **5.1e+02x slower** |
| `linear/2d/max40/eval10000/flipped` | 715 ns | 697 ns | 1.26 µs | 12.6 µs |  |  |  | 18x faster |
| `linear/2d/max200/eval10/ascending` | 7.41 µs | 6.99 µs | 16.6 µs | 6.47 µs | 3.11 µs | 516 ns | 0.68 ns | **1.1e+04x slower** |
| `linear/2d/max200/eval10/flipped` | 14.6 µs | 15.8 µs | 28.4 µs | 12.2 µs |  |  |  | **1.2x slower** |
| `linear/2d/max200/eval100/ascending` | 6.83 µs | 7.44 µs | 22.2 µs | 6.32 µs | 3.2 µs | 518 ns | 0.78 ns | **8.7e+03x slower** |
| `linear/2d/max200/eval100/flipped` | 16.1 µs | 21.9 µs | 29 µs | 13.2 µs |  |  |  | **1.2x slower** |
| `linear/2d/max200/eval1000/ascending` | 7.52 µs | 7.63 µs | 20.9 µs | 6.4 µs | 3.23 µs | 516 ns | 0.78 ns | **9.6e+03x slower** |
| `linear/2d/max200/eval1000/flipped` | 15.6 µs | 15.8 µs | 24.7 µs | 13.1 µs |  |  |  | **1.2x slower** |
| `linear/2d/max200/eval10000/ascending` | 7.59 µs | 7.53 µs | 4.81 µs | 6.4 µs | 3.23 µs | 516 ns | 0.78 ns | **9.7e+03x slower** |
| `linear/2d/max200/eval10000/flipped` | 15.8 µs | 16 µs | 19.1 µs | 12.9 µs |  |  |  | **1.2x slower** |
| `linear/3d/max10/eval10/ascending` | 203 ns | 194 ns | 233 ns | 8.25 µs | 190 ns |  | 1.17 ns | **1.7e+02x slower** |
| `linear/3d/max10/eval10/flipped` | 334 ns | 323 ns | 352 ns | 17.1 µs |  |  |  | 51x faster |
| `linear/3d/max10/eval100/ascending` | 211 ns | 196 ns | 224 ns | 8.27 µs | 188 ns |  | 0.88 ns | **2.4e+02x slower** |
| `linear/3d/max10/eval100/flipped` | 324 ns | 319 ns | 352 ns | 15.9 µs |  |  |  | 49x faster |
| `linear/3d/max10/eval1000/ascending` | 183 ns | 196 ns | 240 ns | 8.35 µs | 190 ns |  | 0.98 ns | **1.9e+02x slower** |
| `linear/3d/max10/eval1000/flipped` | 331 ns | 317 ns | 359 ns | 17.4 µs |  |  |  | 52x faster |
| `linear/3d/max10/eval10000/ascending` | 195 ns | 197 ns | 236 ns | 8.39 µs | 188 ns |  | 0.98 ns | **2e+02x slower** |
| `linear/3d/max10/eval10000/flipped` | 320 ns | 325 ns | 354 ns | 17.1 µs |  |  |  | 53x faster |
| `linear/3d/max40/eval10/ascending` | 2.83 µs | 3.14 µs | 7.84 µs | 7.61 µs | 1.4 µs |  | 0.88 ns | **3.2e+03x slower** |
| `linear/3d/max40/eval10/flipped` | 8.87 µs | 14.4 µs | 13.2 µs | 16.1 µs |  |  |  | 1.8x faster |
| `linear/3d/max40/eval100/ascending` | 2.93 µs | 3.21 µs | 8.06 µs | 8.25 µs | 1.6 µs |  | 0.98 ns | **3e+03x slower** |
| `linear/3d/max40/eval100/flipped` | 9.03 µs | 10.1 µs | 13.3 µs | 18.5 µs |  |  |  | 2x faster |
| `linear/3d/max40/eval1000/ascending` | 3.27 µs | 3.31 µs | 5.39 µs | 8.31 µs | 1.6 µs |  | 0.98 ns | **3.4e+03x slower** |
| `linear/3d/max40/eval1000/flipped` | 8.84 µs | 9.15 µs | 9.83 µs | 17.2 µs |  |  |  | 1.9x faster |
| `linear/3d/max40/eval10000/ascending` | 3.05 µs | 3.08 µs | 3.2 µs | 8.37 µs | 1.57 µs |  | 0.98 ns | **3.1e+03x slower** |
| `linear/3d/max40/eval10000/flipped` | 9.1 µs | 8.82 µs | 9.91 µs | 17.5 µs |  |  |  | 1.9x faster |
| `linear/3d/max200/eval10/ascending` | 2.36 ms | 1.85 ms | 646 µs | 8.33 µs | 1.84 ms |  | 0.98 ns | **2.4e+06x slower** |
| `linear/3d/max200/eval10/flipped` | 3.66 ms | 3.28 ms | 2.46 ms | 17.8 µs |  |  |  | **2.1e+02x slower** |
| `linear/3d/max200/eval100/ascending` | 1.82 ms | 1.81 ms | 661 µs | 8.31 µs | 1.75 ms |  | 0.98 ns | **1.9e+06x slower** |
| `linear/3d/max200/eval100/flipped` | 3.09 ms | 3.16 ms | 1.95 ms | 17.7 µs |  |  |  | **1.7e+02x slower** |
| `linear/3d/max200/eval1000/ascending` | 1.82 ms | 2.02 ms | 667 µs | 8.22 µs | 1.38 ms |  | 0.98 ns | **1.9e+06x slower** |
| `linear/3d/max200/eval1000/flipped` | 3.02 ms | 3.14 ms | 2 ms | 17.5 µs |  |  |  | **1.7e+02x slower** |
| `linear/3d/max200/eval10000/ascending` | 1.8 ms | 1.82 ms | 683 µs | 8.25 µs | 1.34 ms |  | 0.98 ns | **1.8e+06x slower** |
| `linear/3d/max200/eval10000/flipped` | 3.49 ms | 3 ms | 1.81 ms | 17.4 µs |  |  |  | **2e+02x slower** |
| `linear/4d/max10/eval10/ascending` | 255 ns | 250 ns | 250 ns | 10.7 µs |  |  |  | 42x faster |
| `linear/4d/max10/eval10/flipped` | 388 ns | 698 ns | 401 ns | 19.7 µs |  |  |  | 51x faster |
| `linear/4d/max10/eval100/ascending` | 252 ns | 237 ns | 245 ns | 10.8 µs |  |  |  | 43x faster |
| `linear/4d/max10/eval100/flipped` | 381 ns | 360 ns | 400 ns | 19.6 µs |  |  |  | 52x faster |
| `linear/4d/max10/eval1000/ascending` | 243 ns | 244 ns | 250 ns | 10.4 µs |  |  |  | 43x faster |
| `linear/4d/max10/eval1000/flipped` | 369 ns | 367 ns | 405 ns | 19.6 µs |  |  |  | 53x faster |
| `linear/4d/max10/eval10000/ascending` | 234 ns | 249 ns | 270 ns | 10.9 µs |  |  |  | 47x faster |
| `linear/4d/max10/eval10000/flipped` | 657 ns | 366 ns | 394 ns | 33.7 µs |  |  |  | 51x faster |
| `linear/4d/max40/eval10/ascending` | 14.4 µs | 14.4 µs | 36.6 µs | 10.1 µs | 5.99 µs |  | 1.37 ns | **1.1e+04x slower** |
| `linear/4d/max40/eval10/flipped` | 56.5 µs | 57.1 µs | 69.5 µs | 21.5 µs |  |  |  | **2.6x slower** |
| `linear/4d/max40/eval100/ascending` | 14.3 µs | 14.4 µs | 35.5 µs | 9.51 µs | 6.01 µs |  | 1.37 ns | **1e+04x slower** |
| `linear/4d/max40/eval100/flipped` | 56.5 µs | 57.2 µs | 70.2 µs | 21.6 µs |  |  |  | **2.6x slower** |
| `linear/4d/max40/eval1000/ascending` | 14.3 µs | 14.3 µs | 31.5 µs | 10.2 µs | 6 µs |  | 1.37 ns | **1e+04x slower** |
| `linear/4d/max40/eval1000/flipped` | 57.2 µs | 57.4 µs | 70.5 µs | 21.9 µs |  |  |  | **2.6x slower** |
| `linear/4d/max40/eval10000/ascending` | 14.4 µs | 14.7 µs | 30.8 µs | 10.3 µs | 5.65 µs |  | 1.37 ns | **1.1e+04x slower** |
| `linear/4d/max40/eval10000/flipped` | 56.8 µs | 58.4 µs | 73.9 µs | 21.8 µs |  |  |  | **2.6x slower** |
| `linear/4d/max200/eval10/ascending` | 49.7 ms | 50.4 ms | 16.4 ms | 10.3 µs | 40 ms |  | 1.37 ns | **3.6e+07x slower** |
| `linear/4d/max200/eval10/flipped` | 114 ms | 108 ms | 75.8 ms | 21.8 µs |  |  |  | **5.2e+03x slower** |
| `linear/4d/max200/eval100/ascending` | 47.9 ms | 49.9 ms | 15.9 ms | 10.3 µs | 39.3 ms |  | 1.37 ns | **3.5e+07x slower** |
| `linear/4d/max200/eval100/flipped` | 114 ms | 109 ms | 83.3 ms | 22.2 µs |  |  |  | **5.1e+03x slower** |
| `linear/4d/max200/eval1000/ascending` | 52.8 ms | 47.5 ms | 15.7 ms | 10.4 µs | 39 ms |  | 1.37 ns | **3.9e+07x slower** |
| `linear/4d/max200/eval1000/flipped` | 103 ms | 114 ms | 73.7 ms | 22.2 µs |  |  |  | **4.7e+03x slower** |
| `linear/4d/max200/eval10000/ascending` | 49.7 ms | 50.7 ms | 16.1 ms | 10.3 µs | 40.5 ms |  | 1.37 ns | **3.6e+07x slower** |
| `linear/4d/max200/eval10000/flipped` | 100 ms | 102 ms | 71.3 ms | 21.8 µs |  |  |  | **4.6e+03x slower** |
| `nearest/2d/each40/eval10000/ascending` | 652 ns | 639 ns | 1.78 µs | 6.32 µs | 296 ns |  | 0.78 ns | **8.3e+02x slower** |
| `nearest/3d/each40/eval10000/ascending` | 23.1 µs | 23.4 µs | 75.4 µs | 8.25 µs | 9.69 µs |  | 0.98 ns | **2.4e+04x slower** |
| `pchip/2d/each40/eval10000/ascending` | 7.1 µs | 7.08 µs | 7.08 µs | 7.26 µs |  |  |  | 1x faster |
| `pchip/3d/each40/eval10000/ascending` | 518 µs | 839 µs | 509 µs | 9.46 µs |  |  |  | **55x slower** |
| `cubic/2d/each40/eval10000/ascending` | 9.65 µs | 9.71 µs | 9.44 µs | 2.03 ms |  |  | 0.98 ns | **9.9e+03x slower** |
| `cubic/3d/each40/eval10000/ascending` | 378 µs | 324 µs | 319 µs | 164 ms |  |  | 0.98 ns | **3.9e+05x slower** |
| `cubic/2d/each10/eval1000/ascending` | 1.61 µs | 1.61 µs | 1.71 µs | 1.11 ms |  |  | 0.88 ns | **1.8e+03x slower** |
| `cubic/2d/each40/eval1000/ascending` | 9.42 µs | 9.37 µs | 9.16 µs | 2.03 ms |  |  | 0.88 ns | **1.1e+04x slower** |
| `cubic/2d/each100/eval1000/ascending` | 44.6 µs | 48.2 µs | 44.2 µs | 8 ms |  |  | 0.98 ns | **4.6e+04x slower** |
| `cubic/2d/each200/eval1000/ascending` | 197 µs | 208 µs | 218 µs | 30.3 ms |  |  | 0.98 ns | **2e+05x slower** |
| `cubic/2d/each400/eval1000/ascending` | 1.53 ms | 1.5 ms | 1.21 ms | 155 ms |  |  | 0.98 ns | **1.6e+06x slower** |
| `quintic/2d/each10/eval1000/ascending` | 2.11 µs | 2.12 µs | 2.21 µs | 3.79 ms |  |  |  | 1.8e+03x faster |
| `quintic/2d/each40/eval1000/ascending` | 13.4 µs | 13.3 µs | 13.4 µs | 9.39 ms |  |  |  | 7e+02x faster |
| `linear/2d/each400/eval10000/ascending` | 333 µs | 348 µs | 116 µs | 6.57 µs | 298 µs | 1.26 µs | 0.78 ns | **4.3e+05x slower** |
| `linear/2d/each800/eval10000/ascending` | 1.21 ms | 1.34 ms | 443 µs | 6.93 µs | 838 µs | 2.45 µs | 0.78 ns | **1.5e+06x slower** |
| `nearest/2d/auckland-ba32-0401/heldout/ascending` | 16.8 ms | 15.5 ms | 5.34 ms | 7.05 µs | 13.7 ms |  | 0.78 ns | **2.2e+07x slower** |
| `linear/2d/auckland-ba32-0401/heldout/ascending` | 17.2 ms | 17.6 ms | 6.04 ms | 7.58 µs | 13.8 ms | 8.34 µs | 0.78 ns | **2.2e+07x slower** |
| `cubic/2d/auckland-ba32-0401/heldout/ascending` | 78.8 ms | 79.1 ms | 74.3 ms | 17.8 s |  |  | 0.88 ns | **9e+07x slower** |
| `quintic/2d/auckland-ba32-0401/heldout/ascending` | 109 ms | 97.8 ms | 90.7 ms | 23.7 s |  |  |  | 2.2e+02x faster |
| `pchip/2d/auckland-ba32-0401/heldout/ascending` | 78.8 ms | 76.6 ms | 66.2 ms |  |  |  |  |  |

## `NearestNDInterpolator`

`n` points of the unit square, evaluated at the pixels of a square image.

### Evaluation

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: |
| `n10/image50` | 36.7 µs | 36.7 µs | 123 µs | 498 µs | 14x faster |
| `n10/image100` | 140 µs | 59.7 µs | 491 µs | 1.87 ms | 13x faster |
| `n10/image500` | 3.31 ms | 614 µs | 12.7 ms | 49 ms | 15x faster |
| `n100/image50` | 126 µs | 91.7 µs | 207 µs | 644 µs | 5.1x faster |
| `n100/image100` | 468 µs | 141 µs | 793 µs | 2.39 ms | 5.1x faster |
| `n100/image500` | 15.8 ms | 2.32 ms | 19.7 ms | 66 ms | 4.2x faster |
| `n1000/image50` | 358 µs | 123 µs | 355 µs | 881 µs | 2.5x faster |
| `n1000/image100` | 1.28 ms | 325 µs | 1.28 ms | 3.33 ms | 2.6x faster |
| `n1000/image500` | 25.4 ms | 4.86 ms | 26.2 ms | 73.6 ms | 2.9x faster |
| `n10000/image50` | 364 µs | 112 µs | 364 µs | 1.08 ms | 3x faster |
| `n10000/image100` | 1.52 ms | 350 µs | 1.51 ms | 4.24 ms | 2.8x faster |
| `n10000/image500` | 33 ms | 5.99 ms | 33.7 ms | 94.5 ms | 2.9x faster |
| `n100000/image50` | 488 µs | 132 µs | 490 µs | 1.49 ms | 3.1x faster |
| `n100000/image100` | 1.83 ms | 399 µs | 1.97 ms | 5.71 ms | 3.1x faster |
| `n100000/image500` | 48.2 ms | 7.59 ms | 48 ms | 130 ms | 2.7x faster |
| `ghcn-tmax-20250101/heldout` | 172 µs | 72.1 µs | 173 µs | 546 µs | 3.2x faster |

### Construction

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: |
| `n10/image50` | 583 ns | 576 ns | 567 ns | 10.4 µs | 18x faster |
| `n10/image100` | 577 ns | 568 ns | 570 ns | 10.4 µs | 18x faster |
| `n10/image500` | 582 ns | 573 ns | 579 ns | 10.4 µs | 18x faster |
| `n100/image50` | 3.7 µs | 3.62 µs | 3.64 µs | 16.7 µs | 4.5x faster |
| `n100/image100` | 3.66 µs | 3.6 µs | 3.61 µs | 16.4 µs | 4.5x faster |
| `n100/image500` | 4.01 µs | 3.67 µs | 4.04 µs | 16.6 µs | 4.1x faster |
| `n1000/image50` | 38.5 µs | 40 µs | 38.2 µs | 135 µs | 3.5x faster |
| `n1000/image100` | 39.4 µs | 38.1 µs | 38 µs | 130 µs | 3.3x faster |
| `n1000/image500` | 38.1 µs | 38.3 µs | 38.4 µs | 127 µs | 3.3x faster |
| `n10000/image50` | 547 µs | 555 µs | 548 µs | 1.83 ms | 3.3x faster |
| `n10000/image100` | 544 µs | 611 µs | 581 µs | 1.86 ms | 3.4x faster |
| `n10000/image500` | 539 µs | 542 µs | 543 µs | 1.8 ms | 3.3x faster |
| `n100000/image50` | 8.83 ms | 8.78 ms | 8.58 ms | 26.9 ms | 3x faster |
| `n100000/image100` | 8.87 ms | 8.53 ms | 8.98 ms | 26.2 ms | 3x faster |
| `n100000/image500` | 9.12 ms | 8.9 ms | 8.66 ms | 27 ms | 3x faster |
| `ghcn-tmax-20250101/heldout` | 533 µs | 561 µs | 559 µs | 1.84 ms | 3.4x faster |

## `LinearNDInterpolator`

`n` points of the unit square, evaluated at the pixels of a square image (some outside the hull, with no value); `3d`: `n` points of the unit cube, evaluated at 10 000 points of the cube in no order.

### Evaluation

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | spade | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `n10/image50` | 20.1 µs | 20.2 µs | 20.1 µs | 30.2 µs | 106 µs | 1.5x faster |
| `n10/image100` | 81.6 µs | 82 µs | 82.3 µs | 107 µs | 420 µs | 1.3x faster |
| `n10/image500` | 2.08 ms | 437 µs | 1.98 ms | 2.9 ms | 10.6 ms | 1.4x faster |
| `n100/image50` | 58.2 µs | 58.4 µs | 58.1 µs | 74.2 µs | 227 µs | 1.3x faster |
| `n100/image100` | 251 µs | 84.8 µs | 249 µs | 254 µs | 878 µs | 1x faster |
| `n100/image500` | 5.81 ms | 1.16 ms | 5.7 ms | 6.07 ms | 18 ms | 1x faster |
| `n1000/image50` | 62 µs | 63.1 µs | 62 µs | 133 µs | 304 µs | 2.1x faster |
| `n1000/image100` | 321 µs | 97.2 µs | 319 µs | 396 µs | 1.22 ms | 1.2x faster |
| `n1000/image500` | 6.79 ms | 1.35 ms | 6.92 ms | 6.88 ms | 21.8 ms | 1x faster |
| `n10000/image50` | 73.9 µs | 75.5 µs | 75.1 µs | 414 µs | 575 µs | 5.6x faster |
| `n10000/image100` | 382 µs | 110 µs | 383 µs | 1.01 ms | 1.75 ms | 2.6x faster |
| `n10000/image500` | 8.78 ms | 1.58 ms | 8.69 ms | 9.36 ms | 27.2 ms | 1.1x faster |
| `n100000/image50` | 98.6 µs | 105 µs | 102 µs | 1.42 ms | 3.67 ms | 14x faster |
| `n100000/image100` | 533 µs | 152 µs | 500 µs | 3.74 ms | 10.8 ms | 7x faster |
| `n100000/image500` | 12.8 ms | 2.29 ms | 13.1 ms | 22 ms | 52.6 ms | 1.7x faster |
| `3d/n1000/scattered` | 897 µs | 205 µs | 904 µs | 5.29 ms |  | 5.9x faster |
| `3d/n10000/scattered` | 1.5 ms | 294 µs | 1.45 ms | 16.6 ms |  | 11x faster |
| `3d/n100000/scattered` | 4.31 ms | 831 µs | 4.04 ms | 72 ms |  | 17x faster |
| `ghcn-tmax-20250101/heldout` | 99.5 µs | 105 µs | 95.7 µs | 1.19 ms | 1.23 ms | 12x faster |

### Construction

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | spade | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `n10/image50` | 5.81 µs | 5.72 µs | 5.74 µs | 490 µs | 1.91 µs | **3x slower** |
| `n10/image100` | 5.86 µs | 5.73 µs | 5.71 µs | 569 µs | 1.91 µs | **3.1x slower** |
| `n10/image500` | 5.8 µs | 5.01 µs | 5.62 µs | 476 µs | 1.9 µs | **3.1x slower** |
| `n100/image50` | 35.8 µs | 36.2 µs | 35.5 µs | 767 µs | 20.8 µs | **1.7x slower** |
| `n100/image100` | 35.8 µs | 34.6 µs | 35.3 µs | 647 µs | 20.6 µs | **1.7x slower** |
| `n100/image500` | 35.9 µs | 34.5 µs | 35.6 µs | 738 µs | 20.7 µs | **1.7x slower** |
| `n1000/image50` | 614 µs | 615 µs | 615 µs | 2.5 ms | 302 µs | **2x slower** |
| `n1000/image100` | 618 µs | 619 µs | 604 µs | 2.48 ms | 302 µs | **2x slower** |
| `n1000/image500` | 605 µs | 625 µs | 615 µs | 2.81 ms | 308 µs | **2x slower** |
| `n10000/image50` | 7.59 ms | 7.55 ms | 7.62 ms | 25.5 ms | 4.13 ms | **1.8x slower** |
| `n10000/image100` | 7.72 ms | 7.49 ms | 7.8 ms | 26.2 ms | 4.08 ms | **1.9x slower** |
| `n10000/image500` | 8.31 ms | 7.39 ms | 6.99 ms | 25.5 ms | 3.72 ms | **2.2x slower** |
| `n100000/image50` | 76.6 ms | 77.2 ms | 74.6 ms | 385 ms | 48.7 ms | **1.6x slower** |
| `n100000/image100` | 82.6 ms | 76.3 ms | 76.5 ms | 378 ms | 48.8 ms | **1.7x slower** |
| `n100000/image500` | 77.8 ms | 76.8 ms | 78.5 ms | 370 ms | 46.7 ms | **1.7x slower** |
| `3d/n1000/scattered` | 4.63 ms | 4.53 ms | 4.56 ms | 9.47 ms |  | 2x faster |
| `3d/n10000/scattered` | 52.1 ms | 50.7 ms | 50.8 ms | 129 ms |  | 2.5x faster |
| `3d/n100000/scattered` | 579 ms | 557 ms | 567 ms | 1.97 s |  | 3.4x faster |
| `ghcn-tmax-20250101/heldout` | 7.76 ms | 7.55 ms | 7.47 ms | 31.4 ms | 4.6 ms | **1.7x slower** |

## `RBFInterpolator`

`n` points of a square, evaluated at 1000 points in no order (`scattered`) or at the 10 000 pixels of an image (`image`). `global`: all points; `k50`: the 50 nearest.

### Evaluation

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/global/n10/scattered` | 15.8 µs | 16 µs | 15.8 µs | 66.6 µs | 200 µs | 4.2x faster |
| `linear/k50/n10/scattered` | 15.9 µs | 15.9 µs | 15.8 µs | 1.55 ms |  | 98x faster |
| `thin_plate_spline/global/n10/scattered` | 78.2 µs | 78.6 µs | 77.8 µs | 122 µs | 326 µs | 1.6x faster |
| `thin_plate_spline/k50/n10/scattered` | 78.1 µs | 78.2 µs | 78.8 µs | 1.64 ms |  | 21x faster |
| `cubic/global/n10/scattered` | 17.4 µs | 16.6 µs | 17.9 µs | 63.8 µs | 217 µs | 3.7x faster |
| `cubic/k50/n10/scattered` | 17.5 µs | 20.5 µs | 17.6 µs | 1.49 ms |  | 85x faster |
| `quintic/global/n10/scattered` | 25 µs | 25.2 µs | 25.7 µs | 76 µs | 238 µs | 3x faster |
| `quintic/k50/n10/scattered` | 25.4 µs | 24.8 µs | 25.9 µs | 1.56 ms |  | 61x faster |
| `multiquadric/global/n10/scattered` | 30.9 µs | 30.7 µs | 31 µs | 58.8 µs | 264 µs | 1.9x faster |
| `multiquadric/k50/n10/scattered` | 30.5 µs | 30.7 µs | 30.5 µs | 1.55 ms |  | 51x faster |
| `inverse_multiquadric/global/n10/scattered` | 41.1 µs | 41.2 µs | 41.1 µs | 60.5 µs | 294 µs | 1.5x faster |
| `inverse_multiquadric/k50/n10/scattered` | 41.1 µs | 41.2 µs | 41.1 µs | 1.56 ms |  | 38x faster |
| `inverse_quadratic/global/n10/scattered` | 26.2 µs | 26.4 µs | 26.5 µs | 53.6 µs |  | 2x faster |
| `inverse_quadratic/k50/n10/scattered` | 26.2 µs | 26.3 µs | 26.2 µs | 1.56 ms |  | 60x faster |
| `gaussian/global/n10/scattered` | 73.8 µs | 73.9 µs | 73.6 µs | 97.6 µs | 390 µs | 1.3x faster |
| `gaussian/k50/n10/scattered` | 74 µs | 73.9 µs | 73.7 µs | 1.52 ms |  | 21x faster |
| `linear/global/n100/scattered` | 148 µs | 132 µs | 149 µs | 320 µs | 1.89 ms | 2.2x faster |
| `linear/k50/n100/scattered` | 15.3 ms | 7.85 ms | 16 ms | 20.1 ms |  | 1.3x faster |
| `thin_plate_spline/global/n100/scattered` | 765 µs | 202 µs | 771 µs | 956 µs | 3.15 ms | 1.3x faster |
| `thin_plate_spline/k50/n100/scattered` | 25.7 ms | 4.65 ms | 23.4 ms | 31.7 ms |  | 1.2x faster |
| `cubic/global/n100/scattered` | 139 µs | 71.8 µs | 151 µs | 322 µs | 2.05 ms | 2.3x faster |
| `cubic/k50/n100/scattered` | 16.2 ms | 3.2 ms | 16.2 ms | 22.3 ms |  | 1.4x faster |
| `quintic/global/n100/scattered` | 182 µs | 67 µs | 181 µs | 343 µs | 2.25 ms | 1.9x faster |
| `quintic/k50/n100/scattered` | 17.3 ms | 3.33 ms | 17.5 ms | 24.2 ms |  | 1.4x faster |
| `multiquadric/global/n100/scattered` | 295 µs | 94.8 µs | 295 µs | 330 µs | 2.47 ms | 1.1x faster |
| `multiquadric/k50/n100/scattered` | 16.3 ms | 3.37 ms | 16.4 ms | 20.3 ms |  | 1.3x faster |
| `inverse_multiquadric/global/n100/scattered` | 399 µs | 112 µs | 399 µs | 434 µs | 2.82 ms | 1.1x faster |
| `inverse_multiquadric/k50/n100/scattered` | 16.8 ms | 3.16 ms | 17.1 ms | 21 ms |  | 1.2x faster |
| `inverse_quadratic/global/n100/scattered` | 271 µs | 77 µs | 249 µs | 295 µs |  | 1.1x faster |
| `inverse_quadratic/k50/n100/scattered` | 16.9 ms | 3.15 ms | 14.7 ms | 24 ms |  | 1.4x faster |
| `gaussian/global/n100/scattered` | 676 µs | 180 µs | 689 µs | 755 µs | 3.85 ms | 1.1x faster |
| `gaussian/k50/n100/scattered` | 21.6 ms | 4.03 ms | 22.3 ms | 22.9 ms |  | 1.1x faster |
| `linear/global/n1000/scattered` | 1.48 ms | 319 µs | 1.48 ms | 4.06 ms | 19.2 ms | 2.8x faster |
| `linear/k50/n1000/scattered` | 14.7 ms | 2.78 ms | 14.9 ms | 37.1 ms |  | 2.5x faster |
| `thin_plate_spline/global/n1000/scattered` | 9.54 ms | 2.6 ms | 7.89 ms | 14.7 ms | 36.1 ms | 1.5x faster |
| `thin_plate_spline/k50/n1000/scattered` | 24.6 ms | 6.13 ms | 23.6 ms | 50.3 ms |  | 2x faster |
| `cubic/global/n1000/scattered` | 1.7 ms | 385 µs | 1.48 ms | 4.17 ms | 21.8 ms | 2.4x faster |
| `cubic/k50/n1000/scattered` | 17.2 ms | 3.16 ms | 16.8 ms | 39.7 ms |  | 2.3x faster |
| `quintic/global/n1000/scattered` | 1.65 ms | 397 µs | 1.65 ms | 4.12 ms | 23 ms | 2.5x faster |
| `quintic/k50/n1000/scattered` | 18.2 ms | 3.54 ms | 18.4 ms | 43.9 ms |  | 2.4x faster |
| `multiquadric/global/n1000/scattered` | 2.96 ms | 689 µs | 2.97 ms | 4.39 ms | 24.8 ms | 1.5x faster |
| `multiquadric/k50/n1000/scattered` | 17.1 ms | 3.33 ms | 16.8 ms | 44.7 ms |  | 2.6x faster |
| `inverse_multiquadric/global/n1000/scattered` | 4.1 ms | 928 µs | 4.03 ms | 5.53 ms | 30.8 ms | 1.4x faster |
| `inverse_multiquadric/k50/n1000/scattered` | 18.2 ms | 3.74 ms | 18.4 ms | 39.8 ms |  | 2.2x faster |
| `inverse_quadratic/global/n1000/scattered` | 2.45 ms | 561 µs | 2.46 ms | 3.95 ms |  | 1.6x faster |
| `inverse_quadratic/k50/n1000/scattered` | 16.3 ms | 3.1 ms | 16.3 ms | 36.5 ms |  | 2.2x faster |
| `gaussian/global/n1000/scattered` | 6.76 ms | 1.33 ms | 6.73 ms | 8.83 ms | 38.4 ms | 1.3x faster |
| `gaussian/k50/n1000/scattered` | 24 ms | 4.17 ms | 23.4 ms | 45.3 ms |  | 1.9x faster |
| `thin_plate_spline/global/n3000/scattered` | 23.7 ms | 4.34 ms | 23.6 ms | 33.6 ms | 106 ms | 1.4x faster |
| `thin_plate_spline/k50/n10000/scattered` | 24.5 ms | 4.42 ms | 24.3 ms | 48.9 ms |  | 2x faster |
| `thin_plate_spline/k50/n10000/image` | 250 ms | 40.1 ms | 250 ms | 508 ms |  | 2x faster |
| `thin_plate_spline/k50/n100000/scattered` | 32.2 ms | 5 ms | 25.8 ms | 51.4 ms |  | 1.6x faster |
| `thin_plate_spline/k50/n100000/image` | 257 ms | 39.1 ms | 251 ms | 510 ms |  | 2x faster |
| `linear/global/ghcn-tmax-20250101/heldout` | 20.2 ms | 3.78 ms | 20.3 ms | 55.5 ms | 252 ms | 2.7x faster |
| `linear/k50/ghcn-tmax-20250101/heldout` | 17.9 ms | 3.49 ms | 18.2 ms | 45.4 ms |  | 2.5x faster |
| `thin_plate_spline/global/ghcn-tmax-20250101/heldout` | 87.4 ms | 17.1 ms | 94.1 ms | 144 ms | 462 ms | 1.6x faster |
| `thin_plate_spline/k50/ghcn-tmax-20250101/heldout` | 29.3 ms | 5.43 ms | 31.7 ms | 57.2 ms |  | 2x faster |
| `cubic/global/ghcn-tmax-20250101/heldout` | 25.2 ms | 4.51 ms | 25.6 ms | 57.5 ms | 295 ms | 2.3x faster |
| `cubic/k50/ghcn-tmax-20250101/heldout` | 22.4 ms | 3.94 ms | 21.7 ms | 51.3 ms |  | 2.3x faster |
| `quintic/global/ghcn-tmax-20250101/heldout` | 22 ms | 4.27 ms | 22.3 ms | 63.5 ms | 294 ms | 2.9x faster |
| `quintic/k50/ghcn-tmax-20250101/heldout` | 22.7 ms | 4.13 ms | 21.8 ms | 54 ms |  | 2.4x faster |

### Construction

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/global/n10/scattered` | 2.86 µs | 2.85 µs | 2.89 µs | 12.8 µs | 5.95 µs | 2.1x faster |
| `linear/k50/n10/scattered` | 2.89 µs | 2.85 µs | 2.87 µs | 16.3 µs |  | 5.6x faster |
| `thin_plate_spline/global/n10/scattered` | 3.84 µs | 3.66 µs | 3.35 µs | 14.2 µs | 7.14 µs | 1.9x faster |
| `thin_plate_spline/k50/n10/scattered` | 3.64 µs | 3.62 µs | 3.63 µs | 17.3 µs |  | 4.7x faster |
| `cubic/global/n10/scattered` | 3.83 µs | 3.82 µs | 3.81 µs | 14.1 µs | 6.1 µs | 1.6x faster |
| `cubic/k50/n10/scattered` | 3.86 µs | 3.88 µs | 3.8 µs | 17.1 µs |  | 4.4x faster |
| `quintic/global/n10/scattered` | 4.38 µs | 4.31 µs | 4.25 µs | 16.8 µs | 6.41 µs | 1.5x faster |
| `quintic/k50/n10/scattered` | 4.4 µs | 4.28 µs | 4.4 µs | 20.8 µs |  | 4.7x faster |
| `multiquadric/global/n10/scattered` | 3.02 µs | 2.95 µs | 3 µs | 13.5 µs | 7.38 µs | 2.4x faster |
| `multiquadric/k50/n10/scattered` | 2.96 µs | 2.94 µs | 2.91 µs | 17.5 µs |  | 5.9x faster |
| `inverse_multiquadric/global/n10/scattered` | 2.94 µs | 2.92 µs | 2.95 µs | 12.8 µs | 7.56 µs | 2.6x faster |
| `inverse_multiquadric/k50/n10/scattered` | 2.92 µs | 2.92 µs | 2.94 µs | 16.6 µs |  | 5.7x faster |
| `inverse_quadratic/global/n10/scattered` | 2.88 µs | 2.86 µs | 2.9 µs | 13.4 µs |  | 4.7x faster |
| `inverse_quadratic/k50/n10/scattered` | 2.85 µs | 2.85 µs | 2.87 µs | 16.2 µs |  | 5.7x faster |
| `gaussian/global/n10/scattered` | 3.04 µs | 3.07 µs | 3.05 µs | 12.9 µs | 8.55 µs | 2.8x faster |
| `gaussian/k50/n10/scattered` | 3.06 µs | 3.04 µs | 3.11 µs | 16.4 µs |  | 5.4x faster |
| `linear/global/n100/scattered` | 73.8 µs | 118 µs | 77.4 µs | 84.5 µs | 758 µs | 1.1x faster |
| `linear/k50/n100/scattered` | 3.55 µs | 5.9 µs | 3.46 µs | 23 µs |  | 6.5x faster |
| `thin_plate_spline/global/n100/scattered` | 107 µs | 106 µs | 107 µs | 124 µs | 899 µs | 1.2x faster |
| `thin_plate_spline/k50/n100/scattered` | 3.66 µs | 3.49 µs | 3.6 µs | 23.6 µs |  | 6.4x faster |
| `cubic/global/n100/scattered` | 74.8 µs | 78.2 µs | 77.8 µs | 88.8 µs | 776 µs | 1.2x faster |
| `cubic/k50/n100/scattered` | 3.6 µs | 3.6 µs | 3.64 µs | 26.6 µs |  | 7.4x faster |
| `quintic/global/n100/scattered` | 85.8 µs | 84.2 µs | 84.2 µs | 106 µs | 796 µs | 1.2x faster |
| `quintic/k50/n100/scattered` | 3.55 µs | 3.56 µs | 3.61 µs | 25.1 µs |  | 7.1x faster |
| `multiquadric/global/n100/scattered` | 82.8 µs | 82.7 µs | 86 µs | 83.3 µs | 674 µs | 1x faster |
| `multiquadric/k50/n100/scattered` | 3.45 µs | 3.47 µs | 3.54 µs | 22.7 µs |  | 6.6x faster |
| `inverse_multiquadric/global/n100/scattered` | 87.3 µs | 87 µs | 89.5 µs | 89.8 µs | 821 µs | 1x faster |
| `inverse_multiquadric/k50/n100/scattered` | 3.48 µs | 3.54 µs | 3.47 µs | 24.5 µs |  | 7x faster |
| `inverse_quadratic/global/n100/scattered` | 112 µs | 79 µs | 124 µs | 81.1 µs |  | **1.4x slower** |
| `inverse_quadratic/k50/n100/scattered` | 3.75 µs | 3.51 µs | 3.44 µs | 38.2 µs |  | 10x faster |
| `gaussian/global/n100/scattered` | 104 µs | 106 µs | 159 µs | 107 µs | 824 µs | 1x faster |
| `gaussian/k50/n100/scattered` | 3.49 µs | 3.47 µs | 3.47 µs | 22.9 µs |  | 6.5x faster |
| `linear/global/n1000/scattered` | 20.4 ms | 19.6 ms | 20.2 ms | 58.4 ms | 212 ms | 2.9x faster |
| `linear/k50/n1000/scattered` | 39.3 µs | 38.7 µs | 38.6 µs | 145 µs |  | 3.7x faster |
| `thin_plate_spline/global/n1000/scattered` | 26.1 ms | 44.6 ms | 24.4 ms | 63 ms | 240 ms | 2.4x faster |
| `thin_plate_spline/k50/n1000/scattered` | 39.2 µs | 39 µs | 38.7 µs | 148 µs |  | 3.8x faster |
| `cubic/global/n1000/scattered` | 42.3 ms | 21.3 ms | 35.5 ms | 61.9 ms | 246 ms | 1.5x faster |
| `cubic/k50/n1000/scattered` | 39 µs | 37.1 µs | 38.6 µs | 141 µs |  | 3.6x faster |
| `quintic/global/n1000/scattered` | 20.8 ms | 20.9 ms | 21.3 ms | 55.4 ms | 217 ms | 2.7x faster |
| `quintic/k50/n1000/scattered` | 38.9 µs | 38.8 µs | 38.3 µs | 144 µs |  | 3.7x faster |
| `multiquadric/global/n1000/scattered` | 21.5 ms | 22.9 ms | 21.3 ms | 62.2 ms | 218 ms | 2.9x faster |
| `multiquadric/k50/n1000/scattered` | 38.6 µs | 38.8 µs | 38.2 µs | 144 µs |  | 3.7x faster |
| `inverse_multiquadric/global/n1000/scattered` | 24.5 ms | 22.5 ms | 22.2 ms | 63.6 ms | 220 ms | 2.6x faster |
| `inverse_multiquadric/k50/n1000/scattered` | 38.5 µs | 38.6 µs | 38.8 µs | 137 µs |  | 3.5x faster |
| `inverse_quadratic/global/n1000/scattered` | 21.9 ms | 21.9 ms | 21.1 ms | 54 ms |  | 2.5x faster |
| `inverse_quadratic/k50/n1000/scattered` | 38.6 µs | 38.7 µs | 38.1 µs | 136 µs |  | 3.5x faster |
| `gaussian/global/n1000/scattered` | 24.4 ms | 23.3 ms | 23.5 ms | 67.4 ms | 225 ms | 2.8x faster |
| `gaussian/k50/n1000/scattered` | 39.1 µs | 38.6 µs | 38.3 µs | 141 µs |  | 3.6x faster |
| `thin_plate_spline/global/n3000/scattered` | 495 ms | 472 ms | 487 ms | 2.13 s | 7.44 s | 4.3x faster |
| `thin_plate_spline/k50/n10000/scattered` | 782 µs | 755 µs | 672 µs | 1.78 ms |  | 2.3x faster |
| `thin_plate_spline/k50/n10000/image` | 828 µs | 809 µs | 707 µs | 1.76 ms |  | 2.1x faster |
| `thin_plate_spline/k50/n100000/scattered` | 9.67 ms | 9.3 ms | 9.18 ms | 29.3 ms |  | 3x faster |
| `thin_plate_spline/k50/n100000/image` | 8.89 ms | 9.02 ms | 8.73 ms | 27.7 ms |  | 3.1x faster |
| `linear/global/ghcn-tmax-20250101/heldout` | 19.3 s | 18.8 s | 19.1 s | 128 s | 354 s | 6.6x faster |
| `linear/k50/ghcn-tmax-20250101/heldout` | 749 µs | 829 µs | 704 µs | 1.86 ms |  | 2.5x faster |
| `thin_plate_spline/global/ghcn-tmax-20250101/heldout` | 20.5 s | 20.2 s | 21.2 s | 141 s | 354 s | 6.9x faster |
| `thin_plate_spline/k50/ghcn-tmax-20250101/heldout` | 724 µs | 750 µs | 701 µs | 1.84 ms |  | 2.5x faster |
| `cubic/global/ghcn-tmax-20250101/heldout` | 20.5 s | 23 s | 23.8 s | 131 s | 367 s | 6.4x faster |
| `cubic/k50/ghcn-tmax-20250101/heldout` | 715 µs | 837 µs | 724 µs | 1.91 ms |  | 2.7x faster |
| `quintic/global/ghcn-tmax-20250101/heldout` | 18.9 s | 21.4 s | 18.1 s | 161 s | 310 s | 8.5x faster |
| `quintic/k50/ghcn-tmax-20250101/heldout` | 833 µs | 751 µs | 706 µs | 1.74 ms |  | 2.1x faster |

## Memory

Megabytes. "Data read" is the process with the data of the case in memory, before anything is built. The columns after it are what was added to that: the highest while building, with the interpolator built, the highest while evaluating. "Highest" is the highest of the whole run, repeats included.

| case | library | data read | while building | built | while evaluating | highest |
| :--- | :--- | ---: | ---: | ---: | ---: | ---: |
| `oned/cubic_spline/n10000/sorted` | interpx | 4.7 | +0.5 | +0.6 | +0.6 | 5.9 |
| `oned/cubic_spline/n10000/sorted` | interpx, parallel | 3.3 | +2.0 | +2.0 | +2.0 | 6.5 |
| `oned/cubic_spline/n10000/sorted` | SciPy | 71.2 | +1.1 | +1.9 | +1.9 | 73.4 |
| `oned/cubic_spline/n10000/sorted` | scirs2-interpolate | 4.0 | +1.1 | +1.2 | +1.9 | 5.9 |
| `oned/cubic_spline/n10000/sorted` | ndarray-interp | 3.4 | +1.7 | +1.5 | +1.5 | 5.6 |
| `grid/linear/4d/max200/eval10000/ascending` | interpx | 233.3 | +171.8 | +171.8 | +171.8 | 405.1 |
| `grid/linear/4d/max200/eval10000/ascending` | interpx, parallel | 205.3 | +199.8 | +199.8 | +199.8 | 405.4 |
| `grid/linear/4d/max200/eval10000/ascending` | interpx, `from_owned` | 223.4 | +181.6 | +181.3 | +181.3 | 405.0 |
| `grid/linear/4d/max200/eval10000/ascending` | SciPy | 271.2 | +0.3 | +0.3 | +2.0 | 273.5 |
| `grid/linear/4d/max200/eval10000/ascending` | scirs2-interpolate | 209.8 | +187.3 | +195.1 | +195.1 | 404.9 |
| `grid/linear/4d/max200/eval10000/ascending` | interpn | 204.6 | +0.3 | +0.3 | +0.3 | 399.7 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpx | 3.5 | +4.7 | +4.8 | +4.8 | 8.5 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpx, parallel | 5.6 | +2.8 | +2.7 | +2.7 | 8.9 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpx, `from_owned` | 5.6 | +2.8 | +2.7 | +2.7 | 8.5 |
| `grid/cubic/2d/each400/eval1000/ascending` | SciPy | 72.4 | +85.4 | +6.4 | +6.4 | 162.5 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpn | 4.5 | +1.0 | +1.0 | +1.0 | 5.6 |
| `grid/linear/2d/each800/eval10000/ascending` | interpx | 10.0 | +0.0 | +4.8 | +4.8 | 14.9 |
| `grid/linear/2d/each800/eval10000/ascending` | interpx, parallel | 13.1 | +1.9 | +1.9 | +1.9 | 15.4 |
| `grid/linear/2d/each800/eval10000/ascending` | interpx, `from_owned` | 9.7 | +3.3 | +5.2 | +5.2 | 15.0 |
| `grid/linear/2d/each800/eval10000/ascending` | SciPy | 76.3 | +0.2 | +0.2 | +0.8 | 77.1 |
| `grid/linear/2d/each800/eval10000/ascending` | scirs2-interpolate | 10.6 | +0.0 | +4.3 | +4.4 | 15.0 |
| `grid/linear/2d/each800/eval10000/ascending` | ndarray-interp | 13.6 | +0.0 | +0.0 | +0.0 | 13.6 |
| `grid/linear/2d/each800/eval10000/ascending` | interpn | 12.1 | +0.0 | +0.0 | +0.0 | 12.1 |
| `nearest_nd/n100000/image500` | interpx | 13.9 | +1.5 | +0.5 | +2.5 | 17.8 |
| `nearest_nd/n100000/image500` | interpx, parallel | 14.1 | +1.4 | +0.3 | +2.7 | 18.9 |
| `nearest_nd/n100000/image500` | interpx, `from_owned` | 13.3 | +2.6 | +0.9 | +2.9 | 18.2 |
| `nearest_nd/n100000/image500` | SciPy | 77.5 | +4.6 | +4.7 | +12.9 | 92.1 |
| `linear_nd/n100000/image500` | interpx | 13.5 | +24.3 | +21.2 | +23.1 | 41.5 |
| `linear_nd/n100000/image500` | interpx, parallel | 14.8 | +21.9 | +19.9 | +22.4 | 42.4 |
| `linear_nd/n100000/image500` | interpx, `from_owned` | 13.2 | +24.9 | +21.8 | +23.8 | 42.3 |
| `linear_nd/n100000/image500` | SciPy | 77.4 | +67.9 | +14.0 | +25.9 | 145.6 |
| `linear_nd/n100000/image500` | spade | 13.5 | +14.6 | +11.6 | +13.6 | 29.5 |
| `rbf/thin_plate_spline/global/n3000/scattered` | interpx | 3.4 | +77.1 | +4.8 | +4.8 | 80.4 |
| `rbf/thin_plate_spline/global/n3000/scattered` | interpx, parallel | 3.4 | +77.1 | +4.8 | +5.2 | 80.9 |
| `rbf/thin_plate_spline/global/n3000/scattered` | interpx, `from_owned` | 3.3 | +77.2 | +4.8 | +4.9 | 80.6 |
| `rbf/thin_plate_spline/global/n3000/scattered` | SciPy | 71.2 | +73.0 | +0.9 | +9.0 | 144.2 |
| `rbf/thin_plate_spline/global/n3000/scattered` | scirs2-interpolate | 4.3 | +288.3 | +0.2 | +0.2 | 292.5 |
| `rbf/thin_plate_spline/k50/n100000/image` | interpx | 6.6 | +3.6 | +4.8 | +5.1 | 11.9 |
| `rbf/thin_plate_spline/k50/n100000/image` | interpx, parallel | 6.1 | +4.4 | +5.4 | +6.5 | 12.7 |
| `rbf/thin_plate_spline/k50/n100000/image` | interpx, `from_owned` | 6.7 | +6.1 | +4.7 | +5.1 | 13.5 |
| `rbf/thin_plate_spline/k50/n100000/image` | SciPy | 73.7 | +5.2 | +5.4 | +21.8 | 101.1 |

### Over time

The resident set size of each library's process over its run, with the ends of the phases. The time after `evaluated` is the timed repeats.

`oned/cubic_spline/n10000/sorted`

![memory over time, oned/cubic_spline/n10000/sorted](performance/oned-cubic_spline-n10000-sorted.svg)

`grid/linear/4d/max200/eval10000/ascending`

![memory over time, grid/linear/4d/max200/eval10000/ascending](performance/grid-linear-4d-max200-eval10000-ascending.svg)

`grid/cubic/2d/each400/eval1000/ascending`

![memory over time, grid/cubic/2d/each400/eval1000/ascending](performance/grid-cubic-2d-each400-eval1000-ascending.svg)

`grid/linear/2d/each800/eval10000/ascending`

![memory over time, grid/linear/2d/each800/eval10000/ascending](performance/grid-linear-2d-each800-eval10000-ascending.svg)

`nearest_nd/n100000/image500`

![memory over time, nearest_nd/n100000/image500](performance/nearest_nd-n100000-image500.svg)

`linear_nd/n100000/image500`

![memory over time, linear_nd/n100000/image500](performance/linear_nd-n100000-image500.svg)

`rbf/thin_plate_spline/global/n3000/scattered`

![memory over time, rbf/thin_plate_spline/global/n3000/scattered](performance/rbf-thin_plate_spline-global-n3000-scattered.svg)

`rbf/thin_plate_spline/k50/n100000/image`

![memory over time, rbf/thin_plate_spline/k50/n100000/image](performance/rbf-thin_plate_spline-k50-n100000-image.svg)

## Do the libraries compute the same thing

The largest distance of a library's values from SciPy's, over the cases of a family, as a share of the largest value. Under 1e-06 the two compute the same interpolant.

| family | interpx | interpx, parallel | interpx, `from_owned` | scirs2-interpolate | ndarray-interp | interpn | spade |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| One variable | 3e-16 | 3e-16 |  | 6e-02 | 1e-04 |  |  |
| `PPoly` | 9e-16 | 9e-16 |  |  |  |  |  |
| `BSpline` | 4e-16 | 4e-16 |  |  |  |  |  |
| `BSplineCurve` | 4e-16 | 4e-16 |  |  |  |  |  |
| Regular grid | 7e-04 | 7e-04 | 7e-04 | 4e-16 | 6e-16 | 7e-01 |  |
| `NearestNDInterpolator` | 0e+00 | 0e+00 | 0e+00 |  |  |  |  |
| `LinearNDInterpolator` | 6e-15 | 6e-15 | 6e-15 |  |  |  | 6e-15 |
| `RBFInterpolator` | 4e+00 | 4e+00 | 4e+00 | 7e+11 |  |  |  |
