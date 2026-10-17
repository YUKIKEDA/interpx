# Performance: interpx against SciPy and other Rust crates

<!-- Written by compare/perf/report.py from a run of compare/perf/run.py. Do not edit by hand: run it again. -->

Measured on 2026-10-17, at commit `c519bf194b`.

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
| One variable | evaluation | 112 | 111 | 0 | 1 |
| One variable | construction | 112 | 73 | 0 | 39 |
| `PPoly` | evaluation | 4 | 4 | 0 | 0 |
| `PPoly` | construction | 4 | 4 | 0 | 0 |
| `BSpline` | evaluation | 4 | 4 | 0 | 0 |
| `BSpline` | construction | 4 | 4 | 0 | 0 |
| `BSplineCurve` | evaluation | 2 | 2 | 0 | 0 |
| `BSplineCurve` | construction | 2 | 2 | 0 | 0 |
| Regular grid | evaluation | 92 | 74 | 0 | 18 |
| Regular grid | construction | 92 | 29 | 38 | 46 |
| `NearestNDInterpolator` | evaluation | 16 | 16 | 0 | 0 |
| `NearestNDInterpolator` | construction | 16 | 16 | 0 | 0 |
| `LinearNDInterpolator` | evaluation | 19 | 19 | 0 | 0 |
| `LinearNDInterpolator` | construction | 19 | 3 | 0 | 16 |
| `RBFInterpolator` | evaluation | 61 | 61 | 0 | 0 |
| `RBFInterpolator` | construction | 61 | 61 | 0 | 0 |

216 comparisons (a case, an operation, and a library), by method. "Same values" says whether that library's values are SciPy's up to rounding: where they are not, it computes another interpolant under the same name, and its time is not that of the same work.

| method | operation | against | same values | cases | slower by | the largest |
| :--- | :--- | :--- | :--- | ---: | ---: | :--- |
| `grid/cubic` | construction | interpn | no | 8 | 1.89e+03x to 1.53e+08x | `grid/cubic/2d/auckland-ba32-0401/heldout/ascending`: 135 ms against 0.88 ns |
| `grid/linear` | construction | interpn | yes | 35 | 154x to 3.85e+07x | `grid/linear/4d/max200/eval1000/ascending`: 48.8 ms against 1.27 ns |
| `grid/nearest` | construction | interpn | yes | 3 | 775x to 2.12e+07x | `grid/nearest/2d/auckland-ba32-0401/heldout/ascending`: 16.6 ms against 0.78 ns |
| `grid/linear` | construction | SciPy | yes | 35 | 1.17x to 5.08e+03x | `grid/linear/4d/max200/eval10/ascending`: 51.1 ms against 10.1 µs |
| `grid/nearest` | construction | SciPy | yes | 2 | 2.77x to 2.14e+03x | `grid/nearest/2d/auckland-ba32-0401/heldout/ascending`: 16.6 ms against 7.75 µs |
| `grid/linear` | construction | ndarray-interp | yes | 15 | 1.58x to 1.9e+03x | `grid/linear/2d/auckland-ba32-0401/heldout/ascending`: 17.2 ms against 9.05 µs |
| `grid/pchip` | construction | SciPy | yes | 1 | 50.2x | `grid/pchip/3d/each40/eval10000/ascending`: 471 µs against 9.38 µs |
| `oned/akima` | construction | scirs2-interpolate | yes | 8 | 1.48x to 3.13x | `oned/akima/n10000/shuffled`: 141 µs against 45.1 µs |
| `linear_nd` | construction | spade | yes | 16 | 1.41x to 3.11x | `linear_nd/n10/image100`: 5.36 µs against 1.72 µs |
| `grid/linear` | construction | scirs2-interpolate | yes | 28 | 1.11x to 3.01x | `grid/linear/3d/max40/eval100/ascending`: 4.54 µs against 1.51 µs |
| `oned/pchip` | construction | scirs2-interpolate | yes | 12 | 1.38x to 2.83x | `oned/pchip/n10000/shuffled-uneven`: 133 µs against 47 µs |
| `grid/linear` | evaluation | interpn | yes | 17 | 1.1x to 2.72x | `grid/linear/2d/max10/eval1000/ascending`: 27.5 µs against 10.1 µs |
| `grid/nearest` | construction | scirs2-interpolate | yes | 3 | 1.27x to 2.24x | `grid/nearest/3d/each40/eval10000/ascending`: 20.7 µs against 9.25 µs |
| `oned/linear` | construction | ndarray-interp | yes | 12 | 1.27x to 2.13x | `oned/linear/n10/shuffled`: 95.9 ns against 45 ns |
| `oned/linear` | construction | scirs2-interpolate | yes | 14 | 1.16x to 1.76x | `oned/linear/n10000/shuffled-uneven`: 22.7 µs against 12.9 µs |
| `oned/cubic_spline` | construction | scirs2-interpolate | no | 5 | 1.13x to 1.23x | `oned/cubic_spline/n10000/sorted-uneven`: 179 µs against 146 µs |
| `oned/nearest` | evaluation | scirs2-interpolate | yes | 1 | 1.17x | `oned/nearest/n10/shuffled`: 323 ns against 277 ns |
| `grid/nearest` | evaluation | interpn | yes | 1 | 1.11x | `grid/nearest/2d/each40/eval10000/ascending`: 164 µs against 148 µs |

<details><summary>Every comparison where interpx is slower</summary>

"Values apart" is the distance of that library's values from SciPy's, as a share of the largest value.

| case | operation | against | interpx | the other | slower by | values apart |
| :--- | :--- | :--- | ---: | ---: | ---: | ---: |
| `grid/cubic/2d/auckland-ba32-0401/heldout/ascending` | construction | interpn | 135 ms | 0.88 ns | 153134535.11x | 5e-03 |
| `grid/linear/4d/max200/eval1000/ascending` | construction | interpn | 48.8 ms | 1.27 ns | 38477666.46x | 3e-16 |
| `grid/linear/4d/max200/eval10/ascending` | construction | interpn | 51.1 ms | 1.37 ns | 37375341.71x | 2e-16 |
| `grid/linear/4d/max200/eval10000/ascending` | construction | interpn | 46.9 ms | 1.27 ns | 36950882.46x | 4e-16 |
| `grid/linear/4d/max200/eval100/ascending` | construction | interpn | 47.5 ms | 1.37 ns | 34730349.71x | 3e-16 |
| `grid/linear/2d/auckland-ba32-0401/heldout/ascending` | construction | interpn | 17.2 ms | 0.78 ns | 21994624.00x | 0e+00 |
| `grid/nearest/2d/auckland-ba32-0401/heldout/ascending` | construction | interpn | 16.6 ms | 0.78 ns | 21238912.00x | 0e+00 |
| `grid/linear/3d/max200/eval10/ascending` | construction | interpn | 1.82 ms | 0.88 ns | 2072689.78x | 3e-16 |
| `grid/linear/3d/max200/eval100/ascending` | construction | interpn | 1.8 ms | 0.88 ns | 2052096.00x | 3e-16 |
| `grid/linear/3d/max200/eval10000/ascending` | construction | interpn | 1.63 ms | 0.88 ns | 1851960.89x | 5e-16 |
| `grid/linear/3d/max200/eval1000/ascending` | construction | interpn | 1.77 ms | 0.98 ns | 1808384.00x | 3e-16 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | interpn | 1.07 ms | 0.68 ns | 1567012.57x | 3e-16 |
| `grid/cubic/2d/each400/eval1000/ascending` | construction | interpn | 1.33 ms | 0.88 ns | 1510001.78x | 5e-01 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | interpn | 309 µs | 0.68 ns | 452022.86x | 3e-16 |
| `grid/cubic/3d/each40/eval10000/ascending` | construction | interpn | 368 µs | 1.17 ns | 314112.00x | 5e-01 |
| `grid/cubic/2d/each200/eval1000/ascending` | construction | interpn | 174 µs | 0.88 ns | 198257.78x | 2e-01 |
| `grid/cubic/2d/each100/eval1000/ascending` | construction | interpn | 41 µs | 0.88 ns | 46663.11x | 2e-01 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | interpn | 20.7 µs | 0.88 ns | 23562.67x | 0e+00 |
| `grid/linear/4d/max40/eval1000/ascending` | construction | interpn | 14.2 µs | 1.27 ns | 11219.69x | 6e-16 |
| `grid/linear/4d/max40/eval10000/ascending` | construction | interpn | 14.5 µs | 1.37 ns | 10612.57x | 5e-16 |
| `grid/linear/4d/max40/eval10/ascending` | construction | interpn | 14.4 µs | 1.37 ns | 10514.29x | 2e-16 |
| `grid/linear/4d/max40/eval100/ascending` | construction | interpn | 14.2 µs | 1.37 ns | 10409.14x | 2e-16 |
| `grid/linear/2d/max200/eval10/ascending` | construction | interpn | 7.57 µs | 0.78 ns | 9684.00x | 2e-16 |
| `grid/linear/2d/max200/eval10000/ascending` | construction | interpn | 7.54 µs | 0.78 ns | 9650.00x | 3e-16 |
| `grid/cubic/2d/each40/eval10000/ascending` | construction | interpn | 8.39 µs | 0.88 ns | 9550.22x | 7e-01 |
| `grid/cubic/2d/each40/eval1000/ascending` | construction | interpn | 8.39 µs | 0.88 ns | 9550.22x | 7e-01 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | interpn | 7.42 µs | 0.78 ns | 9492.00x | 2e-16 |
| `grid/linear/2d/max200/eval100/ascending` | construction | interpn | 7.41 µs | 0.78 ns | 9484.00x | 1e-16 |
| `grid/linear/4d/max200/eval10/ascending` | construction | SciPy | 51.1 ms | 10.1 µs | 5075.02x |  |
| `grid/linear/4d/max200/eval100/ascending` | construction | SciPy | 47.5 ms | 9.38 µs | 5064.84x |  |
| `grid/linear/4d/max200/eval10000/ascending` | construction | SciPy | 46.9 ms | 9.54 µs | 4916.90x |  |
| `grid/linear/4d/max200/eval100/flipped` | construction | SciPy | 108 ms | 22 µs | 4910.84x |  |
| `grid/linear/4d/max200/eval10000/flipped` | construction | SciPy | 106 ms | 21.7 µs | 4870.24x |  |
| `grid/linear/4d/max200/eval10/flipped` | construction | SciPy | 105 ms | 21.7 µs | 4842.08x |  |
| `grid/linear/4d/max200/eval1000/flipped` | construction | SciPy | 105 ms | 21.8 µs | 4812.10x |  |
| `grid/linear/4d/max200/eval1000/ascending` | construction | SciPy | 48.8 ms | 10.3 µs | 4741.14x |  |
| `grid/linear/3d/max40/eval100/ascending` | construction | interpn | 4.54 µs | 0.98 ns | 4644.00x | 2e-16 |
| `grid/linear/3d/max40/eval1000/ascending` | construction | interpn | 3.4 µs | 0.88 ns | 3871.11x | 2e-16 |
| `grid/linear/3d/max40/eval10000/ascending` | construction | interpn | 3.14 µs | 0.98 ns | 3219.20x | 3e-16 |
| `grid/linear/3d/max40/eval10/ascending` | construction | interpn | 2.94 µs | 0.98 ns | 3006.80x | 2e-16 |
| `grid/linear/2d/auckland-ba32-0401/heldout/ascending` | construction | SciPy | 17.2 ms | 7.01 µs | 2452.57x |  |
| `grid/nearest/2d/auckland-ba32-0401/heldout/ascending` | construction | SciPy | 16.6 ms | 7.75 µs | 2141.02x |  |
| `grid/linear/2d/auckland-ba32-0401/heldout/ascending` | construction | ndarray-interp | 17.2 ms | 9.05 µs | 1898.71x | 0e+00 |
| `grid/cubic/2d/each10/eval1000/ascending` | construction | interpn | 1.66 µs | 0.88 ns | 1892.44x | 3e-01 |
| `grid/nearest/2d/each40/eval10000/ascending` | construction | interpn | 605 ns | 0.78 ns | 775.00x | 0e+00 |
| `grid/linear/2d/max40/eval10000/ascending` | construction | interpn | 396 ns | 0.68 ns | 578.79x | 4e-16 |
| `grid/linear/2d/max40/eval100/ascending` | construction | interpn | 395 ns | 0.68 ns | 577.57x | 2e-16 |
| `grid/linear/2d/max40/eval1000/ascending` | construction | interpn | 398 ns | 0.78 ns | 509.00x | 3e-16 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | ndarray-interp | 1.07 ms | 2.2 µs | 486.56x | 6e-16 |
| `grid/linear/2d/max40/eval10/ascending` | construction | interpn | 375 ns | 0.78 ns | 480.50x | 2e-16 |
| `grid/linear/2d/max10/eval1000/ascending` | construction | interpn | 265 ns | 0.78 ns | 339.06x | 3e-16 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | ndarray-interp | 309 µs | 1.14 µs | 271.09x | 5e-16 |
| `grid/linear/3d/max200/eval1000/ascending` | construction | SciPy | 1.77 ms | 7.8 µs | 226.32x |  |
| `grid/linear/3d/max200/eval10/ascending` | construction | SciPy | 1.82 ms | 8.23 µs | 221.32x |  |
| `grid/linear/3d/max200/eval100/ascending` | construction | SciPy | 1.8 ms | 8.2 µs | 219.95x |  |
| `grid/linear/3d/max200/eval10000/flipped` | construction | SciPy | 3.23 ms | 15.7 µs | 205.33x |  |
| `grid/linear/3d/max200/eval100/flipped` | construction | SciPy | 3.2 ms | 15.8 µs | 202.88x |  |
| `grid/linear/3d/max10/eval10000/ascending` | construction | interpn | 178 ns | 0.88 ns | 202.11x | 4e-16 |
| `grid/linear/3d/max10/eval100/ascending` | construction | interpn | 196 ns | 0.98 ns | 200.90x | 3e-16 |
| `grid/linear/3d/max10/eval10/ascending` | construction | interpn | 195 ns | 0.98 ns | 199.70x | 2e-16 |
| `grid/linear/3d/max200/eval10000/ascending` | construction | SciPy | 1.63 ms | 8.32 µs | 195.59x |  |
| `grid/linear/3d/max200/eval10/flipped` | construction | SciPy | 3.03 ms | 16.1 µs | 188.66x |  |
| `grid/linear/3d/max10/eval1000/ascending` | construction | interpn | 181 ns | 0.98 ns | 185.50x | 5e-16 |
| `grid/linear/3d/max200/eval1000/flipped` | construction | SciPy | 3.19 ms | 17.4 µs | 183.38x |  |
| `grid/linear/2d/max10/eval10000/ascending` | construction | interpn | 137 ns | 0.78 ns | 175.38x | 3e-16 |
| `grid/linear/2d/max10/eval100/ascending` | construction | interpn | 135 ns | 0.78 ns | 172.75x | 2e-16 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | SciPy | 1.07 ms | 6.42 µs | 166.80x |  |
| `grid/linear/2d/max10/eval10/ascending` | construction | interpn | 136 ns | 0.88 ns | 154.33x | 1e-16 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | SciPy | 309 µs | 5.97 µs | 51.73x |  |
| `grid/pchip/3d/each40/eval10000/ascending` | construction | SciPy | 471 µs | 9.38 µs | 50.23x |  |
| `grid/linear/2d/max200/eval10/ascending` | construction | ndarray-interp | 7.57 µs | 515 ns | 14.69x | 2e-16 |
| `grid/linear/2d/max200/eval10000/ascending` | construction | ndarray-interp | 7.54 µs | 515 ns | 14.64x | 5e-16 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | ndarray-interp | 7.42 µs | 516 ns | 14.38x | 4e-16 |
| `grid/linear/2d/max200/eval100/ascending` | construction | ndarray-interp | 7.41 µs | 515 ns | 14.38x | 3e-16 |
| `oned/akima/n10000/shuffled` | construction | scirs2-interpolate | 141 µs | 45.1 µs | 3.13x | 1e-16 |
| `linear_nd/n10/image100` | construction | spade | 5.36 µs | 1.72 µs | 3.11x | 8e-16 |
| `linear_nd/n10/image50` | construction | spade | 5.77 µs | 1.9 µs | 3.04x | 6e-16 |
| `grid/linear/3d/max40/eval100/ascending` | construction | scirs2-interpolate | 4.54 µs | 1.51 µs | 3.01x | 1e-16 |
| `linear_nd/n10/image500` | construction | spade | 5.38 µs | 1.9 µs | 2.83x | 9e-16 |
| `oned/pchip/n10000/shuffled-uneven` | construction | scirs2-interpolate | 133 µs | 47 µs | 2.83x | 3e-16 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | SciPy | 20.7 µs | 7.47 µs | 2.77x |  |
| `grid/linear/4d/max40/eval100/flipped` | construction | SciPy | 54.5 µs | 19.7 µs | 2.77x |  |
| `oned/pchip/n10000/sorted` | construction | scirs2-interpolate | 129 µs | 47 µs | 2.74x | 2e-16 |
| `oned/pchip/n10000/shuffled` | construction | scirs2-interpolate | 128 µs | 46.9 µs | 2.73x | 2e-16 |
| `grid/linear/2d/max10/eval1000/ascending` | evaluation | interpn | 27.5 µs | 10.1 µs | 2.72x | 3e-16 |
| `oned/akima/n10000/sorted` | construction | scirs2-interpolate | 129 µs | 48.5 µs | 2.66x | 1e-16 |
| `grid/linear/4d/max40/eval10/ascending` | construction | scirs2-interpolate | 14.4 µs | 5.45 µs | 2.64x | 2e-16 |
| `grid/linear/4d/max40/eval10000/flipped` | construction | SciPy | 56.1 µs | 21.6 µs | 2.61x |  |
| `oned/pchip/n10000/sorted-uneven` | construction | scirs2-interpolate | 122 µs | 47.2 µs | 2.59x | 3e-16 |
| `grid/linear/4d/max40/eval1000/flipped` | construction | SciPy | 56.1 µs | 21.7 µs | 2.58x |  |
| `grid/linear/4d/max40/eval10/flipped` | construction | SciPy | 53.2 µs | 21.5 µs | 2.48x |  |
| `grid/linear/2d/max200/eval10/ascending` | construction | scirs2-interpolate | 7.57 µs | 3.09 µs | 2.45x | 1e-16 |
| `grid/linear/4d/max40/eval10000/ascending` | construction | scirs2-interpolate | 14.5 µs | 6.03 µs | 2.40x | 4e-16 |
| `grid/linear/4d/max40/eval100/ascending` | construction | scirs2-interpolate | 14.2 µs | 5.96 µs | 2.39x | 3e-16 |
| `grid/linear/2d/max40/eval1000/ascending` | construction | ndarray-interp | 398 ns | 167 ns | 2.39x | 4e-16 |
| `grid/linear/4d/max40/eval1000/ascending` | construction | scirs2-interpolate | 14.2 µs | 5.99 µs | 2.38x | 2e-16 |
| `grid/linear/2d/max40/eval10000/ascending` | construction | ndarray-interp | 396 ns | 166 ns | 2.38x | 4e-16 |
| `grid/linear/2d/max40/eval100/ascending` | construction | ndarray-interp | 395 ns | 167 ns | 2.37x | 2e-16 |
| `grid/linear/2d/max200/eval100/ascending` | construction | scirs2-interpolate | 7.41 µs | 3.2 µs | 2.32x | 2e-16 |
| `grid/linear/3d/max40/eval1000/ascending` | construction | scirs2-interpolate | 3.4 µs | 1.48 µs | 2.30x | 2e-16 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | scirs2-interpolate | 7.42 µs | 3.23 µs | 2.30x | 2e-16 |
| `grid/linear/2d/max40/eval10/ascending` | construction | ndarray-interp | 375 ns | 166 ns | 2.26x | 2e-16 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | scirs2-interpolate | 20.7 µs | 9.25 µs | 2.24x | 0e+00 |
| `oned/linear/n10/shuffled` | construction | ndarray-interp | 95.9 ns | 45 ns | 2.13x | 0e+00 |
| `linear_nd/n1000/image100` | construction | spade | 582 µs | 277 µs | 2.10x | 1e-15 |
| `grid/linear/2d/max200/eval10000/ascending` | construction | scirs2-interpolate | 7.54 µs | 3.59 µs | 2.10x | 2e-16 |
| `grid/nearest/2d/each40/eval10000/ascending` | construction | scirs2-interpolate | 605 ns | 291 ns | 2.08x | 0e+00 |
| `oned/linear/n10/sorted` | construction | ndarray-interp | 94.5 ns | 45.5 ns | 2.08x | 0e+00 |
| `linear_nd/n1000/image50` | construction | spade | 570 µs | 275 µs | 2.07x | 1e-15 |
| `grid/linear/3d/max40/eval10000/ascending` | construction | scirs2-interpolate | 3.14 µs | 1.54 µs | 2.04x | 2e-16 |
| `grid/linear/2d/max40/eval1000/ascending` | construction | scirs2-interpolate | 398 ns | 198 ns | 2.01x | 3e-16 |
| `linear_nd/n1000/image500` | construction | spade | 568 µs | 284 µs | 2.00x | 3e-15 |
| `linear_nd/n100/image50` | construction | spade | 37.6 µs | 18.9 µs | 1.98x | 1e-15 |
| `grid/linear/3d/max40/eval10/ascending` | construction | scirs2-interpolate | 2.94 µs | 1.49 µs | 1.97x | 1e-16 |
| `linear_nd/n100/image100` | construction | spade | 37.4 µs | 19.1 µs | 1.96x | 9e-16 |
| `linear_nd/n10000/image500` | construction | spade | 7.64 ms | 3.99 ms | 1.91x | 3e-15 |
| `grid/linear/2d/max40/eval100/ascending` | construction | scirs2-interpolate | 395 ns | 207 ns | 1.90x | 3e-16 |
| `oned/akima/n1000/sorted` | construction | scirs2-interpolate | 9.14 µs | 4.84 µs | 1.89x | 3e-17 |
| `oned/akima/n1000/shuffled` | construction | scirs2-interpolate | 9.06 µs | 4.84 µs | 1.87x | 3e-17 |
| `oned/pchip/n1000/sorted` | construction | scirs2-interpolate | 8.53 µs | 4.58 µs | 1.86x | 2e-16 |
| `linear_nd/n10000/image50` | construction | spade | 7.37 ms | 3.96 ms | 1.86x | 8e-16 |
| `oned/pchip/n1000/shuffled` | construction | scirs2-interpolate | 8.32 µs | 4.49 µs | 1.85x | 2e-16 |
| `grid/linear/2d/max40/eval10000/ascending` | construction | scirs2-interpolate | 396 ns | 217 ns | 1.82x | 2e-16 |
| `linear_nd/n10000/image100` | construction | spade | 7.22 ms | 4 ms | 1.80x | 8e-16 |
| `oned/linear/n10000/shuffled-uneven` | construction | scirs2-interpolate | 22.7 µs | 12.9 µs | 1.76x | 6e-17 |
| `grid/linear/2d/max40/eval10/ascending` | construction | scirs2-interpolate | 375 ns | 225 ns | 1.67x | 4e-17 |
| `oned/akima/n100/shuffled` | construction | scirs2-interpolate | 967 ns | 579 ns | 1.67x | 1e-16 |
| `linear_nd/n100/image500` | construction | spade | 34.1 µs | 20.6 µs | 1.65x | 1e-15 |
| `oned/linear/n10000/shuffled` | construction | scirs2-interpolate | 21.3 µs | 12.9 µs | 1.65x | 0e+00 |
| `oned/akima/n100/sorted` | construction | scirs2-interpolate | 962 ns | 585 ns | 1.64x | 1e-16 |
| `oned/pchip/n100/sorted` | construction | scirs2-interpolate | 961 ns | 587 ns | 1.64x | 2e-16 |
| `oned/linear/n10000/sorted-uneven` | construction | scirs2-interpolate | 21.3 µs | 13 µs | 1.63x | 6e-17 |
| `oned/pchip/n100/sorted-uneven` | construction | scirs2-interpolate | 962 ns | 591 ns | 1.63x | 2e-16 |
| `oned/pchip/n100/shuffled` | construction | scirs2-interpolate | 961 ns | 593 ns | 1.62x | 2e-16 |
| `oned/linear/n10000/sorted` | construction | scirs2-interpolate | 20.8 µs | 12.9 µs | 1.62x | 0e+00 |
| `oned/pchip/n50/shuffled` | construction | scirs2-interpolate | 541 ns | 335 ns | 1.62x | 1e-16 |
| `grid/linear/2d/max10/eval10000/ascending` | construction | ndarray-interp | 137 ns | 85.5 ns | 1.60x | 4e-16 |
| `grid/linear/2d/max10/eval1000/ascending` | construction | ndarray-interp | 265 ns | 166 ns | 1.60x | 3e-16 |
| `linear_nd/n100000/image100` | construction | spade | 73.2 ms | 46.2 ms | 1.58x | 7e-16 |
| `oned/pchip/n100/shuffled-uneven` | construction | scirs2-interpolate | 962 ns | 609 ns | 1.58x | 2e-16 |
| `grid/linear/2d/max10/eval100/ascending` | construction | ndarray-interp | 135 ns | 85.4 ns | 1.58x | 4e-16 |
| `grid/linear/2d/max10/eval10/ascending` | construction | ndarray-interp | 136 ns | 85.8 ns | 1.58x | 1e-16 |
| `oned/linear/n100/shuffled` | construction | ndarray-interp | 263 ns | 166 ns | 1.58x | 0e+00 |
| `oned/linear/n1000/shuffled` | construction | scirs2-interpolate | 1.99 µs | 1.27 µs | 1.57x | 0e+00 |
| `oned/linear/n10000/shuffled-uneven` | construction | ndarray-interp | 22.7 µs | 14.8 µs | 1.53x | 0e+00 |
| `linear_nd/n100000/image50` | construction | spade | 78.6 ms | 51.7 ms | 1.52x | 8e-16 |
| `grid/linear/4d/max40/eval10/ascending` | construction | SciPy | 14.4 µs | 9.46 µs | 1.52x |  |
| `linear_nd/n100000/image500` | construction | spade | 115 ms | 75.9 ms | 1.51x | 6e-15 |
| `grid/linear/4d/max40/eval100/ascending` | construction | SciPy | 14.2 µs | 9.44 µs | 1.51x |  |
| `oned/linear/n1000/sorted` | construction | scirs2-interpolate | 1.93 µs | 1.28 µs | 1.51x | 0e+00 |
| `oned/akima/n50/shuffled` | construction | scirs2-interpolate | 496 ns | 332 ns | 1.50x | 1e-16 |
| `oned/akima/n50/sorted` | construction | scirs2-interpolate | 497 ns | 336 ns | 1.48x | 1e-16 |
| `oned/linear/n50/sorted` | construction | ndarray-interp | 161 ns | 109 ns | 1.47x | 0e+00 |
| `oned/linear/n50/shuffled` | construction | ndarray-interp | 159 ns | 110 ns | 1.45x | 0e+00 |
| `oned/linear/n100/shuffled-uneven` | construction | ndarray-interp | 265 ns | 182 ns | 1.45x | 0e+00 |
| `oned/linear/n100/sorted-uneven` | construction | ndarray-interp | 264 ns | 183 ns | 1.45x | 0e+00 |
| `oned/linear/n10000/sorted-uneven` | construction | ndarray-interp | 21.3 µs | 14.8 µs | 1.44x | 0e+00 |
| `oned/linear/n1000/shuffled` | construction | ndarray-interp | 1.99 µs | 1.38 µs | 1.44x | 0e+00 |
| `linear_nd/ghcn-tmax-20250101/heldout` | construction | spade | 11.1 ms | 7.88 ms | 1.41x | 7e-16 |
| `grid/linear/4d/max40/eval10000/ascending` | construction | SciPy | 14.5 µs | 10.3 µs | 1.41x |  |
| `grid/linear/2d/max10/eval10/ascending` | evaluation | interpn | 147 ns | 104 ns | 1.41x | 1e-16 |
| `grid/linear/4d/max40/eval1000/ascending` | construction | SciPy | 14.2 µs | 10.2 µs | 1.40x |  |
| `grid/linear/3d/max200/eval1000/ascending` | construction | scirs2-interpolate | 1.77 ms | 1.27 ms | 1.39x | 3e-16 |
| `grid/linear/3d/max10/eval10/ascending` | evaluation | interpn | 215 ns | 155 ns | 1.39x | 2e-16 |
| `grid/linear/3d/max200/eval100/ascending` | construction | scirs2-interpolate | 1.8 ms | 1.3 ms | 1.39x | 3e-16 |
| `oned/pchip/n50/sorted` | construction | scirs2-interpolate | 500 ns | 361 ns | 1.38x | 1e-16 |
| `grid/linear/3d/max10/eval100/ascending` | evaluation | interpn | 2.09 µs | 1.53 µs | 1.36x | 3e-16 |
| `oned/linear/n100/shuffled` | construction | scirs2-interpolate | 263 ns | 193 ns | 1.36x | 0e+00 |
| `grid/linear/3d/max200/eval10/ascending` | construction | scirs2-interpolate | 1.82 ms | 1.34 ms | 1.36x | 3e-16 |
| `grid/linear/2d/max10/eval10000/ascending` | evaluation | interpn | 137 µs | 101 µs | 1.35x | 3e-16 |
| `grid/linear/2d/max10/eval100/ascending` | evaluation | interpn | 1.36 µs | 1.02 µs | 1.34x | 2e-16 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | scirs2-interpolate | 1.07 ms | 810 µs | 1.32x | 2e-16 |
| `oned/linear/n100/sorted` | construction | ndarray-interp | 240 ns | 182 ns | 1.32x | 0e+00 |
| `grid/linear/4d/max200/eval10/ascending` | construction | scirs2-interpolate | 51.1 ms | 38.9 ms | 1.31x | 2e-16 |
| `grid/linear/2d/max200/eval10000/ascending` | construction | SciPy | 7.54 µs | 5.78 µs | 1.31x |  |
| `grid/linear/2d/auckland-ba32-0401/heldout/ascending` | construction | scirs2-interpolate | 17.2 ms | 13.3 ms | 1.29x | 0e+00 |
| `grid/linear/3d/max10/eval10000/ascending` | evaluation | interpn | 199 µs | 155 µs | 1.29x | 4e-16 |
| `grid/nearest/2d/auckland-ba32-0401/heldout/ascending` | construction | scirs2-interpolate | 16.6 ms | 13 ms | 1.27x | 0e+00 |
| `oned/linear/n100/sorted-uneven` | construction | scirs2-interpolate | 264 ns | 207 ns | 1.27x | 3e-17 |
| `oned/linear/n1000/sorted` | construction | ndarray-interp | 1.93 µs | 1.52 µs | 1.27x | 0e+00 |
| `oned/linear/n100/shuffled-uneven` | construction | scirs2-interpolate | 265 ns | 210 ns | 1.26x | 3e-17 |
| `grid/linear/2d/max40/eval10000/ascending` | evaluation | interpn | 166 µs | 133 µs | 1.25x | 4e-16 |
| `grid/linear/3d/max10/eval1000/ascending` | evaluation | interpn | 19 µs | 15.3 µs | 1.24x | 5e-16 |
| `oned/linear/n50/sorted` | construction | scirs2-interpolate | 161 ns | 131 ns | 1.23x | 0e+00 |
| `grid/linear/2d/max200/eval10000/flipped` | construction | SciPy | 16 µs | 12.9 µs | 1.23x |  |
| `oned/linear/n50/shuffled` | construction | scirs2-interpolate | 159 ns | 130 ns | 1.23x | 0e+00 |
| `oned/cubic_spline/n10000/sorted-uneven` | construction | scirs2-interpolate | 179 µs | 146 µs | 1.23x | 7e-06 |
| `grid/linear/3d/max200/eval10000/ascending` | construction | scirs2-interpolate | 1.63 ms | 1.33 ms | 1.22x | 3e-16 |
| `grid/linear/2d/max200/eval1000/flipped` | construction | SciPy | 15.8 µs | 13 µs | 1.22x |  |
| `grid/linear/2d/max200/eval10/flipped` | construction | SciPy | 15.9 µs | 13 µs | 1.22x |  |
| `grid/linear/4d/max200/eval1000/ascending` | construction | scirs2-interpolate | 48.8 ms | 40.3 ms | 1.21x | 3e-16 |
| `oned/cubic_spline/n10000/sorted` | construction | scirs2-interpolate | 176 µs | 145 µs | 1.21x | 4e-05 |
| `grid/linear/2d/max200/eval100/flipped` | construction | SciPy | 15.4 µs | 12.8 µs | 1.20x |  |
| `grid/linear/2d/max200/eval10/ascending` | construction | SciPy | 7.57 µs | 6.3 µs | 1.20x |  |
| `grid/linear/2d/each400/eval10000/ascending` | construction | scirs2-interpolate | 309 µs | 258 µs | 1.20x | 2e-16 |
| `grid/linear/4d/max200/eval10000/ascending` | construction | scirs2-interpolate | 46.9 ms | 39.6 ms | 1.18x | 4e-16 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | SciPy | 7.42 µs | 6.3 µs | 1.18x |  |
| `oned/cubic_spline/n10000/shuffled-uneven` | construction | scirs2-interpolate | 171 µs | 145 µs | 1.18x | 7e-06 |
| `oned/linear/n10/shuffled` | construction | scirs2-interpolate | 95.9 ns | 81.5 ns | 1.18x | 0e+00 |
| `oned/cubic_spline/n10000/shuffled` | construction | scirs2-interpolate | 171 µs | 145 µs | 1.18x | 4e-05 |
| `grid/linear/3d/max40/eval1000/ascending` | evaluation | interpn | 23.4 µs | 20 µs | 1.17x | 2e-16 |
| `grid/linear/3d/max10/eval10/ascending` | construction | scirs2-interpolate | 195 ns | 167 ns | 1.17x | 2e-16 |
| `oned/linear/n10/sorted` | construction | scirs2-interpolate | 94.5 ns | 81 ns | 1.17x | 0e+00 |
| `oned/nearest/n10/shuffled` | evaluation | scirs2-interpolate | 323 ns | 277 ns | 1.17x | 0e+00 |
| `grid/linear/2d/max200/eval100/ascending` | construction | SciPy | 7.41 µs | 6.36 µs | 1.17x |  |
| `grid/linear/3d/max200/eval10/ascending` | evaluation | interpn | 303 ns | 261 ns | 1.16x | 3e-16 |
| `oned/linear/n100/sorted` | construction | scirs2-interpolate | 240 ns | 207 ns | 1.16x | 0e+00 |
| `grid/linear/2d/max40/eval100/ascending` | evaluation | interpn | 1.66 µs | 1.45 µs | 1.15x | 2e-16 |
| `grid/linear/2d/max40/eval1000/ascending` | evaluation | interpn | 16.5 µs | 14.4 µs | 1.15x | 3e-16 |
| `grid/linear/4d/max40/eval10/ascending` | evaluation | interpn | 338 ns | 296 ns | 1.14x | 2e-16 |
| `grid/linear/3d/max200/eval100/ascending` | evaluation | interpn | 2.98 µs | 2.61 µs | 1.14x | 3e-16 |
| `oned/cubic_spline/n1000/shuffled` | construction | scirs2-interpolate | 15.6 µs | 13.8 µs | 1.13x | 4e-04 |
| `grid/linear/4d/max200/eval1000/ascending` | evaluation | interpn | 41.1 µs | 36.4 µs | 1.13x | 3e-16 |
| `grid/nearest/2d/each40/eval10000/ascending` | evaluation | interpn | 164 µs | 148 µs | 1.11x | 0e+00 |
| `grid/linear/4d/max200/eval100/ascending` | construction | scirs2-interpolate | 47.5 ms | 42.7 ms | 1.11x | 1e-16 |
| `grid/linear/3d/max40/eval10/ascending` | evaluation | interpn | 220 ns | 199 ns | 1.10x | 2e-16 |

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

A run (one library on one case, a process of its own) is stopped at 15 minutes. The longest run that finished took 341 s, 38% of that: scirs2-interpolate on `rbf/quintic/global/ghcn-tmax-20250101/heldout`.

No run took more than half the limit (450 s).

## Real-world data

Two data sets, each taken whole in the unit its source gives it out in; what was taken and how it was divided was fixed before anything was measured (`compare/perf/real.py`, docs/design.ja.md section 8.4.7).

- `auckland-ba32-0401`: Auckland Part 1 LiDAR 1m DEM (2024), tile BA32_10000_0401. 4800 by 7200 cells of 1 m. The grid is every second row and column (and the last ones); every other cell is a query, and its elevation is the value held out. Licensed by Regional Software Holdings Limited for re-use under CC-BY-4.0. Hosted by Toitū Te Whenua Land Information New Zealand (https://registry.opendata.aws/nz-elevation).
- `ghcn-tmax-20250101`: GHCN-Daily, maximum temperature of 2025-01-01, every station; longitude and latitude in degrees, as they are. Of every 10 stations in the order of the file the last is a query, and its temperature is the value held out. The RBF kernels are those without a shape parameter. NOAA Global Historical Climatology Network Daily (GHCN-D), from https://registry.opendata.aws/noaa-ghcn (CC0 1.0).

These cases are also in the tables of their families below and in the counts above. Here they are with the error: the largest and the root mean square distance from the values held out (metres, and degrees Celsius), and the number of queries a library gave no value. A library missing from a case does not have it, or did not finish (see the runs that failed).

`grid/nearest/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 16.6 ms | 369 ms | 15.5 | 0.178 | 0 |
| interpx, parallel | 16.3 ms | 64.1 ms | 15.5 | 0.178 | 0 |
| interpx, `from_owned` | 4.98 ms | 363 ms | 15.5 | 0.178 | 0 |
| SciPy | 7.75 µs | 1.37 s | 15.5 | 0.178 | 0 |
| scirs2-interpolate | 13 ms | 158 s | 15.5 | 0.178 | 0 |
| interpn | 0.78 ns | 770 ms | 15.5 | 0.178 | 0 |

`grid/linear/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 17.2 ms | 399 ms | 6.74 | 0.0601 | 0 |
| interpx, parallel | 16.9 ms | 75 ms | 6.74 | 0.0601 | 0 |
| interpx, `from_owned` | 5.77 ms | 401 ms | 6.74 | 0.0601 | 0 |
| SciPy | 7.01 µs | 1.16 s | 6.74 | 0.0601 | 0 |
| scirs2-interpolate | 13.3 ms | 74.8 s | 6.74 | 0.0601 | 0 |
| ndarray-interp | 9.05 µs | 385 ms | 6.74 | 0.0601 | 0 |
| interpn | 0.78 ns | 866 ms | 6.74 | 0.0601 | 0 |

`grid/cubic/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 135 ms | 788 ms | 6.77 | 0.0507 | 0 |
| interpx, parallel | 78.6 ms | 129 ms | 6.77 | 0.0507 | 0 |
| interpx, `from_owned` | 67.6 ms | 789 ms | 6.77 | 0.0507 | 0 |
| SciPy | 15.5 s | 35.8 s | 6.77 | 0.0507 | 0 |
| interpn | 0.88 ns | 2.11 s | 6.76 | 0.0517 | 0 |

`grid/quintic/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 97.4 ms | 1.29 s | 6.77 | 0.0517 | 0 |
| interpx, parallel | 99.3 ms | 241 ms | 6.77 | 0.0517 | 0 |
| interpx, `from_owned` | 85.4 ms | 1.36 s | 6.77 | 0.0517 | 0 |
| SciPy | 28.2 s | 40.3 s | 6.77 | 0.0517 | 0 |

`grid/pchip/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 96.1 ms | 1.07 s | 6.75 | 0.049 | 0 |
| interpx, parallel | 79.6 ms | 190 ms | 6.75 | 0.049 | 0 |
| interpx, `from_owned` | 70.3 ms | 1.09 s | 6.75 | 0.049 | 0 |

`nearest_nd/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 557 µs | 168 µs | 25 | 3.42 | 0 |
| interpx, parallel | 621 µs | 82.4 µs | 25 | 3.42 | 0 |
| interpx, `from_owned` | 582 µs | 171 µs | 25 | 3.42 | 0 |
| SciPy | 1.81 ms | 520 µs | 25 | 3.42 | 0 |

`linear_nd/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 11.1 ms | 136 µs | 22.1 | 2.83 | 1 |
| interpx, parallel | 9.37 ms | 92 µs | 22.1 | 2.83 | 1 |
| interpx, `from_owned` | 10.3 ms | 158 µs | 22.1 | 2.83 | 1 |
| SciPy | 42.7 ms | 1.32 ms | 22.1 | 2.83 | 1 |
| spade | 7.88 ms | 2.06 ms | 22.1 | 2.83 | 1 |

`rbf/linear/global/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 19 s | 19.2 ms | 20.8 | 2.74 | 0 |
| interpx, parallel | 17.4 s | 3.79 ms | 20.8 | 2.74 | 0 |
| interpx, `from_owned` | 18.8 s | 20.6 ms | 20.8 | 2.74 | 0 |
| SciPy | 128 s | 54.7 ms | 20.8 | 2.74 | 0 |
| scirs2-interpolate | 306 s | 256 ms | 20.8 | 2.74 | 0 |

`rbf/linear/k50/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 754 µs | 19.1 ms | 20.8 | 2.75 | 0 |
| interpx, parallel | 860 µs | 3.72 ms | 20.8 | 2.75 | 0 |
| interpx, `from_owned` | 720 µs | 18.9 ms | 20.8 | 2.75 | 0 |
| SciPy | 1.88 ms | 48.4 ms | 20.8 | 2.75 | 0 |

`rbf/thin_plate_spline/global/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 17.6 s | 80 ms | 102 | 5.67 | 0 |
| interpx, parallel | 17.8 s | 20 ms | 102 | 5.67 | 0 |
| interpx, `from_owned` | 18.2 s | 80 ms | 102 | 5.67 | 0 |
| SciPy | 125 s | 121 ms | 102 | 5.67 | 0 |
| scirs2-interpolate | 331 s | 404 ms | 93.7 | 5.41 | 0 |

`rbf/thin_plate_spline/k50/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 774 µs | 32.3 ms | 102 | 5.66 | 0 |
| interpx, parallel | 801 µs | 5.8 ms | 102 | 5.66 | 0 |
| interpx, `from_owned` | 728 µs | 30.7 ms | 102 | 5.66 | 0 |
| SciPy | 1.85 ms | 57.8 ms | 102 | 5.66 | 0 |

`rbf/cubic/global/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 18.8 s | 19.8 ms | 587 | 27.1 | 0 |
| interpx, parallel | 17.5 s | 3.62 ms | 587 | 27.1 | 0 |
| interpx, `from_owned` | 18.8 s | 20.2 ms | 587 | 27.1 | 0 |
| SciPy | 123 s | 55.1 ms | 597 | 24.4 | 0 |
| scirs2-interpolate | 336 s | 282 ms | 100 | 8.48 | 0 |

`rbf/cubic/k50/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 715 µs | 21.5 ms | 567 | 25.3 | 0 |
| interpx, parallel | 1.46 ms | 4.84 ms | 567 | 25.3 | 0 |
| interpx, `from_owned` | 914 µs | 28.7 ms | 567 | 25.3 | 0 |
| SciPy | 1.97 ms | 51.8 ms | 567 | 25.3 | 0 |

`rbf/quintic/global/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 18.6 s | 22.1 ms | 2.55e+05 | 1.36e+04 | 0 |
| interpx, parallel | 39.6 s | 4.08 ms | 2.55e+05 | 1.36e+04 | 0 |
| interpx, `from_owned` | 17.2 s | 19.5 ms | 2.55e+05 | 1.36e+04 | 0 |
| SciPy | 131 s | 54.9 ms | 9.39e+04 | 4.92e+03 | 0 |
| scirs2-interpolate | 339 s | 339 ms | 6.56e+16 | 2.85e+15 | 0 |

`rbf/quintic/k50/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 858 µs | 23.6 ms | 3.74e+03 | 199 | 0 |
| interpx, parallel | 814 µs | 4.79 ms | 3.74e+03 | 199 | 0 |
| interpx, `from_owned` | 753 µs | 22.4 ms | 3.74e+03 | 199 | 0 |
| SciPy | 1.88 ms | 57.1 ms | 3.74e+03 | 200 | 0 |

## One variable

Data of `n` points, evaluated at `4 n` points: in order (`sorted`) and in no order (`shuffled`).

### Evaluation

| case | interpx | interpx, parallel | SciPy | scirs2-interpolate | ndarray-interp | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/n10/sorted` | 123 ns | 190 ns | 5.54 µs | 276 ns | 280 ns | 2.2x faster |
| `linear/n10/shuffled` | 187 ns | 245 ns | 5.55 µs | 300 ns | 274 ns | 1.5x faster |
| `nearest/n10/sorted` | 151 ns | 210 ns | 9.55 µs | 271 ns |  | 1.8x faster |
| `nearest/n10/shuffled` | 323 ns | 360 ns | 9.46 µs | 277 ns |  | **1.2x slower** |
| `previous/n10/sorted` | 149 ns | 204 ns | 11.6 µs |  |  | 78x faster |
| `previous/n10/shuffled` | 283 ns | 354 ns | 11.7 µs |  |  | 41x faster |
| `bspline1/n10/sorted` | 209 ns | 292 ns | 2.67 µs |  |  | 13x faster |
| `bspline1/n10/shuffled` | 527 ns | 656 ns | 2.96 µs |  |  | 5.6x faster |
| `bspline2/n10/sorted` | 272 ns | 373 ns | 3.12 µs |  |  | 11x faster |
| `bspline2/n10/shuffled` | 372 ns | 454 ns | 3.41 µs |  |  | 9.2x faster |
| `bspline3/n10/sorted` | 378 ns | 474 ns | 3.96 µs |  |  | 10x faster |
| `bspline3/n10/shuffled` | 440 ns | 562 ns | 3.98 µs |  |  | 9.1x faster |
| `cubic_spline/n10/sorted` | 157 ns | 213 ns | 3.68 µs | 300 ns | 364 ns | 1.9x faster |
| `cubic_spline/n10/shuffled` | 204 ns | 265 ns | 4.03 µs | 289 ns | 325 ns | 1.4x faster |
| `pchip/n10/sorted` | 155 ns | 194 ns | 3.69 µs | 319 ns |  | 2.1x faster |
| `pchip/n10/shuffled` | 198 ns | 290 ns | 3.62 µs | 293 ns |  | 1.5x faster |
| `akima/n10/sorted` | 156 ns | 207 ns | 4.07 µs | 262 ns |  | 1.7x faster |
| `akima/n10/shuffled` | 225 ns | 283 ns | 4.04 µs | 268 ns |  | 1.2x faster |
| `makima/n10/sorted` | 143 ns | 212 ns | 3.77 µs |  |  | 26x faster |
| `makima/n10/shuffled` | 205 ns | 256 ns | 3.76 µs |  |  | 18x faster |
| `linear/n50/sorted` | 662 ns | 743 ns | 6.18 µs | 1.65 µs | 1.23 µs | 1.9x faster |
| `linear/n50/shuffled` | 891 ns | 946 ns | 8.2 µs | 1.88 µs | 1.38 µs | 1.6x faster |
| `nearest/n50/sorted` | 787 ns | 851 ns | 12.4 µs | 1.81 µs |  | 2.3x faster |
| `nearest/n50/shuffled` | 1.67 µs | 1.74 µs | 12.2 µs | 1.81 µs |  | 1.1x faster |
| `previous/n50/sorted` | 788 ns | 781 ns | 12.7 µs |  |  | 16x faster |
| `previous/n50/shuffled` | 1.48 µs | 1.55 µs | 13 µs |  |  | 8.8x faster |
| `bspline1/n50/sorted` | 1.89 µs | 1.53 µs | 4.26 µs |  |  | 2.3x faster |
| `bspline1/n50/shuffled` | 2.78 µs | 2.86 µs | 5 µs |  |  | 1.8x faster |
| `bspline2/n50/sorted` | 1.49 µs | 1.76 µs | 6.23 µs |  |  | 4.2x faster |
| `bspline2/n50/shuffled` | 2.22 µs | 2.3 µs | 7.47 µs |  |  | 3.4x faster |
| `bspline3/n50/sorted` | 2.01 µs | 1.91 µs | 9.56 µs |  |  | 4.8x faster |
| `bspline3/n50/shuffled` | 2.89 µs | 3.15 µs | 9.35 µs |  |  | 3.2x faster |
| `cubic_spline/n50/sorted` | 683 ns | 794 ns | 5.83 µs | 4.36 µs | 1.64 µs | 2.4x faster |
| `cubic_spline/n50/shuffled` | 990 ns | 1.04 µs | 5.93 µs | 4.07 µs | 1.6 µs | 1.6x faster |
| `pchip/n50/sorted` | 776 ns | 746 ns | 5.23 µs | 4.27 µs |  | 5.5x faster |
| `pchip/n50/shuffled` | 1.07 µs | 1.03 µs | 5.79 µs | 3.82 µs |  | 3.6x faster |
| `akima/n50/sorted` | 694 ns | 752 ns | 5.65 µs | 4.58 µs |  | 6.6x faster |
| `akima/n50/shuffled` | 983 ns | 1.04 µs | 5.94 µs | 4.03 µs |  | 4.1x faster |
| `makima/n50/sorted` | 687 ns | 824 ns | 5.75 µs |  |  | 8.4x faster |
| `makima/n50/shuffled` | 1.06 µs | 1.14 µs | 6.03 µs |  |  | 5.7x faster |
| `linear/n100/sorted` | 1.25 µs | 1.49 µs | 7.94 µs | 3.69 µs | 2.7 µs | 2.2x faster |
| `linear/n100/shuffled` | 1.62 µs | 1.71 µs | 14 µs | 4.04 µs | 2.5 µs | 1.5x faster |
| `nearest/n100/sorted` | 1.51 µs | 1.72 µs | 13.7 µs | 3.48 µs |  | 2.3x faster |
| `nearest/n100/shuffled` | 3.05 µs | 3.19 µs | 13.6 µs | 3.67 µs |  | 1.2x faster |
| `previous/n100/sorted` | 1.64 µs | 1.72 µs | 15 µs |  |  | 9.1x faster |
| `previous/n100/shuffled` | 2.72 µs | 2.81 µs | 15.1 µs |  |  | 5.5x faster |
| `bspline1/n100/sorted` | 2.17 µs | 2.9 µs | 5.33 µs |  |  | 2.5x faster |
| `bspline1/n100/shuffled` | 6.4 µs | 6.33 µs | 12.8 µs |  |  | 2x faster |
| `bspline2/n100/sorted` | 2.72 µs | 3.7 µs | 9.99 µs |  |  | 3.7x faster |
| `bspline2/n100/shuffled` | 4.44 µs | 5.02 µs | 12.6 µs |  |  | 2.9x faster |
| `bspline3/n100/sorted` | 3.64 µs | 3.79 µs | 16.4 µs |  |  | 4.5x faster |
| `bspline3/n100/shuffled` | 6.19 µs | 6.88 µs | 17.9 µs |  |  | 2.9x faster |
| `cubic_spline/n100/sorted` | 1.39 µs | 1.48 µs | 7.07 µs | 19.1 µs | 3.21 µs | 2.3x faster |
| `cubic_spline/n100/shuffled` | 1.95 µs | 2.08 µs | 8.69 µs | 18.2 µs | 3.2 µs | 1.6x faster |
| `pchip/n100/sorted` | 1.43 µs | 1.63 µs | 7.59 µs | 15.7 µs |  | 5.3x faster |
| `pchip/n100/shuffled` | 2.15 µs | 2.26 µs | 8.61 µs | 17.6 µs |  | 4x faster |
| `akima/n100/sorted` | 1.49 µs | 1.48 µs | 7.52 µs | 17.8 µs |  | 5.1x faster |
| `akima/n100/shuffled` | 2.13 µs | 2.25 µs | 8 µs | 19.5 µs |  | 3.8x faster |
| `makima/n100/sorted` | 1.52 µs | 1.62 µs | 7.11 µs |  |  | 4.7x faster |
| `makima/n100/shuffled` | 2.16 µs | 2.24 µs | 8.71 µs |  |  | 4x faster |
| `linear/n1000/sorted` | 13.4 µs | 13.6 µs | 31.7 µs | 110 µs | 26.9 µs | 2x faster |
| `linear/n1000/shuffled` | 16.2 µs | 17.5 µs | 144 µs | 222 µs | 25.1 µs | 1.6x faster |
| `nearest/n1000/sorted` | 15.3 µs | 15.9 µs | 60.2 µs | 151 µs |  | 3.9x faster |
| `nearest/n1000/shuffled` | 30.8 µs | 34 µs | 221 µs | 245 µs |  | 7.2x faster |
| `previous/n1000/sorted` | 15.4 µs | 16.2 µs | 65.3 µs |  |  | 4.2x faster |
| `previous/n1000/shuffled` | 27.1 µs | 29.8 µs | 220 µs |  |  | 8.1x faster |
| `bspline1/n1000/sorted` | 26.8 µs | 27.9 µs | 31.8 µs |  |  | 1.2x faster |
| `bspline1/n1000/shuffled` | 61.7 µs | 65.2 µs | 578 µs |  |  | 9.4x faster |
| `bspline2/n1000/sorted` | 35.9 µs | 36.1 µs | 83.5 µs |  |  | 2.3x faster |
| `bspline2/n1000/shuffled` | 58.1 µs | 63.9 µs | 592 µs |  |  | 10x faster |
| `bspline3/n1000/sorted` | 36.9 µs | 40.3 µs | 140 µs |  |  | 3.8x faster |
| `bspline3/n1000/shuffled` | 86.7 µs | 79.8 µs | 633 µs |  |  | 7.3x faster |
| `cubic_spline/n1000/sorted` | 14 µs | 14.4 µs | 40.4 µs | 1.55 ms | 34.4 µs | 2.5x faster |
| `cubic_spline/n1000/shuffled` | 21.5 µs | 21.9 µs | 223 µs | 1.55 ms | 31.9 µs | 1.5x faster |
| `pchip/n1000/sorted` | 14.8 µs | 15.2 µs | 41.4 µs | 1.32 ms |  | 2.8x faster |
| `pchip/n1000/shuffled` | 21.5 µs | 20.5 µs | 216 µs | 1.31 ms |  | 10x faster |
| `akima/n1000/sorted` | 14.5 µs | 15.4 µs | 40.7 µs | 1.56 ms |  | 2.8x faster |
| `akima/n1000/shuffled` | 21.6 µs | 21.7 µs | 238 µs | 1.55 ms |  | 11x faster |
| `makima/n1000/sorted` | 14.2 µs | 14.1 µs | 39.8 µs |  |  | 2.8x faster |
| `makima/n1000/shuffled` | 19.7 µs | 20.8 µs | 237 µs |  |  | 12x faster |
| `linear/n10000/sorted` | 141 µs | 71.6 µs | 236 µs | 1.3 ms | 272 µs | 1.7x faster |
| `linear/n10000/shuffled` | 208 µs | 87.8 µs | 2.49 ms | 3.19 ms | 285 µs | 1.4x faster |
| `nearest/n10000/sorted` | 159 µs | 86.6 µs | 636 µs | 1.58 ms |  | 4x faster |
| `nearest/n10000/shuffled` | 414 µs | 103 µs | 3.26 ms | 3.61 ms |  | 7.9x faster |
| `previous/n10000/sorted` | 161 µs | 79.8 µs | 691 µs |  |  | 4.3x faster |
| `previous/n10000/shuffled` | 376 µs | 96.2 µs | 3.39 ms |  |  | 9x faster |
| `bspline1/n10000/sorted` | 296 µs | 104 µs | 312 µs |  |  | 1.1x faster |
| `bspline1/n10000/shuffled` | 700 µs | 149 µs | 51.4 ms |  |  | 73x faster |
| `bspline2/n10000/sorted` | 362 µs | 97.2 µs | 816 µs |  |  | 2.3x faster |
| `bspline2/n10000/shuffled` | 893 µs | 215 µs | 50.7 ms |  |  | 57x faster |
| `bspline3/n10000/sorted` | 737 µs | 196 µs | 1.39 ms |  |  | 1.9x faster |
| `bspline3/n10000/shuffled` | 1.22 ms | 292 µs | 52.2 ms |  |  | 43x faster |
| `cubic_spline/n10000/sorted` | 160 µs | 72.3 µs | 399 µs | 157 ms | 349 µs | 2.2x faster |
| `cubic_spline/n10000/shuffled` | 237 µs | 97.2 µs | 3.67 ms | 143 ms | 364 µs | 1.5x faster |
| `pchip/n10000/sorted` | 146 µs | 86.5 µs | 400 µs | 129 ms |  | 2.7x faster |
| `pchip/n10000/shuffled` | 236 µs | 97.1 µs | 3.67 ms | 130 ms |  | 16x faster |
| `akima/n10000/sorted` | 161 µs | 72.3 µs | 398 µs | 153 ms |  | 2.5x faster |
| `akima/n10000/shuffled` | 236 µs | 96.5 µs | 3.68 ms | 152 ms |  | 16x faster |
| `makima/n10000/sorted` | 154 µs | 81.7 µs | 401 µs |  |  | 2.6x faster |
| `makima/n10000/shuffled` | 236 µs | 91 µs | 3.83 ms |  |  | 16x faster |
| `linear/n100/sorted-uneven` | 1.3 µs | 1.4 µs | 8.1 µs | 3.72 µs | 4.13 µs | 2.9x faster |
| `linear/n100/shuffled-uneven` | 2.83 µs | 3.17 µs | 12.7 µs | 3.74 µs | 4.09 µs | 1.3x faster |
| `cubic_spline/n100/sorted-uneven` | 1.48 µs | 1.59 µs | 8.17 µs | 18.6 µs | 4.86 µs | 3.3x faster |
| `cubic_spline/n100/shuffled-uneven` | 3.2 µs | 3.44 µs | 8.75 µs | 19.7 µs | 4.82 µs | 1.5x faster |
| `pchip/n100/sorted-uneven` | 1.41 µs | 1.55 µs | 8.18 µs | 17.2 µs |  | 5.8x faster |
| `pchip/n100/shuffled-uneven` | 3.21 µs | 3.38 µs | 8.69 µs | 17.5 µs |  | 2.7x faster |
| `linear/n10000/sorted-uneven` | 224 µs | 81.4 µs | 266 µs | 1.39 ms | 980 µs | 1.2x faster |
| `linear/n10000/shuffled-uneven` | 681 µs | 174 µs | 2.52 ms | 3.07 ms | 1.91 ms | 2.8x faster |
| `cubic_spline/n10000/sorted-uneven` | 240 µs | 100 µs | 479 µs | 162 ms | 1.04 ms | 2x faster |
| `cubic_spline/n10000/shuffled-uneven` | 801 µs | 189 µs | 3.72 ms | 155 ms | 2.01 ms | 2.5x faster |
| `pchip/n10000/sorted-uneven` | 240 µs | 94.3 µs | 487 µs | 134 ms |  | 2x faster |
| `pchip/n10000/shuffled-uneven` | 788 µs | 182 µs | 3.59 ms | 132 ms |  | 4.6x faster |

### Construction

| case | interpx | interpx, parallel | SciPy | scirs2-interpolate | ndarray-interp | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/n10/sorted` | 94.5 ns | 90 ns | 9.72 µs | 81 ns | 45.5 ns | **2.1x slower** |
| `linear/n10/shuffled` | 95.9 ns | 88.4 ns | 9.62 µs | 81.5 ns | 45 ns | **2.1x slower** |
| `nearest/n10/sorted` | 75.1 ns | 75.1 ns | 10 µs | 81.7 ns |  | 1.1x faster |
| `nearest/n10/shuffled` | 81 ns | 75.6 ns | 10.4 µs | 87.5 ns |  | 1.1x faster |
| `previous/n10/sorted` | 76 ns | 75.7 ns | 10.1 µs |  |  | 1.3e+02x faster |
| `previous/n10/shuffled` | 81.7 ns | 82.8 ns | 9.55 µs |  |  | 1.2e+02x faster |
| `bspline1/n10/sorted` | 103 ns | 104 ns | 15.8 µs |  |  | 1.5e+02x faster |
| `bspline1/n10/shuffled` | 103 ns | 115 ns | 16 µs |  |  | 1.5e+02x faster |
| `bspline2/n10/sorted` | 470 ns | 513 ns | 37.6 µs |  |  | 80x faster |
| `bspline2/n10/shuffled` | 515 ns | 513 ns | 38.8 µs |  |  | 75x faster |
| `bspline3/n10/sorted` | 566 ns | 566 ns | 39.3 µs |  |  | 69x faster |
| `bspline3/n10/shuffled` | 518 ns | 514 ns | 40.1 µs |  |  | 77x faster |
| `cubic_spline/n10/sorted` | 254 ns | 254 ns | 85.9 µs | 410 ns | 516 ns | 1.6x faster |
| `cubic_spline/n10/shuffled` | 231 ns | 231 ns | 86.8 µs | 390 ns | 473 ns | 1.7x faster |
| `pchip/n10/sorted` | 189 ns | 191 ns | 101 µs | 208 ns |  | 1.1x faster |
| `pchip/n10/shuffled` | 182 ns | 190 ns | 103 µs | 197 ns |  | 1.1x faster |
| `akima/n10/sorted` | 190 ns | 190 ns | 82 µs | 196 ns |  | 1x faster |
| `akima/n10/shuffled` | 190 ns | 192 ns | 83.2 µs | 179 ns |  | **1.1x slower** |
| `makima/n10/sorted` | 188 ns | 201 ns | 80.6 µs |  |  | 4.3e+02x faster |
| `makima/n10/shuffled` | 185 ns | 183 ns | 85.2 µs |  |  | 4.6e+02x faster |
| `linear/n50/sorted` | 161 ns | 166 ns | 9.09 µs | 131 ns | 109 ns | **1.5x slower** |
| `linear/n50/shuffled` | 159 ns | 171 ns | 8.94 µs | 130 ns | 110 ns | **1.5x slower** |
| `nearest/n50/sorted` | 115 ns | 116 ns | 10.9 µs | 134 ns |  | 1.2x faster |
| `nearest/n50/shuffled` | 116 ns | 110 ns | 9.77 µs | 130 ns |  | 1.1x faster |
| `previous/n50/sorted` | 108 ns | 117 ns | 10.3 µs |  |  | 96x faster |
| `previous/n50/shuffled` | 118 ns | 120 ns | 10.5 µs |  |  | 89x faster |
| `bspline1/n50/sorted` | 159 ns | 158 ns | 17.6 µs |  |  | 1.1e+02x faster |
| `bspline1/n50/shuffled` | 144 ns | 145 ns | 15.9 µs |  |  | 1.1e+02x faster |
| `bspline2/n50/sorted` | 2 µs | 1.96 µs | 39.9 µs |  |  | 20x faster |
| `bspline2/n50/shuffled` | 2.01 µs | 2.02 µs | 44.6 µs |  |  | 22x faster |
| `bspline3/n50/sorted` | 2.24 µs | 2.24 µs | 42.4 µs |  |  | 19x faster |
| `bspline3/n50/shuffled` | 2.25 µs | 2.31 µs | 42.3 µs |  |  | 19x faster |
| `cubic_spline/n50/sorted` | 829 ns | 853 ns | 94.7 µs | 861 ns | 1.19 µs | 1x faster |
| `cubic_spline/n50/shuffled` | 832 ns | 914 ns | 94.7 µs | 863 ns | 1.13 µs | 1x faster |
| `pchip/n50/sorted` | 500 ns | 495 ns | 102 µs | 361 ns |  | **1.4x slower** |
| `pchip/n50/shuffled` | 541 ns | 496 ns | 112 µs | 335 ns |  | **1.6x slower** |
| `akima/n50/sorted` | 497 ns | 500 ns | 75.2 µs | 336 ns |  | **1.5x slower** |
| `akima/n50/shuffled` | 496 ns | 499 ns | 82.4 µs | 332 ns |  | **1.5x slower** |
| `makima/n50/sorted` | 554 ns | 596 ns | 84.4 µs |  |  | 1.5e+02x faster |
| `makima/n50/shuffled` | 592 ns | 593 ns | 82.8 µs |  |  | 1.4e+02x faster |
| `linear/n100/sorted` | 240 ns | 263 ns | 10.3 µs | 207 ns | 182 ns | **1.3x slower** |
| `linear/n100/shuffled` | 263 ns | 241 ns | 10.1 µs | 193 ns | 166 ns | **1.6x slower** |
| `nearest/n100/sorted` | 160 ns | 161 ns | 9.98 µs | 195 ns |  | 1.2x faster |
| `nearest/n100/shuffled` | 161 ns | 159 ns | 9.96 µs | 190 ns |  | 1.2x faster |
| `previous/n100/sorted` | 176 ns | 177 ns | 10.1 µs |  |  | 58x faster |
| `previous/n100/shuffled` | 161 ns | 167 ns | 10.5 µs |  |  | 65x faster |
| `bspline1/n100/sorted` | 197 ns | 200 ns | 16 µs |  |  | 82x faster |
| `bspline1/n100/shuffled` | 223 ns | 222 ns | 17.5 µs |  |  | 79x faster |
| `bspline2/n100/sorted` | 3.48 µs | 3.82 µs | 47.3 µs |  |  | 14x faster |
| `bspline2/n100/shuffled` | 3.81 µs | 3.82 µs | 42.7 µs |  |  | 11x faster |
| `bspline3/n100/sorted` | 4.03 µs | 3.92 µs | 45.9 µs |  |  | 11x faster |
| `bspline3/n100/shuffled` | 4.11 µs | 3.98 µs | 43.3 µs |  |  | 11x faster |
| `cubic_spline/n100/sorted` | 1.62 µs | 1.58 µs | 87.7 µs | 1.67 µs | 2.05 µs | 1x faster |
| `cubic_spline/n100/shuffled` | 1.58 µs | 1.61 µs | 96.1 µs | 1.51 µs | 1.88 µs | **1x slower** |
| `pchip/n100/sorted` | 961 ns | 962 ns | 101 µs | 587 ns |  | **1.6x slower** |
| `pchip/n100/shuffled` | 961 ns | 984 ns | 111 µs | 593 ns |  | **1.6x slower** |
| `akima/n100/sorted` | 962 ns | 980 ns | 83 µs | 585 ns |  | **1.6x slower** |
| `akima/n100/shuffled` | 967 ns | 975 ns | 76.2 µs | 579 ns |  | **1.7x slower** |
| `makima/n100/sorted` | 1.06 µs | 1.06 µs | 88 µs |  |  | 83x faster |
| `makima/n100/shuffled` | 1.06 µs | 968 ns | 81.5 µs |  |  | 77x faster |
| `linear/n1000/sorted` | 1.93 µs | 2.12 µs | 13.9 µs | 1.28 µs | 1.52 µs | **1.5x slower** |
| `linear/n1000/shuffled` | 1.99 µs | 2.1 µs | 12.8 µs | 1.27 µs | 1.38 µs | **1.6x slower** |
| `nearest/n1000/sorted` | 1.07 µs | 1.07 µs | 15.6 µs | 1.28 µs |  | 1.2x faster |
| `nearest/n1000/shuffled` | 1.07 µs | 1.06 µs | 15.5 µs | 1.27 µs |  | 1.2x faster |
| `previous/n1000/sorted` | 1.08 µs | 1.08 µs | 17.8 µs |  |  | 16x faster |
| `previous/n1000/shuffled` | 1.06 µs | 1.08 µs | 18.1 µs |  |  | 17x faster |
| `bspline1/n1000/sorted` | 1.04 µs | 1.13 µs | 18.1 µs |  |  | 17x faster |
| `bspline1/n1000/shuffled` | 1.12 µs | 1.14 µs | 17.1 µs |  |  | 15x faster |
| `bspline2/n1000/sorted` | 36.6 µs | 36.5 µs | 118 µs |  |  | 3.2x faster |
| `bspline2/n1000/shuffled` | 36.7 µs | 36.6 µs | 117 µs |  |  | 3.2x faster |
| `bspline3/n1000/sorted` | 42.3 µs | 41.4 µs | 126 µs |  |  | 3x faster |
| `bspline3/n1000/shuffled` | 42.2 µs | 41.6 µs | 134 µs |  |  | 3.2x faster |
| `cubic_spline/n1000/sorted` | 16 µs | 16.9 µs | 121 µs | 14.6 µs | 15.7 µs | **1.1x slower** |
| `cubic_spline/n1000/shuffled` | 15.6 µs | 15.6 µs | 113 µs | 13.8 µs | 16 µs | **1.1x slower** |
| `pchip/n1000/sorted` | 8.53 µs | 8.71 µs | 119 µs | 4.58 µs |  | **1.9x slower** |
| `pchip/n1000/shuffled` | 8.32 µs | 8.99 µs | 127 µs | 4.49 µs |  | **1.9x slower** |
| `akima/n1000/sorted` | 9.14 µs | 9.03 µs | 100 µs | 4.84 µs |  | **1.9x slower** |
| `akima/n1000/shuffled` | 9.06 µs | 8.34 µs | 102 µs | 4.84 µs |  | **1.9x slower** |
| `makima/n1000/sorted` | 9.66 µs | 9.9 µs | 99.8 µs |  |  | 10x faster |
| `makima/n1000/shuffled` | 9.44 µs | 9.16 µs | 98.8 µs |  |  | 10x faster |
| `linear/n10000/sorted` | 20.8 µs | 22 µs | 43.4 µs | 12.9 µs | 21.6 µs | **1.6x slower** |
| `linear/n10000/shuffled` | 21.3 µs | 28.3 µs | 50.8 µs | 12.9 µs | 22.1 µs | **1.6x slower** |
| `nearest/n10000/sorted` | 11.6 µs | 11.6 µs | 60 µs | 13 µs |  | 1.1x faster |
| `nearest/n10000/shuffled` | 11.4 µs | 11.9 µs | 60.4 µs | 12.9 µs |  | 1.1x faster |
| `previous/n10000/sorted` | 11.5 µs | 11.6 µs | 83.7 µs |  |  | 7.3x faster |
| `previous/n10000/shuffled` | 11.4 µs | 11.4 µs | 77 µs |  |  | 6.7x faster |
| `bspline1/n10000/sorted` | 11.5 µs | 11.7 µs | 25.9 µs |  |  | 2.3x faster |
| `bspline1/n10000/shuffled` | 11.5 µs | 11.5 µs | 25.7 µs |  |  | 2.2x faster |
| `bspline2/n10000/sorted` | 509 µs | 438 µs | 774 µs |  |  | 1.5x faster |
| `bspline2/n10000/shuffled` | 511 µs | 516 µs | 810 µs |  |  | 1.6x faster |
| `bspline3/n10000/sorted` | 841 µs | 964 µs | 1.02 ms |  |  | 1.2x faster |
| `bspline3/n10000/shuffled` | 486 µs | 584 µs | 1.02 ms |  |  | 2.1x faster |
| `cubic_spline/n10000/sorted` | 176 µs | 175 µs | 711 µs | 145 µs | 183 µs | **1.2x slower** |
| `cubic_spline/n10000/shuffled` | 171 µs | 171 µs | 756 µs | 145 µs | 183 µs | **1.2x slower** |
| `pchip/n10000/sorted` | 129 µs | 134 µs | 253 µs | 47 µs |  | **2.7x slower** |
| `pchip/n10000/shuffled` | 128 µs | 218 µs | 264 µs | 46.9 µs |  | **2.7x slower** |
| `akima/n10000/sorted` | 129 µs | 91.3 µs | 249 µs | 48.5 µs |  | **2.7x slower** |
| `akima/n10000/shuffled` | 141 µs | 92.7 µs | 673 µs | 45.1 µs |  | **3.1x slower** |
| `makima/n10000/sorted` | 99.9 µs | 229 µs | 718 µs |  |  | 7.2x faster |
| `makima/n10000/shuffled` | 140 µs | 104 µs | 799 µs |  |  | 5.7x faster |
| `linear/n100/sorted-uneven` | 264 ns | 266 ns | 10 µs | 207 ns | 183 ns | **1.4x slower** |
| `linear/n100/shuffled-uneven` | 265 ns | 265 ns | 10 µs | 210 ns | 182 ns | **1.5x slower** |
| `cubic_spline/n100/sorted-uneven` | 1.72 µs | 1.73 µs | 87.4 µs | 1.66 µs | 2.07 µs | **1x slower** |
| `cubic_spline/n100/shuffled-uneven` | 1.59 µs | 1.75 µs | 95.1 µs | 1.65 µs | 1.88 µs | 1x faster |
| `pchip/n100/sorted-uneven` | 962 ns | 879 ns | 113 µs | 591 ns |  | **1.6x slower** |
| `pchip/n100/shuffled-uneven` | 962 ns | 964 ns | 111 µs | 609 ns |  | **1.6x slower** |
| `linear/n10000/sorted-uneven` | 21.3 µs | 21.6 µs | 50.9 µs | 13 µs | 14.8 µs | **1.6x slower** |
| `linear/n10000/shuffled-uneven` | 22.7 µs | 21.3 µs | 49.6 µs | 12.9 µs | 14.8 µs | **1.8x slower** |
| `cubic_spline/n10000/sorted-uneven` | 179 µs | 178 µs | 741 µs | 146 µs | 184 µs | **1.2x slower** |
| `cubic_spline/n10000/shuffled-uneven` | 171 µs | 191 µs | 662 µs | 145 µs | 183 µs | **1.2x slower** |
| `pchip/n10000/sorted-uneven` | 122 µs | 222 µs | 275 µs | 47.2 µs |  | **2.6x slower** |
| `pchip/n10000/shuffled-uneven` | 133 µs | 92.9 µs | 255 µs | 47 µs |  | **2.8x slower** |

## `PPoly`

A cubic of 55 intervals, evaluated at points in no order.

### Evaluation

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `intervals55/eval10` | 74.1 ns | 75.1 ns | 3.8 µs | 51x faster |
| `intervals55/eval100` | 715 ns | 795 ns | 4.91 µs | 6.9x faster |
| `intervals55/eval1000` | 7.05 µs | 7.05 µs | 16.2 µs | 2.3x faster |
| `intervals55/eval10000` | 65.5 µs | 72.1 µs | 361 µs | 5.5x faster |

### Construction

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `intervals55/eval10` | 173 ns | 174 ns | 5.63 µs | 32x faster |
| `intervals55/eval100` | 176 ns | 174 ns | 5.66 µs | 32x faster |
| `intervals55/eval1000` | 173 ns | 158 ns | 5.22 µs | 30x faster |
| `intervals55/eval10000` | 173 ns | 159 ns | 5.56 µs | 32x faster |

## `BSpline`

A cubic B-spline of 55 intervals, evaluated at points in no order.

### Evaluation

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `intervals55/eval10` | 166 ns | 163 ns | 2.95 µs | 18x faster |
| `intervals55/eval100` | 1.55 µs | 1.65 µs | 5.93 µs | 3.8x faster |
| `intervals55/eval1000` | 14.2 µs | 15.7 µs | 38.6 µs | 2.7x faster |
| `intervals55/eval10000` | 156 µs | 64.3 µs | 399 µs | 2.6x faster |

### Construction

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `intervals55/eval10` | 176 ns | 179 ns | 10 µs | 57x faster |
| `intervals55/eval100` | 175 ns | 174 ns | 9.95 µs | 57x faster |
| `intervals55/eval1000` | 175 ns | 161 ns | 9.24 µs | 53x faster |
| `intervals55/eval10000` | 175 ns | 174 ns | 9.2 µs | 53x faster |

## `BSplineCurve`

A cubic curve through `n` points of the plane, evaluated at `4 n` parameters in order.

### Evaluation

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `degree3/n100` | 6.56 µs | 4.7 µs | 17.5 µs | 2.7x faster |
| `degree3/n10000` | 493 µs | 234 µs | 1.75 ms | 3.5x faster |

### Construction

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `degree3/n100` | 4.48 µs | 4.7 µs | 50.4 µs | 11x faster |
| `degree3/n10000` | 682 µs | 1.19 ms | 2.48 ms | 3.6x faster |

## Regular grid

Points in no order inside the grid. `maxN`: the first axis has `N` points and each axis half the one before (SciPy's suite); `eachN`: every axis has `N` points; `flipped`: every axis given from its high end.

### Evaluation

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | ndarray-interp | interpn | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/2d/max10/eval10/ascending` | 147 ns | 149 ns | 139 ns | 21.3 µs | 873 ns | 166 ns | 104 ns | **1.4x slower** |
| `linear/2d/max10/eval10/flipped` | 148 ns | 150 ns | 148 ns | 21.9 µs |  |  |  | 1.5e+02x faster |
| `linear/2d/max10/eval100/ascending` | 1.36 µs | 1.44 µs | 1.37 µs | 23.7 µs | 9.29 µs | 1.49 µs | 1.02 µs | **1.3x slower** |
| `linear/2d/max10/eval100/flipped` | 1.37 µs | 1.44 µs | 1.37 µs | 24.4 µs |  |  |  | 18x faster |
| `linear/2d/max10/eval1000/ascending` | 27.5 µs | 13.8 µs | 13.6 µs | 44.7 µs | 92 µs | 27 µs | 10.1 µs | **2.7x slower** |
| `linear/2d/max10/eval1000/flipped` | 13.8 µs | 13.8 µs | 13.6 µs | 44 µs |  |  |  | 3.2x faster |
| `linear/2d/max10/eval10000/ascending` | 137 µs | 59.5 µs | 136 µs | 395 µs | 963 µs | 240 µs | 101 µs | **1.3x slower** |
| `linear/2d/max10/eval10000/flipped` | 137 µs | 60.4 µs | 136 µs | 392 µs |  |  |  | 2.9x faster |
| `linear/2d/max40/eval10/ascending` | 161 ns | 175 ns | 174 ns | 20.2 µs | 1.11 µs | 171 ns | 147 ns | **1.1x slower** |
| `linear/2d/max40/eval10/flipped` | 172 ns | 174 ns | 172 ns | 21.4 µs |  |  |  | 1.2e+02x faster |
| `linear/2d/max40/eval100/ascending` | 1.66 µs | 1.74 µs | 1.66 µs | 24.5 µs | 11 µs | 1.59 µs | 1.45 µs | **1.1x slower** |
| `linear/2d/max40/eval100/flipped` | 1.66 µs | 1.74 µs | 1.66 µs | 24.3 µs |  |  |  | 15x faster |
| `linear/2d/max40/eval1000/ascending` | 16.5 µs | 16.8 µs | 16.6 µs | 50.7 µs | 111 µs | 16.4 µs | 14.4 µs | **1.1x slower** |
| `linear/2d/max40/eval1000/flipped` | 16.5 µs | 16.7 µs | 16.5 µs | 50.6 µs |  |  |  | 3.1x faster |
| `linear/2d/max40/eval10000/ascending` | 166 µs | 69.6 µs | 151 µs | 563 µs | 1.19 ms | 245 µs | 133 µs | **1.2x slower** |
| `linear/2d/max40/eval10000/flipped` | 154 µs | 65.5 µs | 166 µs | 612 µs |  |  |  | 4x faster |
| `linear/2d/max200/eval10/ascending` | 211 ns | 216 ns | 212 ns | 21.4 µs | 2.15 µs | 206 ns | 193 ns | **1.1x slower** |
| `linear/2d/max200/eval10/flipped` | 212 ns | 198 ns | 212 ns | 21.6 µs |  |  |  | 1e+02x faster |
| `linear/2d/max200/eval100/ascending` | 2.06 µs | 2.15 µs | 2.06 µs | 25.3 µs | 21.7 µs | 1.99 µs | 1.89 µs | **1.1x slower** |
| `linear/2d/max200/eval100/flipped` | 1.95 µs | 2.17 µs | 2.07 µs | 25.1 µs |  |  |  | 13x faster |
| `linear/2d/max200/eval1000/ascending` | 20.7 µs | 21 µs | 20.7 µs | 75.4 µs | 228 µs | 24.9 µs | 18.9 µs | **1.1x slower** |
| `linear/2d/max200/eval1000/flipped` | 20.7 µs | 20.9 µs | 20.7 µs | 74.1 µs |  |  |  | 3.6x faster |
| `linear/2d/max200/eval10000/ascending` | 207 µs | 75.4 µs | 207 µs | 886 µs | 2.14 ms | 489 µs | 189 µs | **1.1x slower** |
| `linear/2d/max200/eval10000/flipped` | 207 µs | 81.1 µs | 207 µs | 820 µs |  |  |  | 4x faster |
| `linear/3d/max10/eval10/ascending` | 215 ns | 218 ns | 216 ns | 62.1 µs | 3.42 µs |  | 155 ns | **1.4x slower** |
| `linear/3d/max10/eval10/flipped` | 214 ns | 204 ns | 216 ns | 56.1 µs |  |  |  | 2.6e+02x faster |
| `linear/3d/max10/eval100/ascending` | 2.09 µs | 2.14 µs | 1.9 µs | 67.6 µs | 33.5 µs |  | 1.53 µs | **1.4x slower** |
| `linear/3d/max10/eval100/flipped` | 2.09 µs | 1.97 µs | 2.09 µs | 65.2 µs |  |  |  | 31x faster |
| `linear/3d/max10/eval1000/ascending` | 19 µs | 20.1 µs | 18.9 µs | 123 µs | 315 µs |  | 15.3 µs | **1.2x slower** |
| `linear/3d/max10/eval1000/flipped` | 19.1 µs | 19.5 µs | 20.8 µs | 115 µs |  |  |  | 6x faster |
| `linear/3d/max10/eval10000/ascending` | 199 µs | 72.7 µs | 207 µs | 1.14 ms | 3.49 ms |  | 155 µs | **1.3x slower** |
| `linear/3d/max10/eval10000/flipped` | 207 µs | 70 µs | 209 µs | 1.17 ms |  |  |  | 5.7x faster |
| `linear/3d/max40/eval10/ascending` | 220 ns | 222 ns | 219 ns | 61.5 µs | 3.61 µs |  | 199 ns | **1.1x slower** |
| `linear/3d/max40/eval10/flipped` | 240 ns | 244 ns | 241 ns | 65 µs |  |  |  | 2.7e+02x faster |
| `linear/3d/max40/eval100/ascending` | 2.33 µs | 2.4 µs | 2.33 µs | 72.2 µs | 36.8 µs |  | 2.17 µs | **1.1x slower** |
| `linear/3d/max40/eval100/flipped` | 2.33 µs | 2.41 µs | 2.33 µs | 70.5 µs |  |  |  | 30x faster |
| `linear/3d/max40/eval1000/ascending` | 23.4 µs | 23.4 µs | 21.5 µs | 145 µs | 368 µs |  | 20 µs | **1.2x slower** |
| `linear/3d/max40/eval1000/flipped` | 21.3 µs | 23.5 µs | 21.3 µs | 143 µs |  |  |  | 6.7x faster |
| `linear/3d/max40/eval10000/ascending` | 213 µs | 82.6 µs | 232 µs | 1.56 ms | 3.83 ms |  | 219 µs | 1x faster |
| `linear/3d/max40/eval10000/flipped` | 218 µs | 90.9 µs | 233 µs | 1.55 ms |  |  |  | 7.1x faster |
| `linear/3d/max200/eval10/ascending` | 303 ns | 308 ns | 277 ns | 61.8 µs | 4.93 µs |  | 261 ns | **1.2x slower** |
| `linear/3d/max200/eval10/flipped` | 303 ns | 306 ns | 276 ns | 61.6 µs |  |  |  | 2e+02x faster |
| `linear/3d/max200/eval100/ascending` | 2.98 µs | 3.07 µs | 2.99 µs | 70.9 µs | 52.3 µs |  | 2.61 µs | **1.1x slower** |
| `linear/3d/max200/eval100/flipped` | 2.99 µs | 3.06 µs | 2.98 µs | 69.9 µs |  |  |  | 23x faster |
| `linear/3d/max200/eval1000/ascending` | 29.8 µs | 30.2 µs | 29.9 µs | 190 µs | 557 µs |  | 28.9 µs | **1x slower** |
| `linear/3d/max200/eval1000/flipped` | 29.9 µs | 29.9 µs | 27.2 µs | 205 µs |  |  |  | 6.9x faster |
| `linear/3d/max200/eval10000/ascending` | 272 µs | 95 µs | 297 µs | 1.93 ms | 5.15 ms |  | 264 µs | **1x slower** |
| `linear/3d/max200/eval10000/flipped` | 299 µs | 93.5 µs | 297 µs | 1.88 ms |  |  |  | 6.3x faster |
| `linear/4d/max10/eval10/ascending` | 273 ns | 270 ns | 264 ns | 95.1 µs |  |  |  | 3.5e+02x faster |
| `linear/4d/max10/eval10/flipped` | 275 ns | 267 ns | 261 ns | 103 µs |  |  |  | 3.8e+02x faster |
| `linear/4d/max10/eval100/ascending` | 2.41 µs | 2.51 µs | 2.41 µs | 108 µs |  |  |  | 45x faster |
| `linear/4d/max10/eval100/flipped` | 2.42 µs | 2.59 µs | 2.64 µs | 116 µs |  |  |  | 48x faster |
| `linear/4d/max10/eval1000/ascending` | 26.7 µs | 24.8 µs | 24.6 µs | 210 µs |  |  |  | 7.9x faster |
| `linear/4d/max10/eval1000/flipped` | 24.4 µs | 24.6 µs | 26.6 µs | 209 µs |  |  |  | 8.5x faster |
| `linear/4d/max10/eval10000/ascending` | 243 µs | 88.6 µs | 243 µs | 1.82 ms |  |  |  | 7.5x faster |
| `linear/4d/max10/eval10000/flipped` | 267 µs | 89.5 µs | 267 µs | 1.91 ms |  |  |  | 7.2x faster |
| `linear/4d/max40/eval10/ascending` | 338 ns | 313 ns | 310 ns | 104 µs | 5.91 µs |  | 296 ns | **1.1x slower** |
| `linear/4d/max40/eval10/flipped` | 309 ns | 341 ns | 338 ns | 96 µs |  |  |  | 3.1e+02x faster |
| `linear/4d/max40/eval100/ascending` | 3.19 µs | 3.02 µs | 3.22 µs | 112 µs | 59.6 µs |  | 2.97 µs | **1.1x slower** |
| `linear/4d/max40/eval100/flipped` | 2.96 µs | 3.3 µs | 3.22 µs | 108 µs |  |  |  | 36x faster |
| `linear/4d/max40/eval1000/ascending` | 32.4 µs | 29.9 µs | 32.1 µs | 265 µs | 632 µs |  | 29.6 µs | **1.1x slower** |
| `linear/4d/max40/eval1000/flipped` | 29.4 µs | 32.6 µs | 31.7 µs | 266 µs |  |  |  | 9.1x faster |
| `linear/4d/max40/eval10000/ascending` | 324 µs | 99 µs | 328 µs | 2.46 ms | 6.62 ms |  | 299 µs | **1.1x slower** |
| `linear/4d/max40/eval10000/flipped` | 297 µs | 104 µs | 300 µs | 2.2 ms |  |  |  | 7.4x faster |
| `linear/4d/max200/eval10/ascending` | 392 ns | 419 ns | 385 ns | 95.1 µs | 8.63 µs |  | 392 ns | **1x slower** |
| `linear/4d/max200/eval10/flipped` | 386 ns | 424 ns | 387 ns | 95.3 µs |  |  |  | 2.5e+02x faster |
| `linear/4d/max200/eval100/ascending` | 3.75 µs | 4.22 µs | 3.75 µs | 109 µs | 80.1 µs |  | 3.96 µs | 1.1x faster |
| `linear/4d/max200/eval100/flipped` | 4.09 µs | 3.82 µs | 4.09 µs | 121 µs |  |  |  | 30x faster |
| `linear/4d/max200/eval1000/ascending` | 41.1 µs | 42.1 µs | 37.5 µs | 337 µs | 789 µs |  | 36.4 µs | **1.1x slower** |
| `linear/4d/max200/eval1000/flipped` | 37.8 µs | 42.1 µs | 38.1 µs | 325 µs |  |  |  | 8.6x faster |
| `linear/4d/max200/eval10000/ascending` | 411 µs | 117 µs | 408 µs | 2.99 ms | 7.53 ms |  | 400 µs | **1x slower** |
| `linear/4d/max200/eval10000/flipped` | 375 µs | 116 µs | 401 µs | 3.08 ms |  |  |  | 8.2x faster |
| `nearest/2d/each40/eval10000/ascending` | 164 µs | 57.3 µs | 150 µs | 694 µs | 1.6 ms |  | 148 µs | **1.1x slower** |
| `nearest/3d/each40/eval10000/ascending` | 212 µs | 71 µs | 233 µs | 946 µs | 1.83 ms |  | 208 µs | **1x slower** |
| `pchip/2d/each40/eval10000/ascending` | 544 µs | 135 µs | 543 µs | 1.29 s |  |  |  | 2.4e+03x faster |
| `pchip/3d/each40/eval10000/ascending` | 1.61 ms | 278 µs | 1.61 ms | 3.49 s |  |  |  | 2.2e+03x faster |
| `cubic/2d/each40/eval10000/ascending` | 300 µs | 104 µs | 300 µs | 1.28 ms |  |  | 409 µs | 1.4x faster |
| `cubic/3d/each40/eval10000/ascending` | 663 µs | 174 µs | 724 µs | 3.17 ms |  |  | 856 µs | 1.3x faster |
| `cubic/2d/each10/eval1000/ascending` | 29.1 µs | 29.3 µs | 29.1 µs | 131 µs |  |  | 39.2 µs | 1.3x faster |
| `cubic/2d/each40/eval1000/ascending` | 29.9 µs | 30.2 µs | 32.9 µs | 143 µs |  |  | 44.6 µs | 1.5x faster |
| `cubic/2d/each100/eval1000/ascending` | 32.1 µs | 32.2 µs | 32.1 µs | 180 µs |  |  | 47.1 µs | 1.5x faster |
| `cubic/2d/each200/eval1000/ascending` | 34.6 µs | 34.8 µs | 37.9 µs | 198 µs |  |  | 45.9 µs | 1.3x faster |
| `cubic/2d/each400/eval1000/ascending` | 38.6 µs | 41.8 µs | 38.4 µs | 273 µs |  |  | 49.6 µs | 1.3x faster |
| `quintic/2d/each10/eval1000/ascending` | 48.4 µs | 48.7 µs | 48.4 µs | 227 µs |  |  |  | 4.7x faster |
| `quintic/2d/each40/eval1000/ascending` | 52.6 µs | 48.8 µs | 48.1 µs | 260 µs |  |  |  | 4.9x faster |
| `linear/2d/each400/eval10000/ascending` | 241 µs | 75.4 µs | 220 µs | 949 µs | 4.17 ms | 538 µs | 219 µs | **1.1x slower** |
| `linear/2d/each800/eval10000/ascending` | 239 µs | 83.1 µs | 262 µs | 1.09 ms | 8.14 ms | 595 µs | 241 µs | 1x faster |
| `nearest/2d/auckland-ba32-0401/heldout/ascending` | 369 ms | 64.1 ms | 363 ms | 1.37 s | 158 s |  | 770 ms | 2.1x faster |
| `linear/2d/auckland-ba32-0401/heldout/ascending` | 399 ms | 75 ms | 401 ms | 1.16 s | 74.8 s | 385 ms | 866 ms | **1x slower** |
| `cubic/2d/auckland-ba32-0401/heldout/ascending` | 788 ms | 129 ms | 789 ms | 35.8 s |  |  | 2.11 s | 2.7x faster |
| `quintic/2d/auckland-ba32-0401/heldout/ascending` | 1.29 s | 241 ms | 1.36 s | 40.3 s |  |  |  | 31x faster |
| `pchip/2d/auckland-ba32-0401/heldout/ascending` | 1.07 s | 190 ms | 1.09 s |  |  |  |  |  |

### Construction

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | ndarray-interp | interpn | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/2d/max10/eval10/ascending` | 136 ns | 136 ns | 133 ns | 6.3 µs | 151 ns | 85.8 ns | 0.88 ns | **1.5e+02x slower** |
| `linear/2d/max10/eval10/flipped` | 178 ns | 183 ns | 190 ns | 12.7 µs |  |  |  | 71x faster |
| `linear/2d/max10/eval100/ascending` | 135 ns | 137 ns | 141 ns | 6.33 µs | 146 ns | 85.4 ns | 0.78 ns | **1.7e+02x slower** |
| `linear/2d/max10/eval100/flipped` | 179 ns | 180 ns | 184 ns | 12.7 µs |  |  |  | 71x faster |
| `linear/2d/max10/eval1000/ascending` | 265 ns | 140 ns | 145 ns | 6.37 µs | 145 ns | 166 ns | 0.78 ns | **3.4e+02x slower** |
| `linear/2d/max10/eval1000/flipped` | 178 ns | 182 ns | 190 ns | 12.1 µs |  |  |  | 68x faster |
| `linear/2d/max10/eval10000/ascending` | 137 ns | 139 ns | 154 ns | 6.37 µs | 147 ns | 85.5 ns | 0.78 ns | **1.8e+02x slower** |
| `linear/2d/max10/eval10000/flipped` | 179 ns | 164 ns | 182 ns | 12.6 µs |  |  |  | 70x faster |
| `linear/2d/max40/eval10/ascending` | 375 ns | 399 ns | 819 ns | 6.25 µs | 225 ns | 166 ns | 0.78 ns | **4.8e+02x slower** |
| `linear/2d/max40/eval10/flipped` | 729 ns | 720 ns | 1.14 µs | 12 µs |  |  |  | 17x faster |
| `linear/2d/max40/eval100/ascending` | 395 ns | 413 ns | 868 ns | 6.3 µs | 207 ns | 167 ns | 0.68 ns | **5.8e+02x slower** |
| `linear/2d/max40/eval100/flipped` | 725 ns | 728 ns | 1.19 µs | 12.7 µs |  |  |  | 18x faster |
| `linear/2d/max40/eval1000/ascending` | 398 ns | 399 ns | 858 ns | 5.82 µs | 198 ns | 167 ns | 0.78 ns | **5.1e+02x slower** |
| `linear/2d/max40/eval1000/flipped` | 717 ns | 717 ns | 1.09 µs | 12.8 µs |  |  |  | 18x faster |
| `linear/2d/max40/eval10000/ascending` | 396 ns | 401 ns | 784 ns | 5.78 µs | 217 ns | 166 ns | 0.68 ns | **5.8e+02x slower** |
| `linear/2d/max40/eval10000/flipped` | 661 ns | 713 ns | 1.21 µs | 11.8 µs |  |  |  | 18x faster |
| `linear/2d/max200/eval10/ascending` | 7.57 µs | 7.67 µs | 18.4 µs | 6.3 µs | 3.09 µs | 515 ns | 0.78 ns | **9.7e+03x slower** |
| `linear/2d/max200/eval10/flipped` | 15.9 µs | 14.7 µs | 24.1 µs | 13 µs |  |  |  | **1.2x slower** |
| `linear/2d/max200/eval100/ascending` | 7.41 µs | 7.66 µs | 19.3 µs | 6.36 µs | 3.2 µs | 515 ns | 0.78 ns | **9.5e+03x slower** |
| `linear/2d/max200/eval100/flipped` | 15.4 µs | 15.7 µs | 25.7 µs | 12.8 µs |  |  |  | **1.2x slower** |
| `linear/2d/max200/eval1000/ascending` | 7.42 µs | 7.63 µs | 19.2 µs | 6.3 µs | 3.23 µs | 516 ns | 0.78 ns | **9.5e+03x slower** |
| `linear/2d/max200/eval1000/flipped` | 15.8 µs | 15.6 µs | 26.1 µs | 13 µs |  |  |  | **1.2x slower** |
| `linear/2d/max200/eval10000/ascending` | 7.54 µs | 7.63 µs | 4.82 µs | 5.78 µs | 3.59 µs | 515 ns | 0.78 ns | **9.7e+03x slower** |
| `linear/2d/max200/eval10000/flipped` | 16 µs | 15.9 µs | 20.1 µs | 12.9 µs |  |  |  | **1.2x slower** |
| `linear/3d/max10/eval10/ascending` | 195 ns | 192 ns | 224 ns | 8.11 µs | 167 ns |  | 0.98 ns | **2e+02x slower** |
| `linear/3d/max10/eval10/flipped` | 323 ns | 310 ns | 357 ns | 15.4 µs |  |  |  | 48x faster |
| `linear/3d/max10/eval100/ascending` | 196 ns | 179 ns | 225 ns | 8.15 µs | 181 ns |  | 0.98 ns | **2e+02x slower** |
| `linear/3d/max10/eval100/flipped` | 327 ns | 302 ns | 368 ns | 17 µs |  |  |  | 52x faster |
| `linear/3d/max10/eval1000/ascending` | 181 ns | 200 ns | 234 ns | 7.56 µs | 180 ns |  | 0.98 ns | **1.9e+02x slower** |
| `linear/3d/max10/eval1000/flipped` | 321 ns | 321 ns | 359 ns | 15.6 µs |  |  |  | 49x faster |
| `linear/3d/max10/eval10000/ascending` | 178 ns | 192 ns | 226 ns | 7.69 µs | 182 ns |  | 0.88 ns | **2e+02x slower** |
| `linear/3d/max10/eval10000/flipped` | 322 ns | 322 ns | 363 ns | 15.5 µs |  |  |  | 48x faster |
| `linear/3d/max40/eval10/ascending` | 2.94 µs | 6.4 µs | 7.7 µs | 8.2 µs | 1.49 µs |  | 0.98 ns | **3e+03x slower** |
| `linear/3d/max40/eval10/flipped` | 8.88 µs | 8.98 µs | 13.2 µs | 15.6 µs |  |  |  | 1.8x faster |
| `linear/3d/max40/eval100/ascending` | 4.54 µs | 3.12 µs | 8.04 µs | 10.8 µs | 1.51 µs |  | 0.98 ns | **4.6e+03x slower** |
| `linear/3d/max40/eval100/flipped` | 13.4 µs | 15.6 µs | 11.9 µs | 17.5 µs |  |  |  | 1.3x faster |
| `linear/3d/max40/eval1000/ascending` | 3.4 µs | 2.99 µs | 3.13 µs | 8.47 µs | 1.48 µs |  | 0.88 ns | **3.9e+03x slower** |
| `linear/3d/max40/eval1000/flipped` | 8.1 µs | 8.9 µs | 9.33 µs | 17.3 µs |  |  |  | 2.1x faster |
| `linear/3d/max40/eval10000/ascending` | 3.14 µs | 3.13 µs | 1.95 µs | 8.22 µs | 1.54 µs |  | 0.98 ns | **3.2e+03x slower** |
| `linear/3d/max40/eval10000/flipped` | 8.82 µs | 8.97 µs | 10.2 µs | 16.2 µs |  |  |  | 1.8x faster |
| `linear/3d/max200/eval10/ascending` | 1.82 ms | 1.79 ms | 631 µs | 8.23 µs | 1.34 ms |  | 0.88 ns | **2.1e+06x slower** |
| `linear/3d/max200/eval10/flipped` | 3.03 ms | 3.47 ms | 1.93 ms | 16.1 µs |  |  |  | **1.9e+02x slower** |
| `linear/3d/max200/eval100/ascending` | 1.8 ms | 1.83 ms | 695 µs | 8.2 µs | 1.3 ms |  | 0.88 ns | **2.1e+06x slower** |
| `linear/3d/max200/eval100/flipped` | 3.2 ms | 3.71 ms | 1.92 ms | 15.8 µs |  |  |  | **2e+02x slower** |
| `linear/3d/max200/eval1000/ascending` | 1.77 ms | 1.86 ms | 674 µs | 7.8 µs | 1.27 ms |  | 0.98 ns | **1.8e+06x slower** |
| `linear/3d/max200/eval1000/flipped` | 3.19 ms | 3.09 ms | 1.8 ms | 17.4 µs |  |  |  | **1.8e+02x slower** |
| `linear/3d/max200/eval10000/ascending` | 1.63 ms | 1.79 ms | 643 µs | 8.32 µs | 1.33 ms |  | 0.88 ns | **1.9e+06x slower** |
| `linear/3d/max200/eval10000/flipped` | 3.23 ms | 3.32 ms | 1.92 ms | 15.7 µs |  |  |  | **2.1e+02x slower** |
| `linear/4d/max10/eval10/ascending` | 236 ns | 237 ns | 248 ns | 9.7 µs |  |  |  | 41x faster |
| `linear/4d/max10/eval10/flipped` | 349 ns | 371 ns | 402 ns | 17.9 µs |  |  |  | 51x faster |
| `linear/4d/max10/eval100/ascending` | 224 ns | 217 ns | 248 ns | 9.84 µs |  |  |  | 44x faster |
| `linear/4d/max10/eval100/flipped` | 356 ns | 342 ns | 401 ns | 18.5 µs |  |  |  | 52x faster |
| `linear/4d/max10/eval1000/ascending` | 238 ns | 217 ns | 249 ns | 9.81 µs |  |  |  | 41x faster |
| `linear/4d/max10/eval1000/flipped` | 352 ns | 350 ns | 421 ns | 17.8 µs |  |  |  | 51x faster |
| `linear/4d/max10/eval10000/ascending` | 238 ns | 221 ns | 251 ns | 9.75 µs |  |  |  | 41x faster |
| `linear/4d/max10/eval10000/flipped` | 348 ns | 354 ns | 404 ns | 17.8 µs |  |  |  | 51x faster |
| `linear/4d/max40/eval10/ascending` | 14.4 µs | 13.5 µs | 32.7 µs | 9.46 µs | 5.45 µs |  | 1.37 ns | **1.1e+04x slower** |
| `linear/4d/max40/eval10/flipped` | 53.2 µs | 56.6 µs | 70.2 µs | 21.5 µs |  |  |  | **2.5x slower** |
| `linear/4d/max40/eval100/ascending` | 14.2 µs | 13.4 µs | 35.4 µs | 9.44 µs | 5.96 µs |  | 1.37 ns | **1e+04x slower** |
| `linear/4d/max40/eval100/flipped` | 54.5 µs | 57.3 µs | 67.2 µs | 19.7 µs |  |  |  | **2.8x slower** |
| `linear/4d/max40/eval1000/ascending` | 14.2 µs | 14.3 µs | 28.2 µs | 10.2 µs | 5.99 µs |  | 1.27 ns | **1.1e+04x slower** |
| `linear/4d/max40/eval1000/flipped` | 56.1 µs | 58.2 µs | 68.3 µs | 21.7 µs |  |  |  | **2.6x slower** |
| `linear/4d/max40/eval10000/ascending` | 14.5 µs | 14.1 µs | 32.2 µs | 10.3 µs | 6.03 µs |  | 1.37 ns | **1.1e+04x slower** |
| `linear/4d/max40/eval10000/flipped` | 56.1 µs | 57 µs | 70.2 µs | 21.6 µs |  |  |  | **2.6x slower** |
| `linear/4d/max200/eval10/ascending` | 51.1 ms | 48 ms | 15.4 ms | 10.1 µs | 38.9 ms |  | 1.37 ns | **3.7e+07x slower** |
| `linear/4d/max200/eval10/flipped` | 105 ms | 105 ms | 70.2 ms | 21.7 µs |  |  |  | **4.8e+03x slower** |
| `linear/4d/max200/eval100/ascending` | 47.5 ms | 51.1 ms | 15.4 ms | 9.38 µs | 42.7 ms |  | 1.37 ns | **3.5e+07x slower** |
| `linear/4d/max200/eval100/flipped` | 108 ms | 108 ms | 71.6 ms | 22 µs |  |  |  | **4.9e+03x slower** |
| `linear/4d/max200/eval1000/ascending` | 48.8 ms | 50.9 ms | 16 ms | 10.3 µs | 40.3 ms |  | 1.27 ns | **3.8e+07x slower** |
| `linear/4d/max200/eval1000/flipped` | 105 ms | 103 ms | 72.6 ms | 21.8 µs |  |  |  | **4.8e+03x slower** |
| `linear/4d/max200/eval10000/ascending` | 46.9 ms | 52.4 ms | 15.9 ms | 9.54 µs | 39.6 ms |  | 1.27 ns | **3.7e+07x slower** |
| `linear/4d/max200/eval10000/flipped` | 106 ms | 106 ms | 68.8 ms | 21.7 µs |  |  |  | **4.9e+03x slower** |
| `nearest/2d/each40/eval10000/ascending` | 605 ns | 587 ns | 1.56 µs | 6.36 µs | 291 ns |  | 0.78 ns | **7.8e+02x slower** |
| `nearest/3d/each40/eval10000/ascending` | 20.7 µs | 23 µs | 59.6 µs | 7.47 µs | 9.25 µs |  | 0.88 ns | **2.4e+04x slower** |
| `pchip/2d/each40/eval10000/ascending` | 7.13 µs | 6.74 µs | 6.9 µs | 7.16 µs |  |  |  | 1x faster |
| `pchip/3d/each40/eval10000/ascending` | 471 µs | 542 µs | 487 µs | 9.38 µs |  |  |  | **50x slower** |
| `cubic/2d/each40/eval10000/ascending` | 8.39 µs | 9.15 µs | 8.14 µs | 1.79 ms |  |  | 0.88 ns | **9.6e+03x slower** |
| `cubic/3d/each40/eval10000/ascending` | 368 µs | 290 µs | 290 µs | 147 ms |  |  | 1.17 ns | **3.1e+05x slower** |
| `cubic/2d/each10/eval1000/ascending` | 1.66 µs | 1.65 µs | 1.58 µs | 1.04 ms |  |  | 0.88 ns | **1.9e+03x slower** |
| `cubic/2d/each40/eval1000/ascending` | 8.39 µs | 8.91 µs | 8.93 µs | 1.92 ms |  |  | 0.88 ns | **9.6e+03x slower** |
| `cubic/2d/each100/eval1000/ascending` | 41 µs | 50.4 µs | 39.6 µs | 7.82 ms |  |  | 0.88 ns | **4.7e+04x slower** |
| `cubic/2d/each200/eval1000/ascending` | 174 µs | 278 µs | 191 µs | 25.9 ms |  |  | 0.88 ns | **2e+05x slower** |
| `cubic/2d/each400/eval1000/ascending` | 1.33 ms | 1.42 ms | 1.1 ms | 145 ms |  |  | 0.88 ns | **1.5e+06x slower** |
| `quintic/2d/each10/eval1000/ascending` | 2.16 µs | 2.14 µs | 2.22 µs | 3.4 ms |  |  |  | 1.6e+03x faster |
| `quintic/2d/each40/eval1000/ascending` | 13.5 µs | 12.4 µs | 12.1 µs | 7.73 ms |  |  |  | 5.7e+02x faster |
| `linear/2d/each400/eval10000/ascending` | 309 µs | 319 µs | 97.8 µs | 5.97 µs | 258 µs | 1.14 µs | 0.68 ns | **4.5e+05x slower** |
| `linear/2d/each800/eval10000/ascending` | 1.07 ms | 1.1 ms | 398 µs | 6.42 µs | 810 µs | 2.2 µs | 0.68 ns | **1.6e+06x slower** |
| `nearest/2d/auckland-ba32-0401/heldout/ascending` | 16.6 ms | 16.3 ms | 4.98 ms | 7.75 µs | 13 ms |  | 0.78 ns | **2.1e+07x slower** |
| `linear/2d/auckland-ba32-0401/heldout/ascending` | 17.2 ms | 16.9 ms | 5.77 ms | 7.01 µs | 13.3 ms | 9.05 µs | 0.78 ns | **2.2e+07x slower** |
| `cubic/2d/auckland-ba32-0401/heldout/ascending` | 135 ms | 78.6 ms | 67.6 ms | 15.5 s |  |  | 0.88 ns | **1.5e+08x slower** |
| `quintic/2d/auckland-ba32-0401/heldout/ascending` | 97.4 ms | 99.3 ms | 85.4 ms | 28.2 s |  |  |  | 2.9e+02x faster |
| `pchip/2d/auckland-ba32-0401/heldout/ascending` | 96.1 ms | 79.6 ms | 70.3 ms |  |  |  |  |  |

## `NearestNDInterpolator`

`n` points of the unit square, evaluated at the pixels of a square image.

### Evaluation

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: |
| `n10/image50` | 35.3 µs | 32.6 µs | 114 µs | 907 µs | 26x faster |
| `n10/image100` | 125 µs | 51.4 µs | 451 µs | 1.73 ms | 14x faster |
| `n10/image500` | 2.93 ms | 547 µs | 11.6 ms | 44.3 ms | 15x faster |
| `n100/image50` | 112 µs | 112 µs | 202 µs | 568 µs | 5.1x faster |
| `n100/image100` | 415 µs | 114 µs | 713 µs | 2.31 ms | 5.6x faster |
| `n100/image500` | 9.47 ms | 1.58 ms | 18.2 ms | 56.6 ms | 6x faster |
| `n1000/image50` | 305 µs | 109 µs | 306 µs | 792 µs | 2.6x faster |
| `n1000/image100` | 1.14 ms | 290 µs | 1.14 ms | 2.97 ms | 2.6x faster |
| `n1000/image500` | 23 ms | 4.13 ms | 23.1 ms | 68.4 ms | 3x faster |
| `n10000/image50` | 320 µs | 117 µs | 327 µs | 983 µs | 3.1x faster |
| `n10000/image100` | 1.34 ms | 287 µs | 1.34 ms | 3.83 ms | 2.8x faster |
| `n10000/image500` | 29.7 ms | 5.09 ms | 29.5 ms | 91.1 ms | 3.1x faster |
| `n100000/image50` | 430 µs | 111 µs | 430 µs | 1.45 ms | 3.4x faster |
| `n100000/image100` | 1.95 ms | 351 µs | 1.91 ms | 5.58 ms | 2.9x faster |
| `n100000/image500` | 42.2 ms | 6.81 ms | 42.1 ms | 178 ms | 4.2x faster |
| `ghcn-tmax-20250101/heldout` | 168 µs | 82.4 µs | 171 µs | 520 µs | 3.1x faster |

### Construction

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: |
| `n10/image50` | 509 ns | 554 ns | 522 ns | 9.86 µs | 19x faster |
| `n10/image100` | 513 ns | 513 ns | 515 ns | 9.63 µs | 19x faster |
| `n10/image500` | 540 ns | 507 ns | 517 ns | 9.61 µs | 18x faster |
| `n100/image50` | 3.56 µs | 3.29 µs | 3.57 µs | 15.5 µs | 4.3x faster |
| `n100/image100` | 3.25 µs | 3.55 µs | 3.46 µs | 15.1 µs | 4.7x faster |
| `n100/image500` | 3.34 µs | 3.27 µs | 3.58 µs | 15.7 µs | 4.7x faster |
| `n1000/image50` | 34.3 µs | 36.8 µs | 34.5 µs | 117 µs | 3.4x faster |
| `n1000/image100` | 33.9 µs | 34.2 µs | 34.4 µs | 117 µs | 3.5x faster |
| `n1000/image500` | 34.6 µs | 35.2 µs | 34.9 µs | 117 µs | 3.4x faster |
| `n10000/image50` | 490 µs | 531 µs | 522 µs | 1.67 ms | 3.4x faster |
| `n10000/image100` | 526 µs | 494 µs | 485 µs | 1.64 ms | 3.1x faster |
| `n10000/image500` | 487 µs | 492 µs | 481 µs | 1.63 ms | 3.3x faster |
| `n100000/image50` | 8.62 ms | 8.14 ms | 8.61 ms | 25.7 ms | 3x faster |
| `n100000/image100` | 9.04 ms | 8.36 ms | 8.4 ms | 24.4 ms | 2.7x faster |
| `n100000/image500` | 8.1 ms | 8.03 ms | 8.71 ms | 24.2 ms | 3x faster |
| `ghcn-tmax-20250101/heldout` | 557 µs | 621 µs | 582 µs | 1.81 ms | 3.3x faster |

## `LinearNDInterpolator`

`n` points of the unit square, evaluated at the pixels of a square image (some outside the hull, with no value); `3d`: `n` points of the unit cube, evaluated at 10 000 points of the cube in no order.

### Evaluation

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | spade | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `n10/image50` | 18.6 µs | 17.3 µs | 17 µs | 29.9 µs | 106 µs | 1.6x faster |
| `n10/image100` | 69.2 µs | 70.5 µs | 69.1 µs | 104 µs | 395 µs | 1.5x faster |
| `n10/image500` | 1.74 ms | 366 µs | 1.77 ms | 2.79 ms | 10.2 ms | 1.6x faster |
| `n100/image50` | 55.6 µs | 56.5 µs | 51.1 µs | 72.9 µs | 206 µs | 1.3x faster |
| `n100/image100` | 218 µs | 76.5 µs | 218 µs | 253 µs | 805 µs | 1.2x faster |
| `n100/image500` | 5.17 ms | 934 µs | 5.23 ms | 6.2 ms | 16.9 ms | 1.2x faster |
| `n1000/image50` | 54.1 µs | 55.3 µs | 59.2 µs | 122 µs | 274 µs | 2.3x faster |
| `n1000/image100` | 276 µs | 95.6 µs | 276 µs | 392 µs | 1.11 ms | 1.4x faster |
| `n1000/image500` | 6.22 ms | 1.12 ms | 6.29 ms | 6.45 ms | 21.5 ms | 1x faster |
| `n10000/image50` | 65.8 µs | 61.7 µs | 62.5 µs | 414 µs | 562 µs | 6.3x faster |
| `n10000/image100` | 320 µs | 98.3 µs | 322 µs | 1 ms | 1.72 ms | 3.1x faster |
| `n10000/image500` | 8.19 ms | 1.33 ms | 8.04 ms | 9 ms | 28.1 ms | 1.1x faster |
| `n100000/image50` | 87.3 µs | 98.7 µs | 88.6 µs | 1.28 ms | 3.13 ms | 15x faster |
| `n100000/image100` | 456 µs | 128 µs | 423 µs | 3.45 ms | 6.75 ms | 7.6x faster |
| `n100000/image500` | 12.8 ms | 2.07 ms | 11.4 ms | 22 ms | 111 ms | 1.7x faster |
| `3d/n1000/scattered` | 773 µs | 196 µs | 762 µs | 5.76 ms |  | 7.5x faster |
| `3d/n10000/scattered` | 1.28 ms | 239 µs | 1.18 ms | 15.8 ms |  | 12x faster |
| `3d/n100000/scattered` | 1.91 ms | 531 µs | 2.1 ms | 80.7 ms |  | 42x faster |
| `ghcn-tmax-20250101/heldout` | 136 µs | 92 µs | 158 µs | 1.32 ms | 2.06 ms | 9.6x faster |

### Construction

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | spade | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `n10/image50` | 5.77 µs | 5.23 µs | 5.23 µs | 454 µs | 1.9 µs | **3x slower** |
| `n10/image100` | 5.36 µs | 5.24 µs | 5.2 µs | 463 µs | 1.72 µs | **3.1x slower** |
| `n10/image500` | 5.38 µs | 4.98 µs | 5.23 µs | 453 µs | 1.9 µs | **2.8x slower** |
| `n100/image50` | 37.6 µs | 34.5 µs | 36 µs | 638 µs | 18.9 µs | **2x slower** |
| `n100/image100` | 37.4 µs | 33.2 µs | 37.1 µs | 628 µs | 19.1 µs | **2x slower** |
| `n100/image500` | 34.1 µs | 35.8 µs | 34.8 µs | 672 µs | 20.6 µs | **1.7x slower** |
| `n1000/image50` | 570 µs | 570 µs | 631 µs | 2.23 ms | 275 µs | **2.1x slower** |
| `n1000/image100` | 582 µs | 570 µs | 572 µs | 2.46 ms | 277 µs | **2.1x slower** |
| `n1000/image500` | 568 µs | 577 µs | 564 µs | 2.44 ms | 284 µs | **2x slower** |
| `n10000/image50` | 7.37 ms | 7.27 ms | 7.31 ms | 23.2 ms | 3.96 ms | **1.9x slower** |
| `n10000/image100` | 7.22 ms | 7.35 ms | 7.66 ms | 25 ms | 4 ms | **1.8x slower** |
| `n10000/image500` | 7.64 ms | 7.51 ms | 7.3 ms | 25.1 ms | 3.99 ms | **1.9x slower** |
| `n100000/image50` | 78.6 ms | 75.8 ms | 76.9 ms | 376 ms | 51.7 ms | **1.5x slower** |
| `n100000/image100` | 73.2 ms | 75 ms | 75.9 ms | 363 ms | 46.2 ms | **1.6x slower** |
| `n100000/image500` | 115 ms | 79.3 ms | 79.7 ms | 395 ms | 75.9 ms | **1.5x slower** |
| `3d/n1000/scattered` | 4.35 ms | 4.29 ms | 4.03 ms | 9.78 ms |  | 2.3x faster |
| `3d/n10000/scattered` | 47.1 ms | 47.4 ms | 47.4 ms | 131 ms |  | 2.8x faster |
| `3d/n100000/scattered` | 511 ms | 514 ms | 518 ms | 1.92 s |  | 3.7x faster |
| `ghcn-tmax-20250101/heldout` | 11.1 ms | 9.37 ms | 10.3 ms | 42.7 ms | 7.88 ms | **1.4x slower** |

## `RBFInterpolator`

`n` points of a square, evaluated at 1000 points in no order (`scattered`) or at the 10 000 pixels of an image (`image`). `global`: all points; `k50`: the 50 nearest.

### Evaluation

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/global/n10/scattered` | 15.8 µs | 15.9 µs | 15.8 µs | 64.4 µs | 199 µs | 4.1x faster |
| `linear/k50/n10/scattered` | 15.8 µs | 15.9 µs | 15.8 µs | 1.36 ms |  | 86x faster |
| `thin_plate_spline/global/n10/scattered` | 71.4 µs | 76.4 µs | 71.8 µs | 121 µs | 295 µs | 1.7x faster |
| `thin_plate_spline/k50/n10/scattered` | 77.9 µs | 75 µs | 72.8 µs | 1.53 ms |  | 20x faster |
| `cubic/global/n10/scattered` | 16.9 µs | 16.8 µs | 16.9 µs | 60.8 µs | 215 µs | 3.6x faster |
| `cubic/k50/n10/scattered` | 17.6 µs | 17 µs | 16.9 µs | 1.45 ms |  | 83x faster |
| `quintic/global/n10/scattered` | 23.2 µs | 23.5 µs | 25.9 µs | 67.4 µs | 213 µs | 2.9x faster |
| `quintic/k50/n10/scattered` | 25.9 µs | 25.1 µs | 23.7 µs | 1.46 ms |  | 56x faster |
| `multiquadric/global/n10/scattered` | 27.8 µs | 27.9 µs | 27.8 µs | 52.8 µs | 259 µs | 1.9x faster |
| `multiquadric/k50/n10/scattered` | 27.8 µs | 27.9 µs | 27.8 µs | 1.47 ms |  | 53x faster |
| `inverse_multiquadric/global/n10/scattered` | 41.1 µs | 37.6 µs | 41.1 µs | 54.8 µs | 273 µs | 1.3x faster |
| `inverse_multiquadric/k50/n10/scattered` | 37.4 µs | 37.6 µs | 41.1 µs | 1.35 ms |  | 36x faster |
| `inverse_quadratic/global/n10/scattered` | 23.9 µs | 26.3 µs | 23.8 µs | 48 µs |  | 2x faster |
| `inverse_quadratic/k50/n10/scattered` | 23.9 µs | 26.3 µs | 24 µs | 1.4 ms |  | 58x faster |
| `gaussian/global/n10/scattered` | 73.7 µs | 73.9 µs | 73.7 µs | 96 µs | 391 µs | 1.3x faster |
| `gaussian/k50/n10/scattered` | 73.7 µs | 73.7 µs | 73.7 µs | 1.37 ms |  | 19x faster |
| `linear/global/n100/scattered` | 135 µs | 47.3 µs | 137 µs | 305 µs | 1.89 ms | 2.3x faster |
| `linear/k50/n100/scattered` | 12.9 ms | 2.54 ms | 13.9 ms | 19 ms |  | 1.5x faster |
| `thin_plate_spline/global/n100/scattered` | 701 µs | 188 µs | 699 µs | 943 µs | 3.17 ms | 1.3x faster |
| `thin_plate_spline/k50/n100/scattered` | 24.2 ms | 4.21 ms | 24.9 ms | 23.9 ms |  | **1x slower** |
| `cubic/global/n100/scattered` | 152 µs | 66 µs | 142 µs | 319 µs | 2.07 ms | 2.1x faster |
| `cubic/k50/n100/scattered` | 15.7 ms | 3.09 ms | 14.8 ms | 20.2 ms |  | 1.3x faster |
| `quintic/global/n100/scattered` | 174 µs | 64.7 µs | 180 µs | 302 µs | 2.26 ms | 1.7x faster |
| `quintic/k50/n100/scattered` | 17.9 ms | 3.34 ms | 17.5 ms | 23.4 ms |  | 1.3x faster |
| `multiquadric/global/n100/scattered` | 295 µs | 84.8 µs | 295 µs | 331 µs | 2.47 ms | 1.1x faster |
| `multiquadric/k50/n100/scattered` | 16.1 ms | 3.08 ms | 15.9 ms | 20.5 ms |  | 1.3x faster |
| `inverse_multiquadric/global/n100/scattered` | 399 µs | 111 µs | 399 µs | 403 µs | 2.81 ms | 1x faster |
| `inverse_multiquadric/k50/n100/scattered` | 16.8 ms | 3.17 ms | 16.8 ms | 20.6 ms |  | 1.2x faster |
| `inverse_quadratic/global/n100/scattered` | 224 µs | 75.1 µs | 246 µs | 288 µs |  | 1.3x faster |
| `inverse_quadratic/k50/n100/scattered` | 13.7 ms | 3.31 ms | 14.7 ms | 20.5 ms |  | 1.5x faster |
| `gaussian/global/n100/scattered` | 675 µs | 205 µs | 679 µs | 699 µs | 3.97 ms | 1x faster |
| `gaussian/k50/n100/scattered` | 23.3 ms | 3.9 ms | 21.2 ms | 23 ms |  | **1x slower** |
| `linear/global/n1000/scattered` | 1.47 ms | 340 µs | 1.47 ms | 4.15 ms | 19.2 ms | 2.8x faster |
| `linear/k50/n1000/scattered` | 14.6 ms | 2.74 ms | 13.6 ms | 36.1 ms |  | 2.5x faster |
| `thin_plate_spline/global/n1000/scattered` | 7.66 ms | 1.55 ms | 8.03 ms | 10.9 ms | 32.2 ms | 1.4x faster |
| `thin_plate_spline/k50/n1000/scattered` | 25.4 ms | 4.68 ms | 24.1 ms | 46 ms |  | 1.8x faster |
| `cubic/global/n1000/scattered` | 1.48 ms | 347 µs | 1.38 ms | 4.29 ms | 19.5 ms | 2.9x faster |
| `cubic/k50/n1000/scattered` | 17.7 ms | 3.18 ms | 16.7 ms | 39 ms |  | 2.2x faster |
| `quintic/global/n1000/scattered` | 1.66 ms | 385 µs | 1.76 ms | 4.47 ms | 23.8 ms | 2.7x faster |
| `quintic/k50/n1000/scattered` | 18.6 ms | 3.36 ms | 18.3 ms | 44.3 ms |  | 2.4x faster |
| `multiquadric/global/n1000/scattered` | 2.77 ms | 636 µs | 2.75 ms | 4.5 ms | 25 ms | 1.6x faster |
| `multiquadric/k50/n1000/scattered` | 16.9 ms | 3.26 ms | 16.2 ms | 37.6 ms |  | 2.2x faster |
| `inverse_multiquadric/global/n1000/scattered` | 4.35 ms | 837 µs | 4.17 ms | 5.37 ms | 28.4 ms | 1.2x faster |
| `inverse_multiquadric/k50/n1000/scattered` | 18.1 ms | 3.41 ms | 18.6 ms | 38.8 ms |  | 2.1x faster |
| `inverse_quadratic/global/n1000/scattered` | 2.31 ms | 534 µs | 2.48 ms | 3.92 ms |  | 1.7x faster |
| `inverse_quadratic/k50/n1000/scattered` | 16.1 ms | 3.09 ms | 15.7 ms | 36.9 ms |  | 2.3x faster |
| `gaussian/global/n1000/scattered` | 6.79 ms | 1.46 ms | 6.84 ms | 8.6 ms | 38.6 ms | 1.3x faster |
| `gaussian/k50/n1000/scattered` | 22.1 ms | 4.17 ms | 21.9 ms | 46 ms |  | 2.1x faster |
| `thin_plate_spline/global/n3000/scattered` | 23 ms | 4.32 ms | 22 ms | 32.9 ms | 105 ms | 1.4x faster |
| `thin_plate_spline/k50/n10000/scattered` | 22.1 ms | 4.6 ms | 22.7 ms | 44 ms |  | 2x faster |
| `thin_plate_spline/k50/n10000/image` | 232 ms | 39.3 ms | 233 ms | 460 ms |  | 2x faster |
| `thin_plate_spline/k50/n100000/scattered` | 24.6 ms | 4.7 ms | 24.7 ms | 49.5 ms |  | 2x faster |
| `thin_plate_spline/k50/n100000/image` | 268 ms | 39.4 ms | 241 ms | 492 ms |  | 1.8x faster |
| `linear/global/ghcn-tmax-20250101/heldout` | 19.2 ms | 3.79 ms | 20.6 ms | 54.7 ms | 256 ms | 2.8x faster |
| `linear/k50/ghcn-tmax-20250101/heldout` | 19.1 ms | 3.72 ms | 18.9 ms | 48.4 ms |  | 2.5x faster |
| `thin_plate_spline/global/ghcn-tmax-20250101/heldout` | 80 ms | 20 ms | 80 ms | 121 ms | 404 ms | 1.5x faster |
| `thin_plate_spline/k50/ghcn-tmax-20250101/heldout` | 32.3 ms | 5.8 ms | 30.7 ms | 57.8 ms |  | 1.8x faster |
| `cubic/global/ghcn-tmax-20250101/heldout` | 19.8 ms | 3.62 ms | 20.2 ms | 55.1 ms | 282 ms | 2.8x faster |
| `cubic/k50/ghcn-tmax-20250101/heldout` | 21.5 ms | 4.84 ms | 28.7 ms | 51.8 ms |  | 2.4x faster |
| `quintic/global/ghcn-tmax-20250101/heldout` | 22.1 ms | 4.08 ms | 19.5 ms | 54.9 ms | 339 ms | 2.5x faster |
| `quintic/k50/ghcn-tmax-20250101/heldout` | 23.6 ms | 4.79 ms | 22.4 ms | 57.1 ms |  | 2.4x faster |

### Construction

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/global/n10/scattered` | 2.85 µs | 2.89 µs | 2.85 µs | 11.8 µs | 5.9 µs | 2.1x faster |
| `linear/k50/n10/scattered` | 2.88 µs | 2.86 µs | 2.85 µs | 16.5 µs |  | 5.7x faster |
| `thin_plate_spline/global/n10/scattered` | 3.4 µs | 3.31 µs | 3.61 µs | 13.3 µs | 6.51 µs | 1.9x faster |
| `thin_plate_spline/k50/n10/scattered` | 3.32 µs | 3.32 µs | 3.59 µs | 17.2 µs |  | 5.2x faster |
| `cubic/global/n10/scattered` | 3.46 µs | 3.57 µs | 3.77 µs | 12 µs | 6 µs | 1.7x faster |
| `cubic/k50/n10/scattered` | 3.75 µs | 3.8 µs | 3.79 µs | 17.2 µs |  | 4.6x faster |
| `quintic/global/n10/scattered` | 3.83 µs | 4.25 µs | 4.24 µs | 13.9 µs | 6.21 µs | 1.6x faster |
| `quintic/k50/n10/scattered` | 4.25 µs | 4.25 µs | 3.87 µs | 18.8 µs |  | 4.4x faster |
| `multiquadric/global/n10/scattered` | 2.67 µs | 2.67 µs | 2.68 µs | 10.7 µs | 7.16 µs | 2.7x faster |
| `multiquadric/k50/n10/scattered` | 2.71 µs | 2.69 µs | 2.68 µs | 15.1 µs |  | 5.6x faster |
| `inverse_multiquadric/global/n10/scattered` | 2.93 µs | 2.68 µs | 2.68 µs | 11.7 µs | 7.48 µs | 2.5x faster |
| `inverse_multiquadric/k50/n10/scattered` | 2.68 µs | 2.66 µs | 2.7 µs | 15 µs |  | 5.6x faster |
| `inverse_quadratic/global/n10/scattered` | 2.61 µs | 2.82 µs | 2.6 µs | 12 µs |  | 4.6x faster |
| `inverse_quadratic/k50/n10/scattered` | 2.65 µs | 2.84 µs | 2.85 µs | 16.4 µs |  | 6.2x faster |
| `gaussian/global/n10/scattered` | 3.08 µs | 2.84 µs | 3.05 µs | 12 µs | 8.53 µs | 2.8x faster |
| `gaussian/k50/n10/scattered` | 3.08 µs | 2.8 µs | 3.05 µs | 16.4 µs |  | 5.3x faster |
| `linear/global/n100/scattered` | 67.5 µs | 67.9 µs | 73.6 µs | 76 µs | 704 µs | 1.1x faster |
| `linear/k50/n100/scattered` | 3.21 µs | 3.34 µs | 3.46 µs | 23 µs |  | 7.2x faster |
| `thin_plate_spline/global/n100/scattered` | 96.5 µs | 100 µs | 98.8 µs | 122 µs | 841 µs | 1.3x faster |
| `thin_plate_spline/k50/n100/scattered` | 3.58 µs | 3.35 µs | 3.59 µs | 23.7 µs |  | 6.6x faster |
| `cubic/global/n100/scattered` | 77.9 µs | 78.1 µs | 77.2 µs | 87.9 µs | 741 µs | 1.1x faster |
| `cubic/k50/n100/scattered` | 3.3 µs | 3.5 µs | 3.58 µs | 21.6 µs |  | 6.6x faster |
| `quintic/global/n100/scattered` | 83.7 µs | 84.6 µs | 83.9 µs | 99.7 µs | 809 µs | 1.2x faster |
| `quintic/k50/n100/scattered` | 3.69 µs | 3.67 µs | 3.65 µs | 25.1 µs |  | 6.8x faster |
| `multiquadric/global/n100/scattered` | 81.8 µs | 81.3 µs | 84.2 µs | 83.1 µs | 675 µs | 1x faster |
| `multiquadric/k50/n100/scattered` | 3.49 µs | 3.52 µs | 3.27 µs | 22.9 µs |  | 6.6x faster |
| `inverse_multiquadric/global/n100/scattered` | 85.9 µs | 86.9 µs | 86.6 µs | 85.8 µs | 847 µs | **1x slower** |
| `inverse_multiquadric/k50/n100/scattered` | 3.51 µs | 3.51 µs | 3.5 µs | 21.6 µs |  | 6.1x faster |
| `inverse_quadratic/global/n100/scattered` | 71.4 µs | 79.6 µs | 79 µs | 79.2 µs |  | 1.1x faster |
| `inverse_quadratic/k50/n100/scattered` | 3.56 µs | 3.56 µs | 3.54 µs | 22.6 µs |  | 6.4x faster |
| `gaussian/global/n100/scattered` | 95.5 µs | 107 µs | 115 µs | 106 µs | 839 µs | 1.1x faster |
| `gaussian/k50/n100/scattered` | 3.4 µs | 3.43 µs | 3.49 µs | 20.9 µs |  | 6.2x faster |
| `linear/global/n1000/scattered` | 20 ms | 19.9 ms | 19.7 ms | 58.8 ms | 212 ms | 2.9x faster |
| `linear/k50/n1000/scattered` | 38.6 µs | 35.6 µs | 35 µs | 140 µs |  | 3.6x faster |
| `thin_plate_spline/global/n1000/scattered` | 23.9 ms | 23.4 ms | 23.8 ms | 66.8 ms | 222 ms | 2.8x faster |
| `thin_plate_spline/k50/n1000/scattered` | 38.7 µs | 38.9 µs | 38.4 µs | 142 µs |  | 3.7x faster |
| `cubic/global/n1000/scattered` | 27 ms | 19.7 ms | 20.5 ms | 61.3 ms | 206 ms | 2.3x faster |
| `cubic/k50/n1000/scattered` | 38.8 µs | 39.1 µs | 38.4 µs | 141 µs |  | 3.6x faster |
| `quintic/global/n1000/scattered` | 21.2 ms | 22.2 ms | 22.5 ms | 62.2 ms | 222 ms | 2.9x faster |
| `quintic/k50/n1000/scattered` | 38.9 µs | 38.8 µs | 37.9 µs | 145 µs |  | 3.7x faster |
| `multiquadric/global/n1000/scattered` | 20.7 ms | 21.7 ms | 20.1 ms | 46.2 ms | 209 ms | 2.2x faster |
| `multiquadric/k50/n1000/scattered` | 36.2 µs | 38.9 µs | 38.5 µs | 138 µs |  | 3.8x faster |
| `inverse_multiquadric/global/n1000/scattered` | 21.9 ms | 21.9 ms | 22.2 ms | 59 ms | 225 ms | 2.7x faster |
| `inverse_multiquadric/k50/n1000/scattered` | 35.2 µs | 38.7 µs | 38.4 µs | 147 µs |  | 4.2x faster |
| `inverse_quadratic/global/n1000/scattered` | 20.5 ms | 20.1 ms | 21.1 ms | 59.9 ms |  | 2.9x faster |
| `inverse_quadratic/k50/n1000/scattered` | 38.6 µs | 38.6 µs | 35.2 µs | 138 µs |  | 3.6x faster |
| `gaussian/global/n1000/scattered` | 24.1 ms | 23.3 ms | 24.1 ms | 64.7 ms | 228 ms | 2.7x faster |
| `gaussian/k50/n1000/scattered` | 38.6 µs | 36.6 µs | 34.9 µs | 135 µs |  | 3.5x faster |
| `thin_plate_spline/global/n3000/scattered` | 448 ms | 460 ms | 466 ms | 2.26 s | 6.63 s | 5x faster |
| `thin_plate_spline/k50/n10000/scattered` | 710 µs | 727 µs | 583 µs | 1.7 ms |  | 2.4x faster |
| `thin_plate_spline/k50/n10000/image` | 794 µs | 828 µs | 724 µs | 1.64 ms |  | 2.1x faster |
| `thin_plate_spline/k50/n100000/scattered` | 9.09 ms | 9.56 ms | 9.09 ms | 26.5 ms |  | 2.9x faster |
| `thin_plate_spline/k50/n100000/image` | 9.16 ms | 9.28 ms | 8.88 ms | 26.4 ms |  | 2.9x faster |
| `linear/global/ghcn-tmax-20250101/heldout` | 19 s | 17.4 s | 18.8 s | 128 s | 306 s | 6.8x faster |
| `linear/k50/ghcn-tmax-20250101/heldout` | 754 µs | 860 µs | 720 µs | 1.88 ms |  | 2.5x faster |
| `thin_plate_spline/global/ghcn-tmax-20250101/heldout` | 17.6 s | 17.8 s | 18.2 s | 125 s | 331 s | 7.1x faster |
| `thin_plate_spline/k50/ghcn-tmax-20250101/heldout` | 774 µs | 801 µs | 728 µs | 1.85 ms |  | 2.4x faster |
| `cubic/global/ghcn-tmax-20250101/heldout` | 18.8 s | 17.5 s | 18.8 s | 123 s | 336 s | 6.5x faster |
| `cubic/k50/ghcn-tmax-20250101/heldout` | 715 µs | 1.46 ms | 914 µs | 1.97 ms |  | 2.8x faster |
| `quintic/global/ghcn-tmax-20250101/heldout` | 18.6 s | 39.6 s | 17.2 s | 131 s | 339 s | 7x faster |
| `quintic/k50/ghcn-tmax-20250101/heldout` | 858 µs | 814 µs | 753 µs | 1.88 ms |  | 2.2x faster |

## Memory

Megabytes. "Data read" is the process with the data of the case in memory, before anything is built. The columns after it are what was added to that: the highest while building, with the interpolator built, the highest while evaluating. "Highest" is the highest of the whole run, repeats included.

| case | library | data read | while building | built | while evaluating | highest |
| :--- | :--- | ---: | ---: | ---: | ---: | ---: |
| `oned/cubic_spline/n10000/sorted` | interpx | 4.5 | +0.7 | +0.8 | +0.8 | 6.0 |
| `oned/cubic_spline/n10000/sorted` | interpx, parallel | 3.3 | +2.0 | +2.0 | +2.6 | 6.5 |
| `oned/cubic_spline/n10000/sorted` | SciPy | 71.2 | +1.3 | +1.8 | +1.8 | 73.3 |
| `oned/cubic_spline/n10000/sorted` | scirs2-interpolate | 4.4 | +0.7 | +0.9 | +1.5 | 5.9 |
| `oned/cubic_spline/n10000/sorted` | ndarray-interp | 3.1 | +1.9 | +1.7 | +1.8 | 5.6 |
| `grid/linear/4d/max200/eval10000/ascending` | interpx | 209.2 | +195.6 | +195.6 | +195.6 | 404.9 |
| `grid/linear/4d/max200/eval10000/ascending` | interpx, parallel | 223.3 | +181.8 | +181.5 | +181.7 | 405.4 |
| `grid/linear/4d/max200/eval10000/ascending` | interpx, `from_owned` | 204.6 | +200.1 | +200.1 | +200.1 | 404.9 |
| `grid/linear/4d/max200/eval10000/ascending` | SciPy | 271.3 | +0.3 | +0.4 | +2.2 | 273.6 |
| `grid/linear/4d/max200/eval10000/ascending` | scirs2-interpolate | 214.0 | +185.2 | +191.0 | +191.0 | 404.9 |
| `grid/linear/4d/max200/eval10000/ascending` | interpn | 224.8 | +0.0 | +0.0 | +0.0 | 397.3 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpx | 4.1 | +4.2 | +4.2 | +4.2 | 8.5 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpx, parallel | 5.5 | +2.9 | +2.8 | +2.8 | 8.9 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpx, `from_owned` | 5.6 | +1.8 | +2.7 | +2.7 | 8.5 |
| `grid/cubic/2d/each400/eval1000/ascending` | SciPy | 71.2 | +87.3 | +8.3 | +8.3 | 162.6 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpn | 6.1 | +0.0 | +0.0 | +0.0 | 6.1 |
| `grid/linear/2d/each800/eval10000/ascending` | interpx | 12.2 | +0.0 | +2.7 | +2.7 | 15.0 |
| `grid/linear/2d/each800/eval10000/ascending` | interpx, parallel | 9.8 | +5.2 | +5.2 | +5.2 | 15.5 |
| `grid/linear/2d/each800/eval10000/ascending` | interpx, `from_owned` | 9.5 | +0.4 | +5.5 | +5.5 | 15.0 |
| `grid/linear/2d/each800/eval10000/ascending` | SciPy | 76.5 | +0.2 | +0.2 | +0.8 | 77.3 |
| `grid/linear/2d/each800/eval10000/ascending` | scirs2-interpolate | 10.7 | +0.0 | +4.1 | +4.2 | 15.0 |
| `grid/linear/2d/each800/eval10000/ascending` | ndarray-interp | 9.4 | +0.2 | +0.2 | +0.3 | 9.8 |
| `grid/linear/2d/each800/eval10000/ascending` | interpn | 11.9 | +0.0 | +0.0 | +0.0 | 11.9 |
| `nearest_nd/n100000/image500` | interpx | 14.8 | +0.7 | +0.0 | +1.4 | 18.2 |
| `nearest_nd/n100000/image500` | interpx, parallel | 11.9 | +3.5 | +2.3 | +4.9 | 18.4 |
| `nearest_nd/n100000/image500` | interpx, `from_owned` | 12.5 | +3.2 | +1.7 | +3.7 | 18.2 |
| `nearest_nd/n100000/image500` | SciPy | 77.6 | +4.5 | +4.5 | +11.1 | 92.0 |
| `linear_nd/n100000/image500` | interpx | 12.5 | +24.6 | +15.8 | +17.7 | 41.5 |
| `linear_nd/n100000/image500` | interpx, parallel | 12.6 | +25.6 | +15.6 | +18.1 | 42.4 |
| `linear_nd/n100000/image500` | interpx, `from_owned` | 14.8 | +25.3 | +13.8 | +15.7 | 42.2 |
| `linear_nd/n100000/image500` | SciPy | 75.9 | +68.4 | +15.5 | +27.2 | 145.9 |
| `linear_nd/n100000/image500` | spade | 14.5 | +13.6 | +10.6 | +12.6 | 29.5 |
| `rbf/thin_plate_spline/global/n3000/scattered` | interpx | 3.1 | +77.3 | +5.0 | +5.0 | 80.4 |
| `rbf/thin_plate_spline/global/n3000/scattered` | interpx, parallel | 3.5 | +77.0 | +4.7 | +5.1 | 80.9 |
| `rbf/thin_plate_spline/global/n3000/scattered` | interpx, `from_owned` | 3.3 | +77.2 | +4.8 | +4.8 | 80.5 |
| `rbf/thin_plate_spline/global/n3000/scattered` | SciPy | 71.2 | +73.0 | +0.8 | +9.0 | 144.1 |
| `rbf/thin_plate_spline/global/n3000/scattered` | scirs2-interpolate | 3.3 | +289.2 | +1.2 | +1.2 | 292.5 |
| `rbf/thin_plate_spline/k50/n100000/image` | interpx | 6.9 | +3.5 | +4.5 | +4.8 | 11.9 |
| `rbf/thin_plate_spline/k50/n100000/image` | interpx, parallel | 6.8 | +3.6 | +4.6 | +5.7 | 12.7 |
| `rbf/thin_plate_spline/k50/n100000/image` | interpx, `from_owned` | 6.8 | +5.8 | +4.6 | +4.9 | 13.5 |
| `rbf/thin_plate_spline/k50/n100000/image` | SciPy | 73.7 | +5.4 | +5.6 | +21.9 | 99.0 |

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
| `LinearNDInterpolator` | 4e-15 | 4e-15 | 4e-15 |  |  |  | 6e-15 |
| `RBFInterpolator` | 4e+00 | 4e+00 | 4e+00 | 7e+11 |  |  |  |
