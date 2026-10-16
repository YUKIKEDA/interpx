# Performance: interpx against SciPy and other Rust crates

<!-- Written by compare/perf/report.py from a run of compare/perf/run.py. Do not edit by hand: run it again. -->

Measured on 2026-10-16, at commit `eff8242305`.

## How it was measured

- Machine: Intel64 Family 6 Model 191 Stepping 2, GenuineIntel, 16 logical cores, 48 GB. Windows-11-10.0.26300-SP0.
- rustc 1.99.0 (b940084d7 2026-09-28); Python 3.13.1, NumPy 2.5.4, SciPy 1.18.1.
- The other Rust crates, at the versions in `compare/perf/runners/*/Cargo.toml`: `scirs2-interpolate` (the widest of the Rust crates: one variable, grids, RBF), `ndarray-interp` (the most used for one variable; also bilinear), `interpn` (grids).
- Every library ran each case in a process of its own, pinned to one core, on one thread, on the same data. The libraries took turns case by case, 3 rounds.
- `interpx, from_owned` is interpx built from data that is handed over (`from_owned`) and so not copied. Its column is there for the types that keep their data as it is given (the grid, `NearestNDInterpolator`, `RBFInterpolator`), next to the construction from borrowed data (`new`), which is the one the comparison counts. Each of its samples is given data made just before it, which is not in the nearest cache as the data of the other constructions is: on small grids it reads slower than `new` for that reason.
- `interpx, parallel` is interpx evaluated with `eval_batch_into_par` (the `rayon` feature) on 8 threads, not pinned to a core. Every other column is one thread on one core, and the comparison counts those. A parallel batch times its first points and is shared out only when the rest is worth it, so on small cases the column reads as the serial one does.
- A time is that of one call: construction, or the evaluation of all the queries of the case (`eval_batch_into` for interpx, the call with an array for the others). It is the median of the repeats of a process, and then the median of the rounds.
- A call of two seconds or more is timed once: the construction that made the interpolator, or the evaluation that wrote the values, is its measurement, and it is not made again. A library whose construction or evaluation of a case takes that long is run on that case in the first round only. Its time is then one measurement. (The rule goes by the time, whatever the library.)
- interpx is called slower where its median is more than 1.10 times another library's and its fastest repeat is slower than that library's fastest.
- Memory is the resident set size of the process, sampled every millisecond from outside it.
- Command: `perf/run.py`.

## Where interpx is slower

| family | operation | cases | interpx fastest | slower than SciPy | slower than a Rust crate |
| :--- | :--- | ---: | ---: | ---: | ---: |
| One variable | evaluation | 112 | 112 | 0 | 0 |
| One variable | construction | 112 | 77 | 8 | 32 |
| `PPoly` | evaluation | 4 | 4 | 0 | 0 |
| `PPoly` | construction | 4 | 4 | 0 | 0 |
| `BSpline` | evaluation | 4 | 4 | 0 | 0 |
| `BSpline` | construction | 4 | 4 | 0 | 0 |
| `BSplineCurve` | evaluation | 2 | 2 | 0 | 0 |
| `BSplineCurve` | construction | 2 | 2 | 0 | 0 |
| Regular grid | evaluation | 92 | 74 | 0 | 18 |
| Regular grid | construction | 92 | 27 | 40 | 46 |
| `NearestNDInterpolator` | evaluation | 16 | 16 | 0 | 0 |
| `NearestNDInterpolator` | construction | 16 | 16 | 0 | 0 |
| `RBFInterpolator` | evaluation | 61 | 61 | 0 | 0 |
| `RBFInterpolator` | construction | 61 | 61 | 0 | 0 |

200 comparisons (a case, an operation, and a library), by method. "Same values" says whether that library's values are SciPy's up to rounding: where they are not, it computes another interpolant under the same name, and its time is not that of the same work.

| method | operation | against | same values | cases | slower by | the largest |
| :--- | :--- | :--- | :--- | ---: | ---: | :--- |
| `grid/cubic` | construction | interpn | no | 8 | 1.79e+03x to 9.9e+07x | `grid/cubic/2d/auckland-ba32-0401/heldout/ascending`: 77.3 ms against 0 |
| `grid/linear` | construction | interpn | yes | 35 | 155x to 4.09e+07x | `grid/linear/4d/max200/eval1000/ascending`: 52 ms against 1.27 ns |
| `grid/nearest` | construction | interpn | yes | 3 | 934x to 2.11e+07x | `grid/nearest/2d/auckland-ba32-0401/heldout/ascending`: 16.5 ms against 0 |
| `grid/linear` | construction | SciPy | yes | 37 | 1.17x to 5.79e+03x | `grid/linear/4d/max200/eval1000/flipped`: 127 ms against 21.9 µs |
| `grid/nearest` | construction | SciPy | yes | 2 | 3.04x to 2.14e+03x | `grid/nearest/2d/auckland-ba32-0401/heldout/ascending`: 16.5 ms against 7.69 µs |
| `grid/linear` | construction | ndarray-interp | yes | 15 | 1.55x to 1.78e+03x | `grid/linear/2d/auckland-ba32-0401/heldout/ascending`: 16.1 ms against 9.03 µs |
| `grid/pchip` | construction | SciPy | yes | 1 | 58.9x | `grid/pchip/3d/each40/eval10000/ascending`: 508 µs against 8.63 µs |
| `oned/linear` | construction | ndarray-interp | yes | 14 | 1.38x to 4.75x | `oned/linear/n10000/shuffled`: 92.6 µs against 19.5 µs |
| `oned/pchip` | construction | scirs2-interpolate | yes | 8 | 1.36x to 2.97x | `oned/pchip/n1000/sorted`: 12.8 µs against 4.3 µs |
| `grid/linear` | construction | scirs2-interpolate | yes | 23 | 1.14x to 2.64x | `grid/linear/4d/max40/eval1000/ascending`: 14.6 µs against 5.52 µs |
| `oned/bspline1` | construction | SciPy | yes | 2 | 2.59x to 2.6x | `oned/bspline1/n10000/sorted`: 66.8 µs against 25.7 µs |
| `grid/nearest` | construction | scirs2-interpolate | yes | 3 | 1.29x to 2.56x | `grid/nearest/3d/each40/eval10000/ascending`: 22.8 µs against 8.89 µs |
| `oned/cubic_spline` | construction | scirs2-interpolate | no | 3 | 1.14x to 2.37x | `oned/cubic_spline/n1000/shuffled`: 38.5 µs against 16.2 µs |
| `oned/akima` | construction | scirs2-interpolate | yes | 6 | 1.39x to 2.05x | `oned/akima/n1000/shuffled`: 11.6 µs against 5.63 µs |
| `oned/linear` | construction | SciPy | yes | 4 | 1.85x to 2.02x | `oned/linear/n10000/shuffled`: 92.6 µs against 45.9 µs |
| `oned/linear` | construction | scirs2-interpolate | yes | 14 | 1.18x to 1.69x | `oned/linear/n1000/shuffled`: 2.12 µs against 1.25 µs |
| `grid/linear` | evaluation | interpn | yes | 18 | 1.1x to 1.54x | `grid/linear/3d/max10/eval10/ascending`: 218 ns against 142 ns |
| `oned/nearest` | construction | SciPy | yes | 2 | 1.14x to 1.14x | `oned/nearest/n10000/sorted`: 64 µs against 56 µs |
| `oned/nearest` | construction | scirs2-interpolate | yes | 1 | 1.12x | `oned/nearest/n10000/shuffled`: 67 µs against 59.5 µs |
| `grid/linear` | evaluation | ndarray-interp | yes | 1 | 1.1x | `grid/linear/2d/max40/eval10/ascending`: 172 ns against 156 ns |

<details><summary>Every comparison where interpx is slower</summary>

"Values apart" is the distance of that library's values from SciPy's, as a share of the largest value.

| case | operation | against | interpx | the other | slower by | values apart |
| :--- | :--- | :--- | ---: | ---: | ---: | ---: |
| `grid/cubic/2d/auckland-ba32-0401/heldout/ascending` | construction | interpn | 77.3 ms | 0 | 98994688.00x | 5e-03 |
| `grid/linear/4d/max200/eval1000/ascending` | construction | interpn | 52 ms | 1.27 ns | 40942828.31x | 3e-16 |
| `grid/linear/4d/max200/eval10/ascending` | construction | interpn | 52.5 ms | 1.37 ns | 38401170.29x | 2e-16 |
| `grid/linear/4d/max200/eval100/ascending` | construction | interpn | 52.2 ms | 1.37 ns | 38190372.57x | 3e-16 |
| `grid/linear/4d/max200/eval10000/ascending` | construction | interpn | 49.3 ms | 1.37 ns | 36094537.14x | 4e-16 |
| `grid/linear/2d/auckland-ba32-0401/heldout/ascending` | construction | interpn | 16.1 ms | 0 | 23486025.14x | 0e+00 |
| `grid/nearest/2d/auckland-ba32-0401/heldout/ascending` | construction | interpn | 16.5 ms | 0 | 21102592.00x | 0e+00 |
| `grid/linear/3d/max200/eval10000/ascending` | construction | interpn | 1.82 ms | 0 | 2071779.56x | 5e-16 |
| `grid/linear/3d/max200/eval10/ascending` | construction | interpn | 1.79 ms | 0 | 2036963.56x | 3e-16 |
| `grid/linear/3d/max200/eval1000/ascending` | construction | interpn | 1.7 ms | 0 | 1938090.67x | 3e-16 |
| `grid/linear/3d/max200/eval100/ascending` | construction | interpn | 1.71 ms | 0 | 1750630.40x | 3e-16 |
| `grid/cubic/2d/each400/eval1000/ascending` | construction | interpn | 1.4 ms | 0 | 1596871.11x | 5e-01 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | interpn | 1.09 ms | 0 | 1596854.86x | 3e-16 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | interpn | 319 µs | 0 | 408448.00x | 3e-16 |
| `grid/cubic/3d/each40/eval10000/ascending` | construction | interpn | 324 µs | 1.17 ns | 276565.33x | 5e-01 |
| `grid/cubic/2d/each200/eval1000/ascending` | construction | interpn | 186 µs | 0 | 211512.89x | 2e-01 |
| `grid/cubic/2d/each100/eval1000/ascending` | construction | interpn | 81.7 µs | 0 | 83660.80x | 2e-01 |
| `grid/linear/2d/max200/eval10000/ascending` | construction | interpn | 42.2 µs | 0 | 61696.00x | 3e-16 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | interpn | 22.8 µs | 0 | 23334.40x | 0e+00 |
| `grid/linear/3d/max40/eval1000/ascending` | construction | interpn | 17 µs | 0 | 19345.78x | 2e-16 |
| `grid/linear/3d/max40/eval10000/ascending` | construction | interpn | 14.6 µs | 0 | 16618.67x | 3e-16 |
| `grid/linear/4d/max40/eval1000/ascending` | construction | interpn | 14.6 µs | 1.27 ns | 11475.69x | 6e-16 |
| `grid/linear/4d/max40/eval100/ascending` | construction | interpn | 14.5 µs | 1.37 ns | 10573.71x | 2e-16 |
| `grid/linear/4d/max40/eval10/ascending` | construction | interpn | 14.4 µs | 1.37 ns | 10528.00x | 2e-16 |
| `grid/linear/4d/max40/eval10000/ascending` | construction | interpn | 13.1 µs | 1.37 ns | 9581.71x | 5e-16 |
| `grid/cubic/2d/each40/eval10000/ascending` | construction | interpn | 8.38 µs | 0 | 9534.22x | 7e-01 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | interpn | 7.44 µs | 0 | 9528.00x | 2e-16 |
| `grid/linear/2d/max200/eval10/ascending` | construction | interpn | 7.42 µs | 0 | 9492.00x | 2e-16 |
| `grid/cubic/2d/each40/eval1000/ascending` | construction | interpn | 9.06 µs | 0 | 9273.60x | 7e-01 |
| `grid/linear/2d/max200/eval100/ascending` | construction | interpn | 6.81 µs | 0 | 8720.00x | 1e-16 |
| `grid/linear/3d/max40/eval100/ascending` | construction | interpn | 5.25 µs | 0 | 5978.67x | 2e-16 |
| `grid/linear/4d/max200/eval1000/flipped` | construction | SciPy | 127 ms | 21.9 µs | 5791.47x |  |
| `grid/linear/4d/max200/eval10/ascending` | construction | SciPy | 52.5 ms | 9.56 µs | 5493.06x |  |
| `grid/linear/4d/max200/eval10000/ascending` | construction | SciPy | 49.3 ms | 9.31 µs | 5302.67x |  |
| `grid/linear/4d/max200/eval100/flipped` | construction | SciPy | 108 ms | 20.4 µs | 5259.62x |  |
| `grid/linear/4d/max200/eval1000/ascending` | construction | SciPy | 52 ms | 10.1 µs | 5127.32x |  |
| `grid/linear/4d/max200/eval100/ascending` | construction | SciPy | 52.2 ms | 10.4 µs | 5020.52x |  |
| `grid/linear/4d/max200/eval10000/flipped` | construction | SciPy | 102 ms | 20.4 µs | 5014.69x |  |
| `grid/linear/4d/max200/eval10/flipped` | construction | SciPy | 104 ms | 22 µs | 4732.56x |  |
| `grid/linear/3d/max40/eval10/ascending` | construction | interpn | 2.97 µs | 0 | 3383.11x | 2e-16 |
| `grid/linear/2d/auckland-ba32-0401/heldout/ascending` | construction | SciPy | 16.1 ms | 7.04 µs | 2281.34x |  |
| `grid/nearest/2d/auckland-ba32-0401/heldout/ascending` | construction | SciPy | 16.5 ms | 7.69 µs | 2144.14x |  |
| `grid/cubic/2d/each10/eval1000/ascending` | construction | interpn | 1.58 µs | 0 | 1794.67x | 3e-01 |
| `grid/linear/2d/auckland-ba32-0401/heldout/ascending` | construction | ndarray-interp | 16.1 ms | 9.03 µs | 1778.94x | 0e+00 |
| `grid/nearest/2d/each40/eval10000/ascending` | construction | interpn | 638 ns | 0 | 933.71x | 0e+00 |
| `grid/linear/2d/max40/eval10/ascending` | construction | interpn | 362 ns | 0 | 528.86x | 2e-16 |
| `grid/linear/2d/max40/eval1000/ascending` | construction | interpn | 402 ns | 0 | 514.25x | 3e-16 |
| `grid/linear/2d/max40/eval100/ascending` | construction | interpn | 394 ns | 0 | 503.75x | 2e-16 |
| `grid/linear/2d/max40/eval10000/ascending` | construction | interpn | 392 ns | 0 | 501.38x | 4e-16 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | ndarray-interp | 1.09 ms | 2.44 µs | 446.69x | 6e-16 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | ndarray-interp | 319 µs | 1.17 µs | 271.85x | 5e-16 |
| `grid/linear/3d/max200/eval10000/ascending` | construction | SciPy | 1.82 ms | 7.51 µs | 242.38x |  |
| `grid/linear/3d/max200/eval100/ascending` | construction | SciPy | 1.71 ms | 7.64 µs | 223.84x |  |
| `grid/linear/3d/max200/eval1000/ascending` | construction | SciPy | 1.7 ms | 7.62 µs | 223.67x |  |
| `grid/linear/3d/max10/eval10/ascending` | construction | interpn | 196 ns | 0 | 223.11x | 2e-16 |
| `grid/linear/3d/max200/eval10/ascending` | construction | SciPy | 1.79 ms | 8.23 µs | 217.58x |  |
| `grid/linear/3d/max10/eval100/ascending` | construction | interpn | 189 ns | 0 | 215.33x | 3e-16 |
| `grid/linear/3d/max10/eval10000/ascending` | construction | interpn | 180 ns | 0 | 204.67x | 4e-16 |
| `grid/linear/3d/max200/eval10/flipped` | construction | SciPy | 3.21 ms | 16.1 µs | 199.46x |  |
| `grid/linear/3d/max200/eval10000/flipped` | construction | SciPy | 3.2 ms | 16 µs | 199.35x |  |
| `grid/linear/3d/max200/eval1000/flipped` | construction | SciPy | 3.23 ms | 16.2 µs | 199.26x |  |
| `grid/linear/3d/max10/eval1000/ascending` | construction | interpn | 193 ns | 0 | 197.80x | 5e-16 |
| `grid/linear/2d/max10/eval10/ascending` | construction | interpn | 121 ns | 0 | 177.71x | 1e-16 |
| `grid/linear/2d/max10/eval1000/ascending` | construction | interpn | 121 ns | 0 | 177.29x | 3e-16 |
| `grid/linear/3d/max200/eval100/flipped` | construction | SciPy | 3.08 ms | 17.6 µs | 175.20x |  |
| `grid/linear/2d/each800/eval10000/ascending` | construction | SciPy | 1.09 ms | 6.38 µs | 171.23x |  |
| `grid/linear/2d/max10/eval10000/ascending` | construction | interpn | 132 ns | 0 | 168.62x | 3e-16 |
| `grid/linear/2d/max10/eval100/ascending` | construction | interpn | 121 ns | 0 | 155.38x | 2e-16 |
| `grid/linear/2d/max200/eval10000/ascending` | construction | ndarray-interp | 42.2 µs | 517 ns | 81.61x | 5e-16 |
| `grid/pchip/3d/each40/eval10000/ascending` | construction | SciPy | 508 µs | 8.63 µs | 58.87x |  |
| `grid/linear/2d/each400/eval10000/ascending` | construction | SciPy | 319 µs | 5.97 µs | 53.49x |  |
| `grid/linear/2d/max200/eval10/ascending` | construction | ndarray-interp | 7.42 µs | 472 ns | 15.72x | 2e-16 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | ndarray-interp | 7.44 µs | 517 ns | 14.39x | 4e-16 |
| `grid/linear/2d/max200/eval100/ascending` | construction | ndarray-interp | 6.81 µs | 517 ns | 13.19x | 3e-16 |
| `grid/linear/2d/max200/eval10000/ascending` | construction | SciPy | 42.2 µs | 5.96 µs | 7.07x |  |
| `oned/linear/n10000/shuffled` | construction | ndarray-interp | 92.6 µs | 19.5 µs | 4.75x | 0e+00 |
| `oned/linear/n10000/sorted` | construction | ndarray-interp | 90.5 µs | 19.9 µs | 4.55x | 0e+00 |
| `oned/linear/n10000/shuffled-uneven` | construction | ndarray-interp | 89.9 µs | 20 µs | 4.49x | 0e+00 |
| `grid/linear/2d/max200/eval10000/flipped` | construction | SciPy | 57 µs | 13 µs | 4.40x |  |
| `oned/linear/n10000/sorted-uneven` | construction | ndarray-interp | 92.1 µs | 21.4 µs | 4.29x | 0e+00 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | SciPy | 22.8 µs | 7.5 µs | 3.04x |  |
| `oned/pchip/n1000/sorted` | construction | scirs2-interpolate | 12.8 µs | 4.3 µs | 2.97x | 2e-16 |
| `grid/linear/4d/max40/eval10000/flipped` | construction | SciPy | 56.8 µs | 19.8 µs | 2.87x |  |
| `oned/pchip/n1000/shuffled` | construction | scirs2-interpolate | 11.5 µs | 4.16 µs | 2.77x | 2e-16 |
| `grid/linear/4d/max40/eval100/flipped` | construction | SciPy | 56.5 µs | 20.6 µs | 2.74x |  |
| `grid/linear/4d/max40/eval1000/ascending` | construction | scirs2-interpolate | 14.6 µs | 5.52 µs | 2.64x | 2e-16 |
| `grid/linear/4d/max40/eval10/flipped` | construction | SciPy | 56.4 µs | 21.6 µs | 2.61x |  |
| `oned/bspline1/n10000/sorted` | construction | SciPy | 66.8 µs | 25.7 µs | 2.60x |  |
| `oned/bspline1/n10000/shuffled` | construction | SciPy | 60.8 µs | 23.4 µs | 2.59x |  |
| `grid/linear/4d/max40/eval1000/flipped` | construction | SciPy | 56.4 µs | 21.8 µs | 2.59x |  |
| `grid/linear/2d/max40/eval10000/ascending` | construction | ndarray-interp | 392 ns | 152 ns | 2.57x | 4e-16 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | scirs2-interpolate | 22.8 µs | 8.89 µs | 2.56x | 0e+00 |
| `grid/linear/2d/max200/eval10/ascending` | construction | scirs2-interpolate | 7.42 µs | 2.93 µs | 2.53x | 1e-16 |
| `oned/linear/n10/shuffled` | construction | ndarray-interp | 101 ns | 40.7 ns | 2.48x | 0e+00 |
| `grid/linear/4d/max40/eval10/ascending` | construction | scirs2-interpolate | 14.4 µs | 5.98 µs | 2.41x | 2e-16 |
| `grid/linear/2d/max40/eval1000/ascending` | construction | ndarray-interp | 402 ns | 167 ns | 2.40x | 4e-16 |
| `grid/linear/4d/max40/eval100/ascending` | construction | scirs2-interpolate | 14.5 µs | 6.08 µs | 2.38x | 3e-16 |
| `oned/cubic_spline/n1000/shuffled` | construction | scirs2-interpolate | 38.5 µs | 16.2 µs | 2.37x | 4e-04 |
| `grid/linear/2d/max40/eval100/ascending` | construction | ndarray-interp | 394 ns | 169 ns | 2.33x | 2e-16 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | scirs2-interpolate | 7.44 µs | 3.26 µs | 2.28x | 2e-16 |
| `oned/cubic_spline/n1000/sorted` | construction | scirs2-interpolate | 38.4 µs | 16.9 µs | 2.28x | 4e-04 |
| `grid/linear/4d/max40/eval10000/ascending` | construction | scirs2-interpolate | 13.1 µs | 5.99 µs | 2.19x | 4e-16 |
| `grid/linear/2d/max40/eval10/ascending` | construction | ndarray-interp | 362 ns | 167 ns | 2.16x | 2e-16 |
| `grid/nearest/2d/each40/eval10000/ascending` | construction | scirs2-interpolate | 638 ns | 296 ns | 2.16x | 0e+00 |
| `grid/linear/2d/max200/eval100/ascending` | construction | scirs2-interpolate | 6.81 µs | 3.22 µs | 2.12x | 2e-16 |
| `oned/linear/n10/sorted` | construction | ndarray-interp | 95.1 ns | 45.4 ns | 2.09x | 0e+00 |
| `oned/akima/n1000/shuffled` | construction | scirs2-interpolate | 11.6 µs | 5.63 µs | 2.05x | 3e-17 |
| `oned/linear/n10000/shuffled` | construction | SciPy | 92.6 µs | 45.9 µs | 2.02x |  |
| `grid/linear/3d/max40/eval1000/ascending` | construction | scirs2-interpolate | 17 µs | 8.46 µs | 2.01x | 2e-16 |
| `oned/akima/n1000/sorted` | construction | scirs2-interpolate | 12.6 µs | 6.49 µs | 1.95x | 3e-17 |
| `oned/linear/n10000/sorted-uneven` | construction | SciPy | 92.1 µs | 47.9 µs | 1.93x |  |
| `oned/linear/n10000/sorted` | construction | SciPy | 90.5 µs | 47.8 µs | 1.89x |  |
| `grid/linear/2d/max40/eval1000/ascending` | construction | scirs2-interpolate | 402 ns | 215 ns | 1.87x | 3e-16 |
| `oned/linear/n10000/shuffled-uneven` | construction | SciPy | 89.9 µs | 48.5 µs | 1.85x |  |
| `oned/akima/n100/shuffled` | construction | scirs2-interpolate | 959 ns | 523 ns | 1.83x | 1e-16 |
| `oned/pchip/n100/sorted` | construction | scirs2-interpolate | 968 ns | 538 ns | 1.80x | 2e-16 |
| `oned/pchip/n100/shuffled-uneven` | construction | scirs2-interpolate | 968 ns | 543 ns | 1.78x | 2e-16 |
| `grid/linear/3d/max40/eval1000/flipped` | construction | SciPy | 30.2 µs | 17.3 µs | 1.74x |  |
| `oned/linear/n1000/shuffled` | construction | scirs2-interpolate | 2.12 µs | 1.25 µs | 1.69x | 0e+00 |
| `grid/linear/2d/max10/eval10000/ascending` | construction | ndarray-interp | 132 ns | 78.4 ns | 1.68x | 4e-16 |
| `grid/linear/2d/max40/eval100/ascending` | construction | scirs2-interpolate | 394 ns | 235 ns | 1.68x | 3e-16 |
| `grid/linear/2d/max40/eval10000/ascending` | construction | scirs2-interpolate | 392 ns | 234 ns | 1.67x | 2e-16 |
| `grid/linear/3d/max40/eval10000/flipped` | construction | SciPy | 26.6 µs | 15.9 µs | 1.67x |  |
| `oned/linear/n10000/shuffled` | construction | scirs2-interpolate | 92.6 µs | 56.5 µs | 1.64x | 0e+00 |
| `oned/pchip/n100/shuffled` | construction | scirs2-interpolate | 882 ns | 539 ns | 1.64x | 2e-16 |
| `oned/pchip/n50/shuffled` | construction | scirs2-interpolate | 546 ns | 334 ns | 1.64x | 1e-16 |
| `oned/linear/n1000/sorted` | construction | scirs2-interpolate | 1.94 µs | 1.22 µs | 1.59x | 0e+00 |
| `oned/akima/n100/sorted` | construction | scirs2-interpolate | 890 ns | 562 ns | 1.58x | 1e-16 |
| `oned/linear/n50/shuffled` | construction | ndarray-interp | 156 ns | 98.9 ns | 1.58x | 0e+00 |
| `grid/linear/2d/max10/eval100/ascending` | construction | ndarray-interp | 121 ns | 77.9 ns | 1.56x | 4e-16 |
| `grid/linear/2d/max10/eval1000/ascending` | construction | ndarray-interp | 121 ns | 77.9 ns | 1.56x | 3e-16 |
| `grid/linear/2d/max10/eval10/ascending` | construction | ndarray-interp | 121 ns | 78.4 ns | 1.55x | 1e-16 |
| `oned/akima/n50/shuffled` | construction | scirs2-interpolate | 550 ns | 355 ns | 1.55x | 1e-16 |
| `oned/linear/n50/sorted` | construction | ndarray-interp | 167 ns | 108 ns | 1.55x | 0e+00 |
| `oned/linear/n100/shuffled-uneven` | construction | ndarray-interp | 256 ns | 166 ns | 1.54x | 0e+00 |
| `oned/linear/n10000/sorted-uneven` | construction | scirs2-interpolate | 92.1 µs | 59.8 µs | 1.54x | 6e-17 |
| `oned/linear/n10000/shuffled-uneven` | construction | scirs2-interpolate | 89.9 µs | 58.3 µs | 1.54x | 6e-17 |
| `oned/linear/n10000/sorted` | construction | scirs2-interpolate | 90.5 µs | 58.7 µs | 1.54x | 0e+00 |
| `grid/linear/3d/max10/eval10/ascending` | evaluation | interpn | 218 ns | 142 ns | 1.54x | 2e-16 |
| `oned/linear/n1000/shuffled` | construction | ndarray-interp | 2.12 µs | 1.38 µs | 1.54x | 0e+00 |
| `grid/linear/4d/max40/eval10/ascending` | construction | SciPy | 14.4 µs | 9.42 µs | 1.53x |  |
| `oned/linear/n100/sorted` | construction | ndarray-interp | 252 ns | 166 ns | 1.52x | 0e+00 |
| `grid/linear/2d/max40/eval10/ascending` | construction | scirs2-interpolate | 362 ns | 240 ns | 1.50x | 4e-17 |
| `oned/linear/n100/sorted-uneven` | construction | ndarray-interp | 248 ns | 166 ns | 1.49x | 0e+00 |
| `oned/pchip/n100/sorted-uneven` | construction | scirs2-interpolate | 882 ns | 593 ns | 1.49x | 2e-16 |
| `oned/linear/n100/shuffled` | construction | ndarray-interp | 270 ns | 182 ns | 1.49x | 0e+00 |
| `grid/linear/2d/max10/eval10000/ascending` | evaluation | interpn | 137 µs | 92.4 µs | 1.48x | 3e-16 |
| `grid/linear/3d/max200/eval10000/ascending` | construction | scirs2-interpolate | 1.82 ms | 1.27 ms | 1.43x | 3e-16 |
| `grid/linear/4d/max40/eval1000/ascending` | construction | SciPy | 14.6 µs | 10.2 µs | 1.43x |  |
| `grid/linear/2d/max10/eval10/ascending` | evaluation | interpn | 134 ns | 94.8 ns | 1.42x | 1e-16 |
| `grid/linear/4d/max40/eval100/ascending` | construction | SciPy | 14.5 µs | 10.2 µs | 1.41x |  |
| `grid/linear/4d/max40/eval10000/ascending` | construction | SciPy | 13.1 µs | 9.39 µs | 1.40x |  |
| `oned/akima/n50/sorted` | construction | scirs2-interpolate | 503 ns | 361 ns | 1.39x | 1e-16 |
| `grid/linear/3d/max200/eval10/ascending` | construction | scirs2-interpolate | 1.79 ms | 1.29 ms | 1.38x | 3e-16 |
| `oned/linear/n1000/sorted` | construction | ndarray-interp | 1.94 µs | 1.41 µs | 1.38x | 0e+00 |
| `grid/linear/3d/max10/eval1000/ascending` | evaluation | interpn | 20.9 µs | 15.3 µs | 1.37x | 5e-16 |
| `grid/linear/3d/max10/eval10000/ascending` | evaluation | interpn | 193 µs | 142 µs | 1.36x | 4e-16 |
| `oned/pchip/n50/sorted` | construction | scirs2-interpolate | 496 ns | 365 ns | 1.36x | 1e-16 |
| `grid/linear/2d/max10/eval1000/ascending` | evaluation | interpn | 12.5 µs | 9.23 µs | 1.35x | 3e-16 |
| `grid/linear/3d/max200/eval100/ascending` | construction | scirs2-interpolate | 1.71 ms | 1.27 ms | 1.35x | 3e-16 |
| `grid/linear/2d/max10/eval100/ascending` | evaluation | interpn | 1.25 µs | 928 ns | 1.35x | 2e-16 |
| `grid/linear/2d/max200/eval1000/flipped` | construction | SciPy | 15.7 µs | 11.8 µs | 1.33x |  |
| `oned/linear/n100/shuffled-uneven` | construction | scirs2-interpolate | 256 ns | 193 ns | 1.33x | 3e-17 |
| `oned/linear/n50/sorted` | construction | scirs2-interpolate | 167 ns | 127 ns | 1.32x | 0e+00 |
| `grid/linear/4d/max200/eval10000/ascending` | construction | scirs2-interpolate | 49.3 ms | 37.7 ms | 1.31x | 4e-16 |
| `grid/linear/4d/max200/eval100/ascending` | construction | scirs2-interpolate | 52.2 ms | 40.1 ms | 1.30x | 1e-16 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | scirs2-interpolate | 1.09 ms | 838 µs | 1.30x | 2e-16 |
| `oned/linear/n100/shuffled` | construction | scirs2-interpolate | 270 ns | 209 ns | 1.29x | 0e+00 |
| `grid/nearest/2d/auckland-ba32-0401/heldout/ascending` | construction | scirs2-interpolate | 16.5 ms | 12.8 ms | 1.29x | 0e+00 |
| `oned/linear/n50/shuffled` | construction | scirs2-interpolate | 156 ns | 121 ns | 1.29x | 0e+00 |
| `grid/linear/2d/max200/eval100/flipped` | construction | SciPy | 15.9 µs | 12.3 µs | 1.29x |  |
| `grid/linear/4d/max200/eval1000/ascending` | construction | scirs2-interpolate | 52 ms | 40.6 ms | 1.28x | 3e-16 |
| `grid/linear/3d/max10/eval100/ascending` | evaluation | interpn | 1.95 µs | 1.53 µs | 1.28x | 3e-16 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | SciPy | 7.44 µs | 5.85 µs | 1.27x |  |
| `grid/linear/4d/max200/eval10/ascending` | construction | scirs2-interpolate | 52.5 ms | 41.5 ms | 1.27x | 2e-16 |
| `grid/linear/3d/max200/eval1000/ascending` | construction | scirs2-interpolate | 1.7 ms | 1.35 ms | 1.26x | 3e-16 |
| `grid/linear/2d/max200/eval10/flipped` | construction | SciPy | 15.9 µs | 12.8 µs | 1.24x |  |
| `grid/linear/2d/max200/eval10/ascending` | construction | SciPy | 7.42 µs | 6.06 µs | 1.22x |  |
| `oned/linear/n100/sorted` | construction | scirs2-interpolate | 252 ns | 208 ns | 1.21x | 0e+00 |
| `oned/linear/n100/sorted-uneven` | construction | scirs2-interpolate | 248 ns | 206 ns | 1.20x | 3e-17 |
| `oned/linear/n10/sorted` | construction | scirs2-interpolate | 95.1 ns | 80 ns | 1.19x | 0e+00 |
| `oned/linear/n10/shuffled` | construction | scirs2-interpolate | 101 ns | 85.6 ns | 1.18x | 0e+00 |
| `grid/linear/2d/auckland-ba32-0401/heldout/ascending` | construction | scirs2-interpolate | 16.1 ms | 13.6 ms | 1.18x | 0e+00 |
| `grid/linear/3d/max40/eval1000/ascending` | evaluation | interpn | 23.5 µs | 20 µs | 1.17x | 2e-16 |
| `grid/linear/2d/max200/eval100/ascending` | construction | SciPy | 6.81 µs | 5.82 µs | 1.17x |  |
| `grid/linear/2d/max40/eval10/ascending` | evaluation | interpn | 172 ns | 147 ns | 1.17x | 2e-16 |
| `grid/linear/2d/max40/eval1000/ascending` | evaluation | interpn | 15.1 µs | 13.2 µs | 1.15x | 3e-16 |
| `grid/linear/2d/max40/eval100/ascending` | evaluation | interpn | 1.66 µs | 1.45 µs | 1.15x | 2e-16 |
| `oned/nearest/n10000/sorted` | construction | SciPy | 64 µs | 56 µs | 1.14x |  |
| `grid/linear/3d/max200/eval1000/ascending` | evaluation | interpn | 29.8 µs | 26.1 µs | 1.14x | 3e-16 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | scirs2-interpolate | 319 µs | 280 µs | 1.14x | 2e-16 |
| `oned/cubic_spline/n100/sorted-uneven` | construction | scirs2-interpolate | 1.73 µs | 1.52 µs | 1.14x | 6e-03 |
| `oned/nearest/n10000/shuffled` | construction | SciPy | 67 µs | 58.9 µs | 1.14x |  |
| `grid/linear/2d/max40/eval10000/ascending` | evaluation | interpn | 167 µs | 147 µs | 1.13x | 4e-16 |
| `grid/linear/3d/max200/eval10000/ascending` | evaluation | interpn | 298 µs | 263 µs | 1.13x | 5e-16 |
| `oned/nearest/n10000/shuffled` | construction | scirs2-interpolate | 67 µs | 59.5 µs | 1.12x | 0e+00 |
| `grid/linear/3d/max40/eval10/ascending` | evaluation | interpn | 222 ns | 199 ns | 1.12x | 2e-16 |
| `grid/linear/2d/max200/eval1000/ascending` | evaluation | interpn | 20.7 µs | 18.7 µs | 1.10x | 2e-16 |
| `grid/linear/4d/max40/eval10/ascending` | evaluation | interpn | 328 ns | 298 ns | 1.10x | 2e-16 |
| `grid/linear/2d/max40/eval10/ascending` | evaluation | ndarray-interp | 172 ns | 156 ns | 1.10x | 2e-16 |

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

A run (one library on one case, a process of its own) is stopped at 15 minutes. The longest run that finished took 300 s, 33% of that: scirs2-interpolate on `rbf/cubic/global/ghcn-tmax-20250101/heldout`.

No run took more than half the limit (450 s).

## Real-world data

Two data sets, each taken whole in the unit its source gives it out in; what was taken and how it was divided was fixed before anything was measured (`compare/perf/real.py`, docs/design.ja.md section 8.4.7).

- `auckland-ba32-0401`: Auckland Part 1 LiDAR 1m DEM (2024), tile BA32_10000_0401. 4800 by 7200 cells of 1 m. The grid is every second row and column (and the last ones); every other cell is a query, and its elevation is the value held out. Licensed by Regional Software Holdings Limited for re-use under CC-BY-4.0. Hosted by Toitū Te Whenua Land Information New Zealand (https://registry.opendata.aws/nz-elevation).
- `ghcn-tmax-20250101`: GHCN-Daily, maximum temperature of 2025-01-01, every station; longitude and latitude in degrees, as they are. Of every 10 stations in the order of the file the last is a query, and its temperature is the value held out. The RBF kernels are those without a shape parameter. NOAA Global Historical Climatology Network Daily (GHCN-D), from https://registry.opendata.aws/noaa-ghcn (CC0 1.0).

These cases are also in the tables of their families below and in the counts above. Here they are with the error: the largest and the root mean square distance from the values held out (metres, and degrees Celsius), and the number of queries a library gave no value. A library missing from a case does not have it, or did not finish (see the runs that failed).

`grid/nearest/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 16.5 ms | 364 ms | 15.5 | 0.178 | 0 |
| interpx, parallel | 15.4 ms | 61.8 ms | 15.5 | 0.178 | 0 |
| interpx, `from_owned` | 5.52 ms | 348 ms | 15.5 | 0.178 | 0 |
| SciPy | 7.69 µs | 1.34 s | 15.5 | 0.178 | 0 |
| scirs2-interpolate | 12.8 ms | 158 s | 15.5 | 0.178 | 0 |
| interpn | 0 | 750 ms | 15.5 | 0.178 | 0 |

`grid/linear/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 16.1 ms | 377 ms | 6.74 | 0.0601 | 0 |
| interpx, parallel | 15.5 ms | 65 ms | 6.74 | 0.0601 | 0 |
| interpx, `from_owned` | 5.45 ms | 370 ms | 6.74 | 0.0601 | 0 |
| SciPy | 7.04 µs | 1.09 s | 6.74 | 0.0601 | 0 |
| scirs2-interpolate | 13.6 ms | 73.2 s | 6.74 | 0.0601 | 0 |
| ndarray-interp | 9.03 µs | 369 ms | 6.74 | 0.0601 | 0 |
| interpn | 0 | 772 ms | 6.74 | 0.0601 | 0 |

`grid/cubic/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 77.3 ms | 723 ms | 6.77 | 0.0507 | 0 |
| interpx, parallel | 73.7 ms | 124 ms | 6.77 | 0.0507 | 0 |
| interpx, `from_owned` | 66.6 ms | 736 ms | 6.77 | 0.0507 | 0 |
| SciPy | 14.3 s | 31.8 s | 6.77 | 0.0507 | 0 |
| interpn | 0 | 1.62 s | 6.76 | 0.0517 | 0 |

`grid/quintic/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 90.9 ms | 1.19 s | 6.77 | 0.0517 | 0 |
| interpx, parallel | 95.8 ms | 207 ms | 6.77 | 0.0517 | 0 |
| interpx, `from_owned` | 83.6 ms | 1.2 s | 6.77 | 0.0517 | 0 |
| SciPy | 20.3 s | 35.3 s | 6.77 | 0.0517 | 0 |

`grid/pchip/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 72.1 ms | 1.02 s | 6.75 | 0.049 | 0 |
| interpx, parallel | 72.7 ms | 172 ms | 6.75 | 0.049 | 0 |
| interpx, `from_owned` | 62 ms | 1.03 s | 6.75 | 0.049 | 0 |

`nearest_nd/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 710 µs | 159 µs | 25 | 3.42 | 0 |
| interpx, parallel | 670 µs | 66.5 µs | 25 | 3.42 | 0 |
| interpx, `from_owned` | 721 µs | 171 µs | 25 | 3.42 | 0 |
| SciPy | 1.66 ms | 513 µs | 25 | 3.42 | 0 |

`rbf/linear/global/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 17.6 s | 19 ms | 20.8 | 2.74 | 0 |
| interpx, parallel | 16.6 s | 3.34 ms | 20.8 | 2.74 | 0 |
| interpx, `from_owned` | 17 s | 20.3 ms | 20.8 | 2.74 | 0 |
| SciPy | 121 s | 53.2 ms | 20.8 | 2.74 | 0 |
| scirs2-interpolate | 268 s | 229 ms | 20.8 | 2.74 | 0 |

`rbf/linear/k50/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 821 µs | 17.7 ms | 20.8 | 2.75 | 0 |
| interpx, parallel | 781 µs | 3.35 ms | 20.8 | 2.75 | 0 |
| interpx, `from_owned` | 703 µs | 17.8 ms | 20.8 | 2.75 | 0 |
| SciPy | 1.69 ms | 44.3 ms | 20.8 | 2.75 | 0 |

`rbf/thin_plate_spline/global/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 17 s | 84 ms | 102 | 5.67 | 0 |
| interpx, parallel | 45 s | 79.4 ms | 102 | 5.67 | 0 |
| interpx, `from_owned` | 18.4 s | 85.7 ms | 102 | 5.67 | 0 |
| SciPy | 109 s | 112 ms | 102 | 5.67 | 0 |
| scirs2-interpolate | 282 s | 413 ms | 93.7 | 5.41 | 0 |

`rbf/thin_plate_spline/k50/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 813 µs | 28.8 ms | 102 | 5.66 | 0 |
| interpx, parallel | 808 µs | 5.14 ms | 102 | 5.66 | 0 |
| interpx, `from_owned` | 709 µs | 28.6 ms | 102 | 5.66 | 0 |
| SciPy | 1.73 ms | 54.7 ms | 102 | 5.66 | 0 |

`rbf/cubic/global/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 17.5 s | 19.6 ms | 587 | 27.1 | 0 |
| interpx, parallel | 17.1 s | 3.4 ms | 587 | 27.1 | 0 |
| interpx, `from_owned` | 17.3 s | 17.8 ms | 587 | 27.1 | 0 |
| SciPy | 113 s | 51.2 ms | 597 | 24.4 | 0 |
| scirs2-interpolate | 298 s | 264 ms | 100 | 8.48 | 0 |

`rbf/cubic/k50/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 772 µs | 21.2 ms | 567 | 25.3 | 0 |
| interpx, parallel | 839 µs | 3.61 ms | 567 | 25.3 | 0 |
| interpx, `from_owned` | 717 µs | 19.8 ms | 567 | 25.3 | 0 |
| SciPy | 1.73 ms | 49.1 ms | 567 | 25.3 | 0 |

`rbf/quintic/global/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 18 s | 19.9 ms | 2.55e+05 | 1.36e+04 | 0 |
| interpx, parallel | 17.2 s | 3.99 ms | 2.55e+05 | 1.36e+04 | 0 |
| interpx, `from_owned` | 18.1 s | 22.2 ms | 2.55e+05 | 1.36e+04 | 0 |
| SciPy | 115 s | 58.5 ms | 9.39e+04 | 4.92e+03 | 0 |
| scirs2-interpolate | 282 s | 295 ms | 6.56e+16 | 2.85e+15 | 0 |

`rbf/quintic/k50/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 824 µs | 22.9 ms | 3.74e+03 | 199 | 0 |
| interpx, parallel | 846 µs | 3.91 ms | 3.74e+03 | 199 | 0 |
| interpx, `from_owned` | 783 µs | 23.2 ms | 3.74e+03 | 199 | 0 |
| SciPy | 1.75 ms | 51.7 ms | 3.74e+03 | 200 | 0 |

## One variable

Data of `n` points, evaluated at `4 n` points: in order (`sorted`) and in no order (`shuffled`).

### Evaluation

| case | interpx | interpx, parallel | SciPy | scirs2-interpolate | ndarray-interp | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/n10/sorted` | 112 ns | 187 ns | 4.95 µs | 291 ns | 253 ns | 2.3x faster |
| `linear/n10/shuffled` | 169 ns | 220 ns | 5 µs | 282 ns | 248 ns | 1.5x faster |
| `nearest/n10/sorted` | 176 ns | 202 ns | 9.54 µs | 274 ns |  | 1.6x faster |
| `nearest/n10/shuffled` | 296 ns | 370 ns | 9.8 µs | 279 ns |  | **1.1x slower** |
| `previous/n10/sorted` | 173 ns | 207 ns | 10.6 µs |  |  | 61x faster |
| `previous/n10/shuffled` | 260 ns | 319 ns | 10.6 µs |  |  | 41x faster |
| `bspline1/n10/sorted` | 205 ns | 289 ns | 2.62 µs |  |  | 13x faster |
| `bspline1/n10/shuffled` | 527 ns | 656 ns | 2.72 µs |  |  | 5.2x faster |
| `bspline2/n10/sorted` | 266 ns | 372 ns | 3.08 µs |  |  | 12x faster |
| `bspline2/n10/shuffled` | 342 ns | 422 ns | 3.11 µs |  |  | 9.1x faster |
| `bspline3/n10/sorted` | 342 ns | 433 ns | 3.53 µs |  |  | 10x faster |
| `bspline3/n10/shuffled` | 481 ns | 565 ns | 3.58 µs |  |  | 7.4x faster |
| `cubic_spline/n10/sorted` | 140 ns | 209 ns | 3.72 µs | 297 ns | 364 ns | 2.1x faster |
| `cubic_spline/n10/shuffled` | 224 ns | 273 ns | 3.62 µs | 275 ns | 327 ns | 1.2x faster |
| `pchip/n10/sorted` | 138 ns | 193 ns | 3.62 µs | 291 ns |  | 2.1x faster |
| `pchip/n10/shuffled` | 196 ns | 257 ns | 3.69 µs | 296 ns |  | 1.5x faster |
| `akima/n10/sorted` | 149 ns | 211 ns | 3.64 µs | 250 ns |  | 1.7x faster |
| `akima/n10/shuffled` | 205 ns | 258 ns | 3.61 µs | 260 ns |  | 1.3x faster |
| `makima/n10/sorted` | 140 ns | 211 ns | 3.66 µs |  |  | 26x faster |
| `makima/n10/shuffled` | 204 ns | 257 ns | 3.65 µs |  |  | 18x faster |
| `linear/n50/sorted` | 677 ns | 735 ns | 6.11 µs | 1.92 µs | 1.35 µs | 2x faster |
| `linear/n50/shuffled` | 807 ns | 857 ns | 8.14 µs | 1.89 µs | 1.26 µs | 1.6x faster |
| `nearest/n50/sorted` | 763 ns | 808 ns | 11.2 µs | 1.72 µs |  | 2.3x faster |
| `nearest/n50/shuffled` | 1.55 µs | 1.75 µs | 12.4 µs | 1.73 µs |  | 1.1x faster |
| `previous/n50/sorted` | 761 ns | 809 ns | 12.7 µs |  |  | 17x faster |
| `previous/n50/shuffled` | 1.35 µs | 1.41 µs | 12.7 µs |  |  | 9.4x faster |
| `bspline1/n50/sorted` | 1.11 µs | 1.34 µs | 3.94 µs |  |  | 3.5x faster |
| `bspline1/n50/shuffled` | 3.05 µs | 3.13 µs | 4.86 µs |  |  | 1.6x faster |
| `bspline2/n50/sorted` | 1.36 µs | 1.6 µs | 6.17 µs |  |  | 4.5x faster |
| `bspline2/n50/shuffled` | 2.02 µs | 2.09 µs | 6.61 µs |  |  | 3.3x faster |
| `bspline3/n50/sorted` | 1.82 µs | 1.88 µs | 8.63 µs |  |  | 4.7x faster |
| `bspline3/n50/shuffled` | 2.84 µs | 3.14 µs | 10.4 µs |  |  | 3.6x faster |
| `cubic_spline/n50/sorted` | 729 ns | 757 ns | 5.23 µs | 4.69 µs | 1.61 µs | 2.2x faster |
| `cubic_spline/n50/shuffled` | 986 ns | 1.04 µs | 5.99 µs | 4.46 µs | 1.76 µs | 1.8x faster |
| `pchip/n50/sorted` | 720 ns | 822 ns | 5.59 µs | 4.35 µs |  | 6x faster |
| `pchip/n50/shuffled` | 1.07 µs | 1.14 µs | 5.8 µs | 3.71 µs |  | 3.5x faster |
| `akima/n50/sorted` | 798 ns | 825 ns | 5.21 µs | 4.77 µs |  | 6x faster |
| `akima/n50/shuffled` | 981 ns | 1.04 µs | 5.59 µs | 4.05 µs |  | 4.1x faster |
| `makima/n50/sorted` | 740 ns | 833 ns | 5.6 µs |  |  | 7.6x faster |
| `makima/n50/shuffled` | 973 ns | 1.03 µs | 5.66 µs |  |  | 5.8x faster |
| `linear/n100/sorted` | 1.4 µs | 1.47 µs | 7.32 µs | 4.02 µs | 2.46 µs | 1.8x faster |
| `linear/n100/shuffled` | 1.76 µs | 1.72 µs | 12.8 µs | 4.15 µs | 2.49 µs | 1.4x faster |
| `nearest/n100/sorted` | 1.52 µs | 1.61 µs | 13.3 µs | 3.6 µs |  | 2.4x faster |
| `nearest/n100/shuffled` | 3.06 µs | 3.49 µs | 13.7 µs | 4.11 µs |  | 1.3x faster |
| `previous/n100/sorted` | 1.5 µs | 1.58 µs | 14.9 µs |  |  | 9.9x faster |
| `previous/n100/shuffled` | 2.83 µs | 2.94 µs | 15.4 µs |  |  | 5.4x faster |
| `bspline1/n100/sorted` | 2.3 µs | 2.52 µs | 5.44 µs |  |  | 2.4x faster |
| `bspline1/n100/shuffled` | 6.12 µs | 5.77 µs | 11.8 µs |  |  | 1.9x faster |
| `bspline2/n100/sorted` | 2.68 µs | 3.28 µs | 10.8 µs |  |  | 4x faster |
| `bspline2/n100/shuffled` | 4.43 µs | 4.6 µs | 12.6 µs |  |  | 2.9x faster |
| `bspline3/n100/sorted` | 3.91 µs | 3.77 µs | 15 µs |  |  | 3.8x faster |
| `bspline3/n100/shuffled` | 6.15 µs | 6.28 µs | 18.1 µs |  |  | 2.9x faster |
| `cubic_spline/n100/sorted` | 1.44 µs | 1.54 µs | 6.95 µs | 17.3 µs | 3.2 µs | 2.2x faster |
| `cubic_spline/n100/shuffled` | 1.95 µs | 2.05 µs | 7.85 µs | 17.6 µs | 3.19 µs | 1.6x faster |
| `pchip/n100/sorted` | 1.52 µs | 1.53 µs | 6.91 µs | 17.5 µs |  | 4.6x faster |
| `pchip/n100/shuffled` | 1.95 µs | 2.26 µs | 7.72 µs | 16.2 µs |  | 4x faster |
| `akima/n100/sorted` | 1.46 µs | 1.56 µs | 6.81 µs | 18.4 µs |  | 4.7x faster |
| `akima/n100/shuffled` | 1.95 µs | 2.25 µs | 7.9 µs | 17.6 µs |  | 4.1x faster |
| `makima/n100/sorted` | 1.56 µs | 1.64 µs | 6.77 µs |  |  | 4.3x faster |
| `makima/n100/shuffled` | 1.95 µs | 2.24 µs | 8 µs |  |  | 4.1x faster |
| `linear/n1000/sorted` | 12.7 µs | 13.2 µs | 29.3 µs | 120 µs | 24.5 µs | 1.9x faster |
| `linear/n1000/shuffled` | 16.1 µs | 16.4 µs | 155 µs | 215 µs | 24.7 µs | 1.5x faster |
| `nearest/n1000/sorted` | 16.3 µs | 15.5 µs | 55.9 µs | 148 µs |  | 3.4x faster |
| `nearest/n1000/shuffled` | 31.5 µs | 30.9 µs | 217 µs | 236 µs |  | 6.9x faster |
| `previous/n1000/sorted` | 14.9 µs | 15.3 µs | 60 µs |  |  | 4x faster |
| `previous/n1000/shuffled` | 27.2 µs | 27.4 µs | 204 µs |  |  | 7.5x faster |
| `bspline1/n1000/sorted` | 23.8 µs | 25.2 µs | 30.4 µs |  |  | 1.3x faster |
| `bspline1/n1000/shuffled` | 56.3 µs | 62 µs | 531 µs |  |  | 9.4x faster |
| `bspline2/n1000/sorted` | 34.4 µs | 34.9 µs | 76.2 µs |  |  | 2.2x faster |
| `bspline2/n1000/shuffled` | 58.1 µs | 58.4 µs | 549 µs |  |  | 9.4x faster |
| `bspline3/n1000/sorted` | 35.7 µs | 37.2 µs | 139 µs |  |  | 3.9x faster |
| `bspline3/n1000/shuffled` | 78.5 µs | 86 µs | 577 µs |  |  | 7.4x faster |
| `cubic_spline/n1000/sorted` | 14.4 µs | 16.1 µs | 40.3 µs | 1.55 ms | 31.9 µs | 2.2x faster |
| `cubic_spline/n1000/shuffled` | 19.7 µs | 19.8 µs | 219 µs | 1.42 ms | 32 µs | 1.6x faster |
| `pchip/n1000/sorted` | 14.6 µs | 15.3 µs | 39.4 µs | 1.21 ms |  | 2.7x faster |
| `pchip/n1000/shuffled` | 19.8 µs | 20 µs | 219 µs | 1.2 ms |  | 11x faster |
| `akima/n1000/sorted` | 15.7 µs | 15.7 µs | 37.1 µs | 1.42 ms |  | 2.4x faster |
| `akima/n1000/shuffled` | 19.9 µs | 20.2 µs | 218 µs | 1.42 ms |  | 11x faster |
| `makima/n1000/sorted` | 14.4 µs | 14.8 µs | 37.1 µs |  |  | 2.6x faster |
| `makima/n1000/shuffled` | 19.7 µs | 22 µs | 220 µs |  |  | 11x faster |
| `linear/n10000/sorted` | 140 µs | 61.4 µs | 231 µs | 1.7 ms | 247 µs | 1.7x faster |
| `linear/n10000/shuffled` | 190 µs | 79.8 µs | 2.53 ms | 3.08 ms | 259 µs | 1.4x faster |
| `nearest/n10000/sorted` | 155 µs | 71.9 µs | 584 µs | 1.99 ms |  | 3.8x faster |
| `nearest/n10000/shuffled` | 389 µs | 94.3 µs | 3.19 ms | 3.31 ms |  | 8.2x faster |
| `previous/n10000/sorted` | 155 µs | 72.3 µs | 614 µs |  |  | 4x faster |
| `previous/n10000/shuffled` | 341 µs | 85.2 µs | 2.98 ms |  |  | 8.8x faster |
| `bspline1/n10000/sorted` | 283 µs | 97.8 µs | 288 µs |  |  | 1x faster |
| `bspline1/n10000/shuffled` | 635 µs | 140 µs | 49.4 ms |  |  | 78x faster |
| `bspline2/n10000/sorted` | 360 µs | 85.3 µs | 747 µs |  |  | 2.1x faster |
| `bspline2/n10000/shuffled` | 969 µs | 197 µs | 50 ms |  |  | 52x faster |
| `bspline3/n10000/sorted` | 365 µs | 103 µs | 1.26 ms |  |  | 3.4x faster |
| `bspline3/n10000/shuffled` | 1.22 ms | 241 µs | 49.2 ms |  |  | 40x faster |
| `cubic_spline/n10000/sorted` | 154 µs | 71.4 µs | 364 µs | 147 ms | 318 µs | 2.1x faster |
| `cubic_spline/n10000/shuffled` | 217 µs | 94.2 µs | 3.67 ms | 148 ms | 342 µs | 1.6x faster |
| `pchip/n10000/sorted` | 150 µs | 70.5 µs | 379 µs | 122 ms |  | 2.5x faster |
| `pchip/n10000/shuffled` | 217 µs | 95.7 µs | 3.37 ms | 123 ms |  | 16x faster |
| `akima/n10000/sorted` | 164 µs | 70.9 µs | 376 µs | 151 ms |  | 2.3x faster |
| `akima/n10000/shuffled` | 220 µs | 94.8 µs | 3.62 ms | 146 ms |  | 16x faster |
| `makima/n10000/sorted` | 152 µs | 70.1 µs | 363 µs |  |  | 2.4x faster |
| `makima/n10000/shuffled` | 217 µs | 92.9 µs | 3.41 ms |  |  | 16x faster |
| `linear/n100/sorted-uneven` | 1.17 µs | 1.27 µs | 7.33 µs | 3.76 µs | 3.76 µs | 3.2x faster |
| `linear/n100/shuffled-uneven` | 2.73 µs | 2.85 µs | 12.6 µs | 3.84 µs | 3.75 µs | 1.4x faster |
| `cubic_spline/n100/sorted-uneven` | 1.34 µs | 1.45 µs | 7.4 µs | 17.2 µs | 4.44 µs | 3.3x faster |
| `cubic_spline/n100/shuffled-uneven` | 2.9 µs | 3.32 µs | 7.99 µs | 17.8 µs | 4.41 µs | 1.5x faster |
| `pchip/n100/sorted-uneven` | 1.36 µs | 1.44 µs | 7.28 µs | 16.2 µs |  | 5.3x faster |
| `pchip/n100/shuffled-uneven` | 2.92 µs | 3.03 µs | 8.21 µs | 16.2 µs |  | 2.8x faster |
| `linear/n10000/sorted-uneven` | 206 µs | 86.4 µs | 259 µs | 1.68 ms | 890 µs | 1.3x faster |
| `linear/n10000/shuffled-uneven` | 622 µs | 152 µs | 2.43 ms | 3.13 ms | 1.77 ms | 2.8x faster |
| `cubic_spline/n10000/sorted-uneven` | 240 µs | 94.4 µs | 436 µs | 150 ms | 965 µs | 1.8x faster |
| `cubic_spline/n10000/shuffled-uneven` | 733 µs | 165 µs | 3.41 ms | 142 ms | 1.99 ms | 2.7x faster |
| `pchip/n10000/sorted-uneven` | 222 µs | 94.7 µs | 435 µs | 120 ms |  | 2x faster |
| `pchip/n10000/shuffled-uneven` | 735 µs | 168 µs | 3.38 ms | 123 ms |  | 4.6x faster |

### Construction

| case | interpx | interpx, parallel | SciPy | scirs2-interpolate | ndarray-interp | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/n10/sorted` | 95.1 ns | 93 ns | 8.68 µs | 80 ns | 45.4 ns | **2.1x slower** |
| `linear/n10/shuffled` | 101 ns | 98.2 ns | 8.66 µs | 85.6 ns | 40.7 ns | **2.5x slower** |
| `nearest/n10/sorted` | 82.8 ns | 75.1 ns | 9.67 µs | 88.2 ns |  | 1.1x faster |
| `nearest/n10/shuffled` | 75 ns | 81.8 ns | 9.6 µs | 81.2 ns |  | 1.1x faster |
| `previous/n10/sorted` | 75.2 ns | 82.7 ns | 9.24 µs |  |  | 1.2e+02x faster |
| `previous/n10/shuffled` | 75.6 ns | 75.2 ns | 9.22 µs |  |  | 1.2e+02x faster |
| `bspline1/n10/sorted` | 103 ns | 104 ns | 15.4 µs |  |  | 1.5e+02x faster |
| `bspline1/n10/shuffled` | 103 ns | 104 ns | 15.8 µs |  |  | 1.5e+02x faster |
| `bspline2/n10/sorted` | 475 ns | 512 ns | 37.9 µs |  |  | 80x faster |
| `bspline2/n10/shuffled` | 471 ns | 510 ns | 38 µs |  |  | 81x faster |
| `bspline3/n10/sorted` | 513 ns | 511 ns | 36.3 µs |  |  | 71x faster |
| `bspline3/n10/shuffled` | 512 ns | 507 ns | 36.4 µs |  |  | 71x faster |
| `cubic_spline/n10/sorted` | 240 ns | 244 ns | 95.9 µs | 419 ns | 510 ns | 1.7x faster |
| `cubic_spline/n10/shuffled` | 240 ns | 256 ns | 87.6 µs | 420 ns | 473 ns | 1.8x faster |
| `pchip/n10/sorted` | 184 ns | 177 ns | 103 µs | 194 ns |  | 1.1x faster |
| `pchip/n10/shuffled` | 188 ns | 179 ns | 105 µs | 194 ns |  | 1x faster |
| `akima/n10/sorted` | 201 ns | 183 ns | 77.2 µs | 184 ns |  | **1.1x slower** |
| `akima/n10/shuffled` | 184 ns | 183 ns | 76.5 µs | 192 ns |  | 1x faster |
| `makima/n10/sorted` | 192 ns | 222 ns | 80 µs |  |  | 4.2e+02x faster |
| `makima/n10/shuffled` | 212 ns | 194 ns | 80.9 µs |  |  | 3.8e+02x faster |
| `linear/n50/sorted` | 167 ns | 173 ns | 9.76 µs | 127 ns | 108 ns | **1.5x slower** |
| `linear/n50/shuffled` | 156 ns | 156 ns | 8.86 µs | 121 ns | 98.9 ns | **1.6x slower** |
| `nearest/n50/sorted` | 106 ns | 109 ns | 9.92 µs | 123 ns |  | 1.2x faster |
| `nearest/n50/shuffled` | 108 ns | 118 ns | 10.2 µs | 122 ns |  | 1.1x faster |
| `previous/n50/sorted` | 108 ns | 110 ns | 9.55 µs |  |  | 88x faster |
| `previous/n50/shuffled` | 108 ns | 107 ns | 9.49 µs |  |  | 88x faster |
| `bspline1/n50/sorted` | 144 ns | 161 ns | 15.7 µs |  |  | 1.1e+02x faster |
| `bspline1/n50/shuffled` | 158 ns | 148 ns | 16 µs |  |  | 1e+02x faster |
| `bspline2/n50/sorted` | 2.03 µs | 1.98 µs | 40.8 µs |  |  | 20x faster |
| `bspline2/n50/shuffled` | 1.83 µs | 1.83 µs | 41 µs |  |  | 22x faster |
| `bspline3/n50/sorted` | 2.23 µs | 2.05 µs | 40 µs |  |  | 18x faster |
| `bspline3/n50/shuffled` | 2.04 µs | 2.05 µs | 42.7 µs |  |  | 21x faster |
| `cubic_spline/n50/sorted` | 840 ns | 850 ns | 94.9 µs | 945 ns | 1.08 µs | 1.1x faster |
| `cubic_spline/n50/shuffled` | 920 ns | 843 ns | 97.6 µs | 945 ns | 1.18 µs | 1x faster |
| `pchip/n50/sorted` | 496 ns | 547 ns | 112 µs | 365 ns |  | **1.4x slower** |
| `pchip/n50/shuffled` | 546 ns | 544 ns | 105 µs | 334 ns |  | **1.6x slower** |
| `akima/n50/sorted` | 503 ns | 552 ns | 77.9 µs | 361 ns |  | **1.4x slower** |
| `akima/n50/shuffled` | 550 ns | 518 ns | 77.4 µs | 355 ns |  | **1.6x slower** |
| `makima/n50/sorted` | 599 ns | 599 ns | 87 µs |  |  | 1.5e+02x faster |
| `makima/n50/shuffled` | 551 ns | 550 ns | 81.2 µs |  |  | 1.5e+02x faster |
| `linear/n100/sorted` | 252 ns | 272 ns | 10.1 µs | 208 ns | 166 ns | **1.5x slower** |
| `linear/n100/shuffled` | 270 ns | 248 ns | 9.2 µs | 209 ns | 182 ns | **1.5x slower** |
| `nearest/n100/sorted` | 177 ns | 161 ns | 9.99 µs | 197 ns |  | 1.1x faster |
| `nearest/n100/shuffled` | 165 ns | 176 ns | 9.94 µs | 197 ns |  | 1.2x faster |
| `previous/n100/sorted` | 161 ns | 177 ns | 9.86 µs |  |  | 61x faster |
| `previous/n100/shuffled` | 176 ns | 176 ns | 9.88 µs |  |  | 56x faster |
| `bspline1/n100/sorted` | 216 ns | 216 ns | 16 µs |  |  | 74x faster |
| `bspline1/n100/shuffled` | 214 ns | 215 ns | 15.7 µs |  |  | 73x faster |
| `bspline2/n100/sorted` | 3.82 µs | 3.82 µs | 45.7 µs |  |  | 12x faster |
| `bspline2/n100/shuffled` | 3.48 µs | 3.48 µs | 46.9 µs |  |  | 13x faster |
| `bspline3/n100/sorted` | 4.02 µs | 3.89 µs | 45 µs |  |  | 11x faster |
| `bspline3/n100/shuffled` | 3.9 µs | 3.89 µs | 45.9 µs |  |  | 12x faster |
| `cubic_spline/n100/sorted` | 1.58 µs | 1.58 µs | 98.2 µs | 1.55 µs | 1.87 µs | **1x slower** |
| `cubic_spline/n100/shuffled` | 1.58 µs | 1.73 µs | 89.2 µs | 1.52 µs | 1.9 µs | **1x slower** |
| `pchip/n100/sorted` | 968 ns | 884 ns | 105 µs | 538 ns |  | **1.8x slower** |
| `pchip/n100/shuffled` | 882 ns | 927 ns | 108 µs | 539 ns |  | **1.6x slower** |
| `akima/n100/sorted` | 890 ns | 899 ns | 79 µs | 562 ns |  | **1.6x slower** |
| `akima/n100/shuffled` | 959 ns | 940 ns | 77 µs | 523 ns |  | **1.8x slower** |
| `makima/n100/sorted` | 1.07 µs | 1 µs | 79.2 µs |  |  | 74x faster |
| `makima/n100/shuffled` | 976 ns | 1 µs | 82.1 µs |  |  | 84x faster |
| `linear/n1000/sorted` | 1.94 µs | 1.93 µs | 12.9 µs | 1.22 µs | 1.41 µs | **1.6x slower** |
| `linear/n1000/shuffled` | 2.12 µs | 2.12 µs | 13 µs | 1.25 µs | 1.38 µs | **1.7x slower** |
| `nearest/n1000/sorted` | 984 ns | 988 ns | 14.8 µs | 1.22 µs |  | 1.2x faster |
| `nearest/n1000/shuffled` | 1.08 µs | 986 ns | 16.4 µs | 1.28 µs |  | 1.2x faster |
| `previous/n1000/sorted` | 981 ns | 1 µs | 16.5 µs |  |  | 17x faster |
| `previous/n1000/shuffled` | 1.03 µs | 1.07 µs | 17 µs |  |  | 16x faster |
| `bspline1/n1000/sorted` | 1.03 µs | 1.09 µs | 17 µs |  |  | 17x faster |
| `bspline1/n1000/shuffled` | 1.03 µs | 1.13 µs | 17.5 µs |  |  | 17x faster |
| `bspline2/n1000/sorted` | 34.6 µs | 38.1 µs | 110 µs |  |  | 3.2x faster |
| `bspline2/n1000/shuffled` | 36.6 µs | 37.9 µs | 119 µs |  |  | 3.3x faster |
| `bspline3/n1000/sorted` | 66 µs | 63.7 µs | 138 µs |  |  | 2.1x faster |
| `bspline3/n1000/shuffled` | 68.4 µs | 67.6 µs | 129 µs |  |  | 1.9x faster |
| `cubic_spline/n1000/sorted` | 38.4 µs | 38.6 µs | 125 µs | 16.9 µs | 15.8 µs | **2.4x slower** |
| `cubic_spline/n1000/shuffled` | 38.5 µs | 39.4 µs | 114 µs | 16.2 µs | 17.2 µs | **2.4x slower** |
| `pchip/n1000/sorted` | 12.8 µs | 11.9 µs | 122 µs | 4.3 µs |  | **3x slower** |
| `pchip/n1000/shuffled` | 11.5 µs | 11.7 µs | 123 µs | 4.16 µs |  | **2.8x slower** |
| `akima/n1000/sorted` | 12.6 µs | 12.8 µs | 97.9 µs | 6.49 µs |  | **1.9x slower** |
| `akima/n1000/shuffled` | 11.6 µs | 9.08 µs | 96.8 µs | 5.63 µs |  | **2.1x slower** |
| `makima/n1000/sorted` | 12.5 µs | 12.5 µs | 102 µs |  |  | 8.2x faster |
| `makima/n1000/shuffled` | 12.5 µs | 12.6 µs | 102 µs |  |  | 8.1x faster |
| `linear/n10000/sorted` | 90.5 µs | 79.1 µs | 47.8 µs | 58.7 µs | 19.9 µs | **4.6x slower** |
| `linear/n10000/shuffled` | 92.6 µs | 86.2 µs | 45.9 µs | 56.5 µs | 19.5 µs | **4.8x slower** |
| `nearest/n10000/sorted` | 64 µs | 51.4 µs | 56 µs | 58.3 µs |  | **1.1x slower** |
| `nearest/n10000/shuffled` | 67 µs | 52.8 µs | 58.9 µs | 59.5 µs |  | **1.1x slower** |
| `previous/n10000/sorted` | 62.9 µs | 49.4 µs | 79 µs |  |  | 1.3x faster |
| `previous/n10000/shuffled` | 62.2 µs | 56.5 µs | 76.9 µs |  |  | 1.2x faster |
| `bspline1/n10000/sorted` | 66.8 µs | 66.6 µs | 25.7 µs |  |  | **2.6x slower** |
| `bspline1/n10000/shuffled` | 60.8 µs | 61 µs | 23.4 µs |  |  | **2.6x slower** |
| `bspline2/n10000/sorted` | 492 µs | 521 µs | 818 µs |  |  | 1.7x faster |
| `bspline2/n10000/shuffled` | 527 µs | 487 µs | 758 µs |  |  | 1.4x faster |
| `bspline3/n10000/sorted` | 624 µs | 580 µs | 935 µs |  |  | 1.5x faster |
| `bspline3/n10000/shuffled` | 626 µs | 623 µs | 946 µs |  |  | 1.5x faster |
| `cubic_spline/n10000/sorted` | 236 µs | 336 µs | 739 µs | 341 µs | 360 µs | 1.4x faster |
| `cubic_spline/n10000/shuffled` | 256 µs | 340 µs | 796 µs | 345 µs | 346 µs | 1.3x faster |
| `pchip/n10000/sorted` | 123 µs | 91.6 µs | 244 µs | 154 µs |  | 1.3x faster |
| `pchip/n10000/shuffled` | 119 µs | 90.6 µs | 252 µs | 163 µs |  | 1.4x faster |
| `akima/n10000/sorted` | 94.1 µs | 210 µs | 626 µs | 140 µs |  | 1.5x faster |
| `akima/n10000/shuffled` | 112 µs | 91.2 µs | 618 µs | 130 µs |  | 1.2x faster |
| `makima/n10000/sorted` | 102 µs | 211 µs | 746 µs |  |  | 7.3x faster |
| `makima/n10000/shuffled` | 99.7 µs | 226 µs | 705 µs |  |  | 7.1x faster |
| `linear/n100/sorted-uneven` | 248 ns | 248 ns | 9.12 µs | 206 ns | 166 ns | **1.5x slower** |
| `linear/n100/shuffled-uneven` | 256 ns | 249 ns | 9.1 µs | 193 ns | 166 ns | **1.5x slower** |
| `cubic_spline/n100/sorted-uneven` | 1.73 µs | 1.58 µs | 89.7 µs | 1.52 µs | 1.88 µs | **1.1x slower** |
| `cubic_spline/n100/shuffled-uneven` | 1.58 µs | 1.58 µs | 91.4 µs | 1.52 µs | 1.87 µs | **1x slower** |
| `pchip/n100/sorted-uneven` | 882 ns | 882 ns | 105 µs | 593 ns |  | **1.5x slower** |
| `pchip/n100/shuffled-uneven` | 968 ns | 888 ns | 104 µs | 543 ns |  | **1.8x slower** |
| `linear/n10000/sorted-uneven` | 92.1 µs | 84.8 µs | 47.9 µs | 59.8 µs | 21.4 µs | **4.3x slower** |
| `linear/n10000/shuffled-uneven` | 89.9 µs | 81.8 µs | 48.5 µs | 58.3 µs | 20 µs | **4.5x slower** |
| `cubic_spline/n10000/sorted-uneven` | 249 µs | 329 µs | 713 µs | 287 µs | 354 µs | 1.2x faster |
| `cubic_spline/n10000/shuffled-uneven` | 246 µs | 308 µs | 761 µs | 269 µs | 355 µs | 1.1x faster |
| `pchip/n10000/sorted-uneven` | 109 µs | 96.1 µs | 241 µs | 157 µs |  | 1.4x faster |
| `pchip/n10000/shuffled-uneven` | 113 µs | 203 µs | 240 µs | 168 µs |  | 1.5x faster |

## `PPoly`

A cubic of 55 intervals, evaluated at points in no order.

### Evaluation

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `intervals55/eval10` | 72.2 ns | 79.8 ns | 3.38 µs | 47x faster |
| `intervals55/eval100` | 652 ns | 733 ns | 4.29 µs | 6.6x faster |
| `intervals55/eval1000` | 6.48 µs | 6.9 µs | 15.5 µs | 2.4x faster |
| `intervals55/eval10000` | 66.5 µs | 67.4 µs | 335 µs | 5x faster |

### Construction

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `intervals55/eval10` | 160 ns | 174 ns | 5.23 µs | 33x faster |
| `intervals55/eval100` | 160 ns | 158 ns | 5.21 µs | 33x faster |
| `intervals55/eval1000` | 160 ns | 159 ns | 5.22 µs | 33x faster |
| `intervals55/eval10000` | 159 ns | 173 ns | 5.12 µs | 32x faster |

## `BSpline`

A cubic B-spline of 55 intervals, evaluated at points in no order.

### Evaluation

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `intervals55/eval10` | 146 ns | 151 ns | 2.72 µs | 19x faster |
| `intervals55/eval100` | 1.49 µs | 1.64 µs | 5.77 µs | 3.9x faster |
| `intervals55/eval1000` | 14.1 µs | 15.7 µs | 37.2 µs | 2.6x faster |
| `intervals55/eval10000` | 142 µs | 52.8 µs | 359 µs | 2.5x faster |

### Construction

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `intervals55/eval10` | 176 ns | 161 ns | 9.25 µs | 53x faster |
| `intervals55/eval100` | 161 ns | 161 ns | 9.23 µs | 57x faster |
| `intervals55/eval1000` | 161 ns | 160 ns | 9.27 µs | 58x faster |
| `intervals55/eval10000` | 176 ns | 176 ns | 9.34 µs | 53x faster |

## `BSplineCurve`

A cubic curve through `n` points of the plane, evaluated at `4 n` parameters in order.

### Evaluation

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `degree3/n100` | 4.18 µs | 4.27 µs | 15.9 µs | 3.8x faster |
| `degree3/n10000` | 449 µs | 152 µs | 1.34 ms | 3x faster |

### Construction

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `degree3/n100` | 4.25 µs | 4.25 µs | 48.9 µs | 12x faster |
| `degree3/n10000` | 629 µs | 677 µs | 1.37 ms | 2.2x faster |

## Regular grid

Points in no order inside the grid. `maxN`: the first axis has `N` points and each axis half the one before (SciPy's suite); `eachN`: every axis has `N` points; `flipped`: every axis given from its high end.

### Evaluation

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | ndarray-interp | interpn | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/2d/max10/eval10/ascending` | 134 ns | 136 ns | 138 ns | 19.2 µs | 811 ns | 152 ns | 94.8 ns | **1.4x slower** |
| `linear/2d/max10/eval10/flipped` | 132 ns | 136 ns | 134 ns | 19 µs |  |  |  | 1.4e+02x faster |
| `linear/2d/max10/eval100/ascending` | 1.25 µs | 1.43 µs | 1.24 µs | 21.5 µs | 8.58 µs | 1.36 µs | 928 ns | **1.3x slower** |
| `linear/2d/max10/eval100/flipped` | 1.36 µs | 1.38 µs | 1.26 µs | 21.5 µs |  |  |  | 16x faster |
| `linear/2d/max10/eval1000/ascending` | 12.5 µs | 13.8 µs | 12.5 µs | 40.6 µs | 85 µs | 14.4 µs | 9.23 µs | **1.4x slower** |
| `linear/2d/max10/eval1000/flipped` | 12.5 µs | 12.6 µs | 12.7 µs | 40.1 µs |  |  |  | 3.2x faster |
| `linear/2d/max10/eval10000/ascending` | 137 µs | 47.8 µs | 125 µs | 371 µs | 894 µs | 220 µs | 92.4 µs | **1.5x slower** |
| `linear/2d/max10/eval10000/flipped` | 125 µs | 50.6 µs | 137 µs | 362 µs |  |  |  | 2.9x faster |
| `linear/2d/max40/eval10/ascending` | 172 ns | 160 ns | 172 ns | 21 µs | 1.04 µs | 156 ns | 147 ns | **1.2x slower** |
| `linear/2d/max40/eval10/flipped` | 170 ns | 174 ns | 173 ns | 19.7 µs |  |  |  | 1.2e+02x faster |
| `linear/2d/max40/eval100/ascending` | 1.66 µs | 1.59 µs | 1.66 µs | 22.1 µs | 11.2 µs | 1.59 µs | 1.45 µs | **1.1x slower** |
| `linear/2d/max40/eval100/flipped` | 1.61 µs | 1.74 µs | 1.66 µs | 23.9 µs |  |  |  | 15x faster |
| `linear/2d/max40/eval1000/ascending` | 15.1 µs | 16.9 µs | 16.1 µs | 46.6 µs | 117 µs | 17 µs | 13.2 µs | **1.1x slower** |
| `linear/2d/max40/eval1000/flipped` | 16.6 µs | 16.8 µs | 15.1 µs | 47.8 µs |  |  |  | 2.9x faster |
| `linear/2d/max40/eval10000/ascending` | 167 µs | 56.2 µs | 152 µs | 577 µs | 1.18 ms | 243 µs | 147 µs | **1.1x slower** |
| `linear/2d/max40/eval10000/flipped` | 152 µs | 62.5 µs | 153 µs | 598 µs |  |  |  | 3.9x faster |
| `linear/2d/max200/eval10/ascending` | 212 ns | 196 ns | 212 ns | 20.7 µs | 1.96 µs | 205 ns | 193 ns | **1.1x slower** |
| `linear/2d/max200/eval10/flipped` | 197 ns | 197 ns | 197 ns | 19.2 µs |  |  |  | 98x faster |
| `linear/2d/max200/eval100/ascending` | 1.89 µs | 1.98 µs | 1.89 µs | 22.8 µs | 21.9 µs | 2.18 µs | 1.73 µs | **1.1x slower** |
| `linear/2d/max200/eval100/flipped` | 2.07 µs | 2.17 µs | 2.09 µs | 24.5 µs |  |  |  | 12x faster |
| `linear/2d/max200/eval1000/ascending` | 20.7 µs | 19.3 µs | 20.7 µs | 72.2 µs | 230 µs | 24.2 µs | 18.7 µs | **1.1x slower** |
| `linear/2d/max200/eval1000/flipped` | 20.8 µs | 21.2 µs | 20.8 µs | 70.9 µs |  |  |  | 3.4x faster |
| `linear/2d/max200/eval10000/ascending` | 189 µs | 75.7 µs | 208 µs | 823 µs | 2.16 ms | 483 µs | 173 µs | **1.1x slower** |
| `linear/2d/max200/eval10000/flipped` | 194 µs | 81 µs | 207 µs | 834 µs |  |  |  | 4.3x faster |
| `linear/3d/max10/eval10/ascending` | 218 ns | 204 ns | 198 ns | 56.8 µs | 3.16 µs |  | 142 ns | **1.5x slower** |
| `linear/3d/max10/eval10/flipped` | 198 ns | 203 ns | 199 ns | 56.2 µs |  |  |  | 2.8e+02x faster |
| `linear/3d/max10/eval100/ascending` | 1.95 µs | 1.98 µs | 1.91 µs | 63.4 µs | 33.7 µs |  | 1.53 µs | **1.3x slower** |
| `linear/3d/max10/eval100/flipped` | 1.91 µs | 1.99 µs | 1.98 µs | 63.3 µs |  |  |  | 33x faster |
| `linear/3d/max10/eval1000/ascending` | 20.9 µs | 19.2 µs | 19.1 µs | 116 µs | 339 µs |  | 15.3 µs | **1.4x slower** |
| `linear/3d/max10/eval1000/flipped` | 19.1 µs | 21 µs | 19 µs | 119 µs |  |  |  | 6.2x faster |
| `linear/3d/max10/eval10000/ascending` | 193 µs | 79.9 µs | 191 µs | 1.1 ms | 3.54 ms |  | 142 µs | **1.4x slower** |
| `linear/3d/max10/eval10000/flipped` | 191 µs | 73.3 µs | 210 µs | 1.13 ms |  |  |  | 5.9x faster |
| `linear/3d/max40/eval10/ascending` | 222 ns | 228 ns | 240 ns | 62.3 µs | 3.39 µs |  | 199 ns | **1.1x slower** |
| `linear/3d/max40/eval10/flipped` | 241 ns | 247 ns | 241 ns | 61.7 µs |  |  |  | 2.6e+02x faster |
| `linear/3d/max40/eval100/ascending` | 2.16 µs | 2.42 µs | 2.35 µs | 63.2 µs | 33.3 µs |  | 1.98 µs | **1.1x slower** |
| `linear/3d/max40/eval100/flipped` | 2.14 µs | 2.22 µs | 2.35 µs | 64 µs |  |  |  | 30x faster |
| `linear/3d/max40/eval1000/ascending` | 23.5 µs | 21.7 µs | 23.5 µs | 136 µs | 341 µs |  | 20 µs | **1.2x slower** |
| `linear/3d/max40/eval1000/flipped` | 21.5 µs | 23.6 µs | 23.5 µs | 141 µs |  |  |  | 6.6x faster |
| `linear/3d/max40/eval10000/ascending` | 213 µs | 79.4 µs | 214 µs | 1.47 ms | 3.71 ms |  | 219 µs | 1x faster |
| `linear/3d/max40/eval10000/flipped` | 214 µs | 75.4 µs | 214 µs | 1.43 ms |  |  |  | 6.7x faster |
| `linear/3d/max200/eval10/ascending` | 278 ns | 309 ns | 283 ns | 56.6 µs | 4.5 µs |  | 261 ns | **1.1x slower** |
| `linear/3d/max200/eval10/flipped` | 279 ns | 291 ns | 292 ns | 57 µs |  |  |  | 2e+02x faster |
| `linear/3d/max200/eval100/ascending` | 2.73 µs | 3.08 µs | 3 µs | 69 µs | 50.3 µs |  | 2.86 µs | 1.1x faster |
| `linear/3d/max200/eval100/flipped` | 2.73 µs | 3.07 µs | 2.73 µs | 71.6 µs |  |  |  | 26x faster |
| `linear/3d/max200/eval1000/ascending` | 29.8 µs | 30.2 µs | 29.9 µs | 188 µs | 500 µs |  | 26.1 µs | **1.1x slower** |
| `linear/3d/max200/eval1000/flipped` | 29.8 µs | 30.2 µs | 29.8 µs | 203 µs |  |  |  | 6.8x faster |
| `linear/3d/max200/eval10000/ascending` | 298 µs | 85.4 µs | 299 µs | 1.86 ms | 5.06 ms |  | 263 µs | **1.1x slower** |
| `linear/3d/max200/eval10000/flipped` | 281 µs | 90.8 µs | 273 µs | 1.88 ms |  |  |  | 6.7x faster |
| `linear/4d/max10/eval10/ascending` | 273 ns | 283 ns | 260 ns | 96.4 µs |  |  |  | 3.5e+02x faster |
| `linear/4d/max10/eval10/flipped` | 272 ns | 281 ns | 260 ns | 97.3 µs |  |  |  | 3.6e+02x faster |
| `linear/4d/max10/eval100/ascending` | 2.4 µs | 2.69 µs | 2.62 µs | 108 µs |  |  |  | 45x faster |
| `linear/4d/max10/eval100/flipped` | 2.42 µs | 2.71 µs | 2.42 µs | 107 µs |  |  |  | 44x faster |
| `linear/4d/max10/eval1000/ascending` | 26.8 µs | 26.9 µs | 24.4 µs | 215 µs |  |  |  | 8x faster |
| `linear/4d/max10/eval1000/flipped` | 26.5 µs | 26.9 µs | 25.3 µs | 218 µs |  |  |  | 8.2x faster |
| `linear/4d/max10/eval10000/ascending` | 246 µs | 82 µs | 268 µs | 1.98 ms |  |  |  | 8.1x faster |
| `linear/4d/max10/eval10000/flipped` | 267 µs | 88.3 µs | 245 µs | 1.95 ms |  |  |  | 7.3x faster |
| `linear/4d/max40/eval10/ascending` | 328 ns | 338 ns | 337 ns | 96.3 µs | 6.3 µs |  | 298 ns | **1.1x slower** |
| `linear/4d/max40/eval10/flipped` | 327 ns | 335 ns | 337 ns | 106 µs |  |  |  | 3.2e+02x faster |
| `linear/4d/max40/eval100/ascending` | 2.95 µs | 3.02 µs | 3.21 µs | 118 µs | 65.4 µs |  | 2.97 µs | 1x faster |
| `linear/4d/max40/eval100/flipped` | 2.96 µs | 3.29 µs | 3.21 µs | 111 µs |  |  |  | 38x faster |
| `linear/4d/max40/eval1000/ascending` | 29.6 µs | 32.3 µs | 32 µs | 278 µs | 594 µs |  | 29.5 µs | **1x slower** |
| `linear/4d/max40/eval1000/flipped` | 32.1 µs | 32.4 µs | 29.4 µs | 258 µs |  |  |  | 8.1x faster |
| `linear/4d/max40/eval10000/ascending` | 295 µs | 88.9 µs | 318 µs | 2.19 ms | 6.06 ms |  | 296 µs | 1x faster |
| `linear/4d/max40/eval10000/flipped` | 298 µs | 90.5 µs | 323 µs | 2.4 ms |  |  |  | 8x faster |
| `linear/4d/max200/eval10/ascending` | 389 ns | 384 ns | 415 ns | 105 µs | 8.48 µs |  | 394 ns | 1x faster |
| `linear/4d/max200/eval10/flipped` | 384 ns | 384 ns | 384 ns | 106 µs |  |  |  | 2.7e+02x faster |
| `linear/4d/max200/eval100/ascending` | 3.75 µs | 3.85 µs | 4.11 µs | 120 µs | 79.9 µs |  | 3.97 µs | 1.1x faster |
| `linear/4d/max200/eval100/flipped` | 3.78 µs | 3.83 µs | 3.77 µs | 121 µs |  |  |  | 32x faster |
| `linear/4d/max200/eval1000/ascending` | 41.2 µs | 41.9 µs | 41 µs | 311 µs | 732 µs |  | 39.9 µs | **1x slower** |
| `linear/4d/max200/eval1000/flipped` | 41.7 µs | 38.4 µs | 37.9 µs | 335 µs |  |  |  | 8x faster |
| `linear/4d/max200/eval10000/ascending` | 374 µs | 106 µs | 412 µs | 2.86 ms | 7.32 ms |  | 400 µs | 1.1x faster |
| `linear/4d/max200/eval10000/flipped` | 378 µs | 108 µs | 411 µs | 2.91 ms |  |  |  | 7.7x faster |
| `nearest/2d/each40/eval10000/ascending` | 151 µs | 64.2 µs | 152 µs | 652 µs | 1.74 ms |  | 141 µs | **1.1x slower** |
| `nearest/3d/each40/eval10000/ascending` | 212 µs | 77.3 µs | 233 µs | 981 µs | 2.11 ms |  | 209 µs | **1x slower** |
| `pchip/2d/each40/eval10000/ascending` | 494 µs | 130 µs | 492 µs | 1.28 s |  |  |  | 2.6e+03x faster |
| `pchip/3d/each40/eval10000/ascending` | 1.47 ms | 272 µs | 1.48 ms | 3.49 s |  |  |  | 2.4e+03x faster |
| `cubic/2d/each40/eval10000/ascending` | 325 µs | 104 µs | 301 µs | 1.28 ms |  |  | 410 µs | 1.3x faster |
| `cubic/3d/each40/eval10000/ascending` | 668 µs | 183 µs | 724 µs | 3.03 ms |  |  | 790 µs | 1.2x faster |
| `cubic/2d/each10/eval1000/ascending` | 29.2 µs | 29.4 µs | 29.1 µs | 144 µs |  |  | 39.2 µs | 1.3x faster |
| `cubic/2d/each40/eval1000/ascending` | 32.9 µs | 33.1 µs | 32.9 µs | 155 µs |  |  | 44.6 µs | 1.4x faster |
| `cubic/2d/each100/eval1000/ascending` | 35.2 µs | 32.3 µs | 32.1 µs | 168 µs |  |  | 47.1 µs | 1.3x faster |
| `cubic/2d/each200/eval1000/ascending` | 37.8 µs | 38.1 µs | 37.7 µs | 219 µs |  |  | 50.2 µs | 1.3x faster |
| `cubic/2d/each400/eval1000/ascending` | 41.3 µs | 42.1 µs | 41.1 µs | 276 µs |  |  | 50 µs | 1.2x faster |
| `quintic/2d/each10/eval1000/ascending` | 44.2 µs | 46.3 µs | 44.5 µs | 229 µs |  |  |  | 5.2x faster |
| `quintic/2d/each40/eval1000/ascending` | 48 µs | 51.2 µs | 48.1 µs | 243 µs |  |  |  | 5.1x faster |
| `linear/2d/each400/eval10000/ascending` | 242 µs | 84.1 µs | 242 µs | 967 µs | 4.23 ms | 593 µs | 241 µs | **1x slower** |
| `linear/2d/each800/eval10000/ascending` | 240 µs | 85.2 µs | 262 µs | 1.08 ms | 7.78 ms | 631 µs | 245 µs | 1x faster |
| `nearest/2d/auckland-ba32-0401/heldout/ascending` | 364 ms | 61.8 ms | 348 ms | 1.34 s | 158 s |  | 750 ms | 2.1x faster |
| `linear/2d/auckland-ba32-0401/heldout/ascending` | 377 ms | 65 ms | 370 ms | 1.09 s | 73.2 s | 369 ms | 772 ms | **1x slower** |
| `cubic/2d/auckland-ba32-0401/heldout/ascending` | 723 ms | 124 ms | 736 ms | 31.8 s |  |  | 1.62 s | 2.2x faster |
| `quintic/2d/auckland-ba32-0401/heldout/ascending` | 1.19 s | 207 ms | 1.2 s | 35.3 s |  |  |  | 30x faster |
| `pchip/2d/auckland-ba32-0401/heldout/ascending` | 1.02 s | 172 ms | 1.03 s |  |  |  |  |  |

### Construction

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | ndarray-interp | interpn | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/2d/max10/eval10/ascending` | 121 ns | 132 ns | 135 ns | 5.8 µs | 141 ns | 78.4 ns | 0 | **1.8e+02x slower** |
| `linear/2d/max10/eval10/flipped` | 154 ns | 154 ns | 152 ns | 11.6 µs |  |  |  | 75x faster |
| `linear/2d/max10/eval100/ascending` | 121 ns | 121 ns | 117 ns | 5.74 µs | 137 ns | 77.9 ns | 0 | **1.6e+02x slower** |
| `linear/2d/max10/eval100/flipped` | 170 ns | 164 ns | 168 ns | 11.7 µs |  |  |  | 69x faster |
| `linear/2d/max10/eval1000/ascending` | 121 ns | 133 ns | 127 ns | 5.81 µs | 148 ns | 77.9 ns | 0 | **1.8e+02x slower** |
| `linear/2d/max10/eval1000/flipped` | 156 ns | 167 ns | 152 ns | 11.7 µs |  |  |  | 75x faster |
| `linear/2d/max10/eval10000/ascending` | 132 ns | 123 ns | 129 ns | 6.27 µs | 138 ns | 78.4 ns | 0 | **1.7e+02x slower** |
| `linear/2d/max10/eval10000/flipped` | 168 ns | 169 ns | 165 ns | 11.7 µs |  |  |  | 70x faster |
| `linear/2d/max40/eval10/ascending` | 362 ns | 377 ns | 712 ns | 6.37 µs | 240 ns | 167 ns | 0 | **5.3e+02x slower** |
| `linear/2d/max40/eval10/flipped` | 706 ns | 668 ns | 1.25 µs | 12.8 µs |  |  |  | 18x faster |
| `linear/2d/max40/eval100/ascending` | 394 ns | 359 ns | 646 ns | 6.32 µs | 235 ns | 169 ns | 0 | **5e+02x slower** |
| `linear/2d/max40/eval100/flipped` | 707 ns | 704 ns | 1.16 µs | 12.6 µs |  |  |  | 18x faster |
| `linear/2d/max40/eval1000/ascending` | 402 ns | 390 ns | 644 ns | 5.76 µs | 215 ns | 167 ns | 0 | **5.1e+02x slower** |
| `linear/2d/max40/eval1000/flipped` | 715 ns | 706 ns | 1.17 µs | 11.9 µs |  |  |  | 17x faster |
| `linear/2d/max40/eval10000/ascending` | 392 ns | 356 ns | 641 ns | 5.88 µs | 234 ns | 152 ns | 0 | **5e+02x slower** |
| `linear/2d/max40/eval10000/flipped` | 706 ns | 700 ns | 1.19 µs | 11.7 µs |  |  |  | 17x faster |
| `linear/2d/max200/eval10/ascending` | 7.42 µs | 7.53 µs | 16.4 µs | 6.06 µs | 2.93 µs | 472 ns | 0 | **9.5e+03x slower** |
| `linear/2d/max200/eval10/flipped` | 15.9 µs | 15.6 µs | 26.2 µs | 12.8 µs |  |  |  | **1.2x slower** |
| `linear/2d/max200/eval100/ascending` | 6.81 µs | 6.83 µs | 16.5 µs | 5.82 µs | 3.22 µs | 517 ns | 0 | **8.7e+03x slower** |
| `linear/2d/max200/eval100/flipped` | 15.9 µs | 16.1 µs | 26.7 µs | 12.3 µs |  |  |  | **1.3x slower** |
| `linear/2d/max200/eval1000/ascending` | 7.44 µs | 6.96 µs | 17.5 µs | 5.85 µs | 3.26 µs | 517 ns | 0 | **9.5e+03x slower** |
| `linear/2d/max200/eval1000/flipped` | 15.7 µs | 15.9 µs | 27.3 µs | 11.8 µs |  |  |  | **1.3x slower** |
| `linear/2d/max200/eval10000/ascending` | 42.2 µs | 46.8 µs | 20.5 µs | 5.96 µs | 54.2 µs | 517 ns | 0 | **6.2e+04x slower** |
| `linear/2d/max200/eval10000/flipped` | 57 µs | 60.8 µs | 28.8 µs | 13 µs |  |  |  | **4.4x slower** |
| `linear/3d/max10/eval10/ascending` | 196 ns | 179 ns | 205 ns | 7.48 µs | 191 ns |  | 0 | **2.2e+02x slower** |
| `linear/3d/max10/eval10/flipped` | 298 ns | 296 ns | 337 ns | 16.9 µs |  |  |  | 57x faster |
| `linear/3d/max10/eval100/ascending` | 189 ns | 181 ns | 207 ns | 7.71 µs | 174 ns |  | 0 | **2.2e+02x slower** |
| `linear/3d/max10/eval100/flipped` | 301 ns | 307 ns | 340 ns | 17.1 µs |  |  |  | 57x faster |
| `linear/3d/max10/eval1000/ascending` | 193 ns | 198 ns | 209 ns | 7.58 µs | 190 ns |  | 0 | **2e+02x slower** |
| `linear/3d/max10/eval1000/flipped` | 325 ns | 297 ns | 325 ns | 17.2 µs |  |  |  | 53x faster |
| `linear/3d/max10/eval10000/ascending` | 180 ns | 200 ns | 204 ns | 8.23 µs | 191 ns |  | 0 | **2e+02x slower** |
| `linear/3d/max10/eval10000/flipped` | 301 ns | 323 ns | 341 ns | 15.9 µs |  |  |  | 53x faster |
| `linear/3d/max40/eval10/ascending` | 2.97 µs | 6.67 µs | 5.64 µs | 8.22 µs | 9.23 µs |  | 0 | **3.4e+03x slower** |
| `linear/3d/max40/eval10/flipped` | 8.9 µs | 9.02 µs | 12.5 µs | 17.3 µs |  |  |  | 1.9x faster |
| `linear/3d/max40/eval100/ascending` | 5.25 µs | 6.93 µs | 6.36 µs | 7.55 µs | 6.66 µs |  | 0 | **6e+03x slower** |
| `linear/3d/max40/eval100/flipped` | 8.23 µs | 8.19 µs | 12.6 µs | 15.9 µs |  |  |  | 1.9x faster |
| `linear/3d/max40/eval1000/ascending` | 17 µs | 17.1 µs | 7.22 µs | 8.15 µs | 8.46 µs |  | 0 | **1.9e+04x slower** |
| `linear/3d/max40/eval1000/flipped` | 30.2 µs | 28.3 µs | 12.9 µs | 17.3 µs |  |  |  | **1.7x slower** |
| `linear/3d/max40/eval10000/ascending` | 14.6 µs | 16.6 µs | 5.64 µs | 7.65 µs | 13.5 µs |  | 0 | **1.7e+04x slower** |
| `linear/3d/max40/eval10000/flipped` | 26.6 µs | 28.2 µs | 11.8 µs | 15.9 µs |  |  |  | **1.7x slower** |
| `linear/3d/max200/eval10/ascending` | 1.79 ms | 1.75 ms | 643 µs | 8.23 µs | 1.29 ms |  | 0 | **2e+06x slower** |
| `linear/3d/max200/eval10/flipped` | 3.21 ms | 3.45 ms | 2.05 ms | 16.1 µs |  |  |  | **2e+02x slower** |
| `linear/3d/max200/eval100/ascending` | 1.71 ms | 1.81 ms | 652 µs | 7.64 µs | 1.27 ms |  | 0 | **1.8e+06x slower** |
| `linear/3d/max200/eval100/flipped` | 3.08 ms | 3.47 ms | 1.9 ms | 17.6 µs |  |  |  | **1.8e+02x slower** |
| `linear/3d/max200/eval1000/ascending` | 1.7 ms | 1.8 ms | 653 µs | 7.62 µs | 1.35 ms |  | 0 | **1.9e+06x slower** |
| `linear/3d/max200/eval1000/flipped` | 3.23 ms | 3.36 ms | 1.95 ms | 16.2 µs |  |  |  | **2e+02x slower** |
| `linear/3d/max200/eval10000/ascending` | 1.82 ms | 1.67 ms | 636 µs | 7.51 µs | 1.27 ms |  | 0 | **2.1e+06x slower** |
| `linear/3d/max200/eval10000/flipped` | 3.2 ms | 3.4 ms | 1.85 ms | 16 µs |  |  |  | **2e+02x slower** |
| `linear/4d/max10/eval10/ascending` | 235 ns | 217 ns | 251 ns | 9.56 µs |  |  |  | 41x faster |
| `linear/4d/max10/eval10/flipped` | 375 ns | 370 ns | 388 ns | 17.7 µs |  |  |  | 47x faster |
| `linear/4d/max10/eval100/ascending` | 216 ns | 235 ns | 247 ns | 9.61 µs |  |  |  | 45x faster |
| `linear/4d/max10/eval100/flipped` | 343 ns | 377 ns | 397 ns | 19.5 µs |  |  |  | 57x faster |
| `linear/4d/max10/eval1000/ascending` | 237 ns | 238 ns | 233 ns | 9.57 µs |  |  |  | 40x faster |
| `linear/4d/max10/eval1000/flipped` | 380 ns | 376 ns | 415 ns | 19.3 µs |  |  |  | 51x faster |
| `linear/4d/max10/eval10000/ascending` | 217 ns | 238 ns | 228 ns | 10.8 µs |  |  |  | 50x faster |
| `linear/4d/max10/eval10000/flipped` | 373 ns | 344 ns | 337 ns | 17.8 µs |  |  |  | 48x faster |
| `linear/4d/max40/eval10/ascending` | 14.4 µs | 14.5 µs | 29.5 µs | 9.42 µs | 5.98 µs |  | 1.37 ns | **1.1e+04x slower** |
| `linear/4d/max40/eval10/flipped` | 56.4 µs | 52.1 µs | 69.9 µs | 21.6 µs |  |  |  | **2.6x slower** |
| `linear/4d/max40/eval100/ascending` | 14.5 µs | 14.7 µs | 35.9 µs | 10.2 µs | 6.08 µs |  | 1.37 ns | **1.1e+04x slower** |
| `linear/4d/max40/eval100/flipped` | 56.5 µs | 52.5 µs | 64.8 µs | 20.6 µs |  |  |  | **2.7x slower** |
| `linear/4d/max40/eval1000/ascending` | 14.6 µs | 14.3 µs | 30.8 µs | 10.2 µs | 5.52 µs |  | 1.27 ns | **1.1e+04x slower** |
| `linear/4d/max40/eval1000/flipped` | 56.4 µs | 56.4 µs | 65.2 µs | 21.8 µs |  |  |  | **2.6x slower** |
| `linear/4d/max40/eval10000/ascending` | 13.1 µs | 14.6 µs | 32.2 µs | 9.39 µs | 5.99 µs |  | 1.37 ns | **9.6e+03x slower** |
| `linear/4d/max40/eval10000/flipped` | 56.8 µs | 139 µs | 75.4 µs | 19.8 µs |  |  |  | **2.9x slower** |
| `linear/4d/max200/eval10/ascending` | 52.5 ms | 49.4 ms | 15.6 ms | 9.56 µs | 41.5 ms |  | 1.37 ns | **3.8e+07x slower** |
| `linear/4d/max200/eval10/flipped` | 104 ms | 108 ms | 72.3 ms | 22 µs |  |  |  | **4.7e+03x slower** |
| `linear/4d/max200/eval100/ascending` | 52.2 ms | 48.4 ms | 15.4 ms | 10.4 µs | 40.1 ms |  | 1.37 ns | **3.8e+07x slower** |
| `linear/4d/max200/eval100/flipped` | 108 ms | 109 ms | 79.3 ms | 20.4 µs |  |  |  | **5.3e+03x slower** |
| `linear/4d/max200/eval1000/ascending` | 52 ms | 49.7 ms | 16.7 ms | 10.1 µs | 40.6 ms |  | 1.27 ns | **4.1e+07x slower** |
| `linear/4d/max200/eval1000/flipped` | 127 ms | 104 ms | 71 ms | 21.9 µs |  |  |  | **5.8e+03x slower** |
| `linear/4d/max200/eval10000/ascending` | 49.3 ms | 49.4 ms | 16.5 ms | 9.31 µs | 37.7 ms |  | 1.37 ns | **3.6e+07x slower** |
| `linear/4d/max200/eval10000/flipped` | 102 ms | 116 ms | 74.6 ms | 20.4 µs |  |  |  | **5e+03x slower** |
| `nearest/2d/each40/eval10000/ascending` | 638 ns | 635 ns | 1.5 µs | 5.83 µs | 296 ns |  | 0 | **9.3e+02x slower** |
| `nearest/3d/each40/eval10000/ascending` | 22.8 µs | 21.2 µs | 49.9 µs | 7.5 µs | 8.89 µs |  | 0 | **2.3e+04x slower** |
| `pchip/2d/each40/eval10000/ascending` | 6.5 µs | 7.12 µs | 6.3 µs | 7.26 µs |  |  |  | 1.1x faster |
| `pchip/3d/each40/eval10000/ascending` | 508 µs | 593 µs | 455 µs | 8.63 µs |  |  |  | **59x slower** |
| `cubic/2d/each40/eval10000/ascending` | 8.38 µs | 8.34 µs | 7.98 µs | 1.88 ms |  |  | 0 | **9.5e+03x slower** |
| `cubic/3d/each40/eval10000/ascending` | 324 µs | 471 µs | 327 µs | 149 ms |  |  | 1.17 ns | **2.8e+05x slower** |
| `cubic/2d/each10/eval1000/ascending` | 1.58 µs | 1.58 µs | 1.62 µs | 1.18 ms |  |  | 0 | **1.8e+03x slower** |
| `cubic/2d/each40/eval1000/ascending` | 9.06 µs | 8.38 µs | 8.73 µs | 2.03 ms |  |  | 0 | **9.3e+03x slower** |
| `cubic/2d/each100/eval1000/ascending` | 81.7 µs | 105 µs | 80.7 µs | 7.36 ms |  |  | 0 | **8.4e+04x slower** |
| `cubic/2d/each200/eval1000/ascending` | 186 µs | 186 µs | 227 µs | 32.5 ms |  |  | 0 | **2.1e+05x slower** |
| `cubic/2d/each400/eval1000/ascending` | 1.4 ms | 1.5 ms | 1.19 ms | 156 ms |  |  | 0 | **1.6e+06x slower** |
| `quintic/2d/each10/eval1000/ascending` | 2.08 µs | 2.1 µs | 2.16 µs | 4.16 ms |  |  |  | 2e+03x faster |
| `quintic/2d/each40/eval1000/ascending` | 12.2 µs | 12 µs | 13 µs | 7.63 ms |  |  |  | 6.3e+02x faster |
| `linear/2d/each400/eval10000/ascending` | 319 µs | 338 µs | 106 µs | 5.97 µs | 280 µs | 1.17 µs | 0 | **4.1e+05x slower** |
| `linear/2d/each800/eval10000/ascending` | 1.09 ms | 1.11 ms | 415 µs | 6.38 µs | 838 µs | 2.44 µs | 0 | **1.6e+06x slower** |
| `nearest/2d/auckland-ba32-0401/heldout/ascending` | 16.5 ms | 15.4 ms | 5.52 ms | 7.69 µs | 12.8 ms |  | 0 | **2.1e+07x slower** |
| `linear/2d/auckland-ba32-0401/heldout/ascending` | 16.1 ms | 15.5 ms | 5.45 ms | 7.04 µs | 13.6 ms | 9.03 µs | 0 | **2.3e+07x slower** |
| `cubic/2d/auckland-ba32-0401/heldout/ascending` | 77.3 ms | 73.7 ms | 66.6 ms | 14.3 s |  |  | 0 | **9.9e+07x slower** |
| `quintic/2d/auckland-ba32-0401/heldout/ascending` | 90.9 ms | 95.8 ms | 83.6 ms | 20.3 s |  |  |  | 2.2e+02x faster |
| `pchip/2d/auckland-ba32-0401/heldout/ascending` | 72.1 ms | 72.7 ms | 62 ms |  |  |  |  |  |

## `NearestNDInterpolator`

`n` points of the unit square, evaluated at the pixels of a square image.

### Evaluation

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: |
| `n10/image50` | 35 µs | 32.5 µs | 122 µs | 443 µs | 13x faster |
| `n10/image100` | 124 µs | 44.7 µs | 444 µs | 1.8 ms | 15x faster |
| `n10/image500` | 3.05 ms | 547 µs | 12.2 ms | 47.2 ms | 15x faster |
| `n100/image50` | 116 µs | 112 µs | 204 µs | 592 µs | 5.1x faster |
| `n100/image100` | 424 µs | 118 µs | 781 µs | 2.23 ms | 5.2x faster |
| `n100/image500` | 9.71 ms | 1.56 ms | 18.1 ms | 56.8 ms | 5.9x faster |
| `n1000/image50` | 342 µs | 107 µs | 341 µs | 832 µs | 2.4x faster |
| `n1000/image100` | 1.23 ms | 289 µs | 1.18 ms | 2.99 ms | 2.4x faster |
| `n1000/image500` | 25.2 ms | 4.03 ms | 25.3 ms | 71.1 ms | 2.8x faster |
| `n10000/image50` | 329 µs | 109 µs | 361 µs | 1.08 ms | 3.3x faster |
| `n10000/image100` | 1.48 ms | 281 µs | 1.43 ms | 3.78 ms | 2.5x faster |
| `n10000/image500` | 32.6 ms | 5.25 ms | 32.5 ms | 87.5 ms | 2.7x faster |
| `n100000/image50` | 475 µs | 115 µs | 441 µs | 1.43 ms | 3x faster |
| `n100000/image100` | 1.77 ms | 367 µs | 1.81 ms | 5.6 ms | 3.2x faster |
| `n100000/image500` | 44.4 ms | 6.75 ms | 43.2 ms | 126 ms | 2.8x faster |
| `ghcn-tmax-20250101/heldout` | 159 µs | 66.5 µs | 171 µs | 513 µs | 3.2x faster |

### Construction

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: |
| `n10/image50` | 538 ns | 575 ns | 522 ns | 9.53 µs | 18x faster |
| `n10/image100` | 529 ns | 525 ns | 515 ns | 9.53 µs | 18x faster |
| `n10/image500` | 585 ns | 573 ns | 570 ns | 9.39 µs | 16x faster |
| `n100/image50` | 3.61 µs | 3.44 µs | 3.41 µs | 15.1 µs | 4.2x faster |
| `n100/image100` | 3.68 µs | 3.68 µs | 3.7 µs | 15.3 µs | 4.1x faster |
| `n100/image500` | 3.65 µs | 3.67 µs | 3.63 µs | 15.1 µs | 4.1x faster |
| `n1000/image50` | 58.4 µs | 63.3 µs | 49.1 µs | 127 µs | 2.2x faster |
| `n1000/image100` | 53 µs | 59.2 µs | 40.6 µs | 119 µs | 2.2x faster |
| `n1000/image500` | 55.4 µs | 61.6 µs | 47.6 µs | 120 µs | 2.2x faster |
| `n10000/image50` | 684 µs | 668 µs | 631 µs | 1.72 ms | 2.5x faster |
| `n10000/image100` | 679 µs | 692 µs | 714 µs | 1.73 ms | 2.5x faster |
| `n10000/image500` | 673 µs | 676 µs | 654 µs | 1.84 ms | 2.7x faster |
| `n100000/image50` | 8 ms | 8.14 ms | 8.38 ms | 24.7 ms | 3.1x faster |
| `n100000/image100` | 8.48 ms | 8.4 ms | 8.54 ms | 25.1 ms | 3x faster |
| `n100000/image500` | 8.42 ms | 8.06 ms | 8.18 ms | 26.2 ms | 3.1x faster |
| `ghcn-tmax-20250101/heldout` | 710 µs | 670 µs | 721 µs | 1.66 ms | 2.3x faster |

## `RBFInterpolator`

`n` points of a square, evaluated at 1000 points in no order (`scattered`) or at the 10 000 pixels of an image (`image`). `global`: all points; `k50`: the 50 nearest.

### Evaluation

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/global/n10/scattered` | 15.8 µs | 15.9 µs | 15.8 µs | 63.7 µs | 197 µs | 4x faster |
| `linear/k50/n10/scattered` | 15.8 µs | 15.9 µs | 15.8 µs | 1.35 ms |  | 86x faster |
| `thin_plate_spline/global/n10/scattered` | 77.3 µs | 78.3 µs | 71.2 µs | 112 µs | 324 µs | 1.4x faster |
| `thin_plate_spline/k50/n10/scattered` | 77.8 µs | 78.4 µs | 71.2 µs | 1.55 ms |  | 20x faster |
| `cubic/global/n10/scattered` | 17.9 µs | 16.8 µs | 17.9 µs | 62.2 µs | 215 µs | 3.5x faster |
| `cubic/k50/n10/scattered` | 18.1 µs | 17.4 µs | 16.6 µs | 1.34 ms |  | 74x faster |
| `quintic/global/n10/scattered` | 25.4 µs | 23.7 µs | 23.6 µs | 71.3 µs | 234 µs | 2.8x faster |
| `quintic/k50/n10/scattered` | 23.4 µs | 22.9 µs | 25.3 µs | 1.52 ms |  | 65x faster |
| `multiquadric/global/n10/scattered` | 27.8 µs | 27.9 µs | 27.8 µs | 53 µs | 257 µs | 1.9x faster |
| `multiquadric/k50/n10/scattered` | 27.8 µs | 27.9 µs | 30.5 µs | 1.36 ms |  | 49x faster |
| `inverse_multiquadric/global/n10/scattered` | 37.7 µs | 37.9 µs | 41.1 µs | 54.9 µs | 268 µs | 1.5x faster |
| `inverse_multiquadric/k50/n10/scattered` | 37.4 µs | 37.6 µs | 41.1 µs | 1.37 ms |  | 37x faster |
| `inverse_quadratic/global/n10/scattered` | 23.9 µs | 24.4 µs | 23.8 µs | 47.2 µs |  | 2x faster |
| `inverse_quadratic/k50/n10/scattered` | 23.8 µs | 24.2 µs | 26.2 µs | 1.48 ms |  | 62x faster |
| `gaussian/global/n10/scattered` | 73.6 µs | 67.3 µs | 67.2 µs | 89.8 µs | 356 µs | 1.2x faster |
| `gaussian/k50/n10/scattered` | 67.2 µs | 67.2 µs | 73.7 µs | 1.39 ms |  | 21x faster |
| `linear/global/n100/scattered` | 135 µs | 49.3 µs | 146 µs | 291 µs | 1.77 ms | 2.1x faster |
| `linear/k50/n100/scattered` | 14 ms | 2.46 ms | 13 ms | 18.6 ms |  | 1.3x faster |
| `thin_plate_spline/global/n100/scattered` | 695 µs | 178 µs | 762 µs | 872 µs | 2.92 ms | 1.3x faster |
| `thin_plate_spline/k50/n100/scattered` | 25.6 ms | 4.07 ms | 24.8 ms | 25.1 ms |  | **1x slower** |
| `cubic/global/n100/scattered` | 150 µs | 52.4 µs | 139 µs | 291 µs | 1.89 ms | 1.9x faster |
| `cubic/k50/n100/scattered` | 15.6 ms | 2.78 ms | 14.8 ms | 22.1 ms |  | 1.4x faster |
| `quintic/global/n100/scattered` | 165 µs | 64.6 µs | 181 µs | 301 µs | 2.12 ms | 1.8x faster |
| `quintic/k50/n100/scattered` | 18 ms | 3.03 ms | 17.3 ms | 23.9 ms |  | 1.3x faster |
| `multiquadric/global/n100/scattered` | 269 µs | 80 µs | 269 µs | 331 µs | 2.36 ms | 1.2x faster |
| `multiquadric/k50/n100/scattered` | 16.2 ms | 2.82 ms | 15.6 ms | 19.9 ms |  | 1.2x faster |
| `inverse_multiquadric/global/n100/scattered` | 364 µs | 108 µs | 364 µs | 398 µs | 2.6 ms | 1.1x faster |
| `inverse_multiquadric/k50/n100/scattered` | 16.1 ms | 2.98 ms | 16.3 ms | 19 ms |  | 1.2x faster |
| `inverse_quadratic/global/n100/scattered` | 224 µs | 72.4 µs | 224 µs | 264 µs |  | 1.2x faster |
| `inverse_quadratic/k50/n100/scattered` | 14.2 ms | 2.77 ms | 16 ms | 18.7 ms |  | 1.3x faster |
| `gaussian/global/n100/scattered` | 620 µs | 153 µs | 620 µs | 685 µs | 3.63 ms | 1.1x faster |
| `gaussian/k50/n100/scattered` | 21.4 ms | 3.88 ms | 23.8 ms | 22.7 ms |  | 1.1x faster |
| `linear/global/n1000/scattered` | 1.47 ms | 322 µs | 1.34 ms | 3.7 ms | 17.9 ms | 2.5x faster |
| `linear/k50/n1000/scattered` | 13.6 ms | 2.57 ms | 14 ms | 34.9 ms |  | 2.6x faster |
| `thin_plate_spline/global/n1000/scattered` | 7.25 ms | 1.4 ms | 7.34 ms | 9.85 ms | 30.9 ms | 1.4x faster |
| `thin_plate_spline/k50/n1000/scattered` | 22.4 ms | 3.97 ms | 23.3 ms | 44.1 ms |  | 2x faster |
| `cubic/global/n1000/scattered` | 1.35 ms | 314 µs | 1.36 ms | 3.98 ms | 19.7 ms | 2.9x faster |
| `cubic/k50/n1000/scattered` | 15.9 ms | 2.97 ms | 16.9 ms | 36.8 ms |  | 2.3x faster |
| `quintic/global/n1000/scattered` | 1.52 ms | 328 µs | 1.52 ms | 4.08 ms | 21.9 ms | 2.7x faster |
| `quintic/k50/n1000/scattered` | 18.6 ms | 3.15 ms | 18.6 ms | 43.8 ms |  | 2.4x faster |
| `multiquadric/global/n1000/scattered` | 2.7 ms | 603 µs | 2.76 ms | 4.37 ms | 23.3 ms | 1.6x faster |
| `multiquadric/k50/n1000/scattered` | 16 ms | 2.92 ms | 16 ms | 35.2 ms |  | 2.2x faster |
| `inverse_multiquadric/global/n1000/scattered` | 3.7 ms | 802 µs | 3.7 ms | 4.95 ms | 26.9 ms | 1.3x faster |
| `inverse_multiquadric/k50/n1000/scattered` | 17.1 ms | 3.27 ms | 18.8 ms | 36.4 ms |  | 2.1x faster |
| `inverse_quadratic/global/n1000/scattered` | 2.47 ms | 509 µs | 2.45 ms | 3.83 ms |  | 1.6x faster |
| `inverse_quadratic/k50/n1000/scattered` | 15.3 ms | 3.04 ms | 15.5 ms | 35.3 ms |  | 2.3x faster |
| `gaussian/global/n1000/scattered` | 6.59 ms | 1.19 ms | 6.27 ms | 8.13 ms | 37.9 ms | 1.2x faster |
| `gaussian/k50/n1000/scattered` | 21.2 ms | 3.97 ms | 23.3 ms | 43.2 ms |  | 2x faster |
| `thin_plate_spline/global/n3000/scattered` | 21.4 ms | 4.34 ms | 21.7 ms | 31.8 ms | 90.8 ms | 1.5x faster |
| `thin_plate_spline/k50/n10000/scattered` | 22.6 ms | 4.17 ms | 22.4 ms | 44.9 ms |  | 2x faster |
| `thin_plate_spline/k50/n10000/image` | 231 ms | 37.2 ms | 229 ms | 452 ms |  | 2x faster |
| `thin_plate_spline/k50/n100000/scattered` | 22.9 ms | 4.24 ms | 22.6 ms | 45.4 ms |  | 2x faster |
| `thin_plate_spline/k50/n100000/image` | 231 ms | 36.6 ms | 239 ms | 452 ms |  | 2x faster |
| `linear/global/ghcn-tmax-20250101/heldout` | 19 ms | 3.34 ms | 20.3 ms | 53.2 ms | 229 ms | 2.8x faster |
| `linear/k50/ghcn-tmax-20250101/heldout` | 17.7 ms | 3.35 ms | 17.8 ms | 44.3 ms |  | 2.5x faster |
| `thin_plate_spline/global/ghcn-tmax-20250101/heldout` | 84 ms | 79.4 ms | 85.7 ms | 112 ms | 413 ms | 1.3x faster |
| `thin_plate_spline/k50/ghcn-tmax-20250101/heldout` | 28.8 ms | 5.14 ms | 28.6 ms | 54.7 ms |  | 1.9x faster |
| `cubic/global/ghcn-tmax-20250101/heldout` | 19.6 ms | 3.4 ms | 17.8 ms | 51.2 ms | 264 ms | 2.6x faster |
| `cubic/k50/ghcn-tmax-20250101/heldout` | 21.2 ms | 3.61 ms | 19.8 ms | 49.1 ms |  | 2.3x faster |
| `quintic/global/ghcn-tmax-20250101/heldout` | 19.9 ms | 3.99 ms | 22.2 ms | 58.5 ms | 295 ms | 2.9x faster |
| `quintic/k50/ghcn-tmax-20250101/heldout` | 22.9 ms | 3.91 ms | 23.2 ms | 51.7 ms |  | 2.3x faster |

### Construction

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/global/n10/scattered` | 2.84 µs | 2.84 µs | 2.82 µs | 11.7 µs | 5.9 µs | 2.1x faster |
| `linear/k50/n10/scattered` | 2.83 µs | 2.83 µs | 2.8 µs | 15.4 µs |  | 5.4x faster |
| `thin_plate_spline/global/n10/scattered` | 3.59 µs | 3.54 µs | 3.28 µs | 13.2 µs | 7.11 µs | 2x faster |
| `thin_plate_spline/k50/n10/scattered` | 3.6 µs | 3.58 µs | 3.34 µs | 16.3 µs |  | 4.5x faster |
| `cubic/global/n10/scattered` | 3.82 µs | 3.8 µs | 3.73 µs | 13.1 µs | 5.51 µs | 1.4x faster |
| `cubic/k50/n10/scattered` | 3.66 µs | 3.73 µs | 3.75 µs | 15.7 µs |  | 4.3x faster |
| `quintic/global/n10/scattered` | 3.9 µs | 3.97 µs | 4.3 µs | 15.2 µs | 6.23 µs | 1.6x faster |
| `quintic/k50/n10/scattered` | 3.91 µs | 3.87 µs | 4.24 µs | 18.5 µs |  | 4.7x faster |
| `multiquadric/global/n10/scattered` | 2.68 µs | 2.65 µs | 2.66 µs | 11 µs | 6.6 µs | 2.5x faster |
| `multiquadric/k50/n10/scattered` | 2.64 µs | 2.64 µs | 2.67 µs | 14.8 µs |  | 5.6x faster |
| `inverse_multiquadric/global/n10/scattered` | 2.88 µs | 2.73 µs | 2.63 µs | 10.8 µs | 6.87 µs | 2.4x faster |
| `inverse_multiquadric/k50/n10/scattered` | 2.86 µs | 2.62 µs | 2.64 µs | 16.4 µs |  | 5.7x faster |
| `inverse_quadratic/global/n10/scattered` | 2.82 µs | 2.6 µs | 2.57 µs | 10.7 µs |  | 3.8x faster |
| `inverse_quadratic/k50/n10/scattered` | 2.78 µs | 2.79 µs | 2.8 µs | 15.1 µs |  | 5.4x faster |
| `gaussian/global/n10/scattered` | 2.76 µs | 2.74 µs | 2.74 µs | 11.2 µs | 7.74 µs | 2.8x faster |
| `gaussian/k50/n10/scattered` | 2.75 µs | 2.75 µs | 2.76 µs | 14.9 µs |  | 5.4x faster |
| `linear/global/n100/scattered` | 74.9 µs | 91.4 µs | 82.1 µs | 78.3 µs | 731 µs | 1x faster |
| `linear/k50/n100/scattered` | 3.45 µs | 3.45 µs | 3.4 µs | 21 µs |  | 6.1x faster |
| `thin_plate_spline/global/n100/scattered` | 97.3 µs | 122 µs | 105 µs | 112 µs | 842 µs | 1.2x faster |
| `thin_plate_spline/k50/n100/scattered` | 3.46 µs | 3.29 µs | 3.2 µs | 21.7 µs |  | 6.3x faster |
| `cubic/global/n100/scattered` | 77.6 µs | 72 µs | 89.8 µs | 80.4 µs | 756 µs | 1x faster |
| `cubic/k50/n100/scattered` | 3.32 µs | 3.49 µs | 3.45 µs | 21.9 µs |  | 6.6x faster |
| `quintic/global/n100/scattered` | 85.5 µs | 87 µs | 90.8 µs | 92 µs | 763 µs | 1.1x faster |
| `quintic/k50/n100/scattered` | 3.27 µs | 3.54 µs | 3.19 µs | 23.2 µs |  | 7.1x faster |
| `multiquadric/global/n100/scattered` | 75.7 µs | 96.8 µs | 88.1 µs | 76.7 µs | 621 µs | 1x faster |
| `multiquadric/k50/n100/scattered` | 3.34 µs | 3.45 µs | 3.45 µs | 20.8 µs |  | 6.2x faster |
| `inverse_multiquadric/global/n100/scattered` | 79.6 µs | 84 µs | 90.3 µs | 80.1 µs | 876 µs | 1x faster |
| `inverse_multiquadric/k50/n100/scattered` | 3.22 µs | 3.43 µs | 3.2 µs | 21.2 µs |  | 6.6x faster |
| `inverse_quadratic/global/n100/scattered` | 78.8 µs | 83.7 µs | 84.2 µs | 73 µs |  | **1.1x slower** |
| `inverse_quadratic/k50/n100/scattered` | 3.16 µs | 3.45 µs | 3.43 µs | 22.6 µs |  | 7.1x faster |
| `gaussian/global/n100/scattered` | 96.2 µs | 116 µs | 103 µs | 105 µs | 797 µs | 1.1x faster |
| `gaussian/k50/n100/scattered` | 3.38 µs | 3.36 µs | 3.12 µs | 20.7 µs |  | 6.1x faster |
| `linear/global/n1000/scattered` | 18.9 ms | 19.1 ms | 19.5 ms | 57.9 ms | 207 ms | 3.1x faster |
| `linear/k50/n1000/scattered` | 59.2 µs | 54.5 µs | 55.4 µs | 129 µs |  | 2.2x faster |
| `thin_plate_spline/global/n1000/scattered` | 23.1 ms | 22.1 ms | 22.4 ms | 65.7 ms | 216 ms | 2.8x faster |
| `thin_plate_spline/k50/n1000/scattered` | 61.1 µs | 58 µs | 55.2 µs | 132 µs |  | 2.2x faster |
| `cubic/global/n1000/scattered` | 19.6 ms | 19.5 ms | 19.8 ms | 61.8 ms | 189 ms | 3.2x faster |
| `cubic/k50/n1000/scattered` | 62.4 µs | 56.3 µs | 55.3 µs | 134 µs |  | 2.2x faster |
| `quintic/global/n1000/scattered` | 19.7 ms | 19.7 ms | 19.4 ms | 58.3 ms | 195 ms | 3x faster |
| `quintic/k50/n1000/scattered` | 60.1 µs | 57.4 µs | 58.4 µs | 133 µs |  | 2.2x faster |
| `multiquadric/global/n1000/scattered` | 20.4 ms | 20.8 ms | 20.6 ms | 60.2 ms | 189 ms | 2.9x faster |
| `multiquadric/k50/n1000/scattered` | 60.2 µs | 56.5 µs | 52.8 µs | 130 µs |  | 2.2x faster |
| `inverse_multiquadric/global/n1000/scattered` | 20.7 ms | 20.5 ms | 21.1 ms | 61.6 ms | 192 ms | 3x faster |
| `inverse_multiquadric/k50/n1000/scattered` | 59.1 µs | 56.7 µs | 55.9 µs | 129 µs |  | 2.2x faster |
| `inverse_quadratic/global/n1000/scattered` | 21.1 ms | 20.5 ms | 19.8 ms | 60.4 ms |  | 2.9x faster |
| `inverse_quadratic/k50/n1000/scattered` | 58.4 µs | 51.6 µs | 54.7 µs | 141 µs |  | 2.4x faster |
| `gaussian/global/n1000/scattered` | 24.1 ms | 22.3 ms | 22.8 ms | 63.1 ms | 213 ms | 2.6x faster |
| `gaussian/k50/n1000/scattered` | 61.1 µs | 53.5 µs | 54.5 µs | 130 µs |  | 2.1x faster |
| `thin_plate_spline/global/n3000/scattered` | 453 ms | 438 ms | 442 ms | 2.35 s | 6.05 s | 5.2x faster |
| `thin_plate_spline/k50/n10000/scattered` | 769 µs | 786 µs | 674 µs | 1.76 ms |  | 2.3x faster |
| `thin_plate_spline/k50/n10000/image` | 777 µs | 794 µs | 696 µs | 1.66 ms |  | 2.1x faster |
| `thin_plate_spline/k50/n100000/scattered` | 8.47 ms | 9.18 ms | 8.19 ms | 24.3 ms |  | 2.9x faster |
| `thin_plate_spline/k50/n100000/image` | 8.6 ms | 8.74 ms | 8.36 ms | 24.1 ms |  | 2.8x faster |
| `linear/global/ghcn-tmax-20250101/heldout` | 17.6 s | 16.6 s | 17 s | 121 s | 268 s | 6.9x faster |
| `linear/k50/ghcn-tmax-20250101/heldout` | 821 µs | 781 µs | 703 µs | 1.69 ms |  | 2.1x faster |
| `thin_plate_spline/global/ghcn-tmax-20250101/heldout` | 17 s | 45 s | 18.4 s | 109 s | 282 s | 6.5x faster |
| `thin_plate_spline/k50/ghcn-tmax-20250101/heldout` | 813 µs | 808 µs | 709 µs | 1.73 ms |  | 2.1x faster |
| `cubic/global/ghcn-tmax-20250101/heldout` | 17.5 s | 17.1 s | 17.3 s | 113 s | 298 s | 6.5x faster |
| `cubic/k50/ghcn-tmax-20250101/heldout` | 772 µs | 839 µs | 717 µs | 1.73 ms |  | 2.2x faster |
| `quintic/global/ghcn-tmax-20250101/heldout` | 18 s | 17.2 s | 18.1 s | 115 s | 282 s | 6.4x faster |
| `quintic/k50/ghcn-tmax-20250101/heldout` | 824 µs | 846 µs | 783 µs | 1.75 ms |  | 2.1x faster |

## Memory

Megabytes. "Data read" is the process with the data of the case in memory, before anything is built. The columns after it are what was added to that: the highest while building, with the interpolator built, the highest while evaluating. "Highest" is the highest of the whole run, repeats included.

| case | library | data read | while building | built | while evaluating | highest |
| :--- | :--- | ---: | ---: | ---: | ---: | ---: |
| `oned/cubic_spline/n10000/sorted` | interpx | 3.3 | +1.9 | +2.0 | +2.0 | 6.0 |
| `oned/cubic_spline/n10000/sorted` | interpx, parallel | 3.9 | +1.3 | +1.4 | +1.4 | 6.4 |
| `oned/cubic_spline/n10000/sorted` | SciPy | 71.6 | +0.3 | +1.8 | +1.8 | 73.5 |
| `oned/cubic_spline/n10000/sorted` | scirs2-interpolate | 3.3 | +1.8 | +1.9 | +2.6 | 5.9 |
| `oned/cubic_spline/n10000/sorted` | ndarray-interp | 3.1 | +1.9 | +1.8 | +1.9 | 5.5 |
| `grid/linear/4d/max200/eval10000/ascending` | interpx | 206.2 | +198.9 | +198.9 | +198.9 | 405.1 |
| `grid/linear/4d/max200/eval10000/ascending` | interpx, parallel | 205.4 | +199.7 | +199.4 | +199.4 | 405.4 |
| `grid/linear/4d/max200/eval10000/ascending` | interpx, `from_owned` | 239.8 | +165.3 | +165.3 | +165.3 | 405.1 |
| `grid/linear/4d/max200/eval10000/ascending` | SciPy | 271.3 | +0.1 | +0.1 | +2.0 | 273.4 |
| `grid/linear/4d/max200/eval10000/ascending` | scirs2-interpolate | 224.9 | +178.7 | +180.0 | +180.0 | 405.0 |
| `grid/linear/4d/max200/eval10000/ascending` | interpn | 230.9 | +0.0 | +0.0 | +0.0 | 396.6 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpx | 5.6 | +0.1 | +2.7 | +2.7 | 8.5 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpx, parallel | 5.4 | +0.4 | +3.0 | +3.0 | 8.9 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpx, `from_owned` | 5.6 | +2.7 | +2.7 | +2.7 | 8.5 |
| `grid/cubic/2d/each400/eval1000/ascending` | SciPy | 72.3 | +85.6 | +6.3 | +6.3 | 160.6 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpn | 5.3 | +0.1 | +0.1 | +0.1 | 5.6 |
| `grid/linear/2d/each800/eval10000/ascending` | interpx | 10.9 | +1.9 | +3.9 | +3.9 | 15.0 |
| `grid/linear/2d/each800/eval10000/ascending` | interpx, parallel | 14.6 | +0.0 | +0.5 | +0.5 | 15.4 |
| `grid/linear/2d/each800/eval10000/ascending` | interpx, `from_owned` | 11.2 | +0.3 | +3.8 | +3.8 | 15.0 |
| `grid/linear/2d/each800/eval10000/ascending` | SciPy | 76.0 | +0.3 | +0.4 | +0.7 | 77.0 |
| `grid/linear/2d/each800/eval10000/ascending` | scirs2-interpolate | 12.4 | +0.0 | +2.5 | +2.6 | 15.0 |
| `grid/linear/2d/each800/eval10000/ascending` | ndarray-interp | 13.8 | +0.0 | +0.0 | +0.0 | 13.8 |
| `grid/linear/2d/each800/eval10000/ascending` | interpn | 12.3 | +0.0 | +0.0 | +0.0 | 12.3 |
| `nearest_nd/n100000/image500` | interpx | 11.8 | +3.7 | +2.5 | +4.4 | 17.6 |
| `nearest_nd/n100000/image500` | interpx, parallel | 13.4 | +1.8 | +0.9 | +3.4 | 18.8 |
| `nearest_nd/n100000/image500` | interpx, `from_owned` | 13.3 | +2.6 | +0.9 | +2.9 | 18.2 |
| `nearest_nd/n100000/image500` | SciPy | 76.6 | +5.5 | +5.5 | +14.2 | 91.9 |
| `rbf/thin_plate_spline/global/n3000/scattered` | interpx | 3.0 | +77.4 | +5.1 | +5.1 | 80.4 |
| `rbf/thin_plate_spline/global/n3000/scattered` | interpx, parallel | 3.4 | +77.0 | +4.8 | +5.1 | 80.9 |
| `rbf/thin_plate_spline/global/n3000/scattered` | interpx, `from_owned` | 3.3 | +77.2 | +4.9 | +4.9 | 80.5 |
| `rbf/thin_plate_spline/global/n3000/scattered` | SciPy | 71.2 | +73.0 | +0.9 | +8.9 | 144.3 |
| `rbf/thin_plate_spline/global/n3000/scattered` | scirs2-interpolate | 3.9 | +288.6 | +0.5 | +0.6 | 292.5 |
| `rbf/thin_plate_spline/k50/n100000/image` | interpx | 6.6 | +3.6 | +4.8 | +5.1 | 11.9 |
| `rbf/thin_plate_spline/k50/n100000/image` | interpx, parallel | 5.9 | +4.5 | +5.5 | +6.6 | 12.6 |
| `rbf/thin_plate_spline/k50/n100000/image` | interpx, `from_owned` | 6.5 | +6.4 | +4.9 | +5.3 | 13.5 |
| `rbf/thin_plate_spline/k50/n100000/image` | SciPy | 73.4 | +5.3 | +5.5 | +21.8 | 99.3 |

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

`rbf/thin_plate_spline/global/n3000/scattered`

![memory over time, rbf/thin_plate_spline/global/n3000/scattered](performance/rbf-thin_plate_spline-global-n3000-scattered.svg)

`rbf/thin_plate_spline/k50/n100000/image`

![memory over time, rbf/thin_plate_spline/k50/n100000/image](performance/rbf-thin_plate_spline-k50-n100000-image.svg)

## Do the libraries compute the same thing

The largest distance of a library's values from SciPy's, over the cases of a family, as a share of the largest value. Under 1e-06 the two compute the same interpolant.

| family | interpx | interpx, parallel | interpx, `from_owned` | scirs2-interpolate | ndarray-interp | interpn |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| One variable | 3e-16 | 3e-16 |  | 6e-02 | 1e-04 |  |
| `PPoly` | 9e-16 | 9e-16 |  |  |  |  |
| `BSpline` | 4e-16 | 4e-16 |  |  |  |  |
| `BSplineCurve` | 4e-16 | 4e-16 |  |  |  |  |
| Regular grid | 7e-04 | 7e-04 | 7e-04 | 4e-16 | 6e-16 | 7e-01 |
| `NearestNDInterpolator` | 0e+00 | 0e+00 | 0e+00 |  |  |  |
| `RBFInterpolator` | 4e+00 | 4e+00 | 4e+00 | 7e+11 |  |  |
