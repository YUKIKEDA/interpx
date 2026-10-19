# Performance: interpx against SciPy and other Rust crates

<!-- Written by compare/perf/report.py from a run of compare/perf/run.py. Do not edit by hand: run it again. -->

Measured on 2026-10-18, at commit `a2498de5a6`.

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
| One variable | construction | 112 | 76 | 0 | 36 |
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
| `CloughTocher2DInterpolator` | evaluation | 7 | 7 | 0 | 0 |
| `CloughTocher2DInterpolator` | construction | 7 | 5 | 0 | 2 |
| `RBFInterpolator` | evaluation | 61 | 61 | 0 | 0 |
| `RBFInterpolator` | construction | 61 | 61 | 0 | 0 |

214 comparisons (a case, an operation, and a library), by method. "Same values" says whether that library's values are SciPy's up to rounding: where they are not, it computes another interpolant under the same name, and its time is not that of the same work.

| method | operation | against | same values | cases | slower by | the largest |
| :--- | :--- | :--- | :--- | ---: | ---: | :--- |
| `grid/cubic` | construction | interpn | no | 8 | 906x to 9.42e+07x | `grid/cubic/2d/auckland-ba32-0401/heldout/ascending`: 82.8 ms against 0.88 ns |
| `grid/linear` | construction | interpn | yes | 35 | 155x to 4.08e+07x | `grid/linear/4d/max200/eval10/ascending`: 51.8 ms against 1.27 ns |
| `grid/nearest` | construction | interpn | yes | 3 | 802x to 2.43e+07x | `grid/nearest/2d/auckland-ba32-0401/heldout/ascending`: 16.6 ms against 0.68 ns |
| `grid/linear` | construction | SciPy | yes | 35 | 1.16x to 5.54e+03x | `grid/linear/4d/max200/eval10/ascending`: 51.8 ms against 9.36 µs |
| `grid/nearest` | construction | SciPy | yes | 2 | 2.57x to 2.39e+03x | `grid/nearest/2d/auckland-ba32-0401/heldout/ascending`: 16.6 ms against 6.94 µs |
| `grid/linear` | construction | ndarray-interp | yes | 15 | 1.55x to 2.01e+03x | `grid/linear/2d/auckland-ba32-0401/heldout/ascending`: 16.6 ms against 8.24 µs |
| `grid/pchip` | construction | SciPy | yes | 1 | 59.5x | `grid/pchip/3d/each40/eval10000/ascending`: 507 µs against 8.53 µs |
| `grid/linear` | construction | scirs2-interpolate | yes | 28 | 1.16x to 4.74x | `grid/linear/3d/max40/eval100/ascending`: 7.52 µs against 1.59 µs |
| `linear_nd` | construction | spade | yes | 16 | 1.58x to 3.32x | `linear_nd/n10/image50`: 5.7 µs against 1.72 µs |
| `clough_tocher` | construction | scirs2-interpolate | no | 2 | 1.37x to 2.9x | `clough_tocher/n10/image50`: 9.68 µs against 3.34 µs |
| `linear_nd` | construction | scirs2-interpolate | yes | 3 | 2.83x to 2.85x | `linear_nd/n10/image500`: 5.75 µs against 2.01 µs |
| `oned/linear` | construction | ndarray-interp | yes | 11 | 1.39x to 2.41x | `oned/linear/n10/sorted`: 101 ns against 41.8 ns |
| `grid/nearest` | construction | scirs2-interpolate | yes | 3 | 1.31x to 2.36x | `grid/nearest/3d/each40/eval10000/ascending`: 21 µs against 8.91 µs |
| `oned/akima` | construction | scirs2-interpolate | yes | 8 | 1.45x to 2.21x | `oned/akima/n10000/shuffled`: 106 µs against 47.9 µs |
| `oned/pchip` | construction | scirs2-interpolate | yes | 12 | 1.37x to 2.13x | `oned/pchip/n10000/sorted-uneven`: 102 µs against 48 µs |
| `oned/linear` | construction | scirs2-interpolate | yes | 12 | 1.14x to 1.77x | `oned/linear/n10000/sorted`: 22.7 µs against 12.8 µs |
| `grid/linear` | evaluation | interpn | yes | 17 | 1.11x to 1.49x | `grid/linear/3d/max10/eval100/ascending`: 2.09 µs against 1.4 µs |
| `oned/cubic_spline` | construction | scirs2-interpolate | no | 2 | 1.14x to 1.15x | `oned/cubic_spline/n10000/sorted`: 167 µs against 146 µs |
| `grid/nearest` | evaluation | interpn | yes | 1 | 1.11x | `grid/nearest/2d/each40/eval10000/ascending`: 164 µs against 148 µs |

<details><summary>Every comparison where interpx is slower</summary>

"Values apart" is the distance of that library's values from SciPy's, as a share of the largest value.

| case | operation | against | interpx | the other | slower by | values apart |
| :--- | :--- | :--- | ---: | ---: | ---: | ---: |
| `grid/cubic/2d/auckland-ba32-0401/heldout/ascending` | construction | interpn | 82.8 ms | 0.88 ns | 94164650.67x | 5e-03 |
| `grid/linear/4d/max200/eval10/ascending` | construction | interpn | 51.8 ms | 1.27 ns | 40835229.54x | 2e-16 |
| `grid/linear/4d/max200/eval100/ascending` | construction | interpn | 50.3 ms | 1.27 ns | 39631399.38x | 3e-16 |
| `grid/linear/4d/max200/eval10000/ascending` | construction | interpn | 48 ms | 1.27 ns | 37803401.85x | 4e-16 |
| `grid/linear/4d/max200/eval1000/ascending` | construction | interpn | 48.3 ms | 1.37 ns | 35333485.71x | 3e-16 |
| `grid/nearest/2d/auckland-ba32-0401/heldout/ascending` | construction | interpn | 16.6 ms | 0.68 ns | 24286354.29x | 0e+00 |
| `grid/linear/2d/auckland-ba32-0401/heldout/ascending` | construction | interpn | 16.6 ms | 0.68 ns | 24250514.29x | 0e+00 |
| `grid/linear/3d/max200/eval1000/ascending` | construction | interpn | 2.04 ms | 0.98 ns | 2092646.40x | 3e-16 |
| `grid/linear/3d/max200/eval10000/ascending` | construction | interpn | 1.81 ms | 0.98 ns | 1848524.80x | 5e-16 |
| `grid/linear/3d/max200/eval10/ascending` | construction | interpn | 1.69 ms | 0.98 ns | 1735321.60x | 3e-16 |
| `grid/cubic/2d/each400/eval1000/ascending` | construction | interpn | 1.46 ms | 0.88 ns | 1666332.44x | 5e-01 |
| `grid/linear/3d/max200/eval100/ascending` | construction | interpn | 1.74 ms | 1.07 ns | 1620852.36x | 3e-16 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | interpn | 1.08 ms | 0.68 ns | 1582884.57x | 3e-16 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | interpn | 332 µs | 0.78 ns | 424320.00x | 3e-16 |
| `grid/cubic/2d/each200/eval1000/ascending` | construction | interpn | 193 µs | 0.88 ns | 219363.56x | 2e-01 |
| `grid/cubic/3d/each40/eval10000/ascending` | construction | interpn | 372 µs | 1.76 ns | 211740.44x | 5e-01 |
| `grid/cubic/2d/each100/eval1000/ascending` | construction | interpn | 44.8 µs | 0.98 ns | 45894.40x | 2e-01 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | interpn | 21 µs | 0.88 ns | 23907.56x | 0e+00 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | interpn | 7.44 µs | 0.68 ns | 10889.14x | 2e-16 |
| `grid/linear/2d/max200/eval10/ascending` | construction | interpn | 7.43 µs | 0.68 ns | 10870.86x | 2e-16 |
| `grid/linear/4d/max40/eval10/ascending` | construction | interpn | 14.4 µs | 1.37 ns | 10514.29x | 2e-16 |
| `grid/linear/4d/max40/eval10000/ascending` | construction | interpn | 14.3 µs | 1.37 ns | 10482.29x | 5e-16 |
| `grid/cubic/2d/each40/eval10000/ascending` | construction | interpn | 9.09 µs | 0.88 ns | 10346.67x | 7e-01 |
| `grid/linear/2d/max200/eval100/ascending` | construction | interpn | 7.45 µs | 0.78 ns | 9536.00x | 1e-16 |
| `grid/linear/4d/max40/eval100/ascending` | construction | interpn | 13 µs | 1.37 ns | 9536.00x | 2e-16 |
| `grid/cubic/2d/each40/eval1000/ascending` | construction | interpn | 9.1 µs | 0.98 ns | 9321.60x | 7e-01 |
| `grid/linear/2d/max200/eval10000/ascending` | construction | interpn | 7.05 µs | 0.78 ns | 9028.00x | 3e-16 |
| `grid/linear/3d/max40/eval100/ascending` | construction | interpn | 7.52 µs | 0.98 ns | 7699.20x | 2e-16 |
| `grid/linear/4d/max200/eval10/ascending` | construction | SciPy | 51.8 ms | 9.36 µs | 5537.15x |  |
| `grid/linear/4d/max40/eval1000/ascending` | construction | interpn | 14.3 µs | 2.64 ns | 5430.52x | 6e-16 |
| `grid/linear/4d/max200/eval10/flipped` | construction | SciPy | 109 ms | 21.7 µs | 5053.67x |  |
| `grid/linear/4d/max200/eval100/flipped` | construction | SciPy | 109 ms | 21.8 µs | 5015.30x |  |
| `grid/linear/4d/max200/eval10000/ascending` | construction | SciPy | 48 ms | 9.65 µs | 4971.72x |  |
| `grid/linear/4d/max200/eval10000/flipped` | construction | SciPy | 109 ms | 21.9 µs | 4966.35x |  |
| `grid/linear/4d/max200/eval100/ascending` | construction | SciPy | 50.3 ms | 10.1 µs | 4958.50x |  |
| `grid/linear/4d/max200/eval1000/flipped` | construction | SciPy | 104 ms | 21.8 µs | 4790.02x |  |
| `grid/linear/4d/max200/eval1000/ascending` | construction | SciPy | 48.3 ms | 10.3 µs | 4704.32x |  |
| `grid/linear/3d/max40/eval10000/ascending` | construction | interpn | 2.99 µs | 0.88 ns | 3399.11x | 3e-16 |
| `grid/linear/3d/max40/eval1000/ascending` | construction | interpn | 2.91 µs | 0.98 ns | 2976.80x | 2e-16 |
| `grid/linear/3d/max40/eval10/ascending` | construction | interpn | 2.69 µs | 0.98 ns | 2755.20x | 2e-16 |
| `grid/nearest/2d/auckland-ba32-0401/heldout/ascending` | construction | SciPy | 16.6 ms | 6.94 µs | 2393.08x |  |
| `grid/linear/2d/auckland-ba32-0401/heldout/ascending` | construction | SciPy | 16.6 ms | 7.19 µs | 2305.43x |  |
| `grid/linear/2d/auckland-ba32-0401/heldout/ascending` | construction | ndarray-interp | 16.6 ms | 8.24 µs | 2012.44x | 0e+00 |
| `grid/cubic/2d/each10/eval1000/ascending` | construction | interpn | 1.59 µs | 1.76 ns | 906.22x | 3e-01 |
| `grid/nearest/2d/each40/eval10000/ascending` | construction | interpn | 627 ns | 0.78 ns | 802.50x | 0e+00 |
| `grid/linear/2d/max40/eval10000/ascending` | construction | interpn | 392 ns | 0.68 ns | 574.14x | 4e-16 |
| `grid/linear/2d/max40/eval100/ascending` | construction | interpn | 392 ns | 0.68 ns | 573.00x | 2e-16 |
| `grid/linear/2d/max40/eval10/ascending` | construction | interpn | 393 ns | 0.78 ns | 502.50x | 2e-16 |
| `grid/linear/2d/max40/eval1000/ascending` | construction | interpn | 362 ns | 0.78 ns | 463.06x | 3e-16 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | ndarray-interp | 1.08 ms | 2.44 µs | 443.07x | 6e-16 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | ndarray-interp | 332 µs | 1.25 µs | 264.91x | 5e-16 |
| `grid/linear/3d/max200/eval1000/ascending` | construction | SciPy | 2.04 ms | 8.27 µs | 247.01x |  |
| `grid/linear/3d/max200/eval10000/ascending` | construction | SciPy | 1.81 ms | 7.63 µs | 236.65x |  |
| `grid/linear/3d/max200/eval10/ascending` | construction | SciPy | 1.69 ms | 7.62 µs | 222.52x |  |
| `grid/linear/3d/max10/eval100/ascending` | construction | interpn | 194 ns | 0.88 ns | 221.00x | 3e-16 |
| `grid/linear/3d/max200/eval100/ascending` | construction | SciPy | 1.74 ms | 8.24 µs | 211.37x |  |
| `grid/linear/3d/max10/eval10000/ascending` | construction | interpn | 179 ns | 0.88 ns | 203.89x | 4e-16 |
| `grid/linear/3d/max10/eval1000/ascending` | construction | interpn | 178 ns | 0.88 ns | 203.00x | 5e-16 |
| `grid/linear/3d/max200/eval100/flipped` | construction | SciPy | 3.5 ms | 17.6 µs | 198.48x |  |
| `grid/linear/2d/max10/eval10000/ascending` | construction | interpn | 132 ns | 0.68 ns | 193.43x | 3e-16 |
| `grid/linear/3d/max200/eval1000/flipped` | construction | SciPy | 3.32 ms | 17.2 µs | 193.07x |  |
| `grid/linear/3d/max10/eval10/ascending` | construction | interpn | 179 ns | 0.98 ns | 183.10x | 2e-16 |
| `grid/linear/3d/max200/eval10/flipped` | construction | SciPy | 3.18 ms | 17.5 µs | 181.74x |  |
| `grid/linear/3d/max200/eval10000/flipped` | construction | SciPy | 3.05 ms | 17.1 µs | 177.88x |  |
| `grid/linear/2d/max10/eval10/ascending` | construction | interpn | 135 ns | 0.78 ns | 173.00x | 1e-16 |
| `grid/linear/2d/max10/eval100/ascending` | construction | interpn | 122 ns | 0.78 ns | 156.75x | 2e-16 |
| `grid/linear/2d/max10/eval1000/ascending` | construction | interpn | 136 ns | 0.88 ns | 155.22x | 3e-16 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | SciPy | 1.08 ms | 7.03 µs | 153.82x |  |
| `grid/pchip/3d/each40/eval10000/ascending` | construction | SciPy | 507 µs | 8.53 µs | 59.49x |  |
| `grid/linear/2d/each400/eval10000/ascending` | construction | SciPy | 332 µs | 6.47 µs | 51.27x |  |
| `grid/linear/2d/max200/eval100/ascending` | construction | ndarray-interp | 7.45 µs | 515 ns | 14.46x | 3e-16 |
| `grid/linear/2d/max200/eval10/ascending` | construction | ndarray-interp | 7.43 µs | 515 ns | 14.43x | 2e-16 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | ndarray-interp | 7.44 µs | 516 ns | 14.43x | 4e-16 |
| `grid/linear/2d/max200/eval10000/ascending` | construction | ndarray-interp | 7.05 µs | 516 ns | 13.67x | 5e-16 |
| `grid/linear/3d/max40/eval100/ascending` | construction | scirs2-interpolate | 7.52 µs | 1.59 µs | 4.74x | 1e-16 |
| `linear_nd/n10/image50` | construction | spade | 5.7 µs | 1.72 µs | 3.32x | 6e-16 |
| `linear_nd/n10/image500` | construction | spade | 5.75 µs | 1.73 µs | 3.31x | 9e-16 |
| `linear_nd/n10/image100` | construction | spade | 5.74 µs | 1.74 µs | 3.30x | 8e-16 |
| `linear_nd/n10000/image100` | construction | spade | 11.5 ms | 3.82 ms | 3.01x | 8e-16 |
| `clough_tocher/n10/image50` | construction | scirs2-interpolate | 9.68 µs | 3.34 µs | 2.90x | 1e-01 |
| `grid/linear/4d/max40/eval10/flipped` | construction | SciPy | 56.8 µs | 19.7 µs | 2.89x |  |
| `grid/linear/4d/max40/eval10000/flipped` | construction | SciPy | 56.4 µs | 19.7 µs | 2.86x |  |
| `linear_nd/n10/image500` | construction | scirs2-interpolate | 5.75 µs | 2.01 µs | 2.85x | 9e-16 |
| `linear_nd/n10/image50` | construction | scirs2-interpolate | 5.7 µs | 2 µs | 2.85x | 8e-16 |
| `grid/linear/4d/max40/eval100/flipped` | construction | SciPy | 56.6 µs | 19.9 µs | 2.84x |  |
| `linear_nd/n10/image100` | construction | scirs2-interpolate | 5.74 µs | 2.03 µs | 2.83x | 9e-16 |
| `grid/linear/4d/max40/eval1000/flipped` | construction | SciPy | 56.9 µs | 21.8 µs | 2.61x |  |
| `grid/linear/2d/max40/eval10000/ascending` | construction | ndarray-interp | 392 ns | 152 ns | 2.58x | 4e-16 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | SciPy | 21 µs | 8.17 µs | 2.57x |  |
| `grid/linear/2d/max40/eval100/ascending` | construction | ndarray-interp | 392 ns | 153 ns | 2.56x | 2e-16 |
| `oned/linear/n10/sorted` | construction | ndarray-interp | 101 ns | 41.8 ns | 2.41x | 0e+00 |
| `grid/linear/4d/max40/eval10/ascending` | construction | scirs2-interpolate | 14.4 µs | 5.97 µs | 2.41x | 2e-16 |
| `grid/linear/2d/max200/eval10/ascending` | construction | scirs2-interpolate | 7.43 µs | 3.1 µs | 2.40x | 1e-16 |
| `grid/linear/4d/max40/eval10000/ascending` | construction | scirs2-interpolate | 14.3 µs | 5.99 µs | 2.39x | 4e-16 |
| `grid/linear/2d/max40/eval10/ascending` | construction | ndarray-interp | 393 ns | 166 ns | 2.37x | 2e-16 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | scirs2-interpolate | 21 µs | 8.91 µs | 2.36x | 0e+00 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | scirs2-interpolate | 7.44 µs | 3.24 µs | 2.30x | 2e-16 |
| `grid/linear/2d/max200/eval100/ascending` | construction | scirs2-interpolate | 7.45 µs | 3.25 µs | 2.29x | 2e-16 |
| `oned/linear/n10/shuffled` | construction | ndarray-interp | 101 ns | 45.5 ns | 2.22x | 0e+00 |
| `oned/akima/n10000/shuffled` | construction | scirs2-interpolate | 106 µs | 47.9 µs | 2.21x | 1e-16 |
| `grid/linear/4d/max40/eval100/ascending` | construction | scirs2-interpolate | 13 µs | 5.97 µs | 2.19x | 3e-16 |
| `grid/linear/2d/max200/eval10000/ascending` | construction | scirs2-interpolate | 7.05 µs | 3.24 µs | 2.18x | 2e-16 |
| `grid/nearest/2d/each40/eval10000/ascending` | construction | scirs2-interpolate | 627 ns | 289 ns | 2.17x | 0e+00 |
| `grid/linear/2d/max40/eval1000/ascending` | construction | ndarray-interp | 362 ns | 167 ns | 2.17x | 4e-16 |
| `grid/linear/4d/max40/eval1000/ascending` | construction | scirs2-interpolate | 14.3 µs | 6.61 µs | 2.17x | 2e-16 |
| `oned/pchip/n10000/sorted-uneven` | construction | scirs2-interpolate | 102 µs | 48 µs | 2.13x | 3e-16 |
| `linear_nd/n10000/image50` | construction | spade | 8.42 ms | 3.96 ms | 2.13x | 8e-16 |
| `linear_nd/n10000/image500` | construction | spade | 8 ms | 3.78 ms | 2.12x | 3e-15 |
| `linear_nd/n1000/image50` | construction | spade | 638 µs | 304 µs | 2.10x | 1e-15 |
| `linear_nd/n1000/image500` | construction | spade | 585 µs | 280 µs | 2.09x | 3e-15 |
| `oned/pchip/n10000/shuffled` | construction | scirs2-interpolate | 94.3 µs | 46.7 µs | 2.02x | 2e-16 |
| `oned/akima/n10000/sorted` | construction | scirs2-interpolate | 88.5 µs | 43.9 µs | 2.02x | 1e-16 |
| `oned/pchip/n10000/shuffled-uneven` | construction | scirs2-interpolate | 93.1 µs | 46.8 µs | 1.99x | 3e-16 |
| `oned/akima/n1000/shuffled` | construction | scirs2-interpolate | 8.73 µs | 4.41 µs | 1.98x | 3e-17 |
| `grid/linear/3d/max40/eval10000/ascending` | construction | scirs2-interpolate | 2.99 µs | 1.51 µs | 1.97x | 2e-16 |
| `oned/pchip/n1000/shuffled` | construction | scirs2-interpolate | 8.96 µs | 4.56 µs | 1.97x | 2e-16 |
| `oned/pchip/n1000/sorted` | construction | scirs2-interpolate | 8.98 µs | 4.58 µs | 1.96x | 2e-16 |
| `linear_nd/n1000/image100` | construction | spade | 593 µs | 304 µs | 1.95x | 1e-15 |
| `grid/linear/3d/max40/eval1000/ascending` | construction | scirs2-interpolate | 2.91 µs | 1.49 µs | 1.95x | 2e-16 |
| `oned/pchip/n10000/sorted` | construction | scirs2-interpolate | 89.9 µs | 47.8 µs | 1.88x | 2e-16 |
| `oned/akima/n1000/sorted` | construction | scirs2-interpolate | 8.78 µs | 4.71 µs | 1.86x | 3e-17 |
| `linear_nd/n100000/image50` | construction | spade | 87 ms | 48.4 ms | 1.80x | 8e-16 |
| `oned/pchip/n100/sorted-uneven` | construction | scirs2-interpolate | 979 ns | 545 ns | 1.80x | 2e-16 |
| `grid/linear/2d/max40/eval10/ascending` | construction | scirs2-interpolate | 393 ns | 219 ns | 1.79x | 4e-17 |
| `linear_nd/n100/image100` | construction | spade | 36.7 µs | 20.6 µs | 1.78x | 9e-16 |
| `grid/linear/2d/max40/eval100/ascending` | construction | scirs2-interpolate | 392 ns | 221 ns | 1.77x | 3e-16 |
| `linear_nd/n100/image50` | construction | spade | 36.5 µs | 20.6 µs | 1.77x | 1e-15 |
| `oned/linear/n10000/sorted` | construction | scirs2-interpolate | 22.7 µs | 12.8 µs | 1.77x | 0e+00 |
| `linear_nd/n100/image500` | construction | spade | 36.6 µs | 20.7 µs | 1.77x | 1e-15 |
| `grid/linear/2d/max40/eval10000/ascending` | construction | scirs2-interpolate | 392 ns | 222 ns | 1.77x | 2e-16 |
| `oned/linear/n10000/sorted-uneven` | construction | scirs2-interpolate | 22.5 µs | 12.8 µs | 1.75x | 6e-17 |
| `oned/linear/n10000/shuffled` | construction | scirs2-interpolate | 21.1 µs | 12.3 µs | 1.72x | 0e+00 |
| `grid/linear/3d/max40/eval10/ascending` | construction | scirs2-interpolate | 2.69 µs | 1.58 µs | 1.70x | 1e-16 |
| `linear_nd/n100000/image100` | construction | spade | 78.9 ms | 47 ms | 1.68x | 7e-16 |
| `linear_nd/n100000/image500` | construction | spade | 78.9 ms | 47.3 ms | 1.67x | 6e-15 |
| `oned/linear/n50/sorted` | construction | ndarray-interp | 167 ns | 100 ns | 1.67x | 0e+00 |
| `oned/linear/n1000/shuffled` | construction | scirs2-interpolate | 2.11 µs | 1.27 µs | 1.66x | 0e+00 |
| `grid/linear/2d/max40/eval1000/ascending` | construction | scirs2-interpolate | 362 ns | 219 ns | 1.65x | 3e-16 |
| `oned/linear/n10000/shuffled-uneven` | construction | scirs2-interpolate | 21.1 µs | 12.9 µs | 1.64x | 6e-17 |
| `oned/akima/n100/shuffled` | construction | scirs2-interpolate | 955 ns | 586 ns | 1.63x | 1e-16 |
| `oned/pchip/n100/shuffled-uneven` | construction | scirs2-interpolate | 958 ns | 590 ns | 1.62x | 2e-16 |
| `oned/linear/n100/sorted-uneven` | construction | ndarray-interp | 270 ns | 167 ns | 1.62x | 0e+00 |
| `oned/pchip/n100/sorted` | construction | scirs2-interpolate | 960 ns | 596 ns | 1.61x | 2e-16 |
| `oned/linear/n1000/sorted` | construction | scirs2-interpolate | 2.12 µs | 1.32 µs | 1.61x | 0e+00 |
| `oned/pchip/n100/shuffled` | construction | scirs2-interpolate | 959 ns | 597 ns | 1.61x | 2e-16 |
| `grid/linear/2d/max10/eval1000/ascending` | construction | ndarray-interp | 136 ns | 85.6 ns | 1.59x | 3e-16 |
| `oned/akima/n100/sorted` | construction | scirs2-interpolate | 953 ns | 599 ns | 1.59x | 1e-16 |
| `grid/linear/2d/max10/eval10/ascending` | construction | ndarray-interp | 135 ns | 85.4 ns | 1.58x | 1e-16 |
| `linear_nd/ghcn-tmax-20250101/heldout` | construction | spade | 7.96 ms | 5.04 ms | 1.58x | 7e-16 |
| `grid/linear/2d/max10/eval100/ascending` | construction | ndarray-interp | 122 ns | 77.9 ns | 1.57x | 4e-16 |
| `grid/linear/2d/max10/eval10000/ascending` | construction | ndarray-interp | 132 ns | 85.5 ns | 1.55x | 4e-16 |
| `oned/linear/n50/shuffled` | construction | ndarray-interp | 168 ns | 109 ns | 1.54x | 0e+00 |
| `grid/linear/4d/max40/eval10/ascending` | construction | SciPy | 14.4 µs | 9.38 µs | 1.53x |  |
| `oned/akima/n50/shuffled` | construction | scirs2-interpolate | 548 ns | 364 ns | 1.51x | 1e-16 |
| `oned/pchip/n50/sorted` | construction | scirs2-interpolate | 490 ns | 326 ns | 1.50x | 1e-16 |
| `grid/linear/3d/max10/eval100/ascending` | evaluation | interpn | 2.09 µs | 1.4 µs | 1.49x | 3e-16 |
| `oned/linear/n100/shuffled-uneven` | construction | ndarray-interp | 247 ns | 167 ns | 1.48x | 0e+00 |
| `grid/linear/3d/max200/eval1000/ascending` | construction | scirs2-interpolate | 2.04 ms | 1.38 ms | 1.48x | 3e-16 |
| `grid/linear/2d/max10/eval10000/ascending` | evaluation | interpn | 137 µs | 92.6 µs | 1.47x | 3e-16 |
| `oned/akima/n50/sorted` | construction | scirs2-interpolate | 536 ns | 368 ns | 1.45x | 1e-16 |
| `oned/linear/n100/shuffled` | construction | ndarray-interp | 246 ns | 169 ns | 1.45x | 0e+00 |
| `grid/linear/3d/max10/eval10000/ascending` | evaluation | interpn | 205 µs | 142 µs | 1.45x | 4e-16 |
| `oned/linear/n10000/shuffled` | construction | ndarray-interp | 21.1 µs | 14.8 µs | 1.42x | 0e+00 |
| `grid/linear/4d/max40/eval1000/ascending` | construction | SciPy | 14.3 µs | 10.1 µs | 1.41x |  |
| `grid/linear/2d/max10/eval10/ascending` | evaluation | interpn | 147 ns | 104 ns | 1.41x | 1e-16 |
| `oned/linear/n100/sorted` | construction | ndarray-interp | 258 ns | 183 ns | 1.41x | 0e+00 |
| `grid/linear/4d/max40/eval10000/ascending` | construction | SciPy | 14.3 µs | 10.2 µs | 1.41x |  |
| `oned/linear/n1000/shuffled` | construction | ndarray-interp | 2.11 µs | 1.52 µs | 1.39x | 0e+00 |
| `oned/linear/n1000/sorted` | construction | ndarray-interp | 2.12 µs | 1.53 µs | 1.39x | 0e+00 |
| `grid/linear/2d/max200/eval10000/flipped` | construction | SciPy | 16.2 µs | 11.8 µs | 1.37x |  |
| `oned/pchip/n50/shuffled` | construction | scirs2-interpolate | 486 ns | 355 ns | 1.37x | 1e-16 |
| `clough_tocher/n50/image50` | construction | scirs2-interpolate | 39.1 µs | 28.6 µs | 1.37x | 1e-01 |
| `grid/linear/2d/max10/eval1000/ascending` | evaluation | interpn | 13.7 µs | 10.1 µs | 1.36x | 3e-16 |
| `grid/linear/3d/max10/eval1000/ascending` | evaluation | interpn | 19.1 µs | 14.1 µs | 1.35x | 5e-16 |
| `oned/linear/n100/sorted-uneven` | construction | scirs2-interpolate | 270 ns | 203 ns | 1.33x | 3e-17 |
| `grid/linear/4d/max200/eval10/ascending` | construction | scirs2-interpolate | 51.8 ms | 39.2 ms | 1.32x | 2e-16 |
| `oned/linear/n100/sorted` | construction | scirs2-interpolate | 258 ns | 196 ns | 1.31x | 0e+00 |
| `grid/nearest/2d/auckland-ba32-0401/heldout/ascending` | construction | scirs2-interpolate | 16.6 ms | 12.7 ms | 1.31x | 0e+00 |
| `grid/linear/3d/max200/eval10/ascending` | construction | scirs2-interpolate | 1.69 ms | 1.31 ms | 1.29x | 3e-16 |
| `grid/linear/2d/auckland-ba32-0401/heldout/ascending` | construction | scirs2-interpolate | 16.6 ms | 12.9 ms | 1.29x | 0e+00 |
| `grid/linear/3d/max10/eval10/ascending` | evaluation | interpn | 200 ns | 155 ns | 1.29x | 2e-16 |
| `grid/linear/4d/max40/eval100/ascending` | construction | SciPy | 13 µs | 10.2 µs | 1.27x |  |
| `grid/linear/2d/each800/eval10000/ascending` | construction | scirs2-interpolate | 1.08 ms | 859 µs | 1.26x | 2e-16 |
| `oned/linear/n100/shuffled-uneven` | construction | scirs2-interpolate | 247 ns | 197 ns | 1.26x | 3e-17 |
| `oned/linear/n50/shuffled` | construction | scirs2-interpolate | 168 ns | 134 ns | 1.25x | 0e+00 |
| `oned/linear/n50/sorted` | construction | scirs2-interpolate | 167 ns | 134 ns | 1.25x | 0e+00 |
| `grid/linear/4d/max200/eval100/ascending` | construction | scirs2-interpolate | 50.3 ms | 40.3 ms | 1.25x | 1e-16 |
| `grid/linear/2d/max200/eval100/flipped` | construction | SciPy | 15.8 µs | 12.9 µs | 1.23x |  |
| `grid/linear/2d/max200/eval10/flipped` | construction | SciPy | 14.5 µs | 11.9 µs | 1.23x |  |
| `grid/linear/2d/max10/eval100/ascending` | evaluation | interpn | 1.25 µs | 1.02 µs | 1.22x | 2e-16 |
| `grid/linear/3d/max200/eval100/ascending` | construction | scirs2-interpolate | 1.74 ms | 1.44 ms | 1.21x | 3e-16 |
| `grid/linear/4d/max200/eval1000/ascending` | construction | scirs2-interpolate | 48.3 ms | 40 ms | 1.21x | 3e-16 |
| `grid/linear/2d/max200/eval1000/flipped` | construction | SciPy | 14.2 µs | 11.8 µs | 1.21x |  |
| `grid/linear/3d/max200/eval10000/ascending` | construction | scirs2-interpolate | 1.81 ms | 1.51 ms | 1.20x | 3e-16 |
| `grid/linear/2d/max200/eval10000/ascending` | construction | SciPy | 7.05 µs | 5.9 µs | 1.20x |  |
| `grid/linear/4d/max200/eval10000/ascending` | construction | scirs2-interpolate | 48 ms | 40.5 ms | 1.19x | 4e-16 |
| `grid/linear/2d/max40/eval10/ascending` | evaluation | interpn | 172 ns | 147 ns | 1.17x | 2e-16 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | SciPy | 7.44 µs | 6.38 µs | 1.17x |  |
| `grid/linear/3d/max10/eval100/ascending` | construction | scirs2-interpolate | 194 ns | 167 ns | 1.16x | 2e-16 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | scirs2-interpolate | 332 µs | 286 µs | 1.16x | 2e-16 |
| `grid/linear/2d/max200/eval10/ascending` | construction | SciPy | 7.43 µs | 6.42 µs | 1.16x |  |
| `grid/linear/2d/max200/eval100/ascending` | construction | SciPy | 7.45 µs | 6.43 µs | 1.16x |  |
| `grid/linear/2d/max40/eval10000/ascending` | evaluation | interpn | 167 µs | 145 µs | 1.15x | 4e-16 |
| `grid/linear/2d/max40/eval1000/ascending` | evaluation | interpn | 16.6 µs | 14.4 µs | 1.15x | 3e-16 |
| `oned/cubic_spline/n10000/sorted` | construction | scirs2-interpolate | 167 µs | 146 µs | 1.15x | 4e-05 |
| `oned/linear/n100/shuffled` | construction | scirs2-interpolate | 246 ns | 216 ns | 1.14x | 0e+00 |
| `grid/linear/4d/max40/eval10/ascending` | evaluation | interpn | 310 ns | 273 ns | 1.14x | 2e-16 |
| `oned/cubic_spline/n10000/shuffled` | construction | scirs2-interpolate | 167 µs | 147 µs | 1.14x | 4e-05 |
| `grid/linear/2d/max200/eval1000/ascending` | evaluation | interpn | 19.3 µs | 17.2 µs | 1.12x | 2e-16 |
| `grid/linear/4d/max200/eval10000/ascending` | evaluation | interpn | 408 µs | 365 µs | 1.12x | 4e-16 |
| `grid/nearest/2d/each40/eval10000/ascending` | evaluation | interpn | 164 µs | 148 µs | 1.11x | 0e+00 |
| `grid/linear/4d/max40/eval10000/ascending` | evaluation | interpn | 330 µs | 296 µs | 1.11x | 5e-16 |
| `grid/linear/2d/max200/eval10/ascending` | evaluation | interpn | 195 ns | 176 ns | 1.11x | 2e-16 |
| `grid/linear/3d/max40/eval10/ascending` | evaluation | interpn | 222 ns | 201 ns | 1.11x | 2e-16 |

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

A run (one library on one case, a process of its own) is stopped at 15 minutes. The longest run that finished took 369 s, 41% of that: scirs2-interpolate on `clough_tocher/n100000/image500`.

No run took more than half the limit (450 s).

## Real-world data

Two data sets, each taken whole in the unit its source gives it out in; what was taken and how it was divided was fixed before anything was measured (`compare/perf/real.py`, docs/design.ja.md section 8.4.7).

- `auckland-ba32-0401`: Auckland Part 1 LiDAR 1m DEM (2024), tile BA32_10000_0401. 4800 by 7200 cells of 1 m. The grid is every second row and column (and the last ones); every other cell is a query, and its elevation is the value held out. Licensed by Regional Software Holdings Limited for re-use under CC-BY-4.0. Hosted by Toitū Te Whenua Land Information New Zealand (https://registry.opendata.aws/nz-elevation).
- `ghcn-tmax-20250101`: GHCN-Daily, maximum temperature of 2025-01-01, every station; longitude and latitude in degrees, as they are. Of every 10 stations in the order of the file the last is a query, and its temperature is the value held out. The RBF kernels are those without a shape parameter. NOAA Global Historical Climatology Network Daily (GHCN-D), from https://registry.opendata.aws/noaa-ghcn (CC0 1.0).

These cases are also in the tables of their families below and in the counts above. Here they are with the error: the largest and the root mean square distance from the values held out (metres, and degrees Celsius), and the number of queries a library gave no value. A library missing from a case does not have it, or did not finish (see the runs that failed).

`grid/nearest/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 16.6 ms | 362 ms | 15.5 | 0.178 | 0 |
| interpx, parallel | 16.6 ms | 66.3 ms | 15.5 | 0.178 | 0 |
| interpx, `from_owned` | 5.56 ms | 374 ms | 15.5 | 0.178 | 0 |
| SciPy | 6.94 µs | 1.41 s | 15.5 | 0.178 | 0 |
| scirs2-interpolate | 12.7 ms | 157 s | 15.5 | 0.178 | 0 |
| interpn | 0.68 ns | 755 ms | 15.5 | 0.178 | 0 |

`grid/linear/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 16.6 ms | 381 ms | 6.74 | 0.0601 | 0 |
| interpx, parallel | 16.4 ms | 67.6 ms | 6.74 | 0.0601 | 0 |
| interpx, `from_owned` | 5.5 ms | 386 ms | 6.74 | 0.0601 | 0 |
| SciPy | 7.19 µs | 1.14 s | 6.74 | 0.0601 | 0 |
| scirs2-interpolate | 12.9 ms | 71.2 s | 6.74 | 0.0601 | 0 |
| ndarray-interp | 8.24 µs | 365 ms | 6.74 | 0.0601 | 0 |
| interpn | 0.68 ns | 821 ms | 6.74 | 0.0601 | 0 |

`grid/cubic/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 82.8 ms | 779 ms | 6.77 | 0.0507 | 0 |
| interpx, parallel | 78.1 ms | 135 ms | 6.77 | 0.0507 | 0 |
| interpx, `from_owned` | 72.7 ms | 786 ms | 6.77 | 0.0507 | 0 |
| SciPy | 15.1 s | 33 s | 6.77 | 0.0507 | 0 |
| interpn | 0.88 ns | 1.81 s | 6.76 | 0.0517 | 0 |

`grid/quintic/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 100 ms | 1.27 s | 6.77 | 0.0517 | 0 |
| interpx, parallel | 97.6 ms | 218 ms | 6.77 | 0.0517 | 0 |
| interpx, `from_owned` | 87.6 ms | 1.28 s | 6.77 | 0.0517 | 0 |
| SciPy | 21.3 s | 38.2 s | 6.77 | 0.0517 | 0 |

`grid/pchip/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 79 ms | 1.05 s | 6.75 | 0.049 | 0 |
| interpx, parallel | 77.5 ms | 200 ms | 6.75 | 0.049 | 0 |
| interpx, `from_owned` | 65.3 ms | 1.09 s | 6.75 | 0.049 | 0 |

`nearest_nd/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 524 µs | 157 µs | 25 | 3.42 | 0 |
| interpx, parallel | 558 µs | 72.2 µs | 25 | 3.42 | 0 |
| interpx, `from_owned` | 551 µs | 173 µs | 25 | 3.42 | 0 |
| SciPy | 1.88 ms | 563 µs | 25 | 3.42 | 0 |

`linear_nd/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 7.96 ms | 81.4 µs | 22.1 | 2.83 | 1 |
| interpx, parallel | 7.99 ms | 84.2 µs | 22.1 | 2.83 | 1 |
| interpx, `from_owned` | 7.94 ms | 78.5 µs | 22.1 | 2.83 | 1 |
| SciPy | 30.6 ms | 1.18 ms | 22.1 | 2.83 | 1 |
| scirs2-interpolate | 513 ms | 147 ms | 22.1 | 2.83 | 1 |
| spade | 5.04 ms | 1.36 ms | 22.1 | 2.83 | 1 |

`clough_tocher/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 12.9 ms | 127 µs | 699 | 28 | 1 |
| interpx, parallel | 13.3 ms | 139 µs | 699 | 28 | 1 |
| interpx, `from_owned` | 12.8 ms | 136 µs | 699 | 28 | 1 |
| SciPy | 45.5 ms | 1.22 ms | 699 | 28 | 1 |
| scirs2-interpolate | 666 ms | 144 ms | 22.3 | 2.89 | 1 |

`rbf/linear/global/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 17.3 s | 18.2 ms | 20.8 | 2.74 | 0 |
| interpx, parallel | 17 s | 3.44 ms | 20.8 | 2.74 | 0 |
| interpx, `from_owned` | 17.9 s | 18.5 ms | 20.8 | 2.74 | 0 |
| SciPy | 113 s | 59.1 ms | 20.8 | 2.74 | 0 |
| scirs2-interpolate | 299 s | 239 ms | 20.8 | 2.74 | 0 |

`rbf/linear/k50/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 780 µs | 19.5 ms | 20.8 | 2.75 | 0 |
| interpx, parallel | 762 µs | 3.62 ms | 20.8 | 2.75 | 0 |
| interpx, `from_owned` | 668 µs | 19.9 ms | 20.8 | 2.75 | 0 |
| SciPy | 1.9 ms | 45.8 ms | 20.8 | 2.75 | 0 |

`rbf/thin_plate_spline/global/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 18.8 s | 87.2 ms | 102 | 5.67 | 0 |
| interpx, parallel | 18.2 s | 16.4 ms | 102 | 5.67 | 0 |
| interpx, `from_owned` | 20 s | 88.1 ms | 102 | 5.67 | 0 |
| SciPy | 119 s | 123 ms | 102 | 5.67 | 0 |
| scirs2-interpolate | 287 s | 421 ms | 93.7 | 5.41 | 0 |

`rbf/thin_plate_spline/k50/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 759 µs | 31.5 ms | 102 | 5.66 | 0 |
| interpx, parallel | 810 µs | 5.46 ms | 102 | 5.66 | 0 |
| interpx, `from_owned` | 744 µs | 31.8 ms | 102 | 5.66 | 0 |
| SciPy | 1.96 ms | 60.1 ms | 102 | 5.66 | 0 |

`rbf/cubic/global/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 17.8 s | 18.5 ms | 587 | 27.1 | 0 |
| interpx, parallel | 18.4 s | 3.9 ms | 587 | 27.1 | 0 |
| interpx, `from_owned` | 18.9 s | 19.8 ms | 587 | 27.1 | 0 |
| SciPy | 117 s | 53.2 ms | 597 | 24.4 | 0 |
| scirs2-interpolate | 294 s | 262 ms | 100 | 8.48 | 0 |

`rbf/cubic/k50/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 760 µs | 21.4 ms | 567 | 25.3 | 0 |
| interpx, parallel | 836 µs | 4.29 ms | 567 | 25.3 | 0 |
| interpx, `from_owned` | 680 µs | 21.8 ms | 567 | 25.3 | 0 |
| SciPy | 1.89 ms | 51.5 ms | 567 | 25.3 | 0 |

`rbf/quintic/global/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 21.9 s | 20.3 ms | 2.55e+05 | 1.36e+04 | 0 |
| interpx, parallel | 17.1 s | 3.81 ms | 2.55e+05 | 1.36e+04 | 0 |
| interpx, `from_owned` | 18.2 s | 20.5 ms | 2.55e+05 | 1.36e+04 | 0 |
| SciPy | 126 s | 55.3 ms | 9.39e+04 | 4.92e+03 | 0 |
| scirs2-interpolate | 289 s | 296 ms | 6.56e+16 | 2.85e+15 | 0 |

`rbf/quintic/k50/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 785 µs | 22.5 ms | 3.74e+03 | 199 | 0 |
| interpx, parallel | 814 µs | 4.11 ms | 3.74e+03 | 199 | 0 |
| interpx, `from_owned` | 690 µs | 23.1 ms | 3.74e+03 | 199 | 0 |
| SciPy | 1.76 ms | 55 ms | 3.74e+03 | 200 | 0 |

## One variable

Data of `n` points, evaluated at `4 n` points: in order (`sorted`) and in no order (`shuffled`).

### Evaluation

| case | interpx | interpx, parallel | SciPy | scirs2-interpolate | ndarray-interp | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/n10/sorted` | 126 ns | 193 ns | 5.48 µs | 284 ns | 279 ns | 2.2x faster |
| `linear/n10/shuffled` | 171 ns | 245 ns | 5.6 µs | 282 ns | 274 ns | 1.6x faster |
| `nearest/n10/sorted` | 161 ns | 210 ns | 10.3 µs | 296 ns |  | 1.8x faster |
| `nearest/n10/shuffled` | 333 ns | 401 ns | 10.4 µs | 306 ns |  | **1.1x slower** |
| `previous/n10/sorted` | 166 ns | 214 ns | 11.8 µs |  |  | 71x faster |
| `previous/n10/shuffled` | 290 ns | 322 ns | 11.2 µs |  |  | 39x faster |
| `bspline1/n10/sorted` | 232 ns | 291 ns | 2.92 µs |  |  | 13x faster |
| `bspline1/n10/shuffled` | 597 ns | 600 ns | 2.97 µs |  |  | 5x faster |
| `bspline2/n10/sorted` | 298 ns | 344 ns | 3.38 µs |  |  | 11x faster |
| `bspline2/n10/shuffled` | 340 ns | 416 ns | 3.41 µs |  |  | 10x faster |
| `bspline3/n10/sorted` | 348 ns | 435 ns | 4.03 µs |  |  | 12x faster |
| `bspline3/n10/shuffled` | 440 ns | 566 ns | 3.57 µs |  |  | 8.1x faster |
| `cubic_spline/n10/sorted` | 153 ns | 211 ns | 3.68 µs | 300 ns | 364 ns | 2x faster |
| `cubic_spline/n10/shuffled` | 206 ns | 290 ns | 4.05 µs | 286 ns | 356 ns | 1.4x faster |
| `pchip/n10/sorted` | 155 ns | 197 ns | 3.64 µs | 317 ns |  | 2x faster |
| `pchip/n10/shuffled` | 216 ns | 266 ns | 4 µs | 297 ns |  | 1.4x faster |
| `akima/n10/sorted` | 143 ns | 194 ns | 3.65 µs | 275 ns |  | 1.9x faster |
| `akima/n10/shuffled` | 227 ns | 259 ns | 3.69 µs | 287 ns |  | 1.3x faster |
| `makima/n10/sorted` | 158 ns | 212 ns | 3.87 µs |  |  | 25x faster |
| `makima/n10/shuffled` | 225 ns | 283 ns | 3.66 µs |  |  | 16x faster |
| `linear/n50/sorted` | 624 ns | 689 ns | 6.17 µs | 1.62 µs | 1.26 µs | 2x faster |
| `linear/n50/shuffled` | 816 ns | 938 ns | 8.13 µs | 1.8 µs | 1.38 µs | 1.7x faster |
| `nearest/n50/sorted` | 772 ns | 849 ns | 12.3 µs | 1.85 µs |  | 2.4x faster |
| `nearest/n50/shuffled` | 1.58 µs | 1.65 µs | 12 µs | 1.87 µs |  | 1.2x faster |
| `previous/n50/sorted` | 770 ns | 865 ns | 13.9 µs |  |  | 18x faster |
| `previous/n50/shuffled` | 1.5 µs | 1.46 µs | 13.9 µs |  |  | 9.2x faster |
| `bspline1/n50/sorted` | 1.22 µs | 1.38 µs | 4.04 µs |  |  | 3.3x faster |
| `bspline1/n50/shuffled` | 3.09 µs | 3.16 µs | 5.28 µs |  |  | 1.7x faster |
| `bspline2/n50/sorted` | 1.4 µs | 1.75 µs | 6.69 µs |  |  | 4.8x faster |
| `bspline2/n50/shuffled` | 2.01 µs | 2.28 µs | 6.87 µs |  |  | 3.4x faster |
| `bspline3/n50/sorted` | 1.83 µs | 1.9 µs | 9.2 µs |  |  | 5x faster |
| `bspline3/n50/shuffled` | 3.04 µs | 2.86 µs | 10.2 µs |  |  | 3.4x faster |
| `cubic_spline/n50/sorted` | 724 ns | 815 ns | 5.62 µs | 4.3 µs | 1.74 µs | 2.4x faster |
| `cubic_spline/n50/shuffled` | 997 ns | 1.05 µs | 5.59 µs | 4.44 µs | 1.6 µs | 1.6x faster |
| `pchip/n50/sorted` | 688 ns | 812 ns | 5.2 µs | 3.93 µs |  | 5.7x faster |
| `pchip/n50/shuffled` | 984 ns | 1.14 µs | 5.66 µs | 3.79 µs |  | 3.9x faster |
| `akima/n50/sorted` | 753 ns | 757 ns | 5.69 µs | 4.7 µs |  | 6.2x faster |
| `akima/n50/shuffled` | 1.08 µs | 1.04 µs | 5.9 µs | 4.65 µs |  | 4.3x faster |
| `makima/n50/sorted` | 737 ns | 756 ns | 5.21 µs |  |  | 7.1x faster |
| `makima/n50/shuffled` | 1.07 µs | 1.05 µs | 5.62 µs |  |  | 5.3x faster |
| `linear/n100/sorted` | 1.42 µs | 1.38 µs | 7.94 µs | 3.44 µs | 2.7 µs | 1.9x faster |
| `linear/n100/shuffled` | 1.76 µs | 1.69 µs | 12.8 µs | 3.94 µs | 2.5 µs | 1.4x faster |
| `nearest/n100/sorted` | 1.49 µs | 1.68 µs | 13.5 µs | 3.58 µs |  | 2.4x faster |
| `nearest/n100/shuffled` | 3.21 µs | 3.62 µs | 13.7 µs | 4.13 µs |  | 1.3x faster |
| `previous/n100/sorted` | 1.55 µs | 1.75 µs | 15 µs |  |  | 9.7x faster |
| `previous/n100/shuffled` | 3.01 µs | 2.9 µs | 15.4 µs |  |  | 5.1x faster |
| `bspline1/n100/sorted` | 2.14 µs | 2.51 µs | 5.4 µs |  |  | 2.5x faster |
| `bspline1/n100/shuffled` | 5.94 µs | 6.36 µs | 13.2 µs |  |  | 2.2x faster |
| `bspline2/n100/sorted` | 3 µs | 3.44 µs | 10 µs |  |  | 3.3x faster |
| `bspline2/n100/shuffled` | 4.85 µs | 4.61 µs | 13.9 µs |  |  | 2.9x faster |
| `bspline3/n100/sorted` | 3.98 µs | 4.19 µs | 16.5 µs |  |  | 4.1x faster |
| `bspline3/n100/shuffled` | 6.12 µs | 6.38 µs | 18.1 µs |  |  | 3x faster |
| `cubic_spline/n100/sorted` | 1.54 µs | 1.65 µs | 7.53 µs | 18.4 µs | 3.5 µs | 2.3x faster |
| `cubic_spline/n100/shuffled` | 2.16 µs | 2.26 µs | 8.7 µs | 17.9 µs | 3.52 µs | 1.6x faster |
| `pchip/n100/sorted` | 1.57 µs | 1.69 µs | 7.36 µs | 17.7 µs |  | 4.7x faster |
| `pchip/n100/shuffled` | 2.16 µs | 2.27 µs | 8.84 µs | 17.6 µs |  | 4.1x faster |
| `akima/n100/sorted` | 1.57 µs | 1.68 µs | 7.49 µs | 18.4 µs |  | 4.8x faster |
| `akima/n100/shuffled` | 2.16 µs | 2.27 µs | 8.83 µs | 17.7 µs |  | 4.1x faster |
| `makima/n100/sorted` | 1.56 µs | 1.68 µs | 7.37 µs |  |  | 4.7x faster |
| `makima/n100/shuffled` | 2 µs | 2.26 µs | 7.91 µs |  |  | 4x faster |
| `linear/n1000/sorted` | 14.4 µs | 14.6 µs | 33.3 µs | 136 µs | 27.1 µs | 1.9x faster |
| `linear/n1000/shuffled` | 17.5 µs | 17.7 µs | 161 µs | 242 µs | 27.1 µs | 1.5x faster |
| `nearest/n1000/sorted` | 14.8 µs | 15.9 µs | 63.5 µs | 168 µs |  | 4.3x faster |
| `nearest/n1000/shuffled` | 35.2 µs | 35.4 µs | 218 µs | 253 µs |  | 6.2x faster |
| `previous/n1000/sorted` | 15.4 µs | 16.3 µs | 64.4 µs |  |  | 4.2x faster |
| `previous/n1000/shuffled` | 30.3 µs | 30.5 µs | 224 µs |  |  | 7.4x faster |
| `bspline1/n1000/sorted` | 24.2 µs | 26.8 µs | 34.1 µs |  |  | 1.4x faster |
| `bspline1/n1000/shuffled` | 62.5 µs | 63 µs | 578 µs |  |  | 9.2x faster |
| `bspline2/n1000/sorted` | 34.1 µs | 34.7 µs | 76.7 µs |  |  | 2.2x faster |
| `bspline2/n1000/shuffled` | 63.3 µs | 63.4 µs | 821 µs |  |  | 13x faster |
| `bspline3/n1000/sorted` | 35.9 µs | 39.5 µs | 140 µs |  |  | 3.9x faster |
| `bspline3/n1000/shuffled` | 78.1 µs | 85.8 µs | 633 µs |  |  | 8.1x faster |
| `cubic_spline/n1000/sorted` | 14.6 µs | 15.9 µs | 37.5 µs | 1.56 ms | 32.9 µs | 2.3x faster |
| `cubic_spline/n1000/shuffled` | 21.7 µs | 22.8 µs | 230 µs | 1.56 ms | 35 µs | 1.6x faster |
| `pchip/n1000/sorted` | 16.2 µs | 15.8 µs | 85.6 µs | 1.32 ms |  | 5.3x faster |
| `pchip/n1000/shuffled` | 21.8 µs | 21.9 µs | 238 µs | 1.32 ms |  | 11x faster |
| `akima/n1000/sorted` | 14 µs | 15.8 µs | 41 µs | 1.55 ms |  | 2.9x faster |
| `akima/n1000/shuffled` | 21.7 µs | 21.9 µs | 236 µs | 1.54 ms |  | 11x faster |
| `makima/n1000/sorted` | 15.1 µs | 15.6 µs | 41.7 µs |  |  | 2.8x faster |
| `makima/n1000/shuffled` | 21.7 µs | 33.5 µs | 241 µs |  |  | 11x faster |
| `linear/n10000/sorted` | 149 µs | 73 µs | 241 µs | 1.77 ms | 269 µs | 1.6x faster |
| `linear/n10000/shuffled` | 206 µs | 90.2 µs | 2.5 ms | 2.87 ms | 285 µs | 1.4x faster |
| `nearest/n10000/sorted` | 163 µs | 80.2 µs | 626 µs | 2.1 ms |  | 3.8x faster |
| `nearest/n10000/shuffled` | 428 µs | 114 µs | 3.09 ms | 3.31 ms |  | 7.2x faster |
| `previous/n10000/sorted` | 162 µs | 81.9 µs | 674 µs |  |  | 4.2x faster |
| `previous/n10000/shuffled` | 382 µs | 92.1 µs | 3.23 ms |  |  | 8.5x faster |
| `bspline1/n10000/sorted` | 268 µs | 98.8 µs | 315 µs |  |  | 1.2x faster |
| `bspline1/n10000/shuffled` | 711 µs | 158 µs | 49.9 ms |  |  | 70x faster |
| `bspline2/n10000/sorted` | 353 µs | 101 µs | 812 µs |  |  | 2.3x faster |
| `bspline2/n10000/shuffled` | 978 µs | 225 µs | 51 ms |  |  | 52x faster |
| `bspline3/n10000/sorted` | 400 µs | 114 µs | 1.38 ms |  |  | 3.4x faster |
| `bspline3/n10000/shuffled` | 1.22 ms | 283 µs | 52.6 ms |  |  | 43x faster |
| `cubic_spline/n10000/sorted` | 160 µs | 80.3 µs | 397 µs | 154 ms | 350 µs | 2.2x faster |
| `cubic_spline/n10000/shuffled` | 235 µs | 80.1 µs | 3.67 ms | 158 ms | 366 µs | 1.6x faster |
| `pchip/n10000/sorted` | 154 µs | 70.7 µs | 382 µs | 128 ms |  | 2.5x faster |
| `pchip/n10000/shuffled` | 215 µs | 97.7 µs | 3.5 ms | 119 ms |  | 16x faster |
| `akima/n10000/sorted` | 148 µs | 74.6 µs | 365 µs | 155 ms |  | 2.5x faster |
| `akima/n10000/shuffled` | 215 µs | 95.6 µs | 3.44 ms | 143 ms |  | 16x faster |
| `makima/n10000/sorted` | 177 µs | 69.1 µs | 398 µs |  |  | 2.2x faster |
| `makima/n10000/shuffled` | 235 µs | 92.5 µs | 3.73 ms |  |  | 16x faster |
| `linear/n100/sorted-uneven` | 1.23 µs | 1.42 µs | 8.33 µs | 3.42 µs | 3.87 µs | 2.8x faster |
| `linear/n100/shuffled-uneven` | 2.68 µs | 2.94 µs | 13.1 µs | 3.66 µs | 4.04 µs | 1.4x faster |
| `cubic_spline/n100/sorted-uneven` | 1.38 µs | 1.51 µs | 7.32 µs | 17.2 µs | 4.85 µs | 3.5x faster |
| `cubic_spline/n100/shuffled-uneven` | 3.03 µs | 3.46 µs | 8.76 µs | 19.6 µs | 4.79 µs | 1.6x faster |
| `pchip/n100/sorted-uneven` | 1.49 µs | 1.58 µs | 8.21 µs | 17.1 µs |  | 5.5x faster |
| `pchip/n100/shuffled-uneven` | 3.32 µs | 3.45 µs | 8.64 µs | 17.5 µs |  | 2.6x faster |
| `linear/n10000/sorted-uneven` | 225 µs | 95.3 µs | 280 µs | 1.66 ms | 974 µs | 1.2x faster |
| `linear/n10000/shuffled-uneven` | 680 µs | 161 µs | 2.66 ms | 3.16 ms | 1.93 ms | 2.8x faster |
| `cubic_spline/n10000/sorted-uneven` | 233 µs | 97.7 µs | 481 µs | 153 ms | 1.05 ms | 2.1x faster |
| `cubic_spline/n10000/shuffled-uneven` | 807 µs | 181 µs | 3.7 ms | 168 ms | 2 ms | 2.5x faster |
| `pchip/n10000/sorted-uneven` | 235 µs | 98.6 µs | 486 µs | 133 ms |  | 2.1x faster |
| `pchip/n10000/shuffled-uneven` | 809 µs | 181 µs | 3.69 ms | 132 ms |  | 4.6x faster |

### Construction

| case | interpx | interpx, parallel | SciPy | scirs2-interpolate | ndarray-interp | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/n10/sorted` | 101 ns | 102 ns | 9.62 µs | 93 ns | 41.8 ns | **2.4x slower** |
| `linear/n10/shuffled` | 101 ns | 93.7 ns | 9.05 µs | 101 ns | 45.5 ns | **2.2x slower** |
| `nearest/n10/sorted` | 82.3 ns | 81.9 ns | 10.5 µs | 92.1 ns |  | 1.1x faster |
| `nearest/n10/shuffled` | 76.6 ns | 83.2 ns | 10.7 µs | 89.6 ns |  | 1.2x faster |
| `previous/n10/sorted` | 82.4 ns | 81.2 ns | 10.2 µs |  |  | 1.2e+02x faster |
| `previous/n10/shuffled` | 82.5 ns | 74.2 ns | 10.2 µs |  |  | 1.2e+02x faster |
| `bspline1/n10/sorted` | 112 ns | 114 ns | 17.2 µs |  |  | 1.5e+02x faster |
| `bspline1/n10/shuffled` | 103 ns | 104 ns | 17.1 µs |  |  | 1.7e+02x faster |
| `bspline2/n10/sorted` | 511 ns | 510 ns | 40.7 µs |  |  | 80x faster |
| `bspline2/n10/shuffled` | 508 ns | 472 ns | 37.5 µs |  |  | 74x faster |
| `bspline3/n10/sorted` | 555 ns | 554 ns | 36.7 µs |  |  | 66x faster |
| `bspline3/n10/shuffled` | 513 ns | 554 ns | 38.7 µs |  |  | 75x faster |
| `cubic_spline/n10/sorted` | 236 ns | 246 ns | 86.7 µs | 383 ns | 525 ns | 1.6x faster |
| `cubic_spline/n10/shuffled` | 259 ns | 246 ns | 85.9 µs | 389 ns | 517 ns | 1.5x faster |
| `pchip/n10/sorted` | 181 ns | 184 ns | 101 µs | 203 ns |  | 1.1x faster |
| `pchip/n10/shuffled` | 194 ns | 183 ns | 111 µs | 203 ns |  | 1x faster |
| `akima/n10/sorted` | 181 ns | 198 ns | 75.9 µs | 209 ns |  | 1.2x faster |
| `akima/n10/shuffled` | 196 ns | 193 ns | 81.6 µs | 201 ns |  | 1x faster |
| `makima/n10/sorted` | 191 ns | 207 ns | 86.2 µs |  |  | 4.5e+02x faster |
| `makima/n10/shuffled` | 208 ns | 192 ns | 82.1 µs |  |  | 4e+02x faster |
| `linear/n50/sorted` | 167 ns | 172 ns | 8.82 µs | 134 ns | 100 ns | **1.7x slower** |
| `linear/n50/shuffled` | 168 ns | 161 ns | 9.22 µs | 134 ns | 109 ns | **1.5x slower** |
| `nearest/n50/sorted` | 107 ns | 111 ns | 11 µs | 132 ns |  | 1.2x faster |
| `nearest/n50/shuffled` | 107 ns | 105 ns | 10.7 µs | 132 ns |  | 1.2x faster |
| `previous/n50/sorted` | 116 ns | 118 ns | 10.5 µs |  |  | 91x faster |
| `previous/n50/shuffled` | 117 ns | 115 ns | 10.5 µs |  |  | 89x faster |
| `bspline1/n50/sorted` | 158 ns | 146 ns | 17.1 µs |  |  | 1.1e+02x faster |
| `bspline1/n50/shuffled` | 156 ns | 156 ns | 17 µs |  |  | 1.1e+02x faster |
| `bspline2/n50/sorted` | 2 µs | 2 µs | 40.1 µs |  |  | 20x faster |
| `bspline2/n50/shuffled` | 1.83 µs | 1.85 µs | 42.9 µs |  |  | 23x faster |
| `bspline3/n50/sorted` | 2.03 µs | 2.03 µs | 40.2 µs |  |  | 20x faster |
| `bspline3/n50/shuffled` | 2.22 µs | 2.05 µs | 42.5 µs |  |  | 19x faster |
| `cubic_spline/n50/sorted` | 829 ns | 912 ns | 89.3 µs | 912 ns | 1.18 µs | 1.1x faster |
| `cubic_spline/n50/shuffled` | 907 ns | 831 ns | 87.3 µs | 990 ns | 1.08 µs | 1.1x faster |
| `pchip/n50/sorted` | 490 ns | 486 ns | 109 µs | 326 ns |  | **1.5x slower** |
| `pchip/n50/shuffled` | 486 ns | 491 ns | 102 µs | 355 ns |  | **1.4x slower** |
| `akima/n50/sorted` | 536 ns | 494 ns | 86.9 µs | 368 ns |  | **1.5x slower** |
| `akima/n50/shuffled` | 548 ns | 493 ns | 87 µs | 364 ns |  | **1.5x slower** |
| `makima/n50/sorted` | 580 ns | 579 ns | 85.8 µs |  |  | 1.5e+02x faster |
| `makima/n50/shuffled` | 581 ns | 531 ns | 79.8 µs |  |  | 1.4e+02x faster |
| `linear/n100/sorted` | 258 ns | 246 ns | 10 µs | 196 ns | 183 ns | **1.4x slower** |
| `linear/n100/shuffled` | 246 ns | 246 ns | 9.11 µs | 216 ns | 169 ns | **1.5x slower** |
| `nearest/n100/sorted` | 176 ns | 176 ns | 10.1 µs | 200 ns |  | 1.1x faster |
| `nearest/n100/shuffled` | 173 ns | 176 ns | 10.1 µs | 216 ns |  | 1.2x faster |
| `previous/n100/sorted` | 176 ns | 176 ns | 9.98 µs |  |  | 57x faster |
| `previous/n100/shuffled` | 175 ns | 161 ns | 10.1 µs |  |  | 58x faster |
| `bspline1/n100/sorted` | 215 ns | 196 ns | 16.2 µs |  |  | 75x faster |
| `bspline1/n100/shuffled` | 203 ns | 215 ns | 17.6 µs |  |  | 87x faster |
| `bspline2/n100/sorted` | 3.82 µs | 3.83 µs | 45.4 µs |  |  | 12x faster |
| `bspline2/n100/shuffled` | 3.83 µs | 3.49 µs | 49.3 µs |  |  | 13x faster |
| `bspline3/n100/sorted` | 4.27 µs | 4.31 µs | 48.7 µs |  |  | 11x faster |
| `bspline3/n100/shuffled` | 3.95 µs | 3.95 µs | 46.3 µs |  |  | 12x faster |
| `cubic_spline/n100/sorted` | 1.71 µs | 1.72 µs | 96.9 µs | 1.81 µs | 1.9 µs | 1.1x faster |
| `cubic_spline/n100/shuffled` | 1.72 µs | 1.72 µs | 96.9 µs | 1.66 µs | 2.05 µs | **1x slower** |
| `pchip/n100/sorted` | 960 ns | 961 ns | 113 µs | 596 ns |  | **1.6x slower** |
| `pchip/n100/shuffled` | 959 ns | 959 ns | 118 µs | 597 ns |  | **1.6x slower** |
| `akima/n100/sorted` | 953 ns | 962 ns | 84.1 µs | 599 ns |  | **1.6x slower** |
| `akima/n100/shuffled` | 955 ns | 910 ns | 89.9 µs | 586 ns |  | **1.6x slower** |
| `makima/n100/sorted` | 989 ns | 1.05 µs | 88.4 µs |  |  | 89x faster |
| `makima/n100/shuffled` | 955 ns | 1.05 µs | 85.9 µs |  |  | 90x faster |
| `linear/n1000/sorted` | 2.12 µs | 2.16 µs | 14 µs | 1.32 µs | 1.53 µs | **1.6x slower** |
| `linear/n1000/shuffled` | 2.11 µs | 2.12 µs | 14.6 µs | 1.27 µs | 1.52 µs | **1.7x slower** |
| `nearest/n1000/sorted` | 1.07 µs | 1.08 µs | 15.9 µs | 1.19 µs |  | 1.1x faster |
| `nearest/n1000/shuffled` | 1.06 µs | 1.07 µs | 15.8 µs | 1.27 µs |  | 1.2x faster |
| `previous/n1000/sorted` | 1.07 µs | 1.09 µs | 16.4 µs |  |  | 15x faster |
| `previous/n1000/shuffled` | 1.07 µs | 1.07 µs | 18.1 µs |  |  | 17x faster |
| `bspline1/n1000/sorted` | 1.1 µs | 1.13 µs | 18.6 µs |  |  | 17x faster |
| `bspline1/n1000/shuffled` | 1.13 µs | 1.13 µs | 18.4 µs |  |  | 16x faster |
| `bspline2/n1000/sorted` | 36.5 µs | 37.4 µs | 111 µs |  |  | 3x faster |
| `bspline2/n1000/shuffled` | 34.9 µs | 36.7 µs | 120 µs |  |  | 3.4x faster |
| `bspline3/n1000/sorted` | 38.2 µs | 42.4 µs | 139 µs |  |  | 3.6x faster |
| `bspline3/n1000/shuffled` | 40.9 µs | 42.5 µs | 144 µs |  |  | 3.5x faster |
| `cubic_spline/n1000/sorted` | 15.2 µs | 16.9 µs | 115 µs | 15.9 µs | 15.8 µs | 1x faster |
| `cubic_spline/n1000/shuffled` | 16.6 µs | 22.9 µs | 123 µs | 15.9 µs | 21 µs | **1x slower** |
| `pchip/n1000/sorted` | 8.98 µs | 9 µs | 238 µs | 4.58 µs |  | **2x slower** |
| `pchip/n1000/shuffled` | 8.96 µs | 9.01 µs | 131 µs | 4.56 µs |  | **2x slower** |
| `akima/n1000/sorted` | 8.78 µs | 8.77 µs | 104 µs | 4.71 µs |  | **1.9x slower** |
| `akima/n1000/shuffled` | 8.73 µs | 8.76 µs | 104 µs | 4.41 µs |  | **2x slower** |
| `makima/n1000/sorted` | 8.91 µs | 9.48 µs | 104 µs |  |  | 12x faster |
| `makima/n1000/shuffled` | 9.68 µs | 15.6 µs | 108 µs |  |  | 11x faster |
| `linear/n10000/sorted` | 22.7 µs | 21.2 µs | 50.2 µs | 12.8 µs | 33.7 µs | **1.8x slower** |
| `linear/n10000/shuffled` | 21.1 µs | 21.6 µs | 47 µs | 12.3 µs | 14.8 µs | **1.7x slower** |
| `nearest/n10000/sorted` | 11.3 µs | 11.4 µs | 59.3 µs | 12.9 µs |  | 1.1x faster |
| `nearest/n10000/shuffled` | 11.4 µs | 11.4 µs | 54.1 µs | 12.9 µs |  | 1.1x faster |
| `previous/n10000/sorted` | 11.4 µs | 11.6 µs | 81.1 µs |  |  | 7.1x faster |
| `previous/n10000/shuffled` | 11.4 µs | 11.4 µs | 81.7 µs |  |  | 7.1x faster |
| `bspline1/n10000/sorted` | 11.4 µs | 11.4 µs | 25.2 µs |  |  | 2.2x faster |
| `bspline1/n10000/shuffled` | 11.4 µs | 11.6 µs | 24.8 µs |  |  | 2.2x faster |
| `bspline2/n10000/sorted` | 521 µs | 525 µs | 813 µs |  |  | 1.6x faster |
| `bspline2/n10000/shuffled` | 514 µs | 598 µs | 818 µs |  |  | 1.6x faster |
| `bspline3/n10000/sorted` | 479 µs | 564 µs | 1.01 ms |  |  | 2.1x faster |
| `bspline3/n10000/shuffled` | 592 µs | 569 µs | 1.01 ms |  |  | 1.7x faster |
| `cubic_spline/n10000/sorted` | 167 µs | 177 µs | 706 µs | 146 µs | 183 µs | **1.1x slower** |
| `cubic_spline/n10000/shuffled` | 167 µs | 192 µs | 795 µs | 147 µs | 183 µs | **1.1x slower** |
| `pchip/n10000/sorted` | 89.9 µs | 197 µs | 257 µs | 47.8 µs |  | **1.9x slower** |
| `pchip/n10000/shuffled` | 94.3 µs | 211 µs | 260 µs | 46.7 µs |  | **2x slower** |
| `akima/n10000/sorted` | 88.5 µs | 182 µs | 235 µs | 43.9 µs |  | **2x slower** |
| `akima/n10000/shuffled` | 106 µs | 92 µs | 687 µs | 47.9 µs |  | **2.2x slower** |
| `makima/n10000/sorted` | 147 µs | 227 µs | 751 µs |  |  | 5.1x faster |
| `makima/n10000/shuffled` | 147 µs | 98.4 µs | 778 µs |  |  | 5.3x faster |
| `linear/n100/sorted-uneven` | 270 ns | 262 ns | 10.1 µs | 203 ns | 167 ns | **1.6x slower** |
| `linear/n100/shuffled-uneven` | 247 ns | 256 ns | 9.58 µs | 197 ns | 167 ns | **1.5x slower** |
| `cubic_spline/n100/sorted-uneven` | 1.56 µs | 1.72 µs | 87.7 µs | 1.8 µs | 2.05 µs | 1.1x faster |
| `cubic_spline/n100/shuffled-uneven` | 1.56 µs | 1.72 µs | 92.2 µs | 1.64 µs | 2.05 µs | 1x faster |
| `pchip/n100/sorted-uneven` | 979 ns | 885 ns | 117 µs | 545 ns |  | **1.8x slower** |
| `pchip/n100/shuffled-uneven` | 958 ns | 959 ns | 112 µs | 590 ns |  | **1.6x slower** |
| `linear/n10000/sorted-uneven` | 22.5 µs | 20.9 µs | 48 µs | 12.8 µs | 20.4 µs | **1.7x slower** |
| `linear/n10000/shuffled-uneven` | 21.1 µs | 21.3 µs | 46.2 µs | 12.9 µs | 21.7 µs | **1.6x slower** |
| `cubic_spline/n10000/sorted-uneven` | 167 µs | 177 µs | 641 µs | 159 µs | 183 µs | **1.1x slower** |
| `cubic_spline/n10000/shuffled-uneven` | 168 µs | 176 µs | 776 µs | 159 µs | 184 µs | **1.1x slower** |
| `pchip/n10000/sorted-uneven` | 102 µs | 159 µs | 280 µs | 48 µs |  | **2.1x slower** |
| `pchip/n10000/shuffled-uneven` | 93.1 µs | 92.8 µs | 259 µs | 46.8 µs |  | **2x slower** |

## `PPoly`

A cubic of 55 intervals, evaluated at points in no order.

### Evaluation

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `intervals55/eval10` | 73.3 ns | 81.1 ns | 3.71 µs | 51x faster |
| `intervals55/eval100` | 727 ns | 806 ns | 4.76 µs | 6.5x faster |
| `intervals55/eval1000` | 7.13 µs | 7.39 µs | 16.9 µs | 2.4x faster |
| `intervals55/eval10000` | 71.8 µs | 72.1 µs | 353 µs | 4.9x faster |

### Construction

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `intervals55/eval10` | 160 ns | 174 ns | 5.71 µs | 36x faster |
| `intervals55/eval100` | 174 ns | 174 ns | 5.68 µs | 33x faster |
| `intervals55/eval1000` | 174 ns | 174 ns | 5.17 µs | 30x faster |
| `intervals55/eval10000` | 173 ns | 173 ns | 5.7 µs | 33x faster |

## `BSpline`

A cubic B-spline of 55 intervals, evaluated at points in no order.

### Evaluation

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `intervals55/eval10` | 271 ns | 257 ns | 2.99 µs | 11x faster |
| `intervals55/eval100` | 1.55 µs | 1.63 µs | 5.76 µs | 3.7x faster |
| `intervals55/eval1000` | 14.5 µs | 15.6 µs | 41.6 µs | 2.9x faster |
| `intervals55/eval10000` | 154 µs | 63.8 µs | 398 µs | 2.6x faster |

### Construction

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `intervals55/eval10` | 294 ns | 247 ns | 9.98 µs | 34x faster |
| `intervals55/eval100` | 177 ns | 177 ns | 9.45 µs | 54x faster |
| `intervals55/eval1000` | 177 ns | 176 ns | 10.1 µs | 57x faster |
| `intervals55/eval10000` | 175 ns | 177 ns | 10.1 µs | 57x faster |

## `BSplineCurve`

A cubic curve through `n` points of the plane, evaluated at `4 n` parameters in order.

### Evaluation

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `degree3/n100` | 4.54 µs | 4.67 µs | 16 µs | 3.5x faster |
| `degree3/n10000` | 492 µs | 165 µs | 1.48 ms | 3x faster |

### Construction

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `degree3/n100` | 4.67 µs | 4.67 µs | 47.7 µs | 10x faster |
| `degree3/n10000` | 712 µs | 683 µs | 1.51 ms | 2.1x faster |

## Regular grid

Points in no order inside the grid. `maxN`: the first axis has `N` points and each axis half the one before (SciPy's suite); `eachN`: every axis has `N` points; `flipped`: every axis given from its high end.

### Evaluation

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | ndarray-interp | interpn | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/2d/max10/eval10/ascending` | 147 ns | 150 ns | 137 ns | 21.3 µs | 800 ns | 167 ns | 104 ns | **1.4x slower** |
| `linear/2d/max10/eval10/flipped` | 146 ns | 149 ns | 135 ns | 19.5 µs |  |  |  | 1.3e+02x faster |
| `linear/2d/max10/eval100/ascending` | 1.25 µs | 1.44 µs | 1.25 µs | 24.3 µs | 8.99 µs | 1.36 µs | 1.02 µs | **1.2x slower** |
| `linear/2d/max10/eval100/flipped` | 1.37 µs | 1.43 µs | 1.32 µs | 22.5 µs |  |  |  | 16x faster |
| `linear/2d/max10/eval1000/ascending` | 13.7 µs | 13.8 µs | 13.6 µs | 46.9 µs | 90.4 µs | 15.3 µs | 10.1 µs | **1.4x slower** |
| `linear/2d/max10/eval1000/flipped` | 13.5 µs | 13.8 µs | 12.4 µs | 44 µs |  |  |  | 3.2x faster |
| `linear/2d/max10/eval10000/ascending` | 137 µs | 56.6 µs | 125 µs | 383 µs | 861 µs | 242 µs | 92.6 µs | **1.5x slower** |
| `linear/2d/max10/eval10000/flipped` | 137 µs | 54.6 µs | 137 µs | 391 µs |  |  |  | 2.9x faster |
| `linear/2d/max40/eval10/ascending` | 172 ns | 175 ns | 172 ns | 20.3 µs | 1.06 µs | 158 ns | 147 ns | **1.2x slower** |
| `linear/2d/max40/eval10/flipped` | 160 ns | 178 ns | 158 ns | 21.2 µs |  |  |  | 1.3e+02x faster |
| `linear/2d/max40/eval100/ascending` | 1.52 µs | 1.77 µs | 1.52 µs | 22.3 µs | 10.8 µs | 1.47 µs | 1.45 µs | **1.1x slower** |
| `linear/2d/max40/eval100/flipped` | 1.66 µs | 1.74 µs | 1.66 µs | 31.3 µs |  |  |  | 19x faster |
| `linear/2d/max40/eval1000/ascending` | 16.6 µs | 16.8 µs | 16.6 µs | 50.4 µs | 99.6 µs | 16.8 µs | 14.4 µs | **1.1x slower** |
| `linear/2d/max40/eval1000/flipped` | 16.6 µs | 16.8 µs | 16.6 µs | 51 µs |  |  |  | 3.1x faster |
| `linear/2d/max40/eval10000/ascending` | 167 µs | 66.5 µs | 166 µs | 618 µs | 1.16 ms | 243 µs | 145 µs | **1.1x slower** |
| `linear/2d/max40/eval10000/flipped` | 166 µs | 61 µs | 157 µs | 627 µs |  |  |  | 3.8x faster |
| `linear/2d/max200/eval10/ascending` | 195 ns | 213 ns | 212 ns | 19.6 µs | 2.01 µs | 188 ns | 176 ns | **1.1x slower** |
| `linear/2d/max200/eval10/flipped` | 199 ns | 207 ns | 194 ns | 19.6 µs |  |  |  | 99x faster |
| `linear/2d/max200/eval100/ascending` | 2.06 µs | 2.1 µs | 1.96 µs | 24.8 µs | 21.3 µs | 2.18 µs | 1.89 µs | **1.1x slower** |
| `linear/2d/max200/eval100/flipped` | 1.91 µs | 1.98 µs | 2.07 µs | 25 µs |  |  |  | 13x faster |
| `linear/2d/max200/eval1000/ascending` | 19.3 µs | 20.9 µs | 20.7 µs | 80.2 µs | 205 µs | 26.7 µs | 17.2 µs | **1.1x slower** |
| `linear/2d/max200/eval1000/flipped` | 20.7 µs | 19.3 µs | 18.9 µs | 69.3 µs |  |  |  | 3.3x faster |
| `linear/2d/max200/eval10000/ascending` | 189 µs | 74.7 µs | 204 µs | 857 µs | 2.29 ms | 488 µs | 190 µs | 1x faster |
| `linear/2d/max200/eval10000/flipped` | 207 µs | 75.5 µs | 207 µs | 884 µs |  |  |  | 4.3x faster |
| `linear/3d/max10/eval10/ascending` | 200 ns | 227 ns | 216 ns | 61.2 µs | 3.17 µs |  | 155 ns | **1.3x slower** |
| `linear/3d/max10/eval10/flipped` | 217 ns | 209 ns | 202 ns | 61.9 µs |  |  |  | 2.8e+02x faster |
| `linear/3d/max10/eval100/ascending` | 2.09 µs | 1.97 µs | 1.94 µs | 65.3 µs | 31.6 µs |  | 1.4 µs | **1.5x slower** |
| `linear/3d/max10/eval100/flipped` | 1.9 µs | 2.16 µs | 1.91 µs | 62.8 µs |  |  |  | 33x faster |
| `linear/3d/max10/eval1000/ascending` | 19.1 µs | 20.8 µs | 20.9 µs | 116 µs | 312 µs |  | 14.1 µs | **1.4x slower** |
| `linear/3d/max10/eval1000/flipped` | 19.3 µs | 19.2 µs | 19 µs | 115 µs |  |  |  | 5.9x faster |
| `linear/3d/max10/eval10000/ascending` | 205 µs | 71.8 µs | 191 µs | 1.16 ms | 3.47 ms |  | 142 µs | **1.4x slower** |
| `linear/3d/max10/eval10000/flipped` | 191 µs | 73.8 µs | 190 µs | 1.15 ms |  |  |  | 6x faster |
| `linear/3d/max40/eval10/ascending` | 222 ns | 266 ns | 241 ns | 56.8 µs | 3.35 µs |  | 201 ns | **1.1x slower** |
| `linear/3d/max40/eval10/flipped` | 243 ns | 246 ns | 240 ns | 63.1 µs |  |  |  | 2.6e+02x faster |
| `linear/3d/max40/eval100/ascending` | 2.33 µs | 2.4 µs | 2.33 µs | 69.4 µs | 34.1 µs |  | 2.17 µs | **1.1x slower** |
| `linear/3d/max40/eval100/flipped` | 2.33 µs | 2.4 µs | 2.34 µs | 64 µs |  |  |  | 27x faster |
| `linear/3d/max40/eval1000/ascending` | 21.3 µs | 22.2 µs | 21.2 µs | 134 µs | 360 µs |  | 21.7 µs | 1x faster |
| `linear/3d/max40/eval1000/flipped` | 21.3 µs | 21.4 µs | 21.4 µs | 143 µs |  |  |  | 6.7x faster |
| `linear/3d/max40/eval10000/ascending` | 212 µs | 80.4 µs | 212 µs | 1.41 ms | 3.47 ms |  | 200 µs | **1.1x slower** |
| `linear/3d/max40/eval10000/flipped` | 212 µs | 77.9 µs | 213 µs | 1.43 ms |  |  |  | 6.7x faster |
| `linear/3d/max200/eval10/ascending` | 302 ns | 308 ns | 301 ns | 56.4 µs | 4.85 µs |  | 286 ns | **1.1x slower** |
| `linear/3d/max200/eval10/flipped` | 301 ns | 309 ns | 303 ns | 62.3 µs |  |  |  | 2.1e+02x faster |
| `linear/3d/max200/eval100/ascending` | 2.97 µs | 3.13 µs | 2.97 µs | 78.4 µs | 51.1 µs |  | 2.87 µs | **1x slower** |
| `linear/3d/max200/eval100/flipped` | 2.97 µs | 3.05 µs | 2.98 µs | 72.3 µs |  |  |  | 24x faster |
| `linear/3d/max200/eval1000/ascending` | 29.7 µs | 30.1 µs | 29.7 µs | 200 µs | 498 µs |  | 28.8 µs | **1x slower** |
| `linear/3d/max200/eval1000/flipped` | 29.6 µs | 29.9 µs | 29.6 µs | 199 µs |  |  |  | 6.7x faster |
| `linear/3d/max200/eval10000/ascending` | 275 µs | 94.7 µs | 297 µs | 1.83 ms | 5 ms |  | 289 µs | 1.1x faster |
| `linear/3d/max200/eval10000/flipped` | 297 µs | 103 µs | 297 µs | 1.92 ms |  |  |  | 6.5x faster |
| `linear/4d/max10/eval10/ascending` | 286 ns | 290 ns | 286 ns | 104 µs |  |  |  | 3.6e+02x faster |
| `linear/4d/max10/eval10/flipped` | 286 ns | 295 ns | 285 ns | 104 µs |  |  |  | 3.6e+02x faster |
| `linear/4d/max10/eval100/ascending` | 2.63 µs | 2.71 µs | 2.63 µs | 117 µs |  |  |  | 44x faster |
| `linear/4d/max10/eval100/flipped` | 2.63 µs | 2.71 µs | 2.65 µs | 117 µs |  |  |  | 44x faster |
| `linear/4d/max10/eval1000/ascending` | 26.7 µs | 26.8 µs | 24.5 µs | 223 µs |  |  |  | 8.4x faster |
| `linear/4d/max10/eval1000/flipped` | 24.4 µs | 27.1 µs | 27.2 µs | 228 µs |  |  |  | 9.3x faster |
| `linear/4d/max10/eval10000/ascending` | 275 µs | 80.7 µs | 258 µs | 2 ms |  |  |  | 7.3x faster |
| `linear/4d/max10/eval10000/flipped` | 267 µs | 76.7 µs | 266 µs | 1.93 ms |  |  |  | 7.2x faster |
| `linear/4d/max40/eval10/ascending` | 310 ns | 318 ns | 339 ns | 99.8 µs | 5.88 µs |  | 273 ns | **1.1x slower** |
| `linear/4d/max40/eval10/flipped` | 325 ns | 323 ns | 310 ns | 96.7 µs |  |  |  | 3e+02x faster |
| `linear/4d/max40/eval100/ascending` | 2.99 µs | 3.19 µs | 3.21 µs | 117 µs | 57.9 µs |  | 2.77 µs | **1.1x slower** |
| `linear/4d/max40/eval100/flipped` | 2.96 µs | 3.27 µs | 3.21 µs | 116 µs |  |  |  | 39x faster |
| `linear/4d/max40/eval1000/ascending` | 32.2 µs | 33 µs | 32 µs | 269 µs | 633 µs |  | 42 µs | 1.3x faster |
| `linear/4d/max40/eval1000/flipped` | 32.1 µs | 32.6 µs | 32.3 µs | 271 µs |  |  |  | 8.4x faster |
| `linear/4d/max40/eval10000/ascending` | 330 µs | 106 µs | 322 µs | 2.47 ms | 6.45 ms |  | 296 µs | **1.1x slower** |
| `linear/4d/max40/eval10000/flipped` | 306 µs | 104 µs | 322 µs | 2.29 ms |  |  |  | 7.5x faster |
| `linear/4d/max200/eval10/ascending` | 386 ns | 422 ns | 416 ns | 95.8 µs | 8.38 µs |  | 359 ns | **1.1x slower** |
| `linear/4d/max200/eval10/flipped` | 392 ns | 425 ns | 416 ns | 96 µs |  |  |  | 2.5e+02x faster |
| `linear/4d/max200/eval100/ascending` | 4.06 µs | 4.18 µs | 3.89 µs | 119 µs | 73.3 µs |  | 3.76 µs | **1.1x slower** |
| `linear/4d/max200/eval100/flipped` | 3.88 µs | 4.16 µs | 4.09 µs | 120 µs |  |  |  | 31x faster |
| `linear/4d/max200/eval1000/ascending` | 41.1 µs | 42 µs | 41.1 µs | 309 µs | 776 µs |  | 39.8 µs | **1x slower** |
| `linear/4d/max200/eval1000/flipped` | 41 µs | 41.6 µs | 41.1 µs | 324 µs |  |  |  | 7.9x faster |
| `linear/4d/max200/eval10000/ascending` | 408 µs | 113 µs | 410 µs | 3 ms | 7.98 ms |  | 365 µs | **1.1x slower** |
| `linear/4d/max200/eval10000/flipped` | 412 µs | 115 µs | 414 µs | 3.27 ms |  |  |  | 7.9x faster |
| `nearest/2d/each40/eval10000/ascending` | 164 µs | 64.7 µs | 151 µs | 648 µs | 1.62 ms |  | 148 µs | **1.1x slower** |
| `nearest/3d/each40/eval10000/ascending` | 213 µs | 76.8 µs | 233 µs | 973 µs | 1.97 ms |  | 227 µs | 1.1x faster |
| `pchip/2d/each40/eval10000/ascending` | 528 µs | 139 µs | 532 µs | 1.41 s |  |  |  | 2.7e+03x faster |
| `pchip/3d/each40/eval10000/ascending` | 1.61 ms | 290 µs | 1.61 ms | 3.59 s |  |  |  | 2.2e+03x faster |
| `cubic/2d/each40/eval10000/ascending` | 330 µs | 126 µs | 330 µs | 1.39 ms |  |  | 449 µs | 1.4x faster |
| `cubic/3d/each40/eval10000/ascending` | 726 µs | 220 µs | 731 µs | 3.2 ms |  |  | 864 µs | 1.2x faster |
| `cubic/2d/each10/eval1000/ascending` | 29.2 µs | 29.4 µs | 29.1 µs | 144 µs |  |  | 40.5 µs | 1.4x faster |
| `cubic/2d/each40/eval1000/ascending` | 32.9 µs | 33.1 µs | 32.9 µs | 156 µs |  |  | 44.6 µs | 1.4x faster |
| `cubic/2d/each100/eval1000/ascending` | 35.2 µs | 35.5 µs | 35.2 µs | 182 µs |  |  | 47.1 µs | 1.3x faster |
| `cubic/2d/each200/eval1000/ascending` | 38 µs | 38.1 µs | 38 µs | 218 µs |  |  | 50.2 µs | 1.3x faster |
| `cubic/2d/each400/eval1000/ascending` | 42.1 µs | 42 µs | 42.1 µs | 299 µs |  |  | 52.7 µs | 1.3x faster |
| `quintic/2d/each10/eval1000/ascending` | 48.5 µs | 48.7 µs | 48.4 µs | 248 µs |  |  |  | 5.1x faster |
| `quintic/2d/each40/eval1000/ascending` | 52.6 µs | 52.8 µs | 52.7 µs | 262 µs |  |  |  | 5x faster |
| `linear/2d/each400/eval10000/ascending` | 241 µs | 85.4 µs | 241 µs | 1.04 ms | 4.51 ms | 586 µs | 241 µs | 1x faster |
| `linear/2d/each800/eval10000/ascending` | 262 µs | 101 µs | 263 µs | 1.18 ms | 8.06 ms | 636 µs | 265 µs | 1x faster |
| `nearest/2d/auckland-ba32-0401/heldout/ascending` | 362 ms | 66.3 ms | 374 ms | 1.41 s | 157 s |  | 755 ms | 2.1x faster |
| `linear/2d/auckland-ba32-0401/heldout/ascending` | 381 ms | 67.6 ms | 386 ms | 1.14 s | 71.2 s | 365 ms | 821 ms | **1x slower** |
| `cubic/2d/auckland-ba32-0401/heldout/ascending` | 779 ms | 135 ms | 786 ms | 33 s |  |  | 1.81 s | 2.3x faster |
| `quintic/2d/auckland-ba32-0401/heldout/ascending` | 1.27 s | 218 ms | 1.28 s | 38.2 s |  |  |  | 30x faster |
| `pchip/2d/auckland-ba32-0401/heldout/ascending` | 1.05 s | 200 ms | 1.09 s |  |  |  |  |  |

### Construction

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | ndarray-interp | interpn | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/2d/max10/eval10/ascending` | 135 ns | 129 ns | 130 ns | 6.37 µs | 134 ns | 85.4 ns | 0.78 ns | **1.7e+02x slower** |
| `linear/2d/max10/eval10/flipped` | 167 ns | 169 ns | 165 ns | 12.6 µs |  |  |  | 75x faster |
| `linear/2d/max10/eval100/ascending` | 122 ns | 138 ns | 130 ns | 6.35 µs | 148 ns | 77.9 ns | 0.78 ns | **1.6e+02x slower** |
| `linear/2d/max10/eval100/flipped` | 168 ns | 172 ns | 166 ns | 11.7 µs |  |  |  | 69x faster |
| `linear/2d/max10/eval1000/ascending` | 136 ns | 134 ns | 149 ns | 5.83 µs | 146 ns | 85.6 ns | 0.88 ns | **1.6e+02x slower** |
| `linear/2d/max10/eval1000/flipped` | 171 ns | 171 ns | 165 ns | 12.7 µs |  |  |  | 74x faster |
| `linear/2d/max10/eval10000/ascending` | 132 ns | 134 ns | 124 ns | 5.84 µs | 134 ns | 85.5 ns | 0.68 ns | **1.9e+02x slower** |
| `linear/2d/max10/eval10000/flipped` | 167 ns | 171 ns | 164 ns | 12.7 µs |  |  |  | 76x faster |
| `linear/2d/max40/eval10/ascending` | 393 ns | 391 ns | 773 ns | 6.37 µs | 219 ns | 166 ns | 0.78 ns | **5e+02x slower** |
| `linear/2d/max40/eval10/flipped` | 702 ns | 714 ns | 1.16 µs | 12.8 µs |  |  |  | 18x faster |
| `linear/2d/max40/eval100/ascending` | 392 ns | 392 ns | 690 ns | 6.1 µs | 221 ns | 153 ns | 0.68 ns | **5.7e+02x slower** |
| `linear/2d/max40/eval100/flipped` | 729 ns | 697 ns | 1.27 µs | 12.6 µs |  |  |  | 17x faster |
| `linear/2d/max40/eval1000/ascending` | 362 ns | 394 ns | 784 ns | 5.79 µs | 219 ns | 167 ns | 0.78 ns | **4.6e+02x slower** |
| `linear/2d/max40/eval1000/flipped` | 702 ns | 699 ns | 1.23 µs | 12.7 µs |  |  |  | 18x faster |
| `linear/2d/max40/eval10000/ascending` | 392 ns | 391 ns | 769 ns | 6.32 µs | 222 ns | 152 ns | 0.68 ns | **5.7e+02x slower** |
| `linear/2d/max40/eval10000/flipped` | 699 ns | 696 ns | 1.19 µs | 12.8 µs |  |  |  | 18x faster |
| `linear/2d/max200/eval10/ascending` | 7.43 µs | 7.28 µs | 20.8 µs | 6.42 µs | 3.1 µs | 515 ns | 0.68 ns | **1.1e+04x slower** |
| `linear/2d/max200/eval10/flipped` | 14.5 µs | 14.6 µs | 28.6 µs | 11.9 µs |  |  |  | **1.2x slower** |
| `linear/2d/max200/eval100/ascending` | 7.45 µs | 6.93 µs | 19 µs | 6.43 µs | 3.25 µs | 515 ns | 0.78 ns | **9.5e+03x slower** |
| `linear/2d/max200/eval100/flipped` | 15.8 µs | 15.9 µs | 25.5 µs | 12.9 µs |  |  |  | **1.2x slower** |
| `linear/2d/max200/eval1000/ascending` | 7.44 µs | 7.39 µs | 15.7 µs | 6.38 µs | 3.24 µs | 516 ns | 0.68 ns | **1.1e+04x slower** |
| `linear/2d/max200/eval1000/flipped` | 14.2 µs | 15.8 µs | 26.7 µs | 11.8 µs |  |  |  | **1.2x slower** |
| `linear/2d/max200/eval10000/ascending` | 7.05 µs | 7.68 µs | 4.74 µs | 5.9 µs | 3.24 µs | 516 ns | 0.78 ns | **9e+03x slower** |
| `linear/2d/max200/eval10000/flipped` | 16.2 µs | 16.1 µs | 20.4 µs | 11.8 µs |  |  |  | **1.4x slower** |
| `linear/3d/max10/eval10/ascending` | 179 ns | 195 ns | 217 ns | 8.25 µs | 179 ns |  | 0.98 ns | **1.8e+02x slower** |
| `linear/3d/max10/eval10/flipped` | 316 ns | 291 ns | 323 ns | 16.5 µs |  |  |  | 52x faster |
| `linear/3d/max10/eval100/ascending` | 194 ns | 180 ns | 226 ns | 8.2 µs | 167 ns |  | 0.88 ns | **2.2e+02x slower** |
| `linear/3d/max10/eval100/flipped` | 288 ns | 313 ns | 328 ns | 15.5 µs |  |  |  | 54x faster |
| `linear/3d/max10/eval1000/ascending` | 178 ns | 195 ns | 244 ns | 8.25 µs | 180 ns |  | 0.88 ns | **2e+02x slower** |
| `linear/3d/max10/eval1000/flipped` | 292 ns | 291 ns | 333 ns | 17.1 µs |  |  |  | 58x faster |
| `linear/3d/max10/eval10000/ascending` | 179 ns | 182 ns | 209 ns | 7.63 µs | 179 ns |  | 0.88 ns | **2e+02x slower** |
| `linear/3d/max10/eval10000/flipped` | 291 ns | 314 ns | 344 ns | 15.4 µs |  |  |  | 53x faster |
| `linear/3d/max40/eval10/ascending` | 2.69 µs | 3.04 µs | 8.39 µs | 8.35 µs | 1.58 µs |  | 0.98 ns | **2.8e+03x slower** |
| `linear/3d/max40/eval10/flipped` | 8.8 µs | 14.6 µs | 12.5 µs | 17.8 µs |  |  |  | 2x faster |
| `linear/3d/max40/eval100/ascending` | 7.52 µs | 3.12 µs | 7.71 µs | 8.3 µs | 1.59 µs |  | 0.98 ns | **7.7e+03x slower** |
| `linear/3d/max40/eval100/flipped` | 15.9 µs | 8.87 µs | 13.5 µs | 15.7 µs |  |  |  | **1x slower** |
| `linear/3d/max40/eval1000/ascending` | 2.91 µs | 2.97 µs | 1.9 µs | 8.2 µs | 1.49 µs |  | 0.98 ns | **3e+03x slower** |
| `linear/3d/max40/eval1000/flipped` | 8.65 µs | 8.83 µs | 9.56 µs | 17.1 µs |  |  |  | 2x faster |
| `linear/3d/max40/eval10000/ascending` | 2.99 µs | 2.94 µs | 1.98 µs | 7.61 µs | 1.51 µs |  | 0.88 ns | **3.4e+03x slower** |
| `linear/3d/max40/eval10000/flipped` | 8.03 µs | 8.89 µs | 9.28 µs | 15.6 µs |  |  |  | 1.9x faster |
| `linear/3d/max200/eval10/ascending` | 1.69 ms | 1.74 ms | 685 µs | 7.62 µs | 1.31 ms |  | 0.98 ns | **1.7e+06x slower** |
| `linear/3d/max200/eval10/flipped` | 3.18 ms | 3.2 ms | 1.9 ms | 17.5 µs |  |  |  | **1.8e+02x slower** |
| `linear/3d/max200/eval100/ascending` | 1.74 ms | 2.84 ms | 659 µs | 8.24 µs | 1.44 ms |  | 1.07 ns | **1.6e+06x slower** |
| `linear/3d/max200/eval100/flipped` | 3.5 ms | 3.25 ms | 1.86 ms | 17.6 µs |  |  |  | **2e+02x slower** |
| `linear/3d/max200/eval1000/ascending` | 2.04 ms | 1.86 ms | 664 µs | 8.27 µs | 1.38 ms |  | 0.98 ns | **2.1e+06x slower** |
| `linear/3d/max200/eval1000/flipped` | 3.32 ms | 3.16 ms | 1.89 ms | 17.2 µs |  |  |  | **1.9e+02x slower** |
| `linear/3d/max200/eval10000/ascending` | 1.81 ms | 1.9 ms | 656 µs | 7.63 µs | 1.51 ms |  | 0.98 ns | **1.8e+06x slower** |
| `linear/3d/max200/eval10000/flipped` | 3.05 ms | 3.48 ms | 1.95 ms | 17.1 µs |  |  |  | **1.8e+02x slower** |
| `linear/4d/max10/eval10/ascending` | 241 ns | 243 ns | 249 ns | 10.6 µs |  |  |  | 44x faster |
| `linear/4d/max10/eval10/flipped` | 358 ns | 365 ns | 384 ns | 19.3 µs |  |  |  | 54x faster |
| `linear/4d/max10/eval100/ascending` | 238 ns | 234 ns | 243 ns | 10.5 µs |  |  |  | 44x faster |
| `linear/4d/max10/eval100/flipped` | 353 ns | 358 ns | 382 ns | 19.4 µs |  |  |  | 55x faster |
| `linear/4d/max10/eval1000/ascending` | 236 ns | 229 ns | 229 ns | 10.5 µs |  |  |  | 45x faster |
| `linear/4d/max10/eval1000/flipped` | 357 ns | 370 ns | 382 ns | 19.6 µs |  |  |  | 55x faster |
| `linear/4d/max10/eval10000/ascending` | 243 ns | 216 ns | 252 ns | 10.5 µs |  |  |  | 43x faster |
| `linear/4d/max10/eval10000/flipped` | 355 ns | 354 ns | 397 ns | 17.6 µs |  |  |  | 50x faster |
| `linear/4d/max40/eval10/ascending` | 14.4 µs | 13.5 µs | 36.6 µs | 9.38 µs | 5.97 µs |  | 1.37 ns | **1.1e+04x slower** |
| `linear/4d/max40/eval10/flipped` | 56.8 µs | 56.5 µs | 69.6 µs | 19.7 µs |  |  |  | **2.9x slower** |
| `linear/4d/max40/eval100/ascending` | 13 µs | 13.2 µs | 38 µs | 10.2 µs | 5.97 µs |  | 1.37 ns | **9.5e+03x slower** |
| `linear/4d/max40/eval100/flipped` | 56.6 µs | 57 µs | 68.9 µs | 19.9 µs |  |  |  | **2.8x slower** |
| `linear/4d/max40/eval1000/ascending` | 14.3 µs | 14.4 µs | 32.7 µs | 10.1 µs | 6.61 µs |  | 2.64 ns | **5.4e+03x slower** |
| `linear/4d/max40/eval1000/flipped` | 56.9 µs | 57.3 µs | 68.2 µs | 21.8 µs |  |  |  | **2.6x slower** |
| `linear/4d/max40/eval10000/ascending` | 14.3 µs | 14.3 µs | 31.6 µs | 10.2 µs | 5.99 µs |  | 1.37 ns | **1e+04x slower** |
| `linear/4d/max40/eval10000/flipped` | 56.4 µs | 56.5 µs | 74.5 µs | 19.7 µs |  |  |  | **2.9x slower** |
| `linear/4d/max200/eval10/ascending` | 51.8 ms | 51.1 ms | 16.8 ms | 9.36 µs | 39.2 ms |  | 1.27 ns | **4.1e+07x slower** |
| `linear/4d/max200/eval10/flipped` | 109 ms | 109 ms | 75.7 ms | 21.7 µs |  |  |  | **5.1e+03x slower** |
| `linear/4d/max200/eval100/ascending` | 50.3 ms | 49.1 ms | 15.7 ms | 10.1 µs | 40.3 ms |  | 1.27 ns | **4e+07x slower** |
| `linear/4d/max200/eval100/flipped` | 109 ms | 109 ms | 71.4 ms | 21.8 µs |  |  |  | **5e+03x slower** |
| `linear/4d/max200/eval1000/ascending` | 48.3 ms | 48.4 ms | 15.8 ms | 10.3 µs | 40 ms |  | 1.37 ns | **3.5e+07x slower** |
| `linear/4d/max200/eval1000/flipped` | 104 ms | 102 ms | 79.4 ms | 21.8 µs |  |  |  | **4.8e+03x slower** |
| `linear/4d/max200/eval10000/ascending` | 48 ms | 48.8 ms | 15.9 ms | 9.65 µs | 40.5 ms |  | 1.27 ns | **3.8e+07x slower** |
| `linear/4d/max200/eval10000/flipped` | 109 ms | 103 ms | 69.5 ms | 21.9 µs |  |  |  | **5e+03x slower** |
| `nearest/2d/each40/eval10000/ascending` | 627 ns | 633 ns | 1.45 µs | 5.86 µs | 289 ns |  | 0.78 ns | **8e+02x slower** |
| `nearest/3d/each40/eval10000/ascending` | 21 µs | 22.4 µs | 53.9 µs | 8.17 µs | 8.91 µs |  | 0.88 ns | **2.4e+04x slower** |
| `pchip/2d/each40/eval10000/ascending` | 7.06 µs | 7.03 µs | 6.84 µs | 7.2 µs |  |  |  | 1x faster |
| `pchip/3d/each40/eval10000/ascending` | 507 µs | 569 µs | 494 µs | 8.53 µs |  |  |  | **59x slower** |
| `cubic/2d/each40/eval10000/ascending` | 9.09 µs | 9.25 µs | 8.93 µs | 1.97 ms |  |  | 0.88 ns | **1e+04x slower** |
| `cubic/3d/each40/eval10000/ascending` | 372 µs | 343 µs | 308 µs | 167 ms |  |  | 1.76 ns | **2.1e+05x slower** |
| `cubic/2d/each10/eval1000/ascending` | 1.59 µs | 1.61 µs | 1.71 µs | 1.12 ms |  |  | 1.76 ns | **9.1e+02x slower** |
| `cubic/2d/each40/eval1000/ascending` | 9.1 µs | 12.6 µs | 8.13 µs | 2.05 ms |  |  | 0.98 ns | **9.3e+03x slower** |
| `cubic/2d/each100/eval1000/ascending` | 44.8 µs | 44.7 µs | 42.7 µs | 7.75 ms |  |  | 0.98 ns | **4.6e+04x slower** |
| `cubic/2d/each200/eval1000/ascending` | 193 µs | 202 µs | 198 µs | 28.3 ms |  |  | 0.88 ns | **2.2e+05x slower** |
| `cubic/2d/each400/eval1000/ascending` | 1.46 ms | 1.46 ms | 1.24 ms | 150 ms |  |  | 0.88 ns | **1.7e+06x slower** |
| `quintic/2d/each10/eval1000/ascending` | 2.1 µs | 2.12 µs | 2.19 µs | 3.73 ms |  |  |  | 1.8e+03x faster |
| `quintic/2d/each40/eval1000/ascending` | 13 µs | 13 µs | 12.8 µs | 8.37 ms |  |  |  | 6.4e+02x faster |
| `linear/2d/each400/eval10000/ascending` | 332 µs | 328 µs | 103 µs | 6.47 µs | 286 µs | 1.25 µs | 0.78 ns | **4.2e+05x slower** |
| `linear/2d/each800/eval10000/ascending` | 1.08 ms | 1.36 ms | 422 µs | 7.03 µs | 859 µs | 2.44 µs | 0.68 ns | **1.6e+06x slower** |
| `nearest/2d/auckland-ba32-0401/heldout/ascending` | 16.6 ms | 16.6 ms | 5.56 ms | 6.94 µs | 12.7 ms |  | 0.68 ns | **2.4e+07x slower** |
| `linear/2d/auckland-ba32-0401/heldout/ascending` | 16.6 ms | 16.4 ms | 5.5 ms | 7.19 µs | 12.9 ms | 8.24 µs | 0.68 ns | **2.4e+07x slower** |
| `cubic/2d/auckland-ba32-0401/heldout/ascending` | 82.8 ms | 78.1 ms | 72.7 ms | 15.1 s |  |  | 0.88 ns | **9.4e+07x slower** |
| `quintic/2d/auckland-ba32-0401/heldout/ascending` | 100 ms | 97.6 ms | 87.6 ms | 21.3 s |  |  |  | 2.1e+02x faster |
| `pchip/2d/auckland-ba32-0401/heldout/ascending` | 79 ms | 77.5 ms | 65.3 ms |  |  |  |  |  |

## `NearestNDInterpolator`

`n` points of the unit square, evaluated at the pixels of a square image.

### Evaluation

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: |
| `n10/image50` | 35.2 µs | 35.4 µs | 112 µs | 444 µs | 13x faster |
| `n10/image100` | 135 µs | 54.9 µs | 486 µs | 1.83 ms | 14x faster |
| `n10/image500` | 3.24 ms | 558 µs | 12.2 ms | 45.3 ms | 14x faster |
| `n100/image50` | 111 µs | 112 µs | 202 µs | 625 µs | 5.6x faster |
| `n100/image100` | 452 µs | 116 µs | 777 µs | 2.28 ms | 5x faster |
| `n100/image500` | 10 ms | 1.58 ms | 19.3 ms | 58.1 ms | 5.8x faster |
| `n1000/image50` | 313 µs | 123 µs | 343 µs | 789 µs | 2.5x faster |
| `n1000/image100` | 1.25 ms | 306 µs | 1.24 ms | 3.17 ms | 2.5x faster |
| `n1000/image500` | 23.2 ms | 4.13 ms | 24.3 ms | 69.6 ms | 3x faster |
| `n10000/image50` | 333 µs | 106 µs | 358 µs | 1.04 ms | 3.1x faster |
| `n10000/image100` | 1.42 ms | 272 µs | 1.49 ms | 3.84 ms | 2.7x faster |
| `n10000/image500` | 33.1 ms | 5.22 ms | 30.7 ms | 91.2 ms | 2.8x faster |
| `n100000/image50` | 485 µs | 119 µs | 486 µs | 1.39 ms | 2.9x faster |
| `n100000/image100` | 1.75 ms | 357 µs | 1.81 ms | 5.33 ms | 3x faster |
| `n100000/image500` | 44.9 ms | 6.89 ms | 48.9 ms | 128 ms | 2.8x faster |
| `ghcn-tmax-20250101/heldout` | 157 µs | 72.2 µs | 173 µs | 563 µs | 3.6x faster |

### Construction

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: |
| `n10/image50` | 565 ns | 562 ns | 528 ns | 9.45 µs | 17x faster |
| `n10/image100` | 570 ns | 587 ns | 564 ns | 10.2 µs | 18x faster |
| `n10/image500` | 569 ns | 569 ns | 554 ns | 9.72 µs | 17x faster |
| `n100/image50` | 3.46 µs | 3.4 µs | 3.39 µs | 15 µs | 4.3x faster |
| `n100/image100` | 3.33 µs | 3.43 µs | 3.34 µs | 14.8 µs | 4.4x faster |
| `n100/image500` | 3.41 µs | 3.33 µs | 3.38 µs | 15 µs | 4.4x faster |
| `n1000/image50` | 34.3 µs | 36.5 µs | 37.5 µs | 117 µs | 3.4x faster |
| `n1000/image100` | 37.6 µs | 35.1 µs | 37.1 µs | 127 µs | 3.4x faster |
| `n1000/image500` | 34.8 µs | 35.4 µs | 36.4 µs | 129 µs | 3.7x faster |
| `n10000/image50` | 525 µs | 546 µs | 490 µs | 1.82 ms | 3.5x faster |
| `n10000/image100` | 503 µs | 535 µs | 502 µs | 1.8 ms | 3.6x faster |
| `n10000/image500` | 531 µs | 535 µs | 540 µs | 1.66 ms | 3.1x faster |
| `n100000/image50` | 8.74 ms | 8.24 ms | 8.27 ms | 25.7 ms | 2.9x faster |
| `n100000/image100` | 7.96 ms | 8.34 ms | 8.43 ms | 24.2 ms | 3x faster |
| `n100000/image500` | 9.87 ms | 9.17 ms | 8.83 ms | 26.6 ms | 2.7x faster |
| `ghcn-tmax-20250101/heldout` | 524 µs | 558 µs | 551 µs | 1.88 ms | 3.6x faster |

## `LinearNDInterpolator`

`n` points of the unit square, evaluated at the pixels of a square image (some outside the hull, with no value); `3d`: `n` points of the unit cube, evaluated at 10 000 points of the cube in no order.

### Evaluation

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | spade | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `n10/image50` | 19.3 µs | 18.4 µs | 17.6 µs | 30.2 µs | 56.9 µs | 96.7 µs | 1.6x faster |
| `n10/image100` | 75.9 µs | 69.9 µs | 69 µs | 104 µs | 235 µs | 383 µs | 1.4x faster |
| `n10/image500` | 1.87 ms | 409 µs | 1.67 ms | 2.83 ms | 5.99 ms | 10.3 ms | 1.5x faster |
| `n100/image50` | 53.5 µs | 54.3 µs | 53.1 µs | 73.9 µs | 940 µs | 210 µs | 1.4x faster |
| `n100/image100` | 225 µs | 81.1 µs | 225 µs | 253 µs | 3.8 ms | 881 µs | 1.1x faster |
| `n100/image500` | 5.19 ms | 991 µs | 5.71 ms | 6.1 ms | 80.3 ms | 18.6 ms | 1.2x faster |
| `n1000/image50` | 54.2 µs | 54.9 µs | 50 µs | 135 µs | 20 ms | 305 µs | 2.5x faster |
| `n1000/image100` | 269 µs | 96.5 µs | 253 µs | 391 µs | 55.7 ms | 1.23 ms | 1.5x faster |
| `n1000/image500` | 5.85 ms | 1.19 ms | 6.1 ms | 7.04 ms | 772 ms | 22 ms | 1.2x faster |
| `n10000/image50` | 66.1 µs | 64.5 µs | 65.6 µs | 391 µs | 311 ms | 617 µs | 5.9x faster |
| `n10000/image100` | 306 µs | 107 µs | 305 µs | 1.01 ms | 1.06 s | 1.64 ms | 3.3x faster |
| `n10000/image500` | 7.71 ms | 1.34 ms | 7.82 ms | 9.31 ms | 20.8 s | 29 ms | 1.2x faster |
| `n100000/image50` | 95 µs | 96.7 µs | 86.6 µs | 1.39 ms | 3.28 s | 2.93 ms | 15x faster |
| `n100000/image100` | 436 µs | 132 µs | 426 µs | 3.78 ms | 12.6 s | 8.79 ms | 8.7x faster |
| `n100000/image500` | 11.8 ms | 2.1 ms | 11.3 ms | 21.3 ms | 302 s | 52.7 ms | 1.8x faster |
| `3d/n1000/scattered` | 792 µs | 196 µs | 779 µs | 5.28 ms |  |  | 6.7x faster |
| `3d/n10000/scattered` | 1.29 ms | 269 µs | 1.26 ms | 16.1 ms |  |  | 12x faster |
| `3d/n100000/scattered` | 2.97 ms | 596 µs | 2.79 ms | 75.4 ms |  |  | 25x faster |
| `ghcn-tmax-20250101/heldout` | 81.4 µs | 84.2 µs | 78.5 µs | 1.18 ms | 147 ms | 1.36 ms | 14x faster |

### Construction

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | spade | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `n10/image50` | 5.7 µs | 5.15 µs | 5.33 µs | 493 µs | 2 µs | 1.72 µs | **3.3x slower** |
| `n10/image100` | 5.74 µs | 5.25 µs | 5.15 µs | 474 µs | 2.03 µs | 1.74 µs | **3.3x slower** |
| `n10/image500` | 5.75 µs | 4.88 µs | 5.21 µs | 494 µs | 2.01 µs | 1.73 µs | **3.3x slower** |
| `n100/image50` | 36.5 µs | 37 µs | 36.2 µs | 728 µs | 58.2 µs | 20.6 µs | **1.8x slower** |
| `n100/image100` | 36.7 µs | 35.6 µs | 35.9 µs | 698 µs | 100 µs | 20.6 µs | **1.8x slower** |
| `n100/image500` | 36.6 µs | 35.5 µs | 36.5 µs | 673 µs | 61.2 µs | 20.7 µs | **1.8x slower** |
| `n1000/image50` | 638 µs | 634 µs | 581 µs | 2.91 ms | 4.71 ms | 304 µs | **2.1x slower** |
| `n1000/image100` | 593 µs | 655 µs | 624 µs | 2.53 ms | 4.59 ms | 304 µs | **2x slower** |
| `n1000/image500` | 585 µs | 661 µs | 643 µs | 2.58 ms | 4.24 ms | 280 µs | **2.1x slower** |
| `n10000/image50` | 8.42 ms | 7.87 ms | 8.53 ms | 26.1 ms | 488 ms | 3.96 ms | **2.1x slower** |
| `n10000/image100` | 11.5 ms | 7.82 ms | 7.9 ms | 25.9 ms | 712 ms | 3.82 ms | **3x slower** |
| `n10000/image500` | 8 ms | 7.86 ms | 7.68 ms | 25 ms | 526 ms | 3.78 ms | **2.1x slower** |
| `n100000/image50` | 87 ms | 79 ms | 79.1 ms | 382 ms | 59.8 s | 48.4 ms | **1.8x slower** |
| `n100000/image100` | 78.9 ms | 81.1 ms | 82.8 ms | 374 ms | 57.7 s | 47 ms | **1.7x slower** |
| `n100000/image500` | 78.9 ms | 80.5 ms | 78.7 ms | 381 ms | 60.9 s | 47.3 ms | **1.7x slower** |
| `3d/n1000/scattered` | 4.67 ms | 4.68 ms | 4.44 ms | 9.54 ms |  |  | 2x faster |
| `3d/n10000/scattered` | 52.5 ms | 51.6 ms | 52.6 ms | 128 ms |  |  | 2.4x faster |
| `3d/n100000/scattered` | 557 ms | 563 ms | 558 ms | 1.76 s |  |  | 3.2x faster |
| `ghcn-tmax-20250101/heldout` | 7.96 ms | 7.99 ms | 7.94 ms | 30.6 ms | 513 ms | 5.04 ms | **1.6x slower** |

## `CloughTocher2DInterpolator`

`n` points of the unit square and a smooth function of them, evaluated at the pixels of a square image (some outside the hull, with no value). SciPy at its defaults (`tol=1e-6`, `maxiter=400`); interpx solves for the gradients until the residual is 1e-14 of the right-hand side.

### Evaluation

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `n10/image50` | 56 µs | 56.6 µs | 52 µs | 105 µs | 71.1 µs | 1.3x faster |
| `n50/image50` | 106 µs | 100 µs | 106 µs | 212 µs | 455 µs | 2x faster |
| `n100/image50` | 110 µs | 112 µs | 112 µs | 230 µs | 1 ms | 2.1x faster |
| `n1000/image500` | 13 ms | 2.32 ms | 13 ms | 24.6 ms | 797 ms | 1.9x faster |
| `n10000/image500` | 14.8 ms | 2.74 ms | 14.8 ms | 27.7 ms | 23.4 s | 1.9x faster |
| `n100000/image500` | 21.4 ms | 5.39 ms | 41 ms | 46.9 ms | 296 s | 2.2x faster |
| `ghcn-tmax-20250101/heldout` | 127 µs | 139 µs | 136 µs | 1.22 ms | 144 ms | 9.6x faster |

### Construction

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `n10/image50` | 9.68 µs | 9.63 µs | 9.24 µs | 511 µs | 3.34 µs | **2.9x slower** |
| `n50/image50` | 39.1 µs | 35.8 µs | 38.6 µs | 595 µs | 28.6 µs | **1.4x slower** |
| `n100/image50` | 74.7 µs | 73.3 µs | 74.6 µs | 923 µs | 105 µs | 1.4x faster |
| `n1000/image500` | 1.06 ms | 1.07 ms | 1.05 ms | 3.38 ms | 6.2 ms | 3.2x faster |
| `n10000/image500` | 13.6 ms | 12.8 ms | 12.8 ms | 37.2 ms | 582 ms | 2.7x faster |
| `n100000/image500` | 156 ms | 197 ms | 296 ms | 542 ms | 72.4 s | 3.5x faster |
| `ghcn-tmax-20250101/heldout` | 12.9 ms | 13.3 ms | 12.8 ms | 45.5 ms | 666 ms | 3.5x faster |

## `RBFInterpolator`

`n` points of a square, evaluated at 1000 points in no order (`scattered`) or at the 10 000 pixels of an image (`image`). `global`: all points; `k50`: the 50 nearest.

### Evaluation

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/global/n10/scattered` | 14.4 µs | 15.9 µs | 15.8 µs | 68.3 µs | 185 µs | 4.8x faster |
| `linear/k50/n10/scattered` | 14.4 µs | 14.5 µs | 14.3 µs | 1.4 ms |  | 97x faster |
| `thin_plate_spline/global/n10/scattered` | 79.7 µs | 80 µs | 73.5 µs | 111 µs | 312 µs | 1.4x faster |
| `thin_plate_spline/k50/n10/scattered` | 73.7 µs | 76 µs | 73 µs | 1.43 ms |  | 19x faster |
| `cubic/global/n10/scattered` | 18.1 µs | 18.4 µs | 18.8 µs | 60.8 µs | 215 µs | 3.4x faster |
| `cubic/k50/n10/scattered` | 29 µs | 18.7 µs | 18.8 µs | 3.39 ms |  | 1.2e+02x faster |
| `quintic/global/n10/scattered` | 40 µs | 27.1 µs | 39.4 µs | 73.2 µs | 253 µs | 1.8x faster |
| `quintic/k50/n10/scattered` | 41.1 µs | 27 µs | 40.9 µs | 1.55 ms |  | 38x faster |
| `multiquadric/global/n10/scattered` | 30.5 µs | 28.2 µs | 30.5 µs | 51.5 µs | 256 µs | 1.7x faster |
| `multiquadric/k50/n10/scattered` | 30.5 µs | 30.7 µs | 27.8 µs | 1.39 ms |  | 46x faster |
| `inverse_multiquadric/global/n10/scattered` | 37.5 µs | 37.6 µs | 37.4 µs | 55.1 µs | 264 µs | 1.5x faster |
| `inverse_multiquadric/k50/n10/scattered` | 37.4 µs | 37.6 µs | 41.1 µs | 1.41 ms |  | 38x faster |
| `inverse_quadratic/global/n10/scattered` | 23.8 µs | 29.8 µs | 23.8 µs | 46.6 µs |  | 2x faster |
| `inverse_quadratic/k50/n10/scattered` | 24 µs | 24.1 µs | 26.4 µs | 1.5 ms |  | 63x faster |
| `gaussian/global/n10/scattered` | 74 µs | 74.2 µs | 74.1 µs | 88 µs | 351 µs | 1.2x faster |
| `gaussian/k50/n10/scattered` | 69.5 µs | 73.7 µs | 67.6 µs | 1.49 ms |  | 21x faster |
| `linear/global/n100/scattered` | 147 µs | 57.7 µs | 135 µs | 291 µs | 1.89 ms | 2x faster |
| `linear/k50/n100/scattered` | 13.9 ms | 2.56 ms | 13.8 ms | 19.2 ms |  | 1.4x faster |
| `thin_plate_spline/global/n100/scattered` | 723 µs | 185 µs | 776 µs | 884 µs | 3.21 ms | 1.2x faster |
| `thin_plate_spline/k50/n100/scattered` | 22.8 ms | 4.52 ms | 23.8 ms | 25.5 ms |  | 1.1x faster |
| `cubic/global/n100/scattered` | 139 µs | 56.7 µs | 152 µs | 291 µs | 1.9 ms | 2.1x faster |
| `cubic/k50/n100/scattered` | 15.5 ms | 3.08 ms | 16 ms | 21.2 ms |  | 1.4x faster |
| `quintic/global/n100/scattered` | 172 µs | 65.9 µs | 168 µs | 306 µs | 2.13 ms | 1.8x faster |
| `quintic/k50/n100/scattered` | 16.2 ms | 3.37 ms | 16.6 ms | 22.4 ms |  | 1.4x faster |
| `multiquadric/global/n100/scattered` | 269 µs | 93.5 µs | 269 µs | 331 µs | 2.38 ms | 1.2x faster |
| `multiquadric/k50/n100/scattered` | 15.2 ms | 2.96 ms | 15.2 ms | 19.6 ms |  | 1.3x faster |
| `inverse_multiquadric/global/n100/scattered` | 393 µs | 108 µs | 364 µs | 398 µs | 2.64 ms | 1x faster |
| `inverse_multiquadric/k50/n100/scattered` | 16.6 ms | 3.28 ms | 16.7 ms | 19.6 ms |  | 1.2x faster |
| `inverse_quadratic/global/n100/scattered` | 224 µs | 76.6 µs | 224 µs | 267 µs |  | 1.2x faster |
| `inverse_quadratic/k50/n100/scattered` | 14.1 ms | 2.69 ms | 13.5 ms | 20.4 ms |  | 1.4x faster |
| `gaussian/global/n100/scattered` | 624 µs | 169 µs | 624 µs | 688 µs | 3.82 ms | 1.1x faster |
| `gaussian/k50/n100/scattered` | 20.6 ms | 3.9 ms | 20.2 ms | 21.8 ms |  | 1.1x faster |
| `linear/global/n1000/scattered` | 1.36 ms | 348 µs | 1.49 ms | 4.28 ms | 20.3 ms | 3.1x faster |
| `linear/k50/n1000/scattered` | 14 ms | 2.68 ms | 13.9 ms | 34.3 ms |  | 2.5x faster |
| `thin_plate_spline/global/n1000/scattered` | 7.77 ms | 1.74 ms | 7.33 ms | 10.8 ms | 30.8 ms | 1.4x faster |
| `thin_plate_spline/k50/n1000/scattered` | 23.3 ms | 4.31 ms | 23.1 ms | 45.5 ms |  | 1.9x faster |
| `cubic/global/n1000/scattered` | 1.49 ms | 343 µs | 1.36 ms | 4.16 ms | 19.7 ms | 2.8x faster |
| `cubic/k50/n1000/scattered` | 16.3 ms | 3.15 ms | 16.5 ms | 38.4 ms |  | 2.4x faster |
| `quintic/global/n1000/scattered` | 1.54 ms | 376 µs | 1.54 ms | 4.01 ms | 21.4 ms | 2.6x faster |
| `quintic/k50/n1000/scattered` | 18.5 ms | 3.48 ms | 17.7 ms | 42.7 ms |  | 2.3x faster |
| `multiquadric/global/n1000/scattered` | 2.98 ms | 629 µs | 3.05 ms | 4.3 ms | 25.9 ms | 1.4x faster |
| `multiquadric/k50/n1000/scattered` | 21 ms | 3.14 ms | 16.3 ms | 40.7 ms |  | 1.9x faster |
| `inverse_multiquadric/global/n1000/scattered` | 3.78 ms | 846 µs | 3.72 ms | 5.23 ms | 26.9 ms | 1.4x faster |
| `inverse_multiquadric/k50/n1000/scattered` | 18.7 ms | 3.38 ms | 17.2 ms | 35.8 ms |  | 1.9x faster |
| `inverse_quadratic/global/n1000/scattered` | 2.47 ms | 524 µs | 2.47 ms | 4.09 ms |  | 1.7x faster |
| `inverse_quadratic/k50/n1000/scattered` | 15.6 ms | 3.11 ms | 16 ms | 34.7 ms |  | 2.2x faster |
| `gaussian/global/n1000/scattered` | 6.8 ms | 1.28 ms | 6.33 ms | 8.61 ms | 38.9 ms | 1.3x faster |
| `gaussian/k50/n1000/scattered` | 22.2 ms | 4.17 ms | 23.2 ms | 43.9 ms |  | 2x faster |
| `thin_plate_spline/global/n3000/scattered` | 22.9 ms | 4.33 ms | 22.5 ms | 30.3 ms | 100 ms | 1.3x faster |
| `thin_plate_spline/k50/n10000/scattered` | 24.8 ms | 6.64 ms | 22.9 ms | 47.6 ms |  | 1.9x faster |
| `thin_plate_spline/k50/n10000/image` | 239 ms | 38.4 ms | 240 ms | 450 ms |  | 1.9x faster |
| `thin_plate_spline/k50/n100000/scattered` | 25.6 ms | 4.61 ms | 24.6 ms | 50.9 ms |  | 2x faster |
| `thin_plate_spline/k50/n100000/image` | 240 ms | 38.6 ms | 243 ms | 476 ms |  | 2x faster |
| `linear/global/ghcn-tmax-20250101/heldout` | 18.2 ms | 3.44 ms | 18.5 ms | 59.1 ms | 239 ms | 3.2x faster |
| `linear/k50/ghcn-tmax-20250101/heldout` | 19.5 ms | 3.62 ms | 19.9 ms | 45.8 ms |  | 2.3x faster |
| `thin_plate_spline/global/ghcn-tmax-20250101/heldout` | 87.2 ms | 16.4 ms | 88.1 ms | 123 ms | 421 ms | 1.4x faster |
| `thin_plate_spline/k50/ghcn-tmax-20250101/heldout` | 31.5 ms | 5.46 ms | 31.8 ms | 60.1 ms |  | 1.9x faster |
| `cubic/global/ghcn-tmax-20250101/heldout` | 18.5 ms | 3.9 ms | 19.8 ms | 53.2 ms | 262 ms | 2.9x faster |
| `cubic/k50/ghcn-tmax-20250101/heldout` | 21.4 ms | 4.29 ms | 21.8 ms | 51.5 ms |  | 2.4x faster |
| `quintic/global/ghcn-tmax-20250101/heldout` | 20.3 ms | 3.81 ms | 20.5 ms | 55.3 ms | 296 ms | 2.7x faster |
| `quintic/k50/ghcn-tmax-20250101/heldout` | 22.5 ms | 4.11 ms | 23.1 ms | 55 ms |  | 2.4x faster |

### Construction

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/global/n10/scattered` | 2.57 µs | 2.83 µs | 2.62 µs | 11.8 µs | 5.37 µs | 2.1x faster |
| `linear/k50/n10/scattered` | 2.61 µs | 2.57 µs | 2.57 µs | 15.6 µs |  | 6x faster |
| `thin_plate_spline/global/n10/scattered` | 3.57 µs | 3.56 µs | 3.37 µs | 12.4 µs | 6.44 µs | 1.8x faster |
| `thin_plate_spline/k50/n10/scattered` | 3.62 µs | 3.26 µs | 3.28 µs | 15.9 µs |  | 4.4x faster |
| `cubic/global/n10/scattered` | 3.74 µs | 3.79 µs | 3.83 µs | 13.2 µs | 6.02 µs | 1.6x faster |
| `cubic/k50/n10/scattered` | 5.84 µs | 3.78 µs | 5.55 µs | 33.5 µs |  | 5.7x faster |
| `quintic/global/n10/scattered` | 6.15 µs | 6.06 µs | 6.36 µs | 15.2 µs | 8.83 µs | 1.4x faster |
| `quintic/k50/n10/scattered` | 6.06 µs | 4.23 µs | 6.55 µs | 32.2 µs |  | 5.3x faster |
| `multiquadric/global/n10/scattered` | 2.91 µs | 2.66 µs | 2.64 µs | 10.6 µs | 7.07 µs | 2.4x faster |
| `multiquadric/k50/n10/scattered` | 2.86 µs | 2.9 µs | 2.62 µs | 15.2 µs |  | 5.3x faster |
| `inverse_multiquadric/global/n10/scattered` | 2.63 µs | 2.62 µs | 2.88 µs | 11.3 µs | 6.96 µs | 2.6x faster |
| `inverse_multiquadric/k50/n10/scattered` | 2.78 µs | 2.63 µs | 2.85 µs | 15.2 µs |  | 5.5x faster |
| `inverse_quadratic/global/n10/scattered` | 2.79 µs | 4.47 µs | 2.55 µs | 10.9 µs |  | 3.9x faster |
| `inverse_quadratic/k50/n10/scattered` | 2.57 µs | 2.55 µs | 2.81 µs | 16.2 µs |  | 6.3x faster |
| `gaussian/global/n10/scattered` | 2.95 µs | 2.79 µs | 3.03 µs | 11 µs | 7.71 µs | 2.6x faster |
| `gaussian/k50/n10/scattered` | 3.01 µs | 3.02 µs | 2.77 µs | 16.3 µs |  | 5.4x faster |
| `linear/global/n100/scattered` | 74.1 µs | 68 µs | 70 µs | 78.5 µs | 770 µs | 1.1x faster |
| `linear/k50/n100/scattered` | 3.19 µs | 3.43 µs | 3.45 µs | 22.6 µs |  | 7.1x faster |
| `thin_plate_spline/global/n100/scattered` | 106 µs | 106 µs | 102 µs | 122 µs | 898 µs | 1.2x faster |
| `thin_plate_spline/k50/n100/scattered` | 3.47 µs | 3.26 µs | 3.32 µs | 21.9 µs |  | 6.3x faster |
| `cubic/global/n100/scattered` | 71.4 µs | 78.1 µs | 76.1 µs | 80.7 µs | 733 µs | 1.1x faster |
| `cubic/k50/n100/scattered` | 3.43 µs | 3.22 µs | 3.56 µs | 23.7 µs |  | 6.9x faster |
| `quintic/global/n100/scattered` | 85 µs | 85 µs | 78 µs | 92.1 µs | 765 µs | 1.1x faster |
| `quintic/k50/n100/scattered` | 3.27 µs | 3.51 µs | 3.28 µs | 23.1 µs |  | 7.1x faster |
| `multiquadric/global/n100/scattered` | 76.1 µs | 78.8 µs | 82 µs | 83.7 µs | 665 µs | 1.1x faster |
| `multiquadric/k50/n100/scattered` | 3.16 µs | 3.41 µs | 3.12 µs | 22.9 µs |  | 7.3x faster |
| `inverse_multiquadric/global/n100/scattered` | 86.3 µs | 86.7 µs | 83.4 µs | 80.2 µs | 786 µs | **1.1x slower** |
| `inverse_multiquadric/k50/n100/scattered` | 3.42 µs | 3.24 µs | 3.39 µs | 21.3 µs |  | 6.2x faster |
| `inverse_quadratic/global/n100/scattered` | 71.9 µs | 79.7 µs | 80.1 µs | 74.5 µs |  | 1x faster |
| `inverse_quadratic/k50/n100/scattered` | 3.31 µs | 3.18 µs | 3.47 µs | 22.8 µs |  | 6.9x faster |
| `gaussian/global/n100/scattered` | 95.3 µs | 105 µs | 100 µs | 97.1 µs | 804 µs | 1x faster |
| `gaussian/k50/n100/scattered` | 3.26 µs | 3.41 µs | 3.35 µs | 20.7 µs |  | 6.3x faster |
| `linear/global/n1000/scattered` | 18.5 ms | 21 ms | 25.6 ms | 58.6 ms | 194 ms | 3.2x faster |
| `linear/k50/n1000/scattered` | 38.4 µs | 39 µs | 38 µs | 134 µs |  | 3.5x faster |
| `thin_plate_spline/global/n1000/scattered` | 25 ms | 23.3 ms | 23.2 ms | 61.7 ms | 206 ms | 2.5x faster |
| `thin_plate_spline/k50/n1000/scattered` | 35.2 µs | 38.5 µs | 35 µs | 136 µs |  | 3.8x faster |
| `cubic/global/n1000/scattered` | 21.8 ms | 20.8 ms | 20.7 ms | 60.3 ms | 206 ms | 2.8x faster |
| `cubic/k50/n1000/scattered` | 35.2 µs | 38.6 µs | 38.2 µs | 136 µs |  | 3.9x faster |
| `quintic/global/n1000/scattered` | 19.9 ms | 20.3 ms | 21.6 ms | 60.1 ms | 192 ms | 3x faster |
| `quintic/k50/n1000/scattered` | 38.6 µs | 38.8 µs | 38.2 µs | 133 µs |  | 3.4x faster |
| `multiquadric/global/n1000/scattered` | 20.9 ms | 21.4 ms | 22.8 ms | 57.9 ms | 214 ms | 2.8x faster |
| `multiquadric/k50/n1000/scattered` | 35.6 µs | 38.5 µs | 36.5 µs | 129 µs |  | 3.6x faster |
| `inverse_multiquadric/global/n1000/scattered` | 21.3 ms | 21.7 ms | 22.1 ms | 55.5 ms | 210 ms | 2.6x faster |
| `inverse_multiquadric/k50/n1000/scattered` | 38.5 µs | 38.4 µs | 35 µs | 130 µs |  | 3.4x faster |
| `inverse_quadratic/global/n1000/scattered` | 20.4 ms | 20.3 ms | 21.5 ms | 60.4 ms |  | 3x faster |
| `inverse_quadratic/k50/n1000/scattered` | 35.2 µs | 38.2 µs | 35.1 µs | 129 µs |  | 3.7x faster |
| `gaussian/global/n1000/scattered` | 24 ms | 24.3 ms | 22 ms | 61.8 ms | 214 ms | 2.6x faster |
| `gaussian/k50/n1000/scattered` | 35.1 µs | 35.4 µs | 37.9 µs | 134 µs |  | 3.8x faster |
| `thin_plate_spline/global/n3000/scattered` | 466 ms | 464 ms | 469 ms | 2.08 s | 9.09 s | 4.5x faster |
| `thin_plate_spline/k50/n10000/scattered` | 747 µs | 714 µs | 608 µs | 1.82 ms |  | 2.4x faster |
| `thin_plate_spline/k50/n10000/image` | 754 µs | 787 µs | 736 µs | 1.81 ms |  | 2.4x faster |
| `thin_plate_spline/k50/n100000/scattered` | 9.57 ms | 9.25 ms | 8.47 ms | 25 ms |  | 2.6x faster |
| `thin_plate_spline/k50/n100000/image` | 8.87 ms | 8.97 ms | 8.63 ms | 24.7 ms |  | 2.8x faster |
| `linear/global/ghcn-tmax-20250101/heldout` | 17.3 s | 17 s | 17.9 s | 113 s | 299 s | 6.5x faster |
| `linear/k50/ghcn-tmax-20250101/heldout` | 780 µs | 762 µs | 668 µs | 1.9 ms |  | 2.4x faster |
| `thin_plate_spline/global/ghcn-tmax-20250101/heldout` | 18.8 s | 18.2 s | 20 s | 119 s | 287 s | 6.3x faster |
| `thin_plate_spline/k50/ghcn-tmax-20250101/heldout` | 759 µs | 810 µs | 744 µs | 1.96 ms |  | 2.6x faster |
| `cubic/global/ghcn-tmax-20250101/heldout` | 17.8 s | 18.4 s | 18.9 s | 117 s | 294 s | 6.6x faster |
| `cubic/k50/ghcn-tmax-20250101/heldout` | 760 µs | 836 µs | 680 µs | 1.89 ms |  | 2.5x faster |
| `quintic/global/ghcn-tmax-20250101/heldout` | 21.9 s | 17.1 s | 18.2 s | 126 s | 289 s | 5.7x faster |
| `quintic/k50/ghcn-tmax-20250101/heldout` | 785 µs | 814 µs | 690 µs | 1.76 ms |  | 2.2x faster |

## Memory

Megabytes. "Data read" is the process with the data of the case in memory, before anything is built. The columns after it are what was added to that: the highest while building, with the interpolator built, the highest while evaluating. "Highest" is the highest of the whole run, repeats included.

| case | library | data read | while building | built | while evaluating | highest |
| :--- | :--- | ---: | ---: | ---: | ---: | ---: |
| `oned/cubic_spline/n10000/sorted` | interpx | 4.2 | +1.0 | +1.1 | +1.1 | 6.0 |
| `oned/cubic_spline/n10000/sorted` | interpx, parallel | 3.5 | +1.7 | +1.8 | +1.8 | 6.4 |
| `oned/cubic_spline/n10000/sorted` | SciPy | 71.3 | +1.4 | +1.8 | +1.8 | 73.3 |
| `oned/cubic_spline/n10000/sorted` | scirs2-interpolate | 4.5 | +0.6 | +0.7 | +1.4 | 5.9 |
| `oned/cubic_spline/n10000/sorted` | ndarray-interp | 3.2 | +1.8 | +1.8 | +1.8 | 5.5 |
| `grid/linear/4d/max200/eval10000/ascending` | interpx | 221.1 | +183.7 | +183.7 | +183.7 | 404.9 |
| `grid/linear/4d/max200/eval10000/ascending` | interpx, parallel | 219.8 | +184.9 | +184.9 | +185.0 | 405.4 |
| `grid/linear/4d/max200/eval10000/ascending` | interpx, `from_owned` | 227.4 | +177.3 | +177.3 | +177.4 | 404.9 |
| `grid/linear/4d/max200/eval10000/ascending` | SciPy | 271.3 | +0.3 | +0.3 | +1.2 | 273.5 |
| `grid/linear/4d/max200/eval10000/ascending` | scirs2-interpolate | 229.8 | +170.0 | +175.2 | +175.2 | 404.9 |
| `grid/linear/4d/max200/eval10000/ascending` | interpn | 221.3 | +0.0 | +0.0 | +0.0 | 400.3 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpx | 6.2 | +0.0 | +2.1 | +2.1 | 8.5 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpx, parallel | 5.8 | +2.0 | +2.5 | +2.5 | 9.0 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpx, `from_owned` | 5.6 | +2.6 | +2.7 | +2.7 | 8.5 |
| `grid/cubic/2d/each400/eval1000/ascending` | SciPy | 72.6 | +84.9 | +5.4 | +5.4 | 162.1 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpn | 3.4 | +2.1 | +2.1 | +2.1 | 5.5 |
| `grid/linear/2d/each800/eval10000/ascending` | interpx | 10.5 | +4.4 | +4.4 | +4.4 | 15.0 |
| `grid/linear/2d/each800/eval10000/ascending` | interpx, parallel | 11.9 | +0.0 | +3.0 | +3.0 | 15.5 |
| `grid/linear/2d/each800/eval10000/ascending` | interpx, `from_owned` | 13.5 | +1.3 | +1.3 | +1.3 | 15.0 |
| `grid/linear/2d/each800/eval10000/ascending` | SciPy | 76.6 | +0.1 | +0.2 | +0.7 | 77.4 |
| `grid/linear/2d/each800/eval10000/ascending` | scirs2-interpolate | 10.6 | +0.0 | +4.2 | +4.3 | 15.0 |
| `grid/linear/2d/each800/eval10000/ascending` | ndarray-interp | 11.2 | +0.0 | +0.0 | +0.0 | 11.2 |
| `grid/linear/2d/each800/eval10000/ascending` | interpn | 11.4 | +0.0 | +0.0 | +0.0 | 11.4 |
| `nearest_nd/n100000/image500` | interpx | 14.6 | +0.7 | +0.0 | +1.8 | 18.5 |
| `nearest_nd/n100000/image500` | interpx, parallel | 11.7 | +3.7 | +2.5 | +5.1 | 18.8 |
| `nearest_nd/n100000/image500` | interpx, `from_owned` | 14.2 | +1.8 | +0.4 | +2.4 | 18.6 |
| `nearest_nd/n100000/image500` | SciPy | 77.4 | +4.5 | +4.5 | +12.7 | 91.5 |
| `linear_nd/n100000/image500` | interpx | 12.9 | +26.5 | +22.2 | +24.0 | 41.5 |
| `linear_nd/n100000/image500` | interpx, parallel | 11.6 | +26.5 | +23.5 | +25.9 | 42.6 |
| `linear_nd/n100000/image500` | interpx, `from_owned` | 13.5 | +26.2 | +21.9 | +23.9 | 42.4 |
| `linear_nd/n100000/image500` | SciPy | 77.5 | +67.4 | +13.8 | +25.6 | 145.7 |
| `linear_nd/n100000/image500` | scirs2-interpolate | 14.3 | +12.3 | +8.0 | +10.1 | 26.6 |
| `linear_nd/n100000/image500` | spade | 14.3 | +13.7 | +10.7 | +12.7 | 29.5 |
| `clough_tocher/n100000/image500` | interpx | 14.5 | +47.1 | +27.0 | +29.0 | 64.1 |
| `clough_tocher/n100000/image500` | interpx, parallel | 13.0 | +48.7 | +28.6 | +31.0 | 66.4 |
| `clough_tocher/n100000/image500` | interpx, `from_owned` | 14.2 | +47.8 | +27.3 | +29.3 | 64.2 |
| `clough_tocher/n100000/image500` | SciPy | 77.5 | +67.6 | +17.7 | +29.9 | 146.3 |
| `clough_tocher/n100000/image500` | scirs2-interpolate | 11.9 | +14.8 | +11.6 | +15.5 | 27.7 |
| `rbf/thin_plate_spline/global/n3000/scattered` | interpx | 3.4 | +77.1 | +4.8 | +4.8 | 80.5 |
| `rbf/thin_plate_spline/global/n3000/scattered` | interpx, parallel | 3.4 | +77.1 | +4.8 | +5.2 | 80.9 |
| `rbf/thin_plate_spline/global/n3000/scattered` | interpx, `from_owned` | 3.3 | +77.2 | +4.8 | +4.9 | 80.6 |
| `rbf/thin_plate_spline/global/n3000/scattered` | SciPy | 71.4 | +73.0 | +0.9 | +9.0 | 144.5 |
| `rbf/thin_plate_spline/global/n3000/scattered` | scirs2-interpolate | 2.9 | +289.7 | +1.6 | +1.7 | 292.5 |
| `rbf/thin_plate_spline/k50/n100000/image` | interpx | 6.7 | +3.7 | +4.7 | +5.0 | 11.9 |
| `rbf/thin_plate_spline/k50/n100000/image` | interpx, parallel | 6.1 | +4.2 | +5.3 | +6.5 | 12.7 |
| `rbf/thin_plate_spline/k50/n100000/image` | interpx, `from_owned` | 7.5 | +5.2 | +3.9 | +4.2 | 13.5 |
| `rbf/thin_plate_spline/k50/n100000/image` | SciPy | 73.6 | +5.3 | +5.4 | +23.6 | 100.7 |

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

`clough_tocher/n100000/image500`

![memory over time, clough_tocher/n100000/image500](performance/clough_tocher-n100000-image500.svg)

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
| `LinearNDInterpolator` | 6e-15 | 6e-15 | 6e-15 | inf |  |  | 6e-15 |
| `CloughTocher2DInterpolator` | 1e-03 | 1e-03 | 1e-03 | inf |  |  |  |
| `RBFInterpolator` | 4e+00 | 4e+00 | 4e+00 | 7e+11 |  |  |  |
