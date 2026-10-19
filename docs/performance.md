# Performance: interpx against SciPy and other Rust crates

<!-- Written by compare/perf/report.py from a run of compare/perf/run.py. Do not edit by hand: run it again. -->

Measured on 2026-10-19, at commit `e6651c8742`.

## How it was measured

- Machine: Intel64 Family 6 Model 191 Stepping 2, GenuineIntel, 16 logical cores, 48 GB. Windows-11-10.0.26300-SP0.
- rustc 1.99.1 (58cf8f9c3 2026-10-13); Python 3.13.1, NumPy 2.5.4, SciPy 1.18.1.
- The other Rust crates, at the versions in `compare/perf/runners/*/Cargo.toml`: `scirs2-interpolate` (the widest of the Rust crates: one variable, grids, RBF), `ndarray-interp` (the most used for one variable; also bilinear), `interpn` (grids).
- Every library ran each case in a process of its own, pinned to one core, on one thread, on the same data. The libraries took turns case by case, 3 rounds.
- `interpx, from_owned` is interpx built from data that is handed over (`from_owned`) and so not copied. Its column is there for the types that keep their data as it is given (the grid, `NearestNDInterpolator`, `RBFInterpolator`), next to the construction from borrowed data (`new`), which is the one the comparison counts. Each of its samples is given data made just before it, which is not in the nearest cache as the data of the other constructions is: on small grids it reads slower than `new` for that reason.
- `interpx, parallel` is interpx evaluated with `eval_batch_into_par` (the `rayon` feature) on 8 threads, not pinned to a core. Every other column is one thread on one core, and the comparison counts those. A parallel batch times its first points and is shared out only when the rest is worth it, so on small cases the column reads as the serial one does.
- A time is that of one call: construction, or the evaluation of all the queries of the case (`eval_batch_into` for interpx, the call with an array for the others). It is the median of the repeats of a process, and then the median of the rounds.
- A call of two seconds or more is timed once: the construction that made the interpolator, or the evaluation that wrote the values, is its measurement, and it is not made again. A library whose construction or evaluation of a case was timed so is run on that case in the first round only, and its time is then one measurement; a call near two seconds that was timed several times goes on to the later rounds. (The rule goes by the time, whatever the library.)
- interpx is called slower where its median is more than 1.10 times another library's and its fastest repeat is slower than that library's fastest.
- A process can be slower for a stretch of a few hundred milliseconds, or from start to end, with nothing changed: seen in interpx's construction and in an evaluation that allocates nothing, so not the code's. Such a stretch raises the median of a round, not the fastest repeat, and alone is not counted as slower; a median that differs between rounds is read with that in mind.
- Memory is the resident set size of the process, sampled every millisecond from outside it.
- Command: `perf/run.py`, run in `compare/`.

## Where interpx is slower

| family | operation | cases | interpx fastest | slower than SciPy | slower than a Rust crate |
| :--- | :--- | ---: | ---: | ---: | ---: |
| One variable | evaluation | 112 | 111 | 0 | 1 |
| One variable | construction | 112 | 75 | 0 | 37 |
| `PPoly` | evaluation | 4 | 4 | 0 | 0 |
| `PPoly` | construction | 4 | 4 | 0 | 0 |
| `BSpline` | evaluation | 4 | 4 | 0 | 0 |
| `BSpline` | construction | 4 | 4 | 0 | 0 |
| `BSplineCurve` | evaluation | 2 | 2 | 0 | 0 |
| `BSplineCurve` | construction | 2 | 2 | 0 | 0 |
| Regular grid | evaluation | 92 | 88 | 0 | 4 |
| Regular grid | construction | 92 | 31 | 34 | 46 |
| `NearestNDInterpolator` | evaluation | 16 | 16 | 0 | 0 |
| `NearestNDInterpolator` | construction | 16 | 16 | 0 | 0 |
| `LinearNDInterpolator` | evaluation | 19 | 19 | 0 | 0 |
| `LinearNDInterpolator` | construction | 19 | 3 | 0 | 16 |
| `CloughTocher2DInterpolator` | evaluation | 7 | 7 | 0 | 0 |
| `CloughTocher2DInterpolator` | construction | 7 | 5 | 0 | 2 |
| `RBFInterpolator` | evaluation | 61 | 61 | 0 | 0 |
| `RBFInterpolator` | construction | 61 | 61 | 0 | 0 |

201 comparisons (a case, an operation, and a library), by method. "Same values" says whether that library's values are SciPy's up to rounding: where they are not, it computes another interpolant under the same name, and its time is not that of the same work.

| method | operation | against | same values | cases | slower by | the largest |
| :--- | :--- | :--- | :--- | ---: | ---: | :--- |
| `grid/cubic` | construction | interpn | no | 8 | 1.67e+03x to 8.57e+07x | `grid/cubic/2d/auckland-ba32-0401/heldout/ascending`: 75.3 ms against 0.88 ns |
| `grid/linear` | construction | interpn | yes | 35 | 178x to 3.92e+07x | `grid/linear/4d/max200/eval10000/ascending`: 49.8 ms against 1.27 ns |
| `grid/nearest` | construction | interpn | yes | 3 | 758x to 2.08e+07x | `grid/nearest/2d/auckland-ba32-0401/heldout/ascending`: 16.2 ms against 0.78 ns |
| `grid/linear` | construction | SciPy | yes | 31 | 1.12x to 5.83e+03x | `grid/linear/4d/max200/eval100/flipped`: 119 ms against 20.5 µs |
| `grid/nearest` | construction | SciPy | yes | 2 | 2.75x to 2.09e+03x | `grid/nearest/2d/auckland-ba32-0401/heldout/ascending`: 16.2 ms against 7.76 µs |
| `grid/linear` | construction | ndarray-interp | yes | 15 | 1.56x to 1.95e+03x | `grid/linear/2d/auckland-ba32-0401/heldout/ascending`: 16.1 ms against 8.25 µs |
| `grid/pchip` | construction | SciPy | yes | 1 | 57.6x | `grid/pchip/3d/each40/eval10000/ascending`: 512 µs against 8.88 µs |
| `oned/pchip` | construction | scirs2-interpolate | yes | 12 | 1.35x to 2.66x | `oned/pchip/n10000/shuffled`: 124 µs against 46.7 µs |
| `oned/akima` | construction | scirs2-interpolate | yes | 8 | 1.51x to 2.53x | `oned/akima/n10000/sorted`: 126 µs against 49.7 µs |
| `grid/linear` | construction | scirs2-interpolate | yes | 29 | 1.11x to 2.34x | `grid/linear/4d/max40/eval1000/ascending`: 13.2 µs against 5.63 µs |
| `oned/linear` | construction | ndarray-interp | yes | 13 | 1.29x to 2.2x | `oned/linear/n10/sorted`: 99.8 ns against 45.4 ns |
| `grid/nearest` | construction | scirs2-interpolate | yes | 3 | 1.24x to 2.12x | `grid/nearest/3d/each40/eval10000/ascending`: 20.8 µs against 9.83 µs |
| `clough_tocher` | construction | scirs2-interpolate | no | 2 | 1.13x to 1.81x | `clough_tocher/n10/image50`: 6.05 µs against 3.34 µs |
| `oned/linear` | construction | scirs2-interpolate | yes | 12 | 1.15x to 1.81x | `oned/linear/n1000/shuffled`: 2.11 µs against 1.17 µs |
| `grid/linear` | evaluation | interpn | yes | 4 | 1.13x to 1.74x | `grid/linear/3d/max10/eval10/ascending`: 271 ns against 156 ns |
| `linear_nd` | construction | spade | yes | 16 | 1.23x to 1.52x | `linear_nd/n10000/image50`: 5.62 ms against 3.69 ms |
| `linear_nd` | construction | scirs2-interpolate | yes | 3 | 1.22x to 1.23x | `linear_nd/n10/image500`: 2.5 µs against 2.03 µs |
| `oned/cubic_spline` | construction | scirs2-interpolate | no | 3 | 1.12x to 1.21x | `oned/cubic_spline/n10000/shuffled-uneven`: 176 µs against 145 µs |
| `oned/nearest` | evaluation | scirs2-interpolate | yes | 1 | 1.1x | `oned/nearest/n10/shuffled`: 335 ns against 303 ns |

<details><summary>Every comparison where interpx is slower</summary>

"Values apart" is the distance of that library's values from SciPy's, as a share of the largest value.

| case | operation | against | interpx | the other | slower by | values apart |
| :--- | :--- | :--- | ---: | ---: | ---: | ---: |
| `grid/cubic/2d/auckland-ba32-0401/heldout/ascending` | construction | interpn | 75.3 ms | 0.88 ns | 85723704.89x | 5e-03 |
| `grid/linear/4d/max200/eval10000/ascending` | construction | interpn | 49.8 ms | 1.27 ns | 39189346.46x | 4e-16 |
| `grid/linear/4d/max200/eval10/ascending` | construction | interpn | 49.7 ms | 1.27 ns | 39187140.92x | 2e-16 |
| `grid/linear/4d/max200/eval1000/ascending` | construction | interpn | 49.6 ms | 1.27 ns | 39057644.31x | 3e-16 |
| `grid/linear/4d/max200/eval100/ascending` | construction | interpn | 49.6 ms | 1.37 ns | 36313965.71x | 3e-16 |
| `grid/linear/2d/auckland-ba32-0401/heldout/ascending` | construction | interpn | 16.1 ms | 0.68 ns | 23580379.43x | 0e+00 |
| `grid/nearest/2d/auckland-ba32-0401/heldout/ascending` | construction | interpn | 16.2 ms | 0.78 ns | 20777728.00x | 0e+00 |
| `grid/linear/3d/max200/eval1000/ascending` | construction | interpn | 1.64 ms | 0.88 ns | 1861632.00x | 3e-16 |
| `grid/linear/3d/max200/eval10000/ascending` | construction | interpn | 1.76 ms | 0.98 ns | 1805260.80x | 5e-16 |
| `grid/linear/3d/max200/eval100/ascending` | construction | interpn | 1.65 ms | 0.98 ns | 1689395.20x | 3e-16 |
| `grid/cubic/2d/each400/eval1000/ascending` | construction | interpn | 1.48 ms | 0.88 ns | 1679587.56x | 5e-01 |
| `grid/linear/3d/max200/eval10/ascending` | construction | interpn | 1.6 ms | 0.98 ns | 1640857.60x | 3e-16 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | interpn | 1.11 ms | 0.68 ns | 1623990.86x | 3e-16 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | interpn | 322 µs | 0.68 ns | 471332.57x | 3e-16 |
| `grid/cubic/3d/each40/eval10000/ascending` | construction | interpn | 317 µs | 1.17 ns | 270336.00x | 5e-01 |
| `grid/cubic/2d/each200/eval1000/ascending` | construction | interpn | 169 µs | 0.88 ns | 192227.56x | 2e-01 |
| `grid/cubic/2d/each100/eval1000/ascending` | construction | interpn | 38.5 µs | 0.88 ns | 43832.89x | 2e-01 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | interpn | 20.8 µs | 0.88 ns | 23715.56x | 0e+00 |
| `grid/linear/2d/max200/eval10000/ascending` | construction | interpn | 7.54 µs | 0.68 ns | 11030.86x | 3e-16 |
| `grid/linear/4d/max40/eval1000/ascending` | construction | interpn | 13.2 µs | 1.27 ns | 10382.77x | 6e-16 |
| `grid/linear/4d/max40/eval10000/ascending` | construction | interpn | 13.1 µs | 1.27 ns | 10294.15x | 5e-16 |
| `grid/linear/4d/max40/eval100/ascending` | construction | interpn | 13 µs | 1.27 ns | 10259.69x | 2e-16 |
| `grid/cubic/2d/each40/eval1000/ascending` | construction | interpn | 8.95 µs | 0.88 ns | 10183.11x | 7e-01 |
| `grid/linear/2d/max200/eval10/ascending` | construction | interpn | 6.79 µs | 0.68 ns | 9938.29x | 2e-16 |
| `grid/linear/4d/max40/eval10/ascending` | construction | interpn | 13.1 µs | 1.37 ns | 9549.71x | 2e-16 |
| `grid/cubic/2d/each40/eval10000/ascending` | construction | interpn | 8.2 µs | 0.88 ns | 9326.22x | 7e-01 |
| `grid/linear/2d/max200/eval100/ascending` | construction | interpn | 6.83 µs | 0.78 ns | 8740.00x | 1e-16 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | interpn | 6.82 µs | 0.78 ns | 8728.00x | 2e-16 |
| `grid/linear/4d/max200/eval100/flipped` | construction | SciPy | 119 ms | 20.5 µs | 5832.26x |  |
| `grid/linear/4d/max200/eval10/ascending` | construction | SciPy | 49.7 ms | 9.37 µs | 5311.90x |  |
| `grid/linear/4d/max200/eval100/ascending` | construction | SciPy | 49.6 ms | 9.4 µs | 5279.95x |  |
| `grid/linear/4d/max200/eval10000/ascending` | construction | SciPy | 49.8 ms | 9.52 µs | 5223.32x |  |
| `grid/linear/4d/max200/eval1000/flipped` | construction | SciPy | 103 ms | 20.1 µs | 5116.48x |  |
| `grid/linear/4d/max200/eval10000/flipped` | construction | SciPy | 103 ms | 20.3 µs | 5090.75x |  |
| `grid/linear/4d/max200/eval1000/ascending` | construction | SciPy | 49.6 ms | 10.2 µs | 4852.35x |  |
| `grid/linear/4d/max200/eval10/flipped` | construction | SciPy | 101 ms | 21.7 µs | 4642.18x |  |
| `grid/linear/3d/max40/eval1000/ascending` | construction | interpn | 2.92 µs | 0.88 ns | 3325.33x | 2e-16 |
| `grid/linear/3d/max40/eval10000/ascending` | construction | interpn | 2.86 µs | 0.88 ns | 3253.33x | 3e-16 |
| `grid/linear/3d/max40/eval100/ascending` | construction | interpn | 2.83 µs | 0.98 ns | 2900.00x | 2e-16 |
| `grid/linear/3d/max40/eval10/ascending` | construction | interpn | 2.75 µs | 0.98 ns | 2815.60x | 2e-16 |
| `grid/linear/2d/auckland-ba32-0401/heldout/ascending` | construction | SciPy | 16.1 ms | 7.16 µs | 2250.53x |  |
| `grid/nearest/2d/auckland-ba32-0401/heldout/ascending` | construction | SciPy | 16.2 ms | 7.76 µs | 2092.84x |  |
| `grid/linear/2d/auckland-ba32-0401/heldout/ascending` | construction | ndarray-interp | 16.1 ms | 8.25 µs | 1954.61x | 0e+00 |
| `grid/cubic/2d/each10/eval1000/ascending` | construction | interpn | 1.47 µs | 0.88 ns | 1673.56x | 3e-01 |
| `grid/nearest/2d/each40/eval10000/ascending` | construction | interpn | 592 ns | 0.78 ns | 758.13x | 0e+00 |
| `grid/linear/2d/max40/eval100/ascending` | construction | interpn | 394 ns | 0.68 ns | 576.43x | 2e-16 |
| `grid/linear/2d/max40/eval1000/ascending` | construction | interpn | 391 ns | 0.68 ns | 571.71x | 3e-16 |
| `grid/linear/2d/max40/eval10/ascending` | construction | interpn | 366 ns | 0.68 ns | 534.86x | 2e-16 |
| `grid/linear/2d/max40/eval10000/ascending` | construction | interpn | 358 ns | 0.68 ns | 523.71x | 4e-16 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | ndarray-interp | 1.11 ms | 2.23 µs | 498.94x | 6e-16 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | ndarray-interp | 322 µs | 1.14 µs | 282.48x | 5e-16 |
| `grid/linear/3d/max10/eval10/ascending` | construction | interpn | 195 ns | 0.88 ns | 222.28x | 2e-16 |
| `grid/linear/3d/max200/eval100/ascending` | construction | SciPy | 1.65 ms | 7.53 µs | 218.97x |  |
| `grid/linear/3d/max10/eval1000/ascending` | construction | interpn | 211 ns | 0.98 ns | 216.50x | 5e-16 |
| `grid/linear/3d/max200/eval10000/ascending` | construction | SciPy | 1.76 ms | 8.51 µs | 207.18x |  |
| `grid/linear/3d/max200/eval1000/flipped` | construction | SciPy | 3.23 ms | 16 µs | 202.74x |  |
| `grid/linear/2d/max10/eval1000/ascending` | construction | interpn | 136 ns | 0.68 ns | 199.57x | 3e-16 |
| `grid/linear/2d/max10/eval10/ascending` | construction | interpn | 134 ns | 0.68 ns | 196.29x | 1e-16 |
| `grid/linear/3d/max10/eval10000/ascending` | construction | interpn | 190 ns | 0.98 ns | 194.90x | 4e-16 |
| `grid/linear/3d/max200/eval1000/ascending` | construction | SciPy | 1.64 ms | 8.45 µs | 193.53x |  |
| `grid/linear/3d/max200/eval10/ascending` | construction | SciPy | 1.6 ms | 8.35 µs | 192.01x |  |
| `grid/linear/3d/max200/eval10/flipped` | construction | SciPy | 2.99 ms | 16.2 µs | 184.21x |  |
| `grid/linear/3d/max10/eval100/ascending` | construction | interpn | 178 ns | 0.98 ns | 182.00x | 3e-16 |
| `grid/linear/3d/max200/eval100/flipped` | construction | SciPy | 3.12 ms | 17.2 µs | 181.00x |  |
| `grid/linear/2d/max10/eval100/ascending` | construction | interpn | 122 ns | 0.68 ns | 178.86x | 2e-16 |
| `grid/linear/2d/max10/eval10000/ascending` | construction | interpn | 121 ns | 0.68 ns | 177.57x | 3e-16 |
| `grid/linear/3d/max200/eval10000/flipped` | construction | SciPy | 3.05 ms | 17.4 µs | 175.31x |  |
| `grid/linear/2d/each800/eval10000/ascending` | construction | SciPy | 1.11 ms | 6.47 µs | 171.70x |  |
| `grid/pchip/3d/each40/eval10000/ascending` | construction | SciPy | 512 µs | 8.88 µs | 57.59x |  |
| `grid/linear/2d/each400/eval10000/ascending` | construction | SciPy | 322 µs | 6.02 µs | 53.53x |  |
| `grid/linear/2d/max200/eval10000/ascending` | construction | ndarray-interp | 7.54 µs | 516 ns | 14.61x | 5e-16 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | ndarray-interp | 6.82 µs | 471 ns | 14.47x | 4e-16 |
| `grid/linear/2d/max200/eval10/ascending` | construction | ndarray-interp | 6.79 µs | 470 ns | 14.46x | 2e-16 |
| `grid/linear/2d/max200/eval100/ascending` | construction | ndarray-interp | 6.83 µs | 515 ns | 13.25x | 3e-16 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | SciPy | 20.8 µs | 7.58 µs | 2.75x |  |
| `oned/pchip/n10000/shuffled` | construction | scirs2-interpolate | 124 µs | 46.7 µs | 2.66x | 2e-16 |
| `oned/pchip/n10000/sorted-uneven` | construction | scirs2-interpolate | 123 µs | 46.9 µs | 2.63x | 3e-16 |
| `grid/linear/4d/max40/eval10000/flipped` | construction | SciPy | 51.2 µs | 19.7 µs | 2.59x |  |
| `grid/linear/4d/max40/eval1000/flipped` | construction | SciPy | 51.3 µs | 19.8 µs | 2.59x |  |
| `grid/linear/4d/max40/eval10/flipped` | construction | SciPy | 51.3 µs | 20 µs | 2.57x |  |
| `grid/linear/2d/max40/eval1000/ascending` | construction | ndarray-interp | 391 ns | 153 ns | 2.56x | 4e-16 |
| `oned/akima/n10000/sorted` | construction | scirs2-interpolate | 126 µs | 49.7 µs | 2.53x | 1e-16 |
| `grid/linear/4d/max40/eval100/flipped` | construction | SciPy | 52.4 µs | 21.7 µs | 2.42x |  |
| `grid/linear/2d/max40/eval10/ascending` | construction | ndarray-interp | 366 ns | 153 ns | 2.39x | 2e-16 |
| `grid/linear/2d/max40/eval100/ascending` | construction | ndarray-interp | 394 ns | 165 ns | 2.39x | 2e-16 |
| `grid/linear/2d/max40/eval10000/ascending` | construction | ndarray-interp | 358 ns | 152 ns | 2.36x | 4e-16 |
| `grid/linear/4d/max40/eval1000/ascending` | construction | scirs2-interpolate | 13.2 µs | 5.63 µs | 2.34x | 2e-16 |
| `grid/linear/2d/max200/eval10000/ascending` | construction | scirs2-interpolate | 7.54 µs | 3.36 µs | 2.25x | 2e-16 |
| `grid/linear/2d/max200/eval10/ascending` | construction | scirs2-interpolate | 6.79 µs | 3.09 µs | 2.20x | 1e-16 |
| `oned/linear/n10/sorted` | construction | ndarray-interp | 99.8 ns | 45.4 ns | 2.20x | 0e+00 |
| `oned/akima/n10000/shuffled` | construction | scirs2-interpolate | 95.1 µs | 43.4 µs | 2.19x | 1e-16 |
| `grid/linear/4d/max40/eval10/ascending` | construction | scirs2-interpolate | 13.1 µs | 5.97 µs | 2.19x | 2e-16 |
| `grid/linear/4d/max40/eval100/ascending` | construction | scirs2-interpolate | 13 µs | 5.96 µs | 2.19x | 3e-16 |
| `grid/linear/4d/max40/eval10000/ascending` | construction | scirs2-interpolate | 13.1 µs | 5.99 µs | 2.18x | 4e-16 |
| `grid/linear/2d/max200/eval100/ascending` | construction | scirs2-interpolate | 6.83 µs | 3.18 µs | 2.15x | 2e-16 |
| `oned/pchip/n10000/sorted` | construction | scirs2-interpolate | 91.8 µs | 42.7 µs | 2.15x | 2e-16 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | scirs2-interpolate | 20.8 µs | 9.83 µs | 2.12x | 0e+00 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | scirs2-interpolate | 6.82 µs | 3.24 µs | 2.11x | 2e-16 |
| `oned/pchip/n10000/shuffled-uneven` | construction | scirs2-interpolate | 93.6 µs | 46.9 µs | 1.99x | 3e-16 |
| `grid/nearest/2d/each40/eval10000/ascending` | construction | scirs2-interpolate | 592 ns | 297 ns | 1.99x | 0e+00 |
| `grid/linear/3d/max40/eval1000/ascending` | construction | scirs2-interpolate | 2.92 µs | 1.47 µs | 1.99x | 2e-16 |
| `grid/linear/2d/max40/eval100/ascending` | construction | scirs2-interpolate | 394 ns | 199 ns | 1.98x | 3e-16 |
| `oned/akima/n1000/sorted` | construction | scirs2-interpolate | 8.5 µs | 4.41 µs | 1.93x | 3e-17 |
| `oned/pchip/n1000/sorted` | construction | scirs2-interpolate | 8.28 µs | 4.29 µs | 1.93x | 2e-16 |
| `grid/linear/3d/max40/eval10/ascending` | construction | scirs2-interpolate | 2.75 µs | 1.43 µs | 1.92x | 1e-16 |
| `grid/linear/3d/max40/eval100/ascending` | construction | scirs2-interpolate | 2.83 µs | 1.47 µs | 1.92x | 1e-16 |
| `oned/akima/n1000/shuffled` | construction | scirs2-interpolate | 8.77 µs | 4.59 µs | 1.91x | 3e-17 |
| `grid/linear/3d/max40/eval10000/ascending` | construction | scirs2-interpolate | 2.86 µs | 1.56 µs | 1.83x | 2e-16 |
| `oned/pchip/n1000/shuffled` | construction | scirs2-interpolate | 8.27 µs | 4.56 µs | 1.81x | 2e-16 |
| `clough_tocher/n10/image50` | construction | scirs2-interpolate | 6.05 µs | 3.34 µs | 1.81x | 1e-01 |
| `oned/linear/n1000/shuffled` | construction | scirs2-interpolate | 2.11 µs | 1.17 µs | 1.81x | 0e+00 |
| `grid/linear/2d/max40/eval10/ascending` | construction | scirs2-interpolate | 366 ns | 203 ns | 1.80x | 4e-17 |
| `oned/linear/n10000/sorted-uneven` | construction | scirs2-interpolate | 22.9 µs | 12.9 µs | 1.78x | 6e-17 |
| `grid/linear/2d/max40/eval1000/ascending` | construction | scirs2-interpolate | 391 ns | 220 ns | 1.77x | 3e-16 |
| `grid/linear/2d/max40/eval10000/ascending` | construction | scirs2-interpolate | 358 ns | 202 ns | 1.77x | 2e-16 |
| `grid/linear/3d/max10/eval10/ascending` | evaluation | interpn | 271 ns | 156 ns | 1.74x | 2e-16 |
| `oned/linear/n10/shuffled` | construction | ndarray-interp | 91.4 ns | 53.1 ns | 1.72x | 0e+00 |
| `oned/linear/n10000/shuffled-uneven` | construction | scirs2-interpolate | 21.4 µs | 12.8 µs | 1.67x | 6e-17 |
| `oned/linear/n50/shuffled` | construction | ndarray-interp | 167 ns | 100 ns | 1.67x | 0e+00 |
| `oned/linear/n10000/sorted` | construction | scirs2-interpolate | 21.4 µs | 12.9 µs | 1.66x | 0e+00 |
| `oned/linear/n10000/shuffled` | construction | scirs2-interpolate | 21.3 µs | 12.8 µs | 1.66x | 0e+00 |
| `oned/linear/n1000/sorted` | construction | scirs2-interpolate | 1.92 µs | 1.16 µs | 1.65x | 0e+00 |
| `oned/akima/n100/sorted` | construction | scirs2-interpolate | 963 ns | 584 ns | 1.65x | 1e-16 |
| `oned/akima/n50/sorted` | construction | scirs2-interpolate | 535 ns | 325 ns | 1.65x | 1e-16 |
| `oned/akima/n100/shuffled` | construction | scirs2-interpolate | 882 ns | 536 ns | 1.65x | 1e-16 |
| `oned/pchip/n100/sorted-uneven` | construction | scirs2-interpolate | 965 ns | 590 ns | 1.64x | 2e-16 |
| `oned/pchip/n100/sorted` | construction | scirs2-interpolate | 964 ns | 592 ns | 1.63x | 2e-16 |
| `oned/pchip/n100/shuffled-uneven` | construction | scirs2-interpolate | 966 ns | 595 ns | 1.63x | 2e-16 |
| `oned/pchip/n100/shuffled` | construction | scirs2-interpolate | 961 ns | 594 ns | 1.62x | 2e-16 |
| `grid/linear/2d/max10/eval1000/ascending` | construction | ndarray-interp | 136 ns | 85.4 ns | 1.60x | 3e-16 |
| `grid/linear/2d/max10/eval10/ascending` | construction | ndarray-interp | 134 ns | 85.6 ns | 1.57x | 1e-16 |
| `grid/linear/2d/max10/eval100/ascending` | construction | ndarray-interp | 122 ns | 78.3 ns | 1.56x | 4e-16 |
| `grid/linear/2d/max10/eval10000/ascending` | construction | ndarray-interp | 121 ns | 77.9 ns | 1.56x | 4e-16 |
| `oned/pchip/n50/shuffled` | construction | scirs2-interpolate | 497 ns | 325 ns | 1.53x | 1e-16 |
| `oned/linear/n10000/sorted-uneven` | construction | ndarray-interp | 22.9 µs | 15 µs | 1.52x | 0e+00 |
| `linear_nd/n10000/image50` | construction | spade | 5.62 ms | 3.69 ms | 1.52x | 8e-16 |
| `oned/linear/n50/sorted` | construction | ndarray-interp | 166 ns | 110 ns | 1.51x | 0e+00 |
| `oned/akima/n50/shuffled` | construction | scirs2-interpolate | 535 ns | 354 ns | 1.51x | 1e-16 |
| `oned/linear/n100/sorted-uneven` | construction | ndarray-interp | 251 ns | 166 ns | 1.51x | 0e+00 |
| `oned/linear/n100/shuffled-uneven` | construction | ndarray-interp | 273 ns | 183 ns | 1.49x | 0e+00 |
| `oned/linear/n100/shuffled` | construction | ndarray-interp | 268 ns | 183 ns | 1.47x | 0e+00 |
| `linear_nd/n10000/image100` | construction | spade | 5.79 ms | 3.94 ms | 1.47x | 8e-16 |
| `linear_nd/n100/image100` | construction | spade | 27.3 µs | 18.7 µs | 1.46x | 9e-16 |
| `linear_nd/n1000/image500` | construction | spade | 406 µs | 278 µs | 1.46x | 3e-15 |
| `linear_nd/n10/image500` | construction | spade | 2.5 µs | 1.72 µs | 1.45x | 9e-16 |
| `linear_nd/n1000/image50` | construction | spade | 401 µs | 277 µs | 1.45x | 1e-15 |
| `linear_nd/n10/image100` | construction | spade | 2.49 µs | 1.72 µs | 1.44x | 8e-16 |
| `linear_nd/n1000/image100` | construction | spade | 400 µs | 277 µs | 1.44x | 1e-15 |
| `oned/linear/n10000/sorted` | construction | ndarray-interp | 21.4 µs | 14.8 µs | 1.44x | 0e+00 |
| `linear_nd/n10000/image500` | construction | spade | 5.66 ms | 4 ms | 1.41x | 3e-15 |
| `grid/linear/4d/max40/eval10000/ascending` | construction | SciPy | 13.1 µs | 9.31 µs | 1.40x |  |
| `grid/linear/4d/max40/eval10/ascending` | construction | SciPy | 13.1 µs | 9.31 µs | 1.40x |  |
| `grid/linear/4d/max40/eval100/ascending` | construction | SciPy | 13 µs | 9.31 µs | 1.40x |  |
| `oned/linear/n1000/shuffled` | construction | ndarray-interp | 2.11 µs | 1.51 µs | 1.39x | 0e+00 |
| `oned/linear/n1000/sorted` | construction | ndarray-interp | 1.92 µs | 1.38 µs | 1.39x | 0e+00 |
| `grid/linear/4d/max40/eval1000/ascending` | construction | SciPy | 13.2 µs | 9.57 µs | 1.38x |  |
| `grid/linear/3d/max200/eval1000/ascending` | construction | scirs2-interpolate | 1.64 ms | 1.19 ms | 1.37x | 3e-16 |
| `oned/pchip/n50/sorted` | construction | scirs2-interpolate | 479 ns | 353 ns | 1.35x | 1e-16 |
| `grid/linear/3d/max200/eval10000/ascending` | construction | scirs2-interpolate | 1.76 ms | 1.3 ms | 1.35x | 3e-16 |
| `oned/linear/n100/sorted` | construction | ndarray-interp | 247 ns | 182 ns | 1.35x | 0e+00 |
| `linear_nd/n100/image500` | construction | spade | 27.4 µs | 20.4 µs | 1.34x | 1e-15 |
| `linear_nd/n100/image50` | construction | spade | 27.3 µs | 20.4 µs | 1.34x | 1e-15 |
| `linear_nd/n10/image50` | construction | spade | 2.49 µs | 1.88 µs | 1.33x | 6e-16 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | scirs2-interpolate | 1.11 ms | 842 µs | 1.32x | 2e-16 |
| `linear_nd/n100000/image500` | construction | spade | 57.6 ms | 44.1 ms | 1.31x | 6e-15 |
| `linear_nd/ghcn-tmax-20250101/heldout` | construction | spade | 5.96 ms | 4.56 ms | 1.31x | 7e-16 |
| `oned/linear/n50/shuffled` | construction | scirs2-interpolate | 167 ns | 128 ns | 1.30x | 0e+00 |
| `grid/linear/2d/max200/eval10000/ascending` | construction | SciPy | 7.54 µs | 5.82 µs | 1.29x |  |
| `oned/linear/n10000/shuffled-uneven` | construction | ndarray-interp | 21.4 µs | 16.6 µs | 1.29x | 0e+00 |
| `grid/linear/4d/max200/eval10/ascending` | construction | scirs2-interpolate | 49.7 ms | 38.6 ms | 1.29x | 2e-16 |
| `grid/linear/4d/max200/eval1000/ascending` | construction | scirs2-interpolate | 49.6 ms | 39 ms | 1.27x | 3e-16 |
| `oned/linear/n100/shuffled-uneven` | construction | scirs2-interpolate | 273 ns | 215 ns | 1.27x | 3e-17 |
| `grid/linear/3d/max200/eval100/ascending` | construction | scirs2-interpolate | 1.65 ms | 1.3 ms | 1.27x | 3e-16 |
| `oned/linear/n50/sorted` | construction | scirs2-interpolate | 166 ns | 132 ns | 1.26x | 0e+00 |
| `grid/linear/4d/max200/eval10000/ascending` | construction | scirs2-interpolate | 49.8 ms | 39.5 ms | 1.26x | 4e-16 |
| `linear_nd/n100000/image100` | construction | spade | 58.5 ms | 46.5 ms | 1.26x | 7e-16 |
| `grid/linear/3d/max200/eval10/ascending` | construction | scirs2-interpolate | 1.6 ms | 1.28 ms | 1.25x | 3e-16 |
| `oned/linear/n100/shuffled` | construction | scirs2-interpolate | 268 ns | 215 ns | 1.25x | 0e+00 |
| `grid/nearest/2d/auckland-ba32-0401/heldout/ascending` | construction | scirs2-interpolate | 16.2 ms | 13 ms | 1.24x | 0e+00 |
| `grid/linear/2d/max10/eval1000/ascending` | evaluation | interpn | 11.5 µs | 9.25 µs | 1.24x | 3e-16 |
| `linear_nd/n100000/image50` | construction | spade | 57.8 ms | 46.9 ms | 1.23x | 8e-16 |
| `linear_nd/n10/image500` | construction | scirs2-interpolate | 2.5 µs | 2.03 µs | 1.23x | 9e-16 |
| `linear_nd/n10/image50` | construction | scirs2-interpolate | 2.49 µs | 2.04 µs | 1.22x | 8e-16 |
| `linear_nd/n10/image100` | construction | scirs2-interpolate | 2.49 µs | 2.04 µs | 1.22x | 9e-16 |
| `grid/linear/2d/max10/eval10/ascending` | evaluation | interpn | 115 ns | 95.1 ns | 1.21x | 1e-16 |
| `oned/cubic_spline/n10000/shuffled-uneven` | construction | scirs2-interpolate | 176 µs | 145 µs | 1.21x | 7e-06 |
| `grid/linear/4d/max200/eval100/ascending` | construction | scirs2-interpolate | 49.6 ms | 41.4 ms | 1.20x | 1e-16 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | scirs2-interpolate | 322 µs | 276 µs | 1.17x | 2e-16 |
| `grid/linear/3d/max10/eval1000/ascending` | construction | scirs2-interpolate | 211 ns | 181 ns | 1.17x | 3e-16 |
| `grid/linear/3d/max10/eval10000/ascending` | construction | scirs2-interpolate | 190 ns | 165 ns | 1.16x | 3e-16 |
| `oned/linear/n100/sorted-uneven` | construction | scirs2-interpolate | 251 ns | 217 ns | 1.15x | 3e-17 |
| `grid/linear/2d/max200/eval10/ascending` | construction | SciPy | 6.79 µs | 5.89 µs | 1.15x |  |
| `oned/cubic_spline/n10000/sorted` | construction | scirs2-interpolate | 167 µs | 145 µs | 1.15x | 4e-05 |
| `oned/linear/n100/sorted` | construction | scirs2-interpolate | 247 ns | 215 ns | 1.15x | 0e+00 |
| `grid/linear/2d/max200/eval10/flipped` | construction | SciPy | 14.7 µs | 12.9 µs | 1.13x |  |
| `clough_tocher/n50/image50` | construction | scirs2-interpolate | 32.3 µs | 28.7 µs | 1.13x | 1e-01 |
| `grid/linear/2d/max10/eval10000/ascending` | evaluation | interpn | 104 µs | 92.6 µs | 1.13x | 3e-16 |
| `oned/cubic_spline/n10000/shuffled` | construction | scirs2-interpolate | 167 µs | 149 µs | 1.12x | 4e-05 |
| `grid/linear/2d/max200/eval10000/flipped` | construction | SciPy | 14.5 µs | 13 µs | 1.12x |  |
| `grid/linear/2d/auckland-ba32-0401/heldout/ascending` | construction | scirs2-interpolate | 16.1 ms | 14.5 ms | 1.11x | 0e+00 |
| `oned/nearest/n10/shuffled` | evaluation | scirs2-interpolate | 335 ns | 303 ns | 1.10x | 0e+00 |

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

A run (one library on one case, a process of its own) is stopped at 15 minutes. The longest run that finished took 430 s, 48% of that: scirs2-interpolate on `linear_nd/n100000/image500`.

No run took more than half the limit (450 s).

## Real-world data

Two data sets, each taken whole in the unit its source gives it out in; what was taken and how it was divided was fixed before anything was measured (`compare/perf/real.py`, docs/design.ja.md section 8.4.7).

- `auckland-ba32-0401`: Auckland Part 1 LiDAR 1m DEM (2024), tile BA32_10000_0401. 4800 by 7200 cells of 1 m. The grid is every second row and column (and the last ones); every other cell is a query, and its elevation is the value held out. Licensed by Regional Software Holdings Limited for re-use under CC-BY-4.0. Hosted by Toitū Te Whenua Land Information New Zealand (https://registry.opendata.aws/nz-elevation).
- `ghcn-tmax-20250101`: GHCN-Daily, maximum temperature of 2025-01-01, every station; longitude and latitude in degrees, as they are. Of every 10 stations in the order of the file the last is a query, and its temperature is the value held out. The RBF kernels are those without a shape parameter. NOAA Global Historical Climatology Network Daily (GHCN-D), from https://registry.opendata.aws/noaa-ghcn (CC0 1.0).

These cases are also in the tables of their families below and in the counts above. Here they are with the error: the largest and the root mean square distance from the values held out (metres, and degrees Celsius), and the number of queries a library gave no value. A library missing from a case does not have it, or did not finish (see the runs that failed).

`grid/nearest/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 16.2 ms | 418 ms | 15.5 | 0.178 | 0 |
| interpx, parallel | 15.6 ms | 66.2 ms | 15.5 | 0.178 | 0 |
| interpx, `from_owned` | 6.04 ms | 366 ms | 15.5 | 0.178 | 0 |
| SciPy | 7.76 µs | 3.57 s | 15.5 | 0.178 | 0 |
| scirs2-interpolate | 13 ms | 203 s | 15.5 | 0.178 | 0 |
| interpn | 0.78 ns | 787 ms | 15.5 | 0.178 | 0 |

`grid/linear/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 16.1 ms | 376 ms | 6.74 | 0.0601 | 0 |
| interpx, parallel | 17 ms | 65.8 ms | 6.74 | 0.0601 | 0 |
| interpx, `from_owned` | 5.66 ms | 393 ms | 6.74 | 0.0601 | 0 |
| SciPy | 7.16 µs | 1.1 s | 6.74 | 0.0601 | 0 |
| scirs2-interpolate | 14.5 ms | 71.9 s | 6.74 | 0.0601 | 0 |
| ndarray-interp | 8.25 µs | 357 ms | 6.74 | 0.0601 | 0 |
| interpn | 0.68 ns | 820 ms | 6.74 | 0.0601 | 0 |

`grid/cubic/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 75.3 ms | 775 ms | 6.77 | 0.0507 | 0 |
| interpx, parallel | 73.2 ms | 126 ms | 6.77 | 0.0507 | 0 |
| interpx, `from_owned` | 70.4 ms | 748 ms | 6.77 | 0.0507 | 0 |
| SciPy | 14.6 s | 35.4 s | 6.77 | 0.0507 | 0 |
| interpn | 0.88 ns | 1.7 s | 6.76 | 0.0517 | 0 |

`grid/quintic/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 99 ms | 1.27 s | 6.77 | 0.0517 | 0 |
| interpx, parallel | 94.3 ms | 244 ms | 6.77 | 0.0517 | 0 |
| interpx, `from_owned` | 83.4 ms | 1.28 s | 6.77 | 0.0517 | 0 |
| SciPy | 21.8 s | 49.3 s | 6.77 | 0.0517 | 0 |

`grid/pchip/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 76.6 ms | 1.03 s | 6.75 | 0.049 | 0 |
| interpx, parallel | 74.6 ms | 181 ms | 6.75 | 0.049 | 0 |
| interpx, `from_owned` | 65.5 ms | 1.08 s | 6.75 | 0.049 | 0 |

`nearest_nd/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 539 µs | 160 µs | 25 | 3.42 | 0 |
| interpx, parallel | 516 µs | 74.1 µs | 25 | 3.42 | 0 |
| interpx, `from_owned` | 560 µs | 159 µs | 25 | 3.42 | 0 |
| SciPy | 1.72 ms | 520 µs | 25 | 3.42 | 0 |

`linear_nd/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 5.96 ms | 67.5 µs | 22.1 | 2.83 | 1 |
| interpx, parallel | 5.45 ms | 71.8 µs | 22.1 | 2.83 | 1 |
| interpx, `from_owned` | 5.64 ms | 60.4 µs | 22.1 | 2.83 | 1 |
| SciPy | 30.9 ms | 1.15 ms | 22.1 | 2.83 | 1 |
| scirs2-interpolate | 498 ms | 133 ms | 22.1 | 2.83 | 1 |
| spade | 4.56 ms | 1.21 ms | 22.1 | 2.83 | 1 |

`clough_tocher/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 11.3 ms | 128 µs | 699 | 28 | 1 |
| interpx, parallel | 10.8 ms | 123 µs | 699 | 28 | 1 |
| interpx, `from_owned` | 11.5 ms | 119 µs | 699 | 28 | 1 |
| SciPy | 43.8 ms | 1.18 ms | 699 | 28 | 1 |
| scirs2-interpolate | 604 ms | 139 ms | 22.3 | 2.89 | 1 |

`rbf/linear/global/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 20.5 s | 25.2 ms | 20.8 | 2.74 | 0 |
| interpx, parallel | 16.6 s | 3.64 ms | 20.8 | 2.74 | 0 |
| interpx, `from_owned` | 16.7 s | 18 ms | 20.8 | 2.74 | 0 |
| SciPy | 119 s | 54.5 ms | 20.8 | 2.74 | 0 |
| scirs2-interpolate | 295 s | 241 ms | 20.8 | 2.74 | 0 |

`rbf/linear/k50/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 698 µs | 18.2 ms | 20.8 | 2.75 | 0 |
| interpx, parallel | 811 µs | 3.23 ms | 20.8 | 2.75 | 0 |
| interpx, `from_owned` | 695 µs | 18.9 ms | 20.8 | 2.75 | 0 |
| SciPy | 1.71 ms | 45.2 ms | 20.8 | 2.75 | 0 |

`rbf/thin_plate_spline/global/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 17 s | 79.6 ms | 102 | 5.67 | 0 |
| interpx, parallel | 17 s | 16 ms | 102 | 5.67 | 0 |
| interpx, `from_owned` | 17 s | 79.7 ms | 102 | 5.67 | 0 |
| SciPy | 109 s | 109 ms | 102 | 5.67 | 0 |
| scirs2-interpolate | 285 s | 390 ms | 93.7 | 5.41 | 0 |

`rbf/thin_plate_spline/k50/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 733 µs | 30.8 ms | 102 | 5.66 | 0 |
| interpx, parallel | 735 µs | 5.21 ms | 102 | 5.66 | 0 |
| interpx, `from_owned` | 630 µs | 31.5 ms | 102 | 5.66 | 0 |
| SciPy | 1.92 ms | 55.1 ms | 102 | 5.66 | 0 |

`rbf/cubic/global/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 18.8 s | 18.6 ms | 587 | 27.1 | 0 |
| interpx, parallel | 17.4 s | 3.58 ms | 587 | 27.1 | 0 |
| interpx, `from_owned` | 17.9 s | 19.2 ms | 587 | 27.1 | 0 |
| SciPy | 149 s | 55.9 ms | 597 | 24.4 | 0 |
| scirs2-interpolate | 296 s | 264 ms | 100 | 8.48 | 0 |

`rbf/cubic/k50/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 718 µs | 20.6 ms | 567 | 25.3 | 0 |
| interpx, parallel | 728 µs | 3.68 ms | 567 | 25.3 | 0 |
| interpx, `from_owned` | 599 µs | 20.4 ms | 567 | 25.3 | 0 |
| SciPy | 1.76 ms | 46.4 ms | 567 | 25.3 | 0 |

`rbf/quintic/global/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 16.7 s | 19.9 ms | 2.55e+05 | 1.36e+04 | 0 |
| interpx, parallel | 20.8 s | 4.4 ms | 2.55e+05 | 1.36e+04 | 0 |
| interpx, `from_owned` | 21.3 s | 36.7 ms | 2.55e+05 | 1.36e+04 | 0 |
| SciPy | 111 s | 49.7 ms | 9.39e+04 | 4.92e+03 | 0 |
| scirs2-interpolate | 353 s | 298 ms | 6.56e+16 | 2.85e+15 | 0 |

`rbf/quintic/k50/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 744 µs | 23.4 ms | 3.74e+03 | 199 | 0 |
| interpx, parallel | 790 µs | 3.98 ms | 3.74e+03 | 199 | 0 |
| interpx, `from_owned` | 741 µs | 23.3 ms | 3.74e+03 | 199 | 0 |
| SciPy | 1.81 ms | 53.6 ms | 3.74e+03 | 200 | 0 |

## One variable

Data of `n` points, evaluated at `4 n` points: in order (`sorted`) and in no order (`shuffled`).

### Evaluation

| case | interpx | interpx, parallel | SciPy | scirs2-interpolate | ndarray-interp | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/n10/sorted` | 132 ns | 192 ns | 5.45 µs | 285 ns | 279 ns | 2.1x faster |
| `linear/n10/shuffled` | 173 ns | 248 ns | 5.58 µs | 280 ns | 274 ns | 1.6x faster |
| `nearest/n10/sorted` | 158 ns | 212 ns | 10.4 µs | 294 ns |  | 1.9x faster |
| `nearest/n10/shuffled` | 335 ns | 404 ns | 10.5 µs | 303 ns |  | **1.1x slower** |
| `previous/n10/sorted` | 167 ns | 214 ns | 11.9 µs |  |  | 71x faster |
| `previous/n10/shuffled` | 290 ns | 353 ns | 11.8 µs |  |  | 41x faster |
| `bspline1/n10/sorted` | 232 ns | 318 ns | 2.68 µs |  |  | 12x faster |
| `bspline1/n10/shuffled` | 548 ns | 659 ns | 2.92 µs |  |  | 5.3x faster |
| `bspline2/n10/sorted` | 299 ns | 374 ns | 3.11 µs |  |  | 10x faster |
| `bspline2/n10/shuffled` | 339 ns | 456 ns | 3.14 µs |  |  | 9.3x faster |
| `bspline3/n10/sorted` | 354 ns | 433 ns | 3.9 µs |  |  | 11x faster |
| `bspline3/n10/shuffled` | 440 ns | 537 ns | 3.66 µs |  |  | 8.3x faster |
| `cubic_spline/n10/sorted` | 141 ns | 211 ns | 3.66 µs | 273 ns | 332 ns | 1.9x faster |
| `cubic_spline/n10/shuffled` | 205 ns | 267 ns | 3.99 µs | 279 ns | 357 ns | 1.4x faster |
| `pchip/n10/sorted` | 140 ns | 193 ns | 3.63 µs | 293 ns |  | 2.1x faster |
| `pchip/n10/shuffled` | 219 ns | 266 ns | 3.7 µs | 301 ns |  | 1.4x faster |
| `akima/n10/sorted` | 138 ns | 193 ns | 3.64 µs | 252 ns |  | 1.8x faster |
| `akima/n10/shuffled` | 207 ns | 259 ns | 3.7 µs | 262 ns |  | 1.3x faster |
| `makima/n10/sorted` | 153 ns | 211 ns | 3.63 µs |  |  | 24x faster |
| `makima/n10/shuffled` | 225 ns | 261 ns | 3.91 µs |  |  | 17x faster |
| `linear/n50/sorted` | 657 ns | 721 ns | 6.17 µs | 1.74 µs | 1.23 µs | 1.9x faster |
| `linear/n50/shuffled` | 894 ns | 869 ns | 8.3 µs | 1.65 µs | 1.26 µs | 1.4x faster |
| `nearest/n50/sorted` | 826 ns | 799 ns | 11.2 µs | 1.85 µs |  | 2.2x faster |
| `nearest/n50/shuffled` | 1.6 µs | 1.81 µs | 11.6 µs | 1.7 µs |  | 1.1x faster |
| `previous/n50/sorted` | 838 ns | 818 ns | 14 µs |  |  | 17x faster |
| `previous/n50/shuffled` | 1.51 µs | 1.57 µs | 14 µs |  |  | 9.3x faster |
| `bspline1/n50/sorted` | 1.22 µs | 1.38 µs | 4.24 µs |  |  | 3.5x faster |
| `bspline1/n50/shuffled` | 3.1 µs | 3.16 µs | 5.56 µs |  |  | 1.8x faster |
| `bspline2/n50/sorted` | 1.36 µs | 1.61 µs | 6.62 µs |  |  | 4.9x faster |
| `bspline2/n50/shuffled` | 2 µs | 2.28 µs | 6.68 µs |  |  | 3.3x faster |
| `bspline3/n50/sorted` | 1.85 µs | 1.95 µs | 8.62 µs |  |  | 4.7x faster |
| `bspline3/n50/shuffled` | 2.79 µs | 2.87 µs | 9.35 µs |  |  | 3.4x faster |
| `cubic_spline/n50/sorted` | 762 ns | 816 ns | 5.62 µs | 4.27 µs | 1.6 µs | 2.1x faster |
| `cubic_spline/n50/shuffled` | 988 ns | 1.14 µs | 5.53 µs | 4.11 µs | 1.6 µs | 1.6x faster |
| `pchip/n50/sorted` | 703 ns | 752 ns | 5.16 µs | 4.26 µs |  | 6.1x faster |
| `pchip/n50/shuffled` | 1.07 µs | 1.12 µs | 5.42 µs | 3.69 µs |  | 3.4x faster |
| `akima/n50/sorted` | 764 ns | 807 ns | 5.53 µs | 4.31 µs |  | 5.6x faster |
| `akima/n50/shuffled` | 1.08 µs | 1.15 µs | 6.24 µs | 4.49 µs |  | 4.1x faster |
| `makima/n50/sorted` | 762 ns | 821 ns | 5.67 µs |  |  | 7.4x faster |
| `makima/n50/shuffled` | 1.07 µs | 1.05 µs | 6.1 µs |  |  | 5.7x faster |
| `linear/n100/sorted` | 1.35 µs | 1.53 µs | 7.89 µs | 3.73 µs | 2.7 µs | 2x faster |
| `linear/n100/shuffled` | 1.78 µs | 1.89 µs | 14 µs | 3.95 µs | 2.73 µs | 1.5x faster |
| `nearest/n100/sorted` | 1.54 µs | 1.65 µs | 13.7 µs | 3.91 µs |  | 2.5x faster |
| `nearest/n100/shuffled` | 3.48 µs | 3.4 µs | 14 µs | 4.08 µs |  | 1.2x faster |
| `previous/n100/sorted` | 1.66 µs | 1.76 µs | 16.3 µs |  |  | 9.8x faster |
| `previous/n100/shuffled` | 2.75 µs | 3.15 µs | 16.8 µs |  |  | 6.1x faster |
| `bspline1/n100/sorted` | 2.42 µs | 2.78 µs | 5.81 µs |  |  | 2.4x faster |
| `bspline1/n100/shuffled` | 6.22 µs | 5.92 µs | 12.8 µs |  |  | 2.1x faster |
| `bspline2/n100/sorted` | 2.98 µs | 3.67 µs | 10.9 µs |  |  | 3.6x faster |
| `bspline2/n100/shuffled` | 4.49 µs | 4.62 µs | 13.8 µs |  |  | 3.1x faster |
| `bspline3/n100/sorted` | 3.67 µs | 3.8 µs | 14.9 µs |  |  | 4.1x faster |
| `bspline3/n100/shuffled` | 10.7 µs | 6.8 µs | 18 µs |  |  | 1.7x faster |
| `cubic_spline/n100/sorted` | 1.52 µs | 1.52 µs | 6.92 µs | 17.3 µs | 3.21 µs | 2.1x faster |
| `cubic_spline/n100/shuffled` | 1.96 µs | 2.25 µs | 8.52 µs | 19.4 µs | 3.2 µs | 1.6x faster |
| `pchip/n100/sorted` | 1.58 µs | 1.65 µs | 7.41 µs | 16.4 µs |  | 4.7x faster |
| `pchip/n100/shuffled` | 2.15 µs | 2.07 µs | 8.04 µs | 17.4 µs |  | 3.7x faster |
| `akima/n100/sorted` | 1.59 µs | 1.55 µs | 7.41 µs | 18.4 µs |  | 4.7x faster |
| `akima/n100/shuffled` | 2.15 µs | 2.26 µs | 8.57 µs | 17.6 µs |  | 4x faster |
| `makima/n100/sorted` | 1.56 µs | 1.65 µs | 6.92 µs |  |  | 4.4x faster |
| `makima/n100/shuffled` | 2.15 µs | 2.1 µs | 7.96 µs |  |  | 3.7x faster |
| `linear/n1000/sorted` | 13.4 µs | 15.2 µs | 29.4 µs | 119 µs | 24.6 µs | 1.8x faster |
| `linear/n1000/shuffled` | 17.8 µs | 17.9 µs | 156 µs | 221 µs | 27.1 µs | 1.5x faster |
| `nearest/n1000/sorted` | 16.4 µs | 15.9 µs | 59.6 µs | 172 µs |  | 3.6x faster |
| `nearest/n1000/shuffled` | 32.7 µs | 32.8 µs | 217 µs | 239 µs |  | 6.6x faster |
| `previous/n1000/sorted` | 16.6 µs | 16 µs | 60.7 µs |  |  | 3.7x faster |
| `previous/n1000/shuffled` | 30.3 µs | 27.8 µs | 224 µs |  |  | 7.4x faster |
| `bspline1/n1000/sorted` | 25.4 µs | 26.1 µs | 30.6 µs |  |  | 1.2x faster |
| `bspline1/n1000/shuffled` | 57 µs | 62.7 µs | 531 µs |  |  | 9.3x faster |
| `bspline2/n1000/sorted` | 32.3 µs | 33.4 µs | 83.2 µs |  |  | 2.6x faster |
| `bspline2/n1000/shuffled` | 57.6 µs | 63.7 µs | 598 µs |  |  | 10x faster |
| `bspline3/n1000/sorted` | 36.6 µs | 36.9 µs | 140 µs |  |  | 3.8x faster |
| `bspline3/n1000/shuffled` | 85.6 µs | 80.2 µs | 633 µs |  |  | 7.4x faster |
| `cubic_spline/n1000/sorted` | 15.5 µs | 14.8 µs | 39.8 µs | 1.44 ms | 34.9 µs | 2.2x faster |
| `cubic_spline/n1000/shuffled` | 19.7 µs | 20.1 µs | 221 µs | 1.55 ms | 31.8 µs | 1.6x faster |
| `pchip/n1000/sorted` | 15.7 µs | 15.8 µs | 37.2 µs | 1.26 ms |  | 2.4x faster |
| `pchip/n1000/shuffled` | 19.8 µs | 20.6 µs | 221 µs | 1.23 ms |  | 11x faster |
| `akima/n1000/sorted` | 14.5 µs | 15 µs | 40 µs | 1.55 ms |  | 2.8x faster |
| `akima/n1000/shuffled` | 19.7 µs | 20.3 µs | 218 µs | 1.44 ms |  | 11x faster |
| `makima/n1000/sorted` | 15.6 µs | 15.7 µs | 37.4 µs |  |  | 2.4x faster |
| `makima/n1000/shuffled` | 19.8 µs | 20.1 µs | 228 µs |  |  | 12x faster |
| `linear/n10000/sorted` | 152 µs | 72.2 µs | 245 µs | 1.79 ms | 271 µs | 1.6x faster |
| `linear/n10000/shuffled` | 190 µs | 87.6 µs | 2.34 ms | 3.08 ms | 284 µs | 1.5x faster |
| `nearest/n10000/sorted` | 173 µs | 75.8 µs | 624 µs | 1.91 ms |  | 3.6x faster |
| `nearest/n10000/shuffled` | 430 µs | 99 µs | 3.21 ms | 3.3 ms |  | 7.5x faster |
| `previous/n10000/sorted` | 172 µs | 80.2 µs | 655 µs |  |  | 3.8x faster |
| `previous/n10000/shuffled` | 351 µs | 89.5 µs | 3.23 ms |  |  | 9.2x faster |
| `bspline1/n10000/sorted` | 266 µs | 94.4 µs | 289 µs |  |  | 1.1x faster |
| `bspline1/n10000/shuffled` | 710 µs | 144 µs | 50.2 ms |  |  | 71x faster |
| `bspline2/n10000/sorted` | 366 µs | 93.7 µs | 809 µs |  |  | 2.2x faster |
| `bspline2/n10000/shuffled` | 966 µs | 218 µs | 51.4 ms |  |  | 53x faster |
| `bspline3/n10000/sorted` | 408 µs | 107 µs | 1.3 ms |  |  | 3.2x faster |
| `bspline3/n10000/shuffled` | 1.12 ms | 252 µs | 52.3 ms |  |  | 47x faster |
| `cubic_spline/n10000/sorted` | 162 µs | 70.6 µs | 372 µs | 150 ms | 312 µs | 1.9x faster |
| `cubic_spline/n10000/shuffled` | 236 µs | 97 µs | 3.63 ms | 155 ms | 332 µs | 1.4x faster |
| `pchip/n10000/sorted` | 152 µs | 75.3 µs | 396 µs | 150 ms |  | 2.6x faster |
| `pchip/n10000/shuffled` | 220 µs | 99.4 µs | 3.35 ms | 132 ms |  | 15x faster |
| `akima/n10000/sorted` | 166 µs | 69.1 µs | 368 µs | 164 ms |  | 2.2x faster |
| `akima/n10000/shuffled` | 216 µs | 97.2 µs | 3.42 ms | 150 ms |  | 16x faster |
| `makima/n10000/sorted` | 154 µs | 69.5 µs | 401 µs |  |  | 2.6x faster |
| `makima/n10000/shuffled` | 215 µs | 97.1 µs | 4.14 ms |  |  | 19x faster |
| `linear/n100/sorted-uneven` | 1.34 µs | 1.33 µs | 8.1 µs | 3.8 µs | 4.12 µs | 2.8x faster |
| `linear/n100/shuffled-uneven` | 4.31 µs | 3.17 µs | 13.8 µs | 3.94 µs | 4.22 µs | **1.1x slower** |
| `cubic_spline/n100/sorted-uneven` | 1.35 µs | 1.46 µs | 8.28 µs | 18.8 µs | 4.85 µs | 3.6x faster |
| `cubic_spline/n100/shuffled-uneven` | 3.18 µs | 3.46 µs | 8.78 µs | 19.5 µs | 4.8 µs | 1.5x faster |
| `pchip/n100/sorted-uneven` | 1.48 µs | 1.59 µs | 8.11 µs | 17.1 µs |  | 5.5x faster |
| `pchip/n100/shuffled-uneven` | 3.31 µs | 3.32 µs | 8.69 µs | 17.5 µs |  | 2.6x faster |
| `linear/n10000/sorted-uneven` | 225 µs | 79.2 µs | 284 µs | 1.7 ms | 972 µs | 1.3x faster |
| `linear/n10000/shuffled-uneven` | 676 µs | 170 µs | 2.67 ms | 3.11 ms | 1.97 ms | 2.9x faster |
| `cubic_spline/n10000/sorted-uneven` | 233 µs | 93.8 µs | 440 µs | 159 ms | 1.05 ms | 1.9x faster |
| `cubic_spline/n10000/shuffled-uneven` | 806 µs | 186 µs | 3.69 ms | 147 ms | 1.98 ms | 2.5x faster |
| `pchip/n10000/sorted-uneven` | 219 µs | 92.5 µs | 443 µs | 128 ms |  | 2x faster |
| `pchip/n10000/shuffled-uneven` | 802 µs | 183 µs | 3.68 ms | 128 ms |  | 4.6x faster |

### Construction

| case | interpx | interpx, parallel | SciPy | scirs2-interpolate | ndarray-interp | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/n10/sorted` | 99.8 ns | 103 ns | 9.54 µs | 93.6 ns | 45.4 ns | **2.2x slower** |
| `linear/n10/shuffled` | 91.4 ns | 92.6 ns | 9.53 µs | 98.1 ns | 53.1 ns | **1.7x slower** |
| `nearest/n10/sorted` | 75 ns | 82.6 ns | 10.5 µs | 98.2 ns |  | 1.3x faster |
| `nearest/n10/shuffled` | 82.9 ns | 84.8 ns | 10.5 µs | 103 ns |  | 1.2x faster |
| `previous/n10/sorted` | 82.4 ns | 82.5 ns | 11 µs |  |  | 1.3e+02x faster |
| `previous/n10/shuffled` | 82.8 ns | 82.2 ns | 10.2 µs |  |  | 1.2e+02x faster |
| `bspline1/n10/sorted` | 112 ns | 114 ns | 15.8 µs |  |  | 1.4e+02x faster |
| `bspline1/n10/shuffled` | 112 ns | 113 ns | 17.3 µs |  |  | 1.5e+02x faster |
| `bspline2/n10/sorted` | 524 ns | 526 ns | 42 µs |  |  | 80x faster |
| `bspline2/n10/shuffled` | 496 ns | 523 ns | 37.5 µs |  |  | 76x faster |
| `bspline3/n10/sorted` | 518 ns | 561 ns | 38.9 µs |  |  | 75x faster |
| `bspline3/n10/shuffled` | 526 ns | 561 ns | 36 µs |  |  | 69x faster |
| `cubic_spline/n10/sorted` | 237 ns | 256 ns | 86.4 µs | 381 ns | 500 ns | 1.6x faster |
| `cubic_spline/n10/shuffled` | 238 ns | 244 ns | 88.1 µs | 386 ns | 532 ns | 1.6x faster |
| `pchip/n10/sorted` | 172 ns | 192 ns | 102 µs | 203 ns |  | 1.2x faster |
| `pchip/n10/shuffled` | 189 ns | 173 ns | 102 µs | 218 ns |  | 1.2x faster |
| `akima/n10/sorted` | 181 ns | 173 ns | 75.4 µs | 197 ns |  | 1.1x faster |
| `akima/n10/shuffled` | 175 ns | 178 ns | 77.6 µs | 183 ns |  | 1x faster |
| `makima/n10/sorted` | 187 ns | 202 ns | 78.5 µs |  |  | 4.2e+02x faster |
| `makima/n10/shuffled` | 184 ns | 185 ns | 84.2 µs |  |  | 4.6e+02x faster |
| `linear/n50/sorted` | 166 ns | 158 ns | 8.93 µs | 132 ns | 110 ns | **1.5x slower** |
| `linear/n50/shuffled` | 167 ns | 171 ns | 9.78 µs | 128 ns | 100 ns | **1.7x slower** |
| `nearest/n50/sorted` | 117 ns | 108 ns | 9.8 µs | 120 ns |  | 1x faster |
| `nearest/n50/shuffled` | 116 ns | 119 ns | 10.6 µs | 122 ns |  | 1.1x faster |
| `previous/n50/sorted` | 118 ns | 109 ns | 10.5 µs |  |  | 89x faster |
| `previous/n50/shuffled` | 116 ns | 119 ns | 10.5 µs |  |  | 91x faster |
| `bspline1/n50/sorted` | 157 ns | 145 ns | 15.9 µs |  |  | 1e+02x faster |
| `bspline1/n50/shuffled` | 144 ns | 154 ns | 15.9 µs |  |  | 1.1e+02x faster |
| `bspline2/n50/sorted` | 1.84 µs | 1.84 µs | 44 µs |  |  | 24x faster |
| `bspline2/n50/shuffled` | 1.85 µs | 1.87 µs | 44.2 µs |  |  | 24x faster |
| `bspline3/n50/sorted` | 2.1 µs | 2.02 µs | 39.8 µs |  |  | 19x faster |
| `bspline3/n50/shuffled` | 2.06 µs | 2.1 µs | 39.5 µs |  |  | 19x faster |
| `cubic_spline/n50/sorted` | 914 ns | 836 ns | 88.9 µs | 911 ns | 1.09 µs | **1x slower** |
| `cubic_spline/n50/shuffled` | 834 ns | 849 ns | 86.9 µs | 938 ns | 1.08 µs | 1.1x faster |
| `pchip/n50/sorted` | 479 ns | 490 ns | 102 µs | 353 ns |  | **1.4x slower** |
| `pchip/n50/shuffled` | 497 ns | 524 ns | 101 µs | 325 ns |  | **1.5x slower** |
| `akima/n50/sorted` | 535 ns | 491 ns | 75.1 µs | 325 ns |  | **1.6x slower** |
| `akima/n50/shuffled` | 535 ns | 536 ns | 78.1 µs | 354 ns |  | **1.5x slower** |
| `makima/n50/sorted` | 579 ns | 580 ns | 80.1 µs |  |  | 1.4e+02x faster |
| `makima/n50/shuffled` | 581 ns | 579 ns | 86.1 µs |  |  | 1.5e+02x faster |
| `linear/n100/sorted` | 247 ns | 270 ns | 9.99 µs | 215 ns | 182 ns | **1.4x slower** |
| `linear/n100/shuffled` | 268 ns | 270 ns | 10 µs | 215 ns | 183 ns | **1.5x slower** |
| `nearest/n100/sorted` | 175 ns | 175 ns | 11 µs | 215 ns |  | 1.2x faster |
| `nearest/n100/shuffled` | 176 ns | 173 ns | 10.9 µs | 220 ns |  | 1.3x faster |
| `previous/n100/sorted` | 176 ns | 176 ns | 10.9 µs |  |  | 62x faster |
| `previous/n100/shuffled` | 175 ns | 176 ns | 10.2 µs |  |  | 58x faster |
| `bspline1/n100/sorted` | 215 ns | 217 ns | 15.8 µs |  |  | 74x faster |
| `bspline1/n100/shuffled` | 216 ns | 203 ns | 17.2 µs |  |  | 80x faster |
| `bspline2/n100/sorted` | 3.49 µs | 3.61 µs | 48.1 µs |  |  | 14x faster |
| `bspline2/n100/shuffled` | 3.83 µs | 3.5 µs | 46.4 µs |  |  | 12x faster |
| `bspline3/n100/sorted` | 3.89 µs | 3.91 µs | 48.1 µs |  |  | 12x faster |
| `bspline3/n100/shuffled` | 3.89 µs | 4.26 µs | 44.8 µs |  |  | 12x faster |
| `cubic_spline/n100/sorted` | 1.73 µs | 1.58 µs | 87.7 µs | 1.67 µs | 1.91 µs | **1x slower** |
| `cubic_spline/n100/shuffled` | 1.58 µs | 1.73 µs | 97.5 µs | 1.8 µs | 1.92 µs | 1.1x faster |
| `pchip/n100/sorted` | 964 ns | 966 ns | 110 µs | 592 ns |  | **1.6x slower** |
| `pchip/n100/shuffled` | 961 ns | 880 ns | 109 µs | 594 ns |  | **1.6x slower** |
| `akima/n100/sorted` | 963 ns | 883 ns | 83.7 µs | 584 ns |  | **1.6x slower** |
| `akima/n100/shuffled` | 882 ns | 966 ns | 83.8 µs | 536 ns |  | **1.6x slower** |
| `makima/n100/sorted` | 1.05 µs | 1.05 µs | 79.4 µs |  |  | 76x faster |
| `makima/n100/shuffled` | 1.05 µs | 984 ns | 79.3 µs |  |  | 76x faster |
| `linear/n1000/sorted` | 1.92 µs | 2.12 µs | 12.4 µs | 1.16 µs | 1.38 µs | **1.6x slower** |
| `linear/n1000/shuffled` | 2.11 µs | 1.99 µs | 13.8 µs | 1.17 µs | 1.51 µs | **1.8x slower** |
| `nearest/n1000/sorted` | 1.08 µs | 1.01 µs | 14.4 µs | 1.28 µs |  | 1.2x faster |
| `nearest/n1000/shuffled` | 1.05 µs | 1.07 µs | 14.5 µs | 1.17 µs |  | 1.1x faster |
| `previous/n1000/sorted` | 995 ns | 992 ns | 16.2 µs |  |  | 16x faster |
| `previous/n1000/shuffled` | 1.09 µs | 1.01 µs | 17.9 µs |  |  | 16x faster |
| `bspline1/n1000/sorted` | 1.04 µs | 1.14 µs | 17 µs |  |  | 16x faster |
| `bspline1/n1000/shuffled` | 1.04 µs | 1.12 µs | 16.9 µs |  |  | 16x faster |
| `bspline2/n1000/sorted` | 33.4 µs | 36.7 µs | 119 µs |  |  | 3.6x faster |
| `bspline2/n1000/shuffled` | 33.4 µs | 33.4 µs | 110 µs |  |  | 3.3x faster |
| `bspline3/n1000/sorted` | 38.4 µs | 40 µs | 138 µs |  |  | 3.6x faster |
| `bspline3/n1000/shuffled` | 41.9 µs | 43.1 µs | 127 µs |  |  | 3x faster |
| `cubic_spline/n1000/sorted` | 16.6 µs | 15.2 µs | 112 µs | 15.9 µs | 17.2 µs | **1x slower** |
| `cubic_spline/n1000/shuffled` | 15.2 µs | 15.2 µs | 114 µs | 14.5 µs | 15.7 µs | **1x slower** |
| `pchip/n1000/sorted` | 8.28 µs | 9.06 µs | 125 µs | 4.29 µs |  | **1.9x slower** |
| `pchip/n1000/shuffled` | 8.27 µs | 9.03 µs | 120 µs | 4.56 µs |  | **1.8x slower** |
| `akima/n1000/sorted` | 8.5 µs | 8.79 µs | 95.9 µs | 4.41 µs |  | **1.9x slower** |
| `akima/n1000/shuffled` | 8.77 µs | 8.78 µs | 93.4 µs | 4.59 µs |  | **1.9x slower** |
| `makima/n1000/sorted` | 9.19 µs | 9.7 µs | 104 µs |  |  | 11x faster |
| `makima/n1000/shuffled` | 9.34 µs | 9.26 µs | 101 µs |  |  | 11x faster |
| `linear/n10000/sorted` | 21.4 µs | 21.5 µs | 50.9 µs | 12.9 µs | 14.8 µs | **1.7x slower** |
| `linear/n10000/shuffled` | 21.3 µs | 21.7 µs | 48.5 µs | 12.8 µs | 20.1 µs | **1.7x slower** |
| `nearest/n10000/sorted` | 10.7 µs | 11.5 µs | 57.5 µs | 11.7 µs |  | 1.1x faster |
| `nearest/n10000/shuffled` | 10.9 µs | 11.7 µs | 58.2 µs | 12.5 µs |  | 1.1x faster |
| `previous/n10000/sorted` | 11.7 µs | 10.9 µs | 74.1 µs |  |  | 6.3x faster |
| `previous/n10000/shuffled` | 11.6 µs | 11.7 µs | 74.1 µs |  |  | 6.4x faster |
| `bspline1/n10000/sorted` | 10.6 µs | 11.6 µs | 23.1 µs |  |  | 2.2x faster |
| `bspline1/n10000/shuffled` | 10.6 µs | 11.9 µs | 23 µs |  |  | 2.2x faster |
| `bspline2/n10000/sorted` | 518 µs | 525 µs | 766 µs |  |  | 1.5x faster |
| `bspline2/n10000/shuffled` | 514 µs | 490 µs | 822 µs |  |  | 1.6x faster |
| `bspline3/n10000/sorted` | 524 µs | 562 µs | 962 µs |  |  | 1.8x faster |
| `bspline3/n10000/shuffled` | 448 µs | 566 µs | 1.01 ms |  |  | 2.3x faster |
| `cubic_spline/n10000/sorted` | 167 µs | 184 µs | 609 µs | 145 µs | 167 µs | **1.2x slower** |
| `cubic_spline/n10000/shuffled` | 167 µs | 168 µs | 638 µs | 149 µs | 167 µs | **1.1x slower** |
| `pchip/n10000/sorted` | 91.8 µs | 193 µs | 257 µs | 42.7 µs |  | **2.1x slower** |
| `pchip/n10000/shuffled` | 124 µs | 94.5 µs | 256 µs | 46.7 µs |  | **2.7x slower** |
| `akima/n10000/sorted` | 126 µs | 95.2 µs | 665 µs | 49.7 µs |  | **2.5x slower** |
| `akima/n10000/shuffled` | 95.1 µs | 138 µs | 231 µs | 43.4 µs |  | **2.2x slower** |
| `makima/n10000/sorted` | 112 µs | 206 µs | 752 µs |  |  | 6.7x faster |
| `makima/n10000/shuffled` | 136 µs | 176 µs | 735 µs |  |  | 5.4x faster |
| `linear/n100/sorted-uneven` | 251 ns | 247 ns | 9.97 µs | 217 ns | 166 ns | **1.5x slower** |
| `linear/n100/shuffled-uneven` | 273 ns | 272 ns | 9.05 µs | 215 ns | 183 ns | **1.5x slower** |
| `cubic_spline/n100/sorted-uneven` | 1.61 µs | 1.74 µs | 99.7 µs | 1.8 µs | 2.05 µs | 1.1x faster |
| `cubic_spline/n100/shuffled-uneven` | 1.65 µs | 1.58 µs | 95.9 µs | 1.8 µs | 2.05 µs | 1.1x faster |
| `pchip/n100/sorted-uneven` | 965 ns | 967 ns | 112 µs | 590 ns |  | **1.6x slower** |
| `pchip/n100/shuffled-uneven` | 966 ns | 966 ns | 112 µs | 595 ns |  | **1.6x slower** |
| `linear/n10000/sorted-uneven` | 22.9 µs | 21.1 µs | 50.9 µs | 12.9 µs | 15 µs | **1.8x slower** |
| `linear/n10000/shuffled-uneven` | 21.4 µs | 21.7 µs | 49.5 µs | 12.8 µs | 16.6 µs | **1.7x slower** |
| `cubic_spline/n10000/sorted-uneven` | 153 µs | 269 µs | 718 µs | 159 µs | 183 µs | 1x faster |
| `cubic_spline/n10000/shuffled-uneven` | 176 µs | 174 µs | 692 µs | 145 µs | 183 µs | **1.2x slower** |
| `pchip/n10000/sorted-uneven` | 123 µs | 96.9 µs | 244 µs | 46.9 µs |  | **2.6x slower** |
| `pchip/n10000/shuffled-uneven` | 93.6 µs | 94.5 µs | 258 µs | 46.9 µs |  | **2x slower** |

## `PPoly`

A cubic of 55 intervals, evaluated at points in no order.

### Evaluation

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `intervals55/eval10` | 79 ns | 81.4 ns | 3.64 µs | 46x faster |
| `intervals55/eval100` | 667 ns | 741 ns | 4.4 µs | 6.6x faster |
| `intervals55/eval1000` | 6.52 µs | 7.28 µs | 15.2 µs | 2.3x faster |
| `intervals55/eval10000` | 66.1 µs | 72.8 µs | 360 µs | 5.5x faster |

### Construction

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `intervals55/eval10` | 173 ns | 163 ns | 5.12 µs | 30x faster |
| `intervals55/eval100` | 159 ns | 159 ns | 5.14 µs | 32x faster |
| `intervals55/eval1000` | 158 ns | 173 ns | 5.14 µs | 33x faster |
| `intervals55/eval10000` | 162 ns | 173 ns | 5.18 µs | 32x faster |

## `BSpline`

A cubic B-spline of 55 intervals, evaluated at points in no order.

### Evaluation

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `intervals55/eval10` | 144 ns | 159 ns | 2.7 µs | 19x faster |
| `intervals55/eval100` | 1.4 µs | 1.62 µs | 5.84 µs | 4.2x faster |
| `intervals55/eval1000` | 14 µs | 15.5 µs | 38.7 µs | 2.8x faster |
| `intervals55/eval10000` | 140 µs | 59.8 µs | 393 µs | 2.8x faster |

### Construction

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `intervals55/eval10` | 161 ns | 162 ns | 9.15 µs | 57x faster |
| `intervals55/eval100` | 161 ns | 175 ns | 9.13 µs | 57x faster |
| `intervals55/eval1000` | 161 ns | 169 ns | 9.95 µs | 62x faster |
| `intervals55/eval10000` | 163 ns | 178 ns | 9.19 µs | 56x faster |

## `BSplineCurve`

A cubic curve through `n` points of the plane, evaluated at `4 n` parameters in order.

### Evaluation

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `degree3/n100` | 4.54 µs | 4.69 µs | 17.3 µs | 3.8x faster |
| `degree3/n10000` | 495 µs | 176 µs | 1.47 ms | 3x faster |

### Construction

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `degree3/n100` | 4.65 µs | 4.64 µs | 47.7 µs | 10x faster |
| `degree3/n10000` | 700 µs | 692 µs | 1.5 ms | 2.1x faster |

## Regular grid

Points in no order inside the grid. `maxN`: the first axis has `N` points and each axis half the one before (SciPy's suite); `eachN`: every axis has `N` points; `flipped`: every axis given from its high end.

### Evaluation

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | ndarray-interp | interpn | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/2d/max10/eval10/ascending` | 115 ns | 107 ns | 115 ns | 21.3 µs | 865 ns | 167 ns | 95.1 ns | **1.2x slower** |
| `linear/2d/max10/eval10/flipped` | 106 ns | 117 ns | 115 ns | 19.2 µs |  |  |  | 1.8e+02x faster |
| `linear/2d/max10/eval100/ascending` | 1.04 µs | 1.17 µs | 1.15 µs | 21.7 µs | 8.27 µs | 1.37 µs | 1.02 µs | **1x slower** |
| `linear/2d/max10/eval100/flipped` | 1.04 µs | 1.22 µs | 1.09 µs | 21.8 µs |  |  |  | 21x faster |
| `linear/2d/max10/eval1000/ascending` | 11.5 µs | 10.7 µs | 10.4 µs | 43.9 µs | 90.5 µs | 14.1 µs | 9.25 µs | **1.2x slower** |
| `linear/2d/max10/eval1000/flipped` | 11.4 µs | 10.6 µs | 11.4 µs | 40.2 µs |  |  |  | 3.5x faster |
| `linear/2d/max10/eval10000/ascending` | 104 µs | 55.1 µs | 104 µs | 370 µs | 865 µs | 221 µs | 92.6 µs | **1.1x slower** |
| `linear/2d/max10/eval10000/flipped` | 114 µs | 45.1 µs | 104 µs | 390 µs |  |  |  | 3.4x faster |
| `linear/2d/max40/eval10/ascending` | 150 ns | 152 ns | 151 ns | 19.4 µs | 972 ns | 157 ns | 147 ns | **1x slower** |
| `linear/2d/max40/eval10/flipped` | 150 ns | 151 ns | 150 ns | 21.1 µs |  |  |  | 1.4e+02x faster |
| `linear/2d/max40/eval100/ascending` | 1.5 µs | 1.54 µs | 1.49 µs | 22.2 µs | 9.76 µs | 1.59 µs | 1.45 µs | **1x slower** |
| `linear/2d/max40/eval100/flipped` | 1.37 µs | 1.57 µs | 1.49 µs | 21.8 µs |  |  |  | 16x faster |
| `linear/2d/max40/eval1000/ascending` | 15 µs | 14.1 µs | 15 µs | 49.5 µs | 108 µs | 16.5 µs | 14.4 µs | **1x slower** |
| `linear/2d/max40/eval1000/flipped` | 14.9 µs | 15.1 µs | 14.9 µs | 50.4 µs |  |  |  | 3.4x faster |
| `linear/2d/max40/eval10000/ascending` | 136 µs | 58.2 µs | 136 µs | 569 µs | 1.05 ms | 223 µs | 131 µs | **1x slower** |
| `linear/2d/max40/eval10000/flipped` | 149 µs | 61.8 µs | 136 µs | 575 µs |  |  |  | 3.8x faster |
| `linear/2d/max200/eval10/ascending` | 181 ns | 198 ns | 197 ns | 19.5 µs | 1.91 µs | 206 ns | 192 ns | 1.1x faster |
| `linear/2d/max200/eval10/flipped` | 180 ns | 201 ns | 198 ns | 21.2 µs |  |  |  | 1.2e+02x faster |
| `linear/2d/max200/eval100/ascending` | 1.9 µs | 1.99 µs | 1.91 µs | 24.8 µs | 21.2 µs | 2.18 µs | 1.89 µs | **1x slower** |
| `linear/2d/max200/eval100/flipped` | 1.74 µs | 1.98 µs | 1.91 µs | 25.1 µs |  |  |  | 14x faster |
| `linear/2d/max200/eval1000/ascending` | 17.6 µs | 17.7 µs | 19.2 µs | 79.2 µs | 211 µs | 24.4 µs | 18.9 µs | 1.1x faster |
| `linear/2d/max200/eval1000/flipped` | 17.5 µs | 17.9 µs | 19.1 µs | 73.4 µs |  |  |  | 4.2x faster |
| `linear/2d/max200/eval10000/ascending` | 191 µs | 72.3 µs | 176 µs | 828 µs | 2.28 ms | 496 µs | 194 µs | 1x faster |
| `linear/2d/max200/eval10000/flipped` | 174 µs | 86.8 µs | 176 µs | 841 µs |  |  |  | 4.8x faster |
| `linear/3d/max10/eval10/ascending` | 271 ns | 165 ns | 168 ns | 62.2 µs | 3.46 µs |  | 156 ns | **1.7x slower** |
| `linear/3d/max10/eval10/flipped` | 164 ns | 152 ns | 164 ns | 62.3 µs |  |  |  | 3.8e+02x faster |
| `linear/3d/max10/eval100/ascending` | 1.48 µs | 1.7 µs | 1.48 µs | 70.2 µs | 34.2 µs |  | 1.39 µs | **1.1x slower** |
| `linear/3d/max10/eval100/flipped` | 1.5 µs | 1.7 µs | 1.63 µs | 63.6 µs |  |  |  | 42x faster |
| `linear/3d/max10/eval1000/ascending` | 14.8 µs | 14.9 µs | 14.8 µs | 118 µs | 345 µs |  | 15.3 µs | 1x faster |
| `linear/3d/max10/eval1000/flipped` | 14.7 µs | 16.4 µs | 15.1 µs | 118 µs |  |  |  | 8x faster |
| `linear/3d/max10/eval10000/ascending` | 162 µs | 62.7 µs | 162 µs | 1.15 ms | 3.28 ms |  | 156 µs | **1x slower** |
| `linear/3d/max10/eval10000/flipped` | 148 µs | 59.7 µs | 149 µs | 1.19 ms |  |  |  | 8x faster |
| `linear/3d/max40/eval10/ascending` | 210 ns | 206 ns | 205 ns | 56.1 µs | 3.33 µs |  | 218 ns | 1x faster |
| `linear/3d/max40/eval10/flipped` | 205 ns | 206 ns | 205 ns | 58.2 µs |  |  |  | 2.8e+02x faster |
| `linear/3d/max40/eval100/ascending` | 1.98 µs | 2.04 µs | 2.17 µs | 62.9 µs | 32.6 µs |  | 2 µs | 1x faster |
| `linear/3d/max40/eval100/flipped` | 1.98 µs | 2.04 µs | 1.98 µs | 64.4 µs |  |  |  | 32x faster |
| `linear/3d/max40/eval1000/ascending` | 19.8 µs | 21.8 µs | 21.7 µs | 139 µs | 334 µs |  | 21.9 µs | 1.1x faster |
| `linear/3d/max40/eval1000/flipped` | 19.8 µs | 20.1 µs | 19.8 µs | 150 µs |  |  |  | 7.6x faster |
| `linear/3d/max40/eval10000/ascending` | 199 µs | 69 µs | 198 µs | 1.53 ms | 3.72 ms |  | 200 µs | 1x faster |
| `linear/3d/max40/eval10000/flipped` | 199 µs | 66.1 µs | 199 µs | 1.55 ms |  |  |  | 7.8x faster |
| `linear/3d/max200/eval10/ascending` | 259 ns | 286 ns | 284 ns | 58.6 µs | 4.47 µs |  | 263 ns | 1x faster |
| `linear/3d/max200/eval10/flipped` | 259 ns | 286 ns | 266 ns | 58.3 µs |  |  |  | 2.2e+02x faster |
| `linear/3d/max200/eval100/ascending` | 2.61 µs | 2.62 µs | 2.53 µs | 65.5 µs | 46.8 µs |  | 2.61 µs | 1x faster |
| `linear/3d/max200/eval100/flipped` | 2.77 µs | 2.84 µs | 2.53 µs | 65.3 µs |  |  |  | 24x faster |
| `linear/3d/max200/eval1000/ascending` | 25.3 µs | 28.2 µs | 27.8 µs | 191 µs | 459 µs |  | 26.3 µs | 1x faster |
| `linear/3d/max200/eval1000/flipped` | 27.7 µs | 27.9 µs | 25.3 µs | 185 µs |  |  |  | 6.7x faster |
| `linear/3d/max200/eval10000/ascending` | 277 µs | 85.1 µs | 253 µs | 1.95 ms | 5.13 ms |  | 271 µs | **1x slower** |
| `linear/3d/max200/eval10000/flipped` | 277 µs | 84 µs | 277 µs | 1.98 ms |  |  |  | 7.1x faster |
| `linear/4d/max10/eval10/ascending` | 228 ns | 245 ns | 252 ns | 97.4 µs |  |  |  | 4.3e+02x faster |
| `linear/4d/max10/eval10/flipped` | 235 ns | 248 ns | 240 ns | 99.5 µs |  |  |  | 4.2e+02x faster |
| `linear/4d/max10/eval100/ascending` | 2.11 µs | 2.22 µs | 2.32 µs | 116 µs |  |  |  | 55x faster |
| `linear/4d/max10/eval100/flipped` | 2.14 µs | 2.39 µs | 2.32 µs | 106 µs |  |  |  | 49x faster |
| `linear/4d/max10/eval1000/ascending` | 23.2 µs | 21.5 µs | 23.1 µs | 224 µs |  |  |  | 9.6x faster |
| `linear/4d/max10/eval1000/flipped` | 23.2 µs | 23.1 µs | 21.2 µs | 205 µs |  |  |  | 8.9x faster |
| `linear/4d/max10/eval10000/ascending` | 213 µs | 71.5 µs | 211 µs | 1.92 ms |  |  |  | 9x faster |
| `linear/4d/max10/eval10000/flipped` | 212 µs | 72.7 µs | 211 µs | 1.87 ms |  |  |  | 8.8x faster |
| `linear/4d/max40/eval10/ascending` | 291 ns | 291 ns | 320 ns | 94.8 µs | 5.84 µs |  | 297 ns | 1x faster |
| `linear/4d/max40/eval10/flipped` | 291 ns | 292 ns | 318 ns | 95.8 µs |  |  |  | 3.3e+02x faster |
| `linear/4d/max40/eval100/ascending` | 2.8 µs | 2.84 µs | 3.02 µs | 107 µs | 62.5 µs |  | 2.95 µs | 1.1x faster |
| `linear/4d/max40/eval100/flipped` | 2.79 µs | 2.85 µs | 2.75 µs | 119 µs |  |  |  | 43x faster |
| `linear/4d/max40/eval1000/ascending` | 27.4 µs | 28 µs | 27.4 µs | 267 µs | 610 µs |  | 26.8 µs | **1x slower** |
| `linear/4d/max40/eval1000/flipped` | 27.5 µs | 28.1 µs | 30.2 µs | 253 µs |  |  |  | 9.2x faster |
| `linear/4d/max40/eval10000/ascending` | 282 µs | 87 µs | 301 µs | 2.24 ms | 6.36 ms |  | 270 µs | **1x slower** |
| `linear/4d/max40/eval10000/flipped` | 276 µs | 87.6 µs | 275 µs | 2.13 ms |  |  |  | 7.7x faster |
| `linear/4d/max200/eval10/ascending` | 365 ns | 368 ns | 360 ns | 95.5 µs | 7.87 µs |  | 362 ns | **1x slower** |
| `linear/4d/max200/eval10/flipped` | 398 ns | 364 ns | 363 ns | 97.2 µs |  |  |  | 2.4e+02x faster |
| `linear/4d/max200/eval100/ascending` | 3.65 µs | 4.09 µs | 3.62 µs | 109 µs | 71 µs |  | 3.64 µs | **1x slower** |
| `linear/4d/max200/eval100/flipped` | 3.96 µs | 4.09 µs | 3.92 µs | 119 µs |  |  |  | 30x faster |
| `linear/4d/max200/eval1000/ascending` | 36.6 µs | 40 µs | 36.4 µs | 342 µs | 758 µs |  | 36.4 µs | **1x slower** |
| `linear/4d/max200/eval1000/flipped` | 36.6 µs | 40.6 µs | 39.8 µs | 309 µs |  |  |  | 8.4x faster |
| `linear/4d/max200/eval10000/ascending` | 365 µs | 110 µs | 379 µs | 2.82 ms | 7.56 ms |  | 364 µs | **1x slower** |
| `linear/4d/max200/eval10000/flipped` | 399 µs | 112 µs | 400 µs | 3.02 ms |  |  |  | 7.6x faster |
| `nearest/2d/each40/eval10000/ascending` | 142 µs | 59.1 µs | 138 µs | 653 µs | 1.67 ms |  | 148 µs | 1x faster |
| `nearest/3d/each40/eval10000/ascending` | 195 µs | 72 µs | 199 µs | 954 µs | 1.97 ms |  | 208 µs | 1.1x faster |
| `pchip/2d/each40/eval10000/ascending` | 477 µs | 132 µs | 473 µs | 1.32 s |  |  |  | 2.8e+03x faster |
| `pchip/3d/each40/eval10000/ascending` | 1.58 ms | 274 µs | 1.45 ms | 3.47 s |  |  |  | 2.2e+03x faster |
| `cubic/2d/each40/eval10000/ascending` | 287 µs | 99.6 µs | 301 µs | 1.3 ms |  |  | 409 µs | 1.4x faster |
| `cubic/3d/each40/eval10000/ascending` | 654 µs | 177 µs | 658 µs | 3.03 ms |  |  | 790 µs | 1.2x faster |
| `cubic/2d/each10/eval1000/ascending` | 25 µs | 25.2 µs | 24.9 µs | 143 µs |  |  | 35.8 µs | 1.4x faster |
| `cubic/2d/each40/eval1000/ascending` | 28.6 µs | 31.5 µs | 28.5 µs | 143 µs |  |  | 40.7 µs | 1.4x faster |
| `cubic/2d/each100/eval1000/ascending` | 30.9 µs | 31.5 µs | 33.9 µs | 166 µs |  |  | 43 µs | 1.4x faster |
| `cubic/2d/each200/eval1000/ascending` | 34.5 µs | 34.4 µs | 34.3 µs | 205 µs |  |  | 45.7 µs | 1.3x faster |
| `cubic/2d/each400/eval1000/ascending` | 37.8 µs | 37.2 µs | 37.3 µs | 288 µs |  |  | 48.5 µs | 1.3x faster |
| `quintic/2d/each10/eval1000/ascending` | 42.9 µs | 44.1 µs | 43.8 µs | 228 µs |  |  |  | 5.3x faster |
| `quintic/2d/each40/eval1000/ascending` | 48.1 µs | 53 µs | 48.1 µs | 237 µs |  |  |  | 4.9x faster |
| `linear/2d/each400/eval10000/ascending` | 204 µs | 75.9 µs | 204 µs | 965 µs | 4.11 ms | 542 µs | 220 µs | 1.1x faster |
| `linear/2d/each800/eval10000/ascending` | 223 µs | 77.1 µs | 231 µs | 1.09 ms | 7.36 ms | 584 µs | 242 µs | 1.1x faster |
| `nearest/2d/auckland-ba32-0401/heldout/ascending` | 418 ms | 66.2 ms | 366 ms | 3.57 s | 203 s |  | 787 ms | 1.9x faster |
| `linear/2d/auckland-ba32-0401/heldout/ascending` | 376 ms | 65.8 ms | 393 ms | 1.1 s | 71.9 s | 357 ms | 820 ms | **1.1x slower** |
| `cubic/2d/auckland-ba32-0401/heldout/ascending` | 775 ms | 126 ms | 748 ms | 35.4 s |  |  | 1.7 s | 2.2x faster |
| `quintic/2d/auckland-ba32-0401/heldout/ascending` | 1.27 s | 244 ms | 1.28 s | 49.3 s |  |  |  | 39x faster |
| `pchip/2d/auckland-ba32-0401/heldout/ascending` | 1.03 s | 181 ms | 1.08 s |  |  |  |  |  |

### Construction

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | ndarray-interp | interpn | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/2d/max10/eval10/ascending` | 134 ns | 134 ns | 140 ns | 6.42 µs | 149 ns | 85.6 ns | 0.68 ns | **2e+02x slower** |
| `linear/2d/max10/eval10/flipped` | 176 ns | 161 ns | 183 ns | 11.6 µs |  |  |  | 66x faster |
| `linear/2d/max10/eval100/ascending` | 122 ns | 122 ns | 151 ns | 5.89 µs | 133 ns | 78.3 ns | 0.68 ns | **1.8e+02x slower** |
| `linear/2d/max10/eval100/flipped` | 177 ns | 177 ns | 172 ns | 11.7 µs |  |  |  | 66x faster |
| `linear/2d/max10/eval1000/ascending` | 136 ns | 130 ns | 134 ns | 6.3 µs | 146 ns | 85.4 ns | 0.68 ns | **2e+02x slower** |
| `linear/2d/max10/eval1000/flipped` | 177 ns | 163 ns | 182 ns | 11.6 µs |  |  |  | 66x faster |
| `linear/2d/max10/eval10000/ascending` | 121 ns | 131 ns | 127 ns | 5.94 µs | 145 ns | 77.9 ns | 0.68 ns | **1.8e+02x slower** |
| `linear/2d/max10/eval10000/flipped` | 175 ns | 175 ns | 180 ns | 12.9 µs |  |  |  | 73x faster |
| `linear/2d/max40/eval10/ascending` | 366 ns | 390 ns | 854 ns | 5.89 µs | 203 ns | 153 ns | 0.68 ns | **5.3e+02x slower** |
| `linear/2d/max40/eval10/flipped` | 682 ns | 708 ns | 1.2 µs | 11.8 µs |  |  |  | 17x faster |
| `linear/2d/max40/eval100/ascending` | 394 ns | 387 ns | 804 ns | 5.85 µs | 199 ns | 165 ns | 0.68 ns | **5.8e+02x slower** |
| `linear/2d/max40/eval100/flipped` | 643 ns | 705 ns | 1.23 µs | 11.7 µs |  |  |  | 18x faster |
| `linear/2d/max40/eval1000/ascending` | 391 ns | 391 ns | 822 ns | 6.36 µs | 220 ns | 153 ns | 0.68 ns | **5.7e+02x slower** |
| `linear/2d/max40/eval1000/flipped` | 643 ns | 707 ns | 1.18 µs | 12.8 µs |  |  |  | 20x faster |
| `linear/2d/max40/eval10000/ascending` | 358 ns | 389 ns | 805 ns | 6.15 µs | 202 ns | 152 ns | 0.68 ns | **5.2e+02x slower** |
| `linear/2d/max40/eval10000/flipped` | 699 ns | 645 ns | 1.18 µs | 11.4 µs |  |  |  | 16x faster |
| `linear/2d/max200/eval10/ascending` | 6.79 µs | 6.96 µs | 18.8 µs | 5.89 µs | 3.09 µs | 470 ns | 0.68 ns | **9.9e+03x slower** |
| `linear/2d/max200/eval10/flipped` | 14.7 µs | 16 µs | 25.4 µs | 12.9 µs |  |  |  | **1.1x slower** |
| `linear/2d/max200/eval100/ascending` | 6.83 µs | 7.12 µs | 19.3 µs | 6.26 µs | 3.18 µs | 515 ns | 0.78 ns | **8.7e+03x slower** |
| `linear/2d/max200/eval100/flipped` | 14.3 µs | 16 µs | 29.6 µs | 13.1 µs |  |  |  | **1.1x slower** |
| `linear/2d/max200/eval1000/ascending` | 6.82 µs | 6.95 µs | 20.1 µs | 6.39 µs | 3.24 µs | 471 ns | 0.78 ns | **8.7e+03x slower** |
| `linear/2d/max200/eval1000/flipped` | 14.2 µs | 14.6 µs | 27.3 µs | 12.9 µs |  |  |  | **1.1x slower** |
| `linear/2d/max200/eval10000/ascending` | 7.54 µs | 18.4 µs | 10.9 µs | 5.82 µs | 3.36 µs | 516 ns | 0.68 ns | **1.1e+04x slower** |
| `linear/2d/max200/eval10000/flipped` | 14.5 µs | 22.2 µs | 22.1 µs | 13 µs |  |  |  | **1.1x slower** |
| `linear/3d/max10/eval10/ascending` | 195 ns | 191 ns | 234 ns | 8.51 µs | 180 ns |  | 0.88 ns | **2.2e+02x slower** |
| `linear/3d/max10/eval10/flipped` | 345 ns | 328 ns | 360 ns | 17.4 µs |  |  |  | 50x faster |
| `linear/3d/max10/eval100/ascending` | 178 ns | 193 ns | 263 ns | 7.54 µs | 180 ns |  | 0.98 ns | **1.8e+02x slower** |
| `linear/3d/max10/eval100/flipped` | 302 ns | 330 ns | 407 ns | 17.2 µs |  |  |  | 57x faster |
| `linear/3d/max10/eval1000/ascending` | 211 ns | 180 ns | 213 ns | 7.56 µs | 181 ns |  | 0.98 ns | **2.2e+02x slower** |
| `linear/3d/max10/eval1000/flipped` | 301 ns | 305 ns | 343 ns | 15.6 µs |  |  |  | 52x faster |
| `linear/3d/max10/eval10000/ascending` | 190 ns | 198 ns | 220 ns | 7.69 µs | 165 ns |  | 0.98 ns | **1.9e+02x slower** |
| `linear/3d/max10/eval10000/flipped` | 302 ns | 324 ns | 402 ns | 15.8 µs |  |  |  | 52x faster |
| `linear/3d/max40/eval10/ascending` | 2.75 µs | 2.91 µs | 7.54 µs | 7.95 µs | 1.43 µs |  | 0.98 ns | **2.8e+03x slower** |
| `linear/3d/max40/eval10/flipped` | 8.79 µs | 8.14 µs | 12.3 µs | 15.4 µs |  |  |  | 1.7x faster |
| `linear/3d/max40/eval100/ascending` | 2.83 µs | 2.86 µs | 7.96 µs | 7.69 µs | 1.47 µs |  | 0.98 ns | **2.9e+03x slower** |
| `linear/3d/max40/eval100/flipped` | 8.19 µs | 8.38 µs | 12.6 µs | 15.7 µs |  |  |  | 1.9x faster |
| `linear/3d/max40/eval1000/ascending` | 2.92 µs | 2.93 µs | 4.29 µs | 7.73 µs | 1.47 µs |  | 0.88 ns | **3.3e+03x slower** |
| `linear/3d/max40/eval1000/flipped` | 8.08 µs | 8.16 µs | 9.43 µs | 16.1 µs |  |  |  | 2x faster |
| `linear/3d/max40/eval10000/ascending` | 2.86 µs | 2.88 µs | 3.33 µs | 8.4 µs | 1.56 µs |  | 0.88 ns | **3.3e+03x slower** |
| `linear/3d/max40/eval10000/flipped` | 7.94 µs | 8.1 µs | 10.1 µs | 16.8 µs |  |  |  | 2.1x faster |
| `linear/3d/max200/eval10/ascending` | 1.6 ms | 1.74 ms | 659 µs | 8.35 µs | 1.28 ms |  | 0.98 ns | **1.6e+06x slower** |
| `linear/3d/max200/eval10/flipped` | 2.99 ms | 3.26 ms | 1.88 ms | 16.2 µs |  |  |  | **1.8e+02x slower** |
| `linear/3d/max200/eval100/ascending` | 1.65 ms | 1.76 ms | 618 µs | 7.53 µs | 1.3 ms |  | 0.98 ns | **1.7e+06x slower** |
| `linear/3d/max200/eval100/flipped` | 3.12 ms | 3.31 ms | 1.97 ms | 17.2 µs |  |  |  | **1.8e+02x slower** |
| `linear/3d/max200/eval1000/ascending` | 1.64 ms | 1.94 ms | 657 µs | 8.45 µs | 1.19 ms |  | 0.88 ns | **1.9e+06x slower** |
| `linear/3d/max200/eval1000/flipped` | 3.23 ms | 3.15 ms | 1.81 ms | 16 µs |  |  |  | **2e+02x slower** |
| `linear/3d/max200/eval10000/ascending` | 1.76 ms | 1.81 ms | 619 µs | 8.51 µs | 1.3 ms |  | 0.98 ns | **1.8e+06x slower** |
| `linear/3d/max200/eval10000/flipped` | 3.05 ms | 3.21 ms | 1.86 ms | 17.4 µs |  |  |  | **1.8e+02x slower** |
| `linear/4d/max10/eval10/ascending` | 213 ns | 231 ns | 253 ns | 10.6 µs |  |  |  | 50x faster |
| `linear/4d/max10/eval10/flipped` | 376 ns | 372 ns | 358 ns | 19.7 µs |  |  |  | 53x faster |
| `linear/4d/max10/eval100/ascending` | 212 ns | 229 ns | 248 ns | 10.7 µs |  |  |  | 50x faster |
| `linear/4d/max10/eval100/flipped` | 372 ns | 370 ns | 378 ns | 17.9 µs |  |  |  | 48x faster |
| `linear/4d/max10/eval1000/ascending` | 212 ns | 233 ns | 238 ns | 10.6 µs |  |  |  | 50x faster |
| `linear/4d/max10/eval1000/flipped` | 339 ns | 346 ns | 357 ns | 17.8 µs |  |  |  | 53x faster |
| `linear/4d/max10/eval10000/ascending` | 220 ns | 235 ns | 230 ns | 10.1 µs |  |  |  | 46x faster |
| `linear/4d/max10/eval10000/flipped` | 340 ns | 346 ns | 353 ns | 17.9 µs |  |  |  | 53x faster |
| `linear/4d/max40/eval10/ascending` | 13.1 µs | 13.2 µs | 32.6 µs | 9.31 µs | 5.97 µs |  | 1.37 ns | **9.5e+03x slower** |
| `linear/4d/max40/eval10/flipped` | 51.3 µs | 57 µs | 66.8 µs | 20 µs |  |  |  | **2.6x slower** |
| `linear/4d/max40/eval100/ascending` | 13 µs | 13.3 µs | 35.1 µs | 9.31 µs | 5.96 µs |  | 1.27 ns | **1e+04x slower** |
| `linear/4d/max40/eval100/flipped` | 52.4 µs | 54.3 µs | 70.3 µs | 21.7 µs |  |  |  | **2.4x slower** |
| `linear/4d/max40/eval1000/ascending` | 13.2 µs | 13 µs | 35.4 µs | 9.57 µs | 5.63 µs |  | 1.27 ns | **1e+04x slower** |
| `linear/4d/max40/eval1000/flipped` | 51.3 µs | 56.4 µs | 65.8 µs | 19.8 µs |  |  |  | **2.6x slower** |
| `linear/4d/max40/eval10000/ascending` | 13.1 µs | 13.5 µs | 32.4 µs | 9.31 µs | 5.99 µs |  | 1.27 ns | **1e+04x slower** |
| `linear/4d/max40/eval10000/flipped` | 51.2 µs | 53.1 µs | 77.7 µs | 19.7 µs |  |  |  | **2.6x slower** |
| `linear/4d/max200/eval10/ascending` | 49.7 ms | 48.6 ms | 15.9 ms | 9.37 µs | 38.6 ms |  | 1.27 ns | **3.9e+07x slower** |
| `linear/4d/max200/eval10/flipped` | 101 ms | 105 ms | 69 ms | 21.7 µs |  |  |  | **4.6e+03x slower** |
| `linear/4d/max200/eval100/ascending` | 49.6 ms | 48.2 ms | 16 ms | 9.4 µs | 41.4 ms |  | 1.37 ns | **3.6e+07x slower** |
| `linear/4d/max200/eval100/flipped` | 119 ms | 107 ms | 72.6 ms | 20.5 µs |  |  |  | **5.8e+03x slower** |
| `linear/4d/max200/eval1000/ascending` | 49.6 ms | 49.7 ms | 15.1 ms | 10.2 µs | 39 ms |  | 1.27 ns | **3.9e+07x slower** |
| `linear/4d/max200/eval1000/flipped` | 103 ms | 102 ms | 74.3 ms | 20.1 µs |  |  |  | **5.1e+03x slower** |
| `linear/4d/max200/eval10000/ascending` | 49.8 ms | 48 ms | 16.2 ms | 9.52 µs | 39.5 ms |  | 1.27 ns | **3.9e+07x slower** |
| `linear/4d/max200/eval10000/flipped` | 103 ms | 107 ms | 75 ms | 20.3 µs |  |  |  | **5.1e+03x slower** |
| `nearest/2d/each40/eval10000/ascending` | 592 ns | 579 ns | 1.2 µs | 6.26 µs | 297 ns |  | 0.78 ns | **7.6e+02x slower** |
| `nearest/3d/each40/eval10000/ascending` | 20.8 µs | 21.7 µs | 50.4 µs | 7.58 µs | 9.83 µs |  | 0.88 ns | **2.4e+04x slower** |
| `pchip/2d/each40/eval10000/ascending` | 6.45 µs | 6.44 µs | 6.24 µs | 6.7 µs |  |  |  | 1x faster |
| `pchip/3d/each40/eval10000/ascending` | 512 µs | 506 µs | 459 µs | 8.88 µs |  |  |  | **58x slower** |
| `cubic/2d/each40/eval10000/ascending` | 8.2 µs | 8.93 µs | 7.92 µs | 1.84 ms |  |  | 0.88 ns | **9.3e+03x slower** |
| `cubic/3d/each40/eval10000/ascending` | 317 µs | 279 µs | 283 µs | 147 ms |  |  | 1.17 ns | **2.7e+05x slower** |
| `cubic/2d/each10/eval1000/ascending` | 1.47 µs | 1.47 µs | 1.54 µs | 1.09 ms |  |  | 0.88 ns | **1.7e+03x slower** |
| `cubic/2d/each40/eval1000/ascending` | 8.95 µs | 8.96 µs | 8 µs | 1.82 ms |  |  | 0.88 ns | **1e+04x slower** |
| `cubic/2d/each100/eval1000/ascending` | 38.5 µs | 42.3 µs | 37.5 µs | 7.61 ms |  |  | 0.88 ns | **4.4e+04x slower** |
| `cubic/2d/each200/eval1000/ascending` | 169 µs | 168 µs | 188 µs | 26.3 ms |  |  | 0.88 ns | **1.9e+05x slower** |
| `cubic/2d/each400/eval1000/ascending` | 1.48 ms | 1.37 ms | 1.15 ms | 145 ms |  |  | 0.88 ns | **1.7e+06x slower** |
| `quintic/2d/each10/eval1000/ascending` | 1.92 µs | 1.93 µs | 1.97 µs | 3.58 ms |  |  |  | 1.9e+03x faster |
| `quintic/2d/each40/eval1000/ascending` | 11.7 µs | 11.8 µs | 11.5 µs | 7.39 ms |  |  |  | 6.3e+02x faster |
| `linear/2d/each400/eval10000/ascending` | 322 µs | 327 µs | 99.5 µs | 6.02 µs | 276 µs | 1.14 µs | 0.68 ns | **4.7e+05x slower** |
| `linear/2d/each800/eval10000/ascending` | 1.11 ms | 1.07 ms | 412 µs | 6.47 µs | 842 µs | 2.23 µs | 0.68 ns | **1.6e+06x slower** |
| `nearest/2d/auckland-ba32-0401/heldout/ascending` | 16.2 ms | 15.6 ms | 6.04 ms | 7.76 µs | 13 ms |  | 0.78 ns | **2.1e+07x slower** |
| `linear/2d/auckland-ba32-0401/heldout/ascending` | 16.1 ms | 17 ms | 5.66 ms | 7.16 µs | 14.5 ms | 8.25 µs | 0.68 ns | **2.4e+07x slower** |
| `cubic/2d/auckland-ba32-0401/heldout/ascending` | 75.3 ms | 73.2 ms | 70.4 ms | 14.6 s |  |  | 0.88 ns | **8.6e+07x slower** |
| `quintic/2d/auckland-ba32-0401/heldout/ascending` | 99 ms | 94.3 ms | 83.4 ms | 21.8 s |  |  |  | 2.2e+02x faster |
| `pchip/2d/auckland-ba32-0401/heldout/ascending` | 76.6 ms | 74.6 ms | 65.5 ms |  |  |  |  |  |

## `NearestNDInterpolator`

`n` points of the unit square, evaluated at the pixels of a square image.

### Evaluation

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: |
| `n10/image50` | 32.1 µs | 32.2 µs | 122 µs | 433 µs | 13x faster |
| `n10/image100` | 123 µs | 50.1 µs | 444 µs | 1.72 ms | 14x faster |
| `n10/image500` | 3.02 ms | 536 µs | 11.6 ms | 45.2 ms | 15x faster |
| `n100/image50` | 123 µs | 98 µs | 188 µs | 1.09 ms | 8.8x faster |
| `n100/image100` | 456 µs | 133 µs | 728 µs | 2.44 ms | 5.3x faster |
| `n100/image500` | 9.41 ms | 1.56 ms | 18.4 ms | 55.2 ms | 5.9x faster |
| `n1000/image50` | 344 µs | 120 µs | 311 µs | 798 µs | 2.3x faster |
| `n1000/image100` | 1.14 ms | 294 µs | 1.17 ms | 2.98 ms | 2.6x faster |
| `n1000/image500` | 24.7 ms | 4.14 ms | 23.7 ms | 78.3 ms | 3.2x faster |
| `n10000/image50` | 327 µs | 112 µs | 352 µs | 1.06 ms | 3.2x faster |
| `n10000/image100` | 1.49 ms | 312 µs | 1.47 ms | 4.01 ms | 2.7x faster |
| `n10000/image500` | 33.3 ms | 5.22 ms | 32.5 ms | 95.2 ms | 2.9x faster |
| `n100000/image50` | 439 µs | 116 µs | 484 µs | 1.33 ms | 3x faster |
| `n100000/image100` | 1.86 ms | 353 µs | 1.93 ms | 5.62 ms | 3x faster |
| `n100000/image500` | 47 ms | 6.9 ms | 48.8 ms | 129 ms | 2.7x faster |
| `ghcn-tmax-20250101/heldout` | 160 µs | 74.1 µs | 159 µs | 520 µs | 3.2x faster |

### Construction

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: |
| `n10/image50` | 560 ns | 516 ns | 566 ns | 9.44 µs | 17x faster |
| `n10/image100` | 528 ns | 598 ns | 557 ns | 9.49 µs | 18x faster |
| `n10/image500` | 522 ns | 526 ns | 512 ns | 9.34 µs | 18x faster |
| `n100/image50` | 3.62 µs | 3.85 µs | 3.33 µs | 15 µs | 4.1x faster |
| `n100/image100` | 3.75 µs | 3.31 µs | 3.63 µs | 16.3 µs | 4.4x faster |
| `n100/image500` | 3.65 µs | 3.35 µs | 3.65 µs | 14.9 µs | 4.1x faster |
| `n1000/image50` | 37.9 µs | 34.1 µs | 37.5 µs | 116 µs | 3.1x faster |
| `n1000/image100` | 34.4 µs | 35.3 µs | 34.8 µs | 117 µs | 3.4x faster |
| `n1000/image500` | 35.4 µs | 38.1 µs | 36.8 µs | 128 µs | 3.6x faster |
| `n10000/image50` | 491 µs | 549 µs | 495 µs | 1.8 ms | 3.7x faster |
| `n10000/image100` | 531 µs | 519 µs | 550 µs | 1.82 ms | 3.4x faster |
| `n10000/image500` | 536 µs | 542 µs | 534 µs | 1.79 ms | 3.3x faster |
| `n100000/image50` | 8.49 ms | 8.3 ms | 8.51 ms | 25.3 ms | 3x faster |
| `n100000/image100` | 9.01 ms | 8.74 ms | 9.01 ms | 26.8 ms | 3x faster |
| `n100000/image500` | 8.46 ms | 8.34 ms | 8.42 ms | 26.4 ms | 3.1x faster |
| `ghcn-tmax-20250101/heldout` | 539 µs | 516 µs | 560 µs | 1.72 ms | 3.2x faster |

## `LinearNDInterpolator`

`n` points of the unit square, evaluated at the pixels of a square image (some outside the hull, with no value); `3d`: `n` points of the unit cube, evaluated at 10 000 points of the cube in no order.

### Evaluation

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | spade | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `n10/image50` | 16.4 µs | 16.6 µs | 16.3 µs | 30.3 µs | 57.2 µs | 106 µs | 1.8x faster |
| `n10/image100` | 67.3 µs | 67.5 µs | 61.4 µs | 105 µs | 237 µs | 419 µs | 1.6x faster |
| `n10/image500` | 1.71 ms | 346 µs | 1.59 ms | 2.82 ms | 5.99 ms | 10.3 ms | 1.7x faster |
| `n100/image50` | 45.2 µs | 44.6 µs | 41.3 µs | 73.5 µs | 1.02 ms | 213 µs | 1.6x faster |
| `n100/image100` | 199 µs | 61.5 µs | 186 µs | 254 µs | 2.96 ms | 877 µs | 1.3x faster |
| `n100/image500` | 4.76 ms | 825 µs | 4.51 ms | 5.97 ms | 64.7 ms | 16.9 ms | 1.3x faster |
| `n1000/image50` | 41.1 µs | 44.1 µs | 41.4 µs | 133 µs | 18 ms | 279 µs | 3.2x faster |
| `n1000/image100` | 215 µs | 79.7 µs | 219 µs | 390 µs | 50.6 ms | 1.12 ms | 1.8x faster |
| `n1000/image500` | 5.09 ms | 958 µs | 5.18 ms | 6.8 ms | 748 ms | 20.3 ms | 1.3x faster |
| `n10000/image50` | 49.4 µs | 54.4 µs | 50 µs | 378 µs | 293 ms | 569 µs | 7.7x faster |
| `n10000/image100` | 254 µs | 83.3 µs | 258 µs | 916 µs | 1.05 s | 1.66 ms | 3.6x faster |
| `n10000/image500` | 6.84 ms | 1.24 ms | 6.49 ms | 8.54 ms | 20.6 s | 29.1 ms | 1.2x faster |
| `n100000/image50` | 75.1 µs | 81.8 µs | 75.6 µs | 1.26 ms | 3.4 s | 3.34 ms | 17x faster |
| `n100000/image100` | 395 µs | 115 µs | 359 µs | 3.79 ms | 14.9 s | 7.04 ms | 9.6x faster |
| `n100000/image500` | 10.9 ms | 1.75 ms | 10.1 ms | 20.4 ms | 350 s | 48 ms | 1.9x faster |
| `3d/n1000/scattered` | 778 µs | 190 µs | 787 µs | 5.82 ms |  |  | 7.5x faster |
| `3d/n10000/scattered` | 1.45 ms | 236 µs | 1.23 ms | 15.7 ms |  |  | 11x faster |
| `3d/n100000/scattered` | 2.87 ms | 476 µs | 2.71 ms | 112 ms |  |  | 39x faster |
| `ghcn-tmax-20250101/heldout` | 67.5 µs | 71.8 µs | 60.4 µs | 1.15 ms | 133 ms | 1.21 ms | 17x faster |

### Construction

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | spade | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `n10/image50` | 2.49 µs | 2.31 µs | 2.49 µs | 464 µs | 2.04 µs | 1.88 µs | **1.3x slower** |
| `n10/image100` | 2.49 µs | 2.29 µs | 2.26 µs | 512 µs | 2.04 µs | 1.72 µs | **1.4x slower** |
| `n10/image500` | 2.5 µs | 2.45 µs | 2.46 µs | 473 µs | 2.03 µs | 1.72 µs | **1.5x slower** |
| `n100/image50` | 27.3 µs | 27.5 µs | 25.2 µs | 616 µs | 61.1 µs | 20.4 µs | **1.3x slower** |
| `n100/image100` | 27.3 µs | 25.2 µs | 27.2 µs | 738 µs | 56 µs | 18.7 µs | **1.5x slower** |
| `n100/image500` | 27.4 µs | 25.1 µs | 25.3 µs | 703 µs | 57.2 µs | 20.4 µs | **1.3x slower** |
| `n1000/image50` | 401 µs | 407 µs | 400 µs | 2.72 ms | 4.21 ms | 277 µs | **1.4x slower** |
| `n1000/image100` | 400 µs | 415 µs | 399 µs | 2.4 ms | 4.64 ms | 277 µs | **1.4x slower** |
| `n1000/image500` | 406 µs | 445 µs | 402 µs | 2.71 ms | 4.27 ms | 278 µs | **1.5x slower** |
| `n10000/image50` | 5.62 ms | 5.68 ms | 5.36 ms | 25.4 ms | 453 ms | 3.69 ms | **1.5x slower** |
| `n10000/image100` | 5.79 ms | 5.73 ms | 5.39 ms | 24.9 ms | 460 ms | 3.94 ms | **1.5x slower** |
| `n10000/image500` | 5.66 ms | 5.49 ms | 5.66 ms | 24 ms | 461 ms | 4 ms | **1.4x slower** |
| `n100000/image50` | 57.8 ms | 59.2 ms | 55.8 ms | 366 ms | 65.7 s | 46.9 ms | **1.2x slower** |
| `n100000/image100` | 58.5 ms | 58.2 ms | 57.2 ms | 366 ms | 65.5 s | 46.5 ms | **1.3x slower** |
| `n100000/image500` | 57.6 ms | 62.5 ms | 56.3 ms | 365 ms | 79.5 s | 44.1 ms | **1.3x slower** |
| `3d/n1000/scattered` | 4.82 ms | 4.57 ms | 4.87 ms | 10 ms |  |  | 2.1x faster |
| `3d/n10000/scattered` | 53.3 ms | 48.9 ms | 51.7 ms | 131 ms |  |  | 2.4x faster |
| `3d/n100000/scattered` | 568 ms | 554 ms | 568 ms | 2.32 s |  |  | 4.1x faster |
| `ghcn-tmax-20250101/heldout` | 5.96 ms | 5.45 ms | 5.64 ms | 30.9 ms | 498 ms | 4.56 ms | **1.3x slower** |

## `CloughTocher2DInterpolator`

`n` points of the unit square and a smooth function of them, evaluated at the pixels of a square image (some outside the hull, with no value). SciPy at its defaults (`tol=1e-6`, `maxiter=400`); interpx solves for the gradients until the residual is 1e-14 of the right-hand side.

### Evaluation

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `n10/image50` | 50.7 µs | 50.8 µs | 50.7 µs | 105 µs | 70.2 µs | 1.4x faster |
| `n50/image50` | 95.6 µs | 96.2 µs | 87.8 µs | 214 µs | 424 µs | 2.2x faster |
| `n100/image50` | 98.9 µs | 92 µs | 92.1 µs | 230 µs | 984 µs | 2.3x faster |
| `n1000/image500` | 11.5 ms | 1.94 ms | 12 ms | 24 ms | 754 ms | 2.1x faster |
| `n10000/image500` | 12.5 ms | 2.38 ms | 13.9 ms | 28.2 ms | 21.7 s | 2.3x faster |
| `n100000/image500` | 22.2 ms | 3.73 ms | 20.7 ms | 47.4 ms | 316 s | 2.1x faster |
| `ghcn-tmax-20250101/heldout` | 128 µs | 123 µs | 119 µs | 1.18 ms | 139 ms | 9.2x faster |

### Construction

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `n10/image50` | 6.05 µs | 6.02 µs | 5.99 µs | 523 µs | 3.34 µs | **1.8x slower** |
| `n50/image50` | 32.3 µs | 32.5 µs | 29.7 µs | 635 µs | 28.7 µs | **1.1x slower** |
| `n100/image50` | 65.4 µs | 60.8 µs | 60.2 µs | 779 µs | 104 µs | 1.6x faster |
| `n1000/image500` | 787 µs | 826 µs | 782 µs | 3.46 ms | 5.62 ms | 4.4x faster |
| `n10000/image500` | 10.6 ms | 10.8 ms | 10.7 ms | 34.2 ms | 595 ms | 3.2x faster |
| `n100000/image500` | 136 ms | 132 ms | 131 ms | 515 ms | 79.6 s | 3.8x faster |
| `ghcn-tmax-20250101/heldout` | 11.3 ms | 10.8 ms | 11.5 ms | 43.8 ms | 604 ms | 3.9x faster |

## `RBFInterpolator`

`n` points of a square, evaluated at 1000 points in no order (`scattered`) or at the 10 000 pixels of an image (`image`). `global`: all points; `k50`: the 50 nearest.

### Evaluation

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/global/n10/scattered` | 14.4 µs | 14.5 µs | 14.3 µs | 61.8 µs | 197 µs | 4.3x faster |
| `linear/k50/n10/scattered` | 14.4 µs | 14.5 µs | 14.4 µs | 1.39 ms |  | 97x faster |
| `thin_plate_spline/global/n10/scattered` | 72.9 µs | 72.8 µs | 79.7 µs | 111 µs | 297 µs | 1.5x faster |
| `thin_plate_spline/k50/n10/scattered` | 73 µs | 79.7 µs | 79.8 µs | 1.48 ms |  | 20x faster |
| `cubic/global/n10/scattered` | 17.5 µs | 17.1 µs | 17.2 µs | 61.7 µs | 198 µs | 3.5x faster |
| `cubic/k50/n10/scattered` | 17.2 µs | 17.3 µs | 17.2 µs | 1.38 ms |  | 80x faster |
| `quintic/global/n10/scattered` | 27.2 µs | 24.2 µs | 25.5 µs | 68.5 µs | 214 µs | 2.5x faster |
| `quintic/k50/n10/scattered` | 28.1 µs | 25.5 µs | 26.5 µs | 1.4 ms |  | 50x faster |
| `multiquadric/global/n10/scattered` | 27.8 µs | 30.6 µs | 30.5 µs | 51.9 µs | 256 µs | 1.9x faster |
| `multiquadric/k50/n10/scattered` | 30.5 µs | 28 µs | 29 µs | 1.36 ms |  | 44x faster |
| `inverse_multiquadric/global/n10/scattered` | 41.1 µs | 41.2 µs | 37.4 µs | 54.9 µs | 290 µs | 1.3x faster |
| `inverse_multiquadric/k50/n10/scattered` | 38.4 µs | 37.6 µs | 41.1 µs | 1.39 ms |  | 36x faster |
| `inverse_quadratic/global/n10/scattered` | 26.1 µs | 26.3 µs | 23.8 µs | 46 µs |  | 1.8x faster |
| `inverse_quadratic/k50/n10/scattered` | 23.8 µs | 24.1 µs | 23.8 µs | 1.35 ms |  | 57x faster |
| `gaussian/global/n10/scattered` | 67.5 µs | 67.6 µs | 67.5 µs | 86.8 µs | 353 µs | 1.3x faster |
| `gaussian/k50/n10/scattered` | 67.5 µs | 67.6 µs | 68.2 µs | 1.38 ms |  | 20x faster |
| `linear/global/n100/scattered` | 135 µs | 54.3 µs | 135 µs | 289 µs | 1.73 ms | 2.1x faster |
| `linear/k50/n100/scattered` | 13 ms | 2.52 ms | 12.9 ms | 18.3 ms |  | 1.4x faster |
| `thin_plate_spline/global/n100/scattered` | 701 µs | 190 µs | 709 µs | 959 µs | 2.95 ms | 1.4x faster |
| `thin_plate_spline/k50/n100/scattered` | 24.9 ms | 4.32 ms | 24.8 ms | 24 ms |  | **1x slower** |
| `cubic/global/n100/scattered` | 139 µs | 52.5 µs | 139 µs | 289 µs | 1.93 ms | 2.1x faster |
| `cubic/k50/n100/scattered` | 14.6 ms | 2.83 ms | 14.6 ms | 20.2 ms |  | 1.4x faster |
| `quintic/global/n100/scattered` | 167 µs | 60.7 µs | 167 µs | 300 µs | 2.26 ms | 1.8x faster |
| `quintic/k50/n100/scattered` | 16 ms | 3.06 ms | 16.4 ms | 21.7 ms |  | 1.4x faster |
| `multiquadric/global/n100/scattered` | 269 µs | 94.8 µs | 295 µs | 302 µs | 2.27 ms | 1.1x faster |
| `multiquadric/k50/n100/scattered` | 15.3 ms | 2.92 ms | 14.9 ms | 18.5 ms |  | 1.2x faster |
| `inverse_multiquadric/global/n100/scattered` | 364 µs | 120 µs | 369 µs | 392 µs | 2.67 ms | 1.1x faster |
| `inverse_multiquadric/k50/n100/scattered` | 15.6 ms | 3.02 ms | 16.2 ms | 19.4 ms |  | 1.2x faster |
| `inverse_quadratic/global/n100/scattered` | 224 µs | 77.6 µs | 224 µs | 287 µs |  | 1.3x faster |
| `inverse_quadratic/k50/n100/scattered` | 13.9 ms | 2.82 ms | 15.2 ms | 18.3 ms |  | 1.3x faster |
| `gaussian/global/n100/scattered` | 649 µs | 161 µs | 682 µs | 676 µs | 3.8 ms | 1x faster |
| `gaussian/k50/n100/scattered` | 21.4 ms | 3.74 ms | 22 ms | 22.3 ms |  | 1x faster |
| `linear/global/n1000/scattered` | 1.47 ms | 366 µs | 1.47 ms | 3.68 ms | 18.5 ms | 2.5x faster |
| `linear/k50/n1000/scattered` | 15.3 ms | 2.97 ms | 14.8 ms | 35.5 ms |  | 2.3x faster |
| `thin_plate_spline/global/n1000/scattered` | 7.54 ms | 1.41 ms | 7.93 ms | 10.9 ms | 31.9 ms | 1.4x faster |
| `thin_plate_spline/k50/n1000/scattered` | 22.6 ms | 4.14 ms | 24.4 ms | 45.7 ms |  | 2x faster |
| `cubic/global/n1000/scattered` | 1.4 ms | 339 µs | 1.49 ms | 3.98 ms | 21 ms | 2.9x faster |
| `cubic/k50/n1000/scattered` | 16.5 ms | 3.28 ms | 17 ms | 39.3 ms |  | 2.4x faster |
| `quintic/global/n1000/scattered` | 1.65 ms | 344 µs | 1.65 ms | 4.05 ms | 21.2 ms | 2.5x faster |
| `quintic/k50/n1000/scattered` | 17.8 ms | 3.46 ms | 18.5 ms | 43.8 ms |  | 2.5x faster |
| `multiquadric/global/n1000/scattered` | 2.78 ms | 627 µs | 2.96 ms | 4.21 ms | 25 ms | 1.5x faster |
| `multiquadric/k50/n1000/scattered` | 17.3 ms | 4.14 ms | 19.6 ms | 38.2 ms |  | 2.2x faster |
| `inverse_multiquadric/global/n1000/scattered` | 4.16 ms | 826 µs | 4 ms | 5.18 ms | 29.9 ms | 1.2x faster |
| `inverse_multiquadric/k50/n1000/scattered` | 18.7 ms | 3.42 ms | 17.8 ms | 39.4 ms |  | 2.1x faster |
| `inverse_quadratic/global/n1000/scattered` | 2.45 ms | 529 µs | 2.29 ms | 3.9 ms |  | 1.6x faster |
| `inverse_quadratic/k50/n1000/scattered` | 16.2 ms | 3.34 ms | 16.6 ms | 36.1 ms |  | 2.2x faster |
| `gaussian/global/n1000/scattered` | 8.91 ms | 1.25 ms | 6.83 ms | 9.02 ms | 40.4 ms | 1x faster |
| `gaussian/k50/n1000/scattered` | 21.8 ms | 4.2 ms | 22.3 ms | 44.3 ms |  | 2x faster |
| `thin_plate_spline/global/n3000/scattered` | 22.3 ms | 4.21 ms | 25 ms | 32.5 ms | 103 ms | 1.5x faster |
| `thin_plate_spline/k50/n10000/scattered` | 23 ms | 4.36 ms | 22.3 ms | 48.1 ms |  | 2.1x faster |
| `thin_plate_spline/k50/n10000/image` | 245 ms | 40 ms | 254 ms | 492 ms |  | 2x faster |
| `thin_plate_spline/k50/n100000/scattered` | 25.5 ms | 4.49 ms | 26 ms | 49.1 ms |  | 1.9x faster |
| `thin_plate_spline/k50/n100000/image` | 237 ms | 39.5 ms | 234 ms | 537 ms |  | 2.3x faster |
| `linear/global/ghcn-tmax-20250101/heldout` | 25.2 ms | 3.64 ms | 18 ms | 54.5 ms | 241 ms | 2.2x faster |
| `linear/k50/ghcn-tmax-20250101/heldout` | 18.2 ms | 3.23 ms | 18.9 ms | 45.2 ms |  | 2.5x faster |
| `thin_plate_spline/global/ghcn-tmax-20250101/heldout` | 79.6 ms | 16 ms | 79.7 ms | 109 ms | 390 ms | 1.4x faster |
| `thin_plate_spline/k50/ghcn-tmax-20250101/heldout` | 30.8 ms | 5.21 ms | 31.5 ms | 55.1 ms |  | 1.8x faster |
| `cubic/global/ghcn-tmax-20250101/heldout` | 18.6 ms | 3.58 ms | 19.2 ms | 55.9 ms | 264 ms | 3x faster |
| `cubic/k50/ghcn-tmax-20250101/heldout` | 20.6 ms | 3.68 ms | 20.4 ms | 46.4 ms |  | 2.3x faster |
| `quintic/global/ghcn-tmax-20250101/heldout` | 19.9 ms | 4.4 ms | 36.7 ms | 49.7 ms | 298 ms | 2.5x faster |
| `quintic/k50/ghcn-tmax-20250101/heldout` | 23.4 ms | 3.98 ms | 23.3 ms | 53.6 ms |  | 2.3x faster |

### Construction

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/global/n10/scattered` | 2.55 µs | 2.54 µs | 2.54 µs | 10.9 µs | 5.82 µs | 2.3x faster |
| `linear/k50/n10/scattered` | 2.56 µs | 2.59 µs | 2.56 µs | 16.2 µs |  | 6.3x faster |
| `thin_plate_spline/global/n10/scattered` | 3.27 µs | 3.22 µs | 3.27 µs | 12.2 µs | 6.45 µs | 2x faster |
| `thin_plate_spline/k50/n10/scattered` | 3.51 µs | 3.51 µs | 3.55 µs | 17.3 µs |  | 4.9x faster |
| `cubic/global/n10/scattered` | 3.4 µs | 3.39 µs | 3.7 µs | 12.1 µs | 5.47 µs | 1.6x faster |
| `cubic/k50/n10/scattered` | 3.43 µs | 3.73 µs | 3.43 µs | 15.7 µs |  | 4.6x faster |
| `quintic/global/n10/scattered` | 3.85 µs | 3.82 µs | 3.82 µs | 13.8 µs | 5.69 µs | 1.5x faster |
| `quintic/k50/n10/scattered` | 4.15 µs | 3.94 µs | 4.13 µs | 16.9 µs |  | 4.1x faster |
| `multiquadric/global/n10/scattered` | 2.63 µs | 2.86 µs | 2.87 µs | 10.7 µs | 6.57 µs | 2.5x faster |
| `multiquadric/k50/n10/scattered` | 2.88 µs | 2.86 µs | 2.6 µs | 16.3 µs |  | 5.7x faster |
| `inverse_multiquadric/global/n10/scattered` | 2.86 µs | 2.64 µs | 2.88 µs | 11.8 µs | 7.43 µs | 2.6x faster |
| `inverse_multiquadric/k50/n10/scattered` | 2.67 µs | 2.63 µs | 2.61 µs | 16.5 µs |  | 6.2x faster |
| `inverse_quadratic/global/n10/scattered` | 2.55 µs | 2.55 µs | 2.65 µs | 10.7 µs |  | 4.2x faster |
| `inverse_quadratic/k50/n10/scattered` | 2.53 µs | 2.53 µs | 2.56 µs | 15 µs |  | 5.9x faster |
| `gaussian/global/n10/scattered` | 2.73 µs | 2.73 µs | 2.72 µs | 11.9 µs | 7.65 µs | 2.8x faster |
| `gaussian/k50/n10/scattered` | 2.73 µs | 2.73 µs | 2.67 µs | 14.9 µs |  | 5.5x faster |
| `linear/global/n100/scattered` | 67.7 µs | 68.1 µs | 71.4 µs | 76.4 µs | 684 µs | 1.1x faster |
| `linear/k50/n100/scattered` | 3.1 µs | 3.35 µs | 3.06 µs | 20.7 µs |  | 6.7x faster |
| `thin_plate_spline/global/n100/scattered` | 96.7 µs | 100 µs | 99.5 µs | 122 µs | 831 µs | 1.3x faster |
| `thin_plate_spline/k50/n100/scattered` | 3.16 µs | 3.12 µs | 3.13 µs | 22 µs |  | 6.9x faster |
| `cubic/global/n100/scattered` | 70.9 µs | 71.2 µs | 72.8 µs | 80.3 µs | 709 µs | 1.1x faster |
| `cubic/k50/n100/scattered` | 3.43 µs | 3.12 µs | 3.1 µs | 21.7 µs |  | 6.3x faster |
| `quintic/global/n100/scattered` | 78.3 µs | 78.7 µs | 78.9 µs | 92.2 µs | 748 µs | 1.2x faster |
| `quintic/k50/n100/scattered` | 3.14 µs | 3.17 µs | 3.19 µs | 22.7 µs |  | 7.2x faster |
| `multiquadric/global/n100/scattered` | 74.6 µs | 83.8 µs | 84.7 µs | 76 µs | 607 µs | 1x faster |
| `multiquadric/k50/n100/scattered` | 3.07 µs | 3.39 µs | 3.1 µs | 20.8 µs |  | 6.8x faster |
| `inverse_multiquadric/global/n100/scattered` | 81.8 µs | 79 µs | 90.2 µs | 80.2 µs | 782 µs | **1x slower** |
| `inverse_multiquadric/k50/n100/scattered` | 3.06 µs | 3.07 µs | 3.15 µs | 21.9 µs |  | 7.2x faster |
| `inverse_quadratic/global/n100/scattered` | 71.5 µs | 77.8 µs | 74.6 µs | 79.3 µs |  | 1.1x faster |
| `inverse_quadratic/k50/n100/scattered` | 3.13 µs | 3.38 µs | 3.33 µs | 21.2 µs |  | 6.8x faster |
| `gaussian/global/n100/scattered` | 94.8 µs | 106 µs | 107 µs | 108 µs | 755 µs | 1.1x faster |
| `gaussian/k50/n100/scattered` | 3.37 µs | 3.06 µs | 3.2 µs | 22.7 µs |  | 6.7x faster |
| `linear/global/n1000/scattered` | 20.2 ms | 20.9 ms | 22.7 ms | 60.5 ms | 189 ms | 3x faster |
| `linear/k50/n1000/scattered` | 37.5 µs | 34.2 µs | 38.2 µs | 139 µs |  | 3.7x faster |
| `thin_plate_spline/global/n1000/scattered` | 22.1 ms | 23.7 ms | 24.3 ms | 66.6 ms | 215 ms | 3x faster |
| `thin_plate_spline/k50/n1000/scattered` | 34.1 µs | 35.3 µs | 37 µs | 138 µs |  | 4x faster |
| `cubic/global/n1000/scattered` | 21.1 ms | 21.1 ms | 21.2 ms | 58.4 ms | 203 ms | 2.8x faster |
| `cubic/k50/n1000/scattered` | 34.8 µs | 39.5 µs | 37.2 µs | 145 µs |  | 4.2x faster |
| `quintic/global/n1000/scattered` | 22.1 ms | 20.6 ms | 20.7 ms | 56.8 ms | 195 ms | 2.6x faster |
| `quintic/k50/n1000/scattered` | 34.6 µs | 38 µs | 37.1 µs | 143 µs |  | 4.1x faster |
| `multiquadric/global/n1000/scattered` | 21.2 ms | 21.2 ms | 22.4 ms | 60.3 ms | 197 ms | 2.8x faster |
| `multiquadric/k50/n1000/scattered` | 34.2 µs | 38.1 µs | 37 µs | 141 µs |  | 4.1x faster |
| `inverse_multiquadric/global/n1000/scattered` | 23.6 ms | 22.5 ms | 21.3 ms | 61.7 ms | 202 ms | 2.6x faster |
| `inverse_multiquadric/k50/n1000/scattered` | 37.6 µs | 37.6 µs | 34 µs | 141 µs |  | 3.7x faster |
| `inverse_quadratic/global/n1000/scattered` | 21.9 ms | 20.9 ms | 20.1 ms | 60.2 ms |  | 2.7x faster |
| `inverse_quadratic/k50/n1000/scattered` | 34.1 µs | 37.8 µs | 37.3 µs | 137 µs |  | 4x faster |
| `gaussian/global/n1000/scattered` | 29.3 ms | 23.5 ms | 23.7 ms | 57.7 ms | 222 ms | 2x faster |
| `gaussian/k50/n1000/scattered` | 37.5 µs | 35 µs | 37.2 µs | 142 µs |  | 3.8x faster |
| `thin_plate_spline/global/n3000/scattered` | 461 ms | 498 ms | 486 ms | 2.71 s | 9.4 s | 5.9x faster |
| `thin_plate_spline/k50/n10000/scattered` | 694 µs | 703 µs | 609 µs | 1.67 ms |  | 2.4x faster |
| `thin_plate_spline/k50/n10000/image` | 730 µs | 813 µs | 671 µs | 1.67 ms |  | 2.3x faster |
| `thin_plate_spline/k50/n100000/scattered` | 9.67 ms | 9.43 ms | 8.56 ms | 26.5 ms |  | 2.7x faster |
| `thin_plate_spline/k50/n100000/image` | 8.31 ms | 9.2 ms | 8.21 ms | 27.1 ms |  | 3.3x faster |
| `linear/global/ghcn-tmax-20250101/heldout` | 20.5 s | 16.6 s | 16.7 s | 119 s | 295 s | 5.8x faster |
| `linear/k50/ghcn-tmax-20250101/heldout` | 698 µs | 811 µs | 695 µs | 1.71 ms |  | 2.4x faster |
| `thin_plate_spline/global/ghcn-tmax-20250101/heldout` | 17 s | 17 s | 17 s | 109 s | 285 s | 6.4x faster |
| `thin_plate_spline/k50/ghcn-tmax-20250101/heldout` | 733 µs | 735 µs | 630 µs | 1.92 ms |  | 2.6x faster |
| `cubic/global/ghcn-tmax-20250101/heldout` | 18.8 s | 17.4 s | 17.9 s | 149 s | 296 s | 7.9x faster |
| `cubic/k50/ghcn-tmax-20250101/heldout` | 718 µs | 728 µs | 599 µs | 1.76 ms |  | 2.4x faster |
| `quintic/global/ghcn-tmax-20250101/heldout` | 16.7 s | 20.8 s | 21.3 s | 111 s | 353 s | 6.7x faster |
| `quintic/k50/ghcn-tmax-20250101/heldout` | 744 µs | 790 µs | 741 µs | 1.81 ms |  | 2.4x faster |

## Memory

Megabytes. "Data read" is the process with the data of the case in memory, before anything is built. The columns after it are what was added to that: the highest while building, with the interpolator built, the highest while evaluating. "Highest" is the highest of the whole run, repeats included.

| case | library | data read | while building | built | while evaluating | highest |
| :--- | :--- | ---: | ---: | ---: | ---: | ---: |
| `oned/cubic_spline/n10000/sorted` | interpx | 4.1 | +1.2 | +1.2 | +1.2 | 6.0 |
| `oned/cubic_spline/n10000/sorted` | interpx, parallel | 4.0 | +1.3 | +1.3 | +1.3 | 6.4 |
| `oned/cubic_spline/n10000/sorted` | SciPy | 71.1 | +1.5 | +1.9 | +1.9 | 73.3 |
| `oned/cubic_spline/n10000/sorted` | scirs2-interpolate | 3.6 | +1.5 | +1.6 | +2.2 | 5.9 |
| `oned/cubic_spline/n10000/sorted` | ndarray-interp | 3.2 | +1.8 | +1.8 | +1.8 | 5.5 |
| `grid/linear/4d/max200/eval10000/ascending` | interpx | 223.9 | +181.2 | +181.2 | +181.2 | 405.1 |
| `grid/linear/4d/max200/eval10000/ascending` | interpx, parallel | 227.0 | +177.7 | +177.7 | +178.0 | 405.4 |
| `grid/linear/4d/max200/eval10000/ascending` | interpx, `from_owned` | 235.9 | +168.8 | +168.8 | +168.8 | 404.9 |
| `grid/linear/4d/max200/eval10000/ascending` | SciPy | 271.2 | +0.3 | +0.3 | +2.1 | 273.4 |
| `grid/linear/4d/max200/eval10000/ascending` | scirs2-interpolate | 204.9 | +198.7 | +200.1 | +200.1 | 404.9 |
| `grid/linear/4d/max200/eval10000/ascending` | interpn | 205.3 | +0.0 | +0.0 | +0.0 | 400.4 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpx | 5.6 | +2.3 | +2.7 | +2.7 | 8.5 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpx, parallel | 5.6 | +2.3 | +2.7 | +2.7 | 8.9 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpx, `from_owned` | 4.2 | +3.9 | +4.1 | +4.1 | 8.5 |
| `grid/cubic/2d/each400/eval1000/ascending` | SciPy | 72.6 | +86.1 | +5.4 | +5.4 | 161.9 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpn | 5.4 | +0.1 | +0.1 | +0.1 | 5.5 |
| `grid/linear/2d/each800/eval10000/ascending` | interpx | 12.9 | +0.1 | +1.9 | +1.9 | 15.0 |
| `grid/linear/2d/each800/eval10000/ascending` | interpx, parallel | 10.0 | +4.8 | +4.9 | +4.9 | 15.5 |
| `grid/linear/2d/each800/eval10000/ascending` | interpx, `from_owned` | 11.1 | +0.0 | +3.7 | +3.7 | 15.0 |
| `grid/linear/2d/each800/eval10000/ascending` | SciPy | 76.2 | +0.2 | +0.2 | +0.5 | 77.1 |
| `grid/linear/2d/each800/eval10000/ascending` | scirs2-interpolate | 12.8 | +0.0 | +2.0 | +2.1 | 15.0 |
| `grid/linear/2d/each800/eval10000/ascending` | ndarray-interp | 12.8 | +0.0 | +0.0 | +0.0 | 12.8 |
| `grid/linear/2d/each800/eval10000/ascending` | interpn | 11.6 | +0.0 | +0.0 | +0.0 | 11.6 |
| `nearest_nd/n100000/image500` | interpx | 11.6 | +3.9 | +2.9 | +4.9 | 18.5 |
| `nearest_nd/n100000/image500` | interpx, parallel | 13.0 | +2.5 | +1.2 | +3.7 | 18.4 |
| `nearest_nd/n100000/image500` | interpx, `from_owned` | 12.3 | +3.8 | +2.3 | +4.3 | 18.6 |
| `nearest_nd/n100000/image500` | SciPy | 77.3 | +4.4 | +4.6 | +12.1 | 92.0 |
| `linear_nd/n100000/image500` | interpx | 13.6 | +21.1 | +21.9 | +23.7 | 39.5 |
| `linear_nd/n100000/image500` | interpx, parallel | 14.3 | +20.4 | +21.2 | +23.0 | 38.1 |
| `linear_nd/n100000/image500` | interpx, `from_owned` | 14.8 | +20.2 | +20.2 | +22.2 | 39.0 |
| `linear_nd/n100000/image500` | SciPy | 77.4 | +67.9 | +13.8 | +25.6 | 145.7 |
| `linear_nd/n100000/image500` | scirs2-interpolate | 5.3 | +19.9 | +16.7 | +18.8 | 26.1 |
| `linear_nd/n100000/image500` | spade | 14.0 | +14.1 | +11.1 | +13.1 | 29.5 |
| `clough_tocher/n100000/image500` | interpx | 14.4 | +48.3 | +28.2 | +30.2 | 64.5 |
| `clough_tocher/n100000/image500` | interpx, parallel | 13.2 | +49.5 | +28.3 | +30.8 | 64.9 |
| `clough_tocher/n100000/image500` | interpx, `from_owned` | 13.3 | +49.1 | +28.1 | +30.1 | 64.5 |
| `clough_tocher/n100000/image500` | SciPy | 75.9 | +69.2 | +19.3 | +31.5 | 146.7 |
| `clough_tocher/n100000/image500` | scirs2-interpolate | 14.2 | +12.5 | +9.3 | +11.5 | 27.7 |
| `rbf/thin_plate_spline/global/n3000/scattered` | interpx | 4.3 | +76.2 | +3.9 | +3.9 | 80.5 |
| `rbf/thin_plate_spline/global/n3000/scattered` | interpx, parallel | 4.4 | +76.1 | +3.8 | +4.2 | 80.9 |
| `rbf/thin_plate_spline/global/n3000/scattered` | interpx, `from_owned` | 3.3 | +77.2 | +4.8 | +4.8 | 80.6 |
| `rbf/thin_plate_spline/global/n3000/scattered` | SciPy | 71.2 | +73.0 | +0.9 | +9.0 | 144.2 |
| `rbf/thin_plate_spline/global/n3000/scattered` | scirs2-interpolate | 3.1 | +289.4 | +1.4 | +1.4 | 292.5 |
| `rbf/thin_plate_spline/k50/n100000/image` | interpx | 5.7 | +4.5 | +5.7 | +6.0 | 11.9 |
| `rbf/thin_plate_spline/k50/n100000/image` | interpx, parallel | 6.2 | +4.1 | +5.3 | +6.4 | 12.7 |
| `rbf/thin_plate_spline/k50/n100000/image` | interpx, `from_owned` | 6.9 | +3.3 | +4.5 | +4.9 | 13.5 |
| `rbf/thin_plate_spline/k50/n100000/image` | SciPy | 73.5 | +5.4 | +5.7 | +22.2 | 101.0 |

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
