# Performance: interpx against SciPy and other Rust crates

<!-- Written by compare/perf/report.py from a run of compare/perf/run.py. Do not edit by hand: run it again. -->

Measured on 2026-10-18, at commit `58c8390a56`.

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
| Regular grid | evaluation | 92 | 72 | 0 | 20 |
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
| `grid/cubic` | construction | interpn | no | 8 | 1.73e+03x to 9.54e+07x | `grid/cubic/2d/auckland-ba32-0401/heldout/ascending`: 83.8 ms against 0.88 ns |
| `grid/linear` | construction | interpn | yes | 35 | 172x to 4.44e+07x | `grid/linear/4d/max200/eval10/ascending`: 56.4 ms against 1.27 ns |
| `grid/nearest` | construction | interpn | yes | 3 | 822x to 2.14e+07x | `grid/nearest/2d/auckland-ba32-0401/heldout/ascending`: 16.7 ms against 0.78 ns |
| `grid/linear` | construction | SciPy | yes | 35 | 1.17x to 5.67e+03x | `grid/linear/4d/max200/eval10/flipped`: 115 ms against 20.2 µs |
| `grid/nearest` | construction | SciPy | yes | 2 | 2.59x to 2.17e+03x | `grid/nearest/2d/auckland-ba32-0401/heldout/ascending`: 16.7 ms against 7.68 µs |
| `grid/linear` | construction | ndarray-interp | yes | 15 | 1.44x to 2e+03x | `grid/linear/2d/auckland-ba32-0401/heldout/ascending`: 16.5 ms against 8.24 µs |
| `grid/pchip` | construction | SciPy | yes | 1 | 50.5x | `grid/pchip/3d/each40/eval10000/ascending`: 474 µs against 9.37 µs |
| `linear_nd` | construction | spade | yes | 16 | 1.62x to 5.01x | `linear_nd/n10/image500`: 9.5 µs against 1.9 µs |
| `linear_nd` | construction | scirs2-interpolate | yes | 3 | 2.63x to 4.65x | `linear_nd/n10/image500`: 9.5 µs against 2.04 µs |
| `oned/pchip` | construction | scirs2-interpolate | yes | 12 | 1.33x to 2.83x | `oned/pchip/n10000/sorted-uneven`: 126 µs against 44.4 µs |
| `oned/akima` | construction | scirs2-interpolate | yes | 8 | 1.44x to 2.65x | `oned/akima/n10000/sorted`: 120 µs against 45.4 µs |
| `grid/linear` | construction | scirs2-interpolate | yes | 28 | 1.12x to 2.45x | `grid/linear/2d/max200/eval10/ascending`: 7.36 µs against 3 µs |
| `grid/nearest` | construction | scirs2-interpolate | yes | 3 | 1.24x to 2.23x | `grid/nearest/3d/each40/eval10000/ascending`: 21.1 µs against 9.47 µs |
| `oned/linear` | construction | ndarray-interp | yes | 12 | 1.28x to 2.09x | `oned/linear/n10/sorted`: 86.5 ns against 41.4 ns |
| `oned/linear` | construction | scirs2-interpolate | yes | 12 | 1.21x to 1.78x | `oned/linear/n10000/shuffled-uneven`: 22.9 µs against 12.9 µs |
| `grid/linear` | evaluation | interpn | yes | 18 | 1.1x to 1.48x | `grid/linear/2d/max10/eval1000/ascending`: 13.7 µs against 9.24 µs |
| `oned/cubic_spline` | construction | scirs2-interpolate | no | 2 | 1.18x to 1.47x | `oned/cubic_spline/n10000/sorted-uneven`: 232 µs against 158 µs |
| `grid/linear` | evaluation | ndarray-interp | yes | 2 | 1.12x to 1.12x | `grid/linear/2d/max200/eval10/ascending`: 211 ns against 188 ns |
| `grid/nearest` | evaluation | interpn | yes | 1 | 1.12x | `grid/nearest/2d/each40/eval10000/ascending`: 165 µs against 148 µs |

<details><summary>Every comparison where interpx is slower</summary>

"Values apart" is the distance of that library's values from SciPy's, as a share of the largest value.

| case | operation | against | interpx | the other | slower by | values apart |
| :--- | :--- | :--- | ---: | ---: | ---: | ---: |
| `grid/cubic/2d/auckland-ba32-0401/heldout/ascending` | construction | interpn | 83.8 ms | 0.88 ns | 95370467.56x | 5e-03 |
| `grid/linear/4d/max200/eval10/ascending` | construction | interpn | 56.4 ms | 1.27 ns | 44414188.31x | 2e-16 |
| `grid/linear/4d/max200/eval100/ascending` | construction | interpn | 51.8 ms | 1.27 ns | 40812386.46x | 3e-16 |
| `grid/linear/4d/max200/eval1000/ascending` | construction | interpn | 49.9 ms | 1.27 ns | 39281506.46x | 3e-16 |
| `grid/linear/4d/max200/eval10000/ascending` | construction | interpn | 51.6 ms | 1.37 ns | 37712091.43x | 4e-16 |
| `grid/linear/2d/auckland-ba32-0401/heldout/ascending` | construction | interpn | 16.5 ms | 0.68 ns | 24155574.86x | 0e+00 |
| `grid/nearest/2d/auckland-ba32-0401/heldout/ascending` | construction | interpn | 16.7 ms | 0.78 ns | 21387264.00x | 0e+00 |
| `grid/linear/3d/max200/eval100/ascending` | construction | interpn | 1.86 ms | 0.88 ns | 2120817.78x | 3e-16 |
| `grid/linear/3d/max200/eval10000/ascending` | construction | interpn | 1.78 ms | 0.88 ns | 2023537.78x | 5e-16 |
| `grid/linear/3d/max200/eval1000/ascending` | construction | interpn | 1.77 ms | 0.88 ns | 2009656.89x | 3e-16 |
| `grid/linear/3d/max200/eval10/ascending` | construction | interpn | 1.82 ms | 0.98 ns | 1860505.60x | 3e-16 |
| `grid/cubic/2d/each400/eval1000/ascending` | construction | interpn | 1.48 ms | 0.88 ns | 1681521.78x | 5e-01 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | interpn | 1.04 ms | 0.78 ns | 1330944.00x | 3e-16 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | interpn | 320 µs | 0.78 ns | 410112.00x | 3e-16 |
| `grid/cubic/3d/each40/eval10000/ascending` | construction | interpn | 321 µs | 1.17 ns | 273792.00x | 5e-01 |
| `grid/cubic/2d/each200/eval1000/ascending` | construction | interpn | 190 µs | 0.98 ns | 194969.60x | 2e-01 |
| `grid/cubic/2d/each100/eval1000/ascending` | construction | interpn | 43.1 µs | 0.98 ns | 44160.00x | 2e-01 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | interpn | 21.1 µs | 0.88 ns | 23992.89x | 0e+00 |
| `grid/linear/4d/max40/eval10000/ascending` | construction | interpn | 14.3 µs | 1.27 ns | 11249.23x | 5e-16 |
| `grid/linear/4d/max40/eval1000/ascending` | construction | interpn | 14.2 µs | 1.27 ns | 11190.15x | 6e-16 |
| `grid/linear/2d/max200/eval10/ascending` | construction | interpn | 7.36 µs | 0.68 ns | 10770.29x | 2e-16 |
| `grid/linear/4d/max40/eval10/ascending` | construction | interpn | 13.5 µs | 1.27 ns | 10601.85x | 2e-16 |
| `grid/cubic/2d/each40/eval1000/ascending` | construction | interpn | 9.19 µs | 0.88 ns | 10456.89x | 7e-01 |
| `grid/linear/2d/max200/eval10000/ascending` | construction | interpn | 7.49 µs | 0.78 ns | 9592.00x | 3e-16 |
| `grid/linear/4d/max40/eval100/ascending` | construction | interpn | 13.1 µs | 1.37 ns | 9577.14x | 2e-16 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | interpn | 7.41 µs | 0.78 ns | 9488.00x | 2e-16 |
| `grid/cubic/2d/each40/eval10000/ascending` | construction | interpn | 9.22 µs | 0.98 ns | 9436.80x | 7e-01 |
| `grid/linear/2d/max200/eval100/ascending` | construction | interpn | 6.83 µs | 0.78 ns | 8740.00x | 1e-16 |
| `grid/linear/4d/max200/eval10/flipped` | construction | SciPy | 115 ms | 20.2 µs | 5670.47x |  |
| `grid/linear/4d/max200/eval10/ascending` | construction | SciPy | 56.4 ms | 10.2 µs | 5518.66x |  |
| `grid/linear/4d/max200/eval1000/flipped` | construction | SciPy | 118 ms | 22.4 µs | 5250.84x |  |
| `grid/linear/4d/max200/eval100/flipped` | construction | SciPy | 115 ms | 21.9 µs | 5221.58x |  |
| `grid/linear/4d/max200/eval100/ascending` | construction | SciPy | 51.8 ms | 10.1 µs | 5128.37x |  |
| `grid/linear/4d/max200/eval10000/ascending` | construction | SciPy | 51.6 ms | 10.4 µs | 4973.34x |  |
| `grid/linear/4d/max200/eval1000/ascending` | construction | SciPy | 49.9 ms | 10.2 µs | 4871.22x |  |
| `grid/linear/4d/max200/eval10000/flipped` | construction | SciPy | 105 ms | 22 µs | 4751.13x |  |
| `grid/linear/3d/max40/eval10/ascending` | construction | interpn | 2.96 µs | 0.88 ns | 3364.44x | 2e-16 |
| `grid/linear/3d/max40/eval100/ascending` | construction | interpn | 3.22 µs | 0.98 ns | 3293.60x | 2e-16 |
| `grid/linear/3d/max40/eval1000/ascending` | construction | interpn | 3.12 µs | 0.98 ns | 3199.60x | 2e-16 |
| `grid/linear/3d/max40/eval10000/ascending` | construction | interpn | 2.97 µs | 0.98 ns | 3045.60x | 3e-16 |
| `grid/linear/2d/auckland-ba32-0401/heldout/ascending` | construction | SciPy | 16.5 ms | 7.59 µs | 2175.39x |  |
| `grid/nearest/2d/auckland-ba32-0401/heldout/ascending` | construction | SciPy | 16.7 ms | 7.68 µs | 2174.39x |  |
| `grid/linear/2d/auckland-ba32-0401/heldout/ascending` | construction | ndarray-interp | 16.5 ms | 8.24 µs | 2003.80x | 0e+00 |
| `grid/cubic/2d/each10/eval1000/ascending` | construction | interpn | 1.52 µs | 0.88 ns | 1733.78x | 3e-01 |
| `grid/nearest/2d/each40/eval10000/ascending` | construction | interpn | 642 ns | 0.78 ns | 821.50x | 0e+00 |
| `grid/linear/2d/max40/eval1000/ascending` | construction | interpn | 397 ns | 0.78 ns | 508.25x | 3e-16 |
| `grid/linear/2d/max40/eval10000/ascending` | construction | interpn | 387 ns | 0.78 ns | 495.12x | 4e-16 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | ndarray-interp | 1.04 ms | 2.23 µs | 467.16x | 6e-16 |
| `grid/linear/2d/max40/eval10/ascending` | construction | interpn | 362 ns | 0.78 ns | 464.00x | 2e-16 |
| `grid/linear/2d/max40/eval100/ascending` | construction | interpn | 361 ns | 0.78 ns | 461.75x | 2e-16 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | ndarray-interp | 320 µs | 1.14 µs | 280.51x | 5e-16 |
| `grid/linear/3d/max200/eval10000/ascending` | construction | SciPy | 1.78 ms | 7.88 µs | 225.75x |  |
| `grid/linear/3d/max200/eval100/ascending` | construction | SciPy | 1.86 ms | 8.26 µs | 225.68x |  |
| `grid/linear/3d/max200/eval10/ascending` | construction | SciPy | 1.82 ms | 8.26 µs | 220.06x |  |
| `grid/linear/3d/max200/eval1000/ascending` | construction | SciPy | 1.77 ms | 8.27 µs | 213.45x |  |
| `grid/linear/3d/max10/eval10000/ascending` | construction | interpn | 181 ns | 0.88 ns | 205.39x | 4e-16 |
| `grid/linear/3d/max10/eval10/ascending` | construction | interpn | 196 ns | 0.98 ns | 200.60x | 2e-16 |
| `grid/linear/3d/max10/eval1000/ascending` | construction | interpn | 195 ns | 0.98 ns | 199.40x | 5e-16 |
| `grid/linear/2d/max10/eval10000/ascending` | construction | interpn | 134 ns | 0.68 ns | 196.43x | 3e-16 |
| `grid/linear/3d/max200/eval1000/flipped` | construction | SciPy | 3.11 ms | 15.9 µs | 195.70x |  |
| `grid/linear/3d/max10/eval100/ascending` | construction | interpn | 187 ns | 0.98 ns | 191.90x | 3e-16 |
| `grid/linear/2d/max10/eval1000/ascending` | construction | interpn | 123 ns | 0.68 ns | 179.29x | 3e-16 |
| `grid/linear/3d/max200/eval100/flipped` | construction | SciPy | 3.14 ms | 17.5 µs | 179.17x |  |
| `grid/linear/3d/max200/eval10000/flipped` | construction | SciPy | 3.13 ms | 17.6 µs | 177.98x |  |
| `grid/linear/3d/max200/eval10/flipped` | construction | SciPy | 3.05 ms | 17.4 µs | 175.24x |  |
| `grid/linear/2d/max10/eval100/ascending` | construction | interpn | 135 ns | 0.78 ns | 172.38x | 2e-16 |
| `grid/linear/2d/max10/eval10/ascending` | construction | interpn | 134 ns | 0.78 ns | 171.62x | 1e-16 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | SciPy | 1.04 ms | 6.95 µs | 149.54x |  |
| `grid/pchip/3d/each40/eval10000/ascending` | construction | SciPy | 474 µs | 9.37 µs | 50.55x |  |
| `grid/linear/2d/each400/eval10000/ascending` | construction | SciPy | 320 µs | 6.5 µs | 49.32x |  |
| `grid/linear/2d/max200/eval10000/ascending` | construction | ndarray-interp | 7.49 µs | 475 ns | 15.79x | 5e-16 |
| `grid/linear/2d/max200/eval10/ascending` | construction | ndarray-interp | 7.36 µs | 470 ns | 15.66x | 2e-16 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | ndarray-interp | 7.41 µs | 515 ns | 14.39x | 4e-16 |
| `grid/linear/2d/max200/eval100/ascending` | construction | ndarray-interp | 6.83 µs | 515 ns | 13.26x | 3e-16 |
| `linear_nd/n10/image500` | construction | spade | 9.5 µs | 1.9 µs | 5.01x | 9e-16 |
| `linear_nd/n10/image500` | construction | scirs2-interpolate | 9.5 µs | 2.04 µs | 4.65x | 9e-16 |
| `linear_nd/n10/image50` | construction | spade | 5.7 µs | 1.76 µs | 3.24x | 6e-16 |
| `linear_nd/n10/image50` | construction | scirs2-interpolate | 5.7 µs | 1.88 µs | 3.03x | 8e-16 |
| `grid/linear/4d/max40/eval1000/flipped` | construction | SciPy | 57.2 µs | 19.8 µs | 2.89x |  |
| `oned/pchip/n10000/sorted-uneven` | construction | scirs2-interpolate | 126 µs | 44.4 µs | 2.83x | 3e-16 |
| `linear_nd/n10/image100` | construction | spade | 5.34 µs | 1.91 µs | 2.80x | 8e-16 |
| `oned/pchip/n10000/sorted` | construction | scirs2-interpolate | 124 µs | 46.7 µs | 2.66x | 2e-16 |
| `oned/akima/n10000/sorted` | construction | scirs2-interpolate | 120 µs | 45.4 µs | 2.65x | 1e-16 |
| `grid/linear/4d/max40/eval100/flipped` | construction | SciPy | 53.2 µs | 20.2 µs | 2.64x |  |
| `linear_nd/n10/image100` | construction | scirs2-interpolate | 5.34 µs | 2.03 µs | 2.63x | 9e-16 |
| `grid/linear/2d/max40/eval1000/ascending` | construction | ndarray-interp | 397 ns | 153 ns | 2.60x | 4e-16 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | SciPy | 21.1 µs | 8.16 µs | 2.59x |  |
| `grid/linear/4d/max40/eval10/flipped` | construction | SciPy | 51.5 µs | 20 µs | 2.58x |  |
| `grid/linear/2d/max40/eval10000/ascending` | construction | ndarray-interp | 387 ns | 152 ns | 2.54x | 4e-16 |
| `grid/linear/2d/max200/eval10/ascending` | construction | scirs2-interpolate | 7.36 µs | 3 µs | 2.45x | 1e-16 |
| `grid/linear/4d/max40/eval10000/flipped` | construction | SciPy | 51.9 µs | 21.7 µs | 2.40x |  |
| `grid/linear/2d/max40/eval10/ascending` | construction | ndarray-interp | 362 ns | 152 ns | 2.38x | 2e-16 |
| `grid/linear/4d/max40/eval10000/ascending` | construction | scirs2-interpolate | 14.3 µs | 6.02 µs | 2.37x | 4e-16 |
| `grid/linear/4d/max40/eval1000/ascending` | construction | scirs2-interpolate | 14.2 µs | 5.99 µs | 2.37x | 2e-16 |
| `grid/linear/2d/max40/eval100/ascending` | construction | ndarray-interp | 361 ns | 152 ns | 2.37x | 2e-16 |
| `grid/linear/2d/max200/eval100/ascending` | construction | scirs2-interpolate | 6.83 µs | 2.95 µs | 2.32x | 2e-16 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | scirs2-interpolate | 7.41 µs | 3.23 µs | 2.29x | 2e-16 |
| `linear_nd/n1000/image50` | construction | spade | 627 µs | 275 µs | 2.28x | 1e-15 |
| `grid/linear/4d/max40/eval10/ascending` | construction | scirs2-interpolate | 13.5 µs | 5.97 µs | 2.25x | 2e-16 |
| `grid/linear/2d/max200/eval10000/ascending` | construction | scirs2-interpolate | 7.49 µs | 3.36 µs | 2.23x | 2e-16 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | scirs2-interpolate | 21.1 µs | 9.47 µs | 2.23x | 0e+00 |
| `grid/linear/4d/max40/eval100/ascending` | construction | scirs2-interpolate | 13.1 µs | 5.97 µs | 2.19x | 3e-16 |
| `linear_nd/n1000/image500` | construction | spade | 655 µs | 299 µs | 2.19x | 3e-15 |
| `linear_nd/n1000/image100` | construction | spade | 625 µs | 286 µs | 2.19x | 1e-15 |
| `grid/nearest/2d/each40/eval10000/ascending` | construction | scirs2-interpolate | 642 ns | 295 ns | 2.18x | 0e+00 |
| `oned/pchip/n10000/shuffled-uneven` | construction | scirs2-interpolate | 91.9 µs | 42.7 µs | 2.15x | 3e-16 |
| `oned/akima/n10000/shuffled` | construction | scirs2-interpolate | 93.8 µs | 44.9 µs | 2.09x | 1e-16 |
| `oned/linear/n10/sorted` | construction | ndarray-interp | 86.5 ns | 41.4 ns | 2.09x | 0e+00 |
| `oned/linear/n10/shuffled` | construction | ndarray-interp | 86.6 ns | 41.5 ns | 2.09x | 0e+00 |
| `grid/linear/3d/max40/eval10000/ascending` | construction | scirs2-interpolate | 2.97 µs | 1.44 µs | 2.07x | 2e-16 |
| `grid/linear/3d/max40/eval10/ascending` | construction | scirs2-interpolate | 2.96 µs | 1.45 µs | 2.05x | 1e-16 |
| `grid/linear/3d/max40/eval100/ascending` | construction | scirs2-interpolate | 3.22 µs | 1.58 µs | 2.03x | 1e-16 |
| `linear_nd/n10000/image500` | construction | spade | 7.94 ms | 3.98 ms | 2.00x | 3e-15 |
| `oned/pchip/n1000/shuffled` | construction | scirs2-interpolate | 8.26 µs | 4.18 µs | 1.98x | 2e-16 |
| `grid/linear/3d/max40/eval1000/ascending` | construction | scirs2-interpolate | 3.12 µs | 1.59 µs | 1.97x | 2e-16 |
| `linear_nd/n10000/image50` | construction | spade | 7.95 ms | 4.06 ms | 1.96x | 8e-16 |
| `linear_nd/n100/image100` | construction | spade | 36.8 µs | 18.8 µs | 1.95x | 9e-16 |
| `oned/pchip/n10000/shuffled` | construction | scirs2-interpolate | 83.5 µs | 42.8 µs | 1.95x | 2e-16 |
| `linear_nd/n10000/image100` | construction | spade | 8.05 ms | 4.19 ms | 1.92x | 8e-16 |
| `oned/akima/n1000/shuffled` | construction | scirs2-interpolate | 8.36 µs | 4.39 µs | 1.91x | 3e-17 |
| `oned/akima/n100/shuffled` | construction | scirs2-interpolate | 968 ns | 529 ns | 1.83x | 1e-16 |
| `grid/linear/2d/max40/eval1000/ascending` | construction | scirs2-interpolate | 397 ns | 218 ns | 1.82x | 3e-16 |
| `oned/pchip/n1000/sorted` | construction | scirs2-interpolate | 8.27 µs | 4.55 µs | 1.82x | 2e-16 |
| `linear_nd/n100000/image500` | construction | spade | 85.6 ms | 47.3 ms | 1.81x | 6e-15 |
| `grid/linear/2d/max40/eval10/ascending` | construction | scirs2-interpolate | 362 ns | 203 ns | 1.79x | 4e-17 |
| `linear_nd/n100/image50` | construction | spade | 36.8 µs | 20.7 µs | 1.78x | 1e-15 |
| `oned/linear/n10000/shuffled-uneven` | construction | scirs2-interpolate | 22.9 µs | 12.9 µs | 1.78x | 6e-17 |
| `linear_nd/n100/image500` | construction | spade | 36.8 µs | 20.9 µs | 1.76x | 1e-15 |
| `oned/pchip/n100/shuffled` | construction | scirs2-interpolate | 970 ns | 551 ns | 1.76x | 2e-16 |
| `grid/linear/2d/max40/eval10000/ascending` | construction | scirs2-interpolate | 387 ns | 221 ns | 1.75x | 2e-16 |
| `oned/akima/n1000/sorted` | construction | scirs2-interpolate | 8.33 µs | 4.81 µs | 1.73x | 3e-17 |
| `grid/linear/2d/max10/eval10000/ascending` | construction | ndarray-interp | 134 ns | 78 ns | 1.72x | 4e-16 |
| `linear_nd/n100000/image100` | construction | spade | 81.1 ms | 47.3 ms | 1.72x | 7e-16 |
| `grid/linear/2d/max10/eval10/ascending` | construction | ndarray-interp | 134 ns | 78.1 ns | 1.72x | 1e-16 |
| `oned/linear/n10000/shuffled` | construction | scirs2-interpolate | 21.3 µs | 12.5 µs | 1.70x | 0e+00 |
| `linear_nd/n100000/image50` | construction | spade | 79.3 ms | 47.7 ms | 1.66x | 8e-16 |
| `grid/linear/2d/max40/eval100/ascending` | construction | scirs2-interpolate | 361 ns | 218 ns | 1.66x | 3e-16 |
| `linear_nd/ghcn-tmax-20250101/heldout` | construction | spade | 7.94 ms | 4.9 ms | 1.62x | 7e-16 |
| `oned/linear/n1000/sorted` | construction | scirs2-interpolate | 1.99 µs | 1.23 µs | 1.62x | 0e+00 |
| `oned/linear/n10000/sorted` | construction | scirs2-interpolate | 19.8 µs | 12.4 µs | 1.59x | 0e+00 |
| `oned/linear/n100/sorted` | construction | ndarray-interp | 264 ns | 166 ns | 1.59x | 0e+00 |
| `grid/linear/2d/max10/eval100/ascending` | construction | ndarray-interp | 135 ns | 85.1 ns | 1.58x | 4e-16 |
| `oned/akima/n100/sorted` | construction | scirs2-interpolate | 887 ns | 561 ns | 1.58x | 1e-16 |
| `oned/linear/n1000/shuffled` | construction | scirs2-interpolate | 1.93 µs | 1.24 µs | 1.56x | 0e+00 |
| `oned/pchip/n100/sorted-uneven` | construction | scirs2-interpolate | 883 ns | 580 ns | 1.52x | 2e-16 |
| `grid/linear/4d/max40/eval10000/ascending` | construction | SciPy | 14.3 µs | 9.41 µs | 1.52x |  |
| `oned/linear/n50/shuffled` | construction | ndarray-interp | 151 ns | 99.8 ns | 1.52x | 0e+00 |
| `oned/linear/n10000/sorted-uneven` | construction | scirs2-interpolate | 19.5 µs | 12.9 µs | 1.51x | 6e-17 |
| `oned/linear/n10000/shuffled-uneven` | construction | ndarray-interp | 22.9 µs | 15.2 µs | 1.51x | 0e+00 |
| `oned/linear/n50/sorted` | construction | ndarray-interp | 150 ns | 99.9 ns | 1.50x | 0e+00 |
| `grid/linear/4d/max40/eval1000/ascending` | construction | SciPy | 14.2 µs | 9.48 µs | 1.50x |  |
| `oned/pchip/n100/sorted` | construction | scirs2-interpolate | 883 ns | 591 ns | 1.50x | 2e-16 |
| `grid/linear/2d/max10/eval1000/ascending` | evaluation | interpn | 13.7 µs | 9.24 µs | 1.48x | 3e-16 |
| `oned/akima/n50/sorted` | construction | scirs2-interpolate | 510 ns | 346 ns | 1.48x | 1e-16 |
| `grid/linear/2d/max10/eval10000/ascending` | evaluation | interpn | 136 µs | 92.5 µs | 1.47x | 3e-16 |
| `grid/linear/2d/max10/eval100/ascending` | evaluation | interpn | 1.38 µs | 935 ns | 1.47x | 2e-16 |
| `oned/cubic_spline/n10000/sorted-uneven` | construction | scirs2-interpolate | 232 µs | 158 µs | 1.47x | 7e-06 |
| `grid/linear/3d/max10/eval10000/ascending` | evaluation | interpn | 207 µs | 142 µs | 1.46x | 4e-16 |
| `oned/linear/n100/shuffled-uneven` | construction | ndarray-interp | 265 ns | 182 ns | 1.46x | 0e+00 |
| `oned/linear/n100/sorted-uneven` | construction | ndarray-interp | 242 ns | 167 ns | 1.45x | 0e+00 |
| `oned/linear/n100/shuffled` | construction | ndarray-interp | 241 ns | 166 ns | 1.45x | 0e+00 |
| `oned/akima/n50/shuffled` | construction | scirs2-interpolate | 517 ns | 359 ns | 1.44x | 1e-16 |
| `grid/linear/4d/max200/eval10/ascending` | construction | scirs2-interpolate | 56.4 ms | 39.2 ms | 1.44x | 2e-16 |
| `grid/linear/3d/max200/eval10/ascending` | construction | scirs2-interpolate | 1.82 ms | 1.26 ms | 1.44x | 3e-16 |
| `grid/linear/2d/max10/eval1000/ascending` | construction | ndarray-interp | 123 ns | 85.4 ns | 1.44x | 3e-16 |
| `grid/linear/4d/max40/eval10/ascending` | construction | SciPy | 13.5 µs | 9.44 µs | 1.43x |  |
| `oned/pchip/n100/shuffled-uneven` | construction | scirs2-interpolate | 884 ns | 626 ns | 1.41x | 2e-16 |
| `grid/linear/3d/max10/eval10/ascending` | evaluation | interpn | 216 ns | 155 ns | 1.39x | 2e-16 |
| `grid/linear/3d/max200/eval1000/ascending` | construction | scirs2-interpolate | 1.77 ms | 1.28 ms | 1.38x | 3e-16 |
| `oned/pchip/n50/sorted` | construction | scirs2-interpolate | 501 ns | 364 ns | 1.38x | 1e-16 |
| `grid/linear/3d/max10/eval1000/ascending` | evaluation | interpn | 20.8 µs | 15.2 µs | 1.37x | 5e-16 |
| `oned/linear/n100/sorted` | construction | scirs2-interpolate | 264 ns | 196 ns | 1.34x | 0e+00 |
| `oned/pchip/n50/shuffled` | construction | scirs2-interpolate | 485 ns | 366 ns | 1.33x | 1e-16 |
| `oned/linear/n1000/sorted` | construction | ndarray-interp | 1.99 µs | 1.51 µs | 1.32x | 0e+00 |
| `grid/linear/2d/max10/eval10/ascending` | evaluation | interpn | 137 ns | 104 ns | 1.32x | 1e-16 |
| `grid/linear/4d/max200/eval100/ascending` | construction | scirs2-interpolate | 51.8 ms | 39.7 ms | 1.31x | 1e-16 |
| `grid/linear/3d/max200/eval100/ascending` | construction | scirs2-interpolate | 1.86 ms | 1.44 ms | 1.30x | 3e-16 |
| `grid/linear/3d/max200/eval10000/ascending` | construction | scirs2-interpolate | 1.78 ms | 1.38 ms | 1.29x | 3e-16 |
| `grid/linear/2d/max200/eval10000/ascending` | construction | SciPy | 7.49 µs | 5.81 µs | 1.29x |  |
| `oned/linear/n10000/sorted-uneven` | construction | ndarray-interp | 19.5 µs | 15.2 µs | 1.28x | 0e+00 |
| `oned/linear/n1000/shuffled` | construction | ndarray-interp | 1.93 µs | 1.51 µs | 1.28x | 0e+00 |
| `grid/linear/2d/max200/eval10/ascending` | construction | SciPy | 7.36 µs | 5.78 µs | 1.27x |  |
| `grid/linear/4d/max40/eval100/ascending` | construction | SciPy | 13.1 µs | 10.3 µs | 1.27x |  |
| `grid/linear/2d/max200/eval1000/ascending` | construction | SciPy | 7.41 µs | 5.83 µs | 1.27x |  |
| `oned/linear/n100/shuffled-uneven` | construction | scirs2-interpolate | 265 ns | 210 ns | 1.26x | 3e-17 |
| `grid/linear/2d/max200/eval100/flipped` | construction | SciPy | 15.5 µs | 12.3 µs | 1.26x |  |
| `grid/linear/4d/max200/eval1000/ascending` | construction | scirs2-interpolate | 49.9 ms | 39.6 ms | 1.26x | 3e-16 |
| `grid/linear/3d/max10/eval100/ascending` | evaluation | interpn | 1.91 µs | 1.53 µs | 1.25x | 3e-16 |
| `grid/linear/2d/max40/eval10000/ascending` | evaluation | interpn | 166 µs | 133 µs | 1.25x | 4e-16 |
| `grid/nearest/2d/auckland-ba32-0401/heldout/ascending` | construction | scirs2-interpolate | 16.7 ms | 13.5 ms | 1.24x | 0e+00 |
| `oned/linear/n100/sorted-uneven` | construction | scirs2-interpolate | 242 ns | 195 ns | 1.24x | 3e-17 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | scirs2-interpolate | 1.04 ms | 844 µs | 1.23x | 2e-16 |
| `oned/linear/n100/shuffled` | construction | scirs2-interpolate | 241 ns | 196 ns | 1.23x | 0e+00 |
| `grid/linear/2d/auckland-ba32-0401/heldout/ascending` | construction | scirs2-interpolate | 16.5 ms | 13.6 ms | 1.22x | 0e+00 |
| `grid/linear/2d/max200/eval10000/flipped` | construction | SciPy | 16.1 µs | 13.2 µs | 1.22x |  |
| `oned/linear/n50/sorted` | construction | scirs2-interpolate | 150 ns | 124 ns | 1.22x | 0e+00 |
| `grid/linear/4d/max40/eval10/ascending` | evaluation | interpn | 327 ns | 271 ns | 1.21x | 2e-16 |
| `oned/linear/n50/shuffled` | construction | scirs2-interpolate | 151 ns | 126 ns | 1.21x | 0e+00 |
| `grid/linear/2d/max200/eval100/ascending` | evaluation | interpn | 2.05 µs | 1.72 µs | 1.19x | 1e-16 |
| `grid/linear/3d/max10/eval1000/ascending` | construction | scirs2-interpolate | 195 ns | 163 ns | 1.19x | 3e-16 |
| `grid/linear/2d/max200/eval1000/flipped` | construction | SciPy | 14.5 µs | 12.2 µs | 1.18x |  |
| `oned/cubic_spline/n10000/shuffled` | construction | scirs2-interpolate | 172 µs | 146 µs | 1.18x | 4e-05 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | scirs2-interpolate | 320 µs | 272 µs | 1.18x | 2e-16 |
| `grid/linear/2d/max200/eval10/flipped` | construction | SciPy | 14.5 µs | 12.3 µs | 1.17x |  |
| `grid/linear/2d/max200/eval100/ascending` | construction | SciPy | 6.83 µs | 5.85 µs | 1.17x |  |
| `grid/linear/2d/max40/eval1000/ascending` | evaluation | interpn | 15.4 µs | 13.3 µs | 1.15x | 3e-16 |
| `grid/linear/2d/max40/eval100/ascending` | evaluation | interpn | 1.66 µs | 1.45 µs | 1.14x | 2e-16 |
| `grid/linear/3d/max200/eval1000/ascending` | evaluation | interpn | 29.7 µs | 26.3 µs | 1.13x | 3e-16 |
| `grid/linear/2d/max200/eval10/ascending` | evaluation | ndarray-interp | 211 ns | 188 ns | 1.12x | 2e-16 |
| `grid/linear/3d/max10/eval100/ascending` | construction | scirs2-interpolate | 187 ns | 167 ns | 1.12x | 2e-16 |
| `grid/linear/4d/max200/eval100/ascending` | evaluation | interpn | 4.1 µs | 3.65 µs | 1.12x | 3e-16 |
| `grid/linear/2d/max40/eval100/ascending` | evaluation | ndarray-interp | 1.66 µs | 1.48 µs | 1.12x | 2e-16 |
| `grid/linear/4d/max40/eval100/ascending` | evaluation | interpn | 3.03 µs | 2.71 µs | 1.12x | 2e-16 |
| `grid/nearest/2d/each40/eval10000/ascending` | evaluation | interpn | 165 µs | 148 µs | 1.12x | 0e+00 |
| `grid/linear/3d/max200/eval100/ascending` | evaluation | interpn | 2.97 µs | 2.68 µs | 1.11x | 3e-16 |
| `grid/linear/4d/max200/eval1000/ascending` | evaluation | interpn | 40.7 µs | 37 µs | 1.10x | 3e-16 |

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

A run (one library on one case, a process of its own) is stopped at 15 minutes. The longest run that finished took 374 s, 42% of that: scirs2-interpolate on `linear_nd/n100000/image500`.

No run took more than half the limit (450 s).

## Real-world data

Two data sets, each taken whole in the unit its source gives it out in; what was taken and how it was divided was fixed before anything was measured (`compare/perf/real.py`, docs/design.ja.md section 8.4.7).

- `auckland-ba32-0401`: Auckland Part 1 LiDAR 1m DEM (2024), tile BA32_10000_0401. 4800 by 7200 cells of 1 m. The grid is every second row and column (and the last ones); every other cell is a query, and its elevation is the value held out. Licensed by Regional Software Holdings Limited for re-use under CC-BY-4.0. Hosted by Toitū Te Whenua Land Information New Zealand (https://registry.opendata.aws/nz-elevation).
- `ghcn-tmax-20250101`: GHCN-Daily, maximum temperature of 2025-01-01, every station; longitude and latitude in degrees, as they are. Of every 10 stations in the order of the file the last is a query, and its temperature is the value held out. The RBF kernels are those without a shape parameter. NOAA Global Historical Climatology Network Daily (GHCN-D), from https://registry.opendata.aws/noaa-ghcn (CC0 1.0).

These cases are also in the tables of their families below and in the counts above. Here they are with the error: the largest and the root mean square distance from the values held out (metres, and degrees Celsius), and the number of queries a library gave no value. A library missing from a case does not have it, or did not finish (see the runs that failed).

`grid/nearest/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 16.7 ms | 370 ms | 15.5 | 0.178 | 0 |
| interpx, parallel | 16.7 ms | 63.7 ms | 15.5 | 0.178 | 0 |
| interpx, `from_owned` | 5.26 ms | 373 ms | 15.5 | 0.178 | 0 |
| SciPy | 7.68 µs | 1.38 s | 15.5 | 0.178 | 0 |
| scirs2-interpolate | 13.5 ms | 159 s | 15.5 | 0.178 | 0 |
| interpn | 0.78 ns | 821 ms | 15.5 | 0.178 | 0 |

`grid/linear/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 16.5 ms | 378 ms | 6.74 | 0.0601 | 0 |
| interpx, parallel | 17 ms | 70.2 ms | 6.74 | 0.0601 | 0 |
| interpx, `from_owned` | 5.76 ms | 396 ms | 6.74 | 0.0601 | 0 |
| SciPy | 7.59 µs | 1.14 s | 6.74 | 0.0601 | 0 |
| scirs2-interpolate | 13.6 ms | 73.3 s | 6.74 | 0.0601 | 0 |
| ndarray-interp | 8.24 µs | 395 ms | 6.74 | 0.0601 | 0 |
| interpn | 0.68 ns | 825 ms | 6.74 | 0.0601 | 0 |

`grid/cubic/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 83.8 ms | 787 ms | 6.77 | 0.0507 | 0 |
| interpx, parallel | 81.1 ms | 133 ms | 6.77 | 0.0507 | 0 |
| interpx, `from_owned` | 76.8 ms | 802 ms | 6.77 | 0.0507 | 0 |
| SciPy | 15 s | 34.5 s | 6.77 | 0.0507 | 0 |
| interpn | 0.88 ns | 1.87 s | 6.76 | 0.0517 | 0 |

`grid/quintic/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 96.7 ms | 1.22 s | 6.77 | 0.0517 | 0 |
| interpx, parallel | 96.3 ms | 225 ms | 6.77 | 0.0517 | 0 |
| interpx, `from_owned` | 90.6 ms | 1.27 s | 6.77 | 0.0517 | 0 |
| SciPy | 21.5 s | 35 s | 6.77 | 0.0517 | 0 |

`grid/pchip/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 76 ms | 1.08 s | 6.75 | 0.049 | 0 |
| interpx, parallel | 75.9 ms | 188 ms | 6.75 | 0.049 | 0 |
| interpx, `from_owned` | 66.1 ms | 1.09 s | 6.75 | 0.049 | 0 |

`nearest_nd/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 558 µs | 172 µs | 25 | 3.42 | 0 |
| interpx, parallel | 560 µs | 73.3 µs | 25 | 3.42 | 0 |
| interpx, `from_owned` | 564 µs | 173 µs | 25 | 3.42 | 0 |
| SciPy | 1.87 ms | 558 µs | 25 | 3.42 | 0 |

`linear_nd/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 7.94 ms | 75.6 µs | 22.1 | 2.83 | 1 |
| interpx, parallel | 8.32 ms | 93.7 µs | 22.1 | 2.83 | 1 |
| interpx, `from_owned` | 7.86 ms | 79.5 µs | 22.1 | 2.83 | 1 |
| SciPy | 30 ms | 1.18 ms | 22.1 | 2.83 | 1 |
| scirs2-interpolate | 520 ms | 140 ms | 22.1 | 2.83 | 1 |
| spade | 4.9 ms | 1.32 ms | 22.1 | 2.83 | 1 |

`rbf/linear/global/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 20.6 s | 19.5 ms | 20.8 | 2.74 | 0 |
| interpx, parallel | 17.5 s | 3.7 ms | 20.8 | 2.74 | 0 |
| interpx, `from_owned` | 18 s | 18.5 ms | 20.8 | 2.74 | 0 |
| SciPy | 140 s | 55.4 ms | 20.8 | 2.74 | 0 |
| scirs2-interpolate | 294 s | 236 ms | 20.8 | 2.74 | 0 |

`rbf/linear/k50/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 774 µs | 19.4 ms | 20.8 | 2.75 | 0 |
| interpx, parallel | 828 µs | 3.51 ms | 20.8 | 2.75 | 0 |
| interpx, `from_owned` | 724 µs | 19.5 ms | 20.8 | 2.75 | 0 |
| SciPy | 1.73 ms | 44.4 ms | 20.8 | 2.75 | 0 |

`rbf/thin_plate_spline/global/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 18.4 s | 83 ms | 102 | 5.67 | 0 |
| interpx, parallel | 17.8 s | 16.1 ms | 102 | 5.67 | 0 |
| interpx, `from_owned` | 19 s | 83.2 ms | 102 | 5.67 | 0 |
| SciPy | 117 s | 113 ms | 102 | 5.67 | 0 |
| scirs2-interpolate | 300 s | 402 ms | 93.7 | 5.41 | 0 |

`rbf/thin_plate_spline/k50/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 725 µs | 30 ms | 102 | 5.66 | 0 |
| interpx, parallel | 754 µs | 5.63 ms | 102 | 5.66 | 0 |
| interpx, `from_owned` | 649 µs | 30.4 ms | 102 | 5.66 | 0 |
| SciPy | 1.84 ms | 58.5 ms | 102 | 5.66 | 0 |

`rbf/cubic/global/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 18 s | 18.5 ms | 587 | 27.1 | 0 |
| interpx, parallel | 17.9 s | 3.92 ms | 587 | 27.1 | 0 |
| interpx, `from_owned` | 17.8 s | 18 ms | 587 | 27.1 | 0 |
| SciPy | 124 s | 53.2 ms | 597 | 24.4 | 0 |
| scirs2-interpolate | 307 s | 268 ms | 100 | 8.48 | 0 |

`rbf/cubic/k50/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 752 µs | 21.5 ms | 567 | 25.3 | 0 |
| interpx, parallel | 821 µs | 4.11 ms | 567 | 25.3 | 0 |
| interpx, `from_owned` | 719 µs | 21 ms | 567 | 25.3 | 0 |
| SciPy | 1.86 ms | 49.8 ms | 567 | 25.3 | 0 |

`rbf/quintic/global/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 18.9 s | 21 ms | 2.55e+05 | 1.36e+04 | 0 |
| interpx, parallel | 17.8 s | 4.13 ms | 2.55e+05 | 1.36e+04 | 0 |
| interpx, `from_owned` | 18.6 s | 21.1 ms | 2.55e+05 | 1.36e+04 | 0 |
| SciPy | 122 s | 53.8 ms | 9.39e+04 | 4.92e+03 | 0 |
| scirs2-interpolate | 284 s | 294 ms | 6.56e+16 | 2.85e+15 | 0 |

`rbf/quintic/k50/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 828 µs | 22.6 ms | 3.74e+03 | 199 | 0 |
| interpx, parallel | 820 µs | 4.13 ms | 3.74e+03 | 199 | 0 |
| interpx, `from_owned` | 725 µs | 23.2 ms | 3.74e+03 | 199 | 0 |
| SciPy | 1.87 ms | 53.6 ms | 3.74e+03 | 200 | 0 |

## One variable

Data of `n` points, evaluated at `4 n` points: in order (`sorted`) and in no order (`shuffled`).

### Evaluation

| case | interpx | interpx, parallel | SciPy | scirs2-interpolate | ndarray-interp | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/n10/sorted` | 112 ns | 173 ns | 5.12 µs | 306 ns | 255 ns | 2.3x faster |
| `linear/n10/shuffled` | 187 ns | 243 ns | 5.08 µs | 284 ns | 250 ns | 1.3x faster |
| `nearest/n10/sorted` | 159 ns | 201 ns | 9.61 µs | 273 ns |  | 1.7x faster |
| `nearest/n10/shuffled` | 303 ns | 396 ns | 9.54 µs | 277 ns |  | **1.1x slower** |
| `previous/n10/sorted` | 155 ns | 202 ns | 10.8 µs |  |  | 70x faster |
| `previous/n10/shuffled` | 265 ns | 354 ns | 10.7 µs |  |  | 40x faster |
| `bspline1/n10/sorted` | 211 ns | 293 ns | 2.68 µs |  |  | 13x faster |
| `bspline1/n10/shuffled` | 527 ns | 654 ns | 2.75 µs |  |  | 5.2x faster |
| `bspline2/n10/sorted` | 268 ns | 345 ns | 3.08 µs |  |  | 11x faster |
| `bspline2/n10/shuffled` | 339 ns | 453 ns | 3.11 µs |  |  | 9.2x faster |
| `bspline3/n10/sorted` | 343 ns | 445 ns | 3.62 µs |  |  | 11x faster |
| `bspline3/n10/shuffled` | 440 ns | 520 ns | 3.63 µs |  |  | 8.2x faster |
| `cubic_spline/n10/sorted` | 149 ns | 203 ns | 3.7 µs | 277 ns | 333 ns | 1.9x faster |
| `cubic_spline/n10/shuffled` | 205 ns | 265 ns | 3.67 µs | 289 ns | 325 ns | 1.4x faster |
| `pchip/n10/sorted` | 141 ns | 192 ns | 3.74 µs | 320 ns |  | 2.3x faster |
| `pchip/n10/shuffled` | 204 ns | 265 ns | 3.7 µs | 310 ns |  | 1.5x faster |
| `akima/n10/sorted` | 152 ns | 196 ns | 4 µs | 262 ns |  | 1.7x faster |
| `akima/n10/shuffled` | 227 ns | 259 ns | 3.75 µs | 263 ns |  | 1.2x faster |
| `makima/n10/sorted` | 146 ns | 208 ns | 3.71 µs |  |  | 25x faster |
| `makima/n10/shuffled` | 205 ns | 256 ns | 4.02 µs |  |  | 20x faster |
| `linear/n50/sorted` | 632 ns | 739 ns | 6.2 µs | 1.73 µs | 1.23 µs | 2x faster |
| `linear/n50/shuffled` | 811 ns | 879 ns | 8.32 µs | 1.94 µs | 1.38 µs | 1.7x faster |
| `nearest/n50/sorted` | 828 ns | 829 ns | 12.5 µs | 1.83 µs |  | 2.2x faster |
| `nearest/n50/shuffled` | 1.56 µs | 1.7 µs | 11.7 µs | 1.76 µs |  | 1.1x faster |
| `previous/n50/sorted` | 841 ns | 884 ns | 12.8 µs |  |  | 15x faster |
| `previous/n50/shuffled` | 1.38 µs | 1.5 µs | 12.9 µs |  |  | 9.3x faster |
| `bspline1/n50/sorted` | 1.11 µs | 1.27 µs | 3.93 µs |  |  | 3.5x faster |
| `bspline1/n50/shuffled` | 2.8 µs | 3.14 µs | 4.86 µs |  |  | 1.7x faster |
| `bspline2/n50/sorted` | 1.36 µs | 1.58 µs | 6.68 µs |  |  | 4.9x faster |
| `bspline2/n50/shuffled` | 2.21 µs | 2.09 µs | 6.74 µs |  |  | 3x faster |
| `bspline3/n50/sorted` | 1.83 µs | 2.06 µs | 8.64 µs |  |  | 4.7x faster |
| `bspline3/n50/shuffled` | 2.81 µs | 2.89 µs | 10.2 µs |  |  | 3.6x faster |
| `cubic_spline/n50/sorted` | 718 ns | 758 ns | 5.15 µs | 4.46 µs | 1.61 µs | 2.2x faster |
| `cubic_spline/n50/shuffled` | 986 ns | 1.04 µs | 5.7 µs | 4.47 µs | 1.61 µs | 1.6x faster |
| `pchip/n50/sorted` | 763 ns | 770 ns | 5.32 µs | 4.69 µs |  | 6.1x faster |
| `pchip/n50/shuffled` | 979 ns | 1.04 µs | 5.8 µs | 4.25 µs |  | 4.3x faster |
| `akima/n50/sorted` | 712 ns | 766 ns | 5.26 µs | 4.5 µs |  | 6.3x faster |
| `akima/n50/shuffled` | 1.08 µs | 1.04 µs | 5.56 µs | 3.97 µs |  | 3.7x faster |
| `makima/n50/sorted` | 719 ns | 762 ns | 5.62 µs |  |  | 7.8x faster |
| `makima/n50/shuffled` | 977 ns | 1.05 µs | 6.01 µs |  |  | 6.1x faster |
| `linear/n100/sorted` | 1.32 µs | 1.38 µs | 7.95 µs | 3.66 µs | 2.48 µs | 1.9x faster |
| `linear/n100/shuffled` | 1.62 µs | 1.88 µs | 12.8 µs | 3.88 µs | 2.49 µs | 1.5x faster |
| `nearest/n100/sorted` | 1.59 µs | 1.72 µs | 14.8 µs | 3.97 µs |  | 2.5x faster |
| `nearest/n100/shuffled` | 3.09 µs | 3.28 µs | 14 µs | 3.8 µs |  | 1.2x faster |
| `previous/n100/sorted` | 1.58 µs | 1.7 µs | 15.1 µs |  |  | 9.5x faster |
| `previous/n100/shuffled` | 2.73 µs | 2.85 µs | 16.8 µs |  |  | 6.1x faster |
| `bspline1/n100/sorted` | 2.2 µs | 2.58 µs | 5.77 µs |  |  | 2.6x faster |
| `bspline1/n100/shuffled` | 5.61 µs | 5.77 µs | 11.8 µs |  |  | 2.1x faster |
| `bspline2/n100/sorted` | 2.71 µs | 3.41 µs | 9.78 µs |  |  | 3.6x faster |
| `bspline2/n100/shuffled` | 4.86 µs | 5.07 µs | 12.7 µs |  |  | 2.6x faster |
| `bspline3/n100/sorted` | 3.62 µs | 3.8 µs | 15 µs |  |  | 4.2x faster |
| `bspline3/n100/shuffled` | 6.73 µs | 6.28 µs | 19.1 µs |  |  | 2.8x faster |
| `cubic_spline/n100/sorted` | 1.47 µs | 1.55 µs | 7.47 µs | 18.8 µs | 3.19 µs | 2.2x faster |
| `cubic_spline/n100/shuffled` | 1.96 µs | 2.05 µs | 7.86 µs | 17.3 µs | 3.2 µs | 1.6x faster |
| `pchip/n100/sorted` | 1.47 µs | 1.66 µs | 6.91 µs | 17.6 µs |  | 4.7x faster |
| `pchip/n100/shuffled` | 2.1 µs | 2.07 µs | 7.99 µs | 16.6 µs |  | 3.8x faster |
| `akima/n100/sorted` | 1.47 µs | 1.54 µs | 6.84 µs | 18.3 µs |  | 4.7x faster |
| `akima/n100/shuffled` | 1.96 µs | 2.12 µs | 7.87 µs | 17.5 µs |  | 4x faster |
| `makima/n100/sorted` | 1.47 µs | 1.53 µs | 6.94 µs |  |  | 4.7x faster |
| `makima/n100/shuffled` | 1.96 µs | 2.05 µs | 8.05 µs |  |  | 4.1x faster |
| `linear/n1000/sorted` | 12.6 µs | 13.5 µs | 29.9 µs | 115 µs | 24.8 µs | 2x faster |
| `linear/n1000/shuffled` | 16.1 µs | 18 µs | 143 µs | 220 µs | 24.7 µs | 1.5x faster |
| `nearest/n1000/sorted` | 15 µs | 16.9 µs | 60.3 µs | 162 µs |  | 4x faster |
| `nearest/n1000/shuffled` | 31 µs | 31.3 µs | 201 µs | 251 µs |  | 6.5x faster |
| `previous/n1000/sorted` | 16.6 µs | 15.7 µs | 59.6 µs |  |  | 3.6x faster |
| `previous/n1000/shuffled` | 28.6 µs | 27.7 µs | 206 µs |  |  | 7.2x faster |
| `bspline1/n1000/sorted` | 25.5 µs | 25.6 µs | 31.4 µs |  |  | 1.2x faster |
| `bspline1/n1000/shuffled` | 58.3 µs | 62.1 µs | 530 µs |  |  | 9.1x faster |
| `bspline2/n1000/sorted` | 33.7 µs | 32.1 µs | 77.9 µs |  |  | 2.3x faster |
| `bspline2/n1000/shuffled` | 63.7 µs | 63.8 µs | 543 µs |  |  | 8.5x faster |
| `bspline3/n1000/sorted` | 35.8 µs | 39.6 µs | 129 µs |  |  | 3.6x faster |
| `bspline3/n1000/shuffled` | 86.2 µs | 86.9 µs | 631 µs |  |  | 7.3x faster |
| `cubic_spline/n1000/sorted` | 14.3 µs | 15.1 µs | 38.2 µs | 1.43 ms | 31.5 µs | 2.2x faster |
| `cubic_spline/n1000/shuffled` | 19.8 µs | 21.4 µs | 220 µs | 1.45 ms | 34.9 µs | 1.8x faster |
| `pchip/n1000/sorted` | 14.6 µs | 15.2 µs | 37 µs | 1.2 ms |  | 2.5x faster |
| `pchip/n1000/shuffled` | 19.9 µs | 21.9 µs | 220 µs | 1.22 ms |  | 11x faster |
| `akima/n1000/sorted` | 15.6 µs | 15.1 µs | 39 µs | 1.46 ms |  | 2.5x faster |
| `akima/n1000/shuffled` | 20 µs | 22.1 µs | 222 µs | 1.42 ms |  | 11x faster |
| `makima/n1000/sorted` | 14.4 µs | 16.2 µs | 37.3 µs |  |  | 2.6x faster |
| `makima/n1000/shuffled` | 20.1 µs | 22.1 µs | 219 µs |  |  | 11x faster |
| `linear/n10000/sorted` | 134 µs | 65.7 µs | 235 µs | 1.54 ms | 248 µs | 1.8x faster |
| `linear/n10000/shuffled` | 208 µs | 88.5 µs | 2.33 ms | 2.86 ms | 285 µs | 1.4x faster |
| `nearest/n10000/sorted` | 174 µs | 76.7 µs | 625 µs | 1.96 ms |  | 3.6x faster |
| `nearest/n10000/shuffled` | 384 µs | 96.6 µs | 3.22 ms | 3.05 ms |  | 8x faster |
| `previous/n10000/sorted` | 163 µs | 77.7 µs | 663 µs |  |  | 4.1x faster |
| `previous/n10000/shuffled` | 342 µs | 89.5 µs | 2.98 ms |  |  | 8.7x faster |
| `bspline1/n10000/sorted` | 256 µs | 97.2 µs | 293 µs |  |  | 1.1x faster |
| `bspline1/n10000/shuffled` | 639 µs | 147 µs | 50.3 ms |  |  | 79x faster |
| `bspline2/n10000/sorted` | 324 µs | 92.2 µs | 759 µs |  |  | 2.3x faster |
| `bspline2/n10000/shuffled` | 969 µs | 212 µs | 50.3 ms |  |  | 52x faster |
| `bspline3/n10000/sorted` | 372 µs | 107 µs | 1.27 ms |  |  | 3.4x faster |
| `bspline3/n10000/shuffled` | 1.15 ms | 261 µs | 49.2 ms |  |  | 43x faster |
| `cubic_spline/n10000/sorted` | 150 µs | 75.2 µs | 391 µs | 147 ms | 322 µs | 2.1x faster |
| `cubic_spline/n10000/shuffled` | 239 µs | 94.5 µs | 3.62 ms | 149 ms | 332 µs | 1.4x faster |
| `pchip/n10000/sorted` | 151 µs | 71 µs | 394 µs | 127 ms |  | 2.6x faster |
| `pchip/n10000/shuffled` | 218 µs | 97.4 µs | 3.59 ms | 130 ms |  | 16x faster |
| `akima/n10000/sorted` | 159 µs | 71.3 µs | 403 µs | 142 ms |  | 2.5x faster |
| `akima/n10000/shuffled` | 238 µs | 94.4 µs | 3.64 ms | 154 ms |  | 15x faster |
| `makima/n10000/sorted` | 152 µs | 70.4 µs | 391 µs |  |  | 2.6x faster |
| `makima/n10000/shuffled` | 218 µs | 95.6 µs | 3.62 ms |  |  | 17x faster |
| `linear/n100/sorted-uneven` | 1.19 µs | 1.4 µs | 8.08 µs | 3.66 µs | 3.76 µs | 3.1x faster |
| `linear/n100/shuffled-uneven` | 2.73 µs | 2.98 µs | 12.6 µs | 4.23 µs | 4.06 µs | 1.5x faster |
| `cubic_spline/n100/sorted-uneven` | 1.36 µs | 1.46 µs | 8.06 µs | 17.5 µs | 4.42 µs | 3.2x faster |
| `cubic_spline/n100/shuffled-uneven` | 2.94 µs | 3.35 µs | 8.75 µs | 19 µs | 4.4 µs | 1.5x faster |
| `pchip/n100/sorted-uneven` | 1.34 µs | 1.48 µs | 7.39 µs | 16.1 µs |  | 5.5x faster |
| `pchip/n100/shuffled-uneven` | 2.95 µs | 3.34 µs | 8.06 µs | 16.5 µs |  | 2.7x faster |
| `linear/n10000/sorted-uneven` | 204 µs | 87.1 µs | 280 µs | 1.59 ms | 972 µs | 1.4x faster |
| `linear/n10000/shuffled-uneven` | 678 µs | 162 µs | 2.35 ms | 2.89 ms | 1.91 ms | 2.8x faster |
| `cubic_spline/n10000/sorted-uneven` | 221 µs | 95.5 µs | 440 µs | 149 ms | 964 µs | 2x faster |
| `cubic_spline/n10000/shuffled-uneven` | 796 µs | 196 µs | 3.64 ms | 152 ms | 1.98 ms | 2.5x faster |
| `pchip/n10000/sorted-uneven` | 242 µs | 94.1 µs | 464 µs | 127 ms |  | 1.9x faster |
| `pchip/n10000/shuffled-uneven` | 794 µs | 178 µs | 3.64 ms | 125 ms |  | 4.6x faster |

### Construction

| case | interpx | interpx, parallel | SciPy | scirs2-interpolate | ndarray-interp | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/n10/sorted` | 86.5 ns | 86.7 ns | 8.91 µs | 79.9 ns | 41.4 ns | **2.1x slower** |
| `linear/n10/shuffled` | 86.6 ns | 95 ns | 8.84 µs | 81.1 ns | 41.5 ns | **2.1x slower** |
| `nearest/n10/sorted` | 77 ns | 75.1 ns | 9.61 µs | 80.4 ns |  | 1x faster |
| `nearest/n10/shuffled` | 75.8 ns | 82.1 ns | 9.85 µs | 79.5 ns |  | 1x faster |
| `previous/n10/sorted` | 80.9 ns | 76.2 ns | 9.54 µs |  |  | 1.2e+02x faster |
| `previous/n10/shuffled` | 77.1 ns | 75.6 ns | 9.26 µs |  |  | 1.2e+02x faster |
| `bspline1/n10/sorted` | 105 ns | 106 ns | 17.1 µs |  |  | 1.6e+02x faster |
| `bspline1/n10/shuffled` | 104 ns | 105 ns | 15.7 µs |  |  | 1.5e+02x faster |
| `bspline2/n10/sorted` | 483 ns | 503 ns | 36.6 µs |  |  | 76x faster |
| `bspline2/n10/shuffled` | 490 ns | 517 ns | 40.5 µs |  |  | 83x faster |
| `bspline3/n10/sorted` | 541 ns | 561 ns | 36.5 µs |  |  | 67x faster |
| `bspline3/n10/shuffled` | 526 ns | 536 ns | 36.4 µs |  |  | 69x faster |
| `cubic_spline/n10/sorted` | 251 ns | 255 ns | 85.6 µs | 391 ns | 468 ns | 1.6x faster |
| `cubic_spline/n10/shuffled` | 231 ns | 236 ns | 85.5 µs | 384 ns | 467 ns | 1.7x faster |
| `pchip/n10/sorted` | 174 ns | 177 ns | 109 µs | 191 ns |  | 1.1x faster |
| `pchip/n10/shuffled` | 192 ns | 176 ns | 103 µs | 190 ns |  | **1x slower** |
| `akima/n10/sorted` | 196 ns | 178 ns | 83.1 µs | 195 ns |  | **1x slower** |
| `akima/n10/shuffled` | 196 ns | 194 ns | 77.4 µs | 186 ns |  | **1.1x slower** |
| `makima/n10/sorted` | 191 ns | 192 ns | 82.8 µs |  |  | 4.3e+02x faster |
| `makima/n10/shuffled` | 207 ns | 188 ns | 84.8 µs |  |  | 4.1e+02x faster |
| `linear/n50/sorted` | 150 ns | 163 ns | 9.08 µs | 124 ns | 99.9 ns | **1.5x slower** |
| `linear/n50/shuffled` | 151 ns | 165 ns | 9.02 µs | 126 ns | 99.8 ns | **1.5x slower** |
| `nearest/n50/sorted` | 107 ns | 107 ns | 10.7 µs | 137 ns |  | 1.3x faster |
| `nearest/n50/shuffled` | 108 ns | 117 ns | 9.86 µs | 132 ns |  | 1.2x faster |
| `previous/n50/sorted` | 116 ns | 115 ns | 10.5 µs |  |  | 91x faster |
| `previous/n50/shuffled` | 116 ns | 106 ns | 9.72 µs |  |  | 84x faster |
| `bspline1/n50/sorted` | 144 ns | 158 ns | 15.9 µs |  |  | 1.1e+02x faster |
| `bspline1/n50/shuffled` | 154 ns | 145 ns | 16 µs |  |  | 1e+02x faster |
| `bspline2/n50/sorted` | 1.83 µs | 1.84 µs | 40.6 µs |  |  | 22x faster |
| `bspline2/n50/shuffled` | 1.84 µs | 1.86 µs | 40.8 µs |  |  | 22x faster |
| `bspline3/n50/sorted` | 2.05 µs | 2.24 µs | 39.5 µs |  |  | 19x faster |
| `bspline3/n50/shuffled` | 2.24 µs | 2.04 µs | 39.3 µs |  |  | 18x faster |
| `cubic_spline/n50/sorted` | 831 ns | 836 ns | 87.1 µs | 971 ns | 1.08 µs | 1.2x faster |
| `cubic_spline/n50/shuffled` | 832 ns | 836 ns | 94.4 µs | 938 ns | 1.14 µs | 1.1x faster |
| `pchip/n50/sorted` | 501 ns | 500 ns | 112 µs | 364 ns |  | **1.4x slower** |
| `pchip/n50/shuffled` | 485 ns | 487 ns | 103 µs | 366 ns |  | **1.3x slower** |
| `akima/n50/sorted` | 510 ns | 545 ns | 83 µs | 346 ns |  | **1.5x slower** |
| `akima/n50/shuffled` | 517 ns | 501 ns | 82.7 µs | 359 ns |  | **1.4x slower** |
| `makima/n50/sorted` | 585 ns | 540 ns | 81.1 µs |  |  | 1.4e+02x faster |
| `makima/n50/shuffled` | 538 ns | 589 ns | 85.4 µs |  |  | 1.6e+02x faster |
| `linear/n100/sorted` | 264 ns | 264 ns | 10.1 µs | 196 ns | 166 ns | **1.6x slower** |
| `linear/n100/shuffled` | 241 ns | 241 ns | 9.18 µs | 196 ns | 166 ns | **1.4x slower** |
| `nearest/n100/sorted` | 174 ns | 165 ns | 10 µs | 196 ns |  | 1.1x faster |
| `nearest/n100/shuffled` | 162 ns | 161 ns | 10.5 µs | 195 ns |  | 1.2x faster |
| `previous/n100/sorted` | 161 ns | 166 ns | 10.1 µs |  |  | 63x faster |
| `previous/n100/shuffled` | 161 ns | 161 ns | 11 µs |  |  | 68x faster |
| `bspline1/n100/sorted` | 215 ns | 197 ns | 15.9 µs |  |  | 74x faster |
| `bspline1/n100/shuffled` | 216 ns | 197 ns | 17.3 µs |  |  | 80x faster |
| `bspline2/n100/sorted` | 3.8 µs | 3.84 µs | 43 µs |  |  | 11x faster |
| `bspline2/n100/shuffled` | 3.85 µs | 3.83 µs | 45.8 µs |  |  | 12x faster |
| `bspline3/n100/sorted` | 4.29 µs | 3.92 µs | 44.2 µs |  |  | 10x faster |
| `bspline3/n100/shuffled` | 4.3 µs | 4.29 µs | 47.3 µs |  |  | 11x faster |
| `cubic_spline/n100/sorted` | 1.73 µs | 1.62 µs | 96.7 µs | 1.8 µs | 1.89 µs | 1x faster |
| `cubic_spline/n100/shuffled` | 1.58 µs | 1.59 µs | 88 µs | 1.8 µs | 1.9 µs | 1.1x faster |
| `pchip/n100/sorted` | 883 ns | 972 ns | 102 µs | 591 ns |  | **1.5x slower** |
| `pchip/n100/shuffled` | 970 ns | 892 ns | 102 µs | 551 ns |  | **1.8x slower** |
| `akima/n100/sorted` | 887 ns | 965 ns | 76.6 µs | 561 ns |  | **1.6x slower** |
| `akima/n100/shuffled` | 968 ns | 971 ns | 76.8 µs | 529 ns |  | **1.8x slower** |
| `makima/n100/sorted` | 980 ns | 973 ns | 85.8 µs |  |  | 88x faster |
| `makima/n100/shuffled` | 973 ns | 975 ns | 87.1 µs |  |  | 90x faster |
| `linear/n1000/sorted` | 1.99 µs | 2.02 µs | 12.8 µs | 1.23 µs | 1.51 µs | **1.6x slower** |
| `linear/n1000/shuffled` | 1.93 µs | 1.93 µs | 12.6 µs | 1.24 µs | 1.51 µs | **1.6x slower** |
| `nearest/n1000/sorted` | 1 µs | 1.07 µs | 15.4 µs | 1.22 µs |  | 1.2x faster |
| `nearest/n1000/shuffled` | 994 ns | 1.06 µs | 15.5 µs | 1.28 µs |  | 1.3x faster |
| `previous/n1000/sorted` | 985 ns | 1.02 µs | 16.6 µs |  |  | 17x faster |
| `previous/n1000/shuffled` | 1.03 µs | 1.03 µs | 17.8 µs |  |  | 17x faster |
| `bspline1/n1000/sorted` | 1.09 µs | 1.05 µs | 17.1 µs |  |  | 16x faster |
| `bspline1/n1000/shuffled` | 1.06 µs | 1.14 µs | 18.6 µs |  |  | 18x faster |
| `bspline2/n1000/sorted` | 35.2 µs | 33.5 µs | 116 µs |  |  | 3.3x faster |
| `bspline2/n1000/shuffled` | 36.6 µs | 36.9 µs | 108 µs |  |  | 3x faster |
| `bspline3/n1000/sorted` | 41.9 µs | 42.2 µs | 133 µs |  |  | 3.2x faster |
| `bspline3/n1000/shuffled` | 38.5 µs | 42.1 µs | 136 µs |  |  | 3.5x faster |
| `cubic_spline/n1000/sorted` | 15.5 µs | 15.6 µs | 112 µs | 14.5 µs | 15.7 µs | **1.1x slower** |
| `cubic_spline/n1000/shuffled` | 15.4 µs | 15.5 µs | 121 µs | 14.5 µs | 15.7 µs | **1.1x slower** |
| `pchip/n1000/sorted` | 8.27 µs | 9.06 µs | 118 µs | 4.55 µs |  | **1.8x slower** |
| `pchip/n1000/shuffled` | 8.26 µs | 9.07 µs | 118 µs | 4.18 µs |  | **2x slower** |
| `akima/n1000/sorted` | 8.33 µs | 8.44 µs | 98.3 µs | 4.81 µs |  | **1.7x slower** |
| `akima/n1000/shuffled` | 8.36 µs | 9.09 µs | 94.5 µs | 4.39 µs |  | **1.9x slower** |
| `makima/n1000/sorted` | 9.36 µs | 9.95 µs | 98.7 µs |  |  | 11x faster |
| `makima/n1000/shuffled` | 9.15 µs | 9.96 µs | 101 µs |  |  | 11x faster |
| `linear/n10000/sorted` | 19.8 µs | 21.5 µs | 45.9 µs | 12.4 µs | 19.2 µs | **1.6x slower** |
| `linear/n10000/shuffled` | 21.3 µs | 21.7 µs | 47.4 µs | 12.5 µs | 21.2 µs | **1.7x slower** |
| `nearest/n10000/sorted` | 11.6 µs | 11.4 µs | 58.4 µs | 12.9 µs |  | 1.1x faster |
| `nearest/n10000/shuffled` | 10.7 µs | 11.8 µs | 54.5 µs | 12.4 µs |  | 1.2x faster |
| `previous/n10000/sorted` | 10.8 µs | 11.1 µs | 75.6 µs |  |  | 7x faster |
| `previous/n10000/shuffled` | 10.9 µs | 11 µs | 79 µs |  |  | 7.3x faster |
| `bspline1/n10000/sorted` | 10.8 µs | 11.7 µs | 24.8 µs |  |  | 2.3x faster |
| `bspline1/n10000/shuffled` | 10.6 µs | 11.5 µs | 25.1 µs |  |  | 2.4x faster |
| `bspline2/n10000/sorted` | 448 µs | 399 µs | 752 µs |  |  | 1.7x faster |
| `bspline2/n10000/shuffled` | 517 µs | 490 µs | 782 µs |  |  | 1.5x faster |
| `bspline3/n10000/sorted` | 473 µs | 538 µs | 993 µs |  |  | 2.1x faster |
| `bspline3/n10000/shuffled` | 504 µs | 567 µs | 937 µs |  |  | 1.9x faster |
| `cubic_spline/n10000/sorted` | 171 µs | 176 µs | 642 µs | 158 µs | 183 µs | **1.1x slower** |
| `cubic_spline/n10000/shuffled` | 172 µs | 173 µs | 717 µs | 146 µs | 183 µs | **1.2x slower** |
| `pchip/n10000/sorted` | 124 µs | 95.6 µs | 258 µs | 46.7 µs |  | **2.7x slower** |
| `pchip/n10000/shuffled` | 83.5 µs | 91.1 µs | 255 µs | 42.8 µs |  | **2x slower** |
| `akima/n10000/sorted` | 120 µs | 91 µs | 267 µs | 45.4 µs |  | **2.6x slower** |
| `akima/n10000/shuffled` | 93.8 µs | 96.3 µs | 252 µs | 44.9 µs |  | **2.1x slower** |
| `makima/n10000/sorted` | 137 µs | 99.9 µs | 708 µs |  |  | 5.2x faster |
| `makima/n10000/shuffled` | 137 µs | 100 µs | 275 µs |  |  | 2x faster |
| `linear/n100/sorted-uneven` | 242 ns | 244 ns | 10.1 µs | 195 ns | 167 ns | **1.5x slower** |
| `linear/n100/shuffled-uneven` | 265 ns | 249 ns | 9.27 µs | 210 ns | 182 ns | **1.5x slower** |
| `cubic_spline/n100/sorted-uneven` | 1.58 µs | 1.58 µs | 95.3 µs | 1.67 µs | 1.89 µs | 1.1x faster |
| `cubic_spline/n100/shuffled-uneven` | 1.59 µs | 1.74 µs | 89 µs | 1.65 µs | 2.04 µs | 1x faster |
| `pchip/n100/sorted-uneven` | 883 ns | 885 ns | 103 µs | 580 ns |  | **1.5x slower** |
| `pchip/n100/shuffled-uneven` | 884 ns | 899 ns | 105 µs | 626 ns |  | **1.4x slower** |
| `linear/n10000/sorted-uneven` | 19.5 µs | 21.7 µs | 48.3 µs | 12.9 µs | 15.2 µs | **1.5x slower** |
| `linear/n10000/shuffled-uneven` | 22.9 µs | 19.9 µs | 48.5 µs | 12.9 µs | 15.2 µs | **1.8x slower** |
| `cubic_spline/n10000/sorted-uneven` | 232 µs | 172 µs | 622 µs | 158 µs | 167 µs | **1.5x slower** |
| `cubic_spline/n10000/shuffled-uneven` | 156 µs | 174 µs | 631 µs | 159 µs | 183 µs | 1x faster |
| `pchip/n10000/sorted-uneven` | 126 µs | 90.9 µs | 256 µs | 44.4 µs |  | **2.8x slower** |
| `pchip/n10000/shuffled-uneven` | 91.9 µs | 92.9 µs | 243 µs | 42.7 µs |  | **2.2x slower** |

## `PPoly`

A cubic of 55 intervals, evaluated at points in no order.

### Evaluation

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `intervals55/eval10` | 79.9 ns | 76.3 ns | 3.78 µs | 47x faster |
| `intervals55/eval100` | 646 ns | 793 ns | 4.48 µs | 6.9x faster |
| `intervals55/eval1000` | 6.48 µs | 6.7 µs | 16.3 µs | 2.5x faster |
| `intervals55/eval10000` | 66.6 µs | 67.3 µs | 332 µs | 5x faster |

### Construction

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `intervals55/eval10` | 174 ns | 172 ns | 5.69 µs | 33x faster |
| `intervals55/eval100` | 161 ns | 174 ns | 5.17 µs | 32x faster |
| `intervals55/eval1000` | 160 ns | 159 ns | 5.72 µs | 36x faster |
| `intervals55/eval10000` | 159 ns | 159 ns | 5.21 µs | 33x faster |

## `BSpline`

A cubic B-spline of 55 intervals, evaluated at points in no order.

### Evaluation

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `intervals55/eval10` | 159 ns | 148 ns | 2.95 µs | 18x faster |
| `intervals55/eval100` | 1.42 µs | 1.64 µs | 5.91 µs | 4.1x faster |
| `intervals55/eval1000` | 14.2 µs | 14.4 µs | 38.8 µs | 2.7x faster |
| `intervals55/eval10000` | 143 µs | 57.8 µs | 397 µs | 2.8x faster |

### Construction

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `intervals55/eval10` | 175 ns | 175 ns | 10.3 µs | 59x faster |
| `intervals55/eval100` | 164 ns | 159 ns | 9.58 µs | 58x faster |
| `intervals55/eval1000` | 159 ns | 164 ns | 9.31 µs | 58x faster |
| `intervals55/eval10000` | 159 ns | 159 ns | 10.2 µs | 64x faster |

## `BSplineCurve`

A cubic curve through `n` points of the plane, evaluated at `4 n` parameters in order.

### Evaluation

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `degree3/n100` | 4.42 µs | 4.68 µs | 15.9 µs | 3.6x faster |
| `degree3/n10000` | 490 µs | 155 µs | 1.45 ms | 3x faster |

### Construction

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `degree3/n100` | 4.67 µs | 4.67 µs | 48.2 µs | 10x faster |
| `degree3/n10000` | 687 µs | 690 µs | 1.38 ms | 2x faster |

## Regular grid

Points in no order inside the grid. `maxN`: the first axis has `N` points and each axis half the one before (SciPy's suite); `eachN`: every axis has `N` points; `flipped`: every axis given from its high end.

### Evaluation

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | ndarray-interp | interpn | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/2d/max10/eval10/ascending` | 137 ns | 152 ns | 137 ns | 19.3 µs | 793 ns | 153 ns | 104 ns | **1.3x slower** |
| `linear/2d/max10/eval10/flipped` | 138 ns | 153 ns | 137 ns | 21 µs |  |  |  | 1.5e+02x faster |
| `linear/2d/max10/eval100/ascending` | 1.38 µs | 1.31 µs | 1.37 µs | 23.5 µs | 8.49 µs | 1.49 µs | 935 ns | **1.5x slower** |
| `linear/2d/max10/eval100/flipped` | 1.37 µs | 1.31 µs | 1.37 µs | 21.7 µs |  |  |  | 16x faster |
| `linear/2d/max10/eval1000/ascending` | 13.7 µs | 13 µs | 12.6 µs | 40.3 µs | 92 µs | 14.2 µs | 9.24 µs | **1.5x slower** |
| `linear/2d/max10/eval1000/flipped` | 12.5 µs | 13.8 µs | 12.5 µs | 40.2 µs |  |  |  | 3.2x faster |
| `linear/2d/max10/eval10000/ascending` | 136 µs | 51.1 µs | 125 µs | 361 µs | 876 µs | 240 µs | 92.5 µs | **1.5x slower** |
| `linear/2d/max10/eval10000/flipped` | 125 µs | 54.8 µs | 125 µs | 358 µs |  |  |  | 2.9x faster |
| `linear/2d/max40/eval10/ascending` | 158 ns | 161 ns | 173 ns | 21.1 µs | 1 µs | 157 ns | 147 ns | **1.1x slower** |
| `linear/2d/max40/eval10/flipped` | 157 ns | 176 ns | 159 ns | 19.5 µs |  |  |  | 1.2e+02x faster |
| `linear/2d/max40/eval100/ascending` | 1.66 µs | 1.58 µs | 1.52 µs | 21.8 µs | 10.1 µs | 1.48 µs | 1.45 µs | **1.1x slower** |
| `linear/2d/max40/eval100/flipped` | 1.51 µs | 1.74 µs | 1.52 µs | 24 µs |  |  |  | 16x faster |
| `linear/2d/max40/eval1000/ascending` | 15.4 µs | 15.3 µs | 16.5 µs | 50.2 µs | 106 µs | 15.7 µs | 13.3 µs | **1.2x slower** |
| `linear/2d/max40/eval1000/flipped` | 15.1 µs | 15.3 µs | 16.5 µs | 46.5 µs |  |  |  | 3.1x faster |
| `linear/2d/max40/eval10000/ascending` | 166 µs | 62.7 µs | 153 µs | 567 µs | 1.12 ms | 224 µs | 133 µs | **1.2x slower** |
| `linear/2d/max40/eval10000/flipped` | 151 µs | 64.2 µs | 151 µs | 583 µs |  |  |  | 3.8x faster |
| `linear/2d/max200/eval10/ascending` | 211 ns | 197 ns | 212 ns | 19.4 µs | 1.97 µs | 188 ns | 193 ns | **1.1x slower** |
| `linear/2d/max200/eval10/flipped` | 193 ns | 200 ns | 193 ns | 19.7 µs |  |  |  | 1e+02x faster |
| `linear/2d/max200/eval100/ascending` | 2.05 µs | 1.95 µs | 2.05 µs | 22.9 µs | 19.9 µs | 2.17 µs | 1.72 µs | **1.2x slower** |
| `linear/2d/max200/eval100/flipped` | 2.05 µs | 1.98 µs | 1.88 µs | 23.3 µs |  |  |  | 11x faster |
| `linear/2d/max200/eval1000/ascending` | 18.8 µs | 19.7 µs | 20.6 µs | 73.9 µs | 209 µs | 24 µs | 18.9 µs | 1x faster |
| `linear/2d/max200/eval1000/flipped` | 20.5 µs | 20.7 µs | 19 µs | 72.4 µs |  |  |  | 3.5x faster |
| `linear/2d/max200/eval10000/ascending` | 200 µs | 80.2 µs | 207 µs | 811 µs | 2.32 ms | 454 µs | 188 µs | **1.1x slower** |
| `linear/2d/max200/eval10000/flipped` | 206 µs | 77.6 µs | 206 µs | 891 µs |  |  |  | 4.3x faster |
| `linear/3d/max10/eval10/ascending` | 216 ns | 224 ns | 202 ns | 62.7 µs | 3.38 µs |  | 155 ns | **1.4x slower** |
| `linear/3d/max10/eval10/flipped` | 217 ns | 225 ns | 219 ns | 63.4 µs |  |  |  | 2.9e+02x faster |
| `linear/3d/max10/eval100/ascending` | 1.91 µs | 2.08 µs | 1.93 µs | 66.6 µs | 31.6 µs |  | 1.53 µs | **1.3x slower** |
| `linear/3d/max10/eval100/flipped` | 1.9 µs | 1.96 µs | 2.08 µs | 66.8 µs |  |  |  | 35x faster |
| `linear/3d/max10/eval1000/ascending` | 20.8 µs | 20 µs | 19 µs | 126 µs | 317 µs |  | 15.2 µs | **1.4x slower** |
| `linear/3d/max10/eval1000/flipped` | 19 µs | 20.8 µs | 20.7 µs | 115 µs |  |  |  | 6.1x faster |
| `linear/3d/max10/eval10000/ascending` | 207 µs | 70.9 µs | 209 µs | 1.12 ms | 3.46 ms |  | 142 µs | **1.5x slower** |
| `linear/3d/max10/eval10000/flipped` | 191 µs | 75.8 µs | 194 µs | 1.16 ms |  |  |  | 6x faster |
| `linear/3d/max40/eval10/ascending` | 219 ns | 246 ns | 220 ns | 62.3 µs | 3.34 µs |  | 202 ns | **1.1x slower** |
| `linear/3d/max40/eval10/flipped` | 220 ns | 246 ns | 241 ns | 57.2 µs |  |  |  | 2.6e+02x faster |
| `linear/3d/max40/eval100/ascending` | 2.13 µs | 2.21 µs | 2.13 µs | 64.2 µs | 36.5 µs |  | 2.01 µs | **1.1x slower** |
| `linear/3d/max40/eval100/flipped` | 2.13 µs | 2.21 µs | 2.13 µs | 69.7 µs |  |  |  | 33x faster |
| `linear/3d/max40/eval1000/ascending` | 23.3 µs | 23.5 µs | 23.3 µs | 146 µs | 386 µs |  | 22 µs | **1.1x slower** |
| `linear/3d/max40/eval1000/flipped` | 23.3 µs | 21.5 µs | 23.3 µs | 145 µs |  |  |  | 6.2x faster |
| `linear/3d/max40/eval10000/ascending` | 233 µs | 76 µs | 233 µs | 1.55 ms | 3.63 ms |  | 219 µs | **1.1x slower** |
| `linear/3d/max40/eval10000/flipped` | 233 µs | 83 µs | 233 µs | 1.59 ms |  |  |  | 6.9x faster |
| `linear/3d/max200/eval10/ascending` | 302 ns | 308 ns | 302 ns | 63 µs | 4.91 µs |  | 286 ns | **1.1x slower** |
| `linear/3d/max200/eval10/flipped` | 302 ns | 305 ns | 301 ns | 62.6 µs |  |  |  | 2.1e+02x faster |
| `linear/3d/max200/eval100/ascending` | 2.97 µs | 2.87 µs | 2.97 µs | 65.1 µs | 51.3 µs |  | 2.68 µs | **1.1x slower** |
| `linear/3d/max200/eval100/flipped` | 2.97 µs | 3.05 µs | 2.97 µs | 74 µs |  |  |  | 25x faster |
| `linear/3d/max200/eval1000/ascending` | 29.7 µs | 27.3 µs | 29.7 µs | 198 µs | 495 µs |  | 26.3 µs | **1.1x slower** |
| `linear/3d/max200/eval1000/flipped` | 29.7 µs | 28.1 µs | 29.7 µs | 187 µs |  |  |  | 6.3x faster |
| `linear/3d/max200/eval10000/ascending` | 297 µs | 92 µs | 297 µs | 1.96 ms | 5.15 ms |  | 289 µs | **1x slower** |
| `linear/3d/max200/eval10000/flipped` | 297 µs | 93.2 µs | 273 µs | 2.07 ms |  |  |  | 7x faster |
| `linear/4d/max10/eval10/ascending` | 284 ns | 283 ns | 272 ns | 106 µs |  |  |  | 3.7e+02x faster |
| `linear/4d/max10/eval10/flipped` | 274 ns | 287 ns | 271 ns | 104 µs |  |  |  | 3.8e+02x faster |
| `linear/4d/max10/eval100/ascending` | 2.64 µs | 2.7 µs | 2.64 µs | 116 µs |  |  |  | 44x faster |
| `linear/4d/max10/eval100/flipped` | 2.62 µs | 2.72 µs | 2.64 µs | 107 µs |  |  |  | 41x faster |
| `linear/4d/max10/eval1000/ascending` | 26 µs | 24.4 µs | 26.7 µs | 218 µs |  |  |  | 8.4x faster |
| `linear/4d/max10/eval1000/flipped` | 25.4 µs | 25 µs | 26.7 µs | 225 µs |  |  |  | 8.9x faster |
| `linear/4d/max10/eval10000/ascending` | 245 µs | 89.2 µs | 265 µs | 1.82 ms |  |  |  | 7.4x faster |
| `linear/4d/max10/eval10000/flipped` | 266 µs | 85.5 µs | 267 µs | 1.9 ms |  |  |  | 7.1x faster |
| `linear/4d/max40/eval10/ascending` | 327 ns | 317 ns | 333 ns | 96.6 µs | 6.01 µs |  | 271 ns | **1.2x slower** |
| `linear/4d/max40/eval10/flipped` | 308 ns | 339 ns | 306 ns | 98.4 µs |  |  |  | 3.2e+02x faster |
| `linear/4d/max40/eval100/ascending` | 3.03 µs | 3.27 µs | 3.01 µs | 121 µs | 64.8 µs |  | 2.71 µs | **1.1x slower** |
| `linear/4d/max40/eval100/flipped` | 3.21 µs | 3.03 µs | 3.21 µs | 111 µs |  |  |  | 35x faster |
| `linear/4d/max40/eval1000/ascending` | 29.4 µs | 32.2 µs | 32.2 µs | 253 µs | 629 µs |  | 27 µs | **1.1x slower** |
| `linear/4d/max40/eval1000/flipped` | 31.9 µs | 32.1 µs | 29.4 µs | 249 µs |  |  |  | 7.8x faster |
| `linear/4d/max40/eval10000/ascending` | 306 µs | 105 µs | 321 µs | 2.4 ms | 6.59 ms |  | 296 µs | **1x slower** |
| `linear/4d/max40/eval10000/flipped` | 321 µs | 101 µs | 322 µs | 2.39 ms |  |  |  | 7.4x faster |
| `linear/4d/max200/eval10/ascending` | 415 ns | 420 ns | 383 ns | 105 µs | 7.99 µs |  | 391 ns | **1.1x slower** |
| `linear/4d/max200/eval10/flipped` | 414 ns | 421 ns | 412 ns | 107 µs |  |  |  | 2.6e+02x faster |
| `linear/4d/max200/eval100/ascending` | 4.1 µs | 4.17 µs | 4.1 µs | 120 µs | 80 µs |  | 3.65 µs | **1.1x slower** |
| `linear/4d/max200/eval100/flipped` | 4.08 µs | 4.2 µs | 4.09 µs | 120 µs |  |  |  | 30x faster |
| `linear/4d/max200/eval1000/ascending` | 40.7 µs | 41.3 µs | 40.9 µs | 320 µs | 783 µs |  | 37 µs | **1.1x slower** |
| `linear/4d/max200/eval1000/flipped` | 41 µs | 41.5 µs | 40.9 µs | 329 µs |  |  |  | 8x faster |
| `linear/4d/max200/eval10000/ascending` | 410 µs | 124 µs | 379 µs | 3.2 ms | 8.32 ms |  | 400 µs | **1x slower** |
| `linear/4d/max200/eval10000/flipped` | 411 µs | 115 µs | 414 µs | 3.06 ms |  |  |  | 7.4x faster |
| `nearest/2d/each40/eval10000/ascending` | 165 µs | 64.2 µs | 151 µs | 647 µs | 1.57 ms |  | 148 µs | **1.1x slower** |
| `nearest/3d/each40/eval10000/ascending` | 212 µs | 79.2 µs | 233 µs | 1.02 ms | 1.74 ms |  | 211 µs | **1x slower** |
| `pchip/2d/each40/eval10000/ascending` | 545 µs | 136 µs | 496 µs | 1.27 s |  |  |  | 2.3e+03x faster |
| `pchip/3d/each40/eval10000/ascending` | 1.6 ms | 329 µs | 1.61 ms | 3.39 s |  |  |  | 2.1e+03x faster |
| `cubic/2d/each40/eval10000/ascending` | 332 µs | 107 µs | 331 µs | 1.41 ms |  |  | 449 µs | 1.4x faster |
| `cubic/3d/each40/eval10000/ascending` | 724 µs | 178 µs | 724 µs | 3.25 ms |  |  | 851 µs | 1.2x faster |
| `cubic/2d/each10/eval1000/ascending` | 29.1 µs | 26.8 µs | 29.1 µs | 132 µs |  |  | 39.2 µs | 1.3x faster |
| `cubic/2d/each40/eval1000/ascending` | 30.1 µs | 30.3 µs | 30.2 µs | 146 µs |  |  | 40.7 µs | 1.4x faster |
| `cubic/2d/each100/eval1000/ascending` | 35.3 µs | 35.5 µs | 32.2 µs | 182 µs |  |  | 47.1 µs | 1.3x faster |
| `cubic/2d/each200/eval1000/ascending` | 37.8 µs | 38.1 µs | 36.3 µs | 218 µs |  |  | 50.2 µs | 1.3x faster |
| `cubic/2d/each400/eval1000/ascending` | 41.6 µs | 42.2 µs | 41.6 µs | 300 µs |  |  | 54.2 µs | 1.3x faster |
| `quintic/2d/each10/eval1000/ascending` | 44.1 µs | 44.8 µs | 48.5 µs | 251 µs |  |  |  | 5.7x faster |
| `quintic/2d/each40/eval1000/ascending` | 52.6 µs | 52.9 µs | 52.7 µs | 237 µs |  |  |  | 4.5x faster |
| `linear/2d/each400/eval10000/ascending` | 237 µs | 75.8 µs | 239 µs | 970 µs | 4.48 ms | 592 µs | 227 µs | **1x slower** |
| `linear/2d/each800/eval10000/ascending` | 261 µs | 87.7 µs | 241 µs | 1.17 ms | 8.01 ms | 583 µs | 265 µs | 1x faster |
| `nearest/2d/auckland-ba32-0401/heldout/ascending` | 370 ms | 63.7 ms | 373 ms | 1.38 s | 159 s |  | 821 ms | 2.2x faster |
| `linear/2d/auckland-ba32-0401/heldout/ascending` | 378 ms | 70.2 ms | 396 ms | 1.14 s | 73.3 s | 395 ms | 825 ms | 1x faster |
| `cubic/2d/auckland-ba32-0401/heldout/ascending` | 787 ms | 133 ms | 802 ms | 34.5 s |  |  | 1.87 s | 2.4x faster |
| `quintic/2d/auckland-ba32-0401/heldout/ascending` | 1.22 s | 225 ms | 1.27 s | 35 s |  |  |  | 29x faster |
| `pchip/2d/auckland-ba32-0401/heldout/ascending` | 1.08 s | 188 ms | 1.09 s |  |  |  |  |  |

### Construction

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | ndarray-interp | interpn | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/2d/max10/eval10/ascending` | 134 ns | 123 ns | 142 ns | 6.31 µs | 146 ns | 78.1 ns | 0.78 ns | **1.7e+02x slower** |
| `linear/2d/max10/eval10/flipped` | 167 ns | 171 ns | 179 ns | 12.6 µs |  |  |  | 75x faster |
| `linear/2d/max10/eval100/ascending` | 135 ns | 134 ns | 151 ns | 5.87 µs | 133 ns | 85.1 ns | 0.78 ns | **1.7e+02x slower** |
| `linear/2d/max10/eval100/flipped` | 165 ns | 167 ns | 177 ns | 11.7 µs |  |  |  | 71x faster |
| `linear/2d/max10/eval1000/ascending` | 123 ns | 126 ns | 131 ns | 5.76 µs | 145 ns | 85.4 ns | 0.68 ns | **1.8e+02x slower** |
| `linear/2d/max10/eval1000/flipped` | 157 ns | 171 ns | 177 ns | 11.8 µs |  |  |  | 75x faster |
| `linear/2d/max10/eval10000/ascending` | 134 ns | 135 ns | 131 ns | 5.81 µs | 138 ns | 78 ns | 0.68 ns | **2e+02x slower** |
| `linear/2d/max10/eval10000/flipped` | 158 ns | 162 ns | 162 ns | 11.8 µs |  |  |  | 75x faster |
| `linear/2d/max40/eval10/ascending` | 362 ns | 361 ns | 773 ns | 6.29 µs | 203 ns | 152 ns | 0.78 ns | **4.6e+02x slower** |
| `linear/2d/max40/eval10/flipped` | 655 ns | 651 ns | 1.17 µs | 11.6 µs |  |  |  | 18x faster |
| `linear/2d/max40/eval100/ascending` | 361 ns | 364 ns | 749 ns | 5.85 µs | 218 ns | 152 ns | 0.78 ns | **4.6e+02x slower** |
| `linear/2d/max40/eval100/flipped` | 646 ns | 698 ns | 1.17 µs | 11.8 µs |  |  |  | 18x faster |
| `linear/2d/max40/eval1000/ascending` | 397 ns | 362 ns | 748 ns | 5.82 µs | 218 ns | 153 ns | 0.78 ns | **5.1e+02x slower** |
| `linear/2d/max40/eval1000/flipped` | 653 ns | 687 ns | 1.19 µs | 11.9 µs |  |  |  | 18x faster |
| `linear/2d/max40/eval10000/ascending` | 387 ns | 365 ns | 832 ns | 5.8 µs | 221 ns | 152 ns | 0.78 ns | **5e+02x slower** |
| `linear/2d/max40/eval10000/flipped` | 711 ns | 708 ns | 1.18 µs | 12.8 µs |  |  |  | 18x faster |
| `linear/2d/max200/eval10/ascending` | 7.36 µs | 6.95 µs | 17.5 µs | 5.78 µs | 3 µs | 470 ns | 0.68 ns | **1.1e+04x slower** |
| `linear/2d/max200/eval10/flipped` | 14.5 µs | 14.6 µs | 26.5 µs | 12.3 µs |  |  |  | **1.2x slower** |
| `linear/2d/max200/eval100/ascending` | 6.83 µs | 6.99 µs | 15.4 µs | 5.85 µs | 2.95 µs | 515 ns | 0.78 ns | **8.7e+03x slower** |
| `linear/2d/max200/eval100/flipped` | 15.5 µs | 14.6 µs | 26.5 µs | 12.3 µs |  |  |  | **1.3x slower** |
| `linear/2d/max200/eval1000/ascending` | 7.41 µs | 7.04 µs | 20.2 µs | 5.83 µs | 3.23 µs | 515 ns | 0.78 ns | **9.5e+03x slower** |
| `linear/2d/max200/eval1000/flipped` | 14.5 µs | 15.5 µs | 26.6 µs | 12.2 µs |  |  |  | **1.2x slower** |
| `linear/2d/max200/eval10000/ascending` | 7.49 µs | 11.9 µs | 4.88 µs | 5.81 µs | 3.36 µs | 475 ns | 0.78 ns | **9.6e+03x slower** |
| `linear/2d/max200/eval10000/flipped` | 16.1 µs | 16 µs | 20.7 µs | 13.2 µs |  |  |  | **1.2x slower** |
| `linear/3d/max10/eval10/ascending` | 196 ns | 197 ns | 217 ns | 7.62 µs | 180 ns |  | 0.98 ns | **2e+02x slower** |
| `linear/3d/max10/eval10/flipped` | 299 ns | 317 ns | 348 ns | 17.1 µs |  |  |  | 57x faster |
| `linear/3d/max10/eval100/ascending` | 187 ns | 183 ns | 207 ns | 8.22 µs | 167 ns |  | 0.98 ns | **1.9e+02x slower** |
| `linear/3d/max10/eval100/flipped` | 296 ns | 298 ns | 330 ns | 16.1 µs |  |  |  | 54x faster |
| `linear/3d/max10/eval1000/ascending` | 195 ns | 182 ns | 214 ns | 8.21 µs | 163 ns |  | 0.98 ns | **2e+02x slower** |
| `linear/3d/max10/eval1000/flipped` | 293 ns | 296 ns | 329 ns | 15.7 µs |  |  |  | 53x faster |
| `linear/3d/max10/eval10000/ascending` | 181 ns | 186 ns | 208 ns | 7.53 µs | 179 ns |  | 0.88 ns | **2.1e+02x slower** |
| `linear/3d/max10/eval10000/flipped` | 296 ns | 321 ns | 326 ns | 16.8 µs |  |  |  | 57x faster |
| `linear/3d/max40/eval10/ascending` | 2.96 µs | 2.9 µs | 7.65 µs | 8.24 µs | 1.45 µs |  | 0.88 ns | **3.4e+03x slower** |
| `linear/3d/max40/eval10/flipped` | 8.04 µs | 14.3 µs | 12.5 µs | 15.6 µs |  |  |  | 1.9x faster |
| `linear/3d/max40/eval100/ascending` | 3.22 µs | 3.16 µs | 4.13 µs | 7.48 µs | 1.58 µs |  | 0.98 ns | **3.3e+03x slower** |
| `linear/3d/max40/eval100/flipped` | 8.81 µs | 8.91 µs | 9.17 µs | 15.9 µs |  |  |  | 1.8x faster |
| `linear/3d/max40/eval1000/ascending` | 3.12 µs | 3 µs | 2.72 µs | 7.56 µs | 1.59 µs |  | 0.98 ns | **3.2e+03x slower** |
| `linear/3d/max40/eval1000/flipped` | 8.92 µs | 8.9 µs | 9.97 µs | 17.1 µs |  |  |  | 1.9x faster |
| `linear/3d/max40/eval10000/ascending` | 2.97 µs | 3.17 µs | 4.06 µs | 8.22 µs | 1.44 µs |  | 0.98 ns | **3e+03x slower** |
| `linear/3d/max40/eval10000/flipped` | 8.62 µs | 8.92 µs | 9.39 µs | 17.3 µs |  |  |  | 2x faster |
| `linear/3d/max200/eval10/ascending` | 1.82 ms | 1.81 ms | 631 µs | 8.26 µs | 1.26 ms |  | 0.98 ns | **1.9e+06x slower** |
| `linear/3d/max200/eval10/flipped` | 3.05 ms | 3.21 ms | 1.9 ms | 17.4 µs |  |  |  | **1.8e+02x slower** |
| `linear/3d/max200/eval100/ascending` | 1.86 ms | 1.71 ms | 603 µs | 8.26 µs | 1.44 ms |  | 0.88 ns | **2.1e+06x slower** |
| `linear/3d/max200/eval100/flipped` | 3.14 ms | 3.16 ms | 1.88 ms | 17.5 µs |  |  |  | **1.8e+02x slower** |
| `linear/3d/max200/eval1000/ascending` | 1.77 ms | 1.81 ms | 683 µs | 8.27 µs | 1.28 ms |  | 0.88 ns | **2e+06x slower** |
| `linear/3d/max200/eval1000/flipped` | 3.11 ms | 3.1 ms | 1.83 ms | 15.9 µs |  |  |  | **2e+02x slower** |
| `linear/3d/max200/eval10000/ascending` | 1.78 ms | 1.9 ms | 639 µs | 7.88 µs | 1.38 ms |  | 0.88 ns | **2e+06x slower** |
| `linear/3d/max200/eval10000/flipped` | 3.13 ms | 3.25 ms | 1.95 ms | 17.6 µs |  |  |  | **1.8e+02x slower** |
| `linear/4d/max10/eval10/ascending` | 235 ns | 234 ns | 247 ns | 10.6 µs |  |  |  | 45x faster |
| `linear/4d/max10/eval10/flipped` | 336 ns | 365 ns | 408 ns | 18.2 µs |  |  |  | 54x faster |
| `linear/4d/max10/eval100/ascending` | 237 ns | 235 ns | 240 ns | 10.5 µs |  |  |  | 44x faster |
| `linear/4d/max10/eval100/flipped` | 370 ns | 364 ns | 374 ns | 17.9 µs |  |  |  | 49x faster |
| `linear/4d/max10/eval1000/ascending` | 236 ns | 235 ns | 249 ns | 10.5 µs |  |  |  | 45x faster |
| `linear/4d/max10/eval1000/flipped` | 343 ns | 365 ns | 411 ns | 18.5 µs |  |  |  | 54x faster |
| `linear/4d/max10/eval10000/ascending` | 236 ns | 233 ns | 249 ns | 9.6 µs |  |  |  | 41x faster |
| `linear/4d/max10/eval10000/flipped` | 343 ns | 357 ns | 390 ns | 18.1 µs |  |  |  | 53x faster |
| `linear/4d/max40/eval10/ascending` | 13.5 µs | 13.6 µs | 36.6 µs | 9.44 µs | 5.97 µs |  | 1.27 ns | **1.1e+04x slower** |
| `linear/4d/max40/eval10/flipped` | 51.5 µs | 54.3 µs | 69.3 µs | 20 µs |  |  |  | **2.6x slower** |
| `linear/4d/max40/eval100/ascending` | 13.1 µs | 14.3 µs | 32.8 µs | 10.3 µs | 5.97 µs |  | 1.37 ns | **9.6e+03x slower** |
| `linear/4d/max40/eval100/flipped` | 53.2 µs | 55.9 µs | 69.7 µs | 20.2 µs |  |  |  | **2.6x slower** |
| `linear/4d/max40/eval1000/ascending` | 14.2 µs | 14.3 µs | 30.4 µs | 9.48 µs | 5.99 µs |  | 1.27 ns | **1.1e+04x slower** |
| `linear/4d/max40/eval1000/flipped` | 57.2 µs | 52.4 µs | 69.3 µs | 19.8 µs |  |  |  | **2.9x slower** |
| `linear/4d/max40/eval10000/ascending` | 14.3 µs | 14.8 µs | 34.2 µs | 9.41 µs | 6.02 µs |  | 1.27 ns | **1.1e+04x slower** |
| `linear/4d/max40/eval10000/flipped` | 51.9 µs | 57.3 µs | 74.8 µs | 21.7 µs |  |  |  | **2.4x slower** |
| `linear/4d/max200/eval10/ascending` | 56.4 ms | 50.4 ms | 14.9 ms | 10.2 µs | 39.2 ms |  | 1.27 ns | **4.4e+07x slower** |
| `linear/4d/max200/eval10/flipped` | 115 ms | 106 ms | 68.3 ms | 20.2 µs |  |  |  | **5.7e+03x slower** |
| `linear/4d/max200/eval100/ascending` | 51.8 ms | 48.7 ms | 16.9 ms | 10.1 µs | 39.7 ms |  | 1.27 ns | **4.1e+07x slower** |
| `linear/4d/max200/eval100/flipped` | 115 ms | 111 ms | 73.1 ms | 21.9 µs |  |  |  | **5.2e+03x slower** |
| `linear/4d/max200/eval1000/ascending` | 49.9 ms | 49.3 ms | 18.3 ms | 10.2 µs | 39.6 ms |  | 1.27 ns | **3.9e+07x slower** |
| `linear/4d/max200/eval1000/flipped` | 118 ms | 139 ms | 93.2 ms | 22.4 µs |  |  |  | **5.3e+03x slower** |
| `linear/4d/max200/eval10000/ascending` | 51.6 ms | 50.7 ms | 18.2 ms | 10.4 µs | 46.9 ms |  | 1.37 ns | **3.8e+07x slower** |
| `linear/4d/max200/eval10000/flipped` | 105 ms | 106 ms | 74.5 ms | 22 µs |  |  |  | **4.8e+03x slower** |
| `nearest/2d/each40/eval10000/ascending` | 642 ns | 635 ns | 1.53 µs | 6.03 µs | 295 ns |  | 0.78 ns | **8.2e+02x slower** |
| `nearest/3d/each40/eval10000/ascending` | 21.1 µs | 23.1 µs | 59.2 µs | 8.16 µs | 9.47 µs |  | 0.88 ns | **2.4e+04x slower** |
| `pchip/2d/each40/eval10000/ascending` | 7.05 µs | 7.12 µs | 6.29 µs | 7.2 µs |  |  |  | 1x faster |
| `pchip/3d/each40/eval10000/ascending` | 474 µs | 553 µs | 496 µs | 9.37 µs |  |  |  | **51x slower** |
| `cubic/2d/each40/eval10000/ascending` | 9.22 µs | 8.41 µs | 8.94 µs | 2.02 ms |  |  | 0.98 ns | **9.4e+03x slower** |
| `cubic/3d/each40/eval10000/ascending` | 321 µs | 280 µs | 304 µs | 158 ms |  |  | 1.17 ns | **2.7e+05x slower** |
| `cubic/2d/each10/eval1000/ascending` | 1.52 µs | 1.48 µs | 1.57 µs | 997 µs |  |  | 0.88 ns | **1.7e+03x slower** |
| `cubic/2d/each40/eval1000/ascending` | 9.19 µs | 8.8 µs | 8.18 µs | 1.97 ms |  |  | 0.88 ns | **1e+04x slower** |
| `cubic/2d/each100/eval1000/ascending` | 43.1 µs | 45 µs | 42.2 µs | 7.95 ms |  |  | 0.98 ns | **4.4e+04x slower** |
| `cubic/2d/each200/eval1000/ascending` | 190 µs | 211 µs | 207 µs | 27.2 ms |  |  | 0.98 ns | **1.9e+05x slower** |
| `cubic/2d/each400/eval1000/ascending` | 1.48 ms | 1.43 ms | 1.24 ms | 148 ms |  |  | 0.88 ns | **1.7e+06x slower** |
| `quintic/2d/each10/eval1000/ascending` | 1.94 µs | 2.15 µs | 2.21 µs | 3.66 ms |  |  |  | 1.9e+03x faster |
| `quintic/2d/each40/eval1000/ascending` | 12.5 µs | 13.5 µs | 13.1 µs | 7.74 ms |  |  |  | 6.2e+02x faster |
| `linear/2d/each400/eval10000/ascending` | 320 µs | 326 µs | 105 µs | 6.5 µs | 272 µs | 1.14 µs | 0.78 ns | **4.1e+05x slower** |
| `linear/2d/each800/eval10000/ascending` | 1.04 ms | 1.13 ms | 423 µs | 6.95 µs | 844 µs | 2.23 µs | 0.78 ns | **1.3e+06x slower** |
| `nearest/2d/auckland-ba32-0401/heldout/ascending` | 16.7 ms | 16.7 ms | 5.26 ms | 7.68 µs | 13.5 ms |  | 0.78 ns | **2.1e+07x slower** |
| `linear/2d/auckland-ba32-0401/heldout/ascending` | 16.5 ms | 17 ms | 5.76 ms | 7.59 µs | 13.6 ms | 8.24 µs | 0.68 ns | **2.4e+07x slower** |
| `cubic/2d/auckland-ba32-0401/heldout/ascending` | 83.8 ms | 81.1 ms | 76.8 ms | 15 s |  |  | 0.88 ns | **9.5e+07x slower** |
| `quintic/2d/auckland-ba32-0401/heldout/ascending` | 96.7 ms | 96.3 ms | 90.6 ms | 21.5 s |  |  |  | 2.2e+02x faster |
| `pchip/2d/auckland-ba32-0401/heldout/ascending` | 76 ms | 75.9 ms | 66.1 ms |  |  |  |  |  |

## `NearestNDInterpolator`

`n` points of the unit square, evaluated at the pixels of a square image.

### Evaluation

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: |
| `n10/image50` | 36.7 µs | 36.8 µs | 115 µs | 480 µs | 13x faster |
| `n10/image100` | 140 µs | 56 µs | 492 µs | 1.85 ms | 13x faster |
| `n10/image500` | 3.3 ms | 590 µs | 12.2 ms | 48.1 ms | 15x faster |
| `n100/image50` | 125 µs | 107 µs | 206 µs | 633 µs | 5.1x faster |
| `n100/image100` | 430 µs | 125 µs | 789 µs | 2.37 ms | 5.5x faster |
| `n100/image500` | 9.94 ms | 1.62 ms | 17.3 ms | 57.4 ms | 5.8x faster |
| `n1000/image50` | 351 µs | 124 µs | 353 µs | 797 µs | 2.3x faster |
| `n1000/image100` | 1.17 ms | 305 µs | 1.18 ms | 3.17 ms | 2.7x faster |
| `n1000/image500` | 24.7 ms | 4.46 ms | 22.5 ms | 72.5 ms | 2.9x faster |
| `n10000/image50` | 335 µs | 110 µs | 329 µs | 1.05 ms | 3.1x faster |
| `n10000/image100` | 1.36 ms | 304 µs | 1.45 ms | 4.07 ms | 3x faster |
| `n10000/image500` | 32.5 ms | 5.53 ms | 32.6 ms | 91.4 ms | 2.8x faster |
| `n100000/image50` | 487 µs | 117 µs | 449 µs | 1.46 ms | 3x faster |
| `n100000/image100` | 1.79 ms | 318 µs | 1.86 ms | 5.69 ms | 3.2x faster |
| `n100000/image500` | 43.5 ms | 7.24 ms | 45.3 ms | 175 ms | 4x faster |
| `ghcn-tmax-20250101/heldout` | 172 µs | 73.3 µs | 173 µs | 558 µs | 3.2x faster |

### Construction

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: |
| `n10/image50` | 564 ns | 569 ns | 574 ns | 10.5 µs | 19x faster |
| `n10/image100` | 557 ns | 568 ns | 568 ns | 10.7 µs | 19x faster |
| `n10/image500` | 575 ns | 570 ns | 574 ns | 10.6 µs | 18x faster |
| `n100/image50` | 3.6 µs | 3.69 µs | 3.62 µs | 16.6 µs | 4.6x faster |
| `n100/image100` | 3.28 µs | 3.6 µs | 3.6 µs | 16.5 µs | 5x faster |
| `n100/image500` | 3.61 µs | 3.35 µs | 3.58 µs | 16.6 µs | 4.6x faster |
| `n1000/image50` | 34.9 µs | 35.2 µs | 38.1 µs | 117 µs | 3.4x faster |
| `n1000/image100` | 35.5 µs | 38.3 µs | 34.7 µs | 126 µs | 3.6x faster |
| `n1000/image500` | 36.4 µs | 34.9 µs | 34.7 µs | 127 µs | 3.5x faster |
| `n10000/image50` | 532 µs | 528 µs | 486 µs | 1.63 ms | 3.1x faster |
| `n10000/image100` | 496 µs | 507 µs | 502 µs | 1.67 ms | 3.4x faster |
| `n10000/image500` | 493 µs | 543 µs | 534 µs | 1.79 ms | 3.6x faster |
| `n100000/image50` | 8.52 ms | 8.44 ms | 8.59 ms | 25.8 ms | 3x faster |
| `n100000/image100` | 8.37 ms | 8.59 ms | 8.32 ms | 24.4 ms | 2.9x faster |
| `n100000/image500` | 8.74 ms | 8.73 ms | 8.39 ms | 25.9 ms | 3x faster |
| `ghcn-tmax-20250101/heldout` | 558 µs | 560 µs | 564 µs | 1.87 ms | 3.4x faster |

## `LinearNDInterpolator`

`n` points of the unit square, evaluated at the pixels of a square image (some outside the hull, with no value); `3d`: `n` points of the unit cube, evaluated at 10 000 points of the cube in no order.

### Evaluation

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | spade | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `n10/image50` | 18.5 µs | 18.6 µs | 18.4 µs | 29.9 µs | 52 µs | 95.8 µs | 1.6x faster |
| `n10/image100` | 81.5 µs | 82.1 µs | 81.6 µs | 104 µs | 258 µs | 419 µs | 1.3x faster |
| `n10/image500` | 2.02 ms | 426 µs | 1.95 ms | 2.83 ms | 5.94 ms | 10.3 ms | 1.4x faster |
| `n100/image50` | 58 µs | 52.7 µs | 50.6 µs | 74 µs | 909 µs | 218 µs | 1.3x faster |
| `n100/image100` | 238 µs | 81.2 µs | 233 µs | 254 µs | 2.83 ms | 800 µs | 1.1x faster |
| `n100/image500` | 5.49 ms | 975 µs | 5.44 ms | 5.98 ms | 82.5 ms | 17.7 ms | 1.1x faster |
| `n1000/image50` | 52.1 µs | 56.9 µs | 55.6 µs | 133 µs | 18 ms | 298 µs | 2.5x faster |
| `n1000/image100` | 278 µs | 92.5 µs | 280 µs | 391 µs | 56.7 ms | 1.16 ms | 1.4x faster |
| `n1000/image500` | 6.5 ms | 1.2 ms | 6.32 ms | 6.6 ms | 756 ms | 20.8 ms | 1x faster |
| `n10000/image50` | 65.9 µs | 65.4 µs | 65.7 µs | 411 µs | 308 ms | 610 µs | 6.2x faster |
| `n10000/image100` | 316 µs | 97.5 µs | 316 µs | 1 ms | 1.11 s | 1.76 ms | 3.2x faster |
| `n10000/image500` | 7.75 ms | 1.45 ms | 8.11 ms | 9.21 ms | 22 s | 33.8 ms | 1.2x faster |
| `n100000/image50` | 89.7 µs | 93 µs | 90.6 µs | 1.33 ms | 4 s | 6.52 ms | 15x faster |
| `n100000/image100` | 405 µs | 131 µs | 421 µs | 3.45 ms | 15.3 s | 6.8 ms | 8.5x faster |
| `n100000/image500` | 11.7 ms | 2 ms | 11.6 ms | 21.6 ms | 308 s | 52.5 ms | 1.9x faster |
| `3d/n1000/scattered` | 821 µs | 194 µs | 788 µs | 5.61 ms |  |  | 6.8x faster |
| `3d/n10000/scattered` | 1.25 ms | 274 µs | 1.24 ms | 15.3 ms |  |  | 12x faster |
| `3d/n100000/scattered` | 2.66 ms | 503 µs | 3.22 ms | 72.3 ms |  |  | 27x faster |
| `ghcn-tmax-20250101/heldout` | 75.6 µs | 93.7 µs | 79.5 µs | 1.18 ms | 140 ms | 1.32 ms | 16x faster |

### Construction

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | spade | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `n10/image50` | 5.7 µs | 5.2 µs | 5.16 µs | 440 µs | 1.88 µs | 1.76 µs | **3.2x slower** |
| `n10/image100` | 5.34 µs | 5.69 µs | 5.61 µs | 469 µs | 2.03 µs | 1.91 µs | **2.8x slower** |
| `n10/image500` | 9.5 µs | 4.65 µs | 5.6 µs | 538 µs | 2.04 µs | 1.9 µs | **5x slower** |
| `n100/image50` | 36.8 µs | 35 µs | 36.2 µs | 674 µs | 56.2 µs | 20.7 µs | **1.8x slower** |
| `n100/image100` | 36.8 µs | 32.5 µs | 33.5 µs | 626 µs | 56.2 µs | 18.8 µs | **2x slower** |
| `n100/image500` | 36.8 µs | 32.8 µs | 36.5 µs | 738 µs | 61.6 µs | 20.9 µs | **1.8x slower** |
| `n1000/image50` | 627 µs | 615 µs | 599 µs | 2.41 ms | 4.52 ms | 275 µs | **2.3x slower** |
| `n1000/image100` | 625 µs | 647 µs | 649 µs | 2.44 ms | 4.55 ms | 286 µs | **2.2x slower** |
| `n1000/image500` | 655 µs | 660 µs | 647 µs | 2.46 ms | 4.61 ms | 299 µs | **2.2x slower** |
| `n10000/image50` | 7.95 ms | 7.76 ms | 7.86 ms | 25 ms | 497 ms | 4.06 ms | **2x slower** |
| `n10000/image100` | 8.05 ms | 7.81 ms | 7.98 ms | 25 ms | 486 ms | 4.19 ms | **1.9x slower** |
| `n10000/image500` | 7.94 ms | 8 ms | 7.89 ms | 25.4 ms | 1.14 s | 3.98 ms | **2x slower** |
| `n100000/image50` | 79.3 ms | 80.2 ms | 88 ms | 429 ms | 67 s | 47.7 ms | **1.7x slower** |
| `n100000/image100` | 81.1 ms | 79.8 ms | 78.7 ms | 367 ms | 74.5 s | 47.3 ms | **1.7x slower** |
| `n100000/image500` | 85.6 ms | 78.8 ms | 81.8 ms | 393 ms | 65.7 s | 47.3 ms | **1.8x slower** |
| `3d/n1000/scattered` | 4.71 ms | 4.67 ms | 4.63 ms | 9.61 ms |  |  | 2x faster |
| `3d/n10000/scattered` | 53.4 ms | 53 ms | 52.3 ms | 126 ms |  |  | 2.4x faster |
| `3d/n100000/scattered` | 601 ms | 568 ms | 602 ms | 1.9 s |  |  | 3.2x faster |
| `ghcn-tmax-20250101/heldout` | 7.94 ms | 8.32 ms | 7.86 ms | 30 ms | 520 ms | 4.9 ms | **1.6x slower** |

## `RBFInterpolator`

`n` points of a square, evaluated at 1000 points in no order (`scattered`) or at the 10 000 pixels of an image (`image`). `global`: all points; `k50`: the 50 nearest.

### Evaluation

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/global/n10/scattered` | 15.7 µs | 14.5 µs | 15.7 µs | 64.3 µs | 199 µs | 4.1x faster |
| `linear/k50/n10/scattered` | 15.7 µs | 14.9 µs | 15.8 µs | 1.47 ms |  | 93x faster |
| `thin_plate_spline/global/n10/scattered` | 79.5 µs | 79.4 µs | 79.3 µs | 114 µs | 326 µs | 1.4x faster |
| `thin_plate_spline/k50/n10/scattered` | 73.1 µs | 79.2 µs | 79.2 µs | 1.52 ms |  | 21x faster |
| `cubic/global/n10/scattered` | 18.8 µs | 18 µs | 17.4 µs | 63.8 µs | 217 µs | 3.4x faster |
| `cubic/k50/n10/scattered` | 18.8 µs | 18 µs | 17.3 µs | 1.49 ms |  | 79x faster |
| `quintic/global/n10/scattered` | 27.4 µs | 26.6 µs | 26.5 µs | 69.9 µs | 234 µs | 2.5x faster |
| `quintic/k50/n10/scattered` | 26 µs | 26.7 µs | 26.4 µs | 1.47 ms |  | 56x faster |
| `multiquadric/global/n10/scattered` | 30.5 µs | 30.7 µs | 30.5 µs | 54.9 µs | 257 µs | 1.8x faster |
| `multiquadric/k50/n10/scattered` | 30.5 µs | 27.9 µs | 30.5 µs | 1.48 ms |  | 48x faster |
| `inverse_multiquadric/global/n10/scattered` | 41.1 µs | 41.2 µs | 41.1 µs | 61.1 µs | 267 µs | 1.5x faster |
| `inverse_multiquadric/k50/n10/scattered` | 41.1 µs | 41.3 µs | 41.1 µs | 1.43 ms |  | 35x faster |
| `inverse_quadratic/global/n10/scattered` | 23.8 µs | 26.3 µs | 26.2 µs | 47.3 µs |  | 2x faster |
| `inverse_quadratic/k50/n10/scattered` | 25.7 µs | 24.5 µs | 23.9 µs | 1.35 ms |  | 52x faster |
| `gaussian/global/n10/scattered` | 74.1 µs | 74.3 µs | 68 µs | 88.4 µs | 388 µs | 1.2x faster |
| `gaussian/k50/n10/scattered` | 74 µs | 74.2 µs | 74 µs | 1.45 ms |  | 20x faster |
| `linear/global/n100/scattered` | 148 µs | 48.3 µs | 135 µs | 291 µs | 1.73 ms | 2x faster |
| `linear/k50/n100/scattered` | 13.5 ms | 2.75 ms | 13.1 ms | 19.5 ms |  | 1.4x faster |
| `thin_plate_spline/global/n100/scattered` | 718 µs | 196 µs | 716 µs | 947 µs | 3.05 ms | 1.3x faster |
| `thin_plate_spline/k50/n100/scattered` | 24.5 ms | 4.45 ms | 22.9 ms | 25.7 ms |  | 1x faster |
| `cubic/global/n100/scattered` | 152 µs | 60.3 µs | 140 µs | 321 µs | 1.93 ms | 2.1x faster |
| `cubic/k50/n100/scattered` | 16.2 ms | 3.07 ms | 16 ms | 21.5 ms |  | 1.3x faster |
| `quintic/global/n100/scattered` | 202 µs | 66.8 µs | 168 µs | 329 µs | 2.24 ms | 1.6x faster |
| `quintic/k50/n100/scattered` | 17.5 ms | 3.35 ms | 16.1 ms | 23.6 ms |  | 1.3x faster |
| `multiquadric/global/n100/scattered` | 295 µs | 90 µs | 295 µs | 303 µs | 2.46 ms | 1x faster |
| `multiquadric/k50/n100/scattered` | 15.3 ms | 3.16 ms | 15.8 ms | 20.8 ms |  | 1.4x faster |
| `inverse_multiquadric/global/n100/scattered` | 399 µs | 112 µs | 366 µs | 434 µs | 2.63 ms | 1.1x faster |
| `inverse_multiquadric/k50/n100/scattered` | 16.6 ms | 3.12 ms | 16.4 ms | 21.5 ms |  | 1.3x faster |
| `inverse_quadratic/global/n100/scattered` | 246 µs | 78.5 µs | 246 µs | 291 µs |  | 1.2x faster |
| `inverse_quadratic/k50/n100/scattered` | 14.8 ms | 3.23 ms | 14.4 ms | 20.2 ms |  | 1.4x faster |
| `gaussian/global/n100/scattered` | 676 µs | 173 µs | 676 µs | 749 µs | 3.82 ms | 1.1x faster |
| `gaussian/k50/n100/scattered` | 21.9 ms | 3.91 ms | 21.7 ms | 23 ms |  | 1x faster |
| `linear/global/n1000/scattered` | 1.47 ms | 329 µs | 1.47 ms | 4 ms | 18.2 ms | 2.7x faster |
| `linear/k50/n1000/scattered` | 13.7 ms | 2.84 ms | 14.5 ms | 36.2 ms |  | 2.6x faster |
| `thin_plate_spline/global/n1000/scattered` | 7.16 ms | 1.5 ms | 7.73 ms | 10.5 ms | 32.8 ms | 1.5x faster |
| `thin_plate_spline/k50/n1000/scattered` | 23.5 ms | 4.35 ms | 23.2 ms | 46.1 ms |  | 2x faster |
| `cubic/global/n1000/scattered` | 1.48 ms | 354 µs | 1.49 ms | 3.94 ms | 21.1 ms | 2.7x faster |
| `cubic/k50/n1000/scattered` | 15.9 ms | 3.12 ms | 16 ms | 42.6 ms |  | 2.7x faster |
| `quintic/global/n1000/scattered` | 1.51 ms | 370 µs | 1.54 ms | 3.98 ms | 22.2 ms | 2.6x faster |
| `quintic/k50/n1000/scattered` | 18.1 ms | 3.54 ms | 17.2 ms | 44.2 ms |  | 2.4x faster |
| `multiquadric/global/n1000/scattered` | 2.8 ms | 687 µs | 2.81 ms | 4.06 ms | 25.1 ms | 1.4x faster |
| `multiquadric/k50/n1000/scattered` | 16.4 ms | 3.18 ms | 16.4 ms | 38.4 ms |  | 2.3x faster |
| `inverse_multiquadric/global/n1000/scattered` | 4 ms | 834 µs | 3.96 ms | 5.41 ms | 29.1 ms | 1.4x faster |
| `inverse_multiquadric/k50/n1000/scattered` | 18.6 ms | 3.4 ms | 19.1 ms | 40.1 ms |  | 2.2x faster |
| `inverse_quadratic/global/n1000/scattered` | 2.45 ms | 564 µs | 2.46 ms | 4.01 ms |  | 1.6x faster |
| `inverse_quadratic/k50/n1000/scattered` | 17.1 ms | 3.06 ms | 17.1 ms | 37.6 ms |  | 2.2x faster |
| `gaussian/global/n1000/scattered` | 6.57 ms | 1.37 ms | 6.48 ms | 8.83 ms | 37.6 ms | 1.3x faster |
| `gaussian/k50/n1000/scattered` | 23.4 ms | 4.15 ms | 22.8 ms | 46.4 ms |  | 2x faster |
| `thin_plate_spline/global/n3000/scattered` | 23.4 ms | 4.2 ms | 23.2 ms | 32.5 ms | 99.7 ms | 1.4x faster |
| `thin_plate_spline/k50/n10000/scattered` | 24.8 ms | 4.65 ms | 24.5 ms | 50.5 ms |  | 2x faster |
| `thin_plate_spline/k50/n10000/image` | 253 ms | 39.8 ms | 258 ms | 487 ms |  | 1.9x faster |
| `thin_plate_spline/k50/n100000/scattered` | 24.8 ms | 5.1 ms | 27.6 ms | 51.2 ms |  | 2.1x faster |
| `thin_plate_spline/k50/n100000/image` | 244 ms | 39.6 ms | 251 ms | 499 ms |  | 2x faster |
| `linear/global/ghcn-tmax-20250101/heldout` | 19.5 ms | 3.7 ms | 18.5 ms | 55.4 ms | 236 ms | 2.8x faster |
| `linear/k50/ghcn-tmax-20250101/heldout` | 19.4 ms | 3.51 ms | 19.5 ms | 44.4 ms |  | 2.3x faster |
| `thin_plate_spline/global/ghcn-tmax-20250101/heldout` | 83 ms | 16.1 ms | 83.2 ms | 113 ms | 402 ms | 1.4x faster |
| `thin_plate_spline/k50/ghcn-tmax-20250101/heldout` | 30 ms | 5.63 ms | 30.4 ms | 58.5 ms |  | 1.9x faster |
| `cubic/global/ghcn-tmax-20250101/heldout` | 18.5 ms | 3.92 ms | 18 ms | 53.2 ms | 268 ms | 2.9x faster |
| `cubic/k50/ghcn-tmax-20250101/heldout` | 21.5 ms | 4.11 ms | 21 ms | 49.8 ms |  | 2.3x faster |
| `quintic/global/ghcn-tmax-20250101/heldout` | 21 ms | 4.13 ms | 21.1 ms | 53.8 ms | 294 ms | 2.6x faster |
| `quintic/k50/ghcn-tmax-20250101/heldout` | 22.6 ms | 4.13 ms | 23.2 ms | 53.6 ms |  | 2.4x faster |

### Construction

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/global/n10/scattered` | 2.85 µs | 2.83 µs | 2.85 µs | 12.4 µs | 5.55 µs | 2x faster |
| `linear/k50/n10/scattered` | 2.87 µs | 2.62 µs | 2.87 µs | 16.4 µs |  | 5.7x faster |
| `thin_plate_spline/global/n10/scattered` | 3.65 µs | 3.64 µs | 3.61 µs | 13.7 µs | 7.15 µs | 2x faster |
| `thin_plate_spline/k50/n10/scattered` | 3.66 µs | 3.35 µs | 3.64 µs | 17.4 µs |  | 4.8x faster |
| `cubic/global/n10/scattered` | 3.8 µs | 3.54 µs | 3.72 µs | 14 µs | 5.49 µs | 1.4x faster |
| `cubic/k50/n10/scattered` | 3.81 µs | 3.83 µs | 3.81 µs | 17.2 µs |  | 4.5x faster |
| `quintic/global/n10/scattered` | 4.25 µs | 4.3 µs | 4.27 µs | 15 µs | 6.24 µs | 1.5x faster |
| `quintic/k50/n10/scattered` | 3.9 µs | 4.05 µs | 4.28 µs | 18.8 µs |  | 4.8x faster |
| `multiquadric/global/n10/scattered` | 2.93 µs | 2.94 µs | 2.91 µs | 12.5 µs | 7.18 µs | 2.5x faster |
| `multiquadric/k50/n10/scattered` | 2.95 µs | 2.68 µs | 2.71 µs | 16.4 µs |  | 5.6x faster |
| `inverse_multiquadric/global/n10/scattered` | 2.91 µs | 2.9 µs | 2.91 µs | 12.7 µs | 7.13 µs | 2.4x faster |
| `inverse_multiquadric/k50/n10/scattered` | 2.92 µs | 2.92 µs | 2.9 µs | 15 µs |  | 5.1x faster |
| `inverse_quadratic/global/n10/scattered` | 2.58 µs | 2.84 µs | 2.84 µs | 11.5 µs |  | 4.5x faster |
| `inverse_quadratic/k50/n10/scattered` | 2.6 µs | 2.83 µs | 2.81 µs | 15 µs |  | 5.8x faster |
| `gaussian/global/n10/scattered` | 2.83 µs | 3.05 µs | 2.8 µs | 11.8 µs | 8.47 µs | 3x faster |
| `gaussian/k50/n10/scattered` | 2.79 µs | 3.08 µs | 2.79 µs | 15.1 µs |  | 5.4x faster |
| `linear/global/n100/scattered` | 74 µs | 74.4 µs | 70.3 µs | 83 µs | 730 µs | 1.1x faster |
| `linear/k50/n100/scattered` | 3.26 µs | 3.19 µs | 3.23 µs | 22.6 µs |  | 6.9x faster |
| `thin_plate_spline/global/n100/scattered` | 98.7 µs | 108 µs | 106 µs | 123 µs | 811 µs | 1.2x faster |
| `thin_plate_spline/k50/n100/scattered` | 3.44 µs | 3.48 µs | 3.31 µs | 22.7 µs |  | 6.6x faster |
| `cubic/global/n100/scattered` | 71.7 µs | 75 µs | 78 µs | 87.8 µs | 696 µs | 1.2x faster |
| `cubic/k50/n100/scattered` | 3.4 µs | 3.2 µs | 3.51 µs | 23.6 µs |  | 7x faster |
| `quintic/global/n100/scattered` | 84.6 µs | 82.1 µs | 83.9 µs | 101 µs | 810 µs | 1.2x faster |
| `quintic/k50/n100/scattered` | 3.53 µs | 3.42 µs | 3.23 µs | 25.2 µs |  | 7.1x faster |
| `multiquadric/global/n100/scattered` | 82.2 µs | 83.2 µs | 83.8 µs | 77.7 µs | 678 µs | **1.1x slower** |
| `multiquadric/k50/n100/scattered` | 3.2 µs | 3.44 µs | 3.39 µs | 22.6 µs |  | 7.1x faster |
| `inverse_multiquadric/global/n100/scattered` | 86.4 µs | 87.3 µs | 87.2 µs | 88 µs | 869 µs | 1x faster |
| `inverse_multiquadric/k50/n100/scattered` | 3.47 µs | 3.45 µs | 3.44 µs | 22.6 µs |  | 6.5x faster |
| `inverse_quadratic/global/n100/scattered` | 78.6 µs | 78.8 µs | 79.3 µs | 81.5 µs |  | 1x faster |
| `inverse_quadratic/k50/n100/scattered` | 3.44 µs | 3.46 µs | 3.64 µs | 22.7 µs |  | 6.6x faster |
| `gaussian/global/n100/scattered` | 105 µs | 105 µs | 110 µs | 106 µs | 798 µs | 1x faster |
| `gaussian/k50/n100/scattered` | 3.46 µs | 3.51 µs | 3.5 µs | 22.9 µs |  | 6.6x faster |
| `linear/global/n1000/scattered` | 19.8 ms | 19.6 ms | 19.4 ms | 62.3 ms | 194 ms | 3.2x faster |
| `linear/k50/n1000/scattered` | 38.7 µs | 36.1 µs | 38.4 µs | 137 µs |  | 3.5x faster |
| `thin_plate_spline/global/n1000/scattered` | 23 ms | 23.2 ms | 22.6 ms | 64.4 ms | 208 ms | 2.8x faster |
| `thin_plate_spline/k50/n1000/scattered` | 38.6 µs | 37.3 µs | 35 µs | 135 µs |  | 3.5x faster |
| `cubic/global/n1000/scattered` | 20.6 ms | 19.7 ms | 20.9 ms | 62.2 ms | 201 ms | 3x faster |
| `cubic/k50/n1000/scattered` | 35.3 µs | 38.9 µs | 38.2 µs | 131 µs |  | 3.7x faster |
| `quintic/global/n1000/scattered` | 20.5 ms | 20.3 ms | 20.7 ms | 61.3 ms | 199 ms | 3x faster |
| `quintic/k50/n1000/scattered` | 38.8 µs | 39.2 µs | 39.4 µs | 132 µs |  | 3.4x faster |
| `multiquadric/global/n1000/scattered` | 20.6 ms | 21.1 ms | 21.7 ms | 60.9 ms | 202 ms | 3x faster |
| `multiquadric/k50/n1000/scattered` | 35.1 µs | 38.7 µs | 38.2 µs | 139 µs |  | 4x faster |
| `inverse_multiquadric/global/n1000/scattered` | 21.7 ms | 21.6 ms | 21.7 ms | 60.8 ms | 219 ms | 2.8x faster |
| `inverse_multiquadric/k50/n1000/scattered` | 38.6 µs | 39.4 µs | 35.1 µs | 135 µs |  | 3.5x faster |
| `inverse_quadratic/global/n1000/scattered` | 20.7 ms | 21 ms | 21.4 ms | 60.4 ms |  | 2.9x faster |
| `inverse_quadratic/k50/n1000/scattered` | 39.1 µs | 36.4 µs | 38.8 µs | 132 µs |  | 3.4x faster |
| `gaussian/global/n1000/scattered` | 23.7 ms | 22.9 ms | 23.5 ms | 63.1 ms | 221 ms | 2.7x faster |
| `gaussian/k50/n1000/scattered` | 39 µs | 36.7 µs | 38.2 µs | 136 µs |  | 3.5x faster |
| `thin_plate_spline/global/n3000/scattered` | 491 ms | 448 ms | 477 ms | 1.91 s | 7.54 s | 3.9x faster |
| `thin_plate_spline/k50/n10000/scattered` | 793 µs | 764 µs | 658 µs | 1.81 ms |  | 2.3x faster |
| `thin_plate_spline/k50/n10000/image` | 811 µs | 768 µs | 743 µs | 1.8 ms |  | 2.2x faster |
| `thin_plate_spline/k50/n100000/scattered` | 9.07 ms | 9.84 ms | 9.45 ms | 28 ms |  | 3.1x faster |
| `thin_plate_spline/k50/n100000/image` | 9.36 ms | 9.5 ms | 9.02 ms | 26.4 ms |  | 2.8x faster |
| `linear/global/ghcn-tmax-20250101/heldout` | 20.6 s | 17.5 s | 18 s | 140 s | 294 s | 6.8x faster |
| `linear/k50/ghcn-tmax-20250101/heldout` | 774 µs | 828 µs | 724 µs | 1.73 ms |  | 2.2x faster |
| `thin_plate_spline/global/ghcn-tmax-20250101/heldout` | 18.4 s | 17.8 s | 19 s | 117 s | 300 s | 6.4x faster |
| `thin_plate_spline/k50/ghcn-tmax-20250101/heldout` | 725 µs | 754 µs | 649 µs | 1.84 ms |  | 2.5x faster |
| `cubic/global/ghcn-tmax-20250101/heldout` | 18 s | 17.9 s | 17.8 s | 124 s | 307 s | 6.9x faster |
| `cubic/k50/ghcn-tmax-20250101/heldout` | 752 µs | 821 µs | 719 µs | 1.86 ms |  | 2.5x faster |
| `quintic/global/ghcn-tmax-20250101/heldout` | 18.9 s | 17.8 s | 18.6 s | 122 s | 284 s | 6.5x faster |
| `quintic/k50/ghcn-tmax-20250101/heldout` | 828 µs | 820 µs | 725 µs | 1.87 ms |  | 2.3x faster |

## Memory

Megabytes. "Data read" is the process with the data of the case in memory, before anything is built. The columns after it are what was added to that: the highest while building, with the interpolator built, the highest while evaluating. "Highest" is the highest of the whole run, repeats included.

| case | library | data read | while building | built | while evaluating | highest |
| :--- | :--- | ---: | ---: | ---: | ---: | ---: |
| `oned/cubic_spline/n10000/sorted` | interpx | 4.4 | +0.8 | +0.9 | +0.9 | 5.9 |
| `oned/cubic_spline/n10000/sorted` | interpx, parallel | 3.2 | +2.1 | +2.1 | +2.2 | 6.5 |
| `oned/cubic_spline/n10000/sorted` | SciPy | 71.1 | +0.4 | +1.9 | +1.9 | 73.3 |
| `oned/cubic_spline/n10000/sorted` | scirs2-interpolate | 3.2 | +1.9 | +2.1 | +2.7 | 5.9 |
| `oned/cubic_spline/n10000/sorted` | ndarray-interp | 3.0 | +2.0 | +1.8 | +1.8 | 5.5 |
| `grid/linear/4d/max200/eval10000/ascending` | interpx | 225.7 | +179.0 | +179.0 | +179.0 | 404.9 |
| `grid/linear/4d/max200/eval10000/ascending` | interpx, parallel | 244.2 | +160.9 | +160.9 | +160.9 | 405.4 |
| `grid/linear/4d/max200/eval10000/ascending` | interpx, `from_owned` | 216.6 | +188.5 | +188.5 | +188.5 | 405.1 |
| `grid/linear/4d/max200/eval10000/ascending` | SciPy | 271.3 | +0.3 | +0.4 | +2.0 | 273.6 |
| `grid/linear/4d/max200/eval10000/ascending` | scirs2-interpolate | 210.1 | +193.5 | +194.8 | +194.8 | 404.9 |
| `grid/linear/4d/max200/eval10000/ascending` | interpn | 205.6 | +0.0 | +0.0 | +0.0 | 401.0 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpx | 4.4 | +4.0 | +3.9 | +3.9 | 8.5 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpx, parallel | 5.5 | +2.4 | +2.8 | +2.8 | 8.9 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpx, `from_owned` | 6.3 | +0.0 | +2.0 | +2.0 | 8.5 |
| `grid/cubic/2d/each400/eval1000/ascending` | SciPy | 72.4 | +85.2 | +5.8 | +5.8 | 161.8 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpn | 5.4 | +0.1 | +0.1 | +0.1 | 5.6 |
| `grid/linear/2d/each800/eval10000/ascending` | interpx | 12.4 | +0.7 | +2.4 | +2.4 | 14.9 |
| `grid/linear/2d/each800/eval10000/ascending` | interpx, parallel | 9.5 | +3.5 | +5.5 | +5.5 | 15.5 |
| `grid/linear/2d/each800/eval10000/ascending` | interpx, `from_owned` | 10.4 | +0.0 | +4.5 | +4.5 | 15.0 |
| `grid/linear/2d/each800/eval10000/ascending` | SciPy | 76.4 | +0.0 | +0.1 | +0.5 | 77.1 |
| `grid/linear/2d/each800/eval10000/ascending` | scirs2-interpolate | 11.6 | +0.0 | +3.3 | +3.3 | 15.0 |
| `grid/linear/2d/each800/eval10000/ascending` | ndarray-interp | 10.9 | +0.0 | +0.0 | +0.0 | 10.9 |
| `grid/linear/2d/each800/eval10000/ascending` | interpn | 11.3 | +0.0 | +0.0 | +0.0 | 11.3 |
| `nearest_nd/n100000/image500` | interpx | 12.7 | +2.7 | +1.5 | +3.5 | 18.2 |
| `nearest_nd/n100000/image500` | interpx, parallel | 14.4 | +1.1 | +0.0 | +2.4 | 18.4 |
| `nearest_nd/n100000/image500` | interpx, `from_owned` | 11.8 | +4.0 | +2.4 | +4.3 | 18.2 |
| `nearest_nd/n100000/image500` | SciPy | 77.4 | +4.4 | +4.5 | +13.2 | 92.2 |
| `linear_nd/n100000/image500` | interpx | 13.7 | +23.7 | +21.4 | +23.3 | 41.5 |
| `linear_nd/n100000/image500` | interpx, parallel | 13.5 | +22.8 | +21.5 | +24.0 | 42.9 |
| `linear_nd/n100000/image500` | interpx, `from_owned` | 15.2 | +23.0 | +20.3 | +22.3 | 42.4 |
| `linear_nd/n100000/image500` | SciPy | 77.6 | +66.7 | +13.3 | +25.2 | 145.7 |
| `linear_nd/n100000/image500` | scirs2-interpolate | 13.0 | +13.8 | +8.9 | +11.1 | 26.8 |
| `linear_nd/n100000/image500` | spade | 13.9 | +14.0 | +11.1 | +13.1 | 29.5 |
| `rbf/thin_plate_spline/global/n3000/scattered` | interpx | 3.2 | +77.2 | +4.9 | +4.9 | 80.4 |
| `rbf/thin_plate_spline/global/n3000/scattered` | interpx, parallel | 3.0 | +77.4 | +5.2 | +5.6 | 80.9 |
| `rbf/thin_plate_spline/global/n3000/scattered` | interpx, `from_owned` | 3.3 | +77.2 | +4.8 | +4.8 | 80.5 |
| `rbf/thin_plate_spline/global/n3000/scattered` | SciPy | 71.1 | +73.1 | +0.9 | +9.0 | 144.3 |
| `rbf/thin_plate_spline/global/n3000/scattered` | scirs2-interpolate | 2.9 | +289.7 | +1.6 | +1.6 | 292.5 |
| `rbf/thin_plate_spline/k50/n100000/image` | interpx | 6.7 | +3.5 | +4.7 | +5.0 | 11.9 |
| `rbf/thin_plate_spline/k50/n100000/image` | interpx, parallel | 5.9 | +4.3 | +5.5 | +6.6 | 12.7 |
| `rbf/thin_plate_spline/k50/n100000/image` | interpx, `from_owned` | 6.3 | +6.5 | +5.0 | +5.4 | 13.5 |
| `rbf/thin_plate_spline/k50/n100000/image` | SciPy | 73.8 | +5.2 | +5.3 | +25.7 | 101.3 |

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
| `LinearNDInterpolator` | 6e-15 | 6e-15 | 6e-15 | inf |  |  | 6e-15 |
| `RBFInterpolator` | 4e+00 | 4e+00 | 4e+00 | 7e+11 |  |  |  |
