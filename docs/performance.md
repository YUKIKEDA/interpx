# Performance: interpx against SciPy and other Rust crates

<!-- Written by compare/perf/report.py from a run of compare/perf/run.py. Do not edit by hand: run it again. -->

Measured on 2026-10-13, at commit `526b2f9bf4`.

## How it was measured

- Machine: Intel64 Family 6 Model 191 Stepping 2, GenuineIntel, 16 logical cores, 48 GB. Windows-11-10.0.26300-SP0.
- rustc 1.99.0 (b940084d7 2026-09-28); Python 3.13.1, NumPy 2.5.4, SciPy 1.18.1.
- The other Rust crates, at the versions in `compare/perf/runners/*/Cargo.toml`: `scirs2-interpolate` (the widest of the Rust crates: one variable, grids, RBF), `ndarray-interp` (the most used for one variable; also bilinear), `interpn` (grids).
- Every library ran each case in a process of its own, pinned to one core, on one thread, on the same data. The libraries took turns case by case, 3 rounds.
- `interpx, from_owned` is interpx built from data that is handed over (`from_owned`) and so not copied. Its column is there for the types that keep their data as it is given (the grid, `NearestNDInterpolator`, `RBFInterpolator`), next to the construction from borrowed data (`new`), which is the one the comparison counts. Each of its samples is given data made just before it, which is not in the nearest cache as the data of the other constructions is: on small grids it reads slower than `new` for that reason.
- `interpx, parallel` is interpx evaluated with `eval_batch_into_par` (the `rayon` feature) on 8 threads, not pinned to a core. Every other column is one thread on one core, and the comparison counts those. A parallel batch times its first points and is shared out only when the rest is worth it, so on small cases the column reads as the serial one does.
- A time is that of one call: construction, or the evaluation of all the queries of the case (`eval_batch_into` for interpx, the call with an array for the others). It is the median of the repeats of a process, and then the median of the rounds.
- interpx is called slower where its median is more than 1.10 times another library's and its fastest repeat is slower than that library's fastest.
- Memory is the resident set size of the process, sampled every millisecond from outside it.
- Command: `perf/run.py --rounds 3`.

## Where interpx is slower

| family | operation | cases | interpx fastest | slower than SciPy | slower than a Rust crate |
| :--- | :--- | ---: | ---: | ---: | ---: |
| One variable | evaluation | 112 | 111 | 0 | 1 |
| One variable | construction | 112 | 77 | 8 | 31 |
| `PPoly` | evaluation | 4 | 4 | 0 | 0 |
| `PPoly` | construction | 4 | 4 | 0 | 0 |
| `BSpline` | evaluation | 4 | 4 | 0 | 0 |
| `BSpline` | construction | 4 | 4 | 0 | 0 |
| `BSplineCurve` | evaluation | 2 | 2 | 0 | 0 |
| `BSplineCurve` | construction | 2 | 2 | 0 | 0 |
| Regular grid | evaluation | 87 | 68 | 0 | 19 |
| Regular grid | construction | 87 | 26 | 37 | 43 |
| `NearestNDInterpolator` | evaluation | 15 | 15 | 0 | 0 |
| `NearestNDInterpolator` | construction | 15 | 15 | 0 | 0 |
| `RBFInterpolator` | evaluation | 53 | 53 | 0 | 0 |
| `RBFInterpolator` | construction | 53 | 53 | 0 | 0 |

192 comparisons (a case, an operation, and a library), by method. "Same values" says whether that library's values are SciPy's up to rounding: where they are not, it computes another interpolant under the same name, and its time is not that of the same work.

| method | operation | against | same values | cases | slower by | the largest |
| :--- | :--- | :--- | :--- | ---: | ---: | :--- |
| `grid/linear` | construction | interpn | yes | 34 | 179x to 4.27e+07x | `grid/linear/4d/max200/eval1000/ascending`: 50 ms against 1.17 ns |
| `grid/cubic` | construction | interpn | no | 7 | 1.65e+03x to 1.82e+06x | `grid/cubic/2d/each400/eval1000/ascending`: 1.42 ms against 0 |
| `grid/nearest` | construction | interpn | yes | 2 | 750x to 2.32e+04x | `grid/nearest/3d/each40/eval10000/ascending`: 22.6 µs against 0 |
| `grid/linear` | construction | SciPy | yes | 35 | 1.16x to 5.29e+03x | `grid/linear/4d/max200/eval1000/ascending`: 50 ms against 9.47 µs |
| `grid/linear` | construction | ndarray-interp | yes | 14 | 1.57x to 511x | `grid/linear/2d/each800/eval10000/ascending`: 1.14 ms against 2.22 µs |
| `grid/pchip` | construction | SciPy | yes | 1 | 59.5x | `grid/pchip/3d/each40/eval10000/ascending`: 509 µs against 8.55 µs |
| `oned/linear` | construction | ndarray-interp | yes | 14 | 1.37x to 6.63x | `oned/linear/n10000/sorted`: 99.4 µs against 15 µs |
| `grid/nearest` | construction | SciPy | yes | 1 | 2.95x | `grid/nearest/3d/each40/eval10000/ascending`: 22.6 µs against 7.65 µs |
| `oned/pchip` | construction | scirs2-interpolate | yes | 8 | 1.56x to 2.78x | `oned/pchip/n1000/shuffled`: 14.6 µs against 5.24 µs |
| `oned/bspline1` | construction | SciPy | yes | 2 | 2.61x to 2.69x | `oned/bspline1/n10000/sorted`: 67.3 µs against 25 µs |
| `grid/linear` | construction | scirs2-interpolate | yes | 24 | 1.13x to 2.52x | `grid/linear/2d/max200/eval10/ascending`: 7.47 µs against 2.96 µs |
| `oned/cubic_spline` | construction | ndarray-interp | yes | 1 | 2.47x | `oned/cubic_spline/n1000/shuffled`: 42.7 µs against 17.3 µs |
| `grid/nearest` | construction | scirs2-interpolate | yes | 2 | 2.07x to 2.41x | `grid/nearest/3d/each40/eval10000/ascending`: 22.6 µs against 9.38 µs |
| `oned/cubic_spline` | construction | scirs2-interpolate | no | 3 | 1.12x to 2.31x | `oned/cubic_spline/n1000/shuffled`: 42.7 µs against 18.5 µs |
| `oned/linear` | construction | SciPy | yes | 4 | 1.89x to 2.03x | `oned/linear/n10000/shuffled-uneven`: 97.4 µs against 47.9 µs |
| `oned/linear` | construction | scirs2-interpolate | yes | 12 | 1.15x to 1.67x | `oned/linear/n1000/sorted`: 2.1 µs against 1.26 µs |
| `grid/linear` | evaluation | interpn | yes | 18 | 1.12x to 1.5x | `grid/linear/3d/max10/eval10/ascending`: 215 ns against 143 ns |
| `oned/akima` | construction | scirs2-interpolate | yes | 6 | 1.15x to 1.39x | `oned/akima/n100/shuffled`: 743 ns against 534 ns |
| `oned/nearest` | construction | SciPy | yes | 2 | 1.19x to 1.22x | `oned/nearest/n10000/sorted`: 66.6 µs against 54.6 µs |
| `oned/nearest` | evaluation | scirs2-interpolate | yes | 1 | 1.12x | `oned/nearest/n10/shuffled`: 299 ns against 266 ns |
| `grid/nearest` | evaluation | interpn | yes | 1 | 1.11x | `grid/nearest/2d/each40/eval10000/ascending`: 150 µs against 136 µs |

<details><summary>Every comparison where interpx is slower</summary>

"Values apart" is the distance of that library's values from SciPy's, as a share of the largest value.

| case | operation | against | interpx | the other | slower by | values apart |
| :--- | :--- | :--- | ---: | ---: | ---: | ---: |
| `grid/linear/4d/max200/eval1000/ascending` | construction | interpn | 50 ms | 1.17 ns | 42694314.67x | 3e-16 |
| `grid/linear/4d/max200/eval10/ascending` | construction | interpn | 46.3 ms | 1.17 ns | 39525376.00x | 2e-16 |
| `grid/linear/4d/max200/eval100/ascending` | construction | interpn | 48.3 ms | 1.27 ns | 38038055.38x | 3e-16 |
| `grid/linear/4d/max200/eval10000/ascending` | construction | interpn | 48.2 ms | 1.27 ns | 37992999.38x | 4e-16 |
| `grid/linear/3d/max200/eval10000/ascending` | construction | interpn | 1.75 ms | 0 | 1987356.44x | 5e-16 |
| `grid/linear/3d/max200/eval10/ascending` | construction | interpn | 1.73 ms | 0 | 1967388.44x | 3e-16 |
| `grid/linear/3d/max200/eval100/ascending` | construction | interpn | 1.71 ms | 0 | 1943324.44x | 3e-16 |
| `grid/linear/3d/max200/eval1000/ascending` | construction | interpn | 1.69 ms | 0 | 1928192.00x | 3e-16 |
| `grid/cubic/2d/each400/eval1000/ascending` | construction | interpn | 1.42 ms | 0 | 1822592.00x | 5e-01 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | interpn | 1.14 ms | 0 | 1663414.86x | 3e-16 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | interpn | 312 µs | 0 | 456118.86x | 3e-16 |
| `grid/cubic/3d/each40/eval10000/ascending` | construction | interpn | 423 µs | 0 | 433459.20x | 5e-01 |
| `grid/cubic/2d/each200/eval1000/ascending` | construction | interpn | 179 µs | 0 | 183705.60x | 2e-01 |
| `grid/cubic/2d/each100/eval1000/ascending` | construction | interpn | 84 µs | 0 | 85977.60x | 2e-01 |
| `grid/linear/2d/max200/eval10000/ascending` | construction | interpn | 47.2 µs | 0 | 69046.86x | 3e-16 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | interpn | 22.6 µs | 0 | 23161.60x | 0e+00 |
| `grid/linear/3d/max40/eval10000/ascending` | construction | interpn | 16.3 µs | 0 | 18524.44x | 3e-16 |
| `grid/linear/3d/max40/eval1000/ascending` | construction | interpn | 15.2 µs | 0 | 17251.56x | 2e-16 |
| `grid/linear/4d/max40/eval1000/ascending` | construction | interpn | 14.4 µs | 1.27 ns | 11328.00x | 6e-16 |
| `grid/linear/4d/max40/eval10/ascending` | construction | interpn | 14.4 µs | 1.27 ns | 11323.08x | 2e-16 |
| `grid/linear/2d/max200/eval100/ascending` | construction | interpn | 7.49 µs | 0 | 10957.71x | 1e-16 |
| `grid/linear/2d/max200/eval10/ascending` | construction | interpn | 7.47 µs | 0 | 10934.86x | 2e-16 |
| `grid/cubic/2d/each40/eval1000/ascending` | construction | interpn | 9.12 µs | 0 | 10371.56x | 7e-01 |
| `grid/cubic/2d/each40/eval10000/ascending` | construction | interpn | 9.05 µs | 0 | 10296.89x | 7e-01 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | interpn | 6.85 µs | 0 | 10025.14x | 2e-16 |
| `grid/linear/4d/max40/eval100/ascending` | construction | interpn | 14.3 µs | 1.46 ns | 9779.20x | 2e-16 |
| `grid/linear/4d/max40/eval10000/ascending` | construction | interpn | 13.1 µs | 1.37 ns | 9581.71x | 5e-16 |
| `grid/linear/4d/max200/eval1000/ascending` | construction | SciPy | 50 ms | 9.47 µs | 5285.69x |  |
| `grid/linear/4d/max200/eval10/flipped` | construction | SciPy | 105 ms | 20.1 µs | 5237.92x |  |
| `grid/linear/4d/max200/eval10000/ascending` | construction | SciPy | 48.2 ms | 9.28 µs | 5198.60x |  |
| `grid/linear/4d/max200/eval100/ascending` | construction | SciPy | 48.3 ms | 9.43 µs | 5123.66x |  |
| `grid/linear/4d/max200/eval1000/flipped` | construction | SciPy | 103 ms | 20.4 µs | 5038.62x |  |
| `grid/linear/4d/max200/eval100/flipped` | construction | SciPy | 101 ms | 20.1 µs | 5009.46x |  |
| `grid/linear/4d/max200/eval10000/flipped` | construction | SciPy | 98.6 ms | 20 µs | 4927.06x |  |
| `grid/linear/4d/max200/eval10/ascending` | construction | SciPy | 46.3 ms | 9.51 µs | 4872.46x |  |
| `grid/linear/3d/max40/eval10/ascending` | construction | interpn | 2.79 µs | 0 | 3170.67x | 2e-16 |
| `grid/linear/3d/max40/eval100/ascending` | construction | interpn | 2.78 µs | 0 | 3165.33x | 2e-16 |
| `grid/cubic/2d/each10/eval1000/ascending` | construction | interpn | 1.62 µs | 0 | 1654.00x | 3e-01 |
| `grid/nearest/2d/each40/eval10000/ascending` | construction | interpn | 586 ns | 0 | 750.50x | 0e+00 |
| `grid/linear/2d/max40/eval10/ascending` | construction | interpn | 360 ns | 0 | 526.86x | 2e-16 |
| `grid/linear/2d/max40/eval10000/ascending` | construction | interpn | 358 ns | 0 | 524.00x | 4e-16 |
| `grid/linear/2d/max40/eval1000/ascending` | construction | interpn | 357 ns | 0 | 522.86x | 3e-16 |
| `grid/linear/2d/max40/eval100/ascending` | construction | interpn | 356 ns | 0 | 521.14x | 2e-16 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | ndarray-interp | 1.14 ms | 2.22 µs | 511.24x | 6e-16 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | ndarray-interp | 312 µs | 1.16 µs | 269.21x | 5e-16 |
| `grid/linear/3d/max200/eval10/ascending` | construction | SciPy | 1.73 ms | 7.51 µs | 230.17x |  |
| `grid/linear/3d/max200/eval100/ascending` | construction | SciPy | 1.71 ms | 7.64 µs | 223.54x |  |
| `grid/linear/3d/max10/eval10/ascending` | construction | interpn | 193 ns | 0 | 219.44x | 2e-16 |
| `grid/linear/3d/max10/eval100/ascending` | construction | interpn | 191 ns | 0 | 217.44x | 3e-16 |
| `grid/linear/3d/max200/eval10000/ascending` | construction | SciPy | 1.75 ms | 8.19 µs | 213.22x |  |
| `grid/linear/3d/max200/eval1000/ascending` | construction | SciPy | 1.69 ms | 8.3 µs | 204.18x |  |
| `grid/linear/3d/max10/eval1000/ascending` | construction | interpn | 177 ns | 0 | 201.22x | 5e-16 |
| `grid/linear/3d/max200/eval1000/flipped` | construction | SciPy | 3.21 ms | 15.9 µs | 201.19x |  |
| `grid/linear/3d/max200/eval100/flipped` | construction | SciPy | 3.22 ms | 16 µs | 200.87x |  |
| `grid/linear/3d/max200/eval10000/flipped` | construction | SciPy | 3.2 ms | 16 µs | 199.83x |  |
| `grid/linear/2d/max10/eval1000/ascending` | construction | interpn | 135 ns | 0 | 197.14x | 3e-16 |
| `grid/linear/3d/max10/eval10000/ascending` | construction | interpn | 192 ns | 0 | 197.10x | 4e-16 |
| `grid/linear/2d/max10/eval10/ascending` | construction | interpn | 134 ns | 0 | 195.71x | 1e-16 |
| `grid/linear/2d/max10/eval100/ascending` | construction | interpn | 124 ns | 0 | 182.00x | 2e-16 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | SciPy | 1.14 ms | 6.34 µs | 179.34x |  |
| `grid/linear/2d/max10/eval10000/ascending` | construction | interpn | 122 ns | 0 | 179.00x | 3e-16 |
| `grid/linear/3d/max200/eval10/flipped` | construction | SciPy | 3.04 ms | 17 µs | 178.62x |  |
| `grid/linear/2d/max200/eval10000/ascending` | construction | ndarray-interp | 47.2 µs | 470 ns | 100.36x | 5e-16 |
| `grid/pchip/3d/each40/eval10000/ascending` | construction | SciPy | 509 µs | 8.55 µs | 59.51x |  |
| `grid/linear/2d/each400/eval10000/ascending` | construction | SciPy | 312 µs | 6.03 µs | 51.75x |  |
| `grid/linear/2d/max200/eval10/ascending` | construction | ndarray-interp | 7.47 µs | 470 ns | 15.92x | 2e-16 |
| `grid/linear/2d/max200/eval100/ascending` | construction | ndarray-interp | 7.49 µs | 475 ns | 15.78x | 3e-16 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | ndarray-interp | 6.85 µs | 515 ns | 13.30x | 4e-16 |
| `grid/linear/2d/max200/eval10000/ascending` | construction | SciPy | 47.2 µs | 5.88 µs | 8.03x |  |
| `oned/linear/n10000/sorted` | construction | ndarray-interp | 99.4 µs | 15 µs | 6.63x | 0e+00 |
| `oned/linear/n10000/shuffled` | construction | ndarray-interp | 91.1 µs | 15.3 µs | 5.95x | 0e+00 |
| `grid/linear/2d/max200/eval10000/flipped` | construction | SciPy | 51.7 µs | 11.8 µs | 4.38x |  |
| `oned/linear/n10000/sorted-uneven` | construction | ndarray-interp | 94.1 µs | 23.1 µs | 4.07x | 0e+00 |
| `oned/linear/n10000/shuffled-uneven` | construction | ndarray-interp | 97.4 µs | 31.9 µs | 3.06x | 0e+00 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | SciPy | 22.6 µs | 7.65 µs | 2.95x |  |
| `oned/pchip/n1000/shuffled` | construction | scirs2-interpolate | 14.6 µs | 5.24 µs | 2.78x | 2e-16 |
| `oned/pchip/n1000/sorted` | construction | scirs2-interpolate | 12.8 µs | 4.59 µs | 2.78x | 2e-16 |
| `oned/bspline1/n10000/sorted` | construction | SciPy | 67.3 µs | 25 µs | 2.69x |  |
| `grid/linear/4d/max40/eval1000/flipped` | construction | SciPy | 52.9 µs | 19.8 µs | 2.68x |  |
| `oned/bspline1/n10000/shuffled` | construction | SciPy | 65.3 µs | 25 µs | 2.61x |  |
| `grid/linear/4d/max40/eval10/flipped` | construction | SciPy | 51.4 µs | 19.8 µs | 2.60x |  |
| `grid/linear/4d/max40/eval100/flipped` | construction | SciPy | 57 µs | 22.2 µs | 2.57x |  |
| `grid/linear/2d/max200/eval10/ascending` | construction | scirs2-interpolate | 7.47 µs | 2.96 µs | 2.52x | 1e-16 |
| `oned/cubic_spline/n1000/shuffled` | construction | ndarray-interp | 42.7 µs | 17.3 µs | 2.47x | 2e-16 |
| `grid/linear/4d/max40/eval10/ascending` | construction | scirs2-interpolate | 14.4 µs | 5.86 µs | 2.45x | 2e-16 |
| `grid/linear/4d/max40/eval10000/flipped` | construction | SciPy | 52.7 µs | 21.6 µs | 2.44x |  |
| `grid/linear/2d/max40/eval1000/ascending` | construction | ndarray-interp | 357 ns | 147 ns | 2.43x | 4e-16 |
| `grid/linear/2d/max200/eval100/ascending` | construction | scirs2-interpolate | 7.49 µs | 3.09 µs | 2.42x | 2e-16 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | scirs2-interpolate | 22.6 µs | 9.38 µs | 2.41x | 0e+00 |
| `grid/linear/2d/max40/eval10000/ascending` | construction | ndarray-interp | 358 ns | 150 ns | 2.39x | 4e-16 |
| `grid/linear/4d/max40/eval100/ascending` | construction | scirs2-interpolate | 14.3 µs | 6.02 µs | 2.38x | 3e-16 |
| `grid/linear/4d/max40/eval1000/ascending` | construction | scirs2-interpolate | 14.4 µs | 6.05 µs | 2.38x | 2e-16 |
| `grid/linear/2d/max40/eval100/ascending` | construction | ndarray-interp | 356 ns | 150 ns | 2.38x | 2e-16 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | scirs2-interpolate | 6.85 µs | 2.97 µs | 2.31x | 2e-16 |
| `oned/cubic_spline/n1000/shuffled` | construction | scirs2-interpolate | 42.7 µs | 18.5 µs | 2.31x | 4e-04 |
| `oned/cubic_spline/n1000/sorted` | construction | scirs2-interpolate | 40.9 µs | 18.1 µs | 2.26x | 4e-04 |
| `oned/linear/n10/shuffled` | construction | ndarray-interp | 94.9 ns | 43.3 ns | 2.19x | 0e+00 |
| `grid/linear/2d/max40/eval10/ascending` | construction | ndarray-interp | 360 ns | 164 ns | 2.19x | 2e-16 |
| `grid/linear/4d/max40/eval10000/ascending` | construction | scirs2-interpolate | 13.1 µs | 5.99 µs | 2.19x | 4e-16 |
| `grid/nearest/2d/each40/eval10000/ascending` | construction | scirs2-interpolate | 586 ns | 283 ns | 2.07x | 0e+00 |
| `oned/linear/n10000/shuffled-uneven` | construction | SciPy | 97.4 µs | 47.9 µs | 2.03x |  |
| `oned/linear/n10000/sorted-uneven` | construction | SciPy | 94.1 µs | 48.8 µs | 1.93x |  |
| `oned/linear/n10/sorted` | construction | ndarray-interp | 92.4 ns | 47.9 ns | 1.93x | 0e+00 |
| `oned/linear/n10000/sorted` | construction | SciPy | 99.4 µs | 52.2 µs | 1.90x |  |
| `oned/linear/n10000/shuffled` | construction | SciPy | 91.1 µs | 48.3 µs | 1.89x |  |
| `oned/pchip/n100/shuffled` | construction | scirs2-interpolate | 969 ns | 536 ns | 1.81x | 2e-16 |
| `oned/pchip/n100/sorted` | construction | scirs2-interpolate | 968 ns | 536 ns | 1.80x | 2e-16 |
| `grid/linear/3d/max40/eval1000/flipped` | construction | SciPy | 27.7 µs | 15.6 µs | 1.77x |  |
| `grid/linear/2d/max10/eval10/ascending` | construction | ndarray-interp | 134 ns | 77.8 ns | 1.72x | 1e-16 |
| `grid/linear/2d/max40/eval100/ascending` | construction | scirs2-interpolate | 356 ns | 209 ns | 1.71x | 3e-16 |
| `grid/linear/3d/max40/eval1000/ascending` | construction | scirs2-interpolate | 15.2 µs | 8.89 µs | 1.70x | 2e-16 |
| `grid/linear/2d/max40/eval1000/ascending` | construction | scirs2-interpolate | 357 ns | 212 ns | 1.69x | 3e-16 |
| `grid/linear/2d/max40/eval10/ascending` | construction | scirs2-interpolate | 360 ns | 214 ns | 1.68x | 4e-17 |
| `oned/linear/n1000/sorted` | construction | scirs2-interpolate | 2.1 µs | 1.26 µs | 1.67x | 0e+00 |
| `oned/pchip/n100/shuffled-uneven` | construction | scirs2-interpolate | 966 ns | 581 ns | 1.66x | 2e-16 |
| `grid/linear/2d/max40/eval10000/ascending` | construction | scirs2-interpolate | 358 ns | 215 ns | 1.66x | 2e-16 |
| `oned/linear/n1000/shuffled` | construction | scirs2-interpolate | 1.93 µs | 1.17 µs | 1.65x | 0e+00 |
| `oned/linear/n10000/shuffled-uneven` | construction | scirs2-interpolate | 97.4 µs | 60 µs | 1.62x | 6e-17 |
| `oned/pchip/n50/shuffled` | construction | scirs2-interpolate | 542 ns | 336 ns | 1.61x | 1e-16 |
| `oned/linear/n100/shuffled` | construction | ndarray-interp | 270 ns | 168 ns | 1.60x | 0e+00 |
| `grid/linear/2d/max10/eval100/ascending` | construction | ndarray-interp | 124 ns | 77.8 ns | 1.60x | 4e-16 |
| `oned/pchip/n100/sorted-uneven` | construction | scirs2-interpolate | 883 ns | 554 ns | 1.59x | 2e-16 |
| `oned/linear/n10000/sorted-uneven` | construction | scirs2-interpolate | 94.1 µs | 59.2 µs | 1.59x | 6e-17 |
| `oned/linear/n10000/sorted` | construction | scirs2-interpolate | 99.4 µs | 63 µs | 1.58x | 0e+00 |
| `grid/linear/2d/max10/eval1000/ascending` | construction | ndarray-interp | 135 ns | 85.4 ns | 1.58x | 3e-16 |
| `grid/linear/2d/max10/eval10000/ascending` | construction | ndarray-interp | 122 ns | 78 ns | 1.57x | 4e-16 |
| `oned/pchip/n50/sorted` | construction | scirs2-interpolate | 542 ns | 346 ns | 1.56x | 1e-16 |
| `oned/linear/n50/shuffled` | construction | ndarray-interp | 167 ns | 108 ns | 1.55x | 0e+00 |
| `grid/linear/4d/max40/eval100/ascending` | construction | SciPy | 14.3 µs | 9.51 µs | 1.51x |  |
| `grid/linear/3d/max10/eval10/ascending` | evaluation | interpn | 215 ns | 143 ns | 1.50x | 2e-16 |
| `oned/linear/n50/sorted` | construction | ndarray-interp | 152 ns | 102 ns | 1.49x | 0e+00 |
| `oned/linear/n10000/shuffled` | construction | scirs2-interpolate | 91.1 µs | 61.4 µs | 1.48x | 0e+00 |
| `grid/linear/2d/max10/eval1000/ascending` | evaluation | interpn | 13.7 µs | 9.27 µs | 1.47x | 3e-16 |
| `oned/linear/n100/sorted` | construction | ndarray-interp | 247 ns | 168 ns | 1.47x | 0e+00 |
| `grid/linear/3d/max200/eval10000/ascending` | construction | scirs2-interpolate | 1.75 ms | 1.2 ms | 1.46x | 3e-16 |
| `oned/linear/n100/shuffled` | construction | scirs2-interpolate | 270 ns | 190 ns | 1.42x | 0e+00 |
| `grid/linear/2d/max10/eval10/ascending` | evaluation | interpn | 135 ns | 95.1 ns | 1.42x | 1e-16 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | scirs2-interpolate | 1.14 ms | 805 µs | 1.41x | 2e-16 |
| `grid/linear/4d/max40/eval10/ascending` | construction | SciPy | 14.4 µs | 10.2 µs | 1.41x |  |
| `oned/linear/n100/shuffled-uneven` | construction | ndarray-interp | 259 ns | 184 ns | 1.41x | 0e+00 |
| `grid/linear/4d/max40/eval1000/ascending` | construction | SciPy | 14.4 µs | 10.2 µs | 1.40x |  |
| `oned/linear/n1000/shuffled` | construction | ndarray-interp | 1.93 µs | 1.38 µs | 1.39x | 0e+00 |
| `oned/akima/n100/shuffled` | construction | scirs2-interpolate | 743 ns | 534 ns | 1.39x | 1e-16 |
| `oned/linear/n1000/sorted` | construction | ndarray-interp | 2.1 µs | 1.52 µs | 1.38x | 0e+00 |
| `grid/linear/4d/max40/eval10000/ascending` | construction | SciPy | 13.1 µs | 9.47 µs | 1.38x |  |
| `grid/linear/3d/max200/eval10/ascending` | construction | scirs2-interpolate | 1.73 ms | 1.26 ms | 1.37x | 3e-16 |
| `oned/linear/n100/sorted-uneven` | construction | ndarray-interp | 252 ns | 185 ns | 1.37x | 0e+00 |
| `oned/akima/n50/sorted` | construction | scirs2-interpolate | 469 ns | 347 ns | 1.35x | 1e-16 |
| `grid/linear/3d/max10/eval1000/ascending` | evaluation | interpn | 19 µs | 14.1 µs | 1.35x | 5e-16 |
| `grid/linear/3d/max10/eval10000/ascending` | evaluation | interpn | 192 µs | 142 µs | 1.35x | 4e-16 |
| `grid/linear/2d/max10/eval10000/ascending` | evaluation | interpn | 125 µs | 92.4 µs | 1.35x | 3e-16 |
| `grid/linear/3d/max10/eval100/ascending` | evaluation | interpn | 1.88 µs | 1.39 µs | 1.35x | 3e-16 |
| `oned/akima/n100/sorted` | construction | scirs2-interpolate | 753 ns | 577 ns | 1.31x | 1e-16 |
| `grid/linear/3d/max200/eval1000/ascending` | construction | scirs2-interpolate | 1.69 ms | 1.3 ms | 1.30x | 3e-16 |
| `oned/linear/n100/sorted` | construction | scirs2-interpolate | 247 ns | 190 ns | 1.30x | 0e+00 |
| `grid/linear/3d/max40/eval10000/ascending` | construction | scirs2-interpolate | 16.3 µs | 12.5 µs | 1.30x | 2e-16 |
| `oned/akima/n50/shuffled` | construction | scirs2-interpolate | 437 ns | 337 ns | 1.29x | 1e-16 |
| `grid/linear/2d/max200/eval10/ascending` | construction | SciPy | 7.47 µs | 5.82 µs | 1.28x |  |
| `grid/linear/4d/max200/eval1000/ascending` | construction | scirs2-interpolate | 50 ms | 39.2 ms | 1.27x | 3e-16 |
| `grid/linear/4d/max40/eval10/ascending` | evaluation | interpn | 340 ns | 270 ns | 1.26x | 2e-16 |
| `oned/linear/n100/sorted-uneven` | construction | scirs2-interpolate | 252 ns | 203 ns | 1.25x | 3e-17 |
| `oned/linear/n100/shuffled-uneven` | construction | scirs2-interpolate | 259 ns | 208 ns | 1.25x | 3e-17 |
| `grid/linear/2d/max200/eval1000/flipped` | construction | SciPy | 14.8 µs | 11.9 µs | 1.24x |  |
| `grid/linear/4d/max200/eval10000/ascending` | construction | scirs2-interpolate | 48.2 ms | 39 ms | 1.24x | 4e-16 |
| `grid/linear/2d/max10/eval100/ascending` | evaluation | interpn | 1.25 µs | 1.02 µs | 1.23x | 2e-16 |
| `oned/akima/n1000/sorted` | construction | scirs2-interpolate | 10.9 µs | 8.9 µs | 1.22x | 3e-17 |
| `grid/linear/2d/max200/eval10/flipped` | construction | SciPy | 14.5 µs | 11.9 µs | 1.22x |  |
| `oned/nearest/n10000/sorted` | construction | SciPy | 66.6 µs | 54.6 µs | 1.22x |  |
| `grid/linear/2d/max200/eval100/flipped` | construction | SciPy | 14.5 µs | 11.9 µs | 1.22x |  |
| `grid/linear/2d/max200/eval10000/ascending` | evaluation | interpn | 206 µs | 172 µs | 1.19x | 3e-16 |
| `oned/nearest/n10000/shuffled` | construction | SciPy | 65 µs | 54.4 µs | 1.19x |  |
| `grid/linear/2d/max200/eval100/ascending` | evaluation | interpn | 2.05 µs | 1.72 µs | 1.19x | 1e-16 |
| `grid/linear/4d/max200/eval10/ascending` | construction | scirs2-interpolate | 46.3 ms | 39.1 ms | 1.18x | 2e-16 |
| `oned/linear/n50/shuffled` | construction | scirs2-interpolate | 167 ns | 142 ns | 1.18x | 0e+00 |
| `grid/linear/2d/max200/eval100/ascending` | construction | SciPy | 7.49 µs | 6.37 µs | 1.18x |  |
| `grid/linear/2d/max40/eval10/ascending` | evaluation | interpn | 156 ns | 134 ns | 1.17x | 2e-16 |
| `grid/linear/3d/max10/eval100/ascending` | construction | scirs2-interpolate | 191 ns | 164 ns | 1.17x | 2e-16 |
| `grid/linear/3d/max10/eval10/ascending` | construction | scirs2-interpolate | 193 ns | 166 ns | 1.16x | 2e-16 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | SciPy | 6.85 µs | 5.9 µs | 1.16x |  |
| `grid/linear/4d/max200/eval10/ascending` | evaluation | interpn | 412 ns | 356 ns | 1.16x | 2e-16 |
| `grid/linear/3d/max200/eval10/ascending` | evaluation | interpn | 301 ns | 261 ns | 1.15x | 3e-16 |
| `oned/akima/n1000/shuffled` | construction | scirs2-interpolate | 10.1 µs | 8.84 µs | 1.15x | 3e-17 |
| `oned/linear/n50/sorted` | construction | scirs2-interpolate | 152 ns | 133 ns | 1.15x | 0e+00 |
| `grid/linear/2d/max40/eval10000/ascending` | evaluation | interpn | 151 µs | 132 µs | 1.15x | 4e-16 |
| `grid/linear/4d/max200/eval100/ascending` | construction | scirs2-interpolate | 48.3 ms | 42.2 ms | 1.15x | 1e-16 |
| `grid/linear/2d/max40/eval100/ascending` | evaluation | interpn | 1.51 µs | 1.32 µs | 1.14x | 2e-16 |
| `grid/linear/4d/max40/eval1000/ascending` | evaluation | interpn | 32.5 µs | 28.5 µs | 1.14x | 6e-16 |
| `grid/linear/3d/max200/eval100/ascending` | construction | scirs2-interpolate | 1.71 ms | 1.51 ms | 1.13x | 3e-16 |
| `grid/linear/4d/max200/eval10000/ascending` | evaluation | interpn | 411 µs | 365 µs | 1.12x | 4e-16 |
| `oned/nearest/n10/shuffled` | evaluation | scirs2-interpolate | 299 ns | 266 ns | 1.12x | 0e+00 |
| `oned/cubic_spline/n100/sorted` | construction | scirs2-interpolate | 1.73 µs | 1.55 µs | 1.12x | 5e-03 |
| `grid/nearest/2d/each40/eval10000/ascending` | evaluation | interpn | 150 µs | 136 µs | 1.11x | 0e+00 |

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

| case | interpx | interpx, parallel | SciPy | scirs2-interpolate | ndarray-interp | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/n10/sorted` | 112 ns | 189 ns | 4.91 µs | 262 ns | 277 ns | 2.3x faster |
| `linear/n10/shuffled` | 171 ns | 224 ns | 5.06 µs | 286 ns | 249 ns | 1.5x faster |
| `nearest/n10/sorted` | 148 ns | 213 ns | 9.47 µs | 280 ns |  | 1.9x faster |
| `nearest/n10/shuffled` | 299 ns | 361 ns | 10.2 µs | 266 ns |  | **1.1x slower** |
| `previous/n10/sorted` | 162 ns | 195 ns | 10.6 µs |  |  | 66x faster |
| `previous/n10/shuffled` | 264 ns | 352 ns | 10.7 µs |  |  | 41x faster |
| `bspline1/n10/sorted` | 219 ns | 315 ns | 2.69 µs |  |  | 12x faster |
| `bspline1/n10/shuffled` | 578 ns | 656 ns | 2.74 µs |  |  | 4.7x faster |
| `bspline2/n10/sorted` | 273 ns | 346 ns | 3.29 µs |  |  | 12x faster |
| `bspline2/n10/shuffled` | 340 ns | 416 ns | 3.1 µs |  |  | 9.1x faster |
| `bspline3/n10/sorted` | 346 ns | 436 ns | 3.64 µs |  |  | 11x faster |
| `bspline3/n10/shuffled` | 440 ns | 532 ns | 3.6 µs |  |  | 8.2x faster |
| `cubic_spline/n10/sorted` | 143 ns | 211 ns | 3.58 µs | 271 ns | 323 ns | 1.9x faster |
| `cubic_spline/n10/shuffled` | 204 ns | 290 ns | 4.03 µs | 278 ns | 325 ns | 1.4x faster |
| `pchip/n10/sorted` | 143 ns | 199 ns | 3.65 µs | 298 ns |  | 2.1x faster |
| `pchip/n10/shuffled` | 206 ns | 262 ns | 3.65 µs | 315 ns |  | 1.5x faster |
| `akima/n10/sorted` | 143 ns | 194 ns | 3.67 µs | 274 ns |  | 1.9x faster |
| `akima/n10/shuffled` | 205 ns | 257 ns | 3.67 µs | 285 ns |  | 1.4x faster |
| `makima/n10/sorted` | 143 ns | 195 ns | 3.64 µs |  |  | 25x faster |
| `makima/n10/shuffled` | 204 ns | 271 ns | 3.65 µs |  |  | 18x faster |
| `linear/n50/sorted` | 639 ns | 750 ns | 6.13 µs | 1.53 µs | 1.35 µs | 2.1x faster |
| `linear/n50/shuffled` | 891 ns | 946 ns | 8.21 µs | 1.54 µs | 1.38 µs | 1.6x faster |
| `nearest/n50/sorted` | 747 ns | 813 ns | 11.2 µs | 1.7 µs |  | 2.3x faster |
| `nearest/n50/shuffled` | 1.53 µs | 1.6 µs | 11.3 µs | 1.57 µs |  | 1x faster |
| `previous/n50/sorted` | 731 ns | 836 ns | 12.8 µs |  |  | 18x faster |
| `previous/n50/shuffled` | 1.49 µs | 1.42 µs | 13.9 µs |  |  | 9.3x faster |
| `bspline1/n50/sorted` | 1.1 µs | 1.39 µs | 3.94 µs |  |  | 3.6x faster |
| `bspline1/n50/shuffled` | 2.78 µs | 3.13 µs | 4.87 µs |  |  | 1.8x faster |
| `bspline2/n50/sorted` | 1.36 µs | 1.63 µs | 6.09 µs |  |  | 4.5x faster |
| `bspline2/n50/shuffled` | 2.02 µs | 2.1 µs | 6.78 µs |  |  | 3.4x faster |
| `bspline3/n50/sorted` | 1.83 µs | 1.89 µs | 8.72 µs |  |  | 4.8x faster |
| `bspline3/n50/shuffled` | 3.08 µs | 3.16 µs | 9.31 µs |  |  | 3x faster |
| `cubic_spline/n50/sorted` | 707 ns | 753 ns | 5.26 µs | 4.34 µs | 1.58 µs | 2.2x faster |
| `cubic_spline/n50/shuffled` | 998 ns | 1.14 µs | 5.53 µs | 4.52 µs | 1.75 µs | 1.7x faster |
| `pchip/n50/sorted` | 754 ns | 754 ns | 5.63 µs | 3.98 µs |  | 5.3x faster |
| `pchip/n50/shuffled` | 1.07 µs | 1.03 µs | 5.62 µs | 3.8 µs |  | 3.6x faster |
| `akima/n50/sorted` | 731 ns | 759 ns | 5.21 µs | 4.13 µs |  | 5.7x faster |
| `akima/n50/shuffled` | 979 ns | 1.06 µs | 5.42 µs | 4.06 µs |  | 4.1x faster |
| `makima/n50/sorted` | 712 ns | 809 ns | 5.16 µs |  |  | 7.3x faster |
| `makima/n50/shuffled` | 1.06 µs | 1.08 µs | 5.46 µs |  |  | 5.2x faster |
| `linear/n100/sorted` | 1.33 µs | 1.43 µs | 7.88 µs | 3.61 µs | 2.46 µs | 1.8x faster |
| `linear/n100/shuffled` | 1.62 µs | 1.88 µs | 12.7 µs | 3.36 µs | 2.49 µs | 1.5x faster |
| `nearest/n100/sorted` | 1.5 µs | 1.59 µs | 13.4 µs | 3.33 µs |  | 2.2x faster |
| `nearest/n100/shuffled` | 3.07 µs | 3.51 µs | 13.7 µs | 3.73 µs |  | 1.2x faster |
| `previous/n100/sorted` | 1.65 µs | 1.72 µs | 16.3 µs |  |  | 9.9x faster |
| `previous/n100/shuffled` | 2.97 µs | 3.12 µs | 16.8 µs |  |  | 5.7x faster |
| `bspline1/n100/sorted` | 2.36 µs | 3.06 µs | 5.81 µs |  |  | 2.5x faster |
| `bspline1/n100/shuffled` | 6.12 µs | 6.32 µs | 12.9 µs |  |  | 2.1x faster |
| `bspline2/n100/sorted` | 2.95 µs | 3.63 µs | 10.8 µs |  |  | 3.7x faster |
| `bspline2/n100/shuffled` | 4.86 µs | 4.7 µs | 13.8 µs |  |  | 2.8x faster |
| `bspline3/n100/sorted` | 3.62 µs | 3.94 µs | 16.4 µs |  |  | 4.5x faster |
| `bspline3/n100/shuffled` | 6.17 µs | 6.88 µs | 17.9 µs |  |  | 2.9x faster |
| `cubic_spline/n100/sorted` | 1.47 µs | 1.54 µs | 7.27 µs | 18.9 µs | 3.25 µs | 2.2x faster |
| `cubic_spline/n100/shuffled` | 1.94 µs | 2.08 µs | 7.88 µs | 17.6 µs | 3.2 µs | 1.6x faster |
| `pchip/n100/sorted` | 1.46 µs | 1.54 µs | 7.37 µs | 15.9 µs |  | 5x faster |
| `pchip/n100/shuffled` | 1.94 µs | 2.25 µs | 7.85 µs | 16.3 µs |  | 4x faster |
| `akima/n100/sorted` | 1.49 µs | 1.62 µs | 7.42 µs | 16.6 µs |  | 5x faster |
| `akima/n100/shuffled` | 1.95 µs | 2.25 µs | 7.92 µs | 17.1 µs |  | 4.1x faster |
| `makima/n100/sorted` | 1.46 µs | 1.54 µs | 7.05 µs |  |  | 4.8x faster |
| `makima/n100/shuffled` | 1.94 µs | 2.26 µs | 7.95 µs |  |  | 4.1x faster |
| `linear/n1000/sorted` | 14.5 µs | 14.6 µs | 29.3 µs | 95.9 µs | 25.5 µs | 1.8x faster |
| `linear/n1000/shuffled` | 16.1 µs | 17.8 µs | 144 µs | 202 µs | 25.9 µs | 1.6x faster |
| `nearest/n1000/sorted` | 14.1 µs | 15 µs | 55.8 µs | 130 µs |  | 4x faster |
| `nearest/n1000/shuffled` | 30.9 µs | 31.1 µs | 200 µs | 228 µs |  | 6.5x faster |
| `previous/n1000/sorted` | 14.2 µs | 15.8 µs | 59.7 µs |  |  | 4.2x faster |
| `previous/n1000/shuffled` | 27.3 µs | 29.9 µs | 221 µs |  |  | 8.1x faster |
| `bspline1/n1000/sorted` | 26 µs | 26.8 µs | 31.4 µs |  |  | 1.2x faster |
| `bspline1/n1000/shuffled` | 61.7 µs | 62.1 µs | 532 µs |  |  | 8.6x faster |
| `bspline2/n1000/sorted` | 35.3 µs | 34.8 µs | 83 µs |  |  | 2.4x faster |
| `bspline2/n1000/shuffled` | 63.8 µs | 64.2 µs | 550 µs |  |  | 8.6x faster |
| `bspline3/n1000/sorted` | 39.7 µs | 40 µs | 128 µs |  |  | 3.2x faster |
| `bspline3/n1000/shuffled` | 86.4 µs | 85.2 µs | 631 µs |  |  | 7.3x faster |
| `cubic_spline/n1000/sorted` | 13.9 µs | 16.4 µs | 37.9 µs | 1.55 ms | 34.4 µs | 2.5x faster |
| `cubic_spline/n1000/shuffled` | 21.5 µs | 21.8 µs | 239 µs | 1.56 ms | 34.9 µs | 1.6x faster |
| `pchip/n1000/sorted` | 15.4 µs | 15.6 µs | 40.1 µs | 1.32 ms |  | 2.6x faster |
| `pchip/n1000/shuffled` | 21.5 µs | 21.6 µs | 238 µs | 1.34 ms |  | 11x faster |
| `akima/n1000/sorted` | 13.9 µs | 14 µs | 37.3 µs | 1.55 ms |  | 2.7x faster |
| `akima/n1000/shuffled` | 19.9 µs | 21.6 µs | 219 µs | 1.56 ms |  | 11x faster |
| `makima/n1000/sorted` | 14.7 µs | 15.1 µs | 37.6 µs |  |  | 2.6x faster |
| `makima/n1000/shuffled` | 21.5 µs | 19.9 µs | 240 µs |  |  | 11x faster |
| `linear/n10000/sorted` | 154 µs | 65.5 µs | 237 µs | 1.41 ms | 272 µs | 1.5x faster |
| `linear/n10000/shuffled` | 189 µs | 83.2 µs | 2.38 ms | 3.14 ms | 285 µs | 1.5x faster |
| `nearest/n10000/sorted` | 165 µs | 84.7 µs | 590 µs | 1.63 ms |  | 3.6x faster |
| `nearest/n10000/shuffled` | 416 µs | 94.9 µs | 3.16 ms | 3.14 ms |  | 7.5x faster |
| `previous/n10000/sorted` | 153 µs | 75.4 µs | 611 µs |  |  | 4x faster |
| `previous/n10000/shuffled` | 377 µs | 90.2 µs | 3.21 ms |  |  | 8.5x faster |
| `bspline1/n10000/sorted` | 281 µs | 94.7 µs | 288 µs |  |  | 1x faster |
| `bspline1/n10000/shuffled` | 698 µs | 141 µs | 51 ms |  |  | 73x faster |
| `bspline2/n10000/sorted` | 362 µs | 92.7 µs | 773 µs |  |  | 2.1x faster |
| `bspline2/n10000/shuffled` | 955 µs | 214 µs | 51.8 ms |  |  | 54x faster |
| `bspline3/n10000/sorted` | 402 µs | 109 µs | 1.31 ms |  |  | 3.3x faster |
| `bspline3/n10000/shuffled` | 1.13 ms | 251 µs | 49.6 ms |  |  | 44x faster |
| `cubic_spline/n10000/sorted` | 146 µs | 72.7 µs | 372 µs | 151 ms | 350 µs | 2.4x faster |
| `cubic_spline/n10000/shuffled` | 236 µs | 96.6 µs | 3.52 ms | 145 ms | 333 µs | 1.4x faster |
| `pchip/n10000/sorted` | 148 µs | 71 µs | 371 µs | 125 ms |  | 2.5x faster |
| `pchip/n10000/shuffled` | 216 µs | 98.2 µs | 3.57 ms | 126 ms |  | 17x faster |
| `akima/n10000/sorted` | 159 µs | 71 µs | 380 µs | 150 ms |  | 2.4x faster |
| `akima/n10000/shuffled` | 216 µs | 102 µs | 3.62 ms | 158 ms |  | 17x faster |
| `makima/n10000/sorted` | 156 µs | 75.1 µs | 404 µs |  |  | 2.6x faster |
| `makima/n10000/shuffled` | 216 µs | 96 µs | 3.64 ms |  |  | 17x faster |
| `linear/n100/sorted-uneven` | 1.28 µs | 1.29 µs | 7.35 µs | 3.62 µs | 4.12 µs | 2.8x faster |
| `linear/n100/shuffled-uneven` | 2.75 µs | 3.15 µs | 12.6 µs | 3.65 µs | 4.13 µs | 1.3x faster |
| `cubic_spline/n100/sorted-uneven` | 1.47 µs | 1.49 µs | 8.1 µs | 18.9 µs | 4.83 µs | 3.3x faster |
| `cubic_spline/n100/shuffled-uneven` | 3.36 µs | 3.36 µs | 7.94 µs | 19.4 µs | 4.78 µs | 1.4x faster |
| `pchip/n100/sorted-uneven` | 1.35 µs | 1.57 µs | 7.92 µs | 17.2 µs |  | 5.8x faster |
| `pchip/n100/shuffled-uneven` | 3.22 µs | 3.43 µs | 7.89 µs | 17.8 µs |  | 2.5x faster |
| `linear/n10000/sorted-uneven` | 228 µs | 91.7 µs | 280 µs | 1.28 ms | 901 µs | 1.2x faster |
| `linear/n10000/shuffled-uneven` | 624 µs | 156 µs | 2.4 ms | 2.9 ms | 1.77 ms | 2.8x faster |
| `cubic_spline/n10000/sorted-uneven` | 234 µs | 92.2 µs | 439 µs | 153 ms | 976 µs | 1.9x faster |
| `cubic_spline/n10000/shuffled-uneven` | 730 µs | 180 µs | 3.5 ms | 158 ms | 1.92 ms | 2.6x faster |
| `pchip/n10000/sorted-uneven` | 235 µs | 94.2 µs | 439 µs | 128 ms |  | 1.9x faster |
| `pchip/n10000/shuffled-uneven` | 729 µs | 167 µs | 3.43 ms | 123 ms |  | 4.7x faster |

### Construction

| case | interpx | interpx, parallel | SciPy | scirs2-interpolate | ndarray-interp | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/n10/sorted` | 92.4 ns | 101 ns | 8.8 µs | 88.1 ns | 47.9 ns | **1.9x slower** |
| `linear/n10/shuffled` | 94.9 ns | 92.6 ns | 8.96 µs | 89.6 ns | 43.3 ns | **2.2x slower** |
| `nearest/n10/sorted` | 70.6 ns | 70.2 ns | 9.77 µs | 82.8 ns |  | 1.2x faster |
| `nearest/n10/shuffled` | 70.3 ns | 70.8 ns | 10.4 µs | 92.3 ns |  | 1.3x faster |
| `previous/n10/sorted` | 70.4 ns | 70.9 ns | 10 µs |  |  | 1.4e+02x faster |
| `previous/n10/shuffled` | 64.4 ns | 70.7 ns | 9.25 µs |  |  | 1.4e+02x faster |
| `bspline1/n10/sorted` | 113 ns | 113 ns | 15.6 µs |  |  | 1.4e+02x faster |
| `bspline1/n10/shuffled` | 112 ns | 113 ns | 15.7 µs |  |  | 1.4e+02x faster |
| `bspline2/n10/sorted` | 491 ns | 532 ns | 36.5 µs |  |  | 74x faster |
| `bspline2/n10/shuffled` | 496 ns | 492 ns | 36.6 µs |  |  | 74x faster |
| `bspline3/n10/sorted` | 572 ns | 530 ns | 37.9 µs |  |  | 66x faster |
| `bspline3/n10/shuffled` | 573 ns | 569 ns | 36 µs |  |  | 63x faster |
| `cubic_spline/n10/sorted` | 234 ns | 239 ns | 92 µs | 416 ns | 478 ns | 1.8x faster |
| `cubic_spline/n10/shuffled` | 237 ns | 244 ns | 85.2 µs | 384 ns | 476 ns | 1.6x faster |
| `pchip/n10/sorted` | 181 ns | 195 ns | 101 µs | 190 ns |  | 1.1x faster |
| `pchip/n10/shuffled` | 178 ns | 195 ns | 102 µs | 201 ns |  | 1.1x faster |
| `akima/n10/sorted` | 182 ns | 182 ns | 75.2 µs | 197 ns |  | 1.1x faster |
| `akima/n10/shuffled` | 183 ns | 174 ns | 73.8 µs | 197 ns |  | 1.1x faster |
| `makima/n10/sorted` | 196 ns | 187 ns | 78.2 µs |  |  | 4e+02x faster |
| `makima/n10/shuffled` | 195 ns | 179 ns | 78.1 µs |  |  | 4e+02x faster |
| `linear/n50/sorted` | 152 ns | 169 ns | 9.6 µs | 133 ns | 102 ns | **1.5x slower** |
| `linear/n50/shuffled` | 167 ns | 157 ns | 9.55 µs | 142 ns | 108 ns | **1.5x slower** |
| `nearest/n50/sorted` | 105 ns | 103 ns | 9.7 µs | 146 ns |  | 1.4x faster |
| `nearest/n50/shuffled` | 94.7 ns | 94.2 ns | 10.7 µs | 134 ns |  | 1.4x faster |
| `previous/n50/sorted` | 100 ns | 104 ns | 9.64 µs |  |  | 96x faster |
| `previous/n50/shuffled` | 104 ns | 102 ns | 10.3 µs |  |  | 99x faster |
| `bspline1/n50/sorted` | 157 ns | 145 ns | 15.5 µs |  |  | 98x faster |
| `bspline1/n50/shuffled` | 157 ns | 147 ns | 15.6 µs |  |  | 99x faster |
| `bspline2/n50/sorted` | 1.86 µs | 1.84 µs | 42.3 µs |  |  | 23x faster |
| `bspline2/n50/shuffled` | 1.84 µs | 1.84 µs | 42.5 µs |  |  | 23x faster |
| `bspline3/n50/sorted` | 2.04 µs | 2.23 µs | 38.9 µs |  |  | 19x faster |
| `bspline3/n50/shuffled` | 2.23 µs | 2.24 µs | 39.7 µs |  |  | 18x faster |
| `cubic_spline/n50/sorted` | 839 ns | 836 ns | 86.7 µs | 908 ns | 1.08 µs | 1.1x faster |
| `cubic_spline/n50/shuffled` | 834 ns | 918 ns | 87.4 µs | 934 ns | 1.18 µs | 1.1x faster |
| `pchip/n50/sorted` | 542 ns | 501 ns | 106 µs | 346 ns |  | **1.6x slower** |
| `pchip/n50/shuffled` | 542 ns | 544 ns | 113 µs | 336 ns |  | **1.6x slower** |
| `akima/n50/sorted` | 469 ns | 456 ns | 75.3 µs | 347 ns |  | **1.4x slower** |
| `akima/n50/shuffled` | 437 ns | 427 ns | 75.1 µs | 337 ns |  | **1.3x slower** |
| `makima/n50/sorted` | 468 ns | 471 ns | 78.6 µs |  |  | 1.7e+02x faster |
| `makima/n50/shuffled` | 506 ns | 509 ns | 79.1 µs |  |  | 1.6e+02x faster |
| `linear/n100/sorted` | 247 ns | 247 ns | 9.87 µs | 190 ns | 168 ns | **1.5x slower** |
| `linear/n100/shuffled` | 270 ns | 269 ns | 9.07 µs | 190 ns | 168 ns | **1.6x slower** |
| `nearest/n100/sorted` | 164 ns | 165 ns | 9.94 µs | 191 ns |  | 1.2x faster |
| `nearest/n100/shuffled` | 167 ns | 166 ns | 9.98 µs | 191 ns |  | 1.1x faster |
| `previous/n100/sorted` | 166 ns | 167 ns | 10.1 µs |  |  | 61x faster |
| `previous/n100/shuffled` | 166 ns | 152 ns | 11 µs |  |  | 66x faster |
| `bspline1/n100/sorted` | 215 ns | 214 ns | 17 µs |  |  | 79x faster |
| `bspline1/n100/shuffled` | 218 ns | 215 ns | 17 µs |  |  | 78x faster |
| `bspline2/n100/sorted` | 3.86 µs | 3.85 µs | 48.2 µs |  |  | 12x faster |
| `bspline2/n100/shuffled` | 3.84 µs | 3.85 µs | 47.2 µs |  |  | 12x faster |
| `bspline3/n100/sorted` | 4.06 µs | 4.28 µs | 47.9 µs |  |  | 12x faster |
| `bspline3/n100/shuffled` | 4.28 µs | 4.27 µs | 44.7 µs |  |  | 10x faster |
| `cubic_spline/n100/sorted` | 1.73 µs | 1.58 µs | 94.6 µs | 1.55 µs | 2.06 µs | **1.1x slower** |
| `cubic_spline/n100/shuffled` | 1.58 µs | 1.58 µs | 88.5 µs | 1.66 µs | 1.95 µs | 1.1x faster |
| `pchip/n100/sorted` | 968 ns | 966 ns | 111 µs | 536 ns |  | **1.8x slower** |
| `pchip/n100/shuffled` | 969 ns | 969 ns | 109 µs | 536 ns |  | **1.8x slower** |
| `akima/n100/sorted` | 753 ns | 748 ns | 81.6 µs | 577 ns |  | **1.3x slower** |
| `akima/n100/shuffled` | 743 ns | 750 ns | 76 µs | 534 ns |  | **1.4x slower** |
| `makima/n100/sorted` | 828 ns | 825 ns | 85.8 µs |  |  | 1e+02x faster |
| `makima/n100/shuffled` | 828 ns | 831 ns | 79 µs |  |  | 95x faster |
| `linear/n1000/sorted` | 2.1 µs | 2.02 µs | 12.5 µs | 1.26 µs | 1.52 µs | **1.7x slower** |
| `linear/n1000/shuffled` | 1.93 µs | 2.11 µs | 12.9 µs | 1.17 µs | 1.38 µs | **1.6x slower** |
| `nearest/n1000/sorted` | 1.08 µs | 1.06 µs | 14.8 µs | 1.28 µs |  | 1.2x faster |
| `nearest/n1000/shuffled` | 1.09 µs | 1 µs | 14.6 µs | 1.17 µs |  | 1.1x faster |
| `previous/n1000/sorted` | 997 ns | 1.08 µs | 16.4 µs |  |  | 16x faster |
| `previous/n1000/shuffled` | 1.08 µs | 1.08 µs | 16.5 µs |  |  | 15x faster |
| `bspline1/n1000/sorted` | 1.13 µs | 1.04 µs | 16.9 µs |  |  | 15x faster |
| `bspline1/n1000/shuffled` | 1.12 µs | 1.14 µs | 17.9 µs |  |  | 16x faster |
| `bspline2/n1000/sorted` | 37.4 µs | 36.9 µs | 108 µs |  |  | 2.9x faster |
| `bspline2/n1000/shuffled` | 37.5 µs | 37.1 µs | 112 µs |  |  | 3x faster |
| `bspline3/n1000/sorted` | 72.8 µs | 72.9 µs | 137 µs |  |  | 1.9x faster |
| `bspline3/n1000/shuffled` | 73.6 µs | 73.3 µs | 136 µs |  |  | 1.8x faster |
| `cubic_spline/n1000/sorted` | 40.9 µs | 40.8 µs | 114 µs | 18.1 µs | 17.6 µs | **2.3x slower** |
| `cubic_spline/n1000/shuffled` | 42.7 µs | 40.5 µs | 121 µs | 18.5 µs | 17.3 µs | **2.5x slower** |
| `pchip/n1000/sorted` | 12.8 µs | 9.05 µs | 131 µs | 4.59 µs |  | **2.8x slower** |
| `pchip/n1000/shuffled` | 14.6 µs | 13 µs | 134 µs | 5.24 µs |  | **2.8x slower** |
| `akima/n1000/sorted` | 10.9 µs | 10.5 µs | 97.9 µs | 8.9 µs |  | **1.2x slower** |
| `akima/n1000/shuffled` | 10.1 µs | 10.4 µs | 94.3 µs | 8.84 µs |  | **1.1x slower** |
| `makima/n1000/sorted` | 11.2 µs | 11.7 µs | 98.2 µs |  |  | 8.8x faster |
| `makima/n1000/shuffled` | 12.3 µs | 11.9 µs | 109 µs |  |  | 8.9x faster |
| `linear/n10000/sorted` | 99.4 µs | 87.5 µs | 52.2 µs | 63 µs | 15 µs | **6.6x slower** |
| `linear/n10000/shuffled` | 91.1 µs | 88.6 µs | 48.3 µs | 61.4 µs | 15.3 µs | **6x slower** |
| `nearest/n10000/sorted` | 66.6 µs | 56.6 µs | 54.6 µs | 60.8 µs |  | **1.2x slower** |
| `nearest/n10000/shuffled` | 65 µs | 62.6 µs | 54.4 µs | 62.5 µs |  | **1.2x slower** |
| `previous/n10000/sorted` | 67.7 µs | 58.4 µs | 74.3 µs |  |  | 1.1x faster |
| `previous/n10000/shuffled` | 68.5 µs | 56.3 µs | 75.2 µs |  |  | 1.1x faster |
| `bspline1/n10000/sorted` | 67.3 µs | 64.2 µs | 25 µs |  |  | **2.7x slower** |
| `bspline1/n10000/shuffled` | 65.3 µs | 60.6 µs | 25 µs |  |  | **2.6x slower** |
| `bspline2/n10000/sorted` | 517 µs | 522 µs | 824 µs |  |  | 1.6x faster |
| `bspline2/n10000/shuffled` | 528 µs | 566 µs | 772 µs |  |  | 1.5x faster |
| `bspline3/n10000/sorted` | 638 µs | 654 µs | 1 ms |  |  | 1.6x faster |
| `bspline3/n10000/shuffled` | 630 µs | 650 µs | 956 µs |  |  | 1.5x faster |
| `cubic_spline/n10000/sorted` | 245 µs | 286 µs | 780 µs | 324 µs | 388 µs | 1.3x faster |
| `cubic_spline/n10000/shuffled` | 246 µs | 351 µs | 736 µs | 280 µs | 360 µs | 1.1x faster |
| `pchip/n10000/sorted` | 89.3 µs | 196 µs | 259 µs | 163 µs |  | 1.8x faster |
| `pchip/n10000/shuffled` | 130 µs | 200 µs | 250 µs | 164 µs |  | 1.3x faster |
| `akima/n10000/sorted` | 74.5 µs | 82.4 µs | 709 µs | 144 µs |  | 1.9x faster |
| `akima/n10000/shuffled` | 121 µs | 155 µs | 696 µs | 139 µs |  | 1.1x faster |
| `makima/n10000/sorted` | 85.2 µs | 209 µs | 831 µs |  |  | 9.8x faster |
| `makima/n10000/shuffled` | 83.8 µs | 86.8 µs | 767 µs |  |  | 9.1x faster |
| `linear/n100/sorted-uneven` | 252 ns | 247 ns | 9.75 µs | 203 ns | 185 ns | **1.4x slower** |
| `linear/n100/shuffled-uneven` | 259 ns | 266 ns | 9.21 µs | 208 ns | 184 ns | **1.4x slower** |
| `cubic_spline/n100/sorted-uneven` | 1.73 µs | 1.74 µs | 96.9 µs | 1.66 µs | 2.06 µs | **1x slower** |
| `cubic_spline/n100/shuffled-uneven` | 1.73 µs | 1.58 µs | 95.6 µs | 1.67 µs | 1.89 µs | **1x slower** |
| `pchip/n100/sorted-uneven` | 883 ns | 890 ns | 112 µs | 554 ns |  | **1.6x slower** |
| `pchip/n100/shuffled-uneven` | 966 ns | 967 ns | 111 µs | 581 ns |  | **1.7x slower** |
| `linear/n10000/sorted-uneven` | 94.1 µs | 93.6 µs | 48.8 µs | 59.2 µs | 23.1 µs | **4.1x slower** |
| `linear/n10000/shuffled-uneven` | 97.4 µs | 88 µs | 47.9 µs | 60 µs | 31.9 µs | **3.1x slower** |
| `cubic_spline/n10000/sorted-uneven` | 246 µs | 259 µs | 750 µs | 326 µs | 360 µs | 1.3x faster |
| `cubic_spline/n10000/shuffled-uneven` | 245 µs | 259 µs | 765 µs | 325 µs | 369 µs | 1.3x faster |
| `pchip/n10000/sorted-uneven` | 135 µs | 91.9 µs | 238 µs | 158 µs |  | 1.2x faster |
| `pchip/n10000/shuffled-uneven` | 123 µs | 212 µs | 239 µs | 166 µs |  | 1.3x faster |

## `PPoly`

A cubic of 55 intervals, evaluated at points in no order.

### Evaluation

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `intervals55/eval10` | 71.7 ns | 75.6 ns | 3.45 µs | 48x faster |
| `intervals55/eval100` | 659 ns | 754 ns | 4.34 µs | 6.6x faster |
| `intervals55/eval1000` | 6.54 µs | 7.28 µs | 15.7 µs | 2.4x faster |
| `intervals55/eval10000` | 72.8 µs | 67 µs | 360 µs | 5x faster |

### Construction

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `intervals55/eval10` | 160 ns | 158 ns | 5.15 µs | 32x faster |
| `intervals55/eval100` | 158 ns | 159 ns | 5.22 µs | 33x faster |
| `intervals55/eval1000` | 162 ns | 160 ns | 5.15 µs | 32x faster |
| `intervals55/eval10000` | 173 ns | 160 ns | 5.17 µs | 30x faster |

## `BSpline`

A cubic B-spline of 55 intervals, evaluated at points in no order.

### Evaluation

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `intervals55/eval10` | 146 ns | 147 ns | 2.67 µs | 18x faster |
| `intervals55/eval100` | 1.42 µs | 1.51 µs | 5.73 µs | 4x faster |
| `intervals55/eval1000` | 14.2 µs | 15.3 µs | 37.7 µs | 2.7x faster |
| `intervals55/eval10000` | 142 µs | 52.4 µs | 395 µs | 2.8x faster |

### Construction

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `intervals55/eval10` | 160 ns | 160 ns | 9.18 µs | 57x faster |
| `intervals55/eval100` | 173 ns | 161 ns | 9.27 µs | 54x faster |
| `intervals55/eval1000` | 161 ns | 161 ns | 9.22 µs | 57x faster |
| `intervals55/eval10000` | 160 ns | 174 ns | 10.1 µs | 63x faster |

## `BSplineCurve`

A cubic curve through `n` points of the plane, evaluated at `4 n` parameters in order.

### Evaluation

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `degree3/n100` | 4.14 µs | 4.47 µs | 15.9 µs | 3.8x faster |
| `degree3/n10000` | 447 µs | 152 µs | 1.33 ms | 3x faster |

### Construction

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `degree3/n100` | 4.46 µs | 4.25 µs | 46.6 µs | 10x faster |
| `degree3/n10000` | 638 µs | 689 µs | 1.37 ms | 2.2x faster |

## Regular grid

Points in no order inside the grid. `maxN`: the first axis has `N` points and each axis half the one before (SciPy's suite); `eachN`: every axis has `N` points; `flipped`: every axis given from its high end.

### Evaluation

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | ndarray-interp | interpn | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/2d/max10/eval10/ascending` | 135 ns | 149 ns | 136 ns | 19.2 µs | 782 ns | 151 ns | 95.1 ns | **1.4x slower** |
| `linear/2d/max10/eval10/flipped` | 135 ns | 136 ns | 148 ns | 19.2 µs |  |  |  | 1.4e+02x faster |
| `linear/2d/max10/eval100/ascending` | 1.25 µs | 1.31 µs | 1.25 µs | 21.5 µs | 8.21 µs | 1.38 µs | 1.02 µs | **1.2x slower** |
| `linear/2d/max10/eval100/flipped` | 1.25 µs | 1.31 µs | 1.25 µs | 21.2 µs |  |  |  | 17x faster |
| `linear/2d/max10/eval1000/ascending` | 13.7 µs | 12.6 µs | 12.5 µs | 40 µs | 90.1 µs | 14 µs | 9.27 µs | **1.5x slower** |
| `linear/2d/max10/eval1000/flipped` | 12.5 µs | 12.5 µs | 13.6 µs | 40 µs |  |  |  | 3.2x faster |
| `linear/2d/max10/eval10000/ascending` | 125 µs | 48.9 µs | 125 µs | 362 µs | 846 µs | 220 µs | 92.4 µs | **1.3x slower** |
| `linear/2d/max10/eval10000/flipped` | 124 µs | 49.6 µs | 124 µs | 360 µs |  |  |  | 2.9x faster |
| `linear/2d/max40/eval10/ascending` | 156 ns | 170 ns | 157 ns | 19.1 µs | 987 ns | 169 ns | 134 ns | **1.2x slower** |
| `linear/2d/max40/eval10/flipped` | 157 ns | 159 ns | 172 ns | 19.2 µs |  |  |  | 1.2e+02x faster |
| `linear/2d/max40/eval100/ascending` | 1.51 µs | 1.55 µs | 1.49 µs | 22.3 µs | 9.9 µs | 1.45 µs | 1.32 µs | **1.1x slower** |
| `linear/2d/max40/eval100/flipped` | 1.51 µs | 1.59 µs | 1.5 µs | 22 µs |  |  |  | 15x faster |
| `linear/2d/max40/eval1000/ascending` | 14.8 µs | 15.2 µs | 14.8 µs | 46.1 µs | 102 µs | 14.8 µs | 14.2 µs | **1x slower** |
| `linear/2d/max40/eval1000/flipped` | 15.1 µs | 16.7 µs | 15.1 µs | 46.2 µs |  |  |  | 3.1x faster |
| `linear/2d/max40/eval10000/ascending` | 151 µs | 55.6 µs | 151 µs | 574 µs | 1.08 ms | 225 µs | 132 µs | **1.1x slower** |
| `linear/2d/max40/eval10000/flipped` | 166 µs | 54 µs | 151 µs | 575 µs |  |  |  | 3.5x faster |
| `linear/2d/max200/eval10/ascending` | 193 ns | 193 ns | 193 ns | 19.2 µs | 1.92 µs | 205 ns | 175 ns | **1.1x slower** |
| `linear/2d/max200/eval10/flipped` | 191 ns | 194 ns | 194 ns | 19.2 µs |  |  |  | 1e+02x faster |
| `linear/2d/max200/eval100/ascending` | 2.05 µs | 1.96 µs | 2.05 µs | 24.9 µs | 19.4 µs | 1.99 µs | 1.72 µs | **1.2x slower** |
| `linear/2d/max200/eval100/flipped` | 1.87 µs | 2.14 µs | 2.06 µs | 22.7 µs |  |  |  | 12x faster |
| `linear/2d/max200/eval1000/ascending` | 18.8 µs | 19.4 µs | 20.5 µs | 72.6 µs | 226 µs | 22.3 µs | 17.2 µs | **1.1x slower** |
| `linear/2d/max200/eval1000/flipped` | 19.6 µs | 19 µs | 18.8 µs | 72.8 µs |  |  |  | 3.7x faster |
| `linear/2d/max200/eval10000/ascending` | 206 µs | 71.6 µs | 188 µs | 842 µs | 2.17 ms | 453 µs | 172 µs | **1.2x slower** |
| `linear/2d/max200/eval10000/flipped` | 188 µs | 74.1 µs | 188 µs | 826 µs |  |  |  | 4.4x faster |
| `linear/3d/max10/eval10/ascending` | 215 ns | 205 ns | 214 ns | 55 µs | 3.15 µs |  | 143 ns | **1.5x slower** |
| `linear/3d/max10/eval10/flipped` | 198 ns | 206 ns | 196 ns | 56.3 µs |  |  |  | 2.8e+02x faster |
| `linear/3d/max10/eval100/ascending` | 1.88 µs | 1.97 µs | 1.9 µs | 61 µs | 30.8 µs |  | 1.39 µs | **1.3x slower** |
| `linear/3d/max10/eval100/flipped` | 1.9 µs | 1.97 µs | 1.91 µs | 62.4 µs |  |  |  | 33x faster |
| `linear/3d/max10/eval1000/ascending` | 19 µs | 19 µs | 18.7 µs | 116 µs | 308 µs |  | 14.1 µs | **1.4x slower** |
| `linear/3d/max10/eval1000/flipped` | 18.9 µs | 19.1 µs | 19 µs | 116 µs |  |  |  | 6.1x faster |
| `linear/3d/max10/eval10000/ascending` | 192 µs | 66.2 µs | 190 µs | 1.08 ms | 3.44 ms |  | 142 µs | **1.4x slower** |
| `linear/3d/max10/eval10000/flipped` | 190 µs | 67 µs | 208 µs | 1.12 ms |  |  |  | 5.9x faster |
| `linear/3d/max40/eval10/ascending` | 218 ns | 226 ns | 219 ns | 56.2 µs | 3.36 µs |  | 199 ns | **1.1x slower** |
| `linear/3d/max40/eval10/flipped` | 221 ns | 225 ns | 220 ns | 56.4 µs |  |  |  | 2.5e+02x faster |
| `linear/3d/max40/eval100/ascending` | 2.12 µs | 2.25 µs | 2.12 µs | 62.1 µs | 33.7 µs |  | 1.97 µs | **1.1x slower** |
| `linear/3d/max40/eval100/flipped` | 2.08 µs | 2.26 µs | 2.12 µs | 63.1 µs |  |  |  | 30x faster |
| `linear/3d/max40/eval1000/ascending` | 21.2 µs | 21.4 µs | 20.9 µs | 139 µs | 339 µs |  | 19.8 µs | **1.1x slower** |
| `linear/3d/max40/eval1000/flipped` | 21.3 µs | 21.5 µs | 21.9 µs | 137 µs |  |  |  | 6.4x faster |
| `linear/3d/max40/eval10000/ascending` | 210 µs | 71.2 µs | 212 µs | 1.46 ms | 3.45 ms |  | 200 µs | **1.1x slower** |
| `linear/3d/max40/eval10000/flipped` | 229 µs | 69.9 µs | 212 µs | 1.46 ms |  |  |  | 6.4x faster |
| `linear/3d/max200/eval10/ascending` | 301 ns | 281 ns | 302 ns | 56.4 µs | 4.56 µs |  | 261 ns | **1.2x slower** |
| `linear/3d/max200/eval10/flipped` | 274 ns | 308 ns | 277 ns | 56.6 µs |  |  |  | 2.1e+02x faster |
| `linear/3d/max200/eval100/ascending` | 2.78 µs | 3.05 µs | 2.97 µs | 65.7 µs | 47.1 µs |  | 2.61 µs | **1.1x slower** |
| `linear/3d/max200/eval100/flipped` | 2.95 µs | 3.05 µs | 2.97 µs | 64.9 µs |  |  |  | 22x faster |
| `linear/3d/max200/eval1000/ascending` | 27.1 µs | 30.1 µs | 29.7 µs | 194 µs | 495 µs |  | 27.8 µs | 1x faster |
| `linear/3d/max200/eval1000/flipped` | 29.7 µs | 29.9 µs | 27 µs | 189 µs |  |  |  | 6.4x faster |
| `linear/3d/max200/eval10000/ascending` | 270 µs | 79.1 µs | 270 µs | 1.83 ms | 4.85 ms |  | 264 µs | **1x slower** |
| `linear/3d/max200/eval10000/flipped` | 296 µs | 80.8 µs | 274 µs | 1.87 ms |  |  |  | 6.3x faster |
| `linear/4d/max10/eval10/ascending` | 261 ns | 285 ns | 261 ns | 98.3 µs |  |  |  | 3.8e+02x faster |
| `linear/4d/max10/eval10/flipped` | 260 ns | 281 ns | 261 ns | 104 µs |  |  |  | 4e+02x faster |
| `linear/4d/max10/eval100/ascending` | 2.42 µs | 2.47 µs | 2.62 µs | 115 µs |  |  |  | 48x faster |
| `linear/4d/max10/eval100/flipped` | 2.41 µs | 2.48 µs | 2.41 µs | 106 µs |  |  |  | 44x faster |
| `linear/4d/max10/eval1000/ascending` | 24.5 µs | 26.9 µs | 24.6 µs | 213 µs |  |  |  | 8.7x faster |
| `linear/4d/max10/eval1000/flipped` | 25.2 µs | 24.7 µs | 24.3 µs | 212 µs |  |  |  | 8.4x faster |
| `linear/4d/max10/eval10000/ascending` | 249 µs | 86.7 µs | 266 µs | 1.92 ms |  |  |  | 7.7x faster |
| `linear/4d/max10/eval10000/flipped` | 243 µs | 82 µs | 265 µs | 1.89 ms |  |  |  | 7.8x faster |
| `linear/4d/max40/eval10/ascending` | 340 ns | 321 ns | 310 ns | 95.6 µs | 5.75 µs |  | 270 ns | **1.3x slower** |
| `linear/4d/max40/eval10/flipped` | 305 ns | 314 ns | 314 ns | 97 µs |  |  |  | 3.2e+02x faster |
| `linear/4d/max40/eval100/ascending` | 3.24 µs | 3.13 µs | 3.2 µs | 108 µs | 59.9 µs |  | 2.95 µs | **1.1x slower** |
| `linear/4d/max40/eval100/flipped` | 3.25 µs | 3.32 µs | 3.21 µs | 117 µs |  |  |  | 36x faster |
| `linear/4d/max40/eval1000/ascending` | 32.5 µs | 32.3 µs | 32.4 µs | 274 µs | 639 µs |  | 28.5 µs | **1.1x slower** |
| `linear/4d/max40/eval1000/flipped` | 29.4 µs | 30.8 µs | 32.2 µs | 244 µs |  |  |  | 8.3x faster |
| `linear/4d/max40/eval10000/ascending` | 294 µs | 102 µs | 323 µs | 2.32 ms | 7.22 ms |  | 297 µs | 1x faster |
| `linear/4d/max40/eval10000/flipped` | 293 µs | 101 µs | 319 µs | 2.42 ms |  |  |  | 8.3x faster |
| `linear/4d/max200/eval10/ascending` | 412 ns | 395 ns | 386 ns | 96.2 µs | 7.79 µs |  | 356 ns | **1.2x slower** |
| `linear/4d/max200/eval10/flipped` | 379 ns | 389 ns | 386 ns | 95.2 µs |  |  |  | 2.5e+02x faster |
| `linear/4d/max200/eval100/ascending` | 3.73 µs | 3.83 µs | 3.74 µs | 108 µs | 72 µs |  | 3.6 µs | **1x slower** |
| `linear/4d/max200/eval100/flipped` | 3.72 µs | 3.86 µs | 3.75 µs | 116 µs |  |  |  | 31x faster |
| `linear/4d/max200/eval1000/ascending` | 37.5 µs | 39.8 µs | 37.6 µs | 325 µs | 710 µs |  | 36.2 µs | **1x slower** |
| `linear/4d/max200/eval1000/flipped` | 37.5 µs | 39 µs | 37.4 µs | 311 µs |  |  |  | 8.3x faster |
| `linear/4d/max200/eval10000/ascending` | 411 µs | 114 µs | 373 µs | 2.83 ms | 7.51 ms |  | 365 µs | **1.1x slower** |
| `linear/4d/max200/eval10000/flipped` | 409 µs | 118 µs | 373 µs | 2.75 ms |  |  |  | 6.7x faster |
| `nearest/2d/each40/eval10000/ascending` | 150 µs | 63.7 µs | 150 µs | 653 µs | 1.58 ms |  | 136 µs | **1.1x slower** |
| `nearest/3d/each40/eval10000/ascending` | 212 µs | 75.6 µs | 232 µs | 928 µs | 2.1 ms |  | 227 µs | 1.1x faster |
| `pchip/2d/each40/eval10000/ascending` | 499 µs | 135 µs | 500 µs | 1.2 s |  |  |  | 2.4e+03x faster |
| `pchip/3d/each40/eval10000/ascending` | 1.51 ms | 284 µs | 1.51 ms | 3.28 s |  |  |  | 2.2e+03x faster |
| `cubic/2d/each40/eval10000/ascending` | 330 µs | 97.6 µs | 330 µs | 1.35 ms |  |  | 418 µs | 1.3x faster |
| `cubic/3d/each40/eval10000/ascending` | 727 µs | 177 µs | 722 µs | 3.03 ms |  |  | 859 µs | 1.2x faster |
| `cubic/2d/each10/eval1000/ascending` | 26.6 µs | 27.7 µs | 26.6 µs | 131 µs |  |  | 39.2 µs | 1.5x faster |
| `cubic/2d/each40/eval1000/ascending` | 32.9 µs | 33 µs | 30 µs | 154 µs |  |  | 44.7 µs | 1.4x faster |
| `cubic/2d/each100/eval1000/ascending` | 35.2 µs | 35.4 µs | 32.1 µs | 181 µs |  |  | 47.1 µs | 1.3x faster |
| `cubic/2d/each200/eval1000/ascending` | 35.4 µs | 38.1 µs | 38 µs | 203 µs |  |  | 46.1 µs | 1.3x faster |
| `cubic/2d/each400/eval1000/ascending` | 41.6 µs | 39.9 µs | 42 µs | 297 µs |  |  | 49.5 µs | 1.2x faster |
| `quintic/2d/each10/eval1000/ascending` | 44.6 µs | 48.5 µs | 46 µs | 239 µs |  |  |  | 5.4x faster |
| `quintic/2d/each40/eval1000/ascending` | 48.6 µs | 48.3 µs | 48.8 µs | 238 µs |  |  |  | 4.9x faster |
| `linear/2d/each400/eval10000/ascending` | 218 µs | 83.8 µs | 239 µs | 970 µs | 4.28 ms | 538 µs | 219 µs | 1x faster |
| `linear/2d/each800/eval10000/ascending` | 238 µs | 84.9 µs | 238 µs | 1.08 ms | 7.97 ms | 578 µs | 242 µs | 1x faster |

### Construction

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | ndarray-interp | interpn | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/2d/max10/eval10/ascending` | 134 ns | 137 ns | 134 ns | 5.8 µs | 132 ns | 77.8 ns | 0 | **2e+02x slower** |
| `linear/2d/max10/eval10/flipped` | 153 ns | 158 ns | 166 ns | 11.8 µs |  |  |  | 77x faster |
| `linear/2d/max10/eval100/ascending` | 124 ns | 124 ns | 130 ns | 5.82 µs | 133 ns | 77.8 ns | 0 | **1.8e+02x slower** |
| `linear/2d/max10/eval100/flipped` | 153 ns | 153 ns | 169 ns | 11.9 µs |  |  |  | 78x faster |
| `linear/2d/max10/eval1000/ascending` | 135 ns | 124 ns | 131 ns | 5.87 µs | 146 ns | 85.4 ns | 0 | **2e+02x slower** |
| `linear/2d/max10/eval1000/flipped` | 152 ns | 168 ns | 164 ns | 11.7 µs |  |  |  | 77x faster |
| `linear/2d/max10/eval10000/ascending` | 122 ns | 130 ns | 136 ns | 5.84 µs | 132 ns | 78 ns | 0 | **1.8e+02x slower** |
| `linear/2d/max10/eval10000/flipped` | 152 ns | 156 ns | 152 ns | 11.6 µs |  |  |  | 76x faster |
| `linear/2d/max40/eval10/ascending` | 360 ns | 363 ns | 741 ns | 5.73 µs | 214 ns | 164 ns | 0 | **5.3e+02x slower** |
| `linear/2d/max40/eval10/flipped` | 641 ns | 639 ns | 1.21 µs | 11.7 µs |  |  |  | 18x faster |
| `linear/2d/max40/eval100/ascending` | 356 ns | 356 ns | 619 ns | 5.86 µs | 209 ns | 150 ns | 0 | **5.2e+02x slower** |
| `linear/2d/max40/eval100/flipped` | 638 ns | 638 ns | 1.19 µs | 11.6 µs |  |  |  | 18x faster |
| `linear/2d/max40/eval1000/ascending` | 357 ns | 357 ns | 625 ns | 5.88 µs | 212 ns | 147 ns | 0 | **5.2e+02x slower** |
| `linear/2d/max40/eval1000/flipped` | 639 ns | 639 ns | 1.18 µs | 11.4 µs |  |  |  | 18x faster |
| `linear/2d/max40/eval10000/ascending` | 358 ns | 356 ns | 588 ns | 5.89 µs | 215 ns | 150 ns | 0 | **5.2e+02x slower** |
| `linear/2d/max40/eval10000/flipped` | 697 ns | 637 ns | 1.2 µs | 11.5 µs |  |  |  | 17x faster |
| `linear/2d/max200/eval10/ascending` | 7.47 µs | 6.92 µs | 15.6 µs | 5.82 µs | 2.96 µs | 470 ns | 0 | **1.1e+04x slower** |
| `linear/2d/max200/eval10/flipped` | 14.5 µs | 14.6 µs | 26.8 µs | 11.9 µs |  |  |  | **1.2x slower** |
| `linear/2d/max200/eval100/ascending` | 7.49 µs | 7.37 µs | 16.5 µs | 6.37 µs | 3.09 µs | 475 ns | 0 | **1.1e+04x slower** |
| `linear/2d/max200/eval100/flipped` | 14.5 µs | 15.1 µs | 26.5 µs | 11.9 µs |  |  |  | **1.2x slower** |
| `linear/2d/max200/eval1000/ascending` | 6.85 µs | 6.88 µs | 18.8 µs | 5.9 µs | 2.97 µs | 515 ns | 0 | **1e+04x slower** |
| `linear/2d/max200/eval1000/flipped` | 14.8 µs | 14.5 µs | 26 µs | 11.9 µs |  |  |  | **1.2x slower** |
| `linear/2d/max200/eval10000/ascending` | 47.2 µs | 47.4 µs | 17.5 µs | 5.88 µs | 53.3 µs | 470 ns | 0 | **6.9e+04x slower** |
| `linear/2d/max200/eval10000/flipped` | 51.7 µs | 54.7 µs | 29.3 µs | 11.8 µs |  |  |  | **4.4x slower** |
| `linear/3d/max10/eval10/ascending` | 193 ns | 175 ns | 215 ns | 7.42 µs | 166 ns |  | 0 | **2.2e+02x slower** |
| `linear/3d/max10/eval10/flipped` | 295 ns | 295 ns | 336 ns | 15.7 µs |  |  |  | 53x faster |
| `linear/3d/max10/eval100/ascending` | 191 ns | 176 ns | 203 ns | 7.51 µs | 164 ns |  | 0 | **2.2e+02x slower** |
| `linear/3d/max10/eval100/flipped` | 294 ns | 300 ns | 332 ns | 15.9 µs |  |  |  | 54x faster |
| `linear/3d/max10/eval1000/ascending` | 177 ns | 178 ns | 203 ns | 7.6 µs | 166 ns |  | 0 | **2e+02x slower** |
| `linear/3d/max10/eval1000/flipped` | 295 ns | 298 ns | 320 ns | 15.6 µs |  |  |  | 53x faster |
| `linear/3d/max10/eval10000/ascending` | 192 ns | 176 ns | 185 ns | 7.62 µs | 179 ns |  | 0 | **2e+02x slower** |
| `linear/3d/max10/eval10000/flipped` | 297 ns | 294 ns | 338 ns | 15.8 µs |  |  |  | 53x faster |
| `linear/3d/max40/eval10/ascending` | 2.79 µs | 15.1 µs | 7.21 µs | 7.58 µs | 3.45 µs |  | 0 | **3.2e+03x slower** |
| `linear/3d/max40/eval10/flipped` | 8.08 µs | 8.15 µs | 13.4 µs | 17 µs |  |  |  | 2.1x faster |
| `linear/3d/max40/eval100/ascending` | 2.78 µs | 15.5 µs | 5.59 µs | 7.62 µs | 8.19 µs |  | 0 | **3.2e+03x slower** |
| `linear/3d/max40/eval100/flipped` | 8.33 µs | 8.92 µs | 12.5 µs | 15.6 µs |  |  |  | 1.9x faster |
| `linear/3d/max40/eval1000/ascending` | 15.2 µs | 12.7 µs | 5.77 µs | 7.5 µs | 8.89 µs |  | 0 | **1.7e+04x slower** |
| `linear/3d/max40/eval1000/flipped` | 27.7 µs | 24.9 µs | 12.7 µs | 15.6 µs |  |  |  | **1.8x slower** |
| `linear/3d/max40/eval10000/ascending` | 16.3 µs | 15.2 µs | 6.06 µs | 7.48 µs | 12.5 µs |  | 0 | **1.9e+04x slower** |
| `linear/3d/max40/eval10000/flipped` | 26.1 µs | 24.8 µs | 12.1 µs | 15.9 µs |  |  |  | **1.6x slower** |
| `linear/3d/max200/eval10/ascending` | 1.73 ms | 1.76 ms | 583 µs | 7.51 µs | 1.26 ms |  | 0 | **2e+06x slower** |
| `linear/3d/max200/eval10/flipped` | 3.04 ms | 3.21 ms | 1.94 ms | 17 µs |  |  |  | **1.8e+02x slower** |
| `linear/3d/max200/eval100/ascending` | 1.71 ms | 1.74 ms | 641 µs | 7.64 µs | 1.51 ms |  | 0 | **1.9e+06x slower** |
| `linear/3d/max200/eval100/flipped` | 3.22 ms | 3.09 ms | 1.89 ms | 16 µs |  |  |  | **2e+02x slower** |
| `linear/3d/max200/eval1000/ascending` | 1.69 ms | 1.78 ms | 645 µs | 8.3 µs | 1.3 ms |  | 0 | **1.9e+06x slower** |
| `linear/3d/max200/eval1000/flipped` | 3.21 ms | 3.11 ms | 1.97 ms | 15.9 µs |  |  |  | **2e+02x slower** |
| `linear/3d/max200/eval10000/ascending` | 1.75 ms | 1.73 ms | 619 µs | 8.19 µs | 1.2 ms |  | 0 | **2e+06x slower** |
| `linear/3d/max200/eval10000/flipped` | 3.2 ms | 3.19 ms | 1.91 ms | 16 µs |  |  |  | **2e+02x slower** |
| `linear/4d/max10/eval10/ascending` | 228 ns | 210 ns | 247 ns | 10.5 µs |  |  |  | 46x faster |
| `linear/4d/max10/eval10/flipped` | 361 ns | 351 ns | 362 ns | 18.5 µs |  |  |  | 51x faster |
| `linear/4d/max10/eval100/ascending` | 210 ns | 210 ns | 235 ns | 9.78 µs |  |  |  | 47x faster |
| `linear/4d/max10/eval100/flipped` | 333 ns | 364 ns | 363 ns | 17.9 µs |  |  |  | 54x faster |
| `linear/4d/max10/eval1000/ascending` | 229 ns | 213 ns | 235 ns | 9.73 µs |  |  |  | 43x faster |
| `linear/4d/max10/eval1000/flipped` | 359 ns | 361 ns | 372 ns | 18 µs |  |  |  | 50x faster |
| `linear/4d/max10/eval10000/ascending` | 228 ns | 230 ns | 260 ns | 9.74 µs |  |  |  | 43x faster |
| `linear/4d/max10/eval10000/flipped` | 355 ns | 375 ns | 363 ns | 17.8 µs |  |  |  | 50x faster |
| `linear/4d/max40/eval10/ascending` | 14.4 µs | 14.3 µs | 30.5 µs | 10.2 µs | 5.86 µs |  | 1.27 ns | **1.1e+04x slower** |
| `linear/4d/max40/eval10/flipped` | 51.4 µs | 52.2 µs | 67.1 µs | 19.8 µs |  |  |  | **2.6x slower** |
| `linear/4d/max40/eval100/ascending` | 14.3 µs | 14.8 µs | 37.3 µs | 9.51 µs | 6.02 µs |  | 1.46 ns | **9.8e+03x slower** |
| `linear/4d/max40/eval100/flipped` | 57 µs | 56.9 µs | 68.6 µs | 22.2 µs |  |  |  | **2.6x slower** |
| `linear/4d/max40/eval1000/ascending` | 14.4 µs | 14.4 µs | 30.8 µs | 10.2 µs | 6.05 µs |  | 1.27 ns | **1.1e+04x slower** |
| `linear/4d/max40/eval1000/flipped` | 52.9 µs | 56.4 µs | 63.6 µs | 19.8 µs |  |  |  | **2.7x slower** |
| `linear/4d/max40/eval10000/ascending` | 13.1 µs | 95.8 µs | 31.8 µs | 9.47 µs | 5.99 µs |  | 1.37 ns | **9.6e+03x slower** |
| `linear/4d/max40/eval10000/flipped` | 52.7 µs | 57.1 µs | 72.1 µs | 21.6 µs |  |  |  | **2.4x slower** |
| `linear/4d/max200/eval10/ascending` | 46.3 ms | 48.9 ms | 14.8 ms | 9.51 µs | 39.1 ms |  | 1.17 ns | **4e+07x slower** |
| `linear/4d/max200/eval10/flipped` | 105 ms | 101 ms | 66.7 ms | 20.1 µs |  |  |  | **5.2e+03x slower** |
| `linear/4d/max200/eval100/ascending` | 48.3 ms | 48.6 ms | 15.5 ms | 9.43 µs | 42.2 ms |  | 1.27 ns | **3.8e+07x slower** |
| `linear/4d/max200/eval100/flipped` | 101 ms | 105 ms | 66.6 ms | 20.1 µs |  |  |  | **5e+03x slower** |
| `linear/4d/max200/eval1000/ascending` | 50 ms | 49.7 ms | 15.5 ms | 9.47 µs | 39.2 ms |  | 1.17 ns | **4.3e+07x slower** |
| `linear/4d/max200/eval1000/flipped` | 103 ms | 104 ms | 68.9 ms | 20.4 µs |  |  |  | **5e+03x slower** |
| `linear/4d/max200/eval10000/ascending` | 48.2 ms | 47.9 ms | 15.6 ms | 9.28 µs | 39 ms |  | 1.27 ns | **3.8e+07x slower** |
| `linear/4d/max200/eval10000/flipped` | 98.6 ms | 103 ms | 70 ms | 20 µs |  |  |  | **4.9e+03x slower** |
| `nearest/2d/each40/eval10000/ascending` | 586 ns | 638 ns | 1.21 µs | 5.84 µs | 283 ns |  | 0 | **7.5e+02x slower** |
| `nearest/3d/each40/eval10000/ascending` | 22.6 µs | 20.9 µs | 53.9 µs | 7.65 µs | 9.38 µs |  | 0 | **2.3e+04x slower** |
| `pchip/2d/each40/eval10000/ascending` | 7.05 µs | 6.54 µs | 6.81 µs | 6.7 µs |  |  |  | **1.1x slower** |
| `pchip/3d/each40/eval10000/ascending` | 509 µs | 633 µs | 492 µs | 8.55 µs |  |  |  | **60x slower** |
| `cubic/2d/each40/eval10000/ascending` | 9.05 µs | 9.3 µs | 8.98 µs | 1.91 ms |  |  | 0 | **1e+04x slower** |
| `cubic/3d/each40/eval10000/ascending` | 423 µs | 461 µs | 334 µs | 144 ms |  |  | 0 | **4.3e+05x slower** |
| `cubic/2d/each10/eval1000/ascending` | 1.62 µs | 1.61 µs | 1.68 µs | 995 µs |  |  | 0 | **1.7e+03x slower** |
| `cubic/2d/each40/eval1000/ascending` | 9.12 µs | 9.25 µs | 8.78 µs | 1.92 ms |  |  | 0 | **1e+04x slower** |
| `cubic/2d/each100/eval1000/ascending` | 84 µs | 112 µs | 83.5 µs | 7.41 ms |  |  | 0 | **8.6e+04x slower** |
| `cubic/2d/each200/eval1000/ascending` | 179 µs | 292 µs | 239 µs | 29.8 ms |  |  | 0 | **1.8e+05x slower** |
| `cubic/2d/each400/eval1000/ascending` | 1.42 ms | 1.46 ms | 1.21 ms | 146 ms |  |  | 0 | **1.8e+06x slower** |
| `quintic/2d/each10/eval1000/ascending` | 2.12 µs | 2.12 µs | 2.1 µs | 3.51 ms |  |  |  | 1.7e+03x faster |
| `quintic/2d/each40/eval1000/ascending` | 13.2 µs | 13.2 µs | 12 µs | 7.88 ms |  |  |  | 6e+02x faster |
| `linear/2d/each400/eval10000/ascending` | 312 µs | 354 µs | 101 µs | 6.03 µs | 284 µs | 1.16 µs | 0 | **4.6e+05x slower** |
| `linear/2d/each800/eval10000/ascending` | 1.14 ms | 1.15 ms | 417 µs | 6.34 µs | 805 µs | 2.22 µs | 0 | **1.7e+06x slower** |

## `NearestNDInterpolator`

`n` points of the unit square, evaluated at the pixels of a square image.

### Evaluation

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: |
| `n10/image50` | 119 µs | 88.2 µs | 114 µs | 441 µs | 3.7x faster |
| `n10/image100` | 500 µs | 146 µs | 495 µs | 1.76 ms | 3.5x faster |
| `n10/image500` | 12.9 ms | 2.12 ms | 12.3 ms | 46.5 ms | 3.6x faster |
| `n100/image50` | 186 µs | 84.1 µs | 186 µs | 578 µs | 3.1x faster |
| `n100/image100` | 777 µs | 208 µs | 783 µs | 2.25 ms | 2.9x faster |
| `n100/image500` | 18.3 ms | 3.07 ms | 18 ms | 59 ms | 3.2x faster |
| `n1000/image50` | 333 µs | 121 µs | 308 µs | 796 µs | 2.4x faster |
| `n1000/image100` | 1.14 ms | 292 µs | 1.15 ms | 3.13 ms | 2.7x faster |
| `n1000/image500` | 24.2 ms | 4.13 ms | 24.1 ms | 68.8 ms | 2.8x faster |
| `n10000/image50` | 350 µs | 119 µs | 347 µs | 1.07 ms | 3.1x faster |
| `n10000/image100` | 1.34 ms | 281 µs | 1.46 ms | 4.04 ms | 3x faster |
| `n10000/image500` | 30.7 ms | 5.06 ms | 30.7 ms | 85.9 ms | 2.8x faster |
| `n100000/image50` | 627 µs | 117 µs | 472 µs | 1.44 ms | 2.3x faster |
| `n100000/image100` | 1.86 ms | 364 µs | 1.87 ms | 5.66 ms | 3x faster |
| `n100000/image500` | 43 ms | 6.85 ms | 44 ms | 125 ms | 2.9x faster |

### Construction

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: |
| `n10/image50` | 566 ns | 569 ns | 514 ns | 10.5 µs | 19x faster |
| `n10/image100` | 557 ns | 556 ns | 549 ns | 10.7 µs | 19x faster |
| `n10/image500` | 557 ns | 566 ns | 558 ns | 9.55 µs | 17x faster |
| `n100/image50` | 3.59 µs | 3.9 µs | 3.67 µs | 16.7 µs | 4.6x faster |
| `n100/image100` | 4.25 µs | 3.63 µs | 3.73 µs | 15.8 µs | 3.7x faster |
| `n100/image500` | 3.64 µs | 3.65 µs | 3.57 µs | 16 µs | 4.4x faster |
| `n1000/image50` | 59.7 µs | 59.1 µs | 51.7 µs | 128 µs | 2.1x faster |
| `n1000/image100` | 58.2 µs | 53.9 µs | 42.2 µs | 127 µs | 2.2x faster |
| `n1000/image500` | 60.2 µs | 63.7 µs | 50.7 µs | 124 µs | 2.1x faster |
| `n10000/image50` | 750 µs | 671 µs | 666 µs | 1.8 ms | 2.4x faster |
| `n10000/image100` | 657 µs | 686 µs | 676 µs | 1.67 ms | 2.5x faster |
| `n10000/image500` | 661 µs | 645 µs | 660 µs | 1.68 ms | 2.5x faster |
| `n100000/image50` | 9.54 ms | 8.36 ms | 8.35 ms | 25.8 ms | 2.7x faster |
| `n100000/image100` | 8.48 ms | 8.43 ms | 8.4 ms | 28.3 ms | 3.3x faster |
| `n100000/image500` | 8.26 ms | 9.28 ms | 8.07 ms | 26.1 ms | 3.2x faster |

## `RBFInterpolator`

`n` points of a square, evaluated at 1000 points in no order (`scattered`) or at the 10 000 pixels of an image (`image`). `global`: all points; `k50`: the 50 nearest.

### Evaluation

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/global/n10/scattered` | 14.8 µs | 16 µs | 15.9 µs | 63.2 µs | 183 µs | 4.3x faster |
| `linear/k50/n10/scattered` | 14.5 µs | 15.9 µs | 15.8 µs | 1.41 ms |  | 97x faster |
| `thin_plate_spline/global/n10/scattered` | 78.7 µs | 71.9 µs | 78.6 µs | 122 µs | 326 µs | 1.5x faster |
| `thin_plate_spline/k50/n10/scattered` | 78.5 µs | 78.7 µs | 78.5 µs | 1.6 ms |  | 20x faster |
| `cubic/global/n10/scattered` | 21.3 µs | 20.5 µs | 20.4 µs | 60 µs | 217 µs | 2.8x faster |
| `cubic/k50/n10/scattered` | 19.6 µs | 20 µs | 21.3 µs | 1.4 ms |  | 71x faster |
| `quintic/global/n10/scattered` | 28 µs | 29 µs | 28.1 µs | 73.8 µs | 229 µs | 2.6x faster |
| `quintic/k50/n10/scattered` | 30.6 µs | 29.8 µs | 30.7 µs | 1.46 ms |  | 48x faster |
| `multiquadric/global/n10/scattered` | 30.5 µs | 27.9 µs | 27.8 µs | 52.9 µs | 257 µs | 1.7x faster |
| `multiquadric/k50/n10/scattered` | 30.5 µs | 30.6 µs | 30.5 µs | 1.35 ms |  | 44x faster |
| `inverse_multiquadric/global/n10/scattered` | 38.5 µs | 37.6 µs | 37.4 µs | 55.2 µs | 293 µs | 1.4x faster |
| `inverse_multiquadric/k50/n10/scattered` | 38.6 µs | 41.2 µs | 41.1 µs | 1.52 ms |  | 39x faster |
| `inverse_quadratic/global/n10/scattered` | 23.8 µs | 26.3 µs | 23.8 µs | 48 µs |  | 2x faster |
| `inverse_quadratic/k50/n10/scattered` | 23.8 µs | 26.3 µs | 23.8 µs | 1.39 ms |  | 58x faster |
| `gaussian/global/n10/scattered` | 67.3 µs | 68.4 µs | 68.9 µs | 87.7 µs | 357 µs | 1.3x faster |
| `gaussian/k50/n10/scattered` | 67.2 µs | 74 µs | 69 µs | 1.41 ms |  | 21x faster |
| `linear/global/n100/scattered` | 135 µs | 49.3 µs | 148 µs | 288 µs | 1.88 ms | 2.1x faster |
| `linear/k50/n100/scattered` | 13 ms | 2.63 ms | 13.6 ms | 19.1 ms |  | 1.5x faster |
| `thin_plate_spline/global/n100/scattered` | 769 µs | 190 µs | 700 µs | 880 µs | 3.14 ms | 1.1x faster |
| `thin_plate_spline/k50/n100/scattered` | 24.4 ms | 5.07 ms | 26.9 ms | 25 ms |  | 1x faster |
| `cubic/global/n100/scattered` | 152 µs | 69.4 µs | 153 µs | 321 µs | 2.07 ms | 2.1x faster |
| `cubic/k50/n100/scattered` | 15.3 ms | 2.93 ms | 15.4 ms | 20.6 ms |  | 1.3x faster |
| `quintic/global/n100/scattered` | 167 µs | 64.3 µs | 167 µs | 303 µs | 2.27 ms | 1.8x faster |
| `quintic/k50/n100/scattered` | 16.8 ms | 3.19 ms | 17.2 ms | 22.7 ms |  | 1.3x faster |
| `multiquadric/global/n100/scattered` | 269 µs | 89.6 µs | 269 µs | 305 µs | 2.41 ms | 1.1x faster |
| `multiquadric/k50/n100/scattered` | 16 ms | 3.01 ms | 15.5 ms | 18.9 ms |  | 1.2x faster |
| `inverse_multiquadric/global/n100/scattered` | 378 µs | 110 µs | 364 µs | 395 µs | 2.84 ms | 1x faster |
| `inverse_multiquadric/k50/n100/scattered` | 15.7 ms | 2.94 ms | 15.7 ms | 19.9 ms |  | 1.3x faster |
| `inverse_quadratic/global/n100/scattered` | 246 µs | 83.7 µs | 232 µs | 288 µs |  | 1.2x faster |
| `inverse_quadratic/k50/n100/scattered` | 13.8 ms | 2.7 ms | 13.6 ms | 18.1 ms |  | 1.3x faster |
| `gaussian/global/n100/scattered` | 615 µs | 164 µs | 615 µs | 684 µs | 3.8 ms | 1.1x faster |
| `gaussian/k50/n100/scattered` | 20.6 ms | 3.89 ms | 20.6 ms | 21.7 ms |  | 1.1x faster |
| `linear/global/n1000/scattered` | 1.35 ms | 333 µs | 1.48 ms | 3.99 ms | 19.1 ms | 3x faster |
| `linear/k50/n1000/scattered` | 13.6 ms | 2.76 ms | 15.5 ms | 36.3 ms |  | 2.7x faster |
| `thin_plate_spline/global/n1000/scattered` | 7.09 ms | 1.55 ms | 7.22 ms | 10.4 ms | 29.4 ms | 1.5x faster |
| `thin_plate_spline/k50/n1000/scattered` | 22.8 ms | 4.34 ms | 23.3 ms | 51.7 ms |  | 2.3x faster |
| `cubic/global/n1000/scattered` | 1.38 ms | 334 µs | 1.35 ms | 4.05 ms | 19 ms | 2.9x faster |
| `cubic/k50/n1000/scattered` | 16 ms | 3.12 ms | 15.3 ms | 39.6 ms |  | 2.5x faster |
| `quintic/global/n1000/scattered` | 1.53 ms | 372 µs | 1.52 ms | 4.01 ms | 21.1 ms | 2.6x faster |
| `quintic/k50/n1000/scattered` | 18.3 ms | 3.28 ms | 18.1 ms | 42.8 ms |  | 2.3x faster |
| `multiquadric/global/n1000/scattered` | 2.86 ms | 605 µs | 2.74 ms | 4.18 ms | 24.5 ms | 1.5x faster |
| `multiquadric/k50/n1000/scattered` | 17.5 ms | 3.28 ms | 17.2 ms | 35.9 ms |  | 2x faster |
| `inverse_multiquadric/global/n1000/scattered` | 4 ms | 804 µs | 3.77 ms | 5.17 ms | 27.5 ms | 1.3x faster |
| `inverse_multiquadric/k50/n1000/scattered` | 17.4 ms | 3.3 ms | 17.5 ms | 37 ms |  | 2.1x faster |
| `inverse_quadratic/global/n1000/scattered` | 2.43 ms | 531 µs | 2.27 ms | 3.75 ms |  | 1.5x faster |
| `inverse_quadratic/k50/n1000/scattered` | 16.5 ms | 2.97 ms | 16 ms | 35 ms |  | 2.1x faster |
| `gaussian/global/n1000/scattered` | 6.58 ms | 1.2 ms | 6.75 ms | 8.04 ms | 37.3 ms | 1.2x faster |
| `gaussian/k50/n1000/scattered` | 25.8 ms | 5.24 ms | 23.2 ms | 68.5 ms |  | 2.7x faster |
| `thin_plate_spline/global/n3000/scattered` | 21.3 ms | 4.47 ms | 21.3 ms | 30 ms | 92.8 ms | 1.4x faster |
| `thin_plate_spline/k50/n10000/scattered` | 26.4 ms | 5.12 ms | 24.4 ms | 63.1 ms |  | 2.4x faster |
| `thin_plate_spline/k50/n10000/image` | 234 ms | 39 ms | 234 ms | 545 ms |  | 2.3x faster |
| `thin_plate_spline/k50/n100000/scattered` | 25 ms | 4.54 ms | 24 ms | 46.6 ms |  | 1.9x faster |
| `thin_plate_spline/k50/n100000/image` | 235 ms | 38 ms | 237 ms | 459 ms |  | 2x faster |

### Construction

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/global/n10/scattered` | 1.36 µs | 1.36 µs | 1.39 µs | 11.8 µs | 5.9 µs | 4.3x faster |
| `linear/k50/n10/scattered` | 1.31 µs | 1.37 µs | 1.37 µs | 16.4 µs |  | 13x faster |
| `thin_plate_spline/global/n10/scattered` | 1.96 µs | 1.82 µs | 1.88 µs | 13.3 µs | 7.17 µs | 3.7x faster |
| `thin_plate_spline/k50/n10/scattered` | 1.95 µs | 1.92 µs | 1.9 µs | 17.1 µs |  | 8.8x faster |
| `cubic/global/n10/scattered` | 1.69 µs | 1.69 µs | 1.65 µs | 13.3 µs | 6.01 µs | 3.5x faster |
| `cubic/k50/n10/scattered` | 1.73 µs | 1.68 µs | 1.67 µs | 17 µs |  | 9.9x faster |
| `quintic/global/n10/scattered` | 2.24 µs | 2.15 µs | 2.17 µs | 14.1 µs | 5.74 µs | 2.6x faster |
| `quintic/k50/n10/scattered` | 2.2 µs | 2.15 µs | 2.18 µs | 18.4 µs |  | 8.4x faster |
| `multiquadric/global/n10/scattered` | 1.33 µs | 1.32 µs | 1.32 µs | 10.8 µs | 6.78 µs | 5.1x faster |
| `multiquadric/k50/n10/scattered` | 1.42 µs | 1.44 µs | 1.38 µs | 15 µs |  | 11x faster |
| `inverse_multiquadric/global/n10/scattered` | 1.45 µs | 1.34 µs | 1.39 µs | 10.9 µs | 7.35 µs | 5.1x faster |
| `inverse_multiquadric/k50/n10/scattered` | 1.44 µs | 1.44 µs | 1.44 µs | 14.9 µs |  | 10x faster |
| `inverse_quadratic/global/n10/scattered` | 1.34 µs | 1.29 µs | 1.38 µs | 11.8 µs |  | 8.8x faster |
| `inverse_quadratic/k50/n10/scattered` | 1.37 µs | 1.34 µs | 1.3 µs | 16.3 µs |  | 12x faster |
| `gaussian/global/n10/scattered` | 1.6 µs | 1.5 µs | 1.44 µs | 11.1 µs | 8.29 µs | 5.2x faster |
| `gaussian/k50/n10/scattered` | 1.46 µs | 1.49 µs | 1.48 µs | 16.1 µs |  | 11x faster |
| `linear/global/n100/scattered` | 43.9 µs | 76.4 µs | 53.1 µs | 81.5 µs | 719 µs | 1.9x faster |
| `linear/k50/n100/scattered` | 3.14 µs | 3.44 µs | 3.47 µs | 22.7 µs |  | 7.2x faster |
| `thin_plate_spline/global/n100/scattered` | 79.5 µs | 97.6 µs | 100 µs | 114 µs | 907 µs | 1.4x faster |
| `thin_plate_spline/k50/n100/scattered` | 3.55 µs | 3.54 µs | 3.98 µs | 21.9 µs |  | 6.2x faster |
| `cubic/global/n100/scattered` | 52.2 µs | 87.3 µs | 72.6 µs | 87.8 µs | 798 µs | 1.7x faster |
| `cubic/k50/n100/scattered` | 3.5 µs | 3.6 µs | 3.45 µs | 23.3 µs |  | 6.7x faster |
| `quintic/global/n100/scattered` | 56.4 µs | 80.2 µs | 75.4 µs | 99.8 µs | 775 µs | 1.8x faster |
| `quintic/k50/n100/scattered` | 3.25 µs | 3.55 µs | 3.55 µs | 22.8 µs |  | 7x faster |
| `multiquadric/global/n100/scattered` | 52.7 µs | 84.3 µs | 66.8 µs | 79 µs | 681 µs | 1.5x faster |
| `multiquadric/k50/n100/scattered` | 3.48 µs | 3.47 µs | 3.35 µs | 22.6 µs |  | 6.5x faster |
| `inverse_multiquadric/global/n100/scattered` | 59.6 µs | 80.1 µs | 75.5 µs | 86.7 µs | 872 µs | 1.5x faster |
| `inverse_multiquadric/k50/n100/scattered` | 3.45 µs | 3.45 µs | 3.39 µs | 20.9 µs |  | 6.1x faster |
| `inverse_quadratic/global/n100/scattered` | 52.6 µs | 78.7 µs | 53.2 µs | 74.3 µs |  | 1.4x faster |
| `inverse_quadratic/k50/n100/scattered` | 3.41 µs | 3.48 µs | 3.15 µs | 22.5 µs |  | 6.6x faster |
| `gaussian/global/n100/scattered` | 72.9 µs | 89.7 µs | 91.5 µs | 99.3 µs | 766 µs | 1.4x faster |
| `gaussian/k50/n100/scattered` | 3.23 µs | 3.21 µs | 3.18 µs | 20.7 µs |  | 6.4x faster |
| `linear/global/n1000/scattered` | 15.6 ms | 15.8 ms | 18.1 ms | 59 ms | 190 ms | 3.8x faster |
| `linear/k50/n1000/scattered` | 65 µs | 59.5 µs | 55.1 µs | 140 µs |  | 2.1x faster |
| `thin_plate_spline/global/n1000/scattered` | 21.1 ms | 20 ms | 21.7 ms | 64.4 ms | 210 ms | 3.1x faster |
| `thin_plate_spline/k50/n1000/scattered` | 61.5 µs | 61.5 µs | 57 µs | 211 µs |  | 3.4x faster |
| `cubic/global/n1000/scattered` | 17.6 ms | 17.8 ms | 17.8 ms | 59.4 ms | 191 ms | 3.4x faster |
| `cubic/k50/n1000/scattered` | 61.9 µs | 59.7 µs | 56.4 µs | 131 µs |  | 2.1x faster |
| `quintic/global/n1000/scattered` | 17.4 ms | 18.4 ms | 17.2 ms | 60.7 ms | 191 ms | 3.5x faster |
| `quintic/k50/n1000/scattered` | 62.8 µs | 57.6 µs | 55.1 µs | 144 µs |  | 2.3x faster |
| `multiquadric/global/n1000/scattered` | 18.2 ms | 18.5 ms | 18 ms | 60.7 ms | 195 ms | 3.3x faster |
| `multiquadric/k50/n1000/scattered` | 66.2 µs | 57.2 µs | 59.2 µs | 141 µs |  | 2.1x faster |
| `inverse_multiquadric/global/n1000/scattered` | 19.5 ms | 18.9 ms | 19.2 ms | 60.9 ms | 196 ms | 3.1x faster |
| `inverse_multiquadric/k50/n1000/scattered` | 62.8 µs | 54.2 µs | 55.8 µs | 128 µs |  | 2x faster |
| `inverse_quadratic/global/n1000/scattered` | 17.5 ms | 17 ms | 18.7 ms | 59.6 ms |  | 3.4x faster |
| `inverse_quadratic/k50/n1000/scattered` | 61.1 µs | 55.1 µs | 56.5 µs | 133 µs |  | 2.2x faster |
| `gaussian/global/n1000/scattered` | 19.6 ms | 20.8 ms | 20.6 ms | 60.6 ms | 206 ms | 3.1x faster |
| `gaussian/k50/n1000/scattered` | 69.5 µs | 63.8 µs | 88.3 µs | 211 µs |  | 3x faster |
| `thin_plate_spline/global/n3000/scattered` | 414 ms | 413 ms | 409 ms | 1.89 s | 6.02 s | 4.6x faster |
| `thin_plate_spline/k50/n10000/scattered` | 760 µs | 791 µs | 717 µs | 1.83 ms |  | 2.4x faster |
| `thin_plate_spline/k50/n10000/image` | 824 µs | 801 µs | 746 µs | 1.66 ms |  | 2x faster |
| `thin_plate_spline/k50/n100000/scattered` | 8.87 ms | 9.33 ms | 8.8 ms | 25 ms |  | 2.8x faster |
| `thin_plate_spline/k50/n100000/image` | 9.36 ms | 8.94 ms | 8.81 ms | 25.8 ms |  | 2.8x faster |

## Memory

Megabytes. "Data read" is the process with the data of the case in memory, before anything is built. The columns after it are what was added to that: the highest while building, with the interpolator built, the highest while evaluating. "Highest" is the highest of the whole run, repeats included.

| case | library | data read | while building | built | while evaluating | highest |
| :--- | :--- | ---: | ---: | ---: | ---: | ---: |
| `oned/cubic_spline/n10000/sorted` | interpx | 3.3 | +1.9 | +2.0 | +2.0 | 5.9 |
| `oned/cubic_spline/n10000/sorted` | interpx, parallel | 3.4 | +1.8 | +1.9 | +1.9 | 6.4 |
| `oned/cubic_spline/n10000/sorted` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |
| `oned/cubic_spline/n10000/sorted` | scirs2-interpolate | 3.1 | +2.0 | +2.1 | +2.7 | 5.9 |
| `oned/cubic_spline/n10000/sorted` | ndarray-interp | 2.9 | +2.1 | +2.1 | +2.1 | 5.4 |
| `grid/linear/4d/max200/eval10000/ascending` | interpx | 209.0 | +196.0 | +195.7 | +195.7 | 405.0 |
| `grid/linear/4d/max200/eval10000/ascending` | interpx, parallel | 220.9 | +183.9 | +183.9 | +184.0 | 405.4 |
| `grid/linear/4d/max200/eval10000/ascending` | interpx, `from_owned` | 221.6 | +183.4 | +183.4 | +183.4 | 405.1 |
| `grid/linear/4d/max200/eval10000/ascending` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |
| `grid/linear/4d/max200/eval10000/ascending` | scirs2-interpolate | 204.3 | +198.5 | +200.6 | +200.6 | 404.9 |
| `grid/linear/4d/max200/eval10000/ascending` | interpn | 221.7 | +0.0 | +0.0 | +0.0 | 401.5 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpx | 5.6 | +0.1 | +2.7 | +2.7 | 8.5 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpx, parallel | 5.8 | +0.0 | +2.6 | +2.6 | 8.9 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpx, `from_owned` | 5.6 | +0.1 | +2.7 | +2.7 | 8.5 |
| `grid/cubic/2d/each400/eval1000/ascending` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpn | 3.3 | +2.1 | +2.1 | +2.1 | 5.5 |
| `grid/linear/2d/each800/eval10000/ascending` | interpx | 12.9 | +0.0 | +2.1 | +2.1 | 15.0 |
| `grid/linear/2d/each800/eval10000/ascending` | interpx, parallel | 13.1 | +0.0 | +1.7 | +1.7 | 15.4 |
| `grid/linear/2d/each800/eval10000/ascending` | interpx, `from_owned` | 11.6 | +1.3 | +3.2 | +3.2 | 14.9 |
| `grid/linear/2d/each800/eval10000/ascending` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |
| `grid/linear/2d/each800/eval10000/ascending` | scirs2-interpolate | 12.1 | +0.0 | +2.7 | +2.8 | 14.9 |
| `grid/linear/2d/each800/eval10000/ascending` | ndarray-interp | 12.9 | +0.0 | +0.0 | +0.0 | 12.9 |
| `grid/linear/2d/each800/eval10000/ascending` | interpn | 13.6 | +0.0 | +0.0 | +0.0 | 13.6 |
| `nearest_nd/n100000/image500` | interpx | 14.5 | +0.9 | +0.0 | +1.9 | 17.8 |
| `nearest_nd/n100000/image500` | interpx, parallel | 13.7 | +1.8 | +0.6 | +3.1 | 18.7 |
| `nearest_nd/n100000/image500` | interpx, `from_owned` | 12.6 | +3.5 | +2.0 | +4.0 | 18.6 |
| `nearest_nd/n100000/image500` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |
| `rbf/thin_plate_spline/global/n3000/scattered` | interpx | 3.3 | +77.0 | +4.8 | +4.8 | 80.3 |
| `rbf/thin_plate_spline/global/n3000/scattered` | interpx, parallel | 3.5 | +76.8 | +4.6 | +5.0 | 80.8 |
| `rbf/thin_plate_spline/global/n3000/scattered` | interpx, `from_owned` | 4.3 | +76.0 | +3.8 | +3.8 | 80.4 |
| `rbf/thin_plate_spline/global/n3000/scattered` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |
| `rbf/thin_plate_spline/global/n3000/scattered` | scirs2-interpolate | 3.0 | +289.5 | +1.5 | +1.5 | 292.5 |
| `rbf/thin_plate_spline/k50/n100000/image` | interpx | 5.9 | +4.5 | +5.5 | +5.8 | 11.9 |
| `rbf/thin_plate_spline/k50/n100000/image` | interpx, parallel | 6.7 | +3.5 | +4.7 | +5.8 | 12.7 |
| `rbf/thin_plate_spline/k50/n100000/image` | interpx, `from_owned` | 6.9 | +5.8 | +4.5 | +4.9 | 13.5 |
| `rbf/thin_plate_spline/k50/n100000/image` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |

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
| Regular grid | 4e-04 | 4e-04 | 4e-04 | 4e-16 | 6e-16 | 7e-01 |
| `NearestNDInterpolator` | 0e+00 | 0e+00 | 0e+00 |  |  |  |
| `RBFInterpolator` | 8e-05 | 8e-05 | 8e-05 | 3e+01 |  |  |
