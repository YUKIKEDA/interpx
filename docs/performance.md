# Performance: interpx against SciPy and other Rust crates

<!-- Written by compare/perf/report.py from a run of compare/perf/run.py. Do not edit by hand: run it again. -->

Measured on 2026-10-13, at commit `78050eaa9b` with uncommitted changes.

## How it was measured

- Machine: Intel64 Family 6 Model 191 Stepping 2, GenuineIntel, 16 logical cores, 48 GB. Windows-11-10.0.26300-SP0.
- rustc 1.97.1 (8bab26f4f 2026-07-14); Python 3.13.1, NumPy 2.5.4, SciPy 1.18.1.
- The other Rust crates, at the versions in `compare/perf/runners/*/Cargo.toml`: `scirs2-interpolate` (the widest of the Rust crates: one variable, grids, RBF), `ndarray-interp` (the most used for one variable; also bilinear), `interpn` (grids).
- Every library ran each case in a process of its own, pinned to one core, on one thread, on the same data. The libraries took turns case by case, 3 rounds.
- `interpx, from_owned` is interpx built from data that is handed over (`RegularGridInterpolator::from_owned`) and so not copied. Its column is there for the grids, next to the construction from borrowed data (`new`), which is the one the comparison counts. Each of its samples is given data made just before it, which is not in the nearest cache as the data of the other constructions is: on small grids it reads slower than `new` for that reason.
- A time is that of one call: construction, or the evaluation of all the queries of the case (`eval_batch_into` for interpx, the call with an array for the others). It is the median of the repeats of a process, and then the median of the rounds.
- interpx is called slower where its median is more than 1.10 times another library's and its fastest repeat is slower than that library's fastest.
- Memory is the resident set size of the process, sampled every millisecond from outside it.
- Command: `perf/run.py --rounds 3 --filter ^(oned|grid)/`.

## Where interpx is slower

| family | operation | cases | interpx fastest | slower than SciPy | slower than a Rust crate |
| :--- | :--- | ---: | ---: | ---: | ---: |
| One variable | evaluation | 112 | 111 | 1 | 0 |
| One variable | construction | 112 | 78 | 8 | 30 |
| `PPoly` | evaluation | 4 | 4 | 0 | 0 |
| `PPoly` | construction | 4 | 4 | 0 | 0 |
| `BSpline` | evaluation | 4 | 4 | 0 | 0 |
| `BSpline` | construction | 4 | 4 | 0 | 0 |
| `BSplineCurve` | evaluation | 2 | 2 | 0 | 0 |
| `BSplineCurve` | construction | 2 | 2 | 0 | 0 |
| Regular grid | evaluation | 79 | 50 | 0 | 29 |
| Regular grid | construction | 79 | 18 | 37 | 43 |
| `NearestNDInterpolator` | evaluation | 15 | 15 | 0 | 0 |
| `NearestNDInterpolator` | construction | 15 | 15 | 0 | 0 |
| `RBFInterpolator` | evaluation | 53 | 52 | 1 | 0 |
| `RBFInterpolator` | construction | 53 | 53 | 0 | 0 |

203 comparisons (a case, an operation, and a library), by method. "Same values" says whether that library's values are SciPy's up to rounding: where they are not, it computes another interpolant under the same name, and its time is not that of the same work.

| method | operation | against | same values | cases | slower by | the largest |
| :--- | :--- | :--- | :--- | ---: | ---: | :--- |
| `grid/linear` | construction | interpn | yes | 34 | 175x to 4.06e+07x | `grid/linear/4d/max200/eval10000/ascending`: 47.5 ms against 1.17 ns |
| `grid/cubic` | construction | interpn | no | 7 | 1.88e+03x to 1.6e+06x | `grid/cubic/2d/each400/eval1000/ascending`: 1.41 ms against 0 |
| `grid/nearest` | construction | interpn | yes | 2 | 951x to 2.61e+04x | `grid/nearest/3d/each40/eval10000/ascending`: 22.9 µs against 0 |
| `grid/linear` | construction | SciPy | yes | 35 | 1.17x to 5.62e+03x | `grid/linear/4d/max200/eval100/flipped`: 114 ms against 20.2 µs |
| `grid/linear` | construction | ndarray-interp | yes | 14 | 1.59x to 441x | `grid/linear/2d/each800/eval10000/ascending`: 1.02 ms against 2.32 µs |
| `grid/pchip` | construction | SciPy | yes | 1 | 53.8x | `grid/pchip/3d/each40/eval10000/ascending`: 465 µs against 8.63 µs |
| `oned/linear` | construction | ndarray-interp | yes | 14 | 1.4x to 6.76x | `oned/linear/n10000/sorted-uneven`: 97.6 µs against 14.5 µs |
| `oned/bspline1` | construction | SciPy | yes | 2 | 2.99x to 3.28x | `oned/bspline1/n10000/sorted`: 75.5 µs against 23.1 µs |
| `grid/nearest` | construction | SciPy | yes | 1 | 3.05x | `grid/nearest/3d/each40/eval10000/ascending`: 22.9 µs against 7.51 µs |
| `grid/nearest` | construction | scirs2-interpolate | yes | 2 | 2.4x to 2.69x | `grid/nearest/3d/each40/eval10000/ascending`: 22.9 µs against 8.53 µs |
| `oned/pchip` | construction | scirs2-interpolate | yes | 8 | 1.41x to 2.58x | `oned/pchip/n1000/shuffled`: 11.8 µs against 4.59 µs |
| `grid/linear` | construction | scirs2-interpolate | yes | 25 | 1.11x to 2.54x | `grid/linear/3d/max40/eval1000/ascending`: 15.5 µs against 6.09 µs |
| `oned/cubic_spline` | construction | ndarray-interp | yes | 2 | 2.53x to 2.53x | `oned/cubic_spline/n1000/sorted`: 40.6 µs against 16.1 µs |
| `oned/linear` | construction | SciPy | yes | 4 | 2.08x to 2.22x | `oned/linear/n10000/sorted`: 99.4 µs against 44.7 µs |
| `oned/cubic_spline` | construction | scirs2-interpolate | no | 1 | 2.17x | `oned/cubic_spline/n1000/shuffled`: 39.8 µs against 18.3 µs |
| `oned/akima` | construction | scirs2-interpolate | yes | 6 | 1.59x to 2x | `oned/akima/n1000/shuffled`: 12.7 µs against 6.34 µs |
| `oned/linear` | construction | scirs2-interpolate | yes | 12 | 1.29x to 1.84x | `oned/linear/n1000/shuffled`: 2.12 µs against 1.16 µs |
| `grid/linear` | evaluation | interpn | yes | 28 | 1.11x to 1.78x | `grid/linear/4d/max40/eval10/ascending`: 480 ns against 270 ns |
| `oned/nearest` | construction | SciPy | yes | 2 | 1.14x to 1.18x | `oned/nearest/n10000/sorted`: 61.2 µs against 52 µs |
| `grid/nearest` | evaluation | interpn | yes | 1 | 1.12x | `grid/nearest/2d/each40/eval10000/ascending`: 159 µs against 141 µs |
| `rbf/quintic` | evaluation | SciPy | yes | 1 | 1.12x | `rbf/quintic/global/n10/scattered`: 75 µs against 66.8 µs |
| `oned/bspline1` | evaluation | SciPy | yes | 1 | 1.1x | `oned/bspline1/n10000/sorted`: 314 µs against 285 µs |

<details><summary>Every comparison where interpx is slower</summary>

"Values apart" is the distance of that library's values from SciPy's, as a share of the largest value.

| case | operation | against | interpx | the other | slower by | values apart |
| :--- | :--- | :--- | ---: | ---: | ---: | ---: |
| `grid/linear/4d/max200/eval10000/ascending` | construction | interpn | 47.5 ms | 1.17 ns | 40563968.00x | 4e-16 |
| `grid/linear/4d/max200/eval10/ascending` | construction | interpn | 47.2 ms | 1.17 ns | 40288085.33x | 2e-16 |
| `grid/linear/4d/max200/eval1000/ascending` | construction | interpn | 49.8 ms | 1.27 ns | 39233772.31x | 3e-16 |
| `grid/linear/4d/max200/eval100/ascending` | construction | interpn | 48.8 ms | 1.27 ns | 38400945.23x | 3e-16 |
| `grid/linear/3d/max200/eval100/ascending` | construction | interpn | 1.74 ms | 0 | 1983829.33x | 3e-16 |
| `grid/linear/3d/max200/eval10/ascending` | construction | interpn | 1.67 ms | 0 | 1899975.11x | 3e-16 |
| `grid/linear/3d/max200/eval10000/ascending` | construction | interpn | 1.78 ms | 0 | 1822515.20x | 5e-16 |
| `grid/linear/3d/max200/eval1000/ascending` | construction | interpn | 1.74 ms | 0 | 1777152.00x | 3e-16 |
| `grid/cubic/2d/each400/eval1000/ascending` | construction | interpn | 1.41 ms | 0 | 1602730.67x | 5e-01 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | interpn | 1.02 ms | 0 | 1494601.14x | 3e-16 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | interpn | 302 µs | 0 | 441344.00x | 3e-16 |
| `grid/cubic/3d/each40/eval10000/ascending` | construction | interpn | 436 µs | 1.17 ns | 372437.33x | 5e-01 |
| `grid/cubic/2d/each200/eval1000/ascending` | construction | interpn | 175 µs | 0 | 199452.44x | 2e-01 |
| `grid/cubic/2d/each100/eval1000/ascending` | construction | interpn | 76.5 µs | 0 | 87068.44x | 2e-01 |
| `grid/linear/2d/max200/eval10000/ascending` | construction | interpn | 49.1 µs | 0 | 71789.71x | 3e-16 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | interpn | 22.9 µs | 0 | 26083.56x | 0e+00 |
| `grid/linear/3d/max40/eval1000/ascending` | construction | interpn | 15.5 µs | 0 | 17607.11x | 2e-16 |
| `grid/linear/3d/max40/eval10000/ascending` | construction | interpn | 15.5 µs | 0 | 17600.00x | 3e-16 |
| `grid/linear/4d/max40/eval10/ascending` | construction | interpn | 14.3 µs | 1.17 ns | 12176.00x | 2e-16 |
| `grid/linear/4d/max40/eval10000/ascending` | construction | interpn | 14.4 µs | 1.27 ns | 11328.00x | 5e-16 |
| `grid/linear/4d/max40/eval1000/ascending` | construction | interpn | 14.3 µs | 1.27 ns | 11303.38x | 6e-16 |
| `grid/linear/2d/max200/eval10/ascending` | construction | interpn | 7.54 µs | 0 | 11035.43x | 2e-16 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | interpn | 7.44 µs | 0 | 10880.00x | 2e-16 |
| `grid/linear/2d/max200/eval100/ascending` | construction | interpn | 7.43 µs | 0 | 10861.71x | 1e-16 |
| `grid/linear/4d/max40/eval100/ascending` | construction | interpn | 13.4 µs | 1.27 ns | 10520.62x | 2e-16 |
| `grid/cubic/2d/each40/eval10000/ascending` | construction | interpn | 8.58 µs | 0 | 9760.00x | 7e-01 |
| `grid/cubic/2d/each40/eval1000/ascending` | construction | interpn | 8.57 µs | 0 | 9751.11x | 7e-01 |
| `grid/linear/4d/max200/eval100/flipped` | construction | SciPy | 114 ms | 20.2 µs | 5617.96x |  |
| `grid/linear/4d/max200/eval100/ascending` | construction | SciPy | 48.8 ms | 9.47 µs | 5146.94x |  |
| `grid/linear/4d/max200/eval10000/ascending` | construction | SciPy | 47.5 ms | 9.31 µs | 5105.38x |  |
| `grid/linear/4d/max200/eval10/ascending` | construction | SciPy | 47.2 ms | 9.38 µs | 5033.49x |  |
| `grid/linear/4d/max200/eval1000/flipped` | construction | SciPy | 100 ms | 20 µs | 5004.77x |  |
| `grid/linear/4d/max200/eval10000/flipped` | construction | SciPy | 98.9 ms | 19.9 µs | 4967.83x |  |
| `grid/linear/4d/max200/eval1000/ascending` | construction | SciPy | 49.8 ms | 10.2 µs | 4869.01x |  |
| `grid/linear/4d/max200/eval10/flipped` | construction | SciPy | 106 ms | 22.1 µs | 4800.01x |  |
| `grid/linear/3d/max40/eval10/ascending` | construction | interpn | 2.74 µs | 0 | 2801.60x | 2e-16 |
| `grid/linear/3d/max40/eval100/ascending` | construction | interpn | 2.72 µs | 0 | 2784.40x | 2e-16 |
| `grid/cubic/2d/each10/eval1000/ascending` | construction | interpn | 1.47 µs | 0 | 1880.00x | 3e-01 |
| `grid/nearest/2d/each40/eval10000/ascending` | construction | interpn | 650 ns | 0 | 950.86x | 0e+00 |
| `grid/linear/2d/max40/eval10/ascending` | construction | interpn | 408 ns | 0 | 596.29x | 2e-16 |
| `grid/linear/2d/max40/eval1000/ascending` | construction | interpn | 371 ns | 0 | 542.29x | 3e-16 |
| `grid/linear/2d/max40/eval10000/ascending` | construction | interpn | 370 ns | 0 | 541.43x | 4e-16 |
| `grid/linear/2d/max40/eval100/ascending` | construction | interpn | 406 ns | 0 | 519.25x | 2e-16 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | ndarray-interp | 1.02 ms | 2.32 µs | 441.22x | 6e-16 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | ndarray-interp | 302 µs | 1.14 µs | 265.00x | 5e-16 |
| `grid/linear/3d/max200/eval10000/ascending` | construction | SciPy | 1.78 ms | 7.61 µs | 233.80x |  |
| `grid/linear/3d/max200/eval1000/ascending` | construction | SciPy | 1.74 ms | 7.68 µs | 225.94x |  |
| `grid/linear/3d/max200/eval10/ascending` | construction | SciPy | 1.67 ms | 7.48 µs | 223.12x |  |
| `grid/linear/3d/max200/eval100/ascending` | construction | SciPy | 1.74 ms | 8.28 µs | 210.55x |  |
| `grid/linear/3d/max200/eval10/flipped` | construction | SciPy | 3.27 ms | 16.1 µs | 203.58x |  |
| `grid/linear/3d/max10/eval10000/ascending` | construction | interpn | 197 ns | 0 | 201.40x | 4e-16 |
| `grid/linear/3d/max10/eval1000/ascending` | construction | interpn | 196 ns | 0 | 200.40x | 5e-16 |
| `grid/linear/3d/max10/eval10/ascending` | construction | interpn | 195 ns | 0 | 200.00x | 2e-16 |
| `grid/linear/3d/max10/eval100/ascending` | construction | interpn | 194 ns | 0 | 199.15x | 3e-16 |
| `grid/linear/2d/max10/eval10000/ascending` | construction | interpn | 135 ns | 0 | 198.14x | 3e-16 |
| `grid/linear/3d/max200/eval10000/flipped` | construction | SciPy | 3.15 ms | 16 µs | 197.21x |  |
| `grid/linear/3d/max200/eval1000/flipped` | construction | SciPy | 3.11 ms | 15.9 µs | 195.59x |  |
| `grid/linear/3d/max200/eval100/flipped` | construction | SciPy | 3.06 ms | 15.9 µs | 192.36x |  |
| `grid/linear/2d/max10/eval100/ascending` | construction | interpn | 137 ns | 0 | 175.88x | 2e-16 |
| `grid/linear/2d/max10/eval1000/ascending` | construction | interpn | 137 ns | 0 | 175.50x | 3e-16 |
| `grid/linear/2d/max10/eval10/ascending` | construction | interpn | 137 ns | 0 | 175.38x | 1e-16 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | SciPy | 1.02 ms | 6.33 µs | 161.45x |  |
| `grid/linear/2d/max200/eval10000/ascending` | construction | ndarray-interp | 49.1 µs | 512 ns | 95.83x | 5e-16 |
| `grid/pchip/3d/each40/eval10000/ascending` | construction | SciPy | 465 µs | 8.63 µs | 53.84x |  |
| `grid/linear/2d/each400/eval10000/ascending` | construction | SciPy | 302 µs | 5.99 µs | 50.40x |  |
| `grid/linear/2d/max200/eval10/ascending` | construction | ndarray-interp | 7.54 µs | 514 ns | 14.67x | 2e-16 |
| `grid/linear/2d/max200/eval100/ascending` | construction | ndarray-interp | 7.43 µs | 513 ns | 14.48x | 3e-16 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | ndarray-interp | 7.44 µs | 514 ns | 14.48x | 4e-16 |
| `grid/linear/2d/max200/eval10000/ascending` | construction | SciPy | 49.1 µs | 5.81 µs | 8.45x |  |
| `oned/linear/n10000/sorted-uneven` | construction | ndarray-interp | 97.6 µs | 14.5 µs | 6.76x | 0e+00 |
| `oned/linear/n10000/sorted` | construction | ndarray-interp | 99.4 µs | 21.3 µs | 4.66x | 0e+00 |
| `grid/linear/2d/max200/eval10000/flipped` | construction | SciPy | 58.7 µs | 12.9 µs | 4.57x |  |
| `oned/bspline1/n10000/sorted` | construction | SciPy | 75.5 µs | 23.1 µs | 3.28x |  |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | SciPy | 22.9 µs | 7.51 µs | 3.05x |  |
| `oned/linear/n10000/shuffled` | construction | ndarray-interp | 96.4 µs | 32.2 µs | 2.99x | 0e+00 |
| `oned/bspline1/n10000/shuffled` | construction | SciPy | 75.8 µs | 25.4 µs | 2.99x |  |
| `grid/linear/4d/max40/eval100/flipped` | construction | SciPy | 57.9 µs | 19.8 µs | 2.92x |  |
| `grid/linear/4d/max40/eval10000/flipped` | construction | SciPy | 57.3 µs | 19.9 µs | 2.87x |  |
| `oned/linear/n10000/shuffled-uneven` | construction | ndarray-interp | 94.7 µs | 33.1 µs | 2.86x | 0e+00 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | scirs2-interpolate | 22.9 µs | 8.53 µs | 2.69x | 0e+00 |
| `grid/linear/4d/max40/eval1000/flipped` | construction | SciPy | 57.6 µs | 21.5 µs | 2.69x |  |
| `oned/pchip/n1000/shuffled` | construction | scirs2-interpolate | 11.8 µs | 4.59 µs | 2.58x | 2e-16 |
| `oned/pchip/n1000/sorted` | construction | scirs2-interpolate | 11.8 µs | 4.58 µs | 2.57x | 2e-16 |
| `grid/linear/3d/max40/eval1000/ascending` | construction | scirs2-interpolate | 15.5 µs | 6.09 µs | 2.54x | 2e-16 |
| `oned/cubic_spline/n1000/sorted` | construction | ndarray-interp | 40.6 µs | 16.1 µs | 2.53x | 2e-16 |
| `oned/cubic_spline/n1000/shuffled` | construction | ndarray-interp | 39.8 µs | 15.8 µs | 2.53x | 2e-16 |
| `grid/linear/2d/max40/eval100/ascending` | construction | ndarray-interp | 406 ns | 164 ns | 2.48x | 2e-16 |
| `grid/linear/4d/max40/eval10/flipped` | construction | SciPy | 53 µs | 21.4 µs | 2.47x |  |
| `grid/linear/2d/max40/eval10000/ascending` | construction | ndarray-interp | 370 ns | 150 ns | 2.47x | 4e-16 |
| `grid/linear/2d/max40/eval10/ascending` | construction | ndarray-interp | 408 ns | 165 ns | 2.47x | 2e-16 |
| `grid/linear/4d/max40/eval10/ascending` | construction | scirs2-interpolate | 14.3 µs | 5.91 µs | 2.42x | 2e-16 |
| `grid/nearest/2d/each40/eval10000/ascending` | construction | scirs2-interpolate | 650 ns | 271 ns | 2.40x | 0e+00 |
| `grid/linear/4d/max40/eval10000/ascending` | construction | scirs2-interpolate | 14.4 µs | 5.99 µs | 2.40x | 4e-16 |
| `grid/linear/4d/max40/eval1000/ascending` | construction | scirs2-interpolate | 14.3 µs | 6.05 µs | 2.37x | 2e-16 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | scirs2-interpolate | 7.44 µs | 3.15 µs | 2.36x | 2e-16 |
| `grid/linear/2d/max200/eval100/ascending` | construction | scirs2-interpolate | 7.43 µs | 3.16 µs | 2.35x | 2e-16 |
| `grid/linear/2d/max200/eval10/ascending` | construction | scirs2-interpolate | 7.54 µs | 3.3 µs | 2.29x | 1e-16 |
| `grid/linear/4d/max40/eval100/ascending` | construction | scirs2-interpolate | 13.4 µs | 5.84 µs | 2.29x | 3e-16 |
| `grid/linear/2d/max40/eval1000/ascending` | construction | ndarray-interp | 371 ns | 164 ns | 2.27x | 4e-16 |
| `oned/linear/n10000/sorted` | construction | SciPy | 99.4 µs | 44.7 µs | 2.22x |  |
| `oned/linear/n10/sorted` | construction | ndarray-interp | 94.9 ns | 43.5 ns | 2.18x | 0e+00 |
| `oned/cubic_spline/n1000/shuffled` | construction | scirs2-interpolate | 39.8 µs | 18.3 µs | 2.17x | 4e-04 |
| `oned/linear/n10000/shuffled-uneven` | construction | SciPy | 94.7 µs | 45.6 µs | 2.08x |  |
| `oned/linear/n10000/sorted-uneven` | construction | SciPy | 97.6 µs | 47 µs | 2.08x |  |
| `oned/linear/n10000/shuffled` | construction | SciPy | 96.4 µs | 46.4 µs | 2.08x |  |
| `oned/linear/n10/shuffled` | construction | ndarray-interp | 94.5 ns | 47.1 ns | 2.01x | 0e+00 |
| `oned/akima/n1000/shuffled` | construction | scirs2-interpolate | 12.7 µs | 6.34 µs | 2.00x | 3e-17 |
| `oned/akima/n1000/sorted` | construction | scirs2-interpolate | 12 µs | 6.34 µs | 1.90x | 3e-17 |
| `grid/linear/3d/max40/eval10000/flipped` | construction | SciPy | 28.6 µs | 15.5 µs | 1.84x |  |
| `oned/linear/n1000/shuffled` | construction | scirs2-interpolate | 2.12 µs | 1.16 µs | 1.84x | 0e+00 |
| `oned/linear/n1000/sorted` | construction | scirs2-interpolate | 2.12 µs | 1.17 µs | 1.82x | 0e+00 |
| `grid/linear/4d/max40/eval10/ascending` | evaluation | interpn | 480 ns | 270 ns | 1.78x | 2e-16 |
| `grid/linear/2d/max40/eval100/ascending` | construction | scirs2-interpolate | 406 ns | 229 ns | 1.77x | 3e-16 |
| `grid/linear/2d/max40/eval10/ascending` | construction | scirs2-interpolate | 408 ns | 230 ns | 1.77x | 4e-17 |
| `oned/akima/n50/shuffled` | construction | scirs2-interpolate | 574 ns | 329 ns | 1.74x | 1e-16 |
| `oned/pchip/n100/sorted` | construction | scirs2-interpolate | 956 ns | 554 ns | 1.73x | 2e-16 |
| `oned/akima/n100/sorted` | construction | scirs2-interpolate | 1.03 µs | 602 ns | 1.71x | 1e-16 |
| `oned/akima/n100/shuffled` | construction | scirs2-interpolate | 1.04 µs | 613 ns | 1.70x | 1e-16 |
| `grid/linear/2d/max40/eval1000/ascending` | construction | scirs2-interpolate | 371 ns | 224 ns | 1.66x | 3e-16 |
| `grid/linear/2d/max40/eval10000/ascending` | construction | scirs2-interpolate | 370 ns | 224 ns | 1.65x | 2e-16 |
| `grid/linear/2d/max10/eval10/ascending` | construction | ndarray-interp | 137 ns | 82.9 ns | 1.65x | 1e-16 |
| `oned/linear/n10000/sorted-uneven` | construction | scirs2-interpolate | 97.6 µs | 59.5 µs | 1.64x | 6e-17 |
| `grid/linear/4d/max40/eval100/ascending` | evaluation | interpn | 4.72 µs | 2.93 µs | 1.61x | 2e-16 |
| `oned/linear/n10000/shuffled` | construction | scirs2-interpolate | 96.4 µs | 59.9 µs | 1.61x | 0e+00 |
| `grid/linear/2d/max10/eval100/ascending` | construction | ndarray-interp | 137 ns | 85.5 ns | 1.61x | 4e-16 |
| `oned/linear/n10000/sorted` | construction | scirs2-interpolate | 99.4 µs | 62 µs | 1.60x | 0e+00 |
| `grid/linear/4d/max40/eval1000/ascending` | evaluation | interpn | 44.1 µs | 27.6 µs | 1.60x | 6e-16 |
| `grid/linear/2d/max10/eval1000/ascending` | construction | ndarray-interp | 137 ns | 85.8 ns | 1.60x | 3e-16 |
| `grid/linear/2d/max10/eval10000/ascending` | construction | ndarray-interp | 135 ns | 85.4 ns | 1.59x | 4e-16 |
| `oned/akima/n50/sorted` | construction | scirs2-interpolate | 571 ns | 360 ns | 1.59x | 1e-16 |
| `oned/pchip/n100/sorted-uneven` | construction | scirs2-interpolate | 953 ns | 602 ns | 1.58x | 2e-16 |
| `oned/linear/n100/sorted-uneven` | construction | ndarray-interp | 266 ns | 170 ns | 1.56x | 0e+00 |
| `oned/linear/n10000/shuffled-uneven` | construction | scirs2-interpolate | 94.7 µs | 60.8 µs | 1.56x | 6e-17 |
| `grid/linear/4d/max40/eval10/ascending` | construction | SciPy | 14.3 µs | 9.52 µs | 1.50x |  |
| `oned/pchip/n100/shuffled-uneven` | construction | scirs2-interpolate | 936 ns | 625 ns | 1.50x | 2e-16 |
| `oned/linear/n100/sorted` | construction | scirs2-interpolate | 265 ns | 180 ns | 1.47x | 0e+00 |
| `grid/linear/4d/max40/eval10000/ascending` | evaluation | interpn | 431 µs | 294 µs | 1.47x | 5e-16 |
| `oned/pchip/n50/sorted` | construction | scirs2-interpolate | 530 ns | 365 ns | 1.45x | 1e-16 |
| `grid/linear/3d/max200/eval10000/ascending` | construction | scirs2-interpolate | 1.78 ms | 1.23 ms | 1.45x | 3e-16 |
| `oned/pchip/n100/shuffled` | construction | scirs2-interpolate | 873 ns | 605 ns | 1.44x | 2e-16 |
| `grid/linear/3d/max10/eval10/ascending` | evaluation | interpn | 218 ns | 152 ns | 1.44x | 2e-16 |
| `grid/linear/4d/max200/eval10/ascending` | evaluation | interpn | 515 ns | 359 ns | 1.43x | 2e-16 |
| `oned/linear/n100/shuffled-uneven` | construction | ndarray-interp | 268 ns | 187 ns | 1.43x | 0e+00 |
| `oned/linear/n100/shuffled` | construction | ndarray-interp | 243 ns | 171 ns | 1.43x | 0e+00 |
| `oned/linear/n100/sorted` | construction | ndarray-interp | 265 ns | 187 ns | 1.42x | 0e+00 |
| `grid/linear/4d/max40/eval1000/ascending` | construction | SciPy | 14.3 µs | 10.1 µs | 1.42x |  |
| `oned/linear/n50/shuffled` | construction | ndarray-interp | 163 ns | 115 ns | 1.42x | 0e+00 |
| `oned/pchip/n50/shuffled` | construction | scirs2-interpolate | 522 ns | 370 ns | 1.41x | 1e-16 |
| `oned/linear/n50/sorted` | construction | ndarray-interp | 163 ns | 116 ns | 1.41x | 0e+00 |
| `grid/linear/4d/max40/eval100/ascending` | construction | SciPy | 13.4 µs | 9.5 µs | 1.41x |  |
| `oned/linear/n50/sorted` | construction | scirs2-interpolate | 163 ns | 116 ns | 1.40x | 0e+00 |
| `oned/linear/n1000/sorted` | construction | ndarray-interp | 2.12 µs | 1.52 µs | 1.40x | 0e+00 |
| `oned/linear/n1000/shuffled` | construction | ndarray-interp | 2.12 µs | 1.52 µs | 1.40x | 0e+00 |
| `grid/linear/4d/max40/eval10000/ascending` | construction | SciPy | 14.4 µs | 10.4 µs | 1.39x |  |
| `oned/linear/n100/sorted-uneven` | construction | scirs2-interpolate | 266 ns | 194 ns | 1.37x | 3e-17 |
| `oned/linear/n100/shuffled-uneven` | construction | scirs2-interpolate | 268 ns | 197 ns | 1.36x | 3e-17 |
| `grid/linear/2d/max10/eval10/ascending` | evaluation | interpn | 141 ns | 104 ns | 1.35x | 1e-16 |
| `grid/linear/4d/max200/eval10000/ascending` | evaluation | interpn | 508 µs | 375 µs | 1.35x | 4e-16 |
| `oned/linear/n100/shuffled` | construction | scirs2-interpolate | 243 ns | 180 ns | 1.35x | 0e+00 |
| `grid/linear/3d/max200/eval1000/ascending` | construction | scirs2-interpolate | 1.74 ms | 1.29 ms | 1.35x | 3e-16 |
| `grid/linear/3d/max10/eval10000/ascending` | evaluation | interpn | 196 µs | 146 µs | 1.34x | 4e-16 |
| `grid/linear/2d/max10/eval10000/ascending` | evaluation | interpn | 133 µs | 101 µs | 1.32x | 3e-16 |
| `grid/linear/3d/max10/eval100/ascending` | evaluation | interpn | 1.98 µs | 1.51 µs | 1.31x | 3e-16 |
| `grid/linear/3d/max10/eval1000/ascending` | evaluation | interpn | 19.4 µs | 14.9 µs | 1.31x | 5e-16 |
| `grid/linear/2d/max200/eval100/flipped` | construction | SciPy | 15.4 µs | 12 µs | 1.29x |  |
| `oned/linear/n50/shuffled` | construction | scirs2-interpolate | 163 ns | 127 ns | 1.29x | 0e+00 |
| `grid/linear/2d/max200/eval10/ascending` | construction | SciPy | 7.54 µs | 5.87 µs | 1.28x |  |
| `grid/linear/4d/max200/eval100/ascending` | evaluation | interpn | 5.05 µs | 3.95 µs | 1.28x | 3e-16 |
| `grid/linear/3d/max40/eval10000/ascending` | construction | scirs2-interpolate | 15.5 µs | 12.2 µs | 1.27x | 2e-16 |
| `grid/linear/4d/max200/eval1000/ascending` | evaluation | interpn | 50.6 µs | 39.9 µs | 1.27x | 3e-16 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | SciPy | 7.44 µs | 5.89 µs | 1.26x |  |
| `grid/linear/2d/each800/eval10000/ascending` | construction | scirs2-interpolate | 1.02 ms | 812 µs | 1.26x | 2e-16 |
| `grid/linear/2d/max10/eval100/ascending` | evaluation | interpn | 1.28 µs | 1.02 µs | 1.25x | 2e-16 |
| `grid/linear/2d/max10/eval1000/ascending` | evaluation | interpn | 12.7 µs | 10.1 µs | 1.25x | 3e-16 |
| `grid/linear/2d/max200/eval10/flipped` | construction | SciPy | 14.6 µs | 11.9 µs | 1.23x |  |
| `grid/linear/2d/max40/eval10/ascending` | evaluation | interpn | 179 ns | 146 ns | 1.22x | 2e-16 |
| `grid/linear/3d/max200/eval100/ascending` | construction | scirs2-interpolate | 1.74 ms | 1.44 ms | 1.21x | 3e-16 |
| `grid/linear/4d/max200/eval100/ascending` | construction | scirs2-interpolate | 48.8 ms | 40.7 ms | 1.20x | 1e-16 |
| `grid/linear/3d/max200/eval10/ascending` | construction | scirs2-interpolate | 1.67 ms | 1.39 ms | 1.20x | 3e-16 |
| `grid/linear/4d/max200/eval10000/ascending` | construction | scirs2-interpolate | 47.5 ms | 39.8 ms | 1.20x | 4e-16 |
| `grid/linear/4d/max200/eval1000/ascending` | construction | scirs2-interpolate | 49.8 ms | 41.7 ms | 1.20x | 3e-16 |
| `grid/linear/2d/max40/eval1000/ascending` | evaluation | interpn | 15.6 µs | 13.1 µs | 1.19x | 3e-16 |
| `grid/linear/2d/max40/eval100/ascending` | evaluation | interpn | 1.72 µs | 1.44 µs | 1.19x | 2e-16 |
| `grid/linear/4d/max200/eval10/ascending` | construction | scirs2-interpolate | 47.2 ms | 39.7 ms | 1.19x | 2e-16 |
| `grid/linear/2d/max200/eval1000/flipped` | construction | SciPy | 14.7 µs | 12.5 µs | 1.18x |  |
| `grid/linear/2d/max200/eval1000/ascending` | evaluation | interpn | 20.3 µs | 17.2 µs | 1.18x | 2e-16 |
| `oned/nearest/n10000/sorted` | construction | SciPy | 61.2 µs | 52 µs | 1.18x |  |
| `grid/linear/2d/max200/eval100/ascending` | construction | SciPy | 7.43 µs | 6.37 µs | 1.17x |  |
| `grid/linear/2d/max200/eval100/ascending` | evaluation | interpn | 2.1 µs | 1.8 µs | 1.16x | 1e-16 |
| `grid/linear/3d/max40/eval1000/ascending` | evaluation | interpn | 22.9 µs | 19.7 µs | 1.16x | 2e-16 |
| `grid/linear/3d/max40/eval100/ascending` | evaluation | interpn | 2.5 µs | 2.16 µs | 1.16x | 2e-16 |
| `grid/linear/2d/max40/eval10000/ascending` | evaluation | interpn | 157 µs | 136 µs | 1.15x | 4e-16 |
| `grid/linear/3d/max40/eval10000/ascending` | evaluation | interpn | 230 µs | 201 µs | 1.14x | 3e-16 |
| `oned/nearest/n10000/shuffled` | construction | SciPy | 61.1 µs | 53.8 µs | 1.14x |  |
| `grid/linear/3d/max10/eval100/ascending` | construction | scirs2-interpolate | 194 ns | 171 ns | 1.13x | 2e-16 |
| `grid/linear/3d/max200/eval1000/ascending` | evaluation | interpn | 30.7 µs | 27.2 µs | 1.13x | 3e-16 |
| `grid/nearest/2d/each40/eval10000/ascending` | evaluation | interpn | 159 µs | 141 µs | 1.12x | 0e+00 |
| `grid/linear/3d/max10/eval10000/ascending` | construction | scirs2-interpolate | 197 ns | 175 ns | 1.12x | 3e-16 |
| `rbf/quintic/global/n10/scattered` | evaluation | SciPy | 75 µs | 66.8 µs | 1.12x |  |
| `grid/linear/2d/max200/eval10/ascending` | evaluation | interpn | 214 ns | 191 ns | 1.12x | 2e-16 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | scirs2-interpolate | 302 µs | 271 µs | 1.11x | 2e-16 |
| `grid/linear/2d/max200/eval10000/ascending` | evaluation | interpn | 209 µs | 189 µs | 1.11x | 3e-16 |
| `oned/bspline1/n10000/sorted` | evaluation | SciPy | 314 µs | 285 µs | 1.10x |  |

</details>

## Runs that failed

| case | library | how |
| :--- | :--- | :--- |
| `grid/linear/4d/max10/eval10/ascending` | interpn | exit code 101 |
| `grid/linear/4d/max10/eval10/ascending` | interpx | exit code 101 |
| `grid/linear/4d/max10/eval10/ascending` | interpx, `from_owned` | exit code 101 |
| `grid/linear/4d/max10/eval10/ascending` | scirs2-interpolate | exit code 101 |
| `grid/linear/4d/max10/eval10/flipped` | interpx | exit code 101 |
| `grid/linear/4d/max10/eval10/flipped` | interpx, `from_owned` | exit code 101 |
| `grid/linear/4d/max10/eval100/ascending` | interpn | exit code 101 |
| `grid/linear/4d/max10/eval100/ascending` | interpx | exit code 101 |
| `grid/linear/4d/max10/eval100/ascending` | interpx, `from_owned` | exit code 101 |
| `grid/linear/4d/max10/eval100/ascending` | scirs2-interpolate | exit code 101 |
| `grid/linear/4d/max10/eval100/flipped` | interpx | exit code 101 |
| `grid/linear/4d/max10/eval100/flipped` | interpx, `from_owned` | exit code 101 |
| `grid/linear/4d/max10/eval1000/ascending` | interpn | exit code 101 |
| `grid/linear/4d/max10/eval1000/ascending` | interpx | exit code 101 |
| `grid/linear/4d/max10/eval1000/ascending` | interpx, `from_owned` | exit code 101 |
| `grid/linear/4d/max10/eval1000/ascending` | scirs2-interpolate | exit code 101 |
| `grid/linear/4d/max10/eval1000/flipped` | interpx | exit code 101 |
| `grid/linear/4d/max10/eval1000/flipped` | interpx, `from_owned` | exit code 101 |
| `grid/linear/4d/max10/eval10000/ascending` | interpn | exit code 101 |
| `grid/linear/4d/max10/eval10000/ascending` | interpx | exit code 101 |
| `grid/linear/4d/max10/eval10000/ascending` | interpx, `from_owned` | exit code 101 |
| `grid/linear/4d/max10/eval10000/ascending` | scirs2-interpolate | exit code 101 |
| `grid/linear/4d/max10/eval10000/flipped` | interpx | exit code 101 |
| `grid/linear/4d/max10/eval10000/flipped` | interpx, `from_owned` | exit code 101 |

## One variable

Data of `n` points, evaluated at `4 n` points: in order (`sorted`) and in no order (`shuffled`).

### Evaluation

| case | interpx | SciPy | scirs2-interpolate | ndarray-interp | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: |
| `linear/n10/sorted` | 123 ns | 4.99 µs | 266 ns | 279 ns | 2.2x faster |
| `linear/n10/shuffled` | 184 ns | 5.03 µs | 285 ns | 274 ns | 1.5x faster |
| `nearest/n10/sorted` | 159 ns | 9.4 µs | 258 ns |  | 1.6x faster |
| `nearest/n10/shuffled` | 302 ns | 9.46 µs | 279 ns |  | **1.1x slower** |
| `previous/n10/sorted` | 154 ns | 10.6 µs |  |  | 69x faster |
| `previous/n10/shuffled` | 268 ns | 10.7 µs |  |  | 40x faster |
| `bspline1/n10/sorted` | 296 ns | 2.74 µs |  |  | 9.3x faster |
| `bspline1/n10/shuffled` | 600 ns | 2.69 µs |  |  | 4.5x faster |
| `bspline2/n10/sorted` | 440 ns | 3.06 µs |  |  | 7x faster |
| `bspline2/n10/shuffled` | 477 ns | 3.32 µs |  |  | 7x faster |
| `bspline3/n10/sorted` | 521 ns | 3.66 µs |  |  | 7x faster |
| `bspline3/n10/shuffled` | 605 ns | 3.6 µs |  |  | 6x faster |
| `cubic_spline/n10/sorted` | 140 ns | 3.64 µs | 268 ns | 337 ns | 1.9x faster |
| `cubic_spline/n10/shuffled` | 206 ns | 3.7 µs | 280 ns | 334 ns | 1.4x faster |
| `pchip/n10/sorted` | 138 ns | 3.62 µs | 295 ns |  | 2.1x faster |
| `pchip/n10/shuffled` | 208 ns | 3.65 µs | 291 ns |  | 1.4x faster |
| `akima/n10/sorted` | 141 ns | 3.84 µs | 252 ns |  | 1.8x faster |
| `akima/n10/shuffled` | 207 ns | 3.67 µs | 261 ns |  | 1.3x faster |
| `makima/n10/sorted` | 141 ns | 3.63 µs |  |  | 26x faster |
| `makima/n10/shuffled` | 207 ns | 3.63 µs |  |  | 18x faster |
| `linear/n50/sorted` | 675 ns | 6.65 µs | 1.67 µs | 1.35 µs | 2x faster |
| `linear/n50/shuffled` | 824 ns | 8.83 µs | 1.65 µs | 1.39 µs | 1.7x faster |
| `nearest/n50/sorted` | 861 ns | 12.2 µs | 1.67 µs |  | 1.9x faster |
| `nearest/n50/shuffled` | 1.58 µs | 12.3 µs | 1.66 µs |  | 1x faster |
| `previous/n50/sorted` | 739 ns | 13.6 µs |  |  | 18x faster |
| `previous/n50/shuffled` | 1.39 µs | 12.8 µs |  |  | 9.2x faster |
| `bspline1/n50/sorted` | 1.61 µs | 4.27 µs |  |  | 2.7x faster |
| `bspline1/n50/shuffled` | 3.4 µs | 5.29 µs |  |  | 1.6x faster |
| `bspline2/n50/sorted` | 2.66 µs | 6.33 µs |  |  | 2.4x faster |
| `bspline2/n50/shuffled` | 2.78 µs | 6.64 µs |  |  | 2.4x faster |
| `bspline3/n50/sorted` | 2.61 µs | 8.72 µs |  |  | 3.3x faster |
| `bspline3/n50/shuffled` | 3.73 µs | 9.3 µs |  |  | 2.5x faster |
| `cubic_spline/n50/sorted` | 688 ns | 5.43 µs | 4.53 µs | 1.64 µs | 2.4x faster |
| `cubic_spline/n50/shuffled` | 1.09 µs | 6 µs | 4.12 µs | 1.7 µs | 1.6x faster |
| `pchip/n50/sorted` | 698 ns | 5.74 µs | 3.9 µs |  | 5.6x faster |
| `pchip/n50/shuffled` | 988 ns | 5.49 µs | 3.78 µs |  | 3.8x faster |
| `akima/n50/sorted` | 687 ns | 5.18 µs | 4.41 µs |  | 6.4x faster |
| `akima/n50/shuffled` | 996 ns | 5.96 µs | 4.25 µs |  | 4.3x faster |
| `makima/n50/sorted` | 687 ns | 5.3 µs |  |  | 7.7x faster |
| `makima/n50/shuffled` | 1.03 µs | 5.96 µs |  |  | 5.8x faster |
| `linear/n100/sorted` | 1.3 µs | 7.29 µs | 3.6 µs | 2.52 µs | 1.9x faster |
| `linear/n100/shuffled` | 1.64 µs | 13.3 µs | 3.37 µs | 2.5 µs | 1.5x faster |
| `nearest/n100/sorted` | 1.55 µs | 13.4 µs | 3.34 µs |  | 2.2x faster |
| `nearest/n100/shuffled` | 3.48 µs | 13.7 µs | 3.3 µs |  | **1.1x slower** |
| `previous/n100/sorted` | 1.5 µs | 15 µs |  |  | 10x faster |
| `previous/n100/shuffled` | 3.04 µs | 15.3 µs |  |  | 5.1x faster |
| `bspline1/n100/sorted` | 3.2 µs | 5.37 µs |  |  | 1.7x faster |
| `bspline1/n100/shuffled` | 6.23 µs | 12.9 µs |  |  | 2.1x faster |
| `bspline2/n100/sorted` | 3.9 µs | 9.84 µs |  |  | 2.5x faster |
| `bspline2/n100/shuffled` | 5.96 µs | 12.9 µs |  |  | 2.2x faster |
| `bspline3/n100/sorted` | 5.68 µs | 15 µs |  |  | 2.6x faster |
| `bspline3/n100/shuffled` | 7.94 µs | 19.5 µs |  |  | 2.5x faster |
| `cubic_spline/n100/sorted` | 1.48 µs | 6.86 µs | 18.7 µs | 3.56 µs | 2.4x faster |
| `cubic_spline/n100/shuffled` | 2.16 µs | 8.73 µs | 19 µs | 3.27 µs | 1.5x faster |
| `pchip/n100/sorted` | 1.43 µs | 7.27 µs | 17.1 µs |  | 5.1x faster |
| `pchip/n100/shuffled` | 1.98 µs | 8.64 µs | 17.6 µs |  | 4.4x faster |
| `akima/n100/sorted` | 1.52 µs | 7.6 µs | 16.7 µs |  | 5x faster |
| `akima/n100/shuffled` | 2.19 µs | 8.94 µs | 19.6 µs |  | 4.1x faster |
| `makima/n100/sorted` | 1.54 µs | 7.56 µs |  |  | 4.9x faster |
| `makima/n100/shuffled` | 2 µs | 8.66 µs |  |  | 4.3x faster |
| `linear/n1000/sorted` | 14 µs | 32.1 µs | 120 µs | 27.1 µs | 1.9x faster |
| `linear/n1000/shuffled` | 18 µs | 156 µs | 231 µs | 27.1 µs | 1.5x faster |
| `nearest/n1000/sorted` | 27.8 µs | 101 µs | 197 µs |  | 3.6x faster |
| `nearest/n1000/shuffled` | 37.2 µs | 331 µs | 263 µs |  | 7.1x faster |
| `previous/n1000/sorted` | 16.2 µs | 63.1 µs |  |  | 3.9x faster |
| `previous/n1000/shuffled` | 30.4 µs | 225 µs |  |  | 7.4x faster |
| `bspline1/n1000/sorted` | 33.5 µs | 34.1 µs |  |  | 1x faster |
| `bspline1/n1000/shuffled` | 62.6 µs | 583 µs |  |  | 9.3x faster |
| `bspline2/n1000/sorted` | 44.7 µs | 84.4 µs |  |  | 1.9x faster |
| `bspline2/n1000/shuffled` | 72.2 µs | 543 µs |  |  | 7.5x faster |
| `bspline3/n1000/sorted` | 58.2 µs | 129 µs |  |  | 2.2x faster |
| `bspline3/n1000/shuffled` | 95.3 µs | 594 µs |  |  | 6.2x faster |
| `cubic_spline/n1000/sorted` | 14 µs | 41 µs | 1.56 ms | 32.7 µs | 2.3x faster |
| `cubic_spline/n1000/shuffled` | 20 µs | 220 µs | 1.43 ms | 38.5 µs | 1.9x faster |
| `pchip/n1000/sorted` | 14.5 µs | 37.8 µs | 1.32 ms |  | 2.6x faster |
| `pchip/n1000/shuffled` | 20.1 µs | 227 µs | 1.21 ms |  | 11x faster |
| `akima/n1000/sorted` | 13.6 µs | 37.3 µs | 1.42 ms |  | 2.7x faster |
| `akima/n1000/shuffled` | 19.9 µs | 223 µs | 1.42 ms |  | 11x faster |
| `makima/n1000/sorted` | 13.9 µs | 38.4 µs |  |  | 2.8x faster |
| `makima/n1000/shuffled` | 21.8 µs | 220 µs |  |  | 10x faster |
| `linear/n10000/sorted` | 147 µs | 231 µs | 1.75 ms | 272 µs | 1.6x faster |
| `linear/n10000/shuffled` | 191 µs | 2.35 ms | 3.15 ms | 260 µs | 1.4x faster |
| `nearest/n10000/sorted` | 157 µs | 584 µs | 2.13 ms |  | 3.7x faster |
| `nearest/n10000/shuffled` | 389 µs | 2.93 ms | 3.1 ms |  | 7.5x faster |
| `previous/n10000/sorted` | 163 µs | 608 µs |  |  | 3.7x faster |
| `previous/n10000/shuffled` | 354 µs | 3.01 ms |  |  | 8.5x faster |
| `bspline1/n10000/sorted` | 314 µs | 285 µs |  |  | **1.1x slower** |
| `bspline1/n10000/shuffled` | 706 µs | 51.9 ms |  |  | 74x faster |
| `bspline2/n10000/sorted` | 457 µs | 769 µs |  |  | 1.7x faster |
| `bspline2/n10000/shuffled` | 1.09 ms | 51.2 ms |  |  | 47x faster |
| `bspline3/n10000/sorted` | 613 µs | 1.27 ms |  |  | 2.1x faster |
| `bspline3/n10000/shuffled` | 1.42 ms | 54.7 ms |  |  | 38x faster |
| `cubic_spline/n10000/sorted` | 149 µs | 373 µs | 185 ms | 360 µs | 2.4x faster |
| `cubic_spline/n10000/shuffled` | 223 µs | 3.92 ms | 177 ms | 712 µs | 3.2x faster |
| `pchip/n10000/sorted` | 158 µs | 405 µs | 127 ms |  | 2.6x faster |
| `pchip/n10000/shuffled` | 240 µs | 3.72 ms | 130 ms |  | 15x faster |
| `akima/n10000/sorted` | 160 µs | 408 µs | 161 ms |  | 2.5x faster |
| `akima/n10000/shuffled` | 241 µs | 3.71 ms | 155 ms |  | 15x faster |
| `makima/n10000/sorted` | 160 µs | 399 µs |  |  | 2.5x faster |
| `makima/n10000/shuffled` | 240 µs | 3.68 ms |  |  | 15x faster |
| `linear/n100/sorted-uneven` | 1.26 µs | 8.07 µs | 3.61 µs | 3.78 µs | 2.9x faster |
| `linear/n100/shuffled-uneven` | 2.95 µs | 12.6 µs | 3.62 µs | 4.04 µs | 1.2x faster |
| `cubic_spline/n100/sorted-uneven` | 1.39 µs | 7.4 µs | 17.4 µs | 4.93 µs | 3.6x faster |
| `cubic_spline/n100/shuffled-uneven` | 3.04 µs | 7.96 µs | 19.4 µs | 4.47 µs | 1.5x faster |
| `pchip/n100/sorted-uneven` | 1.47 µs | 7.38 µs | 15.7 µs |  | 5x faster |
| `pchip/n100/shuffled-uneven` | 3.03 µs | 7.91 µs | 16 µs |  | 2.6x faster |
| `linear/n10000/sorted-uneven` | 223 µs | 262 µs | 1.65 ms | 980 µs | 1.2x faster |
| `linear/n10000/shuffled-uneven` | 609 µs | 2.34 ms | 3.09 ms | 1.92 ms | 3.1x faster |
| `cubic_spline/n10000/sorted-uneven` | 217 µs | 458 µs | 148 ms | 1.06 ms | 2.1x faster |
| `cubic_spline/n10000/shuffled-uneven` | 764 µs | 3.62 ms | 145 ms | 1.86 ms | 2.4x faster |
| `pchip/n10000/sorted-uneven` | 234 µs | 482 µs | 120 ms |  | 2.1x faster |
| `pchip/n10000/shuffled-uneven` | 760 µs | 3.37 ms | 129 ms |  | 4.4x faster |

### Construction

| case | interpx | SciPy | scirs2-interpolate | ndarray-interp | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: |
| `linear/n10/sorted` | 94.9 ns | 9.42 µs | 88 ns | 43.5 ns | **2.2x slower** |
| `linear/n10/shuffled` | 94.5 ns | 8.7 µs | 86.9 ns | 47.1 ns | **2x slower** |
| `nearest/n10/sorted` | 64.2 ns | 9.52 µs | 80.5 ns |  | 1.3x faster |
| `nearest/n10/shuffled` | 64.2 ns | 9.47 µs | 80.5 ns |  | 1.3x faster |
| `previous/n10/sorted` | 70.7 ns | 9.49 µs |  |  | 1.3e+02x faster |
| `previous/n10/shuffled` | 64.6 ns | 9.32 µs |  |  | 1.4e+02x faster |
| `bspline1/n10/sorted` | 115 ns | 17 µs |  |  | 1.5e+02x faster |
| `bspline1/n10/shuffled` | 120 ns | 15.5 µs |  |  | 1.3e+02x faster |
| `bspline2/n10/sorted` | 453 ns | 36.2 µs |  |  | 80x faster |
| `bspline2/n10/shuffled` | 452 ns | 39.1 µs |  |  | 87x faster |
| `bspline3/n10/sorted` | 503 ns | 35.9 µs |  |  | 71x faster |
| `bspline3/n10/shuffled` | 505 ns | 35.8 µs |  |  | 71x faster |
| `cubic_spline/n10/sorted` | 253 ns | 85.6 µs | 377 ns | 471 ns | 1.5x faster |
| `cubic_spline/n10/shuffled` | 233 ns | 88.4 µs | 383 ns | 470 ns | 1.6x faster |
| `pchip/n10/sorted` | 175 ns | 101 µs | 193 ns |  | 1.1x faster |
| `pchip/n10/shuffled` | 191 ns | 100 µs | 200 ns |  | 1x faster |
| `akima/n10/sorted` | 191 ns | 74.5 µs | 190 ns |  | **1x slower** |
| `akima/n10/shuffled` | 186 ns | 74.7 µs | 184 ns |  | **1x slower** |
| `makima/n10/sorted` | 208 ns | 76.5 µs |  |  | 3.7e+02x faster |
| `makima/n10/shuffled` | 209 ns | 77.7 µs |  |  | 3.7e+02x faster |
| `linear/n50/sorted` | 163 ns | 9.77 µs | 116 ns | 116 ns | **1.4x slower** |
| `linear/n50/shuffled` | 163 ns | 9.35 µs | 127 ns | 115 ns | **1.4x slower** |
| `nearest/n50/sorted` | 107 ns | 10.8 µs | 127 ns |  | 1.2x faster |
| `nearest/n50/shuffled` | 107 ns | 9.74 µs | 118 ns |  | 1.1x faster |
| `previous/n50/sorted` | 98.2 ns | 10.5 µs |  |  | 1.1e+02x faster |
| `previous/n50/shuffled` | 98.9 ns | 9.62 µs |  |  | 97x faster |
| `bspline1/n50/sorted` | 203 ns | 17.2 µs |  |  | 85x faster |
| `bspline1/n50/shuffled` | 190 ns | 15.9 µs |  |  | 84x faster |
| `bspline2/n50/sorted` | 1.98 µs | 42.1 µs |  |  | 21x faster |
| `bspline2/n50/shuffled` | 1.98 µs | 40.6 µs |  |  | 20x faster |
| `bspline3/n50/sorted` | 2.08 µs | 42.8 µs |  |  | 21x faster |
| `bspline3/n50/shuffled` | 2.2 µs | 41.9 µs |  |  | 19x faster |
| `cubic_spline/n50/sorted` | 834 ns | 93.7 µs | 1 µs | 1.07 µs | 1.2x faster |
| `cubic_spline/n50/shuffled` | 832 ns | 90.8 µs | 992 ns | 1.19 µs | 1.2x faster |
| `pchip/n50/sorted` | 530 ns | 110 µs | 365 ns |  | **1.5x slower** |
| `pchip/n50/shuffled` | 522 ns | 104 µs | 370 ns |  | **1.4x slower** |
| `akima/n50/sorted` | 571 ns | 75.8 µs | 360 ns |  | **1.6x slower** |
| `akima/n50/shuffled` | 574 ns | 82.8 µs | 329 ns |  | **1.7x slower** |
| `makima/n50/sorted` | 598 ns | 84.6 µs |  |  | 1.4e+02x faster |
| `makima/n50/shuffled` | 602 ns | 85 µs |  |  | 1.4e+02x faster |
| `linear/n100/sorted` | 265 ns | 9.06 µs | 180 ns | 187 ns | **1.5x slower** |
| `linear/n100/shuffled` | 243 ns | 9.19 µs | 180 ns | 171 ns | **1.4x slower** |
| `nearest/n100/sorted` | 154 ns | 9.99 µs | 196 ns |  | 1.3x faster |
| `nearest/n100/shuffled` | 167 ns | 9.96 µs | 185 ns |  | 1.1x faster |
| `previous/n100/sorted` | 153 ns | 10.1 µs |  |  | 66x faster |
| `previous/n100/shuffled` | 155 ns | 10.9 µs |  |  | 70x faster |
| `bspline1/n100/sorted` | 311 ns | 15.6 µs |  |  | 50x faster |
| `bspline1/n100/shuffled` | 294 ns | 16.3 µs |  |  | 55x faster |
| `bspline2/n100/sorted` | 3.49 µs | 43.9 µs |  |  | 13x faster |
| `bspline2/n100/shuffled` | 3.8 µs | 43.7 µs |  |  | 11x faster |
| `bspline3/n100/sorted` | 4.26 µs | 47.1 µs |  |  | 11x faster |
| `bspline3/n100/shuffled` | 3.89 µs | 47.6 µs |  |  | 12x faster |
| `cubic_spline/n100/sorted` | 1.58 µs | 90 µs | 1.81 µs | 2.05 µs | 1.2x faster |
| `cubic_spline/n100/shuffled` | 1.73 µs | 95.7 µs | 1.81 µs | 1.87 µs | 1x faster |
| `pchip/n100/sorted` | 956 ns | 103 µs | 554 ns |  | **1.7x slower** |
| `pchip/n100/shuffled` | 873 ns | 113 µs | 605 ns |  | **1.4x slower** |
| `akima/n100/sorted` | 1.03 µs | 83.8 µs | 602 ns |  | **1.7x slower** |
| `akima/n100/shuffled` | 1.04 µs | 84.6 µs | 613 ns |  | **1.7x slower** |
| `makima/n100/sorted` | 1.08 µs | 89.3 µs |  |  | 82x faster |
| `makima/n100/shuffled` | 1.08 µs | 79.9 µs |  |  | 74x faster |
| `linear/n1000/sorted` | 2.12 µs | 14.2 µs | 1.17 µs | 1.52 µs | **1.8x slower** |
| `linear/n1000/shuffled` | 2.12 µs | 14.3 µs | 1.16 µs | 1.52 µs | **1.8x slower** |
| `nearest/n1000/sorted` | 1.76 µs | 16.5 µs | 1.82 µs |  | 1x faster |
| `nearest/n1000/shuffled` | 1.62 µs | 26.2 µs | 1.69 µs |  | 1x faster |
| `previous/n1000/sorted` | 1.09 µs | 16.7 µs |  |  | 15x faster |
| `previous/n1000/shuffled` | 1.1 µs | 18.7 µs |  |  | 17x faster |
| `bspline1/n1000/sorted` | 2.04 µs | 18.8 µs |  |  | 9.2x faster |
| `bspline1/n1000/shuffled` | 2.02 µs | 17.9 µs |  |  | 8.8x faster |
| `bspline2/n1000/sorted` | 37.5 µs | 120 µs |  |  | 3.2x faster |
| `bspline2/n1000/shuffled` | 34.4 µs | 119 µs |  |  | 3.5x faster |
| `bspline3/n1000/sorted` | 72.6 µs | 126 µs |  |  | 1.7x faster |
| `bspline3/n1000/shuffled` | 72.1 µs | 126 µs |  |  | 1.7x faster |
| `cubic_spline/n1000/sorted` | 40.6 µs | 113 µs | 20.3 µs | 16.1 µs | **2.5x slower** |
| `cubic_spline/n1000/shuffled` | 39.8 µs | 113 µs | 18.3 µs | 15.8 µs | **2.5x slower** |
| `pchip/n1000/sorted` | 11.8 µs | 121 µs | 4.58 µs |  | **2.6x slower** |
| `pchip/n1000/shuffled` | 11.8 µs | 126 µs | 4.59 µs |  | **2.6x slower** |
| `akima/n1000/sorted` | 12 µs | 94.9 µs | 6.34 µs |  | **1.9x slower** |
| `akima/n1000/shuffled` | 12.7 µs | 97.3 µs | 6.34 µs |  | **2x slower** |
| `makima/n1000/sorted` | 12.8 µs | 106 µs |  |  | 8.3x faster |
| `makima/n1000/shuffled` | 14.2 µs | 107 µs |  |  | 7.5x faster |
| `linear/n10000/sorted` | 99.4 µs | 44.7 µs | 62 µs | 21.3 µs | **4.7x slower** |
| `linear/n10000/shuffled` | 96.4 µs | 46.4 µs | 59.9 µs | 32.2 µs | **3x slower** |
| `nearest/n10000/sorted` | 61.2 µs | 52 µs | 60.2 µs |  | **1.2x slower** |
| `nearest/n10000/shuffled` | 61.1 µs | 53.8 µs | 56.4 µs |  | **1.1x slower** |
| `previous/n10000/sorted` | 65.4 µs | 74.8 µs |  |  | 1.1x faster |
| `previous/n10000/shuffled` | 61.3 µs | 75.3 µs |  |  | 1.2x faster |
| `bspline1/n10000/sorted` | 75.5 µs | 23.1 µs |  |  | **3.3x slower** |
| `bspline1/n10000/shuffled` | 75.8 µs | 25.4 µs |  |  | **3x slower** |
| `bspline2/n10000/sorted` | 525 µs | 801 µs |  |  | 1.5x faster |
| `bspline2/n10000/shuffled` | 530 µs | 763 µs |  |  | 1.4x faster |
| `bspline3/n10000/sorted` | 633 µs | 960 µs |  |  | 1.5x faster |
| `bspline3/n10000/shuffled` | 672 µs | 1.06 ms |  |  | 1.6x faster |
| `cubic_spline/n10000/sorted` | 249 µs | 778 µs | 376 µs | 470 µs | 1.5x faster |
| `cubic_spline/n10000/shuffled` | 249 µs | 796 µs | 306 µs | 397 µs | 1.2x faster |
| `pchip/n10000/sorted` | 136 µs | 261 µs | 175 µs |  | 1.3x faster |
| `pchip/n10000/shuffled` | 129 µs | 262 µs | 187 µs |  | 1.4x faster |
| `akima/n10000/sorted` | 101 µs | 710 µs | 167 µs |  | 1.7x faster |
| `akima/n10000/shuffled` | 100 µs | 870 µs | 137 µs |  | 1.4x faster |
| `makima/n10000/sorted` | 153 µs | 798 µs |  |  | 5.2x faster |
| `makima/n10000/shuffled` | 102 µs | 797 µs |  |  | 7.8x faster |
| `linear/n100/sorted-uneven` | 266 ns | 9.07 µs | 194 ns | 170 ns | **1.6x slower** |
| `linear/n100/shuffled-uneven` | 268 ns | 10 µs | 197 ns | 187 ns | **1.4x slower** |
| `cubic_spline/n100/sorted-uneven` | 1.73 µs | 88.3 µs | 1.8 µs | 2.05 µs | 1x faster |
| `cubic_spline/n100/shuffled-uneven` | 1.57 µs | 87.5 µs | 1.67 µs | 2.08 µs | 1.1x faster |
| `pchip/n100/sorted-uneven` | 953 ns | 102 µs | 602 ns |  | **1.6x slower** |
| `pchip/n100/shuffled-uneven` | 936 ns | 105 µs | 625 ns |  | **1.5x slower** |
| `linear/n10000/sorted-uneven` | 97.6 µs | 47 µs | 59.5 µs | 14.5 µs | **6.8x slower** |
| `linear/n10000/shuffled-uneven` | 94.7 µs | 45.6 µs | 60.8 µs | 33.1 µs | **2.9x slower** |
| `cubic_spline/n10000/sorted-uneven` | 249 µs | 769 µs | 334 µs | 382 µs | 1.3x faster |
| `cubic_spline/n10000/shuffled-uneven` | 250 µs | 735 µs | 360 µs | 383 µs | 1.4x faster |
| `pchip/n10000/sorted-uneven` | 102 µs | 260 µs | 171 µs |  | 1.7x faster |
| `pchip/n10000/shuffled-uneven` | 121 µs | 261 µs | 169 µs |  | 1.4x faster |

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
| `linear/2d/max10/eval10/ascending` | 141 ns | 141 ns | 21.4 µs | 784 ns | 156 ns | 104 ns | **1.4x slower** |
| `linear/2d/max10/eval10/flipped` | 129 ns | 128 ns | 19.3 µs |  |  |  | 1.5e+02x faster |
| `linear/2d/max10/eval100/ascending` | 1.28 µs | 1.31 µs | 22.7 µs | 9.07 µs | 1.48 µs | 1.02 µs | **1.3x slower** |
| `linear/2d/max10/eval100/flipped` | 1.28 µs | 1.4 µs | 21.6 µs |  |  |  | 17x faster |
| `linear/2d/max10/eval1000/ascending` | 12.7 µs | 13.9 µs | 44.9 µs | 89.1 µs | 16 µs | 10.1 µs | **1.3x slower** |
| `linear/2d/max10/eval1000/flipped` | 13.9 µs | 13.9 µs | 44 µs |  |  |  | 3.2x faster |
| `linear/2d/max10/eval10000/ascending` | 133 µs | 139 µs | 359 µs | 916 µs | 243 µs | 101 µs | **1.3x slower** |
| `linear/2d/max10/eval10000/flipped` | 262 µs | 140 µs | 660 µs |  |  |  | 2.5x faster |
| `linear/2d/max40/eval10/ascending` | 179 ns | 178 ns | 21.2 µs | 1.1 µs | 172 ns | 146 ns | **1.2x slower** |
| `linear/2d/max40/eval10/flipped` | 177 ns | 178 ns | 21.1 µs |  |  |  | 1.2e+02x faster |
| `linear/2d/max40/eval100/ascending` | 1.72 µs | 1.72 µs | 24.5 µs | 11 µs | 1.59 µs | 1.44 µs | **1.2x slower** |
| `linear/2d/max40/eval100/flipped` | 1.73 µs | 1.58 µs | 24.1 µs |  |  |  | 14x faster |
| `linear/2d/max40/eval1000/ascending` | 15.6 µs | 15.6 µs | 50.2 µs | 114 µs | 17.3 µs | 13.1 µs | **1.2x slower** |
| `linear/2d/max40/eval1000/flipped` | 15.6 µs | 15.8 µs | 50.5 µs |  |  |  | 3.2x faster |
| `linear/2d/max40/eval10000/ascending` | 157 µs | 164 µs | 628 µs | 1.11 ms | 246 µs | 136 µs | **1.2x slower** |
| `linear/2d/max40/eval10000/flipped` | 157 µs | 157 µs | 574 µs |  |  |  | 3.7x faster |
| `linear/2d/max200/eval10/ascending` | 214 ns | 214 ns | 21.1 µs | 1.95 µs | 205 ns | 191 ns | **1.1x slower** |
| `linear/2d/max200/eval10/flipped` | 214 ns | 214 ns | 19.4 µs |  |  |  | 91x faster |
| `linear/2d/max200/eval100/ascending` | 2.1 µs | 2.1 µs | 24.7 µs | 20.6 µs | 2.01 µs | 1.8 µs | **1.2x slower** |
| `linear/2d/max200/eval100/flipped` | 1.94 µs | 1.91 µs | 22.8 µs |  |  |  | 12x faster |
| `linear/2d/max200/eval1000/ascending` | 20.3 µs | 20.9 µs | 73.4 µs | 208 µs | 24.6 µs | 17.2 µs | **1.2x slower** |
| `linear/2d/max200/eval1000/flipped` | 19.2 µs | 20.9 µs | 69.3 µs |  |  |  | 3.6x faster |
| `linear/2d/max200/eval10000/ascending` | 209 µs | 209 µs | 818 µs | 2.12 ms | 454 µs | 189 µs | **1.1x slower** |
| `linear/2d/max200/eval10000/flipped` | 210 µs | 210 µs | 840 µs |  |  |  | 4x faster |
| `linear/3d/max10/eval10/ascending` | 218 ns | 219 ns | 61.7 µs | 3.46 µs |  | 152 ns | **1.4x slower** |
| `linear/3d/max10/eval10/flipped` | 218 ns | 219 ns | 56.7 µs |  |  |  | 2.6e+02x faster |
| `linear/3d/max10/eval100/ascending` | 1.98 µs | 1.94 µs | 69.1 µs | 33.4 µs |  | 1.51 µs | **1.3x slower** |
| `linear/3d/max10/eval100/flipped` | 1.94 µs | 2.11 µs | 63.3 µs |  |  |  | 33x faster |
| `linear/3d/max10/eval1000/ascending` | 19.4 µs | 19.4 µs | 116 µs | 345 µs |  | 14.9 µs | **1.3x slower** |
| `linear/3d/max10/eval1000/flipped` | 19.6 µs | 19.4 µs | 116 µs |  |  |  | 5.9x faster |
| `linear/3d/max10/eval10000/ascending` | 196 µs | 195 µs | 1.2 ms | 3.48 ms |  | 146 µs | **1.3x slower** |
| `linear/3d/max10/eval10000/flipped` | 194 µs | 194 µs | 1.1 ms |  |  |  | 5.7x faster |
| `linear/3d/max40/eval10/ascending` | 235 ns | 235 ns | 56.3 µs | 3.38 µs |  | 218 ns | **1.1x slower** |
| `linear/3d/max40/eval10/flipped` | 234 ns | 234 ns | 57.9 µs |  |  |  | 2.5e+02x faster |
| `linear/3d/max40/eval100/ascending` | 2.5 µs | 2.51 µs | 63.2 µs | 37.5 µs |  | 2.16 µs | **1.2x slower** |
| `linear/3d/max40/eval100/flipped` | 2.51 µs | 2.52 µs | 66.5 µs |  |  |  | 26x faster |
| `linear/3d/max40/eval1000/ascending` | 22.9 µs | 22.9 µs | 147 µs | 351 µs |  | 19.7 µs | **1.2x slower** |
| `linear/3d/max40/eval1000/flipped` | 25.1 µs | 22.9 µs | 141 µs |  |  |  | 5.6x faster |
| `linear/3d/max40/eval10000/ascending` | 230 µs | 228 µs | 1.45 ms | 3.49 ms |  | 201 µs | **1.1x slower** |
| `linear/3d/max40/eval10000/flipped` | 229 µs | 230 µs | 1.47 ms |  |  |  | 6.4x faster |
| `linear/3d/max200/eval10/ascending` | 284 ns | 286 ns | 56.5 µs | 4.53 µs |  | 290 ns | 1x faster |
| `linear/3d/max200/eval10/flipped` | 311 ns | 284 ns | 58.6 µs |  |  |  | 1.9e+02x faster |
| `linear/3d/max200/eval100/ascending` | 2.81 µs | 2.81 µs | 64.9 µs | 51.1 µs |  | 2.68 µs | **1x slower** |
| `linear/3d/max200/eval100/flipped` | 2.82 µs | 2.8 µs | 64.3 µs |  |  |  | 23x faster |
| `linear/3d/max200/eval1000/ascending` | 30.7 µs | 28 µs | 184 µs | 509 µs |  | 27.2 µs | **1.1x slower** |
| `linear/3d/max200/eval1000/flipped` | 28.1 µs | 28 µs | 187 µs |  |  |  | 6.7x faster |
| `linear/3d/max200/eval10000/ascending` | 280 µs | 281 µs | 1.94 ms | 4.8 ms |  | 262 µs | **1.1x slower** |
| `linear/3d/max200/eval10000/flipped` | 286 µs | 280 µs | 1.88 ms |  |  |  | 6.6x faster |
| `linear/4d/max10/eval10/ascending` |  |  | 95.9 µs |  |  |  |  |
| `linear/4d/max10/eval10/flipped` |  |  | 95.7 µs |  |  |  |  |
| `linear/4d/max10/eval100/ascending` |  |  | 115 µs |  |  |  |  |
| `linear/4d/max10/eval100/flipped` |  |  | 116 µs |  |  |  |  |
| `linear/4d/max10/eval1000/ascending` |  |  | 224 µs |  |  |  |  |
| `linear/4d/max10/eval1000/flipped` |  |  | 227 µs |  |  |  |  |
| `linear/4d/max10/eval10000/ascending` |  |  | 1.88 ms |  |  |  |  |
| `linear/4d/max10/eval10000/flipped` |  |  | 1.81 ms |  |  |  |  |
| `linear/4d/max40/eval10/ascending` | 480 ns | 481 ns | 97.6 µs | 6.31 µs |  | 270 ns | **1.8x slower** |
| `linear/4d/max40/eval10/flipped` | 440 ns | 456 ns | 104 µs |  |  |  | 2.4e+02x faster |
| `linear/4d/max40/eval100/ascending` | 4.72 µs | 4.73 µs | 107 µs | 63.4 µs |  | 2.93 µs | **1.6x slower** |
| `linear/4d/max40/eval100/flipped` | 4.72 µs | 4.73 µs | 107 µs |  |  |  | 23x faster |
| `linear/4d/max40/eval1000/ascending` | 44.1 µs | 43.2 µs | 267 µs | 630 µs |  | 27.6 µs | **1.6x slower** |
| `linear/4d/max40/eval1000/flipped` | 43.2 µs | 47 µs | 248 µs |  |  |  | 5.7x faster |
| `linear/4d/max40/eval10000/ascending` | 431 µs | 472 µs | 2.47 ms | 6.52 ms |  | 294 µs | **1.5x slower** |
| `linear/4d/max40/eval10000/flipped` | 439 µs | 438 µs | 2.24 ms |  |  |  | 5.1x faster |
| `linear/4d/max200/eval10/ascending` | 515 ns | 516 ns | 95.5 µs | 7.92 µs |  | 359 ns | **1.4x slower** |
| `linear/4d/max200/eval10/flipped` | 565 ns | 560 ns | 108 µs |  |  |  | 1.9e+02x faster |
| `linear/4d/max200/eval100/ascending` | 5.05 µs | 5.56 µs | 110 µs | 72.9 µs |  | 3.95 µs | **1.3x slower** |
| `linear/4d/max200/eval100/flipped` | 5.16 µs | 5.37 µs | 111 µs |  |  |  | 22x faster |
| `linear/4d/max200/eval1000/ascending` | 50.6 µs | 50.6 µs | 327 µs | 750 µs |  | 39.9 µs | **1.3x slower** |
| `linear/4d/max200/eval1000/flipped` | 50.6 µs | 50.5 µs | 308 µs |  |  |  | 6.1x faster |
| `linear/4d/max200/eval10000/ascending` | 508 µs | 557 µs | 2.83 ms | 7.63 ms |  | 375 µs | **1.4x slower** |
| `linear/4d/max200/eval10000/flipped` | 524 µs | 506 µs | 2.78 ms |  |  |  | 5.3x faster |
| `nearest/2d/each40/eval10000/ascending` | 159 µs | 159 µs | 652 µs | 1.6 ms |  | 141 µs | **1.1x slower** |
| `nearest/3d/each40/eval10000/ascending` | 220 µs | 242 µs | 953 µs | 1.88 ms |  | 210 µs | **1x slower** |
| `pchip/2d/each40/eval10000/ascending` | 542 µs | 536 µs | 1.22 s |  |  |  | 2.2e+03x faster |
| `pchip/3d/each40/eval10000/ascending` | 1.59 ms | 1.58 ms | 3.3 s |  |  |  | 2.1e+03x faster |
| `cubic/2d/each40/eval10000/ascending` | 372 µs | 381 µs | 1.3 ms |  |  | 402 µs | 1.1x faster |
| `cubic/3d/each40/eval10000/ascending` | 813 µs | 847 µs | 2.97 ms |  |  | 803 µs | **1x slower** |
| `cubic/2d/each10/eval1000/ascending` | 34 µs | 33.9 µs | 129 µs |  |  | 35.1 µs | 1x faster |
| `cubic/2d/each40/eval1000/ascending` | 36.6 µs | 36.6 µs | 142 µs |  |  | 43.5 µs | 1.2x faster |
| `cubic/2d/each100/eval1000/ascending` | 38.3 µs | 38.5 µs | 166 µs |  |  | 42 µs | 1.1x faster |
| `cubic/2d/each200/eval1000/ascending` | 40.5 µs | 40.1 µs | 200 µs |  |  | 44.9 µs | 1.1x faster |
| `cubic/2d/each400/eval1000/ascending` | 42.6 µs | 42.6 µs | 275 µs |  |  | 48.5 µs | 1.1x faster |
| `quintic/2d/each10/eval1000/ascending` | 60.5 µs | 64.5 µs | 227 µs |  |  |  | 3.8x faster |
| `quintic/2d/each40/eval1000/ascending` | 65.3 µs | 65.5 µs | 237 µs |  |  |  | 3.6x faster |
| `linear/2d/each400/eval10000/ascending` | 220 µs | 221 µs | 980 µs | 4.2 ms | 545 µs | 217 µs | **1x slower** |
| `linear/2d/each800/eval10000/ascending` | 241 µs | 241 µs | 1.09 ms | 8.11 ms | 587 µs | 243 µs | 1x faster |

### Construction

| case | interpx | interpx, `from_owned` | SciPy | scirs2-interpolate | ndarray-interp | interpn | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/2d/max10/eval10/ascending` | 137 ns | 156 ns | 6.44 µs | 146 ns | 82.9 ns | 0 | **1.8e+02x slower** |
| `linear/2d/max10/eval10/flipped` | 170 ns | 189 ns | 12.7 µs |  |  |  | 75x faster |
| `linear/2d/max10/eval100/ascending` | 137 ns | 161 ns | 6.3 µs | 146 ns | 85.5 ns | 0 | **1.8e+02x slower** |
| `linear/2d/max10/eval100/flipped` | 156 ns | 187 ns | 12.6 µs |  |  |  | 81x faster |
| `linear/2d/max10/eval1000/ascending` | 137 ns | 150 ns | 6.55 µs | 144 ns | 85.8 ns | 0 | **1.8e+02x slower** |
| `linear/2d/max10/eval1000/flipped` | 171 ns | 182 ns | 12.7 µs |  |  |  | 74x faster |
| `linear/2d/max10/eval10000/ascending` | 135 ns | 137 ns | 6.29 µs | 137 ns | 85.4 ns | 0 | **2e+02x slower** |
| `linear/2d/max10/eval10000/flipped` | 320 ns | 211 ns | 22.4 µs |  |  |  | 70x faster |
| `linear/2d/max40/eval10/ascending` | 408 ns | 887 ns | 6.35 µs | 230 ns | 165 ns | 0 | **6e+02x slower** |
| `linear/2d/max40/eval10/flipped` | 702 ns | 1.14 µs | 12.8 µs |  |  |  | 18x faster |
| `linear/2d/max40/eval100/ascending` | 406 ns | 755 ns | 6.18 µs | 229 ns | 164 ns | 0 | **5.2e+02x slower** |
| `linear/2d/max40/eval100/flipped` | 1.28 µs | 1.22 µs | 12.8 µs |  |  |  | 10x faster |
| `linear/2d/max40/eval1000/ascending` | 371 ns | 821 ns | 6 µs | 224 ns | 164 ns | 0 | **5.4e+02x slower** |
| `linear/2d/max40/eval1000/flipped` | 643 ns | 1.18 µs | 12.8 µs |  |  |  | 20x faster |
| `linear/2d/max40/eval10000/ascending` | 370 ns | 787 ns | 5.88 µs | 224 ns | 150 ns | 0 | **5.4e+02x slower** |
| `linear/2d/max40/eval10000/flipped` | 698 ns | 1.18 µs | 11.8 µs |  |  |  | 17x faster |
| `linear/2d/max200/eval10/ascending` | 7.54 µs | 19.9 µs | 5.87 µs | 3.3 µs | 514 ns | 0 | **1.1e+04x slower** |
| `linear/2d/max200/eval10/flipped` | 14.6 µs | 29.3 µs | 11.9 µs |  |  |  | **1.2x slower** |
| `linear/2d/max200/eval100/ascending` | 7.43 µs | 17.9 µs | 6.37 µs | 3.16 µs | 513 ns | 0 | **1.1e+04x slower** |
| `linear/2d/max200/eval100/flipped` | 15.4 µs | 27.6 µs | 12 µs |  |  |  | **1.3x slower** |
| `linear/2d/max200/eval1000/ascending` | 7.44 µs | 18.8 µs | 5.89 µs | 3.15 µs | 514 ns | 0 | **1.1e+04x slower** |
| `linear/2d/max200/eval1000/flipped` | 14.7 µs | 27.6 µs | 12.5 µs |  |  |  | **1.2x slower** |
| `linear/2d/max200/eval10000/ascending` | 49.1 µs | 19.2 µs | 5.81 µs | 56.7 µs | 512 ns | 0 | **7.2e+04x slower** |
| `linear/2d/max200/eval10000/flipped` | 58.7 µs | 30.7 µs | 12.9 µs |  |  |  | **4.6x slower** |
| `linear/3d/max10/eval10/ascending` | 195 ns | 229 ns | 8.24 µs | 181 ns |  | 0 | **2e+02x slower** |
| `linear/3d/max10/eval10/flipped` | 333 ns | 348 ns | 17 µs |  |  |  | 51x faster |
| `linear/3d/max10/eval100/ascending` | 194 ns | 226 ns | 8.25 µs | 171 ns |  | 0 | **2e+02x slower** |
| `linear/3d/max10/eval100/flipped` | 307 ns | 339 ns | 15.7 µs |  |  |  | 51x faster |
| `linear/3d/max10/eval1000/ascending` | 196 ns | 241 ns | 7.54 µs | 181 ns |  | 0 | **2e+02x slower** |
| `linear/3d/max10/eval1000/flipped` | 333 ns | 342 ns | 15.7 µs |  |  |  | 47x faster |
| `linear/3d/max10/eval10000/ascending` | 197 ns | 239 ns | 7.53 µs | 175 ns |  | 0 | **2e+02x slower** |
| `linear/3d/max10/eval10000/flipped` | 305 ns | 336 ns | 15.6 µs |  |  |  | 51x faster |
| `linear/3d/max40/eval10/ascending` | 2.74 µs | 7.28 µs | 7.51 µs | 14.6 µs |  | 0 | **2.8e+03x slower** |
| `linear/3d/max40/eval10/flipped` | 8.95 µs | 12.7 µs | 15.6 µs |  |  |  | 1.7x faster |
| `linear/3d/max40/eval100/ascending` | 2.72 µs | 8.34 µs | 7.57 µs | 15.8 µs |  | 0 | **2.8e+03x slower** |
| `linear/3d/max40/eval100/flipped` | 8.92 µs | 12.9 µs | 17.1 µs |  |  |  | 1.9x faster |
| `linear/3d/max40/eval1000/ascending` | 15.5 µs | 5.8 µs | 7.63 µs | 6.09 µs |  | 0 | **1.8e+04x slower** |
| `linear/3d/max40/eval1000/flipped` | 26.9 µs | 12.1 µs | 15.7 µs |  |  |  | **1.7x slower** |
| `linear/3d/max40/eval10000/ascending` | 15.5 µs | 6.97 µs | 7.58 µs | 12.2 µs |  | 0 | **1.8e+04x slower** |
| `linear/3d/max40/eval10000/flipped` | 28.6 µs | 12.3 µs | 15.5 µs |  |  |  | **1.8x slower** |
| `linear/3d/max200/eval10/ascending` | 1.67 ms | 628 µs | 7.48 µs | 1.39 ms |  | 0 | **1.9e+06x slower** |
| `linear/3d/max200/eval10/flipped` | 3.27 ms | 1.84 ms | 16.1 µs |  |  |  | **2e+02x slower** |
| `linear/3d/max200/eval100/ascending` | 1.74 ms | 608 µs | 8.28 µs | 1.44 ms |  | 0 | **2e+06x slower** |
| `linear/3d/max200/eval100/flipped` | 3.06 ms | 1.91 ms | 15.9 µs |  |  |  | **1.9e+02x slower** |
| `linear/3d/max200/eval1000/ascending` | 1.74 ms | 602 µs | 7.68 µs | 1.29 ms |  | 0 | **1.8e+06x slower** |
| `linear/3d/max200/eval1000/flipped` | 3.11 ms | 1.89 ms | 15.9 µs |  |  |  | **2e+02x slower** |
| `linear/3d/max200/eval10000/ascending` | 1.78 ms | 613 µs | 7.61 µs | 1.23 ms |  | 0 | **1.8e+06x slower** |
| `linear/3d/max200/eval10000/flipped` | 3.15 ms | 1.94 ms | 16 µs |  |  |  | **2e+02x slower** |
| `linear/4d/max10/eval10/ascending` |  |  | 9.67 µs |  |  |  |  |
| `linear/4d/max10/eval10/flipped` |  |  | 17.8 µs |  |  |  |  |
| `linear/4d/max10/eval100/ascending` |  |  | 10.6 µs |  |  |  |  |
| `linear/4d/max10/eval100/flipped` |  |  | 19.5 µs |  |  |  |  |
| `linear/4d/max10/eval1000/ascending` |  |  | 9.77 µs |  |  |  |  |
| `linear/4d/max10/eval1000/flipped` |  |  | 19.6 µs |  |  |  |  |
| `linear/4d/max10/eval10000/ascending` |  |  | 9.79 µs |  |  |  |  |
| `linear/4d/max10/eval10000/flipped` |  |  | 17.8 µs |  |  |  |  |
| `linear/4d/max40/eval10/ascending` | 14.3 µs | 38.5 µs | 9.52 µs | 5.91 µs |  | 1.17 ns | **1.2e+04x slower** |
| `linear/4d/max40/eval10/flipped` | 53 µs | 70.7 µs | 21.4 µs |  |  |  | **2.5x slower** |
| `linear/4d/max40/eval100/ascending` | 13.4 µs | 31.7 µs | 9.5 µs | 5.84 µs |  | 1.27 ns | **1.1e+04x slower** |
| `linear/4d/max40/eval100/flipped` | 57.9 µs | 68.4 µs | 19.8 µs |  |  |  | **2.9x slower** |
| `linear/4d/max40/eval1000/ascending` | 14.3 µs | 29.9 µs | 10.1 µs | 6.05 µs |  | 1.27 ns | **1.1e+04x slower** |
| `linear/4d/max40/eval1000/flipped` | 57.6 µs | 71.4 µs | 21.5 µs |  |  |  | **2.7x slower** |
| `linear/4d/max40/eval10000/ascending` | 14.4 µs | 31.6 µs | 10.4 µs | 5.99 µs |  | 1.27 ns | **1.1e+04x slower** |
| `linear/4d/max40/eval10000/flipped` | 57.3 µs | 78.6 µs | 19.9 µs |  |  |  | **2.9x slower** |
| `linear/4d/max200/eval10/ascending` | 47.2 ms | 16.4 ms | 9.38 µs | 39.7 ms |  | 1.17 ns | **4e+07x slower** |
| `linear/4d/max200/eval10/flipped` | 106 ms | 69.6 ms | 22.1 µs |  |  |  | **4.8e+03x slower** |
| `linear/4d/max200/eval100/ascending` | 48.8 ms | 15.1 ms | 9.47 µs | 40.7 ms |  | 1.27 ns | **3.8e+07x slower** |
| `linear/4d/max200/eval100/flipped` | 114 ms | 70.2 ms | 20.2 µs |  |  |  | **5.6e+03x slower** |
| `linear/4d/max200/eval1000/ascending` | 49.8 ms | 16.1 ms | 10.2 µs | 41.7 ms |  | 1.27 ns | **3.9e+07x slower** |
| `linear/4d/max200/eval1000/flipped` | 100 ms | 71.8 ms | 20 µs |  |  |  | **5e+03x slower** |
| `linear/4d/max200/eval10000/ascending` | 47.5 ms | 16.5 ms | 9.31 µs | 39.8 ms |  | 1.17 ns | **4.1e+07x slower** |
| `linear/4d/max200/eval10000/flipped` | 98.9 ms | 68.1 ms | 19.9 µs |  |  |  | **5e+03x slower** |
| `nearest/2d/each40/eval10000/ascending` | 650 ns | 1.1 µs | 5.82 µs | 271 ns |  | 0 | **9.5e+02x slower** |
| `nearest/3d/each40/eval10000/ascending` | 22.9 µs | 52.6 µs | 7.51 µs | 8.53 µs |  | 0 | **2.6e+04x slower** |
| `pchip/2d/each40/eval10000/ascending` | 6.34 µs | 6.14 µs | 6.71 µs |  |  |  | 1.1x faster |
| `pchip/3d/each40/eval10000/ascending` | 465 µs | 450 µs | 8.63 µs |  |  |  | **54x slower** |
| `cubic/2d/each40/eval10000/ascending` | 8.58 µs | 8.96 µs | 1.85 ms |  |  | 0 | **9.8e+03x slower** |
| `cubic/3d/each40/eval10000/ascending` | 436 µs | 309 µs | 153 ms |  |  | 1.17 ns | **3.7e+05x slower** |
| `cubic/2d/each10/eval1000/ascending` | 1.47 µs | 1.67 µs | 1 ms |  |  | 0 | **1.9e+03x slower** |
| `cubic/2d/each40/eval1000/ascending` | 8.57 µs | 8.29 µs | 1.87 ms |  |  | 0 | **9.8e+03x slower** |
| `cubic/2d/each100/eval1000/ascending` | 76.5 µs | 82.2 µs | 7.22 ms |  |  | 0 | **8.7e+04x slower** |
| `cubic/2d/each200/eval1000/ascending` | 175 µs | 218 µs | 27.5 ms |  |  | 0 | **2e+05x slower** |
| `cubic/2d/each400/eval1000/ascending` | 1.41 ms | 1.2 ms | 143 ms |  |  | 0 | **1.6e+06x slower** |
| `quintic/2d/each10/eval1000/ascending` | 1.92 µs | 2.06 µs | 3.46 ms |  |  |  | 1.8e+03x faster |
| `quintic/2d/each40/eval1000/ascending` | 12.4 µs | 12.3 µs | 7.47 ms |  |  |  | 6.1e+02x faster |
| `linear/2d/each400/eval10000/ascending` | 302 µs | 102 µs | 5.99 µs | 271 µs | 1.14 µs | 0 | **4.4e+05x slower** |
| `linear/2d/each800/eval10000/ascending` | 1.02 ms | 414 µs | 6.33 µs | 812 µs | 2.32 µs | 0 | **1.5e+06x slower** |

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
| `nearest_nd/n100000/image500` | interpx | 12.0 | +3.5 | +2.1 | +4.1 | 18.2 |
| `nearest_nd/n100000/image500` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |
| `rbf/thin_plate_spline/global/n3000/scattered` | interpx | 4.3 | +75.9 | +3.8 | +3.8 | 80.3 |
| `rbf/thin_plate_spline/global/n3000/scattered` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |
| `rbf/thin_plate_spline/global/n3000/scattered` | scirs2-interpolate | 2.8 | +289.8 | +1.7 | +1.7 | 292.5 |
| `rbf/thin_plate_spline/k50/n100000/image` | interpx | 7.2 | +3.0 | +4.2 | +4.5 | 11.9 |
| `rbf/thin_plate_spline/k50/n100000/image` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |
| `oned/cubic_spline/n10000/sorted` | interpx | 4.9 | +0.4 | +0.4 | +0.4 | 5.9 |
| `oned/cubic_spline/n10000/sorted` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |
| `oned/cubic_spline/n10000/sorted` | scirs2-interpolate | 3.8 | +1.3 | +1.4 | +2.0 | 5.9 |
| `oned/cubic_spline/n10000/sorted` | ndarray-interp | 3.6 | +1.3 | +1.3 | +1.3 | 5.4 |
| `grid/linear/4d/max200/eval10000/ascending` | interpx | 218.3 | +186.7 | +186.4 | +186.5 | 405.0 |
| `grid/linear/4d/max200/eval10000/ascending` | interpx, `from_owned` | 222.5 | +182.5 | +182.5 | +182.5 | 405.0 |
| `grid/linear/4d/max200/eval10000/ascending` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |
| `grid/linear/4d/max200/eval10000/ascending` | scirs2-interpolate | 209.0 | +193.0 | +195.9 | +195.9 | 404.9 |
| `grid/linear/4d/max200/eval10000/ascending` | interpn | 208.2 | +0.0 | +0.0 | +0.0 | 400.5 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpx | 5.6 | +0.4 | +2.7 | +2.7 | 8.4 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpx, `from_owned` | 5.6 | +2.8 | +2.7 | +2.7 | 8.4 |
| `grid/cubic/2d/each400/eval1000/ascending` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpn | 5.4 | +0.0 | +0.0 | +0.0 | 5.5 |
| `grid/linear/2d/each800/eval10000/ascending` | interpx | 12.5 | +0.7 | +2.4 | +2.4 | 15.0 |
| `grid/linear/2d/each800/eval10000/ascending` | interpx, `from_owned` | 9.9 | +0.0 | +4.9 | +4.9 | 14.9 |
| `grid/linear/2d/each800/eval10000/ascending` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |
| `grid/linear/2d/each800/eval10000/ascending` | scirs2-interpolate | 9.3 | +3.8 | +5.5 | +5.6 | 14.9 |
| `grid/linear/2d/each800/eval10000/ascending` | ndarray-interp | 13.1 | +0.0 | +0.0 | +0.0 | 13.1 |
| `grid/linear/2d/each800/eval10000/ascending` | interpn | 11.8 | +0.0 | +0.0 | +0.0 | 11.8 |

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
| `NearestNDInterpolator` | 0e+00 |  |  |  |  |
| `RBFInterpolator` | 8e-05 |  | 3e+01 |  |  |
