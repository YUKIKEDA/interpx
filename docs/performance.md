# Performance: interpx against SciPy and other Rust crates

<!-- Written by compare/perf/report.py from a run of compare/perf/run.py. Do not edit by hand: run it again. -->

Measured on 2026-10-13, at commit `0902071871`.

## How it was measured

- Machine: Intel64 Family 6 Model 191 Stepping 2, GenuineIntel, 16 logical cores, 48 GB. Windows-11-10.0.26300-SP0.
- rustc 1.99.0 (b940084d7 2026-09-28); Python 3.13.1, NumPy 2.5.4, SciPy 1.18.1.
- The other Rust crates, at the versions in `compare/perf/runners/*/Cargo.toml`: `scirs2-interpolate` (the widest of the Rust crates: one variable, grids, RBF), `ndarray-interp` (the most used for one variable; also bilinear), `interpn` (grids).
- Every library ran each case in a process of its own, pinned to one core, on one thread, on the same data. The libraries took turns case by case, 3 rounds.
- `interpx, from_owned` is interpx built from data that is handed over (`from_owned`) and so not copied. Its column is there for the types that keep their data as it is given (the grid, `NearestNDInterpolator`, `RBFInterpolator`), next to the construction from borrowed data (`new`), which is the one the comparison counts. Each of its samples is given data made just before it, which is not in the nearest cache as the data of the other constructions is: on small grids it reads slower than `new` for that reason.
- A time is that of one call: construction, or the evaluation of all the queries of the case (`eval_batch_into` for interpx, the call with an array for the others). It is the median of the repeats of a process, and then the median of the rounds.
- interpx is called slower where its median is more than 1.10 times another library's and its fastest repeat is slower than that library's fastest.
- Memory is the resident set size of the process, sampled every millisecond from outside it.
- Command: `perf/run.py --rounds 3`.

## Where interpx is slower

| family | operation | cases | interpx fastest | slower than SciPy | slower than a Rust crate |
| :--- | :--- | ---: | ---: | ---: | ---: |
| One variable | evaluation | 112 | 111 | 0 | 1 |
| One variable | construction | 112 | 78 | 8 | 30 |
| `PPoly` | evaluation | 4 | 4 | 0 | 0 |
| `PPoly` | construction | 4 | 4 | 0 | 0 |
| `BSpline` | evaluation | 4 | 4 | 0 | 0 |
| `BSpline` | construction | 4 | 4 | 0 | 0 |
| `BSplineCurve` | evaluation | 2 | 2 | 0 | 0 |
| `BSplineCurve` | construction | 2 | 2 | 0 | 0 |
| Regular grid | evaluation | 87 | 72 | 0 | 15 |
| Regular grid | construction | 87 | 25 | 38 | 43 |
| `NearestNDInterpolator` | evaluation | 15 | 15 | 0 | 0 |
| `NearestNDInterpolator` | construction | 15 | 15 | 0 | 0 |
| `RBFInterpolator` | evaluation | 53 | 53 | 0 | 0 |
| `RBFInterpolator` | construction | 53 | 53 | 0 | 0 |

186 comparisons (a case, an operation, and a library), by method. "Same values" says whether that library's values are SciPy's up to rounding: where they are not, it computes another interpolant under the same name, and its time is not that of the same work.

| method | operation | against | same values | cases | slower by | the largest |
| :--- | :--- | :--- | :--- | ---: | ---: | :--- |
| `grid/linear` | construction | interpn | yes | 34 | 157x to 4.91e+07x | `grid/linear/4d/max200/eval10000/ascending`: 57.6 ms against 1.17 ns |
| `grid/cubic` | construction | interpn | no | 7 | 1.85e+03x to 1.87e+06x | `grid/cubic/2d/each400/eval1000/ascending`: 1.64 ms against 0 |
| `grid/nearest` | construction | interpn | yes | 2 | 900x to 2.35e+04x | `grid/nearest/3d/each40/eval10000/ascending`: 23 µs against 0 |
| `grid/linear` | construction | SciPy | yes | 36 | 1.21x to 6.05e+03x | `grid/linear/4d/max200/eval1000/flipped`: 122 ms against 20.2 µs |
| `grid/linear` | construction | ndarray-interp | yes | 14 | 1.47x to 585x | `grid/linear/2d/each800/eval10000/ascending`: 1.3 ms against 2.23 µs |
| `grid/pchip` | construction | SciPy | yes | 1 | 54.2x | `grid/pchip/3d/each40/eval10000/ascending`: 507 µs against 9.35 µs |
| `oned/linear` | construction | ndarray-interp | yes | 14 | 1.38x to 7.01x | `oned/linear/n10000/sorted`: 104 µs against 14.9 µs |
| `grid/nearest` | construction | SciPy | yes | 1 | 3.06x | `grid/nearest/3d/each40/eval10000/ascending`: 23 µs against 7.5 µs |
| `oned/pchip` | construction | scirs2-interpolate | yes | 8 | 1.37x to 2.87x | `oned/pchip/n1000/sorted`: 12.4 µs against 4.33 µs |
| `oned/bspline1` | construction | SciPy | yes | 2 | 2.7x to 2.87x | `oned/bspline1/n10000/shuffled`: 69 µs against 24.1 µs |
| `grid/linear` | construction | scirs2-interpolate | yes | 23 | 1.18x to 2.67x | `grid/linear/4d/max40/eval10/ascending`: 14.4 µs against 5.37 µs |
| `grid/nearest` | construction | scirs2-interpolate | yes | 2 | 2.1x to 2.49x | `grid/nearest/3d/each40/eval10000/ascending`: 23 µs against 9.22 µs |
| `oned/linear` | construction | SciPy | yes | 4 | 2x to 2.34x | `oned/linear/n10000/shuffled-uneven`: 104 µs against 44.3 µs |
| `oned/cubic_spline` | construction | scirs2-interpolate | no | 2 | 2.31x to 2.31x | `oned/cubic_spline/n1000/shuffled`: 43.5 µs against 18.8 µs |
| `oned/akima` | construction | scirs2-interpolate | yes | 6 | 1.54x to 1.95x | `oned/akima/n1000/sorted`: 13.3 µs against 6.81 µs |
| `oned/linear` | construction | scirs2-interpolate | yes | 12 | 1.17x to 1.67x | `oned/linear/n1000/sorted`: 2.12 µs against 1.27 µs |
| `grid/linear` | evaluation | interpn | yes | 14 | 1.1x to 1.38x | `grid/linear/3d/max10/eval10/ascending`: 198 ns against 144 ns |
| `oned/nearest` | construction | SciPy | yes | 2 | 1.25x to 1.29x | `oned/nearest/n10000/shuffled`: 69.6 µs against 53.9 µs |
| `oned/nearest` | evaluation | scirs2-interpolate | yes | 1 | 1.14x | `oned/nearest/n10/shuffled`: 331 ns against 291 ns |
| `grid/nearest` | evaluation | interpn | yes | 1 | 1.12x | `grid/nearest/2d/each40/eval10000/ascending`: 151 µs against 135 µs |

<details><summary>Every comparison where interpx is slower</summary>

"Values apart" is the distance of that library's values from SciPy's, as a share of the largest value.

| case | operation | against | interpx | the other | slower by | values apart |
| :--- | :--- | :--- | ---: | ---: | ---: | ---: |
| `grid/linear/4d/max200/eval10000/ascending` | construction | interpn | 57.6 ms | 1.17 ns | 49120512.00x | 4e-16 |
| `grid/linear/4d/max200/eval100/ascending` | construction | interpn | 57 ms | 1.27 ns | 44907204.92x | 3e-16 |
| `grid/linear/4d/max200/eval10/ascending` | construction | interpn | 51.7 ms | 1.27 ns | 40725897.85x | 2e-16 |
| `grid/linear/4d/max200/eval1000/ascending` | construction | interpn | 55 ms | 1.37 ns | 40258340.57x | 3e-16 |
| `grid/linear/3d/max200/eval10000/ascending` | construction | interpn | 2.32 ms | 0 | 2375372.80x | 5e-16 |
| `grid/linear/3d/max200/eval100/ascending` | construction | interpn | 2.09 ms | 0 | 2141593.60x | 3e-16 |
| `grid/linear/3d/max200/eval1000/ascending` | construction | interpn | 2.09 ms | 0 | 2140518.40x | 3e-16 |
| `grid/linear/3d/max200/eval10/ascending` | construction | interpn | 2.08 ms | 0 | 2127052.80x | 3e-16 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | interpn | 1.3 ms | 0 | 1904201.14x | 3e-16 |
| `grid/cubic/2d/each400/eval1000/ascending` | construction | interpn | 1.64 ms | 0 | 1867434.67x | 5e-01 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | interpn | 373 µs | 0 | 544987.43x | 3e-16 |
| `grid/cubic/3d/each40/eval10000/ascending` | construction | interpn | 431 µs | 0 | 441241.60x | 5e-01 |
| `grid/cubic/2d/each200/eval1000/ascending` | construction | interpn | 187 µs | 0 | 212878.22x | 2e-01 |
| `grid/cubic/2d/each100/eval1000/ascending` | construction | interpn | 84 µs | 0 | 95573.33x | 2e-01 |
| `grid/linear/2d/max200/eval10000/ascending` | construction | interpn | 51.6 µs | 0 | 66032.00x | 3e-16 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | interpn | 23 µs | 0 | 23532.80x | 0e+00 |
| `grid/linear/3d/max40/eval10000/ascending` | construction | interpn | 18.9 µs | 0 | 19328.00x | 3e-16 |
| `grid/linear/3d/max40/eval1000/ascending` | construction | interpn | 18 µs | 0 | 18470.40x | 2e-16 |
| `grid/cubic/2d/each40/eval1000/ascending` | construction | interpn | 9.22 µs | 0 | 11798.00x | 7e-01 |
| `grid/linear/4d/max40/eval10/ascending` | construction | interpn | 14.4 µs | 1.27 ns | 11308.31x | 2e-16 |
| `grid/linear/4d/max40/eval1000/ascending` | construction | interpn | 14.4 µs | 1.37 ns | 10548.57x | 6e-16 |
| `grid/cubic/2d/each40/eval10000/ascending` | construction | interpn | 9.26 µs | 0 | 10531.56x | 7e-01 |
| `grid/linear/4d/max40/eval10000/ascending` | construction | interpn | 13.3 µs | 1.27 ns | 10500.92x | 5e-16 |
| `grid/linear/4d/max40/eval100/ascending` | construction | interpn | 13.9 µs | 1.37 ns | 10166.86x | 2e-16 |
| `grid/linear/2d/max200/eval10/ascending` | construction | interpn | 7.55 µs | 0 | 9660.00x | 2e-16 |
| `grid/linear/2d/max200/eval100/ascending` | construction | interpn | 7.52 µs | 0 | 9626.00x | 1e-16 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | interpn | 7.48 µs | 0 | 9580.00x | 2e-16 |
| `grid/linear/4d/max200/eval1000/flipped` | construction | SciPy | 122 ms | 20.2 µs | 6053.82x |  |
| `grid/linear/4d/max200/eval10000/ascending` | construction | SciPy | 57.6 ms | 9.59 µs | 6004.95x |  |
| `grid/linear/4d/max200/eval10/flipped` | construction | SciPy | 122 ms | 20.3 µs | 5997.96x |  |
| `grid/linear/4d/max200/eval100/flipped` | construction | SciPy | 120 ms | 20.2 µs | 5957.37x |  |
| `grid/linear/4d/max200/eval1000/ascending` | construction | SciPy | 55 ms | 9.46 µs | 5820.56x |  |
| `grid/linear/4d/max200/eval10000/flipped` | construction | SciPy | 123 ms | 22 µs | 5593.58x |  |
| `grid/linear/4d/max200/eval100/ascending` | construction | SciPy | 57 ms | 10.3 µs | 5526.67x |  |
| `grid/linear/4d/max200/eval10/ascending` | construction | SciPy | 51.7 ms | 9.42 µs | 5491.17x |  |
| `grid/linear/3d/max40/eval10/ascending` | construction | interpn | 2.97 µs | 0 | 3042.80x | 2e-16 |
| `grid/linear/3d/max40/eval100/ascending` | construction | interpn | 2.97 µs | 0 | 3041.60x | 2e-16 |
| `grid/cubic/2d/each10/eval1000/ascending` | construction | interpn | 1.62 µs | 0 | 1847.11x | 3e-01 |
| `grid/nearest/2d/each40/eval10000/ascending` | construction | interpn | 615 ns | 0 | 900.00x | 0e+00 |
| `grid/linear/2d/max40/eval10000/ascending` | construction | interpn | 401 ns | 0 | 586.43x | 4e-16 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | ndarray-interp | 1.3 ms | 2.23 µs | 584.83x | 6e-16 |
| `grid/linear/2d/max40/eval1000/ascending` | construction | interpn | 396 ns | 0 | 579.43x | 3e-16 |
| `grid/linear/2d/max40/eval100/ascending` | construction | interpn | 387 ns | 0 | 565.86x | 2e-16 |
| `grid/linear/2d/max40/eval10/ascending` | construction | interpn | 372 ns | 0 | 544.86x | 2e-16 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | ndarray-interp | 373 µs | 1.25 µs | 297.95x | 5e-16 |
| `grid/linear/3d/max200/eval10000/ascending` | construction | SciPy | 2.32 ms | 8.37 µs | 277.19x |  |
| `grid/linear/3d/max200/eval100/ascending` | construction | SciPy | 2.09 ms | 7.61 µs | 274.73x |  |
| `grid/linear/3d/max200/eval1000/ascending` | construction | SciPy | 2.09 ms | 7.68 µs | 272.25x |  |
| `grid/linear/3d/max200/eval10/ascending` | construction | SciPy | 2.08 ms | 8.34 µs | 248.95x |  |
| `grid/linear/3d/max200/eval10/flipped` | construction | SciPy | 3.61 ms | 16.2 µs | 222.58x |  |
| `grid/linear/3d/max200/eval1000/flipped` | construction | SciPy | 3.55 ms | 16.2 µs | 219.46x |  |
| `grid/linear/3d/max200/eval10000/flipped` | construction | SciPy | 3.64 ms | 17.5 µs | 208.08x |  |
| `grid/linear/2d/each800/eval10000/ascending` | construction | SciPy | 1.3 ms | 6.33 µs | 205.50x |  |
| `grid/linear/3d/max200/eval100/flipped` | construction | SciPy | 3.47 ms | 17.3 µs | 200.59x |  |
| `grid/linear/3d/max10/eval100/ascending` | construction | interpn | 194 ns | 0 | 199.00x | 3e-16 |
| `grid/linear/3d/max10/eval10000/ascending` | construction | interpn | 192 ns | 0 | 196.75x | 4e-16 |
| `grid/linear/3d/max10/eval10/ascending` | construction | interpn | 178 ns | 0 | 182.40x | 2e-16 |
| `grid/linear/3d/max10/eval1000/ascending` | construction | interpn | 178 ns | 0 | 182.30x | 5e-16 |
| `grid/linear/2d/max10/eval100/ascending` | construction | interpn | 133 ns | 0 | 170.75x | 2e-16 |
| `grid/linear/2d/max10/eval10000/ascending` | construction | interpn | 133 ns | 0 | 170.38x | 3e-16 |
| `grid/linear/2d/max10/eval1000/ascending` | construction | interpn | 125 ns | 0 | 160.00x | 3e-16 |
| `grid/linear/2d/max10/eval10/ascending` | construction | interpn | 123 ns | 0 | 157.38x | 1e-16 |
| `grid/linear/2d/max200/eval10000/ascending` | construction | ndarray-interp | 51.6 µs | 516 ns | 99.93x | 5e-16 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | SciPy | 373 µs | 6.01 µs | 61.99x |  |
| `grid/pchip/3d/each40/eval10000/ascending` | construction | SciPy | 507 µs | 9.35 µs | 54.21x |  |
| `grid/linear/2d/max200/eval1000/ascending` | construction | ndarray-interp | 7.48 µs | 489 ns | 15.31x | 4e-16 |
| `grid/linear/2d/max200/eval10/ascending` | construction | ndarray-interp | 7.55 µs | 516 ns | 14.63x | 2e-16 |
| `grid/linear/2d/max200/eval100/ascending` | construction | ndarray-interp | 7.52 µs | 515 ns | 14.59x | 3e-16 |
| `grid/linear/2d/max200/eval10000/ascending` | construction | SciPy | 51.6 µs | 5.86 µs | 8.81x |  |
| `oned/linear/n10000/sorted` | construction | ndarray-interp | 104 µs | 14.9 µs | 7.01x | 0e+00 |
| `grid/linear/2d/max200/eval10000/flipped` | construction | SciPy | 61.9 µs | 13 µs | 4.74x |  |
| `oned/linear/n10000/sorted-uneven` | construction | ndarray-interp | 104 µs | 32.7 µs | 3.17x | 0e+00 |
| `oned/linear/n10000/shuffled-uneven` | construction | ndarray-interp | 104 µs | 32.8 µs | 3.17x | 0e+00 |
| `oned/linear/n10000/shuffled` | construction | ndarray-interp | 102 µs | 33 µs | 3.09x | 0e+00 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | SciPy | 23 µs | 7.5 µs | 3.06x |  |
| `oned/pchip/n1000/sorted` | construction | scirs2-interpolate | 12.4 µs | 4.33 µs | 2.87x | 2e-16 |
| `oned/bspline1/n10000/shuffled` | construction | SciPy | 69 µs | 24.1 µs | 2.87x |  |
| `oned/pchip/n1000/shuffled` | construction | scirs2-interpolate | 12.9 µs | 4.58 µs | 2.82x | 2e-16 |
| `oned/bspline1/n10000/sorted` | construction | SciPy | 68.3 µs | 25.3 µs | 2.70x |  |
| `grid/linear/4d/max40/eval10/ascending` | construction | scirs2-interpolate | 14.4 µs | 5.37 µs | 2.67x | 2e-16 |
| `grid/linear/4d/max40/eval10/flipped` | construction | SciPy | 57.2 µs | 21.5 µs | 2.66x |  |
| `grid/linear/4d/max40/eval1000/flipped` | construction | SciPy | 52.2 µs | 19.8 µs | 2.64x |  |
| `grid/linear/2d/max40/eval1000/ascending` | construction | ndarray-interp | 396 ns | 152 ns | 2.61x | 4e-16 |
| `grid/linear/4d/max40/eval100/flipped` | construction | SciPy | 54.8 µs | 21.2 µs | 2.59x |  |
| `grid/linear/2d/max40/eval100/ascending` | construction | ndarray-interp | 387 ns | 152 ns | 2.55x | 2e-16 |
| `grid/linear/4d/max40/eval1000/ascending` | construction | scirs2-interpolate | 14.4 µs | 5.67 µs | 2.55x | 2e-16 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | scirs2-interpolate | 7.48 µs | 2.98 µs | 2.52x | 2e-16 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | scirs2-interpolate | 23 µs | 9.22 µs | 2.49x | 0e+00 |
| `grid/linear/4d/max40/eval100/ascending` | construction | scirs2-interpolate | 13.9 µs | 5.59 µs | 2.49x | 3e-16 |
| `grid/linear/2d/max40/eval10/ascending` | construction | ndarray-interp | 372 ns | 152 ns | 2.45x | 2e-16 |
| `grid/linear/2d/max40/eval10000/ascending` | construction | ndarray-interp | 401 ns | 164 ns | 2.44x | 4e-16 |
| `grid/linear/2d/max200/eval100/ascending` | construction | scirs2-interpolate | 7.52 µs | 3.1 µs | 2.42x | 2e-16 |
| `grid/linear/4d/max40/eval10000/ascending` | construction | scirs2-interpolate | 13.3 µs | 5.52 µs | 2.42x | 4e-16 |
| `grid/linear/4d/max40/eval10000/flipped` | construction | SciPy | 51.9 µs | 21.5 µs | 2.41x |  |
| `oned/linear/n10000/shuffled-uneven` | construction | SciPy | 104 µs | 44.3 µs | 2.34x |  |
| `grid/linear/2d/max200/eval10/ascending` | construction | scirs2-interpolate | 7.55 µs | 3.25 µs | 2.32x | 1e-16 |
| `oned/cubic_spline/n1000/shuffled` | construction | scirs2-interpolate | 43.5 µs | 18.8 µs | 2.31x | 4e-04 |
| `oned/cubic_spline/n1000/sorted` | construction | scirs2-interpolate | 43.2 µs | 18.7 µs | 2.31x | 4e-04 |
| `oned/linear/n10000/sorted-uneven` | construction | SciPy | 104 µs | 49.3 µs | 2.11x |  |
| `grid/nearest/2d/each40/eval10000/ascending` | construction | scirs2-interpolate | 615 ns | 293 ns | 2.10x | 0e+00 |
| `oned/linear/n10000/shuffled` | construction | SciPy | 102 µs | 50.6 µs | 2.02x |  |
| `oned/linear/n10/shuffled` | construction | ndarray-interp | 95 ns | 47.4 ns | 2.01x | 0e+00 |
| `oned/linear/n10000/sorted` | construction | SciPy | 104 µs | 52 µs | 2.00x |  |
| `oned/linear/n10/sorted` | construction | ndarray-interp | 94.7 ns | 47.6 ns | 1.99x | 0e+00 |
| `oned/akima/n1000/sorted` | construction | scirs2-interpolate | 13.3 µs | 6.81 µs | 1.95x | 3e-17 |
| `grid/linear/3d/max40/eval10000/flipped` | construction | SciPy | 29.4 µs | 15.7 µs | 1.87x |  |
| `grid/linear/2d/max40/eval1000/ascending` | construction | scirs2-interpolate | 396 ns | 215 ns | 1.85x | 3e-16 |
| `grid/linear/3d/max40/eval1000/flipped` | construction | SciPy | 31.3 µs | 17.1 µs | 1.83x |  |
| `oned/pchip/n100/shuffled-uneven` | construction | scirs2-interpolate | 977 ns | 545 ns | 1.79x | 2e-16 |
| `oned/pchip/n100/shuffled` | construction | scirs2-interpolate | 978 ns | 551 ns | 1.77x | 2e-16 |
| `oned/pchip/n100/sorted-uneven` | construction | scirs2-interpolate | 977 ns | 555 ns | 1.76x | 2e-16 |
| `grid/linear/2d/max40/eval10000/ascending` | construction | scirs2-interpolate | 401 ns | 229 ns | 1.75x | 2e-16 |
| `grid/linear/2d/max10/eval100/ascending` | construction | ndarray-interp | 133 ns | 77.8 ns | 1.71x | 4e-16 |
| `grid/linear/2d/max40/eval100/ascending` | construction | scirs2-interpolate | 387 ns | 227 ns | 1.70x | 3e-16 |
| `grid/linear/2d/max10/eval10000/ascending` | construction | ndarray-interp | 133 ns | 78.3 ns | 1.70x | 4e-16 |
| `oned/akima/n100/sorted` | construction | scirs2-interpolate | 989 ns | 582 ns | 1.70x | 1e-16 |
| `oned/akima/n100/shuffled` | construction | scirs2-interpolate | 982 ns | 578 ns | 1.70x | 1e-16 |
| `oned/linear/n1000/sorted` | construction | scirs2-interpolate | 2.12 µs | 1.27 µs | 1.67x | 0e+00 |
| `oned/linear/n1000/shuffled` | construction | scirs2-interpolate | 2.09 µs | 1.28 µs | 1.63x | 0e+00 |
| `grid/linear/2d/max40/eval10/ascending` | construction | scirs2-interpolate | 372 ns | 229 ns | 1.63x | 4e-17 |
| `oned/linear/n10000/shuffled-uneven` | construction | scirs2-interpolate | 104 µs | 64.1 µs | 1.62x | 6e-17 |
| `oned/linear/n10000/sorted-uneven` | construction | scirs2-interpolate | 104 µs | 64.1 µs | 1.62x | 6e-17 |
| `oned/akima/n50/sorted` | construction | scirs2-interpolate | 549 ns | 341 ns | 1.61x | 1e-16 |
| `oned/linear/n50/sorted` | construction | ndarray-interp | 163 ns | 102 ns | 1.60x | 0e+00 |
| `oned/akima/n1000/shuffled` | construction | scirs2-interpolate | 12.9 µs | 8.06 µs | 1.60x | 3e-17 |
| `oned/linear/n10000/shuffled` | construction | scirs2-interpolate | 102 µs | 64.3 µs | 1.58x | 0e+00 |
| `oned/linear/n100/sorted-uneven` | construction | ndarray-interp | 266 ns | 169 ns | 1.58x | 0e+00 |
| `grid/linear/2d/max10/eval10/ascending` | construction | ndarray-interp | 123 ns | 78 ns | 1.58x | 1e-16 |
| `oned/linear/n100/shuffled-uneven` | construction | ndarray-interp | 265 ns | 169 ns | 1.57x | 0e+00 |
| `grid/linear/3d/max200/eval10000/ascending` | construction | scirs2-interpolate | 2.32 ms | 1.48 ms | 1.56x | 3e-16 |
| `oned/linear/n10000/sorted` | construction | scirs2-interpolate | 104 µs | 67.1 µs | 1.55x | 0e+00 |
| `oned/pchip/n100/sorted` | construction | scirs2-interpolate | 900 ns | 584 ns | 1.54x | 2e-16 |
| `oned/akima/n50/shuffled` | construction | scirs2-interpolate | 549 ns | 357 ns | 1.54x | 1e-16 |
| `grid/linear/4d/max40/eval1000/ascending` | construction | SciPy | 14.4 µs | 9.5 µs | 1.52x |  |
| `grid/linear/3d/max200/eval1000/ascending` | construction | scirs2-interpolate | 2.09 ms | 1.38 ms | 1.52x | 3e-16 |
| `grid/linear/4d/max40/eval10/ascending` | construction | SciPy | 14.4 µs | 9.48 µs | 1.52x |  |
| `grid/linear/2d/each800/eval10000/ascending` | construction | scirs2-interpolate | 1.3 ms | 860 µs | 1.51x | 2e-16 |
| `grid/linear/3d/max200/eval10/ascending` | construction | scirs2-interpolate | 2.08 ms | 1.41 ms | 1.47x | 3e-16 |
| `grid/linear/2d/max10/eval1000/ascending` | construction | ndarray-interp | 125 ns | 85.3 ns | 1.47x | 3e-16 |
| `oned/linear/n50/shuffled` | construction | ndarray-interp | 163 ns | 112 ns | 1.46x | 0e+00 |
| `grid/linear/3d/max200/eval100/ascending` | construction | scirs2-interpolate | 2.09 ms | 1.43 ms | 1.46x | 3e-16 |
| `grid/linear/4d/max40/eval100/ascending` | construction | SciPy | 13.9 µs | 9.55 µs | 1.46x |  |
| `oned/linear/n100/sorted` | construction | ndarray-interp | 244 ns | 168 ns | 1.45x | 0e+00 |
| `oned/linear/n1000/sorted` | construction | ndarray-interp | 2.12 µs | 1.47 µs | 1.44x | 0e+00 |
| `oned/pchip/n50/shuffled` | construction | scirs2-interpolate | 529 ns | 368 ns | 1.44x | 1e-16 |
| `oned/linear/n100/shuffled` | construction | ndarray-interp | 264 ns | 184 ns | 1.43x | 0e+00 |
| `grid/linear/4d/max40/eval10000/ascending` | construction | SciPy | 13.3 µs | 9.4 µs | 1.42x |  |
| `grid/linear/3d/max40/eval10000/ascending` | construction | scirs2-interpolate | 18.9 µs | 13.3 µs | 1.41x | 2e-16 |
| `oned/linear/n100/shuffled-uneven` | construction | scirs2-interpolate | 265 ns | 190 ns | 1.39x | 3e-17 |
| `oned/linear/n1000/shuffled` | construction | ndarray-interp | 2.09 µs | 1.52 µs | 1.38x | 0e+00 |
| `oned/linear/n100/shuffled` | construction | scirs2-interpolate | 264 ns | 191 ns | 1.38x | 0e+00 |
| `oned/linear/n100/sorted-uneven` | construction | scirs2-interpolate | 266 ns | 193 ns | 1.38x | 3e-17 |
| `grid/linear/3d/max10/eval10/ascending` | evaluation | interpn | 198 ns | 144 ns | 1.38x | 2e-16 |
| `oned/pchip/n50/sorted` | construction | scirs2-interpolate | 534 ns | 389 ns | 1.37x | 1e-16 |
| `grid/linear/3d/max10/eval1000/ascending` | evaluation | interpn | 20.7 µs | 15.3 µs | 1.36x | 5e-16 |
| `grid/linear/3d/max10/eval10000/ascending` | evaluation | interpn | 193 µs | 143 µs | 1.35x | 4e-16 |
| `grid/linear/2d/max10/eval1000/ascending` | evaluation | interpn | 13.6 µs | 10.1 µs | 1.35x | 3e-16 |
| `grid/linear/2d/max200/eval1000/flipped` | construction | SciPy | 15.7 µs | 12 µs | 1.31x |  |
| `grid/linear/4d/max200/eval100/ascending` | construction | scirs2-interpolate | 57 ms | 43.6 ms | 1.31x | 1e-16 |
| `grid/linear/2d/max200/eval10/flipped` | construction | SciPy | 15.9 µs | 12.2 µs | 1.31x |  |
| `grid/linear/2d/max10/eval10/ascending` | evaluation | interpn | 135 ns | 104 ns | 1.30x | 1e-16 |
| `oned/nearest/n10000/shuffled` | construction | SciPy | 69.6 µs | 53.9 µs | 1.29x |  |
| `grid/linear/4d/max200/eval10000/ascending` | construction | scirs2-interpolate | 57.6 ms | 44.9 ms | 1.28x | 4e-16 |
| `grid/linear/3d/max40/eval1000/ascending` | construction | scirs2-interpolate | 18 µs | 14.2 µs | 1.27x | 2e-16 |
| `grid/linear/2d/max200/eval10/ascending` | construction | SciPy | 7.55 µs | 5.97 µs | 1.27x |  |
| `grid/linear/2d/max200/eval100/ascending` | construction | SciPy | 7.52 µs | 5.95 µs | 1.26x |  |
| `grid/linear/3d/max10/eval100/ascending` | evaluation | interpn | 1.91 µs | 1.52 µs | 1.25x | 3e-16 |
| `oned/nearest/n10000/sorted` | construction | SciPy | 69.4 µs | 55.5 µs | 1.25x |  |
| `grid/linear/2d/max40/eval10000/ascending` | evaluation | interpn | 166 µs | 133 µs | 1.25x | 4e-16 |
| `grid/linear/2d/max10/eval100/ascending` | evaluation | interpn | 1.26 µs | 1.02 µs | 1.24x | 2e-16 |
| `grid/linear/2d/max10/eval10000/ascending` | evaluation | interpn | 125 µs | 101 µs | 1.23x | 3e-16 |
| `oned/linear/n50/shuffled` | construction | scirs2-interpolate | 163 ns | 133 ns | 1.23x | 0e+00 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | SciPy | 7.48 µs | 6.13 µs | 1.22x |  |
| `grid/linear/2d/max200/eval100/flipped` | construction | SciPy | 15.9 µs | 13.1 µs | 1.21x |  |
| `grid/linear/2d/max40/eval100/ascending` | evaluation | interpn | 1.62 µs | 1.33 µs | 1.21x | 2e-16 |
| `grid/linear/4d/max200/eval1000/ascending` | construction | scirs2-interpolate | 55 ms | 45.4 ms | 1.21x | 3e-16 |
| `grid/linear/4d/max200/eval10/ascending` | construction | scirs2-interpolate | 51.7 ms | 43.1 ms | 1.20x | 2e-16 |
| `oned/linear/n50/sorted` | construction | scirs2-interpolate | 163 ns | 137 ns | 1.19x | 0e+00 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | scirs2-interpolate | 373 µs | 315 µs | 1.18x | 2e-16 |
| `oned/linear/n100/sorted` | construction | scirs2-interpolate | 244 ns | 208 ns | 1.17x | 0e+00 |
| `grid/linear/2d/max40/eval1000/ascending` | evaluation | interpn | 15.4 µs | 13.2 µs | 1.17x | 3e-16 |
| `oned/nearest/n10/shuffled` | evaluation | scirs2-interpolate | 331 ns | 291 ns | 1.14x | 0e+00 |
| `grid/linear/4d/max200/eval1000/ascending` | evaluation | interpn | 40.9 µs | 36.5 µs | 1.12x | 3e-16 |
| `grid/linear/4d/max40/eval10/ascending` | evaluation | interpn | 308 ns | 276 ns | 1.12x | 2e-16 |
| `grid/nearest/2d/each40/eval10000/ascending` | evaluation | interpn | 151 µs | 135 µs | 1.12x | 0e+00 |
| `grid/linear/3d/max40/eval10/ascending` | evaluation | interpn | 222 ns | 202 ns | 1.10x | 2e-16 |

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
| `linear/n10/sorted` | 143 ns | 5.39 µs | 266 ns | 277 ns | 1.9x faster |
| `linear/n10/shuffled` | 173 ns | 5.54 µs | 286 ns | 250 ns | 1.4x faster |
| `nearest/n10/sorted` | 166 ns | 10.4 µs | 278 ns |  | 1.7x faster |
| `nearest/n10/shuffled` | 331 ns | 9.59 µs | 291 ns |  | **1.1x slower** |
| `previous/n10/sorted` | 159 ns | 11.8 µs |  |  | 74x faster |
| `previous/n10/shuffled` | 285 ns | 11.7 µs |  |  | 41x faster |
| `bspline1/n10/sorted` | 209 ns | 2.71 µs |  |  | 13x faster |
| `bspline1/n10/shuffled` | 578 ns | 2.95 µs |  |  | 5.1x faster |
| `bspline2/n10/sorted` | 271 ns | 3.33 µs |  |  | 12x faster |
| `bspline2/n10/shuffled` | 340 ns | 3.36 µs |  |  | 9.9x faster |
| `bspline3/n10/sorted` | 375 ns | 3.58 µs |  |  | 9.5x faster |
| `bspline3/n10/shuffled` | 439 ns | 3.63 µs |  |  | 8.3x faster |
| `cubic_spline/n10/sorted` | 140 ns | 3.93 µs | 296 ns | 354 ns | 2.1x faster |
| `cubic_spline/n10/shuffled` | 223 ns | 3.72 µs | 297 ns | 338 ns | 1.3x faster |
| `pchip/n10/sorted` | 140 ns | 3.63 µs | 316 ns |  | 2.3x faster |
| `pchip/n10/shuffled` | 199 ns | 3.92 µs | 317 ns |  | 1.6x faster |
| `akima/n10/sorted` | 139 ns | 3.7 µs | 273 ns |  | 2x faster |
| `akima/n10/shuffled` | 224 ns | 3.7 µs | 285 ns |  | 1.3x faster |
| `makima/n10/sorted` | 140 ns | 3.64 µs |  |  | 26x faster |
| `makima/n10/shuffled` | 207 ns | 3.68 µs |  |  | 18x faster |
| `linear/n50/sorted` | 726 ns | 6.18 µs | 1.67 µs | 1.35 µs | 1.9x faster |
| `linear/n50/shuffled` | 890 ns | 8.21 µs | 1.66 µs | 1.38 µs | 1.6x faster |
| `nearest/n50/sorted` | 791 ns | 12.3 µs | 1.7 µs |  | 2.2x faster |
| `nearest/n50/shuffled` | 1.57 µs | 11.5 µs | 1.71 µs |  | 1.1x faster |
| `previous/n50/sorted` | 776 ns | 13.9 µs |  |  | 18x faster |
| `previous/n50/shuffled` | 1.35 µs | 12.7 µs |  |  | 9.5x faster |
| `bspline1/n50/sorted` | 1.18 µs | 4.02 µs |  |  | 3.4x faster |
| `bspline1/n50/shuffled` | 2.78 µs | 4.82 µs |  |  | 1.7x faster |
| `bspline2/n50/sorted` | 1.37 µs | 6.82 µs |  |  | 5x faster |
| `bspline2/n50/shuffled` | 2.03 µs | 6.83 µs |  |  | 3.4x faster |
| `bspline3/n50/sorted` | 1.86 µs | 9.49 µs |  |  | 5.1x faster |
| `bspline3/n50/shuffled` | 2.82 µs | 9.3 µs |  |  | 3.3x faster |
| `cubic_spline/n50/sorted` | 734 ns | 5.28 µs | 4.31 µs | 1.58 µs | 2.2x faster |
| `cubic_spline/n50/shuffled` | 1.08 µs | 5.52 µs | 4.14 µs | 1.74 µs | 1.6x faster |
| `pchip/n50/sorted` | 742 ns | 5.53 µs | 4.37 µs |  | 5.9x faster |
| `pchip/n50/shuffled` | 985 ns | 5.5 µs | 3.92 µs |  | 4x faster |
| `akima/n50/sorted` | 686 ns | 5.23 µs | 4.53 µs |  | 6.6x faster |
| `akima/n50/shuffled` | 983 ns | 5.51 µs | 4.02 µs |  | 4.1x faster |
| `makima/n50/sorted` | 689 ns | 5.29 µs |  |  | 7.7x faster |
| `makima/n50/shuffled` | 971 ns | 5.67 µs |  |  | 5.8x faster |
| `linear/n100/sorted` | 1.41 µs | 7.32 µs | 3.6 µs | 2.47 µs | 1.7x faster |
| `linear/n100/shuffled` | 1.78 µs | 13.1 µs | 3.33 µs | 2.72 µs | 1.5x faster |
| `nearest/n100/sorted` | 1.52 µs | 13.6 µs | 3.36 µs |  | 2.2x faster |
| `nearest/n100/shuffled` | 3.12 µs | 13.8 µs | 3.63 µs |  | 1.2x faster |
| `previous/n100/sorted` | 1.57 µs | 16.5 µs |  |  | 11x faster |
| `previous/n100/shuffled` | 2.7 µs | 15.4 µs |  |  | 5.7x faster |
| `bspline1/n100/sorted` | 2.33 µs | 5.38 µs |  |  | 2.3x faster |
| `bspline1/n100/shuffled` | 5.72 µs | 12 µs |  |  | 2.1x faster |
| `bspline2/n100/sorted` | 2.7 µs | 10.5 µs |  |  | 3.9x faster |
| `bspline2/n100/shuffled` | 4.86 µs | 12.7 µs |  |  | 2.6x faster |
| `bspline3/n100/sorted` | 3.62 µs | 14.9 µs |  |  | 4.1x faster |
| `bspline3/n100/shuffled` | 6.15 µs | 19.5 µs |  |  | 3.2x faster |
| `cubic_spline/n100/sorted` | 1.41 µs | 6.91 µs | 18.9 µs | 3.44 µs | 2.4x faster |
| `cubic_spline/n100/shuffled` | 2.04 µs | 7.87 µs | 19.3 µs | 3.21 µs | 1.6x faster |
| `pchip/n100/sorted` | 1.42 µs | 6.85 µs | 16 µs |  | 4.8x faster |
| `pchip/n100/shuffled` | 1.95 µs | 7.91 µs | 16.4 µs |  | 4x faster |
| `akima/n100/sorted` | 1.47 µs | 6.98 µs | 16.6 µs |  | 4.7x faster |
| `akima/n100/shuffled` | 2.08 µs | 7.91 µs | 17.2 µs |  | 3.8x faster |
| `makima/n100/sorted` | 1.43 µs | 7 µs |  |  | 4.9x faster |
| `makima/n100/shuffled` | 1.96 µs | 7.98 µs |  |  | 4.1x faster |
| `linear/n1000/sorted` | 14 µs | 29.2 µs | 94.2 µs | 26.8 µs | 1.9x faster |
| `linear/n1000/shuffled` | 16.1 µs | 143 µs | 214 µs | 24.7 µs | 1.5x faster |
| `nearest/n1000/sorted` | 14.9 µs | 61.5 µs | 124 µs |  | 4.1x faster |
| `nearest/n1000/shuffled` | 33.4 µs | 215 µs | 236 µs |  | 6.4x faster |
| `previous/n1000/sorted` | 14.5 µs | 64.8 µs |  |  | 4.5x faster |
| `previous/n1000/shuffled` | 27.1 µs | 218 µs |  |  | 8x faster |
| `bspline1/n1000/sorted` | 26.5 µs | 31 µs |  |  | 1.2x faster |
| `bspline1/n1000/shuffled` | 61.7 µs | 538 µs |  |  | 8.7x faster |
| `bspline2/n1000/sorted` | 32.2 µs | 78.3 µs |  |  | 2.4x faster |
| `bspline2/n1000/shuffled` | 63.7 µs | 593 µs |  |  | 9.3x faster |
| `bspline3/n1000/sorted` | 36.1 µs | 129 µs |  |  | 3.6x faster |
| `bspline3/n1000/shuffled` | 86.5 µs | 631 µs |  |  | 7.3x faster |
| `cubic_spline/n1000/sorted` | 15.1 µs | 40.2 µs | 1.55 ms | 34.4 µs | 2.3x faster |
| `cubic_spline/n1000/shuffled` | 21.6 µs | 237 µs | 1.47 ms | 32 µs | 1.5x faster |
| `pchip/n1000/sorted` | 14.6 µs | 39.7 µs | 1.21 ms |  | 2.7x faster |
| `pchip/n1000/shuffled` | 21.5 µs | 221 µs | 1.31 ms |  | 10x faster |
| `akima/n1000/sorted` | 14.1 µs | 39 µs | 1.54 ms |  | 2.8x faster |
| `akima/n1000/shuffled` | 19.7 µs | 222 µs | 1.55 ms |  | 11x faster |
| `makima/n1000/sorted` | 14.7 µs | 40.3 µs |  |  | 2.7x faster |
| `makima/n1000/shuffled` | 19.7 µs | 225 µs |  |  | 11x faster |
| `linear/n10000/sorted` | 159 µs | 235 µs | 1.34 ms | 270 µs | 1.5x faster |
| `linear/n10000/shuffled` | 208 µs | 2.38 ms | 3.04 ms | 270 µs | 1.3x faster |
| `nearest/n10000/sorted` | 161 µs | 589 µs | 1.54 ms |  | 3.7x faster |
| `nearest/n10000/shuffled` | 385 µs | 3.17 ms | 3.26 ms |  | 8.2x faster |
| `previous/n10000/sorted` | 154 µs | 659 µs |  |  | 4.3x faster |
| `previous/n10000/shuffled` | 349 µs | 3.03 ms |  |  | 8.7x faster |
| `bspline1/n10000/sorted` | 265 µs | 302 µs |  |  | 1.1x faster |
| `bspline1/n10000/shuffled` | 650 µs | 48.7 ms |  |  | 75x faster |
| `bspline2/n10000/sorted` | 337 µs | 818 µs |  |  | 2.4x faster |
| `bspline2/n10000/shuffled` | 971 µs | 49.4 ms |  |  | 51x faster |
| `bspline3/n10000/sorted` | 401 µs | 1.37 ms |  |  | 3.4x faster |
| `bspline3/n10000/shuffled` | 1.23 ms | 53.5 ms |  |  | 44x faster |
| `cubic_spline/n10000/sorted` | 149 µs | 382 µs | 150 ms | 349 µs | 2.3x faster |
| `cubic_spline/n10000/shuffled` | 237 µs | 3.67 ms | 147 ms | 364 µs | 1.5x faster |
| `pchip/n10000/sorted` | 147 µs | 372 µs | 121 ms |  | 2.5x faster |
| `pchip/n10000/shuffled` | 238 µs | 3.41 ms | 122 ms |  | 14x faster |
| `akima/n10000/sorted` | 154 µs | 368 µs | 144 ms |  | 2.4x faster |
| `akima/n10000/shuffled` | 217 µs | 3.42 ms | 146 ms |  | 16x faster |
| `makima/n10000/sorted` | 147 µs | 367 µs |  |  | 2.5x faster |
| `makima/n10000/shuffled` | 217 µs | 3.65 ms |  |  | 17x faster |
| `linear/n100/sorted-uneven` | 1.4 µs | 7.49 µs | 3.39 µs | 3.76 µs | 2.4x faster |
| `linear/n100/shuffled-uneven` | 3.03 µs | 13.8 µs | 3.5 µs | 4.06 µs | 1.2x faster |
| `cubic_spline/n100/sorted-uneven` | 1.38 µs | 8.05 µs | 17.2 µs | 4.41 µs | 3.2x faster |
| `cubic_spline/n100/shuffled-uneven` | 2.95 µs | 8.04 µs | 19.4 µs | 4.65 µs | 1.6x faster |
| `pchip/n100/sorted-uneven` | 1.41 µs | 8.03 µs | 15.9 µs |  | 5.7x faster |
| `pchip/n100/shuffled-uneven` | 3.18 µs | 7.99 µs | 16.4 µs |  | 2.5x faster |
| `linear/n10000/sorted-uneven` | 220 µs | 262 µs | 1.33 ms | 973 µs | 1.2x faster |
| `linear/n10000/shuffled-uneven` | 674 µs | 2.42 ms | 3.03 ms | 1.75 ms | 2.6x faster |
| `cubic_spline/n10000/sorted-uneven` | 236 µs | 476 µs | 149 ms | 967 µs | 2x faster |
| `cubic_spline/n10000/shuffled-uneven` | 796 µs | 3.65 ms | 148 ms | 1.82 ms | 2.3x faster |
| `pchip/n10000/sorted-uneven` | 241 µs | 449 µs | 121 ms |  | 1.9x faster |
| `pchip/n10000/shuffled-uneven` | 804 µs | 3.36 ms | 122 ms |  | 4.2x faster |

### Construction

| case | interpx | SciPy | scirs2-interpolate | ndarray-interp | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: |
| `linear/n10/sorted` | 94.7 ns | 9.45 µs | 89.6 ns | 47.6 ns | **2x slower** |
| `linear/n10/shuffled` | 95 ns | 9.49 µs | 88.5 ns | 47.4 ns | **2x slower** |
| `nearest/n10/sorted` | 70.8 ns | 10.4 µs | 90.2 ns |  | 1.3x faster |
| `nearest/n10/shuffled` | 71.1 ns | 10.5 µs | 89.3 ns |  | 1.3x faster |
| `previous/n10/sorted` | 71.6 ns | 10 µs |  |  | 1.4e+02x faster |
| `previous/n10/shuffled` | 70.6 ns | 9.48 µs |  |  | 1.3e+02x faster |
| `bspline1/n10/sorted` | 113 ns | 16.8 µs |  |  | 1.5e+02x faster |
| `bspline1/n10/shuffled` | 114 ns | 16.8 µs |  |  | 1.5e+02x faster |
| `bspline2/n10/sorted` | 483 ns | 38.3 µs |  |  | 79x faster |
| `bspline2/n10/shuffled` | 517 ns | 39.2 µs |  |  | 76x faster |
| `bspline3/n10/sorted` | 579 ns | 38.4 µs |  |  | 66x faster |
| `bspline3/n10/shuffled` | 524 ns | 38.8 µs |  |  | 74x faster |
| `cubic_spline/n10/sorted` | 257 ns | 93.9 µs | 412 ns | 515 ns | 1.6x faster |
| `cubic_spline/n10/shuffled` | 258 ns | 91.4 µs | 412 ns | 516 ns | 1.6x faster |
| `pchip/n10/sorted` | 194 ns | 108 µs | 193 ns |  | **1x slower** |
| `pchip/n10/shuffled` | 195 ns | 111 µs | 211 ns |  | 1.1x faster |
| `akima/n10/sorted` | 199 ns | 76.1 µs | 190 ns |  | **1x slower** |
| `akima/n10/shuffled` | 199 ns | 74.8 µs | 197 ns |  | **1x slower** |
| `makima/n10/sorted` | 209 ns | 77.8 µs |  |  | 3.7e+02x faster |
| `makima/n10/shuffled` | 209 ns | 83.9 µs |  |  | 4e+02x faster |
| `linear/n50/sorted` | 163 ns | 9.72 µs | 137 ns | 102 ns | **1.6x slower** |
| `linear/n50/shuffled` | 163 ns | 9.56 µs | 133 ns | 112 ns | **1.5x slower** |
| `nearest/n50/sorted` | 104 ns | 10.7 µs | 133 ns |  | 1.3x faster |
| `nearest/n50/shuffled` | 103 ns | 9.92 µs | 133 ns |  | 1.3x faster |
| `previous/n50/sorted` | 104 ns | 10.4 µs |  |  | 1e+02x faster |
| `previous/n50/shuffled` | 95.3 ns | 9.55 µs |  |  | 1e+02x faster |
| `bspline1/n50/sorted` | 158 ns | 16.5 µs |  |  | 1e+02x faster |
| `bspline1/n50/shuffled` | 144 ns | 16 µs |  |  | 1.1e+02x faster |
| `bspline2/n50/sorted` | 2.04 µs | 42.1 µs |  |  | 21x faster |
| `bspline2/n50/shuffled` | 2.03 µs | 40.1 µs |  |  | 20x faster |
| `bspline3/n50/sorted` | 2.24 µs | 41 µs |  |  | 18x faster |
| `bspline3/n50/shuffled` | 2.25 µs | 38.8 µs |  |  | 17x faster |
| `cubic_spline/n50/sorted` | 837 ns | 87.4 µs | 935 ns | 1.09 µs | 1.1x faster |
| `cubic_spline/n50/shuffled` | 916 ns | 89.6 µs | 943 ns | 1.18 µs | 1x faster |
| `pchip/n50/sorted` | 534 ns | 112 µs | 389 ns |  | **1.4x slower** |
| `pchip/n50/shuffled` | 529 ns | 111 µs | 368 ns |  | **1.4x slower** |
| `akima/n50/sorted` | 549 ns | 81.7 µs | 341 ns |  | **1.6x slower** |
| `akima/n50/shuffled` | 549 ns | 75.5 µs | 357 ns |  | **1.5x slower** |
| `makima/n50/sorted` | 543 ns | 84 µs |  |  | 1.5e+02x faster |
| `makima/n50/shuffled` | 595 ns | 79.3 µs |  |  | 1.3e+02x faster |
| `linear/n100/sorted` | 244 ns | 9.11 µs | 208 ns | 168 ns | **1.5x slower** |
| `linear/n100/shuffled` | 264 ns | 9.27 µs | 191 ns | 184 ns | **1.4x slower** |
| `nearest/n100/sorted` | 157 ns | 11 µs | 194 ns |  | 1.2x faster |
| `nearest/n100/shuffled` | 168 ns | 10 µs | 191 ns |  | 1.1x faster |
| `previous/n100/sorted` | 168 ns | 10.7 µs |  |  | 64x faster |
| `previous/n100/shuffled` | 154 ns | 9.91 µs |  |  | 65x faster |
| `bspline1/n100/sorted` | 214 ns | 17 µs |  |  | 79x faster |
| `bspline1/n100/shuffled` | 206 ns | 17 µs |  |  | 83x faster |
| `bspline2/n100/sorted` | 3.85 µs | 46.7 µs |  |  | 12x faster |
| `bspline2/n100/shuffled` | 3.83 µs | 42.9 µs |  |  | 11x faster |
| `bspline3/n100/sorted` | 4.3 µs | 48.6 µs |  |  | 11x faster |
| `bspline3/n100/shuffled` | 4.31 µs | 44.5 µs |  |  | 10x faster |
| `cubic_spline/n100/sorted` | 1.73 µs | 90.1 µs | 1.67 µs | 2.06 µs | **1x slower** |
| `cubic_spline/n100/shuffled` | 1.74 µs | 95.2 µs | 1.67 µs | 1.92 µs | **1x slower** |
| `pchip/n100/sorted` | 900 ns | 110 µs | 584 ns |  | **1.5x slower** |
| `pchip/n100/shuffled` | 978 ns | 112 µs | 551 ns |  | **1.8x slower** |
| `akima/n100/sorted` | 989 ns | 79.2 µs | 582 ns |  | **1.7x slower** |
| `akima/n100/shuffled` | 982 ns | 78.9 µs | 578 ns |  | **1.7x slower** |
| `makima/n100/sorted` | 980 ns | 80.4 µs |  |  | 82x faster |
| `makima/n100/shuffled` | 1.06 µs | 86.9 µs |  |  | 82x faster |
| `linear/n1000/sorted` | 2.12 µs | 13.8 µs | 1.27 µs | 1.47 µs | **1.7x slower** |
| `linear/n1000/shuffled` | 2.09 µs | 12.9 µs | 1.28 µs | 1.52 µs | **1.6x slower** |
| `nearest/n1000/sorted` | 1.05 µs | 15.6 µs | 1.29 µs |  | 1.2x faster |
| `nearest/n1000/shuffled` | 1.09 µs | 15.8 µs | 1.28 µs |  | 1.2x faster |
| `previous/n1000/sorted` | 1.08 µs | 17.9 µs |  |  | 17x faster |
| `previous/n1000/shuffled` | 1.01 µs | 17.9 µs |  |  | 18x faster |
| `bspline1/n1000/sorted` | 1.08 µs | 17.4 µs |  |  | 16x faster |
| `bspline1/n1000/shuffled` | 1.14 µs | 17.1 µs |  |  | 15x faster |
| `bspline2/n1000/sorted` | 34.7 µs | 118 µs |  |  | 3.4x faster |
| `bspline2/n1000/shuffled` | 37.8 µs | 119 µs |  |  | 3.1x faster |
| `bspline3/n1000/sorted` | 72.8 µs | 127 µs |  |  | 1.7x faster |
| `bspline3/n1000/shuffled` | 75.3 µs | 137 µs |  |  | 1.8x faster |
| `cubic_spline/n1000/sorted` | 43.2 µs | 122 µs | 18.7 µs | 16.3 µs | **2.7x slower** |
| `cubic_spline/n1000/shuffled` | 43.5 µs | 122 µs | 18.8 µs | 17.4 µs | **2.5x slower** |
| `pchip/n1000/sorted` | 12.4 µs | 129 µs | 4.33 µs |  | **2.9x slower** |
| `pchip/n1000/shuffled` | 12.9 µs | 125 µs | 4.58 µs |  | **2.8x slower** |
| `akima/n1000/sorted` | 13.3 µs | 102 µs | 6.81 µs |  | **2x slower** |
| `akima/n1000/shuffled` | 12.9 µs | 103 µs | 8.06 µs |  | **1.6x slower** |
| `makima/n1000/sorted` | 12.9 µs | 106 µs |  |  | 8.2x faster |
| `makima/n1000/shuffled` | 13.3 µs | 109 µs |  |  | 8.2x faster |
| `linear/n10000/sorted` | 104 µs | 52 µs | 67.1 µs | 14.9 µs | **7x slower** |
| `linear/n10000/shuffled` | 102 µs | 50.6 µs | 64.3 µs | 33 µs | **3.1x slower** |
| `nearest/n10000/sorted` | 69.4 µs | 55.5 µs | 65.4 µs |  | **1.2x slower** |
| `nearest/n10000/shuffled` | 69.6 µs | 53.9 µs | 63.4 µs |  | **1.3x slower** |
| `previous/n10000/sorted` | 68.3 µs | 77.2 µs |  |  | 1.1x faster |
| `previous/n10000/shuffled` | 65.7 µs | 81.3 µs |  |  | 1.2x faster |
| `bspline1/n10000/sorted` | 68.3 µs | 25.3 µs |  |  | **2.7x slower** |
| `bspline1/n10000/shuffled` | 69 µs | 24.1 µs |  |  | **2.9x slower** |
| `bspline2/n10000/sorted` | 539 µs | 805 µs |  |  | 1.5x faster |
| `bspline2/n10000/shuffled` | 538 µs | 767 µs |  |  | 1.4x faster |
| `bspline3/n10000/sorted` | 658 µs | 1.02 ms |  |  | 1.5x faster |
| `bspline3/n10000/shuffled` | 648 µs | 1.01 ms |  |  | 1.6x faster |
| `cubic_spline/n10000/sorted` | 256 µs | 777 µs | 349 µs | 393 µs | 1.4x faster |
| `cubic_spline/n10000/shuffled` | 256 µs | 782 µs | 298 µs | 392 µs | 1.2x faster |
| `pchip/n10000/sorted` | 94.2 µs | 260 µs | 174 µs |  | 1.8x faster |
| `pchip/n10000/shuffled` | 139 µs | 258 µs | 174 µs |  | 1.2x faster |
| `akima/n10000/sorted` | 96.3 µs | 684 µs | 158 µs |  | 1.6x faster |
| `akima/n10000/shuffled` | 96.4 µs | 695 µs | 141 µs |  | 1.5x faster |
| `makima/n10000/sorted` | 145 µs | 741 µs |  |  | 5.1x faster |
| `makima/n10000/shuffled` | 101 µs | 781 µs |  |  | 7.7x faster |
| `linear/n100/sorted-uneven` | 266 ns | 9.12 µs | 193 ns | 169 ns | **1.6x slower** |
| `linear/n100/shuffled-uneven` | 265 ns | 9.23 µs | 190 ns | 169 ns | **1.6x slower** |
| `cubic_spline/n100/sorted-uneven` | 1.73 µs | 96 µs | 1.58 µs | 1.88 µs | **1.1x slower** |
| `cubic_spline/n100/shuffled-uneven` | 1.7 µs | 96.1 µs | 1.67 µs | 2.05 µs | **1x slower** |
| `pchip/n100/sorted-uneven` | 977 ns | 112 µs | 555 ns |  | **1.8x slower** |
| `pchip/n100/shuffled-uneven` | 977 ns | 106 µs | 545 ns |  | **1.8x slower** |
| `linear/n10000/sorted-uneven` | 104 µs | 49.3 µs | 64.1 µs | 32.7 µs | **3.2x slower** |
| `linear/n10000/shuffled-uneven` | 104 µs | 44.3 µs | 64.1 µs | 32.8 µs | **3.2x slower** |
| `cubic_spline/n10000/sorted-uneven` | 255 µs | 784 µs | 294 µs | 398 µs | 1.2x faster |
| `cubic_spline/n10000/shuffled-uneven` | 256 µs | 795 µs | 296 µs | 386 µs | 1.2x faster |
| `pchip/n10000/sorted-uneven` | 103 µs | 247 µs | 176 µs |  | 1.7x faster |
| `pchip/n10000/shuffled-uneven` | 97.1 µs | 255 µs | 176 µs |  | 1.8x faster |

## `PPoly`

A cubic of 55 intervals, evaluated at points in no order.

### Evaluation

| case | interpx | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: |
| `intervals55/eval10` | 72.4 ns | 3.44 µs | 48x faster |
| `intervals55/eval100` | 664 ns | 4.37 µs | 6.6x faster |
| `intervals55/eval1000` | 6.53 µs | 16 µs | 2.5x faster |
| `intervals55/eval10000` | 66.1 µs | 333 µs | 5x faster |

### Construction

| case | interpx | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: |
| `intervals55/eval10` | 176 ns | 5.61 µs | 32x faster |
| `intervals55/eval100` | 176 ns | 5.2 µs | 30x faster |
| `intervals55/eval1000` | 175 ns | 5.64 µs | 32x faster |
| `intervals55/eval10000` | 175 ns | 5.19 µs | 30x faster |

## `BSpline`

A cubic B-spline of 55 intervals, evaluated at points in no order.

### Evaluation

| case | interpx | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: |
| `intervals55/eval10` | 159 ns | 2.73 µs | 17x faster |
| `intervals55/eval100` | 1.55 µs | 6.25 µs | 4x faster |
| `intervals55/eval1000` | 14.2 µs | 37.7 µs | 2.7x faster |
| `intervals55/eval10000` | 155 µs | 371 µs | 2.4x faster |

### Construction

| case | interpx | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: |
| `intervals55/eval10` | 175 ns | 9.39 µs | 54x faster |
| `intervals55/eval100` | 175 ns | 9.38 µs | 53x faster |
| `intervals55/eval1000` | 162 ns | 10 µs | 62x faster |
| `intervals55/eval10000` | 162 ns | 9.21 µs | 57x faster |

## `BSplineCurve`

A cubic curve through `n` points of the plane, evaluated at `4 n` parameters in order.

### Evaluation

| case | interpx | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: |
| `degree3/n100` | 4.13 µs | 16.4 µs | 4x faster |
| `degree3/n10000` | 445 µs | 1.33 ms | 3x faster |

### Construction

| case | interpx | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: |
| `degree3/n100` | 4.46 µs | 51.1 µs | 11x faster |
| `degree3/n10000` | 707 µs | 1.43 ms | 2x faster |

## Regular grid

Points in no order inside the grid. `maxN`: the first axis has `N` points and each axis half the one before (SciPy's suite); `eachN`: every axis has `N` points; `flipped`: every axis given from its high end.

### Evaluation

| case | interpx | interpx, `from_owned` | SciPy | scirs2-interpolate | ndarray-interp | interpn | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/2d/max10/eval10/ascending` | 135 ns | 147 ns | 19.5 µs | 804 ns | 152 ns | 104 ns | **1.3x slower** |
| `linear/2d/max10/eval10/flipped` | 135 ns | 137 ns | 21.2 µs |  |  |  | 1.6e+02x faster |
| `linear/2d/max10/eval100/ascending` | 1.26 µs | 1.26 µs | 21.9 µs | 8.48 µs | 1.49 µs | 1.02 µs | **1.2x slower** |
| `linear/2d/max10/eval100/flipped` | 1.25 µs | 1.26 µs | 21.9 µs |  |  |  | 17x faster |
| `linear/2d/max10/eval1000/ascending` | 13.6 µs | 12.6 µs | 40.7 µs | 82.6 µs | 14.2 µs | 10.1 µs | **1.3x slower** |
| `linear/2d/max10/eval1000/flipped` | 12.4 µs | 13.6 µs | 41.5 µs |  |  |  | 3.3x faster |
| `linear/2d/max10/eval10000/ascending` | 125 µs | 125 µs | 390 µs | 882 µs | 228 µs | 101 µs | **1.2x slower** |
| `linear/2d/max10/eval10000/flipped` | 129 µs | 125 µs | 358 µs |  |  |  | 2.8x faster |
| `linear/2d/max40/eval10/ascending` | 158 ns | 172 ns | 19.8 µs | 1.06 µs | 169 ns | 147 ns | **1.1x slower** |
| `linear/2d/max40/eval10/flipped` | 159 ns | 162 ns | 20.7 µs |  |  |  | 1.3e+02x faster |
| `linear/2d/max40/eval100/ascending` | 1.62 µs | 1.66 µs | 23.8 µs | 9.94 µs | 1.5 µs | 1.33 µs | **1.2x slower** |
| `linear/2d/max40/eval100/flipped` | 1.52 µs | 1.55 µs | 22.4 µs |  |  |  | 15x faster |
| `linear/2d/max40/eval1000/ascending` | 15.4 µs | 15.1 µs | 46.5 µs | 109 µs | 15.7 µs | 13.2 µs | **1.2x slower** |
| `linear/2d/max40/eval1000/flipped` | 15.2 µs | 16.5 µs | 47.1 µs |  |  |  | 3.1x faster |
| `linear/2d/max40/eval10000/ascending` | 166 µs | 156 µs | 612 µs | 1.14 ms | 224 µs | 133 µs | **1.2x slower** |
| `linear/2d/max40/eval10000/flipped` | 165 µs | 153 µs | 612 µs |  |  |  | 3.7x faster |
| `linear/2d/max200/eval10/ascending` | 193 ns | 211 ns | 19.5 µs | 2.02 µs | 205 ns | 180 ns | **1.1x slower** |
| `linear/2d/max200/eval10/flipped` | 193 ns | 210 ns | 19.5 µs |  |  |  | 1e+02x faster |
| `linear/2d/max200/eval100/ascending` | 1.96 µs | 2.06 µs | 23.2 µs | 20.8 µs | 2.18 µs | 1.89 µs | **1x slower** |
| `linear/2d/max200/eval100/flipped` | 2.06 µs | 1.91 µs | 23 µs |  |  |  | 11x faster |
| `linear/2d/max200/eval1000/ascending` | 18.8 µs | 20.6 µs | 74.6 µs | 206 µs | 24.6 µs | 18.8 µs | 1x faster |
| `linear/2d/max200/eval1000/flipped` | 20.6 µs | 20.6 µs | 69.2 µs |  |  |  | 3.4x faster |
| `linear/2d/max200/eval10000/ascending` | 206 µs | 206 µs | 819 µs | 2.3 ms | 485 µs | 188 µs | **1.1x slower** |
| `linear/2d/max200/eval10000/flipped` | 207 µs | 206 µs | 813 µs |  |  |  | 3.9x faster |
| `linear/3d/max10/eval10/ascending` | 198 ns | 199 ns | 57.7 µs | 3.14 µs |  | 144 ns | **1.4x slower** |
| `linear/3d/max10/eval10/flipped` | 199 ns | 201 ns | 57.1 µs |  |  |  | 2.9e+02x faster |
| `linear/3d/max10/eval100/ascending` | 1.91 µs | 1.91 µs | 67.9 µs | 32 µs |  | 1.52 µs | **1.3x slower** |
| `linear/3d/max10/eval100/flipped` | 1.93 µs | 2.09 µs | 64.9 µs |  |  |  | 34x faster |
| `linear/3d/max10/eval1000/ascending` | 20.7 µs | 20.7 µs | 118 µs | 332 µs |  | 15.3 µs | **1.4x slower** |
| `linear/3d/max10/eval1000/flipped` | 19.1 µs | 20.7 µs | 121 µs |  |  |  | 6.3x faster |
| `linear/3d/max10/eval10000/ascending` | 193 µs | 208 µs | 1.17 ms | 3.35 ms |  | 143 µs | **1.3x slower** |
| `linear/3d/max10/eval10000/flipped` | 193 µs | 209 µs | 1.15 ms |  |  |  | 6x faster |
| `linear/3d/max40/eval10/ascending` | 222 ns | 221 ns | 57.9 µs | 3.38 µs |  | 202 ns | **1.1x slower** |
| `linear/3d/max40/eval10/flipped` | 220 ns | 220 ns | 61.5 µs |  |  |  | 2.8e+02x faster |
| `linear/3d/max40/eval100/ascending` | 2.16 µs | 2.15 µs | 63.7 µs | 36.2 µs |  | 2.17 µs | 1x faster |
| `linear/3d/max40/eval100/flipped` | 2.13 µs | 2.19 µs | 64.8 µs |  |  |  | 30x faster |
| `linear/3d/max40/eval1000/ascending` | 21.4 µs | 22.7 µs | 138 µs | 338 µs |  | 20 µs | **1.1x slower** |
| `linear/3d/max40/eval1000/flipped` | 23.4 µs | 23.4 µs | 141 µs |  |  |  | 6x faster |
| `linear/3d/max40/eval10000/ascending` | 216 µs | 233 µs | 1.51 ms | 3.53 ms |  | 205 µs | **1.1x slower** |
| `linear/3d/max40/eval10000/flipped` | 212 µs | 213 µs | 1.52 ms |  |  |  | 7.2x faster |
| `linear/3d/max200/eval10/ascending` | 300 ns | 301 ns | 61.4 µs | 4.88 µs |  | 287 ns | **1x slower** |
| `linear/3d/max200/eval10/flipped` | 301 ns | 289 ns | 57.3 µs |  |  |  | 1.9e+02x faster |
| `linear/3d/max200/eval100/ascending` | 2.97 µs | 2.97 µs | 66.3 µs | 50.6 µs |  | 2.86 µs | **1x slower** |
| `linear/3d/max200/eval100/flipped` | 2.97 µs | 2.71 µs | 71.1 µs |  |  |  | 24x faster |
| `linear/3d/max200/eval1000/ascending` | 28.2 µs | 27.1 µs | 196 µs | 501 µs |  | 28.9 µs | 1x faster |
| `linear/3d/max200/eval1000/flipped` | 29.7 µs | 29.6 µs | 191 µs |  |  |  | 6.5x faster |
| `linear/3d/max200/eval10000/ascending` | 297 µs | 297 µs | 2.01 ms | 5.12 ms |  | 289 µs | **1x slower** |
| `linear/3d/max200/eval10000/flipped` | 297 µs | 296 µs | 2.08 ms |  |  |  | 7x faster |
| `linear/4d/max10/eval10/ascending` | 283 ns | 284 ns | 96.8 µs |  |  |  | 3.4e+02x faster |
| `linear/4d/max10/eval10/flipped` | 273 ns | 275 ns | 99.6 µs |  |  |  | 3.6e+02x faster |
| `linear/4d/max10/eval100/ascending` | 2.64 µs | 2.64 µs | 106 µs |  |  |  | 40x faster |
| `linear/4d/max10/eval100/flipped` | 2.52 µs | 2.46 µs | 115 µs |  |  |  | 46x faster |
| `linear/4d/max10/eval1000/ascending` | 26.7 µs | 26.7 µs | 227 µs |  |  |  | 8.5x faster |
| `linear/4d/max10/eval1000/flipped` | 24.4 µs | 25.3 µs | 210 µs |  |  |  | 8.6x faster |
| `linear/4d/max10/eval10000/ascending` | 244 µs | 268 µs | 1.95 ms |  |  |  | 8x faster |
| `linear/4d/max10/eval10000/flipped` | 268 µs | 252 µs | 1.95 ms |  |  |  | 7.3x faster |
| `linear/4d/max40/eval10/ascending` | 308 ns | 332 ns | 95.8 µs | 5.95 µs |  | 276 ns | **1.1x slower** |
| `linear/4d/max40/eval10/flipped` | 310 ns | 333 ns | 103 µs |  |  |  | 3.3e+02x faster |
| `linear/4d/max40/eval100/ascending` | 3.16 µs | 3.22 µs | 109 µs | 57.1 µs |  | 2.95 µs | **1.1x slower** |
| `linear/4d/max40/eval100/flipped` | 2.95 µs | 3.19 µs | 117 µs |  |  |  | 40x faster |
| `linear/4d/max40/eval1000/ascending` | 29.7 µs | 32.3 µs | 248 µs | 575 µs |  | 27 µs | **1.1x slower** |
| `linear/4d/max40/eval1000/flipped` | 29.8 µs | 41.7 µs | 253 µs |  |  |  | 8.5x faster |
| `linear/4d/max40/eval10000/ascending` | 295 µs | 319 µs | 2.25 ms | 6.1 ms |  | 272 µs | **1.1x slower** |
| `linear/4d/max40/eval10000/flipped` | 294 µs | 318 µs | 2.32 ms |  |  |  | 7.9x faster |
| `linear/4d/max200/eval10/ascending` | 384 ns | 380 ns | 96.4 µs | 7.79 µs |  | 383 ns | **1x slower** |
| `linear/4d/max200/eval10/flipped` | 385 ns | 379 ns | 96.2 µs |  |  |  | 2.5e+02x faster |
| `linear/4d/max200/eval100/ascending` | 3.75 µs | 3.74 µs | 119 µs | 72.2 µs |  | 3.79 µs | 1x faster |
| `linear/4d/max200/eval100/flipped` | 3.71 µs | 3.73 µs | 109 µs |  |  |  | 29x faster |
| `linear/4d/max200/eval1000/ascending` | 40.9 µs | 40.7 µs | 307 µs | 712 µs |  | 36.5 µs | **1.1x slower** |
| `linear/4d/max200/eval1000/flipped` | 37.5 µs | 37.5 µs | 297 µs |  |  |  | 7.9x faster |
| `linear/4d/max200/eval10000/ascending` | 375 µs | 410 µs | 3.13 ms | 7.35 ms |  | 364 µs | **1x slower** |
| `linear/4d/max200/eval10000/flipped` | 376 µs | 375 µs | 2.95 ms |  |  |  | 7.9x faster |
| `nearest/2d/each40/eval10000/ascending` | 151 µs | 165 µs | 646 µs | 1.62 ms |  | 135 µs | **1.1x slower** |
| `nearest/3d/each40/eval10000/ascending` | 232 µs | 233 µs | 1.01 ms | 1.94 ms |  | 229 µs | **1x slower** |
| `pchip/2d/each40/eval10000/ascending` | 500 µs | 501 µs | 1.21 s |  |  |  | 2.4e+03x faster |
| `pchip/3d/each40/eval10000/ascending` | 1.52 ms | 1.49 ms | 3.15 s |  |  |  | 2.1e+03x faster |
| `cubic/2d/each40/eval10000/ascending` | 302 µs | 302 µs | 1.36 ms |  |  | 417 µs | 1.4x faster |
| `cubic/3d/each40/eval10000/ascending` | 722 µs | 662 µs | 3.21 ms |  |  | 822 µs | 1.1x faster |
| `cubic/2d/each10/eval1000/ascending` | 26.6 µs | 26.6 µs | 142 µs |  |  | 35.8 µs | 1.3x faster |
| `cubic/2d/each40/eval1000/ascending` | 30.1 µs | 30.2 µs | 145 µs |  |  | 41.8 µs | 1.4x faster |
| `cubic/2d/each100/eval1000/ascending` | 32.1 µs | 32.2 µs | 183 µs |  |  | 43.2 µs | 1.3x faster |
| `cubic/2d/each200/eval1000/ascending` | 34.7 µs | 34.9 µs | 200 µs |  |  | 46 µs | 1.3x faster |
| `cubic/2d/each400/eval1000/ascending` | 41.6 µs | 38.9 µs | 274 µs |  |  | 49.6 µs | 1.2x faster |
| `quintic/2d/each10/eval1000/ascending` | 48.5 µs | 46.3 µs | 248 µs |  |  |  | 5.1x faster |
| `quintic/2d/each40/eval1000/ascending` | 48.1 µs | 48 µs | 240 µs |  |  |  | 5x faster |
| `linear/2d/each400/eval10000/ascending` | 241 µs | 220 µs | 965 µs | 4.47 ms | 542 µs | 220 µs | **1.1x slower** |
| `linear/2d/each800/eval10000/ascending` | 262 µs | 263 µs | 1.1 ms | 7.58 ms | 588 µs | 242 µs | **1.1x slower** |

### Construction

| case | interpx | interpx, `from_owned` | SciPy | scirs2-interpolate | ndarray-interp | interpn | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/2d/max10/eval10/ascending` | 123 ns | 134 ns | 5.78 µs | 145 ns | 78 ns | 0 | **1.6e+02x slower** |
| `linear/2d/max10/eval10/flipped` | 158 ns | 165 ns | 12.3 µs |  |  |  | 78x faster |
| `linear/2d/max10/eval100/ascending` | 133 ns | 128 ns | 5.88 µs | 146 ns | 77.8 ns | 0 | **1.7e+02x slower** |
| `linear/2d/max10/eval100/flipped` | 165 ns | 167 ns | 12.4 µs |  |  |  | 75x faster |
| `linear/2d/max10/eval1000/ascending` | 125 ns | 134 ns | 6.33 µs | 145 ns | 85.3 ns | 0 | **1.6e+02x slower** |
| `linear/2d/max10/eval1000/flipped` | 165 ns | 168 ns | 11.6 µs |  |  |  | 70x faster |
| `linear/2d/max10/eval10000/ascending` | 133 ns | 133 ns | 5.83 µs | 145 ns | 78.3 ns | 0 | **1.7e+02x slower** |
| `linear/2d/max10/eval10000/flipped` | 165 ns | 162 ns | 12.5 µs |  |  |  | 75x faster |
| `linear/2d/max40/eval10/ascending` | 372 ns | 687 ns | 6.36 µs | 229 ns | 152 ns | 0 | **5.4e+02x slower** |
| `linear/2d/max40/eval10/flipped` | 697 ns | 1.14 µs | 11.9 µs |  |  |  | 17x faster |
| `linear/2d/max40/eval100/ascending` | 387 ns | 757 ns | 6.33 µs | 227 ns | 152 ns | 0 | **5.7e+02x slower** |
| `linear/2d/max40/eval100/flipped` | 698 ns | 1.03 µs | 11.8 µs |  |  |  | 17x faster |
| `linear/2d/max40/eval1000/ascending` | 396 ns | 667 ns | 6.44 µs | 215 ns | 152 ns | 0 | **5.8e+02x slower** |
| `linear/2d/max40/eval1000/flipped` | 642 ns | 1.12 µs | 11.9 µs |  |  |  | 18x faster |
| `linear/2d/max40/eval10000/ascending` | 401 ns | 724 ns | 6.39 µs | 229 ns | 164 ns | 0 | **5.9e+02x slower** |
| `linear/2d/max40/eval10000/flipped` | 696 ns | 1.05 µs | 12.3 µs |  |  |  | 18x faster |
| `linear/2d/max200/eval10/ascending` | 7.55 µs | 21.2 µs | 5.97 µs | 3.25 µs | 516 ns | 0 | **9.7e+03x slower** |
| `linear/2d/max200/eval10/flipped` | 15.9 µs | 27.1 µs | 12.2 µs |  |  |  | **1.3x slower** |
| `linear/2d/max200/eval100/ascending` | 7.52 µs | 19.9 µs | 5.95 µs | 3.1 µs | 515 ns | 0 | **9.6e+03x slower** |
| `linear/2d/max200/eval100/flipped` | 15.9 µs | 28 µs | 13.1 µs |  |  |  | **1.2x slower** |
| `linear/2d/max200/eval1000/ascending` | 7.48 µs | 20.8 µs | 6.13 µs | 2.98 µs | 489 ns | 0 | **9.6e+03x slower** |
| `linear/2d/max200/eval1000/flipped` | 15.7 µs | 27.5 µs | 12 µs |  |  |  | **1.3x slower** |
| `linear/2d/max200/eval10000/ascending` | 51.6 µs | 21.5 µs | 5.86 µs | 62.1 µs | 516 ns | 0 | **6.6e+04x slower** |
| `linear/2d/max200/eval10000/flipped` | 61.9 µs | 30.1 µs | 13 µs |  |  |  | **4.7x slower** |
| `linear/3d/max10/eval10/ascending` | 178 ns | 227 ns | 8.27 µs | 180 ns |  | 0 | **1.8e+02x slower** |
| `linear/3d/max10/eval10/flipped` | 327 ns | 346 ns | 15.6 µs |  |  |  | 48x faster |
| `linear/3d/max10/eval100/ascending` | 194 ns | 207 ns | 8.33 µs | 180 ns |  | 0 | **2e+02x slower** |
| `linear/3d/max10/eval100/flipped` | 299 ns | 338 ns | 16.8 µs |  |  |  | 56x faster |
| `linear/3d/max10/eval1000/ascending` | 178 ns | 207 ns | 7.66 µs | 169 ns |  | 0 | **1.8e+02x slower** |
| `linear/3d/max10/eval1000/flipped` | 325 ns | 345 ns | 17.1 µs |  |  |  | 52x faster |
| `linear/3d/max10/eval10000/ascending` | 192 ns | 236 ns | 7.59 µs | 180 ns |  | 0 | **2e+02x slower** |
| `linear/3d/max10/eval10000/flipped` | 328 ns | 331 ns | 15.8 µs |  |  |  | 48x faster |
| `linear/3d/max40/eval10/ascending` | 2.97 µs | 8.17 µs | 7.82 µs | 4.71 µs |  | 0 | **3e+03x slower** |
| `linear/3d/max40/eval10/flipped` | 8.86 µs | 12.2 µs | 16.7 µs |  |  |  | 1.9x faster |
| `linear/3d/max40/eval100/ascending` | 2.97 µs | 6.65 µs | 7.7 µs | 7.56 µs |  | 0 | **3e+03x slower** |
| `linear/3d/max40/eval100/flipped` | 8.91 µs | 11.3 µs | 15.8 µs |  |  |  | 1.8x faster |
| `linear/3d/max40/eval1000/ascending` | 18 µs | 8.6 µs | 7.75 µs | 14.2 µs |  | 0 | **1.8e+04x slower** |
| `linear/3d/max40/eval1000/flipped` | 31.3 µs | 14.1 µs | 17.1 µs |  |  |  | **1.8x slower** |
| `linear/3d/max40/eval10000/ascending` | 18.9 µs | 7.02 µs | 7.83 µs | 13.3 µs |  | 0 | **1.9e+04x slower** |
| `linear/3d/max40/eval10000/flipped` | 29.4 µs | 11.5 µs | 15.7 µs |  |  |  | **1.9x slower** |
| `linear/3d/max200/eval10/ascending` | 2.08 ms | 769 µs | 8.34 µs | 1.41 ms |  | 0 | **2.1e+06x slower** |
| `linear/3d/max200/eval10/flipped` | 3.61 ms | 2.08 ms | 16.2 µs |  |  |  | **2.2e+02x slower** |
| `linear/3d/max200/eval100/ascending` | 2.09 ms | 770 µs | 7.61 µs | 1.43 ms |  | 0 | **2.1e+06x slower** |
| `linear/3d/max200/eval100/flipped` | 3.47 ms | 2.18 ms | 17.3 µs |  |  |  | **2e+02x slower** |
| `linear/3d/max200/eval1000/ascending` | 2.09 ms | 766 µs | 7.68 µs | 1.38 ms |  | 0 | **2.1e+06x slower** |
| `linear/3d/max200/eval1000/flipped` | 3.55 ms | 2.08 ms | 16.2 µs |  |  |  | **2.2e+02x slower** |
| `linear/3d/max200/eval10000/ascending` | 2.32 ms | 842 µs | 8.37 µs | 1.48 ms |  | 0 | **2.4e+06x slower** |
| `linear/3d/max200/eval10000/flipped` | 3.64 ms | 2.19 ms | 17.5 µs |  |  |  | **2.1e+02x slower** |
| `linear/4d/max10/eval10/ascending` | 230 ns | 274 ns | 10.5 µs |  |  |  | 46x faster |
| `linear/4d/max10/eval10/flipped` | 365 ns | 367 ns | 19.5 µs |  |  |  | 53x faster |
| `linear/4d/max10/eval100/ascending` | 231 ns | 229 ns | 10.7 µs |  |  |  | 46x faster |
| `linear/4d/max10/eval100/flipped` | 364 ns | 372 ns | 18.8 µs |  |  |  | 52x faster |
| `linear/4d/max10/eval1000/ascending` | 231 ns | 242 ns | 10.6 µs |  |  |  | 46x faster |
| `linear/4d/max10/eval1000/flipped` | 361 ns | 368 ns | 18.6 µs |  |  |  | 51x faster |
| `linear/4d/max10/eval10000/ascending` | 231 ns | 236 ns | 9.76 µs |  |  |  | 42x faster |
| `linear/4d/max10/eval10000/flipped` | 364 ns | 410 ns | 18.1 µs |  |  |  | 50x faster |
| `linear/4d/max40/eval10/ascending` | 14.4 µs | 40.6 µs | 9.48 µs | 5.37 µs |  | 1.27 ns | **1.1e+04x slower** |
| `linear/4d/max40/eval10/flipped` | 57.2 µs | 69 µs | 21.5 µs |  |  |  | **2.7x slower** |
| `linear/4d/max40/eval100/ascending` | 13.9 µs | 37.3 µs | 9.55 µs | 5.59 µs |  | 1.37 ns | **1e+04x slower** |
| `linear/4d/max40/eval100/flipped` | 54.8 µs | 66.8 µs | 21.2 µs |  |  |  | **2.6x slower** |
| `linear/4d/max40/eval1000/ascending` | 14.4 µs | 33.5 µs | 9.5 µs | 5.67 µs |  | 1.37 ns | **1.1e+04x slower** |
| `linear/4d/max40/eval1000/flipped` | 52.2 µs | 66.1 µs | 19.8 µs |  |  |  | **2.6x slower** |
| `linear/4d/max40/eval10000/ascending` | 13.3 µs | 34.2 µs | 9.4 µs | 5.52 µs |  | 1.27 ns | **1.1e+04x slower** |
| `linear/4d/max40/eval10000/flipped` | 51.9 µs | 70.1 µs | 21.5 µs |  |  |  | **2.4x slower** |
| `linear/4d/max200/eval10/ascending` | 51.7 ms | 18 ms | 9.42 µs | 43.1 ms |  | 1.27 ns | **4.1e+07x slower** |
| `linear/4d/max200/eval10/flipped` | 122 ms | 87.1 ms | 20.3 µs |  |  |  | **6e+03x slower** |
| `linear/4d/max200/eval100/ascending` | 57 ms | 19 ms | 10.3 µs | 43.6 ms |  | 1.27 ns | **4.5e+07x slower** |
| `linear/4d/max200/eval100/flipped` | 120 ms | 82.1 ms | 20.2 µs |  |  |  | **6e+03x slower** |
| `linear/4d/max200/eval1000/ascending` | 55 ms | 18.6 ms | 9.46 µs | 45.4 ms |  | 1.37 ns | **4e+07x slower** |
| `linear/4d/max200/eval1000/flipped` | 122 ms | 84.4 ms | 20.2 µs |  |  |  | **6.1e+03x slower** |
| `linear/4d/max200/eval10000/ascending` | 57.6 ms | 18.1 ms | 9.59 µs | 44.9 ms |  | 1.17 ns | **4.9e+07x slower** |
| `linear/4d/max200/eval10000/flipped` | 123 ms | 84.3 ms | 22 µs |  |  |  | **5.6e+03x slower** |
| `nearest/2d/each40/eval10000/ascending` | 615 ns | 1.25 µs | 5.82 µs | 293 ns |  | 0 | **9e+02x slower** |
| `nearest/3d/each40/eval10000/ascending` | 23 µs | 62.8 µs | 7.5 µs | 9.22 µs |  | 0 | **2.4e+04x slower** |
| `pchip/2d/each40/eval10000/ascending` | 6.48 µs | 6.33 µs | 6.74 µs |  |  |  | 1x faster |
| `pchip/3d/each40/eval10000/ascending` | 507 µs | 456 µs | 9.35 µs |  |  |  | **54x slower** |
| `cubic/2d/each40/eval10000/ascending` | 9.26 µs | 8.84 µs | 1.9 ms |  |  | 0 | **1.1e+04x slower** |
| `cubic/3d/each40/eval10000/ascending` | 431 µs | 288 µs | 174 ms |  |  | 0 | **4.4e+05x slower** |
| `cubic/2d/each10/eval1000/ascending` | 1.62 µs | 1.59 µs | 987 µs |  |  | 0 | **1.8e+03x slower** |
| `cubic/2d/each40/eval1000/ascending` | 9.22 µs | 8.9 µs | 1.91 ms |  |  | 0 | **1.2e+04x slower** |
| `cubic/2d/each100/eval1000/ascending` | 84 µs | 84.7 µs | 7.84 ms |  |  | 0 | **9.6e+04x slower** |
| `cubic/2d/each200/eval1000/ascending` | 187 µs | 224 µs | 30.2 ms |  |  | 0 | **2.1e+05x slower** |
| `cubic/2d/each400/eval1000/ascending` | 1.64 ms | 1.28 ms | 157 ms |  |  | 0 | **1.9e+06x slower** |
| `quintic/2d/each10/eval1000/ascending` | 2.14 µs | 2.03 µs | 3.64 ms |  |  |  | 1.7e+03x faster |
| `quintic/2d/each40/eval1000/ascending` | 12.6 µs | 12.2 µs | 7.95 ms |  |  |  | 6.3e+02x faster |
| `linear/2d/each400/eval10000/ascending` | 373 µs | 107 µs | 6.01 µs | 315 µs | 1.25 µs | 0 | **5.4e+05x slower** |
| `linear/2d/each800/eval10000/ascending` | 1.3 ms | 511 µs | 6.33 µs | 860 µs | 2.23 µs | 0 | **1.9e+06x slower** |

## `NearestNDInterpolator`

`n` points of the unit square, evaluated at the pixels of a square image.

### Evaluation

| case | interpx | interpx, `from_owned` | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `n10/image50` | 113 µs | 112 µs | 446 µs | 3.9x faster |
| `n10/image100` | 488 µs | 487 µs | 1.7 ms | 3.5x faster |
| `n10/image500` | 11.6 ms | 11.9 ms | 45.4 ms | 3.9x faster |
| `n100/image50` | 189 µs | 188 µs | 578 µs | 3.1x faster |
| `n100/image100` | 786 µs | 722 µs | 2.19 ms | 2.8x faster |
| `n100/image500` | 18 ms | 18.1 ms | 55.7 ms | 3.1x faster |
| `n1000/image50` | 323 µs | 325 µs | 796 µs | 2.5x faster |
| `n1000/image100` | 1.16 ms | 1.17 ms | 2.95 ms | 2.5x faster |
| `n1000/image500` | 24.9 ms | 23.2 ms | 70.1 ms | 2.8x faster |
| `n10000/image50` | 336 µs | 336 µs | 983 µs | 2.9x faster |
| `n10000/image100` | 1.48 ms | 1.4 ms | 3.85 ms | 2.6x faster |
| `n10000/image500` | 31.6 ms | 32 ms | 85.9 ms | 2.7x faster |
| `n100000/image50` | 478 µs | 477 µs | 1.4 ms | 2.9x faster |
| `n100000/image100` | 1.83 ms | 1.94 ms | 5.63 ms | 3.1x faster |
| `n100000/image500` | 44.7 ms | 45.7 ms | 123 ms | 2.7x faster |

### Construction

| case | interpx | interpx, `from_owned` | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `n10/image50` | 528 ns | 551 ns | 9.48 µs | 18x faster |
| `n10/image100` | 564 ns | 568 ns | 9.52 µs | 17x faster |
| `n10/image500` | 522 ns | 559 ns | 10.3 µs | 20x faster |
| `n100/image50` | 3.36 µs | 3.6 µs | 15 µs | 4.5x faster |
| `n100/image100` | 3.48 µs | 3.34 µs | 16.4 µs | 4.7x faster |
| `n100/image500` | 3.37 µs | 3.65 µs | 15.2 µs | 4.5x faster |
| `n1000/image50` | 59 µs | 47.7 µs | 126 µs | 2.1x faster |
| `n1000/image100` | 58.6 µs | 42.8 µs | 118 µs | 2x faster |
| `n1000/image500` | 62 µs | 51.3 µs | 128 µs | 2.1x faster |
| `n10000/image50` | 730 µs | 650 µs | 1.79 ms | 2.5x faster |
| `n10000/image100` | 677 µs | 658 µs | 1.67 ms | 2.5x faster |
| `n10000/image500` | 639 µs | 682 µs | 1.78 ms | 2.8x faster |
| `n100000/image50` | 9.15 ms | 8.42 ms | 27 ms | 3x faster |
| `n100000/image100` | 8.78 ms | 8.78 ms | 25.8 ms | 2.9x faster |
| `n100000/image500` | 8.91 ms | 8.74 ms | 26.2 ms | 2.9x faster |

## `RBFInterpolator`

`n` points of a square, evaluated at 1000 points in no order (`scattered`) or at the 10 000 pixels of an image (`image`). `global`: all points; `k50`: the 50 nearest.

### Evaluation

| case | interpx | interpx, `from_owned` | SciPy | scirs2-interpolate | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: |
| `linear/global/n10/scattered` | 15.8 µs | 15.7 µs | 63 µs | 197 µs | 4x faster |
| `linear/k50/n10/scattered` | 15.2 µs | 15.7 µs | 1.39 ms |  | 92x faster |
| `thin_plate_spline/global/n10/scattered` | 73.1 µs | 79.8 µs | 111 µs | 297 µs | 1.5x faster |
| `thin_plate_spline/k50/n10/scattered` | 72.9 µs | 73.7 µs | 1.43 ms |  | 20x faster |
| `cubic/global/n10/scattered` | 17.3 µs | 17.2 µs | 60.1 µs | 217 µs | 3.5x faster |
| `cubic/k50/n10/scattered` | 17 µs | 18.6 µs | 1.43 ms |  | 84x faster |
| `quintic/global/n10/scattered` | 25.2 µs | 25.9 µs | 69.6 µs | 219 µs | 2.8x faster |
| `quintic/k50/n10/scattered` | 26.4 µs | 26.5 µs | 1.48 ms |  | 56x faster |
| `multiquadric/global/n10/scattered` | 30.5 µs | 28.2 µs | 52.8 µs | 236 µs | 1.7x faster |
| `multiquadric/k50/n10/scattered` | 27.8 µs | 27.8 µs | 1.42 ms |  | 51x faster |
| `inverse_multiquadric/global/n10/scattered` | 37.4 µs | 37.4 µs | 60.1 µs | 291 µs | 1.6x faster |
| `inverse_multiquadric/k50/n10/scattered` | 41.1 µs | 41.1 µs | 1.36 ms |  | 33x faster |
| `inverse_quadratic/global/n10/scattered` | 23.8 µs | 26.2 µs | 50.9 µs |  | 2.1x faster |
| `inverse_quadratic/k50/n10/scattered` | 26.2 µs | 26.2 µs | 1.51 ms |  | 58x faster |
| `gaussian/global/n10/scattered` | 67.6 µs | 74 µs | 95.2 µs | 354 µs | 1.4x faster |
| `gaussian/k50/n10/scattered` | 68.1 µs | 74 µs | 1.37 ms |  | 20x faster |
| `linear/global/n100/scattered` | 148 µs | 135 µs | 293 µs | 1.86 ms | 2x faster |
| `linear/k50/n100/scattered` | 13.4 ms | 13.8 ms | 20.1 ms |  | 1.5x faster |
| `thin_plate_spline/global/n100/scattered` | 710 µs | 768 µs | 907 µs | 3.12 ms | 1.3x faster |
| `thin_plate_spline/k50/n100/scattered` | 24.8 ms | 23.6 ms | 24.5 ms |  | **1x slower** |
| `cubic/global/n100/scattered` | 140 µs | 152 µs | 292 µs | 1.95 ms | 2.1x faster |
| `cubic/k50/n100/scattered` | 15.3 ms | 15.3 ms | 21.4 ms |  | 1.4x faster |
| `quintic/global/n100/scattered` | 167 µs | 172 µs | 308 µs | 2.25 ms | 1.8x faster |
| `quintic/k50/n100/scattered` | 16.9 ms | 16.9 ms | 23.1 ms |  | 1.4x faster |
| `multiquadric/global/n100/scattered` | 295 µs | 269 µs | 303 µs | 2.46 ms | 1x faster |
| `multiquadric/k50/n100/scattered` | 15.6 ms | 15.6 ms | 20.2 ms |  | 1.3x faster |
| `inverse_multiquadric/global/n100/scattered` | 364 µs | 399 µs | 413 µs | 2.61 ms | 1.1x faster |
| `inverse_multiquadric/k50/n100/scattered` | 17.1 ms | 15.9 ms | 20 ms |  | 1.2x faster |
| `inverse_quadratic/global/n100/scattered` | 224 µs | 224 µs | 268 µs |  | 1.2x faster |
| `inverse_quadratic/k50/n100/scattered` | 14.8 ms | 14.6 ms | 20.5 ms |  | 1.4x faster |
| `gaussian/global/n100/scattered` | 675 µs | 640 µs | 686 µs | 3.8 ms | 1x faster |
| `gaussian/k50/n100/scattered` | 21.5 ms | 21.6 ms | 22.7 ms |  | 1.1x faster |
| `linear/global/n1000/scattered` | 1.47 ms | 1.47 ms | 4.07 ms | 18.6 ms | 2.8x faster |
| `linear/k50/n1000/scattered` | 14.3 ms | 14.7 ms | 34.7 ms |  | 2.4x faster |
| `thin_plate_spline/global/n1000/scattered` | 7.73 ms | 7.77 ms | 10.8 ms | 31.3 ms | 1.4x faster |
| `thin_plate_spline/k50/n1000/scattered` | 23.5 ms | 23.4 ms | 46.1 ms |  | 2x faster |
| `cubic/global/n1000/scattered` | 1.48 ms | 1.47 ms | 4.14 ms | 20.7 ms | 2.8x faster |
| `cubic/k50/n1000/scattered` | 16.9 ms | 17.3 ms | 38.4 ms |  | 2.3x faster |
| `quintic/global/n1000/scattered` | 1.65 ms | 1.65 ms | 4.34 ms | 22.8 ms | 2.6x faster |
| `quintic/k50/n1000/scattered` | 18.3 ms | 18.4 ms | 42.7 ms |  | 2.3x faster |
| `multiquadric/global/n1000/scattered` | 2.95 ms | 2.81 ms | 4.51 ms | 24.7 ms | 1.5x faster |
| `multiquadric/k50/n1000/scattered` | 17.5 ms | 17.7 ms | 37.4 ms |  | 2.1x faster |
| `inverse_multiquadric/global/n1000/scattered` | 4 ms | 4 ms | 5.41 ms | 28.1 ms | 1.4x faster |
| `inverse_multiquadric/k50/n1000/scattered` | 18.7 ms | 18.1 ms | 38.3 ms |  | 2x faster |
| `inverse_quadratic/global/n1000/scattered` | 2.45 ms | 2.45 ms | 4.01 ms |  | 1.6x faster |
| `inverse_quadratic/k50/n1000/scattered` | 16.6 ms | 16.5 ms | 36.1 ms |  | 2.2x faster |
| `gaussian/global/n1000/scattered` | 6.72 ms | 6.21 ms | 8.64 ms | 37.5 ms | 1.3x faster |
| `gaussian/k50/n1000/scattered` | 23 ms | 22.8 ms | 43 ms |  | 1.9x faster |
| `thin_plate_spline/global/n3000/scattered` | 23.2 ms | 22 ms | 31.2 ms | 92.8 ms | 1.3x faster |
| `thin_plate_spline/k50/n10000/scattered` | 23.9 ms | 22.5 ms | 46.1 ms |  | 1.9x faster |
| `thin_plate_spline/k50/n10000/image` | 235 ms | 238 ms | 456 ms |  | 1.9x faster |
| `thin_plate_spline/k50/n100000/scattered` | 25.3 ms | 25.1 ms | 48.3 ms |  | 1.9x faster |
| `thin_plate_spline/k50/n100000/image` | 242 ms | 239 ms | 471 ms |  | 1.9x faster |

### Construction

| case | interpx | interpx, `from_owned` | SciPy | scirs2-interpolate | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: |
| `linear/global/n10/scattered` | 1.34 µs | 1.33 µs | 11.7 µs | 5.89 µs | 4.4x faster |
| `linear/k50/n10/scattered` | 1.25 µs | 1.24 µs | 15.9 µs |  | 13x faster |
| `thin_plate_spline/global/n10/scattered` | 1.87 µs | 1.91 µs | 13.4 µs | 7.14 µs | 3.8x faster |
| `thin_plate_spline/k50/n10/scattered` | 1.83 µs | 1.86 µs | 17.1 µs |  | 9.3x faster |
| `cubic/global/n10/scattered` | 1.58 µs | 1.68 µs | 12.2 µs | 6.03 µs | 3.8x faster |
| `cubic/k50/n10/scattered` | 1.71 µs | 1.65 µs | 16.3 µs |  | 9.5x faster |
| `quintic/global/n10/scattered` | 2.21 µs | 2.16 µs | 15 µs | 5.79 µs | 2.6x faster |
| `quintic/k50/n10/scattered` | 2.07 µs | 2.08 µs | 18.5 µs |  | 8.9x faster |
| `multiquadric/global/n10/scattered` | 1.42 µs | 1.37 µs | 10.7 µs | 6.73 µs | 4.7x faster |
| `multiquadric/k50/n10/scattered` | 1.34 µs | 1.41 µs | 15.4 µs |  | 11x faster |
| `inverse_multiquadric/global/n10/scattered` | 1.33 µs | 1.33 µs | 11.7 µs | 7.53 µs | 5.7x faster |
| `inverse_multiquadric/k50/n10/scattered` | 1.42 µs | 1.43 µs | 14.9 µs |  | 11x faster |
| `inverse_quadratic/global/n10/scattered` | 1.35 µs | 1.32 µs | 11.7 µs |  | 8.7x faster |
| `inverse_quadratic/k50/n10/scattered` | 1.24 µs | 1.3 µs | 16.2 µs |  | 13x faster |
| `gaussian/global/n10/scattered` | 1.52 µs | 1.41 µs | 11.9 µs | 8.47 µs | 5.6x faster |
| `gaussian/k50/n10/scattered` | 1.56 µs | 1.55 µs | 16.3 µs |  | 10x faster |
| `linear/global/n100/scattered` | 52.5 µs | 58.6 µs | 77.1 µs | 780 µs | 1.5x faster |
| `linear/k50/n100/scattered` | 3.49 µs | 3.48 µs | 22.6 µs |  | 6.5x faster |
| `thin_plate_spline/global/n100/scattered` | 78.4 µs | 113 µs | 112 µs | 918 µs | 1.4x faster |
| `thin_plate_spline/k50/n100/scattered` | 3.85 µs | 3.69 µs | 22.6 µs |  | 5.9x faster |
| `cubic/global/n100/scattered` | 46.6 µs | 70.8 µs | 87.5 µs | 772 µs | 1.9x faster |
| `cubic/k50/n100/scattered` | 3.36 µs | 3.53 µs | 21.8 µs |  | 6.5x faster |
| `quintic/global/n100/scattered` | 52.5 µs | 80.8 µs | 97.3 µs | 806 µs | 1.9x faster |
| `quintic/k50/n100/scattered` | 3.57 µs | 3.53 µs | 24.9 µs |  | 7x faster |
| `multiquadric/global/n100/scattered` | 51.4 µs | 63.7 µs | 76.9 µs | 702 µs | 1.5x faster |
| `multiquadric/k50/n100/scattered` | 3.61 µs | 3.46 µs | 20.9 µs |  | 5.8x faster |
| `inverse_multiquadric/global/n100/scattered` | 55.7 µs | 66.9 µs | 80.9 µs | 807 µs | 1.5x faster |
| `inverse_multiquadric/k50/n100/scattered` | 3.48 µs | 3.53 µs | 21.1 µs |  | 6.1x faster |
| `inverse_quadratic/global/n100/scattered` | 52.8 µs | 55.3 µs | 73.9 µs |  | 1.4x faster |
| `inverse_quadratic/k50/n100/scattered` | 3.53 µs | 3.61 µs | 20.9 µs |  | 5.9x faster |
| `gaussian/global/n100/scattered` | 72.9 µs | 89.9 µs | 108 µs | 838 µs | 1.5x faster |
| `gaussian/k50/n100/scattered` | 3.5 µs | 3.46 µs | 21.3 µs |  | 6.1x faster |
| `linear/global/n1000/scattered` | 16.8 ms | 18.2 ms | 61.9 ms | 194 ms | 3.7x faster |
| `linear/k50/n1000/scattered` | 67.3 µs | 60.8 µs | 138 µs |  | 2.1x faster |
| `thin_plate_spline/global/n1000/scattered` | 23 ms | 21.1 ms | 54.5 ms | 218 ms | 2.4x faster |
| `thin_plate_spline/k50/n1000/scattered` | 66.3 µs | 58.8 µs | 140 µs |  | 2.1x faster |
| `cubic/global/n1000/scattered` | 18.8 ms | 19 ms | 64.3 ms | 196 ms | 3.4x faster |
| `cubic/k50/n1000/scattered` | 63.6 µs | 61.2 µs | 140 µs |  | 2.2x faster |
| `quintic/global/n1000/scattered` | 18.9 ms | 19.1 ms | 50.7 ms | 203 ms | 2.7x faster |
| `quintic/k50/n1000/scattered` | 66.5 µs | 58.8 µs | 143 µs |  | 2.2x faster |
| `multiquadric/global/n1000/scattered` | 19.3 ms | 19.6 ms | 57.1 ms | 198 ms | 3x faster |
| `multiquadric/k50/n1000/scattered` | 65.7 µs | 60 µs | 138 µs |  | 2.1x faster |
| `inverse_multiquadric/global/n1000/scattered` | 19.4 ms | 20.3 ms | 63.5 ms | 210 ms | 3.3x faster |
| `inverse_multiquadric/k50/n1000/scattered` | 65 µs | 60.7 µs | 137 µs |  | 2.1x faster |
| `inverse_quadratic/global/n1000/scattered` | 19.2 ms | 18.7 ms | 59.3 ms |  | 3.1x faster |
| `inverse_quadratic/k50/n1000/scattered` | 64.5 µs | 55.7 µs | 140 µs |  | 2.2x faster |
| `gaussian/global/n1000/scattered` | 20.8 ms | 21.2 ms | 62.2 ms | 212 ms | 3x faster |
| `gaussian/k50/n1000/scattered` | 64.3 µs | 60.9 µs | 129 µs |  | 2x faster |
| `thin_plate_spline/global/n3000/scattered` | 442 ms | 426 ms | 2.48 s | 7.17 s | 5.6x faster |
| `thin_plate_spline/k50/n10000/scattered` | 811 µs | 735 µs | 1.77 ms |  | 2.2x faster |
| `thin_plate_spline/k50/n10000/image` | 839 µs | 767 µs | 1.65 ms |  | 2x faster |
| `thin_plate_spline/k50/n100000/scattered` | 9.26 ms | 9.19 ms | 26.4 ms |  | 2.9x faster |
| `thin_plate_spline/k50/n100000/image` | 9.95 ms | 8.92 ms | 26.6 ms |  | 2.7x faster |

## Memory

Megabytes. "Data read" is the process with the data of the case in memory, before anything is built. The columns after it are what was added to that: the highest while building, with the interpolator built, the highest while evaluating. "Highest" is the highest of the whole run, repeats included.

| case | library | data read | while building | built | while evaluating | highest |
| :--- | :--- | ---: | ---: | ---: | ---: | ---: |
| `oned/cubic_spline/n10000/sorted` | interpx | 4.0 | +1.2 | +1.2 | +1.2 | 5.9 |
| `oned/cubic_spline/n10000/sorted` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |
| `oned/cubic_spline/n10000/sorted` | scirs2-interpolate | 4.7 | +0.4 | +0.5 | +1.2 | 5.9 |
| `oned/cubic_spline/n10000/sorted` | ndarray-interp | 4.2 | +0.7 | +0.7 | +0.7 | 5.4 |
| `grid/linear/4d/max200/eval10000/ascending` | interpx | 227.8 | +177.2 | +176.9 | +176.9 | 405.0 |
| `grid/linear/4d/max200/eval10000/ascending` | interpx, `from_owned` | 204.6 | +200.1 | +200.1 | +200.1 | 404.8 |
| `grid/linear/4d/max200/eval10000/ascending` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |
| `grid/linear/4d/max200/eval10000/ascending` | scirs2-interpolate | 228.8 | +170.0 | +176.1 | +176.1 | 404.9 |
| `grid/linear/4d/max200/eval10000/ascending` | interpn | 204.7 | +0.0 | +0.0 | +0.0 | 403.3 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpx | 5.6 | +0.2 | +2.7 | +2.7 | 8.4 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpx, `from_owned` | 4.3 | +4.1 | +4.0 | +4.0 | 8.4 |
| `grid/cubic/2d/each400/eval1000/ascending` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpn | 5.6 | +0.0 | +0.0 | +0.0 | 5.6 |
| `grid/linear/2d/each800/eval10000/ascending` | interpx | 9.7 | +5.2 | +5.2 | +5.2 | 14.9 |
| `grid/linear/2d/each800/eval10000/ascending` | interpx, `from_owned` | 13.6 | +1.1 | +1.1 | +1.1 | 14.9 |
| `grid/linear/2d/each800/eval10000/ascending` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |
| `grid/linear/2d/each800/eval10000/ascending` | scirs2-interpolate | 9.5 | +0.2 | +5.3 | +5.4 | 14.9 |
| `grid/linear/2d/each800/eval10000/ascending` | ndarray-interp | 11.4 | +0.0 | +0.0 | +0.0 | 11.4 |
| `grid/linear/2d/each800/eval10000/ascending` | interpn | 9.9 | +0.0 | +0.0 | +0.0 | 9.9 |
| `nearest_nd/n100000/image500` | interpx | 15.2 | +0.2 | +0.0 | +0.9 | 18.2 |
| `nearest_nd/n100000/image500` | interpx, `from_owned` | 13.8 | +2.0 | +0.7 | +2.7 | 18.6 |
| `nearest_nd/n100000/image500` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |
| `rbf/thin_plate_spline/global/n3000/scattered` | interpx | 3.9 | +76.3 | +4.1 | +4.1 | 80.3 |
| `rbf/thin_plate_spline/global/n3000/scattered` | interpx, `from_owned` | 4.2 | +76.1 | +3.9 | +3.9 | 80.3 |
| `rbf/thin_plate_spline/global/n3000/scattered` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |
| `rbf/thin_plate_spline/global/n3000/scattered` | scirs2-interpolate | 3.3 | +289.2 | +1.2 | +1.2 | 292.5 |
| `rbf/thin_plate_spline/k50/n100000/image` | interpx | 6.7 | +3.8 | +4.7 | +5.0 | 11.9 |
| `rbf/thin_plate_spline/k50/n100000/image` | interpx, `from_owned` | 6.7 | +5.9 | +4.7 | +5.0 | 13.5 |
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

| family | interpx | interpx, `from_owned` | scirs2-interpolate | ndarray-interp | interpn |
| :--- | ---: | ---: | ---: | ---: | ---: |
| One variable | 3e-16 |  | 6e-02 | 1e-04 |  |
| `PPoly` | 9e-16 |  |  |  |  |
| `BSpline` | 4e-16 |  |  |  |  |
| `BSplineCurve` | 4e-16 |  |  |  |  |
| Regular grid | 4e-04 | 4e-04 | 4e-16 | 6e-16 | 7e-01 |
| `NearestNDInterpolator` | 0e+00 | 0e+00 |  |  |  |
| `RBFInterpolator` | 8e-05 | 8e-05 | 3e+01 |  |  |
