# Performance: interpx against SciPy and other Rust crates

<!-- Written by compare/perf/report.py from a run of compare/perf/run.py. Do not edit by hand: run it again. -->

Measured on 2026-10-13, at commit `8d62d6ccbd` with uncommitted changes.

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
| One variable | evaluation | 100 | 96 | 0 | 4 |
| One variable | construction | 100 | 77 | 1 | 23 |
| `PPoly` | evaluation | 4 | 4 | 0 | 0 |
| `PPoly` | construction | 4 | 4 | 0 | 0 |
| `BSpline` | evaluation | 4 | 4 | 0 | 0 |
| `BSpline` | construction | 4 | 4 | 0 | 0 |
| `BSplineCurve` | evaluation | 2 | 2 | 0 | 0 |
| `BSplineCurve` | construction | 2 | 2 | 0 | 0 |
| Regular grid | evaluation | 79 | 56 | 0 | 23 |
| Regular grid | construction | 79 | 20 | 34 | 43 |
| `NearestNDInterpolator` | evaluation | 15 | 15 | 0 | 0 |
| `NearestNDInterpolator` | construction | 15 | 15 | 0 | 0 |
| `RBFInterpolator` | evaluation | 53 | 52 | 1 | 0 |
| `RBFInterpolator` | construction | 53 | 53 | 0 | 0 |

180 comparisons (a case, an operation, and a library), by method. "Same values" says whether that library's values are SciPy's up to rounding: where they are not, it computes another interpolant under the same name, and its time is not that of the same work.

| method | operation | against | same values | cases | slower by | the largest |
| :--- | :--- | :--- | :--- | ---: | ---: | :--- |
| `grid/linear` | construction | interpn | yes | 34 | 171x to 4.26e+07x | `grid/linear/4d/max200/eval1000/ascending`: 49.9 ms against 1.17 ns |
| `grid/cubic` | construction | interpn | no | 7 | 1.68e+03x to 1.6e+06x | `grid/cubic/2d/each400/eval1000/ascending`: 1.4 ms against 0 |
| `grid/nearest` | construction | interpn | yes | 2 | 821x to 2.46e+04x | `grid/nearest/3d/each40/eval10000/ascending`: 21.6 µs against 0 |
| `grid/linear` | construction | SciPy | yes | 32 | 1.14x to 5.27e+03x | `grid/linear/4d/max200/eval10000/flipped`: 107 ms against 20.3 µs |
| `grid/linear` | construction | ndarray-interp | yes | 14 | 1.36x to 451x | `grid/linear/2d/each800/eval10000/ascending`: 1.1 ms against 2.44 µs |
| `grid/pchip` | construction | SciPy | yes | 1 | 53.9x | `grid/pchip/3d/each40/eval10000/ascending`: 472 µs against 8.76 µs |
| `oned/linear` | construction | ndarray-interp | yes | 10 | 1.11x to 5.54x | `oned/linear/n10000/shuffled`: 101 µs against 18.2 µs |
| `grid/nearest` | construction | SciPy | yes | 1 | 2.83x | `grid/nearest/3d/each40/eval10000/ascending`: 21.6 µs against 7.63 µs |
| `oned/linear` | construction | scirs2-interpolate | yes | 9 | 1.12x to 2.74x | `oned/linear/n10000/shuffled`: 101 µs against 36.9 µs |
| `grid/linear` | construction | scirs2-interpolate | yes | 26 | 1.18x to 2.52x | `grid/linear/4d/max40/eval100/ascending`: 14.2 µs against 5.64 µs |
| `grid/nearest` | construction | scirs2-interpolate | yes | 2 | 1.89x to 2.44x | `grid/nearest/3d/each40/eval10000/ascending`: 21.6 µs against 8.86 µs |
| `oned/linear` | evaluation | ndarray-interp | yes | 2 | 1.51x to 2.31x | `oned/linear/n10000/shuffled`: 642 µs against 278 µs |
| `oned/akima` | construction | scirs2-interpolate | yes | 6 | 1.52x to 1.96x | `oned/akima/n1000/shuffled`: 9.66 µs against 4.94 µs |
| `oned/pchip` | construction | scirs2-interpolate | yes | 6 | 1.35x to 1.92x | `oned/pchip/n1000/sorted`: 8.85 µs against 4.61 µs |
| `oned/linear` | construction | SciPy | yes | 1 | 1.9x | `oned/linear/n10000/sorted`: 93.7 µs against 49.3 µs |
| `oned/cubic_spline` | evaluation | ndarray-interp | yes | 2 | 1.23x to 1.84x | `oned/cubic_spline/n10000/shuffled`: 698 µs against 379 µs |
| `grid/linear` | evaluation | interpn | yes | 22 | 1.12x to 1.66x | `grid/linear/4d/max40/eval1000/ascending`: 45.9 µs against 27.6 µs |
| `rbf/quintic` | evaluation | SciPy | yes | 1 | 1.12x | `rbf/quintic/global/n10/scattered`: 75 µs against 66.8 µs |
| `oned/cubic_spline` | construction | scirs2-interpolate | no | 1 | 1.1x | `oned/cubic_spline/n1000/sorted`: 19 µs against 17.2 µs |
| `grid/nearest` | evaluation | interpn | yes | 1 | 1.1x | `grid/nearest/2d/each40/eval10000/ascending`: 153 µs against 139 µs |

<details><summary>Every comparison where interpx is slower</summary>

"Values apart" is the distance of that library's values from SciPy's, as a share of the largest value.

| case | operation | against | interpx | the other | slower by | values apart |
| :--- | :--- | :--- | ---: | ---: | ---: | ---: |
| `grid/linear/4d/max200/eval1000/ascending` | construction | interpn | 49.9 ms | 1.17 ns | 42583637.33x | 3e-16 |
| `grid/linear/4d/max200/eval100/ascending` | construction | interpn | 48.7 ms | 1.17 ns | 41588437.33x | 3e-16 |
| `grid/linear/4d/max200/eval10000/ascending` | construction | interpn | 48.6 ms | 1.17 ns | 41450325.33x | 4e-16 |
| `grid/linear/4d/max200/eval10/ascending` | construction | interpn | 47.4 ms | 1.17 ns | 40475861.33x | 2e-16 |
| `grid/linear/3d/max200/eval10/ascending` | construction | interpn | 1.9 ms | 0 | 2160412.44x | 3e-16 |
| `grid/linear/3d/max200/eval10000/ascending` | construction | interpn | 1.85 ms | 0 | 2107960.89x | 5e-16 |
| `grid/linear/3d/max200/eval100/ascending` | construction | interpn | 1.96 ms | 0 | 2002124.80x | 3e-16 |
| `grid/linear/3d/max200/eval1000/ascending` | construction | interpn | 1.85 ms | 1.07 ns | 1719435.64x | 3e-16 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | interpn | 1.1 ms | 0 | 1607753.14x | 3e-16 |
| `grid/cubic/2d/each400/eval1000/ascending` | construction | interpn | 1.4 ms | 0 | 1595960.89x | 5e-01 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | interpn | 327 µs | 0 | 478573.71x | 3e-16 |
| `grid/cubic/3d/each40/eval10000/ascending` | construction | interpn | 306 µs | 0 | 313702.40x | 5e-01 |
| `grid/cubic/2d/each200/eval1000/ascending` | construction | interpn | 176 µs | 0 | 180121.60x | 2e-01 |
| `grid/cubic/2d/each100/eval1000/ascending` | construction | interpn | 43.3 µs | 0 | 49294.22x | 2e-01 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | interpn | 21.6 µs | 0 | 24586.67x | 0e+00 |
| `grid/linear/4d/max40/eval1000/ascending` | construction | interpn | 14.3 µs | 1.17 ns | 12165.33x | 6e-16 |
| `grid/linear/4d/max40/eval100/ascending` | construction | interpn | 14.2 µs | 1.17 ns | 12106.67x | 2e-16 |
| `grid/linear/4d/max40/eval10000/ascending` | construction | interpn | 13.3 µs | 1.17 ns | 11392.00x | 5e-16 |
| `grid/linear/2d/max200/eval10000/ascending` | construction | interpn | 7.58 µs | 0 | 11081.14x | 3e-16 |
| `grid/linear/4d/max40/eval10/ascending` | construction | interpn | 13.9 µs | 1.27 ns | 10976.00x | 2e-16 |
| `grid/cubic/2d/each40/eval10000/ascending` | construction | interpn | 8.55 µs | 0 | 10944.00x | 7e-01 |
| `grid/linear/2d/max200/eval100/ascending` | construction | interpn | 6.83 µs | 0 | 9997.71x | 1e-16 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | interpn | 6.92 µs | 0 | 8864.00x | 2e-16 |
| `grid/linear/2d/max200/eval10/ascending` | construction | interpn | 6.9 µs | 0 | 8832.00x | 2e-16 |
| `grid/cubic/2d/each40/eval1000/ascending` | construction | interpn | 8.36 µs | 0 | 8563.20x | 7e-01 |
| `grid/linear/4d/max200/eval10000/flipped` | construction | SciPy | 107 ms | 20.3 µs | 5268.01x |  |
| `grid/linear/4d/max200/eval10000/ascending` | construction | SciPy | 48.6 ms | 9.48 µs | 5126.61x |  |
| `grid/linear/4d/max200/eval10/flipped` | construction | SciPy | 99.9 ms | 20.3 µs | 4927.61x |  |
| `grid/linear/4d/max200/eval1000/ascending` | construction | SciPy | 49.9 ms | 10.2 µs | 4885.69x |  |
| `grid/linear/4d/max200/eval10/ascending` | construction | SciPy | 47.4 ms | 9.86 µs | 4809.39x |  |
| `grid/linear/4d/max200/eval100/flipped` | construction | SciPy | 100 ms | 21.3 µs | 4692.58x |  |
| `grid/linear/4d/max200/eval100/ascending` | construction | SciPy | 48.7 ms | 10.5 µs | 4660.99x |  |
| `grid/linear/4d/max200/eval1000/flipped` | construction | SciPy | 99.7 ms | 22 µs | 4540.71x |  |
| `grid/linear/3d/max40/eval1000/ascending` | construction | interpn | 3.13 µs | 0 | 3555.56x | 2e-16 |
| `grid/linear/3d/max40/eval10/ascending` | construction | interpn | 2.96 µs | 0 | 3367.11x | 2e-16 |
| `grid/linear/3d/max40/eval10000/ascending` | construction | interpn | 2.76 µs | 0 | 3136.00x | 3e-16 |
| `grid/linear/3d/max40/eval100/ascending` | construction | interpn | 3.05 µs | 0 | 3118.40x | 2e-16 |
| `grid/cubic/2d/each10/eval1000/ascending` | construction | interpn | 1.48 µs | 0 | 1680.89x | 3e-01 |
| `grid/nearest/2d/each40/eval10000/ascending` | construction | interpn | 561 ns | 0 | 821.00x | 0e+00 |
| `grid/linear/2d/max40/eval10/ascending` | construction | interpn | 367 ns | 0 | 536.86x | 2e-16 |
| `grid/linear/2d/max40/eval1000/ascending` | construction | interpn | 364 ns | 0 | 533.14x | 3e-16 |
| `grid/linear/2d/max40/eval100/ascending` | construction | interpn | 354 ns | 0 | 517.57x | 2e-16 |
| `grid/linear/2d/max40/eval10000/ascending` | construction | interpn | 370 ns | 0 | 473.75x | 4e-16 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | ndarray-interp | 1.1 ms | 2.44 µs | 451.04x | 6e-16 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | ndarray-interp | 327 µs | 1.15 µs | 284.48x | 5e-16 |
| `grid/linear/3d/max200/eval100/ascending` | construction | SciPy | 1.96 ms | 7.71 µs | 253.72x |  |
| `grid/linear/3d/max200/eval10/ascending` | construction | SciPy | 1.9 ms | 7.68 µs | 247.30x |  |
| `grid/linear/3d/max200/eval100/flipped` | construction | SciPy | 3.84 ms | 16.4 µs | 233.67x |  |
| `grid/linear/3d/max200/eval10000/ascending` | construction | SciPy | 1.85 ms | 8.29 µs | 223.39x |  |
| `grid/linear/3d/max200/eval1000/ascending` | construction | SciPy | 1.85 ms | 8.41 µs | 219.52x |  |
| `grid/linear/3d/max200/eval10000/flipped` | construction | SciPy | 3.74 ms | 17.4 µs | 215.74x |  |
| `grid/linear/3d/max200/eval10/flipped` | construction | SciPy | 3.64 ms | 17.5 µs | 208.46x |  |
| `grid/linear/3d/max10/eval10/ascending` | construction | interpn | 183 ns | 0 | 208.00x | 2e-16 |
| `grid/linear/3d/max10/eval10000/ascending` | construction | interpn | 171 ns | 0 | 194.56x | 4e-16 |
| `grid/linear/3d/max10/eval1000/ascending` | construction | interpn | 183 ns | 0 | 187.60x | 5e-16 |
| `grid/linear/2d/max10/eval10/ascending` | construction | interpn | 128 ns | 0 | 187.43x | 1e-16 |
| `grid/linear/3d/max200/eval1000/flipped` | construction | SciPy | 3.23 ms | 17.5 µs | 184.45x |  |
| `grid/linear/3d/max10/eval100/ascending` | construction | interpn | 171 ns | 0 | 175.50x | 3e-16 |
| `grid/linear/2d/max10/eval100/ascending` | construction | interpn | 120 ns | 0 | 175.14x | 2e-16 |
| `grid/linear/2d/max10/eval1000/ascending` | construction | interpn | 119 ns | 0 | 173.57x | 3e-16 |
| `grid/linear/2d/max10/eval10000/ascending` | construction | interpn | 117 ns | 0 | 170.57x | 3e-16 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | SciPy | 1.1 ms | 6.48 µs | 169.74x |  |
| `grid/pchip/3d/each40/eval10000/ascending` | construction | SciPy | 472 µs | 8.76 µs | 53.90x |  |
| `grid/linear/2d/each400/eval10000/ascending` | construction | SciPy | 327 µs | 6.15 µs | 53.21x |  |
| `grid/linear/2d/max200/eval10000/ascending` | construction | ndarray-interp | 7.58 µs | 513 ns | 14.76x | 5e-16 |
| `grid/linear/2d/max200/eval10/ascending` | construction | ndarray-interp | 6.9 µs | 469 ns | 14.71x | 2e-16 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | ndarray-interp | 6.92 µs | 487 ns | 14.23x | 4e-16 |
| `grid/linear/2d/max200/eval100/ascending` | construction | ndarray-interp | 6.83 µs | 512 ns | 13.34x | 3e-16 |
| `oned/linear/n10000/shuffled` | construction | ndarray-interp | 101 µs | 18.2 µs | 5.54x | 0e+00 |
| `oned/linear/n10000/sorted` | construction | ndarray-interp | 93.7 µs | 17.6 µs | 5.34x | 0e+00 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | SciPy | 21.6 µs | 7.63 µs | 2.83x |  |
| `grid/linear/4d/max40/eval1000/flipped` | construction | SciPy | 56.4 µs | 20 µs | 2.82x |  |
| `grid/linear/4d/max40/eval100/flipped` | construction | SciPy | 55.8 µs | 20.1 µs | 2.78x |  |
| `oned/linear/n10000/shuffled` | construction | scirs2-interpolate | 101 µs | 36.9 µs | 2.74x | 0e+00 |
| `grid/linear/4d/max40/eval10000/flipped` | construction | SciPy | 53.5 µs | 20.1 µs | 2.67x |  |
| `grid/linear/4d/max40/eval10/flipped` | construction | SciPy | 56.8 µs | 21.7 µs | 2.62x |  |
| `grid/linear/4d/max40/eval100/ascending` | construction | scirs2-interpolate | 14.2 µs | 5.64 µs | 2.52x | 3e-16 |
| `grid/linear/2d/max40/eval10000/ascending` | construction | ndarray-interp | 370 ns | 151 ns | 2.46x | 4e-16 |
| `grid/linear/2d/max40/eval10/ascending` | construction | ndarray-interp | 367 ns | 150 ns | 2.45x | 2e-16 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | scirs2-interpolate | 21.6 µs | 8.86 µs | 2.44x | 0e+00 |
| `grid/linear/4d/max40/eval1000/ascending` | construction | scirs2-interpolate | 14.3 µs | 5.91 µs | 2.41x | 2e-16 |
| `grid/linear/4d/max40/eval10/ascending` | construction | scirs2-interpolate | 13.9 µs | 5.87 µs | 2.37x | 2e-16 |
| `oned/linear/n10000/sorted` | construction | scirs2-interpolate | 93.7 µs | 39.6 µs | 2.36x | 0e+00 |
| `grid/linear/2d/max40/eval100/ascending` | construction | ndarray-interp | 354 ns | 150 ns | 2.36x | 2e-16 |
| `grid/linear/2d/max200/eval10000/ascending` | construction | scirs2-interpolate | 7.58 µs | 3.25 µs | 2.33x | 2e-16 |
| `oned/linear/n10000/shuffled` | evaluation | ndarray-interp | 642 µs | 278 µs | 2.31x | 0e+00 |
| `grid/linear/2d/max200/eval10/ascending` | construction | scirs2-interpolate | 6.9 µs | 3 µs | 2.30x | 1e-16 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | scirs2-interpolate | 6.92 µs | 3.06 µs | 2.26x | 2e-16 |
| `grid/linear/4d/max40/eval10000/ascending` | construction | scirs2-interpolate | 13.3 µs | 5.95 µs | 2.24x | 4e-16 |
| `grid/linear/3d/max40/eval10/ascending` | construction | scirs2-interpolate | 2.96 µs | 1.33 µs | 2.23x | 1e-16 |
| `grid/linear/2d/max40/eval1000/ascending` | construction | ndarray-interp | 364 ns | 164 ns | 2.23x | 4e-16 |
| `grid/linear/3d/max40/eval1000/ascending` | construction | scirs2-interpolate | 3.13 µs | 1.41 µs | 2.21x | 2e-16 |
| `grid/linear/3d/max40/eval100/ascending` | construction | scirs2-interpolate | 3.05 µs | 1.42 µs | 2.15x | 1e-16 |
| `grid/linear/2d/max200/eval100/ascending` | construction | scirs2-interpolate | 6.83 µs | 3.28 µs | 2.08x | 2e-16 |
| `grid/linear/3d/max40/eval10000/ascending` | construction | scirs2-interpolate | 2.76 µs | 1.34 µs | 2.05x | 2e-16 |
| `oned/linear/n10/shuffled` | construction | ndarray-interp | 89.5 ns | 45.6 ns | 1.96x | 0e+00 |
| `oned/akima/n1000/shuffled` | construction | scirs2-interpolate | 9.66 µs | 4.94 µs | 1.96x | 3e-17 |
| `oned/pchip/n1000/sorted` | construction | scirs2-interpolate | 8.85 µs | 4.61 µs | 1.92x | 2e-16 |
| `oned/pchip/n1000/shuffled` | construction | scirs2-interpolate | 8.79 µs | 4.62 µs | 1.90x | 2e-16 |
| `oned/linear/n10000/sorted` | construction | SciPy | 93.7 µs | 49.3 µs | 1.90x |  |
| `grid/nearest/2d/each40/eval10000/ascending` | construction | scirs2-interpolate | 561 ns | 296 ns | 1.89x | 0e+00 |
| `oned/linear/n10/sorted` | construction | ndarray-interp | 83.3 ns | 45 ns | 1.85x | 0e+00 |
| `oned/cubic_spline/n10000/shuffled` | evaluation | ndarray-interp | 698 µs | 379 µs | 1.84x | 2e-16 |
| `grid/linear/2d/max40/eval10/ascending` | construction | scirs2-interpolate | 367 ns | 208 ns | 1.76x | 4e-17 |
| `oned/akima/n1000/sorted` | construction | scirs2-interpolate | 9.92 µs | 5.67 µs | 1.75x | 3e-17 |
| `grid/linear/2d/max40/eval10000/ascending` | construction | scirs2-interpolate | 370 ns | 212 ns | 1.74x | 2e-16 |
| `oned/akima/n100/shuffled` | construction | scirs2-interpolate | 991 ns | 579 ns | 1.71x | 1e-16 |
| `grid/linear/2d/max40/eval100/ascending` | construction | scirs2-interpolate | 354 ns | 207 ns | 1.71x | 3e-16 |
| `grid/linear/4d/max40/eval1000/ascending` | evaluation | interpn | 45.9 µs | 27.6 µs | 1.66x | 6e-16 |
| `oned/akima/n100/sorted` | construction | scirs2-interpolate | 995 ns | 600 ns | 1.66x | 1e-16 |
| `grid/linear/2d/max40/eval1000/ascending` | construction | scirs2-interpolate | 364 ns | 226 ns | 1.61x | 3e-16 |
| `grid/linear/4d/max40/eval10/ascending` | evaluation | interpn | 466 ns | 295 ns | 1.58x | 2e-16 |
| `oned/linear/n1000/shuffled` | construction | scirs2-interpolate | 1.68 µs | 1.07 µs | 1.57x | 0e+00 |
| `grid/linear/4d/max40/eval10000/ascending` | evaluation | interpn | 420 µs | 269 µs | 1.56x | 5e-16 |
| `oned/akima/n50/sorted` | construction | scirs2-interpolate | 512 ns | 329 ns | 1.55x | 1e-16 |
| `grid/linear/2d/max10/eval100/ascending` | construction | ndarray-interp | 120 ns | 78.3 ns | 1.53x | 4e-16 |
| `oned/pchip/n100/sorted` | construction | scirs2-interpolate | 920 ns | 603 ns | 1.53x | 2e-16 |
| `oned/akima/n50/shuffled` | construction | scirs2-interpolate | 509 ns | 334 ns | 1.52x | 1e-16 |
| `grid/linear/2d/max10/eval1000/ascending` | construction | ndarray-interp | 119 ns | 78.3 ns | 1.51x | 3e-16 |
| `oned/linear/n1000/shuffled` | evaluation | ndarray-interp | 40.7 µs | 27 µs | 1.51x | 0e+00 |
| `oned/pchip/n100/shuffled` | construction | scirs2-interpolate | 923 ns | 613 ns | 1.51x | 2e-16 |
| `grid/linear/2d/max10/eval10/ascending` | construction | ndarray-interp | 128 ns | 85.6 ns | 1.50x | 1e-16 |
| `grid/linear/2d/max200/eval10000/flipped` | construction | SciPy | 18.1 µs | 12.2 µs | 1.48x |  |
| `oned/linear/n1000/sorted` | construction | scirs2-interpolate | 1.7 µs | 1.15 µs | 1.48x | 0e+00 |
| `grid/linear/4d/max40/eval100/ascending` | evaluation | interpn | 4.2 µs | 2.9 µs | 1.45x | 2e-16 |
| `grid/linear/4d/max200/eval10000/ascending` | evaluation | interpn | 537 µs | 371 µs | 1.45x | 4e-16 |
| `grid/linear/3d/max200/eval100/ascending` | construction | scirs2-interpolate | 1.96 ms | 1.39 ms | 1.41x | 3e-16 |
| `oned/pchip/n50/shuffled` | construction | scirs2-interpolate | 473 ns | 338 ns | 1.40x | 1e-16 |
| `oned/linear/n50/sorted` | construction | ndarray-interp | 148 ns | 105 ns | 1.40x | 0e+00 |
| `grid/linear/4d/max40/eval10000/ascending` | construction | SciPy | 13.3 µs | 9.57 µs | 1.39x |  |
| `grid/linear/4d/max200/eval10/ascending` | evaluation | interpn | 543 ns | 390 ns | 1.39x | 2e-16 |
| `grid/linear/4d/max40/eval1000/ascending` | construction | SciPy | 14.3 µs | 10.2 µs | 1.39x |  |
| `grid/linear/4d/max40/eval100/ascending` | construction | SciPy | 14.2 µs | 10.2 µs | 1.39x |  |
| `grid/linear/3d/max200/eval10/ascending` | construction | scirs2-interpolate | 1.9 ms | 1.38 ms | 1.38x | 3e-16 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | scirs2-interpolate | 1.1 ms | 802 µs | 1.37x | 2e-16 |
| `grid/linear/3d/max10/eval1000/ascending` | evaluation | interpn | 18.8 µs | 13.7 µs | 1.37x | 5e-16 |
| `oned/linear/n100/sorted` | construction | ndarray-interp | 234 ns | 171 ns | 1.37x | 0e+00 |
| `oned/linear/n50/shuffled` | construction | ndarray-interp | 148 ns | 109 ns | 1.37x | 0e+00 |
| `grid/linear/4d/max40/eval10/ascending` | construction | SciPy | 13.9 µs | 10.2 µs | 1.36x |  |
| `grid/linear/2d/max10/eval10000/ascending` | construction | ndarray-interp | 117 ns | 85.7 ns | 1.36x | 4e-16 |
| `oned/pchip/n50/sorted` | construction | scirs2-interpolate | 472 ns | 351 ns | 1.35x | 1e-16 |
| `grid/linear/3d/max200/eval1000/ascending` | construction | scirs2-interpolate | 1.85 ms | 1.39 ms | 1.33x | 3e-16 |
| `grid/linear/2d/max10/eval10000/ascending` | evaluation | interpn | 133 µs | 100 µs | 1.33x | 3e-16 |
| `grid/linear/2d/max10/eval10/ascending` | evaluation | interpn | 135 ns | 104 ns | 1.30x | 1e-16 |
| `grid/linear/4d/max200/eval1000/ascending` | construction | scirs2-interpolate | 49.9 ms | 38.5 ms | 1.30x | 3e-16 |
| `grid/linear/3d/max200/eval10000/ascending` | construction | scirs2-interpolate | 1.85 ms | 1.43 ms | 1.30x | 3e-16 |
| `grid/linear/2d/max10/eval100/ascending` | evaluation | interpn | 1.22 µs | 943 ns | 1.29x | 2e-16 |
| `grid/linear/2d/max10/eval1000/ascending` | evaluation | interpn | 12.2 µs | 9.59 µs | 1.28x | 3e-16 |
| `grid/linear/4d/max200/eval100/ascending` | evaluation | interpn | 4.91 µs | 3.93 µs | 1.25x | 3e-16 |
| `oned/linear/n50/shuffled` | construction | scirs2-interpolate | 148 ns | 119 ns | 1.25x | 0e+00 |
| `grid/linear/4d/max200/eval10/ascending` | construction | scirs2-interpolate | 47.4 ms | 38.3 ms | 1.24x | 2e-16 |
| `grid/linear/4d/max200/eval1000/ascending` | evaluation | interpn | 49.3 µs | 39.8 µs | 1.24x | 3e-16 |
| `oned/cubic_spline/n1000/shuffled` | evaluation | ndarray-interp | 44.1 µs | 35.8 µs | 1.23x | 2e-16 |
| `oned/linear/n1000/sorted` | construction | ndarray-interp | 1.7 µs | 1.38 µs | 1.23x | 0e+00 |
| `grid/linear/3d/max40/eval10/ascending` | evaluation | interpn | 247 ns | 201 ns | 1.23x | 2e-16 |
| `grid/linear/4d/max200/eval100/ascending` | construction | scirs2-interpolate | 48.7 ms | 39.9 ms | 1.22x | 1e-16 |
| `grid/linear/4d/max200/eval10000/ascending` | construction | scirs2-interpolate | 48.6 ms | 40.1 ms | 1.21x | 4e-16 |
| `grid/linear/2d/max40/eval10000/ascending` | evaluation | interpn | 166 µs | 137 µs | 1.21x | 4e-16 |
| `grid/linear/2d/max200/eval10/flipped` | construction | SciPy | 14.5 µs | 12 µs | 1.21x |  |
| `grid/linear/3d/max10/eval10/ascending` | evaluation | interpn | 177 ns | 147 ns | 1.20x | 2e-16 |
| `grid/linear/3d/max40/eval100/ascending` | evaluation | interpn | 2.42 µs | 2.04 µs | 1.19x | 2e-16 |
| `grid/linear/2d/max200/eval10000/ascending` | construction | SciPy | 7.58 µs | 6.39 µs | 1.19x |  |
| `oned/linear/n100/shuffled` | construction | scirs2-interpolate | 217 ns | 183 ns | 1.18x | 0e+00 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | scirs2-interpolate | 327 µs | 276 µs | 1.18x | 2e-16 |
| `oned/linear/n100/sorted` | construction | scirs2-interpolate | 234 ns | 198 ns | 1.18x | 0e+00 |
| `grid/linear/2d/max40/eval100/ascending` | evaluation | interpn | 1.55 µs | 1.31 µs | 1.18x | 2e-16 |
| `grid/linear/2d/max200/eval100/ascending` | construction | SciPy | 6.83 µs | 5.81 µs | 1.18x |  |
| `grid/linear/3d/max10/eval10000/ascending` | evaluation | interpn | 175 µs | 149 µs | 1.18x | 4e-16 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | SciPy | 6.92 µs | 5.9 µs | 1.17x |  |
| `oned/linear/n100/shuffled` | construction | ndarray-interp | 217 ns | 187 ns | 1.16x | 0e+00 |
| `oned/linear/n50/sorted` | construction | scirs2-interpolate | 148 ns | 127 ns | 1.16x | 0e+00 |
| `grid/linear/2d/max40/eval1000/ascending` | evaluation | interpn | 16.6 µs | 14.4 µs | 1.15x | 3e-16 |
| `grid/linear/3d/max10/eval100/ascending` | evaluation | interpn | 1.72 µs | 1.49 µs | 1.15x | 3e-16 |
| `grid/linear/2d/max200/eval1000/flipped` | construction | SciPy | 14.8 µs | 13 µs | 1.14x |  |
| `rbf/quintic/global/n10/scattered` | evaluation | SciPy | 75 µs | 66.8 µs | 1.12x |  |
| `grid/linear/3d/max40/eval10000/ascending` | evaluation | interpn | 243 µs | 216 µs | 1.12x | 3e-16 |
| `oned/linear/n10/shuffled` | construction | scirs2-interpolate | 89.5 ns | 79.9 ns | 1.12x | 0e+00 |
| `oned/linear/n1000/shuffled` | construction | ndarray-interp | 1.68 µs | 1.52 µs | 1.11x | 0e+00 |
| `oned/cubic_spline/n1000/sorted` | construction | scirs2-interpolate | 19 µs | 17.2 µs | 1.10x | 4e-04 |
| `grid/nearest/2d/each40/eval10000/ascending` | evaluation | interpn | 153 µs | 139 µs | 1.10x | 0e+00 |

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
| `linear/n10/sorted` | 121 ns | 4.98 µs | 257 ns | 254 ns | 2.1x faster |
| `linear/n10/shuffled` | 182 ns | 5.49 µs | 262 ns | 251 ns | 1.4x faster |
| `nearest/n10/sorted` | 144 ns | 9.39 µs | 250 ns |  | 1.7x faster |
| `nearest/n10/shuffled` | 204 ns | 9.63 µs | 278 ns |  | 1.4x faster |
| `previous/n10/sorted` | 155 ns | 10.7 µs |  |  | 69x faster |
| `previous/n10/shuffled` | 214 ns | 11.7 µs |  |  | 54x faster |
| `bspline1/n10/sorted` | 332 ns | 2.93 µs |  |  | 8.8x faster |
| `bspline1/n10/shuffled` | 364 ns | 2.74 µs |  |  | 7.5x faster |
| `bspline2/n10/sorted` | 407 ns | 3.22 µs |  |  | 7.9x faster |
| `bspline2/n10/shuffled` | 472 ns | 3.12 µs |  |  | 6.6x faster |
| `bspline3/n10/sorted` | 589 ns | 3.65 µs |  |  | 6.2x faster |
| `bspline3/n10/shuffled` | 609 ns | 3.86 µs |  |  | 6.3x faster |
| `cubic_spline/n10/sorted` | 127 ns | 3.7 µs | 277 ns | 338 ns | 2.2x faster |
| `cubic_spline/n10/shuffled` | 212 ns | 3.73 µs | 279 ns | 365 ns | 1.3x faster |
| `pchip/n10/sorted` | 127 ns | 3.86 µs | 299 ns |  | 2.4x faster |
| `pchip/n10/shuffled` | 214 ns | 3.69 µs | 297 ns |  | 1.4x faster |
| `akima/n10/sorted` | 129 ns | 3.8 µs | 251 ns |  | 1.9x faster |
| `akima/n10/shuffled` | 213 ns | 3.7 µs | 269 ns |  | 1.3x faster |
| `makima/n10/sorted` | 127 ns | 3.87 µs |  |  | 31x faster |
| `makima/n10/shuffled` | 231 ns | 3.76 µs |  |  | 16x faster |
| `linear/n50/sorted` | 647 ns | 6.63 µs | 1.54 µs | 1.35 µs | 2.1x faster |
| `linear/n50/shuffled` | 1.29 µs | 8.23 µs | 1.55 µs | 1.27 µs | **1x slower** |
| `nearest/n50/sorted` | 721 ns | 11.3 µs | 1.57 µs |  | 2.2x faster |
| `nearest/n50/shuffled` | 1.45 µs | 11.5 µs | 1.57 µs |  | 1.1x faster |
| `previous/n50/sorted` | 747 ns | 13 µs |  |  | 17x faster |
| `previous/n50/shuffled` | 1.31 µs | 13 µs |  |  | 9.9x faster |
| `bspline1/n50/sorted` | 1.46 µs | 4.17 µs |  |  | 2.9x faster |
| `bspline1/n50/shuffled` | 1.96 µs | 5.29 µs |  |  | 2.7x faster |
| `bspline2/n50/sorted` | 1.97 µs | 6.15 µs |  |  | 3.1x faster |
| `bspline2/n50/shuffled` | 2.78 µs | 6.72 µs |  |  | 2.4x faster |
| `bspline3/n50/sorted` | 2.64 µs | 8.7 µs |  |  | 3.3x faster |
| `bspline3/n50/shuffled` | 3.67 µs | 9.32 µs |  |  | 2.5x faster |
| `cubic_spline/n50/sorted` | 707 ns | 5.14 µs | 4.39 µs | 1.64 µs | 2.3x faster |
| `cubic_spline/n50/shuffled` | 1.34 µs | 5.84 µs | 4.07 µs | 1.65 µs | 1.2x faster |
| `pchip/n50/sorted` | 704 ns | 5.2 µs | 3.93 µs |  | 5.6x faster |
| `pchip/n50/shuffled` | 1.35 µs | 5.5 µs | 3.66 µs |  | 2.7x faster |
| `akima/n50/sorted` | 714 ns | 5.56 µs | 4.41 µs |  | 6.2x faster |
| `akima/n50/shuffled` | 1.41 µs | 5.54 µs | 3.92 µs |  | 2.8x faster |
| `makima/n50/sorted` | 721 ns | 5.31 µs |  |  | 7.4x faster |
| `makima/n50/shuffled` | 1.46 µs | 5.94 µs |  |  | 4.1x faster |
| `linear/n100/sorted` | 1.42 µs | 7.93 µs | 3.6 µs | 2.53 µs | 1.8x faster |
| `linear/n100/shuffled` | 2.69 µs | 13.9 µs | 3.61 µs | 2.61 µs | **1x slower** |
| `nearest/n100/sorted` | 1.57 µs | 14.6 µs | 3.33 µs |  | 2.1x faster |
| `nearest/n100/shuffled` | 2.94 µs | 14.3 µs | 3.27 µs |  | 1.1x faster |
| `previous/n100/sorted` | 1.5 µs | 15.1 µs |  |  | 10x faster |
| `previous/n100/shuffled` | 2.87 µs | 15.4 µs |  |  | 5.4x faster |
| `bspline1/n100/sorted` | 3.19 µs | 5.32 µs |  |  | 1.7x faster |
| `bspline1/n100/shuffled` | 4.32 µs | 12.8 µs |  |  | 3x faster |
| `bspline2/n100/sorted` | 3.97 µs | 10.9 µs |  |  | 2.7x faster |
| `bspline2/n100/shuffled` | 5.91 µs | 13.2 µs |  |  | 2.2x faster |
| `bspline3/n100/sorted` | 5.73 µs | 16.4 µs |  |  | 2.9x faster |
| `bspline3/n100/shuffled` | 7.82 µs | 19.6 µs |  |  | 2.5x faster |
| `cubic_spline/n100/sorted` | 1.52 µs | 7.4 µs | 18.7 µs | 3.31 µs | 2.2x faster |
| `cubic_spline/n100/shuffled` | 3.38 µs | 8.52 µs | 17.9 µs | 3.34 µs | **1x slower** |
| `pchip/n100/sorted` | 1.51 µs | 7.01 µs | 17.2 µs |  | 4.7x faster |
| `pchip/n100/shuffled` | 3.35 µs | 8.53 µs | 17.9 µs |  | 2.5x faster |
| `akima/n100/sorted` | 1.53 µs | 7.49 µs | 18.2 µs |  | 4.9x faster |
| `akima/n100/shuffled` | 3.33 µs | 8.1 µs | 19.2 µs |  | 2.4x faster |
| `makima/n100/sorted` | 1.4 µs | 7.16 µs |  |  | 5.1x faster |
| `makima/n100/shuffled` | 3.06 µs | 8.73 µs |  |  | 2.8x faster |
| `linear/n1000/sorted` | 13.1 µs | 29.2 µs | 111 µs | 26.8 µs | 2x faster |
| `linear/n1000/shuffled` | 40.7 µs | 156 µs | 211 µs | 27 µs | **1.5x slower** |
| `nearest/n1000/sorted` | 14.5 µs | 58 µs | 154 µs |  | 4x faster |
| `nearest/n1000/shuffled` | 41.7 µs | 203 µs | 248 µs |  | 4.9x faster |
| `previous/n1000/sorted` | 14.3 µs | 64.2 µs |  |  | 4.5x faster |
| `previous/n1000/shuffled` | 40.9 µs | 206 µs |  |  | 5x faster |
| `bspline1/n1000/sorted` | 31.3 µs | 30.9 µs |  |  | **1x slower** |
| `bspline1/n1000/shuffled` | 60.4 µs | 579 µs |  |  | 9.6x faster |
| `bspline2/n1000/sorted` | 42.2 µs | 79.4 µs |  |  | 1.9x faster |
| `bspline2/n1000/shuffled` | 78.8 µs | 591 µs |  |  | 7.5x faster |
| `bspline3/n1000/sorted` | 57.9 µs | 140 µs |  |  | 2.4x faster |
| `bspline3/n1000/shuffled` | 95.2 µs | 630 µs |  |  | 6.6x faster |
| `cubic_spline/n1000/sorted` | 15.3 µs | 38.6 µs | 1.46 ms | 32.8 µs | 2.1x faster |
| `cubic_spline/n1000/shuffled` | 44.1 µs | 225 µs | 1.45 ms | 35.8 µs | **1.2x slower** |
| `pchip/n1000/sorted` | 15.3 µs | 37.7 µs | 1.21 ms |  | 2.5x faster |
| `pchip/n1000/shuffled` | 44.1 µs | 220 µs | 1.22 ms |  | 5x faster |
| `akima/n1000/sorted` | 15.4 µs | 38 µs | 1.55 ms |  | 2.5x faster |
| `akima/n1000/shuffled` | 48.6 µs | 239 µs | 1.55 ms |  | 4.9x faster |
| `makima/n1000/sorted` | 15.1 µs | 40.4 µs |  |  | 2.7x faster |
| `makima/n1000/shuffled` | 48.5 µs | 238 µs |  |  | 4.9x faster |
| `linear/n10000/sorted` | 147 µs | 237 µs | 1.66 ms | 251 µs | 1.6x faster |
| `linear/n10000/shuffled` | 642 µs | 2.53 ms | 3.08 ms | 278 µs | **2.3x slower** |
| `nearest/n10000/sorted` | 158 µs | 625 µs | 2.09 ms |  | 3.9x faster |
| `nearest/n10000/shuffled` | 726 µs | 3.21 ms | 3.28 ms |  | 4.4x faster |
| `previous/n10000/sorted` | 166 µs | 657 µs |  |  | 4x faster |
| `previous/n10000/shuffled` | 717 µs | 3.25 ms |  |  | 4.5x faster |
| `bspline1/n10000/sorted` | 330 µs | 307 µs |  |  | **1.1x slower** |
| `bspline1/n10000/shuffled` | 956 µs | 51.3 ms |  |  | 54x faster |
| `bspline2/n10000/sorted` | 451 µs | 825 µs |  |  | 1.8x faster |
| `bspline2/n10000/shuffled` | 1.05 ms | 51 ms |  |  | 48x faster |
| `bspline3/n10000/sorted` | 591 µs | 1.37 ms |  |  | 2.3x faster |
| `bspline3/n10000/shuffled` | 1.39 ms | 51.8 ms |  |  | 37x faster |
| `cubic_spline/n10000/sorted` | 161 µs | 387 µs | 153 ms | 335 µs | 2.1x faster |
| `cubic_spline/n10000/shuffled` | 698 µs | 3.63 ms | 154 ms | 379 µs | **1.8x slower** |
| `pchip/n10000/sorted` | 162 µs | 403 µs | 137 ms |  | 2.5x faster |
| `pchip/n10000/shuffled` | 755 µs | 3.67 ms | 134 ms |  | 4.9x faster |
| `akima/n10000/sorted` | 161 µs | 381 µs | 162 ms |  | 2.4x faster |
| `akima/n10000/shuffled` | 756 µs | 3.76 ms | 153 ms |  | 5x faster |
| `makima/n10000/sorted` | 158 µs | 396 µs |  |  | 2.5x faster |
| `makima/n10000/shuffled` | 756 µs | 3.81 ms |  |  | 5x faster |

### Construction

| case | interpx | SciPy | scirs2-interpolate | ndarray-interp | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: |
| `linear/n10/sorted` | 83.3 ns | 8.66 µs | 86.9 ns | 45 ns | **1.8x slower** |
| `linear/n10/shuffled` | 89.5 ns | 9.52 µs | 79.9 ns | 45.6 ns | **2x slower** |
| `nearest/n10/sorted` | 67.1 ns | 9.52 µs | 86.3 ns |  | 1.3x faster |
| `nearest/n10/shuffled` | 62.5 ns | 9.61 µs | 87.5 ns |  | 1.4x faster |
| `previous/n10/sorted` | 66.4 ns | 9.36 µs |  |  | 1.4e+02x faster |
| `previous/n10/shuffled` | 61.5 ns | 9.35 µs |  |  | 1.5e+02x faster |
| `bspline1/n10/sorted` | 107 ns | 16.1 µs |  |  | 1.5e+02x faster |
| `bspline1/n10/shuffled` | 100 ns | 17.1 µs |  |  | 1.7e+02x faster |
| `bspline2/n10/sorted` | 476 ns | 39.9 µs |  |  | 84x faster |
| `bspline2/n10/shuffled` | 450 ns | 39.5 µs |  |  | 88x faster |
| `bspline3/n10/sorted` | 500 ns | 38.3 µs |  |  | 77x faster |
| `bspline3/n10/shuffled` | 494 ns | 35.8 µs |  |  | 72x faster |
| `cubic_spline/n10/sorted` | 223 ns | 87.2 µs | 417 ns | 497 ns | 1.9x faster |
| `cubic_spline/n10/shuffled` | 224 ns | 89.7 µs | 385 ns | 494 ns | 1.7x faster |
| `pchip/n10/sorted` | 168 ns | 101 µs | 194 ns |  | 1.2x faster |
| `pchip/n10/shuffled` | 183 ns | 101 µs | 195 ns |  | 1.1x faster |
| `akima/n10/sorted` | 176 ns | 76 µs | 181 ns |  | 1x faster |
| `akima/n10/shuffled` | 175 ns | 76 µs | 188 ns |  | 1.1x faster |
| `makima/n10/sorted` | 181 ns | 81.5 µs |  |  | 4.5e+02x faster |
| `makima/n10/shuffled` | 189 ns | 79.5 µs |  |  | 4.2e+02x faster |
| `linear/n50/sorted` | 148 ns | 9.49 µs | 127 ns | 105 ns | **1.4x slower** |
| `linear/n50/shuffled` | 148 ns | 9.66 µs | 119 ns | 109 ns | **1.4x slower** |
| `nearest/n50/sorted` | 84.1 ns | 10.4 µs | 116 ns |  | 1.4x faster |
| `nearest/n50/shuffled` | 83.3 ns | 10 µs | 120 ns |  | 1.4x faster |
| `previous/n50/sorted` | 83.5 ns | 10 µs |  |  | 1.2e+02x faster |
| `previous/n50/shuffled` | 90.6 ns | 10.5 µs |  |  | 1.2e+02x faster |
| `bspline1/n50/sorted` | 172 ns | 16.2 µs |  |  | 94x faster |
| `bspline1/n50/shuffled` | 170 ns | 15.9 µs |  |  | 93x faster |
| `bspline2/n50/sorted` | 1.8 µs | 41.8 µs |  |  | 23x faster |
| `bspline2/n50/shuffled` | 1.79 µs | 39.6 µs |  |  | 22x faster |
| `bspline3/n50/sorted` | 2 µs | 39.1 µs |  |  | 20x faster |
| `bspline3/n50/shuffled` | 2 µs | 41.2 µs |  |  | 21x faster |
| `cubic_spline/n50/sorted` | 825 ns | 88.6 µs | 991 ns | 1.07 µs | 1.2x faster |
| `cubic_spline/n50/shuffled` | 817 ns | 87.2 µs | 929 ns | 1.17 µs | 1.1x faster |
| `pchip/n50/sorted` | 472 ns | 100 µs | 351 ns |  | **1.3x slower** |
| `pchip/n50/shuffled` | 473 ns | 104 µs | 338 ns |  | **1.4x slower** |
| `akima/n50/sorted` | 512 ns | 78 µs | 329 ns |  | **1.6x slower** |
| `akima/n50/shuffled` | 509 ns | 78.4 µs | 334 ns |  | **1.5x slower** |
| `makima/n50/sorted` | 535 ns | 79.1 µs |  |  | 1.5e+02x faster |
| `makima/n50/shuffled` | 579 ns | 86.4 µs |  |  | 1.5e+02x faster |
| `linear/n100/sorted` | 234 ns | 9.87 µs | 198 ns | 171 ns | **1.4x slower** |
| `linear/n100/shuffled` | 217 ns | 10 µs | 183 ns | 187 ns | **1.2x slower** |
| `nearest/n100/sorted` | 127 ns | 10.9 µs | 181 ns |  | 1.4x faster |
| `nearest/n100/shuffled` | 127 ns | 10.6 µs | 195 ns |  | 1.5x faster |
| `previous/n100/sorted` | 138 ns | 9.95 µs |  |  | 72x faster |
| `previous/n100/shuffled` | 129 ns | 10.1 µs |  |  | 78x faster |
| `bspline1/n100/sorted` | 290 ns | 15.7 µs |  |  | 54x faster |
| `bspline1/n100/shuffled` | 292 ns | 17.3 µs |  |  | 59x faster |
| `bspline2/n100/sorted` | 3.79 µs | 43.5 µs |  |  | 11x faster |
| `bspline2/n100/shuffled` | 3.77 µs | 46.8 µs |  |  | 12x faster |
| `bspline3/n100/sorted` | 4.22 µs | 47 µs |  |  | 11x faster |
| `bspline3/n100/shuffled` | 3.9 µs | 46.8 µs |  |  | 12x faster |
| `cubic_spline/n100/sorted` | 1.7 µs | 88.7 µs | 1.8 µs | 2 µs | 1.1x faster |
| `cubic_spline/n100/shuffled` | 1.7 µs | 95.3 µs | 1.8 µs | 2.05 µs | 1.1x faster |
| `pchip/n100/sorted` | 920 ns | 111 µs | 603 ns |  | **1.5x slower** |
| `pchip/n100/shuffled` | 923 ns | 112 µs | 613 ns |  | **1.5x slower** |
| `akima/n100/sorted` | 995 ns | 84.7 µs | 600 ns |  | **1.7x slower** |
| `akima/n100/shuffled` | 991 ns | 80.6 µs | 579 ns |  | **1.7x slower** |
| `makima/n100/sorted` | 954 ns | 81.1 µs |  |  | 85x faster |
| `makima/n100/shuffled` | 959 ns | 85 µs |  |  | 89x faster |
| `linear/n1000/sorted` | 1.7 µs | 13.1 µs | 1.15 µs | 1.38 µs | **1.5x slower** |
| `linear/n1000/shuffled` | 1.68 µs | 13.6 µs | 1.07 µs | 1.52 µs | **1.6x slower** |
| `nearest/n1000/sorted` | 828 ns | 15.7 µs | 1.15 µs |  | 1.4x faster |
| `nearest/n1000/shuffled` | 773 ns | 15.2 µs | 1.16 µs |  | 1.5x faster |
| `previous/n1000/sorted` | 747 ns | 17.5 µs |  |  | 23x faster |
| `previous/n1000/shuffled` | 778 ns | 16.6 µs |  |  | 21x faster |
| `bspline1/n1000/sorted` | 1.74 µs | 17.5 µs |  |  | 10x faster |
| `bspline1/n1000/shuffled` | 1.62 µs | 17.8 µs |  |  | 11x faster |
| `bspline2/n1000/sorted` | 36.2 µs | 119 µs |  |  | 3.3x faster |
| `bspline2/n1000/shuffled` | 36.3 µs | 119 µs |  |  | 3.3x faster |
| `bspline3/n1000/sorted` | 42.8 µs | 139 µs |  |  | 3.2x faster |
| `bspline3/n1000/shuffled` | 42.7 µs | 136 µs |  |  | 3.2x faster |
| `cubic_spline/n1000/sorted` | 19 µs | 117 µs | 17.2 µs | 15.8 µs | **1.2x slower** |
| `cubic_spline/n1000/shuffled` | 17.6 µs | 111 µs | 16.2 µs | 16.1 µs | **1.1x slower** |
| `pchip/n1000/sorted` | 8.85 µs | 124 µs | 4.61 µs |  | **1.9x slower** |
| `pchip/n1000/shuffled` | 8.79 µs | 120 µs | 4.62 µs |  | **1.9x slower** |
| `akima/n1000/sorted` | 9.92 µs | 101 µs | 5.67 µs |  | **1.7x slower** |
| `akima/n1000/shuffled` | 9.66 µs | 104 µs | 4.94 µs |  | **2x slower** |
| `makima/n1000/sorted` | 10.5 µs | 106 µs |  |  | 10x faster |
| `makima/n1000/shuffled` | 10.3 µs | 107 µs |  |  | 10x faster |
| `linear/n10000/sorted` | 93.7 µs | 49.3 µs | 39.6 µs | 17.6 µs | **5.3x slower** |
| `linear/n10000/shuffled` | 101 µs | 52 µs | 36.9 µs | 18.2 µs | **5.5x slower** |
| `nearest/n10000/sorted` | 35.5 µs | 57.2 µs | 36.1 µs |  | 1x faster |
| `nearest/n10000/shuffled` | 37.4 µs | 57.9 µs | 42.2 µs |  | 1.1x faster |
| `previous/n10000/sorted` | 52 µs | 79.2 µs |  |  | 1.5x faster |
| `previous/n10000/shuffled` | 37.7 µs | 80.2 µs |  |  | 2.1x faster |
| `bspline1/n10000/sorted` | 69.6 µs | 24.4 µs |  |  | **2.9x slower** |
| `bspline1/n10000/shuffled` | 54.2 µs | 24.8 µs |  |  | **2.2x slower** |
| `bspline2/n10000/sorted` | 522 µs | 813 µs |  |  | 1.6x faster |
| `bspline2/n10000/shuffled` | 527 µs | 839 µs |  |  | 1.6x faster |
| `bspline3/n10000/sorted` | 654 µs | 1.01 ms |  |  | 1.5x faster |
| `bspline3/n10000/shuffled` | 636 µs | 1.02 ms |  |  | 1.6x faster |
| `cubic_spline/n10000/sorted` | 249 µs | 766 µs | 320 µs | 381 µs | 1.3x faster |
| `cubic_spline/n10000/shuffled` | 242 µs | 752 µs | 325 µs | 380 µs | 1.3x faster |
| `pchip/n10000/sorted` | 125 µs | 271 µs | 177 µs |  | 1.4x faster |
| `pchip/n10000/shuffled` | 107 µs | 274 µs | 169 µs |  | 1.6x faster |
| `akima/n10000/sorted` | 112 µs | 713 µs | 147 µs |  | 1.3x faster |
| `akima/n10000/shuffled` | 115 µs | 740 µs | 152 µs |  | 1.3x faster |
| `makima/n10000/sorted` | 130 µs | 788 µs |  |  | 6.1x faster |
| `makima/n10000/shuffled` | 113 µs | 849 µs |  |  | 7.5x faster |

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
| `linear/2d/max10/eval10/ascending` | 135 ns | 135 ns | 21.4 µs | 857 ns | 166 ns | 104 ns | **1.3x slower** |
| `linear/2d/max10/eval10/flipped` | 135 ns | 123 ns | 21.4 µs |  |  |  | 1.6e+02x faster |
| `linear/2d/max10/eval100/ascending` | 1.22 µs | 1.22 µs | 21.4 µs | 8.99 µs | 1.38 µs | 943 ns | **1.3x slower** |
| `linear/2d/max10/eval100/flipped` | 1.22 µs | 1.22 µs | 21.5 µs |  |  |  | 18x faster |
| `linear/2d/max10/eval1000/ascending` | 12.2 µs | 12.2 µs | 40.8 µs | 84.1 µs | 14.3 µs | 9.59 µs | **1.3x slower** |
| `linear/2d/max10/eval1000/flipped` | 12.3 µs | 12.2 µs | 40.2 µs |  |  |  | 3.3x faster |
| `linear/2d/max10/eval10000/ascending` | 133 µs | 133 µs | 358 µs | 960 µs | 239 µs | 100 µs | **1.3x slower** |
| `linear/2d/max10/eval10000/flipped` | 134 µs | 133 µs | 392 µs |  |  |  | 2.9x faster |
| `linear/2d/max40/eval10/ascending` | 153 ns | 153 ns | 19.4 µs | 997 ns | 157 ns | 146 ns | **1.1x slower** |
| `linear/2d/max40/eval10/flipped` | 156 ns | 154 ns | 21.1 µs |  |  |  | 1.4e+02x faster |
| `linear/2d/max40/eval100/ascending` | 1.55 µs | 1.52 µs | 24.1 µs | 10.1 µs | 1.45 µs | 1.31 µs | **1.2x slower** |
| `linear/2d/max40/eval100/flipped` | 1.66 µs | 1.66 µs | 22.1 µs |  |  |  | 13x faster |
| `linear/2d/max40/eval1000/ascending` | 16.6 µs | 16.6 µs | 49.7 µs | 114 µs | 16.4 µs | 14.4 µs | **1.2x slower** |
| `linear/2d/max40/eval1000/flipped` | 15.4 µs | 16.6 µs | 50.2 µs |  |  |  | 3.3x faster |
| `linear/2d/max40/eval10000/ascending` | 166 µs | 166 µs | 614 µs | 1.11 ms | 244 µs | 137 µs | **1.2x slower** |
| `linear/2d/max40/eval10000/flipped` | 153 µs | 157 µs | 618 µs |  |  |  | 4x faster |
| `linear/2d/max200/eval10/ascending` | 192 ns | 191 ns | 21.2 µs | 1.96 µs | 191 ns | 191 ns | **1x slower** |
| `linear/2d/max200/eval10/flipped` | 191 ns | 191 ns | 21.1 µs |  |  |  | 1.1e+02x faster |
| `linear/2d/max200/eval100/ascending` | 1.86 µs | 1.96 µs | 24.1 µs | 21.6 µs | 1.99 µs | 1.87 µs | 1x faster |
| `linear/2d/max200/eval100/flipped` | 1.98 µs | 2.04 µs | 24.1 µs |  |  |  | 12x faster |
| `linear/2d/max200/eval1000/ascending` | 20.4 µs | 20.4 µs | 74.2 µs | 227 µs | 23.2 µs | 18.9 µs | **1.1x slower** |
| `linear/2d/max200/eval1000/flipped` | 18.7 µs | 19.9 µs | 71.8 µs |  |  |  | 3.8x faster |
| `linear/2d/max200/eval10000/ascending` | 188 µs | 204 µs | 879 µs | 2.3 ms | 459 µs | 176 µs | **1.1x slower** |
| `linear/2d/max200/eval10000/flipped` | 187 µs | 187 µs | 872 µs |  |  |  | 4.7x faster |
| `linear/3d/max10/eval10/ascending` | 177 ns | 189 ns | 56.1 µs | 3.16 µs |  | 147 ns | **1.2x slower** |
| `linear/3d/max10/eval10/flipped` | 176 ns | 174 ns | 56.6 µs |  |  |  | 3.2e+02x faster |
| `linear/3d/max10/eval100/ascending` | 1.72 µs | 1.72 µs | 62.7 µs | 33.4 µs |  | 1.49 µs | **1.2x slower** |
| `linear/3d/max10/eval100/flipped` | 1.87 µs | 1.76 µs | 67.3 µs |  |  |  | 36x faster |
| `linear/3d/max10/eval1000/ascending` | 18.8 µs | 17.8 µs | 119 µs | 332 µs |  | 13.7 µs | **1.4x slower** |
| `linear/3d/max10/eval1000/flipped` | 17.9 µs | 17.5 µs | 123 µs |  |  |  | 6.9x faster |
| `linear/3d/max10/eval10000/ascending` | 175 µs | 176 µs | 1.14 ms | 3.21 ms |  | 149 µs | **1.2x slower** |
| `linear/3d/max10/eval10000/flipped` | 174 µs | 173 µs | 1.17 ms |  |  |  | 6.7x faster |
| `linear/3d/max40/eval10/ascending` | 247 ns | 242 ns | 60.8 µs | 3.43 µs |  | 201 ns | **1.2x slower** |
| `linear/3d/max40/eval10/flipped` | 247 ns | 229 ns | 60.7 µs |  |  |  | 2.5e+02x faster |
| `linear/3d/max40/eval100/ascending` | 2.42 µs | 2.41 µs | 64.2 µs | 35.2 µs |  | 2.04 µs | **1.2x slower** |
| `linear/3d/max40/eval100/flipped` | 2.24 µs | 2.25 µs | 68.1 µs |  |  |  | 30x faster |
| `linear/3d/max40/eval1000/ascending` | 22.4 µs | 23 µs | 138 µs | 347 µs |  | 20.4 µs | **1.1x slower** |
| `linear/3d/max40/eval1000/flipped` | 24.4 µs | 24.4 µs | 141 µs |  |  |  | 5.8x faster |
| `linear/3d/max40/eval10000/ascending` | 243 µs | 243 µs | 1.51 ms | 3.69 ms |  | 216 µs | **1.1x slower** |
| `linear/3d/max40/eval10000/flipped` | 224 µs | 226 µs | 1.51 ms |  |  |  | 6.8x faster |
| `linear/3d/max200/eval10/ascending` | 304 ns | 302 ns | 61.2 µs | 4.84 µs |  | 290 ns | **1x slower** |
| `linear/3d/max200/eval10/flipped` | 303 ns | 302 ns | 59.9 µs |  |  |  | 2e+02x faster |
| `linear/3d/max200/eval100/ascending` | 3 µs | 3 µs | 69.3 µs | 50.7 µs |  | 2.87 µs | **1x slower** |
| `linear/3d/max200/eval100/flipped` | 2.89 µs | 3.01 µs | 66.2 µs |  |  |  | 23x faster |
| `linear/3d/max200/eval1000/ascending` | 28.2 µs | 30 µs | 195 µs | 474 µs |  | 28.8 µs | 1x faster |
| `linear/3d/max200/eval1000/flipped` | 29.2 µs | 30.1 µs | 196 µs |  |  |  | 6.7x faster |
| `linear/3d/max200/eval10000/ascending` | 276 µs | 301 µs | 2.02 ms | 5.15 ms |  | 288 µs | 1x faster |
| `linear/3d/max200/eval10000/flipped` | 300 µs | 300 µs | 2.06 ms |  |  |  | 6.9x faster |
| `linear/4d/max10/eval10/ascending` |  |  | 98.5 µs |  |  |  |  |
| `linear/4d/max10/eval10/flipped` |  |  | 106 µs |  |  |  |  |
| `linear/4d/max10/eval100/ascending` |  |  | 111 µs |  |  |  |  |
| `linear/4d/max10/eval100/flipped` |  |  | 115 µs |  |  |  |  |
| `linear/4d/max10/eval1000/ascending` |  |  | 231 µs |  |  |  |  |
| `linear/4d/max10/eval1000/flipped` |  |  | 230 µs |  |  |  |  |
| `linear/4d/max10/eval10000/ascending` |  |  | 2.06 ms |  |  |  |  |
| `linear/4d/max10/eval10000/flipped` |  |  | 2.1 ms |  |  |  |  |
| `linear/4d/max40/eval10/ascending` | 466 ns | 428 ns | 104 µs | 6.03 µs |  | 295 ns | **1.6x slower** |
| `linear/4d/max40/eval10/flipped` | 468 ns | 466 ns | 104 µs |  |  |  | 2.2e+02x faster |
| `linear/4d/max40/eval100/ascending` | 4.2 µs | 4.6 µs | 112 µs | 64 µs |  | 2.9 µs | **1.4x slower** |
| `linear/4d/max40/eval100/flipped` | 4.59 µs | 4.59 µs | 109 µs |  |  |  | 24x faster |
| `linear/4d/max40/eval1000/ascending` | 45.9 µs | 45.9 µs | 253 µs | 629 µs |  | 27.6 µs | **1.7x slower** |
| `linear/4d/max40/eval1000/flipped` | 42.2 µs | 46 µs | 255 µs |  |  |  | 6.1x faster |
| `linear/4d/max40/eval10000/ascending` | 420 µs | 422 µs | 2.33 ms | 6.43 ms |  | 269 µs | **1.6x slower** |
| `linear/4d/max40/eval10000/flipped` | 427 µs | 436 µs | 2.28 ms |  |  |  | 5.3x faster |
| `linear/4d/max200/eval10/ascending` | 543 ns | 518 ns | 98.5 µs | 8.43 µs |  | 390 ns | **1.4x slower** |
| `linear/4d/max200/eval10/flipped` | 496 ns | 499 ns | 102 µs |  |  |  | 2.1e+02x faster |
| `linear/4d/max200/eval100/ascending` | 4.91 µs | 5.37 µs | 118 µs | 73.5 µs |  | 3.93 µs | **1.2x slower** |
| `linear/4d/max200/eval100/flipped` | 4.91 µs | 4.91 µs | 113 µs |  |  |  | 23x faster |
| `linear/4d/max200/eval1000/ascending` | 49.3 µs | 49.7 µs | 324 µs | 745 µs |  | 39.8 µs | **1.2x slower** |
| `linear/4d/max200/eval1000/flipped` | 49.7 µs | 49.1 µs | 312 µs |  |  |  | 6.3x faster |
| `linear/4d/max200/eval10000/ascending` | 537 µs | 493 µs | 3.01 ms | 7.75 ms |  | 371 µs | **1.4x slower** |
| `linear/4d/max200/eval10000/flipped` | 497 µs | 492 µs | 2.96 ms |  |  |  | 6x faster |
| `nearest/2d/each40/eval10000/ascending` | 153 µs | 151 µs | 668 µs | 1.57 ms |  | 139 µs | **1.1x slower** |
| `nearest/3d/each40/eval10000/ascending` | 213 µs | 213 µs | 944 µs | 1.79 ms |  | 208 µs | **1x slower** |
| `pchip/2d/each40/eval10000/ascending` | 537 µs | 535 µs | 1.24 s |  |  |  | 2.3e+03x faster |
| `pchip/3d/each40/eval10000/ascending` | 1.58 ms | 1.58 ms | 3.27 s |  |  |  | 2.1e+03x faster |
| `cubic/2d/each40/eval10000/ascending` | 365 µs | 368 µs | 1.29 ms |  |  | 404 µs | 1.1x faster |
| `cubic/3d/each40/eval10000/ascending` | 813 µs | 811 µs | 3.03 ms |  |  | 841 µs | 1x faster |
| `cubic/2d/each10/eval1000/ascending` | 33.9 µs | 33.3 µs | 137 µs |  |  | 35.4 µs | 1x faster |
| `cubic/2d/each40/eval1000/ascending` | 36.5 µs | 36.3 µs | 143 µs |  |  | 40.2 µs | 1.1x faster |
| `cubic/2d/each100/eval1000/ascending` | 38.1 µs | 41.4 µs | 173 µs |  |  | 43 µs | 1.1x faster |
| `cubic/2d/each200/eval1000/ascending` | 39.6 µs | 42.9 µs | 219 µs |  |  | 45.3 µs | 1.1x faster |
| `cubic/2d/each400/eval1000/ascending` | 45.6 µs | 42.4 µs | 273 µs |  |  | 48.5 µs | 1.1x faster |
| `quintic/2d/each10/eval1000/ascending` | 60.6 µs | 62.3 µs | 229 µs |  |  |  | 3.8x faster |
| `quintic/2d/each40/eval1000/ascending` | 65.4 µs | 65.8 µs | 257 µs |  |  |  | 3.9x faster |
| `linear/2d/each400/eval10000/ascending` | 214 µs | 215 µs | 967 µs | 4.36 ms | 580 µs | 217 µs | 1x faster |
| `linear/2d/each800/eval10000/ascending` | 235 µs | 236 µs | 1.15 ms | 7.46 ms | 580 µs | 243 µs | 1x faster |

### Construction

| case | interpx | interpx, `from_owned` | SciPy | scirs2-interpolate | ndarray-interp | interpn | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/2d/max10/eval10/ascending` | 128 ns | 132 ns | 6.49 µs | 147 ns | 85.6 ns | 0 | **1.9e+02x slower** |
| `linear/2d/max10/eval10/flipped` | 156 ns | 145 ns | 12.8 µs |  |  |  | 82x faster |
| `linear/2d/max10/eval100/ascending` | 120 ns | 113 ns | 6.54 µs | 146 ns | 78.3 ns | 0 | **1.8e+02x slower** |
| `linear/2d/max10/eval100/flipped` | 143 ns | 143 ns | 11.8 µs |  |  |  | 83x faster |
| `linear/2d/max10/eval1000/ascending` | 119 ns | 136 ns | 5.89 µs | 134 ns | 78.3 ns | 0 | **1.7e+02x slower** |
| `linear/2d/max10/eval1000/flipped` | 144 ns | 160 ns | 11.7 µs |  |  |  | 81x faster |
| `linear/2d/max10/eval10000/ascending` | 117 ns | 138 ns | 5.87 µs | 146 ns | 85.7 ns | 0 | **1.7e+02x slower** |
| `linear/2d/max10/eval10000/flipped` | 156 ns | 166 ns | 12.8 µs |  |  |  | 82x faster |
| `linear/2d/max40/eval10/ascending` | 367 ns | 742 ns | 6.15 µs | 208 ns | 150 ns | 0 | **5.4e+02x slower** |
| `linear/2d/max40/eval10/flipped` | 673 ns | 1.15 µs | 11.8 µs |  |  |  | 17x faster |
| `linear/2d/max40/eval100/ascending` | 354 ns | 638 ns | 6.39 µs | 207 ns | 150 ns | 0 | **5.2e+02x slower** |
| `linear/2d/max40/eval100/flipped` | 670 ns | 1.1 µs | 12.5 µs |  |  |  | 19x faster |
| `linear/2d/max40/eval1000/ascending` | 364 ns | 683 ns | 6.09 µs | 226 ns | 164 ns | 0 | **5.3e+02x slower** |
| `linear/2d/max40/eval1000/flipped` | 667 ns | 1.08 µs | 11.7 µs |  |  |  | 18x faster |
| `linear/2d/max40/eval10000/ascending` | 370 ns | 775 ns | 6.13 µs | 212 ns | 151 ns | 0 | **4.7e+02x slower** |
| `linear/2d/max40/eval10000/flipped` | 622 ns | 1.17 µs | 11.8 µs |  |  |  | 19x faster |
| `linear/2d/max200/eval10/ascending` | 6.9 µs | 21 µs | 6.41 µs | 3 µs | 469 ns | 0 | **8.8e+03x slower** |
| `linear/2d/max200/eval10/flipped` | 14.5 µs | 26.4 µs | 12 µs |  |  |  | **1.2x slower** |
| `linear/2d/max200/eval100/ascending` | 6.83 µs | 15.1 µs | 5.81 µs | 3.28 µs | 512 ns | 0 | **1e+04x slower** |
| `linear/2d/max200/eval100/flipped` | 14.3 µs | 24.8 µs | 13.1 µs |  |  |  | **1.1x slower** |
| `linear/2d/max200/eval1000/ascending` | 6.92 µs | 14.6 µs | 5.9 µs | 3.06 µs | 487 ns | 0 | **8.9e+03x slower** |
| `linear/2d/max200/eval1000/flipped` | 14.8 µs | 23.3 µs | 13 µs |  |  |  | **1.1x slower** |
| `linear/2d/max200/eval10000/ascending` | 7.58 µs | 16.5 µs | 6.39 µs | 3.25 µs | 513 ns | 0 | **1.1e+04x slower** |
| `linear/2d/max200/eval10000/flipped` | 18.1 µs | 28 µs | 12.2 µs |  |  |  | **1.5x slower** |
| `linear/3d/max10/eval10/ascending` | 183 ns | 180 ns | 7.69 µs | 169 ns |  | 0 | **2.1e+02x slower** |
| `linear/3d/max10/eval10/flipped` | 292 ns | 309 ns | 15.7 µs |  |  |  | 54x faster |
| `linear/3d/max10/eval100/ascending` | 171 ns | 189 ns | 8.23 µs | 181 ns |  | 0 | **1.8e+02x slower** |
| `linear/3d/max10/eval100/flipped` | 284 ns | 311 ns | 15.7 µs |  |  |  | 55x faster |
| `linear/3d/max10/eval1000/ascending` | 183 ns | 193 ns | 7.76 µs | 180 ns |  | 0 | **1.9e+02x slower** |
| `linear/3d/max10/eval1000/flipped` | 285 ns | 308 ns | 17.1 µs |  |  |  | 60x faster |
| `linear/3d/max10/eval10000/ascending` | 171 ns | 196 ns | 8.31 µs | 182 ns |  | 0 | **1.9e+02x slower** |
| `linear/3d/max10/eval10000/flipped` | 289 ns | 318 ns | 16.9 µs |  |  |  | 59x faster |
| `linear/3d/max40/eval10/ascending` | 2.96 µs | 6.88 µs | 8.31 µs | 1.33 µs |  | 0 | **3.4e+03x slower** |
| `linear/3d/max40/eval10/flipped` | 9.07 µs | 11.6 µs | 16.8 µs |  |  |  | 1.9x faster |
| `linear/3d/max40/eval100/ascending` | 3.05 µs | 6.12 µs | 7.72 µs | 1.42 µs |  | 0 | **3.1e+03x slower** |
| `linear/3d/max40/eval100/flipped` | 9.39 µs | 11.4 µs | 17 µs |  |  |  | 1.8x faster |
| `linear/3d/max40/eval1000/ascending` | 3.13 µs | 5.86 µs | 8.36 µs | 1.41 µs |  | 0 | **3.6e+03x slower** |
| `linear/3d/max40/eval1000/flipped` | 9.69 µs | 12 µs | 16.8 µs |  |  |  | 1.7x faster |
| `linear/3d/max40/eval10000/ascending` | 2.76 µs | 6.38 µs | 8.21 µs | 1.34 µs |  | 0 | **3.1e+03x slower** |
| `linear/3d/max40/eval10000/flipped` | 9.64 µs | 12.7 µs | 17.2 µs |  |  |  | 1.8x faster |
| `linear/3d/max200/eval10/ascending` | 1.9 ms | 725 µs | 7.68 µs | 1.38 ms |  | 0 | **2.2e+06x slower** |
| `linear/3d/max200/eval10/flipped` | 3.64 ms | 2.09 ms | 17.5 µs |  |  |  | **2.1e+02x slower** |
| `linear/3d/max200/eval100/ascending` | 1.96 ms | 686 µs | 7.71 µs | 1.39 ms |  | 0 | **2e+06x slower** |
| `linear/3d/max200/eval100/flipped` | 3.84 ms | 1.99 ms | 16.4 µs |  |  |  | **2.3e+02x slower** |
| `linear/3d/max200/eval1000/ascending` | 1.85 ms | 735 µs | 8.41 µs | 1.39 ms |  | 1.07 ns | **1.7e+06x slower** |
| `linear/3d/max200/eval1000/flipped` | 3.23 ms | 2 ms | 17.5 µs |  |  |  | **1.8e+02x slower** |
| `linear/3d/max200/eval10000/ascending` | 1.85 ms | 685 µs | 8.29 µs | 1.43 ms |  | 0 | **2.1e+06x slower** |
| `linear/3d/max200/eval10000/flipped` | 3.74 ms | 1.97 ms | 17.4 µs |  |  |  | **2.2e+02x slower** |
| `linear/4d/max10/eval10/ascending` |  |  | 10.7 µs |  |  |  |  |
| `linear/4d/max10/eval10/flipped` |  |  | 19.6 µs |  |  |  |  |
| `linear/4d/max10/eval100/ascending` |  |  | 10.2 µs |  |  |  |  |
| `linear/4d/max10/eval100/flipped` |  |  | 18.3 µs |  |  |  |  |
| `linear/4d/max10/eval1000/ascending` |  |  | 10.6 µs |  |  |  |  |
| `linear/4d/max10/eval1000/flipped` |  |  | 18.1 µs |  |  |  |  |
| `linear/4d/max10/eval10000/ascending` |  |  | 10.6 µs |  |  |  |  |
| `linear/4d/max10/eval10000/flipped` |  |  | 18.1 µs |  |  |  |  |
| `linear/4d/max40/eval10/ascending` | 13.9 µs | 41.3 µs | 10.2 µs | 5.87 µs |  | 1.27 ns | **1.1e+04x slower** |
| `linear/4d/max40/eval10/flipped` | 56.8 µs | 66 µs | 21.7 µs |  |  |  | **2.6x slower** |
| `linear/4d/max40/eval100/ascending` | 14.2 µs | 41.6 µs | 10.2 µs | 5.64 µs |  | 1.17 ns | **1.2e+04x slower** |
| `linear/4d/max40/eval100/flipped` | 55.8 µs | 70.1 µs | 20.1 µs |  |  |  | **2.8x slower** |
| `linear/4d/max40/eval1000/ascending` | 14.3 µs | 32.4 µs | 10.2 µs | 5.91 µs |  | 1.17 ns | **1.2e+04x slower** |
| `linear/4d/max40/eval1000/flipped` | 56.4 µs | 66.4 µs | 20 µs |  |  |  | **2.8x slower** |
| `linear/4d/max40/eval10000/ascending` | 13.3 µs | 31.8 µs | 9.57 µs | 5.95 µs |  | 1.17 ns | **1.1e+04x slower** |
| `linear/4d/max40/eval10000/flipped` | 53.5 µs | 74.4 µs | 20.1 µs |  |  |  | **2.7x slower** |
| `linear/4d/max200/eval10/ascending` | 47.4 ms | 15.5 ms | 9.86 µs | 38.3 ms |  | 1.17 ns | **4e+07x slower** |
| `linear/4d/max200/eval10/flipped` | 99.9 ms | 66.4 ms | 20.3 µs |  |  |  | **4.9e+03x slower** |
| `linear/4d/max200/eval100/ascending` | 48.7 ms | 15.5 ms | 10.5 µs | 39.9 ms |  | 1.17 ns | **4.2e+07x slower** |
| `linear/4d/max200/eval100/flipped` | 100 ms | 67.5 ms | 21.3 µs |  |  |  | **4.7e+03x slower** |
| `linear/4d/max200/eval1000/ascending` | 49.9 ms | 15.2 ms | 10.2 µs | 38.5 ms |  | 1.17 ns | **4.3e+07x slower** |
| `linear/4d/max200/eval1000/flipped` | 99.7 ms | 66.9 ms | 22 µs |  |  |  | **4.5e+03x slower** |
| `linear/4d/max200/eval10000/ascending` | 48.6 ms | 16.3 ms | 9.48 µs | 40.1 ms |  | 1.17 ns | **4.1e+07x slower** |
| `linear/4d/max200/eval10000/flipped` | 107 ms | 70.4 ms | 20.3 µs |  |  |  | **5.3e+03x slower** |
| `nearest/2d/each40/eval10000/ascending` | 561 ns | 1.25 µs | 5.83 µs | 296 ns |  | 0 | **8.2e+02x slower** |
| `nearest/3d/each40/eval10000/ascending` | 21.6 µs | 47.5 µs | 7.63 µs | 8.86 µs |  | 0 | **2.5e+04x slower** |
| `pchip/2d/each40/eval10000/ascending` | 6.87 µs | 6.14 µs | 6.78 µs |  |  |  | **1x slower** |
| `pchip/3d/each40/eval10000/ascending` | 472 µs | 486 µs | 8.76 µs |  |  |  | **54x slower** |
| `cubic/2d/each40/eval10000/ascending` | 8.55 µs | 8.14 µs | 2 ms |  |  | 0 | **1.1e+04x slower** |
| `cubic/3d/each40/eval10000/ascending` | 306 µs | 293 µs | 152 ms |  |  | 0 | **3.1e+05x slower** |
| `cubic/2d/each10/eval1000/ascending` | 1.48 µs | 1.52 µs | 1.01 ms |  |  | 0 | **1.7e+03x slower** |
| `cubic/2d/each40/eval1000/ascending` | 8.36 µs | 8.06 µs | 1.93 ms |  |  | 0 | **8.6e+03x slower** |
| `cubic/2d/each100/eval1000/ascending` | 43.3 µs | 65.9 µs | 8.05 ms |  |  | 0 | **4.9e+04x slower** |
| `cubic/2d/each200/eval1000/ascending` | 176 µs | 175 µs | 28.6 ms |  |  | 0 | **1.8e+05x slower** |
| `cubic/2d/each400/eval1000/ascending` | 1.4 ms | 1.17 ms | 145 ms |  |  | 0 | **1.6e+06x slower** |
| `quintic/2d/each10/eval1000/ascending` | 1.95 µs | 1.98 µs | 3.68 ms |  |  |  | 1.9e+03x faster |
| `quintic/2d/each40/eval1000/ascending` | 12.4 µs | 12.1 µs | 7.95 ms |  |  |  | 6.4e+02x faster |
| `linear/2d/each400/eval10000/ascending` | 327 µs | 113 µs | 6.15 µs | 276 µs | 1.15 µs | 0 | **4.8e+05x slower** |
| `linear/2d/each800/eval10000/ascending` | 1.1 ms | 425 µs | 6.48 µs | 802 µs | 2.44 µs | 0 | **1.6e+06x slower** |

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
| `oned/cubic_spline/n10000/sorted` | interpx | 4.4 | +0.8 | +0.8 | +0.8 | 5.9 |
| `oned/cubic_spline/n10000/sorted` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |
| `oned/cubic_spline/n10000/sorted` | scirs2-interpolate | 4.1 | +1.0 | +1.1 | +1.8 | 5.9 |
| `oned/cubic_spline/n10000/sorted` | ndarray-interp | 3.3 | +1.7 | +1.6 | +1.6 | 5.6 |
| `grid/linear/4d/max200/eval10000/ascending` | interpx | 204.8 | +200.2 | +199.9 | +199.9 | 405.0 |
| `grid/linear/4d/max200/eval10000/ascending` | interpx, `from_owned` | 212.0 | +193.0 | +193.0 | +193.0 | 405.0 |
| `grid/linear/4d/max200/eval10000/ascending` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |
| `grid/linear/4d/max200/eval10000/ascending` | scirs2-interpolate | 212.5 | +186.8 | +192.4 | +192.4 | 404.9 |
| `grid/linear/4d/max200/eval10000/ascending` | interpn | 227.1 | +0.0 | +0.0 | +0.0 | 403.1 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpx | 5.6 | +2.8 | +2.7 | +2.7 | 8.5 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpx, `from_owned` | 5.6 | +2.4 | +2.7 | +2.7 | 8.5 |
| `grid/cubic/2d/each400/eval1000/ascending` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpn | 5.4 | +0.1 | +0.1 | +0.1 | 5.5 |
| `grid/linear/2d/each800/eval10000/ascending` | interpx | 10.0 | +4.9 | +4.9 | +4.9 | 15.0 |
| `grid/linear/2d/each800/eval10000/ascending` | interpx, `from_owned` | 13.1 | +0.0 | +1.7 | +1.7 | 14.9 |
| `grid/linear/2d/each800/eval10000/ascending` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |
| `grid/linear/2d/each800/eval10000/ascending` | scirs2-interpolate | 10.7 | +0.0 | +4.1 | +4.2 | 14.9 |
| `grid/linear/2d/each800/eval10000/ascending` | ndarray-interp | 11.2 | +0.0 | +0.0 | +0.0 | 11.2 |
| `grid/linear/2d/each800/eval10000/ascending` | interpn | 9.7 | +0.0 | +0.0 | +0.0 | 9.8 |

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
| One variable | 3e-16 |  | 6e-02 | 2e-16 |  |
| `PPoly` | 9e-16 |  |  |  |  |
| `BSpline` | 4e-16 |  |  |  |  |
| `BSplineCurve` | 4e-16 |  |  |  |  |
| Regular grid | 4e-04 | 4e-04 | 4e-16 | 6e-16 | 7e-01 |
| `NearestNDInterpolator` | 0e+00 |  |  |  |  |
| `RBFInterpolator` | 8e-05 |  | 3e+01 |  |  |
