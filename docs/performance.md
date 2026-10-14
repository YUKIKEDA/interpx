# Performance: interpx against SciPy and other Rust crates

<!-- Written by compare/perf/report.py from a run of compare/perf/run.py. Do not edit by hand: run it again. -->

Measured on 2026-10-14, at commit `056f1175c3`.

## How it was measured

- Machine: Intel64 Family 6 Model 191 Stepping 2, GenuineIntel, 16 logical cores, 48 GB. Windows-11-10.0.26300-SP0.
- rustc 1.99.0 (b940084d7 2026-09-28); Python 3.13.1, NumPy 2.5.4, SciPy 1.18.1.
- The other Rust crates, at the versions in `compare/perf/runners/*/Cargo.toml`: `scirs2-interpolate` (the widest of the Rust crates: one variable, grids, RBF), `ndarray-interp` (the most used for one variable; also bilinear), `interpn` (grids).
- Every library ran each case in a process of its own, pinned to one core, on one thread, on the same data. The libraries took turns case by case, 3 rounds.
- `interpx, from_owned` is interpx built from data that is handed over (`from_owned`) and so not copied. Its column is there for the types that keep their data as it is given (the grid, `NearestNDInterpolator`, `RBFInterpolator`), next to the construction from borrowed data (`new`), which is the one the comparison counts. Each of its samples is given data made just before it, which is not in the nearest cache as the data of the other constructions is: on small grids it reads slower than `new` for that reason.
- `interpx, parallel` is interpx evaluated with `eval_batch_into_par` (the `rayon` feature) on 8 threads, not pinned to a core. Every other column is one thread on one core, and the comparison counts those. A parallel batch times its first points and is shared out only when the rest is worth it, so on small cases the column reads as the serial one does.
- A time is that of one call: construction, or the evaluation of all the queries of the case (`eval_batch_into` for interpx, the call with an array for the others). It is the median of the repeats of a process, and then the median of the rounds.
- A call of two seconds or more is timed once, not repeated; and a library whose construction or evaluation of a case takes that long is run on that case in the first round only. Its time is then one measurement. (The rule goes by the time, whatever the library.)
- interpx is called slower where its median is more than 1.10 times another library's and its fastest repeat is slower than that library's fastest.
- Memory is the resident set size of the process, sampled every millisecond from outside it.
- Command: `perf/run.py --rounds 3`.

## Where interpx is slower

| family | operation | cases | interpx fastest | slower than SciPy | slower than a Rust crate |
| :--- | :--- | ---: | ---: | ---: | ---: |
| One variable | evaluation | 112 | 111 | 0 | 1 |
| One variable | construction | 112 | 77 | 8 | 32 |
| `PPoly` | evaluation | 4 | 4 | 0 | 0 |
| `PPoly` | construction | 4 | 4 | 0 | 0 |
| `BSpline` | evaluation | 4 | 4 | 0 | 0 |
| `BSpline` | construction | 4 | 4 | 0 | 0 |
| `BSplineCurve` | evaluation | 2 | 2 | 0 | 0 |
| `BSplineCurve` | construction | 2 | 2 | 0 | 0 |
| Regular grid | evaluation | 92 | 77 | 0 | 15 |
| Regular grid | construction | 92 | 28 | 39 | 46 |
| `NearestNDInterpolator` | evaluation | 16 | 16 | 0 | 0 |
| `NearestNDInterpolator` | construction | 16 | 16 | 0 | 0 |
| `RBFInterpolator` | evaluation | 61 | 61 | 0 | 0 |
| `RBFInterpolator` | construction | 61 | 61 | 0 | 0 |

197 comparisons (a case, an operation, and a library), by method. "Same values" says whether that library's values are SciPy's up to rounding: where they are not, it computes another interpolant under the same name, and its time is not that of the same work.

| method | operation | against | same values | cases | slower by | the largest |
| :--- | :--- | :--- | :--- | ---: | ---: | :--- |
| `grid/cubic` | construction | interpn | no | 8 | 1.68e+03x to 9.04e+07x | `grid/cubic/2d/auckland-ba32-0401/heldout/ascending`: 79.5 ms against 0 |
| `grid/linear` | construction | interpn | yes | 35 | 179x to 4.46e+07x | `grid/linear/4d/max200/eval1000/ascending`: 47.9 ms against 1.07 ns |
| `grid/nearest` | construction | interpn | yes | 3 | 926x to 2.24e+07x | `grid/nearest/2d/auckland-ba32-0401/heldout/ascending`: 15.3 ms against 0 |
| `grid/linear` | construction | SciPy | yes | 36 | 1.14x to 5.17e+03x | `grid/linear/4d/max200/eval1000/ascending`: 47.9 ms against 9.27 µs |
| `grid/nearest` | construction | SciPy | yes | 2 | 2.78x to 2.18e+03x | `grid/nearest/2d/auckland-ba32-0401/heldout/ascending`: 15.3 ms against 7.04 µs |
| `grid/linear` | construction | ndarray-interp | yes | 15 | 1.56x to 1.75e+03x | `grid/linear/2d/auckland-ba32-0401/heldout/ascending`: 15.9 ms against 9.05 µs |
| `grid/pchip` | construction | SciPy | yes | 1 | 60.6x | `grid/pchip/3d/each40/eval10000/ascending`: 508 µs against 8.38 µs |
| `oned/linear` | construction | ndarray-interp | yes | 14 | 1.28x to 6.1x | `oned/linear/n10000/shuffled-uneven`: 88.8 µs against 14.6 µs |
| `grid/linear` | construction | scirs2-interpolate | yes | 25 | 1.1x to 2.79x | `grid/linear/3d/max40/eval1000/ascending`: 16.6 µs against 5.96 µs |
| `oned/pchip` | construction | scirs2-interpolate | yes | 8 | 1.5x to 2.73x | `oned/pchip/n1000/shuffled`: 11.5 µs against 4.2 µs |
| `oned/bspline1` | construction | SciPy | yes | 2 | 2.68x to 2.71x | `oned/bspline1/n10000/shuffled`: 62.2 µs against 22.9 µs |
| `grid/nearest` | construction | scirs2-interpolate | yes | 3 | 1.24x to 2.36x | `grid/nearest/3d/each40/eval10000/ascending`: 21 µs against 8.88 µs |
| `oned/cubic_spline` | construction | scirs2-interpolate | no | 3 | 1.13x to 2.23x | `oned/cubic_spline/n1000/shuffled`: 37.5 µs against 16.8 µs |
| `oned/linear` | construction | SciPy | yes | 4 | 1.89x to 1.97x | `oned/linear/n10000/sorted`: 88.7 µs against 45 µs |
| `oned/akima` | construction | scirs2-interpolate | yes | 6 | 1.32x to 1.84x | `oned/akima/n1000/shuffled`: 10.2 µs against 5.54 µs |
| `oned/linear` | construction | scirs2-interpolate | yes | 13 | 1.13x to 1.82x | `oned/linear/n1000/shuffled`: 2.12 µs against 1.17 µs |
| `grid/linear` | evaluation | interpn | yes | 14 | 1.12x to 1.48x | `grid/linear/2d/max10/eval100/ascending`: 1.37 µs against 927 ns |
| `oned/nearest` | construction | SciPy | yes | 2 | 1.17x to 1.17x | `oned/nearest/n10000/shuffled`: 62.4 µs against 53.2 µs |
| `grid/nearest` | evaluation | interpn | yes | 1 | 1.13x | `grid/nearest/2d/each40/eval10000/ascending`: 150 µs against 133 µs |
| `oned/nearest` | evaluation | scirs2-interpolate | yes | 1 | 1.12x | `oned/nearest/n10/shuffled`: 297 ns against 265 ns |
| `oned/nearest` | construction | scirs2-interpolate | yes | 1 | 1.1x | `oned/nearest/n10000/shuffled`: 62.4 µs against 56.5 µs |

<details><summary>Every comparison where interpx is slower</summary>

"Values apart" is the distance of that library's values from SciPy's, as a share of the largest value.

| case | operation | against | interpx | the other | slower by | values apart |
| :--- | :--- | :--- | ---: | ---: | ---: | ---: |
| `grid/cubic/2d/auckland-ba32-0401/heldout/ascending` | construction | interpn | 79.5 ms | 0 | 90446051.56x | 5e-03 |
| `grid/linear/4d/max200/eval1000/ascending` | construction | interpn | 47.9 ms | 1.07 ns | 44580491.64x | 3e-16 |
| `grid/linear/4d/max200/eval10000/ascending` | construction | interpn | 46.4 ms | 1.07 ns | 43178263.27x | 4e-16 |
| `grid/linear/4d/max200/eval100/ascending` | construction | interpn | 45.7 ms | 1.17 ns | 39025408.00x | 3e-16 |
| `grid/linear/4d/max200/eval10/ascending` | construction | interpn | 43.2 ms | 1.17 ns | 36884736.00x | 2e-16 |
| `grid/linear/2d/auckland-ba32-0401/heldout/ascending` | construction | interpn | 15.9 ms | 0 | 23216128.00x | 0e+00 |
| `grid/nearest/2d/auckland-ba32-0401/heldout/ascending` | construction | interpn | 15.3 ms | 0 | 22439936.00x | 0e+00 |
| `grid/linear/3d/max200/eval100/ascending` | construction | interpn | 1.75 ms | 0 | 1987754.67x | 3e-16 |
| `grid/linear/3d/max200/eval10000/ascending` | construction | interpn | 1.63 ms | 0 | 1851505.78x | 5e-16 |
| `grid/linear/3d/max200/eval1000/ascending` | construction | interpn | 1.77 ms | 0 | 1812992.00x | 3e-16 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | interpn | 1.09 ms | 0 | 1593856.00x | 3e-16 |
| `grid/linear/3d/max200/eval10/ascending` | construction | interpn | 1.54 ms | 0 | 1579724.80x | 3e-16 |
| `grid/cubic/2d/each400/eval1000/ascending` | construction | interpn | 1.38 ms | 0 | 1574286.22x | 5e-01 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | interpn | 311 µs | 0 | 454436.57x | 3e-16 |
| `grid/cubic/3d/each40/eval10000/ascending` | construction | interpn | 389 µs | 1.17 ns | 331818.67x | 5e-01 |
| `grid/cubic/2d/each200/eval1000/ascending` | construction | interpn | 172 µs | 0 | 196096.00x | 2e-01 |
| `grid/cubic/2d/each100/eval1000/ascending` | construction | interpn | 81.5 µs | 0 | 92672.00x | 2e-01 |
| `grid/linear/2d/max200/eval10000/ascending` | construction | interpn | 47.1 µs | 0 | 68918.86x | 3e-16 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | interpn | 21 µs | 0 | 21491.20x | 0e+00 |
| `grid/linear/3d/max40/eval1000/ascending` | construction | interpn | 16.6 µs | 0 | 17011.20x | 2e-16 |
| `grid/linear/3d/max40/eval10000/ascending` | construction | interpn | 13.9 µs | 0 | 15854.22x | 3e-16 |
| `grid/linear/4d/max40/eval1000/ascending` | construction | interpn | 13.2 µs | 1.17 ns | 11221.33x | 6e-16 |
| `grid/linear/4d/max40/eval10/ascending` | construction | interpn | 13.1 µs | 1.17 ns | 11178.67x | 2e-16 |
| `grid/linear/4d/max40/eval10000/ascending` | construction | interpn | 13.1 µs | 1.17 ns | 11178.67x | 5e-16 |
| `grid/linear/4d/max40/eval100/ascending` | construction | interpn | 13 µs | 1.17 ns | 11077.33x | 2e-16 |
| `grid/cubic/2d/each40/eval1000/ascending` | construction | interpn | 8.46 µs | 0 | 10830.00x | 7e-01 |
| `grid/linear/2d/max200/eval100/ascending` | construction | interpn | 6.84 µs | 0 | 10006.86x | 1e-16 |
| `grid/linear/2d/max200/eval10/ascending` | construction | interpn | 6.83 µs | 0 | 9993.14x | 2e-16 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | interpn | 6.78 µs | 0 | 9924.57x | 2e-16 |
| `grid/cubic/2d/each40/eval10000/ascending` | construction | interpn | 8.48 µs | 0 | 9653.33x | 7e-01 |
| `grid/linear/4d/max200/eval1000/ascending` | construction | SciPy | 47.9 ms | 9.27 µs | 5168.48x |  |
| `grid/linear/4d/max200/eval10000/ascending` | construction | SciPy | 46.4 ms | 9.23 µs | 5024.55x |  |
| `grid/linear/4d/max200/eval10000/flipped` | construction | SciPy | 99.6 ms | 19.9 µs | 4992.37x |  |
| `grid/linear/4d/max200/eval1000/flipped` | construction | SciPy | 99.1 ms | 20.3 µs | 4883.78x |  |
| `grid/linear/4d/max200/eval10/flipped` | construction | SciPy | 97.6 ms | 20.1 µs | 4845.21x |  |
| `grid/linear/4d/max200/eval100/ascending` | construction | SciPy | 45.7 ms | 9.44 µs | 4842.66x |  |
| `grid/linear/4d/max200/eval100/flipped` | construction | SciPy | 97.5 ms | 20.1 µs | 4840.41x |  |
| `grid/linear/4d/max200/eval10/ascending` | construction | SciPy | 43.2 ms | 9.24 µs | 4679.22x |  |
| `grid/linear/3d/max40/eval100/ascending` | construction | interpn | 3.04 µs | 0 | 3455.11x | 2e-16 |
| `grid/linear/3d/max40/eval10/ascending` | construction | interpn | 2.96 µs | 0 | 3362.67x | 2e-16 |
| `grid/nearest/2d/auckland-ba32-0401/heldout/ascending` | construction | SciPy | 15.3 ms | 7.04 µs | 2178.76x |  |
| `grid/linear/2d/auckland-ba32-0401/heldout/ascending` | construction | SciPy | 15.9 ms | 7.53 µs | 2108.15x |  |
| `grid/linear/2d/auckland-ba32-0401/heldout/ascending` | construction | ndarray-interp | 15.9 ms | 9.05 µs | 1753.64x | 0e+00 |
| `grid/cubic/2d/each10/eval1000/ascending` | construction | interpn | 1.48 µs | 0 | 1680.00x | 3e-01 |
| `grid/nearest/2d/each40/eval10000/ascending` | construction | interpn | 633 ns | 0 | 926.29x | 0e+00 |
| `grid/linear/2d/max40/eval10/ascending` | construction | interpn | 397 ns | 0 | 580.57x | 2e-16 |
| `grid/linear/2d/max40/eval10000/ascending` | construction | interpn | 364 ns | 0 | 532.57x | 4e-16 |
| `grid/linear/2d/max40/eval100/ascending` | construction | interpn | 360 ns | 0 | 526.29x | 2e-16 |
| `grid/linear/2d/max40/eval1000/ascending` | construction | interpn | 359 ns | 0 | 525.86x | 3e-16 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | ndarray-interp | 1.09 ms | 2.22 µs | 489.86x | 6e-16 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | ndarray-interp | 311 µs | 1.12 µs | 277.87x | 5e-16 |
| `grid/linear/3d/max200/eval100/ascending` | construction | SciPy | 1.75 ms | 7.6 µs | 229.78x |  |
| `grid/linear/3d/max200/eval10000/ascending` | construction | SciPy | 1.63 ms | 7.54 µs | 215.85x |  |
| `grid/linear/3d/max200/eval1000/ascending` | construction | SciPy | 1.77 ms | 8.35 µs | 212.12x |  |
| `grid/linear/3d/max10/eval1000/ascending` | construction | interpn | 183 ns | 0 | 208.06x | 5e-16 |
| `grid/linear/3d/max10/eval10000/ascending` | construction | interpn | 180 ns | 0 | 205.22x | 4e-16 |
| `grid/linear/3d/max10/eval100/ascending` | construction | interpn | 180 ns | 0 | 204.89x | 3e-16 |
| `grid/linear/3d/max10/eval10/ascending` | construction | interpn | 180 ns | 0 | 204.44x | 2e-16 |
| `grid/linear/3d/max200/eval1000/flipped` | construction | SciPy | 3.22 ms | 15.8 µs | 204.22x |  |
| `grid/linear/3d/max200/eval10/ascending` | construction | SciPy | 1.54 ms | 7.56 µs | 203.99x |  |
| `grid/linear/2d/max10/eval1000/ascending` | construction | interpn | 135 ns | 0 | 197.14x | 3e-16 |
| `grid/linear/3d/max200/eval10000/flipped` | construction | SciPy | 3.12 ms | 15.8 µs | 197.00x |  |
| `grid/linear/2d/max10/eval100/ascending` | construction | interpn | 134 ns | 0 | 195.86x | 2e-16 |
| `grid/linear/2d/max10/eval10/ascending` | construction | interpn | 134 ns | 0 | 195.43x | 1e-16 |
| `grid/linear/3d/max200/eval10/flipped` | construction | SciPy | 3.02 ms | 16.1 µs | 188.45x |  |
| `grid/linear/2d/max10/eval10000/ascending` | construction | interpn | 122 ns | 0 | 178.57x | 3e-16 |
| `grid/linear/3d/max200/eval100/flipped` | construction | SciPy | 3.12 ms | 17.5 µs | 178.45x |  |
| `grid/linear/2d/each800/eval10000/ascending` | construction | SciPy | 1.09 ms | 6.32 µs | 172.52x |  |
| `grid/linear/2d/max200/eval10000/ascending` | construction | ndarray-interp | 47.1 µs | 463 ns | 101.78x | 5e-16 |
| `grid/pchip/3d/each40/eval10000/ascending` | construction | SciPy | 508 µs | 8.38 µs | 60.61x |  |
| `grid/linear/2d/each400/eval10000/ascending` | construction | SciPy | 311 µs | 5.88 µs | 52.79x |  |
| `grid/linear/2d/max200/eval100/ascending` | construction | ndarray-interp | 6.84 µs | 472 ns | 14.50x | 3e-16 |
| `grid/linear/2d/max200/eval10/ascending` | construction | ndarray-interp | 6.83 µs | 472 ns | 14.48x | 2e-16 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | ndarray-interp | 6.78 µs | 471 ns | 14.40x | 4e-16 |
| `grid/linear/2d/max200/eval10000/ascending` | construction | SciPy | 47.1 µs | 5.75 µs | 8.20x |  |
| `oned/linear/n10000/shuffled-uneven` | construction | ndarray-interp | 88.8 µs | 14.6 µs | 6.10x | 0e+00 |
| `oned/linear/n10000/shuffled` | construction | ndarray-interp | 89.3 µs | 19.9 µs | 4.49x | 0e+00 |
| `grid/linear/2d/max200/eval10000/flipped` | construction | SciPy | 53 µs | 11.9 µs | 4.46x |  |
| `oned/linear/n10000/sorted` | construction | ndarray-interp | 88.7 µs | 21.4 µs | 4.15x | 0e+00 |
| `oned/linear/n10000/sorted-uneven` | construction | ndarray-interp | 88.9 µs | 31.4 µs | 2.83x | 0e+00 |
| `grid/linear/4d/max40/eval10000/flipped` | construction | SciPy | 54.2 µs | 19.3 µs | 2.80x |  |
| `grid/linear/3d/max40/eval1000/ascending` | construction | scirs2-interpolate | 16.6 µs | 5.96 µs | 2.79x | 2e-16 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | SciPy | 21 µs | 7.54 µs | 2.78x |  |
| `oned/pchip/n1000/shuffled` | construction | scirs2-interpolate | 11.5 µs | 4.2 µs | 2.73x | 2e-16 |
| `oned/bspline1/n10000/shuffled` | construction | SciPy | 62.2 µs | 22.9 µs | 2.71x |  |
| `oned/bspline1/n10000/sorted` | construction | SciPy | 61.3 µs | 22.9 µs | 2.68x |  |
| `grid/linear/4d/max40/eval100/flipped` | construction | SciPy | 51.9 µs | 19.5 µs | 2.66x |  |
| `grid/linear/4d/max40/eval10/flipped` | construction | SciPy | 52 µs | 19.6 µs | 2.65x |  |
| `grid/linear/4d/max40/eval1000/flipped` | construction | SciPy | 52.7 µs | 19.9 µs | 2.64x |  |
| `grid/linear/2d/max40/eval10/ascending` | construction | ndarray-interp | 397 ns | 152 ns | 2.61x | 2e-16 |
| `oned/pchip/n1000/sorted` | construction | scirs2-interpolate | 11.5 µs | 4.56 µs | 2.53x | 2e-16 |
| `oned/linear/n10/shuffled` | construction | ndarray-interp | 98.4 ns | 40.5 ns | 2.43x | 0e+00 |
| `grid/linear/4d/max40/eval10/ascending` | construction | scirs2-interpolate | 13.1 µs | 5.4 µs | 2.43x | 2e-16 |
| `grid/linear/2d/max40/eval10000/ascending` | construction | ndarray-interp | 364 ns | 153 ns | 2.38x | 4e-16 |
| `grid/linear/4d/max40/eval10000/ascending` | construction | scirs2-interpolate | 13.1 µs | 5.52 µs | 2.37x | 4e-16 |
| `grid/linear/4d/max40/eval100/ascending` | construction | scirs2-interpolate | 13 µs | 5.48 µs | 2.37x | 3e-16 |
| `grid/linear/3d/max40/eval10000/ascending` | construction | scirs2-interpolate | 13.9 µs | 5.88 µs | 2.37x | 2e-16 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | scirs2-interpolate | 21 µs | 8.88 µs | 2.36x | 0e+00 |
| `grid/linear/4d/max40/eval1000/ascending` | construction | scirs2-interpolate | 13.2 µs | 5.58 µs | 2.36x | 2e-16 |
| `grid/linear/2d/max40/eval1000/ascending` | construction | ndarray-interp | 359 ns | 153 ns | 2.36x | 4e-16 |
| `grid/linear/2d/max40/eval100/ascending` | construction | ndarray-interp | 360 ns | 154 ns | 2.34x | 2e-16 |
| `grid/linear/2d/max200/eval100/ascending` | construction | scirs2-interpolate | 6.84 µs | 2.93 µs | 2.34x | 2e-16 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | scirs2-interpolate | 6.78 µs | 2.96 µs | 2.29x | 2e-16 |
| `oned/cubic_spline/n1000/shuffled` | construction | scirs2-interpolate | 37.5 µs | 16.8 µs | 2.23x | 4e-04 |
| `oned/cubic_spline/n1000/sorted` | construction | scirs2-interpolate | 37.5 µs | 16.8 µs | 2.23x | 4e-04 |
| `grid/nearest/2d/each40/eval10000/ascending` | construction | scirs2-interpolate | 633 ns | 288 ns | 2.20x | 0e+00 |
| `grid/linear/2d/max200/eval10/ascending` | construction | scirs2-interpolate | 6.83 µs | 3.11 µs | 2.20x | 1e-16 |
| `oned/linear/n10/sorted` | construction | ndarray-interp | 93.2 ns | 44.1 ns | 2.11x | 0e+00 |
| `oned/linear/n10000/sorted` | construction | SciPy | 88.7 µs | 45 µs | 1.97x |  |
| `grid/linear/2d/max40/eval10/ascending` | construction | scirs2-interpolate | 397 ns | 205 ns | 1.94x | 4e-17 |
| `oned/linear/n10000/shuffled-uneven` | construction | SciPy | 88.8 µs | 46.8 µs | 1.90x |  |
| `oned/linear/n10000/sorted-uneven` | construction | SciPy | 88.9 µs | 47 µs | 1.89x |  |
| `oned/linear/n10000/shuffled` | construction | SciPy | 89.3 µs | 47.2 µs | 1.89x |  |
| `oned/akima/n1000/shuffled` | construction | scirs2-interpolate | 10.2 µs | 5.54 µs | 1.84x | 3e-17 |
| `grid/linear/3d/max40/eval1000/flipped` | construction | SciPy | 28.9 µs | 15.7 µs | 1.84x |  |
| `oned/linear/n1000/shuffled` | construction | scirs2-interpolate | 2.12 µs | 1.17 µs | 1.82x | 0e+00 |
| `grid/linear/2d/max40/eval10000/ascending` | construction | scirs2-interpolate | 364 ns | 202 ns | 1.80x | 2e-16 |
| `grid/linear/2d/max40/eval100/ascending` | construction | scirs2-interpolate | 360 ns | 200 ns | 1.80x | 3e-16 |
| `grid/linear/2d/max40/eval1000/ascending` | construction | scirs2-interpolate | 359 ns | 201 ns | 1.79x | 3e-16 |
| `grid/linear/2d/max10/eval1000/ascending` | construction | ndarray-interp | 135 ns | 78 ns | 1.73x | 3e-16 |
| `grid/linear/2d/max10/eval10/ascending` | construction | ndarray-interp | 134 ns | 78.1 ns | 1.71x | 1e-16 |
| `oned/pchip/n100/sorted-uneven` | construction | scirs2-interpolate | 887 ns | 535 ns | 1.66x | 2e-16 |
| `oned/pchip/n100/shuffled-uneven` | construction | scirs2-interpolate | 887 ns | 540 ns | 1.64x | 2e-16 |
| `oned/pchip/n100/sorted` | construction | scirs2-interpolate | 886 ns | 541 ns | 1.64x | 2e-16 |
| `oned/linear/n1000/sorted` | construction | scirs2-interpolate | 1.94 µs | 1.2 µs | 1.62x | 0e+00 |
| `oned/pchip/n100/shuffled` | construction | scirs2-interpolate | 868 ns | 539 ns | 1.61x | 2e-16 |
| `oned/pchip/n50/shuffled` | construction | scirs2-interpolate | 535 ns | 333 ns | 1.61x | 1e-16 |
| `oned/linear/n10000/sorted-uneven` | construction | scirs2-interpolate | 88.9 µs | 56.2 µs | 1.58x | 6e-17 |
| `oned/linear/n10000/shuffled` | construction | scirs2-interpolate | 89.3 µs | 56.6 µs | 1.58x | 0e+00 |
| `grid/linear/2d/max10/eval100/ascending` | construction | ndarray-interp | 134 ns | 85.1 ns | 1.57x | 4e-16 |
| `grid/linear/2d/max10/eval10000/ascending` | construction | ndarray-interp | 122 ns | 78.3 ns | 1.56x | 4e-16 |
| `oned/linear/n10000/shuffled-uneven` | construction | scirs2-interpolate | 88.8 µs | 57.2 µs | 1.55x | 6e-17 |
| `oned/linear/n1000/shuffled` | construction | ndarray-interp | 2.12 µs | 1.38 µs | 1.54x | 0e+00 |
| `oned/akima/n1000/sorted` | construction | scirs2-interpolate | 10.2 µs | 6.7 µs | 1.53x | 3e-17 |
| `oned/linear/n10000/sorted` | construction | scirs2-interpolate | 88.7 µs | 58.4 µs | 1.52x | 0e+00 |
| `oned/linear/n100/shuffled-uneven` | construction | ndarray-interp | 248 ns | 164 ns | 1.51x | 0e+00 |
| `oned/linear/n50/shuffled` | construction | ndarray-interp | 153 ns | 101 ns | 1.51x | 0e+00 |
| `oned/pchip/n50/sorted` | construction | scirs2-interpolate | 538 ns | 358 ns | 1.50x | 1e-16 |
| `oned/linear/n50/sorted` | construction | ndarray-interp | 152 ns | 102 ns | 1.50x | 0e+00 |
| `oned/linear/n100/sorted-uneven` | construction | ndarray-interp | 248 ns | 167 ns | 1.49x | 0e+00 |
| `oned/linear/n100/sorted` | construction | ndarray-interp | 246 ns | 167 ns | 1.48x | 0e+00 |
| `grid/linear/2d/max10/eval100/ascending` | evaluation | interpn | 1.37 µs | 927 ns | 1.48x | 2e-16 |
| `oned/linear/n100/shuffled` | construction | ndarray-interp | 246 ns | 167 ns | 1.47x | 0e+00 |
| `grid/linear/3d/max200/eval100/ascending` | construction | scirs2-interpolate | 1.75 ms | 1.19 ms | 1.47x | 3e-16 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | scirs2-interpolate | 1.09 ms | 759 µs | 1.43x | 2e-16 |
| `grid/linear/2d/max10/eval1000/ascending` | evaluation | interpn | 13.2 µs | 9.25 µs | 1.43x | 3e-16 |
| `grid/linear/2d/max10/eval10/ascending` | evaluation | interpn | 135 ns | 94.6 ns | 1.43x | 1e-16 |
| `grid/linear/4d/max40/eval1000/ascending` | construction | SciPy | 13.2 µs | 9.28 µs | 1.42x |  |
| `grid/linear/4d/max40/eval10000/ascending` | construction | SciPy | 13.1 µs | 9.31 µs | 1.41x |  |
| `grid/linear/4d/max40/eval10/ascending` | construction | SciPy | 13.1 µs | 9.32 µs | 1.40x |  |
| `oned/akima/n100/shuffled` | construction | scirs2-interpolate | 737 ns | 528 ns | 1.40x | 1e-16 |
| `grid/linear/4d/max40/eval100/ascending` | construction | SciPy | 13 µs | 9.32 µs | 1.39x |  |
| `grid/linear/3d/max200/eval1000/ascending` | construction | scirs2-interpolate | 1.77 ms | 1.27 ms | 1.39x | 3e-16 |
| `oned/akima/n100/sorted` | construction | scirs2-interpolate | 740 ns | 534 ns | 1.39x | 1e-16 |
| `grid/linear/3d/max10/eval10/ascending` | evaluation | interpn | 195 ns | 141 ns | 1.38x | 2e-16 |
| `grid/linear/3d/max10/eval1000/ascending` | evaluation | interpn | 19 µs | 13.9 µs | 1.37x | 5e-16 |
| `grid/linear/3d/max10/eval100/ascending` | evaluation | interpn | 1.91 µs | 1.4 µs | 1.36x | 3e-16 |
| `oned/akima/n50/sorted` | construction | scirs2-interpolate | 434 ns | 320 ns | 1.36x | 1e-16 |
| `grid/linear/3d/max200/eval10000/ascending` | construction | scirs2-interpolate | 1.63 ms | 1.2 ms | 1.36x | 3e-16 |
| `grid/linear/2d/max10/eval10000/ascending` | evaluation | interpn | 125 µs | 92.3 µs | 1.35x | 3e-16 |
| `grid/linear/3d/max10/eval10000/ascending` | evaluation | interpn | 190 µs | 142 µs | 1.34x | 4e-16 |
| `oned/linear/n100/shuffled` | construction | scirs2-interpolate | 246 ns | 185 ns | 1.33x | 0e+00 |
| `oned/akima/n50/shuffled` | construction | scirs2-interpolate | 425 ns | 323 ns | 1.32x | 1e-16 |
| `oned/linear/n100/shuffled-uneven` | construction | scirs2-interpolate | 248 ns | 189 ns | 1.31x | 3e-17 |
| `oned/linear/n100/sorted` | construction | scirs2-interpolate | 246 ns | 188 ns | 1.31x | 0e+00 |
| `oned/linear/n100/sorted-uneven` | construction | scirs2-interpolate | 248 ns | 190 ns | 1.31x | 3e-17 |
| `grid/linear/3d/max200/eval10/ascending` | construction | scirs2-interpolate | 1.54 ms | 1.2 ms | 1.28x | 3e-16 |
| `oned/linear/n1000/sorted` | construction | ndarray-interp | 1.94 µs | 1.51 µs | 1.28x | 0e+00 |
| `grid/linear/2d/auckland-ba32-0401/heldout/ascending` | construction | scirs2-interpolate | 15.9 ms | 12.6 ms | 1.26x | 0e+00 |
| `grid/nearest/2d/auckland-ba32-0401/heldout/ascending` | construction | scirs2-interpolate | 15.3 ms | 12.4 ms | 1.24x | 0e+00 |
| `grid/linear/2d/max200/eval10/flipped` | construction | SciPy | 14.7 µs | 11.9 µs | 1.23x |  |
| `grid/linear/4d/max200/eval1000/ascending` | construction | scirs2-interpolate | 47.9 ms | 38.9 ms | 1.23x | 3e-16 |
| `grid/linear/2d/max200/eval1000/flipped` | construction | SciPy | 14.3 µs | 11.7 µs | 1.22x |  |
| `grid/linear/4d/max200/eval100/ascending` | construction | scirs2-interpolate | 45.7 ms | 37.7 ms | 1.21x | 1e-16 |
| `grid/linear/4d/max200/eval10/ascending` | construction | scirs2-interpolate | 43.2 ms | 35.9 ms | 1.20x | 2e-16 |
| `grid/linear/4d/max200/eval10000/ascending` | construction | scirs2-interpolate | 46.4 ms | 38.5 ms | 1.20x | 4e-16 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | scirs2-interpolate | 311 µs | 261 µs | 1.19x | 2e-16 |
| `grid/linear/2d/max200/eval100/ascending` | construction | SciPy | 6.84 µs | 5.77 µs | 1.18x |  |
| `oned/nearest/n10000/shuffled` | construction | SciPy | 62.4 µs | 53.2 µs | 1.17x |  |
| `oned/nearest/n10000/sorted` | construction | SciPy | 61.7 µs | 52.8 µs | 1.17x |  |
| `grid/linear/2d/max200/eval1000/ascending` | construction | SciPy | 6.78 µs | 5.82 µs | 1.17x |  |
| `oned/linear/n50/shuffled` | construction | scirs2-interpolate | 153 ns | 131 ns | 1.16x | 0e+00 |
| `grid/linear/2d/max40/eval10/ascending` | evaluation | interpn | 156 ns | 134 ns | 1.16x | 2e-16 |
| `oned/linear/n50/sorted` | construction | scirs2-interpolate | 152 ns | 132 ns | 1.15x | 0e+00 |
| `grid/linear/2d/max40/eval100/ascending` | evaluation | interpn | 1.51 µs | 1.32 µs | 1.14x | 2e-16 |
| `grid/linear/2d/max200/eval100/flipped` | construction | SciPy | 14.4 µs | 12.6 µs | 1.14x |  |
| `grid/linear/2d/max40/eval1000/ascending` | evaluation | interpn | 15.1 µs | 13.2 µs | 1.14x | 3e-16 |
| `grid/linear/2d/max200/eval10/ascending` | construction | SciPy | 6.83 µs | 6 µs | 1.14x |  |
| `oned/linear/n10/shuffled` | construction | scirs2-interpolate | 98.4 ns | 86.9 ns | 1.13x | 0e+00 |
| `grid/nearest/2d/each40/eval10000/ascending` | evaluation | interpn | 150 µs | 133 µs | 1.13x | 0e+00 |
| `oned/cubic_spline/n100/shuffled` | construction | scirs2-interpolate | 1.74 µs | 1.54 µs | 1.13x | 5e-03 |
| `grid/linear/2d/max40/eval10000/ascending` | evaluation | interpn | 151 µs | 134 µs | 1.13x | 4e-16 |
| `grid/linear/3d/max200/eval100/ascending` | evaluation | interpn | 2.95 µs | 2.64 µs | 1.12x | 3e-16 |
| `oned/nearest/n10/shuffled` | evaluation | scirs2-interpolate | 297 ns | 265 ns | 1.12x | 0e+00 |
| `grid/linear/4d/max40/eval1000/ascending` | evaluation | interpn | 29.6 µs | 26.5 µs | 1.12x | 6e-16 |
| `oned/nearest/n10000/shuffled` | construction | scirs2-interpolate | 62.4 µs | 56.5 µs | 1.10x | 0e+00 |
| `grid/linear/3d/max10/eval10000/ascending` | construction | scirs2-interpolate | 180 ns | 164 ns | 1.10x | 3e-16 |

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
| `rbf/thin_plate_spline/global/ghcn-tmax-20250101/heldout` | scirs2-interpolate | stopped at the time limit of 15 minutes |

## Real-world data

Two data sets, each taken whole in the unit its source gives it out in; what was taken and how it was divided was fixed before anything was measured (`compare/perf/real.py`, docs/design.ja.md section 8.4.7).

- `auckland-ba32-0401`: Auckland Part 1 LiDAR 1m DEM (2024), tile BA32_10000_0401. 4800 by 7200 cells of 1 m. The grid is every second row and column (and the last ones); every other cell is a query, and its elevation is the value held out. Licensed by Regional Software Holdings Limited for re-use under CC-BY-4.0. Hosted by Toitū Te Whenua Land Information New Zealand (https://registry.opendata.aws/nz-elevation).
- `ghcn-tmax-20250101`: GHCN-Daily, maximum temperature of 2025-01-01, every station; longitude and latitude in degrees, as they are. Of every 10 stations in the order of the file the last is a query, and its temperature is the value held out. The RBF kernels are those without a shape parameter. NOAA Global Historical Climatology Network Daily (GHCN-D), from https://registry.opendata.aws/noaa-ghcn (CC0 1.0).

These cases are also in the tables of their families below and in the counts above. Here they are with the error: the largest and the root mean square distance from the values held out (metres, and degrees Celsius), and the number of queries a library gave no value. A library missing from a case does not have it, or did not finish (see the runs that failed).

`grid/nearest/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 15.3 ms | 346 ms | 15.5 | 0.178 | 0 |
| interpx, parallel | 16.9 ms | 58.8 ms | 15.5 | 0.178 | 0 |
| interpx, `from_owned` | 5.16 ms | 349 ms | 15.5 | 0.178 | 0 |
| SciPy | 7.04 µs | 1.31 s | 15.5 | 0.178 | 0 |
| scirs2-interpolate | 12.4 ms | 159 s | 15.5 | 0.178 | 0 |
| interpn | 0 | 723 ms | 15.5 | 0.178 | 0 |

`grid/linear/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 15.9 ms | 365 ms | 6.74 | 0.0601 | 0 |
| interpx, parallel | 15.5 ms | 65 ms | 6.74 | 0.0601 | 0 |
| interpx, `from_owned` | 5.19 ms | 365 ms | 6.74 | 0.0601 | 0 |
| SciPy | 7.53 µs | 1.08 s | 6.74 | 0.0601 | 0 |
| scirs2-interpolate | 12.6 ms | 67 s | 6.74 | 0.0601 | 0 |
| ndarray-interp | 9.05 µs | 344 ms | 6.74 | 0.0601 | 0 |
| interpn | 0 | 754 ms | 6.74 | 0.0601 | 0 |

`grid/cubic/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 79.5 ms | 724 ms | 6.77 | 0.0507 | 0 |
| interpx, parallel | 76.5 ms | 132 ms | 6.77 | 0.0507 | 0 |
| interpx, `from_owned` | 65.9 ms | 712 ms | 6.77 | 0.0507 | 0 |
| SciPy | 13.9 s | 32.4 s | 6.77 | 0.0507 | 0 |
| interpn | 0 | 1.6 s | 6.76 | 0.0517 | 0 |

`grid/quintic/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 93 ms | 1.18 s | 6.77 | 0.0517 | 0 |
| interpx, parallel | 90.9 ms | 203 ms | 6.77 | 0.0517 | 0 |
| interpx, `from_owned` | 81.5 ms | 1.18 s | 6.77 | 0.0517 | 0 |
| SciPy | 21.4 s | 34.2 s | 6.77 | 0.0517 | 0 |

`grid/pchip/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 74.8 ms | 1.02 s | 6.75 | 0.049 | 0 |
| interpx, parallel | 74.1 ms | 170 ms | 6.75 | 0.049 | 0 |
| interpx, `from_owned` | 62 ms | 1.03 s | 6.75 | 0.049 | 0 |

`nearest_nd/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 653 µs | 159 µs | 25 | 3.42 | 0 |
| interpx, parallel | 669 µs | 62.6 µs | 25 | 3.42 | 0 |
| interpx, `from_owned` | 667 µs | 172 µs | 25 | 3.42 | 0 |
| SciPy | 1.78 ms | 511 µs | 25 | 3.42 | 0 |

`rbf/linear/global/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 16.6 s | 17.8 ms | 20.8 | 2.74 | 0 |
| interpx, parallel | 16.4 s | 3.69 ms | 20.8 | 2.74 | 0 |
| interpx, `from_owned` | 16.8 s | 18.1 ms | 20.8 | 2.74 | 0 |
| SciPy | 112 s | 52.9 ms | 20.8 | 2.74 | 0 |
| scirs2-interpolate | 272 s | 251 ms | 20.8 | 2.74 | 0 |

`rbf/linear/k50/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 750 µs | 17.8 ms | 20.8 | 2.75 | 0 |
| interpx, parallel | 764 µs | 3.31 ms | 20.8 | 2.75 | 0 |
| interpx, `from_owned` | 698 µs | 17.7 ms | 20.8 | 2.75 | 0 |
| SciPy | 1.71 ms | 41.9 ms | 20.8 | 2.75 | 0 |

`rbf/thin_plate_spline/global/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 18.3 s | 90.2 ms | 102 | 5.67 | 0 |
| interpx, parallel | 17.6 s | 19.6 ms | 102 | 5.67 | 0 |
| interpx, `from_owned` | 18 s | 83.8 ms | 102 | 5.67 | 0 |
| SciPy | 114 s | 126 ms | 102 | 5.67 | 0 |

`rbf/thin_plate_spline/k50/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 780 µs | 28.4 ms | 102 | 5.66 | 0 |
| interpx, parallel | 762 µs | 5.13 ms | 102 | 5.66 | 0 |
| interpx, `from_owned` | 703 µs | 28.3 ms | 102 | 5.66 | 0 |
| SciPy | 1.71 ms | 52.9 ms | 102 | 5.66 | 0 |

`rbf/cubic/global/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 18.7 s | 20.4 ms | 587 | 27.1 | 0 |
| interpx, parallel | 17.5 s | 3.78 ms | 587 | 27.1 | 0 |
| interpx, `from_owned` | 17.7 s | 18.7 ms | 587 | 27.1 | 0 |
| SciPy | 124 s | 57.8 ms | 597 | 24.4 | 0 |
| scirs2-interpolate | 285 s | 275 ms | 100 | 8.48 | 0 |

`rbf/cubic/k50/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 801 µs | 19.7 ms | 567 | 25.3 | 0 |
| interpx, parallel | 803 µs | 3.71 ms | 567 | 25.3 | 0 |
| interpx, `from_owned` | 705 µs | 20.3 ms | 567 | 25.3 | 0 |
| SciPy | 1.69 ms | 45.3 ms | 567 | 25.3 | 0 |

`rbf/quintic/global/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 16.3 s | 20 ms | 2.55e+05 | 1.36e+04 | 0 |
| interpx, parallel | 16.6 s | 4.07 ms | 2.55e+05 | 1.36e+04 | 0 |
| interpx, `from_owned` | 19.2 s | 22.3 ms | 2.55e+05 | 1.36e+04 | 0 |
| SciPy | 119 s | 50.5 ms | 9.39e+04 | 4.92e+03 | 0 |
| scirs2-interpolate | 275 s | 301 ms | 6.56e+16 | 2.85e+15 | 0 |

`rbf/quintic/k50/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 827 µs | 21.7 ms | 3.74e+03 | 199 | 0 |
| interpx, parallel | 758 µs | 3.97 ms | 3.74e+03 | 199 | 0 |
| interpx, `from_owned` | 724 µs | 23.2 ms | 3.74e+03 | 199 | 0 |
| SciPy | 1.7 ms | 50.5 ms | 3.74e+03 | 200 | 0 |

## One variable

Data of `n` points, evaluated at `4 n` points: in order (`sorted`) and in no order (`shuffled`).

### Evaluation

| case | interpx | interpx, parallel | SciPy | scirs2-interpolate | ndarray-interp | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/n10/sorted` | 111 ns | 189 ns | 4.89 µs | 283 ns | 256 ns | 2.3x faster |
| `linear/n10/shuffled` | 187 ns | 224 ns | 5.05 µs | 269 ns | 273 ns | 1.4x faster |
| `nearest/n10/sorted` | 148 ns | 209 ns | 9.55 µs | 256 ns |  | 1.7x faster |
| `nearest/n10/shuffled` | 297 ns | 359 ns | 9.4 µs | 265 ns |  | **1.1x slower** |
| `previous/n10/sorted` | 162 ns | 192 ns | 10.6 µs |  |  | 66x faster |
| `previous/n10/shuffled` | 285 ns | 320 ns | 10.7 µs |  |  | 37x faster |
| `bspline1/n10/sorted` | 242 ns | 319 ns | 2.75 µs |  |  | 11x faster |
| `bspline1/n10/shuffled` | 527 ns | 595 ns | 2.95 µs |  |  | 5.6x faster |
| `bspline2/n10/sorted` | 297 ns | 342 ns | 3.07 µs |  |  | 10x faster |
| `bspline2/n10/shuffled` | 340 ns | 415 ns | 3.06 µs |  |  | 9x faster |
| `bspline3/n10/sorted` | 344 ns | 432 ns | 3.62 µs |  |  | 11x faster |
| `bspline3/n10/shuffled` | 431 ns | 564 ns | 3.56 µs |  |  | 8.3x faster |
| `cubic_spline/n10/sorted` | 143 ns | 193 ns | 3.77 µs | 297 ns | 332 ns | 2.1x faster |
| `cubic_spline/n10/shuffled` | 199 ns | 275 ns | 3.76 µs | 278 ns | 325 ns | 1.4x faster |
| `pchip/n10/sorted` | 143 ns | 192 ns | 3.79 µs | 291 ns |  | 2x faster |
| `pchip/n10/shuffled` | 205 ns | 264 ns | 3.8 µs | 290 ns |  | 1.4x faster |
| `akima/n10/sorted` | 156 ns | 197 ns | 3.81 µs | 265 ns |  | 1.7x faster |
| `akima/n10/shuffled` | 200 ns | 280 ns | 3.83 µs | 287 ns |  | 1.4x faster |
| `makima/n10/sorted` | 144 ns | 194 ns | 3.79 µs |  |  | 26x faster |
| `makima/n10/shuffled` | 204 ns | 257 ns | 3.74 µs |  |  | 18x faster |
| `linear/n50/sorted` | 635 ns | 677 ns | 6.11 µs | 1.56 µs | 1.23 µs | 1.9x faster |
| `linear/n50/shuffled` | 812 ns | 932 ns | 8.07 µs | 1.56 µs | 1.26 µs | 1.6x faster |
| `nearest/n50/sorted` | 741 ns | 759 ns | 11.5 µs | 1.59 µs |  | 2.1x faster |
| `nearest/n50/shuffled` | 1.53 µs | 1.59 µs | 11.4 µs | 1.59 µs |  | 1x faster |
| `previous/n50/sorted` | 723 ns | 771 ns | 12.8 µs |  |  | 18x faster |
| `previous/n50/shuffled` | 1.35 µs | 1.42 µs | 12.8 µs |  |  | 9.5x faster |
| `bspline1/n50/sorted` | 1.1 µs | 1.36 µs | 3.84 µs |  |  | 3.5x faster |
| `bspline1/n50/shuffled` | 2.78 µs | 2.87 µs | 4.83 µs |  |  | 1.7x faster |
| `bspline2/n50/sorted` | 1.36 µs | 1.58 µs | 6.12 µs |  |  | 4.5x faster |
| `bspline2/n50/shuffled` | 2.04 µs | 2.09 µs | 6.74 µs |  |  | 3.3x faster |
| `bspline3/n50/sorted` | 1.83 µs | 1.89 µs | 8.65 µs |  |  | 4.7x faster |
| `bspline3/n50/shuffled` | 2.8 µs | 2.87 µs | 9.25 µs |  |  | 3.3x faster |
| `cubic_spline/n50/sorted` | 708 ns | 749 ns | 5.31 µs | 4.24 µs | 1.58 µs | 2.2x faster |
| `cubic_spline/n50/shuffled` | 987 ns | 1.04 µs | 5.6 µs | 4.04 µs | 1.66 µs | 1.7x faster |
| `pchip/n50/sorted` | 745 ns | 753 ns | 5.38 µs | 3.91 µs |  | 5.2x faster |
| `pchip/n50/shuffled` | 973 ns | 1.13 µs | 5.69 µs | 3.62 µs |  | 3.7x faster |
| `akima/n50/sorted` | 702 ns | 739 ns | 5.16 µs | 4.12 µs |  | 5.9x faster |
| `akima/n50/shuffled` | 980 ns | 1.03 µs | 5.71 µs | 3.91 µs |  | 4x faster |
| `makima/n50/sorted` | 699 ns | 762 ns | 5.36 µs |  |  | 7.7x faster |
| `makima/n50/shuffled` | 964 ns | 1.01 µs | 5.59 µs |  |  | 5.8x faster |
| `linear/n100/sorted` | 1.31 µs | 1.38 µs | 7.12 µs | 3.59 µs | 2.46 µs | 1.9x faster |
| `linear/n100/shuffled` | 1.62 µs | 1.71 µs | 12.7 µs | 3.63 µs | 2.49 µs | 1.5x faster |
| `nearest/n100/sorted` | 1.51 µs | 1.57 µs | 13.2 µs | 3.33 µs |  | 2.2x faster |
| `nearest/n100/shuffled` | 3.02 µs | 3.2 µs | 13.4 µs | 3.36 µs |  | 1.1x faster |
| `previous/n100/sorted` | 1.48 µs | 1.54 µs | 14.9 µs |  |  | 10x faster |
| `previous/n100/shuffled` | 2.71 µs | 2.83 µs | 15.2 µs |  |  | 5.6x faster |
| `bspline1/n100/sorted` | 2.36 µs | 2.71 µs | 5.28 µs |  |  | 2.2x faster |
| `bspline1/n100/shuffled` | 5.53 µs | 5.78 µs | 11.7 µs |  |  | 2.1x faster |
| `bspline2/n100/sorted` | 2.65 µs | 3.3 µs | 9.76 µs |  |  | 3.7x faster |
| `bspline2/n100/shuffled` | 4.43 µs | 4.71 µs | 12.8 µs |  |  | 2.9x faster |
| `bspline3/n100/sorted` | 3.62 µs | 3.77 µs | 14.9 µs |  |  | 4.1x faster |
| `bspline3/n100/shuffled` | 6.12 µs | 6.88 µs | 17.8 µs |  |  | 2.9x faster |
| `cubic_spline/n100/sorted` | 1.44 µs | 1.6 µs | 6.97 µs | 16.9 µs | 3.15 µs | 2.2x faster |
| `cubic_spline/n100/shuffled` | 1.91 µs | 2.12 µs | 8.09 µs | 17.3 µs | 3.13 µs | 1.6x faster |
| `pchip/n100/sorted` | 1.42 µs | 1.51 µs | 7.04 µs | 15.4 µs |  | 5x faster |
| `pchip/n100/shuffled` | 1.95 µs | 2.05 µs | 7.77 µs | 15.8 µs |  | 4x faster |
| `akima/n100/sorted` | 1.4 µs | 1.53 µs | 6.86 µs | 16.5 µs |  | 4.9x faster |
| `akima/n100/shuffled` | 1.91 µs | 2.05 µs | 8.03 µs | 17.4 µs |  | 4.2x faster |
| `makima/n100/sorted` | 1.42 µs | 1.5 µs | 6.98 µs |  |  | 4.9x faster |
| `makima/n100/shuffled` | 1.94 µs | 2.25 µs | 8.05 µs |  |  | 4.1x faster |
| `linear/n1000/sorted` | 13.1 µs | 13.4 µs | 29 µs | 102 µs | 24.5 µs | 1.9x faster |
| `linear/n1000/shuffled` | 16.1 µs | 16.4 µs | 143 µs | 202 µs | 24.7 µs | 1.5x faster |
| `nearest/n1000/sorted` | 14.6 µs | 14.7 µs | 56.8 µs | 129 µs |  | 3.9x faster |
| `nearest/n1000/shuffled` | 30.9 µs | 34.1 µs | 198 µs | 223 µs |  | 6.4x faster |
| `previous/n1000/sorted` | 14.4 µs | 14.8 µs | 62.3 µs |  |  | 4.3x faster |
| `previous/n1000/shuffled` | 27.3 µs | 27.4 µs | 202 µs |  |  | 7.4x faster |
| `bspline1/n1000/sorted` | 26.5 µs | 25.3 µs | 30.8 µs |  |  | 1.2x faster |
| `bspline1/n1000/shuffled` | 56.2 µs | 56.6 µs | 522 µs |  |  | 9.3x faster |
| `bspline2/n1000/sorted` | 31.7 µs | 32 µs | 76.1 µs |  |  | 2.4x faster |
| `bspline2/n1000/shuffled` | 58.1 µs | 57.2 µs | 542 µs |  |  | 9.3x faster |
| `bspline3/n1000/sorted` | 36.4 µs | 40 µs | 129 µs |  |  | 3.5x faster |
| `bspline3/n1000/shuffled` | 79 µs | 79.2 µs | 567 µs |  |  | 7.2x faster |
| `cubic_spline/n1000/sorted` | 14 µs | 14.4 µs | 37.7 µs | 1.42 ms | 31.9 µs | 2.3x faster |
| `cubic_spline/n1000/shuffled` | 19.5 µs | 19.4 µs | 219 µs | 1.42 ms | 31.8 µs | 1.6x faster |
| `pchip/n1000/sorted` | 14 µs | 14.3 µs | 37.6 µs | 1.22 ms |  | 2.7x faster |
| `pchip/n1000/shuffled` | 19.6 µs | 19.9 µs | 214 µs | 1.2 ms |  | 11x faster |
| `akima/n1000/sorted` | 13.9 µs | 14.3 µs | 37.2 µs | 1.39 ms |  | 2.7x faster |
| `akima/n1000/shuffled` | 19.5 µs | 19.7 µs | 215 µs | 1.41 ms |  | 11x faster |
| `makima/n1000/sorted` | 13.9 µs | 14.1 µs | 37.1 µs |  |  | 2.7x faster |
| `makima/n1000/shuffled` | 19.6 µs | 19.9 µs | 217 µs |  |  | 11x faster |
| `linear/n10000/sorted` | 138 µs | 65.8 µs | 235 µs | 1.28 ms | 247 µs | 1.7x faster |
| `linear/n10000/shuffled` | 185 µs | 86.1 µs | 2.31 ms | 2.85 ms | 260 µs | 1.4x faster |
| `nearest/n10000/sorted` | 148 µs | 72 µs | 583 µs | 1.49 ms |  | 3.9x faster |
| `nearest/n10000/shuffled` | 379 µs | 95.5 µs | 2.94 ms | 3.01 ms |  | 7.7x faster |
| `previous/n10000/sorted` | 151 µs | 72.1 µs | 608 µs |  |  | 4x faster |
| `previous/n10000/shuffled` | 344 µs | 91.8 µs | 2.98 ms |  |  | 8.7x faster |
| `bspline1/n10000/sorted` | 264 µs | 98.4 µs | 288 µs |  |  | 1.1x faster |
| `bspline1/n10000/shuffled` | 697 µs | 139 µs | 47.6 ms |  |  | 68x faster |
| `bspline2/n10000/sorted` | 334 µs | 85.2 µs | 745 µs |  |  | 2.2x faster |
| `bspline2/n10000/shuffled` | 892 µs | 195 µs | 48.7 ms |  |  | 55x faster |
| `bspline3/n10000/sorted` | 370 µs | 101 µs | 1.26 ms |  |  | 3.4x faster |
| `bspline3/n10000/shuffled` | 1.13 ms | 245 µs | 48 ms |  |  | 42x faster |
| `cubic_spline/n10000/sorted` | 146 µs | 69.8 µs | 366 µs | 142 ms | 314 µs | 2.2x faster |
| `cubic_spline/n10000/shuffled` | 215 µs | 97 µs | 3.35 ms | 141 ms | 333 µs | 1.5x faster |
| `pchip/n10000/sorted` | 146 µs | 70.7 µs | 371 µs | 117 ms |  | 2.5x faster |
| `pchip/n10000/shuffled` | 215 µs | 93.7 µs | 3.44 ms | 123 ms |  | 16x faster |
| `akima/n10000/sorted` | 145 µs | 72.8 µs | 365 µs | 141 ms |  | 2.5x faster |
| `akima/n10000/shuffled` | 216 µs | 94.8 µs | 3.39 ms | 144 ms |  | 16x faster |
| `makima/n10000/sorted` | 147 µs | 69.9 µs | 370 µs |  |  | 2.5x faster |
| `makima/n10000/shuffled` | 215 µs | 93.8 µs | 3.39 ms |  |  | 16x faster |
| `linear/n100/sorted-uneven` | 1.18 µs | 1.27 µs | 7.34 µs | 3.32 µs | 3.7 µs | 2.8x faster |
| `linear/n100/shuffled-uneven` | 2.75 µs | 2.86 µs | 12.5 µs | 3.66 µs | 3.62 µs | 1.3x faster |
| `cubic_spline/n100/sorted-uneven` | 1.35 µs | 1.47 µs | 7.59 µs | 17 µs | 4.4 µs | 3.3x faster |
| `cubic_spline/n100/shuffled-uneven` | 2.99 µs | 3.12 µs | 8.11 µs | 17.1 µs | 4.42 µs | 1.5x faster |
| `pchip/n100/sorted-uneven` | 1.31 µs | 1.46 µs | 7.4 µs | 15.8 µs |  | 5.6x faster |
| `pchip/n100/shuffled-uneven` | 3.04 µs | 3.16 µs | 7.98 µs | 15.9 µs |  | 2.6x faster |
| `linear/n10000/sorted-uneven` | 208 µs | 88.5 µs | 261 µs | 1.25 ms | 967 µs | 1.3x faster |
| `linear/n10000/shuffled-uneven` | 616 µs | 156 µs | 2.31 ms | 2.81 ms | 1.72 ms | 2.8x faster |
| `cubic_spline/n10000/sorted-uneven` | 217 µs | 91.7 µs | 444 µs | 139 ms | 938 µs | 2x faster |
| `cubic_spline/n10000/shuffled-uneven` | 800 µs | 166 µs | 3.32 ms | 143 ms | 1.82 ms | 2.3x faster |
| `pchip/n10000/sorted-uneven` | 216 µs | 94.8 µs | 444 µs | 118 ms |  | 2.1x faster |
| `pchip/n10000/shuffled-uneven` | 742 µs | 173 µs | 3.43 ms | 119 ms |  | 4.6x faster |

### Construction

| case | interpx | interpx, parallel | SciPy | scirs2-interpolate | ndarray-interp | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/n10/sorted` | 93.2 ns | 92.7 ns | 8.64 µs | 87.3 ns | 44.1 ns | **2.1x slower** |
| `linear/n10/shuffled` | 98.4 ns | 92.7 ns | 8.68 µs | 86.9 ns | 40.5 ns | **2.4x slower** |
| `nearest/n10/sorted` | 81.4 ns | 81.7 ns | 9.48 µs | 81.5 ns |  | 1x faster |
| `nearest/n10/shuffled` | 81.8 ns | 75.5 ns | 9.46 µs | 89.2 ns |  | 1.1x faster |
| `previous/n10/sorted` | 82.5 ns | 75.6 ns | 9.17 µs |  |  | 1.1e+02x faster |
| `previous/n10/shuffled` | 82.2 ns | 75.7 ns | 9.22 µs |  |  | 1.1e+02x faster |
| `bspline1/n10/sorted` | 104 ns | 113 ns | 15.7 µs |  |  | 1.5e+02x faster |
| `bspline1/n10/shuffled` | 102 ns | 103 ns | 16.2 µs |  |  | 1.6e+02x faster |
| `bspline2/n10/sorted` | 509 ns | 508 ns | 36.4 µs |  |  | 72x faster |
| `bspline2/n10/shuffled` | 468 ns | 469 ns | 38.3 µs |  |  | 82x faster |
| `bspline3/n10/sorted` | 555 ns | 504 ns | 34.7 µs |  |  | 63x faster |
| `bspline3/n10/shuffled` | 516 ns | 521 ns | 35.1 µs |  |  | 68x faster |
| `cubic_spline/n10/sorted` | 238 ns | 239 ns | 85.8 µs | 408 ns | 502 ns | 1.7x faster |
| `cubic_spline/n10/shuffled` | 237 ns | 246 ns | 85.8 µs | 410 ns | 496 ns | 1.7x faster |
| `pchip/n10/sorted` | 179 ns | 184 ns | 102 µs | 211 ns |  | 1.2x faster |
| `pchip/n10/shuffled` | 178 ns | 180 ns | 101 µs | 209 ns |  | 1.2x faster |
| `akima/n10/sorted` | 188 ns | 186 ns | 75.8 µs | 193 ns |  | 1x faster |
| `akima/n10/shuffled` | 168 ns | 169 ns | 76.5 µs | 179 ns |  | 1.1x faster |
| `makima/n10/sorted` | 179 ns | 180 ns | 78.4 µs |  |  | 4.4e+02x faster |
| `makima/n10/shuffled` | 179 ns | 180 ns | 77.5 µs |  |  | 4.3e+02x faster |
| `linear/n50/sorted` | 152 ns | 154 ns | 8.87 µs | 132 ns | 102 ns | **1.5x slower** |
| `linear/n50/shuffled` | 153 ns | 154 ns | 8.86 µs | 131 ns | 101 ns | **1.5x slower** |
| `nearest/n50/sorted` | 106 ns | 106 ns | 9.75 µs | 131 ns |  | 1.2x faster |
| `nearest/n50/shuffled` | 116 ns | 107 ns | 10.4 µs | 135 ns |  | 1.2x faster |
| `previous/n50/sorted` | 105 ns | 107 ns | 9.54 µs |  |  | 91x faster |
| `previous/n50/shuffled` | 106 ns | 108 ns | 9.48 µs |  |  | 89x faster |
| `bspline1/n50/sorted` | 157 ns | 143 ns | 15.6 µs |  |  | 1e+02x faster |
| `bspline1/n50/shuffled` | 149 ns | 152 ns | 15.8 µs |  |  | 1.1e+02x faster |
| `bspline2/n50/sorted` | 1.83 µs | 1.84 µs | 39.3 µs |  |  | 22x faster |
| `bspline2/n50/shuffled` | 1.83 µs | 1.84 µs | 39 µs |  |  | 21x faster |
| `bspline3/n50/sorted` | 2.18 µs | 2.04 µs | 38.4 µs |  |  | 18x faster |
| `bspline3/n50/shuffled` | 2.04 µs | 2.04 µs | 39.4 µs |  |  | 19x faster |
| `cubic_spline/n50/sorted` | 837 ns | 838 ns | 86.6 µs | 871 ns | 1.09 µs | 1x faster |
| `cubic_spline/n50/shuffled` | 838 ns | 840 ns | 87.5 µs | 905 ns | 1.07 µs | 1.1x faster |
| `pchip/n50/sorted` | 538 ns | 491 ns | 101 µs | 358 ns |  | **1.5x slower** |
| `pchip/n50/shuffled` | 535 ns | 531 ns | 99.1 µs | 333 ns |  | **1.6x slower** |
| `akima/n50/sorted` | 434 ns | 426 ns | 75 µs | 320 ns |  | **1.4x slower** |
| `akima/n50/shuffled` | 425 ns | 427 ns | 75 µs | 323 ns |  | **1.3x slower** |
| `makima/n50/sorted` | 471 ns | 471 ns | 85 µs |  |  | 1.8e+02x faster |
| `makima/n50/shuffled` | 478 ns | 461 ns | 79 µs |  |  | 1.7e+02x faster |
| `linear/n100/sorted` | 246 ns | 254 ns | 9.03 µs | 188 ns | 167 ns | **1.5x slower** |
| `linear/n100/shuffled` | 246 ns | 248 ns | 9.26 µs | 185 ns | 167 ns | **1.5x slower** |
| `nearest/n100/sorted` | 161 ns | 162 ns | 9.93 µs | 185 ns |  | 1.1x faster |
| `nearest/n100/shuffled` | 161 ns | 169 ns | 9.7 µs | 188 ns |  | 1.2x faster |
| `previous/n100/sorted` | 161 ns | 161 ns | 9.86 µs |  |  | 61x faster |
| `previous/n100/shuffled` | 161 ns | 159 ns | 9.77 µs |  |  | 61x faster |
| `bspline1/n100/sorted` | 196 ns | 199 ns | 15.6 µs |  |  | 80x faster |
| `bspline1/n100/shuffled` | 196 ns | 197 ns | 15.6 µs |  |  | 80x faster |
| `bspline2/n100/sorted` | 3.49 µs | 3.57 µs | 46.6 µs |  |  | 13x faster |
| `bspline2/n100/shuffled` | 3.49 µs | 3.47 µs | 46.7 µs |  |  | 13x faster |
| `bspline3/n100/sorted` | 3.9 µs | 3.9 µs | 43 µs |  |  | 11x faster |
| `bspline3/n100/shuffled` | 3.83 µs | 4.28 µs | 43.2 µs |  |  | 11x faster |
| `cubic_spline/n100/sorted` | 1.58 µs | 1.58 µs | 86.8 µs | 1.51 µs | 1.87 µs | **1x slower** |
| `cubic_spline/n100/shuffled` | 1.74 µs | 1.6 µs | 87.3 µs | 1.54 µs | 1.87 µs | **1.1x slower** |
| `pchip/n100/sorted` | 886 ns | 870 ns | 102 µs | 541 ns |  | **1.6x slower** |
| `pchip/n100/shuffled` | 868 ns | 883 ns | 101 µs | 539 ns |  | **1.6x slower** |
| `akima/n100/sorted` | 740 ns | 738 ns | 75.4 µs | 534 ns |  | **1.4x slower** |
| `akima/n100/shuffled` | 737 ns | 743 ns | 75.8 µs | 528 ns |  | **1.4x slower** |
| `makima/n100/sorted` | 824 ns | 820 ns | 80 µs |  |  | 97x faster |
| `makima/n100/shuffled` | 825 ns | 898 ns | 78.7 µs |  |  | 95x faster |
| `linear/n1000/sorted` | 1.94 µs | 1.94 µs | 12.8 µs | 1.2 µs | 1.51 µs | **1.6x slower** |
| `linear/n1000/shuffled` | 2.12 µs | 1.93 µs | 12.7 µs | 1.17 µs | 1.38 µs | **1.8x slower** |
| `nearest/n1000/sorted` | 1.01 µs | 996 ns | 14.4 µs | 1.17 µs |  | 1.2x faster |
| `nearest/n1000/shuffled` | 984 ns | 992 ns | 14.1 µs | 1.21 µs |  | 1.2x faster |
| `previous/n1000/sorted` | 988 ns | 990 ns | 16.3 µs |  |  | 17x faster |
| `previous/n1000/shuffled` | 984 ns | 996 ns | 16.2 µs |  |  | 16x faster |
| `bspline1/n1000/sorted` | 1.06 µs | 1.04 µs | 17.3 µs |  |  | 16x faster |
| `bspline1/n1000/shuffled` | 1.13 µs | 1.06 µs | 16.9 µs |  |  | 15x faster |
| `bspline2/n1000/sorted` | 36.4 µs | 44.2 µs | 107 µs |  |  | 2.9x faster |
| `bspline2/n1000/shuffled` | 34.4 µs | 38.2 µs | 108 µs |  |  | 3.2x faster |
| `bspline3/n1000/sorted` | 68.6 µs | 69.4 µs | 125 µs |  |  | 1.8x faster |
| `bspline3/n1000/shuffled` | 65.6 µs | 66.8 µs | 125 µs |  |  | 1.9x faster |
| `cubic_spline/n1000/sorted` | 37.5 µs | 40.1 µs | 111 µs | 16.8 µs | 17.5 µs | **2.2x slower** |
| `cubic_spline/n1000/shuffled` | 37.5 µs | 37.4 µs | 111 µs | 16.8 µs | 16 µs | **2.3x slower** |
| `pchip/n1000/sorted` | 11.5 µs | 11.4 µs | 118 µs | 4.56 µs |  | **2.5x slower** |
| `pchip/n1000/shuffled` | 11.5 µs | 11.9 µs | 118 µs | 4.2 µs |  | **2.7x slower** |
| `akima/n1000/sorted` | 10.2 µs | 10.2 µs | 92.4 µs | 6.7 µs |  | **1.5x slower** |
| `akima/n1000/shuffled` | 10.2 µs | 9.16 µs | 92.8 µs | 5.54 µs |  | **1.8x slower** |
| `makima/n1000/sorted` | 12.3 µs | 10.8 µs | 97.4 µs |  |  | 7.9x faster |
| `makima/n1000/shuffled` | 10.7 µs | 12 µs | 98.6 µs |  |  | 9.2x faster |
| `linear/n10000/sorted` | 88.7 µs | 88.2 µs | 45 µs | 58.4 µs | 21.4 µs | **4.1x slower** |
| `linear/n10000/shuffled` | 89.3 µs | 84.4 µs | 47.2 µs | 56.6 µs | 19.9 µs | **4.5x slower** |
| `nearest/n10000/sorted` | 61.7 µs | 54.6 µs | 52.8 µs | 56.3 µs |  | **1.2x slower** |
| `nearest/n10000/shuffled` | 62.4 µs | 55.5 µs | 53.2 µs | 56.5 µs |  | **1.2x slower** |
| `previous/n10000/sorted` | 62.2 µs | 54 µs | 73.1 µs |  |  | 1.2x faster |
| `previous/n10000/shuffled` | 62.3 µs | 54.9 µs | 73.4 µs |  |  | 1.2x faster |
| `bspline1/n10000/sorted` | 61.3 µs | 55.3 µs | 22.9 µs |  |  | **2.7x slower** |
| `bspline1/n10000/shuffled` | 62.2 µs | 60 µs | 22.9 µs |  |  | **2.7x slower** |
| `bspline2/n10000/sorted` | 482 µs | 489 µs | 768 µs |  |  | 1.6x faster |
| `bspline2/n10000/shuffled` | 476 µs | 517 µs | 752 µs |  |  | 1.6x faster |
| `bspline3/n10000/sorted` | 572 µs | 605 µs | 934 µs |  |  | 1.6x faster |
| `bspline3/n10000/shuffled` | 624 µs | 618 µs | 939 µs |  |  | 1.5x faster |
| `cubic_spline/n10000/sorted` | 231 µs | 270 µs | 702 µs | 334 µs | 358 µs | 1.4x faster |
| `cubic_spline/n10000/shuffled` | 246 µs | 258 µs | 714 µs | 252 µs | 346 µs | 1x faster |
| `pchip/n10000/sorted` | 82.8 µs | 101 µs | 239 µs | 154 µs |  | 1.9x faster |
| `pchip/n10000/shuffled` | 87.8 µs | 200 µs | 240 µs | 155 µs |  | 1.8x faster |
| `akima/n10000/sorted` | 115 µs | 84.4 µs | 600 µs | 121 µs |  | 1.1x faster |
| `akima/n10000/shuffled` | 78.3 µs | 69.6 µs | 625 µs | 137 µs |  | 1.7x faster |
| `makima/n10000/sorted` | 118 µs | 196 µs | 696 µs |  |  | 5.9x faster |
| `makima/n10000/shuffled` | 83.6 µs | 198 µs | 679 µs |  |  | 8.1x faster |
| `linear/n100/sorted-uneven` | 248 ns | 249 ns | 8.91 µs | 190 ns | 167 ns | **1.5x slower** |
| `linear/n100/shuffled-uneven` | 248 ns | 248 ns | 8.83 µs | 189 ns | 164 ns | **1.5x slower** |
| `cubic_spline/n100/sorted-uneven` | 1.58 µs | 1.58 µs | 96.6 µs | 1.52 µs | 1.87 µs | **1x slower** |
| `cubic_spline/n100/shuffled-uneven` | 1.58 µs | 1.58 µs | 86.9 µs | 1.5 µs | 1.89 µs | **1.1x slower** |
| `pchip/n100/sorted-uneven` | 887 ns | 889 ns | 103 µs | 535 ns |  | **1.7x slower** |
| `pchip/n100/shuffled-uneven` | 887 ns | 886 ns | 102 µs | 540 ns |  | **1.6x slower** |
| `linear/n10000/sorted-uneven` | 88.9 µs | 86.4 µs | 47 µs | 56.2 µs | 31.4 µs | **2.8x slower** |
| `linear/n10000/shuffled-uneven` | 88.8 µs | 85.1 µs | 46.8 µs | 57.2 µs | 14.6 µs | **6.1x slower** |
| `cubic_spline/n10000/sorted-uneven` | 248 µs | 322 µs | 705 µs | 323 µs | 347 µs | 1.3x faster |
| `cubic_spline/n10000/shuffled-uneven` | 250 µs | 253 µs | 717 µs | 289 µs | 352 µs | 1.2x faster |
| `pchip/n10000/sorted-uneven` | 124 µs | 94.5 µs | 238 µs | 156 µs |  | 1.3x faster |
| `pchip/n10000/shuffled-uneven` | 91.1 µs | 202 µs | 238 µs | 152 µs |  | 1.7x faster |

## `PPoly`

A cubic of 55 intervals, evaluated at points in no order.

### Evaluation

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `intervals55/eval10` | 73.5 ns | 74.2 ns | 3.5 µs | 48x faster |
| `intervals55/eval100` | 659 ns | 725 ns | 4.42 µs | 6.7x faster |
| `intervals55/eval1000` | 6.53 µs | 6.65 µs | 15.5 µs | 2.4x faster |
| `intervals55/eval10000` | 72.9 µs | 68.9 µs | 332 µs | 4.6x faster |

### Construction

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `intervals55/eval10` | 163 ns | 163 ns | 5.07 µs | 31x faster |
| `intervals55/eval100` | 166 ns | 160 ns | 5.22 µs | 32x faster |
| `intervals55/eval1000` | 159 ns | 160 ns | 5.13 µs | 32x faster |
| `intervals55/eval10000` | 158 ns | 165 ns | 5.25 µs | 33x faster |

## `BSpline`

A cubic B-spline of 55 intervals, evaluated at points in no order.

### Evaluation

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `intervals55/eval10` | 142 ns | 147 ns | 2.66 µs | 19x faster |
| `intervals55/eval100` | 1.41 µs | 1.5 µs | 5.76 µs | 4.1x faster |
| `intervals55/eval1000` | 14.2 µs | 14.4 µs | 37.6 µs | 2.7x faster |
| `intervals55/eval10000` | 141 µs | 51.6 µs | 361 µs | 2.6x faster |

### Construction

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `intervals55/eval10` | 162 ns | 162 ns | 9.13 µs | 56x faster |
| `intervals55/eval100` | 162 ns | 165 ns | 9.31 µs | 58x faster |
| `intervals55/eval1000` | 162 ns | 162 ns | 9.25 µs | 57x faster |
| `intervals55/eval10000` | 162 ns | 161 ns | 9 µs | 56x faster |

## `BSplineCurve`

A cubic curve through `n` points of the plane, evaluated at `4 n` parameters in order.

### Evaluation

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `degree3/n100` | 4.11 µs | 4.22 µs | 15.9 µs | 3.9x faster |
| `degree3/n10000` | 474 µs | 158 µs | 1.31 ms | 2.8x faster |

### Construction

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `degree3/n100` | 4.25 µs | 4.23 µs | 45.9 µs | 11x faster |
| `degree3/n10000` | 622 µs | 675 µs | 1.34 ms | 2.1x faster |

## Regular grid

Points in no order inside the grid. `maxN`: the first axis has `N` points and each axis half the one before (SciPy's suite); `eachN`: every axis has `N` points; `flipped`: every axis given from its high end.

### Evaluation

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | ndarray-interp | interpn | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/2d/max10/eval10/ascending` | 135 ns | 136 ns | 135 ns | 19.6 µs | 791 ns | 149 ns | 94.6 ns | **1.4x slower** |
| `linear/2d/max10/eval10/flipped` | 135 ns | 137 ns | 138 ns | 19.7 µs |  |  |  | 1.5e+02x faster |
| `linear/2d/max10/eval100/ascending` | 1.37 µs | 1.32 µs | 1.25 µs | 21.8 µs | 8.55 µs | 1.36 µs | 927 ns | **1.5x slower** |
| `linear/2d/max10/eval100/flipped` | 1.25 µs | 1.34 µs | 1.25 µs | 22.3 µs |  |  |  | 18x faster |
| `linear/2d/max10/eval1000/ascending` | 13.2 µs | 12.6 µs | 12.5 µs | 40.6 µs | 81.9 µs | 14.4 µs | 9.25 µs | **1.4x slower** |
| `linear/2d/max10/eval1000/flipped` | 12.4 µs | 13.8 µs | 12.4 µs | 40.5 µs |  |  |  | 3.3x faster |
| `linear/2d/max10/eval10000/ascending` | 125 µs | 47.7 µs | 125 µs | 361 µs | 858 µs | 231 µs | 92.3 µs | **1.4x slower** |
| `linear/2d/max10/eval10000/flipped` | 124 µs | 46.8 µs | 125 µs | 359 µs |  |  |  | 2.9x faster |
| `linear/2d/max40/eval10/ascending` | 156 ns | 158 ns | 157 ns | 19.5 µs | 968 ns | 154 ns | 134 ns | **1.2x slower** |
| `linear/2d/max40/eval10/flipped` | 160 ns | 159 ns | 157 ns | 19.6 µs |  |  |  | 1.2e+02x faster |
| `linear/2d/max40/eval100/ascending` | 1.51 µs | 1.58 µs | 1.48 µs | 22.5 µs | 10.2 µs | 1.49 µs | 1.32 µs | **1.1x slower** |
| `linear/2d/max40/eval100/flipped` | 1.49 µs | 1.73 µs | 1.51 µs | 22.3 µs |  |  |  | 15x faster |
| `linear/2d/max40/eval1000/ascending` | 15.1 µs | 15.3 µs | 15.1 µs | 46.4 µs | 99.4 µs | 15.4 µs | 13.2 µs | **1.1x slower** |
| `linear/2d/max40/eval1000/flipped` | 15.1 µs | 15.2 µs | 14.8 µs | 46.3 µs |  |  |  | 3.1x faster |
| `linear/2d/max40/eval10000/ascending` | 151 µs | 64.3 µs | 150 µs | 562 µs | 1.03 ms | 222 µs | 134 µs | **1.1x slower** |
| `linear/2d/max40/eval10000/flipped` | 152 µs | 55.8 µs | 155 µs | 570 µs |  |  |  | 3.8x faster |
| `linear/2d/max200/eval10/ascending` | 192 ns | 193 ns | 192 ns | 19.7 µs | 1.96 µs | 187 ns | 176 ns | **1.1x slower** |
| `linear/2d/max200/eval10/flipped` | 192 ns | 194 ns | 192 ns | 19.9 µs |  |  |  | 1e+02x faster |
| `linear/2d/max200/eval100/ascending` | 1.86 µs | 1.98 µs | 1.87 µs | 23 µs | 19.3 µs | 1.96 µs | 1.72 µs | **1.1x slower** |
| `linear/2d/max200/eval100/flipped` | 1.88 µs | 2.14 µs | 1.87 µs | 22.9 µs |  |  |  | 12x faster |
| `linear/2d/max200/eval1000/ascending` | 18.5 µs | 18.7 µs | 18.7 µs | 72.8 µs | 202 µs | 23.4 µs | 17.2 µs | **1.1x slower** |
| `linear/2d/max200/eval1000/flipped` | 18.4 µs | 18.9 µs | 18.7 µs | 68.8 µs |  |  |  | 3.7x faster |
| `linear/2d/max200/eval10000/ascending` | 188 µs | 74.1 µs | 188 µs | 820 µs | 2.09 ms | 453 µs | 189 µs | 1x faster |
| `linear/2d/max200/eval10000/flipped` | 188 µs | 74.5 µs | 189 µs | 812 µs |  |  |  | 4.3x faster |
| `linear/3d/max10/eval10/ascending` | 195 ns | 201 ns | 196 ns | 56.9 µs | 3.07 µs |  | 141 ns | **1.4x slower** |
| `linear/3d/max10/eval10/flipped` | 200 ns | 224 ns | 200 ns | 56.9 µs |  |  |  | 2.8e+02x faster |
| `linear/3d/max10/eval100/ascending` | 1.91 µs | 1.95 µs | 1.89 µs | 62.3 µs | 31.4 µs |  | 1.4 µs | **1.4x slower** |
| `linear/3d/max10/eval100/flipped` | 1.88 µs | 1.96 µs | 1.9 µs | 63.1 µs |  |  |  | 34x faster |
| `linear/3d/max10/eval1000/ascending` | 19 µs | 18.9 µs | 18.6 µs | 119 µs | 320 µs |  | 13.9 µs | **1.4x slower** |
| `linear/3d/max10/eval1000/flipped` | 18.8 µs | 18.9 µs | 18.9 µs | 123 µs |  |  |  | 6.5x faster |
| `linear/3d/max10/eval10000/ascending` | 190 µs | 74.9 µs | 189 µs | 1.07 ms | 3.16 ms |  | 142 µs | **1.3x slower** |
| `linear/3d/max10/eval10000/flipped` | 191 µs | 66.5 µs | 190 µs | 1.07 ms |  |  |  | 5.6x faster |
| `linear/3d/max40/eval10/ascending` | 217 ns | 225 ns | 220 ns | 55.7 µs | 3.27 µs |  | 199 ns | **1.1x slower** |
| `linear/3d/max40/eval10/flipped` | 219 ns | 246 ns | 219 ns | 56 µs |  |  |  | 2.6e+02x faster |
| `linear/3d/max40/eval100/ascending` | 2.13 µs | 2.2 µs | 2.1 µs | 62.1 µs | 33.8 µs |  | 1.99 µs | **1.1x slower** |
| `linear/3d/max40/eval100/flipped` | 2.13 µs | 2.19 µs | 2.09 µs | 63.3 µs |  |  |  | 30x faster |
| `linear/3d/max40/eval1000/ascending` | 21.3 µs | 22 µs | 20.8 µs | 137 µs | 357 µs |  | 19.6 µs | **1.1x slower** |
| `linear/3d/max40/eval1000/flipped` | 21.3 µs | 21.5 µs | 21.3 µs | 136 µs |  |  |  | 6.4x faster |
| `linear/3d/max40/eval10000/ascending` | 213 µs | 70.7 µs | 233 µs | 1.47 ms | 3.35 ms |  | 200 µs | **1.1x slower** |
| `linear/3d/max40/eval10000/flipped` | 208 µs | 69.9 µs | 213 µs | 1.48 ms |  |  |  | 7.1x faster |
| `linear/3d/max200/eval10/ascending` | 275 ns | 280 ns | 274 ns | 56.4 µs | 4.42 µs |  | 286 ns | 1x faster |
| `linear/3d/max200/eval10/flipped` | 274 ns | 278 ns | 273 ns | 57.3 µs |  |  |  | 2.1e+02x faster |
| `linear/3d/max200/eval100/ascending` | 2.95 µs | 2.77 µs | 2.7 µs | 64.9 µs | 46 µs |  | 2.64 µs | **1.1x slower** |
| `linear/3d/max200/eval100/flipped` | 2.95 µs | 2.78 µs | 2.96 µs | 71 µs |  |  |  | 24x faster |
| `linear/3d/max200/eval1000/ascending` | 27.5 µs | 27.3 µs | 27 µs | 184 µs | 462 µs |  | 28.6 µs | 1x faster |
| `linear/3d/max200/eval1000/flipped` | 26.9 µs | 27.2 µs | 26.4 µs | 187 µs |  |  |  | 6.9x faster |
| `linear/3d/max200/eval10000/ascending` | 269 µs | 84.4 µs | 269 µs | 1.82 ms | 4.69 ms |  | 289 µs | 1.1x faster |
| `linear/3d/max200/eval10000/flipped` | 269 µs | 91.9 µs | 272 µs | 1.8 ms |  |  |  | 6.7x faster |
| `linear/4d/max10/eval10/ascending` | 256 ns | 261 ns | 259 ns | 93.9 µs |  |  |  | 3.7e+02x faster |
| `linear/4d/max10/eval10/flipped` | 260 ns | 265 ns | 261 ns | 95.7 µs |  |  |  | 3.7e+02x faster |
| `linear/4d/max10/eval100/ascending` | 2.4 µs | 2.49 µs | 2.41 µs | 104 µs |  |  |  | 43x faster |
| `linear/4d/max10/eval100/flipped` | 2.41 µs | 2.48 µs | 2.41 µs | 104 µs |  |  |  | 43x faster |
| `linear/4d/max10/eval1000/ascending` | 24.4 µs | 24.6 µs | 24.3 µs | 211 µs |  |  |  | 8.6x faster |
| `linear/4d/max10/eval1000/flipped` | 24.1 µs | 24.5 µs | 25.1 µs | 206 µs |  |  |  | 8.6x faster |
| `linear/4d/max10/eval10000/ascending` | 241 µs | 82.2 µs | 255 µs | 1.79 ms |  |  |  | 7.4x faster |
| `linear/4d/max10/eval10000/flipped` | 243 µs | 80.2 µs | 244 µs | 1.82 ms |  |  |  | 7.5x faster |
| `linear/4d/max40/eval10/ascending` | 303 ns | 311 ns | 308 ns | 93.8 µs | 5.74 µs |  | 296 ns | **1x slower** |
| `linear/4d/max40/eval10/flipped` | 302 ns | 341 ns | 337 ns | 94.6 µs |  |  |  | 3.1e+02x faster |
| `linear/4d/max40/eval100/ascending` | 2.94 µs | 3.01 µs | 2.93 µs | 107 µs | 57 µs |  | 2.68 µs | **1.1x slower** |
| `linear/4d/max40/eval100/flipped` | 2.95 µs | 3.05 µs | 2.94 µs | 107 µs |  |  |  | 36x faster |
| `linear/4d/max40/eval1000/ascending` | 29.6 µs | 29.6 µs | 29.4 µs | 243 µs | 571 µs |  | 26.5 µs | **1.1x slower** |
| `linear/4d/max40/eval1000/flipped` | 30.4 µs | 29.4 µs | 29.3 µs | 246 µs |  |  |  | 8.1x faster |
| `linear/4d/max40/eval10000/ascending` | 290 µs | 103 µs | 320 µs | 2.23 ms | 5.88 ms |  | 270 µs | **1.1x slower** |
| `linear/4d/max40/eval10000/flipped` | 293 µs | 91.8 µs | 293 µs | 2.23 ms |  |  |  | 7.6x faster |
| `linear/4d/max200/eval10/ascending` | 376 ns | 382 ns | 386 ns | 94.6 µs | 7.49 µs |  | 356 ns | **1.1x slower** |
| `linear/4d/max200/eval10/flipped` | 379 ns | 383 ns | 379 ns | 95.9 µs |  |  |  | 2.5e+02x faster |
| `linear/4d/max200/eval100/ascending` | 3.73 µs | 3.82 µs | 3.73 µs | 109 µs | 71.4 µs |  | 3.62 µs | **1x slower** |
| `linear/4d/max200/eval100/flipped` | 3.73 µs | 3.82 µs | 3.65 µs | 108 µs |  |  |  | 29x faster |
| `linear/4d/max200/eval1000/ascending` | 37.2 µs | 37.6 µs | 37.3 µs | 314 µs | 707 µs |  | 36.3 µs | **1x slower** |
| `linear/4d/max200/eval1000/flipped` | 37.2 µs | 38 µs | 37.4 µs | 303 µs |  |  |  | 8.1x faster |
| `linear/4d/max200/eval10000/ascending` | 371 µs | 113 µs | 375 µs | 2.74 ms | 7.34 ms |  | 364 µs | **1x slower** |
| `linear/4d/max200/eval10000/flipped` | 376 µs | 110 µs | 372 µs | 2.76 ms |  |  |  | 7.3x faster |
| `nearest/2d/each40/eval10000/ascending` | 150 µs | 54.1 µs | 150 µs | 637 µs | 1.45 ms |  | 133 µs | **1.1x slower** |
| `nearest/3d/each40/eval10000/ascending` | 212 µs | 74.2 µs | 233 µs | 933 µs | 1.7 ms |  | 208 µs | **1x slower** |
| `pchip/2d/each40/eval10000/ascending` | 499 µs | 131 µs | 499 µs | 1.19 s |  |  |  | 2.4e+03x faster |
| `pchip/3d/each40/eval10000/ascending` | 1.51 ms | 278 µs | 1.5 ms | 3.38 s |  |  |  | 2.2e+03x faster |
| `cubic/2d/each40/eval10000/ascending` | 301 µs | 96 µs | 301 µs | 1.29 ms |  |  | 409 µs | 1.4x faster |
| `cubic/3d/each40/eval10000/ascending` | 659 µs | 173 µs | 651 µs | 3.06 ms |  |  | 807 µs | 1.2x faster |
| `cubic/2d/each10/eval1000/ascending` | 26.5 µs | 26.6 µs | 26.3 µs | 130 µs |  |  | 35.8 µs | 1.3x faster |
| `cubic/2d/each40/eval1000/ascending` | 30 µs | 30.2 µs | 30.1 µs | 142 µs |  |  | 40.8 µs | 1.4x faster |
| `cubic/2d/each100/eval1000/ascending` | 31.7 µs | 35.5 µs | 35.2 µs | 172 µs |  |  | 45.3 µs | 1.4x faster |
| `cubic/2d/each200/eval1000/ascending` | 34 µs | 34.2 µs | 34.5 µs | 198 µs |  |  | 45.8 µs | 1.3x faster |
| `cubic/2d/each400/eval1000/ascending` | 38.1 µs | 37.9 µs | 38.2 µs | 271 µs |  |  | 48.4 µs | 1.3x faster |
| `quintic/2d/each10/eval1000/ascending` | 44.1 µs | 44.3 µs | 44.1 µs | 227 µs |  |  |  | 5.1x faster |
| `quintic/2d/each40/eval1000/ascending` | 46.9 µs | 48.1 µs | 49.8 µs | 238 µs |  |  |  | 5.1x faster |
| `linear/2d/each400/eval10000/ascending` | 220 µs | 74.1 µs | 219 µs | 949 µs | 4.07 ms | 532 µs | 219 µs | **1x slower** |
| `linear/2d/each800/eval10000/ascending` | 239 µs | 82.8 µs | 238 µs | 1.07 ms | 7.45 ms | 571 µs | 242 µs | 1x faster |
| `nearest/2d/auckland-ba32-0401/heldout/ascending` | 346 ms | 58.8 ms | 349 ms | 1.31 s | 159 s |  | 723 ms | 2.1x faster |
| `linear/2d/auckland-ba32-0401/heldout/ascending` | 365 ms | 65 ms | 365 ms | 1.08 s | 67 s | 344 ms | 754 ms | **1.1x slower** |
| `cubic/2d/auckland-ba32-0401/heldout/ascending` | 724 ms | 132 ms | 712 ms | 32.4 s |  |  | 1.6 s | 2.2x faster |
| `quintic/2d/auckland-ba32-0401/heldout/ascending` | 1.18 s | 203 ms | 1.18 s | 34.2 s |  |  |  | 29x faster |
| `pchip/2d/auckland-ba32-0401/heldout/ascending` | 1.02 s | 170 ms | 1.03 s |  |  |  |  |  |

### Construction

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | ndarray-interp | interpn | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/2d/max10/eval10/ascending` | 134 ns | 123 ns | 117 ns | 5.78 µs | 132 ns | 78.1 ns | 0 | **2e+02x slower** |
| `linear/2d/max10/eval10/flipped` | 160 ns | 163 ns | 159 ns | 11.8 µs |  |  |  | 74x faster |
| `linear/2d/max10/eval100/ascending` | 134 ns | 129 ns | 118 ns | 5.78 µs | 135 ns | 85.1 ns | 0 | **2e+02x slower** |
| `linear/2d/max10/eval100/flipped` | 175 ns | 163 ns | 155 ns | 11.6 µs |  |  |  | 66x faster |
| `linear/2d/max10/eval1000/ascending` | 135 ns | 124 ns | 121 ns | 5.81 µs | 133 ns | 78 ns | 0 | **2e+02x slower** |
| `linear/2d/max10/eval1000/flipped` | 161 ns | 161 ns | 156 ns | 11.6 µs |  |  |  | 72x faster |
| `linear/2d/max10/eval10000/ascending` | 122 ns | 124 ns | 115 ns | 5.8 µs | 132 ns | 78.3 ns | 0 | **1.8e+02x slower** |
| `linear/2d/max10/eval10000/flipped` | 160 ns | 159 ns | 153 ns | 11.6 µs |  |  |  | 72x faster |
| `linear/2d/max40/eval10/ascending` | 397 ns | 387 ns | 658 ns | 5.82 µs | 205 ns | 152 ns | 0 | **5.8e+02x slower** |
| `linear/2d/max40/eval10/flipped` | 706 ns | 645 ns | 1.21 µs | 11.5 µs |  |  |  | 16x faster |
| `linear/2d/max40/eval100/ascending` | 360 ns | 354 ns | 605 ns | 5.75 µs | 200 ns | 154 ns | 0 | **5.3e+02x slower** |
| `linear/2d/max40/eval100/flipped` | 645 ns | 657 ns | 1.17 µs | 11.5 µs |  |  |  | 18x faster |
| `linear/2d/max40/eval1000/ascending` | 359 ns | 358 ns | 300 ns | 5.84 µs | 201 ns | 153 ns | 0 | **5.3e+02x slower** |
| `linear/2d/max40/eval1000/flipped` | 648 ns | 638 ns | 1.16 µs | 11.4 µs |  |  |  | 18x faster |
| `linear/2d/max40/eval10000/ascending` | 364 ns | 360 ns | 658 ns | 5.73 µs | 202 ns | 153 ns | 0 | **5.3e+02x slower** |
| `linear/2d/max40/eval10000/flipped` | 644 ns | 642 ns | 1.27 µs | 11.7 µs |  |  |  | 18x faster |
| `linear/2d/max200/eval10/ascending` | 6.83 µs | 6.92 µs | 14.9 µs | 6 µs | 3.11 µs | 472 ns | 0 | **1e+04x slower** |
| `linear/2d/max200/eval10/flipped` | 14.7 µs | 14.5 µs | 25.7 µs | 11.9 µs |  |  |  | **1.2x slower** |
| `linear/2d/max200/eval100/ascending` | 6.84 µs | 6.91 µs | 14.6 µs | 5.77 µs | 2.93 µs | 472 ns | 0 | **1e+04x slower** |
| `linear/2d/max200/eval100/flipped` | 14.4 µs | 14.8 µs | 25.4 µs | 12.6 µs |  |  |  | **1.1x slower** |
| `linear/2d/max200/eval1000/ascending` | 6.78 µs | 6.92 µs | 15.2 µs | 5.82 µs | 2.96 µs | 471 ns | 0 | **9.9e+03x slower** |
| `linear/2d/max200/eval1000/flipped` | 14.3 µs | 14.3 µs | 25.9 µs | 11.7 µs |  |  |  | **1.2x slower** |
| `linear/2d/max200/eval10000/ascending` | 47.1 µs | 46.7 µs | 15.7 µs | 5.75 µs | 51 µs | 463 ns | 0 | **6.9e+04x slower** |
| `linear/2d/max200/eval10000/flipped` | 53 µs | 51.8 µs | 28.4 µs | 11.9 µs |  |  |  | **4.5x slower** |
| `linear/3d/max10/eval10/ascending` | 180 ns | 179 ns | 184 ns | 7.52 µs | 177 ns |  | 0 | **2e+02x slower** |
| `linear/3d/max10/eval10/flipped` | 306 ns | 308 ns | 345 ns | 15.7 µs |  |  |  | 51x faster |
| `linear/3d/max10/eval100/ascending` | 180 ns | 179 ns | 185 ns | 7.55 µs | 167 ns |  | 0 | **2e+02x slower** |
| `linear/3d/max10/eval100/flipped` | 334 ns | 304 ns | 321 ns | 15.9 µs |  |  |  | 48x faster |
| `linear/3d/max10/eval1000/ascending` | 183 ns | 181 ns | 189 ns | 7.6 µs | 169 ns |  | 0 | **2.1e+02x slower** |
| `linear/3d/max10/eval1000/flipped` | 309 ns | 306 ns | 318 ns | 17 µs |  |  |  | 55x faster |
| `linear/3d/max10/eval10000/ascending` | 180 ns | 179 ns | 184 ns | 7.53 µs | 164 ns |  | 0 | **2.1e+02x slower** |
| `linear/3d/max10/eval10000/flipped` | 315 ns | 304 ns | 317 ns | 15.5 µs |  |  |  | 49x faster |
| `linear/3d/max40/eval10/ascending` | 2.96 µs | 15 µs | 5.92 µs | 7.55 µs | 11 µs |  | 0 | **3.4e+03x slower** |
| `linear/3d/max40/eval10/flipped` | 8.15 µs | 8.95 µs | 12.3 µs | 15.5 µs |  |  |  | 1.9x faster |
| `linear/3d/max40/eval100/ascending` | 3.04 µs | 15.1 µs | 5.43 µs | 7.51 µs | 6.77 µs |  | 0 | **3.5e+03x slower** |
| `linear/3d/max40/eval100/flipped` | 8.15 µs | 21.6 µs | 12.3 µs | 15.6 µs |  |  |  | 1.9x faster |
| `linear/3d/max40/eval1000/ascending` | 16.6 µs | 15.5 µs | 5.56 µs | 7.56 µs | 5.96 µs |  | 0 | **1.7e+04x slower** |
| `linear/3d/max40/eval1000/flipped` | 28.9 µs | 26.4 µs | 12.6 µs | 15.7 µs |  |  |  | **1.8x slower** |
| `linear/3d/max40/eval10000/ascending` | 13.9 µs | 15 µs | 5.56 µs | 7.43 µs | 5.88 µs |  | 0 | **1.6e+04x slower** |
| `linear/3d/max40/eval10000/flipped` | 27.3 µs | 27.4 µs | 12.4 µs | 15.6 µs |  |  |  | **1.7x slower** |
| `linear/3d/max200/eval10/ascending` | 1.54 ms | 2.39 ms | 595 µs | 7.56 µs | 1.2 ms |  | 0 | **1.6e+06x slower** |
| `linear/3d/max200/eval10/flipped` | 3.02 ms | 3.62 ms | 1.85 ms | 16.1 µs |  |  |  | **1.9e+02x slower** |
| `linear/3d/max200/eval100/ascending` | 1.75 ms | 1.69 ms | 625 µs | 7.6 µs | 1.19 ms |  | 0 | **2e+06x slower** |
| `linear/3d/max200/eval100/flipped` | 3.12 ms | 3.27 ms | 1.96 ms | 17.5 µs |  |  |  | **1.8e+02x slower** |
| `linear/3d/max200/eval1000/ascending` | 1.77 ms | 1.75 ms | 634 µs | 8.35 µs | 1.27 ms |  | 0 | **1.8e+06x slower** |
| `linear/3d/max200/eval1000/flipped` | 3.22 ms | 3.36 ms | 1.86 ms | 15.8 µs |  |  |  | **2e+02x slower** |
| `linear/3d/max200/eval10000/ascending` | 1.63 ms | 1.77 ms | 597 µs | 7.54 µs | 1.2 ms |  | 0 | **1.9e+06x slower** |
| `linear/3d/max200/eval10000/flipped` | 3.12 ms | 3.38 ms | 1.95 ms | 15.8 µs |  |  |  | **2e+02x slower** |
| `linear/4d/max10/eval10/ascending` | 216 ns | 216 ns | 218 ns | 9.73 µs |  |  |  | 45x faster |
| `linear/4d/max10/eval10/flipped` | 351 ns | 339 ns | 329 ns | 17.9 µs |  |  |  | 51x faster |
| `linear/4d/max10/eval100/ascending` | 211 ns | 215 ns | 205 ns | 9.72 µs |  |  |  | 46x faster |
| `linear/4d/max10/eval100/flipped` | 330 ns | 331 ns | 332 ns | 17.8 µs |  |  |  | 54x faster |
| `linear/4d/max10/eval1000/ascending` | 215 ns | 214 ns | 205 ns | 9.63 µs |  |  |  | 45x faster |
| `linear/4d/max10/eval1000/flipped` | 334 ns | 336 ns | 334 ns | 17.7 µs |  |  |  | 53x faster |
| `linear/4d/max10/eval10000/ascending` | 216 ns | 218 ns | 208 ns | 9.65 µs |  |  |  | 45x faster |
| `linear/4d/max10/eval10000/flipped` | 332 ns | 330 ns | 333 ns | 17.9 µs |  |  |  | 54x faster |
| `linear/4d/max40/eval10/ascending` | 13.1 µs | 14.4 µs | 34.8 µs | 9.32 µs | 5.4 µs |  | 1.17 ns | **1.1e+04x slower** |
| `linear/4d/max40/eval10/flipped` | 52 µs | 57.1 µs | 64.2 µs | 19.6 µs |  |  |  | **2.6x slower** |
| `linear/4d/max40/eval100/ascending` | 13 µs | 13.3 µs | 31.3 µs | 9.32 µs | 5.48 µs |  | 1.17 ns | **1.1e+04x slower** |
| `linear/4d/max40/eval100/flipped` | 51.9 µs | 52.6 µs | 63.9 µs | 19.5 µs |  |  |  | **2.7x slower** |
| `linear/4d/max40/eval1000/ascending` | 13.2 µs | 13.5 µs | 25.8 µs | 9.28 µs | 5.58 µs |  | 1.17 ns | **1.1e+04x slower** |
| `linear/4d/max40/eval1000/flipped` | 52.7 µs | 52.4 µs | 63.2 µs | 19.9 µs |  |  |  | **2.6x slower** |
| `linear/4d/max40/eval10000/ascending` | 13.1 µs | 94.7 µs | 29.6 µs | 9.31 µs | 5.52 µs |  | 1.17 ns | **1.1e+04x slower** |
| `linear/4d/max40/eval10000/flipped` | 54.2 µs | 131 µs | 68 µs | 19.3 µs |  |  |  | **2.8x slower** |
| `linear/4d/max200/eval10/ascending` | 43.2 ms | 47.4 ms | 15 ms | 9.24 µs | 35.9 ms |  | 1.17 ns | **3.7e+07x slower** |
| `linear/4d/max200/eval10/flipped` | 97.6 ms | 114 ms | 64.9 ms | 20.1 µs |  |  |  | **4.8e+03x slower** |
| `linear/4d/max200/eval100/ascending` | 45.7 ms | 47.6 ms | 14.7 ms | 9.44 µs | 37.7 ms |  | 1.17 ns | **3.9e+07x slower** |
| `linear/4d/max200/eval100/flipped` | 97.5 ms | 102 ms | 64.6 ms | 20.1 µs |  |  |  | **4.8e+03x slower** |
| `linear/4d/max200/eval1000/ascending` | 47.9 ms | 47.2 ms | 15.4 ms | 9.27 µs | 38.9 ms |  | 1.07 ns | **4.5e+07x slower** |
| `linear/4d/max200/eval1000/flipped` | 99.1 ms | 102 ms | 67.4 ms | 20.3 µs |  |  |  | **4.9e+03x slower** |
| `linear/4d/max200/eval10000/ascending` | 46.4 ms | 47.4 ms | 15.5 ms | 9.23 µs | 38.5 ms |  | 1.07 ns | **4.3e+07x slower** |
| `linear/4d/max200/eval10000/flipped` | 99.6 ms | 99.8 ms | 65.3 ms | 19.9 µs |  |  |  | **5e+03x slower** |
| `nearest/2d/each40/eval10000/ascending` | 633 ns | 584 ns | 1.13 µs | 5.75 µs | 288 ns |  | 0 | **9.3e+02x slower** |
| `nearest/3d/each40/eval10000/ascending` | 21 µs | 23 µs | 54.8 µs | 7.54 µs | 8.88 µs |  | 0 | **2.1e+04x slower** |
| `pchip/2d/each40/eval10000/ascending` | 6.47 µs | 7.09 µs | 6.25 µs | 6.53 µs |  |  |  | 1x faster |
| `pchip/3d/each40/eval10000/ascending` | 508 µs | 582 µs | 450 µs | 8.38 µs |  |  |  | **61x slower** |
| `cubic/2d/each40/eval10000/ascending` | 8.48 µs | 9.15 µs | 8.21 µs | 1.82 ms |  |  | 0 | **9.7e+03x slower** |
| `cubic/3d/each40/eval10000/ascending` | 389 µs | 392 µs | 286 µs | 138 ms |  |  | 1.17 ns | **3.3e+05x slower** |
| `cubic/2d/each10/eval1000/ascending` | 1.48 µs | 1.46 µs | 1.5 µs | 996 µs |  |  | 0 | **1.7e+03x slower** |
| `cubic/2d/each40/eval1000/ascending` | 8.46 µs | 8.47 µs | 8.16 µs | 1.8 ms |  |  | 0 | **1.1e+04x slower** |
| `cubic/2d/each100/eval1000/ascending` | 81.5 µs | 110 µs | 74.1 µs | 6.92 ms |  |  | 0 | **9.3e+04x slower** |
| `cubic/2d/each200/eval1000/ascending` | 172 µs | 173 µs | 268 µs | 33.4 ms |  |  | 0 | **2e+05x slower** |
| `cubic/2d/each400/eval1000/ascending` | 1.38 ms | 1.34 ms | 1.18 ms | 138 ms |  |  | 0 | **1.6e+06x slower** |
| `quintic/2d/each10/eval1000/ascending` | 1.94 µs | 1.94 µs | 1.97 µs | 3.37 ms |  |  |  | 1.7e+03x faster |
| `quintic/2d/each40/eval1000/ascending` | 12.5 µs | 12.4 µs | 12.6 µs | 7.23 ms |  |  |  | 5.8e+02x faster |
| `linear/2d/each400/eval10000/ascending` | 311 µs | 341 µs | 97 µs | 5.88 µs | 261 µs | 1.12 µs | 0 | **4.5e+05x slower** |
| `linear/2d/each800/eval10000/ascending` | 1.09 ms | 1.04 ms | 393 µs | 6.32 µs | 759 µs | 2.22 µs | 0 | **1.6e+06x slower** |
| `nearest/2d/auckland-ba32-0401/heldout/ascending` | 15.3 ms | 16.9 ms | 5.16 ms | 7.04 µs | 12.4 ms |  | 0 | **2.2e+07x slower** |
| `linear/2d/auckland-ba32-0401/heldout/ascending` | 15.9 ms | 15.5 ms | 5.19 ms | 7.53 µs | 12.6 ms | 9.05 µs | 0 | **2.3e+07x slower** |
| `cubic/2d/auckland-ba32-0401/heldout/ascending` | 79.5 ms | 76.5 ms | 65.9 ms | 13.9 s |  |  | 0 | **9e+07x slower** |
| `quintic/2d/auckland-ba32-0401/heldout/ascending` | 93 ms | 90.9 ms | 81.5 ms | 21.4 s |  |  |  | 2.3e+02x faster |
| `pchip/2d/auckland-ba32-0401/heldout/ascending` | 74.8 ms | 74.1 ms | 62 ms |  |  |  |  |  |

## `NearestNDInterpolator`

`n` points of the unit square, evaluated at the pixels of a square image.

### Evaluation

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: |
| `n10/image50` | 33.3 µs | 33.7 µs | 112 µs | 432 µs | 13x faster |
| `n10/image100` | 125 µs | 44.2 µs | 445 µs | 1.73 ms | 14x faster |
| `n10/image500` | 3.01 ms | 524 µs | 11.2 ms | 43.7 ms | 15x faster |
| `n100/image50` | 115 µs | 95.8 µs | 187 µs | 567 µs | 4.9x faster |
| `n100/image100` | 424 µs | 120 µs | 720 µs | 2.22 ms | 5.2x faster |
| `n100/image500` | 9.16 ms | 1.56 ms | 16.9 ms | 53 ms | 5.8x faster |
| `n1000/image50` | 324 µs | 107 µs | 323 µs | 789 µs | 2.4x faster |
| `n1000/image100` | 1.18 ms | 276 µs | 1.17 ms | 3.02 ms | 2.6x faster |
| `n1000/image500` | 23.4 ms | 4.03 ms | 23.4 ms | 67.7 ms | 2.9x faster |
| `n10000/image50` | 330 µs | 106 µs | 328 µs | 974 µs | 2.9x faster |
| `n10000/image100` | 1.36 ms | 278 µs | 1.36 ms | 3.74 ms | 2.7x faster |
| `n10000/image500` | 30.5 ms | 5.24 ms | 30.2 ms | 84.5 ms | 2.8x faster |
| `n100000/image50` | 444 µs | 113 µs | 450 µs | 1.33 ms | 3x faster |
| `n100000/image100` | 1.76 ms | 345 µs | 1.74 ms | 5.22 ms | 3x faster |
| `n100000/image500` | 42.9 ms | 7.04 ms | 43.3 ms | 122 ms | 2.8x faster |
| `ghcn-tmax-20250101/heldout` | 159 µs | 62.6 µs | 172 µs | 511 µs | 3.2x faster |

### Construction

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: |
| `n10/image50` | 522 ns | 510 ns | 507 ns | 9.31 µs | 18x faster |
| `n10/image100` | 506 ns | 508 ns | 504 ns | 9.3 µs | 18x faster |
| `n10/image500` | 532 ns | 530 ns | 508 ns | 9.29 µs | 17x faster |
| `n100/image50` | 3.37 µs | 3.32 µs | 3.25 µs | 15 µs | 4.5x faster |
| `n100/image100` | 3.3 µs | 3.31 µs | 3.28 µs | 15 µs | 4.5x faster |
| `n100/image500` | 3.54 µs | 3.59 µs | 3.25 µs | 15 µs | 4.2x faster |
| `n1000/image50` | 51.8 µs | 59.4 µs | 54.1 µs | 119 µs | 2.3x faster |
| `n1000/image100` | 59.6 µs | 51.8 µs | 42.6 µs | 118 µs | 2x faster |
| `n1000/image500` | 50.4 µs | 51.1 µs | 48.4 µs | 117 µs | 2.3x faster |
| `n10000/image50` | 677 µs | 619 µs | 601 µs | 1.62 ms | 2.4x faster |
| `n10000/image100` | 627 µs | 676 µs | 652 µs | 1.62 ms | 2.6x faster |
| `n10000/image500` | 617 µs | 633 µs | 604 µs | 1.63 ms | 2.6x faster |
| `n100000/image50` | 7.76 ms | 7.8 ms | 7.68 ms | 23.8 ms | 3.1x faster |
| `n100000/image100` | 7.86 ms | 7.7 ms | 7.58 ms | 24.5 ms | 3.1x faster |
| `n100000/image500` | 8.55 ms | 7.83 ms | 7.63 ms | 24.4 ms | 2.9x faster |
| `ghcn-tmax-20250101/heldout` | 653 µs | 669 µs | 667 µs | 1.78 ms | 2.7x faster |

## `RBFInterpolator`

`n` points of a square, evaluated at 1000 points in no order (`scattered`) or at the 10 000 pixels of an image (`image`). `global`: all points; `k50`: the 50 nearest.

### Evaluation

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/global/n10/scattered` | 14.3 µs | 14.5 µs | 15.8 µs | 59.5 µs | 180 µs | 4.1x faster |
| `linear/k50/n10/scattered` | 14.3 µs | 14.2 µs | 14.4 µs | 1.37 ms |  | 96x faster |
| `thin_plate_spline/global/n10/scattered` | 70.9 µs | 70.9 µs | 69.4 µs | 111 µs | 289 µs | 1.6x faster |
| `thin_plate_spline/k50/n10/scattered` | 69.6 µs | 69.8 µs | 69.7 µs | 1.44 ms |  | 21x faster |
| `cubic/global/n10/scattered` | 16.8 µs | 16.5 µs | 16.9 µs | 62.1 µs | 196 µs | 3.7x faster |
| `cubic/k50/n10/scattered` | 17.2 µs | 16.7 µs | 16.8 µs | 1.34 ms |  | 78x faster |
| `quintic/global/n10/scattered` | 23.3 µs | 23.4 µs | 23.8 µs | 68.7 µs | 215 µs | 2.9x faster |
| `quintic/k50/n10/scattered` | 23.2 µs | 23.7 µs | 23.4 µs | 1.38 ms |  | 59x faster |
| `multiquadric/global/n10/scattered` | 27.2 µs | 27.9 µs | 27.8 µs | 52 µs | 236 µs | 1.9x faster |
| `multiquadric/k50/n10/scattered` | 27.8 µs | 27.9 µs | 27.8 µs | 1.35 ms |  | 49x faster |
| `inverse_multiquadric/global/n10/scattered` | 36.6 µs | 41.2 µs | 37.4 µs | 53.9 µs | 267 µs | 1.5x faster |
| `inverse_multiquadric/k50/n10/scattered` | 36.6 µs | 37.4 µs | 41.1 µs | 1.34 ms |  | 36x faster |
| `inverse_quadratic/global/n10/scattered` | 23.8 µs | 24 µs | 23.8 µs | 47.9 µs |  | 2x faster |
| `inverse_quadratic/k50/n10/scattered` | 23.8 µs | 24 µs | 23.8 µs | 1.31 ms |  | 55x faster |
| `gaussian/global/n10/scattered` | 67.2 µs | 65.9 µs | 67.1 µs | 87.8 µs | 356 µs | 1.3x faster |
| `gaussian/k50/n10/scattered` | 67.1 µs | 71.8 µs | 66.9 µs | 1.37 ms |  | 20x faster |
| `linear/global/n100/scattered` | 136 µs | 45.2 µs | 135 µs | 287 µs | 1.72 ms | 2.1x faster |
| `linear/k50/n100/scattered` | 13.7 ms | 2.73 ms | 12.7 ms | 18.3 ms |  | 1.3x faster |
| `thin_plate_spline/global/n100/scattered` | 697 µs | 170 µs | 698 µs | 868 µs | 2.9 ms | 1.2x faster |
| `thin_plate_spline/k50/n100/scattered` | 22.4 ms | 4.04 ms | 22.4 ms | 24 ms |  | 1.1x faster |
| `cubic/global/n100/scattered` | 139 µs | 50.9 µs | 139 µs | 288 µs | 1.9 ms | 2.1x faster |
| `cubic/k50/n100/scattered` | 14.8 ms | 2.89 ms | 14.6 ms | 19.7 ms |  | 1.3x faster |
| `quintic/global/n100/scattered` | 164 µs | 64.9 µs | 165 µs | 301 µs | 2.05 ms | 1.8x faster |
| `quintic/k50/n100/scattered` | 16.4 ms | 2.96 ms | 16.1 ms | 21.3 ms |  | 1.3x faster |
| `multiquadric/global/n100/scattered` | 269 µs | 79.1 µs | 269 µs | 295 µs | 2.26 ms | 1.1x faster |
| `multiquadric/k50/n100/scattered` | 14.9 ms | 2.8 ms | 15.1 ms | 18.2 ms |  | 1.2x faster |
| `inverse_multiquadric/global/n100/scattered` | 364 µs | 101 µs | 364 µs | 391 µs | 2.57 ms | 1.1x faster |
| `inverse_multiquadric/k50/n100/scattered` | 15.8 ms | 3.07 ms | 15.5 ms | 18.9 ms |  | 1.2x faster |
| `inverse_quadratic/global/n100/scattered` | 224 µs | 63.8 µs | 224 µs | 258 µs |  | 1.2x faster |
| `inverse_quadratic/k50/n100/scattered` | 13.6 ms | 2.65 ms | 13.8 ms | 17.7 ms |  | 1.3x faster |
| `gaussian/global/n100/scattered` | 619 µs | 152 µs | 615 µs | 677 µs | 3.47 ms | 1.1x faster |
| `gaussian/k50/n100/scattered` | 20 ms | 3.67 ms | 20.3 ms | 21.1 ms |  | 1.1x faster |
| `linear/global/n1000/scattered` | 1.34 ms | 325 µs | 1.34 ms | 3.69 ms | 17.3 ms | 2.8x faster |
| `linear/k50/n1000/scattered` | 13.7 ms | 2.57 ms | 13.8 ms | 33.3 ms |  | 2.4x faster |
| `thin_plate_spline/global/n1000/scattered` | 6.99 ms | 1.38 ms | 6.98 ms | 9.64 ms | 29.2 ms | 1.4x faster |
| `thin_plate_spline/k50/n1000/scattered` | 22.7 ms | 4.01 ms | 21.7 ms | 42.3 ms |  | 1.9x faster |
| `cubic/global/n1000/scattered` | 1.35 ms | 341 µs | 1.34 ms | 3.68 ms | 18.8 ms | 2.7x faster |
| `cubic/k50/n1000/scattered` | 16 ms | 2.94 ms | 15.6 ms | 35.9 ms |  | 2.2x faster |
| `quintic/global/n1000/scattered` | 1.51 ms | 344 µs | 1.5 ms | 3.77 ms | 20.7 ms | 2.5x faster |
| `quintic/k50/n1000/scattered` | 16.9 ms | 3.13 ms | 17.2 ms | 39.5 ms |  | 2.3x faster |
| `multiquadric/global/n1000/scattered` | 2.72 ms | 611 µs | 2.72 ms | 3.93 ms | 22.5 ms | 1.4x faster |
| `multiquadric/k50/n1000/scattered` | 16.7 ms | 2.96 ms | 16.4 ms | 36.6 ms |  | 2.2x faster |
| `inverse_multiquadric/global/n1000/scattered` | 3.7 ms | 798 µs | 3.64 ms | 4.87 ms | 25.9 ms | 1.3x faster |
| `inverse_multiquadric/k50/n1000/scattered` | 17.5 ms | 3.21 ms | 18 ms | 35.8 ms |  | 2x faster |
| `inverse_quadratic/global/n1000/scattered` | 2.25 ms | 499 µs | 2.24 ms | 3.6 ms |  | 1.6x faster |
| `inverse_quadratic/k50/n1000/scattered` | 15.8 ms | 2.85 ms | 15.9 ms | 33.4 ms |  | 2.1x faster |
| `gaussian/global/n1000/scattered` | 6.12 ms | 1.35 ms | 6.17 ms | 7.93 ms | 34.6 ms | 1.3x faster |
| `gaussian/k50/n1000/scattered` | 21.3 ms | 3.88 ms | 21.2 ms | 40.2 ms |  | 1.9x faster |
| `thin_plate_spline/global/n3000/scattered` | 20.7 ms | 4.19 ms | 21.3 ms | 28.6 ms | 93.3 ms | 1.4x faster |
| `thin_plate_spline/k50/n10000/scattered` | 22.5 ms | 4.03 ms | 22.2 ms | 42.8 ms |  | 1.9x faster |
| `thin_plate_spline/k50/n10000/image` | 224 ms | 36.5 ms | 225 ms | 435 ms |  | 1.9x faster |
| `thin_plate_spline/k50/n100000/scattered` | 23 ms | 4.23 ms | 23.2 ms | 46.3 ms |  | 2x faster |
| `thin_plate_spline/k50/n100000/image` | 228 ms | 37.4 ms | 229 ms | 440 ms |  | 1.9x faster |
| `linear/global/ghcn-tmax-20250101/heldout` | 17.8 ms | 3.69 ms | 18.1 ms | 52.9 ms | 251 ms | 3x faster |
| `linear/k50/ghcn-tmax-20250101/heldout` | 17.8 ms | 3.31 ms | 17.7 ms | 41.9 ms |  | 2.4x faster |
| `thin_plate_spline/global/ghcn-tmax-20250101/heldout` | 90.2 ms | 19.6 ms | 83.8 ms | 126 ms |  | 1.4x faster |
| `thin_plate_spline/k50/ghcn-tmax-20250101/heldout` | 28.4 ms | 5.13 ms | 28.3 ms | 52.9 ms |  | 1.9x faster |
| `cubic/global/ghcn-tmax-20250101/heldout` | 20.4 ms | 3.78 ms | 18.7 ms | 57.8 ms | 275 ms | 2.8x faster |
| `cubic/k50/ghcn-tmax-20250101/heldout` | 19.7 ms | 3.71 ms | 20.3 ms | 45.3 ms |  | 2.3x faster |
| `quintic/global/ghcn-tmax-20250101/heldout` | 20 ms | 4.07 ms | 22.3 ms | 50.5 ms | 301 ms | 2.5x faster |
| `quintic/k50/ghcn-tmax-20250101/heldout` | 21.7 ms | 3.97 ms | 23.2 ms | 50.5 ms |  | 2.3x faster |

### Construction

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/global/n10/scattered` | 1.23 µs | 1.34 µs | 1.26 µs | 10.7 µs | 5.37 µs | 4.4x faster |
| `linear/k50/n10/scattered` | 1.25 µs | 1.24 µs | 1.26 µs | 14.5 µs |  | 12x faster |
| `thin_plate_spline/global/n10/scattered` | 1.92 µs | 1.73 µs | 1.74 µs | 11.9 µs | 6.54 µs | 3.4x faster |
| `thin_plate_spline/k50/n10/scattered` | 1.74 µs | 1.73 µs | 1.72 µs | 15.4 µs |  | 8.8x faster |
| `cubic/global/n10/scattered` | 1.69 µs | 1.56 µs | 1.54 µs | 12.1 µs | 5.48 µs | 3.2x faster |
| `cubic/k50/n10/scattered` | 1.57 µs | 1.58 µs | 1.53 µs | 15.7 µs |  | 10x faster |
| `quintic/global/n10/scattered` | 2.08 µs | 2.08 µs | 2.06 µs | 13.9 µs | 5.74 µs | 2.8x faster |
| `quintic/k50/n10/scattered` | 2.05 µs | 2.01 µs | 2.05 µs | 16.9 µs |  | 8.2x faster |
| `multiquadric/global/n10/scattered` | 1.32 µs | 1.32 µs | 1.32 µs | 11.6 µs | 6.62 µs | 5x faster |
| `multiquadric/k50/n10/scattered` | 1.33 µs | 1.3 µs | 1.3 µs | 15 µs |  | 11x faster |
| `inverse_multiquadric/global/n10/scattered` | 1.3 µs | 1.31 µs | 1.36 µs | 10.8 µs | 6.88 µs | 5.3x faster |
| `inverse_multiquadric/k50/n10/scattered` | 1.32 µs | 1.34 µs | 1.33 µs | 14.6 µs |  | 11x faster |
| `inverse_quadratic/global/n10/scattered` | 1.24 µs | 1.25 µs | 1.26 µs | 10.7 µs |  | 8.6x faster |
| `inverse_quadratic/k50/n10/scattered` | 1.31 µs | 1.25 µs | 1.25 µs | 14.9 µs |  | 11x faster |
| `gaussian/global/n10/scattered` | 1.46 µs | 1.41 µs | 1.42 µs | 10.9 µs | 8.53 µs | 5.9x faster |
| `gaussian/k50/n10/scattered` | 1.43 µs | 1.42 µs | 1.42 µs | 14.9 µs |  | 10x faster |
| `linear/global/n100/scattered` | 49.5 µs | 51.2 µs | 49 µs | 75.2 µs | 722 µs | 1.5x faster |
| `linear/k50/n100/scattered` | 3.14 µs | 3.17 µs | 3.07 µs | 20.9 µs |  | 6.6x faster |
| `thin_plate_spline/global/n100/scattered` | 72.4 µs | 106 µs | 96 µs | 111 µs | 837 µs | 1.5x faster |
| `thin_plate_spline/k50/n100/scattered` | 3.48 µs | 3.23 µs | 3.59 µs | 21.6 µs |  | 6.2x faster |
| `cubic/global/n100/scattered` | 47.4 µs | 70.5 µs | 64.8 µs | 79.6 µs | 740 µs | 1.7x faster |
| `cubic/k50/n100/scattered` | 3.22 µs | 3.47 µs | 3.21 µs | 21.3 µs |  | 6.6x faster |
| `quintic/global/n100/scattered` | 52.4 µs | 78.1 µs | 72.5 µs | 91.3 µs | 752 µs | 1.7x faster |
| `quintic/k50/n100/scattered` | 3.21 µs | 3.32 µs | 3.17 µs | 22.4 µs |  | 7x faster |
| `multiquadric/global/n100/scattered` | 56.5 µs | 73 µs | 57.3 µs | 75.9 µs | 624 µs | 1.3x faster |
| `multiquadric/k50/n100/scattered` | 3.18 µs | 3.39 µs | 3.1 µs | 20.5 µs |  | 6.5x faster |
| `inverse_multiquadric/global/n100/scattered` | 57 µs | 78.9 µs | 61.2 µs | 81 µs | 804 µs | 1.4x faster |
| `inverse_multiquadric/k50/n100/scattered` | 3.18 µs | 3.12 µs | 3.06 µs | 20.7 µs |  | 6.5x faster |
| `inverse_quadratic/global/n100/scattered` | 48.8 µs | 79.1 µs | 54.2 µs | 72.7 µs |  | 1.5x faster |
| `inverse_quadratic/k50/n100/scattered` | 3.1 µs | 3.14 µs | 3.08 µs | 20.5 µs |  | 6.6x faster |
| `gaussian/global/n100/scattered` | 67.7 µs | 93 µs | 79.3 µs | 96.6 µs | 746 µs | 1.4x faster |
| `gaussian/k50/n100/scattered` | 3.13 µs | 3.39 µs | 3.11 µs | 20.6 µs |  | 6.6x faster |
| `linear/global/n1000/scattered` | 15.6 ms | 15.6 ms | 15.3 ms | 59 ms | 182 ms | 3.8x faster |
| `linear/k50/n1000/scattered` | 63.4 µs | 53.4 µs | 54.1 µs | 129 µs |  | 2x faster |
| `thin_plate_spline/global/n1000/scattered` | 18.9 ms | 19.2 ms | 19.4 ms | 63.3 ms | 194 ms | 3.4x faster |
| `thin_plate_spline/k50/n1000/scattered` | 60 µs | 55.5 µs | 55.1 µs | 132 µs |  | 2.2x faster |
| `cubic/global/n1000/scattered` | 16.4 ms | 16.8 ms | 16.5 ms | 58.7 ms | 187 ms | 3.6x faster |
| `cubic/k50/n1000/scattered` | 62.2 µs | 55.6 µs | 55.9 µs | 130 µs |  | 2.1x faster |
| `quintic/global/n1000/scattered` | 16.9 ms | 17.3 ms | 16.9 ms | 58.1 ms | 184 ms | 3.4x faster |
| `quintic/k50/n1000/scattered` | 58.9 µs | 62.2 µs | 55.2 µs | 131 µs |  | 2.2x faster |
| `multiquadric/global/n1000/scattered` | 17.3 ms | 17.4 ms | 17.5 ms | 59.9 ms | 186 ms | 3.5x faster |
| `multiquadric/k50/n1000/scattered` | 59.5 µs | 55.6 µs | 51.9 µs | 131 µs |  | 2.2x faster |
| `inverse_multiquadric/global/n1000/scattered` | 17.7 ms | 17.9 ms | 18.1 ms | 61.4 ms | 192 ms | 3.5x faster |
| `inverse_multiquadric/k50/n1000/scattered` | 57 µs | 54.6 µs | 55 µs | 129 µs |  | 2.3x faster |
| `inverse_quadratic/global/n1000/scattered` | 16.9 ms | 17.1 ms | 17 ms | 59.9 ms |  | 3.6x faster |
| `inverse_quadratic/k50/n1000/scattered` | 58.4 µs | 49.8 µs | 54.1 µs | 128 µs |  | 2.2x faster |
| `gaussian/global/n1000/scattered` | 18.7 ms | 19.1 ms | 19.3 ms | 62.3 ms | 196 ms | 3.3x faster |
| `gaussian/k50/n1000/scattered` | 57 µs | 56 µs | 54.1 µs | 129 µs |  | 2.3x faster |
| `thin_plate_spline/global/n3000/scattered` | 389 ms | 386 ms | 406 ms | 1.64 s | 6.33 s | 4.2x faster |
| `thin_plate_spline/k50/n10000/scattered` | 746 µs | 732 µs | 637 µs | 1.62 ms |  | 2.2x faster |
| `thin_plate_spline/k50/n10000/image` | 731 µs | 746 µs | 671 µs | 1.68 ms |  | 2.3x faster |
| `thin_plate_spline/k50/n100000/scattered` | 8.2 ms | 8.29 ms | 8.13 ms | 24.4 ms |  | 3x faster |
| `thin_plate_spline/k50/n100000/image` | 8.48 ms | 8.55 ms | 8.97 ms | 24 ms |  | 2.8x faster |
| `linear/global/ghcn-tmax-20250101/heldout` | 16.6 s | 16.4 s | 16.8 s | 112 s | 272 s | 6.8x faster |
| `linear/k50/ghcn-tmax-20250101/heldout` | 750 µs | 764 µs | 698 µs | 1.71 ms |  | 2.3x faster |
| `thin_plate_spline/global/ghcn-tmax-20250101/heldout` | 18.3 s | 17.6 s | 18 s | 114 s |  | 6.2x faster |
| `thin_plate_spline/k50/ghcn-tmax-20250101/heldout` | 780 µs | 762 µs | 703 µs | 1.71 ms |  | 2.2x faster |
| `cubic/global/ghcn-tmax-20250101/heldout` | 18.7 s | 17.5 s | 17.7 s | 124 s | 285 s | 6.6x faster |
| `cubic/k50/ghcn-tmax-20250101/heldout` | 801 µs | 803 µs | 705 µs | 1.69 ms |  | 2.1x faster |
| `quintic/global/ghcn-tmax-20250101/heldout` | 16.3 s | 16.6 s | 19.2 s | 119 s | 275 s | 7.3x faster |
| `quintic/k50/ghcn-tmax-20250101/heldout` | 827 µs | 758 µs | 724 µs | 1.7 ms |  | 2.1x faster |

## Memory

Megabytes. "Data read" is the process with the data of the case in memory, before anything is built. The columns after it are what was added to that: the highest while building, with the interpolator built, the highest while evaluating. "Highest" is the highest of the whole run, repeats included.

| case | library | data read | while building | built | while evaluating | highest |
| :--- | :--- | ---: | ---: | ---: | ---: | ---: |
| `oned/cubic_spline/n10000/sorted` | interpx | 3.4 | +1.9 | +1.9 | +1.9 | 5.9 |
| `oned/cubic_spline/n10000/sorted` | interpx, parallel | 3.3 | +1.9 | +2.0 | +2.0 | 6.4 |
| `oned/cubic_spline/n10000/sorted` | SciPy | 71.0 | +0.5 | +2.0 | +2.0 | 73.1 |
| `oned/cubic_spline/n10000/sorted` | scirs2-interpolate | 4.9 | +0.1 | +0.3 | +0.9 | 5.9 |
| `oned/cubic_spline/n10000/sorted` | ndarray-interp | 3.3 | +1.6 | +1.6 | +1.6 | 5.5 |
| `grid/linear/4d/max200/eval10000/ascending` | interpx | 211.8 | +192.9 | +192.9 | +192.9 | 404.8 |
| `grid/linear/4d/max200/eval10000/ascending` | interpx, parallel | 233.9 | +170.8 | +170.8 | +170.8 | 405.3 |
| `grid/linear/4d/max200/eval10000/ascending` | interpx, `from_owned` | 214.1 | +190.9 | +190.9 | +190.9 | 405.0 |
| `grid/linear/4d/max200/eval10000/ascending` | SciPy | 271.1 | +0.1 | +0.1 | +1.7 | 273.1 |
| `grid/linear/4d/max200/eval10000/ascending` | scirs2-interpolate | 210.2 | +190.8 | +194.7 | +194.7 | 404.9 |
| `grid/linear/4d/max200/eval10000/ascending` | interpn | 217.7 | +0.0 | +0.0 | +0.0 | 399.5 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpx | 5.5 | +2.3 | +2.7 | +2.7 | 8.4 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpx, parallel | 4.3 | +3.6 | +4.0 | +4.0 | 8.9 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpx, `from_owned` | 5.6 | +2.9 | +2.7 | +2.7 | 8.5 |
| `grid/cubic/2d/each400/eval1000/ascending` | SciPy | 70.9 | +85.9 | +6.7 | +6.7 | 160.4 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpn | 6.0 | +0.0 | +0.0 | +0.0 | 6.0 |
| `grid/linear/2d/each800/eval10000/ascending` | interpx | 10.3 | +0.0 | +4.5 | +4.5 | 14.9 |
| `grid/linear/2d/each800/eval10000/ascending` | interpx, parallel | 14.0 | +0.0 | +0.8 | +0.8 | 15.4 |
| `grid/linear/2d/each800/eval10000/ascending` | interpx, `from_owned` | 12.2 | +0.7 | +2.7 | +2.7 | 15.0 |
| `grid/linear/2d/each800/eval10000/ascending` | SciPy | 75.2 | +1.2 | +1.2 | +1.2 | 77.1 |
| `grid/linear/2d/each800/eval10000/ascending` | scirs2-interpolate | 11.0 | +0.0 | +3.8 | +3.9 | 14.9 |
| `grid/linear/2d/each800/eval10000/ascending` | ndarray-interp | 12.5 | +0.0 | +0.0 | +0.0 | 12.5 |
| `grid/linear/2d/each800/eval10000/ascending` | interpn | 13.7 | +0.0 | +0.0 | +0.0 | 13.7 |
| `nearest_nd/n100000/image500` | interpx | 13.0 | +2.1 | +1.4 | +3.3 | 18.5 |
| `nearest_nd/n100000/image500` | interpx, parallel | 13.0 | +2.2 | +1.4 | +3.8 | 18.8 |
| `nearest_nd/n100000/image500` | interpx, `from_owned` | 13.3 | +2.5 | +1.3 | +3.3 | 18.6 |
| `nearest_nd/n100000/image500` | SciPy | 77.1 | +4.7 | +4.7 | +13.0 | 91.8 |
| `rbf/thin_plate_spline/global/n3000/scattered` | interpx | 3.3 | +77.0 | +4.8 | +4.8 | 80.3 |
| `rbf/thin_plate_spline/global/n3000/scattered` | interpx, parallel | 3.1 | +77.1 | +5.0 | +5.3 | 80.7 |
| `rbf/thin_plate_spline/global/n3000/scattered` | interpx, `from_owned` | 3.4 | +77.0 | +4.7 | +4.8 | 80.4 |
| `rbf/thin_plate_spline/global/n3000/scattered` | SciPy | 71.0 | +73.1 | +1.0 | +9.1 | 144.2 |
| `rbf/thin_plate_spline/global/n3000/scattered` | scirs2-interpolate | 2.9 | +289.7 | +1.6 | +1.6 | 292.5 |
| `rbf/thin_plate_spline/k50/n100000/image` | interpx | 6.7 | +3.7 | +4.7 | +5.0 | 11.9 |
| `rbf/thin_plate_spline/k50/n100000/image` | interpx, parallel | 6.2 | +4.1 | +5.2 | +6.4 | 12.7 |
| `rbf/thin_plate_spline/k50/n100000/image` | interpx, `from_owned` | 5.9 | +6.7 | +5.5 | +5.8 | 13.5 |
| `rbf/thin_plate_spline/k50/n100000/image` | SciPy | 72.4 | +6.1 | +6.4 | +22.2 | 98.5 |

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
