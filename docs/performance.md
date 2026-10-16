# Performance: interpx against SciPy and other Rust crates

<!-- Written by compare/perf/report.py from a run of compare/perf/run.py. Do not edit by hand: run it again. -->

Measured on 2026-10-17, at commit `c5fe29c7de`.

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
- Command: `perf/run.py`.

## Where interpx is slower

| family | operation | cases | interpx fastest | slower than SciPy | slower than a Rust crate |
| :--- | :--- | ---: | ---: | ---: | ---: |
| One variable | evaluation | 112 | 112 | 0 | 0 |
| One variable | construction | 112 | 77 | 8 | 31 |
| `PPoly` | evaluation | 4 | 4 | 0 | 0 |
| `PPoly` | construction | 4 | 4 | 0 | 0 |
| `BSpline` | evaluation | 4 | 4 | 0 | 0 |
| `BSpline` | construction | 4 | 4 | 0 | 0 |
| `BSplineCurve` | evaluation | 2 | 2 | 0 | 0 |
| `BSplineCurve` | construction | 2 | 2 | 0 | 0 |
| Regular grid | evaluation | 92 | 75 | 0 | 17 |
| Regular grid | construction | 92 | 28 | 38 | 46 |
| `NearestNDInterpolator` | evaluation | 16 | 16 | 0 | 0 |
| `NearestNDInterpolator` | construction | 16 | 16 | 0 | 0 |
| `RBFInterpolator` | evaluation | 61 | 61 | 0 | 0 |
| `RBFInterpolator` | construction | 61 | 61 | 0 | 0 |

197 comparisons (a case, an operation, and a library), by method. "Same values" says whether that library's values are SciPy's up to rounding: where they are not, it computes another interpolant under the same name, and its time is not that of the same work.

| method | operation | against | same values | cases | slower by | the largest |
| :--- | :--- | :--- | :--- | ---: | ---: | :--- |
| `grid/cubic` | construction | interpn | no | 8 | 1.65e+03x to 9.65e+07x | `grid/cubic/2d/auckland-ba32-0401/heldout/ascending`: 89.5 ms against 0 |
| `grid/linear` | construction | interpn | yes | 35 | 156x to 4.18e+07x | `grid/linear/4d/max200/eval100/ascending`: 53.1 ms against 1.27 ns |
| `grid/nearest` | construction | interpn | yes | 3 | 864x to 2.08e+07x | `grid/nearest/2d/auckland-ba32-0401/heldout/ascending`: 16.2 ms against 0 |
| `grid/linear` | construction | SciPy | yes | 35 | 1.17x to 5.69e+03x | `grid/linear/4d/max200/eval10/flipped`: 114 ms against 20.1 µs |
| `grid/linear` | construction | ndarray-interp | yes | 15 | 1.55x to 2.08e+03x | `grid/linear/2d/auckland-ba32-0401/heldout/ascending`: 17.2 ms against 8.26 µs |
| `grid/nearest` | construction | SciPy | yes | 2 | 2.8x to 2.07e+03x | `grid/nearest/2d/auckland-ba32-0401/heldout/ascending`: 16.2 ms against 7.83 µs |
| `grid/pchip` | construction | SciPy | yes | 1 | 60x | `grid/pchip/3d/each40/eval10000/ascending`: 501 µs against 8.35 µs |
| `oned/linear` | construction | ndarray-interp | yes | 14 | 1.4x to 6.6x | `oned/linear/n10000/sorted`: 93.4 µs against 14.1 µs |
| `oned/bspline1` | construction | SciPy | yes | 2 | 2.75x to 2.85x | `oned/bspline1/n10000/sorted`: 66 µs against 23.1 µs |
| `oned/pchip` | construction | scirs2-interpolate | yes | 8 | 1.49x to 2.82x | `oned/pchip/n1000/shuffled`: 11.7 µs against 4.16 µs |
| `grid/linear` | construction | scirs2-interpolate | yes | 26 | 1.12x to 2.5x | `grid/linear/2d/max200/eval1000/ascending`: 7.42 µs against 2.97 µs |
| `oned/cubic_spline` | construction | scirs2-interpolate | no | 3 | 1.13x to 2.29x | `oned/cubic_spline/n1000/shuffled`: 38.6 µs against 16.9 µs |
| `grid/nearest` | construction | scirs2-interpolate | yes | 2 | 2.04x to 2.27x | `grid/nearest/3d/each40/eval10000/ascending`: 20.9 µs against 9.19 µs |
| `oned/cubic_spline` | construction | ndarray-interp | yes | 1 | 2.22x | `oned/cubic_spline/n1000/sorted`: 38.6 µs against 17.4 µs |
| `oned/linear` | construction | SciPy | yes | 4 | 1.85x to 2.13x | `oned/linear/n10000/sorted-uneven`: 97.9 µs against 46 µs |
| `oned/akima` | construction | scirs2-interpolate | yes | 6 | 1.39x to 2.12x | `oned/akima/n1000/shuffled`: 11.8 µs against 5.58 µs |
| `oned/linear` | construction | scirs2-interpolate | yes | 13 | 1.16x to 1.73x | `oned/linear/n1000/shuffled`: 2.12 µs against 1.22 µs |
| `grid/linear` | evaluation | interpn | yes | 16 | 1.13x to 1.48x | `grid/linear/2d/max10/eval1000/ascending`: 13.7 µs against 9.25 µs |
| `oned/nearest` | construction | SciPy | yes | 2 | 1.14x to 1.15x | `oned/nearest/n10000/shuffled`: 67.9 µs against 58.9 µs |
| `grid/nearest` | evaluation | interpn | yes | 1 | 1.11x | `grid/nearest/2d/each40/eval10000/ascending`: 151 µs against 136 µs |

<details><summary>Every comparison where interpx is slower</summary>

"Values apart" is the distance of that library's values from SciPy's, as a share of the largest value.

| case | operation | against | interpx | the other | slower by | values apart |
| :--- | :--- | :--- | ---: | ---: | ---: | ---: |
| `grid/cubic/2d/auckland-ba32-0401/heldout/ascending` | construction | interpn | 89.5 ms | 0 | 96487855.16x | 5e-03 |
| `grid/linear/4d/max200/eval100/ascending` | construction | interpn | 53.1 ms | 1.27 ns | 41829691.08x | 3e-16 |
| `grid/linear/4d/max200/eval1000/ascending` | construction | interpn | 48.4 ms | 1.27 ns | 38160699.08x | 3e-16 |
| `grid/linear/4d/max200/eval10000/ascending` | construction | interpn | 47.4 ms | 1.27 ns | 37355598.77x | 4e-16 |
| `grid/linear/4d/max200/eval10/ascending` | construction | interpn | 46.9 ms | 1.27 ns | 36963800.62x | 2e-16 |
| `grid/linear/2d/auckland-ba32-0401/heldout/ascending` | construction | interpn | 17.2 ms | 0 | 21964544.00x | 0e+00 |
| `grid/nearest/2d/auckland-ba32-0401/heldout/ascending` | construction | interpn | 16.2 ms | 0 | 20752640.00x | 0e+00 |
| `grid/linear/3d/max200/eval10/ascending` | construction | interpn | 1.65 ms | 0 | 1880860.44x | 3e-16 |
| `grid/linear/3d/max200/eval1000/ascending` | construction | interpn | 1.62 ms | 0 | 1846044.44x | 3e-16 |
| `grid/linear/3d/max200/eval10000/ascending` | construction | interpn | 1.75 ms | 0 | 1793792.00x | 5e-16 |
| `grid/cubic/2d/each400/eval1000/ascending` | construction | interpn | 1.39 ms | 0 | 1782656.00x | 5e-01 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | interpn | 1.08 ms | 0 | 1581787.43x | 3e-16 |
| `grid/linear/3d/max200/eval100/ascending` | construction | interpn | 1.67 ms | 1.07 ns | 1558993.45x | 3e-16 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | interpn | 307 µs | 0 | 448512.00x | 3e-16 |
| `grid/cubic/3d/each40/eval10000/ascending` | construction | interpn | 450 µs | 1.27 ns | 354461.54x | 5e-01 |
| `grid/cubic/2d/each200/eval1000/ascending` | construction | interpn | 170 µs | 0 | 193763.56x | 2e-01 |
| `grid/cubic/2d/each100/eval1000/ascending` | construction | interpn | 74.4 µs | 0 | 84679.11x | 2e-01 |
| `grid/linear/2d/max200/eval10000/ascending` | construction | interpn | 56.3 µs | 0 | 72064.00x | 3e-16 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | interpn | 20.9 µs | 0 | 23779.56x | 0e+00 |
| `grid/linear/3d/max40/eval10000/ascending` | construction | interpn | 16.9 µs | 0 | 19200.00x | 3e-16 |
| `grid/linear/3d/max40/eval1000/ascending` | construction | interpn | 15.4 µs | 0 | 17578.67x | 2e-16 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | interpn | 7.42 µs | 0 | 10852.57x | 2e-16 |
| `grid/linear/4d/max40/eval1000/ascending` | construction | interpn | 13.5 µs | 1.27 ns | 10658.46x | 6e-16 |
| `grid/linear/4d/max40/eval100/ascending` | construction | interpn | 13.1 µs | 1.27 ns | 10348.31x | 2e-16 |
| `grid/linear/4d/max40/eval10/ascending` | construction | interpn | 13.1 µs | 1.27 ns | 10328.62x | 2e-16 |
| `grid/linear/4d/max40/eval10000/ascending` | construction | interpn | 13.1 µs | 1.27 ns | 10304.00x | 5e-16 |
| `grid/cubic/2d/each40/eval1000/ascending` | construction | interpn | 8.95 µs | 0 | 10183.11x | 7e-01 |
| `grid/cubic/2d/each40/eval10000/ascending` | construction | interpn | 8.37 µs | 0 | 9520.00x | 7e-01 |
| `grid/linear/2d/max200/eval100/ascending` | construction | interpn | 6.85 µs | 0 | 8764.00x | 1e-16 |
| `grid/linear/2d/max200/eval10/ascending` | construction | interpn | 6.84 µs | 0 | 8754.00x | 2e-16 |
| `grid/linear/4d/max200/eval10/flipped` | construction | SciPy | 114 ms | 20.1 µs | 5692.45x |  |
| `grid/linear/4d/max200/eval100/ascending` | construction | SciPy | 53.1 ms | 9.4 µs | 5649.37x |  |
| `grid/linear/4d/max200/eval1000/ascending` | construction | SciPy | 48.4 ms | 9.34 µs | 5184.88x |  |
| `grid/linear/4d/max200/eval1000/flipped` | construction | SciPy | 104 ms | 20.3 µs | 5115.83x |  |
| `grid/linear/4d/max200/eval10000/flipped` | construction | SciPy | 102 ms | 20.1 µs | 5067.72x |  |
| `grid/linear/4d/max200/eval10/ascending` | construction | SciPy | 46.9 ms | 9.31 µs | 5039.11x |  |
| `grid/linear/4d/max200/eval10000/ascending` | construction | SciPy | 47.4 ms | 9.46 µs | 5013.45x |  |
| `grid/linear/4d/max200/eval100/flipped` | construction | SciPy | 104 ms | 21.9 µs | 4756.66x |  |
| `grid/linear/3d/max40/eval100/ascending` | construction | interpn | 3.8 µs | 0 | 4324.44x | 2e-16 |
| `grid/linear/3d/max40/eval10/ascending` | construction | interpn | 2.84 µs | 0 | 2904.00x | 2e-16 |
| `grid/linear/2d/auckland-ba32-0401/heldout/ascending` | construction | SciPy | 17.2 ms | 7.75 µs | 2215.06x |  |
| `grid/linear/2d/auckland-ba32-0401/heldout/ascending` | construction | ndarray-interp | 17.2 ms | 8.26 µs | 2076.83x | 0e+00 |
| `grid/nearest/2d/auckland-ba32-0401/heldout/ascending` | construction | SciPy | 16.2 ms | 7.83 µs | 2071.54x |  |
| `grid/cubic/2d/each10/eval1000/ascending` | construction | interpn | 1.45 µs | 0 | 1646.22x | 3e-01 |
| `grid/nearest/2d/each40/eval10000/ascending` | construction | interpn | 591 ns | 0 | 864.00x | 0e+00 |
| `grid/linear/2d/max40/eval10/ascending` | construction | interpn | 360 ns | 0 | 526.86x | 2e-16 |
| `grid/linear/2d/max40/eval1000/ascending` | construction | interpn | 359 ns | 0 | 525.71x | 3e-16 |
| `grid/linear/2d/max40/eval100/ascending` | construction | interpn | 389 ns | 0 | 498.00x | 2e-16 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | ndarray-interp | 1.08 ms | 2.23 µs | 485.64x | 6e-16 |
| `grid/linear/2d/max40/eval10000/ascending` | construction | interpn | 359 ns | 0 | 459.00x | 4e-16 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | ndarray-interp | 307 µs | 1.14 µs | 268.43x | 5e-16 |
| `grid/linear/3d/max200/eval10000/ascending` | construction | SciPy | 1.75 ms | 7.87 µs | 222.62x |  |
| `grid/linear/3d/max10/eval1000/ascending` | construction | interpn | 195 ns | 0 | 222.33x | 5e-16 |
| `grid/linear/3d/max200/eval100/ascending` | construction | SciPy | 1.67 ms | 7.6 µs | 220.45x |  |
| `grid/linear/3d/max200/eval10/ascending` | construction | SciPy | 1.65 ms | 7.6 µs | 217.56x |  |
| `grid/linear/3d/max10/eval100/ascending` | construction | interpn | 197 ns | 0 | 201.30x | 3e-16 |
| `grid/linear/3d/max200/eval1000/flipped` | construction | SciPy | 3.17 ms | 15.9 µs | 199.09x |  |
| `grid/linear/3d/max200/eval1000/ascending` | construction | SciPy | 1.62 ms | 8.3 µs | 195.48x |  |
| `grid/linear/2d/max10/eval10000/ascending` | construction | interpn | 133 ns | 0 | 194.14x | 3e-16 |
| `grid/linear/2d/max10/eval1000/ascending` | construction | interpn | 132 ns | 0 | 193.57x | 3e-16 |
| `grid/linear/3d/max200/eval100/flipped` | construction | SciPy | 3.14 ms | 16.3 µs | 193.18x |  |
| `grid/linear/3d/max200/eval10000/flipped` | construction | SciPy | 2.99 ms | 15.9 µs | 188.72x |  |
| `grid/linear/3d/max10/eval10/ascending` | construction | interpn | 179 ns | 0 | 183.60x | 2e-16 |
| `grid/linear/3d/max10/eval10000/ascending` | construction | interpn | 189 ns | 1.07 ns | 175.64x | 4e-16 |
| `grid/linear/3d/max200/eval10/flipped` | construction | SciPy | 3.08 ms | 17.7 µs | 173.38x |  |
| `grid/linear/2d/each800/eval10000/ascending` | construction | SciPy | 1.08 ms | 6.37 µs | 169.87x |  |
| `grid/linear/2d/max10/eval10/ascending` | construction | interpn | 122 ns | 0 | 156.25x | 1e-16 |
| `grid/linear/2d/max10/eval100/ascending` | construction | interpn | 122 ns | 0 | 155.88x | 2e-16 |
| `grid/linear/2d/max200/eval10000/ascending` | construction | ndarray-interp | 56.3 µs | 481 ns | 117.06x | 5e-16 |
| `grid/pchip/3d/each40/eval10000/ascending` | construction | SciPy | 501 µs | 8.35 µs | 59.96x |  |
| `grid/linear/2d/each400/eval10000/ascending` | construction | SciPy | 307 µs | 6 µs | 51.11x |  |
| `grid/linear/2d/max200/eval1000/ascending` | construction | ndarray-interp | 7.42 µs | 462 ns | 16.04x | 4e-16 |
| `grid/linear/2d/max200/eval100/ascending` | construction | ndarray-interp | 6.85 µs | 518 ns | 13.23x | 3e-16 |
| `grid/linear/2d/max200/eval10/ascending` | construction | ndarray-interp | 6.84 µs | 518 ns | 13.21x | 2e-16 |
| `grid/linear/2d/max200/eval10000/ascending` | construction | SciPy | 56.3 µs | 5.81 µs | 9.69x |  |
| `oned/linear/n10000/sorted` | construction | ndarray-interp | 93.4 µs | 14.1 µs | 6.60x | 0e+00 |
| `oned/linear/n10000/shuffled-uneven` | construction | ndarray-interp | 98.7 µs | 19.4 µs | 5.07x | 0e+00 |
| `oned/linear/n10000/sorted-uneven` | construction | ndarray-interp | 97.9 µs | 20.3 µs | 4.83x | 0e+00 |
| `grid/linear/2d/max200/eval10000/flipped` | construction | SciPy | 55.7 µs | 12.2 µs | 4.55x |  |
| `oned/bspline1/n10000/sorted` | construction | SciPy | 66 µs | 23.1 µs | 2.85x |  |
| `grid/linear/4d/max40/eval1000/flipped` | construction | SciPy | 56.2 µs | 19.8 µs | 2.83x |  |
| `oned/pchip/n1000/shuffled` | construction | scirs2-interpolate | 11.7 µs | 4.16 µs | 2.82x | 2e-16 |
| `grid/linear/4d/max40/eval10000/flipped` | construction | SciPy | 56.4 µs | 20 µs | 2.82x |  |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | SciPy | 20.9 µs | 7.46 µs | 2.80x |  |
| `oned/bspline1/n10000/shuffled` | construction | SciPy | 69.5 µs | 25.2 µs | 2.75x |  |
| `oned/linear/n10000/shuffled` | construction | ndarray-interp | 89.9 µs | 32.8 µs | 2.74x | 0e+00 |
| `oned/pchip/n1000/sorted` | construction | scirs2-interpolate | 12.4 µs | 4.56 µs | 2.72x | 2e-16 |
| `grid/linear/4d/max40/eval10/flipped` | construction | SciPy | 51.9 µs | 19.8 µs | 2.62x |  |
| `grid/linear/4d/max40/eval100/flipped` | construction | SciPy | 52.1 µs | 20.4 µs | 2.56x |  |
| `grid/linear/2d/max40/eval100/ascending` | construction | ndarray-interp | 389 ns | 154 ns | 2.53x | 2e-16 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | scirs2-interpolate | 7.42 µs | 2.97 µs | 2.50x | 2e-16 |
| `grid/linear/4d/max40/eval1000/ascending` | construction | scirs2-interpolate | 13.5 µs | 5.51 µs | 2.45x | 2e-16 |
| `grid/linear/4d/max40/eval10/ascending` | construction | scirs2-interpolate | 13.1 µs | 5.51 µs | 2.38x | 2e-16 |
| `grid/linear/4d/max40/eval10000/ascending` | construction | scirs2-interpolate | 13.1 µs | 5.51 µs | 2.37x | 4e-16 |
| `grid/linear/4d/max40/eval100/ascending` | construction | scirs2-interpolate | 13.1 µs | 5.54 µs | 2.37x | 3e-16 |
| `grid/linear/2d/max40/eval1000/ascending` | construction | ndarray-interp | 359 ns | 152 ns | 2.37x | 4e-16 |
| `oned/cubic_spline/n1000/shuffled` | construction | scirs2-interpolate | 38.6 µs | 16.9 µs | 2.29x | 4e-04 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | scirs2-interpolate | 20.9 µs | 9.19 µs | 2.27x | 0e+00 |
| `oned/linear/n10/shuffled` | construction | ndarray-interp | 93.8 ns | 41.5 ns | 2.26x | 0e+00 |
| `oned/linear/n10/sorted` | construction | ndarray-interp | 96 ns | 42.6 ns | 2.25x | 0e+00 |
| `grid/linear/2d/max200/eval10/ascending` | construction | scirs2-interpolate | 6.84 µs | 3.06 µs | 2.23x | 1e-16 |
| `oned/cubic_spline/n1000/sorted` | construction | ndarray-interp | 38.6 µs | 17.4 µs | 2.22x | 2e-16 |
| `grid/linear/2d/max40/eval10/ascending` | construction | ndarray-interp | 360 ns | 167 ns | 2.16x | 2e-16 |
| `grid/linear/2d/max40/eval10000/ascending` | construction | ndarray-interp | 359 ns | 167 ns | 2.15x | 4e-16 |
| `oned/cubic_spline/n1000/sorted` | construction | scirs2-interpolate | 38.6 µs | 18.1 µs | 2.13x | 4e-04 |
| `oned/linear/n10000/sorted-uneven` | construction | SciPy | 97.9 µs | 46 µs | 2.13x |  |
| `grid/linear/2d/max200/eval100/ascending` | construction | scirs2-interpolate | 6.85 µs | 3.23 µs | 2.12x | 2e-16 |
| `oned/akima/n1000/shuffled` | construction | scirs2-interpolate | 11.8 µs | 5.58 µs | 2.12x | 3e-17 |
| `oned/linear/n10000/shuffled-uneven` | construction | SciPy | 98.7 µs | 47.1 µs | 2.10x |  |
| `grid/nearest/2d/each40/eval10000/ascending` | construction | scirs2-interpolate | 591 ns | 290 ns | 2.04x | 0e+00 |
| `oned/linear/n10000/shuffled` | construction | SciPy | 89.9 µs | 45.8 µs | 1.96x |  |
| `grid/linear/3d/max40/eval1000/flipped` | construction | SciPy | 30.8 µs | 15.7 µs | 1.96x |  |
| `oned/akima/n100/shuffled` | construction | scirs2-interpolate | 978 ns | 524 ns | 1.87x | 1e-16 |
| `oned/linear/n10000/sorted` | construction | SciPy | 93.4 µs | 50.5 µs | 1.85x |  |
| `grid/linear/3d/max40/eval10000/ascending` | construction | scirs2-interpolate | 16.9 µs | 9.46 µs | 1.78x | 2e-16 |
| `oned/linear/n1000/shuffled` | construction | scirs2-interpolate | 2.12 µs | 1.22 µs | 1.73x | 0e+00 |
| `grid/linear/2d/max10/eval1000/ascending` | construction | ndarray-interp | 132 ns | 78.1 ns | 1.69x | 3e-16 |
| `grid/linear/2d/max40/eval10/ascending` | construction | scirs2-interpolate | 360 ns | 215 ns | 1.67x | 4e-17 |
| `grid/linear/2d/max40/eval100/ascending` | construction | scirs2-interpolate | 389 ns | 234 ns | 1.66x | 3e-16 |
| `grid/linear/2d/max40/eval10000/ascending` | construction | scirs2-interpolate | 359 ns | 216 ns | 1.66x | 2e-16 |
| `grid/linear/2d/max40/eval1000/ascending` | construction | scirs2-interpolate | 359 ns | 217 ns | 1.66x | 3e-16 |
| `oned/pchip/n100/shuffled` | construction | scirs2-interpolate | 967 ns | 587 ns | 1.65x | 2e-16 |
| `oned/linear/n100/shuffled-uneven` | construction | ndarray-interp | 272 ns | 166 ns | 1.64x | 0e+00 |
| `oned/pchip/n100/sorted-uneven` | construction | scirs2-interpolate | 963 ns | 589 ns | 1.63x | 2e-16 |
| `oned/pchip/n100/sorted` | construction | scirs2-interpolate | 967 ns | 593 ns | 1.63x | 2e-16 |
| `oned/pchip/n50/sorted` | construction | scirs2-interpolate | 545 ns | 334 ns | 1.63x | 1e-16 |
| `oned/linear/n100/shuffled` | construction | ndarray-interp | 270 ns | 166 ns | 1.63x | 0e+00 |
| `oned/linear/n100/sorted` | construction | ndarray-interp | 271 ns | 167 ns | 1.63x | 0e+00 |
| `oned/linear/n10000/shuffled-uneven` | construction | scirs2-interpolate | 98.7 µs | 61.1 µs | 1.62x | 6e-17 |
| `oned/linear/n10000/shuffled` | construction | scirs2-interpolate | 89.9 µs | 56 µs | 1.60x | 0e+00 |
| `oned/linear/n10000/sorted-uneven` | construction | scirs2-interpolate | 97.9 µs | 61.1 µs | 1.60x | 6e-17 |
| `oned/linear/n1000/sorted` | construction | scirs2-interpolate | 2.12 µs | 1.33 µs | 1.60x | 0e+00 |
| `oned/pchip/n50/shuffled` | construction | scirs2-interpolate | 545 ns | 344 ns | 1.58x | 1e-16 |
| `oned/linear/n50/sorted` | construction | ndarray-interp | 156 ns | 99.2 ns | 1.57x | 0e+00 |
| `grid/linear/2d/max10/eval10/ascending` | construction | ndarray-interp | 122 ns | 78 ns | 1.56x | 1e-16 |
| `grid/linear/2d/max10/eval100/ascending` | construction | ndarray-interp | 122 ns | 77.9 ns | 1.56x | 4e-16 |
| `oned/akima/n100/sorted` | construction | scirs2-interpolate | 898 ns | 575 ns | 1.56x | 1e-16 |
| `grid/linear/2d/max10/eval10000/ascending` | construction | ndarray-interp | 133 ns | 85.7 ns | 1.55x | 4e-16 |
| `oned/linear/n50/shuffled` | construction | ndarray-interp | 154 ns | 100 ns | 1.54x | 0e+00 |
| `oned/linear/n1000/shuffled` | construction | ndarray-interp | 2.12 µs | 1.38 µs | 1.54x | 0e+00 |
| `oned/akima/n50/sorted` | construction | scirs2-interpolate | 502 ns | 327 ns | 1.54x | 1e-16 |
| `oned/linear/n10000/sorted` | construction | scirs2-interpolate | 93.4 µs | 61 µs | 1.53x | 0e+00 |
| `oned/akima/n50/shuffled` | construction | scirs2-interpolate | 504 ns | 334 ns | 1.51x | 1e-16 |
| `oned/linear/n100/sorted-uneven` | construction | ndarray-interp | 249 ns | 166 ns | 1.50x | 0e+00 |
| `oned/pchip/n100/shuffled-uneven` | construction | scirs2-interpolate | 884 ns | 592 ns | 1.49x | 2e-16 |
| `grid/linear/2d/max10/eval1000/ascending` | evaluation | interpn | 13.7 µs | 9.25 µs | 1.48x | 3e-16 |
| `grid/linear/4d/max40/eval1000/ascending` | construction | SciPy | 13.5 µs | 9.28 µs | 1.46x |  |
| `grid/linear/3d/max10/eval1000/ascending` | evaluation | interpn | 20.8 µs | 14.4 µs | 1.45x | 5e-16 |
| `grid/linear/4d/max40/eval10/ascending` | construction | SciPy | 13.1 µs | 9.22 µs | 1.42x |  |
| `grid/linear/2d/max10/eval10/ascending` | evaluation | interpn | 134 ns | 94.8 ns | 1.42x | 1e-16 |
| `grid/linear/4d/max40/eval100/ascending` | construction | SciPy | 13.1 µs | 9.35 µs | 1.41x |  |
| `oned/linear/n1000/sorted` | construction | ndarray-interp | 2.12 µs | 1.51 µs | 1.40x | 0e+00 |
| `grid/linear/4d/max40/eval10000/ascending` | construction | SciPy | 13.1 µs | 9.35 µs | 1.40x |  |
| `oned/akima/n1000/sorted` | construction | scirs2-interpolate | 12.2 µs | 8.8 µs | 1.39x | 3e-17 |
| `grid/linear/2d/max10/eval10000/ascending` | evaluation | interpn | 128 µs | 92.3 µs | 1.39x | 3e-16 |
| `oned/linear/n100/shuffled-uneven` | construction | scirs2-interpolate | 272 ns | 197 ns | 1.38x | 3e-17 |
| `oned/linear/n100/shuffled` | construction | scirs2-interpolate | 270 ns | 196 ns | 1.38x | 0e+00 |
| `grid/linear/3d/max10/eval100/ascending` | evaluation | interpn | 2.09 µs | 1.53 µs | 1.36x | 3e-16 |
| `grid/linear/3d/max200/eval10000/ascending` | construction | scirs2-interpolate | 1.75 ms | 1.29 ms | 1.36x | 3e-16 |
| `grid/linear/3d/max200/eval1000/ascending` | construction | scirs2-interpolate | 1.62 ms | 1.2 ms | 1.35x | 3e-16 |
| `grid/linear/2d/max10/eval100/ascending` | evaluation | interpn | 1.25 µs | 931 ns | 1.34x | 2e-16 |
| `grid/linear/3d/max10/eval10000/ascending` | evaluation | interpn | 208 µs | 155 µs | 1.34x | 4e-16 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | scirs2-interpolate | 1.08 ms | 809 µs | 1.34x | 2e-16 |
| `grid/linear/4d/max200/eval100/ascending` | construction | scirs2-interpolate | 53.1 ms | 40.4 ms | 1.31x | 1e-16 |
| `grid/linear/3d/max200/eval10/ascending` | construction | scirs2-interpolate | 1.65 ms | 1.26 ms | 1.31x | 3e-16 |
| `oned/linear/n100/sorted` | construction | scirs2-interpolate | 271 ns | 207 ns | 1.31x | 0e+00 |
| `grid/linear/3d/max200/eval100/ascending` | construction | scirs2-interpolate | 1.67 ms | 1.3 ms | 1.29x | 3e-16 |
| `oned/linear/n50/sorted` | construction | scirs2-interpolate | 156 ns | 121 ns | 1.29x | 0e+00 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | SciPy | 7.42 µs | 5.76 µs | 1.29x |  |
| `grid/linear/3d/max10/eval10/ascending` | evaluation | interpn | 199 ns | 155 ns | 1.29x | 2e-16 |
| `grid/linear/2d/max200/eval10/flipped` | construction | SciPy | 14.7 µs | 11.9 µs | 1.23x |  |
| `grid/linear/4d/max200/eval1000/ascending` | construction | scirs2-interpolate | 48.4 ms | 39.6 ms | 1.22x | 3e-16 |
| `grid/linear/4d/max200/eval10000/ascending` | construction | scirs2-interpolate | 47.4 ms | 38.8 ms | 1.22x | 4e-16 |
| `grid/linear/2d/auckland-ba32-0401/heldout/ascending` | construction | scirs2-interpolate | 17.2 ms | 14.1 ms | 1.22x | 0e+00 |
| `grid/linear/2d/max200/eval100/flipped` | construction | SciPy | 15.9 µs | 13.2 µs | 1.21x |  |
| `oned/linear/n10/sorted` | construction | scirs2-interpolate | 96 ns | 80.1 ns | 1.20x | 0e+00 |
| `grid/linear/2d/max200/eval1000/flipped` | construction | SciPy | 14.3 µs | 11.9 µs | 1.20x |  |
| `grid/linear/2d/max40/eval100/ascending` | evaluation | interpn | 1.57 µs | 1.33 µs | 1.19x | 2e-16 |
| `grid/linear/2d/max40/eval10/ascending` | evaluation | interpn | 158 ns | 134 ns | 1.18x | 2e-16 |
| `oned/linear/n50/shuffled` | construction | scirs2-interpolate | 154 ns | 131 ns | 1.18x | 0e+00 |
| `grid/linear/3d/max40/eval10000/ascending` | evaluation | interpn | 234 µs | 200 µs | 1.17x | 3e-16 |
| `grid/linear/2d/max200/eval10/ascending` | construction | SciPy | 6.84 µs | 5.85 µs | 1.17x |  |
| `oned/linear/n100/sorted-uneven` | construction | scirs2-interpolate | 249 ns | 215 ns | 1.16x | 3e-17 |
| `oned/nearest/n10000/shuffled` | construction | SciPy | 67.9 µs | 58.9 µs | 1.15x |  |
| `grid/linear/3d/max40/eval1000/ascending` | construction | scirs2-interpolate | 15.4 µs | 13.4 µs | 1.15x | 2e-16 |
| `grid/linear/2d/max40/eval1000/ascending` | evaluation | interpn | 15.1 µs | 13.2 µs | 1.15x | 3e-16 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | scirs2-interpolate | 307 µs | 268 µs | 1.14x | 2e-16 |
| `oned/nearest/n10000/sorted` | construction | SciPy | 61.1 µs | 53.5 µs | 1.14x |  |
| `grid/linear/4d/max40/eval10/ascending` | evaluation | interpn | 309 ns | 272 ns | 1.14x | 2e-16 |
| `oned/cubic_spline/n100/shuffled-uneven` | construction | scirs2-interpolate | 1.73 µs | 1.53 µs | 1.13x | 6e-03 |
| `grid/linear/3d/max40/eval10/ascending` | evaluation | interpn | 226 ns | 199 ns | 1.13x | 2e-16 |
| `grid/linear/4d/max200/eval10/ascending` | evaluation | interpn | 414 ns | 365 ns | 1.13x | 2e-16 |
| `grid/linear/4d/max200/eval1000/ascending` | evaluation | interpn | 41.1 µs | 36.3 µs | 1.13x | 3e-16 |
| `grid/linear/3d/max10/eval100/ascending` | construction | scirs2-interpolate | 197 ns | 175 ns | 1.13x | 2e-16 |
| `grid/linear/3d/max10/eval1000/ascending` | construction | scirs2-interpolate | 195 ns | 174 ns | 1.12x | 3e-16 |
| `grid/linear/4d/max200/eval10/ascending` | construction | scirs2-interpolate | 46.9 ms | 41.8 ms | 1.12x | 2e-16 |
| `grid/nearest/2d/each40/eval10000/ascending` | evaluation | interpn | 151 µs | 136 µs | 1.11x | 0e+00 |

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

A run (one library on one case, a process of its own) is stopped at 15 minutes. The longest run that finished took 310 s, 34% of that: scirs2-interpolate on `rbf/thin_plate_spline/global/ghcn-tmax-20250101/heldout`.

No run took more than half the limit (450 s).

## Real-world data

Two data sets, each taken whole in the unit its source gives it out in; what was taken and how it was divided was fixed before anything was measured (`compare/perf/real.py`, docs/design.ja.md section 8.4.7).

- `auckland-ba32-0401`: Auckland Part 1 LiDAR 1m DEM (2024), tile BA32_10000_0401. 4800 by 7200 cells of 1 m. The grid is every second row and column (and the last ones); every other cell is a query, and its elevation is the value held out. Licensed by Regional Software Holdings Limited for re-use under CC-BY-4.0. Hosted by Toitū Te Whenua Land Information New Zealand (https://registry.opendata.aws/nz-elevation).
- `ghcn-tmax-20250101`: GHCN-Daily, maximum temperature of 2025-01-01, every station; longitude and latitude in degrees, as they are. Of every 10 stations in the order of the file the last is a query, and its temperature is the value held out. The RBF kernels are those without a shape parameter. NOAA Global Historical Climatology Network Daily (GHCN-D), from https://registry.opendata.aws/noaa-ghcn (CC0 1.0).

These cases are also in the tables of their families below and in the counts above. Here they are with the error: the largest and the root mean square distance from the values held out (metres, and degrees Celsius), and the number of queries a library gave no value. A library missing from a case does not have it, or did not finish (see the runs that failed).

`grid/nearest/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 16.2 ms | 373 ms | 15.5 | 0.178 | 0 |
| interpx, parallel | 16.8 ms | 62.2 ms | 15.5 | 0.178 | 0 |
| interpx, `from_owned` | 5.27 ms | 385 ms | 15.5 | 0.178 | 0 |
| SciPy | 7.83 µs | 1.39 s | 15.5 | 0.178 | 0 |
| scirs2-interpolate | 16.2 ms | 161 s | 15.5 | 0.178 | 0 |
| interpn | 0 | 760 ms | 15.5 | 0.178 | 0 |

`grid/linear/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 17.2 ms | 389 ms | 6.74 | 0.0601 | 0 |
| interpx, parallel | 16.3 ms | 65.4 ms | 6.74 | 0.0601 | 0 |
| interpx, `from_owned` | 5.14 ms | 421 ms | 6.74 | 0.0601 | 0 |
| SciPy | 7.75 µs | 1.16 s | 6.74 | 0.0601 | 0 |
| scirs2-interpolate | 14.1 ms | 70.4 s | 6.74 | 0.0601 | 0 |
| ndarray-interp | 8.26 µs | 374 ms | 6.74 | 0.0601 | 0 |
| interpn | 0 | 1.57 s | 6.74 | 0.0601 | 0 |

`grid/cubic/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 89.5 ms | 843 ms | 6.77 | 0.0507 | 0 |
| interpx, parallel | 80 ms | 129 ms | 6.77 | 0.0507 | 0 |
| interpx, `from_owned` | 66.8 ms | 779 ms | 6.77 | 0.0507 | 0 |
| SciPy | 14.2 s | 32.9 s | 6.77 | 0.0507 | 0 |
| interpn | 0 | 1.9 s | 6.76 | 0.0517 | 0 |

`grid/quintic/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 98.3 ms | 1.31 s | 6.77 | 0.0517 | 0 |
| interpx, parallel | 92.2 ms | 212 ms | 6.77 | 0.0517 | 0 |
| interpx, `from_owned` | 84.2 ms | 1.21 s | 6.77 | 0.0517 | 0 |
| SciPy | 22.5 s | 35.8 s | 6.77 | 0.0517 | 0 |

`grid/pchip/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 77 ms | 1 s | 6.75 | 0.049 | 0 |
| interpx, parallel | 74.8 ms | 185 ms | 6.75 | 0.049 | 0 |
| interpx, `from_owned` | 61.5 ms | 1.04 s | 6.75 | 0.049 | 0 |

`nearest_nd/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 697 µs | 156 µs | 25 | 3.42 | 0 |
| interpx, parallel | 699 µs | 66.1 µs | 25 | 3.42 | 0 |
| interpx, `from_owned` | 659 µs | 157 µs | 25 | 3.42 | 0 |
| SciPy | 1.7 ms | 519 µs | 25 | 3.42 | 0 |

`rbf/linear/global/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 17.8 s | 18.3 ms | 20.8 | 2.74 | 0 |
| interpx, parallel | 16.7 s | 4.28 ms | 20.8 | 2.74 | 0 |
| interpx, `from_owned` | 17.7 s | 24.4 ms | 20.8 | 2.74 | 0 |
| SciPy | 116 s | 50.1 ms | 20.8 | 2.74 | 0 |
| scirs2-interpolate | 285 s | 234 ms | 20.8 | 2.74 | 0 |

`rbf/linear/k50/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 774 µs | 19.8 ms | 20.8 | 2.75 | 0 |
| interpx, parallel | 789 µs | 3.32 ms | 20.8 | 2.75 | 0 |
| interpx, `from_owned` | 741 µs | 18.3 ms | 20.8 | 2.75 | 0 |
| SciPy | 1.79 ms | 46.7 ms | 20.8 | 2.75 | 0 |

`rbf/thin_plate_spline/global/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 18.1 s | 88.5 ms | 102 | 5.67 | 0 |
| interpx, parallel | 17.7 s | 16.1 ms | 102 | 5.67 | 0 |
| interpx, `from_owned` | 18.6 s | 89.8 ms | 102 | 5.67 | 0 |
| SciPy | 115 s | 116 ms | 102 | 5.67 | 0 |
| scirs2-interpolate | 308 s | 400 ms | 93.7 | 5.41 | 0 |

`rbf/thin_plate_spline/k50/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 788 µs | 29.3 ms | 102 | 5.66 | 0 |
| interpx, parallel | 822 µs | 5.18 ms | 102 | 5.66 | 0 |
| interpx, `from_owned` | 680 µs | 29.9 ms | 102 | 5.66 | 0 |
| SciPy | 1.76 ms | 58.6 ms | 102 | 5.66 | 0 |

`rbf/cubic/global/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 17.4 s | 19.3 ms | 587 | 27.1 | 0 |
| interpx, parallel | 16.5 s | 3.3 ms | 587 | 27.1 | 0 |
| interpx, `from_owned` | 16.9 s | 19 ms | 587 | 27.1 | 0 |
| SciPy | 117 s | 51 ms | 597 | 24.4 | 0 |
| scirs2-interpolate | 277 s | 255 ms | 100 | 8.48 | 0 |

`rbf/cubic/k50/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 764 µs | 20.7 ms | 567 | 25.3 | 0 |
| interpx, parallel | 810 µs | 3.72 ms | 567 | 25.3 | 0 |
| interpx, `from_owned` | 723 µs | 20.1 ms | 567 | 25.3 | 0 |
| SciPy | 1.71 ms | 48.4 ms | 567 | 25.3 | 0 |

`rbf/quintic/global/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 16.1 s | 21.2 ms | 2.55e+05 | 1.36e+04 | 0 |
| interpx, parallel | 16.9 s | 4.13 ms | 2.55e+05 | 1.36e+04 | 0 |
| interpx, `from_owned` | 17.1 s | 19.7 ms | 2.55e+05 | 1.36e+04 | 0 |
| SciPy | 106 s | 48.4 ms | 9.39e+04 | 4.92e+03 | 0 |
| scirs2-interpolate | 281 s | 284 ms | 6.56e+16 | 2.85e+15 | 0 |

`rbf/quintic/k50/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 819 µs | 22.1 ms | 3.74e+03 | 199 | 0 |
| interpx, parallel | 826 µs | 3.86 ms | 3.74e+03 | 199 | 0 |
| interpx, `from_owned` | 754 µs | 21.9 ms | 3.74e+03 | 199 | 0 |
| SciPy | 1.75 ms | 52.9 ms | 3.74e+03 | 200 | 0 |

## One variable

Data of `n` points, evaluated at `4 n` points: in order (`sorted`) and in no order (`shuffled`).

### Evaluation

| case | interpx | interpx, parallel | SciPy | scirs2-interpolate | ndarray-interp | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/n10/sorted` | 112 ns | 171 ns | 5.07 µs | 312 ns | 254 ns | 2.3x faster |
| `linear/n10/shuffled` | 169 ns | 242 ns | 5.5 µs | 309 ns | 249 ns | 1.5x faster |
| `nearest/n10/sorted` | 161 ns | 209 ns | 9.42 µs | 282 ns |  | 1.7x faster |
| `nearest/n10/shuffled` | 290 ns | 359 ns | 9.46 µs | 277 ns |  | **1x slower** |
| `previous/n10/sorted` | 161 ns | 200 ns | 10.6 µs |  |  | 66x faster |
| `previous/n10/shuffled` | 261 ns | 320 ns | 10.9 µs |  |  | 42x faster |
| `bspline1/n10/sorted` | 210 ns | 308 ns | 2.72 µs |  |  | 13x faster |
| `bspline1/n10/shuffled` | 578 ns | 597 ns | 2.77 µs |  |  | 4.8x faster |
| `bspline2/n10/sorted` | 267 ns | 339 ns | 3.05 µs |  |  | 11x faster |
| `bspline2/n10/shuffled` | 338 ns | 413 ns | 3.14 µs |  |  | 9.3x faster |
| `bspline3/n10/sorted` | 343 ns | 436 ns | 3.57 µs |  |  | 10x faster |
| `bspline3/n10/shuffled` | 441 ns | 564 ns | 3.58 µs |  |  | 8.1x faster |
| `cubic_spline/n10/sorted` | 139 ns | 199 ns | 3.83 µs | 284 ns | 333 ns | 2x faster |
| `cubic_spline/n10/shuffled` | 204 ns | 291 ns | 3.87 µs | 276 ns | 326 ns | 1.4x faster |
| `pchip/n10/sorted` | 154 ns | 193 ns | 4.12 µs | 319 ns |  | 2.1x faster |
| `pchip/n10/shuffled` | 196 ns | 261 ns | 3.9 µs | 296 ns |  | 1.5x faster |
| `akima/n10/sorted` | 138 ns | 194 ns | 3.89 µs | 272 ns |  | 2x faster |
| `akima/n10/shuffled` | 205 ns | 283 ns | 3.84 µs | 266 ns |  | 1.3x faster |
| `makima/n10/sorted` | 141 ns | 210 ns | 3.83 µs |  |  | 27x faster |
| `makima/n10/shuffled` | 204 ns | 269 ns | 3.79 µs |  |  | 19x faster |
| `linear/n50/sorted` | 608 ns | 675 ns | 6.17 µs | 1.9 µs | 1.23 µs | 2x faster |
| `linear/n50/shuffled` | 808 ns | 857 ns | 8.77 µs | 1.81 µs | 1.38 µs | 1.7x faster |
| `nearest/n50/sorted` | 768 ns | 809 ns | 12.3 µs | 1.68 µs |  | 2.2x faster |
| `nearest/n50/shuffled` | 1.53 µs | 1.75 µs | 11.3 µs | 1.83 µs |  | 1.2x faster |
| `previous/n50/sorted` | 832 ns | 849 ns | 13.9 µs |  |  | 17x faster |
| `previous/n50/shuffled` | 1.35 µs | 1.55 µs | 12.8 µs |  |  | 9.5x faster |
| `bspline1/n50/sorted` | 1.15 µs | 1.39 µs | 3.93 µs |  |  | 3.4x faster |
| `bspline1/n50/shuffled` | 3.05 µs | 2.85 µs | 4.9 µs |  |  | 1.6x faster |
| `bspline2/n50/sorted` | 1.35 µs | 1.6 µs | 6.22 µs |  |  | 4.6x faster |
| `bspline2/n50/shuffled` | 2.02 µs | 2.1 µs | 6.74 µs |  |  | 3.3x faster |
| `bspline3/n50/sorted` | 1.82 µs | 2.07 µs | 8.96 µs |  |  | 4.9x faster |
| `bspline3/n50/shuffled` | 2.81 µs | 3.15 µs | 9.38 µs |  |  | 3.3x faster |
| `cubic_spline/n50/sorted` | 776 ns | 769 ns | 5.77 µs | 4.35 µs | 1.61 µs | 2.1x faster |
| `cubic_spline/n50/shuffled` | 997 ns | 1.04 µs | 5.63 µs | 4.13 µs | 1.6 µs | 1.6x faster |
| `pchip/n50/sorted` | 720 ns | 753 ns | 5.44 µs | 4.31 µs |  | 6x faster |
| `pchip/n50/shuffled` | 962 ns | 1.03 µs | 5.68 µs | 3.71 µs |  | 3.9x faster |
| `akima/n50/sorted` | 726 ns | 753 ns | 5.27 µs | 4.31 µs |  | 5.9x faster |
| `akima/n50/shuffled` | 980 ns | 1.04 µs | 5.62 µs | 4.27 µs |  | 4.4x faster |
| `makima/n50/sorted` | 780 ns | 817 ns | 5.35 µs |  |  | 6.9x faster |
| `makima/n50/shuffled` | 1.07 µs | 1.13 µs | 5.69 µs |  |  | 5.3x faster |
| `linear/n100/sorted` | 1.28 µs | 1.48 µs | 7.28 µs | 4.02 µs | 2.52 µs | 2x faster |
| `linear/n100/shuffled` | 1.65 µs | 1.7 µs | 12.7 µs | 3.89 µs | 2.72 µs | 1.7x faster |
| `nearest/n100/sorted` | 1.65 µs | 1.59 µs | 13.6 µs | 3.74 µs |  | 2.3x faster |
| `nearest/n100/shuffled` | 3.06 µs | 3.5 µs | 13.9 µs | 4.06 µs |  | 1.3x faster |
| `previous/n100/sorted` | 1.5 µs | 1.72 µs | 16.4 µs |  |  | 11x faster |
| `previous/n100/shuffled` | 2.95 µs | 3.1 µs | 15.7 µs |  |  | 5.3x faster |
| `bspline1/n100/sorted` | 2.31 µs | 2.6 µs | 5.85 µs |  |  | 2.5x faster |
| `bspline1/n100/shuffled` | 5.58 µs | 5.76 µs | 11.9 µs |  |  | 2.1x faster |
| `bspline2/n100/sorted` | 2.65 µs | 3.29 µs | 10 µs |  |  | 3.8x faster |
| `bspline2/n100/shuffled` | 4.44 µs | 4.62 µs | 12.7 µs |  |  | 2.9x faster |
| `bspline3/n100/sorted` | 3.61 µs | 4.12 µs | 15 µs |  |  | 4.1x faster |
| `bspline3/n100/shuffled` | 6.15 µs | 6.29 µs | 18.5 µs |  |  | 3x faster |
| `cubic_spline/n100/sorted` | 1.46 µs | 1.66 µs | 6.98 µs | 18.9 µs | 3.5 µs | 2.4x faster |
| `cubic_spline/n100/shuffled` | 1.95 µs | 2.04 µs | 8.85 µs | 19.4 µs | 3.2 µs | 1.6x faster |
| `pchip/n100/sorted` | 1.59 µs | 1.66 µs | 7.55 µs | 17.5 µs |  | 4.8x faster |
| `pchip/n100/shuffled` | 2.14 µs | 2.26 µs | 8.34 µs | 17.9 µs |  | 3.9x faster |
| `akima/n100/sorted` | 1.47 µs | 1.66 µs | 7.13 µs | 18.3 µs |  | 4.8x faster |
| `akima/n100/shuffled` | 1.95 µs | 2.12 µs | 8.02 µs | 17.6 µs |  | 4.1x faster |
| `makima/n100/sorted` | 1.49 µs | 1.64 µs | 7.76 µs |  |  | 5.2x faster |
| `makima/n100/shuffled` | 1.95 µs | 2.05 µs | 8.12 µs |  |  | 4.2x faster |
| `linear/n1000/sorted` | 13.5 µs | 14.3 µs | 31.2 µs | 122 µs | 27 µs | 2x faster |
| `linear/n1000/shuffled` | 17.6 µs | 16.4 µs | 143 µs | 217 µs | 24.7 µs | 1.4x faster |
| `nearest/n1000/sorted` | 16.3 µs | 16.5 µs | 56.1 µs | 158 µs |  | 3.4x faster |
| `nearest/n1000/shuffled` | 33.7 µs | 34 µs | 206 µs | 259 µs |  | 6.1x faster |
| `previous/n1000/sorted` | 15.2 µs | 16.6 µs | 62.9 µs |  |  | 4.1x faster |
| `previous/n1000/shuffled` | 27.4 µs | 29.6 µs | 203 µs |  |  | 7.4x faster |
| `bspline1/n1000/sorted` | 23.6 µs | 24.7 µs | 33.7 µs |  |  | 1.4x faster |
| `bspline1/n1000/shuffled` | 60.1 µs | 62.1 µs | 585 µs |  |  | 9.7x faster |
| `bspline2/n1000/sorted` | 31.6 µs | 32.2 µs | 85 µs |  |  | 2.7x faster |
| `bspline2/n1000/shuffled` | 58.7 µs | 60 µs | 534 µs |  |  | 9.1x faster |
| `bspline3/n1000/sorted` | 35.9 µs | 36.7 µs | 128 µs |  |  | 3.6x faster |
| `bspline3/n1000/shuffled` | 79.1 µs | 86 µs | 580 µs |  |  | 7.3x faster |
| `cubic_spline/n1000/sorted` | 14.3 µs | 14.6 µs | 37.5 µs | 1.48 ms | 34.8 µs | 2.4x faster |
| `cubic_spline/n1000/shuffled` | 19.7 µs | 21.8 µs | 219 µs | 1.55 ms | 31.8 µs | 1.6x faster |
| `pchip/n1000/sorted` | 14.5 µs | 14.7 µs | 38.8 µs | 1.31 ms |  | 2.7x faster |
| `pchip/n1000/shuffled` | 19.9 µs | 20.3 µs | 237 µs | 1.2 ms |  | 12x faster |
| `akima/n1000/sorted` | 14.3 µs | 15.8 µs | 37.2 µs | 1.42 ms |  | 2.6x faster |
| `akima/n1000/shuffled` | 19.8 µs | 21.9 µs | 219 µs | 1.42 ms |  | 11x faster |
| `makima/n1000/sorted` | 14.4 µs | 16 µs | 37.5 µs |  |  | 2.6x faster |
| `makima/n1000/shuffled` | 21.6 µs | 20.6 µs | 218 µs |  |  | 10x faster |
| `linear/n10000/sorted` | 140 µs | 60.8 µs | 244 µs | 1.65 ms | 247 µs | 1.7x faster |
| `linear/n10000/shuffled` | 190 µs | 84.2 µs | 2.33 ms | 3.02 ms | 283 µs | 1.5x faster |
| `nearest/n10000/sorted` | 155 µs | 71.6 µs | 594 µs | 1.89 ms |  | 3.8x faster |
| `nearest/n10000/shuffled` | 379 µs | 92.3 µs | 2.96 ms | 3.35 ms |  | 7.8x faster |
| `previous/n10000/sorted` | 155 µs | 73.2 µs | 632 µs |  |  | 4.1x faster |
| `previous/n10000/shuffled` | 374 µs | 88.8 µs | 3.2 ms |  |  | 8.5x faster |
| `bspline1/n10000/sorted` | 263 µs | 97.7 µs | 298 µs |  |  | 1.1x faster |
| `bspline1/n10000/shuffled` | 700 µs | 146 µs | 51.7 ms |  |  | 74x faster |
| `bspline2/n10000/sorted` | 359 µs | 91.4 µs | 824 µs |  |  | 2.3x faster |
| `bspline2/n10000/shuffled` | 973 µs | 196 µs | 54.7 ms |  |  | 56x faster |
| `bspline3/n10000/sorted` | 398 µs | 112 µs | 1.38 ms |  |  | 3.5x faster |
| `bspline3/n10000/shuffled` | 1.12 ms | 248 µs | 51 ms |  |  | 45x faster |
| `cubic_spline/n10000/sorted` | 162 µs | 70.5 µs | 356 µs | 154 ms | 320 µs | 2x faster |
| `cubic_spline/n10000/shuffled` | 218 µs | 96.1 µs | 3.68 ms | 148 ms | 332 µs | 1.5x faster |
| `pchip/n10000/sorted` | 153 µs | 70.2 µs | 356 µs | 126 ms |  | 2.3x faster |
| `pchip/n10000/shuffled` | 220 µs | 91.4 µs | 3.71 ms | 127 ms |  | 17x faster |
| `akima/n10000/sorted` | 153 µs | 70.8 µs | 364 µs | 149 ms |  | 2.4x faster |
| `akima/n10000/shuffled` | 217 µs | 94.1 µs | 3.47 ms | 145 ms |  | 16x faster |
| `makima/n10000/sorted` | 153 µs | 71.5 µs | 358 µs |  |  | 2.3x faster |
| `makima/n10000/shuffled` | 239 µs | 96.2 µs | 3.68 ms |  |  | 15x faster |
| `linear/n100/sorted-uneven` | 1.17 µs | 1.29 µs | 7.35 µs | 3.69 µs | 3.77 µs | 3.1x faster |
| `linear/n100/shuffled-uneven` | 2.73 µs | 3.12 µs | 12.5 µs | 4.25 µs | 4.05 µs | 1.5x faster |
| `cubic_spline/n100/sorted-uneven` | 1.34 µs | 1.6 µs | 8.07 µs | 19.4 µs | 4.85 µs | 3.6x faster |
| `cubic_spline/n100/shuffled-uneven` | 2.9 µs | 3.07 µs | 8.06 µs | 17.8 µs | 4.41 µs | 1.5x faster |
| `pchip/n100/sorted-uneven` | 1.35 µs | 1.46 µs | 7.51 µs | 17.6 µs |  | 5.5x faster |
| `pchip/n100/shuffled-uneven` | 2.95 µs | 3.33 µs | 8.05 µs | 17.6 µs |  | 2.7x faster |
| `linear/n10000/sorted-uneven` | 214 µs | 84.1 µs | 262 µs | 1.64 ms | 970 µs | 1.2x faster |
| `linear/n10000/shuffled-uneven` | 679 µs | 170 µs | 2.53 ms | 3.13 ms | 1.75 ms | 2.6x faster |
| `cubic_spline/n10000/sorted-uneven` | 222 µs | 95.5 µs | 439 µs | 153 ms | 1.05 ms | 2x faster |
| `cubic_spline/n10000/shuffled-uneven` | 797 µs | 167 µs | 3.77 ms | 146 ms | 1.82 ms | 2.3x faster |
| `pchip/n10000/sorted-uneven` | 216 µs | 92.9 µs | 444 µs | 122 ms |  | 2.1x faster |
| `pchip/n10000/shuffled-uneven` | 807 µs | 168 µs | 3.54 ms | 123 ms |  | 4.4x faster |

### Construction

| case | interpx | interpx, parallel | SciPy | scirs2-interpolate | ndarray-interp | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/n10/sorted` | 96 ns | 98.4 ns | 9.58 µs | 80.1 ns | 42.6 ns | **2.3x slower** |
| `linear/n10/shuffled` | 93.8 ns | 101 ns | 9.46 µs | 89.1 ns | 41.5 ns | **2.3x slower** |
| `nearest/n10/sorted` | 76.3 ns | 81.7 ns | 9.85 µs | 82.1 ns |  | 1.1x faster |
| `nearest/n10/shuffled` | 77.9 ns | 75.2 ns | 10.4 µs | 81.3 ns |  | 1x faster |
| `previous/n10/sorted` | 84.1 ns | 76.6 ns | 9.24 µs |  |  | 1.1e+02x faster |
| `previous/n10/shuffled` | 76.9 ns | 75.1 ns | 10.1 µs |  |  | 1.3e+02x faster |
| `bspline1/n10/sorted` | 112 ns | 103 ns | 15.7 µs |  |  | 1.4e+02x faster |
| `bspline1/n10/shuffled` | 112 ns | 104 ns | 15.8 µs |  |  | 1.4e+02x faster |
| `bspline2/n10/sorted` | 468 ns | 466 ns | 38.9 µs |  |  | 83x faster |
| `bspline2/n10/shuffled` | 468 ns | 469 ns | 38 µs |  |  | 81x faster |
| `bspline3/n10/sorted` | 525 ns | 557 ns | 41 µs |  |  | 78x faster |
| `bspline3/n10/shuffled` | 519 ns | 511 ns | 37.1 µs |  |  | 71x faster |
| `cubic_spline/n10/sorted` | 235 ns | 242 ns | 88 µs | 395 ns | 465 ns | 1.7x faster |
| `cubic_spline/n10/shuffled` | 237 ns | 261 ns | 89.1 µs | 388 ns | 474 ns | 1.6x faster |
| `pchip/n10/sorted` | 179 ns | 178 ns | 107 µs | 212 ns |  | 1.2x faster |
| `pchip/n10/shuffled` | 194 ns | 180 ns | 103 µs | 194 ns |  | **1x slower** |
| `akima/n10/sorted` | 188 ns | 186 ns | 77.9 µs | 179 ns |  | **1x slower** |
| `akima/n10/shuffled` | 200 ns | 200 ns | 77.6 µs | 196 ns |  | **1x slower** |
| `makima/n10/sorted` | 211 ns | 212 ns | 79.9 µs |  |  | 3.8e+02x faster |
| `makima/n10/shuffled` | 210 ns | 211 ns | 80.3 µs |  |  | 3.8e+02x faster |
| `linear/n50/sorted` | 156 ns | 158 ns | 8.81 µs | 121 ns | 99.2 ns | **1.6x slower** |
| `linear/n50/shuffled` | 154 ns | 160 ns | 9.62 µs | 131 ns | 100 ns | **1.5x slower** |
| `nearest/n50/sorted` | 107 ns | 109 ns | 10.3 µs | 123 ns |  | 1.2x faster |
| `nearest/n50/shuffled` | 108 ns | 110 ns | 9.79 µs | 127 ns |  | 1.2x faster |
| `previous/n50/sorted` | 117 ns | 119 ns | 10.4 µs |  |  | 90x faster |
| `previous/n50/shuffled` | 107 ns | 121 ns | 9.61 µs |  |  | 90x faster |
| `bspline1/n50/sorted` | 157 ns | 153 ns | 17.2 µs |  |  | 1.1e+02x faster |
| `bspline1/n50/shuffled` | 158 ns | 161 ns | 16 µs |  |  | 1e+02x faster |
| `bspline2/n50/sorted` | 2 µs | 1.85 µs | 42.4 µs |  |  | 21x faster |
| `bspline2/n50/shuffled` | 2 µs | 1.95 µs | 41.6 µs |  |  | 21x faster |
| `bspline3/n50/sorted` | 2.05 µs | 2.06 µs | 40.4 µs |  |  | 20x faster |
| `bspline3/n50/shuffled` | 2.23 µs | 2.06 µs | 40.7 µs |  |  | 18x faster |
| `cubic_spline/n50/sorted` | 918 ns | 840 ns | 96.8 µs | 946 ns | 1.08 µs | 1x faster |
| `cubic_spline/n50/shuffled` | 838 ns | 838 ns | 89 µs | 872 ns | 1.08 µs | 1x faster |
| `pchip/n50/sorted` | 545 ns | 503 ns | 105 µs | 334 ns |  | **1.6x slower** |
| `pchip/n50/shuffled` | 545 ns | 498 ns | 103 µs | 344 ns |  | **1.6x slower** |
| `akima/n50/sorted` | 502 ns | 507 ns | 77.8 µs | 327 ns |  | **1.5x slower** |
| `akima/n50/shuffled` | 504 ns | 505 ns | 77.5 µs | 334 ns |  | **1.5x slower** |
| `makima/n50/sorted` | 597 ns | 598 ns | 80.1 µs |  |  | 1.3e+02x faster |
| `makima/n50/shuffled` | 598 ns | 566 ns | 85.7 µs |  |  | 1.4e+02x faster |
| `linear/n100/sorted` | 271 ns | 268 ns | 9.9 µs | 207 ns | 167 ns | **1.6x slower** |
| `linear/n100/shuffled` | 270 ns | 271 ns | 9.17 µs | 196 ns | 166 ns | **1.6x slower** |
| `nearest/n100/sorted` | 176 ns | 177 ns | 10.9 µs | 199 ns |  | 1.1x faster |
| `nearest/n100/shuffled` | 162 ns | 176 ns | 10 µs | 196 ns |  | 1.2x faster |
| `previous/n100/sorted` | 161 ns | 164 ns | 10.8 µs |  |  | 67x faster |
| `previous/n100/shuffled` | 176 ns | 176 ns | 10.4 µs |  |  | 59x faster |
| `bspline1/n100/sorted` | 215 ns | 215 ns | 16.2 µs |  |  | 75x faster |
| `bspline1/n100/shuffled` | 196 ns | 203 ns | 16.2 µs |  |  | 82x faster |
| `bspline2/n100/sorted` | 3.5 µs | 3.42 µs | 45.3 µs |  |  | 13x faster |
| `bspline2/n100/shuffled` | 3.53 µs | 3.6 µs | 44.6 µs |  |  | 13x faster |
| `bspline3/n100/sorted` | 3.89 µs | 4.28 µs | 45.8 µs |  |  | 12x faster |
| `bspline3/n100/shuffled` | 3.89 µs | 3.89 µs | 50.8 µs |  |  | 13x faster |
| `cubic_spline/n100/sorted` | 1.58 µs | 1.73 µs | 90.3 µs | 1.53 µs | 2.05 µs | **1x slower** |
| `cubic_spline/n100/shuffled` | 1.58 µs | 1.58 µs | 90.7 µs | 1.52 µs | 1.88 µs | **1x slower** |
| `pchip/n100/sorted` | 967 ns | 973 ns | 115 µs | 593 ns |  | **1.6x slower** |
| `pchip/n100/shuffled` | 967 ns | 969 ns | 114 µs | 587 ns |  | **1.6x slower** |
| `akima/n100/sorted` | 898 ns | 980 ns | 83.1 µs | 575 ns |  | **1.6x slower** |
| `akima/n100/shuffled` | 978 ns | 899 ns | 85.4 µs | 524 ns |  | **1.9x slower** |
| `makima/n100/sorted` | 1.07 µs | 988 ns | 89.2 µs |  |  | 83x faster |
| `makima/n100/shuffled` | 980 ns | 985 ns | 82.2 µs |  |  | 84x faster |
| `linear/n1000/sorted` | 2.12 µs | 2.12 µs | 14.3 µs | 1.33 µs | 1.51 µs | **1.6x slower** |
| `linear/n1000/shuffled` | 2.12 µs | 1.95 µs | 13.5 µs | 1.22 µs | 1.38 µs | **1.7x slower** |
| `nearest/n1000/sorted` | 1.08 µs | 1.08 µs | 15.7 µs | 1.25 µs |  | 1.2x faster |
| `nearest/n1000/shuffled` | 1.09 µs | 1.08 µs | 16.1 µs | 1.29 µs |  | 1.2x faster |
| `previous/n1000/sorted` | 994 ns | 1.08 µs | 17.1 µs |  |  | 17x faster |
| `previous/n1000/shuffled` | 1.07 µs | 994 ns | 16.6 µs |  |  | 15x faster |
| `bspline1/n1000/sorted` | 1.08 µs | 1.04 µs | 18.7 µs |  |  | 17x faster |
| `bspline1/n1000/shuffled` | 1.13 µs | 1.13 µs | 18.8 µs |  |  | 17x faster |
| `bspline2/n1000/sorted` | 34.6 µs | 35.5 µs | 111 µs |  |  | 3.2x faster |
| `bspline2/n1000/shuffled` | 37.7 µs | 38.1 µs | 111 µs |  |  | 3x faster |
| `bspline3/n1000/sorted` | 71.3 µs | 69.1 µs | 127 µs |  |  | 1.8x faster |
| `bspline3/n1000/shuffled` | 66.5 µs | 67.8 µs | 126 µs |  |  | 1.9x faster |
| `cubic_spline/n1000/sorted` | 38.6 µs | 37.5 µs | 114 µs | 18.1 µs | 17.4 µs | **2.2x slower** |
| `cubic_spline/n1000/shuffled` | 38.6 µs | 40.1 µs | 115 µs | 16.9 µs | 17.6 µs | **2.3x slower** |
| `pchip/n1000/sorted` | 12.4 µs | 11.6 µs | 132 µs | 4.56 µs |  | **2.7x slower** |
| `pchip/n1000/shuffled` | 11.7 µs | 11.9 µs | 132 µs | 4.16 µs |  | **2.8x slower** |
| `akima/n1000/sorted` | 12.2 µs | 19.7 µs | 98.8 µs | 8.8 µs |  | **1.4x slower** |
| `akima/n1000/shuffled` | 11.8 µs | 12 µs | 98.7 µs | 5.58 µs |  | **2.1x slower** |
| `makima/n1000/sorted` | 12.5 µs | 12.8 µs | 105 µs |  |  | 8.4x faster |
| `makima/n1000/shuffled` | 12.7 µs | 13.8 µs | 110 µs |  |  | 8.6x faster |
| `linear/n10000/sorted` | 93.4 µs | 81.2 µs | 50.5 µs | 61 µs | 14.1 µs | **6.6x slower** |
| `linear/n10000/shuffled` | 89.9 µs | 80.5 µs | 45.8 µs | 56 µs | 32.8 µs | **2.7x slower** |
| `nearest/n10000/sorted` | 61.1 µs | 54.4 µs | 53.5 µs | 57.9 µs |  | **1.1x slower** |
| `nearest/n10000/shuffled` | 67.9 µs | 55.3 µs | 58.9 µs | 62.7 µs |  | **1.2x slower** |
| `previous/n10000/sorted` | 63.5 µs | 54.2 µs | 75.1 µs |  |  | 1.2x faster |
| `previous/n10000/shuffled` | 66.7 µs | 56.4 µs | 75 µs |  |  | 1.1x faster |
| `bspline1/n10000/sorted` | 66 µs | 64.9 µs | 23.1 µs |  |  | **2.9x slower** |
| `bspline1/n10000/shuffled` | 69.5 µs | 60 µs | 25.2 µs |  |  | **2.8x slower** |
| `bspline2/n10000/sorted` | 545 µs | 529 µs | 772 µs |  |  | 1.4x faster |
| `bspline2/n10000/shuffled` | 532 µs | 520 µs | 756 µs |  |  | 1.4x faster |
| `bspline3/n10000/sorted` | 630 µs | 634 µs | 969 µs |  |  | 1.5x faster |
| `bspline3/n10000/shuffled` | 614 µs | 602 µs | 954 µs |  |  | 1.6x faster |
| `cubic_spline/n10000/sorted` | 251 µs | 351 µs | 785 µs | 335 µs | 363 µs | 1.3x faster |
| `cubic_spline/n10000/shuffled` | 251 µs | 310 µs | 738 µs | 357 µs | 348 µs | 1.4x faster |
| `pchip/n10000/sorted` | 121 µs | 104 µs | 246 µs | 174 µs |  | 1.4x faster |
| `pchip/n10000/shuffled` | 108 µs | 90.8 µs | 243 µs | 173 µs |  | 1.6x faster |
| `akima/n10000/sorted` | 114 µs | 85 µs | 636 µs | 143 µs |  | 1.2x faster |
| `akima/n10000/shuffled` | 96.9 µs | 203 µs | 642 µs | 143 µs |  | 1.5x faster |
| `makima/n10000/sorted` | 104 µs | 213 µs | 776 µs |  |  | 7.5x faster |
| `makima/n10000/shuffled` | 151 µs | 210 µs | 823 µs |  |  | 5.4x faster |
| `linear/n100/sorted-uneven` | 249 ns | 272 ns | 9.19 µs | 215 ns | 166 ns | **1.5x slower** |
| `linear/n100/shuffled-uneven` | 272 ns | 272 ns | 9.38 µs | 197 ns | 166 ns | **1.6x slower** |
| `cubic_spline/n100/sorted-uneven` | 1.58 µs | 1.74 µs | 99 µs | 1.57 µs | 2.05 µs | **1x slower** |
| `cubic_spline/n100/shuffled-uneven` | 1.73 µs | 1.58 µs | 90.3 µs | 1.53 µs | 1.87 µs | **1.1x slower** |
| `pchip/n100/sorted-uneven` | 963 ns | 882 ns | 105 µs | 589 ns |  | **1.6x slower** |
| `pchip/n100/shuffled-uneven` | 884 ns | 903 ns | 105 µs | 592 ns |  | **1.5x slower** |
| `linear/n10000/sorted-uneven` | 97.9 µs | 82.2 µs | 46 µs | 61.1 µs | 20.3 µs | **4.8x slower** |
| `linear/n10000/shuffled-uneven` | 98.7 µs | 92.8 µs | 47.1 µs | 61.1 µs | 19.4 µs | **5.1x slower** |
| `cubic_spline/n10000/sorted-uneven` | 248 µs | 320 µs | 748 µs | 359 µs | 389 µs | 1.4x faster |
| `cubic_spline/n10000/shuffled-uneven` | 248 µs | 336 µs | 765 µs | 297 µs | 371 µs | 1.2x faster |
| `pchip/n10000/sorted-uneven` | 100 µs | 93 µs | 246 µs | 168 µs |  | 1.7x faster |
| `pchip/n10000/shuffled-uneven` | 108 µs | 92.2 µs | 263 µs | 163 µs |  | 1.5x faster |

## `PPoly`

A cubic of 55 intervals, evaluated at points in no order.

### Evaluation

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `intervals55/eval10` | 79.3 ns | 80.2 ns | 3.53 µs | 45x faster |
| `intervals55/eval100` | 643 ns | 704 ns | 4.49 µs | 7x faster |
| `intervals55/eval1000` | 6.48 µs | 6.66 µs | 15.6 µs | 2.4x faster |
| `intervals55/eval10000` | 66.5 µs | 66.9 µs | 332 µs | 5x faster |

### Construction

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `intervals55/eval10` | 159 ns | 174 ns | 5.25 µs | 33x faster |
| `intervals55/eval100` | 160 ns | 159 ns | 5.21 µs | 33x faster |
| `intervals55/eval1000` | 159 ns | 160 ns | 5.21 µs | 33x faster |
| `intervals55/eval10000` | 174 ns | 159 ns | 5.22 µs | 30x faster |

## `BSpline`

A cubic B-spline of 55 intervals, evaluated at points in no order.

### Evaluation

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `intervals55/eval10` | 146 ns | 161 ns | 2.66 µs | 18x faster |
| `intervals55/eval100` | 1.42 µs | 1.52 µs | 5.79 µs | 4.1x faster |
| `intervals55/eval1000` | 14.2 µs | 14.3 µs | 37.4 µs | 2.6x faster |
| `intervals55/eval10000` | 142 µs | 48.5 µs | 356 µs | 2.5x faster |

### Construction

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `intervals55/eval10` | 175 ns | 161 ns | 9.3 µs | 53x faster |
| `intervals55/eval100` | 160 ns | 161 ns | 10.1 µs | 63x faster |
| `intervals55/eval1000` | 160 ns | 168 ns | 9.46 µs | 59x faster |
| `intervals55/eval10000` | 161 ns | 176 ns | 9.16 µs | 57x faster |

## `BSplineCurve`

A cubic curve through `n` points of the plane, evaluated at `4 n` parameters in order.

### Evaluation

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `degree3/n100` | 4.55 µs | 4.26 µs | 16.5 µs | 3.6x faster |
| `degree3/n10000` | 450 µs | 144 µs | 1.34 ms | 3x faster |

### Construction

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `degree3/n100` | 4.64 µs | 4.24 µs | 48.2 µs | 10x faster |
| `degree3/n10000` | 659 µs | 659 µs | 1.39 ms | 2.1x faster |

## Regular grid

Points in no order inside the grid. `maxN`: the first axis has `N` points and each axis half the one before (SciPy's suite); `eachN`: every axis has `N` points; `flipped`: every axis given from its high end.

### Evaluation

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | ndarray-interp | interpn | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/2d/max10/eval10/ascending` | 134 ns | 137 ns | 135 ns | 19.1 µs | 808 ns | 153 ns | 94.8 ns | **1.4x slower** |
| `linear/2d/max10/eval10/flipped` | 147 ns | 149 ns | 148 ns | 19 µs |  |  |  | 1.3e+02x faster |
| `linear/2d/max10/eval100/ascending` | 1.25 µs | 1.43 µs | 1.25 µs | 21.8 µs | 8.98 µs | 1.36 µs | 931 ns | **1.3x slower** |
| `linear/2d/max10/eval100/flipped` | 1.37 µs | 1.31 µs | 1.37 µs | 23.5 µs |  |  |  | 17x faster |
| `linear/2d/max10/eval1000/ascending` | 13.7 µs | 12.6 µs | 12.5 µs | 40.4 µs | 87.8 µs | 14.2 µs | 9.25 µs | **1.5x slower** |
| `linear/2d/max10/eval1000/flipped` | 12.5 µs | 12.6 µs | 13.6 µs | 40.2 µs |  |  |  | 3.2x faster |
| `linear/2d/max10/eval10000/ascending` | 128 µs | 45.8 µs | 125 µs | 361 µs | 893 µs | 219 µs | 92.3 µs | **1.4x slower** |
| `linear/2d/max10/eval10000/flipped` | 125 µs | 47 µs | 129 µs | 369 µs |  |  |  | 3x faster |
| `linear/2d/max40/eval10/ascending` | 158 ns | 175 ns | 160 ns | 20.9 µs | 1 µs | 169 ns | 134 ns | **1.2x slower** |
| `linear/2d/max40/eval10/flipped` | 157 ns | 159 ns | 157 ns | 19.1 µs |  |  |  | 1.2e+02x faster |
| `linear/2d/max40/eval100/ascending` | 1.57 µs | 1.59 µs | 1.52 µs | 21.9 µs | 11.1 µs | 1.59 µs | 1.33 µs | **1.2x slower** |
| `linear/2d/max40/eval100/flipped` | 1.52 µs | 1.74 µs | 1.67 µs | 21.8 µs |  |  |  | 14x faster |
| `linear/2d/max40/eval1000/ascending` | 15.1 µs | 16.7 µs | 15.1 µs | 46.1 µs | 106 µs | 15.2 µs | 13.2 µs | **1.1x slower** |
| `linear/2d/max40/eval1000/flipped` | 15.1 µs | 16 µs | 15.1 µs | 45.8 µs |  |  |  | 3x faster |
| `linear/2d/max40/eval10000/ascending` | 152 µs | 56.2 µs | 152 µs | 589 µs | 1.13 ms | 229 µs | 144 µs | **1x slower** |
| `linear/2d/max40/eval10000/flipped` | 152 µs | 54.6 µs | 166 µs | 578 µs |  |  |  | 3.8x faster |
| `linear/2d/max200/eval10/ascending` | 193 ns | 197 ns | 194 ns | 19.2 µs | 2 µs | 204 ns | 179 ns | **1.1x slower** |
| `linear/2d/max200/eval10/flipped` | 212 ns | 197 ns | 212 ns | 19.3 µs |  |  |  | 91x faster |
| `linear/2d/max200/eval100/ascending` | 1.89 µs | 2.15 µs | 2.06 µs | 22.6 µs | 21.1 µs | 2.18 µs | 1.72 µs | **1.1x slower** |
| `linear/2d/max200/eval100/flipped` | 2.07 µs | 2.15 µs | 2.06 µs | 25.1 µs |  |  |  | 12x faster |
| `linear/2d/max200/eval1000/ascending` | 18.9 µs | 19.2 µs | 19 µs | 79.1 µs | 207 µs | 22.9 µs | 17.2 µs | **1.1x slower** |
| `linear/2d/max200/eval1000/flipped` | 18.9 µs | 19.3 µs | 20.8 µs | 72.6 µs |  |  |  | 3.8x faster |
| `linear/2d/max200/eval10000/ascending` | 189 µs | 83.4 µs | 192 µs | 835 µs | 2.18 ms | 494 µs | 189 µs | **1x slower** |
| `linear/2d/max200/eval10000/flipped` | 207 µs | 74 µs | 207 µs | 857 µs |  |  |  | 4.1x faster |
| `linear/3d/max10/eval10/ascending` | 199 ns | 217 ns | 217 ns | 57.7 µs | 3.42 µs |  | 155 ns | **1.3x slower** |
| `linear/3d/max10/eval10/flipped` | 198 ns | 204 ns | 217 ns | 57.2 µs |  |  |  | 2.9e+02x faster |
| `linear/3d/max10/eval100/ascending` | 2.09 µs | 1.98 µs | 1.98 µs | 62.7 µs | 31.5 µs |  | 1.53 µs | **1.4x slower** |
| `linear/3d/max10/eval100/flipped` | 1.91 µs | 1.98 µs | 1.91 µs | 68.3 µs |  |  |  | 36x faster |
| `linear/3d/max10/eval1000/ascending` | 20.8 µs | 19.2 µs | 20.2 µs | 131 µs | 307 µs |  | 14.4 µs | **1.4x slower** |
| `linear/3d/max10/eval1000/flipped` | 20.9 µs | 19.4 µs | 21 µs | 118 µs |  |  |  | 5.7x faster |
| `linear/3d/max10/eval10000/ascending` | 208 µs | 68.2 µs | 210 µs | 1.14 ms | 3.19 ms |  | 155 µs | **1.3x slower** |
| `linear/3d/max10/eval10000/flipped` | 192 µs | 68.9 µs | 210 µs | 1.17 ms |  |  |  | 6.1x faster |
| `linear/3d/max40/eval10/ascending` | 226 ns | 225 ns | 223 ns | 57.4 µs | 3.64 µs |  | 199 ns | **1.1x slower** |
| `linear/3d/max40/eval10/flipped` | 223 ns | 246 ns | 241 ns | 57.6 µs |  |  |  | 2.6e+02x faster |
| `linear/3d/max40/eval100/ascending` | 2.14 µs | 2.21 µs | 2.14 µs | 63.8 µs | 36.4 µs |  | 2.17 µs | 1x faster |
| `linear/3d/max40/eval100/flipped` | 2.25 µs | 2.42 µs | 2.34 µs | 63.8 µs |  |  |  | 28x faster |
| `linear/3d/max40/eval1000/ascending` | 23.5 µs | 23.7 µs | 23.4 µs | 137 µs | 350 µs |  | 21.3 µs | **1.1x slower** |
| `linear/3d/max40/eval1000/flipped` | 23.5 µs | 21.6 µs | 23.5 µs | 147 µs |  |  |  | 6.3x faster |
| `linear/3d/max40/eval10000/ascending` | 234 µs | 77.2 µs | 214 µs | 1.63 ms | 3.88 ms |  | 200 µs | **1.2x slower** |
| `linear/3d/max40/eval10000/flipped` | 226 µs | 75.2 µs | 214 µs | 1.48 ms |  |  |  | 6.6x faster |
| `linear/3d/max200/eval10/ascending` | 278 ns | 311 ns | 277 ns | 57.3 µs | 4.5 µs |  | 286 ns | 1x faster |
| `linear/3d/max200/eval10/flipped` | 304 ns | 310 ns | 304 ns | 57.6 µs |  |  |  | 1.9e+02x faster |
| `linear/3d/max200/eval100/ascending` | 2.72 µs | 2.81 µs | 2.98 µs | 65.6 µs | 51.2 µs |  | 2.75 µs | 1x faster |
| `linear/3d/max200/eval100/flipped` | 2.82 µs | 2.81 µs | 2.74 µs | 72 µs |  |  |  | 26x faster |
| `linear/3d/max200/eval1000/ascending` | 27.3 µs | 30.2 µs | 28.6 µs | 190 µs | 478 µs |  | 26.3 µs | **1x slower** |
| `linear/3d/max200/eval1000/flipped` | 29.8 µs | 27.5 µs | 27.6 µs | 189 µs |  |  |  | 6.3x faster |
| `linear/3d/max200/eval10000/ascending` | 272 µs | 82.3 µs | 273 µs | 2 ms | 5.18 ms |  | 289 µs | 1.1x faster |
| `linear/3d/max200/eval10000/flipped` | 272 µs | 93 µs | 278 µs | 1.89 ms |  |  |  | 6.9x faster |
| `linear/4d/max10/eval10/ascending` | 260 ns | 265 ns | 264 ns | 97 µs |  |  |  | 3.7e+02x faster |
| `linear/4d/max10/eval10/flipped` | 259 ns | 278 ns | 255 ns | 96.2 µs |  |  |  | 3.7e+02x faster |
| `linear/4d/max10/eval100/ascending` | 2.41 µs | 2.47 µs | 2.41 µs | 107 µs |  |  |  | 44x faster |
| `linear/4d/max10/eval100/flipped` | 2.5 µs | 2.48 µs | 2.48 µs | 108 µs |  |  |  | 43x faster |
| `linear/4d/max10/eval1000/ascending` | 24.6 µs | 24.5 µs | 24.5 µs | 219 µs |  |  |  | 8.9x faster |
| `linear/4d/max10/eval1000/flipped` | 24.5 µs | 25.6 µs | 24.6 µs | 209 µs |  |  |  | 8.5x faster |
| `linear/4d/max10/eval10000/ascending` | 243 µs | 83.2 µs | 252 µs | 1.93 ms |  |  |  | 8x faster |
| `linear/4d/max10/eval10000/flipped` | 248 µs | 79.2 µs | 248 µs | 1.85 ms |  |  |  | 7.5x faster |
| `linear/4d/max40/eval10/ascending` | 309 ns | 338 ns | 308 ns | 97 µs | 5.83 µs |  | 272 ns | **1.1x slower** |
| `linear/4d/max40/eval10/flipped` | 308 ns | 314 ns | 309 ns | 95.9 µs |  |  |  | 3.1e+02x faster |
| `linear/4d/max40/eval100/ascending` | 2.95 µs | 3.01 µs | 3.14 µs | 119 µs | 59.1 µs |  | 2.71 µs | **1.1x slower** |
| `linear/4d/max40/eval100/flipped` | 2.96 µs | 3.07 µs | 3.2 µs | 110 µs |  |  |  | 37x faster |
| `linear/4d/max40/eval1000/ascending` | 29.2 µs | 32.4 µs | 31.9 µs | 247 µs | 583 µs |  | 27 µs | **1.1x slower** |
| `linear/4d/max40/eval1000/flipped` | 29.5 µs | 30.2 µs | 31.8 µs | 252 µs |  |  |  | 8.5x faster |
| `linear/4d/max40/eval10000/ascending` | 296 µs | 96.5 µs | 295 µs | 2.37 ms | 7.09 ms |  | 296 µs | 1x faster |
| `linear/4d/max40/eval10000/flipped` | 295 µs | 90.3 µs | 318 µs | 2.31 ms |  |  |  | 7.8x faster |
| `linear/4d/max200/eval10/ascending` | 414 ns | 383 ns | 387 ns | 106 µs | 7.71 µs |  | 365 ns | **1.1x slower** |
| `linear/4d/max200/eval10/flipped` | 383 ns | 386 ns | 384 ns | 104 µs |  |  |  | 2.7e+02x faster |
| `linear/4d/max200/eval100/ascending` | 3.97 µs | 3.91 µs | 4.11 µs | 120 µs | 79.2 µs |  | 3.63 µs | **1.1x slower** |
| `linear/4d/max200/eval100/flipped` | 3.77 µs | 3.82 µs | 4.1 µs | 118 µs |  |  |  | 31x faster |
| `linear/4d/max200/eval1000/ascending` | 41.1 µs | 41.4 µs | 41 µs | 315 µs | 724 µs |  | 36.3 µs | **1.1x slower** |
| `linear/4d/max200/eval1000/flipped` | 37.9 µs | 38.8 µs | 37.4 µs | 305 µs |  |  |  | 8.1x faster |
| `linear/4d/max200/eval10000/ascending` | 379 µs | 106 µs | 375 µs | 2.98 ms | 7.99 ms |  | 400 µs | 1.1x faster |
| `linear/4d/max200/eval10000/flipped` | 376 µs | 108 µs | 376 µs | 2.95 ms |  |  |  | 7.8x faster |
| `nearest/2d/each40/eval10000/ascending` | 151 µs | 62.7 µs | 164 µs | 691 µs | 1.61 ms |  | 136 µs | **1.1x slower** |
| `nearest/3d/each40/eval10000/ascending` | 212 µs | 82.8 µs | 213 µs | 947 µs | 2.11 ms |  | 209 µs | **1x slower** |
| `pchip/2d/each40/eval10000/ascending` | 493 µs | 129 µs | 493 µs | 1.32 s |  |  |  | 2.7e+03x faster |
| `pchip/3d/each40/eval10000/ascending` | 1.53 ms | 275 µs | 1.48 ms | 3.26 s |  |  |  | 2.1e+03x faster |
| `cubic/2d/each40/eval10000/ascending` | 305 µs | 94.1 µs | 301 µs | 1.27 ms |  |  | 409 µs | 1.3x faster |
| `cubic/3d/each40/eval10000/ascending` | 726 µs | 171 µs | 660 µs | 3.06 ms |  |  | 835 µs | 1.2x faster |
| `cubic/2d/each10/eval1000/ascending` | 26.6 µs | 26.8 µs | 26.6 µs | 131 µs |  |  | 35.8 µs | 1.3x faster |
| `cubic/2d/each40/eval1000/ascending` | 32.9 µs | 32.9 µs | 30 µs | 142 µs |  |  | 40.7 µs | 1.2x faster |
| `cubic/2d/each100/eval1000/ascending` | 32.1 µs | 35.4 µs | 33.9 µs | 166 µs |  |  | 47.1 µs | 1.5x faster |
| `cubic/2d/each200/eval1000/ascending` | 34.5 µs | 37.9 µs | 34.6 µs | 205 µs |  |  | 45.8 µs | 1.3x faster |
| `cubic/2d/each400/eval1000/ascending` | 37.9 µs | 37.9 µs | 38.2 µs | 273 µs |  |  | 49.6 µs | 1.3x faster |
| `quintic/2d/each10/eval1000/ascending` | 44.2 µs | 44.4 µs | 44.2 µs | 235 µs |  |  |  | 5.3x faster |
| `quintic/2d/each40/eval1000/ascending` | 48 µs | 49.7 µs | 48.1 µs | 238 µs |  |  |  | 5x faster |
| `linear/2d/each400/eval10000/ascending` | 220 µs | 81.2 µs | 220 µs | 961 µs | 4.21 ms | 541 µs | 220 µs | **1x slower** |
| `linear/2d/each800/eval10000/ascending` | 240 µs | 83.7 µs | 262 µs | 1.08 ms | 7.57 ms | 576 µs | 242 µs | 1x faster |
| `nearest/2d/auckland-ba32-0401/heldout/ascending` | 373 ms | 62.2 ms | 385 ms | 1.39 s | 161 s |  | 760 ms | 2x faster |
| `linear/2d/auckland-ba32-0401/heldout/ascending` | 389 ms | 65.4 ms | 421 ms | 1.16 s | 70.4 s | 374 ms | 1.57 s | **1x slower** |
| `cubic/2d/auckland-ba32-0401/heldout/ascending` | 843 ms | 129 ms | 779 ms | 32.9 s |  |  | 1.9 s | 2.2x faster |
| `quintic/2d/auckland-ba32-0401/heldout/ascending` | 1.31 s | 212 ms | 1.21 s | 35.8 s |  |  |  | 27x faster |
| `pchip/2d/auckland-ba32-0401/heldout/ascending` | 1 s | 185 ms | 1.04 s |  |  |  |  |  |

### Construction

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | ndarray-interp | interpn | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/2d/max10/eval10/ascending` | 122 ns | 133 ns | 139 ns | 5.85 µs | 139 ns | 78 ns | 0 | **1.6e+02x slower** |
| `linear/2d/max10/eval10/flipped` | 168 ns | 168 ns | 174 ns | 11.5 µs |  |  |  | 69x faster |
| `linear/2d/max10/eval100/ascending` | 122 ns | 133 ns | 129 ns | 5.85 µs | 144 ns | 77.9 ns | 0 | **1.6e+02x slower** |
| `linear/2d/max10/eval100/flipped` | 167 ns | 155 ns | 178 ns | 12.7 µs |  |  |  | 76x faster |
| `linear/2d/max10/eval1000/ascending` | 132 ns | 123 ns | 129 ns | 5.73 µs | 150 ns | 78.1 ns | 0 | **1.9e+02x slower** |
| `linear/2d/max10/eval1000/flipped` | 152 ns | 154 ns | 177 ns | 11.5 µs |  |  |  | 76x faster |
| `linear/2d/max10/eval10000/ascending` | 133 ns | 121 ns | 130 ns | 5.8 µs | 143 ns | 85.7 ns | 0 | **1.9e+02x slower** |
| `linear/2d/max10/eval10000/flipped` | 169 ns | 169 ns | 166 ns | 11.7 µs |  |  |  | 69x faster |
| `linear/2d/max40/eval10/ascending` | 360 ns | 390 ns | 728 ns | 5.87 µs | 215 ns | 167 ns | 0 | **5.3e+02x slower** |
| `linear/2d/max40/eval10/flipped` | 645 ns | 659 ns | 1.2 µs | 11.6 µs |  |  |  | 18x faster |
| `linear/2d/max40/eval100/ascending` | 389 ns | 358 ns | 745 ns | 5.83 µs | 234 ns | 154 ns | 0 | **5e+02x slower** |
| `linear/2d/max40/eval100/flipped` | 644 ns | 707 ns | 1.2 µs | 11.7 µs |  |  |  | 18x faster |
| `linear/2d/max40/eval1000/ascending` | 359 ns | 395 ns | 709 ns | 5.76 µs | 217 ns | 152 ns | 0 | **5.3e+02x slower** |
| `linear/2d/max40/eval1000/flipped` | 651 ns | 655 ns | 1.1 µs | 11.6 µs |  |  |  | 18x faster |
| `linear/2d/max40/eval10000/ascending` | 359 ns | 356 ns | 781 ns | 5.79 µs | 216 ns | 167 ns | 0 | **4.6e+02x slower** |
| `linear/2d/max40/eval10000/flipped` | 666 ns | 719 ns | 1.19 µs | 11.6 µs |  |  |  | 17x faster |
| `linear/2d/max200/eval10/ascending` | 6.84 µs | 6.97 µs | 17.6 µs | 5.85 µs | 3.06 µs | 518 ns | 0 | **8.8e+03x slower** |
| `linear/2d/max200/eval10/flipped` | 14.7 µs | 14.6 µs | 26.5 µs | 11.9 µs |  |  |  | **1.2x slower** |
| `linear/2d/max200/eval100/ascending` | 6.85 µs | 7.58 µs | 17.6 µs | 6.35 µs | 3.23 µs | 518 ns | 0 | **8.8e+03x slower** |
| `linear/2d/max200/eval100/flipped` | 15.9 µs | 16 µs | 25.3 µs | 13.2 µs |  |  |  | **1.2x slower** |
| `linear/2d/max200/eval1000/ascending` | 7.42 µs | 6.98 µs | 17.8 µs | 5.76 µs | 2.97 µs | 462 ns | 0 | **1.1e+04x slower** |
| `linear/2d/max200/eval1000/flipped` | 14.3 µs | 15.9 µs | 27.3 µs | 11.9 µs |  |  |  | **1.2x slower** |
| `linear/2d/max200/eval10000/ascending` | 56.3 µs | 51.1 µs | 18.8 µs | 5.81 µs | 54.4 µs | 481 ns | 0 | **7.2e+04x slower** |
| `linear/2d/max200/eval10000/flipped` | 55.7 µs | 51.9 µs | 28.8 µs | 12.2 µs |  |  |  | **4.5x slower** |
| `linear/3d/max10/eval10/ascending` | 179 ns | 192 ns | 228 ns | 7.61 µs | 190 ns |  | 0 | **1.8e+02x slower** |
| `linear/3d/max10/eval10/flipped` | 297 ns | 325 ns | 338 ns | 15.4 µs |  |  |  | 52x faster |
| `linear/3d/max10/eval100/ascending` | 197 ns | 195 ns | 221 ns | 7.5 µs | 175 ns |  | 0 | **2e+02x slower** |
| `linear/3d/max10/eval100/flipped` | 323 ns | 301 ns | 376 ns | 15.8 µs |  |  |  | 49x faster |
| `linear/3d/max10/eval1000/ascending` | 195 ns | 197 ns | 231 ns | 7.58 µs | 174 ns |  | 0 | **2.2e+02x slower** |
| `linear/3d/max10/eval1000/flipped` | 326 ns | 326 ns | 369 ns | 16 µs |  |  |  | 49x faster |
| `linear/3d/max10/eval10000/ascending` | 189 ns | 199 ns | 227 ns | 7.48 µs | 175 ns |  | 1.07 ns | **1.8e+02x slower** |
| `linear/3d/max10/eval10000/flipped` | 325 ns | 299 ns | 367 ns | 16.2 µs |  |  |  | 50x faster |
| `linear/3d/max40/eval10/ascending` | 2.84 µs | 6.53 µs | 7.38 µs | 7.59 µs | 13.5 µs |  | 0 | **2.9e+03x slower** |
| `linear/3d/max40/eval10/flipped` | 8.84 µs | 8.39 µs | 12.6 µs | 15.6 µs |  |  |  | 1.8x faster |
| `linear/3d/max40/eval100/ascending` | 3.8 µs | 3.9 µs | 7.32 µs | 7.57 µs | 13.6 µs |  | 0 | **4.3e+03x slower** |
| `linear/3d/max40/eval100/flipped` | 8.91 µs | 8.45 µs | 13.3 µs | 15.9 µs |  |  |  | 1.8x faster |
| `linear/3d/max40/eval1000/ascending` | 15.4 µs | 16.8 µs | 7.91 µs | 8.2 µs | 13.4 µs |  | 0 | **1.8e+04x slower** |
| `linear/3d/max40/eval1000/flipped` | 30.8 µs | 29.1 µs | 13.7 µs | 15.7 µs |  |  |  | **2x slower** |
| `linear/3d/max40/eval10000/ascending` | 16.9 µs | 16.7 µs | 7.39 µs | 8.36 µs | 9.46 µs |  | 0 | **1.9e+04x slower** |
| `linear/3d/max40/eval10000/flipped` | 27.2 µs | 27.7 µs | 12.4 µs | 17.3 µs |  |  |  | **1.6x slower** |
| `linear/3d/max200/eval10/ascending` | 1.65 ms | 1.84 ms | 616 µs | 7.6 µs | 1.26 ms |  | 0 | **1.9e+06x slower** |
| `linear/3d/max200/eval10/flipped` | 3.08 ms | 3.15 ms | 2.02 ms | 17.7 µs |  |  |  | **1.7e+02x slower** |
| `linear/3d/max200/eval100/ascending` | 1.67 ms | 1.62 ms | 628 µs | 7.6 µs | 1.3 ms |  | 1.07 ns | **1.6e+06x slower** |
| `linear/3d/max200/eval100/flipped` | 3.14 ms | 3.15 ms | 1.87 ms | 16.3 µs |  |  |  | **1.9e+02x slower** |
| `linear/3d/max200/eval1000/ascending` | 1.62 ms | 2.02 ms | 622 µs | 8.3 µs | 1.2 ms |  | 0 | **1.8e+06x slower** |
| `linear/3d/max200/eval1000/flipped` | 3.17 ms | 3.26 ms | 1.93 ms | 15.9 µs |  |  |  | **2e+02x slower** |
| `linear/3d/max200/eval10000/ascending` | 1.75 ms | 1.82 ms | 619 µs | 7.87 µs | 1.29 ms |  | 0 | **1.8e+06x slower** |
| `linear/3d/max200/eval10000/flipped` | 2.99 ms | 3.42 ms | 1.91 ms | 15.9 µs |  |  |  | **1.9e+02x slower** |
| `linear/4d/max10/eval10/ascending` | 216 ns | 217 ns | 240 ns | 9.69 µs |  |  |  | 45x faster |
| `linear/4d/max10/eval10/flipped` | 341 ns | 345 ns | 362 ns | 17.8 µs |  |  |  | 52x faster |
| `linear/4d/max10/eval100/ascending` | 216 ns | 218 ns | 229 ns | 9.54 µs |  |  |  | 44x faster |
| `linear/4d/max10/eval100/flipped` | 375 ns | 343 ns | 374 ns | 17.7 µs |  |  |  | 47x faster |
| `linear/4d/max10/eval1000/ascending` | 216 ns | 220 ns | 250 ns | 9.86 µs |  |  |  | 46x faster |
| `linear/4d/max10/eval1000/flipped` | 343 ns | 343 ns | 361 ns | 18 µs |  |  |  | 53x faster |
| `linear/4d/max10/eval10000/ascending` | 216 ns | 220 ns | 239 ns | 10.4 µs |  |  |  | 48x faster |
| `linear/4d/max10/eval10000/flipped` | 347 ns | 369 ns | 341 ns | 17.8 µs |  |  |  | 51x faster |
| `linear/4d/max40/eval10/ascending` | 13.1 µs | 13.4 µs | 32 µs | 9.22 µs | 5.51 µs |  | 1.27 ns | **1e+04x slower** |
| `linear/4d/max40/eval10/flipped` | 51.9 µs | 56.6 µs | 65.4 µs | 19.8 µs |  |  |  | **2.6x slower** |
| `linear/4d/max40/eval100/ascending` | 13.1 µs | 13.7 µs | 33.2 µs | 9.35 µs | 5.54 µs |  | 1.27 ns | **1e+04x slower** |
| `linear/4d/max40/eval100/flipped` | 52.1 µs | 53.3 µs | 67.4 µs | 20.4 µs |  |  |  | **2.6x slower** |
| `linear/4d/max40/eval1000/ascending` | 13.5 µs | 13.4 µs | 28.2 µs | 9.28 µs | 5.51 µs |  | 1.27 ns | **1.1e+04x slower** |
| `linear/4d/max40/eval1000/flipped` | 56.2 µs | 56.5 µs | 62.8 µs | 19.8 µs |  |  |  | **2.8x slower** |
| `linear/4d/max40/eval10000/ascending` | 13.1 µs | 89.3 µs | 30.2 µs | 9.35 µs | 5.51 µs |  | 1.27 ns | **1e+04x slower** |
| `linear/4d/max40/eval10000/flipped` | 56.4 µs | 127 µs | 76.8 µs | 20 µs |  |  |  | **2.8x slower** |
| `linear/4d/max200/eval10/ascending` | 46.9 ms | 47.7 ms | 15.9 ms | 9.31 µs | 41.8 ms |  | 1.27 ns | **3.7e+07x slower** |
| `linear/4d/max200/eval10/flipped` | 114 ms | 105 ms | 77.9 ms | 20.1 µs |  |  |  | **5.7e+03x slower** |
| `linear/4d/max200/eval100/ascending` | 53.1 ms | 48.1 ms | 16.4 ms | 9.4 µs | 40.4 ms |  | 1.27 ns | **4.2e+07x slower** |
| `linear/4d/max200/eval100/flipped` | 104 ms | 104 ms | 72.2 ms | 21.9 µs |  |  |  | **4.8e+03x slower** |
| `linear/4d/max200/eval1000/ascending` | 48.4 ms | 47.5 ms | 15.7 ms | 9.34 µs | 39.6 ms |  | 1.27 ns | **3.8e+07x slower** |
| `linear/4d/max200/eval1000/flipped` | 104 ms | 101 ms | 71.5 ms | 20.3 µs |  |  |  | **5.1e+03x slower** |
| `linear/4d/max200/eval10000/ascending` | 47.4 ms | 46.6 ms | 17.2 ms | 9.46 µs | 38.8 ms |  | 1.27 ns | **3.7e+07x slower** |
| `linear/4d/max200/eval10000/flipped` | 102 ms | 103 ms | 70.5 ms | 20.1 µs |  |  |  | **5.1e+03x slower** |
| `nearest/2d/each40/eval10000/ascending` | 591 ns | 630 ns | 1.43 µs | 5.87 µs | 290 ns |  | 0 | **8.6e+02x slower** |
| `nearest/3d/each40/eval10000/ascending` | 20.9 µs | 21.4 µs | 50.6 µs | 7.46 µs | 9.19 µs |  | 0 | **2.4e+04x slower** |
| `pchip/2d/each40/eval10000/ascending` | 6.49 µs | 6.49 µs | 6.88 µs | 6.58 µs |  |  |  | 1x faster |
| `pchip/3d/each40/eval10000/ascending` | 501 µs | 587 µs | 488 µs | 8.35 µs |  |  |  | **60x slower** |
| `cubic/2d/each40/eval10000/ascending` | 8.37 µs | 8.98 µs | 8.16 µs | 1.79 ms |  |  | 0 | **9.5e+03x slower** |
| `cubic/3d/each40/eval10000/ascending` | 450 µs | 425 µs | 347 µs | 146 ms |  |  | 1.27 ns | **3.5e+05x slower** |
| `cubic/2d/each10/eval1000/ascending` | 1.45 µs | 1.43 µs | 1.52 µs | 1.06 ms |  |  | 0 | **1.6e+03x slower** |
| `cubic/2d/each40/eval1000/ascending` | 8.95 µs | 8.89 µs | 8.71 µs | 1.92 ms |  |  | 0 | **1e+04x slower** |
| `cubic/2d/each100/eval1000/ascending` | 74.4 µs | 102 µs | 80.6 µs | 7.13 ms |  |  | 0 | **8.5e+04x slower** |
| `cubic/2d/each200/eval1000/ascending` | 170 µs | 185 µs | 219 µs | 26.7 ms |  |  | 0 | **1.9e+05x slower** |
| `cubic/2d/each400/eval1000/ascending` | 1.39 ms | 1.37 ms | 1.17 ms | 140 ms |  |  | 0 | **1.8e+06x slower** |
| `quintic/2d/each10/eval1000/ascending` | 1.94 µs | 1.9 µs | 1.98 µs | 3.73 ms |  |  |  | 1.9e+03x faster |
| `quintic/2d/each40/eval1000/ascending` | 12.2 µs | 12.2 µs | 12.2 µs | 8.42 ms |  |  |  | 6.9e+02x faster |
| `linear/2d/each400/eval10000/ascending` | 307 µs | 370 µs | 94.2 µs | 6 µs | 268 µs | 1.14 µs | 0 | **4.5e+05x slower** |
| `linear/2d/each800/eval10000/ascending` | 1.08 ms | 1.08 ms | 442 µs | 6.37 µs | 809 µs | 2.23 µs | 0 | **1.6e+06x slower** |
| `nearest/2d/auckland-ba32-0401/heldout/ascending` | 16.2 ms | 16.8 ms | 5.27 ms | 7.83 µs | 16.2 ms |  | 0 | **2.1e+07x slower** |
| `linear/2d/auckland-ba32-0401/heldout/ascending` | 17.2 ms | 16.3 ms | 5.14 ms | 7.75 µs | 14.1 ms | 8.26 µs | 0 | **2.2e+07x slower** |
| `cubic/2d/auckland-ba32-0401/heldout/ascending` | 89.5 ms | 80 ms | 66.8 ms | 14.2 s |  |  | 0 | **9.6e+07x slower** |
| `quintic/2d/auckland-ba32-0401/heldout/ascending` | 98.3 ms | 92.2 ms | 84.2 ms | 22.5 s |  |  |  | 2.3e+02x faster |
| `pchip/2d/auckland-ba32-0401/heldout/ascending` | 77 ms | 74.8 ms | 61.5 ms |  |  |  |  |  |

## `NearestNDInterpolator`

`n` points of the unit square, evaluated at the pixels of a square image.

### Evaluation

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: |
| `n10/image50` | 32.1 µs | 32.2 µs | 112 µs | 476 µs | 15x faster |
| `n10/image100` | 123 µs | 42.2 µs | 446 µs | 1.67 ms | 14x faster |
| `n10/image500` | 2.98 ms | 536 µs | 11.2 ms | 44 ms | 15x faster |
| `n100/image50` | 112 µs | 107 µs | 186 µs | 575 µs | 5.1x faster |
| `n100/image100` | 414 µs | 116 µs | 712 µs | 2.41 ms | 5.8x faster |
| `n100/image500` | 9.7 ms | 1.54 ms | 17.4 ms | 58.1 ms | 6x faster |
| `n1000/image50` | 312 µs | 116 µs | 314 µs | 798 µs | 2.6x faster |
| `n1000/image100` | 1.17 ms | 278 µs | 1.16 ms | 3.08 ms | 2.6x faster |
| `n1000/image500` | 24.2 ms | 4.04 ms | 24.4 ms | 72.9 ms | 3x faster |
| `n10000/image50` | 327 µs | 95.2 µs | 359 µs | 1.02 ms | 3.1x faster |
| `n10000/image100` | 1.37 ms | 278 µs | 1.47 ms | 4.1 ms | 3x faster |
| `n10000/image500` | 32.1 ms | 5.12 ms | 32.1 ms | 87.7 ms | 2.7x faster |
| `n100000/image50` | 435 µs | 114 µs | 486 µs | 1.51 ms | 3.5x faster |
| `n100000/image100` | 1.91 ms | 338 µs | 2.08 ms | 5.12 ms | 2.7x faster |
| `n100000/image500` | 44.7 ms | 6.8 ms | 41.8 ms | 116 ms | 2.6x faster |
| `ghcn-tmax-20250101/heldout` | 156 µs | 66.1 µs | 157 µs | 519 µs | 3.3x faster |

### Construction

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: |
| `n10/image50` | 568 ns | 527 ns | 517 ns | 9.51 µs | 17x faster |
| `n10/image100` | 527 ns | 532 ns | 523 ns | 9.48 µs | 18x faster |
| `n10/image500` | 528 ns | 579 ns | 563 ns | 9.27 µs | 18x faster |
| `n100/image50` | 3.64 µs | 3.46 µs | 3.49 µs | 15.1 µs | 4.1x faster |
| `n100/image100` | 3.34 µs | 3.4 µs | 3.35 µs | 15.1 µs | 4.5x faster |
| `n100/image500` | 3.65 µs | 3.64 µs | 3.35 µs | 15.1 µs | 4.1x faster |
| `n1000/image50` | 56.4 µs | 57.6 µs | 50.6 µs | 120 µs | 2.1x faster |
| `n1000/image100` | 54.8 µs | 51.1 µs | 42 µs | 120 µs | 2.2x faster |
| `n1000/image500` | 35.4 µs | 61 µs | 50.2 µs | 129 µs | 3.7x faster |
| `n10000/image50` | 724 µs | 681 µs | 658 µs | 1.7 ms | 2.4x faster |
| `n10000/image100` | 625 µs | 635 µs | 690 µs | 1.68 ms | 2.7x faster |
| `n10000/image500` | 665 µs | 681 µs | 658 µs | 1.68 ms | 2.5x faster |
| `n100000/image50` | 8.89 ms | 8.64 ms | 8.5 ms | 25.9 ms | 2.9x faster |
| `n100000/image100` | 8.63 ms | 8.64 ms | 7.98 ms | 26.4 ms | 3.1x faster |
| `n100000/image500` | 8.49 ms | 8.07 ms | 7.88 ms | 26.7 ms | 3.1x faster |
| `ghcn-tmax-20250101/heldout` | 697 µs | 699 µs | 659 µs | 1.7 ms | 2.4x faster |

## `RBFInterpolator`

`n` points of a square, evaluated at 1000 points in no order (`scattered`) or at the 10 000 pixels of an image (`image`). `global`: all points; `k50`: the 50 nearest.

### Evaluation

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/global/n10/scattered` | 14.4 µs | 15.1 µs | 14.4 µs | 59.3 µs | 181 µs | 4.1x faster |
| `linear/k50/n10/scattered` | 14.4 µs | 15.9 µs | 14.4 µs | 1.46 ms |  | 1e+02x faster |
| `thin_plate_spline/global/n10/scattered` | 70.9 µs | 72.2 µs | 73.6 µs | 111 µs | 297 µs | 1.6x faster |
| `thin_plate_spline/k50/n10/scattered` | 71.2 µs | 77.9 µs | 71.1 µs | 1.39 ms |  | 20x faster |
| `cubic/global/n10/scattered` | 17.3 µs | 17.3 µs | 16.7 µs | 60.9 µs | 198 µs | 3.5x faster |
| `cubic/k50/n10/scattered` | 18 µs | 17 µs | 17.6 µs | 1.49 ms |  | 83x faster |
| `quintic/global/n10/scattered` | 23.5 µs | 24.6 µs | 25.4 µs | 69.4 µs | 236 µs | 3x faster |
| `quintic/k50/n10/scattered` | 23.3 µs | 23.3 µs | 22.9 µs | 1.36 ms |  | 58x faster |
| `multiquadric/global/n10/scattered` | 28.2 µs | 27.9 µs | 27.8 µs | 54.1 µs | 236 µs | 1.9x faster |
| `multiquadric/k50/n10/scattered` | 27.8 µs | 27.9 µs | 27.8 µs | 1.41 ms |  | 51x faster |
| `inverse_multiquadric/global/n10/scattered` | 37.4 µs | 37.6 µs | 37.4 µs | 54.8 µs | 265 µs | 1.5x faster |
| `inverse_multiquadric/k50/n10/scattered` | 37.4 µs | 41.2 µs | 37.4 µs | 1.46 ms |  | 39x faster |
| `inverse_quadratic/global/n10/scattered` | 26.2 µs | 24.4 µs | 23.8 µs | 48.6 µs |  | 1.9x faster |
| `inverse_quadratic/k50/n10/scattered` | 26.2 µs | 26.3 µs | 24 µs | 1.48 ms |  | 57x faster |
| `gaussian/global/n10/scattered` | 67.2 µs | 67.6 µs | 73.7 µs | 93.8 µs | 391 µs | 1.4x faster |
| `gaussian/k50/n10/scattered` | 73.7 µs | 73.9 µs | 73.7 µs | 1.44 ms |  | 20x faster |
| `linear/global/n100/scattered` | 148 µs | 50.5 µs | 135 µs | 288 µs | 1.87 ms | 1.9x faster |
| `linear/k50/n100/scattered` | 13.8 ms | 2.58 ms | 13.9 ms | 20.7 ms |  | 1.5x faster |
| `thin_plate_spline/global/n100/scattered` | 702 µs | 183 µs | 764 µs | 870 µs | 2.95 ms | 1.2x faster |
| `thin_plate_spline/k50/n100/scattered` | 23.6 ms | 4.16 ms | 23.2 ms | 24.6 ms |  | 1x faster |
| `cubic/global/n100/scattered` | 138 µs | 54.9 µs | 139 µs | 296 µs | 1.93 ms | 2.1x faster |
| `cubic/k50/n100/scattered` | 15.3 ms | 2.81 ms | 15.3 ms | 22.6 ms |  | 1.5x faster |
| `quintic/global/n100/scattered` | 165 µs | 63.4 µs | 165 µs | 307 µs | 2.06 ms | 1.9x faster |
| `quintic/k50/n100/scattered` | 17.2 ms | 3.09 ms | 16.8 ms | 23.7 ms |  | 1.4x faster |
| `multiquadric/global/n100/scattered` | 269 µs | 76.9 µs | 295 µs | 301 µs | 2.52 ms | 1.1x faster |
| `multiquadric/k50/n100/scattered` | 15.5 ms | 2.86 ms | 16.2 ms | 19.9 ms |  | 1.3x faster |
| `inverse_multiquadric/global/n100/scattered` | 364 µs | 116 µs | 399 µs | 397 µs | 2.81 ms | 1.1x faster |
| `inverse_multiquadric/k50/n100/scattered` | 17.6 ms | 2.96 ms | 17 ms | 21.8 ms |  | 1.2x faster |
| `inverse_quadratic/global/n100/scattered` | 224 µs | 69.2 µs | 224 µs | 291 µs |  | 1.3x faster |
| `inverse_quadratic/k50/n100/scattered` | 14.7 ms | 2.6 ms | 15.3 ms | 20 ms |  | 1.4x faster |
| `gaussian/global/n100/scattered` | 677 µs | 160 µs | 675 µs | 701 µs | 3.59 ms | 1x faster |
| `gaussian/k50/n100/scattered` | 20.4 ms | 3.89 ms | 21.4 ms | 23.2 ms |  | 1.1x faster |
| `linear/global/n1000/scattered` | 1.48 ms | 321 µs | 1.48 ms | 3.98 ms | 18.7 ms | 2.7x faster |
| `linear/k50/n1000/scattered` | 14.7 ms | 2.72 ms | 14.6 ms | 37.1 ms |  | 2.5x faster |
| `thin_plate_spline/global/n1000/scattered` | 7.69 ms | 1.4 ms | 7.4 ms | 10.6 ms | 34.4 ms | 1.4x faster |
| `thin_plate_spline/k50/n1000/scattered` | 24.5 ms | 4.05 ms | 24.1 ms | 50.6 ms |  | 2.1x faster |
| `cubic/global/n1000/scattered` | 1.49 ms | 338 µs | 1.48 ms | 4.01 ms | 20.6 ms | 2.7x faster |
| `cubic/k50/n1000/scattered` | 16.2 ms | 3 ms | 16.8 ms | 40.7 ms |  | 2.5x faster |
| `quintic/global/n1000/scattered` | 1.65 ms | 387 µs | 1.66 ms | 4.07 ms | 23.1 ms | 2.5x faster |
| `quintic/k50/n1000/scattered` | 18.4 ms | 3.37 ms | 18.8 ms | 44.7 ms |  | 2.4x faster |
| `multiquadric/global/n1000/scattered` | 2.95 ms | 602 µs | 2.76 ms | 4.17 ms | 25 ms | 1.4x faster |
| `multiquadric/k50/n1000/scattered` | 18.4 ms | 3.11 ms | 18.4 ms | 40.3 ms |  | 2.2x faster |
| `inverse_multiquadric/global/n1000/scattered` | 3.84 ms | 808 µs | 4.03 ms | 5.07 ms | 27.1 ms | 1.3x faster |
| `inverse_multiquadric/k50/n1000/scattered` | 18.6 ms | 3.44 ms | 18.5 ms | 38.8 ms |  | 2.1x faster |
| `inverse_quadratic/global/n1000/scattered` | 2.46 ms | 522 µs | 2.29 ms | 3.81 ms |  | 1.5x faster |
| `inverse_quadratic/k50/n1000/scattered` | 16.1 ms | 3.2 ms | 16.6 ms | 39.5 ms |  | 2.4x faster |
| `gaussian/global/n1000/scattered` | 6.78 ms | 1.21 ms | 6.82 ms | 8.67 ms | 39.5 ms | 1.3x faster |
| `gaussian/k50/n1000/scattered` | 21.7 ms | 3.94 ms | 21.6 ms | 46.3 ms |  | 2.1x faster |
| `thin_plate_spline/global/n3000/scattered` | 22.5 ms | 3.98 ms | 22.5 ms | 30.2 ms | 93.1 ms | 1.3x faster |
| `thin_plate_spline/k50/n10000/scattered` | 24.5 ms | 4.05 ms | 23.1 ms | 50.2 ms |  | 2x faster |
| `thin_plate_spline/k50/n10000/image` | 244 ms | 37.6 ms | 258 ms | 489 ms |  | 2x faster |
| `thin_plate_spline/k50/n100000/scattered` | 23.8 ms | 4.47 ms | 24.3 ms | 50.3 ms |  | 2.1x faster |
| `thin_plate_spline/k50/n100000/image` | 242 ms | 37.6 ms | 248 ms | 498 ms |  | 2.1x faster |
| `linear/global/ghcn-tmax-20250101/heldout` | 18.3 ms | 4.28 ms | 24.4 ms | 50.1 ms | 234 ms | 2.7x faster |
| `linear/k50/ghcn-tmax-20250101/heldout` | 19.8 ms | 3.32 ms | 18.3 ms | 46.7 ms |  | 2.4x faster |
| `thin_plate_spline/global/ghcn-tmax-20250101/heldout` | 88.5 ms | 16.1 ms | 89.8 ms | 116 ms | 400 ms | 1.3x faster |
| `thin_plate_spline/k50/ghcn-tmax-20250101/heldout` | 29.3 ms | 5.18 ms | 29.9 ms | 58.6 ms |  | 2x faster |
| `cubic/global/ghcn-tmax-20250101/heldout` | 19.3 ms | 3.3 ms | 19 ms | 51 ms | 255 ms | 2.6x faster |
| `cubic/k50/ghcn-tmax-20250101/heldout` | 20.7 ms | 3.72 ms | 20.1 ms | 48.4 ms |  | 2.3x faster |
| `quintic/global/ghcn-tmax-20250101/heldout` | 21.2 ms | 4.13 ms | 19.7 ms | 48.4 ms | 284 ms | 2.3x faster |
| `quintic/k50/ghcn-tmax-20250101/heldout` | 22.1 ms | 3.86 ms | 21.9 ms | 52.9 ms |  | 2.4x faster |

### Construction

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/global/n10/scattered` | 2.82 µs | 2.59 µs | 2.59 µs | 10.8 µs | 5.86 µs | 2.1x faster |
| `linear/k50/n10/scattered` | 2.84 µs | 2.63 µs | 2.58 µs | 14.8 µs |  | 5.2x faster |
| `thin_plate_spline/global/n10/scattered` | 3.29 µs | 3.33 µs | 3.51 µs | 12.3 µs | 6.51 µs | 2x faster |
| `thin_plate_spline/k50/n10/scattered` | 3.33 µs | 3.51 µs | 3.35 µs | 15.6 µs |  | 4.7x faster |
| `cubic/global/n10/scattered` | 3.73 µs | 3.44 µs | 3.48 µs | 12.3 µs | 5.54 µs | 1.5x faster |
| `cubic/k50/n10/scattered` | 3.72 µs | 3.59 µs | 3.5 µs | 17.2 µs |  | 4.6x faster |
| `quintic/global/n10/scattered` | 3.92 µs | 3.89 µs | 4.21 µs | 15.1 µs | 5.73 µs | 1.5x faster |
| `quintic/k50/n10/scattered` | 4.25 µs | 4.28 µs | 4.24 µs | 16.9 µs |  | 4x faster |
| `multiquadric/global/n10/scattered` | 2.64 µs | 2.64 µs | 2.64 µs | 10.7 µs | 7.13 µs | 2.7x faster |
| `multiquadric/k50/n10/scattered` | 2.64 µs | 2.65 µs | 2.62 µs | 15 µs |  | 5.7x faster |
| `inverse_multiquadric/global/n10/scattered` | 2.62 µs | 2.63 µs | 2.56 µs | 10.9 µs | 6.83 µs | 2.6x faster |
| `inverse_multiquadric/k50/n10/scattered` | 2.65 µs | 2.86 µs | 2.63 µs | 15.3 µs |  | 5.8x faster |
| `inverse_quadratic/global/n10/scattered` | 2.77 µs | 2.76 µs | 2.56 µs | 10.9 µs |  | 3.9x faster |
| `inverse_quadratic/k50/n10/scattered` | 2.8 µs | 2.8 µs | 2.8 µs | 15.2 µs |  | 5.4x faster |
| `gaussian/global/n10/scattered` | 2.99 µs | 3 µs | 3 µs | 12.2 µs | 8.51 µs | 2.8x faster |
| `gaussian/k50/n10/scattered` | 3 µs | 3 µs | 3 µs | 16.3 µs |  | 5.4x faster |
| `linear/global/n100/scattered` | 75.2 µs | 77.1 µs | 79.5 µs | 83 µs | 802 µs | 1.1x faster |
| `linear/k50/n100/scattered` | 3.48 µs | 3.47 µs | 3.13 µs | 22.9 µs |  | 6.6x faster |
| `thin_plate_spline/global/n100/scattered` | 97.3 µs | 108 µs | 103 µs | 117 µs | 847 µs | 1.2x faster |
| `thin_plate_spline/k50/n100/scattered` | 3.3 µs | 3.23 µs | 3.3 µs | 23.8 µs |  | 7.2x faster |
| `cubic/global/n100/scattered` | 71.9 µs | 93.8 µs | 84.6 µs | 82.8 µs | 753 µs | 1.2x faster |
| `cubic/k50/n100/scattered` | 3.28 µs | 3.48 µs | 3.18 µs | 21.9 µs |  | 6.7x faster |
| `quintic/global/n100/scattered` | 78.1 µs | 103 µs | 89.8 µs | 93.7 µs | 764 µs | 1.2x faster |
| `quintic/k50/n100/scattered` | 3.27 µs | 3.55 µs | 3.25 µs | 23.1 µs |  | 7x faster |
| `multiquadric/global/n100/scattered` | 82.9 µs | 86.1 µs | 90.8 µs | 78.6 µs | 677 µs | **1.1x slower** |
| `multiquadric/k50/n100/scattered` | 3.47 µs | 3.42 µs | 3.45 µs | 22.6 µs |  | 6.5x faster |
| `inverse_multiquadric/global/n100/scattered` | 80.1 µs | 109 µs | 97.2 µs | 80.1 µs | 840 µs | **1x slower** |
| `inverse_multiquadric/k50/n100/scattered` | 3.17 µs | 3.49 µs | 3.27 µs | 22.7 µs |  | 7.2x faster |
| `inverse_quadratic/global/n100/scattered` | 79.4 µs | 83.1 µs | 85.7 µs | 80.8 µs |  | 1x faster |
| `inverse_quadratic/k50/n100/scattered` | 3.45 µs | 3.47 µs | 3.23 µs | 22.9 µs |  | 6.6x faster |
| `gaussian/global/n100/scattered` | 104 µs | 114 µs | 117 µs | 99.6 µs | 802 µs | **1x slower** |
| `gaussian/k50/n100/scattered` | 3.54 µs | 3.45 µs | 3.48 µs | 22.6 µs |  | 6.4x faster |
| `linear/global/n1000/scattered` | 20.7 ms | 19.8 ms | 20.7 ms | 61.8 ms | 194 ms | 3x faster |
| `linear/k50/n1000/scattered` | 63.8 µs | 58.6 µs | 58.7 µs | 150 µs |  | 2.3x faster |
| `thin_plate_spline/global/n1000/scattered` | 23.9 ms | 24 ms | 23.5 ms | 66.7 ms | 214 ms | 2.8x faster |
| `thin_plate_spline/k50/n1000/scattered` | 63.2 µs | 56.7 µs | 59.9 µs | 150 µs |  | 2.4x faster |
| `cubic/global/n1000/scattered` | 21 ms | 20.5 ms | 21 ms | 59.7 ms | 200 ms | 2.8x faster |
| `cubic/k50/n1000/scattered` | 62.3 µs | 58.8 µs | 61.7 µs | 144 µs |  | 2.3x faster |
| `quintic/global/n1000/scattered` | 21.4 ms | 21.2 ms | 22.1 ms | 59.1 ms | 213 ms | 2.8x faster |
| `quintic/k50/n1000/scattered` | 67.2 µs | 58.8 µs | 60.6 µs | 138 µs |  | 2.1x faster |
| `multiquadric/global/n1000/scattered` | 22 ms | 20.7 ms | 21.3 ms | 62 ms | 197 ms | 2.8x faster |
| `multiquadric/k50/n1000/scattered` | 63.2 µs | 59 µs | 56.1 µs | 132 µs |  | 2.1x faster |
| `inverse_multiquadric/global/n1000/scattered` | 22.3 ms | 21.5 ms | 22.1 ms | 62.3 ms | 204 ms | 2.8x faster |
| `inverse_multiquadric/k50/n1000/scattered` | 66 µs | 59 µs | 58.9 µs | 131 µs |  | 2x faster |
| `inverse_quadratic/global/n1000/scattered` | 21.6 ms | 20.7 ms | 21 ms | 60 ms |  | 2.8x faster |
| `inverse_quadratic/k50/n1000/scattered` | 67.6 µs | 59.7 µs | 59 µs | 147 µs |  | 2.2x faster |
| `gaussian/global/n1000/scattered` | 24.6 ms | 23.4 ms | 24.6 ms | 62.9 ms | 220 ms | 2.6x faster |
| `gaussian/k50/n1000/scattered` | 59.1 µs | 56.2 µs | 56 µs | 135 µs |  | 2.3x faster |
| `thin_plate_spline/global/n3000/scattered` | 444 ms | 433 ms | 443 ms | 1.94 s | 6.33 s | 4.4x faster |
| `thin_plate_spline/k50/n10000/scattered` | 789 µs | 764 µs | 688 µs | 1.67 ms |  | 2.1x faster |
| `thin_plate_spline/k50/n10000/image` | 774 µs | 813 µs | 737 µs | 1.74 ms |  | 2.2x faster |
| `thin_plate_spline/k50/n100000/scattered` | 8.7 ms | 9.16 ms | 8.81 ms | 24.6 ms |  | 2.8x faster |
| `thin_plate_spline/k50/n100000/image` | 8.95 ms | 8.77 ms | 8.4 ms | 25.4 ms |  | 2.8x faster |
| `linear/global/ghcn-tmax-20250101/heldout` | 17.8 s | 16.7 s | 17.7 s | 116 s | 285 s | 6.5x faster |
| `linear/k50/ghcn-tmax-20250101/heldout` | 774 µs | 789 µs | 741 µs | 1.79 ms |  | 2.3x faster |
| `thin_plate_spline/global/ghcn-tmax-20250101/heldout` | 18.1 s | 17.7 s | 18.6 s | 115 s | 308 s | 6.3x faster |
| `thin_plate_spline/k50/ghcn-tmax-20250101/heldout` | 788 µs | 822 µs | 680 µs | 1.76 ms |  | 2.2x faster |
| `cubic/global/ghcn-tmax-20250101/heldout` | 17.4 s | 16.5 s | 16.9 s | 117 s | 277 s | 6.7x faster |
| `cubic/k50/ghcn-tmax-20250101/heldout` | 764 µs | 810 µs | 723 µs | 1.71 ms |  | 2.2x faster |
| `quintic/global/ghcn-tmax-20250101/heldout` | 16.1 s | 16.9 s | 17.1 s | 106 s | 281 s | 6.6x faster |
| `quintic/k50/ghcn-tmax-20250101/heldout` | 819 µs | 826 µs | 754 µs | 1.75 ms |  | 2.1x faster |

## Memory

Megabytes. "Data read" is the process with the data of the case in memory, before anything is built. The columns after it are what was added to that: the highest while building, with the interpolator built, the highest while evaluating. "Highest" is the highest of the whole run, repeats included.

| case | library | data read | while building | built | while evaluating | highest |
| :--- | :--- | ---: | ---: | ---: | ---: | ---: |
| `oned/cubic_spline/n10000/sorted` | interpx | 4.8 | +0.4 | +0.5 | +0.5 | 6.0 |
| `oned/cubic_spline/n10000/sorted` | interpx, parallel | 3.3 | +2.0 | +2.0 | +2.0 | 6.4 |
| `oned/cubic_spline/n10000/sorted` | SciPy | 71.4 | +1.8 | +2.0 | +2.0 | 73.5 |
| `oned/cubic_spline/n10000/sorted` | scirs2-interpolate | 3.1 | +1.9 | +2.1 | +2.7 | 5.9 |
| `oned/cubic_spline/n10000/sorted` | ndarray-interp | 3.8 | +1.1 | +1.1 | +1.1 | 5.5 |
| `grid/linear/4d/max200/eval10000/ascending` | interpx | 230.5 | +174.6 | +174.6 | +174.6 | 405.1 |
| `grid/linear/4d/max200/eval10000/ascending` | interpx, parallel | 227.6 | +177.2 | +177.2 | +177.2 | 405.4 |
| `grid/linear/4d/max200/eval10000/ascending` | interpx, `from_owned` | 225.1 | +179.9 | +179.9 | +179.9 | 405.1 |
| `grid/linear/4d/max200/eval10000/ascending` | SciPy | 271.1 | +0.3 | +0.4 | +2.1 | 273.3 |
| `grid/linear/4d/max200/eval10000/ascending` | scirs2-interpolate | 230.3 | +171.3 | +174.6 | +174.6 | 405.0 |
| `grid/linear/4d/max200/eval10000/ascending` | interpn | 229.4 | +0.0 | +0.0 | +0.0 | 396.4 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpx | 4.3 | +1.4 | +4.0 | +4.0 | 8.5 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpx, parallel | 5.6 | +0.2 | +2.7 | +2.7 | 8.9 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpx, `from_owned` | 3.4 | +5.0 | +4.9 | +4.9 | 8.5 |
| `grid/cubic/2d/each400/eval1000/ascending` | SciPy | 72.2 | +85.2 | +6.3 | +6.3 | 161.3 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpn | 5.5 | +0.0 | +0.0 | +0.0 | 5.6 |
| `grid/linear/2d/each800/eval10000/ascending` | interpx | 12.9 | +1.9 | +1.9 | +1.9 | 15.0 |
| `grid/linear/2d/each800/eval10000/ascending` | interpx, parallel | 12.7 | +0.2 | +2.3 | +2.3 | 15.4 |
| `grid/linear/2d/each800/eval10000/ascending` | interpx, `from_owned` | 12.0 | +2.9 | +2.9 | +2.9 | 15.0 |
| `grid/linear/2d/each800/eval10000/ascending` | SciPy | 76.1 | +0.2 | +0.2 | +0.7 | 76.9 |
| `grid/linear/2d/each800/eval10000/ascending` | scirs2-interpolate | 10.9 | +1.9 | +4.0 | +4.1 | 15.0 |
| `grid/linear/2d/each800/eval10000/ascending` | ndarray-interp | 12.0 | +0.0 | +0.0 | +0.0 | 12.0 |
| `grid/linear/2d/each800/eval10000/ascending` | interpn | 12.2 | +0.0 | +0.0 | +0.0 | 12.2 |
| `nearest_nd/n100000/image500` | interpx | 12.7 | +2.6 | +1.5 | +3.5 | 18.2 |
| `nearest_nd/n100000/image500` | interpx, parallel | 14.6 | +0.6 | +0.0 | +2.1 | 18.8 |
| `nearest_nd/n100000/image500` | interpx, `from_owned` | 11.7 | +3.0 | +2.5 | +4.5 | 18.2 |
| `nearest_nd/n100000/image500` | SciPy | 77.5 | +4.6 | +4.7 | +13.3 | 91.7 |
| `rbf/thin_plate_spline/global/n3000/scattered` | interpx | 3.1 | +77.3 | +5.0 | +5.1 | 80.5 |
| `rbf/thin_plate_spline/global/n3000/scattered` | interpx, parallel | 3.4 | +77.1 | +4.8 | +5.2 | 80.9 |
| `rbf/thin_plate_spline/global/n3000/scattered` | interpx, `from_owned` | 3.4 | +77.1 | +4.7 | +4.7 | 80.5 |
| `rbf/thin_plate_spline/global/n3000/scattered` | SciPy | 71.1 | +73.0 | +0.9 | +8.9 | 144.2 |
| `rbf/thin_plate_spline/global/n3000/scattered` | scirs2-interpolate | 2.7 | +289.9 | +1.8 | +1.8 | 292.6 |
| `rbf/thin_plate_spline/k50/n100000/image` | interpx | 7.3 | +3.0 | +4.1 | +4.4 | 11.9 |
| `rbf/thin_plate_spline/k50/n100000/image` | interpx, parallel | 6.2 | +4.3 | +5.3 | +6.4 | 12.7 |
| `rbf/thin_plate_spline/k50/n100000/image` | interpx, `from_owned` | 6.4 | +6.3 | +5.1 | +5.4 | 13.5 |
| `rbf/thin_plate_spline/k50/n100000/image` | SciPy | 72.7 | +6.1 | +6.2 | +25.3 | 99.7 |

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
