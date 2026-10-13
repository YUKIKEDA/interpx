# Performance: interpx against SciPy and other Rust crates

<!-- Written by compare/perf/report.py from a run of compare/perf/run.py. Do not edit by hand: run it again. -->

Measured on 2026-10-13, at commit `a23222e501`.

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
| One variable | evaluation | 112 | 112 | 0 | 0 |
| One variable | construction | 112 | 77 | 8 | 32 |
| `PPoly` | evaluation | 4 | 4 | 0 | 0 |
| `PPoly` | construction | 4 | 4 | 0 | 0 |
| `BSpline` | evaluation | 4 | 4 | 0 | 0 |
| `BSpline` | construction | 4 | 4 | 0 | 0 |
| `BSplineCurve` | evaluation | 2 | 2 | 0 | 0 |
| `BSplineCurve` | construction | 2 | 2 | 0 | 0 |
| Regular grid | evaluation | 92 | 73 | 0 | 19 |
| Regular grid | construction | 92 | 27 | 40 | 46 |
| `NearestNDInterpolator` | evaluation | 16 | 16 | 0 | 0 |
| `NearestNDInterpolator` | construction | 16 | 16 | 0 | 0 |
| `RBFInterpolator` | evaluation | 61 | 61 | 0 | 0 |
| `RBFInterpolator` | construction | 61 | 61 | 0 | 0 |

205 comparisons (a case, an operation, and a library), by method. "Same values" says whether that library's values are SciPy's up to rounding: where they are not, it computes another interpolant under the same name, and its time is not that of the same work.

| method | operation | against | same values | cases | slower by | the largest |
| :--- | :--- | :--- | :--- | ---: | ---: | :--- |
| `grid/cubic` | construction | interpn | no | 8 | 1.67e+03x to 8.73e+07x | `grid/cubic/2d/auckland-ba32-0401/heldout/ascending`: 76.7 ms against 0 |
| `grid/linear` | construction | interpn | yes | 35 | 182x to 4.12e+07x | `grid/linear/4d/max200/eval100/ascending`: 48.3 ms against 1.17 ns |
| `grid/nearest` | construction | interpn | yes | 3 | 933x to 2.09e+07x | `grid/nearest/2d/auckland-ba32-0401/heldout/ascending`: 16.3 ms against 0 |
| `grid/linear` | construction | SciPy | yes | 37 | 1.14x to 5.23e+03x | `grid/linear/4d/max200/eval10/flipped`: 103 ms against 19.7 µs |
| `grid/nearest` | construction | SciPy | yes | 2 | 2.58x to 2.35e+03x | `grid/nearest/2d/auckland-ba32-0401/heldout/ascending`: 16.3 ms against 6.96 µs |
| `grid/linear` | construction | ndarray-interp | yes | 15 | 1.57x to 2.02e+03x | `grid/linear/2d/auckland-ba32-0401/heldout/ascending`: 16.6 ms against 8.26 µs |
| `grid/pchip` | construction | SciPy | yes | 1 | 54x | `grid/pchip/3d/each40/eval10000/ascending`: 512 µs against 9.48 µs |
| `grid/linear` | construction | scirs2-interpolate | yes | 27 | 1.15x to 5.48x | `grid/linear/3d/max40/eval10000/ascending`: 14.2 µs against 2.58 µs |
| `oned/linear` | construction | ndarray-interp | yes | 14 | 1.41x to 5.12x | `oned/linear/n10000/sorted-uneven`: 90.1 µs against 17.6 µs |
| `oned/pchip` | construction | scirs2-interpolate | yes | 8 | 1.42x to 4.29x | `oned/pchip/n1000/sorted`: 17.9 µs against 4.17 µs |
| `oned/bspline1` | construction | SciPy | yes | 2 | 2.65x to 2.95x | `oned/bspline1/n10000/sorted`: 67.2 µs against 22.8 µs |
| `grid/nearest` | construction | scirs2-interpolate | yes | 3 | 1.24x to 2.33x | `grid/nearest/2d/each40/eval10000/ascending`: 638 ns against 274 ns |
| `oned/cubic_spline` | construction | scirs2-interpolate | no | 3 | 1.13x to 2.26x | `oned/cubic_spline/n1000/shuffled`: 37.5 µs against 16.6 µs |
| `oned/linear` | construction | SciPy | yes | 4 | 1.85x to 1.98x | `oned/linear/n10000/shuffled`: 94.5 µs against 47.8 µs |
| `oned/akima` | construction | scirs2-interpolate | yes | 6 | 1.31x to 1.81x | `oned/akima/n1000/sorted`: 10 µs against 5.55 µs |
| `oned/linear` | construction | scirs2-interpolate | yes | 14 | 1.11x to 1.64x | `oned/linear/n10000/sorted`: 94 µs against 57.2 µs |
| `grid/linear` | evaluation | interpn | yes | 18 | 1.1x to 1.42x | `grid/linear/3d/max10/eval10/ascending`: 201 ns against 142 ns |
| `oned/nearest` | construction | SciPy | yes | 2 | 1.19x to 1.26x | `oned/nearest/n10000/sorted`: 67.2 µs against 53.3 µs |
| `grid/linear` | evaluation | ndarray-interp | yes | 1 | 1.14x | `grid/linear/2d/max40/eval100/ascending`: 1.66 µs against 1.45 µs |
| `grid/nearest` | evaluation | interpn | yes | 1 | 1.12x | `grid/nearest/2d/each40/eval10000/ascending`: 151 µs against 135 µs |
| `oned/nearest` | construction | scirs2-interpolate | yes | 1 | 1.11x | `oned/nearest/n10000/sorted`: 67.2 µs against 60.3 µs |

<details><summary>Every comparison where interpx is slower</summary>

"Values apart" is the distance of that library's values from SciPy's, as a share of the largest value.

| case | operation | against | interpx | the other | slower by | values apart |
| :--- | :--- | :--- | ---: | ---: | ---: | ---: |
| `grid/cubic/2d/auckland-ba32-0401/heldout/ascending` | construction | interpn | 76.7 ms | 0 | 87251740.44x | 5e-03 |
| `grid/linear/4d/max200/eval100/ascending` | construction | interpn | 48.3 ms | 1.17 ns | 41182208.00x | 3e-16 |
| `grid/linear/4d/max200/eval1000/ascending` | construction | interpn | 45.1 ms | 1.17 ns | 38482602.67x | 3e-16 |
| `grid/linear/4d/max200/eval10000/ascending` | construction | interpn | 48.1 ms | 1.27 ns | 37904856.62x | 4e-16 |
| `grid/linear/4d/max200/eval10/ascending` | construction | interpn | 44.4 ms | 1.17 ns | 37893120.00x | 2e-16 |
| `grid/linear/2d/auckland-ba32-0401/heldout/ascending` | construction | interpn | 16.6 ms | 0 | 24351890.29x | 0e+00 |
| `grid/nearest/2d/auckland-ba32-0401/heldout/ascending` | construction | interpn | 16.3 ms | 0 | 20924672.00x | 0e+00 |
| `grid/linear/3d/max200/eval100/ascending` | construction | interpn | 1.76 ms | 0 | 1997596.44x | 3e-16 |
| `grid/linear/3d/max200/eval10/ascending` | construction | interpn | 1.73 ms | 0 | 1972224.00x | 3e-16 |
| `grid/linear/3d/max200/eval1000/ascending` | construction | interpn | 1.61 ms | 0 | 1836487.11x | 3e-16 |
| `grid/linear/3d/max200/eval10000/ascending` | construction | interpn | 1.61 ms | 0 | 1826190.22x | 5e-16 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | interpn | 1.08 ms | 0 | 1577691.43x | 3e-16 |
| `grid/cubic/2d/each400/eval1000/ascending` | construction | interpn | 1.36 ms | 0 | 1549767.11x | 5e-01 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | interpn | 330 µs | 0 | 483181.71x | 3e-16 |
| `grid/cubic/3d/each40/eval10000/ascending` | construction | interpn | 398 µs | 1.17 ns | 339328.00x | 5e-01 |
| `grid/cubic/2d/each200/eval1000/ascending` | construction | interpn | 178 µs | 0 | 202296.89x | 2e-01 |
| `grid/cubic/2d/each100/eval1000/ascending` | construction | interpn | 83.4 µs | 0 | 94876.44x | 2e-01 |
| `grid/linear/2d/max200/eval10000/ascending` | construction | interpn | 44.8 µs | 0 | 65472.00x | 3e-16 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | interpn | 21.1 µs | 0 | 21555.20x | 0e+00 |
| `grid/linear/3d/max40/eval10000/ascending` | construction | interpn | 14.2 µs | 0 | 16103.11x | 3e-16 |
| `grid/linear/3d/max40/eval1000/ascending` | construction | interpn | 14.1 µs | 0 | 16028.44x | 2e-16 |
| `grid/linear/4d/max40/eval1000/ascending` | construction | interpn | 13.7 µs | 1.17 ns | 11677.33x | 6e-16 |
| `grid/linear/4d/max40/eval10000/ascending` | construction | interpn | 13.1 µs | 1.17 ns | 11194.67x | 5e-16 |
| `grid/linear/4d/max40/eval10/ascending` | construction | interpn | 13.1 µs | 1.17 ns | 11152.00x | 2e-16 |
| `grid/linear/4d/max40/eval100/ascending` | construction | interpn | 13.1 µs | 1.17 ns | 11136.00x | 2e-16 |
| `grid/cubic/2d/each40/eval1000/ascending` | construction | interpn | 9.14 µs | 0 | 10400.00x | 7e-01 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | interpn | 6.85 µs | 0 | 10020.57x | 2e-16 |
| `grid/linear/2d/max200/eval100/ascending` | construction | interpn | 6.82 µs | 0 | 9974.86x | 1e-16 |
| `grid/cubic/2d/each40/eval10000/ascending` | construction | interpn | 8.55 µs | 0 | 9728.00x | 7e-01 |
| `grid/linear/2d/max200/eval10/ascending` | construction | interpn | 6.86 µs | 0 | 8776.00x | 2e-16 |
| `grid/linear/4d/max200/eval10/flipped` | construction | SciPy | 103 ms | 19.7 µs | 5234.15x |  |
| `grid/linear/4d/max200/eval10000/ascending` | construction | SciPy | 48.1 ms | 9.3 µs | 5175.21x |  |
| `grid/linear/4d/max200/eval100/ascending` | construction | SciPy | 48.3 ms | 9.37 µs | 5151.21x |  |
| `grid/linear/3d/max40/eval100/ascending` | construction | interpn | 4.93 µs | 0 | 5043.20x | 2e-16 |
| `grid/linear/4d/max200/eval100/flipped` | construction | SciPy | 99.2 ms | 20 µs | 4955.10x |  |
| `grid/linear/4d/max200/eval10000/flipped` | construction | SciPy | 99.9 ms | 20.2 µs | 4940.27x |  |
| `grid/linear/4d/max200/eval1000/flipped` | construction | SciPy | 98.1 ms | 20.1 µs | 4888.96x |  |
| `grid/linear/4d/max200/eval1000/ascending` | construction | SciPy | 45.1 ms | 9.36 µs | 4819.97x |  |
| `grid/linear/4d/max200/eval10/ascending` | construction | SciPy | 44.4 ms | 9.42 µs | 4713.07x |  |
| `grid/linear/3d/max40/eval10/ascending` | construction | interpn | 2.77 µs | 0 | 2832.40x | 2e-16 |
| `grid/linear/2d/auckland-ba32-0401/heldout/ascending` | construction | SciPy | 16.6 ms | 6.92 µs | 2404.41x |  |
| `grid/nearest/2d/auckland-ba32-0401/heldout/ascending` | construction | SciPy | 16.3 ms | 6.96 µs | 2348.98x |  |
| `grid/linear/2d/auckland-ba32-0401/heldout/ascending` | construction | ndarray-interp | 16.6 ms | 8.26 µs | 2016.27x | 0e+00 |
| `grid/cubic/2d/each10/eval1000/ascending` | construction | interpn | 1.47 µs | 0 | 1672.00x | 3e-01 |
| `grid/nearest/2d/each40/eval10000/ascending` | construction | interpn | 638 ns | 0 | 933.00x | 0e+00 |
| `grid/linear/2d/max40/eval100/ascending` | construction | interpn | 395 ns | 0 | 578.00x | 2e-16 |
| `grid/linear/2d/max40/eval1000/ascending` | construction | interpn | 392 ns | 0 | 573.71x | 3e-16 |
| `grid/linear/2d/max40/eval10/ascending` | construction | interpn | 388 ns | 0 | 567.43x | 2e-16 |
| `grid/linear/2d/max40/eval10000/ascending` | construction | interpn | 362 ns | 0 | 528.86x | 4e-16 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | ndarray-interp | 1.08 ms | 2.22 µs | 484.89x | 6e-16 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | ndarray-interp | 330 µs | 1.12 µs | 293.80x | 5e-16 |
| `grid/linear/3d/max200/eval100/ascending` | construction | SciPy | 1.76 ms | 7.62 µs | 230.54x |  |
| `grid/linear/3d/max200/eval10/ascending` | construction | SciPy | 1.73 ms | 7.63 µs | 227.14x |  |
| `grid/linear/3d/max10/eval10/ascending` | construction | interpn | 196 ns | 0 | 222.89x | 2e-16 |
| `grid/linear/3d/max10/eval1000/ascending` | construction | interpn | 195 ns | 0 | 222.33x | 5e-16 |
| `grid/linear/3d/max10/eval10000/ascending` | construction | interpn | 194 ns | 0 | 220.89x | 4e-16 |
| `grid/linear/3d/max200/eval10000/ascending` | construction | SciPy | 1.61 ms | 7.53 µs | 213.30x |  |
| `grid/linear/3d/max200/eval1000/ascending` | construction | SciPy | 1.61 ms | 7.65 µs | 211.08x |  |
| `grid/linear/3d/max200/eval1000/flipped` | construction | SciPy | 3.25 ms | 16.1 µs | 201.58x |  |
| `grid/linear/3d/max200/eval100/flipped` | construction | SciPy | 3.18 ms | 16 µs | 198.02x |  |
| `grid/linear/2d/max10/eval10000/ascending` | construction | interpn | 134 ns | 0 | 196.14x | 3e-16 |
| `grid/linear/2d/max10/eval1000/ascending` | construction | interpn | 133 ns | 0 | 195.29x | 3e-16 |
| `grid/linear/3d/max200/eval10000/flipped` | construction | SciPy | 3.08 ms | 15.8 µs | 194.69x |  |
| `grid/linear/2d/max10/eval100/ascending` | construction | interpn | 133 ns | 0 | 194.21x | 2e-16 |
| `grid/linear/3d/max200/eval10/flipped` | construction | SciPy | 3.11 ms | 16.1 µs | 192.46x |  |
| `grid/linear/3d/max10/eval100/ascending` | construction | interpn | 179 ns | 0 | 183.10x | 3e-16 |
| `grid/linear/2d/max10/eval10/ascending` | construction | interpn | 124 ns | 0 | 181.57x | 1e-16 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | SciPy | 1.08 ms | 6.29 µs | 171.45x |  |
| `grid/linear/2d/max200/eval10000/ascending` | construction | ndarray-interp | 44.8 µs | 486 ns | 92.03x | 5e-16 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | SciPy | 330 µs | 6.04 µs | 54.69x |  |
| `grid/pchip/3d/each40/eval10000/ascending` | construction | SciPy | 512 µs | 9.48 µs | 54.01x |  |
| `grid/linear/2d/max200/eval1000/ascending` | construction | ndarray-interp | 6.85 µs | 469 ns | 14.61x | 4e-16 |
| `grid/linear/2d/max200/eval100/ascending` | construction | ndarray-interp | 6.82 µs | 469 ns | 14.53x | 3e-16 |
| `grid/linear/2d/max200/eval10/ascending` | construction | ndarray-interp | 6.86 µs | 519 ns | 13.20x | 2e-16 |
| `grid/linear/2d/max200/eval10000/ascending` | construction | SciPy | 44.8 µs | 5.75 µs | 7.78x |  |
| `grid/linear/3d/max40/eval10000/ascending` | construction | scirs2-interpolate | 14.2 µs | 2.58 µs | 5.48x | 2e-16 |
| `oned/linear/n10000/sorted-uneven` | construction | ndarray-interp | 90.1 µs | 17.6 µs | 5.12x | 0e+00 |
| `oned/linear/n10000/sorted` | construction | ndarray-interp | 94 µs | 20.2 µs | 4.66x | 0e+00 |
| `oned/linear/n10000/shuffled-uneven` | construction | ndarray-interp | 89.9 µs | 20.1 µs | 4.46x | 0e+00 |
| `grid/linear/2d/max200/eval10000/flipped` | construction | SciPy | 52.3 µs | 11.9 µs | 4.41x |  |
| `oned/linear/n10000/shuffled` | construction | ndarray-interp | 94.5 µs | 21.6 µs | 4.38x | 0e+00 |
| `oned/pchip/n1000/sorted` | construction | scirs2-interpolate | 17.9 µs | 4.17 µs | 4.29x | 2e-16 |
| `oned/bspline1/n10000/sorted` | construction | SciPy | 67.2 µs | 22.8 µs | 2.95x |  |
| `grid/linear/4d/max40/eval10000/flipped` | construction | SciPy | 57.2 µs | 19.8 µs | 2.89x |  |
| `grid/linear/4d/max40/eval100/flipped` | construction | SciPy | 56.4 µs | 19.7 µs | 2.87x |  |
| `grid/linear/4d/max40/eval10/flipped` | construction | SciPy | 53.6 µs | 19.7 µs | 2.72x |  |
| `oned/bspline1/n10000/shuffled` | construction | SciPy | 61.4 µs | 23.2 µs | 2.65x |  |
| `grid/linear/4d/max40/eval1000/flipped` | construction | SciPy | 53.9 µs | 20.5 µs | 2.63x |  |
| `grid/linear/2d/max40/eval100/ascending` | construction | ndarray-interp | 395 ns | 151 ns | 2.62x | 2e-16 |
| `grid/linear/2d/max40/eval1000/ascending` | construction | ndarray-interp | 392 ns | 151 ns | 2.60x | 4e-16 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | SciPy | 21.1 µs | 8.17 µs | 2.58x |  |
| `grid/linear/2d/max40/eval10/ascending` | construction | ndarray-interp | 388 ns | 151 ns | 2.57x | 2e-16 |
| `oned/pchip/n1000/shuffled` | construction | scirs2-interpolate | 11.6 µs | 4.54 µs | 2.56x | 2e-16 |
| `grid/linear/4d/max40/eval1000/ascending` | construction | scirs2-interpolate | 13.7 µs | 5.55 µs | 2.46x | 2e-16 |
| `grid/linear/4d/max40/eval10/ascending` | construction | scirs2-interpolate | 13.1 µs | 5.45 µs | 2.40x | 2e-16 |
| `grid/linear/4d/max40/eval100/ascending` | construction | scirs2-interpolate | 13.1 µs | 5.46 µs | 2.39x | 3e-16 |
| `grid/linear/3d/max40/eval1000/ascending` | construction | scirs2-interpolate | 14.1 µs | 5.99 µs | 2.35x | 2e-16 |
| `grid/nearest/2d/each40/eval10000/ascending` | construction | scirs2-interpolate | 638 ns | 274 ns | 2.33x | 0e+00 |
| `grid/linear/2d/max200/eval100/ascending` | construction | scirs2-interpolate | 6.82 µs | 2.95 µs | 2.31x | 2e-16 |
| `grid/linear/2d/max200/eval10/ascending` | construction | scirs2-interpolate | 6.86 µs | 2.99 µs | 2.29x | 1e-16 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | scirs2-interpolate | 6.85 µs | 3.02 µs | 2.27x | 2e-16 |
| `oned/cubic_spline/n1000/shuffled` | construction | scirs2-interpolate | 37.5 µs | 16.6 µs | 2.26x | 4e-04 |
| `grid/nearest/3d/each40/eval10000/ascending` | construction | scirs2-interpolate | 21.1 µs | 9.37 µs | 2.25x | 0e+00 |
| `oned/cubic_spline/n1000/sorted` | construction | scirs2-interpolate | 36.4 µs | 16.5 µs | 2.21x | 4e-04 |
| `grid/linear/2d/max40/eval10000/ascending` | construction | ndarray-interp | 362 ns | 165 ns | 2.19x | 4e-16 |
| `grid/linear/4d/max40/eval10000/ascending` | construction | scirs2-interpolate | 13.1 µs | 6.04 µs | 2.17x | 4e-16 |
| `oned/linear/n10/sorted` | construction | ndarray-interp | 93 ns | 43.5 ns | 2.14x | 0e+00 |
| `oned/linear/n10/shuffled` | construction | ndarray-interp | 92.4 ns | 43.4 ns | 2.13x | 0e+00 |
| `oned/linear/n10000/shuffled` | construction | SciPy | 94.5 µs | 47.8 µs | 1.98x |  |
| `grid/linear/2d/max40/eval10/ascending` | construction | scirs2-interpolate | 388 ns | 199 ns | 1.95x | 4e-17 |
| `oned/linear/n10000/shuffled-uneven` | construction | SciPy | 89.9 µs | 46.4 µs | 1.94x |  |
| `grid/linear/2d/max40/eval1000/ascending` | construction | scirs2-interpolate | 392 ns | 204 ns | 1.93x | 3e-16 |
| `oned/linear/n10000/sorted-uneven` | construction | SciPy | 90.1 µs | 47.1 µs | 1.91x |  |
| `oned/linear/n10000/sorted` | construction | SciPy | 94 µs | 50.7 µs | 1.85x |  |
| `oned/akima/n1000/sorted` | construction | scirs2-interpolate | 10 µs | 5.55 µs | 1.81x | 3e-17 |
| `grid/linear/2d/max40/eval100/ascending` | construction | scirs2-interpolate | 395 ns | 219 ns | 1.81x | 3e-16 |
| `grid/linear/3d/max40/eval10000/flipped` | construction | SciPy | 27.9 µs | 15.6 µs | 1.78x |  |
| `grid/linear/3d/max40/eval1000/flipped` | construction | SciPy | 27.9 µs | 15.7 µs | 1.78x |  |
| `grid/linear/2d/max10/eval1000/ascending` | construction | ndarray-interp | 133 ns | 77.7 ns | 1.72x | 3e-16 |
| `grid/linear/2d/max10/eval100/ascending` | construction | ndarray-interp | 133 ns | 78.6 ns | 1.69x | 4e-16 |
| `oned/linear/n50/shuffled` | construction | ndarray-interp | 169 ns | 101 ns | 1.67x | 0e+00 |
| `oned/akima/n1000/shuffled` | construction | scirs2-interpolate | 9.47 µs | 5.7 µs | 1.66x | 3e-17 |
| `oned/pchip/n100/shuffled-uneven` | construction | scirs2-interpolate | 892 ns | 537 ns | 1.66x | 2e-16 |
| `oned/pchip/n100/shuffled` | construction | scirs2-interpolate | 973 ns | 588 ns | 1.66x | 2e-16 |
| `grid/linear/2d/max40/eval10000/ascending` | construction | scirs2-interpolate | 362 ns | 219 ns | 1.65x | 2e-16 |
| `oned/linear/n10000/sorted` | construction | scirs2-interpolate | 94 µs | 57.2 µs | 1.64x | 0e+00 |
| `oned/pchip/n100/sorted` | construction | scirs2-interpolate | 898 ns | 549 ns | 1.64x | 2e-16 |
| `oned/linear/n1000/sorted` | construction | scirs2-interpolate | 1.99 µs | 1.21 µs | 1.64x | 0e+00 |
| `oned/pchip/n100/sorted-uneven` | construction | scirs2-interpolate | 891 ns | 545 ns | 1.64x | 2e-16 |
| `oned/linear/n10000/shuffled-uneven` | construction | scirs2-interpolate | 89.9 µs | 55.7 µs | 1.61x | 6e-17 |
| `oned/linear/n10000/sorted-uneven` | construction | scirs2-interpolate | 90.1 µs | 56.2 µs | 1.60x | 6e-17 |
| `grid/linear/2d/max10/eval10/ascending` | construction | ndarray-interp | 124 ns | 77.6 ns | 1.60x | 1e-16 |
| `oned/linear/n1000/shuffled` | construction | scirs2-interpolate | 1.95 µs | 1.23 µs | 1.58x | 0e+00 |
| `grid/linear/2d/max10/eval10000/ascending` | construction | ndarray-interp | 134 ns | 85.3 ns | 1.57x | 4e-16 |
| `oned/linear/n10000/shuffled` | construction | scirs2-interpolate | 94.5 µs | 61.6 µs | 1.53x | 0e+00 |
| `oned/pchip/n50/sorted` | construction | scirs2-interpolate | 506 ns | 334 ns | 1.51x | 1e-16 |
| `grid/linear/3d/max200/eval100/ascending` | construction | scirs2-interpolate | 1.76 ms | 1.18 ms | 1.49x | 3e-16 |
| `oned/linear/n100/sorted-uneven` | construction | ndarray-interp | 249 ns | 169 ns | 1.48x | 0e+00 |
| `oned/linear/n100/shuffled-uneven` | construction | ndarray-interp | 248 ns | 169 ns | 1.47x | 0e+00 |
| `oned/linear/n100/sorted` | construction | ndarray-interp | 247 ns | 169 ns | 1.46x | 0e+00 |
| `oned/linear/n100/shuffled` | construction | ndarray-interp | 247 ns | 169 ns | 1.46x | 0e+00 |
| `grid/linear/4d/max40/eval1000/ascending` | construction | SciPy | 13.7 µs | 9.39 µs | 1.46x |  |
| `oned/akima/n50/sorted` | construction | scirs2-interpolate | 468 ns | 323 ns | 1.45x | 1e-16 |
| `oned/linear/n1000/sorted` | construction | ndarray-interp | 1.99 µs | 1.38 µs | 1.44x | 0e+00 |
| `oned/linear/n50/sorted` | construction | ndarray-interp | 160 ns | 111 ns | 1.44x | 0e+00 |
| `grid/linear/3d/max10/eval10/ascending` | evaluation | interpn | 201 ns | 142 ns | 1.42x | 2e-16 |
| `oned/pchip/n50/shuffled` | construction | scirs2-interpolate | 494 ns | 348 ns | 1.42x | 1e-16 |
| `oned/linear/n50/shuffled` | construction | scirs2-interpolate | 169 ns | 119 ns | 1.42x | 0e+00 |
| `oned/linear/n1000/shuffled` | construction | ndarray-interp | 1.95 µs | 1.38 µs | 1.41x | 0e+00 |
| `oned/akima/n100/shuffled` | construction | scirs2-interpolate | 741 ns | 531 ns | 1.40x | 1e-16 |
| `grid/linear/2d/max10/eval10/ascending` | evaluation | interpn | 132 ns | 95 ns | 1.39x | 1e-16 |
| `grid/linear/4d/max40/eval100/ascending` | construction | SciPy | 13.1 µs | 9.46 µs | 1.38x |  |
| `oned/akima/n100/sorted` | construction | scirs2-interpolate | 741 ns | 538 ns | 1.38x | 1e-16 |
| `grid/linear/4d/max40/eval10/ascending` | construction | SciPy | 13.1 µs | 9.5 µs | 1.38x |  |
| `grid/linear/2d/max10/eval1000/ascending` | evaluation | interpn | 12.5 µs | 9.23 µs | 1.35x | 3e-16 |
| `grid/linear/3d/max200/eval10000/ascending` | construction | scirs2-interpolate | 1.61 ms | 1.19 ms | 1.35x | 3e-16 |
| `grid/linear/3d/max200/eval10/ascending` | construction | scirs2-interpolate | 1.73 ms | 1.29 ms | 1.34x | 3e-16 |
| `grid/linear/2d/max10/eval100/ascending` | evaluation | interpn | 1.25 µs | 927 ns | 1.34x | 2e-16 |
| `grid/linear/3d/max10/eval1000/ascending` | evaluation | interpn | 18.9 µs | 14.1 µs | 1.34x | 5e-16 |
| `grid/linear/3d/max10/eval10000/ascending` | evaluation | interpn | 188 µs | 141 µs | 1.33x | 4e-16 |
| `oned/akima/n50/shuffled` | construction | scirs2-interpolate | 427 ns | 326 ns | 1.31x | 1e-16 |
| `grid/linear/2d/each800/eval10000/ascending` | construction | scirs2-interpolate | 1.08 ms | 825 µs | 1.31x | 2e-16 |
| `grid/linear/2d/max200/eval100/flipped` | construction | SciPy | 15.9 µs | 12.2 µs | 1.31x |  |
| `grid/linear/4d/max200/eval10000/ascending` | construction | scirs2-interpolate | 48.1 ms | 37 ms | 1.30x | 4e-16 |
| `grid/linear/4d/max40/eval10000/ascending` | construction | SciPy | 13.1 µs | 10.2 µs | 1.28x |  |
| `oned/linear/n50/sorted` | construction | scirs2-interpolate | 160 ns | 126 ns | 1.27x | 0e+00 |
| `oned/linear/n100/sorted` | construction | scirs2-interpolate | 247 ns | 195 ns | 1.27x | 0e+00 |
| `oned/linear/n100/shuffled-uneven` | construction | scirs2-interpolate | 248 ns | 196 ns | 1.26x | 3e-17 |
| `oned/nearest/n10000/sorted` | construction | SciPy | 67.2 µs | 53.3 µs | 1.26x |  |
| `grid/linear/4d/max200/eval100/ascending` | construction | scirs2-interpolate | 48.3 ms | 38.4 ms | 1.26x | 1e-16 |
| `oned/linear/n100/sorted-uneven` | construction | scirs2-interpolate | 249 ns | 198 ns | 1.26x | 3e-17 |
| `oned/linear/n100/shuffled` | construction | scirs2-interpolate | 247 ns | 196 ns | 1.26x | 0e+00 |
| `grid/linear/2d/max40/eval100/ascending` | evaluation | interpn | 1.66 µs | 1.32 µs | 1.26x | 2e-16 |
| `grid/linear/2d/auckland-ba32-0401/heldout/ascending` | construction | scirs2-interpolate | 16.6 ms | 13.3 ms | 1.25x | 0e+00 |
| `grid/linear/3d/max10/eval100/ascending` | evaluation | interpn | 1.91 µs | 1.53 µs | 1.25x | 3e-16 |
| `grid/nearest/2d/auckland-ba32-0401/heldout/ascending` | construction | scirs2-interpolate | 16.3 ms | 13.1 ms | 1.24x | 0e+00 |
| `grid/linear/2d/max40/eval1000/ascending` | evaluation | interpn | 16.5 µs | 13.3 µs | 1.24x | 3e-16 |
| `grid/linear/3d/max200/eval1000/ascending` | construction | scirs2-interpolate | 1.61 ms | 1.31 ms | 1.23x | 3e-16 |
| `grid/linear/2d/max10/eval10000/ascending` | evaluation | interpn | 125 µs | 101 µs | 1.23x | 3e-16 |
| `grid/linear/2d/max200/eval1000/flipped` | construction | SciPy | 14.5 µs | 11.9 µs | 1.21x |  |
| `grid/linear/4d/max200/eval1000/ascending` | construction | scirs2-interpolate | 45.1 ms | 37.2 ms | 1.21x | 3e-16 |
| `grid/linear/4d/max200/eval10/ascending` | construction | scirs2-interpolate | 44.4 ms | 36.8 ms | 1.21x | 2e-16 |
| `grid/linear/2d/each400/eval10000/ascending` | construction | scirs2-interpolate | 330 µs | 275 µs | 1.20x | 2e-16 |
| `grid/linear/4d/max40/eval1000/ascending` | evaluation | interpn | 32.2 µs | 27 µs | 1.19x | 6e-16 |
| `grid/linear/2d/max200/eval1000/ascending` | construction | SciPy | 6.85 µs | 5.77 µs | 1.19x |  |
| `grid/linear/2d/max40/eval10/ascending` | evaluation | interpn | 159 ns | 134 ns | 1.19x | 2e-16 |
| `grid/linear/3d/max10/eval10000/ascending` | construction | scirs2-interpolate | 194 ns | 164 ns | 1.19x | 3e-16 |
| `oned/nearest/n10000/shuffled` | construction | SciPy | 61.9 µs | 52.2 µs | 1.19x |  |
| `grid/linear/2d/max200/eval10/ascending` | construction | SciPy | 6.86 µs | 5.84 µs | 1.17x |  |
| `grid/linear/3d/max10/eval10/ascending` | construction | scirs2-interpolate | 196 ns | 169 ns | 1.16x | 2e-16 |
| `grid/linear/3d/max10/eval1000/ascending` | construction | scirs2-interpolate | 195 ns | 169 ns | 1.15x | 3e-16 |
| `grid/linear/4d/max40/eval10/ascending` | evaluation | interpn | 309 ns | 270 ns | 1.15x | 2e-16 |
| `grid/linear/2d/max40/eval100/ascending` | evaluation | ndarray-interp | 1.66 µs | 1.45 µs | 1.14x | 2e-16 |
| `grid/linear/2d/max200/eval10/flipped` | construction | SciPy | 14.5 µs | 12.7 µs | 1.14x |  |
| `grid/linear/2d/max200/eval100/ascending` | construction | SciPy | 6.82 µs | 6 µs | 1.14x |  |
| `oned/linear/n10/sorted` | construction | scirs2-interpolate | 93 ns | 81.9 ns | 1.13x | 0e+00 |
| `grid/linear/2d/max40/eval10000/ascending` | evaluation | interpn | 152 µs | 134 µs | 1.13x | 4e-16 |
| `grid/linear/3d/max200/eval1000/ascending` | evaluation | interpn | 29.6 µs | 26.1 µs | 1.13x | 3e-16 |
| `oned/cubic_spline/n100/sorted-uneven` | construction | scirs2-interpolate | 1.72 µs | 1.53 µs | 1.13x | 6e-03 |
| `grid/linear/2d/max200/eval10000/ascending` | evaluation | interpn | 193 µs | 172 µs | 1.12x | 3e-16 |
| `grid/nearest/2d/each40/eval10000/ascending` | evaluation | interpn | 151 µs | 135 µs | 1.12x | 0e+00 |
| `oned/nearest/n10000/sorted` | construction | scirs2-interpolate | 67.2 µs | 60.3 µs | 1.11x | 0e+00 |
| `oned/linear/n10/shuffled` | construction | scirs2-interpolate | 92.4 ns | 83.3 ns | 1.11x | 0e+00 |
| `grid/linear/3d/max40/eval10/ascending` | evaluation | interpn | 220 ns | 199 ns | 1.10x | 2e-16 |
| `grid/linear/4d/max40/eval10000/ascending` | evaluation | interpn | 298 µs | 270 µs | 1.10x | 5e-16 |

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
| `grid/pchip/2d/auckland-ba32-0401/heldout/ascending` | SciPy | MemoryError: Unable to allocate 695. GiB for an array with shape (25913999, 3601) and data type float64 |

## Real-world data

Two data sets, each taken whole in the unit its source gives it out in; what was taken and how it was divided was fixed before anything was measured (`compare/perf/real.py`, docs/design.ja.md section 8.4.7).

- `auckland-ba32-0401`: Auckland Part 1 LiDAR 1m DEM (2024), tile BA32_10000_0401. 4800 by 7200 cells of 1 m. The grid is every second row and column (and the last ones); every other cell is a query, and its elevation is the value held out. Licensed by Regional Software Holdings Limited for re-use under CC-BY-4.0. Hosted by Toitū Te Whenua Land Information New Zealand (https://registry.opendata.aws/nz-elevation).
- `ghcn-tmax-20250101`: GHCN-Daily, maximum temperature of 2025-01-01, every station; longitude and latitude in degrees, as they are. Of every 10 stations in the order of the file the last is a query, and its temperature is the value held out. The RBF kernels are those without a shape parameter. NOAA Global Historical Climatology Network Daily (GHCN-D), from https://registry.opendata.aws/noaa-ghcn (CC0 1.0).

These cases are also in the tables of their families below and in the counts above. Here they are with the error: the largest and the root mean square distance from the values held out (metres, and degrees Celsius), and the number of queries a library gave no value. A library missing from a case does not have it, or did not finish (see the runs that failed).

`grid/nearest/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 16.3 ms | 346 ms | 15.5 | 0.178 | 0 |
| interpx, parallel | 16.4 ms | 61 ms | 15.5 | 0.178 | 0 |
| interpx, `from_owned` | 5.18 ms | 347 ms | 15.5 | 0.178 | 0 |
| SciPy | 6.96 µs | 1.32 s | 15.5 | 0.178 | 0 |
| scirs2-interpolate | 13.1 ms | 150 s | 15.5 | 0.178 | 0 |
| interpn | 0 | 716 ms | 15.5 | 0.178 | 0 |

`grid/linear/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 16.6 ms | 370 ms | 6.74 | 0.0601 | 0 |
| interpx, parallel | 15.9 ms | 72.6 ms | 6.74 | 0.0601 | 0 |
| interpx, `from_owned` | 5.81 ms | 395 ms | 6.74 | 0.0601 | 0 |
| SciPy | 6.92 µs | 1.07 s | 6.74 | 0.0601 | 0 |
| scirs2-interpolate | 13.3 ms | 75.8 s | 6.74 | 0.0601 | 0 |
| ndarray-interp | 8.26 µs | 373 ms | 6.74 | 0.0601 | 0 |
| interpn | 0 | 841 ms | 6.74 | 0.0601 | 0 |

`grid/cubic/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 76.7 ms | 725 ms | 6.77 | 0.0507 | 0 |
| interpx, parallel | 78.2 ms | 140 ms | 6.77 | 0.0507 | 0 |
| interpx, `from_owned` | 64.9 ms | 748 ms | 6.77 | 0.0507 | 0 |
| SciPy | 14.6 s | 33 s | 6.77 | 0.0507 | 0 |
| interpn | 0 | 1.64 s | 6.76 | 0.0517 | 0 |

`grid/quintic/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 101 ms | 1.3 s | 6.77 | 0.0517 | 0 |
| interpx, parallel | 91.3 ms | 209 ms | 6.77 | 0.0517 | 0 |
| interpx, `from_owned` | 83.4 ms | 1.29 s | 6.77 | 0.0517 | 0 |
| SciPy | 21.3 s | 37.8 s | 6.77 | 0.0517 | 0 |

`grid/pchip/2d/auckland-ba32-0401/heldout/ascending` (8646001 points, 25913999 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 79.2 ms | 1.13 s | 6.75 | 0.049 | 0 |
| interpx, parallel | 71.8 ms | 176 ms | 6.75 | 0.049 | 0 |
| interpx, `from_owned` | 66.4 ms | 1.03 s | 6.75 | 0.049 | 0 |

`nearest_nd/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 675 µs | 159 µs | 25 | 3.42 | 0 |
| interpx, parallel | 721 µs | 70.5 µs | 25 | 3.42 | 0 |
| interpx, `from_owned` | 655 µs | 160 µs | 25 | 3.42 | 0 |
| SciPy | 1.83 ms | 524 µs | 25 | 3.42 | 0 |

`rbf/linear/global/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 18.3 s | 19 ms | 20.8 | 2.74 | 0 |
| interpx, parallel | 19.5 s | 4.24 ms | 20.8 | 2.74 | 0 |
| interpx, `from_owned` | 18.4 s | 18.1 ms | 20.8 | 2.74 | 0 |
| SciPy | 114 s | 55.7 ms | 20.8 | 2.74 | 0 |
| scirs2-interpolate | 314 s | 258 ms | 20.8 | 2.74 | 0 |

`rbf/linear/k50/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 767 µs | 18.2 ms | 20.8 | 2.75 | 0 |
| interpx, parallel | 811 µs | 3.24 ms | 20.8 | 2.75 | 0 |
| interpx, `from_owned` | 762 µs | 19.2 ms | 20.8 | 2.75 | 0 |
| SciPy | 1.77 ms | 44.5 ms | 20.8 | 2.75 | 0 |

`rbf/thin_plate_spline/global/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 18.8 s | 87.8 ms | 102 | 5.67 | 0 |
| interpx, parallel | 37.3 s | 68.6 ms | 102 | 5.67 | 0 |
| interpx, `from_owned` | 18.2 s | 88.1 ms | 102 | 5.67 | 0 |
| SciPy | 116 s | 124 ms | 102 | 5.67 | 0 |
| scirs2-interpolate | 273 s | 412 ms | 93.7 | 5.41 | 0 |

`rbf/thin_plate_spline/k50/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 768 µs | 31 ms | 102 | 5.66 | 0 |
| interpx, parallel | 790 µs | 5.19 ms | 102 | 5.66 | 0 |
| interpx, `from_owned` | 702 µs | 28.6 ms | 102 | 5.66 | 0 |
| SciPy | 1.87 ms | 59.9 ms | 102 | 5.66 | 0 |

`rbf/cubic/global/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 17.2 s | 17.9 ms | 587 | 27.1 | 0 |
| interpx, parallel | 15.9 s | 3.46 ms | 587 | 27.1 | 0 |
| interpx, `from_owned` | 17.8 s | 18.7 ms | 587 | 27.1 | 0 |
| SciPy | 112 s | 53.1 ms | 597 | 24.4 | 0 |
| scirs2-interpolate | 264 s | 270 ms | 100 | 8.48 | 0 |

`rbf/cubic/k50/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 762 µs | 19.4 ms | 567 | 25.3 | 0 |
| interpx, parallel | 767 µs | 3.93 ms | 567 | 25.3 | 0 |
| interpx, `from_owned` | 695 µs | 21.6 ms | 567 | 25.3 | 0 |
| SciPy | 1.68 ms | 44.4 ms | 567 | 25.3 | 0 |

`rbf/quintic/global/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 17.4 s | 21.9 ms | 2.55e+05 | 1.36e+04 | 0 |
| interpx, parallel | 15.8 s | 4.1 ms | 2.55e+05 | 1.36e+04 | 0 |
| interpx, `from_owned` | 16.1 s | 20.5 ms | 2.55e+05 | 1.36e+04 | 0 |
| SciPy | 116 s | 54.4 ms | 9.39e+04 | 4.92e+03 | 0 |
| scirs2-interpolate | 272 s | 279 ms | 6.56e+16 | 2.85e+15 | 0 |

`rbf/quintic/k50/ghcn-tmax-20250101/heldout` (10838 points, 1204 queries)

| library | construction | evaluation | largest error | RMS error | no value |
| :--- | ---: | ---: | ---: | ---: | ---: |
| interpx | 796 µs | 21.2 ms | 3.74e+03 | 199 | 0 |
| interpx, parallel | 773 µs | 4.11 ms | 3.74e+03 | 199 | 0 |
| interpx, `from_owned` | 697 µs | 21.2 ms | 3.74e+03 | 199 | 0 |
| SciPy | 1.7 ms | 49.3 ms | 3.74e+03 | 200 | 0 |

## One variable

Data of `n` points, evaluated at `4 n` points: in order (`sorted`) and in no order (`shuffled`).

### Evaluation

| case | interpx | interpx, parallel | SciPy | scirs2-interpolate | ndarray-interp | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/n10/sorted` | 112 ns | 173 ns | 4.9 µs | 285 ns | 253 ns | 2.3x faster |
| `linear/n10/shuffled` | 170 ns | 224 ns | 4.99 µs | 291 ns | 250 ns | 1.5x faster |
| `nearest/n10/sorted` | 148 ns | 192 ns | 9.35 µs | 268 ns |  | 1.8x faster |
| `nearest/n10/shuffled` | 296 ns | 358 ns | 9.37 µs | 286 ns |  | **1x slower** |
| `previous/n10/sorted` | 151 ns | 192 ns | 10.6 µs |  |  | 70x faster |
| `previous/n10/shuffled` | 267 ns | 319 ns | 10.8 µs |  |  | 40x faster |
| `bspline1/n10/sorted` | 220 ns | 290 ns | 2.67 µs |  |  | 12x faster |
| `bspline1/n10/shuffled` | 527 ns | 595 ns | 2.73 µs |  |  | 5.2x faster |
| `bspline2/n10/sorted` | 272 ns | 335 ns | 3.08 µs |  |  | 11x faster |
| `bspline2/n10/shuffled` | 373 ns | 413 ns | 3.09 µs |  |  | 8.3x faster |
| `bspline3/n10/sorted` | 339 ns | 433 ns | 3.61 µs |  |  | 11x faster |
| `bspline3/n10/shuffled` | 440 ns | 511 ns | 3.57 µs |  |  | 8.1x faster |
| `cubic_spline/n10/sorted` | 141 ns | 193 ns | 3.67 µs | 269 ns | 327 ns | 1.9x faster |
| `cubic_spline/n10/shuffled` | 211 ns | 264 ns | 3.61 µs | 277 ns | 326 ns | 1.3x faster |
| `pchip/n10/sorted` | 142 ns | 211 ns | 3.53 µs | 321 ns |  | 2.3x faster |
| `pchip/n10/shuffled` | 205 ns | 291 ns | 3.58 µs | 317 ns |  | 1.5x faster |
| `akima/n10/sorted` | 151 ns | 208 ns | 3.98 µs | 287 ns |  | 1.9x faster |
| `akima/n10/shuffled` | 205 ns | 279 ns | 3.68 µs | 261 ns |  | 1.3x faster |
| `makima/n10/sorted` | 143 ns | 193 ns | 3.68 µs |  |  | 26x faster |
| `makima/n10/shuffled` | 224 ns | 251 ns | 3.65 µs |  |  | 16x faster |
| `linear/n50/sorted` | 641 ns | 688 ns | 6.64 µs | 1.9 µs | 1.35 µs | 2.1x faster |
| `linear/n50/shuffled` | 811 ns | 861 ns | 8.47 µs | 1.79 µs | 1.24 µs | 1.5x faster |
| `nearest/n50/sorted` | 820 ns | 765 ns | 11.5 µs | 1.82 µs |  | 2.2x faster |
| `nearest/n50/shuffled` | 1.57 µs | 1.59 µs | 11.6 µs | 1.9 µs |  | 1.2x faster |
| `previous/n50/sorted` | 733 ns | 765 ns | 12.8 µs |  |  | 17x faster |
| `previous/n50/shuffled` | 1.37 µs | 1.42 µs | 12.6 µs |  |  | 9.2x faster |
| `bspline1/n50/sorted` | 1.1 µs | 1.36 µs | 3.8 µs |  |  | 3.4x faster |
| `bspline1/n50/shuffled` | 2.77 µs | 2.85 µs | 4.87 µs |  |  | 1.8x faster |
| `bspline2/n50/sorted` | 1.36 µs | 1.59 µs | 6.03 µs |  |  | 4.4x faster |
| `bspline2/n50/shuffled` | 2.02 µs | 2.1 µs | 6.65 µs |  |  | 3.3x faster |
| `bspline3/n50/sorted` | 1.83 µs | 1.89 µs | 8.67 µs |  |  | 4.7x faster |
| `bspline3/n50/shuffled` | 2.81 µs | 2.95 µs | 9.35 µs |  |  | 3.3x faster |
| `cubic_spline/n50/sorted` | 699 ns | 807 ns | 5.11 µs | 4.49 µs | 1.59 µs | 2.3x faster |
| `cubic_spline/n50/shuffled` | 986 ns | 1.03 µs | 5.48 µs | 4.14 µs | 1.6 µs | 1.6x faster |
| `pchip/n50/sorted` | 696 ns | 804 ns | 5.18 µs | 4.29 µs |  | 6.2x faster |
| `pchip/n50/shuffled` | 1.07 µs | 1.03 µs | 5.41 µs | 3.64 µs |  | 3.4x faster |
| `akima/n50/sorted` | 712 ns | 738 ns | 5.3 µs | 4.12 µs |  | 5.8x faster |
| `akima/n50/shuffled` | 980 ns | 1.09 µs | 5.45 µs | 4.14 µs |  | 4.2x faster |
| `makima/n50/sorted` | 711 ns | 749 ns | 5.09 µs |  |  | 7.2x faster |
| `makima/n50/shuffled` | 965 ns | 1.03 µs | 5.37 µs |  |  | 5.6x faster |
| `linear/n100/sorted` | 1.31 µs | 1.49 µs | 7.18 µs | 3.66 µs | 2.46 µs | 1.9x faster |
| `linear/n100/shuffled` | 1.62 µs | 1.71 µs | 12.4 µs | 3.82 µs | 2.49 µs | 1.5x faster |
| `nearest/n100/sorted` | 1.49 µs | 1.58 µs | 13.1 µs | 3.57 µs |  | 2.4x faster |
| `nearest/n100/shuffled` | 3.07 µs | 3.2 µs | 13.3 µs | 3.67 µs |  | 1.2x faster |
| `previous/n100/sorted` | 1.5 µs | 1.71 µs | 14.6 µs |  |  | 9.8x faster |
| `previous/n100/shuffled` | 2.71 µs | 2.83 µs | 15 µs |  |  | 5.5x faster |
| `bspline1/n100/sorted` | 2.16 µs | 2.99 µs | 5.21 µs |  |  | 2.4x faster |
| `bspline1/n100/shuffled` | 5.58 µs | 5.76 µs | 11.6 µs |  |  | 2.1x faster |
| `bspline2/n100/sorted` | 2.94 µs | 3.29 µs | 10.2 µs |  |  | 3.5x faster |
| `bspline2/n100/shuffled` | 4.87 µs | 4.58 µs | 12.6 µs |  |  | 2.6x faster |
| `bspline3/n100/sorted` | 3.64 µs | 3.78 µs | 14.9 µs |  |  | 4.1x faster |
| `bspline3/n100/shuffled` | 6.15 µs | 6.27 µs | 17.8 µs |  |  | 2.9x faster |
| `cubic_spline/n100/sorted` | 1.42 µs | 1.66 µs | 6.92 µs | 19 µs | 3.17 µs | 2.2x faster |
| `cubic_spline/n100/shuffled` | 1.94 µs | 2.06 µs | 7.92 µs | 18.2 µs | 3.2 µs | 1.6x faster |
| `pchip/n100/sorted` | 1.44 µs | 1.53 µs | 6.79 µs | 16 µs |  | 4.7x faster |
| `pchip/n100/shuffled` | 1.94 µs | 2.06 µs | 7.76 µs | 17.7 µs |  | 4x faster |
| `akima/n100/sorted` | 1.41 µs | 1.56 µs | 6.8 µs | 18.1 µs |  | 4.8x faster |
| `akima/n100/shuffled` | 1.95 µs | 2.05 µs | 7.73 µs | 17.3 µs |  | 4x faster |
| `makima/n100/sorted` | 1.41 µs | 1.54 µs | 6.86 µs |  |  | 4.9x faster |
| `makima/n100/shuffled` | 1.94 µs | 2.04 µs | 8.05 µs |  |  | 4.1x faster |
| `linear/n1000/sorted` | 13.1 µs | 13.8 µs | 28.6 µs | 109 µs | 24.5 µs | 1.9x faster |
| `linear/n1000/shuffled` | 16.1 µs | 16.3 µs | 143 µs | 215 µs | 24.6 µs | 1.5x faster |
| `nearest/n1000/sorted` | 14.1 µs | 14.9 µs | 55.4 µs | 147 µs |  | 3.9x faster |
| `nearest/n1000/shuffled` | 30.9 µs | 31.1 µs | 198 µs | 231 µs |  | 6.4x faster |
| `previous/n1000/sorted` | 14.1 µs | 14.7 µs | 59.6 µs |  |  | 4.2x faster |
| `previous/n1000/shuffled` | 27.3 µs | 27.4 µs | 199 µs |  |  | 7.3x faster |
| `bspline1/n1000/sorted` | 25.4 µs | 26 µs | 30.4 µs |  |  | 1.2x faster |
| `bspline1/n1000/shuffled` | 56.3 µs | 56.5 µs | 527 µs |  |  | 9.4x faster |
| `bspline2/n1000/sorted` | 32.4 µs | 32.5 µs | 78.6 µs |  |  | 2.4x faster |
| `bspline2/n1000/shuffled` | 58.1 µs | 58.4 µs | 538 µs |  |  | 9.3x faster |
| `bspline3/n1000/sorted` | 36.3 µs | 36.6 µs | 128 µs |  |  | 3.5x faster |
| `bspline3/n1000/shuffled` | 78.9 µs | 79 µs | 575 µs |  |  | 7.3x faster |
| `cubic_spline/n1000/sorted` | 13.8 µs | 15.5 µs | 37.3 µs | 1.42 ms | 31.6 µs | 2.3x faster |
| `cubic_spline/n1000/shuffled` | 19.5 µs | 19.8 µs | 219 µs | 1.43 ms | 31.9 µs | 1.6x faster |
| `pchip/n1000/sorted` | 13.7 µs | 14.4 µs | 38 µs | 1.2 ms |  | 2.8x faster |
| `pchip/n1000/shuffled` | 19.6 µs | 19.9 µs | 217 µs | 1.21 ms |  | 11x faster |
| `akima/n1000/sorted` | 13.6 µs | 14.7 µs | 36.5 µs | 1.41 ms |  | 2.7x faster |
| `akima/n1000/shuffled` | 19.5 µs | 19.9 µs | 215 µs | 1.41 ms |  | 11x faster |
| `makima/n1000/sorted` | 13.9 µs | 14.9 µs | 37.2 µs |  |  | 2.7x faster |
| `makima/n1000/shuffled` | 19.6 µs | 19.9 µs | 216 µs |  |  | 11x faster |
| `linear/n10000/sorted` | 138 µs | 74.4 µs | 237 µs | 1.57 ms | 247 µs | 1.7x faster |
| `linear/n10000/shuffled` | 190 µs | 86.5 µs | 2.31 ms | 3.01 ms | 284 µs | 1.5x faster |
| `nearest/n10000/sorted` | 149 µs | 71.5 µs | 584 µs | 1.97 ms |  | 3.9x faster |
| `nearest/n10000/shuffled` | 382 µs | 92.8 µs | 2.91 ms | 3.12 ms |  | 7.6x faster |
| `previous/n10000/sorted` | 150 µs | 71.7 µs | 607 µs |  |  | 4.1x faster |
| `previous/n10000/shuffled` | 344 µs | 95.2 µs | 2.93 ms |  |  | 8.5x faster |
| `bspline1/n10000/sorted` | 267 µs | 99.7 µs | 286 µs |  |  | 1.1x faster |
| `bspline1/n10000/shuffled` | 639 µs | 141 µs | 48.6 ms |  |  | 76x faster |
| `bspline2/n10000/sorted` | 333 µs | 86.2 µs | 759 µs |  |  | 2.3x faster |
| `bspline2/n10000/shuffled` | 902 µs | 197 µs | 49.1 ms |  |  | 54x faster |
| `bspline3/n10000/sorted` | 369 µs | 104 µs | 1.25 ms |  |  | 3.4x faster |
| `bspline3/n10000/shuffled` | 1.17 ms | 244 µs | 47.9 ms |  |  | 41x faster |
| `cubic_spline/n10000/sorted` | 143 µs | 69.6 µs | 364 µs | 139 ms | 319 µs | 2.2x faster |
| `cubic_spline/n10000/shuffled` | 211 µs | 94.4 µs | 3.31 ms | 139 ms | 332 µs | 1.6x faster |
| `pchip/n10000/sorted` | 147 µs | 66.3 µs | 370 µs | 116 ms |  | 2.5x faster |
| `pchip/n10000/shuffled` | 215 µs | 94.7 µs | 3.33 ms | 116 ms |  | 15x faster |
| `akima/n10000/sorted` | 147 µs | 70.4 µs | 362 µs | 142 ms |  | 2.5x faster |
| `akima/n10000/shuffled` | 215 µs | 98 µs | 3.33 ms | 142 ms |  | 15x faster |
| `makima/n10000/sorted` | 145 µs | 70.6 µs | 366 µs |  |  | 2.5x faster |
| `makima/n10000/shuffled` | 214 µs | 95.7 µs | 3.33 ms |  |  | 16x faster |
| `linear/n100/sorted-uneven` | 1.18 µs | 1.27 µs | 7.34 µs | 3.64 µs | 3.77 µs | 3.1x faster |
| `linear/n100/shuffled-uneven` | 2.75 µs | 2.91 µs | 12.6 µs | 3.84 µs | 3.68 µs | 1.3x faster |
| `cubic_spline/n100/sorted-uneven` | 1.34 µs | 1.47 µs | 7.38 µs | 17 µs | 4.41 µs | 3.3x faster |
| `cubic_spline/n100/shuffled-uneven` | 3.06 µs | 3.14 µs | 7.9 µs | 18 µs | 4.43 µs | 1.4x faster |
| `pchip/n100/sorted-uneven` | 1.41 µs | 1.46 µs | 7.3 µs | 16.1 µs |  | 5.2x faster |
| `pchip/n100/shuffled-uneven` | 3.04 µs | 3.14 µs | 7.84 µs | 16.2 µs |  | 2.6x faster |
| `linear/n10000/sorted-uneven` | 206 µs | 92.8 µs | 260 µs | 1.6 ms | 894 µs | 1.3x faster |
| `linear/n10000/shuffled-uneven` | 617 µs | 155 µs | 2.31 ms | 2.9 ms | 1.75 ms | 2.8x faster |
| `cubic_spline/n10000/sorted-uneven` | 218 µs | 93.9 µs | 437 µs | 141 ms | 994 µs | 2x faster |
| `cubic_spline/n10000/shuffled-uneven` | 794 µs | 168 µs | 3.44 ms | 140 ms | 1.83 ms | 2.3x faster |
| `pchip/n10000/sorted-uneven` | 238 µs | 95.8 µs | 437 µs | 117 ms |  | 1.8x faster |
| `pchip/n10000/shuffled-uneven` | 730 µs | 173 µs | 3.4 ms | 120 ms |  | 4.7x faster |

### Construction

| case | interpx | interpx, parallel | SciPy | scirs2-interpolate | ndarray-interp | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/n10/sorted` | 93 ns | 93.9 ns | 8.75 µs | 81.9 ns | 43.5 ns | **2.1x slower** |
| `linear/n10/shuffled` | 92.4 ns | 95.4 ns | 8.67 µs | 83.3 ns | 43.4 ns | **2.1x slower** |
| `nearest/n10/sorted` | 75.6 ns | 82.4 ns | 9.48 µs | 80.7 ns |  | 1.1x faster |
| `nearest/n10/shuffled` | 74.8 ns | 75 ns | 9.46 µs | 80.6 ns |  | 1.1x faster |
| `previous/n10/sorted` | 75.1 ns | 74.6 ns | 9.22 µs |  |  | 1.2e+02x faster |
| `previous/n10/shuffled` | 82.4 ns | 76.4 ns | 9.27 µs |  |  | 1.1e+02x faster |
| `bspline1/n10/sorted` | 103 ns | 104 ns | 15.6 µs |  |  | 1.5e+02x faster |
| `bspline1/n10/shuffled` | 103 ns | 105 ns | 15.6 µs |  |  | 1.5e+02x faster |
| `bspline2/n10/sorted` | 462 ns | 459 ns | 36 µs |  |  | 78x faster |
| `bspline2/n10/shuffled` | 466 ns | 485 ns | 36 µs |  |  | 77x faster |
| `bspline3/n10/sorted` | 512 ns | 504 ns | 34.9 µs |  |  | 68x faster |
| `bspline3/n10/shuffled` | 509 ns | 507 ns | 35.4 µs |  |  | 70x faster |
| `cubic_spline/n10/sorted` | 238 ns | 241 ns | 85.7 µs | 410 ns | 488 ns | 1.7x faster |
| `cubic_spline/n10/shuffled` | 241 ns | 237 ns | 86 µs | 424 ns | 492 ns | 1.8x faster |
| `pchip/n10/sorted` | 179 ns | 196 ns | 97.8 µs | 198 ns |  | 1.1x faster |
| `pchip/n10/shuffled` | 178 ns | 195 ns | 102 µs | 209 ns |  | 1.2x faster |
| `akima/n10/sorted` | 187 ns | 186 ns | 83.9 µs | 207 ns |  | 1.1x faster |
| `akima/n10/shuffled` | 169 ns | 167 ns | 75.4 µs | 179 ns |  | 1.1x faster |
| `makima/n10/sorted` | 179 ns | 178 ns | 78.6 µs |  |  | 4.4e+02x faster |
| `makima/n10/shuffled` | 179 ns | 179 ns | 78.1 µs |  |  | 4.4e+02x faster |
| `linear/n50/sorted` | 160 ns | 158 ns | 8.82 µs | 126 ns | 111 ns | **1.4x slower** |
| `linear/n50/shuffled` | 169 ns | 155 ns | 8.86 µs | 119 ns | 101 ns | **1.7x slower** |
| `nearest/n50/sorted` | 117 ns | 108 ns | 10.7 µs | 133 ns |  | 1.1x faster |
| `nearest/n50/shuffled` | 115 ns | 107 ns | 10.7 µs | 132 ns |  | 1.1x faster |
| `previous/n50/sorted` | 106 ns | 108 ns | 9.55 µs |  |  | 90x faster |
| `previous/n50/shuffled` | 117 ns | 109 ns | 9.49 µs |  |  | 81x faster |
| `bspline1/n50/sorted` | 144 ns | 147 ns | 15.7 µs |  |  | 1.1e+02x faster |
| `bspline1/n50/shuffled` | 143 ns | 146 ns | 15.6 µs |  |  | 1.1e+02x faster |
| `bspline2/n50/sorted` | 1.96 µs | 1.84 µs | 39.2 µs |  |  | 20x faster |
| `bspline2/n50/shuffled` | 2 µs | 1.83 µs | 42.7 µs |  |  | 21x faster |
| `bspline3/n50/sorted` | 2.24 µs | 2.04 µs | 39.2 µs |  |  | 18x faster |
| `bspline3/n50/shuffled` | 2.06 µs | 2.23 µs | 39.1 µs |  |  | 19x faster |
| `cubic_spline/n50/sorted` | 837 ns | 843 ns | 87.3 µs | 854 ns | 1.09 µs | 1x faster |
| `cubic_spline/n50/shuffled` | 837 ns | 844 ns | 86.5 µs | 855 ns | 1.07 µs | 1x faster |
| `pchip/n50/sorted` | 506 ns | 504 ns | 101 µs | 334 ns |  | **1.5x slower** |
| `pchip/n50/shuffled` | 494 ns | 497 ns | 100 µs | 348 ns |  | **1.4x slower** |
| `akima/n50/sorted` | 468 ns | 427 ns | 75.1 µs | 323 ns |  | **1.5x slower** |
| `akima/n50/shuffled` | 427 ns | 476 ns | 75.2 µs | 326 ns |  | **1.3x slower** |
| `makima/n50/sorted` | 467 ns | 468 ns | 78.1 µs |  |  | 1.7e+02x faster |
| `makima/n50/shuffled` | 475 ns | 468 ns | 78 µs |  |  | 1.6e+02x faster |
| `linear/n100/sorted` | 247 ns | 248 ns | 9.08 µs | 195 ns | 169 ns | **1.5x slower** |
| `linear/n100/shuffled` | 247 ns | 247 ns | 8.98 µs | 196 ns | 169 ns | **1.5x slower** |
| `nearest/n100/sorted` | 166 ns | 161 ns | 9.89 µs | 195 ns |  | 1.2x faster |
| `nearest/n100/shuffled` | 161 ns | 158 ns | 9.9 µs | 191 ns |  | 1.2x faster |
| `previous/n100/sorted` | 174 ns | 177 ns | 9.94 µs |  |  | 57x faster |
| `previous/n100/shuffled` | 160 ns | 161 ns | 9.85 µs |  |  | 62x faster |
| `bspline1/n100/sorted` | 196 ns | 215 ns | 15.4 µs |  |  | 79x faster |
| `bspline1/n100/shuffled` | 219 ns | 199 ns | 15.8 µs |  |  | 72x faster |
| `bspline2/n100/sorted` | 3.73 µs | 3.48 µs | 45.3 µs |  |  | 12x faster |
| `bspline2/n100/shuffled` | 3.82 µs | 3.48 µs | 43.5 µs |  |  | 11x faster |
| `bspline3/n100/sorted` | 3.91 µs | 3.94 µs | 44.3 µs |  |  | 11x faster |
| `bspline3/n100/shuffled` | 4.03 µs | 3.89 µs | 43.9 µs |  |  | 11x faster |
| `cubic_spline/n100/sorted` | 1.58 µs | 1.62 µs | 87.8 µs | 1.53 µs | 1.89 µs | **1x slower** |
| `cubic_spline/n100/shuffled` | 1.58 µs | 1.73 µs | 94.9 µs | 1.52 µs | 1.88 µs | **1x slower** |
| `pchip/n100/sorted` | 898 ns | 930 ns | 102 µs | 549 ns |  | **1.6x slower** |
| `pchip/n100/shuffled` | 973 ns | 911 ns | 101 µs | 588 ns |  | **1.7x slower** |
| `akima/n100/sorted` | 741 ns | 753 ns | 75.6 µs | 538 ns |  | **1.4x slower** |
| `akima/n100/shuffled` | 741 ns | 737 ns | 76.7 µs | 531 ns |  | **1.4x slower** |
| `makima/n100/sorted` | 807 ns | 857 ns | 79.1 µs |  |  | 98x faster |
| `makima/n100/shuffled` | 821 ns | 833 ns | 79.4 µs |  |  | 97x faster |
| `linear/n1000/sorted` | 1.99 µs | 1.96 µs | 12.6 µs | 1.21 µs | 1.38 µs | **1.6x slower** |
| `linear/n1000/shuffled` | 1.95 µs | 1.95 µs | 13.7 µs | 1.23 µs | 1.38 µs | **1.6x slower** |
| `nearest/n1000/sorted` | 991 ns | 1.08 µs | 14.5 µs | 1.22 µs |  | 1.2x faster |
| `nearest/n1000/shuffled` | 984 ns | 992 ns | 14 µs | 1.22 µs |  | 1.2x faster |
| `previous/n1000/sorted` | 1.01 µs | 996 ns | 16.4 µs |  |  | 16x faster |
| `previous/n1000/shuffled` | 989 ns | 1.06 µs | 16 µs |  |  | 16x faster |
| `bspline1/n1000/sorted` | 1.03 µs | 1.04 µs | 16.8 µs |  |  | 16x faster |
| `bspline1/n1000/shuffled` | 1.09 µs | 1.13 µs | 18.5 µs |  |  | 17x faster |
| `bspline2/n1000/sorted` | 34.1 µs | 35.2 µs | 109 µs |  |  | 3.2x faster |
| `bspline2/n1000/shuffled` | 37.6 µs | 34.8 µs | 107 µs |  |  | 2.9x faster |
| `bspline3/n1000/sorted` | 64.9 µs | 65 µs | 124 µs |  |  | 1.9x faster |
| `bspline3/n1000/shuffled` | 66.3 µs | 65.7 µs | 124 µs |  |  | 1.9x faster |
| `cubic_spline/n1000/sorted` | 36.4 µs | 41.6 µs | 112 µs | 16.5 µs | 15.8 µs | **2.3x slower** |
| `cubic_spline/n1000/shuffled` | 37.5 µs | 36.3 µs | 112 µs | 16.6 µs | 15.8 µs | **2.4x slower** |
| `pchip/n1000/sorted` | 17.9 µs | 24.2 µs | 119 µs | 4.17 µs |  | **4.3x slower** |
| `pchip/n1000/shuffled` | 11.6 µs | 11.9 µs | 117 µs | 4.54 µs |  | **2.6x slower** |
| `akima/n1000/sorted` | 10 µs | 11.1 µs | 93 µs | 5.55 µs |  | **1.8x slower** |
| `akima/n1000/shuffled` | 9.47 µs | 10.6 µs | 92.3 µs | 5.7 µs |  | **1.7x slower** |
| `makima/n1000/sorted` | 10.8 µs | 11 µs | 98.2 µs |  |  | 9.1x faster |
| `makima/n1000/shuffled` | 10.5 µs | 12 µs | 98.6 µs |  |  | 9.4x faster |
| `linear/n10000/sorted` | 94 µs | 86.7 µs | 50.7 µs | 57.2 µs | 20.2 µs | **4.7x slower** |
| `linear/n10000/shuffled` | 94.5 µs | 81.7 µs | 47.8 µs | 61.6 µs | 21.6 µs | **4.4x slower** |
| `nearest/n10000/sorted` | 67.2 µs | 56.4 µs | 53.3 µs | 60.3 µs |  | **1.3x slower** |
| `nearest/n10000/shuffled` | 61.9 µs | 54.1 µs | 52.2 µs | 57.6 µs |  | **1.2x slower** |
| `previous/n10000/sorted` | 61.8 µs | 61.3 µs | 73.2 µs |  |  | 1.2x faster |
| `previous/n10000/shuffled` | 60.9 µs | 53.3 µs | 75.2 µs |  |  | 1.2x faster |
| `bspline1/n10000/sorted` | 67.2 µs | 60.2 µs | 22.8 µs |  |  | **2.9x slower** |
| `bspline1/n10000/shuffled` | 61.4 µs | 56.5 µs | 23.2 µs |  |  | **2.6x slower** |
| `bspline2/n10000/sorted` | 477 µs | 500 µs | 770 µs |  |  | 1.6x faster |
| `bspline2/n10000/shuffled` | 514 µs | 493 µs | 745 µs |  |  | 1.4x faster |
| `bspline3/n10000/sorted` | 571 µs | 610 µs | 989 µs |  |  | 1.7x faster |
| `bspline3/n10000/shuffled` | 594 µs | 596 µs | 929 µs |  |  | 1.6x faster |
| `cubic_spline/n10000/sorted` | 232 µs | 312 µs | 718 µs | 343 µs | 357 µs | 1.5x faster |
| `cubic_spline/n10000/shuffled` | 225 µs | 309 µs | 708 µs | 270 µs | 352 µs | 1.2x faster |
| `pchip/n10000/sorted` | 91 µs | 198 µs | 237 µs | 155 µs |  | 1.7x faster |
| `pchip/n10000/shuffled` | 82.8 µs | 196 µs | 237 µs | 157 µs |  | 1.9x faster |
| `akima/n10000/sorted` | 112 µs | 129 µs | 645 µs | 116 µs |  | 1x faster |
| `akima/n10000/shuffled` | 115 µs | 127 µs | 595 µs | 126 µs |  | 1.1x faster |
| `makima/n10000/sorted` | 76.4 µs | 90.4 µs | 699 µs |  |  | 9.2x faster |
| `makima/n10000/shuffled` | 83.6 µs | 91.2 µs | 693 µs |  |  | 8.3x faster |
| `linear/n100/sorted-uneven` | 249 ns | 248 ns | 8.97 µs | 198 ns | 169 ns | **1.5x slower** |
| `linear/n100/shuffled-uneven` | 248 ns | 252 ns | 9.05 µs | 196 ns | 169 ns | **1.5x slower** |
| `cubic_spline/n100/sorted-uneven` | 1.72 µs | 1.59 µs | 88.4 µs | 1.53 µs | 1.89 µs | **1.1x slower** |
| `cubic_spline/n100/shuffled-uneven` | 1.58 µs | 1.58 µs | 87.4 µs | 1.52 µs | 1.89 µs | **1x slower** |
| `pchip/n100/sorted-uneven` | 891 ns | 912 ns | 102 µs | 545 ns |  | **1.6x slower** |
| `pchip/n100/shuffled-uneven` | 892 ns | 893 ns | 102 µs | 537 ns |  | **1.7x slower** |
| `linear/n10000/sorted-uneven` | 90.1 µs | 87.6 µs | 47.1 µs | 56.2 µs | 17.6 µs | **5.1x slower** |
| `linear/n10000/shuffled-uneven` | 89.9 µs | 85.7 µs | 46.4 µs | 55.7 µs | 20.1 µs | **4.5x slower** |
| `cubic_spline/n10000/sorted-uneven` | 232 µs | 305 µs | 708 µs | 319 µs | 358 µs | 1.4x faster |
| `cubic_spline/n10000/shuffled-uneven` | 226 µs | 260 µs | 698 µs | 320 µs | 352 µs | 1.4x faster |
| `pchip/n10000/sorted-uneven` | 91.1 µs | 99.9 µs | 236 µs | 156 µs |  | 1.7x faster |
| `pchip/n10000/shuffled-uneven` | 91 µs | 82.9 µs | 235 µs | 155 µs |  | 1.7x faster |

## `PPoly`

A cubic of 55 intervals, evaluated at points in no order.

### Evaluation

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `intervals55/eval10` | 73.7 ns | 74.1 ns | 3.4 µs | 46x faster |
| `intervals55/eval100` | 659 ns | 732 ns | 4.38 µs | 6.6x faster |
| `intervals55/eval1000` | 6.53 µs | 6.73 µs | 15.8 µs | 2.4x faster |
| `intervals55/eval10000` | 66.5 µs | 67.2 µs | 329 µs | 5x faster |

### Construction

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `intervals55/eval10` | 164 ns | 167 ns | 5.12 µs | 31x faster |
| `intervals55/eval100` | 164 ns | 161 ns | 5.19 µs | 32x faster |
| `intervals55/eval1000` | 158 ns | 161 ns | 5.15 µs | 33x faster |
| `intervals55/eval10000` | 158 ns | 158 ns | 5.14 µs | 33x faster |

## `BSpline`

A cubic B-spline of 55 intervals, evaluated at points in no order.

### Evaluation

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `intervals55/eval10` | 146 ns | 147 ns | 2.93 µs | 20x faster |
| `intervals55/eval100` | 1.42 µs | 1.51 µs | 5.76 µs | 4.1x faster |
| `intervals55/eval1000` | 14.2 µs | 14.3 µs | 37.1 µs | 2.6x faster |
| `intervals55/eval10000` | 155 µs | 54.1 µs | 360 µs | 2.3x faster |

### Construction

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `intervals55/eval10` | 161 ns | 169 ns | 9.31 µs | 58x faster |
| `intervals55/eval100` | 161 ns | 162 ns | 9.18 µs | 57x faster |
| `intervals55/eval1000` | 175 ns | 163 ns | 9.16 µs | 52x faster |
| `intervals55/eval10000` | 177 ns | 167 ns | 9.26 µs | 52x faster |

## `BSplineCurve`

A cubic curve through `n` points of the plane, evaluated at `4 n` parameters in order.

### Evaluation

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `degree3/n100` | 4.04 µs | 4.23 µs | 15.8 µs | 3.9x faster |
| `degree3/n10000` | 448 µs | 151 µs | 1.34 ms | 3x faster |

### Construction

| case | interpx | interpx, parallel | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: |
| `degree3/n100` | 4.24 µs | 4.24 µs | 46.1 µs | 11x faster |
| `degree3/n10000` | 665 µs | 650 µs | 1.31 ms | 2x faster |

## Regular grid

Points in no order inside the grid. `maxN`: the first axis has `N` points and each axis half the one before (SciPy's suite); `eachN`: every axis has `N` points; `flipped`: every axis given from its high end.

### Evaluation

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | ndarray-interp | interpn | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/2d/max10/eval10/ascending` | 132 ns | 136 ns | 137 ns | 19.1 µs | 782 ns | 152 ns | 95 ns | **1.4x slower** |
| `linear/2d/max10/eval10/flipped` | 134 ns | 139 ns | 147 ns | 19.1 µs |  |  |  | 1.4e+02x faster |
| `linear/2d/max10/eval100/ascending` | 1.25 µs | 1.31 µs | 1.22 µs | 21.3 µs | 8.29 µs | 1.36 µs | 927 ns | **1.3x slower** |
| `linear/2d/max10/eval100/flipped` | 1.25 µs | 1.31 µs | 1.25 µs | 21.3 µs |  |  |  | 17x faster |
| `linear/2d/max10/eval1000/ascending` | 12.5 µs | 12.6 µs | 12.5 µs | 40.3 µs | 90.7 µs | 14 µs | 9.23 µs | **1.3x slower** |
| `linear/2d/max10/eval1000/flipped` | 13.7 µs | 13.6 µs | 13.7 µs | 40.1 µs |  |  |  | 2.9x faster |
| `linear/2d/max10/eval10000/ascending` | 125 µs | 58.5 µs | 128 µs | 360 µs | 872 µs | 238 µs | 101 µs | **1.2x slower** |
| `linear/2d/max10/eval10000/flipped` | 125 µs | 49.1 µs | 125 µs | 394 µs |  |  |  | 3.2x faster |
| `linear/2d/max40/eval10/ascending` | 159 ns | 162 ns | 156 ns | 19.2 µs | 998 ns | 156 ns | 134 ns | **1.2x slower** |
| `linear/2d/max40/eval10/flipped` | 157 ns | 163 ns | 156 ns | 19 µs |  |  |  | 1.2e+02x faster |
| `linear/2d/max40/eval100/ascending` | 1.66 µs | 1.74 µs | 1.51 µs | 23.6 µs | 10.1 µs | 1.45 µs | 1.32 µs | **1.3x slower** |
| `linear/2d/max40/eval100/flipped` | 1.58 µs | 1.74 µs | 1.66 µs | 22 µs |  |  |  | 14x faster |
| `linear/2d/max40/eval1000/ascending` | 16.5 µs | 15 µs | 15.3 µs | 46.4 µs | 105 µs | 15.2 µs | 13.3 µs | **1.2x slower** |
| `linear/2d/max40/eval1000/flipped` | 15.2 µs | 15.2 µs | 16.6 µs | 46.4 µs |  |  |  | 3.1x faster |
| `linear/2d/max40/eval10000/ascending` | 152 µs | 65 µs | 166 µs | 603 µs | 1.14 ms | 229 µs | 134 µs | **1.1x slower** |
| `linear/2d/max40/eval10000/flipped` | 153 µs | 65.1 µs | 152 µs | 567 µs |  |  |  | 3.7x faster |
| `linear/2d/max200/eval10/ascending` | 193 ns | 195 ns | 193 ns | 19.3 µs | 1.96 µs | 205 ns | 193 ns | **1x slower** |
| `linear/2d/max200/eval10/flipped` | 193 ns | 195 ns | 193 ns | 19.4 µs |  |  |  | 1e+02x faster |
| `linear/2d/max200/eval100/ascending` | 1.87 µs | 1.98 µs | 1.88 µs | 23.2 µs | 19.9 µs | 1.99 µs | 1.72 µs | **1.1x slower** |
| `linear/2d/max200/eval100/flipped` | 2.05 µs | 1.95 µs | 1.88 µs | 23.5 µs |  |  |  | 11x faster |
| `linear/2d/max200/eval1000/ascending` | 18.8 µs | 20.7 µs | 18.7 µs | 72 µs | 206 µs | 22.4 µs | 17.2 µs | **1.1x slower** |
| `linear/2d/max200/eval1000/flipped` | 18.9 µs | 19 µs | 18.8 µs | 70.7 µs |  |  |  | 3.7x faster |
| `linear/2d/max200/eval10000/ascending` | 193 µs | 74.8 µs | 189 µs | 834 µs | 2.2 ms | 446 µs | 172 µs | **1.1x slower** |
| `linear/2d/max200/eval10000/flipped` | 189 µs | 74.7 µs | 192 µs | 832 µs |  |  |  | 4.4x faster |
| `linear/3d/max10/eval10/ascending` | 201 ns | 200 ns | 198 ns | 57 µs | 3.18 µs |  | 142 ns | **1.4x slower** |
| `linear/3d/max10/eval10/flipped` | 198 ns | 205 ns | 198 ns | 56.5 µs |  |  |  | 2.8e+02x faster |
| `linear/3d/max10/eval100/ascending` | 1.91 µs | 1.97 µs | 2.09 µs | 62.2 µs | 31.8 µs |  | 1.53 µs | **1.2x slower** |
| `linear/3d/max10/eval100/flipped` | 1.91 µs | 1.97 µs | 1.95 µs | 62.8 µs |  |  |  | 33x faster |
| `linear/3d/max10/eval1000/ascending` | 18.9 µs | 19.2 µs | 20.7 µs | 118 µs | 312 µs |  | 14.1 µs | **1.3x slower** |
| `linear/3d/max10/eval1000/flipped` | 19 µs | 19.2 µs | 18.9 µs | 115 µs |  |  |  | 6.1x faster |
| `linear/3d/max10/eval10000/ascending` | 188 µs | 67.7 µs | 190 µs | 1.11 ms | 3.13 ms |  | 141 µs | **1.3x slower** |
| `linear/3d/max10/eval10000/flipped` | 187 µs | 66.2 µs | 190 µs | 1.11 ms |  |  |  | 5.9x faster |
| `linear/3d/max40/eval10/ascending` | 220 ns | 235 ns | 220 ns | 55.7 µs | 3.29 µs |  | 199 ns | **1.1x slower** |
| `linear/3d/max40/eval10/flipped` | 223 ns | 225 ns | 221 ns | 56.3 µs |  |  |  | 2.5e+02x faster |
| `linear/3d/max40/eval100/ascending` | 2.13 µs | 2.19 µs | 2.13 µs | 62.3 µs | 33.4 µs |  | 1.97 µs | **1.1x slower** |
| `linear/3d/max40/eval100/flipped` | 2.12 µs | 2.21 µs | 2.13 µs | 63 µs |  |  |  | 30x faster |
| `linear/3d/max40/eval1000/ascending` | 21 µs | 21.4 µs | 21.3 µs | 135 µs | 338 µs |  | 19.8 µs | **1.1x slower** |
| `linear/3d/max40/eval1000/flipped` | 21.3 µs | 21.5 µs | 21.4 µs | 134 µs |  |  |  | 6.3x faster |
| `linear/3d/max40/eval10000/ascending` | 213 µs | 70.7 µs | 213 µs | 1.42 ms | 3.83 ms |  | 203 µs | **1x slower** |
| `linear/3d/max40/eval10000/flipped` | 212 µs | 69 µs | 213 µs | 1.39 ms |  |  |  | 6.6x faster |
| `linear/3d/max200/eval10/ascending` | 276 ns | 281 ns | 299 ns | 56.2 µs | 4.9 µs |  | 286 ns | 1x faster |
| `linear/3d/max200/eval10/flipped` | 274 ns | 305 ns | 276 ns | 56.4 µs |  |  |  | 2.1e+02x faster |
| `linear/3d/max200/eval100/ascending` | 2.7 µs | 2.77 µs | 2.7 µs | 64.6 µs | 47 µs |  | 2.61 µs | **1x slower** |
| `linear/3d/max200/eval100/flipped` | 2.69 µs | 3.05 µs | 2.69 µs | 64.4 µs |  |  |  | 24x faster |
| `linear/3d/max200/eval1000/ascending` | 29.6 µs | 27.2 µs | 27 µs | 186 µs | 481 µs |  | 26.1 µs | **1.1x slower** |
| `linear/3d/max200/eval1000/flipped` | 26.9 µs | 27.2 µs | 26.9 µs | 186 µs |  |  |  | 6.9x faster |
| `linear/3d/max200/eval10000/ascending` | 269 µs | 93 µs | 296 µs | 1.82 ms | 4.72 ms |  | 263 µs | **1x slower** |
| `linear/3d/max200/eval10000/flipped` | 269 µs | 91.9 µs | 270 µs | 1.8 ms |  |  |  | 6.7x faster |
| `linear/4d/max10/eval10/ascending` | 273 ns | 287 ns | 283 ns | 106 µs |  |  |  | 3.9e+02x faster |
| `linear/4d/max10/eval10/flipped` | 279 ns | 267 ns | 261 ns | 104 µs |  |  |  | 3.7e+02x faster |
| `linear/4d/max10/eval100/ascending` | 2.37 µs | 2.48 µs | 2.4 µs | 106 µs |  |  |  | 45x faster |
| `linear/4d/max10/eval100/flipped` | 2.55 µs | 2.47 µs | 2.63 µs | 104 µs |  |  |  | 41x faster |
| `linear/4d/max10/eval1000/ascending` | 25.4 µs | 24.6 µs | 24.4 µs | 207 µs |  |  |  | 8.2x faster |
| `linear/4d/max10/eval1000/flipped` | 24.3 µs | 24.4 µs | 25.1 µs | 211 µs |  |  |  | 8.7x faster |
| `linear/4d/max10/eval10000/ascending` | 260 µs | 87.3 µs | 254 µs | 1.78 ms |  |  |  | 6.8x faster |
| `linear/4d/max10/eval10000/flipped` | 244 µs | 83.4 µs | 253 µs | 1.79 ms |  |  |  | 7.3x faster |
| `linear/4d/max40/eval10/ascending` | 309 ns | 315 ns | 336 ns | 93.8 µs | 5.87 µs |  | 270 ns | **1.1x slower** |
| `linear/4d/max40/eval10/flipped` | 310 ns | 311 ns | 308 ns | 95.8 µs |  |  |  | 3.1e+02x faster |
| `linear/4d/max40/eval100/ascending` | 2.92 µs | 3.09 µs | 2.98 µs | 107 µs | 58.6 µs |  | 2.7 µs | **1.1x slower** |
| `linear/4d/max40/eval100/flipped` | 2.95 µs | 3.15 µs | 2.99 µs | 107 µs |  |  |  | 36x faster |
| `linear/4d/max40/eval1000/ascending` | 32.2 µs | 29.8 µs | 29.4 µs | 262 µs | 579 µs |  | 27 µs | **1.2x slower** |
| `linear/4d/max40/eval1000/flipped` | 29.6 µs | 29.5 µs | 30.7 µs | 254 µs |  |  |  | 8.6x faster |
| `linear/4d/max40/eval10000/ascending` | 298 µs | 96.7 µs | 295 µs | 2.25 ms | 6.62 ms |  | 270 µs | **1.1x slower** |
| `linear/4d/max40/eval10000/flipped` | 296 µs | 92.9 µs | 302 µs | 2.22 ms |  |  |  | 7.5x faster |
| `linear/4d/max200/eval10/ascending` | 379 ns | 384 ns | 385 ns | 106 µs | 8.55 µs |  | 358 ns | **1.1x slower** |
| `linear/4d/max200/eval10/flipped` | 421 ns | 422 ns | 415 ns | 95.9 µs |  |  |  | 2.3e+02x faster |
| `linear/4d/max200/eval100/ascending` | 3.75 µs | 3.83 µs | 3.74 µs | 109 µs | 74.2 µs |  | 3.61 µs | **1x slower** |
| `linear/4d/max200/eval100/flipped` | 3.73 µs | 3.93 µs | 3.71 µs | 107 µs |  |  |  | 29x faster |
| `linear/4d/max200/eval1000/ascending` | 37.5 µs | 38.1 µs | 38.1 µs | 310 µs | 742 µs |  | 38.2 µs | 1x faster |
| `linear/4d/max200/eval1000/flipped` | 37.5 µs | 37.8 µs | 37.5 µs | 305 µs |  |  |  | 8.1x faster |
| `linear/4d/max200/eval10000/ascending` | 373 µs | 113 µs | 415 µs | 2.84 ms | 7.51 ms |  | 364 µs | **1x slower** |
| `linear/4d/max200/eval10000/flipped` | 372 µs | 113 µs | 376 µs | 2.78 ms |  |  |  | 7.5x faster |
| `nearest/2d/each40/eval10000/ascending` | 151 µs | 67 µs | 154 µs | 652 µs | 1.49 ms |  | 135 µs | **1.1x slower** |
| `nearest/3d/each40/eval10000/ascending` | 233 µs | 73.1 µs | 220 µs | 1.02 ms | 1.93 ms |  | 228 µs | **1x slower** |
| `pchip/2d/each40/eval10000/ascending` | 500 µs | 138 µs | 500 µs | 1.19 s |  |  |  | 2.4e+03x faster |
| `pchip/3d/each40/eval10000/ascending` | 1.58 ms | 283 µs | 1.52 ms | 3.58 s |  |  |  | 2.3e+03x faster |
| `cubic/2d/each40/eval10000/ascending` | 295 µs | 103 µs | 301 µs | 1.28 ms |  |  | 409 µs | 1.4x faster |
| `cubic/3d/each40/eval10000/ascending` | 723 µs | 176 µs | 668 µs | 2.98 ms |  |  | 794 µs | 1.1x faster |
| `cubic/2d/each10/eval1000/ascending` | 26.5 µs | 27.4 µs | 26.5 µs | 130 µs |  |  | 35.8 µs | 1.4x faster |
| `cubic/2d/each40/eval1000/ascending` | 30.1 µs | 30.9 µs | 30.2 µs | 141 µs |  |  | 40.7 µs | 1.4x faster |
| `cubic/2d/each100/eval1000/ascending` | 35.3 µs | 32.4 µs | 34.2 µs | 181 µs |  |  | 43.2 µs | 1.2x faster |
| `cubic/2d/each200/eval1000/ascending` | 34.6 µs | 34.8 µs | 34.6 µs | 198 µs |  |  | 45.8 µs | 1.3x faster |
| `cubic/2d/each400/eval1000/ascending` | 38.4 µs | 41.9 µs | 38.2 µs | 273 µs |  |  | 49.5 µs | 1.3x faster |
| `quintic/2d/each10/eval1000/ascending` | 44.1 µs | 44.3 µs | 44.1 µs | 227 µs |  |  |  | 5.2x faster |
| `quintic/2d/each40/eval1000/ascending` | 48 µs | 48.3 µs | 47.6 µs | 238 µs |  |  |  | 5x faster |
| `linear/2d/each400/eval10000/ascending` | 218 µs | 72.4 µs | 218 µs | 960 µs | 4.15 ms | 539 µs | 219 µs | 1x faster |
| `linear/2d/each800/eval10000/ascending` | 239 µs | 86.3 µs | 239 µs | 1.1 ms | 7.93 ms | 579 µs | 241 µs | 1x faster |
| `nearest/2d/auckland-ba32-0401/heldout/ascending` | 346 ms | 61 ms | 347 ms | 1.32 s | 150 s |  | 716 ms | 2.1x faster |
| `linear/2d/auckland-ba32-0401/heldout/ascending` | 370 ms | 72.6 ms | 395 ms | 1.07 s | 75.8 s | 373 ms | 841 ms | 1x faster |
| `cubic/2d/auckland-ba32-0401/heldout/ascending` | 725 ms | 140 ms | 748 ms | 33 s |  |  | 1.64 s | 2.3x faster |
| `quintic/2d/auckland-ba32-0401/heldout/ascending` | 1.3 s | 209 ms | 1.29 s | 37.8 s |  |  |  | 29x faster |
| `pchip/2d/auckland-ba32-0401/heldout/ascending` | 1.13 s | 176 ms | 1.03 s |  |  |  |  |  |

### Construction

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | ndarray-interp | interpn | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/2d/max10/eval10/ascending` | 124 ns | 134 ns | 133 ns | 5.81 µs | 145 ns | 77.6 ns | 0 | **1.8e+02x slower** |
| `linear/2d/max10/eval10/flipped` | 160 ns | 176 ns | 154 ns | 11.6 µs |  |  |  | 72x faster |
| `linear/2d/max10/eval100/ascending` | 133 ns | 122 ns | 128 ns | 5.64 µs | 145 ns | 78.6 ns | 0 | **1.9e+02x slower** |
| `linear/2d/max10/eval100/flipped` | 159 ns | 159 ns | 154 ns | 11.5 µs |  |  |  | 72x faster |
| `linear/2d/max10/eval1000/ascending` | 133 ns | 123 ns | 118 ns | 5.81 µs | 132 ns | 77.7 ns | 0 | **2e+02x slower** |
| `linear/2d/max10/eval1000/flipped` | 174 ns | 166 ns | 159 ns | 11.6 µs |  |  |  | 67x faster |
| `linear/2d/max10/eval10000/ascending` | 134 ns | 135 ns | 125 ns | 5.76 µs | 133 ns | 85.3 ns | 0 | **2e+02x slower** |
| `linear/2d/max10/eval10000/flipped` | 175 ns | 162 ns | 157 ns | 12.6 µs |  |  |  | 72x faster |
| `linear/2d/max40/eval10/ascending` | 388 ns | 360 ns | 623 ns | 5.79 µs | 199 ns | 151 ns | 0 | **5.7e+02x slower** |
| `linear/2d/max40/eval10/flipped` | 704 ns | 700 ns | 1.19 µs | 11.5 µs |  |  |  | 16x faster |
| `linear/2d/max40/eval100/ascending` | 395 ns | 395 ns | 600 ns | 5.81 µs | 219 ns | 151 ns | 0 | **5.8e+02x slower** |
| `linear/2d/max40/eval100/flipped` | 702 ns | 709 ns | 1.07 µs | 12.7 µs |  |  |  | 18x faster |
| `linear/2d/max40/eval1000/ascending` | 392 ns | 366 ns | 615 ns | 5.84 µs | 204 ns | 151 ns | 0 | **5.7e+02x slower** |
| `linear/2d/max40/eval1000/flipped` | 647 ns | 704 ns | 1.17 µs | 11.6 µs |  |  |  | 18x faster |
| `linear/2d/max40/eval10000/ascending` | 362 ns | 393 ns | 684 ns | 6.34 µs | 219 ns | 165 ns | 0 | **5.3e+02x slower** |
| `linear/2d/max40/eval10000/flipped` | 712 ns | 705 ns | 1.18 µs | 12.7 µs |  |  |  | 18x faster |
| `linear/2d/max200/eval10/ascending` | 6.86 µs | 6.96 µs | 16.3 µs | 5.84 µs | 2.99 µs | 519 ns | 0 | **8.8e+03x slower** |
| `linear/2d/max200/eval10/flipped` | 14.5 µs | 14.6 µs | 26.4 µs | 12.7 µs |  |  |  | **1.1x slower** |
| `linear/2d/max200/eval100/ascending` | 6.82 µs | 6.95 µs | 14.7 µs | 6 µs | 2.95 µs | 469 ns | 0 | **1e+04x slower** |
| `linear/2d/max200/eval100/flipped` | 15.9 µs | 14.5 µs | 26.8 µs | 12.2 µs |  |  |  | **1.3x slower** |
| `linear/2d/max200/eval1000/ascending` | 6.85 µs | 7.01 µs | 16.3 µs | 5.77 µs | 3.02 µs | 469 ns | 0 | **1e+04x slower** |
| `linear/2d/max200/eval1000/flipped` | 14.5 µs | 14.5 µs | 29.8 µs | 11.9 µs |  |  |  | **1.2x slower** |
| `linear/2d/max200/eval10000/ascending` | 44.8 µs | 45.3 µs | 15.7 µs | 5.75 µs | 51.1 µs | 486 ns | 0 | **6.5e+04x slower** |
| `linear/2d/max200/eval10000/flipped` | 52.3 µs | 54.9 µs | 28.5 µs | 11.9 µs |  |  |  | **4.4x slower** |
| `linear/3d/max10/eval10/ascending` | 196 ns | 194 ns | 189 ns | 7.69 µs | 169 ns |  | 0 | **2.2e+02x slower** |
| `linear/3d/max10/eval10/flipped` | 306 ns | 307 ns | 318 ns | 15.7 µs |  |  |  | 51x faster |
| `linear/3d/max10/eval100/ascending` | 179 ns | 181 ns | 205 ns | 7.61 µs | 164 ns |  | 0 | **1.8e+02x slower** |
| `linear/3d/max10/eval100/flipped` | 318 ns | 333 ns | 326 ns | 16 µs |  |  |  | 50x faster |
| `linear/3d/max10/eval1000/ascending` | 195 ns | 200 ns | 187 ns | 7.67 µs | 169 ns |  | 0 | **2.2e+02x slower** |
| `linear/3d/max10/eval1000/flipped` | 305 ns | 332 ns | 319 ns | 15.6 µs |  |  |  | 51x faster |
| `linear/3d/max10/eval10000/ascending` | 194 ns | 180 ns | 184 ns | 8.27 µs | 164 ns |  | 0 | **2.2e+02x slower** |
| `linear/3d/max10/eval10000/flipped` | 306 ns | 304 ns | 315 ns | 15.7 µs |  |  |  | 51x faster |
| `linear/3d/max40/eval10/ascending` | 2.77 µs | 7.43 µs | 5.39 µs | 7.7 µs | 4.02 µs |  | 0 | **2.8e+03x slower** |
| `linear/3d/max40/eval10/flipped` | 8.18 µs | 8.17 µs | 12.6 µs | 15.5 µs |  |  |  | 1.9x faster |
| `linear/3d/max40/eval100/ascending` | 4.93 µs | 15.3 µs | 5.24 µs | 7.53 µs | 9.4 µs |  | 0 | **5e+03x slower** |
| `linear/3d/max40/eval100/flipped` | 8.04 µs | 27.5 µs | 12.6 µs | 15.4 µs |  |  |  | 1.9x faster |
| `linear/3d/max40/eval1000/ascending` | 14.1 µs | 15.3 µs | 5.96 µs | 7.55 µs | 5.99 µs |  | 0 | **1.6e+04x slower** |
| `linear/3d/max40/eval1000/flipped` | 27.9 µs | 27.9 µs | 12.5 µs | 15.7 µs |  |  |  | **1.8x slower** |
| `linear/3d/max40/eval10000/ascending` | 14.2 µs | 15.1 µs | 5.89 µs | 7.6 µs | 2.58 µs |  | 0 | **1.6e+04x slower** |
| `linear/3d/max40/eval10000/flipped` | 27.9 µs | 27.2 µs | 12.4 µs | 15.6 µs |  |  |  | **1.8x slower** |
| `linear/3d/max200/eval10/ascending` | 1.73 ms | 1.73 ms | 602 µs | 7.63 µs | 1.29 ms |  | 0 | **2e+06x slower** |
| `linear/3d/max200/eval10/flipped` | 3.11 ms | 3.3 ms | 1.94 ms | 16.1 µs |  |  |  | **1.9e+02x slower** |
| `linear/3d/max200/eval100/ascending` | 1.76 ms | 1.71 ms | 593 µs | 7.62 µs | 1.18 ms |  | 0 | **2e+06x slower** |
| `linear/3d/max200/eval100/flipped` | 3.18 ms | 3.55 ms | 1.95 ms | 16 µs |  |  |  | **2e+02x slower** |
| `linear/3d/max200/eval1000/ascending` | 1.61 ms | 1.67 ms | 585 µs | 7.65 µs | 1.31 ms |  | 0 | **1.8e+06x slower** |
| `linear/3d/max200/eval1000/flipped` | 3.25 ms | 3.17 ms | 1.9 ms | 16.1 µs |  |  |  | **2e+02x slower** |
| `linear/3d/max200/eval10000/ascending` | 1.61 ms | 1.85 ms | 645 µs | 7.53 µs | 1.19 ms |  | 0 | **1.8e+06x slower** |
| `linear/3d/max200/eval10000/flipped` | 3.08 ms | 3.32 ms | 1.92 ms | 15.8 µs |  |  |  | **1.9e+02x slower** |
| `linear/4d/max10/eval10/ascending` | 235 ns | 235 ns | 248 ns | 10.6 µs |  |  |  | 45x faster |
| `linear/4d/max10/eval10/flipped` | 331 ns | 333 ns | 336 ns | 17.8 µs |  |  |  | 54x faster |
| `linear/4d/max10/eval100/ascending` | 215 ns | 215 ns | 206 ns | 9.75 µs |  |  |  | 45x faster |
| `linear/4d/max10/eval100/flipped` | 339 ns | 329 ns | 361 ns | 17.7 µs |  |  |  | 52x faster |
| `linear/4d/max10/eval1000/ascending` | 235 ns | 214 ns | 233 ns | 9.63 µs |  |  |  | 41x faster |
| `linear/4d/max10/eval1000/flipped` | 335 ns | 328 ns | 335 ns | 17.9 µs |  |  |  | 54x faster |
| `linear/4d/max10/eval10000/ascending` | 216 ns | 236 ns | 231 ns | 9.58 µs |  |  |  | 44x faster |
| `linear/4d/max10/eval10000/flipped` | 337 ns | 358 ns | 333 ns | 18 µs |  |  |  | 54x faster |
| `linear/4d/max40/eval10/ascending` | 13.1 µs | 14.5 µs | 33.3 µs | 9.5 µs | 5.45 µs |  | 1.17 ns | **1.1e+04x slower** |
| `linear/4d/max40/eval10/flipped` | 53.6 µs | 52.6 µs | 65.2 µs | 19.7 µs |  |  |  | **2.7x slower** |
| `linear/4d/max40/eval100/ascending` | 13.1 µs | 14.3 µs | 32.2 µs | 9.46 µs | 5.46 µs |  | 1.17 ns | **1.1e+04x slower** |
| `linear/4d/max40/eval100/flipped` | 56.4 µs | 51.8 µs | 64.5 µs | 19.7 µs |  |  |  | **2.9x slower** |
| `linear/4d/max40/eval1000/ascending` | 13.7 µs | 13.4 µs | 26.3 µs | 9.39 µs | 5.55 µs |  | 1.17 ns | **1.2e+04x slower** |
| `linear/4d/max40/eval1000/flipped` | 53.9 µs | 52.8 µs | 64.5 µs | 20.5 µs |  |  |  | **2.6x slower** |
| `linear/4d/max40/eval10000/ascending` | 13.1 µs | 95.6 µs | 28.9 µs | 10.2 µs | 6.04 µs |  | 1.17 ns | **1.1e+04x slower** |
| `linear/4d/max40/eval10000/flipped` | 57.2 µs | 132 µs | 75.3 µs | 19.8 µs |  |  |  | **2.9x slower** |
| `linear/4d/max200/eval10/ascending` | 44.4 ms | 48.6 ms | 15.7 ms | 9.42 µs | 36.8 ms |  | 1.17 ns | **3.8e+07x slower** |
| `linear/4d/max200/eval10/flipped` | 103 ms | 105 ms | 69.6 ms | 19.7 µs |  |  |  | **5.2e+03x slower** |
| `linear/4d/max200/eval100/ascending` | 48.3 ms | 46.9 ms | 15.6 ms | 9.37 µs | 38.4 ms |  | 1.17 ns | **4.1e+07x slower** |
| `linear/4d/max200/eval100/flipped` | 99.2 ms | 102 ms | 71.1 ms | 20 µs |  |  |  | **5e+03x slower** |
| `linear/4d/max200/eval1000/ascending` | 45.1 ms | 46.7 ms | 15.4 ms | 9.36 µs | 37.2 ms |  | 1.17 ns | **3.8e+07x slower** |
| `linear/4d/max200/eval1000/flipped` | 98.1 ms | 101 ms | 65.1 ms | 20.1 µs |  |  |  | **4.9e+03x slower** |
| `linear/4d/max200/eval10000/ascending` | 48.1 ms | 47 ms | 15.8 ms | 9.3 µs | 37 ms |  | 1.27 ns | **3.8e+07x slower** |
| `linear/4d/max200/eval10000/flipped` | 99.9 ms | 102 ms | 65.3 ms | 20.2 µs |  |  |  | **4.9e+03x slower** |
| `nearest/2d/each40/eval10000/ascending` | 638 ns | 583 ns | 1.28 µs | 5.98 µs | 274 ns |  | 0 | **9.3e+02x slower** |
| `nearest/3d/each40/eval10000/ascending` | 21.1 µs | 22.8 µs | 55.3 µs | 8.17 µs | 9.37 µs |  | 0 | **2.2e+04x slower** |
| `pchip/2d/each40/eval10000/ascending` | 6.47 µs | 6.48 µs | 6.27 µs | 6.77 µs |  |  |  | 1x faster |
| `pchip/3d/each40/eval10000/ascending` | 512 µs | 625 µs | 456 µs | 9.48 µs |  |  |  | **54x slower** |
| `cubic/2d/each40/eval10000/ascending` | 8.55 µs | 8.46 µs | 8.2 µs | 1.82 ms |  |  | 0 | **9.7e+03x slower** |
| `cubic/3d/each40/eval10000/ascending` | 398 µs | 425 µs | 299 µs | 151 ms |  |  | 1.17 ns | **3.4e+05x slower** |
| `cubic/2d/each10/eval1000/ascending` | 1.47 µs | 1.48 µs | 1.64 µs | 1.08 ms |  |  | 0 | **1.7e+03x slower** |
| `cubic/2d/each40/eval1000/ascending` | 9.14 µs | 9.18 µs | 8.2 µs | 1.81 ms |  |  | 0 | **1e+04x slower** |
| `cubic/2d/each100/eval1000/ascending` | 83.4 µs | 111 µs | 83.9 µs | 7.75 ms |  |  | 0 | **9.5e+04x slower** |
| `cubic/2d/each200/eval1000/ascending` | 178 µs | 190 µs | 209 µs | 26.1 ms |  |  | 0 | **2e+05x slower** |
| `cubic/2d/each400/eval1000/ascending` | 1.36 ms | 1.44 ms | 1.21 ms | 152 ms |  |  | 0 | **1.5e+06x slower** |
| `quintic/2d/each10/eval1000/ascending` | 2.11 µs | 2.12 µs | 1.98 µs | 3.82 ms |  |  |  | 1.8e+03x faster |
| `quintic/2d/each40/eval1000/ascending` | 12.5 µs | 12.4 µs | 12.1 µs | 7.46 ms |  |  |  | 6e+02x faster |
| `linear/2d/each400/eval10000/ascending` | 330 µs | 326 µs | 96.8 µs | 6.04 µs | 275 µs | 1.12 µs | 0 | **4.8e+05x slower** |
| `linear/2d/each800/eval10000/ascending` | 1.08 ms | 1 ms | 419 µs | 6.29 µs | 825 µs | 2.22 µs | 0 | **1.6e+06x slower** |
| `nearest/2d/auckland-ba32-0401/heldout/ascending` | 16.3 ms | 16.4 ms | 5.18 ms | 6.96 µs | 13.1 ms |  | 0 | **2.1e+07x slower** |
| `linear/2d/auckland-ba32-0401/heldout/ascending` | 16.6 ms | 15.9 ms | 5.81 ms | 6.92 µs | 13.3 ms | 8.26 µs | 0 | **2.4e+07x slower** |
| `cubic/2d/auckland-ba32-0401/heldout/ascending` | 76.7 ms | 78.2 ms | 64.9 ms | 14.6 s |  |  | 0 | **8.7e+07x slower** |
| `quintic/2d/auckland-ba32-0401/heldout/ascending` | 101 ms | 91.3 ms | 83.4 ms | 21.3 s |  |  |  | 2.1e+02x faster |
| `pchip/2d/auckland-ba32-0401/heldout/ascending` | 79.2 ms | 71.8 ms | 66.4 ms |  |  |  |  |  |

## `NearestNDInterpolator`

`n` points of the unit square, evaluated at the pixels of a square image.

### Evaluation

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: |
| `n10/image50` | 33.5 µs | 33.5 µs | 123 µs | 437 µs | 13x faster |
| `n10/image100` | 129 µs | 60.5 µs | 446 µs | 1.71 ms | 13x faster |
| `n10/image500` | 3.07 ms | 558 µs | 11.9 ms | 43.4 ms | 14x faster |
| `n100/image50` | 114 µs | 96.8 µs | 190 µs | 568 µs | 5x faster |
| `n100/image100` | 425 µs | 121 µs | 725 µs | 2.18 ms | 5.1x faster |
| `n100/image500` | 9.4 ms | 1.53 ms | 16.9 ms | 54.2 ms | 5.8x faster |
| `n1000/image50` | 324 µs | 115 µs | 323 µs | 802 µs | 2.5x faster |
| `n1000/image100` | 1.17 ms | 279 µs | 1.17 ms | 2.94 ms | 2.5x faster |
| `n1000/image500` | 23.2 ms | 4.16 ms | 23 ms | 67.2 ms | 2.9x faster |
| `n10000/image50` | 334 µs | 108 µs | 330 µs | 983 µs | 2.9x faster |
| `n10000/image100` | 1.49 ms | 281 µs | 1.51 ms | 3.89 ms | 2.6x faster |
| `n10000/image500` | 31 ms | 5.24 ms | 31.5 ms | 83.7 ms | 2.7x faster |
| `n100000/image50` | 449 µs | 114 µs | 448 µs | 1.34 ms | 3x faster |
| `n100000/image100` | 1.77 ms | 358 µs | 1.77 ms | 5.12 ms | 2.9x faster |
| `n100000/image500` | 44.6 ms | 6.92 ms | 42.4 ms | 117 ms | 2.6x faster |
| `ghcn-tmax-20250101/heldout` | 159 µs | 70.5 µs | 160 µs | 524 µs | 3.3x faster |

### Construction

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: |
| `n10/image50` | 509 ns | 508 ns | 510 ns | 9.51 µs | 19x faster |
| `n10/image100` | 554 ns | 510 ns | 521 ns | 9.8 µs | 18x faster |
| `n10/image500` | 509 ns | 554 ns | 534 ns | 9.49 µs | 19x faster |
| `n100/image50` | 3.26 µs | 3.28 µs | 3.54 µs | 15.2 µs | 4.7x faster |
| `n100/image100` | 3.31 µs | 3.3 µs | 3.36 µs | 15 µs | 4.5x faster |
| `n100/image500` | 3.3 µs | 3.41 µs | 3.57 µs | 14.9 µs | 4.5x faster |
| `n1000/image50` | 53 µs | 55.2 µs | 46.9 µs | 119 µs | 2.2x faster |
| `n1000/image100` | 55.2 µs | 54.6 µs | 55 µs | 117 µs | 2.1x faster |
| `n1000/image500` | 53.8 µs | 56.4 µs | 48.2 µs | 126 µs | 2.3x faster |
| `n10000/image50` | 737 µs | 616 µs | 607 µs | 1.66 ms | 2.2x faster |
| `n10000/image100` | 615 µs | 744 µs | 686 µs | 1.64 ms | 2.7x faster |
| `n10000/image500` | 616 µs | 632 µs | 620 µs | 1.77 ms | 2.9x faster |
| `n100000/image50` | 7.91 ms | 7.97 ms | 7.46 ms | 23.8 ms | 3x faster |
| `n100000/image100` | 7.8 ms | 7.9 ms | 7.59 ms | 24.3 ms | 3.1x faster |
| `n100000/image500` | 7.97 ms | 7.87 ms | 8.33 ms | 24.1 ms | 3x faster |
| `ghcn-tmax-20250101/heldout` | 675 µs | 721 µs | 655 µs | 1.83 ms | 2.7x faster |

## `RBFInterpolator`

`n` points of a square, evaluated at 1000 points in no order (`scattered`) or at the 10 000 pixels of an image (`image`). `global`: all points; `k50`: the 50 nearest.

### Evaluation

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/global/n10/scattered` | 14.1 µs | 15.9 µs | 14.4 µs | 61.1 µs | 180 µs | 4.3x faster |
| `linear/k50/n10/scattered` | 14.6 µs | 14.5 µs | 14.4 µs | 1.36 ms |  | 94x faster |
| `thin_plate_spline/global/n10/scattered` | 71.4 µs | 71.1 µs | 71.2 µs | 111 µs | 297 µs | 1.5x faster |
| `thin_plate_spline/k50/n10/scattered` | 71.4 µs | 71.5 µs | 71.4 µs | 1.44 ms |  | 20x faster |
| `cubic/global/n10/scattered` | 16.8 µs | 16.5 µs | 17.2 µs | 58.3 µs | 215 µs | 3.5x faster |
| `cubic/k50/n10/scattered` | 16.9 µs | 16.5 µs | 16.9 µs | 1.31 ms |  | 77x faster |
| `quintic/global/n10/scattered` | 24.4 µs | 25.4 µs | 24.8 µs | 66.3 µs | 215 µs | 2.7x faster |
| `quintic/k50/n10/scattered` | 23.9 µs | 24.6 µs | 25.2 µs | 1.41 ms |  | 59x faster |
| `multiquadric/global/n10/scattered` | 27.8 µs | 28 µs | 27.8 µs | 53.4 µs | 236 µs | 1.9x faster |
| `multiquadric/k50/n10/scattered` | 27.8 µs | 27.9 µs | 27.8 µs | 1.38 ms |  | 50x faster |
| `inverse_multiquadric/global/n10/scattered` | 37.4 µs | 37 µs | 37.4 µs | 54.7 µs | 265 µs | 1.5x faster |
| `inverse_multiquadric/k50/n10/scattered` | 37.4 µs | 41.2 µs | 37.4 µs | 1.39 ms |  | 37x faster |
| `inverse_quadratic/global/n10/scattered` | 23.8 µs | 24 µs | 23.8 µs | 46.4 µs |  | 1.9x faster |
| `inverse_quadratic/k50/n10/scattered` | 23.8 µs | 24 µs | 23.8 µs | 1.34 ms |  | 56x faster |
| `gaussian/global/n10/scattered` | 65.7 µs | 67.2 µs | 67.1 µs | 87.1 µs | 389 µs | 1.3x faster |
| `gaussian/k50/n10/scattered` | 67.9 µs | 73.8 µs | 67.2 µs | 1.4 ms |  | 21x faster |
| `linear/global/n100/scattered` | 132 µs | 57.7 µs | 136 µs | 287 µs | 1.88 ms | 2.2x faster |
| `linear/k50/n100/scattered` | 13.4 ms | 2.53 ms | 12.7 ms | 18.3 ms |  | 1.4x faster |
| `thin_plate_spline/global/n100/scattered` | 696 µs | 189 µs | 695 µs | 869 µs | 2.9 ms | 1.2x faster |
| `thin_plate_spline/k50/n100/scattered` | 23.3 ms | 4.38 ms | 22.7 ms | 23.9 ms |  | 1x faster |
| `cubic/global/n100/scattered` | 139 µs | 50.6 µs | 139 µs | 286 µs | 1.88 ms | 2.1x faster |
| `cubic/k50/n100/scattered` | 14.9 ms | 2.78 ms | 15.5 ms | 21.2 ms |  | 1.4x faster |
| `quintic/global/n100/scattered` | 164 µs | 61.9 µs | 165 µs | 300 µs | 2.06 ms | 1.8x faster |
| `quintic/k50/n100/scattered` | 16.1 ms | 3.12 ms | 16.7 ms | 22.5 ms |  | 1.4x faster |
| `multiquadric/global/n100/scattered` | 265 µs | 79.4 µs | 295 µs | 295 µs | 2.48 ms | 1.1x faster |
| `multiquadric/k50/n100/scattered` | 15.9 ms | 2.83 ms | 16 ms | 19.8 ms |  | 1.2x faster |
| `inverse_multiquadric/global/n100/scattered` | 394 µs | 116 µs | 364 µs | 391 µs | 2.84 ms | **1x slower** |
| `inverse_multiquadric/k50/n100/scattered` | 15.4 ms | 2.98 ms | 15.7 ms | 20.3 ms |  | 1.3x faster |
| `inverse_quadratic/global/n100/scattered` | 224 µs | 63.8 µs | 224 µs | 264 µs |  | 1.2x faster |
| `inverse_quadratic/k50/n100/scattered` | 13.6 ms | 2.86 ms | 14.3 ms | 18.1 ms |  | 1.3x faster |
| `gaussian/global/n100/scattered` | 615 µs | 155 µs | 615 µs | 683 µs | 3.51 ms | 1.1x faster |
| `gaussian/k50/n100/scattered` | 19.8 ms | 3.64 ms | 19.6 ms | 21 ms |  | 1.1x faster |
| `linear/global/n1000/scattered` | 1.34 ms | 308 µs | 1.34 ms | 3.65 ms | 17.6 ms | 2.7x faster |
| `linear/k50/n1000/scattered` | 13.5 ms | 2.69 ms | 13.6 ms | 32.5 ms |  | 2.4x faster |
| `thin_plate_spline/global/n1000/scattered` | 7.06 ms | 1.46 ms | 7.07 ms | 9.57 ms | 29.1 ms | 1.4x faster |
| `thin_plate_spline/k50/n1000/scattered` | 23 ms | 3.97 ms | 22.2 ms | 41.9 ms |  | 1.8x faster |
| `cubic/global/n1000/scattered` | 1.34 ms | 323 µs | 1.35 ms | 3.73 ms | 18.8 ms | 2.8x faster |
| `cubic/k50/n1000/scattered` | 15.4 ms | 3.59 ms | 15.9 ms | 36.5 ms |  | 2.4x faster |
| `quintic/global/n1000/scattered` | 1.51 ms | 340 µs | 1.66 ms | 4.09 ms | 20.8 ms | 2.7x faster |
| `quintic/k50/n1000/scattered` | 16.9 ms | 3.35 ms | 19.9 ms | 40 ms |  | 2.4x faster |
| `multiquadric/global/n1000/scattered` | 2.75 ms | 599 µs | 2.7 ms | 3.98 ms | 23.3 ms | 1.4x faster |
| `multiquadric/k50/n1000/scattered` | 16.3 ms | 3.09 ms | 17.2 ms | 36.3 ms |  | 2.2x faster |
| `inverse_multiquadric/global/n1000/scattered` | 4 ms | 820 µs | 3.77 ms | 4.91 ms | 27.9 ms | 1.2x faster |
| `inverse_multiquadric/k50/n1000/scattered` | 17.2 ms | 3.35 ms | 17.1 ms | 36.1 ms |  | 2.1x faster |
| `inverse_quadratic/global/n1000/scattered` | 2.29 ms | 535 µs | 2.24 ms | 3.59 ms |  | 1.6x faster |
| `inverse_quadratic/k50/n1000/scattered` | 15.4 ms | 2.87 ms | 15.6 ms | 33.5 ms |  | 2.2x faster |
| `gaussian/global/n1000/scattered` | 6.22 ms | 1.23 ms | 6.2 ms | 8.12 ms | 35.5 ms | 1.3x faster |
| `gaussian/k50/n1000/scattered` | 21.7 ms | 3.92 ms | 21.3 ms | 42.4 ms |  | 2x faster |
| `thin_plate_spline/global/n3000/scattered` | 21.1 ms | 3.99 ms | 20.9 ms | 35.5 ms | 104 ms | 1.7x faster |
| `thin_plate_spline/k50/n10000/scattered` | 22.7 ms | 4.07 ms | 22.6 ms | 42.6 ms |  | 1.9x faster |
| `thin_plate_spline/k50/n10000/image` | 221 ms | 38.2 ms | 221 ms | 431 ms |  | 1.9x faster |
| `thin_plate_spline/k50/n100000/scattered` | 23.1 ms | 4.61 ms | 23 ms | 44.9 ms |  | 1.9x faster |
| `thin_plate_spline/k50/n100000/image` | 232 ms | 37.5 ms | 226 ms | 444 ms |  | 1.9x faster |
| `linear/global/ghcn-tmax-20250101/heldout` | 19 ms | 4.24 ms | 18.1 ms | 55.7 ms | 258 ms | 2.9x faster |
| `linear/k50/ghcn-tmax-20250101/heldout` | 18.2 ms | 3.24 ms | 19.2 ms | 44.5 ms |  | 2.4x faster |
| `thin_plate_spline/global/ghcn-tmax-20250101/heldout` | 87.8 ms | 68.6 ms | 88.1 ms | 124 ms | 412 ms | 1.4x faster |
| `thin_plate_spline/k50/ghcn-tmax-20250101/heldout` | 31 ms | 5.19 ms | 28.6 ms | 59.9 ms |  | 1.9x faster |
| `cubic/global/ghcn-tmax-20250101/heldout` | 17.9 ms | 3.46 ms | 18.7 ms | 53.1 ms | 270 ms | 3x faster |
| `cubic/k50/ghcn-tmax-20250101/heldout` | 19.4 ms | 3.93 ms | 21.6 ms | 44.4 ms |  | 2.3x faster |
| `quintic/global/ghcn-tmax-20250101/heldout` | 21.9 ms | 4.1 ms | 20.5 ms | 54.4 ms | 279 ms | 2.5x faster |
| `quintic/k50/ghcn-tmax-20250101/heldout` | 21.2 ms | 4.11 ms | 21.2 ms | 49.3 ms |  | 2.3x faster |

### Construction

| case | interpx | interpx, parallel | interpx, `from_owned` | SciPy | scirs2-interpolate | interpx against the fastest other |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| `linear/global/n10/scattered` | 1.25 µs | 1.26 µs | 1.34 µs | 10.6 µs | 5.49 µs | 4.4x faster |
| `linear/k50/n10/scattered` | 1.35 µs | 1.24 µs | 1.3 µs | 14.8 µs |  | 11x faster |
| `thin_plate_spline/global/n10/scattered` | 1.8 µs | 1.75 µs | 1.94 µs | 12.1 µs | 6.55 µs | 3.7x faster |
| `thin_plate_spline/k50/n10/scattered` | 1.79 µs | 1.8 µs | 1.76 µs | 15.5 µs |  | 8.7x faster |
| `cubic/global/n10/scattered` | 1.71 µs | 1.57 µs | 1.68 µs | 12.8 µs | 6.01 µs | 3.5x faster |
| `cubic/k50/n10/scattered` | 1.69 µs | 1.67 µs | 1.56 µs | 15.7 µs |  | 9.3x faster |
| `quintic/global/n10/scattered` | 2.13 µs | 2.09 µs | 2.23 µs | 13.8 µs | 5.71 µs | 2.7x faster |
| `quintic/k50/n10/scattered` | 2.08 µs | 2.08 µs | 2.12 µs | 17.2 µs |  | 8.3x faster |
| `multiquadric/global/n10/scattered` | 1.48 µs | 1.34 µs | 1.33 µs | 10.7 µs | 6.51 µs | 4.4x faster |
| `multiquadric/k50/n10/scattered` | 1.32 µs | 1.31 µs | 1.35 µs | 15 µs |  | 11x faster |
| `inverse_multiquadric/global/n10/scattered` | 1.3 µs | 1.3 µs | 1.31 µs | 10.5 µs | 6.82 µs | 5.2x faster |
| `inverse_multiquadric/k50/n10/scattered` | 1.3 µs | 1.3 µs | 1.41 µs | 14.8 µs |  | 11x faster |
| `inverse_quadratic/global/n10/scattered` | 1.32 µs | 1.23 µs | 1.25 µs | 10.9 µs |  | 8.3x faster |
| `inverse_quadratic/k50/n10/scattered` | 1.24 µs | 1.23 µs | 1.27 µs | 14.9 µs |  | 12x faster |
| `gaussian/global/n10/scattered` | 1.41 µs | 1.42 µs | 1.43 µs | 10.8 µs | 7.77 µs | 5.5x faster |
| `gaussian/k50/n10/scattered` | 1.42 µs | 1.59 µs | 1.42 µs | 14.9 µs |  | 10x faster |
| `linear/global/n100/scattered` | 46.5 µs | 67.6 µs | 54 µs | 76.1 µs | 727 µs | 1.6x faster |
| `linear/k50/n100/scattered` | 3.42 µs | 3.13 µs | 3.07 µs | 21.8 µs |  | 6.4x faster |
| `thin_plate_spline/global/n100/scattered` | 72.6 µs | 106 µs | 105 µs | 113 µs | 836 µs | 1.6x faster |
| `thin_plate_spline/k50/n100/scattered` | 3.43 µs | 3.43 µs | 3.54 µs | 21.9 µs |  | 6.4x faster |
| `cubic/global/n100/scattered` | 46.9 µs | 56.9 µs | 65 µs | 80.2 µs | 736 µs | 1.7x faster |
| `cubic/k50/n100/scattered` | 3.16 µs | 3.15 µs | 3.12 µs | 22.1 µs |  | 7x faster |
| `quintic/global/n100/scattered` | 52.6 µs | 77.4 µs | 72.2 µs | 91.4 µs | 752 µs | 1.7x faster |
| `quintic/k50/n100/scattered` | 3.18 µs | 3.19 µs | 3.09 µs | 23 µs |  | 7.2x faster |
| `multiquadric/global/n100/scattered` | 51.6 µs | 73 µs | 63.3 µs | 75.8 µs | 625 µs | 1.5x faster |
| `multiquadric/k50/n100/scattered` | 3.11 µs | 3.18 µs | 3.17 µs | 20.8 µs |  | 6.7x faster |
| `inverse_multiquadric/global/n100/scattered` | 55.6 µs | 80.3 µs | 67.9 µs | 86.5 µs | 830 µs | 1.6x faster |
| `inverse_multiquadric/k50/n100/scattered` | 3.44 µs | 3.44 µs | 3.37 µs | 23 µs |  | 6.7x faster |
| `inverse_quadratic/global/n100/scattered` | 53.9 µs | 71.4 µs | 60.4 µs | 73.2 µs |  | 1.4x faster |
| `inverse_quadratic/k50/n100/scattered` | 3.21 µs | 3.21 µs | 3.34 µs | 20.8 µs |  | 6.5x faster |
| `gaussian/global/n100/scattered` | 67.8 µs | 92.4 µs | 79 µs | 107 µs | 734 µs | 1.6x faster |
| `gaussian/k50/n100/scattered` | 3.19 µs | 3.21 µs | 3.36 µs | 20.7 µs |  | 6.5x faster |
| `linear/global/n1000/scattered` | 16 ms | 15.4 ms | 17 ms | 61.7 ms | 180 ms | 3.9x faster |
| `linear/k50/n1000/scattered` | 63.3 µs | 58.9 µs | 52.8 µs | 138 µs |  | 2.2x faster |
| `thin_plate_spline/global/n1000/scattered` | 19.3 ms | 19.4 ms | 20.5 ms | 63.3 ms | 194 ms | 3.3x faster |
| `thin_plate_spline/k50/n1000/scattered` | 59.8 µs | 51.6 µs | 55 µs | 130 µs |  | 2.2x faster |
| `cubic/global/n1000/scattered` | 16.4 ms | 16.9 ms | 17 ms | 59.5 ms | 183 ms | 3.6x faster |
| `cubic/k50/n1000/scattered` | 62.7 µs | 57.1 µs | 61.6 µs | 129 µs |  | 2.1x faster |
| `quintic/global/n1000/scattered` | 17.4 ms | 17.9 ms | 17.1 ms | 59.8 ms | 185 ms | 3.4x faster |
| `quintic/k50/n1000/scattered` | 58.2 µs | 54.5 µs | 53.5 µs | 131 µs |  | 2.2x faster |
| `multiquadric/global/n1000/scattered` | 18.7 ms | 18.2 ms | 17.7 ms | 60.3 ms | 225 ms | 3.2x faster |
| `multiquadric/k50/n1000/scattered` | 61.6 µs | 54.7 µs | 53 µs | 131 µs |  | 2.1x faster |
| `inverse_multiquadric/global/n1000/scattered` | 18.7 ms | 19.5 ms | 17.9 ms | 62.3 ms | 206 ms | 3.3x faster |
| `inverse_multiquadric/k50/n1000/scattered` | 61.9 µs | 55.5 µs | 55.1 µs | 139 µs |  | 2.3x faster |
| `inverse_quadratic/global/n1000/scattered` | 17.2 ms | 17.2 ms | 16.8 ms | 60.1 ms |  | 3.5x faster |
| `inverse_quadratic/k50/n1000/scattered` | 57.9 µs | 54.1 µs | 54 µs | 130 µs |  | 2.3x faster |
| `gaussian/global/n1000/scattered` | 19.5 ms | 19.3 ms | 19.1 ms | 61.9 ms | 202 ms | 3.2x faster |
| `gaussian/k50/n1000/scattered` | 58.4 µs | 49.9 µs | 54.6 µs | 130 µs |  | 2.2x faster |
| `thin_plate_spline/global/n3000/scattered` | 391 ms | 390 ms | 394 ms | 2.09 s | 7.06 s | 5.3x faster |
| `thin_plate_spline/k50/n10000/scattered` | 764 µs | 754 µs | 659 µs | 1.66 ms |  | 2.2x faster |
| `thin_plate_spline/k50/n10000/image` | 741 µs | 736 µs | 724 µs | 1.68 ms |  | 2.3x faster |
| `thin_plate_spline/k50/n100000/scattered` | 8.69 ms | 9.11 ms | 8.6 ms | 23.8 ms |  | 2.7x faster |
| `thin_plate_spline/k50/n100000/image` | 8.46 ms | 8.71 ms | 8.25 ms | 24 ms |  | 2.8x faster |
| `linear/global/ghcn-tmax-20250101/heldout` | 18.3 s | 19.5 s | 18.4 s | 114 s | 314 s | 6.2x faster |
| `linear/k50/ghcn-tmax-20250101/heldout` | 767 µs | 811 µs | 762 µs | 1.77 ms |  | 2.3x faster |
| `thin_plate_spline/global/ghcn-tmax-20250101/heldout` | 18.8 s | 37.3 s | 18.2 s | 116 s | 273 s | 6.2x faster |
| `thin_plate_spline/k50/ghcn-tmax-20250101/heldout` | 768 µs | 790 µs | 702 µs | 1.87 ms |  | 2.4x faster |
| `cubic/global/ghcn-tmax-20250101/heldout` | 17.2 s | 15.9 s | 17.8 s | 112 s | 264 s | 6.5x faster |
| `cubic/k50/ghcn-tmax-20250101/heldout` | 762 µs | 767 µs | 695 µs | 1.68 ms |  | 2.2x faster |
| `quintic/global/ghcn-tmax-20250101/heldout` | 17.4 s | 15.8 s | 16.1 s | 116 s | 272 s | 6.7x faster |
| `quintic/k50/ghcn-tmax-20250101/heldout` | 796 µs | 773 µs | 697 µs | 1.7 ms |  | 2.1x faster |

## Memory

Megabytes. "Data read" is the process with the data of the case in memory, before anything is built. The columns after it are what was added to that: the highest while building, with the interpolator built, the highest while evaluating. "Highest" is the highest of the whole run, repeats included.

| case | library | data read | while building | built | while evaluating | highest |
| :--- | :--- | ---: | ---: | ---: | ---: | ---: |
| `oned/cubic_spline/n10000/sorted` | interpx | 3.4 | +1.9 | +1.9 | +1.9 | 5.9 |
| `oned/cubic_spline/n10000/sorted` | interpx, parallel | 3.1 | +2.0 | +2.1 | +2.2 | 6.4 |
| `oned/cubic_spline/n10000/sorted` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |
| `oned/cubic_spline/n10000/sorted` | scirs2-interpolate | 2.9 | +2.2 | +2.4 | +3.0 | 5.9 |
| `oned/cubic_spline/n10000/sorted` | ndarray-interp | 3.3 | +1.6 | +1.6 | +1.6 | 5.4 |
| `grid/linear/4d/max200/eval10000/ascending` | interpx | 216.3 | +188.4 | +188.4 | +188.5 | 404.8 |
| `grid/linear/4d/max200/eval10000/ascending` | interpx, parallel | 217.5 | +187.6 | +187.6 | +187.8 | 405.4 |
| `grid/linear/4d/max200/eval10000/ascending` | interpx, `from_owned` | 204.6 | +200.1 | +200.1 | +200.1 | 404.8 |
| `grid/linear/4d/max200/eval10000/ascending` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |
| `grid/linear/4d/max200/eval10000/ascending` | scirs2-interpolate | 227.0 | +172.5 | +177.9 | +177.9 | 404.9 |
| `grid/linear/4d/max200/eval10000/ascending` | interpn | 204.7 | +0.0 | +0.0 | +0.0 | 403.6 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpx | 4.3 | +3.6 | +4.0 | +4.0 | 8.5 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpx, parallel | 4.3 | +1.9 | +4.0 | +4.0 | 8.9 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpx, `from_owned` | 5.6 | +0.3 | +2.7 | +2.7 | 8.5 |
| `grid/cubic/2d/each400/eval1000/ascending` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |
| `grid/cubic/2d/each400/eval1000/ascending` | interpn | 5.4 | +0.1 | +0.1 | +0.1 | 5.5 |
| `grid/linear/2d/each800/eval10000/ascending` | interpx | 10.9 | +2.0 | +3.9 | +3.9 | 14.9 |
| `grid/linear/2d/each800/eval10000/ascending` | interpx, parallel | 11.4 | +0.0 | +3.5 | +3.5 | 15.4 |
| `grid/linear/2d/each800/eval10000/ascending` | interpx, `from_owned` | 11.3 | +1.8 | +3.7 | +3.7 | 15.0 |
| `grid/linear/2d/each800/eval10000/ascending` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |
| `grid/linear/2d/each800/eval10000/ascending` | scirs2-interpolate | 9.3 | +0.4 | +5.6 | +5.6 | 14.9 |
| `grid/linear/2d/each800/eval10000/ascending` | ndarray-interp | 9.3 | +0.3 | +0.3 | +0.3 | 9.8 |
| `grid/linear/2d/each800/eval10000/ascending` | interpn | 11.6 | +0.0 | +0.0 | +0.0 | 11.6 |
| `nearest_nd/n100000/image500` | interpx | 12.4 | +2.8 | +2.0 | +4.0 | 17.6 |
| `nearest_nd/n100000/image500` | interpx, parallel | 13.0 | +2.4 | +1.4 | +3.7 | 18.7 |
| `nearest_nd/n100000/image500` | interpx, `from_owned` | 11.9 | +3.7 | +2.7 | +4.6 | 18.6 |
| `nearest_nd/n100000/image500` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |
| `rbf/thin_plate_spline/global/n3000/scattered` | interpx | 4.4 | +75.9 | +3.7 | +3.7 | 80.3 |
| `rbf/thin_plate_spline/global/n3000/scattered` | interpx, parallel | 3.0 | +77.3 | +5.1 | +5.5 | 80.7 |
| `rbf/thin_plate_spline/global/n3000/scattered` | interpx, `from_owned` | 3.4 | +77.0 | +4.7 | +4.7 | 80.4 |
| `rbf/thin_plate_spline/global/n3000/scattered` | SciPy | 4.3 | +0.0 | +0.0 | +0.0 | 4.3 |
| `rbf/thin_plate_spline/global/n3000/scattered` | scirs2-interpolate | 4.0 | +288.6 | +0.5 | +0.5 | 292.5 |
| `rbf/thin_plate_spline/k50/n100000/image` | interpx | 6.7 | +3.7 | +4.7 | +5.0 | 11.9 |
| `rbf/thin_plate_spline/k50/n100000/image` | interpx, parallel | 6.8 | +4.1 | +4.6 | +5.8 | 12.6 |
| `rbf/thin_plate_spline/k50/n100000/image` | interpx, `from_owned` | 5.9 | +6.8 | +5.5 | +5.8 | 13.5 |
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
| Regular grid | 7e-04 | 7e-04 | 7e-04 | 4e-16 | 6e-16 | 7e-01 |
| `NearestNDInterpolator` | 0e+00 | 0e+00 | 0e+00 |  |  |  |
| `RBFInterpolator` | 4e+00 | 4e+00 | 4e+00 | 7e+11 |  |  |
