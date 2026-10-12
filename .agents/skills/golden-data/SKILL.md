---
name: golden-data
description: >-
  Add or regenerate SciPy 1.18 golden data for interpx tests. Use when a new
  interpolator needs numeric-parity tests, when a golden file must change, or
  when the pinned SciPy version moves.
---

# Golden data

Golden data is the committed record of what SciPy 1.18 returns, so that `cargo test` can check interpx without Python. The rules are in `.cursor/rules/scipy-parity.mdc` and `.cursor/rules/numeric-tests.mdc`. The tolerances are in `docs/design.ja.md` §8.2.

## Files

- `compare/pyproject.toml` and `compare/uv.lock`: a `uv` project (`package = false`). SciPy is `1.18.*` and Python is 3.12 or newer (SciPy 1.18 does not support 3.11). The lock records the exact SciPy and NumPy.
- `compare/common/golden.py`: `golden.case(...)` builds a case and `golden.write(name, generator, seed, cases)` checks that SciPy is 1.18, encodes non-finite numbers, and writes `compare/goldens/<name>.json` (a list of scalars on one line, a matrix one row per line). It refuses a file over 1 MB.
- `compare/common/datasets.py`: the data of `docs/design.ja.md` §8.4. Franke F1 to F6, the Genz families and Rosenbrock, the samplers (Halton, uniform, grid with holes, rbox-like geometry), and the adversarial point sets.
- `compare/common/classic_data.py`: Fritsch-Carlson (11), titanium (49), x² with a near point (12), Franke's 100 and 33 node sets, Lawson's 25 node set. Transcribed mechanically from Burkardt's MIT-licensed sources; do not retype them.
- `compare/common/bench_cases.py`: the layer-C benchmark parameters (SciPy's own benchmark suite, and the larger Rust-only sizes). Its runners come with roadmap row M5-5.
- Generators and their files:

  | Generator | File | Reference |
  | --- | --- | --- |
  | `generate_datasets.py` | `datasets.json`, `functions.json` | inputs only; values of the analytic test functions |
  | `generate_oned.py` | `oned.json`, `cubic_spline.json`, `bspline1d.json`, `calculus.json` | `numpy.interp`, `interp1d` (nearest kinds), `CubicSpline`, `PchipInterpolator`, `Akima1DInterpolator`, `CubicHermiteSpline`, `make_interp_spline` |
  | `generate_curve.py` | `curve.json` | `make_interp_spline` with several outputs |
  | `generate_grid.py` | `grid.json`, `grid_pchip.json` | `RegularGridInterpolator`, six methods, 2D to 4D; the cases chosen for `pchip`, with `PchipInterpolator` applied axis by axis where an axis has fewer than four points |
  | `generate_scattered.py` | `scattered.json` | `NearestNDInterpolator` (with the cases chosen for it: rescaled, a bound on the distance, ties marked in `tied`, repeated points), `LinearNDInterpolator`, `CloughTocher2DInterpolator` |
  | `generate_rbf.py` | `rbf.json`, `rbf_system.json` | `RBFInterpolator`, eight kernels; the parts of its system (matrix, right-hand side, shift, scale, powers, evaluation rows) from the internal functions of `_rbfinterp_pythran` and `_rbfinterp_np`, on 12 points at most |
  | `generate_tridiagonal.py` | `tridiagonal.json` | `scipy.linalg.solve_banded` on tridiagonal systems, a dense solve for the cyclic ones |
  | `generate_banded.py` | `banded.json` | LAPACK `dgbsv` on banded systems, a dense solve for the cyclic ones |

- `compare/check_goldens.py`: run last; checks the 8 MB total and that every file has a generator.
- `tests/common/mod.rs`: `load::<Case>(name)` (returns `Result`; tests `.expect` it), `golden_dir`, `Real` (a number or `"NaN"`, `"Infinity"`, `"-Infinity"`), `reals`, `Tol`, `DIRECT_F64`, `close`, `assert_close`. `tests/common/functions.rs` holds Rust copies of the analytic functions.
- `tests/goldens.rs`: checks the files themselves (shape, size, baseline, the interpolation property at data points, the Rust functions against the Python values). A new file is added to its `FILES` list.
- `just gen-goldens`: runs every generator with warnings as errors, then the check. Add one line per new generator.

## What a case stores

Beyond `expected` at the `queries`, the 1D cases store `derivative_1`, `derivative_2`, `antiderivative`, `integrals` (rows of `[a, b, value]`), and the representation itself: `breakpoints` and `coefficients` for a piecewise polynomial (SciPy's layout: one row per power, highest first), or `knots`, `control_points`, and `k` for a B-spline. Use them to test construction separately from evaluation.

## Where SciPy's default is not the reference

`RegularGridInterpolator` with `cubic` or `quintic` uses an iterative sparse solver by default. Its result is about 1e-5 away from the tensor-product spline on uniform 2D grids and up to about 1e-2 away on non-uniform and 3D grids. `expected` for these two methods comes from `cubic_legacy` and `quintic_legacy`, which interpolate along one axis after the other, and the default's result is kept as `expected_scipy_default`. interpx is compared with `expected`. A direct sparse solve (`solver=scipy.sparse.linalg.spsolve`) gives the same spline less accurately, and is only checked to be close (`docs/design.ja.md` §6.4).

`RegularGridInterpolator` with `pchip` refuses an axis of fewer than four points, and interpx takes two. For those cases `expected` is `PchipInterpolator` applied by hand, the last axis first (`fold` in `generate_grid.py`); on the grids SciPy takes, the generator asserts that the two are equal.

When a new reference looks inexact, check the interpolation property first (the reference must return the data at a data point); `tests/goldens.rs` does this for every file.

## Contract

- Run a generator with `uv run --directory compare python -X utf8 <script>` (the recipe does this).
- A golden file records: `meta` (generator, SciPy, NumPy, seed) and `cases`. A case has `name`, `reference` (the SciPy call), `options`, the inputs, and either `expected` or `expected_error` (the error class SciPy raises).
- Random inputs use a NumPy `Generator` with a fixed seed. The file does not record the Python or OS version, so regenerating on another machine gives the same bytes.
- The Rust side reads the JSON in `tests/` through `tests/common`. `cargo test` never runs Python and never writes the files. Floats are written with Python's shortest round-trip repr and read with `serde_json`'s `float_roundtrip`, so the Rust side sees the same `f64`.
- The package `include` list has `/compare/goldens/` and nothing else of `compare/`. The patterns are anchored with a leading `/`, because `tests/` alone also matches `compare/.venv/**/tests/`.
- Check regeneration: stage the files, run `just gen-goldens`, and run `git diff --stat -- compare/goldens`. It must print nothing. This holds on one machine; values that come out of LAPACK (RBF, Clough-Tocher) can differ in the last digits on another platform, so a regeneration elsewhere is reviewed as a change.
- `.gitattributes` keeps `compare/goldens/*.json` and `compare/uv.lock` with LF line endings on every platform.

## Add a case

1. Find the SciPy call that is the reference. Record which class and options in the case.
2. Cover the families `.cursor/rules/spec-checks.mdc` asks for: boundary conditions, extrapolation modes, dtypes where the type is generic, dimensions, degenerate and near-degenerate inputs, and the error cases the design names.
3. Include cases where each term of the algorithm dominates (nearly coincident knots, large aspect ratio, values over many orders of magnitude).
4. Run the generator. Read the diff of the JSON: only the intended cases are added.
5. Write the Rust test with fixed tolerances from §8.2, and do the mutation check from `numeric-tests.mdc`: copy the JSON aside, change one expected value (and, separately, one NaN query), see the test fail, and restore the copy.

## Regenerate

Regeneration is its own change. The PR says why, lists the cases whose expected values moved, and for each says whether SciPy changed or the generator did. Never edit a JSON file by hand.

## Move the SciPy version

A `design` Issue first: update §5 of `docs/design.ja.md` and the differences listed against SciPy. Then regenerate everything, and review the moved cases one by one.
