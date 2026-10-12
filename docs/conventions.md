# Conventions

The directory tree and what each module is for. The meaning of the types is in `docs/design.ja.md`; the rules are in `.cursor/rules/`. A row that adds a module adds its line here in the same PR.

## Tree

```
src/
  lib.rs            crate docs and the public re-exports; no logic
  scalar.rs         Scalar: the sealed bound for f32 and f64; small integers as scalars (crate-private)
  value.rs          Value: what is interpolated, a scalar or [T; M]
  interpolator.rs   Interpolator: eval, and the provided batch methods
  error.rs          InterpolateError and ParamReason
  extrapolate.rs    Extrapolate: behavior outside the domain
  boundary.rs       BoundaryCondition and EndCondition for splines
  validate.rs       input checks shared by the 1D constructors, in the fixed order (crate-private)
  search.rs         Hint, and the interval search PPoly and BSpline share: binary search,
                    hinted search, the per-run decision of the batch, extrapolation, and what
                    each mode makes of a derivative and of a definite integral (crate-private
                    but for Hint)
  tridiagonal.rs    tridiagonal and cyclic systems: Thomas, Sherman-Morrison (crate-private)
  ppoly.rs          PPoly: piecewise polynomial, Horner evaluation, derivatives and integrals
  oned/             interpolators of data in one variable (a private module; mod.rs re-exports)
    cubic_spline.rs CubicSpline: the slopes from a tridiagonal system, then a PPoly of degree 3
    linear.rs       Linear1D: numpy.interp; holds a PPoly of degree 1
    nearest.rs      Nearest1D and NearestKind: nearest, nearest-up, previous, next
  bspline.rs        BSpline: basis functions on the stack, derivatives and integrals;
                    PPoly::from_bspline and BSpline::from_ppoly
tests/
  common/           shared by the integration tests (included with `mod common;`)
    mod.rs          golden loader, Real, Tol, close, assert_close
    functions.rs    Franke, Genz, Rosenbrock: the analytic test functions
  goldens.rs        the golden files themselves: shape, size, baseline, interpolation property
  interpolator.rs   the provided methods of Interpolator; Send + Sync of the public types
  alloc.rs          allocation counts of the evaluation APIs (no harness: a `main` function)
  ppoly.rs          PPoly: golden values, input checks, intervals, extrapolation, properties
  bspline.rs        BSpline and the conversions: golden values (1D, curves, SciPy's coefficients),
                    input checks, repeated knots, extrapolation, properties
  cubic_spline.rs   CubicSpline: SciPy's coefficients and values, input checks, few points,
                    the periodic check, end conditions, properties
  linear.rs         Linear1D: golden values, input checks, exact values at the data points,
                    extrapolation, properties
  nearest.rs        Nearest1D: golden values (bit for bit), input checks, the point each
                    kind returns, extrapolation, properties
  calculus.rs       derivatives and integrals of both types: golden values, the cases where
                    SciPy's two types differ, extrapolation, properties
benches/
  ppoly.rs          time to evaluate, differentiate, and integrate a PPoly (criterion; `just bench`)
  bspline.rs        the same for a BSpline
  linear.rs         time to evaluate a Linear1D, and the same polygon through its PPoly
  nearest.rs        time to evaluate a Nearest1D
  cubic_spline.rs   time to build a CubicSpline (the tridiagonal solvers) and to evaluate it
compare/            Python side: pinned SciPy, data sets, generators, goldens (see the golden-data skill)
docs/               design.ja.md, roadmap.md, conventions.md, adr/
```

## Where things go

- **Public surface.** Every public item is re-exported from `lib.rs`. Modules are private (`mod`, not `pub mod`), so the path of an item is `interpx::Name`.
- **Dependencies between modules.** `scalar` depends on nothing. `value` depends on `scalar`. `error`, `extrapolate` depend on nothing. `interpolator` depends on `error`. `validate` depends on `error`, `scalar`, and `value`. `search` depends on `error`, `extrapolate`, `scalar`, and `value`. `tridiagonal` depends on `scalar` and `value`. `boundary` depends on `scalar` and `value` (an end derivative has the type of a value). An interpolator module depends on these, never the other way round. `bspline` depends on `ppoly` (it holds both conversions), and `ppoly` does not depend on `bspline`.
- **Interpolators.** An interpolator built from data is a file under `src/oned/` (later `curve/`, `grid/`, `scattered/`, as in `docs/design.ja.md` §5.1). These modules are private; `lib.rs` re-exports the types. What several interpolators share (`ppoly`, `bspline`, `search`, `validate`, `tridiagonal`) stays directly under `src/`. An interpolator that is a piecewise polynomial holds a `PPoly` and lends it with `as_ppoly()`; it does not repeat `derivative`, `antiderivative`, `integrate`, or `eval_derivative`.
- **Interval search.** A type that evaluates on a grid of break points or knots calls `search::eval`, `search::eval_hinted`, and `search::eval_batch_into` with a closure for the interval. It does not write its own search, extrapolation, or batch loop. Its `eval_derivative` and `integrate` call `search::eval_derivative` and `search::integrate` in the same way, so that both follow the extrapolation mode as `eval` does. A property of the two types that both must have is tested on both in `tests/calculus.rs`.
- **Unit tests** are `#[cfg(test)] mod tests` at the end of the module they test. **Integration tests** in `tests/` use only the public API.
- **Allocation.** A new interpolator is added to `tests/alloc.rs` in the PR that adds it. The counter is global to the process, so the file is built with `harness = false` and all checks are in its `main`: the test harness has threads of its own, and an allocation on one of them during a measured region is counted (it failed a CI job once). The lints treat this file as library code, so it uses `?`, not `expect`.
- **Thread safety.** A new public type gets a `Send + Sync` assertion in the test file of its module (the core types are in `tests/interpolator.rs`).
- **Input checks.** A 1D constructor calls `validate::coordinates_and_values`; it does not write its own order of checks. Its test file has one input that breaks two adjacent steps and expects the earlier error.
- **Properties.** Property tests (`proptest`) are in the test file of the type. A method with a faster path (a hint, a batch) is compared bit for bit with the plain one.
- **Benchmarks.** A bench file ends with `criterion_group!` and `criterion_main!`, allows `missing_docs` for that macro, and uses no `unwrap`, `expect`, or `panic` (the lints apply to benches).
- **Golden data.** A new golden file is added to `FILES` in `tests/goldens.rs` and to `just gen-goldens`. One file holds 1 MB at most (`tests/goldens.rs`, `compare/common/golden.py`). A row whose cases would take a file over 90% of that splits the file by family first, in the same PR, and does not drop cases to fit (`calculus.json` and `cubic_spline.json` were split from `oned.json` this way; one generator may write several files). A test over every piecewise polynomial reads the files together with `common::load_cases(&common::PIECEWISE)`; a family that gets a file of its own is added to `PIECEWISE`. The unit tests of a crate-private module read a golden file through the same loader: they include it with `#[cfg(test)] #[path = "../tests/common/mod.rs"] mod common;` at the top level of the module file, not inside `mod tests` (an integration test cannot call a crate-private function; inside `mod tests` the path would start from a directory that does not exist, which fails on Linux and macOS).
- **Errors.** A new `InterpolateError` variant comes with the row that first returns it, together with its line in `docs/design.ja.md` §4 and its message test in `src/error.rs`.
