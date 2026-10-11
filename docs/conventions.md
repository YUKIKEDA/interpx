# Conventions

The directory tree and what each module is for. The meaning of the types is in `docs/design.ja.md`; the rules are in `.cursor/rules/`. A row that adds a module adds its line here in the same PR.

## Tree

```
src/
  lib.rs            crate docs and the public re-exports; no logic
  scalar.rs         Scalar: the sealed bound for f32 and f64
  value.rs          Value: what is interpolated, a scalar or [T; M]
  interpolator.rs   Interpolator: eval, and the provided batch methods
  error.rs          InterpolateError
  extrapolate.rs    Extrapolate: behavior outside the domain
  boundary.rs       BoundaryCondition and EndCondition for splines
  validate.rs       input checks shared by the 1D constructors, in the fixed order (crate-private)
  ppoly.rs          PPoly and Hint: piecewise polynomial, interval search, Horner evaluation
tests/
  common/           shared by the integration tests (included with `mod common;`)
    mod.rs          golden loader, Real, Tol, close, assert_close
    functions.rs    Franke, Genz, Rosenbrock: the analytic test functions
  goldens.rs        the golden files themselves: shape, size, baseline, interpolation property
  interpolator.rs   the provided methods of Interpolator; Send + Sync of the public types
  alloc.rs          allocation counts of the evaluation APIs (one test function)
  ppoly.rs          PPoly: golden values, input checks, intervals, extrapolation, properties
benches/
  ppoly.rs          time to evaluate a PPoly (criterion; `just bench`)
compare/            Python side: pinned SciPy, data sets, generators, goldens (see the golden-data skill)
docs/               design.ja.md, roadmap.md, conventions.md, adr/
```

## Where things go

- **Public surface.** Every public item is re-exported from `lib.rs`. Modules are private (`mod`, not `pub mod`), so the path of an item is `interpx::Name`.
- **Dependencies between modules.** `scalar` depends on nothing. `value` depends on `scalar`. `error`, `extrapolate` depend on nothing. `boundary` depends on `scalar`. `interpolator` depends on `error`. `validate` depends on `error`, `scalar`, and `value`. An interpolator module (`ppoly` and those after it) depends on these, never the other way round.
- **Unit tests** are `#[cfg(test)] mod tests` at the end of the module they test. **Integration tests** in `tests/` use only the public API.
- **Allocation.** A new interpolator is added to `tests/alloc.rs` in the PR that adds it. All checks stay in the one test function, because the counter is global to the process.
- **Thread safety.** A new public type gets a `Send + Sync` assertion in the test file of its module (the core types are in `tests/interpolator.rs`).
- **Input checks.** A 1D constructor calls `validate::coordinates_and_values`; it does not write its own order of checks. Its test file has one input that breaks two adjacent steps and expects the earlier error.
- **Properties.** Property tests (`proptest`) are in the test file of the type. A method with a faster path (a hint, a batch) is compared bit for bit with the plain one.
- **Benchmarks.** A bench file ends with `criterion_group!` and `criterion_main!`, allows `missing_docs` for that macro, and uses no `unwrap`, `expect`, or `panic` (the lints apply to benches).
- **Golden data.** A new golden file is added to `FILES` in `tests/goldens.rs` and to `just gen-goldens`.
- **Errors.** A new `InterpolateError` variant comes with the row that first returns it, together with its line in `docs/design.ja.md` §4 and its message test in `src/error.rs`.
