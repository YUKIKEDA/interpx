# 0003: `Interpolator` shape, extrapolation and boundary modes, and errors

- Status: Accepted. Superseded in part by 0009: the multi-output form (`n_outputs`, `eval_into`), the `?Sized` point, the names of the `Extrapolate` and `BoundaryCondition` variants, and the list of error variants.
- Date: 2026-10-11
- Issue: N/A (grill, pre-Issue)

## Context

`type Output = T` cannot express SciPy's `values` of shape `(N, m)` or the vector output of `BSplineCurve<T, D>` (`[T; D]`). `Option` alone cannot express SciPy's `extrapolate` and `bc_type`. A `fill_value` argument on `eval_batch_into` made `None` and `NaN` inconsistent. String payloads in the error type cannot be matched. `NonMonotonic` mixed unsorted and duplicate input.

## Decision

- The trait keeps `type Point: ?Sized` and `type Output`, and adds `n_outputs()` and `eval_into(&p, &mut [T]) -> Result<bool, _>` for run-time output width.
- `eval` returns `None` outside the domain or for a NaN query. `eval_or`, `eval_batch_into` take `fill` at call time. There is no constructor `fill_value`.
- `eval_batch_into` and `eval_into` return `Result` and never `assert`.
- `Extrapolate::{None, Clamp, Extend, Periodic}` for 1D, curve, and grid. `BoundaryCondition::{NotAKnot (default), Natural, Clamped, Periodic, Derivative}` for `CubicSpline` and `BSpline1D`. Scattered types have no extrapolation mode.
- Errors are typed and `#[non_exhaustive]`: `NotSorted` and `DuplicatePoints` are separate; `NonFinite { index }`; `DegenerateGeometry(DegeneracyKind)`; `InvalidParameter { name, reason }`; `SolverFailed(SolverError)`.
- Non-finite input data is `Err(NonFinite)`; a NaN query is `None`.
- `convx::ConvexHullError` is never public. It is converted inside the crate, so `convx`'s version is not part of interpx's semver. Cospherical points are not an error. `InsufficientPoints` and `DegenerateDimension` map to `InsufficientPoints` and `DegenerateGeometry`.
- Derivative, antiderivative, and integral are part of v0.1.0 for `PPoly` and `BSpline`; `roots`/`solve` are after v0.1.0.

## Consequences

SciPy's defaults and modes have a typed counterpart. The error type can grow without a breaking change.
