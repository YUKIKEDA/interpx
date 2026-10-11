# 0009: Value types for multiple outputs, and the concrete trait and its types

- Status: Accepted
- Date: 2026-10-11
- Issue: #9

## Context

ADR 0003 kept a sketch of `Interpolator`: `type Point: ?Sized`, a scalar `eval`, and `n_outputs()` with `eval_into(&p, &mut [T])` for a number of output columns known only at run time. Writing it down for M1-2 showed three problems.

- With more than one output column the scalar `eval` has no meaning. The type can hold a state in which a public method cannot answer, and the only remedies are a rustdoc warning or a run-time error. `.cursor/rules/types.mdc` forbids both.
- A slice of `?Sized` points cannot exist, so the batch methods needed a `where Self::Point: Sized` clause that every caller pays for, while v0.1.0 has only `const D` dimensions (ADR 0002).
- `Extrapolate::None` reads as `Option::None` at the call site, and `BoundaryCondition::{Natural, Clamped, Derivative}` could not express SciPy's per-end form `bc_type=((1, a), (2, b))`, which the golden data already contains.

## Decision

- **A value type carries multiple outputs.** `Value<T>` is a sealed trait implemented for `T` and `[T; M]`, with a zero, `axpy` (self + k * x), and a finiteness check. Interpolators take the value type as a type parameter with the scalar as its default, for example `LinearNDInterpolator<const D: usize, V = f64>`. `Output` is `V`. `n_outputs` and `eval_into` are removed. An output width known only at run time joins the run-time dimension as roadmap row A-4.
- **`Scalar` is a sealed trait** for `f32` and `f64` (`num_traits::Float + Debug + Send + Sync`). It is the one bound on `T`, and with `Value` it makes every interpolator `Send + Sync` by construction.
- **The trait** has `type Point` (Sized) and `type Output`. `eval` is required. `eval_or`, `eval_batch`, and `eval_batch_into` are provided and may be overridden for speed. The trait is not sealed.
- **`Extrapolate::{Forbid, Clamp, Extend, Periodic}`**, with a default per type that follows SciPy 1.18 (`numpy.interp` clamps, `CubicSpline` extends, `Akima1DInterpolator` forbids, and so on).
- **`BoundaryCondition<T>::{NotAKnot, Periodic, Ends(EndCondition<T>, EndCondition<T>)}`** with `EndCondition<T>::{NotAKnot, FirstDerivative(T), SecondDerivative(T)}`. Natural and clamped are constructors. A half-periodic spline cannot be written.
- **Errors.** The enum starts with the variants M1 uses. A later row adds a variant when it first needs one (`#[non_exhaustive]`). There is no dimension-mismatch variant. Input checks run in a fixed order: length, count, finiteness, then order and duplicates.

## Consequences

A multi-output interpolator and a scalar one are different types, so a caller cannot ask a question the type cannot answer. `[T; M]` is contiguous, which matches the `PPoly` coefficient layout `[interval][power][column]`. A caller whose number of outputs is known only at run time has no API in v0.1.0 and must evaluate one column at a time.
