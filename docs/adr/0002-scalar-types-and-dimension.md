# 0002: Scalar types and dimension

- Status: Accepted
- Date: 2026-10-11
- Issue: N/A (grill, pre-Issue)

## Context

`convx` builds Delaunay triangulations from `&[f64]` with a run-time dimension. `kiddo` has a const-generic dimension and supports several scalar types. `faer` supports `f32` and `f64`. A fully generic scalar for scattered data would need silent casts to `f64`.

## Decision

- 1D, curve, and grid code is generic over `T: num_traits::Float` (`f32` and `f64`).
- Scattered code is `f64` only. A caller with `f32` data converts. There is no implicit cast.
- Dimension is `const D: usize` for grid and scattered types. Runtime-dimension wrappers are after v0.1.0 (`A-4`). The `kiddo` documentation targets low dimensions (about 10 or fewer); rustdoc states the recommended range.
- `[[f64; D]]` is passed to `convx` as a flat `&[f64]` without a copy.

## Consequences

A type mismatch is a compile error. Users who want `f32` scattered interpolation convert once at the boundary.
