# 0007: Grid interpolation strategy

- Status: Accepted
- Date: 2026-10-11
- Issue: N/A (grill, pre-Issue)

## Context

SciPy 1.18 `RegularGridInterpolator` supports `linear`, `nearest`, `slinear`, `cubic`, `quintic`, and `pchip`. SciPy rebuilds per-axis splines at evaluation, which is expensive when repeated. Spline interpolation operators along different axes commute; PCHIP derivatives are nonlinear in the values, so axis order changes the result.

## Decision

- `values` is a flat `Vec<T>` in C order (last axis fastest) with `shape: [usize; D]`. `D` is const; a run-time dimension is after v0.1.0 (`A-4`).
- `slinear`, `cubic`, and `quintic` solve the tensor-product B-spline coefficients once at construction. Evaluation reads a local `(k + 1)^D` block. Memory is proportional to the number of grid points.
- `pchip` convolves per axis lazily at evaluation, using local stencils and SciPy's axis order.
- Extrapolation follows `Extrapolate`.

## Consequences

Evaluation is cheap for the spline methods. `pchip` keeps SciPy's results at the cost of per-evaluation work.
