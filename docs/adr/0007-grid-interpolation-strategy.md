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

## Amendments

The decision above is kept as it was written. What later rows changed:

- **M2-2 (#77), 2026-10-12.** `slinear` solves nothing. It is the same function as `linear`, and `GridMethod::Linear` stands for both (`docs/design.ja.md` §6.4).
- **M2-3 (#85), 2026-10-12.** `pchip` keeps the slopes along the last axis, found at construction, beside the values: twice the memory of the grid, as for `cubic`. The first fold, along the last axis, is of the values themselves, so its slopes do not depend on the point. The folds along the earlier axes are of folded values and stay at evaluation, from a local stencil of four points an axis at most. Measured on a prototype, keeping the slopes makes an evaluation 1.6 times as fast in 2D and twice as fast in 3D and 4D (`docs/design.ja.md` §6.4 has the figures). The axis order is SciPy's, and the results are the same.
