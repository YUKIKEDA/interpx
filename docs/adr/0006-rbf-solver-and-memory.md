# 0006: RBF solver and memory

- Status: Accepted
- Date: 2026-10-11
- Issue: N/A (grill, pre-Issue)

## Context

The RBF system with polynomial terms is a saddle-point system, not positive definite, so Cholesky does not apply. A dense system needs O(N²) memory (about 800 MB at N = 10,000). SciPy 1.18 has eight kernels, and the first design listed seven.

## Decision

- Eight kernels as `RbfKernel`: `linear`, `thin_plate_spline`, `cubic`, `quintic`, `gaussian`, `multiquadric`, `inverse_multiquadric`, `inverse_quadratic`.
- `epsilon` is required for the kernels that need it, else `Err(InvalidParameter)`. `degree` defaults to the kernel's minimum; a smaller value is `Err` (SciPy warns). Polynomial terms use shifted and scaled coordinates. Defaults: `thin_plate_spline`, `smoothing = 0`, `neighbors = None`.
- Global mode uses `faer` partial-pivot LU. A singular system is `Err(SolverFailed(Singular))`. A pseudo-inverse fallback is after v0.1.0 (`A-5`).
- The O(N²) memory cost is documented, and `neighbors` is recommended for large N. There is no hard limit, but a size overflow is checked and returned as `Err`.
- `neighbors` mode solves a local problem per query using `kiddo`.

## Consequences

RBF behavior matches SciPy where it can and fails explicitly where SciPy only warns.
