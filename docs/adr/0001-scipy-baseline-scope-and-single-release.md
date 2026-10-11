# 0001: SciPy 1.18 baseline, coverage scope, and a single v0.1.0 release

- Status: Accepted
- Date: 2026-10-11
- Issue: N/A (grill, pre-Issue)

## Context

The first design said "completely cover the SciPy summary table". The scipy.interpolate surface is large (FITPACK, `Barycentric`, `Krogh`, `NdPPoly`, and others), and several behaviors changed across SciPy versions (`RegularGridInterpolator` cubic, quintic, and pchip from 1.13; `interp1d` is legacy; `interp2d` was removed). A staged publication (0.1, 0.2, 0.3) would create compatibility duties on half-finished APIs.

## Decision

- The reference is SciPy **1.18**. Goldens and numeric tests use that version.
- The claim becomes "covers the main classes of `scipy.interpolate`". The list is the coverage matrix and the out-of-scope list in `docs/design.ja.md` §5. `Akima1D` (with makima) and `CubicHermiteSpline` are added; FITPACK, `Barycentric`, `Krogh`, `NdPPoly`, `roots`/`solve`, and derivatives on grid and scattered types are after v0.1.0 (roadmap `A-*`).
- `v0.1.0` is published once, after M1 to M5 are all done. Milestones are internal and never published.

## Consequences

A public release does not constrain the API until the main features exist. Items outside the first release are named in the roadmap, so no deferral is anonymous.
