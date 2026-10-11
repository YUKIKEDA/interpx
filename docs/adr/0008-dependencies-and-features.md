# 0008: Dependencies, features, and the `std` requirement

- Status: Accepted
- Date: 2026-10-11
- Issue: N/A (grill, pre-Issue)

## Context

The first design used `faer` for tridiagonal solves, listed `1d`, `curve`, and `grid` as features, put `convx`, `faer`, and `kiddo` in `default` while `scattered` also depended on them, and mentioned embedded use, while `faer` and `convx` need `std`. `convx` already depends on `faer 0.24` with `default-features = false, features = ["std"]`, and a different `faer` configuration would build two copies.

## Decision

- Tridiagonal systems use a Thomas algorithm written in this crate (periodic case: Sherman-Morrison). `faer` is used only for the dense RBF solve.
- `1d`, `curve`, and `grid` are always built and have no heavy dependency. `scattered` (`convx`, `kiddo`, `faer`) is a default feature; turning it off removes the scattered types and those dependencies.
- `rayon`, `ndarray`, and `serde` are opt-in. There is no `std` feature: the crate requires `std`. `no_std` is after v0.1.0 (`A-8`).
- `faer` is `default-features = false, features = ["std"]`, and its version matches `convx`'s (0.24).
- MSRV 1.89 (from `convx` and `kiddo`), edition 2021, license `MIT OR Apache-2.0`.
- `compare/` (SciPy and other Rust crates for comparison, goldens) is outside the root crate's dependency graph.

## Consequences

The default build is smaller for users who need only 1D and grid. The dependency graph contains one copy of `faer`.
