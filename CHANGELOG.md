# Changelog

Every version of interpx on crates.io, newest first, in the form of [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). A 0.x minor version may break the public API (what the default features export); an entry says so when it does. The MSRV is raised only by a version that names it here.

## [Unreleased]

### Changed

- **The documentation** (the README, the crate documentation, the `description`, and the `keywords`) is rewritten in interpx's own terms, and the README has a Verification section. No public item and no result changes.

## [0.1.0] - 2026-10-19

The first release.

### Computes

- **One variable**: `Linear1D` (`numpy.interp`), `Nearest1D` (the nearest, previous, and next kinds of `interp1d`), `CubicSpline`, `PchipInterpolator`, `Akima1D` (with the `makima` variant), `CubicHermiteSpline`, and `BSpline1D` (`make_interp_spline`).
- **The representations they evaluate through**: `PPoly` and `BSpline`, which can also be built directly, and give derivatives, antiderivatives, and integrals. Every interpolator of one variable but `Nearest1D` lends its own (`as_ppoly`, `as_bspline`).
- **A curve through points of several dimensions**: `BSplineCurve`.
- **Values on a rectilinear grid**: `RegularGridInterpolator`, with the methods `nearest`, `linear`, `cubic`, `quintic`, and `pchip`. `cubic` and `quintic` are solved exactly, as SciPy's `cubic_legacy` and `quintic_legacy`.
- **Scattered data** (the `scattered` feature, on by default): `NearestNDInterpolator`, `LinearNDInterpolator`, `CloughTocher2DInterpolator`, and `RBFInterpolator`, over all points or over the nearest ones.
- **Batches on several threads** (the `rayon` feature): `eval_batch_par` and `eval_batch_into_par`.
- **Storing an interpolator** (the `serde` feature): `Serialize` and `Deserialize` for every interpolator and its parameters, with a format version, and checked when read back.

### Promises

- Where SciPy 1.18 computes the same interpolant, the values are SciPy's to rounding; the tests hold them against values SciPy 1.18 produced. Where interpx computes something else, the README says what and why.
- Every interpolator owns its data and is `Send + Sync + 'static`.
- Constructors check their input and return a typed `InterpolateError`. `eval` returns `None` where there is no value.
- `eval`, `eval_or`, and `eval_batch_into` do not allocate, but for the working space of an `RBFInterpolator` over neighbors that a thread keeps.
- A parallel batch gives the serial batch's results bit for bit, whatever the number of threads.

### Requires

- Rust 1.89 or later, and `std`.
- `num-traits` and `thiserror`; with `scattered`, `convx`, `kiddo`, and `faer`; `rayon` and `serde` with their features.
- No `unsafe` in this crate.

### Measured

A comparison with SciPy 1.18, `scirs2-interpolate`, `ndarray-interp`, `interpn`, and `spade`, time and memory, rows where another library is faster included, is in [docs/performance.md](https://github.com/YUKIKEDA/interpx/blob/v0.1.0/docs/performance.md).
